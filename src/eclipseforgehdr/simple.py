"""The simple stack: calibrate per frame, align, average per tier, merge on
measured exposure ratios. THE FALLBACK -- not the merge.

    read -> bias/dark, flat, demosaic, white balance, colour matrix  (per frame)
         -> one black level for the whole set, from the shortest tier
         -> align -> mean per tier
         -> exposure ratios MEASURED between adjacent tiers (a slope fit)
         -> hat-weighted linear merge
         -> importhdr.build_from_rgb: every enhancement layer, as an import

WHY IT IS THE FALLBACK AND NOT THE MERGE (0.23.5). It was promoted to the
only merge for one build, and on the reference bracket -- 12 tiers at
1.5-3x steps, 2-5 frames each -- it printed the tier boundaries as concentric
isophotes in the corona. The hat weight has hard window edges per exposure
group and nothing feathers across them; the full pipeline's photometric ladder
solve, feather and tier projection exist for exactly that. It had looked clean
on a 23-tier set with 1.2x steps because there the boundaries are too dense to
see. So: the normal path is the merge, and this is what you run when the normal
path gives you something obviously wrong and you want a picture rather than a
diagnosis, or as a control to find out whether a problem is in the merge or in
the data.

WHAT IT HAS SINCE 0.23.5 that it did not before: the per-frame instrumental
steps a Bayer raw needs -- bias, dark, flat when present, the camera's as-shot
white balance and the camera->sRGB matrix -- and raw-saturated photosites
excluded from the tier mean. A 3-plane FITS is used as written.

THE EXPOSURE LADDER IS MEASURED, never read. The header exposure is only a
grouping key. The ratio is the slope of one tier mean against the previous,
not a median of ratios: a ratio is biased by any additive error common to both
tiers, and the black level always carries the shortest tier's sky.

THE BLACK LEVEL IS ONE NUMBER, from the shortest tier's corners, used for every
frame. The sky is LEFT IN: it scales with exposure exactly as the corona does,
so it reaches the merged result as one constant per channel. Per-frame sky
removal (0.23.2-0.23.4, "Remove sky") is gone: measured on a 251-frame set it
cost 26-34% more pixel-to-pixel noise beyond 3 R and shifted the ladder by 15%.

WHAT IT DOES NOT DO, and the report says so: no hot-pixel repair, no per-tier
lunar masking, no ladder solve, no LDIC, no feather, no tier projection.
"""
from __future__ import annotations
import os

import numpy as np


# ---------------------------------------------------------------- reading ---

def _planes_direct(path):
    """A 3-plane FITS read as planes, or None if this file is not one.

    `FitsFrame` re-mosaics an RGB cube into a Bayer CFA so the rest of the app
    can treat every input alike, and the pipeline then interpolates it back.
    That round trip throws three quarters of each channel away and guesses it
    again, and the guess is gradient-directed, so it lands hardest at the lunar
    limb. Siril and PixInsight write RGB cubes. For those the planes are right
    there, and reading them is both more accurate and cheaper.
    """
    if os.path.splitext(path)[1].lower() not in (".fit", ".fits", ".fts"):
        return None
    from .fits import read_fits
    d, h = read_fits(path)
    if d.ndim != 3:
        return None
    if d.shape[0] in (3, 4):
        a = np.asarray(d[:3], np.float32)
    elif d.shape[-1] in (3, 4):
        a = np.ascontiguousarray(np.moveaxis(d[..., :3], -1, 0), np.float32)
    else:
        return None
    ped = h.get("PEDESTAL", h.get("BLKLEVEL", h.get("BLACKLEV")))
    try:
        a = a - float(ped)
    except (TypeError, ValueError):
        pass
    sat = h.get("SATURATE", h.get("DATAMAX", h.get("SATLEVEL")))
    try:
        sat = float(sat)
    except (TypeError, ValueError):
        # No saturation keyword. This used to assume 65535 whatever the file
        # held, so a float FITS scaled 0..1 got a ceiling 65535x too high and
        # every absolute threshold below it (the ratio floor, the hat) was
        # wrong by the same factor. The bit depth is the honest guess, as in
        # fits.py; for float data the data's own maximum is all there is.
        try:
            bp = int(h.get("BITPIX", 16))
        except (TypeError, ValueError):
            bp = 16
        mx = float(np.nanmax(a)) if a.size else 1.0
        sat = (65535.0 if bp == 16 else 255.0 if bp == 8 else max(mx, 1.0))
        sat = max(sat, mx, 1.0)
    return a, sat


def _border(a):
    """Widths (top, bottom, left, right) of the all-zero border a registering
    writer leaves where it shifted a frame (0 in every plane)."""
    z = np.all(a == 0, axis=0)
    rows, cols = z.all(1), z.all(0)

    def run(v):
        return int(v.size) if v.all() else int(np.argmin(v))
    return (run(rows) if rows[0] else 0, run(rows[::-1]) if rows[-1] else 0,
            run(cols) if cols[0] else 0, run(cols[::-1]) if cols[-1] else 0)


def writer_registered(files, n=9):
    """True when the frames were already registered by the program that wrote
    them: in a sample spread over the set, at least half carry a registration
    border, and the borders differ (each frame was shifted by its own amount).
    Only 3-plane FITS can carry one; anything else returns False."""
    if not files:
        return False, []
    idx = sorted(set(np.linspace(0, len(files) - 1, min(n, len(files))).astype(int)))
    found = []
    for k in idx:
        pl = _planes_direct(files[k])
        if pl is None:
            return False, []
        b = _border(pl[0])
        del pl
        if any(b):
            found.append(b)
    return (2 * len(found) >= len(idx) and len(set(found)) > 1), found


def _flag_clipped_planes(a, min_px=64):
    """Clipped pixels of a 3-plane FITS -> NaN in all planes, in place.

    A Bayer raw has its saturation flagged on the sensor's numbers (below); a
    3-plane FITS came here without that, so its clipped core went into the
    tier mean and the alignment as data (2026-10-13). The writer's scaling
    moves the ceiling per plane -- one test set's Siril frames clip at 59889 in red,
    65520 in green and 65535 in blue, and a few hot pixels sit above the red
    one -- so the ceiling is found as a spike in each plane's histogram: the
    most common value in the upper half, held by at least `min_px` pixels
    and 50x more than a typical value there. Everything from 0.5 % below it
    up is flagged (the merge's hat gives it no weight anyway). A frame with
    no clipped core has no such spike and is left alone.
    """
    clip = None
    for c in range(a.shape[0]):
        p = a[c]
        m = float(np.nanmax(p)) if p.size else 0.0
        if not np.isfinite(m) or m <= 0:
            continue
        top = p[p >= 0.5 * m]
        if top.size < min_px:
            continue
        if np.array_equal(top[:1000], np.round(top[:1000])) and m < 2 ** 24:
            lo = int(np.floor(0.5 * m))
            cnt = np.bincount(np.round(top).astype(np.int64) - lo)
            k = int(np.argmax(cnt))
            v, n = float(lo + k), int(cnt[k])
            typ = float(np.median(cnt[cnt > 0]))
        else:
            vals, cnt = np.unique(top, return_counts=True)
            k = int(np.argmax(cnt))
            v, n = float(vals[k]), int(cnt[k])
            typ = float(np.median(cnt))
        if n < min_px or n < 50.0 * typ:
            continue
        hit = p >= 0.995 * v
        clip = hit if clip is None else (clip | hit)
    if clip is not None:
        a[:, clip] = np.nan
    return a


def _load(path, demosaic_method="mhc", calib=None, seconds=None):
    """(3, H, W) float32 scene-linear in sRGB primaries, black level not yet
    removed, plus the saturation level in the same units.

    A 3-plane FITS comes back as it is: whatever wrote it already calibrated and
    colour-managed it. A Bayer raw gets, in this order, the same per-frame
    INSTRUMENTAL steps the full pipeline applies:

        raw saturation mask   decided FIRST, on the sensor's own numbers, at
                              2x2 superpixel resolution -- a flat brightens a
                              corner, a matrix mixes channels, and a threshold
                              tested after either would miss clipped pixels
        - bias - rate*t       additive, before the multiplicative one
        / flat
        demosaic
        * white balance       the camera's as-shot multipliers (pick_wb)
        @ cam2rgb             camera primaries -> sRGB primaries

    Saturated raw pixels come back as NaN in all three channels, so the tier
    mean leaves them out (it counts finite pixels) and no clipped value ever
    enters the merge. Every other stage here already treats NaN as "absent".
    """
    direct = _planes_direct(path)
    if direct is not None:
        a, sat = direct
        _flag_clipped_planes(a)
        return a, sat
    from .raw import open_frame, demosaic, pick_wb
    from .pipeline import tier_headroom
    rf = open_frame(path)
    bay = np.asarray(rf.bayer, np.float32)
    sat_raw = float(rf.sat_level)
    h2, w2 = bay.shape[0] // 2, bay.shape[1] // 2
    _b = bay[: h2 * 2, : w2 * 2].reshape(h2, 2, w2, 2)
    satm = (_b >= sat_raw).any(axis=(1, 3))
    del _b
    if calib:
        bias, rate, flat = calib.get("bias"), calib.get("rate"), calib.get("flat")
        if bias is not None and bias.shape == bay.shape:
            bay = bay - bias
        elif calib.get("black") is not None:
            # ADDITIVE BEFORE MULTIPLICATIVE. A black level that rides through
            # the flat division stops being a constant: 512 ADU divided by a
            # vignette of 0.65 in the corner is 788 there and 512 in the
            # middle, and the one constant subtracted afterwards leaves a
            # gradient. rawpy raws arrive with the camera's black already off;
            # a FITS that reports none gets this one, measured once on the
            # shortest tier's own corners -- see _calibration.
            bay = bay - np.float32(calib["black"])
        if rate is not None and rate.shape == bay.shape and seconds:
            bay = bay - rate * np.float32(seconds)
        if flat is not None and flat.shape == bay.shape:
            bay = bay / flat
    rgb = demosaic(bay, demosaic_method)
    del bay
    wb, _ = pick_wb(rf, "camera")
    rgb *= wb[None, None, :]
    rgb = (rgb.reshape(-1, 3) @ rf.cam2rgb.T).reshape(rgb.shape)
    a = np.ascontiguousarray(rgb.transpose(2, 0, 1), np.float32)
    del rgb
    if satm.any():
        full = np.zeros(a.shape[1:], bool)
        full[: h2 * 2, : w2 * 2] = np.repeat(np.repeat(satm, 2, 0), 2, 1)
        a[:, full] = np.nan
        del full
    del satm
    # the hat's upper shoulder needs the saturation level in OUTPUT units: white
    # balance and the matrix push a raw-saturated pixel above sat_raw
    return a, sat_raw * tier_headroom(wb, rf.cam2rgb)


def _calibration(folder, files, flat_dir, secs, progress, short_file=None):
    """Master bias, dark rate and flat for a Bayer folder, or None.

    Built once, before any light is decoded, through the same cached builders
    the full pipeline uses -- so a folder that already has masterflat.npy /
    masterdark.npz beside its layers reuses them. A 3-plane FITS folder gets
    None: it was calibrated by whatever wrote it.
    """
    if _planes_direct(files[0]) is not None:
        return None, {"note": "3-plane FITS: calibrated by the writer, none applied"}
    from .pipeline import resolve_flat_dir, resolve_calib_dir, workdir as _wd
    from . import flat as _flat, dark as _dark
    from .raw import open_frame
    wd = _wd(folder)
    fd = resolve_flat_dir(folder, flat_dir)
    bd = resolve_calib_dir(folder, "bias")
    dd = resolve_calib_dir(folder, "dark")
    info = {"flat_dir": fd, "bias_dir": bd, "dark_dir": dd,
            "flat_applied": False, "bias_applied": False, "dark_applied": False}
    if not (fd or bd or dd):
        return None, info
    shape = tuple(open_frame(files[0]).bayer.shape)
    out = {}
    if fd:
        try:
            m, fi = _flat.load_or_build(folder, fd, shape, progress, wd)
            if m is not None and tuple(m.shape) == shape:
                out["flat"] = np.asarray(m, np.float32)
                info["flat_applied"] = True
                info["flat"] = fi
        except Exception as e:
            progress.log(f"[warn] flat correction skipped ({e})", None)
            info["flat_error"] = str(e)
    if bd or dd:
        try:
            bias, rate, ci = _dark.load_or_build(
                folder, bd, dd, shape, progress, wd,
                light_iso=None, max_light_seconds=(max(secs) if secs else None))
            ci.pop("defect_map", None)
            if bias is not None and tuple(bias.shape) == shape:
                out["bias"] = np.asarray(bias, np.float32)
                info["bias_applied"] = True
            if rate is not None and tuple(rate.shape) == shape:
                out["rate"] = np.asarray(rate, np.float32)
                info["dark_applied"] = True
            info["calib"] = ci
        except Exception as e:
            progress.log(f"[warn] bias/dark correction skipped ({e})", None)
            info["calib_error"] = str(e)
    # A flat with no bias master: the black level has to come off before the
    # division (see _load). rawpy raws already have it off; a FITS that
    # reported none does not, and its shortest tier's corners are the one
    # place to read it with almost no sky in them.
    if "flat" in out and "bias" not in out and short_file is not None:
        try:
            rf = open_frame(short_file)
            if not getattr(rf, "black_level_reported", True):
                b = np.asarray(rf.bayer, np.float32)
                h, w = b.shape
                cy, cx = max(int(h * 0.04), 8), max(int(w * 0.04), 8)
                blk = float(np.nanmedian(np.concatenate([
                    b[:cy, :cx].ravel(), b[:cy, -cx:].ravel(),
                    b[-cy:, :cx].ravel(), b[-cy:, -cx:].ravel()])))
                out["black"] = blk
                info["black_pre_flat"] = blk
                progress.log("  pre-flat black level %.1f ADU from "
                             "shortest-tier raw corners (no bias master)"
                             % blk, None)
        except Exception as e:
            progress.log(f"[warn] pre-flat black level not measured ({e})", None)
    return (out or None), info


def _pedestal(a, frac=0.04):
    """Per-channel median of the four corners of one frame.

    On the shortest tier this is the black level with almost no sky in it,
    which is the one use it has now: measured once there, applied everywhere.
    """
    h, w = a.shape[1:]
    cy, cx = max(int(h * frac), 8), max(int(w * frac), 8)
    out = np.empty(3, np.float64)
    for c in range(3):
        p = a[c]
        out[c] = np.nanmedian(np.concatenate([
            p[:cy, :cx].ravel(), p[:cy, -cx:].ravel(),
            p[-cy:, :cx].ravel(), p[-cy:, -cx:].ravel()]))
    return out


# -------------------------------------------------------------- alignment ---

def _hipass(x, s=24.0):
    k = int(max(1, round(s)))
    c = np.cumsum(np.pad(x, ((k, k), (k, k)), mode="edge"), axis=0)
    c = (c[2 * k:, :] - c[:-2 * k, :]) / (2 * k)
    c = np.cumsum(c, axis=1)
    c = (c[:, 2 * k:] - c[:, :-2 * k]) / (2 * k)
    return x - c[:x.shape[0], :x.shape[1]]


def _feature(g):
    """Half-resolution band-passed log green for the alignment. NaN in `g`
    marks pixels that hold no usable data.

    HALF RESOLUTION ON PURPOSE: a full-res FFT pair on a 44 Mpx frame costs
    seconds each and there are tens of frames. The parabolic sub-pixel fit
    lands within about 0.1 of a half-res pixel, so 0.2 px at full scale --
    below what this path's own bilinear resampling costs anyway.

    WHAT THE CORRELATION MUST NOT SEE (2026-10-13, a FITS test set):
      * the sensor's fixed pattern -- hot pixels and pixel-scale noise. Phase
        correlation whitens the spectrum, so in a long frame the thousands of
        hot pixels outweigh the corona and the peak lands on the SENSOR's
        offset between the frames. On Siril-registered test frames that is
        the registration shift Siril applied (3, -7 px between two 10 s
        frames, the same in every patch of the sky), not a sky motion. A 3x3
        median and a 1 px blur take it out;
      * clipped pixels. The clipped core grows with the exposure, so its edge
        moves outwards on every side; read as the clip level, its two sides
        gave two opposite peaks of 10-20 px. They are no data now;
      * the black border a registering writer (Siril) leaves where it shifted
        a frame, and the strip a shift leaves empty in a tier stack.
    No-data pixels are filled from their valid neighbours and the feature
    fades to zero over ~16 px towards them, so their edges carry no signal.
    """
    from scipy.ndimage import gaussian_filter, median_filter
    h, w = (g.shape[0] // 2) * 2, (g.shape[1] // 2) * 2
    g2 = g[:h, :w].reshape(h // 2, 2, w // 2, 2).mean((1, 3), dtype=np.float32)
    bad = ~np.isfinite(g2)
    anybad = bool(bad.any())
    g2 = median_filter(np.nan_to_num(g2), 3)
    l = np.log1p(np.clip(g2, 0, None) / 30.0).astype(np.float32)
    del g2
    if anybad:
        bad = _boxmean(bad.astype(np.float64), 4) > 1e-3   # grown by 4 px
        mv = (~bad).astype(np.float32)
        fill = _boxmean(l * mv, 8) / np.maximum(_boxmean(mv, 8), 1e-6)
        l = np.where(bad, fill, l).astype(np.float32)
        del fill
    f = _hipass(gaussian_filter(l, 1.0), 24.0)
    if anybad:
        f *= np.clip((_boxmean(mv, 16) - 0.97) / 0.03, 0.0, 1.0)
    return f.astype(np.float32)


# bumped when the simple stack's result changes, so a cached one is redone
# (2: chained alignment without hot/clipped pixels, merge noise map;
#  3: a set the writer already registered is kept as written;
#  4: clipped pixels of a 3-plane FITS flagged and left out of the tier mean;
#  5: a FITS's black level measured from tier pairs)
SIMPLE_REV = 5
PEAK_MIN = 2.5  # a correlation peak must stand this far above the next one


def _boxmean(x, k):
    # float64: a float32 running sum over a frame loses the small differences
    c = np.cumsum(np.cumsum(np.pad(np.asarray(x, np.float64), ((k + 1, k), (k + 1, k)),
                                   mode="constant"), 0), 1)
    n = 2 * k + 1
    return (c[n:, n:] - c[:-n, n:] - c[n:, :-n] + c[:-n, :-n]) / (n * n)


def _pcorr(a, b, conf=False):
    """Sub-pixel shift taking b onto a. Plain phase correlation. With
    `conf`, also the height of the peak over the next-highest one at least
    3 px away (1 = two equal candidates).

    NOT the pipeline's `_semi_phase_shift`, deliberately: this path is a control
    for the normal one and should share as little with it as possible.
    """
    A = np.fft.rfft2(a)
    B = np.fft.rfft2(b)
    R = A * np.conj(B)
    m = np.abs(R)
    R = R / np.maximum(m, 1e-12 * m.max())
    c = np.fft.irfft2(R, s=a.shape)
    i = np.unravel_index(np.argmax(c), c.shape)
    sub = []
    for ax, ix in enumerate(i):
        n = c.shape[ax]
        im, ip = (ix - 1) % n, (ix + 1) % n
        y0 = c[im, i[1]] if ax == 0 else c[i[0], im]
        y2 = c[ip, i[1]] if ax == 0 else c[i[0], ip]
        d = y0 - 2 * c[i] + y2
        off = 0.0 if abs(d) < 1e-12 else 0.5 * (y0 - y2) / d
        v = ix + off
        if v > n / 2:
            v -= n
        sub.append(v)
    if not conf:
        return np.array(sub)
    c2 = np.roll(np.roll(c, 3 - i[0], 0), 3 - i[1], 1)
    c2[:7, :7] = -np.inf
    return np.array(sub), float(c[i] / max(float(c2.max()), 1e-12))


def _shift(x, dy, dx, out=None):
    """Bilinear shift. Outside the source is NaN so it can be counted out."""
    H, W = x.shape
    yy = np.arange(H, dtype=np.float32) - dy
    xx = np.arange(W, dtype=np.float32) - dx
    y0 = np.floor(yy).astype(np.int32)
    x0 = np.floor(xx).astype(np.int32)
    fy = (yy - y0)[:, None].astype(np.float32)
    fx = (xx - x0)[None, :].astype(np.float32)
    y0c = np.clip(y0, 0, H - 2)
    x0c = np.clip(x0, 0, W - 2)
    r = (x[y0c][:, x0c] * (1 - fy) * (1 - fx)
         + x[y0c][:, x0c + 1] * (1 - fy) * fx
         + x[y0c + 1][:, x0c] * fy * (1 - fx)
         + x[y0c + 1][:, x0c + 1] * fy * fx)
    r[(yy < 0) | (yy > H - 1), :] = np.nan
    r[:, (xx < 0) | (xx > W - 1)] = np.nan
    return r


def _sky_noise_point(g, n_frames):
    """One point of the camera's noise curve from one tier mean: the level and
    the per-frame variance in its darkest 15 % (the sky, where there is little
    structure to mistake for noise). Variance from each pixel's difference to
    its 3x3 neighbourhood mean, squared, median (chi2(1) median 0.455).
    Across the tiers the sky level runs from the black level to bright sky,
    so the points pin both terms of var = a*v + b."""
    v = np.asarray(g[::2, ::2], np.float64)
    fin = np.isfinite(v)
    if fin.sum() < 5000:
        return None
    vf = np.where(fin, v, 0.0)
    k, n = 1, 3
    c = np.cumsum(np.cumsum(np.pad(vf, ((k + 1, k), (k + 1, k))), 0), 1)
    sm = (c[n:, n:] - c[:-n, n:] - c[n:, :-n] + c[:-n, :-n]) / (n * n)
    cf = np.cumsum(np.cumsum(np.pad(fin.astype(np.float64), ((k + 1, k), (k + 1, k))), 0), 1)
    full = (cf[n:, n:] - cf[:-n, n:] - cf[n:, :-n] + cf[:-n, :-n]) > n * n - 0.5
    ok = fin & full
    # selected on the SMOOTHED level: selecting on the pixel itself picks the
    # low tail of the noise where the sky is pure noise, and that tail sits
    # far below its neighbours, which reads as extra variance
    lo = np.quantile(sm[ok], 0.15)
    m = ok & (sm <= lo)
    if m.sum() < 2000:
        return None
    var = float(np.median(((v - sm) ** 2)[m]) * (9.0 / 8.0) / 0.455) * float(n_frames)
    return float(np.median(sm[m])), var


def _fit_noise_curve(points):
    """var_frame = a*v + b through the tiers' sky points; a, b >= 0."""
    if not points:
        return None
    X = np.array([p_[0] for p_ in points], np.float64)
    Y = np.array([p_[1] for p_ in points], np.float64)
    if len(X) >= 3 and np.ptp(X) > 0:
        # robust: least squares, then again without points > 3x off
        a, b = np.polyfit(X, Y, 1)
        r = Y / np.maximum(a * X + b, 1e-12)
        keep = (r > 1 / 3) & (r < 3)
        if keep.sum() >= 3:
            a, b = np.polyfit(X[keep], Y[keep], 1)
    else:
        a, b = 0.0, float(np.median(Y))
    if not (np.isfinite(a) and np.isfinite(b)):
        return None
    if a < 0:
        a, b = 0.0, float(np.median(Y))
    b = max(float(b), 0.05 * float(np.min(Y)))
    return float(a), b


# ------------------------------------------------------------------ merge ---

def _adu_scale(sat):
    """The floor constants below (60, 200, 400) are ADU of a 14-16-bit sensor.
    A float FITS bracket in [0, 1] (sat 1.0) met them as 'everything is below
    the floor' and merged to all zeros -- "could not find the lunar limb"
    (audit 2026-10-10, A15). Below a saturation of 4000 the data are taken as
    scaled, and the constants with them."""
    return 1.0 if sat is None or sat >= 4000.0 else float(sat) / 65535.0


def _ratio(a, b, sat, lo=200.0):
    """The exposure ratio between two tier means: the SLOPE of b against a.

    This was median(b / a) over pixels well exposed in both. A ratio of two
    values is biased by any common additive error -- and there always is one,
    because the black level comes from the shortest tier's corners and carries
    that tier's own sky along with it. Subtract the same offset d from both and
        b / a = (k A - d) / (A - d)  >  k       for d > 0,
    worst where A is small, which is also where most of the qualifying pixels
    are. On a synthetic bracket with a known 5.000x step the median read 5.34 /
    5.45 / 5.10. A straight line b = k a + c absorbs the offset into c and
    leaves k unbiased: the same bracket then measures 5.00x within 1%.

    Least squares on the masked pixels, then once more without the residual
    outliers (misregistration, a star, the odd hot photosite). Measured on
    green and applied to all three channels: on a 23-tier set the channels'
    own ratios agree to under 1.2% per step, so per-channel freedom buys
    nothing and costs a noisier estimate.
    """
    lo = lo * _adu_scale(sat)
    m = (np.isfinite(a) & np.isfinite(b) & (a > lo) & (b > lo)
         & (a < sat * 0.8) & (b < sat * 0.8))
    n = int(m.sum())
    if n < 2000:
        return np.nan, n
    x = a[m].astype(np.float64)
    y = b[m].astype(np.float64)
    if x.size > 400000:                       # one number; a sample is plenty
        idx = np.random.default_rng(0).choice(x.size, 400000, replace=False)
        x, y = x[idx], y[idx]
    for _pass in range(2):
        A = np.vstack([x, np.ones_like(x)]).T
        (k, c), *_ = np.linalg.lstsq(A, y, rcond=None)
        res = y - (k * x + c)
        sd = 1.4826 * float(np.median(np.abs(res - np.median(res)))) + 1e-9
        keep = np.abs(res) < 3.0 * sd
        if keep.sum() < 1000 or keep.all():
            break
        x, y = x[keep], y[keep]
    if not np.isfinite(k) or k <= 0:
        return np.nan, n
    return float(k), n


def _weight(v, sat):
    """Hat: zero at the noise floor, zero towards saturation.

    Weighting by exposure time instead -- what the normal merge does -- gives a
    short tier real weight in the outer field where it holds nothing but noise
    and whatever the pedestal step missed. This does not.

    THE UPPER SHOULDER IS WIDE (0.23.8). It used to run from 0.85 to 0.97 of
    saturation: a tier dropped out of the merge over a brightness factor of
    1.14, i.e. over a dozen pixels of radius in the inner corona. Every tier
    is off the others by a fraction of a percent (ratio error, offset, the
    sensor's own bend below saturation), so each drop-out printed a step of
    ~0.1% along that tier's saturation isophote. Measured on a 23-tier 16-bit
    FITS set: the merged log-luminance, binned by its own level, showed dips
    at exactly the ladder spacing (0.17 dex for the 1.5x steps), and MGN/NAFE
    plus a texture slider made them visible as faint arcs. Now a tier fades
    out from 0.35 to 0.90 of saturation -- a factor 2.6, several tiers at a
    time -- so the same mismatch becomes a gradient no filter can find. The
    upper 10% below saturation is left out entirely.
    """
    f = _adu_scale(sat)
    return (np.clip((v - 60.0 * f) / (400.0 * f), 0, 1)
            * np.clip((sat * 0.90 - v) / (sat * 0.55), 0, 1)).astype(np.float32)


# -------------------------------------------------------------------- run ---

def run(folder, progress, denoise="off", demosaic_method="mhc",
        fnrgf_preset="ours", flat_dir="", partialconv=True, tier_mode="exposure",
        align_mode="auto"):
    """Stack `folder` and build every layer the renderer needs.

    `align_mode`: "auto" keeps a set its writer already registered and aligns
    anything else; "align" always aligns; "keep" never does."""
    from .raw import list_raws, read_exif
    from . import importhdr
    from . import __version__

    try:
        progress.bar_detail = 0.45
    except Exception:
        pass

    files = list_raws(folder)
    if not files:
        raise RuntimeError("no readable frames in this folder")

    # The HEADER exposure is used ONLY as a grouping key. Every photometric
    # number below comes from the pixels.
    tiers = {}
    skipped = 0
    frame_sec, frame_ts = {}, {}
    for p in files:
        try:
            _ex = read_exif(p)
            sec = float(_ex[0])
        except Exception:
            skipped += 1
            continue
        frame_sec[p], frame_ts[p] = sec, _ex[2]
        tiers.setdefault(sec, []).append(p)
    _ginfo = None
    if str(tier_mode or "exposure").strip().lower() == "frame":
        progress.log("Tiers 'Per frame' is not used by the simple stack; "
                     "grouped by shutter speed", None)
    if str(tier_mode or "exposure").strip().lower() == "measured":
        # tiers by measured brightness (grouping.py); the frames keep their
        # header times for the dark current
        from .grouping import by_brightness
        tiers, _ginfo = by_brightness([p for p in frame_sec if frame_sec[p] > 0],
                                      frame_sec, frame_ts, progress)
    secs = sorted(tiers)
    if len(secs) < 2:
        raise RuntimeError(
            "the merge needs at least two exposure groups; this folder "
            "has %d" % len(secs))
    n_used = sum(len(v) for v in tiers.values())
    progress.log(f"{n_used} frames in {len(secs)} exposure groups"
                 + (f" ({skipped} unreadable, skipped)" if skipped else ""), 0.02)

    calib, calib_info = _calibration(folder, files, flat_dir, secs, progress,
                                     short_file=sorted(tiers[secs[0]])[0])
    if calib:
        progress.log("per-frame calibration: "
                     + ", ".join(k for k in ("bias", "rate", "flat") if k in calib)
                     .replace("rate", "dark"), None)
    is_planes = _planes_direct(files[0]) is not None
    # FRAMES THE WRITER ALREADY REGISTERED (2026-10-13, a Siril test set). A
    # registering program (Siril, PixInsight) pads every shifted frame with a
    # black border of its own width. Such a set is kept as written: the
    # writer aligned each frame on its full data, and re-aligning it here only
    # adds error -- on that set the chained correlation moved the 6.5 s frame
    # 6 px off a registration that holds to ~2 px (red/cyan overlay check).
    _am = str(align_mode or "auto").strip().lower()
    if _am not in ("auto", "align", "keep"):
        _am = "auto"
    if _am == "keep":
        prealigned, _pre_b = True, []
        progress.log("alignment off (setting): the frames are stacked as they are", None)
    elif _am == "align":
        prealigned, _pre_b = False, []
    else:
        prealigned, _pre_b = writer_registered(files) if is_planes else (False, [])
    if prealigned and _am == "auto":
        progress.log(f"frames already registered by the program that wrote them "
                     f"({len(_pre_b)} of the sampled frames carry a registration "
                     f"border, widths differ): kept as written, no alignment", None)
    progress.log("simple stack (fallback path)", None)
    progress.log("3-plane FITS, used as written" if is_planes else
                 "Bayer raw: demosaic, camera WB and colour matrix "
                 "per frame", None)
    progress.log("no hot-pixel repair, ladder solve, LDIC or feather; "
                 "ratios from pixels; sky left in", None)

    black = None
    # A 3-PLANE FITS THAT REPORTS NO BLACK LEVEL: measured from pairs of tiers
    # (pipeline.measure_black_from_tiers), as the full pipeline does. The
    # shortest frame's corners carry a little sky and a writer's border can
    # reach them; on a FITS test set they read 6-17 ADU above the measured floors.
    _black_meas = None
    if is_planes:
        try:
            from .raw import open_frame as _of
            if not getattr(_of(files[0]), "black_level_reported", True):
                from .pipeline import measure_black_from_tiers
                _black_meas, _bi = measure_black_from_tiers(tiers, secs)
                if _black_meas is None:
                    progress.log(f"  black level not measurable from the tiers "
                                 f"({_bi.get('why')}): shortest frame's corners used", None)
        except Exception as _e:
            _black_meas = None
            progress.log(f"[odd] black level measurement skipped ({_e})", None)
    ref_feat = None
    # ALIGNMENT IS CHAINED THROUGH THE TIERS (2026-10-13). Every frame used to be
    # tied to the first frame of the shortest tier. A long frame shares almost
    # nothing with it -- its inner corona is clipped where the short frame has
    # all its structure, and the short frame is noise where the long one has
    # its streamers -- so on a FITS test set the longest frame did not correlate at
    # all (shift > 400 px, set to 0) and the long tiers came out 11-25 px rms.
    # Each tier is now aligned to the previous tier's finished stack, which is
    # already in the reference frame's coordinates and is the closest exposure
    # there is. The feature ignores hot pixels and clipped pixels (_feature).
    # A peak that does not stand PEAK_MIN times above the next candidate is
    # not trusted: the frame keeps the previous tier's shift. On
    # Siril-registered test frames that is what the long single frames get (their
    # unclipped field is faint, flat outer corona); on a bracket with
    # structure every tier peaks 3-12x above the next.
    tier_feat = None
    last_shift = (0.0, 0.0)
    n_fallback = 0
    vnum = None
    num = den = None
    prev_g = None
    sat = None
    rel = 1.0
    short_acc, short_n = None, 0
    ladder = []
    noise_models = []
    for i, s in enumerate(secs):
        fs = sorted(tiers[s])
        acc = cnt = None
        clipany = None
        shifts = []
        peaks = []
        for p in fs:
            f, sl = _load(p, demosaic_method, calib=calib,
                          seconds=(frame_sec.get(p) or s) if frame_sec.get(p, 0) > 0 else s)
            if sat is None:
                sat = sl
            if black is None and _black_meas is not None:
                black = np.asarray(_black_meas, np.float32)
                progress.log("  black level (measured from tier pairs, the file "
                             "reports none): R %.1f  G %.1f  B %.1f" % tuple(black), None)
            if black is None:
                # ONE black level for the whole set, from the exposure where
                # the sky contributes least. See the module docstring.
                black = _pedestal(f)
                progress.log("  black level (shortest tier, frame 1): "
                             "R %.1f  G %.1f  B %.1f" % tuple(black), None)
            # a writer's registration border: exactly 0 in every plane before
            # the black level came off (Siril pads a shifted frame with 0).
            # No data -- out of the tier mean and out of the alignment feature
            _zb = np.all(f == 0, axis=0)
            f -= black[:, None, None]
            f[:, _zb] = np.nan
            # clipped and border pixels are NaN by now: no data for the
            # alignment either (see _feature)
            q = None
            # clipped pixels are NaN by now (raw and FITS alike) and stay out
            # of the tier mean; for the ALIGNMENT they read as the clip level,
            # which is what the sensor recorded: the edge of the clipped core
            # is the strongest thing to align on, and on a centred corona it
            # grows evenly with the exposure. Only a writer's border is no data.
            _clip = ~np.isfinite(f[1]) & ~_zb
            ft = None if prealigned else _feature(
                np.where(_clip, np.float32(sat - black[1]), f[1]))
            if prealigned:
                dy = dx = 0.0
                q = None
            elif ref_feat is None:
                ref_feat = ft
                dy = dx = 0.0
            else:
                if ft.shape != ref_feat.shape:
                    raise RuntimeError(
                        "the frames in this folder are not all the same size")
                (dy, dx), q = _pcorr(tier_feat if tier_feat is not None else ref_feat,
                                     ft, conf=True)
                dy, dx = 2.0 * dy, 2.0 * dx
                if (not np.isfinite([dy, dx]).all() or max(abs(dy), abs(dx)) > 400
                        or q < PEAK_MIN):
                    progress.log(f"[warn] {os.path.basename(p)}: no clear correlation "
                                 f"peak ({dy:.1f}, {dx:.1f} px, {q:.2f}x the next); "
                                 f"the previous tier's shift ({last_shift[0]:.1f}, "
                                 f"{last_shift[1]:.1f}) px used", None)
                    dy, dx = last_shift
                    n_fallback += 1
            shifts.append((float(dy), float(dx)))
            if q is not None and ref_feat is not ft:
                peaks.append(q)
            g = np.stack([_shift(f[c], dy, dx) for c in range(3)])
            if not prealigned:
                _sc = _shift(_clip.astype(np.float32), dy, dx)
                # > 0, not > 0.5: the bilinear shift spreads a clipped pixel's NaN
                # one pixel further than half its weight, and a ring of NaN round
                # the core reads as a hole in the alignment feature
                _sc = np.nan_to_num(_sc) > 0.0
                clipany = _sc if clipany is None else (clipany | _sc)
                del _sc
            del f, _zb, _clip
            ok = np.isfinite(g).all(axis=0).astype(np.float32)
            np.nan_to_num(g, copy=False)
            g *= ok[None]
            acc = g if acc is None else acc + g
            cnt = ok if cnt is None else cnt + ok
            del g, ok
        st = (acc / np.maximum(cnt, 1.0)).astype(np.float32)
        st[:, cnt < 0.5] = np.nan
        del acc, cnt
        last_shift = (float(np.median([q[0] for q in shifts])),
                      float(np.median([q[1] for q in shifts])))
        # the next tier's alignment reference: this stack, in the reference
        # frame's coordinates; pixels no frame held are no data (see _feature)
        if prealigned:
            tier_feat = None
        else:
            # where every frame of the tier clipped, the stack reads as the clip
            _g = st[1]
            if clipany is not None:
                _g = np.where(~np.isfinite(_g) & clipany, np.float32(sat - black[1]), _g)
            tier_feat = _feature(_g)
            del _g
        clipany = None

        if prev_g is not None:
            r, npx = _ratio(prev_g, st[1], sat)
            if not np.isfinite(r):
                r = s / secs[i - 1]
                progress.log(f"[warn] {secs[i-1]:.5g}s -> {s:.5g}s: too few shared "
                             f"pixels; header ratio {r:.4f} used", None)
            rel *= r
            ladder.append({"from": secs[i - 1], "to": s,
                           "header": s / secs[i - 1], "measured": r, "px": npx})
        prev_g = st[1].copy()

        if num is None:
            num = np.zeros_like(st)
            den = np.zeros_like(st)
        _pt = _sky_noise_point(st[1], len(fs))
        if _pt is not None:
            noise_models.append({"sec": s, "n": len(fs), "sky": _pt[0], "var_frame": _pt[1]})
        for c in range(3):
            v = st[c]
            okm = np.isfinite(v)
            vv = np.where(okm, v, 0.0).astype(np.float32)
            w = np.where(okm, _weight(vv, sat), np.float32(0.0))
            num[c] += w * vv / np.float32(rel)
            den[c] += w
            if c == 1:
                # the merge's own propagated noise (merge_noise.npy), as the
                # full pipeline writes it: var(merge) = sum w^2 var_t / rel^2,
                # var_t = (a v + b) / n with the camera's noise curve a, b,
                # fitted at the end from the tiers' sky -- so the two terms
                # are accumulated apart and combined then
                _q = w * w / np.float32(rel * rel * len(fs))
                if vnum is None:
                    vnum = [np.zeros_like(vv), np.zeros_like(vv)]
                vnum[0] += _q * vv
                vnum[1] += _q
                del _q

        # the four shortest tiers, in scene units: the inner-corona layer's
        # independent source, which an imported HDR does not have
        if i < 4:
            # a saturated pixel (NaN) used to enter as 0 and print a dark step
            # along each tier's saturation isophote in short_lum (audit A15):
            # it is left out and the mean taken over the tiers that hold it
            sc = st / np.float32(rel)
            lum_s = (0.2126 * sc[0] + 0.7152 * sc[1] + 0.0722 * sc[2])
            okp = np.isfinite(lum_s)
            lum_s = np.where(okp, lum_s, 0.0).astype(np.float32)
            short_acc = lum_s if short_acc is None else short_acc + lum_s
            short_n = okp.astype(np.float32) if short_n is None or isinstance(short_n, int) else short_n + okp
            del sc, lum_s, okp
        del st
        progress.log("  %10.5g s  %2d frames  shift rms %.2f px  max %.2f%s"
                     % (s, len(fs),
                        float(np.sqrt((np.array(shifts) ** 2).sum(1).mean())),
                        float(np.abs(np.array(shifts)).max()),
                        "  peak %.1fx" % min(peaks) if peaks else ""),
                     0.05 + 0.35 * (i + 1) / len(secs))

    if n_fallback:
        progress.log(f"  alignment: {n_fallback} frame(s) without a clear correlation "
                     f"peak kept the previous tier's shift (peak under {PEAK_MIN:g}x "
                     f"the next one)", None)

    hdr = np.where(den > 0, num / np.maximum(den, 1e-9), 0.0)
    relnoise = None
    _curve = _fit_noise_curve([(m_["sky"], m_["var_frame"]) for m_ in noise_models])
    if vnum is not None and _curve is not None:
        _va = np.float32(_curve[0]) * vnum[0] + np.float32(_curve[1]) * vnum[1]
        # 0 where the merge holds nothing, as the full pipeline writes it
        relnoise = np.where(num[1] > 0, np.sqrt(np.maximum(_va, 0)) / np.maximum(num[1], 1e-12), 0.0)
        relnoise = np.clip(np.nan_to_num(relnoise, nan=1.0, posinf=10.0), 0, 10).astype(np.float16)
        del _va
        progress.log(f"merge noise map: camera noise var = {_curve[0]:.3g} v + {_curve[1]:.3g} "
                     f"(from the sky of {len(noise_models)} tiers), propagated through the merge", None)
    elif vnum is not None:
        progress.log("[warn] merge noise map not built: no tier gave a sky noise point", None)
    del num, den, vnum
    rgb = np.ascontiguousarray(np.moveaxis(hdr, 0, -1), np.float32)
    del hdr
    short = ((short_acc / np.maximum(short_n, 1.0)).astype(np.float32)
             if short_acc is not None else None)

    progress.log("exposure ladder: header 1 .. %.4g, "
                 "measured 1 .. %.4g"
                 % (secs[-1] / secs[0], rel), None)
    for L in ladder:
        progress.log("  %10.5g -> %-10.5g header %8.4f   measured %8.4f"
                     % (L["from"], L["to"], L["header"], L["measured"]), None)

    # scale to the same comfortable working range an import gets; every layer
    # downstream is scale-free
    hi = float(np.percentile(rgb, 99.9))
    if hi > 0:
        rgb *= np.float32(50000.0 / hi)
        if short is not None:
            short *= np.float32(50000.0 / hi)

    return importhdr.build_from_rgb(
        folder, rgb, progress, denoise=denoise, short_lum=short,
        merge_noise=relnoise,
        fnrgf_preset=fnrgf_preset, partialconv=partialconv,
        stats={"n_files": n_used, "mode": "simple stack",
               "options": {"simple": True},
               "calibration": calib_info,
               "planes_input": is_planes,
               "simple_ladder": ladder,
               "simple_ladder_span": {"header": secs[-1] / secs[0],
                                      "measured": rel},
               "tier_mode": str(tier_mode or "exposure"),
               "simple_noise_points": noise_models,
               "simple_noise_curve": list(_curve) if _curve else None,
               "simple_align_fallbacks": n_fallback,
               "simple_prealigned": bool(prealigned),
               "align_mode": _am,
               "tier_measured": ({k: v for k, v in _ginfo.items() if k != "frames"}
                                 if _ginfo else None),
               "tiers": [{"sec": s, "n": len(tiers[s])} for s in secs]},
        opts={"mode": "simple", "denoise": denoise,
              "fnrgf_preset": fnrgf_preset,
              "partialconv": bool(partialconv),
              "demosaic": demosaic_method, "n_files": n_used,
              "flat_dir": (calib_info or {}).get("flat_dir", ""),
              "flat_applied": bool((calib_info or {}).get("flat_applied")),
              "bias_applied": bool((calib_info or {}).get("bias_applied")),
              "dark_applied": bool((calib_info or {}).get("dark_applied")),
              "wb_source": "none" if is_planes else "camera",
              "tier_mode": str(tier_mode or "exposure"),
              "simple_rev": SIMPLE_REV,
              "align_mode": _am,
              "prealigned": bool(prealigned),
              "secs": list(secs)})
