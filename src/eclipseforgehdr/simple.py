"""The simple stack: calibrate per frame, align, average per tier, merge on
measured exposure ratios. THE FALLBACK -- not the merge.

    read -> bias/dark, flat, demosaic, white balance, colour matrix  (per frame)
         -> one black level for the whole set, from the shortest tier
         -> align -> mean per tier
         -> exposure ratios MEASURED between adjacent tiers (a slope fit)
         -> hat-weighted linear merge
         -> importhdr.build_from_rgb: every enhancement layer, as an import

WHY IT IS THE FALLBACK AND NOT THE MERGE (0.23.5). It was promoted to the
only merge for one build, and on the 600 mm reference bracket -- 12 tiers at
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
        return direct
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
            progress.log(f"flat correction skipped ({e})", None)
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
            progress.log(f"bias/dark correction skipped ({e})", None)
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
                progress.log("  black level %.1f ADU from the shortest tier's "
                             "raw corners, subtracted before the flat (this "
                             "file reports none and there is no bias master)"
                             % blk, None)
        except Exception as e:
            progress.log(f"pre-flat black level not measured ({e})", None)
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
    """Half-resolution high-passed log green.

    HALF RESOLUTION ON PURPOSE: a full-res FFT pair on a 44 Mpx frame costs
    seconds each and there are tens of frames. The parabolic sub-pixel fit
    lands within about 0.1 of a half-res pixel, so 0.2 px at full scale --
    below what this path's own bilinear resampling costs anyway.
    """
    h, w = (g.shape[0] // 2) * 2, (g.shape[1] // 2) * 2
    g2 = g[:h, :w].reshape(h // 2, 2, w // 2, 2).mean((1, 3))
    return _hipass(np.log1p(np.clip(g2, 0, None) / 30.0), 24.0)


def _pcorr(a, b):
    """Sub-pixel shift taking b onto a. Plain phase correlation.

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
    return np.array(sub)


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


# ------------------------------------------------------------------ merge ---

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
    return (np.clip((v - 60.0) / 400.0, 0, 1)
            * np.clip((sat * 0.90 - v) / (sat * 0.55), 0, 1)).astype(np.float32)


# -------------------------------------------------------------------- run ---

def run(folder, progress, denoise="fine", demosaic_method="mhc",
        fnrgf_preset="ours", flat_dir="", partialconv=True):
    """Stack `folder` and build every layer the renderer needs."""
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
    for p in files:
        try:
            sec = float(read_exif(p)[0])
        except Exception:
            skipped += 1
            continue
        tiers.setdefault(sec, []).append(p)
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
    progress.log("SIMPLE STACK -- the fallback path", None)
    progress.log("3-plane FITS read as planes, used as written" if is_planes else
                 "Bayer raw: demosaic, camera white balance and colour matrix "
                 "per frame, as the full pipeline does", None)
    progress.log("no hot-pixel repair, no ladder solve, no LDIC, no feather; "
                 "the exposure ratios are measured from the pixels and the sky "
                 "is left in for the render to decide about", None)

    black = None
    ref_feat = None
    num = den = None
    prev_g = None
    sat = None
    rel = 1.0
    short_acc, short_n = None, 0
    ladder = []
    for i, s in enumerate(secs):
        fs = sorted(tiers[s])
        acc = cnt = None
        shifts = []
        for p in fs:
            f, sl = _load(p, demosaic_method, calib=calib, seconds=s)
            if sat is None:
                sat = sl
            if black is None:
                # ONE black level for the whole set, from the exposure where
                # the sky contributes least. See the module docstring.
                black = _pedestal(f)
                progress.log("  black level from the shortest tier's first "
                             "frame: R %.1f  G %.1f  B %.1f" % tuple(black), None)
            f -= black[:, None, None]
            # the alignment feature must be finite: a saturated pixel is NaN
            # by now, and reads as the saturation level for this purpose,
            # which is what the sensor recorded there
            ft = _feature(np.where(np.isfinite(f[1]), f[1], np.float32(sat)))
            if ref_feat is None:
                ref_feat = ft
                dy = dx = 0.0
            else:
                if ft.shape != ref_feat.shape:
                    raise RuntimeError(
                        "the frames in this folder are not all the same size")
                dy, dx = 2.0 * _pcorr(ref_feat, ft)
                if not np.isfinite([dy, dx]).all() or max(abs(dy), abs(dx)) > 400:
                    dy = dx = 0.0
            shifts.append((float(dy), float(dx)))
            g = np.stack([_shift(f[c], dy, dx) for c in range(3)])
            del f
            ok = np.isfinite(g).all(axis=0).astype(np.float32)
            np.nan_to_num(g, copy=False)
            g *= ok[None]
            acc = g if acc is None else acc + g
            cnt = ok if cnt is None else cnt + ok
            del g, ok
        st = (acc / np.maximum(cnt, 1.0)).astype(np.float32)
        st[:, cnt < 0.5] = np.nan
        del acc, cnt

        if prev_g is not None:
            r, npx = _ratio(prev_g, st[1], sat)
            if not np.isfinite(r):
                r = s / secs[i - 1]
                progress.log(f"  {secs[i-1]:.5g}s -> {s:.5g}s: too few shared "
                             f"pixels to measure the ratio; falling back to the "
                             f"headers' {r:.4f}", None)
            rel *= r
            ladder.append({"from": secs[i - 1], "to": s,
                           "header": s / secs[i - 1], "measured": r, "px": npx})
        prev_g = st[1].copy()

        if num is None:
            num = np.zeros_like(st)
            den = np.zeros_like(st)
        for c in range(3):
            v = st[c]
            okm = np.isfinite(v)
            vv = np.where(okm, v, 0.0).astype(np.float32)
            w = np.where(okm, _weight(vv, sat), np.float32(0.0))
            num[c] += w * vv / np.float32(rel)
            den[c] += w

        # the four shortest tiers, in scene units: the inner-corona layer's
        # independent source, which an imported HDR does not have
        if i < 4:
            sc = (np.nan_to_num(st) / np.float32(rel))
            lum_s = (0.2126 * sc[0] + 0.7152 * sc[1] + 0.0722 * sc[2])
            short_acc = lum_s if short_acc is None else short_acc + lum_s
            short_n += 1
            del sc, lum_s
        del st
        progress.log("  %10.5g s  %2d frames  shift rms %.2f px  max %.2f"
                     % (s, len(fs),
                        float(np.sqrt((np.array(shifts) ** 2).sum(1).mean())),
                        float(np.abs(np.array(shifts)).max())),
                     0.05 + 0.35 * (i + 1) / len(secs))

    hdr = np.where(den > 0, num / np.maximum(den, 1e-9), 0.0)
    del num, den
    rgb = np.ascontiguousarray(np.moveaxis(hdr, 0, -1), np.float32)
    del hdr
    short = (short_acc / max(short_n, 1)).astype(np.float32) if short_acc is not None else None

    progress.log("exposure ladder, measured from the pixels: header says "
                 "1 .. %.4g, the pixels say 1 .. %.4g"
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
        fnrgf_preset=fnrgf_preset, partialconv=partialconv,
        stats={"n_files": n_used, "mode": "simple stack",
               "options": {"simple": True},
               "calibration": calib_info,
               "planes_input": is_planes,
               "simple_ladder": ladder,
               "simple_ladder_span": {"header": secs[-1] / secs[0],
                                      "measured": rel},
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
              "secs": list(secs)})
