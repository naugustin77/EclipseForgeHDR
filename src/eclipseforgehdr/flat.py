"""Flat-field calibration.

A flat divides out everything multiplicative between the sky and the number in
the file: lens vignetting, the cos^4 falloff, dust shadows on the sensor stack,
filter non-uniformity, and per-photosite sensitivity (PRNU).  For an eclipse it
matters more than for most subjects, because the corona's own radial falloff is
the signal: a 6% vignette across the field is a 6% error in the F-corona
gradient, and every radial filter downstream (MGN, FNRGF, NAFE) then works to
preserve it.

    corrected = (light - black) / F

with F normalised to 1.0 at the centre, so the correction leaves the middle of
the frame -- where the disc, the chromosphere and the calibration all live --
exactly where it was, and only lifts the outer field.

THE NOISE PROBLEM, AND WHY THIS MODULE MEASURES ITS OWN
-------------------------------------------------------
A master flat is not free.  Dividing by it multiplies the light frame by
1/F, and whatever noise F carries is injected into every pixel of every frame,
identically, so stacking cannot average it away.  A flat exposed at 12% of full
well and averaged over 20 frames carries ~1.2% noise per photosite -- which,
applied to correct a 6% vignette, would put more noise in than gradient out.

So the master is smoothed, and the smoothing scale is not a taste setting: the
frames are split into two independent half-stacks, the same processing is
applied to both, and the residual between them IS the noise of the master.  The
Gaussian sigma is chosen as the smallest that brings that measured number under
NOISE_TARGET, then verified by measuring again.  On a clean, well-exposed flat
set the answer is 0 or 1 px and every dust mote survives at full resolution; on
a thin, noisy set it grows until only the vignetting is left.  Either way the
number that ends up in the log is measured on the user's own data, not assumed.

Robustness: frames are normalised per Bayer channel to their own central
median, so exposure and illumination-colour drift between flats cancels and the
master cannot shift the white balance of the lights.  The per-pixel combine is
a min/max-trimmed mean, which removes a cosmic ray, a satellite or a bird in
any one frame.
"""
from __future__ import annotations
import os, json
import numpy as np

# Per-photosite noise the master flat is allowed to inject into every light
# frame, as a fraction.  0.2% is an order of magnitude below the fine structure
# the detail filters work on, and roughly a tenth of a typical vignette.
NOISE_TARGET = 0.002
# Never smooth beyond this fraction of the frame diagonal.  At that point the
# flat is a vignetting model and nothing is gained by going further; it also
# stops a hopeless flat set from turning into a blank field.
MAX_SIGMA_FRAC = 0.02
# Flats outside this window of the saturation level are not usable: too dark is
# noise, too bright is the shoulder of the sensor's response.
MIN_LEVEL_FRAC = 0.02
MAX_LEVEL_FRAC = 0.85

_DIR_NAMES = ("flats", "flat", "flatfield", "flat_field", "flatfields",
              "flat-field", "masterflat")


def find_flat_dir(folder):
    """A conventionally-named flats subfolder of `folder`, or None.

    Case-insensitive, one level down.  Note that list_raws() reads only the
    files directly in a folder, so a flats subfolder can never be mistaken for
    light frames.
    """
    try:
        entries = sorted(os.listdir(folder))
    except OSError:
        return None
    for name in entries:
        p = os.path.join(folder, name)
        if os.path.isdir(p) and name.lower() in _DIR_NAMES:
            return p
    return None


def fingerprint(flat_dir):
    """Identity of the flat files, for cache validation."""
    from .raw import list_raws
    out = []
    if not flat_dir or not os.path.isdir(flat_dir):
        return out
    for p in sorted(list_raws(flat_dir)):
        try:
            st = os.stat(p)
            out.append([os.path.basename(p), int(st.st_size), int(st.st_mtime)])
        except OSError:
            out.append([os.path.basename(p), -1, -1])
    return out


def _central_box(shape, frac=0.15):
    """Even-aligned central box, so the Bayer phase inside it is the frame's."""
    H, W = shape
    hy = max(32, int(0.5 * frac * H)) & ~1
    hx = max(32, int(0.5 * frac * W)) & ~1
    cy = (H // 2) & ~1
    cx = (W // 2) & ~1
    y0, y1 = max(cy - hy, 0), min(cy + hy, H)
    x0, x1 = max(cx - hx, 0), min(cx + hx, W)
    return slice(y0 - y0 % 2, y1), slice(x0 - x0 % 2, x1)


def _channel_medians(a, box=None):
    sub = a[box] if box is not None else a
    return [float(np.median(sub[oy::2, ox::2])) for oy in (0, 1) for ox in (0, 1)]


def _normalise_channels(a, box):
    """Divide each Bayer channel by its own median inside `box`. In place."""
    ok = True
    for i, (oy, ox) in enumerate(((0, 0), (0, 1), (1, 0), (1, 1))):
        m = float(np.median(a[box][oy::2, ox::2]))
        if not np.isfinite(m) or m <= 0:
            ok = False
            break
        a[oy::2, ox::2] /= np.float32(m)
    return ok


def _chain(a, sigma_full, median3=True):
    """Denoise a CFA-sampled flat: per-channel 3x3 median, then Gaussian.

    Both run on each Bayer sub-lattice separately -- the four channels sit at
    different levels and have different shapes, and mixing them would smear the
    colour structure of the flat into a checkerboard.  One sub-lattice pixel is
    two full-res pixels, hence sigma/2.
    """
    from scipy import ndimage
    out = np.empty_like(a)
    s = 0.5 * float(sigma_full)
    for oy in (0, 1):
        for ox in (0, 1):
            c = np.asarray(a[oy::2, ox::2], np.float32)
            if median3:
                c = ndimage.median_filter(c, size=3, mode="nearest")
            if s > 0.05:
                c = ndimage.gaussian_filter(c, s, mode="nearest")
            out[oy::2, ox::2] = c
    return out


def _robust_sd(d):
    m = float(np.median(d))
    return float(1.4826 * np.median(np.abs(d - m)))


def _centre_crop(shape, half):
    """Even-aligned central crop of half-size `half`, so the Bayer phase holds."""
    H, W = shape
    hy = int(min(half, H // 2 - 2)) & ~1
    hx = int(min(half, W // 2 - 2)) & ~1
    cy, cx = (H // 2) & ~1, (W // 2) & ~1
    return slice(cy - hy, cy + hy), slice(cx - hx, cx + hx)


def _measure_noise(mA, mB, sigma, median3=True):
    """Noise of the smoothed master, from the two half-stacks.

    Only a central patch is filtered: the answer is a per-pixel statistic and
    a few megapixels give it to three digits, where filtering the whole frame
    for every trial sigma would cost minutes on a 45 Mpx sensor. The patch is
    grown with sigma and the sd is taken well inside it, so the edge handling
    of the filter never enters the number.
    """
    half = max(512.0, 6.0 * sigma)
    c = _centre_crop(mA.shape, half)
    a = _chain(mA[c].copy(), sigma, median3)
    b = _chain(mB[c].copy(), sigma, median3)
    d = (a - b) * np.float32(0.5)
    k = int(min(3.0 * sigma + 4, 0.4 * min(d.shape))) & ~1
    if k > 0:
        d = d[k:-k, k:-k]
    return _robust_sd(d)


def build_master(flat_dir, shape=None, progress=None, target=NOISE_TARGET):
    """Master flat from every usable frame in `flat_dir`.

    Returns (master, info). `master` is float32 in the same CFA layout and
    shape as a decoded light frame, normalised to 1.0 at the centre of each
    Bayer channel, or None if no usable flat could be built (info["error"]
    then says why, in words meant for the run log).
    """
    from .raw import list_raws, open_frame
    log = (progress.log if progress is not None else (lambda *a, **k: None))
    info = {"dir": flat_dir}
    paths = list_raws(flat_dir) if flat_dir and os.path.isdir(flat_dir) else []
    info["n_found"] = len(paths)
    if not paths:
        info["error"] = f"no readable frames in {flat_dir}"
        return None, info

    sumA = sumB = mn = mx = None
    nA = nB = 0
    used, rejected = [], []
    box = None
    for i, p in enumerate(paths):
        name = os.path.basename(p)
        try:
            rf = open_frame(p)
            b = np.asarray(rf.bayer, np.float32)
            sat = float(rf.sat_level)
            del rf
        except Exception as e:
            rejected.append((name, f"could not be decoded ({e})"))
            continue
        if shape is not None and tuple(b.shape) != tuple(shape):
            rejected.append((name, f"is {b.shape[1]}x{b.shape[0]}, the light "
                                   f"frames are {shape[1]}x{shape[0]}"))
            del b
            continue
        if box is None:
            box = _central_box(b.shape)
            if shape is None:
                shape = b.shape
        meds = _channel_medians(b, box)
        hi = max(meds)
        if not np.isfinite(hi) or hi <= 0:
            rejected.append((name, "has no signal"))
            del b
            continue
        if hi > MAX_LEVEL_FRAC * sat:
            rejected.append((name, f"is exposed to {100 * hi / sat:.0f}% of "
                                   f"saturation — on the shoulder of the "
                                   f"sensor's response, not linear"))
            del b
            continue
        if hi < MIN_LEVEL_FRAC * sat:
            rejected.append((name, f"is exposed to only {100 * hi / sat:.1f}% "
                                   f"of saturation — too dark to calibrate with"))
            del b
            continue
        if not _normalise_channels(b, box):
            rejected.append((name, "has a zero or negative channel median"))
            del b
            continue
        np.clip(b, 0.0, 8.0, out=b)
        if sumA is None:
            sumA = np.zeros_like(b)
            sumB = np.zeros_like(b)
            mn = b.copy()
            mx = b.copy()
        else:
            np.minimum(mn, b, out=mn)
            np.maximum(mx, b, out=mx)
        # Split by parity, not by half: two independent half-stacks that see
        # the same drift in the illumination, so their difference is noise and
        # not a trend.
        if len(used) % 2 == 0:
            sumA += b
            nA += 1
        else:
            sumB += b
            nB += 1
        used.append(name)
        log(f"flat {i + 1}/{len(paths)}: {name}", None)
        del b

    for name, why in rejected:
        log(f"flat rejected: {name} {why}", None)
    info["n_used"] = len(used)
    info["rejected"] = [{"file": n, "why": w} for n, w in rejected]
    n = len(used)
    if n < 2:
        info["error"] = (f"only {n} usable flat frame(s) in {flat_dir} — "
                         f"at least 2 are needed")
        return None, info

    # Built in place, to keep the peak working set at four frame-sized arrays
    # rather than seven -- on a 45 Mpx sensor that is the difference between
    # 0.8 and 1.4 GB.
    master = sumA + sumB
    if n >= 5:
        mn += mx
        master -= mn
        master /= np.float32(n - 2)
        info["combine"] = f"min/max-trimmed mean of {n} frames"
    else:
        master /= np.float32(n)
        info["combine"] = f"mean of {n} frames"
    del mn, mx

    # --- how noisy is it, and how much smoothing does that buy ---
    #
    # var(meanA) = var(meanB) = 2 * var(master), so var((A-B)/2) = var(master):
    # the half-stack difference is not a proxy for the master's noise, it is
    # exactly it -- and it stays exact through any linear or near-linear
    # processing applied identically to both.
    sigma = 0.0
    med3 = False
    if nA > 0 and nB > 0:
        sumA /= np.float32(nA)          # now meanA
        sumB /= np.float32(nB)          # now meanB
        noise0 = _robust_sd(((sumA - sumB) * np.float32(0.5))[box])
        # The half-stacks are plain means; the master is a trimmed one, which
        # averages n-2 samples instead of n. Small (5% at n=20) but it is in
        # the direction of understating the noise, so it is put back.
        trim = float(np.sqrt(n / float(n - 2))) if n >= 5 else 1.0
        noise0 *= trim
        info["noise_raw"] = float(noise0)
        cap = MAX_SIGMA_FRAC * float(np.hypot(*master.shape))
        info["noise_master"] = float(noise0)
        info["sigma_px"] = 0.0
        # A flat set good enough on its own is left completely alone: no
        # smoothing and no median either, so real per-photosite sensitivity
        # (PRNU) is corrected rather than thrown away with the noise.
        if noise0 > target:
            med3 = True
            # White noise through a normalised 2-D Gaussian on the sub-lattice:
            # sd_out = sd_in / (2*sigma_sub*sqrt(pi)); the 3x3 median ahead of
            # it contributes about 0.55.  Used as a starting guess only -- the
            # answer is then measured, and raised until it is met.
            sigma = 2.0 * (0.55 * noise0) / (2.0 * np.sqrt(np.pi) * target)
            sigma = float(min(max(sigma, 0.0), cap))
            # The starting guess assumes white noise. Real sensor noise is not
            # quite white -- on the reference set it needed 2.9x more smoothing
            # than the formula predicted -- so the guess is only a starting
            # point and sigma is raised until the MEASURED noise meets it.
            for _ in range(6):
                nz = trim * _measure_noise(sumA, sumB, sigma, med3)
                info["noise_master"] = float(nz)
                info["sigma_px"] = float(sigma)
                if nz <= target or sigma >= cap:
                    break
                sigma = float(min(max(sigma * 1.7, 1.0), cap))
    info["median3"] = bool(med3)
    del sumA, sumB

    if med3 or sigma > 0:
        master = _chain(master, sigma, med3)

    # Renormalise after smoothing so the centre is exactly 1.0 per channel: the
    # correction must not change the white balance or the overall level of the
    # lights, only their spatial structure.
    _normalise_channels(master, box)
    master[~np.isfinite(master)] = 1.0
    lo, hi = float(master.min()), float(master.max())
    nclip = int((master < 0.2).sum() + (master > 5.0).sum())
    if nclip:
        np.clip(master, 0.2, 5.0, out=master)
    info["clipped"] = nclip
    info["min"] = lo
    info["max"] = hi
    # What the correction actually does, as one number: the ratio between the
    # dimmest and brightest part of the field once the noise is out of it.
    # Percentiles, not min/max, and a 1% margin dropped: a raw decoder that
    # hands back a few masked rows at the frame edge would otherwise report a
    # 30% vignette on a lens that has 8%.
    from scipy import ndimage
    H, W = master.shape
    h, w = H // 2, W // 2
    sp = master[:2 * h, :2 * w].reshape(h, 2, w, 2).mean(axis=(1, 3))
    sm = ndimage.gaussian_filter(sp, max(h, w) / 60.0, mode="nearest")
    my, mx_ = max(int(0.01 * h), 1), max(int(0.01 * w), 1)
    inner = sm[my:h - my, mx_:w - mx_]
    lo = float(np.percentile(inner, 0.5))
    hi_ = float(np.percentile(inner, 99.5))
    info["vignette"] = float(hi_ / max(lo, 1e-6))
    info["corner"] = float(lo / max(hi_, 1e-6))
    return master.astype(np.float32), info


def describe(info):
    """One-line-per-fact summary of a build, for the run log and the report."""
    out = []
    if info.get("error"):
        out.append(f"flat: {info['error']}")
        return out
    out.append(f"flat: {info.get('combine', '')} from "
               f"{os.path.basename(info.get('dir') or '')}/")
    if "noise_raw" in info:
        if not info.get("median3") and not info.get("sigma_px"):
            out.append(f"flat: per-pixel noise {100 * info['noise_raw']:.3f}% "
                       f"— already under the {100 * NOISE_TARGET:.1f}% target, "
                       f"so it is used at full resolution, unsmoothed")
        else:
            out.append(f"flat: per-pixel noise {100 * info['noise_raw']:.3f}% "
                       f"-> {100 * info.get('noise_master', float('nan')):.3f}% "
                       f"after a 3x3 median and a "
                       f"{info.get('sigma_px', 0):.1f} px smooth "
                       f"(target {100 * NOISE_TARGET:.1f}%)")
    out.append(f"flat: corrects a {100 * (info.get('vignette', 1) - 1):.1f}% "
               f"falloff — the dimmest part of the field sits at "
               f"{info.get('corner', 1):.3f} of the brightest")
    if info.get("clipped"):
        out.append(f"flat: {info['clipped']} photosite(s) clipped to the "
                   f"[0.2, 5.0] safety range")
    return out


def superpixel_full(master):
    """Local mean of the four Bayer channels, at full resolution.

    The saturation tests downstream compare demosaiced values against a scalar
    saturation level, and those values are flat-corrected while the level is
    not.  Multiplying by this puts them back in raw units for the comparison,
    which matters at the corners of a strongly vignetted frame: without it a
    long tier's sky reads as clipped where it is merely dim.
    """
    H, W = master.shape
    h, w = H // 2, W // 2
    m = master[:2 * h, :2 * w].reshape(h, 2, w, 2).mean(axis=(1, 3))
    out = np.empty((H, W), np.float32)
    out[:2 * h, :2 * w] = np.repeat(np.repeat(m, 2, 0), 2, 1)
    if H > 2 * h:
        out[2 * h:, :2 * w] = out[2 * h - 1, :2 * w]
    if W > 2 * w:
        out[:, 2 * w:] = out[:, 2 * w - 1:2 * w]
    return out


def load_or_build(folder, flat_dir, shape, progress=None, workdir=None):
    """Cached master flat for a light folder.

    Rebuilds whenever the flat files, the frame shape or this module's
    parameters change; otherwise loads the cached array, which for a 45 Mpx
    sensor saves a couple of minutes on every re-run.
    """
    log = (progress.log if progress is not None else (lambda *a, **k: None))
    fp = fingerprint(flat_dir)
    key = {"files": fp, "shape": list(shape) if shape is not None else None,
           "target": NOISE_TARGET, "fmt": 2}
    npy = os.path.join(workdir, "masterflat.npy") if workdir else None
    js = os.path.join(workdir, "masterflat.json") if workdir else None
    if npy and os.path.exists(npy) and js and os.path.exists(js):
        try:
            cached = json.load(open(js))
            if cached.get("key") == key:
                m = np.load(npy)
                if shape is None or tuple(m.shape) == tuple(shape):
                    info = cached.get("info", {})
                    info["cached"] = True
                    for line in describe(info):
                        log(line + " (cached)", None)
                    return m, info
        except Exception:
            pass
    log(f"building the master flat from {len(fp)} frame(s)...", None)
    master, info = build_master(flat_dir, shape, progress)
    for line in describe(info):
        log(line, None)
    if master is not None and npy:
        try:
            np.save(npy, master)
            json.dump({"key": key, "info": info}, open(js, "w"), indent=1)
        except Exception as e:
            log(f"flat: could not cache the master ({e})", None)
    return master, info
