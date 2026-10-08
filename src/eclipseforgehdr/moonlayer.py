"""THE MOON LAYER (lab 0.24): earthshine structure on the lunar disc.

Built from the long frames of the bracket, after the stack, in about a minute,
and cached as MOON_FILE on the merge grid (0.5 = neutral, the format the old
earth.npy had, so render.py's Earthshine control drives it unchanged).

Recipe (after Project Helion's log on the same frames, re-derived and measured
on Nico's 600 mm set, Oct 2026):
  1. every long frame read from the raw file, calibrated with the cached master
     bias / dark / flat when they fit, superpixel green;
  2. each frame's Moon found by its own limb and the frames stacked ON THE MOON
     (the Moon moves against the Sun and the sensor during totality);
  3. the scattered light on the disc modelled as the merge's corona OUTSIDE the
     Moon convolved with a power-law PSF wing (exponent fitted), plus a smooth
     sky term, fitted on the disc and subtracted;
  4. anything that follows the limb taken off (the median of every ring): no
     lunar feature follows the Moon's own edge, the glare does;
  5. band-passed to the maria scales: the finest scale (a few px) is fixed to
     the SENSOR on this set (measured: early vs late frames correlate 0.41 in
     sensor coordinates, 0.23 on the Moon), so it is left out;
  6. faded out over the last 0.08 R, where the limb's spill owns the light.
Checked against Nico's own PixInsight earthshine of the same frames: same
orientation; correlation 0.54 / 0.67 / 0.61 / 0.27 in the bands 0.02-0.08 /
0.03-0.15 / 0.06-0.3 / 0.1-0.5 R (the largest is where the two scattered-light
models differ most).
"""
import os, json
import numpy as np
from scipy import ndimage
from scipy.signal import fftconvolve

MOON_FILE = "moon_layer.npy"
MOON_META = "moon_layer.json"
MOON_BUILD = 2              # 2: limb band levelled and equalised before the fade (2026-10-07)


def _log(progress, msg, frac=None):
    if progress is not None:
        progress.log(msg, frac)


def needs_build(wd):
    """Missing, built by an older recipe, or older than the merge it sits on."""
    try:
        mp = os.path.join(wd, MOON_FILE)
        if not os.path.exists(mp):
            return True
        if int(json.load(open(os.path.join(wd, MOON_META))).get("build", 0)) < MOON_BUILD:
            return True
        hp = os.path.join(wd, "hdr_rgb.npy")
        return os.path.exists(hp) and os.path.getmtime(mp) < os.path.getmtime(hp)
    except Exception:
        return True


def _fit_limb(G, cy, cx, R0, iters=3):
    """Circle through the steepest rise of log G along 360 rays near R0."""
    L = np.log(np.maximum(ndimage.gaussian_filter(G, 2.0), 1e-3))
    for _ in range(iters):
        pts = []
        rs = np.arange(0.93 * R0, 1.07 * R0, 0.25)
        for t in np.radians(np.arange(0, 360, 1.0)):
            v = ndimage.map_coordinates(L, [cy + rs * np.sin(t), cx + rs * np.cos(t)], order=1, mode="nearest")
            gr = np.gradient(v); i = int(np.argmax(gr))
            if 1 <= i < len(gr) - 1:
                a, b, c = gr[i - 1], gr[i], gr[i + 1]
                den = a - 2 * b + c
                d = 0.5 * (a - c) / den if den != 0 else 0.0
                rr = rs[i] + d * 0.25
                pts.append((cy + rr * np.sin(t), cx + rr * np.cos(t)))
        pts = np.array(pts)
        if len(pts) < 30:
            break
        for _ in range(4):
            A = np.c_[2 * pts[:, 1], 2 * pts[:, 0], np.ones(len(pts))]
            bb = pts[:, 1] ** 2 + pts[:, 0] ** 2
            s = np.linalg.lstsq(A, bb, rcond=None)[0]
            xc, yc = s[0], s[1]; R = float(np.sqrt(s[2] + xc * xc + yc * yc))
            res = np.hypot(pts[:, 1] - xc, pts[:, 0] - yc) - R
            k = np.abs(res) < max(3 * 1.4826 * np.median(np.abs(res)), 0.3)
            pts = pts[k]
        cy, cx, R0 = float(yc), float(xc), R
    return cy, cx, R0


def _masters(wd, shape):
    """Cached master bias / dark rate / flat, each only if it fits the frame."""
    bias = rate = flat = None
    try:
        z = np.load(os.path.join(wd, "masterdark.npz"))
        if "bias" in z.files and z["bias"].shape == shape:
            bias = z["bias"]
        if "rate" in z.files and z["rate"].shape == shape:
            rate = z["rate"]
    except Exception:
        pass
    try:
        f = np.load(os.path.join(wd, "masterflat.npy"), mmap_mode="r")
        if f.shape == shape:
            flat = f
    except Exception:
        pass
    return bias, rate, flat


def build(folder, wd, progress=None, frames=None):
    """Build MOON_FILE / MOON_META. Returns the meta dict, or None."""
    from .raw import open_frame, list_raws, read_exif
    geo = json.load(open(os.path.join(wd, "geometry.json")))
    cy, cx, R = float(geo["cy"]), float(geo["cx"]), float(geo["R"])
    oy, ox = (geo.get("crop_origin") or [0, 0])[:2]
    # the GREEN channel of the merge: the stack below is green, and luminance
    # carries the prominences' H-alpha red, whose scatter the green frames do
    # not see -- with luminance the model over-subtracted under the big
    # prominence (a dark patch at the limb)
    _hr = np.load(os.path.join(wd, "hdr_rgb.npy"), mmap_mode="r")
    H, W = _hr.shape[:2]
    # ---- the frames: the long end of the bracket, where the disc has signal
    if frames is None:
        frames = []
        for p in list_raws(folder):
            try:
                frames.append((p, float(read_exif(p)[0])))
            except Exception:
                continue
    if not frames:
        return None
    tmax = max(t for _, t in frames)
    use = [(p, t) for p, t in frames if t >= tmax / 16.0 - 1e-9]
    # THE LONGEST EXPOSURE IS LEFT OUT when at least two shorter ones remain.
    # Measured on the 600 mm set against Nico's own earthshine: its large-scale
    # pattern is instrumental (4x the others' level in rate units), and with it
    # the maria scale agrees 0.05 with his result, without it 0.27 (0.61 at
    # 0.04-0.2 R instead of 0.53).
    _shorter = sorted(set(t for _, t in use if t < tmax - 1e-9))
    if len(_shorter) >= 2:
        use = [(p, t) for p, t in use if t < tmax - 1e-9]
    if not use:
        return None
    _log(progress, f"Moon layer: stacking {len(use)} long frame(s) "
                   f"({min(t for _, t in use):g}-{max(t for _, t in use):g} s"
                   + (f"; longest {tmax:g} s left out" if len(_shorter) >= 2 else "")
                   + ")...", None)
    # ---- the window: a box about the Moon in frame coordinates, superpixel grid
    Rh = R / 2.0
    half = int(np.ceil(1.35 * R / 2.0)) * 2
    acc = None; wsum = None; ref = None; nused = 0
    for p, t in use:
        try:
            rf = open_frame(p)
            bay = rf.bayer; Hb, Wb = bay.shape
            sat = float(rf.sat_level)
            gy, gx = cy + oy, cx + ox                # the merge's Moon, in frame px
            Y0 = int(gy - half) - int(gy - half) % 2; X0 = int(gx - half) - int(gx - half) % 2
            Y1, X1 = Y0 + 2 * half, X0 + 2 * half
            cfa = np.zeros((2 * half, 2 * half), np.float32)
            yA, yB, xA, xB = max(Y0, 0), min(Y1, Hb), max(X0, 0), min(X1, Wb)
            cfa[yA - Y0:yB - Y0, xA - X0:xB - X0] = bay[yA:yB, xA:xB]
            clip = cfa >= 0.9 * sat
            bias, rate, flat = _masters(wd, (Hb, Wb))
            if bias is not None:
                cfa[yA - Y0:yB - Y0, xA - X0:xB - X0] -= bias[yA:yB, xA:xB]
            if rate is not None:
                cfa[yA - Y0:yB - Y0, xA - X0:xB - X0] -= rate[yA:yB, xA:xB] * np.float32(t)
            if flat is not None:
                cfa[yA - Y0:yB - Y0, xA - X0:xB - X0] /= np.maximum(np.asarray(flat[yA:yB, xA:xB], np.float32), 1e-3)
            del rf, bay
        except Exception as e:
            _log(progress, f"[warn] Moon layer: {os.path.basename(p)} skipped ({e})", None)
            continue
        G = 0.5 * (cfa[0::2, 1::2] + cfa[1::2, 0::2])
        cl = clip.reshape(half, 2, half, 2).any(axis=(1, 3))
        c0 = (gy - Y0) / 2.0 - 0.25; c1 = (gx - X0) / 2.0 - 0.25
        try:
            fy, fx, fR = _fit_limb(G, c0, c1, Rh)
        except Exception:
            continue
        if not (0.9 * Rh < fR < 1.1 * Rh) or np.hypot(fy - c0, fx - c1) > 0.2 * Rh:
            _log(progress, f"[warn] Moon layer: {os.path.basename(p)} skipped, limb not found", None)
            continue
        if ref is None:
            ref = (fy, fx, fR, Y0, X0)
        dy, dx = ref[0] - fy + (Y0 - ref[3]) / 2.0, ref[1] - fx + (X0 - ref[4]) / 2.0
        Gs = ndimage.shift(G, (dy, dx), order=3, mode="nearest") / np.float32(t)
        ok = ndimage.shift(cl.astype(np.float32), (dy, dx), order=1, mode="constant", cval=1.0) < 0.01
        w = np.float32(t) * ok
        acc = Gs * w if acc is None else acc + Gs * w
        wsum = w if wsum is None else wsum + w
        nused += 1
    if acc is None or nused == 0:
        _log(progress, "[warn] Moon layer: no usable frame", None)
        return None
    M = acc / np.maximum(wsum, 1e-9)
    mcy, mcx, mR = ref[0], ref[1], ref[2]
    n = M.shape[0]
    yy, xx = np.mgrid[:n, :n].astype(np.float64)
    r = np.hypot(yy - mcy, xx - mcx) / mR
    # ---- the light outside the Moon, from the merge, on this grid (the grid's
    # Moon placed on the merge's Moon)
    my = cy + (yy - mcy) * 2.0
    mx = cx + (xx - mcx) * 2.0
    _y0, _y1 = max(int(my.min()) - 2, 0), min(int(my.max()) + 3, H)
    _x0, _x1 = max(int(mx.min()) - 2, 0), min(int(mx.max()) + 3, W)
    Gm = np.asarray(_hr[_y0:_y1, _x0:_x1, 1], np.float32)
    I = ndimage.map_coordinates(Gm, [my - _y0, mx - _x0], order=1, mode="constant", cval=0.0)
    del Gm
    I[r < 1.005] = 0.0
    inside = r < 0.97
    sub = np.flatnonzero(inside.ravel())[::2]
    X = (xx - mcx) / mR; Yv = (yy - mcy) / mR
    ky, kx = np.mgrid[-n:n + 1, -n:n + 1].astype(np.float64)
    rho2 = ky * ky + kx * kx
    best = None
    for alpha in (2.0, 2.5, 3.0):
        S = fftconvolve(I, (rho2 + 4.0) ** (-alpha / 2), mode="same")
        A = np.stack([S.ravel()[sub], np.ones(len(sub)), X.ravel()[sub], Yv.ravel()[sub],
                      (X * X + Yv * Yv).ravel()[sub]], -1)
        y = M.ravel()[sub]; keep = np.ones(len(y), bool)
        for _ in range(3):
            sc = np.abs(A[keep]).max(0) + 1e-30
            coef = np.linalg.lstsq(A[keep] / sc, y[keep], rcond=None)[0] / sc
            res = y - A @ coef; s = 1.4826 * np.median(np.abs(res[keep])); keep = np.abs(res) < 3 * s
        if best is None or s < best[0]:
            best = (s, alpha, coef, S)
    s, alpha, coef, S = best
    model = coef[0] * S + coef[1] + coef[2] * X + coef[3] * Yv + coef[4] * (X * X + Yv * Yv)
    E = np.where(inside, M - model, 0.0)
    level = float(np.median(model[r < 0.5]))
    # ---- 4: what follows the limb: the median of every ring
    rb = np.clip((r * mR).astype(np.int32), 0, int(mR) + 2)
    msk = inside.ravel()
    med = np.zeros(rb.max() + 1)
    order = np.argsort(rb.ravel()[msk], kind="stable")
    rv = rb.ravel()[msk][order]; ev = E.ravel()[msk][order]
    cuts = np.searchsorted(rv, np.arange(rb.max() + 2))
    for k in range(rb.max() + 1):
        seg = ev[cuts[k]:cuts[k + 1]]
        if seg.size >= 8:
            med[k] = np.median(seg)
    med = ndimage.gaussian_filter1d(med, 2.0)
    E = np.where(inside, E - med[rb], 0.0)
    # ---- 4b: glare that hugs the limb in some directions only (a bright
    # streamer's spill): per direction, the amplitude of a profile rising
    # toward the limb, b(r) = exp((r - 0.97) / 0.06), fitted on 0.6-0.95 R,
    # smoothed over 10 degrees and taken off
    th = np.arctan2(yy - mcy, xx - mcx)
    bprof = np.exp((r - 0.97) / 0.06)
    zone = inside & (r > 0.6) & (r < 0.95)
    nb = 360
    tb = ((th[zone] + np.pi) / (2 * np.pi) * nb).astype(np.int32) % nb
    ez, bz = E[zone], bprof[zone]
    num = np.bincount(tb, ez * bz, nb); den = np.bincount(tb, bz * bz, nb)
    amp = np.where(den > 0, num / np.maximum(den, 1e-12), 0.0)
    amp = ndimage.gaussian_filter1d(amp, 10.0, mode="wrap")
    ti = ((th + np.pi) / (2 * np.pi) * nb) % nb
    a_map = ndimage.map_coordinates(amp, [ti.ravel()], order=1, mode="grid-wrap").reshape(ti.shape)
    E = np.where(inside, E - a_map * bprof, 0.0)
    # ---- 5: band-pass by normalised convolution inside the disc: above the
    # sensor's fine scale, below 0.3 R
    w = inside.astype(np.float64)
    def nconv(a, sg):
        return ndimage.gaussian_filter(a * w, sg) / np.maximum(ndimage.gaussian_filter(w, sg), 1e-6)
    s_lo = max(3.0, 0.012 * mR); s_hi = 0.3 * mR
    Eb = np.where(inside, nconv(E, s_lo) - nconv(E, s_hi), 0.0)
    # ---- 5b: THE LIMB BAND (build 2, 2026-10-07). Measured on Clifton's
    # 250 mm set and on Nico's 600 mm frames: after the band-pass the structure
    # grew toward the edge (ring std 1.8x the centre's at 0.75-0.86 R) and the
    # rings there sat dark (layer mean 0.44 against 0.50), then both dropped to
    # flat within 0.05 R -- a visible rim. Real maria lose contrast toward the
    # limb (foreshortening), so that growth is residual glare and the
    # normalised convolution's edge bias. Per ring: the median taken off, and
    # the spread held to the inner disc's.
    ring_med = np.zeros(rb.max() + 1); ring_sd = np.zeros(rb.max() + 1)
    ebv = Eb.ravel()[msk][order]
    for k in range(rb.max() + 1):
        seg = ebv[cuts[k]:cuts[k + 1]]
        if seg.size >= 8:
            ring_med[k] = np.median(seg)
            ring_sd[k] = 1.4826 * np.median(np.abs(seg - ring_med[k]))
    ring_med = ndimage.gaussian_filter1d(ring_med, 3.0)
    ring_sd = ndimage.gaussian_filter1d(ring_sd, 3.0)
    _rk = np.arange(rb.max() + 1) / mR
    _ref = ring_sd[(_rk > 0.2) & (_rk < 0.6)]
    sd0 = float(np.median(_ref)) if _ref.size else 0.0
    gain = np.where(ring_sd > 1e-12, np.minimum(1.0, sd0 / np.maximum(ring_sd, 1e-12)), 1.0) if sd0 > 0 else np.ones_like(ring_sd)
    Eb = np.where(inside, (Eb - ring_med[rb]) * gain[rb], 0.0)
    # ---- 6: fade out where the limb's spill owns the light; wider than the
    # 0.87-0.95 R of build 1 so the structure tapers instead of stopping
    fade = np.clip((0.97 - r) / 0.17, 0, 1); fade = fade * fade * (3 - 2 * fade)
    Eb *= fade
    sd = float(1.4826 * np.median(np.abs(Eb[r < 0.8] - np.median(Eb[r < 0.8])))) or 1.0
    e = 0.5 + 0.5 * np.tanh(Eb / (2.5 * sd))
    # ---- onto the merge grid
    Yg = np.arange(H, dtype=np.float64)
    out = np.full((H, W), 0.5, np.float16)
    y0m, y1m = max(int(cy - 1.05 * R), 0), min(int(cy + 1.05 * R) + 1, H)
    x0m, x1m = max(int(cx - 1.05 * R), 0), min(int(cx + 1.05 * R) + 1, W)
    gyy, gxx = np.mgrid[y0m:y1m, x0m:x1m].astype(np.float64)
    py = mcy + (gyy - cy) / 2.0; px = mcx + (gxx - cx) / 2.0
    out[y0m:y1m, x0m:x1m] = ndimage.map_coordinates(e, [py, px], order=1, mode="constant", cval=0.5).astype(np.float16)
    np.save(os.path.join(wd, MOON_FILE), out)
    meta = {"build": MOON_BUILD, "frames": nused, "exposures": sorted(set(round(t, 6) for _, t in use)),
            "psf_alpha": alpha, "disc_level": level, "structure_rms": sd,
            "structure_rel": sd / max(level, 1e-9), "band_px_halfres": [s_lo, s_hi]}
    json.dump(meta, open(os.path.join(wd, MOON_META), "w"), indent=1)
    _log(progress, f"[ok] Moon layer: {nused} frames, scatter model r^-{alpha:g}, "
                   f"maria {100 * sd / max(level, 1e-9):.2f}% of disc light", None)
    return meta
