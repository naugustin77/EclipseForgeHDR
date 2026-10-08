"""Prominence layer: the prominences as a separate RGBA layer for post-processing.

WHY A SEPARATE LAYER. Every attempt to merge the sharp prominence stack into
the composite itself failed on the picture: the renderer normalises and tints
the whole frame, so a deep-red, very bright prominence core either clips or
turns into a flat red blob. Nico's proposal was to keep the prominences out of
the composite's hands and hand them over as their own layer -- same pixel
grid as the 16-bit export, alpha already masked -- to be blended in
Photoshop with Blend If / opacity / screen. This module builds that layer.

TWO SOURCES, EACH FOR WHAT IT DOES WELL (measured on Nico's 600 mm set):
  * the HDR MERGE's red excess has 2-3x less noise than the short stack, but
    is smeared (it holds only 0.33-0.58 of the short stack's peak). It decides
    WHERE a prominence is: the detection gate and the faint parts.
  * the SHORT-EXPOSURE STACK (fast tiers, aligned on the red excess and
    Richardson-Lucy deconvolved, built by build_stack() during the run) is
    sharp. It decides the EDGE of each prominence and carries its bright core.
A cache without the short stack still gets a layer, from the merge alone.

THE CHROMOSPHERE BAND. Along parts of the limb the chromosphere shows as a
continuous pink band 25-50 px deep (5-12 px elsewhere on the same set), and
there the small prominences are as tall as the band, so a height cut takes
them too. The band is instead removed as what it is -- light that is smooth
ALONG the limb: in polar coordinates, subtract the 30th percentile over
+-2.5 deg. Compact prominences stand out of that; prominences that rise well
above the band keep their whole profile down to the limb.

THE EDGE. A fixed-sigma threshold on a bright small prominence lands far out
in its point-spread wings: the merge-only mask was 2.5-5x the area of the
prominence at 25% of its own peak. The edge is therefore taken RELATIVE TO
EACH PROMINENCE'S OWN PEAK in the sharp stack (0.12-0.35 of the local peak),
which brought that ratio to ~1.
"""
from __future__ import annotations
import os, json
import numpy as np
from scipy import ndimage

BOX_R = 1.6          # half-size of the saved crop around the disc, in R
RL_ITERS = 10        # Richardson-Lucy iterations (Nico: "ten iterations look best")
LOCAL_SIG_PX = 20.0  # local-noise window of the prominence gate, px at REF_R
REF_R = 618.4        # the reference set's lunar radius; pixel constants scale from it
P_DEFAULT = (1.0, 0.059, 0.131)   # H-alpha colour (R=1) measured on the reference set
STACK_FILE = "prom_stack.npy"
STACK_META = "prom_stack.json"


def _ss(x):
    return x * x * (3.0 - 2.0 * x)


def disc_box(cy, cx, R, H, W, f=BOX_R):
    y0, x0 = max(int(cy - f * R), 0), max(int(cx - f * R), 0)
    y1, x1 = min(int(cy + f * R) + 1, H), min(int(cx + f * R) + 1, W)
    return y0, y1, x0, x1


def _polar_grid(shape, y0, x0, cy, cx, R):
    h, w = shape
    yy = np.arange(h, dtype=np.float32)[:, None] + (y0 - cy)
    xx = np.arange(w, dtype=np.float32)[None, :] + (x0 - cx)
    rr = np.hypot(yy, xx) / np.float32(R)
    ang = np.degrees(np.arctan2(yy, xx)) % 360.0
    return rr, ang


def _mad(v):
    m = np.median(v)
    return 1.4826 * float(np.median(np.abs(v - m))) + 1e-12, float(m)


def corona_ratio(rgb, rr, ang, nsec=180):
    """Corona R/G per 2 deg sector in the ring 1.03-1.25 R, smoothed over
    ~14 deg: the corona's own colour, so the excess is white-balance and
    sky-gradient independent."""
    ring = (rr > 1.03) & (rr < 1.25)
    r = rgb[..., 0] / np.maximum(rgb[..., 1], 1e-9)
    sec = (ang / (360.0 / nsec)).astype(int) % nsec
    k = np.full(nsec, np.nan)
    for s in range(nsec):
        v = r[ring & (sec == s)]
        if v.size > 50:
            k[s] = np.median(v)
    if np.all(np.isnan(k)):
        return np.full_like(r, np.nanmedian(r[ring]) if ring.any() else 1.0)
    k = np.where(np.isnan(k), np.nanmedian(k), k)
    k = ndimage.uniform_filter1d(k, 7, mode="wrap")
    return k[sec].astype(np.float32)


def prominence_colour(rgb, kR, rr):
    """(1, pG, pB) of the prominence light itself, from the purest cores in the
    merge with the corona under them removed. Falls back to the reference
    set's value when the set has no usable core."""
    ex = rgb[..., 0] - kR * rgb[..., 1]
    zone = (rr > 1.0) & (rr < 1.3)
    if zone.sum() < 1000:
        return np.array(P_DEFAULT)
    thr = np.percentile(ex[zone], 99.8)
    core = zone & (ex > thr)
    if core.sum() < 30:
        return np.array(P_DEFAULT)
    obs = rgb[core].astype(np.float64)
    k = np.median(kR[core])
    kB = np.median((rgb[..., 2] / np.maximum(rgb[..., 1], 1e-9))[(rr > 1.03) & (rr < 1.25)])
    c = np.array([k, 1.0, kB])
    p = np.array(P_DEFAULT)
    for _ in range(5):
        a = (obs[:, 0] - k * obs[:, 1]) / (1 - k * p[1])
        C = obs[:, 1] - a * p[1]
        pr = (obs - C[:, None] * c[None, :]) / np.maximum(a[:, None], 1e-9)
        p = np.median(pr, 0)
    p = p / max(p[0], 1e-9)
    if not (0.0 <= p[1] <= 0.3 and 0.0 <= p[2] <= 0.6):
        return np.array(P_DEFAULT)
    return p


# ----------------------------------------------------------------------------
# the short-exposure stack (built during the run, while the tiers are in memory)
# ----------------------------------------------------------------------------

def _psf_sigma(lum, rr, R):
    """Gaussian sigma of the optics+seeing blur, from the lunar limb: the edge
    of the median radial profile, width of its derivative peak."""
    band = (rr > 1 - 12 / R) & (rr < 1 + 12 / R)
    rp = (rr[band] - 1) * R
    v = lum[band]
    bins = np.arange(-12, 12.01, 0.25)
    idx = np.digitize(rp, bins)
    prof = np.array([np.median(v[idx == i]) if (idx == i).sum() > 20 else np.nan
                     for i in range(1, len(bins))])
    if np.isnan(prof).mean() > 0.3:
        return 2.0
    prof = np.where(np.isnan(prof), np.nanmedian(prof), prof)
    d = np.gradient(ndimage.gaussian_filter1d(prof, 1))
    xc = 0.5 * (bins[1:] + bins[:-1])
    i0 = int(np.argmax(d))
    sl = slice(max(i0 - 24, 0), i0 + 25)
    dd = np.clip(d[sl], 0, None)
    if dd.sum() <= 0:
        return 2.0
    mu = (dd * xc[sl]).sum() / dd.sum()
    sig = float(np.sqrt((dd * (xc[sl] - mu) ** 2).sum() / dd.sum()))
    return float(np.clip(sig, 0.7, 4.0))


def richardson_lucy(img, sigma, iters=RL_ITERS):
    """Plain RL with a Gaussian PSF (symmetric, so the adjoint is itself)."""
    img = np.maximum(img.astype(np.float32), 1e-6)
    est = img.copy()
    for _ in range(iters):
        conv = ndimage.gaussian_filter(est, sigma)
        est *= ndimage.gaussian_filter(img / np.maximum(conv, 1e-6), sigma)
    return est


def _sector_affine(E_ref, E_mov, rr, ang, R, nsec=24, patch=96):
    """Register the moving tier's red excess onto the merge's, sector by sector
    around the limb, and fit a similarity transform to the sector shifts.
    WHY NOT ONE SHIFT: on the reference set the prominence stack was off the
    merge by up to 5 px in different directions at different position angles
    (rotation/scale between the tiers and the merge), which a single shift
    cannot follow. Returns (matrix, offset) for ndimage.affine_transform, and a
    note."""
    from skimage.registration import phase_cross_correlation
    zone = (rr > 1.0) & (rr < 1.3)
    n, m = _mad(E_ref[zone])
    sE = ndimage.gaussian_filter(E_ref, 2)
    sec = (ang / (360.0 / nsec)).astype(int) % nsec
    h, w = E_ref.shape
    hp = patch // 2
    win = np.outer(np.hanning(patch), np.hanning(patch)).astype(np.float32)
    pts, sh = [], []
    for s in range(nsec):
        sel = zone & (sec == s)
        if not sel.any():
            continue
        v = np.where(sel, sE, -np.inf)
        i = int(np.argmax(v))
        py, px = divmod(i, w)
        if (v.flat[i] - m) / n < 20:
            continue
        if py < hp or px < hp or py + hp > h or px + hp > w:
            continue
        a = (E_ref[py - hp:py + hp, px - hp:px + hp].astype(np.float64) - m) * win
        b = (E_mov[py - hp:py + hp, px - hp:px + hp].astype(np.float64) - np.median(E_mov[zone])) * win
        a /= (a.std() + 1e-12); b /= (b.std() + 1e-12)
        try:
            d = phase_cross_correlation(a, b, upsample_factor=10)[0]
        except Exception:
            continue
        # (the returned error is not a usable quality figure across skimage
        # versions; the 20-sigma signal test above and the outlier pass below are)
        if np.hypot(*d) > 12:
            continue
        pts.append((py, px)); sh.append(d)
    c = np.array([h / 2.0, w / 2.0])
    if len(sh) == 0:
        return np.eye(2), np.zeros(2), "no sectors with prominence signal - not aligned"
    P = np.array(pts, float) - c; D = np.array(sh, float)
    if len(sh) < 3:
        t = np.median(D, 0)
        return np.eye(2), -t, f"{len(sh)} sector(s): shift {t[0]:+.2f},{t[1]:+.2f} px"
    # d = t + [[s,-a],[a,s]] p   (y,x ordering)
    A = np.zeros((2 * len(P), 4)); bvec = D.reshape(-1)
    A[0::2, 0] = 1; A[1::2, 1] = 1
    A[0::2, 2] = P[:, 0]; A[0::2, 3] = -P[:, 1]
    A[1::2, 2] = P[:, 1]; A[1::2, 3] = P[:, 0]
    sol, *_ = np.linalg.lstsq(A, bvec, rcond=None)
    ty, tx, sc, ro = sol
    res = bvec - A @ sol
    # drop outliers once
    r2 = np.hypot(res[0::2], res[1::2])
    good = r2 < max(1.0, 3 * np.median(r2))
    if good.sum() >= 3 and not good.all():
        gi = np.repeat(good, 2)
        sol, *_ = np.linalg.lstsq(A[gi], bvec[gi], rcond=None)
        ty, tx, sc, ro = sol
    M = np.array([[sc, -ro], [ro, sc]])
    # out(q) = in(q - d(q)),  d(q) = t + M (q - c)
    matrix = np.eye(2) - M
    offset = M @ c - np.array([ty, tx])
    return matrix, offset, (f"{int(good.sum())}/{len(sh)} sectors: shift {ty:+.2f},{tx:+.2f} px, "
                            f"scale {sc * 1e3:+.2f}e-3, rotation {np.degrees(ro) * 60:+.1f} arcmin")


def build_stack(wd, tiers, tier_box, cy, cx, R, H, W, progress=None):
    """tiers: exposure times to use (fast ones). tier_box(s, box) -> (rgb, valid)
    for the disc box in the cropped frame, in merge units. Writes the
    deconvolved stack and its metadata into the work dir."""
    log = (lambda m: progress.log(m, None)) if progress is not None else print
    y0, y1, x0, x1 = box = disc_box(cy, cx, R, H, W)
    hdr = np.load(os.path.join(wd, "hdr_rgb.npy"), mmap_mode="r")
    ref = np.asarray(hdr[y0:y1, x0:x1], np.float32)
    rr, ang = _polar_grid(ref.shape[:2], y0, x0, cy, cx, R)
    kR = corona_ratio(ref, rr, ang)
    E_ref = ref[..., 0] - kR * ref[..., 1]
    acc = np.zeros_like(ref); accw = np.zeros(ref.shape[:2], np.float32)
    notes = {}
    for s in tiers:
        rgb, valid = tier_box(s, box)
        E = rgb[..., 0] - corona_ratio(rgb, rr, ang) * rgb[..., 1]
        mat, off, note = _sector_affine(E_ref, E, rr, ang, R)
        if not np.allclose(mat, np.eye(2)) or np.any(off != 0):
            rgb = np.stack([ndimage.affine_transform(rgb[..., c], mat, off, order=3,
                                                     mode="nearest") for c in range(3)], -1)
            valid = ndimage.affine_transform(valid.astype(np.float32), mat, off, order=1,
                                             mode="constant", cval=0)
        wgt = np.float32(s) * np.clip(valid, 0, 1)
        acc += rgb * wgt[..., None]; accw += wgt
        notes[f"{s:g}"] = note
        log(f"prominence stack: tier {s:g}s: {note}")
        del rgb, valid, E
    st = acc / np.maximum(accw, 1e-9)[..., None]
    st[accw <= 0] = ref[accw <= 0]
    lum = 0.2126 * st[..., 0] + 0.7152 * st[..., 1] + 0.0722 * st[..., 2]
    sig = _psf_sigma(lum, rr, R)
    off = np.float32(max(0.0, -float(st.min())) + 1e-3)
    st = np.stack([richardson_lucy(st[..., c] + off, sig) - off for c in range(3)], -1)
    np.save(os.path.join(wd, STACK_FILE), st.astype(np.float32))
    json.dump({"box": [y0, y1, x0, x1], "psf_sigma": sig, "rl_iters": RL_ITERS,
               "tiers": [float(t) for t in tiers], "align": notes},
              open(os.path.join(wd, STACK_META), "w"), indent=1)
    log(f"[ok] prominence stack: {len(tiers)} tier(s), PSF sigma {sig:.2f} px, "
        f"{RL_ITERS} RL iterations")


# ----------------------------------------------------------------------------
# the layer
# ----------------------------------------------------------------------------

def _detect(M, kR, rr, ang, R, cy, cx, y0, x0, pG):
    """Where the prominences are, from the merge's red excess, with the
    chromosphere band removed. Returns (sig, resxy, keep, ramp, alpha_soft,
    nprom, a_m) on the box grid."""
    q = R / REF_R
    a_m = (M[..., 0] - kR * M[..., 1]) / (1 - kR * pG)
    zz = (rr > 1.03) & (rr < 1.3)
    moon = np.clip((rr * R - (R - 1)) / 2, 0, 1)
    sm = ndimage.gaussian_filter(a_m, 1.5)
    n_sm, b_sm = _mad(sm[zz]); sig = (sm - b_sm) / n_sm
    # chromosphere band out: polar, background = 30th pct along +-2.5 deg
    NA, NR = 3600, int(300 * q) + 20
    r0 = R - 10 * q
    th = np.radians(np.arange(NA) * 0.1); rad = r0 + np.arange(NR)
    py = (cy - y0) + rad[:, None] * np.sin(th)[None, :]
    px = (cx - x0) + rad[:, None] * np.cos(th)[None, :]
    P = ndimage.map_coordinates(sig, [py, px], order=1, mode="nearest")
    bg = ndimage.percentile_filter(P, 30, size=(1, 51), mode="wrap")
    res = np.maximum(P - bg, 0)
    hh = rad - R
    tall = ((P > 12) & (hh[:, None] > 45 * q)).any(0)
    tall = ndimage.gaussian_filter1d(ndimage.binary_dilation(tall, iterations=3).astype(float),
                                     2, mode="wrap")
    res = np.maximum(res, P * tall[None, :])
    ri = rr * R - r0; ti = (ang * 10) % NA
    resxy = ndimage.map_coordinates(np.concatenate([res, res[:, :1]], 1),
                                    [np.clip(ri, 0, NR - 1), ti], order=1, mode="nearest")
    resxy[ri > NR - 1] = sig[ri > NR - 1]
    core = (resxy > 12) & (rr > 1.0)
    lab, k = ndimage.label(core)
    minpx = max(10, int(40 * q * q))
    sz = ndimage.sum(core, lab, range(1, k + 1)) if k else np.zeros(0)
    keep = np.isin(lab, 1 + np.flatnonzero(sz >= minpx))
    nprom = int((sz >= minpx).sum())
    ramp = _ss(np.clip((resxy - 6) / 14, 0, 1))
    alpha_soft = ndimage.gaussian_filter(ramp * ndimage.binary_dilation(keep, iterations=3), 1.0) * moon
    return sig, resxy, keep, ramp, alpha_soft, nprom, a_m


def build_layer(wd, geo, log=print):
    """Returns (box, rgb, alpha): rgb in display encoding 0..1 and alpha 0..1
    for the disc box (y0,y1,x0,x1) of the cropped frame."""
    cy, cx, R = float(geo["cy"]), float(geo["cx"]), float(geo["R"])
    hdr = np.load(os.path.join(wd, "hdr_rgb.npy"), mmap_mode="r")
    H, W = hdr.shape[:2]
    y0, y1, x0, x1 = disc_box(cy, cx, R, H, W)
    M = np.asarray(hdr[y0:y1, x0:x1], np.float64)
    rr, ang = _polar_grid(M.shape[:2], y0, x0, cy, cx, R)
    q = R / REF_R                                   # pixel constants scale with the disc
    hgt = (rr - 1.0) * R
    kR = corona_ratio(M, rr, ang)
    p = prominence_colour(M, kR, rr)
    pG = p[1]
    zz = (rr > 1.03) & (rr < 1.3)
    moon = np.clip((rr * R - (R - 1)) / 2, 0, 1)
    sig, resxy, keep, ramp, alpha_soft, nprom, a_m = _detect(M, kR, rr, ang, R, cy, cx, y0, x0, pG)
    # --- the sharp stack, when the run made one ---
    n_m, b_m = _mad(a_m[zz]); am = np.clip(a_m - b_m, 0, None)
    sp = os.path.join(wd, STACK_FILE)
    have_stack = False
    if os.path.exists(sp):
        try:
            meta = json.load(open(os.path.join(wd, STACK_META)))
            if list(meta["box"]) == [y0, y1, x0, x1]:
                S = np.load(sp).astype(np.float64)
                have_stack = True
        except Exception:
            have_stack = False
    if have_stack:
        kS = corona_ratio(S, rr, ang)
        a_s = ndimage.gaussian_filter((S[..., 0] - kS * S[..., 1]) / (1 - kS * pG), 0.8)
        n_s, b_s = _mad(a_s[zz]); a_s = np.clip(a_s - b_s, 0, None)
        # same units as the merge: compare both smoothed to the merge's blur
        pr = keep & (rr > 1.0)
        if pr.sum() > 50:
            g1 = ndimage.gaussian_filter(a_s, 4)[pr]; g2 = ndimage.gaussian_filter(am, 4)[pr]
            scale = float(np.median(g2 / np.maximum(g1, 1e-9)))
            if np.isfinite(scale) and 0.2 < scale < 5:
                a_s *= scale; n_s *= scale
        wt = np.clip((ndimage.gaussian_filter(a_s, 2) / n_s - 20) / 40, 0, 1)
        amt = wt * a_s + (1 - wt) * am
        shape = ndimage.gaussian_filter(a_s, 1.2)
    else:
        amt = am
        shape = ndimage.gaussian_filter(am, 1.2)
    # --- tight edge, relative to each prominence's own peak ---
    gate = ndimage.binary_dilation(keep, iterations=4)
    lpk = ndimage.maximum_filter(np.where(gate, shape, 0), size=int(41 * q) | 1)
    rel = np.where(lpk > 0, shape / np.maximum(lpk, 1e-12), 0)
    tr = _ss(np.clip((rel - 0.12) / (0.35 - 0.12), 0, 1))
    alpha_t = ndimage.gaussian_filter(tr * gate * np.clip(ramp * 1.5, 0, 1), 0.7) * moon
    # prominences that rise well above the band keep the soft mask whole: their
    # faint inner parts sit far below their bright core
    lab2, k2 = ndimage.label(ndimage.binary_dilation(keep, iterations=4))
    if k2:
        mxh = ndimage.maximum(np.where(rr > 1, hgt, 0), lab2, range(1, k2 + 1))
        tallxy = ndimage.gaussian_filter(
            np.isin(lab2, 1 + np.flatnonzero(np.asarray(mxh) > 45 * q)).astype(float), 2)
    else:
        tallxy = np.zeros_like(alpha_t)
    alpha = np.clip(np.maximum(alpha_t, alpha_soft * tallxy), 0, 1)
    # --- brightness and colour ---
    inside = alpha > 0.5
    peak = np.percentile(amt[inside], 99.9) if inside.sum() > 20 else max(amt.max(), n_m)
    t = np.clip(np.log1p(amt / n_m) / np.log1p(max(peak, n_m) / n_m), 0, 1)
    cref = np.array([np.median(kR[zz]), 1.0,
                     np.median((M[..., 2] / np.maximum(M[..., 1], 1e-9))[zz])])
    hue = p / cref; hue = np.clip(hue / hue.max(), 0, 1)
    hue_d = hue ** (1 / 2.2)
    rgb = (t[..., None] * hue_d[None, None, :]).astype(np.float32)
    log(f"[ok] prominence layer: {nprom} prominence(s), "
        f"{'sharp stack' if have_stack else 'merge only (no cached stack)'}, "
        f"G/R {p[1]:.3f} B/R {p[2]:.3f}")
    return (y0, y1, x0, x1), rgb, alpha.astype(np.float32), have_stack


def export_layer(wd, geo, out_path, orient="", size="full", log=print):
    """Write the layer as a 16-bit RGBA TIFF on exactly the export's grid."""
    import tifffile
    from . import icc
    from .render import apply_orient
    hdr = np.load(os.path.join(wd, "hdr_rgb.npy"), mmap_mode="r")
    H, W = hdr.shape[:2]
    (y0, y1, x0, x1), rgb, alpha, have_stack = build_layer(wd, geo, log=log)
    full = np.zeros((H, W, 4), np.float32)
    full[y0:y1, x0:x1, :3] = rgb
    full[y0:y1, x0:x1, 3] = alpha
    full = apply_orient(full, orient)
    if size == "half":
        h2, w2 = full.shape[0] // 2 * 2, full.shape[1] // 2 * 2
        f = full[:h2, :w2].reshape(h2 // 2, 2, w2 // 2, 2, 4)
        a = f[..., 3].mean(axis=(1, 3))
        c = (f[..., :3] * f[..., 3:4]).sum(axis=(1, 3)) / np.maximum(
            f[..., 3].sum(axis=(1, 3)), 1e-9)[..., None]
        full = np.concatenate([c, a[..., None]], -1)
    prof = icc.srgb_profile()
    tifffile.imwrite(out_path, (np.clip(full, 0, 1) * 65535 + 0.5).astype(np.uint16),
                     photometric="rgb", extrasamples=[2], compression="zlib",
                     extratags=[(34675, 1, len(prof), prof, False)],
                     description="eclipseforgehdr prominence layer (RGBA, unassociated alpha); "
                                 "same grid as the composite export")
    return have_stack


_DET_CACHE = {}


def detect(wd, geo, shape, log=None):
    """The prominence detection on the full layer grid: a soft 0..1 map of
    where prominences are (band removed, dilated 3 px, feathered), for the
    app's prominence gate and the convolution mask. None when there is no
    merge to read.

    WHY (0.23.10 lab): both older detectors threshold the redness R/((G+B)/2)
    of ONE fast tier. Measured on the 600 mm set after a full run: of the 15
    prominences this detection finds, the older gate covered 4 fully, 2 in
    part (the 225 deg one 76 %) and 9 not at all -- among them the 52 px
    prominence at 68 deg. What the gate had and this lacks: ~100 px of limb
    rim. Callers take the UNION, so nothing either finds is lost."""
    p = os.path.join(wd, "hdr_rgb.npy")
    if not os.path.exists(p):
        return None
    key = (wd, os.path.getmtime(p), tuple(shape))
    if key in _DET_CACHE:
        return _DET_CACHE[key]
    try:
        cy, cx, R = float(geo["cy"]), float(geo["cx"]), float(geo["R"])
        hdr = np.load(p, mmap_mode="r")
        H, W = hdr.shape[:2]
        if (H, W) != tuple(shape):
            return None
        y0, y1, x0, x1 = disc_box(cy, cx, R, H, W)
        M = np.asarray(hdr[y0:y1, x0:x1], np.float64)
        rr, ang = _polar_grid(M.shape[:2], y0, x0, cy, cx, R)
        kR = corona_ratio(M, rr, ang)
        pG = prominence_colour(M, kR, rr)[1]
        _, _, _, _, soft, nprom, _ = _detect(M, kR, rr, ang, R, cy, cx, y0, x0, pG)
        out = np.zeros((H, W), np.float32)
        out[y0:y1, x0:x1] = soft
        if log is not None:
            log(f"  prominence detection (red excess): "
                f"{nprom} prominence(s), {int((soft > 0.3).sum()) / 1e3:.1f}k px")
        _DET_CACHE.clear()
        _DET_CACHE[key] = out
        return out
    except Exception as e:
        if log is not None:
            log(f"[warn] prominence detection failed ({e}); older detectors only")
        return None


# ============================================================================
# 0.24 LAB: THE PROMINENCE LAYER FROM THE SHORT FRAMES, ON A CORONA WITHOUT THEM
# ============================================================================
#
# WHY A NEW DESIGN (Oct 2026, Nico's 600 mm set, prototype in status
# 2026-10-06-prominence-layer-from-short-frames and -open-points).
#
# Everything before this tried to show the prominences THROUGH the corona
# stack: brighten them inside a gate, give them their own detail, exclude them
# from the convolutions. Every one of those draws an edge somewhere -- the gate,
# the exclusion mask, the hole in MGN -- and the partial-convolution masks
# amplify any edge into a line. Nico's list after the Oct 5 run: dark rims,
# smearing, false pink patches, hard edges, no blend, a poor fit with the
# corona. Two causes, measured:
#
#   1. THE STACK SMEARS THEM. The long tiers carry most of the weight and are
#      softened by seeing; the Moon covers one side and uncovers the other over
#      the bracket; and small clipped islands get a weight STEP in the Exact-edge
#      feather (pipeline._feather_weight, fixed there). A single 1/500 s frame
#      shows prominences the stack loses (the 2 o'clock spike).
#   2. MASKS LEAVE SEAMS. Any spatial mask -- even a soft one -- prints in Hill.
#
# So the corona is made WITHOUT prominences by COLOUR ALONE (clean_corona: the
# merge minus its own H-alpha light, a * p, with no mask anywhere), every layer
# is built from that, and the prominences come back as their own layer
# (build_light), stacked from the SHORT frames only: sharp, unclipped and
# colour-true. The render lays it over the finished picture with its own
# transparency (render.py, promLayer). The same layer is what the RGBA export
# writes for Photoshop.

LIGHT_FILE = "prom_light.npy"        # the prominence light, linear RGB, merge units, disc box
LIGHT_META = "prom_light.json"
SHORT_FILE = "prom_short.npy"        # the short-frame stack (before the light is cut out of it)
SHORT_W = "prom_short_w.npy"         # its coverage (sum of weights), 0 = no frame saw the pixel
SHORT_META = "prom_short.json"
CORONA_RGB = "corona_rgb.npy"        # the merge without prominence light, full frame
CORONA_LUM = "corona_lum.npy"
REMOVED_FILE = "prom_removed.npy"    # luminance taken out of the merge, disc box


def layer_mode():
    """The 0.24 prominence layer. ECLIPSEFORGE_PROM_CLASSIC=1 restores the gate."""
    return os.environ.get("ECLIPSEFORGE_PROM_CLASSIC", "") in ("", "0")


def _lum(x):
    return 0.2126 * x[..., 0] + 0.7152 * x[..., 1] + 0.0722 * x[..., 2]


def _height(geo, shape, y0, x0):
    """Height above the MEASURED limb (the per-azimuth mask radius) on a box grid."""
    from .detail import limb_radius_map
    cy, cx, R = float(geo["cy"]), float(geo["cx"]), float(geo["R"])
    h, w = shape
    yy = np.arange(h, dtype=np.float32)[:, None] + np.float32(y0 - cy)
    xx = np.arange(w, dtype=np.float32)[None, :] + np.float32(x0 - cx)
    r = np.sqrt(yy * yy + xx * xx)
    margin = float(geo.get("limb_margin", float(geo.get("Rmask", R + 4.0)) - R))
    prof = geo.get("limb_prof")
    if prof:
        Rm = limb_radius_map(prof, (h, w), cy - y0, cx - x0, margin)
    else:
        Rm = np.float32(R + margin)
    return (r - Rm).astype(np.float32)


def _fit_moon(G, c0, R0):
    """Lunar limb of one frame on a box grid: steepest rise of log G along 360
    rays, circle fit with outlier rejection. Returns (cy, cx, R)."""
    cy, cx = float(c0[0]), float(c0[1])
    Gs = ndimage.gaussian_filter(np.asarray(G, np.float32), 1.5)
    # floored well above the noise inside the disc: a log of near-zero noise
    # has the steepest "rise" in the frame (it found R = 585 for 617 without this): a tenth of the inner corona
    _yy = np.arange(Gs.shape[0], dtype=np.float32)[:, None] - np.float32(cy)
    _xx = np.arange(Gs.shape[1], dtype=np.float32)[None, :] - np.float32(cx)
    _rr = np.sqrt(_yy * _yy + _xx * _xx)
    _an = Gs[(_rr > 1.12 * R0) & (_rr < 1.3 * R0)]
    floor = 0.1 * float(np.median(_an)) if _an.size else 1e-6
    L = np.log10(np.clip(Gs, max(floor, 1e-12), None))
    rs = np.arange(0.9 * R0, 1.1 * R0, 0.5)
    R = R0
    for _ in range(2):
        pts = []
        for a in np.radians(np.arange(0, 360, 1.0)):
            v = ndimage.map_coordinates(L, [cy + rs * np.sin(a), cx + rs * np.cos(a)], order=1)
            i = int(np.argmax(np.gradient(v)))
            pts.append((cy + rs[i] * np.sin(a), cx + rs[i] * np.cos(a)))
        pts = np.array(pts)
        for _ in range(4):
            A = np.c_[2 * pts[:, 1], 2 * pts[:, 0], np.ones(len(pts))]
            b = pts[:, 1] ** 2 + pts[:, 0] ** 2
            s = np.linalg.lstsq(A, b, rcond=None)[0]
            xc, yc = s[0], s[1]
            R = float(np.sqrt(max(s[2] + xc * xc + yc * yc, 1.0)))
            res = np.hypot(pts[:, 1] - xc, pts[:, 0] - yc) - R
            k = np.abs(res) < max(3 * 1.4826 * np.median(np.abs(res)), 0.7)
            if k.sum() < 30:
                break
            pts = pts[k]
        cy, cx = float(yc), float(xc)
    return cy, cx, R


def _red_hp(rgb, kR, moon, R):
    """The positive red excess, band-passed (1.5 - 12 px) and windowed from 6 px
    outside THIS FRAME'S OWN MOON out to 1.4 R: what a frame is registered on.
    The window has to follow the frame's own limb: the chromosphere rim moves
    with the Moon, not the Sun, and inside the merge's window it pulled the
    match toward the Moon (correlation 0.2-0.3 against 0.7-0.97 this way, on
    Nico's short frames)."""
    x = np.clip(rgb[..., 0] - kR * rgb[..., 1], 0, None)
    x = ndimage.gaussian_filter(x, 1.5) - ndimage.gaussian_filter(x, 12.0)
    my, mx, mR = moon
    yy = np.arange(rgb.shape[0], dtype=np.float32)[:, None] - np.float32(my)
    xx = np.arange(rgb.shape[1], dtype=np.float32)[None, :] - np.float32(mx)
    rf = np.sqrt(yy * yy + xx * xx)
    win = np.clip((rf - mR - 6.0) / 6.0, 0, 1) * (rf < 1.4 * R)
    return (x * win).astype(np.float64)


def _xcorr_shift(a, b, maxd):
    """Shift (dy, dx) that moves b onto a, by FFT cross-correlation, searched
    within +-maxd px, sub-pixel by a parabola."""
    H, W = a.shape
    F = np.fft.fftshift(np.fft.irfft2(np.fft.rfft2(a) * np.conj(np.fft.rfft2(b)), s=a.shape))
    cy, cx = H // 2, W // 2
    m = int(maxd)
    sub = F[cy - m:cy + m + 1, cx - m:cx + m + 1]
    i, j = np.unravel_index(int(np.argmax(sub)), sub.shape)
    i += cy - m; j += cx - m

    def par(v0, v1, v2):
        d = v0 - 2 * v1 + v2
        return 0.0 if d == 0 else 0.5 * (v0 - v2) / d
    dy = i - cy + (par(F[i - 1, j], F[i, j], F[i + 1, j]) if 0 < i < H - 1 else 0.0)
    dx = j - cx + (par(F[i, j - 1], F[i, j], F[i, j + 1]) if 0 < j < W - 1 else 0.0)
    nrm = np.sqrt((a * a).sum() * (b * b).sum()) + 1e-30
    return float(dy), float(dx), float(F[i, j] / nrm)


def build_short_stack(wd, items, cy, cx, R, H, W, progress=None):
    """The prominences from the SHORT frames, registered ON THE SUN.

    items: list of dicts {name, sec, get} where get(box) -> (rgb, clipped):
    the frame's pixels on the merge's disc box (linear RGB, per second; tiers
    already in merge units, extra frames not) and a bool map of clipped or
    missing pixels. Roughly placed is enough -- each item is registered here.

    Per item, measured not assumed (all of it was wrong somewhere on Nico's
    set when assumed):
      * REGISTRATION on the red excess (the prominences), not the Moon: the
        Moon moves against the Sun by ~8 px over 22 s of short frames. FFT
        cross-correlation for the coarse shift (frames of the discarded short
        series sat up to 35 px from the merge), three sub-pixel refinements,
        translation only, with each frame's own Moon windowed out (see _red_hp).
      * THE FRAME'S OWN MOON, fitted on the registered frame: a prominence at
        the limb is covered in some frames and not in others, so each frame
        contributes only outside its own disc (ramp 12 px).
      * CALIBRATION ON THE CORONA (1.1-1.4 R), not on the prominences: those
        change through the series by occultation (0.1x to 2x between groups).
      * CLIPPED PIXELS by distance (ramp 8 px), so a small clipped island has
        no weight step around it -- the same fix as the merge's.
    Weight: exposure time x those ramps. Writes SHORT_FILE, SHORT_W, SHORT_META."""
    log = (lambda m: progress.log(m, None)) if progress is not None else print
    y0, y1, x0, x1 = box = disc_box(cy, cx, R, H, W)
    hdr = np.load(os.path.join(wd, "hdr_rgb.npy"), mmap_mode="r")
    ref = np.asarray(hdr[y0:y1, x0:x1], np.float32)
    rr, ang = _polar_grid(ref.shape[:2], y0, x0, cy, cx, R)
    q = R / REF_R
    kref = corona_ratio(ref, rr, ang)
    c0 = (cy - y0, cx - x0)
    E_ref = _red_hp(ref, kref, _fit_moon(ref[..., 1], c0, R), R)
    ring = (rr > 1.1) & (rr < 1.4)
    Gref = ndimage.gaussian_filter(ref[..., 1], 3)
    acc = np.zeros(ref.shape, np.float64); accw = np.zeros(ref.shape[:2], np.float64)
    notes = []

    def ss(t):
        t = np.clip(t, 0, 1)
        return t * t * (3 - 2 * t)
    for it in items:
        try:
            rgb, clipped = it["get"](box)
            rgb = np.asarray(rgb, np.float32); clipped = np.asarray(clipped, bool)
            kR = corona_ratio(rgb, rr, ang)
            mcy, mcx, mR = _fit_moon(rgb[..., 1], c0, R)
            if np.hypot(mcy - c0[0], mcx - c0[1]) > 80 * q or abs(mR - R) > 0.05 * R:
                log(f"[warn] prominence layer: {it['name']} not used (Moon not found)")
                notes.append({"name": it["name"], "sec": it["sec"], "used": False})
                continue
            E = _red_hp(rgb, kR, (mcy, mcx, mR), R)
            if clipped.any():                      # a clipped prominence is not its shape
                E *= ~ndimage.binary_dilation(clipped, iterations=int(round(4 * q)) + 1)
            dy = dx = 0.0; cc = 0.0
            for k_ in range(4):                    # coarse, then three sub-pixel refinements
                Es = E if k_ == 0 else ndimage.shift(E, (dy, dx), order=1, mode="constant")
                a_, b_, cc = _xcorr_shift(E_ref, Es, 80 * q if k_ == 0 else 6)
                dy += a_; dx += b_
            if cc < 0.3:
                log(f"[warn] prominence layer: {it['name']} not used (no match, corr {cc:.2f})")
                notes.append({"name": it["name"], "sec": it["sec"], "used": False, "corr": cc})
                continue
            rgb = np.stack([ndimage.shift(rgb[..., c], (dy, dx), order=3, mode="constant")
                            for c in range(3)], -1)
            inside = ndimage.shift(np.ones(rgb.shape[:2], np.float32), (dy, dx), order=1,
                                   mode="constant", cval=0) > 0.999
            clipped = (ndimage.shift(clipped.astype(np.float32), (dy, dx), order=1,
                                     mode="constant", cval=1) > 0.01) | ~inside
            mcy += dy; mcx += dx
            anote = "translation only"
            yy = np.arange(rgb.shape[0], dtype=np.float32)[:, None] - np.float32(mcy)
            xx = np.arange(rgb.shape[1], dtype=np.float32)[None, :] - np.float32(mcx)
            wmoon = ss((np.sqrt(yy * yy + xx * xx) - mR - 1.0) / (12.0 * q))
            d = ndimage.distance_transform_edt(~clipped) if clipped.any() else np.full(clipped.shape, 1e9)
            wclip = ss(d / (8.0 * q))
            w = wmoon * wclip
            m = ring & (w > 0.99)
            k = 1.0
            if m.sum() > 2000:
                Gs = ndimage.gaussian_filter(rgb[..., 1], 3)
                k = float(np.median(Gref[m] / np.maximum(Gs[m], 1e-12)))
                if not np.isfinite(k) or not (0.3 < k < 3.0):
                    k = 1.0
            wt = float(it["sec"]) * w
            acc += rgb.astype(np.float64) * k * wt[..., None]
            accw += wt
            note = {"name": it["name"], "sec": float(it["sec"]), "used": True,
                    "shift": [round(dy, 2), round(dx, 2)], "corr": round(cc, 3), "fine": anote,
                    "moon": [round(mcy + y0, 2), round(mcx + x0, 2), round(mR, 2)],
                    "corona_cal": round(k, 4)}
            notes.append(note)
            log(f"prominence layer: {it['name']} ({it['sec']:g}s) shift {dy:+.1f},{dx:+.1f} px "
                f"(corr {cc:.2f}), {anote}; corona cal {k:.3f}")
            del rgb, clipped, w, wt, d
        except Exception as e:
            log(f"[warn] prominence layer: {it.get('name')} skipped ({e})")
            notes.append({"name": it.get("name"), "used": False, "error": str(e)})
    if accw.max() <= 0:
        log("[warn] prominence layer: no short frame usable; using the merge")
        return False
    S = acc / np.maximum(accw, 1e-30)[..., None]
    S[accw <= 0] = ref[accw <= 0]
    np.save(os.path.join(wd, SHORT_FILE), S.astype(np.float32))
    np.save(os.path.join(wd, SHORT_W), accw.astype(np.float32))
    json.dump({"box": [y0, y1, x0, x1], "items": notes}, open(os.path.join(wd, SHORT_META), "w"),
              indent=1)
    n = sum(1 for x in notes if x.get("used"))
    log(f"[ok] prominence layer: {n} short frame(s)/tier(s) stacked on the Sun")
    return True


def _fill_core_polar(cor, core, cy, cx, hmap, q=1.0, feather=4.0):
    """Fill `core` (bool) in `cor` (H,W,C linear) from the corona around it,
    on a grid of HEIGHT above the measured limb x angle about the centre
    (cy, cx), in log (0.24, Oct 6 -- replaces the isotropic blur fill, whose
    flat patch and 2 px boundary the detail filters drew as an arc around the
    large prominence):
      A  each row of constant height interpolated across the hole's angular
         span from an 8q px window beside it -- the limb's steep profile and
         the fall-off run straight through;
      T  the corona's fine texture (log, above ~6q px) from beside the hole,
         turned about the centre into it -- from the left near the left edge,
         from the right near the right edge -- so no detail filter sees a
         smooth patch.
    Measured on the 600 mm set's 10 o'clock prominence: detail (|log high-pass|)
    on the fill's edge band 1.12x the surroundings with the old fill, 0.78x
    with this one. Returns the filled image and the blend weight (1 on the
    core, a `feather`q px soft edge outside it)."""
    H, W = core.shape
    ys, xs = np.nonzero(core)
    if ys.size == 0:
        return cor, np.zeros((H, W))
    yy, xx = np.mgrid[:H, :W].astype(np.float64)
    th = np.arctan2(yy - cy, xx - cx)
    t0 = float(np.angle(np.mean(np.exp(1j * th[ys, xs]))))
    rc = hmap[ys, xs].astype(np.float64)          # rows are HEIGHT above the measured limb
    rlimb_map = np.hypot(yy - cy, xx - cx) - hmap
    rmean = float(np.mean(np.hypot(ys - cy, xs - cx)))
    pad = 6.0 * q / rmean; win = 8.0 * q / rmean
    dth_c = np.angle(np.exp(1j * (th[ys, xs] - t0)))
    span = dth_c.max() - dth_c.min() + 2 * pad
    ta, tb = dth_c.min() - pad - win - span, dth_c.max() + pad + win + span
    marg = 12.0 * q
    ra, rb = rc.min() - marg, rc.max() + marg
    dr = 0.5; dt = 0.5 / rmean
    ri = np.arange(ra, rb + dr, dr); tj = np.arange(ta, tb + dt, dt)
    # the limb radius along each column (a function of angle only)
    rl = ndimage.map_coordinates(rlimb_map, [cy + rmean * np.sin(tj + t0), cx + rmean * np.cos(tj + t0)], order=1, mode="nearest")
    HH, TT = np.meshgrid(ri, tj, indexing="ij")
    RR = HH + rl[None, :]
    PY = cy + RR * np.sin(TT + t0); PX = cx + RR * np.cos(TT + t0)
    C = cor.shape[2]
    L = np.log(np.maximum(cor, 1e-6))
    P = np.stack([ndimage.map_coordinates(L[..., k], [PY, PX], order=1, mode="nearest") for k in range(C)], -1)
    grow = max(2, int(round(1.5 * feather * q)))
    hole = ndimage.binary_dilation(core, iterations=grow)
    Pm = ndimage.map_coordinates(hole.astype(np.float64), [PY, PX], order=1, mode="constant") > 0.01
    cols = np.flatnonzero(Pm.any(0)); ja, jb = int(cols.min()), int(cols.max())
    nwin = max(3, int(round(win / dt)))
    sm = 6.0 * q / 0.5                       # 6q px in grid samples
    Ps = ndimage.gaussian_filter(P, (sm, sm, 0))
    # A: per circle, interpolation across [ja, jb] from the windows beside it
    A = Ps.copy()
    # from the SMOOTHED corona: the texture T below carries the rest, so the
    # limb's curvature is not counted twice (it was: +30 % on the limb)
    lw = Ps[:, max(0, ja - nwin):ja]; rw = Ps[:, jb + 1:jb + 1 + nwin]
    vl = np.median(lw, 1); vr = np.median(rw, 1)
    tl = ja - 0.5 * lw.shape[1]; tr = jb + 0.5 * rw.shape[1]
    w = np.clip((np.arange(ja, jb + 1) - tl) / max(tr - tl, 1.0), 0, 1)[None, :, None]
    A[:, ja:jb + 1] = vl[:, None, :] * (1 - w) + vr[:, None, :] * w
    # T: the corona's fine texture turned in from both sides
    T = P - Ps
    nsh = int(jb - ja + 1 + 2 * max(3, int(round(pad / dt))))
    Tl = np.roll(T, nsh, axis=1); Tr = np.roll(T, -nsh, axis=1)
    wt = np.clip((np.arange(len(tj)) - ja) / max(jb - ja, 1), 0, 1)[None, :, None]
    TF = (Tl * (1 - wt) + Tr * wt) / np.sqrt((1 - wt) ** 2 + wt ** 2)
    F = np.where(Pm[..., None], A + TF, P)
    # back to cartesian over the grown hole
    yh, xh = np.nonzero(hole)
    rh = hmap[yh, xh].astype(np.float64); thh = np.angle(np.exp(1j * (np.arctan2(yh - cy, xh - cx) - t0)))
    fi = (rh - ra) / dr; fj = (thh - ta) / dt
    vals = np.stack([np.exp(ndimage.map_coordinates(F[..., k], [fi, fj], order=1, mode="nearest")) for k in range(C)], -1)
    filled = np.zeros_like(cor); filled[yh, xh] = vals
    wc = np.clip(ndimage.gaussian_filter(core.astype(np.float64), feather * q) * 2.0, 0, 1)
    wc = np.maximum(wc, core.astype(float)); wc[~hole] = 0.0
    wc = wc * wc * (3 - 2 * wc)
    return cor * (1 - wc[..., None]) + filled * wc[..., None], wc


def clean_corona(wd, geo, progress=None):
    """The merge without its prominence light, BY COLOUR ALONE, written as
    CORONA_RGB / CORONA_LUM beside the untouched merge.

    corona = merge - a * p, with a = (R - kR G) / (1 - kR pG): the red excess
    over the corona's own colour kR, in units of the prominence colour p. It is
    subtracted everywhere inside 1.45 R (fading to nothing by 1.57 R), signed
    and ungated: no mask, so nothing can draw an edge, and on the corona the
    noise of both signs goes with it (its R becomes its own colour times G).
    Where the prominence is most of the light AND the region is large (the big
    10 o'clock core, >= 2000 px at the reference scale), the colour model is not
    exact (He D3, H-beta stay behind) and that core is filled from around it
    along rows of constant height, with the corona's own texture (see
    _fill_core_polar) -- its edge must not show where the layer thins out."""
    log = (lambda m: progress.log(m, None)) if progress is not None else print
    cy, cx, R = float(geo["cy"]), float(geo["cx"]), float(geo["R"])
    hdr = np.load(os.path.join(wd, "hdr_rgb.npy"), mmap_mode="r")
    H, W = hdr.shape[:2]
    y0, y1, x0, x1 = disc_box(cy, cx, R, H, W)
    M = np.asarray(hdr[y0:y1, x0:x1], np.float64)
    rr, ang = _polar_grid(M.shape[:2], y0, x0, cy, cx, R)
    h = _height(geo, M.shape[:2], y0, x0)
    q = R / REF_R
    kR = corona_ratio(M, rr, ang)
    p = np.asarray(prominence_colour(M, kR, rr), np.float64)
    a = (M[..., 0] - kR * M[..., 1]) / (1.0 - kR * p[1])
    a *= 1.0 - np.clip((rr - 1.45) / 0.12, 0, 1)
    # PIXEL RESOLUTION ONLY WHERE THERE IS A PROMINENCE (2026-10-08). `a` is a
    # per-pixel colour difference, so on the corona -- where it is zero on
    # average -- it is the noise of R and G, and subtracting a * p from every
    # channel put that noise into the luminance every layer is built on:
    # measured at 1.3-1.8 R on the 600 mm set, the pixel noise of corona_lum
    # was 1.66x the merge's and 8-sigma outliers 7x as many. Now the smooth
    # part of `a` (4 px) is used everywhere and the pixel-level part only
    # inside the prominence patches (feathered), where it carries the
    # prominence's own detail.
    _sm = max(2.0, 4.0 * q)
    a_lo = ndimage.gaussian_filter(a, _sm)
    _sh0 = np.clip(_lum(a[..., None] * p[None, None, :]) / np.maximum(_lum(M), 1e-12), 0, 1)
    _pm = ndimage.gaussian_filter(_sh0, max(1.0, 2.0 * q)) > 0.03
    _pm = ndimage.binary_dilation(_pm, iterations=int(round(4 * q)) or 1)
    _pw = ndimage.gaussian_filter(_pm.astype(np.float64), max(1.5, 3.0 * q))
    a = a_lo + (a - a_lo) * np.clip(_pw, 0.0, 1.0)
    del a_lo, _sh0, _pm, _pw
    prom = a[..., None] * p[None, None, :]
    cor = M - prom
    share = np.clip(_lum(prom) / np.maximum(_lum(M), 1e-12), 0, 1)
    core = ndimage.gaussian_filter(share, 1.0) > 0.4
    lab, n = ndimage.label(core)
    if n:
        sz = ndimage.sum(core, lab, range(1, n + 1))
        core = np.isin(lab, 1 + np.flatnonzero(sz >= 2000 * q * q))
    else:
        core = np.zeros_like(core)
    nfill = 0
    if core.any():
        core = ndimage.binary_dilation(core, iterations=2) & (h > -3)
        cor = np.clip(cor, 0, None)
        lab2, n2 = ndimage.label(core)
        bcy, bcx = cy - y0, cx - x0
        for k in range(1, n2 + 1):
            cor, _wc = _fill_core_polar(cor, lab2 == k, bcy, bcx, h.astype(np.float64), q=q)
        nfill = int(core.sum())
    cor = np.clip(cor, 0, None)
    removed = np.clip(_lum(M) - _lum(cor), 0, None)
    # full frame, chunked through a memmap: the merge itself is never touched
    from numpy.lib.format import open_memmap
    out = open_memmap(os.path.join(wd, CORONA_RGB + ".tmp"), mode="w+", dtype=np.float32,
                      shape=(H, W, 3))
    step = 512
    for r0 in range(0, H, step):
        out[r0:r0 + step] = hdr[r0:r0 + step]
    out[y0:y1, x0:x1] = cor.astype(np.float32)
    out.flush(); del out
    os.replace(os.path.join(wd, CORONA_RGB + ".tmp"), os.path.join(wd, CORONA_RGB))
    lum = np.load(os.path.join(wd, "hdr_lum.npy"))
    lum[y0:y1, x0:x1] = _lum(cor).astype(np.float32)
    np.save(os.path.join(wd, CORONA_LUM), lum)
    np.save(os.path.join(wd, REMOVED_FILE), removed.astype(np.float32))
    frac = float(removed[(h > 0) & (rr < 1.3)].sum() / max(_lum(M)[(h > 0) & (rr < 1.3)].sum(), 1e-30))
    log(f"[ok] H-alpha removed by colour "
        f"(1 : {p[1]:.3f} : {p[2]:.3f}): {100 * frac:.1f}% of inner-corona light; "
        f"{nfill} px core filled")
    return {"p": [float(x) for x in p], "removed_frac": frac, "core_filled_px": nfill,
            "box": [y0, y1, x0, x1]}


def build_light(wd, geo, progress=None):
    """The prominence light: H-alpha unmixed from the short-frame stack (or, for
    a cache without one, the older fast-tier stack, or the merge), with a soft
    noise gate. Writes LIGHT_FILE (linear RGB, merge units, disc box)."""
    log = (lambda m: progress.log(m, None)) if progress is not None else print
    cy, cx, R = float(geo["cy"]), float(geo["cx"]), float(geo["R"])
    hdr = np.load(os.path.join(wd, "hdr_rgb.npy"), mmap_mode="r")
    H, W = hdr.shape[:2]
    y0, y1, x0, x1 = disc_box(cy, cx, R, H, W)
    M = np.asarray(hdr[y0:y1, x0:x1], np.float64)
    q = R / REF_R
    rr, ang = _polar_grid(M.shape[:2], y0, x0, cy, cx, R)
    h = _height(geo, M.shape[:2], y0, x0)
    S, cov, source = None, None, "the merge (no short-frame stack in this cache)"
    for fn, cf, nm, mf in ((SHORT_FILE, SHORT_W, "the short frames", SHORT_META),
                          (STACK_FILE, None, "the fast tiers (older stack; re-stack for the short-frame one)",
                           STACK_META)):
        try:
            if not os.path.exists(os.path.join(wd, fn)):
                continue
            meta = json.load(open(os.path.join(wd, mf)))
            if list(meta["box"]) != [y0, y1, x0, x1]:
                continue
            S = np.load(os.path.join(wd, fn)).astype(np.float64)
            cov = (np.load(os.path.join(wd, cf)) > 0) if cf and os.path.exists(os.path.join(wd, cf)) \
                else np.ones(S.shape[:2], bool)
            source = nm
            break
        except Exception:
            S = None
    if S is None:
        S = M; cov = np.ones(M.shape[:2], bool)
    kS = corona_ratio(S, rr, ang)
    p = np.asarray(prominence_colour(S, kS, rr), np.float64)
    if tuple(p) == tuple(P_DEFAULT):
        p = np.asarray(prominence_colour(M, corona_ratio(M, rr, ang), rr), np.float64)
    a = (S[..., 0] - kS * S[..., 1]) / (1.0 - kS * p[1])
    a = np.where(cov, a, 0.0)
    # the noise of the red excess per height above the limb: below 1.5 sigma it
    # is noise (no rectified haze), above 4 sigma all of it; judged on a 2 px
    # smoothed amplitude so the edge is not decided pixel by pixel
    GSM = 2.0
    hb = np.clip(np.round(h).astype(np.int32) + 20, 0, 4000)
    asm = ndimage.gaussian_filter(a, GSM)
    dev = a - ndimage.gaussian_filter(a, 6.0)
    sig = np.zeros(4001)
    zone = cov & (h > -10) & (h < 1500)
    order = np.argsort(hb[zone], kind="stable"); hs = hb[zone][order]; vs = dev[zone][order]
    edges = np.searchsorted(hs, np.arange(4002))
    for b in range(4001):
        lo_, hi_ = edges[b], edges[b + 1]
        if hi_ - lo_ > 200:
            v = vs[lo_:hi_]
            sig[b] = 1.4826 * np.median(np.abs(v - np.median(v)))
    sig = np.maximum(ndimage.maximum_filter1d(sig, 5), 1e-12)
    if (sig > 1e-12).any():
        sig[sig <= 1e-12] = np.median(sig[sig > 1e-12])
    sg = sig[hb]
    # LOCAL NOISE (2026-10-07, Clifton's 250 mm set): one sigma per height
    # under-rates the colour noise in a bright streamer sector, and 67 noise
    # blobs out to 0.8 R passed the gate as faint prominence light (pink
    # speckles at a low Prominence curve). The gate now uses the larger of the
    # ring's sigma and the local one (mean |deviation|, clipped at 4 ring
    # sigma so prominence structure does not raise its own threshold).
    _ad = np.minimum(np.abs(dev), 4.0 * sg)
    sg = np.maximum(sg, 1.2533 * ndimage.gaussian_filter(np.where(cov, _ad, 0.0), LOCAL_SIG_PX * q)
                    / np.maximum(ndimage.gaussian_filter(cov.astype(np.float64), LOCAL_SIG_PX * q), 1e-6))
    del _ad
    # SATURATED CORES (2026-10-07, Clifton's 250 mm set). Pixels saturated in
    # every short frame have no coverage and came out as zero, so the large
    # prominence printed as an arc round a hole. A hole whose rim is
    # prominence light (median > 3 sigma) is filled at that rim's 99th
    # percentile: a lower bound, like the merge's own fill.
    nfill = 0
    hole = (~cov) & (h > -3) & (h < 0.5 * R)
    if hole.any():
        labh, nh = ndimage.label(hole)
        for k in range(1, nh + 1):
            m = labh == k
            if m.sum() < 4:
                continue
            br = ndimage.binary_dilation(m, iterations=3) & ~hole & cov
            if br.sum() < 5 or np.median(a[br] / sg[br]) < 3.0:
                continue
            a[m] = float(np.percentile(a[br], 99))
            nfill += int(m.sum())
        if nfill:
            asm = ndimage.gaussian_filter(a, GSM)
            log(f"[odd] prominence layer: {nfill} px saturated in every short frame, filled from the rim (lower bound)")
    z = np.clip((asm / (sg / GSM) - 1.5) / 2.5, 0, 1); z = z * z * (3 - 2 * z)
    lab, n = ndimage.label(z > 0.5)
    if n:
        szs = ndimage.sum(np.ones_like(z), lab, range(1, n + 1))
        keep = np.isin(lab, 1 + np.flatnonzero(szs >= 15 * q * q))
        z = z * ndimage.gaussian_filter(ndimage.binary_dilation(keep, iterations=3).astype(float), 1.0)
    # THE FRINGE: bright structure keeps its own pixels; where the light is a few
    # sigma only, it is taken from a smoothed copy, so the edge fades instead of
    # breaking into grain
    snr = a / sg
    wsh = np.clip((snr - 3.0) / 5.0, 0, 1); wsh = wsh * wsh * (3 - 2 * wsh)
    a = wsh * a + (1 - wsh) * ndimage.gaussian_filter(a, 1.5)
    a = np.clip(a, 0, None) * np.clip(ndimage.gaussian_filter(z, 1.5), 0, 1)
    # DETACHED FAINT PIECES (2026-10-07): a piece (light above 0.1% of the
    # prominence light near the limb) that starts more than 0.05 R above the
    # limb and peaks below 5% of it is corona colour noise, not a prominence
    # (Clifton's 250 mm set: 68 such pieces, peaks <= 2.6%; Nico's 600 mm set:
    # 1, peak 0.1%). Above 0.05 R only the kept pieces (and 3 px round them)
    # keep their light, so the sub-0.1% residue cannot print at a low
    # Prominence curve either.
    _sel = (h > 0) & (h < 0.2 * R)
    _ref = float(np.percentile(a[_sel], 99.97)) if _sel.any() else 0.0
    ndrop = 0
    if _ref > 0:
        labd, nd = ndimage.label(a > 0.001 * _ref)
        keep = np.zeros(a.shape, bool)
        if nd:
            _ix = np.arange(1, nd + 1)
            _drop = (np.asarray(ndimage.minimum(h, labd, _ix)) > 0.05 * R) & \
                    (np.asarray(ndimage.maximum(a, labd, _ix)) < 0.05 * _ref)
            ndrop = int(_drop.sum())
            keep = np.isin(labd, _ix[~_drop])
        keep = ndimage.binary_dilation(keep, iterations=3) | (h <= 0.05 * R)
        a = np.where(keep, a, 0.0)
    if ndrop:
        log(f"prominence layer: {ndrop} faint detached piece(s) dropped as noise")
    P = (a[..., None] * p[None, None, :]).astype(np.float32)
    np.save(os.path.join(wd, LIGHT_FILE), P)
    json.dump({"box": [y0, y1, x0, x1], "source": source, "p": [float(x) for x in p]},
              open(os.path.join(wd, LIGHT_META), "w"), indent=1)
    lab2, n2 = ndimage.label(_lum(P) > 0.02 * max(float(np.percentile(_lum(P)[(h > 0) & (h < 0.2 * R)], 99.97)), 1e-30))
    log(f"[ok] prominence layer from {source}: colour 1 : {p[1]:.3f} : {p[2]:.3f}, {n2} piece(s)")
    return {"source": source, "p": [float(x) for x in p], "pieces": int(n2)}


def display_maps(wd, geo):
    """What the renderer needs from the layer, on the disc box: L (layer
    luminance over its own 99.97th percentile near the limb), C (its chroma,
    display-encoded, unit luminance, saturation 1), W1/W2 (weights of the
    prominences and of the chromosphere rim, both faded in over the limb) and
    F (where the corona's Hill detail is eased out under the layer). None when
    the cache has no layer."""
    lp_ = os.path.join(wd, LIGHT_FILE)
    if not (layer_mode() and os.path.exists(lp_) and os.path.exists(os.path.join(wd, CORONA_LUM))):
        return None
    meta = json.load(open(os.path.join(wd, LIGHT_META)))
    y0, y1, x0, x1 = [int(v) for v in meta["box"]]
    R = float(geo["R"]); q = R / REF_R
    P = np.load(lp_).astype(np.float32)
    h = _height(geo, P.shape[:2], y0, x0)
    lp = _lum(P)
    sel = (h > 0) & (h < 0.2 * R)
    Lref = max(float(np.percentile(lp[sel], 99.97)) if sel.any() else float(lp.max()), 1e-30)
    Ln = (lp / Lref).astype(np.float32)
    ch = np.where((lp > 0)[..., None], P / np.maximum(lp, 1e-30)[..., None],
                  np.array(meta.get("p", P_DEFAULT), np.float32)[None, None, :])
    ch = np.clip(ch, 0, None) ** (1 / 2.2)
    ch = ch / np.maximum(_lum(ch), 1e-6)[..., None]
    # prominences (pieces that reach above RIMH px) against the chromosphere rim
    RIMH = 5.0 * q
    on = Ln > 0.02
    lab, n = ndimage.label(on)
    if n:
        reach = ndimage.maximum(np.where(on, h, -1e9), lab, range(1, n + 1))
        promm = np.isin(lab, 1 + np.flatnonzero(np.asarray(reach) > RIMH))
    else:
        promm = np.zeros(on.shape, bool)
    promw = np.clip(ndimage.gaussian_filter(ndimage.binary_dilation(promm, iterations=2)
                                            .astype(np.float32), 1.5) * 1.5, 0, 1)
    band = np.clip((RIMH + 3 - h) / 4.0, 0, 1); band = band * band * (3 - 2 * band)
    rimw = (1 - promw) * band
    lw = np.clip((h + 2) / 4.0, 0, 1); lw = lw * lw * (3 - 2 * lw)
    W1 = lw * (promw + (1 - promw) * (1 - band))
    W2 = lw * rimw
    # the Hill fade, at the default strength: the layer's own alpha, its light
    # down into the limb, and wherever light was taken out of the corona
    y = np.clip(Ln, 0, None) ** 0.5
    y = np.where(y < 0.9, y, 0.9 + 0.1 * np.tanh((y - 0.9) / 0.1))
    alpha = np.clip(y * (W1 + W2) * 1.5, 0, 1)
    deep = np.clip((h + 8) / 4.0, 0, 1) * np.maximum(promw, rimw)
    afd = np.clip(y * deep * 1.5, 0, 1)
    rem = np.zeros_like(Ln)
    try:
        rem = np.load(os.path.join(wd, REMOVED_FILE)).astype(np.float32)
        if rem.shape != Ln.shape:
            rem = np.zeros_like(Ln)
    except Exception:
        pass
    rr_ = max(float(np.percentile(rem[sel], 99.97)) if sel.any() else 0.0, 1e-30)
    remw = np.clip(np.sqrt(np.clip(rem / rr_, 0, None)) * 1.5, 0, 1) * deep
    F = np.clip(ndimage.gaussian_filter(np.maximum(np.maximum(alpha, afd), remw), 2.0) * 2.0, 0, 1)
    return (y0, y1, x0, x1), {"L": Ln, "C": ch.astype(np.float32), "W1": W1.astype(np.float32),
                              "W2": W2.astype(np.float32), "F": F.astype(np.float32)}


def overlay(rgb, maps, P):
    """Lay the layer over a finished (display-encoded) picture, in place on the
    slice it is given. Returns alpha. Same formula as the page's promOverlay."""
    s = float(P.get("promLayer", 1.0))
    if s <= 0:
        return None
    gam = float(P.get("promStretch", 0.5)); rim = float(P.get("promRim", 1.0))
    mode = int(P.get("colourMode", 0) or 0)
    sat = 0.0 if mode == 1 else float(P.get("promSat", 1.3))
    peak = float(np.clip(P.get("promPeak", 1.0), 0.3, 1.0))
    y = np.clip(maps["L"], 0, None) ** np.float32(gam)
    y = np.where(y < 0.9, y, 0.9 + 0.1 * np.tanh((y - 0.9) / 0.1)).astype(np.float32)
    alpha = np.clip(y * (maps["W1"] + rim * maps["W2"]) * np.float32(1.5 * s), 0, 1)
    if mode == 2:
        # mono + tints: one picked colour (display RGB, unit luminance) for all
        h = np.clip(np.array([float(P.get("promHueR", 1.0)), float(P.get("promHueG", 0.70)),
                              float(P.get("promHueB", 0.75))], np.float32), 0.02, 1.0)
        C = (h / np.float32(0.2126 * h[0] + 0.7152 * h[1] + 0.0722 * h[2]))[None, None, :]
    else:
        C = maps["C"]
    ch = np.clip(1 + (C - 1) * np.float32(sat), 0, None)
    pdv = y[..., None] * ch
    pdv /= np.maximum(pdv.max(-1), 1.0)[..., None]
    col = pdv / np.maximum(y, 1e-6)[..., None]
    if peak < 1.0:
        # the brightest channel the layer can lay down, hue kept
        col *= np.minimum(1.0, np.float32(peak) / np.maximum(col.max(-1), 1e-6))[..., None]
    rgb *= (1 - alpha)[..., None]
    rgb += col * alpha[..., None]
    np.clip(rgb, 0, 1, out=rgb)
    return alpha, col


def export_layer2(layers, params, out_path, composite, corona, size="full"):
    """The 0.24 layer as a 16-bit RGBA TIFF (unassociated alpha) on the export's
    grid, made so that it laid over the corona export in Normal mode gives the
    composite back EXACTLY.

    `composite` and `corona` are the two rendered pictures (before orientation),
    with and without the layer. Per pixel, colour L and alpha a solve
    composite = corona * (1 - a) + L * a with L inside 0..1: a starts from the
    layer's own alpha (so L is the prominence's own colour wherever that fits),
    and is raised only where the render's colour runs past 1.0 in a channel --
    the faint red fringe, where alpha-over cannot otherwise reach the composite."""
    import tifffile
    from . import icc
    from .render import apply_orient, defaults_for
    P = dict(defaults_for(getattr(layers, "mode", None))); P.update(params or {})
    H, W = composite.shape[:2]
    y0, y1, x0, x1 = layers.pl_box
    O = np.asarray(composite[y0:y1, x0:x1], np.float32)
    B = np.asarray(corona[y0:y1, x0:x1], np.float32)
    sub = np.zeros_like(O)
    r = overlay(sub, layers.pl_full, P)
    a0 = r[0] if r is not None else np.zeros(O.shape[:2], np.float32)
    up = np.where(O > B, (O - B) / np.maximum(1.0 - B, 1e-6), 0.0)
    dn = np.where(O < B, (B - O) / np.maximum(B, 1e-6), 0.0)
    a = np.clip(np.maximum(a0, np.maximum(up, dn).max(-1)), 0, 1)
    a[np.abs(O - B).max(-1) < 0.5 / 65535] = 0.0
    L = np.clip(B + (O - B) / np.maximum(a, 1e-6)[..., None], 0, 1)
    L[a <= 0] = 0.0
    full = np.zeros((H, W, 4), np.float32)
    full[y0:y1, x0:x1, :3] = L
    full[y0:y1, x0:x1, 3] = a
    full = apply_orient(full, P.get("orient", ""))
    if size == "half":
        h2, w2 = full.shape[0] // 2 * 2, full.shape[1] // 2 * 2
        f = full[:h2, :w2].reshape(h2 // 2, 2, w2 // 2, 2, 4)
        aa = f[..., 3].mean(axis=(1, 3))
        c = (f[..., :3] * f[..., 3:4]).sum(axis=(1, 3)) / np.maximum(
            f[..., 3].sum(axis=(1, 3)), 1e-9)[..., None]
        full = np.concatenate([c, aa[..., None]], -1)
    prof = icc.srgb_profile()
    tifffile.imwrite(out_path, (np.clip(full, 0, 1) * 65535 + 0.5).astype(np.uint16),
                     photometric="rgb", extrasamples=[2], compression="zlib",
                     extratags=[(34675, 1, len(prof), prof, False)],
                     description="eclipseforgehdr 0.24 prominence layer (RGBA, unassociated "
                                 "alpha); same grid as the composite export; over the "
                                 "_corona file in Normal mode it gives the composite back")
    return out_path
