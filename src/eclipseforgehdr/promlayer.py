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
        log(f"prominence stack: tier {s:g}s — {note}")
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
    log(f"prominence stack: {len(tiers)} tier(s), PSF sigma {sig:.2f} px from the limb, "
        f"{RL_ITERS} Richardson-Lucy iterations")


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
    log(f"prominence layer: {nprom} prominence(s), "
        f"{'sharp stack' if have_stack else 'merge only (no prominence stack in this cache)'}, "
        f"colour G/R {p[1]:.3f} B/R {p[2]:.3f}")
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

