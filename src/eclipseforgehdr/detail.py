"""Detail-extraction layers: photon-noise-adaptive MGN, robust FNRGF,
inner corona (raw + denoised), prominences, earthshine."""
from __future__ import annotations
import os, json
import numpy as np
from scipy import ndimage


from .nafe import nafe_vn

# NAFE-VN defaults. w and gamma are the published working values -- Habbal,
# Druckmuller & Morgan, IWCIA 2014, LNCS 8466, Fig. 3 caption: "gamma = 2.4,
# w = 0.2, sigma = 12" -- rather than values guessed here. sigma is in units of
# the noise sigma_A and is set adaptively per level. The rest measured on the
# reference set:
#  * eps 0.05 (in rank units) recovers most of the near-limb contrast the
#    unrestricted filter loses -- near/far detail ratio 0.22 -> 0.34 -- while
#    a wider window falls back towards no restriction at all.
#  * K 64 with a Gaussian membership 1.5 bins wide keeps the cumulative
#    histogram smooth; a narrower membership prints concentric contour rings.
NAFE_K = 64
NAFE_W = 0.2
NAFE_GAMMA = 2.4
NAFE_EPS = 0.05
NAFE_GRID = 8


def photon_floor(sig_lin, r, r_sky=(0.75, 0.95)):
    """Sigma of log10 luminance from photon noise: C/sqrt(L), calibrated in the
    structure-free far field (fractions of max radius)."""
    sm = np.maximum(ndimage.gaussian_filter(sig_lin, 8), 1.0)
    L = np.log10(np.clip(sig_lin, 1.0, None))
    hp = L - ndimage.gaussian_filter(L, 2.5)
    rmax = r.max()
    sky = (r > r_sky[0] * rmax) & (r < r_sky[1] * rmax)
    if sky.sum() < 10000:
        sky = r > np.percentile(r, 70)
    mad = 1.4826 * np.median(np.abs(hp[sky] - np.median(hp[sky])))
    C = mad * np.sqrt(np.median(sm[sky]))
    return C / np.sqrt(sm)


def radial_profile_map(L, r, valid, smooth=2.0):
    """Azimuthal mean of L as a function of radius, sampled back onto the image.

    Subtracting this before MGN removes the steep near-limb brightness peak,
    which otherwise survives every local-mean subtraction and shows up as the
    hard bright ring hugging the lunar limb.  Sampled with linear (not nearest)
    interpolation in r, so the steep inner part does not staircase into
    concentric rings."""
    n = int(r.max()) + 2
    ridx = r.astype(np.int32)
    w = valid.astype(np.float32)
    cnt = np.bincount(ridx.ravel(), weights=w.ravel(), minlength=n)[:n]
    ssum = np.bincount(ridx.ravel(), weights=(L * w).ravel().astype(np.float64),
                       minlength=n)[:n]
    prof = ssum / np.maximum(cnt, 1e-6)
    good = cnt > 30
    idx = np.arange(n)
    if good.sum() < 4:
        return np.zeros_like(L)
    prof = np.interp(idx, idx[good], prof[good]).astype(np.float32)
    prof = ndimage.gaussian_filter1d(prof, smooth)
    return ndimage.map_coordinates(prof, [np.clip(r.ravel(), 0, n - 1)],
                                   order=1, mode="nearest").reshape(r.shape)


def fourier_background(L, r, cy, cx, r0, order=2, na=360, smooth=6.0):
    """Low-order-in-azimuth radial background mu(r,theta) of log luminance.

    MGN needs its input flattened first, or the corona's own envelope dominates
    the local statistics. An azimuthal MEAN profile is not enough: the corona
    can be several times brighter on one side, so a mean leaves a large residual
    gradient at every azimuth, which inflates the local sigma and crushes the
    fine structure MGN exists to show. Order 2 in azimuth tracks that envelope
    while being far too smooth to absorb streamers or plumes."""
    H, W = L.shape
    rmax = int(np.hypot(max(cy, H - cy), max(cx, W - cx))) + 8
    nr = max(rmax - r0, 2)
    ang = np.linspace(0, 2 * np.pi, na, endpoint=False)
    cols = [np.ones(na)]
    for m in range(1, order + 1):
        cols += [np.cos(m * ang), np.sin(m * ang)]
    Adm = np.stack(cols, 1)
    sa, ca = np.sin(ang), np.cos(ang)
    nc = 2 * order + 1
    mu = np.zeros((nr, na), np.float32)
    cov_r = np.zeros(nr, np.float32)
    for i in range(nr):
        rad = r0 + i
        ys = cy + rad * sa
        xs = cx + rad * ca
        ok = (ys >= 0) & (ys <= H - 1) & (xs >= 0) & (xs <= W - 1)
        nok = int(ok.sum())
        if nok < 24:
            mu[i] = mu[i - 1]
            cov_r[i] = cov_r[i - 1]
            continue
        # damp the harmonics CONTINUOUSLY as coverage falls, instead of dropping
        # the order in integer steps: a step changes mu discontinuously from one
        # ring to the next and paints a hard concentric circle into the output
        cov = nok / float(na)
        cov_r[i] = cov
        lam = 1e-4 / max(cov, 1e-3) ** 4
        rg = np.diag([0.0] + [lam * ((m + 1) // 2) ** 2 for m in range(1, nc)])
        A = Adm[ok]
        v = ndimage.map_coordinates(L, [ys[ok], xs[ok]], order=1)
        w = np.ones(nok)
        for _ in range(3):
            Aw = A * w[:, None]
            coef = np.linalg.solve(Aw.T @ A + rg * nok, Aw.T @ v)
            res = v - A @ coef
            sg = max(1.4826 * np.median(np.abs(res)), 1e-6)
            w = 1.0 / np.maximum(np.abs(res) / (2.0 * sg), 1.0)
        mu[i] = Adm @ coef
    blend = (1.0 - np.clip((cov_r - 0.3) / 0.4, 0, 1))[:, None].astype(np.float32)
    mu = ((1 - blend) * ndimage.gaussian_filter(mu, (smooth, 0))
          + blend * ndimage.gaussian_filter(mu, (8 * smooth, 0)))
    th = np.arctan2(np.arange(H, dtype=np.float32)[:, None] - cy,
                    np.arange(W, dtype=np.float32)[None, :] - cx)
    tidx = (th % (2 * np.pi)) / (2 * np.pi) * na
    ridx = np.clip(r - r0, 0, nr - 1)
    return ndimage.map_coordinates(np.concatenate([mu, mu[:, :1]], 1),
                                   [ridx.ravel(), tidx.ravel()], order=1,
                                   mode="nearest").reshape(H, W)


def mgn(L, floor_map=None, scales=(1.25, 2.5, 5, 10, 20, 40),
        gains=(0.907, 0.976, 0.994, 0.998, 0.999, 1.0), k=0.7, noise_k=2.0,
        global_wt=0.12, global_gamma=3.2, norm_span=None, valid=None):
    """Multiscale Gaussian Normalization (Morgan & Druckmuller 2014).

    Scales and per-scale gains follow the paper: w = 1.25, 2.5, 5, 10, 20, 40,
    and g_i from its Fig. 4 (the mean local standard deviation of pure noise at
    each kernel width, which is what the gains correct for) -- 0.907 at w=1.25
    rising to ~1 by w=5. The earlier hand-set gains over-suppressed the two
    finest scales.

    `valid`: optional mask.  Where given, the per-scale local mean and local
    standard deviation are computed by normalized convolution, so the occulted
    lunar disc contributes nothing to the statistics of the corona just outside
    it (a flat plateau bleeding into the wide kernels used to produce a bright
    halo band around the limb)."""
    if valid is None:
        lo, hi = (np.percentile(L, 0.5), np.percentile(L, 99.95)) if norm_span is None else norm_span
        m = None
    else:
        lo, hi = (np.percentile(L[valid], 0.5), np.percentile(L[valid], 99.95)) \
            if norm_span is None else norm_span
        m = valid.astype(np.float32)
    xn = np.clip((L - lo) / (hi - lo), 0, 1)
    if m is not None:
        xn *= m
    fl = None if floor_map is None else floor_map / (hi - lo)
    acc = np.zeros_like(xn)
    for wsc, g in zip(scales, gains):
        if m is None:
            B = ndimage.gaussian_filter(xn, wsc)
            S = np.sqrt(np.maximum(ndimage.gaussian_filter((xn - B) ** 2, wsc), 1e-12))
        else:
            Vb = np.maximum(ndimage.gaussian_filter(m, wsc), 1e-3)
            B = ndimage.gaussian_filter(xn, wsc) / Vb
            d = (xn - B) * m
            S = np.sqrt(np.maximum(ndimage.gaussian_filter(d * d, wsc) / Vb, 1e-12))
            del d, Vb
        S = np.maximum(S, 0.004)
        if fl is not None:
            S = np.maximum(S, noise_k * fl)
        acc += g * np.arctan(k * (xn - B) / S)
        del B, S
    acc /= sum(gains)
    out = (0.5 + acc / np.pi)
    if global_wt > 0:
        out = global_wt * xn ** (1 / global_gamma) + (1 - global_wt) * out
    if m is not None:
        out = out * m + 0.5 * (1 - m)
    return out


def fnrgf_robust(lum, r, cy, cx, r0, order=6, na=1440, rmax=None):
    """Fourier Normalizing Radial Gradient Filter (Druckmullerova et al. 2011).

    The Fourier order is reduced as a ring leaves the frame: fitting 13 free
    coefficients to a short covered arc oscillates wildly and jumps from ring to
    ring, which is what produced the concentric arc artifacts in the outer
    field. Coverage-matched order plus ridge damping of the higher harmonics
    keeps consecutive rings consistent."""
    H, W = lum.shape
    if rmax is None:
        rmax = int(np.hypot(max(cy, H - cy), max(cx, W - cx))) + 8
    nr = rmax - r0
    ang = np.linspace(0, 2 * np.pi, na, endpoint=False)
    ca, sa = np.cos(ang), np.sin(ang)
    cols = [np.ones(na)]
    for m in range(1, order + 1):
        cols += [np.cos(m * ang), np.sin(m * ang)]
    Adm = np.stack(cols, axis=1)
    # ridge weights: damp higher harmonics, leave the mean free
    ridge = np.zeros(2 * order + 1)
    for m in range(1, order + 1):
        ridge[2 * m - 1] = ridge[2 * m] = 1e-3 * m * m
    L = np.log10(np.clip(lum, 1.0, None))
    mu_g = np.zeros((nr, na), np.float32)
    sd_g = np.zeros((nr, na), np.float32)
    cov_r = np.zeros(nr, np.float32)
    for i in range(nr):
        rad = r0 + i
        ys = cy + rad * sa
        xs = cx + rad * ca
        ok = (ys >= 0) & (ys <= H - 1) & (xs >= 0) & (xs <= W - 1)
        nok = int(ok.sum())
        cov = nok / float(na)
        cov_r[i] = cov
        if nok < 24:
            mu_g[i] = mu_g[i - 1]; sd_g[i] = sd_g[i - 1]; continue
        # order the covered arc can actually support
        o = int(np.clip(np.floor(order * cov * 1.2), 0, order))
        while o > 0 and nok < 12 * (2 * o + 1):
            o -= 1
        nc = 2 * o + 1
        A = Adm[ok][:, :nc]
        Afull = Adm[:, :nc]
        rg = np.diag(ridge[:nc])
        v = ndimage.map_coordinates(L, [ys[ok], xs[ok]], order=1)
        # IRLS with a soft Huber weight rather than hard sigma-clipping: a hard
        # keep/reject mask flips between neighbouring rings as a streamer drifts
        # across the threshold, and each flip steps the fitted background --
        # which is what drew the thin concentric arcs beside bright streamers.
        w = np.ones(nok)
        coef = None
        sg = 1e-6
        res = np.zeros(nok)
        for _ in range(4):
            Aw = A * w[:, None]
            coef = np.linalg.solve(Aw.T @ A + rg * nok, Aw.T @ v)
            res = v - A @ coef
            sg = max(1.4826 * np.median(np.abs(res)), 1e-6)
            u = np.abs(res) / (2.0 * sg)
            w = 1.0 / np.maximum(u, 1.0)
        mu_g[i] = Afull @ coef
        r2 = np.clip(res, -4 * sg, 4 * sg) ** 2
        Aw = A * w[:, None]
        coef2 = np.linalg.solve(Aw.T @ A + rg * nok, Aw.T @ r2)
        var = np.clip(Afull @ coef2, (0.3 * sg) ** 2, None)
        sd_g[i] = np.sqrt(var)
    # Smooth the background model along radius. It has to be smooth by
    # construction -- real radial structure belongs in the residual, not in the
    # model -- so smooth hard, and harder still further out (where the corona
    # varies slowly) and wherever coverage is poor and the fit is shakier.
    rr_ax = (r0 + np.arange(nr)).astype(np.float32)
    far = np.clip((rr_ax - 2.0 * r0) / (2.0 * r0), 0, 1)
    poor = 1.0 - np.clip((cov_r - 0.35) / 0.35, 0, 1)
    blend = np.maximum(far, poor)[:, None].astype(np.float32)
    mu_g = (1 - blend) * ndimage.gaussian_filter(mu_g, (4, 0)) + \
        blend * ndimage.gaussian_filter(mu_g, (30, 0))
    sd_g = (1 - blend) * ndimage.gaussian_filter(sd_g, (10, 0)) + \
        blend * ndimage.gaussian_filter(sd_g, (40, 0))
    theta = np.arctan2(np.arange(H, dtype=np.float32)[:, None] - cy,
                       np.arange(W, dtype=np.float32)[None, :] - cx)
    tidx = (theta % (2 * np.pi)) / (2 * np.pi) * na
    ridx = np.clip(r - r0, 0, nr - 1)
    mu_map = ndimage.map_coordinates(mu_g, [ridx.ravel(), tidx.ravel()],
                                     order=1, mode="nearest").reshape(H, W)
    sd_map = ndimage.map_coordinates(sd_g, [ridx.ravel(), tidx.ravel()],
                                     order=1, mode="nearest").reshape(H, W)
    return (L - mu_map) / np.maximum(sd_map, 1e-4)


def _ss(x):
    """smootherstep-ish ease, for feathering masks without visible edges"""
    return x * x * (3.0 - 2.0 * x)


def limb_radius_map(prof, r_shape, cy, cx, extra=0.0):
    """Sample a per-azimuth limb radius profile onto the image grid."""
    prof = np.asarray(prof, np.float32)
    na = len(prof)
    H, W = r_shape
    th = np.arctan2(np.arange(H, dtype=np.float32)[:, None] - cy,
                    np.arange(W, dtype=np.float32)[None, :] - cx)
    tidx = (th % (2 * np.pi)) / (2 * np.pi) * na
    pr = np.concatenate([prof, prof[:1]])
    out = ndimage.map_coordinates(pr, [tidx.ravel()], order=1,
                                  mode="nearest").reshape(H, W)
    return out + np.float32(extra)


def _fit(a, shape):
    """Crop or edge-pad a 2-D array to exactly `shape`.

    Binning/decimation steps below round the image size down to a multiple of
    the bin factor; sensors whose height or width is not such a multiple (e.g.
    3708 rows, 3708 % 8 == 4) would otherwise hand back an array a few pixels
    short of the layer grid."""
    h, w = shape
    a = a[:h, :w]
    if a.shape[0] < h or a.shape[1] < w:
        a = np.pad(a, ((0, max(0, h - a.shape[0])), (0, max(0, w - a.shape[1]))),
                   mode="edge")
    return a


def _deband(layer, r, valid, cy=None, cx=None, r0=None, order=6):
    """Remove residual radial trend from a detail layer, allowing the trend to
    vary slowly AROUND the disc as well as with radius.

    A normalized detail layer should have no radial trend, so whatever survives
    is filter residue -- most visibly a rim just outside the limb, where the
    radial pre-flattening is least accurate because the profile is steepest.

    This used to subtract the azimuthal MEAN profile, which removes a perfect
    ring and nothing else. The rim is not a perfect ring: it is produced by the
    local radial gradient, so it is strongest where the corona is brightest and
    weakest on the dark sides. Measured on the reference layer, the mean-only
    version left the rim's variation around the disc completely untouched
    (0.0140, identical before and after) while an azimuthal fit removes it.

    The ORDER is set by where the two things live in azimuth, which is measured
    rather than guessed. On the reference layer the rim is a low-order pattern
    -- m=1 and m=2 alone carry 63% of it, m<=6 carries 82% -- while coronal
    structure at 2R peaks at m=21 (17-degree features) with 80% of its power
    above m=6. They barely overlap, so:

        order  rim removed  streamer power kept
          2       31%             100%
          6       44%             100%
          8       46%              98%
         12       48%              90%

    Order 6 is where the rim stops giving way for free. Order 2 (used in 0.9.1)
    left more than a third of the rim behind, which showed as a halo surviving
    on some sides of the disc and not others -- exactly the m=1/m=2 asymmetry.

    Note this is NOT the order used to pre-flatten the image before MGN: that
    stays at 2, because its job is to remove the corona's envelope without
    touching the structure MGN then has to find.
    """
    base = float(np.mean(layer[valid]))
    if cy is None or cx is None or r0 is None or order < 1:
        trend = radial_profile_map(layer, r, valid)
    else:
        trend = fourier_background(layer, r, cy, cx, int(r0), order=order,
                                   smooth=2.0)
    out = layer - trend + base
    return np.where(valid, out, 0.5).astype(np.float32)


def _soft_norm(x, mask, p_lo=0.5, p_hi=99.7, gain=1.6):
    lo, hi = np.percentile(x[mask], p_lo), np.percentile(x[mask], p_hi)
    return 0.5 + 0.5 * np.tanh(gain * ((x - lo) / max(hi - lo, 1e-6) - 0.5))


def build_layers(wd, progress, denoise="fine", earthshine=False):
    if denoise is True:
        denoise = "fine"
    if denoise is False:
        denoise = "off"
    lstats = {}
    ks = DENOISE_PROFILES.get(denoise, DENOISE_PROFILES["fine"])
    do_dn = any(k > 0 for k in ks)
    geo = json.load(open(os.path.join(wd, "geometry.json")))
    cy, cx, R = geo["cy"], geo["cx"], geo["R"]
    margin = float(geo.get("limb_margin", geo.get("Rmask", R + 4.0) - R))
    prof = geo.get("limb_prof")
    lum = np.load(os.path.join(wd, "hdr_lum.npy"))
    H, W = lum.shape
    yy = np.arange(H, dtype=np.float32)[:, None] - cy
    xx = np.arange(W, dtype=np.float32)[None, :] - cx
    r = np.sqrt(yy * yy + xx * xx)
    disc = r < R - 6
    # per-azimuth mask radius (falls back to a circle for old caches)
    if prof:
        Rmap = limb_radius_map(prof, (H, W), cy, cx, margin)
    else:
        Rmap = np.float32(R + margin)
    disc_m = r < Rmap

    nf = photon_floor(lum, r)
    L = np.log10(np.clip(lum, 1.0, None))
    if do_dn:
        progress.log(f"denoising HDR master (multiscale, profile: {denoise})...", 0.935)
        Ldn = denoise_loglum(L, nf, ks=ks)
        lum_dn = (10.0 ** Ldn).astype(np.float32)
    else:
        Ldn = L
        lum_dn = lum

    progress.log("MGN detail extraction...", 0.94)
    # 1) mask the disc out of the statistics entirely (normalized convolution)
    # 2) subtract the azimuthal radial profile first, so the real brightness
    #    peak at the limb is not re-sharpened into a hard ring
    # 3) normalise the flattened residual on the SAME span as the raw log
    #    luminance, so contrast (and the mgnContrast slider) behaves as before
    valid = r > Rmap
    Lf = (Ldn - fourier_background(Ldn, r, cy, cx, int(R) + 3)).astype(np.float32)
    inner_field = valid & (r < 4 * R)
    sd = float(np.std(Lf[inner_field])) if inner_field.any() else 0.05
    half = max(6.0 * sd, 1e-3)          # gain matched to the residual's own scale
    mgl = mgn(Lf, floor_map=nf, valid=valid, norm_span=(-half, half))
    mgl = _deband(mgl, r, valid, cy, cx, R + margin)
    np.save(os.path.join(wd, "mgn.npy"), mgl.astype(np.float32))
    del nf, L, Lf, mgl

    progress.log("FNRGF detail extraction...", 0.955)
    D = fnrgf_robust(lum_dn, r, cy, cx, int(R) + 4)
    np.save(os.path.join(wd, "fnrgf.npy"), D.astype(np.float32))
    del D

    # --- NAFE with a variable neighbourhood -------------------------------
    # The other two detail layers both need to be told where the Moon is: MGN
    # via `valid`, FNRGF via its radial fit. NAFE-VN needs no geometry at all.
    # Each pixel is ranked against its own neighbourhood, restricted to
    # neighbours of similar brightness, so the dark lunar plateau drops out of
    # the corona's statistics because it is dark, not because a circle was
    # drawn around it. Where the limb fit is imperfect this layer is unaffected.
    progress.log("NAFE (variable neighbourhood)...", 0.958)
    try:
        nv = nafe_vn(Ldn, K=NAFE_K, w=NAFE_W, gamma=NAFE_GAMMA,
                     eps_frac=NAFE_EPS, n_scales=8, grid=NAFE_GRID)
        np.save(os.path.join(wd, "nafe.npy"), nv.astype(np.float32))
        lstats["nafe"] = {"K": NAFE_K, "w": NAFE_W, "eps": NAFE_EPS}
        del nv
    except Exception as e:
        progress.log(f"NAFE layer unavailable ({e})", None)
        np.save(os.path.join(wd, "nafe.npy"), np.full((H, W), 0.5, np.float32))
    del lum_dn, Ldn

    progress.log("inner corona layers...", 0.97)
    shortL = np.load(os.path.join(wd, "short_lum.npy"))

    # --- this layer gets its OWN lunar geometry ---
    # The inner stack is built from the four SHORTEST tiers, which were shot in
    # the first seconds of the bracket; the merged image spans the whole of it.
    # The Moon is therefore in a different place in the two, by ~20 px on the
    # reference set -- the pipeline already measures exactly this offset for the
    # prominence tier and reports it. Masking this layer with the MERGED disc
    # puts the mask off-centre, and an off-centre disc mask prints as a bright
    # arc on one side of the limb and a dark arc on the other. That pair of arcs
    # was the ring in the composite.
    #
    # Each layer is masked by the disc IT sees. The composite then masks by the
    # merged disc as well, so the two exclusions union naturally and no layer
    # contributes light from where its own Moon was.
    cys, cxs, Rs = cy, cx, R
    try:
        _ig = geo.get("inner_geom")
        if _ig:
            cys, cxs, Rs = float(_ig["cy"]), float(_ig["cx"]), float(_ig["R"])
            lstats["inner_geom"] = {"cy": cys, "cx": cxs, "R": Rs,
                                    "offset_px": float(np.hypot(cys - cy, cxs - cx))}
            progress.log(f"inner-stack lunar disc (from the track): "
                         f"({cys:.0f},{cxs:.0f}) R={Rs:.0f}px — "
                         f"{np.hypot(cys - cy, cxs - cx):.0f}px from the merged limb",
                         None)
    except Exception as e:
        progress.log(f"inner-stack geometry unavailable ({e}); using merged", None)
    _yi = np.arange(H, dtype=np.float32)[:, None] - cys
    _xi = np.arange(W, dtype=np.float32)[None, :] - cxs
    r_s = np.sqrt(_yi * _yi + _xi * _xi)
    valid_s = r_s > (Rs + margin)
    inner_field_s = valid_s & (r_s < 4 * Rs)

    nfs = photon_floor(shortL, r_s)
    Ls = np.log10(np.clip(shortL, 1.0, None))
    ann = (r_s > Rs + 2) & (r_s < 1.5 * Rs)
    Lsf = (Ls - fourier_background(Ls, r_s, cys, cxs, int(Rs) + 3)).astype(np.float32)
    shalf = max(6.0 * float(np.std(Lsf[inner_field_s])), 1e-3) \
        if inner_field_s.any() else 0.05
    common = dict(scales=(2.5, 5, 10, 20, 40, 80), k=0.8, global_wt=0.0,
                  valid=valid_s, norm_span=(-shalf, shalf))
    raw = _deband(mgn(Lsf, floor_map=None, gains=(1, 1, 1, 1, 1, 1), **common),
              r_s, valid_s, cys, cxs, Rs + margin)
    np.save(os.path.join(wd, "inner0.npy"), _soft_norm(raw, ann).astype(np.float32))
    Ls_dn = denoise_loglum(Lsf, nfs, ks=tuple(k * 0.85 for k in ks)) if do_dn else Lsf
    clean = _deband(mgn(Ls_dn, floor_map=nfs, gains=(0.8, 0.95, 1, 1, 1, 1), **common),
                r_s, valid_s, cys, cxs, Rs + margin)
    clean = _soft_norm(clean, ann)
    np.save(os.path.join(wd, "inner.npy"), clean.astype(np.float32))
    del nfs, raw, Ls_dn

    # prominence gate v4: Halpha COLOUR detection in a single fast tier.
    # Prominences are deep red; chromosphere glare and corona are not.
    prgb_path = os.path.join(wd, "prom_rgb.npy")
    gate = np.zeros((H, W), np.float32)
    if os.path.exists(prgb_path):
        prgb = np.load(prgb_path).astype(np.float32)   # half-res HxWx3
        if not np.isfinite(prgb).all():
            progress.log("warning: prominence stack has non-finite samples "
                         "(stale cache?) — sanitising", None)
            prgb = np.nan_to_num(prgb, nan=0.0, posinf=0.0, neginf=0.0)
        Rc = ndimage.gaussian_filter(prgb[:, :, 0], 2)
        GB = ndimage.gaussian_filter(0.5 * (prgb[:, :, 1] + prgb[:, :, 2]), 2)
        redness = Rc / np.maximum(GB, 1e-3)
        h2, w2 = redness.shape
        yyh = np.arange(h2, dtype=np.float32)[:, None] - cy / 2
        xxh = np.arange(w2, dtype=np.float32)[None, :] - cx / 2
        # the prominence stack is a single tier, whose Moon sits where THAT
        # tier's Moon sits; use its own limb, not the merged one
        pg = geo.get("prom_geom")
        if pg:
            pcy, pcx = pg["cy"], pg["cx"]
            yyh = np.arange(h2, dtype=np.float32)[:, None] - pcy
            xxh = np.arange(w2, dtype=np.float32)[None, :] - pcx
            rh = np.sqrt(yyh * yyh + xxh * xxh)
            Rh = limb_radius_map(np.asarray(pg["prof"], np.float32),
                                 (h2, w2), pcy, pcx)
        else:
            rh = np.sqrt(yyh * yyh + xxh * xxh)
            Rh = (limb_radius_map(np.asarray(prof, np.float32) / 2.0,
                                  (h2, w2), cy / 2, cx / 2)
                  if prof else np.float32(R / 2))
        # reference redness = the corona's own colour in the ring around the
        # limb, so the gate is white-balance independent; threshold from a
        # robust spread rather than a fixed multiplier
        # The gate window is a fraction of the disc, not a pixel count. The
        # 4/60/8/6/70/25 half-res px it used are 0.013/0.194/0.026/0.019/0.226/
        # 0.081 of the reference lunar radius (310 half-res px), and these
        # coefficients reproduce those numbers to within 0.15 px on that set.
        # Left absolute, on a 150 px disc the same window reaches 0.4 R above
        # the limb: the gate then covers inner-corona loops and promGain
        # brightens them as if they were prominences, while the widened
        # reference annulus inflates its own MAD and pushes the threshold past
        # anything real.
        _Rs = float(np.median(np.asarray(Rh, np.float32)))
        annh = (rh > Rh - 0.013 * _Rs) & (rh < Rh + 0.194 * _Rs)
        lumh = 0.2126 * prgb[:, :, 0] + 0.7152 * prgb[:, :, 1] + 0.0722 * prgb[:, :, 2]
        sel = redness[annh & (lumh > np.percentile(lumh[annh], 20))]
        if sel.size > 1000:
            med = float(np.median(sel))
            mad = 1.4826 * float(np.median(np.abs(sel - med))) + 1e-4
            t0 = max(med + 3.0 * mad, 1.22 * med)
            t1 = max(med + 8.0 * mad, 1.70 * med)
            progress.log(f"prominence colour: corona R/GB {med:.2f}, "
                         f"threshold {t0:.2f}-{t1:.2f}", None)
            lstats["prom"] = {"med_red": med, "t0": t0, "t1": t1}
            g = np.clip((redness - t0) / (t1 - t0), 0, 1)
            g *= _ss(np.clip((rh - (Rh - 0.026 * _Rs)) / (0.019 * _Rs), 0, 1))
            g *= _ss(np.clip(((Rh + 0.226 * _Rs) - rh) / (0.081 * _Rs), 0, 1))
            g = ndimage.gaussian_filter(g, 2)
            gate = _fit(np.clip(g * 1.6, 0, 1).repeat(2, 0).repeat(2, 1), (H, W))
            gate = ndimage.gaussian_filter(gate, 2)
        else:
            progress.log("prominence colour: not enough limb samples", None)
    if "prom" in lstats:
        lstats["prom"]["area_px"] = int((gate > 0.3).sum())
    np.save(os.path.join(wd, "prom.npy"), gate.astype(np.float32))
    del gate, Ls

    progress.log("Pellett layer...", 0.98)
    _pellett(wd, lum, r, cy, cx, R, disc_m)

    ep = os.path.join(wd, "earth.npy")
    if earthshine:
        progress.log("earthshine layer...", 0.99)
        _earthshine(wd, r, cy, cx, R)
    elif os.path.exists(ep):
        os.remove(ep)
    return lstats


def _pellett(wd, lum, r, cy, cx, R, disc, blur_deg=6.0, na=2880):
    """Tangential unsharp (Pellett-style): subtract a rotational blur about the
    disc center from the log luminance. Softer texture than MGN/FNRGF."""
    from skimage.transform import warp_polar
    H, W = lum.shape
    L = np.log10(np.clip(lum, 1.0, None))
    L = L.copy(); L[disc] = np.median(L[disc])
    rmax = int(np.hypot(max(cy, H - cy), max(cx, W - cx))) + 8
    P = warp_polar(L, center=(cy, cx), radius=rmax, output_shape=(na, rmax), order=1)
    V = warp_polar(np.ones_like(L), center=(cy, cx), radius=rmax,
                   output_shape=(na, rmax), order=1)
    sigma_bins = blur_deg / 360.0 * na
    # normalized convolution: out-of-frame samples don't drag the blur down
    Pb = (ndimage.gaussian_filter1d(P * V, sigma_bins, axis=0, mode="wrap") /
          np.maximum(ndimage.gaussian_filter1d(V, sigma_bins, axis=0, mode="wrap"), 1e-3))
    theta = np.arctan2(np.arange(H, dtype=np.float32)[:, None] - cy,
                       np.arange(W, dtype=np.float32)[None, :] - cx)
    tidx = (theta % (2 * np.pi)) / (2 * np.pi) * na
    ridx = np.clip(r, 0, rmax - 1)
    blur = ndimage.map_coordinates(Pb, [tidx.ravel(), ridx.ravel()],
                                   order=1, mode="nearest").reshape(H, W)
    res = L - blur
    outer = r > R + 10
    s = 1.4826 * np.median(np.abs(res[outer] - np.median(res[outer]))) if outer.any() else 0.01
    pel = 0.5 + 0.5 * np.tanh(res / (4.0 * max(s, 1e-4)))
    pel[disc] = 0.5
    np.save(os.path.join(wd, "pellett.npy"), pel.astype(np.float32))


def _earthshine(wd, r, cy, cx, R):
    """Glare-model-subtracted, heavily denoised lunar disc from the longest tiers."""
    H, W = r.shape
    p = os.path.join(wd, "long_lum.npy")
    if not os.path.exists(p):
        np.save(os.path.join(wd, "earth.npy"), np.full((H, W), 0.5, np.float32))
        return
    longL = np.load(p)
    Ll = ndimage.gaussian_filter(np.log10(np.clip(longL, 1.0, None)), 6)
    # low-order Fourier glare model on rings inside the disc
    order, na = 3, 720
    ang = np.linspace(0, 2 * np.pi, na, endpoint=False)
    cols = [np.ones(na)]
    for m in range(1, order + 1):
        cols += [np.cos(m * ang), np.sin(m * ang)]
    Adm = np.stack(cols, axis=1)
    nr = int(R) - 8
    mu_g = np.zeros((nr, na), np.float32)
    sa, ca = np.sin(ang), np.cos(ang)
    for i in range(nr):
        ys = np.clip(cy + i * sa, 0, H - 1); xs = np.clip(cx + i * ca, 0, W - 1)
        v = ndimage.map_coordinates(Ll, [ys, xs], order=1)
        keep = np.ones(na, bool)
        for _ in range(2):
            coef, *_ = np.linalg.lstsq(Adm[keep], v[keep], rcond=None)
            res = v - Adm @ coef
            s = 1.4826 * np.median(np.abs(res[keep]))
            keep = np.abs(res) < 2.5 * max(s, 1e-6)
        mu_g[i] = Adm @ coef
    mu_g = ndimage.gaussian_filter(mu_g, (4, 0))
    theta = np.arctan2(np.arange(H, dtype=np.float32)[:, None] - cy,
                       np.arange(W, dtype=np.float32)[None, :] - cx)
    tidx = (theta % (2 * np.pi)) / (2 * np.pi) * na
    ridx = np.clip(r, 0, nr - 1)
    model = ndimage.map_coordinates(mu_g, [ridx.ravel(), tidx.ravel()],
                                    order=1, mode="nearest").reshape(H, W)
    res = np.where(r < R - 12, Ll - model, 0).astype(np.float32)
    # denoise at coarse scale (earthshine features are large), then normalize
    b = res[: H // 8 * 8, : W // 8 * 8].reshape(H // 8, 8, W // 8, 8).mean(axis=(1, 3))
    b = ndimage.gaussian_filter(b, 1.6)
    up = _fit(ndimage.zoom(b, 8, order=1), (H, W))
    innerm = r < R - 40
    s = up[innerm].std() if innerm.any() else 1.0
    E = 0.5 + 0.5 * np.tanh(up / (2.5 * max(s, 1e-6)))
    feather = np.clip(((R - 25) - r) / 25, 0, 1)
    np.save(os.path.join(wd, "earth.npy"),
            (E * feather + 0.5 * (1 - feather)).astype(np.float32))


DENOISE_PROFILES = {
    "off":    (0.0, 0.0, 0.0, 0.0),
    "fine":   (1.5, 1.0, 0.0, 0.0),
    "medium": (1.8, 1.5, 0.8, 0.0),
    "strong": (2.2, 1.8, 1.2, 0.6),
}


def denoise_loglum(L, noise_map, ks=(1.5, 1.0, 0.0, 0.0), levels=4):
    """A-trous multiscale soft-threshold denoise of log-luminance.
    ks: per-scale-level threshold strengths (finest first; 0 = level untouched).
    noise_map: per-pixel sigma of L (photon model). Dependency-free starlet-style."""
    rng = np.random.default_rng(12345)
    probe = rng.standard_normal((512, 512)).astype(np.float32)
    sigmas = [1.0 * 2 ** i for i in range(levels)]
    # per-level response of unit white noise
    resp = []
    gp_prev = probe
    for s in sigmas:
        gp = ndimage.gaussian_filter(probe, s)
        resp.append(float(np.std(gp_prev - gp)))
        gp_prev = gp
    g_prev = L
    details = []
    for s in sigmas:
        g = ndimage.gaussian_filter(L, s)
        details.append(g_prev - g)
        g_prev = g
    base = g_prev
    rec = base
    for lvl, (d, rsp) in enumerate(zip(details, resp)):
        klvl = ks[lvl] if lvl < len(ks) else 0.0
        if klvl > 0:
            thr = klvl * rsp * noise_map
            d = np.sign(d) * np.maximum(np.abs(d) - thr, 0)
        rec = rec + d
    return rec
