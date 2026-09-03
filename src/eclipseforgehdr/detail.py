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
NAFE_K = 128          # histogram levels; 64 measurably cost corona contrast
NAFE_W = 0.2          # the paper's eq.2 weight; the stored layer is E, so this
                      # is unused here -- nafeMix in render.py is the live w
NAFE_GAMMA = 2.4
NAFE_EPS = 0.10       # value window, in rank units. Swept on the reference set:
                      # 0.02 prints concentric rings (the paper's "fragmentation"),
                      # 0.05 gave high-pass 0.076, 0.10 gives 0.109 at the same
                      # agreement with FNRGF (0.711 -> 0.721), and past 0.15 the
                      # curve is flat and the extra is noise.
NAFE_NOISE_MULT = 4.0 # the paper's sigma in units of sigma_A; their range is 2..12
NAFE_NEIGH_R = 0.13   # neighbourhood sigma as a fraction of the lunar radius.
                      # Measured flat from 0.06 R to 0.51 R, so this is not
                      # critical -- but it has to scale with the disc, or the
                      # neighbourhood means something different on every camera.
NAFE_GRID = 8
# NAFE is the one detail layer that takes no geometry -- it ranks each pixel
# against a neighbourhood restricted in VALUE, so a bad limb fit cannot hurt it.
# That also meant it was the only layer whose input still carried the corona's
# whole radial envelope, and it paid for that twice: 60% of its output range
# went on a large-scale gradient rather than on structure, and the steep falloff
# just outside the limb drove the equalisation into a dark ring.
#
# The envelope is removed with a plain Gaussian high-pass, which needs no circle
# and no limb. Measured on the reference stack against both the shipped version
# and the fitted radial profile MGN uses (which does need the limb, and puts a
# wedge on the disc where it extrapolates inward):
#
#                      large-scale   detail 1.05-1.5 R   ripple outside the mask
#   as shipped            68.3%          0.0611                0.3451
#   fitted radial profile  5.4%          0.1328                0.2543
#   Gaussian, 0.08 R       3.6%          0.1360                0.1210
#
# So the geometry-free option is also the better one on every count: 2.2x the
# near-limb detail and 2.9x less ripple where the ring used to be. Sigma was
# swept over 0.08-0.35 R; smaller is better near the limb and the far field is
# flat across the range. Tied to R, so it means the same on any focal length.
#
# ...but that ring was reduced, not removed, and 0.14.4 got a render showing a
# thick black annulus hugging the disc where FNRGF was clean. The measurement
# above missed it because "ripple" was an RMS, and a broad smooth depression
# barely moves an RMS. It is plainly visible as the MEAN of the layer per radial
# ring, which is how it is measured now.
#
# The cause is the kernel, not the idea. A symmetric Gaussian mean is a bad
# background estimate beside a huge dark hole: near the limb the mean is dragged
# down by the disc, so L - mean overshoots bright at the limb and undershoots
# dark just outside it, ~2 sigma wide -- an unsharp halo, 100 px of black on a
# 622 px disc. The fix is to build the mean by NORMALIZED CONVOLUTION over
# non-disc pixels only, exactly as MGN already does. Measured on a synthetic
# corona carrying a KNOWN fine modulation, so "detail" is a correlation with the
# truth rather than a variance that the artifact itself inflates:
#
#                       ring depth   correlation with the true structure
#   plain Gaussian        0.1684                  0.468
#   normalized conv.      0.0713                  0.814
#
# Two other candidates did better still on a perfect limb fit -- subtracting the
# Fourier radial background reached 0.0047 and 0.938 -- but they need a centre
# and a radius, and with the centre 0.05 R out that option is WORSE than the
# plain Gaussian (0.2267 / 0.542): it manufactures its own artifact. Normalized
# convolution needs only the disc mask, which the renderer already uses to fill
# the disc, so it adds no dependency that is not already there -- and it
# degrades gently: with the centre 0.10 R out it still scores 0.1135 / 0.731,
# better than the plain Gaussian with a perfect fit.
#
# Sigma is left at 0.08 R deliberately. Re-swept with the disc excluded, 0.15 R
# measured a little better (0.0576 / 0.863), but the sigma choice interacts with
# the real scale of coronal structure, which a synthetic of azimuthal cosines
# does not faithfully represent; the kernel correction does not. One change,
# and the sweep is worth repeating on real data.
NAFE_FLATTEN_R = 0.08

# --- progress weighting for the detail stage -------------------------------
# The stage used to split its band of the progress bar evenly by step, and the
# split is nothing like even. Measured wall time (timedetail.py, synthetic
# corona frames, one pass each):
#
#                              10.8 Mpx        43.4 Mpx
#   denoise HDR master            2.9s            18.2s
#   MGN                          26.9s           197.9s
#   FNRGF                         7.6s            22.9s
#   NAFE                         11.1s            56.3s
#   inner corona (all of it)    135.3s          1001.1s
#   prominence colour             0.9s             3.6s
#   Pellett                      11.1s            31.7s
#   ------------------------------------------------------
#   total                       198.7s          1343.9s
#
# Two things follow. The inner-corona block is two thirds to three quarters of
# the stage on its own -- it runs the multiscale normalisation TWICE at full
# resolution, once raw and once denoised (61.0s and 61.7s of its 135.3s at
# 10.8 Mpx) -- so the bar stood still through most of the stage, which is the
# "almost full and needs a long time" of the first field report. And the whole
# stage does not scale with pixel count: 4x the pixels cost 6.8x the time, so
# on a 45 Mpx body it is 22 minutes, not the couple of minutes 6.5% of the bar
# implies. See _BAR_DETAIL in pipeline.py for the band itself.
#
# The weights below are the 43.4 Mpx column (nearest to the cameras in use),
# with the inner block split by its own sub-step measurement. Only ratios
# matter. The MGN-to-inner ratio is 0.199 at 10.8 Mpx and 0.198 at 43.4, so the
# dominant pair is size-stable; the small layers drift a little and are set
# from the large frame. Earthshine is not measured (it was off in both runs)
# and is a placeholder. Every run that finishes now prints its own timing
# summary, which is how these get corrected against real cameras.
_DETAIL_BAND = (0.935, 1.0)
_DETAIL_W = {                 # in pipeline order -- dict order is the order
    "denoise":    18.0,
    "mgn":       198.0,
    "fnrgf":      23.0,
    "nafe":       56.0,
    "inner_bg":   68.0,       # geometry, photon floor, Fourier background
    "inner_raw": 452.0,       # first multiscale pass
    "inner_dn":   25.0,       # denoise of the short stack
    "inner_cln": 457.0,       # second multiscale pass
    "prom":        4.0,
    "promdet":     3.0,
    "pellett":    32.0,
    "earth":      20.0,       # not measured -- placeholder
}


def _detail_fracs():
    lo, hi = _DETAIL_BAND
    tot = float(sum(_DETAIL_W.values()))
    out, c = {}, 0.0
    for k, v in _DETAIL_W.items():
        out[k] = lo + (hi - lo) * c / tot
        c += v
    return out


_DF = _detail_fracs()


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


def resolution_floor(L, r, R, lo=1.3, hi=3.0):
    """The smallest filter scale this image actually resolves, in pixels.

    MGN's finest scales were fixed at 1.25 and 2.5 px. That is a statement about
    a sensor and a lens, not about eclipses: at 1.79 arcsec/px behind 600 mm
    those scales sit near the optical limit, but at 3.19 arcsec/px behind a
    240 mm consumer zoom -- whose real resolution is 3-4 px -- they sit entirely
    below anything the optics delivered, so they can only amplify noise. On that
    dataset the pixel-scale band carried 9.6% of the local mean against 1.9-2.7%
    for the bands where the real structure lives.

    Measured rather than assumed: band-pass sd falls as ~1/s for white noise and
    more slowly for real structure, so the scale where the log-log slope stops
    looking like noise is where the data starts. Returns a floor, never a
    ceiling -- it can only push the ladder coarser than it already is.
    """
    m = (r > lo * R) & (r < hi * R)
    if int(m.sum()) < 20000:
        return 1.25
    ss = np.array([1.0, 1.4, 2.0, 2.8, 4.0, 5.6])
    prev = L
    sd = []
    for s in ss:
        cur = ndimage.gaussian_filter(L, float(s))
        sd.append(float((prev - cur)[m].std()))
        prev = cur
    sd = np.maximum(np.array(sd), 1e-12)
    # local log-log slope; white noise gives about -1, structure is shallower
    sl = np.diff(np.log(sd)) / np.diff(np.log(ss))
    for i, v in enumerate(sl):
        if v > -0.7:
            return float(ss[i])
    return float(ss[-1])


def scale_ladder(top, floor, n=6):
    """n log-spaced scales up to `top`, none finer than `floor`."""
    top = float(max(top, floor * 2.0))
    out, s = [], top
    for _ in range(n):
        out.append(s)
        s /= 2.0
    out = sorted(max(float(x), float(floor)) for x in out)
    # collapse the ones the floor merged together
    ded = [out[0]]
    for x in out[1:]:
        if x > ded[-1] * 1.15:
            ded.append(x)
    return tuple(ded)


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
        progress.log(f"denoising HDR master (multiscale, profile: {denoise})...", _DF["denoise"])
        Ldn = denoise_loglum(L, nf, ks=ks)
        lum_dn = (10.0 ** Ldn).astype(np.float32)
    else:
        Ldn = L
        lum_dn = lum

    progress.log("MGN detail extraction...", _DF["mgn"])
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
    # The ladder: top tied to R so it covers the same coronal structure at any
    # focal length (0.0643 * 622 = 40 px reproduces the reference set exactly),
    # bottom raised to whatever this image actually resolves.
    _fl = resolution_floor(Lf, r, R)
    _sc = scale_ladder(0.0643 * R, _fl)
    lstats["mgn_scales"] = {"px": [round(x, 2) for x in _sc],
                            "resolution_floor_px": round(_fl, 2)}
    progress.log(f"MGN scales {', '.join('%.1f' % x for x in _sc)} px "
                 f"(this image resolves down to {_fl:.1f} px)", None)
    mgl = mgn(Lf, floor_map=nf, valid=valid, norm_span=(-half, half), scales=_sc)
    mgl = _deband(mgl, r, valid, cy, cx, R + margin)
    np.save(os.path.join(wd, "mgn.npy"), mgl.astype(np.float32))
    del nf, L, Lf, mgl

    progress.log("FNRGF detail extraction...", _DF["fnrgf"])
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
    progress.log("NAFE (variable neighbourhood)...", _DF["nafe"])
    try:
        # combine=False: store E, the equalized field, NOT the paper's eq. 2
        # output B = (1-w) T_gamma + w E. B is their final display image and is
        # four fifths gamma transform at w = 0.2 -- as a detail layer it was a
        # second copy of the base image and diluted MGN and FNRGF with it. The
        # eq. 2 mix happens in render.py instead, where the composite envelope
        # is T_gamma and the nafeMix slider is w.
        # Flatten the envelope first -- see NAFE_FLATTEN_R. The local mean is
        # built by normalized convolution over non-disc pixels only: a plain
        # Gaussian straddling the dark disc is what produced the black annulus.
        _s = max(NAFE_FLATTEN_R * R, 4.0)
        _w = (~disc_m).astype(np.float32)
        _num = ndimage.gaussian_filter(Ldn * _w, _s)
        _den = ndimage.gaussian_filter(_w, _s)
        del _w
        np.maximum(_den, 1e-6, out=_den)
        _num /= _den
        del _den
        _Lnf = (Ldn - _num).astype(np.float32)
        del _num
        nv = nafe_vn(_Lnf, K=NAFE_K, gamma=NAFE_GAMMA, combine=False,
                     sigma_sp=max(NAFE_NEIGH_R * R, 8.0) / NAFE_GRID,
                     noise_mult=NAFE_NOISE_MULT,
                     eps_frac=NAFE_EPS, kernel="gauss", grid=NAFE_GRID)
        np.save(os.path.join(wd, "nafe.npy"), nv.astype(np.float32))
        del _Lnf
        lstats["nafe"] = {"K": NAFE_K, "eps": NAFE_EPS, "layer": "E",
                          "flatten_px": round(max(NAFE_FLATTEN_R * R, 4.0), 1),
                          "flatten": "normalized convolution, disc excluded",
                          "noise_mult": NAFE_NOISE_MULT,
                          "neigh_px": round(NAFE_NEIGH_R * R, 1),
                          }
        del nv
    except Exception as e:
        progress.log(f"NAFE layer unavailable ({e})", None)
        np.save(os.path.join(wd, "nafe.npy"), np.full((H, W), 0.5, np.float32))
    del lum_dn, Ldn

    progress.log("inner corona layers...", _DF["inner_bg"])
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
    # same rule, one octave coarser -- 0.1286 * 622 = 80 px on the reference set
    _scs = scale_ladder(0.1286 * Rs * 2.0, max(resolution_floor(Lsf, r_s, Rs), 2.5))
    common = dict(scales=_scs, k=0.8, global_wt=0.0,
                  valid=valid_s, norm_span=(-shalf, shalf))
    # Two full-resolution multiscale passes follow, and together they are the
    # longest thing in the run. Announce each one: without these the progress
    # bar and the log both stand still for the majority of the detail stage.
    progress.log(f"  inner corona: raw pass, scales "
                 f"{', '.join('%.1f' % x for x in _scs)} px", _DF["inner_raw"])
    raw = _deband(mgn(Lsf, floor_map=None, gains=(1, 1, 1, 1, 1, 1), **common),
              r_s, valid_s, cys, cxs, Rs + margin)
    np.save(os.path.join(wd, "inner0.npy"), _soft_norm(raw, ann).astype(np.float32))
    if do_dn:
        progress.log("  inner corona: denoising the short stack", _DF["inner_dn"])
    Ls_dn = denoise_loglum(Lsf, nfs, ks=tuple(k * 0.85 for k in ks)) if do_dn else Lsf
    progress.log("  inner corona: denoised pass", _DF["inner_cln"])
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
            # ...and the prominences' own fine structure, which no corona layer
            # can carry. See prominence_detail.
            try:
                _pd = prominence_detail(prgb[:, :, 0], g)
                _pdf = _fit(_pd.repeat(2, 0).repeat(2, 1), (H, W))
                del _pd
                _pdf = ndimage.gaussian_filter(_pdf, 1.0)
                np.save(os.path.join(wd, "promdet.npy"), _pdf.astype(np.float32))
                lstats["prom"]["detail_layer"] = True
                progress.log("prominence detail: built from the H-alpha tier's "
                             "red channel, scaled inside the gate", None)
                del _pdf
            except Exception as _e:
                progress.log(f"prominence detail layer not built ({_e})", None)
        else:
            progress.log("prominence colour: not enough limb samples", None)
    if "prom" in lstats:
        lstats["prom"]["area_px"] = int((gate > 0.3).sum())
        # ...and how much of it survives the disc mask. A prominence flagged
        # UNDER the mask contributes nothing to the picture, so a report that
        # counts it is telling the user they have prominences they cannot see.
        # Measured on an imported HDR (Val Italo's Siril stack of 327 frames):
        # 60 px flagged, every one of them between 0.97 and 1.01 R -- at or
        # inside the limb -- and 0 px visible. Beyond 1.02 R the highest
        # R/(G+B)/2 anywhere out to 1.3 R is 1.33 against a 1.35 threshold, so
        # there is no H-alpha excess left in that file to find.
        #
        # That is the gate reporting the file correctly, not misfiring. The
        # file's own embedded header says why it could not be otherwise: two
        # GHS stretches (amount 145.65 and 7.57) and an SCNR green subtraction
        # at full strength, all applied before export. What reaches the import
        # is a picture, and a picture no longer carries the chromosphere's
        # colour separately from the corona's. Lowering the threshold would
        # gate on noise; the report should say what happened instead.
        lstats["prom"]["area_visible_px"] = int(((gate > 0.3) & ~disc_m).sum())
    np.save(os.path.join(wd, "prom.npy"), gate.astype(np.float32))
    del gate, Ls

    progress.log("Pellett layer...", _DF["pellett"])
    _pellett(wd, lum, r, cy, cx, R, disc_m)

    ep = os.path.join(wd, "earth.npy")
    if earthshine:
        progress.log("earthshine layer...", _DF["earth"])
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


def prominence_detail(red_half, gate_half, floor_map=None):
    """Fine structure of the prominences, from the H-alpha tier's RED channel.

    WHY A SEPARATE LAYER, AND WHY NOT MGN.

    Prominence interiors came out flat in every corona layer, and both external
    testers said so. Two measured reasons, on the reference bracket's biggest
    prominence (225 deg, R/GB 16.4):

    1. MGN's NORMALISATION WINDOW CLIPS THEM. `hi` is the 99.95th percentile of
       the frame -- the right choice for a corona, where a few hot pixels must
       not set the range. But a prominence is brighter than that: 37% of this
       one sat hard-clipped at xn = 1.0 in the red channel and 21% in the merged
       luminance, with the unclipped remainder squeezed into the top 6% of the
       range. Whatever structure it had was gone before the multiscale filter
       ran. Giving the layer its OWN window, taken inside the gate, drops the
       clipping to 2%.

    2. MGN IS THE WRONG FILTER FOR A COMPACT BRIGHT FEATURE. Its purpose is to
       divide out the local standard deviation so faint structure at any
       brightness comes up equally -- which is exactly what flattens a
       prominence, whose interior variation IS its local sigma. A plain
       multiscale unsharp keeps it. Correlation of the resulting layer's fine
       structure with the red channel's own, which no normalisation can fake:

           MGN, frame window (what the corona layers do)        0.037
           MGN, gate window, fine scales                        0.420
           MGN, gate window, corona scales                      0.375
           plain multiscale unsharp, gate-scaled                0.940   <-- this

       For reference the existing layers score inner 0.317, MGN 0.274, merged
       HDR 0.262 on the same measure.

    THE RED CHANNEL, not luminance. An H-alpha prominence puts most of its
    signal in R, and luminance weights R at 0.2126: going to luminance costs
    more than half the structure before any filter sees it (1.000 -> 0.417),
    which is the single largest loss in the chain.

    Returns a 0..1 layer at the input's resolution, 0.5 where there is nothing.
    """
    L = np.log10(np.clip(np.asarray(red_half, np.float32), 1.0, None))
    acc = np.zeros_like(L)
    # octave spacing with decreasing weight, the usual unsharp ladder; the
    # finest scale carries most of it because prominence structure is fine
    for sig, gain in ((1.0, 1.0), (2.0, 0.7), (4.0, 0.5), (8.0, 0.35)):
        acc += gain * (L - ndimage.gaussian_filter(L, sig))
    del L
    # Scale by the spread INSIDE the prominences, so the layer uses its range on
    # them rather than on whatever else the frame contains. Dilated, so the
    # scaling is not set by the few brightest cores alone.
    sel = ndimage.binary_dilation(np.asarray(gate_half, np.float32) > 0.05,
                                  iterations=6)
    if sel.sum() < 200:
        sel = np.ones_like(acc, bool)
    v = acc[sel]
    sd = 1.4826 * float(np.median(np.abs(v - np.median(v))))
    if not np.isfinite(sd) or sd <= 0:
        sd = float(np.std(acc)) or 1.0
    out = np.clip(0.5 + acc / (6.0 * sd), 0.0, 1.0)
    del acc
    return out.astype(np.float32)


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
