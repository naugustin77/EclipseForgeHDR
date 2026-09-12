"""Detail-extraction layers: photon-noise-adaptive MGN, robust FNRGF,
inner corona (raw + denoised), prominences, earthshine."""
from __future__ import annotations
import os, json
import numpy as np
from scipy import ndimage


from .nafe import nafe_vn

# NAFE-VN defaults. w and gamma are the published working values -- Habbal,
# Druckmuller 2013 (ApJS 207:25) and Druckmuller & Druckmullerova 2014 (LNCS
# 8466 p.262). Published working values, now checked against both primaries:
# w in <0.05, 0.3> with 0.2 used throughout their figures, and sigma in
# <2 sigma_A, 12 sigma_A> -- rather than values guessed here. sigma is in units of
# the noise sigma_A and is set adaptively per level. The rest measured on the
# reference set:
#  * eps 0.05 (in rank units) recovers most of the near-limb contrast the
#    unrestricted filter loses -- near/far detail ratio 0.22 -> 0.34 -- while
#    a wider window falls back towards no restriction at all.
#  * K 64 with a Gaussian membership 1.5 bins wide keeps the cumulative
#    histogram smooth; a narrower membership prints concentric contour rings.
NAFE_K = 128          # histogram levels. NOT A PAPER VALUE -- both papers bin
                      # NOTHING: their histogram runs over the actual discrete
                      # pixel values through a Kronecker delta (2014 eq. 6), on
                      # 14-bit AIA data, i.e. 16384 levels. Our own sweep is the
                      # only justification: 64 measurably cost corona contrast.
                      # (The "several thousand discrete pixel values" line this
                      # once leaned on is real -- 2013 p.2, 2014 p.265 -- but it
                      # argues the CONTINUOUS approximation is safe, which is an
                      # argument for many levels, not for 128.)
NAFE_W = 0.2          # the paper's eq.2 weight; the stored layer is E, so this
                      # is unused here -- nafeMix in render.py is the live w
NAFE_GAMMA = 2.4      # DEAD: only used on nafe_vn's combine=True branch, and the
                      # pipeline always calls combine=False. Left for the day
                      # that branch is used; changing it does nothing today.
NAFE_EPS = 0.10       # value window, in LEVEL units (nafe.py builds a linear
                      # level axis, not a rank axis -- the old "rank units" here
                      # was wrong and contradicted nafe.py's own docstring).
                      # The restriction is CRISP 0/1 (2014 eq. 12), which is what
                      # nafe.py implements. NO NUMERIC VALUE IS PUBLISHED: "The
                      # value of eps must be found experimentally" (2014 p.268).
                      # The paper does support the two failure DIRECTIONS our
                      # sweep found -- too small "cause image fragmentation into
                      # small areas with very high contrast features", too large
                      # "will result in ... nearly identical" to the unrestricted
                      # filter -- but the number 0.10 is ours. Swept on the reference set:
                      # 0.02 prints concentric rings (the paper's "fragmentation"),
                      # 0.05 gave high-pass 0.076, 0.10 gives 0.109 at the same
                      # agreement with FNRGF (0.711 -> 0.721), and past 0.15 the
                      # curve is flat and the extra is noise.
NAFE_NOISE_MULT = 4.0 # the paper's sigma in units of sigma_A; their range is 2..12
# THE NEIGHBOURHOOD SIGMA, AND A BUG THAT TURNED OUT TO BE THE GOOD SETTING.
#
# This was 0.13 and the call site divided it by NAFE_GRID before passing it --
# and nafe_vn() then divides by `grid` AGAIN to reach its decimated working
# grid (nafe.py, `ssp = max(sigma_sp / q, 0.7)`). So the neighbourhood actually
# applied was 0.13 R / 8 = 0.016 R: on the reference frame, 10 px where the
# constant, the docstring and the run report all said 80 px.
#
# Correcting it to a true 0.13 R measures much BETTER and looks worse. Scored
# on the reference bracket as amp*coh:
#
#     shell          as shipped (0.016 R)   "fixed" (0.13 R)
#     1.05-1.30 R    0.194 (coh 0.97)       0.939 (coh 0.99)
#     1.30-1.80 R    0.126 (coh 0.86)       0.430 (coh 0.96)
#     1.80-2.60 R    0.079 (coh 0.53)       0.135 (coh 0.75)
#
# Nearly five times the score at the limb -- and side by side the 0.13 R layer
# is blobby and lobed where the 0.016 R one is filamentary. Same lesson as the
# and the same lesson as the scale ladder: amp*coh rewards large coherent
# structure and cannot see "delicate". The picture decides.
#
# So the double division is removed and the CONSTANT is restated as what was
# always actually running. Behaviour is unchanged; the number now means what it
# says, and the run report stops lying about it. The old sweep that reported
# this "flat from 0.06 R to 0.51 R" ran through the bug, so it really swept
# 0.0075 R to 0.064 R -- 0.13 R was never tested until now.
NAFE_NEIGH_R = 0.016  # neighbourhood sigma as a fraction of the lunar radius.
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
# FNRGF: three published behaviours that were implemented and unreachable until
# 0.23.2. Each was then MEASURED on the 600 mm reference bracket, and the measurements
# disagree with the literature on one of them. Defaults follow the measurement.
#
# THE RING ARTIFACT the first two address is described by the filter's own
# author -- thesis p.81, rings "like on a gramophone disk", caused by impulse
# noise or by "a very high-contrast structure such as a prominence" corrupting
# one ring's polynomial and not its neighbour's.
#
# HOW IT WAS MEASURED: the median of the layer around each inner ring over
# corona pixels, prominence pixels and a 25 px collar excluded from the
# MEASUREMENT in every variant so all are scored on identical pixels. FNRGF
# normalises each ring to zero mean, so a ring's median should sit near zero; a
# prominence pulls the fitted mean up over the whole circle and pushes the
# corona's median down. The gap between rings that carry a prominence and rings
# that do not is therefore the artifact, in the layer's own units:
#
#     shipped 0.23.1      0.0009
#     prominence mask     0.0005     <- halved
#     impulse median      0.0008     <- nothing
#     both                0.0006
#
FNRGF_IMPULSE = 0     # OFF, against the literature, on measurement.
                      # An azimuthal median over the ring samples (fit only).
                      # Thesis p.86 medians the whole input before FNRGF in both
                      # space-based examples. Here it buys nothing -- 0.0009 to
                      # 0.0008 on the gap above -- because the master is already
                      # hot-pixel-repaired (1848 photosites on the reference
                      # bracket) and denoised by the time this runs, so there is
                      # no impulse left to remove. And it is NOT free: it moved
                      # the layer by 29% rms. At na = 1440 on a 3900 px
                      # circumference the samples are 2.7 px apart, so a
                      # 3-sample median spans 8 px of arc on an image that
                      # resolves to 1.4 px -- it is smoothing real azimuthal
                      # structure out of the background model. A large change
                      # with no measurable benefit is a cost, so it ships off.
                      # Kept as a parameter for data that IS impulse-dirty,
                      # which is the case the thesis is describing.
FNRGF_NOISE_FRAC = 0.10   # ON. sqrt(Vn) as a fraction of the measured
                      # outer-field residual (thesis eqs. 6.11-6.13; her band is
                      # 10-15%, matching NAFE's independent tests on AIA data).
                      # Held at her LOWER bound because we denoise before this
                      # runs, so the noise it sees is less than the thesis
                      # assumes. Structure spread against outer-field spread,
                      # higher is better:
                      #
                      #     0.00  0.841     0.10  0.847     0.30  0.893
                      #     0.05  0.843     0.15  0.855
                      #
                      # Monotone and small. Taken because the direction is right
                      # and the floor it supplements is one this file already
                      # describes as the wrong shape. 0 disables it.
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
    "mgn_fine":   99.0,       # same filter, three scales instead of six
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
    halo band around the limb).

    WHERE THIS DEPARTS FROM THE PAPER, stated because the report cites it.
    Eq. 5 there is

        I = h C'_g + (1 - h) (1/n) SUM g_i C'_i

    where C'_g is a globally gamma-normalised version of the whole image, added
    "to give contextual information of the largest scale structure", and the
    paper uses h = 0.7. So in Morgan & Druckmuller most of the output is that
    global term.

    Here `global_wt` is 0.12 by default and the inner-corona passes call with
    0.0. That is deliberate and it follows from a difference in what the layer
    is FOR. Their MGN is the finished picture of an AIA disc, so it has to
    carry its own large-scale context or the disc stops looking like a disc.
    Ours is a DETAIL layer: the renderer builds the brightness envelope
    separately and adds this on top of it, and the radial profile has already
    been removed before normalisation. A large global term would put the
    envelope back and the two would double up.

    `k` is 0.8 at the corona call sites against the paper's 0.7, which is a
    taste setting on the arctan and nothing more.

    Everything else follows the paper: the scale ladder w = 1.25, 2.5, 5, 10,
    20, 40, the per-scale gains from its Fig. 4, and gamma = 3.2."""
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
    # DIVIDE BY THE NUMBER OF SCALES USED, NOT THE SUM OF ALL SUPPLIED GAINS.
    #
    # The paper's Eq. 5 is  I = h*C'_g + ((1-h)/n) * sum_i g_i C'_i,  with n the
    # number of scales. This divided by sum(gains) instead. With the paper's six
    # scales and six gains that is a harmless 2% (5.874 vs 6) -- but `scales` is
    # built per image by scale_ladder(), which COLLAPSES scales the resolution
    # floor has merged and routinely returns fewer than six, while `gains` stays
    # six long. zip() then uses the first n gains and this line still divided by
    # all six, so the whole layer came out under-scaled:
    #
    #   scales used   5      4      3
    #   amplitude   x0.83  x0.66  x0.49
    #
    # A 240 mm lens or any camera whose floor merges the fine end -- exactly the
    # case resolution_floor() was written for -- lost up to half the layer, and
    # `mgnContrast` silently absorbed it.
    _g = list(gains)[:len(list(scales))]
    acc /= max(len(_g), 1)
    out = (0.5 + acc / np.pi)
    if global_wt > 0:
        out = global_wt * xn ** (1 / global_gamma) + (1 - global_wt) * out
    if m is not None:
        out = out * m + 0.5 * (1 - m)
    return out


# ---------------------------------------------------------------------------
# HILL'S ENHANCEMENT CHAIN.
#
#     im_enhanced = im_log + a*M_1 + b*M_2 + c*M_4 + ...
#     M_s (unsharp mask of size s) = im_log - im_blur(s)
#
# straight off his slide, with im_blur the POLAR-ORIENTED, PARTIAL blur of the
# log-mapped image and the amplification factors his example's 100 / 60 / 20 /
# 10 at 2 / 4 / 8 / 16 px.
#
# WHY THIS IS NOT WHAT THE OTHER LAYERS DO, which is the whole point. MGN
# and NAFE all NORMALISE: they divide the local residual by the local standard
# deviation, or replace it by its rank. That is what makes them able to show the
# inner and outer corona at once, and it is also what makes every part of the
# frame end up with the SAME texture amplitude -- faint structure and noise come
# out equally strong, which reads as coarse. Hill's combination is LINEAR: a
# fixed scalar per scale, nothing divided by anything local. Faint stays faint.
# That difference, not the filter kernels, is why his slides look finer than
# our composite did.
#
# WHY THE POLAR BLUR MATTERS, and why the note below was the wrong test for it.
# In the polar image the blur is a plain Gaussian of sigma s in BOTH axes, so
# back in the picture the radial width is s px everywhere while the tangential
# width is s * r / r_max -- narrow near the limb, growing outwards. The kernel
# therefore averages ALONG a streamer and not across it, by a margin that
# increases with radius. Hill: "you'll have poor noise performance in the outer
# corona in particular ... the difference between an image that's sharpened like
# this and just a highpass filter in Photoshop is that polar transformation."
#
# The ACHF note below measured a DIFFERENT kernel -- Druckmullerova's, whose
# arc-length term makes the tangential width a constant physical distance, which
# is what makes it nearly isotropic at our scales. That measurement is correct
# and it does not apply here: Hill's kernel is anisotropic BY CONSTRUCTION and
# by a factor of r_max/r, which on this geometry is about 8x at the limb.
# ---------------------------------------------------------------------------

HILL_SCALES = (2.0, 4.0, 8.0, 16.0, 32.0)
# Hill's own example is 100 / 60 / 20 / 10 at 2 / 4 / 8 / 16 px. Scaled to 1.0
# at the finest and extended by one octave for frames larger than his.
HILL_GAINS = (1.0, 0.6, 0.2, 0.1, 0.05)


def polar_partial_blur(f, w, cy, cx, sigmas, band=384, progress=None):
    """Hill's im_blur: blur f in heliocentric polar coordinates, with the
    convolution PARTIAL so the masked pixels (w = 0) contribute nothing.

    f    log-mapped image
    w    weight image -- 0 on the Moon, prominences and anything else excluded
    sigmas  one blur per scale, in pixels of RADIUS

    Returns one blurred image per sigma, in the same order.

    Only the blurred component makes the polar round trip. The unsharp mask is
    taken in Cartesian space as f - blur, so the interpolation of the warp can
    soften the blur -- which is already smooth -- but can never touch the detail
    the mask is there to carry.

    Worked in radial bands with a 4-sigma overlap, so the peak footprint is one
    band of the polar grid rather than the whole of it. The polar grid is
    sampled at one pixel of RADIUS and one pixel of ARC AT r_max, which is the
    natural choice: it is Nyquist-matched to the picture at the outer edge and
    oversampled everywhere inside it.
    """
    H, W = f.shape
    yy = np.arange(H, dtype=np.float32)[:, None] - cy
    xx = np.arange(W, dtype=np.float32)[None, :] - cx
    rc = np.sqrt(yy * yy + xx * xx)
    tc = np.arctan2(yy, xx)
    rmax = float(rc.max()) + 1.0
    nr = int(np.ceil(rmax))
    nth = int(np.ceil(2.0 * np.pi * rmax))
    dth = 2.0 * np.pi / nth
    tcol = (np.mod(tc, 2.0 * np.pi) / dth).astype(np.float32)
    del tc
    fw = (f * w).astype(np.float32)
    outs = [np.empty((H, W), np.float32) for _ in sigmas]
    # COVERAGE: how much real data each output pixel's kernel actually saw.
    # Needed because a pixel the mask excluded, or one ringed by excluded
    # pixels, gets a blur estimated from almost nothing -- see build_hill.
    covs = [np.zeros((H, W), np.float32) for _ in sigmas]
    pad = int(np.ceil(4.0 * max(sigmas))) + 2
    th = (np.arange(nth, dtype=np.float32) * dth)
    cth, sth = np.cos(th), np.sin(th)
    del th
    nb = (nr + band - 1) // band
    for bi, r0 in enumerate(range(0, nr, band)):
        r1 = min(nr, r0 + band)
        a0, a1 = max(0, r0 - pad), min(nr, r1 + pad)
        ra = np.arange(a0, a1, dtype=np.float32)
        sy = (cy + ra[:, None] * sth[None, :])
        sx = (cx + ra[:, None] * cth[None, :])
        crd = np.stack([sy, sx])
        del sy, sx
        FW = ndimage.map_coordinates(fw, crd, order=1, mode="nearest")
        WW = ndimage.map_coordinates(w, crd, order=1, mode="nearest")
        del crd
        sel = (rc >= r0) & (rc < r1)
        if not sel.any():
            del FW, WW
            continue
        rr = (rc[sel] - a0).astype(np.float32)
        tt = tcol[sel]
        dst = np.stack([rr, tt])
        for k, sg in enumerate(sigmas):
            # 'wrap' in angle: the polar image is periodic, and treating its two
            # edges as a boundary puts a seam along one radius of the picture
            num = ndimage.gaussian_filter(FW, sg, mode=["nearest", "wrap"])
            den = ndimage.gaussian_filter(WW, sg, mode=["nearest", "wrap"])
            B = num / np.maximum(den, 1e-4)
            del num
            outs[k][sel] = ndimage.map_coordinates(B, dst, order=1,
                                                   mode="nearest")
            covs[k][sel] = ndimage.map_coordinates(den, dst, order=1,
                                                   mode="nearest")
            del B, den
        del FW, WW, dst, rr, tt, sel
        if progress is not None:
            progress.log(f"  partial convolution: radial band {bi + 1}/{nb}", None)
    return outs, covs


# ---------------------------------------------------------------------------
# ACHF: TESTED AND REJECTED as a DROP-IN FOR OUR SCALES. This is Druckmullerova's
# kernel, not Hill's -- see the block above before reading it as a rejection of
# polar filtering in general.
#
# Druckmullerova's thesis calls ACHF "the best nowadays used structure
# enhancement technique for images of the solar corona from total solar eclipses
# in white light" (p.68), and it is the one method in that literature we do not
# have. Its kernel (thesis eq. 5.1, p.67) is
#
#     C_{r,phi}(rho, varphi) = exp( -[ (r-rho)^2 + (r(varphi-phi))^2 ] / 2 sigma^2 )
#
# -- a Gaussian in radius crossed with a Gaussian in ARC LENGTH along a
# Sun-centred circle -- applied as f minus its normalized convolution (eq. 5.2),
# with the Moon excluded by the weight function.
#
# Implemented exactly, including the arc-length term (which is what this
# session's earlier polar-convolution attempt got wrong: it used a constant
# kernel in polar-GRID units, i.e. an arc length growing linearly with radius),
# and with only the blurred component making the polar round trip so that
# interpolation cannot contaminate the residual. Then compared against an
# isotropic kernel through an IDENTICAL normalise-and-arctan path with an
# IDENTICAL scale ladder, so the only difference was the kernel shape:
#
#     shell          isotropic        ACHF           change
#     1.02-1.15 R    0.0701           0.0770          +10%
#     1.15-1.40 R    0.0610           0.0607           -1%
#     1.40-1.80 R    0.0438           0.0450           +3%
#     1.80-2.60 R    0.0465           0.0491           +6%
#     2.60-3.40 R    0.0530           0.0554           +4%
#
# and side by side the two images are indistinguishable. Cost: 390s against 21s
# at half resolution, so roughly 26 minutes a run at full resolution.
#
# THE REASON, which is the part worth keeping. The kernel is anisotropic only
# where sigma is a large fraction of r. At our scales it never is:
#
#     sigma 1.25 px at 1.05 R  ->  0.110 deg      sigma 40 px at 1.05 R -> 3.5 deg
#     sigma 1.25 px at 3.40 R  ->  0.034 deg      sigma 40 px at 3.40 R -> 1.1 deg
#
# An isotropic Gaussian of the same sigma covers the same arc to within
# (sigma/r)^2 -- 0.4% at the worst point in our ladder, 4e-6 at the finest. The
# two kernels are numerically almost the same object everywhere we use them.
# ACHF's shape only separates from isotropic at sigma of several hundred pixels,
# which is large-scale structure, not the fine detail this was reached for.
#
# What this does NOT reject: ACHF as Druckmuller ships it, which has 38 tuned
# parameters and a "nonlinearity in the use of filtered images" (thesis p.68)
# that no source specifies. It rejects the kernel SHAPE as a drop-in for ours,
# at our scales, which is what the reading of the literature suggested trying.
# ---------------------------------------------------------------------------


def prominence_mask(wd, geo, shape, r, R, pct=99.0, grow=8, rmax=1.35,
                    progress=None):
    """Pixels to EXCLUDE from the partial convolution, beyond the disc itself.

    Partial (normalized, incomplete) convolution is already how every masked
    filter here estimates a local mean -- (w.f * C) / (w * C), with the mask
    convolved by the same kernel. What was missing is what goes IN the mask.
    Hill's rule is "Moon, proms, stars should be 0"; ours was the Moon alone.

    Why prominences matter more than their size suggests: S, the local standard
    deviation, is MGN's DIVISOR. A prominence is the brightest thing in the
    frame and sits hard against the limb, so every coronal pixel within a kernel
    of one has its contrast divided by a sigma that the prominence set. Masking
    it out gives that annulus its own contrast back.

    MEASURED on the reference bracket, with the limb split by azimuth into rays
    that contain a prominence and rays that do not. The second column is the
    control and it is in the same run: only the near column moves, so this is
    the prominences and not a generic consequence of masking more.

        pct grow  mask %ring   1.02-1.15 near/away   1.15-1.40 near/away
         99   2      1.60%       +6.8% /  +0.2%        -0.0% / +0.2%
         99   5      2.63%      +10.8% /  +0.2%        -0.2% / +0.3%
         99   8      3.73%      +20.2% /  +0.3%        -0.3% / +0.5%   <- shipped
         98   5      5.23%      +24.9% /  +0.3%        +0.6% / +0.4%
         97   5      8.00%      +34.1% /  +0.4%        +4.0% / +0.3%
         99  12      5.27%      +31.8% /  +0.7%        -0.5% / +0.9%

    Two things that table says and a single number would have hidden. The effect
    is confined to 1.02-1.15 R -- past that it is nothing until the mask gets
    large. And it does not saturate: mask more, gain more, with no natural stop.
    So the size is a TRADE, not an optimum -- corona contrast against pixels
    flattened to 0.5 -- and 99/8 is a conservative point on it, not a maximum.
    Anyone is free to move it; the cost of moving it is in the third column.

    STARS ARE NOT MASKED. Hill lists them, and the same argument applies, but
    the 2026 brackets have none to mask yet -- so there is nothing here to
    measure against and nothing is claimed. A star mask belongs with the first
    frame that actually shows stars.

    Returns None when the H-alpha tier was not stacked, in which case callers
    keep the disc-only mask and behave exactly as before.
    """
    p = os.path.join(wd, "prom_rgb.npy")
    if not os.path.exists(p):
        return None
    try:
        prgb = np.nan_to_num(np.load(p).astype(np.float32), nan=0.0,
                             posinf=0.0, neginf=0.0)
        red = ndimage.gaussian_filter(prgb[:, :, 0], 2) / np.maximum(
            ndimage.gaussian_filter(0.5 * (prgb[:, :, 1] + prgb[:, :, 2]), 2), 1e-3)
        del prgb
        # prom_rgb is a half-res single tier; upsample onto the layer grid
        m = _fit(red.repeat(2, 0).repeat(2, 1), shape)
        # Never mask beyond the limb ring: a threshold on redness alone would
        # start eating red-biased sky far out, where there is no prominence to
        # exclude and the corona is all we have.
        ring = (r > R) & (r < rmax * R)
        if ring.sum() < 1000:
            return None
        # AN ABSOLUTE FLOOR AS WELL AS THE PERCENTILE.
        #
        # A bare percentile ALWAYS FIRES. On a bracket with no prominences to
        # speak of it flags the reddest 1% of the ring -- which is noise and
        # ordinary corona -- and the 8 px dilation turns that into a scatter of
        # patches. the 250 mm test set: the real H-alpha gate, which does have
        # an absolute threshold, finds 349 px outside the disc; this mask was
        # punching out 11000, about 32x more. Those holes printed as a row of
        # smooth blobs in an arc below the Moon once Hill's step 4 zeroed the
        # masks inside them (0.22.71) -- and as high-contrast smudges before it.
        #
        # So require the redness to be a genuine outlier against the corona's
        # own, the same way the prominence gate does, AND keep the percentile as
        # a cap so a frame full of prominences cannot mask half the ring. When
        # there is nothing red enough, the floor wins and nothing is masked --
        # which is the correct answer, not a fallback.
        _mr = m[ring]
        _med = float(np.median(_mr))
        _sig = 1.4826 * float(np.median(np.abs(_mr - _med))) + 1e-9
        _floor = _med + 6.0 * _sig
        _pctv = float(np.percentile(_mr, pct))
        out = ring & (m > max(_pctv, _floor))
        if progress is not None:
            progress.log(
                f"  prominence mask: redness median {_med:.2f}, "
                f"floor {_floor:.2f} (median + 6 sigma), "
                f"{pct:g}th percentile {_pctv:.2f} -> "
                f"{'the floor' if _floor > _pctv else 'the percentile'} decides, "
                f"{int(out.sum()) / 1e3:.1f}k px before dilation", None)
        if grow > 0:
            out = ndimage.binary_dilation(out, iterations=int(grow)) & ring
        return out if out.any() else None
    except Exception:
        return None


def fnrgf_robust(lum, r, cy, cx, r0, order=6, na=1440, rmax=None,
                 atte_ave=0.0, atte_dev=0.0, var_floor=0.3, sd_smooth=1.0,
                 noise_var=None, noise_frac=0.0, mask=None, impulse=3):
    """Fourier Normalizing Radial Gradient Filter (Druckmullerova et al. 2011).

    The Fourier order is reduced as a ring leaves the frame: fitting 13 free
    coefficients to a short covered arc oscillates wildly and jumps from ring to
    ring, which is what produced the concentric arc artifacts in the outer
    field. Coverage-matched order plus ridge damping of the higher harmonics
    keeps consecutive rings consistent.

    THE RING ARTIFACT, AND THE AUTHOR'S OWN ACCOUNT OF IT. Druckmullerova's
    doctoral thesis p.81 says this filter produces "visible rings (like on a
    gramophone disk) ... caused by the fact that there is no connection between
    the processing of one Sun-centered ring and the neighboring one. Each circle
    is processed separately. If there is significant impulse noise in the image
    (especially more faulty pixels next to each other) or if there is a very
    high-contrast structure such as a prominence, it affects the trigonometric
    polynomial for the whole circle."

    Two named causes, and `mask` and `impulse` are one each.

    `mask` (True = exclude) is for the prominence. IRLS with a Huber weight only
    DOWN-weights it: the weight is 1/max(u,1), so a 10 sigma outlier still enters
    at 0.2, and a prominence is not one outlier but a contiguous arc of them on
    exactly the inner rings where the artifact shows. We already build
    `prominence_mask` for the partial convolution; this hands it the same one.

    `impulse` is an azimuthal median over the ring SAMPLES, applied to the fit
    only. Her own remedy is a one-pass median on the whole input image, which she
    applies before FNRGF in both space-based examples in the thesis (p.86). Doing
    that here would cost real resolution -- this pipeline resolves to about 1.4 px
    and the input is already hot-pixel-repaired and denoised. Filtering the ring
    samples instead removes an impulse from the BACKGROUND MODEL, which is the
    thing the thesis says it corrupts, and leaves the numerator untouched, so no
    resolution is lost anywhere. 0 disables it."""
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
    # PER-HARMONIC ATTENUATION, applied when the polynomial is EVALUATED, with
    # SEPARATE series for the mean and the standard deviation. That is how
    # Druckmuller's own program does it (FNRGFsoftware/ImgProc.pas,
    # PointAttenFourier) and the thesis says why a hard cutoff will not do:
    # "The necessity to treat information on different spatial frequencies in
    # a different manner was inevitable." It also gives the constraint -- the
    # deviation series must not out-run the mean series -- and the failure mode
    # of pushing the mean series too far: "artificial brightenings in
    # low-contrast parts of the image ... false glimmers of the higher-order
    # sine and cosine functions."
    #
    # `atte_*` is the linear decrement per harmonic; 0 keeps a hard cutoff at
    # `order` with the ridge alone, which is what shipped up to 0.22.28.
    _ka = np.arange(order + 1)
    def _att(step):
        v = np.clip(1.0 - step * _ka, 0.0, 1.0) if step else np.ones(order + 1)
        out = np.ones(2 * order + 1)
        out[0] = v[0]
        for m in range(1, order + 1):
            out[2 * m - 1] = out[2 * m] = v[m]
        return out
    att_a, att_d = _att(atte_ave), _att(atte_dev)
    L = np.log10(np.clip(lum, 1.0, None))
    mu_g = np.zeros((nr, na), np.float32)
    sd_g = np.zeros((nr, na), np.float32)
    cov_r = np.zeros(nr, np.float32)
    n_masked = np.zeros(nr, np.int32)
    v_g = np.zeros((nr, na), np.float32)
    sg_r = np.full(nr, 1e-6, np.float32)
    dead = np.zeros(nr, bool)
    for i in range(nr):
        rad = r0 + i
        ys = cy + rad * sa
        xs = cx + rad * ca
        ok = (ys >= 0) & (ys <= H - 1) & (xs >= 0) & (xs <= W - 1)
        if mask is not None:
            # Sample the mask on the same ring. `> 0.5` after order-0 lookup, so
            # a sample is dropped only when it lands ON a masked pixel.
            _mk = ndimage.map_coordinates(
                mask.astype(np.uint8),
                [np.clip(ys, 0, H - 1), np.clip(xs, 0, W - 1)],
                order=0, mode="nearest")
            _keep = _mk < 1
            # Never let the mask starve a ring: if it would take more than half
            # the covered arc, the fit is worse off without those samples than
            # with them, and the Huber weight is the fallback it always was.
            if int((ok & _keep).sum()) >= max(24, int(0.5 * ok.sum())):
                ok = ok & _keep
                n_masked[i] = int((~_keep).sum())
            del _mk, _keep
        nok = int(ok.sum())
        cov = nok / float(na)
        cov_r[i] = cov
        if nok < 24:
            dead[i] = True
            mu_g[i] = mu_g[i - 1]
            v_g[i] = v_g[i - 1]
            sg_r[i] = sg_r[i - 1]
            continue
        # order the covered arc can actually support
        o = int(np.clip(np.floor(order * cov * 1.2), 0, order))
        while o > 0 and nok < 12 * (2 * o + 1):
            o -= 1
        nc = 2 * o + 1
        A = Adm[ok][:, :nc]
        Afull = Adm[:, :nc]
        rg = np.diag(ridge[:nc])
        v = ndimage.map_coordinates(L, [ys[ok], xs[ok]], order=1)
        if impulse and nok >= 3 * int(impulse):
            # Circular median along the ring, on the FIT SAMPLES only. `wrap` is
            # correct just when the whole ring is covered; on a partial arc the
            # two ends are not neighbours.
            v = ndimage.median_filter(
                v, size=int(impulse),
                mode=("wrap" if nok == na else "nearest"))
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
            res = v - A @ (coef * att_a[:nc])
            sg = max(1.4826 * np.median(np.abs(res)), 1e-6)
            u = np.abs(res) / (2.0 * sg)
            w = 1.0 / np.maximum(u, 1.0)
        mu_g[i] = Afull @ (coef * att_a[:nc])
        r2 = np.clip(res, -4 * sg, 4 * sg) ** 2
        Aw = A * w[:, None]
        coef2 = np.linalg.solve(Aw.T @ A + rg * nok, Aw.T @ r2)
        # THE FLOOR ON THE NORMALISING VARIANCE.
        # Ours is RELATIVE to the local robust scale `sg`, so in the inner
        # corona -- where `sg` is large because the residual there is full of
        # real structure -- the floor rises with the signal and divides that
        # structure back out. Lowering it is not the answer either: a relative
        # floor can collapse toward zero where the fit happens to be very good,
        # and then the division blows up (measured: flat grey with blown
        # blobs). Druckmuller's program does neither. It ADDS an absolute noise
        # variance, estimated once from the outermost rings
        # (EstimateAdditiveNoiseRing), which cannot rise with the signal and
        # cannot collapse to zero.
        #
        # THE MAGNITUDE, because the name is misleading. Vn is NOT the image's
        # noise variance -- it is the variance of a small amount of ARTIFICIAL
        # noise added on top (thesis eqs. 6.11-6.13, sqrt(S^2 + Vn)). Thesis
        # p.83: sqrt(Vn) is "about 10 to 15 percent of the standard deviation of
        # the noise originally contained in the image", which is 1-2% of that
        # noise's VARIANCE, not 100% of it. Estimated as the median of the
        # segment standard deviations over the 10 outermost processed rings. She
        # is also explicit that the opposite -- SUBTRACTING the noise variance,
        # as her own 2011 paper published -- is "a misleading idea" that gives
        # "infinite amplification" wherever the image is noise alone.
        v_g[i] = Afull @ (coef2 * att_d[:nc])
        sg_r[i] = sg
    # THE NORMALISING VARIANCE, chosen here rather than in the loop because the
    # published estimator needs the whole radial range before it can be formed.
    if noise_frac and noise_frac > 0:
        # Her estimator: the median of the segment standard deviations over the
        # 10 OUTERMOST processed rings, with sqrt(Vn) a fixed fraction of it.
        # `sg_r` is the robust scale of the residual after the fitted mean is
        # removed, in the same log10 units as the variance -- in the outer field
        # that residual IS the noise, which is what she is after.
        _use = np.flatnonzero(cov_r >= 0.35)
        _tail = _use[-10:] if _use.size else np.arange(max(nr - 10, 0), nr)
        _sig_n = float(np.median(sg_r[_tail])) if _tail.size else 0.0
        _Vn = (float(noise_frac) * _sig_n) ** 2
        # HER TERM IS AN ADDITION TO OUR FLOOR, NOT A REPLACEMENT FOR IT, and
        # this was measured the wrong way round first. Her S^2 is the variance
        # of PIXEL VALUES in a segment, which is never near zero. Ours is the
        # fitted variance of the RESIDUAL after a robust order-6 fit, which goes
        # to zero wherever the fit is good -- and sqrt(Vn) at 10% of the outer
        # noise is 1% of its variance, nowhere near large enough to floor that.
        # Dropping the floor for her term alone inflated the inner corona's
        # spread ninefold on the reference bracket: the division blowing up,
        # which is exactly what the floor exists to prevent. Keep both.
        var = np.maximum(v_g, (var_floor * sg_r[:, None]) ** 2) + _Vn
        fnrgf_robust.last_noise = (_sig_n, _Vn, len(_tail))
    elif noise_var is not None:
        var = np.maximum(v_g, 0.0) + float(noise_var)
        fnrgf_robust.last_noise = (float("nan"), float(noise_var), 0)
    else:
        var = np.maximum(v_g, (var_floor * sg_r[:, None]) ** 2)
        fnrgf_robust.last_noise = None
    sd_g = np.sqrt(np.maximum(var, 0.0)).astype(np.float32)
    fnrgf_robust.last_masked = int((n_masked > 0).sum())
    fnrgf_robust.last_dead = int(dead.sum())
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
    sd_g = (1 - blend) * ndimage.gaussian_filter(sd_g, (10 * sd_smooth, 0)) + \
        blend * ndimage.gaussian_filter(sd_g, (40 * sd_smooth, 0))
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


HILL_BUILD = 6      # bump to force a rebuild of cached masks; see _deradial


def _radial_median(a, r, valid, nb, smooth=5):
    """Azimuthal median of `a` at each integer radius, lightly smoothed."""
    ri = np.clip(r.astype(np.int32), 0, nb - 1)
    v = ri[valid]
    av = np.asarray(a, np.float32)[valid]
    o = np.argsort(v, kind="stable")
    v, av = v[o], av[o]
    e = np.searchsorted(v, np.arange(nb + 1))
    prof = np.zeros(nb, np.float32)
    last = 0.0
    for i in range(nb):
        s = av[e[i]:e[i + 1]]
        if s.size:
            last = float(np.median(s))
        prof[i] = last
    return ndimage.uniform_filter1d(prof, smooth, mode="nearest")


def _mask_bias(m, r, valid, R, lo=1.05, hi=1.30):
    """The mask's own MEAN over a near-limb shell, in units of its own rms.

    This is the number that measures the featureless collar: a mask is supposed
    to carry azimuthal structure about zero, so a large mean is an offset the
    picture will render as a smooth bright or dark ring with nothing in it.
    """
    sel = valid & (r > lo * R) & (r < hi * R)
    if sel.sum() < 100:
        return 0.0
    v = np.asarray(m, np.float32)[sel]
    sd = float(np.std(v))
    return float(np.mean(v) / sd) if sd > 1e-12 else 0.0


def _deradial(m, r, valid):
    """Remove the mask's azimuthal median at each radius.

    THE FEATURELESS COLLAR AROUND THE MOON, AND WHY IT IS THERE.

    Partial convolution stops the black disc leaking into the blur, which is
    what it is for and it does it. What it cannot fix is that at a radius just
    outside the limb the kernel is ONE-SIDED: every valid sample lies further
    out, where the corona is fainter, so the local mean sits below the centre
    pixel and the residual comes out systematically POSITIVE. On a profile as
    steep as an eclipse corona at 1.05 R that offset is large, smooth, and
    purely a function of radius -- a bright ring with no structure in it.

    Measured on a synthetic with the reference bracket's geometry, as the
    mask's mean over a shell divided by its own rms (so 1.0 means the offset is
    as large as all the real structure put together):

        shell        32 px mask   after this
        1.03-1.08 R     +6.97       +0.55
        1.08-1.15 R     +4.85       -0.01
        1.15-1.30 R     +2.18       -0.04
        1.70-2.50 R     -0.28       -0.10

    and the structure survives: rms kept is 86%, 75%, 88%, 100% in those
    shells, and 100% on the 2 px mask, which barely had the artifact.

    Two other fixes were tried and are worse. Flattening the image by its
    radial profile BEFORE filtering (what MGN does) gets the ratio to about
    -1.8, four times better than nothing but still a visible ring. Odd
    reflection of the polar image at the disc edge, so the kernel is two-sided,
    is far worse at -3.4 to -6.6: it turns the residual structure at the
    boundary into a mirrored copy of itself.

    Why this one is safe: the artifact is a function of radius ALONE, and so is
    what this removes. Hill's masks exist to carry the azimuthal structure --
    streamers, threads -- which is untouched. A genuinely radially-symmetric
    feature is not a streamer, and the radial profile is already the envelope's
    job and the radial-flatten control's.
    """
    nb = int(r.max()) + 2
    prof = _radial_median(m, r, valid, nb)
    return (np.asarray(m, np.float32)
            - np.interp(r, np.arange(nb, dtype=np.float32), prof
                        ).astype(np.float32))


def build_hill(wd, progress, lum_dn=None, disc=None, prom=None,
               cy=None, cx=None, frac=None, denoise="fine"):
    """Hill's unsharp-mask set, written to hill.npy / hill_log.npy.

    Called with everything already in hand during a stack, and with nothing but
    a work directory afterwards -- so a folder stacked before 0.22.64 can gain
    the layer without re-stacking. The rebuild path re-reads hdr_lum.npy and
    redoes the denoise, which is seconds; the polar blur is the whole cost.

    Built on the LOG-MAPPED image, not on log10 luminance: Hill's masks are
    differences of im_log and it is im_log they are added back to, so any other
    base would be a different filter wearing the same name.

    Stored as float16. These are high-pass residuals of a 0..1 image with an
    amplitude of about a thousandth, so float16 carries them to four figures and
    the set costs 2 bytes a pixel a scale instead of 4.

    Returns the stats dict, or None if it could not be built.
    """
    try:
        geo = json.load(open(os.path.join(wd, "geometry.json")))
        if cy is None or cx is None:
            cy, cx = geo["cy"], geo["cx"]
        R = geo["R"]
        # ECLIPSEFORGE_HILL_NODENOISE=1: build the masks on the RAW merged
        # luminance instead of the denoised master.
        #
        # WHY THIS IS A QUESTION AT ALL. Hill builds his masks on the log-mapped
        # HDR straight out of the merge. We build them on `lum_dn`, which the
        # "fine" profile has already soft-thresholded at 1.0-1.5 sigma in its two
        # finest starlet bands -- and then the renderer's Noise threshold slider
        # soft-thresholds the 2 px MASK again. Two soft thresholds in series on
        # the same spatial band do not add up to a stronger one: the first
        # removes the coefficients, the second removes what is left above them,
        # and what survives is sparse and isolated. That is the speckled look,
        # and it is why raising the slider takes structure with it.
        #
        # Not made the default, because "fine" denoise is also what keeps the
        # far field from tearing itself apart, and which way is better is a
        # picture judgement on a real bracket, not a number. One run each.
        _raw_base = os.environ.get("ECLIPSEFORGE_HILL_NODENOISE") == "1"
        if _raw_base:
            lum_dn = None       # force the reload below, ignoring what we were
            progress.log("  partial convolution: building on the RAW merged "
                         "luminance (ECLIPSEFORGE_HILL_NODENOISE=1) — the "
                         "denoise is skipped for the masks only", None)
        if lum_dn is None:
            lum = np.load(os.path.join(wd, "hdr_lum.npy"))
            ks = DENOISE_PROFILES.get(denoise, DENOISE_PROFILES["fine"])
            H, W = lum.shape
            yy = np.arange(H, dtype=np.float32)[:, None] - cy
            xx = np.arange(W, dtype=np.float32)[None, :] - cx
            rr = np.sqrt(yy * yy + xx * xx)
            if _raw_base:
                ks = (0.0, 0.0, 0.0, 0.0)
            if any(k > 0 for k in ks):
                progress.log("  partial convolution: re-deriving the denoised master", None)
                lum_dn = (10.0 ** denoise_loglum(
                    np.log10(np.clip(lum, 1.0, None)),
                    photon_floor(lum, rr), ks=ks)).astype(np.float32)
            else:
                lum_dn = lum
            del lum
        else:
            H, W = lum_dn.shape
            yy = np.arange(H, dtype=np.float32)[:, None] - cy
            xx = np.arange(W, dtype=np.float32)[None, :] - cx
            rr = np.sqrt(yy * yy + xx * xx)
        if disc is None:
            margin = float(geo.get("limb_margin",
                                   geo.get("Rmask", R + 4.0) - R))
            prof = geo.get("limb_prof")
            Rmap = (limb_radius_map(prof, (H, W), cy, cx, margin)
                    if prof else np.float32(R + margin))
            disc = rr < Rmap
        if prom is None:
            prom = prominence_mask(wd, geo, (H, W), rr, R, progress=progress)
        progress.log(f"partial-convolution unsharp masks at "
                     f"{', '.join('%g' % s for s in HILL_SCALES)} px "
                     f"(polar-oriented, partial)...", frac)
        hk = np.float32(10.0 ** float(
            os.environ.get("ECLIPSEFORGE_HILL_LOGK", "6")))
        _pk = max(float(np.percentile(lum_dn, 99.95)), 1e-6)
        xn = np.clip(lum_dn / _pk, 0, 1)
        imlog = (np.log1p(hk * xn) / np.log1p(hk)).astype(np.float32)
        # THE EXPECTED NOISE OF im_log, PER PIXEL, so the renderer can tell a
        # mask coefficient that is structure from one that is noise.
        #
        # photon_floor gives the sigma of log10(luminance) from the photon
        # model. im_log is a different curve, so the sigma has to be carried
        # through it:  d(im_log)/d(log10 L) = ln(10) * xn * k / ((1 + k*xn) *
        # ln(1+k)).  That derivative is the whole reason the far field looks
        # grainy -- it is about 800x larger out there than at the limb at
        # k = 1e3 -- and it is exactly what makes the threshold scale correctly
        # with radius without anyone tuning a radius profile.
        _sl = photon_floor(lum_dn, rr)
        _d = (np.log(10.0) * xn * hk
              / ((1.0 + hk * xn) * np.log1p(hk))).astype(np.float32)
        sig_log = (_sl * _d).astype(np.float32)
        del xn, _sl, _d
        w = (~disc).astype(np.float32)
        if prom is not None:
            w[prom] = 0.0
        blurs, covs = polar_partial_blur(imlog, w, cy, cx, HILL_SCALES,
                                         progress=progress)
        good = ~disc if prom is None else (~disc & ~prom)
        Ms, bias = [], []
        for b, cv in zip(blurs, covs):
            m = (imlog - b).astype(np.float32)
            bias.append(_mask_bias(m, rr, good, R))
            m = _deradial(m, rr, good)
            # HILL'S STEP 4: "replace any pixels restricted by the mask with
            # zero". Leaving it out is what put a row of dark blobs in an arc
            # below the Moon on the 250 mm test set.
            #
            # Where the mask is 0 the partial convolution has no data to work
            # with, so `blur` is whatever num/max(den, 1e-4) happens to produce
            # and the difference from it is meaningless. Measured on a fixture
            # with a prominence mask like the real one, |M| inside the masked
            # patches against |M| over the corona:
            #
            #     2 px  x174    4 px  x91    8 px  x14    16 px  x12    32 px  x9
            #
            # and those patches sit out to 1.35 R, well OUTSIDE the disc mask,
            # so they reach the picture. At 2 px the blobs are two orders of
            # magnitude above the structure the layer is supposed to carry.
            #
            # The taper matters as much as the zeroing: a hard cut leaves the
            # RIM of every patch, where coverage is low but not zero and the
            # division amplifies whatever little it saw. Fading from 0.35 to
            # 0.7 coverage takes the rim with it.
            _t = np.clip((np.asarray(cv, np.float32) - 0.35) / 0.35, 0.0, 1.0)
            m *= (_t * _t)
            m[~good] = 0.0
            Ms.append(m.astype(np.float16))
            del _t, m
        M = np.stack(Ms)
        del blurs, covs, Ms
        bias2 = [_mask_bias(np.asarray(M[i], np.float32), rr, good, R)
                 for i in range(len(HILL_SCALES))]
        # The per-scale response of the polar high-pass to UNIT white noise,
        # measured on a probe rather than assumed, so the threshold is in real
        # sigma units. Same trick denoise_loglum uses for its starlet levels.
        _pr = np.random.default_rng(4242).standard_normal((256, 256)).astype(np.float32)
        _pw = np.ones_like(_pr)
        _pb, _ = polar_partial_blur(_pr, _pw, 128.0, 128.0, HILL_SCALES, band=512)
        resp = [float(np.std((_pr - b)[32:-32, 32:-32])) for b in _pb]
        del _pr, _pw, _pb
        np.save(os.path.join(wd, "hill.npy"), M)
        np.save(os.path.join(wd, "hill_log.npy"), imlog.astype(np.float32))
        np.save(os.path.join(wd, "hill_sigma.npy"), sig_log.astype(np.float32))
        _rms = [float(np.std(np.asarray(M[i], np.float32)[good]))
                for i in range(len(HILL_SCALES))]
        # THE STRUCTURE PART OF EACH MASK'S SPREAD, separated from the noise
        # part, because the renderer divides the Amplification slider by it.
        #
        # `_rms` is the mask's total spread over the corona, and on a short or
        # noisy bracket most of that IS noise -- the 2 px mask especially, which
        # is the one the slider is normalised to. Dividing the gain by a number
        # made of noise means "Amplification 0.1" buys ten percent of the
        # display range PER SIGMA OF NOISE, and the same slider setting delivers
        # a quiet picture on a clean set and a grainy one on a noisy set. That
        # is backwards: the control should hold the amount of STRUCTURE steady.
        #
        # Both terms are already measured. `resp[i]` is this scale's response to
        # unit white noise through the same polar high-pass, and `sig_log` is the
        # per-pixel sigma of im_log from the photon model, so the noise the mask
        # is expected to carry is resp[i] * rms(sig_log), and
        #
        #     structure^2 = total^2 - noise^2
        #
        # with a floor at 25% of the total so a mask that measures as pure noise
        # cannot drive the gain to infinity. Old workdirs have no such figure and
        # the renderer falls back to `rms`, unchanged.
        _sg = float(np.sqrt(np.mean(np.square(
            np.asarray(sig_log, np.float32)[good]))))
        _rms_s = []
        for i, _t in enumerate(_rms):
            _nz = float(resp[i]) * _sg if i < len(resp) else 0.0
            _rms_s.append(float(max((_t * _t - _nz * _nz) ** 0.5
                                    if _t > _nz else 0.0, 0.25 * _t)))
        st = {"scales": [float(s) for s in HILL_SCALES],
              "logK": float(np.log10(float(hk))),
              "build": HILL_BUILD,
              "resp": resp,
              "base": "raw" if _raw_base else "denoised",
              "sigma_rms": _sg,
              "limb_bias_before": bias, "limb_bias_after": bias2,
              "rms": _rms, "rms_struct": _rms_s}
        progress.log("  masks: per-scale rms over the corona "
                     + ", ".join("%g px %.2e" % (s, v)
                                 for s, v in zip(HILL_SCALES, st["rms"])), None)
        progress.log(
            "  of which structure rather than photon noise (the Amplification "
            "slider is normalised to the 2 px figure) "
            + ", ".join("%g px %.0f%%" % (s, 100.0 * b / max(a, 1e-12))
                        for s, a, b in zip(HILL_SCALES, _rms, _rms_s)), None)
        progress.log(
            "  limb bias (the mask's own mean at 1.05-1.3 R, in units of its "
            "structure -- this is the featureless collar) "
            + ", ".join("%g px %+.2f->%+.2f" % (s, a, b)
                        for s, a, b in zip(HILL_SCALES, bias, bias2)), None)
        del M, imlog, w, sig_log
        return st
    except Exception as e:
        progress.log(f"partial-convolution masks skipped ({e})", None)
        for f in ("hill.npy", "hill_log.npy", "hill_sigma.npy"):
            p = os.path.join(wd, f)
            if os.path.exists(p):
                try:
                    os.remove(p)
                except OSError:
                    pass
        return None


def build_layers(wd, progress, denoise="fine", earthshine=False,
                 fnrgf_preset="ours"):
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
    # The partial-convolution mask: the disc, and the prominences too. See
    # prominence_mask -- +27% and +25% of near-limb coronal structure on the
    # rays that have a prominence, nothing measurable on the rays that do not.
    #
    # SCOPED DELIBERATELY. `valid` stays the disc-only mask, because every later
    # user of it -- the deband trend, the inner layers -- was measured
    # against that mask and none of them was measured against this one. Only
    # MGN's own statistics get the tighter mask.
    #
    # Consequence worth stating: a masked pixel's normalized value is undefined
    # (its own brightness is excluded from the mean it would be compared to), so
    # the prominence cores come out flat 0.5 in THIS layer -- Hill's step 4,
    # "replace any pixels restricted by the mask with 0". Prominence structure
    # is carried by the promdet layer and the prominence gate, which exist for
    # exactly that, so nothing is lost; but it does mean the MGN layer viewed on
    # its own now has holes where the prominences are.
    _pm = prominence_mask(wd, geo, (H, W), r, R, progress=progress)
    valid_stat = valid if _pm is None else (valid & ~_pm)
    if _pm is not None:
        progress.log(f"  prominences masked out of the convolution too "
                     f"({_pm.sum() / 1e3:.0f}k px)", None)
        lstats["prom_masked_px"] = int(_pm.sum())
    mgl = mgn(Lf, floor_map=nf, valid=valid_stat, norm_span=(-half, half), scales=_sc)
    mgl = _deband(mgl, r, valid_stat, cy, cx, R + margin)
    np.save(os.path.join(wd, "mgn.npy"), mgl.astype(np.float32))
    del mgl

    # SECOND MGN OVER THE THREE FINE SCALES ONLY.
    #
    # The gains above are Morgan & Druckmuller's NOISE-NORMALISATION constants
    # (their Fig. 4) -- 0.907 at w=1.25 rising to 1.0 by w=5 -- which correct for
    # a small kernel measuring a smaller local sigma on pure noise. They are not
    # amplification factors, and their effect is that all six scales get equal
    # say, tilted very slightly toward the COARSE end.
    #
    # Hill combines his masks as E = 100A + 60B + 20C + 10D (A=2px..D=16px):
    # normalised, 1 : 0.6 : 0.2 : 0.1, the finest scale weighted TEN TIMES the
    # coarsest. That is the opposite emphasis, and side by side on the reference
    # bracket it is the difference between slabs and filaments -- the coarse
    # scales dominate the layer's range, so fine structure arrives already
    # compressed against the rails.
    #
    # amp*coh RATES HILL'S LADDER LOWER in every shell (0.0535 -> 0.0300 at the
    # limb) and it is still the better picture. The score rewards radially
    # coherent structure and the coarse scales carry plenty of it; "delicate" is
    # not a thing it measures. This is an emphasis choice, so it is a slider.
    #
    # Storing the fine group separately makes any ladder that is constant within
    # a group a RENDER-TIME blend, because the layer is a weighted mean of its
    # per-scale terms. Measured, the blend reaches the true ladder: fine/coarse
    # energy ratio 1.09 as shipped, 3.18 at detailScale 1, 3.32 for the real
    # Hill gains -- and the two are indistinguishable side by side.
    progress.log("MGN, fine scales only (detail-balance slider)...",
                 _DF.get("mgn_fine", 99.0))
    _nf = len(_sc) // 2
    mgf = mgn(Lf, floor_map=nf, valid=valid_stat, norm_span=(-half, half),
              scales=_sc[:_nf], gains=(0.907, 0.976, 0.994)[:_nf])
    mgf = _deband(mgf, r, valid_stat, cy, cx, R + margin)
    np.save(os.path.join(wd, "mgn_fine.npy"), mgf.astype(np.float32))
    lstats["mgn_fine_scales"] = [round(x, 2) for x in _sc[:_nf]]
    del nf, L, Lf, mgf

    progress.log("FNRGF detail extraction...", _DF["fnrgf"])
    # EVERY ARGUMENT HERE WAS A DEFAULT UNTIL 0.23.2, INCLUDING THREE THAT
    # IMPLEMENT PUBLISHED BEHAVIOUR WE DOCUMENTED AND NEVER RAN.
    #
    #   mask        the prominence. Thesis p.81 names a prominence as one of two
    #               causes of this filter's own ring artifact. We already build
    #               the mask for the partial convolution.
    #   impulse     the other named cause. An azimuthal median on the FIT
    #               samples; her own remedy is a median on the whole input,
    #               which would cost resolution we do not need to spend.
    #   noise_frac  eqs. 6.11-6.13. sqrt(Vn) at 10-15% of the measured outer
    #               noise, ADDED to the local variance. The alternative our
    #               `else` branch takes is a relative floor that rises with the
    #               signal in the inner corona and divides real structure out.
    #
    # `noise_frac` is held at her lower bound: we denoise the master before this
    # runs, so the residual noise reaching the filter is smaller than the thesis
    # assumes and the 10-15% band is an upper estimate for us, not a target.
    _fnp = str(fnrgf_preset or "ours").strip().lower()
    if _fnp == "published":
        # Thesis p.80: "Setting omega = 30 and the attenuation coefficients Ak
        # linearly decreasing with step 0.05 and Ck by 0.1 is a setting that
        # always works without producing any artifacts."
        _fn_kw = dict(order=30, atte_ave=0.05, atte_dev=0.1)
    else:
        _fn_kw = dict(order=6, atte_ave=0.0, atte_dev=0.0)
    D = fnrgf_robust(lum_dn, r, cy, cx, int(R) + 4,
                     mask=_pm, impulse=FNRGF_IMPULSE,
                     noise_frac=FNRGF_NOISE_FRAC, **_fn_kw)
    _nl = getattr(fnrgf_robust, "last_noise", None)
    _nm = getattr(fnrgf_robust, "last_masked", 0)
    progress.log(
        f"  FNRGF: order {_fn_kw['order']}"
        + (f", attenuation {_fn_kw['atte_ave']:.2f}/{_fn_kw['atte_dev']:.2f} per "
           f"harmonic" if _fn_kw["atte_ave"] else ", hard order cutoff")
        + (f"; prominence samples dropped from {_nm} ring(s)" if _nm else "")
        + (f"; added noise sqrt(Vn) {FNRGF_NOISE_FRAC:.2f}x the outer-field "
           f"residual ({_nl[0]:.4f} dex over {_nl[2]} rings)"
           if _nl and np.isfinite(_nl[0]) else ""), None)
    lstats["fnrgf"] = {"preset": _fnp, "order": _fn_kw["order"],
                       "atte_ave": _fn_kw["atte_ave"],
                       "atte_dev": _fn_kw["atte_dev"],
                       "impulse": FNRGF_IMPULSE,
                       "noise_frac": FNRGF_NOISE_FRAC,
                       "masked_rings": int(_nm),
                       "sigma_outer": (round(_nl[0], 5)
                                       if _nl and np.isfinite(_nl[0]) else None)}
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
        # DEPARTURE FROM THE PAPERS, STATED AS ONE. Both NAFE papers open their
        # method section with "In this image the pixel values have a LINEAR
        # dependence on the intensity of the coronal emission" (2013 p.2, 2014
        # p.265) and never revisit it. We feed LOG luminance, and we high-pass
        # it first. Neither has any counterpart in either paper -- there is no
        # gradient removal, no log, and no lunar masking in their procedure;
        # NRGF and FNRGF appear in their introductions as ALTERNATIVE methods,
        # never as a NAFE pre-step. Both departures are defensible here (see the
        # measurements below and in nafe.py) but they are ours, and they mean no
        # published NAFE parameter transfers to this code by units alone.
        #
        # TESTED. All four combinations, same filter, reference bracket:
        #
        #   input                          1.05-1.30 R     1.30-1.80 R    1.80-2.60 R
        #   log + flatten (this)           0.191 (c0.97)   0.124 (c0.86)  0.079 (c0.53)
        #   log, no flatten                0.151 (c0.95)   0.106 (c0.84)  0.070 (c0.53)
        #   linear + flatten               1.120 (c0.85)   0.704 (c0.83)  0.133 (c0.58)
        #   linear, no flatten (as spec'd) 0.018 (c0.44)   0.000 (c0.81)  0.001 (c0.45)
        #
        # The papers' literal form COLLAPSES here -- 0.018 at the limb and
        # essentially zero beyond it. Their assumption is AIA data, whose
        # dynamic range is a few hundred to one; an eclipse corona spans 7 EV,
        # so a value window of fixed width on a linear axis either contains the
        # whole neighbourhood or none of it, depending on radius.
        #
        # `linear + flatten` scores six times higher than what we ship and is
        # the trap of the day: 1% of its pixels are pinned at the floor and the
        # layer is clipped to black and white. Variance is maximised by
        # binarising, which is why amp*coh rewards it.
        #
        # So both departures are justified, and now measured rather than
        # assumed. Log is what makes eps mean the same thing at every radius.
        #
        # Flatten the envelope first -- see NAFE_FLATTEN_R. The local mean is
        # built by normalized convolution over non-disc pixels only: a plain
        # Gaussian straddling the dark disc is what produced the black annulus.
        #
        # WHERE THE KERNEL CANNOT REACH. `gaussian_filter` truncates at 4 sigma
        # and runs separably, so its support is a SQUARE of half-width 4*_s --
        # 198 px at the reference set's geometry against a disc 618 px in radius. Past that
        # distance from the mask edge, `_den` is not small, it is EXACTLY 0, and
        # so is `_num`. Clamping _den and dividing then handed NAFE the disc
        # interior UNFLATTENED, on which a rank equaliser saturates: measured on
        # the 600 mm run as a plateau at the layer's exact global maximum
        # (0.880797), bounded by the disc eroded by that square -- a rounded
        # square with corners on the axes, which is the "diamond" in the disc.
        # Predicted edge 424 px on the axes / 398 on the diagonals against 408 /
        # 365 measured. Masked out of the picture, but the layer is wrong there
        # and the fallback costs nothing: where the kernel reaches nothing, use
        # the global mean of the non-disc region, which is the best available
        # estimate of what the local mean would have been.
        _s = max(NAFE_FLATTEN_R * R, 4.0)
        _w = (~disc_m).astype(np.float32)
        _num = ndimage.gaussian_filter(Ldn * _w, _s)
        _den = ndimage.gaussian_filter(_w, _s)
        _glob = float(np.sum(Ldn * _w) / max(float(np.sum(_w)), 1.0))
        del _w
        # 1e-3 not 1e-6: below about a thousandth of full weight the quotient is
        # a handful of far-field pixels amplified by three orders of magnitude,
        # which is noise, not a local mean.
        _thin = _den < 1e-3
        np.maximum(_den, 1e-6, out=_den)
        _num /= _den
        del _den
        _nthin = int(_thin.sum())
        if _nthin:
            _num[_thin] = _glob
        del _thin
        _Lnf = (Ldn - _num).astype(np.float32)
        del _num
        if _nthin:
            progress.log(f"  NAFE flatten: {100.0 * _nthin / _Lnf.size:.1f}% of "
                         f"the frame lies further than {4 * _s:.0f}px from open "
                         f"sky (deep inside the disc); the global mean is used "
                         f"there rather than zero", None)
        nv = nafe_vn(_Lnf, K=NAFE_K, gamma=NAFE_GAMMA, combine=False,
                     sigma_sp=max(NAFE_NEIGH_R * R, 1.0),   # nafe_vn divides by `grid` itself
                     noise_mult=NAFE_NOISE_MULT,
                     eps_frac=NAFE_EPS, kernel="gauss", grid=NAFE_GRID)
        np.save(os.path.join(wd, "nafe.npy"), nv.astype(np.float32))
        del _Lnf
        lstats["nafe"] = {"K": NAFE_K, "eps": NAFE_EPS, "layer": "E",
                          "flatten_px": round(max(NAFE_FLATTEN_R * R, 4.0), 1),
                          "flatten": "normalized convolution, disc excluded",
                          "noise_mult": NAFE_NOISE_MULT,
                          "neigh_px": round(max(NAFE_NEIGH_R * R, 1.0), 1),
                          }
        del nv
    except Exception as e:
        progress.log(f"NAFE layer unavailable ({e})", None)
        np.save(os.path.join(wd, "nafe.npy"), np.full((H, W), 0.5, np.float32))
    # HILL'S UNSHARP-MASK SET. Built HERE, before `lum_dn` is released:
    # calling it after that `del` raised UnboundLocalError at the end of an
    # 18-minute stack, which is the worst possible place to find a typo.
    # Everything it needs is also on disk, so
    # it is written as a function that can also be called on an ALREADY STACKED
    # work directory -- see build_hill. The same argument build_hill makes:
    # it would be absurd to charge someone a 13-minute re-stack for a layer that
    # is one transform of data they already have.
    lstats["hill"] = build_hill(wd, progress, lum_dn=lum_dn, disc=disc_m,
                                prom=_pm, cy=cy, cx=cx, frac=_DF["pellett"])

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
        # Measured on an imported HDR (a third tester's Siril stack of 327 frames):
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

    progress.log("tangential filter...", _DF["pellett"])
    _pellett(wd, lum, r, cy, cx, R, disc_m)

    ep = os.path.join(wd, "earth.npy")
    if earthshine:
        progress.log("earthshine layer...", _DF["earth"])
        _earthshine(wd, r, cy, cx, R)
    elif os.path.exists(ep):
        os.remove(ep)
    return lstats


def _pellett(wd, lum, r, cy, cx, R, disc, blur_deg=6.0, na=2880):
    """Tangential unsharp: subtract a rotational blur about the disc centre
    from the log luminance. Softer texture than MGN/FNRGF.

    IT IS BLIND TO TANGENTIAL STRUCTURE, BY CONSTRUCTION, and that is not a
    quirk of this implementation -- it is the documented property of the whole
    method. Druckmullerova's thesis, section 5.1, on exactly this filter:

        "since the filter averages only in the tangential direction (around the
        Sun), it is direction-dependent. Structures that are oriented radially,
        such as plumes, are enhanced the most, structures that are oriented
        tangentially, i.e. parallel with the edge of the Sun (such as tops of
        loops and helmet streamers), are not enhanced at all."

    and Druckmuller, Rusin & Minarovjech (2006) section 4.1 on the same thing:
    "it is blind to tangential structures".

    So a user who raises this layer and sees plumes sharpen while the tops of
    loops do not is seeing the method work as designed, not a bug here.

    THE THESIS IS HARSHER STILL, and it is worth knowing why the layer is kept
    small rather than trusted. A moving average has a sinc frequency response,
    so its amplification is not monotone in frequency -- some scales are
    enhanced and some attenuated -- and its phase spectrum runs over the whole
    of [0, 2pi), which means structures can be SHIFTED. The thesis concludes
    that images filtered this way "are not suitable for further scientific
    interpretation", and records that this filter once had people doubting
    correct magnetic-field models, because it hid the closed loops those models
    predicted: "the eclipse images with the Espenak or similar filter applied
    did not show the reality, they were physically incorrect and the models
    were correct."

    The layer stays because it gives a texture the others do not, and because
    every claim above is about what it FAILS to show rather than about it
    inventing something. It is a look, not a measurement, and the caller should
    treat it that way.

    ATTRIBUTION IS UNVERIFIED. The report credits "Druckmuller 2009, ApJ 706,
    1605, sec. 5.1, after Pellett". ApJ 706, 1605 is the phase-correlation
    paper, and neither the 2006 paper nor the thesis attributes this filter to
    Pellett -- both trace it to Espenak (2000), from Photoshop's Radial Blur /
    Spin. That paper is not in the literature folder, so this is flagged rather
    than corrected.
    """
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
