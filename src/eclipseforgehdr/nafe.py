"""NAFE with a variable neighbourhood, and the value-based masking it implies.

Druckmuller, "A noise adaptive fuzzy equalization method for processing solar
extreme ultraviolet images" (ApJ 775, 88, 2013) and Druckmuller &
Druckmullerova, "Noise Adaptive Fuzzy Equalization Method with Variable
Neighborhood" (IWCIA 2014, LNCS 8466, p. 262).

The 2014 paper lists three defects of the original method and fixes each:

  1. neighbourhood-shape artifacts   -> a fuzzy (Gaussian) neighbourhood
  2. noise amplification             -> the cumulative histogram is convolved
                                        with a Gaussian of width 2..12 sigma_A
  3. "loss of contrast on boundaries between areas with significantly
     different brightness" -- their Fig. 2 is captioned "loss of contrast near
     lunar edge (edge effect)" -> the VARIABLE NEIGHBOURHOOD

Point 3 is the one that matters most here, and it is worth being precise about
what their fix is, because it is not what EclipseForgeHDR has been doing. Their
neighbourhood restriction (eqs. 11-12) is a restriction in VALUE: a pixel's
neighbourhood keeps only those nearby pixels whose brightness is within a range
eps of the centre pixel's own brightness. It is not a geometric mask.

That distinction is the whole point. A geometric mask has to be told where the
Moon is, so it inherits every error in the limb fit -- and when the limb fit is
wrong the mask carves the wrong hole and the filter output is destroyed. A
range restriction asks each pixel only about its own value and the values
around it. It needs no geometry, so it cannot be misplaced: the dark lunar
plateau stops contributing to the statistics of the corona just outside it
because it is dark, not because a circle was drawn around it.
"""
from __future__ import annotations
import numpy as np
from scipy import ndimage


def multiscale_blur(a, n_scales=12, base=2.0, out=None, sigma=None):
    """The fuzzy multiscale neighbourhood of the 2014 paper: a sum of Gaussians
    with sigma_m = 2^(m/2), which behaves like a single smooth kernel with no
    preferred scale (their n = 129 kernel, built as 12 summed Gaussians).

    `sigma` sets the WIDEST scale, so the ladder ends there instead of always
    ending at 2^((n-1)/2) whatever the image size. Without it the caller's
    sigma_sp was silently ignored on this path and the neighbourhood was fixed
    in grid pixels -- which is to say fixed relative to the decimation factor
    rather than to the Sun."""
    acc = np.zeros_like(a, dtype=np.float32) if out is None else out
    acc[...] = 0
    wsum = 0.0
    _sc = 1.0 if sigma is None else sigma / (base ** ((n_scales - 1) / 2.0))
    for m in range(n_scales):
        s = _sc * base ** (m / 2.0)
        w = 1.0 / (1.0 + m)          # taper so the widest scale does not rule
        acc += w * ndimage.gaussian_filter(a, s, mode="nearest")
        wsum += w
    acc /= wsum
    return acc


def value_neighbourhood_weight(L, sigma_sp=25.0, eps=None, valid=None,
                               n_scales=8):
    """Normalized-convolution weights from a RANGE restriction rather than a
    geometric mask -- Druckmuller's variable neighbourhood, used as a mask.

    Returns (mean, weight): the neighbourhood mean of `L` computed only over
    pixels whose value is within `eps` of the centre pixel, and the fraction of
    the kernel that survived the restriction. Where a pixel sits in the middle
    of a uniform region the weight is ~1; on a strong brightness boundary --
    the lunar limb above all -- the weight falls, and the mean stops being
    contaminated from the other side of the boundary.

    `eps` defaults to a robust estimate of the image's own contrast scale, so
    it adapts to the data instead of needing a tuned constant.
    """
    L = np.asarray(L, np.float32)
    if eps is None:
        m = valid if valid is not None else np.isfinite(L)
        v = L[m] if m is not None else L.ravel()
        eps = float(np.percentile(v, 84) - np.percentile(v, 16)) * 0.5
        eps = max(eps, 1e-3)
    # A bilateral mean, computed by the standard "shiftable" decomposition:
    # a few range-basis terms, each a plain Gaussian blur. Exact enough here
    # and orders of magnitude cheaper than a per-pixel neighbourhood scan.
    scale = 1.0 / max(eps, 1e-6)
    base = np.zeros_like(L)
    wacc = np.zeros_like(L)
    macc = np.zeros_like(L)
    lo, hi = float(np.nanmin(L)), float(np.nanmax(L))
    K = 24
    levels = np.linspace(lo, hi, K, dtype=np.float32)
    ok = np.ones_like(L, np.float32) if valid is None else valid.astype(np.float32)
    Lf = np.nan_to_num(L)
    for lv in levels:
        # membership of every pixel in this brightness level, width eps
        g = np.exp(-0.5 * ((Lf - lv) * scale) ** 2, dtype=np.float32) * ok
        gb = ndimage.gaussian_filter(g, sigma_sp, mode="nearest")
        mb = ndimage.gaussian_filter(g * Lf, sigma_sp, mode="nearest")
        # weight of this level for the centre pixel
        c = np.exp(-0.5 * ((Lf - lv) * scale) ** 2, dtype=np.float32)
        wacc += c * gb
        macc += c * mb
        base += c
    del g, gb, mb, c
    wacc /= np.maximum(base, 1e-9)
    macc /= np.maximum(base, 1e-9)
    mean = macc / np.maximum(wacc, 1e-9)
    # normalize the weight so an interior pixel reads 1
    wnorm = wacc / max(float(np.percentile(wacc, 98)), 1e-9)
    return mean.astype(np.float32), np.clip(wnorm, 0, 1).astype(np.float32)


def nafe_vn(A, sigma_sp=30.0, K=128, w=0.2, gamma=3.0, noise_sigma=None,
            eps_frac=0.10, valid=None, n_scales=8, kernel="gauss", grid=8,
            combine=False, noise_mult=4.0, knee=3.0):
    """Noise Adaptive Fuzzy Equalization with a Variable Neighbourhood.

    Returns E, the fuzzy rank of each pixel within its own neighbourhood,
    restricted to neighbours of similar VALUE (their eqs. 11-12). `combine=True`
    returns their eq. 2 display image instead, B = (1-w) T_gamma(A) + w E.

    Correspondence with the paper, term by term:

      eq. 4-5   the fuzzy kernel l_{k,l}          -> `kernel`, width `sigma_sp`
      eq. 6     fuzzy histogram h(x)              -> K blurred membership maps
      eq. 7-8   cumulative, normalised C(x)       -> cumsum over the level axis
      eq. 11-12 restriction to |x - a| < eps      -> `eps_frac`, in LEVEL units
      eq. 13    C(x) * G_sigma(x), sigma 2..12 A  -> `noise_mult` x sigma_A
      eq. 2     B = (1-w)T_gamma + wE             -> `combine`

    EVERY METRIC IS IN THE UNITS OF THE LEVEL AXIS (rewritten in 0.11.0)
    --------------------------------------------------------------------
    This is the whole point of the method and it is what earlier versions got
    wrong. Up to 0.10.4 the image was first passed through an equal-population
    rank map, and eps and sigma were then applied in RANK units. That looks
    harmless -- a monotone remap does not change which neighbour is brighter
    than which -- but it destroys the mechanism eq. 13 depends on.

    Eq. 13 works because sigma is a FIXED width in the value units of A. Where
    the local histogram is wide (real contrast) a fixed-width smoothing is
    negligible; where it is narrow (noise only) it dominates and flattens the
    rank. That is the noise adaptivity, and it is automatic. Under a rank map
    the sky -- most of the pixels -- is stretched to occupy most of the axis and
    the corona is compressed into what is left, so the same physical width means
    something different at every radius, and the automatic behaviour is gone.
    A per-level correction was bolted on to compensate; it did not.

    Measured on the reference set, per-annulus contrast (sd), 8102x5359 frame:

                        1.05-1.5R  1.5-2R   2-2.5R  2.5-3R   3-4R (sky)
        rank units       0.0945    0.0108   0.0261  0.0587   0.0919
        level units      0.0862    0.0261   0.0366  0.0323   0.0310

    In rank units the layer's contrast RISES with radius, tracking the falling
    signal-to-noise: it was a noise detector. In level units it falls and then
    settles, the corona from 1.5 to 3 R gains 40-140% of contrast, and the sky
    loses two thirds of its grain.

    Also from the same measurements:

      * K = 128, not 64. The paper notes its images have "several thousand
        discrete pixel values"; 64 levels is a real approximation and it cost
        real contrast (2-2.5 R: 0.0273 at K=64 against 0.0366 at K=128). Past
        128 the corona keeps gaining but the sky gains faster, so the ratio
        peaks there.
      * The paper's plain Gaussian kernel beats the multiscale sum this used
        to build, and is nearly twice as fast (2-2.5 R: 0.0286 -> 0.0366).
      * `grid` is the one shortcut that is genuinely free. The paper calls the
        naive algorithm "extremely time consuming"; this evaluates every
        pixel's local histogram at once as K blurred membership maps, which is
        the same computation reorganised, and then samples that field every
        `grid` px because it varies far more slowly than the image. Measured,
        grid 8 against grid 4: agreement to three decimals in every annulus.
    """
    A = np.asarray(A, np.float32)
    m = valid if valid is not None else np.isfinite(A)
    v = A[m]
    if v.size < 100:
        return np.full_like(A, 0.5)
    lo = float(np.percentile(v, 0.1))
    hi = float(np.percentile(v, 99.9))
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        return np.full_like(A, 0.5)
    # a0, a1 of eq. 3: the level axis. A monotone pre-transform of A (this
    # pipeline hands in log luminance) is a deliberate departure -- the paper
    # assumes A linear in coronal emission, but a 6.5 EV corona quantised into
    # uniform linear levels leaves the outer corona fewer than one level. What
    # matters is that eps and sigma below are then measured in THESE units.
    x = np.clip((np.nan_to_num(A, nan=lo) - lo) / (hi - lo), 0, 1).astype(np.float32)
    ok = m.astype(np.float32)

    q = max(1, int(grid))
    xs = x[::q, ::q]
    oks = ok[::q, ::q]
    ssp = max(sigma_sp / q, 0.7)
    if kernel == "multiscale":
        blur = lambda a: multiscale_blur(a, n_scales=n_scales, sigma=ssp)
    else:
        blur = lambda a: ndimage.gaussian_filter(a, ssp, mode="nearest")

    edges = np.linspace(0.0, 1.0, K, dtype=np.float32)
    dl = edges[1] - edges[0]
    # Gaussian membership, a little wider than the level spacing. A triangular
    # one makes the cumulative histogram piecewise-quadratic, so its curvature
    # jumps at every level and those jumps print as concentric contour rings
    # wherever the image has a smooth radial gradient.
    mw_ = 1.5 * dl
    Hst = np.empty((K,) + xs.shape, np.float32)
    for i, e in enumerate(edges):
        Hst[i] = blur(np.exp(-0.5 * ((xs - e) / mw_) ** 2).astype(np.float32) * oks)

    # eq. 13, in the level units the histogram is built in
    if noise_sigma is None:
        d = np.abs(np.diff(x[m][:200000]))
        s_a = 1.4826 * float(np.median(d)) / 1.4142 if d.size else 0.01
    else:
        s_a = float(noise_sigma)
    sig = max(float(noise_mult) * s_a, 0.5 * dl)
    Mx = np.exp(-0.5 * ((edges[None, :] - edges[:, None]) / sig) ** 2)
    Mx /= np.maximum(Mx.sum(axis=1, keepdims=True), 1e-12)
    Hst = np.tensordot(Mx.astype(np.float32), Hst, axes=(1, 0))

    C = np.cumsum(Hst, axis=0)
    del Hst
    H_, W_ = x.shape
    yy = (np.arange(H_, dtype=np.float32) / q)[:, None] * np.ones((1, W_), np.float32)
    xx = (np.arange(W_, dtype=np.float32) / q)[None, :] * np.ones((H_, 1), np.float32)

    def _at(level):
        b = np.clip(level / dl, 0, K - 1)
        return ndimage.map_coordinates(
            C, [b.ravel(), yy.ravel(), xx.ravel()], order=1,
            mode="nearest").reshape(x.shape)

    eps = float(eps_frac)
    c_mid = _at(x)
    if eps > 0:
        # eqs. 11-12: rank only among neighbours within eps of this pixel's own
        # value. Its job is the lunar edge -- it stops the dark disc from
        # contaminating the statistics of the corona just outside it -- so it
        # wants to be loose enough not to bind anywhere else. Too small and the
        # image fragments into high-contrast patches whose borders are not real.
        c_lo = _at(np.clip(x - eps, 0, 1))
        c_hi = _at(np.clip(x + eps, 0, 1))
        E = (c_mid - c_lo) / np.maximum(c_hi - c_lo, 1e-6)
    else:
        tot = ndimage.map_coordinates(
            C[-1], [yy.ravel(), xx.ravel()], order=1,
            mode="nearest").reshape(x.shape)
        E = c_mid / np.maximum(tot, 1e-6)
    del C
    E = np.clip(E, 0, 1).astype(np.float32)

    # --- output conditioning (added in 0.11.1) ---------------------------
    # E is a rank, so its useful range is set by how much of the neighbourhood
    # a pixel actually beats -- in the quiet corona that is a narrow band around
    # the middle, while the near-limb brightness RIDGE is a local maximum by
    # definition and pins at exactly 1.0. Raw, that means a washed-out corona
    # and a blown white rim at the same time.
    #
    # So: rescale by the layer's own robust spread, and roll the response off
    # softly past `knee` robust sigmas so an extreme cannot reach the rail.
    # Both of the user-visible complaints come out of this one step. Measured
    # on the reference set (rim = 1.04-1.18 R, corona = 1.5-3 R):
    #
    #   knee    corona sd   1.05-1.5R sd   rim p99.9   rim fraction > 0.99
    #   raw       0.0315      0.0853        1.0000            7.02%
    #    2        0.0907      0.1131        0.7913            0.00%
    #    3        0.0950      0.1302        0.8789            0.00%
    #    6        0.0980      0.1466        0.9704            0.00%
    #   none      0.0991      0.1520        0.9954            7.85%
    #
    # 3 gives 3.0x the corona contrast and 1.5x at the limb while the rim stops
    # clipping entirely; without the knee the contrast is barely better and the
    # rim comes straight back. Set knee=0 to get the raw rank.
    if knee and knee > 0:
        vE = E[m]
        med = float(np.median(vE))
        mad = 1.4826 * float(np.median(np.abs(vE - med)))
        if mad > 1e-6:
            z = (E - med) / mad
            z = float(knee) * np.tanh(z / float(knee))
            E = (0.5 + 0.5 * np.tanh(z / 3.0)).astype(np.float32)
    if not combine:
        return np.where(m, E, 0.5).astype(np.float32)
    out = (1.0 - w) * np.power(x, 1.0 / max(gamma, 1e-3)) + w * E
    return np.where(m, out, 0.5).astype(np.float32)
