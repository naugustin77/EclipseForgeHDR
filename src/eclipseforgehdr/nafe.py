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


def multiscale_blur(a, n_scales=12, base=2.0, out=None):
    """The fuzzy multiscale neighbourhood of the 2014 paper: a sum of Gaussians
    with sigma_m = 2^(m/2), which behaves like a single smooth kernel with no
    preferred scale (their n = 129 kernel, built as 12 summed Gaussians)."""
    acc = np.zeros_like(a, dtype=np.float32) if out is None else out
    acc[...] = 0
    wsum = 0.0
    for m in range(n_scales):
        s = base ** (m / 2.0)
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


def nafe_vn(A, sigma_sp=30.0, K=48, w=0.2, gamma=3.0, noise_sigma=None,
            eps_frac=0.25, valid=None, n_scales=10, fuzzy=True, grid=4,
            combine=False):
    """Noise Adaptive Fuzzy Equalization with a Variable Neighbourhood.

    `A` is a luminance image, ideally already in a log or gamma domain.

    Returns E, the fuzzy rank of each pixel within its own neighbourhood: the
    fraction of the neighbourhood that is darker than it. Because a rank is
    scale-free, faint structure is brought up to the same contrast as bright
    structure without any radial model at all -- which is why NAFE sees the
    inner corona and the outer corona in one pass.

    WHY E AND NOT B (fixed in 0.10.1)
    ---------------------------------
    The paper's output is  B = (1-w) T_gamma(A) + w E_{N,sigma}(A)  (eq. 2)
    with w in 0.05..0.3, and `combine=True` returns exactly that. But B is
    their FINAL DISPLAY IMAGE: T_gamma carries the large-scale brightness and
    E carries the structure, and at w = 0.2 the result is four fifths gamma
    transform by construction.

    This pipeline does not want a display image here. It wants a detail layer
    to mix against MGN and FNRGF, and its composite already supplies the
    large-scale brightness through the envelope -- which is the same role
    T_gamma plays in eq. 2. Returning B therefore added a SECOND copy of the
    base image into the detail term and diluted the other two layers with it.

    Measured on the reference set (decimated x4, K=64, w=0.2, gamma=2.4,
    eps=0.05):

        corr(B, T_gamma)          0.992      <- B is a gamma stretch
        corr(E, T_gamma)          0.562
        high-pass sd of B         0.0180
        high-pass sd of E         0.0675     <- 3.8x the local structure

    So the eq. 2 mix still happens; it happens in render.py, where the
    envelope is T_gamma and the nafeMix slider is w.

    Implementation: the local cumulative histogram is evaluated for every pixel
    at once by quantizing into K levels and blurring each level's membership
    (K blurs, not a per-pixel scan). Then

      * noise adaptivity  -- the cumulative histogram is smoothed ALONG THE
        LEVEL AXIS by `noise_sigma`, so differences smaller than the noise do
        not get stretched (their point 2);
      * variable neighbourhood -- the rank is taken only over levels within
        eps of the centre value, by evaluating the cumulative histogram at
        A-eps and A+eps and renormalizing between them (their eqs. 11-12).
        This is what removes the contrast loss at the lunar edge (their
        point 3), and it needs no knowledge of where the edge is.
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
    x = np.clip((np.nan_to_num(A, nan=lo) - lo) / (hi - lo), 0, 1)
    ok = m.astype(np.float32)

    # Equal-population binning. The corona spans some 12 EV, so on a linear or
    # even a log scale almost all of the outer corona falls inside a single
    # histogram bin and the local rank comes out flat -- the filter then does
    # nothing but a gamma stretch. Mapping through the GLOBAL cumulative
    # distribution first makes every bin equally populated, so the histogram
    # resolves structure at the faint end as finely as at the bright end. The
    # map is monotone, so it does not change what the local rank means.
    qs = np.linspace(0, 100, 512)
    knots = np.percentile(v, qs).astype(np.float64)
    knots = np.maximum.accumulate(knots + np.arange(knots.size) * 1e-9)
    x_r = np.interp(np.nan_to_num(A, nan=lo).ravel(), knots,
                    (qs / 100.0)).reshape(A.shape).astype(np.float32)
    x_lin, x = x, np.clip(x_r, 0, 1)

    # --- local histogram: K level memberships, each spatially blurred ---
    #
    # Built on a COARSE grid. The neighbourhood is tens of pixels across, so
    # the local histogram varies far more slowly than the image does; sampling
    # it every `grid` px and interpolating costs nothing in accuracy and turns
    # K x n_scales full-resolution blurs into K x n_scales small ones. This is
    # the difference between a filter that runs in seconds and one that does
    # not finish.
    q = max(1, int(grid))
    edges = np.linspace(0.0, 1.0, K, dtype=np.float32)
    dl = edges[1] - edges[0]
    xs = x[::q, ::q]
    oks = ok[::q, ::q]
    ssp = sigma_sp / q
    blur = (lambda a: multiscale_blur(a, n_scales=n_scales)) if fuzzy else \
        (lambda a: ndimage.gaussian_filter(a, ssp, mode="nearest"))
    Hst = np.empty((K,) + xs.shape, np.float32)
    # Gaussian, not triangular, membership. A triangular membership makes the
    # cumulative histogram piecewise-quadratic, so its curvature jumps at every
    # bin centre; under a strong stretch those jumps print as concentric
    # contour rings wherever the image has a smooth radial gradient. A Gaussian
    # membership a little wider than the bin spacing makes the whole cumulative
    # histogram smooth, and the rings go with it.
    mw_ = 1.5 * dl
    for i, e in enumerate(edges):
        Hst[i] = blur(np.exp(-0.5 * ((xs - e) / mw_) ** 2).astype(np.float32) * oks)
    # --- noise adaptivity: smooth the histogram along the LEVEL axis ---
    #
    # The paper's sigma is 2..12 times sigma_A, the noise width IN IMAGE VALUE
    # UNITS. It has to be measured there and then carried through the rank map,
    # not measured on the ranks: the rank map is steep wherever the histogram
    # is dense (the sky background), so the same physical noise reads as a huge
    # rank spread and the histogram gets smoothed into uselessness.
    #
    # Getting this wrong is visible two ways at once. Too much smoothing and
    # the local rank goes flat, so the filter degenerates into a plain gamma
    # stretch. Too little and the K-bin cumulative histogram stays a staircase,
    # which prints as concentric contour rings wherever the image has a smooth
    # radial gradient -- exactly the ring artifact this pipeline has been
    # chasing. The correct sigma suppresses the staircase by as much as the
    # real noise justifies and no more.
    # sigma_A is one number in value units, but the rank map's slope varies by
    # orders of magnitude across the frame -- it is steep in the sky, where
    # most pixels live, and shallow in the inner corona. So the smoothing width
    # is computed PER LEVEL by carrying sigma_A through the map at that level,
    # and applied as a single (K x K) mixing matrix over the level axis.
    if noise_sigma is None:
        d = np.abs(np.diff(x_lin[m][:200000]))
        s_a = 1.4826 * float(np.median(d)) / 1.4142 if d.size else 0.01
    else:
        s_a = float(noise_sigma)
    lv_val = np.interp(edges.astype(np.float64), qs / 100.0, knots)   # rank -> value
    lin_lv = (lv_val - lo) / max(hi - lo, 1e-9)
    r_hi = np.interp(np.clip(lin_lv + s_a, 0, 1) * (hi - lo) + lo, knots, qs / 100.0)
    r_lo = np.interp(np.clip(lin_lv - s_a, 0, 1) * (hi - lo) + lo, knots, qs / 100.0)
    sig_lv = np.clip(2.0 * np.abs(r_hi - r_lo) / 2.0, 0.5 * dl, 0.20)  # 2 sigma_A
    Mx = np.exp(-0.5 * ((edges[None, :] - edges[:, None]) / sig_lv[:, None]) ** 2)
    Mx /= np.maximum(Mx.sum(axis=1, keepdims=True), 1e-12)
    Hst = np.tensordot(Mx.astype(np.float32), Hst, axes=(1, 0))
    # --- cumulative histogram, and the rank of each pixel in it ---
    C = np.cumsum(Hst, axis=0)
    del Hst
    H_, W_ = x.shape
    yy = (np.arange(H_, dtype=np.float32) / q)[:, None] * np.ones((1, W_), np.float32)
    xx = (np.arange(W_, dtype=np.float32) / q)[None, :] * np.ones((H_, 1), np.float32)

    def _at(level):
        """Cumulative histogram at a per-pixel level, trilinear over the
        (level, y, x) grid."""
        b = np.clip(level / dl, 0, K - 1)
        return ndimage.map_coordinates(
            C, [b.ravel(), yy.ravel(), xx.ravel()], order=1,
            mode="nearest").reshape(x.shape)

    eps = float(eps_frac)
    c_mid = _at(x)
    if eps > 0:
        # --- variable neighbourhood: rank only among similar-valued pixels ---
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
    # eq. 2: the gamma-stretched original carries the large-scale brightness,
    # the equalized term carries the structure. T_gamma is applied to the
    # LINEAR-in-A version, not the rank-mapped one, so the base image keeps its
    # natural falloff instead of being flattened twice.
    if not combine:
        return np.where(m, E, 0.5).astype(np.float32)
    # eq. 2, for callers that want the paper's standalone display image.
    # T_gamma is applied to the LINEAR-in-A version, not the rank-mapped one,
    # so the base image keeps its natural falloff instead of being flattened
    # twice.
    out = (1.0 - w) * np.power(x_lin, 1.0 / max(gamma, 1e-3)) + w * E
    return np.where(m, out, 0.5).astype(np.float32)
