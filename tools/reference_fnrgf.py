"""A faithful port of Druckmuller's own FNRGF program, to use as an oracle.

Transcribed from the Delphi source shipped in `FNRGFsoftware.zip`
(`ImgProc.pas`, `Settings.ini`), not from the paper. The point is to have the
reference algorithm as something that can be RUN on the same data as ours,
rather than my reading of a description of it.

The reference, step by step:

1.  Statistics are gathered per INTEGER RADIUS r, over the cartesian pixels
    whose rounded radius is r, grouped into `SegmentCount` azimuthal segments
    (`ComputeSumX1SumX2N`). Default SegmentCount = 100.

2.  Per segment: the mean, and a Bessel-corrected standard deviation with an
    additive noise variance folded in before the square root
    (`FourierOnRadius`):

        dev = sqrt( SumX2/(N-1) - SumX1^2/((N-1)N) + Noise_AddVar )

3.  A trigonometric polynomial is fitted to those SEGMENT values -- not to the
    raw pixels -- by the discrete orthogonality sum, at phi = (seg+0.5)*2pi/NS
    (`ComputeFourierCoefs`). Two series: one for the means, one for the
    deviations. Order = (SegmentCount-1) div 2, capped at 80: 49 by default.

4.  Both series are evaluated WITH PER-HARMONIC ATTENUATION and the two series
    have DIFFERENT attenuation profiles (`PointAttenFourier`):

        F(phi) = A_0 a_0 / 2 + SUM_k A_k ( a_k cos k.phi + b_k sin k.phi )

5.  mask = (I - F_ave) / F_dev                        (`ApplyAttenFourier`)

6.  out = norm(I) + MixRatio * norm(mask), MixRatio = 2  (`MixInRatio`)

`Noise_AddVar` is estimated from the data as the median of the per-segment
variances over the outermost few radii (`EstimateAdditiveNoiseRing`).

Two things the reference does NOT do, both of which we do: it does not smooth
the radial profile of either series, and it does not fit robustly -- there is
no sigma-clipping or Huber weighting anywhere in it. Its defence against
outliers is the attenuation and the additive noise variance.

Not ported: the polar-coordinate bookkeeping (we index cartesian directly), the
incremental recompute logic in `ComputeFNRGF`, and the gamma/display code.
"""
import numpy as np


def reference_fnrgf(img, cy, cx, sun_r, init_r=1.0, fin_r=3.0,
                    segment_count=100, order=None,
                    atte_ave_step=0.03, atte_dev_step=0.04,
                    noise_radii=10, mix_ratio=2.0, return_parts=False):
    """`img` is the image the reference would be given (it works on pixel
    values, so hand it the same thing our FNRGF gets: log luminance).

    atte_*_step: the linear decrement per harmonic. The thesis's stated optimum
    is A_k = (1, 0.97, 0.94, ...) and C_k = (1, 0.96, 0.92, ...), i.e. 0.03 and
    0.04, with omega = 50. The shipped Settings.ini is gentler (0.002/0.02).
    """
    img = np.asarray(img, np.float64)
    H, W = img.shape
    if order is None:
        order = min((segment_count - 1) // 2, 80)
    yy = np.arange(H, dtype=np.float64)[:, None] - cy
    xx = np.arange(W, dtype=np.float64)[None, :] - cx
    rad = np.sqrt(yy * yy + xx * xx)
    phi = np.arctan2(yy + 0 * xx, xx + 0 * yy) % (2 * np.pi)
    r_i = np.rint(rad).astype(np.int32)
    seg = np.minimum((phi / (2 * np.pi) * segment_count).astype(np.int32),
                     segment_count - 1)

    r0, r1 = int(round(init_r * sun_r)), int(round(fin_r * sun_r))
    r1 = min(r1, int(min(cy, cx, H - cy, W - cx)) - 2)
    nr = r1 - r0 + 1
    if nr < 8:
        raise ValueError("annulus too thin")

    # --- per (radius, segment) sums, in one pass -------------------------
    inb = (r_i >= r0) & (r_i <= r1)
    flat = ((r_i - r0) * segment_count + seg)[inb]
    v = img[inb]
    n_rs = np.bincount(flat, minlength=nr * segment_count).astype(np.float64)
    s1 = np.bincount(flat, weights=v, minlength=nr * segment_count)
    s2 = np.bincount(flat, weights=v * v, minlength=nr * segment_count)
    n_rs = n_rs.reshape(nr, segment_count)
    s1 = s1.reshape(nr, segment_count)
    s2 = s2.reshape(nr, segment_count)
    safe = np.maximum(n_rs, 2.0)
    ave = s1 / np.maximum(n_rs, 1.0)
    var = s2 / (safe - 1.0) - (s1 * s1) / ((safe - 1.0) * safe)
    var = np.maximum(var, 0.0)

    # --- Noise_AddVar: median of the segment variances on the outer rings -
    tail = var[max(nr - noise_radii, 0):]
    noise_add_var = float(np.median(tail[np.isfinite(tail)])) if tail.size else 0.0
    dev = np.sqrt(var + noise_add_var)

    # --- trigonometric polynomial over the SEGMENTS ----------------------
    ph = (np.arange(segment_count) + 0.5) * (2 * np.pi / segment_count)
    kk = np.arange(1, order + 1)
    C = np.cos(kk[:, None] * ph[None, :])          # order x NS
    S = np.sin(kk[:, None] * ph[None, :])
    nc = 2.0 / segment_count
    A = {}
    for name, arr in (("ave", ave), ("dev", dev)):
        a0 = arr.sum(axis=1) * nc
        ak = (arr @ C.T) * nc                       # nr x order
        bk = (arr @ S.T) * nc
        A[name] = (a0, ak, bk)

    att_a = np.clip(1.0 - atte_ave_step * kk, 0.0, 1.0)
    att_c = np.clip(1.0 - atte_dev_step * kk, 0.0, 1.0)

    # --- evaluate both series at every pixel in the annulus ---------------
    ri = (r_i[inb] - r0)
    pv = phi[inb]
    cos_kp = np.cos(np.outer(pv, kk))               # npix x order
    sin_kp = np.sin(np.outer(pv, kk))
    out = {}
    for name, att in (("ave", att_a), ("dev", att_c)):
        a0, ak, bk = A[name]
        y = 0.5 * a0[ri]
        y += np.einsum("ij,ij->i", cos_kp, ak[ri] * att[None, :])
        y += np.einsum("ij,ij->i", sin_kp, bk[ri] * att[None, :])
        out[name] = y
    del cos_kp, sin_kp

    mask = np.zeros((H, W), np.float64)
    mask[inb] = (v - out["ave"]) / np.where(np.abs(out["dev"]) < 1e-12,
                                            1e-12, out["dev"])
    if return_parts:
        mu = np.zeros((H, W)); sd = np.zeros((H, W))
        mu[inb] = out["ave"]; sd[inb] = out["dev"]
        return mask, mu, sd, noise_add_var, inb

    # --- MixInRatio -------------------------------------------------------
    def nrm(a, m):
        lo, hi = np.min(a[m]), np.max(a[m])
        return (a - lo) / max(hi - lo, 1e-12)
    enh = np.zeros((H, W), np.float64)
    enh[inb] = (nrm(img, inb) + mix_ratio * nrm(mask, inb))[inb]
    return enh
