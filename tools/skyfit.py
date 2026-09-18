"""Prototype: measure the sky as an ADDITIVE per-channel term on a merged HDR.

WHY THIS EXISTS. "Remove sky" subtracts a corner median from every frame before
the merge -- 251 independently estimated constants on a 251-frame set, each
error landing flat on a whole frame, and each one taken before the merge weights
are computed, so the outer corona drops under the hat's noise floor and loses
tiers. Measured on one tester's set, that costs 26-34% more pixel-to-pixel noise
beyond 3 R.

The sky scales with exposure exactly as the corona does, so after each tier is
divided by its exposure ratio the sky arrives as the SAME constant from every
tier. One subtraction on the merged result is therefore the same operation, with
one estimate instead of 251, and it happens after the weighting. That is the
deep-sky order -- calibrate per frame, extract the background after integration.

THE PROBLEM THIS SCRIPT EXISTS TO MEASURE. There is no blank sky in an eclipse
frame: the corona reaches the corners. So the sky has to be separated from the
corona by its RADIAL BEHAVIOUR -- the corona keeps falling, the sky does not.
Fitting I(r) = A r^-k + C per channel does that in principle. In practice A, k
and C trade off against each other once the corona is faint, and an unconstrained
fit absorbs corona into C.

    fit range      C error vs ground truth      k (green)
    k free, 1.6 R        +40%                     5.58   <- unphysical
    k free, 2.0 R        +10%                     3.34
    k free, 2.4 R        -17%                     2.17
    k free, 2.8 R        -73%                     1.30
    k free, 3.2 R       -386%                     0.86   <- degenerate
    k in [2.0,3.5], 1.6 R  -16%                   3.50
    k in [2.0,3.5], 2.0 R  +10%                   3.34
    k in [2.0,3.5], 3.2 R  -22%                   2.00

Constrained, the level lands within about +/-20% wherever the fit starts, and the
COLOUR is recovered to a few percent (R/G 0.723 fitted against 0.734 measured,
B/G 1.159 against 1.260). Ground truth here is the difference between two merges
of the same 251 frames, one with the sky left in and one with it removed per
frame, after matching their scale on the inner corona -- itself imperfect,
because the per-frame corner median eats some corona too.

VERDICT: usable for the colour, good to roughly a fifth for the level. Not yet
good enough to subtract blind. Next: a two-component corona (K + F) instead of
one power law, or anchor the level on the longest frame's own corners and fit
only the shape.

Usage:  python3 tools/skyfit.py <folder>   (the folder you stacked)
"""
import json
import os
import sys

import numpy as np
from scipy.optimize import curve_fit

K_BOUNDS = (2.0, 3.5)      # physical coronal falloff; see the table above
FIT_FROM = 2.0             # R, inner edge of the fit
FIT_TO = 5.2               # R, outer edge


def radial_profile(a, r, lo, hi, nbin=40, minpix=500):
    """Median per channel in log-spaced annuli. A median, so a streamer or a
    prominence occupying part of one azimuth does not move it."""
    edges = np.geomspace(lo, hi, nbin + 1)
    centres = np.sqrt(edges[:-1] * edges[1:])
    out = np.full((nbin, 3), np.nan)
    for i in range(nbin):
        m = (r >= edges[i]) & (r < edges[i + 1])
        if int(m.sum()) >= minpix:
            out[i] = np.median(a[m].reshape(-1, 3), axis=0)
    return centres, out


def fit_sky(a, r, lo=FIT_FROM, hi=FIT_TO, k_bounds=K_BOUNDS):
    """Per-channel additive sky level from I(r) = A r^-k + C."""
    c, p = radial_profile(a, r, lo, hi)
    ok = np.isfinite(p).all(axis=1)
    if ok.sum() < 8:
        raise RuntimeError("too few usable annuli between %.1f and %.1f R" % (lo, hi))

    def model(x, A, k, C):
        return A * x ** (-k) + C

    sky = np.zeros(3)
    kk = np.zeros(3)
    for ch in range(3):
        po, _ = curve_fit(model, c[ok], p[ok, ch], p0=[3000.0, 3.0, 120.0],
                          maxfev=60000,
                          bounds=([0.0, k_bounds[0], -1e5],
                                  [1e9, k_bounds[1], 1e5]))
        sky[ch], kk[ch] = po[2], po[1]
    return sky, kk


def main(folder):
    wd = os.path.join(folder, ".eclipseforgehdr")
    geo = json.load(open(os.path.join(wd, "geometry.json")))
    cy, cx, R = geo["cy"] / 2.0, geo["cx"] / 2.0, geo["R"] / 2.0
    a = np.array(np.load(os.path.join(wd, "hdr_rgb.npy"), mmap_mode="r")[::2, ::2],
                 np.float32)
    h, w = a.shape[:2]
    yy = np.arange(h, dtype=np.float32)[:, None] - cy
    xx = np.arange(w, dtype=np.float32)[None, :] - cx
    r = np.hypot(yy, xx) / R

    sky, kk = fit_sky(a, r)
    print("sky, fitted on the merged HDR (%.1f-%.1f R, k constrained to %.1f-%.1f):"
          % (FIT_FROM, FIT_TO, K_BOUNDS[0], K_BOUNDS[1]))
    for ch, nm in enumerate("RGB"):
        print("   %s  %8.2f    corona falls as r^-%.2f" % (nm, sky[ch], kk[ch]))
    print("   sky colour  R/G %.3f  B/G %.3f" % (sky[0] / sky[1], sky[2] / sky[1]))

    # what subtracting it would do to the far-field colour
    print("\n   radius     before                 after")
    for lo in (2.0, 2.5, 3.0, 3.5, 4.0, 4.5):
        m = (r >= lo) & (r < lo + 0.25)
        v = np.median(a[m].reshape(-1, 3), axis=0)
        u = v - sky
        print("    %.1f R   R/G %5.2f B/G %5.2f     R/G %6.2f B/G %6.2f%s"
              % (lo, v[0] / v[1], v[2] / v[1],
                 u[0] / u[1] if u[1] > 0 else float("nan"),
                 u[2] / u[1] if u[1] > 0 else float("nan"),
                 "   <- over-subtracted" if u.min() <= 0 else ""))


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(os.path.expanduser(sys.argv[1]))
