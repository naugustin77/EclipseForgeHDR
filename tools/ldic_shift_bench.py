#!/usr/bin/env python3
"""Does the LDIC azimuthal fit measure the tiers, or their misregistration?

WHY THIS EXISTS. `_fit_azimuthal_affine` is fitted on ring samples taken from
`stacks_half`, and the comment above that block said the tiers there are
"already in the common frame". They are not: `stacks_half[s]` is intra-tier
aligned only, and `abs_shift[s]` is applied in the merge loop, after the fit.

So every tier is sampled on rings centred at the MID tier's centre. A tier
offset by (dy, dx) has its corona centred somewhere else, and sampling it on
the wrong rings reads its radial profile at r + d*cos(phi - phi0) -- a
cos(phi) modulation whose amplitude is d * |d ln I / d ln r| / r. That is a
DIPOLE in azimuth, which is exactly the shape k(phi) is fitted to represent,
and the azimuthal-mean-1 normalisation removes the monopole, not the dipole.

The test: build tiers that differ ONLY by exposure time and a known shift --
no diffuse light, no cloud, nothing for LDIC to legitimately find. The true
k(phi) is 1 at every azimuth. Anything the fit reports is error.

    python tools/ldic_shift_bench.py
"""
import os
import sys

import numpy as np
from scipy import ndimage

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from eclipseforgehdr.pipeline import (_fit_azimuthal_affine,  # noqa: E402
                                      _LDIC_SEGMENTS)


def corona(H, W, cy, cx, R, seed=3):
    """A corona with a steep limb gradient and real azimuthal structure."""
    yy = np.arange(H, dtype=np.float64)[:, None] - cy
    xx = np.arange(W, dtype=np.float64)[None, :] - cx
    r = np.hypot(yy, xx)
    th = np.arctan2(yy, xx)
    rng = np.random.default_rng(seed)
    # streamers: a few broad lobes plus fine structure, both radius-independent
    thr = 1.0 + 0.30 * np.sin(3 * th + 0.7) + 0.12 * np.sin(11 * th + 2.1)
    img = 3.0e6 * (R / np.maximum(r, R)) ** 6 * thr
    img[r < R] = 0.0
    return img


def sample(stack, ys, xs):
    return stack[ys, xs].astype(np.float64)


def run(shift_px, H=900, W=1300, R=110.0, secs=(1 / 500, 1 / 60, 1 / 8, 1.0)):
    """Return (k spread now, k spread with the shift applied) in percent."""
    cy, cx = H / 2.0, W / 2.0
    truth = corona(H, W, cy, cx, R)

    # per-tier absolute shift, as the alignment network would report it: the
    # mid tier is the reference, the others drift progressively
    ref = len(secs) // 2
    shifts = {s: ((i - ref) * shift_px * 0.8, (i - ref) * shift_px * 0.6)
              for i, s in enumerate(secs)}

    stacks, sats = {}, {}
    for s in secs:
        dy, dx = shifts[s]
        # the tier as the camera recorded it: the scene, moved
        a = ndimage.shift(truth, (dy, dx), order=1, mode="nearest") * s
        stacks[s] = a.astype(np.float32)
        sats[s] = np.zeros((H, W), bool)

    _rs = np.arange(1.00 * R, 2.5 * R, 0.02 * R)
    _NA = 12 * _LDIC_SEGMENTS
    _th = np.linspace(0, 2 * np.pi, _NA, endpoint=False)
    cal = {s: 1.0 for s in secs}
    order = sorted(secs, reverse=True)

    def fit(use_shift, bilinear=False):
        vals, wts = [], []
        for s in order:
            dy, dx = shifts[s] if use_shift else (0.0, 0.0)
            yf = cy + dy + _rs[:, None] * np.sin(_th)
            xf = cx + dx + _rs[:, None] * np.cos(_th)
            if bilinear:
                f = ndimage.map_coordinates(
                    stacks[s].astype(np.float64),
                    [np.clip(yf, 0, H - 1), np.clip(xf, 0, W - 1)],
                    order=1, mode="nearest")
            else:
                ys = np.clip(yf.astype(np.int32), 0, H - 1)
                xs = np.clip(xf.astype(np.int32), 0, W - 1)
                f = sample(stacks[s], ys, xs)
            vals.append(f / (s * cal[s]))
            wts.append(np.ones_like(f))
        return _fit_azimuthal_affine(vals, wts, order)

    def spread(ld):
        # what stats["ldic_k_spread_pct"] reports: the 10-90 percentile span
        return {("%g" % s): 100.0 * float(np.percentile(v[0], 90) -
                                          np.percentile(v[0], 10))
                for s, v in ld.items()}

    return spread(fit(False)), spread(fit(True)), spread(fit(True, True))


def main():
    print("A tier's k(phi) should be FLAT: these tiers differ only by exposure")
    print("time and a rigid shift. The 10-90 percentile span of k, in percent:")
    print()
    print("  tier shift     tier        as now      registered   registered")
    print("  (half-res px)              (unshifted)  (integer)    (bilinear)")
    for sp in (0.0, 2.0, 5.0, 11.0, 25.0):
        now, fixed, bil = run(sp)
        first = True
        for s in sorted(now, key=float):
            lbl = f"  {sp:>5.0f}        " if first else "               "
            print(f"{lbl}{s:>10}  {now[s]:>10.1f}   {fixed[s]:>10.1f}   {bil[s]:>10.1f}")
            first = False
        print()


if __name__ == "__main__":
    main()
