#!/usr/bin/env python3
"""Does the bias measure the read noise, and the dark count what is warm?

TWO CLAIMS, MEASURED SEPARATELY.

1. READ NOISE. `dark._stack` estimates per-frame read noise from the two
   half-stack means: sd(meanA - meanB) = sigma * sqrt(1/nA + 1/nB). Plant a
   known sigma, per channel, and see it come back. This is the number that
   pins the intercept of the hot-pixel noise model, so an error here moves the
   detection threshold on every frame of the run.

2. THE WARM-PIXEL CENSUS. The master dark's rate map is thresholded against
   each photosite's same-colour neighbours to count what is running hot. Plant
   warm photosites at known positions and check two things: how many come
   back, and -- the half that matters -- how many are flagged on a sensor with
   NONE.

   This census is NOT a repair map, and the bench does not compare it against
   the light-frame defect path, because they do not measure the same thing.
   At 1/500 s the worst photosite on the reference sensor (122.9 ADU/s)
   contributes 0.25 ADU: dark current cannot make a detectable outlier in a
   2 ms frame, so what the light tier finds is offset and stuck-response
   defects, not warm pixels. And a warm pixel needs no repair -- its own rate
   is subtracted per photosite, which is the correct fix. A real sensor has no
   warm population to separate either: the rate distribution is a continuum,
   p50 0.234 to p99.9 9.77 to a max of 122.9 ADU/s.

    python tools/smoke_defects.py [-v]
"""
import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, "..", "src"))
sys.path.insert(0, _HERE)

H, W = 400, 600
_READ = (3.1, 2.4, 2.4, 4.6)       # per Bayer offset: R, G, G, B
_GAIN = 0.42                        # ADU per electron-ish, for the light frames
_N_DEFECTS = 260
_RATE = 0.23                        # ADU/s, roughly the S1R II at 30 C


def _rn_map():
    m = np.empty((H, W), np.float32)
    for i, (oy, ox) in enumerate(((0, 0), (0, 1), (1, 0), (1, 1))):
        m[oy::2, ox::2] = _READ[i]
    return m


def _defects(rng):
    """Photosite indices and their excess dark current, in ADU/s."""
    ys = rng.integers(2, H - 2, _N_DEFECTS)
    xs = rng.integers(2, W - 2, _N_DEFECTS)
    keep = {}
    for y, x in zip(ys, xs):
        keep[(int(y), int(x))] = 0.0
    ks = sorted(keep)
    # a realistic spread: mostly mild, a few extreme
    exc = 10.0 ** rng.uniform(0.4, 2.2, len(ks))
    m = np.zeros((H, W), np.float32)
    for (y, x), e in zip(ks, exc):
        m[y, x] = e
    return m


def make_bias(rng, n=50):
    rn = _rn_map()
    return [rng.normal(0.0, rn).astype(np.float32) for _ in range(n)]


def make_darks(rng, defmap, seconds=1.6, n=20):
    rn = _rn_map()
    base = (_RATE + defmap) * seconds
    return [(base + rng.normal(0.0, rn)).astype(np.float32) for _ in range(n)]


def _score(found, truth):
    t = truth > 0
    tp = int((found & t).sum())
    fp = int((found & ~t).sum())
    return tp, fp, 100.0 * tp / max(int(t.sum()), 1)


def main(verbose=False):
    from eclipseforgehdr import dark as D
    from scipy import ndimage
    rng = np.random.default_rng(7)
    defmap = _defects(rng)
    n_def = int((defmap > 0).sum())

    # ---- 1. read noise ---------------------------------------------------
    bias = make_bias(rng)
    nA = len(bias[0::2]); nB = len(bias[1::2])
    mA = np.mean(bias[0::2], axis=0); mB = np.mean(bias[1::2], axis=0)
    f = np.sqrt(1.0 / nA + 1.0 / nB)
    got = [float(D._robust_sd((mA - mB)[oy::2, ox::2]) / f)
           for oy in (0, 1) for ox in (0, 1)]
    print("READ NOISE from %d bias frames" % len(bias))
    print("   offset     planted   measured    error")
    bad = []
    for i, nm in enumerate(("R ", "G1", "G2", "B ")):
        e = 100 * (got[i] / _READ[i] - 1)
        print("     %s      %6.2f     %6.2f  %+7.1f%%" % (nm, _READ[i], got[i], e))
        if abs(e) > 5.0:
            bad.append("read noise %s off by %.1f%%" % (nm, e))

    # ---- 2. defects, from the darks -------------------------------------
    darks = make_darks(rng, defmap)
    t = 1.6
    dA = np.mean(darks[0::2], axis=0); dB = np.mean(darks[1::2], axis=0)
    bA = np.mean(bias[0::2], axis=0); bB = np.mean(bias[1::2], axis=0)
    master = np.mean(darks, axis=0) - np.mean(bias, axis=0)
    # both masters' noise, as dark.build computes `injects` -- the rate map is
    # their difference, so the bias's share is in it too
    def _sig4(xA, xB, yA, yB):
        """Per-Bayer-offset noise of (mean(x) - mean(y)), as dark.build does."""
        return [float(np.hypot(D._robust_sd(((xA - xB) * 0.5)[oy::2, ox::2]),
                               D._robust_sd(((yA - yB) * 0.5)[oy::2, ox::2])))
                for oy in (0, 1) for ox in (0, 1)]
    sig4 = [v / t for v in _sig4(dA, dB, bA, bB)]
    rate = (master / t).astype(np.float32)
    sig = float(np.mean(sig4))
    hot_d = np.zeros((H, W), bool)
    for i, (oy, ox) in enumerate(((0, 0), (0, 1), (1, 0), (1, 1))):
        s = rate[oy::2, ox::2]
        m = ndimage.median_filter(s, size=3, mode="nearest")
        hot_d[oy::2, ox::2] = (s - m) > D.DEFECT_K * sig4[i]
    tp_d, fp_d, rec_d = _score(hot_d, defmap)

    # ---- 3. the threshold, on a sensor with NO defects at all -----------
    # The recall number above is worthless without this. A threshold low
    # enough to find every planted defect will also flag good photosites, and
    # repairing a good photosite is a silent edit to real data -- worse than
    # missing a defect, because nothing downstream can tell it happened.
    clean = make_darks(rng, np.zeros((H, W), np.float32))
    cA = np.mean(clean[0::2], axis=0); cB = np.mean(clean[1::2], axis=0)
    cm = np.mean(clean, axis=0) - np.mean(bias, axis=0)
    csig4 = [v / t for v in _sig4(cA, cB, bA, bB)]
    crate = (cm / t).astype(np.float32)
    fp_clean = 0
    worst = 0.0
    for i, (oy, ox) in enumerate(((0, 0), (0, 1), (1, 0), (1, 1))):
        s = crate[oy::2, ox::2]
        m = ndimage.median_filter(s, size=3, mode="nearest")
        fp_clean += int(((s - m) > D.DEFECT_K * csig4[i]).sum())
        worst = max(worst, float((s - m).max() / csig4[i]))

    print()
    print("WARM PIXELS — %d planted photosites, %d dark frames" % (n_def, len(darks)))
    print("   recall                     %5d found   %.1f%%" % (tp_d, rec_d))
    print("   false positives            %5d" % fp_d)
    print("   sigma of the rate map      %.4f-%.4f ADU/s per channel"
          % (min(sig4), max(sig4)))
    print("   threshold                  %.3f-%.3f ADU/s"
          % (D.DEFECT_K * min(sig4), D.DEFECT_K * max(sig4)))
    print("   faintest planted defect    %.2f ADU/s" % float(defmap[defmap > 0].min()))
    print()
    print("   on a DEFECT-FREE sensor:   %d flagged; largest excursion %.2f sigma"
          % (fp_clean, worst))

    # The light-frame path is deliberately NOT benchmarked against this. It
    # would need a model of what a defect looks like at 1/500 s, and the only
    # honest source for that is a real sensor -- a defect modelled as pure
    # excess dark current contributes 0.2 ADU at that shutter speed and would
    # make the light path look useless for a reason that is an artefact of the
    # fixture. The real comparison is the two counts printed by a real run.
    if rec_d < 60.0:
        bad.append("the dark found only %.0f%% of the planted defects" % rec_d)
    if fp_d > 0.0002 * H * W:
        bad.append("%d false positives against the planted truth" % fp_d)
    if fp_clean > 0:
        bad.append("%d photosites flagged on a sensor that has no defects — "
                   "the threshold is too low and good data is being edited"
                   % fp_clean)
    print()
    if bad:
        for b in bad:
            print("FAIL:", b)
        return 1
    print("OK — read noise recovered per channel to better than 5%%; the master "
          "dark finds %.0f%% of the planted defects, and flags nothing at all "
          "on a defect-free sensor (largest excursion %.2f of the %.0f sigma "
          "threshold)." % (rec_d, worst, D.DEFECT_K))
    return 0





# ---------------------------------------------------------------------------
# 3. DOES PINNING THE INTERCEPT HELP?
#
# On the 600 mm run, feeding the measured read noise into the light-frame
# detector took the count from 3675 to 1848. Fewer is not the same as better,
# and the run itself cannot say which is right -- so plant the truth.
#
# The defects here are OFFSET defects, present in every frame regardless of
# exposure, because that is what the light path actually detects: at 1/500 s
# dark current contributes a quarter of an ADU and cannot make an outlier.
# ---------------------------------------------------------------------------

#: The intercept trial runs on its own, larger grid. On the 400x600 one the
#: corona's CURVATURE across a 3x3 same-colour neighbourhood is ~68 ADU, which
#: swamps every planted defect and floods both paths with limb false positives
#: (523 each, 10.7% recall) -- a fixture artefact, not a finding. What matters
#: is curvature per pixel, not radius in pixels: a linear ramp leaves the
#: median of a symmetric neighbourhood unchanged and contributes nothing. On
#: the real 600 mm frame R is 618 px and that curvature is about 1.3 ADU
#: against 2.7 ADU of read noise, so the grid below is sized to match.
_IH, _IW, _IR = 700, 900, 300.0


def _irn_map():
    m = np.empty((_IH, _IW), np.float32)
    for i, (oy, ox) in enumerate(((0, 0), (0, 1), (1, 0), (1, 1))):
        m[oy::2, ox::2] = _READ[i]
    return m


def _light_tier(rng, offs, n=4, peak=3.0e3):
    """A 1/500 s tier: bright inner corona to near-black outer field.

    A FLAT SKY PROVED NOTHING. With one level everywhere the noise model has
    no signal range, so the fit measures the intercept directly and cannot get
    it wrong -- that version came out 190 against 189. Pinning the intercept
    only matters where gain and read noise trade off, which needs `med` to
    actually span a range.
    """
    rn = _irn_map()
    cy, cx = _IH / 2.0, _IW / 2.0
    yy = np.arange(_IH, dtype=np.float32)[:, None] - cy
    xx = np.arange(_IW, dtype=np.float32)[None, :] - cx
    r = np.hypot(yy, xx)
    corona = peak * (_IR / np.maximum(r, _IR)) ** 3
    corona[r < _IR] = 0.0
    corona += 6.0
    out = []
    for _ in range(n):
        a = corona + offs
        shot = rng.normal(0.0, np.sqrt(np.maximum(a, 0) * _GAIN)) / np.sqrt(_GAIN)
        out.append((a + shot + rng.normal(0.0, rn)).astype(np.float32))
    return out


def intercept_trial(verbose=False):
    from eclipseforgehdr.raw import _outlier_flags
    rng = np.random.default_rng(19)
    # offsets spanning "just detectable" to "obvious", as a real sensor does
    ys = rng.integers(2, _IH - 2, 300); xs = rng.integers(2, _IW - 2, 300)
    offs = np.zeros((_IH, _IW), np.float32)
    mag = 10.0 ** rng.uniform(0.9, 2.3, ys.size)        # 8 .. 200 ADU
    for y, x, m in zip(ys, xs, mag):
        offs[int(y), int(x)] = m
    truth = offs > 0
    n_true = int(truth.sum())
    lights = _light_tier(rng, offs)

    print()
    print("INTERCEPT — %d planted offset defects on a 1/500 s tier" % n_true)
    print("   read noise                 found   recall   false +ve")
    rows = {}
    for label, rn in (("fitted from the tier", None),
                      ("measured from bias  ", list(_READ))):
        cnt = np.zeros((_IH, _IW), np.uint8)
        for b in lights:
            cnt += _outlier_flags(b, 6.0, rn)
        got = cnt >= max(2, int(np.ceil(0.6 * len(lights))))
        tp, fp, rec = _score(got, offs)
        rows[label] = (tp, fp, rec)
        print("   %s     %5d %7.1f%% %10d" % (label, tp, rec, fp))
    a = rows["fitted from the tier"]
    b = rows["measured from bias  "]
    print()
    print("   change: %+d found, %+d false positives" % (b[0] - a[0], b[1] - a[1]))
    bad = []
    if b[1] > a[1]:
        bad.append("pinning the intercept ADDED %d false positives" % (b[1] - a[1]))
    if b[0] < 0.9 * a[0]:
        bad.append("pinning the intercept lost %d real defects (%.0f%% of what "
                   "the fit found)" % (a[0] - b[0], 100 * (1 - b[0] / max(a[0], 1))))
    for x in bad:
        print("   FAIL:", x)
    return 1 if bad else 0


if __name__ == "__main__":
    _v = "-v" in sys.argv
    _rc = main(_v)
    sys.exit(intercept_trial(_v) or _rc)
