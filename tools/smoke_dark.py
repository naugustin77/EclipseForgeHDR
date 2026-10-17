#!/usr/bin/env python3
"""End to end: do bias/ and darks/ beside the lights actually get removed?

THE FIXTURE. A grey corona plus three planted defects, every one of them
additive and fixed in place, which is what makes them invisible to stacking:

    offset      a constant electronic offset under every pixel
    rate        per-pixel dark current, so it grows with the tier's exposure
    hot pixels  a few hundred photosites with a much larger rate

`bias/` gets frames with the offset and read noise only. `darks/` gets frames
at one exposure with the offset, the rate and the hot pixels. The lights get
all three, scaled by their own exposure.

WHAT WOULD BE MISSED WITHOUT THIS. The subtraction can be arithmetically right
and still land in the wrong place -- after the flat instead of before it, or
without the per-tier exposure scaling, or on the wrong Bayer phase. A dark
scaled wrongly is worst on the LONGEST tier, which is the outer field, which is
where the corona work lives. So the test measures the outer field.

    python tools/smoke_dark.py [-v]
"""
import json
import os
import shutil
import sys
import tempfile

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, "..", "src"))
sys.path.insert(0, _HERE)

from smoke_pipeline import _P as Prog                    # noqa: E402
from smoke_channel_floors import write_cube              # noqa: E402

_H, _W, _R = 420, 560, 46.0
#: The ladder ENDS at the dark frames' exposure. A dark is only worth anything
#: on the tiers whose exposure is comparable to it -- the first fixture here
#: topped out at 0.4 s against 1.6 s darks, so the dark signal being corrected
#: was 0.1 ADU under 12 ADU of noise and the correction could only add noise.
#: The module said so ("removes 0.691, injects 0.712") and the measured
#: scatter agreed. That was the fixture being wrong, not the code.
_DARK_T = 1.6             # seconds, the dark frames' exposure
_N_TIERS, _STEP = 8, 24.0 ** (1 / 7.0)
_T0 = _DARK_T / 24.0
_OFFSET = 96.0            # ADU, constant electronic offset
_READ = 3.0               # ADU, read noise


def _defects(seed=5):
    """Per-pixel dark current in ADU/s: a smooth floor plus hot photosites."""
    rng = np.random.default_rng(seed)
    rate = rng.gamma(2.0, 1.5, size=(_H, _W)).astype(np.float32)
    hot = rng.random((_H, _W)) < 0.0008
    rate[hot] += rng.uniform(60.0, 400.0, size=int(hot.sum())).astype(np.float32)
    return rate, int(hot.sum())


def _corona():
    cy, cx = _H / 2.0, _W / 2.0
    yy = np.arange(_H, dtype=np.float64)[:, None] - cy
    xx = np.arange(_W, dtype=np.float64)[None, :] - cx
    r = np.hypot(yy, xx)
    th = np.arctan2(yy, xx)
    thr = 1.0 + 0.30 * np.sin(3 * th + 0.7) + 0.10 * np.sin(11 * th)
    base = 6.0e4 * (_R / np.maximum(r, _R)) ** 4 * thr
    base[r < _R] = 0.0
    return base, r


def build(folder, rate, with_calib, rng, offset=_OFFSET):
    base, _r = _corona()
    gain = np.array([0.70, 0.88, 1.00])[:, None, None]
    os.makedirs(folder, exist_ok=True)
    n = 0
    # Shot noise follows the COLLECTED CHARGE -- corona photons and dark
    # current -- and not the electronic offset, which is added downstream of
    # the well and carries only read noise. The first version of this fixture
    # put sqrt(96) of shot noise on the offset in every frame, inflating the
    # master's measured noise by 4x and making the correction look worthless.
    for i in range(_N_TIERS):
        t = _T0 * _STEP ** i
        for _ in range(2):
            charge = base[None, :, :] * gain * t + rate[None, :, :] * t
            a = charge + offset
            cube = a + rng.normal(0, np.hypot(np.sqrt(np.maximum(charge, 0)),
                                              _READ))
            write_cube(os.path.join(folder, "f%02d.fits" % n), cube, t,
                       "2026-08-12T18:%02d:%02d.0" % (n // 60, n % 60))
            n += 1
    if not with_calib:
        return
    bd = os.path.join(folder, "bias")
    os.makedirs(bd, exist_ok=True)
    for k in range(24):
        a = np.full((3, _H, _W), offset)
        write_cube(os.path.join(bd, "b%02d.fits" % k),
                   a + rng.normal(0, _READ, a.shape), 1 / 8000.0,
                   "2026-08-12T17:%02d:00.0" % (k % 60))
    dd = os.path.join(folder, "darks")
    os.makedirs(dd, exist_ok=True)
    for k in range(24):
        charge = np.broadcast_to(rate[None, :, :] * _DARK_T, (3, _H, _W))
        a = charge + offset
        write_cube(os.path.join(dd, "d%02d.fits" % k),
                   a + rng.normal(0, np.hypot(np.sqrt(np.maximum(charge, 0)),
                                              _READ), a.shape),
                   _DARK_T, "2026-08-12T17:%02d:30.0" % (k % 60))


#: Luminance weights times the per-channel gains. The FITS path has no white
#: balance and an identity camera matrix, so the merged luminance of a planted
#: corona is exactly this constant times `base`.
_LUMW = 0.2126 * 0.70 + 0.7152 * 0.88 + 0.0722 * 1.00


def _truth(lo=3.0, hi=5.0):
    """What the merged outer field must read if everything is removed.

    The merge divides by s*cal, so it returns counts per second: the corona
    rate and nothing else. An uncorrected dark current adds its OWN rate to
    that, at every tier alike, which is why it survives the merge as a bias
    rather than averaging out.
    """
    base, r = _corona()
    m = (r / _R > lo) & (r / _R < hi)
    return _LUMW * float(np.median(base[m]))


def _outer(rgb, cy, cx, R, lo=3.0, hi=5.0):
    """(median, high-frequency scatter) in the annulus.

    The plain scatter of the annulus is dominated by the corona's own radial
    falloff and its azimuthal structure -- real signal, tens of times larger
    than anything a dark frame removes. Subtracting a 5 px median leaves the
    per-pixel part, which is where fixed-pattern noise and hot pixels live and
    where a dark can be seen at all.
    """
    from scipy import ndimage
    lum = (0.2126 * rgb[..., 0] + 0.7152 * rgb[..., 1] + 0.0722 * rgb[..., 2])
    hf = lum - ndimage.median_filter(lum, size=5, mode="nearest")
    yy = np.arange(lum.shape[0])[:, None] - cy
    xx = np.arange(lum.shape[1])[None, :] - cx
    rr = np.hypot(yy, xx) / R
    m = (rr > lo) & (rr < hi) & np.isfinite(lum)
    v, h = lum[m], hf[m]
    return (float(np.median(v)),
            float(1.4826 * np.median(np.abs(h - np.median(h)))))


def run_once(tmp, with_calib, verbose, clean=False):
    from eclipseforgehdr.pipeline import run
    rate, nhot = _defects()
    if clean:
        # No offset, no dark current: what the pipeline returns on data that
        # never had the defects at all. Without this the residual error of the
        # corrected run has nothing to be judged against -- some of it belongs
        # to the pipeline and the fixture's own truth model, not to the
        # correction, and there is no way to tell which from two runs.
        rate = np.zeros_like(rate)
    folder = os.path.join(tmp, "clean" if clean
                          else ("on" if with_calib else "off"))
    build(folder, rate, with_calib, np.random.default_rng(3),
          offset=0.0 if clean else _OFFSET)
    run(folder, Prog(verbose), crop_pc=320)
    wd = os.path.join(folder, ".eclipseforgehdr")
    st = json.load(open(os.path.join(wd, "report.json")))
    geo = json.load(open(os.path.join(wd, "geometry.json")))
    rgb = np.load(os.path.join(wd, "hdr_rgb.npy"))
    return st, _outer(rgb, geo["cy"], geo["cx"], geo["R"]), rate, nhot


def main(verbose=False):
    tmp = tempfile.mkdtemp(prefix="efhdr_dark_")
    try:
        st_off, off, rate, nhot = run_once(tmp, False, verbose)
        st_on, on, _, _ = run_once(tmp, True, verbose)
        _st_cl, clean, _, _ = run_once(tmp, False, verbose, clean=True)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    cal = st_on.get("calib") or {}
    bi = cal.get("bias") or {}
    di = cal.get("dark") or {}
    print("planted    offset %.1f ADU, dark current mean %.4f median %.4f "
          "ADU/s, %d hot photosite(s)"
          % (_OFFSET, float(np.mean(rate)), float(np.median(rate)), nhot))
    print("found      bias dir %s   darks dir %s"
          % (os.path.basename(cal.get("bias_dir") or "-"),
             os.path.basename(cal.get("dark_dir") or "-")))
    print("           bias  used %s frame(s), level %s, noise %s"
          % (bi.get("n_used"), _f(bi.get("level")), _f(bi.get("noise"))))
    print("           dark  used %s frame(s) at %ss, rate mean %s median %s "
          "ADU/s" % (di.get("n_used"), di.get("seconds"),
                     _f(di.get("rate_mean"), 4), _f(di.get("rate_median"), 4)))
    print("applied    bias %s   dark %s"
          % (cal.get("bias_applied"), cal.get("dark_applied")))
    at = di.get("at_longest_tier") or {}
    if at:
        print("           at the longest tier (%.3gs): removes %s, injects %s"
              % (at.get("seconds", 0), _f(at.get("removes"), 3),
                 _f(at.get("injects"), 3)))
    tru = _truth()
    print()
    print("merged outer field, 3-5 R     median   error vs truth   per-px "
          "scatter")
    print("  truth               %11.4g" % tru)
    print("  no defects planted  %11.4g %12.1f%% %14.4g"
          % (clean[0], 100 * (clean[0] / tru - 1), clean[1]))
    print("  calibration OFF     %11.4g %12.1f%% %14.4g"
          % (off[0], 100 * (off[0] / tru - 1), off[1]))
    print("  calibration ON      %11.4g %12.1f%% %14.4g"
          % (on[0], 100 * (on[0] / tru - 1), on[1]))
    print("  ON vs the defect-free run: %+.1f%%"
          % (100 * (on[0] / clean[0] - 1)))

    bad = []
    if not cal.get("bias_applied"):
        bad.append("the bias was not applied")
    if not cal.get("dark_applied"):
        bad.append("the dark was not applied")
    # The MEAN, not the median -- see dark.py's note on why a median across
    # pixels of (skewed truth + symmetric noise) is biased towards the mean.
    true_rate = float(np.mean(rate))
    got = di.get("rate_mean")
    if got is None or abs(got / max(true_rate, 1e-9) - 1.0) > 0.10:
        bad.append("mean dark current came back %s, planted %.4f"
                   % (_f(got, 4), true_rate))
    # THE BIAS, NOT THE SCATTER, IS THE CLAIM. Dark current is added to every
    # tier in proportion to its exposure, so after the merge divides by that
    # exposure it survives as a constant added to the radiance -- identical in
    # every frame, so no amount of stacking touches it. It shows up as the
    # merged outer field reading HIGH against truth, and removing it is the
    # whole point. Per-pixel scatter is reported too but is not asserted on:
    # in this fixture the corona's own photon noise is several times the fixed
    # pattern, so the improvement there is real but small.
    # Judged against the DEFECT-FREE RUN, not against the analytic truth. The
    # pipeline's own pedestal fit, the fake-CFA demosaic and the merge all put
    # a small offset between the analytic corona and what a clean run returns;
    # that offset is not the correction's to answer for, and including it would
    # either hide a real residual or invent one.
    e_off = abs(off[0] / clean[0] - 1)
    e_on = abs(on[0] / clean[0] - 1)
    if not (e_off > 0.10):
        bad.append("the planted defects did not bias the merged outer field "
                   "to begin with (%.1f%%) — the test proves nothing"
                   % (100 * e_off))
    if not (e_on < 0.35 * e_off):
        bad.append("the outer-field bias was not removed: %.1f%% -> %.1f%%"
                   % (100 * e_off, 100 * e_on))
    print()
    if bad:
        for b in bad:
            print("FAIL:", b)
        return 1
    print("OK — bias/ and darks/ were found beside the lights, the mean dark "
          "current came back within %.2f%% of what was planted, and the "
          "merged outer field went from %.1f%% off the defect-free run to "
          "%.1f%%." % (100 * abs(got / true_rate - 1), 100 * e_off, 100 * e_on))
    return 0


def _f(v, n=3):
    return ("%.*f" % (n, v)) if isinstance(v, (int, float)) else str(v)


if __name__ == "__main__":
    sys.exit(main("-v" in sys.argv))
