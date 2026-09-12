#!/usr/bin/env python3
"""Does per-frame grouping rescue a bracket whose headers lie about grouping?

THE ONE THING GROUPING TRUSTS. Everything downstream of the tiers measures
past the exposure metadata -- `cal[s]` is fitted from the data, and
smoke_blind_exposure shows the merged levels survive a wrong header ladder
with one constant and 0.6% of shape. But two frames land in the SAME tier, and
are averaged as though they were the same measurement, purely because their
headers agree. Nothing measures that.

THE FIXTURE. A clean bracket, except that one pair of frames CLAIMS the same
shutter speed while one of the two really got 1.45x more light. Averaging them
produces a tier that is neither frame, at a level no exposure time predicts --
which is what a hand-typed ladder does when two different real exposures get
the same typed value.

    by shutter speed   the two are averaged; the tier is wrong
    per frame          each is calibrated on its own

Measured against a CONTROL whose headers are honest, so "right" is a real
profile and not an argument. The comparison is the merged luminance by radius,
normalised -- an overall scale is meaningless here (see smoke_blind_exposure),
the shape is not.

    python tools/smoke_tier_mode.py [-v]
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

_N_TIERS = 7
_STEP = 2.2
_T0 = 1 / 250.0
#: the frame whose real exposure is longer than the one it claims
_LIAR_TIER = 4
_LIAR_GAIN = 1.45
_EDGES = np.array([1.1, 1.4, 1.8, 2.3, 2.9, 3.6, 4.4])


def _scene(H=380, W=500, R=44.0, peak=7.0e5):
    cy, cx = H / 2.0, W / 2.0
    yy = np.arange(H, dtype=np.float64)[:, None] - cy
    xx = np.arange(W, dtype=np.float64)[None, :] - cx
    r = np.hypot(yy, xx)
    th = np.arctan2(yy, xx)
    thr = 1.0 + 0.30 * np.sin(3 * th + 0.7) + 0.10 * np.sin(11 * th)
    base = peak * (R / np.maximum(r, R)) ** 4 * thr
    base[r < R] = 0.0
    return base


def _write(folder, honest):
    """Two frames per tier. If `honest` is False the second frame of
    _LIAR_TIER really got _LIAR_GAIN more light than its header admits."""
    os.makedirs(folder, exist_ok=True)
    base = _scene()
    gain = np.array([0.70, 0.88, 1.00])[:, None, None]
    rng = np.random.default_rng(23)
    n = 0
    for i in range(_N_TIERS):
        t = _T0 * _STEP ** i
        for k in range(2):
            real = t
            if (not honest) and i == _LIAR_TIER and k == 1:
                real = t * _LIAR_GAIN
            a = base[None, :, :] * gain * real + 4.0
            cube = a + rng.normal(0, np.sqrt(np.maximum(a, 1)))
            write_cube(os.path.join(folder, "f%02d.fits" % n), cube, t,
                       "2026-08-12T18:%02d:%02d.0" % (n // 60, n % 60))
            n += 1
    return n


def _run(folder, mode, verbose):
    from eclipseforgehdr.pipeline import run
    run(folder, Prog(verbose), crop_pc=300, tier_mode=mode)
    wd = os.path.join(folder, ".eclipseforgehdr")
    geo = json.load(open(os.path.join(wd, "geometry.json")))
    rgb = np.load(os.path.join(wd, "hdr_rgb.npy"))
    st = json.load(open(os.path.join(wd, "report.json")))
    return st, rgb, geo


def _profile(rgb, cy, cx, R):
    lum = (0.2126 * rgb[..., 0] + 0.7152 * rgb[..., 1] + 0.0722 * rgb[..., 2])
    H, W = lum.shape
    yy = np.arange(H)[:, None] - cy
    xx = np.arange(W)[None, :] - cx
    rr = np.hypot(yy, xx) / R
    out = []
    for lo, hi in zip(_EDGES[:-1], _EDGES[1:]):
        m = (rr >= lo) & (rr < hi)
        out.append(float(np.median(lum[m])) if m.sum() >= 120 else np.nan)
    return np.asarray(out)


def _shape_err(p, ref):
    """Difference in SHAPE only -- an overall scale carries no information
    here, for the reason smoke_blind_exposure measures."""
    ok = np.isfinite(p) & np.isfinite(ref) & (ref > 0)
    if ok.sum() < 4:
        return float("nan")
    ratio = p[ok] / ref[ok]
    return float(np.max(ratio / np.median(ratio)) - np.min(ratio / np.median(ratio)))


def main(verbose=False):
    tmp = tempfile.mkdtemp(prefix="efhdr_tiermode_")
    try:
        cases = {}
        _write(os.path.join(tmp, "control"), honest=True)
        cases["control"] = _run(os.path.join(tmp, "control"), "exposure", verbose)
        for mode in ("exposure", "frame"):
            d = os.path.join(tmp, "liar_" + mode)
            _write(d, honest=False)
            cases["liar_" + mode] = _run(d, mode, verbose)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    st, _, geo = cases["control"]
    cy, cx, R = geo["cy"], geo["cx"], geo["R"]
    ref = _profile(cases["control"][1], cy, cx, R)

    print("bracket: %d tiers x%.1f, 2 frames each. In the two 'liar' runs the "
          "second frame of tier %d really got %.2fx the light its header claims."
          % (_N_TIERS, _STEP, _LIAR_TIER, _LIAR_GAIN))
    print("tiers formed: control %d | by shutter %d | per frame %d"
          % tuple(len(cases[k][0].get("tiers") or [])
                  for k in ("control", "liar_exposure", "liar_frame")))
    _lt = [t for t in (cases["liar_frame"][0].get("tiers") or [])
           if abs(float(t.get("cal", 1)) - _LIAR_GAIN) < 0.15]
    print("per-frame found the liar: %s"
          % (("cal %.4f on %s" % (_lt[0]["cal"], _lt[0]["best"])) if _lt
             else "NO — no frame came back near %.2f" % _LIAR_GAIN))
    print("clipping: the liar frame is %.0f%% over the well at the limb"
          % (100 * (_LIAR_GAIN * 7.0e5 * _T0 * _STEP ** _LIAR_TIER / 65535 - 1)))
    print()
    print("merged luminance by radius, as a ratio to the honest control")
    print("   r/R        by shutter    per frame")
    pe = _profile(cases["liar_exposure"][1], cy, cx, R)
    pf = _profile(cases["liar_frame"][1], cy, cx, R)
    re_ = pe / ref
    rf = pf / ref
    re_ = re_ / np.nanmedian(re_)
    rf = rf / np.nanmedian(rf)
    for k in range(len(_EDGES) - 1):
        print("  %.1f-%.1f %12.4f %12.4f" % (_EDGES[k], _EDGES[k + 1],
                                             re_[k], rf[k]))
    ee = _shape_err(pe, ref)
    ef = _shape_err(pf, ref)
    print()
    print("shape error against the honest control:  by shutter %.2f%%   "
          "per frame %.2f%%" % (100 * ee, 100 * ef))

    # WHAT THIS ASSERTS, AND WHAT IT DOES NOT.
    #
    # It does NOT assert that per-frame is more accurate. Two fixtures were
    # tried and neither showed a benefit: with the liar in an unclipped tier,
    # grouping cost nothing at all, because the photometric fit simply measured
    # whatever the average came to (cal 1.2258, the mean of 1.00 and 1.45) and
    # the levels stayed right. Moving it to the clipping boundary did not help
    # either -- per-frame came out slightly WORSE, which is the sqrt(2) of SNR
    # each tier loses when it stops averaging two frames, feeding noisier
    # photometric links.
    #
    # Tuning the fixture until it agreed would have been the RHEF mistake over
    # again, so it is left saying what it found. What this asserts is that the
    # mode WORKS: the frames are not grouped, each one is calibrated on its own
    # and the liar is found, and the result does not fall apart. Whether it
    # helps is a question about real data with untrustworthy headers, and the
    # selector exists so that question can be answered by looking.
    bad = []
    _nt = len(cases["liar_frame"][0].get("tiers") or [])
    if _nt != 2 * _N_TIERS:
        bad.append("per-frame mode made %d tiers, not %d — the frames were "
                   "still grouped" % (_nt, 2 * _N_TIERS))
    if not _lt:
        bad.append("per-frame did not find the liar: no frame came back with "
                   "cal near %.2f" % _LIAR_GAIN)
    if not np.isfinite(ee) or not np.isfinite(ef):
        bad.append("a profile could not be measured")
    elif ef > 0.03:
        bad.append("per-frame left %.1f%% of shape error — that is not a cost, "
                   "that is broken" % (100 * ef))
    print()
    if bad:
        for b in bad:
            print("FAIL:", b)
        return 1
    print("OK — per-frame grouping ran over %d frames, calibrated the liar on "
          "its own (cal %.4f against %.2f planted) and left %.2f%% of shape "
          "error against %.2f%% for grouping by shutter speed. It is not more "
          "accurate on this fixture and is not claimed to be."
          % (2 * _N_TIERS, _lt[0]["cal"], _LIAR_GAIN, 100 * ef, 100 * ee))
    return 0


if __name__ == "__main__":
    sys.exit(main("-v" in sys.argv))
