#!/usr/bin/env python3
"""End to end: does a planted per-channel black level come back out?

The unit bench for the estimator is `pedestal_probe.py`. This is the other
half -- the whole pipeline, from FITS on disk to the merged HDR, with a floor
planted in ONE channel and nothing else different about it.

WHAT IT WOULD MISS WITHOUT THIS TEST. The estimator can be right and the
correction still land in the wrong place: subtracted after white balance, or
after the camera matrix, or before the clipping test, or scaled by the exposure
somewhere on the way. Every one of those still produces a number that looks
sensible in the log. The only thing that settles it is the colour of the merged
picture at two radii, so that is what is measured.

THE MEASUREMENT. An additive floor is negligible where the corona is bright and
dominant where it is faint, so it does not shift colour uniformly -- it makes
colour a function of radius. The synthetic corona here is grey by construction:
the same spatial pattern in all three channels, scaled by a fixed per-channel
gain. So the ratio B/G must be the SAME at 1.5 R and at 3.2 R. The test is that
difference, before and after.

    python tools/smoke_channel_floors.py [-v]
"""
import os
import shutil
import sys
import tempfile

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, "..", "src"))
sys.path.insert(0, _HERE)

from smoke_pipeline import _P as Prog                 # noqa: E402


def _card(k, v):
    if isinstance(v, str):
        val = "'%-8s'" % v
    elif isinstance(v, bool):
        val = "T" if v else "F"
    else:
        val = "%20s" % repr(v)
    return ("%-8s= %-70s" % (k, val))[:80]


def write_cube(path, cube, exptime, when):
    """Planes-first (3, H, W) uncompressed 16-bit FITS.

    Not smoke_fits.write_fits: that one writes its own EXPTIME = 0.25 card
    before any extras, and a FITS reader takes the FIRST occurrence of a
    keyword -- so every frame came back as one tier and the fit saw a bracket
    of length 1. Caught by the estimator refusing with "1 tier(s) over 1.0x".
    """
    a = np.asarray(cube)
    d = (np.clip(np.rint(a), 0, 65535) - 32768).astype(">i2")
    cards = [_card("SIMPLE", True), _card("BITPIX", 16), _card("NAXIS", 3),
             _card("NAXIS1", a.shape[2]), _card("NAXIS2", a.shape[1]),
             _card("NAXIS3", a.shape[0]),
             _card("BZERO", 32768), _card("BSCALE", 1),
             _card("EXPTIME", float(exptime)), _card("DATE-OBS", when)]
    head = "".join(cards) + "END" + " " * 77
    head += " " * ((2880 - len(head) % 2880) % 2880)
    body = d.tobytes()
    body += b"\0" * ((2880 - len(body) % 2880) % 2880)
    with open(path, "wb") as f:
        f.write(head.encode("ascii")); f.write(body)



# Gains, so the three channels are not identical and the fit has something to
# separate. Roughly Val's measured rates.
_GAIN = (0.70, 0.88, 1.00)
#: planted floors, ADU. Blue carries an extra 100 on top of the shared 4.
_FLOOR = (0.0, 0.0, 100.0)
_SHARED = 4.0


def bracket(folder, H=420, W=560, R=46.0, n_tiers=8, peak=6.0e4):
    """Eight tiers over 24x, three planes, a grey corona and a blue floor."""
    secs = [(1 / 60.0) * (24.0 ** (i / (n_tiers - 1.0))) for i in range(n_tiers)]
    cy, cx = H / 2.0, W / 2.0
    yy = np.arange(H, dtype=np.float64)[:, None] - cy
    xx = np.arange(W, dtype=np.float64)[None, :] - cx
    r = np.hypot(yy, xx)
    th = np.arctan2(yy, xx)
    thr = 1.0 + 0.30 * np.sin(3 * th + 0.7) + 0.10 * np.sin(11 * th)
    base = peak * (R / np.maximum(r, R)) ** 4 * thr
    base[r < R] = 0.0
    rng = np.random.default_rng(11)
    n = 0
    for s in secs:
        for _ in range(2):
            cube = np.empty((3, H, W), np.float64)
            for c in range(3):
                a = base * _GAIN[c] * s + _SHARED + _FLOOR[c]
                cube[c] = a + rng.normal(0, np.sqrt(np.maximum(a, 1)))
            write_cube(os.path.join(folder, "f%02d.fits" % n), cube, s,
                       "2026-08-12T18:%02d:%02d.0" % (n // 60, n % 60))
            n += 1
    return n, secs


def _colour_by_radius(rgb, cy, cx, R):
    """(B/G at 1.5 R, B/G at 3.2 R). Equal is what grey data must give."""
    H, W = rgb.shape[:2]
    yy = np.arange(H)[:, None] - cy
    xx = np.arange(W)[None, :] - cx
    rr = np.hypot(yy, xx)
    out = []
    for lo, hi in ((1.3, 1.7), (3.0, 3.4)):
        m = (rr > lo * R) & (rr < hi * R)
        if m.sum() < 200:
            out.append(np.nan)
            continue
        g = float(np.median(rgb[..., 1][m]))
        b = float(np.median(rgb[..., 2][m]))
        out.append(b / g if abs(g) > 1e-9 else np.nan)
    return out


def run_once(tmp, disable, verbose):
    from eclipseforgehdr.pipeline import run
    folder = os.path.join(tmp, "off" if disable else "on")
    os.makedirs(folder, exist_ok=True)
    n, secs = bracket(folder)
    os.environ.pop("ECLIPSEFORGE_NO_CHANNEL_FLOORS", None)
    if disable:
        os.environ["ECLIPSEFORGE_NO_CHANNEL_FLOORS"] = "1"
    run(folder, Prog(verbose), crop_pc=320)   # returns None; read report.json
    import json
    wd = os.path.join(folder, ".eclipseforgehdr")
    st = json.load(open(os.path.join(wd, "report.json")))
    geo = json.load(open(os.path.join(wd, "geometry.json")))
    rgb = np.load(os.path.join(wd, "hdr_rgb.npy"))
    return st, _colour_by_radius(rgb, geo["cy"], geo["cx"], geo["R"])


def main(verbose=False):
    tmp = tempfile.mkdtemp(prefix="efhdr_chfloor_")
    try:
        st_off, off = run_once(tmp, True, verbose)
        st_on, on = run_once(tmp, False, verbose)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
        os.environ.pop("ECLIPSEFORGE_NO_CHANNEL_FLOORS", None)

    applied = st_on.get("channel_floors_applied")
    info = ((st_on.get("per_channel_photometry") or {}).get("channel_floors")
            or {})
    true_d = np.asarray(_FLOOR, float) - np.mean(_FLOOR)

    print("planted floors            R %+7.1f  G %+7.1f  B %+7.1f ADU" % _FLOOR)
    print("  as deltas about the mean R %+7.1f  G %+7.1f  B %+7.1f" % tuple(true_d))
    print("fitted and applied        %s" % (
        "  ".join("%+7.1f" % v for v in applied) if applied else "NOTHING"))
    print("  discarded as ladder-like %s" % info.get("deltas_ladder_part"))
    print("  fit residual             %.2f%%   ladder dev %s %%"
          % (100 * info.get("resid_frac", float("nan")),
             info.get("ladder_dev_pct")))
    print()
    print("B/G in the merged image    at 1.5 R   at 3.2 R    difference")
    print("  correction OFF          %9.4f %10.4f %13.1f%%"
          % (off[0], off[1], 100 * (off[1] / off[0] - 1)))
    print("  correction ON           %9.4f %10.4f %13.1f%%"
          % (on[0], on[1], 100 * (on[1] / on[0] - 1)))

    bad = []
    if not applied:
        bad.append("nothing was applied")
    d_off = abs(off[1] / off[0] - 1)
    d_on = abs(on[1] / on[0] - 1)
    if not (d_off > 0.05):
        bad.append("the planted floor did not produce a radial colour cast to "
                   "begin with (%.1f%%) -- the test proves nothing" % (100 * d_off))
    if not (d_on < 0.4 * d_off):
        bad.append("the cast was not reduced: %.1f%% -> %.1f%%"
                   % (100 * d_off, 100 * d_on))
    print()
    if bad:
        for b in bad:
            print("FAIL:", b)
        return 1
    print("OK — a planted %.0f ADU blue floor makes a %.1f%% radial colour cast, "
          "and the fit removes %.0f%% of it"
          % (_FLOOR[2], 100 * d_off, 100 * (1 - d_on / d_off)))
    return 0


if __name__ == "__main__":
    sys.exit(main("-v" in sys.argv))
