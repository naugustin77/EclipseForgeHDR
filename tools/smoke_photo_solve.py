#!/usr/bin/env python3
"""Chain against network: does the extra link actually buy anything?

TWO PARTS, because one bracket cannot answer this.

1. THE INTEGRATION RUN scales one tier by a known factor and runs the whole
   pipeline both ways. It proves the setting reaches the solve, that two-step
   links are measurable on real-shaped data, and that the residuals are
   recorded. It does NOT prove the network is more accurate, and the first
   version of this file wrongly asserted that it would: a uniformly bright
   tier is not a CONTRADICTION. Its two adjacent links read 1.18 and 0.85, the
   two-step link across it reads 1.00, and those three agree -- so there is
   nothing to arbitrate and both paths correctly place the tier high. That is
   the right answer and it says nothing about the solve.

2. THE MONTE CARLO tests the actual claim. What a network buys is not outlier
   rejection but the averaging of independent MEASUREMENT error: the chain
   believes whichever links happen to lie on its path from the reference tier,
   so its error compounds with distance, while the network is constrained from
   both sides and by the two-step links as well. That is a statistical claim
   and needs many draws, not one bracket.

    python tools/smoke_photo_solve.py [-v]
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

_N_TIERS = 9
_STEP = 2.0
_T0 = 1 / 500.0
#: which tier is poisoned, and by how much. 1.18 is inside the 2.5x guard, so
#: the links are kept and have to be arbitrated rather than dropped.
_BAD_TIER = 3
_BAD_GAIN = 1.18


def _frames(H=420, W=560, R=46.0, peak=9.0e5):
    cy, cx = H / 2.0, W / 2.0
    yy = np.arange(H, dtype=np.float64)[:, None] - cy
    xx = np.arange(W, dtype=np.float64)[None, :] - cx
    r = np.hypot(yy, xx)
    th = np.arctan2(yy, xx)
    thr = 1.0 + 0.30 * np.sin(3 * th + 0.7) + 0.10 * np.sin(11 * th)
    base = peak * (R / np.maximum(r, R)) ** 4 * thr
    base[r < R] = 0.0
    gain = np.array([0.70, 0.88, 1.00])[:, None, None]
    rng = np.random.default_rng(5)
    out = []
    for i in range(_N_TIERS):
        t = _T0 * _STEP ** i
        for _ in range(2):
            a = base[None, :, :] * gain * t + 4.0
            if i == _BAD_TIER:
                a = a * _BAD_GAIN
            out.append((a + rng.normal(0, np.sqrt(np.maximum(a, 1))), t))
    return out


def _write(folder, frames):
    os.makedirs(folder, exist_ok=True)
    for n, (cube, t) in enumerate(frames):
        write_cube(os.path.join(folder, "f%02d.fits" % n), cube, t,
                   "2026-08-12T18:%02d:%02d.0" % (n // 60, n % 60))


def _run(folder, mode, verbose):
    from eclipseforgehdr.pipeline import run
    run(folder, Prog(verbose), crop_pc=320, photo_solve=mode)
    wd = os.path.join(folder, ".eclipseforgehdr")
    return json.load(open(os.path.join(wd, "report.json")))


def main(verbose=False):
    tmp = tempfile.mkdtemp(prefix="efhdr_solve_")
    try:
        frames = _frames()
        got = {}
        for mode in ("chain", "network"):
            d = os.path.join(tmp, mode)
            _write(d, frames)
            got[mode] = _run(d, mode, verbose)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    secs = [_T0 * _STEP ** i for i in range(_N_TIERS)]
    print("bracket: %d tiers x%.1f, tier %d (%.4gs) scaled by %.2f"
          % (_N_TIERS, _STEP, _BAD_TIER, secs[_BAD_TIER], _BAD_GAIN))
    net = got["network"].get("photo_network") or {}
    print("solve modes reported: chain=%s network=%s"
          % (got["chain"].get("photo_solve"), got["network"].get("photo_solve")))
    print("network: %s links, %s of them two-step"
          % (net.get("links", "-"), net.get("lag2", "-")))
    if net.get("resid_pct"):
        _r = np.abs(np.asarray(net["resid_pct"], float))
        print("         residual per link: median %.2f%%, worst %.2f%%"
              % (np.median(_r), _r.max()))
    print()

    # cal[] is not in report.json by tier, but the links are -- and the links
    # are what the two paths disagree about. Compare the cumulative factors
    # each path would produce from its own links against the truth.
    rows = {}
    for mode in ("chain", "network"):
        lk = np.asarray(got[mode].get("photometric_links") or [], float)
        rows[mode] = lk
    print("measured adjacent links (1.000 = the exposure ratio predicts it)")
    for mode in ("chain", "network"):
        print("  %-8s %s" % (mode, ", ".join("%.3f" % v for v in rows[mode])))
    print()
    print("  the two links touching tier %d should read about %.3f and %.3f"
          % (_BAD_TIER, _BAD_GAIN, 1.0 / _BAD_GAIN))

    rc = monte_carlo()

    bad = []
    if got["network"].get("photo_solve") != "network":
        bad.append("the network path did not run (%s) — nothing was compared"
                   % got["network"].get("photo_solve"))
    elif not net.get("lag2"):
        bad.append("no two-step links were measurable, so the network reduced "
                   "to the chain and this fixture proves nothing")
    elif not net.get("resid_pct"):
        bad.append("the network ran but recorded no residuals")
    print()
    if bad:
        for b in bad:
            print("FAIL:", b)
        return 1
    print("OK — the network ran over %d links (%d two-step), residuals near "
          "zero because a uniformly bright tier is not a contradiction: its "
          "two adjacent links and the two-step link across it agree. Both "
          "paths place it high, correctly." % (net["links"], net["lag2"]))
    return rc


def monte_carlo(trials=400, n=9, sigma=0.03, seed=3):
    """The claim the integration run above CANNOT test.

    A uniformly bright tier is not a contradiction, so chain and network agree
    on it -- which is the right answer and proves nothing about the solve. What
    the network is actually for is independent MEASUREMENT error on each link:
    the chain believes whichever links lie on its path, the network averages
    all of them. That is a statistical claim, so it needs many draws, not one
    bracket.

    Truth here is logf = 0 for every tier. Each link is measured with Gaussian
    error of `sigma` in the log. The chain gets adjacent links only, because
    that is all it can use; the network gets adjacent and two-step.
    """
    from eclipseforgehdr.pipeline import _solve_photo_network
    rng = np.random.default_rng(seed)
    ref = n // 2
    ec, en = [], []
    for _ in range(trials):
        lag1 = {(i, i + 1): rng.normal(0.0, sigma) for i in range(n - 1)}
        lag2 = {(i, i + 2): rng.normal(0.0, sigma) for i in range(n - 2)}
        logf = np.zeros(n)
        for i in range(ref, n - 1):
            logf[i + 1] = logf[i] + lag1[(i, i + 1)]
        for i in range(ref, 0, -1):
            logf[i - 1] = logf[i] - lag1[(i - 1, i)]
        ec.append(np.abs(logf))
        links = ([(i, j, v, 1e6) for (i, j), v in lag1.items()]
                 + [(i, j, v, 1e6) for (i, j), v in lag2.items()])
        sol, _ = _solve_photo_network(n, links, ref)
        en.append(np.abs(sol))
    ec = np.asarray(ec); en = np.asarray(en)
    print()
    print("MONTE CARLO — %d brackets of %d tiers, every link measured with "
          "%.1f%% error" % (trials, n, 100 * sigma))
    print("   truth is 1.000 for every tier; the table is the error, in %")
    print("   tier      chain    network")
    for k in range(n):
        print("   %4d   %8.2f %10.2f"
              % (k, 100 * (np.exp(ec[:, k].mean()) - 1),
                 100 * (np.exp(en[:, k].mean()) - 1)))
    wc = 100 * (np.exp(ec.mean()) - 1); wn = 100 * (np.exp(en.mean()) - 1)
    print("   mean   %8.2f %10.2f   -> %.0f%% less error"
          % (wc, wn, 100 * (1 - wn / wc)))
    if wn >= wc:
        print("   FAIL: the network is not better (%.2f%% vs %.2f%%)" % (wn, wc))
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main("-v" in sys.argv))
