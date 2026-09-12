#!/usr/bin/env python3
"""Per-channel black level, fitted as slope AND intercept together.

WHY THIS IS NOT THE ESTIMATOR ALREADY IN _tier_colour_check.

`pipeline._fit_pedestal` searches one additive offset that minimises tier-to-
tier disagreement once each tier is divided by its own exposure. Run per
channel, it carries a warning written into the pipeline:

    given a synthetic bracket carrying a pure 3%-per-step COLOUR GAIN and no
    offset at all, it returns -17 and +17 ADU of offset for R and B

which is true and is the reason that path reports rather than corrects. An
offset and a per-step gain error both make the tiers disagree, and a fitter
that is given only an offset to play with will spend it on either.

This module fits the other way round. For one channel, over the tiers, the
outer-field ring median obeys

    m(t) = a * t + b

with `a` the scene rate in ADU/s and `b` the floor left in the data. Both are
free, so a colour GAIN lands in `a` where it belongs and only a genuine floor
lands in `b`. The degeneracy in the warning comes from fitting one of the two;
it does not survive fitting both.

The separation is not free, and the price is stated rather than hidden: a gain
error that COMPOUNDS along the ladder is not a straight line in t, so it bends
the fit and leaks into `b`. That bend is measurable, and it is what the
residual guard below is for. A true floor plus a linear scene fits a straight
line to a fraction of a percent; a compounding gain error does not.

    python3 pedestal_probe.py --self-test
    python3 pedestal_probe.py <folder-of-FITS>
"""
import os
import sys

import numpy as np


# THE GUARD IS CURVATURE, NOT RESIDUAL. This was got wrong once and the
# self-test caught it, so the reasoning is recorded here rather than rederived.
#
# The obvious guard is "does the straight line fit" -- max residual as a
# fraction of the largest median. It does not work. A compounding 3%-per-step
# colour gain, the exact case the pipeline warns about, misses the line by only
# 1.36% and sails through a 2% guard, returning -91 ADU of floor that is not
# there. Tightening the number does not help either: a 1%/step gain misses by
# 0.46%, which is inside the honest scatter of a real bracket.
#
# What separates them is the SHAPE of the miss. A floor plus a linear scene is
# a straight line and its quadratic term is noise. A gain error that compounds
# along the ladder is geometric in t, so it bends, and the bend is large:
#
#     case                      resid    bend, as a fraction of the fitted
#                                        signal at the longest exposure
#     pure floor 107 ADU        0.01%              0.004%
#     Val-like, real noise      0.10%              0.03%
#     1%/step gain, no floor    0.46%              3.7%
#     3%/step gain, no floor    1.36%             13.1%
#
# Two orders of magnitude, where the residual gives less than three. So the fit
# is believed on curvature and the residual is kept only as a second opinion.
_MAX_BEND_FRAC = 0.01           # quadratic term / fitted signal at t_max
_MAX_RESID_FRAC = 0.03          # of the largest median; a loose sanity check
_MIN_TIERS = 6                  # below this, a two-parameter fit is not a fit
_MIN_SPAN = 8.0                 # shortest-to-longest exposure ratio


def fit_floor(secs, medians, max_resid_frac=_MAX_RESID_FRAC):
    """(b, a, info) for one channel.

    `secs`     exposure time per tier, seconds
    `medians`  outer-field ring median per tier, raw ADU, same order

    Returns the fitted floor `b` in ADU, the scene rate `a` in ADU/s, and an
    info dict carrying everything needed to decide whether to believe it.
    Robust: one iteratively reweighted pass with a Huber weight, so a single
    bad tier bends the line less than it would in plain least squares.
    """
    t = np.asarray(secs, float)
    m = np.asarray(medians, float)
    ok = np.isfinite(t) & np.isfinite(m) & (t > 0)
    t, m = t[ok], m[ok]
    info = {"n": int(t.size), "span": float(t.max() / t.min()) if t.size else 0.0}
    if t.size < _MIN_TIERS or info["span"] < _MIN_SPAN:
        info["usable"] = False
        info["why"] = ("only %d tier(s) over %.1fx exposure; a two-parameter fit "
                       "needs %d over %.0fx" % (t.size, info["span"],
                                                _MIN_TIERS, _MIN_SPAN))
        return np.nan, np.nan, info

    w = np.ones_like(t)
    a = b = 0.0
    for _ in range(5):
        W = np.diag(w)
        A = np.vstack([t, np.ones_like(t)]).T
        sol, *_ = np.linalg.lstsq(W @ A, w * m, rcond=None)
        a, b = float(sol[0]), float(sol[1])
        r = m - (a * t + b)
        s = 1.4826 * np.median(np.abs(r - np.median(r)))
        if s <= 0:
            break
        z = np.abs(r) / s
        w = np.where(z <= 2.0, 1.0, 2.0 / np.maximum(z, 1e-9))

    r = m - (a * t + b)
    info["max_resid"] = float(np.abs(r).max())
    info["resid_frac"] = float(np.abs(r).max() / max(np.abs(m).max(), 1e-9))
    info["rate"] = a
    # curvature: does a quadratic in t buy anything? A compounding gain error
    # says yes, a floor says no. Reported as the quadratic term scaled to the
    # span, i.e. how many ADU the bend is worth end to end.
    A2 = np.vstack([t * t, t, np.ones_like(t)]).T
    sol2, *_ = np.linalg.lstsq(A2, m, rcond=None)
    info["bend_adu"] = float(abs(sol2[0]) * (t.max() ** 2))
    info["bend_frac"] = float(info["bend_adu"] / max(abs(a) * t.max(), 1e-9))
    info["usable"] = (info["bend_frac"] <= _MAX_BEND_FRAC
                      and info["resid_frac"] <= max_resid_frac)
    if info["bend_frac"] > _MAX_BEND_FRAC:
        info["why"] = ("the tiers bend away from a straight line by %.1f ADU "
                       "(%.1f%% of the signal at the longest exposure). A floor "
                       "plus a linear scene does not bend; a gain error that "
                       "compounds along the ladder does, and the two cannot be "
                       "told apart once it has. The intercept is not a black "
                       "level and is not used"
                       % (info["bend_adu"], 100 * info["bend_frac"]))
    elif not info["usable"]:
        info["why"] = ("the straight line misses by %.1f ADU (%.1f%% of the "
                       "largest median); too scattered to trust the intercept"
                       % (info["max_resid"], 100 * info["resid_frac"]))
    return b, a, info


def _synth(secs, rate, floor, gain_per_step=1.0, noise=0.0, seed=0):
    """One channel's ring medians. `gain_per_step` COMPOUNDS along the ladder,
    which is the failure mode the pipeline's warning describes."""
    rng = np.random.default_rng(seed)
    out = []
    for i, t in enumerate(secs):
        g = gain_per_step ** i
        out.append(rate * t * g + floor + rng.normal(0, noise))
    return np.asarray(out)


def self_test():
    secs = np.asarray([0.05 * 1.5 ** i for i in range(11)])
    print("%-46s %9s %9s %8s %8s  %s"
          % ("case", "fitted b", "true b", "resid%", "bend%", "verdict"))
    bad = 0

    def row(label, med, true_b, expect_usable, tol):
        nonlocal bad
        b, a, info = fit_floor(secs, med)
        ok = info["usable"] == expect_usable
        if info["usable"] and expect_usable:
            ok = ok and abs(b - true_b) <= tol
        bad += 0 if ok else 1
        print("%-46s %9.2f %9.2f %7.2f%% %7.2f%%  %s%s"
              % (label, b, true_b, 100 * info["resid_frac"],
                 100 * info["bend_frac"],
                 "used" if info["usable"] else "REFUSED",
                 "" if ok else "   <-- WRONG"))

    # 1. the control: a real floor, nothing else
    row("pure floor 107 ADU, clean",
        _synth(secs, 1662, 107.0, noise=0.5), 107.0, True, 3.0)
    # 2. the control in the other direction
    row("no floor at all",
        _synth(secs, 1662, 0.0, noise=0.5), 0.0, True, 3.0)
    # 3. THE WARNING: pure compounding colour gain, no floor whatsoever.
    #    The estimator must not invent a floor here -- or must refuse.
    row("pure 3%/step gain, floor 0  (the warning)",
        _synth(secs, 1662, 0.0, gain_per_step=1.03, noise=0.5), 0.0, False, 3.0)
    row("pure 1%/step gain, floor 0",
        _synth(secs, 1662, 0.0, gain_per_step=1.01, noise=0.5), 0.0, False, 3.0)
    # 4. both at once -- the honest hard case
    row("floor 107 AND 3%/step gain",
        _synth(secs, 1662, 107.0, gain_per_step=1.03, noise=0.5), 107.0, False, 3.0)
    # 5. Val's three channels, as measured
    for name, rate, floor in (("R", 1164.9, 102.73), ("G", 1462.3, 116.23),
                              ("B", 1661.9, 222.76)):
        row("Val-like %s: rate %.0f, floor %.1f" % (name, rate, floor),
            _synth(secs, rate, floor, noise=3.0), floor, True, 6.0)
    # 6. too few tiers
    b, a, info = fit_floor(secs[:4], _synth(secs[:4], 1662, 107.0))
    print("%-46s %9s %9.2f %8s %8s  %s%s"
          % ("only 4 tiers", "n/a", 107.0, "-", "-",
             "REFUSED" if not info["usable"] else "used",
             "" if not info["usable"] else "   <-- WRONG"))
    bad += 0 if not info["usable"] else 1

    print()
    print("SELF TEST %s" % ("PASSED" if bad == 0 else "FAILED (%d)" % bad))
    return bad == 0


if __name__ == "__main__":
    if "--self-test" in sys.argv:
        sys.exit(0 if self_test() else 1)
    print(__doc__)


# ---------------------------------------------------------------------------
# The estimator that survives contact with a real bracket.
#
# fit_floor above refuses Val's set, and it is right to: over 0.05-2.96 s the
# three channels bend away from a straight line by 1.9%, 2.3% and 2.8% of the
# signal. A bend is a gain error compounding along the ladder, and by the
# argument at the top of this file a fitter that cannot see it will spend it on
# the intercept.
#
# But look at what bends. All THREE channels bend by nearly the same amount.
# An exposure-time error is achromatic -- a shutter does not change colour --
# so a common bend is a property of the ladder, not of the sensor, and it
# cannot produce a colour cast. What produces a colour cast is the channels
# having DIFFERENT floors.
#
# So the ladder is made a free parameter instead of a nuisance:
#
#     m_c(i) = a_c * s_i + b_c
#
# with `s_i` one effective exposure per tier SHARED by the three channels, and
# (a_c, b_c) a rate and a floor per channel. Any achromatic error in the stated
# exposure -- a wrong EXPTIME, a shutter that did not honour the request, a
# compounding gain -- lands in `s_i`, where it belongs, and cannot reach `b_c`.
#
# ONLY THE DIFFERENCES BETWEEN CHANNELS ARE RETURNED. The absolute level of the
# floors is constrained only through the channels having different rates, which
# is a weak lever; the differences are constrained directly by the channels
# disagreeing at the faint end, which is a strong one. The common part is
# already fitted, jackknifed and shrunk by `pipeline._fit_pedestal`, and this
# does not touch it.

def fit_floor_deltas(secs, med_rgb, iters=200):
    """Per-channel floor DELTAS about their own mean, in ADU.

    `med_rgb` is (3, n_tiers) of outer-field ring medians in raw ADU.
    Returns (deltas, s_eff, info): deltas sum to zero, `s_eff` is the fitted
    effective exposure ladder normalised to match `secs` in a least-squares
    sense, and info carries the fit quality and the guard verdict.
    """
    t = np.asarray(secs, float)
    M = np.asarray(med_rgb, float)
    info = {"n": int(t.size), "span": float(t.max() / t.min()) if t.size else 0.0}
    if M.shape[0] != 3 or t.size < _MIN_TIERS or info["span"] < _MIN_SPAN:
        info["usable"] = False
        info["why"] = ("needs 3 channels over at least %d tiers and %.0fx "
                       "exposure; got %d tiers over %.1fx"
                       % (_MIN_TIERS, _MIN_SPAN, t.size, info["span"]))
        return np.zeros(3), t.copy(), info

    s = t.copy()
    a = np.zeros(3); b = np.zeros(3)
    for _ in range(iters):
        A = np.vstack([s, np.ones_like(s)]).T
        for c in range(3):
            sol, *_ = np.linalg.lstsq(A, M[c], rcond=None)
            a[c], b[c] = float(sol[0]), float(sol[1])
        den = float((a * a).sum())
        if den <= 0:
            break
        s_new = ((a[:, None] * (M - b[:, None])).sum(axis=0)) / den
        # kill the scale degeneracy: s and a trade off, so pin the scale by
        # least squares against the stated ladder and let the SHAPE be free
        k = float((s_new * t).sum() / max((s_new * s_new).sum(), 1e-30))
        s_new = s_new * k
        if np.max(np.abs(s_new - s)) < 1e-9 * max(1.0, float(t.max())):
            s = s_new
            break
        s = s_new

    # PROJECT OUT WHAT A LADDER ERROR COULD HAVE DONE.
    #
    # (s, b) and (s + d, b - a*d) fit identically for any constant d: a shift
    # in the effective exposure -- shutter lag, a delay the capture software
    # did not report -- is indistinguishable from per-channel floors that are
    # PROPORTIONAL TO EACH CHANNEL'S RATE. The self-test shows this directly:
    # given a pure achromatic gain and no floors at all, the raw deltas come
    # out (+14.2, -1.8, -12.4), and dividing those by the mean-removed rates
    # gives (-0.056, -0.054, -0.056) -- one number, i.e. exactly along `a`.
    #
    # That component is therefore not evidence of a floor and is discarded.
    # What is left is orthogonal to `a`: a difference between the channels
    # that no achromatic error in the exposure can produce. Under-correcting
    # when the true floors happen to lie partly along `a` is the price, and it
    # is the right way round -- this can fail to remove a real cast, but it
    # cannot invent one out of a ladder error.
    R = M - (a[:, None] * s + b[:, None])
    info["max_resid"] = float(np.abs(R).max())
    info["resid_frac"] = float(np.abs(R).max() / max(np.abs(M).max(), 1e-9))
    info["rates"] = a.tolist()
    info["floors_abs"] = b.tolist()
    # how far the fitted ladder had to move from the stated one
    info["ladder_max_dev"] = float(np.max(np.abs(s / np.maximum(t, 1e-12) - 1.0)))
    info["usable"] = info["resid_frac"] <= _MAX_RESID_FRAC
    if not info["usable"]:
        info["why"] = ("even with the ladder free the three channels miss by "
                       "%.1f ADU (%.1f%%); something other than a floor and an "
                       "achromatic exposure error is present"
                       % (info["max_resid"], 100 * info["resid_frac"]))
    braw = b - b.mean()
    arel = a - a.mean()
    den = float((arel * arel).sum())
    proj = (float((braw * arel).sum()) / den) * arel if den > 1e-12 else 0.0 * arel
    info["deltas_raw"] = braw.tolist()
    info["deltas_ladder_part"] = np.asarray(proj).tolist()
    return braw - proj, s, info


def self_test_deltas():
    secs = np.asarray([0.05 * 1.5 ** i for i in range(11)])
    rates = np.array([1140.9, 1428.4, 1615.1])
    print("deltas below are the part ORTHOGONAL to the channel rates --")
    print("the part no achromatic exposure error can explain.\n")
    print("%-48s %22s %22s  %s"
          % ("case", "fitted R G B", "recoverable truth R G B", "verdict"))
    bad = 0

    def row(label, floors, gain_per_step=1.0, noise=1.0, expect=True, tol=4.0):
        nonlocal bad
        M = np.vstack([_synth(secs, rates[c], floors[c],
                              gain_per_step=gain_per_step, noise=noise, seed=c)
                       for c in range(3)])
        d, s, info = fit_floor_deltas(secs, M)
        true = np.asarray(floors, float) - np.mean(floors)
        _ar = rates - rates.mean()
        true = true - (float((true * _ar).sum()) / float((_ar * _ar).sum())) * _ar
        ok = info["usable"] == expect
        if expect and info["usable"]:
            ok = ok and np.max(np.abs(d - true)) <= tol
        bad += 0 if ok else 1
        print("%-48s %22s %22s  %s%s"
              % (label,
                 " ".join("%6.1f" % v for v in d),
                 " ".join("%6.1f" % v for v in true),
                 "used" if info["usable"] else "REFUSED",
                 "" if ok else "   <-- WRONG"))

    row("Val-like floors, honest ladder", [102.7, 116.2, 222.8])
    row("Val-like floors, 3%/step ACHROMATIC gain",
        [102.7, 116.2, 222.8], gain_per_step=1.03)
    row("Val-like floors, 10%/step achromatic gain",
        [102.7, 116.2, 222.8], gain_per_step=1.10)
    row("no per-channel difference at all", [115.0, 115.0, 115.0])
    row("no floors at all", [0.0, 0.0, 0.0])
    row("pure achromatic gain, no floors  (the warning)",
        [0.0, 0.0, 0.0], gain_per_step=1.03)
    print()
    print("DELTA SELF TEST %s" % ("PASSED" if bad == 0 else "FAILED (%d)" % bad))
    return bad == 0
