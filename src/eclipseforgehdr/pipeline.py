"""Full pipeline: quality selection, alignment, calibration, HDR merge, layers.

All intermediate products are cached in <folder>/.eclipseforgehdr/ so re-renders
and exports don't repeat the heavy work.
"""
from __future__ import annotations
import os, json, time
import numpy as np
from scipy import ndimage
from skimage.registration import phase_cross_correlation

from . import align

# A prominence link is trusted only while its anchors agree with each other.
# Measured on the reference set the spread stays under ~1.9 px out to 1/13 s
# and then climbs to 4-9 px as inner-corona glare swallows the prominences,
# so the gate sits just above the good regime rather than at a tier index.
PROM_MAX_SPREAD = 2.0
PROM_WEIGHT = 1.0


def _exp_name(s):
    return f"1/{1 / s:.0f}s" if s < 0.4 else f"{s:g}s"


def input_fingerprint(folder):
    """Identity of the raw files a cached run was built from.

    The cache used to be validated on the render OPTIONS alone, so adding a
    frame to the folder, deleting a bad one, or replacing a file and pressing
    Start (without force) silently reused the previous stack -- and the report
    exported alongside it quoted the previous frame count.
    """
    from .raw import list_raws
    out = []
    for p in sorted(list_raws(folder)):
        try:
            st = os.stat(p)
            out.append([os.path.basename(p), int(st.st_size), int(st.st_mtime)])
        except OSError:
            out.append([os.path.basename(p), -1, -1])
    return out


from .raw import (open_frame, list_raws, read_exif, demosaic_rggb,
                  hot_pixel_map, repair_hot, read_camera_info)
from . import report as _report
from . import __version__


# Where the detail stage begins, in the fractions the pipeline code writes
# (_BAR_PIVOT) and on the bar the user sees (_BAR_DETAIL). Timing on synthetic
# frames (see detail.py) puts the detail stage at 22 minutes on a 43 Mpx frame,
# so handing it the last 6.5% of the bar was never going to look like progress.
# Rather than renumber every frac in this file, the two scales are joined by one
# piecewise-linear map: everything before the pivot is compressed into
# 0.._BAR_DETAIL, everything after is stretched over the rest.
#
# 0.14.1 set _BAR_DETAIL to 0.78 on an assumption -- that stacking costs about
# 3.5x the detail stage -- because the stacking side had never been timed on real
# RAWs. The run-time summary added in the same release then measured it, on a
# 50-frame 14-tier bracket from a 45 Mpx body, 12m27s end to end:
#
#     141s  inner corona: raw pass
#     140s  inner corona: denoised pass
#      31s  MGN
#     ----
#     312s  = 42% of the run in the three NAMED detail steps alone
#
# The other eight detail steps are each below the report's 6th-slowest entry
# (21s), which brackets the whole stage between 42% and 64% of the run. So
# stacking and detail are roughly equal, not 3.5:1, and 0.78 was out by more
# than a factor of two: the bar reached 78% around the halfway mark and crawled
# the rest. 0.52 sits at the low end of the measured bracket, because a bar that
# lags is better than one that arrives at 100% and waits.
#
# One real dataset, so this is an estimate with a range, not a constant. It will
# move with frame count (more frames = more stacking) and with sensor size
# (detail grows faster than pixel count). The report now lists every step over a
# second, so the next few runs can narrow it without guesswork.
#
# A run that IMPORTS a finished HDR skips the stacking entirely, so the detail
# stage is nearly the whole job there; that path lowers bar_detail on its own
# Progress rather than pretending the first half of its bar meant something.
_BAR_PIVOT = 0.935
_BAR_DETAIL = 0.52


class Progress:
    """Log lines plus the two clocks the GUI needs to show it is alive.

    A progress bar alone cannot distinguish "working on a slow step" from
    "hung": the merge sits at one frac for minutes on a big set. t0 and t_line
    let the GUI say how long the run has been going and how long since anything
    last happened, which is the honest liveness signal -- a spinner on its own
    only proves the browser is alive.
    """

    def __init__(self):
        self.lines = []
        self.stamps = []          # seconds since t0, one per line
        self.frac = 0.0
        self.done = False
        self.error = None
        self.t0 = time.time()
        self.t_line = self.t0
        self.bar_pivot = _BAR_PIVOT
        self.bar_detail = _BAR_DETAIL

    def bar(self, f):
        """Authored frac -> the frac the bar shows. See _BAR_DETAIL above."""
        f = float(np.clip(f, 0.0, 1.0))
        p, d = self.bar_pivot, self.bar_detail
        if f <= p:
            return f * (d / p)
        return d + (f - p) * (1.0 - d) / (1.0 - p)

    def log(self, msg, frac=None):
        self.t_line = time.time()
        self.lines.append(msg)
        self.stamps.append(self.t_line - self.t0)
        if frac is not None:
            # never go backwards: a stage that reports a lower frac than one
            # already passed would make the bar jump back and read as a restart
            self.frac = max(self.frac, self.bar(frac))

    def elapsed(self):
        return time.time() - self.t0

    def since(self):
        return time.time() - self.t_line


def _fmt_dur(s):
    s = int(round(float(s)))
    return f"{s // 60}m{s % 60:02d}s" if s >= 60 else f"{s}s"


def _timing_summary(progress, top=5, floor=1.0):
    """Where this run's wall time went, from the progress line timestamps.

    Every step announces itself before it works, so the cost of the step named
    by line i is the gap to line i+1.

    The first version of this reported only the six slowest, which was enough to
    show that the progress bar's split was wrong but not enough to say by how
    much: the named steps accounted for 42% of the run and everything else was
    only known to be "under 21 seconds each", bracketing the answer between 42%
    and 64%. So `steps` now carries every step over a second -- that is what
    makes the next run's weights a measurement instead of another estimate --
    while `slowest` stays short for the one-line log message.
    """
    lines = getattr(progress, "lines", None) or []
    st = getattr(progress, "stamps", None) or []
    if len(st) != len(lines) or len(st) < 3:
        return None
    def _label(m):
        return m.strip().rstrip(".:").split(",")[0][:52]

    durs = sorted(((st[i] - st[i - 1], lines[i - 1]) for i in range(1, len(st))),
                  key=lambda t: -t[0])
    steps = [[round(d, 1), _label(m)] for d, m in durs if d >= floor]
    if not steps:
        return None
    acc = sum(d for d, _ in steps)
    return {"total_s": round(st[-1], 1),
            "slowest": steps[:top],
            "steps": steps,
            # what the listed steps add up to, so a reader can see how much of
            # the run is accounted for rather than having to assume it is all
            "accounted_s": round(acc, 1)}


def workdir(folder):
    d = os.path.join(folder, ".eclipseforgehdr")
    os.makedirs(d, exist_ok=True)
    return d


# ---------- helpers on half-res superpixel luma ----------

def half_luma(bayer):
    h, w = bayer.shape
    b = bayer[: h // 2 * 2, : w // 2 * 2].reshape(h // 2, 2, w // 2, 2)
    return b.mean(axis=(1, 3))


def find_disc(lum, target=520.0):
    """Locate the occulting disc from the LIMB EDGE, not from brightness.

    Returns (cy, cx, R) in the input's own pixels, or None.

    Why not brightness. The old centre estimator took the centroid of the
    brightest 0.05% of pixels. That is only the disc centre while the bright
    inner corona rings the disc evenly. When one sector dominates -- a big
    prominence, an active region, a lopsided inner corona -- every one of those
    pixels sits on the same arc and the centroid lands ON THE LIMB. Measured on
    the reference set: 646 px out on a 620 px disc, i.e. 1.04 R. The half-level
    fit downstream tolerates a seed up to 0.81 R off and fails past ~1.0 R, so
    that estimator was sitting a few percent from a cliff, and a 3% change in
    the limb ring (a flat-field correction) was enough to push it over. The
    result was a disc mask 1000 px from the Moon.

    What replaces it is geometric. The limb is a huge RELATIVE step -- 93x over
    ~20 px on the reference set -- while the corona's own falloff is smooth, so
    in log intensity the limb dominates the gradient no matter how the frame is
    exposed or how bright one side is. Strong-gradient pixels are collected and
    a circle is fitted to them with iterative trimming. Nothing here refers to
    the frame size, so it does not care whether the disc is 3% or 30% of the
    frame: the radius comes out of the data. The one thing demanded of the
    answer is that the inliers go most of the way ROUND, which is what stops a
    streamer edge or a frame border from passing as a limb.
    """
    h, w = lum.shape
    dec = max(1, int(round(min(h, w) / float(target))))
    s = np.asarray(lum[::dec, ::dec], np.float32)
    s = ndimage.gaussian_filter(np.nan_to_num(s, nan=0.0, posinf=0.0, neginf=0.0), 1.5)
    lo = float(np.percentile(s, 1.0))
    span = float(np.percentile(s, 99.9)) - lo
    if not np.isfinite(span) or span <= 0:
        return None
    # log, floored well below the corona: a RELATIVE contrast measure, so it is
    # independent of exposure, of units and of the sky level
    L = np.log(np.maximum(s - lo, span * 1e-4) + span * 1e-4)
    gy, gx = np.gradient(L)
    g = np.hypot(gy, gx)
    hs, ws = s.shape

    # Start from a wide bright mask -- not the extreme tail. Measured across the
    # sweep this lands 0-15 px from the centre, which is plenty to begin with.
    thr = float(np.percentile(s, 98.0))
    ys, xs = np.nonzero(s >= thr)
    if ys.size < 20:
        return None
    cy, cx = float(ys.mean()), float(xs.mean())

    # Then alternate: measure the radius from the azimuthal gradient profile
    # about the current centre, keep only the strong-gradient pixels in a band
    # around THAT radius, and refit the centre to them. Restricting to the band
    # is what makes this work at any disc size -- it is what stops streamers and
    # the corona's own falloff from being fitted along with the limb.
    R = None
    for _ in range(4):
        Rn = _limb_radius_from_gradient(g, cy, cx, R or 0.25 * min(hs, ws))
        if Rn is None or not (2.0 < Rn < 0.48 * min(hs, ws)):
            return None
        R = Rn
        gthr = float(np.percentile(g, 96.0))
        yy = np.arange(hs, dtype=np.float32)[:, None] - cy
        xx = np.arange(ws, dtype=np.float32)[None, :] - cx
        rr = np.hypot(yy, xx)
        band = (g >= gthr) & (np.abs(rr - R) < max(0.25 * R, 3.0))
        py, px = np.nonzero(band)
        if py.size < 40:
            return None
        y = py.astype(np.float64)
        x = px.astype(np.float64)
        for _t in range(4):
            A = np.stack([2 * x, 2 * y, np.ones(x.size)], axis=1)
            try:
                sol, *_ = np.linalg.lstsq(A, x * x + y * y, rcond=None)
            except np.linalg.LinAlgError:
                return None
            ccx, ccy = float(sol[0]), float(sol[1])
            rad = sol[2] + ccx * ccx + ccy * ccy
            if not np.isfinite(rad) or rad <= 0:
                return None
            cR = float(np.sqrt(rad))
            d = np.abs(np.hypot(x - ccx, y - ccy) - cR)
            keep = d <= max(2.0 * float(np.median(d)), 1.0)
            if keep.sum() < 40 or keep.all():
                break
            x, y = x[keep], y[keep]
        if not (0 <= ccy < hs and 0 <= ccx < ws):
            return None
        cy, cx = ccy, ccx

    # Validate. Coverage alone is not enough: where the edge set is noise rather
    # than structure -- the short-exposure stack, whose fast tiers see nothing
    # outside the inner corona -- noise covers every azimuth and a meaningless
    # circle scores a full 360 degrees. The SCATTER of the inliers about the
    # circle is what separates them: measured on the reference set a real limb
    # sits at 0.04 R and the noise-fitted circles at 0.16-0.23 R.
    ang = np.arctan2(y - cy, x - cx)
    bins = np.zeros(36, bool)
    bins[((ang + np.pi) / (2 * np.pi) * 36).astype(int) % 36] = True
    med = float(np.median(np.abs(np.hypot(x - cx, y - cy) - cR)))
    if bins.mean() < 0.75 or med > max(0.10 * cR, 1.5):
        return None
    if not (0.008 * min(hs, ws) < cR < 0.48 * min(hs, ws)):
        return None
    return cy * dec, cx * dec, cR * dec


def _limb_radius_from_gradient(g, cy, cx, R_hint):
    """Radius of the strongest closed edge about (cy, cx): the limb."""
    h, w = g.shape
    rmax = min(cy, h - 1 - cy, cx, w - 1 - cx)
    rmax = float(min(max(rmax, 4.0), 0.48 * min(h, w), 3.0 * max(R_hint, 1.0)))
    if rmax < 6:
        return None
    yy = np.arange(h, dtype=np.float32)[:, None] - cy
    xx = np.arange(w, dtype=np.float32)[None, :] - cx
    rr = np.hypot(yy, xx)
    m = rr < rmax
    idx = rr[m].astype(np.int32)
    tot = np.bincount(idx, weights=g[m].astype(np.float64), minlength=int(rmax) + 1)
    cnt = np.bincount(idx, minlength=int(rmax) + 1).astype(np.float64)
    prof = tot / np.maximum(cnt, 1.0)
    prof[:3] = 0.0
    prof[cnt < 8] = 0.0
    if not prof.any():
        return None
    k = int(np.argmax(prof))
    # parabolic interpolation on the peak, so the answer is not quantised to
    # whole decimated pixels
    if 0 < k < len(prof) - 1:
        a, b, c = prof[k - 1], prof[k], prof[k + 1]
        den = a - 2 * b + c
        if den != 0:
            k = k + 0.5 * (a - c) / den
    return float(k)


def find_center(lum):
    """Centre of the occulting disc. Edge geometry first, brightness only as a
    fallback -- and then over a wide bright mask, never the extreme tail.

    On the reference merged frame the centroid of the brightest 0.05% lands 646
    px from the centre; widening the mask to the brightest 2% brings it to 81
    px, and every value between 1% and 20% lands inside 170 px. The tail was the
    whole problem: it is the part of the distribution a single bright feature
    can own outright.
    """
    d = find_disc(lum)
    if d is not None:
        return d[0], d[1]
    sm = ndimage.gaussian_filter(lum, 8)
    thr = np.percentile(sm, 98.0)
    ys, xs = np.nonzero(sm >= thr)
    if ys.size == 0:
        return lum.shape[0] / 2.0, lum.shape[1] / 2.0
    return float(ys.mean()), float(xs.mean())


def crop_around(img, cy, cx, size):
    h, w = img.shape
    y0 = int(np.clip(cy - size // 2, 0, max(h - size, 0)))
    x0 = int(np.clip(cx - size // 2, 0, max(w - size, 0)))
    return img[y0:y0 + size, x0:x0 + size], y0, x0


def sharpness(crop, satmask):
    valid = ~ndimage.binary_dilation(satmask, iterations=3)
    hf = crop - ndimage.gaussian_filter(crop, 2.0)
    sig = ndimage.gaussian_filter(crop, 4.0)
    lvl = np.percentile(sig, 90)
    m = valid & (sig > 0.25 * lvl)
    if m.sum() < 1000:
        m = valid
    return float(np.mean(hf[m] ** 2)) / max(float(np.mean(sig[m])), 1e-6) ** 2 * 1e4


def prep_pc(crop):
    x = np.log1p(np.clip(crop, 0, None))
    x = x - ndimage.gaussian_filter(x, 30)
    return x * np.hanning(x.shape[0])[:, None] * np.hanning(x.shape[1])[None, :]


def pc_shift(ref, mov, upsample=20):
    sh, err, _ = phase_cross_correlation(prep_pc(ref), prep_pc(mov),
                                         upsample_factor=upsample, normalization=None)
    return sh, err


# NOTE ON THE SIGN OF abs_shift (fixed in 0.8.8)
#
# phase_cross_correlation(ref, mov) returns the shift that must be APPLIED to
# `mov` to register it with `ref`. The network below solves for exactly that
# quantity per tier, so the shift is applied as-is: ndimage.shift(img, +shift).
#
# Every consumer used to negate it. That does not merely fail to align -- it
# leaves a residual of TWICE the true offset, which is worse than doing nothing.
# It went unnoticed because all six consumers negated it consistently, so the
# moon track, the moon-mask trial and the alignment-quality numbers were all
# computed in the same wrong frame and agreed with each other.
#
# Measured on the reference set, flipping to +shift:
#   moon-track scatter about the straight line   23.5 px -> 2.0 px
#   tier-to-tier variance at the limb            0.377  -> 0.073
#   tier-to-tier variance in the corona          0.165  -> 0.136
#   merged limb 20-80% transition                23.0px -> 15.5px
# The Moon's position in a corona-aligned frame is orbital motion and must lie
# on a straight line; 2 px of scatter is that line, 23 px was the bug.
#
# If you touch one of these, touch all of them.

def _moon_track(tier_moon, tier_time, progress):
    """Robust straight-line fit of the Moon's position and radius against time.

    In the corona-aligned frame the Moon's motion is smooth, monotonic and
    linear over one bracket -- it is orbital mechanics, not something that can
    jitter. So the per-tier limb measurements should lie on a line, and the
    scatter about that line is measurement error, not motion.

    Using the raw per-tier values as mask centres is what produced two lunar
    discs in 0.8.2: the scatter is 14-20 px and it clusters, so the exclusion
    discs land in two groups. Replacing each tier's value with the LINE's
    prediction removes that scatter by construction -- one bad measurement can
    no longer place a mask, it can only tug the fit, and the robust weighting
    stops it doing even that.

    Returns {sec: (cy, cx, R)} predicted from the track, plus a stats dict.
    """
    ss = sorted(tier_moon)
    t = np.array([tier_time[x] for x in ss], float)
    if np.ptp(t) <= 0:
        t = np.arange(len(ss), dtype=float)

    def fit(v):
        w = np.ones_like(v)
        for _ in range(6):
            c = np.polyfit(t, v, 1, w=w)
            r = v - np.polyval(c, t)
            sg = 1.4826 * np.median(np.abs(r - np.median(r))) + 1e-6
            w = 1.0 / np.maximum(np.abs(r) / (2.0 * sg), 1.0)
        return c, v - np.polyval(c, t)

    cy, ry = fit(np.array([tier_moon[x][0] for x in ss]))
    cx, rx = fit(np.array([tier_moon[x][1] for x in ss]))
    cr, _ = fit(np.array([tier_moon[x][2] for x in ss]))
    rate = float(np.hypot(cy[0], cx[0]))
    span = rate * float(np.ptp(t))
    info = {"drift_px_per_s": rate, "drift_px_total": span,
            "scatter_y_px": float(np.std(ry)), "scatter_x_px": float(np.std(rx))}
    progress.log(f"lunar track: {rate:.2f} px/s, {span:.0f} px across the "
                 f"bracket; scatter about the line {np.std(ry):.0f}/"
                 f"{np.std(rx):.0f} px", None)
    info["_line"] = (cy, cx, cr)
    return {x: (float(np.polyval(cy, t[i])), float(np.polyval(cx, t[i])),
                float(np.polyval(cr, t[i]))) for i, x in enumerate(ss)}, info


def _moon_mask_helps(stacks_half, sat_half, secs, cal, abs_shift,
                     track, progress, min_gain=0.15):
    """Decide by measurement whether per-tier lunar masking improves the merge.

    Merges the half-res luminance twice and compares. Masking is accepted only
    if ALL of these hold:

      * the merged limb's 20-80% width improves by `min_gain`;
      * the limb FIT stays healthy -- circle-fit rms no worse than 1.5x. This
        is the criterion that discriminates: on the reference set good masking
        costs 1.3x rms while the broken two-disc case costs 3.5x;
      * no pixel just outside the disc loses every contributing tier.

    The fit-health test is the one that matters and the one 0.8.2 lacked. Bad
    masking does not blur the limb, it breaks the disc into pieces -- which
    reads as a SHARPER edge on any width metric while the circle fit falls
    apart. Sharpness alone rewards exactly the failure it should reject.
    """
    if len(track) < 4:
        return False, {"verdict": "too few usable per-tier limb fits"}
    H, W = stacks_half[secs[0]].shape
    yv = np.arange(H, dtype=np.float32)
    xv = np.arange(W, dtype=np.float32)

    def build(mask):
        acc = np.zeros((H, W), np.float32)
        wsum = np.zeros((H, W), np.float32)
        for s in secs:
            a = stacks_half[s].astype(np.float32) / np.float32(s * cal[s])
            ady, adx = abs_shift[s]
            g = (~sat_half[s]).astype(np.float32)
            a = ndimage.shift(a, (ady, adx), order=1, mode="nearest")
            g = ndimage.shift(g, (ady, adx), order=1, mode="nearest", cval=0)
            w = np.float32(s) * ndimage.gaussian_filter(g, 10)
            if mask and s in track:
                cyt, cxt, Rt = (v / 2.0 for v in track[s])
                dt = np.hypot(yv[:, None] - cyt, xv[None, :] - cxt)
                w = w * np.clip((dt - (Rt - 0.5)) / 2.0, 0.0, 1.0)
            acc += w * a
            wsum += w
        return acc / np.maximum(wsum, 1e-9), wsum

    m0, _ = build(False)
    m1, w1 = build(True)
    cy0 = float(np.median([track[s][0] for s in track])) / 2.0
    cx0 = float(np.median([track[s][1] for s in track])) / 2.0
    # the track measured the radius on every tier; use it rather than a frame
    # fraction, which is a statement about focal length
    _r0 = float(np.median([track[s][2] for s in track])) / 2.0
    if not np.isfinite(_r0) or _r0 <= 0:
        _r0 = 0.10 * min(m0.shape)
    f0 = fit_limb_rays(np.clip(m0, 0, None), cy0, cx0, _r0, decim=1)
    f1 = fit_limb_rays(np.clip(m1, 0, None), cy0, cx0, _r0, decim=1)
    if f0 is None or f1 is None:
        return False, {"verdict": "limb fit failed on one of the trial merges"}
    a = align.limb_transition_width(np.clip(m0, 0, None), f0[0], f0[1], f0[2])
    b = align.limb_transition_width(np.clip(m1, 0, None), f1[0], f1[1], f1[2])
    if not a or not b:
        return False, {"verdict": "limb width not measurable"}
    w_off, w_on = a["limb_width_med"] * 2.0, b["limb_width_med"] * 2.0
    gain = (w_off - w_on) / max(w_off, 1e-6)
    rms_ratio = f1[3] / max(f0[3], 1e-6)
    dR = abs(f1[2] - f0[2]) / max(f0[2], 1e-6)
    r = np.hypot(yv[:, None] - f0[0], xv[None, :] - f0[1])
    ring = (r > f0[2] * 1.02) & (r < f0[2] * 1.6)
    hole = float((w1[ring] <= 1e-6).mean())
    reasons = []
    if gain < min_gain:
        reasons.append("no real gain")
    if rms_ratio > 1.5:
        reasons.append(f"limb fit degraded ({rms_ratio:.1f}x rms)")
    # The fitted radius also moves when masking works -- a sharp edge crosses
    # 50% further in than a 23 px ramp does -- and I have no calibrated bound
    # for how much is too much. It is reported, not gated on, rather than
    # invent a threshold that happens to admit the answer I expect.
    if hole >= 0.001:
        reasons.append(f"{100 * hole:.2f}% of the near-limb corona uncovered")
    ok = not reasons
    info = {"limb_width_off_px": w_off, "limb_width_on_px": w_on, "gain": gain,
            "rms_off": f0[3], "rms_on": f1[3], "radius_change": dR,
            "uncovered_frac": hole, "applied": bool(ok),
            "verdict": "applied" if ok else "rejected — " + "; ".join(reasons)}
    progress.log(
        f"moon-mask trial: limb 20-80% {w_off:.1f}px -> {w_on:.1f}px "
        f"({100 * gain:+.0f}%), fit rms {f0[3]:.2f} -> {f1[3]:.2f}, "
        f"radius {100 * dR:+.1f}% -> {info['verdict']}", None)
    return ok, info


def tier_headroom(wb, cam2rgb):
    """Largest output value a raw pixel at saturation can reach, in sat units.

    White balance runs at G=1, so the red and blue gains push those channels
    ABOVE the raw saturation level, and the camera->sRGB matrix can push
    further still -- typically to 3-4x for a strongly red pixel. Normalising by
    sat_level alone therefore clipped strongly coloured highlights (prominence
    H-alpha above all) that the sensor had NOT clipped.

    The cost of dividing by this instead is that scene-white lands at 1/hi
    rather than 1.0, i.e. the file looks about two stops dark. That is the
    right trade for the scene-linear export, where a host is going to rescale
    anyway and losing recorded data is the only real sin. It is the wrong trade
    for the sRGB export, which has to look correct the moment it opens, so
    that one normalises by sat_level and accepts clipping in colours that are
    already blown in the raw for the tiers where it happens.
    """
    m = np.maximum(np.asarray(cam2rgb, np.float64), 0.0)
    return float(max((m @ np.asarray(wb, np.float64)).max(), 1e-6))


_TIER_README = """\
EclipseForgeHDR -- aligned exposure tiers
=========================================

One 16-bit TIFF per exposure tier. Each file is: hot pixels repaired, that
tier's frames aligned to each other and averaged, demosaiced, white-balanced,
converted to sRGB primaries, and shifted into the common cross-tier frame.
NOT tone mapped, NOT merged, NOT stretched.

Every file carries an embedded ICC profile, so Photoshop, Affinity and
PixInsight all interpret it the same way.
{enc_note}

Scale: pixel value 1.0 (65535) corresponds to {hi:.3f} x sensor saturation for
a scene-neutral pixel. Each tier sits at its own real signal level for its own
exposure time, so a short tier looks dark overall and a long tier is
legitimately blown across the inner corona -- that is what the bracket looked
like through the lens. A histogram against the right wall on the long tiers is
the data, not a bug.

HOW TO USE THESE
----------------
Load them as layers and blend them by hand (masks / luminosity masks), or feed
the whole set to an HDR combine (Merge to HDR Pro, PixInsight HDRComposition).

DO NOT mean- or median-stack the set as if it were one exposure. Averaging a
bracket averages blown frames with unblown ones; the result is dominated by the
longest exposure and clips wherever that one clipped. Mean stacking is the
right operation for several frames of the SAME exposure -- which is already
done inside each of these files.
"""


def _write_tier_tiff(folder, sec, cal_s, sat_level, hi, rgb_norm, n_frames,
                     linear, progress):
    """Write one aligned tier as a 16-bit TIFF for third-party stacking.

    The pipeline divides by (exposure x photometric factor) before this point,
    so that is multiplied back out: the file holds the tier at its own real
    signal level, matching its EXIF exposure time.

    Encoding: sRGB by default, because that is what Photoshop and Affinity
    expect and an untagged linear file shown through an sRGB decode is what
    made these look overstretched with a clipped histogram. `linear=True`
    writes scene-linear instead, for PixInsight and for anyone doing their own
    HDR combine. Either way the matching ICC profile is embedded, so no host
    has to guess.
    """
    import tifffile
    from . import icc
    out = os.path.join(folder, "eclipseforge_output", "aligned_tiers")
    os.makedirs(out, exist_ok=True)
    nm = (f"1_{1 / sec:.0f}s" if sec < 0.4 else f"{sec:g}s").replace(".", "p")
    path = os.path.join(out, f"tier_{nm}_{'linear' if linear else 'srgb'}.tif")
    top = sat_level * (hi if linear else 1.0)
    a = np.clip(rgb_norm * np.float32(sec * cal_s / max(top, 1e-9)), 0.0, 1.0)
    if linear:
        prof, enc = icc.linear_profile(), "scene-linear"
    else:
        a = icc.encode_srgb(a)
        prof, enc = icc.srgb_profile(), "sRGB transfer function"
    hi = hi if linear else 1.0
    tifffile.imwrite(path, (a * 65535.0 + 0.5).astype(np.uint16),
                     photometric="rgb",
                     extratags=[(34675, 1, len(prof), prof, False)],
                     description=(f"EclipseForgeHDR aligned tier; exposure {sec}s; "
                                  f"{n_frames} frames; photometric factor "
                                  f"{cal_s:.4f}; {enc}; sRGB primaries; "
                                  f"1.0 = {hi:.4f} x sensor saturation"))
    with open(os.path.join(out, "README.txt"), "w") as fh:
        fh.write(_TIER_README.format(
            hi=hi,
            enc_note=("These are scene-LINEAR (gamma 1.0). Do not let a host\n"
                      "apply an sRGB decode to them -- that is what makes an\n"
                      "untagged linear file look harsh and clipped."
                      if linear else
                      "These carry the sRGB transfer function, so they open\n"
                      "looking correct with no adjustment.")))
    progress.log(f"  wrote {os.path.basename(path)}", None)


def remove_sky_gradient(wd, cy, cx, R, extent_R, stats, progress):
    """Divide out the smooth brightness AND colour gradient across the frame.

    At low solar altitude the sky is neither uniform nor uniformly coloured. On
    the reference set (sun below 7 deg) it varies by 1.20x corner to corner in
    red and 1.34x in blue -- blue steepest, which is what Rayleigh scattering
    with airmass does, and the best evidence that what is being fitted is real
    atmosphere rather than an artefact.

    WHAT IT IS NOT: lens vignetting. That is radial about the FRAME centre, and
    such a model explains 0.0% of the tilt and 2.5% of the curvature here,
    against 54% and 52% for a plane and a quadratic. Flats would not remove it;
    it has to be fitted per run.

    THE ORDER: a plane leaves a curved remainder that reads as darkening on the
    OTHER side and around the corners. Measured on what a plane leaves behind:
    another plane explains 0.8%, a radial model 2.5%, a full quadratic 51.8%.

    PER CHANNEL, not on luminance. One shared correction flattens brightness and
    leaves the colour gradient untouched. Fitting each channel separately takes
    the colour with it -- measured on the reference set, the spread of sky colour
    between frame quadrants:

                  before   after
        R         0.0551   0.0106
        G         0.0123   0.0024
        B         0.0429   0.0072

    WHERE IT IS FITTED matters more than the model. Fitted close in, it absorbs
    the corona's OWN east-west asymmetry, which is real solar structure: the
    fitted amplitude runs 1.22 in the 1.6-2.4 R shell, 0.55 at 2.4-3.2 R, 0.34
    at 3.2-4 R and 0.22 beyond 4 R, converging only once the corona has faded.
    The direction is stable at -13 to -18 deg throughout, which is what says it
    is one real gradient. So the fit is restricted to beyond the MEASURED corona
    extent. After removal the sky shell keeps 15% of its gradient while the
    corona shells keep 93% and 71% of theirs.
    """
    hp = os.path.join(wd, "hdr_rgb.npy")
    if not os.path.exists(hp):
        return
    hdr = np.load(hp, mmap_mode="r")
    H, W, _ = hdr.shape
    d = 6
    S = np.asarray(hdr[::d, ::d], np.float32)
    h, w, _ = S.shape
    yy = np.arange(h, dtype=np.float32)[:, None] - cy / d
    xx = np.arange(w, dtype=np.float32)[None, :] - cx / d
    r = np.hypot(yy, xx)
    Rd = max(R / d, 1e-6)
    rb = np.clip((r / Rd * 8).astype(np.int32), 0, 400)
    # beyond the corona, with a floor so a short bracket cannot fit on the corona
    r_fit = max(float(extent_R or 0.0), 4.0)
    m = (r > r_fit * Rd)
    if m.sum() < 20000:
        progress.log(f"sky gradient: too little sky beyond {r_fit:.1f} R to fit "
                     f"({int(m.sum())} px) — skipped", None)
        return
    # six terms need plenty of sky to be stable; fall back to a plane if not
    quad = int(m.sum()) >= 100000

    def _design(X, Y):
        cols = [np.ones_like(X), X, Y]
        if quad:
            cols += [X * X, Y * Y, X * Y]
        return np.stack(cols, 1)

    Xg = np.ones((h, 1), np.float32) * (xx / w)
    Yg = (yy * np.ones((1, w), np.float32)) / h
    A = _design(Xg[m], Yg[m])
    coeffs, spans = [], []
    for ch in range(3):
        c_ = S[:, :, ch]
        pos = c_[c_ > 0]
        if pos.size < 1000:
            return
        # TRUE log, not log1p: the correction is multiplicative on the linear
        # image, so the fit has to live where a multiplicative change is an
        # additive one. log1p(S/median) compresses exactly where the sky sits,
        # and a fit made there came out at half strength.
        Ls = np.log(np.maximum(c_, max(float(np.percentile(pos, 0.5)), 1e-12)))
        prof = np.array([np.median(Ls[rb == k]) if (rb == k).sum() > 20 else np.nan
                         for k in range(int(rb.max()) + 1)])
        ok = np.flatnonzero(np.isfinite(prof))
        if ok.size < 4:
            return
        prof = np.interp(np.arange(prof.size), ok, prof[ok])
        res = Ls - prof[rb]
        cc, *_ = np.linalg.lstsq(A, res[m], rcond=None)
        mf = A @ cc
        coeffs.append(cc)
        spans.append(float(mf.max() - mf.min()))
    amp = float(max(spans))
    ang = float(np.degrees(np.arctan2(coeffs[1][2], coeffs[1][1])))
    # significance, by bootstrap on the green channel
    Lg = np.log(np.maximum(S[:, :, 1],
                           max(float(np.percentile(S[:, :, 1][S[:, :, 1] > 0], 0.5)), 1e-12)))
    pg = np.array([np.median(Lg[rb == k]) if (rb == k).sum() > 20 else np.nan
                   for k in range(int(rb.max()) + 1)])
    okg = np.flatnonzero(np.isfinite(pg))
    pg = np.interp(np.arange(pg.size), okg, pg[okg])
    resg = (Lg - pg[rb])[m]
    rs = np.random.default_rng(0)
    aa = []
    for _ in range(24):
        i = rs.integers(0, A.shape[0], A.shape[0] // 4)
        c2, *_ = np.linalg.lstsq(A[i], resg[i], rcond=None)
        mm = A @ c2
        aa.append(float(mm.max() - mm.min()))
    sig = amp / max(float(np.std(aa)), 1e-9)
    stats["sky_gradient"] = {"amp_log": amp, "ratio": float(np.exp(amp)),
                             "angle_deg": ang, "sigma": sig,
                             "fitted_beyond_R": r_fit,
                             "order": 2 if quad else 1,
                             "per_channel": [float(np.exp(x)) for x in spans]}
    # a sky model has no business spanning more than a factor of two
    if amp > 0.7:
        progress.log(f"sky gradient: fitted model spans {np.exp(amp):.2f}x — "
                     f"implausible for sky, not applied", None)
        stats["sky_gradient"]["applied"] = False
        return
    if amp < 0.02 or sig < 8.0:
        progress.log(f"sky gradient: {np.exp(amp):.3f}x across the frame "
                     f"({sig:.0f} sigma) — below the threshold, left alone", None)
        stats["sky_gradient"]["applied"] = False
        return
    # Evaluated in row blocks: six full-resolution term arrays per channel would
    # be several gigabytes on a 45 MP frame.
    yf_all = (np.arange(H, dtype=np.float32) - cy) / float(H)
    xf = ((np.arange(W, dtype=np.float32) - cx) / float(W))[None, :]
    out = np.load(hp)
    for y0 in range(0, H, 512):
        y1 = min(y0 + 512, H)
        yf = yf_all[y0:y1, None]
        for ch in range(3):
            cc = coeffs[ch]
            e = cc[1] * xf + cc[2] * yf
            if quad:
                e = e + cc[3] * (xf * xf) + cc[4] * (yf * yf) + cc[5] * (yf * xf)
            out[y0:y1, :, ch] *= np.exp(-e).astype(np.float32)
    np.save(hp, out)
    lum2 = (0.2126 * out[:, :, 0] + 0.7152 * out[:, :, 1]
            + 0.0722 * out[:, :, 2]).astype(np.float32)
    np.save(os.path.join(wd, "hdr_lum.npy"), lum2)
    stats["sky_gradient"]["applied"] = True
    progress.log(f"sky gradient removed per channel: R {np.exp(spans[0]):.3f}x "
                 f"G {np.exp(spans[1]):.3f}x B {np.exp(spans[2]):.3f}x across the "
                 f"frame ({'quadratic' if quad else 'plane'}, tilt {ang:+.0f} deg), "
                 f"fitted beyond {r_fit:.1f} R ({sig:.0f} sigma)", None)


def shift_bayer_even(a, dy, dx):
    """Shift a Bayer mosaic by an EVEN number of pixels without wrapping.

    Even shifts keep every photosite on its own colour, so the mosaic stays
    decodable. np.roll would do that too, but it wraps: a 300 px shift folds a
    300 px strip of the opposite edge into the frame, which then merges as if it
    were real sky. Edge-replication is wrong too, but it is wrong quietly at the
    border instead of loudly in the middle.
    """
    dy = int(round(dy / 2.0)) * 2
    dx = int(round(dx / 2.0)) * 2
    if dy == 0 and dx == 0:
        return a
    out = np.empty_like(a)
    ys_src = slice(max(0, -dy), a.shape[0] - max(0, dy))
    ys_dst = slice(max(0, dy), a.shape[0] - max(0, -dy))
    xs_src = slice(max(0, -dx), a.shape[1] - max(0, dx))
    xs_dst = slice(max(0, dx), a.shape[1] - max(0, -dx))
    out[...] = a[a.shape[0] // 2, a.shape[1] // 2]
    out[ys_dst, xs_dst] = a[ys_src, xs_src]
    # replicate the nearest valid row/column into the vacated border
    if dy > 0:
        out[:dy] = out[dy:dy + 1]
    elif dy < 0:
        out[dy:] = out[dy - 1:dy]
    if dx > 0:
        out[:, :dx] = out[:, dx:dx + 1]
    elif dx < 0:
        out[:, dx:] = out[:, dx - 1:dx]
    return out


# ---------- main pipeline ----------

def resolve_flat_dir(folder, flat_dir=None):
    """Where the flats are: what the caller asked for, or the convention.

    "" / None  -> a flats subfolder of the light folder if there is one
    "off"      -> no flat correction even if such a folder exists
    anything else -> that path, expanded
    """
    from . import flat as _flat
    s = (flat_dir or "").strip()
    if s.lower() in ("off", "none", "no", "-"):
        return None
    if s:
        return os.path.abspath(os.path.expanduser(s))
    return _flat.find_flat_dir(folder)


def run(folder, progress: Progress, crop_pc=1600, denoise="fine",
        earthshine=False, despeckle=True, frames="all", export_tiers=False,
        tier_linear=False, flat_dir=None):
    from . import flat as _flat
    wd = workdir(folder)
    paths = list_raws(folder)
    if len(paths) < 3:
        raise RuntimeError(f"only {len(paths)} raw files found in {folder}")
    progress.log(f"{len(paths)} raw files found", 0.01)
    _flat_dir = resolve_flat_dir(folder, flat_dir)
    stats = {"version": __version__, "folder": folder, "n_files": len(paths),
             "options": {"denoise": denoise, "earthshine": bool(earthshine),
                         "despeckle": bool(despeckle), "frames": frames,
                         "export_tiers": bool(export_tiers),
                         "tier_linear": bool(tier_linear),
                         "flat_dir": _flat_dir},
             "camera_info": read_camera_info(paths[0])}

    # --- metadata & tiers ---
    meta = {}
    for p in paths:
        sec, iso, ts = read_exif(p)
        meta[p] = {"sec": sec, "iso": iso, "ts": ts}
    tiers = {}
    for p in paths:
        tiers.setdefault(meta[p]["sec"], []).append(p)
    secs = sorted(tiers)
    isos = sorted({meta[p]["iso"] for p in paths if meta[p]["iso"]})
    stats["iso"] = ", ".join(str(i) for i in isos) if isos else None
    tss = sorted(x for x in (meta[p]["ts"] for p in paths) if x)
    if tss:
        stats["shot_first"], stats["shot_last"] = tss[0], tss[-1]
    def _epoch(ts):
        try:
            d, tm = str(ts).split(" ")
            Y, Mo, D = (int(x) for x in d.replace("-", ":").split(":"))
            h, mi, se = (int(float(x)) for x in tm.split(":")[:3])
            return (((Y * 12 + Mo) * 31 + D) * 24 + h) * 3600.0 + mi * 60.0 + se
        except Exception:
            return None
    tier_time = {}
    for s_ in secs:
        es = [e for e in (_epoch(meta[p_]["ts"]) for p_ in tiers[s_]) if e is not None]
        tier_time[s_] = float(np.mean(es)) if es else None
    if any(v is None for v in tier_time.values()):
        # no usable timestamps: order the tiers by cumulative frame count, which
        # is the shooting order and close enough to a time axis
        cum, acc_n = {}, 0
        for s_ in secs:
            cum[s_] = acc_n + len(tiers[s_]) / 2.0
            acc_n += len(tiers[s_])
        tier_time = cum
    else:
        t0 = min(tier_time.values())
        tier_time = {kk: vv - t0 for kk, vv in tier_time.items()}

    progress.log(f"{len(secs)} exposure tiers: " +
                 ", ".join(f"{s:g}s x{len(tiers[s])}" for s in secs), 0.03)

    # --- per-tier: decode, quality, intra-align, stack (half-res + full-res bayer) ---
    _edge = {}          # per tier: invalid border (top, bottom, left, right), full-res px
    sat_half = {}
    stacks_half = {}
    stacks_bayer = {}
    quality = {}
    color_info = None
    n_done = 0
    hot = None
    # Master flat, built once and divided out of every frame of every tier.
    # Built BEFORE any light frame is decoded: the build itself holds four
    # frame-sized accumulators, and overlapping that with a tier's worth of
    # decoded frames would double the peak working set for no reason.
    flat_master = None
    stats["flat"] = {"dir": _flat_dir}
    if _flat_dir:
        try:
            flat_master, _fi = _flat.load_or_build(folder, _flat_dir, None,
                                                   progress, wd)
            stats["flat"] = dict(_fi, dir=_flat_dir)
        except Exception as e:
            progress.log(f"flat correction skipped — the master flat could not "
                         f"be built ({e})", None)
            stats["flat"] = {"dir": _flat_dir, "error": str(e)}
            flat_master = None
        if flat_master is None:
            progress.log("continuing without flat correction", None)
    for s in secs:
        files = tiers[s]
        lums, sats, bayers = {}, {}, {}
        raw_bayers = []
        for p in files:
            rf = open_frame(p)
            if color_info is None:
                color_info = {"wb": rf.daylight_wb.tolist(),
                              "cam2rgb": rf.cam2rgb.tolist(),
                              "sat_level": rf.sat_level,
                              "shape": list(rf.shape)}
            raw_bayers.append((p, rf.bayer, rf.sat_level))
            del rf
        # A flat from another body, another crop mode or another orientation is
        # not a flat for these frames. Checked against the frames of THIS tier,
        # so a bracket that changes size partway disables the flat with a
        # sentence instead of raising a broadcast error at the division.
        _lshape = raw_bayers[0][1].shape if raw_bayers else tuple(color_info["shape"])
        if flat_master is not None and flat_master.shape != tuple(_lshape):
            progress.log(f"flat correction DISABLED — the master flat is "
                         f"{flat_master.shape[1]}x{flat_master.shape[0]} px and "
                         f"the {_exp_name(s)} frames are {_lshape[1]}x"
                         f"{_lshape[0]}; flats have to come from the "
                         f"same camera in the same crop mode", None)
            stats["flat"]["error"] = "flat/light frame size mismatch"
            flat_master = None
        stats["flat"]["applied"] = flat_master is not None
        # sensor defects: map them once, on the shortest tier (darkest sky, so
        # real sky objects cannot be mistaken for hot pixels), then repair every
        # frame of every tier with it
        if hot is None and despeckle:
            # ... but only if this tier can actually show a defect. A hot pixel
            # is a pixel far above its neighbours, and in a saturated frame
            # there is no "above": every photosite sits at the clip. One stray
            # overexposed file whose EXIF puts it at the short end becomes the
            # shortest tier, and the map built from it found ZERO defects on a
            # sensor with 434 -- silently turning off hot-pixel repair for the
            # whole run. Measured on a synthetic bracket: 434 -> 0. So a tier
            # that is mostly clipped is skipped and the next one is tried.
            _sf = float(np.mean([float((b >= sl).mean())
                                 for _, b, sl in raw_bayers[:4]]))
            if _sf > 0.5:
                progress.log(f"sensor defect map: not from the {_exp_name(s)} "
                             f"tier — {100 * _sf:.0f}% of it is saturated, which "
                             f"cannot show a hot pixel; trying the next tier",
                             None)
            else:
                hot = hot_pixel_map([b for _, b, _ in raw_bayers[:4]])
                nhot = int(hot.sum()) if hot is not None else 0
                progress.log(f"sensor defect map: {nhot} hot/dead photosites "
                             f"({100.0 * nhot / max(hot.size, 1):.4f}%)", None)
                if nhot > 0.002 * hot.size:  # implausible -> distrust and disable
                    progress.log("defect count implausibly high — skipping repair",
                                 None)
                    hot = np.zeros_like(hot)
                    nhot = 0
                stats["hot_pixels"] = nhot
        for p, bay, sat_level in raw_bayers:
            if hot is not None:
                repair_hot(bay, hot)
            # The clipping test has to be made BEFORE the flat is divided out.
            # A vignetted corner is brightened by the correction, and a pixel
            # measured against the scalar saturation level afterwards would be
            # declared clipped at a fraction of the well it actually filled.
            h2, w2 = bay.shape[0] // 2, bay.shape[1] // 2
            b = bay[: h2 * 2, : w2 * 2].reshape(h2, 2, w2, 2)
            sats[p] = (b >= sat_level).any(axis=(1, 3))
            del b
            if flat_master is not None:
                bay /= flat_master
            lums[p] = half_luma(bay)
            bayers[p] = bay
        del raw_bayers
        # quality
        scores = {}
        centers = {}
        for p in files:
            cy, cx = find_center(lums[p])
            crop, y0, x0 = crop_around(lums[p], cy, cx, min(1200, min(lums[p].shape)))
            scrop = sats[p][y0:y0 + crop.shape[0], x0:x0 + crop.shape[1]]
            scores[p] = sharpness(crop, scrop)
            centers[p] = (cy, cx)
        # --- drop frames that are not totality ---
        #
        # A partial-phase or diamond-ring frame is a different scene: the
        # photosphere is visible and saturates over a large area, where a
        # totality frame at the same shutter speed saturates only on the
        # chromosphere and prominences, if at all. That crescent is then the
        # brightest thing in the frame, so find_center locks onto it instead of
        # the disc, and every geometry measurement downstream inherits the
        # error -- a lunar radius that came out 857 px instead of 451, a disc
        # mask twice the right size, and a merge containing two scenes.
        #
        # Saturated AREA separates them cleanly and needs no threshold tuned to
        # a camera: it is judged against the other frames of the same tier,
        # which were shot seconds apart at the same exposure.
        _sf = {p: float(sats[p].mean()) for p in files}
        _med = float(np.median(list(_sf.values())))
        _lim = max(3.0 * _med, _med + 0.002)
        _drop = [p for p in files if _sf[p] > _lim]
        if _drop and len(_drop) < len(files):
            for p in _drop:
                progress.log(f"{_exp_name(s)}: dropping {os.path.basename(p)} — "
                             f"{100 * _sf[p]:.2f}% of the frame is saturated vs "
                             f"{100 * _med:.2f}% for this tier; this is not a "
                             f"totality frame", None)
            files = [p for p in files if p not in _drop]
            for p in _drop:
                scores.pop(p, None); centers.pop(p, None)
        elif _drop:
            progress.log(f"WARNING: every frame of the {_exp_name(s)} tier looks "
                         f"like partial phase ({100 * _med:.2f}% saturated); "
                         f"keeping them, but the geometry will be unreliable", None)

        ranked = sorted(scores, key=scores.get, reverse=True)
        best = ranked[0]
        if frames == "best":
            use = ranked[:1]
        elif frames == "best50":
            use = ranked[:max(1, (len(ranked) + 1) // 2)]
        else:
            use = ranked
        sc = [scores[q] for q in ranked]
        quality[s] = {"scores": {os.path.basename(k): v for k, v in scores.items()},
                      "best": os.path.basename(best),
                      "used": [os.path.basename(q) for q in use],
                      "spread": float(sc[0] / max(sc[-1], 1e-9))}
        # intra-tier align to best (half-res), stack half & full-res bayer (even shifts)
        # Crop each frame around ITS OWN disc, not around the best frame's.
        #
        # The window used to be fixed at the best frame's centre and reused for
        # every frame in the tier. A frame that moved between captures then sat
        # off-centre inside it, under a Hanning taper that fades the edges to
        # nothing -- so the correlation weakened with offset and finally failed,
        # silently, blurring the tier average. find_center already knows where
        # each frame's disc is; using it makes the residual offset a few px no
        # matter how far the frame jumped, which is what phase correlation is
        # good at. The coarse part comes from the centres, the fine part from
        # the correlation, and the two are simply added.
        cy, cx = centers[best]
        _sz = min(1200, min(lums[best].shape))
        ref_crop, y0, x0 = crop_around(lums[best], cy, cx, _sz)
        acc_h = lums[best].copy()
        acc_b = bayers[best].copy()
        sat_u = sats[best].copy()
        _intra = []
        _edge.setdefault(s, [0.0, 0.0, 0.0, 0.0])
        # Same rule as the cross-tier window: leave it alone unless the frames
        # genuinely move too far for it. Within a tier the exposure is constant,
        # so find_center is at least CONSISTENT here even though it is not a
        # disc locator -- it drifts with exposure, and there is no exposure
        # change inside a tier.
        _fdev = max((np.hypot(centers[p][0] - cy, centers[p][1] - cx)
                     for p in use), default=0.0)
        _fpc = _fdev > 0.15 * _sz
        for p in use:
            if p == best:
                continue
            if _fpc:
                crop, yp, xp = crop_around(lums[p], centers[p][0], centers[p][1], _sz)
                crop = crop[:ref_crop.shape[0], :ref_crop.shape[1]]
            else:
                crop = lums[p][y0:y0 + ref_crop.shape[0], x0:x0 + ref_crop.shape[1]]
                yp, xp = y0, x0
            (dy, dx), err = pc_shift(ref_crop, crop)
            dy += y0 - yp          # put the coarse re-centring back in
            dx += x0 - xp
            _intra.append(float(np.hypot(dy, dx)) * 2.0)
            acc_h += ndimage.shift(lums[p], (dy, dx), order=3, mode="nearest")
            acc_b += shift_bayer_even(bayers[p], dy * 2, dx * 2)
            _edge[s] = [max(_edge[s][0], max(0.0, dy * 2)),
                        max(_edge[s][1], max(0.0, -dy * 2)),
                        max(_edge[s][2], max(0.0, dx * 2)),
                        max(_edge[s][3], max(0.0, -dx * 2))]
            sat_u |= sats[p]
        if _intra:
            quality[s]["intra_shift_px"] = [round(v, 1) for v in _intra]
            # _intra is full-res (x2 at the point it is appended), _sz is the
            # half-res window, so this used to fire at half the intended motion
            if max(_intra) > 0.5 * _sz:
                progress.log(f"WARNING: {_exp_name(s)} frames move up to "
                             f"{max(_intra):.0f}px within the tier — beyond what "
                             f"the {_sz * 2}px alignment window covers reliably",
                             None)
            elif max(_intra) > 40:
                progress.log(f"{_exp_name(s)}: frames move up to "
                             f"{max(_intra):.0f}px within the tier", None)
        stacks_half[s] = acc_h / len(use)
        stacks_bayer[s] = acc_b / len(use)
        sat_half[s] = sat_u
        n_done += 1
        progress.log(f"tier {s:g}s: best {quality[s]['best']}, {len(use)}/{len(files)} "
                     f"frames stacked (sharpness spread x{quality[s]['spread']:.2f})",
                     0.03 + 0.37 * n_done / len(secs))
        del lums, sats, bayers

    # The shot span and the ISO list were taken from every FILE in the folder,
    # before any frame was rejected -- so one stray file dropped as "not a
    # totality frame" still set them. A test frame shot 20 days after the
    # eclipse made the report say the bracket ran
    # "2026:08:12 21:30:12 .. 2026:09:01 13:35:30": three weeks of totality.
    # Recompute from the frames actually stacked.
    _used = {q for s_ in quality for q in quality[s_]["used"]}
    _um = [p for p in paths if os.path.basename(p) in _used]
    if _um:
        _ts = sorted(x for x in (meta[p]["ts"] for p in _um) if x)
        if _ts:
            stats["shot_first"], stats["shot_last"] = _ts[0], _ts[-1]
        _is = sorted({meta[p]["iso"] for p in _um if meta[p]["iso"]})
        stats["iso"] = ", ".join(str(i) for i in _is) if _is else None

    json.dump({str(k): v for k, v in quality.items()},
              open(os.path.join(wd, "quality.json"), "w"), indent=1)

    # --- cross-tier alignment: lag-1 + lag-2 phase correlation, global LS ---
    progress.log("cross-tier alignment...", 0.42)
    mid = secs[len(secs) // 2]
    # Measure the disc ONCE, here, and use its radius as the seed everywhere
    # below. Every one of those seeds used to be a fraction of the frame, which
    # only works while the Moon is a particular size in the frame.
    _dm = find_disc(stacks_half[mid])
    if _dm is not None:
        cym, cxm, _Rseed = _dm
        progress.log(f"disc found on the {_exp_name(mid)} tier: "
                     f"R={_Rseed * 2:.0f}px full-res", None)
    else:
        cym, cxm = find_center(stacks_half[mid])
        _Rseed = 0.10 * min(stacks_half[mid].shape)
        progress.log("disc edge not detected — falling back to a brightness "
                     "centroid and a frame-fraction radius seed", None)
    S = min(crop_pc, min(stacks_half[mid].shape))
    _, y0, x0 = crop_around(stacks_half[mid], cym, cxm, S)

    # Where to put each tier's correlation window.
    #
    # Normally: nowhere special. One window at the middle tier's position, used
    # for every tier, so the phase correlation sees the full relative offset and
    # measures it to a fraction of a pixel. That is what it is good at.
    #
    # Only when the tiers are so far apart that they would fall outside that
    # window is each one re-centred first. That trade is worth making only at
    # large offsets: the crop origin is an integer, so pre-centring quantises
    # away the sub-pixel accuracy. Measured on the reference set, where the
    # tiers sit 22 px apart, the fixed window gives a 0.58 px network residual
    # and pre-centring gives 1.86 px. At 1000 px apart the fixed window fails
    # outright.
    #
    # The coarse position comes from the LIMB FIT, never from find_center.
    # find_center returns the centroid of the brightest 0.05% of pixels, which
    # is a part of the inner corona and moves with exposure -- 197 px of spread
    # across this bracket, on tiers genuinely 22 px apart. Pre-centring on that
    # injected the error it was meant to remove and took the residual to 20 px.
    _tc = {}
    for _s in secs:
        try:
            _f = fit_limb_rays(stacks_half[_s], cym, cxm,
                               _Rseed, decim=1)
        except Exception:
            _f = None
        _tc[_s] = (_f[0], _f[1]) if (_f is not None and _f[3] < 0.05 * _f[2]) \
            else find_center(stacks_half[_s])
    _dev = max(np.hypot(_tc[x][0] - _tc[mid][0], _tc[x][1] - _tc[mid][1])
               for x in secs)
    _percentre = _dev > 0.15 * S
    _org = {}
    for _s in secs:
        if _percentre:
            _, _oy, _ox = crop_around(stacks_half[_s], _tc[_s][0], _tc[_s][1], S)
        else:
            _oy, _ox = y0, x0
        _org[_s] = (_oy, _ox)
    if _percentre:
        progress.log(f"tiers are spread {2 * _dev:.0f}px apart (full-res) — more "
                     f"than the {2 * S}px correlation window handles, so each is "
                     f"windowed on its own disc", None)
    else:
        progress.log(f"tiers sit within {2 * _dev:.0f}px (full-res); one shared "
                     f"correlation window", None)

    def prep_pair(s1, s2):
        p1y, p1x = _org[s1]
        p2y, p2x = _org[s2]
        a = stacks_half[s1][p1y:p1y + S, p1x:p1x + S] / s1
        b = stacks_half[s2][p2y:p2y + S, p2x:p2x + S] / s2
        sat2 = ndimage.binary_dilation(sat_half[s2][p2y:p2y + S, p2x:p2x + S], iterations=4)
        med = np.median(a); sig = 1.4826 * np.median(np.abs(a - med))
        good = (~sat2) & (ndimage.gaussian_filter(a, 3) > med + 5 * sig)
        wgt = ndimage.gaussian_filter(good.astype(np.float32), 5)
        out = []
        for img in (a, b):
            x = np.log1p(np.clip(img, 0, None) / max(med + 5 * sig, 1e-3))
            x -= ndimage.gaussian_filter(x, 25)
            out.append(x * wgt * np.hanning(S)[:, None] * np.hanning(S)[None, :])
        return out

    # A tier whose correlation window carries no usable signal produces an
    # all-zero prepped image -- prep_pair masks out saturated pixels, and when
    # the whole window is saturated the weight map is zero everywhere. Phase
    # correlation on a zero image returns a finite but meaningless shift and an
    # err of NaN, so the weight 1/(err+0.05) is NaN, and ONE such weight poisons
    # a whole row of the design matrix below: LAPACK then reports "DLASCL
    # parameter number 4 had an illegal value" and numpy raises "SVD did not
    # converge in Linear Least Squares". That was a hard crash on any bracket
    # wide enough to blow its longest tier -- 14.3 EV in the report that found
    # it -- and the message named none of it.
    def _usable(x):
        return np.isfinite(x).all() and float(np.abs(x).max()) > 0

    pairs = []
    _dead = {}
    for lag in (1, 2):
        for i in range(len(secs) - lag):
            _a, _b = secs[i], secs[i + lag]
            p1, p2 = prep_pair(_a, _b)
            if not _usable(p1) or not _usable(p2):
                for _t, _p in ((_a, p1), (_b, p2)):
                    if not _usable(_p):
                        _dead[_t] = _dead.get(_t, 0) + 1
                del p1, p2
                continue
            sh, err, _ = phase_cross_correlation(p1, p2, upsample_factor=20,
                                                 normalization=None)
            if not (np.isfinite(sh).all() and np.isfinite(err)):
                continue
            if _percentre:      # add back where the two windows were taken
                sh = np.array([sh[0] + _org[_a][0] - _org[_b][0],
                               sh[1] + _org[_a][1] - _org[_b][1]], float)
            pairs.append((i, i + lag, sh, 1.0 / (err + 0.05)))
    for _t in sorted(_dead):
        progress.log(f"{_exp_name(_t)}: no usable signal in the correlation "
                     f"window (saturated, or nothing above the noise) — this "
                     f"tier cannot be aligned by correlation", None)
    if _dead:
        stats["align_dead_tiers"] = [_exp_name(t) for t in sorted(_dead)]
    n = len(secs)
    ref = n // 2

    # --- prominence links -------------------------------------------------
    # Phase correlation ties tiers together through the corona, which is what
    # the fastest tiers do not have: at 1/4000 s there is a chromosphere ring,
    # a couple of prominences and almost nothing else, so the correlation has
    # little to lock onto and its answer for those tiers is the least reliable
    # one in the network. The prominences are the only high-SNR features those
    # frames share with the rest of the sequence -- and unlike the lunar limb
    # they are SOLAR, so they do not drift against the corona during totality.
    #
    # These go in as extra links in the same weighted least-squares solve
    # rather than overriding the correlation result, so where both are good
    # they average and where one is bad its weight carries the fact.
    _mtrack_line = None
    prom_links, prom_info = [], {}
    try:
        cov = {s: align.signal_coverage(stacks_half[s], sat_half[s],
                                        (slice(y0, y0 + S), slice(x0, x0 + S)))
               for s in secs}
        # one limb fit, on the mid tier, reused for every tier's prominence
        # search: the Moon moves only a few px across the sequence and the
        # search band is 13% of R wide, so it tolerates that easily.
        # Try several seed radii, and several tiers, before giving up.
        #
        # A single 0.10 * min(shape) seed assumes the Moon fills about a fifth of
        # the short side. That is true of one telescope and one sensor; on a
        # compact at 118 mm the same fraction lands somewhere else entirely, and
        # a bracket whose tiers are all short gives find_center a bright ring
        # rather than a disc to centre on. When this returned None the whole
        # prominence path was skipped with one line of log.
        _f = None
        _sd = _Rseed
        _cand = [mid] + [x for x in secs[::max(1, len(secs) // 4)] if x != mid]
        for _tier in _cand:
            for _frac in (1.0, 0.7, 1.4, 2.0):
                try:
                    _t = fit_limb_rays(stacks_half[_tier], cym, cxm,
                                       _frac * _sd, decim=1)
                except Exception:
                    _t = None
                if _t is not None and _t[3] < 0.05 * _t[2]:
                    _f = _t
                    break
            if _f is not None:
                if _tier != mid:
                    progress.log(f"prominence search: limb found on the "
                                 f"{_exp_name(_tier)} tier (the middle tier "
                                 f"would not fit)", None)
                break
        if _f is None:
            raise RuntimeError(f"no limb fit on any of {len(_cand)} tiers "
                               f"(seed radii 0.7-2.0x {_sd:.0f}px half-res)")
        cyp, cxp, R_half = _f[0], _f[1], _f[2]
        # anchor tier: the one whose prominences are clearest. Scored on the
        # data rather than fixed by index, so it adapts to how a given
        # sequence was exposed.
        progress.log(f"prominence search: limb cy {cyp:.0f} cx {cxp:.0f} "
                     f"R {R_half:.0f} (half-res), rms {_f[3]:.2f}", None)
        cand = []
        counts = []
        for i, s in enumerate(secs):
            a = align.find_prominences(stacks_half[s], cyp, cxp, R_half, n_max=6)
            counts.append(f"{_exp_name(s)}:{len(a)}")
            if len(a) >= 2:
                cand.append((len(a), cov[s], i, a))
        progress.log("prominence anchors per tier — " + " ".join(counts), None)
        if not cand:
            progress.log("no tier yielded 2+ prominence anchors; "
                         "alignment falls back to corona correlation alone", None)
        if cand:
            cand.sort(reverse=True)
            _, _, pi, anchors = cand[0]
            p_ref = secs[pi]
            prom_info = {"tier": float(p_ref), "anchors": len(anchors),
                         "links": [], "used": 0}
            raw = []
            for i, s in enumerate(secs):
                if i == pi:
                    continue
                r = align.align_on_prominences(stacks_half[p_ref],
                                               stacks_half[s], anchors,
                                               R=R_half)
                if r is None:
                    continue
                dy, dx, spread, na_ = r
                prom_info["links"].append(
                    {"sec": float(s), "dy": dy, "dx": dx,
                     "spread": spread, "n": na_})
                # A prominence buried in inner-corona glare gives anchors that
                # disagree with each other; that disagreement is the gate.
                # align_on_prominences(ref, tgt) has the same sense as
                # phase_cross_correlation(ref, tgt), so the link reads
                # shift[tgt] - shift[ref], exactly like the correlation pairs.
                if na_ >= 2 and spread <= PROM_MAX_SPREAD:
                    raw.append((pi, i, (dy, dx), 1.0 / (spread + 0.4)))
            if not raw:
                progress.log("prominence links all rejected (anchors disagreed "
                             f"by more than {PROM_MAX_SPREAD:.1f} px); "
                             "corona correlation alone", None)
            if raw:
                # put the two weight families on a common scale so neither
                # silently dominates just because of how its error is defined
                mp = float(np.median([w for *_, w in raw]))
                mc = float(np.median([w for *_, w in pairs]))
                k = (mc / mp) if mp > 0 else 1.0
                prom_links = [(i, j, sh, w * k * PROM_WEIGHT) for i, j, sh, w in raw]
                prom_info["used"] = len(prom_links)
                progress.log(
                    f"prominence anchors: {len(anchors)} on the "
                    f"{_exp_name(p_ref)} tier, {len(prom_links)}/{n - 1} tiers "
                    f"linked by them", None)
    except Exception as e:
        progress.log(f"prominence linking unavailable ({e}); "
                     f"corona correlation only", None)
    stats["prom_align"] = prom_info

    # Every link is checked before it can enter the solve. One non-finite shift
    # or weight makes the whole least-squares problem unsolvable, and the error
    # LAPACK raises names neither the tier nor the reason.
    _links = [(i, j, sh, w) for i, j, sh, w in pairs + prom_links
              if np.isfinite(w) and np.isfinite(sh).all()]
    _drop = len(pairs) + len(prom_links) - len(_links)
    if _drop:
        progress.log(f"{_drop} alignment link(s) discarded as non-finite", None)
    A, by, bx, ws = [], [], [], []
    for i, j, sh, w in _links:
        row = np.zeros(n); row[j] = 1; row[i] = -1
        A.append(row); by.append(sh[0]); bx.append(sh[1]); ws.append(w)
    row = np.zeros(n); row[ref] = 1
    A.append(row); by.append(0); bx.append(0); ws.append(100.0)
    A = np.array(A); ws = np.array(ws)
    Aw = A * ws[:, None]
    try:
        ay = np.linalg.lstsq(Aw, np.array(by) * ws, rcond=None)[0]
        ax = np.linalg.lstsq(Aw, np.array(bx) * ws, rcond=None)[0]
    except np.linalg.LinAlgError as e:
        progress.log(f"WARNING: the alignment network could not be solved ({e}). "
                     f"Falling back to no cross-tier shift — check the limb in "
                     f"the render before trusting it.", None)
        ay = np.zeros(n); ax = np.zeros(n)
        stats["align_failed"] = str(e)
    ay = np.nan_to_num(ay, nan=0.0, posinf=0.0, neginf=0.0)
    ax = np.nan_to_num(ax, nan=0.0, posinf=0.0, neginf=0.0)
    # A tier nothing could link to is not solved for -- least squares would hand
    # back a minimum-norm number for it, which looks like an answer and is not.
    # It takes the shift of the nearest tier that WAS linked: the Moon moves a
    # few px across a whole bracket, so a neighbour's shift is right to within
    # far less than the error of inventing one.
    _linked = sorted({i for i, j, _, _ in _links} | {j for i, j, _, _ in _links})
    if _linked and len(_linked) < n:
        for i in range(n):
            if i not in _linked:
                k = min(_linked, key=lambda t: abs(t - i))
                ay[i], ax[i] = ay[k], ax[k]
                progress.log(f"{_exp_name(secs[i])}: no usable alignment link — "
                             f"taking the shift measured for {_exp_name(secs[k])}",
                             None)
    elif not _linked:
        progress.log("WARNING: no tier could be linked to any other. The tiers "
                     "are merged unshifted; if the mount moved, the limb will "
                     "be smeared.", None)
    abs_shift = {s: (float(ay[i]), float(ax[i])) for i, s in enumerate(secs)}
    # both axes: a network that is perfect in y and 15px inconsistent in x used
    # to report a residual of 0. A[:-1] is empty when every frame shares one
    # shutter speed (no lag pairs, no prominence links), which made .max() raise.
    if len(A) > 1:
        _ry = np.abs(A[:-1] @ ay - np.array(by[:-1]))
        _rx = np.abs(A[:-1] @ ax - np.array(bx[:-1]))
        res = float(np.maximum(_ry, _rx).max())
    else:
        res = 0.0
    progress.log(f"alignment network residual max {res:.2f}px (half-res)", 0.5)
    stats["align_residual"] = res

    # --- photometric calibration ---
    cal = {mid: 1.0}
    logf = np.zeros(n)
    ratios = []
    for i in range(n - 1):
        a = stacks_half[secs[i]] / secs[i]
        b = stacks_half[secs[i + 1]] / secs[i + 1]
        m = ~ndimage.binary_dilation(sat_half[secs[i + 1]], iterations=6)

        # Compare the two tiers only where BOTH of them actually detected
        # something, judged against each tier's own sky noise.
        #
        # This used to take the brightest 20% of the frame, which is corona only
        # when the corona fills the frame. On a short bracket whose corona
        # reaches 1.6 R in a 18 MP frame the top 20% is mostly sky, so the
        # median ratio was measured on noise, came out 0.13-0.39, and every link
        # was rejected as implausible -- leaving every tier at exactly 1.000 and
        # the merge relying on shutter speeds alone.
        def _floor(x):
            sk = (x[m] if m.any() else x.ravel())[::7]
            md = float(np.median(sk))
            return md + 5.0 * 1.4826 * float(np.median(np.abs(sk - md)))
        good = m & (a > _floor(a)) & (b > _floor(b))
        if good.sum() < 1000:
            progress.log(f"photometric link {_exp_name(secs[i])}->"
                         f"{_exp_name(secs[i + 1])}: only {int(good.sum())} px "
                         f"carry signal in both tiers; using exposure time", None)
            ratios.append(0.0); continue
        lr = float(np.log(np.median(b[good] / np.maximum(a[good], 1e-6))))
        if abs(lr) > np.log(2.5):
            # a tier cannot really be 2.5x off its own exposure time; this means
            # the comparison region is noise, not signal (very short tiers) --
            # trust the shutter speed instead of the measurement
            progress.log(f"photometric link {_exp_name(secs[i])}->"
                         f"{_exp_name(secs[i + 1])} rejected (x{np.exp(lr):.2f} "
                         f"over {int(good.sum())} px) — using exposure time", None)
            lr = 0.0
        ratios.append(lr)
    for i in range(ref, n - 1):
        logf[i + 1] = logf[i] + ratios[i]
    for i in range(ref, 0, -1):
        logf[i - 1] = logf[i] - ratios[i - 1]
    cal = {s: float(np.exp(logf[i])) for i, s in enumerate(secs)}
    stats["quality"] = {str(k): v for k, v in quality.items()}
    stats["tiers"] = [{"sec": float(s), "n": len(quality[s]["used"]),
                       "n_avail": len(tiers[s]),
                       "spread": quality[s]["spread"],
                       "best": quality[s]["best"],
                       "sharpness": float(max(quality[s]["scores"].values())),
                       "cal": float(cal[s]),
                       "shift": [float(abs_shift[s][0]), float(abs_shift[s][1])]}
                      for s in secs]
    progress.log("photometric calibration: " +
                 ", ".join(f"{np.exp(l):.3f}" for l in logf), 0.53)
    span_ev = np.log2(max(secs) / min(secs))
    if span_ev < 6:
        progress.log(f"WARNING: the bracket spans only {span_ev:.1f} EV "
                     f"({min(secs):g}s to {max(secs):g}s). A totality corona "
                     "bracket normally spans 10-14 EV; this looks like a "
                     "partial-phase, diamond-ring or beads sequence, which this "
                     "pipeline is not built for.", None)
    bad = [f"{s:g}s x{cal[s]:.2f}" for s in secs if not 0.25 < cal[s] < 4.0]
    if bad:
        progress.log("WARNING: tiers disagree photometrically far beyond their "
                     "exposure ratio (" + ", ".join(bad) + "). Usually means "
                     "clipped highlights, changing cloud/haze, or exposure "
                     "metadata that does not match the actual exposure. "
                     "The merge will be unreliable.", None)

    # --- per-tier lunar limb, in the ALIGNED frame ---
    # The corona alignment puts the Sun in register; the Moon is a different
    # body and moves against it during the bracket. Blending tiers whose Moon
    # sits elsewhere smears the merged limb over tens of px, so each tier must
    # only contribute where IT sees corona.
    tier_moon = {}
    Rms = []
    for s in secs:
        try:
            f = fit_limb_rays(stacks_half[s], cym, cxm,
                              _Rseed, decim=1)
        except Exception:
            f = None
        if f is None or f[3] > 0.05 * f[2]:
            continue
        ady, adx = abs_shift[s]
        # Each tier keeps its OWN radius, not a shared one. The measured limb
        # radius shrinks with exposure -- about 8 px full-res across a 12 EV
        # bracket on the reference set -- because the 50% crossing between disc
        # and near-limb corona moves inward as the corona brightens. Masking a
        # long tier with a short tier's radius would shave real corona.
        tier_moon[s] = ((f[0] + ady) * 2.0, (f[1] + adx) * 2.0, f[2] * 2.0)
        Rms.append(f[2] * 2.0)
    Rmoon = float(np.median(Rms)) if len(Rms) >= 3 else None
    # Rmoon doubles as the moon-mask enable flag and is cleared below whenever
    # the mask is rejected -- which is the routine outcome. Keep the MEASURED
    # radius separately, or anything that only wants the plate scale silently
    # falls back to the reference set's 620px.
    R_measured = Rmoon
    if len(Rms) >= 3:
        # the tiers' own view of the lunar radius, used later to sanity-check
        # the merged fit
        stats["R_consensus"] = float(np.median(Rms))
        _sp = float(np.percentile(Rms, 90) - np.percentile(Rms, 10))
        if _sp > 0.15 * stats["R_consensus"]:
            progress.log(f"WARNING: the tiers disagree about the lunar radius "
                         f"({min(Rms):.0f}-{max(Rms):.0f}px). They should all see "
                         f"the same Moon; this usually means frames of different "
                         f"scenes are mixed in.", None)
        else:
            progress.log(f"lunar radius consensus across tiers: "
                         f"{stats['R_consensus']:.0f}px (spread {_sp:.0f}px)", None)
    if Rmoon is not None:
        cys = np.array([tier_moon[s][0] for s in tier_moon])
        cxs = np.array([tier_moon[s][1] for s in tier_moon])
        my, mx = float(np.median(cys)), float(np.median(cxs))
        # relative to the disc, not an absolute pixel count: 250 px is 0.4 R
        # on this dataset but 1.25 R on a 200 px moon, where a fit that landed
        # on the wrong side of the disc would survive it
        _rej = 0.25 * (Rmoon if Rmoon else 250.0)
        drop = [s for s in list(tier_moon)
                if np.hypot(tier_moon[s][0] - my, tier_moon[s][1] - mx) > _rej]
        for s in drop:
            del tier_moon[s]
        span = float(np.hypot(cys.max() - cys.min(), cxs.max() - cxs.min()))
        stats["moon_drift_px"] = span
        stats["moon_radius_px"] = Rmoon
        progress.log(f"per-tier lunar limb spread: {span:.0f}px (R_moon {Rmoon:.0f}px)",
                     0.54)
        # Masking each tier to its own disc is the right idea and it is what
        # makes the merged limb an edge instead of a 25 px ramp -- but acting on
        # it unconditionally is what destroyed 0.7.0. When a per-tier fit
        # wanders, its exclusion disc lands on real corona and eats it.
        #
        # So the decision is measured, not assumed. Merge the half-res
        # luminance both ways, measure the limb transition on each, and keep
        # the masking only if it actually sharpens the limb without opening
        # holes. A bad fit now loses the comparison instead of wrecking the run.
        # Per-tier lunar masking, decided by measurement.
        #
        # Mask centres come from a robust straight-line fit of the Moon's
        # position against time, never from the raw per-tier values: the raw
        # scatter is 14-20 px and it clusters, which in 0.8.2 punched two
        # separate discs. The line cannot do that.
        #
        # The gate now also checks that the limb FIT survives, because broken
        # masking reads as a sharper edge on any width metric -- that is how
        # 0.8.2 passed a two-moon merge with a "+26% improvement".
        use_moon_mask = False
        try:
            track, tinfo = _moon_track(tier_moon, tier_time, progress)
            _mtrack_line = tinfo.pop("_line", None)
            stats["moon_track"] = tinfo
            use_moon_mask, mm = _moon_mask_helps(
                stacks_half, sat_half, secs, cal, abs_shift, track, progress)
            stats["moon_mask"] = mm
            if use_moon_mask:
                tier_moon = track
        except Exception as e:
            progress.log(f"moon-mask trial skipped ({e})", None)
        if not use_moon_mask:
            Rmoon = None

    # --- full-res merge ---
    wb = np.asarray(color_info["wb"], np.float32)
    cam2rgb = np.asarray(color_info["cam2rgb"], np.float32)
    sat_level = color_info["sat_level"]
    H2, W2 = stacks_bayer[mid].shape
    _yv = np.arange(H2, dtype=np.float32)
    _xv = np.arange(W2, dtype=np.float32)

    def moon_weight(s):
        """0 inside this tier's own lunar disc, 1 outside, feathered by 4px."""
        if Rmoon is None or s not in tier_moon:
            return None
        cyt, cxt, Rt = tier_moon[s]
        dt = np.sqrt((_yv - cyt)[:, None] ** 2 + (_xv - cxt)[None, :] ** 2)
        return np.clip((dt - (Rt - 1.0)) / 4.0, 0.0, 1.0)

    # --- trim the ragged border left by alignment ---
    #
    # Every shift, within a tier and between tiers, vacates a strip at one edge
    # which is filled by replicating the nearest real row or column. That band
    # is not data. It is a few px when the mount held and hundreds when it did
    # not, and downstream it is worse than useless: the detail filters find
    # enormous gradients in it, and the exported tier TIFFs would carry it into
    # whatever registers them next.
    #
    # The size is not guessed -- every shift applied is known, so the fully
    # valid rectangle is exactly the intersection of what each tier contributed.
    # Computed before the merge so the tier TIFFs get the same rectangle.
    _t = _b = _l = _r = 0.0
    for _s in secs:
        _e = _edge.get(_s, [0.0, 0.0, 0.0, 0.0])
        _cy2, _cx2 = 2 * abs_shift[_s][0], 2 * abs_shift[_s][1]
        _t = max(_t, _e[0] + max(0.0, _cy2))
        _b = max(_b, _e[1] + max(0.0, -_cy2))
        _l = max(_l, _e[2] + max(0.0, _cx2))
        _r = max(_r, _e[3] + max(0.0, -_cx2))
    _t, _b = int(np.ceil(_t)), int(np.ceil(_b))
    _l, _r = int(np.ceil(_l)), int(np.ceil(_r))
    _H, _W = H2, W2
    if (_t + _b) > 0.6 * _H or (_l + _r) > 0.6 * _W:
        progress.log(f"WARNING: alignment would trim {_t + _b}x{_l + _r}px, more "
                     f"than half the frame — not trimming; check the frame "
                     f"motion warnings above", None)
        _t = _b = _l = _r = 0
    if _t or _b or _l or _r:
        _t2, _l2 = _t - _t % 2, _l - _l % 2        # keep the Bayer phase
        crop_origin = (_t2, _l2)
        stats["autocrop_px"] = {"top": _t2, "bottom": int(_b), "left": _l2,
                                "right": int(_r),
                                "kept": f"{_W - _l2 - int(_r)}x{_H - _t2 - int(_b)}"}
        progress.log(f"trimming the alignment border: {_t2}/{_b} top/bottom, "
                     f"{_l2}/{_r} left/right -> "
                     f"{_W - _l2 - int(_r)}x{_H - _t2 - int(_b)}px", None)
    else:
        crop_origin = (0, 0)
    _kh = H2 - crop_origin[0] - int(_b)
    _kw = W2 - crop_origin[1] - int(_r)

    _tier_hi = tier_headroom(wb, cam2rgb)
    # 20px is 3.2% of the reference lunar radius. Left absolute it would be a
    # tenth of a 200px disc, feathering the merge weights across a tenth of the
    # Moon; tied to the measured radius it means the same thing everywhere.
    _feather = float(np.clip(0.032 * (R_measured or 620.0), 8.0, 40.0))
    # Puts a flat-corrected value back into raw units for the clipping test
    # below. 1.0 everywhere when no flat was applied.
    _fsat = None
    if flat_master is not None and flat_master.shape == (H2, W2):
        _fsat = _flat.superpixel_full(flat_master)
    acc = np.zeros((H2, W2, 3), np.float32)
    wsum = np.zeros((H2, W2), np.float32)
    for k, s in enumerate(secs):
        rgb = demosaic_rggb(stacks_bayer[s])
        # cmax is the CLIPPING test, so it has to be measured in raw units,
        # before white balance. WB runs at G=1, so the red gain (~2.1x) used to
        # push a red pixel past 0.97*sat_level while its photosite sat at only
        # 0.45 of saturation: every tier declared prominence H-alpha clipped a
        # full stop early, the shortest one included, and with no tier holding
        # weight there the merge filled prominence cores by leakage from tens
        # of px away. Neutral pixels are unaffected either way -- for a neutral
        # subject green is the largest raw channel and WB leaves it alone --
        # which is why this never showed on the corona itself.
        cmax = rgb.max(axis=2)
        if _fsat is not None:
            cmax *= _fsat            # back to raw units -- see _fsat above
        rgb *= wb[None, None, :]
        rgb = (rgb.reshape(-1, 3) @ cam2rgb.T).reshape(H2, W2, 3)
        rgb /= np.float32(s * cal[s])
        ady, adx = abs_shift[s]
        ady, adx = 2 * ady, 2 * adx
        for c in range(3):
            rgb[:, :, c] = ndimage.shift(rgb[:, :, c], (ady, adx), order=1, mode="nearest")
        knee = 0.87 * sat_level
        wsat = 0.5 * (1.0 + np.tanh((knee - cmax) / (0.06 * sat_level)))
        wsat[cmax > 0.97 * sat_level] = 0.0
        wsat = ndimage.shift(wsat, (ady, adx), order=1, mode="nearest", cval=0)
        w = np.float32(s) * ndimage.gaussian_filter(wsat, _feather)
        mw = moon_weight(s)
        if mw is not None:
            w *= mw
            del mw
        if export_tiers:
            _write_tier_tiff(folder, s, cal[s], sat_level, _tier_hi,
                             rgb[crop_origin[0]:crop_origin[0] + _kh,
                                 crop_origin[1]:crop_origin[1] + _kw],
                             len(tiers[s]), tier_linear, progress)
        acc += w[:, :, None] * rgb
        wsum += w
        del rgb, w, wsat, cmax
        progress.log(f"merged tier {s:g}s", 0.55 + 0.3 * (k + 1) / n)
    hdr = acc / np.maximum(wsum[:, :, None], 1e-9)
    # (_kh, _kw) can differ from (H2, W2) with crop_origin still (0,0) when the
    # trim is bottom/right only. The tier TIFFs were already being sliced in
    # that case while the render was not, so the two came out different sizes
    # and the render kept a band of edge-replicated pixels.
    if crop_origin != (0, 0) or (_kh, _kw) != (H2, W2):
        hdr = hdr[crop_origin[0]:crop_origin[0] + _kh,
                  crop_origin[1]:crop_origin[1] + _kw]
    del acc, wsum
    np.save(os.path.join(wd, "hdr_rgb.npy"), hdr.astype(np.float32))
    lum = (0.2126 * hdr[:, :, 0] + 0.7152 * hdr[:, :, 1] + 0.0722 * hdr[:, :, 2]).astype(np.float32)
    np.save(os.path.join(wd, "hdr_lum.npy"), lum)

    # --- short-exposure inner stack ---
    # Which tiers count as "short" is an EXPOSURE question, not a count. Taking
    # secs[:4] spans 4.3 EV on the 14-tier reference bracket but 6.4 EV on a
    # 6-tier one, and since the weight is proportional to exposure time the
    # 4th tier would then carry ~8x the weight of the 1st -- making this layer
    # a blown inner corona instead of the crisp one it exists to be. Capping at
    # 24x the shortest exposure reproduces exactly secs[:4] on the reference
    # set (its 4th tier is 20x; the 5th is 40x).
    inner_secs = [x for x in secs[:4] if x <= 24.0 * secs[0]] or secs[:1]
    if inner_secs != secs[:4]:
        progress.log(f"inner stack: {len(inner_secs)} tier(s) within 24x the "
                     f"shortest exposure (to {_exp_name(inner_secs[-1])})", None)
    accn = None; accw = None
    for s in inner_secs:
        rgb = demosaic_rggb(stacks_bayer[s])
        cmax = rgb.max(axis=2)            # raw units -- see the merge loop
        if _fsat is not None:
            cmax *= _fsat
        rgb *= wb[None, None, :]
        rgb = (rgb.reshape(-1, 3) @ cam2rgb.T).reshape(H2, W2, 3)
        lt = (0.2126 * rgb[:, :, 0] + 0.7152 * rgb[:, :, 1] + 0.0722 * rgb[:, :, 2])
        del rgb
        ady, adx = abs_shift[s]; ady, adx = 2 * ady, 2 * adx
        lt = ndimage.shift(lt / (s * cal[s]), (ady, adx), order=1, mode="nearest")
        wsat = 0.5 * (1.0 + np.tanh((0.82 * sat_level - cmax) / (0.08 * sat_level)))
        del cmax
        wsat = ndimage.shift(wsat, (ady, adx), order=1, mode="nearest", cval=0)
        w = np.float32(s) * ndimage.gaussian_filter(wsat, 1.2 * _feather)
        mw = moon_weight(s)
        if mw is not None:
            w *= mw
            del mw
        accn = lt * w if accn is None else accn + lt * w
        accw = w if accw is None else accw + w
        del lt, w, wsat
    short_lum = (accn / np.maximum(accw, 1e-12)).astype(np.float32)
    if crop_origin != (0, 0) or short_lum.shape != lum.shape:
        short_lum = short_lum[crop_origin[0]:crop_origin[0] + lum.shape[0],
                              crop_origin[1]:crop_origin[1] + lum.shape[1]]
    np.save(os.path.join(wd, "short_lum.npy"), short_lum)
    _fsat = None                      # last user of it; ~180 MB on a 45 Mpx frame
    # Where the Moon is in THIS stack, taken from the fitted track rather than
    # from a limb fit on the stack itself.
    #
    # The limb fit locates the 50% crossing between the disc and the near-limb
    # corona. On the two fastest tiers there is no corona to cross against, so
    # it locks onto the chromosphere ring and its centre scatters by +-10 px in
    # the RAW frames, before any alignment. Two of those tiers build this stack,
    # which is how a fit on it landed 74 px from the merged limb -- further than
    # the Moon travels in the entire bracket. The straight-line track is fitted
    # across all fourteen tiers, so the unstable ones cannot drag it.
    try:
        if "_mtrack_line" in dir() and _mtrack_line is not None:
            _tw = [tier_time[x] for x in inner_secs]
            _t = float(np.mean(_tw))
            _cy = float(np.polyval(_mtrack_line[0], _t))
            _cx = float(np.polyval(_mtrack_line[1], _t))
            _cr = float(np.polyval(_mtrack_line[2], _t))
            stats["inner_geom"] = {"cy": _cy - crop_origin[0],
                                   "cx": _cx - crop_origin[1], "R": _cr}
            progress.log(f"inner-stack lunar disc from the track: "
                         f"({_cy:.0f},{_cx:.0f}) R={_cr:.0f}px", None)
    except Exception as e:
        progress.log(f"inner-stack track geometry unavailable ({e})", None)
    progress.log("short-exposure inner stack done", 0.90)

    # --- disc center + limb radius: coarse fit on HDR luminance (robust),
    # then band-restricted refinement on the short stack (crisp limb) ---
    _d0 = find_disc(lum)
    if _d0 is not None:
        cy0, cx0, _R0 = _d0
    else:
        cy0, cx0 = find_center(lum)
        _R0 = None
    cyA, cxA, RA = fit_limb(lum, cy0, cx0, R0=_R0)
    # Two seeds, judged on their own merits. The gradient fit is only a seed:
    # never reject the half-level fit for disagreeing with it (that guard threw
    # away the good solution whenever the seed was bad, leaving the disc mask
    # tens of px off).
    Hs = min(lum.shape)
    cands = []
    # Fit on the MERGED luminance. short_lum blends the four shortest tiers, and
    # the Moon moves against the corona between tiers, so its limb sits tens of
    # px away from the merged limb the composite actually displays (measured on
    # this dataset: 71 px in y). Masking must follow the image being masked.
    for img in (lum, short_lum):
        for seed in ((cyA, cxA, RA), (cy0, cx0, _R0 or 0.12 * Hs)):
            try:
                f = fit_limb_rays(img, seed[0], seed[1], seed[2])
            except Exception:
                f = None
            if f is None:
                continue
            cyc, cxc, Rc, rms, nk, nt = f[:6]
            if (0.02 * Hs < Rc < 0.45 * Hs and rms < 0.08 * Rc and nk > 0.5 * nt
                    and 0 < cyc < lum.shape[0] and 0 < cxc < lum.shape[1]):
                cands.append((rms, f))
        if cands:
            break                   # merged-luminance fit succeeded; use it
    if cands:
        cands.sort(key=lambda t: t[0])
        fit = cands[0][1]
        cyf, cxf, R, rms, nk, nt, limb_prof = fit
        progress.log(f"limb half-level fit: centre ({cyf:.1f},{cxf:.1f}) R={R:.1f}px, "
                     f"rms {rms:.2f}px over {nk}/{nt} rays "
                     f"(seed fit said R={RA:.0f})", None)
    else:
        progress.log("WARNING: the half-level limb fit failed on every seed and "
                     "image. Falling back to the gradient fit, which is much less "
                     "reliable — check the disc mask before trusting this render.",
                     None)
        cyf, cxf, R = cyA, cxA, RA
        if not (0.02 * Hs < R < 0.25 * Hs):
            R = 0.10 * Hs
            progress.log(f"gradient fit radius implausible too — using R={R:.0f}px; "
                         "the disc mask will need manual trimming", None)
        rms = 0.01 * R
        fit = None
        limb_prof = np.full(720, R, np.float32)
    # Cross-check the merged fit against what the individual tiers measured.
    #
    # The Moon's apparent radius is the same in every tier to within a few px --
    # it is the same object seconds apart. So the tiers give a consensus that a
    # single bad merged fit cannot outvote. On a run where partial-phase frames
    # had contaminated the merge this came out 857 px against a per-tier
    # consensus of 451, and the disc mask was drawn twice the right size with
    # only a soft warning.
    _rc = stats.get("R_consensus")
    if _rc and R > 0 and abs(R - _rc) / _rc > 0.15:
        # This used to say "Using the tiers' value" unconditionally, and then
        # only actually use it past 30% or with no per-azimuth fit -- so between
        # 15% and 30% the log claimed an override that had not happened, on the
        # one line a user reads when the disc mask comes out wrong. Each branch
        # now says what it did.
        _took = fit is None or abs(R - _rc) / _rc > 0.30
        _head = (f"WARNING: the merged limb fit says R={R:.0f}px but the "
                 f"individual tiers agree on R={_rc:.0f}px "
                 f"({100 * (R / _rc - 1):+.0f}%). ")
        if _took:
            progress.log(_head + "Using the tiers' value — the merge probably "
                         "contains frames of different scenes.", None)
            R = float(_rc)
            limb_prof = np.full(720, R, np.float32)
            rms = 0.02 * R
            fit = None
        else:
            progress.log(_head + "KEEPING the merged fit: it is a per-azimuth "
                         "measurement and the disagreement is under 30%, where "
                         "the tiers' single number is not clearly the better "
                         "one. Check the disc mask on the preview — if it is "
                         "the wrong size, the tiers were right.", None)
    progress.log(f"lunar limb: center ({cyf:.1f},{cxf:.1f}) R={R:.1f}px", 0.91)
    # The real limb is not a circle: lunar relief, seeing, and the moon's drift
    # between tiers make it wander a few px about the fitted circle. Masks keyed
    # to the mean radius therefore leave a crescent of real limb visible on
    # whichever side runs large — the dark bow at the trailing edge. Rmask
    # covers the excursion.
    # The mask follows the measured limb per azimuth rather than a circle, so a
    # tight margin suffices: no crescent of unmasked limb where the true edge
    # runs large, and no chromosphere eaten where it runs small.
    # The disc mask has to cover the whole brightness TRANSITION, not just the
    # scatter of the circle fit.
    #
    # This used to be 0.8 * rms capped at 6 px. rms measures how circular the
    # limb is; it says nothing about how wide the edge is. The Moon moves during
    # the bracket, so the merged limb is a ramp -- 25 px on the reference set,
    # 39 px at the 90th percentile -- and a 4.6 px mask leaves ~20 px of half-lit
    # lunar edge inside the "corona" region. MGN and the inner layer then compute
    # their local mean and sigma over a steep gradient that is not corona at all,
    # which prints as a bright band on one side of the mask edge and a dark one
    # on the other. That is the rim.
    #
    # So measure the ramp and cover it. On a well-registered stack the ramp is a
    # couple of px and this stays as tight as before; it only opens up by as much
    # as the data actually demands.
    ramp = 0.0
    try:
        _lw = align.limb_transition_width(np.clip(lum[::2, ::2], 0, None),
                                          cyf / 2, cxf / 2, R / 2)
        if _lw:
            ramp = 2.0 * float(_lw["limb_width_p90"])
    except Exception as _e:
        progress.log(f"limb ramp not measurable ({_e}); disc mask falls back to "
                     f"the circle-fit rms, which may leave a rim", None)
    if ramp <= 0:
        progress.log("limb ramp not measurable; disc mask falls back to the "
                     "circle-fit rms, which may leave a rim", None)
    margin = float(np.clip(max(0.8 * rms, 0.9 * ramp), 1.5, 0.08 * R))
    Rmask = float(R + margin)
    progress.log(f"merged limb ramp {ramp:.0f}px (p90) -> disc mask margin "
                 f"{margin:.1f}px", None)
    stats["W"], stats["H"] = int(lum.shape[1]), int(lum.shape[0])
    stats["geometry"] = {"cy": float(cyf), "cx": float(cxf), "R": float(R),
                         "Rmask": Rmask, "limb_margin": margin,
               "limb_prof": [float(x) for x in limb_prof], "rms": float(rms),
                         "rays_kept": int(fit[4]) if fit else None,
                         "rays": int(fit[5]) if fit else None}
    stats.update(_report.measure_image(lum, cyf, cxf, R))
    # now that the corona's extent is measured, the sky beyond it can be fitted
    try:
        remove_sky_gradient(wd, cyf, cxf, R, stats.get("corona_extent_R"),
                            stats, progress)
        lum = np.load(os.path.join(wd, "hdr_lum.npy"))
    except Exception as e:
        progress.log(f"sky gradient removal skipped ({e})", None)

    # --- how well did the tiers actually land on each other? ---
    # Photoshop's 'Variance' stack mode, as a number. Aligned tiers disagree
    # only where the signal is genuinely different; misaligned ones disagree in
    # a broad ring around the limb, because that is where the image gradient is
    # steepest and a sub-pixel error shows up as a large brightness difference.
    # Saturated pixels are excluded first -- otherwise this measures clipping in
    # the long tiers rather than registration.
    try:
        al = []
        for s in secs:
            a = np.where(ndimage.binary_dilation(sat_half[s], iterations=2),
                         np.nan, stacks_half[s].astype(np.float32))
            a = a / max(float(cal.get(s, 1.0)), 1e-6) / s
            ady, adx = abs_shift[s]
            mk = np.isfinite(a).astype(np.float32)
            aw = ndimage.shift(np.nan_to_num(a), (ady, adx), order=1,
                               mode="constant", cval=0)
            mw = ndimage.shift(mk, (ady, adx), order=1, mode="constant", cval=0)
            al.append(np.where(mw > 0.9, aw / np.maximum(mw, 1e-6), np.nan))
        # al holds the UNCROPPED half-res tiers, so the centre has to have the
        # autocrop added back; lum on the next line is cropped and does not.
        _, aq = align.stack_variance(al, (cyf + crop_origin[0]) / 2,
                                     (cxf + crop_origin[1]) / 2, R / 2)
        lw = align.limb_transition_width(np.clip(lum[::2, ::2], 0, None),
                                         cyf / 2, cxf / 2, R / 2)
        if lw:
            aq.update(lw)
        # report in full-resolution pixels, which is what the user sees
        for kk in ("rim_width_px", "limb_width_med", "limb_width_p90"):
            if kk in aq:
                aq[kk] = float(aq[kk]) * 2.0
        stats["align_quality"] = aq
        if aq:
            progress.log(
                "alignment quality: limb variance %.3f, rim %.0f px, "
                "merged limb 20-80%% %.1f px"
                % (aq.get("cov_limb", float("nan")),
                   aq.get("rim_width_px", float("nan")),
                   aq.get("limb_width_med", float("nan"))), None)
        del al
    except Exception as e:
        progress.log(f"alignment quality not measured ({e})", None)

    progress.log(f"disc mask radius {Rmask:.1f}px (limb {R:.1f} + {Rmask - R:.1f})", None)
    json.dump({"cy": float(cyf), "cx": float(cxf), "R": float(R),
               "inner_geom": stats.get("inner_geom"),
               "Rmask": Rmask, "limb_margin": margin,
               "limb_prof": [float(x) for x in limb_prof],
               "secs": [float(s) for s in secs],
               "cal": {str(k): v for k, v in cal.items()},
               "abs_shift": {str(k): v for k, v in abs_shift.items()},
               # where the layer grid sits inside the SENSOR frame. The master
               # flat is the one cached product that is neither aligned nor
               # trimmed, so the preview needs this to line it up.
               "crop_origin": [int(crop_origin[0]), int(crop_origin[1])]},
              open(os.path.join(wd, "geometry.json"), "w"), indent=1)
    del short_lum

    # --- earthshine stack: mid-long tiers (~0.25-1.0s) carry the best
    # glare-to-earthshine ratio; the longest tiers are deeper into scatter.
    # Off by default: it needs long tiers with real headroom over the scattered
    # glare, which most totality brackets simply do not have. ---
    earth_tiers = [s for s in secs if 0.25 <= s <= 1.0] if earthshine else []
    if earthshine and not earth_tiers:
        earth_tiers = [secs[-1]]
    accn = None; wtot = 0.0
    for s in earth_tiers:
        rgb = demosaic_rggb(stacks_bayer[s])
        rgb *= wb[None, None, :]
        rgb = (rgb.reshape(-1, 3) @ cam2rgb.T).reshape(H2, W2, 3)
        lt = (0.2126 * rgb[:, :, 0] + 0.7152 * rgb[:, :, 1] + 0.0722 * rgb[:, :, 2])
        del rgb
        ady, adx = abs_shift[s]; ady, adx = 2 * ady, 2 * adx
        lt = ndimage.shift(lt / (s * cal[s]), (ady, adx), order=1, mode="nearest")
        accn = lt * s if accn is None else accn + lt * s
        wtot += s
        del lt
    lp = os.path.join(wd, "long_lum.npy")
    if accn is not None:
        # every other layer is in the cropped frame; this one was not, so the
        # earthshine model (built on the cropped grid) could not even broadcast
        # against it -- an outright crash at 99% on any run that trimmed.
        _ll = (accn / wtot).astype(np.float32)
        _ll = _ll[crop_origin[0]:crop_origin[0] + lum.shape[0],
                  crop_origin[1]:crop_origin[1] + lum.shape[1]]
        np.save(lp, _ll)
        del _ll
        progress.log("long-exposure earthshine stack done", 0.925)
    else:
        if os.path.exists(lp):
            os.remove(lp)
        progress.log("earthshine disabled — skipping long-exposure stack", 0.925)

    # --- COLOR stack for Halpha prominence detection ---
    # ONE tier only: the moon drifts against the corona between tiers, so
    # averaging smears the limb and dilutes the prominence colour. Prefer the
    # tier nearest 1/100 s among those <= 1/100 s (chromosphere unsaturated).
    fast = [s for s in secs if s <= 0.0101]
    prom_tier = max(fast) if fast else min(secs)
    s = prom_tier
    rgb = demosaic_rggb(stacks_bayer[s])
    rgb *= wb[None, None, :]
    rgb = (rgb.reshape(-1, 3) @ cam2rgb.T).reshape(H2, W2, 3)
    rgb /= np.float32(s * cal[s])
    ady, adx = abs_shift[s]; ady, adx = 2 * ady, 2 * adx
    for c in range(3):
        rgb[:, :, c] = ndimage.shift(rgb[:, :, c], (ady, adx), order=1, mode="nearest")
    # crop FIRST, then bin: the gate is upsampled x2 and laid over the cropped
    # composite, so an uncropped prominence stack paints every prominence
    # crop_origin px off the limb. crop_origin is even, so the half-res grid
    # stays in phase, and the limb seed below is in cropped coordinates too.
    rgb = rgb[crop_origin[0]:crop_origin[0] + lum.shape[0],
              crop_origin[1]:crop_origin[1] + lum.shape[1]]
    _ph, _pw = rgb.shape[0] // 2, rgb.shape[1] // 2
    prgb = rgb[:_ph * 2, :_pw * 2].reshape(_ph, 2, _pw, 2, 3).mean(axis=(1, 3))
    del rgb
    # float32: these are photometric rates that routinely exceed the float16
    # maximum (65504) -> silent inf, NaN redness, dead prominence gate.
    np.save(os.path.join(wd, "prom_rgb.npy"), prgb.astype(np.float32))
    plum = (0.2126 * prgb[:, :, 0] + 0.7152 * prgb[:, :, 1]
            + 0.0722 * prgb[:, :, 2]).astype(np.float32)
    del prgb
    # The prominence stack is ONE tier, so its Moon sits where that tier's Moon
    # sits — not where the merged limb is. Fit it separately or the gate's
    # radial window is offset by the lunar motion between tiers.
    pf = fit_limb_rays(plum, cyf / 2, cxf / 2, R / 2, decim=1)
    if pf is not None and 0.8 * R / 2 < pf[2] < 1.2 * R / 2:
        pgeo = {"cy": float(pf[0]), "cx": float(pf[1]), "R": float(pf[2]),
                "prof": [float(x) for x in pf[6]]}
        progress.log(f"prominence-tier limb: centre ({pf[0] * 2:.0f},{pf[1] * 2:.0f}) "
                     f"R={pf[2] * 2:.0f}px — {np.hypot(pf[0] * 2 - cyf, pf[1] * 2 - cxf):.0f}px "
                     f"from the merged limb (lunar motion between tiers)", None)
    else:
        pgeo = None
    del plum
    progress.log(f"prominence colour stack from {prom_tier:g}s tier", 0.93)
    del stacks_bayer, stacks_half

    # --- detail layers ---
    # prom_geom must be on disk BEFORE build_layers reads geometry.json. It was
    # written after, so a FIRST run always fell back to the merged limb for the
    # prominence gate and only a rebuild of the same folder used the prominence
    # tier's own limb: the same folder produced two different gates depending
    # on whether the cache was warm.
    if pgeo is not None:
        gj = json.load(open(os.path.join(wd, "geometry.json")))
        gj["prom_geom"] = pgeo
        json.dump(gj, open(os.path.join(wd, "geometry.json"), "w"), indent=1)
    from . import detail
    lstats = detail.build_layers(wd, progress, denoise=denoise, earthshine=earthshine)
    if isinstance(lstats, dict):
        stats.update(lstats)
    json.dump({"export_tiers": bool(export_tiers),
               "tier_linear": bool(tier_linear),
               "inputs": input_fingerprint(folder),
               # the flats are inputs too: adding, replacing or removing one
               # has to invalidate the cache exactly as a light frame does
               "flat_dir": _flat_dir,
               "flat_inputs": _flat.fingerprint(_flat_dir),
               # whether it was actually APPLIED, not just requested: a run
               # that found the flats unusable leaves the cached master on
               # disk, and a contact frame loaded afterwards must not be
               # corrected with a flat the composite never saw
               "flat_applied": flat_master is not None,
               "denoise": denoise, "earthshine": bool(earthshine),
               "despeckle": bool(despeckle), "frames": frames,
               "build": __version__},
              open(os.path.join(wd, "opts.json"), "w"))
    import datetime
    stats["finished"] = datetime.datetime.now().isoformat(timespec="seconds")
    _tm = _timing_summary(progress)
    if _tm:
        stats["timing"] = _tm
        progress.log("run time " + _fmt_dur(_tm["total_s"]) + " — slowest: "
                     + "; ".join(f"{m} {_fmt_dur(d)}" for d, m in _tm["slowest"]),
                     None)
    txt = _report.write(wd, stats)
    for line in txt.split("\n"):
        progress.log(line, None)
    progress.log("pipeline complete — summary above, also in "
                 ".eclipseforgehdr/report.txt", 1.0)
    progress.done = True


def fit_limb_rays(lum, cy0, cx0, R0, n_ang=720, iters=5, decim=2):
    """Limb fit from the 50% crossing between the occulted disc level and the
    near-limb corona level along each ray, fitted robustly as
    r(theta) = R + dx*cos(theta) + dy*sin(theta).

    Returns (cy, cx, R, rms, kept, total, profile) with the per-azimuth radius
    profile measured ABOUT THE RETURNED CENTRE — the fit runs one extra pass
    after convergence for exactly that reason. Building the profile from the
    pre-update pass instead leaves it carrying the last centre correction as a
    spurious sinusoid, which then drives every mask that uses it.

    The gradient estimator this replaced took the maximum of the RAW radial
    gradient, which peaks inside the bright inner corona rather than at the limb
    and does so further out where the corona is brighter — biasing R outward and
    dragging the centre toward the bright side. A per-ray half-level crossing is
    normalized by that ray's own contrast, so azimuthal brightness cannot move
    it."""
    sm = ndimage.gaussian_filter(lum[::decim, ::decim], 2)
    H, W = sm.shape
    cy, cx, R = cy0 / decim, cx0 / decim, R0 / decim
    ang = np.linspace(0, 2 * np.pi, n_ang, endpoint=False)
    sa, ca = np.sin(ang), np.cos(ang)

    def measure(cy, cx, R, wide):
        lo_f, hi_f = (0.35, 1.9) if wide else (0.6, 1.4)
        rr = np.arange(lo_f * R, hi_f * R, 0.25, dtype=np.float32)
        if len(rr) < 20:
            return None
        ys = cy + rr[None, :] * sa[:, None]
        xs = cx + rr[None, :] * ca[:, None]
        inside = ((ys >= 0) & (ys <= H - 1) & (xs >= 0) & (xs <= W - 1)).all(axis=1)
        P = ndimage.map_coordinates(sm, [np.clip(ys, 0, H - 1).ravel(),
                                         np.clip(xs, 0, W - 1).ravel()],
                                    order=1).reshape(n_ang, -1)
        n_in = max(int(0.08 * len(rr)), 8)
        step = float(rr[1] - rr[0])
        A_, r_, i_ = [], [], []
        for i in range(n_ang):
            if not inside[i]:
                continue
            v = P[i]
            lo = float(np.median(v[:n_in])); hi = float(np.percentile(v, 95))
            if hi <= lo * 1.3:
                continue
            half = 0.5 * (lo + hi)
            j = int(np.argmax(v > half))
            if j <= 0:
                continue
            v0, v1 = float(v[j - 1]), float(v[j])
            f = 0.0 if v1 == v0 else (half - v0) / (v1 - v0)
            A_.append(ang[i]); r_.append(float(rr[j - 1]) + f * step); i_.append(i)
        if len(r_) < n_ang // 6:
            return None
        return np.array(A_), np.array(r_), np.array(i_)

    def solve(a, rad):
        M = np.stack([np.ones_like(a), np.cos(a), np.sin(a)], 1)
        keep = np.ones(len(rad), bool)
        c = res = None
        for _ in range(3):
            c, *_ = np.linalg.lstsq(M[keep], rad[keep], rcond=None)
            res = rad - M @ c
            sg = 1.4826 * np.median(np.abs(res[keep]))
            keep = np.abs(res) < 2.5 * max(sg, 1e-3)
        return c, res, keep

    for it in range(iters):
        m = measure(cy, cx, R, wide=(it == 0))
        if m is None:
            return None
        c, res, keep = solve(m[0], m[1])
        cx += c[1]; cy += c[2]; R = c[0]

    # final pass AT the converged centre: this is the frame the profile lives in
    m = measure(cy, cx, R, wide=False)
    if m is None:
        return None
    a, rad, idx = m
    c, res, keep = solve(a, rad)
    R = c[0] + 0.0                      # radius from this same pass
    prof = np.full(n_ang, R, np.float32)
    model = c[0] + c[1] * np.cos(ang) + c[2] * np.sin(ang)
    prof[:] = model
    for jj, rv, kp in zip(idx, rad, keep):
        if kp:
            prof[int(jj)] = rv
    prof = ndimage.gaussian_filter1d(np.concatenate([prof] * 3), 3.0,
                                     mode="nearest")[n_ang:2 * n_ang]
    # residual centre offset should now be ~0; fold whatever is left into
    # the returned centre and leave the profile alone
    cy += c[2]; cx += c[1]
    rms = float(np.std(res[keep])) * decim
    return (cy * decim, cx * decim, R * decim, rms, int(keep.sum()), len(rad),
            (prof * decim).astype(np.float32))


def fit_limb(lum, cy0, cx0, n_ang=720, R0=None):
    """Gradient-maximum limb fit. Only ever a SEED for fit_limb_rays.

    The search band used to be a fixed fraction of the frame -- R had to land
    between a tenth and a third of the short side. That is a statement about
    focal length, not about eclipses: on a 24 MP APS-C body it means the limb is
    only inside the band beyond about 320 mm, and on 24 MP full-frame beyond
    about 510 mm. A 240 mm shot puts the disc at 7.5% of the short side, below
    the floor, so the limb was never even looked at -- the fit returned R = 1005
    px for a 301 px disc and every per-tier limb measurement failed with it,
    taking prominence anchoring and per-tier lunar masking down as well.

    The band now comes from the measured disc instead, so nothing here depends
    on how big the Moon happens to be in the frame.
    """
    sm = ndimage.gaussian_filter(lum[::2, ::2], 2)
    cy, cx = cy0 / 2, cx0 / 2
    if R0 is None:
        d = find_disc(lum)
        if d is not None:
            cy, cx, R0 = d[0] / 2, d[1] / 2, d[2]
    rmax_est = (R0 / 2.0) if R0 else min(sm.shape) / 3
    rmin, rmax = rmax_est * 0.55, rmax_est * 1.8
    for _ in range(4):
        ang = np.linspace(0, 2 * np.pi, n_ang, endpoint=False)
        rr = np.arange(rmin, rmax, 0.5)
        ys = cy + rr[None, :] * np.sin(ang)[:, None]
        xs = cx + rr[None, :] * np.cos(ang)[:, None]
        ys = np.clip(ys, 0, sm.shape[0] - 1); xs = np.clip(xs, 0, sm.shape[1] - 1)
        prof = ndimage.map_coordinates(sm, [ys, xs], order=1)
        grad = np.gradient(prof, axis=1)
        ridx = np.argmax(grad, axis=1)
        redge = rr[ridx]
        strength = grad[np.arange(n_ang), ridx]
        m = strength > np.percentile(strength, 50)
        ye = cy + redge[m] * np.sin(ang[m]); xe = cx + redge[m] * np.cos(ang[m])
        A = np.c_[2 * xe, 2 * ye, np.ones(m.sum())]
        sol, *_ = np.linalg.lstsq(A, xe ** 2 + ye ** 2, rcond=None)
        cx, cy = sol[0], sol[1]
        # A degenerate fit (the true limb outside the [1/6, 1/3] search band, so
        # argmax picks noise inside the disc) makes the radicand negative, R NaN
        # and the NEXT iteration's np.arange(nan, nan, 0.5) raise -- killing the
        # run at 91% instead of reaching the plausibility check downstream, which
        # is written to handle exactly this. Bail with the seed instead.
        _rad = sol[2] + cx ** 2 + cy ** 2
        if not np.isfinite(_rad) or _rad <= 0:
            # rmax_est is half-res, like everything else in this loop; the
            # normal return scales back up the same way. The centre returned is
            # the measured disc's when there is one -- bailing out of the
            # refinement is no reason to throw away a good centre.
            return cy * 2, cx * 2, float(rmax_est * 2)
        R = np.sqrt(_rad)
        if not (0.02 * min(sm.shape) < R < 0.45 * min(sm.shape)):
            # rmax_est is half-res, like everything else in this loop; the
            # normal return scales back up the same way. The centre returned is
            # the measured disc's when there is one -- bailing out of the
            # refinement is no reason to throw away a good centre.
            return cy * 2, cx * 2, float(rmax_est * 2)
        rmin, rmax = R * 0.85, R * 1.15
    return cy * 2, cx * 2, R * 2


def prepare_contact(folder, raw_path, progress):
    """Decode a diamond-ring / Baily's-beads frame, align its lunar disc to the
    composite, tonemap, and cache as a screen-blendable layer."""
    from .raw import open_frame, demosaic_rggb
    wd = workdir(folder)
    geo = json.load(open(os.path.join(wd, "geometry.json")))
    progress.log(f"decoding {os.path.basename(raw_path)}...", 0.1)
    rf = open_frame(raw_path)
    # the contact frame goes through the same optics, so it gets the same flat
    # -- but only if the run that built the composite actually applied one
    _mf = os.path.join(wd, "masterflat.npy")
    _op = os.path.join(wd, "opts.json")
    _used_flat = False
    try:
        _used_flat = bool(json.load(open(_op)).get("flat_applied", False))
    except Exception:
        pass
    if _used_flat and os.path.exists(_mf):
        try:
            _m = np.load(_mf)
            if _m.shape == rf.bayer.shape:
                rf.bayer /= _m
                progress.log("contact frame: flat correction applied", None)
            del _m
        except Exception as e:
            progress.log(f"contact frame: flat not applied ({e})", None)
    rgb = demosaic_rggb(rf.bayer)
    rgb *= rf.daylight_wb[None, None, :]
    rgb = (rgb.reshape(-1, 3) @ rf.cam2rgb.T).reshape(rgb.shape)
    H, W = rgb.shape[:2]
    lum = 0.2126 * rgb[:, :, 0] + 0.7152 * rgb[:, :, 1] + 0.0722 * rgb[:, :, 2]
    progress.log("fitting lunar limb on contact frame...", 0.5)
    cy, cx, R = fit_limb(np.clip(lum, 0, None), geo["cy"], geo["cx"])
    fit = fit_limb_rays(np.clip(lum, 0, None), cy, cx, R)
    if fit is not None and 0.5 * R < fit[2] < 1.5 * R:
        cy, cx, R = fit[0], fit[1], fit[2]
        progress.log(f"contact-frame limb (half-level): R={R:.0f}px, rms {fit[3]:.1f}px", None)
    dy, dx = geo["cy"] - cy, geo["cx"] - cx
    sc = geo["R"] / R if 0.9 < geo["R"] / R < 1.1 else 1.0
    progress.log(f"aligning (shift {dy:+.1f},{dx:+.1f}px, R={R:.0f} vs {geo['R']:.0f}, "
                 f"auto-scale x{sc:.4f})", 0.7)
    for c in range(3):
        rgb[:, :, c] = ndimage.shift(rgb[:, :, c], (dy, dx), order=1, mode="nearest")
    if abs(sc - 1) > 1e-4:
        # scale about the composite disc center so the limbs match
        mat = np.array([[1 / sc, 0], [0, 1 / sc]], np.float64)
        off = [geo["cy"] * (1 - 1 / sc), geo["cx"] * (1 - 1 / sc)]
        for c in range(3):
            rgb[:, :, c] = ndimage.affine_transform(rgb[:, :, c], mat, offset=off,
                                                    order=1, mode="constant", cval=0)
    top = np.percentile(lum, 99.9)
    disp = np.clip(rgb / max(top, 1e-6), 0, 1) ** (1 / 2.2)
    # match the composite frame size if the sensor crop differs slightly
    Hc, Wc = np.load(os.path.join(wd, "hdr_lum.npy"), mmap_mode="r").shape
    out = np.zeros((Hc, Wc, 3), np.float32)
    out[:min(H, Hc), :min(W, Wc)] = disp[:min(H, Hc), :min(W, Wc)]
    np.save(os.path.join(wd, "contact_rgb.npy"), out.astype(np.float16))
    progress.log("contact frame ready", 1.0)


def fit_limb_band(lum, cy0, cx0, R0, band=36, n_ang=720):
    """Refine a limb fit within +/-band px of a known solution (full-res units)."""
    sm = ndimage.gaussian_filter(lum[::2, ::2], 2)
    cy, cx, R = cy0 / 2, cx0 / 2, R0 / 2
    b = band / 2
    for _ in range(3):
        ang = np.linspace(0, 2 * np.pi, n_ang, endpoint=False)
        rr = np.arange(R - b, R + b, 0.5)
        ys = np.clip(cy + rr[None, :] * np.sin(ang)[:, None], 0, sm.shape[0] - 1)
        xs = np.clip(cx + rr[None, :] * np.cos(ang)[:, None], 0, sm.shape[1] - 1)
        prof = ndimage.map_coordinates(sm, [ys, xs], order=1)
        grad = np.gradient(prof, axis=1)
        ridx = np.argmax(grad, axis=1)
        redge = rr[ridx]
        strength = grad[np.arange(n_ang), ridx]
        m = strength > np.percentile(strength, 50)
        ye = cy + redge[m] * np.sin(ang[m]); xe = cx + redge[m] * np.cos(ang[m])
        A = np.c_[2 * xe, 2 * ye, np.ones(m.sum())]
        sol, *_ = np.linalg.lstsq(A, xe ** 2 + ye ** 2, rcond=None)
        cx, cy = sol[0], sol[1]
        R = np.sqrt(sol[2] + cx ** 2 + cy ** 2)
    return cy * 2, cx * 2, R * 2
