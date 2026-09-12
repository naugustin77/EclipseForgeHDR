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
    # A broken or missing ExposureTime reaches here as 0 and used to raise
    # ZeroDivisionError from inside a log line (0.22.78).
    if not s or not np.isfinite(s) or s <= 0:
        return "unknown exposure"
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
                  demosaic as _demosaic,
                  cfa_clip_max, hot_pixel_map, repair_hot, read_camera_info)
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


def load_big(path):
    """np.load, memory-mapped when the filesystem allows it.

    mmap keeps a 45 Mpx three-channel float32 out of RAM, and a plain read is
    the fallback when a filesystem will not map -- network drives and some
    exotic mounts genuinely refuse.

    WHAT THIS WAS ORIGINALLY WRITTEN FOR, AND WHY THAT WAS WRONG. A Windows
    tester reported

        sky gradient removal skipped ([Errno 22] Invalid argument:
        'C:\\Users\\...\\OneDrive\\Desktop\\...\\hdr_rgb.npy')

    and the OneDrive in that path was taken as the cause -- cloud-synced folders
    are known to refuse mmap. Two things later showed that was a coincidence.
    The folder was named that way but was an ordinary local one, not synced at
    all. And this fallback was already shipping when a second Windows run
    produced the same Errno 22 on a path with no OneDrive in it, which on its
    own rules out the load: if opening the map were what failed, the except
    below would have caught it.

    The real fault was at the other end of the same file -- `np.save` over a
    mapping still held open, which Windows refuses and POSIX allows. It is
    fixed where it happens, in `remove_sky_gradient` and
    `neutralise_corona_colour`. This function is kept because the fallback is
    cheap insurance for filesystems that really cannot map, but it never fixed
    the bug it was written for, and one folder's name is not evidence about the
    folder.
    """
    try:
        return np.load(path, mmap_mode="r")
    except (OSError, ValueError):
        return np.load(path)


def workdir(folder):
    from .brand import WORKDIR_NAME
    d = os.path.join(folder, WORKDIR_NAME)
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

    # ...but relative contrast has unbounded variance as the signal goes to
    # zero, and `lo` is the 1st percentile of the WHOLE FRAME. When the disc is
    # darker than the sky and covers more than 1% of the picture -- true of any
    # imported stack whose sky sits well above black -- p1 lands INSIDE the
    # disc, so the flooring never flattens it and the shadow's own noise
    # becomes the strongest edge in the image. Measured on the first imported
    # stack that failed: the disc interior sat 23 units above the floor against
    # 2892 at the limb, yet produced 4548 strong-gradient pixels to the limb's
    # 2092, and the circle fit collapsed to R=37px on a 500px disc.
    #
    # So require SIGNAL as well as contrast: an edge means nothing where there
    # is no light on either side of it. This leaves a bright limb untouched
    # (weight 0.97 there) and takes the dark shadow's noise down fivefold.
    g *= (np.maximum(s - lo, 0.0) /
          (np.maximum(s - lo, 0.0) + 0.002 * span)).astype(np.float32)
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


def _clipped_mean(frames, kappa):
    """Mean of aligned frames with a per-pixel outlier rejected.

    WHY NOT THE SAMPLE SIGMA, which is what kappa-sigma normally uses. With the
    frame counts an eclipse tier actually has, it cannot work: an outlier
    inflates the very sigma it is tested against. For three frames a, a, M the
    ratio |M - median| / sd is at most sqrt(9/2) = 2.12 NO MATTER HOW LARGE M
    IS, so a kappa of 3 rejects nothing, ever, and the option silently does
    nothing. That is how this was first written and the smoke test caught it.

    So the scale is a NOISE MODEL, not a sample statistic. Per tier, fit

        sigma^2(signal) = a * signal + b

    -- photon noise proportional to signal, read noise constant -- by taking a
    robust spread of the frame-to-median deviation in signal bins. A particle
    strike is a handful of pixels in millions, so it cannot move a binned
    median; the fit describes the clean population and the strike then sits far
    outside it. This is the same construction `raw.py::_outlier_flags` already
    uses to find hot photosites, for the same reason.

    Returns (mean, n_rejected).
    """
    st = np.stack([np.asarray(f, np.float32) for f in frames], axis=0)
    med = np.median(st, axis=0)
    dev = st - med[None, :, :]
    # robust sigma per signal bin, from the deviations themselves
    lo = float(np.percentile(med, 1.0))
    hi = float(np.percentile(med, 99.9))
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        return (st.mean(axis=0).astype(np.float32), 0)
    nb = 24
    edges = np.linspace(lo, hi, nb + 1)
    idx = np.clip(np.digitize(med, edges) - 1, 0, nb - 1)
    xs, ys = [], []
    flat_i = idx.ravel()
    flat_d = np.abs(dev).reshape(len(frames), -1)
    for b in range(nb):
        sel = flat_i == b
        n = int(sel.sum())
        if n < 64:
            continue
        # median |deviation| over every frame in this bin -> sigma
        v = flat_d[:, sel]
        m = float(np.median(v))
        if not np.isfinite(m) or m <= 0:
            continue
        xs.append(0.5 * (edges[b] + edges[b + 1]))
        ys.append((1.4826 * m) ** 2)
    del flat_d, flat_i
    if len(xs) >= 3:
        A = np.stack([np.asarray(xs, float), np.ones(len(xs))], axis=1)
        coef, *_ = np.linalg.lstsq(A, np.asarray(ys, float), rcond=None)
        a_, b_ = float(coef[0]), float(coef[1])
    else:
        a_, b_ = 0.0, float(np.median(ys)) if ys else 1.0
    var = np.maximum(a_, 0.0) * np.maximum(med, 0.0) + max(b_, 1e-6)
    # the deviation of a frame from the median of N frames is not the deviation
    # from the truth: the median itself carries noise, and with few frames that
    # is not negligible. Inflate the threshold by the median's own share.
    var = var * (1.0 + 1.0 / max(len(frames) - 1, 1))
    # SELF-CALIBRATE THE CUT, because the sigma above is biased low and the size
    # of the bias depends on the frame count. The deviations are taken from the
    # MEDIAN of the frames, and for an odd count one frame IS the median, so a
    # third of the deviations at n = 3 are exactly zero. A median of |dev| over
    # that mixture under-estimates sigma, the threshold comes out too tight, and
    # the first version of this rejected 7.8% of frame-pixels at n = 3 against
    # the 0.27% a 3 sigma cut should reject -- clipping ordinary noise, which
    # biases the mean it is supposed to protect.
    #
    # Rather than derive the correction per frame count, take it from the data:
    # form t = |dev| / sigma_model and put the cut at whichever is larger,
    # `kappa` or the quantile of t that leaves the intended fraction outside.
    # Any constant error in sigma cancels, and a particle strike -- a handful of
    # pixels in millions, far out in the tail -- is nowhere near that quantile,
    # so it still gets cut. If the frames are clean the quantile lands at kappa
    # and nothing changes.
    _sig = np.sqrt(var)
    t = np.abs(dev) / np.maximum(_sig, 1e-6)[None, :, :]
    # two-sided Gaussian tail beyond kappa
    from math import erfc, sqrt as _sqrt
    _p = erfc(float(kappa) / _sqrt(2.0))
    _q = float(np.quantile(t, max(0.0, 1.0 - _p)))
    _cut = max(float(kappa), _q)
    keep = t <= _cut
    _clipped_mean.last_cut = (_cut, _q, _p)
    del dev, var, t, _sig
    cnt = keep.sum(axis=0)
    keep |= (cnt == 0)[None, :, :]      # never reject every frame at a pixel
    cnt = np.maximum(keep.sum(axis=0), 1)
    out = (np.where(keep, st, 0.0).sum(axis=0) / cnt).astype(np.float32)
    nrej = int(keep.size - int(keep.sum()))
    del st, med, keep, cnt
    return out, nrej


#: Kappa for the per-tier kappa-sigma combine. Loose on purpose: with two or
#: three frames the sample sigma is a poor estimate and a tight kappa starts
#: rejecting signal. The job is a particle strike, not trimming a distribution.
STACK_KAPPA = 3.0

#: Different window shapes on the two images of a correlation pair.
#:
#: OFF, against the literature, on measurement. Druckmullerova's master's thesis
#: sec. 4.2.2 p.57 recommends it -- an identical window is a feature both frames
#: carry in the same place, so it correlates with itself and pulls the peak
#: toward zero shift. The reasoning is sound and the effect is not measurable
#: here: rms error 1.524 -> 1.523 px and the zero-shift pull 0.200 -> 0.208 px,
#: i.e. nothing, or very slightly worse. A Hann window is smooth and makes no
#: sharp structure to correlate on, and on this sensor the zero-pull is
#: dominated by the dust and PRNU in the flat, not by the window. Kept as a
#: switch for data where the window really is the sharpest shared feature.
ASYM_WINDOW = False

#: p and q for the semi-phase correlation, as a fraction of max|F|.
#:
#: Druckmullerova's own value, from the real-data precision test in her master's
#: thesis sec. 4.7.2 that agrees with an independent lunar-centroid reference to
#: 0.40 px: "p = q = 0.01 % max{A_f1, A_f2}".
#:
#: NOT the optimum on our own bench, and taken anyway. Sweeping pq on tiers
#: built from the reference set's own corona and his own master flat (rms error on a known
#: injected shift, px):
#:
#:     pq      1e6*    1e-1   1e-2   1e-3   3e-4   1e-4   1e-5    0
#:     rms     1.414   1.443  1.355  1.248  1.193  1.187  1.174  1.174
#:     zero    0.119   0.089  0.080  0.055  0.041  0.041  0.041  0.041
#:        * pq huge is plain cross-correlation
#:
#: Monotone toward full whitening, so pq = 0 wins there by 1%. It is not taken,
#: because that fixture CANNOT represent the case 0.22.27 measured on real data,
#: where whitening lost by 66%: every tier in the fixture is built from the same
#: merged corona, so even the shortest carries full structure, while the real
#: 1/500 s tier carries almost none and full whitening hands it the same weight
#: as a well-exposed one. pq = 1e-4 sits on the plateau -- 1% off the fixture
#: optimum -- and unlike pq = 0 it has a noise floor, which is the protection
#: the real failure needs. Cheap insurance at the price of one percent.
SEMI_PQ = 1e-4


def _semi_phase_shift(a, b, pq=SEMI_PQ, upsample=20):
    """Shift between two images by SEMI-phase correlation.

    Master's thesis Def. 3.25-3.26:

        Z^{p,q} = F1 conj(F2) / ( (|F1| + p) (|F2| + q) )

    p and q are ADDITIVE constants, not exponents. That distinction is the whole
    value of the method: an exponent is scale-free, so a frequency bin carrying
    only noise gets whitened exactly as hard as one carrying the corona, while
    an additive p is a THRESHOLD in the units of the amplitude spectrum -- bins
    well above it whiten almost fully, bins below it barely at all. It is a
    noise floor on the whitening, which is precisely what our 0.22.27 test found
    missing when full whitening lost by 66%: whitening gave a nearly signal-free
    short tier the same weight as a well-exposed one.

    pq = 0 is full phase correlation; pq large tends to plain cross-correlation.

    Returns (shift, err) with the same sign convention and error scale as
    skimage's phase_cross_correlation, so it drops into the same network.
    """
    F1 = np.fft.fft2(np.asarray(a, np.float64))
    F2 = np.fft.fft2(np.asarray(b, np.float64))
    A1 = np.abs(F1)
    A2 = np.abs(F2)
    p = float(pq) * float(A1.max())
    q = float(pq) * float(A2.max())
    R = (F1 * np.conj(F2)) / ((A1 + p) * (A2 + q) + 1e-30)
    c = np.fft.ifft2(R).real
    del F1, F2, A1, A2, R
    j, i = np.unravel_index(int(np.argmax(c)), c.shape)
    pk = float(c[j, i])
    # parabolic sub-pixel on the wrapped neighbourhood
    def _sub(v0, v1, v2):
        den = v0 - 2 * v1 + v2
        return 0.5 * (v0 - v2) / den if abs(den) > 1e-30 else 0.0
    H_, W_ = c.shape
    dy = _sub(c[(j - 1) % H_, i], pk, c[(j + 1) % H_, i])
    dx = _sub(c[j, (i - 1) % W_], pk, c[j, (i + 1) % W_])
    sy = (j + dy + H_ // 2) % H_ - H_ // 2
    sx = (i + dx + W_ // 2) % W_ - W_ // 2
    # err in the spirit of skimage's: 1 - normalised peak, floored away from 0
    nrm = float(np.sqrt(np.mean(c * c))) * c.size ** 0.0
    err = float(np.clip(1.0 - pk / max(abs(pk) + 8.0 * nrm, 1e-30), 1e-3, 1.0))
    del c
    return np.array([sy, sx], float), err


def sharpness(crop, satmask):
    valid = ~ndimage.binary_dilation(satmask, iterations=3)
    hf = crop - ndimage.gaussian_filter(crop, 2.0)
    sig = ndimage.gaussian_filter(crop, 4.0)
    lvl = np.percentile(sig, 90)
    m = valid & (sig > 0.25 * lvl)
    if m.sum() < 1000:
        m = valid
    # A fully saturated crop leaves `valid` empty, and np.mean of an empty
    # slice is NaN (0.22.78). Every frame in the tier then scored NaN, and
    # "best" and "best half" selected by an ordering over NaN -- which is
    # arbitrary, and looked like a working choice. Zero is the honest score for
    # a frame with nothing measurable in it.
    if not m.any():
        return 0.0
    _v = float(np.mean(hf[m] ** 2)) / max(float(np.mean(sig[m])), 1e-6) ** 2 * 1e4
    return _v if np.isfinite(_v) else 0.0


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
    # SCATTER THAT DWARFS THE TRACK MEANS THERE IS NO TRACK.
    #
    # a second tester's 560mm run: "0.63 px/s, 6 px across the bracket; scatter
    # about the line 26/39 px". The Moon moved six pixels and the per-tier
    # measurements are forty pixels off a straight line -- the fit is describing
    # the alignment error, not the Moon. Printed, and nothing said so.
    #
    # Scaled against the Moon's own radius, which _moon_track has in `cr`: the
    # 560mm scatter is 7.4% of R, where the other three real datasets sit at
    # 0.3%, 0.0% and 0.7%.
    _Rt = float(np.polyval(cr, float(np.median(t)))) if len(cr) else 0.0
    _sc = float(max(np.std(ry), np.std(rx)))
    if _Rt > 0 and _sc > max(0.02 * _Rt, 3.0):
        progress.log(
            f"WARNING: the per-tier lunar positions scatter {_sc:.0f}px about "
            f"the fitted track, which is {100 * _sc / _Rt:.0f}% of the lunar "
            f"radius and {'far more than' if _sc > span else 'comparable to'} "
            f"the {span:.0f}px the Moon actually moved. The track is then "
            f"fitting the cross-tier alignment error rather than the Moon, and "
            f"everything keyed to it -- the inner-stack disc above all -- is "
            f"placed by it.", None)
        info["track_scatter_bad"] = round(_sc, 1)
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


def _fine_structure(img, cy, cx, R, r0, r1):
    """Radially-coherent azimuthal fine structure in a shell, as (amp, coh).

    The point of splitting the two is that grain and corona both raise the
    amplitude of a high-pass, and only one of them is worth having. Real
    coronal structure is radially CONTINUOUS -- a streamer at 1.4 R is still
    there at 1.41 R -- while photon noise is independent from one radial
    sample to the next. So the correlation between adjacent radii separates
    them, and amp * coh is the part that is signal.

    Measured on the reference bracket, this is what says the four shortest
    tiers hold real detail at the limb and nothing but noise further out:

        shell          source        amp     coh     amp*coh
        1.01-1.10 R    merged     0.0326   0.995      0.0325
        1.01-1.10 R    short      0.2101   0.993      0.2086
        1.30-1.80 R    merged     0.0212   0.976      0.0207
        1.30-1.80 R    short      0.0576   0.186      0.0107
        1.80-2.60 R    merged     0.0135   0.874      0.0118
        1.80-2.60 R    short      0.1560   0.013      0.0020
    """
    # SAMPLE AT ONE PIXEL, NOT AT A FIXED COUNT.
    #
    # The coherence term only means anything if adjacent radial samples are
    # independent in the data. A fixed nr over a narrow shell samples far finer
    # than a pixel on a small disc -- 0.125 px at R=60 -- and then bilinear
    # interpolation makes even white noise look perfectly correlated, so a
    # noise-dominated merge scores as if it were full of structure. Caught by
    # the trial's own test case: with the short tiers replaced by pure noise it
    # reported "+174% at the limb" and would have tilted the merge toward them.
    nr = int(np.clip(round((r1 - r0) * R), 8, 256))
    na = int(np.clip(round(2 * np.pi * r1 * R), 256, 4096))
    rr = np.linspace(r0 * R, r1 * R, nr)
    th = np.linspace(0, 2 * np.pi, na, endpoint=False)
    ys = cy + rr[None, :] * np.sin(th)[:, None]
    xs = cx + rr[None, :] * np.cos(th)[:, None]
    P = ndimage.map_coordinates(np.asarray(img, np.float32),
                                [ys.ravel(), xs.ravel()], order=1).reshape(na, nr)
    # NORMALISE BY THE RIGHT KIND OF SCALE FOR THE LAYER IN HAND.
    #
    # This divided every radial column by its azimuthal median, full stop. That
    # is correct for a MULTIPLICATIVE quantity -- luminance, or a 0..1 detail
    # layer whose median sits near 0.5 -- where dividing removes the radial
    # falloff and leaves relative contrast.
    #
    # It is catastrophic for a ZERO-CENTRED one. FNRGF returns (L-mu)/sd, whose
    # azimuthal median is ~0 by construction, so `max(median, 1e-9)` divided by
    # roughly 1e-9 and the metric reported an "amplitude" of 2.6e8. Nothing
    # warned; the number was simply wrong, and it was wrong in the direction of
    # looking like an enormous result.
    #
    # So pick by what the data is. A column whose median dominates its own
    # spread is multiplicative and gets divided; anything else is additive and
    # gets the median SUBTRACTED and is scaled by the shell's robust spread.
    # Within either branch the numbers are comparable, which is all a
    # before/after comparison needs.
    _med = np.median(P, axis=0, keepdims=True)
    _mad = 1.4826 * np.median(np.abs(P - _med), axis=0, keepdims=True)
    _scale = max(float(np.median(_mad)), 1e-12)
    if float(np.median(_med)) > 4.0 * _scale:
        P = P / np.maximum(_med, 1e-9)
    else:
        P = (P - _med) / _scale
    # 2 degrees, expressed in samples, so the scale measured does not change
    # with the disc size either
    _sig = max(2.0, na * 2.0 / 360.0)
    hp = P - ndimage.gaussian_filter1d(P, _sig, axis=0, mode="wrap")
    amp = float(np.std(hp))
    if not np.isfinite(amp) or amp <= 0:
        return 0.0, 0.0
    a, b = hp[:, :-1].ravel(), hp[:, 1:].ravel()
    if a.size < 100 or np.std(a) <= 0 or np.std(b) <= 0:
        return amp, 0.0
    coh = float(np.corrcoef(a, b)[0, 1])
    return amp, (coh if np.isfinite(coh) else 0.0)


def _pick_weight_alpha(stacks_half, sat_half, secs, cal, abs_shift, track,
                       progress, alphas=(1.0, 0.85, 0.7, 0.55)):
    """Choose the exposure exponent in the merge weight, by measurement.

    THE PROBLEM. The merge weights each tier by `s * saturation_rolloff`.
    Weighting by exposure time is statistically optimal for PHOTON NOISE --
    a tier that collected twice the photons deserves twice the say. It is not
    optimal for RESOLUTION, because a long exposure beside a blindingly bright
    edge is degraded by glare and bloom well before it hard-clips, and the
    rolloff does not start until 0.87 of saturation. So at the limb the
    picture goes to the longest not-quite-saturated tier, which is also the
    most smeared one: on the reference bracket that is 67x more weight than
    the sharpest tier carries.

    Measured there, the four shortest tiers hold 1.9x more radially-coherent
    fine structure at 1.01-1.10 R than the merge does -- all the way round the
    disc, 1.38x in the quietest sector and 3.22x in the busiest. That is
    detail the data contains and the picture does not.

    THE KNOB. `w = s**alpha * rolloff`. alpha = 1 is what the merge has always
    done; below 1 it tilts toward the shorter, sharper tiers. Dimensionless,
    so a half-res trial transfers to the full-res merge exactly -- which a
    knee expressed in raw ADU would not.

    WHY IT IS GATED, NOT SET. The same measurement says a blanket tilt would
    be a bad trade: beyond 1.3 R the short tiers are noise (coherence 0.19,
    then 0.013), and giving them weight out there buys limb detail with
    streamers. So the outer shells are guards -- an alpha that improves the
    limb but costs more than 2% of the coherent structure at 1.3-1.8 R or
    1.8-2.6 R does not win. alpha = 1.0 is the default and stays unless
    something measurably beats it.

    This is the harness 0.16.0 should have had: the merge is not changed on
    an argument, it is changed on a number from the data in front of it.
    """
    # `track` IS THE ONLY THING THIS NEEDS, and taking anything else was a bug
    # (0.22.87). The parameter used to include `Rmoon`, which the caller clears
    # whenever the moon MASK is rejected -- a routine outcome, and an unrelated
    # decision. Rejecting the mask therefore switched off this measurement
    # silently, and the run reported "no per-tier lunar track" while holding a
    # perfectly good one.
    #
    # Found on the 600 mm bracket: deleting one already-discarded frame shifted
    # its tier's MEAN TIMESTAMP, which moved a point on the lunar track fit
    # (0.00 -> 0.13 px/s -- the second figure being the physically sensible one
    # for 62 s at this plate scale), which flipped the mask trial from applied
    # to rejected, which cleared Rmoon, which left the exponent at 1.0 instead
    # of 0.55. The picture paid for it: merged limb 20-80% 7.0 -> 13.0 px, limb
    # fit rms 0.40 -> 5.44 px, disc mask margin 8.1 -> 15.2 px, corona 3.9 ->
    # 3.8 R. None of those steps is wrong on its own; the coupling is.
    #
    # Everything below reads `track[s]` for the per-tier disc and takes R0 from
    # its median, and the R0 check two blocks down already rejects a radius
    # that is not usable -- so Rmoon was carrying no information this function
    # did not already have.
    if not track:
        return 1.0, {"verdict": "no per-tier lunar track; left at 1.0"}
    H, W = stacks_half[secs[0]].shape
    yv = np.arange(H, dtype=np.float32)
    xv = np.arange(W, dtype=np.float32)

    def build(alpha):
        acc = np.zeros((H, W), np.float32)
        wsum = np.zeros((H, W), np.float32)
        for s in secs:
            a = stacks_half[s].astype(np.float32) / np.float32(s * cal[s])
            ady, adx = abs_shift[s]
            g = (~sat_half[s]).astype(np.float32)
            a = ndimage.shift(a, (ady, adx), order=1, mode="nearest")
            g = ndimage.shift(g, (ady, adx), order=1, mode="nearest", cval=0)
            w = np.float32(s ** alpha) * ndimage.gaussian_filter(g, 10)
            if s in track:
                cyt, cxt, Rt = (v / 2.0 for v in track[s])
                dt = np.hypot(yv[:, None] - cyt, xv[None, :] - cxt)
                w = w * np.clip((dt - (Rt - 0.5)) / 2.0, 0.0, 1.0)
            acc += w * a
            wsum += w
        return acc / np.maximum(wsum, 1e-9)

    cy0 = float(np.median([track[s][0] for s in track])) / 2.0
    cx0 = float(np.median([track[s][1] for s in track])) / 2.0
    R0 = float(np.median([track[s][2] for s in track])) / 2.0
    if not np.isfinite(R0) or R0 <= 4:
        return 1.0, {"verdict": "lunar radius not usable; left at 1.0"}
    SHELLS = ((1.02, 1.12, "limb"), (1.30, 1.80, "mid"), (1.80, 2.60, "outer"))
    rows, base = {}, None
    for al in alphas:
        try:
            m = np.clip(build(al), 0, None)
            sc = {}
            for r0, r1, nm in SHELLS:
                if r1 * R0 > 0.48 * min(H, W):     # shell off the frame
                    continue
                sc[nm] = _fine_structure(m, cy0, cx0, R0, r0, r1)
            del m
        except Exception as e:
            return 1.0, {"verdict": f"weight trial failed ({e}); left at 1.0"}
        rows[al] = sc
        if base is None:
            base = sc
    if "limb" not in base:
        return 1.0, {"verdict": "limb shell not measurable; left at 1.0"}
    # NO GUARD SHELL, NO CHANGE. On a tight crop -- a long lens filling the
    # frame with the disc -- the mid and outer shells fall off the edge and are
    # skipped, which would leave the limb score to win unopposed. Tilting the
    # merge with nothing watching the outer field is exactly the trade this
    # trial exists to refuse. Caught by its own test case: on a 120 px frame
    # with R=100 it picked 0.55 on pure noise.
    if not any(k in base for k in ("mid", "outer")):
        return 1.0, {"verdict": "no outer shell in frame to guard with; "
                                "left at 1.0"}
    # THE GUARD IS COHERENCE, NOT A DROP IN THE SCORE.
    #
    # Giving weight to a noise-dominated tier RAISES amp in every shell, so a
    # test that only asks "did any score fall" passes it happily -- the trial's
    # own noise case scored +41% at the limb and would have been accepted.
    # What noise cannot fake is radial continuity. Measured on that case, as
    # alpha drops from 1.0 to 0.55 the coherence falls 0.933 -> 0.438 at the
    # limb and 0.398 -> 0.225 in the mid field, monotonically, while in the
    # genuine glare case it barely moves. So an alpha that costs more than 0.05
    # of coherence anywhere is buying grain, not detail, and is refused.
    _score = lambda sc, k: sc[k][0] * max(sc[k][1], 0.0)
    best, best_gain = 1.0, 0.0
    for al, sc in rows.items():
        if al == 1.0 or "limb" not in sc:
            continue
        gain = _score(sc, "limb") / max(_score(base, "limb"), 1e-12) - 1.0
        dcoh = max(base[k][1] - sc[k][1] for k in sc if k in base)
        if gain > max(best_gain, 0.05) and dcoh <= 0.05:
            best, best_gain = al, gain
    info = {"alpha": best, "gain_limb": best_gain,
            "scores": {str(a): {k: [round(v[0], 5), round(v[1], 3)]
                                for k, v in sc.items()}
                       for a, sc in rows.items()}}
    _fmt = "  ".join(
        f"a={a:.2f}:" + "/".join(f"{sc[k][0] * max(sc[k][1], 0):.4f}"
                                 for k in ("limb", "mid", "outer") if k in sc)
        for a, sc in rows.items())
    progress.log(f"merge weight trial (coherent detail limb/mid/outer): {_fmt} "
                 f"-> exposure exponent {best:.2f}"
                 + (f", {100 * best_gain:+.0f}% at the limb" if best != 1.0
                    else " (unchanged — nothing beat it without costing the outer field)"),
                 None)
    return best, info


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


# Low end of the merge weight, as a fraction of sensor saturation. 0 restores
# every build before 0.22.5 exactly. 0.005 is deliberately half the best value
# measured (see the note where it is used): the reconstruction that measured it
# could not separate faint real structure from noise, so this takes the part of
# the gain that is not in doubt.
_MERGE_FLOOR = 0.005

#: Display names for the two feather forms, in ONE place. The log printed the
#: internal value ('plain') while the report printed a UI label ('Detail'), so
#: the same setting appeared under two names in the same run and neither said
#: what it does. These name the MECHANISM: one blurs the weight across each
#: tier's clipping edge, which joins the tiers seamlessly and lets weight leak
#: from clipped pixels into unclipped ones; the other keeps the weight exact
#: at that edge, which is radiometrically correct and leaves the join visible.
FEATHER_NAMES = {"plain": "Blended edge", "taper": "Exact edge",
                 "masked": "Masked edge"}

# WHICH FEATHER SHIPS BY DEFAULT, AND WHY IT IS THE ONE WITH THE KNOWN ERROR.
#
# 0.22.16 replaced the plain blur with a leak-free weight and fixed the radial
# profile exactly. It also made a ring artifact visible that the plain blur had
# been hiding, on every dataset, in MGN, FNRGF and NAFE alike. Fifteen weight
# forms were then rebuilt from the reference set's aligned tiers and scored (tools/):
#
#   variant                              1.02R    ring power
#   none (no feather at all)             1.000       0.94
#   the shipped taper                    1.000       1.00
#   plain blur (<= 0.22.15)              0.747       0.66
#   hard collar exclusion, 2-16 px       0.998    1.01-1.02
#   inverse-variance weighting           0.992       1.15
#   C1 roll-off to zero below the cut    1.002       0.99
#   weight from a smoothed intensity     1.006    0.96-0.97
#   two-pass fill + plain blur           0.993       1.00
#   LDIC azimuthal affine (thesis 4.15)  0.809       0.87
#
# The unfeathered merge -- no blur, no mask, no contour term anywhere -- carries
# the rings at full strength, so no property of the weight creates them. And the
# two-pass fill is decisive about what the plain blur actually does: give the
# clipped collar a correct value instead of the clipped tier's under-report and
# the photometry comes back to 0.993 while the rings return to 1.00. The plain
# blur does not remove the artifact. It covers it with a compensating error of
# the same shape -- both are bounded by each tier's saturation contour.
#
# So the plain feather is a deliberate, documented trade of brightness near the
# limb for pictures that do not show rings -- and the price is NOT a constant.
# 25% of the true level at 1.02 R on the 600 mm reference set; a factor of EIGHT on
# the test set's 360 mm, where it prints as a pink rim. He reported exactly that on
# 0.22.26, which is what killed this as a global default.
#
# Since 0.22.28 the feather is chosen per dataset by measuring that ratio on
# the half-resolution stacks before the merge -- see _pick_feather. This
# constant is only the fallback for when that measurement cannot run.
#
# The artifact itself is NOT understood. What is known: it follows each tier's
# saturation contour; it survives per-tier correction by scale, by colour, by
# radius and by signal level; adjacent tiers disagree by up to 3x within 0-8 px
# outside the longer one's saturated region and by ~1% beyond 8-16 px. See
# TODO 1e and 2.
_FEATHER_DEFAULT = "plain"

# Largest pedestal this will believe, as a fraction of sensor saturation.
# 0.002 is 33 ADU of a 14-bit raw; the three real datasets that show one land at
# 2.9, 5.5 and -2.0 ADU. Anything past this is not a black-level residual and
# the fit has locked onto something else.
_PEDESTAL_MAX = 0.002
# The same ceiling for a file that reported NO black level, so the whole of it
# is still in the pixels rather than a residue of it. 0.02 of saturation is
# 1310 ADU on a 16-bit FITS and 328 on a 14-bit raw: wide enough for a real
# black level -- a third tester's measures 142 -- and still far below anything that
# could be corona. See TODO 22.
_PEDESTAL_UNKNOWN_MAX = 0.02
# How large a per-link offset the linear-fit chain will accept, as a fraction of
# the median signal the link was fitted on. This term models a black-level
# residual, which is a few ADU; anything approaching the signal itself is the
# fit trading slope against offset on data that cannot constrain both. On
# the 560 mm 2024 test set every well-constrained link needed 0.5-3% and the
# one that ran away needed 60%; on the 600 mm reference set the largest is 4%. 0.20
# sits an order of magnitude above the honest ones and well below the failure.
_OFFSET_MAX_FRAC = 0.20


def _pedestal_score(prof, scale, P):
    """Tier-to-tier disagreement in the outer field for a candidate pedestal.

    `prof` is (tiers x shells) of ring medians in RAW units, NaN where a ring is
    clipped or missing; `scale` is s*cal[s] per tier. Subtract P, divide each
    tier by its own scale, and every tier should then report the SAME scene. The
    score is how much they disagree, so the right P minimises it.
    """
    with np.errstate(all="ignore"):
        A = (prof - P) / scale[:, None]
        ref = np.nanmedian(A, axis=0)
        ref = np.where(np.abs(ref) < 1e-12, np.nan, ref)
        sc = np.nanstd(A / ref, axis=0)
    sc = sc[np.isfinite(sc)]
    return float(np.median(sc)) if sc.size else np.inf


def _photometry_requested(folder, requested=None):
    """Is the per-channel linear-fit photometry SELECTED for this folder?

    One predicate, used by opts.json and by server.py's cache key, because they
    used to disagree: the pipeline recorded the mode that was APPLIED and the
    server compared it against the mode that was ASKED FOR, so any fallback made
    the two permanently unequal and forced a full re-stack on every Start.

    ON BY DEFAULT SINCE 0.22.86. It was opt-in through a marker file from
    0.22.68, which is how it should have started -- but the file is invisible
    once dropped, and that cost real work: five of the test set's test runs were
    reported and discussed as linear-fit results when the folder had no marker
    in it, and the reference 600 mm folder turned out to have had one all along
    without anyone remembering. A setting nobody can see is a setting nobody can
    check.

    `requested` is the toolbar's choice, "linfit" or "scalar"; the environment
    variable still overrides both, in either direction, for a bisect. The legacy
    marker files are honoured but no longer needed, and can be deleted.
    """
    _env = os.environ.get("ECLIPSEFORGE_PHOTOMETRY")
    if _env in ("linfit", "hill"):
        return True
    if _env in ("scalar", "off", "none"):
        return False
    if any(os.path.exists(os.path.join(folder, _f)) for _f in
           ("linear_fit_photometry.txt", "hill_photometry.txt")):
        return True
    if requested is not None:
        return str(requested) != "scalar"
    return True


def _ring_sample(stack, sat, yr, xr, shift):
    """Sample one tier on the COMMON frame's rings, in the tier's own frame.

    THE BUG THIS EXISTS TO CLOSE (0.22.76). Three estimators -- the shared
    pedestal, the per-tier radial check and the LDIC azimuthal fit -- sample
    `stacks_half[s]` on rings centred at the mid tier's centre. But
    `stacks_half` is intra-tier aligned only; `abs_shift[s]` is applied in the
    merge loop, after all three have run. A tier offset by d px therefore had
    its radial profile read at r + d*cos(phi - phi0).

    What that costs each of them is different, and worth knowing:

      - LDIC fits k(phi) per azimuth, and a cos(phi) modulation IS a dipole in
        k. It went straight into the fit. tools/ldic_shift_bench.py measures
        a spurious 65% k spread at 5 px of drift where the truth is 0%.
      - the pedestal and the radial check take a MEDIAN over 720 azimuths,
        which cancels the dipole -- but not the second-order term, because the
        profile is convex. Measured on the same fixture, the error in the ring
        median at 11 px of drift runs -5.6% at 1.5 R to -1.4% at 3.5 R. A
        per-tier error that varies with radius is exactly what the pedestal fit
        is built to read as an offset, so it does not cancel there either.

    ndimage.shift(img, +s) gives out[u] = in[u - s] (see the sign note above
    _moon_track), so the point that lands at (y, x) after registration sits at
    (y - ady, x - adx) in the tier's own frame.

    Bilinear for the signal, because the shifts are sub-pixel and an integer
    cast puts a fraction of the same dipole back (4-6% against 1.7% on the
    bench). NEAREST for the saturation mask, because validity is not signal and
    interpolating it invents partial validity along every clipped edge -- the
    same reason the merge resamples its masks at order=1.
    """
    ady, adx = shift
    yy = np.clip(yr - ady, 0, stack.shape[0] - 1)
    xx = np.clip(xr - adx, 0, stack.shape[1] - 1)
    v = ndimage.map_coordinates(stack.astype(np.float64), [yy, xx],
                                order=1, mode="nearest")
    c = ndimage.map_coordinates(sat.astype(np.uint8), [yy, xx],
                                order=0, mode="nearest").astype(bool)
    return v, c


def _fit_pedestal(prof, scale, pmax):
    """One additive offset, in raw units, common to every tier.

    WHY THIS EXISTS. A tier's measured value is `S(r)*s*cal + P`: the scene times
    the exposure, plus whatever constant the black subtraction left behind. Our
    merge divides by `s*cal` and averages, so that leftover P arrives divided by
    the exposure time -- negligible on a long tier, enormous on a short one. It
    is invisible on the bright inner corona and it dominates the outer field.

    MEASURED on the exported tiers of three real datasets, as the scatter
    between tiers over 1.5-3.5 R once each is put back on the same scene scale:

        set                tiers   no offset   one global offset   fitted P
        360 mm test set        12      39.20%          2.87%        +2.87 ADU14
        560 mm 2024 test set   14      36.43%          5.17%        +5.49 ADU14
        250 mm test set         9       2.42%          1.72%        -1.97 ADU14

    Thirteen-fold and sevenfold reductions from ONE number. On the 360 mm set
    the 1/1000 s tier reads 2.76x the scene at 2.9 R before the correction and
    1.18x after. The 250 mm set (a different body) barely has one, which is the
    control: this does not invent an offset where there is none.

    ONE NUMBER, NOT ONE PER TIER. Letting every tier have its own offset only
    moved 2.87% to 2.60% (360 mm) and 5.17% to 4.67% (560 mm), and the long
    tiers' values wandered to -7 ADU because a tier with plenty of signal does
    not constrain its own pedestal at all. A black-level residual is a property
    of the sensor and the session, not of the shutter speed, so the shared
    number is both the better-determined and the more honest model.

    THE PRIOR IS ZERO. P is shrunk by its own significance -- P*P^2/(P^2+sig^2)
    with sig from a leave-one-tier-out jackknife -- so a dataset whose tiers
    disagree for some other reason gets its correction pulled toward nothing
    rather than having one invented for it.

    `prof`, `scale` and `pmax` are all in RAW ADU, which is what stacks_half
    holds; the returned pedestal is in the same units and is subtracted from the
    demosaiced tier before white balance, where it is still one number per
    photosite rather than three per pixel.

    Returns (P_applied, P_raw, sigma, score_before, score_after).
    """
    ok = np.isfinite(prof).sum(axis=1) >= 8
    if ok.sum() < 4:
        return 0.0, 0.0, 0.0, np.nan, np.nan
    prof, scale = prof[ok], scale[ok]
    # 401 points over +-33 ADU is a 0.16 ADU step; over +-1310 it would be
    # 6.6, which is coarse against a black level we are trying to measure to
    # the ADU. Keep the step at or below 1 ADU, capped so the jackknife below
    # (one full re-fit per tier) stays cheap.
    _n = int(min(max(401, 2 * pmax + 1), 4001))
    grid = np.linspace(-pmax, pmax, _n)

    def _best(pr, sc):
        v = [_pedestal_score(pr, sc, p) for p in grid]
        return float(grid[int(np.argmin(v))]), float(np.min(v))

    P, after = _best(prof, scale)
    before = _pedestal_score(prof, scale, 0.0)
    if not np.isfinite(before) or not np.isfinite(after):
        return 0.0, 0.0, 0.0, before, after
    # leave-one-tier-out: a real pedestal does not depend on which tier is in
    n = prof.shape[0]
    jk = []
    for i in range(n):
        k = [j for j in range(n) if j != i]
        jk.append(_best(prof[k], scale[k])[0])
    jk = np.asarray(jk, float)
    sig = float(np.std(jk) * np.sqrt(max(n - 1, 1)))
    shrink = P * P / (P * P + sig * sig) if (P or sig) else 0.0
    return float(P * shrink), float(P), sig, before, after


# Guards for _fit_channel_floors. Both were arrived at by being wrong first;
# tools/pedestal_probe.py is the bench that caught each one.
_CH_MAX_RESID_FRAC = 0.03    # of the largest ring median
_CH_MIN_TIERS = 6            # a two-parameter fit per channel needs room
_CH_MIN_SPAN = 8.0           # shortest-to-longest exposure ratio
_CH_MAX_LADDER_DEV = 0.50    # fitted exposure vs the stated one


def _solve_photo_network(n, links, ref):
    """Least-squares log-factors from a set of measured links.

    `links` is [(i, j, log_ratio, n_pixels)]. Minimises
    sum_k w_k * (logf[j] - logf[i] - lr_k)^2 with logf[ref] pinned, the weight
    being sqrt(n_pixels) -- the standard error of a median falls as 1/sqrt(N).

    Returns (logf, residuals) or (None, None) if the system is not worth
    solving. It is not worth solving when some tier no link reaches: least
    squares would hand back a minimum-norm number for that tier, which looks
    like an answer and is not. Same trap as the alignment network.

    WHEN THIS DIFFERS FROM A CHAIN, which is narrower than it first looks. If
    the links are mutually CONSISTENT the network reproduces the chain exactly
    -- a tier that is genuinely 18% bright makes its two adjacent links read
    1.18 and 0.85 and the two-step link across it read 1.00, and those three
    agree, so there is nothing to arbitrate. What the network buys is the case
    where links CONTRADICT each other: independent measurement error on each,
    or one link measured over a region that misleads it. Then the chain
    believes whichever links lie on its path and the network averages all of
    them.
    """
    if n < 3 or len(links) < n:
        return None, None
    seen = set()
    for i, j, _, _ in links:
        seen.add(int(i)); seen.add(int(j))
    if len(seen) != n:
        return None, None
    wmax = max(np.sqrt(float(w)) for _, _, _, w in links)
    if not np.isfinite(wmax) or wmax <= 0:
        return None, None
    A, b, w = [], [], []
    for i, j, lr, npx in links:
        row = np.zeros(n); row[int(j)] = 1.0; row[int(i)] = -1.0
        A.append(row); b.append(float(lr))
        w.append(float(np.sqrt(float(npx)) / wmax))
    row = np.zeros(n); row[int(ref)] = 1.0
    A.append(row); b.append(0.0); w.append(100.0)
    A = np.asarray(A); w = np.asarray(w); b = np.asarray(b)
    try:
        sol = np.linalg.lstsq(A * w[:, None], b * w, rcond=None)[0]
    except np.linalg.LinAlgError:
        return None, None
    if not np.all(np.isfinite(sol)):
        return None, None
    sol = sol - sol[int(ref)]
    res = [float(sol[int(j)] - sol[int(i)] - float(lr)) for i, j, lr, _ in links]
    return sol, res


def _fit_channel_floors(pmed, secs, iters=200, exposures=None):
    """Per-channel black level, as DIFFERENCES about their own mean, in ADU.

    `pmed` is {sec: [R, G, B]} of outer-field ring medians in raw ADU.
    Returns (deltas, info). `deltas` sums to zero and is subtracted alongside
    the shared scalar pedestal, which keeps the common part.

    WHY NOT THE ESTIMATOR IN `_tier_colour_check`. That one searches a single
    additive offset per channel that minimises tier-to-tier disagreement, and
    its own comment says why it only reports: given a synthetic bracket with a
    pure 3%-per-step colour gain and no offset at all it returns -17 and +17
    ADU for R and B. An offset and a gain error are degenerate when only one of
    them is free.

    So both are free here. Over the tiers, one channel's ring median obeys

        m_c(i) = a_c * s_i + b_c

    with `s_i` ONE effective exposure per tier, SHARED by the three channels.
    A shutter that did not honour the requested time, a wrong EXPTIME, a gain
    compounding along the ladder -- all of these are achromatic, so they land
    in `s_i` and cannot reach `b_c`. On Val's set the fitted ladder comes out
    -15.5% at 1/20 s and -10.4% at 1/13 s against the stated EXPTIME, which is
    the same shape found independently from corner-patch rates; the floors come
    out regardless.

    TWO THINGS THIS DELIBERATELY GIVES UP.

    First, the absolute level. Only the DIFFERENCES between channels are
    returned. The common part is already fitted, jackknifed and shrunk toward
    zero by `_fit_pedestal`, and this does not touch it.

    Second, the part of those differences that a ladder error could have
    produced. `(s, b)` and `(s + d, b - a*d)` fit identically for any constant
    `d`, so floors PROPORTIONAL TO EACH CHANNEL'S RATE cannot be told from a
    shutter-lag shift. The bench shows this directly: given a pure achromatic
    gain and no floors at all, the raw deltas come out (+14.2, -1.8, -12.4),
    and dividing those by the mean-removed rates gives one number. That
    component is projected out. The cost is under-correcting when the true
    floors happen to lie along `a`; the gain is that no ladder error can ever
    be turned into a colour cast. That is the right way round.
    """
    secs = [s for s in secs if s in pmed]
    t = np.asarray(secs, float)
    info = {"n": int(t.size),
            "span": float(t.max() / t.min()) if t.size else 0.0}
    if t.size < _CH_MIN_TIERS or info["span"] < _CH_MIN_SPAN:
        info["usable"] = False
        info["why"] = ("%d tier(s) over %.1fx exposure; needs %d over %.0fx"
                       % (t.size, info["span"], _CH_MIN_TIERS, _CH_MIN_SPAN))
        return np.zeros(3, np.float64), info
    M = np.asarray([[float(pmed[s][c]) for s in secs] for c in range(3)], float)
    if not np.isfinite(M).all():
        info["usable"] = False
        info["why"] = "a ring median is not finite"
        return np.zeros(3, np.float64), info

    # THE LADDER, PINNED RATHER THAN FREE, when a measured one is offered.
    #
    # `sv` free per tier is what the docstring above argues for, and on data
    # whose stated times are roughly right it is fine. On a third tester's FITS it is
    # not: his headers carry a synthetic 1.5^n ladder ("Manual group
    # assignment" in Siril's own comment on the card) where the corner rates
    # measure ~1.25x per step, and the SHORT tiers have almost no sky in the
    # corner to constrain their own exposure. The fit then escapes into the
    # ladder -- +6001% on the shortest tier -- the guard below correctly refuses
    # it, and the floors are not applied on the one dataset that most needs
    # them.
    #
    # `cal` already measures the effective exposure of every tier, on the BRIGHT
    # corona where the signal-to-noise is high, and the whole pipeline already
    # trusts it: `s * cal[s]` is the quantity smoke_blind_exposure proves the
    # merge runs on. Handing that in as the ladder removes the runaway parameter
    # without weakening the protection this function exists to provide -- an
    # achromatic error still lands in the COMMON part of b, which `b - b.mean()`
    # discards, and whatever survives along the rates is still projected out
    # below. So the defence is unchanged; only the unconstrained direction is
    # gone.
    if exposures:
        sv = np.asarray([float(exposures.get(x, x)) for x in secs], float)
        info["ladder"] = "pinned to the measured effective exposure (s * cal)"
    else:
        sv = t.copy()
        info["ladder"] = "free per tier"
    a = np.zeros(3); b = np.zeros(3)
    for _ in range(1 if exposures else iters):
        A = np.vstack([sv, np.ones_like(sv)]).T
        for c in range(3):
            sol, *_ = np.linalg.lstsq(A, M[c], rcond=None)
            a[c], b[c] = float(sol[0]), float(sol[1])
        den = float((a * a).sum())
        if den <= 0:
            break
        sn = ((a[:, None] * (M - b[:, None])).sum(axis=0)) / den
        # s and a trade off; pin the SCALE against the stated ladder and leave
        # the shape free, which is the whole point of the parameter
        sn = sn * float((sn * t).sum() / max((sn * sn).sum(), 1e-30))
        if np.max(np.abs(sn - sv)) < 1e-12 * max(1.0, float(t.max())):
            sv = sn
            break
        sv = sn

    R = M - (a[:, None] * sv + b[:, None])
    info["max_resid"] = float(np.abs(R).max())
    info["resid_frac"] = float(np.abs(R).max() / max(np.abs(M).max(), 1e-9))
    info["rates"] = [round(float(x), 1) for x in a]
    info["ladder_dev_pct"] = [round(100.0 * (float(x) / float(y) - 1.0), 1)
                              for x, y in zip(sv, t)]
    # THE RESIDUAL IS NOT A GUARD ON ITS OWN. `sv` carries one free parameter
    # per tier, so the model can fit almost any set of ring medians and still
    # report a small residual -- the first real run did exactly that, returning
    # a NEGATIVE effective exposure for the shortest tier and a residual that
    # looked fine. An exposure cannot be negative, and a shutter that missed
    # its mark by more than half is not the small achromatic correction this
    # parameter exists to absorb; either way the fit has stopped describing the
    # thing it was built to describe, and the floors that come with it are not
    # black levels.
    _dev = np.abs(sv / np.maximum(t, 1e-12) - 1.0)
    info["ladder_max_dev"] = float(_dev.max())
    # The ladder guard exists because a FREE exposure per tier can hide a bad
    # fit. A pinned ladder has no such freedom -- its deviation from the stated
    # times is a measurement, not a fitted parameter, and on data with wrong
    # headers it is SUPPOSED to be large. Only the residual guard applies then.
    if (sv <= 0).any() or (not exposures
                           and float(_dev.max()) > _CH_MAX_LADDER_DEV):
        info["usable"] = False
        info["why"] = ("the fitted exposure ladder ran to %+.0f%% of the stated "
                       "times (limit %+.0f%%%s); with one free exposure per tier "
                       "a bad fit hides in the ladder rather than in the "
                       "residual, so the floors are not trusted"
                       % (100 * float(_dev.max()), 100 * _CH_MAX_LADDER_DEV,
                          ", and one came out negative" if (sv <= 0).any() else ""))
        return np.zeros(3, np.float64), info
    info["usable"] = info["resid_frac"] <= _CH_MAX_RESID_FRAC
    if not info["usable"]:
        info["why"] = ("even with the ladder free the channels miss by %.1f ADU "
                       "(%.1f%%); something beyond a floor and an achromatic "
                       "exposure error is present"
                       % (info["max_resid"], 100 * info["resid_frac"]))
        return np.zeros(3, np.float64), info

    braw = b - b.mean()
    arel = a - a.mean()
    den = float((arel * arel).sum())
    proj = (float((braw * arel).sum()) / den) * arel if den > 1e-12 \
        else np.zeros(3)
    # THE PROJECTION DEFENDS AGAINST A FREEDOM THAT A PINNED LADDER DOES NOT
    # HAVE. `(s, b)` and `(s + d, b - a*d)` fit identically for any constant d
    # -- that is why the component of the floors parallel to the rates is
    # discarded. But `d` only exists when `s` is a free parameter. Pinned to
    # `s * cal`, measured independently on the bright corona, there is no d to
    # absorb and the component along `a` is as real as any other.
    #
    # MEASURED on a third tester's set, where the true floors happen to lie largely
    # along the rate direction and the projection therefore throws most of the
    # correction away. Floors extrapolated from his two shortest frames, as
    # deltas about the mean, against what each mode returns:
    #
    #     truth        R -46.1  G -29.9  B +76.1
    #     projected    R +15.2  G -40.6  B +25.4   (31% of the blue excess)
    #     full         see tools/smoke_channel_floors.py
    #
    # The projected answer is not noise -- it is the truth minus its own
    # component along `a`, which reproduces to a decimal by hand. It is the
    # method working as designed on data the design did not anticipate.
    if exposures and not _CH_PROJECT_PINNED:
        out = braw
        info["projection"] = "skipped (the ladder is pinned, so nothing is degenerate with it)"
    else:
        out = braw - proj
        info["projection"] = "applied (the ladder is free and degenerate with the floors)"
    info["deltas_raw"] = [round(float(x), 2) for x in braw]
    info["deltas_ladder_part"] = [round(float(x), 2) for x in np.asarray(proj)]
    info["deltas"] = [round(float(x), 2) for x in out]
    return out, info


#: Keep projecting the ladder-parallel component out even when the ladder is
#: pinned to a measured one. The projection exists for a degeneracy that only a
#: FREE ladder has; with `s * cal` handed in there is nothing to be degenerate
#: with, and on data whose floors lie along the rates it discards most of a real
#: correction. True by default would be the conservative choice and is measured
#: against False in tools/smoke_channel_floors.py.
#: Drop a tier from the floor fit once its outer-field level exceeds this
#: multiple of the faintest tier's. Past that the floor is a few percent of what
#: the tier reads and contributes nothing, while any departure from a constant
#: sky rate contributes everything.
_CH_MAX_LEVEL_RATIO = 10.0

_CH_PROJECT_PINNED = False

_LDIC_SEGMENTS = 60
_LDIC_ORDER = 4
_LDIC_TABLE = 16384          # azimuth bins the merge evaluates k and q on


def _ldic_tables(k, q, ck, cq):
    """k(phi) and q(phi) sampled on `_LDIC_TABLE` bins, for the merge to index.

    Both are order-4 trigonometric polynomials, so this is the exact function,
    not an interpolation of the 60 fitted samples. `ck`/`cq` are None only on
    the degenerate path in `smooth()`, where the fit fell back to a constant.
    """
    if ck is None or cq is None:
        return (np.full(_LDIC_TABLE, np.float32(k[0])),
                np.full(_LDIC_TABLE, np.float32(q[0])))
    ang = np.arange(_LDIC_TABLE, dtype=np.float64) / _LDIC_TABLE * 2 * np.pi
    kt = np.full(_LDIC_TABLE, float(ck[0]))
    qt = np.full(_LDIC_TABLE, float(cq[0]))
    for o in range(1, _LDIC_ORDER + 1):
        c, s = np.cos(o * ang), np.sin(o * ang)
        kt += ck[2 * o - 1] * c + ck[2 * o] * s
        qt += cq[2 * o - 1] * c + cq[2 * o] * s
    return kt.astype(np.float32), qt.astype(np.float32)


def _fit_azimuthal_affine(vals, wts, secs_order):
    """Per-tier affine transform that varies with AZIMUTH -- Druckmullerova,
    doctoral thesis eq. 4.15, the LDIC composition.

        g(r,phi) = SUM_i  w(f_i) * ( k_i(phi) * f_i(r,phi) + q_i(phi) )

    k and q are fitted by linear regression in 60 angular segments against the
    composite accumulated so far, starting from the longest exposure, then
    smoothed with a trigonometric polynomial of low order. The thesis says what
    they are for: "to compose images with different distribution of diffuse
    light in the optical system ... or even images that were taken through thin
    clouds." Diffuse light off a large saturated area differs from tier to tier
    and is not axisymmetric, so one scalar per tier cannot express it.

    WHAT WE HAD: one scalar `cal[s]` per tier plus one shared additive
    pedestal. Fitted on the 600 mm reference set, k varies 3-27% around the limb
    depending on the tier and q runs to 9% of the local signal. That is the
    part a scalar cannot reach, and it is azimuthal -- so it changes wherever
    the mix of tiers changes, which is along each tier's saturation contour.

    NORMALISED, DELIBERATELY. The raw transform re-references every tier onto
    the longest exposure and moves the whole radial profile by ~20%, discarding
    the photometric chain and the shared pedestal that are already fitted and
    already measured. Taking the azimuthal MEAN out of k and q leaves exactly
    the part those two cannot express and changes nothing else: measured on the
    bench, ring power 0.87x with the radial profile held inside 4% (the raw
    form scores the same 0.87x and costs 20%).

    Honest about size: 0.87x is a real improvement and not a fix. The ring
    artifact survives it, as it survives every per-tier correction tried so far
    -- by scale, by colour, by radius, by signal level. See TODO 1e.

    `vals[i]` and `wts[i]` are one tier's sampled radiance and weight on a
    (radius x angle) ring grid, angles evenly spaced from 0. Returns
    {sec: (k[NS], q[NS])}.
    """
    NS, NA = _LDIC_SEGMENTS, vals[0].shape[1]
    per = max(NA // NS, 1)
    cen = (np.arange(NS) + 0.5) / NS * 2 * np.pi
    B = [np.ones(NS)]
    for o in range(1, _LDIC_ORDER + 1):
        B += [np.cos(o * cen), np.sin(o * cen)]
    B = np.vstack(B).T

    def smooth(v, fallback):
        """Returns (values at the 60 segment centres, trig coefficients or None).

        The COEFFICIENTS are what the caller wants. k and q are band-limited by
        construction -- an order-4 trigonometric polynomial -- so they can be
        evaluated at any azimuth exactly. Returning only the 60 samples was the
        0.22.26 mistake: the merge then looked the value up by segment index,
        which turns a smooth function into a 60-step staircase and prints a
        radial edge every 6 degrees across the whole frame.
        """
        m = np.isfinite(v)
        if m.sum() < 2 * _LDIC_ORDER + 3:
            return np.full(NS, fallback if not m.any() else np.nanmedian(v)), None
        sol, *_ = np.linalg.lstsq(B[m], v[m], rcond=None)
        return B @ sol, sol

    out = {}
    C = np.zeros_like(vals[0]); Wt = np.zeros_like(vals[0])
    for i, sec in enumerate(secs_order):
        f, w = vals[i], wts[i]
        if i and Wt.max() > 0:
            comp = np.where(Wt > 1e-9, C / np.maximum(Wt, 1e-9), np.nan)
            ks = np.full(NS, np.nan); qs = np.full(NS, np.nan)
            for j in range(NS):
                sl = slice(j * per, (j + 1) * per)
                x = f[:, sl].ravel(); y = comp[:, sl].ravel()
                ok = (np.isfinite(x) & np.isfinite(y) & (w[:, sl].ravel() > 0.5)
                      & (Wt[:, sl].ravel() > 0.05) & (x > 0))
                if ok.sum() < 50:
                    continue
                A = np.vstack([x[ok], np.ones(int(ok.sum()))]).T
                sol, *_ = np.linalg.lstsq(A, y[ok], rcond=None)
                ks[j], qs[j] = sol
            if np.isfinite(ks).sum() >= 20:
                k, ck = smooth(ks, 1.0); q, cq = smooth(qs, 0.0)
                km = float(np.mean(k))
                if np.isfinite(km) and km > 1e-6:
                    k = k / km                       # azimuthal mean 1
                    q = q - float(np.mean(q))        # azimuthal mean 0
                    # the same normalisation on the coefficients, so the
                    # continuous form and the sampled form are the same
                    # function: scale everything by km, then drop the constant
                    # term of q (B's first column is the constant 1).
                    if ck is not None:
                        ck = np.asarray(ck, np.float64) / km
                    if cq is not None:
                        cq = np.asarray(cq, np.float64).copy()
                        cq[0] = 0.0
                    if np.isfinite(k).all() and np.isfinite(q).all():
                        out[sec] = (k.astype(np.float32), q.astype(np.float32),
                                    ck, cq)
        adj = f
        if sec in out:
            kk, qq, ck, cq = out[sec]
            if ck is not None and cq is not None:
                # continuous, exactly as the merge applies it -- so the running
                # composite the NEXT tier is fitted against is the same
                # function, not a staircase version of it
                ang = np.arange(NA) / NA * 2 * np.pi
                kf = np.full(NA, ck[0]); qf = np.full(NA, cq[0])
                for o in range(1, _LDIC_ORDER + 1):
                    c, sn = np.cos(o * ang), np.sin(o * ang)
                    kf += ck[2 * o - 1] * c + ck[2 * o] * sn
                    qf += cq[2 * o - 1] * c + cq[2 * o] * sn
            else:
                kf = np.repeat(kk, per)[:NA]; qf = np.repeat(qq, per)[:NA]
            adj = f * kf[None, :] + qf[None, :]
        good = np.isfinite(adj) & (w > 0)
        C += np.where(good, adj * w, 0.0); Wt += np.where(good, w, 0.0)
    return out


def _bayer_planes(cfa):
    """The three colour planes of an RGGB mosaic at half resolution.

    R is [0::2, 0::2], B is [1::2, 1::2], G is the mean of the two greens. The
    result is the same shape as `stacks_half` and shares its geometry, so the
    ring grids and masks built for the luminance fit index all four identically.
    Each plane is offset from the quad centre by half a photosite; at ring-median
    scale that is far below the measurement.

    Built on demand and thrown away, never held for every tier at once: three
    half-res planes are 0.75x the size of the full-res mosaic they come from, and
    a long bracket at full frame already keeps gigabytes of those.

    CROPPED TO EVEN DIMENSIONS FIRST. A stacked mosaic is not guaranteed to have
    an even number of rows and columns -- the alignment shift trims the border it
    vacated, and it can trim an odd number. On a 4015-row mosaic cfa[0::2] has
    2008 rows and cfa[1::2] has 2007, so averaging the two greens raises
    "operands could not be broadcast together", which is exactly how this failed
    on a 250 mm Canon bracket: (2008,3010) against (2007,3010). Dropping the last
    row and column costs a half-pixel strip at one edge and makes all four quads
    the same shape.
    """
    h = (cfa.shape[0] // 2) * 2
    w = (cfa.shape[1] // 2) * 2
    cfa = cfa[:h, :w]
    return {"R": cfa[0::2, 0::2],
            "G": 0.5 * (cfa[0::2, 1::2] + cfa[1::2, 0::2]),
            "B": cfa[1::2, 1::2]}


_CHANNELS = ("R", "G", "B")


def _linfit_mad(x, y, nmin=200, nmax=100000, rng=None):
    """Straight-line fit by MINIMISING MEAN ABSOLUTE DEVIATION. y = a + b*x.

    This is what PixInsight's LinearFit actually does, and it is NOT least
    squares. From the PCL class reference:

        "Robust straight line fitting by minimization of mean absolute
         deviation."

    and of the intercept:

        "a ... represents a constant additive pedestal present in the whole
         dataset"

    Its outlier resistance comes from the L1 objective itself rather than from
    any rejection step; the class documents no rejection parameters at all. The
    LinearFit *process* adds a plain value-range filter (its reject-low and
    reject-high) BEFORE this fit, which is a different mechanism again -- see
    the caller.

    WHY IT MATTERS THAT THIS IS NOT LEAST SQUARES. 0.22.43 implemented the fit
    as least squares with sigma clipping, on the reasoning that "fit, throw out
    outliers, refit" is the same idea. It is not the same estimator: L1 puts the
    line through the bulk of the data and lets distant points stay distant,
    while clipped L2 decides membership first and then fits the survivors as if
    they were all there was. On a corona bracket the "outliers" are stars,
    prominences and residual alignment error along steep gradients -- real
    structure, not noise -- and the two estimators treat them differently.

    Algorithm is the standard one for this objective (Numerical Recipes 15.7,
    `medfit`): for a given slope the optimal intercept is the MEDIAN of
    (y - b*x), so the problem reduces to one dimension, and the derivative
    sum(x*sign(residual)) is monotonic in b and bracketed and bisected from the
    least-squares solution.

    CHECKED AGAINST THE REAL THING. PixInsight's LinearFit was run on two real
    frames of the reference bracket (1/13 s as reference, 1/8 s as target) and
    this function on the same two files decoded to sensor channels:

        channel   PixInsight    ours     difference
        G          0.625255   0.624724     -0.08%

    Green is the honest comparison, because it is the one channel PixInsight's
    RAW loading leaves closest to the sensor. Agreement to 0.08% on 44 million
    pixels says the ESTIMATOR is right.

    Red and blue did not match, and the reason is not the fit. Its RAW loader
    was set to VNG interpolation with camera white balance, black point
    correction and highlight clipping (read off the preferences dialog), and
    those steps CREATE the colour difference. Same two frames, same fit, same
    0.92 window, varying only what the fit was given:

        sensor channels, no demosaic          R/G 0.9988   B/G 0.9997
        VNG demosaic, no white balance        R/G 1.0069   B/G 0.9819
        VNG demosaic + camera white balance   R/G 1.0333   B/G 0.9835
        PixInsight (VNG + camera WB)          R/G 1.0181   B/G 0.9721

    The sensor channels show no colour difference at all. Demosaicing adds
    some, white balance adds more, and PixInsight lands in the same region --
    not on our exact figure, because its normalisation differs from dcraw's,
    but the mechanism is not in doubt. VNG is gradient-directed, so it
    interpolates R and B differently in two frames whose gradients differ, and
    the two frames of a bracket differ exactly there: around the limb, where
    the longer one is closer to clipping.

    THE CONSEQUENCE FOR WHERE WE APPLY IT. Because the reported per-channel
    slopes depend on the preprocessing, a fit run after a demosaic and a white
    balance is partly measuring those rather than a difference between the two
    exposures. Anyone comparing our per-channel numbers against a LinearFit run
    in PixInsight will find ours smaller, and this table is the reason. `_hill_chain` therefore fits in the sensor
    channel domain, before white balance, where the three channels are
    independent measurements of the same photons. That is a deliberate choice
    and this is the evidence for it.

    Returns (b, a, n_used) -- slope first, to match _linfit_clipped's order.
    """
    m = np.isfinite(x) & np.isfinite(y)
    n = int(m.sum())
    if n < nmin:
        return np.nan, np.nan, n
    X, Y = x[m], y[m]
    if n > nmax:                      # the fit does not need more than this
        r = rng if rng is not None else np.random.default_rng(0)
        k = r.choice(n, nmax, replace=False)
        X, Y = X[k], Y[k]
    N = X.size
    sx = X.sum(); sy = Y.sum(); sxx = float(X @ X); sxy = float(X @ Y)
    den = N * sxx - sx * sx
    if not np.isfinite(den) or den == 0:
        return np.nan, np.nan, n
    bb = (N * sxy - sx * sy) / den                     # least-squares start
    aa = (sxx * sy - sx * sxy) / den
    chi2 = float(np.sum((Y - (aa + bb * X)) ** 2))
    sigb = np.sqrt(chi2 / den) if chi2 > 0 else 0.0

    def rof(b):
        a = float(np.median(Y - b * X))
        d = Y - (a + b * X)
        return float(np.sum(X * np.sign(d))), a

    b1 = bb
    f1, a1 = rof(b1)
    if sigb <= 0 or not np.isfinite(f1):
        return float(bb), float(aa), n
    b2 = bb + np.sign(f1) * 3.0 * sigb
    f2, a2 = rof(b2)
    guard = 0
    while f1 * f2 > 0 and guard < 60:                  # bracket the root
        bb = b2 + 1.6 * (b2 - b1)
        b1, f1 = b2, f2
        b2 = bb
        f2, a2 = rof(b2)
        guard += 1
    if f1 * f2 > 0:
        return float(b1), float(a1), n
    tol = max(abs(sigb) * 0.001, 1e-12)
    a = a1
    while abs(b2 - b1) > tol:                          # bisect
        bb = b1 + 0.5 * (b2 - b1)
        if bb == b1 or bb == b2:
            break
        f, a = rof(bb)
        if f * f1 >= 0:
            f1, b1 = f, bb
        else:
            f2, b2 = f, bb
    b = 0.5 * (b1 + b2)
    a = rof(b)[1]
    return float(b), float(a), n


def _linfit_clipped(x, y, iters=4, clip=3.0, nmin=200):
    """Least squares y = k*x + q with sigma-clipping. PixInsight's LinearFit.

    Both parameters at once, and that is the whole point: a slope and an offset
    fitted one after the other are degenerate -- whichever is fitted first
    absorbs the other's error, in either direction. See the caller.

    The clipping is on the residual, re-estimated from the median absolute
    deviation each round, which is what makes this survive a prominence, a hot
    pixel or a star in the comparison region. Returns (k, q, n_used); k is NaN
    if there was never enough to fit.
    """
    m = np.isfinite(x) & np.isfinite(y)
    k = q = np.nan
    for _ in range(iters + 1):
        nu = int(m.sum())
        if nu < nmin:
            return np.nan, np.nan, nu
        xm, ym = x[m], y[m]
        A = np.empty((nu, 2), np.float64)
        A[:, 0] = xm
        A[:, 1] = 1.0
        try:
            (k, q), *_ = np.linalg.lstsq(A, ym, rcond=None)
        except np.linalg.LinAlgError:
            return np.nan, np.nan, nu
        del A, xm, ym
        r = y - (k * x + q)
        rm = float(np.median(r[m]))
        sd = 1.4826 * float(np.median(np.abs(r[m] - rm)))
        if not np.isfinite(sd) or sd <= 0:
            break
        nm = m & (np.abs(r - rm) < clip * sd)
        if int(nm.sum()) == nu:
            break
        m = nm
    return float(k), float(q), int(m.sum())


# LinearFit's value-range rejection, in fractions of full scale. The PROCESS
# applies these before the fit; the PCL class documents no rejection at all, so
# this is the process's own step and not part of the estimator.
#
# CONFIRMED against a real PixInsight installation: the tool opens with
# Reject low 0.000000 and Reject high 0.920000. Neither documentation page
# states the defaults, so this came from the dialog itself.
#
# 0.92 is worth noting next to our own numbers rather than just copied: the
# merge weight already goes hard to zero above 0.97 of saturation and starts
# its knee at 0.87, so LinearFit's upper limit sits between the two. Three
# independent choices landing in the same band is mild evidence that the band
# is right, not that anyone copied anyone.
_HILL_REJECT_LOW = 0.0
_HILL_REJECT_HIGH = 0.92


def _hill_chain(bayer, sat_half, secs, links_good, cal, sat_level, pedestal,
                progress, stats):
    """Hill's photometric combine: a per-channel LINEAR FIT between neighbours.

    THE PUBLISHED METHOD, implemented as described rather than as re-derived.
    Described in Jonathan Hill's talk "Advanced Solar Eclipse Photography"
    (on YouTube), summarised here in our own words:

      - each exposure is fitted to its NEIGHBOUR, not to a distant reference,
        because adjacent frames overlap best;
      - the chain runs between NEIGHBOURS -- 2->1, 3->2, 4->3 -- because
        adjacent frames overlap best. His reference is the longest exposure and
        he accepts the propagation risk explicitly; we reference the middle
        tier instead, for a reason measured on a 14-tier bracket and set out
        at the chaining step below;
      - the fit is PixInsight `LinearFit`: a slope AND an offset, PER COLOUR
        CHANNEL;

    WHAT LinearFit ACTUALLY IS, from its own documentation rather than from an
    assumption about it. The PCL class reference calls it "robust straight line
    fitting by minimization of mean absolute deviation", describes the intercept
    as "a constant additive pedestal present in the whole dataset", and
    documents NO rejection parameters -- the robustness is the L1 objective
    itself. The process wrapped around it adds a plain value-range filter
    (reject-low / reject-high) before the fit. So there are two separate
    mechanisms, and 0.22.43 had neither of them right: it used least squares
    with sigma clipping for both. Corrected in 0.22.44; see `_linfit_mad`.
      - and the stated reason it matters: "you see a lot of HDR techniques, if
        they don't perform this step, they have ringing ... those ringing
        artifacts if you don't perform this step".

    WHAT WE HAD, AND WHY THIS IS NOT THE SAME. The shipped chain is a special
    case of this one with two constraints bolted on: the slope is a single
    scalar per tier fitted on LUMINANCE and shared by R, G and B, and the offset
    is ONE number shared by every tier as well as every channel. In the shipped
    form a tier's value is

        (raw - pedestal) / (sec * cal[sec])

    which is exactly `K * (raw/sec) + Q` with K constrained to 1/cal and Q
    constrained to -pedestal/(sec*cal). This function lets K and Q be free, per
    tier and per channel, which is what LinearFit does.

    ON THE MEASUREMENTS THAT SAID IT WOULD NOT HELP MUCH. Measured on a 14-tier
    600 mm set, the colour freedom buys nothing resolvable (the three channels
    agree to within the ~1% the measurement itself can distinguish) and the
    per-pair offset improves link agreement most on the LONG tiers, which build
    the outer corona rather than the near-limb band where the rings sit. That is
    an argument about size, not about correctness, and it is measured on one
    bracket with an estimator that has already been wrong twice in this file's
    history. The published method is the defensible default to offer, so it is
    offered, and the measurement stands beside it rather than in place of it.

    ONE PROPERTY OF THE METHOD, NOT A BUG IN THIS CODE. A chain of pairwise
    fits equalises every tier ONTO the reference; it cannot see an offset the
    reference itself carries, because nothing is left to compare it against.
    Verified on a synthetic bracket with +9/+4/-3 ADU injected into every tier:
    the fit recovers it on all tiers except the reference, which comes out at
    exactly 0 by construction. So a black-level residual COMMON to the whole
    bracket survives this mode, where the shipped shared pedestal removed it.
    That is the better failure of the two -- a common offset is a uniform error
    in the composite, while the per-tier part, which this does remove, is the
    one that arrives divided by the exposure time and shows up as tier
    disagreement in the outer field. Worth knowing on a bracket whose fitted
    pedestal is large: it was -0.46 ADU on the 600 mm set (nothing) and +2.87
    ADU on a 360 mm one (not nothing).

    THE DEFAULT SINCE 0.22.86, and a toolbar setting -- "Photometry", with
    "Single scale factor" for the old behaviour. Opt-in through a marker file
    from 0.22.68 until then, which is how it should have started and not how it
    should have stayed: an invisible setting cost five of the test set's test runs,
    reported and discussed as linear-fit results on folders that had no marker
    in them. ECLIPSEFORGE_PHOTOMETRY still overrides the toolbar in either
    direction (`linfit` or `scalar`) for a bisect, and the marker files are
    still honoured.

    Returns ({sec: K[3]}, {sec: Q[3]}) in the convention

        value_on_common_scale = K[sec] * (raw_channel / sec) + Q[sec]

    normalised so the middle tier's mean slope is 1.0, which is the scale the
    rest of the pipeline already expects (cal[mid] == 1). Returns (None, None)
    if the chain could not be built.
    """
    n = len(secs)
    if n < 2:
        return None, None
    mid = n // 2
    rs = np.random.default_rng(0)
    k_link = np.ones((n - 1, 3))
    q_link = np.zeros((n - 1, 3))
    # Which links carried NO usable overlap at all. A link like that tells you
    # nothing about colour, and the chaining below must not carry a colour
    # correction across it -- see the note beside _blind. (0.22.79)
    blind = np.zeros(n - 1, bool)
    nfell = 0
    nqcap = 0
    modes = []
    _prev = None
    for i in range(n - 1):
        pa = _prev if _prev is not None else _bayer_planes(bayer[secs[i]])
        pb = _bayer_planes(bayer[secs[i + 1]])
        _prev = pb
        _qbad = False
        g = links_good[i]
        if g is not None:
            _h = min(g.shape[0], pa["R"].shape[0])
            _w = min(g.shape[1], pa["R"].shape[1])
            g = g[:_h, :_w]
        if g is None or int(g.sum()) < 1000:
            nfell += 1
            blind[i] = True
            del pa
            continue
        _ng = int(g.sum())
        sub = (rs.choice(_ng, 300000, replace=False) if _ng > 300000
               else slice(None))
        _h, _w = g.shape
        # LinearFit's reject-low / reject-high: a plain value-range filter on
        # BOTH images, applied before the fit, in fractions of full scale. Not
        # an outlier test -- the estimator does its own robustness -- so it uses
        # the raw level, and the same pixels for all three channels.
        _lvA = np.zeros_like(pa["G"][:_h, :_w][g][sub])
        _lvB = np.zeros_like(_lvA)
        for c in _CHANNELS:
            _wt = 2.0 if c == "G" else 1.0
            _lvA += _wt * pa[c][:_h, :_w][g][sub]
            _lvB += _wt * pb[c][:_h, :_w][g][sub]
        _sl = max(float(sat_level), 1e-9)
        _keep = ((_lvA / 4.0 / _sl > _HILL_REJECT_LOW)
                 & (_lvA / 4.0 / _sl < _HILL_REJECT_HIGH)
                 & (_lvB / 4.0 / _sl > _HILL_REJECT_LOW)
                 & (_lvB / 4.0 / _sl < _HILL_REJECT_HIGH))
        # HOW MUCH THE LINK CAN ACTUALLY CONSTRAIN. Two free parameters need
        # pixels where BOTH exposures carry real signal, not merely pixels the
        # mask let through. On a 14-tier bracket the short end has almost none:
        # measured on the reference set, the 1/4000->1/2000 link had ZERO
        # pixels above 2% of saturation in both tiers and the 1/2000->1/500
        # link had 23. Fitted anyway, the slope and the offset traded against
        # each other and returned a per-tier factor of 2.529 where the scalar
        # chain read 1.273, and the tiers then agreed three times WORSE at the
        # limb (variance 0.074 -> 0.206).
        #
        # So the model is reduced to what the data supports:
        #   plenty of overlap  -> slope and offset, the full LinearFit
        #   thin overlap       -> slope only; an offset there is fitted noise
        #   almost none        -> the exposure ratio, as the scalar chain does
        _both = ((_lvA / 4.0 / _sl > 0.02) & (_lvB / 4.0 / _sl > 0.02))
        _nboth = int(_both.sum())
        del _lvA, _lvB
        if int(_keep.sum()) < 1000:
            _keep = slice(None)
        _mode = ("full" if _nboth >= 5000
                 else "slope" if _nboth >= 500 else "ratio")
        modes.append(_mode)
        if _mode == "ratio":
            nfell += 3
            del pa
            continue
        for j, c in enumerate(_CHANNELS):
            # The shared pedestal comes off BEFORE the pairwise fit. The two
            # corrections do different jobs and neither replaces the other: the
            # pedestal removes the offset COMMON to the whole bracket, which a
            # chain of pairwise fits is blind to by construction, and the fitted
            # q removes what DIFFERS between two exposures. 0.22.43 switched the
            # shared one off on the assumption that Hill's offsets replaced it;
            # they do not, and the leftover then reached the sky-gradient step
            # and was removed as if it were sky (its blue span went 1.41x ->
            # 1.83x on the reference bracket).
            x = (pa[c][:_h, :_w][g][sub][_keep] - pedestal) / secs[i]
            y = (pb[c][:_h, :_w][g][sub][_keep] - pedestal) / secs[i + 1]
            # MAD minimisation, as PixInsight's LinearFit does it -- not least
            # squares with clipping, which is what 0.22.43 used.
            k, q, _nk = _linfit_mad(x, y, rng=rs)
            if _mode == "slope":
                # Refit through the origin once the offset is dropped, rather
                # than keeping a slope that was fitted alongside one.
                _m = np.isfinite(x) & np.isfinite(y)
                k = (float(np.median(y[_m] / np.maximum(x[_m], 1e-9)))
                     if _m.any() else np.nan)
                q = 0.0
            # Same implausibility guard the luminance chain uses. A link a
            # scalar fit would have thrown away must not enter this one either;
            # falling back to the exposure ratio means k=1, q=0 in these
            # exposure-normalised units.
            #
            # AND A GUARD ON THE OFFSET (0.22.80), which this never had. The
            # slope was checked and q was taken on trust, so a link whose fit
            # ran away in the offset went straight into the chain. the test set's
            # 2024 560 mm set: every tier with real overlap needed |q| of 50-260
            # in these units, and the 0.5s->1s link -- the last one carrying any
            # signal at all -- returned 1149, 1171 and 899. Those propagate to
            # every longer tier, i.e. the whole outer corona.
            #
            # The scale to judge q against is the signal the link was fitted on.
            # An offset worth a large fraction of the median level is not a
            # black-level residual, which is what this term exists to model; it
            # is the fit trading slope against offset on data that cannot
            # constrain both. Falling back to slope-only keeps what the link
            # does know -- the ratio -- and discards what it does not.
            # Judged PER LINK, not per channel: an offset kept on two channels
            # and dropped on the third would itself be a colour error, so the
            # verdict is recorded here and applied to all three below.
            _lvl = float(np.median(np.abs(y[np.isfinite(y)]))) if y.size else 0.0
            if (np.isfinite(q) and _lvl > 0
                    and abs(q) > _OFFSET_MAX_FRAC * _lvl):
                _qbad = True
            if not np.isfinite(k) or k <= 0 or abs(np.log(k)) > np.log(2.5) \
                    or not np.isfinite(q):
                nfell += 1
            else:
                k_link[i, j], q_link[i, j] = k, q
            del x, y
        if _qbad:
            # Drop the offset for the whole link and keep the slopes. The link
            # still knows the RATIO; what it does not know is a black-level
            # residual, and it has just claimed one worth a fifth of the signal.
            q_link[i, :] = 0.0
            modes[-1] = "slope"
            nqcap += 1
        del pa
    del _prev

    # CHAIN FROM THE MIDDLE TIER, NOT THE LONGEST, AND HERE IS WHY.
    #
    # Hill chains from the longest exposure and says plainly that he accepts the
    # risk of propagating an error down the chain. That is a fair trade on HIS
    # bracket, which is four or five frames: nothing is more than a few links
    # from the reference. This app is routinely handed fourteen tiers spanning
    # 12.6 EV, and then the shortest exposure sits THIRTEEN links from the
    # longest -- every one of them contributing its own error, and the last two
    # of them fitted on almost no overlapping signal.
    #
    # Referencing the middle tier halves the worst case: no tier is more than
    # n/2 links away. It is the same choice the scalar chain already makes, for
    # the same reason, and it changes nothing about the fit itself -- only which
    # tier is declared to be 1.000.
    #
    # Measured: chaining from the longest on the reference bracket moved the
    # 1/4000 s tier from 1.273 to 2.529 and tripled the tier-to-tier
    # disagreement at the limb.
    #
    # Composing T_j = T_{j+1} o (k_j, q_j) gives K_j = K_{j+1} k_j and
    # Q_j = K_{j+1} q_j + Q_{j+1}; the same composition runs the other way for
    # tiers above the reference.
    K = np.ones((n, 3))
    Q = np.zeros((n, 3))
    # NO COLOUR ACROSS A BLIND LINK (0.22.79).
    #
    # A failed link leaves k=1, q=0, so the composition below hands the tier its
    # neighbour's ACCUMULATED K and Q unchanged -- including the per-channel
    # part, which was measured on the neighbour's data and says nothing about
    # this tier. Found on the 560 mm 2024 test set, where the 1s->2s and
    # 2s->4s links had literally zero pixels carrying signal in both exposures:
    # the 1s tier's fit (red 9.7% above green, plus its offset) was inherited
    # verbatim by 2s and 4s, which are the tiers that carry the OUTER corona.
    # The merged red channel crossed zero at 3 R and ran to -1234 by 4 R, i.e.
    # 100% of pixels negative, and printed as a cyan ring.
    #
    # Scale may propagate across such a link -- the exposure ratio is still the
    # best estimate of it -- but colour may not. So the accumulated transform is
    # made ACHROMATIC at a blind link: all three channels take the mean, and the
    # offset is dropped, which is exactly what the scalar chain would do.
    #
    # A bracket whose links all fit is untouched: `blind` is all False and this
    # is the arithmetic it always was. Verified on the 600 mm reference set, where
    # every tier's K is distinct and no channel is negative anywhere.
    _flat = []
    for j in range(mid - 1, -1, -1):            # shorter than the reference
        K[j] = K[j + 1] * k_link[j]
        Q[j] = K[j + 1] * q_link[j] + Q[j + 1]
        if blind[j]:
            K[j] = np.full(3, float(np.mean(K[j])))
            Q[j] = np.zeros(3)
            _flat.append(secs[j])
    for j in range(mid + 1, n):                 # longer than the reference
        # invert the link: if y = k x + q maps j-1 onto j, then j onto j-1 is
        # x = (y - q)/k, so K_j = K_{j-1}/k and Q_j = Q_{j-1} - K_{j-1} q/k
        _k = np.where(np.abs(k_link[j - 1]) > 1e-9, k_link[j - 1], 1.0)
        K[j] = K[j - 1] / _k
        Q[j] = Q[j - 1] - K[j - 1] * q_link[j - 1] / _k
        if blind[j - 1]:
            # ACHROMATIC, AND Q DROPS TO ZERO (0.22.79, restored in .82).
            #
            # 0.22.80 made Q carry the neighbour's level instead, reasoning that
            # a blind link teaches us nothing so the neighbour's transform is
            # the honest one. On this bracket that carried the 1s tier's runaway
            # offset of about -1150 into every longer tier and the merge came
            # back 11% negative -- worse than the two thin rings it was meant to
            # remove. Measured on real data, zero is the better of the two.
            #
            # The step that leaves at the 1s/2s crossover is real and is what
            # those two thin rings are. Its CAUSE is the 1s tier's own offset,
            # which is five to eight times its neighbours'; the fix belongs
            # there, not here.
            K[j] = np.full(3, float(np.mean(K[j])))
            Q[j] = np.zeros(3)
            _flat.append(secs[j])

    # Put the result on the scale the rest of the pipeline expects. The shipped
    # convention is cal[mid] == 1, i.e. a mean slope of 1 at the middle tier;
    # Hill's chain is referenced to the longest exposure instead, so the two
    # differ by one constant. Normalising by the mean over the three channels
    # (not per channel) keeps whatever colour the fit found -- discarding it
    # here would silently reimpose the constraint this function exists to lift.
    _nrm = float(np.mean(K[mid]))
    if not np.isfinite(_nrm) or _nrm <= 0:
        progress.log("Hill photometric chain: normalisation failed — falling "
                     "back to the scalar chain", None)
        return None, None
    K /= _nrm
    Q /= _nrm
    if not (np.isfinite(K).all() and np.isfinite(Q).all()):
        progress.log("Hill photometric chain: non-finite terms — falling back "
                     "to the scalar chain", None)
        return None, None

    # 0.22.81 TRIED A NON-NEGATIVITY CONSTRAINT HERE AND IT WAS WRONG.
    #
    # The idea was sound -- K*x + Q maps a radiance to a radiance, so it must
    # not go negative -- but x_lo was taken as the 1st percentile of the whole
    # frame. That includes the occulted disc and the sky, which after pedestal
    # subtraction are noise centred on ZERO, so x_lo came out negative and the
    # constraint demanded a large POSITIVE offset: +41435 on the 1/2000s tier,
    # which is the one carrying the inner corona. The 560 mm render came back
    # entirely red.
    #
    # A constraint like this is only meaningful over the range where the tier
    # actually CONTRIBUTES to the merge -- where its saturation weight is
    # non-zero -- not over every pixel in its frame. Recorded here rather than
    # re-attempted, because three consecutive attempts to patch this from the
    # chain's downstream end each traded one artifact for another.

    Kd = {s: K[i].astype(np.float32) for i, s in enumerate(secs)}
    Qd = {s: Q[i].astype(np.float32) for i, s in enumerate(secs)}
    # How far this moves each tier against the scalar chain it replaces. 1/K is
    # the equivalent of cal[], so this is a like-for-like comparison.
    _dev, _worst = 0.0, None
    for i, s in enumerate(secs):
        c_equiv = 1.0 / max(float(np.mean(K[i])), 1e-9)
        d = abs(c_equiv / max(cal[s], 1e-9) - 1.0)
        if d > _dev:
            _dev, _worst = d, s
    _col = float(np.max(np.abs(K[mid] / max(float(np.mean(K[mid])), 1e-9) - 1.0)))
    stats["hill_chain"] = {
        "K": {("%g" % s): [round(float(v), 5) for v in Kd[s]] for s in secs},
        "Q": {("%g" % s): [round(float(v), 4) for v in Qd[s]] for s in secs},
        "links_fallen_back": nfell,
        "link_model": {f"{_exp_name(secs[i])}->{_exp_name(secs[i + 1])}": m
                       for i, m in enumerate(modes)},
        "worst_dev_vs_scalar_pct": round(100 * _dev, 2),
        "worst_dev_tier": None if _worst is None else ("%g" % _worst),
        "colour_spread_at_mid_pct": round(100 * _col, 2)}
    progress.log(
        f"PHOTOMETRY: per-channel linear-fit chain. "
        f"Slope and offset fitted per colour channel between neighbouring "
        f"exposures, after PixInsight LinearFit as used in Hill's HDR eclipse "
        f"workflow. "
        f"This REPLACES the single luminance scale factor per tier and the one "
        f"shared pedestal. Largest change against the scalar chain: "
        f"{100 * _dev:.1f}% on the {_exp_name(float(_worst)) if _worst else 'n/a'} "
        f"tier. Colour spread at the reference tier {100 * _col:.1f}%. "
        f"{nfell} link-channel fit(s) fell back to the exposure ratio"
        + (", and %d had an offset too large for the signal it was fitted on "
           "and were refitted through the origin" % nqcap if nqcap else "")
        + (". No colour correction was carried across %d blind link(s), so "
           "the %s tier(s) keep the scalar chain's own factor on all three "
           "channels -- a link with no overlapping signal measures no colour"
           % (len(_flat), ", ".join(_exp_name(x) for x in sorted(set(_flat))))
           if _flat else "") + ". "
        f"Chained from the MIDDLE tier, not the longest: on a bracket this "
        f"deep the longest is {n - 1} links from the shortest and the error "
        f"accumulates over every one. Per-link model "
        f"({modes.count('full')} slope+offset, {modes.count('slope')} slope "
        f"only, {modes.count('ratio')} exposure ratio) — a link is only given "
        f"an offset where enough pixels carry signal in both exposures to "
        f"constrain one.", None)
    return Kd, Qd


def _per_channel_photometry(bayer, sat_half, secs, ref, links_good, cal,
                            pedestal, cym, cxm, Rseed, sat_level, progress,
                            stats, wd=None):
    """Does any of our three photometric corrections depend on COLOUR?

    WHY THIS EXISTS. Every photometric correction in this pipeline is fitted on
    luminance and applied identically to R, G and B:

        per-tier scale   cal[s]      one scalar, fitted on the half-res luma
        additive offset  pedestal    one number, shared by every tier
        azimuthal gain   LDIC k,q    fitted on luma, applied to all channels

    The published method this is measured against -- PixInsight's LinearFit, as
    used in Hill's HDR eclipse workflow -- fits a slope AND an offset PER
    CHANNEL, chained from the longest exposure. Ours is more general in space
    (azimuthal, 60 segments, per tier; his is global) and strictly less general
    in colour (one channel; his is three). If the tier-to-tier variation a real
    bracket carries has any colour in it -- thin cloud, changing airmass, the
    sky brightening through totality -- this pipeline cannot express it. It gets
    absorbed into the luminance fit as an average and comes back out as a colour
    error that varies with radius, which is what a corona reported as orange or
    green after merging would look like.

    So: fit all three corrections a second time, per channel, on the SAME pixels
    the luminance fit used, and report the spread. Nothing is corrected and no
    output changes. It is a number, so the decision about a per-channel
    photometric chain is made on four datasets instead of on a video.

    `bayer[s]` is one tier's stacked full-res RGGB mosaic in RAW ADU, before
    white balance -- the WB multiplier is one constant per channel and cancels
    exactly in a ratio between tiers, so the factors below are directly
    comparable across channels. `links_good[i]` is the boolean mask the
    luminance fit used for the link secs[i] -> secs[i+1]; reusing it is the
    point, because then the only thing that differs between the three channels
    is the data.

    Returns a dict for `stats`, or None if there was nothing to measure.
    """
    n = len(secs)
    out = {"scale": {}, "pedestal": {}, "ldic_k_spread_pct": {}}
    _dump = {} if wd else None

    # ONE SHAPE FOR EVERYTHING. The colour planes come from the full-res mosaic
    # and the masks from the half-res luminance stack, and the two need not agree
    # to the last row: an odd-sized mosaic loses its final row to _bayer_planes,
    # while sat_half kept it. Everything below is cropped to the smaller of the
    # two so a plane, a mask and a ring grid can always be indexed together.
    _Hh, _Wh = sat_half[secs[0]].shape
    _ph, _pw = _bayer_planes(bayer[secs[0]])["R"].shape
    _Hh, _Wh = min(_Hh, _ph), min(_Wh, _pw)

    def _fit(a):
        return a[:_Hh, :_Wh]

    # ---- ring samples, in ONE pass over the tiers --------------------------
    # The offset fit and the azimuthal fit both work on ring grids, so a tier's
    # planes are built once, sampled onto both grids, and dropped. What survives
    # the loop is a handful of (rings x angles) arrays per tier, which is
    # nothing.

    def _grid(r0, r1, step, nang, nmin):
        rmax = min(cym, cxm, _Hh - cym, _Wh - cxm) - 2.0
        rs = np.arange(r0 * Rseed, min(r1 * Rseed, rmax), step * Rseed)
        if rs.size < nmin:
            return None
        th = np.linspace(0, 2 * np.pi, nang, endpoint=False)
        return (np.clip((cym + rs[:, None] * np.sin(th)).astype(np.int32),
                        0, _Hh - 1),
                np.clip((cxm + rs[:, None] * np.cos(th)).astype(np.int32),
                        0, _Wh - 1))

    _gp = _grid(1.5, 3.5, 0.05, 720, 8)            # pedestal: the outer field
    _gl = _grid(1.00, 2.5, 0.02, 12 * _LDIC_SEGMENTS, 10)   # LDIC: 1.0-2.5 R
    _sl = max(float(sat_level), 1e-9)
    _pmed = {}          # sec -> {chan: ring medians}, NaN where the ring clips
    _lval = {}          # sec -> {chan: ring samples}
    _lw = {}            # sec -> LDIC weight, from luminance
    for s in secs:
        _pl = _bayer_planes(bayer[s])
        if _gp is not None:
            _ys, _xs = _gp
            _keep = (~sat_half[s][_ys, _xs]).mean(axis=1) > 0.98
            _pmed[s] = {c: np.where(_keep, np.median(_pl[c][_ys, _xs], axis=1),
                                    np.nan) for c in _CHANNELS}
        if _gl is not None:
            _ys, _xs = _gl
            _lval[s] = {c: _pl[c][_ys, _xs].astype(np.float64)
                        for c in _CHANNELS}
            # The weight is built from the LUMINANCE level, as in the shipped
            # LDIC fit, so all three channels are fitted on the same pixels for
            # the same reason and a channel is not selected differently just
            # because it is brighter.
            _fr = (_lval[s]["R"] + 2 * _lval[s]["G"] + _lval[s]["B"]) / 4.0 / _sl
            _u = np.clip((0.85 - _fr) / 0.10, 0, 1)
            _l = np.clip((_fr - 0.004) / 0.008, 0, 1)
            _w = (_u * _u * (3 - 2 * _u)) * (_l * _l * (3 - 2 * _l))
            _lw[s] = np.where(sat_half[s][_ys, _xs], 0.0, _w)
            del _fr, _u, _l, _w
        del _pl

    # ---- 1. the shared additive offset, per channel ------------------------
    # Same outer-field ring medians and the same fitter, one channel at a time.
    # A black-level residual that differs between channels is exactly the thing
    # a single shared pedestal cannot carry.
    #
    # REPORTED, NOT USED to correct anything below. On its own this estimator
    # cannot separate an offset from a gain: given a synthetic bracket carrying
    # a pure 3%-per-step COLOUR GAIN and no offset at all, it returns -17 and
    # +17 ADU of offset for R and B. The gain and the offset are degenerate
    # when each is fitted alone, in both directions, which is exactly why the
    # measurement that matters below fits the two together.
    try:
        if _pmed:
            _sc = np.asarray([s * cal[s] for s in secs], float)
            for c in _CHANNELS:
                _pr = np.asarray([_pmed[s][c] for s in secs], float)
                _P, _Praw, _sig, _b0, _b1 = _fit_pedestal(
                    _pr, _sc, _PEDESTAL_MAX * float(sat_level))
                out["pedestal"][c] = {"applied": round(float(_P), 3),
                                      "fitted": round(float(_Praw), 3),
                                      "scatter_before": round(float(_b0), 5),
                                      "scatter_after": round(float(_b1), 5)}
    except Exception as _e:
        out["pedestal"] = {"error": str(_e)}
    # ---- 1b. the same offsets, with the exposure ladder free ---------------
    # This one IS used, on FITS. See _fit_channel_floors for why it does not
    # inherit the degeneracy that keeps 1. above to reporting.
    # THE FLOOR FIT NEEDS THE FAINT FIELD, AND `_gp` CANNOT REACH IT.
    #
    # `_grid` caps its radius at the nearest frame edge, so on a bracket where
    # the disc fills the frame it stops early -- 2.85 R on a third tester's set. At
    # that radius the corona still swamps a hundred ADU of floor, the fit is
    # ill-conditioned, and with one free exposure per tier it can absorb the
    # mismatch into the ladder instead: the first real run returned an
    # effective exposure of -225% of the stated time, which is a negative
    # exposure, and reported floors that bore no relation to the measured ones.
    #
    # The corners of a frame sit further out than the edge midpoints, so a
    # MASK reaches where a ring grid cannot. Beyond 3 R the corona is faint
    # enough that a floor is a real fraction of the signal, which is the only
    # regime in which these two parameters separate.
    _fmed = {}
    try:
        _yy = np.arange(_Hh, dtype=np.float32)[:, None] - cym
        _xx = np.arange(_Wh, dtype=np.float32)[None, :] - cxm
        _rr = np.hypot(_yy, _xx)
        _fm = (_rr > 3.0 * Rseed)
        if int(_fm.sum()) >= 5000:
            for s in secs:
                _pl = _bayer_planes(bayer[s])
                _ok = _fm & (~sat_half[s])
                if int(_ok.sum()) < 5000:
                    continue
                _fmed[s] = [float(np.median(_pl[c][_ok])) for c in _CHANNELS]
                del _pl
    except Exception:
        _fmed = {}

    try:
        if _fmed:
            _scal = {s: _fmed[s] for s in secs if s in _fmed}
            _scal = {s: v for s, v in _scal.items()
                     if np.isfinite(np.asarray(v, float)).all()}
            # Hand it the MEASURED effective exposure, not the stated one. See
            # the note in _fit_channel_floors: with the ladder free, a bracket
            # whose headers are wrong sends the fit into the ladder instead of
            # the floors, and the guard then refuses the whole thing.
            # ONLY THE TIERS WHERE A FLOOR IS STILL PART OF THE SIGNAL.
            #
            # The model is one rate plus one floor per channel. It breaks on the
            # long end of a third tester's set, where the corner level rises ~2.0x
            # per 1.5x step -- the sky brightened through the end of totality,
            # so a constant rate cannot describe it, and with all 18 tiers in
            # the residual is 3.5% and the fit is refused. Those tiers carry no
            # information about a floor anyway: at 15 s the floor is 0.5% of the
            # corner level.
            #
            # So keep a tier while the floor is at least about a tenth of what
            # it reads. Measured on that set, residual against the cut:
            #
            #     tiers kept      6      7      8     10     11     18
            #     residual    2.53%  1.76%  1.05%  0.25%  0.16%  3.50%
            #
            _lv = {s: float(np.mean(_scal[s])) for s in _scal}
            _lo = min(_lv.values()) if _lv else 0.0
            _keep = [s for s in secs if s in _scal
                     and _lv[s] <= _CH_MAX_LEVEL_RATIO * max(_lo, 1e-9)]
            if len(_keep) >= _CH_MIN_TIERS:
                _scal = {s: _scal[s] for s in _keep}
            _sv = {s: float(s) * float(cal.get(s, 1.0) or 1.0)
                   for s in secs if s in _scal}
            _d, _dinfo = _fit_channel_floors(
                _scal, [s for s in secs if s in _scal], exposures=_sv)
            out["channel_floors"] = _dinfo
            out["channel_floor_deltas"] = [float(x) for x in _d]
    except Exception as _e:
        out["channel_floors"] = {"error": str(_e)}
    del _pmed
    _pv = [v["applied"] for v in out["pedestal"].values()
           if isinstance(v, dict) and np.isfinite(v.get("applied", np.nan))]
    _ped_spread = (max(_pv) - min(_pv)) if len(_pv) == 3 else None
    # ---- 2. slope AND offset together, per link, per channel ---------------
    # This is the measurement the rest of the check rests on, and it is a
    # deliberate copy of the published method being compared against: a linear
    # fit of each exposure against its neighbour, with outlier rejection,
    # yielding a slope and an offset -- PixInsight's LinearFit, as used in the
    # HDR eclipse workflow this pipeline is measured against.
    #
    # JOINTLY, because separately does not work. A synthetic bracket carrying a
    # pure per-channel OFFSET reads as a 39% colour gain if the offset is not
    # removed first; a synthetic carrying a pure per-channel GAIN makes the
    # offset estimator invent +-17 ADU. Fitting one and then the other picks up
    # whichever error was fitted first. A single fit over the link's full
    # dynamic range constrains both, and recovers each cleanly on the same two
    # synthetics.
    #
    # x and y are exposure-normalised, so a slope of 1.000 means the two tiers
    # agree once their shutter speeds are accounted for -- the same convention
    # as cal[].
    _rat = {c: [] for c in _CHANNELS}
    _off = {c: [] for c in _CHANNELS}
    _band_rat = {"lo": {c: [] for c in _CHANNELS},
                 "hi": {c: [] for c in _CHANNELS}}
    _prev = None
    _rs2 = np.random.default_rng(0)
    for i in range(n - 1):
        pa = _prev if _prev is not None else _bayer_planes(bayer[secs[i]])
        pb = _bayer_planes(bayer[secs[i + 1]])
        _prev = pb
        g = None if links_good[i] is None else _fit(links_good[i])
        if g is None or g.sum() < 1000:
            for c in _CHANNELS:
                _rat[c].append(0.0); _off[c].append(np.nan)
                _band_rat["lo"][c].append(np.nan)
                _band_rat["hi"][c].append(np.nan)
            continue
        # A link can select millions of pixels and the fit does not need them; a
        # fixed seed keeps the number reproducible from run to run. The
        # subsample indexes INTO the mask-compacted array rather than into the
        # full plane, because `plane[g]` already compacts in one step and gives
        # the same pixels in the same order for every channel -- whereas
        # ravel()ing a strided colour plane copies the whole thing first.
        _ng = int(g.sum())
        sub = (_rs2.choice(_ng, 300000, replace=False) if _ng > 300000
               else slice(None))
        # SIGNAL BANDS, to tell a real colour gain from sensor non-linearity.
        # A gentle roll-off near the top of the range is not the same size in
        # every channel, because the channels sit at different raw levels for
        # the same scene, so it can MANUFACTURE a colour tilt. A true gain
        # difference cannot know how bright the pixel is. So the same fit is
        # repeated on a faint band and a bright band, defined once from the
        # LONGER tier's luminance so all three channels use identical pixels:
        # if the tilt is the same in both bands it is a gain, and if it is much
        # larger in the bright band it is non-linearity.
        _lb = ((_fit(pb["R"]) + 2 * _fit(pb["G"]) + _fit(pb["B"]))[g][sub]
               / 4.0 / max(float(sat_level), 1e-9))
        _bands = {"lo": (_lb > 0.02) & (_lb < 0.20),
                  "hi": (_lb > 0.35) & (_lb < 0.85)}
        _keepxy = {}
        for c in _CHANNELS:
            x = _fit(pa[c])[g][sub] / secs[i]
            y = _fit(pb[c])[g][sub] / secs[i + 1]
            if _dump is not None:
                # A thinned copy of the actual pixel pairs this link was fitted
                # on. THE POINT: every question asked of this measurement so far
                # -- chain or per link, which signal band, which estimator, how
                # to reject outliers -- has needed another full run of the
                # pipeline to answer. The samples are what the estimator sees,
                # so keeping them lets the next estimator be tried on the same
                # real data without re-stacking anything.
                _st = max(1, x.size // 40000)
                _keepxy[f"{i}_{c}_x"] = x[::_st].astype(np.float32)
                _keepxy[f"{i}_{c}_y"] = y[::_st].astype(np.float32)
            for _bn, _bm in _bands.items():
                # NaN, never 0.0, when the band cannot be fitted. 0.0 is a log
                # ratio of exactly 1.000 -- "the channels agree perfectly" --
                # and a rejected fit reported that way is indistinguishable
                # from a clean result. It happened: two links of a real run
                # printed R/G and B/G of exactly 1.0000 in the bright band
                # because all three channels had been rejected there.
                _bk = np.nan
                if int(_bm.sum()) >= 2000:
                    _bk, _, _ = _linfit_clipped(x[_bm], y[_bm])
                _band_rat[_bn][c].append(
                    np.nan if (not np.isfinite(_bk) or _bk <= 0
                               or abs(np.log(_bk)) > np.log(2.5))
                    else float(np.log(_bk)))
            k, q, _nk = _linfit_clipped(x, y)
            # same implausibility guard as the luminance chain, so a link the
            # luminance fit would have thrown away cannot show up here as a
            # colour effect
            if not np.isfinite(k) or k <= 0 or abs(np.log(k)) > np.log(2.5):
                _rat[c].append(0.0); _off[c].append(np.nan)
            else:
                _rat[c].append(float(np.log(k)))
                # q comes out in exposure-normalised units, where it is
                # dominated by the exposure ratio and unreadable. If both tiers
                # carry the same raw offset P, then
                #     y = k*x + P*(1/s_b - k/s_a)
                # so dividing q by that bracket puts it back into RAW ADU, on
                # the same scale as the shared pedestal, and the two can be
                # compared directly. Verified on a synthetic: +6 ADU injected,
                # +6.00 recovered.
                _d = 1.0 / secs[i + 1] - k / secs[i]
                _off[c].append(float(q / _d) if abs(_d) > 1e-9 else np.nan)
            del x, y
        if _dump is not None:
            _dump.update(_keepxy)
            _dump[f"{i}_lum"] = _lb[::max(1, _lb.size // 40000)].astype(np.float32)
            _dump[f"{i}_secs"] = np.array([secs[i], secs[i + 1]], np.float64)
        del _keepxy
        del pa
    del _prev

    def _chain(rat):
        lf = np.zeros(n)
        for i in range(ref, n - 1):
            lf[i + 1] = lf[i] + rat[i]
        for i in range(ref, 0, -1):
            lf[i - 1] = lf[i] - rat[i - 1]
        return {("%g" % s): round(float(np.exp(lf[i])), 5)
                for i, s in enumerate(secs)}

    for c in _CHANNELS:
        out["scale"][c] = _chain(_rat[c])
    out["link_offset"] = {
        c: {f"{_exp_name(secs[i])}->{_exp_name(secs[i + 1])}":
            (None if not np.isfinite(_off[c][i]) else round(_off[c][i], 4))
            for i in range(n - 1)} for c in _CHANNELS}

    # The number that matters: how far apart are the three chains on each tier,
    # after taking out the fact that the reference tier is 1.000 in all three by
    # construction. Expressed against G, which is where the luminance sits.
    _tilt, _worst_scale, _worst_tier = {}, 0.0, None
    for s in secs:
        k = "%g" % s
        gg = out["scale"]["G"][k]
        if gg > 1e-9:
            _tilt[k] = {"R_over_G": round(out["scale"]["R"][k] / gg, 5),
                        "B_over_G": round(out["scale"]["B"][k] / gg, 5)}
            d = max(abs(_tilt[k]["R_over_G"] - 1.0),
                    abs(_tilt[k]["B_over_G"] - 1.0))
            if d > _worst_scale:
                _worst_scale, _worst_tier = d, k
    out["scale_tilt"] = _tilt

    # The band comparison. A link whose band could not be fitted contributes
    # nothing rather than a zero, so a missing link does not read as agreement:
    # the chain is built only over links all three channels could fit, and the
    # tilt is reported per ADJACENT LINK rather than chained, because a chain
    # with holes in it is not comparable between the two bands.
    _band_tilt = {}
    for _bn in ("lo", "hi"):
        rows = {}
        for i in range(n - 1):
            v = {c: _band_rat[_bn][c][i] for c in _CHANNELS}
            if not all(np.isfinite(x) for x in v.values()):
                continue
            # per-link ratio of the channel slopes, as a percentage
            rows[f"{_exp_name(secs[i])}->{_exp_name(secs[i + 1])}"] = {
                "R_over_G": round(float(np.exp(v["R"] - v["G"])), 5),
                "B_over_G": round(float(np.exp(v["B"] - v["G"])), 5)}
        _band_tilt[_bn] = rows
    out["scale_tilt_by_band"] = _band_tilt

    def _worst_band(rows):
        w = 0.0
        for r in rows.values():
            w = max(w, abs(r["R_over_G"] - 1.0), abs(r["B_over_G"] - 1.0))
        return w

    # Compare only on the links BOTH bands could fit, or the comparison is
    # between different sets of links and says nothing.
    _both = set(_band_tilt["lo"]) & set(_band_tilt["hi"])
    _wlo = _worst_band({k: v for k, v in _band_tilt["lo"].items() if k in _both})
    _whi = _worst_band({k: v for k, v in _band_tilt["hi"].items() if k in _both})
    out["band_links_compared"] = sorted(_both)
    out["worst_tilt_faint_band_pct"] = round(100 * _wlo, 2)
    out["worst_tilt_bright_band_pct"] = round(100 * _whi, 2)

    # ---- 3. the azimuthal gain, per channel --------------------------------
    # k(phi) is the part a scalar cannot reach. Fitting it per channel asks
    # whether the diffuse light that motivates it is grey.
    try:
        if _lval:
            _order = sorted(secs, reverse=True)     # longest exposure first
            _W = [_lw[s] for s in _order]
            for c in _CHANNELS:
                # the SHARED pedestal, exactly as the shipped LDIC fit uses it,
                # so this measures the azimuthal gain under the treatment the
                # pipeline actually applies. The fit runs over 1.0-2.5 R, where
                # a few ADU is parts per thousand of the signal, so the choice
                # does not move it.
                _vals = [(_lval[s][c] - pedestal) / (s * cal[s])
                         for s in _order]
                _f = _fit_azimuthal_affine(_vals, _W, _order)
                out["ldic_k_spread_pct"][c] = {
                    ("%g" % k): round(float(100 * (np.percentile(v[0], 90) -
                                                   np.percentile(v[0], 10))), 1)
                    for k, v in _f.items()}
                del _vals, _f
            del _W
    except Exception as _e:
        out["ldic_k_spread_pct"] = {"error": str(_e)}
    del _lval, _lw

    # ---- verdict -----------------------------------------------------------
    # 1% is not arbitrary: the shipped luminance chain's own links are trusted
    # to a few tenths of a percent where the tiers overlap well, and a 1% colour
    # tilt on a tier is the level at which it becomes visible as a cast in the
    # merged corona rather than as noise in the fit.
    out["worst_scale_tilt_pct"] = round(100 * _worst_scale, 2)
    out["worst_scale_tilt_tier"] = _worst_tier
    if _ped_spread is not None:
        out["pedestal_spread_adu"] = round(float(_ped_spread), 3)
    # the largest azimuthal-gain difference between channels on any tier
    _kd = None
    _kk = out.get("ldic_k_spread_pct") or {}
    if all(isinstance(_kk.get(c), dict) for c in _CHANNELS):
        for t in _kk["G"]:
            v = [_kk[c].get(t) for c in _CHANNELS]
            if all(x is not None for x in v):
                _kd = max(_kd or 0.0, max(v) - min(v))
        if _kd is not None:
            out["ldic_k_spread_diff_pct"] = round(float(_kd), 2)
    stats["per_channel_photometry"] = out

    if _dump:
        try:
            _dp = os.path.join(wd, "per_channel_samples.npz")
            np.savez_compressed(_dp, **_dump)
            out["samples_file"] = os.path.basename(_dp)
            out["samples_mb"] = round(os.path.getsize(_dp) / 1e6, 1)
        except Exception as _e:
            progress.log(f"per-channel samples not written ({_e})", None)
    _dump = None

    # The largest per-link offset any channel needed, in raw ADU. This is the
    # term the shipped pipeline has no per-channel form of at all: one pedestal,
    # shared by every tier and every channel.
    _od = None
    for c in _CHANNELS:
        for v in out["link_offset"][c].values():
            if v is not None:
                _od = max(_od or 0.0, abs(float(v)))
    if _od is not None:
        out["worst_link_offset"] = round(float(_od), 4)

    _pt = (f"; the outer-field offset fitted per channel spreads "
           f"{_ped_spread:.2f} ADU" if _ped_spread is not None else "")
    _ot = (f"; the largest offset any single link needed is {_od:.1f} ADU"
           if _od is not None else "")
    # The band comparison is the part that says WHAT the tilt is, so it goes in
    # the log rather than only in the bundle.
    _bt = ""
    if _both:
        _bt = (f". Over the {len(_both)} link(s) both signal bands could fit, "
               f"the colour difference is {100 * _wlo:.1f}% on faint pixels "
               f"(2-20% of saturation) and {100 * _whi:.1f}% on bright ones "
               f"(35-85%)")
        if _whi > max(2 * _wlo, _wlo + 0.005):
            _bt += (" — much larger where the signal is high, which is what "
                    "sensor non-linearity looks like rather than a gain "
                    "difference, so treat the headline figure as an upper "
                    "bound")
        elif _wlo > max(2 * _whi, _whi + 0.005):
            _bt += (" — larger on faint pixels, where the shorter tier of each "
                    "link is closest to its noise floor; some of that is the "
                    "fit being dragged toward zero by noise in the reference "
                    "rather than a real difference")
        else:
            _bt += (" — the same at both signal levels, which is what a true "
                    "per-channel gain looks like and not what non-linearity "
                    "would do")
    _kt = (f"; the azimuthal gain's swing around the limb differs between "
           f"channels by up to {_kd:.0f} points" if _kd is not None else "")
    if _worst_tier is not None and _worst_scale > 0.01:
        progress.log(
            f"per-channel photometry check: fitting slope and offset together "
            f"per link, the way PixInsight's LinearFit does, the "
            f"{_exp_name(float(_worst_tier))} tier's slope differs between "
            f"colour channels by up to {100 * _worst_scale:.1f}%{_pt}{_ot}{_kt}{_bt}. "
            f"This pipeline fits ONE factor per tier on luminance and applies "
            f"it to R, G and B alike, so a difference of this size cannot be "
            f"corrected and comes back out as a colour cast that changes with "
            f"radius. Recorded, not corrected — the per-channel factors are in "
            f"the diagnostics bundle.", None)
    elif _worst_tier is not None:
        progress.log(
            f"per-channel photometry check: fitting slope and offset together "
            f"per link, the three colour channels agree on every tier's slope "
            f"to {100 * _worst_scale:.1f}%{_pt}{_ot}{_kt}{_bt}. Fitting the "
            f"photometric chain on luminance alone costs this dataset nothing.",
            None)
    return out


def _pick_feather(stacks_half, sat_half, secs, cal, pedestal, sat_level,
                  cym, cxm, Rseed, alpha, sigma_half, progress,
                  advise_only=False):
    """Choose the merge feather from THIS dataset, not from a global default.

    0.22.25 put the plain feather back because it hides the ring artifact.
    That was right for the 600 mm reference set, where the leak it trades for costs
    25% of the true brightness at 1.02 R, and wrong for the 360 mm test set,
    where the same leak costs a factor of EIGHT and prints as a pink rim
    around the limb -- which is exactly what he reported on 0.22.26.

    The severity is a property of the bracket (how much weight a blown tier
    carries, how steep the limb gradient is), it is not knowable in advance,
    and it is cheap to measure: merge the half-resolution stacks both ways and
    compare the near-limb level. A few hundred milliseconds against a run of
    minutes.

    The rule: use the plain feather while its error stays small enough not to
    read as a rim, and fall back to the leak-free taper when it does not. The
    threshold is 25%, which is where the one dataset known to be acceptable
    sits -- so this ships the picture the tester already approved on their data and
    refuses the 8x version on the test set's.

    Returns (mode, ratio) with ratio = plain level / taper level at 1.02 R.
    """
    H, W = stacks_half[secs[0]].shape
    rmax = min(cym, cxm, H - cym, W - cxm) - 2.0
    rs = np.arange(1.00 * Rseed, min(1.30 * Rseed, rmax), 0.01 * Rseed)
    if rs.size < 6:
        return "plain", float("nan")
    th = np.linspace(0, 2 * np.pi, 360, endpoint=False)
    ys = np.clip((cym + rs[:, None] * np.sin(th)).astype(np.int32), 0, H - 1)
    xs = np.clip((cxm + rs[:, None] * np.cos(th)).astype(np.int32), 0, W - 1)
    acc = {"plain": None, "taper": None}
    wsum = {"plain": None, "taper": None}
    for s in secs:
        # `stacks_half` is a LUMA of the Bayer quad, after the flat -- not a
        # per-channel maximum in raw units. 0.22.28 compared it against
        # 0.87*sat_level anyway, so the knee almost never fired, both weights
        # came out nearly identical, and every dataset reported 95-96%
        # whatever its bracket actually did. The saturation mask is the thing
        # that is right at this resolution: `sat_half` is true where ANY
        # photosite in the quad hit the ceiling, measured before the flat.
        # So scale the knee by the luma at which THIS tier clips.
        lum_h = stacks_half[s].astype(np.float32)
        sat = sat_half[s]
        v = (~sat).astype(np.float32)
        _lsat = float(np.median(lum_h[sat])) if sat.any() else 0.0
        if _lsat <= 0:
            x = np.zeros_like(lum_h)          # nothing clips: no knee to apply
        else:
            x = lum_h / np.float32(_lsat)
        q = 0.5 * (1.0 + np.tanh((0.87 - x) / 0.06)).astype(np.float32)
        q[sat] = 0.0
        _lo = _MERGE_FLOOR
        q *= 0.5 * (1.0 + np.tanh((x - _lo) / (0.5 * _lo)))
        q = q.astype(np.float32)
        val = (lum_h - pedestal) / np.float32(s * cal[s])
        sa = np.float32(s ** alpha)
        ws = {"plain": sa * _feather_weight(q, v, sigma_half, "plain"),
              "taper": sa * _feather_weight(q, v, sigma_half, "taper")}
        for k, w in ws.items():
            if acc[k] is None:
                acc[k] = np.zeros_like(val); wsum[k] = np.zeros_like(val)
            acc[k] += w * val; wsum[k] += w
        del lum_h, sat, v, q, val, ws, x
    out = {}
    for k in acc:
        lum = acc[k] / np.maximum(wsum[k], 1e-9)
        out[k] = np.median(lum[ys, xs], axis=1)
    # Over a BAND, not one radius, and with a degeneracy guard. The leak can
    # only ever pull the merged level DOWN, so a ratio above 1 means the
    # comparison is meaningless -- which happens when every tier is clipped
    # across the band and the leak-free merge has no weight left to average.
    # A synthetic bracket blown 100% at the limb returned 8e17 before this.
    band = (rs >= 1.01 * Rseed) & (rs <= 1.15 * Rseed)
    with np.errstate(all="ignore"):
        rr = out["plain"][band] / np.where(out["taper"][band] > 0,
                                           out["taper"][band], np.nan)
    rr = rr[np.isfinite(rr)]
    # the MINIMUM over the band, not the median: a rim is visible where the
    # deficit is worst, and the deficit shrinks quickly with radius (on
    # the test set's 360 mm, 0.13 at 1.02 R but 0.49 by 1.06 R -- a median over the
    # band would blur the two datasets together).
    ratio = float(np.min(rr)) if rr.size >= 3 else float("nan")
    if not np.isfinite(ratio) or ratio > 1.2:
        progress.log("feather trial: no usable comparison at the limb on this "
                     "bracket (every tier clipped there, or no signal) — "
                     "keeping the plain feather", None)
        return "plain", float("nan")
    # THRESHOLD, AND WHAT IT IS CALIBRATED ON. Two real numbers, measured by
    # rebuilding both merges from the exported tiers:
    #     the reference set's 600 mm   0.747  -- he calls this good, no rim
    #     the test set's 360mm 0.126  -- he reported a pink rim on 0.22.26
    # Anything in (0.13, 0.75) separates them; 0.60 means "accept up to a 40%
    # deficit". Two calibration points that far apart do not pin it down, and
    # a third dataset landing between them is what would.
    mode = "plain" if ratio >= 0.60 else "taper"
    if mode == "plain":
        progress.log(f"feather trial: the plain feather reads "
                     f"{100 * ratio:.0f}% of the leak-free level at 1.02 R — "
                     f"small enough to keep, and it suppresses the ring "
                     f"artifact", None)
    else:
        progress.log(f"feather trial: the plain feather reads only "
                     f"{100 * ratio:.0f}% of the leak-free level at 1.02 R on "
                     f"this bracket — that prints as a rim around the limb, so "
                     f"the leak-free weight is used instead. The ring artifact "
                     f"will be visible; it is the smaller error here", None)
    return mode, ratio


def _feather_weight(w, valid, sigma, mode=None):
    """Smooth a merge weight so it neither leaks into, nor steps at, the edge of
    the region where the tier is clipped.

    THREE VERSIONS, TWO OF THEM WRONG, ALL THREE MEASURED. Two-tier merge
    rebuilt from the test set's 360 mm raws (1/125 s and 1/8 s, identical weights,
    only this function changed). "max weight step" is the largest jump one
    tier's weight makes between neighbouring pixels; the profile columns are the
    merged radial profile against the unfeathered reference:

        version                 peak at   1.00R  1.02R  1.06R  1.20R   step
        no feather (reference)   1.025 R  1.000  1.000  1.000  1.000  0.978
        plain blur (<=0.22.15)   1.095 R  3.690  0.674  0.856  0.977  0.027
        blur then mask (0.22.16) 1.025 R  0.999  1.000  1.000  1.000  0.986
        this one                 1.025 R  1.054  1.000  1.000  1.000  0.046

    A PLAIN BLUR goes both ways, so it hands a tier weight in the band where
    that tier is CLIPPED, and a clipped tier under-reports by definition.
    Weight goes as exposure time, so a blown 2 s tier outvotes an unblown
    1/1000 s one two thousand to one: the merged inner corona was dragged down
    over a band just outside the limb, with a bright overshoot at the limb
    itself. That pair was the "pink rim" reported on two datasets.

    NORMALIZED CONVOLUTION AND A HARD MASK fixed the profile exactly and broke
    something else. Restoring the weight's magnitude at the boundary and then
    zeroing it just past that boundary leaves a step of nearly the full weight
    -- 0.986 against the plain blur's 0.027. Every tier's saturation contour
    then prints as its own arc in the merged image. a tester saw it immediately on
    his own 600 mm set: *"strange rims in the corona, in almost all layers.
    Something broke - it was perfect with my data before."*
    It showed on his set and not on the test set's because his bracket runs at
    exposure exponent 0.55 (the alpha trial fires there), so every tier carries
    comparable weight and every tier's boundary is visible; at alpha 1.0 the
    longest unclipped tier dominates and the others' steps do not show.

    THE FIX IS TO TAPER TO ZERO AT THE BOUNDARY FROM INSIDE. `den`, the blurred
    validity mask, is 1 deep inside the usable region and 0.5 on the boundary,
    so a smoothstep on (den-0.5)/0.5 reaches zero exactly where the tier starts
    clipping and 1 where it is safely unclipped. The weight is then continuous
    everywhere, still exactly zero on every clipped pixel, and the merged
    profile is unchanged from 1.02 R outward.

    The cost is real and worth naming: a tier that clips loses a band about one
    sigma wide INSIDE its own clipping contour. That is the conservative
    direction -- a tier within a feather of its ceiling is the one you least
    want to average in -- and a tier that never clips is untouched, because its
    `den` is 1 everywhere.
    """
    if sigma <= 0:
        return np.asarray(w, np.float32)
    # `mode` explicit beats the environment. _pick_feather has to be able to
    # ask for BOTH forms in one process, and reading the variable here made it
    # silently compare a thing with itself -- the trial returned exactly 1.000
    # on brackets that were 100% blown at the limb.
    _mode = (mode or os.environ.get("ECLIPSEFORGE_FEATHER",
                                    _FEATHER_DEFAULT)).lower()
    if _mode == "plain":                       # <= 0.22.15, leaks into the clip
        return ndimage.gaussian_filter(np.asarray(w, np.float32),
                                       sigma).astype(np.float32)
    v = np.asarray(valid, np.float32)
    den = ndimage.gaussian_filter(v, sigma)
    num = ndimage.gaussian_filter(np.asarray(w, np.float32) * v, sigma)
    out = np.where(den > 1e-3, num / np.maximum(den, 1e-3), 0.0)
    if _mode == "masked":                      # 0.22.16, steps at the clip edge
        return (out * v).astype(np.float32)
    t = np.clip((den - 0.5) * 2.0, 0.0, 1.0)
    return (out * t * t * v).astype(np.float32)


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
    from .brand import OUTPUT_NAME
    out = os.path.join(folder, OUTPUT_NAME, "aligned_tiers")
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

    ECLIPSEFORGE_NO_SKYGRAD=1 skips the whole step. It exists because this is a
    PER-CHANNEL division extrapolated over the entire frame from an annulus
    beyond the corona, so any error in it lands as a colour cast on the inner
    corona, where nothing constrains the fit. On the 560 mm 2024 test set the
    fitted spans are R 1.141x G 1.104x B 1.054x -- RED steepest, the opposite of
    the Rayleigh ordering the paragraph above offers as the evidence that what
    is being fitted is atmosphere. The switch turns the question into one run.
    Diagnostic only: it changes nothing in the merge.
    """
    if os.environ.get("ECLIPSEFORGE_NO_SKYGRAD") == "1":
        stats["sky_gradient"] = {"applied": False,
                                 "skipped_by_env": True}
        progress.log("sky gradient: SKIPPED by ECLIPSEFORGE_NO_SKYGRAD=1 — "
                     "diagnostic run, the per-channel division is not applied",
                     None)
        return
    hp = os.path.join(wd, "hdr_rgb.npy")
    if not os.path.exists(hp):
        return
    hdr = load_big(hp)
    H, W, _ = hdr.shape
    d = 6
    # copy=True IS THE WHOLE POINT (0.22.76). hdr_rgb.npy is float32 and
    # np.asarray(view, np.float32) is then a no-copy VIEW of the memmap:
    # np.shares_memory(S, hdr) is True and S.base is still the memmap after the
    # `del` below, so the mapping stayed open through the np.save at the end of
    # this function and Windows refused the write -- the failure this comment
    # already described, with a fix that never took effect.
    S = np.array(hdr[::d, ::d], dtype=np.float32, copy=True)
    # RELEASE THE MAP BEFORE THE WRITE (Windows).
    #
    # `S` is a copy, so `hdr` has no reader left -- but on Windows an open
    # memory mapping LOCKS the file, and the np.save at the end of this function
    # then fails with [Errno 22] Invalid argument on the very path it is trying
    # to write. POSIX allows overwriting a mapped file, which is why this never
    # showed up on macOS. Reported from a Windows run as
    #
    #     sky gradient removal skipped ([Errno 22] Invalid argument:
    #     'C:\\Users\\...\\.eclipseforgehdr\\hdr_rgb.npy')
    #
    # An earlier report of the same error was blamed on OneDrive refusing to
    # memory-map, because the path had OneDrive in it. It was a local folder
    # that merely happened to be named that way, and `load_big`'s fallback --
    # written for that theory -- was already shipping when this one arrived. So
    # the failure was never in OPENING the map. It is in closing it.
    del hdr
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


def neutralise_corona_colour(wd, cy, cx, R, stats, progress):
    """White-balance a stack that arrived without a camera white balance.

    WHY IT IS NEEDED. `RawFile` gets `daylight_whitebalance` and a colour matrix
    from LibRaw. `FitsFrame` and `TiffFrame` have neither -- no FITS convention
    carries them -- so both fall back to wb (1,1,1) and cam2rgb = identity, and
    the merge comes out in RAW SENSOR RGB. A silicon sensor under a CFA is far
    more sensitive in green than in red, so raw-sensor white is not white: on a
    25-tier FITS bracket the corona measured R/((G+B)/2) = 0.56 where a
    colour-managed stack of the same corona reads 1.11, and the picture had a
    blue-cyan rim around the Moon on a brown sky. It also puts the prominence
    gate's reference colour in the wrong place, which is why that run flagged
    nothing.

    WHY THE CORONA IS A LEGITIMATE WHITE REFERENCE, and this is not a taste
    call: the K-corona is photospheric light Thomson-scattered off free
    electrons. Thomson scattering is wavelength-independent, so the inner
    corona has the SUN's spectrum -- it is the one thing in the frame that is
    white by physics rather than by convention. The dusty F-corona is redder and
    takes over further out, so the reference is taken close in, at 1.05-1.6 R,
    and the outward reddening the file actually contains is left alone.

    ONLY when there was no camera white balance to use. A raw bracket keeps the
    body's own numbers; this never runs for it.
    """
    hp = os.path.join(wd, "hdr_rgb.npy")
    if not os.path.exists(hp):
        return
    a = load_big(hp)
    H, W, _ = a.shape
    # Decimate against the DISC, not by a fixed factor. A fixed d=4 leaves a
    # 40 px moon with ~460 px in the reference annulus and the measurement is
    # abandoned; scaling with R keeps roughly 13k-20k px whether the disc is
    # 100 px across or 1000, which is the range these brackets actually span
    # (a 200 mm frame gives R=108 px, a 1000 mm one R=524).
    d = int(np.clip(round(R / 50.0), 1, 8))
    # copy=True, for the reason spelled out in remove_sky_gradient: on a float32
    # file np.asarray is a view and `del a` does not close the mapping. This one
    # happened to survive because `del S` precedes its np.load, which is luck
    # rather than design and would break on any reordering.
    S = np.array(a[::d, ::d], dtype=np.float32, copy=True)
    del a                                  # see remove_sky_gradient: Windows
    h, w, _ = S.shape
    yy = np.arange(h, dtype=np.float32)[:, None] - cy / d
    xx = np.arange(w, dtype=np.float32)[None, :] - cx / d
    r = np.hypot(yy, xx) / max(R / d, 1e-6)
    m = (r > 1.05) & (r < 1.60)
    if m.sum() < 2000:
        progress.log("corona white balance: too few pixels in the 1.05-1.6 R "
                     "annulus to measure — left in raw sensor colour", None)
        return
    med = np.median(S[m].reshape(-1, 3), axis=0).astype(np.float64)
    del S
    if not np.isfinite(med).all() or med.min() <= 0:
        progress.log("corona white balance: the annulus has no usable signal — "
                     "left in raw sensor colour", None)
        return
    lum = 0.2126 * med[0] + 0.7152 * med[1] + 0.0722 * med[2]
    g = (lum / med).astype(np.float32)     # unit-luminance, so brightness holds
    # A sensor's raw R:G:B never needs more than a few times' correction. More
    # than that is a broken channel or a mono frame read as colour, and dividing
    # by it would invent colour rather than correct it.
    if float(g.max() / max(g.min(), 1e-6)) > 8.0:
        progress.log(f"corona white balance: implausible gains "
                     f"{g[0]:.2f}/{g[1]:.2f}/{g[2]:.2f} — not applied", None)
        stats["corona_wb"] = {"gains": [float(x) for x in g], "applied": False}
        return
    out = np.load(hp)
    out *= g[None, None, :]
    np.save(hp, out)
    np.save(os.path.join(wd, "hdr_lum.npy"),
            (0.2126 * out[:, :, 0] + 0.7152 * out[:, :, 1]
             + 0.0722 * out[:, :, 2]).astype(np.float32))
    rgb_before = float(med[0] / max(0.5 * (med[1] + med[2]), 1e-9))
    stats["corona_wb"] = {"gains": [float(x) for x in g], "applied": True,
                          "r_over_gb_before": rgb_before}
    progress.log(f"corona white balance: no camera white balance in these "
                 f"files, so the inner corona (1.05-1.6 R) is used as the white "
                 f"reference — R/GB was {rgb_before:.2f}, gains R {g[0]:.3f} "
                 f"G {g[1]:.3f} B {g[2]:.3f}", None)


def measure_radial_colour(wd, cy, cx, R, stats, progress,
                          r0=1.02, r1=6.0, nbin=96, smooth=5):
    """The corona's colour as a function of radius, stored for the renderer.

    WHY ONE SET OF GAINS IS NOT ENOUGH. `neutralise_corona_colour` above takes
    the K-corona as a white reference, which is sound physics -- Thomson
    scattering is wavelength-independent, so the inner corona carries the Sun's
    own spectrum. It then applies that one measurement everywhere. That is only
    right if the colour of the frame is constant with radius, and on a
    low-altitude eclipse it is not. Measured on a third tester's set, per-channel
    floors removed, on two tiers twelve times apart in exposure:

        radius        R/G     B/G          R/G     B/G
                      (1/4 s tier)         (2.96 s tier)
        1.2-1.5 R    1.301   0.592            clipped
        1.5-2.0 R    1.196   0.710         1.204   0.702
        2.0-2.5 R    1.038   0.893         1.035   0.883
        2.5-3.0 R    0.923   1.009         0.924   0.994
        3.0-3.5 R    0.861   1.071         0.864   1.059
        3.5-4.2 R    0.803   1.130         0.812   1.113
        4.2-5.0 R    0.763   1.164         0.773   1.152

    The two tiers agree to 1% at every radius, so this is the scene as
    recorded -- not a tier, not the ladder, not a black level. R/B swings by 2x
    across the frame: a reddened corona near the limb giving way to a bluer sky
    further out. Neutralising on the inner annulus alone leaves the outer field
    blue, which is the orange-centre-to-grey-edge split these brackets show.

    WHAT THIS IS AND IS NOT. It is a CHOICE, in the same way the existing
    control is: the corona at 7 degrees altitude really is reddened, and
    removing that says you would rather see the K-corona's intrinsic white than
    what the atmosphere did to it. The physically-modelled answer is to fit the
    sky as a separate additive term with its own colour and subtract it; this
    is the cheaper one that flattens the symptom. It is deliberately measured
    and stored rather than applied, so the renderer can dial it live and the
    picture decides.

    Robustness: a MEDIAN per annulus, so a prominence or a bright streamer
    occupying part of one azimuth does not move it. Bins whose luminance falls
    into the noise are not measured, and the last good value is held outward
    rather than extrapolated.
    """
    hp = os.path.join(wd, "hdr_rgb.npy")
    if not os.path.exists(hp):
        return
    try:
        a = load_big(hp)
        H, W, _ = a.shape
        d = int(np.clip(round(R / 50.0), 1, 8))
        S = np.array(a[::d, ::d], dtype=np.float32, copy=True)
        del a
        h, w, _ = S.shape
        yy = np.arange(h, dtype=np.float32)[:, None] - cy / d
        xx = np.arange(w, dtype=np.float32)[None, :] - cx / d
        rr = np.hypot(yy, xx) / max(R / d, 1e-6)
        edges = np.linspace(r0, r1, nbin + 1)
        rg = np.full(nbin, np.nan)
        bg = np.full(nbin, np.nan)
        lv = np.full(nbin, np.nan)
        for i in range(nbin):
            m = (rr >= edges[i]) & (rr < edges[i + 1])
            n = int(m.sum())
            if n < 400:
                continue
            v = np.median(S[m].reshape(-1, 3), axis=0).astype(np.float64)
            if not np.isfinite(v).all() or v[1] <= 0:
                continue
            rg[i], bg[i], lv[i] = v[0] / v[1], v[2] / v[1], v[1]
        del S
        ok = np.isfinite(rg) & np.isfinite(bg)
        # drop the faint tail: once the annulus median stops falling with
        # radius it is sky noise, not corona, and its colour is meaningless
        if ok.sum() >= 8:
            _l = np.where(ok, lv, np.nan)
            _pk = np.nanmax(_l)
            ok &= np.isfinite(_l) & (_l > 1e-4 * _pk)
        if ok.sum() < 8:
            progress.log("radial colour: too few usable annuli to measure", None)
            return
        idx = np.flatnonzero(ok)
        # hold the ends, then smooth: a ratio profile should be smooth in r,
        # and any wiggle left in it would be imprinted on the picture as rings
        rgf = np.interp(np.arange(nbin), idx, rg[idx])
        bgf = np.interp(np.arange(nbin), idx, bg[idx])
        if smooth > 1:
            k = np.ones(int(smooth)) / float(int(smooth))
            pad = int(smooth)
            rgf = np.convolve(np.r_[[rgf[0]] * pad, rgf, [rgf[-1]] * pad],
                              k, "same")[pad:-pad]
            bgf = np.convolve(np.r_[[bgf[0]] * pad, bgf, [bgf[-1]] * pad],
                              k, "same")[pad:-pad]
        centres = 0.5 * (edges[:-1] + edges[1:])
        np.save(os.path.join(wd, "colour_radial.npy"),
                np.vstack([centres, rgf, bgf]).astype(np.float32))
        _sw = float(np.nanmax(rgf / np.maximum(bgf, 1e-6))
                    / max(np.nanmin(rgf / np.maximum(bgf, 1e-6)), 1e-6))
        stats["radial_colour"] = {
            "r_from": float(centres[idx[0]]), "r_to": float(centres[idx[-1]]),
            "rb_swing": round(_sw, 2),
            "rg_inner": round(float(rgf[idx[0]]), 3),
            "rg_outer": round(float(rgf[idx[-1]]), 3),
            "bg_inner": round(float(bgf[idx[0]]), 3),
            "bg_outer": round(float(bgf[idx[-1]]), 3)}
        progress.log(
            f"radial colour measured {centres[idx[0]]:.2f}-{centres[idx[-1]]:.2f} R: "
            f"R/G {rgf[idx[0]]:.3f} -> {rgf[idx[-1]]:.3f}, "
            f"B/G {bgf[idx[0]]:.3f} -> {bgf[idx[-1]]:.3f}, so red over blue "
            f"swings {_sw:.2f}x across the field. One set of white-balance "
            f"gains cannot flatten that; the Neutralise radial slider can. "
            f"This is a CHOICE, not a correction: a corona low in the sky "
            f"really is reddened, and flattening it says you would rather see "
            f"the K-corona's own white.", None)
    except Exception as _e:
        progress.log(f"radial colour not measured ({_e})", None)


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


def resolve_calib_dir(folder, which, arg=None):
    """Where the bias or dark frames are. Same contract as resolve_flat_dir.

    `which` is "bias" or "dark". The conventional subfolder names are `bias/`
    and `darks/` beside the light frames, alongside `flats/`; list_raws() reads
    only the files directly in a folder, so none of the three can ever be
    mistaken for a light frame.
    """
    from . import dark as _dark
    s = (arg or "").strip()
    if s.lower() in ("off", "none", "no", "-"):
        return None
    if s:
        return os.path.abspath(os.path.expanduser(s))
    return (_dark.find_bias_dir(folder) if which == "bias"
            else _dark.find_dark_dir(folder))


def run(folder, progress: Progress, crop_pc=1600, denoise="fine",
        earthshine=False, despeckle=True, frames="all", export_tiers=False,
        tier_linear=False, flat_dir=None, feather="plain",
        wb_source="camera", demosaic_method="mhc",
        photometry="linfit", bias_dir=None, dark_dir=None,
        photo_solve="chain", tier_mode="exposure",
        fnrgf_preset="ours", stack_combine="mean",
        align_filter="isotropic", align_corr="semi"):
    from . import flat as _flat
    from . import dark as _dark
    wd = workdir(folder)
    paths = list_raws(folder)
    if len(paths) < 3:
        raise RuntimeError(f"only {len(paths)} raw files found in {folder}")
    progress.log(f"{len(paths)} raw files found", 0.01)
    # Say out loud which ECLIPSEFORGE_* switches this BUILD actually honours.
    # During a bisect a switch was set on a shell that then launched an older
    # build, which ignored it silently; the run looked like a test and was a
    # repeat of the default, bit for bit. A build that can read the variable
    # says so here, so its absence is the signal that the app was not reinstalled.
    _sw = {k: v for k, v in os.environ.items()
           if k.startswith("ECLIPSEFORGE_") and v not in ("", "0")}
    if _sw:
        progress.log("switches honoured by this build (%s): %s"
                     % (__version__, ", ".join("%s=%s" % kv
                                               for kv in sorted(_sw.items()))),
                     None)
    _flat_dir = resolve_flat_dir(folder, flat_dir)
    _bias_dir = resolve_calib_dir(folder, "bias", bias_dir)
    _dark_dir = resolve_calib_dir(folder, "dark", dark_dir)
    stats = {"version": __version__, "folder": folder, "n_files": len(paths),
             "options": {"denoise": denoise, "earthshine": bool(earthshine),
                         "despeckle": bool(despeckle), "frames": frames,
                         "export_tiers": bool(export_tiers),
                         "tier_linear": bool(tier_linear),
                         "feather": str(feather),
                         "wb_source": str(wb_source),
                         "demosaic": str(demosaic_method),
                         "photo_solve": str(photo_solve),
                         "fnrgf_preset": str(fnrgf_preset),
                         "stack_combine": str(stack_combine),
                         "align_filter": str(align_filter),
                         "align_corr": str(align_corr),
                         "tier_mode": str(tier_mode),
                         "flat_dir": _flat_dir,
                         "bias_dir": _bias_dir, "dark_dir": _dark_dir},
             "camera_info": read_camera_info(paths[0])}

    # --- metadata & tiers ---
    meta = {}
    for p in paths:
        sec, iso, ts = read_exif(p)
        meta[p] = {"sec": sec, "iso": iso, "ts": ts}
    tiers = {}
    # ONE TIER PER FRAME, when asked for.
    #
    # Grouping by EXIF shutter speed is the one place this pipeline TRUSTS the
    # exposure metadata rather than measuring past it. Everything downstream
    # already measures: `cal[s]` is fitted from the data, and a bench over three
    # header ladders showed the merged levels move by one constant and 0.6% of
    # shape when the times are wrong. But two frames only land in the same tier
    # -- and get averaged as though they were the same measurement -- because
    # their headers say so.
    #
    # On a third tester's set the exposure times were typed in by hand from the
    # video, so that is exactly the assumption that does not hold. Druckmüller,
    # Rušin & Minarovjech (2006) do not group at all: every frame is calibrated
    # and enters one weighted average. This is that, reached the cheap way.
    #
    # THE KEY STAYS THE EXPOSURE TIME, nudged by a relative 1e-9 per frame so
    # the dict keeps them apart. That is far below any measurement precision --
    # `s * cal[s]` is unchanged to nine figures, the ordering is unchanged, and
    # `_exp_name(s)` prints the same string -- so nothing downstream has to
    # learn a new kind of key. The alternative, an opaque key plus a parallel
    # seconds lookup, would touch every `s * cal[s]` in this file.
    if str(tier_mode or "exposure").strip().lower() == "frame":
        for k, p in enumerate(sorted(paths, key=lambda q: (meta[q]["sec"], q))):
            tiers[float(meta[p]["sec"]) * (1.0 + k * 1e-9)] = [p]
        progress.log(f"one tier per frame: {len(tiers)} frames, grouped by "
                     f"nothing. The exposure times are still read and still "
                     f"scale the merge, but no two frames are averaged together "
                     f"on the strength of a header agreeing with another "
                     f"header.", None)
    else:
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
            # EXIF writes "2026:08:12 13:47:53"; FITS DATE-OBS is ISO-8601 with
            # a T, "2026-08-12T13:47:53.0" (0.22.78). Splitting on a space alone
            # failed on every FITS bracket, so _epoch returned None for all of
            # them and the lunar track fell back to frame index -- which is a
            # different quantity with the same units and no error anywhere.
            d, tm = str(ts).replace("T", " ").split(" ", 1)
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
    # FITS is stored bottom-up and every other format is not, so say which way
    # up these were read. The first FITS run came out mirrored and the log gave
    # no hint why -- there was nothing about orientation in it at all.
    try:
        from .fits import is_fits as _isf, row_order as _row
        if paths and _isf(paths[0]):
            _ro = _row(paths[0])
            stats["fits_row_order"] = _ro
            progress.log(f"FITS row order: {_ro} — flipped to top-down"
                         if not _ro.startswith("TOP") else
                         f"FITS row order: {_ro} — used as stored", None)
    except Exception:
        pass

    # --- per-tier: decode, quality, intra-align, stack (half-res + full-res bayer) ---
    _edge = {}          # per tier: invalid border (top, bottom, left, right), full-res px
    # Exposure exponent in the merge weight. 1.0 is the historical behaviour and
    # stays unless _pick_weight_alpha measures something better on this data.
    _walpha = 1.0
    sat_half = {}
    stacks_half = {}
    stacks_bayer = {}
    quality = {}
    color_info = None
    n_done = 0
    hot = None
    # Per-tier combine. 'clip' keeps every aligned frame of a tier in memory at
    # once, which is why it is a choice: the peak working set goes up by roughly
    # (frames - 1) x one tier's worth of arrays.
    _clip_stack = str(stack_combine or "mean").strip().lower() == "clip"
    if _clip_stack:
        progress.log(f"per-tier combine: kappa-sigma clipped mean "
                     f"(kappa {STACK_KAPPA:g}) instead of a plain mean", None)
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
    # Master bias and per-second dark rate, built here for the same reason the
    # flat is: the build holds frame-sized accumulators and must not overlap
    # with a tier's worth of decoded lights.
    #
    # ORDER OF OPERATIONS. These are ADDITIVE and the flat is MULTIPLICATIVE,
    # so they come off first: (light - bias - rate*t) / F, not light/F - bias.
    # Applied after the clipping test, which has to see the sensor's own
    # response, and before the flat divides -- the two lines below the test.
    bias_master = dark_rate = None
    #: defects measured on the master dark, and per-Bayer-offset read noise
    #: measured on the master bias. Both stay OUT of `stats`: the first is a
    #: 44-megapixel bool array and report.json is JSON.
    _calib_defects = None
    _calib_read_noise = None
    _warm_logged = False        # the census below sits in the per-tier loop
    _nw = 0
    stats["calib"] = {"bias_dir": _bias_dir, "dark_dir": _dark_dir}
    if _bias_dir or _dark_dir:
        if os.environ.get("ECLIPSEFORGE_NO_DARK") not in (None, "", "0"):
            progress.log("bias/dark correction disabled by "
                         "ECLIPSEFORGE_NO_DARK", None)
            stats["calib"]["error"] = "disabled by switch"
        else:
            try:
                bias_master, dark_rate, _ci = _dark.load_or_build(
                    folder, _bias_dir, _dark_dir, None, progress, wd,
                    light_iso=(isos[0] if len(isos) == 1 else None),
                    max_light_seconds=(max(secs) if secs else None))
                _calib_defects = _ci.pop("defect_map", None)
                _rn = ((_ci.get("bias") or {}).get("read_noise_cfa")
                       if isinstance(_ci.get("bias"), dict) else None)
                if _rn and len(_rn) == 4 and all(v > 0 for v in _rn):
                    _calib_read_noise = [float(v) for v in _rn]
                stats["calib"] = dict(_ci, bias_dir=_bias_dir,
                                      dark_dir=_dark_dir)
            except Exception as e:
                progress.log(f"bias/dark correction skipped — the masters "
                             f"could not be built ({e})", None)
                stats["calib"] = {"bias_dir": _bias_dir, "dark_dir": _dark_dir,
                                  "error": str(e)}
                bias_master = dark_rate = None
    for s in secs:
        files = tiers[s]
        lums, sats, bayers = {}, {}, {}
        raw_bayers = []
        for p in files:
            rf = open_frame(p)
            if color_info is None:
                from .raw import pick_wb
                _wb, _wbname = pick_wb(rf, wb_source)
                color_info = {"wb": _wb.tolist(),
                              "cam2rgb": rf.cam2rgb.tolist(),
                              "sat_level": rf.sat_level,
                              "shape": list(rf.shape)}
                stats["wb_source"] = str(wb_source)
                stats["wb_used"] = _wbname
                stats["wb_gains"] = [float(x) for x in _wb]
                _dwb = getattr(rf, "daylight_wb", np.ones(3, np.float32))
                progress.log(
                    f"white balance: {_wbname} — R {_wb[0]:.4f} B {_wb[2]:.4f} "
                    f"(G = 1). Daylight, for comparison, is R {_dwb[0]:.4f} "
                    f"B {_dwb[2]:.4f}. This is a SETTING; it scales the raw "
                    f"channels before the camera matrix, so it moves the "
                    f"colour of everything in the frame together.", None)
                progress.log(f"demosaic: "
                             f"{'VNG (as PixInsight)' if demosaic_method == 'vng' else 'Malvar-He-Cutler'}",
                             None)
                # FITS and TIFF have no white balance to give (see
                # neutralise_corona_colour). Recorded here because this is the
                # only place that still knows which reader produced the frame.
                stats["no_camera_wb"] = bool(
                    np.allclose(rf.daylight_wb, 1.0, atol=1e-6)
                    and np.allclose(rf.cam2rgb, np.eye(3), atol=1e-6))
                # WHERE THE PER-CHANNEL BLACK LEVEL COMES FROM.
                #
                # A camera raw carries `black_level_per_channel` and `raw.py`
                # subtracts it per Bayer position at load, with a comment on
                # why: a leftover per-channel offset gets multiplied by the
                # white-balance gain and lands on the outer corona, where the
                # signal is weakest. Measured residue after that, from the
                # per-channel check on three real raw brackets: 0.14, 0.26 and
                # 2.14 ADU. So on raw the metadata has already answered this
                # and refitting it from the data would replace a measurement
                # with an inference.
                #
                # FITS has no such header. `fits.py` re-mosaics a 3-plane cube
                # into a CFA, so each plane's own floor survives as a
                # per-Bayer-position offset and nothing downstream removes it.
                # On a third tester's set the three floors are 102.7, 116.2 and
                # 222.8 ADU against one shared scalar of 111.8 -- 107 ADU left
                # in blue, which is nothing against the inner corona and is
                # the entire signal in the outer field. That is the colour
                # cast that changes with radius.
                stats["fit_channel_floors"] = bool(stats["no_camera_wb"])
                # DID THE FILE REPORT A BLACK LEVEL, or did we default to zero?
                # LibRaw always reports one, so on a camera raw the metadata has
                # already answered it. A FITS usually does not -- and since
                # 0.23.3 fits.py no longer accepts Siril's OFFSET card as one,
                # because there OFFSET is the sensor GAIN OFFSET SETTING and is
                # normally 0. The difference decides how far the shared pedestal
                # fit below may travel: a few ADU of residue after a known black
                # level, or a whole unknown one.
                stats["black_level_reported"] = bool(
                    getattr(rf, "black_level_reported", True))
                # Say out loud where the saturation level came from. Getting
                # this wrong is silent and catastrophic -- see RawFile: a body
                # that reports the 14-bit container ceiling rather than its
                # linearity limit put the threshold ABOVE the real clip and fed
                # the flat tops of blown highlights into the merge at full
                # weight. Now it is in the log of every run.
                _lm = getattr(rf, "linear_max", None)
                _ss = getattr(rf, "sat_source", "white_level")
                stats["sat_level"] = float(rf.sat_level)
                stats["sat_from_linear_max"] = _lm is not None
                stats["sat_source"] = _ss
                # Three sources, named apart. 0.22.36 printed "the reported
                # white level" on a run whose number had in fact come from the
                # observed maximum -- literally true, since that is where the
                # chain started, and wrong about the thing the line exists to
                # tell you.
                _sw = {
                    "linear_max":
                        f"the sensor's linearity limit ({_lm:.0f})"
                        if _lm is not None else "the sensor's linearity limit",
                    "observed":
                        "THIS FRAME'S OWN MAXIMUM — the body reported no usable "
                        "linearity limit and its white level sat above the data, "
                        "so the highest photosite in the frame was taken as the "
                        "ceiling. Right whenever something in the frame is "
                        "genuinely blown, which on a totality frame it is; a "
                        "bracket with nothing saturated in its first frame gets "
                        "no such correction",
                    "container":
                        "the file's container range (not a raw sensor file)",
                    # ... and the two a FITS frame can report (0.22.78). Both
                    # fell through to "the reported white level", which is not
                    # what happened and hid the one that is a guess.
                    "SATURATE keyword":
                        "the SATURATE keyword in the FITS header",
                }.get(_ss, ("the FITS header's " + _ss) if "plateau" in _ss
                      else (_ss if _ss.startswith("bit depth")
                            else "the reported white level "
                                 "(no linearity limit given)"))
                progress.log(
                    f"saturation level {rf.sat_level:.0f} ADU above black, from "
                    + _sw, None)
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
        # Same guard, same reason: a bias or dark from another body or another
        # crop mode is not one for these frames, and a sentence here beats a
        # broadcast error at the subtraction.
        for _nm, _arr in (("bias", bias_master), ("dark", dark_rate)):
            if _arr is not None and _arr.shape != tuple(_lshape):
                progress.log(f"{_nm} correction DISABLED — the master is "
                             f"{_arr.shape[1]}x{_arr.shape[0]} px and the "
                             f"{_exp_name(s)} frames are {_lshape[1]}x"
                             f"{_lshape[0]}; calibration frames have to come "
                             f"from the same camera in the same crop mode",
                             None)
                stats["calib"][f"{_nm}_error"] = f"{_nm}/light size mismatch"
                if _nm == "bias":
                    bias_master = None
                else:
                    dark_rate = None
        stats["calib"]["bias_applied"] = bias_master is not None
        stats["calib"]["dark_applied"] = dark_rate is not None
        # sensor defects: map them once, on the shortest tier (darkest sky, so
        # real sky objects cannot be mistaken for hot pixels), then repair every
        # frame of every tier with it
        # THE DARK'S WARM PIXELS ARE NOT THE LIGHT FRAME'S DEFECTS, and this
        # census exists to say so rather than to replace anything.
        #
        # It was briefly the other way round -- the dark's map preferred, on
        # the reasoning that a dark sees a defect with no sky in the way. That
        # reasoning is wrong twice over, and the arithmetic is not close:
        #
        #  * At 1/500 s the worst photosite on the reference sensor, 122.9
        #    ADU/s, contributes 0.25 ADU. Dark current cannot make a
        #    detectable outlier in a 2 ms frame. So the pixels the light tier
        #    flags are NOT warm pixels -- they are offset or stuck-response
        #    defects, which is precisely the population that subtracting a
        #    dark does not fix.
        #  * A warm pixel, conversely, is already corrected: its own rate is
        #    subtracted per photosite. Replacing it with a neighbour median
        #    would discard real data to repair something already repaired.
        #
        # And there is no defect population in a real dark to threshold for.
        # The rate distribution is a continuum with a heavy tail -- p50 0.234,
        # p90 0.942, p99 1.862, p99.9 9.77, max 122.9 ADU/s -- so a cut at six
        # times the read noise takes 0.31% of the sensor, 138907 photosites,
        # and trips the implausibility guard. It was right to trip.
        # ONCE, not once per tier. This sits inside the per-tier loop, and
        # without the guard it printed fourteen identical lines and took
        # fourteen slots in the slowest-steps table -- where a log call is
        # credited with whatever ran before it, so it read as if counting warm
        # pixels cost two minutes. The flag also skips the .sum(), which is a
        # full pass over 44 megapixels, not just the print.
        if _calib_defects is not None and not _warm_logged:
            _warm_logged = True
            _nw = int(_calib_defects.sum())
            stats["warm_pixels"] = _nw
            progress.log(
                f"dark: {_nw} photosites ({100.0 * _nw / _calib_defects.size:.3f}%) "
                f"carry dark current well above their neighbours. NOT repaired "
                f"and not the same thing as a defect: their own rate is "
                f"subtracted per photosite, which is the correct fix. The "
                f"defect map below is for what that cannot fix.", None)
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
                hot = hot_pixel_map([b for _, b, _ in raw_bayers[:4]],
                                    read_noise=_calib_read_noise)
                nhot = int(hot.sum()) if hot is not None else 0
                progress.log(f"sensor defect map: {nhot} hot/dead photosites "
                             f"({100.0 * nhot / max(hot.size, 1):.4f}%)"
                             + (f", threshold from the measured read noise "
                                f"({min(_calib_read_noise):.2f}-"
                                f"{max(_calib_read_noise):.2f} ADU) rather than "
                                f"one fitted from this tier"
                                if _calib_read_noise else ""), None)
                # BOTH NUMBERS, because the bench cannot say which is right on
                # a given sensor. A synthetic tier makes the single-frame fit
                # read HIGH (pinning the intercept found 97 planted defects
                # against 80); this bracket made it read LOW enough to halve
                # the count. Same failure, opposite sign. Printing the two
                # side by side is the only thing that tells us which happened
                # here, and by how much.
                if _calib_read_noise:
                    _ff = getattr(hot_pixel_map.__globals__["_outlier_flags"],
                                  "last_fit_rn", None)
                    if _ff and all(np.isfinite(v) for v in _ff):
                        progress.log(
                            "  read noise: measured from the bias "
                            + "/".join(f"{v:.2f}" for v in _calib_read_noise)
                            + " ADU; the single-frame fit on this tier would "
                              "have used " + "/".join(f"{v:.2f}" for v in _ff)
                            + " ADU. The threshold is 6x these, so the ratio is "
                              "how far the old count was off.", None)
                        stats["read_noise_fitted"] = [round(float(v), 3)
                                                      for v in _ff]
                        stats["read_noise_measured"] = [round(float(v), 3)
                                                        for v in _calib_read_noise]
                stats["hot_pixel_source"] = (
                    "shortest light tier, read noise from the bias"
                    if _calib_read_noise else "shortest light tier")
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
            # Additive corrections, before the multiplicative one. `bias` is
            # the residue left after the black level the raw file reported
            # (raw.py already took that off per Bayer position); `rate` is dark
            # current in ADU per second, so it is scaled by THIS tier's
            # exposure. On the short tiers that term is a small fraction of an
            # ADU and does nothing; on the long ones it is the whole point.
            if bias_master is not None:
                bay -= bias_master
            if dark_rate is not None:
                bay -= dark_rate * np.float32(s)
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
        # KEPT ALIGNED COPIES, only when the combine needs them all at once.
        _clip_h = [acc_h.copy()] if _clip_stack else None
        _clip_b = [acc_b.copy()] if _clip_stack else None
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
            _sh = ndimage.shift(lums[p], (dy, dx), order=3, mode="nearest")
            _sb = shift_bayer_even(bayers[p], dy * 2, dx * 2)
            acc_h += _sh
            acc_b += _sb
            if _clip_stack:
                _clip_h.append(_sh)
                _clip_b.append(_sb)
            del _sh, _sb
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
        # PER-PIXEL REJECTION, OR NONE AT ALL.
        #
        # Up to 0.23.1 this was `acc / len(use)` and nothing else: every frame
        # that passed the SHARPNESS cut was added and the sum divided. Frame
        # level rejection, no pixel-level rejection whatsoever, so a cosmic ray
        # or an uncaught hot pixel in one frame of six survives at 1/6 amplitude
        # in that tier alone -- and then gets amplified by MGN and read as
        # structure, or perturbs one FNRGF ring and prints as a circle.
        #
        # Druckmullerova's thesis sec. 4.1.2 p.38 lists four combines and says
        # exactly this about ours: the arithmetic mean is "the best in noise
        # reduction and is computationally the fastest, but is strongly affected
        # by outliers, i.e. impulse noise present especially in longer exposures
        # caused by hot pixels or when high-energy cosmic particles hit the
        # chip."
        #
        # WHY NOT A MEDIAN, which is the obvious answer and the wrong one here.
        # From the same page: the median "does nothing to smooth out the
        # quantization effects" when combining integer values. Our frames ARE
        # integers out of the raw and the short tiers sit a few ADU above the
        # floor, and contouring is the artifact a tester already reports when
        # several frames are stacked in one tier. A median would make the one
        # thing being complained about worse. Kappa-sigma ends in a MEAN of the
        # survivors, so it keeps the quantisation smoothing and drops the
        # outlier.
        #
        # Kappa is deliberately loose. With two or three frames the sample sigma
        # is a poor estimate and a tight kappa starts rejecting signal; the
        # point here is to catch a particle strike, not to trim a distribution.
        if _clip_stack and len(use) >= 3:
            _n0 = 0
            for _frames, _tag in ((_clip_h, "h"), (_clip_b, "b")):
                _out, _nrej = _clipped_mean(_frames, STACK_KAPPA)
                _n0 += _nrej
                if _tag == "h":
                    stacks_half[s] = _out
                else:
                    stacks_bayer[s] = _out
                del _out
            # acc_h is HALF resolution and acc_b is the full-res mosaic, so
            # the two are not the same size -- counting both as acc_h.size
            # under-counted the denominator by 2.5x and printed a rejection
            # rate 2.5x too high. On the reference bracket that read 0.67%
            # where the cut is calibrated to reject 0.27%, which is the number
            # that made this look wrong when it was right.
            _tot = len(use) * (acc_h.size + acc_b.size)
            quality[s]["clipped_frac"] = float(_n0) / max(_tot, 1)
            if _n0:
                progress.log(f"{_exp_name(s)}: clipped mean rejected "
                             f"{100.0 * _n0 / max(_tot, 1):.5f}% of frame-pixels "
                             f"({_n0} of {_tot})", None)
        else:
            stacks_half[s] = acc_h / len(use)
            stacks_bayer[s] = acc_b / len(use)
        if _clip_stack:
            _clip_h = _clip_b = None
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
        # kept so the alignment-residual check further down has a scale: a shift
        # error only means anything relative to the size of the subject
        stats["R_seed_half"] = float(_Rseed)
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

    # ---- THE THREE ALIGNMENT SELECTORS, and what each measured.
    #
    # All three come from Druckmullerova's registration work and none of them
    # existed here before 0.23.2. Benched on tiers built from the reference set's own merged
    # corona and his own master flat, so the lunar edge, the per-tier saturation
    # edge and the sensor's dust and PRNU are all present, with the shift
    # injected by an exact FFT phase ramp. Two numbers: rms error on a known
    # shift, and the ZERO-PULL -- the shift reported for a pair with no real
    # shift, which is the estimator being dragged by whatever both frames share.
    #
    #     pre-filter   window   correlation     rms err   zero-pull
    #     isotropic    same     cross            1.524      0.200    <- 0.23.1
    #     isotropic    diff     cross            1.523      0.208
    #     tangential   same     cross            1.338      0.196
    #     tangential   diff     cross            1.341      0.190
    #     isotropic    same     semi             1.187      0.041    <- now
    #     tangential   diff     semi             1.161      0.054
    #     isotropic    same     phase (skimage)  1.307      0.193
    #
    # SEMI-PHASE CORRELATION IS THE WIN: 22% off the error and 80% off the
    # zero-pull, which is the predicted mechanism -- an additive p is a noise
    # floor on the whitening, so the fixed pattern that sits in the same place in
    # both frames stops dragging the peak to (0,0). Of that 22%, a control puts
    # 7% on the sub-pixel method rather than the whitening: our parabolic peak
    # fit scores 1.414 against skimage's 20x upsampling at 1.524 on the same
    # correlation surface. The remaining 15% is p and q.
    #
    # THE TANGENTIAL PRE-FILTER IS REAL BUT SMALL ONCE SEMI IS ON: 12% against
    # plain cross-correlation, 2% against semi, which is inside the spread of 33
    # samples. It stays a selector rather than a default -- it is the right tool
    # when the lunar edge or a moving saturation edge dominates, which is a
    # property of the dataset, not of the method.
    #
    # DIFFERENT WINDOWS DO NOTHING HERE; see ASYM_WINDOW.
    _tangential = str(align_filter or "isotropic").strip().lower() == "tangential"
    _corr = str(align_corr or "semi").strip().lower()

    def _arc_blur(x, sigma, org):
        """T_sigma's blur: a Gaussian along a Sun-centred ARC, not a disc.

        Thesis eq. 4.14 / master's thesis eq. 4.1. Subtracting this leaves only
        structure that VARIES ALONG THE ARC. Anything constant around a circle
        -- the lunar limb, the saturation boundary whose radius moves with
        exposure time, a vignette, the radial gradient itself -- cancels by
        construction. That is the point: those are the two features she names as
        hijacking a corona registration, and an isotropic high-pass keeps both.

        Implemented by resampling to polar about the disc centre, blurring along
        the azimuth axis with a kernel whose width in SAMPLES varies as 1/r so
        that its ARC LENGTH is constant (the mistake an earlier polar attempt
        here made was a constant kernel in grid units, i.e. an arc growing with
        radius), then resampling back. Only the BLUR makes the round trip, so
        interpolation error cannot contaminate the residual.
        """
        h, w = x.shape
        _cy = cym - org[0]
        _cx = cxm - org[1]
        rmax = float(np.hypot(max(_cy, h - _cy), max(_cx, w - _cx))) + 2
        nr = int(max(rmax, 8))
        nphi = 1024
        rr = np.linspace(0.0, rmax, nr, dtype=np.float32)
        ph = np.linspace(0, 2 * np.pi, nphi, endpoint=False, dtype=np.float32)
        ys = _cy + rr[:, None] * np.sin(ph)[None, :]
        xs = _cx + rr[:, None] * np.cos(ph)[None, :]
        P = ndimage.map_coordinates(x, [np.clip(ys, 0, h - 1).ravel(),
                                        np.clip(xs, 0, w - 1).ravel()],
                                    order=1, mode="nearest").reshape(nr, nphi)
        # CONSTANT ARC LENGTH: sigma in samples is sigma_px / (r * dphi), so it
        # falls as 1/r. Applied in BANDS of rings whose sigma agrees to 12%,
        # which turns ~850 single-ring calls into ~30 block calls -- 0.29 s for
        # a 1200 px window instead of tens of seconds, for the same result to
        # far under a tenth of a pixel.
        dphi = 2 * np.pi / nphi
        B = np.empty_like(P)
        s_all = sigma / np.maximum(rr * dphi, 1e-3)
        _i = 0
        while _i < nr:
            _s = float(s_all[_i])
            if _s < 0.4:                       # narrower than a sample
                _j = max(int(np.searchsorted(-s_all, -0.4)), _i + 1)
                B[_i:_j] = P[_i:_j]
            elif _s > nphi / 6.0:              # wider than the ring
                _j = max(int(np.searchsorted(-s_all, -(nphi / 6.0))), _i + 1)
                B[_i:_j] = P[_i:_j].mean(axis=1, keepdims=True)
            else:
                _j = max(int(np.searchsorted(-s_all, -(_s * 0.88))), _i + 1)
                B[_i:_j] = ndimage.gaussian_filter1d(P[_i:_j], _s, axis=1,
                                                     mode="wrap")
            _i = _j
        del P
        yy = np.arange(h, dtype=np.float32)[:, None] - _cy
        xx = np.arange(w, dtype=np.float32)[None, :] - _cx
        rq = np.hypot(yy, xx) / max(rmax / (nr - 1), 1e-6)
        tq = (np.arctan2(yy, xx) % (2 * np.pi)) / dphi
        return ndimage.map_coordinates(
            B, [np.clip(rq, 0, nr - 1).ravel(), tq.ravel()],
            order=1, mode="nearest").reshape(h, w).astype(np.float32)

    # TWO DIFFERENT WINDOWS, WHICH LOOKS WRONG AND IS NOT.
    #
    # Master's thesis sec. 4.2.2 p.57: "It is not necessary to apply the same
    # window function on both images to be registered ... If the window
    # functions are different (i.e. in the best case one is rectangular and the
    # other is circular), it assures that they do not bring about similar
    # structures in the images which might lead to incorrect registration."
    #
    # An identical window is a feature BOTH frames carry, in the SAME place, so
    # it correlates with itself and pulls the peak toward zero shift -- the same
    # failure mode as the sensor's fixed pattern, from a source we introduced
    # ourselves. Separable Hann on one, circular Hann on the other.
    _hn = np.hanning(S).astype(np.float32)
    _win_rect = (_hn[:, None] * _hn[None, :])
    _rg = (np.arange(S, dtype=np.float32) - (S - 1) / 2.0) / ((S - 1) / 2.0)
    _rad = np.hypot(_rg[:, None], _rg[None, :])
    _win_circ = np.where(_rad < 1.0,
                         0.5 * (1.0 + np.cos(np.pi * np.clip(_rad, 0, 1))),
                         0.0).astype(np.float32)
    _win = (_win_rect, _win_circ) if ASYM_WINDOW else (_win_rect, _win_rect)

    def prep_pair(s1, s2):
        p1y, p1x = _org[s1]
        p2y, p2x = _org[s2]
        a = stacks_half[s1][p1y:p1y + S, p1x:p1x + S] / s1
        b = stacks_half[s2][p2y:p2y + S, p2x:p2x + S] / s2
        sat2 = ndimage.binary_dilation(sat_half[s2][p2y:p2y + S, p2x:p2x + S], iterations=4)
        med = np.median(a); sig = 1.4826 * np.median(np.abs(a - med))
        good = (~sat2) & (ndimage.gaussian_filter(a, 3) > med + 5 * sig)
        wgt = ndimage.gaussian_filter(good.astype(np.float32), 5)
        # HOW HARD TO HIGH-PASS BEFORE CORRELATING ACROSS EXPOSURES.
        #
        # This was a flat 25 px, which is aggressive: at the reference set's
        # geometry that is 0.08 R, so almost everything the corona actually
        # looks like is removed and the correlation is left with fine texture
        # -- which is where a 1/2000 s tier has nothing but noise.
        #
        # MEASURED. Two aligned tiers are taken from the exported set (so their
        # true relative shift is zero), a known sub-pixel shift is injected into
        # one by an exact FFT phase ramp, and the estimator has to recover it.
        # Adjacent tier pairs only, which is what the network actually links.
        # RMS error in px, over 6 random shifts x 4 pairs:
        #
        #   sigma        0.05R  0.10R  0.20R  0.32R  0.50R  0.75R   none
        #   600mm +/-3    3.01   2.64   2.40   2.35   2.42   2.54   2.18
        #   600mm +/-12   2.99   2.61   2.36   2.31   2.38   2.48   2.15
        #   360mm +/-3    1.01   0.90   0.77   0.72   0.71   0.71   0.71
        #
        # The old 25 px sits at the left end of that. 0.5 R is at the flat
        # optimum on both sets and at both shift magnitudes, and is 18-22%
        # better. Dropping the high-pass entirely is very slightly better again
        # on one set, and is NOT taken: with no high-pass the correlation is
        # driven by the overall brightness distribution, which is exactly what
        # differs between exposures when there is thin cloud or a different
        # distribution of diffuse light -- the case Druckmullerova's thesis
        # names as the hard one. A large but finite high-pass keeps almost all
        # of the gain without that exposure.
        #
        # REVERTED IN 0.22.28. On the 600 mm reference run the network residual
        # went 1.17 -> 2.67 px (half-res), per-tier limb spread 8 -> 12 px,
        # track scatter 1/2 -> 2/4 px, and the step took 1m25s instead of 13s.
        # The table above is real but it was measured on the wrong thing:
        #
        #   * it used the EXPORTED tiers, which are the OUTPUT of alignment --
        #     already registered, resampled and mean-stacked. Their low
        #     frequencies match almost perfectly, which is exactly what a weak
        #     high-pass needs and exactly what the aligner never gets.
        #   * both crops came from the same origin. Real pairs are windowed on
        #     their own per-tier origins with up to ~50 px between them.
        #   * it had no signal-weight mask; prep_pair multiplies by `wgt`.
        #
        # So the harness measured sub-pixel refinement on an easy pair and
        # called it alignment. A synthetic test that cannot fail the way the
        # real thing fails is not evidence. Back to 25 px, which is what
        # produced the 1.17 px residual.
        _hp = 25.0
        out = []
        for _k, img in enumerate((a, b)):
            x = np.log1p(np.clip(img, 0, None) / max(med + 5 * sig, 1e-3))
            if _tangential:
                x = x - _arc_blur(x, _hp, (p1y, p1x) if _k == 0 else (p2y, p2x))
            else:
                x -= ndimage.gaussian_filter(x, _hp)
            out.append(x * wgt * _win[_k])
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
            if _corr == "semi":
                sh, err = _semi_phase_shift(p1, p2)
            elif _corr == "phase":
                sh, err, _ = phase_cross_correlation(p1, p2, upsample_factor=20)
            else:
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
                    # "N/(n-1) tiers linked by them" read as coverage and was
                    # not: it counts LINKS, not tiers, so the test set's 560mm run
                    # said "11/11 tiers linked" on a bracket where five tiers
                    # had zero anchors. Say what the number is.
                    f"prominence anchors: {len(anchors)} on the "
                    f"{_exp_name(p_ref)} tier, giving {len(prom_links)} of "
                    f"{n - 1} possible links", None)
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

    # A RESIDUAL THIS LARGE IS A FAILED ALIGNMENT, NOT A STATISTIC.
    #
    # a second tester's 560 mm set (0.20.1) reported
    #
    #     alignment network residual max 512.09px (half-res)
    #
    # -- 1024 px full-res, about twice the lunar radius -- and the run carried
    # on to "ready" with nothing to say it had failed. Downstream the report
    # then showed a 226 px per-tier limb spread against a 4 px radius spread,
    # a 32.8 px limb-fit rms, a 69 px limb ramp, and shifts of 108 and 118 px on
    # two tiers whose own log lines said they could not be aligned at all.
    #
    # The unlinked-tier fallback above cannot help here: those tiers WERE
    # linked. The links themselves contradict each other, which is what a large
    # residual means -- least squares returns the best compromise between
    # mutually impossible constraints and it looks like an answer.
    #
    # There is no repair to make automatically: which link is wrong is exactly
    # what is unknown. So this says so, unmistakably, and names the worst
    # offenders so the cause is findable. The threshold is a fraction of the
    # lunar radius, because a shift error only matters relative to the subject.
    _rseed = float(stats.get("R_seed_half") or 0.0)
    _tol = max(0.05 * _rseed, 6.0) if _rseed > 0 else 20.0
    if res > _tol and len(A) > 1:
        _bad = np.argsort(np.maximum(_ry, _rx))[::-1][:3]
        _who = []
        for _b in _bad:
            _row = A[int(_b)]
            _ix = [k for k in range(n) if abs(_row[k]) > 0.5]
            if len(_ix) >= 2:
                _who.append("%s<->%s (%.0fpx)"
                            % (_exp_name(secs[_ix[0]]), _exp_name(secs[_ix[-1]]),
                               max(_ry[int(_b)], _rx[int(_b)])))
        progress.log(
            "WARNING: THE CROSS-TIER ALIGNMENT FAILED. The link network is "
            "inconsistent by %.0f px (half-res, %.0f px full-res) against a "
            "tolerance of %.0f px. That is not a measurement error -- the links "
            "contradict each other, and least squares has returned the best "
            "compromise between impossible constraints. Worst: %s. Everything "
            "downstream is built on those shifts: the merged limb, its fitted "
            "radius, the disc mask and every radial filter. Treat this run's "
            "output as unusable and check whether the long tiers have any "
            "corona to correlate on."
            % (res, 2 * res, _tol, "; ".join(_who) or "unidentified"), None)
        stats["align_failed"] = round(float(res), 1)
        stats["align_tolerance"] = round(float(_tol), 1)

    # ONLY REGISTER THE RING SAMPLES IF THE SHIFTS ARE WORTH TRUSTING.
    #
    # The pedestal, the per-tier radial check and the LDIC fit read each tier at
    # its own position from 0.22.76 on -- see _ring_sample, which is the right
    # thing to do and measurably so. But it makes those three estimators DEPEND
    # on abs_shift, where before they ignored it, and abs_shift is not always a
    # measurement: when the link network contradicts itself, least squares still
    # returns a number and the warning above is the only sign of it.
    #
    # Found by measuring rather than by reasoning: on a synthetic bracket whose
    # corona is too smooth for the 25 px correlation to lock onto, the network
    # returned 9 px of residual against a true 1.5 px, and feeding those shifts
    # to the LDIC fit made its spurious k spread WORSE than ignoring them. Where
    # the shifts are noise, the old behaviour is the safer one -- every tier
    # sampled in the same place is at least a consistent error.
    #
    # So the registration follows the alignment's own verdict. A run that prints
    # the alignment-failed warning gets the pre-0.22.76 sampling, and says so.
    _reg_shift = dict(abs_shift)
    if os.environ.get("ECLIPSEFORGE_NO_REGSAMPLE") == "1":
        # THE WAY BACK. This is the one switch that reproduces the pre-0.22.76
        # merge exactly: every tier sampled in the same frame, as it was. It
        # exists because "the ring was not there in the last version" is a
        # claim that needs an A/B on the SAME data, not an argument.
        _reg_shift = {s_: (0.0, 0.0) for s_ in secs}
        progress.log("ring sampling held UNREGISTERED by "
                     "ECLIPSEFORGE_NO_REGSAMPLE=1 — the pedestal, the radial "
                     "check and the azimuthal fit see exactly what they saw "
                     "before 0.22.76", None)
    elif "align_failed" in stats:
        _reg_shift = {s: (0.0, 0.0) for s in secs}
        progress.log("the pedestal, the radial check and the azimuthal fit will "
                     "sample every tier in the SAME frame, because the shifts "
                     "they would otherwise use are the ones that just failed "
                     "their consistency check", None)

    # --- photometric calibration ---
    cal = {mid: 1.0}
    logf = np.zeros(n)
    ratios = []
    # kept so the per-channel check below can be fitted on exactly these pixels
    _links_good = [None] * (n - 1)
    #: (log ratio, px) per adjacent link, or None where it was rejected.
    #: `ratios` cannot carry that distinction -- a rejected link and a
    #: link that measured exactly 1.000 both land on 0.0 there.
    _lag1 = [None] * (n - 1)
    def _measure_link(i, j, keep_good=True):
        """log(ratio) between two tiers' per-second signal, or None.

        Pulled out of the lag-1 loop so the network solve can ask for
        i -> i+2 as well. The measurement is unchanged; what was
        `ratios.append(0.0); continue` is now `return None`, so a
        caller can tell a rejected link from one that measured 0.
        """
        a = stacks_half[secs[i]] / secs[i]
        b = stacks_half[secs[j]] / secs[j]
        m = ~ndimage.binary_dilation(sat_half[secs[j]], iterations=6)

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
        if keep_good:
            _links_good[i] = good
        if good.sum() < 1000:
            progress.log(f"photometric link {_exp_name(secs[i])}->"
                         f"{_exp_name(secs[j])}: only {int(good.sum())} px "
                         f"carry signal in both tiers; using exposure time", None)
            return None
        # REVERTED IN 0.16.2. Read this before changing it again.
        #
        # 0.16.0 replaced this median with sum(b[m])/sum(a[m]) -- no
        # data-dependent selection -- because selecting on a tier's own noise
        # biases that tier upward inside the selection, and the error is
        # one-sided so it compounds down a long chain. That reasoning is still
        # correct, and on synthetic tier pairs whose true ratio is 1.000 the
        # sum estimator returned 1.0015 over nine links against this one's
        # 0.785.
        #
        # It was still wrong, because the synthetic scene was not a real frame.
        # On the reference bracket (8100x5357, corona to 4 R, so most of the
        # picture is sky) the sums are dominated by sky area rather than by
        # corona, and at the short end the sky carries no signal to ratio. The
        # first three links flipped from 0.833/0.910/1.009 to
        # 1.307/1.177/1.091 -- an over-correction in the OPPOSITE direction --
        # and the shortest tier's factor fell from 1.273 to 0.571:
        #
        #                       0.14.5    0.16.1
        #     1/4000s            1.273     0.571
        #     1/2000s            1.061     0.746
        #     1/500s             0.965     0.878
        #     ...
        #     1.6s               1.158     1.326
        #
        #     disagreement rim   2 px      128 px
        #     limb variance      0.072     0.255
        #
        # 128 px of tier disagreement outside the limb is the bright rim, and
        # it is visible in MGN, NAFE and FNRGF. So: the old estimator is biased,
        # the new one is worse on real data, and neither is right.
        #
        # WHAT A REPLACEMENT HAS TO DO. Take the ratio over a region chosen
        # GEOMETRICALLY -- an annulus a little outside the limb, where both
        # tiers of a pair have real signal -- so that the selection depends on
        # neither tier's noise and the sky cannot dominate the sum. And it has
        # to be validated on the reference bracket AND on a wide-field FITS set
        # before it ships, not on synthetics alone. That was the actual mistake
        # here: a change to the core merge went out on simulated evidence.
        lr = float(np.log(np.median(b[good] / np.maximum(a[good], 1e-6))))
        if abs(lr) > np.log(2.5):
            # a tier cannot really be 2.5x off its own exposure time; this means
            # the comparison region is noise, not signal (very short tiers) --
            # trust the shutter speed instead of the measurement
            progress.log(f"photometric link {_exp_name(secs[i])}->"
                         f"{_exp_name(secs[j])} rejected (x{np.exp(lr):.2f} "
                         f"over {int(good.sum())} px) — using exposure time", None)
            return None
        return lr, int(good.sum())

    for i in range(n - 1):
        _r = _measure_link(i, i + 1)
        _lag1[i] = _r
        ratios.append(0.0 if _r is None else _r[0])
    # --- chain, or network -------------------------------------------------
    #
    # THE CHAIN is what this has always done: walk out from the middle tier
    # adding one measured link at a time. It uses each link exactly once, so a
    # link that is off by 3% moves every tier beyond it by 3%, and the errors
    # compound in one direction down the ladder. With only adjacent links
    # measured there is nothing else it could do -- there is no redundancy to
    # average over.
    #
    # THE NETWORK measures i -> i+2 as well and solves all of them together by
    # weighted least squares, minimising
    #
    #     sum_k  w_k * (logf[j_k] - logf[i_k] - lr_k)^2
    #
    # with logf[ref] pinned. Now each tier is constrained by up to four links
    # instead of two, and a single bad one is outvoted rather than propagated.
    # Where the lag-2 links are all rejected this reduces to the chain exactly,
    # which is the honest reason it is offered as a SETTING and not forced: on
    # a bracket with no redundancy the two are the same computation, and on one
    # where the solve misbehaves the chain is a known-good fallback.
    #
    # Same shape as the alignment network above, and for the same reason.
    _solve = str(photo_solve or "chain").strip().lower()
    if _solve not in ("chain", "network"):
        _solve = "chain"
    _net = None
    if _solve == "network" and n >= 3:
        _lag2 = []
        for i in range(n - 2):
            _r = _measure_link(i, i + 2, keep_good=False)
            if _r is not None:
                _lag2.append((i, i + 2, _r[0], _r[1]))
        _all = [(i, i + 1, v[0], v[1]) for i, v in enumerate(_lag1)
                if v is not None] + _lag2
        _sol, _res = _solve_photo_network(n, _all, ref)
        if _sol is None:
            progress.log("photometric network not solvable on this bracket "
                         "(a tier no link reaches, or too few links) — using "
                         "the chain", None)
        else:
            logf = _sol
            _net = {"links": len(_all), "lag2": len(_lag2),
                    "resid_pct": [round(100.0 * (np.exp(r) - 1.0), 3)
                                  for r in _res]}
        if _net is None:
            progress.log(f"photometric solve: network requested but this "
                         f"bracket has no usable redundancy "
                         f"({len(_lag2)} two-step link(s) measured) — using the "
                         f"chain, which is the same answer", None)
            _solve = "chain"
    if _net is None:
        for i in range(ref, n - 1):
            logf[i + 1] = logf[i] + ratios[i]
        for i in range(ref, 0, -1):
            logf[i - 1] = logf[i] - ratios[i - 1]
    cal = {s: float(np.exp(logf[i])) for i, s in enumerate(secs)}
    stats["photo_solve"] = _solve
    if _net is not None:
        stats["photo_network"] = _net
        _rp = np.abs(np.asarray(_net["resid_pct"], float))
        # THE RESIDUAL IS THE POINT. A chain reproduces every link exactly by
        # construction, so residuals near zero here mean the solve found no
        # disagreement to arbitrate and the two paths agree; residuals of a few
        # percent mean the links genuinely disagreed and the network is the one
        # averaging them instead of picking a path through them.
        progress.log(
            f"photometric solve: network over {_net['links']} links "
            f"({_net['lag2']} of them two-step). Residual per link, i.e. how "
            f"far each measurement sits from what the solve settled on: median "
            f"{np.median(_rp):.2f}%, worst {_rp.max():.2f}%. Near zero means "
            f"the links agreed and this is the chain's answer.", None)

    # --- the pedestal every tier shares ---
    # Ring medians over the OUTER field, where the corona is faint enough that a
    # few ADU of leftover black level is a large fraction of the signal and the
    # tiers can be seen to disagree. Inside ~1.5 R the corona swamps it and there
    # is nothing to measure. See _fit_pedestal for the evidence.
    #
    # cal[] above is fitted before this and is therefore itself slightly biased
    # by the pedestal -- but its ratios are taken over pixels selected as being
    # well above each tier's noise floor, where a 3 ADU offset on a 14-bit raw is
    # parts per thousand. Fitting the two jointly would be more correct and is
    # not worth the coupling.
    pedestal = 0.0
    ped_channel = None      # per-channel deltas about the shared one
    # BISECT SWITCH. Three merge changes landed in quick succession (0.22.15
    # pedestal, 0.22.16/.19 feather) and a user reported a texture regression on
    # a set that had been clean. Guessing which one is slower than letting the
    # person with the data turn each off for one run, so:
    #   ECLIPSEFORGE_NO_PEDESTAL=1   restores the pre-0.22.15 merge exactly
    #   ECLIPSEFORGE_FEATHER=plain   hides the rings, leaks at the limb
    #   ECLIPSEFORGE_FEATHER=masked  the 0.22.16 one (steps at the clip edge)
    #   ECLIPSEFORGE_FEATHER=taper   leak-free, rings visible
    # Since 0.22.28 the feather is CHOSEN BY MEASUREMENT per dataset rather
    # than defaulted; setting the variable overrides that. See _pick_feather.
    if os.environ.get("ECLIPSEFORGE_NO_PEDESTAL") == "1":
        progress.log("shared pedestal DISABLED by ECLIPSEFORGE_NO_PEDESTAL=1 — "
                     "this run reproduces the pre-0.22.15 merge", None)
        stats["pedestal_disabled"] = True
    try:
        if os.environ.get("ECLIPSEFORGE_NO_PEDESTAL") == "1":
            raise RuntimeError("disabled")
        _Hh, _Wh = stacks_half[secs[0]].shape
        _rmax = min(cym, cxm, _Hh - cym, _Wh - cxm) - 2.0
        _rs = np.arange(1.5 * _Rseed, min(3.5 * _Rseed, _rmax), 0.05 * _Rseed)
        if _rs.size >= 8:
            _th = np.linspace(0, 2 * np.pi, 720, endpoint=False)
            # the COMMON frame's rings; each tier is read in its own frame --
            # see _ring_sample for what sampling them unregistered cost here
            _yr = cym + _rs[:, None] * np.sin(_th)
            _xr = cxm + _rs[:, None] * np.cos(_th)
            _pr, _sc = [], []
            for s in secs:
                _v, _c = _ring_sample(stacks_half[s], sat_half[s],
                                      _yr, _xr, _reg_shift[s])
                # a ring counts only if essentially none of it is clipped: the
                # median of a partly clipped ring is biased low and would be
                # read as a negative pedestal
                _keep = (~_c).mean(axis=1) > 0.98
                _pr.append(np.where(_keep, np.median(_v, axis=1), np.nan))
                _sc.append(s * cal[s])
            _satadu = float(color_info["sat_level"])
            # HOW FAR THE FIT MAY TRAVEL depends on what the file told us.
            # _PEDESTAL_MAX is sized for the RESIDUE left after a KNOWN black
            # level, and the three real raw sets that show one land at 2.9, 5.5
            # and -2.0 ADU. When the file reported no black level the whole of
            # it is still in the data, and that ceiling makes it unrecoverable:
            # on a third tester's FITS the true level is 142 ADU against a 131 ADU
            # cap, so the fit pinned against its own limit and left ~23 ADU in
            # every channel. The prior is still zero and the jackknife shrinkage
            # is unchanged, so a set with no pedestal does not acquire one --
            # this only stops a real answer being clipped.
            _pmax = _PEDESTAL_MAX * _satadu
            if not stats.get("black_level_reported", True):
                _pmax = _PEDESTAL_UNKNOWN_MAX * _satadu
                progress.log(
                    f"the file reported no black level, so the shared pedestal "
                    f"may range to {_pmax:.0f} ADU rather than "
                    f"{_PEDESTAL_MAX * _satadu:.0f}: what is left in the data "
                    f"is the whole black level, not a residue of one", None)
            pedestal, _praw, _psig, _b0, _b1 = _fit_pedestal(
                np.asarray(_pr, float), np.asarray(_sc, float), _pmax)
            if np.isfinite(_b0) and np.isfinite(_b1):
                stats["pedestal"] = {"applied": float(pedestal),
                                     "fitted": float(_praw),
                                     "sigma": float(_psig),
                                     "scatter_before": float(_b0),
                                     "scatter_after": float(_b1)}
                # sigma comes from a jackknife on a grid and lands at exactly 0
                # whenever every leave-one-out subset picks the same grid point.
                # Dividing by it printed "1238016000000.0 sigma" in a real run.
                _sgtxt = (f"{abs(_praw) / _psig:.1f} sigma" if _psig > 1e-9
                          else "every leave-one-tier-out subset agrees")
                progress.log(
                    f"shared pedestal {pedestal:+.2f} ADU "
                    f"(fitted {_praw:+.2f}, {_sgtxt}): tier-to-tier "
                    f"disagreement in the outer field {100 * _b0:.1f}% -> "
                    f"{100 * _b1:.1f}%. A black level left a few ADU behind "
                    f"arrives divided by the exposure time, so it is nothing on "
                    f"a long tier and everything on a short one.", None)
    except Exception as _e:
        if str(_e) != "disabled":
            progress.log(f"shared pedestal not measurable ({_e}) — the outer "
                         f"field is merged as it stands", None)

    # --- Hill's per-channel linear-fit chain, when selected ---
    # Computed HERE, after cal[] and the shared pedestal, because it replaces
    # both and the log should be able to say by how much. See _hill_chain.
    _hillK = _hillQ = None
    # The marker file's preferred name says what it DOES; the old one is still
    # honoured so a folder set up before 0.22.75 keeps working.
    _hill_on = _photometry_requested(folder, photometry)
    if _hill_on:
        try:
            _hillK, _hillQ = _hill_chain(stacks_bayer, sat_half, secs,
                                         _links_good, cal,
                                         float(color_info["sat_level"]),
                                         pedestal, progress, stats)
        except Exception as _e:
            _hillK = _hillQ = None
            progress.log(f"per-channel linear-fit photometry failed ({_e}) — using the "
                         f"scalar chain", None)
    if _hillK is not None:
        # Everything downstream that works on LUMINANCE -- the feather trial,
        # the per-tier radial check, the LDIC fit, the tier TIFF export -- takes
        # its scale from cal[] and its offset from `pedestal`. Rather than
        # thread a per-channel affine through all of them, give them the
        # channel-mean equivalent of the same transform, which is accurate to
        # the colour spread the fit found (measured at ~1% on the reference
        # bracket). The MERGE, which makes the picture, uses the full
        # per-channel form below.
        stats["photometry_mode"] = "hill"
        _cal_scalar = dict(cal)
        for s in secs:
            _kb = float(np.mean(_hillK[s]))
            if np.isfinite(_kb) and _kb > 1e-9:
                cal[s] = 1.0 / _kb
        # The shared pedestal STAYS. It was already removed before the pairwise
        # fits above, so Q carries only the per-pair part and the two compose
        # rather than duplicate.
        _ped_scalar = pedestal
        progress.log(
            f"the scalar chain is kept for reference only: it read "
            f"{', '.join(f'{_cal_scalar[s]:.3f}' for s in secs)}; the Hill "
            f"chain's equivalent reads "
            f"{', '.join(f'{cal[s]:.3f}' for s in secs)}. The shared pedestal "
            f"of {_ped_scalar:+.2f} ADU is applied as well: it comes off "
            f"before the pairwise fits, so it removes the part common to the "
            f"whole bracket and the fitted offsets remove only what differs "
            f"between neighbours.", None)

    # --- does any tier's error VARY WITH RADIUS? ---
    #
    # The photometric factor is one scalar per tier, fitted where tiers overlap
    # -- mostly the mid field. An error that changes with radius therefore
    # survives calibration untouched, and lands hardest just outside the limb
    # where the filters normalise against a local mean and amplify it.
    #
    # This is measured rather than assumed because a rough version of it, run
    # outside the pipeline on four of the reference raws (one raw per tier, my own
    # centre, the stacked tier's shift applied to a single frame), showed two
    # adjacent tiers disagreeing with the other two by -22% and +64% inside
    # 1.25 R while agreeing to a few percent outside it. That rig was too crude
    # to trust for anything but "look here", so the pipeline now takes the same
    # measurement with the real alignment, the real stacks and the real
    # geometry, on every dataset, and writes it into the diagnostics bundle.
    #
    # It changes nothing. It is a number, so that the next decision about
    # per-tier radial correction is made on five datasets instead of one.
    try:
        _Hh, _Wh = stacks_half[secs[0]].shape
        _rmax = min(cym, cxm, _Hh - cym, _Wh - cxm) - 2.0
        _rs = np.arange(1.00 * _Rseed, min(3.0 * _Rseed, _rmax),
                        0.02 * _Rseed)
        if _rs.size >= 10:
            _th = np.linspace(0, 2 * np.pi, 720, endpoint=False)
            _yr = cym + _rs[:, None] * np.sin(_th)   # see _ring_sample
            _xr = cxm + _rs[:, None] * np.cos(_th)
            _prof = {}
            for s in secs:
                _v, _c = _ring_sample(stacks_half[s], sat_half[s],
                                      _yr, _xr, _reg_shift[s])
                _v = _v - pedestal
                _keep = (~_c).mean(axis=1) > 0.98
                _prof[s] = np.where(_keep, np.median(_v, axis=1),
                                    np.nan) / (s * cal[s])
            _A = np.asarray([_prof[s] for s in secs], float)
            with np.errstate(all="ignore"):
                _ref = np.nanmedian(_A, axis=0)
                _rel = _A / np.where(np.abs(_ref) < 1e-12, np.nan, _ref)
            stats["tier_radial"] = {
                "radius_R": [round(float(r / _Rseed), 3) for r in _rs],
                "rel": {("%g" % s): [None if not np.isfinite(v) else round(float(v), 4)
                                     for v in _rel[i]]
                        for i, s in enumerate(secs)}}
            # worst offender inside 1.3 R, where a scalar factor cannot see it
            _in = _rs < 1.3 * _Rseed
            if _in.sum() > 3:
                _dev = np.nanmax(np.abs(_rel[:, _in] - 1.0), axis=1)
                _k = int(np.nanargmax(_dev))
                _worst = float(_dev[_k])
                stats["tier_radial_worst"] = round(_worst, 3)
                if np.isfinite(_worst) and _worst > 0.15:
                    progress.log(
                        f"per-tier radial check: the {_exp_name(secs[_k])} tier "
                        f"departs from the other tiers by up to "
                        f"{100 * _worst:.0f}% inside 1.3 R while they agree "
                        f"outside it. The photometric factor is ONE number per "
                        f"tier, fitted where the tiers overlap, so an error "
                        f"shaped like this passes through it untouched and the "
                        f"detail filters amplify it just outside the limb. "
                        f"Recorded, not corrected — the profiles are in the "
                        f"diagnostics bundle.", None)
                else:
                    progress.log(
                        f"per-tier radial check: tiers agree to "
                        f"{100 * (_worst if np.isfinite(_worst) else 0):.0f}% "
                        f"inside 1.3 R — no radius-dependent per-tier error",
                        None)
    except Exception as _e:
        progress.log(f"per-tier radial check skipped ({_e})", None)

    # ------------------------------------------------------------------
    # LDIC's per-tier AZIMUTHAL affine transform. See _fit_azimuthal_affine.
    # Fitted here on the half-res stacks and applied in the merge loop below.
    #
    # EACH TIER IS SAMPLED IN ITS OWN FRAME (0.22.76). The comment here used to
    # say the half-res stacks are "already in the common frame". They are not:
    # `stacks_half[s]` is intra-tier aligned only, and `abs_shift[s]` is applied
    # in the merge loop, AFTER this fit. So every tier was sampled on rings
    # centred at the mid tier's centre, and a tier offset by d px had its
    # profile read at r + d*cos(phi - phi0) -- a cos(phi) modulation of the
    # sampled radiance. That is a DIPOLE in azimuth, which is precisely the
    # shape k(phi) is fitted to represent; the azimuthal-mean-1 normalisation
    # removes the monopole and leaves the dipole intact, and the merge then
    # multiplied the REGISTERED tier by it.
    #
    # tools/ldic_shift_bench.py builds tiers that differ only by exposure time
    # and a rigid shift, so the true k(phi) is 1 everywhere and anything the fit
    # reports is error. The 10-90 percentile span of k, in percent:
    #
    #     tier shift    as it was    sampled in the tier's own frame
    #     (half-res)   (unshifted)   integer      bilinear
    #        0 px          0.0         0.0          0.0
    #        2 px         20.9         4.8          1.7
    #        5 px         64.8         0.0          0.0
    #       11 px        130.9         4.4          1.7
    #       25 px        242.7         0.0          0.0
    #
    # For scale: the reference bracket's own fit reports k spreads of 3-27%, so
    # on a set with ten pixels of drift the reported variation was mostly this.
    # The residual 1.7% is bilinear interpolation across the limb gradient and
    # does not fall with a better sampler.
    #
    # ndimage.shift(img, +s) gives out[u] = in[u - s] (see the sign note above),
    # so the point that lands at (y, x) after registration is (y - ady, x - adx)
    # in the tier's own frame. Sampled bilinearly rather than by an integer cast
    # because the shifts are sub-pixel and rounding them re-introduces a
    # fraction of the same dipole -- the middle column above.
    _ldic, _ldic_unreg = {}, {}
    _ldic_rout = 2.5 * _Rseed        # half-res; the fit's outer radius, and the
                                     # radius past which the merge tapers it off
    try:
        if os.environ.get("ECLIPSEFORGE_NO_LDIC") == "1":
            raise RuntimeError("disabled by ECLIPSEFORGE_NO_LDIC=1")
        _Hh, _Wh = stacks_half[secs[0]].shape
        _rmax = min(cym, cxm, _Hh - cym, _Wh - cxm) - 2.0
        _rs = np.arange(1.00 * _Rseed, min(2.5 * _Rseed, _rmax), 0.02 * _Rseed)
        if _rs.size:
            _ldic_rout = float(_rs[-1])
        _NA = 12 * _LDIC_SEGMENTS
        if _rs.size >= 10:
            _th = np.linspace(0, 2 * np.pi, _NA, endpoint=False)
            _yr = cym + _rs[:, None] * np.sin(_th)   # the COMMON frame's rings
            _xr = cxm + _rs[:, None] * np.cos(_th)
            _sl = float(color_info["sat_level"])
            _order = sorted(secs, reverse=True)      # longest exposure first
            _vals, _wts = [], []
            for s in _order:
                _v, _sat = _ring_sample(stacks_half[s], sat_half[s],
                                        _yr, _xr, _reg_shift[s])
                _f = _v - pedestal
                _fr = _v / max(_sl, 1e-9)
                _u = np.clip((0.85 - _fr) / 0.10, 0, 1)
                _l = np.clip((_fr - 0.004) / 0.008, 0, 1)
                _w = (_u * _u * (3 - 2 * _u)) * (_l * _l * (3 - 2 * _l))
                _w = np.where(_sat, 0.0, _w)
                _vals.append(_f / (s * cal[s]))
                _wts.append(_w)
                del _v, _sat
            _ldic = _fit_azimuthal_affine(_vals, _wts, _order)
            del _vals, _wts
            # AND THE SAME FIT WITHOUT REGISTERING, so the log can say how much
            # of this correction was the tiers' own misalignment rather than
            # anything in the sky. Costs one more fit on a ring grid -- seconds
            # -- and it is the only way to see the size of the 0.22.76 change on
            # real data instead of on a bench.
            _ldic_unreg = {}
            if any(abs(v[0]) + abs(v[1]) > 1e-6 for v in _reg_shift.values()):
                try:
                    _v2, _w2 = [], []
                    for s_ in _order:
                        _vv, _ss = _ring_sample(stacks_half[s_], sat_half[s_],
                                                _yr, _xr, (0.0, 0.0))
                        _fr2 = _vv / max(_sl, 1e-9)
                        _u2 = np.clip((0.85 - _fr2) / 0.10, 0, 1)
                        _l2 = np.clip((_fr2 - 0.004) / 0.008, 0, 1)
                        _w2.append(np.where(
                            _ss, 0.0, (_u2 * _u2 * (3 - 2 * _u2))
                            * (_l2 * _l2 * (3 - 2 * _l2))))
                        _v2.append((_vv - pedestal) / (s_ * cal[s_]))
                        del _vv, _ss
                    _ldic_unreg = _fit_azimuthal_affine(_v2, _w2, _order)
                    del _v2, _w2
                except Exception:
                    _ldic_unreg = {}
        if _ldic:
            _sp = {("%g" % k): round(float(100 * (np.percentile(v[0], 90) -
                                                  np.percentile(v[0], 10))), 1)
                   for k, v in _ldic.items()}
            stats["ldic_k_spread_pct"] = _sp
            if _ldic_unreg:
                _spu = {("%g" % k): round(float(100 * (
                            np.percentile(v[0], 90) - np.percentile(v[0], 10))), 1)
                        for k, v in _ldic_unreg.items()}
                stats["ldic_k_spread_pct_unregistered"] = _spu
                _both = [(k, _spu[k], _sp[k]) for k in _sp if k in _spu]
                if _both:
                    _wu = max(_both, key=lambda t: t[1])
                    progress.log(
                        "azimuthal fit, measured both ways: sampling every tier "
                        "in the same frame (as before 0.22.76) makes the gain "
                        "look like it varies by up to %.0f%% around the limb; "
                        "sampling each tier where it actually sits gives %.0f%% "
                        "on that tier. The difference is the tiers' own "
                        "misalignment being read as sky. Per tier: %s"
                        % (_wu[1], _wu[2],
                           ", ".join("%s %.0f%%->%.0f%%"
                                     % (_exp_name(float(k)), u, r)
                                     for k, u, r in _both)), None)
            _wk = max(_sp.items(), key=lambda kv: kv[1])
            progress.log(
                f"azimuthal per-tier correction (Druckmullerova thesis 4.15): "
                f"the gain varies around the limb by up to {_wk[1]:.0f}% "
                f"({_exp_name(float(_wk[0]))} tier), which one scalar per tier "
                f"cannot express. Corrected on {len(_ldic)} of {len(secs)} "
                f"tiers, mean-preserving so the photometric chain and the "
                f"shared pedestal are untouched; applied continuously in "
                f"azimuth and faded out beyond "
                f"{_ldic_rout / max(_Rseed, 1e-6):.1f} R, where it was no "
                f"longer fitted.", None)
        else:
            progress.log("azimuthal per-tier correction: not enough overlap to "
                         "fit; tiers left as they are", None)
    except Exception as _e:
        _ldic = {}
        progress.log(f"azimuthal per-tier correction skipped ({_e})", None)

    # --- is any of the above COLOUR-dependent? ---
    # Diagnostic only: refits the per-tier scale, the shared pedestal and the
    # azimuthal gain per colour channel and records the spread. Nothing here
    # feeds the merge. See _per_channel_photometry.
    try:
        if os.environ.get("ECLIPSEFORGE_NO_CHANNEL_CHECK") != "1":
            _per_channel_photometry(
                stacks_bayer, sat_half, secs, ref, _links_good, cal, pedestal,
                cym, cxm, _Rseed, float(color_info["sat_level"]), progress,
                stats, wd=wd)
            _ch = (stats.get("per_channel_photometry") or {})
            _chd = _ch.get("channel_floor_deltas")
            _chi = _ch.get("channel_floors") or {}
            # VERIFIED END TO END by tools/smoke_channel_floors.py, which
            # plants a known floor in one channel and measures the colour of
            # the MERGED picture at 1.5 R and 3.2 R: a 100 ADU blue floor makes
            # a 70.8% radial cast and this removes 74% of it. The 26% left is
            # the ladder-degenerate component _fit_channel_floors discards on
            # purpose, so it is the expected floor of the method, not a defect.
            #
            # An earlier run of that same bench said the correction made things
            # TEN TIMES WORSE, and the bench was right to say so: the fixture
            # wrote cubes with no ROWORDER card, which tripped a real plane-axis
            # bug in fits._orient and swapped R and B before any of this ran.
            # Recorded because the wrong lesson was nearly drawn from it -- the
            # measurement was sound, the thing being measured was not.
            if (stats.get("fit_channel_floors") and _chd
                    and _chi.get("usable")
                    and os.environ.get("ECLIPSEFORGE_NO_CHANNEL_FLOORS") != "1"):
                ped_channel = np.asarray(_chd, np.float32)
                stats["channel_floors_applied"] = [float(x) for x in ped_channel]
                progress.log(
                    f"per-channel black level: R {_chd[0]:+.1f} G {_chd[1]:+.1f} "
                    f"B {_chd[2]:+.1f} ADU, subtracted on top of the shared "
                    f"pedestal. This format carries no per-channel black level "
                    f"in its header, so each plane's own floor survives into "
                    f"the merge; it is nothing against the inner corona and it "
                    f"is the whole signal in the outer field, which is what "
                    f"makes a colour cast that changes with radius. Fitted with "
                    f"the exposure ladder free, so an achromatic timing error "
                    f"cannot land here; the part indistinguishable from one "
                    f"({_chi.get('deltas_ladder_part')}) is discarded, not "
                    f"applied. Ladder deviation from the stated times: "
                    f"{_chi.get('ladder_dev_pct')} %.", None)
            elif stats.get("fit_channel_floors") and _chi.get("why"):
                progress.log(f"per-channel black level not fitted: "
                             f"{_chi['why']}", None)
    except Exception as _e:
        progress.log(f"per-channel photometry check skipped ({_e})", None)
    _links_good = None      # one half-res boolean per link; nothing needs them now

    stats["quality"] = {str(k): v for k, v in quality.items()}
    stats["tiers"] = [{"sec": float(s), "n": len(quality[s]["used"]),
                       "n_avail": len(tiers[s]),
                       "spread": quality[s]["spread"],
                       "best": quality[s]["best"],
                       "sharpness": float(max(quality[s]["scores"].values())),
                       "cal": float(cal[s]),
                       "shift": [float(abs_shift[s][0]), float(abs_shift[s][1])]}
                      for s in secs]
    # --- is the input actually scene-linear? ---
    #
    # The whole pipeline assumes signal proportional to exposure time. That is
    # true of a raw file and false of anything a raw developer has exported,
    # because a developer applies a tone curve. Inverting sRGB on load (see
    # TiffFrame) undoes the *encoding*; it does not undo Adobe's base curve, and
    # nothing can, because that curve is neither published nor invertible from
    # the file.
    #
    # WHAT DOES AND DOES NOT DETECT IT. The obvious test -- does the ratio
    # between two tiers change with brightness -- FAILS on the commonest case.
    # For a power law y = x^g the ratio is (s_b/s_a)^g at every level, so a
    # gamma is invisible to it; a synthetic check measured a spread of 1.011,
    # indistinguishable from linear. What a gamma does instead is tilt every
    # photometric LINK by the same amount, and the links are already measured
    # above. Fit g from them:
    #
    #     link_i = (s_i+1 / s_i)^(g-1)   ->   g = 1 + mean( ln link / ln step )
    #
    # Measured on the 250 mm test set, the same nine frames both ways:
    #
    #                      tier factors      fitted g   limb fit rms
    #   from the CR2s      0.944 .. 1.015      ~1.00       0.95 px
    #   Lightroom TIFFs    0.251 .. 2.328      ~1.49      10.32 px
    #
    # A factor of 9.3 across the bracket where the raws give 1.08, and a limb
    # that stops being a circle -- "looks like a potato".
    _gam = None
    try:
        _ln = [(ratios[i], np.log(secs[i + 1] / secs[i]))
               for i in range(n - 1) if ratios[i] != 0.0]
        _ln = [(lr, st) for lr, st in _ln if abs(st) > 1e-6]
        if len(_ln) >= 3:
            _gam = 1.0 + float(np.median([lr / st for lr, st in _ln]))
    except Exception:
        _gam = None
    _israw = not str(quality[secs[0]]["best"]).lower().endswith(
        (".tif", ".tiff", ".png", ".jpg", ".jpeg"))
    if _gam is not None:
        stats["linearity_gamma"] = round(_gam, 3)
        if abs(_gam - 1.0) > 0.08:
            progress.log(
                f"WARNING: these frames are NOT scene-linear. The photometric "
                f"links are tilted as if the data carried a gamma of about "
                f"{_gam:.2f} — on linear data every link sits at 1.000 by "
                f"construction and this comes out 1.00. "
                + ("A raw developer's tone curve does exactly this. Inverting "
                   "sRGB on load undoes the encoding, not Adobe's base curve, "
                   "and that curve cannot be inverted from the file: Lightroom "
                   "has no scene-linear TIFF export. Use the raw files. A flat "
                   "built from the same exports carries the same curve and does "
                   "not divide out either, which is why the falloff it reports "
                   "does not match the one measured from raws. "
                   if not _israw else
                   "On raw input this should not happen: check for a black "
                   "level left in the data, or a camera profile applied before "
                   "these files were written. ")
                + "Everything below — the photometric factors, the merge, the "
                  "limb fit — is computed as if the data were linear, so treat "
                  "this run as unreliable.", None)
    progress.log("photometric calibration: " +
                 ", ".join(f"{np.exp(l):.3f}" for l in logf), 0.53)
    span_ev = np.log2(max(secs) / min(secs))
    if span_ev < 6:
        progress.log(f"WARNING: the bracket spans only {span_ev:.1f} EV "
                     f"({min(secs):g}s to {max(secs):g}s). A totality corona "
                     "bracket normally spans 10-14 EV; this looks like a "
                     "partial-phase, diamond-ring or beads sequence, which this "
                     "pipeline is not built for.", None)
    # THE PER-LINK RESIDUALS, NOT ONLY THE CUMULATIVE FACTORS.
    #
    # `cal` is a running product down the chain, so one bad link and a hundred
    # slightly-biased ones look identical in it -- both just end up far from 1.
    # The links themselves tell them apart at a glance, and that distinction is
    # the whole diagnosis: scatter about 1.0 is noise and averages out; a
    # consistent lean to one side is a systematic error and COMPOUNDS.
    #
    # On the run that motivated this, 23 of 24 links sat below 1.0 with a
    # median of 0.798, which is not something anyone would read out of
    # "58.202, 29.812, 15.324, ...".
    _lnk = np.exp(np.asarray(ratios, np.float64)) if ratios else np.ones(0)
    # RECORDED, NOT ONLY LOGGED (0.22.79). These are the numbers that say
    # whether the exposure times predict the signal, and they existed only as a
    # line of text -- so nothing downstream, and no test, could look at them.
    stats["photometric_links"] = [round(float(x), 5) for x in _lnk]
    if _lnk.size:
        progress.log("photometric links (1.000 = the exposure ratio predicts "
                     "this tier exactly): " +
                     ", ".join(f"{x:.3f}" for x in _lnk), None)
        _meas = _lnk[np.abs(np.log(np.maximum(_lnk, 1e-9))) > 1e-9]   # drop fallbacks
        if _meas.size >= 6:
            _lean = float((_meas < 1.0).mean())
            _med = float(np.median(_meas))
            # 1.05 was too tight. The 0.16.1 regression put 12 of 13 links
            # above 1.000 with a median of 1.050 and this stayed silent, which
            # is the one case it most needed to speak up for. A consistent
            # DIRECTION across a dozen links is itself the evidence; the size
            # of the median only has to rule out rounding.
            if (_lean > 0.8 or _lean < 0.2) and abs(np.log(_med)) > np.log(1.02):
                progress.log(
                    f"WARNING: {int(round(_lean * _meas.size)) if _lean > 0.5 else int(round((1 - _lean) * _meas.size))}"
                    f" of {_meas.size} measured links lean the same way "
                    f"(median {_med:.3f}). That is a systematic error, not "
                    f"scatter, and it compounds down the chain to "
                    f"{_med ** _meas.size:.4f}. The usual cause is a black "
                    f"level left in the data — a pedestal does not scale with "
                    f"exposure, so it biases every link in the same direction. "
                    f"Check that the frames are black-subtracted.", None)
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
        # kept for the cross-check against the merged fit further down: how far
        # apart the tiers are is what says whether a given disagreement with the
        # merged fit is large. A flat percentage cannot know that.
        stats["R_consensus_spread"] = _sp
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
        # 226 px of spread on a 525 px Moon (the test set's 560mm run) is 43%: the
        # tiers are not looking at the same place. the 600 mm reference set sits at
        # 1.3%, the test set's own 250mm at 1.0% and his 360mm at 3.3%, so 8% is
        # comfortably clear of anything a working run produces.
        if Rmoon and span > 0.08 * Rmoon:
            progress.log(
                f"WARNING: the tiers' lunar limbs are spread over {span:.0f}px, "
                f"{100 * span / Rmoon:.0f}% of the lunar radius. They are the "
                f"same Moon seconds apart, so this is a cross-tier alignment "
                f"failure, not lunar motion. The merged limb is a smear of "
                f"{span:.0f}px and its fitted radius, the disc mask and every "
                f"radial filter inherit it.", None)
            stats["limb_spread_bad"] = round(float(span), 1)
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
        track = {}          # _moon_track can raise, and the weight trial below
                            # reads `track` -- that was a NameError caught as
                            # "merge weight trial skipped", i.e. a silent revert
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
            # Clears the MASK only. The exposure-exponent trial below is a
            # separate question and reads `track` directly -- it used to be
            # handed Rmoon as well and bailed out here. See _pick_weight_alpha.
            Rmoon = None
        # ...and, with the same half-res merges, how hard the long tiers should
        # be allowed to outvote the short ones. See _pick_weight_alpha.
        # THIRD BISECT SWITCH. ECLIPSEFORGE_WALPHA=1.0 holds the exponent at 1.0
        # and skips the trial. The trial is a 0.22.5-era feature and it fires on
        # some brackets and not others -- on the 600 mm reference set it picks 0.55 and
        # reports "+106% coherent detail at 1.02-1.12 R", which is its purpose:
        # it tilts weight toward the SHORT tiers, which are also the noisiest.
        # Its guard only checks that the mid and outer shells keep their radial
        # coherence, so a fine radial texture near the limb is exactly what it
        # is allowed to add. That makes it the first thing to test once the two
        # changes from today are excluded, and it explains why the same texture
        # is already present in his own 0.22.5 export.
        _wenv = os.environ.get("ECLIPSEFORGE_WALPHA")
        if _wenv:
            try:
                _walpha = float(_wenv)
                progress.log(f"merge weight trial SKIPPED — exposure exponent "
                             f"held at {_walpha:.2f} by ECLIPSEFORGE_WALPHA",
                             None)
                stats["merge_weight"] = {"alpha": _walpha, "forced": True}
            except ValueError:
                _wenv = None
        if not _wenv:
            try:
                _walpha, _winfo = _pick_weight_alpha(
                    stacks_half, sat_half, secs, cal, abs_shift, track,
                    progress)
                stats["merge_weight"] = _winfo
            except Exception as e:
                progress.log(f"merge weight trial skipped ({e})", None)

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
    # THE TRIAL ADVISES; IT DOES NOT DECIDE.
    #
    # 0.22.28 through .30 tried to pick the feather automatically from a
    # half-resolution reconstruction of the merge. Three releases, and it still
    # read 99% on a bracket the offline bench measures at 13% -- the same
    # function, run on arrays rebuilt from that dataset's own exported tiers,
    # returns 0.33, so the estimator is not wrong about the maths, it is wrong
    # about what `stacks_half` and `sat_half` actually contain at this point.
    # Until that is understood the number has not earned a decision.
    #
    # So the feather is a SETTING (GUI: "Merge weight"), and the trial prints
    # what it measured next to it, with its own disagreement stated. A control
    # the user can see and change in a second beats a guess that silently picks
    # wrong -- which is what shipped a pink rim on the test set's data twice.
    _fm = (os.environ.get("ECLIPSEFORGE_FEATHER") or feather or
           _FEATHER_DEFAULT).lower()
    if _fm not in ("plain", "taper", "masked"):
        _fm = _FEATHER_DEFAULT
    stats["feather_mode"] = _fm
    _fratio = float("nan")
    try:
        _adv, _fratio = _pick_feather(
            stacks_half, sat_half, secs, cal, pedestal,
            float(color_info["sat_level"]), cym, cxm, _Rseed, _walpha,
            _feather / 2.0, progress, advise_only=True)
        if np.isfinite(_fratio):
            stats["feather_ratio"] = round(_fratio, 4)
            progress.log(
                f"merge weight: {FEATHER_NAMES.get(_fm, _fm)} (setting). "
                f"The trial measures the "
                f"plain feather at {100 * _fratio:.0f}% of the leak-free level "
                f"at 1.02 R and would suggest '{_adv}' — ADVICE ONLY, and this "
                f"estimator has been wrong before: it read 99% on a bracket "
                f"the offline bench puts at 13%. Trust the picture, not this "
                f"number.", None)
        else:
            progress.log(f"merge weight: {FEATHER_NAMES.get(_fm, _fm)} "
                         f"(setting); the trial could "
                         f"not measure a ratio on this bracket", None)
    except Exception as _e:
        progress.log(f"merge weight: {FEATHER_NAMES.get(_fm, _fm)} (setting); "
                     f"trial skipped ({_e})",
                     None)
    # DO NOT write the choice back into os.environ. 0.22.28 did, so that the
    # merge loop could pick it up -- and the app is ONE LONG-RUNNING PROCESS.
    # The first folder's trial set the variable, every later folder in the same
    # session read it at startup, skipped its own trial, and reported "feather
    # held at 'plain' by ECLIPSEFORGE_FEATHER". the 600 mm reference set chose plain,
    # and the test set's 360 mm -- which needs the leak-free weight -- then inherited
    # it and came out with the pink rim the trial exists to prevent.
    # `_fm` is passed to _feather_weight explicitly below, so nothing needs the
    # environment; writing to it only leaked state between runs.

    acc = np.zeros((H2, W2, 3), np.float32)
    wsum = np.zeros((H2, W2), np.float32)
    _angmap = None
    # Pixels the SHORTEST tier already cannot hold. If the longest tier cannot
    # hold them either then no tier can, and without a fallback they leave the
    # merge as exact zeros -- black. See the fill below the loop.
    _fill_idx = None
    _fill_rgb = None
    for k, s in enumerate(secs):
        rgb = _demosaic(stacks_bayer[s], demosaic_method)
        # cmax is the CLIPPING test, so it has to be measured in raw units,
        # before white balance. WB runs at G=1, so the red gain (~2.1x) used to
        # push a red pixel past 0.97*sat_level while its photosite sat at only
        # 0.45 of saturation: every tier declared prominence H-alpha clipped a
        # full stop early, the shortest one included, and with no tier holding
        # weight there the merge filled prominence cores by leakage from tens
        # of px away. Neutral pixels are unaffected either way -- for a neutral
        # subject green is the largest raw channel and WB leaves it alone --
        # which is why this never showed on the corona itself.
        # THE CLIPPING TEST IS ASKED OF THE MOSAIC, NOT OF THE DEMOSAICED
        # RESULT. See raw.cfa_clip_max: two of every pixel's three channels are
        # interpolated through 5x5 kernels with negative lobes, so a pixel next
        # to a saturated photosite could come out below threshold and enter the
        # merge unflagged, carrying a value partly reconstructed from a
        # photosite that hit the ceiling.
        cmax = cfa_clip_max(stacks_bayer[s])
        if _fsat is not None:
            cmax *= _fsat            # back to raw units -- see _fsat above
        _nclip_before = float((rgb.max(axis=2) *
                               (_fsat if _fsat is not None else 1.0)
                               > 0.97 * sat_level).mean())
        _nclip_after = float((cmax > 0.97 * sat_level).mean())
        if _nclip_after > _nclip_before:
            stats.setdefault("cfa_clip_extra", {})[f"{s:g}"] = \
                round(100.0 * (_nclip_after - _nclip_before), 3)
        # The shared pedestal comes off HERE: after demosaic, before white
        # balance. A black-level residual is one number per photosite, so it is
        # the same in all three channels at this point and becomes three
        # different numbers the moment the WB gains are applied. Subtracting it
        # before the clipping test would be wrong too -- cmax is a question
        # about the sensor's ceiling, not about the scene.
        #
        # KNOWN APPROXIMATION, WITH A FLAT. The flat is divided out before
        # stacks_half is built, so a pedestal that is genuinely constant in RAW
        # units arrives here already divided by the flat, and both the fit and
        # this subtraction treat it as constant instead. The error is the flat's
        # own departure from unity -- a few percent of a few ADU, so hundredths
        # of an ADU, against an effect of 3-6. Not worth fitting in the raw
        # domain, but worth writing down: the one dataset here that carries a
        # master flat is also the one whose fitted pedestal is consistent with
        # zero, so this has never been exercised on data that needs it.
        if _hillK is not None:
            # Hill's per-channel affine, in the raw channel domain and before
            # white balance -- the same place the shared pedestal came off, and
            # for the same reason: this is where a black-level residual is one
            # number per photosite. The white balance and the camera matrix are
            # both linear, so an affine applied here stays affine after them;
            # the matrix mixes the three offsets, which is correct, because an
            # offset in a camera channel IS an offset in a mixture of the
            # output channels.
            if pedestal:
                rgb -= np.float32(pedestal)     # the common part, as always
            if ped_channel is not None:
                # Exact after demosaic, not an approximation: Malvar's
                # gradient-correction terms are differences of neighbouring
                # samples, so they vanish on a field that is constant within
                # each channel. A constant per-channel offset therefore comes
                # out of the demosaic as the same constant per channel.
                rgb -= ped_channel[None, None, :]
            rgb /= np.float32(s)
            rgb *= _hillK[s][None, None, :]     # then the per-pair transform
            rgb += _hillQ[s][None, None, :]
            rgb *= wb[None, None, :]
            rgb = (rgb.reshape(-1, 3) @ cam2rgb.T).reshape(H2, W2, 3)
        else:
            if pedestal:
                rgb -= np.float32(pedestal)
            if ped_channel is not None:
                rgb -= ped_channel[None, None, :]   # see the note above
            rgb *= wb[None, None, :]
            rgb = (rgb.reshape(-1, 3) @ cam2rgb.T).reshape(H2, W2, 3)
            rgb /= np.float32(s * cal[s])
        ady, adx = abs_shift[s]
        ady, adx = 2 * ady, 2 * adx
        # CUBIC, NOT BILINEAR. The shifts are sub-pixel -- they come out of a
        # least-squares network, so they are fractional even though every
        # correlation peak was integer -- which means every tier but the
        # reference is RESAMPLED, and a bilinear resample is a low-pass filter
        # applied before any detail layer ever sees the data.
        #
        # Measured as the amplitude a sinusoid keeps through one sub-pixel
        # shift, averaged over shifts of 0.25, 0.5 and 0.75 px:
        #
        #     feature size   order=1   order=3   order=5
        #        1.4 px       0.698     0.962     0.994
        #        2.8 px       0.562     0.877     0.962
        #        4.0 px       0.763     0.981     0.998
        #        8.0 px       0.937     0.999     1.000
        #
        # 1.4 px is what this image resolves and the finest scale MGN runs at,
        # so bilinear was discarding 30% of it, and 44% at 2.8 px, on eleven of
        # the twelve tiers. That is a lot of the "why is there no fine
        # structure" question answered before any filter is blamed for it.
        #
        # THE COST, both measured. Cubic overshoots a hard step edge by 10%
        # (order=5 by 12%), and the lunar limb is a hard step edge -- but the
        # disc mask sits at R + 8 px and the overshoot is one pixel wide, so it
        # lands inside the mask; saturated regions are already excluded by the
        # mosaic clipping test and its collar. And it is slower: 8.1 s against
        # 1.3 s per full-frame channel here, about five minutes on a twelve-tier
        # 8000 x 5400 stack.
        #
        # The masks below stay at order=1 on purpose: they are validity, not
        # signal, and ringing a mask invents partial validity at every edge.
        for c in range(3):
            rgb[:, :, c] = ndimage.shift(rgb[:, :, c], (ady, adx), order=3,
                                         mode="nearest")
        # LDIC's k_i(phi), q_i(phi) -- applied AFTER the shift, so the azimuth
        # map is the common frame's. Mean-preserving (see
        # _fit_azimuthal_affine), so this only removes the variation around the
        # limb that a scalar cal[s] cannot carry. Same correction on all three
        # channels: k is a gain on the incident light and q is a black-level
        # residual, which is one number per photosite at this point -- exactly
        # the argument that puts the shared pedestal before white balance.
        if _ldic.get(s) is not None:
            if _angmap is None:
                # Azimuth itself, not a segment index. 0.22.26 built a segment
                # INDEX here and looked k and q up in a 60-entry table, which
                # made a smooth function into a 60-step staircase: a hard gain
                # discontinuity every 6 degrees, running from the disc to the
                # frame edge. On the 250 mm test set, where k spans 82%, the
                # steps between neighbouring segments reach ~17% and the whole
                # frame filled with radial wedges. Reported by a tester on 0.22.32.
                _yy = np.arange(H2, dtype=np.float32)[:, None] - 2.0 * cym
                _xx = np.arange(W2, dtype=np.float32)[None, :] - 2.0 * cxm
                # A 16384-entry index, not a per-pixel trig evaluation: eight
                # cos/sin passes over an 8100x5357 frame, once per tier, is
                # gigabytes of temporaries for nothing. At this table size the
                # residual step is 0.06% of k, three orders below the 17% the
                # 60-segment lookup was making.
                _angmap = ((np.arctan2(_yy, _xx) % (2 * np.pi))
                           / (2 * np.pi) * _LDIC_TABLE).astype(np.uint16)
                np.minimum(_angmap, _LDIC_TABLE - 1, out=_angmap)
                # ... and the radius, to keep the correction inside the band it
                # was actually fitted on (1.0-2.5 R). Beyond that it was pure
                # extrapolation applied at full strength out to the corners,
                # which is what put broad fans into the outer field even where
                # the steps were invisible. Cosine taper to identity over the
                # half-radius above the fit's outer edge.
                _rmapf = np.sqrt(_yy * _yy + _xx * _xx, dtype=np.float32)
                del _yy, _xx
                _r_in = np.float32(2.0 * _ldic_rout)          # half-res -> full
                _r_out = np.float32(2.0 * _ldic_rout + _Rseed)
                _t = np.clip((_rmapf - _r_in) / max(float(_r_out - _r_in), 1e-6),
                             0.0, 1.0)
                _ldic_amp = (0.5 * (1.0 + np.cos(np.pi * _t))).astype(np.float32)
                del _t
                # AND THE SAME AT THE INNER EDGE, which this did not have.
                #
                # The fit's band is 1.0-2.5 R. The taper above stopped the
                # extrapolation past 2.5 R and nothing stopped it below 1.0 R,
                # so the full correction ran all the way to r = 0 -- straight
                # through the occulted disc, where there is no signal at all.
                # `q` is an OFFSET, so what it writes there is not a small
                # correction to something: it IS the something. rgb = 0*k + q
                # leaves q(phi) standing on its own, which prints as a
                # pinwheel centred on the Moon, and an order-4 trigonometric
                # polynomial is exactly what a pinwheel looks like.
                #
                # Measured on a third tester's exported 1/4 s tier, azimuthal swing
                # of the median in 24 sectors: 60% at 0.2-0.6 R and 55% at
                # 0.6-0.95 R, inside a disc that should be featureless. The
                # composite hides it under the disc mask; an exported tier does
                # not, and neither does anything that samples across the limb.
                #
                # Faded to identity BELOW 1.0 R rather than at it, so the
                # corona from the limb outward -- the part the fit is actually
                # for, and where the detail layers are most sensitive -- is
                # untouched.
                _ri_hi = np.float32(2.0 * _Rseed)                 # 1.0 R
                _ri_lo = np.float32(2.0 * _Rseed * 0.85)          # 0.85 R
                _ti = np.clip((_ri_hi - _rmapf) / max(float(_ri_hi - _ri_lo), 1e-6),
                              0.0, 1.0)
                _ldic_amp *= (0.5 * (1.0 + np.cos(np.pi * _ti))).astype(np.float32)
                del _rmapf, _ti
            _kk, _qq, _ck, _cq = _ldic[s]
            _ktab, _qtab = _ldic_tables(_kk, _qq, _ck, _cq)
            # gain first, then offset, each freed before the next -- these are
            # full-frame float32 arrays and the merge loop is already the
            # memory high-water mark of the run
            _kf = _ktab[_angmap]
            _kf -= np.float32(1.0)             # taper to identity (k=1, q=0)
            _kf *= _ldic_amp                   # outside the fitted band
            _kf += np.float32(1.0)
            rgb *= _kf[:, :, None]
            del _kf
            _qf = _qtab[_angmap]
            _qf *= _ldic_amp
            rgb += _qf[:, :, None]
            del _qf
        # THE MERGE WEIGHT NEEDS BOTH ENDS, NOT JUST THE TOP.
        #
        # This was a high-end shoulder alone: full weight for every pixel below
        # the knee, including one holding nothing but read noise. Every source
        # says otherwise. Druckmuller, Rusin & Minarovjech 2006, requirement (b)
        # on p.134: "w = 0 in substantially underexposed OR overexposed parts of
        # the corona". Druckmullerova's thesis p.48: the weight is "equal to one
        # on a majority of the dynamic range, equal to zero in the highest AND
        # THE LOWEST part of the dynamic range, with gradual and continuous
        # transitions", because "the lowest part of the dynamic range contains
        # mostly noise. Adding it to the composed image would only increase its
        # noise". Hill's slide shows the same trapezoid.
        #
        # It matters here more than it would elsewhere because _pick_weight_alpha
        # settled on alpha = 0.55, and for s < 1 that gives the SHORT tiers MORE
        # relative weight than alpha = 1 would -- deliberately, since it buys
        # +106% coherent detail at the limb. The cost was paid in the outer
        # field, where those same tiers hold no signal at all. A per-pixel floor
        # is what lets the exponent buy the limb without paying out there; one
        # global exponent cannot separate the two.
        #
        # MEASURED by rebuilding the merge from the reference set's exported tiers (all
        # 14) and looking at azimuthal scatter about a smooth profile, which in
        # a shell that should be smooth is noise:
        #
        #   shell        as shipped   floor 0.002  floor 0.005  floor 0.01
        #   1.2-1.8 R      0.3189       0.3189       0.3189      0.3190
        #   1.8-2.6 R      0.1162       0.1150       0.1140      0.1129
        #   2.6-3.4 R      0.0299       0.0289       0.0282      0.0276
        #   3.4-4.2 R      0.0135       0.0127       0.0123      0.0119
        #
        # Monotone, nothing lost where the corona is bright, and the gain grows
        # with radius exactly as the mechanism predicts: -2.8%, -7.7%, -11.9%.
        #
        # CAVEAT, stated because it is not a clean test: that reconstruction ran
        # on the sRGB-encoded 8-bit-decoded exports at quarter resolution, not on
        # the real linear merge, and its "noise" includes any real azimuthal
        # structure finer than the 9-bin smooth. The DIRECTION is supported by
        # three independent sources and by the measurement; the exact magnitude
        # is not something this test can pin down. Hence a conservative floor.
        knee = 0.87 * sat_level
        wsat = 0.5 * (1.0 + np.tanh((knee - cmax) / (0.06 * sat_level)))
        # Kept separately, and shifted with the weight, because the feather
        # below must not carry weight across it. See _feather_weight.
        _valid = (cmax <= 0.97 * sat_level).astype(np.float32)
        wsat[cmax > 0.97 * sat_level] = 0.0
        _lo = _MERGE_FLOOR * sat_level
        if _lo > 0:
            wsat *= 0.5 * (1.0 + np.tanh((cmax - _lo) / (0.5 * _lo)))
        wsat = ndimage.shift(wsat, (ady, adx), order=1, mode="nearest", cval=0)
        _valid = ndimage.shift(_valid, (ady, adx), order=1, mode="nearest",
                               cval=0)
        # s**alpha, not s: alpha is 1.0 unless the trial above measured that
        # tilting toward the shorter, sharper tiers buys limb detail without
        # costing the outer field. See _pick_weight_alpha.
        _wf = _feather_weight(wsat, _valid, _feather, _fm)
        w = np.float32(s ** _walpha) * _wf
        mw = moon_weight(s)
        if mw is not None:
            w *= mw
        if k == 0:
            # Remember where the shortest tier's weight is EXACTLY zero, and
            # this tier's value there. Exact zero happens only inside a tier's
            # clipping contour (and, with the taper, one feather-width inside
            # it); everywhere else the weight is small but positive. The lunar
            # disc is excluded -- every tier's weight is zero there by design
            # and that region is meant to stay black.
            _hz = _wf <= 0.0
            if mw is not None:
                _hz &= mw > 0.0
            _fill_idx = np.nonzero(_hz.ravel())[0]
            _fill_rgb = (rgb.reshape(-1, 3)[_fill_idx].copy()
                         if _fill_idx.size else None)
            del _hz
        del _wf
        if mw is not None:
            del mw
        if export_tiers:
            _write_tier_tiff(folder, s, cal[s], sat_level, _tier_hi,
                             rgb[crop_origin[0]:crop_origin[0] + _kh,
                                 crop_origin[1]:crop_origin[1] + _kw],
                             len(tiers[s]), tier_linear, progress)
        acc += w[:, :, None] * rgb
        wsum += w
        del rgb, w, wsat, cmax, _valid
        progress.log(f"merged tier {s:g}s", 0.55 + 0.3 * (k + 1) / n)
    # WHERE NO TIER HAS ANY WEIGHT, THE MERGE USED TO PRODUCE BLACK.
    #
    # `hdr = acc / max(wsum, 1e-9)` is 0/0 at every pixel that is clipped in all
    # tiers, and 0 is not a plausible value for the brightest thing in the frame
    # -- it is a hole. hit on the 250 mm test set in 0.22.35: the core
    # of a prominence came out as a black sliver inside the prominence.
    #
    # It appeared with 0.22.35 and was not caused by it. Measured on that set's
    # shortest tier, _MG_4633.CR2 (1/125 s), against the two saturation levels:
    #
    #     sat_level 13978 (<= 0.22.34)   threshold 13558 ADU        0 photosites
    #     sat_level 12925 (0.22.35)      threshold 12537 ADU      227 photosites
    #
    # The data on that body plateaus at 13255 ADU above black, so with the old
    # container-ceiling saturation level the clipping test could not fire AT ALL
    # -- `_valid` was 1 everywhere, the taper never tapered, and no pixel could
    # lose all its weight. The hole is what a working clipping test looks like
    # when the bracket does not reach short enough to hold the prominence core.
    #
    # Why only the Photometric ('taper') and 'masked' weights show it: those
    # multiply by `_valid` and are therefore exactly zero on a clipped pixel.
    # The Detail weight ('plain') is a plain Gaussian blur with no such factor,
    # so weight leaks in from unclipped neighbours and the hole fills itself --
    # with clipped, under-reporting data. That is not the weight being right; it
    # is the same leak that dragged the inner corona down and printed the pink
    # rim, doing something cosmetically useful for once.
    #
    # The fallback: use the SHORTEST tier's value, which is the least saturated
    # estimate that exists. It is a lower bound, not a measurement -- the true
    # value is somewhere above it -- but a lower bound in a prominence core
    # reads as a blown highlight, which is honest, while zero reads as a hole
    # that no amount of stretching can repair. The fix for the picture is a
    # shorter tier at the top of the bracket; the report says so.
    _nfill = 0
    if _fill_idx is not None and _fill_idx.size and _fill_rgb is not None:
        _sel = wsum.ravel()[_fill_idx] <= 0.0
        _nfill = int(_sel.sum())
        if _nfill:
            _hidx = _fill_idx[_sel]
            acc.reshape(-1, 3)[_hidx] = _fill_rgb[_sel]
            wsum.ravel()[_hidx] = 1.0
            del _hidx
        del _sel
    del _fill_idx, _fill_rgb
    stats["all_clipped_px"] = _nfill
    if _nfill:
        progress.log(
            f"{_nfill} px are saturated in EVERY tier, the shortest "
            f"({_exp_name(secs[0])}) included — filled from that tier instead "
            f"of being left at zero, which printed as a black hole inside the "
            f"prominence. Those pixels are a LOWER BOUND on the real "
            f"brightness; only a shorter exposure can measure them", None)
    hdr = acc / np.maximum(wsum[:, :, None], 1e-9)
    # (_kh, _kw) can differ from (H2, W2) with crop_origin still (0,0) when the
    # trim is bottom/right only. The tier TIFFs were already being sliced in
    # that case while the render was not, so the two came out different sizes
    # and the render kept a band of edge-replicated pixels.
    if crop_origin != (0, 0) or (_kh, _kw) != (H2, W2):
        hdr = hdr[crop_origin[0]:crop_origin[0] + _kh,
                  crop_origin[1]:crop_origin[1] + _kw]
    del acc, wsum
    # A RADIANCE CANNOT BE NEGATIVE (0.22.79). Nothing in the merge itself can
    # produce one -- the weights and the data are non-negative -- so a negative
    # channel here always means an OFFSET was subtracted that should not have
    # been, and the only offsets in the chain are the shared pedestal and the
    # linear fit's per-channel q. It is worth saying out loud rather than
    # clipping in silence: on the test set's 560 mm set this was 100% of the red
    # channel beyond 4 R, and it printed as a cyan ring nobody could explain
    # from the log.
    _neg = hdr < 0
    _nfrac = float(_neg.mean())
    if _nfrac > 1e-4:
        _worst = int(np.argmax([float(_neg[:, :, c].mean()) for c in range(3)]))
        progress.log(
            f"WARNING: {100 * _nfrac:.1f}% of the merged samples are NEGATIVE "
            f"(worst channel {'RGB'[_worst]}, {100 * _neg[:, :, _worst].mean():.0f}% "
            f"of it). A merge of non-negative data cannot do that, so an offset "
            f"has been over-subtracted — the shared pedestal ({pedestal:+.2f} "
            f"ADU) or, if the linear-fit chain is on, one of its per-channel "
            f"offsets. Clamped at zero so the picture is not a colour cast, but "
            f"the photometry in that region is not trustworthy.", None)
        stats["merge_negative_frac"] = round(_nfrac, 5)
    np.clip(hdr, 0.0, None, out=hdr)
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
        rgb = _demosaic(stacks_bayer[s], demosaic_method)
        cmax = cfa_clip_max(stacks_bayer[s])   # see the merge loop
        if _fsat is not None:
            cmax *= _fsat
        rgb *= wb[None, None, :]
        rgb = (rgb.reshape(-1, 3) @ cam2rgb.T).reshape(H2, W2, 3)
        lt = (0.2126 * rgb[:, :, 0] + 0.7152 * rgb[:, :, 1] + 0.0722 * rgb[:, :, 2])
        del rgb
        ady, adx = abs_shift[s]; ady, adx = 2 * ady, 2 * adx
        # same interpolation as the merge, or the exported tiers are a softer
        # picture than the composite they are supposed to explain
        lt = ndimage.shift(lt / (s * cal[s]), (ady, adx), order=3, mode="nearest")
        wsat = 0.5 * (1.0 + np.tanh((0.82 * sat_level - cmax) / (0.08 * sat_level)))
        _lvalid = (cmax <= 0.97 * sat_level).astype(np.float32)
        del cmax
        wsat = ndimage.shift(wsat, (ady, adx), order=1, mode="nearest", cval=0)
        _lvalid = ndimage.shift(_lvalid, (ady, adx), order=1, mode="nearest",
                                cval=0)
        # same leak as the main merge, milder because this wsat has no hard
        # zero -- its tanh is already ~0.01 at saturation -- but a blown tier
        # still has no business contributing here either
        w = np.float32(s) * _feather_weight(wsat, _lvalid, 1.2 * _feather)
        mw = moon_weight(s)
        if mw is not None:
            w *= mw
            del mw
        accn = lt * w if accn is None else accn + lt * w
        accw = w if accw is None else accw + w
        del lt, w, wsat, _lvalid
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
        # THE ACCEPTANCE TEST ABOVE ALLOWS rms < 0.08 R, WHICH IS VERY LOOSE.
        # the test set's 560mm fit came in at 32.77 px on R=587 -- 5.6%, comfortably
        # accepted, and the circle it describes is nothing like the Moon. The
        # three usable datasets sit at 0.41%, 0.30% and 1.16%, so there is a
        # wide empty band between "a real limb" and what the filter admits.
        # Warn in it rather than silently taking the fit.
        if R > 0 and rms > 0.02 * R:
            progress.log(
                f"WARNING: the limb fit's rms is {rms:.1f}px, {100 * rms / R:.1f}% "
                f"of the fitted radius. A real lunar limb fits to well under 1%; "
                f"this circle does not describe an edge. Usually the merged limb "
                f"is smeared by a cross-tier alignment error. The disc mask and "
                f"every radial filter are built on this circle.", None)
            stats["limb_fit_rms_bad"] = round(float(rms), 2)
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
    # A FLAT 15% IS THE WRONG TEST, and a real run showed why.
    #
    # a second tester's 360 mm set (0.15.1): the tiers agreed on R = 456 px with a
    # SPREAD OF 2 px, and the merged fit came out 470.2 px -- 3.1%, so this said
    # nothing. But against a 2 px spread, 14 px is seven sigma. The merged limb
    # ramp was 25 px where his other set's was 10, the alignment residual 7.9 px
    # full-res, and the disagreement rim 134 px wide: the merge was smeared, the
    # 50% crossing sat outside the real limb, and R came back too big.
    #
    # That matters far beyond the mask. R and the limb profile set the radial
    # profile MGN divides out, FNRGF's rings, the deband, and the Pellett
    # geometry. Build them on a circle 14 px too large and every one of them
    # prints concentric arcs -- which is exactly the ringing he reported in MGN,
    # FNRGF and NAFE at once. One wrong circle, three layers showing it.
    #
    # So the test is now "large compared with how well the tiers agree", with a
    # 1%-of-R floor so a very tight spread cannot make it fire on nothing.
    # Against the three real datasets available:
    #
    #   set                consensus  spread  merged fit  |diff|  fires?
    #   600 mm ref            617 px    7 px    619.1 px   2.1 px   no
    #   250 mm test set         298 px    7 px    300.2 px   2.2 px   no
    #   360 mm test set         456 px    2 px    470.2 px  14.2 px   YES
    #
    _rc = stats.get("R_consensus")
    _rcs = float(stats.get("R_consensus_spread") or 0.0)
    # R_consensus_spread is a p90-p10 RANGE, not a standard deviation, and the
    # test below wants sigma. Treating the range as sigma made the threshold
    # 2.6x too high, which is how the 560 mm 2024 test set slipped through: 15
    # px range, 28 px disagreement -- under two of the real sigma, but the code
    # asked for four of the range and never fired. He reported it as "the moon
    # mask is too large and covers prominences", which is exactly what a limb 28
    # px large does. For a normal distribution p90-p10 = 2.563 sigma.
    _sig = _rcs / 2.563
    if _rc and R > 0 and abs(R - _rc) > max(4.0 * _sig, 0.01 * _rc) \
            and abs(R - _rc) / _rc <= 0.15:
        progress.log(
            f"WARNING: the merged limb fit says R={R:.0f}px but the tiers agree "
            f"on R={_rc:.0f}px to within {_rcs:.0f}px -- a {abs(R - _rc):.0f}px "
            f"disagreement, and the fit is the one that runs LARGE. "
            f"The 50% crossing sits outside the true limb whenever the edge is "
            f"soft, and the merged limb ramp above says how soft: across four "
            f"real datasets the bias tracks the ramp and nothing else "
            f"(0.3%/8px, 0.8%/9px, 3.2%/21px, 11.9%/69px), matching neither "
            f"the alignment residual nor the tier disagreement. So a wide ramp "
            f"means either a merge smeared by misalignment OR a genuinely soft "
            f"limb -- focus, seeing, or a slow lens -- and the two are told "
            f"apart by the alignment residual and per-tier limb spread above. "
            f"Either way R sets the radial profile MGN divides out, FNRGF's "
            f"rings and the deband, so a circle this size prints concentric "
            f"arcs in all of them, and the disc mask covers real corona.", None)
        stats["limb_fit_disputed"] = round(float(abs(R - _rc)), 1)
        # ...and now correct it, instead of only complaining.
        #
        # Five real datasets, every one biased the same way: the merged
        # half-level fit comes back LARGER than the tiers' own consensus, by an
        # amount that tracks the merged limb ramp and nothing else.
        #
        #   set                 consensus  spread   merged fit   bias   ramp
        #   600 mm ref             617 px    7 px      619.1     +2.1     8
        #   250 mm test set          298 px    7 px      300.3     +2.3     9
        #   360 mm test set          456 px    2 px      470.4    +14.4    21
        #   560 mm 2024 test set     525 px   15 px      553.0    +28.0    28
        #   a second tester 560mm (2024)   525 px  226 px      587.3    +62.3    69
        #
        # The cause is not in dispute: the 50% crossing between disc and
        # near-limb corona sits OUTSIDE the true limb whenever the edge is soft,
        # and the merged edge is soft because the Moon moves against the corona
        # during the bracket. The per-tier fits use the same half-level method
        # but each on a single unsmeared tier, so they carry only their own much
        # narrower ramp -- which is why their median is the better estimate of
        # where the limb actually is.
        #
        # WHY AN OFFSET AND NOT A RESCALE. Both the lunar relief the profile
        # carries and the soft-edge bias are a fixed number of PIXELS, not a
        # fraction of R. Scaling the profile by _rc/R would shrink the relief
        # along with the radius and quietly flatten a real measurement; adding a
        # constant moves the circle and leaves the per-azimuth shape -- the only
        # thing the merged fit measures better than the tiers -- exactly intact.
        #
        # Applied only when the fit runs LARGE. A merged fit that came out SMALL
        # would be the expected direction for a brighter combined corona (the
        # crossing moves inward as the corona brightens, see the per-tier note
        # above), so there is nothing to correct there.
        if fit is not None and R > _rc:
            _off = float(_rc) - float(R)
            limb_prof = np.asarray(limb_prof, np.float32) + _off
            R = float(_rc)
            rms = float(np.hypot(rms, _sig))
            stats["limb_fit_corrected_px"] = round(_off, 1)
            progress.log(
                f"limb fit corrected {_off:+.1f}px to the tiers' consensus "
                f"R={R:.0f}px, keeping the per-azimuth shape. The merged "
                f"half-level crossing runs large on a soft edge; across five "
                f"real datasets the bias is always positive and tracks the "
                f"merged limb ramp. This is the circle MGN's radial profile, "
                f"FNRGF's rings, the deband and the disc mask are all built "
                f"on, so it is corrected here rather than left to the disc "
                f"mask trim slider.", None)
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
    # Before the sky fit, because that fit is per channel: run it on raw sensor
    # colour and it absorbs the sensor's own R:G:B imbalance into what it
    # reports as the sky's colour gradient.
    if stats.get("no_camera_wb"):
        try:
            neutralise_corona_colour(wd, cyf, cxf, R, stats, progress)
        except Exception as e:
            progress.log(f"corona white balance skipped ({e})", None)
    # now that the corona's extent is measured, the sky beyond it can be fitted
    try:
        remove_sky_gradient(wd, cyf, cxf, R, stats.get("corona_extent_R"),
                            stats, progress)
        lum = np.load(os.path.join(wd, "hdr_lum.npy"))
    except Exception as e:
        progress.log(f"sky gradient removal skipped ({e})", None)
    # AFTER the sky fit, for the same reason that fit runs after the white
    # balance: the renderer applies this to the picture the sky fit produced,
    # so the profile has to describe that picture and not the one before it.
    try:
        measure_radial_colour(wd, cyf, cxf, R, stats, progress)
    except Exception as e:
        progress.log(f"radial colour skipped ({e})", None)

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
                "alignment quality: limb variance %.3f, rim %s, "
                "merged limb 20-80%% %.1f px"
                % (aq.get("cov_limb", float("nan")),
                   # "rim nan px" is what a GOOD result printed: with the tiers
                   # agreeing (cov 0.036 on the test set's 360mm) there is no
                   # disagreement rim to measure a width from, and the nan was
                   # the honest answer wearing an alarming face. Say it.
                   ("%.0f px" % aq["rim_width_px"]
                    if np.isfinite(aq.get("rim_width_px", float("nan")))
                    else "none detectable"),
                   aq.get("limb_width_med", float("nan"))), None)
            # TIERS THAT DISAGREE IN VALUE AT THE LIMB ARE A SECOND, SEPARATE
            # CAUSE OF RINGING -- and this number is the one that finds it.
            #
            # a second tester's 250 mm set rings, and none of the alignment guards
            # fire on it: network residual 0.66 px, limb spread 3 px, limb-fit
            # rms 0.95 px, track scatter 0/0. Geometrically it is a clean run.
            # What it has is limb variance 0.793 against 0.075, 0.067 and 0.052
            # on the other three real datasets -- ten times any of them.
            #
            # Measured on the cached layers, ring by ring, as the oscillation of
            # each layer's radial median about a 15-ring smooth, in percent of
            # that layer's own trend:
            #
            #   band          MGN            FNRGF          NAFE           inner
            #   1.02-1.15 R   19.4 / 9.9     10.9 / 2.5     13.9 / 2.7     15.1 / 4.6
            #   1.15-1.40 R    4.5 / 3.7      1.7 / 2.0      2.1 / 0.7      3.2 / 1.5
            #   1.40-1.80 R    3.2 / 3.0      0.6 / 1.5      0.9 / 0.6      1.9 / 1.3
            #                 (250 mm test set / 600 mm ref)
            #
            # The ringing lives in a band 1.02-1.15 R wide and is 2 to 5.3x his,
            # while the merged luminance it is built from oscillates only 1.4%
            # there. So the filters are amplifying a real tier disagreement by
            # about ten, not inventing it.
            #
            # NAFE is the tell: it is 5.3x worse and it is the one layer that
            # does not use the limb fit or the disc mask at all. So this is not
            # geometry.
            #
            # WITHDRAWN, 0.22.17: the conclusion drawn from that -- veiling
            # glare, "a 240 mm f/5.6 zoom flares far more than a 600 mm prime"
            # -- rested on cov_limb 0.793 for that set. Rebuilding the same
            # statistic from that run's own exported tiers gives 0.021. Its
            # tiers agree in the near-limb rim to 0.3-2.1% when compared ring by
            # ring, which is BETTER than the 360 mm set's 1.1-2.8%. See
            # stack_variance: with only two unclipped tiers between 1.00 and
            # 1.10 R, the n>=3 rule was measuring the dim tail of a third,
            # partly clipped one. The ringing in that band is real and still
            # unexplained; the glare story was an artifact of the diagnostic and
            # is no longer offered as the cause.
            _cv = aq.get("cov_limb")
            if _cv is not None and np.isfinite(_cv) and _cv > 0.30:
                progress.log(
                    "WARNING: the tiers disagree by %.2f (coefficient of "
                    "variation) at the limb, where a well-behaved set sits near "
                    "0.05-0.08. They are aligned -- this is a disagreement in "
                    "BRIGHTNESS, not position, in a rim %.0f px wide just "
                    "outside the limb. Every detail filter normalises against a "
                    "local mean, so a disagreement there is amplified into "
                    "concentric rings: measured at 2-5x the reference set in "
                    "1.02-1.15 R. The cause is not established -- an earlier "
                    "build named veiling glare here and that was withdrawn in "
                    "0.22.17, because the number it rested on came from a "
                    "contaminated sample."
                    % (_cv, aq.get("rim_width_px", float("nan"))), None)
                stats["limb_variance_bad"] = round(float(_cv), 3)
            _nt = aq.get("cov_limb_unmeasurable")
            if _nt is not None:
                progress.log(
                    "the tier-agreement test at the limb is not measurable on "
                    "this bracket: only %.0f tier(s) hold unclipped signal "
                    "between 1.00 and 1.10 R and the coefficient of variation "
                    "needs three. Reported as absent rather than as a number "
                    "built on a partly clipped tier -- that is what produced a "
                    "false 0.79 on the 250 mm test set for three releases. A "
                    "shorter tier at the top of the bracket is what would make "
                    "it measurable." % _nt, None)
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
        rgb = _demosaic(stacks_bayer[s], demosaic_method)
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
    rgb = _demosaic(stacks_bayer[s], demosaic_method)
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
    lstats = detail.build_layers(wd, progress, denoise=denoise,
                                 earthshine=earthshine,
                                 fnrgf_preset=fnrgf_preset)
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
               # bias and darks are inputs on exactly the same footing: they
               # change the stacked data, so adding, replacing or removing one
               # has to re-stack. `applied` and not merely `requested`, for the
               # reason given for the flat above.
               "bias_dir": _bias_dir, "dark_dir": _dark_dir,
               "bias_inputs": _dark.fingerprint(_bias_dir),
               "dark_inputs": _dark.fingerprint(_dark_dir),
               "bias_applied": bias_master is not None,
               "dark_applied": dark_rate is not None,
               "denoise": denoise, "earthshine": bool(earthshine),
               "despeckle": bool(despeckle), "frames": frames,
               # The photometric mode changes the merged data itself, so it
               # belongs in the cache key exactly as the merge weight does.
               # Without it, dropping hill_photometry.txt into the folder and
               # pressing Start returns the previous merge and the setting
               # appears to do nothing -- the same trap the feather setting
               # fell into before 0.22.x.
               # the mode that was ASKED FOR, not the one that succeeded
               # (0.22.76). This recorded photometry_mode, which is set only
               # when _hill_chain returns a usable fit; the server compares it
               # against "is the marker file there". So after a fallback -- the
               # chain logs "using the scalar chain" and carries on -- opts said
               # "scalar" while the server expected "hill", the keys never
               # matched again, and every Start re-stacked the whole folder with
               # nothing saying why. What was APPLIED is recorded separately.
               "photometry": ("hill" if _photometry_requested(folder, photometry)
                            else "scalar"),
               "photometry_applied": stats.get("photometry_mode", "scalar"),
               # Both of these change the MERGED data, so both belong in the
               # cache key for exactly the reason the merge weight does: switch
               # one and reuse the previous merge and the setting appears to do
               # nothing. The feather is checked by the server against a key it
               # never wrote -- see the note there.
               "wb_source": str(wb_source),
               "demosaic": str(demosaic_method),
               "feather": str(feather),
               # 0.23.2's four selectors. Three of them change the MERGED data
               # -- the correlation flavour and the pre-filter change where every
               # tier lands, the tier combine changes what each tier IS -- so
               # they belong in the cache key for the same reason the merge
               # weight does. `fnrgf_preset` only changes a detail layer, and it
               # is here anyway: the layer is cached beside the merge and reusing
               # it after a preset change is the same silent staleness.
               "align_corr": str(align_corr),
               "align_filter": str(align_filter),
               "stack_combine": str(stack_combine),
               "fnrgf_preset": str(fnrgf_preset),
               "build": __version__},
              open(os.path.join(wd, "opts.json"), "w"))
    import datetime
    stats["finished"] = datetime.datetime.now().isoformat(timespec="seconds")
    # A step's cost is the gap to the NEXT log line, so the last step of the run
    # has never had one and has never been timed. On the first full report from
    # a real run the whole Pellett layer was simply absent from the list, and
    # the "steps under 1s" remainder quietly absorbed it. Close the interval
    # before measuring it.
    progress.log("assembling the run report", None)
    _tm = _timing_summary(progress)
    if _tm:
        stats["timing"] = _tm
        progress.log("run time " + _fmt_dur(_tm["total_s"]) + " — slowest: "
                     + "; ".join(f"{m} {_fmt_dur(d)}" for d, m in _tm["slowest"]),
                     None)
    txt = _report.write(wd, stats)
    for line in txt.split("\n"):
        progress.log(line, None)
    # Written every run, before the "complete" line, so a tester reporting a
    # problem already has the one file that makes it diagnosable. See
    # diagnostics.py for what is and is not in it -- the short version is
    # profiles and 512px thumbnails of the detail layers, no raw data.
    try:
        from . import diagnostics as _diag
        _diag.write_bundle(wd, folder, stats, progress)
    except Exception:
        pass
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
    from .raw import open_frame, demosaic as _demosaic, pick_wb
    wd = workdir(folder)
    geo = json.load(open(os.path.join(wd, "geometry.json")))
    progress.log(f"decoding {os.path.basename(raw_path)}...", 0.1)
    rf = open_frame(raw_path)
    # the contact frame goes through the same optics, so it gets the same flat
    # -- but only if the run that built the composite actually applied one
    _mf = os.path.join(wd, "masterflat.npy")
    _op = os.path.join(wd, "opts.json")
    _used_flat = False
    # The contact frame is screen-blended onto the composite, so it has to be
    # decoded the SAME way the composite was -- same white balance, same
    # demosaic -- or the ring arrives in a different colour from the corona it
    # sits on. Read from the run's own opts, not from today's defaults.
    _o = {}
    try:
        _o = json.load(open(_op))
        _used_flat = bool(_o.get("flat_applied", False))
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
    rgb = _demosaic(rf.bayer, _o.get("demosaic", "mhc"))
    _cwb, _ = pick_wb(rf, _o.get("wb_source", "camera"))
    rgb *= _cwb[None, None, :]
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
    # THE LIMB FIT CANNOT BE TRUSTED ON A CONTACT FRAME, SO IT IS CHECKED.
    #
    # A totality frame is a dark disc inside a corona, which is what fit_limb
    # was written for. A 2nd/3rd-contact frame is not that. 99.4% of the reference set's
    # P1072722 is below the noise -- the sky is dark, the corona is short-
    # exposed to nothing, and the only structure in the frame is the blazing
    # crescent. A limb finder given that fits the CRESCENT'S arc, because it is
    # the only edge there is.
    #
    # Measured on his frame. The crescent sits 588 px from the composite's
    # lunar centre, and the limb is at R=619 -- so as shot, with no
    # registration at all, it is already within ~30 px of correct, which is
    # what you expect from a tracked sequence. The fit put the lunar centre at
    # (2915, 4248), essentially on top of the crescent, and the resulting shift
    # moved the layer 399 px UP AND LEFT: the crescent ended up 190 px from the
    # centre, well inside the disc where no crescent can be.
    #
    # So: attempt the fit, then check what it is asking for. A tracked sequence
    # needs a small correction or none. A fit demanding a large move has found
    # the wrong thing, and doing nothing is much closer to right than obeying
    # it. This replaces the 0.20.2 reasoning, which blamed a scale mismatch on
    # a different focal length -- the EXIF says 600 mm for both, and the ratio
    # that story rested on came from circle-fitting a crescent, which is not a
    # circle. Wrong diagnosis, wrong fix.
    _rat = geo["R"] / R if R > 1 else 1.0
    _mv = float(np.hypot(dy, dx))
    _bad = (_mv > 0.25 * geo["R"]) or not (0.8 < _rat < 1.25)
    if _bad:
        progress.log(f"contact frame: the limb fit wants to move it {_mv:.0f}px and "
                     f"scale it x{_rat:.3f} (it fitted R={R:.0f}px against the "
                     f"composite's {geo['R']:.0f}px). That is far more than a tracked "
                     f"sequence needs, so the fit found the bead rather than the Moon "
                     f"and is being ignored — the frame is overlaid AS SHOT. Use the "
                     f"ring offset and size sliders if it needs nudging.", None)
        dy, dx, sc = 0.0, 0.0, 1.0
    else:
        sc = _rat
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
    # DID THE ALIGNMENT ACTUALLY WORK? A crescent of photosphere is always
    # OUTSIDE the lunar limb -- the Moon is what is hiding the rest of it. So if
    # the aligned frame's bright arc lands inside the composite's disc, the
    # registration failed, whatever the fit residuals said. This is cheap, it is
    # geometry rather than taste, and it is what was missing when the reported ring
    # came out 327px off centre with nothing in the log to say so.
    try:
        _al = 0.2126 * rgb[:, :, 0] + 0.7152 * rgb[:, :, 1] + 0.0722 * rgb[:, :, 2]
        _s = _al[::4, ::4]
        _m = _s >= np.percentile(_s, 99.98)
        _ys, _xs = np.nonzero(_m)
        if _ys.size > 200:
            _Y, _X = _ys * 4.0, _xs * 4.0
            _rr = np.hypot(_Y - geo["cy"], _X - geo["cx"])
            _med = float(np.median(_rr))
            # algebraic circle fit, valid however far the centre has drifted
            _A = np.stack([_X, _Y, np.ones_like(_X)], axis=1)
            _sol, *_ = np.linalg.lstsq(_A, _X * _X + _Y * _Y, rcond=None)
            _fx, _fy = _sol[0] / 2, _sol[1] / 2
            _fr = float(np.sqrt(max(_sol[2] + _fx * _fx + _fy * _fy, 0.0)))
            progress.log(f"contact frame: its bright arc fits a circle of R={_fr:.0f}px "
                         f"centred {np.hypot(_fy - geo['cy'], _fx - geo['cx']):.0f}px from "
                         f"the composite disc (limb R={geo['R']:.0f}px); arc sits at "
                         f"{_med:.0f}px from that centre", None)
            if _med < 0.90 * geo["R"]:
                progress.log("contact frame: THE RING LANDS INSIDE THE LUNAR DISC. A "
                             "crescent of photosphere cannot be there, so this frame is "
                             "not registered to the composite. The ring sliders cannot "
                             "correct an error this size — the limb fit on the contact "
                             "frame is what needs looking at.", None)
        del _al, _s, _m
    except Exception as _e:
        progress.log(f"contact frame: could not check the ring geometry ({_e})", None)
    top = np.percentile(lum, 99.9)
    disp = np.clip(rgb / max(top, 1e-6), 0, 1) ** (1 / 2.2)
    # match the composite frame size if the sensor crop differs slightly
    Hc, Wc = load_big(os.path.join(wd, "hdr_lum.npy")).shape
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
