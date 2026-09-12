"""Bias and dark calibration.

A bias frame is what the sensor reads with no light and no time: the offset the
electronics put under every pixel, plus read noise.  A dark frame is the same
thing after t seconds, so it also carries the charge that accumulated with no
light -- dark current, which is per-pixel, roughly linear in time, and roughly
doubles every 5-7 C.  Both are additive, so both come off before anything
multiplicative:

    corrected = (light - bias - rate * t) / F

`rate` is per pixel, in ADU per second, and it is measured here as
(master_dark - master_bias) / t_dark.  A bias frame is just the t = 0 case of
the same measurement, which is why one module does both.

WHAT THIS IS AND IS NOT FOR, ON THIS APP
----------------------------------------
On a camera raw, most of the bias is already gone before this module sees
anything.  `raw.py` subtracts `black_level_per_channel` per Bayer position at
load, and the per-channel residue left after that, measured on three real raw
brackets, was 0.14, 0.26 and 2.14 ADU.  So a master bias built from camera raws
is not "the black level" -- the metadata already gave that -- it is the RESIDUE
the metadata got wrong, plus whatever fixed pattern the electronics leave.  It
can legitimately come out near zero, and if it does, that is a result and this
module says so rather than hiding it.

Where the frames earn their keep is the part the metadata cannot describe:
per-pixel structure.  Warm and hot pixels, column or row offsets, and any
fixed spatial pattern that grows with exposure time.  Those live where the
outer corona lives -- tens of ADU above black on a long tier -- and they are
identical in every frame, so stacking the lights cannot average them out.

THE NOISE PROBLEM, AND WHY NOTHING HERE IS SMOOTHED
---------------------------------------------------
Subtracting a master injects that master's noise into every light frame,
identically, so stacking cannot remove it either.  `flat.py` answers this by
smoothing until the injected noise is under a target.  THIS MODULE MUST NOT DO
THAT.  The entire value of a master dark is its per-pixel detail; a Gaussian
over it destroys exactly the hot pixels and column structure it exists to
remove, and leaves only the noise.

So instead of smoothing, this module MEASURES both sides and reports them:

    removes   the fixed pattern actually present in the master, robust sd
    injects   the master's own noise, from two independent half-stacks

If `removes` is not comfortably larger than `injects`, the correction is making
the picture worse, and the log says so in those words.  The half-stack
difference is not a proxy for the master's noise -- var((A-B)/2) is exactly
var(master) -- which is the same argument flat.py rests on.

ISO AND EXPOSURE HAVE TO MATCH
------------------------------
Black level and dark current both move with ISO, and dark current moves with
temperature and time.  A dark shot at a different ISO than the lights is not a
dark for those lights, and this module refuses it rather than subtracting a
number from the wrong sensor state.  Exposure need not match -- that is what
the per-second rate is for -- but it is checked and reported, because
extrapolating a 1.6 s measurement onto a 1/2000 s tier is arithmetically fine
and physically pointless.
"""
from __future__ import annotations
import os, json
import numpy as np

#: A correction is worth applying when what it removes is at least this many
#: times what it injects. Below it the master is mostly its own noise.
WORTH_IT = 2.0

#: Sigma for counting a photosite as WARM in the master dark. Matches the
#: light-frame path in raw._outlier_flags so the two counts are on the same
#: scale -- which is what made it obvious they are not measuring the same
#: population at all. Nothing is repaired from this; see the note in build().
DARK_GATE = True      # subtract a photosite's dark rate only where the master
                      # actually measured something there. On a typical set 99%+
                      # of the master is its own noise rather than dark current,
                      # and subtracting that injects noise while removing
                      # nothing -- which is why a 20-frame master can report
                      # that it removes only 1.1x the noise it adds.
DARK_GATE_KS = (0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0)
DEFECT_K = 6.0

_BIAS_DIR_NAMES = ("bias", "biases", "offset", "offsets", "masterbias",
                   "bias_frames", "bias-frames")
_DARK_DIR_NAMES = ("darks", "dark", "masterdark", "dark_frames",
                   "dark-frames")


def _find_dir(folder, names):
    """A conventionally-named subfolder of `folder`, or None.

    Case-insensitive, one level down. `list_raws()` reads only the files
    directly in a folder, so these can never be mistaken for light frames --
    the same property flats rely on.
    """
    try:
        entries = sorted(os.listdir(folder))
    except OSError:
        return None
    for name in entries:
        p = os.path.join(folder, name)
        if os.path.isdir(p) and name.lower() in names:
            return p
    return None


def find_bias_dir(folder):
    return _find_dir(folder, _BIAS_DIR_NAMES)


def find_dark_dir(folder):
    return _find_dir(folder, _DARK_DIR_NAMES)


def fingerprint(d):
    """Identity of the calibration files, for cache validation."""
    from .raw import list_raws
    out = []
    if not d or not os.path.isdir(d):
        return out
    for p in sorted(list_raws(d)):
        try:
            st = os.stat(p)
            out.append([os.path.basename(p), int(st.st_size), int(st.st_mtime)])
        except OSError:
            out.append([os.path.basename(p), -1, -1])
    return out


def _removes_rms(pattern, per_px_noise):
    """RMS of the REAL fixed pattern in `pattern`, with the master's own noise
    taken back out.

    WHY NOT A ROBUST SD, which is what this used until 0.23.2 and which was
    wrong in both directions.

    Too LOW: once the rate map is gated, 99%+ of it is exactly zero, and a MAD
    of a mostly-zero array is 0. That is what turned a 1.1x verdict into "0.0x"
    on the first gated run -- the correction was fine, the measurement was not.

    Too HIGH, and this one was there all along: on a sensor where only 0.3% of
    photosites carry dark current, the robust spread of the UNGATED map is the
    master's own noise, not the pattern. So the old number compared the noise
    against itself and printed a ratio near 1, which is exactly the "removes
    only 1.1x its own noise" that has been in these reports for months. It was
    not measuring a marginal correction; it was measuring nothing.

    E[pattern^2] = truth^2 + noise^2 in the mean, so subtract the noise power
    once, over the whole frame, and clip at zero. Debiased on the MEAN, never
    per pixel -- E[max(z^2 - 1, 0)] = 0.52 for a standard normal, so per-pixel
    clipping finds signal in pure noise.
    """
    p2 = float(np.mean(np.square(np.asarray(pattern, np.float64))))
    return float(np.sqrt(max(p2 - float(per_px_noise) ** 2, 0.0)))


def _robust_sd(d):
    m = float(np.median(d))
    return float(1.4826 * np.median(np.abs(d - m)))


def _stack(paths, shape, log, label, want_iso=None):
    """Per-pixel master plus two independent half-stacks.

    Frames arrive from `open_frame` with the reported black already subtracted,
    so a master built here is in the SAME units as a decoded light frame and
    subtracts from one directly. That is the whole reason this does not read
    the files itself.
    """
    from .raw import list_raws, open_frame, read_exif
    info = {"n_found": len(paths)}
    sumA = sumB = mn = mx = None
    nA = nB = 0
    used, rejected = [], []
    secs, isos = [], []
    sat0 = None
    for i, p in enumerate(paths):
        name = os.path.basename(p)
        try:
            sec, iso, _ts = read_exif(p)
        except Exception:
            sec, iso = None, None
        try:
            rf = open_frame(p)
            b = np.asarray(rf.bayer, np.float32)
            sat = float(rf.sat_level)
            del rf
        except Exception as e:
            rejected.append((name, f"could not be decoded ({e})"))
            continue
        if shape is not None and tuple(b.shape) != tuple(shape):
            rejected.append((name, f"is {b.shape[1]}x{b.shape[0]}, the light "
                                   f"frames are {shape[1]}x{shape[0]}"))
            del b
            continue
        if want_iso and iso and int(iso) != int(want_iso):
            rejected.append((name, f"was shot at ISO {iso}, the light frames "
                                   f"at ISO {want_iso} — black level and dark "
                                   f"current both move with ISO"))
            del b
            continue
        # A light frame filed in the wrong folder is the failure mode that
        # would do the most damage silently: it would be subtracted from every
        # tier as if it were an offset. A calibration frame is dark by
        # definition, so anything with real signal in it is not one.
        hi = float(np.percentile(b, 99.9))
        if hi > 0.25 * sat:
            rejected.append((name, f"is not a dark frame — its 99.9th "
                                   f"percentile is {100 * hi / sat:.0f}% of "
                                   f"saturation"))
            del b
            continue
        if sat0 is None:
            sat0 = sat
        if sec is not None:
            secs.append(float(sec))
        if iso:
            isos.append(int(iso))
        if sumA is None:
            sumA = np.zeros_like(b)
            sumB = np.zeros_like(b)
            mn = b.copy()
            mx = b.copy()
        else:
            np.minimum(mn, b, out=mn)
            np.maximum(mx, b, out=mx)
        # Parity split, as in flat.py: two half-stacks that saw the same drift,
        # so their difference is noise and not a trend.
        if len(used) % 2 == 0:
            sumA += b
            nA += 1
        else:
            sumB += b
            nB += 1
        used.append(name)
        log(f"{label} {i + 1}/{len(paths)}: {name}", None)
        del b

    for name, why in rejected:
        log(f"{label} rejected: {name} {why}", None)
    info["n_used"] = len(used)
    info["rejected"] = [{"file": n, "why": w} for n, w in rejected]
    n = len(used)
    if n < 2:
        info["error"] = (f"only {n} usable {label} frame(s) — at least 2 are "
                         f"needed")
        return None, None, None, info

    master = sumA + sumB
    if n >= 5:
        mn += mx
        master -= mn
        master /= np.float32(n - 2)
        info["combine"] = f"min/max-trimmed mean of {n} frames"
        trim = float(np.sqrt(n / float(n - 2)))
    else:
        master /= np.float32(n)
        info["combine"] = f"mean of {n} frames"
        trim = 1.0
    del mn, mx
    meanA = (sumA / np.float32(max(nA, 1))) if nA else None
    meanB = (sumB / np.float32(max(nB, 1))) if nB else None
    del sumA, sumB
    if meanA is not None and meanB is not None:
        info["noise"] = float(trim * _robust_sd((meanA - meanB) * np.float32(0.5)))
        # ... and per Bayer offset, because one global spread is a blend of
        # four different sensors. A camera whose blue read noise is twice its
        # green (3.1/2.4/2.4/4.6 ADU on the S1R II fixture) gets a threshold
        # that is too high in green and too low in blue, and the too-low half
        # is the one that edits good data.
        info["noise_cfa"] = [
            float(trim * _robust_sd(
                ((meanA - meanB) * np.float32(0.5))[oy::2, ox::2]))
            for oy in (0, 1) for ox in (0, 1)]
        # PER-FRAME read noise, which is a different number from the master's
        # noise above and the one the rest of the pipeline actually wants.
        #
        # The two half-stack means differ only by noise, so
        #     sd(meanA - meanB) = sigma * sqrt(1/nA + 1/nB)
        # and sigma follows exactly. No trim factor enters: trim describes what
        # min/max rejection did to the MASTER, and these means are untrimmed.
        #
        # Per Bayer offset, not one number. After the per-channel analogue gain
        # a camera applies before the ADC, read noise in ADU is not the same in
        # R, G and B, and the noise model that uses it works per offset anyway.
        _f = float(np.sqrt(1.0 / max(nA, 1) + 1.0 / max(nB, 1)))
        if _f > 0:
            _d = meanA - meanB
            info["read_noise"] = float(_robust_sd(_d) / _f)
            info["read_noise_cfa"] = [
                float(_robust_sd(_d[oy::2, ox::2]) / _f)
                for oy in (0, 1) for ox in (0, 1)]
            del _d
    if secs:
        info["seconds"] = float(np.median(secs))
        info["seconds_spread"] = [float(min(secs)), float(max(secs))]
        if max(secs) > 1.001 * min(secs):
            log(f"{label}: exposures are not all the same "
                f"({min(secs):g}s to {max(secs):g}s) — the rate is scaled from "
                f"the median, {np.median(secs):g}s", None)
    if isos:
        info["iso"] = int(np.median(isos))
    return master.astype(np.float32), meanA, meanB, info


def build(folder, bias_dir, dark_dir, shape=None, progress=None,
          light_iso=None, max_light_seconds=None):
    """Master bias and per-second dark rate for a light folder.

    Returns (bias, rate, info). Either array may be None. `bias` is in decoded
    light-frame units and subtracts directly; `rate` is ADU per second per
    pixel and is multiplied by each tier's exposure before subtracting.
    """
    from .raw import list_raws
    log = (progress.log if progress is not None else (lambda *a, **k: None))
    info = {"bias_dir": bias_dir, "dark_dir": dark_dir}
    bias = rate = None

    if bias_dir and os.path.isdir(bias_dir):
        b, bA, bB, bi = _stack(list_raws(bias_dir), shape, log, "bias",
                               light_iso)
        info["bias"] = bi
        if b is not None:
            bi["removes"] = _robust_sd(b)
            bi["level"] = float(np.median(b))
            bias = b
        del bA, bB

    if dark_dir and os.path.isdir(dark_dir):
        d, dA, dB, di = _stack(list_raws(dark_dir), shape, log, "dark",
                               light_iso)
        info["dark"] = di
        t = di.get("seconds")
        if d is None:
            pass
        elif not t or t <= 0:
            di["error"] = ("no usable exposure time in the dark frames — "
                           "without it the dark cannot be scaled to a tier")
        else:
            # The rate is (dark - bias) / t. Without a bias set the dark's own
            # offset would be subtracted from every tier scaled by exposure,
            # which is wrong for every tier except the dark's own -- so say so
            # and use the dark's median as a stand-in offset rather than
            # silently smearing it across the bracket.
            if bias is not None:
                sig = d - bias
                di["bias_source"] = "master bias"
            else:
                sig = d - np.float32(np.median(d))
                di["bias_source"] = ("the dark's own median (no bias frames) — "
                                     "only the per-pixel structure is used, "
                                     "not the offset")
            rate = (sig / np.float32(t)).astype(np.float32)
            di["rate_median"] = float(np.median(rate))
            # The mean as well as the median, and not for decoration: the
            # master carries per-pixel noise, and a median taken across pixels
            # of (skewed truth + symmetric noise) is pulled towards the MEAN of
            # the truth, not its median. Dark current is strongly right-skewed
            # -- a long tail of warm pixels -- so the two differ by tens of
            # percent and only the mean is unbiased under that noise. Anything
            # comparing this against a known rate has to use the mean.
            di["rate_mean"] = float(np.mean(rate))
            if "noise" in di:
                _bn = (info.get("bias", {}) or {}).get("noise", 0.0) or 0.0
                di["injects"] = float(np.hypot(di["noise"], _bn))
                di["removes"] = _removes_rms(sig, di["injects"])
            else:
                di["removes"] = _robust_sd(sig)
            di["removes_raw_sd"] = _robust_sd(sig)
            di["seconds"] = float(t)
            # THE DEFECT MAP, from the thing defects are defined on.
            #
            # Without darks the pipeline finds hot photosites on the shortest
            # LIGHT tier, as pixels far above their same-colour neighbours. It
            # works, and it is inference: a hot pixel and a star are the same
            # observation there, which is why that path has to use the darkest
            # tier and still cannot separate the two in principle. A dark frame
            # measures the defect directly -- a photosite's own dark current,
            # with no sky in the frame at all -- and that is what the rate map
            # already is.
            #
            # Detected against the same-colour neighbour median, like the light
            # path, but with a KNOWN sigma instead of a fitted one: the master's
            # own noise divided by the exposure, since dark current at these
            # levels contributes no shot noise worth modelling. Per Bayer offset
            # because neighbours of a different colour are a different sensor.
            # `injects`, NOT `noise`. The rate map is (dark master - bias
            # master) / t, so its per-pixel noise is both masters' added in
            # quadrature -- which is exactly what `injects` already is. Using
            # the dark's alone understates it by the bias's share: on 20 darks
            # over 50 bias frames that is a factor 0.845, an 18% low threshold,
            # and it put 62 false positives on a synthetic sensor whose true
            # count is 0. The statistic being thresholded is a pixel minus the
            # median of its eight same-colour neighbours, whose robust spread
            # is 0.86 of the per-pixel sigma and whose largest excursion over
            # 1.9 million clean pixels is 5.6 -- so 6 sigma on the correct
            # scale flags nothing that is not there.
            # A CENSUS, NOT A DEFECT MAP. See the note in pipeline.py where
            # this is consumed: a warm photosite is corrected by subtracting
            # its own rate, and the defects a light frame finds at 1/500 s are
            # a different population entirely -- dark current contributes a
            # quarter of an ADU at that shutter speed. This counts what is warm
            # so the run can say so; nothing is repaired from it.
            _dn = di.get("noise_cfa")
            _bn4 = ((info.get("bias") or {}).get("noise_cfa")
                    or [0.0, 0.0, 0.0, 0.0])
            if _dn and len(_dn) == 4:
                from scipy import ndimage
                _hot = np.zeros(rate.shape, bool)
                _sigs = []
                for _i, (_oy, _ox) in enumerate(((0, 0), (0, 1), (1, 0), (1, 1))):
                    _sig = float(np.hypot(_dn[_i], _bn4[_i])) / float(t)
                    _sigs.append(_sig)
                    _s = rate[_oy::2, _ox::2]
                    _m = ndimage.median_filter(_s, size=3, mode="nearest")
                    _hot[_oy::2, _ox::2] = (_s - _m) > DEFECT_K * _sig
                di["defects"] = int(_hot.sum())
                di["defect_sigma_rate"] = _sigs
                di["defect_k"] = DEFECT_K
                info["defect_map"] = _hot
                del _hot
            # SUBTRACT ONLY WHERE THERE IS SOMETHING TO SUBTRACT.
            #
            # The verdict below compares what the master removes against the
            # noise it injects, and on the reference set's 20 darks that ratio came out 1.1x
            # -- barely worth applying. The census above says why: 0.314% of
            # photosites carry measurable dark current. For the other 99.7% the
            # master's value is not dark current, it is the master's own noise,
            # and subtracting it adds that noise to every frame of every tier
            # while removing nothing.
            #
            # So gate the rate map. A photosite keeps its measured rate where
            # that rate stands above the master's own per-pixel noise, and is
            # set to exactly zero where it does not. This is not a threshold on
            # the DEFECT scale -- that one is 6 sigma against a neighbour median
            # and answers "is this photosite broken". This one is 2 sigma
            # against zero and answers "did we measure anything here at all",
            # which is a much weaker question and the right one: at 2 sigma the
            # pixels kept are those where subtracting helps more than it hurts,
            # and the few real-but-faint ones that fall below it were
            # contributing less than the noise they came with.
            #
            # The gain is arithmetic and large. Noise injected scales as the
            # square root of the FRACTION of pixels left un-gated, so dropping
            # 99% of them takes the injected noise to a tenth of what it was
            # while the fixed pattern removed is almost untouched, because
            # almost all of it lives in the pixels that were kept.
            if _dn and len(_dn) == 4 and DARK_GATE:
                # WHICH THRESHOLD, CHOSEN FROM THE DATA. There is no constant
                # that is right for two different sensors: the trade is dark
                # current left un-subtracted below the cut against master noise
                # injected above it, and which dominates depends entirely on how
                # the dark current is distributed. Measured on four plausible
                # sensors, residual RMS in a corrected frame, best threshold:
                #
                #   sparse + very warm (0.3% warm, like ours)   6 sigma
                #   sparse + mildly warm, tail near the noise   3 sigma
                #   uniform low-level dark current              1.5 sigma
                #   no dark current at all                      any high cut
                #
                # Both terms are estimable from the master alone. Above the cut
                # the error is the master's own noise, sigma^2 per kept
                # photosite. Below it the error is the dark current we declined
                # to subtract -- and that is NOT the measured rate^2 of the
                # zeroed pixels, which carries the same noise, but
                # mean(rate^2) - sigma^2 over the group. Debiased ONCE on the
                # group mean, never per pixel: E[max(z^2 - 1, 0)] = 0.52 for a
                # standard normal, so per-pixel clipping scores 0.52 sigma^2 on
                # a photosite with no dark current at all and the picker then
                # sees signal everywhere. That was the first version's bug and
                # it chose the same threshold for every sensor.
                #
                #     predicted MSE(k) = f * sigma^2
                #                      + (1 - f) * max(mean(rate_zeroed^2) - sigma^2, 0)
                #
                # On the four above this lands on the optimum for two of them
                # and within 3% for the others.
                _sg4 = np.zeros(rate.shape, np.float32)
                for _i, (_oy, _ox) in enumerate(((0, 0), (0, 1), (1, 0), (1, 1))):
                    _sg4[_oy::2, _ox::2] = float(
                        np.hypot(_dn[_i], _bn4[_i])) / float(t)
                _z = np.abs(rate) / np.maximum(_sg4, 1e-9)
                _r2 = (rate.astype(np.float64) ** 2)
                _s2 = (_sg4.astype(np.float64) ** 2)
                _best_k, _best_m = 0.0, None
                for _k in DARK_GATE_KS:
                    _keep = _z > _k if _k > 0 else np.ones(rate.shape, bool)
                    _f = float(_keep.mean())
                    if _f >= 1.0:
                        _m = float(np.mean(_s2))
                    else:
                        _lost = max(float(np.mean(_r2[~_keep]))
                                    - float(np.mean(_s2[~_keep])), 0.0)
                        _m = _f * float(np.mean(_s2)) + (1.0 - _f) * _lost
                    if _best_m is None or _m < _best_m:
                        _best_k, _best_m = _k, _m
                    del _keep
                del _z, _r2, _s2
                _keep = (np.abs(rate) > _best_k * _sg4) if _best_k > 0 \
                    else np.ones(rate.shape, bool)
                _frac = float(_keep.mean())
                _rm_before = float(di["removes"])
                rate = np.where(_keep, rate, np.float32(0.0)).astype(np.float32)
                di["gate_k"] = float(_best_k)
                di["gate_kept_frac"] = _frac
                if "injects" in di:
                    di["injects_ungated"] = float(di["injects"])
                    # the zeroed photosites now carry exactly no noise, so what
                    # the master injects scales with the square root of the
                    # fraction still being subtracted
                    di["injects"] = float(di["injects"]
                                          * np.sqrt(max(_frac, 1e-9)))
                    di["removes"] = _removes_rms(rate * np.float32(t),
                                                 di["injects"])
                else:
                    di["removes"] = _robust_sd(rate * np.float32(t))
                di["removes_ungated"] = _rm_before
                del _keep, _sg4
            # What this actually costs and buys at the longest tier, which is
            # where a dark matters and where it is being extrapolated least.
            if max_light_seconds:
                di["at_longest_tier"] = {
                    "seconds": float(max_light_seconds),
                    "removes": float(di["removes"] * max_light_seconds / t),
                    "injects": float(di.get("injects", 0.0)
                                     * max_light_seconds / t)}
        del dA, dB

    return bias, rate, info


def describe(info):
    """One-line-per-fact summary, for the run log and the report."""
    out = []
    bi = info.get("bias") or {}
    di = info.get("dark") or {}
    if bi:
        if bi.get("error"):
            out.append(f"bias: {bi['error']}")
        else:
            out.append(f"bias: {bi.get('combine', '')} from "
                       f"{os.path.basename(info.get('bias_dir') or '')}/"
                       + (f", ISO {bi['iso']}" if bi.get("iso") else ""))
            if bi.get("read_noise"):
                _c = bi.get("read_noise_cfa") or []
                out.append(
                    f"bias: read noise {bi['read_noise']:.2f} ADU per frame"
                    + (" (R %.2f, G %.2f/%.2f, B %.2f)" % tuple(_c)
                       if len(_c) == 4 else ""))
            rem, inj = bi.get("removes", 0.0), bi.get("noise", 0.0)
            out.append(f"bias: offset {bi.get('level', 0.0):+.3f} ADU on top of "
                       f"the black level the raw file already reported; fixed "
                       f"pattern {rem:.3f} ADU, master noise {inj:.3f} ADU")
            out.append(_verdict("bias", rem, inj))
    if di:
        if di.get("error"):
            out.append(f"dark: {di['error']}")
        elif di.get("rate_median") is not None:
            out.append(f"dark: {di.get('combine', '')} from "
                       f"{os.path.basename(info.get('dark_dir') or '')}/ at "
                       f"{di.get('seconds', 0):g}s"
                       + (f", ISO {di['iso']}" if di.get("iso") else ""))
            out.append(f"dark: offset taken from {di.get('bias_source', '?')}")
            out.append(f"dark: median dark current "
                       f"{di.get('rate_median', 0.0):.4f} ADU/s")
            if di.get("defects") is not None:
                _sg = di.get("defect_sigma_rate") or []
                out.append(
                    f"dark: {di['defects']} photosites carry dark current "
                    f"{di.get('defect_k', DEFECT_K):.0f} sigma above their "
                    f"same-colour neighbours"
                    + (" (sigma %.4f-%.4f ADU/s per channel)"
                       % (min(_sg), max(_sg)) if _sg else "")
                    + " — a census, not defects: each one's rate is subtracted "
                      "individually, which is the correct fix for it")
            if di.get("gate_kept_frac") is not None:
                out.append(
                    f"dark: subtracting the rate on {100.0 * di['gate_kept_frac']:.2f}% "
                    f"of photosites ({di['gate_k']:.1f} sigma above zero, chosen "
                    f"by predicted residual) — on the rest the master measured "
                    f"only its own noise, and subtracting it would add that "
                    f"noise while removing nothing")
            at = di.get("at_longest_tier")
            if at:
                out.append(f"dark: at the longest tier ({at['seconds']:g}s) it "
                           f"removes {at['removes']:.3f} ADU of fixed pattern "
                           f"and injects {at['injects']:.3f} ADU of noise")
                out.append(_verdict("dark", at["removes"], at["injects"]))
            else:
                out.append(_verdict("dark", di.get("removes", 0.0),
                                    di.get("injects", 0.0)))
    return out


def _verdict(what, removes, injects):
    if injects <= 0:
        return (f"{what}: noise of the master could not be measured — treat "
                f"the correction as unverified")
    r = removes / injects
    if r >= WORTH_IT:
        return (f"{what}: removes {r:.1f}x more fixed pattern than it injects "
                f"noise — worth applying")
    return (f"{what}: WARNING — removes only {r:.1f}x its own noise. Below "
            f"{WORTH_IT:.0f}x the master is mostly noise and subtracting it "
            f"puts more in than it takes out. More frames would fix this; "
            f"nothing else will.")


def load_or_build(folder, bias_dir, dark_dir, shape, progress=None,
                  workdir=None, light_iso=None, max_light_seconds=None):
    """Cached masters for a light folder. Same contract as flat.load_or_build."""
    log = (progress.log if progress is not None else (lambda *a, **k: None))
    key = {"bias": fingerprint(bias_dir), "dark": fingerprint(dark_dir),
           "shape": list(shape) if shape is not None else None,
           "iso": int(light_iso or 0),
           # fmt 2 adds the per-offset noise, the read noise and the defect
           # map. A cache written by fmt 1 has none of them, and silently
           # reusing it would leave the defect map empty and the read noise
           # unmeasured on exactly the folders that already have calibration
           # frames -- so it is a different key, not a migration.
           "tmax": round(float(max_light_seconds or 0.0), 6), "fmt": 2}
    npy = os.path.join(workdir, "masterdark.npz") if workdir else None
    js = os.path.join(workdir, "masterdark.json") if workdir else None
    if npy and os.path.exists(npy) and js and os.path.exists(js):
        try:
            cached = json.load(open(js))
            if cached.get("key") == key:
                z = np.load(npy)
                bias = z["bias"] if "bias" in z.files else None
                rate = z["rate"] if "rate" in z.files else None
                if shape is None or all(
                        a is None or tuple(a.shape) == tuple(shape)
                        for a in (bias, rate)):
                    info = cached.get("info", {})
                    info["cached"] = True
                    _ds = info.get("defect_shape")
                    if _ds and "defects" in z.files:
                        _n = int(np.prod(_ds))
                        info["defect_map"] = np.unpackbits(
                            z["defects"])[:_n].astype(bool).reshape(_ds)
                    for line in describe(info):
                        log(line + " (cached)", None)
                    return bias, rate, info
        except Exception:
            pass
    if bias_dir or dark_dir:
        log("building the master bias/dark...", None)
    bias, rate, info = build(folder, bias_dir, dark_dir, shape, progress,
                             light_iso, max_light_seconds)
    for line in describe(info):
        log(line, None)
    if npy and (bias is not None or rate is not None):
        try:
            arrs = {}
            if bias is not None:
                arrs["bias"] = bias
            if rate is not None:
                arrs["rate"] = rate
            # The defect map is an array and rides in the npz; everything else
            # in `info` is JSON. Popping it keeps json.dump from choking on a
            # 44 megapixel bool array, which is how this was found.
            _dm = info.pop("defect_map", None)
            if _dm is not None:
                arrs["defects"] = np.packbits(_dm)
                info["defect_shape"] = list(_dm.shape)
            np.savez_compressed(npy, **arrs)
            json.dump({"key": key, "info": info}, open(js, "w"), indent=1)
            if _dm is not None:
                info["defect_map"] = _dm
        except Exception as e:
            log(f"dark: could not cache the masters ({e})", None)
    return bias, rate, info
