"""TIERS BY MEASURED BRIGHTNESS (2026-10-13).

Grouping by the header's exposure time is the one place the stack trusts the
metadata instead of measuring past it: two frames are averaged as one
measurement only because their headers say they belong together. Usually they
do. On a FITS test set the "exposure times" are labels typed in Siril
("/ Manual group assignment"): the 23 groups step by exactly 1.5x while the
pixels say 1.19x at the short end and 2.0x at the long end, and the group
labelled 2.956 s holds frames at 0.85x, 0.98x, 1.0x and 1.5x of each other.

This module measures every frame's brightness from the pixels and groups
frames that agree. It is a Tiers option ("By measured brightness"), never the
default: on a bracket with honest headers it gives the same tiers, after an
extra read of every frame.

HOW
  * each frame is read once at a quarter of its size (green only);
  * frames are taken in header order (then time, then name), so neighbours
    are close in brightness, and each frame is compared with the one before:
    the slope k of y = k x + c over the pixels both frames hold above 1 % and
    below 70 % of saturation. The fit is made on BIN MEDIANS, with the bins
    taken on the smoothed sum of the two frames. A plain least-squares fit on
    pixels is biased low when x is noisy (regression dilution): on a FITS
    test set's shortest frames it printed a false 5-15 % fall inside every group;
  * the slopes chain into one brightness per frame; sorted, the frames split
    into tiers wherever two neighbours differ by more than GAP (8 %). Real
    steps are 19 % and more (1/4 EV); frame-to-frame scatter inside a tier
    was 0.3-8 % on that set;
  * camera files (exposure from the camera) are only ever SPLIT: a header
    group falls apart where its frames differ by more than GAP, but frames
    with different camera exposure times are never merged (FITS labels are);
  * each pair is placed on each other before it is compared (_align);
  * a tier's key is its median brightness times one scale, chosen so the
    tiers' keys match their frames' header times in the median. The keys stay
    seconds-like, so everything downstream that reads `s` and `s * cal[s]`
    is unchanged; the header times stay on the frames (dark current is scaled
    per frame from them).
"""
import os

import numpy as np

GAP = 1.08          # split a sorted brightness run where neighbours differ by more
LO, HI = 0.01, 0.70  # usable pixel range, fraction of saturation
MIN_PX = 2000


def _green(path):
    """Green at a quarter of the frame size, and the saturation level."""
    from .simple import _planes_direct
    pl = _planes_direct(path)
    if pl is not None:
        a, sat = pl
        g = a[1, 1::4, 1::4]
    else:
        from .raw import open_frame
        rf = open_frame(path)
        g = rf.bayer[0::2, 1::2][::2, ::2]
        sat = rf.sat_level
    return np.ascontiguousarray(g, np.float32), float(sat)


def _box(a, r=2):
    c = np.cumsum(np.cumsum(np.pad(a, ((r + 1, r), (r + 1, r))), 0), 1)
    n = 2 * r + 1
    return (c[n:, n:] - c[:-n, n:] - c[n:, :-n] + c[:-n, :-n]) / (n * n)


#: the frames of a pair are put on top of each other before they are
#: compared (see _align) when the placed pair correlates at least this well;
#: otherwise they are compared where they are
ALIGN_MIN_CORR = 0.5


def _align(x, y, sat):
    """Shift (dy, dx) that puts y on x, and the correlation of the placed pair.

    WHY (2026-11-14). The slope below compares the two frames pixel by
    pixel, so it assumes they show the sky in the same place. Frames shot in
    one burst do; two brackets a minute apart need not. On a test set of two
    brackets whose framing moved about 630 px between them, the unplaced
    pairs gave ratios 6-7 % instead of the 20 % the exposures differ by, three
    pairs of frames from different brackets became one tier each, and each
    such tier averaged two framings into a doubled Moon and corona. So the
    pair is placed first, on the log brightness high-passed (exposure-free),
    by phase correlation; whole pixels at this quarter size are enough for a
    brightness ratio."""
    from scipy import ndimage
    from skimage.registration import phase_cross_correlation

    def prep(a):
        a = np.where(np.isfinite(a), a, 0).astype(np.float64)
        ok = (a > LO * sat) & (a < 0.98 * sat)
        la = np.log(np.maximum(a, LO * sat))
        hp = la - ndimage.uniform_filter(la, 15)
        return np.where(ok, hp, 0.0)

    px, py = prep(x), prep(y)
    try:
        sh, _err, _ph = phase_cross_correlation(px, py, normalization=None)
    except Exception:
        return (0, 0), float("nan")
    dy, dx = int(round(sh[0])), int(round(sh[1]))
    ys = ndimage.shift(py, (dy, dx), order=0, mode="constant", cval=0.0)
    m = (px != 0) & (ys != 0)
    if m.sum() < MIN_PX:
        return (dy, dx), float("nan")
    a, b = px[m], ys[m]
    c = float(np.corrcoef(a, b)[0, 1])
    return (dy, dx), c


def _overlap(x, y, dy, dx):
    """x and y cut to the part they share once y is moved by (dy, dx)."""
    H, W = x.shape
    ys0, ys1 = max(0, -dy), min(H, H - dy)
    xs0, xs1 = max(0, -dx), min(W, W - dx)
    if ys1 - ys0 < 16 or xs1 - xs0 < 16:
        return None, None
    return x[ys0 + dy:ys1 + dy, xs0 + dx:xs1 + dx], y[ys0:ys1, xs0:xs1]


def link(x, y, sat):
    """Brightness of frame y relative to frame x (slope), and the pixel count.
    NaN when the frames share too few usable pixels."""
    if x.shape != y.shape:
        return float("nan"), 0
    # placed only when the placement is sure; otherwise compared as they are
    # (frames of one burst barely move, and a frame whose corona is mostly
    # clipped cannot be placed but still links fine where it is)
    (dy, dx), c = _align(x, y, sat)
    if np.isfinite(c) and c >= ALIGN_MIN_CORR and (dy or dx):
        x, y = _overlap(x, y, dy, dx)
        if x is None:
            return float("nan"), 0
    lo, hi = LO * sat, HI * sat
    m = np.isfinite(x) & np.isfinite(y) & (x > lo) & (y > lo) & (x < hi) & (y < hi)
    if m.sum() < MIN_PX:
        return float("nan"), int(m.sum())
    x0 = np.where(np.isfinite(x), x, 0).astype(np.float64)
    y0 = np.where(np.isfinite(y), y, 0).astype(np.float64)
    sx, sy = _box(x0), _box(y0)
    m &= (sx > lo) & (sy > lo)
    n = int(m.sum())
    if n < MIN_PX:
        return float("nan"), n
    z = (sx / np.median(sx[m]) + sy / np.median(sy[m]))[m]
    q = np.quantile(z, np.linspace(0.02, 0.98, 41))
    b = np.digitize(z, q)
    xv, yv = x0[m], y0[m]
    X, Y = [], []
    for k in range(1, len(q)):
        s = b == k
        if s.sum() > 50:
            X.append(np.median(xv[s]))
            Y.append(np.median(yv[s]))
    if len(X) < 8:
        return float("nan"), n
    k, _c = np.polyfit(np.asarray(X), np.asarray(Y), 1)
    return (float(k) if np.isfinite(k) and k > 0 else float("nan")), n


def by_brightness(paths, sec_of, ts_of=None, progress=None, gap=GAP):
    """{key: [paths]} grouped by measured brightness, and a report dict.

    `sec_of[p]` is the header exposure (used for the order and the key's
    scale only), `ts_of[p]` an optional sortable time stamp."""
    ts_of = ts_of or {}

    def _order(p):
        s = sec_of.get(p)
        return (float(s) if s and s > 0 else float("inf"), str(ts_of.get(p) or ""),
                os.path.basename(p))

    order = sorted(paths, key=_order)
    logb = {}
    prev = prev_p = None
    n_fallback = 0
    for i, p in enumerate(order):
        g, sat = _green(p)
        if prev is None:
            logb[p] = 0.0
        else:
            k, npx = link(prev, g, sat)
            if not np.isfinite(k):
                sa, sb = sec_of.get(prev_p) or 0, sec_of.get(p) or 0
                k = (sb / sa) if (sa > 0 and sb > 0) else 1.0
                n_fallback += 1
                if progress is not None:
                    progress.log(f"[warn] brightness {os.path.basename(prev_p)} -> "
                                 f"{os.path.basename(p)}: "
                                 f"{npx} shared pixels, header ratio {k:.4g} used", None)
            logb[p] = logb[prev_p] + float(np.log(k))
        prev, prev_p = g, p
        if progress is not None and (i + 1) % 25 == 0:
            progress.log(f"  measuring frame brightness: {i + 1}/{len(order)}", None)
    # split the sorted brightness run at the gaps
    #
    # CAMERA HEADERS ARE NOT MERGED ACROSS (2026-11-14). A camera's own
    # exposure time is a measurement of the shutter, not a label: frames that
    # say different times were exposed differently. Merging them is only right
    # for labels typed into FITS headers (the case this module was built
    # for). On a set of two complementary brackets 1/3 EV apart, frames from
    # the two brackets measured within the gap (light and framing changed in
    # the minute between them) and were averaged into one tier. So for camera
    # files the measurement may only SPLIT a header group (frames at another
    # ISO or aperture, a frame that got other light); FITS keep the global
    # grouping.
    from .fits import is_fits

    def _runs(ps):
        srt = sorted(ps, key=lambda p: logb[p])
        out, cur = [], [srt[0]]
        for a, b in zip(srt[:-1], srt[1:]):
            if logb[b] - logb[a] > np.log(gap):
                out.append(cur)
                cur = []
            cur.append(b)
        out.append(cur)
        return out

    fits_ps = [p for p in order if is_fits(p)]
    cam_ps = [p for p in order if not is_fits(p)]
    groups = _runs(fits_ps) if fits_ps else []
    by_sec = {}
    for p in cam_ps:
        by_sec.setdefault(round(float(sec_of.get(p) or 0), 9), []).append(p)
    for _s, ps in sorted(by_sec.items()):
        groups.extend(_runs(ps))
    groups.sort(key=lambda g: float(np.median([logb[p] for p in g])))
    med_b = [float(np.exp(np.median([logb[p] for p in g]))) for g in groups]
    scale = [np.median([sec_of[p] for p in g if (sec_of.get(p) or 0) > 0] or [np.nan]) / mb
             for g, mb in zip(groups, med_b)]
    scale = [s for s in scale if np.isfinite(s) and s > 0]
    A = float(np.median(scale)) if scale else 1.0
    tiers, rep = {}, []
    for g, mb in zip(groups, med_b):
        key = float(A * mb)
        tiers[key] = list(g)
        hdr = sorted({round(float(sec_of.get(p) or 0), 9) for p in g})
        rel = [float(np.exp(logb[p]) / mb) for p in g]
        rep.append({"key": key, "n": len(g), "header": hdr,
                    "spread_pct": round(100 * (max(rel) - min(rel)), 1)})
    # which header groups were split or merged, for the log
    by_hdr = {}
    for t in rep:
        for h in t["header"]:
            by_hdr.setdefault(h, []).append(t["key"])
    split = {h: ks for h, ks in by_hdr.items() if len(ks) > 1}
    merged = [t for t in rep if len(t["header"]) > 1]
    info = {"tiers": rep, "n_header_groups": len(by_hdr), "scale": A,
            "split": {f"{h:g}": [round(k, 6) for k in ks] for h, ks in split.items()},
            "fallback_links": n_fallback,
            "frames": {os.path.basename(p): {"header": float(sec_of.get(p) or 0),
                                             "brightness": round(float(np.exp(logb[p])), 6)}
                       for p in order}}
    if progress is not None:
        progress.log(f"tiers by measured brightness: {len(tiers)} tiers "
                     f"(headers name {len(by_hdr)} groups)", None)
        for h, ks in sorted(split.items()):
            progress.log(f"  header {h:g} s split into {len(ks)} tiers: "
                         + ", ".join(f"{k:.4g}" for k in sorted(ks)), None)
        for t in merged:
            progress.log(f"  tier {t['key']:.4g}: frames labelled "
                         + ", ".join(f"{h:g}" for h in t["header"]) + " s", None)
    return tiers, info
