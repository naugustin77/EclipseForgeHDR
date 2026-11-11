"""Corona in motion (Extras, 0.30.0): the brackets of a repeated set as an
aligned loop, the fine structure brought out so that what moved between
them can be seen.

Needs a run in the stacking mode "HDR per bracket" on a repeated-bracket
set. That run keeps every bracket's own merge as luminance in a box around
the Sun at full size (brackets_lum/win_k.npy, level and transparency matched
to the merge of all frames); older runs kept only a reduced copy of the
whole frame (hdr_k.npy), which still works, at that resolution.

Per bracket: the app's partial convolution (detail.py) -- the log-mapped
luminance minus its blur at the scales 2-32 px, blurred in polar coordinates
with the Moon left out and the radial bias removed -- the masks summed with
the same weights for every bracket and shown on mid grey: brighter than the
surroundings light, darker dark. A loop, a blob or a wave that moves, moves;
a structure that does not, stands still. The Moon is one flat grey disc over
all of them. Without detail.py the filter falls back to one isotropic high
pass of the logarithm.

THE BRACKETS ARE PUT ON THE CORONA FIRST (2026-11-10). They come placed on
one frame, but a real repeated-bracket set showed that this is not yet the sky's: a
star and the corona's own structure both walked 7.9 px across its seven
brackets (78 s), while the grain at the scale of a pixel stayed put (its
correlation between two brackets 0.4-0.5 at zero shift) -- the camera's
pattern, which a registration on fine structure holds on to. The first
version here measured the shift on the fine masks and found 0.24 px. Now:
the log luminance between 5 and 20 px (nothing finer, so no pattern), in
the ring 1.15-2.6 radii, whole-pixel peak then least squares, against the
middle bracket and then the mean; the luminance is shifted (cubic) before
any filter. On that set the star then stays within 0.8 px (0.45 rms, its
own centroid noise) and comes out round (3.0 x 2.8 px against 8.2 x 3.8 px
in the mean of the brackets as placed).

Outputs: an animated GIF and an MP4 of the loop, and the card's still.
Nothing is added or modelled: every pixel is the bracket's own, filtered.
"""
import json
import math
import os

import numpy as np
from PIL import Image
from scipy import ndimage

from . import extras as _ex

OPTIONS = [
    dict(id="scale", label="Detail", type="select", default="look",
         values=[["look", "as in the look"], ["fine", "fine"], ["coarse", "coarse"]],
         title="The partial convolution's scales (2, 4, 8, 16, 32 px) and their weights: as the look in tab 3 "
               "sets them (the five scale sliders and the log stretch), fine (2 and 4 px: thin threads, and "
               "the grain with them), or coarse (8-32 px: broad loops and streamers)."),
    dict(id="grain", label="Grain", type="select", default="soft",
         values=[["keep", "keep"], ["soft", "softened"], ["clean", "clean"]],
         title="Keep: every frame is its bracket alone. Softened: every frame is the mean of its bracket and "
               "the two beside it (the bracket itself counted twice, each weighted by its own noise): about "
               "0.6 of the grain, all the threads, and what moves still moves. Clean: that, and what stands "
               "less than 0.7 of the remaining noise above its surroundings fades to grey (faint threads "
               "can go with it)."),
    dict(id="contrast", label="Contrast", type="select", default="normal",
         values=[["low", "low"], ["normal", "normal"], ["high", "high"]],
         title="How far from mid grey a structure goes: the same for every bracket."),
    dict(id="view", label="View", type="select", default="quad",
         values=[["quad", "four quadrants, 1:1"], ["whole", "the whole corona"],
                 ["0", "zoom at PA 0°"], ["45", "zoom at PA 45°"], ["90", "zoom at PA 90°"],
                 ["135", "zoom at PA 135°"], ["180", "zoom at PA 180°"], ["225", "zoom at PA 225°"],
                 ["270", "zoom at PA 270°"], ["315", "zoom at PA 315°"]],
         title="Four quadrants: four loops, each a quarter of the corona at full resolution with the Moon in "
               "its corner. The whole corona: one loop at the chosen size. A zoom: a 1:1 window beside the limb "
               "at this position angle (from solar north when the star match gives it)."),
    dict(id="turn", label="Solar north up", type="check", default=True,
         title="Turns the frames so that solar north is up (needs the star match or north by hand)."),
    dict(id="pace", label="Pace", type="select", default="normal",
         values=[["slow", "slow (0.8 s per bracket)"], ["normal", "normal (0.45 s)"], ["fast", "fast (0.25 s)"]],
         title="How long each bracket stays in the loop; between two brackets a short dissolve."),
    dict(id="size", label="Size", type="select", default="1200",
         values=[["800", "800 px"], ["1200", "1200 px"], ["1800", "1800 px"]],
         title="The whole corona: the loop's size (square). Quadrants and zooms: the largest size; smaller "
               "when the data has fewer pixels (never enlarged beyond 1:1)."),
    dict(id="steady", label="Steady on the corona", type="check", default=True,
         title="Each bracket moved so that its corona lies on the others': measured on the corona's own "
               "structure of 5 to 20 px, not on the grain. The info says how far apart the brackets lay."),
    dict(id="stamp", label="Time stamp", type="check", default=True,
         title="The seconds from the first bracket, small in a corner."),
]
NEEDS_RERUN = ("the brackets on disk are reduced copies from an older run: stack again in the mode "
               "'HDR per bracket' to keep them at full size")
PACE = {"slow": 0.8, "normal": 0.45, "fast": 0.25}
SIGMA = {"fine": 0.006, "normal": 0.012, "coarse": 0.025}       # high pass, in solar radii
CONTRAST = {"low": 7.0, "normal": 5.0, "high": 3.5}             # +- this many robust sigma to black/white


def available(wd):
    """The brackets the run kept: [(k, path)] -- the full-size boxes when
    there are some, else the reduced whole frames."""
    d = os.path.join(wd, "brackets_lum")
    if not os.path.isdir(d):
        return []
    names = os.listdir(d)
    pre = "win_" if any(n.startswith("win_") for n in names) else "hdr_"
    out = []
    for n in names:
        if n.startswith(pre) and n.endswith(".npy"):
            k = n[len(pre):-4]
            out.append((int(k) if k.isdigit() else k, os.path.join(d, n)))
    out.sort(key=lambda e: (str(type(e[0])), e[0]))
    return out


def _index(wd):
    try:
        return json.load(open(os.path.join(wd, "brackets_lum", "index.json")))
    except Exception:
        return {}


def _hp(lg, valid, sig):
    """High pass of the log luminance: lg minus its normalised blur over the
    valid pixels (outside the Moon, seen by the bracket)."""
    m = valid.astype(np.float32)
    num = ndimage.gaussian_filter(np.where(valid, lg, 0.0).astype(np.float32), sig, mode="nearest")
    den = ndimage.gaussian_filter(m, sig, mode="nearest")
    return np.where(valid, lg - num / np.maximum(den, 1e-3), 0.0).astype(np.float32)


_HD = (1.0, 0.6, 0.2, 0.1, 0.05)                         # the renderer's default ladder, 2..32 px
PC_GAINS = {"fine": (0.5, 1.0, 0.0, 0.0, 0.0),
            "coarse": (0.0, 0.3, 1.0, 0.7, 0.4)}


class _PC:
    """The app's partial convolution (detail.py) on one bracket, with the
    look's settings: the log mapping at the look's stretch (10^hillLogK),
    the polar partial blur per scale with the Moon left out and the radial
    bias removed, the masks summed in im_log units with the look's five
    scale weights (normalised to the default ladder's sum, as the renderer
    does). The mapping's top is taken once, so every bracket gets the same."""

    def __init__(self, mode, rd, cy, cx, look=None):
        from . import detail as _d
        self.d = _d
        look = look or {}
        if mode in PC_GAINS:
            g = PC_GAINS[mode]
        else:
            g = [float(look.get("hill%d" % i, d_)) for i, d_ in enumerate(_HD)]
        ssum = sum(g)
        g = [x * (sum(_HD) / ssum) for x in g] if ssum > 1e-6 else list(_HD)
        self.sc = [(max(0.7, float(s_) / rd), float(gn)) for s_, gn in zip(_d.HILL_SCALES, g) if gn > 1e-6]
        self.full_sc = [float(s_) for s_, gn in zip(_d.HILL_SCALES, g) if gn > 1e-6]
        self.K = 10.0 ** float(look.get("hillLogK", 5.0))
        self.cy, self.cx = cy, cx
        self.pk = None

    def __call__(self, lum, valid, rr, ring):
        if self.pk is None:
            self.pk = max(float(np.percentile(lum[valid], 99.999)) if valid.any() else 1.0, 1e-9)
        xn = np.clip(np.nan_to_num(lum, nan=0.0) / self.pk, 0, None)
        imlog = np.clip(np.log1p(self.K * xn) / np.log1p(self.K), 0.0, 1.0).astype(np.float32)
        w = valid.astype(np.float32)
        blurs, _cov = self.d.polar_partial_blur(imlog, w, self.cy, self.cx, [s_ for s_, _g in self.sc], order=2)
        out = np.zeros_like(imlog)
        for b, (_s, gn) in zip(blurs, self.sc):
            m = self.d._deradial((imlog - b).astype(np.float32), rr, valid)
            out += np.float32(gn) * np.where(valid, m, 0.0).astype(np.float32)
        del blurs, imlog
        return out


def _band(a, valid, s1, s2):
    """The corona's structure between s1 and s2 pixels: a difference of two
    normalised blurs of the log luminance over the valid pixels. Nothing finer
    than s1 is left in it: what a camera prints at the scale of a pixel stays
    on the sensor from bracket to bracket and would hold a registration
    there, while the sky moves."""
    fl = float(np.percentile(a[valid], 0.5)) if valid.any() else 1e-6
    lg = np.log(np.maximum(np.nan_to_num(a, nan=fl), max(fl, 1e-9))).astype(np.float32)
    m = valid.astype(np.float32)

    def nb(s_):
        return (ndimage.gaussian_filter(np.where(valid, lg, 0.0), s_, mode="nearest")
                / np.maximum(ndimage.gaussian_filter(m, s_, mode="nearest"), 1e-3))
    return np.where(valid, nb(s1) - nb(s2), 0.0).astype(np.float32)


def _coarse_shift(a, ref, win, reach):
    """Whole-pixel (dy, dx) that moves `a` onto `ref`: the peak of the
    cross-correlation over the window, within `reach` pixels."""
    A = np.fft.rfft2(a * win)
    B = np.fft.rfft2(ref * win)
    c = np.fft.irfft2(B * np.conj(A), s=a.shape)
    r = int(max(1, min(reach, min(a.shape) // 4)))
    sub = np.concatenate([np.concatenate([c[-r:, -r:], c[-r:, :r + 1]], axis=1),
                          np.concatenate([c[:r + 1, -r:], c[:r + 1, :r + 1]], axis=1)], axis=0)
    i = np.unravel_index(int(np.argmax(sub)), sub.shape)
    return float(i[0] - r), float(i[1] - r)


def _fine_shift(a, ref, m, d0=(0.0, 0.0), it=5):
    """Sub-pixel (dy, dx) that moves `a` onto `ref` over the mask, from d0:
    the least-squares shift on the reference's gradients, repeated."""
    gy, gx = np.gradient(ref)
    gy, gx = gy[m].astype(np.float64), gx[m].astype(np.float64)
    A = np.array([[(gy * gy).sum(), (gy * gx).sum()], [(gy * gx).sum(), (gx * gx).sum()]])
    d = np.array(d0, np.float64)
    if abs(np.linalg.det(A)) < 1e-30:
        return d
    for _ in range(it):
        e = (ndimage.shift(a, d, order=1, mode="nearest") - ref)[m].astype(np.float64)
        st = np.linalg.solve(A, np.array([(gy * e).sum(), (gx * e).sum()]))
        d += np.clip(st, -3.0, 3.0)
        if float(np.hypot(*st)) < 0.01:
            break
    return d


#: the registration's band, stack px: nothing finer than the first, nothing
#: coarser than the second
STEADY_BAND = (5.0, 20.0)


def register(lums, valids, cy, cx, R, r_moon, q=1.0):
    """The shift (dy, dx), working px, that puts every bracket's corona on the
    same place: the brackets' structure between 5 and 20 stack px, in the
    ring from 1.15 to 2.6 lunar radii, each against the middle bracket and
    then against the mean of all; the shifts add up to zero, so the set
    stays where the picture's geometry is. Returns (shifts, left) where
    `left` is each bracket's shift still measured afterwards, working px."""
    n = len(lums)
    h, w = lums[0].shape
    # worked at half size where that still leaves the band two pixels wide
    dec = 2 if (STEADY_BAND[0] / q) >= 4.0 and min(h, w) >= 1200 else 1
    s1 = max(1.5, STEADY_BAND[0] / q / dec)
    s2 = max(4.0 * s1, STEADY_BAND[1] / q / dec)

    def small(a):
        if dec == 1:
            return a
        hh, ww = (a.shape[0] // 2) * 2, (a.shape[1] // 2) * 2
        return a[:hh, :ww].reshape(hh // 2, 2, ww // 2, 2).mean(axis=(1, 3))
    L = [small(np.nan_to_num(np.asarray(a, np.float32), nan=0.0)) for a in lums]
    V = [small(v_.astype(np.float32)) > 0.99 for v_ in valids]
    hh, ww = L[0].shape
    yy = np.arange(hh, dtype=np.float32)[:, None] - (cy / dec - (0.25 if dec == 2 else 0.0))
    xx = np.arange(ww, dtype=np.float32)[None, :] - (cx / dec - (0.25 if dec == 2 else 0.0))
    rr = np.sqrt(yy * yy + xx * xx)
    allv = np.all(V, axis=0)
    zone = (rr > max(1.15 * R, r_moon) / dec) & (rr < 2.6 * R / dec)
    m = zone & ndimage.binary_erosion(allv, iterations=int(math.ceil(s2)))
    if m.sum() < 2000:
        return [(0.0, 0.0)] * n, [(0.0, 0.0)] * n
    Hb = [_band(a, v_, s1, s2) for a, v_ in zip(L, V)]
    win = ndimage.gaussian_filter(m.astype(np.float32), 2.0)
    mid = n // 2
    d = []
    for k in range(n):
        if k == mid:
            d.append(np.zeros(2))
            continue
        c0 = _coarse_shift(Hb[k], Hb[mid], win, reach=40.0 / dec / q + 4)
        d.append(_fine_shift(Hb[k], Hb[mid], m, c0))
    # once more, against the mean of all as registered
    for _pass in range(2):
        ref = np.mean([ndimage.shift(Hb[k], d[k], order=1, mode="nearest") for k in range(n)], axis=0)
        d = [_fine_shift(Hb[k], ref, m, d[k], it=3) for k in range(n)]
    mean_d = np.mean(d, axis=0)
    d = [x - mean_d for x in d]
    ref = np.mean([ndimage.shift(Hb[k], d[k], order=1, mode="nearest") for k in range(n)], axis=0)
    left = [_fine_shift(ndimage.shift(Hb[k], d[k], order=1, mode="nearest"), ref, m, (0.0, 0.0), it=2)
            for k in range(n)]
    return ([(float(x[0]) * dec, float(x[1]) * dec) for x in d],
            [(float(x[0]) * dec, float(x[1]) * dec) for x in left])


def motion(wd, geo, folder, v, out_base, prog=None, preview=False):
    log = (lambda m, f=None: prog.log(m, f) if prog else None)
    br = available(wd)
    if len(br) < 2:
        raise ValueError("needs a run in the stacking mode 'HDR per bracket' on a repeated-bracket set")
    idx = _index(wd)
    full = br[0][1].split(os.sep)[-1].startswith("win_")
    if not full:
        # an older run kept only reduced copies of the brackets: too coarse
        # for fine structure, so no loop is made from them
        raise ValueError(NEEDS_RERUN)
    recs = {int(r["bracket"]): r for r in idx.get("brackets", []) if "bracket" in r}
    if full:
        y0b, y1b, x0b, x1b = idx["win_box"]
        step = 1
    else:
        y0b = x0b = 0
        step = int(idx.get("step") or 1)
    # preview: every 2nd or 4th pixel
    rd = 1 if not preview else (4 if full else 2)
    if full and not preview and v.get("view", "whole") == "whole":
        # the whole corona is shown at the loop's size: worked at about that
        # size, not at full size (the polar blur is the cost, and what is
        # finer than an output pixel would not show)
        rd = max(1, int(min(y1b - y0b, x1b - x0b) / (1.15 * float(v.get("size", 1200)))))
    q = step * rd                                           # stack px per working px
    cy = (float(geo["cy"]) - y0b) / q
    cx = (float(geo["cx"]) - x0b) / q
    R = float(geo["R"]) / q
    Rm = float(geo.get("Rmask", geo["R"])) / q
    ctx = _ex.context(wd, folder, geo, scale=1.0, **_ex.angle_kw(v))
    Rs = ctx["Rs"] / q
    sig = max(1.0, SIGMA.get(v.get("scale", "normal"), 0.012) * Rs)
    if preview:
        br = [br[0], br[len(br) // 2], br[-1]] if len(br) > 2 else br
    # the Moon over all brackets: the largest offset of a bracket's Moon, and
    # the partial convolution's own pad beyond the limb
    offs = [float(r.get("moon_from_target_px") or 0.0) for r in recs.values()]
    r_moon = Rm + (max(offs) if offs else 0.0) / q + 2.0
    try:
        pc = _PC(v.get("scale", "look"), q, cy, cx, look=v.get("_look"))
        r_moon += float(getattr(pc.d, "HILL_LIMB_PAD", 5.0)) / q
        how = ("partial convolution " + "/".join(f"{s_:g}" for s_ in pc.full_sc) + " px"
               + (" (the look's weights)" if v.get("scale", "look") == "look" else ""))
    except Exception:
        pc = None
        how = f"high pass {sig * q:.1f} px"
    log(f"motion: {len(br)} brackets, {how} ...", 0.1)
    lums, ts = [], []
    for k, p in br:
        a = np.load(p, mmap_mode="r")
        a = np.asarray(a[::rd, ::rd], np.float32) if a.ndim == 2 else \
            np.asarray(0.2126 * a[::rd, ::rd, 0] + 0.7152 * a[::rd, ::rd, 1] + 0.0722 * a[::rd, ::rd, 2], np.float32)
        rec = recs.get(int(k)) if str(k).isdigit() else None
        lums.append(a)
        ts.append(float(rec.get("t", 0.0)) if rec else float(len(ts)))
    order = np.argsort(ts)
    lums = [lums[i] for i in order]
    ts = sorted(ts)
    ts = [t_ - ts[0] for t_ in ts]
    t_span = ts[-1] - ts[0] if ts else 0.0
    h, w = lums[0].shape
    yy = np.arange(h, dtype=np.float32)[:, None] - cy
    xx = np.arange(w, dtype=np.float32)[None, :] - cx
    rr = np.sqrt(yy * yy + xx * xx)
    ring = (rr > 1.12 * R) & (rr < 2.6 * R)
    seen = [np.isfinite(a) & (a > 0) for a in lums]
    # steady on the corona. The brackets come placed on one frame, but on
    # that frame the sky may still drift from bracket to bracket; here each
    # bracket is put where its corona is, on the corona's own structure of
    # 5 to 20 px, before anything is filtered.
    shifts = [(0.0, 0.0)] * len(lums)
    left = list(shifts)
    if v.get("steady", True) and len(lums) >= 2:
        log("motion: putting the brackets on the corona ...", 0.12)
        shifts, left = register(lums, [s_ & (rr > r_moon) for s_ in seen], cy, cx, R, r_moon, q)
        for i, s_ in enumerate(shifts):
            if abs(s_[0]) + abs(s_[1]) < 1e-3:
                continue
            fl = np.nan_to_num(lums[i], nan=0.0)
            lums[i] = np.maximum(ndimage.shift(fl, s_, order=3, mode="nearest"), 0.0).astype(np.float32)
            # what was not seen stays out, with the spline's reach round it
            ns = ndimage.shift((~seen[i]).astype(np.float32), s_, order=1, mode="constant", cval=1.0) > 0.01
            seen[i] = ~ndimage.binary_dilation(ns, iterations=2)
    apart = max(float(np.hypot(a_[0] - b_[0], a_[1] - b_[1])) for a_ in shifts for b_ in shifts) * q
    frames = []
    for a, sn in zip(lums, seen):
        valid = sn & (rr > r_moon)
        if pc is not None:
            frames.append((pc(a, valid, rr, ring), valid))
        else:
            fl = float(np.percentile(a[valid], 0.5)) if valid.any() else 1e-6
            lg = np.log(np.maximum(np.nan_to_num(a, nan=fl), max(fl, 1e-9)))
            frames.append((_hp(lg, valid, sig), valid))
            del lg
        log(f"motion: bracket {len(frames)} of {len(lums)} ...", 0.15 + 0.4 * len(frames) / len(lums))
    del lums, seen
    # THE GRAIN. Measured on a real 7-bracket set (2026-11-10): a gain per
    # pixel from its neighbourhood's spread (what "softened" was) took the
    # faint threads down with the grain and left the rest in blotches. What
    # the brackets have instead is each other: the corona hardly changes from
    # one bracket to the next, the grain is new in each. So "softened" is each
    # frame as the mean of its bracket (weight 2) and the two beside it
    # (weight 1 each), every bracket weighted by its own noise: 0.6 of the
    # grain, every thread kept, and what moves still moves (it leaves a
    # quarter of itself one bracket to each side). "Clean" adds the gain per
    # pixel on top, at 0.7 of the noise that is left.
    grain = v.get("grain", "soft")
    rb = np.clip((rr / 4.0).astype(np.int32), 0, None)
    nbin = int(rb.max()) + 1
    n = len(frames)
    allm = np.all([m for _f, m in frames], axis=0)
    grain_note = ""
    if grain in ("soft", "clean") and n >= 2:
        # every bracket's own noise per radius, from the differences between
        # brackets two or more apart (next neighbours still share some of the
        # sensor's pattern after the shift): spread(i - j)^2 = n_i^2 + n_j^2
        pairs = [(i, j) for i in range(n) for j in range(i + 2, n)] or [(i, i + 1) for i in range(n - 1)]
        sub = allm & ring
        bsel = rb[sub]
        order_ = np.argsort(bsel, kind="stable")
        edges = np.searchsorted(bsel[order_], np.arange(nbin + 1))
        Am = np.zeros((len(pairs), n))
        Ym = np.full((len(pairs), nbin), np.nan)
        for pi, (i, j) in enumerate(pairs):
            Am[pi, i] = Am[pi, j] = 1.0
            dv = np.abs(frames[j][0][sub] - frames[i][0][sub])[order_]
            for b0 in range(0, nbin, 4):
                seg = dv[edges[b0]:edges[min(nbin, b0 + 4)]]
                if seg.size > 200:
                    Ym[pi, b0:b0 + 4] = (float(np.median(seg[::max(1, seg.size // 20000)])) / 0.6745) ** 2
            del dv
        ok = np.all(np.isfinite(Ym), axis=0)
        sig2 = None
        if ok.sum() >= 2:
            sol = np.linalg.lstsq(Am, Ym[:, ok], rcond=None)[0]              # n x bins
            tot = Ym[:, ok].mean(axis=0) / 2.0
            sol = np.clip(sol, 0.25 * tot, 4.0 * tot)                         # no bracket without noise
            xs = np.nonzero(ok)[0]
            sig2 = np.stack([ndimage.uniform_filter1d(np.interp(np.arange(nbin), xs, sol[k]), 9, mode="nearest")
                             for k in range(n)]).astype(np.float32)
        if sig2 is not None:
            left2 = sig2.copy()                                               # the noise left in each frame
            if n >= 3:
                W = np.zeros((n, n))
                for k in range(n):
                    if 0 < k < n - 1:
                        W[k, k - 1], W[k, k], W[k, k + 1] = 1.0, 2.0, 1.0
                    elif k == 0:
                        W[k, 0], W[k, 1], W[k, 2] = 2.0, 1.5, 0.5             # the ends lean inward: the same grain
                    else:
                        W[k, n - 1], W[k, n - 2], W[k, n - 3] = 2.0, 1.5, 0.5
                src = [f for f, _m in frames]
                out_f, ts_new = [], []
                for k in range(n):
                    js = [j for j in range(n) if W[k, j] > 0]
                    wj = [np.float32(W[k, j]) / sig2[j] for j in js]          # per radius
                    den = np.sum(wj, axis=0)
                    acc_ = np.zeros_like(src[0])
                    for j, w_ in zip(js, wj):
                        acc_ += (w_ / den)[rb] * src[j]
                    out_f.append((np.where(allm, acc_, 0.0).astype(np.float32), allm))
                    left2[k] = np.sum([(w_ / den) ** 2 * sig2[j] for j, w_ in zip(js, wj)], axis=0)
                    mid_ = nbin // 2
                    wm = np.array([float((w_ / den)[mid_]) for w_ in wj])
                    ts_new.append(float(np.sum(wm * np.array([ts[j] for j in js]))))
                frames = out_f
                ts = ts_new
                del src
                grain_note = ", each frame the mean of its bracket and the two beside it"
            if grain == "clean" or n < 3:
                lam = 0.7 if n >= 3 else (1.0 if grain == "soft" else 2.0)
                sw_ = max(1.5, 1.5 * (pc.sc[0][0] if pc is not None else 1.0))
                out_f = []
                for (f, m), s2k in zip(frames, left2):
                    s2 = (lam * lam * s2k)[rb]
                    V = ndimage.gaussian_filter(f * f, sw_)
                    gk = np.clip((V - s2) / np.maximum(V, 1e-20), 0.0, 1.0)
                    out_f.append((f * ndimage.gaussian_filter(gk, 1.0), m))
                frames = out_f
                grain_note += f", what stands less than {lam:g} x the noise faded"
    # one contrast for all: the 99.5th percentile of the structure in the ring
    pool = np.concatenate([np.abs(f[ring & m])[::5] for f, m in frames])
    p998 = float(np.percentile(pool, 99.8)) if pool.size else 1.0
    lim = max(p998, 1e-9) * {"low": 2.0, "normal": 1.4, "high": 1.0}.get(v.get("contrast", "normal"), 1.4)
    disc = rr <= r_moon + (0.0 if pc is not None else 0.5 * sig)   # the limb's own step under the disc
    # the views, in the frame turned so that solar north is up
    rot = 0.0
    if v.get("turn", True) and ctx.get("north_solar_cw") is not None:
        rot = float(ctx["north_solar_cw"])
    elif v.get("turn", True) and ctx.get("north_cel_cw") is not None:
        rot = float(ctx["north_cel_cw"])
    t_ = math.radians(rot)
    ca, sa = math.cos(t_), math.sin(t_)
    S = int(v.get("size", 1200))
    vw = v.get("view", "quad")
    views = []                                              # (name, centre dy, dx in the turned frame, half)
    if vw == "whole":
        half = min(cy, h - 1 - cy, cx, w - 1 - cx)
        half /= max(1.0, abs(ca) + abs(sa))                 # turned: no empty corners
        views.append(("", 0.0, 0.0, half, S))
    elif vw == "quad":
        hq = 1.1 * Rs                                       # each quarter reaches 2.2 radii
        compass = bool(v.get("turn", True)) and ctx["o"].get("origin") in ("solar", "celestial")
        east_left = ctx["o"].get("sgn", 1) == 1
        for dy, dx, nm in ((-1, -1, "upper-left"), (-1, 1, "upper-right"), (1, -1, "lower-left"), (1, 1, "lower-right")):
            if compass:
                nm = ("N" if dy < 0 else "S") + (("E" if dx < 0 else "W") if east_left else ("W" if dx < 0 else "E"))
            views.append((nm, dy * hq, dx * hq, hq, min(S, int(round(2 * hq)))))
    else:
        pa = float(vw)
        phi = math.radians(ctx["o"]["zero_cw"] - ctx["o"]["sgn"] * pa)
        # the direction in the picture, then into the turned frame
        py_, px_ = -math.cos(phi) * 1.75 * Rs, math.sin(phi) * 1.75 * Rs
        dxo = px_ * ca + py_ * sa
        dyo = -px_ * sa + py_ * ca
        hz = 0.85 * Rs
        views.append((f"pa{int(pa)}", dyo, dxo, hz, min(S, int(round(2 * hz)))))
    if preview:
        views = [(nm, dy, dx, hf, min(Sv, 700 if len(views) == 1 else 350)) for nm, dy, dx, hf, Sv in views]

    # near 1:1 the frames are resampled with a cubic spline (bilinear takes
    # the finest scale down by a different amount at every pixel once the
    # view is turned); the spline's coefficients once per frame
    cubic = (not preview) and all((2 * hf / Sv_) < 1.5 for _n, _dy, _dx, hf, Sv_ in views)
    coef = [ndimage.spline_filter(f, order=3, output=np.float32, mode="mirror") for f, _m in frames] if cubic else None

    def render_view(dy, dx, half, Sv):
        oy, ox = np.mgrid[0:Sv, 0:Sv].astype(np.float32)
        u = dx + (ox - Sv / 2.0) * (2 * half / Sv)
        vv = dy + (oy - Sv / 2.0) * (2 * half / Sv)
        X = cx + u * ca - vv * sa
        Y = cy + u * sa + vv * ca
        dsk = ndimage.map_coordinates(disc.astype(np.float32), [Y, X], order=1, mode="constant", cval=0.0)
        res = []
        fs = max(10, Sv / 40.0)
        for i_, ((f, m), t) in enumerate(zip(frames, ts)):
            if coef is not None:
                g = ndimage.map_coordinates(coef[i_], [Y, X], order=3, mode="mirror", prefilter=False)
            else:
                g = ndimage.map_coordinates(f, [Y, X], order=1, mode="constant", cval=0.0)
            mv = ndimage.map_coordinates(m.astype(np.float32), [Y, X], order=1, mode="constant", cval=0.0)
            img = 0.5 + 0.5 * np.clip(g / lim, -1, 1) * mv
            img = img * (1 - dsk) + 0.5 * dsk               # the Moon: one flat grey disc
            im = Image.fromarray((np.clip(img, 0, 1) * 255 + 0.5).astype(np.uint8)).convert("RGB")
            if v.get("stamp", True):
                d = _ex.Draw(im)
                d.text((Sv - fs, Sv - fs), f"+{t:.0f} s", fill=(235, 235, 235), font=_ex._font(fs), anchor="rs")
            res.append(im)
        return res

    if preview:
        stills = [render_view(dy, dx, hf, Sv)[len(frames) // 2] for _nm, dy, dx, hf, Sv in views]
        if len(stills) == 4:
            Sv = stills[0].width
            mos = Image.new("RGB", (2 * Sv + 6, 2 * Sv + 6), (128, 128, 128))
            for i, st in enumerate(stills):
                mos.paste(st, ((i % 2) * (Sv + 6), (i // 2) * (Sv + 6)))
            if mos.width < 600:                             # a small frame: the card still needs a picture
                kq = int(math.ceil(600.0 / mos.width))
                mos = mos.resize((mos.width * kq, mos.height * kq), Image.LANCZOS)
            return mos, [], f"preview: {len(frames)} of the brackets, four quadrants"
        st0 = stills[0]
        if st0.width < 500:
            kq = int(math.ceil(500.0 / st0.width))
            st0 = st0.resize((st0.width * kq, st0.height * kq), Image.LANCZOS)
        return st0, [], f"preview: {len(frames)} of the brackets"
    # the loops: hold, then a short dissolve into the next, back to the first
    fps = 20
    hold = max(1, int(round(PACE.get(v.get("pace", "normal"), 0.45) * fps)))
    nx = max(1, int(round(0.12 * fps)))
    paths = []
    still = None
    for vi, (nm, dy, dx, hf, Sv) in enumerate(views):
        log(f"motion: writing {('the ' + nm + ' quadrant') if vw == 'quad' else 'the loop'} ...",
            0.6 + 0.35 * vi / len(views))
        pics = render_view(dy, dx, hf, Sv)
        if still is None:
            still = pics[len(pics) // 2]
        seq = []
        for i in range(len(pics)):
            a_, b_ = pics[i], pics[(i + 1) % len(pics)]
            seq += [a_] * hold
            for j in range(1, nx + 1):
                seq.append(Image.blend(a_, b_, j / (nx + 1.0)))
        suf = "_motion" + (f"_{nm}" if nm else "")
        gp = out_base + suf + ".gif"
        kk = 0
        while os.path.exists(gp):
            kk += 1
            gp = out_base + f"{suf}_{kk}.gif"
        # the GIF: every bracket once with its holding time (not a copy per
        # video frame), the dissolve frames between; a quarter of the size
        gfr, gdur = [], []
        for i in range(len(pics)):
            a_, b_ = pics[i], pics[(i + 1) % len(pics)]
            gfr.append(a_.convert("L")); gdur.append(int(1000 * hold / fps))
            for j in range(1, nx + 1):
                gfr.append(Image.blend(a_, b_, j / (nx + 1.0)).convert("L")); gdur.append(int(1000 / fps))
        gfr[0].save(gp, save_all=True, append_images=gfr[1:], duration=gdur, loop=0, optimize=False)
        paths.append(gp)
        try:
            from . import film as _film
            mp = os.path.splitext(gp)[0] + ".mp4"
            E = Sv - Sv % 2
            wr = _film._Writer(mp, E, E, fps)
            if wr.kind == "mp4":
                for _rep in range(2):
                    for im in seq:
                        wr.add(im.crop((0, 0, E, E)))
                paths.append(wr.close())
        except Exception as e:
            log(f"motion: no MP4 ({e})")
    info = (f"{len(frames)} brackets over {t_span:.0f} s, {how}"
            + grain_note
            + (f", brackets put on the corona (they lay up to {apart:.1f} px apart)"
               if v.get("steady", True) else "")
            + (f", {len(views)} loops of {views[0][4]} px" if len(views) > 1 else f", {views[0][4]} px")
            + ("" if full else ", from the reduced copies of an older run (run 'HDR per bracket' again for full size)"))
    return still, paths, info


def make(wd, geo, folder, opts, out_base, prog=None, preview=False):
    """Returns (still PIL image, extra paths, info)."""
    v = {}
    for o in OPTIONS:
        x = (opts or {}).get(o["id"], o["default"])
        if o["type"] == "check":
            x = x if isinstance(x, bool) else str(x).lower() in ("1", "true", "on", "yes")
        v[o["id"]] = x
    v["_look"] = (opts or {}).get("_look") or {}
    return motion(wd, geo, folder, v, out_base, prog=prog, preview=preview)
