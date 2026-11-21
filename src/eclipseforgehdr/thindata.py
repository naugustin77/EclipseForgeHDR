"""THIN DATA, shown as "Limited dataset" (2026-11-18): a set whose light is too thin for the usual looks.

A beginner's bracket -- a compact camera, nothing longer than a few
milliseconds, a bracket of 2-3 stops -- stacks and aligns, but the looks are
made for sets that reach the outer corona. On such a set they lift the sky's
noise into the mid-tones (Clean Steel M4 put the sky at 95 of 255 with a
spread of 18 on the test set) and the detail layers turn that noise into a
mottled band where the corona fades out. The merge itself is smooth there.

`detect()` reads four measurements the run already makes (report.json) and
calls a set thin when at least two of them say so. Measured on ten sets: the
thin test set meets all four, none of the other nine meets any.

    signal                            thin test set   other nine sets
    bracket span (EV)                 2.6             7.0 - 12.9
    corona traced to (R)              1.70            2.10 - 6.54
    merged samples below black        34 %            0 %
    merge noise scale                 5.9             0.47 - 1.68

`look()` gives the Limited dataset look: Clean Steel M4 with the detail layers
(and the FNRGF, partial-convolution and inner-corona mixes), structure and
sharpening off, no sky fill and no radial flatten, and a black
point measured on this set -- the lowest at which the sky's median reaches
black (2 of 255), so the noise goes to black and the corona keeps what is
above it. On the test set that is 0.156: at full size the corona reads 150
at 1.1 R, 87 at 1.25 R, 39 at 1.45 R, the sky 0.6.
"""
from __future__ import annotations
import os, json, math
import numpy as np

SPAN_EV = 5.0          # a bracket narrower than this cannot reach the outer corona
EXTENT_R = 2.0         # the corona traced to less than this
NEG_FRAC = 0.05        # more than this share of merged samples below black
NOISE_SCALE = 3.0      # the merge's noise above the single-exposure model
MIN_SIGNALS = 2

SKY_TARGET = 2.0 / 255.0   # the sky's median at the black point

# Clean Steel M4 (gui.html LOOKS.expM4), then the thin-data changes
_M4 = {"baseLift": 0.45, "bgNeutral": 1, "clarity": 0.1, "containWhites": 1, "coronaLevel": 0.62,
       "detailGain": 1.3, "detailScale": 1, "discSoft": 12, "earthShine": 0.3, "envGamma": 1.1,
       "finalDark": 0.2, "finalLight": 0, "finalMid": 0.35, "finalNoise": 6, "fnCompress": 0.9,
       "fnMix": 0.4, "hillBase": 0, "hillDenoise": 0, "hillGain": 0.02, "hillLogK": 4, "hillMix": 0.22,
       "hillZone": 0.3, "hillZoneR": 2, "hlCompress": 1, "hlDesat": 0, "innerDenoise": 0.4,
       "innerDim": 0, "innerMix": 0.21, "mgnContrast": 0.0, "moonColour": 0.305,
       "moonCurve": [[0, 0], [1, 1]], "moonHueB": 1.0, "moonHueG": 0.766, "moonHueR": 0.696,
       "outGamma": 1.058, "pelGain": 0, "promChroma": 0, "promDetail": 0.7, "promFade": 1,
       "promGain": 0.4, "promHueB": 0.75, "promHueG": 0.7, "promHueR": 1, "promLayer": 1,
       "promPeak": 1, "promRim": 1, "promSat": 0.5, "promStretch": 0.5, "radialFlatten": 0.35,
       "satur": 0.15, "skyColour": 2.477, "skyHueB": 1.0, "skyHueG": 0.802, "skyHueR": 0.477,
       "skyLevel": 0.093, "smoothing": 0.25, "structure": 0.4,
       "toneCurve": [[0, 0], [0.25, 0.42], [0.5, 0.7], [0.75, 0.88], [1, 1]], "bgBlack": 0.1,
       "discLevel": 0.209, "warm": [1.399, 1.041], "innerDetail": 0, "bandDetail": 2.2,
       "hill0": 0.0, "hill1": 1.0, "hill2": 1.0, "hill3": 0.7, "hill4": 0.4, "bandCoarse": 0.0}
# fnMix, hillMix and innerMix off as well (2026-11-19): those three layers
# carried most of the grain into the picture. Measured at full size on the
# thin test set (fine-grain spread, levels of 255, at 1.2 R / 1.4 R / sky):
# 17.5 / 24.9 / 24.4 with them, 8.3 / 11.2 / 6.0 without. Denoise when
# stacking (fine) changed the same numbers by only 3-5 %, because it acts on
# the detail layers, which this look already keeps low.
THIN_CHANGES = {"detailGain": 0.35, "bandDetail": 0, "innerDetail": 0, "structure": 0,
                "finalMid": 0, "finalDark": 0, "skyLevel": 0, "skyColour": 0,
                "radialFlatten": 0.0, "mgnContrast": 0,
                "fnMix": 0.0, "hillMix": 0.0, "innerMix": 0.0}

_MEMO = {}


def detect(wd):
    """{"thin": bool, "why": [lines], "signals": {...}} from report.json."""
    out = {"thin": False, "why": [], "signals": {}}
    try:
        r = json.load(open(os.path.join(wd, "report.json")))
    except Exception:
        return out
    sig = {}
    t = [float(x["sec"]) for x in (r.get("tiers") or []) if x.get("sec")]
    if len(t) >= 2 and min(t) > 0:
        sig["span_ev"] = math.log2(max(t) / min(t))
    for k, v in (("extent_r", r.get("corona_extent_R")),
                 ("neg_frac", r.get("merge_negative_frac")),
                 ("noise_scale", (r.get("merge_noise") or {}).get("scale"))):
        if isinstance(v, (int, float)) and math.isfinite(v):
            sig[k] = float(v)
    why = []
    if sig.get("span_ev", 99) < SPAN_EV:
        why.append(f"the bracket spans only {sig['span_ev']:.1f} EV")
    if sig.get("extent_r", 99) < EXTENT_R:
        why.append(f"the corona is traced only to {sig['extent_r']:.1f} lunar radii")
    if sig.get("neg_frac", 0) > NEG_FRAC:
        why.append(f"{100 * sig['neg_frac']:.0f} % of the sky came out below black")
    if sig.get("noise_scale", 0) > NOISE_SCALE:
        why.append(f"the noise is {sig['noise_scale']:.1f}x that of a single exposure")
    out.update(thin=len(why) >= MIN_SIGNALS, why=why, signals=sig)
    return out


def black_point(layers, params, extent_r=None):
    """The lowest bgBlack at which the sky's median reaches SKY_TARGET, on
    the preview. The sky is everything beyond 2.4 R (or 0.5 R past the traced
    corona when that leaves too little)."""
    from .render import render
    cy, cx, R = layers.geometry(layers.prev_decim)
    H, W = layers.prev["bg"].shape
    yy, xx = np.mgrid[0:H, 0:W]
    rn = np.hypot(yy - cy, xx - cx) / R
    sky = rn > 2.4
    if sky.sum() < 2000 and extent_r:
        sky = rn > float(extent_r) + 0.5
    if sky.sum() < 500:
        return None

    def med(b):
        q = dict(params); q["bgBlack"] = float(b)
        img = render(layers, q, preview=True)
        return float(np.median(img[sky].mean(axis=1)))

    lo, hi = 0.0, 0.4
    if med(lo) <= SKY_TARGET:
        return lo
    if med(hi) > SKY_TARGET:
        return hi
    for _ in range(10):
        mid = 0.5 * (lo + hi)
        if med(mid) > SKY_TARGET:
            lo = mid
        else:
            hi = mid
    return round(hi, 3)


def for_layers(layers):
    """What the page needs: the verdict, and for a thin set the look with its
    measured black point. Remembered per work directory and report."""
    wd = getattr(layers, "wd", None)
    if not wd:
        return {"thin": False}
    try:
        key = (wd, os.path.getmtime(os.path.join(wd, "report.json")))
    except Exception:
        key = (wd, None)
    if key in _MEMO:
        return _MEMO[key]
    d = detect(wd)
    if d["thin"]:
        look = dict(_M4); look.update(THIN_CHANGES)
        try:
            b = black_point(layers, {k: v for k, v in look.items() if k != "warm"},
                            d["signals"].get("extent_r"))
        except Exception:
            b = None
        if b is not None:
            look["bgBlack"] = b
        d["look"] = look
    _MEMO[key] = d
    return d
