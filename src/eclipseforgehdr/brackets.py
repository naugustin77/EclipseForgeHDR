"""Repeated brackets: which frames were shot together (2026-10-14).

A capture program -- SetNC, Umbra, Solar Eclipse Workbench -- or a patient
photographer runs the same exposure ladder again and again through totality.
Typical sets run from a few brackets of a dozen or more exposures in about
a minute and a half to fifteen brackets over about 4 minutes.
EclipseForgeHDR stacks by shutter speed, so until now a 1/8 s frame from the
first bracket and one from the last were averaged as one measurement, although
the sky, the haze and the Moon's position had all changed in between.

This module only FINDS the brackets, from the data the frames carry: their
time stamps and their exposure times, in shooting order. It changes nothing
downstream by itself; the pipeline records the result and reports it, and the
per-bracket steps (photometry, Moon, drift) build on it.

THE RULE. Walk the frames in shooting order. A bracket ends when an exposure
time comes round again that the bracket already holds, or when the gap to the
next frame is far longer than the usual interval. That is order-free -- it
works for ladders shot short-to-long, long-to-short or back and forth -- and
needs no knowledge of the ladder in advance.

A set shot TIER BY TIER (all 1/500 s frames, then all 1/250 s, ...) gives
"brackets" of one exposure each under that rule. That is not a repeated
bracket, so it is recognised and reported as a single pass: repeated brackets
need at least two of them, holding a median of at least three distinct
exposures each.
"""
import os
import re


def _natural(name):
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r"(\d+)", name)]


def _key(sec):
    """Exposure times as written in EXIF (1/125 = 0.008) compared to 0.1 %."""
    return round(float(sec), 6) if sec < 1e-3 else float(f"{float(sec):.4g}")


def find_brackets(paths, sec, epoch, gap_factor=5.0, gap_min_s=10.0):
    """Group frames into brackets.

    paths  : frame paths
    sec    : {path: exposure seconds}
    epoch  : {path: seconds (any origin) or None}
    Returns a dict:
      repeated   True when the set is repeated brackets
      n          number of brackets (1 when not repeated)
      of         {path: bracket index}, in shooting order from 0
      brackets   [{first, last, n_frames, n_exposures, t0, t1}] (t relative to
                 the first frame, None without time stamps)
      span_s     time from the first to the last frame, or None
      reason     one line for the log when not repeated
    """
    have_t = all(epoch.get(p) is not None for p in paths)
    order = sorted(paths, key=lambda p: ((epoch.get(p) if have_t else 0.0),
                                         _natural(os.path.basename(p))))
    if not order:
        return {"repeated": False, "n": 0, "of": {}, "brackets": [],
                "span_s": None, "reason": "no frames"}
    gaps = []
    if have_t:
        gaps = [epoch[b] - epoch[a] for a, b in zip(order, order[1:])]
    med_gap = sorted(gaps)[len(gaps) // 2] if gaps else 0.0
    gap_lim = max(gap_min_s, gap_factor * med_gap)

    groups, cur, seen = [], [], set()
    for i, p in enumerate(order):
        k = _key(sec[p])
        split = k in seen
        if not split and have_t and i > 0 and gaps[i - 1] > gap_lim:
            split = True
        if split and cur:
            groups.append(cur)
            cur, seen = [], set()
        cur.append(p)
        seen.add(k)
    if cur:
        groups.append(cur)

    distinct = sorted(len({_key(sec[p]) for p in g}) for g in groups)
    med_distinct = distinct[len(distinct) // 2]
    t_first = epoch[order[0]] if have_t else None
    span = (epoch[order[-1]] - t_first) if have_t else None
    # at least two REAL brackets: runs that hold at least three exposures and
    # at least half of the fullest one. A single pass plus one re-shot frame,
    # or a ladder shot tier by tier, is not a repeated bracket.
    full = distinct[-1]
    n_real = sum(1 for d in distinct if d >= max(3, 0.5 * full))
    repeated = n_real >= 2 and med_distinct >= 3
    if not repeated:
        why = ("one pass through the exposures" if len(groups) < 2 else
               f"shot exposure by exposure ({med_distinct} distinct per run)"
               if med_distinct < 3 else
               f"one full pass and {len(groups) - 1} short run(s)")
        return {"repeated": False, "n": 1, "of": {p: 0 for p in order},
                "brackets": [{"first": os.path.basename(order[0]),
                              "last": os.path.basename(order[-1]),
                              "n_frames": len(order),
                              "n_exposures": len({_key(sec[p]) for p in order}),
                              "t0": 0.0 if have_t else None,
                              "t1": span}],
                "span_s": span, "reason": why}
    of, info = {}, []
    for b, g in enumerate(groups):
        for p in g:
            of[p] = b
        info.append({"first": os.path.basename(g[0]), "last": os.path.basename(g[-1]),
                     "n_frames": len(g),
                     "n_exposures": len({_key(sec[p]) for p in g}),
                     "t0": (epoch[g[0]] - t_first) if have_t else None,
                     "t1": (epoch[g[-1]] - t_first) if have_t else None})
    return {"repeated": True, "n": len(groups), "of": of, "brackets": info,
            "span_s": span, "reason": ""}


def describe(res):
    """One log line."""
    if not res["n"]:
        return "brackets: no frames"
    span = res["span_s"]
    st = (f" over {int(span // 60)} min {int(round(span % 60))} s" if span and span >= 60
          else f" over {span:.0f} s" if span is not None else "")
    if not res["repeated"]:
        return f"one bracket{st} ({res['reason']})"
    ne = sorted(b["n_exposures"] for b in res["brackets"])
    full = ne[-1]
    short = sum(1 for x in ne if x < full)
    return (f"repeated brackets: {res['n']} brackets of up to {full} exposures{st}"
            + (f", {short} incomplete" if short else ""))
