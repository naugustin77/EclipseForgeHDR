"""NOTES ON A RUN (2026-11-14).

A run that finishes can still carry warnings that matter for the picture:
the log says them once, in orange, among hundreds of lines, and a user who
sees the run carry on reads them as handled. These are the warnings that can
show in the result. After a run they are collected from the log, worded for
the person looking at the picture (what may show, what to try), saved in the
work directory and shown above the result's controls.

The warnings that stop a run are not here: those raise pipeline.RunStopped and
the run ends with the message. The warnings not listed here are fallbacks that
do not change the picture in a way the user can act on; they stay in the log.
"""
import json
import os
import re

NOTES_FILE = "run_notes.json"

# (pattern on the log line, title, what it means for the picture and what to try)
RULES = [
    (r"frames move up to \d+ px, beyond the",
     "Frames moved far within an exposure",
     "Some frames of one exposure moved further than the alignment searches, so they may be "
     "misplaced. Look for a doubled corona or Moon edge; if you see one, check those frames in "
     "1 Images (Frames)."),
    (r"tiers disagree by .* at the limb",
     "Rings at the Moon's edge possible",
     "The run measured the exposures disagreeing at the Moon's edge. Look at the limb and the bright "
     "prominences: if concentric rings show there, Merge weight \"Blended edge\" (new stack) can hide "
     "them, but it can also darken the corona just outside the limb and leave a magenta ring round the "
     "disc. Compare the two results and keep the better one. If no rings show, change nothing.",
     {"id": "feather", "value": "plain", "label": "Try Blended edge (new stack)"}),
    (r"limb fit rms .* limb smeared",
     "Soft Moon edge",
     "The Moon's edge came out smeared: the frames were not placed exactly on each other. Expect a "
     "soft or slightly doubled limb."),
    (r"tier lunar limbs spread .* cross-tier alignment failure",
     "Moon edge not in one place",
     "The exposures put the Moon's edge in slightly different places. The limb may look doubled; "
     "check frames that do not belong in 1 Images (Frames)."),
    (r"no cross-bracket link held",
     "Brackets not tied on the corona",
     "The brackets could not be tied together on the corona. Look for a doubled corona; if you see "
     "one, switch the bracket that does not fit off in 1 Images (Frames)."),
    (r"no usable alignment link, shift taken from",
     "An exposure placed by its neighbour",
     "One exposure could not be aligned on its own and was placed like its neighbour. Look for a "
     "doubled edge in the bright inner corona."),
    (r"frames not scene-linear",
     "Frames not linear",
     "The frames do not respond in proportion to the light, so brightness steps between exposures "
     "may show. Use the camera's raw files if you can."),
    (r"links lean the same way",
     "Exposure ratios drift",
     "The exposures' brightness ratios all drift the same way, which points to a black-level "
     "offset. A gentle brightness step between exposures can show."),
    (r"tiers disagree photometrically beyond exposure ratio",
     "Exposure data look wrong",
     "Some exposures are more than 4x brighter or darker than their exposure times say: the "
     "exposure data are probably wrong (labels, missing data). Try Tiers \"By measured "
     "brightness\" in 2 Processing."),
    (r"bracket spans only",
     "Short bracket",
     "The bracket covers fewer stops than totality needs: the inner corona may be clipped or the "
     "outer corona noisy."),
    (r"outside the predicted totality",
     "Camera clock or site off",
     "Frames fall outside the predicted totality: the camera clock or the observing site is off. "
     "Check the site in 2 Processing; contact frames can set the clock."),
    (r"every frame looks like partial phase",
     "Partial-phase frames",
     "An exposure looks like the partial phase in every frame; the Moon's position and size from it "
     "are unreliable."),
    (r"correction disabled: master",
     "Calibration frames not applied",
     "Flats or darks have another size than the light frames and were not applied."),
    # from 0.5 % on: the log line also appears for a few hundredths of a percent
    # (printed as 0.0%), which shows nowhere in the picture
    (r"(?<![\d.])(?:[1-9]\d*\.\d|0\.[5-9])% of merged samples negative",
     "Sky below black",
     "Part of the sky came out below black and was set to zero; the sky may look flat or blotchy."),
    (r"per-tier radial check: .* departs up to",
     "An exposure's brightness differs",
     "One exposure's brightness profile differs from the others near the Moon; a faint ring may show "
     "where it takes over."),
    (r"prominence layer not built",
     "Prominence layer missing",
     "The prominence layer could not be built: the prominences stay in the corona and their own "
     "controls do nothing."),
    (r"partial convolution: masks not built",
     "Partial convolution missing",
     "The partial-convolution layer could not be built; its controls do nothing in this run."),
    (r"Moon layer not built",
     "Earthshine layer missing",
     "The Moon layer could not be built: the Moon stays smooth and the Earthshine control does "
     "nothing."),
]
# a rule may carry a fourth item: the setting to try ({id, value, label}),
# offered as a button on the note -- only after a run, only when it fired
_RX = [(re.compile(r[0], re.I), r[1], r[2], (r[3] if len(r) > 3 else None)) for r in RULES]


def collect(lines, levels):
    """Notes from a run's log: [{title, text, line}], one per rule."""
    out, seen = [], set()
    for i, ln in enumerate(lines or []):
        lv = levels[i] if levels and i < len(levels) else None
        if lv not in (None, "warn", "fail", "odd"):
            continue
        for rx, title, text, act in _RX:
            if title in seen:
                continue
            if rx.search(str(ln)):
                seen.add(title)
                n = {"title": title, "text": text, "line": str(ln)[:300]}
                if act:
                    n["act"] = dict(act)
                out.append(n)
    return out


def save(wd, lines, levels):
    """Write the notes of this run. A run that did not stack again (layers or
    the Moon layer rebuilt from the cache) keeps the stacking run's notes and
    adds its own."""
    notes = collect(lines, levels)
    if not any("merged tier" in str(l) for l in (lines or [])):
        old = load(wd)
        have = {n["title"] for n in notes}
        notes = [n for n in old if n["title"] not in have] + notes
    try:
        json.dump({"notes": notes}, open(os.path.join(wd, NOTES_FILE), "w"), indent=1)
    except Exception:
        pass
    return notes


def load(wd):
    try:
        return (json.load(open(os.path.join(wd, NOTES_FILE))) or {}).get("notes") or []
    except Exception:
        return []
