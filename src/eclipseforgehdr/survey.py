"""How a folder was shot: the survey behind tab 1 (2026-11-13).

The stack makes assumptions about a shoot -- the framing stays where it was,
every exposure has several frames, all frames are from totality, the corona
fits the frame -- and until now it found out the hard way, minutes into a
run, or not at all. This module looks first. It reads what the frames carry
(exposure, time, ISO, aperture, focal length) and the camera's own small
preview inside each raw file, decodes no raw data, and reports in seconds:

  STAGE 1, from the metadata: the exposure ladder and its steps, repeated or
  complementary brackets, frames per exposure, pauses in the shooting,
  settings that change within the set.

  STAGE 2, from the previews: the Moon's size and its distance to the frame
  edges, whether the framing moved between brackets (measured on the corona,
  both Moons masked), frames with a bright spot at the limb (diamond ring,
  beads), and which exposures are clipped how far out.

Everything is a MEASUREMENT with its unit, in fractions of the Moon's radius
R so it holds for any focal length and sensor; pixels at full size are added
when the raw size can be read. The previews are the camera's 8-bit JPEGs:
they clip earlier than the raw data, so clipping radii are lower limits, and
they are said to be.

WHAT IT IS FOR. The findings show as rows in tab 1 ("How this was shot").
The same record (survey.json in the work folder) carries a proposed role per
frame and proposed stacking settings with their reasons; nothing applies them
yet. `python -m eclipseforgehdr.survey FOLDER` prints the report as text.

Only numpy and PIL are used here, so the survey also runs where the raw
decoder is missing.
"""
from __future__ import annotations
import io
import json
import math
import os
import threading

import numpy as np

SURVEY_BUILD = 9

# The Moon against the Sun during totality, in lunar radii per second: about
# 0.3-0.55 arcsec/s at the observer against a radius of 930-1000 arcsec.
LUNAR_RATE_R = 0.0006
SMALL = 720                 # long side of the analysed preview, px
SMALL_MIN_R = 45.0          # a Moon smaller than this on it: the preview is read at twice the size
GAP_MIN_S = 10.0            # a pause in the shooting (as brackets.find_brackets)
STEP_WIDE_EV = 4.5          # an exposure step from here on is named in the ladder row
MOVE_WARN_R = 0.020         # corona shift between two groups that is a moved framing
DRIFT_INFO_R = 0.020        # total corona drift over the set worth a row
TIGHT_R = 3.0               # nearest frame edge closer than this: the corona runs off
SPOT_IN = 0.0008            # white-clipped area INSIDE the limb (in Moon discs) that marks a bead or a
                            # diamond ring by itself. Clipping outside the limb does not: on two test sets
                            # a large prominence and the brighter side of the inner corona clip one-sidedly
DISC_DEV_GROUP = 0.025      # the dark disc's radius against the other frames at the exposure
DISC_DEV_SET = 0.06         # ... and against the whole set (a long exposure's disc blooms 4 % smaller)
SPOT_OFF = 0.45             # ... when its centroid sits this far (R) from the Moon's centre
SPOT_EXCESS = 0.0008        # ... or this much more than the other frames at the exposure show
SPOT_PEER_S = 3.0           # a shorter exposure this close in time must show the spot too
SPOT_NEAR_S = 1.5           # frames this close in time to a flagged one go with it
SPOT_END_S = 8.0            # ... and all frames beyond it, when it is this close to an end of the set
MAX_PAIRS = 24              # corona links measured between groups
LINK_AGREE_PX = 0.75        # a corona link is trusted when its two band-passes agree to this (preview px)
LINK_PEAK = 0.4             # ... its correlation peak is at least this
                            # ... and the Moon's motion against the corona that it implies is possible


def _ex_txt(t):
    return f"{t:g} s" if t >= 0.4 else f"1/{round(1 / t):g} s"


def _key(sec):
    from .brackets import _key as k
    return k(sec)


# --------------------------------------------------------------- metadata

def frame_meta(path):
    """What a frame says about itself: dict(sec, iso, ts, subsec, fnum, focal,
    model), any of them None. Never raises."""
    m = dict.fromkeys(("sec", "iso", "ts", "subsec", "fnum", "focal", "model"))

    def rat(v):
        try:
            if isinstance(v, (list, tuple)) and len(v) == 2 and not isinstance(v[0], (list, tuple)):
                return float(v[0]) / float(v[1]) if v[1] else None
            if isinstance(v, (list, tuple)):
                return rat(v[0])
            return float(v)
        except Exception:
            return None
    ext = os.path.splitext(path)[1].lower()
    try:
        from .fits import is_fits
        fits = is_fits(path)
    except Exception:
        fits = False
    try:
        if ext == ".cr3":
            from .raw import _cr3_tags
            t = _cr3_tags(path)
            c1, c2 = t.get("CMT1") or {}, t.get("CMT2") or {}
            iso = c2.get(0x8827) or c2.get(0x8833)
            if isinstance(iso, (list, tuple)):
                iso = iso[0] if iso else None
            m.update(sec=rat(c2.get(0x829A)), iso=int(iso) if iso else None,
                     ts=str(c2.get(0x9003) or c1.get(0x0132) or "") or None,
                     subsec=str(c2.get(0x9291) or "").strip() or None,
                     fnum=rat(c2.get(0x829D)), focal=rat(c2.get(0x920A)),
                     model=str(c1.get(0x0110) or "").strip() or None)
        elif not fits:
            try:
                import exifread
            except ImportError:
                exifread = None
            if exifread is not None:
                with open(path, "rb") as f:
                    tg = exifread.process_file(f, details=False)

                def g(*names):
                    for n in names:
                        if n in tg:
                            return tg[n]
                    return None

                def r1(tag):
                    if tag is None or not getattr(tag, "values", None):
                        return None
                    v = tag.values[0]
                    try:
                        return float(v.num) / float(v.den) if v.den else None
                    except AttributeError:
                        return float(v)
                iso = g("EXIF ISOSpeedRatings", "Image ISOSpeedRatings")
                m.update(sec=r1(g("EXIF ExposureTime", "Image ExposureTime")),
                         iso=int(iso.values[0]) if iso is not None and iso.values else None,
                         ts=str(g("EXIF DateTimeOriginal", "Image DateTime") or "") or None,
                         subsec=str(g("EXIF SubSecTimeOriginal", "EXIF SubSecTime") or "").strip() or None,
                         fnum=r1(g("EXIF FNumber")), focal=r1(g("EXIF FocalLength")),
                         model=str(g("Image Model") or "").strip() or None)
            else:
                from PIL import Image
                with Image.open(path) as im:
                    ex = im.getexif()
                    sub = ex.get_ifd(0x8769)
                    m.update(sec=rat(sub.get(0x829A)), iso=(int(sub[0x8827]) if sub.get(0x8827) else None),
                             ts=str(sub.get(0x9003) or ex.get(0x0132) or "") or None,
                             subsec=str(sub.get(0x9291) or "").strip() or None,
                             fnum=rat(sub.get(0x829D)), focal=rat(sub.get(0x920A)),
                             model=str(ex.get(0x0110) or "").strip() or None)
    except Exception:
        pass
    if m["sec"] is None or m["ts"] is None:
        try:
            from .raw import read_exif
            sec, iso, ts = read_exif(path)
            m["sec"] = m["sec"] if m["sec"] is not None else (float(sec) if sec else None)
            m["iso"] = m["iso"] if m["iso"] is not None else (int(iso) if iso else None)
            m["ts"] = m["ts"] or (str(ts) if ts else None)
        except Exception:
            pass
    return m


def _sig(path):
    st = os.stat(path)
    return f"{st.st_size}:{int(st.st_mtime)}"


def _load_cache(wd):
    try:
        c = json.load(open(os.path.join(wd, "survey.json")))
        if c.get("build") == SURVEY_BUILD:
            return c
    except Exception:
        pass
    return {"build": SURVEY_BUILD, "meta": {}, "prev": {}}


def _save_cache(wd, cache):
    try:
        os.makedirs(wd, exist_ok=True)
        tmp = os.path.join(wd, "survey.json.tmp")
        json.dump(cache, open(tmp, "w"))
        os.replace(tmp, os.path.join(wd, "survey.json"))
    except Exception:
        pass


def _raw_width(path):
    """Full width of the raw image in px, or None (needs the raw decoder)."""
    try:
        import rawpy
        with rawpy.imread(path) as r:
            return int(r.sizes.width)
    except Exception:
        return None


# ------------------------------------------------------------- stage 1

def exif_stage(folder, wd, cache=None):
    """The shoot as its metadata tells it. Returns the record (see module
    docstring); record["items"] are the rows for tab 1."""
    from .raw import list_raws
    from .timeline import epoch_of, _exif_table
    from .brackets import find_brackets, _natural
    cache = cache if cache is not None else _load_cache(wd)
    raws = list_raws(folder)
    base = {}
    try:
        base = _exif_table(folder, raws, wd)          # sec, ts: shared with the timeline
    except Exception:
        pass
    frames, unread = [], 0
    for p in raws:
        n = os.path.basename(p)
        try:
            sig = _sig(p)
        except OSError:
            continue
        c = cache["meta"].get(n)
        if not c or c.get("sig") != sig:
            c = dict(frame_meta(p), sig=sig)
            cache["meta"][n] = c
        sec, ts = c.get("sec"), c.get("ts")
        b = base.get(n) or (None, None)
        sec = sec if sec else b[0]
        ts = ts or b[1]
        if not sec:
            unread += 1
            continue
        ep = epoch_of(ts) if ts else None
        if ep is not None and c.get("subsec"):
            try:
                ep += float("0." + "".join(ch for ch in str(c["subsec"]) if ch.isdigit()))
            except ValueError:
                pass
        frames.append({"name": n, "path": p, "sec": float(sec), "ep": ep, "iso": c.get("iso"),
                       "fnum": c.get("fnum"), "focal": c.get("focal"), "model": c.get("model")})
    rec = {"n": len(frames), "unread": unread, "items": [], "suggest": {}, "roles": {}}
    try:
        from .fits import is_fits as _isf
        rec["fits"] = any(_isf(f["path"]) for f in frames)
    except Exception:
        rec["fits"] = False
    items = rec["items"]

    def add(level, title, text, act=None):
        items.append({"level": level, "title": title, "text": text, "act": act})
    if len(frames) < 3:
        rec["frames"] = frames
        return rec
    paths = [f["path"] for f in frames]
    sec = {f["path"]: f["sec"] for f in frames}
    ep = {f["path"]: f["ep"] for f in frames}
    timed = all(v is not None for v in ep.values())
    br = find_brackets(paths, sec, ep)
    order = sorted(frames, key=lambda f: ((f["ep"] if timed else 0.0), _natural(f["name"])))
    t0 = order[0]["ep"] if timed else None
    for i, f in enumerate(order):
        f["t"] = (f["ep"] - t0) if timed else None
        f["bracket"] = int(br["of"].get(f["path"], 0))
        f["order"] = i
    frames = order
    rec["frames"] = frames
    rec["timed"] = timed
    rec["span_s"] = br.get("span_s")

    # the ladder and its steps
    ladder = sorted({_key(f["sec"]) for f in frames})
    per = {k: sum(1 for f in frames if _key(f["sec"]) == k) for k in ladder}
    steps = [math.log2(b / a) for a, b in zip(ladder, ladder[1:])]
    rec.update(ladder=ladder, per_exposure=per, ev_steps=[round(s, 2) for s in steps],
               range_ev=round(math.log2(ladder[-1] / ladder[0]), 1) if len(ladder) > 1 else 0.0)
    # how it was shot
    nb = int(br["n"]) if br["repeated"] else 1
    style, sets = "single", []
    if br["repeated"]:
        sets = [set() for _ in range(nb)]
        for f in frames:
            sets[f["bracket"]].add(_key(f["sec"]))
        jac = [len(a & b) / max(len(a | b), 1) for a, b in zip(sets, sets[1:])]
        style = "complementary" if (jac and sorted(jac)[len(jac) // 2] < 0.34) else "repeated"
    elif "exposure by exposure" in (br.get("reason") or ""):
        style = "by_exposure"
    rec.update(style=style, n_brackets=nb)
    cam = next((f["model"] for f in frames if f.get("model")), None)
    foc = sorted({round(f["focal"]) for f in frames if f.get("focal")})
    fno = sorted({round(f["fnum"], 1) for f in frames if f.get("fnum")})
    iso = sorted({int(f["iso"]) for f in frames if f.get("iso")})
    rec.update(camera=cam, focal=foc, fnum=fno, iso=iso)

    span = br.get("span_s")
    sp = ("" if span is None or span < 1 else (f" over {int(span // 60)} min {int(round(span % 60))} s" if span >= 60
                                                else f" over {span:.0f} s"))
    if style == "repeated":
        full = max(len(s) for s in sets)
        short = sum(1 for s in sets if len(s) < full)
        per_b = None
        if timed and nb > 1:
            starts = [min(f["t"] for f in frames if f["bracket"] == b) for b in range(nb)]
            d = sorted(b - a for a, b in zip(starts, starts[1:]))
            per_b = d[len(d) // 2]
        txt = (f"the same bracket of {full} exposures repeated {nb} times{sp}"
               + (f", one every {per_b:.0f} s" if per_b and per_b >= 1.5 else "")
               + (f"; {short} incomplete" if short else ""))
    elif style == "complementary":
        txt = (f"{nb} brackets{sp} that hold different exposures ("
               + "; ".join(f"{len(s)} from {_ex_txt(min(s))} to {_ex_txt(max(s))}" for s in sets[:4])
               + ("; …" if nb > 4 else "") + "): each exposure exists in one bracket only, "
               "so the brackets must be aligned to each other exactly")
    elif style == "by_exposure":
        txt = f"shot exposure by exposure{sp}: all frames of one exposure, then the next"
    else:
        txt = f"one pass through the exposures{sp}" + (
            f" ({br['reason']})" if br.get("reason") and "one pass" not in br["reason"] else "")
    add("ok", "Shooting pattern", txt)

    def _ev(s):
        return f"{s:.0f}" if abs(s - round(s)) < 0.15 else f"{s:.1f}"
    if len(ladder) > 1:
        wide = [(a, b, s) for a, b, s in zip(ladder, ladder[1:], steps) if s >= STEP_WIDE_EV]
        txt = (f"{len(ladder)} exposures, {_ex_txt(ladder[0])} to {_ex_txt(ladder[-1])} "
               f"({rec['range_ev']:.0f} EV); steps "
               + (f"of {_ev(steps[0])}" if len({_ev(s) for s in steps}) == 1 else
                  ", ".join(_ev(s) for s in steps) if len(steps) <= 8 else
                  f"{_ev(min(steps))} to {_ev(max(steps))}") + " EV")
        if wide:
            a, b, s = max(wide, key=lambda w: w[2])
            rec["wide_step"] = {"short": a, "long": b, "ev": round(s, 1)}
            txt += (f". Wide step {_ex_txt(a)} → {_ex_txt(b)}: where the longer one clips, the shorter "
                    f"one takes over at 1/{2 ** s:.0f} of its full scale (a noisier hand-over)")
        add("info" if wide else "ok", "Exposure ladder", txt)
    else:
        add("warn", "Exposure ladder", f"one exposure only ({_ex_txt(ladder[0])}): nothing to merge into an HDR")
    counts = sorted(per.values())
    med = counts[len(counts) // 2]
    if counts[0] == counts[-1]:
        ctxt = f"{counts[0]} frame{'s' if counts[0] != 1 else ''} per exposure"
    else:
        ctxt = f"{counts[0]} to {counts[-1]} frames per exposure (median {med})"
    rec["frames_per_exposure"] = med
    if med < 3:
        add("info", "Frames per exposure", ctxt + ": too few for a per-pixel rejection "
            "(Tier combine: no rejection is the setting for this)", {"tab": 2})
        rec["suggest"]["stack_combine"] = {"value": "mean", "why": ctxt}
    else:
        add("ok", "Frames per exposure", ctxt)
    # pauses
    if timed and len(frames) > 2:
        gaps = [(b["t"] - a["t"], a, b) for a, b in zip(frames, frames[1:])]
        gs = sorted(g[0] for g in gaps)
        lim = max(GAP_MIN_S, 5.0 * gs[len(gs) // 2])
        big = sorted((g for g in gaps if g[0] > lim), key=lambda g: -g[0])
        rec["pauses"] = [{"s": round(g[0], 1), "after": g[1]["name"], "before": g[2]["name"]} for g in big[:8]]
        if big:
            add("info", "Pauses", (f"{len(big)} pauses in the shooting, the longest {big[0][0]:.0f} s"
                                   if len(big) > 1 else f"one pause of {big[0][0]:.0f} s in the shooting")
                + f" (after {big[0][1]['name']})")
    elif not timed:
        add("warn", "Time stamps", "some frames carry no time stamp: brackets and pauses are read "
            "from the file order only")
    # settings that changed
    def _chg(title, vals, fmt, why):
        if len(vals) > 1:
            cnt = {v: 0 for v in vals}
            key = {"ISO": lambda f: int(f["iso"]) if f.get("iso") else None,
                   "Aperture": lambda f: round(f["fnum"], 1) if f.get("fnum") else None,
                   "Focal length": lambda f: round(f["focal"]) if f.get("focal") else None}[title]
            for f in frames:
                k = key(f)
                if k in cnt:
                    cnt[k] += 1
            add("warn", title, "changes within the set: "
                + ", ".join(f"{fmt(v)} ({cnt[v]} frames)" for v in vals) + ". " + why,
                {"field": "fpanel"})
            rec.setdefault("mixed", {})[title] = {str(v): cnt[v] for v in vals}
            if title in ("ISO", "Aperture"):
                rec["suggest"]["tier_mode"] = {"value": "measured", "why": f"{title} changes within the set"}
    _chg("ISO", iso, lambda v: f"ISO {v}",
         "Frames are grouped by exposure time, so both land in one group. Tiers: By measured brightness "
         "(tab 2) groups them by what they recorded; or switch one ISO off in Frames")
    _chg("Aperture", fno, lambda v: f"f/{v:g}",
         "Frames are grouped by exposure time, so both land in one group. Tiers: By measured brightness "
         "(tab 2) groups them by what they recorded")
    _chg("Focal length", foc, lambda v: f"{v:g} mm",
         "Frames at another image scale cannot be stacked together; switch them off in Frames")
    if cam or foc or fno or iso:
        bits = [cam] if cam else []
        if len(foc) == 1:
            bits.append(f"{foc[0]:g} mm")
        if len(fno) == 1:
            bits.append(f"f/{fno[0]:g}")
        if len(iso) == 1:
            bits.append(f"ISO {iso[0]}")
        if bits:
            add("ok", "Camera", ", ".join(bits))
    return rec


# ----------------------------------------------------------- image tools

def _decode_small(j, small):
    from PIL import Image
    im = Image.open(io.BytesIO(j))
    im.draft("RGB", (small, small))
    im = im.convert("RGB")
    k = max(1, int(round(max(im.size) / small)))
    if k > 1:
        im = im.reduce(k)
    return np.asarray(im, np.uint8)


def _small_rgb(path, small=SMALL):
    """The camera's preview of a raw frame as a small RGB array (long side
    about `small` px, the sensor's orientation), or None."""
    from PIL import Image
    from .frameselect import _embedded_jpeg
    try:
        j = _embedded_jpeg(path)
        if j is not None:
            return _decode_small(j, small)
    except Exception:
        pass
    # the stream offered does not decode (without the raw decoder the search
    # for it can land on raw data): every JPEG stream in the file, largest first
    try:
        b = open(path, "rb").read()
    except OSError:
        return None
    cand, i = [], 0
    while True:
        i = b.find(b"\xff\xd8\xff", i)
        if i < 0:
            break
        try:
            im = Image.open(io.BytesIO(b[i:i + 65536]))
            if im.format == "JPEG":
                cand.append((im.size[0] * im.size[1], i))
        except Exception:
            pass
        i += 3
    for _, i in sorted(cand, reverse=True)[:4]:
        try:
            return _decode_small(b[i:], small)
        except Exception:
            continue
    return None


def _box(a, k):
    h, w = (a.shape[0] // k) * k, (a.shape[1] // k) * k
    return a[:h, :w].reshape(h // k, k, w // k, k).mean((1, 3))


def _bilinear(a, Y, X):
    """a sampled at (Y, X); NaN outside the frame."""
    H, W = a.shape
    ok = (Y >= 0) & (Y <= H - 1.001) & (X >= 0) & (X <= W - 1.001)
    y = np.clip(Y, 0, H - 1.001)
    x = np.clip(X, 0, W - 1.001)
    y0 = y.astype(np.int32)
    x0 = x.astype(np.int32)
    fy, fx = y - y0, x - x0
    v = (a[y0, x0] * (1 - fy) * (1 - fx) + a[y0 + 1, x0] * fy * (1 - fx)
         + a[y0, x0 + 1] * (1 - fy) * fx + a[y0 + 1, x0 + 1] * fy * fx)
    return np.where(ok, v, np.nan)


def _fit_circle(p):
    """Least-squares circle through points (y, x) with three rounds of
    outlier rejection: (cy, cx, R, rms, n) or None."""
    p = np.asarray(p, float)
    y0 = x0 = r0 = None
    for _ in range(4):
        if len(p) < 12:
            return None
        A = np.c_[2 * p[:, 1], 2 * p[:, 0], np.ones(len(p))]
        b = (p ** 2).sum(1)
        c = np.linalg.lstsq(A, b, rcond=None)[0]
        x0, y0 = c[0], c[1]
        r0 = math.sqrt(max(c[2] + x0 * x0 + y0 * y0, 1e-9))
        res = np.hypot(p[:, 0] - y0, p[:, 1] - x0) - r0
        keep = np.abs(res) < max(3 * 1.4826 * np.median(np.abs(res)), 0.6)
        if keep.all():
            break
        p = p[keep]
    res = np.hypot(p[:, 0] - y0, p[:, 1] - x0) - r0
    return float(y0), float(x0), float(r0), float(np.sqrt(np.mean(res ** 2))), int(len(p))


def find_moon(L):
    """The lunar disc on a small luminance image (0-255): (cy, cx, R, rms, n)
    in its pixels, or None.

    Coarse: at about 160 px, the place and radius where a ring just outside
    is brightest against the disc inside it (a matched pair of disc and ring,
    by FFT; a diamond ring adds to a part of the ring only and does not move
    it). Fine: 360 rays, the steepest dark-to-bright step on each, a circle
    through them with outliers rejected, three rounds."""
    k = max(1, int(round(max(L.shape) / 160)))
    s = _box(L, k).astype(np.float64)
    h, w = s.shape
    F = np.fft.rfft2(s)
    yy = np.minimum(np.arange(h), h - np.arange(h))[:, None]
    xx = np.minimum(np.arange(w), w - np.arange(w))[None, :]
    d = np.hypot(yy, xx)
    best = None
    for r in np.geomspace(3.0, 0.48 * min(h, w), 30):
        disc = (d <= 0.85 * r).astype(np.float64)
        ann = ((d >= 1.03 * r) & (d <= 1.3 * r)).astype(np.float64)
        if disc.sum() < 5 or ann.sum() < 5:
            continue
        A = np.fft.irfft2(F * np.conj(np.fft.rfft2(ann / ann.sum())), s=(h, w))
        D = np.fft.irfft2(F * np.conj(np.fft.rfft2(disc / disc.sum())), s=(h, w))
        sc = (A - D) / (A + D + 6.0)
        m = int(math.ceil(0.9 * r))
        if h - 2 * m < 1 or w - 2 * m < 1:
            continue
        sub = sc[m:h - m, m:w - m]
        i = int(np.argmax(sub))
        v = float(sub.flat[i])
        if best is None or v > best[0]:
            y, x = divmod(i, sub.shape[1])
            best = (v, (y + m + 0.5) * k - 0.5, (x + m + 0.5) * k - 0.5, r * k)
    if best is None or best[0] < 0.06:
        return None
    _, cy, cx, R = best
    Lf = L.astype(np.float64)
    th = np.linspace(0, 2 * np.pi, 360, endpoint=False)
    out = None
    for it in range(3):
        lo, hi = (0.72, 1.32) if it == 0 else (0.88, 1.12)
        rr = np.arange(lo * R, hi * R, 0.5)
        if rr.size < 8:
            return None
        Y = cy + rr[None, :] * np.sin(th)[:, None]
        X = cx + rr[None, :] * np.cos(th)[:, None]
        v = _bilinear(Lf, Y, X)
        vs = (v[:, :-2] + v[:, 1:-1] + v[:, 2:]) / 3.0
        g = vs[:, 2:] - vs[:, :-2]
        g = np.where(np.isfinite(g), g, -1e9)
        i = np.argmax(g, axis=1)
        st = g[np.arange(g.shape[0]), i]
        good = st > max(6.0, 0.25 * float(np.median(st[st > -1e8])) if (st > -1e8).any() else 6.0)
        rl = rr[i + 2]
        pts = np.c_[cy + rl * np.sin(th), cx + rl * np.cos(th)][good]
        fit = _fit_circle(pts)
        if fit is None:
            return out
        cy, cx, R = fit[0], fit[1], fit[2]
        out = fit
    # a real limb: most rays agree on one circle
    if out is None or out[4] < 90 or out[3] > 0.04 * out[2]:
        return None
    return out


def measure_frame(path, small=SMALL):
    """What one frame's preview shows: the Moon, the clipping, the level of
    the outer field. Returns a small dict, or {"err": reason}."""
    rgb = _small_rgb(path, small)
    if rgb is None:
        return {"err": "no preview in this file"}
    H, W = rgb.shape[:2]
    L = rgb.astype(np.float32).mean(-1)
    fit = find_moon(L)
    if fit is None:
        top = float(np.percentile(L, 99.9))
        if top < 40.0:
            return {"err": "preview too dark to find the Moon", "w": W, "h": H}
        if float(np.percentile(L, 5.0)) > 215.0:
            return {"err": "preview too bright to find the Moon", "w": W, "h": H}
        return {"err": "no lunar disc found", "w": W, "h": H}
    cy, cx, R, rms, n = fit
    yy, xx = np.ogrid[0:H, 0:W]
    r = np.hypot(yy - cy, xx - cx) / R
    white = (rgb.min(-1) >= 245) & (r > 0.97)          # white-clipped: all three channels
    area = float(white.sum()) / (math.pi * R * R)
    # white INSIDE the limb: the corona never clips there, a bead or a diamond
    # ring blooms onto the disc
    inner = (rgb.min(-1) >= 245) & (r < 0.95)
    area_in = float(inner.sum()) / (math.pi * R * R) if int(inner.sum()) >= 3 else 0.0
    off = 0.0
    if white.any():
        ys, xs = np.nonzero(white)
        off = float(np.hypot(ys.mean() - cy, xs.mean() - cx) / R)
    # per ray: clipped at the limb, and how far out
    th = np.linspace(0, 2 * np.pi, 180, endpoint=False)
    rr = np.arange(1.02, 6.0, 0.02)
    Y = cy + R * rr[None, :] * np.sin(th)[:, None]
    X = cx + R * rr[None, :] * np.cos(th)[:, None]
    inside = (Y >= 0) & (Y <= H - 1) & (X >= 0) & (X <= W - 1)
    wv = white[np.clip(Y, 0, H - 1).astype(np.int32), np.clip(X, 0, W - 1).astype(np.int32)] & inside
    at_limb = wv[:, 1:4].any(1)                         # 1.04-1.08 R
    brk = (~wv[:, 1:]).argmax(1)                        # first unclipped sample past the limb
    ext = np.where(at_limb, rr[1 + brk], 0.0)
    ext = np.where(at_limb & wv[:, 1:].all(1), rr[-1], ext)
    lin = (L / 255.0) ** 2.2
    far = (r > 1.6)
    usable = (L > 12) & (L < 245) & (r > 1.25)
    return {"w": W, "h": H, "cy": round(cy, 2), "cx": round(cx, 2), "R": round(R, 2),
            "rms": round(rms, 2), "n": n, "clip_area": round(area, 4), "clip_off": round(off, 3),
            "clip_in": round(area_in, 5),
            "clip_limb": round(float(at_limb.mean()), 3), "clip_r": round(float(np.median(ext)), 2),
            "far": round(float(lin[far].mean()) if far.any() else 0.0, 5),
            "usable": int(usable.sum())}


def _gauss_fft(a, sigma):
    H, W = a.shape
    fy = np.fft.fftfreq(H)[:, None]
    fx = np.fft.rfftfreq(W)[None, :]
    return np.fft.irfft2(np.fft.rfft2(a) * np.exp(-2 * (np.pi * sigma) ** 2 * (fy * fy + fx * fx)), s=(H, W))


def masked_shift(f1, m1, f2, m2, max_shift):
    """Where the picture in frame 1 sits against frame 2 (dy, dx), by a
    normalised cross-correlation that counts only pixels valid in both frames
    at each trial shift (Padfield 2012): no peak from a mask's own edge, none
    from wrap-around. Returns (dy, dx, peak) or None."""
    H, W = f1.shape
    ph, pw = H + int(max_shift) + 2, W + int(max_shift) + 2
    m1 = m1.astype(np.float64)
    m2 = m2.astype(np.float64)
    f1 = f1 * m1
    f2 = f2 * m2
    F = lambda a: np.fft.rfft2(a, s=(ph, pw))
    I = lambda A: np.fft.irfft2(A, s=(ph, pw))
    F1, F2, M1, M2 = F(f1), F(f2), F(m1), F(m2)
    n = I(M1 * np.conj(M2))
    s1 = I(F1 * np.conj(M2))
    s2 = I(M1 * np.conj(F2))
    nn = np.maximum(n, 1.0)
    num = I(F1 * np.conj(F2)) - s1 * s2 / nn
    d1 = I(F(f1 * f1) * np.conj(M2)) - s1 * s1 / nn
    d2 = I(M1 * np.conj(F(f2 * f2))) - s2 * s2 / nn
    ncc = num / np.sqrt(np.maximum(d1, 1e-12) * np.maximum(d2, 1e-12))
    ok = n > 0.5 * float(n.max())
    sy = np.minimum(np.arange(ph), ph - np.arange(ph))[:, None]
    sx = np.minimum(np.arange(pw), pw - np.arange(pw))[None, :]
    ok &= (sy <= max_shift) & (sx <= max_shift)
    ncc = np.where(ok, ncc, -1.0)
    i = int(np.argmax(ncc))
    py, px = divmod(i, pw)
    peak = float(ncc[py, px])
    if peak <= 0:
        return None

    def sub(a, b, c):
        den = a - 2 * b + c
        return 0.0 if abs(den) < 1e-12 else float(np.clip(0.5 * (a - c) / den, -0.5, 0.5))
    dy = py + sub(ncc[(py - 1) % ph, px], peak, ncc[(py + 1) % ph, px])
    dx = px + sub(ncc[py, (px - 1) % pw], peak, ncc[py, (px + 1) % pw])
    if dy > ph / 2:
        dy -= ph
    if dx > pw / 2:
        dx -= pw
    return dy, dx, peak


def corona_shift(path_a, ma, path_b, mb, small=SMALL):
    """Where the corona in frame A sits against frame B, in lunar radii:
    dict(dy, dx, agree, peak) or None. Both frames' own Moons are masked to
    1.25 R, clipped and empty pixels left out; two band-passes must agree."""
    ra, rb = _small_rgb(path_a, small), _small_rgb(path_b, small)
    if ra is None or rb is None or ra.shape != rb.shape:
        return None
    H, W = ra.shape[:2]
    R = 0.5 * (ma["R"] + mb["R"])
    yy, xx = np.ogrid[0:H, 0:W]
    edge = np.zeros((H, W), bool)
    edge[8:-8, 8:-8] = True
    out = []
    for lo, hi in ((1.5, 6.0), (3.0, 12.0)):
        fs, ms = [], []
        for rgb, m in ((ra, ma), (rb, mb)):
            L = rgb.astype(np.float64).mean(-1)
            good = (rgb.max(-1) < 245) & (L > 10) & edge & (np.hypot(yy - m["cy"], xx - m["cx"]) > 1.25 * m["R"])
            if good.sum() < 4000:
                return None
            ms.append(_gauss_fft(good.astype(np.float64), hi) > 0.97)
            fs.append(_gauss_fft(L, lo) - _gauss_fft(L, hi))
        if ms[0].sum() < 3000 or ms[1].sum() < 3000:
            return None
        r = masked_shift(fs[0], ms[0], fs[1], ms[1], max_shift=0.45 * min(H, W))
        if r is None:
            return None
        out.append(r)
    d = math.hypot(out[0][0] - out[1][0], out[0][1] - out[1][1])
    return {"dy": 0.5 * (out[0][0] + out[1][0]) / R, "dx": 0.5 * (out[0][1] + out[1][1]) / R,
            "agree": d / R, "agree_px": d, "peak": min(out[0][2], out[1][2]), "R": R}


# ------------------------------------------------------------- stage 2

def preview_stage(rec, wd, cache=None, progress=None, stop=None):
    """Adds the rows read from the previews to a stage-1 record."""
    cache = cache if cache is not None else _load_cache(wd)
    frames = rec.get("frames") or []
    items = rec["items"]

    def add(level, title, text, act=None):
        items.append({"level": level, "title": title, "text": text, "act": act})
    if len(frames) < 3:
        return rec
    ext = {os.path.splitext(f["name"])[1].lower() for f in frames}
    from .raw import RAW_EXTS
    if not (ext & RAW_EXTS):
        rec["no_previews"] = True
        add("info", "Preview checks", "need camera raw files with their built-in previews; "
            "framing and clipping are measured in the run instead")
        return rec
    total, dirty = len(frames), 0
    # a small Moon (short focal length) is measured on a larger preview:
    # decided once per set from its middle frames
    small = SMALL
    try:
        probe = [measure_frame(f["path"], SMALL) for f in frames[len(frames) // 2:len(frames) // 2 + 3]]
        rs = [q["R"] for q in probe if not q.get("err")]
        if rs and float(np.median(rs)) < SMALL_MIN_R:
            small = 2 * SMALL
    except Exception:
        pass
    rec["preview_px"] = small
    for i, f in enumerate(frames):
        if stop is not None and stop():
            return rec
        try:
            sig = _sig(f["path"])
        except OSError:
            continue
        c = cache["prev"].get(f["name"])
        if not c or c.get("sig") != sig or c.get("s") != small:
            try:
                c = dict(measure_frame(f["path"], small), sig=sig, s=small)
            except Exception as e:
                c = {"err": f"{type(e).__name__}: {e}", "sig": sig, "s": small}
            cache["prev"][f["name"]] = c
            dirty += 1
            if dirty % 12 == 0:
                _save_cache(wd, cache)
        f["pv"] = c
        if progress is not None:
            progress(i + 1, total)
    if dirty:
        _save_cache(wd, cache)
    good = [f for f in frames if f.get("pv") and not f["pv"].get("err")]
    bad = [f for f in frames if not f.get("pv") or f["pv"].get("err")]
    if len(good) < max(3, 0.5 * len(frames)):
        why = (bad[0]["pv"].get("err") if bad and bad[0].get("pv") else "no preview")
        rec["no_previews"] = True
        add("info", "Preview checks", f"the Moon was found in only {len(good)} of {len(frames)} previews "
            f"({why}); framing and clipping are measured in the run instead")
        return rec
    Rm = float(np.median([f["pv"]["R"] for f in good]))
    timed = rec.get("timed")
    W, H = good[0]["pv"]["w"], good[0]["pv"]["h"]
    rw = _raw_width(good[0]["path"])
    px = (rw / W) if rw else None                    # full-size px per preview px
    rec["geometry"] = {"R_frac_h": round(Rm / H, 4), "R_px": round(Rm * px, 1) if px else None}

    def as_px(v_R):
        return f" ({abs(v_R) * Rm * px:.0f} px)" if px else ""
    # 1. framing: the Moon's size and the nearest edge
    cyc = float(np.median([f["pv"]["cy"] for f in good]))
    cxc = float(np.median([f["pv"]["cx"] for f in good]))
    edges = {"top": cyc / Rm, "bottom": (H - cyc) / Rm, "left": cxc / Rm, "right": (W - cxc) / Rm}
    near = min(edges, key=edges.get)
    rec["geometry"]["edges_R"] = {k: round(v, 2) for k, v in edges.items()}
    txt = (f"the Moon's radius is {100 * Rm / H:.1f} % of the frame height"
           + (f" ({Rm * px:.0f} px, about {975.0 / (Rm * px):.2f} arcsec per px)" if px else "")
           + f"; the frame reaches {min(edges.values()):.1f} R ({near}) to {max(edges.values()):.1f} R")
    if edges[near] < TIGHT_R:
        rec["tight"] = {"R": round(edges[near], 1), "side": near}
        add("info", "Framing", txt + f". Tight: the corona runs off the frame beyond {edges[near]:.1f} R")
    else:
        add("ok", "Framing", txt)

    # 2. frames with a bright spot at the limb, or brighter than their group
    # frames that recorded the same light: same exposure time, ISO and aperture.
    # The key sorts by the light gathered (exposure x ISO / f-number squared).
    iso0 = float(np.median([f["iso"] for f in frames if f.get("iso")] or [100.0]))
    fn0 = float(np.median([f["fnum"] for f in frames if f.get("fnum")] or [8.0]))

    def gkey(f):
        eff = f["sec"] * (float(f["iso"]) / iso0 if f.get("iso") else 1.0) \
            * ((fn0 / float(f["fnum"])) ** 2 if f.get("fnum") else 1.0)
        return (float(f"{eff:.4g}"), _key(f["sec"]), int(f["iso"] or 0), round(float(f["fnum"] or 0.0), 1))
    mixed = len({(k[2], k[3]) for k in map(gkey, frames)}) > 1

    def gname(k):
        return _ex_txt(k[1]) + (f" (ISO {k[2]}" + (f", f/{k[3]:g}" if k[3] else "") + ")" if mixed and k[2] else "")
    by_exp = {}
    for f in good:
        by_exp.setdefault(gkey(f), []).append(f)
    ladder = sorted(by_exp)
    med_area = {k: float(np.median([f["pv"]["clip_area"] for f in by_exp[k]])) for k in ladder}
    med_far = {k: float(np.median([f["pv"]["far"] for f in by_exp[k]])) for k in ladder}
    mad_area = {k: 1.4826 * float(np.median(np.abs(np.array([f["pv"]["clip_area"] for f in by_exp[k]]) - med_area[k])))
                for k in ladder}
    med_R = {k: float(np.median([f["pv"]["R"] for f in by_exp[k]])) for k in ladder}
    flagged = []
    # THE DARK DISC ITSELF. In totality the limb is one full circle of one
    # size; a crescent, a diamond ring or a string of beads breaks it. On a
    # test set's contact frames the disc was not found at all in two of
    # three and came out 10 % too large in the third, while white clipping
    # outside the limb was no larger than a prominence gives on other sets.
    for f in bad:
        if not f.get("pv") or f["pv"].get("err") != "no lunar disc found":
            continue
        k = gkey(f)
        if k in by_exp or any(kk > k for kk in ladder):
            f["flag"] = "no complete dark disc in its preview"
            flagged.append(f)
    # A BEAD IS THERE IN THE SHORT EXPOSURES TOO. A prominence that the Moon
    # is uncovering or covering clips white in the longer exposures only, and
    # more at one end of totality than at the other: on one test set that
    # flagged 18 good frames. The photosphere is brighter by orders of
    # magnitude: a spot counts when the frame is one of the two shortest
    # exposures, or a shorter exposure shot within SPOT_PEER_S shows it too.
    def in_short_too(f, k):
        if ladder.index(k) <= 1:
            return True
        if f.get("t") is None:
            return False
        for g in good:
            if (g is not f and gkey(g) < k and g.get("t") is not None and abs(g["t"] - f["t"]) <= SPOT_PEER_S
                    and g["pv"]["clip_area"] > 0.0003 and g["pv"]["clip_off"] > SPOT_OFF):
                return True
        return False
    for f in good:
        k = gkey(f)
        pv = f["pv"]
        why = None
        one_sided = pv["clip_off"] > SPOT_OFF
        if ((len(by_exp[k]) >= 3 and abs(pv["R"] - med_R[k]) > DISC_DEV_GROUP * med_R[k])
                or abs(pv["R"] - Rm) > DISC_DEV_SET * Rm):
            # near an end of the set: a contact; in mid-totality: a smeared
            # limb (two 2 s frames of one test set, shaken: the disc 3 %
            # smaller, the limb fit three times as rough)
            t_ = f.get("t")
            ends = t_ is None or not timed or min(t_ - frames[0]["t"], frames[-1]["t"] - t_) <= SPOT_END_S
            if ends:
                f["flag"] = "the dark disc is not the size it has in the other frames"
                f["spot"] = True
                flagged.append(f)
            else:
                f["smear"] = True
            continue
        if pv.get("clip_in", 0.0) > SPOT_IN:
            why = "a bright spot at the limb that spills onto the Moon's disc"
        elif (one_sided and len(by_exp[k]) >= 4
              and pv["clip_area"] - med_area[k] > SPOT_EXCESS + 4.0 * mad_area[k]
              and in_short_too(f, k)):
            why = "a bright spot on one side of the limb that the other frames at this exposure do not show"
        elif len(by_exp[k]) >= 3 and pv["clip_area"] > 3.0 * med_area[k] + 0.03:
            why = "far more of it is clipped than in the other frames at this exposure"
        elif len(by_exp[k]) >= 3 and med_far[k] > 0 and pv["far"] > 1.6 * med_far[k] + 0.002:
            why = "its outer field is much brighter than in the other frames at this exposure"
        else:
            j = ladder.index(k)
            if j + 1 < len(ladder) and pv["clip_area"] > 1.5 * med_area[ladder[j + 1]] + 0.03:
                why = "more of it is clipped than of the next longer exposure"
        if why:
            f["flag"] = why
            f["spot"] = why.startswith("a bright spot")
            flagged.append(f)
    # A bead seen in one frame is there in its neighbours too: frames shot
    # within SPOT_NEAR_S of a flagged one go with it, and so does everything
    # between it and the end of the set it is close to (beads fade after second
    # contact and grow toward third). Only from frames flagged for what they
    # show, and only near the ends: a frame flagged in mid-totality is one
    # frame's problem.
    if flagged and timed_all(good):
        t_lo, t_hi = good[0]["t"], good[-1]["t"]
        for f in [x for x in flagged if x.get("spot")]:
            to_end = (t_hi - f["t"]) <= SPOT_END_S
            to_start = (f["t"] - t_lo) <= SPOT_END_S
            for g in good:
                if g.get("flag"):
                    continue
                d = g["t"] - f["t"]
                if abs(d) <= SPOT_NEAR_S or (to_end and d > 0) or (to_start and d < 0):
                    g["flag"] = ("shot within %.1f s of a frame with a bright spot at the limb" % SPOT_NEAR_S
                                 if abs(d) <= SPOT_NEAR_S else "shot between such a frame and the %s of the set"
                                 % ("end" if d > 0 else "start"))
                    flagged.append(g)
    flagged.sort(key=lambda f: f["order"])
    for f in flagged:
        rec["roles"][f["name"]] = {"role": "contact", "why": f["flag"]}
    if flagged:
        names = ", ".join(f["name"] for f in flagged[:5]) + (f" and {len(flagged) - 5} more" if len(flagged) > 5 else "")
        first = next((f for f in flagged if not f["flag"].startswith("shot ")), flagged[0])
        n_a = sum(1 for f in flagged if f["order"] < 0.5 * len(frames))
        n_b = len(flagged) - n_a
        where = (f"{n_a} at the start and {n_b} at the end of the set" if n_a and n_b else
                 "at the start of the set" if n_a else "at the end of the set")
        add("warn", "Not totality?", f"{len(flagged)} frame{'s' if len(flagged) != 1 else ''} ({where}) "
            f"look{'s' if len(flagged) == 1 else ''} like diamond ring or beads. {first['name']}: {first['flag']}. "
            f"Frames: {names}. "
            "The run tests each frame again; switch them off in Frames to be sure", {"field": "fpanel"})
    else:
        add("ok", "Totality", "no frame shows a diamond ring or beads in its preview")
    smeared = [f for f in good if f.get("smear")]
    if smeared:
        add("warn", "Smeared frames?", f"{len(smeared)} frame{'s' if len(smeared) != 1 else ''} in mid-totality "
            f"show{'s' if len(smeared) == 1 else ''} the Moon's edge blurred against the other frames at the "
            "exposure (camera shake, wind?): " + ", ".join(f["name"] for f in smeared[:6])
            + (f" and {len(smeared) - 6} more" if len(smeared) > 6 else "")
            + ". Look at them in Frames and switch them off if they are", {"field": "fpanel"})
        for f in smeared:
            rec["roles"][f["name"]] = {"role": "off", "why": "the Moon's edge is blurred in its preview"}

    # 3. clipping per exposure (previews: lower limits)
    clean = {k: [f for f in by_exp[k] if not f.get("flag")] or by_exp[k] for k in ladder}
    limb = {k: float(np.median([f["pv"]["clip_limb"] for f in clean[k]])) for k in ladder}
    reach = {k: float(np.median([f["pv"]["clip_r"] for f in clean[k]])) for k in ladder}
    rec["clipping"] = [{"sec": k[1], "iso": k[2], "fnum": k[3], "limb": round(limb[k], 2),
                        "to_R": round(reach[k], 2)} for k in ladder]
    free = [k for k in ladder if limb[k] < 0.05]
    clp = [k for k in ladder if limb[k] >= 0.5]
    edge_R = max(edges.values())
    if not free:
        rec["clipped_all"] = gname(ladder[0])
        add("warn", "Clipping", f"even the shortest exposure ({gname(ladder[0])}) is clipped at the limb "
            "in its preview (the raw data hold more): the innermost corona and the prominences may be burnt out")
    else:
        txt = (f"unclipped at the limb up to {gname(free[-1])}" if len(free) < len(ladder) else
               "no exposure is clipped at the limb")
        if clp:
            txt += "; clipped out to " + ", ".join(
                f"{min(reach[k], edge_R):.1f} R at {gname(k)}" for k in clp[-3:])
        add("ok", "Clipping", txt + " (camera previews; the raw data hold more)")

    # 4. did the framing move? corona links between the groups in time
    if rec.get("style") in ("repeated", "complementary"):
        groups = [[f for f in good if f["bracket"] == b] for b in range(rec["n_brackets"])]
    else:
        groups, cur = [], []
        for f in good:
            if cur and timed and f["t"] - cur[-1]["t"] > GAP_MIN_S:
                groups.append(cur)
                cur = []
            cur.append(f)
        if cur:
            groups.append(cur)
        if len(groups) == 1 and len(good) >= 8:
            q = len(good) // 4
            groups = [good[:q], good[q:2 * q], good[2 * q:3 * q], good[3 * q:]]
    groups = [g for g in groups if g]
    links = []
    if len(groups) >= 2:
        gi = list(range(len(groups)))
        if len(gi) > MAX_PAIRS + 1:
            gi = sorted({int(round(v)) for v in np.linspace(0, len(groups) - 1, MAX_PAIRS + 1)})
        pairs = list(zip(gi, gi[1:]))          # a chain, so the links add up over the set

        def pick(ga, gb):
            c = []
            for a in ga:
                for b in gb:
                    if a.get("flag") or b.get("flag") or a.get("smear") or b.get("smear"):
                        continue
                    de = abs(math.log2(gkey(a)[0] / gkey(b)[0]))
                    if de > 2.6:
                        continue
                    c.append(((round(de, 1), -min(a["pv"]["usable"], b["pv"]["usable"])), a, b))
            c.sort(key=lambda t: t[0])
            return c[:3]
        for (i, j) in pairs:
            if stop is not None and stop():
                return rec
            cs = None
            for _, a, b in pick(groups[i], groups[j]):
                try:
                    cs = corona_shift(b["path"], b["pv"], a["path"], a["pv"], small)     # later against earlier
                except Exception:
                    cs = None
                if cs is not None:
                    break
            if cs is None:
                continue
            dt = (b["t"] - a["t"]) if timed else None
            moon = ((b["pv"]["cy"] - a["pv"]["cy"]) / Rm, (b["pv"]["cx"] - a["pv"]["cx"]) / Rm)
            rel = math.hypot(moon[0] - cs["dy"], moon[1] - cs["dx"])     # Moon against the corona
            # the check on the measurement itself: the Moon cannot have moved
            # against the corona by more than its own motion allows
            lim = (LUNAR_RATE_R * abs(dt) if dt is not None else 0.05) * 2.0 + 0.02
            sure = cs["agree_px"] <= LINK_AGREE_PX and cs["peak"] >= LINK_PEAK and rel <= lim
            links.append({"a": i, "b": j, "frame_a": a["name"], "frame_b": b["name"], "dt": dt,
                          "dy": round(cs["dy"], 4), "dx": round(cs["dx"], 4),
                          "shift": round(math.hypot(cs["dy"], cs["dx"]), 4),
                          "moon_rel": round(rel, 4), "peak": round(cs["peak"], 3),
                          "agree_px": round(cs["agree_px"], 2), "sure": bool(sure)})
    rec["links"] = links
    cons = [l for l in links if l["sure"]]
    unsure = [l for l in links if not l["sure"]]
    if cons:
        worst = max(cons, key=lambda l: l["shift"])
        # over the whole set: the links added up along the chain
        total = math.hypot(sum(l["dy"] for l in cons), sum(l["dx"] for l in cons))
        rec["drift_R"] = round(total, 4)
        what = "brackets" if rec.get("style") in ("repeated", "complementary") else "parts of the set"
        if worst["shift"] > MOVE_WARN_R and (len(cons) == 1 or worst["shift"] > 3.0 * float(
                np.median([l["shift"] for l in cons])) + 0.004):
            add("warn", "Framing moved", f"between {what} {worst['a'] + 1} and {worst['b'] + 1} the corona sits "
                f"{100 * worst['shift']:.1f} % of the Moon's radius{as_px(worst['shift'])} apart"
                + (f", {abs(worst['dt']):.0f} s between the two frames compared" if worst["dt"] is not None else "")
                + ". The run has to link those brackets on the corona itself")
            rec["suggest"]["cross_bracket"] = {"value": "corona", "why": f"framing moved by {100 * worst['shift']:.1f} % of R"}
            rec["moved"] = {"R": round(worst["shift"], 4), "px": (round(worst["shift"] * Rm * px) if px else None),
                            "a": worst["a"] + 1, "b": worst["b"] + 1}
        elif total > DRIFT_INFO_R:
            span = rec.get("span_s") or 0
            add("info", "Framing drifts", f"the corona drifts by {100 * total:.1f} % of the Moon's radius"
                f"{as_px(total)} from the first to the last part of the set"
                + (f" ({60 * 100 * total / span:.1f} % per minute)" if span > 20 else "")
                + ": the frames are aligned one by one, and the frame edge cuts a changing part of the corona")
        else:
            add("ok", "Framing steady", f"the corona stays within {100 * worst['shift']:.1f} % of the Moon's radius"
                f"{as_px(worst['shift'])} between {what}"
                + (f", {100 * total:.1f} % over the whole set" if len(cons) > 1 else ""))
    elif rec.get("pauses") or rec.get("style") in ("repeated", "complementary"):
        add("info", "Framing across the set", ("not measured reliably from the previews "
            f"({len(unsure)} comparison{'s' if len(unsure) != 1 else ''} failed its own checks)" if unsure else
            "no two frames of similar exposure to compare from the previews") + "; the run measures it")
    return rec


# ---- SHARPNESS FROM THE RAW FILES (0.32.0, 2026-11-15) -----------------
# The previews are too small to tell a sharp limb from a shaken one at long
# focal lengths: on a test set (Moon radius 1270 px, low Sun, wind) the edge
# of the sharp frames was 6-7 px wide and of the shaken ones 10-26 px, i.e.
# 0.5-2 preview px. So each frame's raw is read once (a few seconds each), the
# limb fitted, and the width of the Moon's edge measured along 360 rays: the
# 20-80 % rise from the dark disc to the corona just outside it. Shake widens
# it in one direction only (two opposite sectors wide, the ones at 90 degrees
# sharp); seeing widens it all round. A frame whose edge is more than
# SHARP_FACTOR times as wide as the set's sharp frames (25th percentile) is
# proposed to be switched off -- never all frames of one exposure.
SHARP_BUILD = 1
SHARP_FACTOR = 1.4
# LARGE DISC ONLY (0.31.2). An exposure whose frames are all soft keeps its
# sharpest one. On a test set with the Moon filling the frame two exposures
# consisted of ONE frame each, both blurred (edge 13 and 15 px against 6-8 px
# in every other frame), and sat between sharp exposures 2.5x apart. Such a
# frame is proposed off as well, leaving the exposure out, when: the disc is
# larger than LONE_DISC_FRAC of the short side; a sharp exposure exists on
# either side and the two are no more than LONE_MAX_GAP apart; and the frame's
# edge is more than SHARP_FACTOR times as wide as theirs. Smaller discs keep
# the old rule unchanged.
LONE_DISC_FRAC = 0.5
LONE_MAX_GAP = 4.0


def _edge_width(L, cy, cx, R):
    """Median 20-80 % width of the limb (px of L) over 360 rays, and the
    median in 12 sectors of 30 degrees (None where no ray was usable)."""
    from scipy import ndimage
    h, w = L.shape
    pad = int(0.08 * R) + 4
    y0, y1 = max(int(cy - R - pad), 0), min(int(cy + R + pad) + 1, h)
    x0, x1 = max(int(cx - R - pad), 0), min(int(cx + R + pad) + 1, w)
    sm = ndimage.gaussian_filter(np.asarray(L[y0:y1, x0:x1], np.float32), 1.0)
    ang = np.linspace(0, 2 * np.pi, 360, endpoint=False)
    rr = np.linspace(0.93 * R, 1.07 * R, 281)
    ys = (cy - y0) + rr[:, None] * np.sin(ang)[None]
    xs = (cx - x0) + rr[:, None] * np.cos(ang)[None]
    ok = (ys.min(0) > 0) & (ys.max(0) < sm.shape[0] - 1) & (xs.min(0) > 0) & (xs.max(0) < sm.shape[1] - 1)
    P = ndimage.map_coordinates(sm, [ys.ravel(), xs.ravel()], order=1).reshape(ys.shape)
    dr = float(rr[1] - rr[0])
    wid = np.full(360, np.nan)
    for k in np.flatnonzero(ok):
        p = P[:, k]
        lo = float(np.median(p[:60]))
        hi = float(np.max(p[140:200]))
        if hi - lo < 10.0 * max(float(np.std(p[:60])), 1.0):
            continue
        q = (p[100:] - lo) / (hi - lo)
        i20, i80 = int(np.argmax(q > 0.2)), int(np.argmax(q > 0.8))
        if i80 > i20:
            wid[k] = (i80 - i20) * dr
    good = wid[np.isfinite(wid)]
    if good.size < 90:
        return None, None
    sec = []
    for i in range(12):
        v = wid[i * 30:(i + 1) * 30]
        v = v[np.isfinite(v)]
        sec.append(round(float(np.median(v)), 2) if v.size >= 5 else None)
    return float(np.median(good)), sec


def measure_raw_edge(path, cy, cx, R):
    """The limb of one raw frame at half size, seeded from its preview:
    {"w": edge width, "sec": per-sector widths, "R": radius} in half-size px."""
    import rawpy
    from .pipeline import fit_limb_rays
    with rawpy.imread(path) as r:
        b = r.raw_image_visible.astype(np.float32)
        blk = float(np.mean(r.black_level_per_channel))
    b -= blk
    h, w = b.shape
    L = b[:h // 2 * 2, :w // 2 * 2].reshape(h // 2, 2, w // 2, 2).mean(axis=(1, 3))
    del b
    f = fit_limb_rays(np.clip(L, 0, None), cy, cx, R)
    if f is None or not np.isfinite(f[2]) or f[3] > 0.05 * f[2] or abs(f[2] - R) > 0.15 * R:
        return {"err": "limb not fitted"}
    wd_, sec = _edge_width(L, float(f[0]), float(f[1]), float(f[2]))
    if wd_ is None:
        return {"err": "limb too faint"}
    return {"w": round(wd_, 3), "sec": sec, "R": round(float(f[2]), 2)}


def sharp_stage(rec, wd, cache=None, progress=None, stop=None):
    """Adds the sharpness rows (raw files) to a record that has been through
    preview_stage. Proposes shaken or soft frames to be switched off."""
    cache = cache if cache is not None else _load_cache(wd)
    frames = rec.get("frames") or []
    if rec.get("no_previews") or len(frames) < 3:
        return rec
    from .raw import RAW_EXTS
    good = [f for f in frames if f.get("pv") and not f["pv"].get("err")
            and os.path.splitext(f["name"])[1].lower() in RAW_EXTS]
    if len(good) < 3:
        return rec
    rw = _raw_width(good[0]["path"])
    if not rw:
        return rec
    k_ = rw / float(good[0]["pv"]["w"]) / 2.0          # half-size px per preview px
    cr = cache.setdefault("raw", {})
    dirty = 0
    for i, f in enumerate(good):
        if stop is not None and stop():
            return rec
        try:
            sig = _sig(f["path"])
        except OSError:
            continue
        c = cr.get(f["name"])
        if not c or c.get("sig") != sig or c.get("b") != SHARP_BUILD:
            pv = f["pv"]
            try:
                c = measure_raw_edge(f["path"], pv["cy"] * k_, pv["cx"] * k_, pv["R"] * k_)
            except Exception as e:
                c = {"err": f"{type(e).__name__}: {e}"}
            c.update(sig=sig, b=SHARP_BUILD)
            cr[f["name"]] = c
            dirty += 1
            if dirty % 8 == 0:
                _save_cache(wd, cache)
        f["sh"] = c
        if progress is not None:
            progress(i + 1, len(good))
    if dirty:
        _save_cache(wd, cache)
    meas = [f for f in good if f.get("sh") and f["sh"].get("w") and not f.get("flag")]
    if len(meas) < max(3, 0.5 * len(good)):
        return rec
    items = rec["items"]
    floor = float(np.percentile([f["sh"]["w"] for f in meas], 25))
    lim = SHARP_FACTOR * floor
    by_exp = {}
    for f in meas:
        by_exp.setdefault((_key(f["sec"]), f.get("iso"), f.get("fnum")), []).append(f)
    # EACH FRAME AGAINST ITS OWN EXPOSURE (2026-11-16). The Moon's edge widens
    # with the exposure: on a seven-bracket test set it was 4.1-4.8 px in every
    # frame up to 1/2 s, 5.0 px in all six 1 s frames and 6.3-6.5 px in all six
    # 2 s frames. Measured against the sharp quarter of the whole set (4.3 px,
    # limit 6.0) five of the six 2 s frames were proposed off as shaken -- 10 of
    # 25 s of exposure, all of it in the outer corona. So the reference for a
    # frame is the sharpest frame of its own exposure or of any longer one
    # (never sharper than the set's floor): only a frame that stands out
    # among frames shot like it is proposed off, and an exposure that is soft
    # all through while the longer ones are sharp is still found.
    _sharpest = {k: min(f["sh"]["w"] for f in fs) for k, fs in by_exp.items()}
    lim_of = {}
    for k in by_exp:
        _base = min(v for k2, v in _sharpest.items() if k2[0] >= k[0])
        lim_of[k] = SHARP_FACTOR * max(floor, _base)
    soft, all_soft = [], []
    for k, fs in by_exp.items():
        bad = [f for f in fs if f["sh"]["w"] > lim_of[k]]
        if bad and len(bad) == len(fs):
            keep = min(fs, key=lambda f: f["sh"]["w"])
            bad = [f for f in bad if f is not keep]
            all_soft.append((k, keep))
        soft += bad
    lone = []
    try:
        _fr = [2.0 * f["pv"]["R"] / max(1.0, min(f["pv"]["w"], f["pv"]["h"])) for f in good
               if f["pv"].get("R") and f["pv"].get("h")]
        if all_soft and _fr and float(np.median(_fr)) > LONE_DISC_FRAC:
            _softk = {k for k, _ in all_soft}
            _softn = {f["name"] for f in soft}
            for k, keep in all_soft:
                if len(by_exp[k]) != 1:
                    continue
                nb = {}
                for k2, fs2 in by_exp.items():
                    if k2 in _softk or k2[1:] != k[1:]:
                        continue
                    ok2 = [f["sh"]["w"] for f in fs2 if f["name"] not in _softn]
                    if not ok2:
                        continue
                    side = "lo" if k2[0] < k[0] else "hi"
                    if side not in nb or abs(np.log(k2[0] / k[0])) < abs(np.log(nb[side][0] / k[0])):
                        nb[side] = (k2[0], ok2)
                if len(nb) == 2 and nb["hi"][0] / nb["lo"][0] <= LONE_MAX_GAP \
                        and keep["sh"]["w"] > SHARP_FACTOR * float(np.median(nb["lo"][1] + nb["hi"][1])):
                    lone.append((k, keep))
    except Exception:
        lone = []
    if lone:
        _lk = {k for k, _ in lone}
        all_soft = [(k, f) for k, f in all_soft if k not in _lk]
        soft += [f for _, f in lone]
    _lone_names = {f["name"] for _, f in lone}
    soft.sort(key=lambda f: f.get("ep") or 0)
    # FRAMES AN EARLIER PRE-CHECK SWITCHED OFF (2026-11-16). A frame switched
    # off on this check's proposal is no longer listed, so it is never looked at
    # again -- and when the rule changes (it did: 2 s frames of a test set were
    # taken for shaken), it stays off for good. Those frames, and only those
    # (their recorded reason is the check's own), are measured here against
    # the current rule; the ones that pass are offered back (rec["sharp_on"]).
    rec["sharp_on"] = {}
    try:
        from .frameselect import load as _fs_load
        _folder = os.path.dirname(good[0]["path"])
        _sel = _fs_load(_folder)
        _why = _sel.get("why") or {}
        _cand = [k for k in sorted(_sel["off"]) if "/" not in k
                 and str(_why.get(k, "")).startswith("survey: ")
                 and any(t in str(_why.get(k, "")) for t in ("shaken", "soft all round", "blurred"))]
        _pv = np.median([[f["pv"]["cy"], f["pv"]["cx"], f["pv"]["R"]] for f in good], axis=0)
        _new = 0
        for k in _cand:
            if stop is not None and stop():
                break
            _p = os.path.join(_folder, k)
            if not os.path.exists(_p):
                continue
            _sg = _sig(_p)
            _m = cache.setdefault("meta", {}).get(k)
            if not _m or _m.get("sig") != _sg:
                _m = dict(frame_meta(_p), sig=_sg)
                cache["meta"][k] = _m
            c = cr.get(k)
            if not c or c.get("sig") != _sg or c.get("b") != SHARP_BUILD:
                try:
                    c = measure_raw_edge(_p, _pv[0] * k_, _pv[1] * k_, _pv[2] * k_)
                except Exception as e:
                    c = {"err": f"{type(e).__name__}: {e}"}
                c.update(sig=_sg, b=SHARP_BUILD)
                cr[k] = c
                _new += 1
            if not c.get("w") or not _m.get("sec"):
                continue
            _kk = (_key(_m["sec"]), _m.get("iso"), _m.get("fnum"))
            if _kk in lim_of and c["w"] <= lim_of[_kk]:
                rec["sharp_on"][k] = round(2 * c["w"], 1)
        if _new:
            _save_cache(wd, cache)
    except Exception:
        rec["sharp_on"] = {}

    def kind(f):
        sec = [v for v in f["sh"]["sec"] if v]
        return "shaken" if sec and max(sec) >= 1.8 * min(sec) else "soft all round"
    # the raw measurement replaces the preview's smear guess
    rec["items"] = items = [it for it in items if it.get("title") != "Smeared frames?"]
    for n, r in list(rec["roles"].items()):
        if r.get("role") == "off" and "preview" in str(r.get("why", "")):
            del rec["roles"][n]
    for f in soft:
        if rec["roles"].get(f["name"], {}).get("role") == "contact":
            continue
        _k = (_key(f["sec"]), f.get("iso"), f.get("fnum"))
        rec["roles"][f["name"]] = {"role": "off", "why": f"{kind(f)}: the Moon's edge is "
                                   f"{2 * f['sh']['w']:.0f} px wide, "
                                   f"{2 * lim_of[_k] / SHARP_FACTOR:.0f} px in the sharp frames shot like it"
                                   + ("; the only frame at this exposure, the exposures next to it cover it"
                                      if f["name"] in _lone_names else "")}
    n_sh = sum(1 for f in soft if kind(f) == "shaken")
    rec["sharpness"] = {"floor_px": round(2 * floor, 2), "limit_px": round(2 * lim, 2),
                        "limit_px_by_exposure": {_ex_txt(k[0]): round(2 * v, 2) for k, v in sorted(lim_of.items())},
                        "frames": {f["name"]: round(2 * f["sh"]["w"], 2) for f in meas},
                        "soft": [f["name"] for f in soft],
                        "all_soft": [_ex_txt(k[0]) for k, _ in all_soft],
                        "lone_out": [_ex_txt(k[0]) for k, _ in lone]}
    base = (f"measured on the raw files: the Moon's edge is {2 * floor:.0f} px wide in the sharp frames "
            f"(the limit set by optics, focus and the air)")
    if soft:
        names = ", ".join(f["name"] for f in soft[:8]) + (f" and {len(soft) - 8} more" if len(soft) > 8 else "")
        txt = (base[0].upper() + base[1:] + f". In {len(soft)} of {len(meas)} frames it is more than "
               f"{SHARP_FACTOR:g} times as wide as in the sharpest frames of the same exposure"
               + (f", {n_sh} of them blurred in one direction (shake)" if n_sh else "")
               + f". Switch them off in Frames: {names}")
        if all_soft:
            txt += (". Every frame is soft at " + ", ".join(_ex_txt(k[0]) for k, _ in all_soft)
                    + ": the sharpest of each stays in")
        if lone:
            txt += (". " + ", ".join(_ex_txt(k[0]) for k, _ in lone)
                    + (" has one frame only and it is blurred; the exposures next to it cover it, so it is "
                       "proposed off too" if len(lone) == 1 else
                       " have one frame each and it is blurred; the exposures next to them cover them, so "
                       "they are proposed off too"))
        items.append({"level": "warn", "title": "Blurred frames", "text": txt, "act": {"field": "fpanel"}})
    else:
        items.append({"level": "ok", "title": "Sharpness", "text": "all frames within "
                      f"{SHARP_FACTOR:g} times the sharpest of their exposure; " + base, "act": None})
    return rec


def timed_all(frames):
    return bool(frames) and all(f.get("t") is not None for f in frames)


def run(folder, wd, progress=None, stop=None, previews=True):
    """Both stages. Returns the record; record["done"] says whether the
    preview stage finished."""
    cache = _load_cache(wd)
    rec = exif_stage(folder, wd, cache)
    _save_cache(wd, cache)
    rec["done"] = False
    if previews:
        n1 = len(rec["items"])
        try:
            preview_stage(rec, wd, cache, progress, stop)
            try:
                sharp_stage(rec, wd, cache, progress, stop)
            except Exception as e:
                rec["items"].append({"level": "info", "title": "Sharpness",
                                     "text": f"not measured ({type(e).__name__}: {e})", "act": None})
            rec["done"] = not (stop is not None and stop())
        except Exception as e:
            del rec["items"][n1:]
            rec["items"].append({"level": "info", "title": "Preview checks",
                                 "text": f"not available ({type(e).__name__}: {e})", "act": None})
            rec["done"] = True
    return rec


def guidance(rec):
    """What the survey means for the user, from the record: the settings to
    change (tab 2 ids and values, with the reason), what to do with single
    frames (counts; the frames themselves are in rec["roles"]), what to expect
    from the result, the limits nothing can change, and "tries": settings to
    try only if the result shows a given problem, each with the way back. The page compares the
    settings with what is set and builds its verdict from that."""
    settings, expect, limits, tries = [], [], [], []
    sug = rec.get("suggest") or {}
    fpe = rec.get("frames_per_exposure")
    if "stack_combine" in sug:
        settings.append({"id": "stackCombine", "value": "mean", "name": "Tier combine",
                         "label": "Mean (no rejection)",
                         "why": (f"{fpe} frame per exposure" if fpe == 1 else f"only {fpe} frames per exposure")
                         + ": a per-pixel rejection has nothing to compare them with and would throw away good pixels"})
    elif fpe and fpe >= 5:
        settings.append({"id": "stackCombine", "value": "clip", "name": "Tier combine",
                         "label": "Clipped mean (kappa-sigma)",
                         "why": f"{fpe} frames per exposure: the clipped mean takes out satellites, planes "
                                "and single hot pixels"})
    # MERGE WEIGHT: NO TIP BEFORE THE RUN (2026-11-15). The conditional
    # "try Blended edge if rings show" with its button read as advice to switch
    # now (a tester: a user would fall for it, get the known Blended-edge
    # failure and blame the suggestion). Whether rings show can only be seen in
    # the result, so the tip moved to 3 Check your result, where it appears
    # only when the run measured the exposures disagreeing at the limb
    # (runnotes.py). Exact edge stays the default.
    if "tier_mode" in sug:
        mx = rec.get("mixed") or {}
        what = " and ".join(k for k in ("ISO", "Aperture") if k in mx) or "ISO"
        settings.append({"id": "tierMode", "value": "measured", "name": "Tiers",
                         "label": "By measured brightness",
                         "why": f"{what} changes within the set: grouped by shutter speed, frames that recorded "
                                "different light would be averaged together",
                         # what to look at after the run, and the way back
                         # (grouping.by_brightness logs every split and merge)
                         "check": "After the run, the log line \"tiers by measured brightness\" and the lines "
                                  "under it show which shutter speeds were split and which frames were put "
                                  "together. Each tier should be one combination of shutter speed, ISO and "
                                  "aperture you shot. If not, set Tiers back to By shutter speed (new stack)"})
    elif rec.get("ladder") and not rec.get("fits"):
        # TIERS FROM THE CAMERA (2026-11-14): the guide also says when a
        # setting left over from another folder is wrong here. "By measured
        # brightness" stayed set after a mixed-ISO test, the next folder (two
        # brackets, honest camera times) was stacked with it, and frames of
        # the two brackets ended up in one tier -- a doubled Moon. With camera
        # files and one ISO and aperture, the shutter speed is the grouping.
        settings.append({"id": "tierMode", "value": "exposure", "name": "Tiers",
                         "label": "By shutter speed",
                         "why": "the camera recorded every frame's exposure, and ISO and aperture do not change: "
                                "frames with the same shutter speed recorded the same light. By measured "
                                "brightness is for mixed ISO or aperture and for missing or wrong exposure data"})
    roles = rec.get("roles") or {}
    n_contact = sum(1 for r in roles.values() if r.get("role") == "contact")
    n_off = sum(1 for r in roles.values() if r.get("role") == "off")
    # what to expect
    if rec.get("style") == "complementary":
        mv = rec.get("moved")
        expect.append({"level": "info", "text":
                       "Each exposure exists in one bracket only, so the brackets are tied together on the "
                       "corona alone"
                       + ((f". The framing moved by {mv['px']} px between them" if mv.get("px") else
                           f". The framing moved by {100 * mv['R']:.1f} % of the Moon's radius between them")
                          + "; the run measures that on the corona and corrects it" if mv else
                          ". The run measures their offset on the corona")
                       + ". After the run, the log line \"cross-bracket links\" shows how they were joined"})
    elif rec.get("moved"):
        mv = rec["moved"]
        expect.append({"level": "info", "text":
                       f"The framing moved between brackets {mv['a']} and {mv['b']}"
                       + (f" by {mv['px']} px" if mv.get("px") else f" by {100 * mv['R']:.1f} % of the Moon's radius")
                       + ": the run measures the offset on the corona and corrects it"})
    if rec.get("tight"):
        t = rec["tight"]
        expect.append({"level": "info", "text":
                       f"Tight framing: the corona runs off the frame beyond {t['R']:.1f} R at the {t['side']}, "
                       "so the outer streamers end at the frame edge there"})
    if fpe == 1:
        expect.append({"level": "info", "text":
                       "One frame per exposure: nothing is averaged, so the result is grainier than from "
                       "repeated brackets. Grain is best removed at the end, after export"})
    if rec.get("wide_step"):
        w = rec["wide_step"]
        expect.append({"level": "info", "text":
                       f"Wide step of {w['ev']:g} EV from {_ex_txt(w['short'])} to {_ex_txt(w['long'])}: "
                       f"where the {_ex_txt(w['long'])} exposure clips, a noisier band in the inner corona"})
    if rec.get("clipped_all"):
        limits.append({"level": "warn", "text":
                       f"Even the shortest exposure ({rec['clipped_all']}) is clipped at the limb in its preview: "
                       "the innermost corona and the prominences may come out burnt (the raw data often hold "
                       "a little more)"})
    if (rec.get("mixed") or {}).get("Focal length"):
        limits.append({"level": "warn", "text":
                       "Frames at another focal length cannot be stacked with the others: switch them off in "
                       "Frames"})
    if rec.get("timed") is False:
        limits.append({"level": "warn", "text":
                       "Some frames carry no time stamp: the brackets are read from the file order, and the "
                       "Moon's motion and the camera clock are not corrected"})
    if len(rec.get("ladder") or []) == 1:
        limits.append({"level": "warn", "text": "One exposure only: there is nothing to merge into an HDR"})
    if rec.get("no_previews"):
        expect.append({"level": "info", "text":
                       "Framing, clipping and contact frames could not be read from previews; the run checks "
                       "them itself"})
    return {"settings": settings, "frames": {"contact": n_contact, "off": n_off},
            "expect": expect, "limits": limits, "tries": tries}


def public(rec):
    """The record without paths and per-frame working data, for the GUI."""
    out = {k: v for k, v in rec.items() if k not in ("frames",)}
    try:
        out["guide"] = guidance(rec)
    except Exception:
        out["guide"] = None
    return out


def report_text(rec):
    """The rows as plain text."""
    mark = {"ok": "ok  ", "info": "note", "warn": "WARN", "missing": "MISS"}
    lines = [f"{rec.get('n', 0)} frames"]
    for it in rec.get("items", []):
        lines.append(f"[{mark.get(it['level'], it['level'])}] {it['title']}: {it['text']}")
    for l in rec.get("links") or []:
        lines.append(f"       link {l['a'] + 1}->{l['b'] + 1}: {l['frame_a']} / {l['frame_b']}  "
                     f"corona ({100 * l['dy']:+.2f}, {100 * l['dx']:+.2f}) % R, Moon against corona "
                     f"{100 * l['moon_rel']:.2f} % R, peak {l['peak']:.2f}, bands agree {l['agree_px']:.2f} px"
                     + ("" if l["sure"] else "  [not trusted]"))
    return "\n".join(lines)


# ------------------------------------------------- background job (server)

_JOB = {"key": None, "thread": None, "rec": None, "done": 0, "total": 0, "gen": 0}
_LOCK = threading.Lock()


def _folder_key(folder):
    from .raw import list_raws
    try:
        raws = list_raws(folder)
        return folder + "|" + str(len(raws)) + "|" + "|".join(
            f"{os.path.basename(p)}:{_sig(p)}" for p in (raws[:3] + raws[-3:])) + f"|{SURVEY_BUILD}"
    except Exception:
        return folder


def proposals(folder):
    """{frame name: {"role": "contact" | "short" | "off", "why": ...}} from
    the finished survey of this folder; {} while it runs or for another folder."""
    try:
        key = _folder_key(folder)
    except Exception:
        return {}
    with _LOCK:
        rec = _JOB["rec"]
        if _JOB["key"] == key and rec and rec.get("done"):
            return dict(rec.get("roles") or {})
    return {}


def reset(wd):
    """"Check again" on tab 1 (0.32.0): forget the survey of this folder,
    rows and cached preview readings alike, so the next status() reads every
    frame afresh. A job still running sees its generation change and stops."""
    with _LOCK:
        _JOB["gen"] += 1
        _JOB.update(key=None, rec=None, done=0, total=0, thread=None)
    _save_cache(wd, {"build": SURVEY_BUILD, "meta": {}, "prev": {}})


def status(folder, wd, start=True):
    """For /api/survey: the rows so far, and whether the preview stage is
    still running. Starts the job when the folder (or its frames) changed."""
    key = _folder_key(folder)
    with _LOCK:
        if _JOB["key"] != key and not start:
            return {"ok": True, "items": [], "pending": True, "waiting": True, "done": 0, "total": 0,
                    "suggest": {}, "roles": {}}
        if _JOB["key"] != key:
            _JOB["gen"] += 1
            gen = _JOB["gen"]
            _JOB.update(key=key, rec=None, done=0, total=0, thread=None)
            if True:
                def work():
                    try:
                        cache = _load_cache(wd)
                        rec = exif_stage(folder, wd, cache)
                        _save_cache(wd, cache)
                        rec["done"] = False
                        with _LOCK:
                            if _JOB["gen"] != gen:
                                return
                            _JOB["rec"] = dict(public(rec), items=list(rec["items"]), done=False)

                        def prog(i, n):
                            _JOB["done"], _JOB["total"] = i, n
                        n1 = len(rec["items"])
                        _JOB["phase"] = "previews"
                        try:
                            preview_stage(rec, wd, cache, prog, lambda: _JOB["gen"] != gen)
                        except Exception as e:
                            del rec["items"][n1:]
                            rec["items"].append({"level": "info", "title": "Preview checks",
                                                 "text": f"not available ({type(e).__name__}: {e})", "act": None})
                        # the preview rows are shown while the raw files are read
                        with _LOCK:
                            if _JOB["gen"] != gen:
                                return
                            _JOB["rec"] = dict(public(rec), items=list(rec["items"]), done=False)
                            _JOB["phase"] = "raw"
                            _JOB["done"], _JOB["total"] = 0, 0
                        try:
                            sharp_stage(rec, wd, cache, prog, lambda: _JOB["gen"] != gen)
                        except Exception as e:
                            rec["items"].append({"level": "info", "title": "Sharpness",
                                                 "text": f"not measured ({type(e).__name__}: {e})", "act": None})
                        rec["done"] = True
                        with _LOCK:
                            if _JOB["gen"] == gen:
                                _JOB["rec"] = public(rec)
                    except Exception as e:
                        with _LOCK:
                            if _JOB["gen"] == gen:
                                _JOB["rec"] = {"items": [{"level": "info", "title": "Survey",
                                                          "text": f"not available ({type(e).__name__}: {e})",
                                                          "act": None}], "done": True}
                t = threading.Thread(target=work, name="survey", daemon=True)
                _JOB["thread"] = t
                t.start()
        rec = _JOB["rec"]
        return {"ok": True, "items": (rec or {}).get("items", []),
                "pending": not (rec or {}).get("done", False),
                "done": _JOB["done"], "total": _JOB["total"], "phase": _JOB.get("phase"),
                "suggest": (rec or {}).get("suggest", {}), "roles": (rec or {}).get("roles", {}),
                "guide": (rec or {}).get("guide"),
                # frames an earlier pre-check switched off that now pass
                # (sharp_stage); the server offers them back as "Switch back on"
                "sharp_on": (rec or {}).get("sharp_on") or {}}


if __name__ == "__main__":
    import sys
    import time
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    if not args:
        print("usage: python -m eclipseforgehdr.survey FOLDER [WORKDIR] [--no-previews] [--json]")
        raise SystemExit(2)
    fol = os.path.abspath(os.path.expanduser(args[0]))
    if len(args) > 1:
        wdir = args[1]
    else:
        from .brand import WORKDIR_NAME
        wdir = os.path.join(fol, WORKDIR_NAME)
    t_ = time.time()
    r_ = run(fol, wdir, previews="--no-previews" not in sys.argv)
    if "--json" in sys.argv:
        print(json.dumps(public(r_), indent=1))
    else:
        print(report_text(r_))
        print(f"({time.time() - t_:.1f} s)")
