"""The shooting timeline of a folder, for the strip under the log (2026-10-15).

What EclipseForgeHDR made of a folder used to be a few lines in the log:
"repeated brackets: 7 brackets of up to 12 exposures over 1 min 26 s". A
picture says it at a glance -- one bar per frame along the time axis, its
height the exposure, so a repeated bracket shows as a saw and a set shot tier
by tier as a staircase -- and shows at once when something was read wrongly:
a frame in the wrong place, a bracket split in two, a gap.

Before a run this needs only the time stamps and exposure times (read once,
kept in the work folder). After a run the same strip also shows what happened
to each frame -- stacked, left out for brightness, left out by the sharpness cut, not
a totality frame -- and the instant whose Moon every tier holds.
"""
import os
import json

from .brackets import find_brackets, _natural


def epoch_of(ts):
    """EXIF "2026:08:12 13:47:53" or FITS "2026-08-12T13:47:53.0" -> seconds
    (any origin), as pipeline._epoch; None when it cannot be read."""
    try:
        d, tm = str(ts).replace("T", " ").split(" ", 1)
        Y, Mo, D = (int(x) for x in d.replace("-", ":").split(":"))
        h, mi, se = (float(x) for x in tm.split(":")[:3])
        return (((Y * 12 + Mo) * 31 + D) * 24 + h) * 3600.0 + mi * 60.0 + se
    except Exception:
        return None


def _exif_table(folder, raws, wd):
    """{basename: (sec, ts)} for every raw, cached in the work folder by name,
    size and modification time."""
    from .raw import read_exif
    cp = os.path.join(wd, "timeline_exif.json")
    try:
        cache = json.load(open(cp))
    except Exception:
        cache = {}
    out, dirty = {}, False
    for p in raws:
        n = os.path.basename(p)
        try:
            stt = os.stat(p)
            sig = f"{stt.st_size}:{int(stt.st_mtime)}"
        except OSError:
            continue
        c = cache.get(n)
        if c and c.get("sig") == sig:
            out[n] = (c.get("sec"), c.get("ts"))
            continue
        try:
            sec, _iso, ts = read_exif(p)
        except Exception:
            sec, ts = None, None
        cache[n] = {"sig": sig, "sec": sec, "ts": ts}
        out[n] = (sec, ts)
        dirty = True
    if dirty:
        try:
            os.makedirs(wd, exist_ok=True)
            json.dump(cache, open(cp, "w"))
        except Exception:
            pass
    return out


def build(folder, wd):
    """The timeline as JSON-able dict (see gui.html, drawTimeline)."""
    from .raw import list_raws
    raws = list_raws(folder)
    tab = _exif_table(folder, raws, wd)
    paths = [p for p in raws if tab.get(os.path.basename(p), (None, None))[0]]
    sec = {p: float(tab[os.path.basename(p)][0]) for p in paths}
    ep = {p: epoch_of(tab[os.path.basename(p)][1]) for p in paths}
    have_t = bool(paths) and all(ep[p] is not None for p in paths)
    res = find_brackets(paths, sec, ep)
    order = sorted(paths, key=lambda p: ((ep[p] if have_t else 0.0),
                                         _natural(os.path.basename(p))))
    # EXIF time stamps are whole seconds and a fast bracket puts several frames
    # in one: spread those evenly across their second, in shooting order
    t = {}
    if have_t:
        t0 = min(ep.values())
        groups = {}
        for p in order:
            groups.setdefault(ep[p], []).append(p)
        for e, g in groups.items():
            for i, p in enumerate(g):
                t[p] = (e - t0) + (i + 0.5) / len(g)
    else:
        for i, p in enumerate(order):
            t[p] = float(i)
    # what the last run did with each frame, when there was one
    status, gain, t_ref, v = {}, {}, None, None
    qp = os.path.join(wd, "quality.json")
    rp = os.path.join(wd, "report.json")
    ran = os.path.exists(qp)
    if ran:
        try:
            q = json.load(open(qp))
            seen = set()
            for tier in q.values():
                used = set(tier.get("used") or [])
                fg = tier.get("frame_gain") or {}
                fl = tier.get("frame_gain_left_out") or {}
                for n in (tier.get("scores") or {}):
                    seen.add(n)
                    if n in used:
                        status[n] = "stacked"
                    elif n in fg and fg[n] is None:
                        status[n] = "brightness"
                        if fl.get(n) is not None:
                            gain[n] = fl[n]
                    else:
                        status[n] = "dropped"
                for n, g in fg.items():
                    if g is not None:
                        gain[n] = g
            for p in paths:
                n = os.path.basename(p)
                if n not in seen:
                    status[n] = "excluded"
        except Exception:
            status, ran = {}, False
        try:
            r = json.load(open(rp))
            pf = r.get("perframe_moon") or {}
            if pf.get("t_ref_s") is not None and have_t:
                t_ref = float(pf["t_ref_s"])
                v = pf.get("v_px_s")
        except Exception:
            pass
    frames = []
    for p in order:
        n = os.path.basename(p)
        frames.append({"name": n, "sec": sec[p], "t": round(t[p], 3),
                       "bracket": int(res["of"].get(p, 0)),
                       "status": status.get(n), "gain": gain.get(n)})
    # THE SHORT EXPOSURES (2026-10-18): the frames kept apart for the
    # prominence layer belong to the shoot too, often the first frames after
    # C2. Drawn on the same time axis, own colour, not part of any bracket.
    n_short = 0
    if have_t:
        try:
            from .pipeline import find_short_dir
            sd = find_short_dir(folder)
            if sd:
                sraws = list_raws(sd)
                stab = _exif_table(sd, sraws, wd)
                t0 = min(ep.values())
                sgrp = {}
                for p in sraws:
                    sc, ts = stab.get(os.path.basename(p), (None, None))
                    e = epoch_of(ts)
                    if sc and e is not None:
                        sgrp.setdefault(e, []).append((p, float(sc)))
                for e, g in sgrp.items():
                    for i, (p, sc) in enumerate(sorted(g, key=lambda q: _natural(os.path.basename(q[0])))):
                        frames.append({"name": os.path.basename(p), "sec": sc,
                                       "t": round((e - t0) + (i + 0.5) / len(g), 3),
                                       "bracket": -1, "status": "short", "gain": None})
                        n_short += 1
                frames.sort(key=lambda f: f["t"])
        except Exception:
            pass
    # C2, mid-totality and C3 at the site, when the site is known (contacts.py)
    contacts = None
    if have_t:
        try:
            from . import observer as _ob, contacts as _ct
            site = _ob.load(folder)
            ct = _ct.shoot_contacts(site, [tab[os.path.basename(p)][1] for p in paths],
                                    contacts=_ct.contact_stamps(folder, wd))
            if ct:
                contacts = {"c2": ct["c2_s"], "mid": ct["mid_s"], "c3": ct["c3_s"],
                            "duration_s": ct["duration_s"], "clock": _ct.offset_text(ct),
                            "outside_s": ct.get("outside_s", 0.0),
                            "contact_c2": ct.get("contact_c2_s"), "contact_c3": ct.get("contact_c3_s")}
        except Exception:
            contacts = None
    return {"ok": True, "n": len(frames) - n_short, "frames": frames,
            "unread": len(raws) - len(paths), "timed": have_t,
            "repeated": bool(res["repeated"]), "n_brackets": int(res["n"]),
            "span_s": res["span_s"], "reason": res.get("reason"),
            "exposures": sorted({float(f"{s:.4g}") for s in sec.values()}),
            "ran": ran, "t_ref": t_ref, "moon_v": v, "contacts": contacts, "n_short": n_short}
