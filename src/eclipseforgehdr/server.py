"""Local web app: folder selection, pipeline run with progress, preview, export."""
from __future__ import annotations
import io, os, json, threading, traceback
import numpy as np
from flask import Flask, request, jsonify, send_file
from PIL import Image

from .brand import APP_NAME as _BRAND_NAME, LAB as _BRAND_LAB
from .pipeline import (run as run_pipeline, Progress, workdir,
                       input_fingerprint as _input_fingerprint,
                       _photometry_requested)
from .render import (Layers, render, export, defaults_for,
                     DEFAULTS as _RENDER_DEFAULTS, classify_settings)

app = Flask(__name__)
STATE = {"folder": None, "progress": None, "layers": None, "thread": None,
         # EVERY background job, not just the pipeline (0.22.78). Only the run
         # was tracked, so Clear cache during a 200 MB export ran rmtree under
         # live memory maps, and Start during one replaced STATE["progress"] so
         # the export's own poll silently followed the pipeline instead.
         "busy": None}


def _short_fingerprint(folder):
    from . import pipeline as _pl
    return _pl.short_fingerprint(folder) if hasattr(_pl, "short_fingerprint") else []


def _env_switches():
    from . import pipeline as _pl
    return _pl.env_switches() if hasattr(_pl, "env_switches") else {}


def _intra_sub():
    """Whether this build stacks the frames of a tier sub-pixel (pipeline
    INTRA_SUBPIXEL); part of the stack-cache key."""
    from . import pipeline as _pl
    return bool(getattr(_pl, "INTRA_SUBPIXEL", False))


def _busy():
    """The name of whatever background job is running, or None."""
    for k in ("thread", "busy"):
        t = STATE.get(k)
        if t is not None and getattr(t, "is_alive", lambda: False)():
            return getattr(t, "name", None) or "a job"
    return None


def _gui_html():
    here = os.path.dirname(__file__)
    return open(os.path.join(here, "gui.html"), encoding="utf-8").read()


@app.get("/")
def index():
    # CROSS-ORIGIN ISOLATION (2026-10-16): the two headers let the page use
    # shared memory, which the preview's worker pool needs (gui.html, POOL).
    # Everything the page loads comes from this server, so nothing is blocked.
    from flask import Response
    r = Response(_gui_html(), mimetype="text/html")
    r.headers["Cross-Origin-Opener-Policy"] = "same-origin"
    r.headers["Cross-Origin-Embedder-Policy"] = "require-corp"
    return r


# THE APP ICON (2026-11-08): the app's icon, also as the page's tab icon
# and beside the name in the tab bar. Same origin, so the page's cross-origin
# isolation (see index) does not block it.
@app.get("/app_icon.png")
@app.get("/favicon.ico")
def app_icon():
    fn = os.path.join(os.path.dirname(os.path.abspath(__file__)), "app_icon.png")
    if not os.path.isfile(fn):
        return "", 404
    return send_file(fn, mimetype="image/png", max_age=86400)


@app.get("/api/version")
def version():
    from . import __version__
    from .brand import APP_NAME, LAB
    return jsonify({"version": __version__, "app": APP_NAME, "lab": LAB})


# ---------- file browser ----------
#
# Pasting an absolute path is fine once you know what the app wants and awkward
# every other time -- early test runs lost time to it, one of them pasting a
# perfectly good path into a field whose Start button could never light up.
#
# This is a browser served from 127.0.0.1, so the obvious answer -- <input
# type="file"> -- is the one thing that cannot work: for security a browser
# hands JavaScript the file's CONTENT and never its path, and the pipeline needs
# a path (it reads 250 raw files off disk, it does not want them uploaded to
# itself). A native OS dialog opened by the server would give a real path, but
# Tk has to own the main thread on macOS and this server is threaded, so that
# trades a paste for a hang.
#
# So the server lists directories and the page draws the picker. It is less
# pretty than a native dialog and it behaves identically on all three platforms,
# which for a tool being tested on machines I cannot reach is the better trade.
_PICK_EXT = {
    "image": {".tif", ".tiff", ".fit", ".fits", ".fts"},
    "settings": {".json"},
    "processing": {".json"},
    "any": None,
}


def _drives():
    r"""Windows drive letters, so the picker can get above C:\Users."""
    if os.name != "nt":
        return []
    out = []
    for c in "ABCDEFGHIJKLMNOPQRSTUVWXYZ":
        d = f"{c}:\\"
        if os.path.exists(d):
            out.append(d)
    return out


@app.get("/api/browse")
def browse():
    from .raw import RAW_EXTS, list_raws
    kind = request.args.get("kind", "dir")
    path = request.args.get("path", "") or ""
    # the page sends where it would like to start and a fallback (the raw
    # folder, say, when its eclipseforge_output/ does not exist yet)
    _cands = [path, request.args.get("fallback", ""), os.path.expanduser("~")]
    path = next((os.path.abspath(os.path.expanduser(c)) for c in _cands
                 if c and os.path.isdir(os.path.expanduser(c))), os.path.expanduser("~"))
    exts = _PICK_EXT.get(kind)
    if kind == "raw":
        exts = set(RAW_EXTS) | _PICK_EXT["image"]
    dirs, files = [], []
    try:
        with os.scandir(path) as it:
            for e in it:
                if e.name.startswith("."):
                    continue
                try:
                    if e.is_dir():
                        dirs.append(e.name)
                    elif kind != "dir" and (exts is None or
                                            os.path.splitext(e.name)[1].lower() in exts):
                        files.append(e.name)
                except OSError:
                    continue
    except PermissionError:
        return jsonify({"ok": False, "error": f"no permission to read {path}"}), 400
    except OSError as e:
        return jsonify({"ok": False, "error": str(e)}), 400
    # How many raw frames each subfolder holds -- the one number that tells you
    # which of twenty folders is the bracket. Capped, so a huge tree stays fast.
    counts = {}
    if kind == "dir" and len(dirs) <= 60:
        for n in dirs:
            try:
                counts[n] = len(list_raws(os.path.join(path, n)))
            except Exception:
                counts[n] = 0
    parent = os.path.dirname(path)
    return jsonify({"ok": True, "path": path,
                    "parent": parent if parent and parent != path else None,
                    "sep": os.sep, "drives": _drives(),
                    "dirs": sorted(dirs, key=str.lower)[:2000],
                    "files": sorted(files, key=str.lower)[:2000],
                    "counts": counts,
                    "here_raws": len(list_raws(path)) if kind == "dir" else 0})


@app.post("/api/folder")
def set_folder():
    _b = _busy()
    if _b:
        return jsonify({"ok": False,
                        "error": f"{_b} is still running — wait for it to "
                                 f"finish before switching folders"}), 400
    folder = os.path.expanduser(request.json.get("path", "").strip())
    if not os.path.isdir(folder):
        return jsonify({"ok": False, "error": f"not a folder: {folder}"}), 400
    from .raw import list_raws
    raws = list_raws(folder)
    STATE["folder"] = folder
    STATE["layers"] = None
    # "cached" has to mean Start will actually reuse it, and Start also checks
    # the build. Reporting a stale-version cache as reusable was a small lie.
    from . import __version__
    _op = os.path.join(workdir(folder), "opts.json")
    cached = os.path.exists(os.path.join(workdir(folder), "geometry.json"))
    if cached:
        try:
            from . import cache_ok
            cached = cache_ok(json.load(open(_op)).get("build"))
        except Exception:
            cached = False
    # offer the conventionally-named flats subfolder, if there is one; the GUI
    # shows it and the user can clear it or point somewhere else
    from .flat import find_flat_dir
    fd = find_flat_dir(folder)
    # Same for bias/ and darks/. These are found and applied with no setting to
    # change, which is exactly why the scan has to SAY it found them: a
    # correction that happens silently is one the user cannot tell apart from a
    # correction that did not happen.
    from .dark import find_bias_dir, find_dark_dir
    bd, dd = find_bias_dir(folder), find_dark_dir(folder)
    # and the short-exposure subfolder (name has "short" or "prom"): its frames
    # feed the prominence layer only -- named here for the same reason
    try:
        from .pipeline import find_short_dir
        sd = find_short_dir(folder)
    except Exception:
        sd = None
    # the user's own frame switches (Frames panel), said here so a folder that
    # stacks fewer frames than it holds does not look like a fault
    try:
        from .frameselect import load as _sel_load
        _sl = _sel_load(folder)
        _sel_n = {"off": len(_sl["off"]), "force": len(_sl["force"])}
    except Exception:
        _sel_n = {"off": 0, "force": 0}
    # the contact-frame folder, and the folders chosen by hand in tab 1
    try:
        from .contacts import find_contact_dir as _fcd
        _cd = _fcd(folder)
    except Exception:
        _cd = None
    try:
        from .frameprod import find_phases_dir as _fpd
        _pd = _fpd(folder)
    except Exception:
        _pd = None
    try:
        from .frameselect import load as _sel_load2
        _dirs_set = dict(_sel_load2(folder).get("dirs") or {})
    except Exception:
        _dirs_set = {}
    # A finished HDR the user might want to import rather than stack: a lone
    # 16-bit TIFF sitting in the folder is almost always exactly that.
    cand = []
    try:
        for n in sorted(os.listdir(folder)):
            if os.path.splitext(n)[1].lower() in (".tif", ".tiff") and \
                    os.path.isfile(os.path.join(folder, n)):
                cand.append(n)
    except OSError:
        pass
    return jsonify({"ok": True, "folder": folder, "raw_count": len(raws),
                    "cached": cached, "flat_dir": fd,
                    "flat_count": len(list_raws(fd)) if fd else 0,
                    "bias_dir": bd,
                    "bias_count": len(list_raws(bd)) if bd else 0,
                    "dark_dir": dd,
                    "dark_count": len(list_raws(dd)) if dd else 0,
                    "short_dir": sd,
                    "short_count": len(list_raws(sd)) if sd else 0,
                    "frames_off": _sel_n["off"], "frames_forced": _sel_n["force"],
                    "contact_dir": _cd, "contact_count": len(list_raws(_cd)) if _cd else 0,
                    "phases_dir": _pd, "phases_count": len(list_raws(_pd)) if _pd else 0,
                    "dirs_set": _dirs_set,
                    "tiffs": cand[:24]})


@app.get("/api/timeline")
def timeline_api():
    """The shooting timeline under the log: one bar per frame (timeline.py)."""
    f = STATE.get("folder")
    if not f:
        return jsonify({"ok": False, "error": "no folder"})
    try:
        from .timeline import build
        return jsonify(build(f, workdir(f)))
    except Exception as e:
        return jsonify({"ok": False, "error": f"timeline: {e}"})


# ---- FRAMES PANEL (2026-11-08) ---------------------------------
# Every frame the app finds, by role, with the camera's preview, and the
# user's switch for each: auto (the app decides), off (never used) or force
# (bracket frames only: used even where an automatic check would leave it
# out). The choice lives in eclipseforge_frames.json in the image folder
# (frameselect.py); nothing is moved or deleted.

def _frame_index(folder, args):
    from . import frameselect as _fs
    idx = {}
    for role, title, forceable, d, ps in _fs.groups(folder, args.get("flat") or None,
                                                    args.get("bias") or None,
                                                    args.get("dark") or None):
        for p in ps:
            idx[_fs.key(p, folder)] = (p, role, forceable)
    return idx


@app.get("/api/frames")
def frames_api():
    f = STATE.get("folder")
    if not f:
        return jsonify({"ok": False, "error": "no folder"})
    try:
        from . import frameselect as _fs
        from .timeline import _exif_table, build as _tl_build
        wd = workdir(f)
        sel = _fs.load(f)
        tl = {}
        try:
            for fr in (_tl_build(f, wd).get("frames") or []):
                tl[fr.get("key")] = fr
        except Exception:
            pass
        # sharpness of the last run, as a fraction of the tier's sharpest frame
        sharp = {}
        try:
            for tier in json.load(open(os.path.join(wd, "quality.json"))).values():
                sc = tier.get("scores") or {}
                mx = max(sc.values()) if sc else 0
                for n, v in sc.items():
                    if mx > 0:
                        sharp[n] = round(float(v) / mx, 3)
        except Exception:
            pass
        # the short exposures the prominence layer took (promlayer.SHORT_META)
        prom = {}
        try:
            from .promlayer import SHORT_META as _SM
            for it in json.load(open(os.path.join(wd, _SM))).get("items") or []:
                n = it.get("name")
                if n and not str(n).startswith("tier "):
                    why = ("" if it.get("used") else
                           ("no match, corr %.2f" % it["corr"]) if it.get("corr") is not None else
                           ("error: " + str(it["error"])) if it.get("error") else "Moon not found")
                    prom[n] = (bool(it.get("used")), why)
        except Exception:
            pass
        groups = []
        for role, title, forceable, d, ps in _fs.groups(f, request.args.get("flat") or None,
                                                        request.args.get("bias") or None,
                                                        request.args.get("dark") or None):
            try:
                tab = _exif_table(d, ps, wd)
            except Exception:
                tab = {}
            frs = []
            # THE CONTACT FRAMES THAT SET THE CLOCK (contacts.shoot_contacts):
            # the last one before the bracket (C2) and the first one after (C3),
            # and only when the site is entered
            clock = {}
            if role == "contact":
                try:
                    from .timeline import epoch_of as _ep
                    from . import observer as _ob
                    site = _ob.load(f)
                    if site and site.get("date") and site.get("lat") is not None:
                        from .raw import list_raws as _lr
                        _lp = _lr(f)
                        _lt = _exif_table(f, _lp, wd)
                        le = [e for e in (_ep(_lt.get(os.path.basename(q), (None, None))[1]) for q in _lp)
                              if e is not None]
                        ce = {p_: _ep(tab.get(os.path.basename(p_), (None, None))[1]) for p_ in ps
                              if _fs.key(p_, f) not in sel["off"]}
                        ce = {p_: e for p_, e in ce.items() if e is not None}
                        if le:
                            bef = [p_ for p_, e in ce.items() if e < min(le)]
                            aft = [p_ for p_, e in ce.items() if e > max(le)]
                            if bef:
                                clock[max(bef, key=ce.get)] = "clock: C2"
                            if aft:
                                clock[min(aft, key=ce.get)] = "clock: C3"
                except Exception:
                    clock = {}
            for p in ps:
                k, n = _fs.key(p, f), os.path.basename(p)
                sec, ts = tab.get(n, (None, None))
                t = tl.get(k) if role in ("lights", "short") else None
                if role == "short" and n in prom and k not in sel["off"]:
                    status, note = ("prom", "") if prom[n][0] else ("promunused", prom[n][1])
                elif role == "contact" and p in clock:
                    status, note = "clock", clock[p]
                else:
                    status, note = ((t or {}).get("status") if role == "lights" else None), ""
                frs.append({"key": k, "name": n, "sec": sec, "ts": ts,
                            "state": ("off" if k in sel["off"] else
                                      "force" if (forceable and k in sel["force"]) else "auto"),
                            "status": status, "note": note,
                            "gain": (t or {}).get("gain") if role == "lights" else None,
                            "sharp": sharp.get(n) if role == "lights" else None})
            groups.append({"role": role, "title": title, "forceable": forceable,
                           "dir": _fs.key(d, f), "frames": frs})
        return jsonify({"ok": True, "groups": groups, "file": _fs.SEL_FILE,
                        "n_off": len(sel["off"]), "n_force": len(sel["force"])})
    except Exception as e:
        return jsonify({"ok": False, "error": f"frames: {e}"})


@app.get("/api/framepreview")
def frame_preview():
    f = STATE.get("folder")
    if not f:
        return jsonify({"ok": False, "error": "no folder"}), 400
    k = request.args.get("k", "")
    kind = request.args.get("s", "thumb")
    if kind not in ("thumb", "big", "crop"):
        kind = "thumb"
    try:
        from . import frameselect as _fs
        hit = _frame_index(f, request.args).get(k)
        if not hit:
            return jsonify({"ok": False, "error": "not a frame of this folder"}), 404
        out = _fs.preview(hit[0], kind, os.path.join(workdir(f), "framepreviews"),
                          mode=request.args.get("m", "strong"))
        return send_file(out, mimetype="image/jpeg", max_age=3600)
    except Exception as e:
        return jsonify({"ok": False, "error": str(e)}), 500


@app.post("/api/frames")
def frames_set():
    """{"set": {key: "auto" | "off" | "force"}} -> saved; the next Start
    re-stacks when the selection changed what it would stack."""
    f = STATE.get("folder")
    if not f:
        return jsonify({"ok": False, "error": "no folder"}), 400
    _b = _busy()
    if _b:
        return jsonify({"ok": False, "error": f"{_b} is still running — change the frames when it has finished"}), 400
    try:
        from . import frameselect as _fs
        data = request.json or {}
        idx = _frame_index(f, data)
        sel = _fs.load(f)
        off, force = set(sel["off"]), set(sel["force"])
        for k, v in (data.get("set") or {}).items():
            if k not in idx:
                continue
            off.discard(k)
            force.discard(k)
            if v == "off":
                off.add(k)
            elif v == "force" and idx[k][2]:
                force.add(k)
        _fs.save(f, off, force)
        return jsonify({"ok": True, "n_off": len(off), "n_force": len(force - off)})
    except Exception as e:
        return jsonify({"ok": False, "error": str(e)}), 500


@app.get("/api/checklist")
def checklist():
    """TAB 1: what the folder holds and what is missing, as a short list of
    findings, each with a level (ok | warn | missing | info), one line of text
    and where to act. Only reads what the other views already read: the time
    stamps of the frames (cached), the folders found, the site."""
    f = STATE.get("folder")
    if not f:
        return jsonify({"ok": False, "error": "no folder"})
    items = []

    def add(level, title, text, act=None):
        items.append({"level": level, "title": title, "text": text, "act": act})
    try:
        from .raw import list_raws, read_exif
        from .timeline import build as _tl_build, _exif_table
        from .dark import find_dark_dir, find_bias_dir, DARK_EXTRAP_WARN, DARK_EXTRAP_MAX, DARK_GROUP_TOL
        from .flat import find_flat_dir
        from .contacts import find_contact_dir
        from .pipeline import find_short_dir
        from . import frameselect as _fs, observer as _ob
        wd = workdir(f)
        lights = list_raws(f)
        tl = {}
        try:
            tl = _tl_build(f, wd)
        except Exception:
            pass
        exps = sorted({float(x) for x in (tl.get("exposures") or [])})
        ex_txt = lambda t: (f"{t:g} s" if t >= 1 else f"1/{round(1 / t):g} s")
        # 1. the bracket
        if len(lights) < 3:
            add("missing", "Corona frames", f"{len(lights)} frame(s) found; at least 3 are needed", {"field": "folder"})
        else:
            br = (f"{tl.get('n_brackets')} repeated brackets" if tl.get("repeated") else "one bracket")
            add("ok", "Corona frames", f"{len(lights)} frames, {len(exps)} exposure(s)"
                + (f" ({ex_txt(exps[0])} to {ex_txt(exps[-1])})" if len(exps) > 1 else
                   (f" ({ex_txt(exps[0])})" if exps else "")) + f", {br}")
        if lights and not tl.get("timed", True):
            add("warn", "Time stamps", "some frames carry no time stamp: no timeline, no clock check, "
                "no Moon-motion correction from the shooting times")
        if tl.get("unread"):
            add("warn", "Exposure times", f"{tl['unread']} frame(s) without an exposure time are left out")
        # 2. darks: which exposures they reach (the run's rule: a dark scaled up
        # more than DARK_EXTRAP_MAX times is not used, above WARN it is noisy)
        dd, bd = find_dark_dir(f), find_bias_dir(f)
        light_iso = None
        try:
            light_iso = int(read_exif(lights[len(lights) // 2])[1] or 0) or None
        except Exception:
            pass
        if dd:
            dps = list_raws(dd)
            dt = _exif_table(dd, dps, wd)
            dsec = sorted(float(v[0]) for v in dt.values() if v and v[0])
            groups = []
            for x in dsec:
                if groups and abs(x - groups[-1]) <= DARK_GROUP_TOL * groups[-1]:
                    continue
                groups.append(x)
            dref = max(dsec) if dsec else 0.0
            over = [t for t in exps if dref and t / dref > DARK_EXTRAP_MAX]
            noisy = [t for t in exps if dref and DARK_EXTRAP_WARN < t / dref <= DARK_EXTRAP_MAX]
            iso_bad = 0
            if light_iso:
                for q in dps[:6]:
                    try:
                        i = int(read_exif(q)[1] or 0)
                        if i and i != light_iso:
                            iso_bad += 1
                    except Exception:
                        pass
            txt = (f"{len(dps)} in {os.path.basename(dd)}/ at "
                   + ", ".join(ex_txt(g) for g in groups[:6]) + ("…" if len(groups) > 6 else ""))
            if over:
                add("warn", "Darks", txt + "; too short for " + ", ".join(ex_txt(t) for t in over)
                    + f" (more than {DARK_EXTRAP_MAX:g}x): those get "
                    + ("the bias only" if bd else "no dark correction"), {"field": "dir_darks"})
            elif noisy:
                add("warn", "Darks", txt + "; scaled up more than "
                    f"{DARK_EXTRAP_WARN:g}x for " + ", ".join(ex_txt(t) for t in noisy)
                    + " (their noise is scaled with them)", {"field": "dir_darks"})
            else:
                add("ok", "Darks", txt + "; they reach every exposure", {"field": "dir_darks"})
            if iso_bad:
                add("warn", "Darks ISO", f"{iso_bad} of the first {min(6, len(dps))} darks were not shot at "
                    f"ISO {light_iso} like the light frames; those are left out", {"field": "dir_darks"})
        else:
            add("info", "Darks", "none: no dark-current correction (hot pixels are still repaired). "
                "Choose a folder if they are elsewhere", {"field": "dir_darks"})
        add("ok" if bd else "info", "Bias",
            f"{len(list_raws(bd))} in {os.path.basename(bd)}/" if bd else
            "none: the black level comes from the camera" + (", and the darks carry it" if dd else ""),
            {"field": "dir_bias"})
        # 3. flats, from the field if one is typed there
        fq = (request.args.get("flat") or "").strip()
        fd = None if fq.lower() in ("off", "none") else (fq if fq and os.path.isdir(fq) else find_flat_dir(f))
        add("ok" if fd else "warn", "Flats",
            f"{len(list_raws(fd))} in {os.path.basename(fd)}/" if fd else
            "none: vignetting and dust are not corrected", {"field": "flatDir"})
        # 4. short exposures and contact frames
        sd = find_short_dir(f)
        add("ok" if sd else "info", "Short exposures",
            f"{len(list_raws(sd))} in {os.path.basename(sd)}/, for the prominence layer" if sd else
            "none: the prominences come from the bracket's shortest exposures", {"field": "dir_short"})
        cd = find_contact_dir(f)
        site = None
        try:
            site = _ob.load(f)
        except Exception:
            pass
        has_site = bool(site and site.get("ut") and site.get("lat") is not None)
        if cd:
            add("ok" if has_site else "info", "Contact frames",
                f"{len(list_raws(cd))} in {os.path.basename(cd)}/"
                + ("; they set the camera clock" if has_site else
                   "; enter the observing site (tab 2) so they can set the camera clock"),
                {"field": "dir_contact"} if has_site else {"tab": 2})
        else:
            add("info", "Contact frames", "none: the camera clock is taken as it is", {"field": "dir_contact"})
        # 5. the site and the user's own frame choices
        add("ok" if has_site else "info", "Observing site",
            "entered" if has_site else
            "not entered: the star match is less exact and the polar map has no solar north", {"tab": 2})
        sel = _fs.load(f)
        if sel["off"] or sel["force"]:
            add("info", "Your frame choices", f"{len(sel['off'])} switched off, {len(sel['force'])} forced in")
        return jsonify({"ok": True, "items": items})
    except Exception as e:
        return jsonify({"ok": False, "error": f"checklist: {e}", "items": items})


@app.post("/api/framedirs")
def frame_dirs():
    """A darks, bias, short-exposure or contact-frame folder chosen by hand:
    {"role": "darks", "value": "/path" | "" (automatic) | "off"}."""
    f = STATE.get("folder")
    if not f:
        return jsonify({"ok": False, "error": "load a folder first"}), 400
    _b = _busy()
    if _b:
        return jsonify({"ok": False, "error": f"{_b} is still running"}), 400
    from . import frameselect as _fs
    data = request.json or {}
    role = str(data.get("role") or "")
    if role not in _fs.DIR_ROLES:
        return jsonify({"ok": False, "error": f"unknown role {role!r}"}), 400
    v = str(data.get("value") or "").strip()
    if v and v.lower() not in ("auto", "off", "none", "-"):
        p = os.path.abspath(os.path.expanduser(v))
        if not os.path.isdir(p):
            return jsonify({"ok": False, "error": f"not a folder: {v}"}), 400
        if os.path.abspath(p) == os.path.abspath(f):
            return jsonify({"ok": False, "error": "that is the image folder itself; choose its subfolder"}), 400
        try:
            from .pipeline import holds_lights as _hl
            if _hl(f, p):
                return jsonify({"ok": False, "error": "that folder holds the light frames"}), 400
        except Exception:
            pass
    _fs.set_dir(f, role, v)
    return jsonify({"ok": True})


@app.post("/api/clearcache")
def clear_cache():
    """Throw away everything the app computed for this folder.

    There was no way to do this from the GUI. The log line said "Start with
    force to redo", but the page never sent force -- that flag was reachable
    only from the command line, so from the GUI a stale cache was permanent.
    Deleting opts.json by hand is not enough either: the next run writes it
    again, so the run after that reports the cache as valid, which is exactly
    what it looked like from the outside.

    Scope: only the app's own `.eclipseforgehdr` subfolder. The raws, the
    flats, and anything exported into the folder are not touched. The name is
    checked before anything is removed -- this deletes a directory tree, and
    the one guard that matters is that it is OUR directory.
    """
    if not request.is_json:            # see /api/run
        return jsonify({"ok": False, "error": "expected a JSON body"}), 415
    _b = _busy()
    if _b:
        return jsonify({"ok": False,
                        "error": f"{_b} is still running — wait for it to "
                                 f"finish before clearing the cache"}), 400
    folder = STATE["folder"]
    if not folder:
        _ip = os.path.expanduser(str((request.json or {}).get("importPath", "")).strip()) \
            if request.is_json else ""
        if _ip and os.path.isfile(_ip):
            folder = os.path.dirname(os.path.abspath(_ip))
    if not folder or not os.path.isdir(folder):
        return jsonify({"ok": False, "error": "choose a folder first"}), 400
    from .brand import WORKDIR_NAME
    wd = os.path.join(folder, WORKDIR_NAME)
    if os.path.basename(wd) != WORKDIR_NAME:
        return jsonify({"ok": False, "error": "refusing to clear that path"}), 400
    if not os.path.isdir(wd):
        STATE["layers"] = None
        return jsonify({"ok": True, "files": 0, "bytes": 0, "path": wd,
                        "note": "nothing cached"})
    n = size = 0
    for root, _dirs, files in os.walk(wd):
        for f in files:
            try:
                size += os.path.getsize(os.path.join(root, f))
                n += 1
            except OSError:
                pass
    # Drop the preview's handles first: Layers holds memory maps into these
    # files, and on Windows an open map makes the delete fail outright.
    STATE["layers"] = None
    import shutil, gc
    gc.collect()
    try:
        shutil.rmtree(wd)
    except OSError as e:
        return jsonify({"ok": False, "error": f"could not clear: {e}"}), 400
    return jsonify({"ok": True, "files": n, "bytes": size, "path": wd})



def _layers_stale(wd):
    """Was this work directory's detail-layer set built by an older recipe?
    (detail.LAYER_BUILD, 0.23.9). The stack is not the question -- only what
    was derived from it -- so this is what lets a layer change reach a cached
    stack in minutes instead of a re-stack."""
    from .detail import LAYER_BUILD
    try:
        _rj = json.load(open(os.path.join(wd, "report.json")))
        have = int((_rj.get("layer_recipe") or {}).get("build", 0))
    except Exception:
        have = 0
    if have >= LAYER_BUILD:
        # the stars are lifted out of the corona while the layers are built
        # (stars.remove in detail.build_layers), so a newer star recipe makes
        # the layers stale too (audit 2026-10-10: stars.json build 6 sat over
        # build-5 patches)
        try:
            from .stars import STARS_BUILD, STARS_FILE
            _sb = int(json.load(open(os.path.join(wd, STARS_FILE))).get("build", 0))
            if _sb < STARS_BUILD:
                return True, have, LAYER_BUILD
        except Exception:
            pass
    return have < LAYER_BUILD, have, LAYER_BUILD


def _star_summary(ly):
    """What stars.json says, for the page: the catalogue solution and the
    lifted stars with their catalogue identity. Never raises."""
    try:
        m = json.load(open(os.path.join(ly.wd, "stars.json")))
        def _c(c, k, d=None):
            return (c.get("catalogue") or {}).get(k, d)
        st = [dict(y=c.get("y"), x=c.get("x"), mag=_c(c, "mag"), name=_c(c, "name", "") or "",
                   label=_c(c, "label", "") or "", hd=_c(c, "hd", 0), hr=_c(c, "hr", 0),
                   sao=_c(c, "sao", 0), hip=_c(c, "hip", 0),
                   ra=_c(c, "ra"), dec=_c(c, "dec"),
                   r=c.get("r"), z=c.get("z"), confirmed=c.get("confirmed", ""))
              for c in m.get("stars", []) if c.get("class") == "star"]
        return {"catalogue": m.get("catalogue") or {}, "stars": st,
                "streamers": m.get("streamers") or [], "shape": m.get("shape"), "moon": m.get("moon"),
                "unconfirmed": sum(1 for c in m.get("stars", []) if c.get("class") == "unconfirmed")}
    except Exception:
        return None


def _build_stars(wd, prog):
    """STAR LIST (0.25.1): stars found in the merge, written to
    stars.json; seconds, rebuilt when missing or older than the merge. Never
    stops a run."""
    try:
        from . import stars as _st
        if _st.needs_build(wd):
            prog.log("finding stars...", None)
            _st.build(wd, prog)
    except Exception as _e:
        if getattr(prog, "cancelled", False):
            raise
        prog.log(f"[warn] star search skipped ({_e})", None)


def _site_changed(wd):
    """The site file is newer than the star list, or was cleared after it."""
    try:
        from . import observer as _ob
        sp = _ob.site_path(os.path.dirname(os.path.abspath(wd)))
        sj = os.path.join(wd, "stars.json")
        if not os.path.exists(sj):
            return os.path.exists(os.path.join(wd, "stars_layer.npz"))
        return os.path.exists(sp) and os.path.getmtime(sp) > os.path.getmtime(sj)
    except Exception:
        return False


def _site_reply(folder):
    from . import observer as _ob
    s = _ob.load(folder) or {}
    cam = None
    try:
        cam = json.load(open(os.path.join(workdir(folder), "report.json"))).get("shot_first")
    except Exception:
        pass
    sun = None
    try:
        aa = _ob.sun_altaz(s)
        if aa is not None:
            Y, M, D, h = s["ut"]
            hh = int(h); mm = int((h - hh) * 60); ss = int(round(((h - hh) * 60 - mm) * 60))
            if ss == 60:
                mm, ss = mm + 1, 0
            sun = {"alt": round(aa[0], 2), "az": round(aa[1], 2),
                   "ut": f"{Y:04d}-{M:02d}-{D:02d} {hh:02d}:{mm:02d}:{ss:02d} UT"}
    except Exception:
        pass
    keep = ("lat_text", "lon_text", "height_m", "date", "time", "utc_offset_text", "lat", "lon")
    # the position against the eclipse's path of totality (eclipses.py), when
    # the date is a bundled eclipse: total there or not, and the local times
    path = None
    try:
        from . import eclipses as _ec
        r = _ec.check_site(s)
        if r is not None:
            off = s.get("utc_offset")
            path = {"date": r["date"], "total": bool(r["total"]), "duration_s": round(r["duration_s"], 1),
                    "sun_alt": round(r["sun_alt"], 1), "mid_ut": _ec.hms(r["mid_ut"]),
                    "mid_local": _ec.hms(r["mid_ut"] + off) if off is not None else None}
            if s.get("ut") and r["total"]:
                path["time_diff_min"] = round((s["ut"][3] - r["mid_ut"]) * 60.0, 1)
    except Exception:
        path = None
    return {"ok": True, "site": {k: s.get(k) for k in keep if k in s}, "camera_time": cam,
            "sun": sun, "path": path, "file": _ob.SITE_FILE}


@app.get("/api/eclipses")
def get_eclipses():
    """The bundled total eclipses and the places in their paths (eclipses.py)."""
    from . import eclipses as _ec
    return jsonify({"ok": True, "eclipses": [{"date": d, "name": e["name"], "cities": _ec.cities(d)}
                                             for d, e in sorted(_ec.ECLIPSES.items())]})


# THE RECIPE CHECK (2026-11-08). Tab 2, above Start: the processing
# options as set in the page against what the folder holds, so a setting that
# cannot help with this data, or works against it, shows before the run and
# not in the picture afterwards. Only settled facts of how each option works
# are checked here (what it needs, what it does with this data); nothing is
# judged on a threshold that was not measured. Never blocks Start: each
# finding carries the suggested value and the page offers a button for it.
@app.post("/api/recipecheck")
def recipecheck():
    import math as _m
    f = STATE.get("folder")
    if not f:
        return jsonify({"ok": False, "error": "no folder"})
    data = request.json if request.is_json else {}
    o = {k: str(v) for k, v in (data.get("processing") or {}).items()}
    rec = data.get("recipe") or None
    items = []

    def add(level, title, text, fix=None):
        items.append({"level": level, "title": title, "text": text, "fix": fix})

    def fix(key, value, label):
        return {"id": _PROCESSING_KEYS[key][1], "value": value, "label": label}
    try:
        from .raw import list_raws, read_camera_info, RAW_EXTS
        from .timeline import build as _tl_build
        from .dark import find_dark_dir
        from . import observer as _ob
        wd = workdir(f)
        lights = list_raws(f)
        tl = {}
        try:
            tl = _tl_build(f, wd)
        except Exception:
            pass
        fr = [x for x in (tl.get("frames") or []) if x.get("bracket", 0) >= 0 and x.get("status") != "off"]
        by = {}
        for x in fr:
            by.setdefault(float(f"{x['sec']:.4g}"), []).append(x)
        exps = sorted(by)
        ex_txt = lambda t: (f"{t:g} s" if t >= 1 else f"1/{round(1 / t):g} s")
        counts = {e: len(by[e]) for e in exps}
        nmin = min(counts.values()) if counts else 0
        steps = [_m.log2(exps[i + 1] / exps[i]) for i in range(len(exps) - 1)]
        maxstep = max(steps) if steps else 0.0
        simple = o.get("simple") == "True"
        frames = o.get("frames", "all")

        # 1. the simple stack: the fallback merge
        if simple and len(exps) > 1:
            add("warn", "Simple stack", "the fallback merge (no ladder, no edge feathering) is on; "
                f"this set has {len(exps)} exposures with steps up to {maxstep:.1f} EV, where it can show "
                "rings. Use it when the normal stack fails", fix("simple", "False", "Use the normal stack"))
        elif not simple:
            add("ok", "Stack", "normal stack")
        # 2. clipped mean: the run applies it only to an exposure with 3 or
        # more frames (pipeline: _clip_stack and _n_added >= 3), the others
        # are a plain mean
        if not simple and o.get("stack_combine") == "clip" and lights and counts:
            used = {e: (1 if frames == "best" else (max(1, c // 2) if frames == "best50" else c))
                    for e, c in counts.items()}
            thin = [e for e in exps if used[e] < 3]
            if frames == "best":
                add("info", "Tier combine", "Best only keeps one frame per exposure, so Clipped mean has "
                    "nothing to compare and every exposure is a plain mean")
            elif len(thin) == len(exps):
                add("info", "Tier combine", f"Clipped mean needs 3 frames of one exposure to find a particle "
                    f"strike; no exposure here has that many (at most {max(used.values())}), so all are a "
                    "plain mean. No harm; strikes stay unless the frames are repeated")
            elif thin:
                add("info", "Tier combine", f"Clipped mean runs on the exposures with 3 or more frames; "
                    + ", ".join(ex_txt(e) for e in thin[:6]) + (" …" if len(thin) > 6 else "")
                    + " have fewer and are a plain mean")
            else:
                add("ok", "Tier combine", "clipped mean")
        # 3. fewer frames than there are
        if not simple and frames in ("best", "best50") and nmin >= 2:
            k = 1 if frames == "best" else max(1, nmin // 2)
            add("info", "Frames", f"{'Best only' if frames == 'best' else 'Best half'} keeps {k} of "
                f"{nmin}+ frames per exposure: noise about {_m.sqrt(nmin / k):.1f}x that of all frames. "
                "If the run's sharpness spread is near 1.0, averaging all costs no detail",
                fix("frames", "all", "Use all frames"))
        # 4. HDR per bracket needs repeated brackets
        if not simple and o.get("bracket_mode") == "hdr":
            if not tl.get("repeated"):
                add("info", "Brackets", "HDR per bracket is meant for repeated brackets; this set is one "
                    "bracket, so the run merges it as usual (Stack per exposure)",
                    fix("bracket_mode", "tiers", "Use Stack per exposure"))
            else:
                add("info", "Brackets", f"HDR per bracket on {tl.get('n_brackets')} brackets: experimental, "
                    "meant for long eclipses where sky or Moon position change between brackets")
        # 5. the Moon's drift inside one exposure's frames, against Moon Align
        site = None
        try:
            site = _ob.load(f)
        except Exception:
            pass
        circ = None
        try:
            from . import eclipses as _ec
            circ = _ec.check_site(site)
        except Exception:
            circ = None
        lock = o.get("intra_lock", "corona")
        if not simple and lock in ("moon", "mixed") and tl.get("timed") and exps:
            spread = {e: (max(x["t"] for x in by[e]) - min(x["t"] for x in by[e])) for e in exps}
            e_w = max(spread, key=spread.get)
            sp = spread[e_w]
            rate_r = None
            try:
                if site and site.get("lat") is not None and site.get("date"):
                    rate_r = _ec.moon_rate(str(site["date"])[:10], float(site["lat"]), float(site["lon"]),
                                           float(site.get("height_m") or 0.0))
            except Exception:
                rate_r = None
            R = None
            try:
                R = float(json.load(open(os.path.join(wd, "geometry.json")))["R"])
            except Exception:
                R = None
            if R and sp > 0:
                px = (rate_r * R if rate_r else 0.5 / (1920.0 / (2.0 * R))) * sp
                lvl = "warn" if (px > 1.0 and lock == "moon") else "info"
                add(lvl, "Intra-tier lock", f"{'Moon Align' if lock == 'moon' else 'Mixed'}: the frames of "
                    f"{ex_txt(e_w)} span {sp:.0f} s, in which the Moon moves about {px:.1f} px against the corona"
                    + ("" if rate_r else " (typical rate; enter the site for this eclipse's)")
                    + (("; locked on the Moon, the corona is smeared by that much" if lock == "moon" else
                        "; exposures with little corona to lock on (usually the fast ones, whose inner corona "
                        "matters most) go on the Moon and smear by up to that much") if px > 1.0 else
                       "; below a pixel, no harm"),
                    fix("intra_lock", "corona", "Use Corona Align") if px > 1.0 else None)
            elif sp > 0:
                add("info", "Intra-tier lock", f"{'Moon Align' if lock == 'moon' else 'Mixed'}: the frames of "
                    f"{ex_txt(e_w)} span {sp:.0f} s; the Moon moves about 0.5 arcsec/s against the corona, "
                    "so locked on the Moon the corona smears by that much (in pixels after the first run)",
                    fix("intra_lock", "corona", "Use Corona Align"))
        elif not simple and lock in ("moon", "mixed"):
            add("info", "Intra-tier lock", f"{'Moon Align' if lock == 'moon' else 'Mixed'}: without time stamps "
                "the Moon's drift inside an exposure cannot be worked out. Corona Align is safe when an "
                "exposure's frames are spread over totality", fix("intra_lock", "corona", "Use Corona Align"))
        elif not simple:
            add("ok", "Intra-tier lock", "Corona Align")
        # 6. colour planes at a low Sun
        if not simple and o.get("colour_planes") == "off":
            if circ and circ.get("sun_alt") is not None and circ["sun_alt"] < 30:
                add("warn", "Colour planes", f"the Sun was {circ['sun_alt']:.0f}\u00b0 high: the air spreads "
                    "red and blue apart (atmospheric dispersion). Align to green measures it and moves them "
                    "only when the offset is 0.2 px or more", fix("colour_planes", "auto", "Align to green"))
            else:
                add("info", "Colour planes", "As recorded: red and blue are not checked against green. "
                    "Align to green is safe at any Sun height (it moves nothing below 0.2 px)",
                    fix("colour_planes", "auto", "Align to green"))
        elif not simple:
            add("ok", "Colour planes", "aligned to green where measured")
        # 7. alignment off on camera raws
        am = o.get("align_mode", "auto")
        if lights and am == "keep":
            raw_in = all(os.path.splitext(p)[1].lower() in RAW_EXTS for p in lights)
            add("warn" if raw_in else "info", "Alignment",
                ("Keep frames as they are: these are camera raws, which no other program has registered, so "
                 "nothing aligns them and any shift between frames stays in the stack") if raw_in else
                ("Keep frames as they are: right only for frames another program already registered; "
                 "Automatic recognises those by their registration border"),
                fix("align_mode", "auto", "Use Automatic"))
        else:
            add("ok", "Alignment", "automatic" if am == "auto" else "always align")
        # 8. hot pixels: no darks and no repair
        if o.get("despeckle") == "False" and not find_dark_dir(f):
            add("warn", "Hot pixels", "no darks and Fix hot pixels is off: hot photosites stay in the picture "
                "as coloured dots", {"id": "despeckle", "value": "True", "label": "Fix hot pixels"})
        # 9. edge weight with big exposure steps
        if not simple and o.get("feather") == "taper" and maxstep >= 2.0:
            add("info", "Merge weight", f"Exact edge with steps up to {maxstep:.1f} EV between exposures: "
                "radiometrically right, but if rings show just outside the limb, Blended edge hides them",
                fix("feather", "plain", "Use Blended edge"))
        # 10. photometric chain over many tiers
        if not simple and o.get("photo_solve") == "chain" and len(exps) >= 4:
            add("info", "Photometric solve", f"Chain over {len(exps)} exposures: one bad link shifts every "
                "exposure after it. Network solves all links together", fix("photo_solve", "network", "Use Network"))
        # 11. the recipe's own data
        if rec and rec.get("camera") is not None:
            cam = {}
            try:
                cam = read_camera_info(lights[len(lights) // 2]) if lights else {}
            except Exception:
                cam = {}
            rc = rec["camera"]
            rname = rc.get("camera") if isinstance(rc, dict) else str(rc)
            rfl = rc.get("focal_mm") if isinstance(rc, dict) else None
            diff = []
            if rname and cam.get("camera") and rname.strip() != cam["camera"].strip():
                diff.append(f"camera {rname} vs {cam['camera']}")
            if rfl and cam.get("focal_mm") and abs(float(rfl) - float(cam["focal_mm"])) > 1:
                diff.append(f"{float(rfl):g} mm vs {float(cam['focal_mm']):g} mm")
            if diff:
                nd = [k for k, (dv, _) in _PROCESSING_KEYS.items() if o.get(k) is not None and o.get(k) != dv]
                add("info", "Recipe", f"{rec.get('file') or 'the recipe'} was saved for other data ("
                    + "; ".join(diff) + ")" + (f"; its non-default settings: {', '.join(nd)}" if nd else ""))
        warn = sum(1 for x in items if x["level"] == "warn")
        return jsonify({"ok": True, "items": items, "warn": warn,
                        "exposures": [ex_txt(e) for e in exps], "per_exposure_min": nmin})
    except Exception as e:
        return jsonify({"ok": False, "error": f"recipe check: {e}", "items": items})


# MAP TILES THROUGH THE APP'S SERVER (2026-11-08). The map picker asked the
# tile servers straight from the page and on the Mac the map stayed empty.
# Fetched here instead: one identifying User-Agent (the tile servers' usage
# policy asks for that), at most two downloads at a time, CA certificates
# from certifi where present (a frozen Python may have none of its own), and
# every tile kept on disk -- a map once seen also works offline in the field.
# When a tile cannot be had the reason is kept and the page can say it.
_TILE_SRC = {"osm": "https://tile.openstreetmap.org/{z}/{x}/{y}.png",
             "carto": "https://{s}.basemaps.cartocdn.com/dark_all/{z}/{x}/{y}.png"}
_TILE_DIR = os.path.join(os.path.expanduser("~"), ".cache", "eclipseforgehdr", "tiles")
_TILE_FRESH_S = 30 * 86400
_TILE_SEM = threading.Semaphore(2)
_TILE_ERR = {}
_TILE_SSL = None


def _tile_ssl():
    global _TILE_SSL
    if _TILE_SSL is None:
        import ssl
        try:
            import certifi
            _TILE_SSL = ssl.create_default_context(cafile=certifi.where())
        except Exception:
            _TILE_SSL = ssl.create_default_context()
    return _TILE_SSL


@app.get("/api/tile/<src>/<int:z>/<int:x>/<int:y>.png")
def map_tile(src, z, x, y):
    import time as _t
    import urllib.request
    from . import __version__
    if src not in _TILE_SRC or not (0 <= z <= 17) or not (0 <= x < 2 ** z) or not (0 <= y < 2 ** z):
        return jsonify({"ok": False, "error": "no such tile"}), 404
    fn = os.path.join(_TILE_DIR, src, str(z), str(x), f"{y}.png")
    if os.path.isfile(fn) and _t.time() - os.path.getmtime(fn) < _TILE_FRESH_S:
        return send_file(fn, mimetype="image/png", max_age=86400)
    url = _TILE_SRC[src].format(s="abcd"[(x + y) % 4], z=z, x=x, y=y)
    req = urllib.request.Request(url, headers={
        "User-Agent": f"EclipseForgeHDR/{__version__} (eclipse image processing; site picker map)"})
    try:
        with _TILE_SEM:
            with urllib.request.urlopen(req, timeout=10, context=_tile_ssl()) as r:
                body = r.read()
                ctype = r.headers.get("Content-Type", "")
        if not body or "image" not in ctype:
            raise ValueError(f"the tile server sent {ctype or 'nothing'} instead of an image")
        os.makedirs(os.path.dirname(fn), exist_ok=True)
        tmp = fn + ".part"
        with open(tmp, "wb") as fh:
            fh.write(body)
        os.replace(tmp, fn)
        _TILE_ERR.pop(src, None)
        return send_file(io.BytesIO(body), mimetype="image/png", max_age=86400)
    except Exception as e:
        msg = str(getattr(e, "reason", "") or e)
        if hasattr(e, "code"):
            msg = f"HTTP {e.code} {msg}"
        _TILE_ERR[src] = msg[:200]
        if os.path.isfile(fn):                     # stale but better than nothing
            return send_file(fn, mimetype="image/png", max_age=600)
        return jsonify({"ok": False, "error": msg}), 502


@app.get("/api/tilestatus")
def map_tile_status():
    return jsonify({"ok": True, "errors": dict(_TILE_ERR), "cache": _TILE_DIR})


@app.get("/api/eclipsepath")
def get_eclipsepath():
    """The path of totality for the map picker: centre line and the two
    limits as [lat, lon] lists (eclipses.path)."""
    from . import eclipses as _ec
    d = (request.args.get("date") or "")[:10]
    if d not in _ec.ECLIPSES:
        return jsonify({"ok": False, "error": "not a bundled eclipse"}), 404
    return jsonify({"ok": True, **_ec.path(d)})


@app.get("/api/eclipselocal")
def get_eclipselocal():
    """Local circumstances at a clicked place: total or not, duration, the
    contacts in UT and local time, the Sun's height. The UTC offset comes
    from the nearest listed place in the path (its time zone on that date)
    when one is within about 300 km, else from the longitude (marked as a
    guess)."""
    import math as _m
    from . import eclipses as _ec
    d = (request.args.get("date") or "")[:10]
    if d not in _ec.ECLIPSES:
        return jsonify({"ok": False, "error": "not a bundled eclipse"}), 404
    try:
        lat = float(request.args.get("lat")); lon = float(request.args.get("lon"))
        h = float(request.args.get("h") or 0.0)
    except (TypeError, ValueError):
        return jsonify({"ok": False, "error": "bad position"}), 400
    lon = (lon + 540.0) % 360.0 - 180.0
    if not (-90.0 <= lat <= 90.0):
        return jsonify({"ok": False, "error": "bad latitude"}), 400
    r = _ec.local(d, lat, lon, h)
    off, src, best = None, "longitude", None
    for c in _ec.cities(d):
        p1, p2 = _m.radians(lat), _m.radians(c["lat"])
        dl = _m.radians(c["lon"] - lon)
        km = 6371.0 * _m.acos(max(-1.0, min(1.0, _m.sin(p1) * _m.sin(p2) + _m.cos(p1) * _m.cos(p2) * _m.cos(dl))))
        if best is None or km < best[0]:
            best = (km, c)
    if best and best[0] < 300.0:
        o = _ec.utc_offset_h(best[1].get("tz"), d, r["mid_ut"]) if best[1].get("tz") else best[1].get("utc_offset")
        if o is not None:
            off, src = o, "near " + best[1]["name"] + (", " + best[1]["country"] if best[1].get("country") else "")
    if off is None:
        off = float(round(lon / 15.0))

    def lt(x):
        return _ec.hms(x + off) if x is not None else None
    out = {"ok": True, "date": d, "lat": round(lat, 5), "lon": round(lon, 5),
           "total": bool(r.get("total")), "duration_s": round(r.get("duration_s") or 0.0, 1),
           "sun_alt": round(r["sun_alt"], 1), "magnitude": round(r.get("magnitude") or 0.0, 4),
           "mid_ut": _ec.hms(r["mid_ut"]), "mid_local": lt(r["mid_ut"]),
           "c2_ut": _ec.hms(r["c2_ut"]) if r.get("c2_ut") is not None else None,
           "c3_ut": _ec.hms(r["c3_ut"]) if r.get("c3_ut") is not None else None,
           "c2_local": lt(r.get("c2_ut")), "c3_local": lt(r.get("c3_ut")),
           "utc_offset": off, "offset_from": src}
    return jsonify(out)


@app.get("/api/site")
def get_site():
    folder = STATE.get("folder")
    if not folder:
        return jsonify({"ok": False, "error": "no folder"}), 400
    return jsonify(_site_reply(folder))


@app.post("/api/site")
def set_site():
    if not request.is_json:
        return jsonify({"ok": False, "error": "expected a JSON body"}), 415
    folder = STATE.get("folder")
    if not folder:
        return jsonify({"ok": False, "error": "choose the folder first"}), 400
    from . import observer as _ob
    data = request.json or {}
    try:
        if data.get("clear"):
            if os.path.exists(_ob.site_path(folder)):
                os.remove(_ob.site_path(folder))
            sj = os.path.join(workdir(folder), "stars.json")
            if os.path.exists(sj):
                os.remove(sj)
        else:
            _ob.save(folder, data)
    except ValueError as e:
        return jsonify({"ok": False, "error": str(e)}), 400
    except OSError as e:
        return jsonify({"ok": False, "error": f"could not write beside the frames ({e})"}), 400
    return jsonify(_site_reply(folder))


def _rebuild_layers(wd, prog, denoise, fnrgf_preset, partialconv, earthshine=False,
                    force=False, why=""):
    """Rebuild every detail layer from the cached merge and record the recipe
    in report.json. Mirrors the tail of the pipeline: build_layers is the last
    stage there, so nothing after it is skipped. `force`: rebuild even when the
    recipe is current, because Denoise changed -- it acts on the layers only,
    so it never needs a re-stack (2026-10-08)."""
    from . import detail
    stale, have, want = _layers_stale(wd)
    # THE OBSERVING SITE AND TIME (observer.py) decide which stars the
    # catalogue confirms, and the stars are lifted out of the corona before the
    # layers are built -- so a changed site means the layers are rebuilt
    if not stale and not force and _site_changed(wd):
        force, why = True, "observing site/time changed, stars re-checked"
    if not stale and not force:
        # only the inner-detail / filament-detail recipe moved: those two
        # need nothing but the cached corona, so they are rebuilt alone
        # (detail.FAST_BUILDS, audit 2026-10-10)
        _fast = detail.fast_stale(wd)
        if not _fast:
            return False
        prog.log(f"detail layers: {', '.join(_fast)} rebuilt from the cached corona "
                 f"(the other layers are current)", 0.5)
        lst = detail.rebuild_fast(wd, prog, _fast)
        _rp = os.path.join(wd, "report.json")
        try:
            _r = json.load(open(_rp)) if os.path.exists(_rp) else {}
            _lr = _r.setdefault("layer_recipe", {})
            _lf = _lr.setdefault("fast", {})
            _lf.update(lst.pop("layer_recipe_fast", {}))
            _r.update(lst)
            from .report import write as _write_report
            _write_report(wd, _r)
        except Exception as _e:
            prog.log(f"[warn] layers rebuilt, not recorded in report ({_e})", None)
        return True
    if stale:
        prog.log(f"detail layers: recipe {have} cached, "
                 f"{want} needed: rebuilding from the cached merge (no re-stack)", 0.5)
    else:
        prog.log(f"detail layers: {why}: rebuilding from the cached merge (no re-stack)", 0.5)
    lst = detail.build_layers(wd, prog, denoise=denoise, earthshine=earthshine,
                              fnrgf_preset=fnrgf_preset, partialconv=partialconv)
    _rp = os.path.join(wd, "report.json")
    try:
        _r = json.load(open(_rp)) if os.path.exists(_rp) else {}
        if isinstance(lst, dict):
            _r.update(lst)
        from .report import write as _write_report
        _write_report(wd, _r)          # report.json + report.txt, as the pipeline does
    except Exception as _e:
        prog.log(f"[warn] detail layers rebuilt, not recorded in report ({_e})", None)
    return True


@app.post("/api/run")
def start_run():
    # JSON ONLY (0.22.78). Both of these accepted a bodiless POST and ran with
    # defaults, so any page open in the same browser could start a run or wipe a
    # cache on 127.0.0.1 with a plain form post. Requiring a JSON body is the
    # cheap half of the fix: a cross-site form cannot set Content-Type.
    if not request.is_json:
        return jsonify({"ok": False, "error": "expected a JSON body"}), 415
    # An import needs no raw folder: its products belong beside the image. The
    # folder was required unconditionally, so pasting a TIFF path and nothing
    # else was a dead end -- reported from the field as "the start button is
    # greyed out". Fall back to the image's own directory.
    # AN IMPORT ALWAYS WORKS BESIDE ITS OWN IMAGE (0.22.78). This used to fall
    # back to the image's directory only when NO folder was loaded -- so with a
    # raw folder open, pasting a finished HDR from somewhere else wrote
    # hdr_lum.npy, geometry.json, report.json and opts.json(mode=import) into
    # that folder's .eclipseforgehdr and destroyed the stack that was already
    # there. Comparing one HDR against a bracket you had just spent a quarter of
    # an hour stacking is exactly when someone does this.
    _ip = (request.json.get("importPath", "") if request.is_json else "") or ""
    _ip = os.path.expanduser(str(_ip).strip())
    _b = _busy()
    if _b:
        # before the folder switch below: a running job keeps its own folder,
        # but an export would otherwise land in the NEW folder (audit 2026-10-10)
        return jsonify({"ok": False, "error": f"{_b} is still running"}), 400
    if _ip and os.path.isfile(_ip):
        STATE["folder"] = os.path.dirname(os.path.abspath(_ip))
    elif STATE["folder"] is None:
        return jsonify({"ok": False,
                        "error": "choose a folder first, or paste the path "
                                 "of one finished HDR to import"}), 400
    force = bool(request.json.get("force", False)) if request.is_json else False
    # DEFAULTS (2026-10-07, tested on two test sets):
    # Exact edge, Network solve, Published FNRGF, Clipped mean. The fallbacks
    # further down that read an OLD cache's opts.json keep the old values --
    # they say what that cache was built with, not what is default now.
    feather = request.json.get("feather", "taper") if request.is_json else "taper"
    if feather not in ("plain", "taper", "masked"):
        feather = "taper"
    denoise = request.json.get("denoise", "off") if request.is_json else "off"
    if denoise is True:
        denoise = "fine"
    if denoise is False:
        denoise = "off"
    # THE OLD EARTHSHINE OPTION IS GONE (2026-10-07): the Moon layer (moonlayer.py)
    # is built after every Start and is what the render uses; the old layer was
    # built and never shown, and as a stacking option it forced a re-stack when
    # toggled. Always off; old caches keep their earth.npy as the fallback.
    earthshine = False
    despeckle = bool(request.json.get("despeckle", True)) if request.is_json else True
    export_tiers = bool(request.json.get("exportTiers", False)) if request.is_json else False
    tier_linear = bool(request.json.get("tierLinear", False)) if request.is_json else False
    frames = (request.json.get("frames", "all") if request.is_json else "all")
    if frames not in ("all", "best50", "best"):
        frames = "all"
    wb_source = (request.json.get("wb", "camera") if request.is_json else "camera")
    if wb_source not in ("camera", "daylight", "none"):
        wb_source = "camera"
    photometry = (request.json.get("photometry", "linfit")
                  if request.is_json else "linfit")
    if photometry not in ("linfit", "scalar"):
        photometry = "linfit"
    photo_solve = (request.json.get("photoSolve", "network")
                   if request.is_json else "network")
    if photo_solve not in ("chain", "network"):
        photo_solve = "network"
    tier_mode = (request.json.get("tierMode", "exposure")
                 if request.is_json else "exposure")
    if tier_mode not in ("exposure", "frame", "measured"):
        tier_mode = "exposure"
    # FRAMES ALREADY REGISTERED (2026-10-13): "auto" decides from the data (a
    # writer's registration borders), "align" always aligns, "keep" never does
    align_mode = (request.json.get("alignMode", "auto")
                  if request.is_json else "auto")
    if align_mode not in ("auto", "align", "keep"):
        align_mode = "auto"
    # WHAT THE CACHE IS COMPARED ON is the decision, not the setting: a stack
    # made before this option existed was aligned, and "auto" on a set nobody
    # registered aligns too -- the same stack. A cached stack is reused only
    # when it was kept/aligned the way this run would do it.
    def _want_prealigned(folder):
        if align_mode == "keep":
            return True
        if align_mode == "align":
            return False
        try:
            from .raw import list_raws as _lr
            from .simple import writer_registered as _wr
            return bool(_wr(_lr(folder))[0])
        except Exception:
            return False
    # FNRGF: our hard order cutoff, or Druckmullerova's published attenuation.
    # Thesis p.80 gives the second as the setting that "always works without
    # producing any artifacts"; ours has never been compared against it.
    fnrgf_preset = (request.json.get("fnrgfPreset", "published")
                    if request.is_json else "published")
    if fnrgf_preset not in ("ours", "published"):
        fnrgf_preset = "published"
    # Per-tier frame combine. 'mean' is what every version up to 0.23.1 did and
    # has no per-pixel rejection at all.
    stack_combine = (request.json.get("stackCombine", "clip")
                     if request.is_json else "clip")
    if stack_combine not in ("mean", "clip"):
        stack_combine = "clip"
    # Repeated brackets: stack each exposure across the brackets, then merge
    # (the default), or merge each bracket on its own and stack the HDRs
    # (bracketmode.py, lab 2026-10-16).
    bracket_mode = (request.json.get("bracketMode", "tiers")
                    if request.is_json else "tiers")
    if bracket_mode not in ("tiers", "hdr"):
        bracket_mode = "tiers"
    # The alignment pre-filter: an isotropic high-pass (ours) or an arc-length
    # one (Druckmullerova's T_sigma, which is what removes the lunar edge and
    # the moving saturation edge).
    align_filter = (request.json.get("alignFilter", "isotropic")
                    if request.is_json else "isotropic")
    if align_filter not in ("isotropic", "tangential"):
        align_filter = "isotropic"
    # Correlation flavour: plain cross-correlation (ours), full whitening, or
    # the semi-phase correlation in between.
    align_corr = (request.json.get("alignCorr", "semi")
                  if request.is_json else "semi")
    if align_corr not in ("cross", "semi", "phase"):
        align_corr = "semi"
    # What the frames WITHIN one tier are aligned on. Separate from the two
    # above, which govern the tier-to-tier network.
    intra_lock = (request.json.get("intraLock", "corona")
                  if request.is_json else "corona")
    if intra_lock not in ("corona", "moon", "mixed"):
        intra_lock = "corona"
    # COLOUR PLANES: move red and blue onto green by the offset measured on the
    # lunar limb (atmospheric dispersion at a low Sun). Changes the merged
    # data, so it is in the cache key below.
    # Default "auto" since 0.23.9: it measures on every run and
    # moves nothing when the offset is below 0.2 px.
    colour_planes = (request.json.get("colourPlanes", "auto")
                     if request.is_json else "auto")
    if colour_planes not in ("off", "auto"):
        colour_planes = "auto"
    # PARTIAL CONVOLUTION ON/OFF. About half the run time on a big merge and
    # nothing else depends on it, so it is turned off while the rest of the
    # settings are being found and on for the final render. NOT part of the
    # stack cache key below: it changes only which detail layers are built,
    # not the stack, so switching it does not force a re-stack -- turning it
    # back on costs the mask build alone, the same path a folder stacked
    # before 0.22.64 takes.
    partialconv = (bool(request.json.get("partialConv", True))
                   if request.is_json else True)
    # VNG IS NO LONGER A TOOLBAR CHOICE. It was added in 0.22.63 so a render
    # could be compared with PixInsight's on equal terms, and that comparison
    # has now been made on the reference bracket: no visible advantage, and the
    # colour runs into the limb slightly less smoothly -- which is exactly what
    # docs/PIXINSIGHT_COMPARISON.md predicted, since VNG is gradient-directed
    # and the steepest gradient in the frame is the limb.
    #
    # The implementation stays (raw.demosaic_vng, verified against LibRaw to
    # 1.5 counts in 65535) and is selected by dropping a file named
    # `vng_demosaic.txt` in the folder of raws, the same way the linear-fit
    # photometry chain is. Off the toolbar, out of the way, still there for the
    # next person who wants to check us against PixInsight.
    # STATE["folder"] rather than `folder`, which is not bound until below
    demosaic_method = "vng" if os.path.exists(os.path.join(
        STATE["folder"] or "", "vng_demosaic.txt")) else "mhc"
    flat_dir = (request.json.get("flatDir", "") if request.is_json else "") or ""
    flat_dir = str(flat_dir).strip()
    # One finished HDR instead of a bracket: skip stacking entirely and run the
    # detail layers on the supplied image.
    import_path = (request.json.get("importPath", "") if request.is_json else "") or ""
    import_path = os.path.expanduser(str(import_path).strip())
    # THE FALLBACK. Align, average, merge on ratios measured from the pixels,
    # and nothing else -- then the same layers every other run gets. It exists
    # so that a bracket the normal path cannot handle still produces a clean
    # stack instead of nothing. See simple.py for what it leaves out.
    simple = bool(request.json.get("simple", False)) if request.is_json else False
    if import_path and not os.path.isfile(import_path):
        return jsonify({"ok": False,
                        "error": f"not a file: {import_path}"}), 400
    prog = Progress()
    STATE["progress"] = prog
    # Bind the folder HERE. work() used to re-read STATE["folder"] at execution
    # time, so picking another folder while a run was in flight processed the
    # new one but loaded the old one's layers.
    folder = STATE["folder"]
    # and drop the previous result: it used to stay live and exportable for the
    # whole re-run, so a run that failed still reported ready and exported the
    # OLD layers under the new settings.
    STATE["layers"] = None

    def work():
        try:
            wd = workdir(folder)
            opts_path = os.path.join(wd, "opts.json")
            from . import __version__, cache_ok as _cache_ok
            from .pipeline import resolve_flat_dir
            from .flat import fingerprint as _flat_fp
            if simple and not import_path:
                from . import simple as _simple
                o = {}
                if os.path.exists(opts_path):
                    try:
                        o = json.load(open(opts_path))
                    except Exception:
                        o = {}
                _fd = resolve_flat_dir(folder, flat_dir)
                same = (o.get("mode") == "simple"
                        and o.get("denoise") == denoise
                        and o.get("fnrgf_preset") == fnrgf_preset
                        # calibration frames change the stacked data, and
                        # nothing the page sends would announce a new flats
                        # folder: compare what was actually applied
                        and o.get("flat_dir", "") == (_fd or "")
                        and o.get("flat_inputs") == _flat_fp(_fd)
                        # ... and whether the partial-convolution layer was
                        # built. Without this, unticking it after a cached
                        # simple run reused the layers and the masks stayed.
                        and bool(o.get("partialconv", True)) == bool(partialconv)
                        # tiers by header or by measured brightness group
                        # the same frames differently
                        and o.get("tier_mode", "exposure") == tier_mode
                        and bool(o.get("prealigned", False)) == _want_prealigned(folder)
                        # a stack made by an older simple stack is redone
                        and int(o.get("simple_rev", 1) or 1)
                        >= getattr(_simple, "SIMPLE_REV", 1)
                        and _cache_ok(o.get("build")))
                have = all(os.path.exists(os.path.join(wd, f))
                           for f in ("prom.npy", "prom_rgb.npy", "pellett.npy"))
                if force or not same or not have:
                    if os.path.exists(opts_path):
                        try:
                            os.remove(opts_path)
                        except OSError:
                            pass
                    _simple.run(folder, prog, denoise=denoise,
                                demosaic_method=demosaic_method,
                                fnrgf_preset=fnrgf_preset,
                                flat_dir=flat_dir,
                                partialconv=partialconv,
                                tier_mode=tier_mode,
                                align_mode=align_mode)
                    try:
                        _o2 = json.load(open(opts_path))
                        _o2["flat_inputs"] = _flat_fp(_fd)
                        json.dump(_o2, open(opts_path, "w"), indent=1)
                    except Exception:
                        pass
                else:
                    prog.log("using cached simple stack", 0.9)
                    _rebuild_layers(wd, prog, denoise, fnrgf_preset, partialconv)
                # THE MOON LAYER ON THE SIMPLE STACK TOO (2026-10-08, a FITS
                # test set): it reads the long frames itself and fits its own
                # scale against the merge, so it needs nothing the full
                # pipeline has and the simple stack lacks. The simple stack
                # sits in the first frame's coordinates (crop origin 0, 0).
                try:
                    from . import moonlayer as _ml
                    if _ml.needs_build(wd):
                        prog.log("building Moon layer (earthshine) from long frames...", None)
                        _ml.build(folder, wd, prog)
                except Exception as _e:
                    if getattr(prog, "cancelled", False):
                        raise
                    prog.log(f"[warn] Moon layer not built ({_e})", None)
                _build_stars(wd, prog)
                prog.log("loading layers for preview...", None)
                STATE["layers"] = Layers(wd)
                prog.log("ready", 1.0)
                prog.done = True
                return
            if import_path:
                from . import importhdr
                o = {}
                if os.path.exists(opts_path):
                    try:
                        o = json.load(open(opts_path))
                    except Exception:
                        o = {}
                same = (o.get("mode") == "import"
                        and o.get("import") == os.path.basename(import_path)
                        and o.get("import_mtime") == int(os.path.getmtime(import_path))
                        and o.get("import_size") == int(os.path.getsize(import_path))
                        and o.get("denoise") == denoise
                        and bool(o.get("partialconv", True)) == bool(partialconv)
                        # NOT the feather (0.22.78): importhdr never writes that
                        # key, so comparing it against the toolbar's Merge weight
                        # made the test fail for every value but 'plain' and
                        # re-imported the same image on every Start. The feather
                        # is a MERGE setting; an import has no merge.
                        and _cache_ok(o.get("build")))
                have = all(os.path.exists(os.path.join(wd, f))
                           for f in ("prom.npy", "prom_rgb.npy", "pellett.npy"))
                if force or not same or not have:
                    if os.path.exists(opts_path):
                        try:
                            os.remove(opts_path)
                        except OSError:
                            pass
                    importhdr.run(folder, import_path, prog, denoise=denoise,
                                  partialconv=partialconv, fnrgf_preset=fnrgf_preset)
                else:
                    prog.log("using cached import", 0.9)
                    _rebuild_layers(wd, prog, denoise, fnrgf_preset, partialconv)
                prog.log("loading layers for preview...", None)
                STATE["layers"] = Layers(wd)
                prog.log("ready", 1.0)
                prog.done = True
                return
            _fd = resolve_flat_dir(folder, flat_dir)
            from .pipeline import resolve_calib_dir as _rcd
            from .dark import fingerprint as _calib_fp
            _bd = _rcd(folder, "bias")
            _dd = _rcd(folder, "dark")
            opts_ok = False
            o = None
            if os.path.exists(opts_path):
                try:
                    o = json.load(open(opts_path))
                except Exception as _e:
                    # a kill mid-write leaves a truncated opts.json; treat it
                    # as stale instead of failing every Start until Clear cache
                    prog.log(f"[warn] opts.json unreadable ({_e}): treating the cache as stale", None)
                    o = None
            if o is not None:
                # Denoise is a LAYER option (2026-10-08): a change rebuilds the
                # detail layers below, never the stack
                opts_ok = (o.get("flat_dir") == _fd
                           and o.get("flat_inputs") == _flat_fp(_fd)
                           # Dropping a bias/ or darks/ folder beside the raws
                           # changes the stacked data, and nothing the page
                           # sends would announce it -- exactly the trap the
                           # flats key exists to close.
                           and o.get("bias_dir") == _bd
                           and o.get("dark_dir") == _dd
                           and o.get("bias_inputs", []) == _calib_fp(_bd)
                           and o.get("dark_inputs", []) == _calib_fp(_dd)
                           and bool(o.get("despeckle", True)) == despeckle
                           and o.get("frames", "all") == frames
                           and bool(o.get("export_tiers", False)) == export_tiers
                           and bool(o.get("tier_linear", False)) == tier_linear
                           # The merge weight changes the merged data itself, so
                           # it belongs in the key. It was checked on the import
                           # path only, which meant switching Detail/Photometric
                           # on a raw folder quietly reused the previous merge
                           # and the setting appeared to do nothing.
                           and o.get("feather", "plain") == feather
                           # Both of these change the MERGED data too: the
                           # white balance scales the raw channels before the
                           # camera matrix, and the demosaic builds the pixels.
                           # The default here is the PRE-0.22.63 behaviour, so
                           # an old workdir compares as stale and restacks
                           # rather than being reused under the new default.
                           and o.get("wb_source", "daylight") == wb_source
                           and o.get("demosaic", "mhc") == demosaic_method
                           # ... and the photometric mode, for the same reason:
                           # it changes the merged data, and it is selected by
                           # a marker file in the raw folder rather than by
                           # anything the page sends, so nothing else here
                           # would notice it appearing or disappearing.
                           # NOTE the stored value is still "hill"/"scalar",
                           # unchanged by 0.22.86 making the fit the default: a
                           # folder already stacked WITH it compares equal and
                           # keeps its stack, and one stacked without it
                           # correctly re-stacks under the new default.
                           and (o.get("photometry", "scalar")
                                == ("hill" if _photometry_requested(
                                        folder, photometry) else "scalar"))
                           # ... and chain vs network, which changes cal[] and
                           # therefore every merged tier. Defaulting the stored
                           # value to "chain" is what lets a folder stacked
                           # before this setting existed keep its cache.
                           and o.get("photo_solve", "chain") == photo_solve
                           # fnrgf_preset is a LAYER option: handled like
                           # denoise below (audit 2026-10-10; it forced a
                           # full re-stack before)
                           and o.get("stack_combine", "mean") == stack_combine
                           and o.get("bracket_mode", "tiers") == bracket_mode
                           and o.get("align_filter", "isotropic") == align_filter
                           and o.get("align_corr", "semi") == align_corr
                           # ... and what the frames within a tier locked on,
                           # which changes every tier's stack. Defaulting the
                           # stored value to "moon" keeps the cache of a folder
                           # stacked before 0.23.8, when that was the only
                           # behaviour.
                           and o.get("intra_lock", "moon") == intra_lock
                           and o.get("colour_planes", "off") == colour_planes
                           # ... and the grouping, which changes what a tier IS
                           and o.get("tier_mode", "exposure") == tier_mode
                           # ... and whether the tiers were aligned or kept
                           and bool(o.get("prealigned", False)) == _want_prealigned(folder)
                           # ... and sub-pixel intra-tier stacking (lab
                           # 0.25.1): a stack made before it compares as stale
                           and bool(o.get("intra_subpixel", False)) == _intra_sub()
                           # ... and the bisect switches the stack was made with
                           and (o.get("switches") or {}) == _env_switches()
                           and _cache_ok(o.get("build"))
                           # ... and the input files themselves. Without this,
                           # adding or removing a frame and pressing Start (not
                           # force) silently reused the previous stack, and the
                           # exported report quoted the previous frame count.
                           and o.get("inputs") == _input_fingerprint(folder)
                           # the short-exposure subfolder (prominence layer);
                           # a cache from before this key is not re-stacked for it
                           and ("short_inputs" not in o
                                or o.get("short_inputs") == _short_fingerprint(folder)))
            have_all = all(os.path.exists(os.path.join(wd, f))
                           for f in ("prom.npy", "prom_rgb.npy", "pellett.npy"))
            if not force and have_all and not opts_ok:
                # SAY WHY. A re-stack that nobody asked for costs half an hour
                # and, without this line, an afternoon of guessing (0.23.9 had
                # three keys missing from opts.json; the next one will be
                # something else). Every key the test above compares, with the
                # cached and the requested value.
                try:
                    _want = {"flat_dir": _fd, "bias_dir": _bd, "dark_dir": _dd,
                             "earthshine": earthshine, "despeckle": despeckle, "frames": frames,
                             "export_tiers": export_tiers, "tier_linear": tier_linear,
                             "feather": feather, "wb_source": wb_source, "demosaic": demosaic_method,
                             "photometry": ("hill" if _photometry_requested(folder, photometry) else "scalar"),
                             "photo_solve": photo_solve, "fnrgf_preset": fnrgf_preset,
                             "stack_combine": stack_combine, "align_filter": align_filter,
                             "bracket_mode": bracket_mode,
                             "align_corr": align_corr, "intra_lock": intra_lock,
                             "colour_planes": colour_planes, "tier_mode": tier_mode,
                             "prealigned": _want_prealigned(folder),
                             "intra_subpixel": _intra_sub(), "switches": _env_switches()}
                    _dflt = {"wb_source": "daylight", "demosaic": "mhc", "photometry": "scalar",
                             "photo_solve": "chain", "fnrgf_preset": "ours", "stack_combine": "mean",
                             "bracket_mode": "tiers",
                             "align_filter": "isotropic", "align_corr": "semi", "intra_lock": "moon",
                             "colour_planes": "off", "tier_mode": "exposure", "feather": "plain",
                             "prealigned": False,
                             "frames": "all", "intra_subpixel": False, "switches": {}}
                    _diff = []
                    _want.pop("fnrgf_preset", None)
                    for _k, _v in _want.items():
                        _c = o.get(_k, _dflt.get(_k, False if isinstance(_v, bool) else None))
                        if isinstance(_v, bool):
                            _c = bool(_c)
                        if _c != _v:
                            _diff.append(f"{_k}: cached {_c!r}, requested {_v!r}")
                    if o.get("flat_inputs") != _flat_fp(_fd):
                        _diff.append("flats folder contents changed")
                    if o.get("bias_inputs", []) != _calib_fp(_bd):
                        _diff.append("bias folder contents changed")
                    if o.get("dark_inputs", []) != _calib_fp(_dd):
                        _diff.append("darks folder contents changed")
                    if not _cache_ok(o.get("build")):
                        from . import __version__ as _ver
                        _diff.append(f"build: cached {o.get('build')!r}, this is {_ver!r}")
                    if o.get("inputs") != _input_fingerprint(folder):
                        _diff.append("frames changed (files added, removed or replaced, or switched in the Frames panel)")
                    if "short_inputs" in o and o.get("short_inputs") != _short_fingerprint(folder):
                        _diff.append("short-exposure frames changed")
                    prog.log("re-stacking, cache mismatch: "
                             + ("; ".join(_diff) if _diff else "no key differs (please report)"), None)
                except Exception as _e:
                    prog.log(f"[odd] re-stacking (cache check failed: {_e})", None)
            if force or not have_all or not opts_ok:
                # A run that dies partway leaves a valid-looking opts.json from
                # the previous run beside a mix of new and old products; clear it
                # first so a later non-forced Start cannot declare that mixture
                # valid and load layers from two different builds.
                if os.path.exists(opts_path):
                    try:
                        os.remove(opts_path)
                    except OSError:
                        pass
                run_pipeline(folder, prog, denoise=denoise, feather=feather,
                             earthshine=earthshine, despeckle=despeckle,
                             frames=frames, export_tiers=export_tiers,
                             tier_linear=tier_linear, flat_dir=flat_dir,
                             wb_source=wb_source,
                             demosaic_method=demosaic_method,
                             photometry=photometry,
                             photo_solve=photo_solve,
                             fnrgf_preset=fnrgf_preset,
                             stack_combine=stack_combine,
                             bracket_mode=bracket_mode,
                             align_filter=align_filter,
                             align_corr=align_corr,
                             intra_lock=intra_lock,
                             colour_planes=colour_planes,
                             partialconv=partialconv,
                             tier_mode=tier_mode,
                             align_mode=align_mode)
            else:
                prog.log("using cached stack", 0.9)
                # 0.23.9: a newer detail-layer recipe is applied to the cached
                # merge here. build_layers also builds the partial-convolution
                # masks when they are on, so the block below then finds them
                # fresh and does nothing.
                try:
                    _oc = json.load(open(opts_path)) if os.path.exists(opts_path) else {}
                except Exception:
                    _oc = {}
                _dn_changed = bool(_oc) and _oc.get("denoise") != denoise
                _fp_changed = bool(_oc) and _oc.get("fnrgf_preset", "ours") != fnrgf_preset
                _why = (f"denoise {_oc.get('denoise')!r} -> {denoise!r}" if _dn_changed
                        else f"FNRGF preset {_oc.get('fnrgf_preset', 'ours')!r} -> {fnrgf_preset!r}")
                _rebuild_layers(wd, prog, denoise, fnrgf_preset, partialconv,
                                earthshine=earthshine, force=_dn_changed or _fp_changed,
                                why=_why)
                if _dn_changed or _fp_changed:
                    try:
                        _oc["denoise"] = denoise
                        _oc["fnrgf_preset"] = fnrgf_preset
                        json.dump(_oc, open(opts_path + ".tmp", "w"), indent=1)
                        os.replace(opts_path + ".tmp", opts_path)
                        # the report prints the options too (audit 2026-10-10)
                        try:
                            _rp = os.path.join(wd, "report.json")
                            _rj = json.load(open(_rp))
                            _rj.setdefault("options", {})["denoise"] = denoise
                            _rj["options"]["fnrgf_preset"] = fnrgf_preset
                            json.dump(_rj, open(_rp + ".tmp", "w"), indent=1)
                            os.replace(_rp + ".tmp", _rp)
                        except Exception:
                            pass
                    except Exception:
                        pass
            # PARTIAL-CONVOLUTION MASKS FROM A CACHED STACK. Everything they need -- the merged
            # luminance, the geometry, the prominence mask -- is already on
            # disk, so a folder stacked before 0.22.64 gains the layer for the
            # cost of the polar blur alone instead of a whole re-stack.
            # SWITCHED OFF: SKIP THIS WHOLE BLOCK, not just the staleness
            # test. Setting _hstale = False here was not enough -- the next
            # branch is `if not _hstale`, which re-derives staleness from
            # report.json, finds no "hill" entry because the run did not write
            # one, reads its build number as 1, and rebuilds. That is exactly
            # what happened on a run of the reference set: the pipeline reported
            # "partial convolution: OFF" and the masks were then built anyway
            # right after "pipeline complete", turning an 18m39s run into 39
            # minutes. Whatever is on disk is also dropped, since it was built
            # on a different merge and nothing downstream could tell.
            if not partialconv:
                for _hf in ("hill.npy", "hill_log.npy", "hill_sigma.npy"):
                    _hfp = os.path.join(wd, _hf)
                    if os.path.exists(_hfp):
                        try:
                            os.remove(_hfp)
                        except OSError:
                            pass
                prog.log("partial convolution off: masks not built, "
                         "cached masks dropped", None)
                _hstale = False
            else:
                _hstale = not os.path.exists(os.path.join(wd, "hill.npy"))
            if partialconv and not _hstale:
                # masks built by an older recipe are rebuilt, not reused: the
                # de-radialisation in 0.22.67 changes what is ON DISK
                from .detail import HILL_BUILD as _HB
                try:
                    _rj = json.load(open(os.path.join(wd, "report.json")))
                    _hj = (_rj.get("hill") or {})
                    _hstale = int(_hj.get("build", 1)) < _HB
                    # ... and rebuilt when the BASE changed, which the build
                    # number cannot see: ECLIPSEFORGE_HILL_NODENOISE picks
                    # between the raw and the denoised master, and without this
                    # the second run of an A/B silently reuses the first run's
                    # masks and reports itself as a test.
                    _want = ("denoised" if (
                        os.environ.get("ECLIPSEFORGE_HILL_DENOISED") == "1"
                        and os.environ.get("ECLIPSEFORGE_HILL_NODENOISE") != "1")
                        else "raw")
                    if not _hstale and _hj.get("base", "denoised") != _want:
                        _hstale = True
                        prog.log("partial-convolution masks: "
                                 "%s master cached, %s needed: "
                                 "rebuilding"
                                 % (_hj.get("base", "denoised"), _want), None)
                    # ... and when the radial smoothing factor changed
                    from .detail import HILL_RADIAL as _HR
                    _wr = float(os.environ.get("ECLIPSEFORGE_HILL_RADIAL", _HR))
                    if not _hstale and abs(float(_hj.get("radial", 0.0)) - _wr) > 1e-6:
                        _hstale = True
                        prog.log("partial-convolution masks: "
                                 "radial smoothing %g cached, "
                                 "%g needed: rebuilding"
                                 % (float(_hj.get("radial", 0.0)), _wr), None)
                    from .detail import HILL_GATE as _HG, HILL_GATE_LEN as _HGL
                    _wg = (float(os.environ.get("ECLIPSEFORGE_HILL_GATE", _HG)),
                           float(os.environ.get("ECLIPSEFORGE_HILL_GATE_LEN", _HGL)))
                    _hg = (float(_hj.get("gate", 0.0)), float(_hj.get("gate_len", 0.0)))
                    if not _hstale and _wr > 0 and (abs(_hg[0] - _wg[0]) > 1e-6 or (
                            _wg[0] > 0 and abs(_hg[1] - _wg[1]) > 1e-6)):
                        _hstale = True
                        prog.log("partial-convolution masks: "
                                 "gate %g/%g cached, "
                                 "%g/%g needed: rebuilding"
                                 % (_hg[0], _hg[1], _wg[0], _wg[1]), None)
                except Exception:
                    _hstale = True
                if _hstale:
                    prog.log("partial-convolution masks: older recipe cached: "
                             "rebuilding", None)
            if partialconv and _hstale:
                from .detail import build_hill
                _hst = build_hill(wd, prog, denoise=denoise)
                if _hst:
                    _rp = os.path.join(wd, "report.json")
                    try:
                        _r = json.load(open(_rp)) if os.path.exists(_rp) else {}
                        _r["hill"] = _hst
                        json.dump(_r, open(_rp, "w"))
                    except Exception as _e:
                        prog.log(f"[warn] partial-convolution masks built, not recorded in "
                                 f"report ({_e})", None)
            # THE MOON LAYER (lab 0.24): built on its own, about a minute, when
            # missing or older than the merge -- no re-stack, no layer rebuild
            try:
                from . import moonlayer as _ml
                if _ml.needs_build(wd):
                    prog.log("building Moon layer (earthshine) from long frames...", None)
                    _ml.build(folder, wd, prog)
            except Exception as _e:
                if getattr(prog, "cancelled", False):
                    raise
                prog.log(f"[warn] Moon layer not built ({_e})", None)
            _build_stars(wd, prog)
            prog.log("loading layers for preview...", None)
            STATE["layers"] = Layers(wd)
            prog.log("ready", 1.0)
            prog.done = True
        except Exception as e:
            if getattr(prog, "cancelled", False):
                prog.error = ("cancelled by you -- nothing from this run is used; "
                              "the next Start runs it again")
            else:
                prog.error = f"{e}\n{traceback.format_exc()}"
            prog.done = True

    t = threading.Thread(target=work, daemon=True, name="the pipeline run")
    STATE["thread"] = t
    t.start()
    return jsonify({"ok": True})


@app.post("/api/cancel")
def cancel_job():
    """Stop the running job (stack, layer build, export) at its next step.
    A cancelled stack leaves no valid cache: opts.json is removed before a
    run starts and written only at its end, so the next Start redoes it."""
    p = STATE.get("progress")
    if p is None or not _busy() or getattr(p, "done", True):
        return jsonify({"ok": False, "error": "nothing is running"})
    p.cancel()
    return jsonify({"ok": True})


@app.get("/api/progress")
def get_progress():
    p = STATE["progress"]
    if p is None:
        return jsonify({"lines": [], "levels": [], "frac": 0, "done": False, "error": None,
                        "elapsed": 0.0, "since": 0.0, "step": "",
                        "ready": STATE["layers"] is not None})
    # "step" is the last log line: it is what the run is actually doing, and
    # pairing it with "since" tells the user whether a long wait is a slow step
    # or a dead one.
    _lv = getattr(p, "levels", None) or []
    return jsonify({"lines": p.lines[-250:], "frac": p.frac, "done": p.done,
                    "levels": _lv[-250:] if len(_lv) == len(p.lines) else [],
                    "error": p.error,
                    "cancelled": bool(getattr(p, "cancelled", False)),
                    "cancelling": bool(getattr(p, "cancel_requested", False)) and not p.done,
                    "elapsed": round(p.elapsed(), 1),
                    "since": round(p.since(), 1),
                    "step": p.lines[-1] if p.lines else "",
                    "ready": STATE["layers"] is not None})


@app.get("/api/cachestate")
def cache_state():
    """What the next Start would do with the loaded folder, for the GUI's Start
    label and the amber marks on the options: the cached run's options (minus
    the input lists), whether the stack is reusable at all, and whether the
    layer set or the partial-convolution masks are behind the build."""
    folder = STATE.get("folder")
    if not folder:
        return jsonify({"ok": False, "error": "no folder"}), 400
    wd = workdir(folder)
    out = {"ok": True, "cached": False, "opts": None, "layers_stale": False,
           "hill_stale": False, "mode": None}
    _op = os.path.join(wd, "opts.json")
    if os.path.exists(os.path.join(wd, "geometry.json")) and os.path.exists(_op):
        try:
            o = json.load(open(_op))
            from . import cache_ok
            out["cached"] = bool(cache_ok(o.get("build")))
            out["mode"] = o.get("mode")
            out["opts"] = {k: v for k, v in o.items()
                           if k not in ("inputs", "flat_inputs", "bias_inputs", "dark_inputs")}
            out["inputs_changed"] = (o.get("inputs") != _input_fingerprint(folder))
            # the other frame folders, which the Frames panel can switch too
            try:
                from .flat import fingerprint as _ffp
                from .dark import fingerprint as _dfp
                out["inputs_changed"] = bool(
                    out["inputs_changed"]
                    or ("short_inputs" in o and o.get("short_inputs") != _short_fingerprint(folder))
                    or ("flat_inputs" in o and o.get("flat_applied") and o.get("flat_inputs") != _ffp(o.get("flat_dir")))
                    or ("bias_inputs" in o and o.get("bias_applied") and o.get("bias_inputs") != _dfp(o.get("bias_dir")))
                    or ("dark_inputs" in o and o.get("dark_applied") and o.get("dark_inputs") != _dfp(o.get("dark_dir"))))
            except Exception:
                pass
        except Exception as e:
            out["error"] = str(e)
    if out["cached"]:
        try:
            from .detail import fast_stale as _fs
            out["layers_stale"] = bool(_layers_stale(wd)[0] or _fs(wd))
        except Exception:
            pass
        try:
            from .detail import HILL_BUILD as _HB
            if os.path.exists(os.path.join(wd, "hill.npy")):
                _rj = json.load(open(os.path.join(wd, "report.json")))
                out["hill_stale"] = int(((_rj.get("hill") or {}).get("build", 1))) < _HB
            else:
                out["hill_stale"] = True
        except Exception:
            out["hill_stale"] = True
    return jsonify(out)


@app.get("/api/geometry")
def get_geometry():
    ly = STATE["layers"]
    if ly is None:
        return jsonify({"ok": False}), 400
    from . import __version__
    cy, cx, R = ly.geometry(ly.prev_decim)
    h, w = ly.prev["bg"].shape
    try:                       # the export's size, for the Export page
        fh, fw = ly.full["bg"].shape[:2]
    except Exception:
        fh, fw = h * ly.prev_decim, w * ly.prev_decim
    return jsonify({"ok": True, "W": w, "H": h, "cy": cy, "cx": cx, "R": R,
                    "fullW": int(fw), "fullH": int(fh),
                    "defaults": defaults_for(getattr(ly, "mode", None)),
                    "version": __version__,
                    "app": _BRAND_NAME, "lab": _BRAND_LAB,
                    "Rmask": ly.mask_radius(ly.prev_decim),
                    "limbProf": ([float(x) / ly.prev_decim for x in ly.limb_prof]
                                 if ly.limb_prof else None),
                    "limbMargin": ly.limb_margin / ly.prev_decim,
                    "bgChroma": [float(x) for x in ly.bg_chroma],
                    "rprofile": [round(float(x), 5) for x in ly.rprof_1d],
                    # the Sky pedestal (2026-10-19): the sky light measured
                    # in xn units, and the profile's scale -- see render.Layers
                    "skyP3": float(getattr(ly, "sky_p3", 0.0) or 0.0),
                    "profMax": float(getattr(ly, "prof_max", 1.0) or 1.0),
                    "coronaGain": [float(x) for x in ly.corona_gain],
                    "coronaRoverGB": ly.corona_r_over_gb,
                    "hasHill": bool(getattr(ly, "has_hill", False)),
                    "nHill": int(getattr(ly, "n_hill", 0)),
                    "hillScales": [float(x) for x in getattr(ly, "hill_scales", [])],
                    "hillRms": [float(x) for x in getattr(ly, "hill_rms", [])],
                    # the STRUCTURE rms, which is what render.py normalises the
                    # gain to since 0.22.84. The page used hillRms (total) until
                    # 0.23.7: equal on a denoised set, 4x apart on a raw one,
                    # so the export came out 4x stronger than the preview.
                    "hillRmsStruct": [float(x) for x
                                      in getattr(ly, "hill_rms_struct", [])],
                    # 0.23.9: the scale the export actually divides by (render.py,
                    # Layers.hill_ref / HILL_REF_C). The page must use the same
                    # number or the preview and the export part again.
                    "hillR0": _hill_r0(ly),
                    "hillLogKBuild": float(getattr(ly, "hill_logk", 6.0)),
                    "hillResp": [float(x) for x in getattr(ly, "hill_resp", [])],
                    # the page renders at preview scale, so it needs the same
                    # per-scale threshold correction render(preview=True) uses
                    "hillPrevAtt": [float(x) for x
                                    in getattr(ly, "hill_prev_att", [])],
                    "hasHillSigma": bool(getattr(ly, "has_hill_sigma", False)),
                    "hillSigmaScale": (
                        float(np.percentile(np.asarray(ly.prev["hillsigma"],
                                                       np.float32), 99.9))
                        if "hillsigma" in getattr(ly, "prev", {}) else 1.0),
                    "hillSpan": [_hill_span(ly, i)
                                 for i in range(int(getattr(ly, "n_hill", 0)))],
                    "decim": ly.prev_decim, "has_contact": ly.has_contact,
                    "has_flat": getattr(ly, "has_flat", False),
                    "has_promdet": getattr(ly, "has_promdet", False),
                    "hasPromLayer": bool(getattr(ly, "has_promlayer", False)),
                    "hasStars": bool(getattr(ly, "has_stars", False)),
                    "hasInnerDetail": bool(getattr(ly, "has_inner_detail", False)),
                    "hasBandDetail": bool(getattr(ly, "has_band_detail", False)),
                    "hasBandFine": bool(getattr(ly, "has_band_fine", False)),
                    "nStars": len(getattr(ly, "st_full", []) or []),
                    # per star on the preview grid: centre, gain ring, the
                    # corona's grain there, patch box and colour (stars.overlay)
                    "starCat": _star_summary(ly),
                    "stars": [{"cy": float(t[0]), "cx": float(t[1]), "rin": float(t[2]),
                               "rout": float(t[3]), "sig": float(t[4]), "y0": int(t[5]),
                               "x0": int(t[6]), "h": int(t[7].shape[0]), "w": int(t[7].shape[1]),
                               "c": [float(v) for v in np.asarray(t[9]).reshape(-1)[:3]]}
                              for t in (getattr(ly, "st_prev", []) or [])],
                    "stLMax": float(getattr(ly, "st_lmax", 1.0)),
                    "has_nafe": getattr(ly, "has_nafe", True),
                    "has_mgn_fine": getattr(ly, "has_mgn_fine", False),
                    "flat_range": list(getattr(ly, "flat_range", []) or []),
                    "flat_error": getattr(ly, "flat_error", None),
                    "has_earth": getattr(ly, "has_earth", False)})


@app.post("/api/contact")
def load_contact():
    _b = _busy()
    if _b:
        return jsonify({"ok": False, "error": f"{_b} is still running"}), 400
    ly = STATE["layers"]
    if ly is None:
        return jsonify({"ok": False, "error": "run the pipeline first"}), 400
    path = os.path.expanduser(request.json.get("path", "").strip())
    if not os.path.isfile(path):
        return jsonify({"ok": False, "error": f"not a file: {path}"}), 400
    prog = Progress()
    STATE["progress"] = prog

    def work():
        try:
            from .pipeline import prepare_contact
            prepare_contact(STATE["folder"], path, prog)
            ly.reload_contact()
            prog.done = True
        except Exception as e:
            prog.error = str(e)
            prog.done = True

    _t = threading.Thread(target=work, daemon=True, name="the contact-frame load")
    STATE["busy"] = _t
    _t.start()
    return jsonify({"ok": True})


def _hill_r0(ly):
    """The Amplification slider's divisor, exactly as render.render uses it."""
    try:
        from .render import HILL_REF_C
        if (os.environ.get("ECLIPSEFORGE_HILL_OLDNORM") != "1"
                and float(getattr(ly, "hill_ref", 0.0)) > 0.0):
            return float(ly.hill_ref) / HILL_REF_C
    except Exception:
        pass
    return 0.0          # the page falls back to hillRmsStruct / hillRms


def _hill_span(ly, i, k=8.0):
    """Half-range of the transport encoding for partial-convolution mask i, in its own units.
    8 rms clips a hard prominence edge and nothing else -- measured on the
    reference bracket at 0.02% of coronal pixels on the finest scale."""
    rms = getattr(ly, "hill_rms", None)
    v = float(rms[i]) if rms and i < len(rms) else 0.0
    return max(v, 1e-9) * k


@app.get("/api/layer/<name>")
def get_layer(name):
    ly = STATE["layers"]
    if ly is None or name not in ly.prev:
        return "no layer", 404
    arr = ly.prev[name]
    if name == "ratio":
        # same clip as render.py's np.clip(ratio, 0.2, 3.0), so the preview's
        # saturation does not diverge from the export in deep prominence red
        img = Image.fromarray((np.clip(arr, 0.2, 3.0) / 3.0 * 255).astype(np.uint8))
    elif name == "bg":
        # SENT AT 16 BITS, packed high byte / low byte into R and G. The log
        # envelope has to CUBE this layer to get back to luminance, which
        # multiplies the quantisation error by three and does it where the far
        # field lives: at 8 bits the outer corona sits near byte 16, where one
        # step is 19% in xn. At 16 bits it is 0.005%.
        q = np.clip(arr, 0, 1)
        u = np.rint(q * 65535.0).astype(np.uint16)
        rgb = np.zeros(u.shape + (3,), np.uint8)
        rgb[:, :, 0] = (u >> 8).astype(np.uint8)
        rgb[:, :, 1] = (u & 0xFF).astype(np.uint8)
        img = Image.fromarray(rgb)
    elif name == "hilllog" or name.startswith("hill"):
        # 16-BIT, PACKED, AND SIGNED FOR THE MASKS. These carry the finest
        # structure in the picture; sending them 8-bit would quantise a
        # thousandth-amplitude residual into steps as large as the residual.
        # The masks are symmetric about zero, so they go out mapped from
        # [-hillSpan*rms, +hillSpan*rms] onto 0..65535 and the page inverts it
        # with the span it is given in the geometry.
        if name == "hilllog":
            q = np.clip(arr, 0, 1)
        elif name == "hillsigma":
            # sigma of im_log; small and positive. Sent over its own 99.9th
            # percentile so the 16 bits are spent where the values are.
            sc = float(np.percentile(np.asarray(arr, np.float32), 99.9)) or 1.0
            q = np.clip(np.asarray(arr, np.float32) / sc, 0.0, 1.0)
        else:
            i = int(name[4:])
            sc = _hill_span(ly, i)
            q = np.clip(np.asarray(arr, np.float32) / sc, -1.0, 1.0) * 0.5 + 0.5
        u = np.rint(q * 65535.0).astype(np.uint16)
        rgb = np.zeros(u.shape + (3,), np.uint8)
        rgb[:, :, 0] = (u >> 8).astype(np.uint8)
        rgb[:, :, 1] = (u & 0xFF).astype(np.uint8)
        img = Image.fromarray(rgb)
    elif name == "pl_L":
        # the prominence layer's luminance over its own 99.97th percentile:
        # cores run above 1, so it goes out over 0..4, 16-bit packed like bg
        q = np.clip(np.asarray(arr, np.float32) / 4.0, 0, 1)
        u = np.rint(q * 65535.0).astype(np.uint16)
        rgb = np.zeros(u.shape + (3,), np.uint8)
        rgb[:, :, 0] = (u >> 8).astype(np.uint8)
        rgb[:, :, 1] = (u & 0xFF).astype(np.uint8)
        img = Image.fromarray(rgb)
    elif name == "pl_C":
        # its chroma at unit luminance: red runs to ~2.5, so 0..4
        img = Image.fromarray((np.clip(np.asarray(arr, np.float32) / 4.0, 0, 1) * 255 + 0.5)
                              .astype(np.uint8))
    elif name == "st_L":
        # star light in base-envelope units, over its own maximum, 16-bit packed
        sc = float(getattr(ly, "st_lmax", 1.0)) or 1.0
        q = np.clip(np.asarray(arr, np.float32) / sc, 0, 1)
        u = np.rint(q * 65535.0).astype(np.uint16)
        rgb = np.zeros(u.shape + (3,), np.uint8)
        rgb[:, :, 0] = (u >> 8).astype(np.uint8)
        rgb[:, :, 1] = (u & 0xFF).astype(np.uint8)
        img = Image.fromarray(rgb)
    elif name in ("indet", "bdet", "bdetf"):
        # the inner / filament detail layer, signed, over +-SPAN, 16-bit packed
        from .render import INNER_DETAIL_SPAN, BAND_DETAIL_SPAN
        _sp = INNER_DETAIL_SPAN if name == "indet" else BAND_DETAIL_SPAN
        q = np.clip(np.asarray(arr, np.float32) / _sp, -1.0, 1.0) * 0.5 + 0.5
        u = np.rint(q * 65535.0).astype(np.uint16)
        rgb = np.zeros(u.shape + (3,), np.uint8)
        rgb[:, :, 0] = (u >> 8).astype(np.uint8)
        rgb[:, :, 1] = (u & 0xFF).astype(np.uint8)
        img = Image.fromarray(rgb)
    elif name == "st_C":
        img = Image.fromarray((np.clip(np.asarray(arr, np.float32) / 4.0, 0, 1) * 255 + 0.5)
                              .astype(np.uint8))
    elif name == "contact":
        img = Image.fromarray((np.clip(arr, 0, 1) * 255).astype(np.uint8))
    else:
        img = Image.fromarray((np.clip(arr, 0, 1) * 255).astype(np.uint8))
    buf = io.BytesIO()
    img.save(buf, format="PNG")
    buf.seek(0)
    return send_file(buf, mimetype="image/png")


# ---------- settings files ----------
#
# A look is worth carrying from one bracket to another, and worth sending to
# somebody else -- which is a different problem from saving your own work. An
# exchanged file arrives on a machine whose sliders are somewhere else entirely,
# written by a build that may not be this one, from a camera that is not this
# one. So three rules, all of them about the recipient rather than the author:
#
#   1. A KEY THE FILE DOES NOT MENTION IS RESET TO ITS DEFAULT, not left at
#      whatever the recipient happened to have. Otherwise the same file lands
#      differently on two machines and neither of them is wrong, which is the
#      opposite of exchangeable.
#   2. A KEY THIS BUILD DOES NOT KNOW IS IGNORED, quietly. A file from a later
#      version has to load on an earlier one; the alternative is that a preset
#      library rots every release.
#   3. THE SETTINGS TIED TO ONE IMAGE'S GEOMETRY ARE SKIPPED unless the file
#      came from the folder it is being loaded into. See render._LOCAL_KEYS.
#
# The file records the version that wrote it because at least one slider has
# changed meaning: hillGain was renormalised in 0.22.84 and a value from before
# that means something else. The load says so rather than silently obeying.
_SETTINGS_EXT = ".efsettings.json"

# THE PROCESSING SETTINGS, which are a different animal from the render ones.
#
# A render setting moves the picture as you drag it. Every one of these decides
# how the data was STACKED, so changing one is worth a re-run and nothing else.
# That difference is the whole reason they are a separate block in the file and
# are never applied silently on load: the page sets the controls and says a
# re-stack is needed, and the user presses Start when they mean to.
#
# Saved from what the RUN ACTUALLY USED, read back out of the work directory,
# not from where the toolbar happens to be sitting. If someone changes a
# selector and saves without re-running, the toolbar and the picture disagree,
# and the file should describe the picture -- that is what makes it worth
# reloading. `from` in the saved block says which of the two it got.
#
# id -> (default, the element id on the page)
_PROCESSING_KEYS = {
    "align_corr":    ("semi",      "alignCorr"),
    "intra_lock":    ("corona",    "intraLock"),
    "colour_planes": ("auto",      "colourPlanes"),
    "partialconv":   ("True",      "partialConv"),
    "align_filter":  ("isotropic", "alignFilter"),
    "stack_combine": ("mean",      "stackCombine"),
    "bracket_mode":  ("tiers",     "bracketMode"),
    "fnrgf_preset":  ("ours",      "fnrgfPreset"),
    "photo_solve":   ("chain",     "photoSolve"),
    "tier_mode":     ("exposure",  "tierMode"),
    "align_mode":    ("auto",      "alignMode"),
    "feather":       ("plain",     "feather"),
    "wb_source":     ("camera",    "wb"),
    "photometry":    ("linfit",    "photometry"),
    "frames":        ("all",       "frames"),
    "denoise":       ("off",       "denoise"),
    # THE SIMPLE STACK FLAG. Without it a saved recipe could not turn the
    # fallback back on: the recipe is built from what the run reported, and the
    # simple path reported only denoise and fnrgf_preset.
    "simple":        ("False",     "simpleStack"),
}


def _run_processing():
    """What the loaded folder was actually stacked with, if anything was."""
    try:
        o = (json.load(open(os.path.join(workdir(STATE["folder"]),
                                         "report.json"))) or {}).get("options")
        if isinstance(o, dict):
            return {k: o[k] for k in _PROCESSING_KEYS if k in o}, "the last run"
    except Exception:
        pass
    return {}, None


def _folder_camera():
    """The camera line from this folder's run report, for provenance only."""
    try:
        rp = os.path.join(workdir(STATE["folder"]), "report.json")
        return (json.load(open(rp)) or {}).get("camera")
    except Exception:
        return None


@app.post("/api/settings/save")
def settings_save():
    if STATE["folder"] is None:
        return jsonify({"ok": False, "error": "load a folder first"}), 400
    from . import __version__
    data = request.json if request.is_json else {}
    params = data.get("params") or {}
    name = str(data.get("name") or "settings").strip() or "settings"
    name = os.path.basename(name)
    if name.endswith(_SETTINGS_EXT):
        name = name[:-len(_SETTINGS_EXT)]
    out = os.path.join(STATE["folder"], "eclipseforge_output")
    os.makedirs(out, exist_ok=True)
    path = os.path.join(out, name + _SETTINGS_EXT)
    # Only keys this build actually renders with. A stray key from the page
    # would otherwise be preserved for ever and mean nothing to anyone.
    keep = {k: params[k] for k in _RENDER_DEFAULTS if k in params}
    _por, _loc = classify_settings(keep)
    ly = STATE.get("layers")
    doc = {
        "eclipseforge_settings": 1,
        "version": __version__,
        "source": {
            "folder": STATE["folder"],
            "camera": _folder_camera(),
            "note": data.get("note") or "",
        },
        "portable": {k: keep[k] for k in sorted(_por)},
        "local": {k: keep[k] for k in sorted(_loc)},
    }

    try:
        if ly is not None:
            doc["source"]["lunar_radius_px"] = round(float(ly.R), 1)
    except Exception:
        pass
    try:
        json.dump(doc, open(path, "w", encoding="utf-8"), indent=1)
    except OSError as e:
        return jsonify({"ok": False, "error": str(e)}), 400
    return jsonify({"ok": True, "path": path,
                    "n_portable": len(doc["portable"]),
                    "n_local": len(doc["local"])})


@app.post("/api/settings/load")
def settings_load():
    from . import __version__
    data = request.json if request.is_json else {}
    path = os.path.expanduser(str(data.get("path") or "").strip())
    if not path or not os.path.isfile(path):
        return jsonify({"ok": False, "error": f"no such file: {path}"}), 400
    try:
        doc = json.load(open(path, encoding="utf-8"))
    except Exception as e:
        return jsonify({"ok": False, "error": f"not readable as JSON ({e})"}), 400
    if not isinstance(doc, dict):
        return jsonify({"ok": False, "error": "not a settings file"}), 400
    if "eclipseforge_processing" in doc:
        # a processing file has no render keys: it used to come back as "ok"
        # with every slider reset (audit 2026-10-10)
        return jsonify({"ok": False, "error": "that is a processing-options file, "
                                              "not a render-settings file (use the other Load button)"}), 400

    mode = "import" if (STATE.get("layers") is not None
                        and getattr(STATE["layers"], "mode", "") == "import")         else "raw"
    base = dict(defaults_for(mode))

    # Two shapes are accepted: this file, and the flat .params.json sidecar the
    # export has always written beside an image -- which is the same dict without
    # the provenance, and which people already have lying next to every render
    # they liked. Refusing to read it would be pedantry.
    if "eclipseforge_settings" in doc:
        src = dict(doc.get("portable") or {})
        loc = dict(doc.get("local") or {})
        wrote = str(doc.get("version") or "?")
        origin = (doc.get("source") or {}).get("folder") or ""
        note = (doc.get("source") or {}).get("note") or ""
        shape = "settings file"
    else:
        _p, _l = classify_settings(doc)
        src = {k: doc[k] for k in _p if k in doc}
        loc = {k: doc[k] for k in _l if k in doc}
        wrote, origin, note = "?", "", ""
        shape = "exported .params.json"
        # a sidecar in this folder's own eclipseforge_output describes this
        # very stack, so its disc, black point and ring settings apply too
        if STATE.get("folder") and os.path.normcase(os.path.abspath(os.path.dirname(path))) == \
                os.path.normcase(os.path.abspath(os.path.join(STATE["folder"], "eclipseforge_output"))):
            origin = STATE["folder"]

    same_stack = bool(origin) and STATE["folder"] and \
        os.path.normcase(os.path.abspath(origin)) == \
        os.path.normcase(os.path.abspath(STATE["folder"]))

    applied, unknown = {}, []
    for k, v in list(src.items()) + (list(loc.items()) if same_stack else []):
        if k == "orient":
            continue
        if k not in base:
            unknown.append(k)
            continue
        if isinstance(base[k], list):          # the tone curve: [[in, out], ...]
            try:
                pts = [[min(1.0, max(0.0, float(a))), min(1.0, max(0.0, float(b)))] for a, b in v]
                applied[k] = pts if pts else list(base[k])
            except (TypeError, ValueError):
                unknown.append(k)
            continue
        try:
            applied[k] = float(v)
        except (TypeError, ValueError):
            unknown.append(k)
    out = dict(base)          # rule 1: anything unmentioned goes back to default
    out.update(applied)

    notes = []
    if shape == "exported .params.json":
        notes.append("read from an export sidecar, which carries no version — "
                     "check the result rather than trusting it")
    if wrote not in ("?", __version__):
        notes.append(f"written by {wrote}, this is {__version__}")
        if _older_than(wrote, "0.22.84") and "hillGain" in applied:
            notes.append("Amplification was renormalised in 0.22.84 (against "
                         "the mask's structure rather than its total spread), "
                         "so this file's value means something else here — "
                         "expect to re-set it")
    if loc and not same_stack:
        notes.append("%d setting(s) tied to the source image were skipped: %s"
                     % (len(loc), ", ".join(sorted(loc))))
    if unknown:
        notes.append("%d unknown key(s) ignored: %s"
                     % (len(unknown), ", ".join(sorted(unknown)[:6])))
    return jsonify({"ok": True, "params": out, "notes": notes, "note": note,
                    "same_stack": same_stack, "n": len(applied),
                    "file": os.path.basename(path)})




# ---- block 3b: the Versions strip ------------------------------------------
_VERSIONS_NOTES = ".efversions.json"      # {image file name: note}, in eclipseforge_output


def _outdir():
    return os.path.join(STATE["folder"], "eclipseforge_output") if STATE.get("folder") else None


def _thumb_path(img):
    return os.path.splitext(img)[0] + "_thumb.jpg"


def _make_thumb(img, thumb):
    """A thumbnail for an export made before thumbnails existed: read the
    file back (16-bit TIFF through tifffile, the rest through PIL)."""
    from .render import thumbnail
    ext = os.path.splitext(img)[1].lower()
    if ext in (".tif", ".tiff"):
        import tifffile
        a = tifffile.imread(img)
    else:
        from PIL import Image
        a = np.asarray(Image.open(img))
    a = np.asarray(a)
    if a.ndim == 3 and a.shape[2] == 4:
        a = a[:, :, :3]
    scale = 65535.0 if a.dtype == np.uint16 else 255.0
    thumbnail(a.astype(np.float32) / scale, thumb)


def _read_notes(out):
    try:
        d = json.load(open(os.path.join(out, _VERSIONS_NOTES), encoding="utf-8"))
        return d if isinstance(d, dict) else {}
    except Exception:
        return {}


@app.get("/api/exports")
def list_exports():
    """The exports of the loaded folder, newest first: every image with a
    .params.json sidecar beside it (the sidecar is what makes it loadable)."""
    out = _outdir()
    if not out or not os.path.isdir(out):
        return jsonify({"ok": True, "items": []})
    notes = _read_notes(out)
    items = []
    for fn in os.listdir(out):
        if not fn.endswith(".params.json"):
            continue
        img = os.path.join(out, fn[:-len(".params.json")])
        if not os.path.isfile(img):
            continue
        try:
            st = os.stat(img)
        except OSError:
            continue
        base = os.path.basename(img)
        view = "composite"
        for v in ("_mgn", "_fnrgf", "_nafe", "_inner", "_pellett", "_promgate", "_hill"):
            if v in base:
                view = v[1:]
        items.append({"file": base, "mtime": st.st_mtime, "size": st.st_size,
                      "thumb": os.path.isfile(_thumb_path(img)),
                      "note": notes.get(base, ""), "view": view,
                      "sidecar": os.path.join(out, fn)})
    items.sort(key=lambda d: -d["mtime"])
    return jsonify({"ok": True, "items": items[:60], "dir": out})


@app.get("/api/thumb")
def get_thumb():
    out = _outdir()
    fn = os.path.basename(request.args.get("f", ""))
    if not out or not fn:
        return jsonify({"ok": False, "error": "no folder"}), 400
    img = os.path.join(out, fn)
    if not os.path.isfile(img):
        return jsonify({"ok": False, "error": "no such export"}), 404
    th = _thumb_path(img)
    if not os.path.isfile(th):
        try:
            _make_thumb(img, th)
        except Exception as e:
            return jsonify({"ok": False, "error": str(e)}), 500
    return send_file(th, mimetype="image/jpeg", max_age=0)


@app.post("/api/exports/note")
def set_export_note():
    out = _outdir()
    data = request.json if request.is_json else {}
    fn = os.path.basename(str(data.get("file") or ""))
    if not out or not fn:
        return jsonify({"ok": False, "error": "no folder"}), 400
    notes = _read_notes(out)
    note = str(data.get("note") or "").strip()
    if note:
        notes[fn] = note
    else:
        notes.pop(fn, None)
    try:
        os.makedirs(out, exist_ok=True)
        json.dump(notes, open(os.path.join(out, _VERSIONS_NOTES), "w", encoding="utf-8"), indent=1)
    except OSError as e:
        return jsonify({"ok": False, "error": str(e)}), 400
    return jsonify({"ok": True})


_PROCESSING_EXT = ".efprocess.json"


@app.post("/api/processing/save")
def processing_save():
    """The stack recipe, in its own file. See _PROCESSING_KEYS for why."""
    if STATE["folder"] is None:
        return jsonify({"ok": False, "error": "load a folder first"}), 400
    from . import __version__
    data = request.json if request.is_json else {}
    name = str(data.get("name") or "processing").strip() or "processing"
    name = os.path.basename(name)
    if name.endswith(_PROCESSING_EXT):
        name = name[:-len(_PROCESSING_EXT)]
    proc, src = _run_processing()
    if not proc:
        sent = data.get("processing") or {}
        proc = {k: sent[k] for k in _PROCESSING_KEYS if k in sent}
        src = "the toolbar (this folder has not been stacked yet)"
    if not proc:
        return jsonify({"ok": False,
                        "error": "nothing to save: no run in this folder and "
                                 "no settings sent"}), 400
    out = os.path.join(STATE["folder"], "eclipseforge_output")
    os.makedirs(out, exist_ok=True)
    path = os.path.join(out, name + _PROCESSING_EXT)
    doc = {"eclipseforge_processing": 1,
           "version": __version__,
           "source": {"folder": STATE["folder"], "camera": _folder_camera(),
                      "from": src, "note": data.get("note") or ""},
           "processing": {k: str(proc[k]) for k in sorted(proc)}}
    try:
        json.dump(doc, open(path, "w", encoding="utf-8"), indent=1)
    except OSError as e:
        return jsonify({"ok": False, "error": str(e)}), 400
    return jsonify({"ok": True, "path": path, "n": len(doc["processing"]),
                    "from": src})


@app.post("/api/processing/load")
def processing_load():
    """Return the recipe and say which keys differ from this folder's last run.

    NOT APPLIED, and not defaulted. The look file's rule -- a key it does not
    mention goes back to its default -- is right there and wrong here: these
    decide how the data is STACKED, so quietly resetting one and having the next
    Start re-stack differently is the sort of thing nobody notices until the
    picture has changed. A key the file omits is left where it is.
    """
    from . import __version__
    data = request.json if request.is_json else {}
    path = os.path.expanduser(str(data.get("path") or "").strip())
    if not path or not os.path.isfile(path):
        return jsonify({"ok": False, "error": f"no such file: {path}"}), 400
    try:
        doc = json.load(open(path, encoding="utf-8"))
    except Exception as e:
        return jsonify({"ok": False, "error": f"not readable as JSON ({e})"}), 400
    if not isinstance(doc, dict) or "eclipseforge_processing" not in doc:
        return jsonify({"ok": False,
                        "error": "not a processing file (a look is saved and "
                                 "loaded with the Settings buttons above)"}), 400
    proc, unknown = {}, []
    for k, v in (doc.get("processing") or {}).items():
        (proc.__setitem__(k, str(v)) if k in _PROCESSING_KEYS
         else unknown.append(k))
    notes = []
    wrote = str(doc.get("version") or "?")
    if wrote not in ("?", __version__):
        notes.append(f"written by {wrote}, this is {__version__}")
    if unknown:
        notes.append("%d key(s) this build does not know, ignored: %s"
                     % (len(unknown), ", ".join(sorted(unknown)[:6])))
    ran, _ = _run_processing()
    changed = sorted(k for k, v in proc.items()
                     if v != str(ran.get(k, _PROCESSING_KEYS[k][0])))
    return jsonify({"ok": True, "processing": proc, "notes": notes,
                    "note": (doc.get("source") or {}).get("note") or "",
                    "from": (doc.get("source") or {}).get("from"),
                    "changed": changed, "stacked": bool(ran),
                    "ids": {k: _PROCESSING_KEYS[k][1] for k in proc},
                    "n": len(proc), "file": os.path.basename(path),
                    "camera": (doc.get("source") or {}).get("camera")})


def _older_than(a, b):
    """Version compare on dotted integers; unparseable sorts as 'older'."""
    def t(v):
        try:
            return tuple(int(x) for x in str(v).split("."))
        except ValueError:
            return ()
    ta, tb = t(a), t(b)
    return (ta < tb) if (ta and tb) else True


# THE 1:1 CROP (0.25.1, 2026-10-11). The preview is decimated, and
# decimation folds the finest things -- the 1.5 and 3 px sharpening bands,
# the grain, the inner-detail and filament layers -- into one noise texture;
# until now the only way to see the export's pixels was the export. This
# renders a window of the full-resolution picture with the current sliders
# through the export path (crop_layers + render), in a second or two. The
# picture-wide numbers (inner corona level's ring median, contain whites,
# Moon disc median, the structure slider's ring means) come from a preview
# render of the whole frame (`pin`), cached per slider set, so the crop is
# the export's pixels and not a crop-normalised picture.
CROP_MARGIN = 320            # px of context rendered around the window, cropped off after


@app.post("/api/crop")
def do_crop():
    if not request.is_json:
        return jsonify({"ok": False, "error": "expected a JSON body"}), 415
    _b = _busy()
    if _b:
        return jsonify({"ok": False, "error": f"{_b} is still running — wait for it to finish"}), 400
    ly = STATE["layers"]
    if ly is None:
        return jsonify({"ok": False, "error": "run the pipeline first"}), 400
    data = request.json or {}
    params = dict(data.get("params") or {})
    view = str(data.get("view") or "composite")
    try:
        x, y = int(round(float(data["x"]))), int(round(float(data["y"])))
    except Exception:
        return jsonify({"ok": False, "error": "x and y (full-resolution px) are needed"}), 400
    size = int(min(max(int(data.get("size", 600) or 600), 128), 1200))
    H, W = ly.full["bg"].shape
    x = min(max(x, 0), W - 1); y = min(max(y, 0), H - 1)
    from .render import crop_layers as _crop, apply_orient as _ao
    import time as _t
    t0 = _t.time()
    try:
        # the picture-wide numbers, from a preview render of the whole frame;
        # one cached entry, keyed on the sliders that reach the render
        key = json.dumps({k: v for k, v in params.items() if k != "orient"}, sort_keys=True, default=str)
        cache = STATE.setdefault("crop_pin", {})
        if cache.get("key") != key or cache.get("lyid") != id(ly):
            pin = {}
            render(ly, params, preview=True, view="composite", pin=pin)
            cache.clear(); cache.update(key=key, pin=pin, lyid=id(ly))
        pin = dict(cache["pin"])
        h = size // 2
        y0, y1 = y - h - CROP_MARGIN, y - h + size + CROP_MARGIN
        x0, x1 = x - h - CROP_MARGIN, x - h + size + CROP_MARGIN
        C = _crop(ly, y0, y1, x0, x1)
        cy0, cy1, cx0, cx1 = C.crop_box
        rgb = render(C, params, preview=False, view=view, pin=pin)
        wy0, wx0 = max(y - h, 0) - cy0, max(x - h, 0) - cx0
        rgb = rgb[wy0:wy0 + size, wx0:wx0 + size]
        if rgb.ndim == 2:
            rgb = np.repeat(rgb[:, :, None], 3, axis=2)
        rgb = _ao(np.ascontiguousarray(rgb, np.float32), str(params.get("orient", "") or ""))
        a8 = (np.clip(rgb, 0, 1) * 255.0 + 0.5).astype(np.uint8)
        from PIL import Image
        buf = io.BytesIO()
        Image.fromarray(a8, "RGB").save(buf, format="PNG", compress_level=1)
        buf.seek(0)
        resp = send_file(buf, mimetype="image/png")
        resp.headers["X-Crop-Box"] = json.dumps([int(max(y - h, 0)), int(max(x - h, 0)),
                                                 int(a8.shape[0]), int(a8.shape[1])])
        resp.headers["X-Crop-Seconds"] = f"{_t.time() - t0:.1f}"
        return resp
    except Exception as e:
        return jsonify({"ok": False, "error": f"crop failed: {e}"}), 500


def _polar_angles(po):
    """The by-hand solar north (deg cw from image up, or None) and the fine
    turn (deg) from the polar options; blank or bad text means not set."""
    def num(v):
        try:
            v = str(v).strip().replace(",", ".")
            return float(v) if v else None
        except (TypeError, ValueError):
            return None
    m = num(po.get("manual")); off = num(po.get("offset")) or 0.0
    return m, off


# EXPORT THUMBNAILS (tab 4, 2026-11-08): small pictures of what each export
# product will look like, rendered from the PREVIEW layers (a quarter of the
# resolution) with the page's current settings, so they cost a fraction of an
# export. Kept in memory by request, a few dozen at most.
_XTHUMB = {}


def _thumb_jpeg(arr, width=520, rgba_over=None):
    a = np.clip(np.asarray(arr, np.float32), 0, 1)
    if a.ndim == 3 and a.shape[2] == 1:
        a = a[:, :, 0]
    img = Image.fromarray((a * 255 + 0.5).astype(np.uint8))
    if img.mode != "RGB":
        img = img.convert("RGB")
    if img.width > width:
        img = img.resize((width, max(1, int(round(img.height * width / img.width)))), Image.LANCZOS)
    b = io.BytesIO()
    img.save(b, "JPEG", quality=85)
    return b.getvalue()


@app.post("/api/xthumb")
def export_thumb():
    import hashlib
    from .render import render, apply_orient, plain_stack
    ly = STATE["layers"]
    if ly is None:
        return jsonify({"ok": False, "error": "no result yet"}), 400
    data = request.json if request.is_json else {}
    kind = str(data.get("kind") or "composite")
    params = dict(data.get("params") or {})
    _or = str(data.get("orient", "") or "").strip().lower()
    params["orient"] = _or if _or in ("flipv", "fliph", "180", "cw", "ccw") else ""
    key = hashlib.sha1(json.dumps([STATE.get("folder"), id(ly), kind, params, data.get("polar")],
                                  sort_keys=True, default=str).encode()).hexdigest()
    if key in _XTHUMB:
        return send_file(io.BytesIO(_XTHUMB[key]), mimetype="image/jpeg", max_age=0)
    try:
        wd = workdir(STATE["folder"])
        q = ly.prev_decim
        if kind == "plain":
            geo = json.load(open(os.path.join(wd, "geometry.json")))
            arr = apply_orient(plain_stack(wd, geometry=geo, step=max(2, q)), params["orient"])
            jpg = _thumb_jpeg(arr)
        elif kind == "hdr32":
            # AS THE FILE OPENS, NO STRETCH (2026-11-08): the merge scaled as
            # the 32-bit export scales it (1.0 = the corona's 99.99th
            # percentile outside the Moon), clipped at 1, and shown the way a
            # colour-managed program shows a file with a linear profile --
            # sRGB-encoded for the screen, no tone curve on top
            from .pipeline import load_big
            geo = json.load(open(os.path.join(wd, "geometry.json")))
            st = max(2, q)
            a = np.asarray(load_big(os.path.join(wd, "hdr_rgb.npy"))[::st, ::st], np.float32)
            a = np.nan_to_num(np.maximum(a, 0.0), nan=0.0, posinf=0.0)
            Y = 0.2126 * a[:, :, 0] + 0.7152 * a[:, :, 1] + 0.0722 * a[:, :, 2]
            yy = np.arange(a.shape[0], dtype=np.float32)[:, None] - np.float32(float(geo["cy"]) / st)
            xx = np.arange(a.shape[1], dtype=np.float32)[None, :] - np.float32(float(geo["cx"]) / st)
            sel = (Y > 0) & ((yy * yy + xx * xx) > np.float32((1.02 * float(geo["R"]) / st) ** 2))
            ref = float(np.percentile(Y[sel], 99.99)) if sel.any() else (float(Y.max()) or 1.0)
            x = np.clip(a / np.float32(max(ref, 1e-12)), 0.0, 1.0)
            x = np.where(x <= 0.0031308, 12.92 * x, 1.055 * np.power(x, 1 / 2.4) - 0.055)
            jpg = _thumb_jpeg(apply_orient(x.astype(np.float32), params["orient"]))
        elif kind.startswith("view:"):
            v = kind[5:]
            arr = render(ly, params, preview=True, view=v)
            arr = apply_orient(arr, params["orient"])
            jpg = _thumb_jpeg(arr if v in ("composite", "moon") else arr[..., 0] if arr.ndim == 3 else arr)
        elif kind in ("polar", "polarlayers"):
            from . import polarmap as _pm
            po = data.get("polar") or {}
            rgb = render(ly, params, preview=True, view="composite")
            cy, cx, R = ly.geometry(q)
            pmap, meta, _o = _pm.polar_of(
                rgb, cy, cx, R, wd, STATE["folder"],
                north="celestial" if po.get("north") == "celestial" else "solar",
                outer="frame" if po.get("outer") == "frame" else "circle",
                radial="log" if po.get("radial") == "log" else "linear",
                centre=(str(po.get("centre") or "S").upper() if str(po.get("centre") or "S").upper()
                        in ("N", "E", "S", "W") else "S"),
                mirrored=bool(po.get("mirrored", False)),
                gray=int(params.get("colourMode", 0) or 0) == 1, unit_scale=float(q),
                manual_cw=_polar_angles(po)[0], offset_deg=_polar_angles(po)[1])
            # THE EXPORT'S SHAPE: a linear map has one row per pixel of radius,
            # so on the preview (q times fewer pixels) it came out q times
            # flatter than the file (2026-11-08, "too slim"). Its columns are
            # averaged in blocks of q instead, which gives the file's shape. The
            # log map is conformal and already has it.
            if po.get("radial") != "log" and q > 1 and pmap.shape[1] >= 2 * q:
                _w = pmap.shape[1] // q * q
                pmap = pmap[:, :_w].reshape(pmap.shape[0], _w // q, q, pmap.shape[2]).mean(axis=2)
            if kind == "polar" and data.get("format") in ("tif16", "tif8"):
                jpg = _thumb_jpeg(pmap)                      # the TIFF is the clean map
            else:
                im = _pm.annotate(pmap, meta, grid=bool(po.get("grid", True)),
                                  earth=bool(po.get("earth", True)))
                jpg = _thumb_jpeg(np.asarray(im, np.float32) / 255.0)
        else:                                                # the picture
            arr = render(ly, params, preview=True, view="composite")
            jpg = _thumb_jpeg(apply_orient(arr, params["orient"]))
    except Exception as e:
        return jsonify({"ok": False, "error": f"thumbnail: {e}"}), 500
    if len(_XTHUMB) > 40:
        _XTHUMB.pop(next(iter(_XTHUMB)))
    _XTHUMB[key] = jpg
    return send_file(io.BytesIO(jpg), mimetype="image/jpeg", max_age=0)


@app.post("/api/export")
def do_export():
    # BUSY GUARD (audit 2026-10-10): without it a second export orphaned the
    # first (its progress and cancel went to the new one), and a Start or
    # Clear cache then rewrote the .npy files under the first export's memory
    # maps -- a SIGBUS that kills the whole server.
    _b = _busy()
    if _b:
        return jsonify({"ok": False, "error": f"{_b} is still running — wait for it to finish"}), 400
    ly = STATE["layers"]
    if ly is None:
        return jsonify({"ok": False, "error": "run the pipeline first"}), 400
    data = request.json
    fmt = data.get("format", "tif16")
    view = data.get("view", "composite")
    # the extra products of tab 4 (2026-11-08): "plain" is the merge with one
    # neutral stretch (render.plain_stack), "polarlayers" the annotated polar
    # map as a base layer and a transparent overlay (polarmap.export_layers)
    product = str(data.get("product") or "picture")
    if product not in ("picture", "plain", "polarlayers"):
        return jsonify({"ok": False, "error": f"unknown product {product!r}"}), 400
    size = data.get("size", "full")
    prom_layer = bool(data.get("promLayer", False))
    moon_layer = bool(data.get("moonLayer", False))
    _po = data.get("polar") or {}
    polar = bool(_po.get("on", False))
    polar_north = "celestial" if _po.get("north") == "celestial" else "solar"
    polar_outer = "frame" if _po.get("outer") == "frame" else "circle"
    polar_grid = bool(_po.get("grid", True))
    polar_earth = bool(_po.get("earth", True))
    polar_radial = "log" if _po.get("radial") == "log" else "linear"
    polar_centre = str(_po.get("centre") or "S").upper()
    polar_centre = polar_centre if polar_centre in ("N", "E", "S", "W") else "S"
    polar_mirror = bool(_po.get("mirrored", False))
    polar_only = polar and bool(_po.get("only", False))
    polar_manual, polar_offset = _polar_angles(_po)
    # carousel tiles of the map (2026-11-10): n tiles of a social size
    try:
        _tn = int(_po.get("tiles") or 0)
    except (TypeError, ValueError):
        _tn = 0
    polar_tiles = ({"n": _tn, "fmt": str(_po.get("tileFmt") or "3:4")}
                   if 2 <= _tn <= 12 else None)
    params = data.get("params", {})
    # Purely cosmetic and applied last, so it needs no re-run and no re-render
    # of anything cached -- see render.apply_orient.
    _or = str(data.get("orient", "") or "").strip().lower()
    params = dict(params, orient=(_or if _or in ("flipv", "fliph", "180", "cw", "ccw") else ""))
    if fmt not in ("tif16", "tif8", "png", "jpg", "hdr32"):
        return jsonify({"ok": False, "error": f"unknown format {fmt!r}"}), 400
    ext = {"tif16": ".tif", "tif8": ".tif", "png": ".png", "jpg": ".jpg", "hdr32": ".tif"}[fmt]
    name = os.path.basename(str(data.get("name") or "").strip()) or "eclipseforge_render"
    if fmt == "hdr32":
        # the merge, not the render: the view, the look and the layers do not apply
        name = name.replace("render", "merge") if "render" in name else name + "_merge"
        view = "composite"
        prom_layer = moon_layer = polar = False
        product = "picture"
    if product != "picture":
        view = "composite"
        prom_layer = moon_layer = polar = polar_only = False
        name += "_plainstack" if product == "plain" else ""
    if view != "composite":
        name += "_" + view
    if size == "half":
        name += "_50pct"
    outdir = os.path.join(STATE["folder"], "eclipseforge_output")
    os.makedirs(outdir, exist_ok=True)
    base = os.path.join(outdir, name)
    _suf = {"tif16": "_16bit", "hdr32": "_32bit_linear"}.get(fmt, "")
    path = base + _suf + ext
    i = 1
    while os.path.exists(path):
        path = base + _suf + f"_{i}" + ext
        i += 1
    prog = Progress()
    STATE["progress"] = prog
    prog.log(f"exporting {fmt} ...", 0.1)

    if product == "polarlayers":
        _k, _sfx = 0, ""
        while (os.path.exists(base + "_polar_base_16bit" + _sfx + ".tif")
               or os.path.exists(base + "_polar_overlay_16bit" + _sfx + ".tif")):
            _k += 1
            _sfx = f"_{_k}"
        path = base + "_polar_base_16bit" + _sfx + ".tif"
        _ovp = base + "_polar_overlay_16bit" + _sfx + ".tif"

    def work():
        if product == "plain":
            try:
                from .render import plain_stack, write_picture, apply_orient
                prog.log("writing the plain HDR stack (the merge, one neutral stretch, no look) ...", 0.3)
                _wd = workdir(STATE["folder"])
                _geo = json.load(open(os.path.join(_wd, "geometry.json")))
                _a = apply_orient(plain_stack(_wd, size=size, geometry=_geo), params.get("orient", ""))
                write_picture(_a, fmt, path, gray=int(params.get("colourMode", 0) or 0) == 1,
                              description="eclipseforgehdr plain HDR stack: the merge with one "
                                          "logarithmic luminance stretch; no detail layers, no look")
                try:
                    from .render import thumbnail
                    thumbnail(_a, _thumb_path(path))
                except Exception:
                    pass
                json.dump(params, open(path + ".params.json", "w"), indent=1)
                prog.log(f"[ok] saved {path}", 1.0)
            except Exception as e:
                prog.error = f"plain stack: {e}"
            prog.done = True
            return
        if product == "polarlayers":
            try:
                from . import polarmap as _pm
                prog.log("rendering the picture for the polar map ...", 0.3)
                _rgb = render(ly, params, preview=False, view="composite")
                _gcy, _gcx, _gR = ly.geometry(1)
                prog.log("unwrapping ...", 0.6)
                pmap, meta, _o = _pm.polar_of(_rgb, _gcy, _gcx, _gR, workdir(STATE["folder"]),
                                              STATE["folder"], north=polar_north, outer=polar_outer,
                                              step=2.0 if size == "half" else 1.0, radial=polar_radial,
                                              centre=polar_centre, mirrored=polar_mirror, manual_cw=polar_manual, offset_deg=polar_offset,
                                              gray=int(params.get("colourMode", 0) or 0) == 1)
                del _rgb
                prog.log("writing the base and the overlay ...", 0.85)
                W, H = _pm.export_layers(pmap, meta, path, _ovp, grid=polar_grid, earth=polar_earth)
                json.dump(params, open(path + ".params.json", "w"), indent=1)
                prog.log(f"[ok] saved {path} and {os.path.basename(_ovp)} ({W} x {H}; put the "
                         f"overlay over the base, blend mode Normal)", 1.0)
            except Exception as e:
                prog.error = f"polar layers: {e}"
            prog.done = True
            return
        if fmt == "hdr32":
            try:
                from .render import export_hdr32
                prog.log("writing the merged HDR (32-bit float, linear) ...", 0.3)
                _wd = workdir(STATE["folder"])
                _geo = json.load(open(os.path.join(_wd, "geometry.json")))
                _ref = export_hdr32(_wd, path, orient=params.get("orient", ""), size=size,
                                    geometry=_geo)
                prog.log(f"[ok] saved {path} (scene-linear; 1.0 = the corona's 99.99th "
                         f"luminance percentile, {_ref:.4g} in the stack's units)", 1.0)
                prog.done = True
            except Exception as e:
                prog.error = str(e)
                prog.done = True
            return
        if polar_only and view == "composite":
            # ONLY THE POLAR MAP (2026-11-08): the picture is rendered at full
            # resolution as for an export but not written; the map pair gets
            # its own number so it never overwrites an earlier one
            try:
                from . import polarmap as _pm
                prog.log("rendering the picture for the polar map ...", 0.3)
                _rgb = render(ly, params, preview=False, view="composite")
                # ONE FILE, in the format chosen (2026-11-08: the polar-only
                # export wrote a TIFF beside the JPG). TIFF: the clean map,
                # data only; JPG / PNG: the annotated map with axes and labels.
                _tail = {"tif16": "_polar_16bit.tif", "tif8": "_polar_8bit.tif",
                         "png": "_polar.png", "jpg": "_polar.jpg"}[fmt]
                _k, _out = 0, base + _tail
                while os.path.exists(_out):
                    _k += 1
                    _out = base + _tail.replace(".", f"_{_k}.", 1)
                _pt = _out if fmt in ("tif16", "tif8") else None
                _pj = None if _pt else _out
                prog.log("writing polar map ...", 0.8)
                _gcy, _gcx, _gR = ly.geometry(1)
                _line = _pm.export(_rgb, _gcy, _gcx, _gR, workdir(STATE["folder"]),
                                   STATE["folder"], _pt, _pj, north=polar_north,
                                   outer=polar_outer, grid=polar_grid, earth=polar_earth,
                                   radial=polar_radial, centre=polar_centre, mirrored=polar_mirror, manual_cw=polar_manual, offset_deg=polar_offset,
                                   step=2.0 if size == "half" else 1.0,
                                   gray=int(params.get("colourMode", 0) or 0) == 1,
                                   tif_bits=8 if fmt == "tif8" else 16, tiles=polar_tiles)
                json.dump(params, open(_out + ".params.json", "w"), indent=1)
                prog.log(f"[ok] saved {_out} ({_line})", 1.0)
                prog.done = True
            except Exception as e:
                prog.error = f"polar map: {e}"
                prog.done = True
            return
        try:
            _notes = []
            # the rendered composite is kept for the prominence layer and the quality grading
            _keep = {} if view == "composite" else None
            export(ly, params, fmt, path, view=view, size=size, notes=_notes,
                   thumb=_thumb_path(path), keep=_keep)
            json.dump(params, open(path + ".params.json", "w"), indent=1)
            # PICTURE QUALITY (lab 0.24): graded measurements of what was just written
            _pq = None
            if _keep is not None and _keep.get("rgb") is not None:
                try:
                    from . import picquality as _pqm
                    prog.log("grading picture ...", 0.75)
                    _gcy, _gcx, _gR = ly.geometry(1)
                    _pq = _pqm.grade(_keep["rgb"], _gcy, _gcx, _gR, getattr(ly, "Rmask", None))
                    _sm = _pqm.summary(_pq)
                    _lv = _pqm.level(_pq) if hasattr(_pqm, "level") else "info"
                    _notes.append({"warn": "[warn] ", "odd": "[odd] ", "ok": "[ok] "}.get(_lv, "")
                                  + _sm + " (details in report)")
                except Exception as _e:
                    _notes.append(f"[warn] picture quality not measured ({_e})")
            # the prominences as their own RGBA layer, on the same grid
            if prom_layer and view == "composite":
                try:
                    from . import promlayer as _pl
                    _wd = workdir(STATE["folder"])
                    _geo = json.load(open(os.path.join(_wd, "geometry.json")))
                    # render_16bit_30.tif -> render_prominences_16bit_30.tif (pairs by number)
                    _st = os.path.splitext(path)[0]
                    _pp = ("_prominences_16bit".join(_st.rsplit("_16bit", 1))
                           if "_16bit" in os.path.basename(_st) else _st + "_prominences_16bit") + ".tif"
                    prog.log("building prominence layer ...", 0.8)
                    if getattr(ly, "has_promlayer", False):
                        # THE 0.24 LAYER: the RGBA file holds exactly what the
                        # preview lays over the picture, and the corona it goes
                        # on is written beside it, so the pair rebuilds the
                        # composite in Photoshop (layer on top, Normal)
                        _cb = ("_corona_16bit".join(_st.rsplit("_16bit", 1))
                               if "_16bit" in os.path.basename(_st) else _st + "_corona") + ext
                        _kc = {}
                        if _keep.get("corona") is not None:
                            # kept by render(): no second render (audit 2026-10-10)
                            export(ly, dict(params, promLayer=0.0), fmt, _cb, view=view, size=size,
                                   keep=_kc, rgb_in=_keep["corona"])
                        else:
                            export(ly, dict(params, promLayer=0.0), fmt, _cb, view=view, size=size,
                                   keep=_kc)
                        prog.log(f"[ok] saved {_cb} (corona without prominence layer)", 0.85)
                        _pl.export_layer2(ly, params, _pp, _keep.get("rgb"), _kc.get("rgb"),
                                          size=size)
                        _keep.pop("corona", None)
                        del _kc
                        _sharp = True
                    else:
                        _sharp = _pl.export_layer(_wd, _geo, _pp, orient=params.get("orient", ""),
                                                  size=size, log=lambda m: prog.log(m, None))
                    prog.log(f"[ok] saved {_pp}", 0.95)
                    if not _sharp:
                        _notes.append("[warn] prominence layer from merge only (no cached "
                                      "stack); Clear cache and re-run for "
                                      "sharp edges")
                except Exception as _e:
                    _notes.append(f"[warn] prominence layer not written ({_e})")
            # the Moon as its own RGBA layer, on the same grid (lab 0.24)
            if moon_layer and view == "composite" and _keep.get("rgb") is not None:
                try:
                    from .render import export_moon_layer
                    _st = os.path.splitext(path)[0]
                    _mp = ("_moon_16bit".join(_st.rsplit("_16bit", 1))
                           if "_16bit" in os.path.basename(_st) else _st + "_moon_16bit") + ".tif"
                    prog.log("writing Moon layer ...", 0.9)
                    export_moon_layer(ly, params, _mp, _keep["rgb"], size=size)
                    prog.log(f"[ok] saved {_mp}", 0.95)
                except Exception as _e:
                    _notes.append(f"[warn] Moon layer not written ({_e})")
            # the picture unwrapped around the Sun: position angle x radius (0.30.0)
            if polar and view == "composite" and _keep.get("rgb") is not None:
                try:
                    from . import polarmap as _pm
                    # ONE FILE BESIDE THE PICTURE, in the picture's format
                    # (2026-11-08): TIFF = the clean map, JPG / PNG = the
                    # annotated map. Named after the picture; numbered on if
                    # a polar-only export already took the name.
                    _st = os.path.splitext(path)[0]
                    if fmt == "tif16" and "_16bit" in os.path.basename(_st):
                        _o0 = "_polar_16bit".join(_st.rsplit("_16bit", 1)) + ".tif"
                    else:
                        _o0 = _st + {"tif16": "_polar_16bit.tif", "tif8": "_polar_8bit.tif",
                                     "png": "_polar.png", "jpg": "_polar.jpg"}[fmt]
                    _out, _k = _o0, 0
                    while os.path.exists(_out):
                        _k += 1
                        _out = os.path.splitext(_o0)[0] + f"_{_k}" + os.path.splitext(_o0)[1]
                    _pt = _out if fmt in ("tif16", "tif8") else None
                    _pj = None if _pt else _out
                    prog.log("writing polar map ...", 0.92)
                    _gcy, _gcx, _gR = ly.geometry(1)
                    _line = _pm.export(_keep["rgb"], _gcy, _gcx, _gR, workdir(STATE["folder"]),
                                       STATE["folder"], _pt, _pj, north=polar_north,
                                       outer=polar_outer, grid=polar_grid, earth=polar_earth, radial=polar_radial, centre=polar_centre,
                                       mirrored=polar_mirror, manual_cw=polar_manual, offset_deg=polar_offset,
                                       step=2.0 if size == "half" else 1.0,
                                       gray=int(params.get("colourMode", 0) or 0) == 1,
                                       tif_bits=8 if fmt == "tif8" else 16, tiles=polar_tiles)
                    prog.log(f"[ok] saved {os.path.basename(_out)} ({_line})", 0.96)
                except Exception as _e:
                    _notes.append(f"[warn] polar map not written ({_e})")
            try:
                from . import report as _report
                rp = os.path.join(workdir(STATE["folder"]), "report.json")
                st = json.load(open(rp)) if os.path.exists(rp) else {}
                st["params"] = params
                st["export"] = {"file": os.path.basename(path), "format": fmt,
                                "view": view, "size": size}
                # what Neutralise sky cast divides by: the far field's colour,
                # normalised to unit luminance, as the Layers measured it
                try:
                    st["bg_chroma"] = [round(float(x), 3) for x in ly.bg_chroma]
                except Exception:
                    pass
                if _pq:
                    st["picture_quality"] = _pq
                _txt = _report.build(st)
                if _pq:
                    from . import picquality as _pqm
                    _txt += "\n\n" + _pqm.text(_pq)
                open(os.path.splitext(path)[0] + "_report.txt", "w").write(_txt + "\n")
            except Exception:
                pass
            prog.log(f"[ok] saved {path}", 1.0)
            # Logged AFTER the save line so it is the one the page shows: a
            # warning nobody reads is the same as no warning.
            for _n in _notes:
                prog.log(_n, 1.0)
            prog.done = True
        except Exception as e:
            prog.error = str(e)
            prog.done = True

    t = threading.Thread(target=work, daemon=True, name="the export")
    STATE["busy"] = t
    t.start()
    return jsonify({"ok": True, "path": (base + "_polar." + ("tif" if fmt in ("tif16", "tif8") else fmt)
                                         + " (polar map only)")
                    if (polar_only and view == "composite") else path})


# EXTRAS (tab 5, 0.30.0): products drawn from the finished picture and its
# geometry -- see extras.py. The page builds its cards from /api/extras; a
# preview is the product made on the preview layers (1/4), a generate is the
# same on the full render, written to eclipseforge_output/extras/.
_XPREV = {}


def _extra_angles(data):
    po = data.get("polar") or {}
    m, off = _polar_angles(po)
    return ("celestial" if po.get("north") == "celestial" else "solar"), m, off


@app.get("/api/extras")
def extras_list():
    from . import extras as _ex
    if not STATE.get("folder"):
        return jsonify({"ok": True, "items": _ex.catalogue("", "", False)})
    wd = workdir(STATE["folder"])
    return jsonify({"ok": True, "items": _ex.catalogue(wd, STATE["folder"], STATE["layers"] is not None)})


@app.post("/api/extra/preview")
def extra_preview():
    import hashlib
    from . import extras as _ex
    ly = STATE["layers"]
    if ly is None:
        return jsonify({"ok": False, "error": "no result yet"}), 400
    data = request.json if request.is_json else {}
    kind = str(data.get("kind") or "")
    opts = dict(data.get("opts") or {})
    params = dict(data.get("params") or {})
    if kind == "activity":
        opts["_look"] = params            # corona in motion: the look's partial-convolution settings
    # the poster's map in a second look (the standard variant of a soft look)
    map_params = dict(data.get("mapParams") or {}) if kind == "poster" else {}
    north, man, off = _extra_angles(data)
    key = hashlib.sha1(json.dumps([STATE.get("folder"), id(ly), kind, opts, params, map_params, north, man, off],
                                  sort_keys=True, default=str).encode()).hexdigest()
    if key in _XPREV:
        return send_file(io.BytesIO(_XPREV[key]), mimetype="image/jpeg", max_age=0)
    try:
        wd = workdir(STATE["folder"])
        q = ly.prev_decim
        geo = json.load(open(os.path.join(wd, "geometry.json")))
        rgb = render(ly, params, preview=True, view="composite")
        rgb_m = render(ly, map_params, preview=True, view="composite") if map_params else None
        res = _ex.make(kind, rgb, geo, wd, STATE["folder"], opts, scale=float(q), north=north,
                       manual_cw=man, offset_deg=off, preview=True, rgb_map=rgb_m)
        im = res["image"]
        if im.width > 1600:
            im = im.resize((1600, max(1, int(round(im.height * 1600 / im.width)))), Image.LANCZOS)
        b = io.BytesIO()
        im.save(b, "JPEG", quality=86)
        jpg = b.getvalue()
    except Exception as e:
        return jsonify({"ok": False, "error": f"preview: {e}"}), 500
    if len(_XPREV) > 30:
        _XPREV.pop(next(iter(_XPREV)))
    _XPREV[key] = jpg
    return send_file(io.BytesIO(jpg), mimetype="image/jpeg", max_age=0)


@app.post("/api/extra/generate")
def extra_generate():
    from . import extras as _ex
    _b = _busy()
    if _b:
        return jsonify({"ok": False, "error": f"{_b} is still running — wait for it to finish"}), 400
    ly = STATE["layers"]
    if ly is None:
        return jsonify({"ok": False, "error": "run the pipeline first"}), 400
    data = request.json if request.is_json else {}
    kind = str(data.get("kind") or "")
    if _ex.entry(kind) is None:
        return jsonify({"ok": False, "error": f"unknown product {kind!r}"}), 400
    opts = dict(data.get("opts") or {})
    params = dict(data.get("params") or {})
    if kind == "activity":
        opts["_look"] = params            # corona in motion: the look's partial-convolution settings
    map_params = dict(data.get("mapParams") or {}) if kind == "poster" else {}
    north, man, off = _extra_angles(data)
    opts["_angles"] = [north, man, off]          # film and sound build their own sky context
    name = os.path.basename(str(data.get("name") or "").strip()) or "eclipseforge_render"
    outdir = os.path.join(STATE["folder"], "eclipseforge_output", "extras")
    os.makedirs(outdir, exist_ok=True)
    base = os.path.join(outdir, f"{name}_{kind}")
    prog = Progress()
    STATE["progress"] = prog
    prog.log(f"{kind}: rendering the picture at full size ...", 0.1)

    def work():
        try:
            wd = workdir(STATE["folder"])
            if kind == "film":
                from . import film as _film
                from .render import plain_stack
                path, secs, vk = _film.make(ly, params, wd, STATE["folder"], opts, base, prog=prog,
                                            render=render, plain_stack=plain_stack)
                json.dump({"product": kind, "options": opts, "params": params},
                          open(path + ".params.json", "w"), indent=1)
                prog.log(f"[ok] saved {path} ({secs:.0f} s" + (", animated GIF: no ffmpeg found" if vk == "gif" else "") + ")", 1.0)
                prog.done = True
                return
            geo = json.load(open(os.path.join(wd, "geometry.json")))
            rgb = render(ly, params, preview=False, view="composite")
            if kind == "sound":
                from . import sound as _sound
                paths, info = _sound.make(rgb, geo, wd, STATE["folder"], opts, base, prog=prog)
                json.dump({"product": kind, "options": opts, "params": params, "info": info},
                          open(paths[0] + ".params.json", "w"), indent=1)
                more = (" and " + ", ".join(os.path.basename(p) for p in paths[1:])) if len(paths) > 1 else ""
                prog.log(f"[ok] saved {paths[0]}{more} ({info})", 1.0)
                prog.done = True
                return
            rgb_m = None
            if map_params:
                prog.log(f"{kind}: rendering the map's look at full size ...", 0.4)
                rgb_m = render(ly, map_params, preview=False, view="composite")
            prog.log(f"{kind}: drawing ...", 0.55)
            res = _ex.make(kind, rgb, geo, wd, STATE["folder"], opts, scale=1.0, north=north,
                           manual_cw=man, offset_deg=off, prog=prog, out_base=base, rgb_map=rgb_m)
            del rgb, rgb_m
            prog.log(f"{kind}: writing ...", 0.85)
            paths = _ex.write(res, base) + list(res.get("extra_paths") or [])
            json.dump(dict({"product": kind, "options": opts, "params": params, "info": res.get("info", "")},
                           **({"map_params": map_params} if map_params else {})),
                      open(paths[0] + ".params.json", "w"), indent=1)
            more = (" and " + ", ".join(os.path.basename(p) for p in paths[1:])) if len(paths) > 1 else ""
            prog.log(f"[ok] saved {paths[0]}{more} ({res.get('info', '')})", 1.0)
        except Exception as e:
            prog.error = f"{kind}: {e}"
        prog.done = True

    t = threading.Thread(target=work, daemon=True, name="the extras")
    STATE["busy"] = t
    t.start()
    return jsonify({"ok": True, "base": base})


def main(folder=None, port=None, open_browser=True):
    from .brand import DEFAULT_PORT, APP_NAME
    if port is None:
        port = DEFAULT_PORT
    if folder:
        STATE["folder"] = os.path.abspath(os.path.expanduser(folder))
    if open_browser:
        import webbrowser, threading as th
        th.Timer(1.0, lambda: webbrowser.open(f"http://127.0.0.1:{port}")).start()
    app.run(host="127.0.0.1", port=port, debug=False)
