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
    return _gui_html()


@app.get("/api/version")
def version():
    from . import __version__
    from .brand import APP_NAME, LAB
    return jsonify({"version": __version__, "app": APP_NAME, "lab": LAB})


# ---------- file browser ----------
#
# Pasting an absolute path is fine once you know what the app wants and awkward
# every other time -- both testers lost time to it, one of them pasting a
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
                    "tiffs": cand[:24]})


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
    return have < LAYER_BUILD, have, LAYER_BUILD


def _rebuild_layers(wd, prog, denoise, fnrgf_preset, partialconv, earthshine=False,
                    force=False, why=""):
    """Rebuild every detail layer from the cached merge and record the recipe
    in report.json. Mirrors the tail of the pipeline: build_layers is the last
    stage there, so nothing after it is skipped. `force`: rebuild even when the
    recipe is current, because Denoise changed -- it acts on the layers only,
    so it never needs a re-stack (2026-10-08)."""
    from . import detail
    stale, have, want = _layers_stale(wd)
    if not stale and not force:
        return False
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
    if _ip and os.path.isfile(_ip):
        STATE["folder"] = os.path.dirname(os.path.abspath(_ip))
    elif STATE["folder"] is None:
        return jsonify({"ok": False,
                        "error": "choose a folder first, or paste the path "
                                 "of one finished HDR to import"}), 400
    _b = _busy()
    if _b:
        return jsonify({"ok": False, "error": f"{_b} is still running"}), 400
    force = bool(request.json.get("force", False)) if request.is_json else False
    # DEFAULTS (2026-10-07, Nico, tested on his 600 mm and Clifton's 250 mm sets):
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
    if tier_mode not in ("exposure", "frame"):
        tier_mode = "exposure"
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
    # Default "auto" since 0.23.9 (Nico, 21 Sep): it measures on every run and
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
    # `vng_demosaic.txt` in the folder of raws, the same way the Hill
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
                                partialconv=partialconv)
                    try:
                        _o2 = json.load(open(opts_path))
                        _o2["flat_inputs"] = _flat_fp(_fd)
                        json.dump(_o2, open(opts_path, "w"), indent=1)
                    except Exception:
                        pass
                else:
                    prog.log("using cached simple stack", 0.9)
                    _rebuild_layers(wd, prog, denoise, fnrgf_preset, partialconv)
                # THE MOON LAYER ON THE SIMPLE STACK TOO (2026-10-08, Val's
                # FITS): it reads the long frames itself and fits its own
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
            if os.path.exists(opts_path):
                o = json.load(open(opts_path))
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
                           and o.get("fnrgf_preset", "ours") == fnrgf_preset
                           and o.get("stack_combine", "mean") == stack_combine
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
                           and _cache_ok(o.get("build"))
                           # ... and the input files themselves. Without this,
                           # adding or removing a frame and pressing Start (not
                           # force) silently reused the previous stack, and the
                           # exported report quoted the previous frame count.
                           and o.get("inputs") == _input_fingerprint(folder))
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
                             "align_corr": align_corr, "intra_lock": intra_lock,
                             "colour_planes": colour_planes, "tier_mode": tier_mode}
                    _dflt = {"wb_source": "daylight", "demosaic": "mhc", "photometry": "scalar",
                             "photo_solve": "chain", "fnrgf_preset": "ours", "stack_combine": "mean",
                             "align_filter": "isotropic", "align_corr": "semi", "intra_lock": "moon",
                             "colour_planes": "off", "tier_mode": "exposure", "feather": "plain",
                             "frames": "all"}
                    _diff = []
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
                        _diff.append("raw files changed (name, size or date)")
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
                             align_filter=align_filter,
                             align_corr=align_corr,
                             intra_lock=intra_lock,
                             colour_planes=colour_planes,
                             partialconv=partialconv,
                             tier_mode=tier_mode)
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
                _rebuild_layers(wd, prog, denoise, fnrgf_preset, partialconv,
                                earthshine=earthshine, force=_dn_changed,
                                why=f"denoise {_oc.get('denoise')!r} -> {denoise!r}")
                if _dn_changed:
                    try:
                        _oc["denoise"] = denoise
                        json.dump(_oc, open(opts_path, "w"), indent=1)
                    except Exception:
                        pass
            # HILL MASKS FROM A CACHED STACK. Everything they need -- the merged
            # luminance, the geometry, the prominence mask -- is already on
            # disk, so a folder stacked before 0.22.64 gains the layer for the
            # cost of the polar blur alone instead of a whole re-stack.
            # SWITCHED OFF: SKIP THIS WHOLE BLOCK, not just the staleness
            # test. Setting _hstale = False here was not enough -- the next
            # branch is `if not _hstale`, which re-derives staleness from
            # report.json, finds no "hill" entry because the run did not write
            # one, reads its build number as 1, and rebuilds. That is exactly
            # what happened on Nico's 600 mm run: the pipeline reported
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
    layer set or the Hill masks are behind the build."""
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
        except Exception as e:
            out["error"] = str(e)
    if out["cached"]:
        try:
            out["layers_stale"] = bool(_layers_stale(wd)[0])
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
    return jsonify({"ok": True, "W": w, "H": h, "cy": cy, "cx": cx, "R": R,
                    "defaults": defaults_for(getattr(ly, "mode", None)),
                    "version": __version__,
                    "app": _BRAND_NAME, "lab": _BRAND_LAB,
                    "Rmask": ly.mask_radius(ly.prev_decim),
                    "limbProf": ([float(x) / ly.prev_decim for x in ly.limb_prof]
                                 if ly.limb_prof else None),
                    "limbMargin": ly.limb_margin / ly.prev_decim,
                    "bgChroma": [float(x) for x in ly.bg_chroma],
                    "rprofile": [round(float(x), 5) for x in ly.rprof_1d],
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
                    "has_nafe": getattr(ly, "has_nafe", True),
                    "has_mgn_fine": getattr(ly, "has_mgn_fine", False),
                    "flat_range": list(getattr(ly, "flat_range", []) or []),
                    "flat_error": getattr(ly, "flat_error", None),
                    "has_earth": getattr(ly, "has_earth", False)})


@app.post("/api/contact")
def load_contact():
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
    """Half-range of the transport encoding for Hill mask i, in its own units.
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
    "fnrgf_preset":  ("ours",      "fnrgfPreset"),
    "photo_solve":   ("chain",     "photoSolve"),
    "tier_mode":     ("exposure",  "tierMode"),
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
                    "n": len(proc), "file": os.path.basename(path)})


def _older_than(a, b):
    """Version compare on dotted integers; unparseable sorts as 'older'."""
    def t(v):
        try:
            return tuple(int(x) for x in str(v).split("."))
        except ValueError:
            return ()
    ta, tb = t(a), t(b)
    return (ta < tb) if (ta and tb) else True


@app.post("/api/export")
def do_export():
    ly = STATE["layers"]
    if ly is None:
        return jsonify({"ok": False, "error": "run the pipeline first"}), 400
    data = request.json
    fmt = data.get("format", "tif16")
    view = data.get("view", "composite")
    size = data.get("size", "full")
    prom_layer = bool(data.get("promLayer", False))
    moon_layer = bool(data.get("moonLayer", False))
    params = data.get("params", {})
    # Purely cosmetic and applied last, so it needs no re-run and no re-render
    # of anything cached -- see render.apply_orient.
    _or = str(data.get("orient", "") or "").strip().lower()
    params = dict(params, orient=(_or if _or in ("flipv", "fliph", "180", "cw", "ccw") else ""))
    ext = {"tif16": ".tif", "tif8": ".tif", "png": ".png", "jpg": ".jpg"}[fmt]
    name = data.get("name") or "eclipseforge_render"
    if view != "composite":
        name += "_" + view
    if size == "half":
        name += "_50pct"
    outdir = os.path.join(STATE["folder"], "eclipseforge_output")
    os.makedirs(outdir, exist_ok=True)
    base = os.path.join(outdir, name)
    path = base + ("_16bit" if fmt == "tif16" else "") + ext
    i = 1
    while os.path.exists(path):
        path = base + ("_16bit" if fmt == "tif16" else "") + f"_{i}" + ext
        i += 1
    prog = Progress()
    STATE["progress"] = prog
    prog.log(f"exporting {fmt} ...", 0.1)

    def work():
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
                    _notes.append(("[warn] " if "BAD" in _sm else "[odd] " if "check" in _sm else "[ok] ")
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
                        export(ly, dict(params, promLayer=0.0), fmt, _cb, view=view, size=size,
                               keep=_kc)
                        prog.log(f"[ok] saved {_cb} (corona without prominence layer)", 0.85)
                        _pl.export_layer2(ly, params, _pp, _keep.get("rgb"), _kc.get("rgb"),
                                          size=size)
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
    return jsonify({"ok": True, "path": path})


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
