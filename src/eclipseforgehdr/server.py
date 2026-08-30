"""Local web app: folder selection, pipeline run with progress, preview, export."""
from __future__ import annotations
import io, os, json, threading, traceback
import numpy as np
from flask import Flask, request, jsonify, send_file
from PIL import Image

from .pipeline import (run as run_pipeline, Progress, workdir,
                       input_fingerprint as _input_fingerprint)
from .render import Layers, render, export, DEFAULTS

app = Flask(__name__)
STATE = {"folder": None, "progress": None, "layers": None, "thread": None}


def _gui_html():
    here = os.path.dirname(__file__)
    return open(os.path.join(here, "gui.html"), encoding="utf-8").read()


@app.get("/")
def index():
    return _gui_html()


@app.get("/api/version")
def version():
    from . import __version__
    return jsonify({"version": __version__})


@app.post("/api/folder")
def set_folder():
    if STATE["thread"] and STATE["thread"].is_alive():
        return jsonify({"ok": False,
                        "error": "a run is in progress — wait for it to finish "
                                 "before switching folders"}), 400
    folder = os.path.expanduser(request.json.get("path", "").strip())
    if not os.path.isdir(folder):
        return jsonify({"ok": False, "error": f"not a folder: {folder}"}), 400
    from .raw import list_raws
    raws = list_raws(folder)
    STATE["folder"] = folder
    STATE["layers"] = None
    cached = os.path.exists(os.path.join(workdir(folder), "geometry.json"))
    return jsonify({"ok": True, "folder": folder, "raw_count": len(raws),
                    "cached": cached})


@app.post("/api/run")
def start_run():
    if STATE["folder"] is None:
        return jsonify({"ok": False, "error": "choose a folder first"}), 400
    if STATE["thread"] and STATE["thread"].is_alive():
        return jsonify({"ok": False, "error": "already running"}), 400
    force = bool(request.json.get("force", False)) if request.is_json else False
    denoise = request.json.get("denoise", "fine") if request.is_json else "fine"
    if denoise is True:
        denoise = "fine"
    if denoise is False:
        denoise = "off"
    earthshine = bool(request.json.get("earthshine", False)) if request.is_json else False
    despeckle = bool(request.json.get("despeckle", True)) if request.is_json else True
    export_tiers = bool(request.json.get("exportTiers", False)) if request.is_json else False
    tier_linear = bool(request.json.get("tierLinear", False)) if request.is_json else False
    frames = (request.json.get("frames", "all") if request.is_json else "all")
    if frames not in ("all", "best50", "best"):
        frames = "all"
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
            from . import __version__
            opts_ok = False
            if os.path.exists(opts_path):
                o = json.load(open(opts_path))
                opts_ok = (o.get("denoise") == denoise
                           and bool(o.get("earthshine", False)) == earthshine
                           and bool(o.get("despeckle", True)) == despeckle
                           and o.get("frames", "all") == frames
                           and bool(o.get("export_tiers", False)) == export_tiers
                           and bool(o.get("tier_linear", False)) == tier_linear
                           and o.get("build") == __version__
                           # ... and the input files themselves. Without this,
                           # adding or removing a frame and pressing Start (not
                           # force) silently reused the previous stack, and the
                           # exported report quoted the previous frame count.
                           and o.get("inputs") == _input_fingerprint(folder))
            have_all = all(os.path.exists(os.path.join(wd, f))
                           for f in ("prom.npy", "prom_rgb.npy", "pellett.npy"))
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
                run_pipeline(folder, prog, denoise=denoise,
                             earthshine=earthshine, despeckle=despeckle,
                             frames=frames, export_tiers=export_tiers,
                             tier_linear=tier_linear)
            else:
                prog.log("using cached pipeline products (Start with force to redo)", 0.9)
            prog.log("loading layers for preview...", None)
            STATE["layers"] = Layers(wd)
            prog.log("ready", 1.0)
            prog.done = True
        except Exception as e:
            prog.error = f"{e}\n{traceback.format_exc()}"
            prog.done = True

    t = threading.Thread(target=work, daemon=True)
    STATE["thread"] = t
    t.start()
    return jsonify({"ok": True})


@app.get("/api/progress")
def get_progress():
    p = STATE["progress"]
    if p is None:
        return jsonify({"lines": [], "frac": 0, "done": False, "error": None,
                        "ready": STATE["layers"] is not None})
    return jsonify({"lines": p.lines[-250:], "frac": p.frac, "done": p.done,
                    "error": p.error, "ready": STATE["layers"] is not None})


@app.get("/api/geometry")
def get_geometry():
    ly = STATE["layers"]
    if ly is None:
        return jsonify({"ok": False}), 400
    from . import __version__
    cy, cx, R = ly.geometry(ly.prev_decim)
    h, w = ly.prev["bg"].shape
    return jsonify({"ok": True, "W": w, "H": h, "cy": cy, "cx": cx, "R": R,
                    "defaults": DEFAULTS, "version": __version__,
                    "Rmask": ly.mask_radius(ly.prev_decim),
                    "limbProf": ([float(x) / ly.prev_decim for x in ly.limb_prof]
                                 if ly.limb_prof else None),
                    "limbMargin": ly.limb_margin / ly.prev_decim,
                    "bgChroma": [float(x) for x in ly.bg_chroma],
                    "decim": ly.prev_decim, "has_contact": ly.has_contact,
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

    threading.Thread(target=work, daemon=True).start()
    return jsonify({"ok": True})


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
    elif name == "contact":
        img = Image.fromarray((np.clip(arr, 0, 1) * 255).astype(np.uint8))
    else:
        img = Image.fromarray((np.clip(arr, 0, 1) * 255).astype(np.uint8))
    buf = io.BytesIO()
    img.save(buf, format="PNG")
    buf.seek(0)
    return send_file(buf, mimetype="image/png")


@app.post("/api/export")
def do_export():
    ly = STATE["layers"]
    if ly is None:
        return jsonify({"ok": False, "error": "run the pipeline first"}), 400
    data = request.json
    fmt = data.get("format", "tif16")
    view = data.get("view", "composite")
    size = data.get("size", "full")
    params = data.get("params", {})
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
            export(ly, params, fmt, path, view=view, size=size)
            json.dump(params, open(path + ".params.json", "w"), indent=1)
            try:
                from . import report as _report
                rp = os.path.join(workdir(STATE["folder"]), "report.json")
                st = json.load(open(rp)) if os.path.exists(rp) else {}
                st["params"] = params
                st["export"] = {"file": os.path.basename(path), "format": fmt,
                                "view": view, "size": size}
                open(os.path.splitext(path)[0] + "_report.txt", "w").write(
                    _report.build(st) + "\n")
            except Exception:
                pass
            prog.log(f"saved {path}", 1.0)
            prog.done = True
        except Exception as e:
            prog.error = str(e)
            prog.done = True

    t = threading.Thread(target=work, daemon=True)
    t.start()
    return jsonify({"ok": True, "path": path})


def main(folder=None, port=8765, open_browser=True):
    if folder:
        STATE["folder"] = os.path.abspath(os.path.expanduser(folder))
    if open_browser:
        import webbrowser, threading as th
        th.Timer(1.0, lambda: webbrowser.open(f"http://127.0.0.1:{port}")).start()
    app.run(host="127.0.0.1", port=port, debug=False)
