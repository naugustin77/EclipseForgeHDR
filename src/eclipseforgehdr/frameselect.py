"""The user's own say over which frames are used (2026-11-08).

Every frame the app finds -- the bracket, the short exposures, the contact
frames, flats, darks, bias -- can be switched OFF, and a bracket frame can be
FORCED IN past the automatic checks (not-totality saturation test, sharpness
choice "best"/"best half", the 15 % brightness test). Nothing is moved or
deleted: the choice is a small file in the image folder,
eclipseforge_frames.json, so it survives a cleared cache and is the same for
every build that opens the folder.

    {"version": 1, "off": ["frame_0012.raw", "darks/dark_0007.raw"],
                   "force": ["frame_0058.raw"]}

Keys are paths relative to the image folder, "/"-separated; a calibration
folder outside it is keyed by its absolute path.

WHERE IT ACTS. raw.list_raws() drops the switched-off frames of whatever
folder it lists (the selection file of that folder or of its parent is
consulted, which covers the image folder and its subfolders). Every reader of
frames goes through list_raws, and so do the cache fingerprints of the lights,
the short exposures, the flats, the darks and the bias: switching a frame off
changes the fingerprint and the next Start re-stacks exactly as for a removed
file. FORCED frames are added to the light fingerprint for the same reason.
"""
import hashlib
import io
import json
import os

SEL_FILE = "eclipseforge_frames.json"
_CACHE = {}
_BYPASS = [0]          # > 0: list_raws lists everything (the panel's own scan)


def sel_path(folder):
    return os.path.join(folder, SEL_FILE)


def load(folder):
    """{"off": set, "force": set} for an image folder; empty sets when there
    is no file. Cached on the file's modification time: list_raws runs often."""
    p = sel_path(folder)
    try:
        mt = os.stat(p).st_mtime_ns
    except OSError:
        return {"off": set(), "force": set(), "dirs": {}}
    c = _CACHE.get(p)
    if c and c[0] == mt:
        return c[1]
    try:
        j = json.load(open(p))
        v = {"off": set(map(str, j.get("off") or [])), "force": set(map(str, j.get("force") or [])),
             "dirs": {str(k): str(d) for k, d in (j.get("dirs") or {}).items() if d}}
    except Exception:
        v = {"off": set(), "force": set(), "dirs": {}}
    _CACHE[p] = (mt, v)
    return v


def save(folder, off, force, dirs=None):
    p = sel_path(folder)
    off, force = sorted(set(off)), sorted(set(force) - set(off))
    if dirs is None:
        dirs = load(folder).get("dirs") or {}
    dirs = {k: v for k, v in dirs.items() if v}
    if not off and not force and not dirs:
        try:
            os.remove(p)
        except OSError:
            pass
    else:
        tmp = p + ".tmp"
        doc = {"version": 1, "off": off, "force": force}
        if dirs:
            doc["dirs"] = dict(sorted(dirs.items()))
        json.dump(doc, open(tmp, "w"), indent=1)
        os.replace(tmp, p)
    _CACHE.pop(p, None)


# FOLDERS CHOSEN BY HAND (2026-11-08): darks, bias, short exposures and
# contact frames are found by their folder names; when that fails, or finds the
# wrong one, the user names the folder in tab 1. Stored in the same file under
# "dirs": {"darks": "path", "short": "off", ...}; a path inside the image
# folder is kept relative. "off": use none even if one would be found.
DIR_ROLES = ("darks", "bias", "short", "contact", "phases")


def dir_override(folder, role):
    """("auto", None) | ("off", None) | ("path", absolute path)."""
    try:
        v = (load(folder).get("dirs") or {}).get(role)
    except Exception:
        v = None
    if not v:
        return "auto", None
    if v == "off":
        return "off", None
    p = v if os.path.isabs(v) else os.path.join(folder, v)
    return ("path", p) if os.path.isdir(p) else ("auto", None)


def set_dir(folder, role, value):
    """value: "" (auto), "off", or a folder path."""
    sel = load(folder)
    dirs = dict(sel.get("dirs") or {})
    v = (value or "").strip()
    if not v or v.lower() == "auto":
        dirs.pop(role, None)
    elif v.lower() in ("off", "none", "-"):
        dirs[role] = "off"
    else:
        dirs[role] = key(os.path.expanduser(v), folder)
    save(folder, sel["off"], sel["force"], dirs)


def key(path, folder):
    a, b = os.path.abspath(path), os.path.abspath(folder)
    try:
        r = os.path.relpath(a, b)
    except ValueError:              # another drive (Windows)
        return a.replace(os.sep, "/")
    if r.startswith(".."):
        return a.replace(os.sep, "/")
    return r.replace(os.sep, "/")


def drop_off(paths, listing_dir):
    """`paths` without the frames switched off, judged by the selection file
    of the listed folder and of its parent."""
    if not paths or _BYPASS[0]:
        return paths
    d = os.path.abspath(listing_dir)
    for base in (d, os.path.dirname(d)):
        off = load(base)["off"]
        if off:
            paths = [p for p in paths if key(p, base) not in off]
    return paths


def forced(folder):
    """Keys of the frames forced in (bracket frames, relative to `folder`)."""
    return set(load(folder)["force"])


def fingerprint_extra(folder):
    """What the forced set adds to the light fingerprint; [] when nothing is forced."""
    f = sorted(forced(folder))
    return [["#forced", f]] if f else []


# ------------------------------------------------------------------ the groups

ROLES = (("lights", "Corona brackets", True), ("short", "Short exposures (prominences)", False),
         ("contact", "Contact frames", False), ("phases", "Partial phases (Extras)", False),
         ("flats", "Flats", False), ("darks", "Darks", False), ("bias", "Bias", False))


def groups(folder, flat_dir=None, bias_dir=None, dark_dir=None):
    """Every frame the app finds for this folder, by role, switched-off frames
    included: [(role, title, forceable, dir, [paths])]."""
    from .raw import list_raws
    from .pipeline import find_short_dir, resolve_flat_dir, resolve_calib_dir
    from .contacts import find_contact_dir
    dirs = {"lights": folder}
    # the folders are found with nothing switched off: a short-exposure or
    # darks folder whose frames are all off must still be listed here, or it
    # could never be switched back on
    _BYPASS[0] += 1
    try:
        _found = _find_dirs(folder, flat_dir, bias_dir, dark_dir, find_short_dir,
                            find_contact_dir, resolve_flat_dir, resolve_calib_dir)
    finally:
        _BYPASS[0] -= 1
    dirs.update(_found)
    out = []
    for role, title, forceable in ROLES:
        d = dirs.get(role)
        if not d or not os.path.isdir(d):
            continue
        try:
            ps = list_raws(d, selected=False)
        except Exception:
            ps = []
        if ps:
            out.append((role, title, forceable, d, ps))
    return out


def _find_dirs(folder, flat_dir, bias_dir, dark_dir, find_short_dir, find_contact_dir,
               resolve_flat_dir, resolve_calib_dir):
    dirs = {}
    from .frameprod import find_phases_dir
    for role, fn in (("short", lambda: find_short_dir(folder)),
                     ("contact", lambda: find_contact_dir(folder)),
                     ("phases", lambda: find_phases_dir(folder)),
                     ("flats", lambda: resolve_flat_dir(folder, flat_dir)),
                     ("darks", lambda: resolve_calib_dir(folder, "dark", dark_dir)),
                     ("bias", lambda: resolve_calib_dir(folder, "bias", bias_dir))):
        try:
            dirs[role] = fn()
        except Exception:
            dirs[role] = None
    return dirs


# ------------------------------------------------------------------ previews

def _embedded_jpeg(path):
    """The camera's own preview inside a raw file, as JPEG bytes, or None."""
    try:
        import rawpy
        with rawpy.imread(path) as r:
            t = r.extract_thumb()
            if t.format == rawpy.ThumbFormat.JPEG:
                return bytes(t.data)
    except Exception:
        pass
    # without rawpy (or a format it will not hand out): the largest JPEG
    # stream in the file. Lossless-JPEG raw data does not open in PIL and is
    # passed over.
    from PIL import Image
    try:
        b = open(path, "rb").read()
    except OSError:
        return None
    best, i = None, 0
    while True:
        i = b.find(b"\xff\xd8\xff", i)
        if i < 0:
            break
        try:
            im = Image.open(io.BytesIO(b[i:]))
            if im.format == "JPEG":
                n = im.size[0] * im.size[1]
                if best is None or n > best[1]:
                    best = (i, n)
        except Exception:
            pass
        i += 3
    return b[best[0]:] if best else None


STRETCHES = ("strong", "gentle", "exposed", "soft")


def _stretch(a, mode="strong", white=None):
    """A linear frame to 8-bit for a look, one curve for all channels.

    strong:  asinh (k 25) from the sky's level (0.5th percentile of the
             luminance) to the frame's brightest (99.95th): faint corona
             shows, only what is truly saturated prints white.
    gentle:  the same range, asinh k 4: closer to how the frame looks.
    exposed: no normalising at all: the data against the sensor's full scale
             (`white`), with the sRGB curve: a short exposure looks dark, a
             dark frame almost black, as they were shot.
    soft:    as shot like `exposed`, but with an asinh curve (k 6) from the
             frame's own dark level: the highlights roll off instead of
             printing flat white -- for the film of the frames (Extras)."""
    import numpy as np
    a = np.asarray(a, np.float32)
    if a.ndim == 3 and a.shape[0] in (3, 4) and a.shape[-1] not in (3, 4):
        a = np.moveaxis(a, 0, -1)
    if a.ndim == 3:
        a = a[..., :3]
    if mode in ("exposed", "soft"):
        if not white:
            mx = float(np.max(a)) if a.size else 1.0
            white = 1.0 if mx <= 1.0 else (255.0 if mx <= 255 else (4095.0 if mx <= 4095 else
                                                                      (16383.0 if mx <= 16383 else 65535.0)))
        if mode == "soft":
            L = a.mean(-1) if a.ndim == 3 else a
            lo = float(np.percentile(L[::4, ::4], 0.3))
            x = np.clip((a - lo) / max(white - lo, 1e-9), 0, 1)
            x = np.arcsinh(6.0 * x) / np.arcsinh(6.0)
            return (x * 255 + 0.5).astype(np.uint8)
        x = np.clip(a / white, 0, 1)
        x = np.where(x <= 0.0031308, 12.92 * x, 1.055 * np.power(x, 1 / 2.4) - 0.055)
        return (x * 255 + 0.5).astype(np.uint8)
    L = a.mean(-1) if a.ndim == 3 else a
    s = L[::4, ::4]
    lo, hi = float(np.percentile(s, 0.5)), float(np.percentile(s, 99.95))
    x = np.clip((a - lo) / max(hi - lo, 1e-9), 0, 1)
    k = 4.0 if mode == "gentle" else 25.0
    x = np.arcsinh(k * x) / np.arcsinh(k)
    return (x * 255 + 0.5).astype(np.uint8)


import threading
_DECODE = threading.Semaphore(2)     # raw decodes at a time (memory, cores)


def _from_raw(path, full=False, mode="strong"):
    """The raw data itself, linear, camera white balance, stretched: unlike the
    camera's 8-bit JPEG it shows the whole range the sensor caught. Half size
    unless `full` (the 1:1 window)."""
    import rawpy
    with _DECODE:
        with rawpy.imread(path) as r:
            # user_flip=0: as the sensor recorded it, the way the stack reads
            # the raw (raw_image_visible). Left to the file, the decoder turns
            # a frame the camera's orientation sensor flagged as portrait by
            # 90 deg, and the frame no longer matches the stack
            rgb = r.postprocess(half_size=not full, output_bps=16, gamma=(1, 1),
                                no_auto_bright=True, use_camera_wb=True, user_flip=0)
    return _stretch(rgb, mode, white=65535.0)


def _pil(path, full=False, mode="strong"):
    from PIL import Image
    from .fits import is_fits
    ext = os.path.splitext(path)[1].lower()
    if is_fits(path):
        from .fits import read_fits
        d, _ = read_fits(path)
        return Image.fromarray(_stretch(d, mode))
    if ext in (".tif", ".tiff"):
        import tifffile
        return Image.fromarray(_stretch(tifffile.imread(path), mode))
    try:
        return Image.fromarray(_from_raw(path, full=full, mode=mode))
    except Exception:
        pass                       # no rawpy, or a file it will not read
    j = _embedded_jpeg(path)
    if j is None:
        raise ValueError("no preview in this file")
    im = Image.open(io.BytesIO(j))
    if not full:
        im.draft("RGB", (1600, 1600))
    return im.convert("RGB")


def preview(path, kind, cache_dir, mode="strong"):
    """JPEG file for a frame: kind "thumb" (240 px), "big" (1600 px) or "crop"
    (a 1:1 window of 900 px on the disc). Cached by path, size and date."""
    st = os.stat(path)
    h = hashlib.sha1(f"{os.path.abspath(path)}|{st.st_size}|{st.st_mtime_ns}|{kind}|{mode}|4".encode()).hexdigest()[:20]
    os.makedirs(cache_dir, exist_ok=True)
    out = os.path.join(cache_dir, h + ".jpg")
    if os.path.exists(out):
        return out
    if mode not in STRETCHES:
        mode = "strong"
    im = _pil(path, full=(kind == "crop"), mode=mode)
    if kind == "crop":
        import numpy as np
        from scipy import ndimage
        g = np.asarray(im.convert("L").reduce(8), np.float32)
        g = ndimage.uniform_filter(g, 9)
        # centred on the brightest smoothed area: the inner corona at the
        # limb, where sharpness, cloud and trailing show first
        yy, xx = np.unravel_index(int(np.argmax(g)), g.shape)
        cy, cx = yy * 8 + 4, xx * 8 + 4
        W, H = im.size
        s = 900
        x0 = int(min(max(cx - s // 2, 0), max(W - s, 0)))
        y0 = int(min(max(cy - s // 2, 0), max(H - s, 0)))
        im = im.crop((x0, y0, x0 + min(s, W), y0 + min(s, H)))
    else:
        n = 240 if kind == "thumb" else 1600
        im.thumbnail((n, n))
    im.save(out, quality=92)
    return out
