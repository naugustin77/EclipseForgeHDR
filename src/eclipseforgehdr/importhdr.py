"""Import a finished HDR and run the enhancement layers on it.

For people who already have a merged corona image -- from Siril, PixInsight,
Photoshop, AstroSurface, or from this app's own aligned tier TIFFs -- and want
MGN, FNRGF, NAFE-VN, the inner-corona layer, Pellett and the prominence gate
without going through stacking again. It also sidesteps every stage that can
fail on an awkward bracket, which is why it exists.

WHAT IT DOES NOT DO, and the report says so:

  * No photometry, no alignment, no per-tier lunar masking. Those describe a
    stack, and there is no stack here.
  * The inner-corona layer stops being an INDEPENDENT source. Normally it is a
    separate MGN of the four shortest tiers, which see the inner corona
    unsaturated; from one image it is a second view of the same pixels at
    different scales. Still useful, no longer independent.
  * The prominence gate is measurably weaker. It keys on H-alpha redness in a
    fast tier where the chromosphere is not blown; a merged image has usually
    tone-compressed exactly that. Measured on the reference stack: the corona's
    own R/GB reads 1.84 against 3.02 from the real fast tier.
  * No earthshine layer -- that needs the longest tiers on their own.

LINEARITY is the one thing that has to be right. MGN, FNRGF and NAFE all work
on log luminance and assume the value is proportional to coronal brightness. A
display-gamma image has a corona whose falloff is roughly half as steep as the
truth -- measured on a Photoshop stack of this app's own sRGB tier exports, the
log-log radial slope read -1.72 against -3.38 for the same scene linear, a
ratio of 0.510 where sRGB predicts ~0.45.

So the transfer function is not guessed. It is read from the file's embedded
ICC profile, which is where the answer actually is, and only falls back to a
declaration when there is no profile. Applying the inverse to that same stack
brought the slope back to -3.17, within 6% of the truth.
"""
from __future__ import annotations
import os, json, struct
import numpy as np


class ImportError_(RuntimeError):
    pass


# ---------- transfer function ----------

def _icc_gamma(icc):
    """Read the tone response out of an ICC profile.

    Returns ("linear", None), ("gamma", g), ("srgb", None) or (None, None) when
    the profile does not say. Only the red TRC is inspected: a profile whose
    channels disagree is not one we can use anyway.
    """
    try:
        n = struct.unpack(">I", icc[128:132])[0]
        tag = None
        for i in range(n):
            sig, off, size = struct.unpack(">4sII", icc[132 + 12 * i:144 + 12 * i])
            if sig == b"rTRC":
                tag = icc[off:off + size]
                break
        if tag is None or tag[:4] != b"curv":
            return None, None
        cnt = struct.unpack(">I", tag[8:12])[0]
        if cnt == 0:
            return "linear", 1.0                      # identity curve
        if cnt == 1:                                  # u8Fixed8 gamma
            g = struct.unpack(">H", tag[12:14])[0] / 256.0
            return ("linear", 1.0) if abs(g - 1.0) < 0.02 else ("gamma", g)
        # a sampled curve: identify it by shape rather than by name
        t = np.frombuffer(tag[12:12 + 2 * cnt], dtype=">u2").astype(np.float64) / 65535.0
        x = np.linspace(0.0, 1.0, cnt)
        if float(np.max(np.abs(t - x))) < 0.01:
            return "linear", 1.0
        srgb = np.where(x <= 0.04045, x / 12.92, ((x + 0.055) / 1.055) ** 2.4)
        if float(np.max(np.abs(t - srgb))) < 0.02:
            return "srgb", None
        # otherwise fit a pure power law to it
        m = (x > 0.05) & (x < 0.95) & (t > 1e-6)
        if m.sum() > 8:
            g = float(np.polyfit(np.log(x[m]), np.log(t[m]), 1)[0])
            if 0.5 < g < 4.0:
                return "gamma", g
    except Exception:
        pass
    return None, None


def _linearise(a, kind, g, progress):
    log = (progress.log if progress is not None else (lambda *x, **k: None))
    if kind == "linear":
        log("import: the file declares a linear tone curve — used as it is", None)
        return a
    if kind == "srgb":
        log("import: the file declares the sRGB transfer function — inverting it", None)
        return np.where(a <= 0.04045, a / 12.92, ((a + 0.055) / 1.055) ** 2.4)
    if kind == "gamma":
        log(f"import: the file declares gamma {g:.3f} — inverting it", None)
        return np.power(np.maximum(a, 0.0), g)
    raise ImportError_(
        "this image carries no colour profile, so whether it is scene-linear "
        "cannot be read from the file. Re-export it with a profile embedded, or "
        "put 'linear' in the filename to declare that it already is.")


# ---------- reading ----------

def read_image(path, progress=None, assume=None):
    """Read one HDR image and return it scene-linear, float32 HxWx3.

    `assume` overrides the file: "linear", "srgb", or a numeric gamma.
    """
    log = (progress.log if progress is not None else (lambda *x, **k: None))
    ext = os.path.splitext(path)[1].lower()
    icc = None
    if ext in (".fit", ".fits", ".fts"):
        from .fits import read_fits
        a, _hdr = read_fits(path)
        a = np.asarray(a, np.float32)
        kind, g = "linear", 1.0            # FITS is linear by construction
        log("import: FITS — linear by construction", None)
    else:
        import tifffile
        with tifffile.TiffFile(path) as tf:
            pg = tf.pages[0]
            bits = pg.tags["BitsPerSample"].value
            bits = bits[0] if isinstance(bits, (tuple, list)) else bits
            if "InterColorProfile" in pg.tags:
                icc = pg.tags["InterColorProfile"].value
            a = pg.asarray()
        if bits < 16 and a.dtype == np.uint8:
            raise ImportError_(
                "this is an 8-bit image. The corona spans several thousand to "
                "one; 8 bits cannot hold it and the detail layers would be "
                "working on quantisation steps. Re-export at 16 bits or more.")
        a = np.asarray(a, np.float32)
        if a.dtype != np.float32 or a.max() > 1.001:
            a = a / (65535.0 if a.max() > 255.5 else 255.0)
        kind, g = (None, None)
        if icc is not None:
            kind, g = _icc_gamma(icc)
            if kind:
                log("import: colour profile found in the file", None)
    if assume:
        try:
            g2 = float(assume)
            kind, g = "gamma", g2
        except (TypeError, ValueError):
            kind, g = str(assume).lower(), None
        log(f"import: tone curve overridden by the caller ({assume})", None)
    elif kind is None and ("linear" in os.path.basename(path).lower()
                           or os.environ.get("ECLIPSEFORGE_TIFF_LINEAR") == "1"):
        kind, g = "linear", 1.0
    a = np.clip(np.nan_to_num(np.asarray(a, np.float32)), 0.0, None)
    if a.ndim == 2:
        a = np.stack([a] * 3, axis=-1)
    if a.ndim != 3 or a.shape[2] < 3:
        raise ImportError_(f"expected a 2-D or 3-channel image, got {a.shape}")
    a = a[:, :, :3]
    a = _linearise(a, kind, g, progress)
    # scale to a comfortable working range; every layer downstream is scale-free
    hi = float(np.percentile(a, 99.9))
    if hi <= 0:
        raise ImportError_("this image has no signal")
    return (a * (50000.0 / hi)).astype(np.float32), (kind, g)


# ---------- the run ----------

def run(folder, image_path, progress, denoise="fine", assume=None):
    """Build every cached product the renderer needs, from one image."""
    from .pipeline import (workdir, find_disc, fit_limb_rays, remove_sky_gradient,
                           _exp_name)
    from . import detail, report as _report
    from . import __version__
    import datetime

    wd = workdir(folder)
    # No stacking on this path: reading the file, finding the disc and fitting
    # the sky take seconds, and build_layers takes the rest of the run. Leaving
    # the bar's default split would have shown 78% for that first handful of
    # seconds and then crawled -- exactly backwards.
    try:
        progress.bar_detail = 0.12
    except Exception:
        pass
    progress.log(f"importing {os.path.basename(image_path)} ...", 0.05)
    rgb, (kind, g) = read_image(image_path, progress, assume=assume)
    H, W, _ = rgb.shape
    progress.log(f"{W}x{H}, scene-linear, "
                 f"{np.log2(float(np.percentile(rgb, 99.99)) / max(float(np.percentile(rgb, 1)), 1e-6)):.1f} EV "
                 f"between the 1st and 99.99th percentile", 0.15)

    lum = (0.2126 * rgb[:, :, 0] + 0.7152 * rgb[:, :, 1]
           + 0.0722 * rgb[:, :, 2]).astype(np.float32)
    progress.log("locating the lunar disc ...", 0.2)
    d = find_disc(lum)
    fit = fit_limb_rays(lum, *d) if d is not None else None
    if fit is None:
        raise ImportError_(
            "could not find the lunar limb in this image. It has to be a "
            "totality frame with the Moon fully inside the field.")
    cy, cx, R, rms, nk, nt, prof = fit
    progress.log(f"lunar limb: centre ({cy:.1f},{cx:.1f}) R={R:.1f}px, "
                 f"rms {rms:.2f}px over {nk}/{nt} rays", 0.25)

    stats = {"version": __version__, "folder": folder, "imported": image_path,
             "import_tone": kind, "import_gamma": g, "n_files": 1,
             "W": W, "H": H, "mode": "imported HDR",
             "options": {"denoise": denoise},
             "geometry": {"cy": cy, "cx": cx, "R": R, "rms": rms,
                          "rays_kept": nk, "rays": nt}}

    np.save(os.path.join(wd, "hdr_rgb.npy"), rgb)
    np.save(os.path.join(wd, "hdr_lum.npy"), lum)
    # No separate short-exposure stack exists, so the inner layer runs on the
    # same pixels. It is still a different filter at different scales, but it
    # is no longer an independent measurement and the report says so.
    np.save(os.path.join(wd, "short_lum.npy"), lum)
    h2, w2 = H // 2, W // 2
    np.save(os.path.join(wd, "prom_rgb.npy"),
            rgb[:2 * h2, :2 * w2].reshape(h2, 2, w2, 2, 3)
            .mean(axis=(1, 3)).astype(np.float32))
    for gone in ("long_lum.npy",):
        p = os.path.join(wd, gone)
        if os.path.exists(p):
            os.remove(p)

    margin = max(4.0, 0.042 * R)
    json.dump({"cy": cy, "cx": cx, "R": R, "Rmask": R + margin,
               "limb_margin": margin,
               "limb_prof": [float(x) for x in prof] if prof is not None else None,
               "inner_geom": {"cy": cy, "cx": cx, "R": R},
               "prom_geom": ({"cy": cy / 2, "cx": cx / 2, "R": R / 2,
                              "prof": [float(x) / 2 for x in prof]}
                             if prof is not None else None),
               "secs": [], "cal": {}, "abs_shift": {}},
              open(os.path.join(wd, "geometry.json"), "w"), indent=1)

    # The sky fit lives in the stacking path, but it works on the merged image
    # alone -- so an imported HDR gets it too. Without it the sky's colour
    # gradient goes straight into the composite and into NAFE.
    progress.log("fitting the sky gradient ...", 0.35)
    try:
        st = _report.measure_image(lum, cy, cx, R)
        stats.update(st)
        remove_sky_gradient(wd, cy, cx, R, st.get("corona_extent_R"), stats, progress)
    except Exception as e:
        progress.log(f"sky gradient not removed ({e})", None)

    progress.log("building the detail layers ...", 0.45)
    lstats = detail.build_layers(wd, progress, denoise=denoise, earthshine=False)
    if isinstance(lstats, dict):
        stats.update(lstats)
    try:
        stats.update(_report.measure_image(np.load(os.path.join(wd, "hdr_lum.npy")),
                                           cy, cx, R))
    except Exception:
        pass
    stats["finished"] = datetime.datetime.now().isoformat(timespec="seconds")
    json.dump({"import": os.path.basename(image_path),
               "import_mtime": int(os.path.getmtime(image_path)),
               "import_size": int(os.path.getsize(image_path)),
               "denoise": denoise, "assume": assume,
               "mode": "import", "build": __version__},
              open(os.path.join(wd, "opts.json"), "w"))
    txt = _report.write(wd, stats)
    for line in txt.split("\n"):
        progress.log(line, None)
    progress.log("import complete", 1.0)
    return stats
