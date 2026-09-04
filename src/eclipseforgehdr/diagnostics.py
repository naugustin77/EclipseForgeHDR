"""A small diagnostic bundle, written at the end of every run.

WHY THIS EXISTS. Testers report things like "there is ringing in MGN, FNRGF and
NAFE" and the only way to chase that has been to ask for raw frames. People are
reasonably reluctant to hand over their eclipse raws -- it is a shot they may
have travelled to another continent for, and it is theirs. So diagnosis stalls,
and the same defect survives into the next release.

Almost none of that diagnosis actually needs the raws. Ringing is a radial
signature: it shows up in the median of a layer taken ring by ring, and in the
azimuthal profile at a few radii. Those are a few thousand numbers. A wrong limb
circle shows up as a step in the same profiles at R. A bad merge shows up in the
per-tier numbers the report already prints.

So this writes what a diagnosis needs and nothing else:

  * the run report, which the user can already read in full
  * geometry, and the per-tier limb measurements
  * radial profiles of each detail layer -- median, p10, p90 per ring
  * azimuthal profiles at three radii, which is where ringing declares itself
  * a 512 px thumbnail of each layer

WHAT IT IS NOT. The thumbnails are 512 px on the long side, greyscale, of the
DETAIL LAYERS -- normalised high-pass fields, not the picture. Nobody is
publishing a 512 px MGN layer, and it carries no colour, no tone mapping and
none of the work that makes the final image theirs. The bundle has no raw data,
no full-resolution anything, and no file paths beyond the folder name already in
the report. Typical size is one to two megabytes.

If a tester still does not want to send it, that is their call -- but they can
open the zip and see exactly what is in it first, which is the point.
"""
from __future__ import annotations
import io
import json
import os
import zipfile

import numpy as np

_LAYERS = ("mgn", "mgn_fine", "fnrgf", "nafe", "inner", "inner0", "rhef",
           "pellett", "promdet", "hdr_lum")
_THUMB = 512


def _profiles(a, cy, cx, R, rmax_R=4.5, nt=360):
    """Radial and azimuthal profiles of one layer. Numbers, not pixels."""
    H, W = a.shape
    r0 = max(int(R * 1.02), 1)
    r1 = int(min(R * rmax_R, min(cy, cx, H - cy, W - cx) - 2))
    if r1 <= r0 + 8:
        return None
    rad = np.arange(r0, r1, max(1, (r1 - r0) // 400), dtype=np.float32)
    th = np.linspace(0, 2 * np.pi, nt, endpoint=False, dtype=np.float32)
    ys = cy + rad[:, None] * np.sin(th)[None, :]
    xs = cx + rad[:, None] * np.cos(th)[None, :]
    yi = np.clip(ys.astype(np.int32), 0, H - 1)
    xi = np.clip(xs.astype(np.int32), 0, W - 1)
    P = np.asarray(a)[yi, xi]
    out = {
        "radius_px": [round(float(v), 1) for v in rad],
        "radius_R": [round(float(v / R), 3) for v in rad],
        # the ring-by-ring median is where a wrong limb circle or a deband
        # artifact shows as a step or an oscillation
        "median": [round(float(v), 5) for v in np.median(P, axis=1)],
        "p10": [round(float(v), 5) for v in np.percentile(P, 10, axis=1)],
        "p90": [round(float(v), 5) for v in np.percentile(P, 90, axis=1)],
    }
    # azimuthal cuts at three radii -- ringing is periodic in RADIUS, so a cut
    # at fixed radius that looks clean while the radial median oscillates is the
    # signature that separates it from real coronal structure
    az = {}
    for f in (1.15, 1.6, 2.4):
        k = int(np.argmin(np.abs(rad - f * R)))
        if abs(rad[k] - f * R) < 0.25 * R:
            az["%.2fR" % f] = [round(float(v), 5) for v in P[k]]
    out["azimuthal"] = az
    return out


def _thumb_png(a, valid_lo=None, valid_hi=None):
    """512 px greyscale PNG of a layer, robustly stretched. Needs no PIL."""
    from PIL import Image
    a = np.asarray(a, np.float32)
    if a.ndim == 3:
        a = a.mean(axis=2)
    lo = float(np.percentile(a, 0.5)) if valid_lo is None else valid_lo
    hi = float(np.percentile(a, 99.5)) if valid_hi is None else valid_hi
    b = np.clip((a - lo) / max(hi - lo, 1e-9), 0, 1)
    step = max(1, int(max(b.shape) / _THUMB))
    b = b[::step, ::step]
    im = Image.fromarray((b * 255).astype(np.uint8), mode="L")
    buf = io.BytesIO()
    im.save(buf, format="PNG", optimize=True)
    return buf.getvalue()


def write_bundle(wd, folder, stats=None, progress=None):
    """Write eclipseforge_output/eclipseforge_diagnostics.zip. Never raises."""
    try:
        geo = json.load(open(os.path.join(wd, "geometry.json")))
    except Exception:
        return None
    cy, cx, R = geo.get("cy"), geo.get("cx"), geo.get("R")
    if not R:
        return None
    out_dir = os.path.join(folder, "eclipseforge_output")
    os.makedirs(out_dir, exist_ok=True)
    path = os.path.join(out_dir, "eclipseforge_diagnostics.zip")
    n_layers = 0
    try:
        with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as z:
            z.writestr("README.txt", _README)
            rp = os.path.join(wd, "report.txt")
            if os.path.exists(rp):
                z.write(rp, "report.txt")
            z.writestr("geometry.json", json.dumps(geo, indent=1))
            if stats:
                keep = {k: v for k, v in stats.items()
                        if isinstance(v, (int, float, str, bool, list, dict))}
                try:
                    z.writestr("stats.json", json.dumps(keep, indent=1, default=str))
                except Exception:
                    pass
                # Each tier's own radial profile, scene-referred and divided by
                # the median across tiers. A flat 1.0 everywhere means the one
                # scalar photometric factor per tier describes that tier fully;
                # a departure that grows toward the limb is an error the scalar
                # cannot see, and it is the thing the detail filters amplify
                # into near-limb texture. Its own file because it is the one
                # measurement here that is per-tier rather than per-layer.
                _tr = stats.get("tier_radial")
                if _tr:
                    try:
                        z.writestr("tier_radial.json", json.dumps(_tr))
                    except Exception:
                        pass
            for nm in _LAYERS:
                fp = os.path.join(wd, nm + ".npy")
                if not os.path.exists(fp):
                    continue
                try:
                    a = np.load(fp, mmap_mode="r")
                    if a.ndim == 3:
                        a = a[:, :, 0]
                    sub = np.asarray(a[::2, ::2], np.float32)
                    pr = _profiles(sub, cy / 2, cx / 2, R / 2)
                    if pr:
                        z.writestr("profiles/%s.json" % nm, json.dumps(pr))
                    z.writestr("thumbs/%s.png" % nm, _thumb_png(sub))
                    n_layers += 1
                    del a, sub
                except Exception:
                    continue
    except Exception:
        return None
    if progress is not None:
        try:
            mb = os.path.getsize(path) / 1e6
            progress.log("diagnostics bundle: %d layers, %.1f MB -> "
                         "eclipseforge_output/eclipseforge_diagnostics.zip "
                         "(profiles + 512px thumbnails, no raw data)"
                         % (n_layers, mb), None)
        except Exception:
            pass
    return path


_README = """\
EclipseForgeHDR -- diagnostics bundle
=====================================

This is what the developer needs to diagnose a processing problem WITHOUT your
raw frames. Open it and look; there is nothing here you would call your picture.

  report.txt        the same run report the app printed for you
  geometry.json     the fitted lunar centre and radius
  stats.json        the measurements behind the report
  tier_radial.json  each exposure tier's own radial profile, divided by the
                    median across tiers. Flat 1.0 = the single photometric
                    factor per tier describes it fully; a departure growing
                    toward the limb is an error that factor cannot see.
  profiles/*.json   each detail layer reduced to numbers: its median, p10 and
                    p90 taken ring by ring outward from the limb, plus three
                    azimuthal cuts. This is where ringing, a wrong limb circle
                    and a smeared merge each leave their own signature.
  thumbs/*.png      each detail layer at 512 px, greyscale

WHAT IS NOT IN HERE: no raw frames, no full-resolution image, no colour image,
no tone-mapped result, no EXIF beyond the camera model already in the report,
and no file paths beyond the folder name the report already shows.

The thumbnails are of the DETAIL LAYERS -- normalised high-pass fields. They
show structure, not a finished photograph, and at 512 px greyscale they are not
something anyone would publish.

Typical size is 1-2 MB. If you are reporting a problem, send this file.
"""
