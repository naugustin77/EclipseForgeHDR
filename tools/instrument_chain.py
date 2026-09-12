#!/usr/bin/env python3
"""Dump the per-frame chain, stage by stage, so the ring artifact can be traced
to the step that creates it.

WHY THIS EXISTS. a tester mean-stacked two sets of the same nine raws in Photoshop:
EclipseForgeHDR's exported aligned tiers, and Lightroom's exports. The EFHDR
stack rings; the Lightroom stack does not. So something in our per-frame
handling makes it. Reproducing that chain from outside the app was tried and
produced three confounded results in a row (no white balance, a shift applied
to the mosaic instead of the channels). The fix is not a better reproduction --
it is to stop reproducing and instrument the real code.

This script imports the app's own modules and walks one tier through the same
operations the merge loop performs, saving the MEAN over all tiers after each
stage. The mean is what the tester stacked, so the artifact appears in exactly the
same form; whichever stage it first appears in is the one that makes it.

Stages, in the order pipeline.run() applies them:

    0_bayer        raw, black-subtracted, rolled to RGGB      (raw.RawFile)
    1_flat         after the master flat is divided out
    2_demosaic     after demosaic_rggb                        (luminance)
    3_pedestal     after the shared pedestal subtraction
    4_wb           after the white-balance gains
    5_matrix       after the camera -> sRGB colour matrix
    6_shift        after the sub-pixel alignment shift
    7_export       clipped and sRGB-encoded, as _write_tier_tiff writes it

Everything is cropped to a window about the disc so the files stay small, and
only the mean over the nine tiers is written -- 8 files, ~10 MB each.

RUN IT WITH THE SAME PYTHON THAT RUNS THE APP (it needs rawpy):

    python3 tools/instrument_chain.py "/path/to/the/raw/folder"

It reads geometry.json and masterflat.npy from that folder's .eclipseforgehdr,
so run it on a folder the app has already processed. Nothing is modified: the
output goes to <folder>/eclipseforge_output/chain_stages/.

NOT reproduced here, deliberately, and worth remembering when reading the
result: hot-pixel repair (needs the map built across the four shortest tiers)
and the LDIC azimuthal correction (needs the fit against the running
composite). Stage 6 is therefore the tier as it stands just BEFORE LDIC, while
the tier TIFF the app writes has LDIC applied. If the arcs are already in an
earlier stage, neither matters.
"""
import json
import os
import sys

import numpy as np

HALF = 800          # crop half-width about the disc centre, full-res px


def main(folder):
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                    "..", "src"))
    from scipy import ndimage
    from eclipseforgehdr.raw import (RawFile, demosaic_rggb, list_raws,
                                     read_exif)

    wd = os.path.join(folder, ".eclipseforgehdr")
    geo = json.load(open(os.path.join(wd, "geometry.json")))
    flat = None
    fp = os.path.join(wd, "masterflat.npy")
    if os.path.exists(fp):
        flat = np.load(fp)
        print(f"master flat {flat.shape}")
    cal = {float(k): float(v) for k, v in geo["cal"].items()}
    shifts = {float(k): v for k, v in geo["abs_shift"].items()}
    oy, ox = geo.get("crop_origin", [0, 0])
    cy, cx = float(geo["cy"]) + oy, float(geo["cx"]) + ox
    print(f"disc centre in the untrimmed frame: ({cx:.1f}, {cy:.1f})")

    # The shared pedestal is fitted during the run and not cached; it is a few
    # tenths of an ADU and cannot make a contour, so it is reported and skipped
    # rather than guessed at.
    pedestal = float(os.environ.get("EFHDR_PEDESTAL", "0.0"))

    files = sorted(list_raws(folder))
    print(f"{len(files)} raw files")

    stages = ["0_bayer", "1_flat", "2_demosaic", "3_pedestal", "4_wb",
              "5_matrix", "6_shift", "7_export"]
    acc = {k: None for k in stages}
    n = 0
    # pipeline.py takes wb, cam2rgb and sat_level from the FIRST frame and uses
    # them for every tier (color_info, set once at line ~1791). Do the same, so
    # this is the app's arithmetic and not a per-frame variant of it.
    cinfo = {}

    def add(key, img):
        c = crop(img)
        if acc[key] is None:
            acc[key] = c.astype(np.float64)
        else:
            acc[key] += c

    def crop(a):
        y0, x0 = int(round(cy)) - HALF, int(round(cx)) - HALF
        y0 = max(y0 - y0 % 2, 0)
        x0 = max(x0 - x0 % 2, 0)
        return np.ascontiguousarray(a[y0:y0 + 2 * HALF, x0:x0 + 2 * HALF])

    for path in files:
        sec = float(read_exif(path)[0])
        rf = RawFile(path)
        # geometry.json keys the tiers by the exposure as a float; match on the
        # nearest key so 1/60 s does not miss 0.016666666666666666
        key = min(cal, key=lambda k: abs(k - sec)) if cal else sec
        if abs(key - sec) > 0.02 * max(sec, 1e-9):
            print(f"    !! no tier within 2% of {sec:g}s -- using cal 1.0")
            key = None
        c = cal.get(key, 1.0)
        dy, dx = shifts.get(key, (0.0, 0.0))
        print(f"  {os.path.basename(path)}  {sec:g}s  cal {c:.4f}  "
              f"shift ({2*dy:+.2f},{2*dx:+.2f}) px", flush=True)

        if not cinfo:
            cinfo = {"wb": np.asarray(rf.daylight_wb, np.float32),
                     "cam2rgb": np.asarray(rf.cam2rgb, np.float32),
                     "sat_level": float(rf.sat_level)}
            print(f"    wb {cinfo['wb']}  sat_level {cinfo['sat_level']:.1f}")
        bay = rf.bayer.astype(np.float32)
        add("0_bayer", bay)
        if flat is not None and flat.shape == bay.shape:
            bay = bay / flat
        add("1_flat", bay)

        rgb = demosaic_rggb(bay)
        del bay
        add("2_demosaic", rgb.mean(axis=2))
        rgb -= np.float32(pedestal)
        add("3_pedestal", rgb.mean(axis=2))
        rgb *= cinfo["wb"][None, None, :]
        add("4_wb", rgb.mean(axis=2))
        h, w = rgb.shape[:2]
        rgb = (rgb.reshape(-1, 3) @ cinfo["cam2rgb"].T).reshape(h, w, 3)
        add("5_matrix", rgb.mean(axis=2))
        rgb /= np.float32(sec * c)
        for ch in range(3):
            rgb[:, :, ch] = ndimage.shift(rgb[:, :, ch], (2 * dy, 2 * dx),
                                          order=1, mode="nearest")
        add("6_shift", rgb.mean(axis=2))
        # what _write_tier_tiff writes, minus the LDIC step it cannot know here
        a = np.clip(rgb * np.float32(sec * c / max(cinfo["sat_level"], 1e-9)), 0.0, 1.0)
        add("7_export", a.mean(axis=2))
        del rgb, a, rf
        n += 1

    out = os.path.join(folder, "eclipseforge_output", "chain_stages")
    os.makedirs(out, exist_ok=True)
    meta = {"n": n, "half": HALF, "R": float(geo["R"]),
            "cy": float(cy) - (int(round(cy)) - HALF),
            "cx": float(cx) - (int(round(cx)) - HALF),
            "pedestal_used": pedestal, "folder": folder}
    for k in stages:
        np.save(os.path.join(out, f"mean_{k}.npy"),
                (acc[k] / n).astype(np.float32))
    json.dump(meta, open(os.path.join(out, "meta.json"), "w"), indent=1)
    print(f"\nwrote {len(stages)} stage means to {out}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(os.path.expanduser(sys.argv[1]))
