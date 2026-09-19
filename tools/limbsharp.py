"""How sharp is the lunar limb -- in ONE RAW FRAME, and in the merge built
from them?

The question this answers: when a set looks soft at 1:1, is the softness in
the data (optics, focus, seeing, shutter shake) or does the pipeline lose it
(alignment, intra-tier averaging, the HDR merge)? Both ends are measured the
same way, in the same units, so the two numbers can be put side by side.

THE MEASURE is the 20-80% transition width of the lunar limb along 720 rays:
the distance over which the profile climbs from 20% to 80% of the step between
the occulted disc's own level and the near-limb corona. The Moon's edge is the
sharpest thing in the frame -- an occulting disc with no atmosphere of its own
-- so its width IS the system's point spread, everything included. Reported
per ray as the median and the 90th percentile, in full-resolution pixels and
in arcseconds (the plate scale comes from the fitted lunar radius against the
Moon's apparent diameter, the same assumption the run report uses).

    python tools/limbsharp.py FRAME.CR2 [FRAME2.CR2 ...]
    python tools/limbsharp.py --folder /path/to/bracket     # + the merge
    python tools/limbsharp.py --folder /path --lab          # lab workdir

A single frame is read at its own Bayer green (half resolution, no demosaic,
no calibration): a demosaic interpolates across the very edge being measured,
which would soften it by about a pixel before anything else got the chance.
The half-res width is doubled so both numbers are in full-resolution pixels.

WHAT THE COMPARISON MEANS. If one frame and the merge agree, the pipeline is
carrying the sharpness it was given and the limit is the data. If the merge is
wider, the loss is between them -- alignment residual, lunar motion inside a
tier, or the merge weight -- and the run report's own numbers (network
residual, lunar track, tier shifts) say which.
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "src"))
from eclipseforgehdr.pipeline import find_disc, fit_limb_rays   # noqa: E402

MOON_ARCSEC = 1920.0            # apparent diameter, +/-4% over a lunation


def green_half(path):
    """The Bayer green channel of one raw frame, at half resolution."""
    import rawpy
    with rawpy.imread(path) as r:
        cfa = r.raw_image_visible.astype(np.float32)
        black = float(np.mean(r.black_level_per_channel))
        cols = r.raw_colors_visible[:2, :2]
    g = [(i, j) for i in range(2) for j in range(2) if cols[i, j] == 1]
    if len(g) != 2:
        g = [(0, 1), (1, 0)]
    H, W = cfa.shape
    H -= H % 2
    W -= W % 2
    a = cfa[g[0][0]:H:2, g[0][1]:W:2]
    b = cfa[g[1][0]:H:2, g[1][1]:W:2]
    return np.maximum(0.5 * (a + b) - black, 0.0)


def limb_width(lum, cy, cx, R, n_ang=720):
    """20-80% transition width per ray, in the input's own pixels.

    Sampled at 1/8 px along each ray and normalised by that ray's OWN step, so
    a bright sector cannot widen or narrow the number: this is a width, not a
    contrast.
    """
    ang = np.linspace(0, 2 * np.pi, n_ang, endpoint=False)
    rr = np.arange(0.75 * R, 1.35 * R, 0.125, dtype=np.float32)
    ys = cy + rr[None, :] * np.sin(ang)[:, None]
    xs = cx + rr[None, :] * np.cos(ang)[:, None]
    from scipy import ndimage
    prof = ndimage.map_coordinates(lum, np.stack([ys, xs]), order=1,
                                   mode="nearest")
    lo = np.median(prof[:, rr < 0.9 * R], axis=1)          # inside the disc
    hi = np.median(prof[:, (rr > 1.05 * R) & (rr < 1.12 * R)], axis=1)
    out = []
    for k in range(n_ang):
        if not np.isfinite(hi[k]) or hi[k] - lo[k] < 1e-6:
            continue
        p = (prof[k] - lo[k]) / (hi[k] - lo[k])
        i20 = np.argmax(p > 0.2)
        i80 = np.argmax(p > 0.8)
        if i80 <= i20 or i20 == 0:
            continue
        out.append((rr[i80] - rr[i20]))
    return np.array(out)


def report(name, lum, scale_px, seed=None):
    d = find_disc(lum) if seed is None else seed
    if d is None:
        print(f"{name}: no disc found")
        return
    cy, cx, R = d
    cy, cx, R, rms, kept, tot, _ = fit_limb_rays(lum, cy, cx, R)
    w = limb_width(lum, cy, cx, R)
    if not w.size:
        print(f"{name}: limb found (R={R:.1f}) but no usable rays")
        return
    asec = MOON_ARCSEC / (2.0 * R)          # arcsec per pixel of THIS image
    med, p90 = float(np.median(w)), float(np.percentile(w, 90))
    print(f"{name}:")
    print(f"   limb fit   R = {R * scale_px:7.1f} full-res px, "
          f"rms {rms * scale_px:.2f} px over {kept}/{tot} rays")
    print(f"   20-80%     {med * scale_px:5.2f} px (p90 {p90 * scale_px:5.2f})"
          f"   = {med * asec:5.2f}\" (p90 {p90 * asec:5.2f}\")")


def main(argv):
    folder = None
    lab = False
    frames = []
    i = 0
    while i < len(argv):
        a = argv[i]
        if a == "--folder":
            i += 1
            folder = argv[i]
        elif a == "--lab":
            lab = True
        else:
            frames.append(a)
        i += 1
    if folder:
        from eclipseforgehdr.raw import list_raws
        frames = frames or list_raws(folder)
    if not frames and not folder:
        print(__doc__)
        return 1
    for f in frames:
        try:
            g = green_half(f)
        except Exception as e:                    # noqa: BLE001
            print(f"{os.path.basename(f)}: cannot read ({e})")
            continue
        # half-res image, so every length doubles to be a full-res figure
        report(os.path.basename(f), g, 2.0)
    if folder:
        wd = os.path.join(folder, ".eclipseforgehdr-lab" if lab
                          else ".eclipseforgehdr")
        p = os.path.join(wd, "hdr_lum.npy")
        if os.path.exists(p):
            print()
            report("MERGED (" + os.path.basename(wd) + ")",
                   np.load(p, mmap_mode="r")[:], 1.0)
        else:
            print(f"\nno merged luminance in {wd} — stack the folder first")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
