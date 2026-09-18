"""Regression test for the merge, on a synthetic bracket. No rawpy needed.

Writes a small 2-D Bayer FITS bracket (four exposure groups, known 5.000x steps,
a reddish r^-3 corona, a blue sky that scales with exposure, a black level,
per-frame shifts, saturation on the long tiers), optionally with a vignetting
flat set beside it, runs simple.run with the layer builder stubbed out, and
checks the merged HDR against the scene it was made from.

    python3 tools/synthetic_bracket_test.py            # from the repo root

A 2-D FITS goes through FitsFrame -> rf.bayer, so this exercises the Bayer
path of simple._load -- saturation masking, demosaic, white balance (unity for
FITS), colour matrix (identity for FITS), calibration -- without a camera raw.
"""
import os
import shutil
import sys

import numpy as np
from scipy.ndimage import shift as _sh

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))
WORK = os.path.join(HERE, "_synthetic_bracket")


def write_fits(path, data_u16, exptime, sat=65535):
    """Minimal FITS writer: BITPIX 16 with BZERO 32768 (unsigned convention)."""
    H, W = data_u16.shape
    cards = ["SIMPLE  =                    T", "BITPIX  =                   16",
             "NAXIS   =                    2", "NAXIS1  = %20d" % W,
             "NAXIS2  = %20d" % H, "BZERO   =                32768",
             "BSCALE  =                    1", "EXPTIME = %20.6f" % exptime,
             "SATURATE= %20d" % sat, "BAYERPAT= 'RGGB    '", "END"]
    hdr = "".join(c.ljust(80) for c in cards)
    hdr += " " * ((2880 - len(hdr) % 2880) % 2880)
    d = (data_u16.astype(np.int32) - 32768).astype(">i2").tobytes()
    d += b"\0" * ((2880 - len(d) % 2880) % 2880)
    with open(path, "wb") as f:
        f.write(hdr.encode("ascii") + d)


def make_scene(folder, vignette=None, seed=1):
    H, W = 480, 640
    cy, cx, R = 240.0, 320.0, 40.0
    yy, xx = np.mgrid[0:H, 0:W].astype(np.float64)
    r = np.hypot(yy - cy, xx - cx) / R
    with np.errstate(divide="ignore"):
        corona = np.where(r > 1.0, 4000.0 * np.maximum(r, 1e-3) ** -3.0, 0.0)
    corona_rgb = np.stack([corona * 1.3, corona, corona * 0.6])
    sky = np.array([12.0, 17.0, 22.0])
    black = 512.0
    tiers = [(0.002, 4), (0.010, 4), (0.050, 3), (0.250, 2)]
    rng = np.random.default_rng(seed)
    shutil.rmtree(folder, ignore_errors=True)
    os.makedirs(folder)
    k = 0
    for t, n in tiers:
        for _ in range(n):
            dy, dx = rng.uniform(-3, 3, 2)
            img = np.empty((H, W))
            for (oy, ox, c) in ((0, 0, 0), (0, 1, 1), (1, 0, 1), (1, 1, 2)):
                plane = _sh((corona_rgb[c] + sky[c]) * t * 1000.0, (dy, dx),
                            order=1, mode="nearest")
                img[oy::2, ox::2] = plane[oy::2, ox::2]
            if vignette is not None:
                img *= vignette
            img += black + rng.normal(0, 3.0, img.shape)
            img = np.clip(np.round(img), 0, 65535).astype(np.uint16)
            write_fits(os.path.join(folder, "frame_%03d.fit" % k), img, t)
            k += 1
    return dict(r=r, xx=xx, cx=cx, corona_rgb=corona_rgb, sky=sky, black=black,
                yy=yy, cy=cy, rng=rng)


def run_merge(folder, flat_dir):
    from eclipseforgehdr import simple, importhdr
    cap = {}
    importhdr.build_from_rgb = lambda folder, rgb, progress, **kw: cap.update(rgb=rgb, kw=kw)

    class P:
        def log(self, m, frac=None):
            pass
    simple.run(folder, P(), flat_dir=flat_dir)
    return cap["rgb"], cap["kw"]["stats"], cap["kw"]["opts"]


def main():
    fails = 0

    def check(cond, msg):
        nonlocal fails
        print(("  ok    " if cond else "  FAIL  ") + msg)
        fails += 0 if cond else 1

    print("1. plain bracket")
    g = make_scene(WORK)
    rgb, st, op = run_merge(WORK, "nonexistent")
    lad = [x["measured"] for x in st["simple_ladder"]]
    check(np.isfinite(rgb).all(), "merged HDR is finite everywhere")
    check(all(abs(x - 5.0) < 0.05 for x in lad),
          "ladder %s against a true 5.000 per step" % ["%.3f" % x for x in lad])
    r, sky, cr = g["r"], g["sky"], g["corona_rgb"]
    m = (r > 3.0) & (r < 4.0)
    far = rgb[m].mean(0) / rgb[m].mean(0)[1]
    exp = (cr[:, m].mean(1) + sky) / (cr[1, m].mean() + sky[1])
    check(abs(far[0] - exp[0]) < 0.03 and abs(far[2] - exp[2]) < 0.03,
          "far-field colour R/G %.3f B/G %.3f, scene %.3f %.3f (sky left in)"
          % (far[0], far[2], exp[0], exp[2]))
    check(op["wb_source"] == "camera" and not op["flat_applied"],
          "opts record white balance and no flat")

    print("2. vignetted bracket with a flat set beside it")
    yy, xx, cy, cx = g["yy"], g["xx"], g["cy"], g["cx"]
    vig = 1.0 - 0.35 * ((yy - cy) ** 2 + (xx - cx - 150) ** 2) / (cy ** 2 + cx ** 2)
    m1 = (r > 3.0) & (r < 3.6) & (xx < cx)
    m2 = (r > 3.0) & (r < 3.6) & (xx >= cx)
    base = rgb[m1][:, 1].mean() / rgb[m2][:, 1].mean()
    g2 = make_scene(WORK, vignette=vig)
    rgb_v, _, _ = run_merge(WORK, "nonexistent")
    lr_v = rgb_v[m1][:, 1].mean() / rgb_v[m2][:, 1].mean()
    os.makedirs(os.path.join(WORK, "flats"))
    for i in range(5):
        write_fits(os.path.join(WORK, "flats", "flat_%d.fit" % i),
                   np.clip(np.round(30000 * vig + g2["black"]
                                    + g2["rng"].normal(0, 20, vig.shape)),
                           0, 65535).astype(np.uint16), 0.01)
    shutil.rmtree(os.path.join(WORK, ".eclipseforgehdr"), ignore_errors=True)
    rgb_f, st_f, op_f = run_merge(WORK, "")
    lr_f = rgb_f[m1][:, 1].mean() / rgb_f[m2][:, 1].mean()
    check(op_f["flat_applied"], "flat built and applied")
    check(st_f["calibration"].get("black_pre_flat") is not None,
          "black level taken off before the flat (this FITS reports none)")
    check(abs(lr_f - base) < 0.02,
          "far-field left/right: %.3f unvignetted, %.3f vignetted, %.3f with "
          "the flat" % (base, lr_v, lr_f))
    shutil.rmtree(WORK, ignore_errors=True)
    print("\n%s" % ("ALL PASSED" if not fails else "%d FAILED" % fails))
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
