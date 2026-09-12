"""Does each tier read HIGH just outside its own saturated region?

Every correction so far treated a tier as globally wrong (a scale, a colour, a
nonlinearity). None of them can touch an error whose SHAPE is the saturated
region itself -- charge spill from saturated photosites, or veiling glare off a
large saturated area. Both raise the pixels just outside the saturation
contour, by an amount that decays with distance from it. That is a halo bounded
by exactly the curve the rings follow, it lives in the tier data, and no
reweighting can remove it.

Test: for each adjacent tier pair, the longer tier's radiance estimate divided
by the shorter one's, binned by DISTANCE OUTSIDE the longer tier's saturated
region. Flat 1.0 = no spill.
"""
import json, os, re, sys
import numpy as np, cv2
cv2.setNumThreads(4)
F = sys.argv[1]
def srgb(x): return np.where(x <= 0.04045, x/12.92, ((x+0.055)/1.055)**2.4)
def secs(n):
    t = re.match(r"tier_(.+?)_srgb\.tif$", n).group(1)
    return 1.0/float(t[2:].rstrip("s").replace("p",".")) if t.startswith("1_") \
        else float(t.rstrip("s").replace("p","."))
geo = json.load(open(os.path.join(F, ".eclipseforgehdr", "geometry.json")))
cal = {float(k): float(v) for k, v in geo["cal"].items()}
tdir = os.path.join(F, "eclipseforge_output", "aligned_tiers")
files = sorted((secs(f), f) for f in os.listdir(tdir) if f.endswith("_srgb.tif"))
D = np.array([0, 2, 4, 8, 16, 32, 64, 128, 256])
print(f"{'pair (long/short)':22s} " +
      " ".join(f"{a:>3d}-{b:<3d}" for a, b in zip(D[:-1], D[1:])) + "   px outside")
prev = None
for s, fn in files:
    im = cv2.imread(os.path.join(tdir, fn), cv2.IMREAD_UNCHANGED)
    a = srgb(im.astype(np.float32)/65535.0); del im
    cur = (s, a.max(axis=2), a.mean(axis=2)/np.float32(s*cal.get(s, 1.0)))
    del a
    if prev is not None:
        (s0, c0, l0), (s1, c1, l1) = prev, cur       # s1 longer, s0 shorter
        sat = (c1 > 0.97).astype(np.uint8)
        if 5000 < sat.sum() < sat.size*0.9:
            # distance from the saturated region, measured outside it
            dist = cv2.distanceTransform(1-sat, cv2.DIST_L2, 3)
            m = (c1 <= 0.90) & (c0 >= 0.02) & (l0 > 0) & (l1 > 0)
            row = []
            for lo, hi in zip(D[:-1], D[1:]):
                k = m & (dist >= lo) & (dist < hi)
                row.append(np.median(l1[k]/l0[k]) if k.sum() > 2000 else np.nan)
            print(f"{s1:.5g}s / {s0:.5g}s".ljust(22) +
                  " ".join("  --   " if not np.isfinite(x) else f"{x:7.3f}"
                           for x in row), flush=True)
            del dist, m
        del sat
    prev = cur
print("\n>1 near the contour, falling with distance = the longer tier reads high")
print("just outside its own saturated area. No reweighting can remove that.")
