#!/usr/bin/env python3
"""Rebuild the merge from exported aligned tiers under each feather variant and
score the result BOTH ways: radial leak and azimuthal contour structure.

Every previous test of _feather_weight measured a radial profile or a local
weight step. The artifact Nico sees is azimuthal -- it traces each tier's
saturation contour, which follows the corona, so it is fingered and lobed, not
circular. This scores that directly.

Input: <folder>/eclipseforge_output/aligned_tiers/tier_*_srgb.tif
       <folder>/.eclipseforgehdr/geometry.json   (R, cal, secs)

The sRGB export holds a = clip(raw / sat_level, 0, 1) per channel, so after an
sRGB decode every threshold in the merge becomes dimensionless: knee 0.87,
hard clip 0.97, floor 0.005. Nothing about sat_level needs to be recovered.

Not modelled: moon_weight (acts inside the disc, not in 1.05-1.5 R) and the
pedestal (a constant, identical across variants). Both cancel in a comparison.
"""
import json, os, re, sys, time
import numpy as np
import cv2

cv2.setNumThreads(4)

FOLDER = sys.argv[1]
OUT = sys.argv[2] if len(sys.argv) > 2 else "/tmp/feather_bench.json"
MERGE_FLOOR = 0.005
ALPHAS = (1.0, 0.55)
VARIANTS = ("none", "plain", "norm", "masked", "taper", "smooth")


def log(*a):
    print(f"[{time.strftime('%H:%M:%S')}]", *a, flush=True)


def srgb_to_linear(x):
    return np.where(x <= 0.04045, x / 12.92, ((x + 0.055) / 1.055) ** 2.4)


def blur(a, s):
    k = int(2 * round(3 * s) + 1)
    return cv2.GaussianBlur(a, (k, k), s, borderType=cv2.BORDER_REPLICATE)


def secs_from_name(n):
    m = re.match(r"tier_(.+?)_(srgb|linear)\.tif$", n)
    t = m.group(1)
    if t.startswith("1_"):
        return 1.0 / float(t[2:].rstrip("s").replace("p", "."))
    return float(t.rstrip("s").replace("p", "."))


def feather(mode, w, cmax, sigma):
    """The five weight forms. `w` is wsat as the pipeline builds it; `cmax` is
    the per-pixel channel max in units of saturation."""
    if mode == "none":
        return w
    if mode == "plain":                       # <= 0.22.15
        return blur(w, sigma)
    if mode == "smooth":                      # candidate: no mask anywhere
        # A weight that is a pointwise function of a SMOOTHED intensity, rather
        # than a smoothed function of intensity. Dilating first guarantees that
        # any pixel within sigma of a clipped one sees a clipped value, so the
        # weight is exactly zero there -- no leak -- and nothing binary is ever
        # multiplied in, so no contour is written into the weight.
        rad = max(int(round(sigma)), 1)
        ker = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (2 * rad + 1,) * 2)
        cs = blur(cv2.dilate(cmax, ker), sigma)
        u = np.clip((0.95 - cs) / 0.15, 0.0, 1.0)     # exact 0 at 0.95 sat
        out = (u * u * (3 - 2 * u)).astype(np.float32)
        lo = MERGE_FLOOR
        out *= 0.5 * (1.0 + np.tanh((blur(cmax, sigma) - lo) / (0.5 * lo)))
        return out.astype(np.float32)
    v = (cmax <= 0.97).astype(np.float32)
    den = blur(v, sigma)
    num = blur(w * v, sigma)
    out = np.where(den > 1e-3, num / np.maximum(den, 1e-3), 0.0).astype(np.float32)
    if mode == "norm":                        # isolates the 1/den gain alone
        return out
    if mode == "masked":                      # 0.22.16
        return out * v
    t = np.clip((den - 0.5) * 2.0, 0.0, 1.0)  # 0.22.19, current
    return out * t * t * v


geo = json.load(open(os.path.join(FOLDER, ".eclipseforgehdr", "geometry.json")))
R = float(geo["R"])
cy, cx = float(geo["cy"]), float(geo["cx"])
oy, ox = geo.get("crop_origin", [0, 0])
cal = {float(k): float(v) for k, v in geo["cal"].items()}
sigma = float(np.clip(0.032 * R, 8.0, 40.0))
log(f"R={R:.1f}  feather sigma={sigma:.1f}  centre=({cx:.1f},{cy:.1f}) crop={oy},{ox}")

tdir = os.path.join(FOLDER, "eclipseforge_output", "aligned_tiers")
files = sorted([f for f in os.listdir(tdir) if f.endswith("_srgb.tif")])
log(f"{len(files)} tiers")

acc = {}   # (variant, alpha) -> [sum w*lum, sum w]  at half resolution
HW = None

for i, fn in enumerate(files):
    s = secs_from_name(fn)
    c = cal.get(s, 1.0)
    img = cv2.imread(os.path.join(tdir, fn), cv2.IMREAD_UNCHANGED)
    if img is None:
        log(f"  !! unreadable {fn}"); continue
    a = srgb_to_linear(img.astype(np.float32) / 65535.0)
    del img
    cmax = a.max(axis=2)
    lum = a.mean(axis=2) / np.float32(s * c)          # relative radiance
    del a
    if HW is None:
        H, W = cmax.shape
        HW = (W // 2, H // 2)
        log(f"  full {W}x{H} -> half {HW[0]}x{HW[1]}")
    # wsat exactly as the pipeline builds it, in units of saturation
    wsat = 0.5 * (1.0 + np.tanh((0.87 - cmax) / 0.06)).astype(np.float32)
    wsat[cmax > 0.97] = 0.0
    wsat *= 0.5 * (1.0 + np.tanh((cmax - MERGE_FLOOR) / (0.5 * MERGE_FLOOR)))
    wsat = wsat.astype(np.float32)
    for mode in VARIANTS:
        fw = feather(mode, wsat, cmax, sigma)
        for al in ALPHAS:
            w = (s ** al) * fw
            k = (mode, al)
            if k not in acc:
                acc[k] = [np.zeros(HW[::-1], np.float32),
                          np.zeros(HW[::-1], np.float32)]
            acc[k][0] += cv2.resize(w * lum, HW, interpolation=cv2.INTER_AREA)
            acc[k][1] += cv2.resize(w, HW, interpolation=cv2.INTER_AREA)
            del w
        del fw
    del wsat, cmax, lum
    log(f"  merged {fn}  ({i+1}/{len(files)})")

# ---- scoring -------------------------------------------------------------
hcy, hcx = (cy - oy) / 2.0, (cx - ox) / 2.0
hR = R / 2.0
NA = 1024
rr = np.arange(0.95 * hR, 2.30 * hR, 0.5, dtype=np.float32)
th = np.arange(NA, dtype=np.float32) * (2 * np.pi / NA)
mapx = (hcx + rr[:, None] * np.cos(th)[None, :]).astype(np.float32)
mapy = (hcy + rr[:, None] * np.sin(th)[None, :]).astype(np.float32)

SHELLS = ((1.05, 1.25), (1.25, 1.50), (1.50, 2.00))
BANDS = {"m20_80": (20, 80), "m80_250": (80, 250)}
RADII = (1.00, 1.02, 1.06, 1.20, 1.50)

res = {}
for (mode, al), (num, den) in acc.items():
    lum = num / np.maximum(den, 1e-9)
    pol = cv2.remap(lum, mapx, mapy, cv2.INTER_LINEAR,
                    borderMode=cv2.BORDER_REPLICATE)
    med = np.median(pol, axis=1, keepdims=True)
    norm = pol / np.maximum(med, 1e-12)
    F = np.abs(np.fft.rfft(norm - norm.mean(axis=1, keepdims=True), axis=1))
    P = (F * F) / NA
    e = {}
    for sn, (lo, hi) in zip(("s1", "s2", "s3"), SHELLS):
        sel = (rr >= lo * hR) & (rr < hi * hR)
        for bn, (m0, m1) in BANDS.items():
            e[f"{sn}_{bn}"] = float(P[sel, m0:m1].mean())
    prof = {}
    for f in RADII:
        j = int(np.argmin(np.abs(rr - f * hR)))
        prof[f"{f:.2f}R"] = float(med[j, 0])
    res[f"{mode}@{al}"] = {"az": e, "prof": prof}
    log(f"scored {mode}@{al}")

json.dump({"R": R, "sigma": sigma, "folder": FOLDER, "res": res},
          open(OUT, "w"), indent=1)
log("wrote", OUT)

# ---- report --------------------------------------------------------------
for al in ALPHAS:
    ref = res[f"none@{al}"]["prof"]
    pl = res[f"plain@{al}"]["az"]
    print(f"\n=== alpha {al} ===")
    print(f"{'variant':9s} " + " ".join(f"{k:>7s}" for k in ref) +
          "   " + " ".join(f"{k:>10s}" for k in pl))
    for mode in VARIANTS:
        r = res[f"{mode}@{al}"]
        pr = " ".join(f"{r['prof'][k]/max(ref[k],1e-12):7.3f}"
                      for k in ref)
        az = " ".join(f"{r['az'][k]/max(pl[k],1e-30):10.2f}" for k in pl)
        print(f"{mode:9s} {pr}   {az}")
print("\nprofile columns: merged level / unfeathered reference (1.000 = no radial leak)")
print("azimuthal columns: ring-normalised power / the plain blur's (1.00 = no contour structure)")
