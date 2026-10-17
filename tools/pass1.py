"""Fit k_i(r) per tier from rings each tier leaves >=98% unclipped. Save it."""
import json, os, re, sys
import numpy as np, cv2
cv2.setNumThreads(4)
F, OUT = sys.argv[1], sys.argv[2]
def srgb(x): return np.where(x <= 0.04045, x/12.92, ((x+0.055)/1.055)**2.4)
def secs(n):
    t = re.match(r"tier_(.+?)_srgb\.tif$", n).group(1)
    return 1.0/float(t[2:].rstrip("s").replace("p",".")) if t.startswith("1_") \
        else float(t.rstrip("s").replace("p","."))
geo = json.load(open(os.path.join(F, ".eclipseforgehdr", "geometry.json")))
R = float(geo["R"]); cy = float(geo["cy"]); cx = float(geo["cx"])
oy, ox = geo.get("crop_origin", [0, 0])
cal = {float(k): float(v) for k, v in geo["cal"].items()}
tdir = os.path.join(F, "eclipseforge_output", "aligned_tiers")
files = sorted(f for f in os.listdir(tdir) if f.endswith("_srgb.tif"))
im0 = cv2.imread(os.path.join(tdir, files[0]), cv2.IMREAD_UNCHANGED)
H, W = im0.shape[:2]; del im0
yy = np.arange(H, dtype=np.float32)[:, None]-(cy-oy)
xx = np.arange(W, dtype=np.float32)[None, :]-(cx-ox)
rad = np.sqrt(yy*yy+xx*xx).ravel(); del yy, xx
E = np.arange(1.00, 2.60, 0.01, dtype=np.float32)*R
NB = len(E)-1
idx = (np.digitize(rad, E)-1).astype(np.int32); del rad
inb = (idx >= 0) & (idx < NB); idxc = np.clip(idx, 0, NB-1); del idx
prof = {}
for fn in files:
    im = cv2.imread(os.path.join(tdir, fn), cv2.IMREAD_UNCHANGED)
    a = srgb(im.astype(np.float32)/65535.0); del im
    cm = a.max(axis=2).ravel(); s = secs(fn)
    v = (a.mean(axis=2).ravel()/np.float32(s*cal.get(s, 1.0))).astype(np.float64)
    del a
    tot = np.bincount(idxc[inb], minlength=NB).astype(np.float64)
    clp = np.bincount(idxc[inb & (cm > 0.97)], minlength=NB).astype(np.float64)
    clean = (tot > 400) & (clp/np.maximum(tot, 1) < 0.02)
    ok = inb & (cm >= 0.02) & clean[idxc]
    cnt = np.bincount(idxc[ok], minlength=NB).astype(np.float64)
    sm = np.bincount(idxc[ok], weights=v[ok], minlength=NB)
    prof[s] = np.where(cnt > 400, sm/np.maximum(cnt, 1), np.nan).tolist()
    print("profiled", fn, int(np.isfinite(prof[s]).sum()), "rings", flush=True)
    del cm, v, ok, cnt, sm
P = np.array([prof[s] for s in sorted(prof)])
ref = np.nanmedian(P, axis=0)
K = {}
for s in sorted(prof):
    p = np.array(prof[s])
    k = np.where(np.isfinite(p) & np.isfinite(ref) & (p > 0),
                 ref/np.maximum(p, 1e-12), 1.0)
    k = np.clip(np.nan_to_num(k, nan=1.0), 0.25, 4.0)
    kk = np.convolve(k, np.ones(9)/9.0, mode="same")
    kk[:4] = kk[4]; kk[-4:] = kk[-5]
    K[str(s)] = kk.tolist()
    print(f"  tier {s:g}s  max |k-1| inside 1.5R = "
          f"{100*np.max(np.abs(kk[:50]-1)):.0f}%", flush=True)
json.dump({"R": R, "NB": NB, "K": K}, open(OUT, "w"))
print("wrote", OUT)
