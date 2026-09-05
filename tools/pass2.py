"""Merge under plain / taper / taper+k and render MGN-like views to look at."""
import json, os, re, sys
import numpy as np, cv2
cv2.setNumThreads(4)
F, KJ, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
def srgb(x): return np.where(x <= 0.04045, x/12.92, ((x+0.055)/1.055)**2.4)
def blur(a, s):
    k = int(2*round(3*s)+1)
    return cv2.GaussianBlur(a, (k, k), s, borderType=cv2.BORDER_REPLICATE)
def secs(n):
    t = re.match(r"tier_(.+?)_srgb\.tif$", n).group(1)
    return 1.0/float(t[2:].rstrip("s").replace("p",".")) if t.startswith("1_") \
        else float(t.rstrip("s").replace("p","."))
geo = json.load(open(os.path.join(F, ".eclipseforgehdr", "geometry.json")))
R = float(geo["R"]); cy = float(geo["cy"]); cx = float(geo["cx"])
oy, ox = geo.get("crop_origin", [0, 0])
cal = {float(k): float(v) for k, v in geo["cal"].items()}
sigma = float(np.clip(0.032*R, 8.0, 40.0))
KD = json.load(open(KJ)); NB = KD["NB"]
K = {float(k): np.array(v, np.float32) for k, v in KD["K"].items()}
ALPHA = 0.55            # this set's measured exposure exponent
tdir = os.path.join(F, "eclipseforge_output", "aligned_tiers")
files = sorted(f for f in os.listdir(tdir) if f.endswith("_srgb.tif"))
im0 = cv2.imread(os.path.join(tdir, files[0]), cv2.IMREAD_UNCHANGED)
H, W = im0.shape[:2]; del im0
yy = np.arange(H, dtype=np.float32)[:, None]-(cy-oy)
xx = np.arange(W, dtype=np.float32)[None, :]-(cx-ox)
rad = np.sqrt(yy*yy+xx*xx); del yy, xx
E = np.arange(1.00, 2.60, 0.01, dtype=np.float32)*R
ib = np.clip(np.digitize(rad.ravel(), E)-1, 0, NB-1).reshape(H, W)
outside = (rad < 1.00*R) | (rad >= 2.59*R)
HW = (W//2, H//2)
NAMES = ("plain", "taper", "taper_k")
acc = {n: [np.zeros(HW[::-1], np.float32), np.zeros(HW[::-1], np.float32)]
       for n in NAMES}
for fn in files:
    im = cv2.imread(os.path.join(tdir, fn), cv2.IMREAD_UNCHANGED)
    a = srgb(im.astype(np.float32)/65535.0); del im
    cm = a.max(axis=2); s = secs(fn)
    lum = a.mean(axis=2)/np.float32(s*cal.get(s, 1.0)); del a
    ws = 0.5*(1.0+np.tanh((0.87-cm)/0.06)).astype(np.float32)
    ws[cm > 0.97] = 0.0
    ws *= 0.5*(1.0+np.tanh((cm-0.005)/0.0025))
    ws = ws.astype(np.float32)
    v = (cm <= 0.97).astype(np.float32)
    den = blur(v, sigma); num = blur(ws*v, sigma)
    o = np.where(den > 1e-3, num/np.maximum(den, 1e-3), 0.0).astype(np.float32)
    t = np.clip((den-0.5)*2.0, 0.0, 1.0)
    wt = (s**ALPHA)*o*t*t*v
    kk = K.get(s)
    lk = lum*np.where(outside, np.float32(1.0), kk[ib]) if kk is not None else lum
    for nm, w, L in (("plain", (s**ALPHA)*blur(ws, sigma), lum),
                     ("taper", wt, lum), ("taper_k", wt, lk)):
        acc[nm][0] += cv2.resize(w*L, HW, interpolation=cv2.INTER_AREA)
        acc[nm][1] += cv2.resize(w,   HW, interpolation=cv2.INTER_AREA)
    del cm, lum, ws, v, den, num, o, t, wt, lk
    print("merged", fn, flush=True)
hcy, hcx, hR = (cy-oy)/2.0, (cx-ox)/2.0, R/2.0
Hh, Wh = acc["taper"][0].shape
yy = np.arange(Hh, dtype=np.float32)[:, None]-hcy
xx = np.arange(Wh, dtype=np.float32)[None, :]-hcx
radh = np.sqrt(yy*yy+xx*xx); del yy, xx
disc = radh < 1.01*hR
def mgn(lum):
    x = np.log(np.maximum(lum, 1e-7)).astype(np.float32)
    fill = np.where(disc, np.median(x[~disc]), x).astype(np.float32)
    out = np.zeros_like(fill)
    for k in (4, 8, 16, 32, 64):
        m = blur(fill, k); d = fill-m
        sd = np.sqrt(np.maximum(blur(d*d, k), 1e-12))
        out += np.arctan(3.0*d/sd)
    return out/5.0
half = 1100
y0, x0 = max(int(hcy)-half, 0), max(int(hcx)-half, 0)
sl = (slice(y0, y0+2*half), slice(x0, x0+2*half))
imgs, prof = {}, {}
for nm in NAMES:
    lum = acc[nm][0]/np.maximum(acc[nm][1], 1e-9)
    imgs[nm] = mgn(lum)[sl]
    prof[nm] = [float(np.median(lum[(radh >= f*hR-1) & (radh < f*hR+1)]))
                for f in (1.02, 1.06, 1.20, 1.50)]
print("\nradial profile, relative to taper (1.02 1.06 1.20 1.50 R):")
for nm in NAMES:
    print(f"  {nm:9s} " + " ".join(f"{a/b:6.3f}" for a, b in zip(prof[nm], prof["taper"])))
os.makedirs(OUT, exist_ok=True)
def png(a, p, lo=-1.2, hi=1.2):
    cv2.imwrite(p, (np.clip((a-lo)/(hi-lo), 0, 1)*255).astype(np.uint8))
for nm in NAMES:
    png(imgs[nm], os.path.join(OUT, f"nico_{nm}.png"))
png(imgs["taper"]-imgs["plain"], os.path.join(OUT, "nico_taper_minus_plain.png"), -0.5, 0.5)
png(imgs["taper_k"]-imgs["taper"], os.path.join(OUT, "nico_k_minus_taper.png"), -0.5, 0.5)
rr = radh[sl]
for nm, d in (("taper-plain", imgs["taper"]-imgs["plain"]),
              ("k-taper", imgs["taper_k"]-imgs["taper"])):
    s_ = " ".join(f"{lo:.2f}-{hi:.2f}R {np.sqrt(np.nanmean(d[(rr>=lo*hR)&(rr<hi*hR)]**2)):.3f}"
                  for lo, hi in ((1.02,1.10),(1.10,1.25),(1.25,1.50)))
    print(f"  rms {nm:12s} {s_}", flush=True)
print("wrote", OUT)
