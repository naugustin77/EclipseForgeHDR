"""Stop optimising a number. Build the two merges that differ ON SCREEN --
plain and the shipped taper -- push both through an MGN-like stretch, and look
at where they differ. No assumed metric.
"""
import json, os, re, sys
import numpy as np, cv2
cv2.setNumThreads(4)
F, OUT = sys.argv[1], sys.argv[2]
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
tdir = os.path.join(F, "eclipseforge_output", "aligned_tiers")
files = sorted(f for f in os.listdir(tdir) if f.endswith("_srgb.tif"))
HW = None; acc = {"plain": None, "taper": None, "cont": None}
for fn in files:
    im = cv2.imread(os.path.join(tdir, fn), cv2.IMREAD_UNCHANGED)
    a = srgb(im.astype(np.float32)/65535.0); del im
    cm = a.max(axis=2); s = secs(fn)
    lum = a.mean(axis=2)/np.float32(s*cal.get(s, 1.0)); del a
    if HW is None:
        H, W = cm.shape; HW = (W//2, H//2)
        for k in acc: acc[k] = [np.zeros(HW[::-1], np.float32),
                                np.zeros(HW[::-1], np.float32)]
    ws = 0.5*(1.0+np.tanh((0.87-cm)/0.06)).astype(np.float32)
    ws[cm > 0.97] = 0.0
    ws *= 0.5*(1.0+np.tanh((cm-0.005)/0.0025))
    ws = ws.astype(np.float32)
    v = (cm <= 0.97).astype(np.float32)
    den = blur(v, sigma); num = blur(ws*v, sigma)
    o = np.where(den > 1e-3, num/np.maximum(den, 1e-3), 0.0).astype(np.float32)
    t = np.clip((den-0.5)*2.0, 0.0, 1.0)
    # CANDIDATE: roll the weight to EXACT zero, C1, at 0.90 sat -- well below
    # the 0.97 validity threshold. The taper machinery still stops the spatial
    # leak, but by the time it multiplies by t*t*v the weight there is already
    # zero, so there is no step left on the contour for it to print. Cost: the
    # brightest 7% of each tier is discarded, and a shorter tier always covers
    # it.
    u = np.clip((0.90 - cm)/0.10, 0.0, 1.0)
    wc = (u*u*(3-2*u)).astype(np.float32)
    wc *= 0.5*(1.0+np.tanh((cm-0.005)/0.0025))
    dc = blur(v, sigma); nc = blur(wc*v, sigma)
    oc = np.where(dc > 1e-3, nc/np.maximum(dc, 1e-3), 0.0).astype(np.float32)
    tc = np.clip((dc-0.5)*2.0, 0.0, 1.0)
    for nm, w in (("plain", s*blur(ws, sigma)), ("taper", s*o*t*t*v),
                  ("cont", s*oc*tc*tc*v)):
        acc[nm][0] += cv2.resize(w*lum, HW, interpolation=cv2.INTER_AREA)
        acc[nm][1] += cv2.resize(w,   HW, interpolation=cv2.INTER_AREA)
    del cm, lum, ws, v, den, num, o, t
    print("merged", fn, flush=True)

hcy, hcx, hR = (cy-oy)/2.0, (cx-ox)/2.0, R/2.0
Hh, Wh = acc["plain"][0].shape
yy = np.arange(Hh, np.float32)[:, None]-hcy if False else (np.arange(Hh, dtype=np.float32)[:, None]-hcy)
xx = np.arange(Wh, dtype=np.float32)[None, :]-hcx
rad = np.sqrt(yy*yy+xx*xx)
disc = rad < 1.01*hR

def mgn(lum):
    """MGN in miniature: local contrast normalised at several scales."""
    x = np.log(np.maximum(lum, 1e-7))
    x[disc] = np.nan
    out = np.zeros_like(x); n = 0
    fill = np.where(disc, np.nanmedian(x), x).astype(np.float32)
    for k in (4, 8, 16, 32, 64):
        m = blur(fill, k)
        d = fill - m
        sd = np.sqrt(np.maximum(blur(d*d, k), 1e-12))
        out += np.arctan(3.0*d/sd); n += 1
    return out/n

half = 900
y0, x0 = int(hcy)-half, int(hcx)-half
sl = (slice(max(y0,0), y0+2*half), slice(max(x0,0), x0+2*half))
imgs = {}
prof = {}
for nm in ("plain", "taper", "cont"):
    lum = acc[nm][0]/np.maximum(acc[nm][1], 1e-9)
    imgs[nm] = mgn(lum)[sl]
    pol = np.median(lum[np.newaxis], axis=0)
    prof[nm] = [float(np.median(pol[(rad >= f*hR-1) & (rad < f*hR+1)]))
                for f in (1.02, 1.06, 1.20, 1.50)]
for nm in ("plain", "taper", "cont"):
    print(nm, "  profile vs taper: " + " ".join(
        f"{a/b:6.3f}" for a, b in zip(prof[nm], prof["taper"])), flush=True)
d = imgs["taper"] - imgs["plain"]
print("difference: rms %.4f   p99 %.4f   max %.4f" % (
    float(np.sqrt((d*d).mean())), float(np.percentile(np.abs(d), 99)),
    float(np.abs(d).max())), flush=True)
dm = np.where(disc[sl], np.nan, d)
for lo, hi in ((1.02,1.10),(1.10,1.25),(1.25,1.50),(1.50,2.00)):
    m = (rad[sl] >= lo*hR) & (rad[sl] < hi*hR)
    print(f"  {lo:.2f}-{hi:.2f} R   rms diff {np.sqrt(np.nanmean(d[m]**2)):.4f}", flush=True)

def png(a, path, lo=-1.2, hi=1.2):
    u = np.clip((a-lo)/(hi-lo), 0, 1)
    cv2.imwrite(path, (u*255).astype(np.uint8))
os.makedirs(OUT, exist_ok=True)
png(imgs["taper"], os.path.join(OUT, "mgn_taper.png"))
png(imgs["plain"], os.path.join(OUT, "mgn_plain.png"))
png(d, os.path.join(OUT, "mgn_difference.png"), -0.5, 0.5)
png(imgs["cont"], os.path.join(OUT, "mgn_cont.png"))
dc2 = imgs["cont"] - imgs["plain"]
png(dc2, os.path.join(OUT, "mgn_cont_minus_plain.png"), -0.5, 0.5)
for lo, hi in ((1.02,1.10),(1.10,1.25),(1.25,1.50)):
    m = (rad[sl] >= lo*hR) & (rad[sl] < hi*hR)
    print(f"  cont-plain {lo:.2f}-{hi:.2f} R  rms {np.sqrt(np.nanmean(dc2[m]**2)):.4f}",
          flush=True)
print("wrote", OUT, flush=True)
