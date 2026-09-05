"""Full-resolution merge under one weight form, scaled onto the pipeline's own
hdr_lum units, so it can be fed to the real build_layers() and judged by the
real MGN instead of my stand-in."""
import json, os, re, sys
import numpy as np, cv2
cv2.setNumThreads(4)
F, MODE, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
def blur(a, s):
    k = int(2*round(3*s)+1)
    return cv2.GaussianBlur(a, (k, k), s, borderType=cv2.BORDER_REPLICATE)
def secs(n):
    t = re.match(r"tier_(.+?)_srgb\.tif$", n).group(1)
    return 1.0/float(t[2:].rstrip("s").replace("p",".")) if t.startswith("1_") \
        else float(t.rstrip("s").replace("p","."))
W = os.path.join(F, ".eclipseforgehdr")
geo = json.load(open(os.path.join(W, "geometry.json")))
R = float(geo["R"]); oy, ox = geo.get("crop_origin", [0, 0])
cy, cx = float(geo["cy"])-oy, float(geo["cx"])-ox
cal = {float(k): float(v) for k, v in geo["cal"].items()}
sigma = float(np.clip(0.032*R, 8.0, 40.0)); ALPHA = 0.55
tdir = os.path.join(F, "eclipseforge_output", "aligned_tiers")
files = sorted(f for f in os.listdir(tdir) if f.endswith("_srgb.tif"))
A = B = None
for fn in files:
    im = cv2.imread(os.path.join(tdir, fn), cv2.IMREAD_UNCHANGED)
    s = secs(fn); c = cal.get(s, 1.0)
    # channel by channel: never hold a 3-channel float32 copy
    cm = None; acc3 = None
    for ch in range(3):
        x = im[:, :, ch].astype(np.float32)/65535.0
        x = np.where(x <= 0.04045, x/12.92, ((x+0.055)/1.055)**2.4).astype(np.float32)
        cm = x.copy() if cm is None else np.maximum(cm, x)
        acc3 = x if acc3 is None else acc3+x
        del x
    del im
    lum = acc3/np.float32(3.0*s*c); del acc3
    if A is None:
        A = np.zeros(cm.shape, np.float32); B = np.zeros(cm.shape, np.float32)
    q = 0.5*(1.0+np.tanh((0.87-cm)/0.06)).astype(np.float32)
    q[cm > 0.97] = 0.0
    q *= 0.5*(1.0+np.tanh((cm-0.005)/0.0025))
    q = q.astype(np.float32)
    if MODE == "plain":
        w = np.float32(s**ALPHA)*blur(q, sigma)
    else:
        v = (cm <= 0.97).astype(np.float32)
        den = blur(v, sigma); num = blur(q*v, sigma)
        o = np.where(den > 1e-3, num/np.maximum(den, 1e-3), 0.0).astype(np.float32)
        t = np.clip((den-0.5)*2.0, 0.0, 1.0)
        w = np.float32(s**ALPHA)*o*t*t*v
        del v, den, num, o, t
    A += w*lum; B += w
    del cm, lum, q, w
    print("merged", fn, flush=True)
out = A/np.maximum(B, 1e-9); del A, B
ref = np.load(os.path.join(W, "hdr_lum.npy"), mmap_mode="r")
yy = np.arange(out.shape[0], dtype=np.float32)[:, None]-cy
xx = np.arange(out.shape[1], dtype=np.float32)[None, :]-cx
rr = np.sqrt(yy*yy+xx*xx); del yy, xx
band = (rr > 1.3*R) & (rr < 1.6*R)
k = float(np.median(np.asarray(ref)[band])/max(np.median(out[band]), 1e-30))
out *= np.float32(k)
print(f"scale onto hdr_lum units: {k:.6g}")
for f in (1.02, 1.10, 1.30, 1.60, 2.00):
    m = (rr >= f*R-2) & (rr < f*R+2)
    print(f"  {f:.2f}R  mine {np.median(out[m]):12.4f}   pipeline "
          f"{np.median(np.asarray(ref)[m]):12.4f}", flush=True)
np.save(OUT, out.astype(np.float32))
print("wrote", OUT, out.shape)
