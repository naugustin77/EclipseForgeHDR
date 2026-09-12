"""Hard-exclude each tier's contaminated collar.

The collar measurement: adjacent tiers disagree by up to 3x within 0-8 px
outside the longer tier's saturated region, and by ~1% beyond 8-16 px. The
shipped taper does not exclude that collar, it only fades it -- den just
outside a large saturated area is 0.655 at 8 px, so t*t = 0.096. Ten percent of
full weight on a 3x error is a 30% contribution, printed along the contour.

So: v = 1 - dilate(saturated, Rc). Exactly zero throughout the collar, for the
tier that is saturated, at several radii.
"""
import json, os, re, sys
import numpy as np, cv2
cv2.setNumThreads(4)
F = sys.argv[1]
def srgb(x): return np.where(x <= 0.04045, x/12.92, ((x+0.055)/1.055)**2.4)
def blur(a, s):
    k = int(2*round(3*s)+1)
    return cv2.GaussianBlur(a, (k, k), s, borderType=cv2.BORDER_REPLICATE)
def secs(n):
    t = re.match(r"tier_(.+?)_srgb\.tif$", n).group(1)
    return 1.0/float(t[2:].rstrip("s").replace("p",".")) if t.startswith("1_") \
        else float(t.rstrip("s").replace("p","."))
geo = json.load(open(os.path.join(F, ".eclipseforgehdr", "geometry.json")))
oy, ox = geo.get("crop_origin", [0, 0])
R = float(geo["R"])/2.0
cy, cx = (float(geo["cy"])-oy)/2.0, (float(geo["cx"])-ox)/2.0
cal = {float(k): float(v) for k, v in geo["cal"].items()}
sigma = float(np.clip(0.032*R*2, 8.0, 40.0))/2.0
ALPHA = 0.55
tdir = os.path.join(F, "eclipseforge_output", "aligned_tiers")
files = sorted(f for f in os.listdir(tdir) if f.endswith("_srgb.tif"))
# collar radii in HALF-res px; double for full res
RC = (0, 2, 4, 8, 16)
NAMES = ["cur", "plain"] + [f"c{r}" for r in RC]
acc = None
for fn in files:
    im = cv2.imread(os.path.join(tdir, fn), cv2.IMREAD_UNCHANGED)
    a = srgb(im.astype(np.float32)/65535.0); del im
    H, W = a.shape[:2]
    a = cv2.resize(a, (W//2, H//2), interpolation=cv2.INTER_AREA)
    cm = a.max(axis=2); s = secs(fn)
    lum = a.mean(axis=2)/np.float32(s*cal.get(s, 1.0)); del a
    if acc is None:
        acc = {n: [np.zeros(cm.shape, np.float32), np.zeros(cm.shape, np.float32)]
               for n in NAMES}
    q = 0.5*(1.0+np.tanh((0.87-cm)/0.06)).astype(np.float32)
    q[cm > 0.97] = 0.0
    q *= 0.5*(1.0+np.tanh((cm-0.005)/0.0025))
    q = q.astype(np.float32)
    sat = (cm > 0.97).astype(np.uint8)
    sa = np.float32(s**ALPHA)
    for r in RC:
        if r:
            k = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (2*r+1,)*2)
            v = 1.0 - cv2.dilate(sat, k).astype(np.float32)
        else:
            v = 1.0 - sat.astype(np.float32)
        den = blur(v, sigma); num = blur(q*v, sigma)
        o = np.where(den > 1e-3, num/np.maximum(den, 1e-3), 0.0).astype(np.float32)
        t = np.clip((den-0.5)*2.0, 0.0, 1.0)
        w = sa*o*t*t*v
        n = f"c{r}"
        acc[n][0] += w*lum; acc[n][1] += w
        if r == 0:
            acc["cur"][0] += w*lum; acc["cur"][1] += w
        del v, den, num, o, t, w
    wp = sa*blur(q, sigma)
    acc["plain"][0] += wp*lum; acc["plain"][1] += wp
    del cm, lum, q, sat, wp
    print("merged", fn, flush=True)
yy = np.arange(acc["cur"][0].shape[0], dtype=np.float32)[:, None]-cy
xx = np.arange(acc["cur"][0].shape[1], dtype=np.float32)[None, :]-cx
rad = np.sqrt(yy*yy+xx*xx); disc = rad < 1.01*R; del yy, xx
def mgn(l):
    x = np.log(np.maximum(l, 1e-7)).astype(np.float32)
    f = np.where(disc, np.median(x[~disc]), x).astype(np.float32)
    o = np.zeros_like(f)
    for k in (2, 4, 8, 16, 32):
        d = f-blur(f, k)
        o += np.arctan(3.0*d/np.sqrt(np.maximum(blur(d*d, k), 1e-12)))
    return o/5.0
NA = 1024
rr = np.arange(1.02, 1.60, 0.004, dtype=np.float32)*R
th = np.arange(NA, dtype=np.float32)*(2*np.pi/NA)
mx = (cx + rr[:, None]*np.cos(th)[None, :]).astype(np.float32)
my = (cy + rr[:, None]*np.sin(th)[None, :]).astype(np.float32)
print(f"\n{'variant':8s} {'1.02R':>7s} {'1.06R':>7s} {'1.20R':>7s} {'1.50R':>7s}"
      f" {'ring power':>11s} {'vs cur':>8s}")
res = {}
for n in NAMES:
    l = acc[n][0]/np.maximum(acc[n][1], 1e-9)
    p = [float(np.median(l[(rad >= f*R-1) & (rad < f*R+1)])) for f in (1.02,1.06,1.20,1.50)]
    pol = cv2.remap(mgn(l), mx, my, cv2.INTER_LINEAR, borderMode=cv2.BORDER_REPLICATE)
    hp = pol - cv2.GaussianBlur(pol, (1, 21), 0, borderType=cv2.BORDER_REPLICATE)
    Fp = np.abs(np.fft.rfft(hp-hp.mean(axis=1, keepdims=True), axis=1))**2/NA
    res[n] = (p, float(Fp[:, 1:20].mean()))
bp, br = res["cur"]
for n in NAMES:
    p, r = res[n]
    print(f"{n:8s} " + " ".join(f"{x/y:7.3f}" for x, y in zip(p, bp)) +
          f" {r:11.5f} {r/br:7.2f}x")
print("\ncN = weight exactly zero within N half-res px outside the saturated region.")
