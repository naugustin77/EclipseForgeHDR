"""Keep the smooth weight; remove the reason it hurts.

Only the plain blur reduces the rings (0.66x, and a radial-rescale control
proves that is not just its photometry). Its one distinguishing property is
that weight crosses each tier's clipping boundary -- which is exactly why it
was removed in 0.22.16, because a CLIPPED TIER UNDER-REPORTS and drags the
merge down.

Both can be had. The leak is only harmful because the clipped pixel carries a
wrong value. So: merge once conservatively to get an estimate, substitute that
estimate wherever a tier is clipped, then merge again with a smooth weight.
Nothing under-reports, and the weight never has to know where the contour is.
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

def load(fn):
    im = cv2.imread(os.path.join(tdir, fn), cv2.IMREAD_UNCHANGED)
    a = srgb(im.astype(np.float32)/65535.0); del im
    H, W = a.shape[:2]
    a = cv2.resize(a, (W//2, H//2), interpolation=cv2.INTER_AREA)
    cm = a.max(axis=2); s = secs(fn)
    return s, cm, a.mean(axis=2)/np.float32(s*cal.get(s, 1.0))

def wsat(cm):
    q = 0.5*(1.0+np.tanh((0.87-cm)/0.06)).astype(np.float32)
    q[cm > 0.97] = 0.0
    q *= 0.5*(1.0+np.tanh((cm-0.005)/0.0025))
    return q.astype(np.float32)

# ---- pass 1: the shipped conservative merge, to get an estimate ----------
A = B = None
for fn in files:
    s, cm, lum = load(fn)
    if A is None:
        A = np.zeros(cm.shape, np.float32); B = np.zeros(cm.shape, np.float32)
    q = wsat(cm); v = (cm <= 0.97).astype(np.float32)
    den = blur(v, sigma); num = blur(q*v, sigma)
    o = np.where(den > 1e-3, num/np.maximum(den, 1e-3), 0.0).astype(np.float32)
    t = np.clip((den-0.5)*2.0, 0.0, 1.0)
    w = np.float32(s**ALPHA)*o*t*t*v
    A += w*lum; B += w
    del cm, lum, q, v, den, num, o, t, w
est = A/np.maximum(B, 1e-9)
print("pass 1 done", flush=True)

# ---- pass 2 ---------------------------------------------------------------
# MGN divides by the LOCAL standard deviation over windows of 4-64 px. In the
# inner corona the noise is tiny, so a 1% systematic step at a tier hand-over
# is many sigma there -- a visible ring. And the current hand-over happens over
# about one feather sigma, ~20 px, which sits right in the middle of MGN's
# sensitive band. So the transition has to be made WIDE compared with 64 px,
# not merely smooth. Filling the clipped pixels first is what lets the feather
# be widened without the leak that widening would otherwise cause.
MULT = (1.0, 3.0, 6.0)
NAMES = ["cur", "plain"] + [f"fill{m:g}" for m in MULT]
acc = {n: [np.zeros(est.shape, np.float32), np.zeros(est.shape, np.float32)]
       for n in NAMES}
for fn in files:
    s, cm, lum = load(fn)
    q = wsat(cm); v = (cm <= 0.97).astype(np.float32)
    den = blur(v, sigma); num = blur(q*v, sigma)
    o = np.where(den > 1e-3, num/np.maximum(den, 1e-3), 0.0).astype(np.float32)
    t = np.clip((den-0.5)*2.0, 0.0, 1.0)
    sa = np.float32(s**ALPHA)
    # substitute the estimate wherever this tier is at or near its ceiling
    g = np.clip((cm-0.85)/0.10, 0.0, 1.0).astype(np.float32)
    filled = lum*(1.0-g) + est*g
    # a weight that no longer needs to know where the contour is
    qs = 0.5*(1.0+np.tanh((0.87-cm)/0.06)).astype(np.float32)
    qs *= 0.5*(1.0+np.tanh((cm-0.005)/0.0025))
    qs = qs.astype(np.float32)                 # NO hard cut, NO mask
    acc["cur"][0] += (sa*o*t*t*v)*lum; acc["cur"][1] += sa*o*t*t*v
    wp = sa*blur(q, sigma); acc["plain"][0] += wp*lum; acc["plain"][1] += wp
    del wp
    for m in MULT:
        w = sa*blur(qs, m*sigma)
        acc[f"fill{m:g}"][0] += w*filled; acc[f"fill{m:g}"][1] += w
        del w
    del cm, lum, q, v, den, num, o, t, g, filled, qs
print("pass 2 done", flush=True)

yy = np.arange(est.shape[0], dtype=np.float32)[:, None]-cy
xx = np.arange(est.shape[1], dtype=np.float32)[None, :]-cx
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
print(f"\n{'variant':11s} {'1.02R':>7s} {'1.06R':>7s} {'1.20R':>7s} {'1.50R':>7s}"
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
    print(f"{n:11s} " + " ".join(f"{x/y:7.3f}" for x, y in zip(p, bp)) +
          f" {r:11.5f} {r/br:7.2f}x")
