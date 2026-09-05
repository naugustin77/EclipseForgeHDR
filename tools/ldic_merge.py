"""LDIC composition (Druckmullerova thesis eq. 4.15) against ours.

  g = SUM_i w(f_i) * ( k_i(phi) f_i + q_i(phi) )

k_i, q_i by linear regression against the composite accumulated so far, in 60
azimuth segments, longest exposure first, smoothed by a trigonometric
polynomial of order 4. Scored against the shipped weight with the same ring
metric (polar, high-pass in radius, power at m<20).
"""
import json, os, re, sys
import numpy as np, cv2
cv2.setNumThreads(4)
F = sys.argv[1]
FULL = len(sys.argv) > 2 and sys.argv[2] == "full"
OUT = sys.argv[3] if len(sys.argv) > 3 else None
NS, ORDER = 60, 4
def srgb(x): return np.where(x <= 0.04045, x/12.92, ((x+0.055)/1.055)**2.4)
def blur(a, s):
    k = int(2*round(3*s)+1)
    return cv2.GaussianBlur(a, (k, k), s, borderType=cv2.BORDER_REPLICATE)
def secs(n):
    t = re.match(r"tier_(.+?)_srgb\.tif$", n).group(1)
    return 1.0/float(t[2:].rstrip("s").replace("p",".")) if t.startswith("1_") \
        else float(t.rstrip("s").replace("p","."))
Wd = os.path.join(F, ".eclipseforgehdr")
geo = json.load(open(os.path.join(Wd, "geometry.json")))
oy, ox = geo.get("crop_origin", [0, 0])
D = 1.0 if FULL else 2.0
R = float(geo["R"])/D
cy, cx = (float(geo["cy"])-oy)/D, (float(geo["cx"])-ox)/D
cal = {float(k): float(v) for k, v in geo["cal"].items()}
sigma = float(np.clip(0.032*R*D, 8.0, 40.0))/D
ALPHA = 0.55
tdir = os.path.join(F, "eclipseforge_output", "aligned_tiers")
flist = sorted(((secs(f), f) for f in os.listdir(tdir) if f.endswith("_srgb.tif")),
               reverse=True)
def load(fn):
    im = cv2.imread(os.path.join(tdir, fn), cv2.IMREAD_UNCHANGED)
    a = srgb(im.astype(np.float32)/65535.0); del im
    if not FULL:
        H, W = a.shape[:2]
        a = cv2.resize(a, (W//2, H//2), interpolation=cv2.INTER_AREA)
    s = secs(fn)
    return a.max(axis=2), (a.mean(axis=2)/np.float32(s*cal.get(s, 1.0)))
cm, _ = load(flist[0][1]); H, W = cm.shape; del cm
yy = np.arange(H, dtype=np.float32)[:, None]-cy
xx = np.arange(W, dtype=np.float32)[None, :]-cx
rad = np.sqrt(yy*yy+xx*xx)
phi = np.arctan2(yy+0*xx, xx+0*yy).astype(np.float32)
seg = np.clip(np.floor((phi+np.pi)/(2*np.pi)*NS).astype(np.int32), 0, NS-1)
del yy, xx
def trig_fit(vals, fallback):
    c = (np.arange(NS)+0.5)/NS*2*np.pi - np.pi
    m = np.isfinite(vals)
    if m.sum() < 2*ORDER+3:
        return np.full(NS, fallback if not m.any() else np.nanmedian(vals))
    B = [np.ones(NS)]
    for o in range(1, ORDER+1):
        B += [np.cos(o*c), np.sin(o*c)]
    B = np.vstack(B).T
    sol, *_ = np.linalg.lstsq(B[m], vals[m], rcond=None)
    return B @ sol
G = np.zeros((H, W), np.float32); WT = np.zeros((H, W), np.float32)
Gc = np.zeros((H, W), np.float32); WTc = np.zeros((H, W), np.float32)
for s, fn in flist:
    cm, lum = load(fn)
    # ours, unchanged
    q0 = 0.5*(1.0+np.tanh((0.87-cm)/0.06)).astype(np.float32)
    q0[cm > 0.97] = 0.0
    q0 *= 0.5*(1.0+np.tanh((cm-0.005)/0.0025))
    v = (cm <= 0.97).astype(np.float32)
    den = blur(v, sigma); num = blur(q0.astype(np.float32)*v, sigma)
    o = np.where(den > 1e-3, num/np.maximum(den, 1e-3), 0.0).astype(np.float32)
    t = np.clip((den-0.5)*2.0, 0.0, 1.0)
    wc = np.float32(s**ALPHA)*o*t*t*v
    Gc += wc*lum; WTc += wc
    del q0, v, den, num, o, t, wc
    # LDIC
    u = np.clip((0.85-cm)/0.10, 0, 1)*np.clip((cm-0.004)/0.008, 0, 1)
    w = (u*u*(3-2*u)).astype(np.float32); del u
    if WT.max() > 0:
        comp = G/np.maximum(WT, 1e-9)
        ok = (WT > 0.05) & (w > 0.5) & (rad > 1.0*R) & (rad < 2.5*R) & (lum > 0)
        ks = np.full(NS, np.nan); qs = np.full(NS, np.nan)
        for j in range(NS):
            m = ok & (seg == j)
            if m.sum() < 500: continue
            x = lum[m].astype(np.float64); y = comp[m].astype(np.float64)
            sol, *_ = np.linalg.lstsq(np.vstack([x, np.ones_like(x)]).T, y, rcond=None)
            ks[j], qs[j] = sol
        kk = trig_fit(ks, 1.0).astype(np.float32)
        qq = trig_fit(qs, 0.0).astype(np.float32)
        adj = lum*kk[seg] + qq[seg]
        del comp, ok
    else:
        adj = lum
    G += w*adj; WT += w
    del cm, lum, w, adj
    print("composed", fn, flush=True)
gl = G/np.maximum(WT, 1e-9); gc = Gc/np.maximum(WTc, 1e-9)
if OUT:
    ref = np.load(os.path.join(Wd, "hdr_lum.npy"), mmap_mode="r")
    band = (rad > 1.3*R) & (rad < 1.6*R)
    k = float(np.median(np.asarray(ref)[band])/max(np.median(gl[band]), 1e-30))
    np.save(OUT, (gl*np.float32(k)).astype(np.float32))
    print("wrote", OUT, "scale", k)
disc = rad < 1.01*R
def mgn(l):
    x = np.log(np.maximum(l, 1e-7)).astype(np.float32)
    f = np.where(disc, np.median(x[~disc]), x).astype(np.float32)
    o = np.zeros_like(f)
    for kk_ in (2, 4, 8, 16, 32):
        d = f-blur(f, kk_)
        o += np.arctan(3.0*d/np.sqrt(np.maximum(blur(d*d, kk_), 1e-12)))
    return o/5.0
NA = 1024
rr = np.arange(1.02, 1.60, 0.004, dtype=np.float32)*R
th = np.arange(NA, dtype=np.float32)*(2*np.pi/NA)
mx = (cx + rr[:, None]*np.cos(th)[None, :]).astype(np.float32)
my = (cy + rr[:, None]*np.sin(th)[None, :]).astype(np.float32)
print(f"\n{'variant':8s} {'1.02R':>7s} {'1.06R':>7s} {'1.20R':>7s} {'1.50R':>7s}"
      f" {'ring power':>11s} {'vs ours':>8s}")
res = {}
for n, l in (("ours", gc), ("LDIC", gl)):
    p = [float(np.median(l[(rad >= f*R-1) & (rad < f*R+1)])) for f in (1.02,1.06,1.20,1.50)]
    pol = cv2.remap(mgn(l), mx, my, cv2.INTER_LINEAR, borderMode=cv2.BORDER_REPLICATE)
    hp = pol - cv2.GaussianBlur(pol, (1, 21), 0, borderType=cv2.BORDER_REPLICATE)
    Fp = np.abs(np.fft.rfft(hp-hp.mean(axis=1, keepdims=True), axis=1))**2/NA
    res[n] = (p, float(Fp[:, 1:20].mean()))
bp, br = res["ours"]
for n in ("ours", "LDIC"):
    p, r = res[n]
    print(f"{n:8s} " + " ".join(f"{x/y:7.3f}" for x, y in zip(p, bp)) +
          f" {r:11.5f} {r/br:7.2f}x")
