"""Which TERM of the feather prints the rings? Difference images, not metrics.

taper = blur(w*v)/blur(v) * t^2 * v. Three separable suspects:
  1/den  -- renormalisation. Boosts the weight by up to 2x in a sigma-wide band
            just INSIDE each tier's clipping contour, where wsat is nowhere near
            zero because the corona gradient is steep. A 2x weight modulation
            shaped like the contour is exactly a ring.
  t^2    -- the taper ramp, also keyed to the contour.
  v      -- the hard mask.
plain has none of them and shows no rings; `none` has the hard wsat step only.
"""
import json, os, re, sys
import numpy as np, cv2
cv2.setNumThreads(4)
F, WANT, TAG = sys.argv[1], sys.argv[2].split(","), sys.argv[3]
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
sigma = float(np.clip(0.032*R, 8.0, 40.0)); ALPHA = 0.55
tdir = os.path.join(F, "eclipseforge_output", "aligned_tiers")
files = sorted(f for f in os.listdir(tdir) if f.endswith("_srgb.tif"))
HW = None; acc = {n: None for n in WANT}
for fn in files:
    im = cv2.imread(os.path.join(tdir, fn), cv2.IMREAD_UNCHANGED)
    a = srgb(im.astype(np.float32)/65535.0); del im
    cm = a.max(axis=2); s = secs(fn)
    lum = a.mean(axis=2)/np.float32(s*cal.get(s, 1.0)); del a
    if HW is None:
        H, W = cm.shape; HW = (W//2, H//2)
        for n in acc: acc[n] = [np.zeros(HW[::-1], np.float32),
                                np.zeros(HW[::-1], np.float32)]
    ws = 0.5*(1.0+np.tanh((0.87-cm)/0.06)).astype(np.float32)
    ws[cm > 0.97] = 0.0
    ws *= 0.5*(1.0+np.tanh((cm-0.005)/0.0025))
    ws = ws.astype(np.float32)
    v = (cm <= 0.97).astype(np.float32)
    den = blur(v, sigma); num = blur(ws*v, sigma)
    o = np.where(den > 1e-3, num/np.maximum(den, 1e-3), 0.0).astype(np.float32)
    t = np.clip((den-0.5)*2.0, 0.0, 1.0)
    # a weight already zero at 0.90 sat, C1, so the mask has nothing to cut
    u = np.clip((0.90-cm)/0.10, 0.0, 1.0)
    wc = (u*u*(3-2*u)).astype(np.float32)
    wc *= 0.5*(1.0+np.tanh((cm-0.005)/0.0025))
    # A tier's contribution currently ends over a signal range of 0.06*sat --
    # about a fifth of a stop. Tiers differ hugely in SNR at the same radius
    # (28x measured), so where one drops out that fast, the merged NOISE LEVEL
    # steps. MGN divides by local sigma, so a step in noise prints as a ring
    # just as readily as a step in brightness. Widening the roll-off makes the
    # hand-over gradual in SIGNAL space -- no spatial blur, so no leak and no
    # mask needed for it.
    def knee(c0, wd):
        q = 0.5*(1.0+np.tanh((c0-cm)/wd)).astype(np.float32)
        q[cm > 0.97] = 0.0
        q *= 0.5*(1.0+np.tanh((cm-0.005)/0.0025))
        q = q.astype(np.float32)
        nu = blur(q*v, sigma)          # den and t are already computed above
        oo = np.where(den > 1e-3, nu/np.maximum(den, 1e-3), 0.0).astype(np.float32)
        return oo*t*t*v
    W_ = {"none":     ws,
          "wide15":   knee(0.87, 0.15),
          "wide25":   knee(0.85, 0.25),
          "wide40":   knee(0.80, 0.40),
          "plain":    blur(ws, sigma),
          "taper":    o*t*t*v,
          "nonorm":   num*v,                 # masked blur, NO renormalisation
          "nonorm_t": num*t*t*v,             # taper ramp, NO renormalisation
          "cont_pl":  blur(wc, sigma),        # C1 roll-off + plain blur, no mask
          # Inverse variance. The radiance estimate is x/(s*cal), so its
          # variance goes as (x + x_r)/(s*cal)^2 -- shot noise plus a read
          # floor. s**alpha is a scene-blind stand-in for that; at alpha 0.55
          # it hands a nearly signal-free short tier almost as much weight as a
          # well exposed long one, and where that short tier's contribution
          # starts and stops -- its saturation contour -- the NOISE LEVEL of
          # the merge changes abruptly. MGN prints a change in noise level as
          # a ring just as readily as a change in brightness.
          "iv":       (o*t*t*v*np.float32((s*cal.get(s,1.0))**2)
                       / (cm + np.float32(3e-4)) / np.float32(s**ALPHA)),
          "iv_plain": (blur(ws, sigma)*np.float32((s*cal.get(s,1.0))**2)
                       / (cm + np.float32(3e-4)) / np.float32(s**ALPHA))}
    for n in WANT:
        w = np.float32(s**ALPHA)*W_[n]
        acc[n][0] += cv2.resize(w*lum, HW, interpolation=cv2.INTER_AREA)
        acc[n][1] += cv2.resize(w,   HW, interpolation=cv2.INTER_AREA)
    del cm, lum, ws, v, den, num, o, t, u, wc, W_
    print("merged", fn, flush=True)
hcy, hcx, hR = (cy-oy)/2.0, (cx-ox)/2.0, R/2.0
Hh, Wh = acc[WANT[0]][0].shape
yy = np.arange(Hh, dtype=np.float32)[:, None]-hcy
xx = np.arange(Wh, dtype=np.float32)[None, :]-hcx
radh = np.sqrt(yy*yy+xx*xx); disc = radh < 1.01*hR
def mgn(l):
    x = np.log(np.maximum(l, 1e-7)).astype(np.float32)
    f = np.where(disc, np.median(x[~disc]), x).astype(np.float32)
    o = np.zeros_like(f)
    for k in (4, 8, 16, 32, 64):
        d = f-blur(f, k)
        o += np.arctan(3.0*d/np.sqrt(np.maximum(blur(d*d, k), 1e-12)))
    return o/5.0
half = 1100
y0, x0 = max(int(hcy)-half, 0), max(int(hcx)-half, 0)
sl = (slice(y0, y0+2*half), slice(x0, x0+2*half))
for n in WANT:
    l = acc[n][0]/np.maximum(acc[n][1], 1e-9)
    np.save(f"{TAG}_{n}_mgn.npy", mgn(l)[sl])
    np.save(f"{TAG}_{n}_prof.npy", np.array(
        [np.median(l[(radh >= f*hR-1) & (radh < f*hR+1)]) for f in
         (1.02, 1.06, 1.20, 1.50)]))
np.save(f"{TAG}_rad.npy", radh[sl]/hR)
print("saved", WANT, flush=True)
