"""The merge weight has no SNR term. Test whether adding one fixes it.

single.py showed the azimuthal power at 1.05-1.25 R scales with how SHORT the
tier is -- 28x for 1/1000 s against 1/60 s -- which is the signature of shot
noise, not of corona structure. wsat only asks whether a pixel is UNCLIPPED; a
1/1000 s frame at 1.15 R is unclipped and nearly signal-free, so it gets full
weight and its noise goes straight into the merge.

Inverse-variance: the radiance estimate is x/(s*cal) with x the signal in units
of saturation, so its variance goes as (x + xr)/(s*cal)^2 -- shot noise in x
plus a read floor. Weighting by the inverse of that is the textbook estimator
(Mann & Picard 1995; Robertson 1999; Granados 2010) and is what s**alpha is a
crude stand-in for.
"""
import json, os, re, sys, time
import numpy as np, cv2
cv2.setNumThreads(4)
F = sys.argv[1]
def log(*a): print(f"[{time.strftime('%H:%M:%S')}]", *a, flush=True)
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
log(f"R={R:.1f} sigma={sigma:.1f} {len(files)} tiers")

VAR = ("cur_a1", "cur_a055", "floor02", "iv_3e4", "iv_1e3")
HW = None; acc = {}
for fn in files:
    im = cv2.imread(os.path.join(tdir, fn), cv2.IMREAD_UNCHANGED)
    a = srgb(im.astype(np.float32)/65535.0); del im
    cm = a.max(axis=2); s = secs(fn); c = cal.get(s, 1.0)
    lum = a.mean(axis=2)/np.float32(s*c); del a
    if HW is None:
        H, W = cm.shape; HW = (W//2, H//2)
    base = 0.5*(1.0+np.tanh((0.87-cm)/0.06)).astype(np.float32)
    base[cm > 0.97] = 0.0
    v = (cm <= 0.97).astype(np.float32)
    den = blur(v, sigma); num = blur(base*v, sigma)
    o = np.where(den > 1e-3, num/np.maximum(den, 1e-3), 0.0).astype(np.float32)
    t = np.clip((den-0.5)*2.0, 0.0, 1.0)
    fw = (o*t*t*v).astype(np.float32)          # the shipped taper, SNR-free
    def lowend(lo):
        return 0.5*(1.0+np.tanh((cm-lo)/(0.5*lo))).astype(np.float32)
    W_ = {}
    W_["cur_a1"]   = (s**1.0) * fw * lowend(0.005)
    W_["cur_a055"] = (s**0.55) * fw * lowend(0.005)
    W_["floor02"]  = (s**1.0) * fw * lowend(0.02)
    for nm, xr in (("iv_3e4", 3e-4), ("iv_1e3", 1e-3)):
        W_[nm] = ((s*c)**2 * fw / (cm + np.float32(xr))).astype(np.float32)
    for nm in VAR:
        w = W_[nm]
        if nm not in acc:
            acc[nm] = [np.zeros(HW[::-1], np.float32), np.zeros(HW[::-1], np.float32)]
        acc[nm][0] += cv2.resize(w*lum, HW, interpolation=cv2.INTER_AREA)
        acc[nm][1] += cv2.resize(w, HW, interpolation=cv2.INTER_AREA)
    del cm, lum, base, v, den, num, o, t, fw, W_
    log(f"  merged {fn}")

hR = R/2.0; NA = 1024
rr = np.arange(0.95*hR, 2.30*hR, 0.5, dtype=np.float32)
th = np.arange(NA, dtype=np.float32)*(2*np.pi/NA)
mx = ((cx-ox)/2.0 + rr[:, None]*np.cos(th)[None, :]).astype(np.float32)
my = ((cy-oy)/2.0 + rr[:, None]*np.sin(th)[None, :]).astype(np.float32)
print(f"\n{'variant':10s} {'1.02R':>7s} {'1.06R':>7s} {'1.20R':>7s} {'1.50R':>7s} {'2.00R':>7s}"
      f" {'s1_m20_80':>10s} {'s1_m80_250':>11s} {'s2_m20_80':>10s} {'s3_m20_80':>10s}")
b = None; rp = None
for nm in VAR:
    num, den = acc[nm]
    lum = num/np.maximum(den, 1e-9)
    pol = cv2.remap(lum, mx, my, cv2.INTER_LINEAR, borderMode=cv2.BORDER_REPLICATE)
    med = np.median(pol, axis=1, keepdims=True)
    n = pol/np.maximum(med, 1e-12)
    P = np.abs(np.fft.rfft(n - n.mean(axis=1, keepdims=True), axis=1))**2/NA
    row = []
    for lo, hi, m0, m1 in ((1.05,1.25,20,80),(1.05,1.25,80,250),
                           (1.25,1.50,20,80),(1.50,2.00,20,80)):
        sel = (rr >= lo*hR) & (rr < hi*hR)
        row.append(float(P[sel, m0:m1].mean()))
    pr = [float(med[int(np.argmin(np.abs(rr-f*hR))), 0]) for f in (1.02,1.06,1.20,1.50,2.00)]
    if b is None: b, rp = list(row), list(pr)
    print(f"{nm:10s} " + " ".join(f"{x/y:7.3f}" for x, y in zip(pr, rp)) +
          " " + " ".join(f"{x/y:10.2f}" for x, y in zip(row, b)))
print("\nall columns relative to cur_a1 (the shipped 0.22.24 weight).")
print("profile 1.000 = same photometry; azimuthal <1 = less fine-scale variance.")
