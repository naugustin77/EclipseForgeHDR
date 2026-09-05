"""Control: score INDIVIDUAL tiers with the same azimuthal metric.

If one unmerged tier scores like the merge, the metric is measuring the corona
itself and the merge adds nothing -- and the plain blur's low score means only
that a blur removes fine detail, not that it removes an artifact.
"""
import json, os, re, sys, time
import numpy as np, cv2
cv2.setNumThreads(4)
F = sys.argv[1]
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
hR = R/2.0; NA = 1024
rr = np.arange(0.95*hR, 2.30*hR, 0.5, dtype=np.float32)
th = np.arange(NA, dtype=np.float32)*(2*np.pi/NA)
mx = ((cx-ox)/2.0 + rr[:, None]*np.cos(th)[None, :]).astype(np.float32)
my = ((cy-oy)/2.0 + rr[:, None]*np.sin(th)[None, :]).astype(np.float32)
def score(lum):
    pol = cv2.remap(lum, mx, my, cv2.INTER_LINEAR, borderMode=cv2.BORDER_REPLICATE)
    med = np.median(pol, axis=1, keepdims=True)
    n = pol/np.maximum(med, 1e-12)
    P = np.abs(np.fft.rfft(n - n.mean(axis=1, keepdims=True), axis=1))**2/NA
    out = []
    for lo, hi, m0, m1 in ((1.05,1.25,20,80),(1.05,1.25,80,250),(1.25,1.50,20,80)):
        sel = (rr >= lo*hR) & (rr < hi*hR)
        out.append(float(P[sel, m0:m1].mean()))
    return out
print(f"{'tier':10s} {'clip% 1.05-1.25R':>16s} {'s1_m20_80':>10s} {'s1_m80_250':>11s} {'s2_m20_80':>10s}")
rows = {}
for fn in files:
    im = cv2.imread(os.path.join(tdir, fn), cv2.IMREAD_UNCHANGED)
    a = srgb(im.astype(np.float32)/65535.0); del im
    cm = a.max(axis=2); s = secs(fn)
    lum = a.mean(axis=2)/np.float32(s*cal.get(s,1.0)); del a
    H, W = cm.shape
    yy = np.arange(H, dtype=np.float32)[:,None]-(cy-oy)
    xx = np.arange(W, dtype=np.float32)[None,:]-(cx-ox)
    rad = np.sqrt(yy*yy+xx*xx); del yy, xx
    band = (rad >= 1.05*R) & (rad < 1.25*R)
    cf = float((cm[band] > 0.97).mean())*100
    half = cv2.resize(lum, (W//2, H//2), interpolation=cv2.INTER_AREA)
    rows[fn] = (cf, score(half))
    del cm, lum, rad, band, half
    print(f"{fn[5:-10]:10s} {cf:15.1f}% " +
          " ".join(f"{v:10.4g}" for v in rows[fn][1]), flush=True)
clean = [f for f in files if rows[f][0] < 1.0]
if clean:
    b = rows[clean[-1]][1]
    print(f"\nrelative to the longest tier that is unclipped in 1.05-1.25R "
          f"({clean[-1][5:-10]}):")
    for f in files:
        print(f"  {f[5:-10]:10s} " + " ".join(f"{v/bb:8.2f}" for v, bb in zip(rows[f][1], b)))
