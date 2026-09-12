"""Does a per-tier RADIAL correction remove the azimuthal structure?

feather_bench showed the unfeathered merge carries it too, so the weight's
shape is not the source. The remaining candidate is the tiers disagreeing with
each other as a function of radius (TODO 1a), printing wherever the merge hands
over from one tier to the next. This fits k_i(r) from the tiers themselves and
re-scores.
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
log(f"R={R:.1f} sigma={sigma:.1f}  {len(files)} tiers")

def load(fn):
    im = cv2.imread(os.path.join(tdir, fn), cv2.IMREAD_UNCHANGED)
    a = srgb(im.astype(np.float32)/65535.0)
    cm = a.max(axis=2)
    s = secs(fn)
    return s, cm, (a.mean(axis=2)/np.float32(s*cal.get(s, 1.0)))

im0 = cv2.imread(os.path.join(tdir, files[0]), cv2.IMREAD_UNCHANGED)
H, W = im0.shape[:2]; del im0
yy = (np.arange(H, dtype=np.float32)[:, None] - (cy-oy))
xx = (np.arange(W, dtype=np.float32)[None, :] - (cx-ox))
rad = np.sqrt(yy*yy + xx*xx); del yy, xx

# ---- pass 1: each tier's radial profile, unclipped pixels only -----------
EDGES = np.arange(1.00, 2.60, 0.01, dtype=np.float32) * R
NB = len(EDGES) - 1
idx = (np.digitize(rad.ravel(), EDGES) - 1).astype(np.int32)
inb = (idx >= 0) & (idx < NB)
idxc = np.clip(idx, 0, NB - 1)
prof = {}
for fn in files:
    s_, cm, lum = load(fn)
    cmr = cm.ravel()
    tot = np.bincount(idxc[inb], minlength=NB).astype(np.float64)
    clp = np.bincount(idxc[inb & (cmr > 0.97)], minlength=NB).astype(np.float64)
    # A ring is usable for this tier only if the tier is unclipped essentially
    # all the way round it. Fitting from the unclipped pixels of a PARTLY
    # clipped ring samples only the faint azimuths, so the tier looks darker
    # than it is and k comes out large -- precisely on the tiers that matter.
    clean = (tot > 400) & (clp / np.maximum(tot, 1) < 0.02)
    ok = inb & (cmr >= 0.02) & clean[idxc]
    v = lum.ravel().astype(np.float64)
    cnt = np.bincount(idxc[ok], minlength=NB).astype(np.float64)
    sm = np.bincount(idxc[ok], weights=v[ok], minlength=NB)
    p = np.where(cnt > 400, sm / np.maximum(cnt, 1), np.nan)
    prof[s_] = p
    log(f"  profiled {fn}  ({int(np.isfinite(p).sum())}/{NB} clean rings)")
    del cm, cmr, lum, v, ok, cnt, sm

ref = np.nanmedian(np.stack([prof[s] for s in prof]), axis=0)
K = {}
for s, p in prof.items():
    k = np.where(np.isfinite(p) & np.isfinite(ref) & (p > 0), ref/np.maximum(p, 1e-12), 1.0)
    k = np.clip(k, 0.25, 4.0)
    kk = np.convolve(np.nan_to_num(k, nan=1.0), np.ones(9)/9.0, mode="same")
    kk[:4] = kk[4]; kk[-4:] = kk[-5]
    K[s] = kk.astype(np.float32)
    d = np.nanmax(np.abs(kk[:30] - 1.0))
    log(f"  tier {s:g}s  max |k-1| inside 1.3R = {100*d:.0f}%")

ib = idxc.reshape(H, W)
inside = (rad < 1.00*R) | (rad >= 2.59*R)

# ---- pass 2: merge with and without the correction -----------------------
HW = (W//2, H//2)
acc = {k: [np.zeros(HW[::-1], np.float32), np.zeros(HW[::-1], np.float32)]
       for k in ("taper", "taper_k", "none", "none_k")}
for fn in files:
    s, cm, lum = load(fn)
    ws = 0.5*(1.0+np.tanh((0.87-cm)/0.06)).astype(np.float32)
    ws[cm > 0.97] = 0.0
    ws *= 0.5*(1.0+np.tanh((cm-0.005)/0.0025))
    ws = ws.astype(np.float32)
    v = (cm <= 0.97).astype(np.float32)
    den = blur(v, sigma); num = blur(ws*v, sigma)
    o = np.where(den > 1e-3, num/np.maximum(den, 1e-3), 0.0).astype(np.float32)
    t = np.clip((den-0.5)*2.0, 0.0, 1.0)
    wt = o*t*t*v
    kmap = np.where(inside, np.float32(1.0), K[s][ib])
    lk = lum*kmap
    for nm, w, L in (("none", ws, lum), ("none_k", ws, lk),
                     ("taper", wt, lum), ("taper_k", wt, lk)):
        acc[nm][0] += cv2.resize(w*L, HW, interpolation=cv2.INTER_AREA)
        acc[nm][1] += cv2.resize(w,   HW, interpolation=cv2.INTER_AREA)
    del cm, lum, ws, v, den, num, o, t, wt, kmap, lk
    log(f"  merged {fn}")

# ---- score ---------------------------------------------------------------
hR = R/2.0; NA = 1024
rr = np.arange(0.95*hR, 2.30*hR, 0.5, dtype=np.float32)
th = np.arange(NA, dtype=np.float32)*(2*np.pi/NA)
mx = ((cx-ox)/2.0 + rr[:, None]*np.cos(th)[None, :]).astype(np.float32)
my = ((cy-oy)/2.0 + rr[:, None]*np.sin(th)[None, :]).astype(np.float32)
print(f"\n{'variant':9s} {'1.02R':>7s} {'1.06R':>7s} {'1.20R':>7s}"
      f" {'s1_m20_80':>10s} {'s1_m80_250':>11s} {'s2_m20_80':>10s}")
base = {}
for nm in ("none", "none_k", "taper", "taper_k"):
    num, den = acc[nm]
    lum = num/np.maximum(den, 1e-9)
    pol = cv2.remap(lum, mx, my, cv2.INTER_LINEAR, borderMode=cv2.BORDER_REPLICATE)
    med = np.median(pol, axis=1, keepdims=True)
    n = pol/np.maximum(med, 1e-12)
    P = np.abs(np.fft.rfft(n - n.mean(axis=1, keepdims=True), axis=1))**2/NA
    row = []
    for lo, hi, m0, m1 in ((1.05, 1.25, 20, 80), (1.05, 1.25, 80, 250),
                           (1.25, 1.50, 20, 80)):
        sel = (rr >= lo*hR) & (rr < hi*hR)
        row.append(float(P[sel, m0:m1].mean()))
    if nm == "none": base = list(row)
    pr = [float(med[int(np.argmin(np.abs(rr-f*hR))), 0]) for f in (1.02, 1.06, 1.20)]
    if nm == "none": ref_pr = list(pr)
    print(f"{nm:9s} " + " ".join(f"{a/b:7.3f}" for a, b in zip(pr, ref_pr)) +
          " " + " ".join(f"{a/b:10.2f}" for a, b in zip(row, base)))
print("\nprofile: level / uncorrected unfeathered.  azimuthal: power / uncorrected unfeathered.")
