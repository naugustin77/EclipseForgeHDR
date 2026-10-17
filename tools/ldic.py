"""Measure the term LDIC has and we do not: a per-tier AFFINE transform that
varies with AZIMUTH.

Druckmullerova, Doctoral thesis, eq. 4.15:
    g(r,phi) = SUM_i  w(f_i(r,phi)) * ( k_i(phi)*f_i(r,phi) + q_i(phi) )
k and q are fitted by linear regression in each of 60 angular segments, against
the composite accumulated SO FAR, starting from the longest exposure, then
smoothed with a low-order trigonometric polynomial. The thesis says plainly
what they are for: "to compose images with different distribution of diffuse
light in the optical system ... or even images that were taken through thin
clouds."

EclipseForgeHDR uses ONE SCALAR per tier (cal[s]) and one shared additive
pedestal. If k_i varies with azimuth, that variation is a real error we are not
correcting, and it changes wherever the tier mix changes -- which is along the
saturation contours, exactly where the rings are.

This does not build the composer. It fits k_i(phi), q_i(phi) against the
running composite and reports how far from (1, 0) they are.
"""
import json, os, re, sys
import numpy as np, cv2
cv2.setNumThreads(4)
F = sys.argv[1]
NS = 60
def srgb(x): return np.where(x <= 0.04045, x/12.92, ((x+0.055)/1.055)**2.4)
def secs(n):
    t = re.match(r"tier_(.+?)_srgb\.tif$", n).group(1)
    return 1.0/float(t[2:].rstrip("s").replace("p",".")) if t.startswith("1_") \
        else float(t.rstrip("s").replace("p","."))
geo = json.load(open(os.path.join(F, ".eclipseforgehdr", "geometry.json")))
oy, ox = geo.get("crop_origin", [0, 0])
R = float(geo["R"])/2.0
cy, cx = (float(geo["cy"])-oy)/2.0, (float(geo["cx"])-ox)/2.0
cal = {float(k): float(v) for k, v in geo["cal"].items()}
tdir = os.path.join(F, "eclipseforge_output", "aligned_tiers")
files = sorted(((secs(f), f) for f in os.listdir(tdir) if f.endswith("_srgb.tif")),
               reverse=True)                      # longest exposure first
def load(fn):
    im = cv2.imread(os.path.join(tdir, fn), cv2.IMREAD_UNCHANGED)
    a = srgb(im.astype(np.float32)/65535.0); del im
    H, W = a.shape[:2]
    a = cv2.resize(a, (W//2, H//2), interpolation=cv2.INTER_AREA)
    s = secs(fn)
    return a.max(axis=2), a.mean(axis=2)/np.float32(s*cal.get(s, 1.0))
cm0, _ = load(files[0][1])
H, W = cm0.shape; del cm0
yy = np.arange(H, dtype=np.float32)[:, None]-cy
xx = np.arange(W, dtype=np.float32)[None, :]-cx
rad = np.sqrt(yy*yy+xx*xx)
seg = np.floor((np.arctan2(yy+0*xx, xx+0*yy)+np.pi)/(2*np.pi)*NS).astype(np.int32)
seg = np.clip(seg, 0, NS-1); del yy, xx
G = np.zeros((H, W), np.float32); Wt = np.zeros((H, W), np.float32)
print(f"{'tier':10s} {'segs':>5s} {'k median':>9s} {'k spread(%)':>12s} "
      f"{'q/median signal(%)':>19s}")
for s, fn in files:
    cm, lum = load(fn)
    w = np.clip((0.85-cm)/0.10, 0, 1)*np.clip((cm-0.004)/0.008, 0, 1)
    w = (w*w*(3-2*w)).astype(np.float32)          # LDIC's w: reject top and bottom
    if Wt.max() > 0:
        comp = G/np.maximum(Wt, 1e-9)
        ok = (Wt > 0.05) & (w > 0.5) & (rad > 1.0*R) & (rad < 2.5*R) & (lum > 0)
        ks, qs, n = [], [], 0
        for j in range(NS):
            m = ok & (seg == j)
            if m.sum() < 500:
                ks.append(np.nan); qs.append(np.nan); continue
            x = lum[m].astype(np.float64); y = comp[m].astype(np.float64)
            A = np.vstack([x, np.ones_like(x)]).T
            sol, *_ = np.linalg.lstsq(A, y, rcond=None)
            ks.append(sol[0]); qs.append(sol[1]); n += 1
        ks = np.array(ks, float); qs = np.array(qs, float)
        if n >= 10:
            km = np.nanmedian(ks)
            spread = 100*(np.nanpercentile(ks, 90)-np.nanpercentile(ks, 10))/max(km, 1e-9)
            msig = np.nanmedian(comp[ok])
            qsp = 100*np.nanmedian(np.abs(qs))/max(msig, 1e-9)
            print(f"{s:<10.5g} {n:5d} {km:9.4f} {spread:12.1f} {qsp:19.2f}", flush=True)
    G += w*lum; Wt += w
    del cm, lum, w
print("\nk spread = 10th-to-90th percentile of k across the 60 azimuth segments,")
print("as a percent of that tier's median k. A single scalar cal[s] can only")
print("express the median; the spread is the part we are not correcting.")
