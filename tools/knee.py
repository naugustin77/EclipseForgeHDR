"""Does WIDENING the roll-off in signal space kill the rings?

A tier currently stops contributing over a signal range of 0.06*sat, about a
fifth of a stop. Tiers differ by up to 28x in SNR at the same radius, so where
one drops out that fast the merged NOISE LEVEL steps -- and MGN divides by
local sigma, so a step in noise prints as a ring exactly like a step in
brightness. Widening the knee makes the hand-over gradual in SIGNAL space: no
spatial blur, so no leak, and the mask has nothing new to bite on.

Half resolution throughout (sigma scaled with it) so every variant fits in one
pass; the comparison between variants is what matters, not the absolute level.
"""
import json, os, re, sys
import numpy as np, cv2
cv2.setNumThreads(4)
F, TAG = sys.argv[1], sys.argv[2]
def srgb(x): return np.where(x <= 0.04045, x/12.92, ((x+0.055)/1.055)**2.4)
def blur(a, s):
    k = int(2*round(3*s)+1)
    return cv2.GaussianBlur(a, (k, k), s, borderType=cv2.BORDER_REPLICATE)
def secs(n):
    t = re.match(r"tier_(.+?)_srgb\.tif$", n).group(1)
    return 1.0/float(t[2:].rstrip("s").replace("p",".")) if t.startswith("1_") \
        else float(t.rstrip("s").replace("p","."))
geo = json.load(open(os.path.join(F, ".eclipseforgehdr", "geometry.json")))
R = float(geo["R"])/2.0; cy = (float(geo["cy"])-geo.get("crop_origin",[0,0])[0])/2.0
cx = (float(geo["cx"])-geo.get("crop_origin",[0,0])[1])/2.0
cal = {float(k): float(v) for k, v in geo["cal"].items()}
sigma = float(np.clip(0.032*R*2, 8.0, 40.0))/2.0
ALPHA = 0.55
tdir = os.path.join(F, "eclipseforge_output", "aligned_tiers")
files = sorted(f for f in os.listdir(tdir) if f.endswith("_srgb.tif"))
KNEES = {"cur": (0.87, 0.06)}
# Roll to EXACT zero, C1, at a level below the 0.97 hard cut, so the cut
# removes nothing and there is no step on the saturation contour at all.
# (The previous widening test was wrong: it widened the tanh but left the cut
# in place, which makes the step LARGER -- w40 still jumped 0.30 at 0.97.)
ZERO = {"z90": 0.90, "z80": 0.80, "z70": 0.70}
# Only spatial smoothing of the WEIGHT reduced the rings, and removing the
# hard cut did nothing. So the source is that w(x) = f(cm(x)) inherits every
# bit of the image's own fine structure: the tier MIXTURE is modulated by the
# corona's brightness pattern, and wherever two tiers disagree at all, that
# disagreement is printed along the isophotes. The cure is a mixture that
# varies smoothly in POSITION. Derive w from a smoothed intensity, and dilate
# first by the same scale so no pixel within it of a clipped one gets weight --
# that is what stops the leak, with no mask and no contour term.
SMOOTH = {"sm1": 1.0, "sm2": 2.0, "sm4": 4.0}
NAMES = ["none", "plain"] + list(KNEES) + ["z80_nf"] + list(SMOOTH)
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
    v = (cm <= 0.97).astype(np.float32)
    den = blur(v, sigma)
    t = np.clip((den-0.5)*2.0, 0.0, 1.0); t *= t
    floor = 0.5*(1.0+np.tanh((cm-0.005)/0.0025)).astype(np.float32)
    sa = np.float32(s**ALPHA)
    for n in NAMES:
        if n in SMOOTH:
            S = SMOOTH[n]*sigma
            rd = max(int(round(S)), 1)
            ker = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (2*rd+1,)*2)
            cs = blur(cv2.dilate(cm, ker), S)
            u = np.clip((0.90-cs)/0.45, 0.0, 1.0)
            w = sa*(u*u*(3-2*u)).astype(np.float32)*blur(floor, S)
            del cs, u
        elif n in ZERO or n == "z80_nf":
            hi = ZERO.get(n, 0.80)
            u = np.clip((hi-cm)/(0.5*hi), 0.0, 1.0)
            q = (u*u*(3-2*u)).astype(np.float32)*floor
            if n == "z80_nf":                     # no spatial feather at all
                w = sa*q
            else:
                nu = blur(q*v, sigma)
                w = sa*np.where(den > 1e-3, nu/np.maximum(den, 1e-3), 0.0
                                ).astype(np.float32)*t*v
                del nu
            del q, u
        elif n in KNEES:
            c0, wd = KNEES[n]
            q = 0.5*(1.0+np.tanh((c0-cm)/wd)).astype(np.float32)
            q[cm > 0.97] = 0.0; q *= floor
            nu = blur(q*v, sigma)
            w = sa*np.where(den > 1e-3, nu/np.maximum(den, 1e-3), 0.0
                            ).astype(np.float32)*t*v
            del q, nu
        else:
            q = 0.5*(1.0+np.tanh((0.87-cm)/0.06)).astype(np.float32)
            q[cm > 0.97] = 0.0; q *= floor
            w = sa*(q if n == "none" else blur(q, sigma))
            del q
        acc[n][0] += w*lum; acc[n][1] += w
        del w
    del cm, lum, v, den, t, floor
    print("merged", fn, flush=True)
yy = np.arange(acc["none"][0].shape[0], dtype=np.float32)[:, None]-cy
xx = np.arange(acc["none"][0].shape[1], dtype=np.float32)[None, :]-cx
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
# CONTROL. `plain` is the only variant with both a changed radial profile and
# a lower score. A purely RADIAL rescaling cannot create or destroy a ring that
# follows an isophote, so if applying plain's radial profile to the shipped
# merge reproduces plain's score, the metric is measuring the photometry, not
# the rings, and every number in this table is worthless.
lc = acc["cur"][0]/np.maximum(acc["cur"][1], 1e-9)
lp = acc["plain"][0]/np.maximum(acc["plain"][1], 1e-9)
rb = np.clip((rad/R-1.00)/0.006, 0, 264).astype(np.int32)
pc = np.bincount(rb.ravel(), lc.ravel(), 265)/np.maximum(np.bincount(rb.ravel(), None, 265), 1)
pp = np.bincount(rb.ravel(), lp.ravel(), 265)/np.maximum(np.bincount(rb.ravel(), None, 265), 1)
ratio = np.where(pc > 0, pp/np.maximum(pc, 1e-12), 1.0).astype(np.float32)
acc["ctrl"] = [lc*ratio[rb], np.ones_like(lc)]
NAMES = NAMES + ["ctrl"]
out = {}
for n in NAMES:
    l = acc[n][0]/np.maximum(acc[n][1], 1e-9)
    p = [float(np.median(l[(rad >= f*R-1) & (rad < f*R+1)])) for f in (1.02,1.06,1.20,1.50)]
    pol = cv2.remap(mgn(l), mx, my, cv2.INTER_LINEAR, borderMode=cv2.BORDER_REPLICATE)
    hp = pol - cv2.GaussianBlur(pol, (1, 21), 0, borderType=cv2.BORDER_REPLICATE)
    Fp = np.abs(np.fft.rfft(hp-hp.mean(axis=1, keepdims=True), axis=1))**2/NA
    out[n] = (p, float(Fp[:, 1:20].mean()))
bp, bring = out["cur"]
for n in NAMES:
    p, ring = out[n]
    print(f"{n:8s} " + " ".join(f"{x/y:7.3f}" for x, y in zip(p, bp)) +
          f" {ring:11.5f} {ring/bring:7.2f}x")
print("\nprofile relative to `cur` (the shipped weight); ring power likewise.")
