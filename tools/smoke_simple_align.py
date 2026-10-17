"""Simple stack: chained tier alignment and the merge noise map.

A synthetic bracket (8 tiers x 2 frames, x2.4 steps, so the long tiers clip
the inner corona and the short ones are noise in the outer field) whose
frames drift by a known shift through the bracket (0 -> 30 px), plus fine
streamers so there is something to align on. The merged picture is compared
with the truth shifted to the reference frame: an alignment error shows as a
doubled / blurred outer field (high-pass correlation with the truth).
"""
import os, sys, json, tempfile, shutil, re, importlib
import numpy as np
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, "..", "src")); sys.path.insert(0, _HERE)
from smoke_pipeline import _P as Prog
from smoke_channel_floors import write_cube
from scipy import ndimage

H, W, R = 520, 700, 60.0
NT, STEP, T0 = 8, 2.4, 1 / 500.0

def scene():
    cy, cx = H / 2, W / 2
    yy, xx = np.mgrid[:H, :W].astype(np.float64)
    r = np.hypot(yy - cy, xx - cx); th = np.arctan2(yy - cy, xx - cx)
    rng = np.random.default_rng(5)
    st = sum(rng.uniform(0.1, 0.35) * np.cos(k * th + rng.uniform(0, 6)) for k in (5, 9, 17, 29, 41, 63))
    base = 2.5e6 * (R / np.maximum(r, R)) ** 4.5 * (1 + 0.5 * st)
    base[r < R] = 0
    return base

def write(folder, base):
    os.makedirs(folder, exist_ok=True)
    gain = np.array([0.7, 0.88, 1.0])[:, None, None]
    rng = np.random.default_rng(11)
    n = 0; truth = {}
    for i in range(NT):
        t = T0 * STEP ** i
        for k in range(2):
            dy, dx = 30.0 * n / (2 * NT - 1), -15.0 * n / (2 * NT - 1)
            img = ndimage.shift(base, (dy, dx), order=1, mode="constant")
            a = img[None] * gain * t + 300.0 + 40.0 * t   # black + a little sky
            cube = a + rng.normal(0, np.sqrt(np.maximum(a, 1)) * 2.0)
            write_cube(os.path.join(folder, "f%02d.fits" % n), cube, t, "2026-08-12T18:29:%02d.0" % n)
            truth["f%02d.fits" % n] = (dy, dx); n += 1
    return truth

def run(mod_name, folder, mode):
    import eclipseforgehdr.simple as S
    if mod_name != "eclipseforgehdr.simple":
        spec = importlib.util.spec_from_file_location("eclipseforgehdr._simple_old", mod_name)
        S = importlib.util.module_from_spec(spec); S.__package__ = "eclipseforgehdr"; spec.loader.exec_module(S)
    P = Prog(False)
    S.run(folder, P, tier_mode=mode, partialconv=False)
    wd = os.path.join(folder, ".eclipseforgehdr")
    return P.lines, np.load(os.path.join(wd, "hdr_lum.npy")), wd

def score(lum, base):
    # the merge is in the first frame's coordinates (shift 0)
    hp = lambda a: np.log1p(np.maximum(a, 0)) - ndimage.gaussian_filter(np.log1p(np.maximum(a, 0)), 3)
    yy, xx = np.mgrid[:H, :W]; r = np.hypot(yy - H / 2, xx - W / 2) / R
    m = (r > 1.6) & (r < 3.0)
    a, b = hp(lum)[m], hp(base)[m]
    return float(np.corrcoef(a, b)[0, 1])

if __name__ == "__main__":
    old = sys.argv[1] if len(sys.argv) > 1 else None
    base = scene()
    tmp = tempfile.mkdtemp(prefix="efhdr_salign_")
    try:
        d = os.path.join(tmp, "b"); truth = write(d, base)
        for label, mod in (("old", old), ("new", "eclipseforgehdr.simple")):
            if mod is None: continue
            shutil.rmtree(os.path.join(d, ".eclipseforgehdr"), ignore_errors=True)
            lines, lum, wd = run(mod, d, "exposure")
            sh = [l for l in lines if "shift rms" in l]
            print(label, "outer-field detail correlation with truth: %.3f" % score(lum, base),
                  "| merge_noise.npy:", os.path.exists(os.path.join(wd, "merge_noise.npy")))
            for l in sh: print("   ", l.strip())
            fil = [l for l in lines if "filament detail" in l.lower()]
            print("   ", fil[-1] if fil else "")
            nm = [l for l in lines if "merge noise" in l.lower() or "noise floor" in l.lower()]
            for l in nm: print("   ", l)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
