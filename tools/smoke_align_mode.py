"""Alignment setting (auto / align / keep) on both paths.

A synthetic bracket written twice: once registered by a 'writer' (black
borders of different widths, frames already aligned), once drifting with no
borders. Checks what each setting does on each."""
import os, sys, tempfile, shutil, json
import numpy as np
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, "..", "src")); sys.path.insert(0, _HERE)
from smoke_pipeline import _P as Prog
from smoke_channel_floors import write_cube
import smoke_simple_align as T
from eclipseforgehdr import simple as S, pipeline as PL

def make(folder, prereg):
    base = T.scene(); gain = np.array([0.7, 0.88, 1.0])[:, None, None]
    rng = np.random.default_rng(3); n = 0
    from scipy import ndimage
    for i in range(T.NT):
        t = T.T0 * T.STEP ** i
        for k in range(2):
            img = base if prereg else ndimage.shift(base, (30.0 * n / 15, -15.0 * n / 15), order=1)
            a = img[None] * gain * t + 300.0 + 40.0 * t
            cube = a + rng.normal(0, np.sqrt(np.maximum(a, 1)) * 2.0)
            if prereg:
                bot, left = int(rng.integers(3, 20)), int(rng.integers(2, 15))
                cube[:, -bot:, :] = 0; cube[:, :, :left] = 0
            write_cube(os.path.join(folder, "f%02d.fits" % n), cube, t, "2026-08-12T18:29:%02d.0" % n)
            n += 1

ok = True
tmp = tempfile.mkdtemp(prefix="efhdr_am_")
try:
    for prereg in (True, False):
        d = os.path.join(tmp, "pre" if prereg else "drift"); os.makedirs(d); make(d, prereg)
        for mode in ("auto", "align", "keep"):
            shutil.rmtree(os.path.join(d, ".eclipseforgehdr"), ignore_errors=True)
            P = Prog(False)
            S.run(d, P, tier_mode="exposure", partialconv=False, align_mode=mode)
            sh = [float(l.split("shift rms")[1].split("px")[0]) for l in P.lines if "shift rms" in l]
            kept = any("kept as written" in l or "alignment off" in l for l in P.lines)
            want_kept = (mode == "keep") or (mode == "auto" and prereg)
            good = (kept == want_kept) and ((max(sh) == 0.0) == want_kept or (prereg and not want_kept))
            ok &= good
            print("simple  %-5s %-5s kept=%-5s max tier rms %5.2f  %s" % ("pre" if prereg else "drift", mode, kept, max(sh), "ok" if good else "FAIL"))
    # full pipeline on the registered set, auto and align
    d = os.path.join(tmp, "pre")
    for mode in ("auto", "align"):
        shutil.rmtree(os.path.join(d, ".eclipseforgehdr"), ignore_errors=True)
        P = Prog(False)
        try:
            PL.run(d, P, align_mode=mode, partialconv=False)
            lines = P.lines
            kept = any("kept as written" in l for l in lines)
            cross = [l for l in lines if "cross-tier alignment" in l or "network residual" in l]
            good = kept == (mode == "auto")
            ok &= good
            print("full    pre   %-5s kept=%-5s %s | %s" % (mode, kept, "ok" if good else "FAIL", " / ".join(c.strip() for c in cross)))
        except Exception as e:
            ok = False; print("full    pre   %-5s CRASH %s" % (mode, e))
            import traceback; traceback.print_exc()
    print("OK" if ok else "FAIL")
finally:
    shutil.rmtree(tmp, ignore_errors=True)
