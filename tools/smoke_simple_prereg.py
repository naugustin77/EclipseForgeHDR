"""Simple stack on a set a writer already registered (Siril-style black
borders of different widths): the stack must keep the frames as written."""
import os, sys, tempfile, shutil, json
import numpy as np
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, "..", "src")); sys.path.insert(0, _HERE)
from smoke_pipeline import _P as Prog
from smoke_channel_floors import write_cube
import smoke_simple_align as T
from eclipseforgehdr import simple as S

base = T.scene()
tmp = tempfile.mkdtemp(prefix="efhdr_prereg_")
try:
    gain = np.array([0.7, 0.88, 1.0])[:, None, None]
    rng = np.random.default_rng(3); n = 0
    for i in range(T.NT):
        t = T.T0 * T.STEP ** i
        for k in range(2):
            a = base[None] * gain * t + 300.0 + 40.0 * t
            cube = a + rng.normal(0, np.sqrt(np.maximum(a, 1)) * 2.0)
            bot, left = int(rng.integers(3, 20)), int(rng.integers(2, 15))
            cube[:, -bot:, :] = 0; cube[:, :, :left] = 0
            write_cube(os.path.join(tmp, "f%02d.fits" % n), cube, t, "2026-08-12T18:29:%02d.0" % n)
            n += 1
    P = Prog(False)
    S.run(tmp, P, tier_mode="exposure", partialconv=False)
    reg = [l for l in P.lines if "already registered" in l]
    sh = [l for l in P.lines if "shift rms" in l]
    rep = json.load(open(os.path.join(tmp, ".eclipseforgehdr", "report.json")))
    txt = open(os.path.join(tmp, ".eclipseforgehdr", "report.txt")).read()
    ok = (len(reg) >= 1 and all("shift rms 0.00" in l for l in sh)
          and rep.get("simple_prealigned") is True and "Alignment: none." in txt)
    print(reg[0] if reg else "NO registered line")
    for l in sh: print("  ", l.strip())
    print("report simple_prealigned:", rep.get("simple_prealigned"),
          "| methods:", [l for l in txt.splitlines() if l.startswith("* Alignment")][0][:90])
    print("OK" if ok else "FAIL")
finally:
    shutil.rmtree(tmp, ignore_errors=True)
