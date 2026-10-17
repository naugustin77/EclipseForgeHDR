import os, sys, json, tempfile, shutil, numpy as np
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, "..", "src")); sys.path.insert(0, _HERE)
import smoke_tier_mode as S
from smoke_pipeline import _P as Prog
from eclipseforgehdr import simple
tmp = tempfile.mkdtemp(prefix="efhdr_simple_meas_")
out = {}
try:
    for name, honest, mode in (("control", True, "exposure"), ("liar_exposure", False, "exposure"), ("liar_measured", False, "measured")):
        d = os.path.join(tmp, name); S._write(d, honest=honest)
        simple.run(d, Prog(False), tier_mode=mode, partialconv=False)
        wd = os.path.join(d, ".eclipseforgehdr")
        out[name] = (np.load(os.path.join(wd, "hdr_rgb.npy")), json.load(open(os.path.join(wd, "geometry.json"))), json.load(open(os.path.join(wd, "report.json"))))
finally:
    shutil.rmtree(tmp, ignore_errors=True)
g = out["control"][1]; ref = S._profile(out["control"][0], g["cy"], g["cx"], g["R"])
for k in ("liar_exposure", "liar_measured"):
    p = S._profile(out[k][0], g["cy"], g["cx"], g["R"]); r = p / ref; r = r / np.nanmedian(r)
    print(k, "tiers", len(out[k][2].get("tiers") or []), "shape err %.2f%%" % (100 * S._shape_err(p, ref)), np.round(r, 4))
