#!/usr/bin/env python3
"""Repeated brackets end to end (2026-10-14).

Three brackets of four exposures, shot round-robin like SetNC or Umbra, the
second and third bracket through 10 % less and 10 % more light (haze). The
pipeline must find 3 brackets from the time stamps and exposure order, and the
per-bracket brightness it reports from the frame gains must follow the injected
transparency. Run: python tools/smoke_brackets.py
"""
import os, sys, json, shutil, tempfile
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
sys.path.insert(0, os.path.dirname(__file__))
from smoke_pipeline import write_fits, _P

def brackets(folder, H=420, W=560, R=46.0, secs=(1/60, 1/15, 1/4, 1.0),
             gains=(1.0, 0.9, 1.1), pedestal=4.0, peak=6.0e5):
    cy, cx = H / 2.0, W / 2.0
    yy = np.arange(H, dtype=np.float64)[:, None] - cy
    xx = np.arange(W, dtype=np.float64)[None, :] - cx
    r = np.hypot(yy, xx); th = np.arctan2(yy, xx)
    base = peak * (R / np.maximum(r, R)) ** 4 * (1.0 + 0.30 * np.sin(3 * th + 0.7) + 0.10 * np.sin(11 * th))
    base[r < R] = 0.0
    rng = np.random.default_rng(5); n = 0; t = 0
    for b, g in enumerate(gains):
        for s in secs:
            a = base * s * g + pedestal
            a = a + rng.normal(0, np.sqrt(np.maximum(a, 1)))
            write_fits(os.path.join(folder, "f%02d.fits" % n), a, s,
                       "2026-08-12T18:%02d:%02d.0" % (t // 60, t % 60))
            n += 1; t += 2
        t += 4
    return n

def main():
    from eclipseforgehdr.pipeline import run, workdir
    folder = tempfile.mkdtemp(prefix="efhdr_brk_")
    try:
        nf = brackets(folder)
        prog = _P(False)
        run(folder, prog, crop_pc=320, denoise="off")
        st = json.load(open(os.path.join(workdir(folder), "report.json")))
        bk = st.get("brackets") or {}
        for l in prog.lines:
            if "bracket" in l: print("  |", l)
        assert bk.get("repeated") and bk.get("n") == 3, bk
        rows = bk.get("brightness") or []
        longs = [r["long"] for r in rows]
        print("brightness, long tiers:", longs, " injected (median-normalised): [1.0, 0.9, 1.1]")
        assert len(rows) == 3 and all(x is not None for x in longs)
        assert abs(longs[1] - 0.9) < 0.03 and abs(longs[2] - 1.1) < 0.03, longs
        txt = open(os.path.join(workdir(folder), "report.txt")).read()
        assert "repeated brackets" in txt
        print("OK - %d frames, 3 brackets found, transparency followed" % nf)
    finally:
        shutil.rmtree(folder, ignore_errors=True)

if __name__ == "__main__":
    main()
