#!/usr/bin/env python3
"""Run detail.build_layers end to end on a synthetic work directory.

WHY THIS EXISTS. 0.22.67 shipped a call to build_hill placed AFTER
`del lum_dn`, so every layer built correctly and then the run died with an
UnboundLocalError on the last step -- at the end of an 18-minute stack, on a
tester's machine. Each piece had been tested on its own; the function that
strings them together had not been run at all.

This builds a small fake HDR, runs the whole layer chain on it, and fails
loudly if anything raises. Seconds, not minutes. Run it before shipping
anything that touches detail.py:

    python tools/smoke_layers.py
"""
import json
import os
import shutil
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))


class _P:
    def __init__(self, verbose):
        self.verbose = verbose

    def log(self, m, f=None):
        if self.verbose:
            print("  |", m[:160])


def fixture(wd, H=520, W=760, R=55.0):
    cy, cx = H / 2.0, W / 2.0
    yy = np.arange(H, dtype=np.float32)[:, None] - cy
    xx = np.arange(W, dtype=np.float32)[None, :] - cx
    r = np.hypot(yy, xx)
    th = np.arctan2(yy, xx)
    rng = np.random.default_rng(9)
    thr = 1 + 0.25 * np.sin(37 * th) + 0.12 * np.sin(97 * th + 1)
    lum = (3e6 * (R / np.maximum(r, R)) ** 6 * thr).astype(np.float32)
    lum += (rng.normal(0, 1, (H, W)) * np.sqrt(np.maximum(lum, 1))).astype(np.float32)
    lum[r < R] = 0.0
    rgb = np.stack([lum * 1.15, lum, lum * 0.85], -1).astype(np.float32)
    for n, a in (("hdr_lum", lum), ("short_lum", lum)):
        np.save(os.path.join(wd, n + ".npy"), a)
    for n in ("hdr_rgb", "short_rgb"):
        np.save(os.path.join(wd, n + ".npy"), rgb)
    json.dump({"cy": cy, "cx": cx, "R": R, "Rmask": R + 6, "limb_margin": 6.0,
               "limb_prof": [float(R)] * 720, "rms": 1.0,
               "inner_cy": cy, "inner_cx": cx, "inner_R": R},
              open(os.path.join(wd, "geometry.json"), "w"))
    return lum


def main(verbose=True):
    from eclipseforgehdr import detail
    from eclipseforgehdr.render import Layers, render, DEFAULTS
    wd = tempfile.mkdtemp(prefix="efhdr_smoke_")
    try:
        fixture(wd)
        st = detail.build_layers(wd, _P(verbose), denoise="fine") or {}
        need = ["mgn", "fnrgf", "nafe", "rhef", "inner", "inner0", "prom",
                "pellett", "hill", "hill_log"]
        missing = [n for n in need
                   if not os.path.exists(os.path.join(wd, n + ".npy"))]
        if missing:
            raise SystemExit("layers missing after build_layers: %s" % missing)
        json.dump({k: v for k, v in st.items()},
                  open(os.path.join(wd, "report.json"), "w"), default=float)

        # and the render, both paths, with the Hill chain on and off
        ly = Layers(wd)
        for mix in (0.0, 1.0):
            P = dict(DEFAULTS)
            P["hillMix"] = mix
            for prev in (True, False):
                out = render(ly, P, preview=prev)
                if not np.isfinite(out).all():
                    raise SystemExit(
                        "render produced non-finite values (hillMix=%g, "
                        "preview=%s)" % (mix, prev))
        # the Hill VIEW: grayscale, in range, and actually moved by the
        # sliders. It is exported as a 16-bit TIFF, so a view that silently
        # fell back to something else would be shipped as the wrong picture.
        P = dict(DEFAULTS)
        hv = render(ly, P, preview=False, view="partialconv")
        if not np.isfinite(hv).all():
            raise SystemExit("partial-convolution view produced non-finite values")
        if not (np.allclose(hv[:, :, 0], hv[:, :, 1])
                and np.allclose(hv[:, :, 0], hv[:, :, 2])):
            raise SystemExit("partial-convolution view is not monochrome")
        P2 = dict(P); P2["hillGain"] = P["hillGain"] * 4
        if float(np.max(np.abs(render(ly, P2, preview=False, view="partialconv") - hv))) < 1e-4:
            raise SystemExit("partial-convolution view does not respond to hillGain")

        hb = st.get("hill") or {}
        if hb:
            worst = max(abs(x) for x in hb.get("limb_bias_after", [0.0]))
            print("hill limb bias after de-radialising, worst scale: %+.2f"
                  % worst)
            if worst > 1.0:
                raise SystemExit(
                    "the near-limb collar is back: a mask's own mean is as "
                    "large as its structure (%.2f)" % worst)
        print("OK — build_layers and render complete, %d layers written" % len(need))
    finally:
        shutil.rmtree(wd, ignore_errors=True)


if __name__ == "__main__":
    main(verbose="-q" not in sys.argv)
