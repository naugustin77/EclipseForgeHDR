#!/usr/bin/env python3
"""One Moon for a single pass: the whole pipeline on a synthetic bracket shot
once through, where the Moon moves against a fixed corona from tier to tier.

What it checks (lab, 2026-10-15):
  * the pipeline picks "one Moon for all tiers" on such a set;
  * the merged lunar disc sits where the Moon was at that instant (the
    geometry the report and the layers use), within 1.5 px;
  * the merged limb stays round (fit rms under 1 px; the same frames with
    ECLIPSEFORGE_ONE_MOON=0 are run and printed beside it);
  * no hole opens in the corona just outside the disc.

    python tools/smoke_onemoon.py
"""
import json
import os
import shutil
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
sys.path.insert(0, os.path.dirname(__file__))
from smoke_pipeline import write_fits, _P      # noqa: E402

H, W, R = 640, 840, 130.0
SECS = (1 / 500, 1 / 125, 1 / 30, 1 / 8, 1 / 2)
DT = 8.0                       # s between tiers (one pass)
V = (-0.30, -0.18)             # Moon against the corona, px/s (11 px over the bracket, 9 % of R: enough smear for the moon-mask trial to want a mask)


def moon_at(t):
    return H / 2.0 + 20.0 + V[0] * t, W / 2.0 + 15.0 + V[1] * t


def bracket(folder):
    cy, cx = H / 2.0 + 20.0, W / 2.0 + 15.0          # the Sun
    yy = np.arange(H, dtype=np.float64)[:, None]
    xx = np.arange(W, dtype=np.float64)[None, :]
    r = np.hypot(yy - cy, xx - cx)
    th = np.arctan2(yy - cy, xx - cx)
    thr = 1.0 + 0.30 * np.sin(3 * th + 0.7) + 0.10 * np.sin(11 * th)
    # bright enough that the long tiers clip at the limb, as on a real bracket
    corona = 6.0e5 * (R / np.maximum(r, R * 0.98)) ** 5 * thr
    rng = np.random.default_rng(5)
    # fine structure fixed to the Sun, for the moon-free corona links to lock on
    from scipy import ndimage
    tex = ndimage.gaussian_filter(rng.normal(size=(H, W)), 2.0)
    corona = corona * (1.0 + 0.25 * tex / tex.std())
    n = 0
    for i, s in enumerate(SECS):
        for k in range(2):
            t = i * DT + k * 0.5
            my, mx = moon_at(t)
            rm = np.hypot(yy - my, xx - mx)
            edge = np.clip((rm - (R - 1.0)) / 2.0, 0.0, 1.0)      # a 2 px limb
            a = corona * edge * s + 4.0 + 30.0 * (1 - edge) * s   # dim earthshine
            a = a + rng.normal(0, np.sqrt(np.maximum(a, 1)))
            write_fits(os.path.join(folder, "f%02d.fits" % n), a, s,
                       "2026-08-12T18:30:%02d.%d" % (int(t), int(10 * (t % 1))))
            n += 1
    return n


def one(env, verbose):
    from eclipseforgehdr.pipeline import run, workdir
    folder = tempfile.mkdtemp(prefix="efhdr_onemoon_")
    old = {k: os.environ.get(k) for k in env}
    try:
        os.environ.update(env)
        bracket(folder)
        prog = _P(verbose)
        run(folder, prog, crop_pc=320, denoise="fine")
        wd = workdir(folder)
        stats = json.load(open(os.path.join(wd, "report.json")))
        geo = json.load(open(os.path.join(wd, "geometry.json")))
        hdr = np.load(os.path.join(wd, "hdr_lum.npy"))
        return stats, geo, hdr, prog.lines
    finally:
        for k, v in old.items():
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v
        shutil.rmtree(folder, ignore_errors=True)


def _limb(hdr, geo):
    from eclipseforgehdr.pipeline import fit_limb_rays
    f = fit_limb_rays(np.clip(hdr, 0, None), geo["cy"], geo["cx"], geo["R"], decim=1)
    yy, xx = np.mgrid[:hdr.shape[0], :hdr.shape[1]]
    rr = np.hypot(yy - f[0], xx - f[1])
    ring = (rr > f[2] * 1.02) & (rr < f[2] * 1.15)
    hole = float((hdr[ring] <= 0).mean())
    return f, hole


def main(verbose=False):
    on, geo1, h1, lines1 = one({}, verbose)
    off, geo0, h0, _ = one({"ECLIPSEFORGE_ONE_MOON": "0"}, verbose)
    om = on.get("one_moon") or {}
    for ln in lines1:
        if "one Moon" in ln or "moon-mask trial" in ln or "lunar track" in ln:
            print("  |", ln[:160])
    assert om.get("used"), ("one Moon not used", om, on.get("moon_mask"))
    # tier_time counts from the first tier's mean time (0.25 s in); FITS rows
    # are stored bottom-up and the reader flips them
    ty, tx = moon_at(float(om["t"]) + 0.25)
    ty = H - 1 - ty
    oy, ox = geo1.get("crop_origin", (0, 0)) if isinstance(geo1.get("crop_origin"), (list, tuple)) else (0, 0)
    dy, dx = geo1["cy"] + oy - ty, geo1["cx"] + ox - tx
    f1, hole1 = _limb(h1, geo1)
    f0, hole0 = _limb(h0, geo0)
    print(f"one Moon at t={om['t']} s: merged disc ({geo1['cy'] + oy:.1f}, {geo1['cx'] + ox:.1f}) "
          f"vs the Moon then ({ty:.1f}, {tx:.1f}): off by ({dy:+.1f}, {dx:+.1f}) px")
    print(f"merged limb fit rms: one Moon {f1[3]:.2f} px, each tier its own {f0[3]:.2f} px; "
          f"holes 1.02-1.15 R: {100 * hole1:.2f} % vs {100 * hole0:.2f} %")
    assert np.hypot(dy, dx) < 1.5, (dy, dx)
    # Not "better than each tier its own": on this synthetic set the own-disc
    # mask is rejected and the unmasked merge's limb is essentially one tier's,
    # round to 0.5 px. One Moon leaves a 0.6 px oval across the track (the
    # limb band built from fewer tiers along it). On the reference test set, with
    # hard clipping, one Moon gave 2.4 px against 3.6 px for the own discs.
    assert f1[3] < 1.0, f1[3]
    assert hole1 < 0.001, hole1
    print("OK - one Moon for all tiers: the disc sits at one instant, the limb is round, no holes")


if __name__ == "__main__":
    main(verbose="-v" in sys.argv)
