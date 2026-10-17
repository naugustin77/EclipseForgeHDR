#!/usr/bin/env python3
"""Repeated brackets with the Moon moving across the corona (2026-10-14).

Five brackets of four exposures, round-robin, over 54 s. The corona is fixed
to the Sun (a radial falloff with fine texture, what an alignment holds on
to); the Moon slides across it by a fifth of its radius and sits exactly on
the Sun at the middle of the shoot; the mount drifts 3 px, Moon and corona
together. All frames get the SAME light, so:

  * every per-frame brightness gain must come out 1.00 -- the moving limb
    zone must not be read as transparency;
  * the merged corona, registered onto the truth, must follow it all round
    the limb: the per-sector ratio merged / truth in 24 sectors, just outside
    the limb and further out, as flat as with the Moon held still.

Run with the per-frame Moon on (automatic: repeated brackets) and off
(ECLIPSEFORGE_PERFRAME_MOON=0), printed side by side; "on" is asserted, and
asserted to beat "off" at the limb.
Run: python tools/smoke_moonmotion.py
"""
import os, sys, json, shutil, tempfile
import numpy as np
from scipy import ndimage
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
sys.path.insert(0, os.path.dirname(__file__))
from smoke_pipeline import write_fits, _P

H, W, R = 420, 560, 46.0
MOVE = 0.2 * R                       # Moon against the corona, full-res px over the shoot
DRIFT = 3.0                          # mount drift, Moon and corona together


def corona():
    cy, cx = H / 2.0, W / 2.0
    yy = np.arange(H, dtype=np.float64)[:, None]
    xx = np.arange(W, dtype=np.float64)[None, :]
    r = np.hypot(yy - cy, xx - cx)
    th = np.arctan2(yy - cy, xx - cx)
    tex = ndimage.gaussian_filter(np.random.default_rng(3).normal(size=(H, W)), 2.0)
    tex /= tex.std()
    return 6.0e5 * (R / np.maximum(r, R)) ** 4 * (1.0 + 0.04 * np.sin(23 * th)) * (1.0 + 0.15 * tex)


def frames(folder, secs=(1 / 60, 1 / 15, 1 / 4, 1.0), nb=5, pedestal=4.0):
    cy, cx = H / 2.0, W / 2.0
    yy = np.arange(H, dtype=np.float64)[:, None]
    xx = np.arange(W, dtype=np.float64)[None, :]
    base = corona()
    times, t = [], 0
    for b in range(nb):
        for s in secs:
            times.append(t); t += 2
        t += 4
    t_mid = 0.5 * (times[0] + times[-1])
    rng = np.random.default_rng(7)
    for i, s in enumerate(secs * nb):
        tt = times[i]
        f = (tt - t_mid) / (times[-1] - times[0])
        rm = np.hypot(yy - cy, xx - (cx + MOVE * f))
        a = np.where(rm < R, 1200.0, base) * s
        a = ndimage.shift(a, (0.0, DRIFT * f), order=1, mode="nearest") + pedestal
        a = a + rng.normal(0, np.sqrt(np.maximum(a, 1)))
        write_fits(os.path.join(folder, "f%02d.fits" % i), a, s,
                   "2026-08-12T18:%02d:%02d.0" % (tt // 60, tt % 60))
    # the truth at mid-shoot, as the pipeline will hold it (FITS is bottom-up)
    truth = base.copy()
    truth[np.hypot(yy - cy, xx - cx) < R] = 1200.0
    return truth[::-1].copy(), (H - 1 - cy, cx)


def sectors(img, c, lo, hi, n=24):
    yy, xx = np.mgrid[0:img.shape[0], 0:img.shape[1]]
    rad = np.hypot(yy - c[0], xx - c[1]) / R
    s = ((np.arctan2(yy - c[0], xx - c[1]) + np.pi) / (2 * np.pi) * n).astype(int) % n
    ring = (rad > lo) & (rad < hi)
    return np.array([np.median(img[ring & (s == k)]) for k in range(n)])


def against_truth(L, truth, c):
    """Register the merge onto the truth on the outer corona, then the
    max/min of merged / truth over 24 sectors in three rings."""
    from skimage.registration import phase_cross_correlation
    T = truth[:L.shape[0], :L.shape[1]]
    yy, xx = np.mgrid[0:T.shape[0], 0:T.shape[1]]
    r = np.hypot(yy - c[0], xx - c[1]) / R
    out = (r > 1.5) & (r < 3.0)
    lt = np.log(np.maximum(T, 1.0))
    ll = np.log(np.maximum(L * np.median(T[out]) / max(np.median(L[out]), 1e-9), 1.0))
    w = ((r > 1.3) & (r < 3.0)).astype(float)
    hpf = lambda x: (x - ndimage.gaussian_filter(x, 4)) * w
    sh = phase_cross_correlation(hpf(lt), hpf(ll), upsample_factor=20, normalization=None)[0]
    Lr = ndimage.shift(L, sh, order=1)
    return {f"{lo}-{hi}": float(np.max(q) / np.min(q)) for lo, hi in ((1.0, 1.04), (1.1, 1.25), (1.25, 1.6))
            for q in [sectors(Lr, c, lo, hi) / sectors(T, c, lo, hi)]}


def one(env):
    from eclipseforgehdr import pipeline as PP
    folder = tempfile.mkdtemp(prefix="efhdr_mm_")
    old = os.environ.get("ECLIPSEFORGE_PERFRAME_MOON")
    try:
        truth, c = frames(folder)
        if env is None:
            os.environ.pop("ECLIPSEFORGE_PERFRAME_MOON", None)
        else:
            os.environ["ECLIPSEFORGE_PERFRAME_MOON"] = env
        prog = _P(False)
        PP.run(folder, prog, crop_pc=320, denoise="off")
        wd = PP.workdir(folder)
        q = json.load(open(os.path.join(wd, "quality.json")))
        gains = [g for v in q.values() for g in (v.get("frame_gain") or {}).values()]
        L = np.load(os.path.join(wd, "hdr_lum.npy"))
        if L.ndim == 3:
            L = L.mean(axis=2)
        st = json.load(open(os.path.join(wd, "report.json")))
        return {"gains": [g for g in gains if g is not None],
                "dropped": sum(1 for g in gains if g is None),
                "spread": against_truth(L, truth, c), "pfm": st.get("perframe_moon"),
                "lines": [l for l in prog.lines if "Moon" in l or "lunar track" in l]}
    finally:
        if old is None:
            os.environ.pop("ECLIPSEFORGE_PERFRAME_MOON", None)
        else:
            os.environ["ECLIPSEFORGE_PERFRAME_MOON"] = old
        shutil.rmtree(folder, ignore_errors=True)


def main():
    res = {"off": one("0"), "on": one(None)}
    for tag in ("off", "on"):
        r = res[tag]
        print(f"{tag:>3}: frame gains {min(r['gains']):.3f}..{max(r['gains']):.3f} "
              f"({r['dropped']} dropped); merged/truth sector spread "
              + ", ".join(f"{k} R {v:.3f}" for k, v in r["spread"].items()))
    for l in res["on"]["lines"]:
        print("  |", l)
    on, off = res["on"], res["off"]
    v = (on["pfm"] or {}).get("v_px_s")
    print(f"Moon velocity measured {np.hypot(*v):.3f} px/s, put in {MOVE / 54.0:.3f}" if v else "no velocity")
    assert on["pfm"] and on["pfm"].get("all"), on["pfm"]
    assert v and abs(np.hypot(*v) / (MOVE / 54.0) - 1.0) < 0.15, v
    assert on["dropped"] == 0, on
    assert min(on["gains"]) > 0.99 and max(on["gains"]) < 1.01, on["gains"]
    # was < 0.85 x off. With prominence links dropped for tiers that hold a
    # moon-free corona link (2026-10-15) "off" fell 1.98 -> 1.42 and "on"
    # 1.54 -> 1.36: much of the old gap was the prominence links pulling
    # "off" toward the Moon. "on" must still win, and stay at its new level.
    assert on["spread"]["1.0-1.04"] < 0.98 * off["spread"]["1.0-1.04"], (on["spread"], off["spread"])
    assert on["spread"]["1.0-1.04"] < 1.45, on["spread"]
    assert on["spread"]["1.1-1.25"] < 1.10 and on["spread"]["1.25-1.6"] < 1.10, on["spread"]
    print("OK - Moon placed per frame: gains flat, corona follows the truth at the limb")


if __name__ == "__main__":
    main()
