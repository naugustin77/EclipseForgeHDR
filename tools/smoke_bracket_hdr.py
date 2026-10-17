#!/usr/bin/env python3
"""HDR per bracket against stacking per exposure (2026-10-16).

Five brackets of four exposures over 54 s, the Moon sliding across a fixed
corona by a fifth of its radius, the mount drifting 3 px, each bracket through
its own transparency (haze), and a bright streak (a plane) across the corona
in ONE frame. Both modes run on the same frames:

  * "tiers": every exposure stacked across the brackets, then one merge;
  * "hdr":   every bracket merged on its own into the main frame, then the
             bracket HDRs combined with rejection.

Checked for "hdr": all five brackets merged and placed; the transparency it
reports follows the injected one (relative to the mean); the streak does not
print; the corona against the truth, sector by sector around the limb and
further out, is as flat as with "tiers".
Run: python tools/smoke_bracket_hdr.py
"""
import os, sys, json, shutil, tempfile, time
import numpy as np
from scipy import ndimage
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
sys.path.insert(0, os.path.dirname(__file__))
from smoke_pipeline import write_fits, _P
from smoke_moonmotion import corona, sectors, against_truth, H, W, R, MOVE, DRIFT

TRANS = (1.00, 0.92, 1.05, 0.97, 1.03)
STREAK = (2, 2)            # bracket 3, its third exposure (1/4 s)


def frames(folder, secs=(1 / 60, 1 / 15, 1 / 4, 1.0), pedestal=4.0):
    nb = len(TRANS)
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
    rng = np.random.default_rng(11)
    i = 0
    streak = None
    for b in range(nb):
        for k, s in enumerate(secs):
            tt = times[i]
            f = (tt - t_mid) / (times[-1] - times[0])
            rm = np.hypot(yy - cy, xx - (cx + MOVE * f))
            a = np.where(rm < R, 1200.0, base) * s * TRANS[b]
            if (b, k) == STREAK:
                # a 3 px wide line across the corona, +50 % of the corona under
                # it (below saturation, so the merge cannot drop it as clipped)
                d = np.abs((yy - cy) * 0.6 - (xx - cx) * 0.8 - 1.6 * R)
                streak = (d < 1.5)
                a = a + streak * 0.5 * np.where(rm < R, 0.0, base) * s * TRANS[b]
            a = ndimage.shift(a, (0.0, DRIFT * f), order=1, mode="nearest") + pedestal
            a = a + rng.normal(0, np.sqrt(np.maximum(a, 1)))
            write_fits(os.path.join(folder, "f%02d.fits" % i), a, s,
                       "2026-08-12T18:%02d:%02d.0" % (tt // 60, tt % 60))
            i += 1
    truth = base.copy()
    truth[np.hypot(yy - cy, xx - cx) < R] = 1200.0
    return truth[::-1].copy(), (H - 1 - cy, cx), streak[::-1].copy()


def one(mode, folder, truth, c, streak):
    from eclipseforgehdr import pipeline as PP
    wd = PP.workdir(folder)
    shutil.rmtree(wd, ignore_errors=True)
    prog = _P(False)
    t0 = time.time()
    PP.run(folder, prog, crop_pc=320, denoise="off", bracket_mode=mode)
    dt = time.time() - t0
    L = np.load(os.path.join(wd, "hdr_lum.npy"))
    st = json.load(open(os.path.join(wd, "report.json")))
    # the streak: merged / truth on its pixels against its neighbours, with
    # the merge registered onto the truth (as against_truth does), where the
    # 1/4 s frame carries the merge (1.4-2 R; further out the 1 s tier does)
    from skimage.registration import phase_cross_correlation
    T = truth[:L.shape[0], :L.shape[1]]
    yy, xx = np.mgrid[0:T.shape[0], 0:T.shape[1]]
    rr = np.hypot(yy - c[0], xx - c[1]) / R
    w = ((rr > 1.3) & (rr < 3.0)).astype(float)
    hpf = lambda x: (x - ndimage.gaussian_filter(x, 4)) * w
    sh = phase_cross_correlation(hpf(np.log(np.maximum(T, 1.0))), hpf(np.log(np.maximum(L, 1e-9))),
                                 upsample_factor=20, normalization=None)[0]
    q = ndimage.shift(L, sh, order=1) / np.maximum(T, 1e-9)
    s_ = streak[:L.shape[0], :L.shape[1]]
    band = (rr > 1.4) & (rr < 2.0)
    near = ndimage.binary_dilation(s_, iterations=6) & ~ndimage.binary_dilation(s_, iterations=3) & band
    ratio_streak = float(np.median(q[s_ & band]) / np.median(q[near]))
    return {"spread": against_truth(L, truth, c), "streak": ratio_streak, "seconds": dt,
            "bk": st.get("bracket_hdr"), "lines": prog.lines}


def main():
    folder = tempfile.mkdtemp(prefix="efhdr_bh_")
    try:
        truth, c, streak = frames(folder)
        res = {m: one(m, folder, truth, c, streak) for m in ("tiers", "hdr")}
    finally:
        shutil.rmtree(folder, ignore_errors=True)
    for m, r in res.items():
        print(f"{m:>5}: {r['seconds']:.0f} s; streak/neighbours {r['streak']:.3f}; merged/truth sector spread "
              + ", ".join(f"{k} R {v:.3f}" for k, v in r["spread"].items()))
    for l in res["hdr"]["lines"]:
        if "bracket" in l.lower():
            print("  |", l)
    bk = res["hdr"]["bk"] or {}
    assert bk.get("used") and bk.get("n") == len(TRANS), bk
    tr = np.array([b["transparency"] for b in bk["brackets"]])
    inj = np.array(TRANS) / np.mean(TRANS)
    print("transparency measured", np.round(tr / tr.mean(), 3), "injected", np.round(inj, 3))
    assert np.max(np.abs(tr / tr.mean() - inj)) < 0.02, (tr, inj)
    offs = [b["offset_px"] for b in bk["brackets"]]
    print("bracket offsets (full-res px)", offs)
    # the streak is +50 % in one frame; stacking per exposure clips it within
    # its tier, HDR per bracket rejects it across the brackets
    assert res["hdr"]["streak"] < 1.02, res["hdr"]["streak"]
    for key in res["tiers"]["spread"]:
        assert res["hdr"]["spread"][key] < res["tiers"]["spread"][key] * 1.05 + 0.02, (key, res)
    print("OK - HDR per bracket: all brackets placed, transparency followed, streak rejected, "
          "corona as flat as stacking per exposure")


if __name__ == "__main__":
    main()
