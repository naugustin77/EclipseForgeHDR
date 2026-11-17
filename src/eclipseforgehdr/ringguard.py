"""THE RING GUARD (2026-11-16): concentric rings in the picture as it
will be exported, measured before it is saved.

Why. A picture with rings round the Moon, posted with the app's name on it, is
the failure users meet first. Until now rings were graded only on export,
after the file was written, and only beyond 1.2 R -- where they were never
found; the rings people see sit from the limb to about 1.6 R.

How. The picture is rendered at full resolution round the Moon, its display
brightness averaged in 0.5 px rings in eight 45-degree sectors, and each
sector's profile band-passed along the radius: "fine" 0.6-4 px (rings from
the detail layers and sharpening), "coarse" 4-30 px (hand-overs between
exposures, scaled down on small discs). A ring is what the sectors have in
common: the median of the eight band-passed profiles, its rms in % of full
scale ("amp"), and how well each sector follows it ("coh"). A prominence, a
streamer or a star sits in one or two sectors and drops out of the median.

Calibrated (2026-11-16) on full-resolution exports: one with fine rings at
the limb that were visible and later fixed (fine limb 0.78 %, coh 0.95; inner
0.14 %, coh 0.85), the same set fixed (0.21 / 0.03 %), a set with tier rings
and dust (fine inner 0.54 %, coh 0.72) and five good exports of other sets
(fine limb at most 0.23 %, fine inner at most 0.07 %, coh at most 0.30 there).
The limits sit about halfway in between.
"""
import numpy as np
from scipy import ndimage

ZONES = (("limb", 1.02, 1.15), ("inner", 1.15, 1.6), ("outer", 1.6, 2.2))
#: the coarse band reads the corona's own steep fall-off as a "ring" next to
#: the limb, so it starts further out
COARSE_FROM = 1.25
BANDS = (("fine", 0.6, 4.0), ("coarse", 4.0, 30.0))      # gaussian sigmas in px, along the radius


def measure(Y, cy, cx, R, rmask=None, nsec=8, step=0.5):
    """Y: display luminance (0..1), full grid. Returns {band: {zone: dict(amp, coh, at)}}:
    amp = rms of the ring profile common to the sectors (in % of full scale),
    coh = mean correlation of each sector's profile with that common profile,
    at = radius (R) of the largest common departure."""
    H, W = Y.shape
    r0 = max(float(rmask or 0.0) + 3.0, 1.01 * R)
    r1 = 2.4 * R          # beyond the last zone, so the band-pass has no edge there
    y0, y1 = int(max(cy - r1, 0)), int(min(cy + r1 + 1, H))
    x0, x1 = int(max(cx - r1, 0)), int(min(cx + r1 + 1, W))
    sub = np.asarray(Y[y0:y1, x0:x1], np.float64)
    yy, xx = np.mgrid[y0:y1, x0:x1]
    rr = np.hypot(yy - cy, xx - cx)
    th = (np.arctan2(yy - cy, xx - cx) + np.pi) / (2 * np.pi)
    sel = (rr >= r0) & (rr < r1)
    rb = ((rr[sel] - r0) / step).astype(np.int64)
    sb = np.minimum((th[sel] * nsec).astype(np.int64), nsec - 1)
    nb = int((r1 - r0) / step) + 1
    v = sub[sel]
    idx = sb * nb + rb
    s1 = np.bincount(idx, v, nsec * nb).reshape(nsec, nb)
    n1 = np.bincount(idx, None, nsec * nb).reshape(nsec, nb)
    prof = np.where(n1 > 3, s1 / np.maximum(n1, 1), np.nan)
    rad = (r0 + (np.arange(nb) + 0.5) * step) / R
    out = {}
    for bname, sa, sb_ in BANDS:
        if bname == "coarse":
            # in R: on a small disc 30 px is a third of the radius and reads
            # the limb's own fall-off
            sa, sb_ = max(1.5, min(sa, 0.01 * R)), max(6.0, min(sb_, 0.06 * R))
        bp = np.full_like(prof, np.nan)
        for k in range(nsec):
            p = prof[k]; ok = np.isfinite(p)
            if ok.sum() < 20:
                continue
            q = np.interp(np.arange(nb), np.flatnonzero(ok), p[ok])
            b = ndimage.gaussian_filter1d(q, sa / step) - ndimage.gaussian_filter1d(q, sb_ / step)
            b[~ok] = np.nan
            bp[k] = b
        out[bname] = {}
        for zname, za, zb in ZONES:
            if bname == "coarse":
                if zname == "limb":
                    continue
                za = max(za, COARSE_FROM)
            m = (rad >= za) & (rad < zb)
            if m.sum() < 8:
                continue
            Z = bp[:, m]
            good = np.isfinite(Z).mean(axis=1) > 0.8
            if good.sum() < 4:
                continue
            Z = Z[good]
            com = np.nanmedian(Z, axis=0)
            comf = np.nan_to_num(com)
            cs = []
            for z in Z:
                f = np.isfinite(z)
                if f.sum() > 8 and np.std(z[f]) > 0 and np.std(comf[f]) > 0:
                    cs.append(np.corrcoef(z[f], comf[f])[0, 1])
            amp = float(np.sqrt(np.nanmean(com ** 2)))
            out[bname][zname] = {"amp": round(100 * amp, 3), "coh": round(float(np.mean(cs)) if cs else 0.0, 2),
                                 "at": round(float(rad[m][np.nanargmax(np.abs(com))]), 3), "n": int(good.sum())}
    return out


#: (band, zone) -> (amp %, coh) above which the rings are taken as visible
LIMITS = {("fine", "limb"): (0.45, 0.80), ("fine", "inner"): (0.10, 0.60),
          ("fine", "outer"): (0.10, 0.60), ("coarse", "inner"): (0.50, 0.60),
          ("coarse", "outer"): (0.50, 0.60)}

#: the look's detail turned off, to tell rings from the look from rings in the stack
DETAIL_OFF = {"detailGain": 0.0, "hillMix": 0.0, "structure": 0.0, "bandDetail": 0.0,
              "innerDetail": 0.0, "clarity": 0.0}


def found(res):
    """[(band, zone, amp, coh, at)] over the limits, strongest first."""
    out = []
    for (b, z), (la, lc) in LIMITS.items():
        v = (res.get(b) or {}).get(z)
        if v and v["amp"] >= la and v["coh"] >= lc:
            out.append((b, z, v["amp"], v["coh"], v["at"], v["amp"] / la))
    return sorted(out, key=lambda t: -t[5])
