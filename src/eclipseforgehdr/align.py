"""Cross-tier alignment and an objective measure of how well it worked.

Two published findings drive this module:

* Enhance, then correlate. The 2021 MNRAS eclipse-polarimetry pipeline
  (blind deconvolution -> NAFE -> cross-correlation) correlates ENHANCED
  images, not raw ones. On a short tier the corona is a few ADU above the
  read noise, so a raw correlation locks onto noise.
* The lunar edge is fine for coarse registration but "cannot be used as a
  reference feature for precise registration because it is dynamic during
  the TSE" (same paper). The Moon and the Sun are different targets.

Which means the fastest tiers, where the corona is essentially absent, have
to be tied to the sequence by the only high-SNR SOLAR features they contain:
the prominences.
"""
from __future__ import annotations
import numpy as np
from scipy import ndimage


# ---------- enhancement ----------

def enhance(img, sigma_lo=25.0, sigma_hi=2.0, eps=1e-3):
    """Local-contrast-normalized log image: structure at every brightness is
    brought to comparable amplitude, which is what phase correlation needs."""
    L = np.log1p(np.clip(img, 0, None) / max(float(np.median(np.abs(img))), eps))
    hp = L - ndimage.gaussian_filter(L, sigma_lo)
    sc = ndimage.gaussian_filter(np.abs(hp), sigma_lo)
    out = hp / np.maximum(sc, np.percentile(sc, 20) + eps)
    return ndimage.gaussian_filter(out, sigma_hi).astype(np.float32)


def signal_coverage(img, sat=None, region=None):
    """Fraction of the correlation region carrying real structure above the
    noise. Measured inside the crop the correlation actually uses, not over
    the whole frame (which is mostly empty sky). Below about 0.35 a tier has
    too little corona to phase-correlate on."""
    sub = img if region is None else img[region]
    med = float(np.median(sub))
    sig = 1.4826 * float(np.median(np.abs(sub - med)))
    good = ndimage.gaussian_filter(sub, 3) > med + 5 * sig
    if sat is not None:
        good &= ~(sat if region is None else sat[region])
    return float(good.mean())


# ---------- prominence anchors ----------

def find_prominences(img, cy, cx, R, n_max=4, na=1440, k=4.0, debug=False):
    """Compact bright features sitting on the limb. Returns patch centres
    (y, x) in image coordinates, brightest first.

    Prominences are SOLAR features, so they are valid alignment anchors; the
    chromosphere ring is not, because its visible extent is cut by the Moon.

    A prominence is a NARROW peak in azimuth confined to a thin band just
    outside the limb. Coronal streamers are the opposite -- broad in azimuth
    and radially extended to the frame edge -- so a detector based on "how far
    out does this azimuth stay bright" finds streamers, not prominences. What
    separates them is the azimuthal scale: subtract a broad azimuthal median
    from the near-limb brightness and only the compact features survive.
    """
    H, W = img.shape
    ang = np.linspace(0, 2 * np.pi, na, endpoint=False)
    rr = np.arange(R * 1.005, R * 1.13, 0.5, dtype=np.float32)
    ys = cy + rr[None, :] * np.sin(ang)[:, None]
    xs = cx + rr[None, :] * np.cos(ang)[:, None]
    ok = ((ys >= 0) & (ys <= H - 1) & (xs >= 0) & (xs <= W - 1)).all(axis=1)
    if ok.sum() < na // 4:
        return []
    P = ndimage.map_coordinates(np.log1p(np.clip(img, 0, None)),
                                [np.clip(ys, 0, H - 1).ravel(),
                                 np.clip(xs, 0, W - 1).ravel()],
                                order=1).reshape(na, -1)
    # brightness of the near-limb band per azimuth, robust to a single hot px
    b = np.percentile(P, 80, axis=1).astype(np.float32)
    b[~ok] = np.nan
    # broad azimuthal background: median over +/-15 deg, wrapped
    wide = max(9, int(round(na * 15.0 / 360.0))) | 1
    b3 = np.concatenate([b, b, b])
    fill = np.nanmedian(b) if np.isfinite(b).any() else 0.0
    b3f = np.where(np.isfinite(b3), b3, fill)
    bg = ndimage.median_filter(b3f, size=wide, mode="nearest")[na:2 * na]
    ex = np.where(np.isfinite(b), b - bg, 0.0).astype(np.float32)
    ex = ndimage.gaussian_filter1d(np.concatenate([ex] * 3), na / 720.0)[na:2 * na]
    # sigma over the IN-FRAME azimuths only. Out-of-frame ones were filled with
    # exact zeros above; with the disc near a frame edge (up to 75% of the
    # annulus may be outside and still pass the guard) those zeros drag the MAD
    # toward zero, the threshold follows, and every in-frame azimuth is flagged
    # as a prominence -- handing the NCC matcher arbitrary points on the limb.
    _e = ex[ok] if ok.any() else ex
    sig = 1.4826 * float(np.median(np.abs(_e - np.median(_e))))
    if not np.isfinite(sig) or sig <= 0:
        return []
    thr = k * sig
    prom = ex > thr
    if debug:
        print(f"    [prom] sigma {sig:.4f} thr {thr:.4f} "
              f"max excess {ex.max():.4f} ({ex.max() / sig:.1f} sigma), "
              f"{prom.sum()} of {na} azimuths flagged")
    if not prom.any():
        return []
    # contiguous sectors on the circle
    idx = np.flatnonzero(prom)
    groups, cur = [], [idx[0]]
    for a, bb in zip(idx[:-1], idx[1:]):
        if bb - a <= 3:
            cur.append(bb)
        else:
            groups.append(cur); cur = [bb]
    groups.append(cur)
    if len(groups) > 1 and prom[0] and prom[-1]:
        groups[0] = groups[-1] + groups[0]; groups.pop()
    min_w = max(3, int(round(na * 1.5 / 360.0)))       # at least ~1.5 deg wide
    scored = []
    for g in groups:
        if len(g) < min_w:
            continue
        gi = np.array(g) % na
        w = ex[gi]
        # brightness-weighted azimuth centroid, unwrapped so a group that
        # straddles 0 deg does not average to the opposite side of the disc
        a_un = np.unwrap(ang[gi])
        a_mid = float(np.sum(a_un * w) / np.sum(w))
        # radial centroid of the peak azimuth, so the patch sits ON the feature
        p = P[gi[int(np.argmax(w))]]
        p = np.clip(p - np.percentile(p, 10), 0, None)
        rad = float(np.sum(rr * p) / max(np.sum(p), 1e-9)) if p.sum() > 0 else R * 1.05
        scored.append((float(w.max()), a_mid, rad))
    scored.sort(reverse=True)
    out = []
    for _, a_mid, rad in scored[:n_max]:
        out.append((float(cy + rad * np.sin(a_mid)), float(cx + rad * np.cos(a_mid))))
    return out


def _subpixel_peak(c):
    j, i = np.unravel_index(int(np.argmax(c)), c.shape)
    dy = dx = 0.0
    if 0 < j < c.shape[0] - 1:
        a, b, d = c[j - 1, i], c[j, i], c[j + 1, i]
        den = a - 2 * b + d
        if abs(den) > 1e-12:
            dy = 0.5 * (a - d) / den
    if 0 < i < c.shape[1] - 1:
        a, b, d = c[j, i - 1], c[j, i], c[j, i + 1]
        den = a - 2 * b + d
        if abs(den) > 1e-12:
            dx = 0.5 * (a - d) / den
    return j + dy, i + dx, float(c[j, i])


def align_on_prominences(ref, tgt, anchors, half=60, search=80, min_peak=0.3,
                         debug=False, R=None):
    """Normalized cross-correlation of prominence patches. NCC is invariant to
    the brightness scaling between tiers, so raw (log) data is fine.

    Returns (dy, dx, spread_px, n_used) -- the shift that maps tgt onto ref,
    in the same sense as phase_cross_correlation(ref, tgt).
    """
    from skimage.feature import match_template
    # 60/80 px are 0.19 R / 0.26 R at the reference lunar radius (310 half-res
    # px) -- a patch that holds one prominence. Left absolute, at R=150 the same
    # patch spans 0.65-1.45 R and is dominated by the lunar edge, which is the
    # one feature that MOVES between tiers. Worse, every anchor then contains
    # the same limb, so they all agree and the spread gate cannot see it: the
    # tiers get registered to the Moon instead of the corona.
    if R:
        half = int(np.clip(round(0.19 * R), 20, 120))
        search = int(np.clip(round(0.26 * R), 28, 160))
    H, W = ref.shape
    A = np.log1p(np.clip(ref, 0, None)).astype(np.float32)
    B = np.log1p(np.clip(tgt, 0, None)).astype(np.float32)
    offs = []
    for (py, px) in anchors:
        y0, x0 = int(round(py)) - half, int(round(px)) - half
        y1, x1 = y0 + 2 * half, x0 + 2 * half
        wy0, wx0 = y0 - search, x0 - search
        wy1, wx1 = y1 + search, x1 + search
        if y0 < 0 or x0 < 0 or y1 > H or x1 > W:
            continue
        if wy0 < 0 or wx0 < 0 or wy1 > H or wx1 > W:
            continue
        tmpl = A[y0:y1, x0:x1]
        if float(tmpl.std()) < 1e-6:
            continue
        win = B[wy0:wy1, wx0:wx1]
        c = match_template(win, tmpl)
        py_, px_, peak = _subpixel_peak(c)
        if peak < min_peak:                # no believable match
            continue
        # template sits at (search, search) in the window when aligned
        offs.append((search - py_, search - px_, peak))
    if len(offs) < 1:
        return None
    o = np.array(offs)
    if len(o) >= 3:
        # drop anchors that disagree with the consensus -- a prominence can be
        # cut off by the frame edge or swamped by glare in the slow tiers
        med = np.median(o[:, :2], axis=0)
        d = np.hypot(o[:, 0] - med[0], o[:, 1] - med[1])
        mad = 1.4826 * float(np.median(d)) + 0.5
        keep = d < max(3.0 * mad, 2.0)
        if keep.sum() >= 2:
            o = o[keep]
    dy, dx = float(np.median(o[:, 0])), float(np.median(o[:, 1]))
    spread = float(np.hypot(o[:, 0].std(), o[:, 1].std())) if len(o) > 1 else 0.0
    if debug:
        for r in o:
            print(f"      anchor dy {r[0]:+6.2f} dx {r[1]:+6.2f} ncc {r[2]:.3f}")
    return dy, dx, spread, len(o)


def choose_shift(prom_result, corona_shift, max_spread=2.0, min_anchors=2):
    """Pick between the prominence-anchored and corona-correlated shift for one
    tier, and say why.

    Neither method is right everywhere, and which one is right is measurable
    rather than assumed. The prominences are the only solar features a very
    short tier contains, so they win there; in a long tier they are buried in
    inner-corona glare and their anchors start disagreeing with each other.
    That disagreement -- the spread across anchors -- is the gate. On the test
    set it stays under 1.9 px out to 1/13 s and then climbs to 4-9 px, which is
    exactly where the prominence answer stops being trustworthy.

    Returns (dy, dx, source).
    """
    if prom_result is None:
        return corona_shift[0], corona_shift[1], "corona"
    dy, dx, spread, n = prom_result
    if n < min_anchors or spread > max_spread:
        return corona_shift[0], corona_shift[1], "corona"
    return dy, dx, "prominences"


# ---------- alignment quality ----------

def stack_variance(tiers_aligned, cy, cx, R):
    """Coefficient of variation across photometrically-normalized aligned
    tiers -- the same thing Photoshop's 'Variance' stack mode shows. A tight,
    narrow rim at the limb means the tiers agree; a wide one means they do not.

    Returns (cov_map, stats dict).
    """
    if len(tiers_aligned) < 3:
        return None, {}
    S = np.stack(tiers_aligned, axis=0)
    good = np.isfinite(S) & (S > 0)
    n = good.sum(axis=0)
    Sm = np.where(good, S, np.nan)
    with np.errstate(invalid="ignore"):
        mu = np.nanmean(Sm, axis=0)
        sd = np.nanstd(Sm, axis=0)
        cov = np.where((n >= 3) & (mu > 0), sd / np.maximum(mu, 1e-9), np.nan)
    H, W = cov.shape
    yy = np.arange(H, dtype=np.float32)[:, None] - cy
    xx = np.arange(W, dtype=np.float32)[None, :] - cx
    r = np.hypot(yy, xx)
    out = {}
    # NOTE: every measurement here is taken OUTSIDE the limb. Inside the disc
    # there is no signal, so the coefficient of variation is large by
    # construction and says nothing about alignment; a rim measured across a
    # band that includes the disc interior mostly reports that plateau.
    # HOW MANY TIERS ACTUALLY CARRY THE NUMBER. This is not bookkeeping; it is
    # what makes the number mean anything, and leaving it out produced a false
    # alarm that stood for three releases.
    #
    # Clifton's 250 mm set reported cov_limb 0.793 against 0.037 on his 360 mm
    # set of the same eclipse, and that 20x gap was read as veiling glare --
    # written into the report, the changelog and the backlog. Rebuilding the
    # same statistic from that run's OWN exported tiers gives **0.021**. The
    # 360 mm set rebuilds at 0.039 against the pipeline's 0.037, so the method
    # is sound; the 250 mm number is not.
    #
    # The cause is the `n >= 3` line above meeting a bracket whose shortest tier
    # is 1/125 s. Between 1.00 and 1.10 R that set has only TWO tiers with
    # unclipped signal. Every pixel where a third one survived the clipping cut
    # is a pixel where a partly-clipped tier happened to be dim -- the exclusion
    # keeps precisely the dim tail of a blown tier -- so the median is taken
    # over a biased remnant instead of over agreeing tiers. The 360 mm set has
    # five or six clean tiers there and never enters that regime.
    #
    # So: count the clean tiers, and refuse to report rather than report a
    # number built on two of them plus contamination.
    for lo, hi, tag in ((1.00, 1.10, "limb"), (1.3, 2.2, "corona")):
        band = (r > lo * R) & (r < hi * R)
        m = band & np.isfinite(cov)
        if band.sum() > 500:
            out[f"tiers_{tag}"] = float(np.median(n[band]))
        if m.sum() > 500:
            if out.get(f"tiers_{tag}", 0.0) >= 3.0:
                out[f"cov_{tag}"] = float(np.median(cov[m]))
            else:
                # the statistic is not measurable here, and saying so is the
                # honest output -- a number from the contaminated remnant is
                # worse than no number
                out[f"cov_{tag}_unmeasurable"] = float(
                    out.get(f"tiers_{tag}", 0.0))
    # Width of the high-variance rim just outside the limb: how far the CoV
    # excess over the coronal baseline persists. This is the number the user
    # reads off a Photoshop 'Variance' layer -- a misaligned stack smears the
    # limb disagreement outward, an aligned one keeps it to a few px.
    m = np.isfinite(cov) & (r >= R) & (r < 1.60 * R)
    if m.sum() > 500:
        ri = np.round(r[m]).astype(np.int32)
        v = cov[m].astype(np.float64)
        n0 = ri.max() + 1
        cnt = np.bincount(ri, minlength=n0).astype(np.float64)
        s = np.bincount(ri, weights=v, minlength=n0)
        prof = np.where(cnt > 20, s / np.maximum(cnt, 1), np.nan)
        rad = np.arange(n0, dtype=np.float32)
        band = np.isfinite(prof) & (rad >= R) & (rad < 1.60 * R)
        if band.sum() > 20:
            pr, ra = prof[band], rad[band]
            base = float(np.median(pr[ra > 1.35 * R])) if (ra > 1.35 * R).sum() > 5 \
                else float(np.min(pr))
            pk = float(pr[0])          # the rim starts at the limb by definition
            if pk > base * 1.15:
                half = 0.5 * (pk + base)
                below = np.flatnonzero(pr <= half)
                out["rim_width_px"] = float(ra[below[0]] - ra[0]) if len(below) \
                    else float(ra[-1] - ra[0])
                out["rim_peak_cov"] = pk
                out["rim_base_cov"] = base
    return cov, out


def limb_transition_width(lum, cy, cx, R, na=72, ref_px=5.0):
    """20-80% width of the merged limb, per azimuth. A sharp lunar edge should
    be a few px; tens of px means the tiers were blended out of register.

    THE REFERENCE LEVEL HAS TO SIT JUST OUTSIDE THE EDGE (fixed 0.16.3).

    This used `hi = percentile(v, 95)` over 0.75-1.35 R. The inner corona is
    still climbing steeply tens of px beyond the limb -- on the reference
    bracket it peaks at R+27 -- so that percentile is the near-limb CORONA
    PEAK, and "80% of it" lands far outside the lunar edge. The number then
    measures the corona's own radial gradient rather than the edge.

    The proof is that it scaled with wherever the reference was taken, on the
    reference bracket's own merged luminance:

        reference at   R+5   R+10   R+15   R+20   R+30
        median width   6.5    9.8   12.5   16.0   20.0 px
        p90 width      9.0   11.5   16.0   21.0   28.0 px

    A real edge width converges as the reference moves out. This one grows
    without bound, which is what says it is not an edge measurement.

    It matters because `margin = 0.9 * ramp` sets the disc mask, so an
    inflated width hides real corona: 25.1 px of margin on that bracket where
    the edge itself is 6.5 px wide. Both testers reported losing prominences,
    and the largest one sat inside that annulus.
    """
    H, W = lum.shape
    ang = np.linspace(0, 2 * np.pi, na, endpoint=False)
    rr = np.arange(0.75 * R, 1.35 * R, 0.5, dtype=np.float32)
    ys = cy + rr[None, :] * np.sin(ang)[:, None]
    xs = cx + rr[None, :] * np.cos(ang)[:, None]
    ok = ((ys >= 0) & (ys <= H - 1) & (xs >= 0) & (xs <= W - 1)).all(axis=1)
    P = ndimage.map_coordinates(lum, [np.clip(ys, 0, H - 1).ravel(),
                                      np.clip(xs, 0, W - 1).ravel()],
                                order=1).reshape(na, -1)
    ws = []
    for i in range(na):
        if not ok[i]:
            continue
        v = P[i]
        lo = float(np.median(v[:30]))
        # the corona immediately outside the edge, not the near-limb peak
        _c = int(round((R - 0.75 * R) / 0.5))
        _k = int(round(ref_px / 0.5))
        hi = float(np.median(v[max(_c + _k - 3, 0):_c + _k + 4]))
        if hi <= lo * 1.3:
            continue
        r20 = rr[int(np.argmax(v > lo + 0.2 * (hi - lo)))]
        r80 = rr[int(np.argmax(v > lo + 0.8 * (hi - lo)))]
        if r80 > r20:
            ws.append(float(r80 - r20))
    if not ws:
        return None
    return {"limb_width_med": float(np.median(ws)),
            "limb_width_p90": float(np.percentile(ws, 90))}
