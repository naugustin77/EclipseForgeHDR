"""HDR PER BRACKET (lab, 2026-10-16).

The default for repeated brackets stacks every exposure across all brackets
(frame gains and the Moon placed per frame), then merges the stacked tiers
into one HDR. This is the other order, the one some published workflows use: one HDR per
bracket, then the HDRs aligned and averaged.

  1. Each bracket is merged ON ITS OWN by the normal pipeline (pipeline.run
     with `_sub`): its own cross-tier alignment, photometric ladder, merge
     weights and Moon -- every frame of it shot within seconds, so the ladder
     carries no change of sky or haze, and the Moon barely moves.
  2. It is merged straight into the MAIN run's frame, so no second resample
     softens it: the bracket keeps its tiers' placement relative to each
     other, and the whole bracket moves by one offset, measured moon-free
     against the main run's tier stacks of the same exposures
     (`bracket_offset`).
  3. The bracket HDRs are scaled to the main merge (one factor each: the
     transparency of that bracket), and combined per pixel: weighted by the
     merge weight each HDR has there (its effective exposure, so a bracket
     missing its long tiers counts less in the faint corona), near the limb
     only the HDRs whose Moon matches the main merge's Moon (the per-frame
     Moon's weight), and with a per-pixel rejection across brackets against a
     noise model, so a bracket hit by haze, a plane or a cosmic ray does not
     print.

The main run's per-exposure merge is still made: it is the frame the brackets
are placed in, the reference they are scaled to, and the source of everything
built from the tier stacks after the merge (the inner short stack, the
prominence stack). Only the merged HDR itself is replaced.
"""
from __future__ import annotations
import os, time, shutil
import numpy as np

BRACKET_KAPPA = 3.0        # rejection threshold across brackets, in noise sigma
BRACKET_MIN_FRAMES = 3     # a bracket with fewer frames is not merged on its own
BRACKET_BLOCK = 192        # rows per block in the combine (memory)


def _subprogress(parent, tag, f0, f1):
    """A Progress for one bracket's run: its lines are kept, warnings and
    failures go to the main log with the bracket named, the main bar moves
    from f0 to f1, and a cancel in the main run stops it."""
    from .pipeline import Progress, split_level, Cancelled

    class _Sub(Progress):
        def log(self, msg, frac=None, level=None):
            m, lev = split_level(msg, level)
            self.lines.append(m)
            self.levels.append(lev)
            self.t_line = time.time()
            self.stamps.append(self.t_line - self.t0)
            try:
                parent.t_line = self.t_line          # the main run is alive
                if frac is not None and hasattr(parent, "bar"):
                    parent.frac = max(parent.frac, parent.bar(f0 + (f1 - f0) * float(np.clip(frac, 0, 1))))
            except Exception:
                pass
            if lev in ("warn", "fail"):
                parent.log(f"[{lev}] {tag}: {m}", None)
            if getattr(parent, "cancel_requested", False):
                parent.cancelled = True
                self.cancelled = True
                raise Cancelled("cancelled by the user")

    s = _Sub()
    return s


def bracket_offset(stacks_half, sat_half, abs_shift, secs, ref, progress):
    """One offset (half-res dy, dx) that puts this bracket's merge frame onto
    the main run's: for every tier with a main tier of the same exposure, the
    moon-free corona shift between this bracket's frame and the main stack,
    plus the main tier's placement, less this tier's own placement. The
    median over the tiers; returns (offset or None, info)."""
    from .pipeline import (corona_shift_moonfree, _frame_moon_fit, crop_around,
                           MOONFREE_MIN_Q)
    rsecs = list(ref["secs"])
    mcache = ref.setdefault("_moon_half", {})
    got = []
    for s in secs:
        sm = min(rsecs, key=lambda x: abs(np.log(float(x) / float(s))))
        if abs(np.log(float(sm) / float(s))) > 0.05:
            continue
        A = ref["stacks_half"].get(sm)
        B = stacks_half.get(s)
        if A is None or B is None or A.shape != B.shape:
            continue
        if sm not in mcache:
            mcache[sm] = _frame_moon_fit(A)
        ma, mb = mcache[sm], _frame_moon_fit(B)
        if ma is None or mb is None:
            continue
        R = 0.5 * (ma[2] + mb[2])
        S = int(min(1200, min(A.shape)))
        ca, y0, x0 = crop_around(A, 0.5 * (ma[0] + mb[0]), 0.5 * (ma[1] + mb[1]), S)
        cb = B[y0:y0 + ca.shape[0], x0:x0 + ca.shape[1]]
        h, w = min(ca.shape[0], cb.shape[0]), min(ca.shape[1], cb.shape[1])
        ca, cb = ca[:h, :w], cb[:h, :w]
        _sa = ref["sat_half"].get(sm)
        _sb = sat_half.get(s)
        sa = _sa[y0:y0 + h, x0:x0 + w] if _sa is not None else None
        sb = _sb[y0:y0 + h, x0:x0 + w] if _sb is not None else None
        try:
            r = corona_shift_moonfree(ca, cb, sa, sb, (ma[0] - y0, ma[1] - x0),
                                      (mb[0] - y0, mb[1] - x0), R, with_q=True)
        except Exception:
            r = None
        if r is None or r[2] < MOONFREE_MIN_Q:
            continue
        dy, dx, q = r
        gy = dy + ref["abs_shift"][sm][0] - abs_shift[s][0]
        gx = dx + ref["abs_shift"][sm][1] - abs_shift[s][1]
        got.append((float(s), gy, gx, q))
    if not got:
        return None, {"why": "no tier matched the main run's stacks on the corona"}
    gy = float(np.median([g[1] for g in got]))
    gx = float(np.median([g[2] for g in got]))
    sc = float(np.median([np.hypot(g[1] - gy, g[2] - gx) for g in got])) if len(got) > 1 else 0.0
    info = {"dy": 2.0 * gy, "dx": 2.0 * gx, "tiers": len(got),
            "scatter_px": round(2.0 * sc, 3),
            "q_median": round(float(np.median([g[3] for g in got])), 1)}
    return (gy, gx), info


def _half_lum(rgb):
    h2, w2 = rgb.shape[0] // 2, rgb.shape[1] // 2
    L = (0.2126 * rgb[:2 * h2, :2 * w2, 0] + 0.7152 * rgb[:2 * h2, :2 * w2, 1]
         + 0.0722 * rgb[:2 * h2, :2 * w2, 2])
    return np.asarray(L, np.float32).reshape(h2, 2, w2, 2).mean(axis=(1, 3))


def _half_rgb(rgb):
    h2, w2 = rgb.shape[0] // 2, rgb.shape[1] // 2
    return np.asarray(rgb[:2 * h2, :2 * w2], np.float32).reshape(h2, 2, w2, 2, 3).mean(axis=(1, 3))


def _moon(Lh, seed=None):
    """(cy, cx, R) of the Moon on a half-res luminance, full-res px."""
    from .pipeline import find_disc, fit_limb_rays
    d = seed if seed is not None else find_disc(Lh)
    if d is None:
        return None
    f = fit_limb_rays(np.clip(Lh, 0, None), d[0], d[1], d[2], decim=1)
    if f is None:
        return None
    return float(f[0]) * 2.0, float(f[1]) * 2.0, float(f[2]) * 2.0


def _moon_w_rows(r0, r1, W, cy, cx, ty, tx, R, margin, soft):
    """The per-frame Moon's weight (pipeline._pfm_weight) for rows r0..r1 of
    a full-res frame: 1 away from both limbs, and in the limb band a Gaussian
    in the difference of the distances to this HDR's limb and the target's."""
    from .pipeline import _PFM_EPS
    w = np.ones((r1 - r0, W), np.float32)
    pad = R + margin + 4
    y0 = int(max(r0, np.floor(min(cy, ty) - pad)))
    y1 = int(min(r1, np.ceil(max(cy, ty) + pad) + 1))
    x0 = int(max(0, np.floor(min(cx, tx) - pad)))
    x1 = int(min(W, np.ceil(max(cx, tx) + pad) + 1))
    if y1 <= y0 or x1 <= x0:
        return w
    yy, xx = np.mgrid[y0:y1, x0:x1].astype(np.float32)
    m = float(max(margin, 1.0))
    dp = np.clip(np.hypot(yy - cy, xx - cx) - R, -m, m)
    dt = np.clip(np.hypot(yy - ty, xx - tx) - R, -m, m)
    ww = np.exp(-0.5 * ((dp - dt) / float(soft)) ** 2)
    w[y0 - r0:y1 - r0, x0:x1] = np.maximum(ww, _PFM_EPS)
    return w


def hdr_per_bracket(folder, wd, hdr, bracket_of, sec_of, epoch, run_kw, ref,
                    progress, stats, f0=0.86, f1=0.93):
    """Merge every bracket on its own into the main frame, then combine.
    `hdr` is the main run's merge (cropped, scene-linear); returns the
    combined HDR of the same shape, or `hdr` itself when fewer than two
    brackets could be merged."""
    from .pipeline import run as _run, Cancelled, _exp_name
    groups = {}
    for p, b in bracket_of.items():
        groups.setdefault(int(b), []).append(p)
    order = sorted(groups)
    bdir = os.path.join(wd, "brackets")
    shutil.rmtree(bdir, ignore_errors=True)
    os.makedirs(bdir, exist_ok=True)
    (cy0, cx0), kh, kw = ref["crop"]
    ep = [e for e in (epoch.get(p) for p in bracket_of) if e is not None]
    tfirst = min(ep) if ep else None
    progress.log(f"HDR per bracket (experimental): {len(order)} brackets, each merged on its own "
                 f"into this frame, then aligned brackets combined", None)
    info = {"used": False, "brackets": []}
    stats["bracket_hdr"] = info
    items = []
    t_start = time.time()
    for k, b in enumerate(order):
        ps = sorted(groups[b], key=lambda p: (epoch.get(p) or 0.0, os.path.basename(p)))
        nexp = len({round(float(np.log(max(sec_of.get(p) or 1e-9, 1e-9))), 2) for p in ps})
        rec = {"bracket": b + 1, "frames": len(ps), "exposures": nexp}
        tt = [epoch.get(p) for p in ps if epoch.get(p) is not None]
        if tt and tfirst is not None:
            rec["t"] = round(float(np.mean(tt)) - tfirst, 1)
        if len(ps) < BRACKET_MIN_FRAMES or nexp < BRACKET_MIN_FRAMES:
            rec["skipped"] = f"{len(ps)} frames, {nexp} exposures"
            info["brackets"].append(rec)
            progress.log(f"bracket {b + 1}: {len(ps)} frames, {nexp} exposures: "
                         f"too few to merge on its own, left out", None)
            continue
        sub = _subprogress(progress, f"bracket {b + 1}",
                           f0 + (f1 - f0) * k / len(order), f0 + (f1 - f0) * (k + 1) / len(order))
        t0 = time.time()
        try:
            res = _run(folder, sub, **run_kw, bracket_mode="tiers",
                       _sub={"paths": ps, "wd": os.path.join(bdir, f"b{b + 1}"),
                             "masters_wd": wd, "ref": ref})
        except Cancelled:
            raise
        except Exception as e:
            rec["failed"] = str(e)[:200]
            info["brackets"].append(rec)
            progress.log(f"[warn] bracket {b + 1}: not merged ({e}); left out", None)
            continue
        h = res["hdr"][cy0:cy0 + kh, cx0:cx0 + kw]
        w = np.array(res["wsum"][cy0:cy0 + kh, cx0:cx0 + kw], np.float32)
        # the bracket's own trimmed rectangle (its tiers' shifts), in this frame
        (by0, bx0), bkh, bkw = res["crop"]
        ry0, ry1 = max(0, by0 - cy0), min(kh, by0 + bkh - cy0)
        rx0, rx1 = max(0, bx0 - cx0), min(kw, bx0 + bkw - cx0)
        valid = np.zeros((kh, kw), bool)
        if ry1 > ry0 and rx1 > rx0:
            valid[ry0:ry1, rx0:rx1] = True
        w[~valid] = 0.0
        hp = os.path.join(bdir, f"hdr_{b + 1}.npy")
        wp = os.path.join(bdir, f"w_{b + 1}.npy")
        np.save(hp, np.asarray(h, np.float32))
        np.save(wp, w)
        off = res.get("offset") or {}
        rec.update({"offset_px": [round(off.get("dy", 0.0), 2), round(off.get("dx", 0.0), 2)],
                    "offset_tiers": off.get("tiers"), "offset_scatter_px": off.get("scatter_px"),
                    "valid_frac": round(float(valid.mean()), 4),
                    "seconds": round(time.time() - t0, 1),
                    "photometry": res.get("photometry")})
        items.append({"b": b + 1, "h": hp, "w": wp, "rec": rec})
        info["brackets"].append(rec)
        progress.log(f"bracket {b + 1}/{len(order)}: {len(ps)} frames, merged and placed "
                     f"({off.get('dy', 0.0):+.1f}, {off.get('dx', 0.0):+.1f}) px on this frame "
                     f"from {off.get('tiers')} tier(s), scatter {off.get('scatter_px', 0):.2f} px, "
                     f"{time.time() - t0:.0f} s", None)
        del res, h, w, valid
    if len(items) < 2:
        progress.log("[warn] HDR per bracket: fewer than two bracket HDRs, "
                     "the per-exposure merge is kept", None)
        info["why"] = "fewer than two bracket HDRs"
        shutil.rmtree(bdir, ignore_errors=True)
        return hdr

    # ---- the Moon of each HDR, and of the main merge (the target) ----
    Lm = _half_lum(hdr)
    tgt = _moon(Lm)
    _Xm = None
    for it in items:
        Lb = _half_lum(np.load(it["h"], mmap_mode="r"))
        it["moon"] = _moon(Lb, (tgt[0] / 2, tgt[1] / 2, tgt[2] / 2) if tgt else None)
        # ---- transparency and level: this HDR against the main merge ----
        # One gain (the bracket's transparency) from the brighter corona,
        # where an offset does not count, and one offset per channel from the
        # faintest part, where it is all that counts: each bracket's own
        # pedestal fit comes from four to fourteen single frames, and the sky
        # and haze differ from bracket to bracket. Measured on the synthetic
        # bracket: without the offset the brackets disagree in the outer field
        # by about 5 noise sigma.
        wb = np.load(it["w"], mmap_mode="r")[:2 * Lb.shape[0]:2, :2 * Lb.shape[1]:2]
        g, off = None, np.zeros(3, np.float32)
        if tgt is not None and it["moon"] is not None:
            yy, xx = np.ogrid[0:Lb.shape[0], 0:Lb.shape[1]]
            rt = np.hypot(yy - tgt[0] / 2, xx - tgt[1] / 2) / (tgt[2] / 2)
            rb = np.hypot(yy - it["moon"][0] / 2, xx - it["moon"][1] / 2) / (it["moon"][2] / 2)
            m = (rt > 1.3) & (rb > 1.3) & (wb > 0)
            if m.sum() > 4000:
                if _Xm is None:
                    _Xm = _half_rgb(hdr)
                Xm_m = _Xm[m]
                Xb = _half_rgb(np.load(it["h"], mmap_mode="r"))[m]
                lm = Lm[m]
                bright = lm > np.percentile(lm, 60)
                faint = lm < np.percentile(lm, 30)
                lb = 0.2126 * Xb[:, 0] + 0.7152 * Xb[:, 1] + 0.0722 * Xb[:, 2]
                ol = 0.0
                for _ in range(3):
                    g = float(np.median((lb[bright] - ol) / np.maximum(lm[bright], 1e-30)))
                    off = np.array([np.median(Xb[faint, c] - g * Xm_m[faint, c]) for c in range(3)],
                                   np.float32)
                    ol = float(0.2126 * off[0] + 0.7152 * off[1] + 0.0722 * off[2])
                del Xb, Xm_m
        it["gain"] = g if (g is not None and np.isfinite(g) and g > 0) else 1.0
        it["offset"] = off if np.all(np.isfinite(off)) else np.zeros(3, np.float32)
        it["rec"]["level_offset"] = [float(v) for v in it["offset"]]
        it["rec"]["transparency"] = round(it["gain"], 4)
        if it["moon"] is not None and tgt is not None:
            it["rec"]["moon_from_target_px"] = round(float(np.hypot(it["moon"][0] - tgt[0],
                                                                    it["moon"][1] - tgt[1])), 2)
        del Lb, wb
    if tgt is None:
        progress.log("[warn] HDR per bracket: no lunar limb on the main merge; "
                     "brackets combined without the Moon's weight", None)
    gl = ", ".join(f"{it['b']}: {it['gain']:.3f}" for it in items)
    progress.log(f"bracket transparency against the merge of all frames: {gl}", None)
    _ol = [float(0.2126 * it["offset"][0] + 0.7152 * it["offset"][1] + 0.0722 * it["offset"][2])
           for it in items]
    _sky = float(np.percentile(Lm[Lm > 0], 5)) if (Lm > 0).any() else 0.0
    progress.log("bracket level offsets against the merge of all frames, in its faint field's "
                 "5th-percentile units: " + ", ".join(f"{it['b']}: {o / max(_sky, 1e-30):+.2f}"
                                                      for it, o in zip(items, _ol)), None)
    ml = ", ".join(f"{it['b']}: {it['rec'].get('moon_from_target_px', float('nan')):.1f}"
                   for it in items)
    progress.log(f"bracket Moons from the target Moon (px): {ml}", None)

    # ---- the combine, by rows: weights, rejection, mean ----
    n = len(items)
    H, W = kh, kw
    Rt = tgt[2] if tgt else 0.0
    margin, soft = max(8.0, 0.02 * Rt), 2.0
    hm = [np.load(it["h"], mmap_mode="r") for it in items]
    wm = [np.load(it["w"], mmap_mode="r") for it in items]
    gains = np.array([it["gain"] for it in items], np.float32)
    offs = np.stack([it["offset"] for it in items]).astype(np.float32)

    def _block(r0, r1):
        X = ((np.stack([np.asarray(h[r0:r1], np.float32) for h in hm]) - offs[:, None, None, :])
             / gains[:, None, None, None])
        Wt = np.stack([np.asarray(w[r0:r1], np.float32) for w in wm])
        if tgt is not None:
            MW = np.stack([_moon_w_rows(r0, r1, W, it["moon"][0], it["moon"][1], tgt[0], tgt[1],
                                        Rt, margin, soft) if it["moon"] is not None
                           else np.ones((r1 - r0, W), np.float32) for it in items])
        else:
            MW = np.ones_like(Wt)
        L = 0.2126 * X[..., 0] + 0.7152 * X[..., 1] + 0.0722 * X[..., 2]
        return X, Wt, MW, L

    # THE NOISE MODEL FOR THE REJECTION, from the data (as _clipped_mean does
    # for the frames of one tier): the deviation of each bracket from the
    # median of the brackets, scaled to a common merge weight, has a spread
    # sigma^2(m) = a m + b + (e m)^2 in the median m -- photon noise, a floor
    # (read noise, sky), and a part proportional to the signal: the brackets
    # are separate merges, with their own ladder, placement and resampling,
    # and in the bright inner corona those differ by more than the photon
    # noise (measured on the synthetic bracket: 0.1 px of placement is 1 %
    # at 1.1 R). Fitted on binned robust spreads over a sample of rows. Only
    # pixels clear of every bracket's limb band are tested: there the
    # brackets differ because the Moon moved, which is the Moon's weight's
    # business, not an outlier. A quantile guard, as in _clipped_mean, keeps
    # a model that is still too tight from rejecting more than kappa implies.
    Wref = float(np.median([np.median(np.asarray(w[::7, ::7])[np.asarray(w[::7, ::7]) > 0])
                            for w in wm]))

    def _test_stats(X, Wt, MW, L):
        clear = (MW >= 0.9999).all(axis=0)
        ok = (Wt > 0)
        nv = ok.sum(axis=0)
        Lq = np.where(ok, L, np.nan)
        with np.errstate(all="ignore"):
            med = np.nanmedian(Lq, axis=0)
        test = clear & (nv >= 3) & np.isfinite(med)
        dev = (L - np.nan_to_num(med)[None]) * np.sqrt(np.maximum(Wt, 0) / Wref)
        return ok, test, np.nan_to_num(med), dev

    meds, devs = [], []
    for r0 in range(0, H, max(BRACKET_BLOCK, H // 12)):
        r1 = min(H, r0 + 32)
        X, Wt, MW, L = _block(r0, r1)
        ok, test, med, dev = _test_stats(X, Wt, MW, L)
        sel = ok & test[None]
        if sel.any():
            meds.append(np.broadcast_to(med[None], L.shape)[sel])
            devs.append(np.abs(dev[sel]))
        del X, Wt, MW, L
    model = None
    if meds and sum(m.size for m in meds) > 2000 and n >= 3:
        mm = np.concatenate(meds); dd = np.concatenate(devs)
        lo, hi = np.percentile(mm, 1), np.percentile(mm, 99.9)
        edges = np.unique(np.percentile(mm, np.linspace(1, 99.9, 25)))
        xs, ys = [], []
        for i in range(len(edges) - 1):
            sel = (mm >= edges[i]) & (mm < edges[i + 1])
            if sel.sum() < 64:
                continue
            xs.append(float(np.median(mm[sel])))
            ys.append(float((1.4826 * np.median(dd[sel])) ** 2))
        if len(xs) >= 3:
            xs = np.asarray(xs); ys = np.asarray(ys)
            A = np.stack([xs, np.ones_like(xs), xs * xs], 1)
            # relative fit: every bin counts alike, bright and faint
            coef = np.linalg.lstsq(A / ys[:, None], np.ones_like(ys), rcond=None)[0]
            coef = np.maximum(coef, 0.0)
            if coef[1] <= 0:
                coef[1] = float(np.min(ys)) * 0.5
            model = [float(c) for c in coef]
    from math import erfc, sqrt as _sqrt
    _p = erfc(BRACKET_KAPPA / _sqrt(2.0))
    _infl = float(np.sqrt(1.0 + 1.0 / max(n - 1, 1)))

    def _tval(med, dev):
        a_, b_, c_ = model
        sig = np.sqrt(a_ * np.maximum(med, 0) + b_ + c_ * med * med)[None] * _infl
        return np.abs(dev) / np.maximum(sig, 1e-30)

    cut = BRACKET_KAPPA
    if model is not None:
        ts = []
        for r0 in range(0, H, max(BRACKET_BLOCK, H // 12)):
            r1 = min(H, r0 + 32)
            X, Wt, MW, L = _block(r0, r1)
            ok, test, med, dev = _test_stats(X, Wt, MW, L)
            sel = ok & test[None]
            if sel.any():
                ts.append(_tval(med, dev)[sel])
            del X, Wt, MW, L
        if ts:
            q = float(np.quantile(np.concatenate(ts), max(0.0, 1.0 - _p)))
            cut = max(BRACKET_KAPPA, q)
    info["noise_model"] = {"a": model[0], "b": model[1], "e": float(np.sqrt(model[2]))} if model else None
    info["cut_sigma"] = round(cut, 3)

    out = np.empty((H, W, 3), np.float32)
    _dbg = os.environ.get("ECLIPSEFORGE_BRACKET_DEBUG", "") not in ("", "0")
    rejmap = np.zeros((H, W), np.uint8) if _dbg else None
    tmap = np.zeros((H, W), np.float32) if _dbg else None
    nrej = 0
    nsamp = 0
    for r0 in range(0, H, BRACKET_BLOCK):
        r1 = min(H, r0 + BRACKET_BLOCK)
        X, Wt, MW, L = _block(r0, r1)
        wgt = Wt * MW
        keep = Wt > 0
        if model is not None:
            ok, test, med, dev = _test_stats(X, Wt, MW, L)
            t = _tval(med, dev)
            rej = ok & test[None] & (t > cut)
            # never reject every bracket at a pixel
            allr = rej.sum(axis=0) >= ok.sum(axis=0)
            rej &= ~allr[None]
            keep &= ~rej
            if _dbg:
                rejmap[r0:r1] = rej.sum(axis=0)
                tmap[r0:r1] = np.where(ok & test[None], t, 0).max(axis=0)
            nrej += int(rej.sum())
            nsamp += int((ok & test[None]).sum())
            del med, dev, t, rej
        ww = wgt * keep
        s = ww.sum(axis=0)
        acc = (ww[..., None] * X).sum(axis=0)
        blk = acc / np.maximum(s, 1e-30)[..., None]
        nodat = s <= 0
        if nodat.any():
            blk[nodat] = hdr[r0:r1][nodat]
        out[r0:r1] = blk
        del X, Wt, MW, L, wgt, ok, keep, ww, s, acc, blk
    info.update({"used": True, "n": n, "rejected_frac": round(nrej / max(nsamp, 1), 6),
                 "kappa": BRACKET_KAPPA,
                 "target_moon": [round(v, 2) for v in tgt] if tgt else None,
                 "seconds": round(time.time() - t_start, 1)})
    progress.log(f"[ok] HDR per bracket: {n} bracket HDRs combined; "
                 f"{100.0 * nrej / max(nsamp, 1):.3f}% of samples rejected "
                 f"(> {cut:.1f} sigma from the median of the brackets, outside the limb band); "
                 f"{time.time() - t_start:.0f} s", None)
    del hm, wm
    if _dbg:
        np.save(os.path.join(wd, "bracket_rej.npy"), rejmap)
        np.save(os.path.join(wd, "bracket_tmax.npy"), tmap)
    else:
        shutil.rmtree(bdir, ignore_errors=True)
    return out
