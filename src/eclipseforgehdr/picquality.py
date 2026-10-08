"""Picture quality: graded measurements of the FINISHED picture (lab 0.24).

Everything here is measured on the display-encoded composite as it is written,
scale-free (radii in R), so the same thresholds hold for any set and size. Each
check returns a value, a grade (ok / check / bad) and one line saying what it
means and what to reach for. Thresholds were set by running the same code on
published renders (Druckmüller 2026, Lefaudeux 2017/2024, Project Helion) and on
our own renders with known faults.

The list follows what Project Helion's grading names (ring_break, inner_halo,
grain_level, corona_acutance, limb_width, colour_step, corona_neutrality,
limb_blue, clipped_black/white, posterization, flat_blocks), each re-derived
here from first principles rather than copied.
"""
import numpy as np
from scipy import ndimage

OK, CHECK, BAD = "ok", "check", "bad"


def _Y(a):
    return 0.2126 * a[..., 0] + 0.7152 * a[..., 1] + 0.0722 * a[..., 2]


def _grade(v, lo, hi):
    """lo/hi: the ok and bad thresholds of a 'lower is better' value."""
    return OK if v <= lo else (BAD if v >= hi else CHECK)


def _profile(v, rr, edges, w=None):
    """per radial bin: median of v (sampled), nan where too few pixels"""
    idx = np.digitize(rr, edges) - 1
    out = np.full(len(edges) - 1, np.nan)
    ok = (idx >= 0) & (idx < len(edges) - 1)
    if w is not None:
        ok &= w
    iv, vv = idx[ok], v[ok]
    order = np.argsort(iv, kind="stable")
    iv, vv = iv[order], vv[order]
    cuts = np.searchsorted(iv, np.arange(len(edges)))
    for k in range(len(edges) - 1):
        seg = vv[cuts[k]:cuts[k + 1]]
        if seg.size >= 30:
            out[k] = float(np.median(seg))
    return out


def grade(rgb, cy, cx, R, Rmask=None):
    """rgb: HxWx3 float 0..1, display-encoded. Returns {key: {value, grade, unit, note}}."""
    rgb = np.clip(np.asarray(rgb, np.float32), 0, 1)
    H, W = rgb.shape[:2]
    Rmask = float(Rmask or R)
    # work on a copy binned so R is ~300-600 px: plenty for every radial check
    b = max(1, int(round(R / 300.0)))
    if b > 1:
        H2, W2 = H // b * b, W // b * b
        sm = rgb[:H2, :W2].reshape(H2 // b, b, W2 // b, b, 3).mean(axis=(1, 3))
    else:
        sm = rgb
    cyb, cxb, Rb, Rmb = cy / b, cx / b, R / b, Rmask / b
    h, w = sm.shape[:2]
    yy = np.arange(h, dtype=np.float32)[:, None] - cyb
    xx = np.arange(w, dtype=np.float32)[None, :] - cxb
    rr = np.sqrt(yy * yy + xx * xx) / Rb
    Y = _Y(sm)
    out = {}
    def put(k, v, g, unit, note):
        out[k] = {"value": round(float(v), 4), "grade": g, "unit": unit, "note": note}

    rmax = float(min(np.max(rr), 3.4))
    # THE LAST FULL CIRCLE (lab 0.24): beyond the nearest frame edge a ring is
    # only part of a circle, and the azimuthal median jumps where a side of the
    # corona drops out -- on Val Italo's FITS that read as a 4 % "ring" at
    # 2.87 R, exactly the bottom edge. The ring and colour-step profiles stop there.
    rfull = float(min(cyb, h - cyb, cxb, w - cxb) / Rb) - 0.03
    rprof = min(rmax, rfull)
    # ---- rings: bumps in the azimuthal median of log Y against its own smooth fall-off.
    # From 1.2 R (the limb peak and its fall-off are no ring), local cubic fit over
    # 0.12 R as the trend, per-bin noise taken out with a 2-bin smooth.
    step = 1.0 / Rb
    edges = np.arange(1.15, rprof, max(step, 0.004))
    if edges.size > 60:
        from scipy.signal import savgol_filter
        L = np.log(np.maximum(Y, 1e-4))
        p = _profile(L, rr, edges)
        good = np.isfinite(p)
        if good.sum() > 60:
            x = 0.5 * (edges[1:] + edges[:-1])
            pf = np.interp(x, x[good], p[good])
            win = int(0.12 / (edges[1] - edges[0])) | 1
            res = ndimage.gaussian_filter1d(pf - savgol_filter(pf, max(win, 7), 3, mode="interp"), 2.0)
            sel = x > 1.2
            amp = float(np.max(np.abs(res[sel]))) if sel.any() else 0.0
            at = float(x[sel][np.argmax(np.abs(res[sel]))]) if sel.any() else 0.0
            put("rings", 100 * amp, _grade(100 * amp, 0.8, 2.0), "%",
                "concentric bands: largest departure of the azimuthal median from its own "
                "smooth fall-off beyond 1.2 R, at %.2f R (a ring or a tier seam shows here)" % at)
    # ---- dark rim / halo at the limb, per direction: beyond the limb peak the corona
    # should only fall; a darker band with brighter corona further out is a rim
    # (prominence feet, a cut-out, a merge edge). Prominence pixels (red) left out.
    nth = 360
    r0 = Rmb / Rb + 1.0 / Rb
    rs = np.arange(r0, 1.15, 0.5 / Rb)
    if rs.size > 10:
        th = np.radians(np.arange(nth) + 0.5)
        py = cyb + np.outer(np.sin(th), rs * Rb); px = cxb + np.outer(np.cos(th), rs * Rb)
        P = ndimage.map_coordinates(ndimage.gaussian_filter(Y, 1.0), [py, px], order=1, mode="nearest")
        Gp = ndimage.map_coordinates(sm[..., 1], [py, px], order=1, mode="nearest")
        Rp = ndimage.map_coordinates(sm[..., 0], [py, px], order=1, mode="nearest")
        red = Rp > 1.25 * np.maximum(Gp, 1e-4)
        # 7 degree running median across directions: the dark lanes between single
        # streamers are structure, a rim along a stretch of limb is not
        Pn = np.where(red, np.nan, P)
        with np.errstate(all="ignore"):
            import warnings
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                P = np.nanmedian(np.stack([np.roll(Pn, k, axis=0) for k in range(-3, 4)]), axis=0)
        P = np.where(np.isfinite(P), P, 0.0)
        red = ndimage.maximum_filter(red, size=(7, 1), mode="wrap")
        dips = []
        kmax = max(3, int(0.06 / (rs[1] - rs[0])))
        for t in range(nth):
            if red[t].mean() > 0.05:
                continue
            p_ = P[t]; k = int(np.argmax(p_[:kmax]))
            tail = p_[k:]
            if tail.size < 4:
                continue
            fmax = np.maximum.accumulate(tail[::-1])[::-1]
            dips.append(float(np.max((fmax - tail) / np.maximum(fmax, 1e-3))))
        if len(dips) > 30:
            d95 = 100 * float(np.percentile(dips, 95)); share = 100 * float(np.mean(np.array(dips) > 0.05))
            put("limb_dark_rim", d95, _grade(d95, 6.0, 15.0), "%",
                "dark band just outside the limb with brighter corona beyond it, worst 5 %% of "
                "directions (%.0f %% of directions dip more than 5 %%)" % share)
    # ---- limb width: 10-90 % rise from the disc level to the corona peak, in px of the full image
    e3 = np.arange(0.9, 1.15, 0.5 / Rb)
    p3 = _profile(Y, rr, e3)
    if np.isfinite(p3).sum() > 10:
        x3 = 0.5 * (e3[1:] + e3[:-1]); g3 = np.isfinite(p3)
        x3, p3 = x3[g3], p3[g3]
        k = int(np.argmax(p3)); lo_, hi_ = float(np.min(p3[:max(k, 1)])), float(p3[k])
        if hi_ - lo_ > 0.05 and k > 2:
            seg = p3[:k + 1]
            t10 = np.interp(lo_ + 0.1 * (hi_ - lo_), seg, x3[:k + 1])
            t90 = np.interp(lo_ + 0.9 * (hi_ - lo_), seg, x3[:k + 1])
            wpx = (t90 - t10) * R
            put("limb_width", (t90 - t10) * 1000, _grade((t90 - t10) * 1000, 30.0, 50.0), "mR",
                "the Moon's edge, 10-90 %% rise: %.1f px at this size (Disc edge softness "
                "and the limb itself)" % wpx)
    # ---- colour: neutrality of the inner corona, steps, limb fringe
    ann = (rr > 1.05) & (rr < 1.6)
    if ann.sum() > 200:
        v = sm[ann]; yv = _Y(v); sel = (yv > np.percentile(yv, 10)) & (v[:, 0] < 1.25 * np.maximum(v[:, 1], 1e-4))
        c = np.median(v[sel], 0); rg, bg = c[0] / max(c[1], 1e-6), c[2] / max(c[1], 1e-6)
        dev = 100 * max(abs(np.log(rg)), abs(np.log(bg)))
        put("corona_neutral", dev, _grade(dev, 6.0, 12.0), "%",
            "inner corona (1.05-1.6 R) colour off white: R/G %.3f, B/G %.3f "
            "(a preset or Pick white point sets it)" % (rg, bg))
    e4 = np.arange(1.05, rprof, 0.02)
    if e4.size > 10:
        notprom = sm[..., 0] < 1.25 * np.maximum(sm[..., 1], 1e-4)
        lg = np.log(np.maximum(sm[..., 1], 1e-4))
        prg = _profile(np.log(np.maximum(sm[..., 0], 1e-4)) - lg, rr, e4, notprom)
        pbg = _profile(np.log(np.maximum(sm[..., 2], 1e-4)) - lg, rr, e4, notprom)
        g4 = np.isfinite(prg) & np.isfinite(pbg)
        if g4.sum() > 8:
            # a step is a change over one 0.02 R bin that is not part of the smooth trend
            def stp(p):
                p = p[g4]; return np.abs(np.diff(p) - ndimage.median_filter(np.diff(p), 5, mode="nearest"))
            s = 100 * float(max(stp(prg).max(), stp(pbg).max()))
            put("colour_step", s, _grade(s, 2.0, 5.0), "%",
                "largest sudden colour change between neighbouring 0.02 R rings (a seam between "
                "tiers or a correction that ends at a radius)")
    rim = (rr > Rmb / Rb + 1.0 / Rb) & (rr < Rmb / Rb + 1.0 / Rb + 0.025)
    inner = (rr > 1.06) & (rr < 1.15)
    if rim.sum() > 100 and inner.sum() > 100:
        def ch(m):
            v = sm[m]; v = v[v[:, 0] < 1.25 * np.maximum(v[:, 1], 1e-4)]
            c = np.median(v, 0) if len(v) > 50 else np.array([1, 1, 1.0]); return c / max(c[1], 1e-6)
        a_, b_ = ch(rim), ch(inner)
        db = 100 * float(np.log(a_[2] / b_[2]))
        put("limb_blue", db, _grade(abs(db), 4.0, 10.0), "%",
            "blue at the limb against the inner corona (+ = bluer rim: a colour-plane offset "
            "or the disc colour bleeding out)")
    # ---- grain: pixel noise measured ALONG THE RADIUS. Coronal structure is
    # mostly radial, so the radial derivative carries the noise and little of the
    # fine structure (a plain high-pass rates Lefaudeux's hair-fine streamers as
    # grain). Central difference minus its local mean, robust sigma, x sqrt(2) for
    # the difference, in 8-bit levels; per zone.
    zones = [("1.2-1.6 R", 1.2, 1.6), ("1.6-2.4 R", 1.6, 2.4), ("2.4-3.2 R", 2.4, 3.2)]
    # on the binned copy (R ~ 300 px): grain as it reads at a normal viewing size,
    # the same for a 24 Mpx export and a web-size reference
    sf = b
    yf = Y
    gy = np.zeros_like(yf); gx = np.zeros_like(yf)
    gy[1:-1] = 0.5 * (yf[2:] - yf[:-2]); gx[:, 1:-1] = 0.5 * (yf[:, 2:] - yf[:, :-2])
    yyf = (np.arange(yf.shape[0], dtype=np.float32)[:, None] * sf - cy)
    xxf = (np.arange(yf.shape[1], dtype=np.float32)[None, :] * sf - cx)
    rrf = np.sqrt(yyf * yyf + xxf * xxf)
    gr = (gy * yyf + gx * xxf) / np.maximum(rrf, 1.0)
    gr = gr - ndimage.uniform_filter(gr, 7)
    rrf /= R
    gz = {}
    for nm, a, b2 in zones:
        m = (rrf > a) & (rrf < b2)
        if m.sum() > 2000:
            v = gr[m][::5]
            gz[nm] = 255 * np.sqrt(2.0) * 1.4826 * float(np.median(np.abs(v - np.median(v))))
    del gy, gx, gr, yyf, xxf
    if gz:
        gmax = max(gz.values())
        put("grain", gmax, _grade(gmax, 2.5, 5.0), "levels",
            "pixel noise along the radius at R = 300 px, robust sigma in 8-bit levels, worst zone (" +
            ", ".join("%s %.1f" % kv for kv in gz.items()) + ")")
        if len(gz) > 1:
            ratio = max(gz.values()) / max(min(gz.values()), 1e-3)
            put("grain_gap", ratio, _grade(ratio, 2.5, 5.0), "x",
                "grain of the noisiest zone against the calmest (the texture changes with radius)")
    # ---- detail (acutance): fine structure at 1.3-2 R, log Y high-pass at R/150 scales
    s = Rb / 150.0
    Lg = np.log(np.maximum(Y, 1e-3))
    hpd = ndimage.gaussian_filter(Lg, 1 * s) - ndimage.gaussian_filter(Lg, 8 * s)
    m = (rr > 1.3) & (rr < 2.0)
    if m.sum() > 500:
        d = float(np.std(hpd[m]))
        put("detail", d, OK if d >= 0.035 else CHECK, "",
            "fine-structure contrast at 1.3-2 R (published renders: 0.05-0.07; below 0.035 "
            "the corona reads soft)")
        # flat blocks: tiles of ~R/40 with almost no structure, inside 1.1-2.4 R
        t = max(4, int(Rb / 40))
        hh, ww = h // t * t, w // t * t
        hb = np.abs(hpd[:hh, :ww]).reshape(hh // t, t, ww // t, t).mean(axis=(1, 3))
        rt = rr[:hh, :ww].reshape(hh // t, t, ww // t, t).mean(axis=(1, 3))
        zt = (rt > 1.1) & (rt < 2.4)
        if zt.sum() > 50:
            ref = float(np.median(hb[zt]))
            fb = 100 * float(np.mean(hb[zt] < 0.15 * ref))
            put("flat_blocks", fb, _grade(fb, 1.0, 5.0), "%",
                "share of the 1.1-2.4 R corona in patches with no structure at all (filled or "
                "smoothed-over areas)")
    # ---- clipping, outside the Moon
    outm = rr > Rmb / Rb + 2.0 / Rb
    if outm.sum() > 1000:
        v = rgb[::b, ::b][:h, :w][outm]
        cw = 100 * float(np.mean(v.max(1) >= 254.5 / 255))
        ck = 100 * float(np.mean(v.max(1) <= 0.5 / 255))
        put("clipped_white", cw, _grade(cw, 0.5, 2.0), "%",
            "corona at full white in any channel (lose structure there: Contain whites, "
            "Inner corona level or Highlight compression)")
        put("clipped_black", ck, _grade(ck, 0.5, 5.0), "%",
            "outside the Moon at pure black (Black point too high, or the sky lift off)")
    # ---- posterisation in the smooth outer field: empty histogram bins at 12-bit
    far = rr > min(2.6, rmax - 0.2)
    if far.sum() > 2000:
        q = np.round(_Y(rgb[::b, ::b][:h, :w][far]) * 4095).astype(np.int32)
        lo, hi = np.percentile(q, [5, 95]).astype(int)
        if hi - lo >= 8:
            qq = q[(q >= lo) & (q <= hi)] - lo
            hist = np.bincount(qq, minlength=hi - lo + 1)[:hi - lo + 1]
            empty = 100 * float(np.mean(hist == 0))
            put("posterisation", empty, _grade(empty, 10.0, 30.0), "%",
                "empty levels in the outer field's histogram at 12 bits (banding)")
    return out


def summary(q):
    bad = [k for k, v in q.items() if v["grade"] == BAD]
    chk = [k for k, v in q.items() if v["grade"] == CHECK]
    if not bad and not chk:
        return "picture quality: all %d checks ok" % len(q)
    s = "picture quality: "
    if bad:
        s += "BAD " + ", ".join(bad)
    if chk:
        s += ("; " if bad else "") + "check " + ", ".join(chk)
    return s


def text(q):
    lines = ["PICTURE QUALITY (measured on the exported picture)"]
    for k, v in q.items():
        lines.append("  %-15s %-6s %9.3f %-6s %s" % (k, v["grade"].upper(), v["value"], v["unit"], v["note"]))
    return "\n".join(lines)
