"""Picture quality: graded measurements of the FINISHED picture (lab 0.24).

Everything here is measured on the display-encoded composite as it is written,
scale-free (radii in R), so the same thresholds hold for any set and size. Each
check returns a value, a grade (ok / check / bad) and one line saying what it
means and what to reach for. Thresholds were set by running the same code on
published renders and on
our own renders with known faults.

The list follows common grading names (ring_break, inner_halo,
grain_level, corona_acutance, limb_width, colour_step, corona_neutrality,
limb_blue, clipped_black/white, posterization, flat_blocks), each re-derived
here from first principles rather than copied.
"""
import numpy as np
from scipy import ndimage
from .parallel import gaussian_filter as _pgauss   # scipy's, on several threads

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
    # corona drops out -- on a FITS test set that read as a 4 % "ring" at
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
            out["_rings_at"] = at
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
        P = ndimage.map_coordinates(_pgauss(Y, 1.0), [py, px], order=1, mode="nearest")
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
    # fine structure (a plain high-pass rates a published render's hair-fine streamers as
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
        out["_grain_zone"] = max(gz, key=gz.get)
        put("grain", gmax, _grade(gmax, 2.5, 5.0), "levels",
            "pixel noise along the radius at R = 300 px, robust sigma in 8-bit levels, worst zone (" +
            ", ".join("%s %.1f" % kv for kv in gz.items()) + "). Measured before any denoising: "
            "Grain smoothing here, or a denoiser after export, brings it down")
        if len(gz) > 1:
            ratio = max(gz.values()) / max(min(gz.values()), 1e-3)
            put("grain_gap", ratio, _grade(ratio, 2.5, 5.0), "x",
                "grain of the noisiest zone against the calmest (the texture changes with radius)")
    # ---- detail (acutance): fine structure at 1.3-2 R, log Y high-pass at R/150 scales
    s = Rb / 150.0
    Lg = np.log(np.maximum(Y, 1e-3))
    hpd = _pgauss(Lg, 1 * s) - _pgauss(Lg, 8 * s)
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


# OUR OWN RATING (2026-10-13). Thirteen measurements averaged told a user
# little: clipping and banding are nearly always excellent and lifted the mean,
# and "weakest grain (weak)" said neither what it means nor what to do. The
# checks are now read as five things a viewer sees -- detail, noise,
# artefacts, colour, tones -- each graded by its worst check on seven steps
# (excellent to 0.5x the reference, very good 0.75x, good 1x, fair 1.5x,
# weak 2.5x, poor 4x, then very poor; a common grading scale). The picture's grade is their average, at most two steps above the
# weakest. Tones count only when they fail. Each category says whether what
# limits it is the FRAMES (detail, noise: light, focal length, seeing) or the
# PROCESSING (artefacts, colour, tones: settings or a fault), and the
# headline gives one sentence and one thing to try for the weakest.
# References: detail 0.05 is the low end of the published renders
# (0.05-0.07) -- "good" means that, not "passes".
REF = {"rings": 0.8, "limb_dark_rim": 6.0, "limb_width": 30.0, "corona_neutral": 6.0,
       "colour_step": 2.0, "limb_blue": 4.0, "grain": 2.5, "grain_gap": 2.5,
       "detail": 0.05, "flat_blocks": 1.0, "clipped_white": 0.5, "clipped_black": 0.5,
       "posterisation": 10.0}
HIGHER_BETTER = {"detail"}
STEPS = ((0.5, "excellent"), (0.75, "very good"), (1.0, "good"), (1.5, "fair"),
         (2.5, "weak"), (4.0, "poor"))
WORDS = ["very poor", "poor", "weak", "fair", "good", "very good", "excellent"]
CATEGORIES = (
    ("detail", "Detail", "frames", ("detail", "limb_width")),
    ("noise", "Noise", "frames", ("grain", "grain_gap")),
    ("artefacts", "Artefacts", "processing", ("rings", "limb_dark_rim", "flat_blocks", "colour_step")),
    ("colour", "Colour", "processing", ("corona_neutral", "limb_blue")),
    ("tones", "Tones", "processing", ("clipped_white", "clipped_black", "posterisation")),
)


def step(key, value):
    """(score 1..7, word) of one check against its reference."""
    ref = REF.get(key)
    if ref is None:
        return None
    v = abs(float(value))
    r = (ref / max(v, 1e-9)) if key in HIGHER_BETTER else v / ref
    for k, (lim, word) in enumerate(STEPS):
        if r <= lim:
            return 7 - k, word
    return 1, "very poor"


def _advice(key, q):
    """One plain sentence and one thing to try, for the check that limits a
    category."""
    v = q[key]["value"]
    if key == "grain":
        z = q.get("_grain_zone", "")
        if z.startswith("1.2"):
            # the inner corona has light to spare: grain there is what the
            # detail layers lifted
            return ("inner corona grainy (%s): the detail settings lift grain there" % z,
                    "lower Detail gain, or Grain smoothing 0.5")
        return ("outer corona grainy (%s): the frames hold little light there" % (z or "far field"),
                "Grain smoothing 0.5, or denoise after export")
    if key == "grain_gap":
        return ("grain changes with distance from the Moon", "Grain smoothing, or denoise after export")
    if key == "detail":
        return ("fine structure is soft: what the frames resolve (focal length, seeing, "
                "number of frames)", "Detail gain or Structure lift it, and the grain with it")
    if key == "limb_width":
        return ("the Moon's edge is soft", "Disc edge softness")
    if key == "rings":
        return ("a faint ring around the Sun at %.2f R (a seam between exposures)"
                % q.get("_rings_at", 0.0), "lower Structure or Radial flatten; if it stays, report it")
    if key == "limb_dark_rim":
        return ("a dark band just outside the Moon's edge", "Inner corona level or Disc edge softness")
    if key == "flat_blocks":
        return ("patches of corona without structure (smoothed over or filled)",
                "lower Grain smoothing; if they stay, report it")
    if key == "colour_step":
        return ("the colour changes abruptly at one distance from the Moon", "if it shows, report it")
    if key == "corona_neutral":
        return ("the inner corona is not neutral white", "Pick white point, or a colour preset")
    if key == "limb_blue":
        return ("a coloured fringe at the Moon's edge", "check Colour planes / Disc edge softness")
    if key == "clipped_white":
        return ("%.1f %% of the corona is burned to white" % v,
                "Inner corona level or Highlight compression")
    if key == "clipped_black":
        return ("%.1f %% of the sky is pure black" % v, "lower Black point or raise the sky lift")
    if key == "posterisation":
        return ("the smooth outer field shows banding", "export 16-bit")
    return (key, "")


def categories(q):
    """[(id, name, source, score, word, limiting check)] for the categories
    that have measurements."""
    out = []
    for cid, name, src, keys in CATEGORIES:
        sc = [(step(k, q[k]["value"]), k) for k in keys if k in q and isinstance(q[k], dict)]
        sc = [(s, k) for s, k in sc if s is not None]
        if not sc:
            continue
        (score, word), k = min(sc, key=lambda t: t[0][0])
        if k == "grain" and str(q.get("_grain_zone", "")).startswith("1.2"):
            src = "processing"   # inner-corona grain is what the detail layers lifted
        out.append((cid, name, src, score, word, k))
    return out


def overall(q):
    """(score 1..7, word, categories). Average of the categories (tones only
    when below good), at most two steps above the weakest."""
    cats = categories(q)
    use = [c for c in cats if c[0] != "tones" or c[3] < 5]
    if not use:
        return None, "", cats
    mean = float(np.mean([c[3] for c in use]))
    sc = int(min(np.round(mean), min(c[3] for c in use) + 2))
    sc = int(np.clip(sc, 1, 7))
    return sc, WORDS[sc - 1], cats


def level(q):
    """log level of the headline: warn below fair, odd when a PROCESSING
    category is weak or worse, ok otherwise."""
    sc, _, cats = overall(q)
    if sc is None:
        return "info"
    if sc < 4:
        return "warn"
    if any(c[2] == "processing" and c[3] <= 3 for c in cats):
        return "odd"
    return "ok"


def summary(q):
    sc, word, cats = overall(q)
    if sc is None:
        return "picture quality: not graded"
    s = "picture quality: %s" % word
    weak = sorted([c for c in cats if c[3] <= 4 and (c[0] != "tones" or c[3] < 5)],
                  key=lambda c: (c[3], c[2] != "processing"))
    if not weak:
        return s + ", nothing to fix"
    c = weak[0]
    what, todo = _advice(c[5], q)
    s += ". %s%s" % (what[0].upper(), what[1:])
    if todo:
        s += " \u2192 " + todo
    return s


def text(q):
    sc, word, cats = overall(q)
    lines = ["PICTURE QUALITY (measured on the exported picture, before any denoising)"]
    if sc is not None:
        lines.append("  overall: %s (%d of 7) -- the average of the categories below, at most "
                     "two steps above the weakest" % (word, sc))
        lines.append("")
        for cid, name, src, score, w, k in cats:
            what, todo = _advice(k, q)
            limited = "set by the frames" if src == "frames" else "set by the processing"
            lines.append("  %-10s %-10s %s" % (name, w, limited))
            if score <= 4:
                lines.append("  %-10s %-10s %s%s" % ("", "", what, ("; try: " + todo) if todo else ""))
        lines += ["", "  scale: each check against its reference -- excellent to 0.5x, very good to "
                      "0.75x, good to 1x,", "  fair to 1.5x, weak to 2.5x, poor to 4x, beyond that "
                      "very poor. The category takes its worst check.", "", "  checks:"]
    for k, v in q.items():
        if not isinstance(v, dict):
            continue
        st = step(k, v["value"])
        lines.append("  %-15s %-10s %9.3f %-6s %s" % (k, st[1] if st else v["grade"],
                                                      v["value"], v["unit"], v["note"]))
    return "\n".join(lines)
