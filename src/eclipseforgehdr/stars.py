"""STAR DETECTION (lab 0.25.1, step 2 of the star layer).

Finds stars in the merge without knowing in advance where they are, and writes
them to STARS_FILE in the work directory. Nothing is drawn or rendered yet;
this is the list the star layer will be built from.

How a candidate becomes a star:
  1. The merge is binned so the lunar radius is about 300 px, which puts every
     set on the same scale whatever the focal length.
  2. A star-sized band-pass (Gaussian minus a Gaussian three times wider) at
     four widths. Each response is divided by its own local noise, measured
     robustly in 32 px blocks, so the threshold means the same thing in the
     bright inner corona and in the dark sky. The best width is kept.
  3. THE SAME SEARCH ON THE INVERTED IMAGE. Nothing real is a dark point, so
     every dark "detection" is noise or residual structure, and their count
     says how many bright ones noise alone produces. The star threshold is set
     just above the strongest dark one (never below 6 sigma).
  4. Shape: round (second-moment ellipticity < 0.4), and not concentrated in
     one pixel (a leftover hot pixel or cosmic ray puts more than 12 % of the
     central 7x7 light into its peak pixel; the stars measured hold 3-6 %).
  5. Not inside coronal structure: the band-pass in a ring 8-25 binned px
     around the candidate must look like sky (clutter ratio below 1.5 times
     the image's typical value). A compact knot in a streamer is round and
     strong but sits in clutter 2-6 times the sky's. Such candidates are
     listed as "corona", not dropped silently.
Measured on the reference stack (0.25.1): the two stars confirmed by hand
(6.9 R, z 13; 2.6 R, z 6.7) are the only two found, and a coronal knot at
1.46 R (z 7.0, clutter 6.5x, corona-coloured) is set aside. On another test
stack: the three stars marked by eye (6.8, 7.0, 6.7 R) are found, nothing
else. Cross-matched between the two sets through the anchor star (rotation
319.3 deg, scale 1.015): the shared star at 2.6 R lands 18-37 px from its
counterpart.
Candidates between z 5 and the star threshold are kept as "possible" with the
number the inverted image produced in the same range, which is the number to
expect from noise alone.
"""
import os, json
import numpy as np
from scipy import ndimage

STARS_FILE = "stars.json"
# A change here that alters WHICH stars are found must also bump detail.LAYER_BUILD:
# the stars are lifted out of the corona while the layers are built (detail.py).
STARS_BUILD = 6                # 6: two-star catalogue check with the site known, faint (4.5-5 sigma) candidates confirmable (2026-10-10). 5: sensor artefacts rejected by colour and profile (2026-10-10). 4: catalogue numbers and the main streamers for the sky map. 3: observing site and time used when given (observer.py). 2: checked against the bundled ecliptic catalogue (starcat.py)
SIGMAS = (0.6, 1.0, 1.5, 2.25)        # band-pass widths, binned px (R ~ 300 px)
RMIN = 1.4                             # inner limit, lunar radii
Z_FLOOR = 6.0
Z_POSSIBLE = 5.0
# lab 0.25.1 (2026-10-10): candidates between Z_VERIFY and Z_POSSIBLE ("faint")
# go to the catalogue check too, but only as verification: one is kept, as a
# star, only when a solution puts a catalogue star under it; it never seeds a
# solution, and unconfirmed ones are dropped. A mag-7.7 star at z 5.09 in one
# stack fell to just under 5 in the next, and the match lost its third star.
Z_VERIFY = 4.5
ELL_MAX = 0.4
PEAKFRAC_MAX = 0.12
CLUTTER_MAX = 1.5
FORCED_ADD = False             # catalogue stars missed by the search are NOT added (verify only)
# ARTEFACT TEST (2026-10-10, a test stack: a "star" with no
# catalogue star under it was a hot green photosite -- light in G only, R and
# B BELOW the background, and twice as compact as the stack's real stars).
# Starlight is broadband and spread by the optics; a sensor defect is neither.
ART_GREEN_ONLY = 0.25          # R and B both under this fraction of G ...
ART_SNR = 8.0                  # ... with G at least this significant
ART_NEG_SIGMA = 3.0            # or a channel this far below the background while another is >= ART_SNR
ART_COMPACT = 1.8              # peak/flux above this multiple of the real stars' median
ART_COMPACT_SNR = 8.0          # (only judged when the aperture flux is this significant)


def _log(progress, msg, frac=None):
    if progress is not None:
        progress.log(msg, frac)


def needs_build(wd):
    """Missing, built by an older recipe, or older than the merge."""
    try:
        sp = os.path.join(wd, STARS_FILE)
        if not os.path.exists(sp):
            return True
        if int(json.load(open(sp)).get("build", 0)) < STARS_BUILD:
            return True
        hp = os.path.join(wd, "hdr_lum.npy")
        if os.path.exists(hp) and os.path.getmtime(sp) < os.path.getmtime(hp):
            return True
        from .observer import site_path
        st = site_path(os.path.dirname(os.path.abspath(wd)))
        return os.path.exists(st) and os.path.getmtime(sp) < os.path.getmtime(st)
    except Exception:
        return True


def load(wd):
    return json.load(open(os.path.join(wd, STARS_FILE)))


def _local_rms(D, blk=32):
    H, W = D.shape
    h, w = max(1, H // blk), max(1, W // blk)
    a = np.abs(D[:h * blk, :w * blk]).reshape(h, blk, w, blk).transpose(0, 2, 1, 3).reshape(h, w, -1)
    m = ndimage.median_filter(np.median(a, axis=2) * 1.4826, 3)
    out = ndimage.zoom(m, (H / h, W / w), order=1)
    if out.shape[0] < H or out.shape[1] < W:
        out = np.pad(out, ((0, max(0, H - out.shape[0])), (0, max(0, W - out.shape[1]))), mode="edge")
    return np.maximum(out[:H, :W], 1e-12)


def _shape(img, y, x, s):
    r = int(max(4, 3 * s))
    w = img[y - r:y + r + 1, x - r:x + r + 1].astype(np.float64)
    bg = np.median(img[max(0, y - 3 * r):y + 3 * r + 1, max(0, x - 3 * r):x + 3 * r + 1])
    yy, xx = np.mgrid[-r:r + 1, -r:r + 1]
    v = np.clip(w - bg, 0, None) * np.exp(-(yy ** 2 + xx ** 2) / (2 * (1.5 * s) ** 2))
    t = v.sum() + 1e-12
    my, mx = (v * yy).sum() / t, (v * xx).sum() / t
    syy = (v * (yy - my) ** 2).sum() / t
    sxx = (v * (xx - mx) ** 2).sum() / t
    sxy = (v * (yy - my) * (xx - mx)).sum() / t
    tr, det = syy + sxx, syy * sxx - sxy ** 2
    q = np.sqrt(max(tr * tr / 4 - det, 0))
    return float(2 * q / (tr + 1e-12))


def _artefact_metrics(L, RGB, y, x, f):
    """Full-resolution look at a candidate: the peak found again within the
    bin, peak / aperture flux (r <= 4.5 px, background from 9-15 px), and the
    flux in each colour (r <= 3 px) with its significance."""
    H, W = L.shape
    if y < 24 or x < 24 or y >= H - 24 or x >= W - 24:
        return None
    w = np.asarray(L[y - f - 2:y + f + 3, x - f - 2:x + f + 3], np.float64)
    dy, dx = np.unravel_index(int(np.argmax(w)), w.shape)
    y, x = y - f - 2 + dy, x - f - 2 + dx
    P = np.asarray(L[y - 16:y + 17, x - 16:x + 17], np.float64)
    a, b = np.mgrid[-16:17, -16:17]
    r = np.hypot(a, b)
    ann = (r >= 9) & (r <= 15)
    ap = r <= 4.5
    bg = float(np.median(P[ann]))
    sd = 1.4826 * float(np.median(np.abs(P[ann] - bg))) + 1e-30
    flux = float((P[ap] - bg).sum())
    out = dict(pof=float((P[16, 16] - bg) / flux) if flux > 0 else 1.0,
               flux_snr=flux / (sd * np.sqrt(ap.sum())))
    if RGB is not None:
        Q = np.asarray(RGB[y - 16:y + 17, x - 16:x + 17], np.float64)
        ap3 = r <= 3
        fl, sn = [], []
        for k in range(3):
            q = Q[..., k]
            qb = float(np.median(q[ann]))
            qs = 1.4826 * float(np.median(np.abs(q[ann] - qb))) * np.sqrt(ap3.sum()) + 1e-30
            fk = float((q[ap3] - qb).sum())
            fl.append(fk); sn.append(fk / qs)
        out.update(rgb=fl, rgb_snr=sn)
    return out


def _colour_artefact(m):
    """Light in one colour only, or a colour below the background next to a
    strong one: not starlight. Returns the reason or None."""
    if not m or "rgb" not in m:
        return None
    R, G, B = m["rgb"]
    sR, sG, sB = m["rgb_snr"]
    if sG >= ART_SNR and R < ART_GREEN_ONLY * G and B < ART_GREEN_ONLY * G:
        return "light in green only"
    if max(m["rgb_snr"]) >= ART_SNR and min(m["rgb_snr"]) <= -ART_NEG_SIGMA:
        return "a colour below the background"
    return None


def streamers(Lb, cy, cx, R, nmax=12):
    """The main streamers, for the sketch map in the Stars group: angle
    (clockwise from up on the picture, deg) and outer reach (R). Excess of the
    log luminance over the corona's broad azimuthal shape (25 deg smoothing),
    in units of that ring's spread, averaged over 1.25-2.8 R; peaks above 1
    and 8 deg apart; the reach is the first gap in the excess along the ray.
    Measured on the reference stack: the streamers at 281, 353, 195, 113 and
    147 deg -- the left, top, bottom and the two lower-right ones in that
    picture.""" 
    H, W = Lb.shape
    rad = R * np.linspace(1.25, 6.0, 48)
    ang = np.radians(np.arange(0, 360, 0.5))
    yy = cy - np.cos(ang)[None, :] * rad[:, None]
    xx = cx + np.sin(ang)[None, :] * rad[:, None]
    ok = (yy > 2) & (yy < H - 3) & (xx > 2) & (xx < W - 3)
    v = ndimage.map_coordinates(np.log(np.maximum(ndimage.gaussian_filter(Lb, 1.5), 1e-12)), [yy.ravel(), xx.ravel()], order=1).reshape(yy.shape)
    v = np.where(ok, v, np.nan)
    # excess over the corona's broad azimuthal shape (25 deg), so rays and
    # streamers stand out of the polar-hole / equator contrast
    vf = np.where(ok, v, np.nanmedian(v, axis=1, keepdims=True))
    base = ndimage.gaussian_filter1d(vf, 25, axis=1, mode="wrap")
    d = np.where(ok, vf - base, np.nan)
    mad = 1.4826 * np.nanmedian(np.abs(d - np.nanmedian(d, axis=1, keepdims=True)), axis=1, keepdims=True) + 1e-9
    e = d / mad
    e = ndimage.gaussian_filter1d(np.nan_to_num(e, nan=0.0), 3, axis=1, mode="wrap")   # 1.5 deg
    e = np.where(ok, e, np.nan)
    inner = slice(0, 16)                                  # 1.25-2.8 R
    score = np.nanmean(e[inner], axis=0)
    score = np.nan_to_num(score, nan=-9)
    pk = (score == ndimage.maximum_filter1d(score, 17, mode="wrap")) & (score > 1.0)  # 8 deg apart
    out = []
    for j in np.flatnonzero(pk):
        col = ndimage.gaussian_filter1d(np.nan_to_num(e[:, j], nan=-9), 1.5)
        good = np.flatnonzero(col > 0.6)
        if not len(good):
            continue
        # the run from the inner edge outwards
        run = good[0]
        for k in good[1:]:
            if k == run + 1: run = k
            else: break
        r1 = rad[run] / R
        out.append(dict(angle=round(float(np.degrees(ang[j])), 1), r1=round(float(r1), 2), score=round(float(score[j]), 2)))
    out.sort(key=lambda d: -d["score"])
    return out[:nmax]


def find(wd, progress=None):
    g = json.load(open(os.path.join(wd, "geometry.json")))
    L = np.load(os.path.join(wd, "hdr_lum.npy"), mmap_mode="r")
    f = max(1, int(round(float(g["R"]) / 300.0)))
    H, W = L.shape
    H2, W2 = H // f, W // f
    Lb = np.asarray(L[:f * H2, :f * W2], np.float32).reshape(H2, f, W2, f).mean(axis=(1, 3))
    cy, cx, R = float(g["cy"]) / f, float(g["cx"]) / f, float(g["R"]) / f
    # a stack smaller than the R ~ 300 px scale (a decimated set, a short lens)
    # keeps its own pixels; every width below is scaled to it instead
    ksc = min(1.0, R / 300.0)
    sigmas = [max(0.5, s * ksc) for s in SIGMAS]
    yy, xx = np.ogrid[:H2, :W2]
    rr = np.hypot(yy - cy, xx - cx) / R
    edge = 30
    ok = rr > RMIN
    ok[:edge] = False; ok[-edge:] = False; ok[:, :edge] = False; ok[:, -edge:] = False
    Z = np.empty((len(sigmas), H2, W2), np.float32)
    for i, s in enumerate(sigmas):
        D = ndimage.gaussian_filter(Lb, s) - ndimage.gaussian_filter(Lb, 3 * s)
        Z[i] = D / _local_rms(D)
    Dc = ndimage.gaussian_filter(Lb, max(0.5, ksc)) - ndimage.gaussian_filter(Lb, 3.0 * max(0.5, ksc))
    hf = Lb - ndimage.uniform_filter(Lb, 3)

    def clutter(y, x, r0=max(4, int(round(8 * ksc))), r1=max(12, int(round(25 * ksc)))):
        if y - r1 < 0 or x - r1 < 0 or y + r1 >= H2 or x + r1 >= W2:
            return float("nan")
        w = Dc[y - r1:y + r1 + 1, x - r1:x + r1 + 1]
        h = hf[y - r1:y + r1 + 1, x - r1:x + r1 + 1]
        a, b = np.mgrid[-r1:r1 + 1, -r1:r1 + 1]
        m = np.hypot(a, b) > r0
        return float(np.std(w[m]) / (1.4826 * np.median(np.abs(h[m])) + 1e-12))

    res = {}
    Zpos = Z.max(0)
    for sign, name in ((1.0, "pos"), (-1.0, "neg")):
        Zs = sign * Z
        Zm = Zs.max(0)
        Si = Zs.argmax(0)
        pk = (Zm == ndimage.maximum_filter(Zm, max(5, int(9 * ksc) | 1))) & (Zm > Z_POSSIBLE - 0.5) & ok
        rows = []
        img = sign * Lb
        for y, x in zip(*np.nonzero(pk)):
            s = sigmas[Si[y, x]]
            if 3 * int(max(4, 3 * s)) >= min(y, x, H2 - y, W2 - x):
                continue
            Y, X = int(y) * f + f // 2, int(x) * f + f // 2
            w = sign * np.asarray(L[max(0, Y - 3):Y + 4, max(0, X - 3):X + 4], np.float64)
            bw = np.median(sign * np.asarray(L[max(0, Y - 20):Y + 21, max(0, X - 20):X + 21], np.float64))
            v = np.clip(w - bw, 0, None)
            rows.append(dict(y=Y, x=X, z=float(Zm[y, x]), sigma=float(s * f),
                             ell=_shape(img, int(y), int(x), s),
                             peakfrac=float(v.max() / (v.sum() + 1e-12)),
                             clutter=clutter(int(y), int(x)),
                             r=float(rr[y, x]),
                             angle=float(np.degrees(np.arctan2(-(y - cy), x - cx)) % 360)))
        res[name] = rows
    allc = [c["clutter"] for c in res["pos"] + res["neg"] if np.isfinite(c["clutter"])]
    cref = float(np.median(allc)) if allc else 1.0

    # a much coarser grid (a decimated stack, R well under 300 px) concentrates
    # a real star into fewer pixels too; between 0.6 and 1 the limit stays,
    # because hot-pixel residue sits at 0.12-0.15 (a test stack)
    pf_max = PEAKFRAC_MAX if ksc >= 0.6 else min(0.3, PEAKFRAC_MAX / ksc)

    def shape_ok(c):
        return c["ell"] < ELL_MAX and c["peakfrac"] < pf_max

    def sky(c):
        return np.isfinite(c["clutter"]) and c["clutter"] < CLUTTER_MAX * cref

    negz = [c["z"] for c in res["neg"] if shape_ok(c) and sky(c)]
    t = max(Z_FLOOR, (max(negz) + 0.25) if negz else Z_FLOOR)
    n_chance = sum(1 for z in negz if Z_POSSIBLE <= z < t)
    out = []
    for c in sorted(res["pos"], key=lambda c: -c["z"]):
        if c["z"] < Z_VERIFY or not shape_ok(c):
            continue
        if c["z"] < Z_POSSIBLE:
            if not sky(c):
                continue
            cls = "faint"
        elif not sky(c):
            cls = "corona"
        elif c["z"] >= t:
            cls = "star"
        else:
            cls = "possible"
        c = {k: (round(v, 3) if isinstance(v, float) else v) for k, v in c.items()}
        c["clutter"] = round(c["clutter"] / cref, 2)
        c["class"] = cls
        out.append(c)
    # ARTEFACT TEST, part 1 (colour): before the catalogue, so a sensor defect
    # takes no part in the match
    try:
        _rp = os.path.join(wd, "hdr_rgb.npy")
        RGB = np.load(_rp, mmap_mode="r") if os.path.exists(_rp) else None
    except Exception:
        RGB = None
    for c in out:
        m = _artefact_metrics(L, RGB, c["y"], c["x"], f)
        if m is None:
            continue
        c["pof"] = round(m["pof"], 4)
        c["flux_snr"] = round(m["flux_snr"], 1)
        if "rgb_snr" in m:
            c["rgb_snr"] = [round(v, 1) for v in m["rgb_snr"]]
        why = _colour_artefact(m) if c["class"] in ("star", "possible", "faint") else None
        if why:
            c["class"] = "artefact"; c["why"] = why
    for c in out:
        if c["class"] == "faint":
            c["faint"] = True
    cat = _catalogue_check(wd, g, out, Zpos, f, (H, W), ok, progress)
    # faint candidates the catalogue did not confirm are not kept
    out = [c for c in out if c["class"] != "faint"]
    for c in out:
        c.pop("faint", None)
    # ARTEFACT TEST, part 2 (profile): sharper than the stack's own stars --
    # the catalogue-confirmed ones, else the strongest stars found -- is a
    # sensor defect, not light through the lens. Confirmed stars are kept.
    refs = [c["pof"] for c in out if c.get("catalogue") and "pof" in c and c["class"] == "star"]
    if len(refs) < 2:
        refs = [c["pof"] for c in sorted(out, key=lambda c: -c["z"])
                if c["class"] == "star" and "pof" in c][:5]
    if len(refs) >= 2:
        pref = float(np.median(refs))
        for c in out:
            if (c["class"] in ("star", "possible") and not c.get("catalogue") and "pof" in c
                    and c.get("flux_snr", 0) >= ART_COMPACT_SNR and c["pof"] > ART_COMPACT * pref):
                c["class"] = "artefact"
                c["why"] = f"sharper than the stars ({c['pof'] / pref:.1f}x their peak/flux)"
    n_art = sum(1 for c in out if c["class"] == "artefact")
    if n_art:
        _log(progress, f"stars: {n_art} candidate(s) rejected as sensor artefacts (" +
             "; ".join(sorted({c['why'].split(' (')[0] for c in out if c['class'] == 'artefact'})) + ")", None)
    # streamers are not drawn on the sky map any more (2026-10-10): the
    # luminance-based finder did not match the finished picture by eye
    strm = []
    meta = dict(build=STARS_BUILD, bin=f, threshold=round(t, 2), possible_from=Z_POSSIBLE,
                expected_by_chance=n_chance, rmin=RMIN,
                n_star=sum(1 for c in out if c["class"] == "star"),
                n_possible=sum(1 for c in out if c["class"] == "possible"),
                n_corona=sum(1 for c in out if c["class"] == "corona"),
                n_artefact=sum(1 for c in out if c["class"] == "artefact"),
                catalogue=cat, streamers=strm, shape=[int(H), int(W)],
                moon=[float(g["cy"]), float(g["cx"]), float(g["R"])], stars=out)
    return meta


def _catalogue_check(wd, g, out, Zpos, f, shape, sky_ok, progress=None):
    """Match the candidates to the bundled ecliptic catalogue (starcat.py) and
    act on a clear solution, in place on `out`:
      * a "possible" (or "corona") candidate sitting on a catalogue star
        becomes a star -- the catalogue is independent evidence;
      * a "star" with no catalogue star under it, although its signal says it
        is brighter than the catalogue's limit by a margin, is demoted to
        "unconfirmed" and not lifted -- it is something else;
      * (off, FORCED_ADD) a catalogue star the transform puts on open sky,
        with 4 sigma of signal there but missed by the search, added as a star.
    Without a timestamp, without the catalogue file, or without a clear
    solution nothing changes and the reason is recorded."""
    try:
        from . import starcat as _sc
        rep = {}
        try:
            rep = json.load(open(os.path.join(wd, "report.json")))
        except Exception:
            pass
        ts = rep.get("shot_first")
        from . import observer as _ob
        site = _ob.load(os.path.dirname(os.path.abspath(wd)))
        if not ts and not (site and site.get("ut")):
            return {"status": "no timestamp in report.json and no observing time given"}
        H2, W2 = Zpos.shape

        def zfun(y, x, w=3):
            yb, xb = int(round(y / f)), int(round(x / f))
            if not (0 <= yb < H2 and 0 <= xb < W2):
                return 0.0
            return float(np.nanmax(Zpos[max(0, yb - w):yb + w + 1, max(0, xb - w):xb + w + 1]))

        dets = [c for c in out if c["class"] != "artefact"]
        sol = _sc.match(dets, g, ts, shape, zfun=zfun, progress=progress, site=site)
        # A WRONG POSITION IS WORSE THAN NONE (2026-10-10: a test site entered
        # with E instead of W longitude put the Sun 1.8 deg up instead of ~10, the
        # refraction model then bent the predictions and 5 matches fell to 2).
        # With a place given, also match on the time alone; when that finds
        # clearly more stars, the place does not fit this sky: use the time-only
        # solution and say so.
        site_warning = None
        if site and site.get("ut") and site.get("lat") is not None:
            _nm = lambda s_: 0 if (not s_ or s_.get("ambiguous")) else int(s_.get("matches", 0))
            st_t = {k: v for k, v in site.items() if k not in ("lat", "lon", "lat_text", "lon_text", "height_m")}
            sol_t = _sc.match(dets, g, ts, shape, zfun=zfun, progress=progress, site=st_t)
            if _nm(sol_t) >= _nm(sol) + 2:
                site_warning = (f"the position given fits the stars worse than the time alone "
                                f"({_nm(sol)} vs {_nm(sol_t)} matches); position not used - check it, "
                                f"west longitudes and south latitudes are negative")
                _log(progress, "stars: " + site_warning, None)
                sol = sol_t
        if sol is None:
            return {"status": ("no match (fewer than 2 candidates on catalogue stars)"
                               if site and site.get("ut") and site.get("lat") is not None else
                               "no match (fewer than 3 candidates on catalogue stars)"),
                    "site_warning": site_warning}
        if sol.get("ambiguous"):
            return {"status": "ambiguous", "site_warning": site_warning, **{k: v for k, v in sol.items() if k != "ambiguous"}}
        pairs = sol["pairs"]
        for i, c in enumerate(dets):
            if i in pairs and c["class"] != "artefact":
                c["catalogue"] = pairs[i]
                if c["class"] != "star":
                    c["class"] = "star"
                    c["confirmed"] = "catalogue"
        # magnitude from signal, calibrated on the matched stars: z ~ flux / noise
        zm = [(dets[i]["z"], p["mag"]) for i, p in pairs.items()]
        m0 = float(np.median([m + 2.5 * np.log10(max(z, 1e-3)) for z, m in zm]))
        demoted = 0
        for i, c in enumerate(dets):
            if c["class"] == "star" and i not in pairs:
                est = m0 - 2.5 * np.log10(max(c["z"], 1e-3))
                c["mag_est"] = round(float(est), 2)
                if est <= _sc.MAG_SCORE - 1.0:
                    c["class"] = "unconfirmed"
                    demoted += 1
        added = 0
        have = [(c["y"], c["x"]) for c in dets]
        # The catalogue VERIFIES what the search found; it does not add.
        # Adding a catalogue star from a signal below the search's own list
        # (FORCED_ADD) stays available but off: the catalogue is meant as an
        # independent check, never as a source of stars.
        if not FORCED_ADD:
            sol["predicted"] = []
        R = float(g["R"])
        for p in sol["predicted"]:
            y, x = p["y"], p["x"]
            yb, xb = int(round(y / f)), int(round(x / f))
            if not (0 <= yb < H2 and 0 <= xb < W2) or not sky_ok[yb, xb]:
                continue
            if any(np.hypot(y - hy, x - hx) < 0.05 * R for hy, hx in have):
                continue
            if zfun(y, x) < float(sol.get("forced_z", 4.0)):
                continue
            w = 3
            blk = Zpos[max(0, yb - w):yb + w + 1, max(0, xb - w):xb + w + 1]
            k = np.unravel_index(np.nanargmax(blk), blk.shape)
            Y = (max(0, yb - w) + k[0]) * f + f // 2
            X = (max(0, xb - w) + k[1]) * f + f // 2
            cy, cx = float(g["cy"]), float(g["cx"])
            out.append(dict(y=int(Y), x=int(X), z=round(float(np.nanmax(blk)), 3), sigma=float(f),
                            ell=0.0, peakfrac=0.0, clutter=1.0,
                            r=round(float(np.hypot(Y - cy, X - cx) / R), 3),
                            angle=round(float(np.degrees(np.arctan2(-(Y - cy), X - cx)) % 360), 3),
                            **{"class": "star", "confirmed": "catalogue (forced)",
                               "catalogue": {k2: p[k2] for k2 in ("mag", "ra", "dec")}}))
            added += 1
        info = {"status": "solved", "matches": sol["matches"], "mirror": sol["mirror"],
                "arcsec_per_px": round(sol["arcsec_per_px"], 4),
                "rotation_deg": round(sol["rotation_deg"], 2),
                "north_deg_cw_from_up": round(sol["north_deg_cw_from_up"], 1),
                "east_deg_cw_from_up": round(sol["east_deg_cw_from_up"], 1),
                "fov_deg": [round(v, 3) for v in sol["fov_deg"]],
                "moon_sd_arcmin": round(sol["moon_sd_arcmin"], 2),
                "sun_ra": round(sol["sun_ra"], 4), "sun_dec": round(sol["sun_dec"], 4),
                "date_source": sol.get("date_source", "frames"),
                "site_used": sol.get("site_used", "none"),
                "zenith_deg_cw_from_up": (round(sol["zenith_deg_cw_from_up"], 1)
                                          if "zenith_deg_cw_from_up" in sol else None),
                "sun_alt": round(sol["sun_alt"], 2) if "sun_alt" in sol else None,
                "sun_az": round(sol["sun_az"], 2) if "sun_az" in sol else None,
                "site_warning": site_warning,
                "confirmed": sum(1 for c in out if c.get("confirmed") == "catalogue"),
                "added": added, "demoted": demoted}
        _log(progress, f"stars: catalogue match on {sol['matches']} stars "
                       f"({sol['arcsec_per_px']:.2f}\"/px, north {info['north_deg_cw_from_up']:.1f} deg "
                       f"clockwise from up{', mirrored' if sol['mirror'] else ''}); "
                       f"{info['confirmed']} confirmed, {added} added, {demoted} demoted", None)
        return info
    except Exception as e:
        return {"status": f"failed ({e})"}


def build(wd, progress=None):
    meta = find(wd, progress)
    with open(os.path.join(wd, STARS_FILE), "w") as fh:
        json.dump(meta, fh, indent=1)
    _log(progress, f"stars: {meta['n_star']} found (z >= {meta['threshold']:.1f}), "
                   f"{meta['n_possible']} possible (z {Z_POSSIBLE:.0f}-{meta['threshold']:.1f}; "
                   f"about {meta['expected_by_chance']} expected from noise alone), "
                   f"{meta['n_corona']} set aside inside coronal structure", None)
    return meta


# --------------------------------------------------------------------------
# THE STAR LAYER (step 3). Real data, not drawn points.
#
# Each star found above is lifted out of the corona the detail layers are built
# on (corona_rgb / corona_lum, the merge itself is never touched) and kept as
# its own small patch. Per star:
#   * the background under it is a quadratic surface fitted, per channel and
#     with two rounds of 3-sigma clipping, to a ring of 1-2 star radii around
#     it (other stars masked out);
#   * the star's light is the data minus that surface, smoothed at half the
#     star's width and sharpened back by three Van Cittert steps: the star's
#     shape is taken out, the pixel noise stays behind in the starless picture
#     (no smooth "hole" in the grain);
#   * its colour is the aperture sum of each channel over the luminance,
#     shrunk towards the corona's colour where the star is faint (signal to
#     noise 5 gives half and half), so a faint star is not painted a random
#     colour;
#   * a cosine taper over the last 3 px, so nothing has an edge.
# The render puts the stars back over the finished picture (overlay below),
# like the prominence layer, so no detail filter ever sharpens them into rings.

LAYER_FILE = "stars_layer.npz"
LAYER_BUILD = 1


def _ring_surface(img, ra, mask_other):
    """Quadratic surface per channel fitted on the ring ra+1 .. 2 ra."""
    h, w = img.shape[:2]
    c = (h - 1) / 2.0
    yy, xx = np.mgrid[0:h, 0:w].astype(np.float64)
    yy -= c; xx -= c
    rr = np.hypot(yy, xx)
    ring = (rr >= ra + 1) & (rr <= 2 * ra) & ~mask_other
    A = np.stack([np.ones_like(yy), yy, xx, yy * yy, xx * xx, yy * xx], -1)
    out = np.empty(img.shape, np.float64)
    for ch in range(img.shape[2]):
        v = img[..., ch].astype(np.float64)
        m = ring.copy()
        coef = None
        for _ in range(3):
            if m.sum() < 30:
                break
            coef, *_ = np.linalg.lstsq(A[m], v[m], rcond=None)
            res = v - A @ coef
            s = 1.4826 * np.median(np.abs(res[m])) + 1e-30
            m = ring & (np.abs(res) < 3 * s)
        out[..., ch] = (A @ coef) if coef is not None else np.median(v[ring]) if ring.any() else 0.0
    return out, rr


def remove(wd, progress=None, classes=("star",)):
    """Find the stars, lift them out of corona_rgb / corona_lum in place, and
    write their light to LAYER_FILE. Returns a stats dict."""
    meta = find(wd, progress)
    with open(os.path.join(wd, STARS_FILE), "w") as fh:
        json.dump(meta, fh, indent=1)
    sel = [c for c in meta["stars"] if c["class"] in classes]
    cr = os.path.join(wd, "corona_rgb.npy"); cl = os.path.join(wd, "corona_lum.npy")
    if not (os.path.exists(cr) and os.path.exists(cl)):
        _log(progress, "star layer: no corona files to lift the stars from, skipped", None)
        return {"stars_layer": 0}
    rgb = np.load(cr, mmap_mode="r+")
    lum = np.load(cl, mmap_mode="r+")
    H, W = lum.shape
    boxes, patches, info, centres = [], [], [], []
    taken = [(c["y"], c["x"], float(np.clip(round(3.0 * c["sigma"] + 2), 6, 40))) for c in sel]
    for c in sel:
        ra = float(np.clip(round(3.0 * c["sigma"] + 2), 6, 40))
        hb = int(2 * ra + 3)
        y, x = int(c["y"]), int(c["x"])
        if y - hb < 0 or x - hb < 0 or y + hb >= H or x + hb >= W:
            continue
        img = np.asarray(rgb[y - hb:y + hb + 1, x - hb:x + hb + 1], np.float64)
        n = 2 * hb + 1
        oy, ox = np.mgrid[0:n, 0:n] - hb
        other = np.zeros((n, n), bool)
        for (ty, tx, tr) in taken:
            if (ty, tx) != (c["y"], c["x"]):
                other |= np.hypot(oy + y - ty, ox + x - tx) <= tr + 1
        bg, rr = _ring_surface(img, ra, other)
        ex = img - bg
        E = 0.2126 * ex[..., 0] + 0.7152 * ex[..., 1] + 0.0722 * ex[..., 2]
        # re-centre on the star's own light
        Es0 = ndimage.gaussian_filter(E, max(1.0, c["sigma"] / 2.0))
        k = np.unravel_index(np.argmax(np.where(rr <= ra / 2 + 1, Es0, -np.inf)), Es0.shape)
        rr = np.hypot(oy - (k[0] - hb), ox - (k[1] - hb))
        Es = ndimage.gaussian_filter(E, max(0.8, c["sigma"] / 2.0))
        w = np.where(rr <= ra, 1.0, np.where(rr <= ra + 3, 0.5 * (1 + np.cos(np.pi * (rr - ra) / 3.0)), 0.0))
        ap = rr <= ra
        ring = (rr >= ra + 1) & (rr <= 2 * ra) & ~other
        sE = 1.4826 * np.median(np.abs(E[ring] - np.median(E[ring]))) + 1e-30
        # Smoothing alone flattens a compact core (at sigma 1.1 px it left half
        # of the peak behind) and taking the data pixel for pixel takes the
        # grain with it (a grainless disc 24 px wide at 6.9 R). So the smoothed
        # light is sharpened back by three Van Cittert steps, S += Es - G*S:
        # the star's own shape returns at the scales the smoothing kept, the
        # pixel noise it removed does not.
        ss = max(0.8, c["sigma"] / 2.0)
        S = Es.copy()
        for _ in range(3):
            S += Es - ndimage.gaussian_filter(S, ss)
        # NOT CLIPPED AT ZERO (audit 2026-10-10, A8). 21-37 % of the aperture
        # is negative noise; clipping kept the positive half and lifted
        # 20-65 % too much light (8240 -> 8969 counts on one star), leaving a
        # 32 px dark disc of -0.13 to -0.28 % under every star for the detail
        # layers to sharpen, and a speckle halo in the patch. Instead the
        # sharpened light is tapered by its own significance: pixels at the
        # noise level stay in the corona (both signs), the core is lifted.
        sS = 1.4826 * np.median(np.abs(S[ring] - np.median(S[ring]))) + 1e-30
        from scipy.special import erf as _erf
        gsig = _erf(np.clip(S, 0, None) / (2.0 * sS))
        prof = np.clip(S, 0, None) * gsig * w
        tot = float(E[ap].sum())
        snr = tot / (sE * np.sqrt(ap.sum()))
        bgl = 0.2126 * bg[hb, hb, 0] + 0.7152 * bg[hb, hb, 1] + 0.0722 * bg[hb, hb, 2]
        bgcol = bg[hb, hb] / max(bgl, 1e-30)
        col = np.array([ex[..., i][ap].sum() for i in range(3)]) / max(tot, 1e-30)
        t = snr * snr / (snr * snr + 25.0)
        col = np.clip(t * col + (1 - t) * bgcol, 0.05, None)
        col /= 0.2126 * col[0] + 0.7152 * col[1] + 0.0722 * col[2]
        star = prof[..., None] * col[None, None, :]
        rgb[y - hb:y + hb + 1, x - hb:x + hb + 1] = (img - star).astype(np.float32)
        lum[y - hb:y + hb + 1, x - hb:x + hb + 1] = (
            np.asarray(lum[y - hb:y + hb + 1, x - hb:x + hb + 1], np.float64) - prof).astype(np.float32)
        boxes.append((y - hb, x - hb, n, n))
        centres.append((y - hb + float(k[0]), x - hb + float(k[1]), ra))
        patches.append(prof.astype(np.float32))
        info.append(dict(y=c["y"], x=c["x"], z=c["z"], r=c["r"], radius=ra,
                         flux=round(tot, 3), snr=round(float(snr), 1),
                         colour=[round(float(v), 3) for v in col], peak=float(prof.max())))
    rgb.flush(); lum.flush(); del rgb, lum
    np.savez(os.path.join(wd, LAYER_FILE), build=LAYER_BUILD,
             boxes=np.array(boxes, np.int32).reshape(-1, 4),
             colours=np.array([i["colour"] for i in info], np.float32).reshape(-1, 3),
             peaks=np.array([i["peak"] for i in info], np.float32),
             centres=np.array(centres, np.float32).reshape(-1, 3),
             **{f"p{i}": p for i, p in enumerate(patches)})
    _log(progress, f"star layer: {len(info)} star(s) lifted out of the corona "
                   f"({', '.join('%.1f R snr %.0f' % (i['r'], i['snr']) for i in info) or 'none'}); "
                   f"they come back over the finished picture", None)
    return {"stars_layer": len(info), "stars_layer_info": info,
            "stars_found": meta["n_star"], "stars_possible": meta["n_possible"]}


def load_layer(wd):
    """(boxes, patches, colours, peaks, centres) or None; centres (y, x, radius)."""
    p = os.path.join(wd, LAYER_FILE)
    if not os.path.exists(p):
        return None
    z = np.load(p)
    n = len(z["boxes"])
    if n == 0 or "centres" not in z.files:
        return None
    return (z["boxes"], [z[f"p{i}"] for i in range(n)], z["colours"], z["peaks"], z["centres"])


def rings(ra, q):
    """The ring the display gain is measured in, on a grid decimated by q."""
    return max(ra / q + 2.0, 3.0), max((2.0 * ra + 4.0) / q, 10.0)


def _hp_mad(Y, cy, cx, rin, rout):
    """Robust spread of the 3x3 high-pass of Y in the ring rin..rout about
    (cy, cx): the grain amplitude there."""
    H, W = Y.shape
    r0 = int(np.ceil(rout)) + 1
    y0, y1 = max(1, int(cy) - r0), min(H - 1, int(cy) + r0 + 1)
    x0, x1 = max(1, int(cx) - r0), min(W - 1, int(cx) + r0 + 1)
    if y1 - y0 < 3 or x1 - x0 < 3:
        return 0.0
    w = np.asarray(Y[y0 - 1:y1 + 1, x0 - 1:x1 + 1], np.float64)
    m = sum(w[1 + a:w.shape[0] - 1 + a, 1 + b:w.shape[1] - 1 + b]
            for a in (-1, 0, 1) for b in (-1, 0, 1)) / 9.0
    hp = w[1:-1, 1:-1] - m
    yy, xx = np.mgrid[y0:y1, x0:x1]
    rr = np.hypot(yy - cy, xx - cx)
    v = hp[(rr >= rin) & (rr <= rout)]
    if v.size < 20:
        return 0.0
    return float(1.4826 * np.median(np.abs(v - np.median(v))))


def overlay(rgb, stars_, P, corona_gain=None):
    """Screen the stars over a finished (display-encoded) picture, in place.

    stars_: per star (cy, cx, rin, rout, sig, y0, x0, L, N, C) on this
    picture's grid -- L the star's light in the merge's luminance units, N its
    profile over its peak (for the size control), C its colour at unit
    luminance, sig the pixel noise of the corona there (same units, same grid).

    THE STAR IS SHOWN WITH THE PICTURE'S OWN LOCAL GAIN. The display
    amplification at the star's place is measured on the finished picture
    itself, as the grain amplitude in a ring around the star over the data's
    grain amplitude there, and the star's light is multiplied by it. So at
    strength 1 a star stands out of the grain as far as it stood out of the
    noise in the stack, whatever the detail layers and curves did to that part
    of the picture -- the base envelope alone puts a 4-sigma star at half a
    percent of the display, invisible. Same formula as the page's starOverlay.
    """
    s = float(P.get("starStrength", 1.0))
    if s <= 0 or not stars_:
        return
    sat = float(P.get("starColour", 1.0))
    if int(P.get("colourMode", 0) or 0) == 1:
        sat = 0.0
    Y = 0.2126 * rgb[..., 0] + 0.7152 * rgb[..., 1] + 0.0722 * rgb[..., 2]
    gains = [(_hp_mad(Y, st[0], st[1], st[2], st[3]) / max(st[4], 1e-30)) for st in stars_]
    H, W = rgb.shape[:2]
    for st, g in zip(stars_, gains):
        cy, cx, rin, rout, sig, y0, x0, L, N, C = st
        h, w = L.shape
        ya, xa = max(0, y0), max(0, x0)
        yb, xb = min(H, y0 + h), min(W, x0 + w)
        if yb <= ya or xb <= xa or g <= 0:
            continue
        Ls = L[ya - y0:yb - y0, xa - x0:xb - x0]
        d = Ls * np.float32(g * s)
        col = np.asarray(C, np.float32).reshape(-1)[:3].copy()
        if corona_gain is not None:
            col *= np.asarray(corona_gain, np.float32)
            col /= max(float(0.2126 * col[0] + 0.7152 * col[1] + 0.0722 * col[2]), 1e-6)
        col = np.clip(1.0 + (col - 1.0) * np.float32(sat), 0, None)
        a = np.clip(d[..., None] * col[None, None, :], 0, 1)
        sub = rgb[ya:yb, xa:xb]
        sub[...] = 1.0 - (1.0 - sub) * (1.0 - a)
