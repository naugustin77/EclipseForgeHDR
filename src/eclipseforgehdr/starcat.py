"""STARS AGAINST A CATALOGUE.

The stars the search finds are matched to a star catalogue, so a faint star
that the search could only call "possible" is confirmed when a catalogue star
sits exactly there, and a "star" that matches nothing keeps its class but says
so. Nothing here needs a location or an exact clock:

  * the catalogue (CAT_FILE) is every star to mag 10 within 12 deg of the
    ecliptic -- the Sun is always on it -- J2000 positions with proper motion,
    built by tools/build_star_catalogue.py from the CDS files of the
    Hipparcos and Tycho catalogues (ESA 1997, I/239), with HD numbers from
    Hipparcos/Tycho, HR numbers and Bayer/Flamsteed names from the Bright Star
    Catalogue (V/50) and SAO numbers from the BSC and the SAO catalogue
    (I/131A); `label` is what a planetarium would print: Bayer/Flamsteed,
    else HR, else HD, else SAO, else HIP or TYC;
  * the Sun's position comes from the frames' date and clock, but the clock is
    treated as unknown to +-14 h (cameras keep local time or no time zone at
    all; one test set read 21:30 for 18:30 UT). Over that span the Sun moves along
    the ecliptic only, so the unknown is a shift ALONG the ecliptic, and the
    true Sun must sit under the Moon's centre at totality;
  * the scale comes from the lunar radius in pixels: at any total eclipse the
    Moon's apparent semi-diameter is 14.7-16.9 arcmin;
  * the camera's rotation, and whether the picture is mirrored, are free.

A hypothesis is made from two detections and two catalogue stars (their
separation fixes scale and rotation), kept only if its scale is a Moon's and
it puts the Sun on the ecliptic under the Moon, and scored by how many of the
detections land on catalogue stars. Atmospheric refraction (the Sun was 7 deg
up on one test set; the frame's top and bottom are lifted differently) and lens
distortion are not modelled; the match tolerance grows with the distance from
the centre to cover them: on the reference test set the three matched stars sit
12-37 px from their predicted places at 2.6-6.9 R.

Measured on the reference stack: the two stars and the faint "possible" at
4.3 R match catalogue stars of mag 6.31, 6.83 and 7.66; the coronal knot at
1.46 R and the seven other possibles match nothing.
"""
import os, json
import numpy as np

CAT_FILE = os.path.join(os.path.dirname(__file__), "star_catalogue.npz")
D2R = np.pi / 180.0
SD_RANGE = (14.6 / 60.0, 17.0 / 60.0)    # lunar semi-diameter at a total eclipse, deg
CLOCK_SPAN_H = 14.0
MAG_PAIRS = 8.0                           # catalogue stars that may seed a hypothesis
MAG_SCORE = 8.5                           # ... and that score a hypothesis
# lab 0.28.1 (2026-10-18): the catalogue goes to V 10 (61k stars), but the
# solution is still found and ranked on the stars to MAG_SCORE, so that the
# denser faint catalogue cannot add chance alignments; only once the solution
# is chosen, the remaining detections are matched to the stars to MAG_CONFIRM
# at the fit tolerance. A detection that the search could only call "faint"
# or "possible" is kept when a catalogue star sits under it; the long high
# eclipses that reach past mag 8.5 are what this is for.
MAG_CONFIRM = 10.0
LINE_TOL_R = 0.25                         # Sun to the ecliptic line through the Moon, lunar radii
# Two stars are enough when the place and time are given (lab 0.25.1,
# 2026-10-10): the pair fixes scale, rotation and offset exactly, and is
# accepted only if it puts the Sun within this many lunar radii of the Moon's
# centre, and the scale on a Moon of 14.6-17 arcmin. At totality the Sun sits
# under the Moon (about 0.03 R off at most in 2026), but the fit also carries
# the unmodelled distortion: one test set's 3-star solution puts it 0.062 R off, its
# true pair 0.067 R; the chance pairs in the same search 0.15-0.16 R.
SUN_OFFSET_MAX_R = 0.10
ZP_SPREAD = 0.8                           # mag: photometric consistency of a solution
# Total (and hybrid) solar eclipses, from NASA's eclipse catalogue (F. Espenak).
# Tried when the frames' own date gives no solution -- a camera whose clock
# was never set (one test set's frames say 13 January 2026).
ECLIPSE_DATES = ["2008:08:01", "2009:07:22", "2010:07:11", "2012:11:13", "2013:11:03",
                 "2015:03:20", "2016:03:09", "2017:08:21", "2019:07:02", "2020:12:14",
                 "2021:12:04", "2023:04:20", "2024:04:08", "2026:08:12", "2027:08:02",
                 "2028:07:22", "2030:11:25", "2031:11:14", "2033:03:30", "2034:03:20",
                 "2035:09:02", "2037:07:13", "2038:12:26", "2039:12:15", "2041:04:30",
                 "2042:04:20", "2043:04:09", "2044:08:23", "2045:08:12"]


def _jd(Y, M, D, h):
    if M <= 2:
        Y -= 1; M += 12
    A = Y // 100
    B = 2 - A + A // 4
    return int(365.25 * (Y + 4716)) + int(30.6001 * (M + 1)) + D + B - 1524.5 + h / 24.0


def _precess_matrix(JD):
    T = (JD - 2451545.0) / 36525.0
    z = (2306.2181 * T + 1.09468 * T * T) / 3600 * D2R
    ze = (2306.2181 * T + 0.30188 * T * T) / 3600 * D2R
    th = (2004.3109 * T - 0.42665 * T * T) / 3600 * D2R
    cz, sz, cze, sze, ct, st = np.cos(z), np.sin(z), np.cos(ze), np.sin(ze), np.cos(th), np.sin(th)
    return np.array([[cze * ct * cz - sze * sz, -sze * ct * cz - cze * sz, -st * cz],
                     [cze * ct * sz + sze * cz, -sze * ct * sz + cze * cz, -st * sz],
                     [cze * st, -sze * st, ct]])


def sun_j2000(JD):
    """Apparent Sun (low precision, ~0.01 deg) referred to the J2000 frame."""
    T = (JD - 2451545.0) / 36525.0
    L0 = 280.46646 + 36000.76983 * T
    M = (357.52911 + 35999.05029 * T) * D2R
    C = (1.914602 - 0.004817 * T) * np.sin(M) + 0.019993 * np.sin(2 * M) + 0.000289 * np.sin(3 * M)
    lon = (L0 + C - 0.00569) * D2R
    eps = (23.439291 - 0.0130042 * T) * D2R
    ra = np.arctan2(np.cos(eps) * np.sin(lon), np.cos(lon))
    de = np.arcsin(np.sin(eps) * np.sin(lon))
    v = _precess_matrix(JD).T @ np.array([np.cos(de) * np.cos(ra), np.cos(de) * np.sin(ra), np.sin(de)])
    return float(np.degrees(np.arctan2(v[1], v[0])) % 360), float(np.degrees(np.arcsin(v[2])))


def _gnomonic(ra, de, ra0, de0):
    ra, de, a0, d0 = ra * D2R, de * D2R, ra0 * D2R, de0 * D2R
    cc = np.sin(d0) * np.sin(de) + np.cos(d0) * np.cos(de) * np.cos(ra - a0)
    xi = np.cos(de) * np.sin(ra - a0) / cc
    eta = (np.cos(d0) * np.sin(de) - np.sin(d0) * np.cos(de) * np.cos(ra - a0)) / cc
    return np.degrees(xi), np.degrees(eta), cc


def parse_time(ts):
    """EXIF '2026:08:12 21:30:28' or ISO '2026-08-12T21:30:28' -> (Y, M, D, hours)."""
    try:
        d, t = str(ts).replace("T", " ").split(" ", 1)
        Y, M, D = (int(x) for x in d.replace("-", ":").split(":")[:3])
        hh = [float(x) for x in t.split(":")[:3]]
        return Y, M, D, hh[0] + hh[1] / 60 + (hh[2] if len(hh) > 2 else 0) / 3600
    except Exception:
        return None


def refraction_deg(alt):
    """Standard atmospheric refraction (Saemundsson), degrees, 1010 hPa 10 C."""
    alt = np.asarray(alt, np.float64)
    return 1.02 / np.tan(np.radians(alt + 10.3 / (alt + 5.11))) / 60.0


def altaz_of_date_j2000(ra, dec, JD, lat, lon):
    """Geometric altitude and azimuth (deg, azimuth from north through east) of
    J2000 positions at JD (UT) for a site; precession applied, nutation and
    aberration ignored (well under the match tolerance)."""
    ra = np.atleast_1d(np.asarray(ra, np.float64)); dec = np.atleast_1d(np.asarray(dec, np.float64))
    v = np.stack([np.cos(dec * D2R) * np.cos(ra * D2R), np.cos(dec * D2R) * np.sin(ra * D2R), np.sin(dec * D2R)])
    w = _precess_matrix(JD) @ v
    rd = np.arctan2(w[1], w[0]); dd = np.arcsin(np.clip(w[2], -1, 1))
    T = (JD - 2451545.0) / 36525.0
    gmst = (280.46061837 + 360.98564736629 * (JD - 2451545.0) + 0.000387933 * T * T) % 360.0
    Hh = (gmst + lon) * D2R - rd
    p = lat * D2R
    alt = np.arcsin(np.sin(p) * np.sin(dd) + np.cos(p) * np.cos(dd) * np.cos(Hh))
    az = np.arctan2(-np.cos(dd) * np.sin(Hh), np.sin(dd) * np.cos(p) - np.cos(dd) * np.sin(p) * np.cos(Hh))
    out_alt, out_az = np.degrees(alt), np.degrees(az) % 360.0
    if out_alt.size == 1:
        return float(out_alt[0]), float(out_az[0])
    return out_alt, out_az


def load_catalogue():
    if not os.path.exists(CAT_FILE):
        return None
    z = np.load(CAT_FILE)
    return {k: z[k] for k in z.files}


def _match(dets, geo, timestamp, shape, zfun=None, progress=None, site=None):
    """dets: list of dicts with y, x, z (full-res merge pixels). Returns a dict
    with the solution and, per detection index, the catalogue star it sits on;
    or None when no solution stands out."""
    cat = load_catalogue()
    have_time = bool(site and site.get("ut"))
    have_place = have_time and site.get("lat") is not None and site.get("lon") is not None
    tm = tuple(site["ut"]) if have_time else parse_time(timestamp)
    if cat is None or tm is None or len(dets) < 2:
        return None
    cy, cx, R = float(geo["cy"]), float(geo["cx"]), float(geo["R"])
    JD = _jd(tm[0], tm[1], tm[2], tm[3])
    ra0, de0 = sun_j2000(JD)
    yrs = (JD - 2451545.0) / 365.25
    ra = cat["ra"] + cat["pmra"] / 3.6e6 * yrs / np.cos(np.radians(cat["dec"]))
    de = cat["dec"] + cat["pmdec"] / 3.6e6 * yrs
    xi, eta, cc = _gnomonic(ra, de, ra0, de0)
    H, W = shape
    if have_place:
        from .observer import refraction_factor
        rf = refraction_factor(site.get("height_m"))
        a0, z0 = altaz_of_date_j2000(ra0, de0, JD, site["lat"], site["lon"])
        a0 = a0 + rf * float(refraction_deg(a0))

    def project(r_, d_):
        """Sky plane about the Sun, as seen (no mirror): complex x + i y in
        degrees. With the site known: the horizon frame, refraction applied,
        x towards increasing azimuth, y up. Without: the equatorial frame,
        x = -xi (east to the left), y = eta (north up)."""
        if have_place:
            al, az = altaz_of_date_j2000(r_, d_, JD, site["lat"], site["lon"])
            al, az = np.atleast_1d(al).astype(np.float64), np.atleast_1d(az).astype(np.float64)
            al = al + rf * refraction_deg(np.maximum(al, -1.0))
            A, Z, A0, Z0 = al * D2R, az * D2R, a0 * D2R, z0 * D2R
            c2 = np.sin(A0) * np.sin(A) + np.cos(A0) * np.cos(A) * np.cos(Z - Z0)
            x = np.cos(A) * np.sin(Z - Z0) / c2
            y = (np.cos(A0) * np.sin(A) - np.sin(A0) * np.cos(A) * np.cos(Z - Z0)) / c2
            return np.degrees(x) + 1j * np.degrees(y)
        x_, e_, _ = _gnomonic(np.atleast_1d(r_), np.atleast_1d(d_), ra0, de0)
        return -x_ + 1j * e_
    reach = np.hypot(H, W) / R * SD_RANGE[1] + 1.0          # deg: frame half-diagonal + clock span
    near = (cc > 0) & (np.hypot(xi, eta) < reach)
    xi, eta, mag = xi[near], eta[near], cat["mag"][near]
    names = cat["name"][near] if "name" in cat else np.array([""] * int(near.sum()))
    nn = int(near.sum())
    ids = {k: (cat[k][near] if k in cat else (np.zeros(nn, np.int32) if k != "label" else np.array([""] * nn)))
           for k in ("label", "hd", "hr", "sao", "hip")}

    def _ids(c):
        return dict(label=str(ids["label"][c]), hd=int(ids["hd"][c]), hr=int(ids["hr"][c]),
                    sao=int(ids["sao"][c]), hip=int(ids["hip"][c]))
    cra, cde = ra[near], de[near]
    # the ecliptic's direction at the Sun, from the Sun a day later
    ra1, de1 = sun_j2000(JD + 1.0)
    e1 = project(np.array([ra1]), np.array([de1]))[0] - project(np.array([ra0]), np.array([de0]))[0]
    edir = e1 / abs(e1)
    # deg the Sun can be off along the ecliptic: the whole clock span, or
    # with the time given only the Moon-to-Sun offset and a few minutes
    span = (CLOCK_SPAN_H / 24.0 * 0.99) if not have_time else 0.15 * SD_RANGE[1]
    sky = project(cra, cde)
    D = np.array([complex(d["x"] - cx, -(d["y"] - cy)) for d in dets])
    # "faint" candidates (below the search's own floor, stars.Z_VERIFY) may be
    # confirmed by a solution but never seed one
    order = np.array([i for i in np.argsort([-d["z"] for d in dets])
                      if not dets[i].get("faint")][:8], dtype=int)
    seed = np.flatnonzero(mag <= MAG_PAIRS)
    scale_lo, scale_hi = R / SD_RANGE[1], R / SD_RANGE[0]  # px per deg
    score_cat = np.flatnonzero(mag <= MAG_SCORE)
    rD = np.abs(D)
    # Refraction and distortion are not modelled, so the tolerance grows with
    # the distance from the centre: two stars fix a hypothesis exactly and a
    # third, 2.6-6.9 R out, then sits up to ~80 px off on one test frame (the
    # 7-deg-high Sun lifts the frame's top and bottom by different amounts).
    # Loose for finding, tighter once the fit has used every match.
    if have_place:
        # refraction modelled: on that frame 12-16 px at 2.6-6.9 R
        tol_find = 12.0 + 0.015 * rD
        tol_fit = 8.0 + 0.008 * rD
    else:
        tol_find = 12.0 + 0.035 * rD
        tol_fit = 10.0 + 0.02 * rD
    min_match = 2 if have_place else 3

    def evaluate(a, b, mir, tol):
        S = np.conj(sky[score_cat]) if mir else sky[score_cat]
        P = a * S + b
        dist = np.abs(D[:, None] - P[None, :])
        j = np.argmin(dist, axis=1)
        dmin = dist[np.arange(len(D)), j]
        used, pairs = set(), {}
        for i in np.argsort(dmin):
            if dmin[i] < tol[i] and j[i] not in used:
                used.add(j[i]); pairs[int(i)] = int(score_cat[j[i]])
        return pairs

    def missed(a, b, mir, pairs):
        """Catalogue stars at least as bright as the faintest matched one that
        the transform puts on open sky in the frame, with no detection there
        and (when the caller can look) no signal either. A real solution has
        none: if a mag-7.7 star was seen, a mag-6.5 star in the same sky would
        have been too. This is what tells a chance alignment from the sky."""
        if not pairs:
            return 99
        lim = max(mag[c] for c in pairs.values()) - 0.3
        S = np.conj(sky) if mir else sky
        k = np.flatnonzero(mag <= lim)
        P = a * S[k] + b
        py, px = cy - P.imag, cx + P.real
        m = 0.06 * min(H, W)
        ins = (px > m) & (px < W - m) & (py > m) & (py < H - m) & (np.abs(P) > 1.6 * R)
        n = 0
        mc = set(pairs.values())
        for kk, y_, x_ in zip(k[ins], py[ins], px[ins]):
            if kk in mc:
                continue
            if zfun is not None and zfun(y_, x_) >= 3.0:
                continue
            n += 1
        return n

    hyps = []
    for mir in (False, True):
        Sk = np.conj(sky) if mir else sky
        for ii in range(len(order)):
            for jj in range(len(order)):
                if ii == jj:
                    continue
                i, j = order[ii], order[jj]
                dD = D[j] - D[i]
                if abs(dD) < 0.5 * R:
                    continue
                # all seed pairs (p, q) at once
                Cp = Sk[seed][:, None]; Cq = Sk[seed][None, :]
                dC = Cq - Cp
                with np.errstate(divide="ignore", invalid="ignore"):
                    a = dD / dC
                s = np.abs(a)
                good = (s > scale_lo) & (s < scale_hi)
                if not good.any():
                    continue
                pp, qq = np.nonzero(good)
                a = a[pp, qq]
                b = D[i] - a * Sk[seed][pp]
                # the Moon's centre (image 0) back on the sky: must lie on the
                # ecliptic through the nominal Sun, within the clock span
                zm = -b / a
                ed = np.conj(edir) if mir else edir
                along = (zm * np.conj(ed)).real
                across = (zm * np.conj(ed)).imag
                sd = R / np.abs(a)                            # implied lunar radius, deg
                keep = (np.abs(across) < LINE_TOL_R * sd) & (np.abs(along) < span)
                for a_, b_ in zip(a[keep], b[keep]):
                    hyps.append((a_, b_, mir))
    if not hyps:
        return None
    def refine(pr, mir):
        di = np.array([D[i] for i in pr]); ci = np.array([sky[c] for c in pr.values()])
        if mir:
            ci = np.conj(ci)
        if len(pr) < 3:
            if have_place:
                # two stars with the site known: the exact similarity through
                # both; the Sun-offset and scale checks below decide
                a_ = complex((di[1] - di[0]) / (ci[1] - ci[0]))
                return a_, complex(di[0] - a_ * ci[0])
            # two stars: the Sun is held under the Moon's centre (offset 0)
            # so the pair still has something to disagree with
            a_ = complex(np.vdot(ci, di) / np.vdot(ci, ci))
            return a_, 0j
        A = np.stack([ci, np.ones_like(ci)], 1)
        sol = np.linalg.lstsq(A, di, rcond=None)[0]
        return sol[0], sol[1]

    scored, seen = [], []
    for a, b, mir in hyps:
        pr = evaluate(a, b, mir, tol_find)
        if len(pr) < min_match:
            continue
        a, b = refine(pr, mir)
        pr = evaluate(a, b, mir, tol_fit)
        if len(pr) < min_match:
            continue
        a, b = refine(pr, mir)
        pr = evaluate(a, b, mir, tol_fit)
        if len(pr) < min_match:
            continue
        if len(pr) == 2 and have_place:
            if abs(b) > SUN_OFFSET_MAX_R * R or not (scale_lo < abs(a) < scale_hi):
                continue
        key = (mir, round(float(np.degrees(np.angle(a))) / 2.0), tuple(sorted(pr)))
        if key in seen:
            continue
        seen.append(key)
        # PHOTOMETRIC CONSISTENCY: the search's z is the star's signal over
        # the local noise, so mag + 2.5 log10(z) is roughly one number for
        # every matched star. On the true solutions it spreads 0.1-0.55 mag
        # (two test sets); the chance three-star alignments found when
        # the clock is given a wrong eclipse date spread 0.95-1.3.
        zp = sorted(mag[c] + 2.5 * np.log10(max(dets[i]["z"], 1e-3)) for i, c in pr.items())
        spread = (zp[-1] - zp[0]) if len(zp) <= 3 else (zp[-2] - zp[1])
        if spread > ZP_SPREAD:
            continue
        # a mirrored picture (a star diagonal) is the exception: it needs four
        if mir and len(pr) < 4:
            continue
        miss = missed(a, b, mir, pr)
        sc = len(pr) - 1.5 * miss - (1.0 if mir else 0.0)
        if len(pr) == 2:
            # TWO STARS (audit 2026-10-10): every pair scored exactly 2.0, so
            # a chance pair inside the Sun-offset window tied with the true one
            # and the match was dropped as ambiguous (4 of 40 trials). Rank the
            # pairs by how far they put the Sun from the Moon's centre, by the
            # magnitude consistency, and by whether the two candidates are
            # stars in their own right (class "star") rather than "possible".
            _weak = sum(1 for i in pr if dets[i].get("class") not in ("star", None))
            sc -= 4.0 * abs(b) / R + spread / ZP_SPREAD + 0.5 * _weak
        scored.append((sc, len(pr), miss, a, b, mir, pr))
    if not scored:
        return None
    scored.sort(key=lambda t: -t[0])
    best = scored[0]
    sc, nmatch, miss, a, b, mir, pr = best
    rivals = [t for t in scored[1:] if abs(np.angle(t[3] / a)) > np.radians(2.0) or t[5] != mir]
    _amb = 0.5 if nmatch > 2 else 0.12          # two-star scores are continuous
    if miss > 0 or (rivals and rivals[0][0] >= sc - _amb):
        return {"ambiguous": True, "matches": nmatch, "missed_bright": miss,
                "rival_matches": rivals[0][1] if rivals else 0}
    rot = np.angle(a)
    best = (None, None, None, None, mir)
    S = np.conj(sky) if best[4] else sky
    # CONFIRMATION PASS: the solution stands; the detections it did not use
    # are now tried against the catalogue to MAG_CONFIRM, nearest star within
    # the fit tolerance, each catalogue star at most once
    faint_cat = np.flatnonzero((mag > MAG_SCORE) & (mag <= MAG_CONFIRM))
    confirm_lim, confirm_miss = MAG_CONFIRM, None
    if len(faint_cat):
        P = a * S[faint_cat] + b
        used = set(pr.values())
        free = [i for i in range(len(D)) if i not in pr]
        # HOW DEEP THE PICTURE GOES: the same logic as missed() above, now for
        # the faint stars. Walk the catalogue stars the solution puts on open
        # sky (inside the frame's margins, clear of the corona) from bright
        # to faint; the first one with no signal at its place sets the limit.
        # A candidate on a star fainter than that is not confirmed: if a
        # mag-9.2 star left nothing in the stack, a detection on a mag-9.8
        # star is noise that happens to sit there (2026-10-18).
        py, px = cy - P.imag, cx + P.real
        m = 0.06 * min(H, W)
        ins = (px > m) & (px < W - m) & (py > m) & (py < H - m) & (np.abs(P) > 1.6 * R)
        for q in np.argsort(mag[faint_cat]):
            if not ins[q]:
                continue
            seen = zfun is not None and zfun(py[q], px[q]) >= 3.0
            if not seen and len(D):
                seen = bool(np.min(np.abs(D - P[q])) < (tol_fit.max() if np.ndim(tol_fit) else tol_fit))
            if not seen:
                confirm_miss = float(mag[faint_cat[q]])
                confirm_lim = min(MAG_CONFIRM, confirm_miss + 0.3)
                break
        if free:
            dist = np.abs(D[free][:, None] - P[None, :])
            j = np.argmin(dist, axis=1)
            dmin = dist[np.arange(len(free)), j]
            for q in np.argsort(dmin):
                i = free[q]
                c = int(faint_cat[j[q]])
                if dmin[q] < tol_fit[i] and c not in used and mag[c] <= confirm_lim:
                    used.add(c); pr[int(i)] = c
    resid = {i: float(abs(D[i] - (a * S[c] + b))) for i, c in pr.items()}
    out = {"ambiguous": False, "matches": nmatch, "confirmed_faint": len(pr) - nmatch,
           "confirm_limit_mag": round(float(confirm_lim), 2),
           "confirm_missing_mag": (round(confirm_miss, 2) if confirm_miss is not None else None),
           "mirror": bool(best[4]), "score": float(sc),
           "scale_px_per_deg": float(abs(a)), "arcsec_per_px": float(3600.0 / abs(a)),
           "rotation_deg": float(np.degrees(rot)),
           "moon_sd_arcmin": float(R / abs(a) * 60.0),
           "sun_ra": ra0, "sun_dec": de0,
           "pairs": {int(i): dict(mag=float(mag[c]), ra=float(cra[c]), dec=float(cde[c]),
                                   name=str(names[c]), resid_px=round(resid[i], 1), **_ids(c))
                     for i, c in pr.items()}}
    # where north and east (and, with the site known, the zenith) point in
    # the picture, clockwise from "up", from small steps on the sky
    def _dir(v):
        w = a * (np.conj(v) if best[4] else v)
        return float(np.degrees(np.arctan2(w.real, w.imag)) % 360)
    p0 = project(np.array([ra0]), np.array([de0]))[0]
    pn = project(np.array([ra0]), np.array([de0 + 0.2]))[0] - p0
    pe = project(np.array([ra0 + 0.2 / np.cos(de0 * D2R)]), np.array([de0]))[0] - p0
    out["north_deg_cw_from_up"] = _dir(pn)
    out["east_deg_cw_from_up"] = _dir(pe)
    out["site_used"] = "time and place" if have_place else ("time" if have_time else "none")
    if have_place:
        out["zenith_deg_cw_from_up"] = _dir(1j)
        out["sun_alt"], out["sun_az"] = float(a0), float(z0)
    out["forced_z"] = 3.5 if have_place else 4.0
    out["fov_deg"] = [float(W / abs(a)), float(H / abs(a))]
    # catalogue stars bright enough to be expected, predicted in the frame:
    # positions on the merge, for a forced look by the caller
    P = a * S + b
    py, px = cy - P.imag, cx + P.real
    inside = (px > 0) & (px < W) & (py > 0) & (py < H)
    out["predicted"] = [dict(y=float(py[k]), x=float(px[k]), mag=float(mag[k]),
                             ra=float(cra[k]), dec=float(cde[k]), name=str(names[k]), **_ids(k))
                        for k in np.flatnonzero(inside & (mag <= MAG_SCORE))]
    return out


def match(dets, geo, timestamp, shape, zfun=None, progress=None, site=None):
    """dets: list of dicts with y, x, z (full-res merge pixels). Returns the
    solution dict, {"ambiguous": True, ...}, or None.

    First with the frames' own date and clock. If that gives no solution, each
    total eclipse in ECLIPSE_DATES is tried at noon UT (+-14 h covers the day);
    a date is taken only with at least five matches and two more than any
    other date -- with three, a wrong date gave a unique-looking chance fit in
    2 of 20 trials on one test stack."""
    if site and site.get("ut"):
        r = _match(dets, geo, timestamp, shape, zfun, progress, site=site)
        if r is not None and not r.get("ambiguous"):
            r["date_source"] = "observing site and time"
        return r
    r = _match(dets, geo, timestamp, shape, zfun, progress)
    if r is not None and not r.get("ambiguous"):
        r["date_source"] = "frames"
        return r
    found = []
    for d in ECLIPSE_DATES:
        if timestamp and str(timestamp).replace("-", ":")[:10] == d:
            continue
        rr = _match(dets, geo, d + " 12:00:00", shape, zfun, progress)
        if rr is not None and not rr.get("ambiguous"):
            found.append((rr["matches"], d, rr))
    if not found:
        return r
    found.sort(key=lambda t: -t[0])
    top = found[0][2]

    def same_sky(x):
        # the same day of the year puts the Sun in the same place: 2026-08-12
        # and 2045-08-12 give one and the same solution
        return (abs(((x["sun_ra"] - top["sun_ra"] + 180) % 360) - 180) < 1.0
                and abs(x["sun_dec"] - top["sun_dec"]) < 1.0)
    group = [t for t in found if same_sky(t[2])]
    others = [t for t in found if not same_sky(t[2])]
    if found[0][0] >= 5 and (not others or others[0][0] <= found[0][0] - 2):
        try:
            yr = int(str(timestamp)[:4])
        except Exception:
            yr = 2026
        n, d, rr = min(group, key=lambda t: (-t[0], abs(int(t[1][:4]) - yr)))
        rr["date_source"] = ("eclipse table (" + d + "; the frames say " + str(timestamp) + ")")
        return rr
    return r
