"""Polar map: the finished picture unwrapped around the Sun (0.30.0).

Position angle across, distance from the centre up -- the corona's streamers
stand as vertical curtains, and their position angles can be read off.

ANGLE. Position angle, counter-clockwise on the sky through east, from solar
north by default. Celestial north comes from the star match (stars.json:
north and east as degrees clockwise from up in the stack frame, which also
gives the mirror sense); solar north is turned from it by the P angle of the
date (Meeus, Astronomical Algorithms, ch. 29). Without a star match the angle
starts at image up, counter-clockwise, and the axis label says so.

RADIUS. In solar radii when the scale is known (star match: arcsec per pixel;
solar semi-diameter from the date), else in Moon radii. The map starts at the
Moon's limb.

SAMPLING. One row per pixel of radius at the picture's scale, 0.1 deg per
column, six sub-samples per column averaged so the outer rows, where a column
spans many pixels, do not alias. Bilinear. Where the frame does not reach,
neutral grey.
"""
import json
import os

import numpy as np
from scipy import ndimage

NO_DATA = np.array([0.16, 0.16, 0.17], np.float32)
NT_FULL = 3600          # columns: 0.1 deg
SUBS = 6                # sub-samples per column


def _jd(Y, M, D, h):
    if M <= 2:
        Y -= 1
        M += 12
    A = Y // 100
    B = 2 - A + A // 4
    return int(365.25 * (Y + 4716)) + int(30.6001 * (M + 1)) + D + B - 1524.5 + h / 24.0


def _sun(JD):
    """Apparent ecliptic longitude (deg), true obliquity (deg), distance (AU)."""
    T = (JD - 2451545.0) / 36525.0
    L0 = 280.46646 + 36000.76983 * T + 0.0003032 * T * T
    M = np.radians(357.52911 + 35999.05029 * T - 0.0001537 * T * T)
    e = 0.016708634 - 0.000042037 * T
    C = ((1.914602 - 0.004817 * T) * np.sin(M) + (0.019993 - 0.000101 * T) * np.sin(2 * M)
         + 0.000289 * np.sin(3 * M))
    v = M + np.radians(C)
    R = 1.000001018 * (1 - e * e) / (1 + e * np.cos(v))
    Om = np.radians(125.04 - 1934.136 * T)
    lam = L0 + C - 0.00569 - 0.00478 * np.sin(Om)
    eps = 23 + (26 + (21.448 - 46.815 * T) / 60) / 60 + 0.00256 * np.cos(Om)
    return float(lam), float(eps), float(R)


def p_angle(JD):
    """Position angle of the Sun's rotation axis from celestial north, degrees,
    positive towards east (Meeus ch. 29; his example 1992 Oct 13.0 gives 26.27)."""
    lam, eps, _ = _sun(JD)
    K = 73.6667 + 1.3958333 * (JD - 2396758.0) / 36525.0
    x = np.degrees(np.arctan(-np.cos(np.radians(lam)) * np.tan(np.radians(eps))))
    y = np.degrees(np.arctan(-np.cos(np.radians(lam - K)) * np.tan(np.radians(7.25))))
    return float(x + y)


def solar_sd_arcsec(JD):
    return 959.63 / _sun(JD)[2]


def _time_jd(wd, folder):
    """JD of the shoot: the observing time entered, else the first frame's
    timestamp (a camera clock off by hours moves P by hundredths of a degree)."""
    try:
        from . import observer as _ob
        site = _ob.load(folder)
        if site and site.get("ut"):
            return _jd(*tuple(site["ut"])[:4])
    except Exception:
        pass
    try:
        from .starcat import parse_time
        rep = json.load(open(os.path.join(wd, "report.json")))
        tm = parse_time(rep.get("shot_first"))
        if tm:
            return _jd(*tm)
    except Exception:
        pass
    return None


def orientation(wd, folder, north="solar", manual_cw=None, offset_deg=0.0):
    """Where the angle starts and which way it runs, in the stack frame.

    `manual_cw`: solar north entered by hand, degrees clockwise from image up;
    used instead of the star match (or without one). `offset_deg`: a small
    turn added to whatever angle is used, for fine adjustment.

    Returns dict: zero_cw (deg clockwise from up where theta = 0), sgn (+1 when
    theta runs counter-clockwise in the picture), px_per_unit, unit ("sun" |
    "moon"), origin ("solar" | "celestial" | "up"), p (P angle used), note."""
    o = dict(zero_cw=0.0, sgn=1, px_per_unit=None, unit="moon", origin="up", p=0.0, note="")
    try:
        cat = (json.load(open(os.path.join(wd, "stars.json"))).get("catalogue") or {})
    except Exception:
        cat = {}
    JD = _time_jd(wd, folder)
    if cat.get("status") == "solved":
        n_cw = float(cat["north_deg_cw_from_up"])
        e_cw = float(cat["east_deg_cw_from_up"])
        o["sgn"] = 1 if ((e_cw - n_cw) % 360.0) > 180.0 else -1
        o["zero_cw"] = n_cw
        o["origin"] = "celestial"
        if north == "solar" and JD is not None:
            o["p"] = p_angle(JD)
            o["zero_cw"] = n_cw - o["sgn"] * o["p"]
            o["origin"] = "solar"
        if JD is not None and cat.get("arcsec_per_px"):
            o["px_per_unit"] = solar_sd_arcsec(JD) / float(cat["arcsec_per_px"])
            o["unit"] = "sun"
        o["note"] = "north from %d catalogue stars" % int(cat.get("matches", 0))
    else:
        o["note"] = "no star match: angle from image up"
    if manual_cw is not None:
        # by hand: solar north where the user says. The sense (east
        # counter-clockwise from north) is the sky's as a camera records it;
        # a star match, when there is one, still decides it.
        # North by hand counts as much as a star match: celestial north
        # follows from it through the P angle of the date, the same relation
        # the star match uses the other way round.
        sol = float(manual_cw) % 360.0
        o["zero_cw"] = sol
        o["origin"] = "solar"
        o["p"] = 0.0
        o["note"] = "solar north entered by hand (%.1f deg from image up)" % sol
        if north == "celestial" and JD is not None:
            o["p"] = p_angle(JD)
            o["zero_cw"] = (sol + o["sgn"] * o["p"]) % 360.0
            o["origin"] = "celestial"
            o["note"] = ("celestial north from solar north entered by hand (%.1f deg from image up) "
                         "and the P angle %+.1f deg" % (sol, o["p"]))
    if offset_deg:
        o["zero_cw"] = (o["zero_cw"] + float(offset_deg)) % 360.0
        o["note"] += ", turned %+.1f deg" % float(offset_deg)
    return o


def unwrap(rgb, cy, cx, R_moon, o, outer="frame", step=1.0, radial="linear", centre="S",
           mirrored=False):
    """The picture (H, W, C float 0-1, stack frame) as a polar map.

    Returns (map (rows, cols, C) float32, meta). Rows run outward from the
    top: row 0 is the outermost radius, as the map is read with r up."""
    H, W = rgb.shape[:2]
    unit = o["px_per_unit"] or float(R_moon)
    r0 = float(R_moon) / unit
    rfull = min(cy, H - 1 - cy, cx, W - 1 - cx) / unit
    rcorner = max(np.hypot(y, x) for y in (cy, H - 1 - cy) for x in (cx, W - 1 - cx)) / unit
    r1 = rfull if outer == "circle" else rcorner
    nt = int(NT_FULL / step)
    if radial == "log":
        # CONFORMAL: equal steps in ln r, one step = one column in radians, so
        # a small round thing (a star) stays round at every radius
        rs = np.exp(np.arange(np.log(r0), np.log(r1), 2.0 * np.pi / nt))[::-1]
    else:
        rs = np.arange(r0, r1, step / unit)[::-1]
    # NOT MIRRORED: with r up, the angle has to run clockwise on the sky from
    # left to right, so the map reads like the sky at the limb (2026-11-08: a
    # large prominence looked mirrored with the angle rising to the right). Column 0 is 360 deg, the last column 0 deg.
    # MIDDLE: N, E, S or W sits at the centre column. The angle normally
    # FALLS to the right (sky-true, see above); `mirrored` lets it rise to the
    # right instead -- the classic layout, a mirror image of the sky.
    dirn = 1.0 if mirrored else -1.0
    mid = {"N": 0.0, "E": 90.0, "S": 180.0, "W": 270.0}.get(str(centre).upper(), 180.0)
    start = mid - dirn * 180.0                       # angle at the left edge
    if start <= 0 and dirn < 0:
        start += 360.0
    th = (start + dirn * (np.arange(nt * SUBS) + 0.5) / (nt * SUBS) * 360.0) % 360.0
    phi = np.radians(o["zero_cw"] - o["sgn"] * th)          # clockwise from up
    sp, cp = np.sin(phi), np.cos(phi)
    C = rgb.shape[2] if rgb.ndim == 3 else 1
    chans = [np.ascontiguousarray(rgb[..., c] if rgb.ndim == 3 else rgb, np.float32) for c in range(C)]
    out = np.empty((len(rs), nt, C), np.float32)
    cov = np.empty((len(rs), nt), np.float32)
    for a in range(0, len(rs), 96):
        rr = rs[a:a + 96, None] * unit
        yy = cy - rr * cp
        xx = cx + rr * sp
        ok = (yy >= 0) & (yy <= H - 1) & (xx >= 0) & (xx <= W - 1)
        n = ok.reshape(len(rr), nt, SUBS).sum(2)
        for c in range(C):
            v = ndimage.map_coordinates(chans[c], [yy, xx], order=1, mode="nearest")
            v = np.where(ok, v, 0.0).reshape(len(rr), nt, SUBS).sum(2)
            out[a:a + 96, :, c] = v / np.maximum(n, 1)
        cov[a:a + 96] = n / float(SUBS)
    nd = NO_DATA[:C] if C == 3 else NO_DATA[:1].mean(keepdims=True)
    out = np.where((cov > 0.5)[..., None], out, nd)
    return np.clip(out, 0.0, 1.0), dict(r0=r0, r1=float(rs[0]) if len(rs) else r1,
                                        rfull=rfull, unit=o["unit"], origin=o["origin"],
                                        radial=radial, start=start, dirn=dirn)


# ---------------------------------------------------------------- annotated JPG

def _font(size):
    from PIL import ImageFont
    try:
        return ImageFont.load_default(size=size)
    except TypeError:                     # Pillow < 10.1: fixed bitmap font
        return ImageFont.load_default()


EARTH_DIAM_RSUN = 12742.0 / 695700.0   # Earth's mean diameter in solar radii


def _canvas(h0, w0, width, s):
    """Canvas geometry for a map of h0 x w0 drawn `width` wide, margins scaled by s."""
    pw = int(width)
    ph = max(int(240 * s), int(round(pw * h0 / w0)))   # the map's own shape: no stretch on top
    L, T, Rm, B = (int(round(v * s)) for v in (110, 56, 36, 100))
    return pw, ph, L, T, Rm, B


def _overlay(meta, size, pw, ph, L, T, s, grid=True, earth=True):
    """Everything drawn over the map -- guide lines, frame, ticks, numbers,
    N/E/S/W, axis titles, the Earth -- on a transparent canvas of `size`."""
    from PIL import Image, ImageDraw
    ov = Image.new("RGBA", size, (0, 0, 0, 0))
    d = ImageDraw.Draw(ov)
    w1, w2 = max(1, int(round(s))), max(2, int(round(2 * s)))
    r0, r1 = meta["r0"], meta["r1"]
    lg = meta.get("radial") == "log"
    if lg:
        yof = lambda r: T + (np.log(r1) - np.log(r)) / (np.log(r1) - np.log(r0)) * ph
    else:
        yof = lambda r: T + (r1 - r) / (r1 - r0) * ph
    S0 = meta.get("start", 360.0)                     # angle at the left edge; falls to the right
    xk = lambda k: L + k / 12.0 * pw                   # 12 steps of 30 deg
    D0 = meta.get("dirn", -1.0)

    def ang(k):
        raw = S0 + D0 * 30 * k
        v = int(round(raw % 360))
        return 360 if (v == 0 and raw > 0 and k in (0, 12)) else v
    span = r1 - r0
    rstep = 0.5 if span <= 8 else (1.0 if span <= 16 else 2.0)
    if lg:
        rt = [v for v in (1.0, 1.25, 1.5, 2, 2.5, 3, 4, 5, 6, 8, 10, 12, 15, 20) if r0 <= v <= r1]
    else:
        rt = list(np.arange(np.ceil(r0 / rstep) * rstep, r1 + 1e-9, rstep))
    if grid:
        for k in range(1, 12):
            t = ang(k)
            d.line([(xk(k), T), (xk(k), T + ph)], fill=(255, 255, 255, 56 if t % 90 == 0 else 34),
                   width=w2 if t % 90 == 0 else w1)
        for r in rt:
            whole = abs(r - round(r)) < 1e-6
            d.line([(L, yof(r)), (L + pw, yof(r))], fill=(255, 255, 255, 52 if whole else 30),
                   width=w2 if whole else w1)
    if meta["r1"] > meta["rfull"] + 0.05:
        y = yof(meta["rfull"])
        x = L
        while x < L + pw:
            d.line([(x, y), (min(x + 12 * s, L + pw), y)], fill=(255, 255, 255, 120), width=w1)
            x += 24 * s
    f, fs = _font(int(round(22 * s))), _font(int(round(18 * s)))
    tc, lc = (200, 204, 214, 255), (150, 156, 170, 255)
    d.rectangle([L, T, L + pw, T + ph], outline=(85, 85, 102, 255), width=w1)
    for k in range(13):
        d.line([(xk(k), T + ph), (xk(k), T + ph + 8 * s)], fill=lc, width=w1)
        d.text((xk(k), T + ph + 12 * s), str(ang(k)), fill=lc, font=fs, anchor="mt")
    for r in rt:
        d.line([(L - 8 * s, yof(r)), (L, yof(r))], fill=lc, width=w1)
        d.text((L - 12 * s, yof(r)), ("%g" % r), fill=lc, font=fs, anchor="rm")
    if meta["origin"] in ("solar", "celestial"):
        for k in range(13):
            s_ = {0: "N", 90: "E", 180: "S", 270: "W", 360: "N"}.get(ang(k))
            if s_:
                d.text((xk(k), T - 10 * s), s_, fill=(235, 238, 245, 255), font=f, anchor="mb")
    if meta["r1"] > meta["rfull"] + 0.05:
        d.text((L + pw - 10 * s, yof(meta["rfull"]) - 6 * s), "full circle in the frame up to here",
               fill=lc, font=fs, anchor="rb")
    if earth and meta["unit"] == "sun":
        ex, ey = L + 34 * s, T + 34 * s
        if lg:      # conformal: to scale in both directions, for the radius it sits at
            rh = np.exp(np.log(r1) - (ey - T) / ph * (np.log(r1) - np.log(r0)))
            rad = max(1.5, 0.5 * EARTH_DIAM_RSUN / rh * ph / (np.log(r1) - np.log(r0)))
            lab = "Earth, to scale at this distance"
        else:
            rad = max(1.5, 0.5 * EARTH_DIAM_RSUN * ph / (r1 - r0))
            lab = "Earth, to the radial scale"
        d.ellipse([ex - rad, ey - rad, ex + rad, ey + rad], fill=(255, 255, 255, 255))
        d.text((ex + rad + 10 * s, ey), lab, fill=lc, font=fs, anchor="lm")
    xl = {"solar": "position angle (deg) from solar north through east",
          "celestial": "position angle (deg) from celestial north through east",
          "up": "angle (deg) from image up, counter-clockwise (no star match)"}[meta["origin"]]
    if D0 > 0:
        xl += "  (mirrored view)"
    d.text((L + pw / 2, T + ph + 48 * s), xl, fill=tc, font=f, anchor="mt")
    yl = ("r / solar radius" if meta["unit"] == "sun" else "r / Moon radius") + (" (log)" if lg else "")
    tw, th = int(600 * s), int(40 * s)
    tmp = Image.new("RGBA", (tw, th), (0, 0, 0, 0))
    ImageDraw.Draw(tmp).text((tw / 2, th / 2), yl, fill=tc, font=f, anchor="mm")
    tmp = tmp.rotate(90, expand=True)
    ov.alpha_composite(tmp, (int(8 * s), max(0, int(T + ph / 2 - tw / 2))))
    return ov


_BG = (11, 13, 18)


def annotate(pmap, meta, grid=True, earth=True, width=2000):
    """The map with axes, N/E/S/W, (optionally) faint guide lines and, in solar
    radii, the Earth as a dot to the radial scale."""
    from PIL import Image
    a = (np.clip(pmap, 0, 1) * 255 + 0.5).astype(np.uint8)
    if a.shape[2] == 1:
        a = np.repeat(a, 3, axis=2)
    h0, w0 = a.shape[:2]
    s = 1.0
    pw, ph, L, T, Rm, B = _canvas(h0, w0, width, s)
    img = Image.fromarray(a).resize((pw, ph), Image.LANCZOS)
    can = Image.new("RGBA", (pw + L + Rm, ph + T + B), _BG + (255,))
    can.paste(img, (L, T))
    can.alpha_composite(_overlay(meta, can.size, pw, ph, L, T, s, grid=grid, earth=earth))
    return can.convert("RGB")


def export_layers(pmap, meta, base_path, overlay_path, grid=True, earth=True):
    """THE ANNOTATED MAP AS TWO LAYERS (2026-11-08), for finishing elsewhere:
    the base is the map at its own resolution (no resampling) on the labelled
    canvas, 16-bit, no labels; the overlay is everything drawn on top of it --
    guide lines, frame, numbers, N/E/S/W, the Earth -- on transparency, a
    16-bit RGBA TIFF on the same canvas. Denoise or retouch the base, put the
    overlay on top (Normal), and the pair is the annotated map."""
    import tifffile
    from . import icc
    h0, w0 = pmap.shape[:2]
    s = w0 / 2000.0
    pw, ph, L, T, Rm, B = _canvas(h0, w0, w0, s)
    if ph != h0:                       # a very flat map: keep its pixels, no stretch
        ph = h0
    W, H = pw + L + Rm, ph + T + B
    g = pmap.shape[2] == 1
    base = np.empty((H, W) if g else (H, W, 3), np.float32)
    if g:
        base[:] = (0.2126 * _BG[0] + 0.7152 * _BG[1] + 0.0722 * _BG[2]) / 255.0
        base[T:T + ph, L:L + pw] = np.clip(pmap[..., 0], 0, 1)
    else:
        base[:] = np.array(_BG, np.float32) / 255.0
        base[T:T + ph, L:L + pw] = np.clip(pmap, 0, 1)
    _prof = icc.srgb_profile()
    tifffile.imwrite(base_path, (base * 65535 + 0.5).astype(np.uint16), compression="zlib",
                     photometric="minisblack" if g else "rgb",
                     extratags=[] if g else [(34675, 1, len(_prof), _prof, False)],
                     description="eclipseforgehdr polar map, base layer (no labels)")
    ov = np.asarray(_overlay(meta, (W, H), pw, ph, L, T, s, grid=grid, earth=earth), np.uint16) * 257
    tifffile.imwrite(overlay_path, ov, photometric="rgb", extrasamples=[2], compression="zlib",
                     extratags=[(34675, 1, len(_prof), _prof, False)],
                     description="eclipseforgehdr polar map, overlay layer (labels and guide lines "
                                 "on transparency; put it over the base layer)")
    return W, H


def polar_of(rgb, cy, cx, R_moon, wd, folder, north="solar", outer="frame", step=1.0,
             gray=False, radial="linear", centre="S", mirrored=False, unit_scale=1.0,
             manual_cw=None, offset_deg=0.0):
    """The unwrapped map and its meta, as export() makes it (no files).
    `unit_scale`: for a picture smaller than the stack (the preview, 1/4), the
    factor its pixels are smaller by, so the radii in solar radii stay right."""
    o = orientation(wd, folder, north=north, manual_cw=manual_cw, offset_deg=offset_deg)
    if unit_scale != 1.0 and o.get("px_per_unit"):
        o = dict(o, px_per_unit=o["px_per_unit"] / float(unit_scale))
    src = rgb
    if gray and rgb.ndim == 3:
        src = (0.2126 * rgb[..., 0] + 0.7152 * rgb[..., 1] + 0.0722 * rgb[..., 2])[..., None]
    pmap, meta = unwrap(src, cy, cx, R_moon, o, outer=outer, step=step, radial=radial, centre=centre,
                         mirrored=mirrored)
    return pmap, meta, o


# CAROUSEL TILES (2026-11-10). The map cut into n tiles of a social size
# that lie side by side as one strip: swiped through, it runs on without a
# seam, and as it spans exactly 360 degrees the last tile's right edge is the
# first tile's left edge. The map keeps its shape: scaled so that n tiles
# hold its full width; when it is then taller than the tile leaves room for,
# its outer part (the faintest) is cut off the top, else it sits in the middle
# with black above and below. The angles are along the foot; a label that
# would fall on a cut between two tiles is set beside its tick instead.
TILE_SIZES = {"3:4": (1080, 1440), "4:5": (1080, 1350), "1:1": (1080, 1080), "9:16": (1080, 1920)}


def tile_strip(strip, th, n, fmt, stem, origin="solar", quality=93):
    """`strip` (rows, cols, C) float 0-1, row 0 outermost; `th` the angle of
    each column (deg). Writes n tiles as stem_tile_01_of_0n.jpg ...; returns
    the paths."""
    from PIL import Image, ImageDraw
    from . import extras as _ex
    n = max(2, int(n))
    Tw, Th = TILE_SIZES.get(fmt, TILE_SIZES["3:4"])
    a = np.clip(np.asarray(strip, np.float32), 0, 1)
    if a.ndim == 2:
        a = a[..., None]
    if a.shape[2] == 1:
        a = np.repeat(a, 3, axis=2)
    rows, cols = a.shape[:2]
    Wt = n * Tw
    band = int(Tw * 0.085)                       # the angles under the map
    top_pad = int(Tw * 0.03)
    room = Th - band - top_pad
    k = Wt / float(cols)
    hs = int(round(rows * k))
    if hs > room:                                # too tall: the outer rows go
        keep = int(round(room / k))
        a = a[rows - keep:]
        hs = room
    im = Image.fromarray((a * 255 + 0.5).astype(np.uint8)).resize((Wt, hs), Image.LANCZOS)
    can = Image.new("RGB", (Wt, Th), _ex.BG)
    y0 = top_pad + (room - hs) // 2
    can.paste(im, (0, y0))
    d = _ex.Draw(can)
    yb = y0 + hs
    fs = Tw / 44.0
    f_n, f_l = _ex._font(fs), _ex._font(fs * 1.15, True)
    step = 30 if n <= 6 else 15
    thc = np.asarray(th, np.float64)
    # column of each whole label angle (the angle falls or rises steadily)
    xs = (np.arange(cols) + 0.5) * k
    for lab in range(0, 360, step):
        dd = ((thc - lab + 180.0) % 360.0) - 180.0
        j = np.nonzero(np.sign(dd[:-1]) != np.sign(dd[1:]))[0]
        cand = [i for i in j if abs(dd[i]) < 5 and abs(dd[i + 1]) < 5]
        if abs(dd[0]) < 0.06:
            cand.append(0)
        for i in cand:
            x = float(xs[i] + (xs[min(i + 1, cols - 1)] - xs[i]) * (dd[i] / (dd[i] - dd[min(i + 1, cols - 1)])
                                                                    if dd[i] != dd[min(i + 1, cols - 1)] else 0))
            nesw = {0: "N", 90: "E", 180: "S", 270: "W"}.get(lab) if origin in ("solar", "celestial") else None
            t_ = nesw or f"{lab}\u00b0"
            f_ = f_l if nesw else f_n
            w_ = d.textlength(t_, font=f_)
            # the cut nearest this tick: keep the label on one tile
            cut = round(x / Tw) * Tw
            anc = "mt"
            if abs(x - cut) < w_ / 2 + 4:
                anc = "lt" if x >= cut else "rt"
            d.line([(x, yb), (x, yb + fs * 0.45)], fill=(165, 165, 165), width=1)
            d.text((x + (3 if anc == "lt" else (-3 if anc == "rt" else 0)), yb + fs * 0.7), t_,
                   fill=_ex.GOLD if nesw else (205, 205, 205), font=f_, anchor=anc)
    paths = []
    for i in range(n):
        p_ = f"{stem}_tile_{i + 1:02d}_of_{n:02d}.jpg"
        can.crop((i * Tw, 0, (i + 1) * Tw, Th)).save(p_, quality=quality, subsampling=0)
        paths.append(p_)
    return paths


def export(rgb, cy, cx, R_moon, wd, folder, tif_path, jpg_path, north="solar", outer="frame",
           grid=True, earth=True, step=1.0, gray=False, radial="linear", centre="S",
           mirrored=False, tif_bits=16, manual_cw=None, offset_deg=0.0, tiles=None):
    """Writes the clean TIFF (16-bit, or 8 with tif_bits=8) and the annotated
    picture (JPG, or PNG when jpg_path ends in .png). Either path may be None
    to skip that file. Returns a log line."""
    import tifffile
    from . import icc
    o = orientation(wd, folder, north=north, manual_cw=manual_cw, offset_deg=offset_deg)
    src = rgb
    if gray and rgb.ndim == 3:
        src = (0.2126 * rgb[..., 0] + 0.7152 * rgb[..., 1] + 0.0722 * rgb[..., 2])[..., None]
    pmap, meta = unwrap(src, cy, cx, R_moon, o, outer=outer, step=step, radial=radial, centre=centre,
                         mirrored=mirrored)
    g = pmap.shape[2] == 1
    _prof = icc.srgb_profile()
    if tif_path:
        _mx, _dt = (255, np.uint8) if tif_bits == 8 else (65535, np.uint16)
        tifffile.imwrite(tif_path, ((np.clip(pmap[..., 0] if g else pmap, 0, 1)) * _mx + 0.5).astype(_dt),
                       compression="zlib", photometric="minisblack" if g else "rgb",
                       extratags=[] if g else [(34675, 1, len(_prof), _prof, False)],
                       description="eclipseforgehdr polar map: columns 360-0 deg left to right (%s), rows r %.3f-%.3f %s (%s)"
                       % (meta["origin"], meta["r1"], meta["r0"], "solar radii" if meta["unit"] == "sun"
                          else "Moon radii", radial))
    if jpg_path:
        _im = annotate(pmap, meta, grid=grid, earth=earth)
        if str(jpg_path).lower().endswith(".png"):
            _im.save(jpg_path, icc_profile=_prof)
        else:
            _im.save(jpg_path, quality=92, icc_profile=_prof)
    tile_note = ""
    if tiles and tiles.get("n"):
        stem = os.path.splitext(jpg_path or tif_path)[0]
        nt_ = pmap.shape[1]
        thc = (meta["start"] + meta["dirn"] * (np.arange(nt_) + 0.5) / nt_ * 360.0) % 360.0
        tp = tile_strip(pmap, thc, int(tiles["n"]), str(tiles.get("fmt") or "3:4"), stem, origin=meta["origin"])
        tile_note = f"; {len(tp)} carousel tiles {tiles.get('fmt') or '3:4'}"
    where = {"solar": ("from solar north (P %+.2f deg)" % o["p"]) if manual_cw is None else "from solar north (by hand)",
             "celestial": "from celestial north", "up": "from image up"}[o["origin"]]
    return (("polar map %d x %d, angle %s, r %.2f-%.2f %s (" + radial + "); %s")
            % (pmap.shape[1], pmap.shape[0], where, meta["r0"], meta["r1"],
               "solar radii" if meta["unit"] == "sun" else "Moon radii", o["note"]) + tile_note)
