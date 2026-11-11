"""TOTAL SOLAR ECLIPSES: LOCAL CIRCUMSTANCES AND A CITY LIST (0.25.1;
2017-2037 added 2026-11-08, each checked against NASA's greatest-eclipse
duration: within 0.1 s, mid-eclipse within 0.2 s).

For the observing-site form (observer.py): pick the eclipse and a city and
the form fills in latitude, longitude, height, date, local time of
mid-totality and the UTC offset; or type a position and it is checked
against the path of totality.

THE ELEMENTS. Besselian elements of each eclipse, as published by NASA:
"Eclipse Predictions by Fred Espenak, NASA's GSFC" (eclipse.gsfc.nasa.gov,
SEsearch/SEdata.php). Polynomials in t = t1 - t0 (hours, TDT); the cubic
terms of d, l1, l2 and mu are not printed there and are 0 here (they are
below 1e-7 for every eclipse of the century).

THE CALCULATION follows the Explanatory Supplement to the Astronomical
Almanac (1961, ch. 9) as used by Espenak's and O'Byrne's local-circumstance
codes: the observer's place in the fundamental plane (xi, eta, zeta) from
geodetic latitude, longitude and height; mid-eclipse where the shadow axis
passes closest (Newton steps); total where that distance is inside the
umbra's radius on the observer's plane, |L2'|; second and third contact from
the same quadratic; the Sun's altitude from the declination and hour angle.
Lunar limb profile is not modelled (NASA: limits move 1-3 km, durations
1-3 s). Checked against NASA's central-line table (eclipses_test in the
status notes): mid-totality to the second, durations within 1 s.

THE CITIES (eclipse_cities.json, built by tools/build_eclipse_cities.py):
GeoNames (geonames.org, CC BY 4.0) populated places, as Stellarium ships
them (data/base_locations.bin.gz): name, region, population, latitude,
longitude, elevation, IANA time zone; country from the time zone
(tzdata zone.tab). Every place inside a path of totality with at least
POP_MIN inhabitants, plus a few smaller places added by hand.
"""
import math, json, os

ECLIPSES = {
    # 2026 Aug 12, Saros 126, greatest eclipse 17:47:06 TDT, Delta T 75.4 s
    "2026-08-12": dict(t0=18.0, dT=75.4,
                       x=(0.4755140, 0.5189249, -0.0000773, -0.0000080),
                       y=(0.7711830, -0.2301680, -0.0001246, 0.0000038),
                       d=(14.7966700, -0.0120650, -0.0000030),
                       l1=(0.5379550, 0.0000939, -0.0000121),
                       l2=(-0.0081420, 0.0000935, -0.0000121),
                       mu=(88.747787, 15.003090, 0.0),
                       tf1=0.0046141, tf2=0.0045911,
                       name="12 August 2026 (Greenland, Iceland, Spain)"),
    # 2027 Aug 02, Saros 136, greatest eclipse 10:07:50 TDT, Delta T 76.0 s
    "2027-08-02": dict(t0=10.0, dT=76.0,
                       x=(-0.0197720, 0.5447123, -0.0000446, -0.0000092),
                       y=(0.1600610, -0.2111582, -0.0001217, 0.0000038),
                       d=(17.7624702, -0.0101810, -0.0000040),
                       l1=(0.5305960, 0.0000138, -0.0000128),
                       l2=(-0.0154640, 0.0000137, -0.0000128),
                       mu=(328.422546, 15.002100, 0.0),
                       tf1=0.0046064, tf2=0.0045834,
                       name="2 August 2027 (Spain, North Africa, Egypt, Saudi Arabia, Yemen, Somalia)"),
    # The ones below were added 2026-11-08 from NASA's Besselian-element pages
    # (eclipse.gsfc.nasa.gov/SEbeselm, "Eclipse Predictions by Fred Espenak,
    # NASA's GSFC"); there the cubic terms of d, l1, l2 and mu are blank (0).
    # 2017 Aug 21, Saros 145, greatest eclipse 18:26:40.3 TDT, Delta T 68.4 s
    "2017-08-21": dict(t0=18.0, dT=68.4,
                       x=(-0.129576, 0.5406409, -0.0000293, -0.0000081),
                       y=(0.485417, -0.1416394, -0.0000905, 0.0000021),
                       d=(11.86697, -0.013622, -0.000002),
                       l1=(0.542112, 0.0001241, -0.0000118),
                       l2=(-0.004025, 0.0001235, -0.0000117),
                       mu=(89.24545, 15.003937, 0.0),
                       tf1=0.0046222, tf2=0.0045992,
                       name="21 August 2017 (USA)"),
    # 2019 Jul 02, Saros 127, greatest eclipse 19:24:07.5 TDT, Delta T 69.0 s
    "2019-07-02": dict(t0=19.0, dT=69.0,
                       x=(-0.215634, 0.5662070, 0.0000274, -0.0000088),
                       y=(-0.650709, 0.0106401, -0.0001272, -0.0000003),
                       d=(23.01295, -0.003187, -0.000005),
                       l1=(0.537650, -0.0000898, -0.0000120),
                       l2=(-0.008465, -0.0000894, -0.0000120),
                       mu=(103.97974, 14.999507, 0.0),
                       tf1=0.0045984, tf2=0.0045755,
                       name="2 July 2019 (Chile, Argentina)"),
    # 2020 Dec 14, Saros 142, greatest eclipse 16:14:39.4 TDT, Delta T 69.5 s
    "2020-12-14": dict(t0=16.0, dT=69.5,
                       x=(-0.181827, 0.5633549, 0.0000217, -0.0000089),
                       y=(-0.269642, -0.0858121, 0.0001884, 0.0000015),
                       d=(-23.25776, -0.001986, 0.000006),
                       l1=(0.543881, 0.0000970, -0.0000126),
                       l2=(-0.002266, 0.0000966, -0.0000125),
                       mu=(61.26593, 14.996498, 0.0),
                       tf1=0.0047502, tf2=0.0047266,
                       name="14 December 2020 (Chile, Argentina)"),
    # 2021 Dec 04, Saros 152, greatest eclipse 07:34:38.1 TDT, Delta T 69.8 s
    "2021-12-04": dict(t0=8.0, dT=69.8,
                       x=(0.025194, 0.5683012, 0.0000400, -0.0000092),
                       y=(-0.983559, -0.1315136, 0.0002212, 0.0000024),
                       d=(-22.27472, -0.005178, 0.000006),
                       l1=(0.537804, -0.0000160, -0.0000131),
                       l2=(-0.008292, -0.0000160, -0.0000131),
                       mu=(302.45218, 14.997278, 0.0),
                       tf1=0.0047434, tf2=0.0047198,
                       name="4 December 2021 (Antarctica)"),
    # 2023 Apr 20 (hybrid), Saros 129, greatest eclipse 04:17:55.1 TDT, Delta T 70.3 s
    "2023-04-20": dict(t0=4.0, dT=70.3,
                       x=(0.026967, 0.4950173, 0.0000137, -0.0000072),
                       y=(-0.427322, 0.2441985, -0.0000494, -0.0000037),
                       d=(11.41179, 0.013741, -0.000002),
                       l1=(0.546803, 0.0001216, -0.0000116),
                       l2=(0.000662, 0.0001210, -0.0000115),
                       mu=(240.24289, 15.003411, 0.0),
                       tf1=0.0046549, tf2=0.0046318,
                       name="20 April 2023 (Western Australia, Timor-Leste, Indonesia; hybrid)"),
    # 2024 Apr 08, Saros 139, greatest eclipse 18:18:29.0 TDT, Delta T 70.6 s
    "2024-04-08": dict(t0=18.0, dT=70.6,
                       x=(-0.318157, 0.5117105, 0.0000326, -0.0000085),
                       y=(0.219747, 0.2709586, -0.0000594, -0.0000047),
                       d=(7.58620, 0.014844, -0.000002),
                       l1=(0.535813, 0.0000618, -0.0000128),
                       l2=(-0.010274, 0.0000615, -0.0000127),
                       mu=(89.59122, 15.004084, 0.0),
                       tf1=0.0046683, tf2=0.0046450,
                       name="8 April 2024 (Mexico, USA, Canada)"),
    # 2028 Jul 22, Saros 146, greatest eclipse 02:56:39.0 TDT, Delta T 72.1 s
    "2028-07-22": dict(t0=3.0, dT=72.1,
                       x=(-0.154300, 0.5449941, -0.0000226, -0.0000095),
                       y=(-0.586380, -0.1746077, -0.0001022, 0.0000029),
                       d=(20.18231, -0.007974, -0.000005),
                       l1=(0.535236, -0.0000859, -0.0000123),
                       l2=(-0.010847, -0.0000854, -0.0000122),
                       mu=(223.37866, 15.001018, 0.0),
                       tf1=0.0046016, tf2=0.0045786,
                       name="22 July 2028 (Australia, New Zealand)"),
    # 2030 Nov 25, Saros 133, greatest eclipse 06:51:36.5 TDT, Delta T 77.3 s
    "2030-11-25": dict(t0=7.0, dT=77.3,
                       x=(0.044200, 0.5787731, 0.0000157, -0.0000099),
                       y=(-0.392726, -0.0551896, 0.0001744, 0.0000008),
                       d=(-20.76100, -0.007989, 0.000005),
                       l1=(0.538213, -0.0000379, -0.0000130),
                       l2=(-0.007886, -0.0000377, -0.0000130),
                       mu=(288.27457, 14.998360, 0.0),
                       tf1=0.0047361, tf2=0.0047125,
                       name="25 November 2030 (Namibia, South Africa, Lesotho)"),
    # 2033 Mar 30, Saros 120, greatest eclipse 18:02:35.3 TDT, Delta T 78.9 s
    "2033-03-30": dict(t0=18.0, dT=78.9,
                       x=(-0.318785, 0.5554225, 0.0000227, -0.0000094),
                       y=(0.924678, 0.1756610, -0.0000800, -0.0000029),
                       d=(4.09368, 0.015719, -0.000001),
                       l1=(0.534941, 0.0000276, -0.0000130),
                       l2=(-0.011141, 0.0000275, -0.0000129),
                       mu=(88.92808, 15.004456, 0.0),
                       tf1=0.0046807, tf2=0.0046574,
                       name="30 March 2033 (Alaska, Russian Far East)"),
    # 2034 Mar 20, Saros 130, greatest eclipse 10:18:44.9 TDT, Delta T 79.5 s
    "2034-03-20": dict(t0=10.0, dT=79.5,
                       x=(-0.259593, 0.5481665, 0.0000230, -0.0000094),
                       y=(0.220828, 0.1755788, -0.0000080, -0.0000028),
                       d=(-0.05513, 0.016042, -0.000000),
                       l1=(0.538630, -0.0000665, -0.0000127),
                       l2=(-0.007470, -0.0000662, -0.0000126),
                       mu=(328.13922, 15.004398, 0.0),
                       tf1=0.0046951, tf2=0.0046718,
                       name="20 March 2034 (Nigeria, Cameroon, Chad, Saudi Arabia, Kuwait, Iran, Afghanistan, Pakistan, India)"),
    # 2035 Sep 02, Saros 145, greatest eclipse 01:56:45.9 TDT, Delta T 80.6 s
    "2035-09-02": dict(t0=2.0, dT=80.6,
                       x=(0.134339, 0.5377713, -0.0000358, -0.0000080),
                       y=(0.348988, -0.1584657, -0.0000595, 0.0000023),
                       d=(8.01771, -0.014783, -0.000002),
                       l1=(0.541919, 0.0001104, -0.0000119),
                       l2=(-0.004197, 0.0001098, -0.0000118),
                       mu=(210.02995, 15.004648, 0.0),
                       tf1=0.0046328, tf2=0.0046097,
                       name="2 September 2035 (China, North Korea, Japan)"),
    # 2037 Jul 13, Saros 127, greatest eclipse 02:40:35.3 TDT, Delta T 81.9 s
    "2037-07-13": dict(t0=3.0, dT=81.9,
                       x=(0.141589, 0.5635982, -0.0000003, -0.0000087),
                       y=(-0.733774, -0.0318211, -0.0001131, 0.0000003),
                       d=(21.78243, -0.006047, -0.000005),
                       l1=(0.538383, -0.0001102, -0.0000120),
                       l2=(-0.007716, -0.0001096, -0.0000119),
                       mu=(223.55016, 15.000225, 0.0),
                       tf1=0.0045993, tf2=0.0045764,
                       name="13 July 2037 (Australia, New Zealand)"),
}
CITY_FILE = os.path.join(os.path.dirname(__file__), "eclipse_cities.json")
D2R = math.pi / 180.0


def _poly(c, t):
    return sum(a * t ** k for k, a in enumerate(c))


def _dpoly(c, t):
    return sum(k * a * t ** (k - 1) for k, a in enumerate(c) if k)


def _state(E, t, rsp, rcp, lon_e):
    """Shadow geometry relative to the observer at t (hours from t0, TDT)."""
    x, y = _poly(E["x"], t), _poly(E["y"], t)
    dx, dy = _dpoly(E["x"], t), _dpoly(E["y"], t)
    d = _poly(E["d"], t) * D2R
    dd = _dpoly(E["d"], t) * D2R
    mu = _poly(E["mu"], t)
    dmu = _dpoly(E["mu"], t) * D2R
    # hour angle of the shadow axis at the observer; mu is on the ephemeris
    # time scale, the 1.002738 * 15 deg/h * dT/3600 term brings it to UT
    H = (mu + lon_e - 0.00417807 * E["dT"]) * D2R
    xi = rcp * math.sin(H)
    eta = rsp * math.cos(d) - rcp * math.cos(H) * math.sin(d)
    zeta = rsp * math.sin(d) + rcp * math.cos(H) * math.cos(d)
    dxi = dmu * rcp * math.cos(H)
    deta = dmu * xi * math.sin(d) - zeta * dd
    u, v = x - xi, y - eta
    a, b = dx - dxi, dy - deta
    L1 = _poly(E["l1"], t) - zeta * E["tf1"]
    L2 = _poly(E["l2"], t) - zeta * E["tf2"]
    return u, v, a, b, L1, L2, zeta, d, H


def local(date, lat, lon, height_m=0.0):
    """Local circumstances at a place (geodetic lat, lon east-positive, deg).
    Returns dict: total (bool), mid_ut (hours UT on the eclipse date),
    c2_ut, c3_ut, duration_s (0 when not total), sun_alt (deg, geometric,
    at mid-eclipse), magnitude, obscured (bool: any eclipse at all)."""
    E = ECLIPSES[date]
    phi = lat * D2R
    U = math.atan(0.99664719 * math.tan(phi))
    h = (height_m or 0.0) / 6378140.0
    rsp = 0.99664719 * math.sin(U) + h * math.sin(phi)
    rcp = math.cos(U) + h * math.cos(phi)
    t = 0.0
    for _ in range(8):                   # mid-eclipse: closest approach
        u, v, a, b, L1, L2, zeta, d, H = _state(E, t, rsp, rcp, lon)
        n2 = a * a + b * b
        tau = -(u * a + v * b) / n2
        t += tau
        if abs(tau) < 1e-7:
            break
    u, v, a, b, L1, L2, zeta, d, H = _state(E, t, rsp, rcp, lon)
    m = math.hypot(u, v)
    sun_alt = math.degrees(math.asin(max(-1.0, min(1.0, math.sin(d) * math.sin(phi)
                                                    + math.cos(d) * math.cos(phi) * math.cos(H)))))
    out = dict(mid_ut=t + E["t0"] - E["dT"] / 3600.0, sun_alt=sun_alt,
               magnitude=(L1 - m) / (L1 + L2) if (L1 + L2) else 0.0,
               obscured=m < L1 and zeta > 0, total=False, duration_s=0.0,
               c2_ut=None, c3_ut=None, c1_ut=None, c4_ut=None)
    if m < L1 and zeta > 0:
        # first and fourth contacts: |shadow axis - observer| = L1 (the
        # penumbra's edge), the same refinement as for the umbra below
        cs = []
        for sgn in (-1.0, 1.0):
            tc = t
            for _ in range(8):
                u, v, a, b, L1c, L2c, zc, dc, Hc = _state(E, tc, rsp, rcp, lon)
                n = math.hypot(a, b)
                S = (a * v - u * b) / (n * L1c)
                tau = -(u * a + v * b) / (n * n) + sgn * abs(L1c / n * math.sqrt(max(0.0, 1 - S * S)))
                tc += tau
                if abs(tau) < 1e-8:
                    break
            cs.append(tc)
        out.update(c1_ut=cs[0] + E["t0"] - E["dT"] / 3600.0, c4_ut=cs[1] + E["t0"] - E["dT"] / 3600.0)
    if L2 < 0 and m < abs(L2) and zeta > 0:
        # second and third contacts: |shadow axis - observer| = |L2'|, refined
        cs = []
        for sgn in (-1.0, 1.0):
            tc = t
            for _ in range(6):
                u, v, a, b, L1c, L2c, zc, dc, Hc = _state(E, tc, rsp, rcp, lon)
                n = math.hypot(a, b)
                S = (a * v - u * b) / (n * L2c)
                tau = -(u * a + v * b) / (n * n) + sgn * abs(L2c / n * math.sqrt(max(0.0, 1 - S * S)))
                tc += tau
                if abs(tau) < 1e-8:
                    break
            cs.append(tc)
        c2, c3 = cs
        out.update(total=True, c2_ut=c2 + E["t0"] - E["dT"] / 3600.0,
                   c3_ut=c3 + E["t0"] - E["dT"] / 3600.0,
                   duration_s=(c3 - c2) * 3600.0,
                   mid_ut=0.5 * (c2 + c3) + E["t0"] - E["dT"] / 3600.0)
    return out


def _obs(lat):
    phi = lat * D2R
    U = math.atan(0.99664719 * math.tan(phi))
    return 0.99664719 * math.sin(U), math.cos(U)


def _ground(E, t, off):
    """The point on the ground (geodetic lat, lon; sea level) where the shadow
    axis, moved by `off` (fraction of |L2'| across the track: 0 = centre line,
    +1 / -1 = the two limits), meets the Earth at time t; None where it does
    not, where the Sun is below the horizon, or where the shadow is not
    umbral. Newton on the same geometry as local()."""
    x, y = _poly(E["x"], t), _poly(E["y"], t)
    if off:   # start from the offset point in the fundamental plane
        dx, dy = _dpoly(E["x"], t), _dpoly(E["y"], t)
        n0 = math.hypot(dx, dy) or 1.0
        r0 = off * abs(_poly(E["l2"], t))
        x, y = x + r0 * dy / n0, y - r0 * dx / n0
    if x * x + y * y >= 0.9999:
        return None
    d = _poly(E["d"], t) * D2R
    mu = _poly(E["mu"], t)
    # start: the spherical Earth
    z = math.sqrt(1 - x * x - y * y)
    sp = y * math.cos(d) + z * math.sin(d)
    H = math.atan2(x, z * math.cos(d) - y * math.sin(d))
    lat = math.degrees(math.atan(math.tan(math.asin(max(-1, min(1, sp)))) / 0.99330546))
    lon = (math.degrees(H) - mu + 0.00417807 * E["dT"] + 540.0) % 360.0 - 180.0

    def f(la, lo):
        rsp, rcp = _obs(la)
        u, v, a, b, L1, L2, zeta, dd, HH = _state(E, t, rsp, rcp, lo)
        n = math.hypot(a, b) or 1.0
        nx, ny = -b / n, a / n
        r = off * abs(L2)
        return u - r * nx, v - r * ny, L2, zeta
    for _ in range(25):
        fu, fv, L2, zeta = f(lat, lon)
        if abs(fu) < 1e-9 and abs(fv) < 1e-9:
            break
        e = 1e-4
        a1, b1, _, _ = f(lat + e, lon)
        a2, b2, _, _ = f(lat, lon + e)
        j11, j21 = (a1 - fu) / e, (b1 - fv) / e
        j12, j22 = (a2 - fu) / e, (b2 - fv) / e
        det = j11 * j22 - j12 * j21
        if abs(det) < 1e-14:
            return None
        dla = (fu * j22 - fv * j12) / det
        dlo = (fv * j11 - fu * j21) / det
        big = max(abs(dla), abs(dlo))
        if big > 3.0:             # grazing near sunrise/sunset: damp the step
            dla, dlo = dla * 3.0 / big, dlo * 3.0 / big
        lat -= dla
        lon -= dlo
        if abs(lat) > 89.9:
            return None
    fu, fv, L2, zeta = f(lat, lon)
    if abs(fu) > 1e-6 or abs(fv) > 1e-6 or zeta <= 0 or L2 >= 0:
        return None
    return round(lat, 4), round((lon + 540.0) % 360.0 - 180.0, 4)


_PATHS = {}


def path(date, step_min=1.0):
    """The path of totality for the map: {"central": [[lat, lon], ...],
    "north": [...], "south": [...]} in time order, sea level. Limits follow
    the umbra's edge across the track at each instant."""
    if date in _PATHS:
        return _PATHS[date]
    E = ECLIPSES[date]
    out = {"central": [], "north": [], "south": [], "date": date, "name": E["name"]}
    # each line on the 1-minute grid, then refined where neighbours lie more
    # than ~4 km apart: near sunrise and sunset the shadow races over the
    # ground (several km/s at a low Sun) and a minute is hundreds of km there
    ts = _times(E, step_min)
    tc = [t for t in ts if _ground(E, t, 0.0)]
    out["central"] = _track(E, 0.0, tc)
    a_, b_ = _track(E, 1.0, tc), _track(E, -1.0, tc)
    # which limit is the northern one: the side with the higher mean latitude
    ma = sum(p[0] for p in a_) / len(a_) if a_ else 0
    mb = sum(p[0] for p in b_) / len(b_) if b_ else 0
    north_off = 1.0 if ma >= mb else -1.0
    out["north"], out["south"] = (a_, b_) if ma >= mb else (b_, a_)
    # the two ends of the band, where the shadow first touches and last leaves
    # the Earth (sunrise and sunset along the path): for offsets across the
    # track, the first and the last moment the point exists, to about 1 s.
    # Ordered from the north limit to the south limit.
    out["start_cap"], out["end_cap"] = [], []
    if tc:
        for k in range(21):
            f = north_off * (1.0 - k / 10.0)
            for first, key in ((True, "start_cap"), (False, "end_cap")):
                t_edge = _edge(E, f, tc[0] if first else tc[-1], first, step_min)
                if t_edge is not None:
                    g = _ground(E, t_edge, f)
                    if g:
                        out[key].append(list(g))
    _PATHS[date] = out
    return out


def _km(p, q):
    la1, la2 = math.radians(p[0]), math.radians(q[0])
    dl = math.radians(q[1] - p[1])
    c = math.sin(la1) * math.sin(la2) + math.cos(la1) * math.cos(la2) * math.cos(dl)
    return 6371.0 * math.acos(max(-1.0, min(1.0, c)))


def _track(E, off, ts, dev_km=0.15, max_km=40.0):
    """One line of the band (offset `off`) at the times `ts`, with points
    added by halving the time step wherever the line bends away from the
    straight segment by more than `dev_km`, or two neighbours are more than
    `max_km` apart (down to a fraction of a second)."""
    pts = [(t, _ground(E, t, off)) for t in ts]
    pts = [(t, g) for t, g in pts if g]
    if not pts:
        return []

    def fill(t1, g1, t2, g2, depth):
        if depth >= 10:
            return []
        L = _km(g1, g2)
        if L <= 0.5:
            return []
        tm = 0.5 * (t1 + t2)
        gm = _ground(E, tm, off)
        if not gm:
            return []
        dl = (g2[1] - g1[1] + 540.0) % 360.0 - 180.0
        mid = (0.5 * (g1[0] + g2[0]), g1[1] + 0.5 * dl)
        if L <= max_km and _km(gm, mid) <= dev_km:
            return []
        return fill(t1, g1, tm, gm, depth + 1) + [gm] + fill(tm, gm, t2, g2, depth + 1)
    out = []
    for (t1, g1), (t2, g2) in zip(pts, pts[1:]):
        out.append(list(g1))
        out.extend(list(g) for g in fill(t1, g1, t2, g2, 0))
    out.append(list(pts[-1][1]))
    return out


def _times(E, step_min):
    t, ts = -3.5, []
    while t <= 3.5:
        ts.append(t)
        t += step_min / 60.0
    return ts


def _edge(E, f, t_c, first, step_min):
    """The first (or last) moment, near t_c, at which the offset-f point of the
    band is on the sunlit Earth; bisection to about a second."""
    st = step_min / 60.0
    grid = [t_c + i * st for i in range(-40, 41)]
    ok = [_ground(E, t, f) is not None for t in grid]
    if not any(ok):
        return None
    if first:
        i = ok.index(True)
        if i == 0:
            return grid[0]
        lo, hi = grid[i - 1], grid[i]          # lo invalid, hi valid
    else:
        i = len(ok) - 1 - ok[::-1].index(True)
        if i == len(ok) - 1:
            return grid[-1]
        lo, hi = grid[i + 1], grid[i]          # lo invalid, hi valid
    for _ in range(7):
        m = (lo + hi) / 2.0
        if _ground(E, m, f) is not None:
            hi = m
        else:
            lo = m
    return hi


K_MOON = 0.2725076        # the Moon's radius in Earth radii (the elements' unit)


def moon_rate(date, lat, lon, height_m=0.0):
    """How fast the Moon moves against the Sun as seen from this place at
    mid-eclipse, in MOON RADII per second, or None when the eclipse is not
    in the list or not seen there.

    On the fundamental plane the observer sits (u, v) from the shadow axis,
    which is the Sun-Moon line, and the Sun-Moon separation seen from there is
    that distance over the observer's distance to the Moon -- as is the Moon's
    own radius K_MOON over the same distance. So in Moon radii the separation
    is simply hypot(u, v) / K_MOON, and its rate hypot(u', v') / K_MOON, with
    the observer's own motion (Earth's rotation) already in u', v'. Times the
    Moon's radius in pixels this is the lunar track's expected speed on the
    sensor, for any lens and camera."""
    E = ECLIPSES.get(str(date))
    if E is None:
        return None
    try:
        loc = local(date, lat, lon, height_m or 0.0)
    except Exception:
        return None
    if not loc.get("obscured"):
        return None
    phi = lat * D2R
    U = math.atan(0.99664719 * math.tan(phi))
    h = (height_m or 0.0) / 6378140.0
    rsp = 0.99664719 * math.sin(U) + h * math.sin(phi)
    rcp = math.cos(U) + h * math.cos(phi)
    t = loc["mid_ut"] - E["t0"] + E["dT"] / 3600.0
    u, v, a, b, L1, L2, zeta, d, H = _state(E, t, rsp, rcp, lon)
    return math.hypot(a, b) / 3600.0 / K_MOON


def utc_offset_h(tz, date, ut_h):
    """UTC offset (hours) of an IANA zone at that moment, or None."""
    try:
        from zoneinfo import ZoneInfo
        import datetime as _dt
        Y, M, D = (int(x) for x in date.split("-"))
        when = _dt.datetime(Y, M, D, tzinfo=_dt.timezone.utc) + _dt.timedelta(hours=ut_h)
        return when.astimezone(ZoneInfo(tz)).utcoffset().total_seconds() / 3600.0
    except Exception:
        return None


def hms(h):
    s = int(round((h % 24.0) * 3600.0))
    return f"{s // 3600:02d}:{(s // 60) % 60:02d}:{s % 60:02d}"


def cities(date):
    """The bundled city list for one eclipse (sorted west to east along the path)."""
    try:
        return json.load(open(CITY_FILE)).get(date, [])
    except Exception:
        return []


def check_site(site):
    """For a saved site (observer.normalise): the eclipse on its date, and
    whether the place sees totality, with local times. None when the date is
    not a bundled eclipse or there is no position."""
    if not site or site.get("lat") is None or not site.get("date"):
        return None
    date = str(site["date"])[:10]
    if date not in ECLIPSES:
        return None
    r = local(date, float(site["lat"]), float(site["lon"]), float(site.get("height_m") or 0.0))
    r["date"] = date
    return r
