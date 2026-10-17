"""TOTAL SOLAR ECLIPSES: LOCAL CIRCUMSTANCES AND A CITY LIST (lab 0.25.1).

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
               c2_ut=None, c3_ut=None)
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
