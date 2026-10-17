"""THE OBSERVING SITE AND TIME (lab 0.25.1), optional.

Where and when the bracket was shot: latitude, longitude, height, and the date
and local clock time of totality with its offset from UTC. Every eclipse
chaser knows these, and they make the star catalogue match (starcat.py)
exact instead of searched:

  * the time puts the Sun -- and so the catalogue -- at its true place, so
    the +-14 h clock search and the eclipse-date table are not needed;
  * time and place give the Sun's altitude, and with it the atmospheric
    refraction that lifts the bottom of a low-Sun frame more than its top
    (one test set with the Sun 7 deg up: the stars sat 27-37 px from their
    catalogue places without it, 12-16 px with it).

Saved beside the frames as SITE_FILE, so it survives a cleared cache and
travels with the data. Coordinates are accepted as decimal degrees, degrees
and decimal minutes, or degrees minutes seconds, with a decimal point or a
decimal comma, a sign or a hemisphere letter:
    41.855   41,855   41 51.3   41 51,3   41 51 18   41°51'18"N   N 41 51 18
    -1.922   1 55.32 W   1°55'19" W   W1.922
"""
import os, re, json, math

SITE_FILE = "eclipseforge_site.json"


def parse_angle(text, kind="lat"):
    """Decimal degrees from DD.dddd / DD MM.mmmm / DD MM SS.ss (point or comma
    decimals, sign or N/S/E/W). Raises ValueError with a sentence."""
    s = str(text or "").strip()
    if not s:
        raise ValueError(f"no {'latitude' if kind == 'lat' else 'longitude'} given")
    u = s.upper().replace("\u2212", "-")
    hemi = None
    m1 = re.match(r"^\s*([NSEW])\s*(.*)$", u)
    m2 = re.match(r"^(.*?)\s*([NSEW])\s*$", u)
    if m1:
        hemi, u = m1.group(1), m1.group(2)
    elif m2:
        u, hemi = m2.group(1), m2.group(2)
    if re.search(r"[NSEW]", u):
        raise ValueError(f"'{s}': only one hemisphere letter, at the start or the end")
    sign = 1.0
    u = u.strip()
    if u.startswith("-"):
        sign, u = -1.0, u[1:]
    elif u.startswith("+"):
        u = u[1:]
    for ch in "°º'′’\"″”:;":
        u = u.replace(ch, " ")
    toks = u.split()
    if not 1 <= len(toks) <= 3:
        raise ValueError(f"'{s}': use DD.dddd, DD MM.mmmm or DD MM SS")
    vals = []
    for k, t in enumerate(toks):
        if "," in t and "." not in t:
            t = t.replace(",", ".")
        if not re.fullmatch(r"\d+(\.\d*)?|\.\d+", t):
            raise ValueError(f"'{s}': '{toks[k]}' is not a number")
        if "." in t and k < len(toks) - 1:
            raise ValueError(f"'{s}': only the last part may have decimals")
        vals.append(float(t))
    deg = vals[0]
    if len(vals) >= 2:
        if vals[1] >= 60:
            raise ValueError(f"'{s}': minutes must be below 60")
        deg += vals[1] / 60.0
    if len(vals) == 3:
        if vals[2] >= 60:
            raise ValueError(f"'{s}': seconds must be below 60")
        deg += vals[2] / 3600.0
    if hemi:
        if kind == "lat" and hemi not in "NS":
            raise ValueError(f"'{s}': a latitude is N or S")
        if kind == "lon" and hemi not in "EW":
            raise ValueError(f"'{s}': a longitude is E or W")
        if sign < 0:
            raise ValueError(f"'{s}': a minus sign and a hemisphere letter together")
        if hemi in "SW":
            sign = -1.0
    deg *= sign
    lim = 90.0 if kind == "lat" else 180.0
    if abs(deg) > lim:
        raise ValueError(f"'{s}': beyond {lim:g} degrees")
    return deg


def parse_offset(text):
    """UTC offset in hours from '+2', '-5', '+5:30', '5,5', 'UTC+2'."""
    s = str(text or "").strip().upper().replace("UTC", "").replace("GMT", "").replace("−", "-").strip()
    if not s:
        return 0.0
    m = re.fullmatch(r"([+-]?)(\d{1,2})(?:[:h](\d{2}))?|([+-]?)(\d{1,2}[.,]\d+)", s)
    if not m:
        raise ValueError(f"UTC offset '{text}': use e.g. +2, -5 or +5:30")
    if m.group(2) is not None:
        v = int(m.group(2)) + (int(m.group(3)) / 60.0 if m.group(3) else 0.0)
        v = -v if m.group(1) == "-" else v
    else:
        v = float(m.group(5).replace(",", "."))
        v = -v if m.group(4) == "-" else v
    if abs(v) > 14:
        raise ValueError(f"UTC offset '{text}': beyond 14 hours")
    return v


def ut_from_local(date, time, offset_h):
    """('2026-08-12', '20:29:30', 2.0) -> (Y, M, D, hours UT), day rolled over."""
    d = [int(x) for x in str(date).replace(":", "-").replace(".", "-").split("-")[:3]]
    t = [float(x.replace(",", ".")) for x in str(time).split(":")]
    if len(d) != 3 or not 1 <= len(t) <= 3:
        raise ValueError("date as YYYY-MM-DD and time as HH:MM or HH:MM:SS")
    h = t[0] + (t[1] if len(t) > 1 else 0) / 60.0 + (t[2] if len(t) > 2 else 0) / 3600.0 - offset_h
    import datetime as _dt
    base = _dt.datetime(d[0], d[1], d[2]) + _dt.timedelta(hours=h)
    return base.year, base.month, base.day, base.hour + base.minute / 60.0 + (base.second + base.microsecond / 1e6) / 3600.0


def split_pair(text):
    """'41.864906, -1.978736' (as Google Maps copies it), '41.86 -1.98',
    '41,86; -1,98' or '41°51'54"N 1°58'43"W' -> (lat_text, lon_text), or None
    when the text holds one coordinate only."""
    t = str(text or "").strip()
    for pat in (r"\s*;\s*", r",\s+", r"(?<=\d),(?=[-+]?\d+\.)"):
        parts = [p for p in re.split(pat, t) if p.strip()]
        if len(parts) == 2:
            return parts[0].strip(), parts[1].strip()
    m = re.match(r"^(.*?[NS])\s+(.*[EW])$", t.upper())
    if m:
        return m.group(1).strip(), m.group(2).strip()
    toks = t.split()
    if len(toks) == 2 and all(re.fullmatch(r"[-+]?\d+([.,]\d+)?", x) for x in toks):
        return toks[0], toks[1]
    return None


def site_path(folder):
    return os.path.join(folder, SITE_FILE)


def load(folder):
    """The saved site as entered plus parsed numbers, or None."""
    try:
        s = json.load(open(site_path(folder)))
        return normalise(s)
    except Exception:
        return None


def normalise(raw):
    """Parse what the user typed; returns the dict to save (raises ValueError).
    Each part is optional, but latitude and longitude come together, and a time
    needs its date."""
    out = {"lat_text": str(raw.get("lat_text", raw.get("lat", "")) or "").strip(),
           "lon_text": str(raw.get("lon_text", raw.get("lon", "")) or "").strip(),
           "height_m": raw.get("height_m", raw.get("height", "")),
           "date": str(raw.get("date", "") or "").strip(),
           "time": str(raw.get("time", "") or "").strip(),
           "utc_offset_text": str(raw.get("utc_offset_text", raw.get("utc_offset", "")) or "").strip()}
    if out["lat_text"] and not out["lon_text"]:
        pr = split_pair(out["lat_text"])          # a pasted "lat, lon" pair
        if pr:
            out["lat_text"], out["lon_text"] = pr
    if bool(out["lat_text"]) != bool(out["lon_text"]):
        raise ValueError("give both latitude and longitude, or neither")
    if out["lat_text"]:
        out["lat"] = parse_angle(out["lat_text"], "lat")
        out["lon"] = parse_angle(out["lon_text"], "lon")
    h = out["height_m"]
    if h in ("", None):
        out["height_m"] = None
    else:
        try:
            out["height_m"] = float(str(h).replace(",", "."))
        except ValueError:
            raise ValueError(f"height '{h}': a number of metres")
        if not -500 <= out["height_m"] <= 9000:
            raise ValueError(f"height {out['height_m']:g} m is out of range")
    if out["time"] and not out["date"]:
        raise ValueError("a time needs its date")
    if out["date"]:
        off = parse_offset(out["utc_offset_text"])
        out["utc_offset"] = off
        if out["time"]:
            out["ut"] = list(ut_from_local(out["date"], out["time"], off))
    return out


def save(folder, raw):
    s = normalise(raw)
    with open(site_path(folder), "w") as fh:
        json.dump(s, fh, indent=1)
    return s


def refraction_factor(height_m, temp_c=15.0):
    """Pressure and temperature scaling of the standard refraction formula
    (1010 hPa, 10 C); pressure from the height with an 8.4 km scale height."""
    p = 1010.0 * math.exp(-(height_m or 0.0) / 8400.0)
    return (p / 1010.0) * (283.0 / (273.0 + temp_c))


def sun_altaz(site):
    """Sun altitude (refracted) and azimuth for a site with time and place."""
    if not site or "ut" not in site or "lat" not in site:
        return None
    from . import starcat as _sc
    Y, M, D, h = site["ut"]
    JD = _sc._jd(Y, M, D, h)
    alt, az = _sc.altaz_of_date_j2000(*_sc.sun_j2000(JD), JD, site["lat"], site["lon"])
    alt_r = alt + _sc.refraction_deg(alt) * refraction_factor(site.get("height_m"))
    return float(alt_r), float(az)
