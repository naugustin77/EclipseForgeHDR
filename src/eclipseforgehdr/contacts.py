"""WHERE THE SHOOT SITS IN TOTALITY (lab 0.28.1, 2026-10-18).

With the observing site known (observer.py), eclipses.local gives second and
third contact and mid-totality at that place to the second. The frames carry
the camera's clock, which is local time, UTC or whatever the owner set, so the
offset between the camera and UT is found from the frames themselves: the
clock offset, in whole hours (the site's own UTC offset and UTC tried first),
that puts the middle of the shoot inside totality. There is at most one such
offset, totality being minutes long.

What it is for:
  * the timeline strip marks C2, mid-totality and C3 beside the frames, so a
    sequence is seen at once to cover the start, the middle or the end;
  * on a repeated-bracket set the Moon of ONE instant is stacked into every
    tier (pipeline, per-frame Moon); that instant is mid-totality when the
    shoot spans it, the middle of the shoot otherwise. Mid-totality is the
    symmetric choice -- the prominences on both limbs as the eye saw them in
    the middle of the eclipse -- and it does not move when a bracket at the
    start or the end is added or left out.
A camera clock is rarely right to the second, and the marks are only as good
as the clock. CONTACT FRAMES fix that: the diamond ring and Baily's beads are
shot at second and third contact, by everyone, and a subfolder whose name has
"diamond", "ring", "bead" or "contact" in it is read for them. The last of
them before the bracket is taken as C2, the first after it as C3, their middle
as mid-totality on the camera's clock, and the predicted C2/C3 are placed
around that middle (the ring frames themselves are shot a few seconds outside
totality, so the predicted duration is kept). With one contact frame only,
that contact sets the clock. Without any, the whole hours are searched as
above and the marks carry the camera clock's error, which the log states when
frames fall outside the predicted totality.
"""
import os
import datetime as _dt

CONTACT_DIR_WORDS = ("diamond", "ring", "bead", "contact")


def find_contact_dir(folder):
    """A subfolder of contact frames (diamond ring, Baily's beads), or None."""
    try:
        for name in sorted(os.listdir(folder)):
            p = os.path.join(folder, name)
            n = name.lower()
            if (os.path.isdir(p) and any(w in n for w in CONTACT_DIR_WORDS)
                    and not n.startswith("eclipseforge")):
                return p
    except Exception:
        pass
    return None


def contact_stamps(folder, wd):
    """Time stamps of the contact frames, [] when there are none."""
    d = find_contact_dir(folder)
    if not d:
        return []
    try:
        from .raw import list_raws
        from .timeline import _exif_table
        raws = list_raws(d)
        tab = _exif_table(d, raws, wd)
        return [tab[os.path.basename(p)][1] for p in raws
                if tab.get(os.path.basename(p), (None, None))[1]]
    except Exception:
        return []


def _cam_hours(ts, date):
    """Camera clock of an EXIF/FITS stamp as hours from 00:00 of `date`
    (negative or beyond 24 when the clock is on the previous or next day)."""
    try:
        d, tm = str(ts).replace("T", " ").split(" ", 1)
        Y, Mo, D = (int(x) for x in d.replace("-", ":").split(":"))
        h, mi, se = (float(x) for x in tm.split(":")[:3])
        y0, m0, d0 = (int(x) for x in str(date).replace(":", "-").split("-")[:3])
        days = (_dt.date(Y, Mo, D) - _dt.date(y0, m0, d0)).days
        return days * 24.0 + h + mi / 60.0 + se / 3600.0
    except Exception:
        return None


def shoot_contacts(site, stamps, contacts=None, slack_s=90.0):
    """site: observer.load() dict (date, lat, lon, height_m, utc_offset);
    stamps: the frames' time stamps; contacts: the contact frames' stamps.
    Returns None, or a dict:
      c2_s, mid_s, c3_s   seconds from the first frame
      duration_s          totality at the site
      offset_h            camera clock minus UT, hours
      offset_from         how it was found ("contacts", "site", "utc", "search")
      contact_c2_s, contact_c3_s, contact_span_s   the contact frames used
      first_ut_h          UT of the first frame, hours on the eclipse date
    """
    if not site or site.get("lat") is None or not site.get("date"):
        return None
    try:
        from . import eclipses as _ecl
        loc = _ecl.local(site["date"], site["lat"], site["lon"], site.get("height_m") or 0.0)
    except Exception:
        return None
    if not loc.get("total") or loc.get("c2_ut") is None:
        return None
    cam = [v for v in (_cam_hours(t, site["date"]) for t in stamps) if v is not None]
    if not cam:
        return None
    cam.sort()
    t_first, t_mid = cam[0], 0.5 * (cam[0] + cam[-1])
    c2, c3 = float(loc["c2_ut"]), float(loc["c3_ut"])
    sl = slack_s / 3600.0
    # the contact frames, when there are any: the clock from them
    cc = sorted(v for v in (_cam_hours(t, site["date"]) for t in (contacts or [])) if v is not None)
    before = [v for v in cc if v < cam[0]]
    after = [v for v in cc if v > cam[-1]]
    c2cam = max(before) if before else None
    c3cam = min(after) if after else None
    if c2cam is not None or c3cam is not None:
        if c2cam is not None and c3cam is not None:
            off = 0.5 * (c2cam + c3cam) - 0.5 * (c2 + c3)
        elif c2cam is not None:
            off = c2cam - c2
        else:
            off = c3cam - c3
        ut_first = t_first - off
        over = max(0.0, (c2 - ut_first) * 3600.0, (cam[-1] - off - c3) * 3600.0)
        return {"outside_s": round(over, 1),
                "c2_s": round((c2 - ut_first) * 3600.0, 1),
                "mid_s": round((0.5 * (c2 + c3) - ut_first) * 3600.0, 1),
                "c3_s": round((c3 - ut_first) * 3600.0, 1),
                "duration_s": round((c3 - c2) * 3600.0, 1),
                "shoot_s": round((cam[-1] - cam[0]) * 3600.0, 1),
                "offset_h": round(off, 6), "offset_from": "contacts",
                "contact_c2_s": (round((c2cam - t_first) * 3600.0, 1) if c2cam is not None else None),
                "contact_c3_s": (round((c3cam - t_first) * 3600.0, 1) if c3cam is not None else None),
                "contact_span_s": (round((c3cam - c2cam) * 3600.0, 1)
                                   if c2cam is not None and c3cam is not None else None),
                "first_ut_h": round(ut_first, 6)}
    cands = []
    so = site.get("utc_offset")
    if so is not None:
        cands.append((float(so), "site"))
    cands.append((0.0, "utc"))
    for k in sorted(range(-12, 15), key=lambda k: abs(k - (so or 0.0))):
        cands.append((float(k), "search"))
    seen = set()
    for off, how in cands:
        if off in seen:
            continue
        seen.add(off)
        ut_mid = t_mid - off
        if c2 - sl <= ut_mid <= c3 + sl:
            ut_first = t_first - off
            # frames outside the predicted totality: the camera clock (or the
            # site) is off by about that much; said, not corrected
            over = max(0.0, (c2 - ut_first) * 3600.0, (cam[-1] - off - c3) * 3600.0)
            return {"outside_s": round(over, 1),
                    "c2_s": round((c2 - ut_first) * 3600.0, 1),
                    "mid_s": round((0.5 * (c2 + c3) - ut_first) * 3600.0, 1),
                    "c3_s": round((c3 - ut_first) * 3600.0, 1),
                    "duration_s": round((c3 - c2) * 3600.0, 1),
                    "shoot_s": round((cam[-1] - cam[0]) * 3600.0, 1),
                    "offset_h": off, "offset_from": how,
                    "first_ut_h": round(ut_first, 6)}
    return None


def _hms(h):
    sgn = "-" if h < 0 else "+"
    s = int(round(abs(h) * 3600.0))
    return f"{sgn}{s // 3600}:{(s % 3600) // 60:02d}:{s % 60:02d}"


def offset_text(ct):
    off = ct["offset_h"]
    s = f"UTC{off:+g}" if off else "UTC"
    if ct["offset_from"] == "contacts":
        n = sum(1 for k in ("contact_c2_s", "contact_c3_s") if ct.get(k) is not None)
        return (f"camera clock UTC{_hms(off)}, set from {n} contact frame{'s' if n > 1 else ''}"
                + (f" ({ct['contact_span_s']:.0f} s apart)" if ct.get("contact_span_s") else ""))
    return {"site": f"camera clock on the site's time ({s})",
            "utc": "camera clock on UTC",
            "search": f"camera clock {s}, found from the frames"}[ct["offset_from"]]
