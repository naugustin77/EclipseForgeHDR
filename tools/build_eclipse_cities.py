"""Build src/eclipseforgehdr/eclipse_cities.json: the places inside each
bundled path of totality (eclipses.ECLIPSES), for the site form's city list.

Source: GeoNames (geonames.org, CC BY 4.0) populated places as Stellarium
ships them, data/base_locations.bin.gz (a Qt QDataStream: QMap<QString,
StelLocation>). Country from the IANA time zone via tzdata's zone.tab and
iso3166.tab.

Kept: every place where eclipses.local() says total, with population >=
POP_MIN, every regional or national capital, plus the TOP_PER_COUNTRY largest of each country (so small
countries in the path are not empty), plus places asked for by name (ALWAYS).

    python3 tools/build_eclipse_cities.py /Applications/Stellarium.app/Contents/Resources/data/base_locations.bin.gz
"""
import gzip, struct, json, sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from eclipseforgehdr import eclipses as E

POP_MIN = 20000                 # 50000 until 2026-10-10: Soria (39 000, 1:42 of totality) was missing
TOP_PER_COUNTRY = 3
ALWAYS = {("Thuwal", "Asia/Riyadh"), ("Al Lith", "Asia/Riyadh")}   # 2026-10-10


def read_stellarium(path):
    b = gzip.open(path).read()
    p = [0]
    def take(n):
        v = b[p[0]:p[0] + n]; p[0] += n; return v
    def u32(): return struct.unpack(">I", take(4))[0]
    def i32(): return struct.unpack(">i", take(4))[0]
    def f64(): return struct.unpack(">d", take(8))[0]
    def qs():
        n = u32()
        return None if n == 0xFFFFFFFF else take(n).decode("utf-16-be")
    out = []
    for _ in range(u32()):
        qs(); name, state, region = qs(), qs(), qs()
        role = take(2).decode("utf-16-be")
        pop, lat, lon, alt, _lp = i32(), f64(), f64(), i32(), i32()
        tz, planet, _land = qs(), qs(), qs(); take(1)
        out.append(dict(name=name, state=state or "", region=region or "", role=role, pop=pop,
                        lat=lat, lon=lon, alt=alt, tz=tz or "", planet=planet or ""))
    return out


def countries():
    zi = "/usr/share/zoneinfo"
    names = {}
    for line in open(os.path.join(zi, "iso3166.tab"), encoding="utf-8"):
        if line.strip() and not line.startswith("#"):
            cc, nm = line.rstrip("\n").split("\t")[:2]; names[cc] = nm
    tz2cc = {}
    for fn in ("zone.tab", "zone1970.tab"):
        for line in open(os.path.join(zi, fn), encoding="utf-8"):
            if line.strip() and not line.startswith("#"):
                f = line.rstrip("\n").split("\t")
                for cc in f[0].split(","):
                    tz2cc.setdefault(f[2], cc)
    return names, tz2cc


def main(src):
    L = [o for o in read_stellarium(src) if o["planet"] == "Earth" and o["pop"] > 0 and o["tz"] != "UTC"]
    names, tz2cc = countries()
    out = {}
    for date in E.ECLIPSES:
        inpath = []
        for o in L:
            r = E.local(date, o["lat"], o["lon"], max(o["alt"], 0))
            if not r["total"]:
                continue
            cc = tz2cc.get(o["tz"], "")
            off = E.utc_offset_h(o["tz"], date, r["mid_ut"])
            inpath.append(dict(name=o["name"], state=o["state"], country=names.get(cc, cc), cc=cc,
                               lat=round(o["lat"], 5), lon=round(o["lon"], 5), height_m=int(o["alt"]),
                               tz=o["tz"], pop=o["pop"], utc_offset=off,
                               mid_ut=E.hms(r["mid_ut"]), mid_local=E.hms(r["mid_ut"] + (off or 0.0)),
                               duration_s=round(r["duration_s"], 1), sun_alt=round(r["sun_alt"], 1),
                               _t=r["mid_ut"], _cap=o["role"] in ("R", "C")))
        keep = set()
        bycc = {}
        for c in inpath:
            bycc.setdefault(c["cc"], []).append(c)
        for cc, lst in bycc.items():
            for c in sorted(lst, key=lambda c: -c["pop"])[:TOP_PER_COUNTRY]:
                keep.add(id(c))
        sel = [c for c in inpath if c["pop"] >= POP_MIN or id(c) in keep or c.get("_cap")
               or (c["name"], c["tz"]) in ALWAYS]
        # countries in the order the shadow reaches them, places by name
        first = {cc: min(c["_t"] for c in lst) for cc, lst in bycc.items()}
        sel.sort(key=lambda c: (first[c["cc"]], c["country"], c["name"]))
        for c in sel:
            c.pop("_t"); c.pop("_cap")
        out[date] = sel
        print(date, len(inpath), "in the path,", len(sel), "kept")
    dst = os.path.join(os.path.dirname(E.__file__), "eclipse_cities.json")
    json.dump(out, open(dst, "w"), ensure_ascii=False, separators=(",", ":"))
    print("written", dst, os.path.getsize(dst), "bytes")


if __name__ == "__main__":
    main(sys.argv[1])
