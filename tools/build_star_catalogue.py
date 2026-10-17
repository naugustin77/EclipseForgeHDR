#!/usr/bin/env python3
"""Rebuild src/eclipseforgehdr/star_catalogue.npz from the CDS source files.

The bundled catalogue holds every star to V 8.5 within 12 degrees of the
ecliptic (the Sun is always on it, so these are the only stars that can stand
in an eclipse picture), with its J2000 position and proper motion, V magnitude
and catalogue numbers. Sources, all from the CDS (cdsarc.cds.unistra.fr):

  I/239  hip_main.dat   The Hipparcos Main Catalogue (ESA 1997): position
                        (ICRS, epoch J1991.25), proper motion, V, HD number
  V/50   catalog(.gz)   The Bright Star Catalogue, 5th rev. ed. (Hoffleit &
                        Warren 1991): HR number, Bayer/Flamsteed name, SAO
  I/131A sao.dat(.gz)   The SAO Star Catalog (SAO Staff 1966): SAO number by HD
  I/239  tyc_main.dat   The Tycho Catalogue (ESA 1997): the stars to V 8.5 that
                        Hipparcos did not observe (optional)

Usage:
  python3 tools/build_star_catalogue.py HIP_MAIN BSC_CATALOG [SAO_DAT] [--tycho TYC_MAIN] OUT.npz

tyc_main.dat is 350 MB; only its lines with V <= 8.5 are used, so a file
holding just those lines works the same.

Gzipped or plain files both work. Columns are read at the byte positions of
each catalogue's ReadMe.
"""
import gzip
import sys

import numpy as np

MAG_LIMIT = 8.5
ECL_LAT_MAX = 12.0                  # degrees from the ecliptic
OBLIQUITY = np.radians(23.4392911)  # J2000
EPOCH_HIP = 1991.25

GREEK = {"Alp": "α", "Bet": "β", "Gam": "γ", "Del": "δ", "Eps": "ε", "Zet": "ζ",
         "Eta": "η", "The": "θ", "Iot": "ι", "Kap": "κ", "Lam": "λ", "Mu": "μ",
         "Nu": "ν", "Xi": "ξ", "Omi": "ο", "Pi": "π", "Rho": "ρ", "Sig": "σ",
         "Tau": "τ", "Ups": "υ", "Phi": "φ", "Chi": "χ", "Psi": "ψ", "Ome": "ω"}
SUP = {"1": "¹", "2": "²", "3": "³", "4": "⁴", "5": "⁵", "6": "⁶", "7": "⁷", "8": "⁸", "9": "⁹"}


def _open(path):
    with open(path, "rb") as f:
        gz = f.read(2) == b"\x1f\x8b"
    return gzip.open(path, "rt", encoding="latin-1") if gz else open(path, "rt", encoding="latin-1")


def _f(s):
    s = s.strip()
    return float(s) if s else None


def _i(s):
    s = s.strip()
    return int(s) if s.isdigit() else 0


def read_hipparcos(path):
    """hip_main.dat, bytes per the I/239 ReadMe (1-based, inclusive)."""
    rows = []
    with _open(path) as f:
        for line in f:
            if len(line) < 396 or line[0] != "H":
                continue
            v, ra, de = _f(line[41:46]), _f(line[51:63]), _f(line[64:76])
            if v is None or ra is None or de is None or v > MAG_LIMIT:
                continue
            pmra, pmde = _f(line[87:95]) or 0.0, _f(line[96:104]) or 0.0
            rows.append((_i(line[8:14]), ra, de, pmra, pmde, v, _i(line[390:396])))
    return rows


def read_tycho(path):
    """tyc_main.dat, same layout as hip_main.dat for the astrometry (I/239
    ReadMe): TYC 3-14, Vmag 42-46, RAdeg 52-63, DEdeg 65-76, pmRA 88-95,
    pmDE 97-104, HIP 211-216, HD 310-315. Only stars without a HIP number
    (the Hipparcos entry is used for those)."""
    rows = []
    with _open(path) as f:
        for line in f:
            if len(line) < 315 or line[0] != "T":
                continue
            if line[210:216].strip():
                continue
            v, ra, de = _f(line[41:46]), _f(line[51:63]), _f(line[64:76])
            if v is None or ra is None or de is None or v > MAG_LIMIT:
                continue
            pmra, pmde = _f(line[87:95]) or 0.0, _f(line[96:104]) or 0.0
            tyc = "-".join(line[2:14].split())
            rows.append((0, ra, de, pmra, pmde, v, _i(line[309:315]), tyc))
    return rows


def read_bsc(path):
    """catalog, bytes per the V/50 ReadMe: HR 1-4, Name 5-14, HD 26-31, SAO 32-37."""
    by_hd = {}
    with _open(path) as f:
        for line in f:
            if len(line) < 37:
                continue
            hr, hd, sao = _i(line[0:4]), _i(line[25:31]), _i(line[31:37])
            if hd:
                by_hd[hd] = (hr, line[4:14], sao)
    return by_hd


def read_sao(path):
    """sao.dat, bytes per the I/131A ReadMe: SAO 1-6, deleted flag 7, Vmag
    81-84, HD 118-123. One SAO number per HD (the brightest if several)."""
    by_hd = {}
    with _open(path) as f:
        for line in f:
            if len(line) < 123 or line[6] == "D":
                continue
            sao, hd = _i(line[0:6]), _i(line[117:123])
            if not (sao and hd):
                continue
            v = _f(line[80:84])
            v = 99.9 if v is None else v
            if hd not in by_hd or v < by_hd[hd][1]:
                by_hd[hd] = (sao, v)
    return {hd: s for hd, (s, _) in by_hd.items()}


def bayer_flamsteed(name10):
    """BSC Name field -> ('alp Leo' style key, display label) or ('', '')."""
    if not name10.strip():
        return "", ""
    flam, bay, sup, con = name10[0:3].strip(), name10[3:6].strip(), name10[6:7].strip(), name10[7:10].strip()
    if bay:
        g = GREEK.get(bay, bay)
        key = (bay.lower() + sup + " " + con).strip()
        return key, g + SUP.get(sup, sup) + " " + con
    if flam:
        return flam + " " + con, flam + " " + con
    return "", ""


def ecliptic_latitude(ra, de):
    ra, de = np.radians(ra), np.radians(de)
    s = np.sin(de) * np.cos(OBLIQUITY) - np.cos(de) * np.sin(OBLIQUITY) * np.sin(ra)
    return np.degrees(np.arcsin(np.clip(s, -1, 1)))


def build(hip_path, bsc_path, sao_path=None, tyc_path=None):
    rows = [r + ("",) for r in read_hipparcos(hip_path)]
    if tyc_path:
        rows += read_tycho(tyc_path)
    tycid = np.array([r[7] for r in rows], dtype="<U16")
    hip = np.array([r[0] for r in rows], np.int32)
    ra = np.array([r[1] for r in rows])
    de = np.array([r[2] for r in rows])
    pmra = np.array([r[3] for r in rows])
    pmde = np.array([r[4] for r in rows])
    mag = np.array([r[5] for r in rows], np.float32)
    hd = np.array([r[6] for r in rows], np.int32)
    # epoch J1991.25 -> J2000 along the proper motion (pmRA is mu_alpha*cos(delta))
    dt = 2000.0 - EPOCH_HIP
    de2 = de + pmde * dt / 3.6e6
    ra2 = (ra + pmra * dt / 3.6e6 / np.cos(np.radians(de))) % 360.0
    keep = np.abs(ecliptic_latitude(ra2, de2)) <= ECL_LAT_MAX
    bsc = read_bsc(bsc_path)
    sao_hd = read_sao(sao_path) if sao_path else {}
    n = int(keep.sum())
    tycid = tycid[keep]
    out = dict(ra=ra2[keep], dec=de2[keep], pmra=pmra[keep].astype(np.float32),
               pmdec=pmde[keep].astype(np.float32), mag=mag[keep], hip=hip[keep], hd=hd[keep],
               hr=np.zeros(n, np.int32), sao=np.zeros(n, np.int32))
    names, labels = [], []
    for k in range(n):
        h = int(out["hd"][k])
        hr, nm, sao = bsc.get(h, (0, "", 0)) if h else (0, "", 0)
        sao = sao or sao_hd.get(h, 0)
        out["hr"][k], out["sao"][k] = hr, sao
        key, lab = bayer_flamsteed(nm)
        names.append(key)
        labels.append(lab or ("HR %d" % hr if hr else "HD %d" % h if h else
                              "SAO %d" % sao if sao else "HIP %d" % out["hip"][k] if out["hip"][k] else
                              "TYC " + tycid[k]))
    out["name"] = np.array(names, dtype="<U40")
    out["label"] = np.array(labels, dtype="<U40")
    order = np.argsort(out["mag"], kind="stable")
    return {k: v[order] for k, v in out.items()}


def main(argv):
    args = list(argv[1:])
    tyc_path = None
    if "--tycho" in args:
        i = args.index("--tycho")
        tyc_path = args[i + 1]
        del args[i:i + 2]
    if len(args) not in (3, 4):
        print(__doc__)
        return 2
    hip_path, bsc_path = args[0], args[1]
    sao_path = args[2] if len(args) == 4 else None
    cat = build(hip_path, bsc_path, sao_path, tyc_path)
    np.savez_compressed(args[-1], **cat)
    print("%d stars to V %.1f within %.0f deg of the ecliptic (%d Hipparcos, %d Tycho); "
          "%d with HD, %d with HR, %d with SAO"
          % (len(cat["hip"]), MAG_LIMIT, ECL_LAT_MAX, (cat["hip"] > 0).sum(), (cat["hip"] == 0).sum(),
             (cat["hd"] > 0).sum(), (cat["hr"] > 0).sum(), (cat["sao"] > 0).sum()))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
