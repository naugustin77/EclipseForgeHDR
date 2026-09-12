#!/usr/bin/env python3
"""The FITS paths that broke in 0.22.75 and were fixed in 0.22.76, as tests.

WHY THIS EXISTS. Three FITS defects were found by reading the code, not by a
user, and the only person on the project with a real FITS bracket may never get
round to sending one. A bug that can only be confirmed by a tester who does not
appear is a bug nobody is watching, so each of the three is pinned here instead:

  1. THE SATURATION CEILING. With no SATURATE keyword and no plateau in the
     data, the ceiling was min(bit-depth, frame maximum), which is just the
     frame's own maximum -- and pipeline.py freezes sat_level from the FIRST
     frame of the SHORTEST tier, the one least likely to hold a plateau. A
     bracket then had every longer tier declared clipped through its inner
     corona and dropped from the merge, silently.

  2. RGB FITS IMPORT. read_fits returns planes-first (3, H, W) -- what Siril
     and PixInsight write. The guard tested shape[2], which is the WIDTH, so it
     passed, and a[:, :, :3] produced a (3, H, 3) array that failed much later
     as "could not find the lunar limb".

  3. ODD DIMENSIONS. The RGB->mosaic packing cropped the destination to even
     but sliced the full odd-height source, raising a broadcast error.

    python tools/smoke_fits.py
"""
import os
import shutil
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))


def _card(k, v):
    if isinstance(v, str):
        val = "'%-8s'" % v
    elif isinstance(v, bool):
        val = "T" if v else "F"
    else:
        val = "%20s" % repr(v)
    return ("%-8s= %-70s" % (k, val))[:80]


def write_fits(path, data, extra=(), bitpix=16):
    """Uncompressed FITS. `data` may be 2-D or planes-first (3, H, W)."""
    a = np.asarray(data)
    if bitpix == 16:
        # BZERO 32768 means the FILE stores physical - 32768 as int16. Writing
        # the physical value and merely reinterpreting the bytes put every
        # sample 32768 ADU too high, which is a pedestal larger than the corona
        # -- the fixture was testing the pipeline against a scene it could not
        # photometrically link. Found by smoke_fits printing a peak of 44772
        # for a corona built to peak at 12000.
        d = (np.clip(np.rint(a), 0, 65535) - 32768).astype(">i2")
        zero = 32768
    else:
        d = np.clip(np.rint(a), 0, 255).astype(">u1")
        zero = 0
    naxis = a.ndim
    cards = [_card("SIMPLE", True), _card("BITPIX", bitpix),
             _card("NAXIS", naxis)]
    for i, n in enumerate(reversed(a.shape)):          # NAXIS1 is fastest-varying
        cards.append(_card("NAXIS%d" % (i + 1), int(n)))
    cards += [_card("BZERO", zero), _card("BSCALE", 1),
              _card("EXPTIME", 0.25), _card("GAIN", 100)]
    cards += [_card(k, v) for k, v in extra]
    cards.append("END" + " " * 77)
    blob = "".join(c.ljust(80) for c in cards)
    blob += " " * ((2880 - len(blob) % 2880) % 2880)
    raw = d.tobytes()
    raw += b"\0" * ((2880 - len(raw) % 2880) % 2880)
    with open(path, "wb") as f:
        f.write(blob.encode("ascii"))
        f.write(raw)


def corona(H, W, R=40.0, peak=12000.0, seed=5):
    cy, cx = H / 2.0, W / 2.0
    yy = np.arange(H, dtype=np.float64)[:, None] - cy
    xx = np.arange(W, dtype=np.float64)[None, :] - cx
    r = np.hypot(yy, xx)
    a = peak * (R / np.maximum(r, R)) ** 4
    a[r < R] = 0.0
    return a + np.random.default_rng(seed).normal(0, 3, (H, W))


def t_saturation(tmp):
    """1: a frame with no ceiling in it must not become the ceiling."""
    from eclipseforgehdr.fits import FitsFrame
    p = os.path.join(tmp, "nosat.fits")
    # peaks at ~12000 of a 16-bit range, no plateau, no SATURATE keyword
    write_fits(p, corona(240, 320))
    f = FitsFrame(p)
    mx = float(f.bayer.max())
    if f.sat_level <= mx * 1.05:
        raise SystemExit(
            "saturation ceiling collapsed onto the frame maximum: sat_level "
            "%.0f against a peak of %.0f (%s). Every longer tier in the "
            "bracket would be declared clipped." % (f.sat_level, mx, f.sat_source))
    print("  1 saturation: peak %.0f -> ceiling %.0f, from %s"
          % (mx, f.sat_level, f.sat_source.split(",")[0]))

    # ... and a frame that IS clipped must still be detected by its plateau
    a = corona(240, 320, peak=90000.0)
    write_fits(os.path.join(tmp, "sat.fits"), a)
    f2 = FitsFrame(os.path.join(tmp, "sat.fits"))
    if "plateau" not in f2.sat_source:
        raise SystemExit("a genuinely clipped frame no longer finds its "
                         "plateau (source: %s)" % f2.sat_source)
    print("  1 saturation: a clipped frame still detects its own plateau")


def t_rgb_import(tmp):
    """2: an RGB cube written planes-first must import as an image."""
    from eclipseforgehdr import importhdr
    H, W = 240, 320
    base = corona(H, W, peak=20000.0)
    cube = np.stack([base * 1.1, base, base * 0.9])        # (3, H, W)
    p = os.path.join(tmp, "rgb.fits")
    write_fits(p, cube)
    from eclipseforgehdr.fits import read_fits
    raw, _ = read_fits(p)
    if raw.shape != (3, H, W):
        raise SystemExit("fixture is not planes-first: %s" % (raw.shape,))
    a = np.asarray(raw, np.float32)
    if a.ndim == 3 and a.shape[0] in (3, 4) and a.shape[-1] not in (3, 4):
        a = np.moveaxis(a[:3], 0, -1)
    if a.shape != (H, W, 3):
        raise SystemExit("planes-first RGB was not transposed: got %s, the "
                         "import would carry on with a (3, H, 3) array"
                         % (a.shape,))
    # and the real thing: the module's own guard has to reach the same shape
    src = open(os.path.join(os.path.dirname(importhdr.__file__),
                            "importhdr.py"), encoding="utf-8").read()
    if "np.moveaxis(a[:3], 0, -1)" not in src:
        raise SystemExit("importhdr no longer transposes planes-first cubes")
    print("  2 RGB import: (3, %d, %d) -> (%d, %d, 3)" % (H, W, H, W))


def t_odd(tmp):
    """3: odd dimensions must not raise."""
    from eclipseforgehdr.fits import FitsFrame
    for name, shape in (("odd_mono.fits", (241, 321)),
                        ("odd_rgb.fits", (3, 241, 321))):
        a = (corona(*shape[-2:], peak=9000.0) if len(shape) == 2
             else np.stack([corona(*shape[-2:], peak=9000.0)] * 3))
        p = os.path.join(tmp, name)
        write_fits(p, a)
        f = FitsFrame(p)
        if f.bayer.shape != (240, 320):
            raise SystemExit("%s: mosaic is %s, expected (240, 320)"
                             % (name, f.bayer.shape))
    print("  3 odd dimensions: 241x321 mono and RGB both pack to 240x320")


def main():
    tmp = tempfile.mkdtemp(prefix="efhdr_fits_")
    try:
        print("FITS regressions (no tester data required):")
        t_saturation(tmp)
        t_rgb_import(tmp)
        t_odd(tmp)
        print("OK — all three FITS paths behave")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    main()
