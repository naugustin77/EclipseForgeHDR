"""Minimal ICC v2 profile builder (sRGB primaries, choosable tone curve).

Why this exists: a 16-bit TIFF with no colour information is a guess for every
program that opens it.  Photoshop assigns its working space (sRGB), Affinity
does something similar, PixInsight assumes linear.  We write *scene-linear*
data, so a host that applies an sRGB decode to it is wrong by the whole
transfer function -- which is exactly the "overstretched, clipped histogram"
symptom you get from an untagged linear TIFF.

Rather than take a dependency on a colour-management library, this builds the
~1 KB profile by hand.  Both profiles use the sRGB/Rec.709 primaries and D65
white, Bradford-adapted to the D50 PCS, which is what the reference sRGB
profile does; only the TRC differs:

    srgb_profile()    the piecewise sRGB transfer function (1024-entry curve)
    linear_profile()  gamma 1.0

Reference: ICC.1:2001-04 (ICC v2), clauses 6.1 (header), 6.2 (tag table),
6.5.3 (curveType), 6.5.10 (textDescriptionType), 6.5.14 (XYZType).
"""
from __future__ import annotations
import struct
import numpy as np

# sRGB primaries and D65 white, Bradford-adapted to the D50 PCS.  These are the
# values carried by the reference sRGB profile, quoted to the s15Fixed16
# precision the format stores anyway.
_PRIMARIES_D50 = {
    "rXYZ": (0.43607, 0.22249, 0.01392),
    "gXYZ": (0.38515, 0.71687, 0.09708),
    "bXYZ": (0.14307, 0.06061, 0.71410),
}
_D50 = (0.96420, 1.00000, 0.82491)


def _s15f16(x):
    return struct.pack(">i", int(round(x * 65536.0)))


def _xyz_tag(xyz):
    return b"XYZ " + b"\0" * 4 + b"".join(_s15f16(v) for v in xyz)


def _text_tag(s):
    return b"text" + b"\0" * 4 + s.encode("ascii", "replace") + b"\0"


def _desc_tag(s):
    b = s.encode("ascii", "replace") + b"\0"
    return (b"desc" + b"\0" * 4 + struct.pack(">I", len(b)) + b
            + b"\0" * 8            # unicode language code + count
            + b"\0" * 3            # scriptcode code (u16) + count (u8)
            + b"\0" * 67)          # scriptcode description


def _curve_gamma(g):
    """curveType with a single u8Fixed8 gamma value."""
    return b"curv" + b"\0" * 4 + struct.pack(">I", 1) + struct.pack(">H", int(round(g * 256.0)))


def _curve_srgb(n=1024):
    """curveType as an n-entry LUT sampling the sRGB transfer function.

    The table maps *device* value -> linear light, i.e. it is the sRGB EOTF,
    which is the direction ICC TRC tags are defined in.
    """
    x = np.linspace(0.0, 1.0, n)
    y = np.where(x <= 0.04045, x / 12.92, ((x + 0.055) / 1.055) ** 2.4)
    t = np.clip(np.round(y * 65535.0), 0, 65535).astype(">u2")
    return b"curv" + b"\0" * 4 + struct.pack(">I", n) + t.tobytes()


def _pad4(b):
    return b + b"\0" * ((4 - len(b) % 4) % 4)


def _build(desc, trc):
    tags = [(b"desc", _desc_tag(desc)),
            (b"wtpt", _xyz_tag(_D50)),
            (b"rXYZ", _xyz_tag(_PRIMARIES_D50["rXYZ"])),
            (b"gXYZ", _xyz_tag(_PRIMARIES_D50["gXYZ"])),
            (b"bXYZ", _xyz_tag(_PRIMARIES_D50["bXYZ"])),
            (b"rTRC", trc), (b"gTRC", trc), (b"bTRC", trc),
            (b"cprt", _text_tag("Public domain"))]

    # rTRC/gTRC/bTRC are byte-identical, so point all three at one blob.
    blobs, offsets, data = {}, [], b""
    base = 128 + 4 + 12 * len(tags)
    for sig, blob in tags:
        if blob not in blobs:
            blobs[blob] = (base + len(data), len(blob))
            data += _pad4(blob)
        offsets.append((sig,) + blobs[blob])

    table = struct.pack(">I", len(tags)) + b"".join(
        sig + struct.pack(">II", off, size) for sig, off, size in offsets)

    size = 128 + len(table) + len(data)
    hdr = (struct.pack(">I", size) + b"\0" * 4 + struct.pack(">I", 0x02100000)
           + b"mntr" + b"RGB " + b"XYZ "
           + b"\0" * 12                       # date/time -- deliberately zero
           + b"acsp" + b"\0" * 4 * 4          # platform, flags, manufacturer, model
           + b"\0" * 8                        # attributes
           + struct.pack(">I", 0)             # rendering intent: perceptual
           + b"".join(_s15f16(v) for v in _D50)
           + b"\0" * 4 + b"\0" * 16 + b"\0" * 28)
    assert len(hdr) == 128, len(hdr)
    return hdr + table + data


def srgb_profile():
    """sRGB primaries, sRGB transfer function."""
    return _build("sRGB IEC61966-2.1 (EclipseForgeHDR)", _curve_srgb())


def linear_profile():
    """sRGB primaries, gamma 1.0 -- scene-linear."""
    return _build("Linear sRGB primaries (EclipseForgeHDR)", _curve_gamma(1.0))


def encode_srgb(a):
    """Scene-linear float in [0,1] -> sRGB-encoded float in [0,1]."""
    a = np.clip(a, 0.0, 1.0)
    return np.where(a <= 0.0031308, a * 12.92,
                    1.055 * np.power(np.maximum(a, 1e-12), 1.0 / 2.4) - 0.055)
