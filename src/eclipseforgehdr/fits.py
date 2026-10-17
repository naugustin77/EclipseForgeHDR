"""FITS input, for capture software that writes it rather than camera raw.

INDI/EKOS, SharpCap, N.I.N.A., FireCapture and most observatory control systems
save FITS. The format matters here because a FITS frame is already linear and
already black-subtracted by the driver -- so it arrives without the two things
`raw.py` normally gets from LibRaw, a black level and a saturation ceiling, and
both have to be recovered from the file.

READING IT
----------
`astropy.io.fits` if it is installed, then `fitsio`, then a built-in reader.
The built-in one handles plain uncompressed FITS, which is what every camera
writes; astropy is worth having for tile-compressed files and the stranger
corners of the standard, and the pipeline says so if it meets one it cannot
read. Nothing here is a hard dependency -- `pip install eclipseforgehdr[fits]`
adds astropy for anyone who wants the full coverage.

WHAT THE PIPELINE NEEDS FROM THE HEADER
---------------------------------------
    EXPTIME / EXPOSURE      the exposure in seconds. Without it there are no
                            tiers, and the run stops with that said plainly.
    BAYERPAT / BAYERPATTERN 'RGGB' / 'BGGR' / 'GRBG' / 'GBRG'. Absent, the
                            frame is treated as monochrome.
    XBAYROFF / YBAYROFF     odd Bayer origin offsets, applied before the roll.
    DATE-OBS                for the run report and the lunar track.
    GAIN / EGAIN / ISOSPEED reported only.
    PEDESTAL / BLKLEVEL /   a pedestal the driver added, subtracted if present.
    OFFSET / BZERO-style
    SATURATE / DATAMAX      the saturation ceiling, if the writer records one.

There is no colour matrix and no white balance in a FITS header, so both are
identity here and the frames stay in raw sensor colour through the merge. That
much has not changed, and for the reason it always had: inventing a white
balance before the photometry would put a colour error where it cannot be seen.

What changed is that leaving it there was worse. A silicon sensor under a CFA
is much more sensitive in green than in red, so raw-sensor white is not white:
on a 25-tier FITS bracket the corona measured R/((G+B)/2) = 0.56 against 1.11
for a colour-managed stack of the same corona, which is a blue-cyan rim on a
brown sky, and it also put the prominence gate's reference colour in the wrong
place. Warmth and Tint can chase that by eye, but they cannot know the answer.

So the balance is now measured, once, AFTER the merge and the photometry, from
the inner corona -- see `pipeline.neutralise_corona_colour`. The K-corona is
photospheric light Thomson-scattered off free electrons, and Thomson scattering
is wavelength-independent, so the inner corona carries the Sun's own spectrum:
it is white by physics, not by convention. That makes it a reference rather
than a guess, and taking it after the merge keeps the photometry clean.
"""
from __future__ import annotations
import os
import numpy as np

FITS_EXTS = {".fits", ".fit", ".fts", ".fts.fz", ".fits.fz"}

_BITPIX = {8: ">u1", 16: ">i2", 32: ">i4", 64: ">i8", -32: ">f4", -64: ">f8"}


def is_fits(path):
    n = os.path.basename(path).lower()
    return any(n.endswith(e) for e in FITS_EXTS)


# ---------- readers ----------

def _parse_header(f, path):
    """The primary header's cards from an open file positioned at 0; the file
    is left at the first data byte."""
    hdr = {}
    if f.read(6) != b"SIMPLE":
        raise ValueError(f"{os.path.basename(path)} is not a FITS file")
    f.seek(0)
    while True:
            block = f.read(2880)
            if not block or len(block) < 2880:
                raise ValueError(f"{os.path.basename(path)}: truncated FITS header")
            end = False
            for i in range(0, 2880, 80):
                card = block[i:i + 80].decode("latin-1")
                key = card[:8].strip()
                if key == "END":
                    end = True
                    break
                if not key or card[8:10] != "= ":
                    continue
                val = card[10:].split("/")[0].strip()
                if val.startswith("'"):
                    val = val.strip().strip("'").strip()
                elif val in ("T", "F"):
                    val = (val == "T")
                else:
                    try:
                        val = int(val)
                    except ValueError:
                        try:
                            val = float(val.replace("D", "E"))
                        except ValueError:
                            pass
                hdr.setdefault(key, val)
            if end:
                break
    return hdr


def _read_builtin_header(path):
    with open(path, "rb") as f:
        return _parse_header(f, path)


def _read_builtin(path):
    """Minimal reader for uncompressed FITS: 2880-byte blocks of 80-char cards,
    then big-endian data scaled by BZERO/BSCALE."""
    with open(path, "rb") as f:
        hdr = _parse_header(f, path)
        naxis = int(hdr.get("NAXIS", 0))
        if naxis < 2:
            raise ValueError(f"{os.path.basename(path)}: NAXIS={naxis}, not an image")
        dims = [int(hdr[f"NAXIS{i}"]) for i in range(1, naxis + 1)]
        bp = int(hdr.get("BITPIX", 16))
        if bp not in _BITPIX:
            raise ValueError(f"{os.path.basename(path)}: unsupported BITPIX {bp}")
        n = int(np.prod(dims))
        buf = f.read(n * abs(bp) // 8)
        if len(buf) < n * abs(bp) // 8:
            raise ValueError(f"{os.path.basename(path)}: truncated FITS data")
        data = np.frombuffer(buf, dtype=_BITPIX[bp]).reshape(dims[::-1])
    bz, bs = float(hdr.get("BZERO", 0.0)), float(hdr.get("BSCALE", 1.0))
    # SCALE IN FLOAT64 (audit 2026-10-10, A16): a 32-bit integer image with
    # BZERO 2^31 cast to float32 FIRST loses its low 8 bits (24-bit mantissa
    # against values near 2^31), and the pixels came out quantised in steps
    # of 256 ADU. float32 only after the offset is applied.
    if bs != 1.0 or bz != 0.0:
        data = data.astype(np.float64)
        if bs != 1.0:
            data *= bs
        if bz != 0.0:
            data += bz
    data = data.astype(np.float32)
    return data, hdr


def _orient(data, hdr):
    """Put a FITS image the same way up as every other image the pipeline sees.

    FITS row 1 is the BOTTOM of the frame (the standard's origin is lower-left);
    every raster format the rest of this code touches has row 0 at the top. Read
    straight through, a FITS bracket therefore comes out vertically mirrored --
    which is exactly what the first person to feed it FITS reported: "the image
    came out upside down".

    Astro capture software knows this is a trap and most of it writes ROWORDER.
    Honour that when it is there, and otherwise follow the standard.

    THE BAYER TRAP. Flipping an image with an EVEN number of rows swaps the
    CFA's row parity -- RGGB becomes GBRG -- because the new row 0 was the old
    row H-1, which is odd. Getting the flip right and the pattern wrong would
    trade an upside-down picture for a miscoloured one, so the parity shift is
    returned with the data and added to the Bayer y-offset below.

    Returns (data, what it did, the row-parity shift the flip introduced).
    """
    # Whatever this decides, it must be BAYER-correct -- that is the only part
    # of orientation that matters here, because a wrong CFA parity is a colour
    # error and cannot be undone downstream. Which way up the picture ends up is
    # settled at the other end of the pipeline, losslessly and without a
    # re-stack: see render.apply_orient.
    ro = str(hdr.get("ROWORDER", hdr.get("ROW_ORDER", "")) or "").strip().upper()
    if ro.startswith("TOP"):          # already top-down; leave it alone
        return data, "TOP-DOWN (ROWORDER)", 0
    why = "BOTTOM-UP (ROWORDER)" if ro.startswith("BOTTOM") else \
          "BOTTOM-UP (FITS default, no ROWORDER)"
    # FLIP THE ROW AXIS, NOT AXIS 0. A mono or CFA frame is (H, W) and axis 0
    # is rows, but a colour cube is stored PLANES FIRST -- (3, H, W) -- and
    # there axis 0 is the plane index, so `data[::-1]` reversed R, G, B into
    # B, G, R and left the picture upside down as well. Green sits in the
    # middle and is unmoved, which is what makes it hard to see: the corona
    # keeps a plausible colour and only the red/blue balance is inverted, and
    # the corona white balance downstream then scales it into looking almost
    # right. Found by planting 1000/2000/3000 in the three planes of a cube
    # with no ROWORDER card and reading back R=3000, G=2000, B=1000.
    #
    # Only files that reach this branch were ever affected: a cube carrying
    # ROWORDER='TOP-DOWN' returns above and never flipped. Siril writes that
    # card, so its exports are clean; the FITS default is bottom-up, so a
    # writer that omits ROWORDER was not.
    # the ROW axis of a 3-D image: planes-first (3, H, W) is the FITS
    # convention, but an (H, W, 3) array (NAXIS1 = 3) was flipped along its
    # columns here (audit 2026-10-10, A16)
    if data.ndim == 3:
        ax = 0 if (data.shape[-1] in (3, 4) and data.shape[0] not in (3, 4)) else -2
    else:
        ax = 0
    idx = [slice(None)] * data.ndim
    idx[ax] = slice(None, None, -1)
    return data[tuple(idx)].copy(), why, (data.shape[ax] - 1) % 2


# Where the parity shift rides from the reader to the Bayer code. Private, and
# in the header dict so it cannot be separated from the data it belongs to.
_YPAR = "_EFH_YPARITY"
_ROW = "_EFH_ROWORDER"


def read_fits(path):
    """(data, header dict), the data top-down. See _orient."""
    d, h = _read_fits_raw(path)
    d, why, par = _orient(d, h)
    h = dict(h)
    h[_YPAR], h[_ROW] = par, why
    return d, h


def row_order(path):
    """Which way up the file said it was; for the log (header only)."""
    try:
        hdr = _read_fits_header(path)
        ro = str(hdr.get("ROWORDER", hdr.get("ROW_ORDER", "")) or "").strip().upper()
        if ro.startswith("TOP"):
            return "TOP-DOWN (ROWORDER)"
        return "BOTTOM-UP (ROWORDER)" if ro.startswith("BOTTOM") else \
            "BOTTOM-UP (FITS default, no ROWORDER)"
    except Exception:
        return "unknown"


def _read_fits_raw(path):
    """(data, header dict), exactly as stored. astropy, fitsio, then built-in."""
    try:
        from astropy.io import fits as _af
        with _af.open(path, memmap=False) as hd:
            for h in hd:
                d = getattr(h, "data", None)
                if d is not None and getattr(d, "ndim", 0) >= 2:
                    return (np.asarray(d, np.float32),
                            {k: h.header[k] for k in h.header if k})
        raise ValueError(f"{os.path.basename(path)}: no image data in any HDU")
    except ImportError:
        pass
    try:
        import fitsio
        with fitsio.FITS(path) as ff:
            for h in ff:
                try:
                    d = h.read()
                except Exception:
                    continue
                if d is not None and getattr(d, "ndim", 0) >= 2:
                    hh = h.read_header()
                    return (np.asarray(d, np.float32),
                            {k: hh[k] for k in hh.keys()})
        raise ValueError(f"{os.path.basename(path)}: no image data in any HDU")
    except ImportError:
        pass
    try:
        return _read_builtin(path)
    except ValueError as e:
        raise ValueError(
            f"{e}. If this is a tile-compressed or otherwise unusual FITS, "
            f"install astropy (pip install 'eclipseforgehdr[fits]') and try again."
        ) from None


# ---------- metadata ----------

def _first(hdr, *keys):
    for k in keys:
        if k in hdr and hdr[k] not in ("", None):
            return hdr[k]
    return None


def _read_fits_header(path):
    """The first image HDU's header, without reading its pixels (audit
    2026-10-10, A16: every header read loaded the whole image -- minutes on
    a 250-frame bracket just to list exposure times)."""
    try:
        from astropy.io import fits as _af
        with _af.open(path, memmap=True, lazy_load_hdus=True) as hd:
            for h in hd:
                if int(h.header.get("NAXIS", 0)) >= 2:
                    return {k: h.header[k] for k in h.header if k}
        raise ValueError(f"{os.path.basename(path)}: no image data in any HDU")
    except ImportError:
        pass
    try:
        import fitsio
        with fitsio.FITS(path) as ff:
            for h in ff:
                hh = h.read_header()
                if int(hh.get("NAXIS", 0) or 0) >= 2:
                    return {k: hh[k] for k in hh.keys()}
        raise ValueError(f"{os.path.basename(path)}: no image data in any HDU")
    except ImportError:
        pass
    return _read_builtin_header(path)


def fits_exif(path):
    """(exposure_seconds, gain, timestamp) in the shape read_exif returns."""
    hdr = _read_fits_header(path)          # header only; the pixels stay on disk
    return _exif_from_header(hdr, path)


def _exif_from_header(hdr, path):
    sec = _first(hdr, "EXPTIME", "EXPOSURE", "EXPOSURE_TIME", "ITIME")
    if sec is None:
        raise ValueError(
            f"no EXPTIME or EXPOSURE in {os.path.basename(path)} — the exposure "
            f"time is what groups frames into tiers, so it cannot be guessed. "
            f"Most capture software writes it; if yours does not, the header can "
            f"be edited before processing.")
    gain = _first(hdr, "GAIN", "EGAIN", "ISOSPEED", "ISO") or 0
    ts = str(_first(hdr, "DATE-OBS", "DATE_OBS", "DATE") or "")
    try:
        gain = int(round(float(gain)))
    except (TypeError, ValueError):
        gain = 0
    return float(sec), gain, ts


def plane_clip_mask(rgb, min_px=64):
    """Pixels clipped in ANY plane of a 3-plane (H, W, 3) image, and the
    ceiling found per plane (None where the plane shows none).

    A writer that scales the planes (white balance, Siril's debayer) moves the
    clip level per plane: one FITS test set's Siril frames clip at 59889 in red, 65520
    in green and 65535 in blue, and a few hot pixels sit above the red one. So
    one ceiling for the frame (the bit depth, 65535) misses the red clip, and
    clipped red entered the merge as data (2026-10-13). The ceiling is found as
    a spike in each plane's histogram: the most common value in the upper
    half, held by at least `min_px` pixels and 50x more than a typical value
    there; everything from 0.5 % below it is flagged. A frame with no clipped
    core has no such spike and nothing is flagged."""
    clip = None
    ceil = []
    for c in range(rgb.shape[-1]):
        pl = rgb[..., c]
        m = float(pl.max()) if pl.size else 0.0
        if not np.isfinite(m) or m <= 0:
            ceil.append(None)
            continue
        top = pl[pl >= 0.5 * m]
        if top.size < min_px:
            ceil.append(None)
            continue
        if m < 2 ** 24 and np.array_equal(top[:1000], np.round(top[:1000])):
            lo = int(np.floor(0.5 * m))
            cnt = np.bincount(np.round(top).astype(np.int64) - lo)
            k = int(np.argmax(cnt))
            v, n = float(lo + k), int(cnt[k])
            typ = float(np.median(cnt[cnt > 0]))
        else:
            vals, cnt = np.unique(top, return_counts=True)
            k = int(np.argmax(cnt))
            v, n = float(vals[k]), int(cnt[k])
            typ = float(np.median(cnt))
        if n < min_px or n < 50.0 * typ:
            ceil.append(None)
            continue
        ceil.append(v)
        hit = pl >= 0.995 * v
        clip = hit if clip is None else (clip | hit)
    return clip, ceil


class FitsFrame:
    """A FITS bracket frame, presented to the pipeline as an RGGB mosaic.

    Colour cameras write a CFA frame plus BAYERPAT; that is rolled to RGGB the
    same way `raw.py` rolls a camera raw. Monochrome frames (no BAYERPAT) are
    laid into a mosaic with every photosite equal, so the rest of the pipeline
    runs unchanged and the result is a greyscale corona -- which is what a mono
    camera recorded. A 3-plane FITS is taken as already-debayered RGB.
    """

    def __init__(self, path):
        self.path = path
        self.plane_ceilings = None
        self.clipped_px = 0
        data, hdr = read_fits(path)
        self.header = hdr
        self.exposure, self.gain, self.timestamp = _exif_from_header(hdr, path)
        if data.ndim == 3:
            # FITS stores planes first
            if data.shape[0] in (3, 4):
                data = np.moveaxis(data[:3], 0, -1)
            elif data.shape[-1] in (3, 4):
                data = data[:, :, :3]
            else:
                raise ValueError(
                    f"{os.path.basename(path)}: {data.shape} is not an image "
                    f"this pipeline can use (expected 2-D CFA or 3-plane RGB)")
        data = np.asarray(data, np.float32)

        # a pedestal the driver added, if it says so.
        #
        # OFFSET IS NOT ON THIS LIST, and was until 0.23.3. In Siril, INDI and
        # SharpCap headers OFFSET is the camera's GAIN OFFSET SETTING -- a
        # configuration number, often 0 -- and not a pedestal in ADU. Siril says
        # so in the card's own comment: `OFFSET = 0 / Sensor gain offset`.
        # Reading it as a black level told the rest of the pipeline that the
        # black level was KNOWN and equal to zero, when in fact it was unknown
        # and, on a FITS test set, 142 ADU. Everything downstream is sized for a
        # few ADU of residue (_PEDESTAL_MAX caps the shared fit at 0.002 of
        # saturation, 131 ADU there), so a full black level cannot be recovered
        # by the fit that is supposed to mop up what this missed. See TODO 22.
        ped = _first(hdr, "PEDESTAL", "BLKLEVEL", "BLACKLEV")
        reported = ped is not None
        if ped is None and "OFFSET" in hdr:
            # RESCUE, so a writer that genuinely puts an ADU pedestal in OFFSET
            # is not silently dropped: take it only when the card's own comment
            # does NOT describe a camera setting. Siril writes
            # "Sensor gain offset"; INDI and SharpCap write similar. A writer
            # meaning a pedestal says "pedestal", "bias" or "black", or says
            # nothing at all -- and a bare OFFSET with no comment is still not
            # trusted, because the setting is the commoner meaning.
            try:
                _c = str(hdr.comments["OFFSET"]).lower()
            except Exception:
                _c = ""
            if any(w in _c for w in ("pedestal", "bias", "black")):
                ped = hdr["OFFSET"]
                reported = ped is not None
        try:
            ped = float(ped)
        except (TypeError, ValueError):
            ped = 0.0
            reported = False
        # PEDESTAL is conventionally the value ADDED, so subtract it; a negative
        # value in the header means the writer already recorded it as a
        # correction, so respect the sign it gave
        if ped:
            data = data - ped
        self.black_level = float(ped)
        # Whether the FILE said so, as opposed to us defaulting to zero. A raw
        # from LibRaw always reports one; a FITS usually does not. The
        # difference matters because "black level is zero" and "black level is
        # unknown" call for different behaviour downstream, and only the first
        # of those is safe to assume.
        self.black_level_reported = bool(reported)
        data = np.maximum(data, 0.0)

        if data.ndim == 3:
            H, W = data.shape[:2]
            H -= H % 2; W -= W % 2
            # CROP THE SOURCE TOO (0.22.76). Shrinking H/W to even sized the
            # DESTINATION but the slices below still ran over the full odd
            # height, so an odd-height frame raised "could not broadcast input
            # array from shape (51,100) into shape (50,100)" and killed the run.
            # The mono branch below has always done this; this one did not.
            data = data[:H, :W]
            _clip, self.plane_ceilings = plane_clip_mask(data)
            cfa = np.empty((H, W), np.float32)
            cfa[0::2, 0::2] = data[0::2, 0::2, 0]
            cfa[0::2, 1::2] = data[0::2, 1::2, 1]
            cfa[1::2, 0::2] = data[1::2, 0::2, 1]
            cfa[1::2, 1::2] = data[1::2, 1::2, 2]
            # a pixel clipped in ANY plane reads as the frame's ceiling at its
            # photosite, so the saturation test downstream (one level for the
            # whole bracket) cannot miss a plane that clips lower
            if _clip is not None:
                _bp = int(hdr.get("BITPIX", 16))
                _top = (65535.0 if _bp == 16 else 255.0 if _bp == 8
                        else float(np.nanmax(data)))
                cfa[_clip] = max(_top - self.black_level, float(cfa.max()), 1.0)
                self.clipped_px = int(_clip.sum())
            del _clip
            self.bayer = cfa
            self.mono = False
            self.pattern = "RGB (already debayered)"
        else:
            pat = _first(hdr, "BAYERPAT", "BAYERPATTERN", "COLORTYP")
            pat = str(pat).strip().upper() if pat else ""
            H, W = data.shape
            H -= H % 2; W -= W % 2
            d = data[:H, :W]
            if pat in ("RGGB", "BGGR", "GRBG", "GBRG"):
                # Bayer origin offsets, if the writer records them
                ox = int(_first(hdr, "XBAYROFF", "XBAYERPAT") or 0) % 2
                # + the parity the top-down flip introduced. Addition mod 2 is
                # correct whichever direction each shift ran in.
                oy = (int(_first(hdr, "YBAYROFF", "YBAYERPAT") or 0)
                      + int(hdr.get(_YPAR, 0))) % 2
                if ox or oy:
                    d = np.roll(np.roll(d, -oy, axis=0), -ox, axis=1)
                # roll so the pattern starts on R, exactly as raw.py does
                dy, dx = {"RGGB": (0, 0), "GRBG": (0, 1),
                          "GBRG": (1, 0), "BGGR": (1, 1)}[pat]
                if dy or dx:
                    d = np.roll(np.roll(d, -dy, axis=0), -dx, axis=1)
                self.bayer = np.ascontiguousarray(d, np.float32)
                self.mono = False
                self.pattern = pat
            else:
                # monochrome: same value at every photosite, so demosaic returns
                # a neutral image and nothing downstream has to special-case it
                self.bayer = np.ascontiguousarray(d, np.float32)
                self.mono = True
                self.pattern = "mono"

        # SATURATION CEILING. LibRaw hands raw.py a white level; a FITS header
        # usually does not. Taken in order: a keyword if the writer recorded one;
        # otherwise a saturation PLATEAU in the data -- a real ceiling shows up
        # as many pixels sharing the maximum, where a merely bright frame has a
        # handful; otherwise the bit-depth ceiling. Guessing high is the
        # dangerous direction (clipped pixels then merge as if valid), so the
        # plateau test is preferred over the nominal depth.
        sat = _first(hdr, "SATURATE", "DATAMAX", "SATLEVEL", "WHITELEV")
        try:
            sat = float(sat) - self.black_level
        except (TypeError, ValueError):
            sat = None
        mx = float(self.bayer.max()) if self.bayer.size else 1.0
        how = "SATURATE keyword"
        if sat is None or not np.isfinite(sat) or sat <= 0:
            plateau = int((self.bayer >= mx * 0.999).sum())
            # a real ceiling is a MINORITY of pixels sharing the maximum. If
            # most of the frame sits there it is a flat or uniform image, not a
            # saturated one, and taking its maximum as the ceiling would tell
            # the merge that a perfectly good frame is clipped everywhere.
            frac = plateau / max(self.bayer.size, 1)
            if plateau >= 50 and frac < 0.20:
                sat, how = mx, (f"saturation plateau ({plateau} px, "
                                f"{100 * frac:.2f}% of the frame)")
            else:
                # NOT min(nominal, mx) (0.22.76). With no plateau in the frame
                # that min() is just `mx`, the frame's own maximum -- and
                # pipeline.py freezes sat_level from the FIRST frame of the
                # SHORTEST tier, the one least likely to hold a plateau. A
                # 16-bit corona frame peaking at 12000 ADU then set an 11700 ADU
                # ceiling for the whole bracket, and every longer tier had its
                # inner corona declared clipped and dropped from the merge, with
                # nothing in the log to say so. Guessing high risks merging a
                # clipped pixel as valid; guessing low throws away good data in
                # every tier but one. The bit-depth ceiling is the honest guess
                # when the data shows no ceiling of its own; `mx` only raises it
                # if the data somehow exceeds it.
                bp = int(hdr.get("BITPIX", 16))
                nominal = (65535.0 if bp == 16 else
                           255.0 if bp == 8 else max(mx, 1.0))
                sat = max(nominal - self.black_level, max(mx, 1.0), 1.0)
                how = ("bit depth (%d-bit), no saturation plateau in the data "
                       "-- if this frame is clipped, add SATURATE to the header"
                       % bp)
        if data.ndim == 3:
            how += ("; clipped pixels found per plane and per frame, at each "
                    "plane's own ceiling")
        self.white_level = float(sat)
        self.sat_level = max(float(sat) * 0.975, 1.0)
        self.sat_source = how
        # no colour science in a FITS header -- see the module docstring
        self.daylight_wb = np.ones(3, np.float32)
        self.cam2rgb = np.eye(3, dtype=np.float32)

    @property
    def shape(self):
        return self.bayer.shape
