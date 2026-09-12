"""Raw decoding and metadata for exposure-bracketed eclipse sequences."""
from __future__ import annotations
import os
import numpy as np

RAW_EXTS = {".rw2", ".raf", ".nef", ".cr2", ".cr3", ".arw", ".orf", ".dng", ".pef"}


def list_raws(folder):
    """The frames of a bracket, in name order.

    TIFFs COUNT ONLY WHEN THERE IS NOTHING ELSE (0.22.78). A folder holding a
    bracket may also hold a finished HDR -- the server offers exactly such a
    file in the import box, and someone comparing an external stack against
    their own drops one in. It was then also counted as a frame, and Start
    either aborted at read_exif ("TIFF exports must keep EXIF") or, if the file
    had inherited EXIF from its source raw, merged a gamma-encoded finished
    image into that exposure tier. Anything Start writes goes to
    eclipseforge_output/, so this is only about files the user put there.

    A TIFF-only folder still works: that is the TIFF-bracket path, and it is
    reached whenever there is no raw or FITS file to prefer.
    """
    from .fits import is_fits
    raws, tiffs = [], []
    for name in sorted(os.listdir(folder)):
        ext = os.path.splitext(name)[1].lower()
        p = os.path.join(folder, name)
        if ext in RAW_EXTS or is_fits(name):
            raws.append(p)
        elif ext in {".tif", ".tiff"}:
            tiffs.append(p)
    return raws if raws else tiffs


def _exif_via_tifffile(path):
    """Fallback for TIFFs: read ExposureTime/ISO from the main or Exif IFD."""
    import tifffile
    with tifffile.TiffFile(path) as tf:
        tags = {t.name: t.value for t in tf.pages[0].tags.values()}
        exif = tags.get("ExifTag")
        if isinstance(exif, dict):
            tags.update(exif)
    exp = tags.get("ExposureTime")
    if exp is None:
        raise ValueError(f"no ExposureTime metadata in {path} — "
                         "TIFF exports must keep EXIF (e.g. Lightroom: include metadata)")
    sec = float(exp[0]) / float(exp[1]) if isinstance(exp, tuple) else float(exp)
    iso = tags.get("ISOSpeedRatings") or tags.get("PhotographicSensitivity") or 0
    if isinstance(iso, (tuple, list)):
        iso = iso[0]
    return sec, int(iso), str(tags.get("DateTime", ""))


def read_camera_info(path):
    """Best-effort camera/lens strings for the run report. Never raises."""
    info = {}
    try:
        import exifread
        with open(path, "rb") as f:
            tags = exifread.process_file(f, details=False)
        def g(*names):
            for n in names:
                if n in tags:
                    return str(tags[n]).strip()
            return None
        info["camera"] = " ".join(x for x in (g("Image Make"), g("Image Model")) if x) or None
        info["lens"] = g("EXIF LensModel", "MakerNote LensModel")
        fl = tags.get("EXIF FocalLength")
        if fl is not None:
            v = fl.values[0]
            info["focal_mm"] = float(v.num) / float(v.den)
        fn = tags.get("EXIF FNumber")
        if fn is not None:
            v = fn.values[0]
            info["f_number"] = float(v.num) / float(v.den)
    except Exception:
        pass
    return {k: v for k, v in info.items() if v}


def read_exif(path):
    """Return (exposure_seconds, iso, timestamp_str)."""
    from .fits import is_fits, fits_exif
    if is_fits(path):
        return fits_exif(path)
    try:
        import exifread
    except ImportError:
        return _exif_via_tifffile(path)
    with open(path, "rb") as f:
        tags = exifread.process_file(f, details=False)
    if not any(k in tags for k in ("EXIF ExposureTime", "Image ExposureTime")) \
            and os.path.splitext(path)[1].lower() in {".tif", ".tiff"}:
        return _exif_via_tifffile(path)
    exp = tags.get("EXIF ExposureTime") or tags.get("Image ExposureTime")
    if exp is None:
        raise ValueError(f"no ExposureTime metadata in {path} — "
                         "TIFF exports must keep EXIF (e.g. Lightroom: include metadata)")
    v = exp.values[0]
    sec = float(v.num) / float(v.den)
    iso_tag = tags.get("EXIF ISOSpeedRatings") or tags.get("Image ISOSpeedRatings")
    iso = int(iso_tag.values[0]) if iso_tag else 0
    ts = str(tags.get("EXIF DateTimeOriginal") or tags.get("Image DateTime") or "")
    return sec, iso, ts


def _rawpy_help(err):
    """Turn rawpy's import failure into something a user can act on.

    rawpy is a compiled extension, and when it will not load, Windows says
    only `DLL load failed while importing _rawpy: The specified module could
    not be found.` -- which names neither the missing thing nor the fix. Both
    known causes have simple answers, and neither is guessable from that line.
    """
    import platform
    msg = [f"the raw decoder (rawpy) could not be loaded: {err}", ""]
    if os.name == "nt":
        if platform.machine().upper() in ("ARM64", "AARCH64"):
            msg += [
                "This Python is the ARM64 build of Windows Python, and rawpy",
                "publishes no ARM64 wheel -- so it can never load here. Install",
                "the x86-64 build of Python from python.org instead (Windows",
                "on ARM runs it under emulation), then reinstall this app.",
            ]
        else:
            msg += [
                "On Windows this is nearly always the missing Microsoft Visual",
                "C++ Redistributable, which rawpy's compiled part needs and a",
                "fresh Windows does not have. Install it from",
                "  https://aka.ms/vs/17/release/vc_redist.x64.exe",
                "reboot, and try again.",
            ]
    else:
        msg += ["Reinstalling the app usually rebuilds it: pipx install --force ."]
    msg += ["", "FITS input does not use rawpy and still works."]
    return "\n".join(msg)


class RawFile:
    """Decoded bayer data normalized to RGGB layout, plus color info."""

    def __init__(self, path):
        try:
            import rawpy
        except ImportError as e:
            raise ImportError(_rawpy_help(e)) from None
        self.path = path
        with rawpy.imread(path) as raw:
            bayer = raw.raw_image_visible.astype(np.float32)
            # Per-channel black, not the mean of the four. Where a body really
            # does report different pedestals per colour, the leftover offset is
            # then multiplied by the white-balance gain (~2x in red) and lands on
            # the outer corona, which sits only tens of ADU above black -- an
            # exposure-dependent colour cast exactly where the signal is weakest.
            # On bodies that report four equal values this is a no-op.
            blp = np.asarray(raw.black_level_per_channel, np.float32)[:4]
            black = float(np.mean(blp))
            self.white_level = float(raw.white_level)
            # Read INSIDE the `with`: rawpy closes the file on exit and every
            # attribute then raises. Reading it after the block (0.22.34) meant
            # the guard below always fell through to "no linearity limit given"
            # on a body that reports one perfectly well.
            try:
                _lm_raw = raw.camera_white_level_per_channel
            except Exception:
                _lm_raw = None
            pat0 = raw.raw_pattern
            if (pat0 is not None and pat0.shape == (2, 2)
                    and float(blp.max() - blp.min()) > 0.5):
                self.bayer = bayer.copy()
                for yy in range(2):
                    for xx in range(2):
                        self.bayer[yy::2, xx::2] -= float(blp[int(pat0[yy, xx])])
            else:
                self.bayer = bayer - black
            # roll pattern to RGGB
            pat = raw.raw_pattern
            if pat is None or pat.shape != (2, 2) or set(pat.ravel()) != {0, 1, 2, 3}:
                raise RuntimeError(
                    "this sensor is not a 2x2 Bayer mosaic (X-Trans, Foveon and "
                    "monochrome sensors are not supported); its colour would be "
                    "decoded wrongly rather than failing visibly")
            pat = pat.copy()  # 2x2 of color indices, 0=R 1=G 2=B 3=G2
            dy = dx = 0
            for yy in range(2):
                for xx in range(2):
                    if pat[yy, xx] == 0:
                        dy, dx = yy, xx
            if dy or dx:
                self.bayer = np.roll(np.roll(self.bayer, -dy, axis=0), -dx, axis=1)
            self.daylight_wb = _norm_wb(raw.daylight_whitebalance)
            # AS SHOT, the multipliers the photographer's own setting produced.
            # Kept alongside daylight rather than instead of it; which one the
            # merge uses is a setting (see pick_wb).
            self.camera_wb = _norm_wb(raw.camera_whitebalance)
            cm = np.asarray(raw.rgb_xyz_matrix, np.float32)[:3, :3]  # XYZ -> cam
            self.cam2rgb = _cam2srgb(cm)
        # THE CONTAINER CEILING IS NOT THE SATURATION LEVEL.
        #
        # LibRaw reports two different ceilings and they are not the same thing.
        # `white_level` is the container's -- 16383 for 14-bit Canon data.
        # `camera_white_level_per_channel` is LibRaw's linear_max: the level at
        # which the sensor stops responding linearly, which is where the data
        # actually stops. On the test set's Canon 1500D they differ by 8.6%:
        #
        #     _MG_4637.CR2   white_level 16383   linear_max 15092   black 2047
        #
        # Taking the container ceiling put sat_level at 13977 when the data
        # clips at 13045, so every photosite on the flat top of a blown
        # highlight was declared VALID and merged at full weight. That is a
        # contour-shaped error at each tier's true clip radius, and it is what
        # printed the rings: 230% tier-to-tier disagreement at the limb, a 74 px
        # disagreement rim, and concentric arcs sitting within 0.6 px of a clip
        # contour. Adobe writes the true value into a DNG (white_level 15092
        # there), which is why the same code gave a clean run on DNGs of the
        # same nine frames -- the reference set's test, 0.22.33, and the thing that found it.
        #
        # Guarded, because linear_max is not always populated: use it only when
        # every channel reports a positive value below the container ceiling and
        # not implausibly far below it (a stray small value would collapse the
        # range and make the whole frame read as clipped).
        self.linear_max = None
        # Which of the three sources the saturation level ends up coming from.
        # Reported in the log, because "the reported white level" was printed on
        # a run whose number actually came from the observed maximum -- true and
        # misleading at the same time, on exactly the line that exists to make
        # this diagnosable. the reference set's 600 mm S1R II run, 0.22.36.
        self.sat_source = "white_level"
        _lm = _lm_raw
        if _lm is not None:
            _lm = np.asarray([v for v in _lm if v is not None], np.float64)
            if (_lm.size and np.isfinite(_lm).all() and _lm.min() > 0
                    and _lm.min() < self.white_level
                    and _lm.min() > 0.5 * self.white_level
                    and _lm.min() > black + 1.0):
                self.linear_max = float(_lm.min())
                self.white_level = self.linear_max
                self.sat_source = "linear_max"
        # A fraction of the usable range, not a literal 400 ADU margin. 400 is
        # 2.8% of a 14-bit range but 11% of a 12-bit one, and on a camera whose
        # reported white level is already black-subtracted it can go negative --
        # which makes every pixel read as saturated and the merge come out black.
        self.sat_level = max((self.white_level - black) * 0.975, 1.0)
        # Cross-check the reported white level against the data. LibRaw
        # occasionally reports 0, an already-black-subtracted value, or a 12-bit
        # ceiling for 14-bit data. Any of those makes sat_level far too low, so
        # the merge treats the whole inner corona as clipped in EVERY tier,
        # weights go to zero and the result is a black frame -- cached, with no
        # exception raised. The observed maximum is a hard lower bound on where
        # saturation can be: nothing can read above it.
        #
        # AND IT IS A MUCH WEAKER SAFETY NET THAN IT LOOKS, because of WHICH
        # frame reaches it. pipeline.py takes color_info from the FIRST frame it
        # opens, and it walks the tiers shortest-first -- so this cross-check is
        # always applied to the shortest exposure in the bracket, which is the
        # one least likely to contain a saturated photosite at all.
        #
        # Measured on the reference S1R II bracket, on the 1/4000 s frame that
        # sets color_info for the whole run: raw maximum 1766 ADU. Nothing
        # within four stops of saturation, so this branch cannot fire on that
        # set. It is the mechanism that is unguarded, not this dataset.
        #
        # THE BLACK LEVEL ON THAT BODY IS 512, AND THE CEILING IS 16319.
        # An earlier version of this comment said 576 and 16383, on the theory
        # that LibRaw adds its own black of 64 to the cblack of 512 that
        # raw-identify prints. That was wrong. Both readings give the same
        # usable range -- 16383 - 576 and 16319 - 512 are both 15807, which is
        # why sat_level is 15411.8 either way and why the report cannot tell
        # them apart. The data can:
        #
        #   - decoded with LibRaw's unprocessed_raw, which modifies nothing, the
        #     sky in five widely separated patches reads a median of exactly
        #     512.0 on the 1/4000 s frame, 512-513 at 1/500 s, 516-517 at
        #     1/100 s and 526-532 at 1/30 s. A floor of 512 with real sky
        #     accumulating linearly on top of it: 0.113 ADU per 1/4000 s unit,
        #     consistent across the last two tiers.
        #   - with 576 subtracted the sky would be NEGATIVE on every tier,
        #     -64 ADU at the short end and -49 ADU at 1/30 s.
        #   - the shared pedestal settles it. Refitting it here on those frames
        #     gives -0.16 ADU with tier disagreement 45.1% -> 13.2% for a black
        #     of 512, against a value PEGGED at the +-0.002*sat clamp and
        #     disagreement 582% -> 571% for 576. The real 0.22.38 run on that
        #     dataset fitted -0.46 ADU with 44.6% -> 8.4%. Only 512 reproduces
        #     the run.
        #
        # So rawpy is reporting white_level 16319 -- the linearity limit, as its
        # own field rather than through camera_white_level_per_channel -- and
        # black 512, and the pipeline is subtracting the right number. Read the
        # frame, not the metadata summary, and check the fit against a real run
        # before rewriting this again.
        _obs = float(self.bayer.max()) if self.bayer.size else 0.0
        if _obs > self.sat_level:
            self.sat_level = _obs * 0.975
            self.white_level = _obs + black
            self.sat_source = "observed"

    @property
    def shape(self):
        return self.bayer.shape


def _norm_wb(v):
    """LibRaw multipliers -> a gain triple normalised to G = 1.

    Some bodies report zeros here; without the guard the whole merge becomes
    inf/nan and gets cached that way.
    """
    w = np.asarray(v, np.float32).ravel()[:3]
    if w.size < 3 or not np.isfinite(w).all() or w[1] <= 0 or (w <= 0).any():
        return np.ones(3, np.float32)
    return (w / w[1]).astype(np.float32)


WB_SOURCES = ("camera", "daylight", "none")


def pick_wb(rf, source="camera"):
    """Which white balance the merge should apply, and what to call it.

    'camera' is AS SHOT -- the multipliers the photographer's own setting
    produced -- and is the default since 0.22.63.

    Daylight was the default before that, on the argument that it is a property
    of the camera model rather than of a menu setting and so is reproducible
    between photographers. That argument is weaker than it looked: for a body
    LibRaw has no table entry for, `daylight_whitebalance` is DERIVED from the
    colour matrix, not measured. On the Panasonic DC-S1RM2 it comes out at
    R 2.258 / B 1.232 against the camera's own Fine Weather preset of
    R 2.227 / B 1.707 -- a 39% difference in blue, all of it in the direction
    that renders the corona too orange.

    Falls back to daylight when a body reports no usable as-shot value, and to
    unity for readers that carry no white balance at all (FITS, TIFF).
    """
    if source == "none":
        return np.ones(3, np.float32), "none (1:1:1)"
    cam = getattr(rf, "camera_wb", None)
    day = getattr(rf, "daylight_wb", np.ones(3, np.float32))
    if source == "camera":
        if cam is not None and not np.allclose(cam, 1.0, atol=1e-6):
            return np.asarray(cam, np.float32), "camera (as shot)"
        if np.allclose(day, 1.0, atol=1e-6):
            # no white balance in this file type at all
            return np.ones(3, np.float32), "none (this format carries no white balance)"
        return np.asarray(day, np.float32), \
            "daylight (this file reports no usable as-shot value)"
    return np.asarray(day, np.float32), "daylight"


def _cam2srgb(xyz2cam):
    """Standard DNG-style matrix: normalize cam_from_sRGB rows, invert."""
    srgb2xyz = np.array([[0.4124564, 0.3575761, 0.1804375],
                         [0.2126729, 0.7151522, 0.0721750],
                         [0.0193339, 0.1191920, 0.9503041]], np.float32)
    if not np.isfinite(xyz2cam).all() or abs(xyz2cam).sum() < 1e-6:
        return np.eye(3, dtype=np.float32)
    cam_from_srgb = xyz2cam @ srgb2xyz
    cam_from_srgb /= cam_from_srgb.sum(axis=1, keepdims=True)
    return np.linalg.inv(cam_from_srgb).astype(np.float32)


# Malvar-He-Cutler demosaic (RGGB)
_KG = np.array([[0, 0, -1, 0, 0], [0, 0, 2, 0, 0], [-1, 2, 4, 2, -1],
                [0, 0, 2, 0, 0], [0, 0, -1, 0, 0]], np.float32) / 8
_KR = np.array([[0, 0, 0.5, 0, 0], [0, -1, 0, -1, 0], [-1, 4, 5, 4, -1],
                [0, -1, 0, -1, 0], [0, 0, 0.5, 0, 0]], np.float32) / 8
_KC = _KR.T.copy()
_KD = np.array([[0, 0, -1.5, 0, 0], [0, 2, 0, 2, 0], [-1.5, 0, 6, 0, -1.5],
                [0, 2, 0, 2, 0], [0, 0, -1.5, 0, 0]], np.float32) / 8


def cfa_clip_max(cfa):
    """The largest PHOTOSITE value that feeds each demosaiced pixel.

    THIS IS THE CLIPPING TEST, AND IT HAS TO BE ASKED OF THE MOSAIC.

    Until 0.22.26 the merge asked it of the demosaiced result --
    `demosaic_rggb(cfa).max(axis=2)` -- which is a different question and
    silently the wrong one. Two of a demosaiced pixel's three channels are
    interpolated from photosites up to two away through the kernels above, and
    those kernels have negative lobes: an interpolated channel next to a
    saturated photosite can come out BELOW the clipping threshold. The pixel is
    then declared valid and enters the merge at full weight carrying a value
    reconstructed, in part, from a photosite that hit the ceiling.

    Measured on the 600 mm reference set: adjacent tiers disagree by up to 3x within
    0-8 px outside the longer tier's saturated region and by ~1% beyond
    8-16 px. This is one of the two mechanisms that can produce that collar --
    the other, charge spill or veiling glare off the saturated area, is real
    physics and is not fixed here.

    The footprint is the union of the four kernels' non-zero taps, not a full
    5x5 square: the corners are zero in every one of them, so a corner
    photosite genuinely does not reach the centre pixel and excluding it would
    be needlessly conservative.

    A single saturated photosite in an otherwise dim field, ceiling 16000 ADU,
    background 300:

        old test (demosaiced max)      new test (mosaic max)
          0.3   0.3   0.3               0.3  16.0   0.3
          0.3   4.2   8.2   4.2         16.0 16.0  16.0
          0.3   8.2  16.0   8.2         16.0 16.0  16.0     (x1000 ADU)

    It contaminates 13 pixels and the old test flags 1. The nearest neighbours
    read 8200 and 4200 against a 16000 ceiling -- half-built from a clipped
    photosite and comfortably under any threshold -- and the diamond tips at
    distance 2 read exactly the background, reaching the centre only through a
    negative lobe, undetectable by any test applied after the demosaic.

    Honesty about what this buys: hard-excluding a collar of 2-16 px around
    each tier's saturated region was measured on the bench (tools/collar.py)
    and moved the ring artifact by 1-2%. This is a correctness fix -- clipped
    photosites no longer enter the merge unflagged -- not a cure for that.
    """
    from scipy import ndimage
    fp = ((_KG != 0) | (_KR != 0) | (_KC != 0) | (_KD != 0))
    fp[2, 2] = True                     # the pixel's own photosite
    return ndimage.maximum_filter(cfa, footprint=fp, mode="nearest")




# ---------------------------------------------------------------------------
# VNG: "Interpolation using a Threshold-based Variable Number of Gradients"
# (Chang, Cheng & Ward 1999). This is the algorithm PixInsight's RAW module
# uses, so it is what someone comparing the two renders is comparing against.
#
# The two tables below were EXTRACTED from LibRaw's own misc_demosaic.cpp
# rather than retyped, so they cannot have drifted in transcription. The
# reading of them follows LibRaw exactly:
#
#   1. a weighted bilinear fill first (orthogonal 2, diagonal 1) -- the
#      gradients are measured on THAT, not on the mosaic;
#   2. eight directional gradients per pixel, numbered clockwise from NW;
#   3. threshold at gmin + gmax/2, keep the directions at or below it;
#   4. the measured channel is passed through untouched and the other two are
#      reconstructed from the COLOUR DIFFERENCE averaged over those directions.
#
# Step 4 is why VNG holds an edge better than a plain bilinear and why it can
# also smear a one-pixel feature: the difference is averaged over however many
# directions survived the threshold, which in flat noise is all eight.
#
# NOT THE DEFAULT. Malvar-He-Cutler stays the default -- it has lower
# interpolation error on the same data, and it is what every render up to
# 0.22.62 was made with. This exists so a render can be compared with
# PixInsight's on equal terms.
_VNG_CHOOD = ((-1, -1), (-1, 0), (-1, 1), (0, 1),
              (1, 1), (1, 0), (1, -1), (0, -1))

_VNG_TERMS = (
    (-2,-2,0,-1,0,1), (-2,-2,0,0,1,1), (-2,-1,-1,0,0,1),
    (-2,-1,0,-1,0,2), (-2,-1,0,0,0,3), (-2,-1,0,1,1,1),
    (-2,0,0,-1,0,6), (-2,0,0,0,1,2), (-2,0,0,1,0,3),
    (-2,1,-1,0,0,4), (-2,1,0,-1,1,4), (-2,1,0,0,0,6),
    (-2,1,0,1,0,2), (-2,2,0,0,1,4), (-2,2,0,1,0,4),
    (-1,-2,-1,0,0,-128), (-1,-2,0,-1,0,1), (-1,-2,1,-1,0,1),
    (-1,-2,1,0,1,1), (-1,-1,-1,1,0,-120), (-1,-1,1,-2,0,64),
    (-1,-1,1,-1,0,34), (-1,-1,1,0,0,51), (-1,-1,1,1,1,17),
    (-1,0,-1,2,0,8), (-1,0,0,-1,0,68), (-1,0,0,1,0,17),
    (-1,0,1,-2,1,64), (-1,0,1,-1,0,102), (-1,0,1,0,1,34),
    (-1,0,1,1,0,51), (-1,0,1,2,1,16), (-1,1,1,-1,1,68),
    (-1,1,1,0,0,102), (-1,1,1,1,0,34), (-1,1,1,2,0,16),
    (-1,2,0,1,0,4), (-1,2,1,0,1,4), (-1,2,1,1,0,4),
    (0,-2,0,0,1,-128), (0,-1,0,1,1,-120), (0,-1,1,-2,0,64),
    (0,-1,1,0,0,17), (0,-1,2,-2,0,64), (0,-1,2,-1,0,32),
    (0,-1,2,0,0,48), (0,-1,2,1,1,16), (0,0,0,2,1,8),
    (0,0,2,-2,1,64), (0,0,2,-1,0,96), (0,0,2,0,1,32),
    (0,0,2,1,0,48), (0,0,2,2,1,16), (0,1,1,0,0,68),
    (0,1,1,2,0,16), (0,1,2,-1,1,64), (0,1,2,0,0,96),
    (0,1,2,1,0,32), (0,1,2,2,0,16), (1,-2,1,0,0,-128),
    (1,-1,1,1,0,-120), (1,0,1,2,0,8), (1,0,2,-1,0,64),
    (1,0,2,1,0,16),
)


def demosaic_rggb(cfa):
    from scipy import ndimage
    conv = lambda k: ndimage.convolve(cfa, k, mode="mirror")
    gm, rrow, rcol, rdiag = conv(_KG), conv(_KR), conv(_KC), conv(_KD)
    R = np.empty_like(cfa); G = np.empty_like(cfa); B = np.empty_like(cfa)
    G[:] = gm
    G[0::2, 1::2] = cfa[0::2, 1::2]
    G[1::2, 0::2] = cfa[1::2, 0::2]
    R[0::2, 0::2] = cfa[0::2, 0::2]
    R[0::2, 1::2] = rrow[0::2, 1::2]
    R[1::2, 0::2] = rcol[1::2, 0::2]
    R[1::2, 1::2] = rdiag[1::2, 1::2]
    B[1::2, 1::2] = cfa[1::2, 1::2]
    B[1::2, 0::2] = rrow[1::2, 0::2]
    B[0::2, 1::2] = rcol[0::2, 1::2]
    B[0::2, 0::2] = rdiag[0::2, 0::2]
    return np.stack([R, G, B], axis=-1)


def _fcol(y, x):
    """RGGB: 0=R 1=G 2=B. LibRaw's fcol() for this pattern. The +144 offsets
    below are LibRaw's own trick for keeping negative indices positive; they
    are multiples of the 2x2 period, so they change nothing."""
    return ((0, 1), (1, 2))[y & 1][x & 1]


def _lin_interp(cfa):
    """LibRaw's lin_interpolate(): the weighted bilinear fill VNG measures its
    gradients on. Orthogonal neighbours weight 2, diagonal 1, and only
    neighbours of a colour OTHER than the centre's contribute."""
    H, W = cfa.shape
    P = np.pad(cfa, 1, mode="edge")
    out = np.empty((H, W, 3), np.float32)
    for py in (0, 1):
        for px in (0, 1):
            f = _fcol(py, px)
            nr, nc = len(range(py, H, 2)), len(range(px, W, 2))
            acc, wsum = [None, None, None], [0, 0, 0]
            for dy in (-1, 0, 1):
                for dx in (-1, 0, 1):
                    c = _fcol(py + dy + 48, px + dx + 48)
                    if c == f:
                        continue
                    w = 1 << ((dy == 0) + (dx == 0))
                    v = P[py + 1 + dy: py + 1 + dy + 2 * nr: 2,
                          px + 1 + dx: px + 1 + dx + 2 * nc: 2]
                    acc[c] = v * w if acc[c] is None else acc[c] + v * w
                    wsum[c] += w
            for c in range(3):
                sl = (slice(py, py + 2 * nr, 2), slice(px, px + 2 * nc, 2), c)
                out[sl] = cfa[py::2, px::2] if c == f else acc[c] / wsum[c]
    return out


def demosaic_vng(cfa, band=1024):
    """VNG demosaic, RGGB. Same result shape and dtype as demosaic_rggb.

    Banded over rows so peak memory stays near the Malvar path's: the eight
    gradient planes are a quarter-frame each, and holding them for a whole
    8000x5400 frame at once costs about 350 MB on top of everything else.
    """
    H, W = cfa.shape
    img = _lin_interp(np.asarray(cfa, np.float32))
    out = np.empty_like(img)
    for r0 in range(0, H, band):
        _vng_band(img, out, r0, min(H, r0 + band), H, W)
    return out


def _vng_band(img, out, r0, r1, H, W):
    y0, y1 = max(0, r0 - 2), min(H, r1 + 2)
    P = np.pad(img[y0:y1], ((2 - (r0 - y0), 2 - (y1 - r1)), (2, 2), (0, 0)),
               mode="edge")
    for py in (0, 1):
        rr0 = r0 + ((py - r0) % 2)          # first absolute row of this phase
        if rr0 >= r1:
            continue
        nr = len(range(rr0, r1, 2))
        oy = (rr0 - r0) + 2                 # its row inside P
        for px in (0, 1):
            nc = len(range(px, W, 2))
            ox = px + 2

            def nb(dy, dx, c):
                return P[oy + dy: oy + dy + 2 * nr: 2,
                         ox + dx: ox + dx + 2 * nc: 2, c]

            f = _fcol(rr0, px)
            gval = np.zeros((8, nr, nc), np.float32)
            for (Y1, X1, Y2, X2, wt, gr) in _VNG_TERMS:
                color = _fcol(rr0 + Y1 + 144, px + X1 + 144)
                if _fcol(rr0 + Y2 + 144, px + X2 + 144) != color:
                    continue
                # `diag` follows the TERM's colour, not the centre's -- green
                # sits on both diagonals of a Bayer cell and red and blue do not
                diag = 2 if (_fcol(rr0 + 144, px + 145) == color
                             and _fcol(rr0 + 145, px + 144) == color) else 1
                if abs(Y1 - Y2) == diag and abs(X1 - X2) == diag:
                    continue
                d = np.abs(nb(Y1, X1, color) - nb(Y2, X2, color))
                if wt:
                    d = d * float(1 << wt)
                m = gr & 0xFF               # the table stores this as a signed char
                for g in range(8):
                    if m & (1 << g):
                        gval[g] += d
            # LibRaw's threshold is gmin + (gmax >> 1) on integers; on float
            # data the shift is a halving and that is what this is
            thold = gval.min(axis=0) + 0.5 * gval.max(axis=0)
            good = gval <= thold[None]
            num = good.sum(axis=0).astype(np.float32)
            centre = nb(0, 0, f)
            ssum = np.zeros((3, nr, nc), np.float32)
            for g, (dy, dx) in enumerate(_VNG_CHOOD):
                # where the 1-step neighbour is a different colour but the
                # 2-step one is the centre's, LibRaw averages the centre with
                # that 2-step pixel instead of reading the neighbour
                two = (_fcol(rr0 + dy + 144, px + dx + 144) != f
                       and _fcol(rr0 + 2 * dy + 144, px + 2 * dx + 144) == f)
                gm = good[g]
                for c in range(3):
                    v = (0.5 * (centre + nb(2 * dy, 2 * dx, f))
                         if (c == f and two) else nb(dy, dx, c))
                    ssum[c] += np.where(gm, v, 0.0)
            for c in range(3):
                sl = (slice(rr0, r1, 2), slice(px, W, 2), c)
                out[sl] = centre if c == f else \
                    centre + (ssum[c] - ssum[f]) / num


DEMOSAICS = ("mhc", "vng")


def demosaic(cfa, method="mhc"):
    """The one entry point. `method` is a SETTING; see DEMOSAICS."""
    return demosaic_vng(cfa) if method == "vng" else demosaic_rggb(cfa)


TIFF_EXTS = {".tif", ".tiff"}


class TiffFrame:
    """16-bit (or 8-bit) TIFF bracket frame. Assumes display-gamma (sRGB) encoding
    unless ECLIPSEFORGE_TIFF_LINEAR=1 is set or the filename contains 'linear'.
    The RGB data is converted to a synthetic RGGB mosaic so the rest of the
    pipeline (stacking, saturation weighting, demosaic) runs unchanged."""

    SCALE = 15000.0

    def __init__(self, path):
        import tifffile
        self.path = path
        arr = tifffile.imread(path)
        if arr.ndim == 2:
            arr = np.stack([arr] * 3, axis=-1)
        if arr.shape[2] > 3:
            arr = arr[:, :, :3]
        # SCALE BY THE DTYPE, NOT BY THE DATA (0.22.76). This used to divide by
        # 65535 or 255 depending on whether the brightest pixel exceeded 255.5 --
        # so a 16-bit frame whose brightest pixel was below 256, which is exactly
        # what a very short totality frame looks like, was divided by 255 and
        # arrived 257x too bright, then clipped to 1.0: the whole tier read as
        # saturated. importhdr.py documents fixing the same bug; this copy was
        # missed. Floats are left alone (already 0..1 by convention) and only
        # clipped.
        _int = np.issubdtype(arr.dtype, np.integer)
        _full = float(np.iinfo(arr.dtype).max) if _int else 1.0
        arr = arr.astype(np.float32)
        if _int:
            arr /= _full
        arr = np.clip(arr, 0, 1)
        linear = (os.environ.get("ECLIPSEFORGE_TIFF_LINEAR") == "1"
                  or "linear" in os.path.basename(path).lower())
        if not linear:  # inverse sRGB
            a = arr
            arr = np.where(a <= 0.04045, a / 12.92, ((a + 0.055) / 1.055) ** 2.4)
        arr *= self.SCALE
        H, W = arr.shape[:2]
        H -= H % 2; W -= W % 2
        arr = arr[:H, :W]        # the slices below run over the SOURCE; an odd
                                 # height raised a broadcast error. See fits.py.
        cfa = np.empty((H, W), np.float32)
        cfa[0::2, 0::2] = arr[0::2, 0::2, 0]   # R
        cfa[0::2, 1::2] = arr[0::2, 1::2, 1]   # G
        cfa[1::2, 0::2] = arr[1::2, 0::2, 1]   # G
        cfa[1::2, 1::2] = arr[1::2, 1::2, 2]   # B
        self.bayer = cfa
        self.white_level = self.SCALE
        self.sat_level = 0.97 * self.SCALE
        self.linear_max = None          # a TIFF carries no sensor linearity limit
        self.sat_source = "container"   # ... nor a white level: this is 1.0
        self.daylight_wb = np.ones(3, np.float32)
        self.cam2rgb = np.eye(3, dtype=np.float32)

    @property
    def shape(self):
        return self.bayer.shape


def open_frame(path):
    from .fits import is_fits, FitsFrame
    if is_fits(path):
        return FitsFrame(path)
    if os.path.splitext(path)[1].lower() in TIFF_EXTS:
        return TiffFrame(path)
    return RawFile(path)


# ---------- hot / dead pixel repair ----------

def _outlier_flags(bayer, k=6.0, read_noise=None):
    """Pixels deviating from the median of their same-colour neighbours by more
    than k sigma, with sigma from a photon+read noise model.

    `read_noise` is a per-Bayer-offset sequence in ADU, measured from bias
    frames. When it is given the intercept of the model is KNOWN instead of
    fitted. That matters because the fit estimates gain and read noise together
    from one frame, and the two trade off against each other: the shortest
    tier of an eclipse bracket has almost no signal range to separate them
    with, so a poorly determined intercept moves the threshold everywhere. A
    bias frame measures the intercept with the shutter closed, which is the
    only condition under which it is measurable on its own.
    """
    from scipy import ndimage
    flag = np.zeros(bayer.shape, bool)
    for oy in (0, 1):
        for ox in (0, 1):
            sub = np.asarray(bayer[oy::2, ox::2], np.float32)
            med = ndimage.median_filter(sub, size=3, mode="nearest")
            d = sub - med
            lo, hi = np.percentile(med, 5), np.percentile(med, 99)
            g = rn2 = None
            _rn = None
            xs = ys = None
            if read_noise is not None:
                try:
                    _rn = float(read_noise[oy * 2 + ox])
                except (TypeError, IndexError, ValueError):
                    _rn = None
            if _rn is not None and _rn > 0 and hi > lo:
                # THE MEDIAN FILTER'S OWN NOISE IS IN `d` TOO. d = pixel minus
                # the median of its 3x3 same-colour neighbourhood, so its
                # variance is the pixel's plus the median's, and for a 3x3
                # median of independent samples the latter is about 1/3 of one
                # sample's. Ignoring it would make sigma ~15% too small and
                # over-flag everywhere. The gain is still fitted -- only the
                # intercept is pinned.
                _c = 1.0 + 1.0 / 3.0
                bins = np.linspace(lo, hi, 12)
                idx = np.digitize(med.ravel(), bins)
                dr, mr = d.ravel(), med.ravel()
                xs, ys = [], []
                for b in range(1, len(bins)):
                    m = idx == b
                    if m.sum() > 2000:
                        xs.append(float(mr[m].mean()))
                        ys.append(float((1.4826 * np.median(np.abs(dr[m]))) ** 2))
                if len(xs) >= 3:
                    xs = np.asarray(xs); ys = np.asarray(ys)
                    # one free parameter: ys - c*rn^2 = g * c * xs
                    _y = ys - _c * _rn * _rn
                    _x = _c * np.maximum(xs, 0.0)
                    _den = float(_x @ _x)
                    g = max(float(_x @ _y) / _den, 0.0) if _den > 0 else 0.0
                    rn2 = _c * _rn * _rn
                    g = g * _c
                    # WHAT THE FIT WOULD HAVE SAID, recorded so a run can
                    # print both. The bench cannot settle which number is
                    # right on a given sensor: on a synthetic tier the fit
                    # comes out HIGH (pinning found 97 defects against 80),
                    # on the 600 mm reference set it must have come out LOW (the
                    # count fell 3675 -> 1848). Same failure, opposite sign,
                    # so the only way to know is to look at the two numbers
                    # for the data in hand.
                    try:
                        _A = np.stack([np.asarray(xs), np.ones(len(xs))], 1)
                        _sol = np.linalg.lstsq(_A, np.asarray(ys), rcond=None)[0]
                        _fitrn = float(np.sqrt(max(float(_sol[1]), 0.0) / _c))
                    except Exception:
                        _fitrn = float("nan")
                    _outlier_flags.last_fit_rn[oy * 2 + ox] = _fitrn
            if g is None and hi > lo:
                bins = np.linspace(lo, hi, 12)
                idx = np.digitize(med.ravel(), bins)
                dr, mr = d.ravel(), med.ravel()
                xs, ys = [], []
                for b in range(1, len(bins)):
                    m = idx == b
                    if m.sum() > 2000:
                        xs.append(float(mr[m].mean()))
                        ys.append(float((1.4826 * np.median(np.abs(dr[m]))) ** 2))
                if len(xs) >= 3:
                    A = np.stack([np.array(xs), np.ones(len(xs))], 1)
                    sol = np.linalg.lstsq(A, np.array(ys), rcond=None)[0]
                    g, rn2 = max(float(sol[0]), 0.0), max(float(sol[1]), 1e-6)
            if g is None:
                g, rn2 = 0.0, max(float(np.var(d)), 1e-6)
            sig = np.sqrt(np.maximum(g * np.maximum(med, 0) + rn2, 1e-6))
            flag[oy::2, ox::2] = np.abs(d) > k * sig
    return flag


#: What the single-frame fit put the read noise at, per Bayer offset, on the
#: last call that was ALSO given a measured value. Diagnostic only.
_outlier_flags.last_fit_rn = [float("nan")] * 4


def hot_pixel_map(bayers, k=6.0, frac=0.6, read_noise=None):
    """Sensor defects sit at the same photosite in every frame. Detected on the
    SHORTEST exposure tier, where the sky is essentially black, so real sky
    objects (stars, planets) are far below the noise and cannot be flagged."""
    n = len(bayers)
    if n == 0:
        return None
    # A frame shot in a crop mode, or a stray file from another body, would
    # otherwise raise a bare broadcast error here or an IndexError later in
    # repair_hot, with nothing in the message about frame sizes.
    shp = bayers[0].shape
    if any(b.shape != shp for b in bayers):
        raise RuntimeError(
            "the frames in this tier are not all the same size "
            f"({', '.join(sorted({str(b.shape) for b in bayers}))}) -- mixing "
            "crop modes, sensor sizes or cameras in one bracket is not supported")
    _outlier_flags.last_fit_rn = [float("nan")] * 4
    kk = k if n > 1 else k + 2.0      # single frame: no cross-frame vote, be stricter
    cnt = np.zeros(bayers[0].shape, np.uint8)
    for b in bayers:
        cnt += _outlier_flags(b, kk, read_noise)
    need = 1 if n == 1 else max(2, int(np.ceil(frac * n)))
    return cnt >= need


def repair_hot(bayer, hot):
    """Replace flagged photosites by the median of their same-colour neighbours."""
    from scipy import ndimage
    if hot is None or not hot.any():
        return bayer
    if hot.shape != bayer.shape:
        return bayer                      # defect map built on a different size
    for oy in (0, 1):
        for ox in (0, 1):
            hs = hot[oy::2, ox::2]
            if not hs.any():
                continue
            sub = bayer[oy::2, ox::2]
            med = ndimage.median_filter(np.asarray(sub, np.float32), size=3, mode="nearest")
            sub[hs] = med[hs]
    return bayer
