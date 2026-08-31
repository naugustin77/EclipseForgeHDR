"""Raw decoding and metadata for exposure-bracketed eclipse sequences."""
from __future__ import annotations
import os
import numpy as np

RAW_EXTS = {".rw2", ".raf", ".nef", ".cr2", ".cr3", ".arw", ".orf", ".dng", ".pef"}


def list_raws(folder):
    from .fits import is_fits
    out = []
    for name in sorted(os.listdir(folder)):
        ext = os.path.splitext(name)[1].lower()
        if ext in RAW_EXTS or ext in {".tif", ".tiff"} or is_fits(name):
            out.append(os.path.join(folder, name))
    return out


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


class RawFile:
    """Decoded bayer data normalized to RGGB layout, plus color info."""

    def __init__(self, path):
        import rawpy
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
            self.daylight_wb = np.asarray(raw.daylight_whitebalance[:3], np.float32)
            if not np.isfinite(self.daylight_wb).all() or self.daylight_wb[1] <= 0:
                # some bodies report zeros here; without this the whole merge
                # becomes inf/nan and gets cached that way
                self.daylight_wb = np.ones(3, np.float32)
            else:
                self.daylight_wb = self.daylight_wb / self.daylight_wb[1]
            cm = np.asarray(raw.rgb_xyz_matrix, np.float32)[:3, :3]  # XYZ -> cam
            self.cam2rgb = _cam2srgb(cm)
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
        _obs = float(self.bayer.max()) if self.bayer.size else 0.0
        if _obs > self.sat_level:
            self.sat_level = _obs * 0.975
            self.white_level = _obs + black

    @property
    def shape(self):
        return self.bayer.shape


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
        arr = arr.astype(np.float32)
        if arr.max() > 1.001:
            arr /= 65535.0 if arr.max() > 255.5 else 255.0
        arr = np.clip(arr, 0, 1)
        linear = (os.environ.get("ECLIPSEFORGE_TIFF_LINEAR") == "1"
                  or "linear" in os.path.basename(path).lower())
        if not linear:  # inverse sRGB
            a = arr
            arr = np.where(a <= 0.04045, a / 12.92, ((a + 0.055) / 1.055) ** 2.4)
        arr *= self.SCALE
        H, W = arr.shape[:2]
        H -= H % 2; W -= W % 2
        cfa = np.empty((H, W), np.float32)
        cfa[0::2, 0::2] = arr[0::2, 0::2, 0]   # R
        cfa[0::2, 1::2] = arr[0::2, 1::2, 1]   # G
        cfa[1::2, 0::2] = arr[1::2, 0::2, 1]   # G
        cfa[1::2, 1::2] = arr[1::2, 1::2, 2]   # B
        self.bayer = cfa
        self.white_level = self.SCALE
        self.sat_level = 0.97 * self.SCALE
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

def _outlier_flags(bayer, k=6.0):
    """Pixels deviating from the median of their same-colour neighbours by more
    than k sigma, with sigma from a fitted photon+read noise model."""
    from scipy import ndimage
    flag = np.zeros(bayer.shape, bool)
    for oy in (0, 1):
        for ox in (0, 1):
            sub = np.asarray(bayer[oy::2, ox::2], np.float32)
            med = ndimage.median_filter(sub, size=3, mode="nearest")
            d = sub - med
            lo, hi = np.percentile(med, 5), np.percentile(med, 99)
            g = rn2 = None
            if hi > lo:
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


def hot_pixel_map(bayers, k=6.0, frac=0.6):
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
    kk = k if n > 1 else k + 2.0      # single frame: no cross-frame vote, be stricter
    cnt = np.zeros(bayers[0].shape, np.uint8)
    for b in bayers:
        cnt += _outlier_flags(b, kk)
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
