"""Recombination render from cached layers, preview generation, exports."""
from __future__ import annotations
import os, sys, json
import numpy as np
from scipy import ndimage
from scipy.special import erf
from PIL import Image

# NAFE REMOVED (2026-10-08, Nico). With Denoise Off as the default the rank
# equaliser turned unsmoothed noise near the limb into 4-sd blobs (18 between
# 1.02 and 1.47 R on the 600 mm set) and added nothing measurable to the
# streamers. The layer is no longer built, loaded or mixed; True brings it back.
NAFE_ENABLED = False

# Starting points, not truths. These are the settings the reference bracket was
# worked to by eye once the layers were behaving (0.12.0); every one of them is
# a taste call and every one is a slider.
DEFAULTS = {
    "fnMix": 0.17, "detailGain": 1.0, "baseLift": 0.255, "envGamma": 1.0,
    "logK": 0.0,
    # Hill's chain. hillMix 0 renders exactly as every build before 0.22.64;
    # hill0..hill4 are his 100:60:20:10 ratios at 2:4:8:16 px, extended one
    # octave, and hillGain is the one master that scales the set.
    # hillDenoise 3.0 (was 0): measured on the reference set at threshold 0 the
    # masks beyond 2 R are pure photon noise, amplified into visible grain; at
    # 3 sigma they are zero there and the inner corona keeps its structure.
    # 0.23.7 defaults are what the 600 mm reference set settled on, viewed as
    # Hill views his (masks on mid grey, half size): raw base, no threshold,
    # the 2 px mask off -- on this set it carries noise only -- and the
    # ladder starting at 4 px.
    "hillMix": 0.0, "hillGain": 0.025, "hillLogK": 3.0, "hillDenoise": 0.0,
    # hillBase: how much of the log-mapped base the Partial-conv view shows
    # under the masks. 1 = Hill's im_enhanced; 0 = the masks alone on mid
    # grey, which is his "E" panel, the layer the blend overlays, and the
    # honest way to judge them -- the default since 0.23.7.
    "hillBase": 0.0,
    "hill0": 0.0, "hill1": 1.0, "hill2": 0.6, "hill3": 0.2, "hill4": 0.1,
    "radialFlatten": 0.5, "mgnContrast": 0.15, "fnCompress": 0.9,
    "structure": 0.0, "clarity": 0.39, "smoothing": 0.25, "pelGain": 0.15,
    "nafeMix": 0.0, "innerMix": 0.31, "innerDenoise": 0.4, "innerDim": 0.05, "promGain": 0.4,
    # How much of the detail inside the prominence gate comes from the
    # prominence's own layer rather than from the corona filters. 0 reproduces
    # every build before 0.18 exactly. See detail.prominence_detail.
    "promDetail": 0.7,
    # Where the MGN layer sits between all six scales weighted alike (0, every
    # build before 0.22) and the three FINE scales alone (1). See the long note
    # in detail.build_layers: our gains are noise-normalisation constants, not
    # Hill's 1:0.6:0.2:0.1 amplification ladder, and the difference is slabs
    # against filaments.
    #
    # DEFAULT 1.0 -- chosen by Nico from a blind-ish four-way contact sheet at
    # matched contrast, against my metric, which rates it 44% WORSE at the limb.
    # amp*coh cannot see "delicate": it rewards radially coherent structure and
    # the coarse scales carry plenty. Drag to 0 for the old look.
    "detailScale": 1.0,
    # Prominence structure carried in green/blue instead of luminance -- see
    # the long note where it is applied. OFF by default: the measurement is
    # sound and the mechanism works, but at +0.6 it makes dense prominence
    # material go pink, and pink prominences look wrong. Negative values run it
    # the other way (dense material goes deeper red), which is the direction
    # worth trying before this is called a dead end.
    "promChroma": 0.0,
    # THE 0.24 PROMINENCE LAYER (promlayer.py). Used when the work directory has
    # one; the three prominence controls above are then inactive, because the
    # corona they acted on no longer contains the prominences.
    #   promLayer   strength of the layer (opacity scale); 0 = no prominences
    #   promStretch the layer's own curve (lower = faint parts brighter)
    #   promRim     the chromosphere rim along the limb, 0..1
    #   promSat     the layer's colour saturation
    #   promFade    how far the corona's Hill detail is eased out under the
    #               layer (1 = fully where the layer covers it)
    "promLayer": 1.0, "promStretch": 0.5, "promRim": 1.0, "promSat": 1.3, "promFade": 1.0,
    # RAW-PATH WHITE BALANCE, SOLVED RATHER THAN EYEBALLED.
    #
    # These were 0.9 / 1.205, a by-eye correction from 0.12.0, and tint at 1.205
    # is a straight 20% GREEN MULTIPLY. After the unit-luminance renormalisation
    # below (green carries 0.7152 of luminance) the net is green +5% and red and
    # blue -13% relative -- the green-grey sky Nico reported, and the same cast
    # the IMPORT_DEFAULTS note below had already measured and fixed for imports
    # while the raw path was left alone. Half a fix for two years.
    #
    # Solvable exactly, because the far sky is neutral before these are applied:
    # `ratio` is driven to 1.0 where the signal is near the noise floor, and the
    # bgNeutral division is weighted by `cconf`, which measures 0.015 out there.
    # So the rendered sky chroma is just (temp, tint, 1/temp) renormalised, and
    #     R/G = temp/tint        B/G = 1/(temp*tint)
    #
    # Nico's own PixInsight version -- "lowered the green curve a bit and raised
    # the blue curve a little" -- measured off his JPEG, linearised, normalised
    # to green: sky R/G 1.064, B/G 1.884. Inverting the two expressions gives
    # temp 0.752, tint 0.706, which reproduces that sky exactly.
    #
    #              temp   tint    sky R/G   sky B/G
    #   was        0.900  1.205     0.747     0.922      green-dominant
    #   now        0.752  0.706     1.064     1.884      his target, exactly
    #
    # Relative to green this is R/G x1.425 and B/G x2.043 -- his two curve moves,
    # as numbers. It cools the corona too, which is why `satur` is the other half
    # of what he did: he pushed saturation to restore the warm/cool separation.
    #
    # BOTH BACK TO 1.0 IN 0.22.23, at Nico's request, and the derivation above
    # deserves the correction rather than quiet deletion.
    #
    # It rests on "the far sky is neutral before these are applied". Measured on
    # hdr_rgb, the far sky is not neutral on ANY of the five datasets, and the
    # five do not agree with each other:
    #
    #     set                    far sky R/G   B/G
    #     Clifton 2024 560mm         1.244    0.908
    #     Clifton 2026 360mm         1.206    0.880
    #     Clifton 2026 250mm         1.116    0.757
    #     Nico Lumix 600mm           0.992    0.647
    #     Nico Sony                  0.553    1.328
    #
    # So a pair of numbers solved against one photographer's rendering of one
    # eclipse cannot transfer: they landed near his target on his set because
    # his own sky chroma was folded into the fit, and they render violet on a
    # Canon whose sky chroma differs. 1.0/1.0 applies no cast at all and leaves
    # the decision where it belongs -- with the two sliders and the person
    # looking at the picture.
    #
    # The principled fix is a corona-referenced white balance: the K-corona is
    # Thomson-scattered photospheric light and is white by physics, so the
    # 1.05-1.6 R annulus is a defensible white reference. `corona_white_balance`
    # already does exactly this, and is gated to files with no camera WB, so it
    # never runs on a raw bracket. Opening that gate is the real answer and is
    # not this change.
    "temp": 1.0, "tint": 1.0, "bgNeutral": 1.0, "satur": 1.0, "hlCompress": 0.0, "hlDesat": 0.0,
    "outGamma": 1.0, "bgBlack": 0.02,
    # Final sharpening (lab, Oct 2026): Nimmervoll's Photoshop finish with a
    # noise gate. All off by default -- it is a finish, not a default look.
    "finalMid": 0.0, "finalDark": 0.0, "finalLight": 0.0, "finalNoise": 3.0,
    "discLevel": 0.045, "discTrim": 0.0, "earthShine": 0.0,
    # block 3a: zoned Hill gains, tone curve, contained whites, disc edge
    "hillZone": 0.0, "hillZoneR": 1.8,
    "toneCurve": [[0.0, 0.0], [1.0, 1.0]], "containWhites": 1.0,
    "discSoft": 12.0,
    "ringBlend": 0.0, "ringScale": 1.0, "ringDX": 0.0, "ringDY": 0.0,
    # the look controls (lab 0.24), all off -- see look()
    "coronaLevel": 0.0, "skyLevel": 0.0, "skyColour": 0.0, "moonColour": 0.0,
    # THE MOON GROUP (lab 0.24): the Moon's own tone curve -- the main curve no
    # longer reaches the disc -- and the colour its tint heads towards, picked
    # on the page (display RGB; the default is the look's blue, LOOK_BLUE).
    "moonCurve": [[0.0, 0.0], [1.0, 1.0]],
    "moonHueR": 0.579, "moonHueG": 0.756, "moonHueB": 1.0,
    # the sky lift's colour (lab 2026-10-08): the colour the lift heads towards
    # at Sky lift colour 1, picked as for the Moon; the default is the old blue
    # axis exactly, so every look before it renders unchanged
    "skyHueR": 0.579, "skyHueG": 0.756, "skyHueB": 1.0,
    # COLOUR MODE (lab 2026-10-08): 0 = colour from the data (as before),
    # 1 = monochrome (no colour anywhere; the composite exports as one
    # channel), 2 = mono + tints (grey corona; sky, Moon and prominences take
    # the picked colours). promHue*: the prominence colour used in mode 2.
    # promPeak: the prominences' brightest channel is held to this (1 = off).
    "colourMode": 0, "promHueR": 1.0, "promHueG": 0.70, "promHueB": 0.75, "promPeak": 1.0,
}

# An IMPORTED image has already been through somebody's colour management --
# that is what makes it an import rather than a stack -- so `temp` and `tint`
# start neutral there.
#
# WHY: on the raw path the merge applies the camera's daylight white balance
# and colour matrix, and 0.9/1.205 is the by-eye correction that sat on top of
# THAT, for THAT camera, under a sun 7 degrees up. It is a taste call about one
# specific pipeline, not a property of coronae. Applied a second time to a file
# that is already balanced it is simply a tint.
#
# MEASURED on Val Italo's 16-bit sRGB stack (4680x3132), median chroma
# normalised to unit luminance, source against the rendered composite:
#
#                       source           as shipped        temp/tint 1.0
#   limb  1.02-1.15 R   1.058 .990 .925  0.832 1.062 .885  1.046 .996 .901
#   inner 1.15-1.5  R   1.088 .985 .889  0.853 1.058 .853  1.074 .991 .868
#   mid   1.5 -2.5  R   1.155 .973 .811  0.902 1.050 .791  1.133 .980 .803
#
# Green becomes the strongest channel in every shell -- the teal cast he
# reported. Neutral reproduces the source to within 2%, the remainder being
# `bgNeutral` (1.5%). The per-channel sky-gradient division, the other
# suspect, was measured on the same file at 0.1% or less: its tilt is +89 deg,
# so it cancels in a radial median.
IMPORT_DEFAULTS = dict(DEFAULTS, temp=1.0, tint=1.0)


def defaults_for(mode):
    """Slider starting points for a workdir built by `mode` ('import' or not)."""
    return IMPORT_DEFAULTS if mode == "import" else DEFAULTS


# WHICH SETTINGS TRAVEL TO A DIFFERENT STACK, AND WHICH DO NOT.
#
# A settings file is for carrying a look from one bracket to another, and for
# sending one to somebody else. That only works for the settings that mean the
# same thing on data they have never seen.
#
# LOCAL settings are tied to the geometry or the level of the one image they
# were set on. Applied blind to another stack they do not carry a look across,
# they carry a mistake across, and the result is worse than the defaults -- so
# they are skipped unless the file came from the very folder being loaded into.
#
#   discLevel, discTrim   the occulting disc: radius and margin in pixels
#   ring*                 the diamond ring, registered to one contact frame
#   bgBlack               a black point in this image's own level units
#
# Everything else is dimensionless or is normalised against something the app
# measures per dataset, so the number means the same amount of the same thing
# anywhere. Two are worth naming because they were NOT portable until recently:
#
#   hillGain    normalised since 0.22.84 against the 2 px mask's STRUCTURE rms
#               rather than its total rms. Before that the same value meant
#               "this much noise" on a noisy set and "this much detail" on a
#               clean one -- 0.63 against 1.00 between two of the test sets --
#               so it did not transfer at all. A settings file written by an
#               older build therefore carries a hillGain that means something
#               else, which is why the file records the version that wrote it.
#   coronaNeutral, bgNeutral
#               these scale a correction that is re-measured on every dataset,
#               so the SLIDER travels even though the correction it scales does
#               not, which is exactly what portability means here.
#
# baseLift, logK, hillLogK and envGamma are the awkward ones. They are stretch
# parameters, and the coronal range runs 6.7 to 8.2 EV across the test sets, so
# the same stretch lands differently. They travel -- a look is mostly stretch
# and refusing to carry them would make the feature pointless -- but they are
# the first thing to reach for when a loaded look is close but not right.
_LOCAL_KEYS = frozenset({
    "discLevel", "discTrim", "bgBlack",
    "ringBlend", "ringScale", "ringDX", "ringDY",
})


def classify_settings(keys=None):
    """({portable}, {local}) over the render parameters."""
    ks = set(DEFAULTS) if keys is None else set(keys)
    return frozenset(ks - _LOCAL_KEYS), frozenset(ks & _LOCAL_KEYS)




# Orientation is the one property of the output that no measurement depends on.
#
# The limb fit is a circle, MGN, FNRGF, NAFE and Pellett are all radial or
# tangential about that circle, and the Bayer decode was settled at read time.
# Turn the finished picture any way up and every number in the report is the
# number it was. So this belongs at the very last step, applied to the array on
# its way out -- never at read time, where it would change what is cached and
# cost a full re-stack to undo something purely cosmetic.
#
# It exists because a FITS file cannot say which way up the camera was. ROWORDER
# describes row order, not camera rotation, and there is no orientation keyword
# in FITS at all -- a portrait-shot frame arrives on its side and nothing in a
# corona can tell you so. A rotationally symmetric subject has no up.
_ORIENT = {"": lambda a: a,
           "flipv": lambda a: a[::-1],
           "fliph": lambda a: a[:, ::-1],
           "180":   lambda a: a[::-1, ::-1],
           "cw":    lambda a: np.rot90(a, -1),
           "ccw":   lambda a: np.rot90(a, 1)}


def apply_orient(a, orient):
    """Turn a finished image. Lossless: every one of these is a view permutation."""
    f = _ORIENT.get(str(orient or "").lower())
    return a if f is None else np.ascontiguousarray(f(a))


def _decim(a, q):
    """Decimate by q with an area average. See the note at the call site: plain
    subsampling of a fine detail layer aliases into coarse mottle."""
    if q <= 1:
        return np.asarray(a, np.float32)
    a = np.asarray(a, np.float32)
    H, W = a.shape[:2]
    h, w = H // q * q, W // q * q
    if h < q or w < q:
        return a[::q, ::q]
    out = a[:h, :w].reshape(h // q, q, w // q, q, *a.shape[2:]).mean(axis=(1, 3))
    return np.ascontiguousarray(out, np.float32)


def _fill_disc(a, cy, cx, Rm):
    """Continue the layer across the disc edge by reflecting it inward.

    The detail layers are a flat 0.5 inside the disc. Any blur taken across that
    boundary -- which is what Clarity and Grain smoothing subtract -- averages
    the plateau together with the corona and comes back too low just outside the
    limb, so the unsharp difference prints a bright rim. Measured on the
    reference layer, Clarity at 0.39 amplified the rim by 25% and at 1.0 by 64%;
    reflecting first removes most of that (0.0588 -> 0.0494 and 0.0772 -> 0.0532).

    Only the blurred VARIANTS are built from the filled copy. The layer itself
    is untouched, so nothing invented inside the disc reaches the output.
    """
    H, W = a.shape
    yy = np.arange(H, dtype=np.float32)[:, None] - cy
    xx = np.arange(W, dtype=np.float32)[None, :] - cx
    r = np.sqrt(yy * yy + xx * xx)
    inside = r < Rm
    if not inside.any():
        return a
    rr = np.where(inside, 2.0 * Rm - r, r)
    th = np.arctan2(yy, xx)
    ys = np.clip(cy + rr * np.sin(th), 0, H - 1)
    xs = np.clip(cx + rr * np.cos(th), 0, W - 1)
    out = a.copy()
    out[inside] = ndimage.map_coordinates(a, [ys[inside], xs[inside]],
                                          order=1, mode="nearest")
    return out


# Amplification slider reference: see Layers.hill_ref. The ratio, on the 600 mm
# reference stack, between the 32 px mask's spread at Rmask-1.5 R and the old
# normaliser (0.25 x the 2 px mask's spread over the frame).
HILL_REF_C = 6.25


class Layers:
    """Cached full-res layers, plus a decimated copy for interactive previews."""

    def __init__(self, wd, preview_decim=4):
        self.wd = wd
        # How these products were built. Only the colour defaults read it, but
        # they have to: an import is already white-balanced (see IMPORT_DEFAULTS).
        self.mode = None
        try:
            self.mode = json.load(open(os.path.join(wd, "opts.json"))).get("mode")
        except Exception:
            pass
        geo = json.load(open(os.path.join(wd, "geometry.json")))
        self.cy, self.cx, self.R = geo["cy"], geo["cx"], geo["R"]
        self.Rmask = float(geo.get("Rmask", self.R + 4.0))
        self.limb_prof = geo.get("limb_prof")
        self.limb_margin = float(geo.get("limb_margin", self.Rmask - self.R))
        lum = np.load(os.path.join(wd, "hdr_lum.npy"))
        # THE 0.24 PROMINENCE LAYER: the picture is built on the corona without
        # its prominence light -- but its WHITE POINT stays the merge's, so taking
        # the brightest thing in the frame out does not re-stretch the corona.
        from . import promlayer as _plm
        self.has_promlayer = False
        _cl = os.path.join(wd, _plm.CORONA_LUM)
        _use_corona = (_plm.layer_mode() and os.path.exists(_cl)
                       and os.path.exists(os.path.join(wd, _plm.CORONA_RGB))
                       and os.path.exists(os.path.join(wd, _plm.LIGHT_FILE)))
        _hi_merge = float(np.percentile(lum, 99.97))
        if _use_corona:
            lum = np.load(_cl)
        self.shape = lum.shape
        H, W = lum.shape
        yy = np.arange(H, dtype=np.float32)[:, None] - self.cy
        xx = np.arange(W, dtype=np.float32)[None, :] - self.cx
        r = np.sqrt(yy * yy + xx * xx)
        # BLACK POINT FROM THE SKY, NOT FROM A PERCENTILE (fixed in 0.11.5)
        #
        # This used lo = percentile(lum, 2). On a wide field the sky IS most of
        # the frame, so that lands the sky within a noise sigma of the clip
        # point -- and then `xn` is a small difference between two nearly equal
        # numbers, which turns a tiny real variation in the sky into an
        # enormous one on screen. Measured on the reference set:
        #
        #                         corner   mid-edge   ratio   sky clipped
        #   the data itself                            1.030
        #   lo = p2          Bg   0.0478   0.0637      1.332      2.4%
        #   lo = sky - 5 sig Bg   0.1253   0.1284      1.025      0.0%
        #
        # A 3% brightness difference across the frame was being shown as 33%,
        # and 2.4% of the sky was crushed to pure black. That -- not any detail
        # layer -- is the "vignetting" that survived every layer being set to
        # zero.
        #
        # It also un-pins the radial flatten control: with the sky that close to
        # the floor, `rprof` sat on its 0.12 clamp everywhere in the outer field
        # (corner and mid-edge both exactly 0.1200), so radial flattening did
        # nothing out there. After the fix it reads 0.165-0.170 and works.
        # Net through both: 1.332 -> 1.013, against 1.030 in the data.
        hi = _hi_merge
        _far = r > 2.5 * max(float(self.R), 1.0)
        lo = float(np.percentile(lum, 2))
        if _far.sum() > 20000:
            _sv = lum[_far]
            _sm = float(np.median(_sv))
            _ss = 1.4826 * float(np.median(np.abs(_sv - _sm)))
            if np.isfinite(_sm) and _ss > 0:
                # never ABOVE the 1st percentile, so this can only ever clip
                # less than the old rule, never more
                lo = min(_sm - 5.0 * _ss, float(np.percentile(lum, 1)))
        del _far
        if not np.isfinite(lo) or hi <= lo:
            lo = float(np.percentile(lum, 2)); hi = _hi_merge
        self.black_point = lo
        xn = np.clip((lum - lo) / (hi - lo), 0, 1)
        Bg = xn ** (1 / 3.0)
        # Sized from the image, not capped at a literal 6000 px: on a larger
        # sensor, or with the disc off-centre, every radius past the cap
        # collapsed into one bin and the radial-flatten control quietly stopped
        # working in the outer field.
        #
        # Built by sorting once instead of scanning the whole frame 1500 times.
        # The old loop was ~2 minutes on a 45 MP frame, paid on every load.
        nr = int(r.max()) + 2
        rid = np.clip(r.astype(np.int32), 0, nr - 1)
        nb = (nr + 3) // 4
        bins = rid.ravel() // 4
        order = np.argsort(bins, kind="stable")
        bs = bins[order]
        vs = Bg.ravel()[order]
        edges = np.searchsorted(bs, np.arange(nb + 1))
        prof = np.zeros(nr, np.float32)
        for kb in range(nb):
            a, b = edges[kb], edges[kb + 1]
            if b - a > 300:
                prof[kb * 4:min(kb * 4 + 4, nr)] = np.median(vs[a:b])
        del order, bs, vs, bins
        ok = prof > 0
        if ok.any():
            idx = np.flatnonzero(ok)
            prof[~ok] = np.interp(np.flatnonzero(~ok), idx, prof[idx])
        prof = ndimage.gaussian_filter1d(prof, 25)
        prof /= max(prof.max(), 1e-6)
        # INTERPOLATE, do not round the radius. `prof[rid]` with
        # rid = r.astype(int) makes the divisor piecewise-constant inside each
        # one-pixel-wide annulus, so radial flattening divides the picture by a
        # staircase and prints concentric rings. Measured on a 600 mm frame,
        # the step between neighbouring integer radii is 0.28% of the local
        # level near the limb -- small, but perfectly circular and coherent,
        # which is what both the eye and the detail filters pick out. Reported
        # as "radial flatten introduces some radial pattern; with it off, gone".
        rprof = np.maximum(
            np.interp(r, np.arange(nr, dtype=np.float32), prof),
            0.12).astype(np.float32)
        # Kept so the page can evaluate the profile itself. Sending it as a
        # full-resolution 8-BIT image was the second ring source: in the outer
        # field only 1-9% of radii differed at all, leaving wide flat annuli
        # separated by jumps of up to 2%. It is one number per radius, so there
        # is no reason to rasterise it.
        self.rprof_1d = prof.astype(np.float32)

        mg = np.load(os.path.join(wd, "mgn.npy"))
        mg = np.clip((mg - np.percentile(mg, 1)) /
                     (np.percentile(mg, 99.7) - np.percentile(mg, 1)), 0, 1)
        _mfp = os.path.join(wd, "mgn_fine.npy")
        self.has_mgn_fine = os.path.exists(_mfp)
        if self.has_mgn_fine:
            mgf = np.load(_mfp)
            mgf = np.clip((mgf - np.percentile(mgf, 1)) /
                          (np.percentile(mgf, 99.7) - np.percentile(mgf, 1)), 0, 1)
        D = np.load(os.path.join(wd, "fnrgf.npy"))
        nfp = os.path.join(wd, "nafe.npy")
        # ABSENT, NOT FLAT (0.22.78). A workdir stacked before NAFE existed got
        # a constant 0.5 here, which is not "no NAFE" -- the default nafeMix of
        # 0.15 then blended 15% of a flat field into mg and fn and took ~5% of
        # the local contrast out of the picture, with the slider still moving
        # and nothing in the log. Leaving the key out makes the blend below skip
        # it, the same way a missing promdet is handled.
        self.has_nafe = NAFE_ENABLED and os.path.exists(nfp)
        nfl = np.load(nfp) if self.has_nafe else None
        fnl = np.clip(D / 8.0 + 0.5, 0, 1)     # raw sigma units, ±4σ window
        del D
        inner = np.load(os.path.join(wd, "inner.npy"))
        inner0 = np.load(os.path.join(wd, "inner0.npy"))
        # THE MOON LAYER (moonlayer.py, lab 0.24) when the work directory has
        # one; the old Earthshine-option layer otherwise
        _mlp = os.path.join(wd, "moon_layer.npy")
        ep = _mlp if os.path.exists(_mlp) else os.path.join(wd, "earth.npy")
        self.has_earth = os.path.exists(ep)
        earth = (np.load(ep).astype(np.float32) if self.has_earth
                 else np.full(lum.shape, 0.5, np.float32))
        if earth.shape != lum.shape:
            self.has_earth = False
            earth = np.full(lum.shape, 0.5, np.float32)
        gate = np.load(os.path.join(wd, "prom.npy"))
        pelp = os.path.join(wd, "pellett.npy")
        pel = np.load(pelp) if os.path.exists(pelp) else np.full(lum.shape, 0.5, np.float32)
        from .pipeline import load_big
        hdr = load_big(os.path.join(wd, _plm.CORONA_RGB if _use_corona else "hdr_rgb.npy"))
        Ls = ndimage.gaussian_filter(lum, 6)
        # Colour is only meaningful where something was actually detected.
        #
        # The divisor used to be floored at a literal 1e-3 in linear luminance
        # units, which means nothing across cameras and exposures. On a short,
        # noisy bracket the far field sits near zero, each channel's ratio slams
        # into its [0.2, 3.0] clip in a random direction, and the render paints
        # the sky as fully saturated RGB speckle. Floor at the sky's own noise
        # instead, and fade the chroma to neutral where the signal is not above
        # it -- an undetected sky has no colour to report.
        _sk = Ls[r > 0.80 * float(r.max())]
        if _sk.size > 5000:
            _bg = float(np.median(_sk))
            _nz = 1.4826 * float(np.median(np.abs(_sk - _bg))) + 1e-9
        else:
            _bg, _nz = 0.0, max(float(np.median(Ls)) * 1e-3, 1e-9)
        _floor = max(_bg + 2.0 * _nz, 1e-9)
        # SMOOTHSTEP, AND OVER A WIDER SPAN. This used to be
        # np.clip((Ls - _floor) / (4 * _nz), 0, 1) -- a linear ramp with the
        # corners left on, so the derivative jumped at both ends of it.
        #
        # In the far field the luminance falls slowly with radius, so four
        # sigmas of it occupy a narrow annulus, and a chroma that is fading in a
        # straight line and then stops dead draws a RING there. Nico found it by
        # setting the white balance to None, which puts a strong green inside
        # against the exactly-neutral far field and makes the boundary obvious;
        # it is present at every white balance, just quieter. His check is the
        # one that settles what it is: the ring is not in the single raws even
        # pushed hard in curves, and it cannot be -- this fade exists only here.
        #
        # 8 sigmas instead of 4, with a smoothstep, so both knees are C1. The
        # fade still reaches exactly neutral at the floor, which is the point of
        # it: an undetected sky has no colour to report, and letting its chroma
        # noise through is what this was written to stop.
        _t = np.clip((Ls - _floor) / (8.0 * _nz), 0.0, 1.0)
        _conf = (_t * _t * (3.0 - 2.0 * _t)).astype(np.float32)
        del _t
        # WHERE THE FADE GOES WHEN THE SKY WAS SUBTRACTED (lab 0.24). With the
        # sky in the merge, the far field's colour IS the sky's and fading it to
        # neutral is harmless. With the sky subtracted it is the corona's own
        # colour (and a floor in that colour) all the way out -- on Nico's 600 mm
        # set R:G:B 1.59 : 0.88 : 0.44 of luminance, a low Sun's reddened
        # corona -- and fading THAT to neutral put a hue edge around the whole
        # corona at 3-4.5 R: a brown ring against grey under the camera's
        # balance, a grey ring against blue under the corona white. So when the
        # pipeline reports the sky subtracted, the fade goes to the far field's
        # measured colour instead, and the colour is continuous through the edge.
        _fade_to = np.ones(3, np.float32)
        try:
            _ssub = (json.load(open(os.path.join(wd, "report.json"))).get("sky_sub") or {})
        except Exception:
            _ssub = {}
        if _ssub.get("applied"):
            _fm = r > 0.72 * float(r.max())
            if _fm.sum() > 10000:
                _fc = np.median(np.asarray(hdr[::4, ::4], np.float32)[_fm[::4, ::4]].reshape(-1, 3), axis=0)
                _fl = 0.2126 * _fc[0] + 0.7152 * _fc[1] + 0.0722 * _fc[2]
                if np.isfinite(_fc).all() and _fc.min() > 0 and _fl > 0:
                    _fade_to = np.clip(_fc / _fl, 0.2, 3.0).astype(np.float32)
            del _fm
        self.colour_fade_to = _fade_to
        ratio = np.empty(lum.shape + (3,), np.float32)
        for c in range(3):
            rc = ndimage.gaussian_filter(
                np.ascontiguousarray(hdr[:, :, c]), 6) / np.maximum(Ls, _floor)
            ratio[:, :, c] = _fade_to[c] + _conf * (rc - _fade_to[c])
            del rc
        self.colour_floor = float(_floor)
        self.colour_conf_frac = float((_conf > 0.5).mean())
        # THE CORONA AS A WHITE REFERENCE, measured here so it costs no re-run.
        #
        # The K-corona is photospheric light Thomson-scattered off free
        # electrons, and Thomson scattering is wavelength-independent. The inner
        # corona therefore carries the SUN's spectrum: it is the one thing in
        # the frame that is white by physics rather than by convention. The
        # dusty F-corona is redder and takes over further out, so the reference
        # is taken close in, at 1.05-1.6 R, and the genuine outward reddening
        # the file contains is left alone.
        #
        # WHY IT IS MEASURED FOR EVERY DATASET, not only for the ones that
        # arrived without a camera white balance. `neutralise_corona_colour`
        # has always done this in the pipeline, gated on there being no camera
        # white balance to use -- on the theory that a raw bracket keeps the
        # body's own numbers and is therefore already right. Measured on a
        # 600 mm raw bracket that HAS camera white balance, the merged corona
        # reads R/G 1.95 and B/G 0.28 at 1.05 R: twice the red and a quarter of
        # the blue. The pipeline's own prominence gate reports the same thing
        # from the other side ("corona R/GB 3.02") and has been calibrated
        # around it for as long as it has existed. So the gate was wrong, and
        # the gains belong to every dataset as a CONTROL rather than to a few as
        # an automatic fallback.
        #
        # Stored, never applied here. `coronaNeutral` in the render decides how
        # much of it to use, so the default picture is unchanged and the
        # correction is adjustable without re-stacking.
        self.corona_gain = np.ones(3, np.float32)
        self.corona_r_over_gb = None
        try:
            _cm = (r > 1.05 * self.R) & (r < 1.60 * self.R)
            if _cm.sum() >= 2000:
                _md = np.median(hdr[_cm].reshape(-1, 3),
                                axis=0).astype(np.float64)
                if np.isfinite(_md).all() and _md.min() > 0:
                    _l = 0.2126 * _md[0] + 0.7152 * _md[1] + 0.0722 * _md[2]
                    _gg = _l / _md
                    # A few times' correction is a colour cast; more than that
                    # is a broken channel or a mono frame read as colour, and
                    # dividing by it would invent colour rather than correct it.
                    if float(_gg.max() / max(_gg.min(), 1e-6)) <= 8.0:
                        self.corona_gain = _gg.astype(np.float32)
                        self.corona_r_over_gb = float(
                            _md[0] / max(0.5 * (_md[1] + _md[2]), 1e-9))
            del _cm
        except Exception:
            pass
        # Keep the confidence map. Neutralising the sky cast has to be applied
        # WITH it: `ratio` above is already faded to exactly neutral wherever
        # the signal is below the noise, so dividing that region by the sky's
        # colour a second time does not neutralise anything -- it tips an
        # already-grey sky blue. Weighted by confidence the correction is
        # consistent: full where there is real chroma to correct, absent where
        # the chroma was discarded. (0.11.4)
        self._cconf = _conf
        # colour of the sky far from the corona. At low sun altitude extinction
        # crushes blue, so the background carries a real yellow-green cast that
        # is atmosphere, not corona; dividing it out neutralises the sky while
        # leaving the corona's own (very different) colour recognisable.
        # Measured on the HDR itself, not on `ratio` (fixed in 0.11.4).
        # `ratio` carries the confidence fade above, which drives it to exactly
        # 1.0 wherever the signal is near the noise floor -- which is precisely
        # this region. Mean confidence out here measures 0.015, so the old
        # measurement returned R 1.000 G 1.000 B 1.000 on a sky whose real
        # colour is R 0.985 G 1.037 B 0.681, and the Neutralise sky cast slider
        # did nothing at any setting.
        rmx = float(r.max())
        farm = r > 0.72 * rmx
        if farm.sum() > 10000:
            # subsampled: hdr is a memmap and this region is tens of millions
            # of pixels; a median over every 4th row and column is the same
            # number for a fraction of the memory
            _fs = farm[::4, ::4]
            _hs = np.asarray(hdr[::4, ::4], np.float32)[_fs].reshape(-1, 3)
            bc = (np.median(_hs, axis=0).astype(np.float32) if _hs.shape[0] > 2000
                  else np.ones(3, np.float32))
            del _hs, _fs
            if not np.isfinite(bc).all() or bc.min() <= 0:
                bc = np.ones(3, np.float32)
        else:
            bc = np.ones(3, np.float32)
        bl = 0.2126 * bc[0] + 0.7152 * bc[1] + 0.0722 * bc[2]
        self.bg_chroma = np.clip(bc / max(float(bl), 1e-6), 0.3, 3.0).astype(np.float32)
        # No sky cast is left to neutralise once the sky was subtracted: the far
        # field is the corona's own colour (see _fade_to above). Dividing by it
        # whitened the corona where the colour is trusted and not beyond, which
        # is the same hue edge the other way round (lab 0.24).
        if bool((_ssub or {}).get("applied")):
            self.bg_chroma = np.ones(3, np.float32)
        # The prominence detail layer is optional: a workdir built before 0.18
        # does not have one, and a run whose gate found nothing does not write
        # one. When it is missing the key is ABSENT, not filled with 0.5 --
        # a flat 0.5 is NOT the identity here, because the blend pulls `det`
        # TOWARDS the layer, so a flat one would erase detail inside the gate
        # rather than do nothing. Caught by the fallback test, which is the
        # only reason this comment is not still claiming otherwise.
        _pdp = os.path.join(wd, "promdet.npy")
        self.has_promdet = os.path.exists(_pdp)
        # HILL'S CHAIN: the log-mapped base and its set of unsharp masks. Both
        # optional -- a workdir stacked before 0.22.64 has neither, and the
        # renderer falls back to the multiplicative path with hillMix ignored.
        _hp, _hl = os.path.join(wd, "hill.npy"), os.path.join(wd, "hill_log.npy")
        self.has_hill = os.path.exists(_hp) and os.path.exists(_hl)
        self.hill_scales, self.hill_rms, self.hill_rms_struct = [], [], []
        if self.has_hill:
            _hm = np.load(_hp, mmap_mode="r")
            self.n_hill = int(_hm.shape[0])
            try:
                _rep = json.load(open(os.path.join(wd, "report.json")))
                _hs = (_rep.get("hill") or {})
                self.hill_scales = [float(x) for x in _hs.get("scales", [])]
                self.hill_rms = [float(x) for x in _hs.get("rms", [])]
                self.hill_rms_struct = [float(x)
                                        for x in _hs.get("rms_struct", [])]
            except Exception:
                pass
            if len(self.hill_rms) != self.n_hill:
                # Measure it here rather than render with an unknown scale.
                # ROBUSTLY: a plain std over the whole plane is dominated by the
                # occulted disc, where the masks are the difference between a
                # black plateau and a partial convolution that has no data to
                # work with. On the reference fixture that read 3.4e-2 against a
                # true 5.7e-3 -- a factor of six, straight into the slider scale.
                #
                # Excluding the disc, which is the SAME estimator build_hill
                # records, so a slider means one thing whichever path measured
                # it. A robust spread was tried instead and reads half the value
                # on this data: the masks are peaky, most of the frame carries
                # nothing, and MAD is not std on a distribution like that.
                _keep = (r[::8, ::8] > self.Rmask)
                self.hill_rms = []
                for i in range(self.n_hill):
                    _s = np.asarray(_hm[i, ::8, ::8], np.float32)[_keep]
                    _s = _s[np.isfinite(_s)]
                    self.hill_rms.append(float(np.std(_s)) if _s.size else 1.0)
                del _keep
            if not self.hill_scales:
                from .detail import HILL_SCALES as _HS
                self.hill_scales = [float(x) for x in _HS[:self.n_hill]]
            # the k the masks were BUILT at; the render can move away from it
            self.hill_resp = []
            try:
                _hj = (json.load(open(os.path.join(
                    wd, "report.json"))).get("hill") or {})
                self.hill_logk = float(_hj.get("logK", 6.0))
                self.hill_resp = [float(x) for x in _hj.get("resp", [])]
            except Exception:
                self.hill_logk = 6.0
            # THE AMPLIFICATION SLIDER'S REFERENCE (0.23.9): the spread of the
            # COARSEST mask just outside the limb, Rmask to 1.5 R. Measured on
            # the 600 mm reference stack that shell of the 32 px mask is 99%
            # structure, while the figure used before (rms_struct[0], the 2 px
            # mask over the whole frame) was a fixed 25% of a spread made of
            # far-field noise -- the photon model reads 13-24% high at 1.5-4 R,
            # so total^2 - noise^2 went negative on every scale and every
            # scale sat on the floor. HILL_REF_C = 6.25 is that stack's ratio
            # between the two, so any slider value gives the picture it gave
            # before on the reference set; other stacks now get the same
            # amount of STRUCTURE per slider unit, not the same amount of noise.
            # Computed here from hill.npy, so no mask rebuild is needed.
            self.hill_ref = 0.0
            try:
                _top = self.n_hill - 1
                _q = 4
                _mt = np.asarray(_hm[_top, ::_q, ::_q], np.float32)
                _ry = (np.arange(_mt.shape[0], dtype=np.float32) * _q - self.cy)[:, None]
                _rx = (np.arange(_mt.shape[1], dtype=np.float32) * _q - self.cx)[None, :]
                _rq = np.sqrt(_ry * _ry + _rx * _rx)
                _sh = (_rq > self.Rmask) & (_rq < 1.5 * float(self.R)) & (_mt != 0) & np.isfinite(_mt)
                if _sh.sum() > 1000:
                    self.hill_ref = float(np.std(_mt[_sh]))
                del _mt, _ry, _rx, _rq, _sh
            except Exception:
                self.hill_ref = 0.0
            # the per-pixel sigma of im_log; absent on masks built before
            # 0.22.74, in which case the noise threshold is simply unavailable
            _hs = os.path.join(wd, "hill_sigma.npy")
            self.has_hill_sigma = (os.path.exists(_hs)
                                   and len(self.hill_resp) == self.n_hill)
        self.full = {"bg": Bg, "rprof": rprof, "mgn": mg, "fnrgf": fnl,
                     "inner": inner, "inner0": inner0, "earth": earth,
                     "prom": gate, "pel": pel, "ratio": ratio,
                     "cconf": self._cconf}
        if self.has_nafe:
            self.full["nafe"] = nfl
        if self.has_mgn_fine:
            self.full["mgn_fine"] = mgf
        if self.has_promdet:
            self.full["promdet"] = np.load(_pdp, mmap_mode="r")
        if self.has_hill:
            self.full["hilllog"] = np.load(_hl, mmap_mode="r")
            if self.has_hill_sigma:
                self.full["hillsigma"] = np.load(_hs, mmap_mode="r")
            for i in range(self.n_hill):
                self.full[f"hill{i}"] = _hm[i]
        q = preview_decim
        # THE PREVIEW IS AREA-AVERAGED, NOT SUBSAMPLED.
        #
        # This used to be `v[::q, ::q]`, which throws away 15 of every 16 pixels
        # with no prefilter. On a detail layer whose finest scale is 1.4 px that
        # is aliasing, and aliasing of fine radial threads does not look like
        # missing detail -- it looks like COARSE detail, because the energy folds
        # down into low frequencies as mottle. The preview was therefore a
        # systematically rougher picture than the export it was supposed to
        # stand for, and every judgement about "too coarse" was being made on
        # it. Averaging costs nothing and makes the preview an honest miniature.
        self.prev = {k: _decim(v, q) for k, v in self.full.items()}
        # clarity/smoothing variants (preview scale; full-res computed on demand)
        for key in (("mgn", "fnrgf") + (("nafe",) if self.has_nafe else ())
                    + (("mgn_fine",) if self.has_mgn_fine else ())):
            f = _fill_disc(self.prev[key], self.cy / q, self.cx / q, self.Rmask / q)
            self.prev[key + "_lo"] = ndimage.gaussian_filter(f, 8.0 / q * 2)
            self.prev[key + "_sm"] = ndimage.gaussian_filter(f, 0.6)
            del f
        self.prev_decim = q
        # THE 0.24 PROMINENCE LAYER: its maps on the disc box at full resolution
        # (the full-frame arrays would cost ~0.7 GB), and full-frame at preview
        # scale, area-averaged on the same q grid as everything else.
        self.pl_box, self.pl_full = None, None
        if _use_corona:
            try:
                _r = _plm.display_maps(wd, geo)
                if _r is not None:
                    self.pl_box, self.pl_full = _r
                    self.has_promlayer = True
                    y0, y1, x0, x1 = self.pl_box
                    Y0, X0 = y0 // q * q, x0 // q * q
                    Y1, X1 = min(-(-y1 // q) * q, H), min(-(-x1 // q) * q, W)
                    hp_, wp_ = self.prev["bg"].shape
                    for _k, _v in self.pl_full.items():
                        _sh = (Y1 - Y0, X1 - X0) + _v.shape[2:]
                        _tmp = np.zeros(_sh, np.float32)
                        if _k == "C":
                            _tmp[...] = 1.0
                        _tmp[y0 - Y0:y1 - Y0, x0 - X0:x1 - X0] = _v
                        _d = _decim(_tmp, q)
                        _full = np.zeros((hp_, wp_) + _v.shape[2:], np.float32)
                        if _k == "C":
                            _full[...] = 1.0
                        _ty, _tx = Y0 // q, X0 // q
                        _hh = min(_d.shape[0], hp_ - _ty); _ww = min(_d.shape[1], wp_ - _tx)
                        _full[_ty:_ty + _hh, _tx:_tx + _ww] = _d[:_hh, :_ww]
                        self.prev["pl_" + _k] = _full
                    del _tmp, _d, _full
            except Exception as _e:
                print(f"prominence layer maps not loaded ({_e})")
                self.has_promlayer = False
        # HOW MUCH OF EACH MASK THE PREVIEW'S AREA AVERAGE REMOVES.
        #
        # Averaging 4x4 is the honest way to shrink a picture, and for every
        # LINEAR term the preview is then a true miniature of the export: the
        # amplification is a weighted sum of the masks, so averaging the masks
        # and summing gives the same answer as summing and averaging.
        #
        # The NOISE THRESHOLD is not linear, and that is where the preview and
        # the export came apart. It compares |M| against k*resp*sigma. A 2 px
        # unsharp mask is zero-mean over a 4 px box, so the average knocks its
        # amplitude down to about a third -- while `hill_sigma` is a smooth
        # photon-noise map that the same average leaves untouched. Measured on
        # the smoke fixture, preview std / full-res std per scale:
        #
        #     2 px 0.29   4 px 0.45   8 px 0.70   16 px 0.90   32 px 0.96
        #
        # So at one slider setting the threshold was cutting the 2 px mask at
        # roughly 3.5x its intended level on screen and at the intended level in
        # the TIFF: the preview looked cleaner than the file it stood for, which
        # is the same class of dishonesty the area average was introduced to fix
        # and the reason the export could look grainier than what was tuned.
        #
        # The fix is NOT to subsample the masks instead. That would restore the
        # amplitude by aliasing the finest layers in the app -- the 2 px mask is
        # the worst case there is for it -- and buy an honest threshold with a
        # dishonest picture. Scale the THRESHOLD by the same factor the average
        # took out, and the preview cuts the same FRACTION of coefficients the
        # export does while staying a true miniature.
        #
        # Measured per work directory rather than assumed: the ratio depends on
        # the plate scale, which is why it is not a constant in the code.
        self.hill_prev_att = []
        if self.has_hill:
            _rr = np.hypot(
                np.arange(self.full["bg"].shape[0], dtype=np.float32)[:, None]
                - self.cy,
                np.arange(self.full["bg"].shape[1], dtype=np.float32)[None, :]
                - self.cx)
            _mf = _rr[::8, ::8] > self.Rmask          # corona only; the disc is 0
            _mp = _rr[::q, ::q][:self.prev["hill0"].shape[0],
                                :self.prev["hill0"].shape[1]] > self.Rmask
            del _rr
            for i in range(self.n_hill):
                try:
                    _sf = float(np.std(np.asarray(
                        self.full[f"hill{i}"][::8, ::8], np.float32)[_mf]))
                    _sp = float(np.std(np.asarray(
                        self.prev[f"hill{i}"], np.float32)[_mp]))
                    _a = _sp / _sf if _sf > 1e-12 else 1.0
                except Exception:
                    _a = 1.0
                self.hill_prev_att.append(
                    float(_a) if np.isfinite(_a) and 0 < _a <= 1.0 else 1.0)
            del _mf, _mp
        _geom_for_fill["g"] = (self.cy, self.cx, self.Rmask)   # full-res variants
        self.flat_range = None
        self.flat_error = None
        self._load_flat(geo)
        self.reload_contact()

    def _load_flat(self, geo):
        """The master flat as a QC view.

        It is the one cached product that is neither aligned nor trimmed, so it
        has to be cut down to the layer grid before it can be shown beside the
        others -- crop_origin says where that grid sits in the sensor frame.
        Displayed on its own p0.5-p99.5, because a 16% falloff shown linearly
        over 0..1 is invisible, which is exactly the kind of thing this view
        exists to catch.

        The layer grid is FULL sensor resolution (the merge demosaics to it),
        while the master flat is a Bayer array on that same grid. The first
        version of this reduced the flat to superpixels -- which does kill the
        Bayer checkerboard, but leaves it at half the layer grid, so the crop
        could never match and the view silently never appeared. Reduce to
        superpixels and then put it BACK on the full grid: the checkerboard is
        gone and the geometry lines up with everything else.
        """
        p = os.path.join(self.wd, "masterflat.npy")
        if not os.path.exists(p):
            return False
        try:
            m = np.load(p).astype(np.float32)
            h2, w2 = m.shape[0] // 2, m.shape[1] // 2
            sp = m[:2 * h2, :2 * w2].reshape(h2, 2, w2, 2).mean(axis=(1, 3))
            del m
            H, W = self.full["bg"].shape
            # crop_origin is in full-res px and even (Bayer phase), so it maps
            # exactly onto the superpixel grid; crop there, then expand, so the
            # full-res copy that gets made is only the size of the view.
            oy, ox = geo.get("crop_origin") or (0, 0)
            oy, ox = int(oy) // 2, int(ox) // 2
            nh, nw = (H + 1) // 2, (W + 1) // 2
            if sp.shape[0] < oy + nh or sp.shape[1] < ox + nw:   # old cache/odd size
                oy = max((sp.shape[0] - nh) // 2, 0)
                ox = max((sp.shape[1] - nw) // 2, 0)
            sp = sp[oy:oy + nh, ox:ox + nw]
            if sp.shape != (nh, nw):
                self.flat_error = (f"master flat {sp.shape} does not cover the "
                                   f"{nh}x{nw} view grid")
                return False
            f = np.repeat(np.repeat(sp, 2, axis=0), 2, axis=1)[:H, :W]
            del sp
            lo, hi = np.percentile(f, [0.5, 99.5])
            self.flat_range = (float(f.min()), float(f.max()), float(lo), float(hi))
            f = np.clip((f - lo) / max(hi - lo, 1e-6), 0, 1).astype(np.float32)
            self.full["flat"] = f
            self.prev["flat"] = f[::self.prev_decim, ::self.prev_decim]
            return True
        except Exception as e:
            # Never take the render down for a QC view -- but say so, because a
            # silently missing button is what hid the bug above for a release.
            self.flat_error = f"{type(e).__name__}: {e}"
            sys.stderr.write(f"master-flat preview unavailable ({e})\n")
            return False

    @property
    def has_flat(self):
        return "flat" in self.full

    def reload_contact(self):
        p = os.path.join(self.wd, "contact_rgb.npy")
        if os.path.exists(p):
            c = np.load(p).astype(np.float32)
            self.full["contact"] = c
            self.prev["contact"] = c[::self.prev_decim, ::self.prev_decim]
            return True
        return False

    @property
    def has_contact(self):
        return "contact" in self.full

    def geometry(self, decim=1):
        return self.cy / decim, self.cx / decim, self.R / decim

    def mask_radius(self, decim=1):
        return self.Rmask / decim

    def mask_radius_map(self, shape, decim=1):
        """Per-azimuth mask radius on the given grid; scalar if no profile."""
        if not self.limb_prof:
            return self.Rmask / decim
        from .detail import limb_radius_map
        return limb_radius_map(np.asarray(self.limb_prof, np.float32) / decim,
                               shape, self.cy / decim, self.cx / decim,
                               self.limb_margin / decim)


_geom_for_fill = {}


def _variant(src, key, kind, preview):
    k = key + "_" + kind
    if k in src:
        return src[k]
    sigma = 16.0 if kind == "lo" else 0.6 * 4   # full-res equivalents
    g = _geom_for_fill.get("g")
    if g is not None:
        return ndimage.gaussian_filter(_fill_disc(src[key], *g), sigma)
    return ndimage.gaussian_filter(src[key], sigma)


# THE MGN VIEW shows the MGN layer itself, stretched for inspection (2026-10-07):
# it used to show the composite's MGN term -- the layer times MGN contrast
# (0.15 by default) with the NAFE mix blended in -- which left about one display
# level of MGN in the outer corona and read as a posterised, blotchy layer.
VIEW_MGN_STRETCH = 2.0


def _detail_layers(src, P, preview=True, views=False):
    """Runtime-transformed detail layers (shared by composite render and layer views)."""
    # Detail balance: blend the all-scale layer with the fine-only one. Because
    # an MGN layer is a weighted mean of its per-scale terms, this IS a ladder
    # change, not an approximation of one -- measured fine/coarse ratio 1.09 at
    # 0, 3.18 at 1, against 3.32 for Hill's true gains.
    _ds = float(P.get("detailScale", 0.0)) if "mgn_fine" in src else 0.0

    def _mgsrc(suffix=None):
        a = src["mgn"] if suffix is None else _variant(src, "mgn", suffix, preview)
        if _ds <= 0:
            return a
        b = src["mgn_fine"] if suffix is None else _variant(src, "mgn_fine", suffix, preview)
        return (1.0 - _ds) * a + _ds * b

    mgr = _mgsrc()
    if P["clarity"] > 0:
        mgr = mgr + P["clarity"] * (mgr - _mgsrc("lo"))
    if P["smoothing"] > 0:
        mgr = (1 - P["smoothing"]) * mgr + P["smoothing"] * _mgsrc("sm")
    mg = np.clip(0.5 + (mgr - 0.5) * P["mgnContrast"], 0, 1)
    fnr = src["fnrgf"]
    if P["clarity"] > 0:
        fnr = fnr + P["clarity"] * (fnr - _variant(src, "fnrgf", "lo", preview))
    if P["smoothing"] > 0:
        fnr = (1 - P["smoothing"]) * fnr + P["smoothing"] * _variant(src, "fnrgf", "sm", preview)
    D = (fnr - 0.5) * 8.0                     # back to sigma units
    s_hi = 2.5 / max(P["fnCompress"], 1e-4)   # 0 = FNRGF off (flat 0.5)
    fn = (np.where(D >= 0, np.tanh(D / s_hi), np.tanh(D / (1.5 * s_hi))) + 1) / 2
    fn0 = fn
    # NAFE-VN rides with the other two rather than replacing them: it sees the
    # faint outer structure they flatten away, and because it needs no disc
    # geometry it stays clean at the limb where they are most fragile.
    #
    # The stored layer is E, the equalized field -- not the paper's eq. 2
    # output B = (1-w) T_gamma(A) + w E. That combination happens HERE and one
    # level up: the composite's envelope plays the role of T_gamma, and
    # nafeMix is w. Their w runs 0.05..0.3, so the useful part of this slider
    # is the bottom third; past that the rank field starts to overwhelm the
    # envelope's own falloff.
    nf_ = src.get("nafe")
    if NAFE_ENABLED and nf_ is not None and P.get("nafeMix", 0.0) > 0:
        a = float(np.clip(P["nafeMix"], 0.0, 1.0))
        nfd = nf_
        if P["clarity"] > 0:
            nfd = nfd + P["clarity"] * (nfd - _variant(src, "nafe", "lo", preview))
        if P["smoothing"] > 0:
            nfd = (1 - P["smoothing"]) * nfd + P["smoothing"] * _variant(src, "nafe", "sm", preview)
        nfd = np.clip(0.5 + (nfd - float(np.median(nfd))) * 1.0, 0, 1)
        mg = np.clip((1 - a) * mg + a * nfd, 0, 1)
        fn = np.clip((1 - a) * fn + a * nfd, 0, 1)
    inner_eff = (1 - P["innerDenoise"]) * src["inner0"] + P["innerDenoise"] * src["inner"]
    if views:
        return mg, fn, inner_eff, mgr, fn0
    return mg, fn, inner_eff


# FINAL SHARPENING (lab, Oct 2026). Nico's Photoshop finish, read out of the
# action file "Nimmervoll's Aktionen V1.2" and reproduced here, plus a noise
# gate Photoshop cannot have.
#
# The action: five copies of the finished image, High Pass at 1.5 / 3 / 6 /
# 12 / 24 px, all in Soft Light; the five stacked into a smart object and
# combined three ways -- the MEDIAN across the scales ("Mittel Toene", 80 %),
# the MINIMUM ("Tiefen", 60 %: only the darkening excursions) and the MAXIMUM
# ("Lichter", 40 %: only the brightening ones) -- grouped and pasted back.
# Soft light of a 0.5-grey high pass is unsharp masking that eases off in the
# shadows and highlights; the min/max split lets the dark and the bright side
# of every edge be weighted separately.
#
# What is different here, and why:
#   * The radii follow the solar radius (1.5 px at R = 620, the 600 mm
#     reference bracket), so the same angular scales are sharpened on any
#     camera. Hill's ladder is fixed in pixels; this one is not.
#   * A noise gate. Measured on the reference set at full resolution the
#     action raises the 1 px band by x1.7, the 2 px band by x1.5 and the
#     8-16 px structure by x1.1-1.2: most of what it adds is grain. The
#     gate takes each scale's high pass through erf(|d| / (k sigma_s)) --
#     the WOW attenuation (Auchere et al. 2023) -- where sigma_s is the
#     expected noise amplitude of that high pass. sigma_s comes from the
#     picture itself: the robust spread of the finest band in 32 px tiles
#     (the finest band of a finished render is noise -- the set resolves to
#     1.4 px), carried to the other scales by each filter's measured white-
#     noise response. No model, no camera table; k is the "Noise gate"
#     slider, 0 switches it off and gives exactly the Photoshop result.
#   * Applied to the FINISHED display image, after tone and colour, before
#     the ring overlay -- where the action runs in Photoshop.
#   * Row bands with a margin of three times the widest radius, so a 44 Mpx
#     export does not hold five full-size high-pass planes at once.


# ---- block 3a: tone curve, contained whites, zoned Hill gains ---------------
def tone_lut_pts(pts, n=1024):
    """Monotone cubic (Fritsch-Carlson) through the user's points [[in, out],
    ...] on 0..1, as a LUT. Points are sorted by `in`; duplicate or reversed
    `out` values are clamped so the curve never turns back on itself."""
    P = sorted((float(x), float(y)) for x, y in pts)
    if not P or P[0][0] > 1e-6:
        P.insert(0, (0.0, 0.0))
    if P[-1][0] < 1 - 1e-6:
        P.append((1.0, 1.0))
    xs = np.array([p[0] for p in P]); ys = np.array([p[1] for p in P])
    keep = np.concatenate([[True], np.diff(xs) > 1e-4]); xs, ys = xs[keep], ys[keep]
    ys = np.maximum.accumulate(np.clip(ys, 0, 1))
    k = xs.size
    t = np.linspace(0, 1, n)
    if k < 2:
        return t.astype(np.float32)
    h = np.diff(xs); d = np.diff(ys) / h
    m = np.zeros(k); m[0] = d[0]; m[-1] = d[-1]
    for i in range(1, k - 1):
        m[i] = 0.0 if d[i - 1] * d[i] <= 0 else 2.0 / (1.0 / d[i - 1] + 1.0 / d[i])
    out = np.empty(n)
    for i in range(k - 1):
        sel = (t >= xs[i]) & (t <= xs[i + 1]); u = (t[sel] - xs[i]) / h[i]
        h00 = 2 * u ** 3 - 3 * u ** 2 + 1; h10 = u ** 3 - 2 * u ** 2 + u
        h01 = -2 * u ** 3 + 3 * u ** 2; h11 = u ** 3 - u ** 2
        out[sel] = h00 * ys[i] + h10 * h[i] * m[i] + h01 * ys[i + 1] + h11 * h[i] * m[i + 1]
    return np.clip(out, 0, 1).astype(np.float32)


def curve_is_identity(pts):
    try:
        return all(abs(float(x) - float(y)) < 1e-6 for x, y in pts)
    except Exception:
        return True


def _lut_apply(rgb, pts):
    lut = tone_lut_pts(pts)
    idx = np.clip(rgb * (lut.size - 1), 0, lut.size - 1)
    lo = np.floor(idx).astype(np.int32); fr = (idx - lo).astype(np.float32)
    hi_i = np.minimum(lo + 1, lut.size - 1)
    return (lut[lo] * (1 - fr) + lut[hi_i] * fr).astype(np.float32)


def apply_tone(rgb, pts, ceiling, moon_pts=None, eg=None):
    """The tone curve (user points) on every channel, then 'contain whites':
    the picture is scaled down so the 99.99th percentile of its luminance
    sits at `ceiling` (1.0 = off). Never scales up.

    THE MOON HAS ITS OWN CURVE (lab 0.24, Nico): the disc sits low on the
    display scale (Moon disc level ~0.1), exactly where a mid-tone lift bites
    hardest, so the main curve moved the Moon a lot and the corona hardly at
    all. With `eg` (1 outside the Moon, 0 on the disc, the disc's soft edge in
    between) the main curve acts outside the Moon and `moon_pts` on it."""
    g_on = bool(pts) and not curve_is_identity(pts)
    m_on = eg is not None and bool(moon_pts) and not curve_is_identity(moon_pts)
    if eg is None:
        if g_on:
            rgb = _lut_apply(rgb, pts)
    elif g_on or m_on:
        G = _lut_apply(rgb, pts) if g_on else rgb
        M = _lut_apply(rgb, moon_pts) if m_on else rgb
        e = np.asarray(eg, np.float32)[:, :, None]
        rgb = (G * e + M * (1 - e)).astype(np.float32)
        del G, M, e
    if ceiling < 0.999:
        Y = 0.2126 * rgb[:, :, 0] + 0.7152 * rgb[:, :, 1] + 0.0722 * rgb[:, :, 2]
        p = float(np.percentile(Y[::4, ::4], 99.99))
        if p > ceiling:
            rgb = (rgb * np.float32(ceiling / p)).astype(np.float32)
    return rgb


# THE LOOK CONTROLS (lab 0.24). Four controls that act on the finished picture
# and set what a viewer reads as the "look" of a published eclipse, measured off
# Druckmüller, Lefaudeux and Project Helion renders (looks/refs.json):
#   coronaLevel  the median display brightness of the inner corona, 1.05-1.15 R.
#                Helion and Lefaudeux sit at 0.53-0.69, i.e. the corona is NOT
#                white there; ours was 0.9+. 0 = off (picture left as it is).
#   skyLevel     a screen lift of everything outside the Moon towards the sky
#                colour: the sky becomes ~skyLevel, the white corona is untouched.
#   skyColour    the lift's hue on one axis: 0 neutral grey, 1 Helion's sky blue
#                (R/G 0.766, B/G 1.323), 2 Druckmüller's deeper blue.
#   moonColour   the Moon disc's hue on the same axis (0 = neutral grey).
# All four off by default, so every render before them is unchanged.
LOOK_BLUE = np.array([0.766, 1.0, 1.323])


def look_chroma(t, base=None):
    b = LOOK_BLUE if base is None else np.asarray(base, np.float64)
    c = b ** float(t)
    return (c / (0.2126 * c[0] + 0.7152 * c[1] + 0.0722 * c[2])).astype(np.float32)


def sky_base(P):
    """The sky lift's direction from the picked colour (display RGB), green as
    the reference channel; the default pick is LOOK_BLUE / 1.323, the old axis."""
    h = np.array([float(P.get("skyHueR", 0.579)), float(P.get("skyHueG", 0.756)),
                  float(P.get("skyHueB", 1.0))], np.float64)
    h = np.clip(h, 0.02, 1.0)
    return h / h[1]


def moon_base(P):
    """The Moon tint's direction from the picked colour (display RGB), with
    green as the reference channel as LOOK_BLUE has it -- so the default pick
    (LOOK_BLUE / 1.323) gives exactly the old blue axis."""
    h = np.array([float(P.get("moonHueR", 0.579)), float(P.get("moonHueG", 0.756)),
                  float(P.get("moonHueB", 1.0))], np.float64)
    h = np.clip(h, 0.02, 1.0)
    return h / h[1]


def corona_white(layers, P, lo=1.05, hi=1.6):
    """Temperature and Tint that render the 1.05-1.6 R ring neutral -- the
    page's solveWhite over that ring, same chain (ratio, sky cast, corona gain)."""
    src = layers.prev if layers.prev is not None else layers.full
    decim = layers.prev_decim if layers.prev is not None else 1
    cy, cx, R = layers.geometry(decim)
    H, W = src["bg"].shape
    yy = np.arange(H, dtype=np.float32)[:, None] - cy
    xx = np.arange(W, dtype=np.float32)[None, :] - cx
    rr = np.sqrt(yy * yy + xx * xx)
    Rm = max(layers.mask_radius(decim), lo * R)
    m = (rr >= Rm) & (rr < hi * R)
    a = np.clip(np.asarray(src["ratio"], np.float32)[m], 0.2, 3.0)
    if P.get("bgNeutral", 0) > 0:
        a = a / (1.0 + np.asarray(src["cconf"], np.float32)[m][:, None]
                 * (layers.bg_chroma[None, :] ** P["bgNeutral"] - 1.0))
    if P.get("coronaNeutral", 0) > 0:
        a = a * layers.corona_gain[None, :] ** P["coronaNeutral"]
    m0, m1, m2 = (float(np.median(a[:, k])) for k in range(3))
    temp = float(np.sqrt(m2 / max(m0, 1e-6)))
    return temp, float(m0 * temp / max(m1, 1e-6))


def look(rgb, P, r, R, eg):
    cl = float(P.get("coronaLevel", 0) or 0)
    sl = float(P.get("skyLevel", 0) or 0)
    mc = float(P.get("moonColour", 0) or 0)
    if cl <= 0 and sl <= 0 and mc == 0:
        return rgb
    rgb = np.asarray(rgb, np.float32)
    if cl > 0:
        m = (r >= 1.05 * R) & (r < 1.15 * R) & (eg > 0.999)
        if m.sum() > 50:
            v = rgb[m]
            med = float(np.median(0.2126 * v[:, 0] + 0.7152 * v[:, 1] + 0.0722 * v[:, 2]))
            if med > 1e-4:
                k = np.float32(cl / med)
                # outside the Moon only: the disc keeps the level Moon disc level gives it
                rgb = rgb * (1 + (k - 1) * eg)[:, :, None]
                _hc = float(P.get("hlCompress", 0.0) or 0.0)
                if k > 1 or _hc > 0:     # hue-preserving shoulder, as the main knee
                    # HIGHLIGHT COMPRESSION ACTS HERE (2026-10-08). The main
                    # knee sits before this level step, on a picture that is
                    # still too dark to reach it; the step then lifted
                    # everything back, so the slider did nothing in any look
                    # that sets an inner corona level (all of them). The knee
                    # is now this shoulder's: 0.9 at 0 (unchanged), 0.55 at 1.
                    kn = np.float32(0.9 - 0.35 * _hc)
                    mx = rgb.max(axis=2)
                    ms = np.where(mx <= kn, mx, kn + (1 - kn) * np.tanh((mx - kn) / (1 - kn)))
                    rgb *= np.where(mx > 1e-6, ms / np.maximum(mx, 1e-6), 1.0)[:, :, None]
        # CONTAIN WHITES after the level step too, for the same reason:
        # applied before it, the step scaled it straight back out.
        cw = float(P.get("containWhites", 1.0) or 1.0)
        if cw < 0.999:
            Yc = 0.2126 * rgb[:, :, 0] + 0.7152 * rgb[:, :, 1] + 0.0722 * rgb[:, :, 2]
            pc = float(np.percentile(Yc[::4, ::4], 99.99))
            del Yc
            if pc > cw:
                rgb = (rgb * np.float32(cw / pc)).astype(np.float32)
    if sl > 0:
        lift = np.clip(sl * look_chroma(P.get("skyColour", 0) or 0, sky_base(P))[None, None, :]
                       * eg[:, :, None], 0, 0.95)
        rgb = 1 - (1 - np.clip(rgb, 0, 1)) * (1 - lift)
    if mc != 0:
        # the tint is the DISC's own level times the hue shift, so the bright limb
        # inside the edge ramp is not multiplied into a coloured ring
        dm = eg <= 0
        if dm.sum() > 20:
            v = rgb[dm]
            yd = float(np.median(0.2126 * v[:, 0] + 0.7152 * v[:, 1] + 0.0722 * v[:, 2]))
            # per pixel (the maria keep their tint in proportion), capped at
            # twice the disc's level so the bright limb in the edge ramp is not
            # multiplied into a coloured ring
            ypx = np.minimum(0.2126 * rgb[:, :, 0] + 0.7152 * rgb[:, :, 1] + 0.0722 * rgb[:, :, 2], 2.0 * yd)
            rgb = rgb + ((1 - eg) * ypx)[:, :, None] * (look_chroma(mc, moon_base(P))[None, None, :] - 1)
    return np.clip(rgb, 0, 1).astype(np.float32)


# Zoned Hill gains: a radial emphasis per mask scale, fine scales inside the
# zone radius, coarse ones outside, with a 0.4 R feather. zone 0 = as before.
HILL_ZONE_SIGN = (1.0, 1.0, 0.0, -1.0, -1.0)      # 2, 4, 8, 16, 32 px
HILL_ZONE_FEATHER_R = 0.4


def hill_zone_t(r, R, zone_r):
    t = np.clip((r - zone_r * R) / (HILL_ZONE_FEATHER_R * R) + 0.5, 0, 1)
    return (t * t * (3 - 2 * t)).astype(np.float32)        # 0 inside, 1 outside


FINAL_RADII = (1.5, 3.0, 6.0, 12.0, 24.0)     # px at R = 620
FINAL_REF_R = 620.0
# the gate per scale, relative to the slider: the finest scale is nearly all
# noise and gets the full k, the coarse ones carry structure and little noise
# (WOW uses n_s = 5, 3, 1 on its first three scales and nothing beyond)
FINAL_GATE = (1.0, 0.6, 0.35, 0.2, 0.1)


def _soft_light(a, b):
    """Photoshop soft light, a = base, b = blend, both in 0..1"""
    d = np.where(a <= 0.25, ((16 * a - 12) * a + 4) * a, np.sqrt(np.clip(a, 0, None)))
    return np.where(b <= 0.5, a - (1 - 2 * b) * a * (1 - a), a + (2 * b - 1) * (d - a))


def _noise_response(radii):
    """rms of each BAND (finest: n - G_r0 n; then G_r(i-1) n - G_ri n) for unit
    white noise n -- the factor that carries a sigma measured on the finest
    band to the other bands"""
    rng = np.random.default_rng(7)
    n = rng.standard_normal((256, 256)).astype(np.float32)
    g = [ndimage.gaussian_filter(n, s) for s in radii]
    out = [float((n - g[0]).std())]
    for i in range(1, len(radii)):
        out.append(float((g[i - 1] - g[i]).std()))
    return np.array(out)


def final_sharpen(rgb, R, mid=0.0, dark=0.0, light=0.0, k=2.0, tile=32, band_rows=768):
    """rgb: finished display image, float32 in 0..1, (H, W, 3). R: solar radius
    on THIS grid (already divided by the preview decimation).

    The five Photoshop high passes are rebuilt from BANDS: b0 = I - G_1.5 I,
    b_i = G_r(i-1) I - G_ri I, and HP_i = b0 + ... + b_i. With the gate off the
    sums are exactly the high passes. With it on, each band is attenuated
    against its own noise level, so the pixel noise that every high pass
    carries (the 24 px one as much as the 1.5 px one) is taken out at the
    source instead of being let through wherever a coarse plane is strong."""
    if max(mid, dark, light) <= 1e-6:
        return rgb
    H, W = rgb.shape[:2]
    radii = [max(0.6, s * (float(R) / FINAL_REF_R)) for s in FINAL_RADII]
    resp = _noise_response(radii); resp = resp / max(resp[0], 1e-9)
    Y = (0.2126 * rgb[:, :, 0] + 0.7152 * rgb[:, :, 1] + 0.0722 * rgb[:, :, 2]).astype(np.float32)
    # the noise amplitude of the finest band, per tile, smoothed, back on the grid
    sig = None
    if k > 0:
        b1 = Y - ndimage.gaussian_filter(Y, radii[0])
        ty, tx = max(H // tile, 1), max(W // tile, 1)
        bt = b1[:ty * tile, :tx * tile].reshape(ty, tile, tx, tile)
        s = 1.4826 * np.median(np.abs(bt), axis=(1, 3)).astype(np.float32)
        del b1, bt
        s = ndimage.gaussian_filter(s, 1.0)
        sig = ndimage.zoom(s, (H / s.shape[0], W / s.shape[1]), order=1).astype(np.float32)
        sig = sig[:H, :W]
        if sig.shape != (H, W):
            sig = np.pad(sig, ((0, H - sig.shape[0]), (0, W - sig.shape[1])), mode="edge")
    margin = int(np.ceil(3 * radii[-1])) + 2
    out = np.empty_like(rgb)
    for y0 in range(0, H, band_rows):
        y1 = min(H, y0 + band_rows)
        a0, a1 = max(0, y0 - margin), min(H, y1 + margin)
        blk = rgb[a0:a1]; Yb = Y[a0:a1]
        gates = None
        if sig is not None:
            gates = []
            prev = Yb
            for i, s in enumerate(radii):
                cur = ndimage.gaussian_filter(Yb, s)
                gates.append(erf(np.abs(prev - cur) / np.maximum(k * FINAL_GATE[i] * resp[i] * sig[a0:a1], 1e-6)).astype(np.float32))
                prev = cur
            del prev, cur
        hp = np.empty((len(radii),) + blk.shape, np.float32)
        for c in range(3):
            prev = blk[:, :, c]; acc = np.zeros(prev.shape, np.float32)
            for i, s in enumerate(radii):
                cur = ndimage.gaussian_filter(blk[:, :, c], s)
                band = prev - cur
                if gates is not None:
                    band = band * gates[i]
                acc = acc + band
                hp[i, :, :, c] = acc + 0.5
                prev = cur
        med = np.median(hp, axis=0); mn = hp.min(axis=0); mx = hp.max(axis=0)
        del hp, gates
        o = blk.astype(np.float32)
        if mid > 0:
            o = o + mid * (_soft_light(o, np.clip(med, 0, 1)) - o)
        if dark > 0:
            o = o + dark * (_soft_light(o, np.clip(mn, 0, 1)) - o)
        if light > 0:
            o = o + light * (_soft_light(o, np.clip(mx, 0, 1)) - o)
        out[y0:y1] = np.clip(o[y0 - a0:y1 - a0], 0, 1)
        del med, mn, mx, o
    return out


# STRUCTURE (2026-10-08, Nico: "Clarity does nothing, a Structure slider would
# be cool"). Local contrast on the FINISHED picture's brightness, between
# 0.01 R and 0.08 R (6-50 px on the 600 mm set) -- the width of the streamers
# and the gaps between them. Clarity works on the MGN/FNRGF layers only, which
# the looks mix in at 0.1-0.2, so it barely shows; this acts on everything.
#  * band = NC(ln Y, 0.01 R) - NC(ln Y, 0.08 R), normalized convolution with
#    the Moon left out, so the disc's black does not leak into the averages;
#  * the band's mean per 1 px ring is taken out, so the corona's steep radial
#    falloff is not turned into a bright or dark ring at the limb;
#  * weight: midtones (4Y(1-Y)), 0 on the Moon, ramping in over 1.01-1.04 R;
#  * applied as a brightness ratio, so the colour is kept.
# 0 = off and the picture is untouched. Page (structureBoost) and export share it.
STRUCT_R = (0.01, 0.08)
STRUCT_K = 1.5
STRUCT_LIMB = (1.01, 1.04)


def structure_boost(rgb, cy, cx, R, s):
    if s <= 1e-6:
        return rgb
    H, W = rgb.shape[:2]
    yy = np.arange(H, dtype=np.float32)[:, None] - np.float32(cy)
    xx = np.arange(W, dtype=np.float32)[None, :] - np.float32(cx)
    rp = np.sqrt(yy * yy + xx * xx)
    del yy, xx
    Y = (0.2126 * rgb[:, :, 0] + 0.7152 * rgb[:, :, 1] + 0.0722 * rgb[:, :, 2]).astype(np.float32)
    l = np.log(np.maximum(Y, 1e-4)).astype(np.float32)
    m = (rp > STRUCT_LIMB[0] * R).astype(np.float32)
    lm = l * m
    del l
    s1, s2 = max(STRUCT_R[0] * R, 0.5), max(STRUCT_R[1] * R, 1.0)
    band = (ndimage.gaussian_filter(lm, s1) / np.maximum(ndimage.gaussian_filter(m, s1), 1e-4)).astype(np.float32)
    # the wide blur on a block-averaged copy (sigma 49 px at full size), then
    # back up: the same field to ~1e-3, a fraction of the time and memory
    q = int(s2 // 6) if s2 > 12 else 1
    if q > 1:
        H2, W2 = -(-H // q), -(-W // q)
        def _dn(a):
            p = np.pad(a, ((0, H2 * q - H), (0, W2 * q - W)), mode="edge")
            return p.reshape(H2, q, W2, q).mean((1, 3)).astype(np.float32)
        lo = ndimage.gaussian_filter(_dn(lm), s2 / q) / np.maximum(ndimage.gaussian_filter(_dn(m), s2 / q), 1e-4)
        yi = (np.arange(H, dtype=np.float32) + 0.5) / q - 0.5
        xi = (np.arange(W, dtype=np.float32) + 0.5) / q - 0.5
        for y0 in range(0, H, 1024):
            y1 = min(H, y0 + 1024)
            YY, XX = np.meshgrid(yi[y0:y1], xi, indexing="ij")
            band[y0:y1] -= ndimage.map_coordinates(lo, [YY, XX], order=1, mode="nearest")
        del lo, YY, XX
    else:
        band -= (ndimage.gaussian_filter(lm, s2) / np.maximum(ndimage.gaussian_filter(m, s2), 1e-4))
    del lm
    rb = rp.astype(np.int32)
    sel = m > 0
    n_ = np.bincount(rb[sel], minlength=int(rb.max()) + 1)
    s_ = np.bincount(rb[sel], weights=band[sel].astype(np.float64), minlength=int(rb.max()) + 1)
    prof = np.where(n_ > 0, s_ / np.maximum(n_, 1), 0.0).astype(np.float32)
    band -= prof[rb]
    del rb, sel
    w = np.clip(4.0 * Y * (1.0 - Y), 0.0, 1.0) * np.clip(
        (rp / R - STRUCT_LIMB[0]) / (STRUCT_LIMB[1] - STRUCT_LIMB[0]), 0.0, 1.0)
    g = np.exp(np.float32(s * STRUCT_K) * band * w).astype(np.float32)
    return np.clip(rgb * g[:, :, None], 0, 1).astype(np.float32)


def colour_mode(P, params=None):
    """Apply the Colour mode to a parameter set (in place, returned). Monochrome
    zeroes every colour control; mono + tints zeroes the corona's own colour and
    keeps the sky, Moon and prominence tints. Mode 0 changes nothing."""
    P.update(params or {})
    m = int(P.get("colourMode", 0) or 0)
    if m == 1:
        P.update(satur=0.0, skyColour=0.0, moonColour=0.0, promSat=0.0)
    elif m == 2:
        P["satur"] = 0.0
    return P


def render(layers: Layers, params, preview=False, view="composite"):
    P = colour_mode(dict(defaults_for(getattr(layers, "mode", None))), params)
    if getattr(layers, "has_promlayer", False):
        # THE 0.24 PROMINENCE LAYER: the corona no longer holds the prominences,
        # so the gate's three terms would only act on the corona where they were
        P["promGain"] = 0.0; P["promDetail"] = 0.0; P["promChroma"] = 0.0
    src = layers.prev if preview else layers.full
    decim = layers.prev_decim if preview else 1
    cy, cx, R = layers.geometry(decim)
    H, W = src["bg"].shape
    if view in ("mgn", "fnrgf"):
        _, _, _, mgr, fn0 = _detail_layers(src, P, preview=preview, views=True)
        v = np.clip(0.5 + (mgr - 0.5) * VIEW_MGN_STRETCH, 0, 1) if view == "mgn" else fn0
        return np.repeat(np.asarray(v, np.float32)[:, :, None], 3, axis=2)
    mg, fn, inner_eff = _detail_layers(src, P, preview=preview)
    if view == "inner":
        return np.repeat(inner_eff[:, :, None], 3, axis=2)
    if view == "prom":
        return np.repeat(src["prom"][:, :, None], 3, axis=2)
    if view == "tangential":
        return np.repeat(src["pel"][:, :, None], 3, axis=2)
    if view == "nafe":
        return np.repeat(np.clip(src.get("nafe", np.full_like(mg, 0.5)), 0, 1)[:, :, None], 3, axis=2)
    if view == "flat":
        f = src.get("flat")
        if f is None:
            return np.full((H, W, 3), 0.5, np.float32)
        return np.repeat(np.clip(f, 0, 1)[:, :, None], 3, axis=2)

    yy = np.arange(H, dtype=np.float32)[:, None] - cy
    xx = np.arange(W, dtype=np.float32)[None, :] - cx
    r = np.sqrt(yy * yy + xx * xx)
    Re = layers.mask_radius_map((H, W), decim) + P["discTrim"] / decim
    # the edge ramp ends 2 px outside the mask radius; its width is a slider
    _dsw = max(float(P.get("discSoft", 12.0)), 1.0) / decim
    edge = np.clip((r - (Re + 2.0 / decim - _dsw)) / _dsw, 0, 1)
    def _ss(x):
        x = np.clip(x, 0, 1)
        return x * x * (3 - 2 * x)
    wf = _ss((r - 1.02 * R) / (0.55 * R))
    wI = edge * _ss((1.45 * R - r) / (0.40 * R))
    # Glare dim gets its OWN profile. It used to share wI, whose smoothstep
    # window closes at 1.45 R -- and at full strength that is a 3.3x brightness
    # ramp ending at a definite radius, which prints as a ring. Instrumental
    # glare does not end at a radius; it is a broad wing off the limb. An
    # exponential decay from the mask edge has no boundary to see: it is
    # monotone, never reaches zero, and its log-slope changes smoothly
    # everywhere. Scale length 0.6 R puts it at 0.47 where the old window shut
    # off and 0.11 by 3 R.
    wG = edge * np.exp(-np.maximum(r - Re, 0.0) / (0.60 * R))

    # THE ENVELOPE CURVE: cube root, or Hill's log stretch. `bg` holds
    # xn**(1/3), so xn comes back by cubing it -- exactly here, and to 0.005% in
    # the browser, which is why that layer is sent at 16 bits.
    #
    # Hill uses log(1 + k*I) with k about a million, and says it is what makes
    # the inner and outer corona visible at once. The reason it works is that a
    # log gives every factor of two the same display range wherever it sits.
    # Measured on the reference bracket, output units per doubling:
    #
    #     radius        1.1 R    2 R     3 R     4 R     6 R
    #     cube root     0.197   0.057   0.044   0.041   0.039
    #     log, k=1e6    0.050   0.050   0.050   0.050   0.050
    #
    # OFF BY DEFAULT and exactly the old curve at 0. The control is log10(k);
    # Hill's own value is 6. Turn radial flatten down as this goes up -- they do
    # the same job from opposite ends and both at once inverts the picture.
    # REMOVED 2026-10-08 (Nico): at Hill's 6 his picture went white and no look
    # used it; old settings that carry logK render as 0.
    _lk = 0.0
    if _lk > 0:
        _K = np.float32(10.0 ** _lk)
        _env = (np.log1p(_K * np.clip(src["bg"], 0, 1) ** 3)
                / np.log1p(_K)).astype(np.float32)
    else:
        _env = src["bg"]
    B = (_env ** P["envGamma"]) * (1 - P["innerDim"] * wG)
    del _env
    B = B / (src["rprof"] ** P["radialFlatten"])
    B *= (1 - P["radialFlatten"] * 0.5)       # keep overall level roughly stable
    det = (1 - P["fnMix"] * wf) * mg + P["fnMix"] * wf * fn
    det = det + P["pelGain"] * (src["pel"] - 0.5)
    det = det * (1 - P["innerMix"] * wI) + P["innerMix"] * wI * inner_eff
    # PROMINENCE DETAIL, inside the gate and nowhere else.
    #
    # The corona filters cannot carry it: MGN's normalisation window clips the
    # prominence (37% of the reference one sat hard against xn = 1.0) and its
    # local-sigma division is precisely what flattens a compact bright feature.
    # Measured as correlation with the red channel's own fine structure, which
    # normalisation cannot fake: the prominence layer scores 0.940 against
    # inner 0.317, MGN 0.274 and the merged HDR 0.262.
    #
    # `prom` is the same gate that decides where prominences are brightened, so
    # this can only act where the H-alpha colour test already fired -- it cannot
    # reach the corona at any setting.
    if P.get("promDetail", 0) > 0 and "promdet" in src:
        _pw = P["promDetail"] * np.clip(src["prom"], 0, 1)
        det = det * (1 - _pw) + _pw * src["promdet"]
        del _pw
    Y = B * (P["baseLift"] + P["detailGain"] * det)
    # HILL'S CHAIN, as an alternative to the line above.
    #
    #     im_enhanced = im_log + a*M_1 + b*M_2 + c*M_4 + ...
    #
    # Two differences from everything else here, and they are the point:
    #
    #   ADDITIVE, not multiplicative. The line above is envelope x detail, so a
    #   detail layer's contrast is proportional to the local brightness. Hill
    #   ADDS the masks to the log-mapped image, so a thread of a given contrast
    #   in the log domain renders the same whether it sits at 1.1 R or 6 R.
    #
    #   LINEAR, not normalised. MGN and NAFE each divide by a LOCAL
    #   statistic, which is what lets them show the whole corona at once and
    #   also what gives every part of the frame the same texture amplitude.
    #   Hill's weights are five scalars. Faint structure stays faint, and noise
    #   in the far field stays noise instead of being lifted to signal level.
    #
    # The amplification factors are Hill's own 100:60:20:10 ratios at 2:4:8:16
    # px, extended one octave. They are divided by the rms of the FINEST mask --
    # one number for the whole set, so the ratios between scales are exactly as
    # he gives them -- which is what makes a slider mean the same thing on a
    # different camera.
    _hmix = float(P.get("hillMix", 0.0))
    _hview = (view == "partialconv")
    # Where the log-mapped base sits in the Partial-conv view, leaving the rest
    # of the range for the masks. Hill's own base runs 0.13 to 0.75. See the
    # note where it is applied.
    _HILL_HEADROOM = 0.75
    if _hview and not (getattr(layers, "has_hill", False) and "hilllog" in src):
        return np.full((H, W, 3), 0.5, np.float32)
    if (_hmix > 0 or _hview) and getattr(layers, "has_hill", False) and "hilllog" in src:
        # THE BASE STRETCH IS LIVE, not the one the masks were built at.
        #
        # Hill's k = 1e6 is calibrated for a normalisation whose outer corona
        # sits near 1e-6. Ours does not. On the reference bracket, xn 1.0 at the
        # limb down to 5e-4 in the far field maps to im_log 1.00 down to 0.45 --
        # so the WHOLE corona lands in the top third of the display range, the
        # inner region flattens into a solid blob and the sky lifts to grey.
        # That is what a k too large for the data looks like, and it is the
        # first thing to turn down.
        #
        # The stored layer is invertible, so the base can be re-mapped without
        # rebuilding anything: xn = expm1(im_log * ln(1+Kb)) / Kb exactly
        # recovers the normalised linear image.
        #
        # The masks come along by a scalar. Wherever K*xn >> 1,
        #     log(1+K*xn)/log(1+K) = [ln K + ln xn] / ln(1+K)
        # so an unsharp mask of im_log is an unsharp mask of ln(xn) divided by
        # ln(1+K) -- i.e. the SAME mask at any K, up to that constant. Scaling
        # by ln(1+Kb)/ln(1+Kr) therefore moves the masks to the new base
        # exactly in the regime that carries the corona, and leaves them
        # bounded in the far field, where they were built bounded.
        _kb = float(getattr(layers, "hill_logk", 6.0))
        _kr = float(P.get("hillLogK", _kb))
        _base = np.asarray(src["hilllog"], np.float32)
        _mscale = 1.0
        if abs(_kr - _kb) > 1e-6:
            _Kb, _Kr = 10.0 ** _kb, 10.0 ** _kr
            _xn = np.expm1(np.clip(_base, 0, 1) * np.log1p(_Kb)) / _Kb
            _base = (np.log1p(_Kr * np.clip(_xn, 0, None))
                     / np.log1p(_Kr)).astype(np.float32)
            del _xn
            _mscale = float(np.log1p(_Kb) / np.log1p(_Kr))
        # THE LADDER IS NORMALISED, so the five scale sliders set BALANCE and
        # the master sets STRENGTH. Without this they do both at once: all five
        # at 2 is a sum of 10 against Hill's 1.95, so dragging them up to
        # "see the effect better" multiplied the whole chain by five and drove
        # the composite into clipping -- which looks like a blown-out inner
        # corona, not like too much amplification.
        _HD = (1.0, 0.6, 0.2, 0.1, 0.05)
        _hg = [float(P.get("hill%d" % i, d)) for i, d in enumerate(_HD)]
        _n = min(layers.n_hill, len(_hg))
        _sum = sum(_hg[:_n])
        _norm = (sum(_HD[:_n]) / _sum) if _sum > 1e-6 else 0.0
        _hg = [g * _norm for g in _hg]
        # NORMALISE THE GAIN TO THE MASK'S STRUCTURE, NOT ITS TOTAL SPREAD.
        #
        # `hill_rms[0]` is the 2 px mask's spread over the whole corona, and on
        # a noisy bracket most of that is photon noise, so the slider was
        # effectively calibrated in units of noise: the same setting gave a
        # quiet picture on a clean set and a grainy one on a noisy set, and
        # turning it down to kill the grain took the structure with it at the
        # same rate. build_hill now separates the two (total^2 = structure^2 +
        # resp^2 * sigma^2) and this uses the structure part, so the slider
        # holds the amount of real detail steady and lets the grain fall where
        # the data puts it. Workdirs built before 0.22.84 carry no such figure
        # and fall back to the old normalisation, unchanged.
        _rl = (layers.hill_rms_struct
               if len(getattr(layers, "hill_rms_struct", [])) == layers.n_hill
               else layers.hill_rms)
        _r0 = float(_rl[0]) if _rl else 1.0
        if (os.environ.get("ECLIPSEFORGE_HILL_OLDNORM") != "1"
                and getattr(layers, "hill_ref", 0.0) > 0.0):
            _r0 = layers.hill_ref / HILL_REF_C
        _k = float(P.get("hillGain", 0.06)) * _mscale / max(_r0, 1e-9)
        # SOFT THRESHOLD AGAINST THE PIXEL'S OWN EXPECTED NOISE.
        #
        # An unsharp mask of a signal-free region is not flat: f - blur(f) gives
        # every noise speck a negative ring one to three pixels out. Measured on
        # pure white noise, the mask's autocorrelation reads +1.000, -0.061,
        # -0.027, -0.009 at lags 0..3 while the input noise reads +1.000,
        # -0.001, +0.001 -- so each grain acquires a dark surround, and a bright
        # point with a dark halo is what a lit bump looks like. That is the
        # embossed "3D pattern" in the background, and it is what the method
        # does wherever there is nothing but noise to work on.
        #
        # So Donoho's soft threshold, at k times the sigma this pixel is
        # expected to carry: below it the coefficient is consistent with noise
        # and goes to zero, above it the coefficient keeps its own amplitude
        # minus the threshold. Hill's linearity survives where there is signal,
        # which is the only place it means anything.
        #
        # AT 0 THIS IS EXACTLY THE OLD BEHAVIOUR, coefficient for coefficient.
        _hth = float(P.get("hillDenoise", 0.0))
        _sig = (np.asarray(src["hillsigma"], np.float32)
                if (_hth > 0 and "hillsigma" in src) else None)
        # ... and, in the preview, scaled by what the area average took out of
        # each mask, so the same FRACTION of coefficients falls here as in the
        # export. See the note beside hill_prev_att.
        _att = (layers.hill_prev_att if (preview and
                len(getattr(layers, "hill_prev_att", [])) == layers.n_hill)
                else None)
        _E = np.zeros_like(Y)
        _zs = float(P.get("hillZone", 0.0))
        _zt = hill_zone_t(r, R, float(P.get("hillZoneR", 1.8))) if _zs > 1e-6 else None
        for i in range(_n):
            if not _hg[i]:
                continue
            _mi = np.asarray(src["hill%d" % i], np.float32)
            if _zt is not None and i < len(HILL_ZONE_SIGN) and HILL_ZONE_SIGN[i]:
                # fine masks weigh 1+zone inside and 1-zone outside, coarse the reverse
                _mi = _mi * (1.0 + np.float32(_zs * HILL_ZONE_SIGN[i]) * (1.0 - 2.0 * _zt))
            if _sig is not None and i < len(layers.hill_resp):
                _t = np.float32(_hth * layers.hill_resp[i]
                                * (_att[i] if _att else 1.0)) * _sig
                _mi = np.sign(_mi) * np.maximum(np.abs(_mi) - _t, 0.0)
                del _t
            _E += np.float32(_hg[i]) * _mi
            del _mi
        del _sig
        _E *= np.float32(_k)                    # the masks, in im_log units
        Yh = _base + _E
        _base_keep = _base if _hview else None
        del _base
        # THE HILL RESULT ON ITS OWN, as a monochrome view -- Hill's
        # im_enhanced and nothing else: no envelope, no radial flatten, no
        # prominence term, no disc fill, and NO level match (there is nothing
        # to match it to here). It moves with hillLogK, hillGain, hillDenoise
        # and the five scale sliders, so it can be tuned in the preview and
        # then written out through the normal grayscale export path as a
        # 16-bit TIFF, the same way MGN/FNRGF/NAFE already are.
        if _hview:
            # HEADROOM, so the masks have somewhere to go.
            #
            # Our log map is normalised to put the 99.95th percentile of the
            # luminance at xn = 1, so `_base` reaches 1.000 at the limb and the
            # WHOLE display range is spent before a single mask is added. Every
            # positive coefficient there then clipped, which is why the bright
            # collar just outside the Moon came out as a flat white band with
            # the detail missing -- the one place the masks are strongest.
            #
            # Hill's own log-mapped HDR (his stretching slide) runs from about
            # 0.13 in the sky to 0.75 at the limb, and he says so out loud on
            # that slide: "notice that the inner corona and prominences are not
            # clipped". The room above 0.75 is where a*M1 + b*M2 + ... lands.
            #
            # Scaling base and masks TOGETHER, so this is a display scaling and
            # nothing else: it cannot change the balance between them, and the
            # composite path below is untouched (it level-matches to the
            # multiplicative render's own annulus mean, which already sets its
            # exposure).
            # BASE WEIGHT (0.23.7). Yh = base + k*E. The view shows
            #     hillBase * 0.75 * base + (1 - hillBase) * 0.5 + 0.75 * k * E
            # so at 0 the masks sit on mid grey with the same amplitude they
            # have at 1: the slider moves the base under them and nothing else.
            _hb = float(np.clip(P.get("hillBase", 1.0), 0.0, 1.0))
            _Y = np.float32(_HILL_HEADROOM) * Yh
            del Yh
            if _hb < 1.0:
                _Y += np.float32((1.0 - _hb)) * (np.float32(0.5)
                                                - np.float32(_HILL_HEADROOM) * _base_keep)
            # THE DISC IS PAINTED OUT TO THE MASKS' OWN EDGE (0.23.7). The
            # masks are zero inside Rmask + HILL_LIMB_PAD by construction, so
            # that band showed the bare log base at its brightest -- a bright
            # rim with nothing on it. Hill's disc runs flat to the corona
            # edge; so does this one now, at the base's own value at the
            # Moon's centre (the composite is untouched -- its edge blend
            # covers the same band).
            from .detail import HILL_LIMB_PAD as _HLP
            _pad = (Re - P["discTrim"] / decim) + (_HLP + 1.0) / decim
            _cv = float(_Y[int(round(cy)), int(round(cx))]) if (
                0 <= int(round(cy)) < H and 0 <= int(round(cx)) < W) else 0.5
            _Y = np.where(r < _pad, np.float32(_cv), _Y).astype(np.float32)
            return np.repeat(np.clip(_Y, 0, 1)[:, :, None], 3, axis=2)
        # OVERLAY ONTO THE COMPOSITE (0.23.7), not crossfaded and not
        # multiplied.
        #
        # The crossfade mixed two images with different tone curves (the
        # log-mapped one is flat, so matched in the mean it was brighter at
        # the limb and darker outside: blow-out and lost outer detail). A
        # multiplication by exp(k E) -- Hill's straight addition in the log
        # domain -- is right for LINEAR data, but Y here is already
        # tone-mapped, so it clipped the inner corona just the same.
        #
        # Overlay of a mid-grey layer is bounded by construction: below 0.5
        # the base is multiplied, above 0.5 it is screened, and no mask can
        # push a pixel past white; the contrast it adds tapers toward either
        # end of the range. It is what a high-pass layer does in Photoshop
        # and what Hill points to ("soft light blending mode"). The grey
        # layer is EXACTLY what the Partial-conv view draws at Base weight
        # 0: 0.5 + 0.75 k E. So the view shows the layer, and the slider is
        # its opacity.
        _d = np.clip(np.float32(_hmix * _HILL_HEADROOM) * _E, -0.5, 0.5)
        # THE 0.24 PROMINENCE LAYER: under the layer the corona's Hill detail is
        # eased out, so the trace the colour removal leaves at a prominence's
        # foot does not show through the layer's soft edge
        if getattr(layers, "has_promlayer", False) and float(P.get("promFade", 1.0)) > 0:
            _f = float(P.get("promFade", 1.0))
            if preview:
                _d *= (1.0 - np.float32(_f) * src["pl_F"][:_d.shape[0], :_d.shape[1]])
            else:
                _y0, _y1, _x0, _x1 = layers.pl_box
                _d[_y0:_y1, _x0:_x1] *= (1.0 - np.float32(_f) * layers.pl_full["F"])
        _lo = Y < 0.5
        Y = np.where(_lo, Y * (1.0 + 2.0 * _d),
                     1.0 - (1.0 - Y) * (1.0 - 2.0 * _d)).astype(np.float32)
        del Yh, _E, _d, _lo
    # prominence: local-contrast modulation inside the gate, with a small
    # positive bias so a detected prominence gains presence, not just texture.
    #
    # The contrast term used to be driven by the inner-corona layer, which was
    # the best thing available before 0.18. It is not any more. Scored as
    # correlation with the H-alpha red channel's own fine structure, measured on
    # the bright core of each prominence clear of the disc mask:
    #
    #     prominence   inner-driven   promdet-driven
    #     az 226           0.714          0.828
    #     az  44           0.357          0.361
    #     az  10           0.143          0.284
    #
    # The bias stays at 0.30. Lowering it to 0.10 scores 0.828 against 0.807 on
    # the big prominence, which is not worth changing how bright every user's
    # prominences render.
    #
    # Scoped: `prom` is 0 outside the gate, so the whole term is exactly 1
    # there, and a workdir with no promdet falls back to the old driver and
    # renders bit-identically. Both are checkable, which is the point.
    _pdrv = src["promdet"] if "promdet" in src else inner_eff
    Y = Y * (1 + P["promGain"] * src["prom"] * (0.30 + 1.3 * (_pdrv - 0.5)))
    del _pdrv
    Yd = P["discLevel"] * (1 + P["earthShine"] * (2 * src["earth"] - 1))
    Y = Y * edge + Yd * (1 - edge)
    del B, det, inner_eff, Yd, wf, wI, wG, mg, fn

    # NOT raised to the Saturation exponent here. Saturation is applied AFTER
    # Temperature, Tint and the sky-cast correction below, so that at 0 the
    # chroma collapses to white as balanced, not to the balance gains
    # themselves. Measured before 0.23.8: Saturation at its floor left the
    # whole frame the flat colour of the WB gains -- blue with the corona as
    # the white reference, warm with the sky.
    a = np.clip(src["ratio"], 0.2, 3.0)
    if P.get("bgNeutral", 0) > 0:
        # weighted by the same confidence that built `ratio` -- see the note
        # where _cconf is stored
        _bn = (layers.bg_chroma[None, None, :] ** P["bgNeutral"]) - 1.0
        a = a / (1.0 + src["cconf"][:, :, None] * _bn)
    # THE CORONA AS A WHITE REFERENCE, applied FLAT. This is 0.22.48's line,
    # restored unchanged after ten versions of trying to improve it.
    #
    # Everything added to it afterwards was a WEIGHT -- fade it by the corona
    # fraction (0.22.50), by the chroma confidence, by the sky fraction -- and
    # every one of those weights drew a boundary in the picture, because a
    # weight that varies with radius applied to a gain of about 4.7 in blue IS a
    # visible edge. The tester's verdict on each of them was the same, and the
    # verdict on this one was "0.22.48 was still good".
    #
    # What it does out in the field is not an accident to be corrected: the far
    # field's chroma has already been faded to neutral, and a flat gain turns
    # neutral into the sky-blue that makes the corona read as white against a
    # sky. That is the picture that works.
    # coronaNeutral (0.22.48-0.23.4) was a flat gain from the 1.05-1.6 R
    # annulus at a chosen strength. Since 0.23.5 the page's "Auto: corona"
    # button samples that same ring through this chain and sets Temperature
    # and Tint instead -- one mechanism for white, and one the user can also
    # aim by hand. The key is still honoured here for an old settings file.
    if P.get("coronaNeutral", 0) > 0:
        a *= (layers.corona_gain[None, None, :] ** P["coronaNeutral"])
    # Neutralise radial (0.22.x-0.23.4) lived here: a stored per-annulus colour
    # profile divided out at render time. Removed in 0.23.5. It was measured at
    # stack time on the full pipeline only, stored in the workdir rather than
    # the settings file, never implemented in the browser preview, and
    # unbounded -- a channel median through zero became a gain of a million.
    a[:, :, 0] *= P["temp"]
    a[:, :, 2] /= P["temp"]
    a[:, :, 1] *= P.get("tint", 1.0)
    # Saturation: the balanced chroma's departure from white, raised to a
    # power. 1 leaves it alone, 0 is monochrome, 2 doubles it in the log.
    if abs(P["satur"] - 1.0) > 1e-9:
        a = a ** np.float32(P["satur"])
    # renormalise to unit luminance so colour moves never change brightness
    al = (0.2126 * a[:, :, 0] + 0.7152 * a[:, :, 1] + 0.0722 * a[:, :, 2])
    a /= np.maximum(al, 1e-6)[:, :, None]
    del al
    a = a * edge[:, :, None] + (1 - edge)[:, :, None]   # neutral moon disc
    # PROMINENCE STRUCTURE, CARRIED IN COLOUR RATHER THAN BRIGHTNESS.
    #
    # Measured on the reference bracket: inside a prominence the red channel is
    # at 255 in 100% of the bright core, while green sits near 59 and blue near
    # 51. Red is the max channel everywhere in there, so the hue-preserving knee
    # below sets every channel from `ms` and the chroma ratio -- which means the
    # luminance detail we spent 0.18 and 0.19 improving cannot reach the picture
    # at all. Before the knee a promDetail change is 59 levels; after it, 13.
    #
    # Green and blue have ~200 unused levels. Putting the prominence's own
    # structure there costs nothing in red and shows immediately:
    #
    #                     G spread (p10-p90)   agreement with H-alpha in G
    #   promChroma 0.0        49 - 79 levels             0.368
    #   promChroma 0.6        57 - 105 levels            0.596
    #
    # and on the other two prominences 0.565 -> 0.762 and 0.740 -> 0.827. On
    # screen it is a mean of 13.6 levels, p90 27, against 1.5 for the 0.19.0
    # luminance change. Physically it reads as the denser material going pinker
    # rather than redder, which is what more continuum through more material
    # actually looks like.
    #
    # Scoped like everything else in the gate: `prom` is 0 outside it, so the
    # factor is exactly 1 there (measured residual 1.2e-07, float rounding).
    if abs(P.get("promChroma", 0)) > 1e-6 and "promdet" in src:
        _f = 1.0 + P["promChroma"] * np.clip(src["prom"], 0, 1) * (
            2.0 * np.asarray(src["promdet"], np.float32) - 1.0)
        if abs(P["satur"] - 1.0) > 1e-9:      # follows Saturation, so 0 is grey here too
            _f = np.clip(_f, 1e-3, None) ** np.float32(P["satur"])
        a[:, :, 1] *= _f
        a[:, :, 2] *= _f
        del _f
    rgb = Y[:, :, None] * a
    del a, Y, edge
    # hue-preserving highlight shoulder (parametric knee) + optional white rolloff
    m = rgb.max(axis=2)
    knee = 0.9 - 0.35 * P["hlCompress"]          # 0.9 (off) .. 0.55 (strong)
    ms = np.where(m <= knee, m, knee + (1 - knee) * np.tanh((m - knee) / (1 - knee)))
    scale = np.where(m > 1e-6, ms / np.maximum(m, 1e-6), 1.0)
    rgb *= scale[:, :, None]
    t = np.clip((ms - 0.9) / 0.1, 0, 1) * P["hlDesat"]
    rgb = rgb * (1 - t[:, :, None]) + (ms * t)[:, :, None]
    del m, ms, scale, t
    rgb = np.clip(rgb, 0, 1) ** (1 / P["outGamma"])
    rgb = np.clip((rgb - P["bgBlack"]) / (1 - P["bgBlack"]), 0, 1)
    _eg = np.clip((r - (Re + 2.0 / decim - _dsw)) / _dsw, 0, 1)
    # with an inner corona level, contain whites moves into look() (see there)
    rgb = apply_tone(np.ascontiguousarray(rgb, np.float32), P.get("toneCurve") or [],
                     1.0 if float(P.get("coronaLevel", 0) or 0) > 0 else float(P.get("containWhites", 1.0)),
                     moon_pts=P.get("moonCurve") or [], eg=_eg)
    # the look controls -- see look()
    rgb = look(rgb, P, r, R, _eg)
    del _eg
    # STRUCTURE, then FINAL SHARPENING, on the finished picture, before the ring overlay.
    if view in ("composite", "moon") and float(P.get("structure", 0.0)) > 1e-6:
        rgb = structure_boost(np.ascontiguousarray(rgb, np.float32), cy, cx, R,
                              float(P["structure"]))
    if max(float(P.get("finalMid", 0)), float(P.get("finalDark", 0)),
           float(P.get("finalLight", 0))) > 1e-6:
        rgb = final_sharpen(np.ascontiguousarray(rgb, np.float32), R,
                            float(P["finalMid"]), float(P["finalDark"]),
                            float(P["finalLight"]), float(P.get("finalNoise", 0)))
    # THE 0.24 PROMINENCE LAYER, over the finished picture (after the final
    # sharpening, so its soft edge is not sharpened into a line)
    if view == "composite" and getattr(layers, "has_promlayer", False):
        from .promlayer import overlay as _pov
        if preview:
            _m = {k: src["pl_" + k] for k in ("L", "C", "W1", "W2", "F")}
            _hh, _ww = rgb.shape[:2]
            _m = {k: v[:_hh, :_ww] for k, v in _m.items()}
            _pov(rgb, _m, P)
        else:
            _y0, _y1, _x0, _x1 = layers.pl_box
            _sub = np.ascontiguousarray(rgb[_y0:_y1, _x0:_x1])
            _pov(_sub, layers.pl_full, P)
            rgb[_y0:_y1, _x0:_x1] = _sub
    if P["ringBlend"] > 0 and "contact" in src:
        c = src["contact"]
        sc = P["ringScale"]
        rdy, rdx = P["ringDY"] / decim, P["ringDX"] / decim
        if abs(sc - 1) > 1e-4 or abs(rdy) > 1e-3 or abs(rdx) > 1e-3:
            # out(y,x) = in(cy + (y-cy-rdy)/sc, cx + (x-cx-rdx)/sc)
            mat = np.array([[1 / sc, 0], [0, 1 / sc]], np.float64)
            off = [cy - (cy + rdy) / sc, cx - (cx + rdx) / sc]
            ct = np.empty_like(c)
            for ch in range(3):
                ct[:, :, ch] = ndimage.affine_transform(
                    c[:, :, ch], mat, offset=off, order=1, mode="constant", cval=0)
            c = ct
        rgb = 1 - (1 - rgb) * (1 - P["ringBlend"] * c)
    # THE MOON VIEW (lab 0.24): the finished Moon on its own -- the composite's
    # disc (level, earthshine, colour, its own curve), black outside it, the
    # disc's soft edge as the fade. What the Moon-layer export carries.
    if view == "moon":
        rgb = rgb * (1.0 - moon_alpha(layers, P, decim, rgb.shape[:2]))[:, :, None]
    return rgb


def moon_alpha(layers, P, decim, shape):
    """1 on the Moon's disc, 0 outside, the disc's soft edge between -- the
    same ramp the render uses to put the disc in (1 - edge)."""
    H, W = shape
    cy, cx, R = layers.geometry(decim)
    yy = np.arange(H, dtype=np.float32)[:, None] - cy
    xx = np.arange(W, dtype=np.float32)[None, :] - cx
    r = np.sqrt(yy * yy + xx * xx)
    Re = layers.mask_radius_map((H, W), decim) + float(P.get("discTrim", 0.0)) / decim
    _dsw = max(float(P.get("discSoft", 12.0)), 1.0) / decim
    return (1.0 - np.clip((r - (Re + 2.0 / decim - _dsw)) / _dsw, 0, 1)).astype(np.float32)


def export_moon_layer(layers, params, out_path, composite, size="full"):
    """THE MOON AS ITS OWN LAYER (lab 0.24, Nico): a 16-bit RGBA TIFF on the
    export's grid -- the composite's own disc pixels, alpha = the disc with its
    soft edge. Laid over the composite (Normal) it changes nothing; it is there
    to be edited on its own."""
    import tifffile
    from . import icc
    P = dict(defaults_for(getattr(layers, "mode", None))); P.update(params or {})
    a = moon_alpha(layers, P, 1, composite.shape[:2])
    full = np.concatenate([np.asarray(composite, np.float32) * (a[:, :, None] > 0),
                           a[:, :, None]], -1)
    full = apply_orient(full, P.get("orient", ""))
    if size == "half":
        h2, w2 = full.shape[0] // 2 * 2, full.shape[1] // 2 * 2
        f = full[:h2, :w2].reshape(h2 // 2, 2, w2 // 2, 2, 4)
        aa = f[..., 3].mean(axis=(1, 3))
        c = (f[..., :3] * f[..., 3:4]).sum(axis=(1, 3)) / np.maximum(
            f[..., 3].sum(axis=(1, 3)), 1e-9)[..., None]
        full = np.concatenate([c, aa[..., None]], -1)
    prof = icc.srgb_profile()
    tifffile.imwrite(out_path, (np.clip(full, 0, 1) * 65535 + 0.5).astype(np.uint16),
                     photometric="rgb", extrasamples=[2], compression="zlib",
                     extratags=[(34675, 1, len(prof), prof, False)])


def _range_note(a, gray):
    """What fraction of the container this view actually uses.

    A detail view is rendered at the slider setting it will contribute to the
    composite, not at a setting that fills the histogram: measured on the
    reference render, MGN occupies 14% of the 16-bit range where FNRGF occupies
    81%. Someone compositing the layers by hand needs to know that before they
    stretch one in another program, so it goes in the file.
    """
    try:
        lo, hi = np.percentile(a, [1, 99])
        return "p1 %.4f p99 %.4f (%.0f%% of full scale)" % (lo, hi, 100 * (hi - lo))
    except Exception:
        return "unknown"


def _png_add_iccp(path, prof):
    """Insert an iCCP chunk into a PNG that was written without one.

    OpenCV writes the only 16-bit PNG we can produce and has no way to embed a
    profile, so this file went out untagged -- the exact case the sRGB tagging
    was added to prevent, and the one that makes PixInsight read it as linear.
    Rewriting the container is cheap and lossless: iCCP goes after IHDR, which
    is where the spec requires it (PNG 1.2, 4.2.2.4).
    """
    import struct
    import zlib
    try:
        with open(path, "rb") as f:
            blob = f.read()
        if blob[:8] != b"\x89PNG\r\n\x1a\n" or b"iCCP" in blob[:4096]:
            return
        ihdr_len = struct.unpack(">I", blob[8:12])[0]
        cut = 8 + 12 + ihdr_len                  # end of the IHDR chunk
        body = b"sRGB\0\0" + zlib.compress(prof, 9)
        chunk = (struct.pack(">I", len(body)) + b"iCCP" + body
                 + struct.pack(">I", zlib.crc32(b"iCCP" + body) & 0xFFFFFFFF))
        with open(path, "wb") as f:
            f.write(blob[:cut] + chunk + blob[cut:])
    except Exception:
        pass                                     # a missing tag is not fatal


def thumbnail(arr, path, width=480):
    """A small sRGB JPEG of a finished picture (float 0..1, HxW or HxWx3) for
    the Versions strip: box-binned to about `width` px wide, so a 44 Mpx
    export costs a few hundred ms and 30 kB."""
    a = np.asarray(arr, np.float32)
    if a.ndim == 2:
        a = np.repeat(a[:, :, None], 3, axis=2)
    f = max(1, int(round(a.shape[1] / float(width))))
    H2, W2 = a.shape[0] // f * f, a.shape[1] // f * f
    if f > 1:
        a = a[:H2, :W2].reshape(H2 // f, f, W2 // f, f, 3).mean(axis=(1, 3))
    Image.fromarray((np.clip(a, 0, 1) * 255 + 0.5).astype(np.uint8)).save(
        path, quality=86, optimize=True)
    return path


def export(layers: Layers, params, fmt, out_path, view="composite", size="full",
           notes=None, thumb=None, keep=None):
    """fmt: tif16 | tif8 | png (16-bit when OpenCV present, else 8-bit) | jpg.
    notes: an optional list the caller passes in to collect warnings about the
    file that was just written -- see the black-point check below.
    view: composite | mgn | fnrgf | nafe | inner | prom | tangential |
    partialconv | flat (detail views export grayscale).
    size: full | half (half = 2x2 binned, ~2x better SNR)."""
    rgb = render(layers, params, preview=False, view=view)
    if keep is not None:          # the rendered picture, for the prominence-layer export
        keep["rgb"] = rgb
    # THE BLACK POINT CLIPS EACH CHANNEL ON ITS OWN, and used to say nothing
    # when it did. Where green and blue fall under it and red does not, the
    # outer field prints as flat saturated red with the structure still visible
    # in it -- a tester's export did exactly that, and looked like a stacking
    # fault rather than a slider. Count it and say so.
    if notes is not None and view == "composite":
        try:
            _bb = float((params or {}).get("bgBlack", 0.0) or 0.0)
            if _bb > 0:
                _h, _w = rgb.shape[:2]
                _dy = np.arange(_h, dtype=np.float32)[:, None] - layers.cy
                _dx = np.arange(_w, dtype=np.float32)[None, :] - layers.cx
                _out = (_dy * _dy + _dx * _dx) > (layers.Rmask + 8.0) ** 2
                _n = float(_out.sum())
                if _n > 1000.0:
                    _flat = (_out & (rgb[:, :, 1] <= 0.0) & (rgb[:, :, 2] <= 0.0)
                             & (rgb[:, :, 0] > 0.2))
                    _f = float(_flat.sum()) / _n
                    if _f > 0.01:
                        notes.append(
                            "WARNING: Black point %.3f drives green and blue to "
                            "zero on %.0f%% of the corona while red stays up, so "
                            "that area prints as flat red. Lower Black point, or "
                            "correct the colour before raising it." % (_bb, 100.0 * _f))
        except Exception:
            pass
    rgb = apply_orient(rgb, (params or {}).get("orient", ""))
    if size == "half":
        H2, W2 = rgb.shape[0] // 2 * 2, rgb.shape[1] // 2 * 2
        rgb = rgb[:H2, :W2].reshape(H2 // 2, 2, W2 // 2, 2, -1).mean(axis=(1, 3))
        if rgb.shape[-1] == 1:
            rgb = rgb[:, :, 0]
        rgb = np.ascontiguousarray(rgb)
    gray = view not in ("composite", "moon")
    arr16 = (rgb[:, :, 0] if gray else rgb)
    # MONOCHROME composite: one channel, the picture's luminance
    if view == "composite" and int((params or {}).get("colourMode", 0) or 0) == 1 and rgb.ndim == 3:
        gray = True
        arr16 = np.ascontiguousarray(0.2126 * rgb[:, :, 0] + 0.7152 * rgb[:, :, 1] + 0.0722 * rgb[:, :, 2])
    # Tag the file. The render is already in display encoding (the browser
    # preview shows these exact bytes as sRGB), so sRGB is what it IS -- but it
    # was going out untagged, leaving every host to guess. Photoshop guessed
    # right by accident; PixInsight assumes linear and got it wrong.
    from . import icc
    _prof = icc.srgb_profile()
    # ONLY ON RGB (0.22.78). This is an RGB sRGB profile, and it was being
    # written into minisblack TIFFs and mode-L JPEGs, where the profile class
    # does not match the data -- some hosts flag it, others apply it. A detail
    # view is a single channel of measurement, not colour, so it goes untagged.
    _tags = [] if gray else [(34675, 1, len(_prof), _prof, False)]
    if fmt == "tif16":
        import tifffile
        tifffile.imwrite(out_path, (arr16 * 65535 + 0.5).astype(np.uint16),
                         compression="zlib", extratags=_tags,
                         photometric="minisblack" if gray else "rgb",
                         description="eclipseforgehdr %s, range %s, params: %s"
                         % (view, _range_note(arr16, gray), json.dumps(params)))
    elif fmt == "tif8":
        import tifffile
        tifffile.imwrite(out_path, (arr16 * 255 + 0.5).astype(np.uint8),
                         compression="zlib", extratags=_tags,
                         photometric="minisblack" if gray else "rgb")
    elif fmt == "png":
        arr = (arr16 * 65535 + 0.5).astype(np.uint16)
        try:
            import cv2
            cv2.imwrite(out_path, arr if gray else arr[:, :, ::-1])
            if not gray:
                _png_add_iccp(out_path, _prof)   # cv2 writes no profile at all
        except Exception:
            Image.fromarray((arr16 * 255 + 0.5).astype(np.uint8)).save(
                out_path, icc_profile=(None if gray else _prof))
    elif fmt == "jpg":
        Image.fromarray((arr16 * 255 + 0.5).astype(np.uint8)).save(
            out_path, quality=92, optimize=True,
            icc_profile=(None if gray else _prof))
    else:
        raise ValueError(fmt)
    if thumb:
        try:
            thumbnail(arr16, thumb)
        except Exception:
            pass
    return out_path
