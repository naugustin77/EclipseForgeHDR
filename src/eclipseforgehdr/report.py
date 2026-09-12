"""Run report: what went in, what was measured, what was done to it."""
from __future__ import annotations
import os, json, math
import numpy as np

# Apparent lunar diameter during a total eclipse. It has to exceed the Sun's
# 1890"-1950", and varies with the Moon's distance; 1920" +/- 4% covers the
# realistic range. Used only to quote an approximate plate scale.
MOON_ARCSEC = 1920.0
MOON_ARCSEC_TOL = 0.04

METHODS = [
    ("Frame selection",
     "per-frame sharpness score (high-pass energy over mean signal, saturation-masked); "
     "best frame per tier keeps the reference, the rest are aligned to it and averaged"),
    ("Sensor defects",
     "hot/dead photosites mapped on the shortest tier against a fitted photon+read "
     "noise model, repaired by the median of same-colour neighbours"),
    ("Alignment",
     "SEMI-PHASE correlation of a gradient-flattened log corona, high-passed at "
     "25 px, with lag-1 and lag-2 links solved together by weighted least "
     "squares. The cross-power spectrum is divided by the two amplitude spectra "
     "plus a constant p = q = 1e-4 of the peak amplitude (Druckmullerova, "
     "Phase-correlation based image registration, 2010, Def. 3.25) -- an "
     "ADDITIVE constant, not an exponent, so it acts as a noise floor on the "
     "whitening: frequencies well above it are flattened, frequencies carrying "
     "only noise are left alone. That is the middle ground between the two ends "
     "this pipeline had tested. Plain cross-correlation (0.23.1 and earlier) and "
     "full whitening are both still selectable. Measured on tiers built from a "
     "real bracket and its own master flat, with the shift injected by an exact "
     "phase ramp: 1.524 px rms error and a 0.200 px pull toward zero shift for "
     "plain cross-correlation, against 1.187 px and 0.041 px for semi-phase -- "
     "the pull being the estimator locking onto what both frames share, which "
     "is the sensor's own dust and fixed pattern"),
    ("Prominence anchors",
     "the fastest tiers hold almost no corona to correlate on, so they are also tied in "
     "by normalized cross-correlation of prominence patches -- solar features, unlike the "
     "lunar limb, which drifts against the corona during totality. Anchors are found as "
     "compact azimuthal peaks on the limb; a link is used only while its anchors agree "
     "with each other, and enters the same least-squares network as an extra link"),
    ("Demosaic",
     "Malvar-He-Cutler gradient-corrected bilinear interpolation; the clipping "
     "test is asked of the MOSAIC, so a pixel interpolated in part from a "
     "saturated photosite is flagged (the kernels have negative lobes, so such "
     "a pixel can otherwise read below threshold)"),
    ("Photometric calibration",
     "per-tier scale factors from median ratios of overlapping unsaturated signal, "
     "chained to the middle tier; plus a per-tier affine transform varying with "
     "azimuth (k_i(phi), q_i(phi)) fitted in 60 angular segments against the "
     "running composite and smoothed by a trigonometric polynomial, after "
     "Druckmullerova's LDIC (doctoral thesis eq. 4.15), applied mean-preserving "
     "so it corrects only what a scalar cannot"),
    ("HDR merge",
     "weighted average of the calibrated tiers in LINEAR camera-native colour, "
     "then to sRGB primaries. The weight is the tier's exposure time raised to "
     "the merge exponent, times a roll-off that reaches zero at saturation and "
     "again at a floor of 0.5% of it. This is the Debevec and Malik (1997) "
     "framework -- a weighted average of exposures put on a common scale -- but "
     "NOT their weighting: eq. 4 there is a hat on the PIXEL VALUE with no "
     "exposure term, and eq. 6 averages in the log domain. Their hat exists "
     "because an unknown nonlinear film or 8-bit response is least reliable at "
     "its toe and shoulder; raw sensor data has no toe, so the only reason to "
     "distrust a low value is noise. The estimation-theoretic weight for that "
     "case is Robertson, Borman and Stevenson (2003) eq. 7, which weights each "
     "tier by its certainty times the exposure time SQUARED -- reducing to an "
     "exponent of 1 where photon noise dominates and 2 where read noise does. "
     "This merge's exponent is chosen by a detail trial instead, and can sit "
     "below that range: it is buying sharpness, not signal-to-noise"),
    ("Limb fit",
     "50% crossing between disc and near-limb corona level along 720 rays, "
     "robust fit of r = R + dx cos t + dy sin t, re-centred over five iterations"),
    ("MGN detail",
     "Multiscale Gaussian Normalization (Morgan & Druckmuller 2014, Sol. Phys. 289, 2945), "
     "photon-noise-adaptive, disc excluded by normalized convolution, radial profile "
     "removed before normalization. Their scale ladder, their per-scale gains "
     "and their gamma, but NOT their eq. 5 global term: the paper adds a "
     "globally gamma-normalised copy of the whole image at h = 0.7 so an AIA "
     "disc keeps its own large-scale context, and here that weight is 0.12 "
     "(0 on the inner-corona passes) because this is a detail layer and the "
     "renderer builds the brightness envelope separately -- carrying it twice "
     "would double it"),
    ("FNRGF detail",
     "Fourier Normalizing Radial Gradient Filter (Morgan, Habbal & Woo 2006; "
     "Druckmullerova et al. 2011, ApJS 194, 25) with coverage-matched Fourier order "
     "and Huber-weighted robust fitting"),
    ("Tangential filter",
     "tangential (rotational) unsharp about the disc centre, normalized convolution "
     "in polar space. BLIND TO TANGENTIAL STRUCTURE by construction: it "
     "averages around the Sun, so radially oriented features such as plumes "
     "are enhanced most and tangentially oriented ones -- the tops of loops "
     "and helmet streamers -- are not enhanced at all (Druckmullerova, thesis "
     "sec. 5.1; Druckmuller, Rusin & Minarovjech 2006, sec. 4.1). A moving "
     "average also has a sinc frequency response, so its gain is not monotone "
     "in frequency, and its phase spectrum is not flat, so structures can be "
     "shifted. Treat this layer as a look rather than a measurement. Method "
     "traced to Espenak (2000), from Photoshop's Radial Blur / Spin; an "
     "earlier attribution here to Pellett via Druckmuller 2009 is unverified"),
    ("Partial convolution",
     "unsharp masks of the log-stretched image at 2, 4, 8, 16 and 32 px, added back "
     "with weights 1 / 0.6 / 0.2 / 0.1 / 0.05. The blur behind each mask is taken in "
     "polar coordinates, so it averages along a streamer rather than across it, and "
     "it is PARTIAL: the occulted disc and the prominences are excluded from both the "
     "numerator and the normalisation, so nothing is smeared out of them into the "
     "corona. Additive and linear, where MGN and NAFE are multiplicative and "
     "locally normalised. Masks are soft-thresholded against a per-pixel photon-noise "
     "model carried through the derivative of the log map. After Jonathan Hill, "
     "\"Advanced Solar Eclipse Photography\""),
    ("Inner corona",
     "separate MGN of a saturation-weighted stack of the shortest tiers -- up to "
     "four, and only those within 24x of the shortest, so a bracket whose fourth "
     "tier is already long does not drag a bright frame into it"),
    ("Prominences",
     "H-alpha redness R/((G+B)/2) in a single fast tier, thresholded against the "
     "robust spread of the corona's own colour in the limb annulus"),
    ("NAFE detail",
     "Noise Adaptive Fuzzy Equalization with a Variable Neighbourhood (Druckmuller 2013, "
     "ApJS 207:25; Druckmuller & Druckmullerova, IWCIA 2014, LNCS 8466, 262-271): each pixel is "
     "ranked within a fuzzy multiscale neighbourhood restricted in VALUE rather than by a "
     "geometric mask, which is what removes the contrast loss at the lunar edge and makes "
     "this the one detail layer that does not depend on the limb fit"),
    ("Denoise",
     "a-trous (starlet) multiscale soft thresholding of log luminance against a "
     "per-pixel photon-noise sigma"),
]


def _fmt_exp(s):
    return f"1/{1 / s:.0f}s" if s < 0.4 else f"{s:g}s"


def build(stats):
    """stats -> plain-text report."""
    L = []
    A = L.append
    A("EclipseForgeHDR run report")
    A("=" * 60)
    A(f"version      : {stats.get('version', '?')}")
    A(f"folder       : {stats.get('folder', '?')}")
    if stats.get("finished"):
        A(f"processed    : {stats['finished']}")
    _tm = stats.get("timing") or {}
    if _tm.get("total_s"):
        def _d(s):
            s = int(round(float(s)))
            return f"{s // 60}m{s % 60:02d}s" if s >= 60 else f"{s}s"
        A(f"run time     : {_d(_tm['total_s'])}")
        # Every step over a second, not just the slowest few: the point of this
        # block is to make the progress bar's weights measurable, and a list
        # that stops at the top four leaves the rest of the run unaccounted for.
        _steps = _tm.get("steps") or _tm.get("slowest") or []
        _tot = float(_tm["total_s"]) or 1.0
        for d, m in _steps:
            A(f"               {_d(d):>7}  {100 * float(d) / _tot:4.1f}%  {m}")
        if _tm.get("accounted_s"):
            _un = float(_tm["total_s"]) - float(_tm["accounted_s"])
            A(f"               {_d(_un):>7}  {100 * _un / _tot:4.1f}%  "
              f"(steps under 1s, not listed)")
    cam = stats.get("camera_info") or {}
    if cam:
        bits = [cam.get("camera"), cam.get("lens")]
        if cam.get("focal_mm"):
            bits.append(f"{cam['focal_mm']:g} mm")
        if cam.get("f_number"):
            bits.append(f"f/{cam['f_number']:g}")
        A("camera       : " + "  ".join(b for b in bits if b))
    if stats.get("iso"):
        A(f"ISO          : {stats['iso']}")
    if stats.get("shot_first"):
        A(f"frames taken : {stats['shot_first']}  ..  {stats.get('shot_last', '')}")

    if stats.get("mode") == "imported HDR":
        A("")
        A("IMPORTED HDR")
        A("-" * 60)
        A(f"source       : {os.path.basename(str(stats.get('imported', '?')))}")
        _tc = {"srgb": "sRGB transfer function (inverted on import)",
               "linear": "linear (used as it is)",
               "gamma": f"gamma {stats.get('import_gamma') or '?'} (inverted on import)"
               }.get(stats.get("import_tone"), "not declared")
        A(f"tone curve   : {_tc}")
        A("This image was not stacked here, so there is nothing to report about")
        A("alignment, photometry or per-tier lunar masking, and no earthshine")
        A("layer was built. Two layers are weaker than they would be from a")
        A("bracket, and it is worth knowing which:")
        A("  * inner corona: normally a separate MGN of the shortest tiers, which")
        A("    see the inner corona unsaturated. From one image it is a second")
        A("    view of the same pixels -- a different filter, not new data.")
        A("  * prominences : the gate keys on H-alpha redness in a fast tier. A")
        A("    merged image has usually compressed exactly that, so the gate has")
        A("    less to work with and flags less.")
    _imported = stats.get("mode") == "imported HDR"
    tiers = stats.get("tiers", [])
    if not (_imported and not tiers):
        A("")
        A("EXPOSURE STACK")
        A("-" * 60)
        # THE SIMPLE PATH FILLS NONE OF THESE COLUMNS. It scores no frame, fits
        # no per-tier factor and solves no shift network, so printing 0.0 /
        # 1.00 / 1.000 / 0.0 for every tier reads as "measured, and perfect"
        # when the truth is "never measured". Show what it does have instead.
        _simple_tab = stats.get("mode") == "simple stack"
        if _simple_tab:
            A(f"{'exposure':>10}  {'used':>6}   {'relative exposure, MEASURED from the pixels':<44}")
            _lad = {round(float(x['to']), 9): x for x in (stats.get('simple_ladder') or [])}
            _rel = 1.0
            for t in tiers:
                _L = _lad.get(round(float(t['sec']), 9))
                if _L:
                    _rel *= float(_L['measured'])
                    _txt = (f"x{_rel:<11.4g} step {_L['measured']:.4f}"
                            f"  (header said {_L['header']:.4f})")
                else:
                    _txt = f"x{_rel:<11.4g} reference tier"
                A(f"{_fmt_exp(t['sec']):>10}  "
                  f"{str(t['n']) + '/' + str(t.get('n_avail', t['n'])):>6}   {_txt}")
        else:
            A(f"{'exposure':>10}  {'used':>6}  {'best frame':<18} {'sharp':>7} {'spread':>6}  "
              f"{'photom.':>7}  {'shift px':>9}")
            for t in tiers:
                sh = t.get("shift", [0, 0])
                A(f"{_fmt_exp(t['sec']):>10}  "
                  f"{str(t['n']) + '/' + str(t.get('n_avail', t['n'])):>6}  "
                  f"{t.get('best', ''):<18} {t.get('sharpness', 0):>7.1f} "
                  f"{t.get('spread', 1):>6.2f}  {t.get('cal', 1):>7.3f}  "
                  f"{math.hypot(sh[0], sh[1]) * 2:>9.1f}")
    if tiers:
        secs = [t["sec"] for t in tiers]
        n_used = sum(t["n"] for t in tiers)
        integ = sum(t["sec"] * t["n"] for t in tiers)
        A("")
        A(f"tiers        : {len(tiers)}   frames stacked: {n_used}"
          f"   files found: {stats.get('n_files', n_used)}")
        A(f"exposure span: {_fmt_exp(min(secs))} .. {_fmt_exp(max(secs))}"
          f"   = {math.log2(max(secs) / min(secs)):.1f} EV")
        A(f"integration  : {integ:.2f} s total open-shutter time")
        _span = stats.get("simple_ladder_span")
        if stats.get("mode") == "simple stack" and _span:
            _h, _m = float(_span["header"]), float(_span["measured"])
            A(f"photometric  : no tier factors are fitted on this path. The "
              f"ladder is MEASURED:")
            A(f"             : the headers span 1 .. {_h:.4g}, the pixels span "
              f"1 .. {_m:.4g}"
              + (f"  ({_h / _m:.2f}x apart)" if _m > 0 else ""))
            if _m > 0 and (_h / _m > 1.05 or _m / _h > 1.05):
                A(f"             : so the exposure times in these headers do "
                  f"NOT predict the signal, and")
                A(f"             : the merge used the measured ratios. Check "
                  f"how the files were written.")
        else:
            cals = [t.get("cal", 1) for t in tiers]
            A(f"photometric  : tier factors {min(cals):.3f} .. {max(cals):.3f} "
              f"(1.000 = exposure time exactly predicts signal)")
    if stats.get("align_residual") is not None:
        A(f"alignment    : {stats['align_residual']:.2f} px max network residual "
          f"(half-res) = {stats['align_residual'] * 2:.2f} px full-res")
    pa = stats.get("prom_align") or {}
    if pa.get("used"):
        A(f"prominence links: {pa.get('anchors', 0)} anchors on the "
          f"{_fmt_exp(pa['tier'])} tier tie in {pa['used']} tier(s)")
        sp = [l["spread"] for l in pa.get("links", []) if l.get("spread") is not None]
        if sp:
            A(f"             : anchor spread {min(sp):.2f}..{max(sp):.2f} px "
              f"(links above {2.0:.1f} px are not used)")
    ig = stats.get("inner_geom") or {}
    if ig.get("offset_px") is not None:
        A(f"inner stack  : its own lunar limb sits {ig['offset_px']:.0f} px from the "
          f"merged one (the short tiers were shot earlier); this layer is masked "
          f"with its own disc, not the merged one")
    mt = stats.get("moon_track") or {}
    if mt:
        A(f"lunar motion : {mt['drift_px_per_s']:.2f} px/s, "
          f"{mt['drift_px_total']:.0f} px across the bracket "
          f"(scatter about the straight-line track "
          f"{mt.get('scatter_y_px', 0):.0f}/{mt.get('scatter_x_px', 0):.0f} px)")
    mm = stats.get("moon_mask") or {}
    if mm.get("verdict"):
        A(f"moon masking : {mm['verdict']}")
        if "limb_width_off_px" in mm:
            A(f"             : merged limb {mm['limb_width_off_px']:.1f} px "
              f"unmasked vs {mm['limb_width_on_px']:.1f} px masked; "
              f"circle-fit rms {mm.get('rms_off', 0):.2f} -> {mm.get('rms_on', 0):.2f}")
    aq = stats.get("align_quality") or {}
    if aq:
        A("alignment check:")
        if "cov_limb" in aq:
            A(f"             : tier-to-tier variance {aq['cov_limb']:.3f} at the limb, "
              f"{aq.get('cov_corona', float('nan')):.3f} in the corona "
              f"(lower = the tiers agree; over "
              f"{aq.get('tiers_limb', float('nan')):.0f} unclipped tiers at the "
              f"limb)")
        elif "cov_limb_unmeasurable" in aq:
            A(f"             : tier-to-tier variance at the limb NOT MEASURABLE "
              f"— only {aq['cov_limb_unmeasurable']:.0f} tier(s) hold unclipped "
              f"signal between 1.00 and 1.10 R, and the coefficient of "
              f"variation needs three. A number built on a partly clipped tier "
              f"is worse than none: that is what produced a false 0.79 here for "
              f"three releases. A shorter tier at the top of the bracket would "
              f"make it measurable")
        if "rim_width_px" in aq:
            import math as _m
            A("             : disagreement rim "
              + (f"{aq['rim_width_px']:.0f} px wide just outside the limb"
                 if _m.isfinite(aq["rim_width_px"])
                 else "none detectable — the tiers agree at the limb"))
        if "limb_width_med" in aq:
            A(f"             : merged limb 20-80% transition "
              f"{aq['limb_width_med']:.1f} px (p90 "
              f"{aq.get('limb_width_p90', float('nan')):.1f} px)")
    _ls = stats.get("ldic_k_spread_pct")
    if _ls:
        _w = max(_ls.items(), key=lambda kv: kv[1])
        A(f"azimuthal fit: per-tier gain and offset varying with azimuth "
          f"(Druckmullerova thesis eq. 4.15), fitted in 60 segments against "
          f"the running composite and mean-preserved. Worst tier varied "
          f"{_w[1]:.0f}% around the limb ({_fmt_exp(float(_w[0]))}); a single "
          f"photometric scalar per tier cannot express that. Corrected on "
          f"{len(_ls)} tiers")
    _ce = stats.get("cfa_clip_extra")
    if _ce:
        _w = max(_ce.items(), key=lambda kv: kv[1])
        A(f"clipping test: asked of the MOSAIC since 0.22.26, not of the "
          f"demosaiced result — a pixel built in part from a saturated "
          f"photosite is now flagged. This run marks up to {_w[1]:.2f}% more of "
          f"the frame invalid (worst tier {_w[0]}s); those pixels were entering "
          f"the merge at full weight carrying interpolated clipped data")
    _ac = stats.get("all_clipped_px")
    if _ac:
        A(f"bracket reach: {_ac} px are saturated in EVERY tier, the shortest "
          f"one included — nothing in this bracket measures them. They are "
          f"filled from the shortest tier, so they are a LOWER BOUND on the "
          f"real brightness, not a measurement; before 0.22.36 they came out "
          f"as an unrecoverable black hole. This is almost always a prominence "
          f"core. One more stop at the short end of the bracket fixes it at "
          f"the camera, and nothing else can")
    _fm = stats.get("feather_mode"); _fr = stats.get("feather_ratio")
    if _fm:
        from .pipeline import FEATHER_NAMES as _FN
        _n = _FN.get(_fm, _fm)
        A(f"merge weight : {_n} — a SETTING, not a measurement. Blended edge "
          f"blurs the merge weight across each tier's clipping boundary, so "
          f"the tiers join without a seam and a ring artifact just outside the "
          f"limb is suppressed. It is suppressed by covering it with an error "
          f"of the same shape: weight leaks from clipped pixels into unclipped "
          f"ones and the corona reads low there — on some brackets a few "
          f"percent, on others enough to show as a coloured rim. Exact edge "
          f"keeps the weight exact at that boundary, which is radiometrically "
          f"correct and leaves the rings visible. Change it in the toolbar")
        if _ac and _fm == "plain":
            A(f"             : note this weight also papers over the "
              f"all-tiers-clipped pixels above, using leaked weight from "
              f"unclipped neighbours. It looks better there and the values are "
              f"no more real than the fallback's")
        if _fr:
            A(f"             : the trial measured the Detail weight at "
              f"{100 * _fr:.0f}% of the Photometric level at 1.02 R. ADVICE "
              f"ONLY — this estimator read 99% on a bracket the offline bench "
              f"puts at 13%, so it does not choose and should not be trusted "
              f"over the picture")
    _lg = stats.get("linearity_gamma")
    if _lg is not None and abs(_lg - 1.0) > 0.08:
        A(f"linearity    : *** NOT scene-linear — the photometric links are "
          f"tilted as if the data carried a gamma of about {_lg:.2f}, where "
          f"linear data gives 1.00. A raw developer's tone curve does this, "
          f"and inverting sRGB on load does not undo it. Use the raw files; a "
          f"flat exported the same way carries the same curve and does not "
          f"divide out. Every number in this report assumes linearity")
    _pd = stats.get("pedestal")
    if _pd:
        # sigma comes from a leave-one-tier-out jackknife on a grid, so it can
        # land at exactly 0 when every subset agrees. Say "well determined"
        # rather than printing a number the grid cannot support.
        _sg = (f"{abs(_pd['fitted']) / _pd['sigma']:.0f} sigma"
               if _pd["sigma"] > 1e-9 else
               "every leave-one-tier-out subset agrees")
        if _pd["applied"]:
            A(f"pedestal     : {_pd['applied']:+.2f} ADU subtracted from every "
              f"tier (fitted {_pd['fitted']:+.2f}; {_sg}) — a black level left "
              f"a few ADU behind arrives divided by the exposure time, so it is "
              f"nothing on a long tier and everything on a short one")
            A(f"             : tier-to-tier disagreement in the outer field "
              f"{100 * _pd['scatter_before']:.1f}% -> "
              f"{100 * _pd['scatter_after']:.1f}%")
        else:
            A(f"pedestal     : none applied. The best fit was "
              f"{_pd['fitted']:+.2f} ADU ({_sg}) and was shrunk to zero by its "
              f"own uncertainty; the tiers already agree to "
              f"{100 * _pd['scatter_before']:.1f}% in the outer field")
    for _k, _msg in (
        ("limb_spread_bad",
         "the tiers' lunar limbs are spread over %s px -- they are the same Moon "
         "seconds apart, so this is a cross-tier alignment failure"),
        ("limb_fit_rms_bad",
         "the limb fit's rms is %s px; a real lunar limb fits to well under 1%% "
         "of its radius, so this circle does not describe an edge"),
        ("limb_variance_bad",
         "the tiers disagree by %s (coefficient of variation) at the limb where "
         "a well-behaved set sits near 0.05-0.08 — a disagreement in BRIGHTNESS, "
         "not position, which the detail filters amplify into concentric rings. "
         "The CAUSE is not established: 0.22.17 withdrew the veiling-glare "
         "explanation an earlier build printed here"),
        ("track_scatter_bad",
         "the per-tier lunar positions scatter %s px about the fitted track, so "
         "the track is fitting the alignment error rather than the Moon"),
    ):
        _v = stats.get(_k)
        if _v:
            A("             : *** " + (_msg % _v))
    _af = stats.get("align_failed")
    if _af:
        A(f"             : *** CROSS-TIER ALIGNMENT FAILED — the link network is "
          f"inconsistent by {_af:.0f} px half-res "
          f"({2 * _af:.0f} px full-res) against a tolerance of "
          f"{stats.get('align_tolerance', 0):.0f} px. The links contradict each "
          f"other; the solved shifts are a compromise between impossible "
          f"constraints. The merged limb, its radius, the disc mask and every "
          f"radial filter are built on those shifts. This run's output is not "
          f"usable — check whether the long tiers hold any corona to correlate on.")
    _d = stats.get("limb_fit_disputed")
    _dc = stats.get("limb_fit_corrected_px")
    if _d and _dc:
        A(f"limb correction: the merged half-level fit ran {_d:.0f} px LARGER "
          f"than the tiers' own radius, far more than they disagree with each "
          f"other, so it was moved {_dc:+.1f} px onto their consensus with its "
          f"per-azimuth shape kept. The 50% crossing sits outside the true limb "
          f"whenever the edge is soft, and the merged edge is soft because the "
          f"Moon moves against the corona during the bracket; the per-tier fits "
          f"use the same method on unsmeared single tiers. R sets the radial "
          f"profile MGN divides out, FNRGF's rings, the deband and the disc "
          f"mask, so this circle had to be right before any of them.")
    elif _d:
        A(f"             : *** the merged limb fit disagrees with the tiers' own "
          f"radius by {_d:.0f} px, far more than they disagree with each other, "
          f"and it was NOT corrected — the fit has no per-azimuth profile to "
          f"move. R sets the radial profile MGN divides out, FNRGF's rings and "
          f"the deband, so a wrong circle prints concentric arcs in all of them "
          f"at once. Check the disc mask size on the preview.")
    ac = stats.get("autocrop_px") or {}
    if ac:
        A(f"alignment trim: {ac['top']}/{ac['bottom']} top/bottom, "
          f"{ac['left']}/{ac['right']} left/right removed — the border each "
          f"shift vacated; {ac['kept']} px kept")
    _mw = stats.get("merge_weight") or {}
    if _mw.get("alpha") is not None:
        _al = float(_mw["alpha"])
        if abs(_al - 1.0) < 1e-9:
            A(f"merge weight : exposure exponent 1.00 (unchanged) — no tilt "
              f"toward the short tiers measurably beat it")
        else:
            A(f"merge weight : exposure exponent {_al:.2f} — tilted toward the "
              f"shorter, sharper tiers, {100 * _mw.get('gain_limb', 0):+.0f}% "
              f"coherent detail at 1.02-1.12 R")
            A(f"             : accepted only because the mid and outer shells "
              f"kept their radial coherence; that guard is what separates "
              f"detail from grain")
    _ish = [max(q.get("intra_shift_px") or [0]) for q in (stats.get("quality") or {}).values()]
    if _ish and max(_ish) > 0:
        A(f"frame motion : up to {max(_ish):.0f} px between frames within a tier "
          f"(each frame is windowed on its own disc, so this is corrected, not "
          f"tolerated)")
    if stats.get("hot_pixels") is not None:
        A(f"sensor defects: {stats['hot_pixels']} photosites repaired")
    # THE TWO DECODE SETTINGS THE REPORT NEVER STATED (0.22.78). Both are in
    # stats and both change every pixel: the white balance is a SETTING the
    # user picks, and the saturation level is a measurement that has been wrong
    # before and is silent when it is. A report that omits them cannot be used
    # to explain a colour cast or a black frame after the fact.
    if stats.get("wb_used"):
        _g = stats.get("wb_gains") or []
        A(f"white balance: {stats['wb_used']}"
          + (f" — R {_g[0]:.4f} B {_g[2]:.4f}, G = 1" if len(_g) == 3 else ""))
    if stats.get("sat_level"):
        A(f"saturation   : {stats['sat_level']:.0f} ADU above black, from "
          f"{stats.get('sat_source', 'the reported white level')}")
    _fl = stats.get("flat") or {}
    if _fl.get("dir") and _fl.get("applied") and _fl.get("n_used"):
        A(f"flat field   : {_fl['combine']} from "
          f"{os.path.basename(_fl['dir'])}/, corrects a "
          f"{100 * (_fl.get('vignette', 1) - 1):.1f}% falloff")
        if "noise_master" in _fl:
            A(f"             : master flat noise {100 * _fl['noise_raw']:.3f}% "
              f"per photosite, {100 * _fl['noise_master']:.3f}% after a "
              f"{_fl.get('sigma_px', 0):.1f} px smooth — that is what the "
              f"division injects into every frame")
        if _fl.get("rejected"):
            A(f"             : {len(_fl['rejected'])} flat frame(s) rejected "
              f"({_fl['rejected'][0]['file']} {_fl['rejected'][0]['why']})")
    elif _fl.get("dir"):
        A(f"flat field   : NOT applied — {_fl.get('error', 'unavailable')}")

    A("")
    A("MEASURED")
    A("-" * 60)
    g = stats.get("geometry") or {}
    if g:
        R = g.get("R", 0)
        A(f"image        : {stats.get('W', '?')} x {stats.get('H', '?')} px")
        A(f"lunar limb   : centre ({g.get('cy', 0):.1f}, {g.get('cx', 0):.1f}) px, "
          f"R = {R:.1f} px  (diameter {2 * R:.1f} px)")
        if g.get("rms") is not None:
            A(f"limb fit     : {g['rms']:.2f} px rms over {g.get('rays_kept', '?')}"
              f"/{g.get('rays', '?')} rays; disc mask at R+{g.get('Rmask', R) - R:.1f} px")
        if R > 0 and stats.get("W"):
            scale = MOON_ARCSEC / (2 * R)
            A(f"plate scale  : ~{scale:.2f} arcsec/px  (assuming a {MOON_ARCSEC:.0f}\" "
              f"lunar disc, +/-{MOON_ARCSEC_TOL * 100:.0f}%)")
            fov_w = stats["W"] * scale / 3600.0
            fov_h = stats["H"] * scale / 3600.0
            A(f"field of view: {fov_w:.2f} x {fov_h:.2f} deg "
              f"= {stats['W'] / (2 * R):.1f} x {stats['H'] / (2 * R):.1f} lunar diameters")
            if cam.get("focal_mm"):
                pitch = cam["focal_mm"] * 1000.0 * scale / 206265.0
                A(f"implied pitch: {pitch:.2f} um/px at {cam['focal_mm']:g} mm "
                  f"(cross-check against your sensor)")
    if stats.get("corona_extent_R"):
        A(f"corona traced: out to {stats['corona_extent_R']:.1f} lunar radii "
          f"before the signal drops into the sky noise")
    _tr = stats.get("tone_reach")
    if _tr and _tr.get("contrast_at_extent") is not None:
        _pc = 100.0 * _tr["contrast_at_extent"]
        A(f"              : at that radius the render shows structure with "
          f"{_pc:.0f}% of the contrast it gives just above the limb — the "
          f"detail layers are flat in radius by construction, and Radial "
          f"flatten ({_tr['radial_flatten']:.2f}) sets how much of the "
          f"brightness gradient is put back on top of them. This is a look, "
          f"not an error: raise it to show the far field, lower it to keep the "
          f"inner corona's modelling")
    if stats.get("hdr_range_ev"):
        A(f"coronal range: {stats['hdr_range_ev']:.1f} EV between the inner corona "
          f"at the limb and the outer background")
    sg = stats.get("sky_gradient")
    if sg:
        if sg.get("applied"):
            pc = sg.get("per_channel")
            span = (f"R {pc[0]:.3f}x G {pc[1]:.3f}x B {pc[2]:.3f}x" if pc
                    else f"{sg['ratio']:.3f}x")
            A(f"sky gradient : {span} across the frame, removed per channel "
              f"({'quadratic' if sg.get('order', 1) >= 2 else 'plane'}, tilt "
              f"{sg['angle_deg']:+.0f} deg, fitted beyond "
              f"{sg['fitted_beyond_R']:.1f} R, {sg['sigma']:.0f} sigma)")
        elif sg.get("skipped_by_env"):
            A("sky gradient : SKIPPED by ECLIPSEFORGE_NO_SKYGRAD=1 — the "
              "per-channel division was not applied (diagnostic run)")
        elif sg.get("ratio") is not None:
            A(f"sky gradient : {sg['ratio']:.3f}x across the frame — below the "
              f"threshold, left alone")
        else:
            # The fit can bail out before it has a ratio to report (too little
            # sky beyond the corona, an implausible model). Say so rather than
            # raising KeyError at the end of a 3-minute run -- which is exactly
            # what the env switch above did on its first test.
            A("sky gradient : not applied")
    p = stats.get("prom") or {}
    if p:
        A(f"prominences  : corona colour R/GB = {p.get('med_red', 0):.2f}, "
          f"gate threshold {p.get('t0', 0):.2f}-{p.get('t1', 0):.2f}, "
          f"{p.get('area_px', 0)} px flagged")
        # Flagged and visible are different numbers, and the difference is the
        # answer to "where are my prominences?". Only said when it matters.
        if p.get("detail_layer"):
            A(f"             : a prominence detail layer was built from the "
              f"H-alpha tier's red channel — the corona filters cannot carry "
              f"this structure (MGN clips it and divides out its local sigma)")
        _vis = p.get("area_visible_px")
        if _vis is not None and _vis < p.get("area_px", 0):
            if _vis == 0:
                A(f"             : none of it outside the disc mask — the gate "
                  f"found redness only at or inside the limb, so the "
                  f"prominence slider has nothing to act on. Prominences may "
                  f"still be visible as brightness; this layer only adds the "
                  f"ones it can identify by colour")
            else:
                A(f"             : {_vis} px of that outside the disc mask; the "
                  f"rest sits under it and does not reach the picture")

    A("")
    A("PROCESSING")
    A("-" * 60)
    o = stats.get("options") or {}
    # SAY THE MODE FIRST, LOUDLY. The selectors below are printed as the user
    # left them, which is right -- but on the simple path most of them were
    # never consulted, and a reader who sees "correlation: semi-phase" against a
    # run that did plain phase correlation in simple.py has been told something
    # false by omission. One line at the top fixes that for all of them.
    _simple = stats.get("mode") == "simple stack"
    if _simple:
        A("MODE         : SIMPLE STACK — the fallback path")
        A("             : align, average each exposure group, measure the "
          "exposure ratios")
        A("             : from the pixels, merge. No dark, no flat, no "
          "hot-pixel repair, no")
        A("             : ladder solve, no LDIC, no feather, no tier "
          "projection.")
        A("             : The selectors marked [bypassed] below were NOT "
          "consulted.")
        A("")
    A(f"denoise      : {o.get('denoise', '?')}")
    A(f"earthshine   : {'on' if o.get('earthshine') else 'off'}")
    # EVERY SELECTOR THAT CHANGES THE RESULT, NAMED. Two of these used to be
    # readable only by inference -- the FNRGF preset from the order it logged,
    # the merge weight from its own line -- and the other three not at all. A
    # run whose settings cannot be read off its own report cannot be compared
    # with another run, which is the whole point of having selectors.
    _sel = [("correlation  ", "align_corr", "semi",
             {"cross": "cross-correlation",
              "semi": "semi-phase (amplitude floor)",
              "phase": "full phase (whitened)"}),
            ("align filter ", "align_filter", "isotropic",
             {"isotropic": "isotropic high-pass",
              "tangential": "tangential (radial removed)"}),
            ("tier combine ", "stack_combine", "mean",
             {"mean": "mean (no per-pixel rejection)",
              "clip": "clipped mean (kappa-sigma)"}),
            ("FNRGF        ", "fnrgf_preset", "ours",
             {"ours": "EFHDR (order 6, hard cutoff)",
              "published": "published (order 30, attenuated)"}),
            ("photom solve ", "photo_solve", "chain",
             {"chain": "chain (from the middle tier)",
              "network": "network (least squares)"}),
            ("tiers        ", "tier_mode", "exposure",
             {"exposure": "by shutter speed",
              "frame": "per frame (grouping ignored)"})]
    # Which of these the simple path actually consults. FNRGF does -- it is a
    # detail-layer setting and the detail layers run normally. The rest are
    # merge and alignment settings that simple.py does not read at all.
    _used_by_simple = {"fnrgf_preset"}
    for _lab, _key, _dflt, _names in _sel:
        _v = str(o.get(_key, _dflt))
        _note = "" if _v == _dflt else "   [not the default]"
        if _simple and _key not in _used_by_simple:
            _note = "   [bypassed — simple stack]"
        A(f"{_lab}: {_names.get(_v, _v)}{_note}")
    # NEITHER OF THESE HAPPENS ON AN IMPORT (0.22.78). An imported HDR has no
    # raw frames: nothing was despeckled and nothing was selected, but the
    # report printed "repaired" and "all frames per tier (best SNR)" all the
    # same -- stating as done two steps the run never reached.
    if _simple:
        A("hot pixels   : not repaired — the simple path does no defect "
          "correction")
        A("frames used  : every frame of every group, averaged")
        A("merge weight : n/a — a hat on the pixel value, no feather, not the "
          "toolbar setting")
        A("white balance: n/a — the simple path applies none")
    elif stats.get("mode") == "imported HDR":
        A("hot pixels   : n/a — an imported image has no raw frames to repair")
        A("frames used  : n/a — one finished image, imported")
    else:
        A(f"hot pixels   : "
          f"{'repaired' if o.get('despeckle', True) else 'left as shot'}")
        _fm = {"all": "all frames per tier (best SNR)",
               "best50": "sharpest half per tier",
               "best": "sharpest frame only (max detail)"}
        A(f"frames used  : "
          f"{_fm.get(o.get('frames', 'all'), o.get('frames', 'all'))}")

    A("")
    A("METHODS")
    A("-" * 60)
    # In import mode two of these would describe something that did not happen:
    # there are no tiers, so "a stack of the four shortest tiers" and "a single
    # fast tier" are not true of this run. The layers are still built, so the
    # entries stay -- restated to say what was actually done.
    _imp = stats.get("mode") == "imported HDR"
    _sim = stats.get("mode") == "simple stack"
    _skip = ({"Frame selection", "Sensor defects", "Alignment", "Prominence anchors",
              "Demosaic", "Photometric calibration", "HDR merge"} if _imp else set())
    # THE SIMPLE PATH RUNS NONE OF THESE EITHER, and unlike the selectors above
    # these are not settings -- they are claims about what the code did. Three
    # are replaced with what simple.py actually does; the rest are dropped.
    if _sim:
        _skip = {"Frame selection", "Sensor defects", "Prominence anchors",
                 "Photometric calibration"}
    _swap = {
        "Inner corona":
            "MGN of the imported image at coarser scales. With a bracket this is "
            "a separate stack of the shortest tiers, which see the inner corona "
            "unsaturated; here it is the same pixels through a different filter",
        "Prominences":
            "H-alpha redness R/((G+B)/2) measured on the imported image, "
            "thresholded against the robust spread of the corona's own colour in "
            "the limb annulus. With a bracket this reads a fast tier where the "
            "chromosphere is not blown, so an imported image flags less",
    } if _imp else {}
    if _sim:
        _swap = {
            "Alignment":
                "plain phase correlation of a half-resolution high-passed log "
                "green channel, every frame tied directly to the first frame of "
                "the shortest tier. No link network, no prominence anchors, no "
                "lag-2 solve -- deliberately sharing nothing with the normal "
                "path, so that this run is an independent check on it",
            "HDR merge":
                "weighted average of the tier means divided by exposure ratios "
                "MEASURED between adjacent tiers on pixels well exposed in both, "
                "never from the header. The weight is a hat on the pixel value "
                "alone -- zero at the noise floor, zero towards saturation -- "
                "with no exposure-time term, so a short tier carries no weight "
                "in the outer field where it holds only noise",
            "Demosaic":
                "a 3-plane FITS is read as planes, with no mosaic/demosaic round "
                "trip at all; a Bayer raw takes Malvar-He-Cutler as usual",
            "Inner corona":
                "MGN of the mean of the four shortest tiers, put on the merged "
                "scene's scale. Independent of the merge, as on the normal path",
        }
    # The Demosaic entry used to be a constant, so a run made with VNG still
    # reported Malvar-He-Cutler in its own methods section -- the report
    # describing a pipeline the run had not used.
    if (stats.get("options") or {}).get("demosaic") == "vng":
        _swap["Demosaic"] = (
            "VNG, Variable Number of Gradients (Chang, Cheng & Ward 1999), from "
            "the gradient tables in LibRaw's own misc_demosaic.cpp -- the same "
            "algorithm PixInsight's RAW module uses, so a render can be compared "
            "with one of its on equal terms. Agrees with LibRaw's output to 1.5 "
            "counts in 65535 over 44 million pixels once LibRaw's own 16-bit "
            "clip is applied. The clipping test is asked of the MOSAIC, so a "
            "pixel interpolated in part from a saturated photosite is flagged")
    for name, desc in METHODS:
        if name in _skip:
            continue
        A(f"* {name}: {_swap.get(name, desc)}")
    A("")
    A("All methods are published and their patents (where any existed) expired;")
    A("see README for citations.")
    return "\n".join(L)


def write(wd, stats):
    """Write report.json + report.txt into the workdir; return the text."""
    txt = build(stats)
    try:
        json.dump(stats, open(os.path.join(wd, "report.json"), "w"), indent=1, default=str)
        open(os.path.join(wd, "report.txt"), "w").write(txt + "\n")
    except Exception:
        pass
    return txt


def measure_image(lum, cy, cx, R):
    """Corona extent and achieved dynamic range, from the merged luminance."""
    out = {}
    try:
        H, W = lum.shape
        d = lum[::4, ::4]
        yy = np.arange(d.shape[0], dtype=np.float32)[:, None] - cy / 4
        xx = np.arange(d.shape[1], dtype=np.float32)[None, :] - cx / 4
        r = np.hypot(yy, xx)
        Rq = R / 4
        ri = r.astype(np.int32)
        n = int(r.max()) + 1
        cnt = np.bincount(ri.ravel(), minlength=n)[:n]
        s = np.bincount(ri.ravel(), weights=d.ravel().astype(np.float64), minlength=n)[:n]
        prof = s / np.maximum(cnt, 1)
        sky = (r > 0.8 * r.max())
        bg = float(np.median(d[sky]))
        noise = 1.4826 * float(np.median(np.abs(d[sky] - bg)))
        thr = bg + 3 * noise
        idx = np.flatnonzero((prof > thr) & (np.arange(n) > Rq))
        if len(idx):
            k = idx.max()
            # first radius beyond the limb where it drops for good
            drop = np.flatnonzero(prof[int(Rq):] <= thr)
            k = int(Rq) + int(drop[0]) if len(drop) else k
            out["corona_extent_R"] = float(k / Rq)
        # HOW MUCH OF THAT REACH THE TONE CURVE ACTUALLY SHOWS.
        #
        # The extent above is a property of the DATA. What the viewer sees is
        # the data times the render's envelope, and render.py builds the picture
        # as a PRODUCT: Y = B * (baseLift + detailGain * det), with
        # B = bg / rprof ** radialFlatten. The detail layers leave MGN and FNRGF
        # with roughly radius-independent amplitude -- that is what a
        # normalising filter is for, and FNRGF's author says so outright: "After
        # the processing, all circles have identical mean pixel value and
        # identical standard deviation of pixel values" -- so multiplying by an
        # envelope that still falls with radius puts the falloff straight back
        # into the rendered contrast.
        #
        # `radialFlatten` at its 0.5 default removes half the gradient in the
        # exponent, which leaves the other half in the picture. Measured on the
        # reference bracket, display contrast per unit of detail amplitude,
        # normalised to the innermost shell:
        #
        #     shell      1.2R   1.5R   2.0R   2.6R   3.2R   3.9R
        #     rf 0.50    0.513  0.225  0.140  0.118  0.108  0.103
        #     rf 0.75    0.716  0.474  0.375  0.344  0.329  0.320
        #     rf 1.00    1.000  1.000  1.000  1.000  1.000  1.000
        #
        # -- and on the same bracket the signal at 3.2-3.9 R is still 8x the sky
        # noise. So structure that is plainly in the data is being shown at a
        # tenth of the contrast it gets at the limb, which is why the corona
        # "stops" well inside where it was traced.
        #
        # This is a TONE CHOICE and it is not made here: a flat envelope shows
        # the far field at the price of the inner corona's modelling, and which
        # is wanted is a matter of taste and of the picture. What was missing is
        # that the choice was invisible. The number below makes it a decision.
        # The report is written at stack time, before any tone has been picked,
        # so this is against the DEFAULT -- which is the number that matters,
        # since it is what the first look at the picture will have used.
        try:
            from .render import DEFAULTS as _RD
            _rf = float(_RD.get("radialFlatten", 0.5))
        except Exception:
            _rf = 0.5
        if out.get("corona_extent_R", 0) > 1.5 and prof[int(Rq)] > 0:
            _rr = np.arange(n, dtype=np.float64)
            _in = max(int(Rq * 1.05), int(Rq) + 1)
            _out_r = min(int(Rq * out["corona_extent_R"]), n - 1)
            if _out_r > _in and prof[_in] > bg and prof[_out_r] > bg:
                # B falls as (bg/rprof)^rf, and with bg == rprof by construction
                # that is rprof^(1-rf): the SURVIVING fraction of the gradient.
                _ratio = ((prof[_out_r] - bg) / (prof[_in] - bg)) ** (1.0 - _rf)
                out["tone_reach"] = {
                    "radial_flatten": _rf,
                    "contrast_at_extent": float(max(_ratio, 0.0)),
                    "extent_R": float(out["corona_extent_R"])}
        inner = (r > Rq + 0.5) & (r < Rq + 2.5)
        if inner.any() and bg > 0:
            peak = float(np.median(d[inner]))
            if peak > bg:
                out["hdr_range_ev"] = float(np.log2(peak / max(bg, 1e-6)))
    except Exception:
        pass
    return out
