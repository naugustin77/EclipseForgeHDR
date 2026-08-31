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
     "phase correlation on gradient-flattened log corona (Druckmuller 2009, ApJ 706, 1605); "
     "lag-1 and lag-2 links solved together by weighted least squares"),
    ("Prominence anchors",
     "the fastest tiers hold almost no corona to correlate on, so they are also tied in "
     "by normalized cross-correlation of prominence patches -- solar features, unlike the "
     "lunar limb, which drifts against the corona during totality. Anchors are found as "
     "compact azimuthal peaks on the limb; a link is used only while its anchors agree "
     "with each other, and enters the same least-squares network as an extra link"),
    ("Demosaic",
     "Malvar-He-Cutler gradient-corrected bilinear interpolation"),
    ("Photometric calibration",
     "per-tier scale factors from median ratios of overlapping unsaturated signal, "
     "chained to the middle tier"),
    ("HDR merge",
     "saturation-weighted linear-colour merge (Debevec-style weighting) in camera "
     "native colour, then to sRGB primaries"),
    ("Limb fit",
     "50% crossing between disc and near-limb corona level along 720 rays, "
     "robust fit of r = R + dx cos t + dy sin t, re-centred over three iterations"),
    ("MGN detail",
     "Multiscale Gaussian Normalization (Morgan & Druckmuller 2014, Sol. Phys. 289, 2945), "
     "photon-noise-adaptive, disc excluded by normalized convolution, radial profile "
     "removed before normalization"),
    ("FNRGF detail",
     "Fourier Normalizing Radial Gradient Filter (Morgan, Habbal & Woo 2006; "
     "Druckmullerova et al. 2011, ApJS 194, 25) with coverage-matched Fourier order "
     "and Huber-weighted robust fitting"),
    ("Pellett layer",
     "tangential (rotational) unsharp about the disc centre, normalized convolution "
     "in polar space"),
    ("Inner corona",
     "separate MGN of a saturation-weighted stack of the four shortest tiers"),
    ("Prominences",
     "H-alpha redness R/((G+B)/2) in a single fast tier, thresholded against the "
     "robust spread of the corona's own colour in the limb annulus"),
    ("NAFE detail",
     "Noise Adaptive Fuzzy Equalization with a Variable Neighbourhood (Druckmuller 2013, "
     "ApJ 775, 88; Druckmuller & Druckmullerova, IWCIA 2014, LNCS 8466, 262): each pixel is "
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

    A("")
    A("EXPOSURE STACK")
    A("-" * 60)
    tiers = stats.get("tiers", [])
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
              f"(lower = the tiers agree)")
        if "rim_width_px" in aq:
            A(f"             : disagreement rim {aq['rim_width_px']:.0f} px wide "
              f"just outside the limb")
        if "limb_width_med" in aq:
            A(f"             : merged limb 20-80% transition "
              f"{aq['limb_width_med']:.1f} px (p90 "
              f"{aq.get('limb_width_p90', float('nan')):.1f} px)")
    ac = stats.get("autocrop_px") or {}
    if ac:
        A(f"alignment trim: {ac['top']}/{ac['bottom']} top/bottom, "
          f"{ac['left']}/{ac['right']} left/right removed — the border each "
          f"shift vacated; {ac['kept']} px kept")
    _ish = [max(q.get("intra_shift_px") or [0]) for q in (stats.get("quality") or {}).values()]
    if _ish and max(_ish) > 0:
        A(f"frame motion : up to {max(_ish):.0f} px between frames within a tier "
          f"(each frame is windowed on its own disc, so this is corrected, not "
          f"tolerated)")
    if stats.get("hot_pixels") is not None:
        A(f"sensor defects: {stats['hot_pixels']} photosites repaired")
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
        else:
            A(f"sky gradient : {sg['ratio']:.3f}x across the frame — below the "
              f"threshold, left alone")
    p = stats.get("prom") or {}
    if p:
        A(f"prominences  : corona colour R/GB = {p.get('med_red', 0):.2f}, "
          f"gate threshold {p.get('t0', 0):.2f}-{p.get('t1', 0):.2f}, "
          f"{p.get('area_px', 0)} px flagged")

    A("")
    A("PROCESSING")
    A("-" * 60)
    o = stats.get("options") or {}
    A(f"denoise      : {o.get('denoise', '?')}")
    A(f"earthshine   : {'on' if o.get('earthshine') else 'off'}")
    A(f"hot pixels   : {'repaired' if o.get('despeckle', True) else 'left as shot'}")
    _fm = {"all": "all frames per tier (best SNR)",
           "best50": "sharpest half per tier",
           "best": "sharpest frame only (max detail)"}
    A(f"frames used  : {_fm.get(o.get('frames', 'all'), o.get('frames', 'all'))}")

    A("")
    A("METHODS")
    A("-" * 60)
    for name, desc in METHODS:
        A(f"* {name}: {desc}")
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
        inner = (r > Rq + 0.5) & (r < Rq + 2.5)
        if inner.any() and bg > 0:
            peak = float(np.median(d[inner]))
            if peak > bg:
                out["hdr_range_ev"] = float(np.log2(peak / max(bg, 1e-6)))
    except Exception:
        pass
    return out
