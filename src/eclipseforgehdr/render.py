"""Recombination render from cached layers, preview generation, exports."""
from __future__ import annotations
import os, json
import numpy as np
from scipy import ndimage
from PIL import Image

# Starting points, not truths. These are the settings the reference bracket was
# worked to by eye once the layers were behaving (0.12.0); every one of them is
# a taste call and every one is a slider.
DEFAULTS = {
    "fnMix": 0.17, "detailGain": 1.0, "baseLift": 0.255, "envGamma": 1.0,
    "radialFlatten": 0.5, "mgnContrast": 0.04, "fnCompress": 0.9,
    "clarity": 0.39, "smoothing": 0.25, "pelGain": 0.0,
    "nafeMix": 0.15, "innerMix": 0.31, "innerDenoise": 0.4, "innerDim": 0.05, "promGain": 0.4,
    "temp": 0.9, "tint": 1.205, "bgNeutral": 1.0, "satur": 1.0, "hlCompress": 0.1, "hlDesat": 0.0,
    "outGamma": 1.0, "bgBlack": 0.005,
    "discLevel": 0.045, "discTrim": 0.0, "earthShine": 0.0,
    "ringBlend": 0.0, "ringScale": 1.0, "ringDX": 0.0, "ringDY": 0.0,
}



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


class Layers:
    """Cached full-res layers, plus a decimated copy for interactive previews."""

    def __init__(self, wd, preview_decim=4):
        self.wd = wd
        geo = json.load(open(os.path.join(wd, "geometry.json")))
        self.cy, self.cx, self.R = geo["cy"], geo["cx"], geo["R"]
        self.Rmask = float(geo.get("Rmask", self.R + 4.0))
        self.limb_prof = geo.get("limb_prof")
        self.limb_margin = float(geo.get("limb_margin", self.Rmask - self.R))
        lum = np.load(os.path.join(wd, "hdr_lum.npy"))
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
        hi = float(np.percentile(lum, 99.97))
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
            lo = float(np.percentile(lum, 2)); hi = float(np.percentile(lum, 99.97))
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
        rprof = np.maximum(prof[rid], 0.12).astype(np.float32)

        mg = np.load(os.path.join(wd, "mgn.npy"))
        mg = np.clip((mg - np.percentile(mg, 1)) /
                     (np.percentile(mg, 99.7) - np.percentile(mg, 1)), 0, 1)
        D = np.load(os.path.join(wd, "fnrgf.npy"))
        nfp = os.path.join(wd, "nafe.npy")
        nfl = np.load(nfp) if os.path.exists(nfp) else np.full_like(D, 0.5)
        fnl = np.clip(D / 8.0 + 0.5, 0, 1)     # raw sigma units, ±4σ window
        del D
        inner = np.load(os.path.join(wd, "inner.npy"))
        inner0 = np.load(os.path.join(wd, "inner0.npy"))
        ep = os.path.join(wd, "earth.npy")
        self.has_earth = os.path.exists(ep)
        earth = np.load(ep) if self.has_earth else np.full(lum.shape, 0.5, np.float32)
        gate = np.load(os.path.join(wd, "prom.npy"))
        pelp = os.path.join(wd, "pellett.npy")
        pel = np.load(pelp) if os.path.exists(pelp) else np.full(lum.shape, 0.5, np.float32)
        hdr = np.load(os.path.join(wd, "hdr_rgb.npy"), mmap_mode="r")
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
        _conf = np.clip((Ls - _floor) / (4.0 * _nz), 0.0, 1.0).astype(np.float32)
        ratio = np.empty(lum.shape + (3,), np.float32)
        for c in range(3):
            rc = ndimage.gaussian_filter(
                np.ascontiguousarray(hdr[:, :, c]), 6) / np.maximum(Ls, _floor)
            ratio[:, :, c] = 1.0 + _conf * (rc - 1.0)
            del rc
        self.colour_floor = float(_floor)
        self.colour_conf_frac = float((_conf > 0.5).mean())
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
        self.full = {"bg": Bg, "rprof": rprof, "mgn": mg, "fnrgf": fnl,
                     "inner": inner, "inner0": inner0, "earth": earth,
                     "prom": gate, "pel": pel, "ratio": ratio,
                     "nafe": nfl, "cconf": self._cconf}
        q = preview_decim
        self.prev = {k: v[::q, ::q] for k, v in self.full.items()}
        # clarity/smoothing variants (preview scale; full-res computed on demand)
        for key in ("mgn", "fnrgf", "nafe"):
            f = _fill_disc(self.prev[key], self.cy / q, self.cx / q, self.Rmask / q)
            self.prev[key + "_lo"] = ndimage.gaussian_filter(f, 8.0 / q * 2)
            self.prev[key + "_sm"] = ndimage.gaussian_filter(f, 0.6)
            del f
        self.prev_decim = q
        _geom_for_fill["g"] = (self.cy, self.cx, self.Rmask)   # full-res variants
        self.reload_contact()

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


def _detail_layers(src, P, preview=True):
    """Runtime-transformed detail layers (shared by composite render and layer views)."""
    mgr = src["mgn"]
    if P["clarity"] > 0:
        mgr = mgr + P["clarity"] * (mgr - _variant(src, "mgn", "lo", preview))
    if P["smoothing"] > 0:
        mgr = (1 - P["smoothing"]) * mgr + P["smoothing"] * _variant(src, "mgn", "sm", preview)
    mg = np.clip(0.5 + (mgr - 0.5) * P["mgnContrast"], 0, 1)
    fnr = src["fnrgf"]
    if P["clarity"] > 0:
        fnr = fnr + P["clarity"] * (fnr - _variant(src, "fnrgf", "lo", preview))
    if P["smoothing"] > 0:
        fnr = (1 - P["smoothing"]) * fnr + P["smoothing"] * _variant(src, "fnrgf", "sm", preview)
    D = (fnr - 0.5) * 8.0                     # back to sigma units
    s_hi = 2.5 / max(P["fnCompress"], 1e-4)   # 0 = FNRGF off (flat 0.5)
    fn = (np.where(D >= 0, np.tanh(D / s_hi), np.tanh(D / (1.5 * s_hi))) + 1) / 2
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
    if nf_ is not None and P.get("nafeMix", 0.0) > 0:
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
    return mg, fn, inner_eff


def render(layers: Layers, params, preview=False, view="composite"):
    P = dict(DEFAULTS); P.update(params or {})
    src = layers.prev if preview else layers.full
    decim = layers.prev_decim if preview else 1
    cy, cx, R = layers.geometry(decim)
    H, W = src["bg"].shape
    mg, fn, inner_eff = _detail_layers(src, P, preview=preview)
    if view == "mgn":
        return np.repeat(mg[:, :, None], 3, axis=2)
    if view == "fnrgf":
        return np.repeat(fn[:, :, None], 3, axis=2)
    if view == "inner":
        return np.repeat(inner_eff[:, :, None], 3, axis=2)
    if view == "prom":
        return np.repeat(src["prom"][:, :, None], 3, axis=2)
    if view == "pellett":
        return np.repeat(src["pel"][:, :, None], 3, axis=2)
    if view == "nafe":
        return np.repeat(np.clip(src.get("nafe", np.full_like(mg, 0.5)), 0, 1)[:, :, None], 3, axis=2)

    yy = np.arange(H, dtype=np.float32)[:, None] - cy
    xx = np.arange(W, dtype=np.float32)[None, :] - cx
    r = np.sqrt(yy * yy + xx * xx)
    Re = layers.mask_radius_map((H, W), decim) + P["discTrim"] / decim
    edge = np.clip((r - (Re - 10 / decim)) / (12.0 / decim), 0, 1)
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

    B = (src["bg"] ** P["envGamma"]) * (1 - P["innerDim"] * wG)
    B = B / (src["rprof"] ** P["radialFlatten"])
    B *= (1 - P["radialFlatten"] * 0.5)       # keep overall level roughly stable
    det = (1 - P["fnMix"] * wf) * mg + P["fnMix"] * wf * fn
    det = det + P["pelGain"] * (src["pel"] - 0.5)
    det = det * (1 - P["innerMix"] * wI) + P["innerMix"] * wI * inner_eff
    Y = B * (P["baseLift"] + P["detailGain"] * det)
    # prominence: local-contrast modulation inside the gate, with a small
    # positive bias so a detected prominence gains presence, not just texture
    Y = Y * (1 + P["promGain"] * src["prom"] * (0.30 + 1.3 * (inner_eff - 0.5)))
    Yd = P["discLevel"] * (1 + P["earthShine"] * (2 * src["earth"] - 1))
    Y = Y * edge + Yd * (1 - edge)
    del B, det, inner_eff, Yd, wf, wI, wG, mg, fn

    a = np.clip(src["ratio"], 0.2, 3.0) ** P["satur"]
    if P.get("bgNeutral", 0) > 0:
        # weighted by the same confidence that built `ratio` -- see the note
        # where _cconf is stored
        _bn = (layers.bg_chroma[None, None, :] ** P["bgNeutral"]) - 1.0
        a = a / (1.0 + src["cconf"][:, :, None] * _bn)
    a[:, :, 0] *= P["temp"]
    a[:, :, 2] /= P["temp"]
    a[:, :, 1] *= P.get("tint", 1.0)
    # renormalise to unit luminance so colour moves never change brightness
    al = (0.2126 * a[:, :, 0] + 0.7152 * a[:, :, 1] + 0.0722 * a[:, :, 2])
    a /= np.maximum(al, 1e-6)[:, :, None]
    del al
    a = a * edge[:, :, None] + (1 - edge)[:, :, None]   # neutral moon disc
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
    return rgb


def export(layers: Layers, params, fmt, out_path, view="composite", size="full"):
    """fmt: tif16 | tif8 | png (16-bit when OpenCV present, else 8-bit) | jpg.
    view: composite | mgn | fnrgf | nafe | inner | prom | pellett (detail views export grayscale).
    size: full | half (half = 2x2 binned, ~2x better SNR)."""
    rgb = render(layers, params, preview=False, view=view)
    if size == "half":
        H2, W2 = rgb.shape[0] // 2 * 2, rgb.shape[1] // 2 * 2
        rgb = rgb[:H2, :W2].reshape(H2 // 2, 2, W2 // 2, 2, -1).mean(axis=(1, 3))
        if rgb.shape[-1] == 1:
            rgb = rgb[:, :, 0]
        rgb = np.ascontiguousarray(rgb)
    gray = view != "composite"
    arr16 = (rgb[:, :, 0] if gray else rgb)
    # Tag the file. The render is already in display encoding (the browser
    # preview shows these exact bytes as sRGB), so sRGB is what it IS -- but it
    # was going out untagged, leaving every host to guess. Photoshop guessed
    # right by accident; PixInsight assumes linear and got it wrong.
    from . import icc
    _prof = icc.srgb_profile()
    _tags = [(34675, 1, len(_prof), _prof, False)]
    if fmt == "tif16":
        import tifffile
        tifffile.imwrite(out_path, (arr16 * 65535 + 0.5).astype(np.uint16),
                         compression="zlib", extratags=_tags,
                         photometric="minisblack" if gray else "rgb",
                         description="eclipseforgehdr %s, params: %s" % (view, json.dumps(params)))
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
        except Exception:
            Image.fromarray((arr16 * 255 + 0.5).astype(np.uint8)).save(out_path)
    elif fmt == "jpg":
        Image.fromarray((arr16 * 255 + 0.5).astype(np.uint8)).save(
            out_path, quality=92, optimize=True, icc_profile=_prof)
    else:
        raise ValueError(fmt)
    return out_path
