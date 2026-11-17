"""THE MOON LAYER (0.24): earthshine structure on the lunar disc.

Built from the long frames of the bracket, after the stack, in about a minute,
and cached as MOON_FILE on the merge grid (0.5 = neutral, the format the old
earth.npy had, so render.py's Earthshine control drives it unchanged).

Recipe (derived and measured
on the reference test set, Oct 2026):
  1. every long frame read from the raw file, calibrated with the cached master
     bias / dark / flat when they fit, superpixel green;
  2. each frame's Moon found by its own limb and the frames stacked ON THE MOON
     (the Moon moves against the Sun and the sensor during totality);
  3. the scattered light on the disc modelled as the merge's corona OUTSIDE the
     Moon convolved with a power-law PSF wing (exponent fitted), plus a smooth
     sky term, fitted on the disc and subtracted;
  4. anything that follows the limb taken off (the median of every ring): no
     lunar feature follows the Moon's own edge, the glare does;
  5. band-passed to the maria scales: the finest scale (a few px) is fixed to
     the SENSOR on this set (measured: early vs late frames correlate 0.41 in
     sensor coordinates, 0.23 on the Moon), so it is left out;
  6. faded out over the last 0.08 R, where the limb's spill owns the light.
Checked against an independent earthshine processing of the same frames: same
orientation; correlation 0.54 / 0.67 / 0.61 / 0.27 in the bands 0.02-0.08 /
0.03-0.15 / 0.06-0.3 / 0.1-0.5 R (the largest is where the two scattered-light
models differ most).
"""
import os, json
import numpy as np
from scipy import ndimage
from .parallel import gaussian_filter as _pgauss   # scipy's, on several threads
from scipy.signal import fftconvolve

MOON_FILE = "moon_layer.npy"
MOON_META = "moon_layer.json"
MOON_BUILD = 9              # 2: limb band levelled and equalised before the fade (2026-10-07)
                            # 3: scaled by the maria agreement of two half-stacks (2026-10-28)
                            # 4: the structure checked against a map of the real Moon (lunarcheck.py, 2026-11-14)
                            # 5: the Force earthshine variant built alongside (2026-11-14)
                            # 9: blobs of the edge pass faded to the plain disc (2026-11-17)
#: the variant the Moon group's "Force earthshine" box switches to, and the
#: small file that remembers the box for this folder
MOON_FORCE_FILE = "moon_layer_forced.npy"
FORCE_STATE = "moon_force.json"
#: below this agreement between the two half-stacks the structure is taken as
#: noise or glare and the disc is left smooth
MARIA_MIN_AGREE = 0.10


def _log(progress, msg, frac=None):
    if progress is not None:
        progress.log(msg, frac)


def forced_on(wd):
    """Is Force earthshine ticked for this folder?"""
    try:
        return bool(json.load(open(os.path.join(wd, FORCE_STATE))).get("on", False))
    except Exception:
        return False


def set_forced(wd, on):
    json.dump({"on": bool(on)}, open(os.path.join(wd, FORCE_STATE), "w"))


def layer_path(wd):
    """The Moon layer to load: the forced variant when the box is ticked and
    one was built, the normal layer otherwise."""
    fp = os.path.join(wd, MOON_FORCE_FILE)
    if forced_on(wd) and os.path.exists(fp):
        return fp
    return os.path.join(wd, MOON_FILE)


def force_info(wd):
    """What the box can do here, for the page: available, on, source, note."""
    try:
        f = dict(json.load(open(os.path.join(wd, MOON_META))).get("force") or {})
    except Exception:
        f = {}
    f["available"] = bool(f.get("available")) and os.path.exists(os.path.join(wd, MOON_FORCE_FILE))
    f["on"] = forced_on(wd)
    return f


MOON_DUMP = "moon_stack_dump.npz"       # lab diagnostics (ECLIPSEFORGE_MOON_DUMP=1)


def needs_build(wd):
    """Missing, built by an older recipe, or older than the merge it sits on."""
    try:
        if (os.environ.get("ECLIPSEFORGE_MOON_DUMP", "") not in ("", "0")
                and not os.path.exists(os.path.join(wd, MOON_DUMP))):
            return True
        mp = os.path.join(wd, MOON_FILE)
        if not os.path.exists(mp):
            return True
        if int(json.load(open(os.path.join(wd, MOON_META))).get("build", 0)) < MOON_BUILD:
            return True
        hp = os.path.join(wd, "hdr_rgb.npy")
        return os.path.exists(hp) and os.path.getmtime(mp) < os.path.getmtime(hp)
    except Exception:
        return True


def _fit_limb(G, cy, cx, R0, iters=3):
    """Circle through the steepest rise of log G along 360 rays near R0."""
    L = np.log(np.maximum(_pgauss(G, 2.0), 1e-3))
    for _ in range(iters):
        pts = []
        rs = np.arange(0.93 * R0, 1.07 * R0, 0.25)
        for t in np.radians(np.arange(0, 360, 1.0)):
            v = ndimage.map_coordinates(L, [cy + rs * np.sin(t), cx + rs * np.cos(t)], order=1, mode="nearest")
            gr = np.gradient(v); i = int(np.argmax(gr))
            if 1 <= i < len(gr) - 1:
                a, b, c = gr[i - 1], gr[i], gr[i + 1]
                den = a - 2 * b + c
                d = 0.5 * (a - c) / den if den != 0 else 0.0
                rr = rs[i] + d * 0.25
                pts.append((cy + rr * np.sin(t), cx + rr * np.cos(t)))
        pts = np.array(pts)
        if len(pts) < 30:
            break
        for _ in range(4):
            A = np.c_[2 * pts[:, 1], 2 * pts[:, 0], np.ones(len(pts))]
            bb = pts[:, 1] ** 2 + pts[:, 0] ** 2
            s = np.linalg.lstsq(A, bb, rcond=None)[0]
            xc, yc = s[0], s[1]; R = float(np.sqrt(s[2] + xc * xc + yc * yc))
            res = np.hypot(pts[:, 1] - xc, pts[:, 0] - yc) - R
            k = np.abs(res) < max(3 * 1.4826 * np.median(np.abs(res)), 0.3)
            pts = pts[k]
        cy, cx, R0 = float(yc), float(xc), R
    return cy, cx, R0


def _masters(wd, shape):
    """Cached master bias / dark rate / flat, each only if it fits the frame."""
    bias = rate = flat = None
    try:
        z = np.load(os.path.join(wd, "masterdark.npz"))
        if "bias" in z.files:
            _b = z["bias"]          # NpzFile inflates on every access: read once
            if _b.shape == shape:
                bias = _b
        if "rate" in z.files:
            _r = z["rate"]
            if _r.shape == shape:
                rate = _r
    except Exception:
        pass
    try:
        f = np.load(os.path.join(wd, "masterflat.npy"), mmap_mode="r")
        if f.shape == shape:
            flat = f
    except Exception:
        pass
    return bias, rate, flat


def build(folder, wd, progress=None, frames=None):
    """Build MOON_FILE / MOON_META. Returns the meta dict, or None."""
    _mcache = {}
    from .raw import open_frame, list_raws, read_exif
    geo = json.load(open(os.path.join(wd, "geometry.json")))
    cy, cx, R = float(geo["cy"]), float(geo["cx"]), float(geo["R"])
    oy, ox = (geo.get("crop_origin") or [0, 0])[:2]
    # the GREEN channel of the merge: the stack below is green, and luminance
    # carries the prominences' H-alpha red, whose scatter the green frames do
    # not see -- with luminance the model over-subtracted under the big
    # prominence (a dark patch at the limb)
    _hr = np.load(os.path.join(wd, "hdr_rgb.npy"), mmap_mode="r")
    H, W = _hr.shape[:2]
    # ---- the frames: the long end of the bracket, where the disc has signal
    if frames is None:
        frames = []
        for p in list_raws(folder):
            try:
                frames.append((p, float(read_exif(p)[0])))
            except Exception:
                continue
    if not frames:
        return None
    tmax = max(t for _, t in frames)
    use = [(p, t) for p, t in frames if t >= tmax / 16.0 - 1e-9]
    # THE LONGEST EXPOSURE IS LEFT OUT when at least two shorter ones remain.
    # Measured on the reference set against that independent earthshine: its large-scale
    # pattern is instrumental (4x the others' level in rate units), and with it
    # the maria scale agrees 0.05 with that result, without it 0.27 (0.61 at
    # 0.04-0.2 R instead of 0.53).
    _shorter = sorted(set(t for _, t in use if t < tmax - 1e-9))
    if len(_shorter) >= 2:
        use = [(p, t) for p, t in use if t < tmax - 1e-9]
    if not use:
        return None
    _log(progress, f"Moon layer: stacking {len(use)} long frame(s) "
                   f"({min(t for _, t in use):g}-{max(t for _, t in use):g} s"
                   + (f"; longest {tmax:g} s left out" if len(_shorter) >= 2 else "")
                   + ")...", None)
    # ---- the window: a box about the Moon in frame coordinates, superpixel grid
    Rh = R / 2.0
    half = int(np.ceil(1.35 * R / 2.0)) * 2
    acc = None; wsum = None; ref = None; nused = 0
    hacc = [None, None]; hw = [None, None]      # the two halves (see MARIA AGREEMENT)
    _hn = [0, 0]
    for _fi, (p, t) in enumerate(use):
        try:
            rf = open_frame(p)
            bay = rf.bayer; Hb, Wb = bay.shape
            sat = float(rf.sat_level)
            gy, gx = cy + oy, cx + ox                # the merge's Moon, in frame px
            Y0 = int(gy - half) - int(gy - half) % 2; X0 = int(gx - half) - int(gx - half) % 2
            Y1, X1 = Y0 + 2 * half, X0 + 2 * half
            cfa = np.zeros((2 * half, 2 * half), np.float32)
            yA, yB, xA, xB = max(Y0, 0), min(Y1, Hb), max(X0, 0), min(X1, Wb)
            cfa[yA - Y0:yB - Y0, xA - X0:xB - X0] = bay[yA:yB, xA:xB]
            clip = cfa >= 0.9 * sat
            # the masters once per build, not per frame: each call inflated
            # the 44 MP compressed dark four times (audit 2026-10-10, ~5 s/frame)
            if _mcache.get("shape") != (Hb, Wb):
                _mcache.update(shape=(Hb, Wb), m=_masters(wd, (Hb, Wb)))
            bias, rate, flat = _mcache["m"]
            if bias is not None:
                cfa[yA - Y0:yB - Y0, xA - X0:xB - X0] -= bias[yA:yB, xA:xB]
            if rate is not None:
                cfa[yA - Y0:yB - Y0, xA - X0:xB - X0] -= rate[yA:yB, xA:xB] * np.float32(t)
            if flat is not None:
                cfa[yA - Y0:yB - Y0, xA - X0:xB - X0] /= np.maximum(np.asarray(flat[yA:yB, xA:xB], np.float32), 1e-3)
            del rf, bay
        except Exception as e:
            _log(progress, f"[warn] Moon layer: {os.path.basename(p)} skipped ({e})", None)
            continue
        G = 0.5 * (cfa[0::2, 1::2] + cfa[1::2, 0::2])
        cl = clip.reshape(half, 2, half, 2).any(axis=(1, 3))
        c0 = (gy - Y0) / 2.0 - 0.25; c1 = (gx - X0) / 2.0 - 0.25
        try:
            fy, fx, fR = _fit_limb(G, c0, c1, Rh)
        except Exception:
            continue
        if not (0.9 * Rh < fR < 1.1 * Rh) or np.hypot(fy - c0, fx - c1) > 0.2 * Rh:
            _log(progress, f"[warn] Moon layer: {os.path.basename(p)} skipped, limb not found", None)
            continue
        if ref is None:
            ref = (fy, fx, fR, Y0, X0)
        dy, dx = ref[0] - fy + (Y0 - ref[3]) / 2.0, ref[1] - fx + (X0 - ref[4]) / 2.0
        Gs = ndimage.shift(G, (dy, dx), order=3, mode="nearest") / np.float32(t)
        ok = ndimage.shift(cl.astype(np.float32), (dy, dx), order=1, mode="constant", cval=1.0) < 0.01
        w = np.float32(t) * ok
        acc = Gs * w if acc is None else acc + Gs * w
        wsum = w if wsum is None else wsum + w
        _h = 0 if _fi < len(use) / 2.0 else 1
        hacc[_h] = Gs * w if hacc[_h] is None else hacc[_h] + Gs * w
        hw[_h] = w if hw[_h] is None else hw[_h] + w
        _hn[_h] += 1
        nused += 1
    if acc is None or nused == 0:
        _log(progress, "[warn] Moon layer: no usable frame", None)
        return None
    M = acc / np.maximum(wsum, 1e-9)
    mcy, mcx, mR = ref[0], ref[1], ref[2]

    def _layer(M):
        """Steps 3-6 on one stacked disc: the band-passed, faded structure Eb
        (before the display normalisation), plus what the meta records."""
        n = M.shape[0]
        yy, xx = np.mgrid[:n, :n].astype(np.float64)
        r = np.hypot(yy - mcy, xx - mcx) / mR
        # ---- the light outside the Moon, from the merge, on this grid (the grid's
        # Moon placed on the merge's Moon)
        my = cy + (yy - mcy) * 2.0
        mx = cx + (xx - mcx) * 2.0
        _y0, _y1 = max(int(my.min()) - 2, 0), min(int(my.max()) + 3, H)
        _x0, _x1 = max(int(mx.min()) - 2, 0), min(int(mx.max()) + 3, W)
        Gm = np.asarray(_hr[_y0:_y1, _x0:_x1, 1], np.float32)
        I = ndimage.map_coordinates(Gm, [my - _y0, mx - _x0], order=1, mode="constant", cval=0.0)
        del Gm
        I[r < 1.005] = 0.0
        # ---- 3: the scattered light, fitted inside 0.97 R
        sub = np.flatnonzero((r < 0.97).ravel())[::2]
        X = (xx - mcx) / mR; Yv = (yy - mcy) / mR
        ky, kx = np.mgrid[-n:n + 1, -n:n + 1].astype(np.float64)
        rho2 = ky * ky + kx * kx
        best = None
        for alpha in (2.0, 2.5, 3.0):
            S = fftconvolve(I, (rho2 + 4.0) ** (-alpha / 2), mode="same")
            A = np.stack([S.ravel()[sub], np.ones(len(sub)), X.ravel()[sub], Yv.ravel()[sub],
                          (X * X + Yv * Yv).ravel()[sub]], -1)
            y = M.ravel()[sub]; keep = np.ones(len(y), bool)
            for _ in range(3):
                sc = np.abs(A[keep]).max(0) + 1e-30
                coef = np.linalg.lstsq(A[keep] / sc, y[keep], rcond=None)[0] / sc
                res = y - A @ coef; s = 1.4826 * np.median(np.abs(res[keep])); keep = np.abs(res) < 3 * s
            if best is None or s < best[0]:
                best = (s, alpha, coef, S)
        s, alpha, coef, S = best
        model = coef[0] * S + coef[1] + coef[2] * X + coef[3] * Yv + coef[4] * (X * X + Yv * Yv)
        level = float(np.median(model[r < 0.5]))
        th = np.arctan2(yy - mcy, xx - mcx)
        nb = 360
        ti = ((th + np.pi) / (2 * np.pi) * nb) % nb
        rb = np.clip((r * mR).astype(np.int32), 0, int(mR) + 2)
        s_lo = max(3.0, 0.012 * mR); s_hi = 0.3 * mR

        def _structure(lim, along_limb, sd_ref=None):
            """Steps 4-5b inside `lim` R. `along_limb`: the edge pass -- also take
            off what runs along the limb in the edge zone (step 4c) and keep the
            fine scale only (build 8); its spread is held to `sd_ref`, the inner
            pass's. Returns the structure and the inner disc's ring spread."""
            inside = r < lim
            E = np.where(inside, M - model, 0.0)
            # ---- 4: what follows the limb: the median of every ring
            msk = inside.ravel()
            med = np.zeros(rb.max() + 1)
            order = np.argsort(rb.ravel()[msk], kind="stable")
            rv = rb.ravel()[msk][order]; ev = E.ravel()[msk][order]
            cuts = np.searchsorted(rv, np.arange(rb.max() + 2))
            for k in range(rb.max() + 1):
                seg = ev[cuts[k]:cuts[k + 1]]
                if seg.size >= 8:
                    med[k] = np.median(seg)
            med = ndimage.gaussian_filter1d(med, 2.0)
            E = np.where(inside, E - med[rb], 0.0)
            # ---- 4b: glare that hugs the limb in some directions only (a bright
            # streamer's spill): per direction, the amplitude of a profile rising
            # toward the limb, b(r) = exp((r - 0.97) / 0.06), fitted on 0.6-0.95 R,
            # smoothed over 10 degrees and taken off
            bprof = np.exp((r - 0.97) / 0.06)
            zone = inside & (r > 0.6) & (r < 0.95)
            tb = ((th[zone] + np.pi) / (2 * np.pi) * nb).astype(np.int32) % nb
            ez, bz = E[zone], bprof[zone]
            num = np.bincount(tb, ez * bz, nb); den = np.bincount(tb, bz * bz, nb)
            amp = np.where(den > 0, num / np.maximum(den, 1e-12), 0.0)
            amp = ndimage.gaussian_filter1d(amp, 10.0, mode="wrap")
            a_map = ndimage.map_coordinates(amp, [ti.ravel()], order=1, mode="grid-wrap").reshape(ti.shape)
            E = np.where(inside, E - a_map * bprof, 0.0)
            if along_limb:
                # ---- 4c: WHAT RUNS ALONG THE LIMB (build 7). From 0.82 R outward,
                # at every radius, the running median over 20 degrees of arc is
                # taken off. Leftover glare follows the limb for long stretches;
                # lunar features are compact and stay.
                nth = 1440
                rr_ = np.arange(0.82 * mR - 4, lim * mR + 1, 1.0)
                tt = np.linspace(-np.pi, np.pi, nth, endpoint=False)
                Pp = ndimage.map_coordinates(E, [mcy + rr_[:, None] * np.sin(tt)[None],
                                                 mcx + rr_[:, None] * np.cos(tt)[None]],
                                             order=1, mode="nearest")
                kk = max(3, int(round(20.0 / 360.0 * nth)) | 1)
                Bp = ndimage.median_filter(Pp, size=(1, kk), mode="wrap")
                Bp = ndimage.uniform_filter1d(Bp, 9, axis=1, mode="wrap")
                Bp = ndimage.gaussian_filter1d(Bp, 2.0, axis=0)
                Bc = ndimage.map_coordinates(Bp, [np.clip(r * mR - rr_[0], 0, len(rr_) - 1),
                                                  ((th + np.pi) / (2 * np.pi) * nth) % nth],
                                             order=1, mode="nearest")
                wz = np.clip((r - 0.82) / 0.04, 0, 1); wz = wz * wz * (3 - 2 * wz)
                E = np.where(inside, E - wz * Bc, 0.0)
            # ---- 5: band-pass by normalised convolution inside the disc: above the
            # sensor's fine scale, below 0.3 R
            w = inside.astype(np.float64)

            def nconv(a, sg):
                return _pgauss(a * w, sg) / np.maximum(_pgauss(w, sg), 1e-6)
            if along_limb:
                # THE EDGE PASS KEEPS THE FINE SCALE ONLY (build 8, 2026-11-16).
                # After 4c nothing large is left from 0.86 R outward (sector
                # means within 0.2 of the inner disc's spread on a weak test
                # set), but further in the large patches are still there. The
                # wide low-pass averaged them into the rim and its subtraction
                # put the OPPOSITE sign there: under bright patches of +2 the rim
                # sat at -1.0 to -1.4, under dark ones of -0.7 at +0.8 to +1.0.
                # So the wide low-pass is not subtracted here.
                Eb = np.where(inside, nconv(E, s_lo), 0.0)
            else:
                Eb = np.where(inside, nconv(E, s_lo) - nconv(E, s_hi), 0.0)
            # ---- 5b: THE LIMB BAND (build 2, 2026-10-07). Measured on two
            # test sets (one of them the reference set): after the band-pass the structure
            # grew toward the edge (ring std 1.8x the centre's at 0.75-0.86 R) and the
            # rings there sat dark (layer mean 0.44 against 0.50), then both dropped to
            # flat within 0.05 R -- a visible rim. Real maria lose contrast toward the
            # limb (foreshortening), so that growth is residual glare and the
            # normalised convolution's edge bias. Per ring: the median taken off, and
            # the spread held to the inner disc's.
            ring_med = np.zeros(rb.max() + 1); ring_sd = np.zeros(rb.max() + 1)
            ebv = Eb.ravel()[msk][order]
            for k in range(rb.max() + 1):
                seg = ebv[cuts[k]:cuts[k + 1]]
                if seg.size >= 8:
                    ring_med[k] = np.median(seg)
                    ring_sd[k] = 1.4826 * np.median(np.abs(seg - ring_med[k]))
            ring_med = ndimage.gaussian_filter1d(ring_med, 3.0)
            ring_sd = ndimage.gaussian_filter1d(ring_sd, 3.0)
            _rk = np.arange(rb.max() + 1) / mR
            _ref = ring_sd[(_rk > 0.2) & (_rk < 0.6)]
            sd0 = float(np.median(_ref)) if _ref.size else 0.0
            if sd_ref is not None:
                sd0 = float(sd_ref)       # the edge pass: held to the inner pass's spread
            gain = (np.where(ring_sd > 1e-12, np.minimum(1.0, sd0 / np.maximum(ring_sd, 1e-12)), 1.0)
                    if sd0 > 0 else np.ones_like(ring_sd))
            return np.where(inside, (Eb - ring_med[rb]) * gain[rb], 0.0), sd0

        # THE INNER DISC AS IT WAS, THE EDGE FROM A SECOND PASS (build 7,
        # 2026-11-16). Measured on the reference set, ring by ring against the
        # real Moon (band 0.02-0.12 R): the structure matched 0.5-0.6 out to
        # 0.90 R but was faded from 0.80 R, and beyond 0.90 R the match fell to
        # 0.1-0.2 while the two halves of the frames agreed MORE -- leftover
        # glare, steady, not lunar. Build 6 took that glare off with a narrow
        # limb profile per direction; it kept the detail (0.56 / 0.48 at
        # 0.90-0.95 / 0.95-0.985 R) but carved a dark band along the rim in two
        # quadrants and changed the inner disc (0.87 against build 5). Here the
        # inner disc is build 5's, untouched; from 0.86 to 0.92 R it hands over
        # to a second pass that reaches 0.995 R and has what runs along the limb
        # taken off (4c): 0.63 / 0.57 / 0.27 in the same rings, and a quarter of
        # build 6's limb-parallel banding.
        # BUILD 8 (2026-11-16): the edge pass carries the fine scale only (see
        # step 5) and the hand-over runs from 0.80 to 0.92 R, so the large
        # patches of the inner disc run out softly instead of ending in a band
        # of the opposite sign. Weak test set: rim sector means within 0.2 of
        # the inner spread (build 7: -1.4 to +1.0). Reference set: whole-disc
        # match with the real Moon 0.646 (build 7: 0.638).
        Ea, _sd_in = _structure(0.97, False)
        Ee, _ = _structure(0.995, True, sd_ref=_sd_in)
        wb = np.clip((r - 0.80) / 0.12, 0, 1); wb = wb * wb * (3 - 2 * wb)
        # BLOBS AT THE EDGE FADED OUT (build 9, 2026-11-17). On a 2023 test set
        # the edge pass carried compact patches of up to -11 times the inner
        # disc's spread at 0.90-0.97 R (and a few bright ones), where the stack
        # itself had nothing of the kind: leftover of the glare model beside
        # bright limb features. Shown, they were black holes along the limb.
        # Maria are not that strong near the limb (foreshortening). Where the
        # edge pass, smoothed over 0.012 R, stands more than 1.5 times the inner
        # spread off zero, it fades to the plain disc (gone at 3 times), the
        # fade widened over the patch's shoulders and softened over 0.02 R, so
        # what remains is the smooth disc, not a hole with a rim. Measured:
        # that set 22 % of the 0.80-0.985 R zone faded, a high-Sun reference set
        # 3 % (no visible change), an extreme close-up set 18 % (its limb
        # blobs gone). The inner disc (below 0.80 R) is not touched.
        if _sd_in > 0:
            _a = _pgauss(np.abs(Ee) / _sd_in, 0.012 * mR)
            _c = np.clip((_a - 1.5) / 1.5, 0, 1); _c = 1.0 - _c * _c * (3 - 2 * _c)
            _c = ndimage.minimum_filter(_c, size=max(3, int(round(0.03 * mR)) | 1))
            _c = np.clip(_pgauss(_c, 0.02 * mR), 0, 1)
            Ee = Ee * _c
        Eb = Ea * (1.0 - wb) + Ee * wb
        # ---- 6: fade out at the very edge, 0.95-0.985 R
        fade = np.clip((0.985 - r) / 0.035, 0, 1); fade = fade * fade * (3 - 2 * fade)
        Eb *= fade
        return Eb, r, level, alpha, s_lo, s_hi

    if os.environ.get("ECLIPSEFORGE_MOON_DUMP", "") not in ("", "0"):
        # lab diagnostics: the stacked disc (all frames and the two halves) and
        # the merge's green light on the same grid, as _layer receives them
        try:
            _n = M.shape[0]
            _yy, _xx = np.mgrid[:_n, :_n].astype(np.float64)
            _my = cy + (_yy - mcy) * 2.0; _mx = cx + (_xx - mcx) * 2.0
            _a0, _a1 = max(int(_my.min()) - 2, 0), min(int(_my.max()) + 3, H)
            _b0, _b1 = max(int(_mx.min()) - 2, 0), min(int(_mx.max()) + 3, W)
            _I = ndimage.map_coordinates(np.asarray(_hr[_a0:_a1, _b0:_b1, 1], np.float32),
                                         [_my - _a0, _mx - _b0], order=1, mode="constant", cval=0.0)
            np.savez(os.path.join(wd, MOON_DUMP), M=M.astype(np.float32),
                     A=(hacc[0] / np.maximum(hw[0], 1e-9)).astype(np.float32) if hacc[0] is not None else M * 0,
                     B=(hacc[1] / np.maximum(hw[1], 1e-9)).astype(np.float32) if hacc[1] is not None else M * 0,
                     I=_I.astype(np.float32), m=np.array([mcy, mcx, mR]), g=np.array([cy, cx, R, H, W]),
                     n=np.array([nused] + list(_hn)), t=np.array(sorted(t for _, t in use)))
            del _yy, _xx, _my, _mx, _I
            _log(progress, "Moon layer: stack written to " + MOON_DUMP + " (ECLIPSEFORGE_MOON_DUMP)", None)
        except Exception as _e:
            _log(progress, f"[warn] Moon layer: dump failed ({_e})", None)
    Eb, r, level, alpha, s_lo, s_hi = _layer(M)
    # MARIA AGREEMENT (2026-10-28). The same layer from the first and the
    # second half of the frames (in shooting order), each stacked on its own.
    # Maria are fixed to the Moon and are in both; noise is not, and leftover
    # glare is fixed to the Sun, which the Moon moves across during totality, so
    # it lines up less well between an early and a late half. The correlation c
    # of the two halves inside 0.8 R is the share of signal in either half;
    # the whole stack has half the noise variance, so its signal share is
    # 2c / (1 + c). The layer is scaled by that share: earthshine as strong as
    # the data supports, and a smooth disc where the maria are not there.
    _agree = None; _share = 1.0; _EA = _EB = None
    if min(_hn) >= 1 and all(a is not None for a in hacc):
        try:
            _EA = _layer(hacc[0] / np.maximum(hw[0], 1e-9))[0]
            _EB = _layer(hacc[1] / np.maximum(hw[1], 1e-9))[0]
            _m = r < 0.8
            _a = _EA[_m] - _EA[_m].mean(); _b = _EB[_m] - _EB[_m].mean()
            _agree = float(np.dot(_a, _b) / max(np.sqrt(np.dot(_a, _a) * np.dot(_b, _b)), 1e-30))
            _c = max(_agree, 0.0)
            _share = 2.0 * _c / (1.0 + _c)
            if _c < MARIA_MIN_AGREE:
                _share = 0.0
            del _a, _b, _m
        except Exception as _e:
            _log(progress, f"[warn] Moon layer: maria agreement not measured ({_e})", None)
            _agree = None; _share = 1.0; _EA = _EB = None
    # THE REAL MOON (2026-11-14, lunarcheck.py). The halves test above is
    # internal and can be fooled both ways: it kept 22 % of maria that are
    # there on one test set, and 89 % of a pattern that is not the Moon on a
    # set of four frames. So the
    # structure is also compared with a map of the real Moon's surface, in
    # every orientation, against random patterns of the same texture. The map
    # decides only whether, and how strongly, the user's own structure goes
    # into the result; nothing from the map is drawn.
    _share_halves = _share
    _chk = None
    try:
        from . import lunarcheck as _lc
        _chk = _lc.check(_lc.disc_from(Eb, mcy, mcx, mR))
        if _chk["verdict"] == "confirmed":
            _share = max(_share, _lc.share_for(_chk["best"]))
        elif _chk["verdict"] == "rejected":
            _share = 0.0
    except Exception as _e:
        _log(progress, f"[warn] Moon layer: the check against the lunar map did not run ({_e}); "
                       "the halves test alone decides", None)
        _chk = None
    # FORCE EARTHSHINE (2026-11-14). People like earthshine, and a strict
    # check takes it away from some sets where it may be there. So a second
    # variant is built for the Moon group's "Force earthshine" box. It still
    # has to beat chance against the real Moon; it only accepts what the
    # default does not:
    #   - all frames match above chance but not by CONFIRM_SD: the structure
    #     at share_for(best) (at least 50 %) instead of the halves test's share;
    #   - all frames do not beat chance, but one half of the frames does: that
    #     half, with about 1.4x the noise. The better of two halves gives
    #     chance two tries, and the random patterns get two tries as well
    #     (lunarcheck.check_halves).
    # Nothing that is at chance level against the real Moon is forced in.
    _Ef = None; _share_f = None
    if _chk is None:
        _force = {"available": False, "why": "the check against the lunar map did not run"}
    elif _chk["verdict"] == "confirmed":
        _force = {"available": False, "why": "not needed: earthshine confirmed"}
    elif _chk["verdict"] == "undecided":
        _sf = float(_lc.share_for(_chk["best"]))
        if _sf > _share + 0.02:
            _Ef, _share_f = Eb, _sf
            _force = {"available": True, "source": "all frames", "share": round(_sf, 3), "check": _chk}
        else:
            _force = {"available": False, "why": "not needed: it would not add earthshine"}
    elif _EA is None or _EB is None:
        _force = {"available": False, "why": "nothing to use: all frames do not match the real "
                                             "Moon better than chance, and there are no two halves to try"}
    else:
        try:
            _ch = _lc.check_halves(_lc.disc_from(_EA, mcy, mcx, mR), _lc.disc_from(_EB, mcy, mcx, mR))
            if _ch["verdict"] != "rejected":
                _Ef = _EA if _ch["half"] == 0 else _EB
                _share_f = float(_lc.share_for(_ch["best"]))
                _force = {"available": True, "share": round(_share_f, 3), "check": _ch,
                          "source": "the " + ("first" if _ch["half"] == 0 else "second") + " half of the frames"}
            else:
                _force = {"available": False, "check": _ch,
                          "why": "nothing to use: neither all frames nor either half of them "
                                 "matches the real Moon better than chance"}
        except Exception as _e:
            _force = {"available": False, "why": f"the halves could not be checked ({_e})"}
    del _EA, _EB
    # ---- onto the merge grid
    Yg = np.arange(H, dtype=np.float64)

    def _to_grid(Ex, sh):
        sdx = float(1.4826 * np.median(np.abs(Ex[r < 0.8] - np.median(Ex[r < 0.8])))) or 1.0
        e = 0.5 + np.float64(0.5 * sh) * np.tanh(Ex / (2.5 * sdx))
        return _onto_merge(e), sdx

    def _onto_merge(e):
        out = np.full((H, W), 0.5, np.float16)
        y0m, y1m = max(int(cy - 1.05 * R), 0), min(int(cy + 1.05 * R) + 1, H)
        x0m, x1m = max(int(cx - 1.05 * R), 0), min(int(cx + 1.05 * R) + 1, W)
        gyy, gxx = np.mgrid[y0m:y1m, x0m:x1m].astype(np.float64)
        _dy, _dx = gyy - cy, gxx - cx
        _uq = geo.get("unsquash")
        if _uq and _uq.get("S"):
            # THE MERGE IS STRETCHED ROUND (pipeline.fit_limb_oval): the raw
            # frames this layer is stacked from are not, so the merge's offsets
            # from the disc centre are taken back through S^-1 before they index
            # the stack -- the Moon's oval then lands exactly on the merge's
            # round disc instead of up to 1.5 % of R off along one axis.
            _Si = np.linalg.inv(np.asarray(_uq["S"], float))
            _dy, _dx = (_Si[0, 0] * _dy + _Si[0, 1] * _dx), (_Si[1, 0] * _dy + _Si[1, 1] * _dx)
        py = mcy + _dy / 2.0; px = mcx + _dx / 2.0
        out[y0m:y1m, x0m:x1m] = ndimage.map_coordinates(e, [py, px], order=1, mode="constant", cval=0.5).astype(np.float16)
        return out

    out, sd = _to_grid(Eb, _share)
    np.save(os.path.join(wd, MOON_FILE), out)
    del out
    _fp = os.path.join(wd, MOON_FORCE_FILE)
    if _Ef is not None:
        np.save(_fp, _to_grid(_Ef, _share_f)[0])
    elif os.path.exists(_fp):
        os.remove(_fp)
    del _Ef
    meta = {"build": MOON_BUILD, "frames": nused, "exposures": sorted(set(round(t, 6) for _, t in use)),
            "psf_alpha": alpha, "disc_level": level, "structure_rms": sd,
            "structure_rel": sd / max(level, 1e-9), "band_px_halfres": [s_lo, s_hi],
            "half_frames": list(_hn), "maria_agreement": _agree, "maria_share": _share,
            "maria_share_halves": _share_halves, "maria_check": _chk, "force": _force}
    json.dump(meta, open(os.path.join(wd, MOON_META), "w"), indent=1)
    _log(progress, f"[ok] Moon layer: {nused} frames, scatter model r^-{alpha:g}, "
                   f"maria {100 * sd / max(level, 1e-9):.2f}% of disc light", None)
    if _agree is None:
        _log(progress, "Moon layer: maria agreement not measured (needs frames in both "
                       "halves of the sequence)" + ("" if _chk else "; earthshine layer used as built"), None)
    elif _chk is not None:
        _log(progress, f"Moon layer: the first and second half of the frames agree {_agree:.2f} on the disc "
                       f"(on its own that would give {100 * _share_halves:.0f} %)", None)
    elif _share <= 0.0:
        _log(progress, f"[warn] Moon layer: the two halves of the frames agree only "
                       f"{_agree:.2f} on the disc -- no maria found in the data, only noise "
                       f"or glare; the Moon is left smooth", None)
    else:
        _log(progress, f"Moon layer: maria agreement {_agree:.2f} between the first and "
                       f"second half of the frames -- earthshine structure at "
                       f"{100 * _share_halves:.0f} %", None)
    if _chk is not None:
        _txt = (f"best match {_chk['best']:.2f} against {_chk['random_mean']:.2f} +- {_chk['random_sd']:.2f} "
                f"for random patterns of the same texture ({_chk['sd']:+.1f} sd)")
        if _chk["verdict"] == "confirmed":
            _log(progress, f"[ok] Moon layer: the structure on the disc matches the real Moon's maria "
                           f"({_txt}): earthshine at {100 * _share:.0f} %. The lunar map is used for this "
                           "check only; the picture shows your own data", None)
        elif _chk["verdict"] == "rejected":
            _log(progress, f"[warn] Moon layer: the structure on the disc does not match the real Moon's "
                           f"maria ({_txt}): glare or noise, the Moon is left smooth", None)
        else:
            _log(progress, f"Moon layer: the match with the real Moon's maria is not clear ({_txt}); "
                           f"earthshine at {100 * _share:.0f} % from the halves test", None)
    if _force.get("available"):
        _c = _force["check"]
        _log(progress, f"Moon layer: \"Force earthshine\" (Moon group) can use the earthshine from "
                       f"{_force['source']} (best match {_c['best']:.2f} against {_c['random_mean']:.2f} "
                       f"+- {_c['random_sd']:.2f} for random patterns, {_c['sd']:+.1f} sd): above chance but "
                       f"not confirmed, so it may be luck, and "
                       + ("half the frames carry more noise" if "half" in _force["source"]
                          else "it goes in at " + f"{100 * _force['share']:.0f} %"), None)
    elif _chk is not None and _chk["verdict"] == "rejected":
        _log(progress, f"Moon layer: \"Force earthshine\" not available here: {_force['why']}", None)
    return meta
