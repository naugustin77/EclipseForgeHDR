"""The simple stack: align, average, merge on measured ratios. The fallback.

WHY THIS EXISTS. The normal path fits a photometric ladder, a shared pedestal,
per-channel floors, an azimuthal affine per tier and a feather. Each is there
for a measured reason, and when the result is wrong none of them can be blamed
individually because there is no baseline to compare against. This is the
baseline, and it is offered to the user as a fallback: when nothing else gives a
clean stack, this gives a stack.

    read -> per-frame per-channel pedestal -> align -> mean per tier
         -> exposure ratios MEASURED between tiers -> hat-weighted linear merge

Then it hands the merge to `importhdr.build_from_rgb`, which is what an imported
HDR already goes through, so the result opens in the preview with MGN, FNRGF,
NAFE-VN, the inner-corona layer, Pellett and the prominence gate exactly like
any other run.

WHAT IT DELIBERATELY DOES NOT DO, and the report says so: no dark, no flat, no
hot-pixel repair, no lunar masking per tier, no ladder solve, no LDIC, no
feather, no tier projection. Anything that can fail on an awkward bracket is
absent, which is the point.

THE ONE THING IT DOES THAT THE NORMAL PATH DOES NOT. The per-frame pedestal is
the median of the frame's own corners, per channel -- which is black level PLUS
sky. The normal path removes one constant shared by every tier, correct for a
black level and wrong for sky, because sky scales with exposure exactly like the
corona and so survives the merge as a constant added to every pixel. On Val
Italo's set that constant is 2.4x the corona in red and 6.8x in blue at 2.85 R,
with the opposite colour. Removing it is most of why this path looks clean.

ITS LIMIT, stated because a fallback that hides its own failure is worse than no
fallback. One number per frame per channel cannot follow a sky that varies
ACROSS the frame. On a small field it works. On the 600 mm reference set the sky varies
about 28% corner to corner and a constant leaves broad colour blobs beyond ~4 R
-- the failure documented in docs/SKY_SUBTRACTION.md, which no version of this
has escaped. Judge it by eye.
"""
from __future__ import annotations
import os

import numpy as np


# ---------------------------------------------------------------- reading ---

def _planes_direct(path):
    """A 3-plane FITS read as planes, or None if this file is not one.

    `FitsFrame` re-mosaics an RGB cube into a Bayer CFA so the rest of the app
    can treat every input alike, and the pipeline then interpolates it back.
    That round trip throws three quarters of each channel away and guesses it
    again, and the guess is gradient-directed, so it lands hardest at the lunar
    limb. Siril and PixInsight write RGB cubes. For those the planes are right
    there, and reading them is both more accurate and cheaper.
    """
    if os.path.splitext(path)[1].lower() not in (".fit", ".fits", ".fts"):
        return None
    from .fits import read_fits
    d, h = read_fits(path)
    if d.ndim != 3:
        return None
    if d.shape[0] in (3, 4):
        a = np.asarray(d[:3], np.float32)
    elif d.shape[-1] in (3, 4):
        a = np.ascontiguousarray(np.moveaxis(d[..., :3], -1, 0), np.float32)
    else:
        return None
    ped = h.get("PEDESTAL", h.get("BLKLEVEL", h.get("BLACKLEV")))
    try:
        a = a - float(ped)
    except (TypeError, ValueError):
        pass
    sat = h.get("SATURATE", h.get("DATAMAX", h.get("SATLEVEL")))
    try:
        sat = float(sat)
    except (TypeError, ValueError):
        sat = 65535.0
    return a, sat


def _load(path, demosaic_method="mhc"):
    """(3, H, W) float32 raw ADU, pedestal not yet removed, and sat_level."""
    direct = _planes_direct(path)
    if direct is not None:
        return direct
    from .raw import open_frame, demosaic
    rf = open_frame(path)
    rgb = demosaic(rf.bayer, demosaic_method)
    return (np.ascontiguousarray(rgb.transpose(2, 0, 1), np.float32),
            float(rf.sat_level))


def _pedestal(a, frac=0.04):
    """Per-channel median of the four corners: black level AND sky, together.

    Taken per channel because a colour cast that changes with radius is exactly
    a per-channel additive term, and one shared number cannot express it.
    """
    h, w = a.shape[1:]
    cy, cx = max(int(h * frac), 8), max(int(w * frac), 8)
    out = np.empty(3, np.float64)
    for c in range(3):
        p = a[c]
        out[c] = np.median(np.concatenate([
            p[:cy, :cx].ravel(), p[:cy, -cx:].ravel(),
            p[-cy:, :cx].ravel(), p[-cy:, -cx:].ravel()]))
    return out


# -------------------------------------------------------------- alignment ---

def _hipass(x, s=24.0):
    k = int(max(1, round(s)))
    c = np.cumsum(np.pad(x, ((k, k), (k, k)), mode="edge"), axis=0)
    c = (c[2 * k:, :] - c[:-2 * k, :]) / (2 * k)
    c = np.cumsum(c, axis=1)
    c = (c[:, 2 * k:] - c[:, :-2 * k]) / (2 * k)
    return x - c[:x.shape[0], :x.shape[1]]


def _feature(g):
    """Half-resolution high-passed log green.

    HALF RESOLUTION ON PURPOSE: a full-res FFT pair on a 44 Mpx frame costs
    seconds each and there are tens of frames. The parabolic sub-pixel fit
    lands within about 0.1 of a half-res pixel, so 0.2 px at full scale --
    below what this path's own bilinear resampling costs anyway.
    """
    h, w = (g.shape[0] // 2) * 2, (g.shape[1] // 2) * 2
    g2 = g[:h, :w].reshape(h // 2, 2, w // 2, 2).mean((1, 3))
    return _hipass(np.log1p(np.clip(g2, 0, None) / 30.0), 24.0)


def _pcorr(a, b):
    """Sub-pixel shift taking b onto a. Plain phase correlation.

    NOT the pipeline's `_semi_phase_shift`, deliberately: this path is a control
    for the normal one and should share as little with it as possible.
    """
    A = np.fft.rfft2(a)
    B = np.fft.rfft2(b)
    R = A * np.conj(B)
    m = np.abs(R)
    R = R / np.maximum(m, 1e-12 * m.max())
    c = np.fft.irfft2(R, s=a.shape)
    i = np.unravel_index(np.argmax(c), c.shape)
    sub = []
    for ax, ix in enumerate(i):
        n = c.shape[ax]
        im, ip = (ix - 1) % n, (ix + 1) % n
        y0 = c[im, i[1]] if ax == 0 else c[i[0], im]
        y2 = c[ip, i[1]] if ax == 0 else c[i[0], ip]
        d = y0 - 2 * c[i] + y2
        off = 0.0 if abs(d) < 1e-12 else 0.5 * (y0 - y2) / d
        v = ix + off
        if v > n / 2:
            v -= n
        sub.append(v)
    return np.array(sub)


def _shift(x, dy, dx, out=None):
    """Bilinear shift. Outside the source is NaN so it can be counted out."""
    H, W = x.shape
    yy = np.arange(H, dtype=np.float32) - dy
    xx = np.arange(W, dtype=np.float32) - dx
    y0 = np.floor(yy).astype(np.int32)
    x0 = np.floor(xx).astype(np.int32)
    fy = (yy - y0)[:, None].astype(np.float32)
    fx = (xx - x0)[None, :].astype(np.float32)
    y0c = np.clip(y0, 0, H - 2)
    x0c = np.clip(x0, 0, W - 2)
    r = (x[y0c][:, x0c] * (1 - fy) * (1 - fx)
         + x[y0c][:, x0c + 1] * (1 - fy) * fx
         + x[y0c + 1][:, x0c] * fy * (1 - fx)
         + x[y0c + 1][:, x0c + 1] * fy * fx)
    r[(yy < 0) | (yy > H - 1), :] = np.nan
    r[:, (xx < 0) | (xx > W - 1)] = np.nan
    return r


# ------------------------------------------------------------------ merge ---

def _ratio(a, b, sat, lo=200.0):
    """Median b/a over pixels well exposed in BOTH: the measured ratio.

    Measured on green and applied to all three channels. On a third tester's 23
    tiers the channels' own ratios agree to under 1.2% per step and diverge
    6.8% over the whole ladder, so per-channel freedom buys nothing and costs a
    noisier estimate.
    """
    m = (np.isfinite(a) & np.isfinite(b) & (a > lo) & (b > lo)
         & (a < sat * 0.8) & (b < sat * 0.8))
    if m.sum() < 2000:
        return np.nan, int(m.sum())
    return float(np.median(b[m] / a[m])), int(m.sum())


def _weight(v, sat):
    """Hat: zero at the noise floor, zero towards saturation.

    Weighting by exposure time instead -- what the normal merge does -- gives a
    short tier real weight in the outer field where it holds nothing but noise
    and whatever the pedestal step missed. This does not.
    """
    return (np.clip((v - 60.0) / 400.0, 0, 1)
            * np.clip((sat * 0.97 - v) / (sat * 0.12), 0, 1)).astype(np.float32)


# -------------------------------------------------------------------- run ---

def run(folder, progress, denoise="fine", demosaic_method="mhc",
        fnrgf_preset="ours"):
    """Stack `folder` the simple way and build every layer the renderer needs."""
    from .raw import list_raws, read_exif
    from . import importhdr
    from . import __version__

    try:
        progress.bar_detail = 0.45
    except Exception:
        pass

    files = list_raws(folder)
    if not files:
        raise RuntimeError("no readable frames in this folder")

    # The HEADER exposure is used ONLY as a grouping key. Every photometric
    # number below comes from the pixels -- which matters, because on Val
    # Italo's set the headers are Siril "manual group assignment" labels
    # stepping a uniform 1.5x where the pixels measure 1.19 / 1.49 / 1.99.
    tiers = {}
    skipped = 0
    for p in files:
        try:
            sec = float(read_exif(p)[0])
        except Exception:
            skipped += 1
            continue
        tiers.setdefault(sec, []).append(p)
    secs = sorted(tiers)
    if len(secs) < 2:
        raise RuntimeError(
            "the simple stack needs at least two exposure groups; this folder "
            "has %d" % len(secs))
    n_used = sum(len(v) for v in tiers.values())
    progress.log(f"simple stack: {n_used} frames in {len(secs)} exposure groups"
                 + (f" ({skipped} unreadable, skipped)" if skipped else ""), 0.02)
    progress.log("no dark, no flat, no hot-pixel repair, no ladder solve, no "
                 "LDIC, no feather — this path is the fallback and does the "
                 "least a program can do and still merge a bracket", None)

    ref_feat = None
    num = den = None
    prev_g = None
    sat = None
    rel = 1.0
    short_acc, short_n = None, 0
    ladder = []
    for i, s in enumerate(secs):
        fs = sorted(tiers[s])
        acc = cnt = None
        shifts = []
        for p in fs:
            f, sl = _load(p, demosaic_method)
            if sat is None:
                sat = sl
            f -= _pedestal(f)[:, None, None]
            ft = _feature(f[1])
            if ref_feat is None:
                ref_feat = ft
                dy = dx = 0.0
            else:
                if ft.shape != ref_feat.shape:
                    raise RuntimeError(
                        "the frames in this folder are not all the same size")
                dy, dx = 2.0 * _pcorr(ref_feat, ft)
                if not np.isfinite([dy, dx]).all() or max(abs(dy), abs(dx)) > 400:
                    dy = dx = 0.0
            shifts.append((float(dy), float(dx)))
            g = np.stack([_shift(f[c], dy, dx) for c in range(3)])
            del f
            ok = np.isfinite(g[0]).astype(np.float32)
            np.nan_to_num(g, copy=False)
            acc = g if acc is None else acc + g
            cnt = ok if cnt is None else cnt + ok
            del g, ok
        st = (acc / np.maximum(cnt, 1.0)).astype(np.float32)
        st[:, cnt < 0.5] = np.nan
        del acc, cnt

        if prev_g is not None:
            r, npx = _ratio(prev_g, st[1], sat)
            if not np.isfinite(r):
                r = s / secs[i - 1]
                progress.log(f"  {secs[i-1]:.5g}s -> {s:.5g}s: too few shared "
                             f"pixels to measure the ratio; falling back to the "
                             f"headers' {r:.4f}", None)
            rel *= r
            ladder.append({"from": secs[i - 1], "to": s,
                           "header": s / secs[i - 1], "measured": r, "px": npx})
        prev_g = st[1].copy()

        if num is None:
            num = np.zeros_like(st)
            den = np.zeros_like(st)
        for c in range(3):
            v = st[c]
            okm = np.isfinite(v)
            vv = np.where(okm, v, 0.0).astype(np.float32)
            w = np.where(okm, _weight(vv, sat), np.float32(0.0))
            num[c] += w * vv / np.float32(rel)
            den[c] += w

        # the four shortest tiers, in scene units: the inner-corona layer's
        # independent source, which an imported HDR does not have
        if i < 4:
            sc = (np.nan_to_num(st) / np.float32(rel))
            lum_s = (0.2126 * sc[0] + 0.7152 * sc[1] + 0.0722 * sc[2])
            short_acc = lum_s if short_acc is None else short_acc + lum_s
            short_n += 1
            del sc, lum_s
        del st
        progress.log("  %10.5g s  %2d frames  shift rms %.2f px  max %.2f"
                     % (s, len(fs),
                        float(np.sqrt((np.array(shifts) ** 2).sum(1).mean())),
                        float(np.abs(np.array(shifts)).max())),
                     0.05 + 0.35 * (i + 1) / len(secs))

    hdr = np.where(den > 0, num / np.maximum(den, 1e-9), 0.0)
    del num, den
    rgb = np.ascontiguousarray(np.moveaxis(hdr, 0, -1), np.float32)
    del hdr
    short = (short_acc / max(short_n, 1)).astype(np.float32) if short_acc is not None else None

    progress.log("exposure ladder, measured from the pixels: header says "
                 "1 .. %.4g, the pixels say 1 .. %.4g"
                 % (secs[-1] / secs[0], rel), None)
    for L in ladder:
        progress.log("  %10.5g -> %-10.5g header %8.4f   measured %8.4f"
                     % (L["from"], L["to"], L["header"], L["measured"]), None)

    # scale to the same comfortable working range an import gets; every layer
    # downstream is scale-free
    hi = float(np.percentile(rgb, 99.9))
    if hi > 0:
        rgb *= np.float32(50000.0 / hi)
        if short is not None:
            short *= np.float32(50000.0 / hi)

    return importhdr.build_from_rgb(
        folder, rgb, progress, denoise=denoise, short_lum=short,
        fnrgf_preset=fnrgf_preset,
        stats={"n_files": n_used, "mode": "simple stack",
               "simple_ladder": ladder,
               "simple_ladder_span": {"header": secs[-1] / secs[0],
                                      "measured": rel},
               "tiers": [{"sec": s, "n": len(tiers[s])} for s in secs]},
        opts={"mode": "simple", "denoise": denoise,
              "fnrgf_preset": fnrgf_preset,
              "demosaic": demosaic_method, "n_files": n_used,
              "secs": list(secs)})
