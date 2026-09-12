#!/usr/bin/env python3
"""The dumb path: align, average, merge on MEASURED exposure ratios. Nothing else.

WHY THIS EXISTS. The shipped pipeline fits a photometric ladder, a shared
pedestal, per-channel floors, an azimuthal affine per tier, a feather and a
detail chain. Each of those is there for a measured reason. But when the result
is wrong there is no way to tell WHICH of them is wrong, because there is no
baseline to compare against. This is the baseline: the least a program can do
and still produce an HDR corona.

    read -> per-frame per-channel pedestal -> align -> mean per tier
         -> exposure ratios MEASURED between tiers -> hat-weighted linear merge

That is all. No ladder solve, no LDIC, no FNRGF, no tier projection, no feather.

WHAT THE PEDESTAL STEP ACTUALLY REMOVES, and why it is the interesting part.
The median of the frame's own corners is black level PLUS sky. Subtracting it
per frame removes both, and the sky is the term the shipped chain has no model
for: `_fit_pedestal` removes ONE constant shared by every tier, which is right
for a black level and wrong for sky, because sky scales with the exposure
exactly like the corona and so survives the merge as a constant added to every
pixel. On a third tester's set that constant is 2.4x the corona in red and 6.8x in
blue at 2.85 R, with the opposite colour -- which is the warm-inside,
blue-outside picture he reported.

THE LIMIT, stated plainly. One number per frame per channel cannot follow a sky
that varies across the frame. On Val's set (2.9 R to the frame edge) it works.
On the 600 mm reference set the sky's brightness varies about 28% corner to corner, and
subtracting a constant leaves broad colour blobs beyond ~4 R -- the failure
documented in docs/SKY_SUBTRACTION.md, which no version of this has escaped.
So: judge the output by eye, expect it to win on some sets and not others, and
do not read the far field of a set where it loses.

    python tools/simple_stack.py FOLDER [--scale 2] [--out DIR] [--max-frames N]
"""
import argparse
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "src"))


# ---------------------------------------------------------------- reading ---

def _planes_direct(path):
    """A 3-plane FITS read as planes, or None if this file is not one.

    WHY THIS EXISTS. `FitsFrame` re-mosaics an RGB cube into a Bayer CFA so the
    rest of the app can treat every input the same way, and `demosaic` then
    interpolates it back. That round trip is lossless for nothing: it throws
    three quarters of each channel away and guesses it again, and the guess is
    gradient-directed, so it lands hardest exactly at the lunar limb. Siril and
    PixInsight write RGB cubes, which is what most FITS anyone imports will be,
    and for those the planes are right there. Reading them directly is more
    accurate AND avoids pulling scipy in for a demosaic that should not happen.
    """
    if os.path.splitext(path)[1].lower() not in (".fit", ".fits", ".fts"):
        return None
    from eclipseforgehdr.fits import read_fits
    d, h = read_fits(path)
    if d.ndim != 3:
        return None                      # a mono/CFA FITS: the normal path
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


def _load(path, scale, demosaic_method="mhc"):
    """(3, H, W) float32 in raw ADU, pedestal NOT yet removed, plus sat_level."""
    direct = _planes_direct(path)
    if direct is not None:
        a, sat = direct
    else:
        from eclipseforgehdr.raw import open_frame, demosaic
        rf = open_frame(path)
        rgb = demosaic(rf.bayer, demosaic_method)      # H, W, 3
        a = np.ascontiguousarray(rgb.transpose(2, 0, 1), np.float32)
        sat = float(rf.sat_level)
    if scale > 1:
        h = (a.shape[1] // scale) * scale
        w = (a.shape[2] // scale) * scale
        a = a[:, :h, :w].reshape(3, h // scale, scale,
                                 w // scale, scale).mean((2, 4))
    return a, sat


def _pedestal(a, frac=0.04):
    """Per-channel median of the four corners: black level AND sky, together.

    Taken on the DEMOSAICED frame rather than the CFA so each channel gets its
    own number -- which is the whole point, since a colour cast that changes
    with radius is exactly a per-channel additive term.
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
    """Box high-pass. Cheap, separable, and good enough to align on."""
    k = int(max(1, round(s)))
    c = np.cumsum(np.pad(x, ((k, k), (k, k)), mode="edge"), axis=0)
    c = (c[2 * k:, :] - c[:-2 * k, :]) / (2 * k)
    c = np.cumsum(c, axis=1)
    c = (c[:, 2 * k:] - c[:, :-2 * k]) / (2 * k)
    return x - c[:x.shape[0], :x.shape[1]]


def _feature(g, floor=30.0):
    return _hipass(np.log1p(np.clip(g, 0, None) / floor), 24.0)


def _pcorr(a, b):
    """Sub-pixel shift taking b onto a. Plain phase correlation.

    NOT the pipeline's `_semi_phase_shift`, deliberately: this file is supposed
    to share as little as possible with the thing it is a control for.
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


def _shift(x, dy, dx):
    """Bilinear shift; outside the source becomes NaN so it can be counted out."""
    H, W = x.shape
    yy = np.arange(H, dtype=np.float32) - dy
    xx = np.arange(W, dtype=np.float32) - dx
    y0 = np.floor(yy).astype(np.int32)
    x0 = np.floor(xx).astype(np.int32)
    fy = (yy - y0)[:, None].astype(np.float32)
    fx = (xx - x0)[None, :].astype(np.float32)
    y0c = np.clip(y0, 0, H - 2)
    x0c = np.clip(x0, 0, W - 2)
    out = (x[y0c][:, x0c] * (1 - fy) * (1 - fx)
           + x[y0c][:, x0c + 1] * (1 - fy) * fx
           + x[y0c + 1][:, x0c] * fy * (1 - fx)
           + x[y0c + 1][:, x0c + 1] * fy * fx)
    out[(yy < 0) | (yy > H - 1), :] = np.nan
    out[:, (xx < 0) | (xx > W - 1)] = np.nan
    return out


# ------------------------------------------------------------------ merge ---

def _ratio(a, b, sat, lo=200.0):
    """Median b/a over pixels well exposed in BOTH. The measured exposure ratio.

    Measured on green, and the SAME number is applied to all three channels. On
    Val's 23 tiers the three channels' own ratios agree to under 1.2% per step,
    so the colour freedom buys nothing and costs a noisier estimate.
    """
    m = (np.isfinite(a) & np.isfinite(b) & (a > lo) & (b > lo)
         & (a < sat * 0.8) & (b < sat * 0.8))
    if m.sum() < 2000:
        return np.nan, int(m.sum())
    return float(np.median(b[m] / a[m])), int(m.sum())


def _weight(v, sat):
    """Hat: down at the noise floor, down towards saturation.

    Weighting by exposure time instead (what the shipped merge does) gives a
    short tier real weight in the outer field where it holds nothing but noise
    and whatever the pedestal step missed. This does not.
    """
    return (np.clip((v - 60.0) / 400.0, 0, 1)
            * np.clip((sat * 0.97 - v) / (sat * 0.12), 0, 1))


# ------------------------------------------------------------------- main ---

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("folder")
    ap.add_argument("--scale", type=int, default=2,
                    help="block-average by this factor first (default 2; "
                         "1 = full resolution, slower and much hungrier)")
    ap.add_argument("--out", default=None,
                    help="output directory (default FOLDER/simple_stack)")
    ap.add_argument("--max-frames", type=int, default=0,
                    help="use at most N frames per tier (0 = all)")
    ap.add_argument("--demosaic", default="mhc", choices=("mhc", "vng"))
    a = ap.parse_args(argv)

    from eclipseforgehdr.raw import list_raws, read_exif
    # expanduser BEFORE abspath: a quoted "~/..." reaches us verbatim, because
    # the shell only expands an unquoted tilde, and abspath would then glue it
    # onto the current directory and report a baffling missing path.
    folder = os.path.abspath(os.path.expanduser(a.folder))
    out = (os.path.abspath(os.path.expanduser(a.out)) if a.out
           else os.path.join(folder, "simple_stack"))
    os.makedirs(out, exist_ok=True)

    files = list_raws(folder)
    if not files:
        print("no readable frames in", folder)
        return 2

    # ---- group by exposure. The HEADER value is used ONLY as a grouping key;
    # every photometric number below comes from the pixels. On a third tester's set
    # the headers are Siril "manual group assignment" labels stepping a uniform
    # 1.5x where the pixels measure 1.19 / 1.49 / 1.99 -- so as a label they are
    # perfectly good and as an exposure they are fiction.
    tiers = {}
    for p in files:
        try:
            sec = float(read_exif(p)[0])
        except Exception as e:
            print("  skipping %s (%s)" % (os.path.basename(p), e))
            continue
        tiers.setdefault(sec, []).append(p)
    secs = sorted(tiers)
    print("%d frames, %d exposure groups" % (len(files), len(secs)))
    if len(secs) < 2:
        print("need at least two exposure groups to merge")
        return 2

    t0 = time.time()
    ref_feat = None
    stacks, sat = [], None
    for i, s in enumerate(secs):
        fs = sorted(tiers[s])
        if a.max_frames:
            fs = fs[:a.max_frames]
        acc = cnt = None
        shifts = []
        for p in fs:
            f, sl = _load(p, a.scale, a.demosaic)
            sat = sl if sat is None else sat
            f -= _pedestal(f)[:, None, None]
            ft = _feature(f[1])
            if ref_feat is None:
                ref_feat = ft
                dy = dx = 0.0
            else:
                dy, dx = _pcorr(ref_feat, ft)
                if not np.isfinite([dy, dx]).all() or max(abs(dy), abs(dx)) > 200:
                    dy = dx = 0.0
            shifts.append((float(dy), float(dx)))
            g = np.stack([_shift(f[c], dy, dx) for c in range(3)])
            ok = np.isfinite(g[0]).astype(np.float32)
            g = np.nan_to_num(g)
            acc = g if acc is None else acc + g
            cnt = ok if cnt is None else cnt + ok
        st = acc / np.maximum(cnt, 1.0)
        st[:, cnt < 0.5] = np.nan
        stacks.append(st.astype(np.float32))
        m = np.array(shifts)
        print("  %10.5f s  %2d frames  shift rms %.2f px  max %.2f  [%.0fs]"
              % (s, len(fs), float(np.sqrt((m ** 2).sum(1).mean())),
                 float(np.abs(m).max()), time.time() - t0))
        sys.stdout.flush()

    # ---- the ladder, measured ------------------------------------------
    print()
    print("exposure ratios MEASURED between adjacent groups (green):")
    print("  from        to          header    measured     px")
    rel = [1.0]
    for i in range(len(stacks) - 1):
        r, n = _ratio(stacks[i][1], stacks[i + 1][1], sat)
        hdr = secs[i + 1] / secs[i]
        if not np.isfinite(r):
            print("  %9.5f -> %9.5f  %8.3f   too few shared px" % (secs[i], secs[i + 1], hdr))
            r = hdr
        else:
            print("  %9.5f -> %9.5f  %8.3f  %8.3f  %9d"
                  % (secs[i], secs[i + 1], hdr, r, n))
        rel.append(rel[-1] * r)
    rel = np.array(rel)
    print()
    print("  ladder: header says 1 .. %.4g, the pixels say 1 .. %.4g"
          % (secs[-1] / secs[0], rel[-1]))

    # ---- merge ----------------------------------------------------------
    num = np.zeros_like(stacks[0], np.float32)
    den = np.zeros_like(stacks[0], np.float32)
    for i, st in enumerate(stacks):
        for c in range(3):
            v = st[c]
            ok = np.isfinite(v)
            w = np.where(ok, _weight(np.nan_to_num(v), sat), 0.0).astype(np.float32)
            num[c] += w * np.nan_to_num(v) / np.float32(rel[i])
            den[c] += w
    hdr = np.where(den > 0, num / np.maximum(den, 1e-9), np.nan)
    np.save(os.path.join(out, "hdr_simple.npy"), hdr)
    cov = float((den > 0).all(0).mean())
    print("  merged; %.2f%% of pixels have at least one usable tier" % (100 * cov))

    _preview(hdr, os.path.join(out, "simple_preview.png"))
    fits_path = os.path.join(out, "simple_linear.fits")
    off = _write_fits(hdr, fits_path)
    print()
    print("wrote", out)
    print("  hdr_simple.npy      linear, (3, H, W) float32, TRUE zero point")
    print("  simple_preview.png  asinh stretch, for the eye only")
    if off is not None:
        print("  simple_linear.fits  the same data + %.3f, for EFHDR's "
              "\"import one HDR\" field" % off)
        print()
        print("TO USE IT IN THE APP: in EFHDR, leave the raw folder empty, put "
              "this FITS in")
        print("the \"or import one HDR\" field, and Start. You get MGN, FNRGF, "
              "NAFE, Pellett")
        print("and every render control, applied to this merge instead of the "
              "app's own.")
    return 0


def _write_fits(hdr, path):
    """Write the merge as a linear 3-plane FITS for `importhdr`.

    WITH A SMALL POSITIVE OFFSET, deliberately. Once the sky is subtracted the
    far field straddles zero, and importhdr clips negatives:

        a = np.clip(np.nan_to_num(np.asarray(a, np.float32)), 0.0, None)

    docs/SKY_SUBTRACTION.md is explicit about what that does -- "Rectifying the
    noise gives each channel a mean proportional to its own sigma, which invents
    colour out of nothing." So we lift the whole image by a few times the far
    field's own noise before writing, which costs nothing against a corona
    thousands of times larger and keeps the noise symmetric where it matters.
    The .npy keeps the true zero point; this file is for the app.
    """
    try:
        from astropy.io import fits as _af
    except ImportError:
        print("  (no astropy, skipping the FITS for import)")
        return None
    a = np.nan_to_num(np.asarray(hdr, np.float32), nan=0.0)
    lum = a.mean(0)
    lo = lum[np.isfinite(lum)]
    sig = 1.4826 * float(np.median(np.abs(lo - np.median(lo)))) if lo.size else 0.0
    off = float(4.0 * sig)
    if not np.isfinite(off) or off <= 0:
        off = float(max(-a.min(), 0.0)) + 1.0
    _af.PrimaryHDU((a + off).astype(np.float32)).writeto(path, overwrite=True)
    return off


def _preview(hdr, path):
    """An asinh stretch and nothing else. This is a look, not a render."""
    try:
        from PIL import Image
    except ImportError:
        print("  (no PIL, skipping the preview)")
        return
    lum = np.nanmean(hdr, axis=0)
    hi = np.nanpercentile(lum, 99.8)
    if not np.isfinite(hi) or hi <= 0:
        return
    f = np.arcsinh(np.clip(lum / hi, 0, None) * 60) / np.arcsinh(60)
    ch = [hdr[c] / np.maximum(lum, 1e-9) for c in range(3)]
    rgb = np.clip(np.nan_to_num(np.stack([f * ch[c] for c in range(3)])), 0, 1)
    Image.fromarray((rgb.transpose(1, 2, 0) * 255 + 0.5).astype(np.uint8)).save(path)


if __name__ == "__main__":
    sys.exit(main())
