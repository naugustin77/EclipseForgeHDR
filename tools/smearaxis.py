#!/usr/bin/env python3
"""Is the merged corona smeared along ONE direction, and which?

The question this answers: a merge built from tiers that do not sit exactly on
each other is not blurred evenly -- the tiers are displaced along a particular
axis, so the merge loses structure along that axis and keeps it across. That
prints as a soft WEDGE in the picture, at two azimuths 180 degrees apart, and
it is easy to mistake for a coronal hole or for lunar motion.

THE MEASURE is cross-streamer structure over radial structure, per azimuth.
The corona is resampled to polar about the fitted disc, mildly smoothed to put
pixel noise below the streamers, and at each azimuth the mean |d/dphi| is
divided by the mean |d/dr|, both normalised by the local level. Streamers are
radial, so |d/dphi| is the contrast ACROSS them -- the thing a smear destroys
-- while |d/dr| is a reference that noise raises equally. The ratio is free of
the render's stretch, which is why it can be taken off a finished PNG.

    python tools/smearaxis.py --folder /path/to/bracket [--lab]
    python tools/smearaxis.py --folder /path --against /other/path
    python tools/smearaxis.py a.png b.png        # two renders of the same set

ONE RUN gives the azimuth profile plus, from the work directory, the principal
axis of the per-tier shifts and how elongated they are. A set whose tiers are
displaced along one axis is a set whose merge can only be smeared along it, so
that axis is the prediction and the profile is the evidence.

TWO RUNS is the clean test, and the one to use when judging a change. The
corona's own structure is not uniform in azimuth -- there are genuinely more
streamers in some directions -- and that swamps a single profile. Dividing one
run's profile by the other's cancels it exactly, because both runs are the same
corona. What survives is the difference in smear, and a one-axis smear appears
in it as a cos(2*phi): destroyed where the smear runs across the streamers,
apparently IMPROVED 90 degrees away, where the same smear removes radial
structure instead and lifts the ratio. The tool fits that cos(2*phi) and
reports its amplitude and axis. Amplitude near zero means the two runs are
smeared the same way; a large amplitude means one of them is smeared along the
reported axis, and the sign says which.

WORKED EXAMPLE (a 360 mm round-robin bracket, 0.23.8). Moon Align against Corona
Align: the two runs measure the same overall (+0.1%), but the ratio fits a
cos(2*phi) of amplitude 18.8% about an axis of 74 degrees, and the per-tier
shifts in the work directory lie along 82 degrees with 46x elongation. Same
axis: Moon Align smears each tier's corona along the Moon's drift, Corona Align
does not. The soft wedge at 330-350 degrees survives in BOTH runs, weaker in
Corona Align, because what remains there is the cross-tier residual -- 9.65 px
full-res along that same axis -- which neither lock touches. That residual is
what a cross-tier change has to reduce, and this tool is how to tell whether
it did, on real data rather than on a synthetic.

Needs numpy, scipy and Pillow -- run it with the same python the app uses.
"""
import os
import sys
import json

try:
    import numpy as np
    from scipy import ndimage, optimize
except ImportError as _e:                        # noqa: BLE001
    # The app is installed with pipx, so its dependencies live in pipx's venv
    # and NOT in the system python. Name the interpreter that has them rather
    # than leaving a bare ImportError: on macOS `python` does not exist at all
    # and `python3` is the system one, which has no scipy.
    _v = os.path.expanduser("~/.local/pipx/venvs/eclipseforgehdr/bin/python")
    sys.stderr.write(
        f"{__file__.split(os.sep)[-1]}: {_e}\n\n"
        "This needs numpy, scipy and Pillow. They are in the venv pipx made\n"
        "for the app, not in the system python. Run it with:\n\n"
        f"    {_v} {os.path.relpath(__file__)} ...\n\n"
        "If that path does not exist, ask pipx where it put them:\n\n"
        "    ls \"$(pipx environment --value PIPX_LOCAL_VENVS)\"\n")
    raise SystemExit(2)

R_IN, R_OUT = 1.03, 1.45          # the annulus that carries streamer structure
NA = 720                          # azimuth samples before binning
BIN = 10                          # degrees per reported bin


def _load_image(path):
    """Luminance from a render, or the merged luminance from a work directory."""
    if path.endswith(".npy"):
        return np.asarray(np.load(path, mmap_mode="r")[:], np.float32)
    from PIL import Image
    a = np.asarray(Image.open(path).convert("RGB"), np.float32)
    return 0.2126 * a[:, :, 0] + 0.7152 * a[:, :, 1] + 0.0722 * a[:, :, 2]


def _fit_disc(L):
    """Centre and radius of the occulted disc, fitted to its boundary.

    Fitted rather than centroided because the disc often runs off the frame in
    a crop, and a centroid of a clipped circle is not its centre. Boundary
    points on the frame edge are dropped for the same reason.

    The threshold is taken in the LOG, which is what makes the same rule work
    on a finished render and on linear merged data. On the linear merge a
    threshold set as a fraction of the range is meaningless -- the corona spans
    seven stops, so a few percent of the maximum still sits far above the
    faint outer field, and a first version of this returned R = 994 px where
    the pipeline's own fit says 456. Prefer geometry.json whenever there is one
    (see _resolve); this is the fallback for a bare image.
    """
    Ln = np.log1p(np.maximum(L, 0) / max(np.percentile(L, 60), 1e-6))
    lo, hi = np.percentile(Ln, 2), np.percentile(Ln, 99)
    dark = ndimage.binary_opening(Ln < lo + 0.10 * (hi - lo), np.ones((9, 9)))
    lab, n = ndimage.label(dark)
    if n == 0:
        raise RuntimeError("no occulted disc found")
    sz = ndimage.sum(dark, lab, range(1, n + 1))
    disc = lab == (int(np.argmax(sz)) + 1)
    H, W = L.shape
    edge = disc & ~ndimage.binary_erosion(disc)
    ey, ex = np.nonzero(edge)
    k = (ey > 2) & (ey < H - 3) & (ex > 2) & (ex < W - 3)
    ey, ex = ey[k], ex[k]
    if ey.size < 50:
        raise RuntimeError("the disc boundary is too short to fit")
    ys, xs = np.nonzero(disc)
    p = optimize.leastsq(lambda q: np.hypot(ey - q[0], ex - q[1]) - q[2],
                         [ys.mean(), xs.mean(), np.sqrt(disc.sum() / np.pi)])[0]
    return float(p[0]), float(p[1]), float(p[2])


def profile(L, cy, cx, R):
    """Cross-streamer over radial structure, per azimuth bin (0 deg = up, cw)."""
    rr = np.linspace(R_IN * R, R_OUT * R, 300)
    ph = np.radians(np.linspace(0, 360, NA, endpoint=False))
    P = ndimage.map_coordinates(
        L, [cy - rr[:, None] * np.cos(ph)[None, :],
            cx + rr[:, None] * np.sin(ph)[None, :]],
        order=1, mode="constant", cval=np.nan)
    ok = np.isfinite(P) & (P > 0.02 * np.nanpercentile(P, 99))
    Ps = ndimage.gaussian_filter(np.nan_to_num(P), (1.0, 1.0))
    gt = np.abs(np.gradient(Ps, axis=1))
    gr = np.abs(np.gradient(Ps, axis=0))
    out = np.full(360 // BIN, np.nan)
    for i, a0 in enumerate(range(0, 360, BIN)):
        s = slice(int(a0 / 360 * NA), int((a0 + BIN) / 360 * NA))
        m = ok[:, s]
        if m.sum() < 400:
            continue
        lvl = np.maximum(Ps[:, s][m], 1e-6)
        out[i] = (gt[:, s][m] / lvl).mean() / (gr[:, s][m] / lvl).mean()
    return out


def fit_two_fold(ratio):
    """Fit A*cos(2*(phi - axis)) to a log ratio. Returns (amplitude %, axis deg).

    Two-fold because a smear has no head or tail: it destroys structure the
    same way at phi and at phi+180. Fitting in the log makes the amplitude a
    symmetric percentage rather than one that depends on which run is on top.
    """
    a0 = np.radians(np.arange(0, 360, BIN) + BIN / 2.0)
    y = np.log(ratio)
    good = np.isfinite(y)
    if good.sum() < 8:
        return float("nan"), float("nan")
    A = np.stack([np.ones(good.sum()), np.cos(2 * a0[good]), np.sin(2 * a0[good])], 1)
    c = np.linalg.lstsq(A, y[good], rcond=None)[0]
    amp = float(np.hypot(c[1], c[2]))
    axis = float((np.degrees(np.arctan2(c[2], c[1])) / 2.0) % 180.0)
    return 100.0 * (np.exp(amp) - 1.0), axis


def tier_axis(wd):
    """Principal axis of the per-tier shifts, from the work directory.

    None when the input was a render rather than a folder -- the renders carry
    no shift table, so the prediction simply is not available and the profile
    has to stand on its own.
    """
    if not wd:
        return None
    p = os.path.join(wd, "geometry.json")
    if not os.path.exists(p):
        return None
    sh = (json.load(open(p)).get("abs_shift") or {})
    if len(sh) < 3:
        return None
    V = np.array([[sh[k][0], sh[k][1]] for k in sorted(sh, key=float)], float)
    V = V - V.mean(0)
    w, U = np.linalg.eigh(np.cov(V.T))
    pa = U[:, int(np.argmax(w))]
    axis = float((np.degrees(np.arctan2(pa[1], -pa[0]))) % 180.0)
    elong = float(np.sqrt(max(w) / max(min(w), 1e-9)))
    span = float(np.hypot(*(V.max(0) - V.min(0))))
    return axis, elong, span


def _disc_of(L, wd):
    """The disc the pipeline measured, or a fit when there is no work directory.

    Always prefer geometry.json: that circle was fitted on 720 rays and then
    corrected onto the tiers' own consensus, which is a better number than
    anything this tool would recover from a finished image.
    """
    if wd:
        p = os.path.join(wd, "geometry.json")
        if os.path.exists(p):
            g = json.load(open(p))
            return float(g["cy"]), float(g["cx"]), float(g["R"])
    return _fit_disc(L)


def _resolve(path, lab):
    """A folder becomes its merged luminance; a file is taken as given."""
    if os.path.isdir(path):
        wd = os.path.join(path, ".eclipseforgehdr-lab" if lab
                          else ".eclipseforgehdr")
        f = os.path.join(wd, "hdr_lum.npy")
        if not os.path.exists(f):
            raise RuntimeError(f"no merged luminance in {wd} — stack it first")
        return f, wd
    return path, None


def report_one(name, prof, wd):
    print(f"\n{name}")
    ta = tier_axis(wd) if wd else None
    if ta:
        axis, elong, span = ta
        print(f"   per-tier shifts lie along {axis:.0f} deg, {elong:.0f}x elongated, "
              f"spanning {span:.1f} px (half-res)")
        print(f"   -> if the merge is smeared it is smeared along {axis:.0f} deg, "
              f"which costs structure near {(axis + 90) % 180:.0f} and "
              f"{(axis + 270) % 360:.0f} deg")
    print(f"   {'azimuth':>9}{'cross/radial':>14}")
    med = np.nanmedian(prof)
    for i, a0 in enumerate(range(0, 360, BIN)):
        if not np.isfinite(prof[i]):
            continue
        d = 100 * (prof[i] / med - 1)
        bar = "#" * int(round(28 * prof[i] / (2 * med)))
        print(f"   {a0:>6}-{a0 + BIN:<3}{prof[i]:10.2f}{d:+7.0f}%  {bar}")
    print(f"   median {med:.2f}")


def main(argv):
    lab = False
    against = None
    paths = []
    i = 0
    while i < len(argv):
        a = argv[i]
        if a == "--lab":
            lab = True
        elif a in ("--folder", "--against"):
            i += 1
            if a == "--against":
                against = argv[i]
            else:
                paths.append(argv[i])
        elif a in ("-h", "--help"):
            print(__doc__)
            return 0
        else:
            paths.append(a)
        i += 1
    if against:
        paths.append(against)
    if not paths:
        print(__doc__)
        return 1

    runs = []
    for p in paths:
        f, wd = _resolve(p, lab)
        L = _load_image(f)
        cy, cx, R = _disc_of(L, wd)
        runs.append((os.path.basename(p.rstrip("/")) or p, profile(L, cy, cx, R), wd, R))
        print(f"{runs[-1][0]}: disc R = {R:.0f} px in {L.shape[1]}x{L.shape[0]}")

    for nm, prof, wd, _ in runs:
        report_one(nm, prof, wd)

    if len(runs) >= 2:
        (n1, p1, wd1, _), (n2, p2, wd2, _) = runs[0], runs[1]
        ratio = p2 / p1
        amp, axis = fit_two_fold(ratio)
        print(f"\n{n2} AGAINST {n1}")
        print(f"   overall          {100 * (np.nanmean(p2) / np.nanmean(p1) - 1):+.1f}%"
              f"   (a one-axis smear moves structure between directions, so this"
              f" stays near zero even when one run is smeared)")
        print(f"   two-fold term    {amp:.1f}% about an axis of {axis:.0f} deg")
        print(f"   {'azimuth':>9}{'ratio':>9}")
        for i, a0 in enumerate(range(0, 360, BIN)):
            if not np.isfinite(ratio[i]):
                continue
            print(f"   {a0:>6}-{a0 + BIN:<3}{100 * (ratio[i] - 1):+8.0f}%")
        ta = tier_axis(wd1) or tier_axis(wd2)
        if ta:
            d = abs(((ta[0] - axis + 90) % 180) - 90)
            print(f"\n   the per-tier shifts lie along {ta[0]:.0f} deg; the fitted "
                  f"smear axis is {axis:.0f} deg — {d:.0f} deg apart")
            print("   " + ("SAME AXIS: the difference between these runs is smear "
                           "along the direction the tiers are displaced"
                           if d < 20 else
                           "DIFFERENT AXES: whatever separates these runs is not "
                           "the tier displacement"))
        print("\n   Lower this two-fold term and the wedge goes away; that is the "
              "number a cross-tier change has to move.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
