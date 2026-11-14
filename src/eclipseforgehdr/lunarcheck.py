"""Is the structure on the lunar disc the Moon's own? (2026-11-14)

The Moon layer (moonlayer.py) lifts the earthshine: the light and dark of the
lunar surface, lit by the Earth, on a disc that is otherwise scattered light.
What it lifts has to be the Moon's surface and not glare, noise or the
sensor. Until now the only test was internal: the layer built from the first
and from the second half of the frames had to agree. On one test set (the
Moon moving 12 px between the halves) that test kept 22 % of maria that
are clearly there; on another (four frames, two brackets) it accepted 89 % of
a pattern that is not the Moon at all.

THIS MODULE CHECKS THE PATTERN AGAINST A MAP OF THE REAL MOON. The map is
used for that check and for nothing else: nothing from it is drawn into the
picture, and the layer that is shown is always the user's own data. It only
decides whether the user's own structure is let into the result, and how
strongly.

THE MAP. lunar_albedo.png: a grayscale reduction (720 x 360 px, plate carree,
longitude -180..180 left to right, latitude +90..-90 top to bottom) of the
Lunar Reconnaissance Orbiter Camera colour map of the "CGI Moon Kit" of
NASA's Scientific Visualization Studio (lroc_color_poles_1k). NASA imagery,
public domain.

THE TEST.
  1. The map is drawn as the near side seen from the Earth, for every
     position angle (3 deg steps), both ways round (a mirrored picture) and
     librations of -6, 0 and +6 deg in longitude and latitude.
  2. Each drawing and the user's layer are band-passed alike (the maria
     scales) and compared inside 0.85 R by correlation; the best match counts.
  3. The same search is run on random patterns that have exactly the user's
     layer's texture (its power spectrum, random phases). That is how high a
     best match climbs by chance over so many trials.
  4. Confirmed: the best match stands at least CONFIRM_SD standard deviations
     above the random patterns and above every one of them. Rejected: it
     does not rise above their 95th percentile. In between: not decided.
Measured on three test sets (their Moon layers as built before this check):
     test set A   best 0.76, random 0.32 +- 0.06 (max 0.45)  -> 6.9 sd, confirmed
     test set B   best 0.67, random 0.31 +- 0.06 (max 0.41)  -> 6.3 sd, confirmed
     test set C   best 0.43, random 0.35 +- 0.06 (max 0.46)  -> 1.4 sd, rejected
and the matched position angle of test set B lies 16 deg from the
direction to celestial north measured from its stars, within the lunar axis'
own tilt.
"""
import os

import numpy as np
from scipy import ndimage

MAP_FILE = "lunar_albedo.png"
MAP_SOURCE = "LRO camera colour map, NASA Scientific Visualization Studio (CGI Moon Kit)"
GRID = 121                 # px across the disc for the comparison
INNER = 0.85               # compared inside this radius (R)
STEP_DEG = 3
LIBRATION = (-6, 0, 6)
N_RANDOM = 60
CONFIRM_SD = 4.0
_BANK = {}


def _load_map():
    from PIL import Image
    p = os.path.join(os.path.dirname(__file__), MAP_FILE)
    return np.asarray(Image.open(p).convert("L"), np.float64) / 255.0


def _render(A, n, l_deg, b_deg, mirror):
    """The map seen as a disc n x n px: north up, the IAU east limb to the
    right, the sub-observer point at latitude b, longitude l."""
    H, W = A.shape
    yy, xx = np.mgrid[0:n, 0:n].astype(np.float64)
    c = (n - 1) / 2.0
    y = -(yy - c) / c
    x = (xx - c) / c
    if mirror:
        x = -x
    rho2 = x * x + y * y
    inside = rho2 < 1
    z = np.sqrt(np.clip(1 - rho2, 0, 1))
    b, l = np.radians(b_deg), np.radians(l_deg)
    Y1 = y * np.cos(b) + z * np.sin(b)
    Z1 = -y * np.sin(b) + z * np.cos(b)
    X2 = x * np.cos(l) + Z1 * np.sin(l)
    Z2 = -x * np.sin(l) + Z1 * np.cos(l)
    lat = np.degrees(np.arcsin(np.clip(Y1, -1, 1)))
    lon = np.degrees(np.arctan2(X2, Z2))
    py = (90 - lat) / 180 * H - 0.5
    px = (lon + 180) / 360 * W - 0.5
    v = ndimage.map_coordinates(A, [py, px % W], order=1, mode="nearest")
    v[~inside] = float(np.mean(v[inside]))
    return v


def _geometry(n=GRID):
    yy, xx = np.mgrid[0:n, 0:n]
    r = np.hypot(yy - (n - 1) / 2, xx - (n - 1) / 2) / ((n - 1) / 2)
    return r, r < INNER


def _prep(d, r, m, n=GRID):
    d = d.copy()
    d[r >= 0.95] = float(np.mean(d[m]))
    x = ndimage.gaussian_filter(d, 1.2) - ndimage.gaussian_filter(d, n * 0.15)
    x = x[m]
    s = float(x.std())
    return (x - x.mean()) / (s if s > 0 else 1.0)


def _bank():
    """Every drawing of the map, band-passed: (matrix, keys). Built once."""
    if "M" in _BANK:
        return _BANK["M"], _BANK["K"]
    A = _load_map()
    r, m = _geometry()
    T, K = [], []
    for mir in (False, True):
        for l in LIBRATION:
            for b in LIBRATION:
                base = _render(A, GRID, l, b, mir)
                for P in range(0, 360, STEP_DEG):
                    t = ndimage.rotate(base, P, reshape=False, order=1, mode="nearest")
                    T.append(_prep(t, r, m))
                    K.append((P, l, b, mir))
    _BANK["M"], _BANK["K"] = np.stack(T).astype(np.float32), K
    return _BANK["M"], _BANK["K"]


def disc_from(layer, cy, cx, R):
    """The layer's disc resampled onto the comparison grid."""
    n = GRID
    yy, xx = np.mgrid[0:n, 0:n].astype(np.float64)
    u = (yy - (n - 1) / 2) / ((n - 1) / 2)
    v = (xx - (n - 1) / 2) / ((n - 1) / 2)
    k = max(R / n * 0.7, 0.5)
    a = ndimage.gaussian_filter(np.asarray(layer, np.float64), k)
    return ndimage.map_coordinates(a, [cy + u * R, cx + v * R], order=1, mode="nearest")


def _scores(disc, seed):
    """The best match of a disc and the best matches of N_RANDOM random
    patterns with its texture: (best, index of the best drawing, array)."""
    M, _ = _bank()
    r, m = _geometry()
    c = M @ _prep(disc, r, m).astype(np.float32) / m.sum()
    i = int(np.argmax(c))
    rng = np.random.default_rng(seed)
    base = np.asarray(disc, np.float64).copy()
    base[r >= 0.95] = float(np.mean(base[m]))
    F = np.abs(np.fft.rfft2(base))
    rand = []
    for _ in range(N_RANDOM):
        ph = np.exp(1j * rng.uniform(0, 2 * np.pi, F.shape))
        ph[0, 0] = 1
        s = np.fft.irfft2(F * ph, s=base.shape)
        rand.append(float(np.max(M @ _prep(s, r, m).astype(np.float32) / m.sum())))
    return float(c[i]), i, np.asarray(rand)


def _verdict(best, i, rand):
    _, K = _bank()
    mu, sd = float(rand.mean()), float(rand.std()) or 1e-9
    z = (best - mu) / sd
    p95 = float(np.percentile(rand, 95))
    if z >= CONFIRM_SD and best > rand.max():
        verdict = "confirmed"
    elif best <= p95:
        verdict = "rejected"
    else:
        verdict = "undecided"
    P, l, b, mir = K[i]
    return {"best": round(best, 3), "pa_deg": P, "libration": [l, b], "mirrored": bool(mir),
            "random_mean": round(mu, 3), "random_sd": round(sd, 3),
            "random_max": round(float(rand.max()), 3), "random_p95": round(p95, 3),
            "sd": round(z, 1), "verdict": verdict, "map": MAP_SOURCE}


def check(disc, seed=0):
    """The verdict on a disc (from disc_from): dict with best (correlation),
    pa_deg / libration / mirrored (the best drawing), random_mean,
    random_sd, random_max, random_p95, sd (how far the best stands out) and
    verdict "confirmed" | "rejected" | "undecided"."""
    return _verdict(*_scores(disc, seed))


def check_halves(disc_a, disc_b, seed=0):
    """The better of two half-stacks, judged fairly. Taking the better of two
    gives chance two tries, so each random trial counts the better of its two
    random patterns as well (one with each half's texture). Returns the
    verdict dict of the better half plus "half": 0 or 1, and each half's own
    best match ("best_halves")."""
    ba, ia, ra = _scores(disc_a, seed)
    bb, ib, rb = _scores(disc_b, seed + 1)
    h, best, i = (0, ba, ia) if ba >= bb else (1, bb, ib)
    out = _verdict(best, i, np.maximum(ra, rb))
    out["half"] = h
    out["best_halves"] = [round(ba, 3), round(bb, 3)]
    return out


def share_for(best):
    """How strongly confirmed maria go into the result (0..1): the match
    against the real Moon, 0.75 and better counting as full strength."""
    return float(np.clip(best / 0.75, 0.5, 1.0))
