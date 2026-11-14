"""Smoke test: the per-tier clipped mean must not leave a trail where an
object crossed the field (pipeline._grow_rejection, lab 2026-11-14).

Twelve synthetic frames of sky noise; in each, a soft blob (core far above
the noise, a wide faint halo inside it) at another place along a line. The
clipped mean without growth leaves the halos; with growth the residual along
the path has to be at the noise level, and on frames without any outlier the
result has to equal the plain mean.
"""
import sys, os
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from eclipseforgehdr import pipeline as P

rng = np.random.default_rng(1)
H, W, N = 600, 400, 12
yy, xx = np.mgrid[:H, :W]
frames = []
for i in range(N):
    f = 1000 + rng.normal(0, 30, (H, W))
    cy, cx = 60 + i * 42, 120 + i * 12
    r2 = (yy - cy) ** 2 + (xx - cx) ** 2
    f += 20000 * np.exp(-r2 / (2 * 2.0 ** 2)) + 70 * np.exp(-r2 / (2 * 12.0 ** 2))
    frames.append(f.astype(np.float32))
path = np.zeros((H, W), bool)
for i in range(N):
    cy, cx = 60 + i * 42, 120 + i * 12
    path[max(cy - 30, 0):cy + 30, max(cx - 30, 0):cx + 30] = True
from scipy import ndimage
# judged smoothed over a few px, as the detail layers see it
_ref = ndimage.gaussian_filter(np.mean([(1000 + rng.normal(0, 30, (H, W))) for _ in range(N)], axis=0), 3.0)
noise = float(_ref.std())


def worst(out):
    return float(np.abs(ndimage.gaussian_filter(out.astype(np.float64), 3.0)[path] - 1000).max()) / noise


seed = P.GROW_SEED
P.GROW_SEED = 1e9
old, _ = P._clipped_mean(frames, 3.0)
P.GROW_SEED = seed
new, nrej = P._clipped_mean(frames, 3.0)
print(f"path residual: without growth {worst(old):.1f} sd, with growth {worst(new):.1f} sd")
assert worst(old) > 6.0, "the test object should show without growth"
assert worst(new) < 5.0, "growth should take the object out"

clean = [(1000 + rng.normal(0, 30, (H, W))).astype(np.float32) for _ in range(N)]
c_out, _ = P._clipped_mean(clean, 3.0)
d = float(np.abs(c_out - np.mean(clean, axis=0)).max()) / (30 / np.sqrt(N))
print(f"clean frames: largest difference from the plain mean {d:.2f} sd")
assert d < 3.0
print("ok")
