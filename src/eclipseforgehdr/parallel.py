"""Threads for the independent parts of a stack (lab 0.28.x, 2026-10-18).

The heavy steps are numpy and scipy.ndimage calls, which release the GIL, so
plain threads give most of the speed of several cores with NO copy of the
arrays: the tiers' demosaic and shift, the per-frame hot-pixel test, the raw
decode, and the detail layers' bands and scales all run in a small pool.

Measured on two cores (cloud, 2026-10-18): demosaic + shift of two tiers in
4.8 s against 7.9 s one after the other; a sigma-12 Gaussian blur 1.0 s
against 1.4 s.

How many: ECLIPSEFORGE_THREADS when set, otherwise the smaller of the core
count and one worker per 8 GB of memory, at most 4 -- the memory cap because
every worker holds a tier or a band of its own in flight (a 45 Mpx tier
demosaiced is 540 MB), and the README's budget of 2 GB per 4 tiers has to
stay true on an 8 GB machine.  1 means: exactly the old sequential code.
"""
import os
import numpy as np
from concurrent.futures import ThreadPoolExecutor

_CACHE = {}


def _total_ram_gb():
    try:
        return os.sysconf("SC_PHYS_PAGES") * os.sysconf("SC_PAGE_SIZE") / 2.0 ** 30
    except (ValueError, OSError, AttributeError):
        return None


def workers(cap=4, gb_per_worker=8.0):
    """Pool size for this process; cached. `gb_per_worker` 0 = the step
    holds nothing extra per thread (a strip of a blur), so only the cores
    count."""
    key = ("w", cap, gb_per_worker)
    if key in _CACHE:
        return _CACHE[key]
    v = os.environ.get("ECLIPSEFORGE_THREADS", "").strip()
    if v:
        try:
            n = max(1, int(v))
        except ValueError:
            n = 1
    else:
        n = os.cpu_count() or 1
        ram = _total_ram_gb()
        if ram and gb_per_worker > 0:
            n = min(n, max(1, int(ram // gb_per_worker)))
        n = max(1, min(cap, n))
    _CACHE[key] = n
    return n


def tmap(fn, items, cap=4):
    """map() over a thread pool, results in order; sequential when 1 worker."""
    items = list(items)
    n = min(workers(cap), len(items))
    if n <= 1:
        return [fn(x) for x in items]
    with ThreadPoolExecutor(max_workers=n) as ex:
        return list(ex.map(fn, items))


class Prefetch:
    """Compute items ahead of a sequential consumer, `depth` at a time.

    for s in secs: use(pf.get(s))  -- get() returns item s and starts the
    next one(s) in the background. depth 0 = compute on demand."""

    def __init__(self, fn, keys, depth=None, cap=4):
        self.fn, self.keys = fn, list(keys)
        d = (workers(cap) - 1) if depth is None else depth
        self.depth = max(0, min(d, 2, len(self.keys) - 1))
        self.ex = ThreadPoolExecutor(max_workers=max(1, self.depth)) if self.depth else None
        self.fut, self.next = {}, 0
        self._fill()

    def _fill(self):
        while self.ex is not None and self.next < len(self.keys) and len(self.fut) < self.depth + 1:
            k = self.keys[self.next]; self.next += 1
            self.fut[k] = self.ex.submit(self.fn, k)

    def get(self, k):
        if self.ex is None:
            return self.fn(k)
        f = self.fut.pop(k, None)
        if f is None:
            f = self.ex.submit(self.fn, k)
        v = f.result()
        self._fill()
        return v

    def close(self):
        if self.ex is not None:
            for f in self.fut.values():
                f.cancel()
            self.ex.shutdown(wait=False)
            self.ex = None


def gaussian_filter(a, sigma, order=0, output=None, mode="reflect", cval=0.0,
                    truncate=4.0, radius=None):
    """scipy.ndimage.gaussian_filter, the same arithmetic, on several threads.

    scipy runs one 1-D pass per axis, and every line of a pass is independent
    of the others, so each pass is cut into strips along the OTHER axis and
    the strips run in parallel. Each element goes through exactly the 1-D
    correlations scipy would give it, in the same order and dtype, so the
    result is bit-identical to scipy's. Falls back to scipy for anything
    unusual (not 2-D, order > 0, a 1-worker pool, or a small image)."""
    from scipy import ndimage
    a = np.asarray(a)
    n = workers(8, 0.0)
    if a.ndim != 2 or order != 0 or n <= 1 or a.size < 1_000_000 or radius is not None:
        return ndimage.gaussian_filter(a, sigma, order=order, output=output, mode=mode,
                                       cval=cval, truncate=truncate)
    sig = sigma if np.iterable(sigma) else (sigma, sigma)
    if len(sig) != 2:
        return ndimage.gaussian_filter(a, sigma, order=order, output=output, mode=mode,
                                       cval=cval, truncate=truncate)
    modes = mode if (isinstance(mode, (list, tuple))) else (mode, mode)
    if output is None:
        out = np.empty(a.shape, dtype=a.dtype)          # scipy: output dtype = input dtype
    else:
        out = output
    src = a
    for axis in (0, 1):
        s = float(sig[axis])
        if s <= 0.0:
            if src is not out:
                out[...] = src
                src = out
            continue
        other = 1 - axis
        L = src.shape[other]
        k = min(n, max(1, L // 256))
        edges = np.linspace(0, L, k + 1).astype(int)
        def _strip(i, _axis=axis, _s=s, _src=src, _mode=modes[axis]):
            sl = [slice(None), slice(None)]
            sl[other] = slice(int(edges[i]), int(edges[i + 1]))
            sl = tuple(sl)
            ndimage.gaussian_filter1d(_src[sl], _s, axis=_axis, output=out[sl],
                                      mode=_mode, cval=cval, truncate=truncate)
        if k > 1:
            with ThreadPoolExecutor(max_workers=k) as ex:
                list(ex.map(_strip, range(k)))
        else:
            _strip(0)
        src = out
    return out
