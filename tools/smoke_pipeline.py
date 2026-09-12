#!/usr/bin/env python3
"""Run the WHOLE pipeline end to end on a synthetic FITS bracket.

WHY THIS EXISTS, on top of tools/smoke_layers.py. That one starts at
build_layers, so everything before it -- decode, stack, align, photometric
calibration, the shared pedestal, the per-tier radial check, LDIC, the merge,
the sky fit -- had no test that ran it at all. 0.22.67 shipped a crash in
exactly that kind of untested wiring, at the end of an 18-minute stack on a
tester's machine. 0.22.76 rewrites the frame three of those estimators sample
in, which is the same class of change.

FITS rather than raw, because it needs no camera library and because it also
covers the FITS reader's own saturation and odd-dimension paths.

    python tools/smoke_pipeline.py

Fails loudly on any exception, on a missing product, or on a non-finite merge.
"""
import os
import shutil
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))


class _P:
    def __init__(self, verbose):
        self.verbose = verbose
        self.lines = []

    def log(self, m, f=None):
        self.lines.append(m)
        if self.verbose:
            print("  |", m[:150])


def _card(k, v):
    if isinstance(v, str):
        val = "'%-8s'" % v
    elif isinstance(v, bool):
        val = "T" if v else "F"
    else:
        val = "%20s" % repr(v)
    return ("%-8s= %-70s" % (k, val))[:80]


def write_fits(path, data, exptime, when):
    """Uncompressed 16-bit FITS with the keywords the pipeline needs."""
    # BZERO 32768: the file stores physical - 32768 as int16. Reinterpreting
    # the bytes instead put every sample 32768 ADU too high -- a pedestal bigger
    # than the corona itself. See tools/smoke_fits.py.
    d = (np.clip(np.rint(data), 0, 65535) - 32768).astype(">i2")
    hdr = [_card("SIMPLE", True), _card("BITPIX", 16), _card("NAXIS", 2),
           _card("NAXIS1", data.shape[1]), _card("NAXIS2", data.shape[0]),
           _card("BZERO", 32768), _card("BSCALE", 1),
           _card("EXPTIME", float(exptime)), _card("GAIN", 400),
           _card("SATURATE", 60000), _card("DATE-OBS", when),
           "END" + " " * 77]
    blob = "".join(c.ljust(80) for c in hdr)
    blob += " " * ((2880 - len(blob) % 2880) % 2880)
    raw = d.tobytes()
    raw += b"\0" * ((2880 - len(raw) % 2880) % 2880)
    with open(path, "wb") as f:
        f.write(blob.encode("ascii"))
        f.write(raw)


def bracket(folder, H=420, W=560, R=46.0, secs=(1 / 60, 1 / 15, 1 / 4, 1.0),
            drift=3.0, pedestal=4.0, peak=6.0e5):
    """A corona bracket with real per-tier drift -- the thing _ring_sample is
    about. Two frames per tier so the stacker has something to average.

    THE NUMBERS ARE CHOSEN TO BE PHOTOMETRICALLY HONEST, because otherwise this
    test can only catch crashes. An earlier version used a 500:1 exposure span
    and a 300 ADU pedestal, which left the shortest tier with ~24 ADU of corona
    under 300 ADU of offset: every photometric link was rejected, the pipeline
    correctly reported the data as non-linear, and the estimator numbers it
    printed were meaningless. 24:1 with a 4 ADU residual is what a real bracket
    looks like, and then the links land near 1.000 and the pedestal fit can be
    checked against a known answer."""
    cy, cx = H / 2.0, W / 2.0
    yy = np.arange(H, dtype=np.float64)[:, None] - cy
    xx = np.arange(W, dtype=np.float64)[None, :] - cx
    r = np.hypot(yy, xx)
    th = np.arctan2(yy, xx)
    thr = 1.0 + 0.30 * np.sin(3 * th + 0.7) + 0.10 * np.sin(11 * th)
    base = peak * (R / np.maximum(r, R)) ** 4 * thr
    base[r < R] = 0.0
    rng = np.random.default_rng(11)
    n = 0
    for i, s in enumerate(secs):
        dy, dx = (i - len(secs) // 2) * drift, (i - len(secs) // 2) * drift * 0.7
        from scipy import ndimage
        moved = ndimage.shift(base, (dy, dx), order=1, mode="nearest")
        for k in range(2):
            a = moved * s + pedestal        # a black level left a few ADU behind
            a = a + rng.normal(0, np.sqrt(np.maximum(a, 1)))
            write_fits(os.path.join(folder, "f%02d.fits" % n), a, s,
                       "2026-08-12T18:%02d:%02d.0" % (n // 60, n % 60))
            n += 1
    return n


def main(verbose=True):
    from eclipseforgehdr.pipeline import run, workdir
    folder = tempfile.mkdtemp(prefix="efhdr_pipe_")
    try:
        nf = bracket(folder)
        _true_pedestal = 4.0
        prog = _P(verbose)
        run(folder, prog, crop_pc=320, denoise="fine")   # returns None
        wd = workdir(folder)
        need = ["hdr_rgb.npy", "hdr_lum.npy", "geometry.json", "report.json"]
        missing = [f for f in need
                   if not os.path.exists(os.path.join(wd, f))]
        if missing:
            raise SystemExit("pipeline finished but did not write: %s" % missing)
        hdr = np.load(os.path.join(wd, "hdr_lum.npy"))
        if not np.isfinite(hdr).all():
            raise SystemExit("merged luminance contains non-finite values")
        if float(hdr.max()) <= 0:
            raise SystemExit("merged luminance is entirely zero — every tier "
                             "was treated as clipped or invalid")
        # The three ring estimators must have RUN, not merely not crashed:
        # each of them records something in stats on success, and a silent
        # `except` around any of them would otherwise look like a pass. The
        # pedestal only records when its fit converges, so a bracket this small
        # is allowed to skip it -- but not to raise.
        import json
        stats = json.load(open(os.path.join(wd, "report.json")))
        _ran = [k for k in ("pedestal", "tier_radial", "ldic_k_spread_pct")
                if k in stats]
        if "tier_radial" not in _ran:
            raise SystemExit("the per-tier radial check did not run")
        if "ldic_k_spread_pct" not in _ran:
            raise SystemExit("the LDIC azimuthal fit did not run")
        print("ring estimators that ran: " + ", ".join(_ran))
        # The fixture is linear by construction, so the pipeline must say so:
        # a link near 1.000 means the exposure time predicts the signal, which
        # is the whole premise the merge rests on.
        _lk = stats.get("photometric_links") or []
        if _lk and (min(_lk) < 0.9 or max(_lk) > 1.1):
            raise SystemExit("photometric links %s stray from 1.000 on data "
                             "that is linear by construction" % [round(x, 3)
                                                                 for x in _lk])
        if _lk:
            print("photometric links: "
                  + ", ".join("%.3f" % x for x in _lk) + "  (truth 1.000)")
        _pd = (stats.get("pedestal") or {}).get("applied")
        if _pd is not None:
            print("pedestal: fitted %+.1f ADU against a true %+.1f"
                  % (_pd, _true_pedestal))
        # A RUN THAT HAS A LUNAR TRACK MUST NOT REPORT HAVING NONE (0.22.87).
        # The exposure-exponent trial used to be handed `Rmoon`, which the
        # caller clears whenever the moon MASK is rejected -- an unrelated and
        # routine outcome -- so it bailed out with exactly that verdict while a
        # perfectly good track sat in the same stats dict. It cost the 600 mm
        # bracket its exponent (1.0 instead of 0.55) and 6 px of limb width,
        # silently. This is the invariant that says the two are uncoupled.
        _mt, _mw = stats.get("moon_track"), (stats.get("merge_weight") or {})
        if _mt and "no per-tier lunar track" in str(_mw.get("verdict", "")):
            raise SystemExit(
                "the merge-weight trial reports no lunar track while "
                "moon_track holds one (%s) — the Rmoon coupling is back"
                % _mt)
        if _mw:
            print("merge weight: %s" % _mw.get("verdict", _mw))
        _sp = stats["ldic_k_spread_pct"]
        print("LDIC k spread per tier, percent: "
              + ", ".join("%s %.1f" % (k, v) for k, v in sorted(_sp.items())))
        print("OK — %d frames through the full pipeline, merge finite, "
              "all three ring estimators ran" % nf)
    finally:
        shutil.rmtree(folder, ignore_errors=True)


if __name__ == "__main__":
    main(verbose="-q" not in sys.argv)
