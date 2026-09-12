#!/usr/bin/env python3
"""Rebuild ONLY the partial-convolution masks in an existing work directory.

WHY THIS EXISTS. Pressing Start re-runs every detail layer -- MGN, RHEF, FNRGF,
NAFE and the two inner-corona passes -- even when the stack itself is cached.
On the 600 mm bracket that is about five minutes, of which the masks are two,
and an A/B of one mask setting pays it twice for no reason. This does the polar
blur and nothing else, then writes the new stats into report.json so the
renderer picks up the scale it should normalise the Amplification slider to.

    python tools/rebuild_hill.py <folder or work directory> [denoise profile]

    ECLIPSEFORGE_HILL_NODENOISE=1 python tools/rebuild_hill.py <...>

The denoise profile must match the one the folder was stacked with, because the
masks are built on that master -- it defaults to "fine", which is the app's own
default. Reload the folder in the app afterwards; the masks are read from disk.

Run it with the SAME interpreter the app is installed under, so the masks are
built by the code that will render them:

    "$HOME/Library/Application Support/pipx/venvs/eclipseforgehdr/bin/python" \
        tools/rebuild_hill.py ~/Pictures/.../RAWs
"""
import json
import os
import sys
import time

# Prefer the installed package; fall back to the source tree beside this file.
try:
    from eclipseforgehdr.detail import build_hill
except ImportError:
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
    from eclipseforgehdr.detail import build_hill


class _P:
    def log(self, m, frac=None):
        print("  |", m[:170])
        sys.stdout.flush()


def main():
    if len(sys.argv) < 2:
        raise SystemExit(__doc__)
    p = os.path.abspath(os.path.expanduser(sys.argv[1]))
    denoise = sys.argv[2] if len(sys.argv) > 2 else "fine"
    # Accept either the folder of raws or the work directory inside it, and
    # decide by what is IN the directory rather than by its name -- a work
    # directory copied or renamed is still a work directory.
    wd = p if os.path.exists(os.path.join(p, "hdr_lum.npy")) \
        else os.path.join(p, ".eclipseforgehdr")
    need = ("hdr_lum.npy", "geometry.json")
    missing = [f for f in need if not os.path.exists(os.path.join(wd, f))]
    if missing:
        raise SystemExit("%s is not a stacked work directory (no %s)"
                         % (wd, ", ".join(missing)))
    base = "raw" if os.environ.get("ECLIPSEFORGE_HILL_NODENOISE") == "1" \
        else "denoised"
    print("rebuilding the partial-convolution masks in %s\n  base: %s   "
          "denoise profile: %s" % (wd, base, denoise))
    t0 = time.time()
    st = build_hill(wd, _P(), denoise=denoise)
    if not st:
        raise SystemExit("build_hill did not produce a mask set — see above")
    rp = os.path.join(wd, "report.json")
    try:
        rep = json.load(open(rp)) if os.path.exists(rp) else {}
    except Exception:
        rep = {}
    rep["hill"] = st
    json.dump(rep, open(rp, "w"), indent=1)
    print("done in %.0fs — base=%s, recipe build %d"
          % (time.time() - t0, st.get("base"), st.get("build", 0)))
    print("  reload the folder in the app to see them (no re-stack needed)")


if __name__ == "__main__":
    main()
