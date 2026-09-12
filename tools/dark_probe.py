#!/usr/bin/env python3
"""Measure what a set of bias and dark frames actually says about a sensor.

STANDALONE ON PURPOSE. Imports only numpy and rawpy, not the eclipseforgehdr
package, so it can be run with whatever interpreter already has rawpy -- on a
pipx install that is:

    ~/.local/pipx/venvs/eclipseforgehdr/bin/python tools/dark_probe.py \
        --bias ~/Pictures/eclipse/bias --dark ~/Pictures/eclipse/dark

EFHDR has no dark-frame support and this script does not add any. It is a
measurement, and it exists to answer questions that cannot be answered from
the eclipse frames themselves, because in those the thing being measured is
inseparable from real corona.

WHAT THE BIAS FRAMES ANSWER

  1. Does the black level the raw file REPORTS match the black that is
     actually there? `raw.py` subtracts `black_level_per_channel` and trusts
     it. If the reported value is off by a few ADU per channel, that error
     lands on the outer corona -- which sits only tens of ADU above black --
     and it is then multiplied by the white-balance gain, so it arrives as a
     colour cast exactly where the signal is weakest.

  2. IS THE DISTRIBUTION CLIPPED AT ZERO? This is the one that matters most
     and the one that cannot be guessed. A bias frame is noise about a mean;
     its histogram must be symmetric. If the camera clamps at the black level
     before writing the file, the left half is missing, the measured mean sits
     ABOVE the true black, and every pedestal estimate downstream is fitting
     something that has already been truncated. The test is the histogram
     shape, and the giveaway is a spike in the lowest occupied bin plus a
     median that sits above the mode.

  3. Read noise per channel, in ADU. Needed to judge whether anything else
     measured here is above the noise at all.

WHAT THE DARK FRAMES ANSWER

  4. Dark current in ADU/s per channel, at the temperature they were shot at.
     Dark current roughly doubles every 5-7 C, so this number is only valid
     near the ambient of the shoot -- which is why they are worth taking at
     30 C rather than indoors.

  5. Hot pixels: how many, and how bright at the longest tier's exposure.

  6. FIXED SPATIAL STRUCTURE. Anything that sits in the same place in every
     frame and grows with exposure time lands where the outer-field gradients
     live, and stacking cannot average it out. Amp glow is the classic case
     and a tester reports the S1R II has none, which is a prediction this makes
     checkable rather than a reason not to look: the script reports
     centre-versus-corner and the worst edge, and writes a downsampled map, so
     "there is nothing there" becomes a number instead of an assumption.

NOTHING HERE IS A CORRECTION. It prints numbers.
"""
import argparse
import glob
import os
import sys

import numpy as np

_RAW_EXT = (".rw2", ".raw", ".cr2", ".cr3", ".nef", ".arw", ".dng", ".orf",
            ".raf", ".pef", ".srw", ".rwl")


def _files(folder):
    out = []
    for p in sorted(glob.glob(os.path.join(folder, "*"))):
        if os.path.splitext(p)[1].lower() in _RAW_EXT:
            out.append(p)
    return out


def _open(path):
    """Raw mosaic BEFORE black subtraction, plus what the file claims."""
    import rawpy
    with rawpy.imread(path) as raw:
        bayer = raw.raw_image_visible.astype(np.float32)
        pat = raw.raw_pattern
        meta = {"black": np.asarray(raw.black_level_per_channel, float)[:4],
                "white": float(raw.white_level),
                "pattern": None if pat is None else np.asarray(pat).copy()}
        return bayer.copy(), meta


def _plane_index(pattern):
    """Map each of the four mosaic positions to a colour name."""
    if pattern is None or np.asarray(pattern).shape != (2, 2):
        return {(0, 0): "c0", (0, 1): "c1", (1, 0): "c2", (1, 1): "c3"}
    names = {0: "R", 1: "G", 2: "B", 3: "G2"}
    out = {}
    for y in range(2):
        for x in range(2):
            out[(y, x)] = names.get(int(pattern[y, x]), "c?")
    return out


def _stack(paths, label, limit=None):
    """Per-pixel mean, plus per-frame per-channel scalars. Streams."""
    if limit:
        paths = paths[:limit]
    acc = None
    meta0 = None
    per_frame = []
    hist = np.zeros(1 << 16, np.int64)
    for i, p in enumerate(paths):
        bayer, meta = _open(p)
        if acc is None:
            acc = np.zeros(bayer.shape, np.float64)
            meta0 = meta
        elif bayer.shape != acc.shape:
            print("  skipping %s: %s, expected %s"
                  % (os.path.basename(p), bayer.shape, acc.shape))
            continue
        acc += bayer
        idx = _plane_index(meta["pattern"])
        row = {}
        for (y, x), nm in idx.items():
            sub = bayer[y::2, x::2]
            row[nm] = (float(sub.mean()), float(sub.std()))
        per_frame.append(row)
        v = np.clip(bayer, 0, 65535).astype(np.uint16).ravel()
        hist += np.bincount(v, minlength=1 << 16)
        sys.stdout.write("\r  %s: %d/%d" % (label, i + 1, len(paths)))
        sys.stdout.flush()
    print()
    if acc is None:
        return None
    return {"master": (acc / len(per_frame)).astype(np.float32),
            "meta": meta0, "per_frame": per_frame, "hist": hist,
            "n": len(per_frame)}


def _chan_stats(per_frame):
    """Mean of the per-frame means, and read noise from the per-frame stds."""
    keys = list(per_frame[0].keys())
    out = {}
    for k in keys:
        mu = np.array([f[k][0] for f in per_frame])
        sd = np.array([f[k][1] for f in per_frame])
        out[k] = {"mean": float(mu.mean()), "mean_spread": float(mu.std()),
                  "noise": float(np.median(sd))}
    return out


def _clip_report(hist, black_mean):
    """Is the left tail of the bias histogram intact, or cut off?"""
    occ = np.nonzero(hist)[0]
    if occ.size == 0:
        return None
    lo, hi = int(occ[0]), int(occ[-1])
    vals = np.arange(lo, hi + 1)
    cnt = hist[lo:hi + 1].astype(float)
    tot = cnt.sum()
    mode = int(vals[int(np.argmax(cnt))])
    csum = np.cumsum(cnt)
    med = float(vals[int(np.searchsorted(csum, 0.5 * tot))])
    mean = float((vals * cnt).sum() / tot)
    # How much of the distribution sits in the single lowest occupied value?
    # In an unclipped Gaussian that bin holds a vanishing fraction; in a
    # clipped one it holds the whole missing tail.
    spike = float(cnt[0] / tot)
    # Symmetry: counts an equal distance either side of the mode.
    d = max(1, int(round(0.5 * (mode - lo))))
    a = float(cnt[max(0, (mode - lo) - d):(mode - lo)].sum())
    b = float(cnt[(mode - lo) + 1:(mode - lo) + 1 + d].sum())
    return {"lowest": lo, "highest": hi, "mode": mode, "median": med,
            "mean": mean, "spike_frac": spike, "left": a, "right": b,
            "asym": (b - a) / max(a + b, 1.0), "d": d}


def _spatial(master, label, out_png=None, ny=24, nx=32):
    """Block-average the master frame down to a small map."""
    H, W = master.shape
    by, bx = H // ny, W // nx
    m = master[:ny * by, :nx * bx].reshape(ny, by, nx, bx).mean(axis=(1, 3))
    m = m - float(np.median(m))
    cy0, cx0 = ny // 4, nx // 4
    centre = float(np.median(m[cy0:ny - cy0, cx0:nx - cx0]))
    corners = {"TL": float(np.median(m[:cy0, :cx0])),
               "TR": float(np.median(m[:cy0, -cx0:])),
               "BL": float(np.median(m[-cy0:, :cx0])),
               "BR": float(np.median(m[-cy0:, -cx0:]))}
    edges = {"top": float(np.median(m[:2, :])),
             "bottom": float(np.median(m[-2:, :])),
             "left": float(np.median(m[:, :2])),
             "right": float(np.median(m[:, -2:]))}
    if out_png:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            v = float(np.nanmax(np.abs(m))) or 1.0
            plt.figure(figsize=(6, 4.5))
            plt.imshow(m, cmap="coolwarm", vmin=-v, vmax=v)
            plt.colorbar(label="ADU above the frame median")
            plt.title("%s — fixed spatial structure" % label)
            plt.tight_layout()
            plt.savefig(out_png, dpi=110)
            plt.close()
        except Exception as e:
            print("  (no PNG: %s)" % e)
            out_png = None
    return {"map": m, "centre": centre, "corners": corners, "edges": edges,
            "png": out_png}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--bias", required=True, help="folder of bias frames")
    ap.add_argument("--dark", default=None, help="folder of dark frames")
    ap.add_argument("--dark-seconds", type=float, default=None,
                    help="exposure of the dark frames, for the ADU/s rate")
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--outdir", default=".")
    a = ap.parse_args()

    try:
        import rawpy                                        # noqa: F401
    except ImportError:
        print("rawpy is not available to this interpreter. On a pipx install:\n"
              "  ~/.local/pipx/venvs/eclipseforgehdr/bin/python %s ..."
              % sys.argv[0])
        return 2

    bf = _files(a.bias)
    if not bf:
        print("no raw files in %s" % a.bias)
        return 2
    print("BIAS  %d frame(s) in %s" % (len(bf), a.bias))
    B = _stack(bf, "bias", a.limit)
    bs = _chan_stats(B["per_frame"])
    blk = B["meta"]["black"]
    idx = _plane_index(B["meta"]["pattern"])
    order = [idx[(0, 0)], idx[(0, 1)], idx[(1, 0)], idx[(1, 1)]]

    print()
    print("1. REPORTED BLACK vs MEASURED BLACK")
    print("   the file says black_level_per_channel = %s, white_level = %g"
          % (np.array2string(blk, precision=1), B["meta"]["white"]))
    print("   channel   measured    reported    difference   read noise")
    pat = B["meta"]["pattern"]
    for (y, x), nm in idx.items():
        rep = float(blk[int(pat[y, x])]) if pat is not None else float(blk.mean())
        me = bs[nm]["mean"]
        print("   %-8s %9.2f %11.2f %13.2f %11.2f"
              % (nm, me, rep, me - rep, bs[nm]["noise"]))
    dmax = max(abs(bs[nm]["mean"]
                   - (float(blk[int(pat[y, x])]) if pat is not None
                      else float(blk.mean())))
               for (y, x), nm in idx.items())
    nmin = min(bs[nm]["noise"] for nm in order)
    print("   worst channel is off by %.2f ADU, against a read noise of %.2f "
          "ADU (%.2f sigma of a single pixel, but the master of %d frames "
          "measures the mean to about %.3f ADU)"
          % (dmax, nmin, dmax / max(nmin, 1e-6), B["n"],
             nmin / max(np.sqrt(B["n"] * B["master"].size / 4.0), 1)))

    print()
    print("2. IS THE BIAS CLIPPED?")
    cr = _clip_report(B["hist"], float(blk.mean()))
    if cr is None:
        print("   no data")
    else:
        print("   occupied range %d .. %d ADU, mode %d, median %.1f, mean %.2f"
              % (cr["lowest"], cr["highest"], cr["mode"], cr["median"],
                 cr["mean"]))
        print("   lowest occupied value holds %.4f%% of all pixels"
              % (100 * cr["spike_frac"]))
        print("   counts within %d ADU below the mode: %.0f ; above: %.0f "
              "(asymmetry %+.3f)" % (cr["d"], cr["left"], cr["right"],
                                     cr["asym"]))
        verdict = []
        if cr["spike_frac"] > 0.002:
            verdict.append("a spike in the lowest bin")
        if cr["asym"] > 0.10:
            verdict.append("a missing left tail")
        if cr["lowest"] >= float(blk.min()) - 0.5:
            verdict.append("nothing below the reported black level")
        if verdict:
            print("   CLIPPED, on the evidence of: " + ", ".join(verdict))
            print("   Consequence: the measured black sits ABOVE the true "
                  "black, and every pedestal estimate downstream is fitting a "
                  "truncated distribution.")
        else:
            print("   NOT clipped — the tail extends below the reported black "
                  "and the distribution is symmetric about the mode. Pedestal "
                  "estimates downstream are working on intact data.")

    if a.dark:
        df = _files(a.dark)
        if not df:
            print("\nno raw files in %s" % a.dark)
            return 2
        print()
        print("DARK  %d frame(s) in %s" % (len(df), a.dark))
        D = _stack(df, "dark", a.limit)
        ds = _chan_stats(D["per_frame"])
        print()
        print("4. DARK CURRENT")
        t = a.dark_seconds
        print("   channel   dark mean   minus bias" +
              ("      ADU/s" if t else ""))
        for nm in order:
            d = ds[nm]["mean"] - bs[nm]["mean"]
            line = "   %-8s %11.2f %12.3f" % (nm, ds[nm]["mean"], d)
            if t:
                line += " %10.4f" % (d / t)
            print(line)
        if not t:
            print("   (pass --dark-seconds to get a rate)")

        print()
        print("5. HOT PIXELS")
        if D["master"].shape == B["master"].shape:
            res = D["master"] - B["master"]
            sd = float(np.median(np.abs(res - np.median(res)))) * 1.4826
            for k in (10, 20, 50):
                n = int((res > k * max(sd, 1e-6)).sum())
                print("   above %2dx the frame sigma (%.2f ADU): %8d px "
                      "(%.5f%% of the sensor)"
                      % (k, k * sd, n, 100.0 * n / res.size))
            print("   brightest pixel in the dark, above bias: %.1f ADU"
                  % float(res.max()))
        else:
            print("   bias and dark frame sizes differ — skipped")

        print()
        print("6. FIXED SPATIAL STRUCTURE")
        for label, S in (("bias", B), ("dark", D)):
            png = os.path.join(a.outdir, "spatial_%s.png" % label)
            sp = _spatial(S["master"], label, png)
            c = sp["corners"]
            print("   %s: centre %+.2f ADU, corners TL %+.2f TR %+.2f "
                  "BL %+.2f BR %+.2f" % (label, sp["centre"], c["TL"],
                                         c["TR"], c["BL"], c["BR"]))
            e = sp["edges"]
            print("         edges top %+.2f bottom %+.2f left %+.2f right %+.2f"
                  % (e["top"], e["bottom"], e["left"], e["right"]))
            worst = max(abs(v - sp["centre"]) for v in
                        list(c.values()) + list(e.values()))
            print("         worst departure from centre: %.2f ADU%s"
                  % (worst, "   -> %s" % sp["png"] if sp["png"] else ""))
        print()
        print("   For scale: how much does this matter to the corona? Compare "
              "the worst departure above against the outer-field signal in a "
              "1.6 s eclipse frame, which on a 600 mm set runs tens of ADU "
              "above black. A structure of a few ADU is a few percent there "
              "and nothing at all inside 2 R.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
