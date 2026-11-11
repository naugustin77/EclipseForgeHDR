"""The corona as sound (Extras, 0.30.0).

A radial line sweeps once around the Sun. At every moment the brightness
along that line, relative to the median corona at the same distance, sets the
loudness of one partial per distance band: near the limb low notes, far out
high notes. A streamer sounds, a gap goes quiet, a prominence rings a bell
(its tone from its shape and size: thin and tall high, broad and low deep).
Nothing but the picture's own numbers: the brightness is read from the
finished picture as tab 3 shows it; the prominences are found on the stack
itself, so they are the same in every look.

Writes a WAV and, when ffmpeg is there, an MP4 with the picture, the
sweeping line and a bar per partial.
"""
import math
import os

import numpy as np
from PIL import Image, ImageDraw
from scipy import ndimage
from scipy.io import wavfile

from . import extras as _ex
from . import film as _film

OPTIONS = [
    dict(id="turn", label="One turn", type="select", default="24",
         values=[["16", "16 s"], ["24", "24 s"], ["40", "40 s"], ["60", "60 s"]],
         title="How long the line takes to go once around the Sun; the clip is one turn."),
    dict(id="pitch", label="Pitch range", type="select", default="mid",
         values=[["low", "low (60–1200 Hz)"], ["mid", "middle (90–1800 Hz)"], ["high", "high (150–3000 Hz)"]],
         title="The partials from the limb (lowest) to the outer corona (highest)."),
    dict(id="law", label="Loudness", type="select", default="normal",
         values=[["gentle", "gentle"], ["normal", "normal"], ["strong", "strong"]],
         title="Always louder where the picture is brighter; gentle follows the brightness evenly, strong lets the streamers ring out and the gaps fall silent."),
    dict(id="outer", label="Out to", type="select", default="3",
         values=[["2", "2 radii"], ["3", "3 radii"], ["4", "4 radii"]],
         title="The outermost distance that sounds (where the picture still shows corona)."),
    dict(id="bells", label="Prominence bells", type="check", default=True,
         title="A bell when the line passes a prominence, its tone from the prominence's shape and size: a tall, "
               "thin one rings high and bright, a broad, low one deep and long; the bigger, the louder."),
    dict(id="aspect", label="Format", type="select", default="16:9",
         values=_ex.SOCIAL_VALUES),
    _ex.SAFE_OPT,
    _ex.rotate_opt("Turns the picture in the clip around the Sun. The sound stays the same: the line always starts "
                   "at the angle origin and goes through east. North and zenith need the star match (or north by hand)."),
    dict(id="video", label="Video", type="check", default=True,
         title="Also write the clip (MP4) with the sweeping line; off: only the WAV."),
]
PITCH = {"low": (60.0, 1200.0), "mid": (90.0, 1800.0), "high": (150.0, 3000.0)}
LAW = {"gentle": 1.0, "normal": 1.8, "strong": 3.0}
NB = 48
NA = 720          # half a degree per cell
KF_MIN = 0.6      # the smallest lettering in a safe zone, of its size without one


def _opts(opts):
    v = {}
    for o in OPTIONS:
        x = (opts or {}).get(o["id"], o["default"])
        v[o["id"]] = (x if isinstance(x, bool) else str(x).lower() in ("1", "true", "on", "yes")) \
            if o["type"] == "check" else str(x if x is not None else o["default"])
    if v["rotate"] not in [r[0] for r in _ex.ROTATE_VALUES]:
        v["rotate"] = "none"
    return v


def _rot_of(ctx, v):
    th = _ex.turn_angle(ctx, v.get("rotate", "none"))
    return 0.0 if th is None else float(th) % 360.0


# THE BELLS. Each prominence its own bell, from its shape and its size as
# found on the stack: the slenderness (height over its foot along the limb)
# and the size (the geometric mean of the two). Tall and thin rings high,
# bright and short; broad and low rings deep, with the hum of a large bell,
# a softer strike, a slow beat and a long tail; the bigger, the louder and
# the longer. The partials are a bell's (hum, prime, tierce, quint, nominal
# and above), their weights following the shape.
BELL_RATIOS = (0.5, 1.0, 1.2, 1.5, 2.0, 2.5, 2.667, 3.0, 4.0)


def _bell_of(p, ctx):
    """dict(f, slender, size_km, tau, loud) for a prominence; None when too small."""
    hk = p["height_km"] if p["height_km"] is not None else p["height"] * _ex.KM_PER_RSUN
    if hk < 5000:
        return None
    rm = ctx["Rm"] / ctx["Rs"]
    base = max(2000.0, math.radians(float(p.get("width_deg") or 1.0)) * rm * _ex.KM_PER_RSUN)
    slender = hk / base
    size = math.sqrt(hk * base)
    lf = 9.3 + 0.55 * math.log2(slender) - 0.3 * math.log2(size / 15000.0)
    f = float(np.clip(2.0 ** lf, 130.0, 1600.0))
    shape = float(np.clip((math.log2(slender) + 3.0) / 4.5, 0.0, 1.0))      # 0 broad and low, 1 tall and thin
    tau = float(np.clip(0.35 + 1.1 * math.sqrt(size / 50000.0) * (1.25 - 0.5 * shape), 0.3, 2.6))
    loud = 0.3 + 0.7 * float(np.clip(size / 60000.0, 0.0, 1.0)) ** 0.7
    return dict(f=f, slender=slender, size_km=size, shape=shape, tau=tau, loud=loud, hk=hk)


def _ring(b, t, t0, T):
    """The bell `b` struck at t0, ringing on across the loop."""
    s = b["shape"]
    dt = (t - t0) % T
    amps = (0.1 + 0.9 * (1 - s), 1.0, 0.7, 0.35, 0.5 + 0.5 * s, 0.15 + 0.45 * s,
            0.1 + 0.35 * s, 0.1 + 0.35 * s, 0.05 + 0.3 * s)
    stretch = 1.0 + 0.012 * s                       # the thin bells a little glassy
    t_att = 0.002 + 0.02 * (1 - s)                   # the broad ones a softer strike
    att = 1.0 - np.exp(-dt / t_att)
    out = np.zeros_like(t)
    for k, (r, a) in enumerate(zip(BELL_RATIOS, amps)):
        fk = b["f"] * r * (stretch ** (r - 1.0))
        if fk > 16000.0:
            continue
        tk = b["tau"] * (1.6 if r < 1 else r ** -0.7)
        e = np.exp(-dt / tk)
        tone = np.sin(2 * np.pi * fk * t + k)
        if s < 0.6 and r <= 1.0:                     # the slow beat of a large bell
            tone = 0.5 * tone + 0.5 * np.sin(2 * np.pi * fk * (1.0 + 0.004 * (1 - s)) * t + k + 1.0)
        out += a * e * tone
    return out * att / sum(amps)


def _analyse(rgb, ctx, outer):
    """Loudness per partial and angle (NB x NA), the band edges and the
    prominences, from the picture. Every cell is the mean of the picture
    over its whole wedge (8 rays per half degree, a sample every 1.5 px out),
    so no streamer between two rays is missed."""
    src = np.asarray(rgb, np.float32)
    r0 = ctx["Rm"] / ctx["Rs"] + 0.02
    H, W = src.shape[:2]
    rmax = min(ctx["cx"], W - 1 - ctx["cx"], ctx["cy"], H - 1 - ctx["cy"]) / ctx["Rs"]
    r1 = min(float(outer), rmax)
    sub = 8
    th = (np.arange(NA * sub) + 0.5) * 360.0 / (NA * sub)        # the sweep runs with the angle
    nr = max(20, int(round((r1 - r0) * ctx["Rs"] / 1.5)))
    rr = r0 + (np.arange(nr) + 0.5) / nr * (r1 - r0)
    phi = np.radians(ctx["o"]["zero_cw"] - ctx["o"]["sgn"] * th)[None, :]
    lum = 0.2126 * src[..., 0] + 0.7152 * src[..., 1] + 0.0722 * src[..., 2]
    edges = np.geomspace(r0, r1, NB + 1)
    band = np.clip(np.searchsorted(edges, rr, side="right") - 1, 0, NB - 1)
    acc = np.zeros((NB, NA * sub), np.float64)
    cnt = np.zeros(NB)
    for k in range(0, nr, 64):
        rk = rr[k:k + 64, None]
        Y = ctx["cy"] - np.cos(phi) * rk * ctx["Rs"]
        X = ctx["cx"] + np.sin(phi) * rk * ctx["Rs"]
        Pk = ndimage.map_coordinates(lum, [Y, X], order=1, mode="nearest")
        Pk = Pk / np.maximum(np.median(Pk, axis=1, keepdims=True), 1e-6)   # vs the median corona there
        np.add.at(acc, band[k:k + 64], Pk)
        np.add.at(cnt, band[k:k + 64], 1)
    amp = (acc / np.maximum(cnt, 1)[:, None]).reshape(NB, NA, sub).mean(2).astype(np.float32)
    strip, th_s, rr_s = _ex._limb_strip(src, ctx, 1.35)
    proms = _ex.find_prominences(src, ctx, strip, th_s, rr_s)
    return amp, edges, proms, r1


def _gain(amp, law):
    """Loudness 0.02..1 from the brightness ratio, monotonic over its whole
    range: a brighter part is always louder, a darker gap always quieter."""
    lo, hi = np.percentile(amp, 0.5), np.percentile(amp, 99.5)
    x = np.clip((amp - lo) / max(hi - lo, 1e-6), 0.0, 1.0)
    return (0.02 + 0.98 * x ** law).astype(np.float32)


def _at(row, pa):
    """A row of the map (NA cells round the Sun) at angle(s) pa."""
    return np.interp(pa, (np.arange(NA) + 0.5) * 360.0 / NA, row, period=360)


def make(rgb, geo, wd, folder, opts, out_base, prog=None, scale=1.0):
    """Writes <out_base>.wav (and .mp4). Returns (paths, info)."""
    v = _opts(opts)
    log = (lambda m, f=None: prog.log(m, f) if prog else None)
    ctx = _ex.context(wd, folder, geo, scale=scale, **_ex.angle_kw(opts))
    T = float(v["turn"])
    f_lo, f_hi = PITCH.get(v["pitch"], PITCH["mid"])
    law = LAW.get(v["law"], 1.8)
    log("sound: reading the corona around the Sun ...", 0.1)
    amp, edges, proms, r1 = _analyse(rgb, ctx, v["outer"])
    gain = _gain(amp, law)
    freqs = np.geomspace(f_lo, f_hi, NB)
    sr = 44100
    n = int(sr * T)
    t = np.arange(n, dtype=np.float64) / sr
    th_t = (t / T) * 360.0
    log("sound: the partials ...", 0.3)
    audio = np.zeros(n)
    sm = max(64, int(sr * T / 360.0 * 0.5))       # half a degree: no clicks, no blur
    for i in range(NB):
        gi = _at(gain[i], th_t)
        env = ndimage.uniform_filter1d(gi, sm, mode="wrap")
        audio += env * np.sin(2 * np.pi * freqs[i] * t + i) / NB
    bells = []
    if v["bells"]:
        for p in proms:
            b = _bell_of(p, ctx)
            if b is None:
                continue
            audio += _ring(b, t, p["pa"] / 360.0 * T, T) * b["loud"] * 0.9
            bells.append(dict(pa=p["pa"], height_km=round(b["hk"]), width_deg=p.get("width_deg"),
                              f=round(b["f"]), slender=round(b["slender"], 2)))
    audio /= max(1e-9, np.abs(audio).max()) * 1.05
    wav = out_base + ".wav"
    k = 0
    while os.path.exists(wav):
        k += 1
        wav = out_base + f"_{k}.wav"
    wavfile.write(wav, sr, (audio * 32767).astype(np.int16))
    paths = [wav]
    info = (f"{T:.0f} s, {NB} partials {f_lo:.0f}–{f_hi:.0f} Hz, out to {r1:.1f} {_ex._unit(ctx)}"
            + (f", {len(bells)} bells {min(b['f'] for b in bells):.0f}–{max(b['f'] for b in bells):.0f} Hz"
               if len(bells) > 1 else (f", 1 bell {bells[0]['f']:.0f} Hz" if bells else "")))
    if v["video"]:
        log("sound: the clip ...", 0.5)
        mp4 = os.path.splitext(wav)[0] + ".mp4"
        try:
            _video(rgb, ctx, gain, edges, mp4, wav, T, v["aspect"], f_lo, f_hi, r1, log, v["safe"], _rot_of(ctx, v))
            paths.append(mp4)
        except Exception as e:
            log(f"sound: no clip ({e}); the WAV is there")
    return paths, info


def still(rgb, geo, wd, folder, opts, scale=1.0):
    """One frame of the clip, for the card."""
    v = _opts(opts)
    ctx = _ex.context(wd, folder, geo, scale=scale, **_ex.angle_kw(opts))
    f_lo, f_hi = PITCH.get(v["pitch"], PITCH["mid"])
    amp, edges, _p, r1 = _analyse(rgb, ctx, v["outer"])
    gain = _gain(amp, LAW.get(v["law"], 1.8))
    draw = _layout(rgb, ctx, gain, edges, float(v["turn"]), v["aspect"], f_lo, f_hi, r1, v["safe"], _rot_of(ctx, v))
    return draw(118.0)


def _turned_square(src, ctx, th, half, side):
    """The square of half-width `half` (picture px) around the Sun with the
    picture direction `th` (deg clockwise from up) standing up, at `side` px."""
    n = int(min(max(side, 2 * half), 4096))
    g = (np.arange(n, dtype=np.float32) + 0.5) / n * 2 * half - half
    v, u = np.meshgrid(g, g, indexing="ij")
    t = math.radians(th)
    ca, sa = math.cos(t), math.sin(t)
    X = ctx["cx"] + u * ca - v * sa
    Y = ctx["cy"] + u * sa + v * ca
    a = np.stack([ndimage.map_coordinates(src[..., c], [Y, X], order=1, mode="constant", cval=0.0)
                  for c in range(3)], -1)
    im = Image.fromarray(_ex._to8(a))
    return im if n == side else im.resize((side, side), Image.LANCZOS)


def _layout(rgb, ctx, gain, edges, T, aspect, f_lo, f_hi, r1, safe="off", rot=0.0):
    """Returns draw(pa) -> PIL frame of the clip at position angle `pa`."""
    W, H = _film.SIZES.get(aspect, _film.SIZES["16:9"])
    src = np.asarray(rgb, np.float32)
    s_top, s_right, s_bot, s_left, s_cw, s_ch = _ex.safe_margins(W, H, safe)
    rx = None                                    # the right edge of the lettering (upright layouts)
    hb = H                                       # the foot of the lettering
    kf = 1.0                                     # the lettering's size (smaller in a tight safe zone)
    # layout: the picture square, the bars beside (16:9) or below (the rest)
    if aspect == "16:9":
        side = H
        px, py = (W - side) // 2, 0
        bx0, bx1, by0, by1 = px + side + 40, W - 40, 80, H - 60
        lx, ly = 40, 40
    elif aspect in ("1:1", "9:16") and not s_top:
        # the two layouts as they were first made
        fs0 = max(14, int(min(W, H) / 40))
        if aspect == "1:1":
            side = int(H * 0.66)
            px, py = (W - side) // 2, 0
            bx0, bx1, by0, by1 = 60, W - 60, side + int(fs0 * 4.2), H - int(fs0 * 5.4)
            lx, ly = 60, side + int(fs0 * 0.6)
        else:
            side = W
            px, py = 0, int(H * 0.12)
            bx0, bx1, by0, by1 = 80, W - 80, py + side + int(fs0 * 4.2), H - int(fs0 * 5.0)
            lx, ly = 80, py + side + int(fs0 * 0.6)
    else:
        # upright, any shape, and in a safe zone. The picture keeps the size
        # and place it has without a zone (it may lie under the app's bars:
        # only the lettering has to stay clear of them); it moves up only
        # when the lettering would not fit below it, and the lettering is
        # set smaller to fit the room left inside the zone. Only when even
        # the smallest lettering does not fit (a reel's zone on a square)
        # does the picture give way, just enough.
        fs0 = max(14, int(min(W, H) / 40))
        need = fs0 * 15.2                         # 0.6 + 3.6 above the bars, 6 for them, 5 below
        if aspect == "1:1":
            side, py = int(H * 0.66), 0
        elif aspect == "9:16":
            side, py = W, int(H * 0.12)
        else:
            py = int(H * 0.04)
            side = max(int(min(W, H - py - need)), int(W * 0.4))
        hb = H - s_bot
        if py + side + need > hb:
            py = max(0, int(hb - need - side))
        kf = min(1.0, (hb - max(py + side, s_top)) / need)
        if kf < KF_MIN:
            # the smallest lettering; the bars give up height first (down to
            # 2.4 font sizes from 6), the picture only after that
            kf = KF_MIN
            need_min = fs0 * kf * (0.6 + 3.6 + 2.4 + 5.0)
            if hb - max(py + side, s_top) < need_min:
                side = max(int(W * 0.4), int(hb - need_min - py))
        px = (W - side) // 2
        fsk = fs0 * kf
        ly = max(py + side, s_top) + int(fsk * 0.6)
        lx = max(60, s_left)
        rx = W - max(60, s_right)
        if s_cw and ly - fsk > H - s_ch:                 # the buttons' column beside the lettering
            rx = min(rx, W - s_cw)
        bx0, bx1, by0, by1 = lx, rx, ly + int(fsk * 3.6), hb - int(fsk * 5.0)
    half = (r1 + 0.1) * ctx["Rs"]
    rot = float(rot or 0.0) % 360.0
    Hs, Ws = src.shape[:2]
    if rot or ctx["cx"] - half < 0 or ctx["cy"] - half < 0 or ctx["cx"] + half > Ws or ctx["cy"] + half > Hs:
        # turned, or the square reaches past the picture's edge (a frame
        # tight around the corona): sampled, black beyond the edge
        pic = _turned_square(src, ctx, rot, half, side)
        if rot:
            o = dict(ctx["o"])
            o["zero_cw"] = (o["zero_cw"] - rot) % 360.0      # the line follows the turned picture
            ctx = dict(ctx, o=o)
    else:
        pic = Image.fromarray(_ex._to8(src)).resize(
            (side, side), Image.LANCZOS, box=(ctx["cx"] - half, ctx["cy"] - half, ctx["cx"] + half, ctx["cy"] + half))
    sc = side / (2 * half)
    cxp, cyp = px + side / 2.0, py + side / 2.0
    fs = max(12, int(round(min(W, H) / 40 * kf)))
    f_s = _ex._font(fs * 0.8)
    f_cap = _ex._font(fs * 0.56, True)                       # labels: letter-spaced capitals
    f_it = _ex._serif(fs * 0.9)                              # sentences: serif italic
    unit = _ex._unit(ctx)
    horizontal = aspect == "16:9"
    org = ctx["o"].get("origin")
    ref = {"solar": "from solar north", "celestial": "from celestial north"}.get(org, "from the picture\u2019s top")
    lo_lab, hi_lab = f"{edges[0]:.1f} {unit} \u00b7 {f_lo:.0f} Hz", f"{r1:.1f} {unit} \u00b7 {f_hi / 1000:.1f} kHz"
    key = [("one turn", f"{T:.0f} seconds"), ("low notes", "near the limb"), ("high notes", "far out"),
           ("loud", "a streamer"), ("quiet", "a gap"), ("high bell", "a thin prominence"), ("deep bell", "a broad one")]
    # everything that does not move, drawn once (the sheet's voice: a kicker
    # in the accent colour, labels in letter-spaced capitals, sentences in
    # italic); the frame adds the line, the notes, the bars and the angle
    static = Image.new("RGBA", (W, H), (0, 0, 0, 0))
    ds = _ex.Draw(static)
    if horizontal:
        ds.rectangle([0, 0, px, H], fill=_ex.BG + (255,))
        ds.tracked((lx, ly), "THE CORONA AS SOUND", fill=_ex.GOLD + (255,), font=_ex._font(fs * 0.62, True),
                   anchor="la", track=0.2)
        ang_xy, ang_f, ang_an = (lx, ly + fs * 1.5), _ex._font(fs * 3.4), "la"
        ds.tracked((lx, ly + fs * 6.2), ref.upper(), fill=_ex.MUTED + (255,), font=f_cap, anchor="la", track=0.16)
        if org in ("solar", "celestial"):
            ds.tracked((lx, ly + fs * 7.2), "THROUGH EAST", fill=_ex.MUTED + (255,), font=f_cap, anchor="la", track=0.16)
        yk = ly + fs * 10.2
        ds.line([(lx, yk - fs * 0.9), (px - lx, yk - fs * 0.9)], fill=(48, 48, 48, 255), width=max(1, fs // 14))
        for i, (lab, val) in enumerate(key):
            ds.tracked((lx, yk + i * fs * 1.55 + fs * 0.28), lab.upper(), fill=_ex.MUTED + (255,), font=f_cap,
                       anchor="la", track=0.16)
            ds.text((lx + fs * 5.6, yk + i * fs * 1.55), val, fill=_ex.TEXT + (255,), font=f_it, anchor="la")
        ds.text((lx, H - fs * 4.6), "From the picture\u2019s brightness,\nagainst the median corona\nat the same distance.",
                fill=(140, 140, 140, 255), font=_ex._serif(fs * 0.8))
        ds.tracked((bx0, by0 - fs * 2.6), "PLAYING NOW", fill=_ex.MUTED + (255,), font=f_cap, anchor="la", track=0.16)
        ds.text((bx0, by1 + 8), lo_lab, fill=(135, 135, 135, 255), font=f_s)
        ds.text((bx0, by0 - fs * 1.3), hi_lab, fill=(135, 135, 135, 255), font=f_s)
    else:
        rx_ = rx if rx is not None else W - lx
        ds.tracked((lx, ly + fs * 0.35), "THE CORONA AS SOUND", fill=_ex.GOLD + (255,), font=_ex._font(fs * 0.66, True),
                   anchor="la", track=0.2)
        ang_xy, ang_f, ang_an = (rx_, ly - fs * 0.1), _ex._font(fs * 1.9), "ra"
        ds.tracked((rx_, ly + fs * 2.3), ref.upper(), fill=_ex.MUTED + (255,), font=f_cap, anchor="ra", track=0.16)
        ds.tracked((bx0, by0 - fs * 1.5), "PLAYING NOW", fill=_ex.MUTED + (255,), font=f_cap, anchor="la", track=0.16)
        ds.text((bx0, by1 + 8), lo_lab, fill=(135, 135, 135, 255), font=f_s)
        ds.text((bx1, by1 + 8), hi_lab, fill=(135, 135, 135, 255), font=f_s, anchor="ra")
        # the two sentences, smaller where the lettering is narrower
        l1 = f"One turn in {T:.0f} seconds. Low notes near the limb, high notes far out;"
        l2 = "loud is a streamer, quiet a gap; a prominence rings, thin ones high, broad ones deep."
        f_it2 = f_it
        wmax = rx_ - lx
        wt = max(ds.textlength(l1, font=f_it), ds.textlength(l2, font=f_it))
        if wt > wmax:
            f_it2 = _ex._serif(fs * 0.9 * wmax / wt)
        ds.text((lx, hb - fs * 2.6), l1, fill=_ex.TEXT + (255,), font=f_it2)
        ds.text((lx, hb - fs * 1.45), l2, fill=_ex.TEXT + (255,), font=f_it2)

    def draw(pa):
        fr = Image.new("RGB", (W, H), _ex.BG)
        fr.paste(pic, (px, py))
        d = _ex.Draw(fr, "RGBA")
        phi = math.radians(ctx["o"]["zero_cw"] - ctx["o"]["sgn"] * pa)
        ux, uy = math.sin(phi), -math.cos(phi)
        r_in, r_out = ctx["Rm"] / ctx["Rs"] * ctx["Rs"] * sc, r1 * ctx["Rs"] * sc
        d.line([(cxp + ux * r_in, cyp + uy * r_in), (cxp + ux * r_out, cyp + uy * r_out)],
               fill=(255, 220, 150, 230), width=max(2, fs // 7))
        for i in range(NB):
            rmid = math.sqrt(edges[i] * edges[i + 1]) * ctx["Rs"] * sc
            gi = float(_at(gain[i], pa))
            rad = 1.5 + 0.35 * fs * gi
            d.ellipse([cxp + ux * rmid - rad, cyp + uy * rmid - rad, cxp + ux * rmid + rad, cyp + uy * rmid + rad],
                      fill=(255, 200, 100, 140))
        fr.paste(static, (0, 0), static)
        d = _ex.Draw(fr, "RGBA")
        # the bars: one per partial, low notes at the bottom (16:9) or left
        for i in range(NB):
            gi = float(_at(gain[i], pa))
            a = int(70 + 180 * gi)
            if horizontal:
                yb = by1 - (i + 0.5) * (by1 - by0) / NB
                wbar = 6 + (bx1 - bx0 - 6) * gi
                d.rectangle([bx0, yb - 3, bx0 + wbar, yb + 3], fill=(232, 163, 61, a))
            else:
                xb = bx0 + (i + 0.5) * (bx1 - bx0) / NB
                hbar = 6 + (by1 - by0 - 6) * gi
                d.rectangle([xb - 4, by1 - hbar, xb + 4, by1], fill=(232, 163, 61, a))
        # the one figure that moves: the angle
        d.text(ang_xy, f"{pa:.0f}\u00b0", fill=_ex.INK + (255,), font=ang_f, anchor=ang_an)
        return fr
    return draw


def _video(rgb, ctx, gain, edges, mp4, wav, T, aspect, f_lo, f_hi, r1, log, safe="off", rot=0.0):
    W, H = _film.SIZES.get(aspect, _film.SIZES["16:9"])
    fps = 25
    NF = int(T * fps)
    draw = _layout(rgb, ctx, gain, edges, T, aspect, f_lo, f_hi, r1, safe, rot)
    wr = _film._Writer(mp4, W, H, fps)
    if wr.kind != "mp4":
        raise RuntimeError("no ffmpeg")
    for kf in range(NF):
        wr.add(draw(kf / NF * 360.0))
    tmp = wr.close()
    # mux the sound in
    import shutil
    import subprocess
    exe = wr.exe
    muxed = os.path.splitext(mp4)[0] + "_a.mp4"
    subprocess.run([exe, "-y", "-loglevel", "error", "-i", tmp, "-i", wav, "-c:v", "copy", "-c:a", "aac",
                    "-b:a", "160k", "-shortest", muxed], check=True)
    shutil.move(muxed, mp4)
