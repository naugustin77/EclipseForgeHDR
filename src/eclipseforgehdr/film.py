"""The film (Extras, 0.30.0): from the single frames to the finished
picture, which then unrolls into its polar map.

STAGES. 1 The camera frames, one per exposure from the shortest to the
longest, each dissolving into the next (the best frame of each tier, as the
run picked it, developed as shot with a soft shoulder, aligned to the stack
by the tier's own shift). 2 A dissolve into the plain merged stack (the merge
with one neutral stretch, render.plain_stack), then into the same stack brought to
the finished picture's brightness and colour ring by ring, so the
last dissolve into the finished picture, held, shows only what the detail
layers add. 4 The picture unrolls: the ring around the Sun opens at the
seam and straightens into the strip, the scale unchanged, the seam feathered
so nothing tears. 5 The flat map is moved into its final place: in 16:9 and
1:1 it grows to fill the width; in 9:16 it is larger than the frame and
slides slowly from one end to the other. 6 A hold on the map.

GEOMETRY. One solar radius is a sixth of the frame's short side, the same
through the dissolves, so nothing jumps between stages. Everything is drawn
from the finished picture itself (the composite at full resolution), the
stack and the frames in the folder; no modelled content.

OUTPUT. MP4 (H.264) through ffmpeg when it is on the system or bundled with
imageio-ffmpeg; else an animated GIF, and the log says so.
"""
import json
import math
import os
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor

import numpy as np
from PIL import Image, ImageDraw
from scipy import ndimage

from . import extras as _ex

SIZES = dict(_ex.SOCIAL_SIZES)
PACE = {"brisk": 0.75, "calm": 1.0, "slow": 1.35}
BGC = np.array(_ex.BG, np.float32) / 255.0

OPTIONS = [
    dict(id="aspect", label="Format", type="select", default="16:9",
         values=_ex.SOCIAL_VALUES,
         title="The frame, at the size the apps show it unscaled. In 9:16 the unrolled map is wider than the "
               "frame and slides across it."),
    _ex.SAFE_OPT,
    _ex.rotate_opt("Turns the whole film around the Sun: the frames, the stack and the finished picture stand "
                   "this way up. North and zenith need the star match (or north by hand). The map at the end "
                   "always has the angle origin in its middle."),
    dict(id="pace", label="Pace", type="select", default="calm",
         values=[["brisk", "brisk"], ["calm", "calm"], ["slow", "slow"]],
         title="How long each stage takes; the holds on the finished picture and the map scale with it."),
    dict(id="frames", label="Single frames", type="check", default=True,
         title="Start with the camera frames, one per exposure, dissolving into each other."),
    dict(id="stack", label="Merged stack", type="check", default=True,
         title="Pass through the plain merged stack before the finished picture."),
    dict(id="unroll", label="Unroll into the map", type="check", default=True,
         title="End with the picture unrolling into its polar map."),
    dict(id="captions", label="Captions", type="check", default=True,
         title="Top left over each stage: what it is, its figure (the exposure, the number of frames, the date) "
               "and one sentence; off: the pictures alone."),
    dict(id="fps", label="Frame rate", type="select", default="25",
         values=[["25", "25"], ["30", "30"], ["50", "50"]]),
]


def _ex_txt(t):
    return f"{t:g} s" if t >= 1 else f"1/{round(1 / t):g} s"


def _ease(t):
    return 0.5 - 0.5 * math.cos(math.pi * min(max(t, 0.0), 1.0))


def _resample(src, cy, cx, scale, W, H, dy=0.0, dx=0.0, rot=0.0):
    """`src` (h, w, 3) float around (cy, cx) at `scale` output px per source px,
    the centre landing at (W/2 + dx, H/2 + dy) of a W x H frame; the picture
    direction `rot` (deg clockwise from up) stands up."""
    yy, xx = np.mgrid[0:H, 0:W].astype(np.float32)
    u = (xx - (W / 2.0 + dx)) / scale
    v = (yy - (H / 2.0 + dy)) / scale
    if rot % 360.0:
        t = math.radians(rot)
        ca, sa = math.cos(t), math.sin(t)
        X = cx + u * ca - v * sa
        Y = cy + u * sa + v * ca
    else:
        X, Y = cx + u, cy + v
    out = np.empty((H, W, 3), np.float32)
    for c in range(3):
        out[..., c] = ndimage.map_coordinates(src[..., c], [Y, X], order=1, mode="constant", cval=BGC[c])
    return out


def _pil_frame(a, W, H):
    return Image.fromarray(_ex._to8(a))


class _Writer:
    """Frames to a video file: ffmpeg (system or imageio-ffmpeg), else GIF."""

    def __init__(self, path, W, H, fps):
        self.W, self.H, self.fps = W, H, fps
        self.path = path
        self.proc = None
        self.gif = None
        exe = shutil.which("ffmpeg") or next((p for p in ("/opt/homebrew/bin/ffmpeg", "/usr/local/bin/ffmpeg")
                                               if os.path.exists(p)), None)
        if exe is None:
            try:
                import imageio_ffmpeg
                exe = imageio_ffmpeg.get_ffmpeg_exe()
            except Exception:
                exe = None
        self.exe = exe
        if exe:
            self.proc = subprocess.Popen(
                [exe, "-y", "-loglevel", "error", "-f", "rawvideo", "-pix_fmt", "rgb24", "-s", f"{W}x{H}",
                 "-r", str(fps), "-i", "-", "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "19",
                 "-movflags", "+faststart", path], stdin=subprocess.PIPE, stderr=subprocess.PIPE)
        else:
            self.gif = []
        self.n = 0

    @property
    def kind(self):
        return "mp4" if self.proc is not None else "gif"

    def add(self, im):
        if self.proc is not None:
            self.proc.stdin.write(np.asarray(im.convert("RGB"), np.uint8).tobytes())
        else:
            if self.n % 2 == 0:                                  # GIF at half the rate, half size
                self.gif.append(im.convert("RGB").reduce(2))
        self.n += 1

    def close(self):
        if self.proc is not None:
            self.proc.stdin.close()
            err = self.proc.stderr.read().decode(errors="replace")
            rc = self.proc.wait()
            if rc != 0:
                raise RuntimeError(f"ffmpeg failed: {err.strip()[-300:]}")
            return self.path
        p = os.path.splitext(self.path)[0] + ".gif"
        if self.gif:
            self.gif[0].save(p, save_all=True, append_images=self.gif[1:], duration=int(2000 / self.fps),
                             loop=0, optimize=False)
        return p


def _frames_of_run(wd, folder):
    """One frame per tier, the run's best, shortest exposure first:
    [(sec, path, (cy, cx) in the frame's full-size pixels), ...]."""
    try:
        rep = json.load(open(os.path.join(wd, "report.json")))
        geo = json.load(open(os.path.join(wd, "geometry.json")))
    except Exception:
        return []
    co = geo.get("crop_origin") or [0, 0]
    ash = geo.get("abs_shift") or {}
    out = []
    for t in rep.get("tiers") or []:
        name = t.get("best")
        if not name:
            continue
        p = os.path.join(folder, name)
        if not os.path.exists(p):
            # the run may have read the lights from a subfolder
            hits = [os.path.join(r, name) for r, _d, fs in os.walk(folder) if name in fs]
            if not hits:
                continue
            p = hits[0]
        sec = float(t["sec"])
        key = next((k for k in ash if abs(float(k) - sec) < 1e-9), None)
        sh = ash.get(key, t.get("shift") or [0.0, 0.0])
        # stack pixel = aligned pixel - crop origin; aligned = frame shifted by +shift
        cy = float(geo["cy"]) + float(co[0]) - float(sh[0])
        cx = float(geo["cx"]) + float(co[1]) - float(sh[1])
        out.append((sec, p, (cy, cx)))
    out.sort(key=lambda e: e[0])
    return out


def _find_moon(a, r_px):
    """The Moon's centre (y, x) in a developed frame `a`, found in the frame
    itself: a dark disc of radius `r_px` with brighter light just outside it,
    matched over the whole frame on the log brightness. Returns (y, x, score);
    score is the brightness step at the edge in natural-log units."""
    g = np.log(np.maximum(a.mean(2), 1e-4)).astype(np.float32)
    k = max(1, int(math.ceil(max(g.shape) / 700.0)))
    if k > 1:
        h, w = (g.shape[0] // k) * k, (g.shape[1] // k) * k
        g = g[:h, :w].reshape(h // k, k, w // k, k).mean((1, 3))
    r = r_px / k
    n = int(math.ceil(1.3 * r)) + 2
    yy, xx = np.mgrid[-n:n + 1, -n:n + 1].astype(np.float32)
    rr = np.hypot(yy, xx)
    inner = (rr < 0.92 * r).astype(np.float32)
    outer = ((rr > 1.04 * r) & (rr < 1.25 * r)).astype(np.float32)
    t = outer / outer.sum() - inner / inner.sum()
    from scipy.signal import fftconvolve
    c = fftconvolve(g - g.mean(), t[::-1, ::-1], mode="same")
    # only centres whose disc lies in the frame
    m = int(math.ceil(0.5 * r))
    c[:m] = c[-m:] = -np.inf
    c[:, :m] = c[:, -m:] = -np.inf
    iy, ix = np.unravel_index(int(np.argmax(c)), c.shape)
    # the score: how much brighter just outside the disc than inside it,
    # in natural-log units (1 = e times), from the match itself
    return (iy + 0.5) * k - 0.5, (ix + 0.5) * k - 0.5, float(c[iy, ix])


def _moon_in(a, r0):
    """The Moon in a developed frame, its size not trusted: the best disc
    over radii 0.55..1.8 r0, refined by a parabola in log radius.
    Returns (y, x, r, score)."""
    rs = r0 * np.geomspace(0.55, 1.8, 15)
    res = [_find_moon(a, r) for r in rs]
    sc = np.array([t[2] for t in res])
    i = int(np.argmax(sc))
    r = rs[i]
    if 0 < i < len(rs) - 1:
        l0, l1, l2 = sc[i - 1], sc[i], sc[i + 1]
        den = l0 - 2 * l1 + l2
        if den < 0:
            r = rs[i] * (rs[1] / rs[0]) ** float(np.clip(0.5 * (l0 - l2) / den, -0.5, 0.5))
    y, x, s = _find_moon(a, r)
    return y, x, float(r), float(s)


def _place_frames(runs, w_unc, R, rep, sample, log):
    """The camera frames laid on the film: where the stack says the Moon is
    in each frame, at the stack's scale -- checked on the frames themselves.
    A decoder that hands the frames back at another size or framing than the
    run read them at (or a frame the run moved) would otherwise put the Sun
    off the middle and the Moon at the wrong size. The size is corrected
    once for all frames (the median over the frames where the Moon is
    clear), the place frame by frame."""
    tiers = ((rep.get("moon_track") or {}).get("tiers") or {})
    got = []
    for sec, path, (fcy, fcx) in runs:
        try:
            a, fs = _develop(path, w_unc)
        except Exception:
            continue
        Rt = R
        for k_, t_ in tiers.items():
            try:
                if abs(float(k_) - sec) < 1e-6 * max(1.0, sec) + 1e-7 and t_.get("R"):
                    Rt = float(t_["R"])
            except Exception:
                pass
        try:
            my, mx, mr, msc = _moon_in(a, Rt * fs)
        except Exception:
            my = mx = mr = None
            msc = 0.0
        got.append(dict(sec=sec, a=(np.clip(a, 0, 1) * 255).astype(np.uint8), fs=fs, py=fcy * fs, px=fcx * fs,
                        R=Rt, my=my, mx=mx, q=(mr / (Rt * fs) if mr else None), sc=msc))
        del a
    if not got:
        return []
    h0, w0 = got[0]["a"].shape[:2]
    qs = [g["q"] for g in got if g["q"] and g["sc"] > 0.25]
    q = float(np.median(qs)) if len(qs) >= max(2, len(got) // 3) else 1.0
    if abs(q - 1.0) < 0.06:
        q = 1.0
    moved = 0
    out = []
    for g in got:
        fs = g["fs"] * q
        py_, px_ = g["py"] * q, g["px"] * q
        if g["my"] is not None and g["sc"] > 0.25 and math.hypot(g["my"] - py_, g["mx"] - px_) > 0.15 * g["R"] * fs:
            py_, px_ = g["my"], g["mx"]
            moved += 1
        out.append((g["sec"], sample(g["a"].astype(np.float32) / 255.0, py_, px_, fs)))
        g["a"] = None
    msg = f"film: {len(out)} frames, decoded at {w0}\u00d7{h0}"
    if q != 1.0:
        msg += f"; the Moon in them is {q:.2f}\u00d7 the size the run's frame width gives, laid at the stack's size"
    if moved:
        msg += f"; {moved} placed on the Moon found in them (the stack's position did not match)"
    log(msg)
    return out


def _develop(path, rep_w_uncropped):
    """A frame as shot with a soft shoulder, float 0-1, and the factor its
    pixels are smaller than the run's full-size frame."""
    from . import frameselect as _fs
    im = _fs._pil(path, full=False, mode="soft")
    a = np.asarray(im.convert("RGB"), np.float32) / 255.0
    return a, (a.shape[1] / float(rep_w_uncropped)) if rep_w_uncropped else 1.0


def storyboard(rgb, geo, wd, folder, opts, scale=1.0):
    """The card's preview: three stills in the chosen format side by side --
    a frame stage (the picture, dim), the picture, the map."""
    v = _opts(opts)
    W, H = SIZES[v["aspect"]]
    k = 300.0 / max(W, H)
    w, h = int(W * k), int(H * k)
    rpx = min(w, h) / (5.0 if H > W else 6.0)
    # the picture reduced first (area average) to about twice the stills'
    # scale: sampled straight from the full size it aliased into grain
    src = np.asarray(rgb, np.float32)
    ctx0 = _ex.context(wd, folder, geo, scale=scale, **_ex.angle_kw(opts))
    q = max(1.0, ctx0["Rs"] / (2.0 * rpx))
    if q > 1.05:
        im_ = Image.fromarray(_ex._to8(src))
        im_ = im_.resize((max(2, int(round(im_.width / q))), max(2, int(round(im_.height / q)))), Image.LANCZOS)
        q = src.shape[1] / float(im_.width)
        src = np.asarray(im_, np.float32) / 255.0
    else:
        q = 1.0
    ctx = _ex.context(wd, folder, geo, scale=scale * q, **_ex.angle_kw(opts))
    rot = _rot_of(ctx, v)
    if rot:
        rpx *= _zoom_needed(src, ctx, w, h, rpx, rot)
    sc = rpx / ctx["Rs"]
    pic = _resample(src, ctx["cy"], ctx["cx"], sc, w, h, rot=rot)
    try:
        from .render import plain_stack
        st = max(2, int(round(scale)))
        ps = np.asarray(plain_stack(wd, geometry=geo, step=st), np.float32)
        plain = _resample(ps, float(geo["cy"]) / st, float(geo["cx"]) / st, rpx / (ctx["Rs"] * scale / st), w, h,
                          rot=rot)
    except Exception:
        plain = np.clip(pic, 0, 1) ** 1.6 * 0.9
    # the end of the film: solar north up and the map grown to fill the width
    # (or, in the reel, larger than the frame, as it slides)
    hs, ws = src.shape[:2]
    r_max = max(1.5, min(ctx["cx"], ws - 1 - ctx["cx"], ctx["cy"], hs - 1 - ctx["cy"]) / ctx["Rs"] - 0.05)
    r_end = max(rpx * 1.6 if v["aspect"] == "9:16" else w / (2 * math.pi) * 0.995, rpx)
    r_end = max(r_end, (h * 0.92) / max(0.5, r_max - 1.2))      # the map's top inside the picture's faded edge
    s_top, s_right, s_bot, s_left, s_cw, s_ch = _ex.safe_margins(w, h, v["safe"])
    flat = _flat_frame(src, ctx, w, h, r_end, h - max(int(h * 0.08), s_bot), 0.0,
                       rot=(None if ctx["o"]["origin"] in ("solar", "celestial") else rot))
    pad = 10
    sheet = Image.new("RGB", (3 * w + 4 * pad, h + 2 * pad + 26), _ex.BG)
    d = _ex.Draw(sheet)
    for i, (a, lab) in enumerate(((plain, "frames to the stack"), (pic, "the finished picture"), (flat, "unrolled"))):
        x = pad + i * (w + pad)
        sheet.paste(Image.fromarray(_ex._to8(a)), (x, pad))
        if s_top or s_bot:
            # where the lettering may go: the app covers the rest
            _safe_outline(d, x, pad, w, h, (s_top, s_right, s_bot, s_left, s_cw, s_ch))
        d.tracked((x + w / 2, h + pad + 7), lab.upper(), fill=(165, 165, 165), font=_ex._font(11, True),
                  anchor="mt", track=0.18)
    return sheet


def _safe_outline(d, x, y, w, h, m):
    """A faint dashed line round the safe zone of a w x h still at (x, y)."""
    top, right, bot, left, cw, ch = m
    x0, y0, x1, y1 = x + left, y + top, x + w - right, y + h - bot
    pts = [(x0, y0), (x1, y0)]
    if cw and ch:
        yc = y + h - ch
        pts += [(x1, yc), (x + w - cw, yc), (x + w - cw, y1)]
    else:
        pts += [(x1, y1)]
    pts += [(x0, y1), (x0, y0)]
    for (ax, ay), (bx, by) in zip(pts[:-1], pts[1:]):
        n = max(1, int(math.hypot(bx - ax, by - ay) / 6))
        for k in range(0, n, 2):
            t0, t1 = k / n, min(1.0, (k + 1) / n)
            d.line([(ax + (bx - ax) * t0, ay + (by - ay) * t0), (ax + (bx - ax) * t1, ay + (by - ay) * t1)],
                   fill=(232, 163, 61), width=1)


def _opts(opts):
    v = {}
    for o in OPTIONS:
        x = (opts or {}).get(o["id"], o["default"])
        v[o["id"]] = (x if isinstance(x, bool) else str(x).lower() in ("1", "true", "on", "yes")) \
            if o["type"] == "check" else str(x if x is not None else o["default"])
    if v["aspect"] not in SIZES:
        v["aspect"] = "16:9"
    if v["safe"] not in _ex.SAFE_MARGINS:
        v["safe"] = "off"
    if v["pace"] not in PACE:
        v["pace"] = "calm"
    if v["rotate"] not in [r[0] for r in _ex.ROTATE_VALUES]:
        v["rotate"] = "none"
    return v


def _rot_of(ctx, v, log=None):
    """The picture direction (deg clockwise from up) the film stands up."""
    th = _ex.turn_angle(ctx, v.get("rotate", "none"))
    if th is None:
        if log:
            log(f"film: {v['rotate']} is not known here (no star match): as captured")
        return 0.0
    return float(th) % 360.0


def _unroll(src, ctx, W, H, t, y0, rpx, dx=0.0, rot=None):
    """The picture with the ring of radius `rpx` (film px) opened by `t`
    (0 the picture, 1 the flat strip whose base line is at `y0`); the
    picture direction `rot` (deg clockwise from up in the picture; default
    solar north) stands up on the canvas, so it ends in the middle of the
    strip with the corona upward. `dx` moves the whole thing sideways. The
    seam at the bottom is feathered into the background as it opens."""
    yy, xx = np.mgrid[0:H, 0:W].astype(np.float32)
    k = max((1.0 - t) / rpx, 1e-7)
    Cy = y0 + 1.0 / k
    vx, vy = xx - (W / 2.0 + dx), yy - Cy
    rho = np.hypot(vx, vy)
    alpha = np.arctan2(vx, -vy)                               # clockwise from up at the centre of curvature
    s_ = alpha / k                                            # arc length along the strip
    rr = rpx + (rho - 1.0 / k)                                # radius in film px (rpx at the base line)
    phi = s_ / rpx
    sc = ctx["Rs"] / rpx
    z0 = math.radians(ctx["o"]["zero_cw"] if rot is None else rot)
    ang = z0 + phi                                            # the picture, turned, as it is
    X = ctx["cx"] + np.sin(ang) * rr * sc
    Y = ctx["cy"] - np.cos(ang) * rr * sc
    h, w = src.shape[:2]
    out = np.empty((H, W, 3), np.float32)
    for c in range(3):
        out[..., c] = ndimage.map_coordinates(src[..., c], [Y, X], order=1, mode="nearest")
    if t > 1e-6:
        fe = max(1.5, 0.06 * rpx * math.sqrt(t))
        a_seam = np.clip((math.pi * rpx - np.abs(s_)) / fe, 0.0, 1.0)
    else:
        a_seam = 1.0
    # the picture's own edges fade out over a fifth of a radius instead of
    # cutting (they show as the ring opens)
    fe_p = max(2.0, 0.2 * rpx * sc)
    a_in = (np.clip((X + 1) / fe_p, 0, 1) * np.clip((w - X) / fe_p, 0, 1)
            * np.clip((Y + 1) / fe_p, 0, 1) * np.clip((h - Y) / fe_p, 0, 1))
    a = (a_seam * a_in * (rr >= 0))[..., None]
    return out * a + BGC * (1 - a)


def _flat_frame(src, ctx, W, H, rpx, y0, dx, rot=None):
    return _unroll(src, ctx, W, H, 1.0, y0, rpx, dx=dx, rot=rot)


def _zoom_needed(src, ctx, W, H, rpx, rot):
    """The factor the picture must be enlarged by so that, turned by `rot`
    (deg) around the Sun at the frame's centre, it still covers the frame."""
    h, w = src.shape[:2]
    sc = ctx["Rs"] / rpx
    z = 1.0
    for u, v in ((-W / 2.0, -H / 2.0), (W / 2.0, -H / 2.0), (-W / 2.0, H / 2.0), (W / 2.0, H / 2.0)):
        rr = math.hypot(u, v)
        ang = math.radians(rot) + math.atan2(u, -v)
        dx_, dy_ = math.sin(ang) * rr * sc, -math.cos(ang) * rr * sc
        for d, room in ((dx_, (w - 1 - ctx["cx"]) if dx_ > 0 else ctx["cx"]),
                        (dy_, (h - 1 - ctx["cy"]) if dy_ > 0 else ctx["cy"])):
            if abs(d) > 1e-6 and room > 1:
                z = max(z, abs(d) / room)
    return min(z, 1.8)


def _run_ordered(fn, specs, sink, workers=None):
    """fn(spec) -> frame, computed a few at a time on threads (map_coordinates
    and the big ufuncs release the GIL), handed to `sink` in order."""
    n = workers or max(1, min(4, (os.cpu_count() or 2) - 1))
    with ThreadPoolExecutor(max_workers=n) as ex:
        pending = []
        for sp in specs:
            pending.append(ex.submit(fn, sp))
            if len(pending) >= 2 * n:
                sink(pending.pop(0).result())
        for f in pending:
            sink(f.result())



def _match_tones(a, ref, r_in=0.0):
    """`a` brought to the radial brightness and colour of `ref`: both frames
    have the Sun in the middle; per ring of 2 px and per channel, a is scaled
    so its median equals ref's median there. The structure (every streamer,
    every gap) stays a's own; only how bright each distance is, and its
    colour, become the finished picture's.

    Only where BOTH hold picture. The film's frame is larger than the picture
    on one side, and a ring that runs out of the picture used to count the
    empty part as black: where about half of a ring lay outside, each median
    fell to its picture's darkest pixels there, their ratio jumped, and the
    stretched stack got a dark circle at that distance. Rings with less than
    60 % of their pixels in both pictures take the gain of the last ring
    that had them."""
    H, W = a.shape[:2]
    yy, xx = np.mgrid[0:H, 0:W]
    ring = (np.hypot(yy - H / 2.0, xx - W / 2.0) / 2.0).astype(np.int32)
    n = int(ring.max()) + 1
    idx = np.arange(n)
    both = (a.max(axis=-1) > 1e-4) & (ref.max(axis=-1) > 1e-4)
    k0 = min(n - 1, int(r_in / 2.0) + 2)                 # inside the Moon: the gain at its edge
    both |= ring < k0
    cnt = np.bincount(ring.ravel(), minlength=n).astype(np.float64)
    cov = np.bincount(ring.ravel(), weights=both.ravel().astype(np.float64), minlength=n) / np.maximum(cnt, 1.0)
    ok = (cov >= 0.6) & (cnt * cov >= 40)
    lab = np.where(both, ring, -1)
    out = np.empty_like(a)
    for c in range(3):
        ma = np.asarray(ndimage.median(a[..., c], lab, idx), np.float64)
        mr = np.asarray(ndimage.median(ref[..., c], lab, idx), np.float64)
        lg = np.log(np.maximum(mr, 1e-4) / np.maximum(ma, 1e-4))
        lg[~np.isfinite(lg)] = 0.0
        lg[:k0] = lg[k0]
        good = ok.copy()
        good[:k0 + 1] = True
        last = np.maximum.accumulate(np.where(good, idx, 0))
        lg = lg[last]                                    # a ring without enough picture: the last one that had it
        # close to the limb the brightness falls fast and the gain follows
        # ring by ring; further out it is taken over 10 rings, or single
        # rings print as faint circles in the smooth outer field
        g2 = ndimage.gaussian_filter1d(lg, 2.0, mode="nearest")
        g10 = ndimage.gaussian_filter1d(lg, 10.0, mode="nearest")
        wf = np.clip((idx - 1.3 * k0) / max(1.0, 0.7 * k0), 0.0, 1.0)
        g = g2 * (1.0 - wf) + g10 * wf
        out[..., c] = np.clip(a[..., c] * np.exp(g)[ring], 0.0, 1.0)
    return out


def make(ly, params, wd, folder, opts, out_base, prog=None, render=None, plain_stack=None):
    """Writes the film next to `out_base`. Returns (path, seconds)."""
    v = _opts(opts)
    W, H = SIZES[v["aspect"]]
    fps = int(v["fps"])
    pace = PACE[v["pace"]]
    log = (lambda m, f=None: prog.log(m, f) if prog else None)
    geo = json.load(open(os.path.join(wd, "geometry.json")))
    rep = {}
    try:
        rep = json.load(open(os.path.join(wd, "report.json")))
    except Exception:
        pass
    ctx = _ex.context(wd, folder, geo, scale=1.0, **_ex.angle_kw(opts))
    rpx = min(W, H) / (5.0 if H > W else 6.0)                 # one unit radius in film px
    sc = rpx / ctx["Rs"]
    # ---- the finished picture, full resolution, kept at about twice the film scale
    log("film: rendering the finished picture ...", 0.05)
    full = np.asarray(render(ly, params, preview=False, view="composite"), np.float32)
    q = max(1.0, 1.0 / (2.0 * sc))                             # source at twice the film scale
    if q > 1.05:
        im = Image.fromarray(_ex._to8(full))
        im = im.resize((max(2, int(round(im.width / q))), max(2, int(round(im.height / q)))), Image.LANCZOS)
        q = full.shape[1] / float(im.width)
        src = np.asarray(im, np.float32) / 255.0
    else:
        q = 1.0
        src = full
    del full
    ctxs = _ex.context(wd, folder, geo, scale=float(q), **_ex.angle_kw(opts))
    rot = _rot_of(ctxs, v, log)
    if rot:
        # turned, the picture grows just enough to keep the frame covered
        rpx *= _zoom_needed(src, ctxs, W, H, rpx, rot)
    yA = H / 2.0 - rpx
    final = _unroll(src, ctxs, W, H, 0.0, yA, rpx, rot=rot)      # the picture as it is (or turned)
    # ---- the plain stack
    plain = None
    if v["stack"] and plain_stack is not None:
        try:
            log("film: the plain merged stack ...", 0.15)
            ps = np.asarray(plain_stack(wd, size="half", geometry=geo), np.float32)
            plain = _resample(ps, ctx["cy"] / 2.0, ctx["cx"] / 2.0, rpx / (ctx["Rs"] / 2.0), W, H, rot=rot)
            del ps
        except Exception as e:
            log(f"film: no plain stack ({e})")
    # ---- the frames
    cam = []
    if v["frames"]:
        log("film: developing the frames ...", 0.2)
        ac = (rep.get("autocrop_px") or {})
        w_unc = float(rep.get("W") or 0) + float(ac.get("left") or 0) + float(ac.get("right") or 0)
        cam = _place_frames(_frames_of_run(wd, folder), w_unc, float(geo.get("R", ctx["R"])), rep,
                            lambda a, cy_, cx_, fs: _resample(a, cy_, cx_, rpx / (ctx["Rs"] * fs), W, H, rot=rot),
                            log)
    n_fr = int(rep.get("n_files") or 0)
    n_exp = len(rep.get("tiers") or [])
    # ---- write
    out_path = out_base + ".mp4"
    kk = 0
    while os.path.exists(out_path) or os.path.exists(os.path.splitext(out_path)[0] + ".gif"):
        kk += 1
        out_path = out_base + f"_{kk}.mp4"
    wr = _Writer(out_path, W, H, fps)
    fs_cap = max(14, int(min(W, H) / 26))
    caps = v["captions"]
    total = [0.0]
    # THE CAPTION, in the sheet's voice: a kicker (letter-spaced capitals, the
    # accent colour), the figure or name large and light, one sentence in
    # italic; top left over a soft dark gradient instead of a bar. A caption
    # is drawn once and laid over every frame that carries it.
    # the safe zone (reels, stories, posts): the lettering starts below the
    # app's header and stays clear of its sides; the darkening grows with it
    s_top, s_right, s_bot, s_left, _cw, _ch = _ex.safe_margins(W, H, v["safe"])
    mx, y0 = max(int(fs_cap * 0.9), s_left), max(int(fs_cap * 0.75), s_top)
    mr = max(int(fs_cap * 0.9), s_right)
    scrim_h = int(y0 + fs_cap * 6.25)
    ramp = 1.0 - np.arange(scrim_h, dtype=np.float32) / scrim_h          # 1 at the top, 0 at its foot
    scrim = Image.fromarray(np.repeat((ramp * ramp * (3 - 2 * ramp) * 0.5 * 255).astype(np.uint8)[:, None], W, 1), "L")
    f_k, f_t, f_d = _ex._font(fs_cap * 0.5, True), _ex._font(fs_cap * 1.25), _ex._serif(fs_cap * 0.68)
    cap_cache = {}

    def cap_layer(label, sub, bar):
        key = (label, sub, bar)
        if key not in cap_cache:
            if len(cap_cache) > 64:
                cap_cache.clear()
            kick, title = label if isinstance(label, tuple) else ("", label)
            lay = Image.new("RGBA", (W, scrim_h), (0, 0, 0, 0))
            d = _ex.Draw(lay)
            halo = {}        # no dark outline: fading in, it would show before the letters do
            if kick:
                d.tracked((mx, y0), kick.upper(), fill=_ex.GOLD + (255,), font=f_k, anchor="la", track=0.22)
            ft, fd = f_t, f_d
            if title and d.textlength(title, font=ft) > W - mx - mr:
                ft = _ex._font(fs_cap * 1.25 * (W - mx - mr) / d.textlength(title, font=ft))
            sub_lines = [sub] if sub else []
            if sub and d.textlength(sub, font=fd) > W - mx - mr:
                # too long for one line (a narrow safe zone): two lines at
                # the same size, split at the space nearest the middle that
                # fits; smaller only when even that does not
                words = sub.split(" ")
                best = None
                for k in range(1, len(words)):
                    a_, b_ = " ".join(words[:k]), " ".join(words[k:])
                    wmax = max(d.textlength(a_, font=fd), d.textlength(b_, font=fd))
                    if best is None or wmax < best[0]:
                        best = (wmax, [a_, b_])
                if best is not None:
                    sub_lines = best[1]
                    if best[0] > W - mx - mr:
                        fd = _ex._serif(fs_cap * 0.68 * (W - mx - mr) / best[0])
                else:
                    fd = _ex._serif(fs_cap * 0.68 * (W - mx - mr) / d.textlength(sub, font=fd))
            if title:
                d.text((mx, y0 + fs_cap * 0.95), title, fill=_ex.INK + (255,), font=ft, anchor="la", **halo)
            for i_, ln in enumerate(sub_lines):
                d.text((mx, y0 + fs_cap * (2.7 + 1.05 * i_)), ln, fill=(215, 215, 215, 255), font=fd, anchor="la",
                       **halo)
            if bar is not None and len(cam) > 1:
                # which exposure of the run this is: one hairline per frame
                bw, gp, bh = max(2, int(fs_cap * 0.1)), max(3, int(fs_cap * 0.2)), int(fs_cap * 0.55)
                x0 = W - mr - len(cam) * (bw + gp) + gp
                for j in range(len(cam)):
                    x = x0 + j * (bw + gp)
                    d.rectangle([x, y0, x + bw, y0 + bh], fill=_ex.GOLD + (255,) if j <= bar else (90, 90, 90, 255))
            cap_cache[key] = lay
        return cap_cache[key]

    def emit(a, label="", sub="", alpha=1.0, bar=None):
        im = Image.fromarray(_ex._to8(a))
        if caps and (label or sub) and alpha > 0:
            lay = cap_layer(label, sub, bar)
            if alpha >= 0.999:
                im.paste((0, 0, 0), (0, 0, W, scrim_h), scrim)
                im.paste(lay, (0, 0), lay)
            else:
                im.paste((0, 0, 0), (0, 0, W, scrim_h), scrim.point(lambda p_: int(p_ * alpha)))
                im.paste(lay, (0, 0), lay.getchannel("A").point(lambda p_: int(p_ * alpha)))
        wr.add(im)
        total[0] += 1.0 / fps

    def hold(a, sec, *cap):
        for _ in range(int(round(fps * sec))):
            emit(a, *cap)

    def xfade(a, b, sec, cap_a, cap_b):
        n = max(1, int(round(fps * sec)))
        for j in range(n):
            t = _ease((j + 1) / n)
            emit(a * (1 - t) + b * t, *(cap_b if t > 0.5 else cap_a))

    log("film: writing the frames ...", 0.35)
    # kicker, title, one sentence. The finished picture carries the sheet's
    # own title: what it is and the date
    kind = _ex._eclipse_kind(ctx)
    date = _ex._long_date(ctx["date"]) if ctx.get("date") else ""
    FI = ((kind if date else "the finished picture", date or "Detail layers and the look"),
          "The finished picture: the detail layers and the look." if date else "", 1.0, None)
    PL = (("the merged stack", f"{n_fr} frames, {n_exp} exposures" if n_fr and n_exp else "Aligned and merged"),
          "Aligned and merged into one; no look yet.", 1.0, None)
    if cam:
        CAM = lambda i: ((f"single frame {i + 1} of {len(cam)}", _ex_txt(cam[i][0])),
                         "Straight from the camera, one exposure after the other.", 1.0, i)
        hold(cam[0][1], 0.8 * pace, *CAM(0))
        for i in range(len(cam) - 1):
            xfade(cam[i][1], cam[i + 1][1], 0.45 * pace, CAM(i), CAM(i + 1))
            hold(cam[i + 1][1], 0.28 * pace, *CAM(i + 1))
        last, lastcap = cam[-1][1], CAM(len(cam) - 1)
    else:
        last, lastcap = None, None
    if plain is not None:
        if last is not None:
            xfade(last, plain, 1.6 * pace, lastcap, PL)
        hold(plain, 2.5 * pace, *PL)
        last, lastcap = plain, PL
        # the stack stretched to the picture's tones (and colour) before the
        # detail layers come in, so the last dissolve shows only what the
        # layers and the look add
        ST = (("the merged stack, stretched", "The picture\u2019s tones"),
              "Brightness and colour of the finished picture, no detail layers yet.", 1.0, None)
        stretched = _match_tones(plain, final, rpx * ctx["Rm"] / ctx["Rs"])
        xfade(plain, stretched, 1.4 * pace, PL, ST)
        hold(stretched, 1.3 * pace, *ST)
        last, lastcap = stretched, ST
    if last is not None:
        xfade(last, final, 1.8 * pace, lastcap, FI)
    hold(final, 3.0 * pace, *FI)
    nf = int(round(fps * 0.7))
    for j in range(nf):
        emit(final, FI[0], FI[1], 1.0 - (j + 1) / nf, None)
    hold(final, 0.6 * pace)
    if v["unroll"]:
        log("film: unrolling ...", 0.6)
        # the turn: the picture turns until the angle origin (solar north)
        # stands up, the shorter way round, so the strip opens at its south
        # (Without a star match or an angle entered by hand there is no north:
        # the picture stays as it is shown, turned or not, and that up is the
        # middle of the map.)
        known = ctxs["o"]["origin"] in ("solar", "celestial")
        rot_m = None if known else rot                       # the direction in the map's middle
        z = float(ctxs["o"]["zero_cw"]) % 360.0
        delta = ((z - rot + 180.0) % 360.0) - 180.0 if known else 0.0
        if abs(delta) >= 1.0:
            TU = (("turned", ("Solar" if ctx["o"]["origin"] == "solar" else "Celestial") + " north up"),
                  f"Turned {abs(delta):.0f}\u00b0 {'counter-clockwise' if delta > 0 else 'clockwise'} about the Sun.",
                  1.0, None)
            Nr = int(round(fps * (1.2 + 1.2 * abs(delta) / 90.0) * pace))
            # the picture grows as it turns, just enough to keep the frame
            # covered (a turned rectangle leaves corners bare)
            tts = [_ease((j + 1) / Nr) for j in range(Nr)]
            zn = [_zoom_needed(src, ctxs, W, H, rpx, rot + delta * tt) for tt in tts]
            z_end = max(zn) if zn else 1.0
            zs, zm = [], 1.0
            for tt, zq in zip(tts, zn):
                zm = max(zm, zq, 1.0 + (z_end - 1.0) * tt)
                zs.append(zm)
            _run_ordered(lambda a: (_unroll(src, ctxs, W, H, 0.0, H / 2.0 - rpx * a[1], rpx * a[1], rot=rot + delta * a[0]), a[0]),
                         list(zip(tts, zs)),
                         lambda r: emit(r[0], TU[0], TU[1], min(1.0, 4.0 * r[1]), None))
            rpx = rpx * z_end
            yA = H / 2.0 - rpx
            turned = _unroll(src, ctxs, W, H, 0.0, yA, rpx, rot=rot + delta)
            hold(turned, 0.8 * pace, *TU)
            nf = int(round(fps * 0.5))
            for j in range(nf):
                emit(turned, TU[0], TU[1], 1.0 - (j + 1) / nf, None)
        # stage 4: the ring opens at the picture's own scale; the base line
        # moves from the top of the Sun down to a line near the bottom
        yB = H - max(60, int(H * 0.08), s_bot)
        N = int(round(fps * 3.4 * pace))
        _run_ordered(lambda t: _unroll(src, ctxs, W, H, t, yA + (yB - yA) * t, rpx, rot=rot_m),
                     [_ease((j + 1) / N) for j in range(N)], emit)
        # stage 5: into place. 16:9 and 1:1: grow until the strip fills the
        # width. 9:16: grow to 1.6x and slide from one end to the other.
        wide = v["aspect"] == "9:16"
        yB2 = H - max(60, int(H * 0.07), s_bot)
        # the map's top must stay inside the picture: the frame reaches
        # r = yB2 / rpx_end at the top, the picture r_max all round
        hs, ws = src.shape[:2]
        r_max = max(1.5, min(ctxs["cx"], ws - 1 - ctxs["cx"], ctxs["cy"], hs - 1 - ctxs["cy"]) / ctxs["Rs"] - 0.05)
        rpx_end = (rpx * 1.6) if wide else max(rpx, W / (2 * math.pi) * 0.995)
        # r = 1 at the base line; the picture's edges fade out over a fifth of
        # a radius, and no bar hides the map's top now: it stays inside that
        rpx_end = max(rpx_end, yB2 / max(0.5, r_max - 1.2))
        if not wide:
            wide = 2 * math.pi * rpx_end > W * 1.15          # too wide for the frame: slide as the reel does
        N2 = int(round(fps * 2.0 * pace))
        half_span = (2 * math.pi * rpx_end - W) / 2.0 if wide else 0.0
        _run_ordered(lambda t: _flat_frame(src, ctxs, W, H, rpx + (rpx_end - rpx) * t, yB + (yB2 - yB) * t,
                                           dx=half_span * t, rot=rot_m),
                     [_ease((j + 1) / N2) for j in range(N2)], emit)
        org = ctx["o"]["origin"]
        UN = (("unrolled", "The same picture, as a map"),
              ("Position angle across, distance from the Sun upward; "
               + ("solar north" if org == "solar" else ("celestial north" if org == "celestial" else
                  ("the picture\u2019s top" if not rot else "the top of the turned picture")))
               + " in the middle."), 1.0, None)
        if wide:
            log("film: the slide ...", 0.8)
            N3 = int(round(fps * 7.0 * pace))
            _run_ordered(lambda t: (_flat_frame(src, ctxs, W, H, rpx_end, yB2, dx=half_span - 2 * half_span * t,
                                                rot=rot_m), t),
                         [_ease((j + 1) / N3) for j in range(N3)],
                         lambda r: emit(r[0], UN[0], UN[1], min(1.0, 3.0 * r[1]), None))
            lastf = _flat_frame(src, ctxs, W, H, rpx_end, yB2, dx=-half_span, rot=rot_m)
        else:
            lastf = _flat_frame(src, ctxs, W, H, rpx_end, yB2, dx=0.0, rot=rot_m)
            nf = int(round(fps * 0.6))
            for j in range(nf):
                emit(lastf, UN[0], UN[1], (j + 1) / nf, None)
        hold(lastf, 3.5 * pace, *UN)
    else:
        hold(final, 1.5 * pace)
    log("film: encoding ...", 0.92)
    path = wr.close()
    return path, total[0], wr.kind
