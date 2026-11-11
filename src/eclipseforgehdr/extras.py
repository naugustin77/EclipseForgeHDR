"""Extras (tab 5, 0.30.0): products made from the finished picture.

Every product here is drawn from what the run produced -- the composite,
its geometry (Moon centre and radius), the star match when there is one (sky
orientation, scale, Sun altitude), the site and the frame times. Nothing is
modelled or added; where a figure is not in the data (no star match: no
scale in solar radii) the product says so or leaves it out.

THE CATALOGUE lists the products with their options; the page builds its
cards from it, so a new product is one entry and one function here.

THE FRAME. Products work in the stack's own frame, as the layers hold it;
the Orientation setting of the Export tab does not apply to them.

SCALE. `make(..., scale=q)` draws on a picture q times smaller than the stack
(the page's preview layers are 1/4); the geometry is divided by q and every
size in the drawing follows the picture's width, so the preview shows the
product as the file will look.

OUTPUT. Presentation products: 8-bit JPG or PNG, and transparent PNG for the
layers. The 16-bit path for the picture itself is the Export tab.
"""
import json
import math
import os

import numpy as np
from PIL import Image, ImageDraw, ImageFont
from scipy import ndimage

from . import polarmap as _pm

GOLD = (232, 163, 61)
BG = (0, 0, 0)
KM_PER_RSUN = 695700.0
EARTH_DIAM_KM = 12742.0

# ---------------------------------------------------------------- catalogue
# type: check | select | text ; values: [[value, label], ...] for select
ASPECTS = {"3:2": (3, 2), "4:3": (4, 3), "5:4": (5, 4), "16:9": (16, 9), "1:1": (1, 1),
           "4:5": (4, 5), "3:4": (3, 4), "2:3": (2, 3), "9:16": (9, 16)}
FILES_OPT = dict(id="files", label="Files", type="select", default="one",
                 values=[["one", "the picture"], ["layers", "picture + layers + SVG"]],
                 title="Picture + layers: besides the finished file, the picture without lettering, the lettering "
                       "alone as a transparent PNG, and an SVG with the picture linked under the lettering as "
                       "vectors and live text (Inter), to open in a vector drawing app and work on there (a PDF/X-4 "
                       "for print is then exported from that app).")
# THE SOCIAL FORMATS (2026-11-10). The sizes the apps show without
# scaling, 1080 px wide, and the parts of the frame their own buttons and text
# cover. Measured off the published layout charts, in px of a 1080-wide frame:
# top, right, bottom, left, and the column of buttons at the right edge
# (its width from the edge, its height from the bottom). "reel" is the union
# of the Instagram reel, TikTok and Shorts charts (the widest margin of each
# side): Instagram 250 top / 320 bottom / 120 sides, column 768 high; TikTok
# 160 / 480 / 120, column 1080 high and 300 from the edge; Shorts 240 / 380 /
# 60-120, column 1200 high. A reel shown in the feed is cut to its middle 3:4
# (240 px off top and bottom), which the 250 px top already clears. "story":
# 250 px top and bottom. "post": posts and carousels, 120 px at the right for
# the counter and the buttons, 180 at the bottom.
SOCIAL_SIZES = {"16:9": (1920, 1080), "1:1": (1080, 1080), "4:5": (1080, 1350), "3:4": (1080, 1440),
                "9:16": (1080, 1920)}
SOCIAL_VALUES = [["16:9", "16:9 (1920 × 1080)"], ["1:1", "1:1 post (1080 × 1080)"],
                 ["4:5", "4:5 post (1080 × 1350)"], ["3:4", "3:4 post (1080 × 1440)"],
                 ["9:16", "9:16 reel / story (1080 × 1920)"]]
SAFE_MARGINS = {"reel": (250, 120, 480, 120, 300, 1200), "story": (250, 60, 250, 60, 0, 0),
                "post": (60, 120, 180, 60, 0, 0)}
SAFE_OPT = dict(id="safe", label="Safe zone", type="select", default="off",
                values=[["off", "off"], ["reel", "reels (Instagram, TikTok, Shorts)"], ["story", "stories"],
                        ["post", "posts and carousels"]],
                title="Keeps all lettering clear of what the app lays over the frame: its buttons, the account "
                      "line and the caption. Reels: 250 px at the top, 480 at the bottom, 120 at the sides and "
                      "the column of buttons at the lower right. Stories: 250 at top and bottom. Posts and "
                      "carousels: the counter and buttons at the right. For the square and upright formats.")


def _aspect(frame):
    """(aw, ah) of a frame name ("4:5", or any "W:H"), None for "as"/natural."""
    if frame in ASPECTS:
        return ASPECTS[frame]
    try:
        a_, b_ = str(frame).split(":")
        a_, b_ = float(a_), float(b_)
        return (a_, b_) if a_ > 0 and b_ > 0 else None
    except Exception:
        return None


SOCIAL_OPT = dict(id="social", label="For", type="select", default="full",
                  values=[["full", "full resolution"]] + [x for x in SOCIAL_VALUES],
                  title="Full resolution: the picture's own pixels, in the shape set by Frame. A social format: "
                        "exactly that size (1080 px wide), the lettering laid out for it and large enough to read "
                        "on a phone; Frame is then set by the format.")
#: body lettering of a still at a social size, px (about 8 pt on a phone)
SOCIAL_FS = 17.0


def _scale_ctx(ctx, f):
    """The context of the same picture enlarged by `f` (every pixel figure)."""
    c = dict(ctx)
    for k in ("cx", "cy", "R", "Rm", "Rs"):
        c[k] = ctx[k] * f
    c["scale"] = ctx["scale"] / f
    if ctx.get("km_px"):
        c["km_px"] = ctx["km_px"] * f
    c["stars"] = [dict(st, x=st["x"] * f, y=st["y"] * f) for st in ctx["stars"]]
    return c


def _shrink(rgb, ow, oh):
    """The picture resampled to ow x oh (area average when smaller)."""
    im = Image.fromarray(_to8(np.asarray(rgb, np.float32)), "RGB") if rgb.dtype != np.uint8 else Image.fromarray(rgb)
    # 16-bit detail is not needed at 1080 px: 8-bit after the look, as every file
    return np.asarray(im.resize((int(ow), int(oh)), Image.LANCZOS), np.float32) / 255.0


def safe_margins(W, H, kind):
    """(top, right, bottom, left, col_w, col_h) in px of a W x H frame for the
    safe zone `kind`; zeros when off or for a landscape frame."""
    m = SAFE_MARGINS.get(str(kind or "off"))
    if not m or W > H:
        return (0, 0, 0, 0, 0, 0)
    k = W / 1080.0
    return tuple(int(round(x * k)) for x in m)


# the frame choices on the cards: landscape, square, portrait
FRAME_VALUES = [["3:2", "3:2"], ["4:3", "4:3"], ["5:4", "5:4"], ["16:9", "16:9"], ["1:1", "1:1"],
                ["4:5", "4:5 portrait"], ["3:4", "3:4 portrait"], ["2:3", "2:3 portrait"], ["9:16", "9:16 portrait"]]

CATALOGUE = [
    dict(id="overlay", title="Sky & Sun overlay", group="picture",
         out="JPG/PNG + transparent PNG layer",
         desc="The finished picture with the sky drawn on it: radius rings, the position-angle ring, "
              "the solar axis and equator, a hairline to every prominence with its height, the catalogue "
              "stars by name, a compass with celestial north and the zenith, a scale bar with the Earth to "
              "scale, and the title with the facts of the shot. The layer alone is written too, to lay over "
              "your own export.",
         options=[
             dict(id="rings", label="Radius rings", type="check", default=True,
                  title="Faint circles at 1.5, 2 and 2.5 radii, labelled."),
             dict(id="paring", label="Angle ring", type="check", default=True,
                  title="A ring with ticks every 5 degrees and labels every 30, position angle from north."),
             dict(id="axis", label="Solar axis", type="check", default=True,
                  title="The Sun's rotation axis (solid) and equator (dashed), with N, E, S, W at the ring."),
             dict(id="stars", label="Stars", type="check", default=True,
                  title="The stars matched to the catalogue, circled and named with their magnitude."),
             dict(id="proms", label="Prominences", type="check", default=True,
                  title="A hairline from every prominence found along the limb to its height above the Moon's edge, "
                        "in kilometres and in Earths (needs the star match for the scale)."),
             dict(id="compass", label="Compass", type="check", default=True,
                  title="Celestial north and east, and the direction to the zenith, in a corner."),
             dict(id="scale", label="Scale bar", type="check", default=True,
                  title="A 500 000 km bar and the Earth at the same scale."),
             dict(id="info", label="Facts", type="check", default=True,
                  title="Date, mid-totality time, the Sun's altitude, the site's coordinates, and how north was found."),
             dict(id="turn", label="Turn", type="select", default="none",
                  values=[["none", "as captured"], ["solar", "solar north up"], ["celestial", "celestial north up"],
                          ["zenith", "zenith up"], ["cw90", "90\u00b0 clockwise"], ["ccw90", "90\u00b0 counter-clockwise"],
                          ["r180", "180\u00b0"]],
                  title="Turns the picture around the Sun before drawing. North and zenith need the star match "
                        "(or north by hand); the layer then fits this turned picture, not the plain export."),
             dict(id="frame", label="Frame", type="select", default="as",
                  values=[["as", "as captured"]] + FRAME_VALUES,
                  title="Crops to this shape with the Sun in the middle, the largest crop that fits inside the "
                        "(turned) picture. As captured keeps the whole picture when it is not turned."),
             dict(id="size", label="Lettering", type="select", default="1",
                  values=[["0.8", "small"], ["1", "normal"], ["1.3", "large"]],
                  title="The size of the lettering relative to the picture."),
             SOCIAL_OPT,
             SAFE_OPT,
             FILES_OPT,
             dict(id="fmt", label="File", type="select", default="jpg",
                  values=[["jpg", "JPG"], ["png", "PNG"]], title="The format of the composite picture."),
         ]),
    dict(id="poster", title="Poster", group="picture", out="JPG/PNG",
         desc="The round picture above, its own unrolled polar map below, joined by a soft fade: the corona "
              "as seen and as a curtain of streamers with their position angles.",
         options=[
             dict(id="height", label="Map up to", type="select", default="2.7",
                  values=[["2", "2 radii"], ["2.7", "2.7 radii"], ["3.5", "3.5 radii"], ["frame", "the frame"]],
                  title="How far out the unrolled map reaches."),
             dict(id="turn", label="Turn", type="select", default="solar",
                  values=[["solar", "solar north up"], ["celestial", "celestial north up"], ["zenith", "zenith up"],
                          ["none", "as captured"]],
                  title="Turns the picture so this direction is up; the map has the same direction in its middle, so "
                        "the strip reads as the picture opened at the bottom. Needs the star match (or north by hand)."),
             dict(id="frame", label="Frame", type="select", default="as",
                  values=[["as", "natural"]] + [f_ for f_ in FRAME_VALUES if f_[0] in ("1:1", "4:5", "3:4", "2:3", "9:16")],
                  title="The poster's shape (the picture over its map is upright by nature). A taller frame shows "
                        "more corona above the picture first; what is still missing is black, the picture's top "
                        "fading into it."),
             dict(id="fade", label="Join", type="select", default="soft",
                  values=[["soft", "soft fade"], ["hard", "hard edge"]],
                  title="How the picture passes into the map."),
             dict(id="caption", label="Lettering", type="check", default=True,
                  title="Position angles under the map, the title (what, when, one line of facts) and a small "
                        "line of what the sheet shows and was made from; off: the picture alone."),
             dict(id="maplook", label="Map look", type="select", default="same",
                  values=[["same", "as the picture"], ["standard", "standard (the picture soft)"]],
                  title="With a soft look in tab 3, the map can come from its standard variant: the round picture "
                        "bright and soft (with the earthshine), the map below with the inner corona's full detail "
                        "and contrast. Only the settings in which the two variants differ change; the colour and "
                        "white stay. With a look that is not a soft one, the map uses the picture's look."),
             SOCIAL_OPT,
             SAFE_OPT,
             FILES_OPT,
             dict(id="fmt", label="File", type="select", default="jpg",
                  values=[["jpg", "JPG"], ["png", "PNG"]]),
         ]),
    dict(id="limb", title="Limb panorama", group="picture", out="JPG/PNG",
         desc="The inner corona unwrapped along the Moon's edge at the picture's own resolution, in rows, "
              "position angle along the bottom: every prominence and the base of every streamer in one strip.",
         options=[
             dict(id="outer", label="Out to", type="select", default="1.32",
                  values=[["1.2", "1.2 radii"], ["1.32", "1.32 radii"], ["1.5", "1.5 radii"]],
                  title="The outer edge of the strip."),
             dict(id="rows", label="Rows", type="select", default="4",
                  values=[["2", "2 (180° each)"], ["3", "3 (120° each)"], ["4", "4 (90° each)"]],
                  title="The circle cut into this many rows."),
             dict(id="zoom", label="Zoom", type="select", default="1",
                  values=[["1", "1:1"], ["2", "2:1"]], title="Enlargement of the strip."),
             dict(id="marks", label="Mark prominences", type="check", default=True,
                  title="A marker above every prominence found along the limb."),
             FILES_OPT,
             dict(id="fmt", label="File", type="select", default="jpg",
                  values=[["jpg", "JPG"], ["png", "PNG"]]),
         ]),
]

COMING = []


# THE FEW OPTIONS A CARD SHOWS (0.30.0, the review): everything else sits
# under "More options" on the card, with its default chosen so the product
# looks right without opening it.
BASIC = {
    "overlay": ["turn", "frame", "social", "safe"],
    "poster": ["turn", "frame", "maplook", "social", "safe"],
    "limb": ["outer", "rows"],
    "activity": ["view", "scale", "grain"],
    "film": ["aspect", "safe", "rotate", "pace"],
    "sound": ["turn", "aspect", "safe", "rotate"],
}


def _mark_basic(e):
    keep = BASIC.get(e["id"])
    if keep is None:
        return e
    e = dict(e)
    e["options"] = [dict(o, adv=(o["id"] not in keep)) for o in e.get("options", [])]
    return e


def entries():
    return [_mark_basic(e) for e in _entries()]


def _entries():
    """CATALOGUE plus the products that live in their own modules."""
    from . import film as _film
    from . import sound as _sound
    from . import activity as _act
    return list(CATALOGUE) + [
        dict(id="activity", title="Corona in motion", group="motion", out="GIF + MP4",
             desc="The brackets of a repeated set as a loop, steadied on the corona, every bracket through the "
                  "app's partial convolution (the same masks and weights for all) on mid grey, the Moon one flat "
                  "disc: loops, blobs and waves that moved between the brackets move. Needs the stacking mode "
                  "'HDR per bracket'.",
             options=list(_act.OPTIONS)),
        dict(id="film", title="Film", group="motion", out="MP4 (H.264)",
             desc="From the single frames through the merged stack to the finished picture, which turns so that "
                  "north is up and then unrolls into its polar map. 16:9, 1:1 or a 9:16 reel, where the map "
                  "slides across the frame. The card shows three stills; the clip takes a minute or two to make.",
             options=list(_film.OPTIONS)),
        dict(id="sound", title="Corona sound", group="motion", out="WAV + MP4",
             desc="A line sweeps once around the Sun; the brightness along it, relative to the median corona at the "
                  "same distance, is the loudness of one note per distance: low near the limb, high far out. "
                  "Streamers sound, gaps go quiet, prominences ring a bell. Nothing but the picture's numbers.",
             options=list(_sound.OPTIONS))]


def entry(kind):
    return next((g for g in entries() if g["id"] == kind), None)


def catalogue(wd, folder, have_layers):
    """The products with `ready` and a note, for the page."""
    out = []
    try:
        cat = (json.load(open(os.path.join(wd, "stars.json"))).get("catalogue") or {})
    except Exception:
        cat = {}
    solved = cat.get("status") == "solved"
    for g in entries():
        e = dict(g, ready=bool(have_layers), coming=False)
        if not have_layers:
            e["why"] = "needs a result (tab 2)"
        elif g["id"] == "activity":
            from .activity import available
            n = len(available(wd))
            if n < 2:
                e["ready"] = False
                e["why"] = "needs a run in the stacking mode 'HDR per bracket' on a repeated-bracket set"
            else:
                from .activity import _index, NEEDS_RERUN
                if _index(wd).get("win_box"):
                    e["why"] = f"{n} brackets merged on their own in the work folder"
                else:
                    e["ready"] = False
                    e["why"] = NEEDS_RERUN
        elif not solved:
            e["why"] = ("no star match: angles from image up, radii in Moon radii"
                        if g["id"] != "poster" else "no star match: angles from image up")
        else:
            e["why"] = ""
        out.append(e)
    for g in COMING:
        out.append(dict(g, options=[], out="", ready=False, coming=True, why="next round"))
    return out


def option_values(kind, opts):
    """The options of `kind` with defaults filled in and values typed."""
    g = entry(kind)
    if g is None:
        raise ValueError(f"unknown product {kind!r}")
    v = {}
    for o in g["options"]:
        x = (opts or {}).get(o["id"], o["default"])
        if o["type"] == "check":
            v[o["id"]] = (x if isinstance(x, bool) else str(x).lower() in ("1", "true", "on", "yes"))
        else:
            v[o["id"]] = str(x if x is not None else o["default"])
    return v


# ---------------------------------------------------------------- drawing helpers
_FONT_CANDIDATES = [
    "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
    "/Library/Fonts/Arial Unicode.ttf",
    "/System/Library/Fonts/Supplemental/Arial Unicode.ttf",
    "C:/Windows/Fonts/arial.ttf",
]
_FONT_PATH = next((p for p in _FONT_CANDIDATES if os.path.exists(p)), None)


# THE TYPEFACE (0.30.0): Inter, bundled (fonts/, SIL Open Font License), so
# every product looks the same on every machine: Light for text, Medium for
# the few headings. It has no sun symbol; Draw() below draws that by hand.
_HERE = os.path.dirname(os.path.abspath(__file__))
_INTER = {False: os.path.join(_HERE, "fonts", "Inter-Light.otf"),
          True: os.path.join(_HERE, "fonts", "Inter-Medium.otf")}
# A second voice, used sparingly (0.30.0, the design pass): Lora Italic
# (bundled, SIL Open Font License) for the few lines that are read as a
# sentence -- a deck under a title, a caption, the name of a thing on the
# picture. Labels, figures and everything technical stay in Inter.
_SERIF = os.path.join(_HERE, "fonts", "Lora-Italic.ttf")
_FONTS = {}


def _serif(size):
    size = max(6, int(round(size)))
    key = (size, "serif")
    if key not in _FONTS:
        try:
            _FONTS[key] = ImageFont.truetype(_SERIF, size)
        except Exception:
            _FONTS[key] = _font(size)
    return _FONTS[key]


def _font(size, bold=False):
    size = max(6, int(round(size)))
    key = (size, bool(bold))
    if key in _FONTS:
        return _FONTS[key]
    f = None
    for p in ([_INTER[bool(bold)]] + ([_FONT_PATH] if _FONT_PATH else [])):
        try:
            if p and os.path.exists(p):
                f = ImageFont.truetype(p, size)
                break
        except Exception:
            continue
    if f is None:
        try:
            f = ImageFont.load_default(size=size)
        except Exception:
            f = ImageFont.load_default()
    _FONTS[key] = f
    return f


SUN = "\u2609"                                  # the unit symbol; drawn by Draw() where the font lacks it
INK = (232, 232, 232)                           # headings
TEXT = (205, 205, 205)                          # labels
MUTED = (150, 150, 150)                         # secondary lines


class _SunDraw(ImageDraw.ImageDraw):
    """ImageDraw whose text() draws the sun symbol (a circle with a dot) where
    the typeface has none; everything else is the plain ImageDraw."""

    def text(self, xy, text, fill=None, font=None, anchor=None, *args, **kw):
        t = str(text)
        if SUN not in t or "\n" in t or font is None:
            return super().text(xy, t.replace(SUN, "sun") if "\n" in t else t, fill, font, anchor, *args, **kw)
        anchor = anchor or "la"
        h_a, v_a = anchor[0], anchor[1]
        parts = t.split(SUN)
        sz = getattr(font, "size", 12)
        gw = 0.62 * sz                               # the symbol's advance
        widths = [self.textlength(p_, font=font) for p_ in parts]
        total = sum(widths) + gw * (len(parts) - 1)
        x, y = xy
        x0 = x - (total / 2.0 if h_a == "m" else (total if h_a == "r" else 0.0))
        # the baseline for the vertical anchor asked for
        try:
            b1 = font.getbbox("Hg", anchor="l" + v_a)
            b2 = font.getbbox("Hg", anchor="ls")
            yb = y + (b1[1] - b2[1])
        except Exception:
            yb = y + 0.8 * sz
        cx = x0
        for i, p_ in enumerate(parts):
            if p_:
                self._plain((cx, yb), p_, fill, font, "ls", *args, **kw)
            cx += widths[i]
            if i < len(parts) - 1:
                r = 0.24 * sz
                ccx, ccy = cx + gw / 2.0, yb - 0.36 * sz
                lw = max(1, int(round(sz / 14.0)))
                self.ellipse([ccx - r, ccy - r, ccx + r, ccy + r], outline=fill, width=lw)
                rd = max(0.8, 0.055 * sz)
                self.ellipse([ccx - rd, ccy - rd, ccx + rd, ccy + rd], fill=fill)
                cx += gw

    def _plain(self, xy, text, fill, font, anchor, *args, **kw):
        ImageDraw.ImageDraw.text(self, xy, text, fill, font, anchor, *args, **kw)

    @staticmethod
    def _baseline(y, font, v_a):
        """The baseline's y for a text anchored at y with vertical anchor v_a."""
        try:
            b1 = font.getbbox("Hg", anchor="l" + v_a)
            b2 = font.getbbox("Hg", anchor="ls")
            return y + (b1[1] - b2[1])
        except Exception:
            return y + 0.8 * getattr(font, "size", 12)

    def tracked_w(self, text, font, track=0.14):
        t = str(text)
        sp = track * getattr(font, "size", 12)
        return sum(self.textlength(ch, font=font) for ch in t) + sp * max(0, len(t) - 1)

    def tracked(self, xy, text, fill=None, font=None, anchor="la", track=0.14, **kw):
        """Letter-spaced text (labels in capitals): drawn letter by letter,
        `track` em between the letters. Returns its width."""
        t = str(text)
        sp = track * getattr(font, "size", 12)
        w = self.tracked_w(t, font, track)
        x, y = xy
        x0 = x - (w / 2.0 if anchor[0] == "m" else (w if anchor[0] == "r" else 0.0))
        yb = self._baseline(y, font, anchor[1])
        self._tracked_rec(x0, yb, t, fill, font, sp)
        cx = x0
        for ch in t:
            ImageDraw.ImageDraw.text(self, (cx, yb), ch, fill, font, "ls", **kw)
            cx += self.textlength(ch, font=font) + sp
        return w

    def _tracked_rec(self, x0, yb, t, fill, font, sp):
        pass


# ---------------------------------------------------------------- layers and vector lettering
# With "Files: picture + layers" a product draws its lettering not on the
# picture but on a transparent layer of the same size, and every line,
# circle, polygon and text is also written down as SVG. The picture stays
# clean (written as its own file), the layer is a transparent PNG, and the
# SVG has the picture linked under the lettering as vectors (live text in
# the bundled Inter), to work on in a vector drawing app.
import threading
_CAP = threading.local()


def _svg_col(c):
    if c is None:
        return "none", 1.0
    if isinstance(c, int):
        c = (c, c, c)
    c = tuple(int(x) for x in c)
    a = (c[3] / 255.0) if len(c) > 3 else 1.0
    return f"#{c[0]:02x}{c[1]:02x}{c[2]:02x}", a


def _esc(t):
    return str(t).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


class _RecDraw(_SunDraw):
    """_SunDraw that also writes every drawing call as an SVG element."""

    def __init__(self, im, rec, mode=None):
        super().__init__(im, mode)
        self._rec = rec

    def _stroke(self, col, w):
        h, a = _svg_col(col)
        return f'stroke="{h}" stroke-opacity="{a:.3f}" stroke-width="{max(1, w or 1)}"'

    def line(self, xy, fill=None, width=0, joint=None):
        pts = [tuple(p) for p in xy] if xy and not isinstance(xy[0], (int, float)) else \
            [(xy[i], xy[i + 1]) for i in range(0, len(xy), 2)]
        self._rec.append(f'<polyline points="{" ".join(f"{x:.2f},{y:.2f}" for x, y in pts)}" fill="none" '
                         f'{self._stroke(fill, width)} stroke-linecap="butt"/>')
        return super().line(xy, fill=fill, width=width, joint=joint)

    def ellipse(self, xy, fill=None, outline=None, width=1):
        x0, y0, x1, y1 = [float(v) for v in (list(xy[0]) + list(xy[1]) if len(xy) == 2 else xy)]
        fh, fa = _svg_col(fill)
        st = self._stroke(outline, width) if outline is not None else 'stroke="none"'
        self._rec.append(f'<ellipse cx="{(x0 + x1) / 2:.2f}" cy="{(y0 + y1) / 2:.2f}" rx="{(x1 - x0) / 2:.2f}" '
                         f'ry="{(y1 - y0) / 2:.2f}" fill="{fh}" fill-opacity="{fa:.3f}" {st}/>')
        return super().ellipse(xy, fill=fill, outline=outline, width=width)

    def polygon(self, xy, fill=None, outline=None, width=1):
        pts = [tuple(p) for p in xy]
        fh, fa = _svg_col(fill)
        st = self._stroke(outline, width) if outline is not None else 'stroke="none"'
        self._rec.append(f'<polygon points="{" ".join(f"{x:.2f},{y:.2f}" for x, y in pts)}" fill="{fh}" '
                         f'fill-opacity="{fa:.3f}" {st}/>')
        return super().polygon(xy, fill=fill, outline=outline, width=width)

    def rectangle(self, xy, fill=None, outline=None, width=1):
        x0, y0, x1, y1 = [float(v) for v in (list(xy[0]) + list(xy[1]) if len(xy) == 2 else xy)]
        fh, fa = _svg_col(fill)
        st = self._stroke(outline, width) if outline is not None else 'stroke="none"'
        self._rec.append(f'<rect x="{x0:.2f}" y="{y0:.2f}" width="{x1 - x0:.2f}" height="{y1 - y0:.2f}" '
                         f'fill="{fh}" fill-opacity="{fa:.3f}" {st}/>')
        return super().rectangle(xy, fill=fill, outline=outline, width=width)

    @staticmethod
    def _face(font):
        try:
            fam, sty = font.getname()
        except Exception:
            fam, sty = "Inter", "Light"
        if fam == "Lora":
            return 'font-family="Lora, Georgia, serif" font-style="italic" font-weight="400"'
        return ('font-family="Inter, Helvetica, Arial, sans-serif" font-weight="'
                + ("500" if "Medium" in str(sty) else "300") + '"')

    def _rec_text(self, x, yb, t, fill, font, ta="start", extra=""):
        # the baseline is worked out here (alphabetic, the one every vector
        # app agrees on) rather than left to dominant-baseline
        h, a = _svg_col(fill if fill is not None else (255, 255, 255))
        self._rec.append(f'<text x="{x:.2f}" y="{yb:.2f}" {self._face(font)} font-size="{getattr(font, "size", 12)}" '
                         f'fill="{h}" fill-opacity="{a:.3f}" text-anchor="{ta}"{extra}>{_esc(t)}</text>')

    def _plain(self, xy, text, fill, font, anchor, *args, **kw):
        # a part of a text with the sun symbol in it: anchored "ls"
        self._rec_text(xy[0], xy[1], text, fill, font)
        ImageDraw.ImageDraw.text(self, xy, text, fill, font, anchor, *args, **kw)

    def _tracked_rec(self, x0, yb, t, fill, font, sp):
        self._rec_text(x0, yb, t, fill, font, extra=f' letter-spacing="{sp:.2f}"')

    def text(self, xy, text, fill=None, font=None, anchor=None, *args, **kw):
        t = str(text)
        if SUN in t and "\n" not in t and font is not None:
            # the parts and the symbol's circles are written down one by one
            return _SunDraw.text(self, xy, t, fill, font, anchor, *args, **kw)
        if t.strip() and font is not None:
            size = getattr(font, "size", 12)
            an = anchor or "la"
            ta = {"l": "start", "m": "middle", "r": "end"}.get(an[0], "start")
            x, y = xy
            yb = self._baseline(y, font, an[1])
            for i, ln in enumerate(t.split("\n")):
                self._rec_text(x, yb + i * size * 1.2, ln, fill, font, ta)
        return super().text(xy, text, fill, font, anchor, *args, **kw)


def Draw(im, mode=None):
    cap = getattr(_CAP, "cap", None)
    if cap is None:
        return _SunDraw(im, mode)
    if im.mode == "RGBA" and cap.get("own_layer"):
        return _RecDraw(im, cap["rec"], mode)
    # the lettering goes on a transparent layer the size of this picture
    key = id(im)
    if key not in cap["layers"]:
        cap["layers"][key] = (im, Image.new("RGBA", im.size, (0, 0, 0, 0)))
    return _RecDraw(cap["layers"][key][1], cap["rec"], "RGBA")


def _cap_follow(old, new, how):
    """A captured picture replaced by `new` (padded, cropped): its layer
    follows by `how(layer) -> layer`, and the SVG coordinates by the offset
    the caller adds."""
    cap = getattr(_CAP, "cap", None)
    if cap is None or id(old) not in cap["layers"]:
        return
    _o, lay = cap["layers"].pop(id(old))
    cap["layers"][id(new)] = (new, how(lay))


def _svg_doc(W, H, href, rec, dx=0.0, dy=0.0):
    body = "\n".join(rec)
    g = f'<g transform="translate({dx:.2f},{dy:.2f})">' if (dx or dy) else "<g>"
    return (f'<?xml version="1.0" encoding="UTF-8"?>\n'
            f'<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" '
            f'width="{W}" height="{H}" viewBox="0 0 {W} {H}">\n'
            f'<!-- the picture is linked (keep it next to this file); the lettering is live text in Inter '
            f'(SIL Open Font License, in the app\'s fonts folder) -->\n'
            f'<image x="0" y="0" width="{W}" height="{H}" xlink:href="{_esc(href)}" href="{_esc(href)}"/>\n'
            f'{g}\n{body}\n</g>\n</svg>\n')


def _unit(ctx):
    return ("R" + SUN) if ctx["unit"] == "sun" else "R moon"


def _to8(a):
    a = np.clip(np.asarray(a, np.float32), 0, 1)
    if a.ndim == 2:
        a = np.repeat(a[:, :, None], 3, 2)
    elif a.shape[2] == 1:
        a = np.repeat(a, 3, 2)
    return (a * 255 + 0.5).astype(np.uint8)


def _sample(src, Y, X, order=1):
    """src (H,W,3) float, sampled at rows Y / cols X (arrays)."""
    return np.stack([ndimage.map_coordinates(src[..., c], [Y, X], order=order, mode="nearest")
                     for c in range(src.shape[2])], -1)


def _text_w(d, t, f):
    try:
        return d.textlength(t, font=f)
    except Exception:
        return len(t) * f.size * 0.55


def _dashed(d, p0, p1, col, w, dash, gap):
    x0, y0 = p0
    x1, y1 = p1
    L = math.hypot(x1 - x0, y1 - y0)
    n = 0.0
    while n < L:
        a = n / L
        b = min(1.0, (n + dash) / L)
        d.line([(x0 + (x1 - x0) * a, y0 + (y1 - y0) * a), (x0 + (x1 - x0) * b, y0 + (y1 - y0) * b)],
               fill=col, width=w)
        n += dash + gap


# ---------------------------------------------------------------- the context
def angle_kw(opts):
    """The sky angles a product was asked with (north, north by hand, the
    fine turn), as context() keywords: carried in opts["_angles"] to the
    products that build their own context (film, sound, corona in motion),
    so north entered by hand replaces the star match there too."""
    a = (opts or {}).get("_angles") or []
    try:
        north = "celestial" if (len(a) > 0 and a[0] == "celestial") else "solar"
        man = float(a[1]) if len(a) > 1 and a[1] is not None else None
        off = float(a[2]) if len(a) > 2 and a[2] else 0.0
    except Exception:
        return {}
    return dict(north=north, manual_cw=man, offset_deg=off)


def context(wd, folder, geo, scale=1.0, north="solar", manual_cw=None, offset_deg=0.0):
    """What the drawings need, at the picture's scale: centre, Moon radius, the
    unit radius in pixels, the sky orientation, the star match, the site."""
    o = _pm.orientation(wd, folder, north=north, manual_cw=manual_cw, offset_deg=offset_deg)
    s = float(scale)
    cy, cx = float(geo["cy"]) / s, float(geo["cx"]) / s
    R = float(geo["R"]) / s
    Rm = float(geo.get("Rmask", geo["R"])) / s
    Rs = (float(o["px_per_unit"]) / s) if o.get("px_per_unit") else R
    try:
        sj = json.load(open(os.path.join(wd, "stars.json")))
    except Exception:
        sj = {}
    cat = sj.get("catalogue") or {}
    solved = cat.get("status") == "solved"
    stars = []
    for st in sj.get("stars", []) or []:
        c = st.get("catalogue")
        if st.get("class") == "star" and c:
            stars.append(dict(x=float(st["x"]) / s, y=float(st["y"]) / s, mag=c.get("mag"),
                              name=(c.get("name") or c.get("label") or "")))
    site = {}
    try:
        from . import observer as _ob
        site = _ob.load(folder) or {}
    except Exception:
        site = {}
    JD = _pm._time_jd(wd, folder)
    date = _date_of(JD, site)
    o_sol = _pm.orientation(wd, folder, north="solar", manual_cw=manual_cw, offset_deg=offset_deg)
    o_cel = _pm.orientation(wd, folder, north="celestial", manual_cw=manual_cw, offset_deg=offset_deg)
    north_cel = o_cel["zero_cw"] if o_cel["origin"] == "celestial" else None
    if north_cel is None and manual_cw is not None and JD is not None:
        # north by hand stands in for the star match everywhere: celestial
        # north follows from it through the P angle of the date, as it does
        # the other way round with stars
        try:
            north_cel = (o_sol["zero_cw"] + o_sol.get("sgn", 1) * _pm.p_angle(JD)) % 360.0
        except Exception:
            north_cel = None
    return dict(wd=wd, cy=cy, cx=cx, R=R, Rm=Rm, Rs=Rs, o=o, unit=o["unit"], solved=solved, cat=cat,
                stars=stars, site=site, JD=JD, date=date, scale=s, turned=0.0,
                north_solar_cw=(o_sol["zero_cw"] if o_sol["origin"] == "solar" else None),
                north_cel_cw=north_cel,
                km_px=(Rs / KM_PER_RSUN) if o["unit"] == "sun" else None)


def _date_of(JD, site):
    if site.get("date"):
        return str(site["date"])
    if JD is None:
        return ""
    # Meeus ch. 7: JD to calendar date
    Z = int(JD + 0.5)
    F = JD + 0.5 - Z
    if Z < 2299161:
        A = Z
    else:
        al = int((Z - 1867216.25) / 36524.25)
        A = Z + 1 + al - al // 4
    B = A + 1524
    C = int((B - 122.1) / 365.25)
    D = int(365.25 * C)
    E = int((B - D) / 30.6001)
    day = B - D - int(30.6001 * E)
    month = E - 1 if E < 14 else E - 13
    year = C - 4716 if month > 2 else C - 4715
    return f"{year:04d}-{month:02d}-{day:02d}"


def _nice_date(iso):
    try:
        y, m, d = (int(x) for x in str(iso)[:10].split("-"))
        return f"{d} {['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'][m - 1]} {y}"
    except Exception:
        return str(iso)


def _dir(ctx, th):
    """Position angle (deg, from the angle origin, east-positive on the sky)
    -> unit vector in the picture (x right, y down)."""
    phi = math.radians(ctx["o"]["zero_cw"] - ctx["o"]["sgn"] * th)
    return math.sin(phi), -math.cos(phi)


def _pt(ctx, th, r):
    ux, uy = _dir(ctx, th)
    return ctx["cx"] + ux * r * ctx["Rs"], ctx["cy"] + uy * r * ctx["Rs"]


def _eclipse_kind(ctx):
    try:
        from . import eclipses as _ec
        d = ctx.get("date") or ""
        E = _ec.ECLIPSES.get(d)
        if E is not None:
            nm = str(E.get("name", "")).lower()
            if "hybrid" in nm:
                return "Hybrid solar eclipse"
            if "annular" in nm:
                return "Annular solar eclipse"
            return "Total solar eclipse"
    except Exception:
        pass
    return "Solar eclipse"


def _mid_totality_local(ctx):
    """Mid-totality at the site from the eclipse elements, as local clock
    time (the site's UTC offset), or None when the site, the date or the
    offset is not known or the site saw no totality."""
    site = ctx.get("site") or {}
    try:
        from . import eclipses as _ec
        if site.get("lat") is None or not ctx.get("date") or site.get("utc_offset") is None:
            return None
        if ctx["date"] not in _ec.ECLIPSES:
            return None
        loc = _ec.local(ctx["date"], site["lat"], site["lon"], site.get("height_m") or 0.0)
        if not loc.get("total"):
            return None
        h = (float(loc["mid_ut"]) + float(site["utc_offset"])) % 24.0
        hh = int(h)
        mm = int((h - hh) * 60)
        ss = int(round(((h - hh) * 60 - mm) * 60))
        if ss == 60:
            ss, mm = 0, mm + 1
        if mm == 60:
            mm, hh = 0, (hh + 1) % 24
        return f"{hh:02d}:{mm:02d}:{ss:02d}"
    except Exception:
        return None


def _facts(ctx):
    """The lines of the facts block, from the data only."""
    lines = []
    head = _eclipse_kind(ctx)
    if ctx["date"]:
        head += " \u00b7 " + _nice_date(ctx["date"])
    lines.append(head)
    site, cat = ctx["site"], ctx["cat"]
    l2 = []
    mid = _mid_totality_local(ctx)
    if mid:
        l2.append(f"mid-totality {mid} local")
    elif site.get("time"):
        l2.append(f"{site['time']} local")
    if ctx["solved"] and cat.get("sun_alt") is not None:
        l2.append(f"Sun {float(cat['sun_alt']):.1f}\u00b0 high")
    if l2:
        lines.append("  \u00b7  ".join(l2))
    if site.get("lat") is not None and site.get("lon") is not None:
        la, lo = float(site["lat"]), float(site["lon"])
        t = f"{abs(la):.4f}\u00b0 {'N' if la >= 0 else 'S'}  {abs(lo):.4f}\u00b0 {'E' if lo >= 0 else 'W'}"
        if site.get("height_m") is not None:
            t += f"  \u00b7  {float(site['height_m']):.0f} m"
        lines.append(t)
    if ctx["solved"]:
        t = f"north from {int(cat.get('matches', 0))} catalogue stars"
        if ctx["o"]["origin"] == "solar":
            t += f"  \u00b7  P {ctx['o']['p']:+.2f}\u00b0"
        lines.append(t)
    elif ctx["o"]["origin"] == "solar":
        lines.append("solar north entered by hand")
    else:
        lines.append("no star match: angles from the top of the picture")
    return lines


# ---------------------------------------------------------------- the sheet's voice
# One hierarchy for every lettered product (0.30.0, the design pass):
#   kicker   what it is, small capitals letter-spaced, the one accent colour
#   title    the date, large and light
#   deck     one sentence of facts, serif italic
#   labels   capitals letter-spaced and muted over their value
#   callout  a hairline from the thing itself to a figure and a few words
#   sources  what the sheet was made from, small, at the foot
# Everything said is from the run's own data; a figure that is not known is
# left out, never rounded up into a claim.
MONTHS = ["January", "February", "March", "April", "May", "June", "July", "August", "September", "October",
          "November", "December"]


def _tick_lab(ctx, a):
    """The label of an angle tick: N, E, S, W at the quarters when north is
    known (star match or by hand), else the bare angle from the picture's top."""
    if ctx["o"]["origin"] in ("solar", "celestial"):
        return {0: "N", 90: "E", 180: "S", 270: "W"}.get(int(a) % 360, f"{int(a)}\u00b0")
    return f"{int(a)}\u00b0"


def _origin_words(ctx):
    o = ctx["o"]["origin"]
    return ("position angle from solar north" if o == "solar" else
            ("position angle from celestial north" if o == "celestial" else "angle from the picture\u2019s top (no star match)"))


def _long_date(iso):
    try:
        y, m, d = (int(x) for x in str(iso)[:10].split("-"))
        return f"{d} {MONTHS[m - 1]} {y}"
    except Exception:
        return str(iso or "")


def _totality_s(ctx):
    """Seconds of totality at the site from the eclipse elements, or None."""
    try:
        from . import eclipses as _ec
        site = ctx["site"]
        if site.get("lat") is not None and ctx.get("date"):
            loc = _ec.local(ctx["date"], site["lat"], site["lon"], site.get("height_m") or 0.0)
            if loc.get("total"):
                return float(loc["duration_s"])
    except Exception:
        pass
    return None


def _deck(ctx):
    """One line of facts under a title, from the data only."""
    parts = []
    mid = _mid_totality_local(ctx)
    if mid:
        parts.append(f"mid-totality at {mid} local time")
    dur = _totality_s(ctx)
    if dur:
        parts.append(f"{dur:.0f} seconds of totality")
    if ctx["solved"] and ctx["cat"].get("sun_alt") is not None:
        parts.append(f"the Sun {float(ctx['cat']['sun_alt']):.0f}° above the horizon")
    line = ", ".join(parts)
    return (line[0].upper() + line[1:]) if line else ""


def _sources(ctx):
    """What the picture was made from, as short parts: frames and exposures
    from the run's report, how north was found."""
    parts = []
    try:
        rep = json.load(open(os.path.join(ctx["wd"], "report.json")))
        n = int(rep.get("n_files") or 0)
        secs = [float(t["sec"]) for t in (rep.get("tiers") or []) if t.get("sec")]
        if n:
            parts.append(f"{n} frames")
        if len(secs) > 1:
            parts.append(f"{len(secs)} exposures, {_fmt_s(min(secs))} to {_fmt_s(max(secs))}")
    except Exception:
        pass
    if ctx["solved"]:
        parts.append(f"north from {int(ctx['cat'].get('matches', 0))} catalogue stars")
    elif ctx["o"]["origin"] == "solar":
        parts.append("solar north entered by hand")
    return parts


def _wrap(text, font, wmax, measure):
    """`text` in as few lines as fit `wmax`, broken at spaces, the lines kept
    about even (a long line over a short word reads badly)."""
    if wmax is None or measure(text, font) <= wmax:
        return [text]
    words = text.split(" ")
    for n in range(2, 6):
        # the even split into n lines: each break at the space nearest its share
        tot = measure(text, font)
        lines, cur = [], []
        for w_ in words:
            cand = " ".join(cur + [w_])
            if cur and measure(cand, font) > tot / n * 1.08 and len(lines) < n - 1:
                lines.append(" ".join(cur))
                cur = [w_]
            else:
                cur.append(w_)
        lines.append(" ".join(cur))
        if max(measure(l_, font) for l_ in lines) <= wmax:
            return lines
    return lines


def _deck_lines(ctx, fs, wmax):
    line = _deck(ctx)
    if not line:
        return []
    m = Draw(Image.new("L", (4, 4)))
    return _wrap(line, _serif(fs * 1.3), wmax, lambda t, f: m.textlength(t, font=f))


def _plate_title(d, x, y, fs, ctx, align="l", deck=True, big=3.4, wmax=None):
    """Kicker, date and (optionally) the deck at (x, y); returns the y below.
    `wmax`: the width the lettering must keep within (the deck wraps, the
    date gets smaller)."""
    an = align
    d.tracked((x, y), _eclipse_kind(ctx).upper(), fill=GOLD, font=_font(fs * 0.95, True), anchor=an + "a",
              track=0.22)
    y += fs * 2.0
    if ctx.get("date"):
        fd = _font(fs * big)
        t_ = _long_date(ctx["date"])
        if wmax and d.textlength(t_, font=fd) > wmax:
            fd = _font(fs * big * wmax / d.textlength(t_, font=fd))
        d.text((x, y), t_, fill=INK, font=fd, anchor=an + "a")
        y += fs * big * 1.32
    if deck:
        f_dk = _serif(fs * 1.3)
        lines = _wrap(_deck(ctx), f_dk, wmax, lambda t, f: d.textlength(t, font=f)) if _deck(ctx) else []
        for i, line in enumerate(lines):
            d.text((x, y + i * fs * 1.7), line, fill=TEXT, font=f_dk, anchor=an + "a")
        if lines:
            y += fs * 2.3 + (len(lines) - 1) * fs * 1.7
    return y


def _source_line(d, xc, y, wmax, parts, fs, col=(120, 120, 120)):
    """The small line at the foot: the parts in letter-spaced capitals,
    centred at xc, on as many lines as the width asks for. Returns the y below."""
    f = _font(fs, True)
    sep = "  \u00b7  "
    lines, cur = [], []
    for p_ in parts:
        t = sep.join(cur + [p_]).upper()
        if cur and d.tracked_w(t, f, 0.16) > wmax:
            lines.append(sep.join(cur).upper())
            cur = [p_]
        else:
            cur.append(p_)
    if cur:
        lines.append(sep.join(cur).upper())
    for i, t in enumerate(lines):
        d.tracked((xc, y + i * fs * 1.9), t, fill=col, font=f, anchor="ma", track=0.16)
    return y + len(lines) * fs * 1.9


def _round_km(km):
    """A prominence's height as it can honestly be said: its top is found
    to a few pixels, so to the nearest 5 000 km (1 000 km under 10 000)."""
    return float(round(km / 5000.0) * 5000) if km >= 10000 else float(round(km / 1000.0) * 1000)


def _earths(km):
    n = km / EARTH_DIAM_KM
    return f"{n:.0f} Earths" if n >= 2.5 else (f"{n:.1f} Earths" if n >= 0.75 else "smaller than Earth")


def _callout(d, p0, p1, head, sub, fs, lw=1, halo=True):
    """A hairline from the thing (a small open ring at p0) to p1, and there a
    figure (head) over a few words (sub, serif italic). The text stands on the
    far side of p1 from p0."""
    x0, y0 = p0
    x1, y1 = p1
    r = max(2.0, fs * 0.16)
    L = math.hypot(x1 - x0, y1 - y0) or 1.0
    ux, uy = (x1 - x0) / L, (y1 - y0) / L
    d.line([(x0 + ux * r, y0 + uy * r), (x1, y1)], fill=(255, 255, 255, 150), width=lw)
    d.ellipse([x0 - r, y0 - r, x0 + r, y0 + r], outline=(255, 255, 255, 200), width=lw)
    left = ux < -0.25
    centred = abs(ux) <= 0.25
    an = "r" if left else ("m" if centred else "l")
    tx = x1 + (0 if centred else (-fs * 0.45 if left else fs * 0.45))
    up = centred and uy < 0
    ty = y1 - (fs * 2.9 if up else (-fs * 0.5 if centred else fs * 0.75))
    kw = dict(stroke_width=max(1, int(fs * 0.10)), stroke_fill=(0, 0, 0, 120)) if halo else {}
    d.text((tx, ty), head, fill=(245, 245, 245, 255), font=_font(fs * 1.05, True), anchor=an + "a", **kw)
    if sub:
        d.text((tx, ty + fs * 1.4), sub, fill=(215, 215, 215, 235), font=_serif(fs * 0.98), anchor=an + "a", **kw)


def _prom_callouts(d, ctx, proms, fs, lw, pt, gap=0.13, maxn=6):
    """Callouts for the prominences found: from each one's top outward. `pt`
    maps (position angle, radius in units) to the drawing's pixels."""
    Rm = ctx["Rm"] / ctx["Rs"]
    sel = sorted(sorted(proms, key=lambda p: -p["score"])[:maxn], key=lambda p: p["pa"])
    last_pa, last_long = None, False
    for p in sel:
        close = last_pa is not None and abs(p["pa"] - last_pa) < 16.0
        long_ = close and not last_long
        r0 = Rm + p["height"] + 0.012
        r1 = r0 + gap + (0.17 if long_ else 0.0)
        if p["height_km"] is not None and p["height_km"] >= 1000:
            hk = _round_km(p["height_km"])
            head = f"\u2248 {hk:,.0f} km".replace(",", " ")
            sub = f"prominence \u00b7 {_earths(hk)}"
        else:
            head = "Prominence"
            sub = f"at {p['pa']:.0f}°"
        _callout(d, pt(p["pa"], r0), pt(p["pa"], r1), head, sub, fs, lw)
        last_pa, last_long = p["pa"], long_


# ---------------------------------------------------------------- turn and frame



ROTATE_VALUES = [["none", "as captured"], ["solar", "solar north up"], ["celestial", "celestial north up"],
                 ["zenith", "zenith up"], ["cw90", "90\u00b0 clockwise"], ["ccw90", "90\u00b0 counter-clockwise"],
                 ["r180", "180\u00b0"]]


def rotate_opt(title):
    return dict(id="rotate", label="Picture up", type="select", default="none", values=ROTATE_VALUES, title=title)


def turn_angle(ctx, turn):
    """Degrees (clockwise from up in the picture) of the direction that is to
    stand up, or None when it is not known."""
    if turn == "solar":
        return ctx.get("north_solar_cw")
    if turn == "celestial":
        return ctx.get("north_cel_cw")
    if turn == "zenith":
        z = ctx["cat"].get("zenith_deg_cw_from_up") if ctx["solved"] else None
        return float(z) if z is not None else None
    return {"cw90": 270.0, "ccw90": 90.0, "r180": 180.0}.get(turn, 0.0)


def reframe(rgb, ctx, turn="none", frame="as", fy=0.5):
    """The picture turned so that `turn` stands up and cropped to `frame`
    around the Sun; the context follows (centre, angle origin, stars, the
    compass angles). Returns (rgb, ctx). Without a turn and with the frame as
    captured, the picture is returned as it is."""
    th = turn_angle(ctx, turn)
    th = 0.0 if th is None else float(th) % 360.0
    H, W = rgb.shape[:2]
    if th == 0.0 and _aspect(frame) is None and fy == 0.5:
        return rgb, ctx
    aw, ah = _aspect(frame) or (W, H)
    cx, cy = ctx["cx"], ctx["cy"]
    # the largest centred rectangle of this shape whose corners, turned back
    # into the picture, stay inside it
    t = math.radians(th)
    f = None
    # the Sun at `fy` of the height from the top (0.5: the middle)
    for u, v in ((-aw, -2 * ah * fy), (aw, -2 * ah * fy), (-aw, 2 * ah * (1 - fy)), (aw, 2 * ah * (1 - fy))):
        # output offset (u, v) -> picture offset: a direction at angle a (cw
        # from up) in the output was at a + th in the picture
        a = math.atan2(u, -v) + t
        r = math.hypot(u, v)
        dx_, dy_ = math.sin(a) * r, -math.cos(a) * r
        for d, room in ((dx_, (W - 1 - cx) if dx_ > 0 else cx), (dy_, (H - 1 - cy) if dy_ > 0 else cy)):
            if abs(d) > 1e-9:
                f = room / abs(d) if f is None else min(f, room / abs(d))
    f = f or 1.0
    ow, oh = int(2 * f * aw), int(2 * f * ah)
    ow, oh = max(2, ow - ow % 2), max(2, oh - oh % 2)
    yy, xx = np.mgrid[0:oh, 0:ow].astype(np.float32)
    u, v = xx - ow / 2.0, yy - oh * fy
    if th == 0.0:
        X, Y = cx + u, cy + v
    else:
        ca, sa = math.cos(t), math.sin(t)
        # rotate the output offset by +th (cw) into the picture
        X = cx + u * ca + (-v) * sa          # x' = u cos + (-v) sin  (y up positive = -v)
        Y = cy - (-(u) * sa + (-v) * ca)
    out = _sample(np.asarray(rgb, np.float32), Y, X, order=1)
    c2 = dict(ctx)
    c2["cx"], c2["cy"] = ow / 2.0, oh * fy
    c2["turned"] = th
    o = dict(ctx["o"])
    o["zero_cw"] = (o["zero_cw"] - th) % 360.0
    c2["o"] = o
    for k in ("north_solar_cw", "north_cel_cw"):
        if ctx.get(k) is not None:
            c2[k] = (ctx[k] - th) % 360.0
    cat = dict(ctx["cat"])
    for k in ("north_deg_cw_from_up", "east_deg_cw_from_up", "zenith_deg_cw_from_up"):
        if cat.get(k) is not None:
            cat[k] = (float(cat[k]) - th) % 360.0
    c2["cat"] = cat
    stars = []
    for st in ctx["stars"]:
        dx_, dy_ = st["x"] - cx, st["y"] - cy
        a = math.atan2(dx_, -dy_) - t
        r = math.hypot(dx_, dy_)
        x, y = c2["cx"] + math.sin(a) * r, c2["cy"] - math.cos(a) * r
        if 0 <= x < ow and 0 <= y < oh:
            stars.append(dict(st, x=x, y=y))
    c2["stars"] = stars
    return out.astype(np.float32), c2


# ---------------------------------------------------------------- 1 overlay
def overlay(rgb, ctx, v):
    """Returns (composite PIL RGB, layer PIL RGBA)."""
    H, W = rgb.shape[:2]
    cx, cy, Rs = ctx["cx"], ctx["cy"], ctx["Rs"]
    k = float(v["size"])
    # sizes follow the picture's area, not its width: a portrait frame is
    # looked at as large as a landscape one (W / 165 on a 3:2 picture)
    G = math.sqrt(W * H) * 1.2333
    fs = max(7.0, G / 165.0) * k           # body lettering
    if v.get("_social"):
        fs = SOCIAL_FS * k                  # read on a phone: the size in px, not a share of the picture
    lw = max(1, int(round(G / 2000.0)))     # line width
    # the safe zone (px of this canvas): the corner lettering stays inside it
    s_top, s_right, s_bot, s_left, s_cw, s_ch = v.get("_safe_px") or (0, 0, 0, 0, 0, 0)
    halo = dict(stroke_width=max(1, int(fs * 0.10)), stroke_fill=(0, 0, 0, 120))
    ov = Image.new("RGBA", (W, H), (0, 0, 0, 0))
    d = Draw(ov)
    white = (255, 255, 255)
    unit = _unit(ctx)
    # the largest ring that fits the frame (in a safe zone: inside it, with
    # room for the ring's labels; the reel's column of buttons is left out of
    # this, or the ring would shrink to a third)
    rmax = min(cx, W - 1 - cx, cy, H - 1 - cy) / Rs
    if s_top or s_left:
        rmax = min(rmax, (min(cx - s_left, W - s_right - cx, cy - s_top, H - s_bot - cy) - fs * 1.6) / Rs)
    if v.get("_social") and H >= W:
        rmax = min(rmax, (H - s_bot - fs * 13.0 - cy) / Rs)      # clear of the lettering at the foot
    if v["rings"]:
        for r in (1.5, 2.0, 2.5):
            if r + 0.1 > rmax:
                break
            rr = r * Rs
            for kk in range(0, 360, 2):
                a0, a1 = math.radians(kk), math.radians(kk + 1.2)
                d.line([(cx + rr * math.cos(a0), cy + rr * math.sin(a0)),
                        (cx + rr * math.cos(a1), cy + rr * math.sin(a1))], fill=white + (70,), width=lw)
            x, y = _pt(ctx, 222, r)
            d.text((x, y), f"{r:g} {unit}", fill=white + (170,), font=_font(fs), anchor="mm", **halo)
    RR = min(3.2, rmax - 0.25)
    if v["paring"] and RR > 1.3:
        d.ellipse([cx - RR * Rs, cy - RR * Rs, cx + RR * Rs, cy + RR * Rs], outline=white + (110,), width=lw)
        for th in range(0, 360, 5):
            big = th % 30 == 0
            x0, y0 = _pt(ctx, th, RR)
            x1, y1 = _pt(ctx, th, RR + (0.07 if big else 0.035))
            d.line([(x0, y0), (x1, y1)], fill=white + (200 if big else 120,), width=lw + (1 if big else 0))
            if big:
                x2, y2 = _pt(ctx, th, RR + 0.14)
                d.text((x2, y2), f"{th}\u00b0", fill=white + (210,), font=_font(fs * 1.1), anchor="mm")
    if v["axis"]:
        r_in, r_out = ctx["Rm"] / Rs + 0.03, (RR - 0.05 if RR > 1.3 else min(rmax - 0.1, 2.5))
        for th in (0, 180):
            d.line([_pt(ctx, th, r_in), _pt(ctx, th, r_out)], fill=GOLD + (230,), width=lw + 1)
        for th in (90, 270):
            _dashed(d, _pt(ctx, th, r_in), _pt(ctx, th, r_out), (255, 225, 170, 170), lw,
                    dash=G / 230.0, gap=G / 375.0)
        # north, east, south and west only when north is known (star match
        # or entered by hand); else the lines are the picture's own axes
        for th, lab in (((0, "N"), (180, "S"), (90, "E"), (270, "W"))
                        if ctx["o"]["origin"] in ("solar", "celestial") else ()):
            x, y = _pt(ctx, th, r_out - 0.27)
            d.text((x, y), lab, fill=GOLD + (255,), font=_font(fs * 2.0, True), anchor="mm")
        x, y = _pt(ctx, 10, min(2.3, r_out - 0.6))
        if ctx["o"]["origin"] == "solar":
            lab = f"solar axis  (P = {ctx['o']['p']:+.1f}\u00b0)" if ctx["solved"] else "solar axis (by hand)"
        elif ctx["o"]["origin"] == "celestial":
            lab = "celestial north"
        else:
            lab = "picture up (no star match)"
        d.text((x, y), lab, fill=GOLD + (235,), font=_font(fs * 1.05), anchor="mm")
    if v.get("proms", True):
        try:
            strip, th_s, rr_s = _limb_strip(np.asarray(rgb, np.float32), ctx, 1.35)
            _prom_callouts(d, ctx, find_prominences(None, ctx, strip, th_s, rr_s), fs, lw, lambda pa, r: _pt(ctx, pa, r))
        except Exception:
            pass
    if v["stars"]:
        rs = max(3, G / 230.0)
        f_st, f_mag = _serif(fs * 1.12), _font(fs * 0.95)
        circles = [(RR, 2)] if (v["paring"] and RR > 1.3) else []
        circles += [(r_, 1) for r_ in (1.5, 2.0, 2.5) if v["rings"] and r_ + 0.1 <= rmax]

        def crosses(x0, x1, yc):
            # how badly a drawn circle runs through the label's box: the
            # solid angle ring counts twice a dashed radius ring
            ds = [math.hypot(xx - cx, yy - cy) / Rs for xx in (x0, x1) for yy in (yc - fs * 0.7, yc + fs * 0.7)]
            return sum(w_ for c_, w_ in circles if min(ds) < c_ < max(ds))
        for st in ctx["stars"]:
            x, y = st["x"], st["y"]
            if not (0 <= x < W and 0 <= y < H):
                continue
            if (s_top or s_left) and not (s_left + fs <= x <= W - s_right - fs and s_top + fs <= y <= H - s_bot - fs):
                continue                    # under the app's lettering: its name would be too
            if v.get("_social") and (y > H - s_bot - fs * 13 or (x < fs * 14 + s_left and y < fs * 14 + s_top)):
                continue                    # a small frame: the corners hold the compass, scale and facts
            d.ellipse([x - rs, y - rs, x + rs, y + rs], outline=(140, 200, 255, 230), width=lw)
            # the name as a thing on the sky (italic), its brightness as a figure
            nm = st["name"].strip()
            mag = f"{float(st['mag']):.1f} mag" if st.get("mag") is not None else ""
            wn = _text_w(d, nm, f_st)
            wl = wn + (fs * 0.6 + _text_w(d, mag, f_mag) if mag else 0)
            c_l = crosses(x - rs * 1.6 - wl, x - rs * 1.6, y) + (10 if x - rs * 1.6 - wl < fs else 0)
            c_r = crosses(x + rs * 1.6, x + rs * 1.6 + wl, y) + (10 if x + rs * 1.6 + wl > W - fs else 0)
            on_left = c_l <= c_r
            x0 = x - rs * 1.6 - wl if on_left else x + rs * 1.6
            d.text((x0, y), nm, fill=(170, 215, 255, 240), font=f_st, anchor="lm", **halo)
            if mag:
                d.text((x0 + wn + fs * 0.6, y), mag, fill=(170, 215, 255, 200), font=f_mag, anchor="lm", **halo)
    if v["compass"] and ctx["solved"]:
        cat = ctx["cat"]
        L = G / 29.0 if not v.get("_social") else fs * 3.4
        bx, by = max(L * 2.0, s_left + L * 1.7), max(L * 2.0, s_top + L * 1.7)
        d.ellipse([bx - L * 1.3, by - L * 1.3, bx + L * 1.3, by + L * 1.3], outline=white + (55,), width=lw)

        def arrow(a_cw, col, lab, dash=False, lab_r=1.24):
            a = math.radians(a_cw)
            ux, uy = math.sin(a), -math.cos(a)
            x1, y1 = bx + ux * L, by + uy * L
            if dash:
                _dashed(d, (bx, by), (x1, y1), col, lw + 1, dash=L / 8.0, gap=L / 13.0)
            else:
                d.line([(bx, by), (x1, y1)], fill=col, width=lw + 1)
            ah = L / 8.0
            an = math.atan2(uy, ux)
            d.polygon([(x1, y1), (x1 - ah * math.cos(an - 0.4), y1 - ah * math.sin(an - 0.4)),
                       (x1 - ah * math.cos(an + 0.4), y1 - ah * math.sin(an + 0.4))], fill=col)
            d.text((bx + ux * (L * lab_r), by + uy * (L * lab_r)), lab, fill=col, font=_font(fs * 1.15, True),
                   anchor="mm")
        arrow(float(cat["north_deg_cw_from_up"]), (235, 235, 235, 255), "N")
        arrow(float(cat["east_deg_cw_from_up"]), (235, 235, 235, 255), "E")
        if cat.get("zenith_deg_cw_from_up") is not None:
            zc = float(cat["zenith_deg_cw_from_up"])
            near = min(abs(((zc - float(cat[k_])) + 180.0) % 360.0 - 180.0)
                       for k_ in ("north_deg_cw_from_up", "east_deg_cw_from_up"))
            arrow(zc, (150, 200, 150, 255), "zenith", True, 1.62 if near < 40.0 else 1.24)
        d.tracked((bx, by + L * 1.55), "CELESTIAL", fill=(170, 170, 170, 220), font=_font(fs * 0.8, True),
                  anchor="mm", track=0.2)
    sc_box = None                           # the scale block's extent, for the facts below
    if v["scale"]:
        x0, y0 = max(G / 46.0, s_left), min(H - G / 19.5, H - s_bot - fs * 2.6)
        if ctx["km_px"]:
            km_px = ctx["km_px"]
            bar_km = 500000.0
            bar = bar_km * km_px
            while bar > W * 0.28 and bar_km > 50000:
                bar_km /= 2
                bar = bar_km * km_px
            er = EARTH_DIAM_KM / 2 * km_px
            d.text((x0, y0 - fs * 3.5),
                   f"1 {unit} = 695 700 km = {Rs * ctx['scale']:.0f} px \u00b7 "
                   f"{float(ctx['cat']['arcsec_per_px']):.2f}\u2033/px", fill=(205, 205, 205, 230),
                   font=_font(fs), anchor="lm")
            if v.get("proms", True):
                d.text((x0, y0 - fs * 5.3), "Prominence heights are measured from the Moon\u2019s edge.",
                       fill=(190, 190, 190, 230), font=_serif(fs), anchor="lm")
            d.line([(x0, y0), (x0 + bar, y0)], fill=white + (230,), width=lw + 1)
            for kk in range(6):
                xx = x0 + kk * bar / 5
                d.line([(xx, y0 - fs * 0.4), (xx, y0 + fs * 0.4)], fill=white + (230,), width=lw)
            d.text((x0, y0 + fs * 0.9), "0", fill=(228, 228, 228, 230), font=_font(fs), anchor="mt")
            d.text((x0 + bar, y0 + fs * 0.9), f"{bar_km:,.0f} km".replace(",", " "), fill=(228, 228, 228, 230),
                   font=_font(fs), anchor="mt")
            ex, ey = x0 + bar + er + fs * 6.5, y0
            _wt = _text_w(d, f"1 {unit} = 695 700 km = {Rs * ctx['scale']:.0f} px \u00b7 "
                             f"{float(ctx['cat']['arcsec_per_px']):.2f}\u2033/px", _font(fs))
            if v.get("proms", True):
                _wt = max(_wt, _text_w(d, "Prominence heights are measured from the Moon\u2019s edge.", _serif(fs)))
            sc_box = (x0, y0 - fs * (6.0 if v.get("proms", True) else 4.2), max(ex + er + fs * 2.5, x0 + _wt),
                      y0 + fs * 2.4)
            d.ellipse([ex - er, ey - er, ex + er, ey + er], fill=(120, 170, 255, 255))
            d.text((ex, ey + max(er, fs * 0.5) + fs * 0.4), "Earth", fill=(170, 200, 255, 240), font=_font(fs),
                   anchor="mt")
        else:
            sc_box = (x0, y0 - fs * 1.4, x0 + fs * 22, y0 + fs * 1.4)
            d.text((x0, y0 - fs * 0.6), f"1 {unit} = {Rs * ctx['scale']:.0f} px",
                   fill=(205, 205, 205, 230), font=_font(fs), anchor="lm")
            d.text((x0, y0 + fs * 0.7), "no star match: no scale in kilometres", fill=(150, 150, 150, 220),
                   font=_font(fs * 0.9), anchor="lm")
    if v["info"]:
        # the title block, bottom right: what, when, then the facts
        lines = _facts(ctx)[1:]
        fn = _font(fs * 1.02)
        lh = fs * 1.7
        tx = W - max(G / 50.0, s_right)
        block = fs * 2.0 + fs * 2.2 * 1.32 + lh * len(lines)
        foot = min(H - G / 26.0, H - s_bot)
        if s_cw and foot > H - s_ch:
            tx = min(tx, W - s_cw)          # beside the app's column of buttons
        ty = foot - block
        if sc_box is not None and (s_top or v.get("_social")):
            # a narrow frame: the facts and the scale would meet; the facts
            # then stand above the scale
            wf = max([_text_w(d, t, fn) for t in lines] + [fs * 14]
                     + ([_text_w(d, _long_date(ctx["date"]), _font(fs * 2.2))] if ctx.get("date") else []))
            if tx - wf < sc_box[2] and ty + block > sc_box[1]:
                ty = sc_box[1] - fs * 1.2 - block
        yb_ = _plate_title(d, tx, ty, fs, ctx, align="r", deck=False, big=2.2)
        for i, t in enumerate(lines):
            d.text((tx, yb_ + i * lh), t, fill=(215, 215, 215, 235), font=fn, anchor="ra")
    comp = Image.fromarray(_to8(rgb)).convert("RGBA")
    comp.alpha_composite(ov)
    return comp.convert("RGB"), ov


# ---------------------------------------------------------------- shared unrolling
def _angle_cols(th, lab):
    """The columns (float, pixel centres) where the angle `lab` falls on a map
    whose columns carry the angles `th`: both ends when it is at the seam."""
    nt = len(th)
    step = 360.0 / nt
    d = np.abs(((th - lab) + 180.0) % 360.0 - 180.0)
    out = []
    for i in np.nonzero(d <= step * 0.51)[0]:
        if not out or i - out[-1] > 2:
            out.append(int(i))
    cols = [float(i) for i in out]
    # an angle exactly at the seam sits on the outer edge of the first or last column
    if 0 in out:
        cols[0] = 0.0
    if nt - 1 in out:
        cols[-1] = float(nt - 1)
    return cols


def _pa_of_up(ctx, up_cw):
    """The position angle that points to `up_cw` (deg clockwise from the
    picture's top)."""
    return ((ctx["o"]["zero_cw"] - up_cw) * ctx["o"]["sgn"]) % 360.0


def _turned_window(src, ctx, rot, above, below, half_w):
    """The part of the picture around the Sun with the direction `rot` (deg
    clockwise from up) turned to the top: `above` px up, `below` px down,
    `half_w` px to each side; black where the picture does not reach."""
    H, W = src.shape[:2]
    oh, ow = int(round(above + below)), int(round(2 * half_w))
    yy, xx = np.mgrid[0:oh, 0:ow].astype(np.float32)
    u, v = xx - half_w, yy - above                     # output offsets (right, down)
    t = math.radians(rot)
    ca, sa = math.cos(t), math.sin(t)
    X = ctx["cx"] + u * ca - v * sa
    Y = ctx["cy"] + u * sa + v * ca
    out = _sample(src, Y, X)
    out[~((X >= 0) & (X <= W - 1) & (Y >= 0) & (Y <= H - 1))] = 0.0
    return out


def _flat_map(src, ctx, r0, r1, centre="N", width=None, mid_deg=None):
    """The picture unrolled around the centre with the arc length kept at one
    unit radius: width = 2 pi Rs columns (or `width`), rows from r1 (top) down
    to r0. The angle falls to the right, `centre` in the middle (sky-true, as
    the polar map). Returns float32 (rows, cols, 3) and the angle per column."""
    Rs = ctx["Rs"]
    nt = int(width or round(2 * math.pi * Rs))
    nr = max(2, int(round((r1 - r0) * Rs)))
    mid = float(mid_deg) if mid_deg is not None else \
        {"N": 0.0, "E": 90.0, "S": 180.0, "W": 270.0}.get(str(centre).upper(), 0.0)
    th = (mid + 180.0 - (np.arange(nt) + 0.5) / nt * 360.0) % 360.0        # falls to the right
    rr = r1 - (np.arange(nr) + 0.5) / nr * (r1 - r0)                        # row 0 outermost
    phi = np.radians(ctx["o"]["zero_cw"] - ctx["o"]["sgn"] * th)[None, :]
    Y = ctx["cy"] - np.cos(phi) * rr[:, None] * Rs
    X = ctx["cx"] + np.sin(phi) * rr[:, None] * Rs
    H, W = src.shape[:2]
    ok = (X >= 0) & (X <= W - 1) & (Y >= 0) & (Y <= H - 1)
    out = _sample(src, Y, X)
    out[~ok] = _pm.NO_DATA
    return out.astype(np.float32), th, rr


# ---------------------------------------------------------------- 2 poster
def poster(rgb, ctx, v, rgb_map=None):
    """`rgb_map`: the same picture in another look for the map below (the
    soft look above, its standard variant below); None: the picture's."""
    src = np.asarray(rgb, np.float32)
    src_m = np.asarray(rgb_map, np.float32) if rgb_map is not None and np.shape(rgb_map) == src.shape else src
    H, W = src.shape[:2]
    cx, cy, Rs = ctx["cx"], ctx["cy"], ctx["Rs"]
    bgc = np.array(BG, np.float32) / 255.0
    rmax = min(cx, W - 1 - cx, cy, H - 1 - cy) / Rs
    top = rmax if v["height"] == "frame" else min(float(v["height"]), rmax)
    r_bot = max(0.0, ctx["Rm"] / Rs - 0.05)
    # the picture turned so that the chosen direction is up; the map then has
    # that same direction in its middle, so the strip is the picture opened
    # at the bottom (the unrolling read straight down)
    rot = turn_angle(ctx, v.get("turn", "solar"))
    rot = 0.0 if rot is None else float(rot)
    mp, th, rr = _flat_map(src_m, ctx, r_bot, top, mid_deg=_pa_of_up(ctx, rot))
    wm = mp.shape[1]
    # the round picture: as wide as the map, up to 2.5 R above the centre and
    # 2.3 R below (the fade starts in the lower corona)
    up = min(2.5, rmax + 1.5) * Rs
    fade = int(0.9 * Rs) if v["fade"] == "soft" else 1
    fs = float(v["_fs"]) if v.get("_fs") else max(7.0, wm / 150.0)
    # under the map: the angle ticks, then the title block (kicker, date,
    # deck) and one small line of what the sheet shows and was made from
    cap_h = int(fs * 17.5) if v["caption"] else int(fs * 0.8)
    # where the title block goes (a social frame may keep it left of the
    # app's buttons); a deck too long for it takes more lines
    lx0, lx1 = v.get("_letter_x") or (0.0, float(wm))
    n_dk = len(_deck_lines(ctx, fs, (lx1 - lx0) * 0.96)) if v["caption"] else 0
    if n_dk > 1:
        cap_h += int((n_dk - 1) * fs * 1.7)
    if v.get("_letter_x"):
        cap_h += int(fs * 2.8)                 # the narrower source line takes two more lines
    down = min(2.3, rmax + 1.5) * Rs
    Hn = int(round(up + down)) + mp.shape[0] - fade + cap_h
    aw, ah = _aspect(v.get("frame", "as")) or (wm, Hn)
    want_h = int(round(wm * ah / aw))
    kt = 1.0                                   # the title block's size against the lettering
    if want_h > Hn:
        # a taller frame: more of the corona above first, as far as the
        # picture reaches straight up from the Sun; then a larger title
        up += min(want_h - Hn, max(0.0, rmax * Rs - up))
        left_ = want_h - (int(round(up + down)) + mp.shape[0] - fade + cap_h)
        if v["caption"] and left_ > 0:
            kt = min(2.2, 1.0 + 0.5 * left_ / (fs * 12.9))
            cap_h = int(fs * (4.6 + 12.9 * kt) + max(0, n_dk - 1) * fs * 1.7 * kt)
    rnd = _turned_window(src, ctx, rot, up, down, wm / 2.0)
    Hc = rnd.shape[0] + mp.shape[0] - fade + cap_h
    # what is still missing to the frame's height is black: the picture's top
    # then fades into it (a cut from sky to black would read as a mistake)
    extra = max(0, want_h - Hc)
    y_off = int(extra * 0.4)
    if extra > 0:
        ft_ = int(min(1.8 * Rs, 0.4 * rnd.shape[0]))
        ramp = np.clip(np.arange(rnd.shape[0], dtype=np.float32) / max(1, ft_), 0, 1)
        rnd = rnd * (ramp * ramp * (3 - 2 * ramp))[:, None, None]
    Hp = Hc + extra
    out = np.empty((Hp, wm, 3), np.float32)
    out[:] = bgc
    out[y_off:y_off + rnd.shape[0]] = rnd
    ys = y_off + rnd.shape[0] - fade
    a = (0.5 - 0.5 * np.cos(np.pi * np.clip(np.arange(mp.shape[0]) / max(fade, 1), 0, 1)))[:, None, None]
    seg = out[ys:ys + mp.shape[0]].copy()
    seg[fade:] = bgc
    out[ys:ys + mp.shape[0]] = seg * (1 - a) + mp * a
    im = Image.fromarray(_to8(out))
    if v["caption"]:
        d = Draw(im)
        yb = ys + mp.shape[0]
        for th_lab in range(0, 360, 30):
            for col in _angle_cols(th, th_lab):
                d.line([(col, yb), (col, yb + fs * 0.5)], fill=(165, 165, 165), width=1)
                lab = _tick_lab(ctx, th_lab)
                anc = "lt" if col < fs * 1.5 else ("rt" if col > wm - fs * 1.5 else "mt")
                d.text((col, yb + fs * 0.8), lab, fill=GOLD if lab in "NESW" and len(lab) == 1 else (205, 205, 205),
                       font=_font(fs * (1.15 if len(lab) == 1 else 1.0), len(lab) == 1), anchor=anc)
        xc_, wl_ = (lx0 + lx1) / 2.0, (lx1 - lx0)
        y_t = _plate_title(d, xc_, yb + fs * (3.6 + kt), fs * kt, ctx, align="m", deck=True, big=3.0,
                           wmax=wl_ * 0.96)
        unit_w = "solar radii" if ctx["unit"] == "sun" else "Moon radii"
        src_line = [_origin_words(ctx), f"map from the Moon\u2019s edge to {top:.1f} {unit_w}"] + _sources(ctx)
        _source_line(d, xc_, y_t + fs * 0.9 * kt, wl_ * 0.94, src_line, fs * 0.72 * min(kt, 1.5))
    return _pad_to(im, v.get("frame", "as"))


def _poster_social(rgb, ctx, v, soc, preview=False, rgb_map=None):
    """The poster at a social size: laid out for the part of the frame the app
    leaves free (the safe zone, or all of it), at that size, and set on black
    in the frame."""
    Wt, Ht = SOCIAL_SIZES[soc]
    if preview:
        Wt, Ht = int(Wt * 0.6), int(Ht * 0.6)
    top, right, bot, left, _cw, _ch = safe_margins(Wt, Ht, v.get("safe", "off"))
    # the poster's lettering is all at its foot and along the map: it keeps
    # clear of the bottom and the sides; its top is picture, which may run
    # under the app's header (cut there it showed a hard edge across the
    # corona), so it reaches the top of the frame
    top = 0
    bw, bh = Wt - left - right, Ht - top - bot
    v2 = dict(v, frame=f"{bw}:{bh}", _fs=SOCIAL_FS * (0.6 if preview else 1.0))
    if _cw:
        # the title block centred in what the column of buttons leaves free
        v2["_letter_x"] = (0.0, float(Wt - _cw - left))
    # the map is 2 pi unit radii wide: the picture enlarged or reduced so that
    # it comes out at the box's width (twice: a frame wider than the poster
    # pads its sides, and the second pass takes that out)
    f = bw / (2 * math.pi * ctx["Rs"])
    im = None
    for _ in range(2):
        h_, w_ = rgb.shape[:2]
        src = _shrink(rgb, max(2, round(w_ * f)), max(2, round(h_ * f)))
        src_m = _shrink(rgb_map, src.shape[1], src.shape[0]) if rgb_map is not None else None
        im = poster(src, _scale_ctx(ctx, f), v2, rgb_map=src_m)
        if abs(im.width - bw) <= 2:
            break
        f *= bw / float(im.width)
    im = im.resize((bw, bh), Image.LANCZOS) if im.size != (bw, bh) else im
    out = Image.new("RGB", (Wt, Ht), BG)
    out.paste(im, (left, top))
    return out


def _pad_to(im, frame):
    """The picture centred on black to the frame's shape (black is the
    background of every product, so nothing is added that is not there)."""
    if _aspect(frame) is None:
        return im
    aw, ah = _aspect(frame)
    W, H = im.size
    if abs(W * ah - H * aw) <= max(aw, ah):
        return im
    W2, H2 = (W, int(round(W * ah / aw))) if W * ah > H * aw else (int(round(H * aw / ah)), H)
    out = Image.new(im.mode, (W2, H2), BG if im.mode == "RGB" else (0, 0, 0, 0))
    ox, oy = (W2 - W) // 2, (H2 - H) // 2
    out.paste(im, (ox, oy))

    def _lay(lay):
        l2 = Image.new("RGBA", (W2, H2), (0, 0, 0, 0))
        l2.paste(lay, (ox, oy))
        return l2
    _cap_follow(im, out, _lay)
    cap = getattr(_CAP, "cap", None)
    if cap is not None:
        cap["off"][0] += ox
        cap["off"][1] += oy
    return out


# ---------------------------------------------------------------- shared: the limb strip
def _limb_strip(src, ctx, r1):
    """1:1 unwrapped strip from the Moon's edge (bottom) out to r1 (top), one
    column per pixel of circumference at the limb, the angle falling to the
    right with south in the middle. Returns strip, th, rr (in units)."""
    Rs = ctx["Rs"]
    r0 = ctx["Rm"] / Rs - 0.01
    # a multiple of 12, so the panorama's 2, 3 or 4 rows each start and end
    # exactly on N, E, S or W (or every 120 deg)
    nt = max(12, int(round(2 * math.pi * Rs * r0 / 12.0)) * 12)
    return _flat_map(src, ctx, r0, r1, centre="S", width=nt)


def prominences(strip, th, rr, ctx, sens=5.0):
    """Prominences along the limb, from the strip: red excess (R/G) over a
    running median along the limb, grouped, with height from the outward walk.
    Returns a list of dicts: pa, width_deg, height (units), height_km (or
    None), score; by position angle."""
    nt = strip.shape[1]
    Rm = ctx["Rm"] / ctx["Rs"]
    win = max(9, int(0.078 * nt) | 1)
    rg = strip[..., 0] / np.maximum(strip[..., 1], 1e-3)
    rows = (rr > Rm + 0.004) & (rr < Rm + 0.04)
    if rows.sum() < 1:
        return []
    prof = rg[rows].max(0)
    exc = prof - ndimage.median_filter(prof, size=win, mode="wrap")
    noise = 1.4826 * np.median(np.abs(exc - np.median(exc)))
    if not np.isfinite(noise) or noise <= 0:
        return []
    lab, n = ndimage.label(exc > sens * noise)
    min_cols = max(2, int(0.4 / 360.0 * nt))
    groups = []
    for i in range(1, n + 1):
        idx = np.where(lab == i)[0]
        if len(idx) < min_cols:
            continue
        # angular extent, wrap-safe through the column index
        groups.append([int(idx.min()), int(idx.max()), float(exc[idx].max() / noise)])
    groups.sort()
    merged = []
    gap_cols = int(8.0 / 360.0 * nt)
    for a, b, sc in groups:
        if merged and a - merged[-1][1] < gap_cols:
            merged[-1] = [merged[-1][0], b, max(merged[-1][2], sc)]
        else:
            merged.append([a, b, sc])
    if len(merged) > 1 and merged[0][0] + nt - merged[-1][1] < gap_cols:   # across the seam
        a, b, sc = merged.pop()
        merged[0] = [a - nt, merged[0][1], max(sc, merged[0][2])]
    ex2 = rg - ndimage.median_filter(rg, size=(1, win), mode="wrap")
    nz2 = 1.4826 * np.median(np.abs(ex2[rows]))
    out = []
    for a, b, sc in merged:
        cols = np.arange(a, b + 1) % nt
        e = ndimage.uniform_filter1d(ex2[:, cols].max(1), 5)
        order = np.where(rr > Rm)[0][::-1]
        rtop = Rm
        for kk in order:
            if e[kk] > 0.9 * sens * nz2:
                rtop = rr[kk]
            elif rr[kk] > Rm + 0.01:
                break
        pa_mid = float(th[(a + (b - a) // 2) % nt])
        width = float((b - a + 1) * 360.0 / nt)
        h = float(max(0.0, rtop - Rm))
        out.append(dict(pa=round(pa_mid, 1), width_deg=round(width, 1), height=h,
                        height_km=(round(h * KM_PER_RSUN / 1000.0) * 1000.0 if ctx["unit"] == "sun" else None),
                        score=round(sc, 1), col=(a + (b - a) // 2) % nt))
    out.sort(key=lambda p: p["pa"])
    return out


_DPROM = {}


def data_prominences(ctx, sens=5.0):
    """The prominences found on the stack itself (hdr_rgb.npy, linear), not
    on the rendered look: a mono look with tints, or a warm or a cold one,
    would otherwise move the red excess the search reads and find more or
    fewer of them. The position angles are the picture's (they follow its
    turn); `col` is left out, each drawing finds its own column. None when
    the stack is not there."""
    wd = ctx.get("wd")
    hp = os.path.join(wd or "", "hdr_rgb.npy")
    gp = os.path.join(wd or "", "geometry.json")
    if not wd or not os.path.exists(hp) or not os.path.exists(gp):
        return None
    s = float(ctx.get("scale", 1.0))
    zero = (float(ctx["o"]["zero_cw"]) + float(ctx.get("turned", 0.0))) % 360.0
    Rs = float(ctx["Rs"]) * s
    key = (hp, os.path.getmtime(hp), round(zero, 3), round(Rs, 2), int(ctx["o"]["sgn"]), ctx["unit"], float(sens))
    if key in _DPROM:
        return [dict(p) for p in _DPROM[key]]
    geo = json.load(open(gp))
    from .pipeline import load_big
    hdr = load_big(hp)
    H, W = hdr.shape[:2]
    cy, cx = float(geo["cy"]), float(geo["cx"])
    m = int(math.ceil(1.4 * Rs)) + 2
    y0, y1 = max(0, int(cy) - m), min(H, int(cy) + m + 1)
    x0, x1 = max(0, int(cx) - m), min(W, int(cx) + m + 1)
    crop = np.asarray(hdr[y0:y1, x0:x1, :3], np.float32)
    del hdr
    o = dict(ctx["o"])
    o["zero_cw"] = zero
    cd = dict(ctx, cy=cy - y0, cx=cx - x0, R=float(geo["R"]), Rm=float(geo.get("Rmask", geo["R"])),
              Rs=Rs, o=o, scale=1.0, turned=0.0)
    crop = np.maximum(crop, 0.0)
    strip, th, rr = _limb_strip(crop, cd, 1.35)
    out = prominences(strip, th, rr, cd, sens=sens)
    for p in out:
        p.pop("col", None)
    if len(_DPROM) > 8:
        _DPROM.pop(next(iter(_DPROM)))
    _DPROM[key] = [dict(p) for p in out]
    return out


def find_prominences(src, ctx, strip=None, th=None, rr=None, sens=5.0):
    """The prominences for a drawing: from the stack when it is there (the
    same ones in every look), else from the picture. `col` is the nearest
    column of this drawing's own strip `th` (or of a 1:1 limb strip)."""
    try:
        out = data_prominences(ctx, sens=sens)
    except Exception:
        out = None
    if strip is None:
        strip, th, rr = _limb_strip(src, ctx, 1.35)
    if out is None:
        return prominences(strip, th, rr, ctx, sens=sens)
    tha = np.asarray(th, np.float64)
    for p in out:
        dd = np.abs((tha - p["pa"] + 180.0) % 360.0 - 180.0)
        p["col"] = int(np.argmin(dd))
    return out


# ---------------------------------------------------------------- 3 limb panorama
def limb(rgb, ctx, v):
    src = np.asarray(rgb, np.float32)
    r1 = float(v["outer"])
    strip, th, rr = _limb_strip(src, ctx, r1)
    proms = find_prominences(src, ctx, strip, th, rr) if v["marks"] else []
    z = int(v["zoom"])
    nr, nt = strip.shape[:2]
    rows = int(v["rows"])
    q = nt // rows
    s8 = Image.fromarray(_to8(strip))
    if z > 1:
        s8 = s8.resize((nt * z, nr * z), Image.LANCZOS)
    sh, sw = nr * z, q * z
    fs = max(7.0, sw / 95.0)
    left = int(fs * 4.5)
    head = int(fs * 4.6)
    gap = int(fs * 2.8)
    pano = Image.new("RGB", (left + sw + int(fs), head + rows * (sh + gap) + int(fs)), BG)
    d = Draw(pano)
    unit = _unit(ctx)
    ft = _font(fs * 1.7, True)
    d.text((fs, fs * 0.9), "Limb panorama", fill=INK, font=ft)
    sub = (f"{ctx['Rm'] / ctx['Rs']:.2f}\u2013{r1:.2f} {unit} unwrapped at {z}:1 \u00b7 "
           + (_origin_words(ctx) + ", falling to the right as on the sky"
              if ctx["o"]["origin"] in ("solar", "celestial") else _origin_words(ctx))
           + (" \u00b7 \u25bc prominences found, their height above the Moon\u2019s edge" if proms else ""))
    d.text((fs + _text_w(d, "Limb panorama", ft) + fs * 1.2, fs * 1.45), sub, fill=(150, 150, 150),
           font=_font(fs * 0.95))
    for kk in range(rows):
        y = head + kk * (sh + gap)
        pano.paste(s8.crop((kk * sw, 0, (kk + 1) * sw, sh)), (left, y))
        for pa in range(0, 360, 10):
          for col in _angle_cols(th, pa):
            xx = left + (col + 0.5) * z - kk * sw
            if left - 1 <= xx <= left + sw + 1:
                xx = min(max(xx, left), left + sw)
                d.line([(xx, y + sh), (xx, y + sh + fs * 0.4)], fill=(165, 165, 165), width=1)
                if pa % 30 == 0:
                    lab = _tick_lab(ctx, pa)
                    anc = "lt" if xx < left + fs else ("rt" if xx > left + sw - fs else "mt")
                    d.text((xx, y + sh + fs * 0.55), lab, fill=GOLD if len(lab) == 1 else (205, 205, 205),
                           font=_font(fs * 0.9, len(lab) == 1), anchor=anc)
        for rv in np.arange(1.1, r1 - 0.01, 0.1):
            yy = y + (r1 - rv) / (r1 - rr[-1]) * sh
            d.text((left - fs * 0.4, yy), f"{rv:.1f}", fill=(135, 135, 135), font=_font(fs * 0.8), anchor="rm")
    for p in proms:
        xg = p["col"] * z
        kk = int(xg // sw)
        xx = left + xg - kk * sw
        y = head + kk * (sh + gap)
        d.polygon([(xx - fs * 0.4, y - fs * 0.7), (xx + fs * 0.4, y - fs * 0.7), (xx, y - fs * 0.1)],
                  fill=(255, 140, 170))
        if p["height_km"] is not None and p["height_km"] >= 1000:
            lab = f"\u2248 {_round_km(p['height_km']):,.0f} km".replace(",", " ")
            right = xx + fs * 0.8 + _text_w(d, lab, _font(fs * 0.85)) < left + sw
            d.text((xx + fs * 0.8 if right else xx - fs * 0.8, y - fs * 0.45), lab, fill=(255, 185, 205),
                   font=_font(fs * 0.85), anchor="lm" if right else "rm")
    return pano, proms


def _fmt_s(t):
    return f"{t:g} s" if t >= 1 else f"1/{round(1 / t):g} s"


# ---------------------------------------------------------------- entry
def make(kind, rgb, geo, wd, folder, opts, scale=1.0, north="solar", manual_cw=None, offset_deg=0.0, prog=None,
         preview=False, out_base=None, rgb_map=None):
    """One product from the composite `rgb` (float 0-1, H x W x 3, the stack
    frame at 1/scale of the stack's size). Returns dict(image=PIL RGB,
    layer=PIL RGBA or None, fmt, info=str, data=dict)."""
    opts = dict(opts or {}, _angles=[north, manual_cw, offset_deg])
    v = option_values(kind, opts)
    if preview and rgb is not None and kind in ("overlay", "poster", "limb", "sound") \
            and rgb.shape[1] < 1400:
        # a small preview picture (a small sensor, or the 1/4 layers of a
        # small frame) drawn on as it is gets lettering far too big for it:
        # enlarged first, so the card shows the proportions of the file
        k = int(math.ceil(1400.0 / rgb.shape[1]))
        im_ = Image.fromarray(_to8(rgb)).resize((rgb.shape[1] * k, rgb.shape[0] * k), Image.BILINEAR)
        rgb = np.asarray(im_, np.float32) / 255.0
        if rgb_map is not None:
            im_ = Image.fromarray(_to8(rgb_map)).resize((rgb.shape[1], rgb.shape[0]), Image.BILINEAR)
            rgb_map = np.asarray(im_, np.float32) / 255.0
        scale = float(scale) / k
    ctx = context(wd, folder, geo, scale=scale, north=north, manual_cw=manual_cw, offset_deg=offset_deg)
    if kind in LAYERED and v.get("files") == "layers" and not preview \
            and not (kind == "poster" and v.get("social", "full") in SOCIAL_SIZES):
        # (a poster at a social size is set into its frame afterwards: one file)
        _CAP.cap = dict(layers={}, rec=[], own_layer=(kind == "overlay"), off=[0, 0])
        try:
            res = _make(kind, rgb, ctx, geo, wd, folder, opts, v, scale, prog, preview, out_base, rgb_map)
            cap = _CAP.cap
        finally:
            _CAP.cap = None
        im = res["image"]
        if kind == "overlay":
            base, layer = res["base"], res["layer"]
        else:
            base = im
            layer = cap["layers"].get(id(im), (None, Image.new("RGBA", im.size, (0, 0, 0, 0))))[1]
            comp = im.convert("RGBA")
            comp.alpha_composite(layer)
            res["image"] = comp.convert("RGB")
            res["layer"] = layer
        res["base"] = base
        res["svg"] = (base.width, base.height, cap["rec"], cap["off"][0], cap["off"][1])
        return res
    return _make(kind, rgb, ctx, geo, wd, folder, opts, v, scale, prog, preview, out_base, rgb_map)


LAYERED = ("overlay", "poster", "limb")


def _make(kind, rgb, ctx, geo, wd, folder, opts, v, scale, prog, preview, out_base, rgb_map=None):
    if rgb is not None and (rgb.ndim == 2 or rgb.shape[2] == 1):
        rgb = np.repeat(np.asarray(rgb, np.float32).reshape(rgb.shape[0], rgb.shape[1], 1), 3, 2)
    soc = v.get("social", "full") if kind in ("overlay", "poster") else "full"
    if soc in SOCIAL_SIZES:
        v = dict(v, frame=soc, _social=True)
    if kind == "overlay":
        fy = 0.5
        if soc in SOCIAL_SIZES and SOCIAL_SIZES[soc][1] >= SOCIAL_SIZES[soc][0]:
            # an upright or square social frame: the facts and the scale take a band at
            # the foot; the Sun sits in the middle of what is left above it
            Wt, Ht = SOCIAL_SIZES[soc]
            t_, _r, b_, _l, _c, _h = safe_margins(Wt, Ht, v.get("safe", "off"))
            band = SOCIAL_FS * float(v["size"]) * 13.0
            fy = min(0.5, (t_ + (Ht - b_ - band)) / 2.0 / Ht)
        rgb, ctx = reframe(rgb, ctx, v.get("turn", "none"), v.get("frame", "as"), fy=fy)
        if soc in SOCIAL_SIZES:
            # laid out at the size it is shown at, not shrunk afterwards
            # (hairlines of a 6000 px sheet vanish at 1080)
            Wt, Ht = SOCIAL_SIZES[soc]
            if preview:
                Wt, Ht = int(Wt * 0.6), int(Ht * 0.6)
            f = Wt / float(rgb.shape[1])
            rgb = _shrink(rgb, Wt, Ht)
            ctx = _scale_ctx(ctx, f)
            if preview:
                v = dict(v, size=str(float(v["size"]) * 0.6))
        v = dict(v, _safe_px=safe_margins(rgb.shape[1], rgb.shape[0], v.get("safe", "off")))
        im, layer = overlay(rgb, ctx, v)
        return dict(image=im, layer=layer, fmt=v["fmt"], base=Image.fromarray(_to8(rgb)),
                    info=f"overlay on {im.width} x {im.height}; " + ctx["o"]["note"], data={})
    if kind == "poster":
        if v.get("maplook", "same") == "same":
            rgb_map = None
        how = ("; the map in the standard look" if rgb_map is not None else
               ("; the map in the picture's look (tab 3's look is not a soft one)"
                if v.get("maplook", "same") == "standard" else ""))
        if soc in SOCIAL_SIZES:
            im = _poster_social(rgb, ctx, v, soc, preview, rgb_map=rgb_map)
            return dict(image=im, layer=None, fmt=v["fmt"],
                        info=f"poster {im.width} x {im.height}{how}; " + ctx["o"]["note"], data={})
        im = poster(rgb, ctx, v, rgb_map=rgb_map)
        return dict(image=im, layer=None, fmt=v["fmt"], info=f"poster {im.width} x {im.height}{how}; " + ctx["o"]["note"],
                    data={})
    if kind == "limb":
        im, proms = limb(rgb, ctx, v)
        return dict(image=im, layer=None, fmt=v["fmt"],
                    info=f"limb panorama {im.width} x {im.height}; {len(proms)} prominences marked; " + ctx["o"]["note"],
                    data={"prominences": [{k: p[k] for k in ("pa", "width_deg", "height_km", "score")} for p in proms]})
    if kind == "film":
        from . import film as _film
        im = _film.storyboard(rgb, geo, wd, folder, opts, scale=scale)
        return dict(image=im, layer=None, fmt="jpg", info="film storyboard", data={})
    if kind == "sound":
        from . import sound as _sound
        im = _sound.still(rgb, geo, wd, folder, opts, scale=scale)
        return dict(image=im, layer=None, fmt="jpg", info="sound clip still", data={})
    if kind == "activity":
        from . import activity as _act
        im, extra, info = _act.make(wd, geo, folder, dict(v, _look=(opts or {}).get("_look") or {},
                                                          _angles=(opts or {}).get("_angles")),
                                    out_base or "", prog=prog, preview=preview)
        return dict(image=im, layer=None, fmt=v.get("fmt", "jpg"), info=info, data={}, extra_paths=extra)
    raise ValueError(f"unknown product {kind!r}")


def write(result, base):
    """Writes the product's files next to `base` (path without extension).
    Returns the list of paths."""
    paths = []
    fmt = result.get("fmt", "jpg")
    p = base + ("." + ("png" if fmt == "png" else "jpg"))
    k = 0
    while os.path.exists(p):
        k += 1
        p = base + f"_{k}." + ("png" if fmt == "png" else "jpg")
    if fmt == "png":
        result["image"].save(p, "PNG", optimize=False)
    else:
        result["image"].save(p, "JPEG", quality=93, subsampling=0)
    paths.append(p)
    stem = base + (f"_{k}" if k else "")
    if result.get("layer") is not None:
        lp = stem + "_layer.png"
        result["layer"].save(lp, "PNG")
        paths.append(lp)
    if result.get("svg") is not None and result.get("base") is not None:
        # the picture without lettering, and the lettering as vectors over it
        bp = stem + "_picture." + ("png" if fmt == "png" else "jpg")
        if fmt == "png":
            result["base"].save(bp, "PNG", optimize=False)
        else:
            result["base"].save(bp, "JPEG", quality=95, subsampling=0)
        paths.append(bp)
        W_, H_, rec, dx, dy = result["svg"]
        sp = stem + "_lettering.svg"
        with open(sp, "w", encoding="utf-8") as fh:
            fh.write(_svg_doc(W_, H_, os.path.basename(bp), rec, dx, dy))
        paths.append(sp)
    if result.get("data"):
        with open(paths[0] + ".json", "w") as f:
            json.dump(result["data"], f, indent=1)
    return paths
