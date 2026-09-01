# Changelog — EclipseForgeHDR

Newest first. Entries from 0.6.1 onward were written at the time. The 0.7.2 –
0.9.8 block was reconstructed afterwards from the code comments that name their
own version and from the development record; where a change cannot be pinned to
an exact version it is filed under the release it is known to precede.

## 0.14.5 — The black annulus in the NAFE layer

Reported from a real render: NAFE showed a thick black ring hugging the disc
and extending well outside it, with a bright rim at the limb, where FNRGF was
clean.

The cause is the envelope flattening added in 0.14.0, and specifically its
kernel rather than its idea. A symmetric Gaussian mean is a poor background
estimate beside a large dark hole: near the limb the mean is pulled down by the
occulted disc, so `L − mean` overshoots bright at the limb and undershoots dark
just outside it, roughly 2σ wide. At 0.08 R on a 622 px disc that is about
100 px of black — an unsharp halo, textbook.

0.14.0's own measurement missed it because "ripple" was an RMS, and a broad
smooth depression barely moves an RMS. It is obvious in the **mean of the layer
per radial ring**, which is how it is measured now.

The local mean is now built by **normalized convolution over non-disc pixels
only** — the same treatment MGN already gives its statistics. Measured on a
synthetic corona carrying a *known* fine modulation, so that "detail" is a
correlation with the truth rather than a variance the artifact itself inflates:

| flattening | ring depth | correlation with the true structure |
|---|---:|---:|
| plain Gaussian (0.14.0–0.14.4) | 0.1684 | 0.468 |
| normalized convolution | **0.0713** | **0.814** |

Mean layer value at 1.15 R went from 0.333 to 0.445 against a far-field 0.50 —
the visible depth of the ring, more than halved — and the streamers now run
cleanly down to the limb instead of stopping at a dark band.

### Why not the option that scored better

Subtracting the Fourier radial background (what FNRGF and the inner corona use)
scored far better still on a perfect limb fit — 0.0047 and 0.938 — but it needs
a centre and a radius, and it fails hard when they are wrong:

| centre error | ring depth | correlation |
|---|---:|---:|
| exact | 0.0047 | 0.938 |
| 0.02 R | 0.1347 | 0.666 |
| 0.05 R | 0.2267 | 0.542 |
| 0.10 R | 0.2533 | 0.416 |

Past about 0.05 R it is **worse than what it replaces**: it manufactures its own
artifact. Normalized convolution needs only the disc mask, which the renderer
already uses to fill the disc, so it adds no dependency that was not already
there — and it degrades gently rather than inverting. With the centre 0.10 R out
it still scores 0.1135 / 0.731, better than the plain Gaussian with a *perfect*
fit. Across every error tested it was never worse than the version it replaces.

σ stays at 0.08 R. Re-swept with the disc excluded, 0.15 R measured slightly
better (0.0576 / 0.863), but σ interacts with the real scale of coronal
structure, which a synthetic built from azimuthal cosines does not faithfully
represent, while the kernel correction does not. One change at a time; the sweep
is worth repeating on real data.

Re-running is required — the NAFE layer itself changes.

## 0.14.4 — Progress weights from a measured run; report span excludes rejects

### The progress split, measured instead of assumed

0.14.1 gave the detail stage 22% of the bar on an assumption — that stacking
costs about 3.5× the detail stage — because the stacking side had never been
timed on real RAWs. The run-time summary added in that same release then
measured it, on a 50-frame 14-tier bracket from a 45 Mpx body, 12m27s end to
end:

| step | time | share of run |
|---|---:|---:|
| inner corona: raw pass | 141 s | 19% |
| inner corona: denoised pass | 140 s | 19% |
| MGN | 31 s | 4% |
| **three named detail steps** | **312 s** | **42%** |

The remaining detail steps are each below the report's 6th entry (21 s), which
brackets the whole stage between 42% and 64% of the run. So stacking and detail
are roughly equal, not 3.5:1, and the bar was reaching 78% around the halfway
mark and then crawling. The detail stage now starts at 52% of the bar rather
than 78% — the low end of the measured bracket, since a bar that lags is better
than one that arrives at 100% and waits.

This is one dataset, so it is an estimate with a range rather than a constant,
and it will move with frame count (more frames, more stacking) and sensor size
(detail grows faster than pixel count). To let the next runs settle it without
extrapolation, the report now lists **every step over a second** with its
percentage, plus a line for what the listed steps do not account for — instead
of stopping at the slowest four and leaving the rest of the run unexplained.

### The report's shot span counted rejected frames

The bracket's time span and ISO list were taken from every *file* in the
folder, before any frame was rejected. A test frame shot three weeks after the
eclipse and correctly dropped as "not a totality frame" still set the span, so
the report read:

    frames taken : 2026:08:12 21:30:12  ..  2026:09:01 13:35:30

Three weeks of totality. Both are now computed from the frames actually
stacked. Verified against a control: with the stray frame dropped, the span
matches the clean run exactly (10:20:00–10:23:14 rather than 10:20:00–10:59:00);
where a stray frame is genuinely used, it still counts.

## 0.14.3 — A stray overexposed frame in the raw folder

Found by testing what happens when one blown file is dropped into an otherwise
healthy folder. The 0.14.0 guard covers a saturated *tier*; a stray frame is
the same problem in three positions, and where it lands decides what it can
damage. All three finish without raising, both before and after this release —
the crash guard holds — but two of them were quietly doing the wrong thing.

### A saturated shortest tier disabled hot-pixel repair

Hot pixels are mapped once, on the shortest tier, because its sky is darkest.
A hot pixel is a photosite far above its neighbours — and in a saturated frame
there is no "above": everything sits at the clip. A stray overexposed file
whose EXIF puts it at the short end becomes the shortest tier on its own, and
the map built from it found **zero** defects on a sensor carrying 434. The
existing guard only caught an implausibly *high* count, so this passed silently
and hot-pixel repair was off for the entire run.

A tier that is more than half clipped is now skipped for this purpose and the
next tier is used, with a line saying so.

| defect map built from | before | after |
|---|---:|---:|
| healthy shortest tier (control) | 434 | 434 |
| stray blown frame as the shortest tier | **0** | 434 |

### The limb-radius warning claimed an override it had not made

When the merged limb fit disagrees with the per-tier consensus by more than
15%, the log said:

> ... the individual tiers agree on R=80px (+22%). **Using the tiers' value** —
> the merge probably contains frames of different scenes.

The code only actually adopts the tiers' value past **30%**, or when there is
no per-azimuth fit. Between 15% and 30% it printed that sentence and then kept
the merged fit. This is the one line a user reads when the disc mask comes out
the wrong size, so it was misdirecting exactly when it mattered.

The behaviour is unchanged — a per-azimuth measurement is not obviously worse
than the tiers' single number at 22% — but each branch now says what it did,
and the keeping branch points at what to check:

> ... **KEEPING the merged fit**: it is a per-azimuth measurement and the
> disagreement is under 30%, where the tiers' single number is not clearly the
> better one. Check the disc mask on the preview — if it is the wrong size, the
> tiers were right.

Whether 30% is the right place for that line is a separate question, and not
one a synthetic frame can answer.

### Unchanged, and confirmed

A stray blown frame *inside* a healthy tier was already handled correctly and
still is — it is identified against its own tier and dropped by name:

    1/30s: dropping ZSTRAY.rw2 — 100.00% of the frame is saturated vs 5.11%
    for this tier; this is not a totality frame

Cached products from 0.14.2 and earlier are **not** reused by this release: the
defect map and the limb radius both feed the merge, so the products genuinely
differ. Re-running is required.

## 0.14.2 — Flat view rendered the prominence gate

0.14.1 got the button back and then drew the wrong layer into it. The server
was serving the flat correctly — `/api/layer/flat` and the export path were
both right — but the browser's preview picked its single-layer source with a
ternary chain that ended in a bare fallback to the prominence gate:

    VIEW==="pellett" ? L.pel[j]/255 : pr[j]/255

Every view spelled out in that chain worked, and `prom` "worked" only by being
the fallback. Any view added later — Flat, when it arrived — silently drew the
prominence gate instead. Measured, before and after, as the mean absolute
difference between what each button draws:

| | before | after |
|---|---:|---:|
| Flat vs Prom gate | **0.00** (identical) | 160.69 |
| Flat vs what `/api/layer/flat` serves | — | r = 1.0000 |

Each view now names its own source, and a view whose layer is missing falls
back to the composite with a console warning rather than to whichever layer
happened to sit last in the chain. All eight views were checked pairwise; none
of them draws the same thing as any other.

This is the second bug in three days where the flat preview failed by showing
something plausible instead of nothing. Both were invisible for the same
reason: no test compared what the user sees against what the pipeline built.
There is one now (`viewtest.py`), and it drives the real page in a real browser
rather than the functions underneath it.

## 0.14.1 — Master-flat preview and progress reporting

### Master-flat preview failed to load

0.14.0 added the master flat as a quality-control preview and the button for it
never showed up on the first real run. The loader reduced the flat to
superpixels — correctly, to kill the Bayer checkerboard — and then cropped that
half-resolution array to the **full-resolution** layer grid. On the reporting
frame that meant asking for 5358×8100 out of a 2716×4076 array, which can never
match, so the size check failed and the loader returned `False`. It failed
silently, which is why it shipped: a missing button looks exactly like a run
with no flats.

Fixed by putting the reduced flat back on the full grid (crop on the superpixel
grid first, then expand, so the temporary is the size of the view and not of
the sensor). Verified on the reporting shapes and three others — uncropped, an
odd-sized layer grid, and an old cache with no `crop_origin`:

| case | before | after |
|---|---|---|
| 5432×8152 sensor, 5358×8100 grid, crop (32,34) | no button | loads, dust visible, no checkerboard |
| no autocrop | loads | loads |
| odd layer grid | loads | loads |
| old cache, no `crop_origin` | loads (centred) | loads (centred) |

And it no longer fails silently: the reason goes to the log and to
`/api/geometry`, so the next one of these is visible rather than invisible.

### Progress reporting reweighted by measured time

Reported from the field: the bar sits nearly full for a long time, and the only
sign the run is alive is the terminal. Both halves of that are now addressed.

**Where the time actually goes.** The detail stage was timed step by step on
synthetic corona frames at two sizes:

| step | 10.8 Mpx | 43.4 Mpx |
|---|---:|---:|
| denoise HDR master | 2.9 s | 18.2 s |
| MGN | 26.9 s | 197.9 s |
| FNRGF | 7.6 s | 22.9 s |
| NAFE | 11.1 s | 56.3 s |
| **inner corona** | **135.3 s** | **1001.1 s** |
| prominence colour | 0.9 s | 3.6 s |
| Pellett | 11.1 s | 31.7 s |
| **total** | **198.7 s** | **1343.9 s** |

Two things came out of that. The inner-corona block is two thirds to three
quarters of the stage on its own — it runs the multiscale normalisation twice
at full resolution, once raw and once denoised (61.0 s and 61.7 s of its 135.3 s
at 10.8 Mpx) — and it had no logging inside it at all, so the bar and the log
both stood still through most of the stage. And the stage does not scale with
pixel count: 4× the pixels cost 6.8× the time, so on a 45 Mpx body the detail
layers alone are 22 minutes. The bar was giving all of that 6.5% of its width.

So: the detail stage now gets 22% of the bar instead of 6.5%, divided by
measured wall time rather than evenly, and the inner-corona block reports its
four sub-steps. Its two long passes get 7.4% of the bar each, where before the
pair got about 1%. The 22% is still an estimate — it assumes stacking takes
about 3.5× the detail stage, and the stacking side has never been timed on real
RAWs — so it errs towards lagging, because a bar that arrives at 100% and waits
is worse than one that arrives late.

**A run that says how long it took.** Every finished run now prints its own
timing summary — total, and the slowest steps by name — to the log and to
`report.txt`. That is what will replace the estimate above with real numbers
from real cameras.

**Signs of life.** Beside the bar there is now a spinner, an elapsed clock
counting from the server's own start time, and the name of the current step.
Once a step has been running for more than 45 seconds the line says which step
it is and for how long, so a long wait is legible instead of alarming. The bar
itself carries moving stripes while a run is live (and doesn't, when it isn't).
The importing and diamond-ring paths get the same treatment; import lowers the
detail band to 12%, since it skips stacking and the detail layers are nearly
the whole job there.

### Cache compatibility

0.14.1 changes nothing that is written into the work directory, so a folder
already processed with 0.14.0 is still reusable and does not need re-stacking
for the sake of a version number. Builds are declared interchangeable
explicitly, in `CACHE_COMPAT`, rather than by relaxing the check.

## 0.14.0 — Alignment crash guard, HDR import, NAFE input flattening

### Alignment failure on a fully saturated tier

Two people hit `SVD did not converge in Linear Least Squares` partway through a
run, on Windows, and the theory going round was that it was a Windows problem.
It is not, and it is not memory either.

`prep_pair` masks out saturated pixels before correlating a pair of tiers. When
a tier's whole correlation window is saturated the weight map is zero
everywhere, so the prepped image is exactly zero. Phase correlation on a zero
image returns a finite but meaningless shift and an **`err` of NaN** — and the
weight `1/(err+0.05)` is therefore NaN. One NaN weight poisons a full row of the
design matrix, and both `lstsq` calls die:

| input | returned shift | returned `err` | weight |
|---|---|---|---|
| two real images | [0, 0] | 0.0009 | 19.6 |
| one all-zero | [−0.75, −0.75] | **nan** | **NaN** |

Reproduced with a synthetic bracket whose longest tier is blown — same warning,
same `DLASCL parameter number 4 had an illegal value`, same exception, on Linux.
It takes a bracket wide enough to blow one end; the report that found it spans
14.3 EV with single frames at 4.2, 11.9 and 20 s.

Now: degenerate windows are detected before correlating, non-finite shifts and
weights can never enter the solve, the solve itself is wrapped, and a tier that
nothing could link to takes the shift of its nearest linked neighbour rather
than a minimum-norm number that looks like an answer. Every one of those says so
in the log, naming the tier.

### Import of a finished HDR

A folder of raws is no longer the only way in. Point the new box at one 16-bit
TIFF or FITS — from Siril, PixInsight, Photoshop, or this app's own aligned tier
exports — and it runs the disc fit, the sky-gradient fit and every detail layer
on that image. All seven views export exactly as they do from a stack.

**The tone curve is read from the file, not guessed.** MGN, FNRGF and NAFE all
work on log luminance and assume the value is proportional to coronal
brightness. Measured on a Photoshop mean stack of this app's own sRGB tier
exports: the corona's log-log radial slope read **−1.72** against **−3.38** for
the same scene linear, a ratio of 0.510 where sRGB predicts ~0.45. Applying the
inverse the embedded ICC profile declares brought it to **−3.17**, within 6% of
the truth. Files with no profile must declare themselves; 8-bit input is
refused, because the corona spans several thousand to one.

The report says plainly what an imported image gives up — no alignment,
photometry or per-tier lunar masking, no earthshine, the inner-corona layer is a
second view of the same pixels rather than an independent measurement, and the
prominence gate is weaker because a merged image has usually compressed the
H-alpha contrast it keys on (R/GB 1.84 against 3.02 from a real fast tier).

### NAFE input flattening

NAFE was the only detail layer whose input still carried the full radial
falloff — MGN and FNRGF both subtract it first. It paid twice: **60% of its
output range** went on a large-scale gradient instead of structure, and the
steep falloff outside the limb drove the equalisation into a dark ring.

The envelope now comes off with a plain Gaussian high-pass at 0.08 R, which
needs no circle and no limb — NAFE keeps the independence that makes it the
layer to trust when the limb fit is shaky. Measured three ways on the same
image, against the fitted radial profile MGN uses:

| | large-scale range | detail 1.05–1.5 R | ripple outside the mask |
|---|---|---|---|
| as shipped | 68.3% | 0.0611 | 0.3451 |
| fitted radial profile (needs the limb) | 5.4% | 0.1328 | 0.2543 |
| **Gaussian 0.08 R (no geometry)** | **3.6%** | **0.1360** | **0.1210** |

**2.2× the near-limb detail and 2.9× less ripple** where the dark ring used to
be — and the geometry-free option turned out to be the better one on every
count. The fitted profile also extrapolates inward and leaves a wedge across the
disc, which would have shown in an exported NAFE layer.

Worth recording what does *not* work, since it is the obvious thing to try:
subtracting a smooth model from the finished layer flattens it (68% → 9%) but
recovers **exactly zero** detail — identical to four decimals. The damage is
done by a nonlinear equalisation; a linear correction afterwards cannot undo it.

### Resolution-adaptive MGN scales

The scale ladders were fixed pixel counts. Now the top is tied to R — `0.0643 *
622 = 40 px` reproduces the reference set exactly — and the bottom is raised to
whatever the image actually resolves, measured from the band-pass falloff (white
noise falls as ~1/s, structure more slowly). At 3.19 arcsec/px behind a 240 mm
consumer zoom the old 1.25 and 2.5 px scales sat entirely below the optics: the
pixel-scale band carried 9.6% of the local mean against 1.9–2.7% for the bands
where the real structure lives. The reference ladder is unchanged.

### Other changes

- **Master flat preview.** Its own view button, cut to the layer grid via a new
  `crop_origin` in geometry.json, stretched to its own 0.5–99.5 percentile —
  a 16% falloff shown linearly over 0..1 is invisible, which is the point. Min,
  max and the stretch bounds are in the tooltip.
- **The master flat repairs its own outer border.** Some decoders hand back a
  row of masked photosites inside what they call the visible area; on the
  reference sensor row 0 reads 0.644 against 0.936 four rows in, which would
  brighten the top rows of every light frame by up to 55%. Edge lines are judged
  against their inward neighbours, never against a global median — the first
  version of this compared globally and erased 7–9 lines of perfectly good
  vignetting on every edge.
- Exported detail views carry their own p1/p99 range in the TIFF description.
  A view is rendered at the setting it contributes to the composite, not one
  that fills a histogram: MGN occupies 14% of the 16-bit range where FNRGF
  occupies 81%, and anyone compositing by hand needs to know before stretching.
- "Cached pipeline products found" now checks the build version, as Start does.

## 0.13.1 — Edge-based lunar disc detection

Fixes two failures that both came from the same wrong assumption: that you can
guess where the Moon is, and how big it is, from brightness and from the frame
size. One of them put a second black disc 1000 px from the Moon in a 0.13.0
render; the other has been silently disabling features on every short-focal-
length dataset since the beginning.

**The centre.** `find_center` took the centroid of the brightest **0.05%** of
pixels. That is the disc centre only while the bright inner corona rings the
disc evenly. When one sector dominates — a big prominence, an active region, a
lopsided inner corona — every one of those pixels sits on the same arc and the
centroid lands *on the limb*. Measured on the reference set: **646 px out on a
620 px disc, i.e. 1.04 R.**

I measured what the half-level fit downstream can actually absorb:

| seed centre offset | result |
|---|---|
| 0 – 500 px (0.81 R) | R 622.3, rms 5.00, 720/720 rays, **0.2 px from truth** |
| 620 px (1.00 R) | **fails** |

So the estimator had been sitting a few percent from a cliff, and a 3% change in
the limb ring — which is all the new flat-field correction does there — pushed it
over. Everything alarming in that run was downstream of it: the 67 px limb ramp
(was 21), `rim nan px`, "coronal range 1.2 EV" (was 6.5), the sky gradient
refitted about a point 1000 px from the Moon.

**The radius.** `fit_limb` searched for the limb between a tenth and a third of
the short side. That is a statement about focal length, not about eclipses: on a
24 MP APS-C body the limb is only inside that band beyond about **320 mm**, and
on 24 MP full-frame beyond about **510 mm**. A 240 mm frame puts the disc at
7.5% of the short side — below the floor — so the limb was never looked at. It
returned R = 1005 px for a 301 px disc, and every per-tier limb measurement
failed with it, taking prominence anchoring and per-tier lunar masking down too.

**What replaces both.** A new `find_disc()` locates the disc from the limb
*edge*. The limb is a huge relative step — 93× over ~20 px on the reference set
— while the corona's own falloff is smooth, so in log intensity it dominates the
gradient regardless of exposure or how bright one side is. A wide bright mask
gives a starting centre, the radius comes from the peak of the azimuthally
averaged log gradient, strong-gradient pixels in a band around that radius are
fitted with a circle, and the two alternate. Nothing in it refers to the frame
size, so the radius comes out of the data.

Two gates decide whether the answer is a limb at all: the inliers must go at
least 270° round, **and** their scatter about the circle must be tight. Coverage
alone is not enough — on the short-exposure stack, where the fast tiers see
nothing outside the inner corona, noise covers every azimuth and a meaningless
circle scores a full 360°. A real limb sits at 0.04 R of scatter; those noise
circles at 0.16–0.23 R.

Measured across disc sizes from 3% to 30% of the short side, on 3:2, square and
2.25:1 frames, with the inner corona lopsided up to 4:1 and prominences up to 3×
the limb brightness — 18 + 7 cases:

| | before | after |
|---|---|---|
| `fit_limb` R at 7.5% of the short side (240 mm) | **+344%** | **−1.2%** |
| `fit_limb` R at 3% of the short side | +1232% | −5.2% |
| centre, reference merged HDR | 646 px (1.04 R) | **62 px** |
| centre, short-exposure stack | 611 px | **1 px** |
| centre, lopsided inner corona | up to 1.03 R — past the cliff | 0–15 px |
| `fit_limb_rays` from that seed | failed | 720/720 rays, **0.2 px from truth** |
| seed cost | — | 0.08 s/frame, 4 s per 49-frame bracket |

The measured radius is now also carried into every place that used to seed a
limb fit with a frame fraction — the per-tier centres for the correlation
windows, the prominence-search ladder (rungs are now 0.7–2.0× the measured
radius instead of 6–24% of the frame), the moon-mask trial, and the merged-image
fallback. `prepare_contact` gets it through `fit_limb` for free.

Expect the knock-on effects to go with it: the same bad per-tier centres are why
that run claimed the tiers were "spread 607px apart" while the lunar track said
the Moon moved 1 px across the whole bracket, and why the alignment residual went
0.80 → 2.13 px half-res between two runs on identical files.

## 0.13.0 — Flat-field calibration

Put your flat frames in a subfolder of the raw folder called `flats/` and they
are found and used. No flats, nothing changes — the whole feature is inert when
the folder is absent. A path box next to the folder box takes flats from
somewhere else, or the word `off` to ignore a `flats/` folder that is there.

A master flat is built once per folder and cached, then divided out of every
frame of every tier before anything else touches it. It removes lens
vignetting, the cos⁴ falloff, dust shadows on the sensor stack and
per-photosite sensitivity. For an eclipse this matters more than for most
subjects, because the corona's own radial falloff *is* the signal: a 6%
vignette is a 6% error in the F-corona gradient, and every radial filter
downstream then works to preserve it.

**The master flat measures its own noise, and that sets how much it is
smoothed.** A flat is not free: dividing by it injects whatever noise it
carries into every frame identically, so stacking cannot average it away. The
frames are split into two independent half-stacks; because `var(mean of half) =
2 × var(mean of all)`, the difference between the two halves *is* the master's
noise rather than a proxy for it. The Gaussian σ is then raised until the
measured noise falls under 0.2% per photosite, and the numbers go in the log
and the run report:

```
flat: min/max-trimmed mean of 20 frames from flats/
flat: per-pixel noise 1.318% -> 0.183% after a 5.7 px smooth (target 0.2%)
flat: corrects a 8.4% falloff — the dimmest part of the field sits at 0.923
      of the brightest
```

So a clean, well-exposed flat set keeps its dust motes at full resolution (σ
lands at 0 or 1 px), and a thin or under-exposed one degrades gracefully to a
vignetting model instead of adding more noise than it removes. The first guess
at σ comes from a white-noise formula, but sensor noise is not white — on the
reference flats it needed 2.9× more smoothing than the formula predicted — so
the guess is only a starting point and the answer is measured.

Other things that fall out of doing this properly:

- **Per Bayer channel, normalised to the centre of the frame.** Each flat is
  divided by its own central median *per channel*, so exposure and
  illumination-colour drift between flats cancels and the master cannot shift
  the white balance or the overall level of the lights — only their spatial
  structure. Measured end to end on a synthetic bracket: the central level
  moves by 0.4%, and the correction the merged HDR actually received matches
  the vignette that was baked into the frames to 0.9% median over 60% of the
  field.
- **Min/max-trimmed per-pixel combine** (5 frames or more), so a cosmic ray, a
  satellite or a bird in any one flat is dropped rather than averaged in at
  1/N. Tested with a synthetic trail: 15% contamination becomes 1.4%, which is
  the master's own noise floor.
- **The clipping test is made before the flat is divided out.** This is the
  part that would have been a silent regression. A vignetted corner is
  brightened by the correction, and a pixel compared against the scalar
  saturation level *afterwards* reads as clipped at a fraction of the well it
  actually filled. With a 30% vignette, a corner pixel that deserves full
  weight 1.000 in the merge would have been given 0.269 — the longest tier
  would have lost three-quarters of its contribution in the corners. Per-frame
  saturation masks are now taken from the uncorrected data, and in the merge
  and the inner stack the flat-corrected value is put back into raw units for
  the comparison.
- **Flats that are not usable say so, one line each**: exposed above 85% of
  saturation (on the shoulder of the response), below 2% (too dark to
  calibrate with), or a different frame size from the lights. A flat set that
  cannot produce a master leaves a message in the log and the run continues
  uncorrected rather than failing.
- The flat files are part of the cache key, so adding, replacing or removing
  one re-runs the stack instead of silently reusing the old one; a contact
  frame loaded later gets the same flat, but only if the run that built the
  composite actually applied one.

This does not replace the sky-gradient fit from 0.11.x. Vignetting is radial
about the frame centre and fixed by the optics; the low-altitude sky gradient
is a tilted plane fixed by the atmosphere, and on the reference set a radial
model explains 0.0% of it. Both run, in that order, and each is reported
separately.

## 0.12.0 — FITS input

Folders of FITS frames work as input, for capture software that writes it
rather than camera raw — INDI/EKOS, SharpCap, N.I.N.A., FireCapture. Colour
(CFA + `BAYERPAT`), monochrome, and already-debayered 3-plane cubes; the Bayer
pattern is rolled to RGGB the same way a camera raw is, and `XBAYROFF` /
`YBAYROFF` are respected.

**No new dependency.** `astropy.io.fits` is used if installed, then `fitsio`,
then a built-in reader for plain uncompressed FITS — which is what cameras
actually write. astropy earns its place on tile-compressed files and the odd
corners of the standard, so it is an optional extra
(`pip install 'eclipseforgehdr[fits]'`) and the error message says to install it
if a file turns up that the built-in reader cannot open.

A FITS frame arrives without the two things LibRaw normally supplies, so both
are recovered from the file:

- **Exposure** from `EXPTIME` / `EXPOSURE`, and it is required rather than
  guessed — it is what groups frames into tiers. A file without it stops the
  run with that said plainly.
- **The saturation ceiling** from `SATURATE` / `DATAMAX` if the writer records
  one; otherwise from a saturation *plateau* in the data, since a real ceiling
  shows as many pixels sharing the maximum. That test now also requires the
  plateau to be under 20% of the frame — caught by a test of my own where a
  uniform frame tripped it and would have told the merge that a perfectly good
  frame was clipped everywhere. Bit depth is the last fallback. The report says
  which was used, because guessing high is the dangerous direction: clipped
  pixels then merge as if they were valid.

There is no colour matrix or white balance in a FITS header, so both are
identity and the README says so. A colour-camera frame comes out green-dominant
as captured; Warmth, Tint and Neutralise sky cast are the controls. Inventing a
white balance would bury a colour error inside the photometry where nothing can
see it.

Tested against hand-written FITS covering int16-with-BZERO, float32, RGGB and
BGGR (which must decode identically — they do), mono, 3-plane cubes, pedestal
subtraction, and a missing exposure time.

**New defaults**, from the settings the reference bracket was worked to: radial
flatten 0.5, FNRGF share 0.17 / strength 0.9, MGN contrast 0.04, NAFE-VN mix
0.15, Pellett off, short-exposure detail 0.31, glare dim 0.05, warmth 0.9, tint
1.205, neutralise sky cast 1.0, saturation 1.0, highlight compression 0.1,
output gamma 1.0, black point 0.005.

## 0.11.5 — Black-point correction, previously misdiagnosed as vignetting

Setting every structure slider to zero and still seeing it was the decisive
clue: it was never NAFE, or any detail layer. It was the base envelope.

`Bg` normalised the merged luminance with `lo = percentile(lum, 2)`. On a wide
field the sky **is** most of the frame, so that lands the black point within a
noise sigma of the sky itself — and `xn` becomes a small difference between two
nearly equal numbers, which turns a tiny real variation into an enormous one on
screen. Measured on the reference set, corner against mid-edge:

| | corner | mid-edge | ratio | sky clipped to black |
|---|---|---|---|---|
| the data itself | | | **1.030** | |
| `lo = p2` (before) | 0.0478 | 0.0637 | **1.332** | 2.4% |
| `lo = sky − 5σ` (now) | 0.1253 | 0.1284 | 1.025 | 0.0% |

A 3% brightness difference across the frame was being displayed as 33%, and
2.4% of the sky was crushed to pure black. The black point is now measured from
the sky itself — median and MAD beyond 2.5 R — and clamped never to sit above
the 1st percentile, so it can only ever clip less than the old rule, never more.

It also **un-pinned the radial flatten control**. With the sky that close to the
floor, `rprof` sat on its 0.12 clamp everywhere in the outer field — corner and
mid-edge both read exactly 0.1200 — so radial flattening did nothing out there
at any setting. It now reads 0.165–0.170 and works.

End to end, rendered composite with every detail layer at zero: corner-to-edge
ratio **1.332 → 1.000**. With the default layers on, 1.017.

**The background now starts lighter**, because the crush is gone rather than
hidden. `Black point` is the control, and it now behaves uniformly instead of
distorting shape. Measured on the reference set:

| bgBlack | sky | corona | corona/sky |
|---|---|---|---|
| 0.02 (default) | 0.185 | 0.251 | 1.36 |
| 0.08 | 0.131 | 0.202 | 1.54 |
| 0.11 | 0.102 | 0.179 | 1.75 |

About 0.11 reproduces the old background darkness — but note the corona-to-sky
contrast *improves* as you crush, which it could not do before, because the
crush used to take the shape with it.

## 0.11.4 — Sky colour neutralisation; inoperative slider corrected

**Neutralise sky cast had no effect at any setting, and never had.** It divides
the chroma field by the measured background colour — but that colour was being
measured on `ratio`, which carries a confidence fade that drives it to exactly
1.0 wherever the signal is near the noise floor. That is precisely the region it
was measured in: mean confidence there is **0.015**. So it returned R 1.000
G 1.000 B 1.000 on a sky whose real colour is R 0.985 G 1.037 B 0.681, and
dividing by unity does nothing. It is now measured on the HDR itself.

That exposed a second problem. With a real value in hand, dividing the whole
chroma field by it tipped the *far sky* blue — because `ratio` had already
forced that region neutral, so the correction was applied twice there and once
in the corona. The correction is now weighted by the same confidence that built
`ratio`: full where there is real chroma to correct, absent where the chroma was
discarded. Measured across the slider:

| bgNeutral | far sky B/R | corona B/R |
|---|---|---|
| 0.00 | 1.000 | 0.166 |
| 0.50 | 1.000 | 0.200 |
| 1.00 | 1.000 | 0.239 |

The sky stays exactly neutral at every setting while the corona's atmospheric
cast comes off monotonically — which is what the label promises.

**The sky gradient is now removed per channel.** One shared correction flattens
brightness and leaves the colour gradient behind. Fitting each channel takes the
colour with it, and the fitted spans are R 1.202× G 1.262× **B 1.335×** — blue
steepest, which is what Rayleigh scattering with airmass does, and the best
evidence that what is being fitted is real atmosphere. Measured on the sky:

| | before | after |
|---|---|---|
| brightness spread across quadrants | 1.122× | **1.019×** |
| colour spread, R | 0.0551 | **0.0106** |
| colour spread, G | 0.0123 | **0.0024** |
| colour spread, B | 0.0429 | **0.0072** |

**New defaults.** The starting slider positions are now the settings the
reference bracket was worked to by eye, rather than the ones inherited from
before the layers behaved: MGN contrast 0.4 → 0.1, FNRGF strength 1.5 → 0.5,
NAFE-VN mix 0 → 0.2, short-exposure detail 0.24 → 0.6, glare dim 0.35 → 0.15,
radial flatten 0.2 → 0.35, base lift 0.18 → 0.255, grain smoothing 0 → 0.25,
neutralise sky cast 0.7 → 0.5, and a few smaller moves. They are starting
points, not truths — every one is a taste call and every one is a slider.

## 0.11.3 — Second-order sky gradient model

0.11.2 took out the tilt and left a curved remainder, which then read as
darkening on the *other* side and around the corners. Measured on what the
plane left behind, beyond 4.1 R:

| model fitted to the remainder | explains |
|---|---|
| another plane | 0.8% ← the tilt really was gone |
| radial about the frame centre (vignetting) | 2.5% ← still not vignetting |
| **full quadratic in x, y** | **51.8%** |

So the sky curves across a 3.5° field near the horizon, and a plane cannot
follow it. The model is now a full quadratic, still fitted only beyond the
measured corona extent:

| shell | before | after | retained |
|---|---|---|---|
| 1.6–2.4 R | 0.1735 | 0.1608 | 93% ← corona survives |
| 2.4–3.2 R | 0.0810 | 0.0576 | 71% ← corona survives |
| 4.0–6.0 R | 0.0594 | 0.0088 | **15%** ← sky collapses |

Six terms need plenty of sky to be stable, so it falls back to the plane below
100k fitted pixels. A model spanning more than 2× across the frame is rejected
outright — that is not sky. And the full-resolution evaluation is done in row
blocks, since six term arrays at 45 MP would be a gigabyte.

## 0.11.2 — Sky gradient removal

A broad dark gradient across one side of the frame, strongest at low solar
altitude. Measured on the reference set it is **1.20× corner to corner**.

**It is not vignetting.** A radial model about the frame centre — what lens
vignetting looks like — explains **0.0%** of it, against 54% for a plane in x
and y. So flats would not have removed it, and it has to be fitted per run
rather than calibrated once.

**Where you fit it matters more than the model.** Fitted close in, the plane
absorbs the corona's own east–west asymmetry, which is real solar structure.
The fitted amplitude runs 1.22 in the 1.6–2.4 R shell, 0.55 at 2.4–3.2 R, 0.34
at 3.2–4 R and 0.22 beyond 4 R, converging only once the corona has faded into
the sky — while the *direction* is stable at −13° to −18° everywhere, which is
what says it is one real gradient rather than a fitting accident. So the fit is
restricted to beyond the measured `corona_extent_R`, and only what it finds
there is removed:

| shell | before | after | retained |
|---|---|---|---|
| 1.6–2.4 R | 1.2235 | 1.0381 | 85% ← corona, survives |
| 2.4–3.2 R | 0.5530 | 0.3676 | 66% ← corona, survives |
| 3.2–4.0 R | 0.3425 | 0.1572 | 46% |
| 4.0–6.0 R | 0.2184 | 0.0331 | **15%** ← sky, collapses |

The correction is multiplicative and identical in all three channels, so colour
is untouched. It is gated: below 2% amplitude or 8σ it is measured, reported
and left alone. The run report prints what was found either way.

One thing this cost me an iteration: the fit must be in **true log**, not
`log1p(S/median)`. The correction is multiplicative on the linear image, so the
fit has to live in the space where a multiplicative change is an additive one —
and `log1p` compresses exactly where the sky sits. Fitted that way the
correction came out at half strength (sky residual 0.105 → 0.063 instead of
→ 0.017).

## 0.11.1 — Limb rim artefact and composite contrast

Both of these turned out to be the same thing.

E is a rank, so its useful range is set by how much of its neighbourhood a
pixel actually beats. In the quiet corona that is a narrow band around the
middle — hence the flat look — while the near-limb brightness **ridge** is a
local maximum by definition and pins at exactly 1.0 — hence the blown white
rim. Raw, you get a washed-out corona and a saturated rim at the same time, and
they are two faces of one problem: the layer was being used unscaled.

The layer is now rescaled by its own robust spread, with the response rolled
off softly past 3 robust sigmas so an extreme cannot reach the rail. Measured
at full resolution (rim = 1.04–1.18 R):

| | 1.05–1.5 R | 1.5–2 R | 2–2.5 R | 2.5–3 R | 3–4 R | rim p99.9 | rim >0.99 |
|---|---|---|---|---|---|---|---|
| 0.11.0 | 0.0747 | 0.0272 | 0.0268 | 0.0276 | 0.0301 | 1.0000 | 6.75% |
| **0.11.1** | **0.1225** | **0.0882** | **0.0843** | **0.0856** | 0.0940 | **0.8792** | **0.00%** |

Corona contrast 3.1–3.3×, limb 1.6×, and the rim stops clipping entirely.
Repeatability under an independent σ_A is unchanged (+0.96 at the limb, +0.77
to +0.81 outside), so this is a rescale of the same structure, not new
structure invented out of noise.

The knee position was chosen by measurement: at 2 the rim is darkest but the
corona gives up contrast; past 6 the rim starts clipping again; without it the
contrast barely improves *and* the rim returns. `knee=0` returns the raw rank
for anyone who wants the paper's E untouched.

Note the sky looks grainier than in 0.11.0 — everything is stretched three
times, sky included. That is honest, and the mix slider plus the composite's
own envelope decide how much of it reaches the image.

## 0.11.0 — NAFE-VN in the published parameter units

You were right to send the paper. The variable neighbourhood was implemented
correctly — eqs. 11–12 are there, and the ratio this code computes,
`(C(a) − C(a−ε)) / (C(a+ε) − C(a−ε))`, is exactly the restricted rank the paper
defines. This *was* NAFE-VN and not plain NAFE. What was wrong was the units
everything was measured in, and that broke the one thing that makes the method
work on noisy data.

**The bug.** The image was passed through an equal-population rank map first,
and ε and σ were then applied in *rank* units. That looks harmless — a monotone
remap does not change which neighbour is brighter than which. But eq. 13's noise
adaptivity depends on σ being a **fixed width in the value units of A**: where
the local histogram is wide (real contrast) a fixed smoothing is negligible,
where it is narrow (noise only) it dominates and flattens the rank. Under a rank
map the sky — most of the pixels — is stretched to fill most of the axis and the
corona is squeezed into what is left, so the same physical width means something
different at every radius. The automatic behaviour is gone. A per-level
correction had been bolted on to compensate; it did not.

Per-annulus contrast (sd) on the reference frame, full resolution:

| | 1.05–1.5 R | 1.5–2 R | 2–2.5 R | 2.5–3 R | 3–4 R (sky) |
|---|---|---|---|---|---|
| 0.10.4, rank units | 0.0944 | 0.0105 | 0.0253 | 0.0602 | 0.1041 |
| **0.11.0, level units** | 0.0747 | **0.0272** | **0.0268** | 0.0276 | **0.0301** |

In rank units the layer's contrast *rose* with radius, tracking the falling
signal-to-noise — it was a noise detector, which is why the sky was grainy and
the corona flat. In level units it is nearly constant with radius, which is what
a scale-free rank filter should produce: the corona from 1.5 to 3 R gains
40–160% of contrast and the sky loses two thirds of its grain. Repeatability
under an independent σ_A of added noise is +0.99 at the limb and +0.77 to +0.82
everywhere else.

Three other things the re-reading turned up:

- **K = 128, not 64.** The paper notes its images have "several thousand
  discrete pixel values". 64 levels is a real approximation and it cost real
  contrast (2–2.5 R: 0.0273 at K=64 against 0.0366 at K=128). Past 128 the
  corona keeps gaining but the sky gains faster.
- **The paper's plain Gaussian kernel beats the multiscale sum** this had been
  building for the fuzzy weights l_{k,l} of eqs. 4–5 — more corona contrast and
  nearly twice as fast.
- **The speed is legitimate.** The paper calls the naive per-pixel algorithm
  "extremely time consuming"; this evaluates every pixel's local histogram at
  once as K blurred membership maps, which is the same computation reorganised
  rather than an approximation. The one genuine shortcut is sampling that field
  every `grid` px — and grid 8 against grid 4 agrees to three decimals in every
  annulus, so it is free. K=64 was the shortcut that was not free.

## 0.10.4 — Revert of 0.10.3

**0.10.3 broke the NAFE layer and should not be used.** Reverted; this build is
0.10.2 behaviour.

Weighting the rank map toward r < 1.5 R starved the *mid* corona, which is
neither in the favoured region nor negligible. Measured at full resolution on
the reference set, the field from 1.5 to 3 R collapsed to sd 0.044 with p1–p99
of 0.46–0.75 — a flat grey plateau with a hard edge where the weighting ran
out — while inside 1.5 R it clipped at 1.000.

It should not have shipped, and the reason it did is worth recording:

- The metric was **high-pass sd**, which only sees pixel-scale variation. The
  streamers are mid-scale, so the collapse was invisible to it. The right
  metric is per-annulus contrast, and by that measure 0.10.3 fails immediately.
- The check image was **cropped to 2.4 R**, just inside the plateau's edge.

Both checks were run on a decimated probe rather than the full-resolution code
path, which is also where the σ_A error in 0.10.3 hid. Any future attempt here
must be judged on per-annulus contrast, out to at least 4 R, on a picture of
the whole frame, computed on the shipped path.

Verified after the revert: 0.10.4 matches 0.10.1 to within 0.002 sd in every
annulus from 1.05 R to 4 R.

The 0.10.2 changes (ε, noise σ, a working `sigma_sp`) are kept — they were
measured on the same footing and are small.

## 0.10.3 — NAFE histogram range (reverted in 0.10.4)

Testing the "crop the field" hypothesis from 0.10.2 turned into a real fix.

The rank map decides how finely NAFE can resolve brightness, and it was built
from the whole frame. A 600 mm full-frame bracket is 6.5 x 4.4 lunar diameters
and overwhelmingly sky, so the map spent its resolution on sky noise and had
almost none left for the inner corona — exactly where a rank filter should be
strongest.

Cropping the working field to 2 lunar diameters did triple the corona detail,
confirming the diagnosis. But cropping throws away real corona, and a hard mask
on the statistics was worse still: it pinned 61% of the frame at rank 0 or 1 and
blew the outer corona out. The answer is a **weighted** rank map — pixels inside
1.5 R count fully, everything outside counts at 0.10 — so the map still spans
the whole value range and nothing clamps.

Measured on the reference set. "repeat" is the correlation between two runs of
the same image differing by one independent σ_A of added noise, which separates
reproducible structure from amplified noise:

| | inner detail (1.05–1.8 R) | repeat | outer detail (2.2–3.2 R) | repeat |
|---|---|---|---|---|
| 0.10.1 | 0.00559 | +0.906 | 0.04850 | +0.663 |
| 0.10.3 | **0.01552** | **+0.961** | 0.04246 | +0.656 |

Inner-corona detail **2.8×**, and it is *more* reproducible than before, not
less — so it is structure, not amplified grain. The outer corona pays 12%.

Two things this test also caught:

- σ_A must still be measured over the **whole** frame. Measuring it inside
  1.5 R, where the radial gradient is steep, made consecutive-pixel differences
  signal-dominated, σ_A came out too large, the level smoothing over-fired and
  the outer corona lost 38% of its detail for nothing.
- The FNRGF-correlation metric used in 0.10.2 is contaminated: the two layers
  share a noise realisation, so part of that agreement was agreement on noise.
  The added-noise repeatability test replaces it.

## 0.10.2 — NAFE parameters derived from measurement

Two of the three NAFE constants had never been swept. Both moved, both measured
on the reference set (corona = 1.05–2.2 R, sky = beyond 3.2 R):

- **ε 0.05 → 0.10.** The value window was in the steep part of its curve. At
  0.02 it prints concentric rings — the paper's "fragmentation" — at 0.05 the
  high-pass structure is 0.076, at 0.10 it is 0.109 with the agreement against
  FNRGF unchanged (0.711 → 0.721), and past 0.15 the curve is flat and the
  extra is noise.
- **Noise σ: 2 σ_A → 4 σ_A.** The paper's range is 2–12 and this sat at the
  bottom of it, with a cap that clipped even that.

  | σ / σ_A | corona detail | sky grain | detail/grain |
  |---|---|---|---|
  | 2 | 0.00843 | 0.12421 | 0.068 |
  | **4** | **0.00820** | **0.06621** | **0.124** |
  | 6 | 0.00782 | 0.04488 | 0.174 |
  | 12 | 0.00696 | 0.03085 | 0.226 |

  Sky grain halves from 2 to 4 for 3% of the corona detail. It keeps falling
  after that, but a broad dark halo appears around the disc as the level
  smoothing lets the radial gradient back in, and by 12 the halo has swallowed
  the streamers. 4 is the last value with no visible cost.

- **`sigma_sp` now does something.** On the fuzzy path it was silently ignored
  and the neighbourhood was fixed in grid pixels — i.e. fixed relative to the
  decimation factor rather than to the Sun. It is now the widest scale of the
  multiscale kernel and defaults to 0.13 R. Measured flat from 0.06 R to
  0.51 R, so the value is not critical; scaling with the disc is.

End to end on the reference set, detail-to-sky-grain improves 0.094 → 0.124.
That is real but modest, and it does not close the gap to the paper's Fig. 5 —
see the note in `nafe.py` on why.

## 0.10.1 — NAFE output term corrected

- **The NAFE layer was 99% a gamma stretch.** The 2014 paper's output is
  `B = (1−w)·T_γ(A) + w·E` (eq. 2) with w in 0.05–0.3, and that is what was
  being stored — but B is their *final display image*: T_γ carries the
  large-scale brightness and E carries the structure, so at w = 0.2 the result
  is four fifths gamma transform by construction. As a detail layer that was a
  second copy of the base image, and mixing it in diluted MGN and FNRGF rather
  than adding to them. It looked, correctly, like a greyscale version of the
  plain stack.
  The layer is now **E**, the equalized field. Measured on the reference set
  (×4 decimated, K=64, γ=2.4, ε=0.05):

  | | stored before (B) | stored now (E) |
  |---|---|---|
  | correlation with a plain gamma stretch | 0.992 | 0.562 |
  | high-pass structure (sd) | 0.0180 | 0.0675 |

  Equation 2 has not been abandoned — it now happens where it belongs, one
  level up: the composite's envelope is T_γ and the **nafeMix slider is w**.
  Since the paper's w runs 0.05–0.3, the useful part of that slider is its
  bottom third. `nafe_vn(..., combine=True)` still returns the paper's B for
  anyone who wants the standalone image.

## 0.10.0 — Full code audit

A line-by-line review of every module, with a specific eye on datasets other
than the two this was built against. Nothing here changes how the reference
render looks except where noted.

**Wrong pixels**

- **The clipping test was made after white balance.** `cmax` was the maximum of
  the *white-balanced* channels but was compared against the raw saturation
  level. White balance runs at G=1, so the ~2.1× red gain pushed a red pixel
  past the clipping threshold while its photosite sat at only 0.45 of
  saturation: every tier declared prominence H-alpha clipped a full stop early
  — the shortest tier included — and with no tier holding weight there, the
  merge filled prominence cores by leakage from tens of pixels away. Neutral
  pixels are unaffected either way (for a neutral subject green is the largest
  raw channel and white balance leaves it alone), which is why this never
  showed on the corona itself. **This is the most likely reason prominences
  have been hard to recover.**
- **`long_lum.npy` was saved in the uncropped frame** while everything else was
  cropped, so the earthshine model could not even broadcast against it: an
  outright crash at 99% on any run with `earthshine` on and a non-zero
  alignment trim.
- **`prom_rgb.npy` was saved in the uncropped frame too.** The prominence gate
  is upsampled ×2 and laid over the cropped composite, so every prominence was
  painted `crop_origin` pixels off the limb. A few pixels on the reference set;
  the whole gate off the limb on a set with real drift.
- **`prom_geom` was written to `geometry.json` after the layers were built**, so
  a first run always fell back to the merged limb for the prominence gate and
  only a *rebuild* of the same folder used the prominence tier's own limb. The
  same folder produced two different gates depending on whether the cache was
  warm.
- **The autocrop was skipped on the render when only the bottom and right
  needed trimming** (`crop_origin` stays `(0,0)` in that case). The tier TIFFs
  were being sliced, the render was not — different sizes, and the render kept
  a band of edge-replicated pixels for MGN to find gradients in.

**Preview vs export**

The browser preview *is* the preview — there is no server-side preview render —
and its JavaScript had drifted from the numpy it mirrors. Five divergences, all
now resolved in favour of the numpy:

- `discTrim` was applied in preview pixels against a mask radius that was
  already decimated, so the slider moved the mask **4× too far** in the preview
  and correctly in the export.
- The FNRGF blend ramp was linear over 1.05–1.40 R in the preview and a
  smootherstep over 1.02–1.57 R in the export (0.43 vs 0.25 of the mix at
  1.2 R).
- The inner-corona weight was linear over 1.38 R/0.33 R against a smoothstep
  over 1.45 R/0.40 R.
- The prominence lift used `2·(inner−0.5)` in the preview against
  `0.30 + 1.3·(inner−0.5)` in the export — the "presence" bias existed only on
  export, so the preview showed no lift at all at neutral.
- The chroma ratio was clipped to [0, 2.5] in the preview and [0.2, 3.0] in the
  export, so deep prominence red looked more saturated in the preview than it
  came out.
- The NAFE median was sampled with a stride of 7 over an RGBA byte array, so one
  sample in four was the alpha byte (always 255) and the "median" was really the
  ~67th percentile — the NAFE-mixed field rendered darker in the preview than in
  the export.

**Colour management**

- **The composite export carried no ICC profile.** It is already in display
  encoding, so Photoshop guessed right by accident; PixInsight assumes linear
  and got it wrong. 16-bit TIFF, 8-bit TIFF and JPEG exports now embed the sRGB
  profile.
- **Black level is now subtracted per channel** rather than as the mean of the
  four. On bodies that report four equal values this is a no-op; where they
  genuinely differ, the leftover offset was being multiplied by the white
  balance gain and landing on the outer corona, which sits only tens of ADU
  above black.

**Constants that encoded the reference dataset**

Each of these reproduces the reference numbers to within a fraction of a pixel
and now means the same thing at any plate scale:

- The prominence NCC patch and search radius were 60/80 absolute pixels —
  0.19 R/0.26 R on the reference disc, but **0.65–1.45 R on a 300 px disc**,
  where the patch is dominated by the lunar edge, the one feature that *moves*
  between tiers. Worse, every anchor then contains the same limb, so they all
  agree and the spread gate cannot see it: the tiers get registered to the Moon
  instead of the corona. Now fractions of the measured radius.
- The prominence gate window was 4/60/8/6/70/25 half-res pixels — 0.19 R above
  the limb on the reference set, **0.4 R on a 300 px disc**, where the gate
  covers inner-corona loops and `promGain` brightens them as if they were
  prominences. Now fractions of the measured radius.
- The near-limb saturation feather was 20 pixels — 3.2% of the reference lunar
  radius, a tenth of a 200 px one. Now tied to the measured radius.
- The short-exposure inner stack took `secs[:4]` — a tier *count*. That spans
  4.3 EV on the 14-tier reference bracket but 6.4 EV on a 6-tier one, and since
  the weight is proportional to exposure time the 4th tier would carry ~8× the
  weight of the 1st, making this layer a blown inner corona instead of the crisp
  one it exists to be. Now capped at 24× the shortest exposure, which selects
  exactly `secs[:4]` on the reference set.

**Crashes and silent failures**

- `fit_limb` could produce a negative radicand (true limb outside its search
  band), hence a NaN radius, hence `np.arange(nan, nan, 0.5)` raising on the
  next iteration — killing the run at 91% instead of reaching the plausibility
  check written to handle exactly that.
- The alignment-network residual reduced over an empty array when every frame
  shares one shutter speed, and it only ever measured the y axis: a network
  perfect in y and 15 px inconsistent in x reported 0.
- `sat_level` is now cross-checked against the data. LibRaw occasionally reports
  a white level of 0, an already-black-subtracted value, or a 12-bit ceiling for
  14-bit data; any of those made the merge treat the whole inner corona as
  clipped in *every* tier, and the result was a black frame — cached, with no
  exception raised.
- Frames of differing size in one tier (a crop-mode shot, a stray file from
  another body) raised a bare broadcast error from the hot-pixel mapper. Now a
  message that says what is wrong.
- `find_prominences` computed its detection sigma over an array where
  out-of-frame azimuths had been filled with exact zeros. With the disc near a
  frame edge those zeros dragged the MAD toward zero, the threshold followed,
  and every in-frame azimuth was flagged as a prominence — handing the matcher
  arbitrary points on the limb. Now measured over the in-frame azimuths only.
- `stack_variance` was handed a cropped-frame centre for uncropped tiers, so on
  a run with real drift the "limb variance" and "rim" numbers in the report
  described mid-corona. Diagnostic only, but it is the number used to judge a
  run.
- The intra-tier motion warning compared full-res displacements against a
  half-res window, so it fired at half the intended motion.

**Cache and session**

- The cache was validated on the render options alone. Adding a frame to the
  folder, deleting a bad one or replacing a file and pressing Start (without
  force) silently reused the previous stack, and the report exported alongside
  it quoted the previous frame count. The raw file list, sizes and mtimes are
  now part of the key.
- A run that dies partway left a valid-looking `opts.json` beside a mix of new
  and old products; a later non-forced Start would declare that mixture valid.
  The options file is now cleared before a run rather than after.
- The previous result stayed live and exportable for the whole of a re-run, so a
  failed run still reported ready and exported the *old* layers under the new
  settings.
- The run thread re-read the selected folder at execution time, so picking
  another folder mid-run processed the new one and loaded the old one's layers.
  The folder is now bound when the run starts, and switching folders during a
  run is refused.

**Known and deliberately not changed** (each would move the reference render,
and each is documented in place): MGN's `global_wt` is 0.12 where Morgan &
Druckmüller Eq. 5 uses h = 0.7; MGN's scales and the FNRGF smoothing lengths are
in absolute pixels; `nafe.py`'s `sigma_sp` argument is a no-op in the default
fuzzy path; and the radial flattening that precedes MGN/FNRGF starts at R+3
rather than at the disc mask radius, so a few rings are fitted on the limb ramp.

## 0.9.10

- **Glare dim no longer prints a ring.** It shared its radial weight with the
  short-exposure detail layer, whose smoothstep window closes at 1.45 R. At full
  strength that is a 3.3× brightness ramp ending at a definite radius, and a
  definite radius is a visible circle. Glare is a broad wing off the limb, not
  something that stops somewhere, so it now has its own exponential profile
  (scale length 0.6 R) with no boundary at all. Measured outside the limb
  (r > 1.1 R): peak curvature of the multiplier drops from 21.5 at r = 1.42 R
  to 1.54, peak slope from 2.61 to 0.97.

## 0.9.9

- Tier TIFF export now carries an **embedded ICC profile** — sRGB-encoded by
  default (opens looking correct in Photoshop and Affinity), or scene-linear
  with a linear profile via the checkbox, for PixInsight. Untagged linear data
  was being decoded as sRGB by the host, which is what made the exports look
  harsh with clipped histograms.
- Fixed the export normaliser: it divided by sensor saturation *after* white
  balance and the colour matrix had already pushed red and blue above it, so
  strongly coloured highlights the sensor had not clipped were clipped in the
  file. The linear export now divides by the measured headroom instead.
- A `README.txt` is written alongside the tier TIFFs explaining the scale and
  why the set must not be mean-stacked as if it were one exposure.

## 0.9.8

- **Fixed the alignment regression introduced in 0.9.5.** `find_center` was
  being used as the per-tier disc locator for cross-tier pre-centring; it
  returns the centroid of the brightest 0.05% of pixels, which drifts **197 px**
  with exposure on tiers that are genuinely 22 px apart. The coarse centre now
  comes from the limb fit, and pre-centring is gated on the measured spread
  (> 15% of the correlation window) so small offsets keep the shared fixed
  window. Measured: fixed window 0.58 px, limb pre-centred 1.86 px,
  `find_center` pre-centred 20.01 px.

## 0.9.5

- Autocrop of the alignment border, computed before the merge from the recorded
  per-tier edge vacancies and the cross-tier shifts, applied to the render and
  to each tier TIFF; even origin so the Bayer phase survives, and a 60% guard
  against trimming most of the frame.
- Large frame-to-frame offsets: pre-centring each tier's correlation window on
  its own disc position, for sets where the disc moves hundreds of pixels
  between brackets. (This is the change 0.9.8 had to gate — it was right for the
  Sony set and wrong for the Lumix one.)
- Credit line changed to "© 2026 N. Augustin with ClaudeAI".

## 0.9.4

- Renamed from CoronaForge to **EclipseForgeHDR** (branding collisions), with
  the subtitle "High-Dynamic-Range Solar Eclipse Image Processing". Version
  number added to the build filename, `efhdr` short alias added, and the version
  made dynamic from `__init__.py` — it had been written in two places and the
  two had drifted (`pipx list` said 0.7.1 while the code said 0.8.0).
- **Export of the aligned exposure tiers** as 16-bit TIFFs, for stacking or
  blending in PixInsight, Photoshop or Affinity.

## 0.9.3

- Saturation-based rejection of non-totality frames: a frame whose saturated
  area exceeds three times the tier median is dropped. Added after a Sony set
  gained partial-phase frames with a bright crescent, which pushed the fitted
  lunar radius from 451 px to 857.
- Cross-tier radius consensus: the median of the per-tier limb radii, used to
  override the merged fit when it deviates by more than 15%.
- The disc mask margin is derived from the measured limb ramp rather than a
  fixed number.
- Photometric calibration regions are now selected by signal-to-noise against
  each tier's own sky noise, replacing a "brightest 20% of the frame" rule that
  measured noise on a short bracket and rejected every link.

## 0.9.2

- **MGN halo fixed.** Two wrong hypotheses first: moving the background start
  radius made it worse (0.0331 → 0.0615), and reflecting the layer across the
  mask edge gained 7%. The actual cause was `_deband` subtracting an azimuthal
  *mean*, which leaves the rim's variation *around* the disc untouched — that
  residual measured 0.0140 before and after. Now debanded at order 6, chosen
  from a measured order-vs-rim/streamer trade recorded in the code.

## 0.9.1

- NAFE-VN added as a mixable detail layer with its own view and slider, after
  three failed attempts: 360 full-resolution blurs (too slow, fixed with a
  coarse histogram grid), a flat output (fixed with equal-population rank
  binning), and contour rings (fixed with Gaussian membership and a noise sigma
  correctly carried through the rank map per level).
- Chroma reconstruction gained a noise-relative floor and a confidence fade, so
  colour is not invented where there is no signal.
- The radial profile was rewritten by sorting rather than binning: 21× faster
  and no 6000 px cap.

## 0.9.0

- Prominence anchors: prominences are located automatically on a reference tier
  and matched between tiers by normalised cross-correlation, adding independent
  hard links to the alignment network. Two rewrites — the first two attempts
  found coronal streamers instead of prominences, because a detector based on
  "how far out does this azimuth stay bright" finds exactly the wrong thing.
- Alignment quality metrics: per-tier variance at the limb and in the corona,
  the disagreement rim, and the limb 20–80% transition width, all in the report.
- Value-based neighbourhood replacing geometric masking in the equalisation.

## 0.8.8

- **The cross-tier shift sign was inverted.** `phase_cross_correlation(ref,
  mov)` returns the shift to apply to `mov`; every one of the six consumers
  negated it, which does not merely fail to align — it leaves a residual of
  *twice* the true offset. It went unnoticed because all six negated it
  consistently, so the moon track, the moon-mask trial and the alignment quality
  numbers were all computed in the same wrong frame and agreed with each other.
  Measured on the reference set after flipping to `+shift`: moon-track scatter
  about the straight line 23.5 → 2.0 px, tier-to-tier variance at the limb
  0.377 → 0.073, in the corona 0.165 → 0.136, merged limb 20–80% transition
  23.0 → 15.5 px. This one bug explains the "125 px lunar spread", the 0.89 px/s
  apparent drift against a physical 0.35, the 74 px inner-stack offset and the
  0.8.2 mask failure below.

## 0.8.2

- **Shipped a merge with two Moons in it.** The per-tier lunar mask was gated on
  "is the limb sharper" and "are there holes" — and a fragmented disc scores
  *sharper*. Fixed by constraining the mask centres to a robust straight-line
  lunar track (the Moon's motion in a corona-aligned frame is orbital and must
  lie on a line) and adding a circle-fit-rms gate: 1.3× for good masking, 3.5×
  for the broken version.

## 0.8.0 – 0.8.6

- Redundant-link alignment network: lag-1 *and* lag-2 links between tiers,
  solved globally by weighted least squares, so no single bad pair can drag the
  chain.
- Robust straight-line fitting of the lunar track across all tiers, used to
  predict each tier's disc position and radius.
- Photometric cross-calibration rewritten with Huber IRLS.
- Hot-pixel mapping across frames with a voting threshold.

## 0.7.2 – 0.7.9

- Sky-cast neutralisation, warmth/tint separation and luminance-preserving
  colour moves consolidated.
- TIFF bracket input.
- Report generation.

## Changes in 0.9.10

- **Glare dim no longer prints a ring.** It shared its radial weight with the
  short-exposure detail layer, whose smoothstep window closes at 1.45 R. At
  full strength that is a 3.3x brightness ramp ending at a definite radius,
  and a definite radius is a visible circle. Glare is a broad wing off the
  limb, not something that stops somewhere, so it now has its own exponential
  profile (scale length 0.6 R) with no boundary at all. Measured outside the
  limb (r > 1.1 R): peak curvature of the multiplier drops from 21.5 at
  r = 1.42 R to 1.54, and peak slope from 2.61 to 0.97.

## Changes in 0.9.9

- Tier TIFF export now carries an **embedded ICC profile** — sRGB-encoded by
  default (opens looking correct in Photoshop and Affinity), or scene-linear
  with a linear profile via the checkbox, for PixInsight. Untagged linear data
  was being decoded as sRGB by the host, which is what made the exports look
  harsh with clipped histograms.
- Fixed the export normaliser: it divided by sensor saturation *after* white
  balance and the colour matrix had already pushed red and blue above it, so
  strongly coloured highlights the sensor had not clipped were clipped in the
  file. The linear export now divides by the measured headroom instead.
- A `README.txt` is written alongside the tier TIFFs explaining the scale and
  why the set must not be mean-stacked as if it were one exposure.

## Changes in 0.7.1

- **Reverted the per-tier lunar masking added in 0.7.0.** It destroyed the merge.
  The idea depends on every tier's limb fit being trustworthy, and the shortest
  tiers are noise-dominated: on the Lumix run their fits scattered by 124 px, so
  the union of exclusion discs ate most of the inner corona. The limb fit then
  failed on the wrecked image, fell back to the gradient fit (R=1376 against a
  true 625), and the disc mask became a giant circle. The Moon cannot actually
  move far in a corona-aligned frame across one bracket, so a large measured
  spread means bad fits rather than real motion — either way, acting on it is
  unsafe. The spread is still measured and reported; it is no longer applied.
- The gradient fallback now refuses an implausible radius instead of writing it
  into geometry.json, and says loudly in the log when the good fit failed.
- The limb-profile fix from 0.7.0 is kept — that one was right.

## Changes in 0.7.0

- **The limb profile was carrying a spurious 118 px sinusoid — this is what put
  the disc mask out of register.** `fit_limb_rays` iterates: it measures the
  limb about the current centre, solves for the correction, then applies it. The
  per-azimuth profile was being built from the LAST measurement pass, i.e. in
  the pre-correction frame, so it kept whatever centre offset that pass still
  had. On the Lumix run the seed was bad (the gradient fit said R=1319 against a
  true 625, so the first search window did not even contain the limb) and the
  final correction was still ~59 px — which appeared in the stored profile as a
  sinusoid of 118 px peak-to-peak. The mask then behaved like a circle offset by
  59 px along the up-left/down-right axis: eating the chromosphere on one side,
  leaving dark limb on the other.
  The fit now runs an extra pass AT the converged centre and builds the profile
  there. Measured after the fix: profile swing 20 px (real non-circularity), and
  the residual sinusoid in the profile is **0.65 px**, down from ~59. The search
  window also opens wide on the first iteration and the fit runs 5 passes, so it
  converges to the identical answer from a seed 70 px and 2x off.
- **Each tier is now masked to its own lunar disc during the merge.** The corona
  alignment registers the Sun; the Moon is a different body and moves against it
  during the bracket, so blending tiers whose Moon sits elsewhere smears the
  merged limb. Every tier's limb is fitted in the aligned frame and that tier
  contributes only where it actually sees corona. The run logs the total lunar
  motion across the bracket.

## Changes in 0.6.10

- **Sky-cast neutralisation** ("Neutralise sky cast", default 0.7). At low sun
  altitude, extinction crushes blue: measured on the Lumix set the far sky is
  R 0.95 / G 1.05 / B 0.70 — a real yellow-green, atmosphere rather than corona.
  The background chroma is now measured automatically beyond 0.72 of the frame
  radius and divided out at the chosen strength. The corona's own colour is far
  from the sky's (R 1.66 / G 0.88 / B 0.23), so it survives essentially intact:
  the sky goes neutral while the orange stays. At 0 nothing is applied, at 1.0
  the sky is fully neutral, and the slider runs to 1.2 if you want to push past.
- **Tint slider** (green/magenta), the axis Warmth could not reach — Warmth
  trades R against B, so a green cast was previously uncorrectable.
- Colour moves now renormalise to unit luminance, so Warmth, Tint, Saturation
  and neutralisation change colour without also changing brightness. Saturation
  in particular no longer darkens as you raise it.

## Changes in 0.6.9

- **The limb is now fitted on the merged luminance, not on short_lum.** The Moon
  moves against the corona between tiers, so the four shortest tiers' blended
  limb sits far from the merged limb the composite actually displays — measured
  on the Lumix set: short_lum 2758.5, prominence tier 2707.0, merged 2687.6, a
  71 px spread in y. The fit was correct; it was fitting the wrong image. Masks
  now follow the image being masked.
- **The prominence gate uses the prominence tier's own limb**, fitted separately
  and stored as `prom_geom`, because that stack is a single tier whose Moon sits
  ~20 px from the merged one. The run logs the offset.
- **MGN concentric circles fixed.** The Fourier background dropped its azimuthal
  order in integer steps as rings left the frame; each step moved mu
  discontinuously between adjacent rings and painted a hard circle into the
  output. Harmonics are now damped continuously with coverage and the model is
  smoothed harder along radius where coverage is poor. Ring-step rms drops to
  0.0014 across the whole field, inside and outside 4R.
- **Composite preview fixed.** `L.rmask` is a per-pixel array but was indexed
  with the RGBA byte offset, so three quarters of the frame went NaN and
  rendered black — the "strip" preview. Detail views were unaffected, which is
  why only the composite looked broken.
- **Frame selection is now a choice** (Frames: all / best half / best only).
  All frames per tier averaged (default) gains sqrt(N) in signal-to-noise;
  best-only keeps maximum sharpness. Each tier now logs its sharpness spread so
  the trade-off is visible — on the Lumix set it is 1.06-1.22 for most tiers,
  meaning the frames are near-identical in sharpness and averaging is close to
  free.
- Photometric links that disagree with their own exposure time by more than
  2.5x are rejected and the shutter speed is trusted instead. On very short
  tiers the comparison region is read noise, not signal, which is where the
  16.98x factor came from.

## Changes in 0.6.8

- **Limb fit no longer discarded when it was right.** 0.6.5 accepted the new
  half-level fit only if it agreed with the old gradient fit to within 25% in R
  — i.e. it gated the reliable estimator on the unreliable one. On the Lumix run
  the seed was bad enough to trip that guard, so the stored centre fell back to
  the gradient fit and sat 70 px off in y: the disc mask missed the top of the
  limb entirely, which is the dark crescent between the arrows. The fit now runs
  from two independent seeds, is judged only on its own residuals (rms, rays
  kept, plausible radius), and the better of the two wins. Verified to converge
  to the same answer from a deliberately 70 px-wrong seed.
- **The disc mask follows the measured limb per azimuth.** geometry.json now
  stores the 720-point limb profile r(theta), and the disc mask, the detail-layer
  masks and the prominence gate's radial window all use it instead of a circle.
  The real limb wanders ~12 px peak-to-peak about the best circle, so a circular
  mask had to be either too big (eating the chromosphere) or too small (leaving
  the bow). The margin is now ~2 px instead of ~9.
- **MGN rebuilt — the earlier fix was cosmetic and this one is not.** MGN was
  being flattened by an azimuthally-AVERAGED radial profile. The corona is
  several times brighter on one side, so that left a large residual gradient at
  every azimuth, which inflated the per-scale local sigma and crushed exactly
  the fine structure MGN exists to show. It is now flattened by a low-order
  (order 2) Fourier background mu(r,theta) — the local envelope, too smooth to
  absorb streamers — and normalized on a span matched to the residual's own
  scale rather than the full luminance range. Measured: near-limb detail
  contrast +52%, structure-to-noise 0.85 -> 0.94, no rim, no halo.
  **Your saved mgnContrast will read much stronger — start lower.**

## Changes in 0.6.7

- **Run report.** Every run now ends with a summary: the exposure stack (each
  tier's exposure, frame count, chosen best frame, sharpness score, photometric
  factor and alignment shift), total EV span and integration time, alignment
  residual, sensor defects repaired, the measured lunar limb and its fit rms,
  an approximate plate scale and field of view derived from the lunar disc,
  how far out the corona was traced before it drops into the sky noise, the
  coronal brightness range, the prominence-gate thresholds, the processing
  options in force, and a list of every method used with its citation. It is
  printed into the progress log, written to `.eclipseforgehdr/report.txt` (plus
  `report.json` for machine use), and written again as
  `<name>_report.txt` beside every export — that copy also records the exact
  slider values used for that render.

## Changes in 0.6.6

- **MGN and FNRGF can be switched off.** "MGN contrast" and "FNRGF strength"
  now go down to 0 (both were floored at 0.4; 0 divided by zero in the FNRGF
  compressor, which is why). At 0 each layer collapses to a flat 0.5 and drops
  out of the mix entirely. "Disc mask trim" now spans +/-40 px instead of +/-10.
- **FNRGF ring/block artifacts fixed.** Rings that leave the frame were still
  being fitted with all 13 Fourier coefficients over whatever short arc
  remained, and the hard sigma-clip flipped between neighbouring rings whenever
  a bright streamer crossed the threshold — each flip stepping the background
  model. The order is now matched to the covered arc, the higher harmonics are
  ridge-damped, rejection is a soft Huber weight instead of a hard mask, and the
  model is smoothed hard along radius (harder further out and wherever coverage
  is poor). The hard-edged bright block beside the west streamer is gone.
- **MGN rim fixed.** Any residual azimuthally-averaged radial trend is now
  subtracted from the finished detail layers — by construction a normalized
  detail layer has none, so whatever survives is filter residue. The bright
  spike a few px outside the limb drops from 0.60 to 0.54 (background 0.536),
  and the residual radial trend over R..4R falls 3x.
- **Dark bow at the trailing limb.** The real limb is not a circle — lunar
  relief, seeing, and the moon's drift against the corona between tiers move it
  a few px about the fitted circle (measured: 3.6 px in y, 4.3 px in x across
  the tier sequence, half-res). Masks keyed to the mean radius leave a crescent
  of real limb visible wherever the true edge runs large. geometry.json now
  carries `Rmask = R + 2.5*rms`, used for the disc mask and the detail-layer
  masks, while `R` stays the true limb for the prominence gate.

## Changes in 0.6.5

- **Limb fit rewritten — this was the real cause of the missed prominences.**
  The old estimator took the maximum of the raw radial brightness gradient
  along each ray. That maximum lies inside the bright inner corona, not at the
  limb, and lies further out where the corona is brighter — so R came out too
  large and the centre was dragged toward the bright side. On the Lumix set it
  stored centre (2686.0, 3959.2) R=644.0 where the true limb is
  (2686.9, 3975.0) R=625.2: up to ~35 px of error on one side. Every mask keyed
  to R inherited it, and the prominence gate's radial window sat outside the
  limb on that side, cutting off exactly the bases of the two biggest
  prominences. It now finds the 50% crossing between the disc level and the
  near-limb corona level along each ray — normalized per ray, so azimuthal
  brightness cannot bias it — and fits r(t) = R + dx*cos t + dy*sin t robustly,
  re-centring over three iterations. Measured on the real data: 3.5 px rms over
  720/720 rays. The run logs the rms and ray count, and falls back to the old
  fit if the new one is implausible.
- The diamond-ring/contact frame is fitted the same way, so the ring overlay
  should now land on the limb without manual nudging.

## Changes in 0.6.4

- **Hot/dead pixel repair** (on by default, "Fix hot pixels" next to Denoise).
  Sensor defects sit at the same photosite in every frame, so they are mapped
  once on the *shortest* exposure tier — where the sky is essentially black and
  real sky objects are far below the noise, so a star or planet can never be
  mistaken for a defect — then repaired in every frame of every tier by the
  median of the same-colour neighbours. Outliers are judged against a fitted
  photon+read noise model, not a fixed threshold. The run logs how many
  photosites were found and refuses to repair if the count is implausible.
  Dust *shadows* are not addressed by this: they need flat fields.
- Evaluated ACHF-style radius-adaptive kernel scaling in MGN and did **not**
  ship it — measured on real data it changed the structure-to-noise ratio by
  about 2% (0.85 -> 0.87) and was visually indistinguishable. MGN's per-scale
  local-sigma normalization plus the photon-noise floor already supply the
  adaptivity that ACHF's variable kernel provides.

## Changes in 0.6.3

- **Earthshine is now opt-in and off by default.** There is an Earthshine
  checkbox next to Denoise; leave it unticked and the long-exposure stack and
  the earthshine layer are never built at all (faster run, ~350 MB less cache),
  and the Earthshine slider is hidden. Tick it before Start if a dataset does
  have long tiers with real headroom over the scattered glare. Changing the
  checkbox invalidates the cache, so the next Start rebuilds.

## Changes in 0.6.2

- Fixed an IndexError in the earthshine layer on sensors whose height or width
  is not a multiple of 8 (e.g. 3708 rows). Binning/upsampling steps are now
  size-safe throughout; the prominence gate had the same latent bug for odd
  dimensions.
- The pipeline now warns when a bracket spans less than 6 EV (a totality corona
  bracket normally spans 10-14 EV) and when tiers disagree photometrically far
  beyond their exposure ratio, which points at clipping, cloud, or wrong
  exposure metadata rather than at the merge.

## Changes in 0.6.1

- **MGN inner ring fixed.** The bright band hugging the lunar limb was the
  corona's real brightness peak being re-sharpened, plus the flat occulted disc
  bleeding into the wide MGN kernels. MGN now (a) excludes the disc from its
  local statistics by normalized convolution and (b) subtracts the azimuthal
  radial brightness profile before normalizing, on the same contrast span as
  before. Detail now runs continuously from the limb outwards. The same fix is
  applied to the inner-corona layers.
- **Prominence detection fixed.** The colour stack was written as float16;
  its photometric values exceed the float16 limit (65504), so most of it became
  inf and the redness test returned NaN — hence the near-empty gate. It is now
  float32, built from a **single** fast tier (≤ 1/100 s, the moon drifts between
  tiers so averaging smeared the limb), and the threshold is derived from the
  robust spread of the corona's own colour instead of a fixed multiplier.
- Prominence boost now carries a small positive bias, so a detected prominence
  gains presence rather than only texture.
- The layer cache is invalidated when the app version changes, so an old
  workdir is rebuilt automatically instead of silently reusing broken layers.

## Method references

- Alignment: phase correlation on gradient-flattened log corona
  (Druckmüller 2009, ApJ 706, 1605)
- MGN: Morgan & Druckmüller 2014, Sol. Phys. 289, 2945
- NRGF/FNRGF: Morgan, Habbal & Woo 2006; Druckmüllerová, Morgan &
  Druckmüller 2011, ApJS 194, 25
- Earthshine visualization: inverted-disc local-contrast normalization
  (after Adam Block's HDRMT technique), implemented as a glare-model
  subtraction plus multiscale normalization of the longest exposures.

### Alignment (v0.8.0)

Cross-tier registration uses two independent sources of information, solved
together in one weighted least-squares network:

* **Corona phase correlation** on gradient-flattened log luminance
  (Druckmuller 2009, ApJ 706, 1605), lag-1 and lag-2 links.
* **Prominence anchors.** The fastest tiers contain almost no corona to
  correlate on, so they are also tied in by normalized cross-correlation of
  prominence patches. Prominences are solar features; the lunar limb is not,
  and drifts against the corona during totality
  (cf. MNRAS 503, 5715, 2021: the lunar edge "cannot be used as a reference
  feature for precise registration because it is dynamic during the TSE").
  A prominence link is used only while its anchors agree with each other.

Alignment quality is measured and reported: tier-to-tier coefficient of
variation in the limb and corona annuli (the same view as Photoshop's
"Variance" stack mode), the width of the disagreement rim just outside the
limb, and the 20-80% transition width of the merged limb.

### NAFE with a variable neighbourhood (v0.8.0)

Druckmuller 2013 (ApJ 775, 88) and Druckmuller & Druckmullerova,
"Noise Adaptive Fuzzy Equalization Method with Variable Neighborhood"
(IWCIA 2014, LNCS 8466, p. 262).

Each pixel is ranked within a fuzzy multiscale neighbourhood, and that
neighbourhood is restricted **in value** -- to neighbours of similar
brightness -- rather than by a geometric mask. This is the paper's fix for
"loss of contrast on boundaries between areas with significantly different
brightness", their Fig. 2 being captioned "loss of contrast near lunar edge".

It matters here for a second reason. MGN needs a `valid` mask and FNRGF needs
a radial origin, so both inherit any error in the limb fit. NAFE-VN needs no
geometry at all: the dark lunar plateau drops out of the corona's statistics
because it is dark, not because a circle was drawn around it. It is therefore
the one detail layer that stays correct when the limb fit is not.
