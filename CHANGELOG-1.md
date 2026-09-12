# Changelog — EclipseForgeHDR

Newest first. Entries from 0.6.1 onward were written at the time. The 0.7.2 –
0.9.8 block was reconstructed afterwards from the code comments that name their
own version and from the development record; where a change cannot be pinned to
an exact version it is filed under the release it is known to precede.

## 0.22.66

**Two bugs in how the Hill chain was wired into the composite.** Render-only:
a 0.22.65 work directory, masks included, is reused as it stands.

### The base stretch was fixed at Hill's k, which is wrong for our data

`im_log = log(1 + k*xn)/log(1 + k)` with k = 1e6 is calibrated for a
normalisation whose outer corona sits near 1e-6. Ours does not. Measured on the
reference bracket:

    xn 1.0e+00  (limb)      ->  im_log 1.000
    xn 2.0e-01  (1.5 R)     ->  im_log 0.884
    xn 5.0e-02  (2 R)       ->  im_log 0.783
    xn 1.0e-02  (3 R)       ->  im_log 0.667
    xn 3.0e-03  (4 R)       ->  im_log 0.580
    xn 5.0e-04  (far field) ->  im_log 0.450

The whole corona lands in the top third of the display range and the sky lifts
off the floor. That is the solid orange inner region and the grey background,
and no amount of amplification fixes it, because the masks are being added to a
base that has no room left.

`log10(k)` is now a slider, default 3, and it moves WITHOUT rebuilding
anything. The stored layer is invertible -- `xn = expm1(im_log*ln(1+Kb))/Kb`
recovers the normalised linear image exactly -- and the masks follow by a
scalar: wherever `K*xn >> 1`,

    log(1+K*xn)/log(1+K) = [ln K + ln xn] / ln(1+K)

so an unsharp mask of `im_log` is an unsharp mask of `ln xn` divided by
`ln(1+K)`, i.e. the same mask at any k up to that constant. Scaling by
`ln(1+Kb)/ln(1+Kr)` moves them exactly in the regime that carries the corona,
and leaves them bounded in the far field where they were built bounded. The
page does the same re-map through a 65536-entry table, which is exact because
the layer is transported at 16 bits.

Measured, fine-scale rms in the render: 0.018 at k=1e6, 0.026 at 1e4, 0.035 at
1e3, 0.067 at 1e2.

### The crossfade was comparing brightness, not structure

`Y` on the existing path is an envelope times a detail modulation; `Y_hill` is
a log-mapped 0..1 image. On the reference fixture their means differ by about
3x, so every intermediate `hillMix` was mostly just making the picture
brighter -- which is exactly what a blown-out inner corona looks like.

The two are now matched on their mean over the 1.15-3 R annulus before the
crossfade, so the slider compares structure at constant brightness. The page
accumulates both sums inside the pixel loop it is already running and applies
the ratio on the next frame, which it schedules itself; neither sum depends on
the ratio, so there is no feedback and it settles in one frame. A mean rather
than a median so the page does not have to sort a million floats per slider
move, and so the two implementations agree.

After the fix, mean brightness at hillMix 0 / 0.3 / 1.0 reads 0.147 / 0.142 /
0.130 against 0.147 / 0.30 / 0.44 before.

## 0.22.65

**The Hill masks can be built from a work directory that is already stacked.**
Press Start on a stacked folder and, if `hill.npy` is missing, it is built from
the merged luminance and the geometry that are already on disk -- the polar blur
alone, not a 13-minute re-stack. The same argument the RHEF layer makes.

Also: the renderer's fallback measurement of the mask scale (used when the run
report carries no `hill` section) measured over the WHOLE plane including the
occulted disc, where the masks are the difference between a black plateau and a
partial convolution with nothing to work with. On the test fixture that read
3.4e-2 against a true 5.7e-3 -- a factor of six, straight into what the
amplification slider means. It now excludes the disc, which is the same
estimator the builder records, and the two agree to under 1%.

## 0.22.64

**Hill's enhancement chain, built as he describes it**, plus a preview that no
longer aliases. A 0.22.63 work directory is reused as it stands; re-stack to
build the Hill masks (about 5 minutes on an 8000 x 5400 frame).

### The chain

    im_enhanced = im_log + a*M_1 + b*M_2 + c*M_4 + ...
    M_s (unsharp mask of size s) = im_log - im_blur(s)

straight off his slide. `im_blur` is the polar-oriented, PARTIAL blur of the
log-mapped image, with the Moon and the prominences excluded from the
convolution. Scales 2, 4, 8, 16, 32 px; amplification factors default to his
own 100 : 60 : 20 : 10 ratios, extended one octave, under one master gain.
Off by default (`hillMix` 0), and at 0 the render is bit-identical to 0.22.63.

**Why this is not a rearrangement of what was already here.** MGN, RHEF and
NAFE all NORMALISE -- they divide the local residual by a local statistic, or
replace it by its rank. That is what lets them show the inner and outer corona
at once, and it is also why every part of the frame ends up with the same
texture amplitude: faint structure and noise come out equally strong. Hill's
combination is LINEAR, five scalars, nothing divided by anything local, and it
is ADDED to the log image rather than multiplied into an envelope. Faint stays
faint. That difference, not the choice of kernel, is the gap between his slides
and our composite.

### The polar blur, and a correction

The note in detail.py said polar filtering had been "tested and rejected with
numbers". It had not. What was tested was DRUCKMULLEROVA's ACHF kernel, whose
arc-length term holds the tangential width at a constant physical distance --
which is exactly what makes it nearly isotropic at our scales, and that
measurement stands.

Hill's kernel is a different object. It is a plain Gaussian in the POLAR IMAGE,
so back in the picture the radial width is s px everywhere while the tangential
width is s * r / r_max: narrow near the limb, growing outwards. It averages
along a streamer and not across it, by a factor of r_max/r -- on this geometry
about 8x at the limb. Measured on an impulse at two radii, sigma_radial 8.00 px
at both, sigma_tangential 1.07 px at 0.12 r_max and 5.92 px at 0.74 r_max,
against 0.98 and 5.90 predicted. That is the anisotropy the earlier note said
was not there, and it is why Hill says the polar transform is the whole
difference between his result and a Photoshop high-pass.

Implementation notes: partial convolution as (w.f * C) / (w * C) with the same
kernel on the mask, per his slide; only the blurred component makes the polar
round trip, so interpolation can soften the blur but can never touch the detail
the mask carries; worked in radial bands with a 4-sigma overlap so the peak
footprint is one band rather than the whole polar grid; the angular axis wraps,
checked -- the step across the seam is smaller than the largest step elsewhere
on the same ring.

### The preview was aliasing

Preview layers were decimated with `v[::q, ::q]` -- plain subsampling, no
prefilter, 15 of every 16 pixels discarded. On a detail layer whose finest
scale is 1.4 px that is aliasing, and aliased fine threads do not read as
missing detail, they read as COARSE detail, because the energy folds down into
low frequencies as mottle. So the preview was a systematically rougher picture
than the export it stood for, and judgements about the render being too coarse
were being made on it.

Now area-averaged. Measured against the export downsampled properly: fine-scale
rms was 12% high and is now 3% high, and the rms error against that reference
halves. Exports never went through this path and do not change.

## 0.22.63

**White balance is now camera as shot, and VNG is available as a demosaic.**
Both are toolbar settings; both change the merged data, so both are in the
cache key and switching either one re-stacks. A work directory from any
earlier version has to be re-stacked as well.

### White balance: as shot, by default

The old default was LibRaw's `daylight_whitebalance`, chosen because a daylight
value is a property of the camera model rather than of whatever the
photographer had dialled in, and so is reproducible between people.

That premise does not hold for a body LibRaw carries no table entry for: the
value is DERIVED from the colour matrix, not measured. On the Panasonic
DC-S1RM2, normalised to G = 1:

    as shot (544/256/462)           R 2.125   B 1.805   R/B 1.177
    daylight (LibRaw pre_mul)       R 2.258   B 1.232   R/B 1.833
    the camera's own Fine Weather   R 2.227   B 1.707   R/B 1.305

The derived daylight value is 46% short of blue against as-shot and 39% short
against the camera's own daylight preset, so it rendered the corona 56% redder
in R/B than the camera saw it. That is a large part of why the app needed a
warmth control below 1 and a corona white balance at 0.66 to look right.

Daylight and None (1:1:1) remain available in the toolbar. Bodies that report
no usable as-shot value fall back to daylight; FITS and TIFF carry no white
balance and are unaffected. The run log prints which was used and what the
other one would have been.

### VNG interpolation

Implemented from the gradient tables in LibRaw's own `misc_demosaic.cpp` --
extracted from that source rather than retyped -- and following its reading of
them exactly: weighted bilinear first, eight directional gradients measured on
that, threshold at gmin + gmax/2, then the measured channel passed through
untouched with the other two rebuilt from the colour difference averaged over
the surviving directions.

Checked against LibRaw's own VNG output on a real frame of the reference
bracket. After applying LibRaw's 16-bit clip -- the only thing that ever
differed -- the largest disagreement over 44 million pixels is 1.5 counts in
65535, which is integer rounding against our float arithmetic. Banded over rows
so peak memory stays near the Malvar path's; the banding is exact (a
single-band and a 128-row-band render are bit-identical). Costs about 1.7x the
Malvar path's time on a full frame.

Malvar-He-Cutler stays the default: it has lower interpolation error, and it is
what every render up to 0.22.62 was made with. VNG is there so a render can be
compared with PixInsight's on equal terms.

### Also

- `feather` was checked in the server's cache key but never written into
  opts.json, so a non-plain merge weight compared against a default of "plain"
  and re-ran every time. It is written now.
- The contact (diamond ring) frame is decoded with the same white balance and
  demosaic the composite was built with, read from the run's own opts, rather
  than with today's defaults.

## 0.22.62

**Neutralise corona is back, exactly as 0.22.48 had it.** Render only — a
0.22.61 work directory is reused as it stands.

The tester's verdict, after everything: *"0.22.48 was still good."* His settings
in that build were corona 0.66, warmth 0.895, sky cast 1, and the picture is a
deep blue field with a white corona and a magenta prominence, with no edge
anywhere in it.

0.22.48's entire colour application was one line:

```python
if P.get("coronaNeutral", 0) > 0:
    a *= (layers.corona_gain[None, None, :] ** P["coronaNeutral"])
```

Flat. No mask, no fade, no per-pixel weight. Every version after it added a
weight to that line — fade by the corona fraction (0.22.50), by the chroma
confidence, by the sky fraction — and **every one of those weights drew a
boundary**, because a weight that varies with radius, applied to a gain of 4.7
in blue, IS a visible edge. Twelve releases of adding weights, and the verdict
on each was the same.

The thing they were all trying to correct is not an error. Out in the field the
chroma has already been faded to neutral, and a flat gain turns that neutral
into sky-blue — which is what makes the corona read as white against a sky
rather than as an orange smear on grey. 0.22.50 called that a green ring and
"fixed" it. It was the picture working.

Measured with the tester's settings, from 1.5 R to 7.5 R: R/G falls 1.18 → 0.58
and B/G rises 1.20 → 2.73, smoothly, with no step at any radius. That is why
there is no edge — the transition is monotone, so it reads as a corona against
a sky rather than as a boundary.

## 0.22.61

**The colour work is removed.** Everything from 0.22.49 to 0.22.60 that touched
colour is gone: the sky subtraction and its checkbox, `Neutralise corona`, the
corona white-balance gains, and all the chroma machinery built to hold them up.
Needs a re-stack, because the merge goes back to what 0.22.48 produced.

Ten releases of patching a thing that never once produced a picture the tester
would accept. The measurements underneath it are sound and the mechanism was
not, and continuing to repair it was the wrong call several versions ago.

**Kept**, because none of it is colour: the Hill photometric chain, the fitted
tier pedestal, the radial-flatten ring fix, and Hill's log stretch — which is a
tone control, not a colour one, and which the tester asked for by name.

**What survives as knowledge** is written up in `docs/SKY_SUBTRACTION.md`: the
two-component test that says the corona's colour is constant with radius (red
against green is a straight line to 1.4% over a hundredfold range in signal),
the 800-count additive sky that nothing in the pipeline removes, why every shape
model tried either invented colour or left a gradient, and the four traps found
on the way — the occulted disc owning the black point, clipping at zero
inventing colour, a grey pedestal diluting the corona, and the order of the fade
against the gains.

And the honest reason it took ten tries: every check was on luminance, and this
class of error is invisible there. The renders were never looked at until the
end. Any future attempt at this is judged on a rendered image first.

## 0.22.60

**A plane, not a quadratic** — and this one was decided by looking at the
render, not by a number. Supersedes 0.22.59, which was never usable: it fixed
the shared shape and kept the quadratic, and the quadratic is the other half of
the same problem. Needs a re-stack.

0.22.57 added the quadratic on the strength of `remove_sky_gradient`'s note that
one explains 51.8% of what a plane leaves behind. That measurement is about ITS
fit — multiplicative, in log space, one shared shape for brightness. This fit is
additive and per channel, and the extra freedom goes somewhere else: the
amplitudes the three channels want diverge, and the difference between them is a
colour pattern with the quadratic's own shape. Rendered on the reference bracket
with the quadratic: a blue region across the upper right, a warm one on the
left, plainly visible. With a plane, same data, same settings: grey and even.

The far field is where the corona is faintest but never absent, so every term
added here is another way to absorb corona and call it sky. A plane can tilt; it
cannot sculpt.

## 0.22.59

**The sky subtraction was painting colour onto the background.** Broad pink and
blue regions across the frame, over an image whose luminance measured flat.
Superseded by 0.22.60 — this half of the fix is necessary and was not
sufficient.

The surface was fitted independently in each channel. Free to do that, each one
absorbs a different part of whatever is left in the fitting annulus — and what
is left there is corona, which is very coloured (R/G about 2). The three fitted
surfaces came out with genuinely different SHAPES. Each normalised to its own x
coefficient:

```
        x       y      xx      yy      xy
  R   1.00   -0.52   -0.07    0.58    1.99
  G   1.00   -0.69   -1.16   -0.26    1.20
  B   1.00   -0.74   -1.21   -0.41    0.69
```

Red has no xx curvature where green and blue have a great deal. Subtract three
surfaces differing like that and you have subtracted a colour map: a pink region
on one side, a blue corner on the other.

The sky has one shape. Its brightness varies across the frame, and its colour
varies a little along that same gradient — more airmass lower down means more
reddening — which a per-channel AMPLITUDE on a shared shape expresses exactly.
Fifteen free parameters express that too, and also express things that are not
sky. Five plus three cannot invent a colour pattern.

The shape is now fitted once, on luminance, and each channel gets one amplitude.
On the reference bracket the straight-line residual of the corona model drops
from 2.5% to 0.3% at the same time — the shared shape describes the data better,
not worse.

**How this was found, which matters more than the fix.** Every check I had was
on numbers, and the numbers said the background was flat: far-field luminance
even to 15 counts across the frame after 859 were removed. It was flat. The
error was entirely in colour, and nothing measured colour at large scale. It
took rendering the image and looking at it. Renders are now part of checking a
change to the colour path, not just profiles.

## 0.22.58

**In 0.22.57 the sky subtraction never ran.** One line in the log said so:

```
sky subtraction skipped (output array is read-only)
```

`hdr` is a memory map opened read-only, and `np.asarray(hdr[::4, ::4],
np.float32)` does **not** copy when the dtype already matches — it returns a
read-only view. 0.22.57 changed the fit to take the sky's shape out of its
input in place, which on that view raises. The exception was caught one level
up and the run fell back to the old gradient division, so the sky level was
never removed and the picture looked as though the subtraction simply had not
helped.

The array is now copied explicitly, in both the caller and the fit. Verified by
reproducing the exact call — `load_big`, the same slice, the same dtype — and
asserting the view is read-only before running the fit through it.

The fallback message is no longer quiet either. It now says the subtraction
FAILED, that the sky level was not removed, and what ran instead.

**My own test was the reason this shipped.** I checked the new fit by handing it
`A.copy()`, which is exactly the case that cannot fail. Tests of a changed code
path go through the real caller from now on.

Verified end to end on the reference frames through the pipeline entry point:
far-field median by quadrant +54.8 / +49.2 / +63.3 / +50.0 on a 24-count
pedestal, disc at exactly 0, corona colour R/G 1.963 B/G 0.243 with a
straight-line residual of 1.8%.

Its own cache family: a 0.22.57 work directory contains the old gradient
division, not this.

## 0.22.57

**The sky subtraction was subtracting a flat level from a tilted sky.**
Needs a re-stack — this is the one that changes the merged data.

Your sky is not flat. On the 600 mm reference bracket it runs 153 counts in x
and 95 in y across the frame, on an 853-count sky — 18%, which is what a sky
does at 6.8° altitude across a 2.3° field.

`fit_airlight` measured that correctly and then threw it away. It carried a
guard reading "a sky tilt of more than 30% of its own level is not sky", the
two components summed past 30%, so the shape was zeroed and only the level was
removed. The result, measured on that run: one quadrant of the far field left
at −90 counts and another at +89, 17% of the frame below the black point and
rendering as solid black, the rest a coloured wash. With the corona
neutralisation up it showed as a blue quadrant and a red one.

Two things were wrong, and one more was missing:

- **The guard.** It now limits the fitted model's *span* against the sky level
  and allows a factor of two before refusing — the same shape of test
  `remove_sky_gradient` has always used — instead of a third.
- **The model.** A plane, now a quadratic where there is enough sky to support
  one. `remove_sky_gradient`'s own note has said for many versions that on what
  a plane leaves behind, another plane explains 0.8% and a full quadratic 51.8%.
  The sky near the horizon is not linear.
- **The black point** is now never above the 0.5th percentile of the luminance.
  If the sky model is right this changes nothing; if it is wrong, it is what
  stops a quadrant going solid black before anyone can see why.

Measured on the reference frames, far-field median by quadrant after the
subtraction:

```
              before        after
upper-left     +8.2         +47.1
upper-right   +89.0         +41.8
lower-left    −90.1         +57.3
lower-right   +75.8         +44.6
```

A spread of 180 counts becomes 15, on a pedestal of 23 — and what is left is
real corona, which at 5 R is about 25 counts of it. Pixels below the black
point: 17% before, none after (the 2.8% the test still shows is the occulted
disc, which 0.22.53 restores to zero).

**Every work directory stacked with Subtract sky by 0.22.52–0.22.56 has the
flat-sky error baked into it.** The cache key cannot distinguish them — the
option is the same — so this release opens its own cache family and everything
re-stacks.

## 0.22.56

**The blue band at 4.5–5 R is gone.** Render-only — no re-stack.

I had called it the measurement floor. It was not; it was two mistakes of mine
stacked on each other, and the data was never the limit.

**The grey pedestal was diluting the colour.** `subtract_airlight` deliberately
leaves a small neutral pedestal so the far field stays clear of the log floor
the detail layers use. Grey cannot print a coloured ring by itself — but it
dilutes, and the corona gains then amplify the dilution. Blue-to-green as the
renderer measured it, with a 23-count pedestal:

```
radius      2R      4R      5R      6R
measured   0.257   0.320   0.434   0.675
corona      0.24    0.24    0.24    0.24
```

At 6 R the corona contributes about 20 counts of green and 5 of blue, and the
pedestal adds 23 to each — dragging the ratio most of the way to neutral. Blue's
gain is 4.7, so that drift came out as a band. The pedestal is now subtracted
again before the chroma is formed, which is exact because the stack records it.
The luminance path keeps it, which is what it was for.

**The fade to neutral was in the wrong place.** `ratio` was faded towards
neutral where the signal is weak, and the gains were applied afterwards. Those
two do not commute: a blend towards neutral RAISES blue, which sits far below
neutral in a corona, and the gains then multiply the raised value. At 5 R, where
confidence is 0.75:

```
measured      B/L 0.21   G/L 0.87    B/G 0.24  (correct)
faded first   B/L 0.41   G/L 0.90    B/G 0.46
then gains    B/L 1.94   G/L 1.04    B/G 1.87  <- the ring
```

And it was a *ring* rather than a drift because confidence falls to zero further
out: the error rises, peaks where confidence is partial, and collapses again.
Two monotone curves multiplied, one rising and one falling.

The fade now happens once, at the end, in log chroma — `(rc·g)^conf`. At
confidence 1 that is the corrected colour, at 0 exactly neutral, and in between
the same correction applied less.

Result, blue-to-green with the corona neutralisation up, 1.1 R to 7 R:

```
before   1.51  1.63  1.79  1.74  1.59  1.44  1.32  1.22  1.13  1.04  1.01
after    1.06  1.03  1.03  1.03  1.02  1.01  1.01  1.01  1.01  1.01  1.03
```

Flat to 1% across the whole field.

What is left is red: a broad rise to R/G 1.16 at 4 R, easing back to 0.96 by
7 R. It is not a level error — nudging the fitted levels to flatten it makes
blue worse by more than it makes red better, in both directions — and the
F-corona genuinely reddens outward, so some of it is the picture rather than the
pipeline. It is gradual over five solar radii rather than a band.

## 0.22.55

**Fixes the orange haze over the whole frame** that 0.22.52 put on a
sky-subtracted stack. Render-only — no re-stack.

0.22.52 made the chroma fall back to the corona's own measured colour, rather
than to neutral, where it could not be measured. That is true of the data and
wrong for the picture: the fallback is decided when the layers load, while the
thing that would neutralise it — the corona gains — is a slider set afterwards.
With **Neutralise corona** at its default of 0, the entire outer field was
painted with the corona's uncorrected colour, which at 6.8° of sun altitude is
R/G 1.98. A deep orange haze over everything.

A default must not depend on a slider agreeing with it. The fallback is neutral
again, always.

The problem it was reaching for is real and is now solved where the slider is
known: the corona gains are weighted by the chroma confidence as well as by the
corona's share of the pixel. These gains are a real transform — blue ×4.7 — so
applied to a pixel holding only noise they turn grey into blue, and the far
field is not black in the render because radial flattening lifts it. Measured
beyond 6 R: B/G 4.2 without the weight, 1.07 with it.

What is left, honestly: a blue band at 4.5–5 R reaching B/G 1.9. There the
corona is about 6% of the sky that was removed, so the 1% the level is good to
is a 30% error in blue. It is the same measurement floor described in 0.22.52,
now visible in one ring rather than spread over the bright corona. The log
stretch will make it more visible, not less.

## 0.22.54

**Hill's log stretch**, as an alternative to the cube-root envelope. New
slider **Log stretch (Hill)** under STRUCTURE, 0 = off and the picture is
bit-identical to before.

`log(1 + k·I)`, with the slider setting log10(k); Hill's own value is 6. He
uses it because it "closely resembles how the corona looks to our eyes", and
the reason it does is that a logarithm gives every factor of two in brightness
the same amount of display range, wherever it sits. Measured on the
sky-subtracted reference bracket, output units per doubling:

```
radius       1.1R    2R     3R     4R     6R
cube root   0.197  0.057  0.044  0.041  0.039
log, k=1e6  0.050  0.050  0.050  0.050  0.050
```

This became worth having only now. With the sky still in the picture the
corona's own range is about 130:1 and any sensible curve shows all of it; with
the sky subtracted it is 1900:1, and the cube root spends most of its output on
the inner corona. Rendered luminance with radial flatten off:

```
             1.1R    1.5R    2R      3R      4R      6R
cube root   0.546   0.254   0.143   0.087   0.069   0.043
log, k=1e5  0.671   0.537   0.435   0.345   0.302   0.209
```

**Turn radial flatten down when you turn this up.** They do the same job from
opposite ends, and both at once inverts the picture — at k=1e6 with the
default flatten of 0.5 the field at 4 R renders brighter than the corona at
1.1 R. The pairing that works on the reference set is log stretch 5–6 with
radial flatten at 0.

The `bg` layer now goes to the browser at 16 bits, packed high/low byte like
the sky layer. The envelope has to cube it to recover luminance, which
multiplies the quantisation error by three and does it exactly where the far
field lives: at 8 bits the outer corona sits near byte 16, where one step is
19% in normalised luminance — visible banding in the region the log stretch
exists to show. Measured transport error at k=1e6: 8.6e-5 against 2.2e-2.

## 0.22.53

**The sky subtraction was being handed straight back by the display mapping.**
Found while reading Hill's stretching section, not from a symptom.

The renderer takes its black point from the luminance, by two rules: a low
percentile, or the sky median minus five sigma. With the sky subtracted both
are wrong, and wrong by about the amount that was removed.

The merge sets the occulted disc to exactly zero — it is masked, not measured
— so subtracting an 800-count sky drove it to −800, and at 3% of the frame it
owned the low percentile. The sigma rule failed differently: it takes the MAD
of everything beyond 2.5 R, which after the subtraction is the corona's own
falloff rather than sky plus noise, so 81 counts of structure were read as 81
counts of noise and the rule returned −309 on a background sitting at 24.

Either way the black point put the pedestal back. Measured, `xn` came out
**identical at every radius** with and without the subtraction — the entire
gain in outer contrast, cancelled.

The stack already records what it left behind: the grey pedestal and the sky's
own pixel noise, the latter from neighbouring-pixel differences so the corona's
gradient cannot inflate it. The black point now comes from those two numbers,
three sigma below the pedestal. Same shape of rule as PixInsight's STF
autostretch, which clips 1.25 sigma below the median; ours is looser because
nothing here may crush sky to black.

What it buys, on the reference bracket — normalised luminance by radius:

```
            1.1R      2R        3R        4R        6R
before   4.33e-1   1.03e-2   4.68e-3   3.78e-3   3.30e-3
after    4.31e-1   7.21e-3   1.62e-3   7.11e-4   2.32e-4
```

The corona's own range goes from 130:1 to 1900:1. Before, everything past 3 R
was sitting on a pedestal fourteen times larger than the corona there.

The subtraction also no longer touches the disc: zero stays zero. That part
changes the merged data, so it applies to the next stack — a 0.22.52 work
directory renders correctly as it stands, because the new black point does not
look at percentiles at all.

**Consequence worth knowing:** with a real 1900:1 range the cube-root envelope
now leaves the outer corona much darker than it used to. That is the honest
picture rather than a fault, but it is what Hill's log stretch exists to solve
— see the note in the code.

## 0.22.52

The sky is now **subtracted from the image** rather than corrected in the
colour. New checkbox, **Subtract sky**, off by default; switching it on
costs a re-stack because it changes the merged data and every layer built
on it.

0.22.51 fitted the sky correctly and then did the wrong thing with the
answer: it corrected the sky's colour at render time and left the sky's
*brightness* in the picture. That cannot work, and the arithmetic says why.
The correction divides by the corona's share of each pixel, so the 1% the
fitted level is good to comes back divided by the same number — 27% in blue
where the corona is a seventh of the sky. Painted across a large, still-bright,
smooth area, where the eye reads 2% as a cast. Measured, it was a 13% blue
deficit over half the frame; reported as an olive disc. Tuning the weights
could not have fixed it.

Subtract the sky properly and that error lands on a part of the picture that
is now nearly black. This is why PixInsight's DBE is subtractive and needs no
masks, and this follows it: a smooth model per channel, fitted beyond the
corona, taken out of the merged image before any layer is built.

Measured on the reference bracket, R/G and B/G with the corona neutralisation
full up:

```
                1.1R   2.0R   2.5R   3.0R   4.0R
before  R/G     1.17   1.20   1.19   1.16   1.10
        B/G     1.51   1.79   1.74   1.59   1.32
after   R/G     0.99   1.01   1.03   1.05   1.04
        B/G     1.06   1.06   1.11   1.15   1.33
```

The ring over the bright inner corona is gone, and with the sky out of the way
the corona's own colour is the same number at every radius — so the corona
gains apply flat, with no radial fade.

Three things fall out of it, all measured rather than assumed:

- **Nothing is clipped at zero.** Rectifying the noise looks harmless and
  invents colour: each channel's rectified mean is proportional to its own
  sigma. Clipped, the far field came out R/G 0.55 and B/G 0.97 where the
  corona has neither, and the gains turned that into a blue outer field.
- **A grey pedestal of two sky-noise sigma is left behind** — 27 counts against
  the 806 removed. Subtracting to exactly zero puts much of the far field under
  the `log10(clip(lum, 1))` floor every detail layer uses.
- **The chroma falls back to the corona's own measured colour**, not to neutral,
  where it cannot be measured. Neutral is right while the sky is in the picture;
  with the sky gone it hands a grey pixel to gains that are a real transform
  (blue ×4.7 here) and turns it blue.

**Neutralise sky cast** becomes inert on such a stack — there is no sky cast
left, and its default of fully-on was multiplying blue by 1.47 across the whole
picture.

What remains: a blue drift beyond about 4 R, reaching B/G 1.5 by 5 R. That is
the 1% limit on measuring an 845-count sky, in a region where the corona is a
fifteenth of it. Two defensible estimators for the level disagree by exactly
that much and move it in opposite directions, which is what says it is the
measurement limit and not a bug. It sits where the corona is faint rather than
over the bright inner corona, which is the trade this change makes.

The render-time **Remove airlight** slider from 0.22.51 is gone. It was the
wrong mechanism and leaving it as an option would only invite tuning it.

## 0.22.51

The sky's own light is now separated from the corona instead of being
tinted along with it. New control, **Remove airlight**, off by default.

The atmosphere does two different things to a low-sun eclipse. It
*removes* corona light on the way down — multiplicative, uniform over a
field this small, strongly wavelength-dependent, and the reason a corona
photographed with the sun 6.8° up comes out orange. And it *adds*
scattered sunlight along the line of sight, which has its own colour and
its own gradient. The first is what "Neutralise corona" corrects. Nothing
in the app had ever removed the second: the sky-gradient step divides,
which is the wrong shape for something additive, and it takes out only
the tilt — the level has always been left in the picture.

That level is why the corona neutralisation printed a ring. On the
reference bracket the airlight is 845 counts in green, against a corona
that has fallen to 170 counts by 3.5 R. Out there most of the pixel is
sky, so gains measured on the inner corona over-correct by an amount that
grows with radius. 0.22.50 faded them by the corona fraction, which
turned an over-correction into a correction that varies with radius — a
ring either way, as the tester reported.

The separation is measured, not assumed. Writing each channel as one
shared corona shape plus a per-channel sky level predicts that red
plotted against green must be a straight line, whatever the corona looks
like. On the reference bracket, radial medians from 1.2 to 6.9 R — a
hundredfold range in signal:

    R = 1.926 * G - 844      median residual 1.4%
    B = 0.244 * G + 344      median residual 0.8%

So the corona's colour does not change with radius at all; the entire
apparent outward shift (R/G 2.05 → 0.91, B/G 0.31 → 0.64) is the airlight
taking over. The fitted level, 781 / 844 / 550, agrees within 3% with the
far-field median measured independently at r > 5.5 R, 785 / 866 / 558 —
two routes sharing no arithmetic.

What it changes in the rendered picture, measured as R/G and B/G from
1.1 to 7 R with the corona neutralisation full up:

    before   R/G 0.99–1.16   B/G 1.04–1.95, peaking at 2 R
    after    R/G 0.99–1.06   B/G 0.86–1.12, no peak

The peak is the ring; it is gone, and the blue swing is a quarter of what
it was. With the airlight out, the corona's colour is the same number
everywhere, so the extinction gains apply flat and need no radial fade.

Where the corona is a small fraction of the sky the subtraction divides
by that fraction, and the 1% the fitted level is good to comes back
divided by it as well — 27% in blue by 4 R. So the correction fades to
neutral below a corona fraction of 0.45, in log chroma rather than
linearly, because the gains are multiplicative and large. Far out, where
the corona's colour cannot be measured at all, the answer is grey rather
than something invented.

**Neutralise sky cast** now fades out as this comes up. It divides by the
sky's colour; this removes the sky. Doing both was correcting the same
thing twice, and it was over-correcting the far field — B/G 2.56 at 5 R
with the sky already subtracted.

The fit refuses datasets it cannot describe — too little field beyond the
corona, a level that disagrees with the far field, or channels that are
not linear in each other — and the control is then absent rather than
wrong. It runs at load time on the merged HDR already on disk, so no
re-stack: an existing work directory gains the control.

## 0.22.50

The corona neutralisation no longer tints the outer field green.

The correction is measured on the inner corona, which is white by
physics, but it was applied at full strength across the whole frame. The
colour cast it corrects is not uniform: it fades outward, from roughly
twice as much red as green at the limb to nearly balanced at four solar
radii. Applied out there it therefore overshot, leaving green as the
strongest channel — a green disc around the sun, reported by a tester.

Turning the sky-cast control up hid it, which is two wrong things
cancelling rather than either being right.

The correction now fades with how much of each pixel is corona rather
than sky, using the sky level the renderer already measures. At the limb
it is unchanged; in the far field, where the sky dominates and is
genuinely blue, it barely applies at all. The sky keeps its own colour,
which it should — nothing about it is white by physics.

## 0.22.49

The corona neutralisation no longer discolours prominences.

Neutralising the corona means multiplying blue by about four, because
the corona arrives with barely a quarter the blue it should have. A
prominence is hydrogen emitting at a single red wavelength and has
almost no blue to begin with, so quadrupling it turned crimson into
magenta.

The correction now fades out inside the prominence gate the app already
computes, so prominences keep their own colour while everything else is
corrected. The fade is gradual rather than a switch, so there is no
colour seam at a prominence's edge.

Measured on the reference colours: a prominence went from red 3.30,
green 0.39, blue 0.33 to red 2.35, green 0.54, blue 1.64 — blue ending
up well above green, which is what magenta is. It now keeps its
original colour.

## 0.22.48

**A control that neutralises the corona's colour.** New "Neutralise
corona" slider, in the Color group, off by default so existing pictures
are unchanged.

The inner corona is photospheric light scattered off free electrons, and
that scattering does not depend on wavelength — so the inner corona
carries the sun's own spectrum and is the one thing in the frame that is
white by physics rather than by convention. The app has always known
this and has always used it, but only for files that arrived with no
camera white balance at all, on the assumption that a raw bracket keeps
the camera's own numbers and is therefore already right.

Measured on a 600 mm raw bracket that does have camera white balance,
the merged corona reads twice as much red as green and barely a quarter
the blue. The app's own prominence detection reports the same thing from
the other side and has been calibrated around it for as long as it has
existed. So the restriction was wrong, and the correction now belongs to
every dataset as an adjustable control.

It is measured when the picture is loaded, so it needs no re-processing
and costs nothing until the slider is moved. At 1 the inner corona is
exactly neutral; in between it is applied partially. Brightness never
changes — only colour.

**Radial flatten no longer prints rings.** Two causes, both measured on
real data after a tester reported that the control introduced a radial
pattern which vanished when it was turned off.

The radial brightness profile was being read at whole-pixel radius, so
the divisor was constant inside each one-pixel ring and jumped between
them — dividing the picture by a staircase. It is now interpolated,
which reduces the steps by a factor of two hundred.

Worse, the on-screen preview received that profile as an 8-bit image.
Across the outer field only a few dozen radii differed at all, leaving
wide flat rings separated by jumps of up to 3%. The page now receives
the profile as numbers — one per radius, a few thousand of them — and
evaluates it directly. An older cached run without them falls back to
the previous behaviour.

## 0.22.47

The optional Hill mode now applies the shared offset as well as its own
per-pair ones. Default runs are unaffected.

When the mode was added, the app's single shared offset was switched off
on the assumption that the per-pair offsets replaced it. They don't. A
chain of pairwise fits matches each exposure to its neighbour, so it can
correct what DIFFERS between two exposures but is blind to anything the
whole bracket has in common — there is nothing left to compare that
against.

The leftover then reached the sky-removal step and was taken out as if
it were sky: on a test bracket that step's blue correction rose from
1.41x to 1.83x, doing work that wasn't its own.

The shared offset is now removed before the pairwise fits, so the two
compose instead of one standing in for the other. Verified by injecting
an offset common to every exposure plus a different one per exposure:
the shared correction takes the common part and the fitted offsets
recover the per-exposure part to a hundredth of an ADU, with the
brightness factors unchanged.

## 0.22.46

Two corrections to the optional Hill photometric mode, both from its
first run on a real bracket. Default runs remain unaffected.

That run made the exposures agree markedly *worse* with each other —
tier-to-tier variance at the limb rose from 0.074 to 0.206 — while
looking identical. Two causes, both about bracket depth rather than
about the method.

**A link with nothing in it was still given two free parameters.** On a
fourteen-exposure bracket the shortest end has almost no overlap: the
two fastest exposures had zero and twenty-three pixels respectively
carrying real signal in both frames of the pair. Fitted anyway, the
slope and the offset traded against each other and returned a
correction of 2.53x where the previous method read 1.27x, and that
error then travelled down the whole chain. Each pair is now given only
as many free parameters as its overlap can support: slope and offset
where there is plenty, slope alone where there is little, and the plain
exposure ratio where there is none.

**The chain now starts from the middle exposure, not the longest.**
The published method chains from the longest and accepts that errors
accumulate along the way — reasonable for the four or five frames it
was described for. With fourteen, the shortest exposure sits thirteen
links from the longest. Starting in the middle halves the worst case
and changes nothing else; it is also where the app's original method
has always started.

The log now states which model each pair received and how many of each.

## 0.22.45

Fixes how an imported image is scaled. Affects the "import a finished
image" path only; processing raw frames is untouched.

The scale was decided from the brightest pixel in the file: above 255 it
was treated as 16-bit, otherwise as 8-bit. That misreads two real cases.
A 16-bit image whose brightest pixel happens to be below 256 — a faint
corona exported at 16 bits easily is — was treated as 8-bit and arrived
257 times too bright. A floating-point image carrying values above 1,
which some software writes for linear data, was divided by 255 and
arrived almost black.

Neither showed up as a wrong overall brightness, because the import
renormalises at the end. They showed up in the *shape* of the image: an
sRGB curve has to be undone before the data is linear, and that
undoing is only correct if the values are in the range it expects.
Measured on a dark 16-bit test image with a known answer, the recovered
image was off by 91%; it is now off by 11%, the rest being the
coarseness of the file itself.

The scale now comes from the file's declared type — full scale for its
bit depth for integers, and for floating point the established 0-to-1
convention, falling back to the image's own maximum if it exceeds that.

Files written by PixInsight were never affected: its TIFF and FITS
readers normalise floating-point data to 0–1 on the way out, so what it
writes always took the correct branch. If you have previously imported a
floating-point image, or a very dark 16-bit one, clear that folder's
cache once so it is re-read.

## 0.22.44

Corrects how the new photometric mode fits its line, after reading
PixInsight's own documentation rather than assuming.

`LinearFit` is **not** least squares. Its developer reference describes it
as "robust straight line fitting by minimization of mean absolute
deviation", calls the intercept "a constant additive pedestal present in
the whole dataset", and lists no rejection settings at all — the
resistance to outliers is that objective, not a rejection step. Separately,
the tool wrapped around it discards pixels outside a brightness range
before fitting.

0.22.43 used least squares with sigma clipping, on the reasoning that
"fit, discard outliers, refit" amounts to the same thing. It does not.
The two treat bright oddities differently, and on a corona those
oddities are stars, prominences and alignment residue along steep edges
— real structure rather than noise. Both mechanisms are now implemented
as documented.

The brightness limits are 0 and 0.92 of full scale, confirmed against a
real PixInsight installation rather than taken from documentation —
neither documentation page states them. For comparison, the merge weight
here already falls to zero above 0.97 of saturation and begins its
rolloff at 0.87, so the published limit sits between the two.

Verified again on the synthetic bracket: injected per-channel gains are
recovered exactly, and the fit reproduces the slope and intercept from
the worked example in the source material (3.0500 and −0.010) on data
built to match it.

## 0.22.43

Adds Hill's photometric combine as a selectable alternative to the one
the app has always used. **Off by default — a normal run is unchanged.**

The app corrects each exposure with a single brightness factor, plus one
additive offset shared by the whole bracket. The established method, as
used in PixInsight and described in Hill's HDR eclipse workflow, instead
fits a straight line between each pair of neighbouring exposures — a
slope *and* an offset, for each colour channel separately, with outlier
rejection — and chains those from the longest exposure. Hill states that
omitting this step is what produces ringing in the corona, which is an
artifact reported on several datasets here.

Measurements taken on one 14-tier bracket suggested the gain would be
small: the colour freedom bought nothing above the roughly 1% the
measurement itself could resolve, and the per-pair offset helped most on
the long exposures, which build the outer corona rather than the
near-limb band where the rings appear. That is an argument about size,
not about correctness, on one dataset, with an estimator that has been
wrong before. The published method is the defensible one to offer, so it
is offered, and the measurement is recorded beside it.

To switch it on, put a file named `hill_photometry.txt` in the folder of
raw frames, then clear the cache and run. (`ECLIPSEFORGE_PHOTOMETRY=hill`
does the same from a command line.) The log states plainly which mode
produced the image and what the two chains differ by.

One property of the method worth knowing: a chain of pairwise fits
equalises every exposure onto the reference and therefore cannot see an
offset the reference itself carries. A black level common to the whole
bracket survives this mode, where the shared offset removed it. The
per-exposure part — the part that arrives divided by the shutter speed
and shows up as disagreement in the outer field — is removed.

Verified against a synthetic bracket with a known per-channel gain drift
and known per-exposure offsets: the gains are recovered exactly and the
offsets on every exposure but the reference.

## 0.22.42

Two faults in the per-channel measurement, both found on a real run.

A link whose signal band could not be fitted was recorded as a ratio of
exactly 1.000 — "the three colour channels agree perfectly" — which is
indistinguishable from a clean result. Two links of a 14-tier run
reported exactly that. Such a link is now left out of the comparison
instead of counted as agreement.

The run also stores the pixel pairs each fit was made from, in
`.eclipseforgehdr/per_channel_samples.npz` (a few MB, thinned). Every
question asked of this measurement so far — chained or per link, which
signal band, which estimator, how to reject outliers — has required
another full re-stack to answer. The samples are exactly what the
estimator sees, so the next estimator can be tried against the same real
data without re-processing anything.

Still diagnostic only.

## 0.22.41

The per-channel check now repeats each fit in two signal bands — faint
pixels (2–20% of saturation) and bright ones (35–85%) — and reports both.

This tells apart the two things that produce the same headline number. A
true per-channel gain difference does not know how bright a pixel is, so
it reads the same in both bands. Sensor non-linearity near the top of the
range does depend on brightness, and because the three channels sit at
different raw levels for the same scene it hits them unequally — which
manufactures a colour difference that is not one. On test brackets: an
injected 3% colour gain reads 3.0% in both bands; an injected roll-off
reads 0.2% faint against 5.8% bright.

The bands are defined once from the longer exposure's brightness, so all
three channels are compared on identical pixels, and only links that both
bands could fit are compared.

Still diagnostic only.

## 0.22.40

The per-channel photometry check added in 0.22.39 failed on brackets whose
stacked mosaic has an odd number of rows or columns, logging

```
per-channel photometry check skipped (operands could not be broadcast
together with shapes (2008,3010) (2007,3010) )
```

and producing no measurement. The alignment trim can remove an odd number
of border rows, after which the two green sub-grids of the mosaic differ by
one row and cannot be averaged. The colour planes are now cropped to even
dimensions, and the planes, the masks and the ring grids are all held to a
single common shape.

Diagnostic only, as before: nothing else in the run is affected, and a run
that hit this still produced a correct image.

## 0.22.39

Every photometric correction in the merge is fitted on brightness alone
and applied to red, green and blue alike: one scale factor per exposure
tier, one shared additive offset, and one azimuthal gain profile. If the
variation between exposures carries any colour — thin cloud, changing
air mass, the sky brightening through totality — that colour cannot be
expressed, and it leaves the merge as a cast that changes with radius.

Established methods for this step fit a slope **and** an offset **per
channel** for each pair of exposures. Published descriptions also name
ringing in the corona as the direct consequence of omitting the step,
which is the artifact reported on several datasets.

This release measures whether it matters, and changes no output. Each
run now refits all three corrections per colour channel and records the
spread in the report and the run log: the per-tier slope, the outer-field
offset, and the azimuthal gain. Slope and offset are fitted **together**,
because fitting them one after the other does not work — on a test
bracket carrying a pure per-channel offset, measuring the slope first
reports a 39% colour difference that is not there, and measuring the
offset first on a bracket carrying a pure colour gain invents ±17 ADU of
offset. A single joint fit recovers both.

Nothing is corrected. The numbers decide whether a per-channel chain is
worth building, and the check can be turned off with
`ECLIPSEFORGE_NO_CHANNEL_CHECK=1`.

The merged image is unchanged, so an existing work directory is still
reused and no re-stack is forced. The measurement is taken during the
merge, though, so a reused run does not produce it: **to get these
numbers on a dataset that has already been processed, run it once with
force.** The result is the same picture.

Also corrected: the 0.22.37 note stating the Panasonic S1R II black
level as 576 with a 16383 ceiling. The correct figures are 512 and
16319. Both give the same usable range, so no run could distinguish
them, but only 512 is consistent with the frames themselves and with
the shared pedestal that real runs fit. No code behaviour changes — the
value has always been read from the file.

## 0.22.38

Exporting a 16-bit TIFF from a 600 mm Panasonic S1R II dataset left the
app showing "exporting tif16 ..." forever, though the 207 MB export had
finished on disk.

The polling loops in `gui.html` rescheduled themselves at the end of an
unprotected async body; any thrown error ended the chain uncaught,
freezing the button and status bar.

`pollAgain()` now retries with backoff to 4 s, resets on success, and
reports after ~6 s of silence instead of appearing frozen.

Also fixed: the orientation listener was registered inside the export
click handler, doing nothing until the first export, then duplicating
each export after.

## 0.22.37

A 600 mm Panasonic S1R II run on 0.22.36 logged:

```
saturation level 15412 ADU above black, from the reported white level
(no linearity limit given)
```

Correct, a no-op here: that body's black level is 512 and its ceiling is
16319, giving (16319 − 512) × 0.975 = 15412.

`color_info` comes from the first (shortest) frame, so the cross-check
never fires here: raw maximum 1766 ADU, nowhere near saturation.

(Corrected in 0.22.39. This entry first read "black is 576 (64 + 512
cblack)" with a ceiling of 16383, which gives the same 15412 because
16383 − 576 and 16319 − 512 are both 15807 — so no report can tell the
two apart. Measurement on the frames can, and says 512: decoded without
modification the sky reads a median of exactly 512.0 at 1/4000 s and
rises linearly with exposure above it, and only a black of 512
reproduces the shared pedestal that real runs on that dataset fit.)

`RawFile` now records and logs `sat_source` (`linear_max`, `observed`,
`white_level`, or `container`). A 0.22.36 stack is reused as-is.

## 0.22.36

On a 250 mm Canon EOS Rebel T7 dataset with Photometric, the core of a
prominence rendered as a black sliver; with Detail it did not.

`hdr = acc / max(wsum, 1e-9)` evaluates to 0/0 at a pixel no tier can
hold — a hole no stretching repairs. Present since 0.22.35. Measured on
the shortest tier (1/125 s):

```
sat_level 13978 (<= 0.22.34)   threshold 13558 ADU        0 photosites clipped
sat_level 12925 (0.22.35)      threshold 12537 ADU      227 photosites clipped
```

This body's data plateaus at 13255 ADU above black, so the old
container-ceiling level could never trigger clipping. The 227
photosites form one sliver, 45 x 17 px.

Only Photometric and `masked` show it, since both zero a clipped
pixel's weight; Detail's plain Gaussian blur leaks weight in from
unclipped neighbours instead — the same leak behind the pink rim
elsewhere.

**Fix.** Where no tier has weight, the merge falls back to the shortest
tier's value (flagged as a lower bound); the lunar disc is excluded.
The report gains a `bracket reach` line with the pixel count — the real
fix is one more stop at the short end of the bracket.

## 0.22.35

0.22.34's saturation fix never executed: it read
`camera_white_level_per_channel` after the `rawpy.imread` block had
closed the file, the attribute raised, the exception was swallowed, and
every frame fell back to the container ceiling. A 250 mm Canon EOS
Rebel T7 run logged saturation level 13978 ADU above black on a body whose
linearity limit rawpy reports correctly when read inside the block.

Now read next to `white_level` while the file is open. 0.22.34's
products are identical to 0.22.33's, so neither is reusable.

## 0.22.34

Identifies the ring artifact's source: a 250 mm Canon EOS Rebel T7 dataset
(the same body is sold as the EOS 2000D and EOS 1500D; LibRaw normalises all
three to "EOS 1500D") run through 0.22.33 from CR2s and Lightroom DNGs of the
same nine frames:

```
                                     CR2        DNG
tier-to-tier variance at the limb   0.793      0.013
disagreement rim                    74 px      none detectable
worst per-tier radial departure      230%        18%
azimuthal gain spread                 81%        27%
alignment residual                 0.66 px    0.37 px
coronal range                      5.8 EV     7.5 EV
```

**Cause.** LibRaw reports two ceilings: `white_level` (container) and
`camera_white_level_per_channel` (linear_max, the true clip point):

```
CR2 (raw)      white_level 16383   linear_max 15092   black 2047
DNG (Adobe)    white_level 15092   linear_max 15092   black 2048
```

`RawFile` used `white_level`: `sat_level = (16383 - 2047) * 0.975 =
13977` against a real clip point of **13045** — 7% high, so every
blown-highlight photosite was declared valid and merged at full weight,
producing the ring at each tier's clip radius, matching within 0.6 px.
Adobe writes the true value into a DNG.

**Fix.** The linearity limit is now preferred when reported (guarded
positive, below the container ceiling, above half of it, and above
black); the log states the saturation level and ceiling used.

Closes most of TODO 1c and 1e: rings were not caused by the merge
weight, feather, flat, demosaic, filters, tier set, frame count, or
noise floor.

## 0.22.33

A 250 mm Canon EOS Rebel T7 dataset run on 0.22.32 produced broad radial
wedges in FNRGF, RHEF, Pellett, and faintly MGN — 0.22.26's code, first
visible here.

**Bug.** LDIC's `k_i(phi)`/`q_i(phi)` are fitted in 60 azimuthal
segments and smoothed by an order-4 trigonometric polynomial, so gain
should be smooth by construction. The merge instead did
`rgb *= _kk[_segmap]` — the raw segment index, a hard discontinuity
every 6 degrees, scaling with fit amplitude: a synthetic 30% gain gives
a 7% jump between segments; here, where k spans 82%, about 17%.

**Fix.** The merge now evaluates the polynomial via a 16384-bin table,
reducing the residual step from 7-17% to 0.026% of k. Verified: the
table matches the fitted samples at segment centres to 2.4e-4, keeps the
azimuthal mean at exactly 1.000 for k and 0 for q, and is exact identity
on identical tiers.

**Also fixed: extrapolation.** Fitted over 1.0-2.5 R but applied at full
strength to the frame corners (4.3 R here) — the other cause of the
outer-field fan. It now fades to identity over the half lunar radius
above the fit's outer edge.

Both present since 0.22.26, invisible when tiers agree closely.

## 0.22.32

**Clear cache button**, next to Start. Deletes the app's own
`.eclipseforgehdr` subfolder only — raws, flats, exports stay — so the
next Start redoes the pipeline from raw files. Requires two clicks, four
seconds apart, since it discards tens of gigabytes of computed products.

Previously nothing in the GUI could trigger this — the log referenced a
`force` flag that existed only on the command line, and deleting
`opts.json` manually was undone by the next run. The log now names the
button.

**The merge weight is now part of the cache key.** `opts.json` recorded
it since 0.22.31; the import path compared it but the stacking path
didn't, so switching Detail/Photometric and pressing Start reused the
old merge. 0.22.31 cache remains valid; reuse is now conditional on
weight matching.

## 0.22.31

Three releases (0.22.28-0.22.30) attempted to choose the feather
automatically without a reliable measurement — on a 360 mm Canon EOS
80D dataset the trial still read **99%** against a bench figure of
**13%** (**0.33** run on the dataset's own exported tiers).

**The trial no longer decides the outcome.** The merge weight is now a
toolbar control:

- **Detail (hides rings)** — the plain feather; default.
- **Photometric (no rim)** — the leak-free weight, rings visible.

The trial still runs and reports its measurement as advice. `opts.json`
records the setting and the cache check compares it, so switching
triggers a re-run. `ECLIPSEFORGE_FEATHER` still overrides.

**Rationale.** No default serves every bracket; a labeled switch beats
an estimator wrong by 7x.

## 0.22.30

0.22.29 fixed the leak that had stopped the feather trial running. On a
360 mm Canon EOS 80D dataset it then reported **96%**, chose plain, and
the pink rim stayed. A 600 mm S1R II dataset had reported 95% — two
brackets six times apart in the measured quantity, both reading ~95%:
the estimator wasn't measuring the leak at all.

**Bug.** `_pick_feather` treated `stacks_half` as a per-channel raw
maximum compared against `0.87 * sat_level`; it's actually post-flat
Bayer-quad luma, well below `sat_level` even at ceiling, so the knee
and hard cut almost never fired.

`sat_half` is correct: true where any photosite in the quad hit the
ceiling, measured before the flat. The trial now uses it for the hard
cut and validity, and scales the soft knee by the luma at which each
tier actually clips.

**Validated against the offline bench**:

```
                estimator   bench   decision
600mm S1R II      0.809     0.747   plain
360mm EOS 80D     0.329     0.126   taper
```

Magnitudes still differ from bench figures, but the ordering is correct
and both now sit clear of the 0.60 threshold.

## 0.22.29

Processing a 600 mm S1R II dataset followed by a 360 mm Canon EOS 80D
dataset in the same session left the second with the pink rim the
0.22.28 trial exists to prevent:

```
switches honoured by this build (0.22.28): ECLIPSEFORGE_FEATHER=plain
feather held at 'plain' by ECLIPSEFORGE_FEATHER
```

The trial never ran: 0.22.28 wrote its choice back into `os.environ`,
and since the app is a single long-running process, every folder
processed afterward inherited the first dataset's `plain` choice at
startup. The write was also unnecessary, since `_fm` is passed to
`_feather_weight` explicitly; removed.

Anyone who processed multiple folders under 0.22.28 must re-run all but
the first.

## 0.22.28

On 0.22.26, one dataset was clean, another showed a pink rim again —
both from 0.22.25's global plain-feather default:

```
plain feather's level at 1.02 R, against the leak-free merge
  600mm S1R II      0.747      no visible rim
  360mm EOS 80D     0.126      a pink rim around the limb
```

The feather is no longer defaulted: `_pick_feather` rebuilds both
weights on the half-resolution stacks before the merge, takes the worst
near-limb ratio over 1.01-1.15 R, and picks plain if its deficit stays
under 40%, else taper. `ECLIPSEFORGE_FEATHER` still overrides.

The threshold is 0.60, calibrated on the two datasets above; anything in
(0.13, 0.75) separates them, so it is loosely pinned.

### 0.22.27's alignment change reverted

A real 600 mm S1R II run, 0.22.26 against 0.22.27:

```
                                   0.22.26   0.22.27
alignment network residual (half)   1.17 px   2.67 px
per-tier lunar limb spread          8 px      12 px
track scatter about the line        1/2 px    2/4 px
that step's run time                13 s      1m25s
```

Reverted to a flat 25 px. The 0.22.27 sweep measured the wrong thing:
already-registered, resampled, mean-stacked tiers sharing one origin
(real pairs use per-tier origins up to ~50 px apart), with no
signal-weight mask despite `prep_pair` multiplying by `wgt`.

Also fixed in the trial: `_feather_weight` read its mode from the
environment, comparing plain against itself and returning 1.000 on
fully blown brackets (mode is now explicit); and with every tier
clipped, the leak-free merge had no weight to average, giving a ratio
of 8e17 (now degenerate above 1.2).

## 0.22.27

0.22.1 left open whether to enable phase correlation. Now measured: no.
Two already-aligned tiers (true shift zero) get a known sub-pixel shift
via FFT phase ramp; each estimator must recover it with the same
parabolic refinement. Self-test accuracy: ~0.1 px.

RMS error over 4 adjacent pairs x 8 shifts on a 600 mm S1R II dataset:

```
ours, plain cross-correlation                     3.02 px
best of 36 regularised phase-correlation settings 5.00 px
the same band-pass, WITHOUT the whitening         2.91 px
```

`H` is worth ~0.1 px; whitening costs 2.1 px. Textbook phase correlation
(no `H`, no regularization) is worse still: 12.5 px against 5.9 px on
the full pair set. The tangential-blur high-pass (`naavis/eclipsetools`)
measures 5.40 px; not adopted.

### High-pass filter was too aggressive

`prep_pair` high-passed at a flat 25 px (0.08 R at the reference
geometry), removing nearly everything the corona looks like. Sweeping
the radius, RMS error in px:

```
sigma        0.05R  0.10R  0.20R  0.32R  0.50R  0.75R   none
600mm +/-3    3.01   2.64   2.40   2.35   2.42   2.54   2.18
600mm +/-12   2.99   2.61   2.36   2.31   2.38   2.48   2.15
360mm +/-3    1.01   0.90   0.77   0.72   0.71   0.71   0.71
```

Set to `0.5 R`, clamped to [25 px, S/6]: 18-22% lower error on both
datasets. Dropping it entirely scores slightly better on one dataset
but isn't adopted: correlation then follows overall brightness, which
differs under thin cloud. Only the cross-exposure filter changed;
`prep_pc` remains at 30 px.

## 0.22.26

Two changes, from the collar measurement and the reference thesis.

### Clipping test targeted the wrong data

`cmax` came from the demosaiced result: two of three channels per pixel
pass through 5x5 Malvar-He-Cutler kernels with negative lobes, so a
pixel beside a saturated photosite can read below threshold though
partly built from it. One saturated photosite in a dim field, ceiling
16000 ADU, background 300:

```
      old test (demosaiced max)          new test (mosaic max)
        0.3    0.3    0.3    0.3          0.3   16.0    0.3   0.3
        0.3    4.2    8.2    4.2         16.0   16.0   16.0  16.0    x1000 ADU
        0.3    8.2   16.0    8.2         16.0   16.0   16.0  16.0
```

This contaminates 13 pixels; the old test flags only 1.
`raw.cfa_clip_max` now evaluates the mosaic directly, over the four
kernels' non-zero taps — a correctness fix, not a ring fix: hard
exclusion of a 2-16 px collar moved the ring by only 1-2%.

### The merge gains a thesis term

Druckmullerova, doctoral thesis eq. 4.15 — the LDIC composition this
pipeline adapts:

```
g(r,phi) = SUM_i  w(f_i) * ( k_i(phi) * f_i(r,phi) + q_i(phi) )
```

`k` and `q` are a per-image affine transform varying with azimuth,
fitted in 60 angular segments and smoothed by a trigonometric
polynomial, replacing one scalar per tier and one shared pedestal.
Fitted on the S1R II dataset, gain varies 3-27% around the limb by
tier, offset runs to 9% of local signal; median gain sits below 1 for
eleven of twelve tiers.

Applied mean-preserving, addressing only variation a scalar can't carry.
Bench: ring power 0.87x, radial profile held within 4%. Synthetic test: a
10% cosine gain recovered to 0.1% in amplitude and phase; identical
tiers give exactly k = 1, q = 0. `ECLIPSEFORGE_NO_LDIC=1` disables it.

The ring survives both changes, every per-tier correction tried, and the
unfeathered merge.

## 0.22.25

**The default merge weight reverts to the plain feather (pre-0.22.16).**
The 600 mm S1R II dataset was clean before 0.22.16 and hasn't been
since; a deliberate trade, stated on every run.

**Bench results across fifteen variants**, scored by ring power on the
S1R II dataset at its own exposure exponent (0.55):

```
variant                              1.02R    ring power
none (no feather at all)             1.000       0.94
the shipped taper (.19-.24)          1.000       1.00
plain blur (<= 0.22.15)              0.747       0.66
hard collar exclusion, 2-16 px       0.998    1.01-1.02
inverse-variance weighting           0.992       1.15
C1 roll-off to zero below the cut    1.002       0.99
weight from a smoothed intensity     1.006    0.96-0.97
two-pass fill + plain blur           0.993       1.00
LDIC azimuthal affine (thesis 4.15)  0.809       0.87
```

`none` shows the rings at full strength, so the weight's shape doesn't
create them. The two-pass fill shows what plain blur does: a correct
collar value restores photometry to 0.993 while rings return to 1.00 —
a compensating error of the same shape.

0.22.25 knowingly ships two errors that cancel. Cost: the corona just
outside the limb reads ~25% low at 1.02 R here, worse elsewhere (8x on
a 360 mm Canon EOS 80D dataset). `ECLIPSEFORGE_FEATHER=taper` gives
correct brightness with the artifact visible.

**Ruled out, each by measurement:** the renderer's fixed-radius layer
blends; the feather; short-tier noise; per-tier correction by scale,
colour, radius, signal level; hard exclusion of the collar.

**Also found, not yet fixed:** adjacent tiers disagree by up to 3x
within 0-8 px outside the saturated region, ~1% beyond 8-16 px — likely
charge spill or veiling glare — and the merge lacks Druckmullerova's
LDIC term (thesis eq. 4.15), added in 0.22.26.

Render Warmth/Tint defaults remain at 1.0 (0.22.23).

## 0.22.24

TODO 1a was reopened with a rough rig, showing two tiers from the S1R II
dataset departing by -22% and +64% inside 1.25 R — too rough to build
on.

The pipeline now measures this itself: each tier's radial profile,
scene-referred, divided by the median across tiers, 1.0 to 3.0 R,
written to `tier_radial.json`. Verified on synthetic tiers with a
planted +60% near-limb error on tier 2, correctly named at 60%.

**Changes nothing in the output** — no correction applied, no layer
moves; a 0.22.19+ work directory is still reused. The goal is gathering
this across five datasets before deciding whether a radial per-tier
correction is warranted.

## 0.22.23

**Warmth and Tint default to 1.0.** 0.22.9 set them to 0.752 / 0.706,
assuming the far sky is neutral before these apply — not true on any
dataset measured:

```
set                     far sky R/G   B/G
560mm Canon EOS 450D        1.244    0.908
360mm Canon EOS 80D         1.206    0.880
250mm Canon EOS Rebel T7       1.116    0.757
600mm Panasonic S1R II      0.992    0.647
Sony                        0.553    1.328
```

The coefficients were fit to one eclipse, matching that target only
because its own sky chroma was folded in — they render violet
elsewhere. 1.0 / 1.0 leaves the choice to the sliders.

### Bisect: the feather is the sole cause

Three switches tested on the S1R II dataset:

```
change                       texture at 1.06 R (NAFE)   elsewhere
pedestal off                        1.00                 1.000
exposure exponent 1.0 (was 0.55)    1.03                 0.97-1.05
feather plain (pre-0.22.16)         0.36                 0.95-1.05
```

Only the feather moves the result, at 1.06 R; it also softens the limb
(ramp 6 to 8 px, fit rms 2.37 to 2.62, R 617.4 to 619.4 vs. tier
consensus 617).

```
m80-250 coherence at 15 px lag   1.06R  1.10R  1.25R  1.40R  1.90R
hdr_lum                           0.33   0.49   0.72   0.85   0.77
MGN                               0.30   0.41   0.54   0.68   0.72
NAFE                              0.25   0.44   0.62   0.70   0.65
```

Beyond 1.4 R the structure is radially coherent, consistent with real
corona; at 1.06 R, not.

### Per-tier radial error

Four raws, aligned with the pipeline's own shifts, scene-referred,
compared ring by ring:

```
tier / median-of-tiers   1.03R  1.06R  1.10R  1.15R  1.20R  1.30R  1.60R
1/60s                    1.000  1.008  0.999  0.994  1.035  0.999  0.996
1/30s                    0.947  0.912  0.912  0.911  0.775  1.001  1.050
1/20s                    1.216  1.216  1.226  1.359  1.640  0.925  0.978
1/13s                    1.000  0.992  1.001  1.006  0.965  1.008  1.004
```

Two tiers disagree with the others by −22% and +64% inside 1.25 R,
agreeing outside it. Saturation accounts for only 0-13% of the ring —
a per-tier error varying with radius, what TODO #1 was meant to fix
and 0.22.16 prematurely closed as unnecessary.

One raw per tier, manually estimated centre/radius; too large to be a
setup artifact, but not final.

## 0.22.22

A bisect run set `ECLIPSEFORGE_WALPHA=1.0` against a build predating
the switch. Silently ignored, it appeared valid, returning
bit-identical to default — every band 1.000, alpha 0.55, coherence
unchanged to four decimals — costing twenty minutes on an unintended
reproducibility check.

The run now logs any `ECLIPSEFORGE_*` variables it can see, with the
version reading them — a build too old to recognize a switch can't
print the line, signaling an outdated build. Running `efhdr --version`
first is the other half of this check.

## 0.22.21

The pedestal term is ruled out by measurement:
`ECLIPSEFORGE_NO_PEDESTAL=1` on the S1R II dataset, before and after
against the same rings:

```
                       with pedestal   without    ratio
hdr_lum m80-250 amp        43.938      43.936     1.000
   coherence 5/15/40    0.88/0.74/0.47  0.88/0.74/0.47
merged level 1.1-2.3 R                             1.000 .. 0.999
```

Identical in every layer and band. The texture is also not new — earlier
full-resolution exports already show it:

```
export             version   m80-250 rel-amp   coherence 5/15/40
render_16bit_2      0.21.4        0.047          0.29/0.23/0.10
render_16bit_5      0.22.5        0.075          0.46/0.36/0.13
```

Points to the exposure-exponent trial (0.22.5): it picks 0.55 here,
reporting "+106% coherent detail at 1.02-1.12 R," and tilts merge weight
toward the short, noisiest tiers. Its guard only checks mid/outer-shell
radial coherence, so it can add fine limb texture freely.

`ECLIPSEFORGE_WALPHA=1.0` holds the exponent at 1.0 and skips the trial.
Nothing changes by default; the report states when the switch is active.

## 0.22.20

0.22.19 removed the per-tier arcs; texture amplitude in the detail
layers is now 0.024–0.030, with fine structure left but no rings.

Measured on a 600 mm Panasonic S1R II dataset at 0.22.19:

* not the flat — merged texture correlates **−0.001** with the flat's
  own texture (7 px shift control: +0.001)
* not a per-tier boundary — no sharp radial features past 1.1 R
* not the filters — already in `hdr_lum`; FNRGF's version correlates
  **+0.999** with the merged luminance
* radially coherent, not noise — m80–250 coherence 0.897/0.760/0.483 at
  5/15/40 px lag vs 0.075/0.041/0.017 for m250–700

No measurement separates this from real coronal structure. Three merge
changes landed in a row — the pedestal (0.22.15), two feathers (0.22.16,
0.22.19) — too close to bisect, so this release adds switches instead:

```
ECLIPSEFORGE_NO_PEDESTAL=1     restores the pre-0.22.15 merge exactly
ECLIPSEFORGE_FEATHER=plain     the pre-0.22.16 feather (leaks into the clip)
ECLIPSEFORGE_FEATHER=masked    the 0.22.16 one (steps at the clip edge)
ECLIPSEFORGE_FEATHER=taper     the default, 0.22.19
```

Verified: `plain` leaks 0.46 of full weight into the clip at step 0.080;
`masked` leaks nothing at step 1.000; `taper` leaks nothing at step
0.133.

## 0.22.19

On a 600 mm Panasonic S1R II dataset at 0.22.18, rims appeared in the
corona in almost all layers: 0.22.16 fixed one thing and broke another.

A two-tier merge rebuilt from a 560 mm Canon EOS 450D dataset's raw
frames, identical weights, only the feather changed. "Step" is the
largest jump one tier's weight makes between neighbours:

```
version                  peak at   1.00R  1.02R  1.06R  1.20R   step
no feather (reference)    1.025 R  1.000  1.000  1.000  1.000  0.978
plain blur (<= 0.22.15)   1.095 R  3.690  0.674  0.856  0.977  0.027
blur then mask (0.22.16)  1.025 R  0.999  1.000  1.000  1.000  0.986
this one                  1.025 R  1.054  1.000  1.000  1.000  0.046
```

0.22.16 restored the weight at the boundary then zeroed it just past —
step **0.986 against the old 0.027** — so every tier's saturation
contour printed as its own arc, visible at exposure exponent 0.55
where every tier's boundary carries weight.

**The fix.** `den`, the blurred validity mask, is 1 inside the usable
region and 0.5 on the boundary, so a smoothstep on `(den−0.5)·2` reaches
zero exactly where the tier clips — continuous, still zero on clipped
pixels, profile unchanged from 1.02 R out. Step down to 0.046. Cost: a
clipping tier loses a feather-sigma-wide band inside its contour.

## 0.22.18

Reprocessing a 250 mm Canon EOS Rebel T7 dataset from Lightroom TIFF
exports instead of raw files produced a non-round sun/moon and
non-functioning flats — one cause, unflagged through nine pages of
plausible numbers.

Same nine frames, both ways:

```
                     tier factors      photometric links   limb fit rms
from the CR2s        0.944 .. 1.015    no systematic tilt     0.95 px
Lightroom TIFFs      0.251 .. 2.328    7 of 7 leaning up     10.32 px
```

A factor of 9.3 across the bracket where the raws give 1.08, and a limb
fit rms 3.3% of its own radius.

Everything assumes signal proportional to exposure time; a raw file
has that, an export does not, since Adobe applies an unpublished tone
curve `TiffFrame`'s sRGB inversion can't undo. Flats fail the same
way: falloff disagrees, 184.0% (TIFF) vs 105.6% (raw), same lens, same
30 frames.

A ratio-between-tiers test is blind to a gamma; instead, since a gamma
tilts every photometric link equally, `g` is fit from the links:

```
link_i = (s_i+1 / s_i)^(g-1)    ->    g = 1 + median( ln link / ln step )
```

Exact on planted values (1.00 → 1.000, 1.49 → 1.490, 0.80 → 0.800). The
TIFF run reads **g ≈ 1.46**. Past 8% from unity the run now flags the
frames as not scene-linear, names the tone curve, and notes the flat
and every number below assume linearity.

## 0.22.17

0.22.16 removed the pink rim on a 250 mm Canon EOS Rebel T7 dataset, but
the rings did not go, and limb variance came back 0.791 against 0.793 —
unchanged, as predicted. That statistic does not hold up.

Rebuilding the same coefficient of variation from the 250 mm dataset's
exported tiers, over 1.00–1.10 R, scene-referred, saturated pixels
excluded:

```
set          pipeline says   rebuilt from its own tiers
250 mm           0.791                 0.021
360 mm           0.037                 0.039
```

The 360 mm dataset reproduces to within 0.002; the 250 mm number does
not. Ring by ring, the 250 mm tiers agree in the near-limb rim to
**0.3–2.1%** — *better* than the 360 mm dataset's 1.1–2.8%.

`stack_variance` needs three tiers per pixel; between 1.00–1.10 R
**only two of the 250 mm bracket's nine tiers hold unclipped signal**,
so excluding clipped pixels keeps the dim tail of a blown tier and the
median is a biased remnant. The 360 mm dataset has five or six clean
tiers there, and never enters that regime.

The statistic now counts clean tiers and, below three, reports **not
measurable** instead of a contaminated number, naming the fix: a
shorter top-bracket tier.

The veiling-glare explanation from 0.22.11, and the 0.22.13 entry
calling the 360 mm dataset a confirming experiment, are both
withdrawn: both rested on 0.793 vs 0.037, a tier-count artifact, not
optics. The source comment stays, marked withdrawn — its ringing
measurements (2–5.3x the reference dataset in 1.02–1.15 R) are still
real and unexplained.

**The rings are real and their cause is once again open.** Excluded so
far: the merge weight (0.22.16), the limb fit (0.22.14), the pedestal
(0.22.15), tier-to-tier brightness disagreement (this entry).

## 0.22.16

On a 360 mm Canon EOS 80D dataset and a 560 mm Canon EOS 450D dataset
(2024), a pink rim appeared around the limb — a defect 0.22.14 had
merely stopped hiding. `gaussian_filter(wsat, 0.032*R)` blurs both
ways, handing weight to **clipped** pixels: weight goes as exposure
time, so a blown 2 s tier outvotes an unblown 1/1000 s one two
thousand to one.

A controlled two-tier merge rebuilt from the 560 mm dataset's raw
frames (1/125 s and 1/8 s, LibRaw-decoded, identical weights, only the
feather changed):

```
mode              peak at    1.00R   1.02R   1.04R   1.06R   1.09R
no feather         1.025 R   1.000   1.000   1.000   1.000   1.000
0.22.15 as shipped 1.095 R   3.690   0.674   0.748   0.856   0.900
0.22.16            1.025 R   0.999   1.000   1.000   1.000   1.000
```

Both raw frames peak at 1.02–1.03 R, falling monotonically outward; the
shipped merge peaked at **1.095 R** — brightness rising away from the
limb — with a bright overshoot and a suppressed band outside: the pink
rim.

Radial peak in the merged luminance, which should sit at the limb, in
every dataset: 560 mm (2024) 1.093 R / 5.4x suppression, 360 mm
1.083 R / 3.2x, 250 mm 1.057 R / 1.2x, 600 mm 1.030 R / 1.05x, Sony
1.020 R / 1.00x — longstanding, uncovered when 0.22.14 moved the disc
mask onto the real limb instead of 28 px outside it.

**The fix**: normalized convolution, already used to keep the lunar
disc out of MGN. Smooth the weight and validity mask together, divide,
re-apply — smooth where usable, full magnitude at the boundary, zero
wherever clipped. Applied to both feathered weights.

Also fixed: the pedestal progress line printed `1238016000000.0 sigma`
at a zero-jackknife; now states that in words.

## 0.22.15

A backlog item asked for a per-tier photometric correction varying
with radius. On three real datasets the thing that varies is a
**single additive constant** from the black subtraction — one number
fixes it.

Put each aligned tier back on the scene scale (divide by exposure and
photometric factor); over 1.5–3.5 R they should coincide and do not:

```
set                  tiers   scatter   with one offset   fitted P
360 mm                 12     39.20%        2.93%        +2.95 ADU (14-bit)
560 mm (2024)          14     36.43%        5.24%        +5.57 ADU
250 mm                  9      2.42%        1.91%        -1.97 ADU, shrunk to -1.31
```

Thirteenfold and sevenfold, from one number per dataset; the 250 mm
dataset, a different camera body, barely has one — the control. The
merge divides by `s·cal`, so a leftover offset P arrives divided by
exposure time — overwhelming on a short tier, invisible on a long one.

One shared number, not one per tier: fitting each its own offset moved
2.93% → 2.60% (360 mm) and 5.24% → 4.67% (560 mm) while the long tiers'
values wandered to −7 ADU — a black-level residual belongs to the
sensor and session, not the shutter speed.

P is shrunk toward zero by a leave-one-tier-out jackknife (synthetic
no-pedestal tiers apply −0.06; a planted +3.00 applies +3.10), and
refuses below four usable tiers. Applied after demosaic and **before
white balance**, and before the clipping test, since `cmax` is about
the sensor ceiling, not the scene.

Still open: the near-limb tier disagreement on the 250 mm dataset (limb
variance 0.793).

0.22.14 stated the merge amplifies near-limb azimuthal structure 1.5x
at 1.3 R. **Wrong**: that compared a compressed 1/8 s frame and a
2x2-binned raw against a full-resolution merge. Repeated on aligned
tiers, linearised, merged/best-tier runs **0.97/1.02/1.03/1.08/1.08**
at 1.15–2.80 R — no near-limb amplification. The 250 mm dataset is
different — merged/best-tier runs **0.25/0.29/0.35/0.46/0.60** at
1.05–1.40 R, losing three quarters of the real structure there, still
unexplained.

## 0.22.14

On a 560 mm dataset (2024), the moon mask was reported too large,
covering prominences: merged half-level fit **R = 553.0 px** against a
per-tier consensus of **525 px**, disc mask at 578.2 — twenty-eight
pixels of real corona under the mask.

Fifth dataset, every one biased the same way:

```
set                            consensus  range   merged fit    bias    ramp
600 mm                             617 px    7 px      619.1     +2.1      8
250 mm                              298 px    7 px      300.3     +2.3      9
360 mm                              456 px    2 px      470.4    +14.4     21
560 mm (2024)                       525 px   15 px      553.0    +28.0     28
560 mm (2024, 0.20.1 run)           525 px  226 px      587.3    +62.3     69
```

Always positive, tracking the merged limb ramp: a soft edge's 50%
crossing sits **outside** the true limb, and the merged edge is soft
since the Moon moves mid-bracket. Per-tier fits, on unsmeared tiers,
give the better estimate.

**The fix moves the merged fit's per-azimuth shape onto the tiers'
radius** — an offset, not a rescale, since the bias is a fixed pixel
count, not a fraction of R.

`R_consensus_spread` is a **p90–p10 range**, but the 0.22.6 test used
it as a standard deviation; p90−p10 = 2.563 σ, so "four sigma" was
really ten, and the 2024 dataset (15 px range, 28 px disagreement, 4.8
real sigma) passed, found only by inspection. Converted:

```
set                            diff   threshold   fires
600 mm                          2.1        10.9   no
250 mm                          2.3        10.9   no
360 mm                          14.4         4.6   YES -> R=456
560 mm (2024)                   28.0        23.4   YES -> R=525
560 mm (2024, 0.20.1 run)       62.3       352.7   no
```

The last row stays quiet correctly: those tiers disagree by 226 px, so
their consensus is not evidence of anything. R sets the radial profile
MGN, FNRGF, deband and disc mask all use, so it had to be right first;
0.22.14 opens a new cache family for it.

Also settled: the radial spokes are real, not a pipeline artifact — a
raw frame and `hdr_lum` share the same azimuthal structure (coherence
+0.958 vs +0.985, ring correlation +0.782 to +0.963), except **1.3 R,
where the merged modulation is 1.5x the raw frame's** (8.56% vs 5.79%).

## 0.22.13

A 360 mm Canon EOS 80D dataset, same eclipse as the 250 mm dataset, is
the natural experiment 0.22.11 was waiting on.

**Limb variance 0.036, against the 250 mm dataset's 0.793** — 22x
lower, below even the reference dataset's 0.075: the 250 mm
disagreement is optics-dependent, as predicted, and the photometric
correction backlog item is unblocked.

0.22.6 stated a disputed limb fit means the merge is "probably smeared
by a cross-tier alignment error." True of the 560 mm dataset, **false
of the 360 mm dataset**, which aligns to 1.36 px with 3 px of limb
spread and is still 14 px large. The bias tracks the **merged limb
ramp** alone across all four datasets:

```
set                consensus   fit     bias    ramp   align resid   limb spread
600 mm                617     619.1   +0.3%     8px      1.17px        8px
250 mm                298     300.3   +0.8%     9px      0.66px        3px
360 mm                456     470.4   +3.2%    21px      1.36px        3px
560 mm (2024, 0.20.1 run)   525   587.3  +11.9%    69px    512.09px      226px
```

**On a soft edge the 50% crossing sits outside the true limb** — a wide
ramp can come from a smeared merge *or* a genuinely soft limb, and
alignment residual and limb spread now tell those apart instead of the
warning picking one.

The disc mask sits at 489.3 px on a Moon whose tiers say 456 — 19 px of
real corona underneath; logged with three candidate fixes. And
`rim nan px` on this run, correct but alarming (no rim to measure when
tiers agree), now prints "none detectable — the tiers agree at the
limb."

## 0.22.12

Ten open items, scattered across changelog entries, source comments
and discussion threads, are now collected in `TODO.md` with the
evidence, uncertainty, and what would settle each — ordered by
expected value.

Top of the list: a **radially-resolved per-tier photometric
correction**. Current calibration is one scalar per tier, right only
if a tier's error is constant across the frame; a 250 mm Canon EOS
Rebel T7 dataset shows it need not be (limb variance 0.793 against
0.075/0.067/0.052 elsewhere, a 62 px disagreement rim against 4 px).
The structure exists in the literature: Druckmüllerová's LDIC (thesis
§4.1.4, eq. 4.15) fits an affine `k_i(φ)`, `q_i(φ)` per segment.

Written up but deliberately **not built yet** — it rests on a glare
hypothesis with one supporting observation, and a 350 mm dataset is
about to test it.

Rule for the file: an item leaves either by shipping, or by being
**rejected with a measurement**, filed beside the thing it rejects.

## 0.22.11

A 250 mm Canon EOS Rebel T7 dataset rings, and **not one alignment guard
fires on it**: network residual 0.66 px, limb spread 3 px, limb-fit rms
0.95 px, track scatter 0/0 — geometrically clean, so the 360 mm and
250 mm ringing are two different faults.

Measured on the cached layers, ring by ring (250 mm dataset / 600 mm
reference dataset, % oscillation of radial median about a 15-ring
smooth): at 1.02–1.15 R, MGN 19.4/9.9, FNRGF 10.9/2.5, NAFE 13.9/2.7,
inner 15.1/4.6 — **2 to 5.3× the reference dataset**, dropping to
near-parity by 1.8–2.6 R. The ringing lives **almost entirely in
1.02–1.15 R**, while the merged luminance itself oscillates only
**1.4%** — the filters amplify a real disagreement by about ten, not
invent one. **NAFE is the tell**: 5.3× worse despite using neither the
limb fit nor the disc mask — not geometry.

**Limb variance 0.793**, against 0.075/0.067/0.052 elsewhere, and a
62 px disagreement rim against 4 px: the tiers disagree about
*brightness*. Likely cause: veiling glare around a very bright limb —
a 240 mm f/5.6 kit zoom, not a 600 mm prime, which no flat removes.
Warned above 0.30, firing on the 250 mm dataset alone:

```
250 mm dataset   cov 0.793  rim  62 px   WARNS   [rings]
360 mm dataset   cov 0.067  rim 134 px   quiet   [rings — wrong limb circle, 0.22.6]
560 mm dataset   cov 0.052  rim  14 px   quiet   [broken alignment, 0.22.8]
600 mm dataset   cov 0.075  rim   4 px   quiet   [clean]
```

Three ringing reports, three different causes, each now identified by
its own measurement.

## 0.22.10

Continued monitoring turned up three more numbers in a 560 mm Canon
EOS 450D dataset's log, each diagnostic of total failure and unguarded.

**Per-tier limb spread**: 226 px on a 525 px Moon, 43% — warned above
8% of R.

**Limb-fit rms**: 32.77 px on R=587.3, 5.6% — the old filter allowed
anything under **8% of R**, so it passed a circle nothing like a lunar
limb. Warned above 2%.

**Lunar track scatter**: 0.63 px/s, 6 px across the bracket, scatter
about the line 26/39 px — describing the alignment error, not the
Moon's motion. Warned above 2% of R.

Against the four real datasets (limb spread/rms/track scatter): 560 mm
43.0/5.58/7.4% all WARN; 360 mm 3.3/1.15/0.7% ok; 250 mm 1.0/0.30/0.0%
ok; 600 mm 1.3/0.41/0.3% ok — each fires on the broken run alone. With
0.22.6 and 0.22.8 that is five independent checks on one failure that
otherwise reports "ready."

## 0.22.9

The `tint = 1.205` default is a straight **20% green multiply**; after
unit-luminance renormalisation the net is green +5%, red/blue −13%
relative — the green-grey sky cast in renders. `IMPORT_DEFAULTS`
already corrected this for imports; the raw path kept the by-eye
0.9/1.205 values from 0.12.0.

The far sky is neutral before `temp`/`tint` apply, so sky chroma
reduces to `R/G = temp/tint`, `B/G = 1/(temp·tint)`. A manually
corrected version measured off a JPEG export: sky R/G 1.064, B/G
1.884, giving `temp=0.752, tint=0.706` against the old 0.900/1.205 —
R/G ×1.425, B/G ×2.043 relative to green. It cools the corona too, why
`satur` is the other half; a taste call, now backed by a reference,
taking effect next render with no re-stack.

Adding 0.22.9 to the flat `CACHE_COMPAT` set silently made **every
pre-0.22.5 stack reusable again**, since `cache_ok` only checked both
builds were in the one set. Replaced with `CACHE_FAMILIES`, three
frozensets: 0.22.8 and 0.22.5 stacks stay reusable by 0.22.9; 0.22.4
and 0.20.1 do not, since the merge changed at 0.22.5.

## 0.22.8

A 560 mm dataset (2024 eclipse, processed at 0.20.1) failed completely
and reported "ready":

```
alignment network residual max 512.09px (half-res)      = 1024 px full-res
per-tier lunar limb spread: 226px (R_moon 525px)        vs a 4 px radius spread
limb half-level fit: R=587.3px, rms 32.77px             vs consensus 525px
merged limb ramp 69px -> disc mask margin 47.0px
1/4s  shift 108.4 px   "no usable signal in the correlation window"
0.5s  shift 117.8 px   "no usable signal in the correlation window"
prominence anchors: 1/15s:0 1/8s:0 1/4s:0 0.5s:0 1s:0
...
pipeline complete
ready
```

A residual of 1024 px full-res is **twice the lunar radius**. Five
tiers had no prominence anchors, three had no correlation signal, two
of those still producing 108 and 118 px shifts.

**Three fixes.** 1. A large residual is now a named failure: the links
contradict each other, and least squares returns the best compromise
between impossible constraints — no automatic repair, since which
link is wrong is unknown; names the three worst. Tolerance is 5% of
the seed lunar radius: this run's 512.09 px against 15.1 px FAILS
LOUDLY, versus 3.95/12.8, 0.88/7.1 and 1.17/17.3 px quiet elsewhere.

2. 0.22.6's limb check, built from the 360 mm dataset, is validated
elsewhere: here the merged fit is 587.3 px against a tier consensus of
525 px ± 4 — **62 px out, but only 11.9%**, so the old flat 15% rule
said nothing. The new relative test fires at 62.3 px against 16 px.

3. A misleading log line — "prominence anchors: 6 on the 1/60s tier,
11/11 tiers linked by them" — counted links, not tiers. Now reads
"giving 11 of 11 possible links."

## 0.22.7

Diagnosis this week repeatedly stalled on requests for raw eclipse
data — rarely something photographers want to share — letting the same
defects survive into the next release. Almost none of it needs the
raws: ringing, a wrong limb circle and a smeared merge are all visible
in per-tier numbers and radial profiles — a few thousand numbers, not
an image.

So every run now writes `eclipseforge_output/eclipseforge_diagnostics.zip`:

```
README.txt        what is in here, in plain words, for the user to read first
report.txt        the run report they can already read
geometry.json     the fitted lunar centre and radius
stats.json        the measurements behind the report
profiles/*.json   each layer as numbers -- median, p10, p90 ring by ring from
                  the limb outward, plus azimuthal cuts at 1.15, 1.6 and 2.4 R
thumbs/*.png      each detail layer at 512 px, greyscale
```

**0.48 MB** with three layers cached, about 1.5 MB for a full run.
Excludes raw frames, full-resolution output, colour images,
tone-mapped results, EXIF, and file paths beyond the folder name.
Verified: on a well-fitted run MGN's ring-by-ring median oscillates a
smooth 6.3% of its trend.

## 0.22.6

Logs and exports for two 0.15.1 datasets reported ringing artifacts in
both, worse at 250 mm, with the 360 mm details appearing stretched as
if misaligned. The 360 mm log had every number needed to catch it:

```
lunar radius consensus across tiers: 456px (spread 2px)
limb half-level fit: ... R=470.2px, rms 5.43px   (seed fit said R=482)
merged limb ramp 25px -> disc mask margin 22.5px
alignment network residual max 3.95px (half-res) = 7.90 px full-res
disagreement rim 134 px wide just outside the limb
```

Tiers agree on the lunar radius **to within 2 px**, the merged fit came
back **14 px larger**; the old flat 15% threshold said nothing at 3.1%,
but against a 2 px spread, 14 px is seven sigma — the merge is smeared
by a cross-tier alignment error, inflating R.

**Very likely the cause of the reported ringing**: a circle 14 px too
big prints concentric arcs into everything built on it (MGN, FNRGF, the
deband, Pellett geometry), matching ringing reports in all three at
once.

The test is now "large compared with how well the tiers agree," with a
1%-of-R floor: 600 mm (diff 2.1, threshold 28) and 250 mm (2.2, 28)
stay quiet, 360 mm (14.2, 8) fires — discriminating on all three real
datasets. It warns rather than overrides, since the right repair is
the alignment, not the circle.

Also fixed: both runs show `sky gradient removal skipped ([Errno 22]
Invalid argument: ...hdr_rgb.npy)`, the Windows memory-map file lock
fixed after 0.15.1. The 250 mm run also hit the MGN divisor bug fixed
in 0.22.2 (**5 scales**, 1.4/2.4/4.8/9.7/19.3 px), so that layer came
out at 0.83x amplitude.

## 0.22.5

Exported aligned tiers unblocked two items untestable from cached
products alone.

`wsat` was a high-end shoulder alone — full weight below the knee,
even for pixels holding only read noise, against Druckmüller, Rušin &
Minarovjech 2006 requirement (b): "w = 0 in substantially
underexposed **or** overexposed parts." It matters more since
`_pick_weight_alpha`'s α = 0.55 gives **short** tiers more relative
weight, buying +106% coherent detail at the limb, paid for in the
outer field where those tiers hold no signal — a per-pixel floor buys
the limb without paying out there.

Measured by rebuilding the merge from all 14 exported tiers:

```
shell        as shipped   floor 0.002  floor 0.005  floor 0.01
1.2-1.8 R      0.3189       0.3189       0.3189      0.3190
1.8-2.6 R      0.1162       0.1150       0.1140      0.1129
2.6-3.4 R      0.0299       0.0289       0.0282      0.0276
3.4-4.2 R      0.0135       0.0127       0.0123      0.0119
```

Nothing lost where the corona is bright; −2.8%, −7.7%, −11.9% going
outward. Shipped at **half** the best measured value, since
quarter-resolution sRGB exports can't separate faint structure from
noise. `_MERGE_FLOOR = 0` restores earlier builds exactly.

A LinearFit-style intercept for the photometric calibration, previously
recommended on Hill's slide, tested against the tiers, pair by pair:

```
pair             predicted   measured k   k/pred
1/20 -> 1/13       1.597       1.5856      0.993
1/13 -> 1/8        1.611       1.5967      0.991
1/8  -> 1/5        1.602       1.5821      0.988
1/5  -> 0.3125     1.619       1.6016      0.989
0.3125 -> 0.5      1.630       1.6128      0.989
0.5  -> 1.6        3.412       3.4429      1.009
```

**Gain-only calibration reproduces every tier ratio to within 1.2%**;
the affine intercept buys almost nothing here. No pedestal in this
data — an earlier "large pedestal" reading was **PIL misreading the
16-bit RGB TIFFs** (256 distinct values where `imageio` reads 65,530
from the same file), caught by a second library.

## 0.22.4

Two promising leads from a literature audit were built and measured.
Neither ships; the negative results are recorded in the code so they
are not retried.

The thesis calls ACHF the best structure-enhancement technique for
white-light eclipse corona, the one method this codebase lacked.
Implemented exactly — a Gaussian in radius crossed with a Gaussian in
**arc length** along a Sun-centred circle. Against an isotropic kernel,
identical normalisation and ladder, only kernel shape varying:

```
shell          isotropic     ACHF      change
1.02-1.15 R      0.0701     0.0770      +10%
1.15-1.40 R      0.0610     0.0607       -1%
1.40-1.80 R      0.0438     0.0450       +3%
1.80-2.60 R      0.0465     0.0491       +6%
2.60-3.40 R      0.0530     0.0554       +4%
```

Indistinguishable side by side, and 390 s against 21 s at half
resolution. **Why**: the kernel is anisotropic only where σ is a large
fraction of r; at every scale in the ladder the two are almost the
same object, separating only at several hundred pixels — large-scale
structure, not the fine detail ACHF targets. Rejects the kernel shape
as a drop-in here, not ACHF as Druckmüller ships it.

A per-annulus gate fading RHEF out where an annulus holds no signal,
via radial coherence, was built and tested. Coherence is **0.99 at
1.5 R, 0.94 at 2.0 R, 0.78 at 2.5 R, still 0.55 at 3.0 R** — open where
the grain is worst. Wrong diagnosis: a rank transform **discards
amplitude**, stretching a 1% real modulation and a 1% noise
fluctuation across the same range; no presence gate fixes a signal
that is merely faint. RHEF stays opt-in at 0.

## 0.22.3

Two papers the previous audit could not check are now available:
Druckmüller 2013 (ApJS 207:25) and, inside the IWCIA proceedings at
printed p. 262, the variable-neighbourhood paper itself. Every
previously unverifiable claim now has a verdict.

**Confirmed, verbatim:** σ's band is ⟨2σ_A, 12σ_A⟩ (2013 p.3; 2014
p.268); w runs ⟨0.05, 0.3⟩ (2013 p.2; 2014 p.265); the n = 129 kernel
is twelve summed Gaussians, σ_m = 2^(m/2), m = 1…12, c_m = 1 (2014 eqs.
14–16); the Fig. 2 "loss of contrast near lunar edge" caption is real;
every equation nafe.py maps is correct.

**Ours, not theirs, now labelled:**
- `NAFE_K = 128` — neither paper bins anything (a Kronecker delta over
  actual pixel values, 2014 eq. 6); the cited "several thousand
  discrete pixel values" line argues *continuous* is safe, not that
  128 is right. An internal sweep only.
- `NAFE_EPS = 0.10` — crisp 0/1 in the paper (2014 eq. 12), but **no
  numeric value published**: "must be found experimentally."
- **Log input and `NAFE_FLATTEN_R`** — both papers assume linear
  dependence on coronal emission, no gradient removal, no log, no
  masking specified. Both departures now stated at the call site.

**Citation corrected:** 2013 is **ApJS 207:25**, not ApJ 775, 88. Fixed
in nafe.py, report.py, README.

A previous entry called every headline detail layer designed for a
different regime than white-light eclipses. **Wrong about NAFE-VN**:
the 2014 paper targets this data in its abstract, uses an eclipse
image (Fig. 2), and publishes a worked white-light result (Fig. 5) —
the earlier claim came from the survey paper's taxonomy filing *NAFE*
under EUV tools. ACHF remains the one method the literature calls best
for white light that this codebase lacks.

## 0.22.2

A literature audit checked the code against primary sources: the MGN
and NAFE/ACHF papers, the FNRGF chapter of Druckmüllerová's thesis, the
Hill slides, and the alignment chapters.

Eq. 5 divides MGN's multiscale sum by **n**, the scale count; the
implementation divided by `sum(gains)` — harmless at six scales (5.874
vs 6), but `scale_ladder()` returns fewer scales below the resolution
floor while `gains` stays six long, so the divisor stayed 5.874,
under-scaling the layer:

```
scales used     5       4       3
amplitude     x0.83   x0.66   x0.49
```

A 240 mm lens lost up to half the MGN layer, absorbed silently by
`mgnContrast`; six-scale datasets were unaffected.

`detail.py` divided `NAFE_NEIGH_R * R` by `NAFE_GRID`, and `nafe_vn()`
divided by `grid` **again**: applied neighbourhood was 0.13 R / 8 =
0.016 R — 10 px where docs said 80 px. **The bug was the better
setting**: a true 0.13 R scores higher (0.194 → 0.939 at the limb,
coherence 0.53 → 0.75) but looks blobby where 0.016 R was filamentary.
Double division removed, constant restated as 0.016; behaviour
unchanged to within 2%. The sweep that "justified" 0.13 ran through
the bug and actually swept 0.0075–0.064 R.

Labels only, no behaviour change: phase correlation is unused
(`normalization=None`) though the report cited Druckmüller 2009 for
it; `NAFE_GAMMA` is **dead**, and `multiscale_blur`/
`value_neighbourhood_weight` have no callers; `NAFE_EPS`'s "rank
units" comment now reads "level units."

Verified, not changed: MGN's `k=0.7`, `γ=3.2`, kernels and Fig. 4 gains
are faithful to the paper, though the source comment's rationale was
backwards. The paper says nothing about masking an occulting disc; the
normalized convolution here is an extension.
</content>

## 0.22.1

The detail balance slider did nothing when a work directory had no
`mgn_fine.npy`: correctly guarded off, but with no indication it was
inert — the same failure mode as 0.21.0's blank preview.

Sliders depending on an optional cached layer are now greyed, labelled
"— re-stack needed", and forced to 0: `detailScale` (needs
`mgn_fine.npy`) and `rhefMix` (needs `rhef.npy`).

`mgn_fine.npy` could in principle be rebuilt from `hdr_lum.npy` without
a re-stack, but only by reproducing the normalisation span and
prominence mask exactly, or the blend is silently wrong. Not
implemented.

## 0.22.0

The rendered corona looked bolder and rougher than intended, lacking
delicate structure — traced to an unquestioned line in the per-scale
gains.

The gains, `0.907, 0.976, 0.994, 0.998, 0.999, 1.0`, are Morgan &
Druckmüller's noise-normalisation constants (Fig. 4), not amplification
factors — all six scales get equal say, tilted toward the coarse end.
Hill's talk instead combines masks `1 : 0.6 : 0.2 : 0.1`, finest
weighted ten times the coarsest, and won a four-way comparison despite
scoring worse on amp×coh in every shell (0.0535 → 0.0300 at the limb,
−44%), since the metric rewards radially coherent structure, which the
coarse scales carry. Visual assessment took precedence.

### How it is exposed

Storing the three fine MGN scales as their own layer turns any ladder
constant within a group into a render-time blend. Measured fine/coarse
energy ratio:

```
as shipped            1.09
detailScale 0.50      1.72
detailScale 1.00      3.18
Hill's true gains     3.32     <- indistinguishable from 1.00 side by side
```

New slider, **Detail balance (fine →)** (`detailScale`, default 1.0); at
0, every earlier build is reproduced exactly, costing one extra MGN pass
and one cached layer.

A re-stack is required, since `mgn_fine.npy` cannot be reconstructed
from the cache; older work directories still load and render as before.

Still open: the background is green-grey rather than neutral, and the
corona stops reading at about 2 R — tone, not detail, and not addressed
here.

## 0.21.4

RHEF now defaults to 0. Shipping it at 1.0 in 0.21.0 made the picture
worse, visible immediately on review.

RHEF forces a flat histogram in every annulus — correct where an
annulus holds corona, but past about 3 R in a 600 mm frame, where it
holds only sky, it stretches noise across the full tonal range and
fills the outer field with grain manufactured from nothing. The amp×coh
score missed this: rank-equalised noise plus the layer's own 1.6 px
smooth still reads as moderately coherent, so score rose while the
image degraded — Gilly & Cranmer's method targets EUV disk imagery
where every annulus carries signal, but an eclipse frame runs out of
corona long before it runs out of frame.

The layer is still worth a drag between 1.6 and 3 R, where it measured
+28% against MGN. Default-on needs a per-annulus SNR gate zeroing the
weight once an annulus is noise.

### Amplification factors: an open question the metric can't settle

Hill's ladder (`1 : 0.6 : 0.2 : 0.1`, finest weighted ten times the
coarsest) is the opposite emphasis to the gains in use here (`0.907,
0.976, 0.994, 0.998, 0.999, 1.0`, Morgan & Druckmüller's
noise-normalisation constants, Fig. 4). Measured on the reference
dataset:

```
ladder                  1.05-1.30  1.30-1.80  1.80-2.60   fine/coarse
ours (M&D noise norm)     0.0535     0.0602     0.0425        1.09
Hill 1:.6:.2:.1           0.0300     0.0195     0.0197        3.32
Hill, gentler             0.0425     0.0328     0.0257        2.11
steeper than Hill         0.0241     0.0146     0.0175        4.26
```

Hill's ladder triples the fine/coarse ratio (1.09 → 3.32) but lowers
amp×coh in every shell, since the score rewards radially coherent
structure, which the coarse scales carry — an emphasis choice the
metric can't adjudicate. Not implemented: the gains apply when
`mgn.npy` is built, so a slider needs a second MGN layer or the
per-scale terms stored separately.

## 0.21.3

`rhefMix` was registered under the slider group `"detail"`, but the
valid groups are `structure`, `inner`, `prom`, `color`, `tone`, `ring` —
there is no `detail`.

The builder does `$("g-"+grp).appendChild(div)`, so the wrong string
threw a TypeError before any sliders were built — no controls, no
preview, no visible error, though the server logged normally, since the
failure was entirely client-side.

Two changes: the group is now `structure`, next to FNRGF share and
Clarity, and the builder warns and skips an unknown group instead of
throwing. A build check now confirms every group id used by a slider
exists as a div — a syntax check can't catch this, since the file
parses fine and `getElementById` returning null is a runtime fact.

No pipeline change; a re-stack is not needed. Reload the page.

## 0.21.2

Druckmullerova's own Delphi FNRGF implementation (`.pas` sources) was
reviewed against `fnrgf_robust`, line by line. No change needed: same
core, `(I - Ave(r,phi)) / Dev(r,phi)`. Differences from the reference
are equally or better motivated:

| | reference (`ImgProc.pas`) | ours |
|---|---|---|
| azimuthal fit | mean/sd in `SegmentCount` bins, Fourier fitted to the bin values | Huber IRLS directly on 1440 samples — no binning quantisation |
| high orders | `Atte[series,k]` per order, user-set, **default 0** (i.e. plain NRGF) | ridge `1e-3·m²`, always on, coverage-matched order |
| output | `norm(Input) + MixRatio·norm(Mask)` | `fnMix` in the renderer, same job |
| noise | `Noise_AddVar` added inside the deviation | *(absent)* |

The noise term targets a divisor sigma that is really grain, but
measured, it gives no improvement: coherence is unchanged to three
decimals in every shell (0.506, 0.452, 0.533, 0.624 before and after),
while amplitude falls 25% in the outer shell — a uniform rescale, and
here `fnCompress`/`fnMix` already own that role. Not added.

### A defect found in `_fine_structure`

The metric divided every radial column by its azimuthal median —
correct for luminance, catastrophic for FNRGF's `(L-mu)/sd`, whose
median is ~0 by construction, so the guard `max(median, 1e-9)` divided
by about 1e-9 and reported an amplitude of 2.6e8, unflagged.

The normalisation branch is now chosen by the data: a column whose
median dominates its own spread is divided as before; otherwise the
median is subtracted and scaled by the shell's robust spread. Every
layer already measured is bit-for-bit unchanged — the four MGN shells
from 0.21.0 still read 0.0535, 0.0602, 0.0425, 0.0380. A measurement
fix only; no pipeline behaviour changes.

## 0.21.1

Hill's guidance says the Moon, prominences and stars should read zero
in the mask. Partial (normalised, incomplete) convolution already
estimates every masked filter's local mean here, which is why the limb
has never shown Hill's black ring; what was missing was what went in
the mask — previously the Moon alone.

Prominences matter more than their area suggests: S, the local standard
deviation, is MGN's divisor, and a prominence — brightest thing in the
frame, hard against the limb — sets the sigma nearby coronal pixels get
divided by.

Measured with the limb split by azimuth into rays with and without a
prominence (control, same run):

```
pct grow  mask %ring   1.02-1.15 near/away   1.15-1.40 near/away
 99   2      1.60%       +6.8% /  +0.2%        -0.0% / +0.2%
 99   5      2.63%      +10.8% /  +0.2%        -0.2% / +0.3%
 99   8      3.73%      +20.2% /  +0.3%        -0.3% / +0.5%   <- shipped
 98   5      5.23%      +24.9% /  +0.3%        +0.6% / +0.4%
 97   5      8.00%      +34.1% /  +0.4%        +4.0% / +0.3%
 99  12      5.27%      +31.8% /  +0.7%        -0.5% / +0.9%
```

Only the near column moves, confirming this is prominences, not masking
more in general. It does not saturate, so mask size is a trade, not an
optimum, and 99/8 is conservative: an earlier reading of +26.9%/+24.9%,
using an ad-hoc mask about four times larger, does not reproduce below
8% of the ring.

Only MGN's statistics use the tighter mask; prominence cores come out
flat 0.5 there, and the gate (100% of the masked pixels) supplies
structure via promdet instead. Stars aren't masked — the same argument
applies, but no available dataset contains one.

## 0.21.0

Addresses grain outside the corona, based on Gilly & Cranmer,
*Visualization of High Dynamic Range Solar Imagery and the Radial
Histogram Equalizing Filter*, Sol. Phys. 300, 174 (2025): every pixel
becomes its percentile rank within a one-pixel annulus (eq. 1) —

```
I_out[A_i] = rank(I_in[A_i]) / N_{A_i}
```

— no kernel, no scales, no parameters.

With no spatial kernel it can't compete with MGN near the limb; with no
local sigma to divide by, it beats MGN badly further out, where that
sigma mostly measures grain. Both given the same mild smooth, scored as
amp×coh:

| shell | MGN | RHEF | |
|---|---|---|---|
| 1.05–1.30 R | 0.0531 | 0.0318 | −40% |
| 1.30–1.80 R | 0.0628 | 0.0372 | −41% |
| 1.80–2.60 R | 0.0492 | 0.0629 | **+28%** |
| 2.60–3.40 R | 0.0424 | 0.0570 | **+34%** |

and radial coherence alone: 0.489 → 0.752 at 1.8–2.6 R, 0.349 → 0.470 at
2.6–3.4 R. The renderer crosses over on radius rather than choosing;
weight and crossover were picked by measurement:

```
rhefMix   1.05-1.30  1.30-1.80  1.80-2.60  2.60-3.40
  0.4        +0%        -0%        +8%       +11%
  0.8        +0%        -1%       +32%       +36%
  1.0        +0%        -1%       +46%       +50%

crossover from 1.3 R    -0%       -15%       +48%      +50%
crossover from 1.6 R    +0%        -1%       +46%      +50%
crossover from 1.9 R    +0%        +0%       +32%      +50%
```

New slider **Outer field (RHEF)** (`rhefMix`, default 1.0) and an RHEF
layer view; at 0 every earlier build is reproduced, crossover fixed at
1.6–2.2 R.

No re-stack needed: RHEF is a single `lexsort` of luminance already on
disk, cached in seconds on first load, excluding the occulted disc from
each annulus's ranking. Measured on one dataset.

### What this replaced

Polar-domain convolution was tried first and rejected: an earlier
+155%/+344% gain, re-scored as amp×coh with proper azimuthal sampling,
turned out to be mostly the warp's own interpolation smoothing — the
polar kernel loses in every shell (−46%, −42%, −27%, −10%).

## 0.20.3

0.20.2 had the wrong diagnosis: it attributed a scale mismatch to the
contact frame being shot at a different focal length, but EXIF shows
600 mm for both — the ratio behind it came from circle-fitting a
crescent, which is not a circle. Retracted.

Measured on the third-contact frame against the composite it was loaded
into:

| | distance from the composite lunar centre |
|---|---|
| the crescent, in the frame as shot | **588 px** |
| the composite's lunar limb | 619 px |
| the crescent, after `prepare_contact` | **190 px** |

As shot, with no registration at all, the frame is already within ~30
px of right, as a tracked sequence should give; `prepare_contact` then
moved it 399 px up and left, putting the crescent well inside the lunar
disc. The cause: `fit_limb` was written for a totality frame, a dark
disc inside a corona. A contact frame is not that — 99.4% of it is
below the noise, dark sky with one blazing crescent — so a limb finder
fits the crescent's arc, the only edge in the picture.

The fit is now checked before it is obeyed: past a quarter of the
lunar radius, or a rescale beyond 0.8–1.25, the frame is overlaid as
shot instead, and the log says so and points at the ring sliders.

Kept from 0.20.2: the log no longer reports `auto-scale x1.0000` on a
silent give-up, alignment is checked against the lunar limb, and the
wider ring sliders stay (0.80–1.25, ±200 px).

Reload the contact frame to pick this up; the composite stack is
untouched.

## 0.20.2

A composite export with a diamond ring showed the ring in the wrong
place. The export was not at fault — the diamond lands within 1 px of
where the ring parameters put it — the cached contact layer was already
wrong.

Measured on the stored `contact_rgb.npy`: the ring arc circle-fits to R
= 553 px centred 327 px from the composite disc, whose limb is at R =
619 px — the crescent sits inside the lunar disc, impossible since the
Moon hides the rest of it. A registration failure, not a taste
question.

`prepare_contact` computed the scale it needed, `geo.R / R`, applying it
only when it fell in 0.9–1.1 and silently falling back to 1.0 outside
that band. The ratio needed was 1.120, just outside; 553 x 1.120 = 619,
so the refused scale accounts for the error exactly.

Three changes: the band is now **0.5–2.0**, since a different focal
length is normal; the log states any refusal in full; and after
alignment the bright arc's distance from the disc centre is
circle-fitted and reported, flagging the frame as unregistered if
inside the lunar limb.

The ring sliders were also too narrow: **Ring disc size** now spans
0.80–1.25 (was 0.96–1.04) and the two offsets ±200 px (were ±30).

Reload the contact frame to pick this up; the composite stack is
untouched.

## 0.20.1

Review of 0.20.0 found the pink prominences visually wrong. The
measurement stands — green and blue are the only channels with
headroom, and structure does appear there — but a technically correct
mechanism can still look wrong, and visual judgment takes precedence
here.

`promChroma` now defaults to 0, range -1.5 to +1.5: positive values make
dense material go pinker (the previously shipped, rejected behaviour),
negative values push it toward deeper red. Skipped entirely at zero.

## 0.20.0

Setting Prominence detail to 1.0 after 0.19.0 produced no visible
change. Measurement found the cause.

Inside a prominence the red channel sits at 255 across 100% of the
bright core (green ~59, blue ~51) — the maximum channel everywhere, so
the hue-preserving highlight knee sets all three channels from `ms` and
the chroma ratio, and luminance detail cannot reach the render at all.
The same promDetail change measures 59 levels before the knee, 13
after; the 0.19.0 driver swap measures 1.5 levels on screen.

Green and blue have roughly 200 unused levels, so the prominence's own
structure can go there at no cost in red:

| | G spread (p10–p90) | agreement with H-alpha in G |
|---|---|---|
| promChroma 0.0 | 49–79 levels | 0.368 |
| promChroma 0.6 | 57–105 levels | **0.596** |

and on the two smaller prominences 0.565 → 0.762 and 0.740 → 0.827. On
screen, a mean of 13.6 levels and p90 of 27, against 1.5 for luminance.

New slider, **Prominence colour detail** (`promChroma`, default 0.6). At
0 every earlier build is reproduced exactly; outside the gate the
factor is exactly 1 (residual 1.2e-07, float rounding). The effect
reads as denser prominence material going pinker rather than redder —
stylistic, hence the slider. Measured on one dataset, three
prominences; a 0.18.0 cache is reused as-is.

## 0.19.0

0.18.0 built a prominence detail layer and blended it into `det`, but
left the contrast term — which gives a detected prominence its texture
— driven by the inner-corona layer, the best option before the new
layer existed. Scored as correlation with the H-alpha red channel's own
fine structure on each prominence's bright core (600 mm Panasonic S1R
II dataset, 14 tiers):

| prominence | inner-driven | promdet-driven |
|---|---|---|
| az 226 (the large one) | 0.714 | **0.828** |
| az 44 | 0.357 | **0.361** |
| az 10 | 0.143 | **0.284** |

No prominence scores worse. The positive bias stays at 0.30: lowering
it to 0.10 scores 0.828 against 0.807 on the large prominence, not
enough to justify changing overall prominence brightness.

Scoped by construction: `prom` is zero outside the gate, so a work
directory with no `promdet.npy` falls back to the old driver and
renders bit-identically. Measured on one dataset. A 0.18.0 cache is
reused as-is — this release touches only `render.py` and `gui.html`.

Two things unfixed: the disc mask is not covering the prominences
(mask weight runs 0.73–0.92; an earlier 83%-hidden reading sampled the
whole gate, mostly inside the disc) — and 83–91% of the large
prominence's bright core still renders above the highlight knee, since
the envelope pins at the 99.97th percentile and a prominence is only
5,870 px of 43 million. Raising highlight compression does not recover
it.

## 0.18.0

Prominence interiors rendered flat, though short exposures show the
structure is present. Two causes, on the reference dataset's largest
prominence (225 deg, R/GB 16.4):

MGN's normalisation window clips them — `hi` is the frame's 99.95th
percentile, correct for a corona but too low for a prominence:

    red channel        37% of the prominence hard-clipped at xn = 1.0
    merged luminance   21% clipped, the rest squeezed into the top 6% of range

Giving the layer its own window, taken inside the gate, drops the
clipping to 2%. MGN is also the wrong filter for a compact bright
feature: dividing out the local standard deviation flattens a
prominence whose interior variation is its own local sigma. Correlation
with the red channel's fine structure:

    MGN, frame window (what the corona layers do)        0.037
    MGN, gate window, corona scales                       0.375
    MGN, gate window, fine scales                          0.420
    plain multiscale unsharp, gate-scaled                  0.940   <-- shipped

    existing layers on the same measure:  inner 0.317, MGN 0.274, merged 0.262

It also uses red, not luminance, which weights R at 0.2126, costing
more than half the structure before any filter sees it (1.000 to
0.417). So `promdet.npy` is a multiscale unsharp of log(red) from the
H-alpha tier, scaled by the prominence gate's spread, blended only
where `prom` fires. New slider, **Prominence detail**, default 0.7; at
0 the render matches 0.17.0.

Measured end to end, correlation of the render's detail with ground
truth inside the prominence:

    promDetail 0.0 (0.17.0)   0.368
    promDetail 0.4            0.535
    promDetail 0.7            0.633
    promDetail 1.0            0.682

Two approaches that looked correct were ruled out by measurement:
partial convolution with the prominences masked out (Druckmuller's
method) makes things worse (0.0850 to 0.0561, 0.0519 dilated), since
excluding a bright feature drags the local mean down and the arctan
saturates; correlation against a fixed reference is used instead of an
earlier, invalid ratio comparison.

An older work directory with no `promdet.npy` omits the key entirely
rather than substitute a flat 0.5, which would erase detail instead of
doing nothing.

## 0.17.0

Short-exposure frames appeared to carry detail lost in the merge.
Measured by comparing `short_lum` (the four shortest tiers) with the
merged HDR, separating structure from grain by radial coherence:

    shell          source        amp     coh    amp*coh
    1.01-1.10 R    merged     0.0326   0.995     0.0325
    1.01-1.10 R    short      0.2101   0.993     0.2086
    1.30-1.80 R    merged     0.0212   0.976     0.0207
    1.30-1.80 R    short      0.0576   0.186     0.0107
    1.80-2.60 R    merged     0.0135   0.874     0.0118
    1.80-2.60 R    short      0.1560   0.013     0.0020

At the limb, short tiers carry 1.9x more coherent fine structure round
the disc (1.38x quietest sector, 3.22x busiest). The merge loses it
because `w = s * rolloff` weights by exposure time — optimal for photon
noise, blind to glare-smearing near a bright edge — so the longest
not-quite-saturated tier there carries 67x the weight of the sharpest.

Streamers are different: beyond 1.3 R the short tiers are noise
(coherence 0.186, then 0.013), nothing to recover; the a-trous denoise
is not responsible (coherent structure retained is 99.9%, 100.7%,
107.4% across the outer shells).

The merge now has one knob, `w = s**alpha * rolloff` (alpha=1.0
default), set by a measurement trial — four alphas, half-res merges,
seconds — used unless another value wins. The guard is coherence, not
raw score, since a noise-dominated tier raises amplitude in every shell
(the noise case scored +41% at the limb): as alpha falls to 0.55 there,
coherence goes 0.933 -> 0.438 at the limb and 0.398 -> 0.225 mid-field,
while the genuine glare case barely moves. An alpha costing more than
0.05 of coherence anywhere is refused.

Seven cases, all passing, two failing first:

    glare-smeared long tiers            -> tilts to 0.70, +25% at the limb
    every tier equally sharp            -> stays at 1.00
    short tiers replaced by pure noise  -> refused  (failed first: chose 0.55)
    no lunar track                      -> stays at 1.00
    two-tier bracket                    -> stays at 1.00
    shells off the frame                -> refused  (failed first: chose 0.55)
    alpha 1.0 vs the old code           -> bit-identical over all 14 tiers

The off-frame case matters: a long lens filling the frame leaves no
outer shell to guard with, so the trial declines rather than deciding
blind. (The fine-structure metric itself was fixed to sample at one
pixel rather than a fixed count, which had let bilinear interpolation
make noise look coherent on a small disc.)

## 0.16.3

Prominences were being lost, including the largest one. On the
reference dataset the mask sat at R+25.1 px; recovering them required
setting Disc mask trim to -40.

The margin came from a metric that cannot measure an edge:
`limb_transition_width` took its 100% reference as `percentile(v, 95)`
over 0.75-1.35 R, but the inner corona is still climbing well beyond
the limb — it peaks at R+27 here — so "80% of it" lands far outside
the lunar edge. Proof — the answer scaled with the reference:

    reference at   R+5   R+10   R+15   R+20   R+30
    median width   6.5    9.8   12.5   16.0   20.0 px
    p90 width      9.0   11.5   16.0   21.0   28.0 px

A real edge converges as the reference moves outward; this one grows
without bound — it was measuring the corona's radial gradient, not the
Moon.

The reference is now taken just outside the edge, at R+5 px. Measured:

                        shipped     corrected
    ramp (p90)           28.0 px      11.0 px
    disc mask margin     25.1 px       9.9 px

15.2 px of annulus recovered all the way round — 23 arcsec at this
plate scale — exactly where the prominences were. On an imported stack
whose limb is genuinely sharp, the corrected metric returns the same
margin as before (9.9 px).

Still outstanding: the fitted R is itself about 8 px too large (merged
disc floor at R-7.8 px; per-tier consensus 617 px against the merged
fit's 622.9), likely the same bias in `fit_limb_rays`. Left for a
separate release; per-tier fits, tighter by 7 px, suggest fitting the
disc on a short exposure instead.

## 0.16.2

0.16.0 replaced the per-tier photometric link with a ratio of sums and
no data-dependent selection. On the reference dataset this made things
much worse and reintroduced a bright rim:

                       0.14.5    0.16.1
    1/4000s             1.273     0.571
    1/2000s             1.061     0.746
    1/500s              0.965     0.878
    1.6s                1.158     1.326

    disagreement rim    2 px      128 px
    limb variance       0.072     0.255

That 128 px of tier disagreement outside the limb is visible in MGN,
NAFE and FNRGF. The first three links flipped from 0.833 / 0.910 /
1.009 to 1.307 / 1.177 / 1.091 — an over-correction opposite the bias
it was meant to remove.

The validation was flawed: the synthetic tier pairs had a corona
filling most of the frame, but a real wide-field dataset is mostly sky
(8100x5357 px, corona reaching 4 R), so the sums were dominated by sky
area, with nothing at the short end to take a ratio of. Neither
estimator is correct; thresholding a tier against its own noise biases
it upward, one-sidedly and compounding — why a 25-tier FITS dataset
reads 23 of 24 links below 1.000.

What stays: the per-link residual line and systematic-lean warning from
0.16.0, diagnostic only; the FITS colour balance, memmap fix and
orientation control are untouched.

A replacement must take the ratio over an annulus just outside the
limb, chosen geometrically, and validate on both datasets before
shipping. Cached products from 0.16.0 and 0.16.1 are not reused.

## 0.16.1

0.16.0 placed the FITS row-order override at read time, making it part
of the cache key, so changing orientation cost a full re-stack — six
minutes on the 200 mm dataset — for something no measurement depends
on. Since a flip is trivially reversible, it is now applied at the end.

Reasoning: the limb fit is a circle, and MGN, FNRGF, NAFE and Pellett
are all radial or tangential about it, so turning the finished picture
any way up changes nothing downstream — only the Bayer decode, settled
at read time, is a real, undoable error.

So the read-time override is gone and an **Orientation** control sits
with the export settings instead — as captured, flip vertical/
horizontal, 180 degrees, 90 CW, 90 CCW — applying to preview and
export at no re-run cost, and lossless: verified as a pure pixel
permutation in all six states.

Applied at the final canvas blit, since the composite is computed in
the layers' own coordinate frame — turning it sooner would move the
disc out from under every radial weight. Rationale: FITS has no
orientation keyword, and ROWORDER describes row order, not which way
the camera was held, so a portrait-shot dataset arrives on its side
with nothing to reveal that.

## 0.16.0

A 25-tier 200 mm FITS dataset, run on 0.15.5, produced a log with five
separate faults. Four are fixed here; the fifth is now
user-overridable.

The photometric chain was biased by construction: its factors ran
58.202 down to 0.043, tripping the "tiers disagree" warning. Per-link
residuals, which should each sit at 1.000, showed 23 of 24 below it,
median 0.798 — one error made 24 times (0.798^24 = 0.0044).

The cause is the selection, not the arithmetic: the link was
`median(b/a)` over pixels passing `a > floor(a) & b > floor(b)`, and
thresholding a tier against its own noise biases it upward inside the
selection — Eddington bias. Measured on synthetic pairs (true ratio
1.000, nine links compounded):

    median(b/a), select on both (shipped)   0.785
    median(b/a), select on b only           1.366
    sum(b)/sum(a), select on b only         1.246
    Huber slope b~a, select on b            0.815
    sum(b)/sum(a), NO data selection        1.0015   <-- now shipped

Only removing the selection is unbiased; the saturation mask stays,
keying on saturation, not noise. Caveat: the black level must be
subtracted first, since a pedestal biases every link equally — 64 ADU
left in reads 0.005 for sums, 0.077 for the old medians, both wrong.
The warning now shows per-link residuals and flags when most lean the
same way.

FITS came out in raw sensor colour, since `FitsFrame` has no white
balance or colour matrix to give. A CFA sensor is far more sensitive in
green than red — the corona measured R/((G+B)/2) = 0.56 against 1.11
for a colour-managed stack, a blue-cyan rim on brown sky, leaving the
prominence gate's reference colour meaningless (0 px flagged). The
balance is now measured once from the inner corona (1.05-1.6 R, after
merge and photometry), since Thomson-scattered K-corona light is
wavelength-independent and carries the Sun's own spectrum. Gains beyond
8x between channels are refused; raw datasets are untouched.

The Windows sky-gradient failure was the wrong bug. It reported

    sky gradient removal skipped ([Errno 22] Invalid argument:
    '.eclipseforgehdr\hdr_rgb.npy')

— the opposite end of the file: `remove_sky_gradient` memory-maps
`hdr_rgb.npy` and saves over it while the mapping is open, which
Windows locks against. The map is now released before the write.
Verified unchanged — 1.064 / 1.049 / 1.025 per channel, as in 0.15.3.
An earlier theory blaming a "OneDrive" path is retired as coincidence;
`load_big`'s fallback would already have caught a real mmap failure,
and its docstring now reflects the evidence.

FITS row order is now overridable by hand — `ROWORDER` is a convention,
not a guarantee — via a toolbar control ("from the header" default,
"bottom-up", "top-down"), part of the cache key. The white-balance
reference now decimates against the disc rather than a fixed factor,
so a 108 px moon gets the same ~15k reference pixels as a 524 px one.

Cached products from 0.15.x are not reused: the photometric factors
changed.

## 0.15.5

0.15.4's report said no light reached the picture, but a prominence still
showed in the composite: the gate found 0 px there, while the mask
covers everything, not just gate signal.

Measured on the 12 px annulus 0.15.4 stopped covering:

    median luminance, 160-200 deg      55253
    median luminance, everywhere else  43950
    excess in that sector              1.26x

That 26% is prominence light the gate doesn't cover. The report now
states what the gate measured, not what the picture shows:

    prominences  : corona colour R/GB = 1.11, gate threshold 1.35-1.88, 60 px flagged
                 : none of it outside the disc mask - the gate found redness only
                   at or inside the limb, so the prominence slider has nothing to
                   act on. Prominences may still be visible as brightness; this
                   layer only adds the ones it can identify by colour

Report text only; 0.15.4's cached products are reused.

## 0.15.4

An imported 16-bit sRGB stack got a teal cast absent from the source;
three causes were measured.

**Colour.** Median chroma, normalised to luminance:

                       source           as shipped        temp/tint 1.0
    limb  1.02-1.15 R  1.058 .990 .925  0.832 1.062 .885  1.046 .996 .901
    inner 1.15-1.5  R  1.088 .985 .889  0.853 1.058 .853  1.074 .991 .868
    mid   1.5 -2.5  R  1.155 .973 .811  0.902 1.050 .791  1.133 .980 .803

Sky-gradient division moves colour ≤0.1%; `bgNeutral` accounts for 1.5%.
The rest is `temp` 0.9 and `tint` 1.205, a by-eye white-balance
correction redundant for an already colour-managed import. Both now
start at 1.0 for imports; neutral reproduces the source within 2%. RAW
output is unchanged.

**Disc mask.** The import path used `max(4.0, 0.042 * R)` — 22.0 px on a
524 px Moon. It now measures the 20-80% brightness transition instead,
since mask width depends on seeing, not disc size: 11.0 px here at the
90th percentile, so 9.9 sufficed where 22.0 was asked.

**Report.** `stats["geometry"]` never carried `Rmask` on import, so it
printed `disc mask at R+0.0 px` when it was actually R+22.0.

**Prominence gate.** Unfixed by the above: 60 px flagged, all between
0.97-1.01 R, under the mask, invisible. Beyond 1.02 R the highest
R/(G+B)/2 to 1.3 R is 1.33 against a 1.35 threshold — no excess left.
The embedded Siril header explains why: 327 frames, background
neutralization, colour calibration, two GHS stretches (pivot 0.001
amount 145.65, pivot 0.079 amount 7.57), two background subtractions,
full-strength SCNR. The report now states this instead of implying 60
invisible prominences:

    prominences  : corona colour R/GB = 1.11, gate threshold 1.35-1.88, 60 px flagged
                 : none of it outside the disc mask - the gate found redness only
                   at or inside the limb, so nothing reaches the picture

**Linearity.** Held despite the two GHS stretches. Log-log falloff over
1.1-2.4 R:

    as stored, treated as linear        -1.52     (reference at display gamma: -1.72)
    after the ICC sRGB inversion        -3.19     (reference scene-linear:     -3.38)

The embedded profile put this file within 6% of reference slope; face
value would be off 2x. Nothing clipped: 0.000% at both ends.

Cached products from 0.15.3 aren't reused: the disc-mask margin feeds
`geometry.json`, which every detail layer builds against.

## 0.15.3

An imported stack failed with "could not find the lunar limb" — it had
actually found a 37 px disc on a 500 px Moon.

`find_disc` uses the log-gradient of the image, since relative contrast
is exposure-independent, but has unbounded variance near zero; the
guarding floor was the 1st percentile of the whole frame. When the disc
is darker than sky and covers over 1% of the picture, that percentile
lands inside the disc and the shadow's own noise becomes the strongest
edge:

| region | value above the floor | median &#124;∇log&#124; |
|---|---:|---:|
| inside the disc | 23 | 0.059 |
| at the limb | 2892 | 0.504 |
| sky | 1447 | 0.011 |

42% of the disc's interior pixels cleared the strong-gradient threshold —
4548 votes inside the shadow against 2092 on the limb — collapsing the
circle fit inward.

The fix weights the gradient by signal as well as contrast:
`(s-lo)/(s-lo + 0.002*span)`, 0.97 at the limb and 0.19 in the shadow.

| | before | after |
|---|---|---|
| the failing stack | no disc found (R=37 px) | centre (1523,2390) R=524 px, 720/720 rays, rms 2.34 px |
| synthetic sweep, R/short 0.030-0.300 | all pass | all pass, unchanged |

## 0.15.2

Windows reports a compiled extension failing to load as `DLL load failed
while importing _rawpy: The specified module could not be found.` —
naming neither cause nor fix. The app now detects which of two known
causes applies: a missing Microsoft Visual C++ Redistributable, or an
ARM64 Python, for which rawpy publishes no wheel — both now in the
README.

## 0.15.1

Every field wanting a location now has a Browse... button: raw folder,
flats folder, single-HDR import, diamond-ring frame. Folders list frame
counts; pickers offer only files of the right type;
Escape or a click outside closes it; picking the raw folder loads it
immediately — typed paths had cost real testing time, including a valid
path that left Start permanently disabled.

A native `<input type="file">` can't work: it hands JavaScript only the
file's content, never its location, and a server-side dialog needs Tk to
own the main thread on macOS while this server is threaded. The server
instead lists directories and the page draws the picker.

## 0.15.0

The first external run on FITS from a 249-frame bracket surfaced four
defects, all now fixed, plus two diagnosed but not yet fixed.

**FITS came out upside down.** FITS row 1 is the bottom of the frame,
while every other format here has row 0 at the top, so a bracket read
straight through comes out mirrored. It now reads `ROWORDER` when
present, else follows the standard. Flipping an even row count also
swaps Bayer row parity — RGGB becomes GBRG — so the parity shift now
travels with the data. Tested at even/odd heights with `ROWORDER`
absent, `BOTTOM-UP` and `TOP-DOWN`: correct in all four cases.

**One OneDrive folder silently disabled three features.**

    sky gradient removal skipped ([Errno 22] Invalid argument:
    'hdr_rgb.npy')

`np.load(mmap_mode="r")` fails on a cloud-synced Windows folder, used in
three places — sky-gradient fit, contact-frame loader, colour estimate —
so one folder took out three features, only one reported. It now falls
back to a plain read.

**The diamond ring did nothing.** The contact-frame path field had no
Enter handler, unlike the folder field above it. Enter now loads it.

**Importing one HDR needed a raw folder that an import does not use.**
Start was enabled only by a successful folder load, so a pasted TIFF
path left the button greyed out. Typing an import path now enables
Start.

**Diagnosed, not yet fixed.** The merged limb fit came out R = 542 px
against a per-tier consensus of 526 — 3% large, below the 15%
tier-consensus check — landing the mask 41 px beyond the real Moon and
eating the prominences; shrinking the mask reveals the merged limb's own
26 px transition. Separately, photometric calibration ran away — factors
9.40, 7.08, 5.36, 4.07 on the four shortest tiers, a smooth geometric
progression rather than noise — flagged by the pipeline as unreliable,
as intended. Cause unknown.

## 0.14.6

A step's cost is the gap to the next log line, so a run's final step was
never measured — Pellett was quietly absorbed into "steps under 1s". The
gap is now closed before it is measured.

The same report confirmed the 0.14.4 weighting: with every step over a
second listed, the detail stage on a 50-frame 45 Mpx run comes to 343 s
of 734 s — 47% — against the 48% the bar allocates. The earlier estimate
could only bracket it between 42% and 64%.

## 0.14.5

NAFE showed a thick black ring hugging the disc, with a bright rim at
the limb, where FNRGF was clean.

Cause: the envelope-flattening kernel added in 0.14.0. A symmetric
Gaussian mean is a poor background estimate beside a large dark hole:
pulled down near the limb by the occulted disc, so `L - mean` overshoots
bright at the limb and undershoots dark outside it, ~100 px of black at
0.08 R on a 622 px disc. 0.14.0's RMS-based check missed it, since a
broad depression barely moves an RMS; it's obvious in mean layer value
per radial ring, now the measure used.

The local mean is now built by normalized convolution over non-disc
pixels only, as MGN already does. Measured on a synthetic corona with
known fine modulation:

| flattening | ring depth | correlation with the true structure |
|---|---:|---:|
| plain Gaussian (0.14.0-0.14.4) | 0.1684 | 0.468 |
| normalized convolution | 0.0713 | 0.814 |

Mean layer value at 1.15 R went from 0.333 to 0.445 against far-field
0.50 — ring depth more than halved — and streamers now run cleanly to
the limb.

**Why not Fourier subtraction.** It scores better on a perfect limb fit
(0.0047/0.938) but needs a centre/radius and fails hard when wrong:

| centre error | ring depth | correlation |
|---|---:|---:|
| exact | 0.0047 | 0.938 |
| 0.02 R | 0.1347 | 0.666 |
| 0.05 R | 0.2267 | 0.542 |
| 0.10 R | 0.2533 | 0.416 |

Normalized convolution degrades gently instead: 0.10 R off centre it
still scores 0.1135/0.731, beating the plain Gaussian's best case
(0.1684/0.468). σ stays at 0.08 R (0.15 R scored 0.0576/0.863 on the
sweep). Re-running is required.

## 0.14.4

**Progress split, measured instead of assumed.** 0.14.1 gave the detail
stage 22% of the bar on an assumption — stacking costs about 3.5x the
detail stage. Measured instead, on a 50-frame 14-tier bracket from a
45 Mpx body, 12m27s end to end:

| step | time | share of run |
|---|---:|---:|
| inner corona: raw pass | 141 s | 19% |
| inner corona: denoised pass | 140 s | 19% |
| MGN | 31 s | 4% |
| three named detail steps | 312 s | 42% |

Remaining steps bracket the stage between 42% and 64% of the run —
roughly equal to stacking, not 3.5:1 — and the bar was reaching 78% at
halfway then crawling. Detail now starts at 52%, the low end.

**The report's shot span counted rejected frames.** Time span and ISO
list came from every file in the folder, before rejection — a test
frame shot three weeks after the eclipse and correctly dropped still set
the span:

    frames taken : 2026:08:12 21:30:12  ..  2026:09:01 13:35:30

Both are now computed from frames actually stacked: dropping the stray
frame makes the span match the clean run exactly (10:20:00-10:23:14 vs
10:20:00-10:59:00); a genuinely used stray frame still counts.

## 0.14.3

Found by dropping a blown file into a healthy folder: the 0.14.0 guard's
problem, in three positions. Two were quietly wrong.

**A saturated shortest tier disabled hot-pixel repair.** Hot pixels map
once on the shortest (darkest) tier, where a saturated frame has no
"above the neighbours" to detect. A stray overexposed file at the short
end becomes the shortest tier alone, and its map found zero defects on a
sensor carrying 434. A tier more than half clipped is now skipped, next
tier used, logged.

| defect map built from | before | after |
|---|---:|---:|
| healthy shortest tier (control) | 434 | 434 |
| stray blown frame as the shortest tier | 0 | 434 |

**The limb-radius warning claimed an override it had not made.** Past
15% disagreement, the log said:

> ... the individual tiers agree on R=80px (+22%). Using the tiers' value —
> the merge probably contains frames of different scenes.

The code only adopts the tiers' value past 30%; 15-30% kept the merged
fit but printed this anyway. Each branch now states what it did:

> ... KEEPING the merged fit: it's a per-azimuth measurement and the
> disagreement is under 30%. Check the disc mask on the preview — if
> it's the wrong size, the tiers were right.

Cached 0.14.2-and-earlier products aren't reused: the defect map and
limb radius both feed the merge.

## 0.14.2

0.14.1 restored the button but drew the wrong layer into it. The server
served the flat correctly, but the browser's preview picked its layer
source with a ternary chain ending in a bare fallback to the prominence
gate:

    VIEW==="pellett" ? L.pel[j]/255 : pr[j]/255

`prom` "worked" only by being the fallback — any later view, Flat
included, silently drew it instead. Mean absolute difference between
what each button draws:

| | before | after |
|---|---:|---:|
| Flat vs Prom gate | 0.00 (identical) | 160.69 |
| Flat vs what `/api/layer/flat` serves | — | r = 1.0000 |

Each view now names its own source, falling back to the composite with a
console warning if its layer is missing. A new test (`viewtest.py`) now
drives the real page in a real browser to catch this class of bug.

## 0.14.1

**Master-flat preview failed to load.** 0.14.0 added the master flat as
a QC preview; its button never appeared on the first real run. The
loader reduced the flat to superpixels, then cropped that half-res array
to the full-res layer grid — asking for 5358x8100 out of a 2716x4076
array, which can never match, so the loader returned `False` silently.
Fixed by cropping on the superpixel grid first, then expanding. Verified
on the reporting shapes and three others:

| case | before | after |
|---|---|---|
| 5432x8152 sensor, 5358x8100 grid, crop (32,34) | no button | loads, dust visible, no checkerboard |
| no autocrop | loads | loads |
| odd layer grid | loads | loads |
| old cache, no `crop_origin` | loads (centred) | loads (centred) |

The reason now goes to the log and to `/api/geometry` instead of failing
silently.

**Progress reporting reweighted by measured time.** The bar sat nearly
full for a long time with no other sign a run was alive. The detail
stage was timed step by step at two sizes:

| step | 10.8 Mpx | 43.4 Mpx |
|---|---:|---:|
| denoise HDR master | 2.9 s | 18.2 s |
| MGN | 26.9 s | 197.9 s |
| FNRGF | 7.6 s | 22.9 s |
| NAFE | 11.1 s | 56.3 s |
| inner corona | 135.3 s | 1001.1 s |
| prominence colour | 0.9 s | 3.6 s |
| Pellett | 11.1 s | 31.7 s |
| total | 198.7 s | 1343.9 s |

The inner-corona block, two-thirds to three-quarters of the stage,
doesn't scale with pixel count: 4x the pixels cost 6.8x the time, so at
45 Mpx the detail layers alone take 22 minutes — the bar gave that 6.5%
of its width. It now gets 22% instead, divided by measured wall time,
with the inner-corona block's two long passes at 7.4% each, up from
about 1%. The 22% remains an estimate, assumed at 3.5x stacking. Beside
the bar there is now a spinner, an elapsed clock, and the current
step's name; import lowers the detail band to 12%. Every finished run
now prints its own timing summary to the log and `report.txt`.

**Cache compatibility.** 0.14.1 changes nothing on disk, so a 0.14.0
folder remains reusable, via an explicit `CACHE_COMPAT` declaration.

## 0.14.0

**Alignment failure on a fully saturated tier.** `SVD did not converge
in Linear Least Squares` partway through a run: `prep_pair` masks
saturated pixels before correlating tiers, and a fully saturated window
prepares to exactly zero, so phase correlation returns `err` NaN,
poisoning a row of the design matrix:

| input | returned shift | returned `err` | weight |
|---|---|---|---|
| two real images | [0, 0] | 0.0009 | 19.6 |
| one all-zero | [-0.75, -0.75] | nan | NaN |

Reproduced with a synthetic bracket spanning 14.3 EV, longest tier
blown. Degenerate windows are now caught before correlating, and an
unlinkable tier takes its nearest linked neighbour's shift instead of a
meaningless answer.

**Import of a finished HDR.** Pointing a new box at one 16-bit TIFF or
FITS runs the disc fit, sky-gradient fit and every detail layer,
exporting exactly as a stack does. The tone curve is read from the file,
not guessed: a Photoshop mean stack of this app's own exports measured
corona log-log slope -1.72 against -3.38 scene-linear (ratio 0.510 vs
~0.45 for sRGB); inverting the embedded ICC profile brought it to -3.17,
within
6% of truth. No-profile files must declare themselves; 8-bit is refused.
The report states what an import gives up: alignment, photometry,
per-tier lunar masking, earthshine, an independent inner-corona
measurement, and prominence-gate strength (R/GB 1.84 vs 3.02 from a real
fast tier).

**NAFE input flattening.** NAFE alone still carried the full radial
falloff before processing: 60% of its output range went to the gradient
instead of structure, driving the equalisation into a dark ring near the
limb. It now comes off with a plain Gaussian high-pass at 0.08 R,
needing no circle or limb fit:

| | large-scale range | detail 1.05-1.5 R | ripple outside the mask |
|---|---|---|---|
| as shipped | 68.3% | 0.0611 | 0.3451 |
| fitted radial profile (needs the limb) | 5.4% | 0.1328 | 0.2543 |
| Gaussian 0.08 R (no geometry) | 3.6% | 0.1360 | 0.1210 |

2.2x the near-limb detail and 2.9x less ripple — the geometry-free
option won on every count. A smooth model subtracted from the finished
layer instead flattens it (68% -> 9%) but recovers zero detail: the
nonlinear equalisation had already done the damage.

**Resolution-adaptive MGN scales.** The top scale is now tied to R
(`0.0643 * 622 = 40 px`), the bottom raised to whatever the image
actually resolves. At 3.19 arcsec/px behind a 240 mm zoom, the old
1.25/2.5 px scales sat below the optics: 9.6% of the local mean against
1.9-2.7% where real structure lives.

**Other changes.** Master flat preview gets its own view button, cut to
the layer grid and stretched to its own 0.5-99.5 percentile. The flat
now repairs its own outer border (masked photosites: 0.644 vs 0.936 four
rows in). Exported detail views carry their own p1/p99 range in the
TIFF description — MGN occupies 14% of the 16-bit range where FNRGF
occupies 81%. "Cached pipeline products found" now checks build version.

## 0.13.1

Fixes two failures from one wrong assumption: that Moon location and
size can be guessed from brightness and frame size.

**Centre.** `find_center` took the centroid of the brightest 0.05% of
pixels — valid only while inner corona rings the disc evenly. When one
sector dominates, those pixels sit on one arc and the centroid lands on
the limb: 646 px out on a 620 px disc, 1.04 R, on the reference set.

| seed centre offset | result |
|---|---|
| 0 – 500 px (0.81 R) | R 622.3, rms 5.00, 720/720 rays, 0.2 px from truth |
| 620 px (1.00 R) | fails |

A 3% change in the limb ring from the new flat-field correction pushed
the estimator past its cliff: a 67 px limb ramp (was 21), `rim nan px`,
"coronal range 1.2 EV" (was 6.5), and a sky gradient refitted 1000 px
from the Moon.

**Radius.** `fit_limb` searched between a tenth and a third of the short
side — a statement about focal length, not eclipses. A 240 mm frame puts
the disc at 7.5%, below the floor: it returned R = 1005 px for a 301 px
disc.

**Replacement.** A new `find_disc()` locates the disc from the limb
edge — a 93x relative step over ~20 px that dominates the log-intensity
gradient regardless of exposure. A wide bright mask seeds a centre, the
radius comes from the peak of the azimuthal log gradient, and
strong-gradient pixels near it are fitted with a circle, alternating.
Two gates confirm a real limb: inliers span 270+ degrees, and scatter is
tight (0.04 R for a real limb; 0.16-0.23 R for noise).

Measured across disc sizes 3-30% of the short side, three aspect
ratios, lopsided up to 4:1, prominences up to 3x limb brightness —
18 + 7 cases:

| | before | after |
|---|---|---|
| `fit_limb` R at 7.5% of the short side (240 mm) | +344% | -1.2% |
| `fit_limb` R at 3% of the short side | +1232% | -5.2% |
| centre, reference merged HDR | 646 px (1.04 R) | 62 px |
| centre, short-exposure stack | 611 px | 1 px |
| centre, lopsided inner corona | up to 1.03 R — past the cliff | 0–15 px |
| `fit_limb_rays` from that seed | failed | 720/720 rays, 0.2 px from truth |
| seed cost | — | 0.08 s/frame, 4 s per 49-frame bracket |

The measured radius now seeds everywhere that used to guess from a frame
fraction.

## 0.13.0

Flat frames in a `flats/` subfolder are found and used automatically; a
path box takes flats from elsewhere, or `off` ignores a present folder.

A master flat is built once per folder, cached, and divided out of every
frame first — removing vignetting, dust shadows and per-photosite
sensitivity. This matters more for an eclipse: the corona's own radial
falloff is the signal, so a 6% vignette is a 6% error in it.

**The master flat measures its own noise, which sets how much it is
smoothed.** Dividing by a flat injects its noise into every frame
identically, so stacking can't average it away. Frames split into two
half-stacks; since `var(mean of half) = 2 x var(mean of all)`, their
difference is the master's noise directly. Gaussian σ is raised until
measured noise falls under 0.2% per photosite:

```
flat: min/max-trimmed mean of 20 frames from flats/
flat: per-pixel noise 1.318% -> 0.183% after a 5.7 px smooth (target 0.2%)
flat: corrects a 8.4% falloff — the dimmest part of the field sits at 0.923
      of the brightest
```

A clean flat set keeps its dust motes at full resolution; a thin one
degrades gracefully to a vignetting model instead of adding noise. The
white-noise starting guess for σ needed 2.9x more smoothing in practice,
since sensor noise isn't white.

Other results: per-Bayer-channel normalisation to frame centre cancels
drift between flats (central level moves 0.4%, matching the baked-in
vignette to 0.9% median over 60% of the field); min/max-trimmed
per-pixel combine (5+ frames) drops a cosmic ray or satellite instead of
averaging it in (15% contamination becomes 1.4%); the clipping test now
runs before the flat is divided out (a 30% vignette would otherwise give
a corner pixel deserving weight 1.000 only 0.269); unusable flats say so
(over 85% saturation, under 2%, or wrong frame size) and the run
continues uncorrected; flat files are part of the cache key.

This does not replace the sky-gradient fit from 0.11.x: vignetting is
radial and fixed by the optics, the sky gradient a tilted plane fixed by
the atmosphere; a radial model explains 0.0% of it on the reference set.

## 0.12.0

Folders of FITS frames now work as input, for capture software that
writes it instead of camera raw — INDI/EKOS, SharpCap, N.I.N.A.,
FireCapture. Colour (CFA + `BAYERPAT`), monochrome, and debayered
3-plane cubes are supported; Bayer pattern rolls to RGGB as with camera
raw, and `XBAYROFF`/`YBAYROFF` are respected. No new dependency:
`astropy.io.fits` is used if installed, then `fitsio`, then a built-in
reader for plain uncompressed FITS; astropy is an optional extra for
tile-compressed files.

A FITS frame arrives without two things LibRaw normally supplies, both
recovered from the file: exposure, from `EXPTIME`/`EXPOSURE` (required;
a file without it stops the run, since it groups frames into tiers), and
the saturation ceiling, from `SATURATE`/`DATAMAX` if present, else a
saturation plateau under 20% of the frame, else bit depth — the report
states which was used.

There is no colour matrix or white balance in a FITS header, so both are
identity — a colour-camera frame comes out green-dominant as captured,
with Warmth/Tint/Neutralise sky cast as the controls.

New defaults:
radial flatten 0.5, FNRGF share 0.17 / strength 0.9, MGN contrast 0.04,
NAFE-VN mix 0.15, Pellett off, short-exposure detail 0.31, glare dim
0.05, warmth 0.9, tint 1.205, neutralise sky cast 1.0, saturation 1.0,
highlight compression 0.1, output gamma 1.0, black point 0.005.

## 0.11.5

Setting every structure slider to zero and still seeing the artifact
traced it to the base envelope, not any detail layer.

`Bg` normalised merged luminance with `lo = percentile(lum, 2)`. On a
wide field the sky is most of the frame, so the black point lands within
a noise sigma of the sky, and `xn` becomes a small difference between
nearly equal numbers, turning tiny real variation into an enormous one:

| | corner | mid-edge | ratio | sky clipped to black |
|---|---|---|---|---|
| the data itself | | | 1.030 | |
| `lo = p2` (before) | 0.0478 | 0.0637 | 1.332 | 2.4% |
| `lo = sky - 5σ` (now) | 0.1253 | 0.1284 | 1.025 | 0.0% |

A 3% brightness difference was displayed as 33%, with 2.4% of the sky
crushed to black. The black point is now measured from the sky itself
(median and MAD beyond 2.5 R), clamped never above the 1st percentile.

It also un-pinned radial flatten: `rprof` had sat on its 0.12 clamp
everywhere outward, doing nothing at any setting. It now reads
0.165-0.170 and works.

End to end, composite with every detail layer at zero: corner-to-edge
ratio 1.332 -> 1.000; with default layers on, 1.017. The background now
starts lighter:

| bgBlack | sky | corona | corona/sky |
|---|---|---|---|
| 0.02 (default) | 0.185 | 0.251 | 1.36 |
| 0.08 | 0.131 | 0.202 | 1.54 |
| 0.11 | 0.102 | 0.179 | 1.75 |

About 0.11 reproduces the old background darkness — but corona-to-sky
contrast now improves as you crush, which it could not before.

## 0.11.4

**Neutralise sky cast had no effect at any setting, and never had.** It
divides chroma by the measured background colour, but that colour was
measured on `ratio`, whose confidence fade drives it to exactly 1.0
where signal is near the noise floor — precisely where it was measured
(mean confidence 0.015). It returned R 1.000 G 1.000 B 1.000 on a sky
whose real colour is R 0.985 G 1.037 B 0.681; it is now measured on the
HDR itself. That exposed dividing the whole field by the real value
tipping the far sky blue, since `ratio` had already forced that region
neutral — now weighted by the same confidence that built `ratio`:

| bgNeutral | far sky B/R | corona B/R |
|---|---|---|
| 0.00 | 1.000 | 0.166 |
| 0.50 | 1.000 | 0.200 |
| 1.00 | 1.000 | 0.239 |

**The sky gradient is now removed per channel**, taking colour with it.
Fitted spans: R 1.202x G 1.262x B 1.335x, blue steepest:

| | before | after |
|---|---|---|
| brightness spread across quadrants | 1.122x | 1.019x |
| colour spread, R | 0.0551 | 0.0106 |
| colour spread, G | 0.0123 | 0.0024 |
| colour spread, B | 0.0429 | 0.0072 |

**New defaults**, the settings the reference bracket was worked to by
eye: MGN contrast 0.4 -> 0.1, FNRGF strength 1.5 -> 0.5, NAFE-VN mix 0 ->
0.2, short-exposure detail 0.24 -> 0.6, glare dim 0.35 -> 0.15, radial
flatten 0.2 -> 0.35, base lift 0.18 -> 0.255, grain smoothing 0 -> 0.25,
neutralise sky cast 0.7 -> 0.5, and a few smaller moves.

## 0.11.3

0.11.2 removed the tilt but left a curved remainder, measured beyond 4.1 R:

| model fitted to the remainder | explains |
|---|---|
| another plane | 0.8% ← the tilt really was gone |
| radial about the frame centre (vignetting) | 2.5% ← still not vignetting |
| **full quadratic in x, y** | **51.8%** |

The sky curves near the horizon; a plane cannot follow it. Model is now a full
quadratic, fitted beyond the measured corona extent:

| shell | before | after | retained |
|---|---|---|---|
| 1.6–2.4 R | 0.1735 | 0.1608 | 93% ← corona survives |
| 2.4–3.2 R | 0.0810 | 0.0576 | 71% ← corona survives |
| 4.0–6.0 R | 0.0594 | 0.0088 | **15%** ← sky collapses |

Falls back to a plane below 100k fitted pixels; a fit spanning more than 2× the
frame is rejected. Full-resolution evaluation now runs in row blocks.

## 0.11.2

A broad dark gradient crosses one side of the frame, measured at 1.20× corner to
corner. Not vignetting (radial model explains 0.0% vs 54% for a plane); flats do
not remove it, so it is fitted per run.

Fitted close in, the model absorbs real corona asymmetry instead: amplitude 1.22
(1.6–2.4 R) to 0.22 (beyond 4 R), direction stable −13° to −18° — one real
gradient. Fit restricted to beyond the measured corona extent:

| shell | before | after | retained |
|---|---|---|---|
| 1.6–2.4 R | 1.2235 | 1.0381 | 85% ← corona, survives |
| 2.4–3.2 R | 0.5530 | 0.3676 | 66% ← corona, survives |
| 3.2–4.0 R | 0.3425 | 0.1572 | 46% |
| 4.0–6.0 R | 0.2184 | 0.0331 | **15%** ← sky, collapses |

Multiplicative and identical across channels: colour untouched. Below 2%
amplitude or 8σ, left uncorrected.

Must fit in true log, not `log1p(S/median)`, which compresses where the sky
sits: in log1p space the correction came out at half strength (sky residual
0.105 → 0.063 instead of → 0.017).

## 0.11.1

E is a rank: the quiet corona sits in a narrow flat band, and the near-limb
ridge, a local max, pins at 1.0, blowing out the rim.

Layer now rescaled by its own robust spread, rolled off past 3σ. Measured at
full resolution (rim = 1.04–1.18 R):

| | 1.05–1.5 R | 1.5–2 R | 2–2.5 R | 2.5–3 R | 3–4 R | rim p99.9 | rim >0.99 |
|---|---|---|---|---|---|---|---|
| 0.11.0 | 0.0747 | 0.0272 | 0.0268 | 0.0276 | 0.0301 | 1.0000 | 6.75% |
| **0.11.1** | **0.1225** | **0.0882** | **0.0843** | **0.0856** | 0.0940 | **0.8792** | **0.00%** |

Corona contrast improves 3.1–3.3×, limb 1.6×, rim stops clipping. Repeatability
(independent σ_A) unchanged (+0.96 limb, +0.77–0.81 outside): a rescale, not new
structure. Rolloff knee of 3 chosen by measurement: at 2 the corona loses
contrast, past 6 the rim clips again.

Sky looks grainier; the mix slider and composite envelope control how much
reaches the final image.

## 0.11.0

Variable neighbourhood was correct (eqs. 11–12, NAFE-VN); the bug was applying ε
and σ in rank units, not A's value units (eq. 13 needs σ fixed there) — under a
rank map the sky filled most of the axis, so the same width meant something
different at every radius.

Per-annulus contrast (sd), full resolution:

| | 1.05–1.5 R | 1.5–2 R | 2–2.5 R | 2.5–3 R | 3–4 R (sky) |
|---|---|---|---|---|---|
| 0.10.4, rank units | 0.0944 | 0.0105 | 0.0253 | 0.0602 | 0.1041 |
| **0.11.0, level units** | 0.0747 | **0.0272** | **0.0268** | 0.0276 | **0.0301** |

In rank units contrast rose with radius, tracking falling SNR — a noise
detector. In level units it's nearly constant: corona 1.5–3 R gains 40–160%, sky
loses two thirds of its grain. Repeatability (independent σ_A): +0.99 limb,
+0.77–0.82 elsewhere.

- **K = 128, not 64.** 64 levels cost contrast (2–2.5 R: 0.0273 vs 0.0366 at
  K=128); past 128 the sky gains faster than the corona.
- Paper's plain Gaussian kernel beats the multiscale sum used for fuzzy weights:
  more corona contrast, nearly 2× faster.
- Speed is legitimate: evaluating every pixel's histogram as K blurred maps is
  the same computation reorganised. Grid 8 vs 4 agrees to three decimals, so it
  is free; K=64 was the shortcut that was not.

## 0.10.4

0.10.3 broke the NAFE layer; reverted to 0.10.2 behaviour.

Weighting the rank map toward r < 1.5 R starved the mid corona: 1.5–3 R
collapsed to sd 0.044 (p1–p99 0.46–0.75), a flat plateau with a hard edge;
inside 1.5 R it clipped at 1.000.

Shipped because the check metric, high-pass sd, misses mid-scale collapse that
per-annulus contrast catches; the check image was cropped to 2.4 R and both
checks ran on a decimated probe.

Verified after revert: 0.10.4 matches 0.10.1 to within 0.002 sd in every annulus
from 1.05 R to 4 R. The 0.10.2 changes (ε, noise σ, working sigma_sp) are kept.

## 0.10.3

Testing the crop-the-field hypothesis from 0.10.2: the rank map spends its
resolution on sky — a 600 mm bracket is 6.5 × 4.4 lunar diameters,
overwhelmingly sky.

Cropping to 2 lunar diameters tripled corona detail but throws away real corona;
a hard mask was worse (61% pinned at rank 0 or 1). Fix: a weighted rank map —
full weight inside 1.5 R, 0.10 outside.

"repeat": correlation between two runs differing by one independent σ_A of added
noise:

| | inner detail (1.05–1.8 R) | repeat | outer detail (2.2–3.2 R) | repeat |
|---|---|---|---|---|
| 0.10.1 | 0.00559 | +0.906 | 0.04850 | +0.663 |
| 0.10.3 | **0.01552** | **+0.961** | 0.04246 | +0.656 |

Inner-corona detail improves 2.8× and is more reproducible; outer corona pays
12%. σ_A must also be measured over the whole frame — inside 1.5 R it made
differences signal-dominated, over-fired the level smoothing, and cost the outer
corona 38% of its detail. The FNRGF-correlation metric used in 0.10.2 is
contaminated by shared noise; the added-noise repeatability test replaces it.

## 0.10.2

Two of three NAFE constants had never been swept (corona = 1.05–2.2 R, sky =
beyond 3.2 R):

- **ε 0.05 → 0.10.** At 0.02 it prints concentric rings; at 0.05 high-pass
  structure is 0.076; at 0.10 it is 0.109, FNRGF agreement unchanged (0.711 →
  0.721); past 0.15 the extra is noise.
- **Noise σ: 2 σ_A → 4 σ_A**, within the paper's 2–12 range:

  | σ / σ_A | corona detail | sky grain | detail/grain |
  |---|---|---|---|
  | 2 | 0.00843 | 0.12421 | 0.068 |
  | **4** | **0.00820** | **0.06621** | **0.124** |
  | 6 | 0.00782 | 0.04488 | 0.174 |
  | 12 | 0.00696 | 0.03085 | 0.226 |

  Sky grain halves from 2 to 4 for 3% of the corona detail; past 4 a halo
  swallows the streamers by 12.
- **`sigma_sp` now does something** (was ignored on the fuzzy path): the widest
  multiscale kernel scale, default 0.13 R; flat from 0.06 R to 0.51 R, so
  scaling with the disc matters, not the value.

Detail-to-sky-grain improves 0.094 → 0.124 — real but modest.

## 0.10.1

The NAFE layer was 99% a gamma stretch: the paper's B = (1−w)·T_γ(A) + w·E (eq.
2, w in 0.05–0.3) was stored as the detail layer — at w = 0.2, four fifths gamma
transform, diluting MGN and FNRGF.

Layer is now E, the equalised field. Measured (×4 decimated, K=64, γ=2.4,
ε=0.05):

| | stored before (B) | stored now (E) |
|---|---|---|
| correlation with a plain gamma stretch | 0.992 | 0.562 |
| high-pass structure (sd) | 0.0180 | 0.0675 |

Eq. 2 now happens one level up: the composite's envelope is T_γ and nafeMix is
w. `nafe_vn(..., combine=True)` still returns the paper's B standalone.

## 0.10.0

A line-by-line review of every module. Nothing here changes the reference render
except where noted.

**Wrong pixels**

- Clipping test compared white-balanced max against raw saturation level: the
  ~2.1× red gain pushed prominence H-alpha past threshold a stop early, filling
  cores by leakage.
- The long-exposure luminance and prominence-colour layers were both saved
  uncropped: the former crashed earthshine runs with a non-zero alignment trim
  at 99%, the latter's ×2-upsampled gate painted prominences off the limb by the
  crop offset.
- Prominence geometry cached after layers were built: a first run used the
  merged limb; only a rebuild used the prominence tier's own.
- Autocrop skipped when only bottom/right needed trimming: tier TIFFs were
  sliced, the render was not, leaving edge-replicated pixels for MGN.

**Preview vs export** (JS had drifted from the numpy it mirrors, now resolved):

- Disc-trim slider moved the mask 4× too far in preview.
- FNRGF ramp and inner-corona weight both used mismatched curves: FNRGF linear
  1.05–1.40 R vs smootherstep 1.02–1.57 R (0.43 vs 0.25 of the mix at 1.2 R);
  inner-corona linear 1.38 R/0.33 R vs smoothstep 1.45 R/0.40 R.
- Prominence lift's "presence" bias existed only on export, so preview showed
  none at neutral; chroma ratio was clipped to [0, 2.5] in preview vs [0.2, 3.0]
  in export, understating saturation.
- NAFE median sampled at stride 7 over RGBA hit the alpha byte one sample in
  four (~67th percentile, not median), darkening NAFE in preview.

**Colour management**

- Composite export carried no ICC profile; PixInsight decoded it wrong.
  16-bit/8-bit TIFF and JPEG exports now embed sRGB.
- Black level now subtracted per channel, not as a four-channel mean; the
  leftover offset had landed on the outer corona after white-balance gain.

**Constants that encoded the reference dataset** (now fractions of the radius,
not absolute pixels):

- Prominence patch/search radius (60/80 px, 0.19/0.26 R reference) was 0.65–1.45
  R on a 300 px disc, registering tiers to the lunar edge.
- Prominence gate window (0.19 R reference) was 0.4 R on a 300 px disc, letting
  promGain brighten inner-corona loops as prominences.
- Near-limb saturation feather (20 px, 3.2% of reference radius) was a tenth of
  a 200 px disc's radius.
- Short-exposure inner stack took the first four tiers by count — 4.3 EV
  reference, 6.4 EV on a 6-tier bracket, where the 4th tier would carry ~8× the
  1st's weight. Now capped at 24× the shortest exposure.

**Crashes and silent failures**

- Negative radicand in the limb fit produced a NaN radius, crashing at 91%
  before the plausibility check meant to catch it.
- Alignment-network residual reduced over an empty array when every frame shared
  one shutter speed and only measured y: 15 px inconsistent in x reported 0.
- Saturation level now cross-checked, since LibRaw can report white level 0, an
  already-black-subtracted value, or a 12-bit ceiling for 14-bit data — any
  silently produced a cached black frame. Frames of differing size in one tier
  now raise a clear message, not a bare broadcast error.
- Prominence-detection sigma computed over zero-filled out-of-frame azimuths
  dragged the threshold to zero near a frame edge, flagging every azimuth as a
  prominence. Now measured over in-frame azimuths only.
- Variance diagnostic used a cropped-frame centre for uncropped tiers, so drift
  runs reported mid-corona numbers as "limb variance."
- Intra-tier motion warning compared full-res displacements against a half-res
  window, firing at half the intended motion.

**Cache and session**

- Cache key was render options only, so replacing a frame without forcing reused
  the old stack; file list, sizes and mtimes are now part of the key.
- A run that died partway left a valid options file a later non-forced Start
  would accept; now cleared before a run, not after.
- A failed run still reported ready and exported old layers under new settings,
  since the previous result stayed exportable through a re-run.
- Run thread re-read the selected folder at execution time; now bound at start,
  switching mid-run refused.

**Known, not changed** (each would move the reference render): MGN's global
weight is 0.12 where Morgan & Druckmüller Eq. 5 uses h = 0.7; MGN's scales and
FNRGF's smoothing lengths are absolute pixels; NAFE sigma_sp is a no-op in the
default fuzzy path; radial flattening before MGN/FNRGF starts at R+3, not the
disc mask radius.

## 0.9.10

- **Glare dim no longer prints a ring.** It shared radial weight with the short-
  exposure layer's smoothstep, closing at 1.45 R, reading as a circle. Glare now
  has its own exponential profile (scale length 0.6 R), no boundary. Outside the
  limb (r > 1.1 R): peak curvature drops from 21.5 at r = 1.42 R to 1.54, peak
  slope from 2.61 to 0.97.

## 0.9.9

- Tier TIFF export now embeds an ICC profile — sRGB by default, or scene-linear
  via a PixInsight checkbox. Untagged linear data was decoded as sRGB, clipping
  histograms.
- Fixed the export normaliser: it divided by sensor saturation after white
  balance, after the colour matrix pushed red/blue above it, clipping
  highlights. Now divides by measured headroom.
- A README now ships with the tier TIFFs: don't mean-stack the set as one
  exposure.

## 0.9.8

- **Fixed the alignment regression from 0.9.5.** The brightest-pixel centroid
  finder, the per-tier disc locator for pre-centring, drifts 197 px with
  exposure on tiers 22 px apart. Coarse centre now comes from the limb fit,
  gated on measured spread (>15% of the correlation window). Measured: fixed
  window 0.58 px, limb pre-centred 1.86 px, centroid pre-centred 20.01 px.

## 0.9.5

- Autocrop of the alignment border, from per-tier edge vacancies and cross-tier
  shifts, applied to render and tier TIFFs; even origin preserves Bayer phase,
  60% guard against over-trimming.
- Per-tier pre-centring of the correlation window on its own disc position, for
  datasets where the disc moves hundreds of pixels between brackets, later gated
  by 0.9.8.

## 0.9.4

- Renamed from CoronaForge to EclipseForgeHDR (branding collision). Version
  added to the build filename, an efhdr alias added, version now dynamic from
  one source (had drifted: 0.7.1 vs 0.8.0).
- Export of aligned exposure tiers as 16-bit TIFFs, for stacking or blending in
  PixInsight, Photoshop or Affinity.

## 0.9.3

- Saturation-based rejection of non-totality frames (saturated area >3× tier
  median dropped), after a Sony dataset's partial-phase crescents pushed the
  fitted lunar radius from 451 px to 857.
- Cross-tier radius consensus overrides the merged limb fit when it deviates
  >15%.
- Disc mask margin is derived from the measured limb ramp, not a fixed number.
- Photometric calibration regions now selected by signal-to-noise against each
  tier's sky noise, replacing a "brightest 20%" rule that rejected every link on
  a short bracket.

## 0.9.2

- **MGN halo fixed.** Two wrong hypotheses first (moving the background start
  radius: 0.0331 → 0.0615, worse; reflecting the layer: +7%). Cause was the
  deband step's azimuthal-mean subtraction leaving the rim's variation untouched
  (residual 0.0140 before and after). Now debanded at order 6, chosen from a
  measured order-vs-rim/streamer trade.

## 0.9.1

- NAFE-VN added as a mixable detail layer, after three failed attempts — 360
  full-res blurs, a flat output, contour rings — fixed respectively by a coarse
  histogram grid, equal-population rank binning, and Gaussian membership with a
  per-level noise sigma.
- Chroma reconstruction gained a noise-relative floor and confidence fade, so
  colour is not invented where there is no signal.
- Radial profile rewritten by sorting rather than binning: 21× faster, no 6000
  px cap.

## 0.9.0

- Prominence anchors: located on a reference tier, matched between tiers by
  normalised cross-correlation, adding hard links to the alignment network. Two
  rewrites — a "how far out does this azimuth stay bright" detector found
  coronal streamers instead.
- Alignment quality metrics added to the report: per-tier limb/corona variance,
  the disagreement rim, and the limb 20–80% transition width.
- Value-based neighbourhood replacing geometric masking in the equalisation.

## 0.8.8

**The cross-tier shift sign was inverted.** All six consumers negated the cross-
correlation shift, leaving a residual of twice the true offset — unnoticed since
all six agreed in the same wrong frame. Measured after the fix: moon-track
scatter 23.5 → 2.0 px, tier-to-tier limb variance 0.377 → 0.073, corona variance
0.165 → 0.136, merged limb 20–80% transition 23.0 → 15.5 px. Explains the "125
px lunar spread," the 0.89 px/s apparent drift against a physical 0.35, the 74
px inner-stack offset, and the 0.8.2 mask failure.

## 0.8.2

- Fixed a merge defect producing two lunar discs in the composite. Mask
  centres are now constrained to a robust straight-line lunar track, with
  a circle-fit-rms gate (1.3x accepted, 3.5x broken).

## 0.8.0 – 0.8.6

- Redundant-link alignment network (lag-1/lag-2 links between tiers,
  solved globally by weighted least squares) and a robust straight-line
  fit of the lunar track predicting each tier's disc position/radius.
- Photometric cross-calibration rewritten with Huber IRLS.
- Hot-pixel mapping across frames with a voting threshold.

## 0.7.2 – 0.7.9

- Sky-cast neutralisation, warmth/tint separation and luminance-preserving
  colour moves consolidated.
- TIFF bracket input.
- Report generation.

## 0.9.10

- Fixed a ring artifact in glare dimming, from sharing radial weight
  with the short-exposure detail layer's smoothstep window (closes at
  1.45 R, a 3.3x ramp). Glare now uses its own exponential profile
  (scale length 0.6 R) with no boundary. Outside the limb (r > 1.1 R),
  peak curvature drops from 21.5 at r = 1.42 R to 1.54, and peak slope
  from 2.61 to 0.97.

## 0.9.9

- Tier TIFF export now embeds an ICC profile: sRGB by default, or
  scene-linear via a checkbox. Untagged linear data was previously
  decoded as sRGB, producing clipped-looking exports.
- Fixed the export normaliser clipping strongly coloured highlights; it
  now divides by measured headroom instead of sensor saturation.
- A `README.txt` now ships alongside tier TIFFs explaining the scale.

## 0.7.1

- Reverted the per-tier lunar masking added in 0.7.0: the noise-dominated
  shortest tiers scattered fits by 124 px on a Lumix dataset, and the
  exclusion-disc union removed most of the inner corona; the limb fit
  then fell back to the gradient fit (R=1376 against a true 625). Spread
  is still measured and reported but no longer applied.
- The gradient fallback now rejects an implausible radius instead of
  writing it to geometry.json and logs when the primary fit fails; the
  limb-profile fix from 0.7.0 is retained.

## 0.7.0

- Fixed a spurious 118 px sinusoid in the limb profile that misaligned
  the disc mask: `fit_limb_rays` used the stale pre-correction centre
  for each pass's profile. On a Lumix dataset the seed was bad (R=1319
  vs true 625) and the correction was still about 59 px, appearing as a
  118 px peak-to-peak sinusoid. The fit now re-profiles at the converged
  centre over 5 passes. After: profile swing is 20 px, residual sinusoid
  is 0.65 px (down from about 59), converging to the same answer from a
  seed 70 px and 2x off.
- Each tier is now masked to its own lunar disc during merge, since the
  Moon moves against the corona during the bracket; every tier's limb is
  fitted in the aligned frame and contributes only where it sees corona.

## 0.6.10

- Sky-cast neutralisation ("Neutralise sky cast", default 0.7): at low
  sun altitude, extinction crushes blue (Lumix dataset: far sky R 0.95 /
  G 1.05 / B 0.70). Chroma beyond 0.72 of frame radius is measured and
  divided out at the chosen strength; the corona's own colour (R 1.66 /
  G 0.88 / B 0.23) survives intact. 0 applies nothing, 1.0 fully
  neutralises, slider runs to 1.2.
- Tint slider (green/magenta), an axis Warmth could not reach.
- Colour moves now renormalise to unit luminance, so Warmth, Tint,
  Saturation and neutralisation change colour without changing
  brightness; Saturation no longer darkens as it increases.

## 0.6.9

- The limb is now fitted on merged luminance rather than short_lum,
  whose blended limb sat far off (Lumix dataset: short_lum 2758.5,
  prominence tier 2707.0, merged 2687.6, a 71 px spread in y); the
  prominence gate keeps its own limb as `prom_geom` (Moon ~20 px off).
- Fixed concentric-ring artifacts in MGN by damping harmonics
  continuously with coverage instead of dropping azimuthal order in
  integer steps. Ring-step rms drops to 0.0014 across the field, inside
  and outside 4R.
- Fixed the composite preview: `L.rmask` was indexed with the RGBA byte
  offset, blacking out three quarters of the frame.
- Frame selection is now a choice (Frames: all / best half / best
  only); averaging (default) gains sqrt(N) in signal-to-noise (Lumix
  sharpness spread: 1.06-1.22 for most tiers).
- Photometric links disagreeing with exposure time by more than 2.5x are
  rejected in favor of shutter speed, fixing a spurious 16.98x factor on
  short tiers.

## 0.6.8

- Limb fit no longer discarded when correct: 0.6.5 required a new fit to
  agree with the old gradient fit within 25% in R; a bad seed tripped
  that guard on a Lumix dataset, falling back to the gradient fit 70 px
  off in y. The fit now runs from two independent seeds and keeps the
  better by its own residuals.
- The disc mask now follows the measured limb per azimuth via a
  720-point profile r(theta) in geometry.json. Margin is now about 2 px,
  down from about 9.
- MGN rebuilt: azimuthally-averaged flattening left a residual gradient
  that suppressed fine structure. It is now flattened by a low-order
  (order 2) Fourier background mu(r,theta). Near-limb detail contrast
  improves 52%, structure-to-noise rises from 0.85 to 0.94, no rim or
  halo. Saved mgnContrast values will read stronger; start lower.

## 0.6.7

- Run report: every run now ends with a summary (exposure stack, EV
  span, alignment residual, sensor defects, limb fit rms, plate scale,
  corona trace depth, brightness range, prominence-gate thresholds, and
  methods used with citations). Written to the progress log,
  `.eclipseforgehdr/report.txt` (plus `report.json`), and as
  `<name>_report.txt` beside every export, recording that render's
  slider values.

## 0.6.6

- MGN and FNRGF can now be switched off, down to 0 (previously floored
  at 0.4); each layer then collapses to a flat 0.5. "Disc mask trim" now
  spans +/-40 px instead of +/-10.
- Fixed FNRGF ring/block artifacts: order now matches the covered arc as
  rings leave the frame, higher harmonics are ridge-damped, rejection
  uses a soft Huber weight instead of a hard sigma-clip, and the model
  is smoothed harder along radius.
- Fixed an MGN rim artifact by subtracting residual azimuthally-averaged
  radial trend from finished detail layers. Bright spike outside the
  limb drops from 0.60 to 0.54 (background 0.536); residual radial trend
  over R to 4R falls 3x.
- Fixed a dark bow at the trailing limb: the real limb is not a circle
  (measured: 3.6 px in y, 4.3 px in x, half-res). geometry.json now
  carries `Rmask = R + 2.5*rms` for disc/detail masks, while `R` remains
  the true limb for the prominence gate.

## 0.6.5

- Limb fit rewritten — the real cause of missed prominences: the old
  estimator's raw gradient maximum lies inside the bright inner corona,
  so R came out too large. On a Lumix dataset it stored centre (2686.0,
  3959.2), R=644.0, where the true limb is (2686.9, 3975.0), R=625.2 —
  up to about 35 px of error, cutting off the bases of the two largest
  prominences. The fit now finds the 50% crossing between disc and
  near-limb corona level per ray and fits r(t) = R + dx*cos t + dy*sin t
  robustly (3.5 px rms over 720/720 rays), falling back to the old fit
  if implausible.
- The diamond-ring/contact frame is fitted the same way.

## 0.6.4

- Hot/dead pixel repair (on by default, "Fix hot pixels" next to
  Denoise): defects are mapped once on the shortest exposure tier, then
  repaired per frame by the median of same-colour neighbours, judged
  against a fitted photon-plus-read-noise model. The run logs how many
  photosites were found. Dust shadows are not addressed; they require
  flat fields.
- Evaluated ACHF-style radius-adaptive kernel scaling for MGN and did
  not ship it: it changed structure-to-noise by about 2% (0.85 to 0.87),
  visually indistinguishable.

## 0.6.3

- Earthshine is now opt-in and off by default: an Earthshine checkbox
  sits next to Denoise; leaving it unticked skips the long-exposure
  stack and earthshine layer (faster run, about 350 MB less cache).
  Changing the checkbox invalidates the cache.

## 0.6.2

- Fixed an IndexError in the earthshine layer on sensors whose height or
  width is not a multiple of 8 (e.g. 3708 rows); binning/upsampling is
  now size-safe throughout (the prominence gate had the same latent
  bug).
- The pipeline now warns when a bracket spans less than 6 EV (normal:
  10-14 EV) or when tiers disagree photometrically beyond their exposure
  ratio.

## 0.6.1

- Fixed an MGN inner-ring artifact from the corona's brightness peak
  being re-sharpened, combined with the flat occulted disc bleeding into
  wide MGN kernels. MGN now excludes the disc from its local statistics
  via normalized convolution and subtracts the azimuthal radial
  brightness profile before normalizing.
- Fixed prominence detection: the colour stack was written as float16,
  exceeding its limit (65504), so the redness test returned NaN and the
  gate came out near-empty. It is now float32, built from a single fast
  tier (<=1/100 s), with the threshold derived from the robust spread of
  the corona's own colour.
- Prominence boost now carries a small positive bias, and the layer
  cache is invalidated when the app version changes.

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
