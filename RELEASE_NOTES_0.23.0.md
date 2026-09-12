# EclipseForgeHDR 0.23.0

The first tagged release since 0.22.33, and a large one. The short version: the
photometric calibration was rebuilt and is now fitted per colour channel, a new
detail layer after Jonathan Hill's partial convolution was added, three
estimators were found to be reading each exposure in the wrong frame, and a
merge measurement turned out to have been switched off by an unrelated
decision on every dataset in the project.

Everything below was measured on real brackets: a 600 mm Panasonic S1R II set,
and 250 mm / 360 mm Canon sets from a tester. Where a number is quoted it comes
from one of those runs, not from a synthetic.

## Photometry

**Per-channel linear fit, now the default.** Each pair of neighbouring
exposures is fitted with a slope and an offset, separately in R, G and B, after
PixInsight's LinearFit as used in Hill's HDR eclipse workflow. This replaces
one luminance scale factor per exposure. It is a dropdown in the toolbar —
"Per-channel (linear fit)" or "Single scale factor" — because both are worth
having and because a setting nobody can see is a setting nobody can check.

The chain is built outward from the middle exposure rather than from the
longest, so on a twelve-stop bracket the error accumulates over six links
instead of eleven.

**A shared pedestal**, fitted in the outer field and subtracted from every
exposure before the pairwise fits. A black level left a few ADU behind arrives
divided by the exposure time, so it is nothing on a long frame and everything
on a short one. On the 360 mm set it takes tier-to-tier disagreement in the
outer field from 27.6 % to 2.9 %.

**A per-exposure azimuthal correction** — gain and offset varying around the
disc, fitted in 60 segments and smoothed by a trigonometric polynomial, after
Druckmüllerová's LDIC (doctoral thesis eq. 4.15), applied mean-preserving so it
corrects only what a single scalar cannot express. On the 600 mm set the worst
exposure varies 34 % around the limb.

## A new detail layer: partial convolution

After Jonathan Hill, "Advanced Solar Eclipse Photography". Unsharp masks of the
log-mapped image at 2, 4, 8, 16 and 32 px, combined linearly. Two things make
it different from the filters already here:

- the blur behind each mask is taken in polar coordinates, so it averages
  **along** a streamer rather than across it;
- it is **partial** — the occulted disc and the prominences are excluded from
  both the numerator and the normalisation, so nothing is smeared out of them
  into the corona.

MGN, NAFE and RHEF all normalise locally, which is what lets them show the
inner and outer corona at once and also what gives every part of the frame the
same texture amplitude. This one is linear: faint structure stays faint.

It has its own monochrome view and exports as a 16-bit TIFF like the other
layers, a soft threshold against a per-pixel photon-noise model, and five scale
sliders that set balance while a master sets strength.

It is **experimental**. Known issues are listed at the end.

## Colour

**Neutralise corona** and **Neutralise sky cast** sliders. The first makes the
corona white by physics rather than by convention, taking its reference at
1.05–1.6 R and leaving the genuine outward reddening alone. Both are
multiplicative, which is their limitation — see the known issues.

**White balance is a toolbar setting**: camera as-shot (the default), daylight,
or none. It scales the raw channels before the camera matrix, so it moves the
colour of the whole frame together.

## Fixes worth naming

**Three estimators were reading each exposure in the wrong frame.** The shared
pedestal, the per-exposure radial check and the LDIC azimuthal fit all sampled
rings centred on the middle exposure's centre, but the per-exposure alignment
shift is applied later, in the merge. An exposure offset by *d* px had its
radial profile read at `r + d·cos(φ − φ₀)` — which for LDIC is a dipole in the
very quantity it fits. A bench with known shifts measures a spurious 65 % gain
spread at 5 px of drift where the truth is 0 %.

**A merge measurement was dead code on every real dataset.** The trial that
chooses how hard long exposures may outvote short ones was reading a variable
the caller cleared whenever the *moon mask* was rejected — an unrelated and
routine outcome. With a correct lunar track, no bracket in this project accepts
that mask, so the trial had never run on any of them. Fixing it recovers, on
the three sets: exposure exponent 0.55 / 0.70 / 0.55, +87 % / +459 % / +65 %
coherent detail just outside the limb, and merged limb transitions of 8.0 px,
3.0 px and 8.0 px.

**The merge was resampling every exposure with bilinear interpolation**, and
the clipping test was asked of the demosaiced result rather than of the mosaic
— so a pixel reconstructed in part from a saturated photosite could read below
threshold and enter the merge at full weight.

**A featureless bright collar around the Moon.** Just outside the limb the
partial-convolution kernel is one-sided: every valid sample lies further out,
where the corona is fainter, so the residual comes out systematically positive.
Measured at 7× the layer's own structure at 1.03–1.08 R, now 0.55.

**FITS**: the saturation ceiling no longer collapses onto the frame's own
maximum, planes-first RGB cubes (what Siril and PixInsight write) import
correctly, and odd dimensions no longer raise. All three are pinned as tests.

**Windows**: a memory-mapped array kept a file open across the write that
followed, which failed as `[Errno 22] Invalid argument` on the path being
written.

## Removed

**The sky subtraction**, built across 0.22.51–0.22.60 and taken out again.
Every failure was in the *shape* of the sky, not its level or its colour: any
model with enough freedom to describe a real horizon sky has enough freedom to
eat the outer corona.

**VNG demosaic** comes off the toolbar. It was added so a render could be
compared with PixInsight's on equal terms; that comparison has been made and
showed no advantage.

## Known issues

**The sky's own light is still in the picture.** `remove_sky_gradient` divides
out the smooth gradient across the frame; it never subtracts a level. Measured
on the 600 mm and 560 mm sets, one corona of constant colour plus one
background of constant colour and level reproduces the radial colour behaviour
of both to 1–2 %. So the outward colour change is the background taking over,
not the corona changing — and the multiplicative sliders cannot remove an
additive term. Whether removing it belongs in this app at all is an open
question: fitting a real sky needs a measurement no eclipse bracket contains,
and every attempt so far has cost more outer corona than it bought.

**Partial convolution leaves a residue at the prominence boundary.** On the
600 mm set the prominence gate covers azimuths 44°, 220–230° and 306°, and the
narrow negative lines in the masks sit at 46°, 233° and 307° — every one at a
prominence, at every scale, growing with scale. Hill's own slides show the same
thing. The mechanism is not yet understood and no fix is guessed at.

**Per-channel photometry is measured but not applied to the merge.** The merge
still uses one factor per exposure on luminance. The per-channel factors are
recorded in the diagnostics bundle.

**One tester set is not solved.** A 2024 560 mm bracket shot through thin cloud
keeps a colour cast that none of the above removes. It is not in the validated
set for this release.

## Upgrading

A 0.22.87 work directory is reused as it stands. Anything older re-stacks — the
merge itself changed several times across this span, and reusing those products
under this build would show an old merge under a new report.

```
pipx install --force .
```

## Methods and citations

Every method here is published and any patents that existed have expired. MGN
(Morgan & Druckmüller 2014), RHEF (Gilly & Cranmer 2025), FNRGF (Morgan, Habbal
& Woo 2006; Druckmüllerová et al. 2011), NAFE (Druckmüller 2013), the
tangential filter (Druckmüller 2009), LDIC (Druckmüllerová, doctoral thesis),
and the partial-convolution chain after Jonathan Hill, "Advanced Solar Eclipse
Photography". Full references in the README.

Thanks to the testers who ran their own brackets through it and reported back
with the failures rather than the successes — most of the fixes above started
as somebody's screenshot.
