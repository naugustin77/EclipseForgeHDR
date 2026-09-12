# How our raw handling compares with PixInsight's

PixInsight is the reference most eclipse photographers will measure this app
against, so the differences should be stated by us rather than discovered by
someone else. Its settings below were read off a real installation
(PixInsight 1.9.4 / PCL 2.10.3, RAW module on LibRaw 0.22.1), not from
documentation.

Everything marked **differs** is a choice with a reason. Where the reason is
weak, it says so.

## Option by option

| RAW option | PixInsight | EclipseForgeHDR | |
|---|---|---|---|
| Interpolation | VNG | Malvar–He–Cutler | **differs**; VNG available, see below |
| White balance | Camera (as shot) | **Camera (as shot)** | same since 0.22.63 |
| Black point correction | on | on, per channel | ours is finer |
| Highlights clipping | on | not clipped — flagged and excluded | **differs** |
| Auto crop | on | on (visible area only) | same |
| Auto flip | on | not applied | **differs**, cosmetic |
| Interpolate RGB as 4 colours | off | off | same |
| Super-pixel / half-size | off | off | same |
| FBDD / wavelet noise reduction | off | off | same |

## The three that matter

### White balance: as shot (changed in 0.22.63)

**Superseded.** The argument below was the reason for the old default and is
kept because the reasoning is still worth having; the measurement that
overturned it follows.

This is the one that visibly changes colour, and it is deliberate.

Camera white balance is whatever the photographer had dialled in — auto,
cloudy, tungsten, a custom value. Two people shooting the same eclipse with
different settings would get different colour out of an identical pipeline,
and neither result would be reproducible by anyone else.

Daylight white balance is a property of the camera model, published in its own
metadata and identical for every copy of that body. Using it means the colour
of the output depends on the sensor and the scene and not on a menu. It is also
the physically apt choice here: the corona is photospheric light scattered by
free electrons, which is close to the solar spectrum the daylight coefficients
describe.

Consequence for anyone comparing: our colour will not match a PixInsight
render of the same frames unless their camera was set to daylight.

**What was wrong with that.** The premise is that daylight white balance is a
property of the camera model, published in its metadata. For a body LibRaw has
no table entry for, it is not: `daylight_whitebalance` is DERIVED from the
colour matrix. On the Panasonic DC-S1RM2, normalised to G = 1:

```
                                      R        B      R/B
as shot (544/256/462)               2.125    1.805    1.177
daylight (LibRaw pre_mul)           2.258    1.232    1.833
the camera's own Fine Weather       2.227    1.707    1.305
```

The derived value is 39% short of blue against the camera's own daylight
preset, and 46% short against as-shot. That renders the corona 56% redder in
R/B than as shot, which is a large part of why the app needed a warmth control
below 1 and a corona white balance at 0.66 to look right.

So as-shot is the default from 0.22.63, and the old behaviour stays available
as a setting. The reproducibility argument was real but it was buying
reproducibility of a number that was itself an approximation.

### Highlight clipping: never

PixInsight clips highlights to solid white. For a single image that is a
reasonable default. For a bracket it destroys the thing the bracket exists for:
a photosite saturated in one exposure carries real signal in a shorter one, and
clipping it to white throws that away and then merges the white.

We do the opposite. Saturation is detected on the **mosaic**, before
interpolation, so a pixel built partly from a saturated photosite is caught —
the interpolation kernels have negative lobes, so such a pixel can otherwise
read below any threshold applied afterwards. Those pixels then get zero weight
in the merge and their value comes from a shorter exposure instead.

### Interpolation: Malvar–He–Cutler, with VNG available for comparison

VNG is implemented from the tables in LibRaw's own `misc_demosaic.cpp` and
checked against LibRaw's output on a real frame: after applying LibRaw's 16-bit
clip, the largest disagreement over 44 million pixels is 1.5 counts in 65535,
which is integer rounding. It costs about 1.7x the Malvar path's time.

**It was tried on the reference bracket and it is not the default.** A full VNG
stack against a matched Malvar one showed no visible advantage, and the colour
runs into the limb slightly less smoothly -- which is what the measurement below
predicts, since VNG is gradient-directed and the steepest gradient in the frame
is the limb. It was a toolbar choice in 0.22.63-0.22.69; from 0.22.70 it is
selected by dropping a file named `vng_demosaic.txt` in the folder of raws, so
that anyone checking this app against PixInsight can still match its
interpolation exactly without the option cluttering the toolbar. The run report
says which was used.

Malvar–He–Cutler remains the default, for the reasons below.

VNG dates from the 1990s; Malvar–He–Cutler (2004) has lower interpolation error
on the same data. That alone would be a weak reason to differ, since both are
respectable.

The stronger reason is that in this pipeline the demosaic must not influence the
photometry at all, and with a per-channel linear fit it otherwise does. Measured
on two exposures of the same bracket, fitting the sensor channels directly
against each other and then repeating it after interpolation:

```
what the fit was given                  R/G      B/G
sensor channels, no demosaic           0.9988   0.9997
VNG demosaic, no white balance         1.0069   0.9819
VNG demosaic + camera white balance    1.0333   0.9835
PixInsight (VNG + camera WB)           1.0181   0.9721
```

The sensor channels carry no colour difference between the two exposures. VNG
introduces one, and white balance adds more. VNG is gradient-directed, so it
interpolates red and blue differently in two frames whose gradients differ —
and two exposures of a bracket differ exactly at the limb, where the longer one
is closer to clipping.

So our photometric fit runs on the sensor channels, before any interpolation or
white balance, where the three channels are independent measurements of the
same photons. The demosaic choice then affects only the rendered detail, never
the photometry.

**This predicts that our per-channel numbers will be smaller than a LinearFit
run inside PixInsight on the same frames, and the table above is why.**

### Auto flip

Not applied. Every geometric step here works in the sensor frame, and output
orientation is a render-time choice. Cosmetic only.

## What was checked directly against PixInsight

`LinearFit` was run on two real frames of the reference bracket (1/13 s as
reference, 1/8 s as target) and our implementation on the same two files:

```
channel   PixInsight     ours      difference
G          0.625255    0.624724      -0.08%
```

Green is the honest comparison — the channel PixInsight's raw loading leaves
closest to the sensor. Agreement to 0.08% over 44 million pixels says the
estimator matches. Red and blue differ for the reason given above, and that
difference is a property of the preprocessing, not of the fit.

Our fit follows LinearFit's documented behaviour: minimisation of mean absolute
deviation rather than least squares, with a value-range filter — reject low 0,
reject high 0.92 — applied before it. Both were taken from PixInsight's own
class reference and preferences dialog rather than assumed.
