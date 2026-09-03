# Changelog — EclipseForgeHDR

Newest first. Entries from 0.6.1 onward were written at the time. The 0.7.2 –
0.9.8 block was reconstructed afterwards from the code comments that name their
own version and from the development record; where a change cannot be pinned to
an exact version it is filed under the release it is known to precede.

## 0.21.2 — read Druckmullerova's FNRGF source; fixed the measuring stick instead

Nico supplied `FNRGFsoftware.zip` — Druckmullerova's own Delphi implementation,
with `.pas` sources. Read against `fnrgf_robust`, line by line.

**Our FNRGF needs no change.** Same core, `(I - Ave(r,phi)) / Dev(r,phi)`, with
both terms azimuthal trigonometric polynomials fitted per radius. Where we
differ we are equal or better motivated:

| | reference (`ImgProc.pas`) | ours |
|---|---|---|
| azimuthal fit | mean/sd in `SegmentCount` bins, Fourier fitted to the bin values | Huber IRLS directly on 1440 samples — no binning quantisation |
| high orders | `Atte[series,k]` per order, user-set, **default 0** (i.e. plain NRGF) | ridge `1e-3·m²`, always on, coverage-matched order |
| output | `norm(Input) + MixRatio·norm(Mask)` | `fnMix` in the renderer, same job |
| noise | `Noise_AddVar` added inside the deviation | *(absent)* |

That last row looked like the find of the day — an additive noise **variance**
inside the divisor, estimated by `EstimateAdditiveNoiseRing` as the median
per-segment variance over the outermost rings. It is aimed at exactly the
complaint this whole session has been about: not dividing by a sigma that is
really grain.

Implemented and measured, it does nothing for us. Coherence is unchanged to
three decimals in every shell (0.506, 0.452, 0.533, 0.624 before and after);
only the amplitude falls, by 25% in the outer shell. Which is what the algebra
says it must be: our fitted variance at 3.0 R is 2.99e-5 and their estimator
gives 2.06e-5, so adding it scales the outer field by about 1/sqrt(2) — and a
uniform scaling cannot change coherence. In their pipeline that rescale matters
because the mask is min-max normalised before mixing; in ours `fnCompress` and
`fnMix` already own that decision. Not added.

### The actual defect this turned up, in `_fine_structure`

The metric normalised every radial column by dividing by its azimuthal median.
Correct for a MULTIPLICATIVE quantity — luminance, or a 0..1 detail layer with a
median near 0.5. Catastrophic for a ZERO-CENTRED one: FNRGF returns `(L-mu)/sd`,
whose azimuthal median is ~0 by construction, so the guard `max(median, 1e-9)`
divided by about 1e-9 and the metric reported an amplitude of **2.6e8**.

Nothing warned. The number was simply wrong, and wrong in the direction of
looking like an enormous result — which is the dangerous direction, and it is
the third time in one session that a first number was an artifact.

Now the branch is chosen by what the data is: a column whose median dominates
its own spread is divided as before; anything else has the median subtracted and
is scaled by the shell's robust spread. **Verified that every layer the metric
was already used on takes the first branch and is bit-for-bit unchanged** —
the four MGN shells 0.21.0 was built on still read 0.0535, 0.0602, 0.0425,
0.0380.

No pipeline behaviour changes in this release. It is a measurement fix, which on
this project is the same kind of thing as a bug fix.

## 0.21.1 — prominences join the partial-convolution mask

From Jonathan Hill's slide: *"Moon, proms, stars should be 0 in the mask."*

Partial (normalized, incomplete) convolution was already how every masked filter
here estimates a local mean — `(w·f * C) / (w * C)`, the mask convolved by the
same kernel, which is why our limb has never shown the black ring in the left
half of his comparison. What was missing was what goes *in* the mask: ours held
the Moon alone.

Prominences matter more than their area suggests because **S, the local standard
deviation, is MGN's divisor**. A prominence is the brightest thing in the frame
and sits hard against the limb, so every coronal pixel within a kernel of one has
its contrast divided by a sigma the prominence set.

Measured with the limb split by azimuth into rays that contain a prominence and
rays that do not — the second column is the control, in the same run:

```
pct grow  mask %ring   1.02-1.15 near/away   1.15-1.40 near/away
 99   2      1.60%       +6.8% /  +0.2%        -0.0% / +0.2%
 99   5      2.63%      +10.8% /  +0.2%        -0.2% / +0.3%
 99   8      3.73%      +20.2% /  +0.3%        -0.3% / +0.5%   <- shipped
 98   5      5.23%      +24.9% /  +0.3%        +0.6% / +0.4%
 97   5      8.00%      +34.1% /  +0.4%        +4.0% / +0.3%
 99  12      5.27%      +31.8% /  +0.7%        -0.5% / +0.9%
```

Only the near column moves, so this is prominences and not a general consequence
of masking more. Two things the table says that a single figure would have
hidden: the effect is **confined to 1.02–1.15 R**, and it **does not saturate** —
mask more, gain more, with no natural stopping point. So the mask size is a trade
(corona contrast against pixels flattened to 0.5), not an optimum, and 99/8 is a
conservative point on it. The third column is the cost of moving it.

An earlier version of this measurement, quoted mid-session as +26.9% and +24.9%,
used an ad-hoc mask about four times larger in area than the one that shipped.
The 1.15–1.40 R half of it does not reproduce at any mask size below 8% of the
ring. Treat the table above as the result.

Scope: only MGN's own statistics use the tighter mask. `valid` stays disc-only
for RHEF, the deband trend and the inner layers, none of which were measured
against it. Prominence cores therefore come out flat 0.5 in the MGN layer —
Hill's step 4 — and the prominence gate covers 100% of the masked pixels, so
promdet supplies the structure there.

**Stars are not masked.** Hill lists them and the same argument applies, but the
2026 brackets have none, so there is nothing to measure and nothing is claimed.

## 0.21.0 — RHEF carries the outer field

Nico's standing complaint about grain outside the corona, answered from a paper
he put in the reading folder: Gilly & Cranmer, *Visualization of High Dynamic
Range Solar Imagery and the Radial Histogram Equalizing Filter*, Sol. Phys. 300,
174 (2025). RHEF is the whole of their eq. 1 —

```
I_out[A_i] = rank(I_in[A_i]) / N_{A_i}
```

— every pixel becomes its percentile rank within a one-pixel annulus. No kernel,
no scales, no parameters.

Because it has no spatial kernel it cannot compete with MGN near the limb, where
the picture is made of fine multiscale structure. Because it has no local sigma
to divide by, it beats MGN badly further out, where that sigma is mostly
measuring grain. Both layers given the same mild smooth, scored as amp×coh:

| shell | MGN | RHEF | |
|---|---|---|---|
| 1.05–1.30 R | 0.0531 | 0.0318 | −40% |
| 1.30–1.80 R | 0.0628 | 0.0372 | −41% |
| 1.80–2.60 R | 0.0492 | 0.0629 | **+28%** |
| 2.60–3.40 R | 0.0424 | 0.0570 | **+34%** |

and the radial coherence on its own, which is the grain question directly:
**0.489 → 0.752** at 1.8–2.6 R, **0.349 → 0.470** at 2.6–3.4 R.

So they are complementary, and the renderer crosses over on radius rather than
choosing. Weight and crossover were both picked by measurement:

```
rhefMix   1.05-1.30  1.30-1.80  1.80-2.60  2.60-3.40
  0.4        +0%        -0%        +8%       +11%
  0.8        +0%        -1%       +32%       +36%
  1.0        +0%        -1%       +46%       +50%

crossover from 1.3 R    -0%       -15%       +48%      +50%
crossover from 1.6 R    +0%        -1%       +46%      +50%
crossover from 1.9 R    +0%        +0%       +32%      +50%
```

New slider **Outer field (RHEF)** (`rhefMix`, default 1.0) and an RHEF layer
view. At 0 every earlier build is reproduced. The crossover is fixed at
1.6–2.2 R because that is where it stops costing anything inside while still
collecting everything outside.

**No re-stack.** RHEF is one `lexsort` of luminance you already have on disk, so
a work directory from 0.18 onward builds it on first load, in seconds, and
caches it.

The one eclipse-specific adaptation is excluding the occulted disc from each
annulus's ranking — left in, those pixels are a solid block of identical dark
values that would own the bottom of every inner annulus's distribution.

Measured on one dataset, and implemented from the paper rather than from their
`sunkit_image` code, which this machine could not reach. Worth checking against
their implementation when someone can.

### What this replaced

Polar-domain convolution, the last untried item from the Astro Imaging Channel
talk, was tested first and **rejected**. The +155%/+344% figure quoted earlier
was amplitude only, on shells where isotropic MGN's coherence is 0.489 and
0.349 — more than half of that "detail" was grain. Re-scored as amp×coh with
proper azimuthal sampling and against a control that takes the identical
warp/unwarp round trip and then applies the *ordinary* kernel, the polar kernel
loses in every shell (−46%, −42%, −27%, −10%). The apparent gain was the warp's
own interpolation smoothing. Caveat: the test reused the paper's scale list as
polar-grid pixels, which mean something different there, so this rejects the
implementation rather than the idea.

## 0.20.3 — The limb fit on a contact frame found the bead, not the Moon

0.20.2 had the wrong diagnosis. It blamed a scale mismatch on the contact frame
having been shot at a different focal length. Nico's EXIF says 600 mm for both,
and the ratio that story rested on came from circle-fitting a crescent — which
is not a circle. Retracted.

Here is what actually happens, measured on P1072722 (3rd contact) against the
composite it was loaded into:

| | distance from the composite lunar centre |
|---|---|
| the crescent, in the frame as shot | **588 px** |
| the composite's lunar limb | 619 px |
| the crescent, after `prepare_contact` | **190 px** |

As shot, with no registration at all, the frame is already within ~30 px of
right — which is what a tracked sequence should give. `prepare_contact` then
moved it **399 px up and left**, putting the crescent well inside the lunar
disc. Back-solving the shift, the fit placed the lunar centre at (2915, 4248):
on top of the crescent.

The cause is that `fit_limb` was written for a totality frame — a dark disc
inside a corona. A contact frame is not that. 99.4% of this one is below the
noise: dark sky, no corona at that exposure, and one blazing crescent. Given
that, a limb finder fits the crescent's arc, because it is the only edge in the
picture.

So the fit is now **checked before it is obeyed**. If it asks to move the frame
by more than a quarter of the lunar radius, or to rescale it beyond 0.8–1.25, it
has found the wrong thing: the frame is overlaid **as shot**, and the log says
so and points at the ring sliders. For a tracked sequence doing nothing is far
closer to right than obeying a bad fit.

Kept from 0.20.2, on their own merits: the log no longer reports `auto-scale
x1.0000` when it silently gave up, and after alignment the bright arc is checked
against the lunar limb — a crescent of photosphere is always outside the Moon,
so a crescent inside it is a registration failure that should announce itself.
The wider ring sliders stay too (size 0.80–1.25, offsets ±200 px): widened for a
bad reason, but a manual nudge needs more room than ±30 px.

Reload the contact frame to pick this up. The composite stack is untouched.

## 0.20.2 — The contact frame could fail to register and say nothing

Nico exported a composite with a diamond ring and found the ring in the wrong
place. The export was not at fault: the diamond lands within 1 px of where the
ring parameters put it. The cached contact layer was already wrong.

Measured on the stored `contact_rgb.npy`: the ring arc circle-fits to
**R = 553 px centred 327 px from the composite disc**, whose limb is at
R = 619 px. The crescent therefore sits *inside* the lunar disc — and a crescent
of photosphere cannot be there, because the Moon is what hides the rest of it.
That is a registration failure, not a taste question.

`prepare_contact` computed the scale it needed, `geo.R / R`, and then applied it
only when it fell in 0.9–1.1, falling back to 1.0 outside that band — silently,
with the log still printing `auto-scale x1.0000`. The ratio this frame needed
was **1.120**, just outside. 553 x 1.120 = 619, so the refused scale accounts for
the radius error exactly.

Three changes:

* The band is now **0.5–2.0**. A contact frame shot at another focal length is
  an ordinary thing to do — people zoom out for the diamond ring — and refusing
  to scale it is worse than scaling it.
* Outside even that band the log **says so**, in full, instead of reporting a
  scale of 1.0000 as though nothing happened.
* After alignment the frame is **checked**: the bright arc is circle-fitted and
  its distance from the composite disc centre reported. If it lands inside the
  lunar limb the log says the frame is not registered and that the sliders
  cannot fix it. Geometry, not taste — a crescent is always outside the Moon.

The ring sliders were also too narrow to correct anything real: **Ring disc size**
now spans 0.80–1.25 (was 0.96–1.04) and the two offsets ±200 px (were ±30).

Reload the contact frame to pick this up; the composite stack is untouched.

## 0.20.1 — Prominence colour detail is off by default

Nico's verdict on 0.20.0 was immediate: pink prominences look wrong. The
measurement stands — green and blue really are the only channels with room, and
the structure really does show up there — but a correct mechanism pointed in an
ugly direction is still ugly, and taste beats a correlation coefficient on a
question of how a picture should look.

`promChroma` now defaults to **0**, so nothing changes unless it is asked for,
and its range runs **-1.5 to +1.5**. Positive makes dense material go pinker
(what was shipped and rejected); negative makes it go deeper red, which is the
direction worth trying before the idea is abandoned. The whole term is skipped
when the slider sits at zero.

## 0.20.0 — Prominence structure, carried in colour

Nico set Prominence detail to 1.0 after 0.19.0 and reported no visible change.
He was right, and measuring why found the real obstacle.

**Inside a prominence the red channel is at 255 in 100% of the bright core.**
Green sits near 59, blue near 51. Red is the maximum channel everywhere in
there, so the hue-preserving highlight knee sets all three channels from `ms`
and the chroma ratio — which means luminance detail cannot reach the picture at
all, however good it is. The same promDetail change measures 59 levels before
the knee and 13 after it. The 0.19.0 driver swap measures 1.5 levels on screen.
That is what "no visible change" was.

Green and blue have roughly 200 unused levels. Putting the prominence's own
structure there costs nothing in red:

| | G spread (p10–p90) | agreement with H-alpha in G |
|---|---|---|
| promChroma 0.0 | 49–79 levels | 0.368 |
| promChroma 0.6 | 57–105 levels | **0.596** |

and on the two smaller prominences 0.565 → 0.762 and 0.740 → 0.827. On screen
it is a mean of 13.6 levels and p90 of 27, against 1.5 for the luminance route.

New slider, **Prominence colour detail** (`promChroma`, default 0.6). At 0 every
build before this one is reproduced exactly. Outside the gate `prom` is 0 so the
factor is exactly 1 — measured residual 1.2e-07, which is float rounding.

It reads as the denser prominence material going pinker rather than redder.
That is roughly what more continuum through more material looks like, but it is
a look, and the slider is there because it is a taste call.

Measured on one dataset, three prominences. **A 0.18.0 cache is reused as it
stands.**

## 0.19.0 — Let the prominence layer drive the prominence contrast

0.18.0 built a prominence detail layer and blended it into `det`, but left the
*contrast* term — the one that gives a detected prominence its texture — driven
by the inner-corona layer, which was the best thing available before the new
layer existed. It is not any more.

Scored as correlation with the H-alpha red channel's own fine structure, which
normalisation cannot fake, on each prominence's bright core clear of the disc
mask (Nico's 12 Aug 2026 bracket, 600 mm, 14 tiers):

| prominence | inner-driven | promdet-driven |
|---|---|---|
| az 226 (the large one) | 0.714 | **0.828** |
| az 44 | 0.357 | **0.361** |
| az 10 | 0.143 | **0.284** |

No prominence gets worse. The positive bias stays at 0.30: lowering it to 0.10
scores 0.828 against 0.807 on the large prominence, which does not justify
changing how bright everyone's prominences render.

The change is scoped by construction. `prom` is zero outside the gate, so the
whole term is exactly 1 there and nothing in the corona moves; a work directory
with no `promdet.npy` falls back to the old driver and renders bit-identically.

Measured on one dataset. It wants a second bracket before it is trusted as a
general result.

**A 0.18.0 cache is reused as it stands** — this release touches `render.py` and
`gui.html` only, so nobody pays for another 16-minute stack.

Two things this release does NOT fix, both measured while looking for the cause:

* **The disc mask is not covering the prominences.** Where they are actually
  bright, none are hidden and the mask weight runs 0.73–0.92. An earlier reading
  of 83% hidden came from sampling the whole gate, most of which sits inside the
  disc where the colour test fires on noise. The mask is fine.
* **83–91% of the large prominence's bright core still renders above the
  highlight knee**, in every setting tried. The envelope pins at the frame's
  99.97th percentile and a prominence is 5,870 px out of 43 million, so its own
  structure lands where there is no room left. Raising the highlight compression
  does not recover it — measured across 0.1 to 1.0, fine structure inside the
  prominence goes slightly *down*. This is the next real problem.

## 0.18.0 — Prominences get their own detail layer

Both testers said prominence interiors came out flat, and Nico's short-exposure
frames show why that is a loss: the structure is plainly there. Two measured
reasons, on the reference bracket's largest prominence (225 deg, R/GB 16.4).

**MGN's normalisation window clips them.** `hi` is the 99.95th percentile of the
frame — correct for a corona, where a few hot pixels must not set the range —
and a prominence is brighter than that:

    red channel        37% of the prominence hard-clipped at xn = 1.0
    merged luminance   21% clipped, the rest squeezed into the top 6% of range

Whatever structure it had was gone before the multiscale filter ran. Giving the
layer its own window, taken inside the gate, drops the clipping to 2%.

**And MGN is the wrong filter for a compact bright feature.** Its purpose is to
divide out the local standard deviation so faint structure comes up equally at
any brightness — which is exactly what flattens a prominence, whose interior
variation IS its local sigma. Measured as correlation with the red channel's own
fine structure, which no normalisation can fake:

    MGN, frame window (what the corona layers do)        0.037
    MGN, gate window, corona scales                      0.375
    MGN, gate window, fine scales                        0.420
    plain multiscale unsharp, gate-scaled                0.940   <-- shipped

    existing layers on the same measure:  inner 0.317, MGN 0.274, merged 0.262

**And it uses the RED channel, not luminance.** An H-alpha prominence puts most
of its signal in R, and luminance weights R at 0.2126, so the conversion costs
more than half the structure before any filter sees it — 1.000 to 0.417, the
single largest loss in the chain.

So `promdet.npy` is a multiscale unsharp of log(red) from the H-alpha tier,
scaled by the spread inside the prominence gate. It is blended into the detail
only where `prom` fires — the same gate that already decides where prominences
are brightened — so it cannot touch the corona at any setting. New slider,
**Prominence detail**, default 0.7; at 0 the render is identical to 0.17.0.

Measured end to end on the reference bracket, correlation of the render's detail
with the truth inside the prominence:

    promDetail 0.0 (0.17.0)   0.368
    promDetail 0.4            0.535
    promDetail 0.7            0.633
    promDetail 1.0            0.682

**Two things that looked right and were not**, both killed by measurement rather
than argument:

*Partial convolution with the prominences masked out.* Druckmuller's method, and
the obvious candidate — exclude the bright feature so it does not inflate the
local statistics. It makes things WORSE (0.0850 to 0.0561, and 0.0519 dilated),
because excluding a bright feature drags the local mean down, `(x - B)` is then
large everywhere inside it, and the arctan saturates. Partial convolution solves
the DARK plateau problem the Moon creates. Different problem, different tool.

*My own comparison metric.* The first version of this analysis compared
`rms(highpass)/mean(lowpass)` across stages — but that denominator is the local
brightness for an IMAGE and about 0.5 for a normalised LAYER, so "red channel
0.2932 against inner 0.1513" was never a comparison. Every number above uses
correlation against a fixed reference instead, which is invariant to whatever
normalisation a stage applies.

An old workdir has no `promdet.npy`, and the renderer omits the key entirely
rather than substituting a flat 0.5. That distinction is not cosmetic: the blend
pulls `det` TOWARDS the layer, so a flat one would ERASE detail inside the gate
rather than do nothing. The first version of this release did exactly that, and
the fallback test caught it — on a 0.17.0 cache, turning the new slider up would
have flattened the prominences it was meant to sharpen. 0.17.0 caches load fine;
they simply have nothing to blend until the folder is re-run.

Verified: with no layer the slider has no effect at all; with a layer,
promDetail 0 is bit-identical to 0.17.0, and the change is exactly zero wherever
the gate is zero.

## 0.17.0 — Letting the data say how much the short exposures are worth

Nico: "All the detail is in the short exposures. So it is there but somehow got
lost. Is there a way of weighting the short exposures against the longer ones?"

Measured on his own run, comparing `short_lum` (the four shortest tiers) with
the merged HDR, radial profile divided out, and separating structure from grain
by radial coherence -- real corona is continuous from one radius to the next,
photon noise is not:

    shell          source        amp     coh    amp*coh
    1.01-1.10 R    merged     0.0326   0.995     0.0325
    1.01-1.10 R    short      0.2101   0.993     0.2086
    1.30-1.80 R    merged     0.0212   0.976     0.0207
    1.30-1.80 R    short      0.0576   0.186     0.0107
    1.80-2.60 R    merged     0.0135   0.874     0.0118
    1.80-2.60 R    short      0.1560   0.013     0.0020

He was right about the limb: 1.9x more coherent fine structure in the short
tiers, all the way round the disc -- 1.38x in the quietest sector, 3.22x in the
busiest -- and coherence 0.993 says it is structure, not grain. And he was
right that the merge loses it. `w = s * rolloff` weights by exposure time,
which is optimal for photon noise and blind to the fact that a long exposure
beside a blindingly bright edge is glare-smeared long before it clips. At the
limb the longest not-quite-saturated tier carries 67x the weight of the
sharpest one.

He was NOT right about the streamers, and the same table says so: beyond 1.3 R
the short tiers are noise (coherence 0.186, then 0.013). Nothing to recover
there, and a blanket tilt would trade streamers for chromosphere. Also
measured, and also worth knowing: the a-trous denoise is innocent of the
streamer question. Coherent structure retained is 99.9% at 1.1-1.3 R, 100.7% at
1.3-1.8 R, and 107.4% at 1.8-2.6 R -- out there it IMPROVES the ratio by
stripping incoherent grain. So more MGN or NAFE knobs would not have found
anything.

**So the merge now has one dimensionless knob, `w = s**alpha * rolloff`, and a
trial that sets it by measurement.** Four alphas, half-res merges, seconds --
the same pattern `_moon_mask_helps` has used since 0.8. alpha = 1.0 is the
historical behaviour and the default; nothing else is used unless it wins.

**The guard is coherence, not score.** This is the part worth keeping. Giving
weight to a noise-dominated tier RAISES the amplitude in every shell, so a test
that asks "did any score drop" waves it through -- the trial's own noise case
scored +41% at the limb. What noise cannot fake is radial continuity: as alpha
falls from 1.0 to 0.55 on that case, coherence goes 0.933 -> 0.438 at the limb
and 0.398 -> 0.225 in the mid field, monotonically, while in the genuine glare
case it barely moves. So an alpha costing more than 0.05 of coherence anywhere
is buying grain and is refused.

Seven cases, all passing, and two of them were found failing first:

    glare-smeared long tiers            -> tilts to 0.70, +25% at the limb
    every tier equally sharp            -> stays at 1.00
    short tiers replaced by pure noise  -> refused  (failed first: chose 0.55)
    no lunar track                      -> stays at 1.00
    two-tier bracket                    -> stays at 1.00
    shells off the frame                -> refused  (failed first: chose 0.55)
    alpha 1.0 vs the old code           -> bit-identical over all 14 tiers

The off-frame case matters for exactly the range Nico is worried about: a long
lens filling the frame with the disc leaves no outer shell to guard with, and
tilting the merge with nothing watching the outer field is the trade this trial
exists to refuse. It now declines rather than deciding blind.

The fine-structure metric also had to learn to sample at one pixel rather than
a fixed count -- over a narrow shell on a small disc a fixed nr samples at
0.125 px, and bilinear interpolation then makes white noise look perfectly
coherent. That is what let the noise case through the first time.

This is the harness 0.16.0 should have had. The merge is not changed on an
argument; it is changed on a number taken from the data in front of it, with a
guard that refuses the trade when the number is grain.

## 0.16.3 — The disc mask was measuring the corona, not the Moon

Both testers reported losing prominences, the largest one included. On the
reference bracket the mask sat at R+25.1 px and Nico had to set Disc mask trim
to -40 to get them back.

**The margin came from a metric that cannot measure an edge.**
`limb_transition_width` took its 100% reference as `percentile(v, 95)` over
0.75-1.35 R. But the inner corona is still climbing steeply well beyond the
limb -- on this bracket it peaks at R+27 -- so that percentile IS the near-limb
corona peak, and "80% of it" lands far outside the lunar edge.

The proof is that the answer scaled with wherever the reference was put, on the
same merged luminance:

    reference at   R+5   R+10   R+15   R+20   R+30
    median width   6.5    9.8   12.5   16.0   20.0 px
    p90 width      9.0   11.5   16.0   21.0   28.0 px

A real edge width converges as the reference moves outward. This one grows
without bound, which is what says it was measuring the corona's own radial
gradient and not the Moon at all.

The reference is now taken just outside the edge, at R+5 px. Measured:

                        shipped     corrected
    ramp (p90)           28.0 px      11.0 px
    disc mask margin     25.1 px       9.9 px

15.2 px of annulus recovered all the way round -- 23 arcsec at this plate scale
-- which is exactly where the prominences were. On an imported stack whose limb
is genuinely sharp the corrected metric returns the SAME margin it did before
(9.9 px), so the case that was already right is undisturbed.

**Still outstanding, and deliberately not fixed here.** The fitted R is itself
about 8 px too large: the merged image leaves the disc floor at R-7.8 px
(p10 -8.8, p90 -7.3, so it is consistent all round) and the per-tier lunar
radius consensus reads 617 px against the merged fit's 622.9. `fit_limb_rays`
uses a half-level crossing with the same kind of reference, so it is likely the
same bias. That would account for the remaining gap between this mask and the
Moon's true edge. It is left for its own release because R feeds every radial
weight in the pipeline, and after 0.16.1 two changes to the geometry in one
version is not a trade worth making.

Nico's proposal -- find the disc on a short or median exposure rather than the
merged image -- is the right instinct and the numbers support it: the per-tier
fits already agree to within 7 px while the merged fit sits 6 px outside them.
A short tier has a sharp, high-contrast edge and no bright near-limb corona to
inflate the reference. The existing code fits the merged image on purpose,
because the Moon moves against the corona between tiers and the mask must
follow the image being masked -- but the pipeline already fits a linear lunar
track, so it knows where the Moon is at any epoch and can have both.

## 0.16.2 — Reverting the photometric estimator: right reasoning, wrong evidence

0.16.0 replaced the per-tier photometric link with a ratio of sums and no
data-dependent selection. On the reference bracket that made things much worse,
and it put the bright rim back:

                       0.14.5    0.16.1
    1/4000s             1.273     0.571
    1/2000s             1.061     0.746
    1/500s              0.965     0.878
    1.6s                1.158     1.326

    disagreement rim    2 px      128 px
    limb variance       0.072     0.255

128 px of tier disagreement outside the limb is that rim, visible in MGN, NAFE
and FNRGF and faintly in FNRGF's neighbours. The first three links flipped from
0.833 / 0.910 / 1.009 to 1.307 / 1.177 / 1.091 — an over-correction in the
opposite direction to the bias it was meant to remove.

**Why the evidence was bad.** The synthetic tier pairs it was validated on had
a corona filling most of the frame. A real wide-field bracket is mostly sky —
8100x5357 px with the corona reaching 4 R — so the sums are dominated by sky
area rather than by corona, and at the short end the sky carries nothing to take
a ratio of. The estimator was measuring the wrong thing, and the synthetic
scene could not show that.

The original reasoning stands: thresholding a tier against its own noise biases
that tier upward inside the selection, one-sidedly, and it compounds. That is
still what makes a 25-tier FITS bracket read 23 of 24 links below 1.000. So the
old estimator is biased, the new one was worse on real data, and neither is
right.

**What stays.** The per-link residual line and the systematic-lean warning
introduced in 0.16.0 — those are diagnostic only, they cost nothing, and they
are what made this regression legible in one glance instead of a reconstruction
by hand. The FITS colour balance, the Windows memmap fix and the orientation
control are untouched.

**What a replacement has to do**, written into the code so the next attempt
starts from it: take the ratio over a region chosen GEOMETRICALLY — an annulus
just outside the limb where both tiers of a pair carry real signal — so that
the selection depends on neither tier's noise and sky area cannot dominate. And
be validated on the reference bracket AND a wide-field FITS set before it ships.
The mistake was not the idea; it was shipping a change to the core merge on
simulated evidence alone.

Cached products from 0.16.0 and 0.16.1 are not reused — their merged tiers came
from the reverted estimator.

## 0.16.1 — Orientation belongs at the end, not the beginning

0.16.0 put the FITS row-order override at read time, where it was part of the
cache key. That meant changing which way up the picture is cost a full re-stack
— six minutes on the 200 mm bracket — to alter something no measurement depends
on. Nico's read was better: a flip is trivially reversible, so just turn it at
the end.

He is right, and the reason is worth stating. The limb fit is a circle. MGN,
FNRGF, NAFE and Pellett are all radial or tangential about that circle. The
Bayer decode is settled at read time and stays there, because a wrong CFA parity
IS a real error and cannot be undone downstream. Everything else is invariant:
turn the finished picture any way up and every number in the report is the
number it was.

So the read-time override is gone and an **Orientation** control sits with the
export settings — as captured, flip vertical, flip horizontal, 180 degrees,
90 CW, 90 CCW. It applies to the preview and to the export, costs no re-run, and
is lossless: verified as a pure permutation of pixels in all six states, with
cw→ccw and flipv→flipv both exact round trips.

In the preview it is applied at the final canvas blit and nowhere earlier. The
composite is computed in the layers' own coordinate frame against the disc
centre and radius, so turning it any sooner would move the disc out from under
every radial weight.

WHY IT IS NEEDED AT ALL, restated now that the reason is clearer: FITS has no
orientation keyword. ROWORDER describes row order, not which way the camera was
held. A portrait-shot bracket therefore arrives on its side, nothing in the data
can reveal that — a corona is rotationally symmetric enough to have no up — and
so the answer has to come from whoever is looking at it.

## 0.16.0 — Five bugs one FITS log was carrying at once

A 25-tier 200 mm FITS bracket from Val Italo, run on 0.15.5, produced a log with
five separate faults in it. Four are fixed here and the fifth is now something
the user can override.

**The photometric chain was biased by construction, and the bias compounded.**
The log's factors ran 58.202 down to 0.043 and tripped the "tiers disagree
photometrically" warning. Decomposed into per-link residuals — each one is
supposed to sit at 1.000, because the exposure ratio is already divided out —
23 of 24 fell below it, with a median of 0.798. That is not scatter. It is one
error made 24 times, and 0.798^24 is 0.0044.

The cause is the selection, not the arithmetic. The link was
`median(b/a)` over pixels passing `a > floor(a) & b > floor(b)`, and
thresholding a tier against its own noise biases that tier upward inside the
selection — Eddington bias — which drags the ratio down. Measured on synthetic
tier pairs whose true ratio is 1.000 by construction, compounded over nine
links:

    median(b/a), select on both (shipped)   0.785
    median(b/a), select on b only           1.366
    sum(b)/sum(a), select on b only         1.246
    Huber slope b~a, select on b            0.815
    sum(b)/sum(a), NO data selection        1.0015   <-- now shipped

Every data-dependent selection leans toward the tier it selects on; only
removing the selection is unbiased. It holds at read noise 12 e- (0.960) and
with the corona filling a quarter of the frame (0.999). The saturation mask
stays — it keys on saturation, not on noise.

One caveat is now in the code and in the warning: this estimator needs the
black level subtracted. A pedestal does not scale with exposure, so it biases
every link the same way; with 64 ADU left in, sums read 0.005 and the old
medians 0.077. Both are wrong, sums more obviously so.

**The warning now shows the per-link residuals**, not only the running product.
A running product cannot distinguish one bad link from a hundred slightly
biased ones, and that distinction is the entire diagnosis. When most links lean
the same way it says so, names the compounded result, and points at the black
level. Verified to fire on the real 25-tier numbers and to stay quiet on
simulated healthy chains and on a chain with a single 1.5x outlier.

**FITS came out in raw sensor colour.** `FitsFrame` has no white balance and no
colour matrix to give — no FITS convention carries them — so both were identity
and the merge stayed as the sensor saw it. A CFA sensor is far more sensitive in
green than red, so that is not neutral: the corona measured R/((G+B)/2) = 0.56
against 1.11 for a colour-managed stack of the same corona, which is the
blue-cyan rim on a brown sky in that render, and it put the prominence gate's
reference colour somewhere meaningless (0 px flagged).

The balance is now measured from the inner corona, at 1.05-1.6 R, once, AFTER
the merge and the photometry. That is a reference and not a guess: the K-corona
is photospheric light Thomson-scattered off free electrons, and Thomson
scattering is wavelength-independent, so the inner corona carries the Sun's own
spectrum. It is the one thing in the frame that is white by physics. The F-corona
is redder and takes over further out, so the reference stays close in and the
outward reddening the file really contains is left alone. Gains beyond 8x
between channels are refused rather than applied. Raw brackets are untouched —
this runs only where there was no camera white balance to use.

**The Windows sky-gradient failure was the wrong bug.** It reported

    sky gradient removal skipped ([Errno 22] Invalid argument:
    'C:\Users\...\.eclipseforgehdr\hdr_rgb.npy')

It is the opposite end of the same file: `remove_sky_gradient` memory-maps
`hdr_rgb.npy`, copies what it needs, and then np.saves over it while the mapping
is still open. POSIX permits that; Windows locks the file and raises Errno 22.
The map is released before the write now. Verified unchanged on real data — the
same 1.064 / 1.049 / 1.025 per channel as 0.15.3.

This also retires the explanation that has been in `load_big` since it was
written. An earlier Windows report of the identical error carried OneDrive in
its path, and that was taken as the cause — cloud-synced folders do refuse
mmap. It was a coincidence twice over: the folder was named that way but was an
ordinary local one, and `load_big`'s fallback was already shipping in 0.15.5
when this second report arrived on a path with no OneDrive in it. If opening
the map were what failed, that fallback would have caught it. So there was one
bug, reported twice, and the fix added the first time addressed a cause that
did not exist. `load_big` stays — a plain-read fallback is cheap insurance for
filesystems that genuinely cannot map — but its docstring now says what the
evidence actually supports.

**And FITS row order can now be set by hand.** `ROWORDER` is honoured by
default and that has not changed, but the keyword is a convention rather than a
guarantee: a file re-saved by a tool that flipped the data without rewriting the
keyword arrives declaring the opposite of what it holds, and nothing in the file
can catch it — an upside-down corona is still a plausible corona. A FITS rows
control in the toolbar takes "from the header" (default), "bottom-up" or
"top-down". It is part of the cache key, so changing it re-runs, and the Bayer
row-parity shift is applied on the manual flip exactly as on the automatic one.

Also: the corona white-balance reference now decimates against the disc rather
than by a fixed factor, so a 108 px moon gets the same ~15k reference pixels as
a 524 px one instead of being abandoned for want of samples.

Cached products from 0.15.x are not reused: the photometric factors changed, so
every merged tier and every layer built on them differs.

## 0.15.5 — The report was right and still misleading

0.15.4's new prominence line reads "none of it outside the disc mask ... so
nothing reaches the picture". On the run that prompted it, a prominence then
appeared in the composite that had not been there before. Both facts are true
and the sentence was still wrong.

The mask covers everything, not just gate signal. Measured on the 12 px of
annulus 0.15.4 stopped covering, in the sector where the H-alpha excess sits:

    median luminance, 160-200 deg      55253
    median luminance, everywhere else  43950
    excess in that sector              1.26x

That 26% is prominence LIGHT, and uncovering it is what put the prominence on
screen. The gate contributed nothing — it flagged 0 px there, exactly as
reported. So the line now says what is actually true of the gate, and says that
prominences can be visible without it:

    prominences  : corona colour R/GB = 1.11, gate threshold 1.35-1.88, 60 px flagged
                 : none of it outside the disc mask - the gate found redness only
                   at or inside the limb, so the prominence slider has nothing to
                   act on. Prominences may still be visible as brightness; this
                   layer only adds the ones it can identify by colour

Report text only. 0.15.4's cached products are reused — nothing that writes into
the work directory changed.

## 0.15.4 — Three things the import path was doing to somebody else's colour

An imported HDR came back green. Val Italo's 16-bit sRGB stack ran end to end on
0.15.3 and the composite arrived with a teal cast that is not in the file.

**The colour.** Three suspects, measured on that file rather than argued about.
Median chroma, normalised to unit luminance:

                       source           as shipped        temp/tint 1.0
    limb  1.02-1.15 R  1.058 .990 .925  0.832 1.062 .885  1.046 .996 .901
    inner 1.15-1.5  R  1.088 .985 .889  0.853 1.058 .853  1.074 .991 .868
    mid   1.5 -2.5  R  1.155 .973 .811  0.902 1.050 .791  1.133 .980 .803

The per-channel sky-gradient division was the obvious suspect and is innocent:
its tilt is +89 deg, so it cancels in a radial median and moves the corona's
colour by 0.1% or less. `bgNeutral` accounts for 1.5%. The whole of the rest is
`temp` 0.9 and `tint` 1.205 — 21% off red and 6% onto green once the renderer
renormalises to unit luminance, which turns a source corona that is warm in
every shell into one where green is the strongest channel.

Those two numbers are a by-eye correction that sat on top of the camera's
daylight white balance and colour matrix, for one camera, under a sun 7 degrees
up. They are a taste call about the RAW path, not a property of coronae — and an
imported image has already been through somebody's colour management, which is
what makes it an import. So `temp` and `tint` now start at 1.0 for an import and
are untouched for a stack. Neutral reproduces the source to within 2%.

The RAW path renders bit-for-bit identically: verified by rendering the same
layers with the mode marker removed and differencing (max |diff| 0.0).

**The disc mask, on the import path, was guessed.** It was
`max(4.0, 0.042 * R)` — 22.0 px on a 524 px Moon, with no reference to the
image. The stacking path stopped guessing in 0.13: it measures the 20-80%
brightness transition and covers that, because what the mask has to hide is the
half-lit lunar edge, and how wide that is depends on registration and seeing,
not on the disc's size. An import has a *better* limb than a stack, not a worse
one, so the blind rule was backwards. Measured on Val's file the transition is
11.0 px at the 90th percentile, so the same rule the stack uses asks for 9.9 px
where the blind one took 22.0. That is a 12 px ring around the whole limb, about
22 arcsec at his 1.83 arcsec/px, that was being covered for no reason.

**And the report was lying about it.** `stats["geometry"]` on the import path
never carried `Rmask`, so the report fell back to `R` and printed `disc mask at
R+0.0 px` on every import ever run, while the mask was actually at R+22.0.

**Where are my prominences?** Not fixed by any of the above, and the honest
answer is measurable. On Val's file the gate flags 60 px and every one of them
lies between 0.97 and 1.01 R — at or inside the limb, under the mask, invisible.
Beyond 1.02 R the highest R/(G+B)/2 anywhere out to 1.3 R is 1.33 against a 1.35
threshold: there is no H-alpha excess left in that file to find.

That is the gate reading the file correctly rather than misfiring, and the file
says so itself. It carries its Siril header in the TIFF's ImageDescription: 327
frames stacked in Siril 1.4.4, then Background neutralization, Color
Calibration, two GHS stretches (pivot 0.001 amount 145.65, and pivot 0.079
amount 7.57), two background subtractions, and SCNR average-neutral at amount
1.00 — a full-strength green subtraction. What arrives at the import is a
finished picture, and a finished picture no longer carries the chromosphere's
colour separately from the corona's. Lowering the threshold would gate on noise.
So the report says what happened instead of implying 60 prominences he cannot
see:

    prominences  : corona colour R/GB = 1.11, gate threshold 1.35-1.88, 60 px flagged
                 : none of it outside the disc mask - the gate found redness only
                   at or inside the limb, so nothing reaches the picture

**The linearity survived, and that was worth checking.** A history with two GHS
stretches in it is exactly the case the import path warns about, since MGN,
FNRGF and NAFE all assume the value is proportional to coronal brightness. It
holds anyway. Log-log corona falloff over 1.1–2.4 R, measured on the file:

    as stored, treated as linear        -1.52     (reference at display gamma: -1.72)
    after the ICC sRGB inversion        -3.19     (reference scene-linear:     -3.38)

So reading the transfer function out of the embedded profile, rather than
guessing it, put this file within 6% of the reference slope — and had it been
taken at face value it would have been off by a factor of two. Nothing is
clipped at either end either: 0.000% of pixels at zero, 0.000% at saturation.

Cached products from 0.15.3 are not reused: the disc-mask margin is written into
`geometry.json` and every detail layer is built against it.

## 0.15.3 — Disc detection outvoted by noise in the Moon's shadow

An imported stack failed with "could not find the lunar limb in this image" on
two different machines. It was not failing to find a disc — it was finding a
37 px one on a 500 px Moon.

`find_disc` works on the gradient of the LOG of the image, deliberately: a
relative-contrast measure is independent of exposure, of units and of the sky
level, which is what lets one detector work from a 42 px disc to a 420 px one.
But relative contrast has unbounded variance as the signal approaches zero, and
the floor that is supposed to protect against that was the **1st percentile of
the whole frame**. When the disc is darker than the sky and covers more than 1%
of the picture — true of any stack whose sky sits well above black — that
percentile lands *inside the disc*, the flooring never flattens it, and the
shadow's own noise becomes the strongest edge in the image:

| region | value above the floor | median &#124;∇log&#124; |
|---|---:|---:|
| inside the disc | 23 | 0.059 |
| at the limb | 2892 | 0.504 |
| sky | 1447 | 0.011 |

42% of the disc's interior pixels cleared the detector's own strong-gradient
threshold, giving **4548 votes inside the shadow against 2092 on the limb** —
so the circle fit collapsed inward.

The fix asks for signal as well as contrast: an edge means nothing where there
is no light on either side of it. The gradient is now weighted by
`(s−lo)/(s−lo + 0.002·span)`, which is 0.97 at the limb and 0.19 in the shadow.

| | before | after |
|---|---|---|
| the failing stack | no disc found (R=37 px) | centre (1523,2390) R=524 px, **720/720 rays, rms 2.34 px** |
| synthetic sweep, R/short 0.030–0.300 | all pass | all pass, unchanged |

The whole import now runs through to a finished set of layers.

## 0.15.2 — `DLL load failed while importing _rawpy` now says what to do

Windows reports a compiled extension that will not load as `DLL load failed
while importing _rawpy: The specified module could not be found.` — naming
neither the missing thing nor the fix. Both known causes have one-line answers
and neither is guessable from that message, so the app now detects which
applies and says so: the missing Microsoft Visual C++ Redistributable (absent
from a fresh Windows), or an ARM64 Python, for which rawpy publishes no wheel
at all. Both are also in the README's install section now.

Found by trying to run it in a Windows VM on an Apple-silicon Mac, which hit
the second one.

## 0.15.1 — Browse buttons, so nobody has to type a path

Every field that wants a location now has a **Browse…** next to it: the raw
folder, the flats folder, the single-HDR import and the diamond-ring frame.
Folders list how many frames each one holds, so the bracket is obvious among
twenty subfolders; the file pickers offer only files of the right type. Escape
or a click outside closes it. Picking the raw folder loads it straight away.

Both testers lost time to typed paths — one pasted a perfectly good path into a
field whose Start button could never light up — so this is worth more than it
looks.

**Why it is not a native OS dialog.** The obvious answer, `<input type="file">`,
is the one thing that cannot work here: for security a browser hands JavaScript
the file's *content* and never its location, and the pipeline needs the path —
it reads 250 raw files off the disk, it does not want them uploaded to itself. A
native dialog opened by the server would give a real path, but Tk must own the
main thread on macOS while this server is threaded, so that trades a paste for a
hang. The server therefore lists directories and the page draws the picker:
plainer than a native dialog, and identical on all three platforms, which for a
tool being tested on machines its author cannot reach is the better trade.

## 0.15.0 — First external FITS run: four defects it found

The first tester to complete a run did it on FITS from a 249-frame bracket, and
found four things at once. All four are fixed; two more are diagnosed and not.

### FITS came out upside down

FITS row 1 is the **bottom** of the frame — the standard's origin is lower-left
— and every other format the pipeline touches has row 0 at the top. Read
straight through, a FITS bracket comes out vertically mirrored. It now reads
`ROWORDER` when the capture software wrote one and follows the standard when it
did not, and says which in the log.

The fix has a trap in it, and the trap is worse than the bug. Flipping an image
with an **even** number of rows swaps the Bayer row parity — RGGB becomes GBRG —
so the obvious fix trades an upside-down picture for a miscoloured one. The
parity shift now travels with the data and is added to the Bayer y-offset.
Tested at even and odd heights, with `ROWORDER` absent, `BOTTOM-UP` and
`TOP-DOWN`: bar at the top and red on the red photosites in all four.

### One OneDrive folder silently disabled three features

    sky gradient removal skipped ([Errno 22] Invalid argument:
    'C:\Users\...\OneDrive\Desktop\...\hdr_rgb.npy')

`np.load(mmap_mode="r")` fails on a cloud-synced Windows folder. Memory-mapping
keeps a 45 Mpx three-channel float32 out of RAM, which is worth having — but the
same call sat in **three** places: the sky-gradient fit, the contact-frame
loader, and the renderer's colour estimate. So one folder could take out three
unrelated features, and only one of them said so. It now falls back to a plain
read, which costs memory and keeps the feature.

That also explains the intermittent contact frame: the same error, on a code
path whose failure the GUI reported only as nothing happening.

### The diamond ring did nothing — for both testers

The contact-frame path field had no Enter handler. The folder field directly
above it loads on Enter, so pasting a path and pressing Enter is the obvious
move, and it did precisely nothing. Both testers reported this independently,
in the same words, and neither had a broken file. Enter now loads it, and stays
inert while the button is disabled.

### Importing one HDR needed a raw folder that an import does not use

Start was enabled only by a successful folder load, so pasting a TIFF path and
nothing else left the button greyed out forever — reported verbatim. An import
needs no raw folder: its products belong beside the image. Typing an import path
now enables Start, and the server falls back to the image's own directory.

### Diagnosed, not yet fixed

**Prominences cut off by the disc mask.** The merged limb fit came out
R = 542 px against a per-tier consensus of 526 — 3% large — and the mask sits at
R + 25, so it lands about 41 px beyond the real Moon and eats the prominences.
The tier-consensus cross-check exists for exactly this but only fires past 15%.
Shrinking the mask reveals a grey ring, which is the merged limb's own 26 px
transition; that transition is the real fault and 3% is a symptom, not the
cause.

**Photometric calibration ran away on this set** — factors 9.40, 7.08, 5.36,
4.07 on the four shortest tiers, a smooth geometric progression rather than
noise. The pipeline detected it and said the merge would be unreliable, which is
the check working. What it means for a FITS bracket from this capture chain is
not yet known, and a wide limb transition is exactly what a bad photometric
merge produces — so this and the mask are probably one problem, not two.

## 0.14.6 — The last step of a run was never timed

A step's cost is the gap to the next log line, so the final step of a run has
never had one — and has therefore never been measured. On the first full
step-by-step report from a real run, the entire Pellett layer was simply absent
from the list and the "steps under 1s" remainder quietly absorbed it. Since the
whole point of that report is to let real runs replace estimated progress
weights, a hole in it is worth closing before it collects any more data. The
interval is now closed before it is measured.

The same report confirmed the 0.14.4 weighting from the other direction. With
every step over a second listed, the detail stage on a 50-frame 45 Mpx run comes
to 343 s of 734 s — **47%** — against the 48% the bar now allocates. The earlier
top-six estimate could only bracket it between 42% and 64%.

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
