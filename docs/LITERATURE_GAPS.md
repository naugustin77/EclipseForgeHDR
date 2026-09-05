# What the literature has that we do not

A pass over the papers, theses and reference source in `Literature and
resources/`, plus a search for comparable projects, looking for capability we
are missing rather than for confirmation. Ordered by expected value.

Everything here is either a quotation from a source or a measurement on Nico's
600 mm set. Where a measurement did not support the idea, that is said.

---

## 1. Phase correlation — TESTED AND REJECTED (0.22.27). I was wrong.

**Status: measured on two datasets. Turning phase normalisation on is 66%
WORSE. What the test did find was a different, real 18-22% win.**

I ranked this highest. It was the wrong call, and the measurement says so:

```
ours, plain cross-correlation                     3.02 px
best of 36 regularised phase-correlation settings 5.00 px
the same band-pass, WITHOUT the whitening         2.91 px
```

The band-pass `H` is worth ~0.1 px; the amplitude whitening costs 2.1 px.
Whitening weights every spatial frequency equally, and in a 1/2000 s frame most
frequencies hold only read noise. The thesis recommends phase correlation for
registering images of *comparable* quality — a 14-stop bracket is a different
problem, and `normalization=None` turns out to be right.

The tangential-blur high-pass from `eclipsetools` (item 5 below) was tested in
the same run: 5.40 px against our 3.02. Also rejected.

**What the test seemed to find, and did not:** that our high-pass was too
aggressive, and 0.5 R was 18-22% better. Shipped in 0.22.27 and **reverted in
0.22.28** — on Nico's real run the network residual went 1.17 -> 2.67 px
(half-res) and the per-tier limb spread 8 -> 12 px. The sweep ran on the
*exported* tiers, which are the output of alignment: already registered and
mean-stacked, from a common origin, with no signal-weight mask. It measured
sub-pixel refinement on an easy pair, not alignment.

The phase-correlation result above still stands: that comparison put two
estimators on the same data, so the relative answer survives even though the
absolute numbers came from an easy case.

The original text of this item is kept below, because the reasoning looked
sound and was not.

---

### (original, refuted)

**Status: real gap, well specified, not yet built. The best-supported item here.**

Our own run report already admits it:

> Alignment: cross-correlation of a gradient-flattened log corona, after
> Druckmuller 2009 (ApJ 706, 1605) — but **WITHOUT that paper's phase
> normalisation**: skimage is called with `normalization=None`, so the
> amplitude spectrum is not divided out and this is not, strictly, phase
> correlation.

Druckmüllerová's thesis §4.1.3 states exactly the property we are giving up:

> The phase correlation proved to be a powerful tool [...] **It can register
> images taken with different exposure times, different distribution of diffuse
> light**, can be extended to subpixel precision.

Different exposure times and differing diffuse light is precisely our bracket —
and "differing distribution of diffuse light" is the same phrase the thesis
uses to justify LDIC's `k_i(φ), q_i(φ)`. Plain cross-correlation is dominated
by the brightest low-frequency content, which is exactly what changes between a
1/2000 s and a 1.6 s frame.

Current cost, from the last run: **2.34 px max network residual** at full
resolution, limb fit **2.37 px rms** over 720/720 rays. That is 0.4% of R and
not sub-pixel.

The thesis gives the practical regularised form (eq. 4.13) rather than the
textbook one, because the textbook version divides by an amplitude that can be
near zero:

```
P(x,y) = F^-1 {  H(ξ,η) · F1(ξ,η) F2*(ξ,η) / ( (|F1(ξ,η)|+p) · (|F2(ξ,η)|+q) )  }
```

with `H` a bounded even function (a band-pass window) and `p, q > 0`. It also
covers rotation and scale via the polar amplitude spectrum (Reddy & Chatterji
1996), which we do not need — the thesis says parallactic rotation over a few
minutes is negligible — and sub-pixel extension, which we do need.

**What to do:** implement eq. 4.13 as the link estimator, keep the existing
weighted-least-squares network over lag-1/lag-2 links, and compare network
residual and limb-fit rms directly against the current numbers. This is a
like-for-like test with an existing baseline.

---

## 2. FNRGF: we use order 6 with a hard cutoff; the reference uses ~50 with a taper

**Status: real divergence, measured, and the measurement does NOT clearly favour
the change. Worth doing properly, not worth assuming.**

The zip in the literature folder is the **Delphi source of Druckmüller's own
FNRGF program**, not a description of it. From `ImgProc.pas` and `Settings.ini`:

- `SegmentCount 100`, and `FourOrder := (SegmentCount-1) div 2` → **order 49**.
- Two separate per-harmonic attenuation tables, `Atte[Ave,k]` and `Atte[Dev,k]`,
  applied when the polynomial is evaluated (`PointAttenFourier`) — a smooth
  taper, not a truncation.
- The shipped defaults decay linearly: the mean series to 0.51 at k=49, the
  standard-deviation series to 0.02.
- The segment standard deviation is Bessel-corrected and carries an **additive
  noise variance** inside the square root, estimated from the data by
  `EstimateAdditiveNoiseRing` (median of segment deviations at the outermost
  radii).

Thesis §6.1.2 explains the design and gives an optimum:

> The middle image has about optimal setting of attenuation coefficients (A_k
> set to (1, 0.97, 0.94, ...), C_k set to (1, 0.96, 0.92, ...), **ω = 50**)

and states the constraint plainly:

> using a high order of the trigonometric polynomial **for standard deviations
> which is not in accordance with the order for the averages gives completely
> wrong results**. Using a high order for averages not followed by the standard
> deviations, on the other hand, is not such a big mistake.

and the failure mode of pushing A too far:

> If the A_k s are set too high, it causes **artificial brightenings in
> low-contrast parts of the image. They are false glimmers of the higher-order
> sine and cosine functions.**

**Ours:** `fnrgf_robust(order=6)`, hard cutoff, and the *same* `1e-3·m²` ridge
on the mean and the deviation. By the thesis's own account that is the "too low
for both" regime — safe, but "they only do not make use of the full advantage
of the FNRGF".

**Measured** on Nico's 600 mm set, half resolution, everything standardised per
shell so the comparison is scale-free (amplitude is not comparable between
variants, each normalises by its own σ):

```
variant                                1.05-1.3R      1.3-1.8R      1.8-2.6R
                                       coh / struct   coh / struct  coh / struct
ours: fnrgf_robust, order 6            0.748 / 69.4   0.977 / 52.2  0.916 / 3.3
attenuated, order 20, A -.03, C -.04   0.661 / 84.8   0.965 / 33.1  0.788 / 3.2
attenuated, order 40, A -.03, C -.04   0.657 / 85.4   0.961 / 32.1  0.784 / 3.4
attenuated, order 49, A -.02, C -.02   0.577 / 69.3   0.926 / 24.0  0.622 / 3.0
```

`coh` is radial coherence at lag 5 px, the codebase's own separator of real
radial streamers from texture; `struct` is azimuthal power m 40-250 over m
250-500 on ring-median-normalised data — structure against noise.

The one clear signal is **+23% structure-to-noise in the inner shell** at order
20-40. Everything else is worse, and coherence drops in every shell.

**The test is not clean, and that matters more than the numbers.** The
attenuated variants were written as plain least squares; they do not have our
IRLS/Huber robust fitting, our coverage-matched order, or our ridge. So this
compares "attenuated high order, non-robust" against "low order, robust" and
cannot separate the two effects. A fair test keeps the robust fit and changes
only the order and the attenuation.

**What to do:** add the two attenuation series to `fnrgf_robust` itself,
keeping the IRLS and the coverage matching, and re-run this table. If the inner
shell keeps its +23% without losing coherence, ship it.

---

## 3. The noise floor on the normalising σ is a different quantity from ours

The reference adds an **estimated additive noise variance** to the segment
variance before the square root, and estimates it from the outermost rings of
the actual image. Ours (`fnrgf_robust`) clips the variance at `(0.3·sg)²` — a
floor *relative to the local robust scale* — and MGN clips at an absolute
`0.004` plus a `photon_floor` map.

These are not the same thing. A relative floor scales with whatever the local
scatter happens to be, including where that scatter is real structure; an
absolute noise variance estimated once from the faint field does not. Cheap to
try, and the reference's estimator is fully specified in
`EstimateAdditiveNoiseRing`.

---

## 4. MGN's global term: the paper uses h = 0.7, we use 0.12

Morgan & Druckmüller 2014 §2:

> C_g is included to give contextual information of the largest scale
> structure. We used **h = 0.7** in this work.

We use `global_wt = 0.12`, a factor of six lower. For an eclipse image the
global gamma-transformed term reinstates exactly the radial gradient we work to
remove, so a lower value is defensible — but it is a large departure from the
published value and there is no measurement in the tree justifying the number.
Worth one sweep to find out whether 0.12 was chosen or inherited.

Everything else in our MGN matches the paper: scales 1.25-40 px, k = 0.7, gains
rising to 1 with kernel width.

---

## 5. Comparable projects

Very little is published; Nico is right that people do not share. Two useful
finds.

**`naavis/eclipsetools`** (Python, on GitHub) is the closest comparable
pipeline. It uses:

- **Phase correlation** for alignment, with a specific preprocessing: mask the
  Moon and edges, then high-pass by subtracting a *tangentially blurred*
  version — which removes the radial gradient without touching azimuthal
  structure. Our flattening is radial-profile subtraction; a tangential blur is
  a different and arguably better matched filter for this. Independent
  confirmation of item 1.
- **Linear fits of each exposure to a reference** for the HDR combine — the
  same family as LDIC and as Hill's PixInsight LinearFit, and confirmation that
  a per-image affine transform is the standard approach rather than one scalar.
- **Enhanced unsharp masking with partial convolution**, kernels excluding
  Moon-contaminated pixels — we do this in MGN already.
- One trick we do **not** use: **fit and remove a linear trend before the
  convolution, then reapply it afterwards**, to stop the bright lunar edge
  contaminating a wide kernel. That is a cheap, targeted fix for exactly the
  near-limb band where our artifacts live.

**Viladrich, RCE 2024** (2024 eclipse, 585 frames) is worth knowing for one
structural choice: **star-based registration with ephemeris drift correction**,
rather than registering on the corona or the limb at all. It sidesteps the
lunar drift problem entirely. He also reports a much simpler enhancement that
works — divide the HDR image by a Gaussian blur with radially varying σ plus a
constant — and uses FNRGF via `sunkit-image`, which is a second independent
implementation we could diff against ours.

---

## 6. Confirmation that our merge is now at or beyond Hill's

Hill's talk (transcript, ~line 1250) describes his combine as PixInsight
`LinearFit` of each exposure to a reference, plus a trapezoid weight — a
**global** k and q per image. Since 0.22.26 we fit `k_i(φ), q_i(φ)` per azimuth
segment, which is the thesis version and strictly more general.

One detail of his that we already get right, worth recording because it is easy
to get wrong:

> you determine the weights based on the images **before** you do a linear fit
> [...] because when you do a linear fit to your images the exposure values will
> go greater than one so you'll end up rejecting most of your picture

Our `wsat` is computed from raw `cmax` before any photometric transform, and
0.22.26's azimuthal correction is applied to `rgb` after the weight is already
fixed. Correct as it stands — do not move that order.

---

## Ranked, with what each is worth

| # | item | evidence | expected value |
|---|---|---|---|
| ~~1~~ | ~~Regularised phase correlation~~ | **REJECTED 0.22.27**: 5.00 px against our 3.02 | none — the whitening costs 2.1 px |
| ~~1b~~ | ~~Relax the alignment high-pass to 0.5 R~~ | **REVERTED 0.22.28**: real residual 1.17 -> 2.67 px; the sweep ran on already-aligned tiers | none |
| ~~5~~ | ~~Tangential-blur high-pass~~ | **REJECTED 0.22.27**: 5.40 px against our 3.02 | none |
| 2 | FNRGF attenuation series, robust fit kept | thesis §6.1.2 with stated optima; +23% inner structure/noise, other shells worse | medium, needs the clean test |
| 5 | De-trend before convolution, re-trend after | `eclipsetools`; targets the near-limb band | medium, cheap |
| 3 | Additive noise variance on σ | reference `EstimateAdditiveNoiseRing` | medium |
| 4 | MGN global weight sweep | paper says 0.7, we use 0.12, no measurement on file | low, one sweep |

None of these is aimed at the ring artifact. That remains unexplained, and
nothing in the literature describes it — which is itself worth noting: the
published pipelines do not report it, and they differ from us most in the
alignment step.
