# EclipseForgeHDR — open items

Things known to be wrong, unfinished, or untested, with the evidence for each and
what would settle it. Kept here rather than in changelog prose so they stop
getting lost between releases.

Ordered by expected value, not by effort. An item moves out of this file when it
is either shipped (into CHANGELOG.md) or **rejected with a measurement** — a
rejection is a result and belongs in the source next to the thing it rejects, as
ACHF and the RHEF gate are.

---

## 1. Per-tier photometric correction — RESOLVED (0.22.15 + 0.22.16)

Both halves are shipped, and the second half was not what this item thought.

**The additive half (0.22.15).** A single pedestal shared by every tier, fitted
from tier overlap in the outer field and shrunk by its own jackknife
uncertainty. Real runs: 360 mm 26.8% -> 3.2%, 2024 560 mm 10.6% -> 2.7%,
250 mm 3.2% -> 2.0%.

**The near-limb half (0.22.16) was never photometric at all.** It was the merge
weight's Gaussian feather leaking into each tier's clipped region, so blown
tiers dragged the merged inner corona down over a band just outside the limb —
5.4x on the 2024 560 mm set, 3.2x on the 360 mm, with a bright overshoot at the
limb itself. A controlled two-tier merge from the raws pinned it and the fix
restores the unfeathered profile to within 0.1%. See CHANGELOG 0.22.16.

No `k_i(r)`, `q_i(r)` per radial band was needed for either. The radially
varying affine transform this item was designed around is not built and, on the
evidence of five datasets, is not warranted.

**What this leaves open.** The 250 mm set's limb variance is still 0.793 with a
62 px disagreement rim, and 0.22.16 does not touch it — that rim is a
tier-to-tier disagreement in brightness, measured on the tiers themselves, not
an artifact of how they are combined. Whether it survives the corrected merge is
the first thing to check on the next run of that folder. If it does, the veiling
glare hypothesis stands and wants its own item.

## 1c. The rings in the 250 mm set, cause unknown

**Status:** four candidates excluded by measurement, none left standing.

Concentric arcs in 1.02–1.15 R on Clifton's 250 mm set, 2–5.3x the reference
set's, worst in NAFE — the one layer that uses neither the limb fit nor the disc
mask. Still present at 0.22.17.

Excluded so far, each with a measurement:

* **the limb fit** (0.22.14) — R corrected onto the tiers' consensus; that set's
  fit was only 2.3 px large and the correction did not fire.
* **the pedestal** (0.22.15) — fitted at −1.26 ADU, 1.4 sigma, shrunk to −0.83.
  Outer field only; nothing near the limb.
* **the merge weight** (0.22.16) — the feather leak is fixed and the rings
  remain. A reconstruction at knees 0.87 → 0.25 moves the near-limb azimuthal
  modulation by under 4%.
* **tier disagreement** (0.22.17) — the 0.79 was a contaminated statistic; the
  tiers agree to 0.3–2.1% ring by ring, better than the 360 mm set.

What has NOT been tested, in order of what I would try next:

1. **The flat.** This is the only one of the five datasets with a master flat,
   and the only one with the rings. It corrects a 105.6% falloff — the dimmest
   part of the field at 0.486 of the brightest — and the Moon sits within a few
   hundred px of the frame centre, so the flat's own radial structure is very
   nearly *concentric with the disc*. A flat divides every tier equally, so it
   cannot show up in any tier-to-tier statistic, which is exactly why every test
   above is blind to it. Re-run the folder with `flats/` renamed away: if the
   rings go, this is it, and the master flat's own radial profile is the next
   thing to look at.
2. **The demosaic** near the clipping ceiling, where one channel clips and the
   others do not.
3. **The starlet denoise** — its per-pixel sigma is a photon-noise model, and
   the flat's division changes the noise per pixel without changing the signal.

## 1a. A per-tier error that varies with radius — MEASURED on four datasets

0.22.24 put the measurement in the pipeline. All four sets, run by Nico:

```
set                    worst tier departure inside 1.3 R   rings seen?
Clifton 250mm            1/60s   +230%                     strong
Nico 600mm               1/2000s  +28%                     yes, "artifact back"
Clifton 360mm            tiers agree to 12%                some, "likely from Inner"
Clifton 2024 560mm       tiers agree to 10%                some, plus a colour shift
```

**The effect is real and it is enormous where it appears.** 230% on the 250 mm
set is not a subtlety, it is on the set with the worst rings, and the
photometric factor — one scalar per tier, fitted where tiers overlap in the mid
field — cannot see it. That is this item, confirmed, and it is worth building
the correction for.

**But it does not explain everything, and that is the important half.** Two sets
whose tiers agree to 10-12% still show concentric rings. So a per-tier radial
correction would fix the 250 mm set and leave the others as they are. See 1d.

**Design, when it is built:** fit `k_i(r)` per tier over radial bands against
the median across tiers, with a strong prior toward 1.0, only where a tier has
signal and is unclipped, smoothed in radius. The profiles already computed by
0.22.24 are exactly the input; nothing new needs measuring.

---

## 1d. REJECTED: the render's fixed-radius blends are not the artifact

`wI` does fade the inner layer out across an absolute radius (1.05-1.45 R) and
`wf` fades FNRGF in across another, and both are the failure mode the codebase
already documents for `wG`. But they cannot be what Nico is reporting, and it
took him telling me twice before I checked: `render()` returns at

```python
if view == "mgn":   return ...        # before wf and wI are ever computed
if view == "fnrgf": return ...
```

so the layer views bypass both blends entirely. The artifact is visible in the
MGN and FNRGF views, therefore it is upstream of the renderer — in `hdr_lum`,
i.e. in the merge. Nico's own bisect said so before this item was written:
`ECLIPSEFORGE_FEATHER=plain` was *"basically gone"*, and that switch touches
nothing but the merge weight.

Leaving the two blends alone anyway is not defensible on this evidence, only
un-prioritised: they are still boundary-shaped, and if a ring ever survives in
the composite that is absent from every layer view, this is the first suspect.

---

## 1e. The merge bench — what has been excluded, and what is still open

**Status:** measurement infrastructure built; cause not yet found.

`tools/feather_bench.py` and its siblings rebuild the merge from the exported
aligned-tier TIFFs, so any weight form can be tried in about twenty seconds
with no re-stack. This should have existed before 0.22.16 shipped. Run on
Clifton's 360 mm set, all figures against the unfeathered merge:

```
variant     1.02R   1.06R   1.20R    azimuthal power, 1.05-1.25 R (m20-80)
none        1.000   1.000   1.000        9.08
plain       0.126   0.486   0.996        1.00
masked      1.001   1.000   1.000        9.16
taper       1.002   0.994   1.000        9.38
```

**Excluded as the cause, each by its own measurement:**

- *The feather.* `none` — no blur, no mask, no contour term of any kind —
  scores the same as the shipped taper. The weight's shape is not the source.
- *Short-tier noise.* Single tiers score 28x apart in this band, scaling with
  how short the exposure is, which is shot noise. But a proper inverse-variance
  weight, `w = (s·cal)²/(x + x_r)`, changes the merged result by under 1%
  everywhere: `s¹` already suppresses the short tiers there. See `invvar.py`.
- *My own candidate.* A weight rolled to exact zero, C1, at 0.90 sat — below
  the 0.97 validity threshold, so the mask has nothing left to step on — is
  indistinguishable from the taper (rms 0.3493 against 0.3492). Rejected
  before shipping, which is the entire point of building the bench.

**What survived:** the per-tier radial error of item 1a accounts for about 26%
of the structure in that shell (1.00 -> 0.74 with `k_i(r)` applied), fitted only
from rings where a tier is ≥98% unclipped. Real, worth building, not the whole
story.

**What `plain` actually buys.** Not a cure — a smother. It suppresses the
metric by leaking weight to the long tiers near the limb, and it costs a factor
of eight in inner-corona photometry: 0.126 of the true level at 1.02 R. On the
360 mm set the visible difference between plain and taper is almost entirely
one bright halo ring at the limb.

### On Nico's own 600 mm set (14 tiers, alpha 0.55) — the rings are located

The tier export exists now, so the bench runs on the set that actually shows
the artifact. The difference image between the plain blur and the shipped taper
is unmistakable: **nested closed contours, one per tier, following the corona's
own isophotes.** They are each tier's saturation contour, printed.

Ring metric — polar transform, high-pass along RADIUS, power at m < 20. A ring
is sharp in radius and smooth in azimuth; a streamer is the opposite, so this
separates them. All figures against `none`, the exact unfeathered merge:

```
variant     1.02R   1.06R   1.20R   1.50R   ring power
none        1.000   1.000   1.000   1.000       1.00
plain       0.752   0.932   0.997   1.000       0.63
taper       1.006   1.001   1.000   1.000       1.11
nonorm      1.002   1.003   1.000   1.000       0.97
cont_pl     0.753   0.934   0.998   1.000       0.63
iv          0.992   0.963   0.999   0.983       1.15
```

**`none` already scores 1.00.** The rings are in the merge with no feather of
any kind. Nothing about the weight's shape creates them; blurring the weight
softens them by 37%, which is the whole of what `plain` buys, and it costs 25%
of the true brightness at 1.02 R.

**Ruled out, each by measurement on this set:**

| candidate | result |
| --- | --- |
| per-tier radial photometry, `k_i(r)` | rms 0.016 against the rings' 0.272 |
| per-tier colour (one scalar `cal[s]`) | ≤ 8%, and flat where the rings are |
| the `1/den` renormalisation | 13% of the ring power (see below) |
| inverse-variance weighting | 1.15x — worse, not better |
| a C1 weight rolled to zero at 0.90 sat | indistinguishable from the taper |
| per-tier nonlinearity vs signal level | ≤ 1-3% for every adjacent pair |

### FOUND: a large photometric error in a collar around each tier's saturation

Every correction above treated a tier as *globally* wrong — a scale, a colour,
a nonlinearity. None of them can touch an error whose SHAPE is the saturated
region itself. There is one, and it is big. Adjacent tier pairs on Nico's
600 mm set, the longer tier's radiance estimate over the shorter one's, binned
by distance OUTSIDE the longer tier's saturated region:

```
pair                 0-2px   2-4px   4-8px  8-16px  16-32px  ...  128-256px
0.333s / 0.2s        2.964   2.044   1.210   0.991    0.991         0.993
0.125s / 0.0769s     1.888   1.958   1.019   0.991    0.988         0.991
0.2s   / 0.125s      1.103   1.088   1.014   0.994    0.990         0.993
0.5s   / 0.333s      0.797   0.946   1.000   0.989    0.993         1.001
1.6s   / 0.5s        0.348   0.403   0.667   0.969    0.983         1.013
```

Beyond about 8-16 px every pair agrees to ~1%. Inside it they disagree by up to
**3x, and the sign differs between tiers.** High is what charge spill or veiling
glare off a large saturated area does; low is what a demosaic does when its
kernel reaches into saturated photosites — `cmax` is measured after demosaic, so
a pixel whose own photosite was fine but whose neighbour clipped is silently
interpolated from a clipped value and is never flagged.

This is contour-shaped by construction, it lives in the tier data, and no
reweighting can remove it — which is precisely why eight weight forms all
failed. It is the first candidate whose shape, amplitude and location all match
the artifact.

**Next:** establish which of the two mechanisms dominates (they are separable —
spill scales with the saturated area, demosaic contamination does not), then
decide whether the fix is to widen the invalid collar, or to measure `cmax`
before demosaic on the Bayer plane where clipping actually happens.

### 2. THE MERGE FORMULA ITSELF IS MISSING A TERM (Nico's question, and he is right)

Druckmullerova, *Doctoral thesis* §4.1.4, eq. 4.15 — the LDIC composition that
Hill says works best:

```
g(r,phi) = SUM_i  w(f_i(r,phi)) * ( k_i(phi) * f_i(r,phi) + q_i(phi) )
```

`k_i` and `q_i` are an AFFINE transform per image that VARIES WITH AZIMUTH,
fitted by linear regression in 60 angular segments against the composite
accumulated so far — longest exposure first — then smoothed with a
trigonometric polynomial of order <= 4. The thesis states what they are for:
*"to compose images with different distribution of diffuse light in the optical
system ... or even images that were taken through thin clouds."*

EclipseForgeHDR has **one scalar `cal[s]` per tier and one shared additive
pedestal.** That expresses the median of `k_i` and nothing else. Fitted on
Nico's 600 mm set exactly as the thesis describes:

```
tier        k median   k spread across azimuth   |q| as % of local signal
0.5s          0.923            23.8%                      9.28
0.33s         0.967             3.9%                      4.62
0.125s        0.980             5.0%                      3.22
0.0167s       0.994            10.6%                      1.49
0.002s        0.998            13.5%                      1.41
0.0005s       1.023            27.3%                      2.34
```

k varies with azimuth by **3-27%** depending on the tier and q runs to **9%**
of the local signal. Both are real terms we do not correct, both are azimuthal,
and an azimuthal error changes wherever the tier mix changes — which is along
the saturation contours where the rings are. Note also that the median k sits
BELOW 1 for eleven of twelve tiers, so `cal[s]` is biased as well.

**Also note our item 1a had this as a RADIAL correction `k_i(r)`.** The
reference method makes it azimuthal. That is very likely why fitting `k_i(r)`
bought 26% on one set and nothing on this one.

**Implemented and measured (`tools/ldic_merge.py`):** the full LDIC
composition, sequential from the longest exposure, 60 segments, order-4
trigonometric smoothing. Ring power **0.87x** against the shipped merge. Real,
worth having, not a cure — and it shifts the radial profile by ~20%, since it
references everything to the longest exposure rather than to exposure time.

**Instrument note.** These numbers are now trustworthy: the bench reproduces
the pipeline's own `hdr_lum` to 0.7-2.5% at every radius, and running the
pipeline's real MGN path on the bench's output gives plain/cur = 0.70x against
the crude metric's 0.66x. The crude metric was not the problem.

**Alignment** (Nico asked): 2.34 px max network residual at full resolution,
limb fit 2.37 px rms over 720/720 rays — 0.4% of R. Not sub-pixel, worth
improving, but a shift error is not ring-shaped and does not explain this.

**Still untested:** per-tier PSF differences,
and per-tier flat residual. A tier-pair ratio binned by signal level does show
an additive residual on the two longest tiers (1.09 at low signal falling to
0.97 at high), which the single shared pedestal of 0.22.15 cannot absorb — that
is a real finding but it lives in the faint outer field, not where the rings are.

**One free improvement, measured and shippable:** dropping the `1/den`
renormalisation (`nonorm`) is better on BOTH axes than the shipped taper — ring
power 0.97 against 1.11, radial profile 1.002 against 1.006. It is not the fix
and must not be sold as one.

---

## 1b. The limb fit on a soft edge — MOSTLY SOLVED BY 0.22.16

The soft merged edge this item was about was largely **the feather leak**, not
the optics. Clifton's 2024 560 mm set, same frames, before and after:

```
                       0.22.15    0.22.16
merged limb ramp p90     28 px      2 px
limb fit rms           0.96 px   0.46 px
fitted R vs consensus   +28 px    +0.7 px   (the correction no longer fires)
disc mask margin       25.2 px    1.8 px
```

A blown tier bleeding weight into the band where it clips was smearing the
merged limb by more than an order of magnitude in ramp width, and the half-level
crossing was running large on the smear it produced. With that gone the fit
lands on the tiers' consensus by itself.

**What remains:** the 0.22.14 correction is still the right guard — it did fire
on the 360 mm set at 0.22.16 (−14.2 px) — but the bias table in its changelog
entry was measured on data that carried the feather leak, so those ramp numbers
describe the old merge, not a property of the optics. Worth re-measuring across
the sets at 0.22.16+ before anyone leans on it again.

## 2. `bgNeutral` does almost nothing to the sky

**Status:** confirmed defect, unfixed. Cheap.

The "Neutralise sky cast" slider divides by the measured sky chroma weighted by
`cconf`: `a / (1 + cconf * (bg_chroma**bgNeutral - 1))`. Mean confidence in the
far sky **measures 0.015** — the code says so itself, in the comment right above
the measurement. So the correction is applied at ~1.5% strength exactly where the
sky is, and at full strength where the corona is.

An earlier fix (0.11.4) corrected how `bg_chroma` is *measured* — it had been
measured on `ratio`, which is driven to 1.0 out there — and left the weighting
alone. So the slider still does not do what its label says.

**What is unclear:** whether the intent is "remove the atmospheric cast from where
there is signal" (in which case the label is wrong) or "neutralise the sky" (in
which case the weighting is wrong). Both are defensible; they are different
features. Decide, then fix one of the two.

---

## 3. The corona reads to ~2 R where it is traced to 3.7 R

**Status:** unmeasured. Nico's standing complaint, and the half of "finer detail
and a more neutral background" that 0.22.9 did not address.

The report says the corona is traced out to 3.7 lunar radii before the signal
drops into the sky noise, and the rendered picture stops reading at about 2 R.
That is tone, not detail extraction: `radialFlatten` (0.5) and `baseLift` (0.255)
are doing the work, and neither has ever been measured against how far the corona
is actually detectable.

**What would settle it:** the radial median of the *rendered* composite against
the radial median of `hdr_lum`, per shell, over a sweep of `radialFlatten`. The
question is at what radius the render's contrast per decade of real signal falls
below what the eye resolves. Cheap — cached data only.

---

## 4. Detail balance is a two-layer blend, not a full basis

**Status:** shipped and working (0.22.0), but coarser than it could be.

`detailScale` blends the all-scale MGN with a fine-scales-only layer. That reaches
Hill's ladder (fine/coarse energy 3.18 against his 3.32, indistinguishable side by
side) but only along one axis: any ladder constant within each group.

Storing the six per-scale terms would make **any** ladder a render-time linear
combination, since an MGN layer is a weighted mean of its per-scale terms. Cost is
~260 MB as uint8 per work directory. Probably not worth it unless someone wants a
ladder the current axis cannot reach.

---

## 5. `mgn_fine.npy` cannot be rebuilt without a re-stack

**Status:** known, deliberate, still annoying.

RHEF rebuilds itself on first load when its file is absent. `mgn_fine` cannot,
because it needs the flattened denoised luminance that only exists during the
stack. So a pre-0.22 work directory has the Detail balance slider greyed out
(0.22.1 makes that visible) until the folder is re-run.

It *is* reconstructible from `hdr_lum.npy` — denoise, subtract the Fourier
background, run the three fine scales — but it has to reproduce `norm_span` and
the prominence mask **exactly**, or the two MGN layers stop being commensurable
and the blend is quietly wrong. Worth doing carefully; not worth doing quickly.

---

## 6. Stars are not in the partial-convolution mask

**Status:** correct to defer. No data yet.

Hill's rule is "Moon, proms, stars should be 0 in the mask". 0.21.1 added the
prominences (+20% near-limb coronal structure on the rays that have one, nothing
on the rays that do not). Stars are listed and not masked, because the 2026
brackets have none to mask — so there is nothing to measure against and nothing
to claim. Belongs with the first dataset that actually shows stars.

---

## 7. We call it phase correlation and do not run phase correlation

**Status:** label fixed (0.22.2), behaviour not decided.

Both call sites pass `normalization=None` to `phase_cross_correlation`, which
disables the amplitude normalisation — that is plain cross-correlation. The
report cited Druckmüller 2009 for a method the flag switches off; the report text
now says what we actually do.

Whether to turn it **on** is a separate, measurable question, and worth asking on
Clifton's 560 mm set, where the alignment fails outright. Unnormalised
cross-correlation is defensible in high noise (skimage's own docs say so) — but it
is a choice nobody here has made deliberately.

---

## 8. Speed: the inner corona is a third of the run

**Status:** measured, untouched.

From the reference run's own timing: inner corona raw pass 2m25s and denoised
pass 2m24s out of 13m19s — **36%**, for one layer, because it runs the multiscale
normalisation twice at full resolution. Nothing else comes close.

If speed is ever the goal, that is the only place worth looking. FFT-based
Gaussians instead of `ndimage.gaussian_filter`, or the GPU. Not a correctness
issue, so it stays below everything above it.

---

## 9. A diamond artifact inside the lunar disc in NAFE

**Status:** cosmetic, masked, unexplained.

The NAFE layer shows a clear rotated-square (L1-ball) pattern inside the lunar
disc, visible in both the shipped and the corrected-neighbourhood versions. It
sits well inside the disc mask, so it never reaches the picture — but a diamond is
the signature of a Chebyshev/Manhattan neighbourhood somewhere, and nothing in
`nafe_vn` should produce one. Worth understanding before trusting that function in
a region that *does* reach the picture.

---

## 10. RHEF: no amplitude-preserving variant

**Status:** rejected as-is (0.22.4), possibly revivable.

RHEF is off by default because a rank transform discards amplitude: inside one
annulus a 1% real modulation and a 1% noise fluctuation both stretch to the full
range. The SNR gate does not help — the data's radial coherence is still 0.55 at
3.0 R, so the gate stays open exactly where the grain is.

What might work is a variant that keeps rank *ordering* but rescales to the
annulus's own robust amplitude. That is close to reinventing MGN, which is why it
is at the bottom of this list rather than in the source.
