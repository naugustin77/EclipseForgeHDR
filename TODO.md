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

## 1d. The render blends layers across fixed radii, and that is ring-shaped

**Status:** identified in the source, not yet tested. Costs a slider move.

Nico, on two sets whose tiers agree to 10-12%: *"some concentric rings in the
composite, likely from Inner."* The render has two blends keyed to absolute
radius, both landing in the 1.0-1.6 R band where the rings are:

```
wI = smoothstep((1.45R - r) / 0.40R)      inner corona: fades OUT, gone at 1.45 R
wf = smoothstep((r - 1.02R) / 0.55R)      FNRGF share: fades IN,  full at 1.57 R
```

`wI` blends a DIFFERENT layer in and out again over 1.05-1.45 R. Wherever the
inner layer and MGN disagree, that blend traces the disagreement as an annulus —
by construction, not by accident.

**The codebase already knows this failure mode.** The glare-dim term used to
share `wI`, and the comment where it was given its own profile says it plainly:
*"its smoothstep window closes at 1.45 R -- and at full strength that is a 3.3x
brightness ramp ending at a definite radius, which prints as a ring."* `wG` was
fixed. `wI` was not.

**The test costs nothing and needs no re-run:** set Short-exposure detail
(`innerMix`) to 0 and look; then FNRGF share (`fnMix`) to 0. Both re-render from
cached layers in seconds. If the rings go with `innerMix`, the fix is to give
the inner layer a profile with no boundary — an exponential decay from the mask
edge, exactly as `wG` was given — rather than a window that shuts at a radius.

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
