# Changelog — EclipseForgeHDR

Newest first. Entries from 0.6.1 onward were written at the time. The 0.7.2 –
0.9.8 block was reconstructed afterwards from the code comments that name their
own version and from the development record; where a change cannot be pinned to
an exact version it is filed under the release it is known to precede.

## 0.10.0 — full code audit

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
