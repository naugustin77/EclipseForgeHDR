# PixInsight's denoise tools, and what of them we could use

Notes taken from the TGVDenoise interface and the MultiscaleLinearTransform
reference, supplied by a tester in September 2026. Kept because the useful part is
the STRUCTURE of these tools, not their numbers -- the numbers are in
PixInsight's own units and do not transfer.

Nothing here is implemented yet. This is the reading, with a verdict on each
item, so the next person does not have to do it again.

## What we already do

`detail.denoise_loglum` is an a-trous (starlet) multiscale decomposition with
Donoho soft thresholding, per-scale threshold strengths (`DENOISE_PROFILES`),
and a per-pixel threshold that scales with the photon-noise sigma at that
pixel. That is the same algorithm family MultiscaleLinearTransform implements,
and the per-pixel noise map does the job PixInsight's "linear mask" does --
arguably better, since it comes from a noise model rather than from an
amplified copy of the image.

## MultiscaleLinearTransform

Same decomposition (starlet / a trous, or Gaussian-pyramid variant), same
Donoho soft thresholding, same per-layer thresholds. Three things it has that
we do not:

**Amount below 1.** Their soft threshold is followed by a strength control, and
their own guidance is emphatic about it: "the objective of denoise is not to
completely remove noise, which would create an artificial effect". Our
threshold sets sub-threshold coefficients to zero, i.e. Amount fixed at 1.
This is the cheapest item on the list and the one most likely to matter here,
because "blotchy, artificially looking" is close to the complaint that our
renders read as over-processed.

**Iterations.** They recommend a smaller threshold with Amount below 1 applied
two or three times over one hard single pass, for the same reason.

**Deringing, dark and bright.** For when a layer is BIASED rather than
thresholded -- Gibbs ringing around high-contrast edges. Not relevant to our
denoise, directly relevant to any sharpening that boosts a layer, which is what
the Hill masks do.

Their "linear mask" (amplification 50-200, Gaussian smoothness 1-3 px,
inverted) is the same idea as our photon-noise map and we should not copy it.

## TGVDenoise

Total Generalised Variation (Bredies, Kunisch & Pock) -- second-order TV, which
avoids the staircasing that makes first-order TV unusable on luminance.
Iterative primal-dual, 100 iterations by default. Not in scipy or scikit-image;
`skimage.restoration.denoise_tv_chambolle` is first-order only.

The transferable part is not the algorithm, it is the SPLIT: it denoises
lightness and chrominance separately and hits chrominance harder. Off a real
installation:

    channel        strength   edge protection   smoothness   iterations
    lightness        5.00         0.0020           2.00          100
    chrominance      7.00         0.0030           2.00          100

MultiscaleLinearTransform makes the same split available through its Target
parameter (Lightness / Luminance / Chrominance / RGB-K), and its own docs
recommend exactly this: "a softer noise reduction to the Lightness channel that
contains details and a stronger one to Chrominance".

## The one that matters here

An earlier draft of this file said we do not denoise chroma at all. That is
wrong and the correction is the interesting part. `render.Layers` builds the
chroma ratio as

    rc = gaussian_filter(hdr[:, :, c], 6) / max(Ls, floor)

so each channel IS denoised -- by a fixed 6 px Gaussian -- and then gated by the
confidence fade, which forces the result to exactly neutral below the sky noise
floor. The gap is not the absence of a denoise. It is that ours is a blunt
instrument of fixed width followed by a hard gate.

The width is the part with a measurable cost:

    set       plate scale   sigma = 6 px is   FWHM
    600 mm     1.55 "/px        9.3 "         ~22 "
    250 mm     3.12 "/px       18.7 "         ~44 "

Prominences are typically 20-60 " tall, so the colour of a prominence is being
averaged over a large fraction of the prominence. That is very likely why they
read flat in colour, and why `promChroma` exists at all -- it puts back, in
green and blue, structure the chroma blur had removed.

And the fade is what drew the ring a tester found at White balance None; 0.22.72
smoothstepped it, which softened the symptom and not the cause.

So the proposal both PixInsight tools point at, adjusted for what we actually
do: replace the fixed blur with an EDGE-AWARE chroma denoise, and the fade can
then be relaxed rather than removed. Two gains to expect and to check for --
prominence and inner-corona colour stays put instead of smearing, and the ring
goes at its cause. One honest limit: where the sky has no colour above the
noise there is none to recover, so some fade is still needed out there.

Chroma is piecewise-smooth by nature, so first-order TV -- `denoise_tv_chambolle`,
which IS in scikit-image -- should carry most of it without the staircasing that
rules it out for luminance. The starlet path we already have, run on the chroma
channels with the luminance noise map, is the other candidate.
