# The sky subtraction: what was measured, and why it was reverted

Built across 0.22.51–0.22.60 and removed in 0.22.61. The machinery is gone; the
measurements are not, because they are correct and the next attempt should start
from them rather than repeat them.

## What is true

**The atmosphere does two separate things.** It removes corona light on the way
down — multiplicative, uniform over a 3.5 degree field, strongly wavelength
dependent, and the reason a corona shot at 6.8 degrees altitude comes out
orange. And it adds scattered sunlight along the line of sight, which is
additive and has its own colour and gradient.

**The additive part is large and nothing removes it.** On the 600 mm reference
bracket it is about 800 counts per channel, against a corona that is down to 170
counts by 3.5 R and 55 by 5 R. `remove_sky_gradient` divides, which is the wrong
shape for it, and it fits the residual after a radial profile, so only the TILT
ever reaches the image — the level has always stayed in the picture.

**The two-component model is real, and it was tested before it was used.**
Writing each channel as one shared corona shape plus a per-channel sky level
predicts that red plotted against green must be a straight line, whatever the
corona looks like. Measured, radial medians from 1.2 to 6.9 R — a hundredfold
range in signal:

    R = 1.926 * G - 844      median residual 1.4%
    B = 0.244 * G + 344      median residual 0.8%

So the corona's colour does not change with radius, and the entire apparent
outward shift (R/G 2.05 -> 0.91) is the sky taking over. The fitted level agrees
within 3% with the far-field median measured independently.

**Removing it does buy what it promised.** With the level gone the corona's own
range goes from 130:1 to 1900:1, and the outer corona stops sitting on a
pedestal fourteen times larger than itself.

## Why it still failed

Every failure was in the SHAPE of the sky, not in its level or its colour.

* **Fitted per channel, it invents colour.** Each channel's surface is free to
  absorb a different part of what is left in the fitting annulus, and what is
  left there is corona, which is very coloured. The three fitted surfaces came
  out with different shapes (xx term -0.07 / -1.16 / -1.21 in R / G / B), and
  subtracting those is subtracting a colour map: broad pink and blue regions
  over a background whose LUMINANCE measures flat to 15 counts.
* **One shared shape with three amplitudes is better and not enough.** It cannot
  invent a colour pattern, but a quadratic still lets the three amplitudes
  diverge, and their difference is a colour pattern shaped like the quadratic.
* **A plane cannot describe the sky either.** After the best plane fit, the four
  quadrants of the far field sat at R+122/G+86/B+46, R+65/G+37/B+23,
  R+85/G+45/B+23, R+128/G+84/B+41 — a factor of two, on a 23-count pedestal.

The bind: the far field is where the corona is faintest but never absent, so
every term added to the sky model is another way to absorb corona and call it
sky, and every term removed leaves more real gradient behind.

## Things worth carrying forward

* **The occulted disc is not sky.** The merge sets it to exactly zero, so
  subtracting 800 counts from it drives it to -800, and at 3% of the frame it
  owns any low percentile — including the renderer's black point, which then
  added the whole pedestal back and cancelled the subtraction exactly.
* **Do not clip the result at zero.** Rectifying the noise gives each channel a
  mean proportional to its own sigma, which invents colour out of nothing.
* **A grey pedestal is not neutral in effect.** It dilutes the corona's colour
  outward, and a corona white balance then amplifies the dilution. If a pedestal
  is kept for the detail layers' log floor, the colour path must subtract it
  again.
* **Order matters in the colour chain.** Fading chroma toward neutral and then
  applying gains do not commute: the fade raises blue, which sits far below
  neutral in a corona, and the gains multiply the raised value. If both are
  needed, fade last, in log chroma.
* **Luminance checks cannot see this class of error.** Every numeric test passed
  while the picture was plainly wrong. Any future attempt is checked on a
  rendered image.

## The measurement that reframes all of this (found after the revert)

The renderer has faded chroma to neutral wherever the smoothed luminance falls
below the sky's noise since 0.11.4 -- `ratio = 1 + conf*(rc - 1)`. On the
reference bracket that confidence runs 1.000 at 2 R, 0.608 at 4 R and 0.000 at
6 R, and the rendered colour follows it exactly: R/G 1.76, 1.11, 1.00. So the
corona's colour does not fade outward in the picture because the corona changes.
It stops because the fade reaches zero, and the ramp between draws a visible
edge -- reported by the tester, with arrows, as "a clear sharp edge between
orange and more neutral".

Switch the fade off and there is no edge. What appears instead is the sky's own
colour, olive-green, filling the entire frame. The fade has been covering up
exactly the problem the subtraction was built to solve.

**And the sky's colour is very nearly ONE NUMBER.** Measured beyond 5 R on the
un-subtracted data:

    quadrant        level     R/G     B/G
    upper-left      828.0    0.977   0.631
    upper-right     920.2    0.921   0.671
    lower-left      720.6    1.014   0.615
    lower-right     899.8    0.964   0.646

The BRIGHTNESS varies 28% across the frame. The COLOUR varies about 10%, and
part of that 10% is residual corona rather than sky -- the radial series from
5 R to 7.5 R shows R/G falling 1.00 -> 0.90 as the corona thins out, which is
the corona leaving, not the sky changing.

That splits the problem in two, and only one half is hard:

* **The sky's colour** needs ONE gain vector for the whole frame. No surface, no
  per-channel spatial fit, no masks -- and therefore no way to invent a colour
  pattern, which was the failure mode of every version above.
* **The sky's brightness** is a single shared shape across the three channels,
  which is what a luminance-only fit produces and which by construction cannot
  create colour.

Everything that went wrong came from fitting colour and shape together. They
should not have been one fit.

## What a next attempt should probably do differently

Fit the sky where there is no corona at all rather than where it is merely
faint — which on a 3.5 degree field may mean it cannot be fitted from the
science frames, and has to come from somewhere else: a frame taken off the sun,
the partial phases, or a measurement the photographer supplies. A model with
enough freedom to describe a real horizon sky has enough freedom to eat the
outer corona, and no amount of care in fitting escapes that.
