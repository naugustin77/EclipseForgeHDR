# Contributing

This is a small, single-author project. Issues and pull requests are welcome;
please read this first, because the constraints here are unusual.

## The thing that matters most: measure it

Almost every parameter in this codebase was arrived at by measuring, not by
taste, and the reasoning is written down next to the code — often with the
numbers that settled it. Several changes that *looked* obviously right made the
image measurably worse, and are documented as such so nobody repeats them.

So: if you change a constant or an algorithm, say what you measured. The run
report already prints alignment residuals, per-tier variance at the limb and in
the corona, the limb transition width, lunar drift and the photometric factors.
"Before / after" on those numbers, from the same dataset, is the currency here.
A change justified only by a screenshot is very hard to accept, because a
corona image can be made to look better while being less true.

## Datasets

The pipeline has been developed against two sets: a good one (Panasonic
DC-S1R II, 600 mm f/8, 49 frames, 14 tiers, 12.6 EV, lunar radius ~620 px) and a
weak one (Sony DSC-HX99, 118 mm, low SNR). That is not many. **The most useful
contribution is a test on a third dataset**, especially one that differs in the
ways the code is least sure about:

- a lunar radius much smaller than ~600 px
- an untracked or alt-az mount, where field rotation is present (the aligner
  currently solves translation only, so this is expected to fail — but knowing
  *how* it fails is valuable)
- fewer than 8 tiers, or an EV span under 8
- a camera whose sensor is not RGGB-phase Bayer, or with unusual black/white
  levels

Please include the run report (`eclipseforge_output/..._report.txt`) with any
issue. It records everything the pipeline decided and why.

## Scope

Things that fit: alignment, merging, the enhancement layers, colour, robustness
across cameras, and the honest reporting of what the pipeline did.

Things that do not: turning this into a general astro stacker. There are good
ones already, and the whole reason this exists is that the corona is not a star
field.

## Style

- Comments explain *why*, and cite the measurement or the paper. There are a lot
  of them and that is deliberate.
- No new dependencies without a strong reason.
- Published methods are cited in the code where they are implemented; if you
  deviate from a paper's constants, say so and say why.

## Licensing of contributions

The project is source-available and not yet licensed (see
[LICENSE](LICENSE)). By opening a pull request you agree that your contribution
may be included and released under whatever licence the project eventually
adopts, alongside the rest of the code. If that is not acceptable to you, please
open an issue describing the change instead — a good bug report is worth as much
as a patch here.

## Attribution

Most of the code was written by Claude (Anthropic) in collaboration with the
author, iterating against real data. That is stated so nobody is surprised by
the density of the comments or the occasional confession of a bug in them.
