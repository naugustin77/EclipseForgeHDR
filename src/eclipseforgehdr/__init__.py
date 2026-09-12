__version__ = "0.23.3"

# 0.23.3 changes the SHARED PEDESTAL for any file whose header reports no black
# level (every FITS, normally): fits.py stops reading OFFSET as one, and the fit
# may then range to 0.02 of saturation instead of 0.002. That moves the merged
# levels, so hdr_lum.npy and every layer built from it are different data and a
# cached workdir from 0.23.2 or earlier is NOT reusable. Deliberately not listed
# below. See TODO 22.
# Builds whose cached pipeline products are interchangeable with this one's.
# A release that only changes the interface should not cost the user another
# full stack -- 0.14.1 touched progress reporting and the master-flat preview,
# neither of which writes anything the renderer reads differently.
#
# Only add a version here after checking that NOTHING which writes into the
# work directory changed between them: the merge, the layers, geometry.json.
# When in doubt leave it out; a needless re-run costs time, a wrong reuse costs
# a wrong picture.
# 0.22.5 changes the MERGE WEIGHT, so hdr_lum.npy and every layer built from it
# are different data. A cached workdir from any earlier build is not reusable
# and a full re-stack is correct -- exactly the case the "a wrong reuse costs a
# wrong picture" rule exists for. It opens the third family below.
# 0.22.0 adds mgn_fine.npy and CANNOT rebuild it from a finished work directory
# -- it needs the flattened, denoised luminance that only exists during the
# stack. An older workdir is still listed here because it loads and renders
# exactly as it did: `detailScale` is guarded on the layer's presence and is
# inert without it. So this is the 0.17.0 situation ("loads" is not "has the
# feature") with the sharp edge removed -- nobody gets a wrong picture, they
# just do not get the slider until they re-run.
# 0.21.0 adds a cached product (rhef.npy), and would normally be excluded on
# that ground alone. It is listed anyway because the renderer BUILDS the layer
# on first load when the file is absent -- it is one sort of the luminance
# already on disk, seconds rather than the sixteen minutes a re-stack costs --
# and then caches it. So an older work directory gains the feature rather than
# merely tolerating its absence, which is the distinction that kept 0.17.0 out.
# 0.20.0 and 0.19.0 change render.py and gui.html only -- which layer drives the
# prominence contrast term. Nothing that writes into the work directory moved,
# so a 0.18.0 stack is reused as it stands and nobody pays for another run.
# 0.18.0 writes a new cached product (promdet.npy), so a 0.17.0 workdir is
# missing it. The key is then ABSENT from the layer dict rather than filled
# with a flat 0.5 -- 0.5 is not the identity here, the blend pulls towards the
# layer -- so an old cache still LOADS; it just has nothing to blend.
# 0.17.0 is deliberately not listed, because "loads" and "has the feature"
# differ.
# 0.16.0 and 0.16.1 are deliberately NOT listed, and were removed from this set
# in 0.16.2: their photometric factors came from an estimator that has since
# been reverted, so their merged tiers are wrong. Reusing them would show the
# bad merge under the corrected build.
# 0.16.0 changes the photometric factors, so every merged tier and every
# layer built on them differs. Deliberately not listed: reusing a 0.15.x
# merge here would show the old merge under the new report.
# 0.14.3 changes the defect map and the limb override, both of which write into
# the merge -- so it is deliberately NOT listed above: its products differ.
# 0.15.4 likewise: it changes the disc-mask margin an import writes into
# geometry.json, which every detail layer is built against.


# One flat set cannot express this. 0.22.5 changed the MERGE, so anything it
# writes is different data from anything before it -- but the moment 0.22.5+
# were added to the same set as the older builds, every pre-merge stack became
# reusable again and the fix would have been silently withheld from anyone with
# an old work directory. Families, not one bag.
CACHE_FAMILIES = (
    frozenset({"0.14.0", "0.14.1", "0.14.2", "0.15.4", "0.15.5"}),
    frozenset({"0.18.0", "0.19.0", "0.20.0", "0.20.1", "0.20.2", "0.20.3",
               "0.21.0", "0.21.1", "0.21.2", "0.21.3", "0.21.4",
               "0.22.0", "0.22.1", "0.22.2", "0.22.3", "0.22.4"}),
    # 0.22.5 changed the merge weight; .6-.9 add guards, the diagnostics bundle
    # and render defaults, none of which is cached data.
    frozenset({"0.22.5", "0.22.6", "0.22.7", "0.22.8", "0.22.9",
               "0.22.10", "0.22.11", "0.22.12", "0.22.13"}),
    # 0.22.14 moves the fitted lunar radius onto the tiers' consensus when the
    # merged half-level fit runs large. R is what MGN's radial profile, FNRGF's
    # rings, the deband and the disc mask are all built on, so every cached
    # layer from an affected run is different data. Its own family.
    frozenset({"0.22.14"}),
    # 0.22.15 subtracts a fitted pedestal from every tier before the merge, so
    # hdr_lum and every layer built on it are different data -- most of all in
    # the outer field, which is the whole point. Its own family.
    frozenset({"0.22.15"}),
    # 0.22.16 stops the merge weight leaking into each tier's clipped region.
    # That changes hdr_lum everywhere just outside the limb -- the whole point --
    # so nothing earlier is reusable. Its own family.
    # 0.22.17 and .18 change DIAGNOSTICS only -- how the tier-agreement
    # statistic is reported, and a linearity check on the input -- and write
    # nothing the renderer reads differently, so a 0.22.16 work directory is
    # reused as it stands.
    frozenset({"0.22.16", "0.22.17", "0.22.18"}),
    # 0.22.19 changes the merge weight again -- 0.22.16's hard mask left a step
    # at every tier's clipping contour, which printed one arc per tier. Cached
    # products from .16-.18 carry those arcs. Its own family.
    # 0.22.20 adds two environment switches and changes nothing by default, so
    # a 0.22.19 work directory is reused as it stands. A run made WITH a switch
    # set is deliberately not distinguishable here -- that is the point of a
    # bisect, and such a directory should be re-run before it is trusted.
    # 0.22.21 adds a third switch and changes nothing by default.
    # 0.22.22 adds a startup log line only.
    # 0.22.23 changes render defaults only -- nothing cached moves.
    # 0.22.24 adds a measurement and a file in the diagnostics bundle; it
    # changes nothing the renderer reads, so a .19+ work directory is reused.
    frozenset({"0.22.19", "0.22.20", "0.22.21", "0.22.22", "0.22.23",
               "0.22.24"}),
    # 0.22.25 changes the DEFAULT merge weight back to the plain feather, so
    # hdr_lum and every layer built on it are different data from .19-.24 --
    # deliberately so. Its own family. A run made with ECLIPSEFORGE_FEATHER set
    # is not distinguishable here, which is the point of a bisect switch; such
    # a directory should be re-run before it is trusted.
    frozenset({"0.22.25"}),
    # 0.22.26 moves the clipping test onto the mosaic, so more pixels are
    # correctly marked invalid and every tier's contribution to the merge
    # changes. hdr_lum and every layer differ. Its own family.
    frozenset({"0.22.26"}),
    # 0.22.27 changes the cross-exposure alignment filter, so every tier's
    # shift can differ and everything downstream of the stack is different
    # data. Its own family.
    frozenset({"0.22.27"}),
    # 0.22.28 chooses the merge feather per dataset, so hdr_lum can differ from
    # a .27 run on the same folder. Its own family.
    # 0.22.29 fixes a leak that let one folder's feather choice apply to the
    # next folder in the same session, so a .28 work directory may have been
    # merged with the wrong weight. Its own family.
    frozenset({"0.22.28"}),
    frozenset({"0.22.29"}),
    # 0.22.30 fixes the feather trial, which was measuring nothing and so
    # returned "plain" on every dataset. hdr_lum differs wherever the trial now
    # picks the other weight. Its own family.
    # 0.22.31 makes the feather a setting rather than an automatic choice, so
    # a work directory carries whichever weight was asked for; opts.json now
    # records it and the cache check compares it. Its own family.
    frozenset({"0.22.30"}),
    # 0.22.32 adds a Clear cache button and puts the merge weight into the
    # cache key (it was compared on the import path only). Nothing that writes
    # into the work directory moved, and opts.json already recorded the weight,
    # so a 0.22.31 directory is reused as it stands -- correctly now, which is
    # the fix.
    frozenset({"0.22.31", "0.22.32"}),
    # 0.22.33 fixes how the azimuthal per-tier correction is APPLIED: 0.22.26
    # looked k and q up by 60-segment index, which made a smooth function into
    # a staircase and printed a radial edge every 6 degrees, and it applied the
    # fit at full strength far outside the radii it was fitted on. Every tier's
    # contribution to the merge changes, so hdr_lum and every layer differ.
    # Its own family.
    frozenset({"0.22.33"}),
    # 0.22.34 corrects the SATURATION LEVEL: a body that reports the container
    # ceiling rather than its linearity limit had the clipping threshold set
    # 7% too high, so the flat tops of blown highlights were merged as valid
    # data. Every tier's valid mask changes and with it hdr_lum and every
    # layer. Its own family, and no earlier cache is trustworthy.
    # 0.22.34's fix never fired: it read camera_white_level_per_channel after
    # rawpy's `with` block had closed the file, so the attribute raised and the
    # guard fell through to the container ceiling every time. 0.22.35 reads it
    # inside the block. Products from .34 are byte-identical to .33's and both
    # carry the wrong valid mask, so neither is reusable here.
    frozenset({"0.22.34"}),
    # 0.22.36 fills the pixels that are clipped in every tier from the shortest
    # tier instead of leaving them at zero. hdr_rgb and every layer built on it
    # differ wherever such a pixel exists, so .35's products are not reusable.
    frozenset({"0.22.35"}),
    # 0.22.37 changes a log string and adds a stats key. Nothing that writes
    # into the work directory moves, so a 0.22.36 stack is reused as it stands.
    # 0.22.38 changes gui.html only. Nothing cached moves.
    # 0.22.39 adds a per-channel photometry MEASUREMENT and corrects a comment.
    # It writes one new block into the report and nothing else; the merge, the
    # layers and geometry.json are byte-for-byte what 0.22.36 produced.
    # 0.22.40 fixes a crash in that measurement on an odd-sized mosaic. Still
    # nothing cached moves.
    # 0.22.57 CHANGES WHAT THE SKY SUBTRACTION SUBTRACTS -- a quadratic surface
    # where .52-.56 zeroed the shape entirely on any dataset with a real sky
    # gradient, and subtracted a flat level from a tilted sky. A work directory
    # made by any of those with `Subtract sky` on carries that error baked into
    # hdr_rgb and every layer. The cache key cannot tell the two apart -- the
    # option is the same in both -- so this gets its OWN family and everyone
    # re-stacks. Exactly the "a wrong reuse costs a wrong picture" case.
    # 0.22.56 is render-only: the pedestal comes off the colour path and the
    # fade to neutral moves after the gains. No re-stack.
    # 0.22.55 is render-only: the chroma fallback goes back to neutral.
    # 0.22.54 adds Hill's log envelope and sends the `bg` layer at 16 bits.
    # Render and transport only -- nothing cached moves.
    # 0.22.53 fixes the black point on a sky-subtracted stack -- a render-side
    # change, so a 0.22.52 work directory is reused as it stands and gains the
    # fix. It also stops the subtraction reaching the occulted disc, which DOES
    # change hdr_rgb; that part only applies to a stack made after it, and a
    # .52 directory renders correctly without it because the new black point no
    # longer looks at percentiles at all.
    # 0.22.52 adds the SKY SUBTRACTION, which does change the merged image and
    # every layer built on it -- but only when it is switched on, and the switch
    # is part of the cache key (opts.json "remove_sky"), so a work directory
    # made without it is reused exactly as before and one made with it can never
    # be confused for one made without. That is what lets it stay in this family
    # rather than opening another.
    # 0.22.51 separates the sky's additive airlight from the corona. It is
    # measured at LOAD time from the merged HDR already on disk and applied in
    # the render, so a work directory from any build in this family gains the
    # control without a re-stack -- the 0.21.0 situation, not the 0.18.0 one.
    frozenset({"0.22.36", "0.22.37", "0.22.38", "0.22.39", "0.22.40", "0.22.41", "0.22.42", "0.22.43", "0.22.44", "0.22.45", "0.22.46", "0.22.47", "0.22.48", "0.22.49", "0.22.50", "0.22.51", "0.22.52", "0.22.53", "0.22.54", "0.22.55", "0.22.56"}),
    # 0.22.58: in 0.22.57 the subtraction never ran at all (it threw on a
    # read-only memory map and the run silently fell back to the old gradient
    # division), so a .57 work directory holds a DIFFERENT image again. Its own
    # family for the same reason .57 has one.
    frozenset({"0.22.57"}),
    # 0.22.59 fits ONE sky shape with a per-channel amplitude instead of three
    # independent surfaces, which changes what is subtracted. Own family.
    frozenset({"0.22.58"}),
    frozenset({"0.22.59"}),
    # 0.22.60 drops the quadratic back to a plane, so what is subtracted differs
    # from .59 again. Own family.
    frozenset({"0.22.60"}),
    # 0.22.61 REMOVES the sky subtraction and the colour machinery built around
    # it. The merge goes back to what 0.22.48 produced, which is different data
    # from anything .52-.60 wrote. Own family.
    # 0.22.62 restores 0.22.48's corona white balance -- render only, nothing
    # cached moves, so a 0.22.61 work directory is reused as it stands.
    frozenset({"0.22.61", "0.22.62"}),
    # 0.22.63 changes the DEFAULT white balance to camera as-shot and adds VNG
    # as a second demosaic. Both scale or rebuild the merged data, so 0.22.63
    # stands alone: every earlier work directory has to be restacked. It is
    # also its own family so that switching either setting invalidates the
    # cache through opts.json rather than through the version.
    # 0.22.64 ADDS the Hill unsharp-mask set to what a stack produces and
    # changes nothing that was already cached, so a 0.22.63 work directory is
    # reused as it stands -- it simply has no Hill layers and the controls say
    # so. Re-stack to build them.
    # 0.22.65 only adds a way to BUILD the Hill masks from an already-stacked
    # work directory, so it reuses 0.22.63 and 0.22.64 exactly as they are.
    # 0.22.66 is render-only -- the base stretch and the level match both act
    # at render time, so a 0.22.65 work directory (Hill masks included) is
    # reused exactly as it stands and nothing needs rebuilding.
    # 0.22.67 changes the Hill MASKS, but they carry their own recipe version
    # (detail.HILL_BUILD) and the server rebuilds them on their own when it is
    # stale -- so the expensive stack is still reused and nothing else moves.
    frozenset({"0.22.63", "0.22.64", "0.22.65", "0.22.66", "0.22.67",
               "0.22.68"}),
    # 0.22.69 changes the INTERPOLATION the merge resamples every tier with, so
    # the merged data itself is different and there is nothing to reuse.
    # 0.22.70 only takes VNG off the toolbar and makes the run report say which
    # demosaic it used -- the default path is bit-identical to 0.22.69.
    # 0.22.71 changes the Hill masks only; they carry their own recipe version
    # (detail.HILL_BUILD) and are rebuilt on their own, so the stack is reused.
    # 0.22.72 is render-only: the chroma fade is built when the layers are
    # loaded, so nothing on disk moves and nothing is rebuilt.
    frozenset({"0.22.69", "0.22.70", "0.22.71", "0.22.72"}),
    # 0.22.73 changes the prominence mask, which MGN and the Hill masks are
    # both built with, so the layers themselves differ. Its own family.
    # 0.22.74 adds a per-pixel noise map beside the Hill masks and a monochrome
    # Hill view. The map is written by build_hill, which carries its own recipe
    # version (detail.HILL_BUILD), so the masks are rebuilt on their own and
    # the expensive stack is reused; the view is render-only.
    # 0.22.75 renames things on screen and in the run log only -- no product
    # on disk changes, so 0.22.74's cache is reused as it stands.
    frozenset({"0.22.73", "0.22.74", "0.22.75"}),
    # 0.22.76 fixes the frame the pedestal, the per-tier radial check and the
    # LDIC azimuthal fit sample each tier in -- see _ring_sample. All three feed
    # the merge, so hdr_rgb itself differs and there is nothing to reuse. Its
    # own family, and a full re-stack.
    # 0.22.77 adds a diagnostic and an escape hatch; the default merge is
    # 0.22.76's, bit for bit, so its products are reused as they stand.
    # 0.22.78 is input handling, the server, the render fallbacks and report
    # text. Nothing it changes can alter a merged image, so 0.22.76's products
    # are reused as they stand -- with one exception the user has to opt into:
    # a folder holding both raws and a stray TIFF no longer counts that TIFF as
    # a frame, and such a folder's input fingerprint therefore changes, which
    # the fingerprint check catches on its own and re-stacks.
    # 0.22.79 adds tests and one recorded stat; no computation changes.
    frozenset({"0.22.76", "0.22.77", "0.22.78", "0.22.79"}),
    # 0.22.80 and 0.22.81 each changed the photometric chain, so each wrote a
    # different merge and neither is reusable by anything. No family.
    # 0.22.82 reverted that chain; 0.22.83 adds ECLIPSEFORGE_NO_SKYGRAD, which
    # is off unless it is set, so with it unset the two write identical
    # products. 0.22.84 touches the partial-convolution MASKS (their own recipe
    # version, detail.HILL_BUILD, rebuilds them on their own in seconds), the
    # renderer and the server's own cache check -- nothing that can alter a
    # merged image. So a 0.22.82 or 0.22.83 work directory is reused as it
    # stands.
    #
    # WITHOUT THIS ENTRY A VERSION IN NO FAMILY RE-STACKS, which is the safe
    # default and was also simply wrong here: 0.22.84 asked for a 14-minute
    # re-stack of the 600 mm bracket to rebuild masks that take two minutes,
    # and an A/B of one mask setting would have cost half an hour of it.
    # 0.22.85 is 0.22.84 renumbered. Two different file sets went out as
    # 0.22.84 -- the first without the family entry below, so it re-stacked --
    # and a version string that does not identify the code is the one thing
    # this project cannot afford: it is how a switch gets set on a shell that
    # launches a build which ignores it. The products are identical, so a
    # 0.22.84 work directory is reused as it stands.
    # 0.22.86 makes the per-channel linear fit the DEFAULT and puts it on
    # the toolbar. That changes the merge -- but only for a folder that was
    # stacked WITHOUT it, and opts.json already carries the mode, so the
    # per-folder check catches exactly those and leaves the rest alone. A
    # folder stacked with the marker file compares equal and keeps its
    # stack, which is the whole reason the stored value stayed "hill".
    frozenset({"0.22.82", "0.22.83", "0.22.84", "0.22.85", "0.22.86"}),
    # 0.22.87 uncouples the exposure-exponent trial from the moon-MASK
    # verdict. On any bracket where the mask was rejected -- the routine
    # outcome -- the trial had been bailing out, so the exponent stayed at
    # 1.0 where it should have been fitted. Where it now fires, the merge
    # weights differ and hdr_lum is different data. Its own family.
    # 0.23.0 IS 0.22.87 with a version number. Nothing that writes into a
    # work directory changed, so a 0.22.87 stack is reused as it stands --
    # a release should not cost every user a re-stack for a number.
    # 0.23.1 adds settings files: render, server and page only. Nothing that
    # writes into a work directory moved, so a 0.23.0 stack is reused.
    #
    # 0.23.2 IS NOT IN THIS SET, deliberately. It changes the alignment
    # estimator (semi-phase correlation is now the default), the dark rate map
    # (gated to photosites where something was actually measured) and the FNRGF
    # layer (prominence mask and the published noise term). Every one of those
    # changes what is written into the work directory -- where each tier lands,
    # what each tier IS, and what the detail layer contains -- so a cached stack
    # from 0.23.1 is different data and reusing it would show the old alignment
    # under the new report. A full re-stack is correct here, and it is exactly
    # the case the "a wrong reuse costs a wrong picture" rule exists for.
    frozenset({"0.22.87", "0.23.0", "0.23.1"}),

)


def cache_ok(build):
    """Can products written by build `build` be reused by this build?"""
    if build == __version__:
        return True
    return any(build in f and __version__ in f for f in CACHE_FAMILIES)
