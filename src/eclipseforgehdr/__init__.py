__version__ = "0.22.26"

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
)


def cache_ok(build):
    """Can products written by build `build` be reused by this build?"""
    if build == __version__:
        return True
    return any(build in f and __version__ in f for f in CACHE_FAMILIES)
