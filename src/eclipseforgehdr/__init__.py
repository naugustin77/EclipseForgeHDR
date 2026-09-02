__version__ = "0.16.2"

# Builds whose cached pipeline products are interchangeable with this one's.
# A release that only changes the interface should not cost the user another
# full stack -- 0.14.1 touched progress reporting and the master-flat preview,
# neither of which writes anything the renderer reads differently.
#
# Only add a version here after checking that NOTHING which writes into the
# work directory changed between them: the merge, the layers, geometry.json.
# When in doubt leave it out; a needless re-run costs time, a wrong reuse costs
# a wrong picture.
CACHE_COMPAT = frozenset({"0.14.0", "0.14.1", "0.14.2",
                          "0.15.4", "0.15.5"})
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


def cache_ok(build):
    """Can products written by build `build` be reused by this build?"""
    return build == __version__ or (build in CACHE_COMPAT
                                    and __version__ in CACHE_COMPAT)
