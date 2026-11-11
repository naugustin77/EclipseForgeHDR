"""Frame folders beside the bracket: the partial phases.

A subfolder of partial-phase frames is listed in the Frames panel and left
out of the stack; it can also be set by hand in tab 1 ("off" for none).
"""
import os

PHASES_WORDS = ("phase", "partial", "sequence", "crescent")


def find_phases_dir(folder):
    """A subfolder of partial-phase frames, or None (a folder set by hand in
    tab 1 first; "off" means none)."""
    try:
        from .frameselect import dir_override
        st, p = dir_override(folder, "phases")
        if st == "off":
            return None
        if st == "path":
            return p
    except Exception:
        pass
    try:
        for name in sorted(os.listdir(folder)):
            p = os.path.join(folder, name)
            n = name.lower()
            if os.path.isdir(p) and any(w in n for w in PHASES_WORDS) and not n.startswith("eclipseforge"):
                return p
    except Exception:
        pass
    return None
