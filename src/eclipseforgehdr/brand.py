"""Everything that has to differ between the released app and a test build.

WHY THIS EXISTS. A second install of this app next to a working one is not a
rename. Three things collide if only the package name changes, and every one of
them costs the user something they cannot get back cheaply:

  * THE WORK DIRECTORY. Both builds would write `.eclipseforgehdr/` inside the
    user's raw folder. The cache key does not know which build filled it, so a
    test run silently invalidates a finished release run, and the next release
    run re-stacks from the raws. On a 45 MP twelve-tier bracket that is a
    quarter of an hour, and it happens without a word.
  * THE PORT. Both bind 127.0.0.1:8765, so the second one to start either
    fails or attaches to the first one's server and reports the wrong version.
  * THE OUTPUT FOLDER. Exports and diagnostics bundles from the two builds land
    on each other under the same names, and the user cannot tell them apart
    afterwards.

So a lab build gets its own of each. It is driven by an environment variable
rather than a build flag because the value has to be visible to the pipeline,
the server and the GUI alike, and because a stale hardcoded copy in one of the
three is exactly the bug this file exists to prevent -- there is one definition
and everything reads it.

`ECLIPSEFORGE_LAB=1` is set by the lab build's own console script (see
`cli.main_lab`), so the user never types it. Setting it by hand on a normal
install gives the same separation, which is a legitimate way to try a risky
setting against a folder whose release cache you want to keep.
"""
import os

LAB = os.environ.get("ECLIPSEFORGE_LAB") == "1"

#: Shown in the GUI title, the header and the credit line.
APP_NAME = "EclipseForgeHDR Lab" if LAB else "EclipseForgeHDR"

#: Per-folder cache directory, inside the user's raw folder.
WORKDIR_NAME = ".eclipseforgehdr-lab" if LAB else ".eclipseforgehdr"

#: Where exports and the diagnostics bundle are written, beside the raws.
OUTPUT_NAME = "eclipseforge_lab_output" if LAB else "eclipseforge_output"

#: Default loopback port for the GUI server.
DEFAULT_PORT = 8766 if LAB else 8765
