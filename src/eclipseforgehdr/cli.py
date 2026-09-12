"""eclipseforgehdr command-line entry point."""
import argparse
import os


def main():
    from .brand import APP_NAME, DEFAULT_PORT, LAB
    prog = "eclipseforgehdr-lab" if LAB else "eclipseforgehdr"
    ap = argparse.ArgumentParser(
        prog=prog,
        description=f"{APP_NAME} — High-Dynamic-Range Solar Eclipse Image Processing")
    ap.add_argument("folder", nargs="?", default=None,
                    help="folder with the bracketed raw files (can also be chosen in the GUI)")
    ap.add_argument("--port", type=int, default=DEFAULT_PORT)
    ap.add_argument("--no-browser", action="store_true",
                    help="don't open the browser automatically")
    from . import __version__
    ap.add_argument("--version", action="version",
                    version=f"{prog} {__version__}")
    args = ap.parse_args()
    from .server import main as serve
    serve(folder=args.folder, port=args.port, open_browser=not args.no_browser)


def main_lab():
    """Entry point for the side-by-side test build.

    Sets the flag BEFORE anything imports `brand`, which is what keeps the two
    installs off each other's work directory, port and output folder. Doing it
    here rather than in the package means a normal install is untouched: the
    released build has no code path that can turn this on by accident.
    """
    os.environ["ECLIPSEFORGE_LAB"] = "1"
    main()


if __name__ == "__main__":
    main()
