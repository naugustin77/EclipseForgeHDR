"""eclipseforgehdr command-line entry point."""
import argparse


def main():
    ap = argparse.ArgumentParser(
        prog="eclipseforgehdr",
        description="EclipseForgeHDR — High-Dynamic-Range Solar Eclipse Image Processing")
    ap.add_argument("folder", nargs="?", default=None,
                    help="folder with the bracketed raw files (can also be chosen in the GUI)")
    ap.add_argument("--port", type=int, default=8765)
    ap.add_argument("--no-browser", action="store_true",
                    help="don't open the browser automatically")
    from . import __version__
    ap.add_argument("--version", action="version",
                    version=f"eclipseforgehdr {__version__}")
    args = ap.parse_args()
    from .server import main as serve
    serve(folder=args.folder, port=args.port, open_browser=not args.no_browser)


if __name__ == "__main__":
    main()
