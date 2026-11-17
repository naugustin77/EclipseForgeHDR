"""EclipseForgeHDR as a desktop app (Mac .app, Windows .exe).

The same program as the pip install, frozen with PyInstaller, plus this
launcher:

  * starts the app's own image server in a thread on 127.0.0.1;
  * shows the interface in a window of its own (pywebview: the system's
    WebKit on macOS, Edge WebView2 on Windows) instead of a browser tab;
  * when the window closes, the whole program ends, a run in progress with it.

A double-clicked app has no terminal, so everything the server and the
pipeline print goes to a log file:
  macOS    ~/Library/Logs/EclipseForgeHDR.log
  Windows  %LOCALAPPDATA%\\EclipseForgeHDR\\EclipseForgeHDR.log
(the previous start's log is kept as .log.1).

The port is fixed (8765, as the pip install) rather than picked freely because
the page keeps a few conveniences in the browser storage of its address (the
last folders, run durations), and a new port every start would lose them. If
the port is taken by EclipseForgeHDR itself, the app opens a window on that
server; if another program holds it, start the app with EFHDR_APP_PORT=8775.

Build: packaging/build_app.py (run by .github/workflows/release-apps.yml).
"""
import json
import os
import socket
import sys
import threading
import time
import traceback
import urllib.request

# EFHDR_APP_LAB=1: a test build of this app made from the lab tree keeps the
# lab's own work directory, output folder and port, so it never touches the
# caches of a released install (set before anything imports brand)
if os.environ.get("EFHDR_APP_LAB") == "1":
    os.environ["ECLIPSEFORGE_LAB"] = "1"

APP_TITLE = "EclipseForgeHDR"
PORT = int(os.environ.get("EFHDR_APP_PORT", "8768" if os.environ.get("ECLIPSEFORGE_LAB") == "1" else "8765"))
URL = f"http://127.0.0.1:{PORT}/"
HOME = os.path.expanduser("~")
if sys.platform == "darwin":
    LOG = os.path.join(HOME, "Library", "Logs", f"{APP_TITLE}.log")
    STORAGE = os.path.join(HOME, "Library", "Application Support", APP_TITLE)
elif os.name == "nt":
    _base = os.environ.get("LOCALAPPDATA") or os.path.join(HOME, "AppData", "Local")
    LOG = os.path.join(_base, APP_TITLE, f"{APP_TITLE}.log")
    STORAGE = os.path.join(_base, APP_TITLE, "webview")
else:
    LOG = os.path.join(HOME, ".local", "state", APP_TITLE, f"{APP_TITLE}.log")
    STORAGE = os.path.join(HOME, ".local", "share", APP_TITLE, "webview")


def _log_to_file():
    """No terminal: both output streams go to the log."""
    if sys.stdout is not None and getattr(sys.stdout, "isatty", lambda: False)():
        return                                  # started from a terminal: keep it there
    try:
        os.makedirs(os.path.dirname(LOG), exist_ok=True)
        if os.path.exists(LOG):
            os.replace(LOG, LOG + ".1")
        f = open(LOG, "a", buffering=1, encoding="utf-8")
        sys.stdout = sys.stderr = f
    except OSError:
        pass


def _port_open():
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        s.settimeout(0.3)
        return s.connect_ex(("127.0.0.1", PORT)) == 0


def _ours():
    """True when what answers on PORT is an EclipseForgeHDR server."""
    try:
        with urllib.request.urlopen(URL + "api/version", timeout=2) as r:
            v = json.loads(r.read().decode("utf-8"))
        return "EclipseForgeHDR" in str(v.get("app", ""))
    except Exception:
        return False


def _serve():
    try:
        from eclipseforgehdr.server import main as serve
        serve(folder=None, port=PORT, open_browser=False)
    except BaseException:
        traceback.print_exc()
        sys.stdout.flush()


def _wait_for_server(timeout=180.0):
    t0 = time.time()
    while time.time() - t0 < timeout:
        if _port_open():
            return True
        time.sleep(0.2)
    return False


WINDOW_TITLE = APP_TITLE


def _page(msg):
    return ("<html><body style='font:15px system-ui,-apple-system,Segoe UI,sans-serif;"
            "background:#111;color:#ddd;padding:40px'><h2>" + WINDOW_TITLE + "</h2><p>" + msg
            + "</p><p style='color:#888'>Log: " + LOG + "</p></body></html>")


def main():
    _log_to_file()
    from eclipseforgehdr import __version__
    global WINDOW_TITLE
    WINDOW_TITLE = ("EclipseForgeHDR Lab App" if os.environ.get("ECLIPSEFORGE_LAB") == "1" else APP_TITLE) \
        + " " + __version__                     # the version in the window title
    print(f"{time.strftime('%Y-%m-%d %H:%M:%S')} {APP_TITLE} {__version__} starting, port {PORT}", flush=True)

    # EFHDR_APP_NOWINDOW=1: the server only, for the build's own start test
    if os.environ.get("EFHDR_APP_NOWINDOW") == "1":
        _serve()
        return

    import webview
    own_server = True
    if _port_open():
        if _ours():
            own_server = False                  # a second start: another window on the running one
            print("an EclipseForgeHDR server already answers on this port; opening a window on it", flush=True)
        else:
            webview.create_window(WINDOW_TITLE, html=_page(
                f"Port {PORT} is taken by another program. Quit that program, or start "
                f"{APP_TITLE} with the environment variable EFHDR_APP_PORT=8775."),
                width=720, height=320)
            webview.start()
            return
    if own_server:
        threading.Thread(target=_serve, name="efhdr-server", daemon=True).start()

    # the window opens at once with a short note and switches to the interface
    # when the server answers (the first start unpacks numpy, scipy and the
    # rest, which can take some seconds)
    win = webview.create_window(WINDOW_TITLE, html=_page("Starting …"),
                                width=1600, height=1000, min_size=(1000, 700))

    def _load():
        if _wait_for_server():
            win.load_url(URL)
        else:
            win.load_html(_page("The image server did not start. The log says why."))

    try:
        webview.settings["ALLOW_DOWNLOADS"] = True
    except Exception:
        pass
    os.makedirs(STORAGE, exist_ok=True)
    webview.start(_load, private_mode=False, storage_path=STORAGE)
    print(f"{time.strftime('%Y-%m-%d %H:%M:%S')} window closed, quitting", flush=True)
    sys.stdout.flush()
    os._exit(0)                                 # ends the server and any run in progress


if __name__ == "__main__":
    import multiprocessing
    multiprocessing.freeze_support()            # a frozen Windows program that starts processes needs this
    main()
