"""Build the EclipseForgeHDR desktop app on this machine (macOS, Windows, or
Linux for a test), then start the frozen program once to check it serves.

    python packaging/build_app.py            build, test, package
    python packaging/build_app.py --no-test  skip the start test
    python packaging/build_app.py --lab      a test build from the lab tree: named
                                             "EclipseForgeHDR Lab App", it keeps
                                             the lab's work folder, output folder
                                             and its own port (8768)

Needs: the package installed in this Python (pip install .), plus
pyinstaller and pywebview. The GitHub workflow release-apps.yml runs this on
a Mac and on Windows for every release tag and attaches the results:
  macOS    EclipseForgeHDR-<version>-macOS-<arm64|x86_64>.dmg
  Windows  EclipseForgeHDR-<version>-Windows-Setup.exe (an installer: Start menu
           entry, uninstall in Windows Settings), and the same as a .zip
The apps are not signed: macOS and Windows ask once whether to open them.
"""
import os
import platform
import shutil
import subprocess
import sys
import time
import urllib.request
import zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
LAB = "--lab" in sys.argv
NAME = "EclipseForgeHDR Lab App" if LAB else "EclipseForgeHDR"
DIST = os.path.join(ROOT, "dist")
WORK = os.path.join(ROOT, "build", "app")


def version():
    import eclipseforgehdr
    return eclipseforgehdr.__version__


def freeze():
    import eclipseforgehdr
    icon = os.path.join(os.path.dirname(eclipseforgehdr.__file__), "app_icon.png")
    args = [sys.executable, "-m", "PyInstaller", "--noconfirm", "--clean", "--windowed",
            "--log-level", "WARN", "--name", NAME,
            "--distpath", DIST, "--workpath", WORK, "--specpath", WORK,
            "--collect-data", "eclipseforgehdr",
            "--collect-submodules", "eclipseforgehdr",
            "--collect-all", "rawpy",
            "--collect-all", "webview",
            "--collect-data", "certifi",
            "--collect-all", "imageio_ffmpeg",
            "--collect-submodules", "skimage",
            "--collect-data", "skimage"]
    if os.path.exists(icon):
        args += ["--icon", icon]
    if LAB:                         # the lab switch baked in: a double-clicked app has no environment
        os.makedirs(WORK, exist_ok=True)
        hook = os.path.join(WORK, "lab_hook.py")
        open(hook, "w").write('import os\nos.environ["EFHDR_APP_LAB"] = "1"\n')
        args += ["--runtime-hook", hook]
    if sys.platform == "darwin":
        args += ["--osx-bundle-identifier", "org.eclipseforgehdr.labapp2" if LAB else "org.eclipseforgehdr.app"]
    args.append(os.path.join(HERE, "app_main.py"))
    print("== freezing:", " ".join(args[2:6]), "...", flush=True)
    subprocess.check_call(args)


def app_bundle():
    """The Mac app as named for the user, with its version: "EclipseForgeHDR 0.32.0.app"."""
    return os.path.join(DIST, f"{NAME} {version()}.app")


def name_mac_app():
    """PyInstaller writes NAME.app; it is renamed to carry the version, and the
    version goes into its Info.plist too (Finder's Get Info shows it)."""
    import plistlib
    src, dst = os.path.join(DIST, NAME + ".app"), app_bundle()
    if os.path.exists(dst):
        shutil.rmtree(dst)
    os.rename(src, dst)
    pl = os.path.join(dst, "Contents", "Info.plist")
    with open(pl, "rb") as f:
        info = plistlib.load(f)
    info["CFBundleShortVersionString"] = info["CFBundleVersion"] = version()
    info["CFBundleDisplayName"] = f"{NAME} {version()}"
    with open(pl, "wb") as f:
        plistlib.dump(info, f)


def exe_path():
    if sys.platform == "darwin":
        return os.path.join(app_bundle(), "Contents", "MacOS", NAME)
    if os.name == "nt":
        return os.path.join(DIST, NAME, NAME + ".exe")
    return os.path.join(DIST, NAME, NAME)


def start_test(port=8799, timeout=240):
    """The frozen program, without its window, must answer /api/version."""
    env = dict(os.environ, EFHDR_APP_NOWINDOW="1", EFHDR_APP_PORT=str(port))
    p = subprocess.Popen([exe_path()], env=env)
    try:
        t0 = time.time()
        while time.time() - t0 < timeout:
            try:
                with urllib.request.urlopen(f"http://127.0.0.1:{port}/api/version", timeout=2) as r:
                    body = r.read().decode("utf-8")
                print("== start test ok:", body, flush=True)
                with urllib.request.urlopen(f"http://127.0.0.1:{port}/", timeout=5) as r:
                    assert b"EclipseForgeHDR" in r.read(), "the page did not come back"
                return
            except Exception:
                if p.poll() is not None:
                    raise SystemExit(f"the frozen program ended with code {p.returncode}")
                time.sleep(1)
        raise SystemExit("the frozen program did not answer in time")
    finally:
        p.kill()


def package():
    v = version().replace(" ", "")
    os.makedirs(os.path.join(ROOT, "release"), exist_ok=True)
    if sys.platform == "darwin":
        arch = platform.machine() or "arm64"
        out = os.path.join(ROOT, "release", f"{NAME.replace(' ', '')}-{v}-macOS-{arch}.dmg")
        stage = os.path.join(WORK, "dmg")
        shutil.rmtree(stage, ignore_errors=True)
        os.makedirs(stage)
        subprocess.check_call(["ditto", app_bundle(), os.path.join(stage, os.path.basename(app_bundle()))])
        os.symlink("/Applications", os.path.join(stage, "Applications"))
        if os.path.exists(out):
            os.remove(out)
        subprocess.check_call(["hdiutil", "create", "-volname", f"{NAME} {v}", "-srcfolder", stage,
                               "-ov", "-format", "UDZO", out])
    else:
        tag = "Windows" if os.name == "nt" else "Linux"
        out = os.path.join(ROOT, "release", f"{NAME.replace(' ', '')}-{v}-{tag}.zip")
        src = os.path.join(DIST, NAME)
        with zipfile.ZipFile(out, "w", zipfile.ZIP_DEFLATED) as z:
            for d, _, fs in os.walk(src):
                for f in fs:
                    full = os.path.join(d, f)
                    z.write(full, os.path.join(NAME, os.path.relpath(full, src)))
    print("== packaged:", out, f"({os.path.getsize(out) / 1e6:.0f} MB)", flush=True)
    if os.name == "nt":
        setup = windows_installer(v)
        if setup:
            print("== installer:", setup, f"({os.path.getsize(setup) / 1e6:.0f} MB)", flush=True)
    return out


ISS = r"""
[Setup]
AppId={{6E3C2B0A-9C1E-4F4B-9C7E-EF4D0A1C2B31}
AppName=EclipseForgeHDR
AppVersion=%(v)s
AppPublisher=EclipseForgeHDR
DefaultDirName={autopf}\EclipseForgeHDR
DefaultGroupName=EclipseForgeHDR
DisableProgramGroupPage=yes
PrivilegesRequiredOverridesAllowed=dialog
OutputDir=%(out)s
OutputBaseFilename=EclipseForgeHDR-%(v)s-Windows-Setup
Compression=lzma2/max
SolidCompression=yes
ArchitecturesInstallIn64BitMode=x64compatible
SetupIconFile=%(ico)s
UninstallDisplayIcon={app}\EclipseForgeHDR.exe
WizardStyle=modern

[Tasks]
Name: "desktopicon"; Description: "Create a desktop shortcut"; Flags: unchecked

[Files]
Source: "%(src)s\*"; DestDir: "{app}"; Flags: recursesubdirs createallsubdirs ignoreversion

[Icons]
Name: "{group}\EclipseForgeHDR"; Filename: "{app}\EclipseForgeHDR.exe"
Name: "{group}\Uninstall EclipseForgeHDR"; Filename: "{uninstallexe}"
Name: "{autodesktop}\EclipseForgeHDR"; Filename: "{app}\EclipseForgeHDR.exe"; Tasks: desktopicon

[Run]
Filename: "{app}\EclipseForgeHDR.exe"; Description: "Start EclipseForgeHDR"; Flags: nowait postinstall skipifsilent
"""


def windows_installer(v):
    """An Inno Setup installer around the frozen folder, when Inno Setup is
    there (the GitHub Windows machines have it)."""
    iscc = shutil.which("iscc") or next((p for p in (
        r"C:\Program Files (x86)\Inno Setup 6\ISCC.exe", r"C:\Program Files\Inno Setup 6\ISCC.exe")
        if os.path.exists(p)), None)
    if not iscc:
        print("== no Inno Setup: the .zip only", flush=True)
        return None
    ico = os.path.join(WORK, "app.ico")
    try:
        import eclipseforgehdr
        from PIL import Image
        Image.open(os.path.join(os.path.dirname(eclipseforgehdr.__file__), "app_icon.png")).save(
            ico, sizes=[(256, 256), (64, 64), (48, 48), (32, 32), (16, 16)])
    except Exception as e:
        print("== icon not converted:", e, flush=True)
        ico = ""
    iss = os.path.join(WORK, "setup.iss")
    text = ISS % {"v": v, "out": os.path.join(ROOT, "release"), "src": os.path.join(DIST, NAME), "ico": ico}
    if not ico:
        text = text.replace("SetupIconFile=\n", "")
    open(iss, "w", encoding="utf-8").write(text)
    subprocess.check_call([iscc, "/Q", iss])
    out = os.path.join(ROOT, "release", f"EclipseForgeHDR-{v}-Windows-Setup.exe")
    return out if os.path.exists(out) else None


if __name__ == "__main__":
    freeze()
    if sys.platform == "darwin":
        name_mac_app()
    if "--no-test" not in sys.argv:
        start_test()
    package()
