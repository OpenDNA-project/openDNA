#!/usr/bin/env bash
# Builds openDNA for macOS (dist/openDNA.app, zipped as dist/openDNA-macos.zip) or Linux (the folder dist/openDNA, packed as
# dist/openDNA-linux.tar.gz) from packaging/openDNA.spec, smoke-tests it (scripts/smoke_unix.py) and prints its size and SHA-256.
# The macOS / Linux counterpart of scripts/build_exe.ps1.
#
#   pip install -r requirements-dev.txt   (and Pillow on macOS: PyInstaller makes the .icns icon with it)
#   maturin develop -m core/Cargo.toml --release
#   bash scripts/build_unix.sh [--no-smoke]
#
# Uses $PYTHON (default: python3), which must be the one where opendna_core and PyInstaller are installed.
set -euo pipefail
cd "$(dirname "$0")/.."
PYTHON="${PYTHON:-python3}"
SMOKE=1
[ "${1:-}" = "--no-smoke" ] && SMOKE=0
step() { printf '\n== %s\n' "$*"; }
sha256() { if command -v sha256sum >/dev/null; then sha256sum "$1" | cut -d' ' -f1; else shasum -a 256 "$1" | cut -d' ' -f1; fi; }

case "$(uname -s)" in
    Darwin) OS=macos ;;
    Linux) OS=linux ;;
    *) echo "build_unix.sh is for macOS and Linux; on Windows use scripts/build_exe.ps1" >&2; exit 1 ;;
esac

step "Checking the setup ($OS)"
"$PYTHON" -c "import sys; print('python', sys.version.split()[0], sys.executable)"
"$PYTHON" -c "import PyInstaller; print('PyInstaller', PyInstaller.__version__)"
"$PYTHON" -c "import PySide6; print('PySide6', PySide6.__version__)"
"$PYTHON" -c "import opendna_core; print('opendna_core', opendna_core.__file__)" ||
    { echo "opendna_core is not installed in this Python. Run: maturin develop -m core/Cargo.toml --release" >&2; exit 1; }
[ "$OS" = macos ] && { "$PYTHON" -c "import PIL" || { echo "Pillow is needed for the macOS icon: pip install pillow" >&2; exit 1; }; }

step "Building (PyInstaller, about 1-2 minutes)"
rm -rf dist/openDNA dist/openDNA.app dist/openDNA-macos.zip dist/openDNA-linux.tar.gz
"$PYTHON" -m PyInstaller --noconfirm --clean --log-level WARN --workpath build/pyinstaller --distpath dist packaging/openDNA.spec

if [ "$OS" = macos ]; then
    PROGRAM=dist/openDNA.app/Contents/MacOS/openDNA
    [ -x "$PROGRAM" ] || { echo "PyInstaller finished but $PROGRAM is not there" >&2; exit 1; }
    rm -rf dist/openDNA   # (the folder the .app was made from)
else
    PROGRAM=dist/openDNA/openDNA
    [ -x "$PROGRAM" ] || { echo "PyInstaller finished but $PROGRAM is not there" >&2; exit 1; }
    # Next to the program, where the user sees them: how to run it, the licence, the notices.
    cp packaging/LINUX_README.txt dist/openDNA/README.txt
    cp LICENSE dist/openDNA/LICENSE
    for f in THIRD_PARTY_NOTICES.md TRADEMARKS.md; do [ -f "$f" ] && cp "$f" dist/openDNA/; done
    ls dist/openDNA/_internal/PySide6/Qt/lib/libQt6Core.so.6 >/dev/null 2>&1 || ls dist/openDNA/_internal/libQt6Core.so.6 >/dev/null ||
        { echo "libQt6Core is not a plain file in the build, so Qt could not be replaced" >&2; exit 1; }
fi

if [ "$SMOKE" = 1 ]; then
    step "Smoke test"
    "$PYTHON" scripts/smoke_unix.py "$PROGRAM"
else
    echo "(smoke test skipped)"
fi

# Packed last, after the smoke test, so the archive is of what was tested.
step "Packing"
if [ "$OS" = macos ]; then
    ARCHIVE=dist/openDNA-macos.zip
    ditto -c -k --sequesterRsrc --keepParent dist/openDNA.app "$ARCHIVE"   # keeps the bundle's symlinks and attributes (zip -r does not)
else
    ARCHIVE=dist/openDNA-linux.tar.gz
    tar -C dist -czf "$ARCHIVE" openDNA   # keeps the executable bit, which a zip loses
fi
echo "$ARCHIVE  $(du -h "$ARCHIVE" | cut -f1)"
echo "SHA-256 $(sha256 "$ARCHIVE")"
