# -*- mode: python ; coding: utf-8 -*-
"""PyInstaller recipe for dist/openDNA.exe (one file, no console window) and, with OPENDNA_BUILD=folder, dist/openDNA-folder/.

Build it with  powershell scripts/build_exe.ps1  (that script checks the setup and smoke-tests the result).
By hand, from the repository root:  python -m PyInstaller --noconfirm --clean --workpath build/pyinstaller --distpath dist packaging/openDNA.spec

Why a committed spec: the data files the program reads at run time are listed here, so the bundle never depends on what happened
to be lying around. The program finds them with Path(__file__).with_name(...) (flat modules in the bundle root), so each file is
bundled next to where the module sits: app/parts/ -> parts/, app/vendor/3Dmol-min.js -> vendor/3Dmol-min.js, and so on.
If you add a file the app reads (a new .tsv/.json/.toml, an image), add it to DATA below, otherwise the exe fails on that feature
only (not at startup), which is the classic frozen-app bug. scripts/smoke_exe.ps1 checks that every file listed here is in the bundle.

Name and version of the Windows file properties come from app/project.py (NAME, VERSION, AUTHOR), so a release has one place to bump.
"""
import importlib.util
import os
from pathlib import Path

from PyInstaller.utils.win32.versioninfo import (
    FixedFileInfo, StringFileInfo, StringStruct, StringTable, VarFileInfo, VarStruct, VSVersionInfo,
)

ROOT = Path(SPECPATH).resolve().parent          # SPECPATH (set by PyInstaller) is this file's folder: packaging/
APP = ROOT / "app"

# app/project.py is plain constants with no imports: load it by path (the spec runs before `app` is on sys.path).
_spec = importlib.util.spec_from_file_location("ge_project", APP / "project.py")
project = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(project)

# (source file, folder inside the bundle). Required files: the build stops when one is missing.
DATA = [
    (APP / "parts", "parts"),                  # parts_library.py: the bundled collections (one file per source)
    (APP / "tools.toml", "."),                 # runner.py: recipes for external tools
    (APP / "sponsors.json", "."),              # sponsors.py: startup page
    (ROOT / "demo" / "pBAD30_eilR.gb", "demo"),   # sponsors.py: DEMO, "Open the Demo Plasmid" on the startup page
    (APP / "vendor" / "3Dmol-min.js", "vendor"),   # structure.py: VENDOR, the protein viewer's library
    (APP / "vendor" / "3Dmol.LICENSE", "vendor"),  # BSD-3-Clause: must travel with the library
    (ROOT / "LICENSE", "."),                   # AGPL-3.0
]
# main.py: window icon (all sizes); sponsors.py: the wordmark on the welcome page. Drawn by packaging/make_icon.py.
DATA += [(APP / "logo" / f"icon_{size}.png", "logo") for size in (16, 24, 32, 48, 64, 128, 256)]
DATA += [(APP / "logo" / "wordmark.png", "logo"), (APP / "logo" / "wordmark_dark.png", "logo")]
# Optional files: bundled when present (a notice file written later must not break an older checkout).
for name in ("THIRD_PARTY_NOTICES.md", "TRADEMARKS.md"):   # (Help > About names both)
    if (ROOT / name).exists():
        DATA.append((ROOT / name, "."))
missing = [str(src) for src, _ in DATA if not src.exists()]
if missing:
    raise SystemExit("openDNA.spec: missing data file(s): " + ", ".join(missing))
datas = [(str(src), dest) for src, dest in DATA]

# Windows "Details" tab of the exe: file version 0.1.0.0 from VERSION "0.1.0".
parts = [int(x) for x in project.VERSION.split(".")[:3]]
vtuple = tuple(parts + [0] * (4 - len(parts)))
version_info = VSVersionInfo(
    ffi=FixedFileInfo(filevers=vtuple, prodvers=vtuple, mask=0x3F, flags=0x0, OS=0x40004, fileType=0x1, subtype=0x0),
    kids=[
        StringFileInfo([StringTable("040904B0", [   # US English, Unicode
            StringStruct("CompanyName", project.AUTHOR),
            StringStruct("FileDescription", f"{project.NAME} - offline DNA editor"),
            StringStruct("FileVersion", project.VERSION + ".0"),
            StringStruct("InternalName", project.NAME),
            StringStruct("LegalCopyright", f"Copyright (C) 2026 {project.AUTHOR}. Licensed under AGPL-3.0."),
            StringStruct("OriginalFilename", f"{project.NAME}.exe"),
            StringStruct("ProductName", project.NAME),
            StringStruct("ProductVersion", project.VERSION),
            StringStruct("Comments", project.REPO_URL),
        ])]),
        VarFileInfo([VarStruct("Translation", [0x0409, 1200])]),
    ],
)

a = Analysis(
    [str(APP / "main.py")],
    pathex=[str(APP)],                          # the app's modules import each other by flat name (from editor import ...)
    binaries=[],
    datas=datas,
    hiddenimports=[],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=["tkinter", "pytest", "Bio", "numpy"],   # never used by the app; keeps stray installs of them out of the bundle
    noarchive=False,
    optimize=0,
)
pyz = PYZ(a.pure)

# Two layouts from this one recipe. Default: one file, dist/openDNA.exe. OPENDNA_BUILD=folder: dist/openDNA-folder/ (the exe
# plus an _internal folder where the Qt DLLs are plain files). The folder exists so that a user can replace Qt, which the LGPL asks
# for and a one-file exe cannot offer (docs/PUBLISHING.md, "Qt (LGPL)"). scripts/build_exe.ps1 runs PyInstaller once for each.
FOLDER = os.environ.get("OPENDNA_BUILD") == "folder"

exe = EXE(
    pyz,
    a.scripts,
    [] if FOLDER else a.binaries,
    [] if FOLDER else a.datas,
    [],
    exclude_binaries=FOLDER,
    name=project.NAME,
    icon=str(ROOT / "packaging" / "openDNA.ico"),
    version=version_info,
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=False,          # UPX-packed exes are flagged by antivirus far more often, and the saving is small
    upx_exclude=[],
    runtime_tmpdir=None,
    console=False,      # a window program: no black console box
    disable_windowed_traceback=False,   # a crash shows a dialog instead of vanishing
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
)

if FOLDER:
    coll = COLLECT(
        exe,
        a.binaries,
        a.datas,
        strip=False,
        upx=False,
        upx_exclude=[],
        name=project.NAME + "-folder",
    )
