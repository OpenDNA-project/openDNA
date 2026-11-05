<#
Builds dist\OpenGenEdit.exe (one file) and dist\OpenGenEdit-folder.zip (the same program as a folder, where the Qt files can be replaced:
see docs/PUBLISHING.md, "Qt (LGPL)") from packaging\OpenGenEdit.spec, smoke-tests both and prints their size and SHA-256.

  pip install -r requirements-dev.txt
  maturin develop -m core/Cargo.toml --release
  powershell -File scripts/build_exe.ps1 [-NoSmoke] [-NoFolder] [-Python path\to\python.exe]

Uses the Python on PATH (or -Python, or $env:PYTHON); it must be the one where genomeedit_core and PyInstaller are installed.
Run it from anywhere: it works from the repository root. The exe is not bit-for-bit reproducible (PyInstaller stamps the build time).
#>
param(
    [switch]$NoSmoke,
    [switch]$NoFolder,
    [string]$Python = $(if ($env:PYTHON) { $env:PYTHON } else { "python" })
)
$ErrorActionPreference = "Stop"
$root = Split-Path $PSScriptRoot -Parent
Set-Location $root

function Step($msg) { Write-Host "`n== $msg" -ForegroundColor Cyan }
function Run-Python { & $Python @args; if ($LASTEXITCODE -ne 0) { throw "python $($args -join ' ') failed (exit $LASTEXITCODE)" } }

Step "Checking the setup"
Run-Python -c "import sys; print('python', sys.version.split()[0], sys.executable)"
Run-Python -c "import PyInstaller; print('PyInstaller', PyInstaller.__version__)"
Run-Python -c "import PySide6; print('PySide6', PySide6.__version__)"
try { Run-Python -c "import genomeedit_core; print('genomeedit_core', genomeedit_core.__file__)" }
catch { throw "genomeedit_core is not installed in this Python. Run: maturin develop -m core/Cargo.toml --release" }
foreach ($f in "packaging\OpenGenEdit.spec", "packaging\OpenGenEdit.ico", "packaging\FOLDER_README.txt", "app\main.py", "app\parts\ncbi_vectors.tsv", "app\parts\engineering_proteins.tsv", "app\vendor\3Dmol-min.js") {
    if (-not (Test-Path $f)) { throw "missing $f (is this a complete checkout?)" }
}
if (-not (Test-Path THIRD_PARTY_NOTICES.md)) { Write-Host "note: THIRD_PARTY_NOTICES.md not found, so it is not bundled" -ForegroundColor Yellow }

Step "Building dist\OpenGenEdit.exe (PyInstaller, about 1-2 minutes)"
if (Test-Path dist\OpenGenEdit.exe) { Remove-Item dist\OpenGenEdit.exe -Force }
# Warnings only; the full log is build\pyinstaller\OpenGenEdit\warn-OpenGenEdit.txt.
Run-Python -m PyInstaller --noconfirm --clean --log-level WARN --workpath build\pyinstaller --distpath dist packaging\OpenGenEdit.spec
if (-not (Test-Path dist\OpenGenEdit.exe)) { throw "PyInstaller finished but dist\OpenGenEdit.exe is not there" }

$exe = Get-Item dist\OpenGenEdit.exe
$hash = (Get-FileHash $exe.FullName -Algorithm SHA256).Hash.ToLower()
Step "Built"
"{0}  {1:N1} MB ({2:N0} bytes)" -f $exe.FullName, ($exe.Length / 1MB), $exe.Length
"SHA-256 $hash"
$vi = $exe.VersionInfo
"File properties: $($vi.ProductName) $($vi.FileVersion) | $($vi.CompanyName) | $($vi.LegalCopyright)"

if (-not $NoFolder) {
    Step "Building dist\OpenGenEdit-folder (the same recipe, folder layout)"
    if (Test-Path dist\OpenGenEdit-folder) { Remove-Item dist\OpenGenEdit-folder -Recurse -Force }
    if (Test-Path dist\OpenGenEdit-folder.zip) { Remove-Item dist\OpenGenEdit-folder.zip -Force }
    $env:GENOMEEDIT_BUILD = "folder"
    try { Run-Python -m PyInstaller --noconfirm --clean --log-level WARN --workpath build\pyinstaller-folder --distpath dist packaging\OpenGenEdit.spec }
    finally { Remove-Item Env:GENOMEEDIT_BUILD -ErrorAction SilentlyContinue }
    if (-not (Test-Path dist\OpenGenEdit-folder\OpenGenEdit.exe)) { throw "PyInstaller finished but dist\OpenGenEdit-folder\OpenGenEdit.exe is not there" }
    # Next to the exe, where the user sees them: how to swap Qt, the licence, the notices.
    Copy-Item packaging\FOLDER_README.txt dist\OpenGenEdit-folder\README.txt
    Copy-Item LICENSE dist\OpenGenEdit-folder\LICENSE
    if (Test-Path THIRD_PARTY_NOTICES.md) { Copy-Item THIRD_PARTY_NOTICES.md dist\OpenGenEdit-folder\THIRD_PARTY_NOTICES.md }
    if (Test-Path TRADEMARKS.md) { Copy-Item TRADEMARKS.md dist\OpenGenEdit-folder\TRADEMARKS.md }
    if (-not (Get-ChildItem dist\OpenGenEdit-folder -Recurse -Filter "Qt6Core.dll")) { throw "Qt6Core.dll is not a plain file in the folder build, so Qt could not be replaced" }
}

if ($NoSmoke) { Write-Host "`n(smoke test skipped)" }
else {
    Step "Smoke test (one file)"
    & powershell -NoProfile -ExecutionPolicy Bypass -File "$PSScriptRoot\smoke_exe.ps1" -Exe $exe.FullName
    if ($LASTEXITCODE -ne 0) { throw "smoke test failed" }
    if (-not $NoFolder) {
        Step "Smoke test (folder)"
        & powershell -NoProfile -ExecutionPolicy Bypass -File "$PSScriptRoot\smoke_exe.ps1" -Exe (Resolve-Path dist\OpenGenEdit-folder\OpenGenEdit.exe).Path
        if ($LASTEXITCODE -ne 0) { throw "smoke test of the folder build failed" }
    }
}

if (-not $NoFolder) {
    # Zipped last, after the smoke test: running the program writes nothing into the folder, but a zip must be of what was tested.
    # Python's zipfile writes the forward-slash paths every unzip tool expects (Compress-Archive of Windows PowerShell 5.1 does not).
    Step "Zipping dist\OpenGenEdit-folder.zip"
    Run-Python -c "import shutil; shutil.make_archive('dist/OpenGenEdit-folder', 'zip', 'dist', 'OpenGenEdit-folder')"
    $zip = Get-Item dist\OpenGenEdit-folder.zip
    "{0}  {1:N1} MB" -f $zip.FullName, ($zip.Length / 1MB)
    "SHA-256 $((Get-FileHash $zip.FullName -Algorithm SHA256).Hash.ToLower())"
}
