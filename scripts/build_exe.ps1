<#
Builds the Windows program from packaging\openDNA.spec as a folder (dist\openDNA-folder, where the Qt files can be replaced: see
docs/PUBLISHING.md, "Qt (LGPL)"), smoke-tests it, packs it two ways and smoke-tests the installer:
  dist\openDNA-setup.exe     the installer (packaging\openDNA.iss, Inno Setup): the main download, and what updates itself
  dist\openDNA-folder.zip    the same program as a portable folder; it does not update itself

  pip install -r requirements-dev.txt
  maturin develop -m core/Cargo.toml --release
  powershell -File scripts/build_exe.ps1 [-NoSmoke] [-NoInstaller] [-Python path\to\python.exe]

Needs Inno Setup 6 (ISCC.exe on PATH or in its usual folder: winget install JRSoftware.InnoSetup) unless -NoInstaller.
Uses the Python on PATH (or -Python, or $env:PYTHON); it must be the one where opendna_core and PyInstaller are installed.
Run it from anywhere: it works from the repository root. The build is not bit-for-bit reproducible (PyInstaller stamps the build time).
#>
param(
    [switch]$NoSmoke,
    [switch]$NoInstaller,
    [string]$Python = $(if ($env:PYTHON) { $env:PYTHON } else { "python" })
)
$ErrorActionPreference = "Stop"
$root = Split-Path $PSScriptRoot -Parent
Set-Location $root

function Step($msg) { Write-Host "`n== $msg" -ForegroundColor Cyan }
function Run-Python { & $Python @args; if ($LASTEXITCODE -ne 0) { throw "python $($args -join ' ') failed (exit $LASTEXITCODE)" } }
function Run-Script($script, $scriptArgs) {
    & powershell -NoProfile -ExecutionPolicy Bypass -File "$PSScriptRoot\$script" @scriptArgs
    if ($LASTEXITCODE -ne 0) { throw "$script failed" }
}

Step "Checking the setup"
Run-Python -c "import sys; print('python', sys.version.split()[0], sys.executable)"
Run-Python -c "import PyInstaller; print('PyInstaller', PyInstaller.__version__)"
Run-Python -c "import PySide6; print('PySide6', PySide6.__version__)"
try { Run-Python -c "import opendna_core; print('opendna_core', opendna_core.__file__)" }
catch { throw "opendna_core is not installed in this Python. Run: maturin develop -m core/Cargo.toml --release" }
foreach ($f in "packaging\openDNA.spec", "packaging\openDNA.ico", "packaging\openDNA.iss", "packaging\FOLDER_README.txt", "app\main.py", "app\parts\ncbi_vectors.tsv", "app\parts\engineering_proteins.tsv", "app\vendor\3Dmol-min.js") {
    if (-not (Test-Path $f)) { throw "missing $f (is this a complete checkout?)" }
}
if (-not (Test-Path THIRD_PARTY_NOTICES.md)) { Write-Host "note: THIRD_PARTY_NOTICES.md not found, so it is not bundled" -ForegroundColor Yellow }
$iscc = $null
if (-not $NoInstaller) {
    $iscc = (Get-Command ISCC.exe -ErrorAction SilentlyContinue).Source
    if (-not $iscc) {
        $iscc = @("$env:LOCALAPPDATA\Programs\Inno Setup 6\ISCC.exe", "${env:ProgramFiles(x86)}\Inno Setup 6\ISCC.exe", "$env:ProgramFiles\Inno Setup 6\ISCC.exe") |
            Where-Object { Test-Path $_ } | Select-Object -First 1
    }
    if (-not $iscc) { throw "Inno Setup 6 (ISCC.exe) was not found. Install it (winget install JRSoftware.InnoSetup) or pass -NoInstaller." }
    Write-Host "Inno Setup: $iscc"
}

Step "Building dist\openDNA-folder (PyInstaller, about 1-2 minutes)"
if (Test-Path dist\openDNA-folder) { Remove-Item dist\openDNA-folder -Recurse -Force }
foreach ($f in "dist\openDNA-folder.zip", "dist\openDNA-setup.exe") { if (Test-Path $f) { Remove-Item $f -Force } }
$env:OPENDNA_BUILD = "folder"
# Warnings only; the full log is build\pyinstaller\openDNA\warn-openDNA.txt.
try { Run-Python -m PyInstaller --noconfirm --clean --log-level WARN --workpath build\pyinstaller --distpath dist packaging\openDNA.spec }
finally { Remove-Item Env:OPENDNA_BUILD -ErrorAction SilentlyContinue }
$exe = Get-Item dist\openDNA-folder\openDNA.exe -ErrorAction SilentlyContinue
if (-not $exe) { throw "PyInstaller finished but dist\openDNA-folder\openDNA.exe is not there" }
# Next to the exe, where the user sees them: how to swap Qt, the licence, the notices.
Copy-Item packaging\FOLDER_README.txt dist\openDNA-folder\README.txt
Copy-Item LICENSE dist\openDNA-folder\LICENSE
if (Test-Path THIRD_PARTY_NOTICES.md) { Copy-Item THIRD_PARTY_NOTICES.md dist\openDNA-folder\THIRD_PARTY_NOTICES.md }
if (Test-Path TRADEMARKS.md) { Copy-Item TRADEMARKS.md dist\openDNA-folder\TRADEMARKS.md }
if (-not (Get-ChildItem dist\openDNA-folder -Recurse -Filter "Qt6Core.dll")) { throw "Qt6Core.dll is not a plain file in the folder build, so Qt could not be replaced" }
$vi = $exe.VersionInfo
Step "Built"
"{0}  {1:N1} MB" -f $exe.FullName, ($exe.Length / 1MB)
"File properties: $($vi.ProductName) $($vi.FileVersion) | $($vi.CompanyName) | $($vi.LegalCopyright)"

if ($NoSmoke) { Write-Host "`n(smoke test skipped)" }
else {
    Step "Smoke test (folder)"
    Run-Script "smoke_exe.ps1" @("-Exe", $exe.FullName)
}

if ($iscc) {
    Step "Building dist\openDNA-setup.exe (Inno Setup)"
    $version = (Select-String -Path app\project.py -Pattern '^VERSION\s*=\s*"([^"]+)"').Matches[0].Groups[1].Value
    & $iscc /Qp "/DAppVersion=$version" packaging\openDNA.iss
    if ($LASTEXITCODE -ne 0) { throw "ISCC failed (exit $LASTEXITCODE)" }
    $setup = Get-Item dist\openDNA-setup.exe
    "{0}  {1:N1} MB" -f $setup.FullName, ($setup.Length / 1MB)
    if (-not $NoSmoke) {
        Step "Smoke test (installer: install, start, uninstall)"
        Run-Script "smoke_installer.ps1" @("-Setup", $setup.FullName)
    }
}

# Zipped last, after the smoke test: running the program writes nothing into the folder, but a zip must be of what was tested.
# Python's zipfile writes the forward-slash paths every unzip tool expects (Compress-Archive of Windows PowerShell 5.1 does not).
Step "Zipping dist\openDNA-folder.zip"
Run-Python -c "import shutil; shutil.make_archive('dist/openDNA-folder', 'zip', 'dist', 'openDNA-folder')"
$zip = Get-Item dist\openDNA-folder.zip
"{0}  {1:N1} MB" -f $zip.FullName, ($zip.Length / 1MB)
foreach ($f in $zip, $setup) {
    if ($f) { "SHA-256 $((Get-FileHash $f.FullName -Algorithm SHA256).Hash.ToLower())  $($f.Name)" }
}
