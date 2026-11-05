<#
Smoke test of the built exe: start it the way a user would (with a file on the command line) on a stripped-down environment, check that
it shows its window with the file loaded, stays alive, shows no error dialog, bundles its data files and loads nothing from the
build machine.

  powershell -File scripts/smoke_exe.ps1 [-Exe dist\OpenGenEdit.exe] [-File demo\demo_vector.gb] [-Seconds 5]

Exit code 0 = passed. It is what scripts/build_exe.ps1 runs after the build (and so what the release workflow runs).
A one-file exe starts TWO processes (a small launcher that unpacks the bundle into %TEMP%\_MEIxxxx, and the program itself), so the
checks look at both. A windowed exe has no console, so a crash cannot be read from stderr: it shows up as a dead process, or as an
error dialog (PyInstaller's "Unhandled exception" box), and both are checked.
#>
param(
    [string]$Exe = "dist\OpenGenEdit.exe",
    [string]$File = "demo\demo_vector.gb",
    [int]$Seconds = 5
)
$ErrorActionPreference = "Stop"
Set-Location (Split-Path $PSScriptRoot -Parent)
$Exe = (Resolve-Path $Exe).Path
$File = (Resolve-Path $File).Path
$procName = [IO.Path]::GetFileNameWithoutExtension($Exe)

# What must be inside the bundle (the same list as DATA in packaging/OpenGenEdit.spec).
$expected = @("parts\ncbi_vectors.tsv", "parts\engineering_proteins.tsv", "tools.toml", "sponsors.json", "vendor\3Dmol-min.js", "vendor\3Dmol.LICENSE", "LICENSE", "logo\icon_256.png", "logo\wordmark.png", "logo\wordmark_dark.png")
if (Test-Path THIRD_PARTY_NOTICES.md) { $expected += "THIRD_PARTY_NOTICES.md" }   # bundled by the spec when it exists
if (Test-Path TRADEMARKS.md) { $expected += "TRADEMARKS.md" }

if (-not ("WinList" -as [type])) {
    Add-Type -TypeDefinition @"
using System;
using System.Collections.Generic;
using System.Runtime.InteropServices;
using System.Text;
public static class WinList {
    delegate bool EnumProc(IntPtr h, IntPtr l);
    [DllImport("user32.dll")] static extern bool EnumWindows(EnumProc p, IntPtr l);
    [DllImport("user32.dll")] static extern uint GetWindowThreadProcessId(IntPtr h, out uint pid);
    [DllImport("user32.dll", CharSet = CharSet.Unicode)] static extern int GetWindowText(IntPtr h, StringBuilder s, int n);
    [DllImport("user32.dll")] static extern bool IsWindowVisible(IntPtr h);
    // Titles of the visible top-level windows that belong to one of the given process ids.
    public static List<string> Titles(uint[] pids) {
        var ids = new HashSet<uint>(pids);
        var found = new List<string>();
        EnumWindows(delegate (IntPtr h, IntPtr l) {
            uint pid; GetWindowThreadProcessId(h, out pid);
            if (ids.Contains(pid) && IsWindowVisible(h)) {
                var sb = new StringBuilder(512); GetWindowText(h, sb, 512);
                if (sb.Length > 0) found.Add(sb.ToString());
            }
            return true;
        }, IntPtr.Zero);
        return found;
    }
}
"@
}

function Exe-Pids { @(Get-Process -Name $procName -ErrorAction SilentlyContinue | ForEach-Object { [uint32]$_.Id }) }
function Window-Titles { @([WinList]::Titles([uint32[]](Exe-Pids))) }

$failures = New-Object System.Collections.Generic.List[string]
function Fail($msg) { $script:failures.Add($msg); Write-Host "  FAIL: $msg" -ForegroundColor Red }
function Ok($msg) { Write-Host "  ok:   $msg" }

# A new folder to run from, a PATH with Windows only (no Python, no conda, no Rust), and a private data folder, so the test neither
# helps itself to the dev machine nor touches the user's primer database.
$work = Join-Path ([IO.Path]::GetTempPath()) ("ge-smoke-" + [guid]::NewGuid().ToString("N").Substring(0, 8))
New-Item -ItemType Directory -Path $work | Out-Null
$savedPath = $env:PATH
$savedData = $env:GENOMEEDIT_DATA_DIR
$env:PATH = "$env:SystemRoot\System32;$env:SystemRoot"
Remove-Item Env:PYTHONPATH, Env:PYTHONHOME, Env:VIRTUAL_ENV, Env:CONDA_PREFIX -ErrorAction SilentlyContinue
$env:GENOMEEDIT_DATA_DIR = Join-Path $work "data"

Write-Host "Smoke test: $Exe"
Write-Host "  opening $File, running from $work"
$launcher = $null
try {
    $launcher = Start-Process -FilePath $Exe -ArgumentList ('"' + $File + '"') -WorkingDirectory $work -PassThru
    $fileName = [IO.Path]::GetFileName($File)
    $title = $null
    $deadline = (Get-Date).AddSeconds(60)
    # Wait for the main window to show the file's name (the title reads "file.gb (N bp) - OpenGenEdit").
    while ((Get-Date) -lt $deadline) {
        Start-Sleep -Milliseconds 500
        if ($launcher.HasExited) { break }
        $title = Window-Titles | Where-Object { $_ -like "*$fileName*" } | Select-Object -First 1
        if ($title) { break }
    }
    if ($launcher.HasExited) {
        Fail "the exe exited by itself with code $($launcher.ExitCode) before showing the file"
    } elseif (-not $title) {
        Fail "no window with the file's name appeared within 60 s (windows seen: $((Window-Titles) -join ' | '))"
    } else {
        Ok "window shown: $title"
        # Stay alive a few seconds: a late crash (the sponsors page, a timer, the first paint) would show now.
        Start-Sleep -Seconds $Seconds
        if ($launcher.HasExited) {
            Fail "the exe died within $Seconds s of starting (exit code $($launcher.ExitCode))"
        } else {
            Ok "still running after $Seconds s"
        }
        $titles = Window-Titles
        $bad = $titles | Where-Object { $_ -match "exception|traceback|error|failed|fatal" }
        if ($bad) { Fail "an error window is showing: $($bad -join ' | ')" } else { Ok "no error dialog (windows: $($titles -join ' | '))" }

        # The bundle: the unpack folder is where python*.dll was loaded from.
        $mei = $null
        $stray = New-Object System.Collections.Generic.List[string]
        $injected = New-Object System.Collections.Generic.List[string]
        $win = $env:SystemRoot.ToLower()
        foreach ($p in @(Get-Process -Name $procName -ErrorAction SilentlyContinue)) {
            try { $mods = @($p.Modules) } catch { continue }
            $py = $mods | Where-Object { $_.ModuleName -like "python3*.dll" } | Select-Object -First 1
            if ($py) {
                $mei = Split-Path $py.FileName -Parent
                foreach ($m in $mods) {
                    $f = $m.FileName.ToLower()
                    if (-not $f.StartsWith($mei.ToLower()) -and -not $f.StartsWith($win) -and $f -ne $Exe.ToLower()) {
                        # Other software injects its own DLLs into every window (display managers, overlays): only a DLL from a
                        # developer tool means the exe leans on the build machine.
                        if ($f -match "conda|python|site-packages|\\\.cargo|rustup|\\msys|\\mingw") { $stray.Add($m.FileName) }
                        else { $injected.Add($m.FileName) }
                    }
                }
            }
        }
        if (-not $mei) {
            Fail "could not find the unpacked bundle (no python3*.dll loaded)"
        } else {
            Ok "bundle unpacked to $mei"
            $lost = @($expected | Where-Object { -not (Test-Path (Join-Path $mei $_)) })
            if ($lost.Count) { Fail "data file(s) missing from the bundle: $($lost -join ', ')" } else { Ok "all $($expected.Count) data files are in the bundle" }
            # The folder build exists so that Qt can be replaced: its Qt DLL must be a plain file next to the exe, not packed.
            if ((Test-Path (Join-Path (Split-Path $Exe -Parent) "_internal")) -and -not (Test-Path (Join-Path $mei "PySide6\Qt6Core.dll"))) {
                Fail "folder build: _internal\PySide6\Qt6Core.dll is not a plain file"
            }
            if (-not (Get-ChildItem $mei -Recurse -Filter "genomeedit_core*.pyd" -ErrorAction SilentlyContinue)) {
                Fail "genomeedit_core*.pyd is not in the bundle"
            } else { Ok "Rust module genomeedit_core is in the bundle" }
            if ($stray.Count) { Fail "loads DLLs from a developer tool outside the bundle: $($stray -join '; ')" } else { Ok "no DLL loaded from Python, conda, Rust or MSYS" }
            if ($injected.Count) { Write-Host "  note: other software injected DLLs into the window (ignored): $(($injected | Split-Path -Leaf) -join ', ')" }
        }
    }
} finally {
    $env:PATH = $savedPath
    if ($savedData) { $env:GENOMEEDIT_DATA_DIR = $savedData } else { Remove-Item Env:GENOMEEDIT_DATA_DIR -ErrorAction SilentlyContinue }
    if ($launcher -and -not $launcher.HasExited) { & taskkill /PID $launcher.Id /T /F | Out-Null }
    Get-Process -Name $procName -ErrorAction SilentlyContinue | Where-Object { $_.Path -eq $Exe } | Stop-Process -Force -ErrorAction SilentlyContinue
    Start-Sleep -Milliseconds 800
    Remove-Item $work -Recurse -Force -ErrorAction SilentlyContinue
}
if ($failures.Count) { Write-Host "Smoke test FAILED ($($failures.Count) problem(s))" -ForegroundColor Red; exit 1 }
Write-Host "Smoke test passed" -ForegroundColor Green
