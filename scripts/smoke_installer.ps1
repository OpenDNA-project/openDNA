<#
Smoke test of the Windows installer: install silently into a new folder (the way the updater runs it, plus /DIR), check what the
installer promises, start the installed program (scripts\smoke_exe.ps1), uninstall, check that nothing is left.

  powershell -File scripts/smoke_installer.ps1 [-Setup dist\openDNA-setup.exe]

Exit code 0 = passed. Run by scripts/build_exe.ps1 (and so by the release workflow). It writes the file-type registration into the
current user's registry and removes it again with the uninstall.
#>
param([string]$Setup = "dist\openDNA-setup.exe")
$ErrorActionPreference = "Stop"
Set-Location (Split-Path $PSScriptRoot -Parent)
$Setup = (Resolve-Path $Setup).Path
$dir = Join-Path ([IO.Path]::GetTempPath()) ("ge-install-" + [guid]::NewGuid().ToString("N").Substring(0, 8))
$assoc = "HKCU:\Software\Classes\.gb\OpenWithProgids"
$failures = New-Object System.Collections.Generic.List[string]
function Fail($msg) { $script:failures.Add($msg); Write-Host "  FAIL: $msg" -ForegroundColor Red }
function Ok($msg) { Write-Host "  ok:   $msg" }
function Registered { $null -ne (Get-ItemProperty $assoc -Name "openDNA.SequenceFile" -ErrorAction SilentlyContinue) }

Write-Host "Installer smoke test: $Setup"
try {
    # The updater's own arguments (app/updater.py), minus /RELAUNCH (it would start a window), plus the folder and a log.
    $p = Start-Process -FilePath $Setup -Wait -PassThru -ArgumentList "/VERYSILENT", "/SUPPRESSMSGBOXES", "/NORESTART", "/CLOSEAPPLICATIONS", "/DIR=`"$dir`"", "/LOG=`"$dir.log`""
    if ($p.ExitCode -ne 0) { Fail "the installer exited with code $($p.ExitCode) (log: $dir.log)" }
    elseif (-not (Test-Path "$dir\openDNA.exe")) { Fail "openDNA.exe is not in $dir after the install" }
    else {
        Ok "installed to $dir"
        # unins000.exe beside the exe is how the program knows it is an installed copy (updater.installed_copy)
        if (Test-Path "$dir\unins000.exe") { Ok "uninstaller present" } else { Fail "unins000.exe is missing: the program would not update itself" }
        if (Registered) { Ok "file types registered" } else { Fail "the .gb file type was not registered" }
        & powershell -NoProfile -ExecutionPolicy Bypass -File "$PSScriptRoot\smoke_exe.ps1" -Exe "$dir\openDNA.exe"
        if ($LASTEXITCODE -ne 0) { Fail "the installed program failed its smoke test" } else { Ok "installed program runs" }
        $u = Start-Process -FilePath "$dir\unins000.exe" -Wait -PassThru -ArgumentList "/VERYSILENT", "/SUPPRESSMSGBOXES", "/NORESTART"
        if ($u.ExitCode -ne 0) { Fail "the uninstaller exited with code $($u.ExitCode)" }
        # The uninstaller finishes by a helper process that deletes the last files: give it a moment.
        for ($i = 0; $i -lt 20 -and (Test-Path "$dir\openDNA.exe"); $i++) { Start-Sleep -Milliseconds 500 }
        if (Test-Path "$dir\openDNA.exe") { Fail "openDNA.exe is still there after the uninstall" } else { Ok "uninstalled" }
        if (Registered) { Fail "the file-type registration is still there after the uninstall" } else { Ok "file types unregistered" }
    }
} finally {
    Remove-Item $dir, "$dir.log" -Recurse -Force -ErrorAction SilentlyContinue
}
if ($failures.Count) { Write-Host "Installer smoke test FAILED ($($failures.Count) problem(s))" -ForegroundColor Red; exit 1 }
Write-Host "Installer smoke test passed" -ForegroundColor Green
