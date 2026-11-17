; Inno Setup script for the Windows installer (dist\openDNA-setup.exe), built by scripts\build_exe.ps1 from the folder build
; dist\openDNA-folder (packaging\openDNA.spec with OPENDNA_BUILD=folder). By hand:
;   ISCC /DAppVersion=0.7.0 packaging\openDNA.iss
;
; Per-user install (no administrator rights, nothing outside the user's own profile and registry hive), so the program can
; replace its own files when it updates itself: app\updater.py runs this installer silently with /RELAUNCH=1. The presence of
; unins000.exe beside openDNA.exe is how the program knows it is an installed copy (updater.installed_copy). User data is in
; %APPDATA%\openDNA, not here: an update or an uninstall never touches it.
#ifndef AppVersion
  #error Pass the version: ISCC /DAppVersion=0.7.0 packaging\openDNA.iss
#endif

[Setup]
AppId={{8E0B7C4A-3F52-4D6B-9A1E-5C2D7F90B3A6}
AppName=openDNA
AppVersion={#AppVersion}
AppPublisher=openDNA project
AppPublisherURL=https://github.com/OpenDNA-project/openDNA
DefaultDirName={autopf}\openDNA
DisableProgramGroupPage=yes
PrivilegesRequired=lowest
ArchitecturesAllowed=x64compatible
ArchitecturesInstallIn64BitMode=x64compatible
MinVersion=10.0
SetupIconFile=openDNA.ico
UninstallDisplayIcon={app}\openDNA.exe
OutputDir=..\dist
OutputBaseFilename=openDNA-setup
Compression=lzma2
SolidCompression=yes
WizardStyle=modern
CloseApplications=yes
RestartApplications=no
VersionInfoVersion={#AppVersion}

[Tasks]
Name: "assoc"; Description: "Offer openDNA for sequence files (.gb, .fasta, ...) in ""Open with"""; GroupDescription: "File types:"

[InstallDelete]
; an update must not leave files the new build no longer ships (a stale Qt file beside the new ones)
Type: filesandordirs; Name: "{app}\_internal"

[Files]
Source: "..\dist\openDNA-folder\*"; DestDir: "{app}"; Flags: recursesubdirs ignoreversion

[Icons]
Name: "{autoprograms}\openDNA"; Filename: "{app}\openDNA.exe"

[Run]
Filename: "{app}\openDNA.exe"; Parameters: "--register-file-types"; Flags: runhidden; Tasks: assoc
Filename: "{app}\openDNA.exe"; Description: "Start openDNA"; Flags: nowait postinstall skipifsilent
; the updater's silent run: start the new version once the files are in place
Filename: "{app}\openDNA.exe"; Flags: nowait runasoriginaluser; Check: Relaunching

[UninstallRun]
Filename: "{app}\openDNA.exe"; Parameters: "--unregister-file-types"; Flags: runhidden; RunOnceId: "UnregisterFileTypes"

[Code]
function Relaunching: Boolean;
begin
  Result := ExpandConstant('{param:RELAUNCH|0}') = '1';
end;
