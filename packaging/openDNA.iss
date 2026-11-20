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
; Unchecked, and offered only when WSL2 already runs (see WslReady): a silent update (app\updater.py) never runs it. WSL itself is never
; installed here: that takes an administrator and a restart, so the user does it (the page says how).
Name: "wsltools"; Description: "Set up pLannotate and the Terminal's tools in WSL2 now (a console window shows the progress)"; GroupDescription: "Optional tools (several hundred MB, a few minutes; you can also do it later from Detect Features or with oDNA-cli wsltools install):"; Flags: unchecked; Check: WslReady

[InstallDelete]
; an update must not leave files the new build no longer ships (a stale Qt file beside the new ones)
Type: filesandordirs; Name: "{app}\_internal"

[Files]
Source: "..\dist\openDNA-folder\*"; DestDir: "{app}"; Flags: recursesubdirs ignoreversion

[Icons]
Name: "{autoprograms}\openDNA"; Filename: "{app}\openDNA.exe"

[Run]
Filename: "{app}\openDNA.exe"; Parameters: "--register-file-types"; Flags: runhidden; Tasks: assoc
; the cores page's answer, saved in settings.ini (not on the updater's silent run: the user's value stays)
Filename: "{app}\openDNA.exe"; Parameters: "--set-cores {code:ChosenCores}"; Flags: runhidden; Check: not WizardSilent
Filename: "{app}\openDNA.exe"; Description: "Start openDNA"; Flags: nowait postinstall skipifsilent
; the optional tools: oDNA-cli.exe is a console program, so its text (conda's downloads, pLannotate's databases) shows in a terminal window
Filename: "{app}\oDNA-cli.exe"; Parameters: "wsltools install"; StatusMsg: "Setting up pLannotate and the Terminal's tools in WSL2 (a few minutes)..."; Flags: waituntilterminated skipifdoesntexist; Tasks: wsltools
; the updater's silent run: start the new version once the files are in place
Filename: "{app}\openDNA.exe"; Flags: nowait runasoriginaluser; Check: Relaunching

[UninstallRun]
Filename: "{app}\openDNA.exe"; Parameters: "--unregister-file-types"; Flags: runhidden; RunOnceId: "UnregisterFileTypes"

[Code]
var
  WslKnown, WslOk: Boolean;
  WslNote: TNewStaticText;
  CoresPage: TInputQueryWizardPage;

{ CPU cores of this computer (Windows sets NUMBER_OF_PROCESSORS). }
function Cpus: Integer;
begin
  Result := StrToIntDef(GetEnv('NUMBER_OF_PROCESSORS'), 1);
  if Result < 1 then Result := 1;
end;

{ The cores openDNA may use: the saved value (a reinstall), else min(8, cores - 1), at least 1 (app\tool_settings.py DEFAULT_CORES). }
function DefaultCores: Integer;
begin
  Result := Cpus - 1;
  if Result > 8 then Result := 8;
  if Result < 1 then Result := 1;
  Result := StrToIntDef(GetIniString('Performance', 'cores', IntToStr(Result), ExpandConstant('{userappdata}\openDNA\settings.ini')), Result);
end;

function ChosenCores(Param: String): String;
begin
  Result := IntToStr(StrToIntDef(Trim(CoresPage.Values[0]), DefaultCores));
end;

function NextButtonClick(CurPageID: Integer): Boolean;
var
  N: Integer;
begin
  Result := True;
  if CurPageID = CoresPage.ID then
  begin
    N := StrToIntDef(Trim(CoresPage.Values[0]), 0);
    if (N < 1) or (N > Cpus) then
    begin
      MsgBox('Give a whole number from 1 to ' + IntToStr(Cpus) + '.', mbError, MB_OK);
      Result := False;
    end;
  end;
end;

function Relaunching: Boolean;
begin
  Result := ExpandConstant('{param:RELAUNCH|0}') = '1';
end;

{ A WSL distribution is installed and starts (wsl -e true succeeds). Asked once: it can take a few seconds to start WSL. }
function WslReady: Boolean;
var
  Code: Integer;
begin
  if not WslKnown then
  begin
    WslOk := Exec(ExpandConstant('{sys}\wsl.exe'), '-e true', '', SW_HIDE, ewWaitUntilTerminated, Code) and (Code = 0);
    WslKnown := True;
  end;
  Result := WslOk;
end;

{ Without WSL the task is not listed: say why, and what the user does (this installer never installs WSL: an administrator and a restart). }
procedure InitializeWizard;
begin
  WslNote := TNewStaticText.Create(WizardForm);
  WslNote.Parent := WizardForm.SelectTasksPage;
  WslNote.WordWrap := True;
  WslNote.Left := WizardForm.TasksList.Left;
  WslNote.Width := WizardForm.TasksList.Width;
  WslNote.Height := ScaleY(48);
  WslNote.Top := WizardForm.TasksList.Top + WizardForm.TasksList.Height - WslNote.Height;
  WslNote.Caption := 'The optional tools (pLannotate, and the Terminal''s bioinformatics tools) run in WSL2, which is not set up on this computer. ' +
    'To use them, install WSL yourself in an administrator terminal (wsl --install, then restart), and afterwards run oDNA-cli wsltools install, ' +
    'or press Install pLannotate in Detect Features.';
  WslNote.Visible := False;
  CoresPage := CreateInputQueryPage(wpSelectTasks, 'CPU cores',
    'How many of this computer''s ' + IntToStr(Cpus) + ' CPU cores may openDNA use?',
    'Analyses such as Local Align and Folder Search run on several cores at once. Leaving one core free keeps the window and ' +
    'the rest of the computer responsive. You can change it later in Tools > Settings > General > Performance.');
  CoresPage.Add('Cores (1 to ' + IntToStr(Cpus) + '):', False);
  CoresPage.Values[0] := IntToStr(DefaultCores);
end;

procedure CurPageChanged(CurPageID: Integer);
begin
  if (CurPageID = wpSelectTasks) and not WslNote.Visible and not WslReady then   { (going back to the page must not shrink the list again) }
  begin
    WslNote.Visible := True;
    WizardForm.TasksList.Height := WizardForm.TasksList.Height - WslNote.Height - ScaleY(6);
  end;
end;
