# Publishing openDNA on GitHub

A checklist for taking the project further. **Done so far:** the repository is public (since 0.7.0) under the AGPL-3.0 ([LICENSE](../LICENSE)), with CI and releases (v0.8.2 to v0.8.4). The sections below are the how-to for each step; read the ones about the past as the record of why.

## 1. Decide first

- **Licence.** Chosen: AGPL-3.0 ([LICENSE](../LICENSE)), kept as the licence of the project. It keeps the source open and
  obliges every fork or hosted copy to stay open under the same licence, so changes can be pulled back; it does not prevent
  forks, and is not meant to. Notes: (1) on a **public** GitHub repository anyone can fork it whatever the licence says.
  (2) To keep one project, the other levers: the **name policy** ([TRADEMARKS.md](../TRADEMARKS.md): forks are welcome under a
  name and an icon of their own; the AGPL's section 7(e) lets the author decline trademark rights), welcoming pull requests,
  and releases that work. A registered trademark on "openDNA" / "oDNA" would back the policy up; none is registered yet.
  (3) Contributions: decided. Contributors assign their copyright to the project owner, who stays the sole copyright holder and
  can change the licence later ([CONTRIBUTING.md](../CONTRIBUTING.md#copyright-of-contributions); the pull-request template asks
  them to tick that they agree). A signed agreement (for example through CLA Assistant) would be firmer than a checkbox.
  (4) The dependencies are compatible: PySide6 is LGPL, `gb-io` is MIT (its `LICENSE` is in `core/vendor/gb-io`), PyO3 is
  MIT/Apache-2.0, 3Dmol.js is BSD-3-Clause. The PyInstaller build bundles Qt (LGPL), which asks that users can
  replace the Qt libraries: the installed folder keeps Qt as plain files for that (see "Qt (LGPL)" below). (5) If the work was done
  as a student or an employee, check the institution's intellectual-property policy. This is a pointer, not legal advice.
- **Public or private first?** Creating it private, checking the page, then flipping it public costs nothing.
- **The Sponsors startup page.** `app/sponsors.json` has three empty "Your ad here" slots shown on startup (labelled as sponsored,
  shown at the first start of each day, can be turned off). Decide that you want that in a public tool; if not, remove the page or leave the slots empty
  and the `contact` blank (as now).
- **Honest status.** The README says "early / personal project". [test-report/README.md](test-report/README.md) lists every test; what is known to be wrong is in
  [KNOWN_ISSUES.md](KNOWN_ISSUES.md). Keep both in step.

## 2. Already checked (as of this commit)

- History has one author with a GitHub no-reply address (`...@users.noreply.github.com`), so no private email is exposed.
- No keys or tokens were ever committed. Your own `.dna` / `.gb` / `.fa` files in the repo root are git-ignored
  (`.gitignore`), and so are `core/target/`, `build/`, `dist/`, `__pycache__/` and `tests/fake_reads/`.
- The data files shipped are `demo/demo_*.gb` (made-up random sequences from `scripts/make_demo_files.py`), the example
  sequences in `demo/` that you made yourself and own (`pBAD30_eilR.gb`, `pFG051_Chlor.gb`, `pFG051_Chlor.dna` and the
  27 MB `BW25113_2026.gb`, which the tests also use), and `app/parts.tsv` (the common-feature library), plus `app/vendor/3Dmol-min.js` (3Dmol.js 2.4.2, BSD-3-Clause, from cdnjs; its licence is
  `app/vendor/3Dmol.LICENSE`; it must be bundled in the exe and keep its licence file).

Re-run the history check before pushing, since it takes seconds:

```bash
git log --format='%an <%ae>' | sort -u
git log --all --name-only --format= | sort -u | grep -i -E '\.(dna|gb|gbk|bam|fa|fasta|ini|env)$|secret|token|key'
```

## 3. Prove a clean clone builds

The most common first issue on any repo is "it does not build". Do this once, in a throwaway folder:

```powershell
git clone . ..\openDNA-clean
cd ..\openDNA-clean
uv venv --python 3.12 .venv; .venv\Scripts\activate
uv pip install -r requirements-dev.txt
maturin develop -m core/Cargo.toml --release
python scripts/make_demo_files.py
pytest -q tests
```

`cargo test` needs to find Python: set `PYO3_PYTHON` to that environment's `python.exe` (see the README).
Tests that need your own files skip themselves when the files are absent, so a clean clone should be green.

## 4. Create the repository

With the [GitHub CLI](https://cli.github.com/) (`gh auth login` once):

```powershell
gh repo create openDNA --private --source . --remote origin --description "Local DNA editor with Rust core and PySide6 GUI"
git push -u origin main
```

Add topics (`bioinformatics`, `dna`, `plasmid`, `golden-gate`, `gibson-assembly`, `pyside6`, `rust`, `pyo3`) in the
repository settings, then make it public when ready: `gh repo edit --visibility public --accept-visibility-change-consequences`.

## 5. Make it look good

- The README is a showcase page: one section per feature, each with a picture. `python scripts/make_screenshots.py` makes
  every picture again from the real program and the sequences in `demo/`; redo them whenever the interface changes
  (`CONTRIBUTING.md` says how).
- The README download links point at the releases page (`releases/latest`), because the file names carry the version (`oDNA_beta_X_Y_Z-setup.exe`).
- Optional: a code of conduct. Skip it until someone asks (`CONTRIBUTING.md` and the issue templates are in place).

## 6. Continuous integration

`.github/workflows/ci.yml` runs on every push to `main` and on pull requests (Windows, about 5 minutes): it builds the Rust core, runs
`cargo test`, `pyflakes` and `pytest` (the tests that need your own `.dna` / `.fa` files skip themselves on the runner). The workflow uses `actions/checkout@v5` and `actions/setup-python@v6`; bump them when newer majors appear. The cargo cache keeps runs short.

## 7. Building the Windows installer

One command builds the program as a folder, tests it, wraps it in the installer `dist\openDNA-setup.exe` (Inno Setup, per user, no
administrator rights), and tests that (installs it into a temporary folder, starts it, uninstalls it). Install Inno Setup 6 first (`winget install JRSoftware.InnoSetup`); no Python is needed on the PC that runs the result:

```powershell
pip install -r requirements-dev.txt
maturin develop -m core/Cargo.toml --release
powershell -NoProfile -ExecutionPolicy Bypass -File scripts\build_exe.ps1
```

(`-ExecutionPolicy Bypass` is needed on a PC where PowerShell scripts are blocked, which is the Windows default for a normal user.
Options: `-NoSmoke` skips the tests, `-NoInstaller` skips the installer (and the need for Inno Setup), `-Python <path>` or `$env:PYTHON` picks the interpreter; it must be the one that has
`opendna_core` and PyInstaller.)

What is where:

| File | What it does |
|---|---|
| `packaging/openDNA.spec` | The PyInstaller recipe. Lists every data file the app reads (`parts.tsv`, `tools.toml`, `sponsors.json`, `vendor/3Dmol-min.js` and its licence, `LICENSE`, and `THIRD_PARTY_NOTICES.md` when present). Name, version, author and copyright of the Windows file properties come from `app/project.py`. No UPX packing (antivirus programs dislike it). |
| `packaging/openDNA.ico`, `packaging/make_icon.py` | The icon (an open plasmid ring with a promoter and a terminator, on a light tile, 16 to 256 px) and the script that draws it with Qt. The same script writes the window icons and wordmarks in `app/logo/` and the SVGs in `docs/images/`. Re-run `python packaging/make_icon.py` only to change the picture. |
| `scripts/build_exe.ps1` | Checks the setup, runs PyInstaller on the spec (folder layout), runs the smoke test, builds the installer, runs the installer smoke test, prints sizes and SHA-256. |
| `packaging/openDNA.iss` | The Inno Setup script: per-user install in `%LOCALAPPDATA%\Programs\openDNA`, Start Menu entry, file types registered at install (`openDNA.exe --register-file-types`) and removed at uninstall. The updater runs it silently with `/RELAUNCH=1`, which starts the new version at the end. `unins000.exe` beside the exe is how a running copy knows it was installed (`updater.installed_copy`). User data is in `%APPDATA%\openDNA` and is never touched. |
| `scripts/smoke_installer.ps1` | Installs `dist\openDNA-setup.exe` silently into a temporary folder, checks the uninstaller and the file-type registration, runs `smoke_exe.ps1` on the installed exe, uninstalls, checks that the exe and the registration are gone. |
| `scripts/smoke_exe.ps1` | Starts the exe on a PATH that has only Windows on it (no Python, conda or Rust), with `demo/demo_vector.gb` on the command line. Passes when the window shows the file's name, the process is still alive 5 seconds later, no error dialog is showing, every bundled data file is inside the unpacked bundle, `opendna_core` is in it, and no DLL was loaded from Python, conda, Rust or MSYS. |
| `tests/test_project.py` | Fails if the version differs between `app/project.py`, `core/Cargo.toml` and `core/pyproject.toml`, if a file the app reads with `Path(__file__).with_name(...)` is not in the spec, or if the icon is damaged. |

If you add a file the app reads at run time, add it to `DATA` in the spec (and to `$expected` in `smoke_exe.ps1`). Forgetting it is the
classic frozen-app bug: the exe starts fine and only the feature that needs the file fails. The test above catches the common case.

Known limits: the build is not bit-for-bit reproducible (PyInstaller stamps the build time), so two builds of the same tag have
different hashes; the release's `.sha256` is for the file you upload. The smoke test proves the exe starts and that its data is
there; it does not click through the features (the pytest suite does that, on the source).

## 8. Releases

A release is one tag. Pushing a tag `v*` starts `.github/workflows/release.yml`: it checks the tag against `VERSION`, then builds on
three machines. Windows runs the Rust and Python tests and the script above (the installer, smoke-tested); macOS
(Apple Silicon and Intel) and Linux run `scripts/build_unix.sh`, which smoke-tests what it builds. Each file gets a `.sha256`
while building; the release job joins them into one `SHA256SUMS` (`<hex>  <file name>` per line, the format `sha256sum -c` and the in-app updater read) and all go into a **draft** GitHub release:

| System | Files | What the program updates itself with |
|---|---|---|
| Windows | `oDNA_beta_X_Y_Z-setup.exe` | the setup file, run silently |
| macOS | `-macos-arm64.dmg` / `-macos-x86_64.dmg` (first install), `-macos-<arch>.zip` | the zip: the `.app` is replaced in place |
| Linux | `-linux-x86_64.AppImage` | the AppImage: the file named by `$APPIMAGE` is replaced |

`updater.release_asset` and `release.yml` must agree on these names. Copies run from source never
replace themselves: they say where to download the new version. The release never publishes by itself: you do that by hand.

1. **Bump the version in three places** to the same value (say 0.1.1): `VERSION` in `app/project.py`, `version` in
   `core/Cargo.toml`, `version` in `core/pyproject.toml`. Building changes `core/Cargo.lock` too: commit it.
   `pytest tests/test_project.py` fails if the three differ.
2. **Try it locally.** Run the build script (section 7) and install `dist\openDNA-setup.exe` yourself once.
3. **Commit, push to `main`, wait for CI to go green** (section 6).
4. **Tag and push the tag.** The tag must be `v` plus the version exactly, or the workflow stops at its first step:

   ```powershell
   git tag v0.1.1
   git push origin v0.1.1
   ```

5. **Watch the run** (Actions > release, about 10 to 15 minutes). When it is green, a draft appears under Releases.
6. **Review the draft.** Every file in the table above is attached, with one `SHA256SUMS`, the notes say what changed (edit them; the
   workflow only writes a default text). Put what changed BEFORE the **Downloads.** line: the program's update window shows only the text above it. Download the setup file from the draft, check `Get-FileHash` against its line of `SHA256SUMS`, and
   install it on a PC that has no Python. To test the update itself, install the previous release first, then publish this one.
7. **Publish** by pressing *Publish release* on the draft page, or `gh release edit v0.1.1 --draft=false --latest`. Only now does the
   release become public and the "latest" one that the in-app update check reads (`UPDATES_ENABLED` in
   `app/project.py`, on since 0.7.0 now that the repository is public).

To try the pipeline without a tag: Actions > release > *Run workflow* (from any branch). It builds the same files and creates a draft
named after the version in the code, with no tag until published. Delete the draft afterwards (`gh release delete v0.1.1`).
Something wrong after step 4: `gh release delete v0.1.1 --cleanup-tag` removes the draft and the tag on GitHub, then `git tag -d v0.1.1`
locally, fix, tag again. Never move a tag that was published.

The release workflow has run for v0.8.2 to v0.8.4.

### Code signing and SmartScreen

The exe is **not signed**, and nothing here pretends otherwise. What users will see:

- **Windows SmartScreen** ("Windows protected your PC", *More info* > *Run anyway*) on first run of a downloaded unsigned exe. It goes
  away only as the file builds reputation (many downloads over time), and a new build restarts that. A code-signing certificate
  helps, but is not an instant fix: since 2024 even EV certificates no longer give immediate reputation. The release notes say this.
- **A signature also protects the updater.** The `.sha256` file catches a corrupted download, but not a stolen GitHub account (an
  attacker would replace both files). Only a signature, checked by the updater, closes that gap (see the note at the top of
  `app/updater.py`).

Options when it matters (prices and eligibility change, check them): *SignPath Foundation* signs open-source projects for free, with
conditions (public repository, an OSI licence, a build made in CI, which this workflow is); Microsoft's *Trusted Signing* (now called
Artifact Signing) is a monthly fee, with identity checks and country limits; a classic OV certificate from a certificate authority
costs a few hundred dollars a year and the key must be on a hardware token. Signing would be one more step in `release.yml` after the
build and before the checksum. Do not put a certificate or its password in the repository.

### Windows Defender and other antivirus false positives

Unsigned PyInstaller one-file exes are a well-known source of wrong detections (the bootloader is shared by thousands of programs,
some of them malware). It can happen to any release, and to some users and not others. What lowers the risk (all already done or
cheap): no UPX packing (done), a complete version resource and icon (done), building in CI from a public tag so the build is
traceable, and, if it keeps happening, building PyInstaller's bootloader from source (`pip install pyinstaller --no-binary pyinstaller`,
needs a C compiler) so it differs from the common one.

When Defender flags a release:

1. Note the detection name (Windows Security > Protection history).
2. Submit the file to Microsoft: <https://www.microsoft.com/en-us/wdsi/filesubmission>. Sign in, choose *Software developer*, upload
   `openDNA.exe`, tick that it is a false positive (the detection name goes in the box), and say it is open source with the
   repository link and the tag. Answers usually come in a few days, and the fix reaches users through definition updates.
3. For other vendors use their false-positive pages (VirusTotal lists which ones flagged the file and links to them).
4. Put a line in the release notes while it lasts: how to verify the hash, and that the source and the build recipe are public.

Uploading to VirusTotal shares the file with the antivirus vendors: do it for a public release, not for a build you keep private.

### Qt (LGPL): the installed folder

openDNA's licence is AGPL-3.0, but the Qt libraries inside the exe (PySide6) are LGPL-3.0, and the LGPL asks that a user can
**replace the Qt libraries with their own build** and still run the program. This is a pointer to that obligation, not legal advice;
if the project grows, ask someone qualified. [THIRD_PARTY_NOTICES.md](../THIRD_PARTY_NOTICES.md) says how to get the Qt source
and how to replace Qt.

A PyInstaller one-file exe packs the Qt DLLs inside the executable, so a user cannot just drop in another DLL. The Windows releases
therefore no longer contain a one-file exe: the installer carries the one-folder layout (`packaging/openDNA.spec` with
`OPENDNA_BUILD=folder`, which `scripts/build_exe.ps1` sets), where `_internal\PySide6` and `_internal\shiboken6` hold plain files:

`-setup.exe` installs that folder per user, so Qt stays replaceable: the installed `_internal` folder (an update puts the release's Qt
back). No portable zip is released: the installer is small, and the LGPL is met by the plain Qt files it installs, by the AppImage on Linux
(also a folder build) and by the public source.

The installed folder also holds `README.txt` (how to swap Qt), `LICENSE` and `THIRD_PARTY_NOTICES.md` next to the exe. The smoke test
runs on both layouts and, for the folder, checks that `Qt6Core.dll` is a plain file.
Help > About names the Qt and PySide6 versions in use. Tested once by hand (version 0.1.0, built with Qt 6.12.0): ran the installed folder, followed the
`README.txt` steps with the PySide6-Essentials and shiboken6 6.10.1 wheels, and the program started and passed `smoke_exe.ps1` with Qt 6.10.1 loaded.
Redo that after a PyInstaller or PySide6 upgrade; the smoke test only proves that the Qt DLL is a plain file, not that a swap still works.
The updater only looks for the file its system installs (table in section 8) and `SHA256SUMS` in a release; the other files are ignored by it.
