# Publishing openDNA on GitHub

A checklist for taking the project further. **Done so far:** the code is in a private GitHub repository under the AGPL-3.0 ([LICENSE](../LICENSE)). The steps below are what is left before opening it up, building releases or adding CI.

## 1. Decide first

- **Licence.** Chosen: AGPL-3.0 ([LICENSE](../LICENSE)), kept as the licence of the project. It keeps the source open and
  obliges every fork or hosted copy to stay open under the same licence, so changes can be pulled back; it does not prevent
  forks, and is not meant to. Notes: (1) on a **public** GitHub repository anyone can fork it whatever the licence says.
  (2) To keep one project, the other levers: the **name policy** ([TRADEMARKS.md](../TRADEMARKS.md): forks are welcome under a
  name and an icon of their own; the AGPL's section 7(e) lets the author decline trademark rights), welcoming pull requests,
  and releases that work. A registered trademark on "openDNA" / "oDNA" would back the policy up; none is registered yet.
  (3) Before accepting pull requests from other people, decide whether to ask for a contributor agreement (it keeps the
  option of changing the licence later); without one, each contributor keeps the copyright of their part, under the AGPL.
  (4) The dependencies are compatible: PySide6 is LGPL, `gb-io` is MIT (its `LICENSE` is in `core/vendor/gb-io`), PyO3 is
  MIT/Apache-2.0, 3Dmol.js is BSD-3-Clause. The PyInstaller one-file `.exe` bundles Qt (LGPL), which asks that users can
  replace the Qt libraries: releases also carry a one-folder zip for that (see "Qt (LGPL)" below). (5) If the work was done
  as a student or an employee, check the institution's intellectual-property policy. This is a pointer, not legal advice.
- **Public or private first?** Creating it private, checking the page, then flipping it public costs nothing.
- **The Sponsors startup page.** `app/sponsors.json` has three empty "Your ad here" slots shown on startup (labelled as sponsored,
  shown at the first start of each day, can be turned off). Decide that you want that in a public tool; if not, remove the page or leave the slots empty
  and the `contact` blank (as now).
- **Honest status.** The README says "early / personal project". [TESTING.md](TESTING.md) lists what is still
  untested (two items: a comparison with a .dna file's own primer list on a real file that has primers, and the location forms listed under Known problems in TODO.md). Keep both in step.

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
conda create -n ge-clean python=3.13 -y; conda activate ge-clean
pip install -r requirements-dev.txt
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

- Add 2-3 screenshots to `docs/images/` and reference them from the README. `GE_SHOTS=<folder>` makes the GUI tests
  save screenshots (see the README); the Assembly panel (the junction cards) is the most distinctive.
- Pin the README's *Quick start* to a tested path (section 3 does exactly that).
- Optional: `CONTRIBUTING.md`, issue templates, a code of conduct. Skip them until someone asks.

## 6. Continuous integration

`.github/workflows/ci.yml` runs on every push to `main` and on pull requests (Windows, about 5 minutes): it builds the Rust core, runs
`cargo test`, `pyflakes` and `pytest` (the tests that need your own `.dna` / `.fa` files skip themselves on the runner). The first run
passed. GitHub warns that `actions/checkout@v4` and `actions/setup-python@v5` use Node 20: bump them when newer majors appear.
Private-repo minutes are limited (Windows counts double); the cargo cache keeps later runs short.

## 7. Building the exe

One command builds `dist\openDNA.exe` (one file, about 39 MB, no Python needed on the PC that runs it) and then tests it:

```powershell
pip install -r requirements-dev.txt
maturin develop -m core/Cargo.toml --release
powershell -NoProfile -ExecutionPolicy Bypass -File scripts\build_exe.ps1
```

(`-ExecutionPolicy Bypass` is needed on a PC where PowerShell scripts are blocked, which is the Windows default for a normal user.
Options: `-NoSmoke` skips the test, `-Python <path>` or `$env:PYTHON` picks the interpreter; it must be the one that has
`opendna_core` and PyInstaller.)

What is where:

| File | What it does |
|---|---|
| `packaging/openDNA.spec` | The PyInstaller recipe. Lists every data file the app reads (`parts.tsv`, `tools.toml`, `sponsors.json`, `vendor/3Dmol-min.js` and its licence, `LICENSE`, and `THIRD_PARTY_NOTICES.md` when present). Name, version, author and copyright of the Windows file properties come from `app/project.py`. No UPX packing (antivirus programs dislike it). |
| `packaging/openDNA.ico`, `packaging/make_icon.py` | The icon (an open plasmid ring with a promoter and a terminator, on a light tile, 16 to 256 px) and the script that draws it with Qt. The same script writes the window icons and wordmarks in `app/logo/` and the SVGs in `docs/images/`. Re-run `python packaging/make_icon.py` only to change the picture. |
| `scripts/build_exe.ps1` | Checks the setup, runs PyInstaller on the spec, prints the size, SHA-256 and file properties, then runs the smoke test. |
| `scripts/smoke_exe.ps1` | Starts the exe on a PATH that has only Windows on it (no Python, conda or Rust), with `demo/demo_vector.gb` on the command line. Passes when the window shows the file's name, the process is still alive 5 seconds later, no error dialog is showing, every bundled data file is inside the unpacked bundle, `opendna_core` is in it, and no DLL was loaded from Python, conda, Rust or MSYS. |
| `tests/test_project.py` | Fails if the version differs between `app/project.py`, `core/Cargo.toml` and `core/pyproject.toml`, if a file the app reads with `Path(__file__).with_name(...)` is not in the spec, or if the icon is damaged. |

If you add a file the app reads at run time, add it to `DATA` in the spec (and to `$expected` in `smoke_exe.ps1`). Forgetting it is the
classic frozen-app bug: the exe starts fine and only the feature that needs the file fails. The test above catches the common case.

Known limits: the build is not bit-for-bit reproducible (PyInstaller stamps the build time), so two builds of the same tag have
different hashes; the release's `.sha256` is for the file you upload. The smoke test proves the exe starts and that its data is
there; it does not click through the features (the pytest suite does that, on the source).

## 8. Releases

A release is one tag. Pushing a tag `v*` starts `.github/workflows/release.yml` (Windows, same setup as CI): it checks the tag
against `VERSION`, runs the Rust and Python tests, builds the exe with the script above (which smoke-tests it), writes
`openDNA.exe.sha256` (`<hex>  openDNA.exe`, the format `sha256sum -c` and the in-app updater read), and creates a **draft**
GitHub release holding both files. It never publishes: you do that by hand.

1. **Bump the version in three places** to the same value (say 0.1.1): `VERSION` in `app/project.py`, `version` in
   `core/Cargo.toml`, `version` in `core/pyproject.toml`. Building changes `core/Cargo.lock` too: commit it.
   `pytest tests/test_project.py` fails if the three differ.
2. **Try it locally.** Run the build script (section 7) and start `dist\openDNA.exe` yourself once.
3. **Commit, push to `main`, wait for CI to go green** (section 6).
4. **Tag and push the tag.** The tag must be `v` plus the version exactly, or the workflow stops at its first step:

   ```powershell
   git tag v0.1.1
   git push origin v0.1.1
   ```

5. **Watch the run** (Actions > release, about 10 to 15 minutes). When it is green, a draft appears under Releases.
6. **Review the draft.** Both files are attached, the size is about 39 MB, the notes say what changed (edit them; the workflow only
   writes a default text). Download the exe from the draft, check `Get-FileHash openDNA.exe` against the `.sha256` file, and run it on
   a PC that has no Python.
7. **Publish** by pressing *Publish release* on the draft page, or `gh release edit v0.1.1 --draft=false --latest`. Only now does the
   release become public and the "latest" one that the in-app update check reads (once `UPDATES_ENABLED` is switched on in
   `app/project.py`, which must wait for the repository to be public).

To try the pipeline without a tag: Actions > release > *Run workflow* (from any branch). It builds the same files and creates a draft
named after the version in the code, with no tag until published. Delete the draft afterwards (`gh release delete v0.1.1`).
Something wrong after step 4: `gh release delete v0.1.1 --cleanup-tag` removes the draft and the tag on GitHub, then `git tag -d v0.1.1`
locally, fix, tag again. Never move a tag that was published.

The workflow file has only been checked by reading and by running each of its PowerShell steps locally; it has not run on GitHub yet.
The first real run may need a small fix (most likely candidates: the smoke test on the runner's desktop, which has no GPU, and
`gh release create` permissions in a private repository).

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

### Qt (LGPL): the one-file exe and the one-folder zip

openDNA's licence is AGPL-3.0, but the Qt libraries inside the exe (PySide6) are LGPL-3.0, and the LGPL asks that a user can
**replace the Qt libraries with their own build** and still run the program. This is a pointer to that obligation, not legal advice;
if the project grows, ask someone qualified. [THIRD_PARTY_NOTICES.md](../THIRD_PARTY_NOTICES.md) says how to get the Qt source
and how to replace Qt.

In a PyInstaller one-file exe the Qt DLLs are packed inside the executable and unpacked to a temporary folder at each start, so a user
cannot just drop in another DLL. So every release carries **two downloads built from the same recipe** (`packaging/openDNA.spec`;
the environment variable `OPENDNA_BUILD=folder` selects the folder layout, and `scripts/build_exe.ps1` builds both, `-NoFolder` skips the folder):

| | What it is | Qt replaceable | Updates itself |
|---|---|---|---|
| `openDNA.exe` | one file, about 40 MB | only by rebuilding from the public source | yes, once `UPDATES_ENABLED` is on |
| `openDNA-folder.zip` | the exe plus an `_internal` folder, about 40 MB zipped | yes: `_internal\PySide6` and `_internal\shiboken6` hold plain files that can be swapped (steps in the `README.txt` inside the zip) | no: `updater.is_folder_build()` makes it say "download the new zip", because swapping one exe would leave the user's Qt files out of step with it |

The zip also holds `README.txt` (how to swap Qt), `LICENSE` and `THIRD_PARTY_NOTICES.md` next to the exe. `release.yml` attaches the zip
and its `.sha256` to the draft release; the smoke test runs on both layouts and, for the folder, checks that `Qt6Core.dll` is a plain file.
Help > About names the Qt and PySide6 versions in use. Tested once by hand (version 0.1.0, built with Qt 6.12.0): unzipped the zip, followed the
`README.txt` steps with the PySide6-Essentials and shiboken6 6.10.1 wheels, and the program started and passed `smoke_exe.ps1` with Qt 6.10.1 loaded.
Redo that after a PyInstaller or PySide6 upgrade; the smoke test only proves that the Qt DLL is a plain file, not that a swap still works.
The updater only looks for `openDNA.exe` and `openDNA.exe.sha256` in a release; the zip is ignored by it.
