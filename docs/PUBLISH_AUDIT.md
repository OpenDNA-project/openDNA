# Pre-publication audit

Checklist for making the repository public for testers. Kept up to date with the tree, so that what it says can be
trusted at a glance; the ranked items below are the work that is still open.

This is an engineering checklist, **not legal advice**. Where it says a licence is compatible, read "as far as we can tell".

Ranking: **BLOCKER** = fix before flipping the repository public. **SHOULD-FIX** = fix before or just after the first
release. **NICE** = polish.

## Summary

- The git history is clean: one author (a GitHub no-reply address), no keys or tokens, no real sequence files, no personal
  paths, nothing deleted that still sits in history. No blocker comes from the history.
- The blockers are about **what a release contains** (licence texts, the 3D viewer file), the **Qt/LGPL decision for the
  one-file exe**, and **a working private channel for security reports**.
- The program contacts the internet only when you use a network feature (list in section 5). There is no telemetry.

---

## BLOCKERS

### B1. A release must carry the licence texts and the 3D viewer, and the old build command does not

- The README's old PyInstaller command (`README.md`, the block under "Build from source") bundles only `parts.tsv`,
  `tools.toml`, `sponsors.json`. It does **not** bundle `app/vendor/3Dmol-min.js`, which `app/structure.py:39` (`VENDOR`) reads
  when the 3D page is made (`structure.py:251`). An exe built with it would fail on *Open 3D view* (only that feature, not at
  start-up), and would not carry `app/vendor/3Dmol.LICENSE` (BSD-3-Clause requires the notice with binary redistributions).
- Evidence: the `dist/openDNA.exe` that was in the working tree on 2026-11-01 (dated 2026-10-29) contained no `3Dmol`
  string (it has since been removed from `dist/`).
- The new, untracked `packaging/openDNA.spec` lists the vendor files, `LICENSE` and `THIRD_PARTY_NOTICES.md`. **To do:**
  commit it, build the release only with it, then check the exe on a PC without Python, including *Open 3D view in the
  browser* (the smoke test mentioned in the spec should assert this).
- The zip that users download should hold: the exe, `LICENSE`, `THIRD_PARTY_NOTICES.md`, `README.md`, and a link to the
  matching source tag (AGPL: the source must be offered for what you convey; a link to the tag on the same site is the usual way).

### B2. Qt (LGPL) inside a PyInstaller one-file exe: decide, and say what you decided

PySide6 6.12.0 is licensed `LGPL-3.0-only OR GPL-2.0-only OR GPL-3.0-only`; openDNA uses it under the LGPL-3.0. In plain terms
the LGPL asks, for an app that bundles Qt: (1) say that Qt is used and under which licence, and ship the licence text;
(2) give the corresponding Qt source (or a written offer / link: Qt's source is public at download.qt.io; PySide6's at
code.qt.io); (3) let users replace the Qt libraries with their own build and run your program with them (and allow
reverse engineering for debugging such a change).

A one-file exe packs the Qt DLLs inside itself and unpacks them to a temporary folder at each start, so a user cannot simply swap
a DLL. Recommendation (not legal advice), cheapest first: (a) also publish a **one-folder** build (zip of `dist\openDNA\`),
where the Qt DLLs are ordinary files that can be replaced; or (b) keep the one-file exe and state, in the README and
`THIRD_PARTY_NOTICES.md`, that anyone can rebuild it from source with another PySide6 version (the build is scripted and the source
is AGPL). `THIRD_PARTY_NOTICES.md` already says (b) and gives the Qt source locations; if you choose (a), add one line.
**Resolved: both.** Each release now also carries `openDNA-folder.zip` (a one-folder build with a `README.txt` on swapping Qt), built by
`scripts/build_exe.ps1` and `release.yml`; see `docs/PUBLISHING.md`, "Qt (LGPL)".
Also ship the LGPL-3.0 and GPL-3.0 texts or link them (the notices file links them; Qt's own wheel metadata carries them).
GPL-3.0 and AGPL-3.0 are expressly allowed to be combined (section 13 of each), so the AGPL project and LGPL Qt do not conflict.

### B3. Security reports need a channel that works before the repository is public

`SECURITY.md` (new) tells people to use GitHub's private vulnerability reporting. That only works once it is switched on:
repository **Settings > Code security > Private vulnerability reporting > Enable** (it is off by default). Also fill in
`CONTACT_EMAIL` in `app/project.py:9` and the `<project email>` placeholders in `SECURITY.md` and `PRIVACY.md` once a
project-only address exists, or delete those lines. Do not use a personal address.

### B4. Commit `app/project.py` (and what imports it) together — done

`app/project.py` is untracked, but `app/sponsors.py` (modified) and `app/main.py` (modified) import it, so a clean clone of the
last commit would not match what you tested, and a commit that forgets the file would not start. Its docstring says
`tests/test_project.py` checks that `VERSION` equals the version in `core/Cargo.toml` and `core/pyproject.toml`; **that test
file does not exist yet** (see S8). Run the clean-clone check of `docs/PUBLISHING.md` section 3 on the final commit.

---

## SHOULD-FIX

### S1. 3D page: text from other servers is put inside a `<script>` without escaping `</script>` — done

`app/structure.py` (`viewer_html`, the line `const pdb = {json.dumps(structure.pdb)}`) writes the PDB text, which comes from RCSB,
the ESM Atlas, AlphaFold DB or a local command, straight into an inline script. `json.dumps` does not escape `</script>` or
`<!--`, so a PDB `REMARK` line containing `</script><script>...` would end the script and run code in the local page (opened
from `file://` in the browser). The 3Dmol library text next to it *is* escaped (`.replace("</script", ...)`). The risk is low
(HTTPS to well-known services) but it is a two-line fix: apply the same `replace("</", "<\\/")` to the JSON. Add a test.

### S1b. Self-updater of an unsigned exe (new code): make the trust model explicit and keep it off until you mean it

`app/updater.py` replaces the running exe with a download from GitHub. Good: HTTPS only, GitHub hosts only (each redirect hop checked), SHA-256 compared, user must click.
Limits it states itself: the checksum file lives in the same release as the exe, so a hijacked release (stolen GitHub account or token) would be accepted; the exe is unsigned. For a
test program with a handful of testers this is a normal trade-off, but (1) say it in the README and `SECURITY.md` (done), (2) protect the GitHub account with 2-factor
authentication and a hardware or passkey login, and enable "require signed commits / immutable releases" if offered, (3) leave `UPDATES_ENABLED = False` until the repository is public and
the first release exists (flipping it earlier shows a message and is harmless), (4) make the SHA-256 of the exe visible in the release notes too, so a human can compare it
without the program, (5) consider a signature (minisign or Sigstore) of `openDNA.exe.sha256` with a public key stored in the repository, which would defeat a replaced pair of files.

### S2. AlphaFold DB: the model URL taken from the server's answer is not checked — done

`structure.py:125`: `fetch(entry["pdbUrl"])`. Python's `urlopen` also opens `file://` and `ftp://` URLs, so a wrong or hostile
answer could make the program read a local file and place its text in the page (`from_alphafold` does not call
`looks_like_pdb`). Require `https://` and a host ending `ebi.ac.uk` (or the AlphaFold download host), and call
`looks_like_pdb` on the result like the other sources do. Same idea for `fetch` in general: refuse non-`https` schemes.

### S3. The data folder is not what the docs say, and settings are split in three places

`app/main.py` never calls `setOrganizationName` / `setApplicationName`, so Qt names the data folder after the executable:
`%APPDATA%\openDNA\` for `openDNA.exe`, but `%APPDATA%\python\` when run from source with `python app/main.py`
(both folders exist on the author's PC). `app/primer_db.py:3` and `docs/USER_GUIDE.md:17` say
`%APPDATA%\openDNA\openDNA\`, which does not exist. Consequences: documentation is wrong, and a user who tries the source
version and the exe sees two different primer databases. Fix: set organization/application name once at start-up (then
decide whether to migrate the old folder), and correct the two docs. `PRIVACY.md` states today's behaviour.

Also: `QSettings("openDNA", "openDNA")` (`app/primer_db.py:28`, `app/primer_panel.py:34`) writes to the **Windows registry**
(`HKCU\Software\openDNA\openDNA`: PCR conditions, primer options including the path of the chosen reference genome), while
every other setting goes to `.ini` files in the data folder. Moving these two to the same `.ini` style would keep "everything is in
one folder" true and make uninstalling clean. **Resolved:** both now use `primers.ini` in the data folder (`primer_db.settings()`); the
registry holds only the optional file-association entries. Older values under the registry key are not copied over.

### S4. Public wording that overclaims (README fixed; others remain)

- `README.md` (old): "82 Python tests, 64 Rust tests" was stale (430 Python tests are collected; about 80 `#[test]` in
  `core/src`). The rewritten README no longer quotes counts that go out of date.
- `README.md` (old): "Edits that cannot silently break annotations" and "provably unique in the product": softened in the rewrite
  ("designed so that", "checked").
- `PLAN.md` used to describe the program as a replacement for another commercial plasmid editor and as "-like" it. The README links
  PLAN.md, and a public "replacement for / compatible with <product>" statement invites a trademark or endorsement complaint and
  testers will hold you to it. Reworded: the product name now appears nowhere in the repository except the trademark line in
  `THIRD_PARTY_NOTICES.md` and the file signature the `.dna` reader checks.
- `docs/TODO.md:75` says "about 50" items remain in `TESTING.md`; `docs/PUBLISHING.md:19` says "about 60". Actually 2 are unticked
  (lines 87 and 138, one of them a known bug marked xfail) and 149 are ticked. Update both so the public docs agree with TESTING.md.
- `app/sponsors.py` text and `About` (`app/main.py:431-441`) are fine ("NO WARRANTY", "keep a copy of your original files"). Keep
  "early test version" in both until the open items are done.

### S5. The bundled part library (`app/parts.tsv`): provenance and terms

94 sequences, each with an accession in its `source` column, built by `scripts/build_parts.py` from 22 public GenBank/INSDC records
of well-known vectors (pBR322 J01749, pUC19 L09137, pEGFP-N1 U55762, pKD46 AY048746, pCas9 KM099231, pCAMBIA-1300 AF234296 ...).
Problems: (1) the script's docstring says "NCBI data is public domain". NCBI's own policy is more careful: NCBI places no
restrictions on the *database*, but submitters may claim copyright or patent rights in a record, and NCBI cannot say they do not.
Replace the sentence in `scripts/build_parts.py:2` with that wording. (2) Nothing in `README.md` or Help > Citations lists the
parts library or its source records. Add a Citations entry (the project rule: new panels built on datasets go in `app/citations.py`) and
keep the accession column visible. (3) Several of these are commercial vectors (pEGFP-N1, pCAMBIA, pGEM-luc) whose sequence is
published; as far as we can tell short functional elements (promoters, origins, resistance genes) are facts, but the list should be
reviewed once by the author, record by record, before it ships in a public binary. `THIRD_PARTY_NOTICES.md` states the current
position honestly.

### S6. Name of the Sidewinder reference is unverified

`app/citations.py:28`: "Robinson et al. Sidewinder DNA assembly. Nature (2026)", no DOI, and the entry itself says the full reference
was not recorded. The README and `docs/USER_GUIDE.md` cite it. Before the public release, look the paper up, put the real citation
(authors, title, volume, DOI) in `citations.py`, and re-check the strand layout (the code and the docs both say it is simplified).
The README rewrite keeps the existing hedge. This is the one citation that could be wrong.

### S7. Python version

`app/runner.py` imports `tomllib` (Python 3.11+), while `core/pyproject.toml` says `requires-python = ">=3.10"` and the PyO3 build
targets `abi3-py310`. The README says Python 3.13 (what CI and the author use). State "3.11 or newer" or pin 3.13 in
`core/pyproject.toml`'s docs. Minor, but a tester on 3.10 would get a traceback.

### S8. Tests and docs referenced but missing

- ~~`tests/test_project.py` (promised by the docstring of `app/project.py`) does not exist.~~ Written: 9 tests.
- The `README` tests line `pytest -s tests` and `docs/PUBLISHING.md` `pytest -q tests` differ; both work. `CONTRIBUTING.md` uses `-q`.
- `docs/PUBLISHING.md:24` says the history check was done "as of this commit"; re-run it at the final commit (commands in section 1 below).

### S9. A tester's own files reach git by accident — done

`.gitignore` ignores `/*.dna /*.gb /*.gbk /*.fa /*.fasta` in the repository root, and `*.bam`, `*.bai`, `*.fastq`
anywhere in the tree. Sequence files are tracked only where they are meant to be, in `demo/`. Simulated reads live in
`tests/fake_reads/`, which is ignored by path. Still worth a glance at `git status` before committing, since a real
sequence dropped into `demo/` would be picked up by `git add .`.

---

## NICE

- N1. `LICENSE` is the unmodified AGPL-3.0 text (good). The source files have no per-file licence/SPDX header (only `app/main.py` mentions the
  licence in the About text). The AGPL's "how to apply" suggests a short header; `# SPDX-License-Identifier: AGPL-3.0-or-later`
  decision is yours (is it "only" or "or later"? README says just "AGPL-3.0": say which). The README rewrite says "AGPL-3.0".
- N2. Copyright line: the old README said "Copyright (C) 2026 Ant-Champ" (a GitHub login), `About` says "Antoine Champie". The rewrite uses the
  author's name; make the About, the exe's `LegalCopyright` (`packaging/openDNA.spec`) and the README agree.
- N3. The recent commit messages carry a `Co-Authored-By: Claude ...` trailer. Nothing wrong; if you want a plain
  statement of AI assistance, put one line in the README.
- N4. `PLAN.md:219` says "Ubuntu is already installed here" (the author's PC). Harmless.
- N5. `.github/workflows/ci.yml` uses `actions/checkout@v4`, `setup-python@v5`, `rust-cache@v2`, `dtolnay/rust-toolchain@stable` by tag. For a public repo
  consider pinning actions to commit hashes (supply-chain hygiene) and giving the workflow `permissions: contents: read`.
- N6. Name collision check for "openDNA" (other software, GitHub projects, trademark databases) was not done. A name policy for forks
  ("do not call a fork openDNA") is mentioned in `docs/PUBLISHING.md`; no trademark is registered.
- N7. `app/sponsors.json` is empty (`contact` blank, three empty slots) and the sponsors page now shows only the welcome text when no slot is
  filled. Fine. If you ever add a sponsor: it must stay labelled, and `PRIVACY.md` already says the program loads nothing for it.
- N8. The 3D page written to the data folder embeds 3Dmol.js (BSD-3-Clause) without the licence text next to it (the library's own header points
  to a `3Dmol-min.js.LICENSE.txt` that is not shipped). The notice is in `app/vendor/3Dmol.LICENSE` and `THIRD_PARTY_NOTICES.md`; putting the
  licence as an HTML comment in the generated page would be tidy.

---

## 1. Tree and history

Commands used (re-run before flipping public):

```bash
git log --all --format='%an <%ae>' | sort -u          # Ant-Champ <35179122+Ant-Champ@users.noreply.github.com>
git log --all --format='%cn <%ce>' | sort -u          # same
git log --all --name-only --format= | sort -u         # every path ever committed
git log --all --diff-filter=D --name-only --format=   # nothing deleted
git ls-files                                          # what is in the tree now
```

Findings:

- 64 commits, one branch (`main`), no tags, one remote (the private repository). Author and committer: the GitHub no-reply address only.
- No email address anywhere in the tree except `me@example.org` in `tests/test_ncbi.py:67,70` (a fake). No `C:\Users\...`, `V:\HM_Softwares` or similar absolute
  path in code, docs or tests (searched case-insensitively for `antoi`, `HM_Softwares`, `C:\`, `/Users/`, `champie`, `usherbrooke`): the only hits are
  generic (`V:\data\x.fq` in a docstring and `PLAN.md`). Note that `app/project.py` (untracked) names the author and a blank email field.
- No keys, tokens or passwords committed (searched the tree and every revision for `ghp_`, `AKIA`, private-key headers, `api_key`, `secret`, `password`).
  The words `token`/`api_key` in `app/` are an internal counter and the NCBI API-key field the *user* types (stored in the user's own `ncbi.ini`, see `PRIVACY.md`).
- Sequence data: the files shipped are `demo/demo_*.gb` (random, seeded sequences written by `scripts/make_demo_files.py`, "made-up sequences, not real
  parts"), the author's own example sequences in `demo/` (`pBAD30_eilR.gb`, `pFG051_Chlor.gb`, `pFG051_Chlor.dna` and `BW25113_2026.gb`, which the tests
  also read, through `tests/sample_files.py`), `demo/vibrio_cholerae_N16961.gb` (NCBI RefSeq NC_002505 / NC_002506, public domain US government data,
  read by `tests/test_multi_genbank.py`), and `app/parts.tsv` (see S5). The author made the four examples and owns them; they go out under the
  project's own licence (`THIRD_PARTY_NOTICES.md`). No `.bam`, `.bai` or `.fastq` is tracked anywhere.
- `docs/images/*.png` (7 screenshots, made by `scripts/make_screenshots.py`) show only sequences shipped in `demo/` (the author's own
  plasmids, the BW25113 genome, the RefSeq V. cholerae genome) and, in one picture, simulated reads made by `tests/fake_reads.py`.
  No personal path, user name or outside data is visible; the script writes its settings to a temporary folder, not the author's.
- `app/sponsors.json`: `{"contact": "", "slots": [{}, {}, {}]}`; nothing in it.
- Repository size: about 21 MB, most of it `demo/BW25113_2026.gb` (27 MB of text) and `demo/vibrio_cholerae_N16961.gb` (10 MB of text),
  both of which pack down well. No other large binaries in history.

## 2. Licensing

Everything below is "as far as we can tell"; see `THIRD_PARTY_NOTICES.md` for the per-component table (all 71 third-party Rust crates are listed there with
their licences, read from `cargo metadata --locked`).

| Piece | Licence | Compatible with AGPL-3.0? | Notice duty |
|---|---|---|---|
| PySide6-Essentials / shiboken6 6.12.0 (Qt 6) | LGPL-3.0 (or GPL-2/3) | Yes (dynamic use under the LGPL; GPL-3 option also combines with AGPL-3) | licence text, Qt source link, replaceability: see B2 |
| PyO3 0.23.5 | MIT or Apache-2.0 | Yes (Apache-2.0 is compatible with GPL-3/AGPL-3, not GPL-2) | notice in notices file |
| maturin 1.15.0 | MIT or Apache-2.0 | build tool, not shipped | none |
| `gb-io` 0.9.0 patched, `core/vendor/gb-io` | MIT | Yes | keep `core/vendor/gb-io/LICENSE`; the patches are marked "openDNA patch" in the source |
| `noodles` 0.118 and its 7 sub-crates, `quick-xml` 0.42 | MIT | Yes | notice |
| other Rust crates (71 third-party crates in total) | MIT, Apache-2.0, "MIT OR Apache-2.0", Unlicense OR MIT, Zlib (`zlib-rs`), Unicode-3.0 (`unicode-ident`), Apache-2.0 WITH LLVM-exception (`target-lexicon`, build time) | Yes | notice (MIT and Apache ask that the copyright and licence text travel with the binary: the notices file points to each crate's repository) |
| 3Dmol.js 2.4.2 (version as stated in `docs/PUBLISHING.md`), `app/vendor/3Dmol-min.js` | BSD-3-Clause (includes parts of GLmol, Three.js, jQuery) | Yes | the full notice in `app/vendor/3Dmol.LICENSE` must travel with it (B1) |
| PyInstaller 6.22.3 (build tool, bootloader ends up in the exe) | GPL-2.0-or-later with a special exception that allows non-free and other programs | Yes (the exception covers the bootloader) | none beyond mention |
| Python interpreter and its bundled libraries inside the exe | PSF-2.0 and the licences of OpenSSL, zlib, libffi, expat, etc. | Yes | mention; exact list depends on the build: check the built exe |
| `app/parts.tsv` | records of the public sequence databases | see S5 | attribution by accession (already in the file) |
| Enzyme list `app/enzymes.py` (47 enzymes: recognition site and cut offsets) | facts about enzymes | as far as we can tell | none; it was written by hand from common knowledge. Not taken from REBASE or a vendor catalogue (see `docs/TODO.md`: the NEB/REBASE catalogue is deliberately not bundled until permission) |
| Ladders `app/tools.py:280-288` | "typical band sets", explicitly not one vendor's list | n/a | the names ("1 kb", "1 kb Plus", "100 bp Plus") resemble product names; the docs say "typical sets" |
| Overhang data `app/overhangs.py` | 25 overhangs and junction sets from Potapov et al. 2018 (published supplementary tables) | facts from a paper | cited in `citations.py` |
| Test-only: pytest, Biopython, pyflakes | MIT, BSD-3-Clause / Biopython licence, MIT | not shipped | none |

Not bundled but contacted (their data terms are their own): NCBI records (public, submitters' rights may apply), AlphaFold DB models (CC BY 4.0),
RCSB PDB files (CC0), ESM Atlas predictions (Meta's terms). The structure window already says where a sequence goes; Help > Citations cites each.

## 3. Names and trademarks

Used nominatively (to say what it reads, what a method is, or what a service is). None is used in the program's own name, logo or colours.

| Name | Where | Comment |
|---|---|---|
| Name of the `.dna` format's owner | removed from the README, UI strings, comments and `PLAN.md`; files are called ".dna files" | Avoid "replacement for", "compatible with" and "-like" wording naming a product in public text (S4). No logo, colour scheme or product file is shipped. The colour qualifiers `ApEinfo_fwdcolor` are a de-facto format shared by several tools. |
| BLAST | none (the tool was renamed Local Align, `docs/TODO.md`) | good: the README says Local Align is not BLAST |
| Addgene, NEB, REBASE | only `docs/TODO.md` (planned, blocked for permission) and none in the program | fine |
| AlphaFold, ESMFold, PDB | `app/structure.py`, docs, citations | service names, cited; the window says what is sent to whom. No endorsement is implied by the README wording. |
| Biopython, Qt, PySide6, Python, GitHub, Windows, WSL, SmartScreen | docs | nominative |
| Ladder names ("1 kb Plus" ...) | `app/tools.py` | generic descriptions, disclaimed in the docs |

Claims about endorsement or compatibility: none found in the program's UI. `PLAN.md` is the only place (S4).

## 4. Claims in the docs (honesty check)

- README (rewritten): "early test version", the two unticked TESTING.md items named, no test counts, no "provably", no "cannot".
- `docs/TESTING.md`: honest and detailed; two items open (one is a known bug, xfail).
- `docs/USER_GUIDE.md`: states limits (Cas9 score is a checklist, ladder sets are typical, Sidewinder simplified, 3D model is a model).
  It has the wrong data-folder path (S3).
- `app/citations.py`: the Sidewinder entry is explicitly incomplete (S6).
- Test coverage claims in `docs/TESTING.md` such as "Agreement with the primer_bind stored in the real .dna file" rely on **files that are not in
  the repository** (they skip on a clean clone and in CI). Fine, but say so (the README does: "some tests need your own files and skip themselves").

## 5. Network behaviour (everything that contacts the internet)

Searched `app/` for `urllib`, `http`, `socket`, `QDesktopServices`, `QNetwork`, `webbrowser`, `openUrl`, `subprocess`, `QProcess`, `openExternalLinks`.
In the committed code only one function makes HTTP requests: `fetch` in `app/structure.py:60` (also used by `app/ncbi.py`). The working tree adds
`app/updater.py` (its own `fetch`, HTTPS and GitHub hosts only, redirect hops checked) and `app/feedback.py` (opens the browser only). There is no other network code.
These three files were not committed and not yet wired into `app/main.py` when this was written (`grep updater app/main.py` is empty): re-run the search
above on the final commit, and make sure `PRIVACY.md` still matches.

| Feature | Host | What is sent | When |
|---|---|---|---|
| File > Import from NCBI | `eutils.ncbi.nlm.nih.gov` (HTTPS) | your search words or accession; the tool name `openDNA` (also as User-Agent); your optional API key if you typed it (as URL parameters) | when you press Search / Open |
| 3D structure: AlphaFold DB | `alphafold.ebi.ac.uk` (then the model URL it returns) | the UniProt accession only | Get structure |
| 3D structure: PDB | `search.rcsb.org`, `files.rcsb.org` | the **protein sequence** (search), then a PDB code | Get structure (also in Auto) |
| 3D structure: ESMFold | `api.esmatlas.com` (Meta) | the **protein sequence** (up to 400 aa) | Get structure (also in Auto) |
| 3D structure: Local command | none by openDNA | the command *you* typed runs on your PC (`subprocess.run`, no shell); what it does is up to it | Get structure |
| Open 3D view in the browser | none: it opens a local `.html` file (`QDesktopServices.openUrl(file)`) in your default browser; the page has 3Dmol.js and the data inline and loads nothing | nothing | button |
| Tools > Install Mapping Tools | conda inside WSL downloads bwa, samtools, minimap2 from conda-forge/bioconda | what conda sends (openDNA sends nothing itself) | you run it |
| Links in About, Citations, Welcome (doi.org, GitHub) | your browser | nothing from openDNA | when you click |
| Update check (`app/updater.py`, new and untracked at the time of writing) | `api.github.com`, then GitHub release hosts for the download | a GET for "latest release" with User-Agent `openDNA/<version>`; the exe is downloaded only after the user presses the button | **off while `project.UPDATES_ENABLED` is False**; once true: once a day at start-up (setting *Tools > Settings > Updates*) and on *Check for Updates* |
| Help > Report a Problem / Request a Feature (`app/feedback.py`, new) | opens `github.com/OpenDNA-project/openDNA/issues/new?...` in the browser | version, Windows / Python / Qt versions and the last 30 log lines, **scrubbed** of sequences, file names, paths and the user name; nothing is sent until the user presses Submit on GitHub | on the menu item |
| Telemetry, analytics, ads | none | none | none |
| Unexpected-error dialog (`app/feedback.py`) | none by itself | the "Report on GitHub" button opens the prefilled issue page as above | on the user's click |

The 3D window says before *Get structure* which service the sequence goes to (`NOTICE`, `structure.py:26-33`). The NCBI dialog has no similar
sentence about the email/key; consider one line in the dialog. `scripts/build_parts.py` also downloads from NCBI but is a developer script,
not part of the program.

## 6. What the program writes, and where

All writes found by searching `app/` for `open(`, `write_text`, `write_bytes`, `os.replace`, `unlink`, `mkdir`, `QSettings`, `tempfile`:

- In the **data folder** (`primer_db.data_dir()`: `%APPDATA%\openDNA\` for the exe, `%APPDATA%\python\` from source, or `OPENDNA_DATA_DIR`):
  `primers.tsv` (+ a transient `primers.tmp`), `parts.tsv` (your part library), `*.ini` (`primers` [PCR conditions, primer options, reference-genome path], `settings`, `theme`, `workspace`, `recent`, `search`, `ncbi`
  [may hold your NCBI API key in plain text], `structure` [the local command text]), `ncbi\` (downloaded records), `structures\` (3D cache and the generated page).
- **Windows registry, file association** (`app/file_assoc.py`, new): only after the user answers Yes to a question, the exe writes `HKCU\Software\Classes\openDNA.SequenceFile`,
  `...\Classes\Applications\openDNA.exe`, one `OpenWithProgids` value per extension, `HKCU\Software\openDNA\Capabilities` and `HKCU\Software\RegisteredApplications\openDNA`. All under the
  user's own hive (no admin rights); it does not touch `UserChoice`. The question is asked once. This is the only code that writes outside the data folder in the
  registry; `PRIVACY.md` lists it and how to remove it.
- **Log files** (`app/feedback.py`): `<data folder>\logs\opendna.log`, 3 x 500 KB, scrubbed.
- **Self-update** (`app/updater.py`): writes the new exe beside the running one, renames the running one to `openDNA.exe.old`, relaunches. That is a write to the program's own folder,
  but only after the user presses the button and the SHA-256 matched. (This replaces the earlier statement that nothing is written in the program's own folder.)
- **Next to the user's file**: `<file>.autosave.gb` recovery copy every minute (deleted on save or clean close; for a `.dna` file the copy goes beside it and the
  `.dna` itself is never written), and Save / Save As where the user chooses. Save is atomic (temp file then rename, `app/primer_db.py:48-52` and the editor).
- **Next to the output chosen in Map Reads**: `<out>_run\` with `ref.fa`, `command.txt`, `log.txt` and the BAM (`app/runner.py:105-118`, `app/main.py:681-683`).
- `%TEMP%\opendna_fold_*`: scratch folder for a Local command, removed afterwards (`structure.py:170`).
- Apart from the self-update above, nothing is written into the program's own folder, and nothing is written outside the places above. Tests write to a temporary data folder (`OPENDNA_DATA_DIR`).
- It runs external programs only on request: the user's local structure command, and `wsl -e bash -lc ...` for read mapping and tool installation (recipes in
  `app/tools.toml`; file paths are shell-quoted, `app/runner.py:55`).
- Opening files: the parsers are Rust (`gb-io`, `noodles`, an XML reader for `.dna`); malformed files are tested to give an error, not a crash
  (`tests/test_roundtrip.py`), but they have not been fuzzed. Mention that in SECURITY.md (done).
