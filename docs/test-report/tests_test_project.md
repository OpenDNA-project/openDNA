# `tests/test_project.py`

[← all test files](README.md) · source: [`tests/test_project.py`](../../tests/test_project.py)

## `test_release_metadata_and_exe_recipe_are_consistent`

**Result:** ✅ passed in 0.01 s · [line 23](../../tests/test_project.py#L23) · **Area:** Release packaging · **Checked against:** property/invariant · **Strength:** 🟡 adequate

**Why it matters:** Without it, a release could ship with mismatched version numbers (breaking the updater) or a Windows exe missing data files the app needs (parts library, tool settings, 3D viewer), which only shows up for users of the frozen build.

Checks project.py's version equals core/Cargo.toml and core/pyproject.toml and is major.minor.patch; URLs derive from REPO; name, author, update switch and file-type extensions are sane; every Path(__file__).with_name(...) data file is bundled by packaging/openDNA.spec and the spec's files exist; the icon is a valid multi-size .ico.

| # | Input | Expected output |
|---|---|---|
| 1 | project.VERSION, Cargo.toml version, pyproject.toml version | all equal and match \d+.\d+.\d+ |
| 2 | REPO and derived URLs | owner/name; REPO_URL = https://github.com/REPO; ISSUES_URL, RELEASES_URL, LATEST_RELEASE_API built from it; all https, no spaces |
| 3 | NAME, AUTHOR, UPDATES_ENABLED | NAME == 'openDNA'; AUTHOR not blank; UPDATES_ENABLED is bool |
| 4 | FILE_TYPES | non-empty; each key starts with '.', lower case, no spaces, label non-blank |
| 5 | Scan app/*.py for Path(__file__).with_name("x") | finds at least parts, sponsors.json, tools.toml, vendor; each name appears quoted in the spec; spec bundles vendor/3Dmol-min.js and its LICENSE |
| 6 | Paths named in the spec as APP/... or ROOT/...; fixed list of data files | all exist |
| 7 | packaging/openDNA.ico header and directory | type 1 icon; each entry square, 32 bpp, inside the file; sizes include 16, 32, 48, 256 |

**Review notes:**

- Not a test of DNA editing; it is a release-hygiene check. Relevant for packaging but could be marked as such.
- The bundling check is text matching on the spec ('"name"' appears anywhere), so a name mentioned but not actually bundled would pass.
- Only catches data files reached through the exact with_name("...") pattern.

## `test_old_data_folder_is_copied_once`

**Result:** ✅ passed in 0.01 s · [line 84](../../tests/test_project.py#L84) · **Area:** Data folder migration · **Checked against:** hand-computed · **Strength:** 🟡 adequate

**Why it matters:** Without it, users upgrading from the program's old name could lose their primers and settings, or a later start could overwrite their new data with the old copy.

Checks copy_old_data copies the old-name data folder to the new location when the new one does not exist, leaves the old folder in place, never copies over an existing new folder, and returns False when there is no old folder.

| # | Input | Expected output |
|---|---|---|
| 1 | Old folder with primers.tsv and logs/, no new folder | returns True; new primers.tsv has the same text and logs/ exists; old primers.tsv still there |
| 2 | New primers.tsv changed, call again | returns False; new file still 'changed' |
| 3 | Target with no old sibling folder | returns False |

**Review notes:**

- Lives in test_project.py (release metadata) though it is about primer_db's data migration.
- The copy-error path (half copy returns False) is not tested.
