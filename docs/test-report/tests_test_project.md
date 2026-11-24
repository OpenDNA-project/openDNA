# `tests/test_project.py`

[← all test files](README.md) · source: [`tests/test_project.py`](../../tests/test_project.py)

## `test_release_metadata_and_exe_recipe_are_consistent`

**Result:** ✅ passed in 0.95 s · [line 23](../../tests/test_project.py#L23) · **Area:** Release packaging · **Checked against:** property/invariant · **Strength:** 🟡 adequate

**Why it matters:** Without it, a release could ship with mismatched version numbers (breaking the updater) or a Windows exe missing data files the app needs (parts library, tool settings, 3D viewer), which only shows up for users of the frozen build; or a missing or malformed changelog would reach the start page and the release description.

Checks project.py's version equals core/Cargo.toml and core/pyproject.toml and is major.minor.patch; URLs derive from REPO; name, author, update switch and file-type extensions are sane; every Path(__file__).with_name(...) data file is bundled by packaging/openDNA.spec and the spec's files exist; the icon is a valid multi-size .ico. Also checks CHANGELOG.md: the first section is Unreleased, a section exists for project.VERSION, every other section is '## x.y.z (date)' with Features, Bugfixes, UI/UX labels in that order, each with '- ' lines, and at most 3 Features lines.

| # | Input | Expected output |
|---|---|---|
| 1 | project.VERSION, Cargo.toml version, pyproject.toml version | all equal and match \d+.\d+.\d+ |
| 2 | REPO and derived URLs | owner/name; REPO_URL = https://github.com/REPO; ISSUES_URL, RELEASES_URL, LATEST_RELEASE_API built from it; all https, no spaces |
| 3 | NAME, AUTHOR, UPDATES_ENABLED | NAME == 'openDNA'; AUTHOR not blank; UPDATES_ENABLED is bool |
| 4 | FILE_TYPES | non-empty; each key starts with '.', lower case, no spaces, label non-blank |
| 5 | Scan app/*.py for Path(__file__).with_name("x") | finds at least parts, sponsors.json, tools.toml, vendor; each name appears quoted in the spec; spec bundles vendor/3Dmol-min.js and its LICENSE |
| 6 | Paths named in the spec as APP/... or ROOT/...; fixed list of data files | all exist |
| 7 | packaging/openDNA.ico header and directory | type 1 icon; each entry square, 32 bpp, inside the file; sizes include 16, 32, 48, 256 |
| 8 | CHANGELOG.md headings and section bodies | first heading Unreleased; project.VERSION has a section; others '## x.y.z (YYYY-MM-DD)'; labels from Features/Bugfixes/UI/UX in order, each with '- ' lines; Features at most 3 lines |

**Review notes:**

- Not a test of DNA editing; it is a release-hygiene check. Relevant for packaging but could be marked as such.
- The bundling check is text matching on the spec ('"name"' appears anywhere), so a name mentioned but not actually bundled would pass.
- Only catches data files reached through the exact with_name("...") pattern.
- The changelog check is about shape only: it cannot know that a line is true or that a user-visible change was written down.

**Fixed by the audit of November 2026:**

- 13 Nov 2026: gap closed. The proposed checks were added, pass on the current code and fail under the mutation that slipped through before. Was: 'every data file the app reads is bundled by the spec' is a substring test ('"name"' in SPEC): a name in a comment, a hiddenimports list or an excludes list satisfies it, so a frozen exe missing tools.toml would ship.

## `test_old_data_folder_is_copied_once`

**Result:** ✅ passed in 0.01 s · [line 109](../../tests/test_project.py#L109) · **Area:** Data folder migration · **Checked against:** hand-computed · **Strength:** 🟡 adequate

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

## `test_welcome_slots_hold_four_and_expire`

**Result:** ✅ passed in 0.03 s · [line 123](../../tests/test_project.py#L123) · **Area:** Welcome page · **Checked against:** hand-computed · **Strength:** 🟡 adequate

**Why it matters:** Without it, the naming vote on the welcome page could ship with a broken link, a missing picture or a fifth slot silently dropped, or stay up after its end date; or the page could open at every start (or never) instead of once a day.

Checks sponsors.load returns four filled slots with an https link, a button, a 'sends' text and existing pictures, empties slots after their 'until' date, and keeps the old behaviour (padded to three, never expiring, SPONSORED heading) for a file without those keys; then that the page is due at startup once a day by default and at every start when the user asks.

| # | Input | Expected output |
|---|---|---|
| 1 | load on 2026-12-31 (the last day) | 4 slots, each with name, https url without spaces, button, sends, and image files that exist |
| 2 | load on 2027-01-01 | [{}, {}, {}, {}] |
| 3 | A file with one slot {name: A} and no heading | ([{name: A}, {}, {}], '', 'SPONSORED') |
| 4 | due() with no saved state; then last_shown = today; then every_start = True | True (every_start False); False; True |

**Review notes:**

- Does not check the card drawing or the click on the picture (no widget tests); the page was only rendered once offscreen by hand.
- Does not detect a placeholder form link (FORM_ID) left in sponsors.json.
