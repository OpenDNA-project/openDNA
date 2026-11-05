# `tests/test_project.py`

[← all test files](README.md) · source: [`tests/test_project.py`](../../tests/test_project.py)

## `test_project_metadata_is_consistent`

**Result:** ✅ passed in 0.00 s · [line 19](../../tests/test_project.py#L19) · **Area:** Project metadata / release · **Checked against:** property/invariant · **Strength:** 🟢 strong

A release must bump the version in app/project.py, core/Cargo.toml and core/pyproject.toml together, and links, name, author, update switch and file-type list must be well-formed.

| # | Input | Expected output |
|---|---|---|
| 1 | project.VERSION, Cargo.toml package.version, pyproject.toml project.version | all equal; matches r'\d+\.\d+\.\d+' |
| 2 | project.REPO and URLs | owner/name; REPO_URL == 'https://github.com/'+REPO; ISSUES/RELEASES = +'/issues', +'/releases'; API 'https://api.github.com/repos/<REPO>/releases/latest'; all https, no spaces |
| 3 | NAME, AUTHOR, UPDATES_ENABLED | NAME == 'openDNA'; AUTHOR non-blank; UPDATES_ENABLED is bool |
| 4 | FILE_TYPES | non-empty; each key starts with '.', lowercase, no space, length>1; label non-blank |

## `test_exe_recipe_bundles_every_data_file_the_app_reads`

**Result:** ✅ passed in 0.01 s · [line 50](../../tests/test_project.py#L50) · **Area:** Packaging (PyInstaller spec) · **Checked against:** property/invariant · **Strength:** 🟡 adequate

Every data file the app reads via Path(__file__).with_name(...) must be bundled by packaging/openDNA.spec, the files the spec names must exist, and the icon must be a multi-size .ico.

| # | Input | Expected output |
|---|---|---|
| 1 | Regex scan of app/*.py for Path(__file__).with_name("x") | finds at least parts, sponsors.json, tools.toml, vendor; each name appears quoted in the spec |
| 2 | Spec text | contains APP / "vendor" / "3Dmol-min.js" and 3Dmol.LICENSE |
| 3 | Every (APP / "…") and (ROOT / "…") path in the spec; 6 listed files | exists under app/ or repo root |
| 4 | packaging/openDNA.ico header | reserved 0, type 1; each entry square, 32 bpp, data inside file; sizes include 16, 32, 48, 256 |

**Review notes:**

- L60: membership is a substring test ('"name"' anywhere in the spec), so a name mentioned in a spec comment or another path would satisfy it without being bundled.
- L57: only the exact double-quoted with_name("x") form is scanned; other ways of reaching data files would be missed (none exist today besides frozen-aware primer_panel.py:35).
