# Working on openDNA

## Python: use the repo-local `.venv`

The tests need PySide6 **and** `opendna_core`, the compiled Rust core. Only the `.venv` at the repo
root has both (Python 3.12, created with `uv`). On this machine it is:

```
D:\GitHub_Projects\openDNA\.venv\Scripts\python.exe
```

Call that interpreter by its full path (from the repo root: `.venv\Scripts\python.exe`). A bare
`python` on PATH is the plain uv-managed 3.12 with none of the project's packages, so imports fail
in ways that look like broken code rather than a missing environment.

```bash
# tests (what CI runs)
.venv\Scripts\python.exe -m pytest -q tests

# one file, or one test
.venv\Scripts\python.exe -m pytest -q tests/test_assembly.py

# lint (also in CI)
.venv\Scripts\python.exe -m pyflakes app tests

# the Rust side (from core/, as CI does); PYO3_PYTHON must name the venv's interpreter
cd core && PYO3_PYTHON="D:\GitHub_Projects\openDNA\.venv\Scripts\python.exe" cargo test
```

The suite is about 95 tests of logic, one role each, grouped by area (no UI tests: nothing builds a
window, dialog or widget). It takes about 35 seconds. Keep it that way: extend the test whose role
fits rather than adding a new one, and do not add tests that drive widgets.

`scripts/test_changed.py` runs test files in parallel and needs `pytest-xdist`, which is installed
in `.venv`.

### Recreating the environment

```bash
uv venv --python 3.12 .venv
uv pip install --python .venv\Scripts\python.exe -r requirements.txt pytest pytest-timeout pytest-xdist pyflakes biopython
```

Then rebuild the core (below). Do not run `uv init` in the repo: it creates a stub `pyproject.toml`
pinned to a newer Python.

### Rebuilding `opendna_core`

It is compiled from `core/` into the environment, so after any change under `core/src` it must be
rebuilt or Python keeps importing the old module and the change appears to do nothing:

```bash
.venv\Scripts\python.exe -m maturin develop -m core/Cargo.toml --release
```

`maturin develop` insists on a virtualenv: set `VIRTUAL_ENV` to the `.venv` path first (that is all
`openDNA-test.bat` does before calling it). CI builds a wheel instead — see
`.github/workflows/ci.yml`.

`D:\GitHub_Projects\openDNA-test.bat` — one directory **above** the repo, so not version-controlled —
is the dev launcher: it rebuilds the core and opens the app on `demo/BW25113_2026.gb`. It locates the
repo and `.venv` relative to itself; keep it and this file in agreement.

## The test report: keep it in step with the tests

`docs/test-report/report.html` is the owner's way of checking that every test is relevant. It lists
each test (Python and Rust) with why it matters, its inputs, the output it expects, what that answer
is checked against, a strength rating and review notes. The text comes from
`docs/test-report/tests.json`; `scripts/test_report.py` builds the HTML page and the Markdown pages
beside it from that file plus the latest run.

Whenever you add, rename, merge, delete or substantially change a test (Python `def test_*` or Rust
`#[test]`), in the same change:

1. Edit its entry in `tests.json` (`id` is `tests/<file>.py::<name>` or `core/src/<file>.rs::tests::<name>`).
   Fields: `area`; `role` (one sentence: what would break unnoticed for a user without this test);
   `purpose`; `cases` (list of `{"input", "expected"}`, concrete values, in the order the test runs
   them); `oracle` (one of `hand-computed`, `independent implementation`, `in-test reference model`,
   `property/invariant`, `round-trip/self-consistency`, `regression value`, `smoke`); `strength`
   (`strong`, `adequate`, `weak`; be honest); `issues` (review notes, may be empty).
   Remove the entries of deleted tests.
2. Rebuild the report and commit the generated files with the change:

```bash
QT_QPA_PLATFORM=offscreen .venv\Scripts\python.exe -m pytest -q tests --junitxml=junit.xml
(cd core && PYO3_PYTHON="D:\GitHub_Projects\openDNA\.venv\Scripts\python.exe" cargo test > ../cargo.txt 2>&1)
.venv\Scripts\python.exe scripts/test_report.py --junit junit.xml --cargo cargo.txt
```

The script prints how many tests have a written description and lists entries whose test no
longer exists: both must be clean (every test, Python and Rust, described; no stale entries). It
does not delete the Markdown page of a removed test file: delete `docs/test-report/tests_<file>.md`
yourself. Do not commit `junit.xml` or `cargo.txt`.

`docs/test-report/AUDIT.md` (the report's "Audit" tab) is generated too: every test on one page with
its role, reference, strength and review notes. Never edit it by hand; change `tests.json` and
rebuild. When a review note in `tests.json` names a test that has since been renamed or merged,
update the note.

## A Qt failure that is the version, not the code

The assembly panel's Actions menu works on PySide6 6.12 (what `.venv` has) and failed on 6.8.3 with
`RuntimeError: Internal C++ object (QMenu) already deleted` — a garbage-collection difference in Qt,
not a bug in the app. `requirements.txt` only asks for `>=6.8`, so a fresh install can land on a
version that fails. Check the PySide6 version before chasing that error in the code.

## Before pushing a PR

Run `/simplify` and then `/code-review` on the change before pushing any PR, and act on what they find.

## Agent skills

### Issue tracker

Issues and specs are GitHub issues on `OpenDNA-project/openDNA`, via the `gh` CLI. See `docs/agents/issue-tracker.md`.

### Triage labels

The five default labels: `needs-triage`, `needs-info`, `ready-for-agent`, `ready-for-human`, `wontfix`. See `docs/agents/triage-labels.md`.

### Domain docs

Single-context: one `GLOSSARY.md` and `docs/adr/` at the repo root, created when first needed. See `docs/agents/domain.md`.
