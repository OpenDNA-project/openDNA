# Working on openDNA

## Python: use the `opendna` conda environment

The tests need PySide6 **and** `opendna_core`, the compiled Rust core. Only the `opendna` conda
environment has both. On this machine it is:

```
C:\Users\Ant\anaconda3\envs\opendna\python.exe
```

Call that interpreter by its full path. Do **not** use a bare `python`: the first one on PATH is
Inkscape's bundled interpreter (no pytest), and the base Anaconda has pytest but neither PySide6 nor
`opendna_core`. Both fail in ways that look like broken code rather than a missing environment.

```bash
# tests (what CI runs)
"C:\Users\Ant\anaconda3\envs\opendna\python.exe" -m pytest -q tests

# one file, or one test
"C:\Users\Ant\anaconda3\envs\opendna\python.exe" -m pytest -q tests/test_assembly.py

# lint (also in CI)
"C:\Users\Ant\anaconda3\envs\opendna\python.exe" -m pyflakes app tests

# the Rust side (from core/, as CI does)
cd core && cargo test
```

Qt tests open real windows. Set `QT_QPA_PLATFORM=offscreen` to keep them off the screen.

`scripts/test_changed.py` runs test files in parallel and needs `pytest-xdist`, which is **not**
installed in the environment: `python -m pip install pytest-xdist` before using it.

### Rebuilding `opendna_core`

It is compiled from `core/` into the environment, so after any change under `core/src` it must be
rebuilt or Python keeps importing the old module and the change appears to do nothing:

```bash
"C:\Users\Ant\anaconda3\envs\opendna\python.exe" -m maturin develop -m core/Cargo.toml --release
```

`maturin develop` insists on a virtualenv. The environment is a conda env, so set `VIRTUAL_ENV` to
its path first (that is all `openDNA-test.bat` does before calling it). CI builds a wheel instead —
see `.github/workflows/ci.yml`.

`C:\Users\Ant\code\openDNA-test.bat` — one directory **above** the repo, so not version-controlled —
is the dev launcher: it rebuilds the core and opens the app on `demo/BW25113_2026.gb`. It is the
reference for which interpreter is the dev one; keep it and this file in agreement.

## A Qt failure that is the version, not the code

`tests/test_assembly.py::test_actions_menu_pointer_product_box_and_sidewinder_panel` passes on
PySide6 6.12 (what the `opendna` environment has) and fails on 6.8.3 with
`RuntimeError: Internal C++ object (QMenu) already deleted` — a garbage-collection difference in Qt,
not a bug in the app. `requirements.txt` only asks for `>=6.8`, so a fresh install can land on a
version that fails. Check the PySide6 version before chasing that error in the code.
