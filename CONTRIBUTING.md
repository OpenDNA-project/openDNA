# Contributing

Thank you for testing OpenGenEdit. It is an early test version made by one person, so reports from real use are the most useful thing you can give.

## Bugs, questions and feature requests

Use **GitHub Issues** with the templates: https://github.com/Ant-Champ/GenomeEdit/issues/new/choose

For a bug, say: the version (Help > About), what you did, what you expected, what happened, and a small file that shows it if you can (never a file
you cannot share). Security problems go through [SECURITY.md](SECURITY.md), not an issue.

## Code changes

Small pull requests are welcome; for a big change, open an issue first so we agree on it. By sending a change you agree that it is licensed under the
same licence as the project, the **GNU AGPL v3** ([LICENSE](LICENSE)). Do not add files you do not have the right to share (real sequence files, other
people's data, a library without a compatible licence: add new third-party pieces to [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md)). If a new panel is
based on a paper or a dataset, add it to Help > Citations (`app/citations.py`).

**Your own fork** is welcome too, under the same licence; give it a name and an icon of its own, as [TRADEMARKS.md](TRADEMARKS.md)
explains (the names "OpenGenEdit" and "oDNA" stay with this project).

### Run the tests

You need Rust, the Visual Studio C++ Build Tools and a conda environment (not `base`; see the README).

```powershell
conda create -n genomeedit python=3.13
conda activate genomeedit
pip install -r requirements-dev.txt
maturin develop -m core/Cargo.toml --release     # rebuild after every change in core/
pytest -q tests
```

The Rust tests need to find Python:

```powershell
$env:PYO3_PYTHON = (Get-Command python).Source   # the genomeedit environment's python.exe
cd core; cargo test; cd ..
```

Some tests use your own `.dna` / `.gb` files in the repository root (not in git) and skip themselves when the files are absent, so a clean clone should be
green. [docs/TESTING.md](docs/TESTING.md) lists what is covered and what is not; a fix should come with a test, and the test plan should be updated.
CI runs `cargo test`, `pyflakes app tests` and `pytest` on every pull request.

### Style (what the code around it does)

- Comments explain **why**, not what the next line does. Look at a neighbouring file and match it.
- **Rust is commented heavily**, on purpose: the author does not write Rust, so each non-obvious construct is explained in plain words (see
  [docs/RUST_PRIMER.md](docs/RUST_PRIMER.md)). Keep that when you change `core/`. Keep the logic in plain Rust modules; `core/src/lib.rs` is the only
  file that knows about Python.
- Python: PySide6, no extra dependencies without a reason. Reuse what exists: new tools embed the existing sequence view and open as workspace tabs,
  they do not get a copy of it.
- Anything that talks to the internet must be listed in [PRIVACY.md](PRIVACY.md) and tell the user what is sent before it is sent.
- `core/vendor/gb-io` is a patched copy of a third-party crate: mark changes "OpenGenEdit patch".
- No emojis in code, docs or messages.
