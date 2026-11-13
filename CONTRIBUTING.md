# Contributing

Thank you for testing openDNA. It is an early test version made by one person, so reports from real use are the most useful thing you can give.

## Bugs, questions and feature requests

Use **GitHub Issues** with the templates: https://github.com/OpenDNA-project/openDNA/issues/new/choose

For a bug, say: the version (Help > About), what you did, what you expected, what happened, and a small file that shows it if you can (never a file
you cannot share). Security problems go through [SECURITY.md](SECURITY.md), not an issue.

## Code changes

Small pull requests are welcome; for a big change, open an issue first so we agree on it. Read **Copyright of contributions** below
first: sending a change means you accept it. Do not add files you do not have the right to share (real sequence files, other
people's data, a library without a compatible licence: add new third-party pieces to [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md)). If a new panel is
based on a paper or a dataset, add it to Help > Citations (`app/citations.py`).

### Copyright of contributions

openDNA has one copyright holder, Antoine Champie (the project owner), so that the project can change its licence later without
having to find and ask every contributor. By submitting a contribution (a pull request, a patch, or code, text, images or data sent
in any other way), you:

1. **assign** to the project owner all copyright and related rights in your contribution, worldwide, for the full term of those
   rights and at no charge;
2. **waive**, or where a waiver is not possible agree not to assert, any right that cannot be assigned (such as moral rights in some
   countries) against the project owner and anyone who receives openDNA from the project owner;
3. agree that the project owner may publish, license and **relicense** openDNA, your contribution included, under any terms,
   including a licence other than the AGPL, without asking you or telling you first;
4. confirm that the contribution is your own work and that you have the right to make this assignment. If an employer or a
   university may own what you write, get its written permission first. Code or data by others may be included only under a
   licence compatible with the project's, named in [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).

Your contribution is still released under the project's licence, the **GNU AGPL v3** ([LICENSE](LICENSE)), so you keep, like
everyone, the rights that licence gives. If you do not agree with these terms, do not send the contribution. The pull-request
template asks you to confirm them.

**Unfinished features** carry an orange **WIP** tag. `FEATURES` in `app/wip.py` is the switch: one line per feature, `True` while it
is a work in progress. Tag each place the feature shows with the helpers of that file (`wip.Badge("key")` beside a title or a choice,
`wip.mark_action(action, "key")` on a menu item, `wip.Delegate` on a list); when the feature is done, set its line to `False` and
every tag of it goes. `SHOW_TAGS = False` hides them all.

**Your own fork** is welcome too, under the same licence; give it a name and an icon of its own, as [TRADEMARKS.md](TRADEMARKS.md)
explains (the names "openDNA" and "oDNA" stay with this project).

### Run from source

Windows is the primary platform. Install [Rust](https://rustup.rs/), the Visual Studio C++ Build Tools and
[uv](https://docs.astral.sh/uv/), then, from the repository root:

```powershell
uv venv --python 3.12 .venv
.venv\Scripts\activate
uv pip install -r requirements.txt
maturin develop -m core/Cargo.toml --release
python app/main.py demo/pBAD30_eilR.gb
```

Use a dedicated virtual environment such as this `.venv`, not a Python that other tools share (a conda `base`, for one, puts its own
ICU 73 on the DLL path; PySide6's `Qt6Core.dll` then loads it instead of Windows' ICU and fails with "DLL load failed … The
specified procedure could not be found"). Any Python 3.10 or newer with `venv` works in place of uv.

The sequences in `demo/` are real: two plasmids, the *E. coli* BW25113 genome and *Vibrio cholerae* N16961 (both chromosomes).
`python scripts/make_demo_files.py` writes the three small made-up fragments used by the tutorial.

### Build the exe

A standalone `dist\openDNA.exe` (one file, about 39 MB, runs on a PC without Python): after
`pip install -r requirements-dev.txt` and `maturin develop -m core/Cargo.toml --release`, run

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File scripts/build_exe.ps1
```

It builds from the committed recipe `packaging/openDNA.spec`, prints the size and SHA-256, and smoke-tests the exe: it starts it,
opens a demo file and checks that it stays alive (`-NoSmoke` skips that). It also builds `dist\openDNA-folder.zip`, the same program
as a folder with the Qt files replaceable (`-NoFolder` skips it). Releases are built from a `v*` tag by
`.github/workflows/release.yml` into a draft GitHub release, under the name in `app/project.py` (`EXE_NAME`); see
[docs/PUBLISHING.md](docs/PUBLISHING.md).

### The README pictures

`python scripts/make_screenshots.py` makes every picture in `docs/images/` again, from the real program and the real sequences in
`demo/`. Run it with real fonts (do not set `QT_QPA_PLATFORM`); its settings go to a temporary folder, never yours. Name one or
more pictures to redo only those (`editor`, `map`, `gibson`, `primers`, `cas9`, `align`, `tracks`). The `tracks` picture maps
simulated reads, so it needs bwa and samtools the way *Tools > Map Reads* does. Redo the pictures when you change the interface.

### Run the tests

You need Rust, the Visual Studio C++ Build Tools and a virtual environment (see above).

```powershell
uv venv --python 3.12 .venv
.venv\Scripts\activate
uv pip install -r requirements-dev.txt
maturin develop -m core/Cargo.toml --release     # rebuild after every change in core/
pytest -q tests
```

The Rust tests need to find Python:

```powershell
$env:PYO3_PYTHON = (Get-Command python).Source   # the .venv's python.exe
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
- `core/vendor/gb-io` is a patched copy of a third-party crate: mark changes "openDNA patch".
- No emojis in code, docs or messages.
