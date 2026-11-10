# Phase 7 — what changed, and what a reviewer must check

> **Merged.** Everything below is on `main` as of `c687874`: the seven Phase 7 branches, the file-writing fix, and the
> alignment, primer-origin and protein-view work that was still loose in the working tree. Measured on merged `main`:
> **417 passed, 2 skipped**, pyflakes clean, `cargo test` 102 passed, and the app opens the 4.6 Mb demo genome with no
> missing bundled files. The merge turned up two collisions no branch could see on its own — Alt+R claimed by both the
> Repeats and the Recode tools, and a folder-search test whose "random" query was `demo/msa_demo.fasta`'s first 60
> bases, both generated from `random.seed(7)` — fixed in `c687874`. The `claude/*` branches and their worktrees are
> still there; nothing is pushed.

Phase 7 of [PLAN.md](PLAN.md) is being built by several Claude sessions, each in its own git worktree under
`.claude/worktrees/` and on its own branch. This file is the ledger: what each one touched, what is still unverified, and
what a reviewer should check before any of it is merged. **Keep it current** — a session that reports back gets its row
filled in, and anything it says it did not finish goes in the "To verify" column, not into a commit message where it will
be missed.

Nothing here has been merged into `feature/settings-and-codons` or `main`.

## Running order

One session at a time, by the owner's instruction. `degenerate-primers` holds the shared compiled core, so it goes
first; the rest were stopped mid-turn and stay open, to be resumed in this order, each only once the one before is idle:

1. **degenerate-primers** — resumed; holds the Rust core and has nothing committed yet.
2. **restriction-cloning** — stopped mid-turn, nothing committed.
3. **msa-view** — stopped mid-turn, nothing committed.
4. **library-search** — stopped mid-turn, nothing committed.
5. **codon-tools** and **repeats-dotplot** — committed already; they need re-verification against a core built from
   their own branch, not more feature work.
6. **oligo-order** — committed already; needs one clean full-suite run.

A stopped session keeps its worktree and its uncommitted edits, but it was interrupted, so its files may be half-written:
the first thing each resumed session should do is re-read what it had changed and re-run its own tests.

## Verification pass (in progress)

Run from the main checkout, one branch at a time, nothing else building. Only two distinct `core/` trees exist among the
committed branches, so two rebuilds cover all three:

| core tree | built from | verifies |
|---|---|---|
| `bd486677` (base `b40dac3`, feature/settings-and-codons) | `repeats-dotplot` worktree | codon-tools, repeats-dotplot |
| `785a57da` (base `e5ab2df`, main) | `oligo-order-export` worktree | oligo-order |

`degenerate-primers` has its own, third core (uncommitted changes to `core/src/primers.rs`, `lib.rs`, `search.rs`), so it
needs its own build when its turn comes.

### Results

| Branch | Python suite | pyflakes | cargo test | Verdict |
|---|---|---|---|---|
| `codon-tools` (`a695d6a`) | **377 passed, 2 skipped** | clean | 96 passed | verified |
| `repeats-dotplot` (`f528b72`) | **363 passed, 2 skipped** | clean | 96 passed | verified — the 2 failures it reported were the wrong core, as it suspected |
| `degenerate-primers` (`dab7f11`) | **326 passed, 2 skipped** with `tests/test_record_list.py` ignored (its base `main` carries 50 test files, not 52) | clean | 99 passed | verified, as far as its base allows |
| `restriction-cloning` (`587ea13`) | **367 passed, 2 skipped** | clean | — | verified; its changes to the three existing assembly modes still need a human read |
| `msa-view` (`2c5c1f4`) | **362 passed, 2 skipped** | clean | — | verified; **only the built-in aligner was ever run** — the MAFFT/Clustal recipes in `tools.toml` are untested, and the view has never been opened in the real app |
| `library-search` (`db64bde`) | **362 passed, 2 skipped** | clean | — | verified; it renamed a menu item after its own last full run, so this run is the only complete one on that commit |
| `oligo-order-export` (`453402b`) | aborts on its own base; **364 passed, 2 skipped** when applied to `feature/settings-and-codons` | clean | 95 passed | code verified; the branch's base is broken, see below |

### `main` hung in the record-list test — nothing to do with Phase 7, and now fixed

The oligo branch's abort traces to `tests/test_record_list.py::test_a_line_dragged_in_the_map_reorders_the_file`, which
hangs inside `app/main.py` `_write_record` → `save` (the Rust write) until `pytest-timeout` fires and the process dies
with a Windows access violation.

**Checked on `e5ab2df` (`main`) itself, with no Phase 7 code at all: the same test hangs.** It ran past 400 s and was
killed. The branch before it in the same file passes. So this is a bug on `main`, which the oligo branch merely
inherited by being cut from there; the two branches based on `feature/settings-and-codons` run the same test fine.

It is worse than one test. With only `test_a_line_dragged_in_the_map_reorders_the_file` deselected, the suite still dies
with an access violation, and the dump holds no frames from `app/` or `tests/` at all — it crashes at interpreter level,
after or between tests. The whole file has to be left out (`--ignore=tests/test_record_list.py`) for a run on a
`main`-based branch to complete. So at least two tests in that file leave the process in a state it cannot survive.

Test counts differ by base, which explains why `main`-based branches report fewer: `main` carries 50 test files,
`feature/settings-and-codons` 52.

**Diagnosed and fixed** on `claude/safe-file-writes` (`c1821aa`, cut from `main`). The save was not hanging in the Rust
write: it was failing. On Windows a rename over a file is refused while any thread still has that file open, and the app
has such readers — the background search for primer sites reads the other records of the open file. `_write_record`
caught the error and opened `QMessageBox.critical`, a modal box that under `QT_QPA_PLATFORM=offscreen` nobody can ever
close. That is the hang, and the access violation is Qt being torn down with a modal dialog still up.

The fix already existed on `feature/settings-and-codons` and was never on `main`: `app/fileio.py` retries the rename for
~5 s while the target is a writable file and raises at once for anything else, `core/src/document.rs` does the same for
the Rust save, and the four call sites (`main.py`, `multigb.py` ×3, `primer_db.py`) go through it. Picked onto `main` as
one commit, with `tests/test_fileio.py`.

Measured on `c1821aa`: **329 passed, 2 skipped** (the whole suite, nothing ignored), pyflakes clean, `cargo test` 95
passed. `main` had not completed a suite run at all before this.

This was a real user-facing bug, not a test artefact: saving while a primer search was running failed outright, with a
"Save Failed" box, on every GenBank write the app makes.

To judge the export code on its own, `453402b` was applied (uncommitted) onto a worktree at `feature/settings-and-codons`
and the suite run there: **364 passed, 2 skipped, pyflakes clean** (357 on that base plus its 7 new tests). The export
work is sound. What is not sound is the base it was committed onto.

**So `claude/oligo-order-export-538e2b` should be replayed onto `feature/settings-and-codons` before it is merged
anywhere** — a cherry-pick of `453402b` applies cleanly, as the run above proves. Doing it needs a git identity: this
repository has none set, which is also why that commit is authored "Ant" with the owner's email.

## What no test has covered

Every Phase 7 result above comes from offscreen tests. **Not one of these features has been seen in the running app.**
For a tool whose value is what a biologist sees on screen — a dot plot, an MSA view, four assembly cards that should
read alike — green suites do not answer the question that matters. Open each one in the real window before merging.

Also untested: the MAFFT and Clustal Omega recipes (neither tool is installed on this machine), and every judgment
change made to code that already worked — `587ea13`'s rework of the Gibson, Golden Gate and Sidewinder modes above all.

## All seven are committed and verified — what is left

1. ~~**Fix the crash on `main`.**~~ Done: `claude/safe-file-writes` (`c1821aa`), 329 passed on `main` with nothing
   ignored. The two branches cut from `main` (`degenerate-primers`, `oligo-order`) can show a green suite once they sit
   on top of it. Not merged into `main`.
2. **Reconcile the two bases.** Five branches sit on `feature/settings-and-codons`, two on `main`. Replaying the two
   onto the feature branch is the simpler direction; `453402b` has already been proved to cherry-pick cleanly.
3. **Set a git identity** in this repository. There is none configured, so authorship is inconsistent: some commits
   are "Ant" with the owner's email (`a695d6a`, `453402b`), others "Claude <noreply@anthropic.com>" (`dab7f11`,
   `c1821aa`, and every commit on `main` and `feature/settings-and-codons`).
4. **Look at all of it in the running app**, which no test has done.
5. **Read the judgment changes**, above all `587ea13`'s rework of the three existing assembly modes.

## Branches

| Session | Branch | Item | State | Head |
|---|---|---|---|---|
| oligo-order | `claude/oligo-order-export-538e2b` | I — oligo order export | committed, suite never completed | `453402b` |
| degenerate-primers | `claude/degenerate-primers-4a570f` | A — IUPAC primers | **stopped mid-turn, nothing committed** — its work exists only as uncommitted files in its worktree | `e5ab2df` (= `main`) |
| codon-tools | `claude/codon-tools-f1cbb6` | E — codon optimisation, cross-code translation | committed; reported 377 passed / 2 skipped, but against its own private core build | `a695d6a` |
| restriction-cloning | `claude/restriction-ligation-60c032` | F — restriction–ligation + assembly UX | committed, verified | `587ea13` |
| msa-view | `claude/msa-view-ff687e` | G — multiple sequence alignment view | committed, verified | `2c5c1f4` |
| library-search | `claude/library-search-1ded43` | H — search a folder of plasmids | committed, verified | `db64bde` |
| repeats-dotplot | `claude/repeats-dotplot-c01ecf` | J — repeat finder and dot plot | committed; 2 suite failures, attributed to a stale core | `f528b72` |

## Every Python test result so far is void

The seven worktrees share one conda environment, and `opendna_core` is a single compiled module in its
`site-packages` — not an editable link, a built copy that each `maturin develop` **overwrites for every worktree at
once**. They ran concurrently, so they overwrote each other's core while each other's suites were running.

Measured directly: the module installed right now exports `tm_range` and finds binding sites for a primer containing
`N`, so it is **degenerate-primers' build**, and it does *not* contain the `align.rs` end-trimming fix from the main
checkout. That single fact explains all three oddities the sessions reported — the oligo suite crashing inside the test
loop, codon-tools' "private core build" caveat, and repeats-dotplot's two failures in
`test_codon_usage.py` and `test_primer_binding.py`, neither of which it had touched.

So: **no branch has a trustworthy full-suite result yet, including the main checkout.** Each branch must be verified on
its own, one at a time, with `maturin develop` run from that branch's worktree immediately before its suite, and nothing
else building meanwhile. Until then, treat the sessions' green reports as "it compiled and their own new tests passed".

Bases differ: `degenerate-primers` was cut from `main` (`e5ab2df`), everything else from `feature/settings-and-codons`
(`b40dac3`), which carries 7 further commits (codon editor, protein view, demo genome). A reviewer comparing two branches
must diff against the right base, or the codon and protein-view work will look like someone's deletion.

## Files changed, per branch

### `claude/oligo-order-export-538e2b` — committed as `453402b`

`app/primer_db.py` (+43), `app/primer_panel.py` (+58), `app/main.py` (+8), `app/db_tabs.py` (+5),
`docs/USER_GUIDE.md` (+17), `tests/test_primer_order_export.py` (new, 66 lines).

**To verify:** the session committed before the full Python suite finished and stopped it, so only its own 7 new tests are
known to pass; `pyflakes` was clean. A full run of the suite on this branch **has not yet produced a trustworthy result**:
the attempt exited 0 but printed a faulthandler stack from inside `pytest_runtestloop`, most likely because another
worktree rebuilt the shared `opendna_core` underneath it. Re-run it when no other session is building, keep the whole
output, and do not accept exit code 0 alone as a pass — read the `N passed` line (`conftest.py` deliberately ends the
process by hand after the summary, so a stack trace near the end is normal and a crash *during* the run is not). Also check the plate well order (down the
columns, A1→H1 then A2), the roll-over past 96 oligos, and that a primer name holding a comma or a quote survives the
CSV. The commit was authored as "Ant" with the user's email because git had no identity set in that worktree — the user
may want it amended.

### The other six

Working; rows to be filled in as each reports. Expected areas, from the briefs:

- **degenerate-primers** — `core/src/primers.rs` (IUPAC matching in `binding_sites`, new `tm_range`), `core/src/lib.rs`,
  `app/primer_panel.py`, possibly `app/primer_db.py` and `app/primer_design.py`. The only session allowed to rebuild the
  compiled core.
- **codon-tools** — new `app/codon_tools.py`, menu hooks in `app/main.py`, reuse of `app/genetic_code.py`,
  `app/codon_usage.py`, `app/enzymes.py`.
- **restriction-cloning** — `app/assembly.py`, `app/assembly_panel.py`, `app/cut_sites.py`, `app/overhangs.py`,
  `tests/test_assembly.py`. This one also changes the three existing assembly modes for consistency: that part is
  judgment and needs a human read, not just green tests.
- **msa-view** — new `app/msa_view.py`, `app/tools.toml` recipes, view registration (`app/view_picker.py`).
- **library-search** — new tool beside `app/local_align.py`, reusing `app/fileio.py` and `app/multigb.py`.
- **repeats-dotplot** — new `app/repeats.py` plus a drawn view in the style of `app/map_view.py`.

## Uncommitted in the main checkout (`C:/Users/Ant/code/openDNA`)

Two unrelated sets of changes sit in the working tree, on `feature/settings-and-codons`:

**From the alignment and primer audit** (tested: 358 Python tests, 98 Rust tests, `pyflakes` clean):

- `core/src/align.rs` — a mismatch within 11 bases of an alignment end was trimmed away and the hit reported shorter at
  100 % identity; it is now kept when ≥ 5 exact bases lie beyond it (`KEEP_RUN`). Short queries (40–80 bases) with 3–5
  mismatches were reported as a stub or missed entirely; k-mer size now drops to ≤ 9 under 200 bases and a trim that
  would cut more than a quarter of the alignment is refused. Two new Rust tests.
- `app/primer_design.py` — primer design now works at and across the origin of a circular sequence (three-lap search,
  `raw` coordinates for the product size); the dialog's "Across the origin is not supported" refusal is gone.
- `tests/test_primer_origin.py` — one new test covering design at base 1, at the last base, and across the origin.
- `PLAN.md` — the Phase 7 section.

**Not from these sessions** (someone is editing the 3D structure view in this checkout): `app/structure.py`,
`app/structure_view.py`, `tests/test_structure.py` — ligand atoms gain a `group` naming the molecule they belong to.
Left untouched; it must not be swept into a Phase 7 commit.

## How to check a branch

Per [CLAUDE.md](CLAUDE.md) — the full interpreter path matters, a bare `python` is Inkscape's:

```bash
"C:\Users\Ant\anaconda3\envs\opendna\python.exe" -m pytest -q tests        # QT_QPA_PLATFORM=offscreen
"C:\Users\Ant\anaconda3\envs\opendna\python.exe" -m pyflakes app tests
cd core && cargo test                                                       # PYO3_PYTHON must point at that interpreter
```

**The one real hazard:** every worktree shares the single `opendna` conda environment, and `opendna_core` is compiled
into it from whichever checkout last ran `maturin develop`. So a Python test run in one worktree can be testing another
branch's Rust core. Before trusting a Rust-dependent result, rebuild from the branch under test, and do not run two
rebuilds at once. A sudden `ImportError: opendna_core` means a rebuild is in progress elsewhere, not broken code.
