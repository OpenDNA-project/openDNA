# Audit of the unit tests (November 2026)

[← unit-test report](README.md)

**The question:** do the tests run, do they check what they claim to check, and what is still not tested at all?

**The short answer:** the suite is in good shape. Every Rust test passes. Every Python test passes apart from one that depends on the font and failed on Linux (now fixed). The scientific core (Tm, ORFs, translation, restriction cuts, primer binding, GenBank round trips, pileups) is mostly checked against independent references: Biopython, hand-worked values, or a plain model written inside the test. Not against numbers copied from the program. The audit found **one test that checked nothing**, a handful of checks that could not fail, and some untested code. That code includes `app/variants.py`, which writes read variants into the user's DNA and was 13 % covered. These are fixed or listed below.

## 1. How the audit was done

| Step | What |
|---|---|
| Ran everything | `cargo test` (Rust) and `pytest tests` (Python, Qt offscreen) on Linux, Python 3.13, the Rust core built with maturin, as in CI. `pyflakes app tests scripts` is clean. |
| Static scan | Every Python test parsed: none has no assertion, none catches exceptions broadly (`except:` / `except Exception`), none is marked `xfail` or `skip`; 8 call `pytest.skip()` at run time (listed in §4). |
| History check | Commit `a758ace` merged 866 small tests into 218 topical ones. The number of `assert` statements went *up* (2743 → 2785), so no check was lost in the merge. |
| Test-by-test review | Each of the 375 tests was read alongside the code it tests. Its inputs and expected outputs were written down ([the report](README.md)), and each test was rated strong, adequate or weak. Coordinates, reverse complements, enzyme cuts, codon tables and Tm values in the expected answers were spot-checked by hand. **No biological or arithmetic error was found in any expected value.** |
| Claims checked | Every problem below marked *fixed* was reproduced first (e.g. the empty primer list, the overwritten origin site, the PAM pixel count) and the fixed test was run. |
| Coverage | Python line coverage with `pytest-cov`. Rust coverage with LLVM source-based coverage (`-C instrument-coverage`), from `cargo test` and from the Python suite running the instrumented core. |

Two tests have no written entry, only the code extract the report generator makes: `test_app_logic.py::test_history_texts` and `test_cut_ui.py::test_spacer_finder_automatic_check`. Their pages show their docstrings and every `assert` line as written.

## 2. Results of the run

| Suite | Tests | Passed | Failed | Skipped |
|---|---|---|---|---|
| Rust (`cargo test`) | 101 | 101 | 0 | 0 |
| Python, before the audit | 267 | 263 | 1 | 3 |
| Python, after the fixes and new tests (bwa, samtools and minimap2 installed) | 274 | 272 | 0 | 2 |
| Python, after removing duplicates (§4b) | 271 | 269 | 0 | 2 |
| Rust, after removing duplicates (§4b) | 95 | 95 | 0 | 0 |

- **The failure:** `test_cut_sites.py::test_a_saved_guide_site_is_drawn`. The drawing is right (the PAM box is about 27 px wide). But the test counted PAM-coloured pixels along one pixel row, and that row runs through the black "PAM" label inside the box. With the font used on Linux (DejaVu Sans), the letters cover most of the row and only 11 pixels were left, under the 18 needed. Windows fonts are narrower, so CI passed. *Fixed:* box widths are now measured column by column over the whole lane height.
- **The skips** (they skip in CI too, so these parts are never checked automatically):
  - the read-mapping half of `test_gui.py::test_read_colours_local_align_and_map_reads`: mapping, variant calling, coverage, applying variants and cancelling. It needs bwa, samtools and minimap2, which the Windows CI does not have. Installed for this audit, it **passes**;
  - `test_map_view.py::test_human_genome_chromosomes_listed_drawn_and_zoomed` and `test_multi_genbank.py::test_vibrio_cholerae_two_chromosomes`: they need a human genome and a *V. cholerae* file at `V:/HM_Softwares/TestData/...`, a drive outside the repository.

## 3. Problems found and fixed

| Test | Problem | Fix |
|---|---|---|
| `test_primers.py::test_primers_bind_where_the_dna_file_says` | **Checked nothing.** It imported `sample_files.PLASMID` (a GenBank file) under the name `PLASMID_DNA`. That file's primer list is empty, so the loop never ran and the test always passed. Pointed at the real `.dna` file, its exact-span comparison would also have failed for oAC812 / oAC813, whose stored site includes a 5′ tail. | Replaced by `test_every_primer_the_dna_file_annotates_is_found_on_its_strand_with_the_same_3_prime_end`. It checks all 15 stored primers, asserts the list is not empty, and names the two primers stored with a tail. This is the direction `test_primer_binding.py` did not cover (it checks that every site we find is stored, not that every stored site is found). |
| `test_primer_binding.py::test_our_sites_are_the_binding_sites_stored_in_each_real_dna_file` | Same risk in the other direction: its loops assert nothing if they find nothing. | Counts the sites it compares: at least one per stored primer (15). |
| `test_core_functions.py::test_digest_with_every_enzyme_planted_in_a_random_circle_equals_biopython` | The EcoRI site planted across the origin was overwritten by the BsaI site planted across the same origin. So the "EcoRI across the origin" case never existed in any of the 5 circles. The final check passed because of an EcoRI site elsewhere. | EcoRI and BsaI take turns, one per circle. The cut across the origin is asserted at its exact position. |
| `test_cut_sites.py::test_a_saved_guide_site_is_drawn` | Font-dependent failure (above). | Measure over the lane height. |
| `test_map_view.py::test_map_layout_crowding_clicks_scale_and_geometry` | `assert hit is None or hit in kept` also passes when `_hit()` never finds anything. | The test now counts the hits too: more than 40 of the 320 clicks must land on a feature (147 do). |
| `test_search.py::test_search_bar` | The docstring says the options are applied, but only the stored values were read back. | After *fuzzy = 2* the test asserts that the search now finds the GAATTC sites, which are 2 mismatches from GGATCC. |
| `test_search.py::test_sequence_search_engine` | One comment said the opposite of its assertion ("1 mismatch in 4 bases is refused"; it is allowed). Another comment was garbled. | Comments corrected. |
| `test_updater.py::test_job_gc_handling` | "Collector back on only after the last of several jobs" was not shown, because both jobs were released together. | Separate gates: after the first job ends the collector must still be off. |
| `test_updater.py::test_release_json` | `error is None or "not understood" in error` accepted either outcome for every tag. | Only the `v0.3.0-beta1` tag may give the message. |
| `test_enzymes.py::test_the_digest_dialog_filter_gel_ladders_and_lanes` | The docstring came after an `import`, so the test had no `__doc__`. | Moved. |
| `core/src/tracks.rs::gff_bed_and_queries` | `assert_eq!(t.query(..,2), (t.query(..,2).0, true))` compared the result with itself. Only the "cut" flag was really checked. | Asserts the two items kept (lacZ, big,gene) and that no cut is flagged when 3 items fit a limit of 3. |

## 4. New tests

| Test | Why |
|---|---|
| `tests/test_variants.py` (6 tests) | `app/variants.py` writes variants found in reads into the sequence. It was 13 % covered, and its only test path (the mapping section above) never runs in CI. The new tests compare every kind of change (SNP, deletion, insertion, at the first and last base and after the end) with the same edit on a plain string. They also check: several changes applied right to left in either order, one undo step each, the DNA lock, unusable insertions skipped, and the `variation` feature locations (1-based), `/replace` and read counts. They need no external tools. |
| `test_primers.py::test_tm_of_self_complementary_primers_matches_biopython` | The Tm formula has a separate branch for primers that are their own reverse complement (whole concentration instead of C/4, plus a symmetry term). The Biopython comparison skips palindromes and the Rust test uses none, so the branch was untested everywhere. Six palindromes at three buffer conditions now agree with Biopython's `selfcomp` mode within 0.01 °C. The test also asserts that the other formula would be more than 1 °C off. |
| `test_primers.py::test_every_primer_the_dna_file_annotates_...` | See §3. |

## 4b. Fewer tests to review

To make the suite quicker to check by hand, every pair of tests that looked alike was compared line by line. Most were not duplicates: they use different inputs or references, or each holds a check found nowhere else. Those were kept. The tests removed below were each fully covered by another test, and nothing biological was lost. **Result: 375 → 366 tests, about 200 fewer lines to read.**

| Removed | Its checks are now in |
|---|---|
| `core/src/annotate.rs::translates_with_the_standard_code` | `code.rs::the_standard_code_and_the_known_differences` (the identical assertion on the same function) |
| `core/src/document.rs::random_edits_then_full_undo_is_lossless` | `random_edits_undo_redo_follow_a_snapshot_and_string_model` (600 edits per topology, every step compared) |
| `core/src/document.rs::redo_is_dropped_by_a_new_edit_and_replays_in_order` | the same snapshot test (redo gone after every edit, every redo replayed) |
| `core/src/document.rs::deleting_a_whole_feature_removes_it` | `deleting_exactly_the_feature_or_more_removes_it_and_undo_restores_it_whole` (stricter) |
| `core/src/code.rs::alternative_start_codons_only_count_when_asked` | `bio.rs::alternative_starts_and_other_tables`, and `test_genetic_code.py` against Biopython |
| `core/src/locations.rs::insert_rules` | `document.rs::insert_at_the_exact_start_or_end_of_a_feature_does_not_grow_it` and `edit_moves_features_and_undo_restores_exactly`, with the bases compared |
| `test_app_logic.py::test_autosave_and_recovery` | merged into `test_roundtrip.py::test_autosave_and_recovery`, with its own checks added (60 s timer, no rewrite without a change, a feature added by the user recovered, the `*` title) |
| `test_gui.py::test_about_dialog_and_light_dark_themes` | merged into `test_menus_inputs_and_online.py::test_about_dialog_and_view_theme` |
| `test_real_files.py::test_genome_opens_and_round_trips` | merged into `test_roundtrip.py::test_real_file_round_trips`. The genome joins its byte-identical save loop, with its size, topology and feature count checked. |

Two repeated sections inside tests were also removed. The fuzzy-location block of `test_bugfixes.py`, which `test_core_functions.py::test_locations_in_genbank_files` checks more strongly. And the re-run of every self-test through the panel in `test_selftest.py`; the panel itself is tested with fake checks in `test_app_logic.py`.

**Not removed, by choice: Rust tests repeated by Python tests through the bindings.** About 17 Rust tests in `align.rs`, `primers.rs`, `annotate.rs`, `code.rs` and `bio.rs` are also checked, more strongly, by Python tests calling the same Rust functions. CI runs both suites, so removing them would lose nothing there. But a developer who runs only `cargo test` would lose every check of `primers.rs` and `annotate.rs`, and most of `align.rs`. That is a decision for the maintainers. Two assertions would have to move to Python first: `tm("A") is None`, and the full-duplex ΔG < −20 kcal/mol.

## 5. Problems found and not fixed

These are listed per test in the report under **Review notes** (search the pages for them). The ones that matter most:

**Checks weaker than their names** (a wrong program could pass):
- `core/src/align.rs::max_hits_still_limits_the_number_of_pieces_and_keeps_the_best` checks the counts only, not that the best pieces are kept.
- `core/src/document.rs::deleting_exactly_the_feature_or_more_removes_it_and_undo_restores_it_whole` has a "keep what is left" branch that none of its 5 ranges reaches. Partial deletion is covered by the test just before it.
- `core/src/document.rs::edited_ranges_basics_origin_and_save` checks the width of the edited range across the origin, not where it is.
- `gc_bins_exact_and_from_blocks` tests the cached-block path only on all-GC and all-AT sequences.
- `core/src/primers.rs::tm_is_sane`, `dimers_and_hairpins` and `code.rs::every_table_has_64_codons_and_known_ids` use loose windows or lengths only. Exact Tm, ΔG and codon tables are checked in the Python tests (`test_primers.py`, `test_core_functions.py`, `test_genetic_code.py`), so these are weak only as stand-alone Rust tests.
- `test_gui.py::test_genome_is_responsive` prints timings but asserts no time limit.
- `test_gui.py::test_display_filters_translation_menus_history_and_new_file` (L1103-1104): "Find also looks for the bases" also passes when the hit is the feature found by name.
- `test_cut_ui.py::test_guide_look_and_interface_highlight` (L145, L147): two of its "points up" comparisons are always true; the "hollow" outline is not checked.
- `test_core_functions.py::test_genome_fasta_fetch_and_motif_search` compares `Fasta.search_motif` with `find_motif`, which it is built on.
- `test_core_functions.py::test_orfs_agree_with_a_brute_force_scan_and_wrap_the_origin` checks the ORF across the origin by count and strand only. The Rust test `bio.rs::finds_an_orf_that_crosses_the_origin_of_a_ring` pins its coordinates.
- `test_protein.py::test_mass_pi_charge_and_layout` checks pI and charge with inequalities only (basic peptide pI > 10, acidic < 4), so a wrong pKa table could pass. In `test_follows_the_selection_...` the after-an-edit section accepts "nothing shown".
- `test_copy_and_search_strand.py::test_copy_and_copy_reverse_complement` claims each copied feature keeps its bases but checks only its strand. Its across-the-origin selection carries no feature at all.
- `test_checkup.py::test_origins_and_prediction_thresholds` checks the promoter/terminator thresholds on one random sequence with a loose bound (≤ 8), counting with the functions under test.
- `test_cut_sites.py::test_enzyme_tooltips_by_site_position` checks, for every enzyme, only that the tooltip contains its name and the labels. The exact coordinates are checked for EcoRI, KpnI, EcoRV and BsaI in other tests.
- `test_cut_ui.py::test_right_click_adds_the_same_feature_as_the_button` (L285): "add 2 rows" re-adds the guide of row 0 and checks only that the count rises by 2.
- `test_primer_origin.py` (L153): "the temporary site is not provisional any more" also passes if the site is gone.

**Robustness (could break or flake, but checks the right thing):**
- Shared state: several tests change the module-wide primer database, settings or `theme.is_dark` without `try/finally` (`test_gui.py`, `test_ui_tools.py`, `test_db_tabs.py`, `test_menus_inputs_and_online.py`). One failure can then cascade into later tests. `test_history.py` shows the isolated pattern to copy.
- Order dependence: `test_enzymes.py` L97 checks EcoRV only if EcoRV cuts, which depends on the seeded random generator having been used by the previous test.
- Timing: background-thread results are polled with 5–10 s wall-clock loops (`test_ncbi.py`, `test_structure.py`, `test_reference.py`, `test_import_export.py`, `test_multi_genbank.py`, `test_updater.py`, `test_db_tabs.py`). Some end without saying what they waited for. `test_ncbi.py` (L93-100) asserts a real elapsed time under 0.6 s.
- Rust tests write fixed file names (`ge_test.gff3`, `ge_test_*.gb`, ...) into the shared temp directory and mostly leave them there, so two runs at once could collide.
- Silent partial skips: some sections run only `if <demo file>.exists()`, with no skip reported (`test_roundtrip.py`, `test_import_export.py`). The demo files are committed, so these sections do run today.

## 6. What is not tested

**Python line coverage: 89.0 %** (16 423 of 18 453 statements in `app/`, before the new tests). The least covered modules, and the untested functions that matter:

| Module | Covered | Not tested |
|---|---|---|
| `variants.py` | 13 % → covered by the new tests | (was: applying variants, the `variation` features) |
| `runner.py` | 51 % | the mapping-tools dialog and running bwa/minimap2 (`__init__`, `_accept`, `run_folder`): needs the tools |
| `alignments.py` | 53 % | the background variant scan (`start_scan`, `work`), `flag_level`, `describe` |
| `genome_view.py` | 65 % | the Alignment view's BAM/variant tracks and their right-click menus (the GUI path to `variants.py`) |
| `main.py` | 78 % | menus and dialogs not opened by any test |
| `file_assoc.py` | 80 % | the Windows-registry functions (`_open`, `get`, `exists`, `delete_value`): not reachable on Linux, so check the Windows CI before calling them untested |
| `feature_table.py` | 81 % | the right-click menu (`_menu`), group as gene, delete from the table |
| `enzymes.py` | 84 % | `Enzyme.summary` (the enzyme's one-line description) |
| `ncbi.py` | 85 % | `download_genome`, `open_url` (network; the rest uses fakes) |
| `structure.py` | 87 % | `fetch` (the real HTTP call; every source above it is tested with fakes) |
| `updater.py` | 89 % | the real `fetch`: its cancel, its size limits (`MAX_TEXT`, `MAX_EXE`) and its early stop. The tests use their own fake. |

**Rust coverage from `cargo test` alone:** 84 % of lines and 71 % of functions in `core/src`. These figures include the test code that sits in the same files. `lib.rs` (the Python bindings, 0 % here) is reached only through the Python tests. Measuring it from the Python suite running an instrumented core was attempted, but it failed: the process aborts at exit, in the same Qt teardown crash that `tests/conftest.py` avoids with `os._exit`, before the coverage data is written. So the combined figure is not given. Least covered: `locations.rs` 83 % of lines (no Rust test edits a sequence carrying an `a^b`, `order(...)`, `bond(...)`, `one-of(...)` or remote location), `search.rs` 88 %, `code.rs` 89 %, `dna_file.rs` 90 %, `fasta.rs` 90 %.

**Behaviours with no test (and how to close them):**
1. **Read mapping end to end in CI.** For this audit, bwa, samtools and minimap2 were installed from the Ubuntu packages and the skipped section of `test_read_colours_local_align_and_map_reads` was run: it **passes**. It just never runs in CI. Installing the three packages in a Linux CI job (or adding a small BAM they produced as test data) would keep it checked.
2. **Multi-chromosome genomes.** Add a small two-record GenBank and a small multi-sequence FASTA to `demo/`, so `test_vibrio_cholerae_two_chromosomes` and `test_human_genome_...` have a stand-in that is in the repository. Today they only run on the author's machine.
3. **The real download path of the updater** (`updater.fetch`): its size limits and cancel, with a local HTTP server.
4. **A Linux job in CI.** CI runs only on Windows. The font failure above shows that the pixel tests depend on the platform. A Linux job would catch that (and it is how this audit ran).

## 7. Documentation corrected

- `docs/TEST_CATALOGUE.md` said "810 tests" and "100 tests" (true before the merge commit) and "Nothing is skipped". It now gives the current counts, names the three tests that skip, and links this report.
- `docs/TESTING.md` links this report.
