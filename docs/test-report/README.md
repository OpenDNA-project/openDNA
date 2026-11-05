# Unit-test report

What every automated test of OpenGenEdit feeds in, what it expects to come out, what the answer is checked against, and
how much a pass is worth. Written so that a user can check, test by test, that the tests check relevant things and check
them well. **The findings of the audit (what was wrong, what was fixed, what is still not tested) are in [AUDIT.md](AUDIT.md).**

## How to read a test entry

- **Result**: the outcome of the last run this report was built from (✅ passed, ❌ failed, ⏭️ skipped, · not run).
- **Checked against**: where the expected answer comes from. *hand-computed*: worked out on paper and written in the test; *independent implementation*: another program (mostly Biopython) computes it; *in-test reference model*: a plain model written inside the test (e.g. the sequence as a Python string); *property/invariant*: a rule that must always hold (e.g. undo restores the original); *round-trip/self-consistency*: save then reload gives the same thing; *GUI state*: what the window shows; *smoke*: only that nothing crashes.
- **Strength**: 🟢 strong = exact values or an independent check, edge cases included; 🟡 adequate = the main behaviour, some edges or exact values missing; 🔴 weak = a clearly wrong program could still pass.
- **Input / Expected output**: the concrete cases, in the order the test runs them. Many tests are *topical*: one test function holds several sections, each a row.
- **Review notes**: what the reviewer found could be stronger. **Fixed by the audit**: what was changed.

## Totals (366 tests: 271 Python, 95 Rust)

Last run: ✅ 364 passed · ❌ 0 failed · ⏭️ 2 skipped · · 0 not run

Review: 🟡 adequate: 102 · not reviewed: 2 · 🟢 strong: 255 · 🔴 weak: 7

## Test files

| Test file | Tests | Passed | Failed | Skipped | Not run | Strength (reviewed) |
|---|---|---|---|---|---|---|
| [`tests/test_align_structural.py`](tests_test_align_structural.md) | 2 | 2 | 0 | 0 | 0 | 🟢 strong 2 |
| [`tests/test_app_logic.py`](tests_test_app_logic.md) | 12 | 12 | 0 | 0 | 0 | not reviewed 1 · 🟢 strong 11 |
| [`tests/test_assembly.py`](tests_test_assembly.md) | 9 | 9 | 0 | 0 | 0 | 🟢 strong 9 |
| [`tests/test_bugfixes.py`](tests_test_bugfixes.md) | 2 | 2 | 0 | 0 | 0 | 🟡 adequate 1 · 🟢 strong 1 |
| [`tests/test_cas9.py`](tests_test_cas9.md) | 4 | 4 | 0 | 0 | 0 | 🟢 strong 4 |
| [`tests/test_checkup.py`](tests_test_checkup.md) | 8 | 8 | 0 | 0 | 0 | 🟡 adequate 4 · 🟢 strong 4 |
| [`tests/test_copy_and_search_strand.py`](tests_test_copy_and_search_strand.md) | 2 | 2 | 0 | 0 | 0 | 🟡 adequate 1 · 🟢 strong 1 |
| [`tests/test_core_functions.py`](tests_test_core_functions.md) | 23 | 23 | 0 | 0 | 0 | 🟡 adequate 4 · 🟢 strong 19 |
| [`tests/test_cut_sites.py`](tests_test_cut_sites.md) | 8 | 8 | 0 | 0 | 0 | 🟡 adequate 2 · 🟢 strong 6 |
| [`tests/test_cut_ui.py`](tests_test_cut_ui.md) | 9 | 9 | 0 | 0 | 0 | 🟡 adequate 5 · not reviewed 1 · 🟢 strong 3 |
| [`tests/test_db_tabs.py`](tests_test_db_tabs.md) | 6 | 6 | 0 | 0 | 0 | 🟡 adequate 1 · 🟢 strong 5 |
| [`tests/test_edited_bases.py`](tests_test_edited_bases.md) | 3 | 3 | 0 | 0 | 0 | 🟡 adequate 2 · 🟢 strong 1 |
| [`tests/test_editing.py`](tests_test_editing.md) | 1 | 1 | 0 | 0 | 0 | 🟢 strong 1 |
| [`tests/test_enzymes.py`](tests_test_enzymes.md) | 2 | 2 | 0 | 0 | 0 | 🟢 strong 2 |
| [`tests/test_feedback.py`](tests_test_feedback.md) | 6 | 6 | 0 | 0 | 0 | 🟢 strong 6 |
| [`tests/test_file_assoc.py`](tests_test_file_assoc.md) | 6 | 6 | 0 | 0 | 0 | 🟢 strong 6 |
| [`tests/test_file_history.py`](tests_test_file_history.md) | 3 | 3 | 0 | 0 | 0 | 🟡 adequate 1 · 🟢 strong 2 |
| [`tests/test_genetic_code.py`](tests_test_genetic_code.md) | 6 | 6 | 0 | 0 | 0 | 🟢 strong 6 |
| [`tests/test_gui.py`](tests_test_gui.md) | 10 | 10 | 0 | 0 | 0 | 🟢 strong 9 · 🔴 weak 1 |
| [`tests/test_history.py`](tests_test_history.md) | 2 | 2 | 0 | 0 | 0 | 🟢 strong 2 |
| [`tests/test_import_export.py`](tests_test_import_export.md) | 7 | 7 | 0 | 0 | 0 | 🟡 adequate 5 · 🟢 strong 2 |
| [`tests/test_input_validation.py`](tests_test_input_validation.md) | 5 | 5 | 0 | 0 | 0 | 🟡 adequate 5 |
| [`tests/test_local_align_pieces.py`](tests_test_local_align_pieces.md) | 3 | 3 | 0 | 0 | 0 | 🟢 strong 3 |
| [`tests/test_map_view.py`](tests_test_map_view.md) | 4 | 3 | 0 | 1 | 0 | 🟡 adequate 3 · 🔴 weak 1 |
| [`tests/test_menus_inputs_and_online.py`](tests_test_menus_inputs_and_online.md) | 9 | 9 | 0 | 0 | 0 | 🟡 adequate 2 · 🟢 strong 7 |
| [`tests/test_multi_genbank.py`](tests_test_multi_genbank.md) | 3 | 2 | 0 | 1 | 0 | 🟡 adequate 1 · 🟢 strong 1 · 🔴 weak 1 |
| [`tests/test_ncbi.py`](tests_test_ncbi.md) | 2 | 2 | 0 | 0 | 0 | 🟡 adequate 1 · 🟢 strong 1 |
| [`tests/test_parts_library.py`](tests_test_parts_library.md) | 3 | 3 | 0 | 0 | 0 | 🟡 adequate 1 · 🟢 strong 2 |
| [`tests/test_primer_binding.py`](tests_test_primer_binding.md) | 17 | 17 | 0 | 0 | 0 | 🟡 adequate 1 · 🟢 strong 16 |
| [`tests/test_primer_origin.py`](tests_test_primer_origin.md) | 3 | 3 | 0 | 0 | 0 | 🟡 adequate 1 · 🟢 strong 2 |
| [`tests/test_primers.py`](tests_test_primers.md) | 3 | 3 | 0 | 0 | 0 | 🟢 strong 3 |
| [`tests/test_project.py`](tests_test_project.md) | 2 | 2 | 0 | 0 | 0 | 🟡 adequate 1 · 🟢 strong 1 |
| [`tests/test_protein.py`](tests_test_protein.md) | 3 | 3 | 0 | 0 | 0 | 🟡 adequate 3 |
| [`tests/test_real_files.py`](tests_test_real_files.md) | 1 | 1 | 0 | 0 | 0 | 🟢 strong 1 |
| [`tests/test_reference.py`](tests_test_reference.md) | 4 | 4 | 0 | 0 | 0 | 🟢 strong 4 |
| [`tests/test_roundtrip.py`](tests_test_roundtrip.md) | 14 | 14 | 0 | 0 | 0 | 🟢 strong 14 |
| [`tests/test_search.py`](tests_test_search.md) | 3 | 3 | 0 | 0 | 0 | 🟡 adequate 1 · 🟢 strong 2 |
| [`tests/test_selftest.py`](tests_test_selftest.md) | 1 | 1 | 0 | 0 | 0 | 🟡 adequate 1 |
| [`tests/test_structure.py`](tests_test_structure.md) | 3 | 3 | 0 | 0 | 0 | 🟡 adequate 1 · 🟢 strong 2 |
| [`tests/test_ui_stability.py`](tests_test_ui_stability.md) | 9 | 9 | 0 | 0 | 0 | 🟡 adequate 4 · 🟢 strong 5 |
| [`tests/test_ui_tools.py`](tests_test_ui_tools.md) | 11 | 11 | 0 | 0 | 0 | 🟡 adequate 6 · 🟢 strong 5 |
| [`tests/test_ui_views.py`](tests_test_ui_views.md) | 12 | 12 | 0 | 0 | 0 | 🟡 adequate 7 · 🟢 strong 5 |
| [`tests/test_updater.py`](tests_test_updater.md) | 14 | 14 | 0 | 0 | 0 | 🟡 adequate 3 · 🟢 strong 11 |
| [`tests/test_variants.py`](tests_test_variants.md) | 6 | 6 | 0 | 0 | 0 | 🟢 strong 6 |
| [`tests/test_view_edges.py`](tests_test_view_edges.md) | 3 | 3 | 0 | 0 | 0 | 🟡 adequate 2 · 🟢 strong 1 |
| [`tests/test_view_help.py`](tests_test_view_help.md) | 2 | 2 | 0 | 0 | 0 | 🟡 adequate 2 |
| [`core/src/align.rs`](core_src_align.md) | 11 | 11 | 0 | 0 | 0 | 🟡 adequate 3 · 🟢 strong 7 · 🔴 weak 1 |
| [`core/src/alignments.rs`](core_src_alignments.md) | 1 | 1 | 0 | 0 | 0 | 🟢 strong 1 |
| [`core/src/annotate.rs`](core_src_annotate.md) | 1 | 1 | 0 | 0 | 0 | 🟢 strong 1 |
| [`core/src/bio.rs`](core_src_bio.md) | 3 | 3 | 0 | 0 | 0 | 🟢 strong 3 |
| [`core/src/cas9.rs`](core_src_cas9.md) | 8 | 8 | 0 | 0 | 0 | 🟡 adequate 4 · 🟢 strong 4 |
| [`core/src/code.rs`](core_src_code.md) | 7 | 7 | 0 | 0 | 0 | 🟡 adequate 4 · 🟢 strong 2 · 🔴 weak 1 |
| [`core/src/dna_file.rs`](core_src_dna_file.md) | 2 | 2 | 0 | 0 | 0 | 🟡 adequate 1 · 🟢 strong 1 |
| [`core/src/document.rs`](core_src_document.md) | 43 | 43 | 0 | 0 | 0 | 🟡 adequate 12 · 🟢 strong 31 |
| [`core/src/fasta.rs`](core_src_fasta.md) | 4 | 4 | 0 | 0 | 0 | 🟡 adequate 1 · 🟢 strong 3 |
| [`core/src/locations.rs`](core_src_locations.md) | 3 | 3 | 0 | 0 | 0 | 🟡 adequate 2 · 🟢 strong 1 |
| [`core/src/primers.rs`](core_src_primers.md) | 5 | 5 | 0 | 0 | 0 | 🟢 strong 3 · 🔴 weak 2 |
| [`core/src/search.rs`](core_src_search.md) | 6 | 6 | 0 | 0 | 0 | 🟡 adequate 2 · 🟢 strong 4 |
| [`core/src/tracks.rs`](core_src_tracks.md) | 1 | 1 | 0 | 0 | 0 | 🟡 adequate 1 |

## Rebuilding this report

The descriptions are in [tests.json](tests.json); the pages are written by `scripts/test_report.py`, which also lists any
test that has no description yet (with its docstring and asserts taken from the code), so a new test never goes missing:

```bash
QT_QPA_PLATFORM=offscreen python -m pytest tests -q --junitxml=junit.xml
(cd core && cargo test > ../cargo.txt 2>&1)
python scripts/test_report.py --junit junit.xml --cargo cargo.txt
```
