# Unit-test report

What every automated test of openDNA feeds in, what it expects to come out, what the answer is checked against, and
how much a pass is worth. Written so that a user can check, test by test, that the tests check relevant things and check
them well. **[AUDIT.md](AUDIT.md) lists every test on one page with why it matters, its strength and the review notes.**

## How to read a test entry

- **Result**: the outcome of the last run this report was built from (✅ passed, ❌ failed, ⏭️ skipped, · not run).
- **Checked against**: where the expected answer comes from. *hand-computed*: worked out on paper and written in the test; *independent implementation*: another program (mostly Biopython) computes it; *in-test reference model*: a plain model written inside the test (e.g. the sequence as a Python string); *property/invariant*: a rule that must always hold (e.g. undo restores the original); *round-trip/self-consistency*: save then reload gives the same thing; *GUI state*: what the window shows; *smoke*: only that nothing crashes.
- **Strength**: 🟢 strong = exact values or an independent check, edge cases included; 🟡 adequate = the main behaviour, some edges or exact values missing; 🔴 weak = a clearly wrong program could still pass.
- **Why it matters**: what would break unnoticed for a user if the test did not exist. A test whose reason you cannot agree with is a candidate for removal.
- **Input / Expected output**: the concrete cases, in the order the test runs them. Many tests are *topical*: one test function holds several sections, each a row.
- **Review notes**: what the reviewer found could be stronger. **Fixed by the audit**: what was changed.

## Totals (159 tests: 94 Python, 65 Rust)

Last run: ✅ 159 passed · ❌ 0 failed · ⏭️ 0 skipped · · 0 not run

Review: 🟡 adequate: 47 · 🟢 strong: 111 · 🔴 weak: 1

## Test files

| Test file | Tests | Passed | Failed | Skipped | Not run | Strength (reviewed) |
|---|---|---|---|---|---|---|
| [`tests/test_assembly.py`](tests_test_assembly.md) | 8 | 8 | 0 | 0 | 0 | 🟡 adequate 1 · 🟢 strong 7 |
| [`tests/test_checkup.py`](tests_test_checkup.md) | 6 | 6 | 0 | 0 | 0 | 🟡 adequate 2 · 🟢 strong 4 |
| [`tests/test_codons.py`](tests_test_codons.md) | 4 | 4 | 0 | 0 | 0 | 🟢 strong 4 |
| [`tests/test_core_alignment.py`](tests_test_core_alignment.md) | 4 | 4 | 0 | 0 | 0 | 🟡 adequate 2 · 🟢 strong 2 |
| [`tests/test_core_formats.py`](tests_test_core_formats.md) | 4 | 4 | 0 | 0 | 0 | 🟢 strong 4 |
| [`tests/test_core_sequence.py`](tests_test_core_sequence.md) | 6 | 6 | 0 | 0 | 0 | 🟡 adequate 1 · 🟢 strong 5 |
| [`tests/test_documents.py`](tests_test_documents.md) | 8 | 8 | 0 | 0 | 0 | 🟡 adequate 3 · 🟢 strong 5 |
| [`tests/test_enzymes_and_guides.py`](tests_test_enzymes_and_guides.md) | 5 | 5 | 0 | 0 | 0 | 🟢 strong 5 |
| [`tests/test_feedback.py`](tests_test_feedback.md) | 3 | 3 | 0 | 0 | 0 | 🟡 adequate 1 · 🟢 strong 2 |
| [`tests/test_file_assoc.py`](tests_test_file_assoc.md) | 3 | 3 | 0 | 0 | 0 | 🟢 strong 3 |
| [`tests/test_genes_and_history.py`](tests_test_genes_and_history.md) | 2 | 2 | 0 | 0 | 0 | 🟡 adequate 1 · 🟢 strong 1 |
| [`tests/test_genetic_code.py`](tests_test_genetic_code.md) | 2 | 2 | 0 | 0 | 0 | 🟢 strong 2 |
| [`tests/test_input_parsing.py`](tests_test_input_parsing.md) | 1 | 1 | 0 | 0 | 0 | 🟢 strong 1 |
| [`tests/test_multi_record.py`](tests_test_multi_record.md) | 2 | 2 | 0 | 0 | 0 | 🟢 strong 2 |
| [`tests/test_online.py`](tests_test_online.md) | 4 | 4 | 0 | 0 | 0 | 🟡 adequate 2 · 🟢 strong 2 |
| [`tests/test_parts_library.py`](tests_test_parts_library.md) | 2 | 2 | 0 | 0 | 0 | 🟡 adequate 1 · 🟢 strong 1 |
| [`tests/test_primer_binding.py`](tests_test_primer_binding.md) | 5 | 5 | 0 | 0 | 0 | 🟡 adequate 1 · 🟢 strong 4 |
| [`tests/test_primer_design.py`](tests_test_primer_design.md) | 4 | 4 | 0 | 0 | 0 | 🟡 adequate 1 · 🟢 strong 3 |
| [`tests/test_project.py`](tests_test_project.md) | 2 | 2 | 0 | 0 | 0 | 🟡 adequate 2 |
| [`tests/test_protein.py`](tests_test_protein.md) | 1 | 1 | 0 | 0 | 0 | 🟢 strong 1 |
| [`tests/test_reference.py`](tests_test_reference.md) | 2 | 2 | 0 | 0 | 0 | 🟡 adequate 1 · 🟢 strong 1 |
| [`tests/test_repeats_msa.py`](tests_test_repeats_msa.md) | 2 | 2 | 0 | 0 | 0 | 🟡 adequate 2 |
| [`tests/test_roundtrip.py`](tests_test_roundtrip.md) | 6 | 6 | 0 | 0 | 0 | 🟡 adequate 1 · 🟢 strong 5 |
| [`tests/test_search.py`](tests_test_search.md) | 1 | 1 | 0 | 0 | 0 | 🟡 adequate 1 |
| [`tests/test_selftest.py`](tests_test_selftest.md) | 1 | 1 | 0 | 0 | 0 | 🟡 adequate 1 |
| [`tests/test_updater.py`](tests_test_updater.md) | 5 | 5 | 0 | 0 | 0 | 🟡 adequate 1 · 🟢 strong 4 |
| [`tests/test_variants.py`](tests_test_variants.md) | 1 | 1 | 0 | 0 | 0 | 🟢 strong 1 |
| [`core/src/align.rs`](core_src_align.md) | 12 | 12 | 0 | 0 | 0 | 🟡 adequate 3 · 🟢 strong 9 |
| [`core/src/alignments.rs`](core_src_alignments.md) | 1 | 1 | 0 | 0 | 0 | 🟢 strong 1 |
| [`core/src/annotate.rs`](core_src_annotate.md) | 1 | 1 | 0 | 0 | 0 | 🟢 strong 1 |
| [`core/src/bio.rs`](core_src_bio.md) | 3 | 3 | 0 | 0 | 0 | 🟢 strong 3 |
| [`core/src/cas9.rs`](core_src_cas9.md) | 8 | 8 | 0 | 0 | 0 | 🟡 adequate 4 · 🟢 strong 4 |
| [`core/src/code.rs`](core_src_code.md) | 7 | 7 | 0 | 0 | 0 | 🟡 adequate 5 · 🟢 strong 2 |
| [`core/src/dna_file.rs`](core_src_dna_file.md) | 2 | 2 | 0 | 0 | 0 | 🟡 adequate 1 · 🟢 strong 1 |
| [`core/src/document.rs`](core_src_document.md) | 12 | 12 | 0 | 0 | 0 | 🟡 adequate 2 · 🟢 strong 10 |
| [`core/src/fasta.rs`](core_src_fasta.md) | 4 | 4 | 0 | 0 | 0 | 🟡 adequate 1 · 🟢 strong 3 |
| [`core/src/locations.rs`](core_src_locations.md) | 3 | 3 | 0 | 0 | 0 | 🟡 adequate 2 · 🟢 strong 1 |
| [`core/src/primers.rs`](core_src_primers.md) | 5 | 5 | 0 | 0 | 0 | 🟡 adequate 1 · 🟢 strong 3 · 🔴 weak 1 |
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
