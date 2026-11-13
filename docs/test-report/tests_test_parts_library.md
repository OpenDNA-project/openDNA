# `tests/test_parts_library.py`

[← all test files](README.md) · source: [`tests/test_parts_library.py`](../../tests/test_parts_library.py)

## `test_collections_sources_and_switches`

**Result:** ✅ passed in 0.13 s · [line 60](../../tests/test_parts_library.py#L60) · **Area:** Features Database (parts library) · **Checked against:** hand-computed · **Strength:** 🟢 strong

**Why it matters:** Without it, the parts library could show parts from a source the user switched off, lose source/measurement metadata, let users delete shipped parts, or crash on one bad row.

Using temporary collection files, checks collection metadata and measured-strength parsing, describe() text, switching a source off and on, the user's own 'Local' parts (dedup by sequence ignoring case and against shipped parts, its own switch), that only Local parts can be removed with the header kept, and that broken rows are skipped.

| # | Input | Expected output |
|---|---|---|
| 1 | measured.tsv (2 terminators with TE 88 and 3.5, title/reference/doi) and vectors.tsv (AmpR CDS) | collections {measured, vectors}; title, doi, reference parsed; counts {terminator:2}, measured 2; vectors {CDS:1}, measured 0 |
| 2 | load_parts / describe | Tstrong strength 88, units TE, organism E. coli, collection measured; describe 'TE 88 (E. coli)'; AmpR describe '' |
| 3 | set_on('measured', False) then True | only AmpR loaded; is_on false; collections switches [(measured, False), (vectors, True)]; back to 3 parts |
| 4 | add_parts: 'mine', 'again' (same sequence other case), 'AmpR copy' (same as shipped AmpR) | returns 1; Local collection title LOCAL, counts {promoter:1}; parts {AmpR, mine}; Local switched off -> {AmpR} |
| 5 | add mine (my.gb) and other (b.gb); remove mine and AmpR | mine has source_name LOCAL and source 'my.gb'; AmpR source_name 'Cloning vectors (NCBI records)'; remove returns 1; remaining {AmpR, other}; user file still starts '#title: <LOCAL>' |
| 6 | odd.tsv: row with no sequence, protein row 'MKPG', terminator with strength 'very strong' | only 'bad number' loaded; strength None; describe '' |

**Review notes:**

- The first two sections rewrite the same two files without clearing; harmless but redundant.
- Five subtopics in one test (the docstring calls it a merged test).

## `test_shipped_parts_meet_the_data_quality_rules`

**Result:** ✅ passed in 0.02 s · [line 135](../../tests/test_parts_library.py#L135) · **Area:** Features Database (shipped data) · **Checked against:** property/invariant · **Strength:** 🟡 adequate

**Why it matters:** Without it, the shipped parts library could gain junk entries (primer sites, cloning-site slices, unnamed types), parts too short to annotate meaningfully, duplicates, or parts without a citation, cluttering every user's automatic annotation.

Loads the real app/parts collections and checks each file parses to a titled non-empty source, no part name matches build_parts.NEVER, no vague types, every part meets the length floor for its type, no repeated names or sequences (raw rows included), and every part cites a source record.

| # | Input | Expected output |
|---|---|---|
| 1 | pl.collections() over app/parts | non-empty; every collection titled with at least one part; total loaded parts equals the sum of counts |
| 2 | Each shipped part | name not matched by build_parts.NEVER; kind not misc_feature/regulatory/sig_peptide; build_parts.FLOORS has a floor for its kind and length >= floor; source not blank |
| 3 | All part names | no duplicates |
| 4 | Raw rows of every file, sequence upper-cased | no sequence stored more than once across or within collections |

**Review notes:**

- A data-lint check on shipped files rather than a code test; it fails when data changes, not when code breaks. Still useful as a guard.
- The rules (NEVER, FLOORS) come from scripts/build_parts.py, so a loosened rule there would pass silently.
- Duplicate check is exact sequence only; reverse-complement duplicates or near-duplicates are not caught.
