# `tests/test_parts_library.py`

[← all test files](README.md) · source: [`tests/test_parts_library.py`](../../tests/test_parts_library.py)

## `test_collections_and_switches`

**Result:** ✅ passed in 0.02 s · [line 64](../../tests/test_parts_library.py#L64) · **Area:** Features Database: collections · **Checked against:** hand-computed · **Strength:** 🟢 strong

Each collection file carries its title, reference, DOI and measured strengths; a source switched off is not used; the user's own parts form the 'Local' collection and skip duplicate sequences; broken rows are skipped.

| # | Input | Expected output |
|---|---|---|
| 1 | measured.tsv (2 terminators, TE 88 and 3.5) and vectors.tsv (AmpR CDS) | ids {measured, vectors}; title/doi/reference parsed; counts {terminator: 2}, measured 2; vectors {CDS: 1}, measured 0 |
| 2 | load_parts() | Tstrong strength 88, units TE, organism E. coli, collection measured; describe 'TE 88 (E. coli)'; AmpR describe '' |
| 3 | set_on('measured', False) then True | only ['AmpR']; is_on false; collections [(measured,False),(vectors,True)]; back to 3 parts |
| 4 | add_parts mine (lowercase), again (same seq upper), 'AmpR copy' (same seq as AmpR) | 1 added; Local collection {promoter: 1}; parts {AmpR, mine}; Local off -> {AmpR} |
| 5 | odd.tsv: empty sequence, protein 'MKPG', strength 'very strong' | only 'bad number' kept, strength None, describe '' |

## `test_the_bundled_collections_are_readable_and_the_database_lists_them`

**Result:** ✅ passed in 0.14 s · [line 112](../../tests/test_parts_library.py#L112) · **Area:** Features Database: shipped files · **Checked against:** round-trip/self-consistency · **Strength:** 🟡 adequate

Every collection file shipped in app/parts parses with a title and at least one part, and the Features Database tab shows every part with one Source value per collection.

| # | Input | Expected output |
|---|---|---|
| 1 | BUNDLED_DIR = app/parts; collections() | non-empty, each with a title and > 0 parts |
| 2 | DatabaseTab(FEATURES) | row count == total of collection counts; Source filter values == collection titles |

**Review notes:**

- Line 120: the expected row count comes from pl.collections() counts, the same parser that feeds the table, so a parsing bug that drops rows from both would not be seen; no known part or count is pinned.

## `test_local_features_are_a_source_and_can_be_removed`

**Result:** ✅ passed in 0.01 s · [line 125](../../tests/test_parts_library.py#L125) · **Area:** Features Database: user's own parts · **Checked against:** hand-computed · **Strength:** 🟢 strong

Parts the user adds belong to the source 'Local' and remember their file; only Local parts can be removed, and the user file keeps its header and remaining rows.

| # | Input | Expected output |
|---|---|---|
| 1 | vectors.tsv + add_parts mine (promoter, my.gb), other (RBS, b.gb) | mine: (LOCAL, 'my.gb'); AmpR source_name 'Cloning vectors (NCBI records)' |
| 2 | remove_parts([mine, AmpR]) | returns 1; remaining {AmpR, other}; user file starts '#title: Local' header |
