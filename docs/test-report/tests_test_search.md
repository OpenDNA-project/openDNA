# `tests/test_search.py`

[← all test files](README.md) · source: [`tests/test_search.py`](../../tests/test_search.py)

## `test_sequence_search_engine`

**Result:** ✅ passed in 0.01 s · [line 46](../../tests/test_search.py#L46) · **Area:** Sequence search: DNA / IUPAC / mismatches / protein · **Checked against:** hand-computed · **Strength:** 🟢 strong

Search finds DNA on both strands with IUPAC codes and a bounded number of mismatches, wraps across a circular origin, and finds protein motifs in translations on both strands.

| # | Input | Expected output |
|---|---|---|
| 1 | 2000 bp random seq (seed 5) with GAATTC at 300, ACGTACGTAC at 700, rc('GGATCCAAAT') at 1200; query 'GGATCCAAAT' | [(1200, 1210, -1)]; with both=False: [] |
| 2 | 'GRATTC' with IUPAC on/off; 'ACGTNCGTAC' | (300,306) among hits / []; (700,710) among hits |
| 3 | SEQ[700:710] with base 4 mutated; fuzzy 0 and 1 | (700,710) not found; found with mismatches == 1 |
| 4 | allowed_mismatches(fuzzy=3) for length 4 and 30 | 1 and 3 |
| 5 | circular 'TTCAAAAAAAAAAAAAAGA' (19 bp), 'GATTC' forward only | [(17, 22)]; search_spans split [(0,3),(17,19)]; hits_between(0,5)==[(0,3)], (10,16)==[] |
| 6 | 'ATGGCTAAATGGCATTAAAAAA' (M A K W H *), 'AKW' translations on/off | [(3, 12, 'protein', 1)]; [] with translations off |
| 7 | reverse complement of that seq, 'AXW' (X = any) | 1 hit on strand -1 |

**Review notes:**

- Line 82: the reverse-strand protein hit checks only count and strand, not its coordinates (expected 10..19).

**Fixed by the audit of November 2026:**

- The two misleading comments (lines 66 and 69) were corrected.

## `test_names_scope_and_file_search`

**Result:** ✅ passed in 0.36 s · [line 85](../../tests/test_search.py#L85) · **Area:** Sequence search: names, primers, enzymes, scope, file · **Checked against:** hand-computed · **Strength:** 🟢 strong

Search finds features by name (notes only with 'details'), primers by name where they bind, enzymes by name with that option, respects the selection scope, and searches every sequence of a multi-FASTA file.

| # | Input | Expected output |
|---|---|---|
| 1 | CDS 101..400 gene=abcD product=kinase db_xref=needle; queries 'abcd', 'needle' (details off/on) | [(100,400)]; []; [(100,400)] |
| 2 | primer 'fwdprimer' = SEQ[500:520]; query 'fwdprim' primers on/off | [(500,520)]; [] |
| 3 | 'EcoRI' with enzymes=True / default | [(300, 306, 'enzyme')]; [] |
| 4 | scope=selection, 'GAATTC', selection 250-350 / 400-500 / bare cursor at 10 | [(300,306)] no note; []; [(300,306)] with note 'Nothing is selected' |
| 5 | g.fa chrA ('ACGT'*50+'GGTACC'+'TTTT'*20), chrB ('TTTT'*30+'GGTACC'+'CCATGG'+'ACAC'*10); search_file 'GGTACC' | returns 'GGTACC'; chrA [(200,206)], chrB [(120,126)] (palindrome: one hit per site) |
| 6 | search_file for 'kinase' | None (not a sequence motif) |

## `test_search_bar`

**Result:** ✅ passed in 0.08 s · [line 146](../../tests/test_search.py#L146) · **Area:** GUI: search bar · **Checked against:** hand-computed · **Strength:** 🟡 adequate

The search bar counts hits, steps next/previous with wrap-around (also via Edit menu/F3), clears stale hits after an edit, and applies and remembers its options.

| # | Input | Expected output |
|---|---|---|
| 1 | seq with GAATTC at 100, 206, 312; type 'GAATTC' | search_spans [(100,106),(206,212),(312,318)], counter '3 hits' |
| 2 | cursor 0, step(+1) x4, step(-1), find_next, find_previous | selects 100 ('1 of 3'), 312 ('3 of 3'), wraps to 100, back to 312, 100, 312 |
| 3 | cursor 210 then step(+1); cursor 206 then step(+1) | 312; 206 (hit at the cursor counts) |
| 4 | type a base after searching; then 'ATATATAT'; empty text; escape() | spans cleared, counter ''; 'No hits'; '' and []; [] |
| 5 | 'GGATCC'; fuzzy on with 2; both off; scope selection; new SearchBar | 'No hits'; options().fuzzy==2 and settings fuzzy_n==2, both True; both False saved; scope 'selection'; new bar options equal |
| 6 | fixture: sys.excepthook collected during the test | no uncaught exceptions (e.g. while painting hits) |

**Fixed by the audit of November 2026:**

- After fuzzy = 2 the test now asserts the search finds hits (GGATCC is 2 bases from the GAATTC sites), so the option is shown to be applied, not only stored.
