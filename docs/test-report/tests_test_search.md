# `tests/test_search.py`

[← all test files](README.md) · source: [`tests/test_search.py`](../../tests/test_search.py)

## `test_sequence_search_engine`

**Result:** ✅ passed in 0.02 s · [line 44](../../tests/test_search.py#L44) · **Area:** Sequence search (Find) · **Checked against:** hand-computed · **Strength:** 🟡 adequate

**Why it matters:** Without it, the Find box could miss DNA matches on the other strand or with IUPAC/mismatches, fail across the origin, ignore its options for names, primers, enzymes and translations, or search outside the selection.

Runs search.search_sequence on a 2000 bp random sequence with planted sites and checks DNA, IUPAC, fuzzy, circular, protein, name/detail/primer/enzyme, scope and FASTA-file searches, plus the editor's search-hit bookkeeping.

| # | Input | Expected output |
|---|---|---|
| 1 | GGATCCAAAT planted as rc at 1200 | [(1200,1210,-1)]; none with both=False |
| 2 | GRATTC with/without IUPAC; ACGTNCGTAC | (300,306) found / nothing; (700,710) found |
| 3 | probe at 700 with 1 mismatch, fuzzy 0 / 1 | not found / found with 1 mismatch |
| 4 | allowed_mismatches(fuzzy=3) for 4 and 30 bases | 1 and 3 (at most a third) |
| 5 | GATTC across origin of a 19 bp circle | [(17,22)]; search_spans split to [(0,3),(17,19)]; hits_between works |
| 6 | AKW in ATG GCT AAA TGG CAT TAA with translations on/off | [(3,12,'protein',+1)] / [] |
| 7 | CDS gene abcD, db_xref needle | 'abcd' finds (100,400); 'needle' only with details=True |
| 8 | primer 'fwdprimer' at 500..520 | 'fwdprim' finds (500,520); not with primers=False |
| 9 | 'EcoRI' with/without enzymes option | [(300,306,'enzyme')] / [] |
| 10 | selection 250..350 / 400..500 / no selection | (300,306), no note / [] / searches all with 'Nothing is selected' note |
| 11 | FASTA file chrA/chrB, GGTACC; 'kinase' | chrA (200,206), chrB (120,126) one hit per palindrome; non-DNA query returns None |
| 12 | set_search_hits with strands | -1 kept; no strand or both strands -> 0; empty clears |

**Review notes:**

- The reverse-strand protein search with 'AXW' computes hits but asserts nothing.
- IUPAC and fuzzy checks use 'in'/'any', so extra spurious hits would not be noticed.
- One long test mixing engine, editor bookkeeping and file search; a failure early hides the rest.
