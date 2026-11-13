# `tests/test_repeats_msa.py`

[← all test files](README.md) · source: [`tests/test_repeats_msa.py`](../../tests/test_repeats_msa.py)

## `test_repeats_are_found_where_planted_and_not_in_random_dna`

**Result:** ✅ passed in 0.55 s · [line 36](../../tests/test_repeats_msa.py#L36) · **Area:** Repeat finder · **Checked against:** hand-computed · **Strength:** 🟡 adequate

**Why it matters:** Without it, the repeat finder could miss planted direct, inverted or tandem repeats, flood random DNA with chance hits, or hang on long low-complexity input.

Plants repeats in random DNA and checks repeats.find_repeats and pair_repeats return them at their places (within 2 bases), that random DNA gives nothing at the defaults, thresholds/mismatches/origin/capping behave, and pair_dots returns both dot sets.

| # | Input | Expected output |
|---|---|---|
| 1 | 30 bp x 3 direct, 16 bp inverted pair, GATTACAGCT x 5 | one direct (length 30-32, 3 copies), one inverted (16-18, strands +1/-1), one tandem (unit 10, 5 copies); not capped |
| 2 | random 20 kb linear, 5 kb circular | no repeats |
| 3 | random 20 kb with max_loop=None | some repeats (chance stems) |
| 4 | 14 bp direct pair | not listed by default; listed with min_direct=12 |
| 5 | 30 bp pair with one mismatch | not listed by default; listed with mismatches=1 |
| 6 | 25 bp direct over the origin of a circle | one direct, copies near 115 and 240 |
| 7 | ACGT x 200 000 | capped with 'skipped' in message |
| 8 | two sequences sharing 40 bp direct and 30 bp inverted | direct near (50,100,40), inverted near (150,30,30); two random 3 kb: none |
| 9 | sequence against itself | the repeated pair once, not the main diagonal |
| 10 | pair_dots | non-empty direct and inverted |

**Review notes:**

- Positions are checked within 2 bases and lengths within a range, so small off-by-one errors pass.
- pair_dots is checked only for being non-empty.
- In the two-sequence check `near` uses zip without a length check.

## `test_multiple_alignment_consensus_identity_and_parsing`

**Result:** ✅ passed in 0.05 s · [line 99](../../tests/test_repeats_msa.py#L99) · **Area:** Multiple sequence alignment · **Checked against:** hand-computed · **Strength:** 🟡 adequate

**Why it matters:** Without it, the built-in multiple aligner could lose bases, misplace a gap, or give wrong consensus, identity and column/position mapping in the alignment view.

Aligns four copies of a 200 bp sequence (two substitutions, one 3-base deletion) and checks rows, gap placement, consensus, identities and conservation; checks column/position mapping on a 2-row toy; checks aligned-FASTA parsing and refusals.

| # | Input | Expected output |
|---|---|---|
| 1 | ref, sub50, sub120, del80 (200 bp) | rows without gaps equal inputs; width 200; del80 row has '---' at 80..83; row 0 is the reference |
| 2 | consensus; identity rows 0,1,3 | consensus = reference; 100, 100*199/200, 100*197/200 |
| 3 | reference None | identities() equal identity(i, 0) |
| 4 | conservation at 0, 50, 81 | 1.0, 0.75, 0.75 |
| 5 | AC--GT vs ACTTGT | positions [0,1,2,2,2,3]; column_of(0,2)=4, column_of(1,2)=2 |
| 6 | >a AC- / GT; >b desc ACTG / T- | names ['a','b desc'], rows 'AC-GT', 'ACTGT-' |
| 7 | one sequence; sequence over MAX_BUILTIN | ValueError |

**Review notes:**

- identity(1) for sub120 (row 2) is not checked.
- Parsing accepts rows of different widths (5 and 6) without comment; the test does not say whether that should be refused.
