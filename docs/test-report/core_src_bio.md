# `core/src/bio.rs`

[← all test files](README.md) · source: [`core/src/bio.rs`](../../core/src/bio.rs)

## `finds_orfs_on_both_strands`

**Result:** ✅ passed · [line 77](../../core/src/bio.rs#L77) · **Area:** ORF finder · **Checked against:** hand-computed · **Strength:** 🟢 strong

An ATG..stop ORF is found with exact top-strand coordinates, the minimum length is applied, and the same ORF written on the bottom strand is reported reverse.

| # | Input | Expected output |
|---|---|---|
| 1 | 'CCATGAAAAAATAACC', min 4 codons, linear | [Orf{start 2, end 14, forward}] |
| 2 | same, min 5 codons | empty |
| 3 | revcomp of the sequence, min 4 | [Orf{2, 14, forward=false}] |

## `alternative_starts_and_other_tables`

**Result:** ✅ passed · [line 88](../../core/src/bio.rs#L88) · **Area:** ORF finder / genetic codes · **Checked against:** hand-computed · **Strength:** 🟢 strong

Alternative start codons count only when enabled and only if the chosen table lists them; table 4 reads TGA as Trp.

| # | Input | Expected output |
|---|---|---|
| 1 | 'CCGTGAAAAAATAACC', table 11, alt starts off / on | empty / [Orf{2,14,fwd}] |
| 2 | same, table 1, alt on | empty (GTG not a table-1 start) |
| 3 | 'CCTTGAAAGTGAAATAACC', table 1, alt on | [Orf{2,17,fwd}] (longest ORF from first start TTG) |
| 4 | 'ATGTGAAAATAA', min 4, table 4 vs table 1 | 1 ORF vs empty |

**Review notes:**

- Table-4 case (line 102) checks only len()==1, not the coordinates 0..12.

## `finds_an_orf_that_crosses_the_origin_of_a_ring`

**Result:** ✅ passed · [line 107](../../core/src/bio.rs#L107) · **Area:** ORF finder · **Checked against:** hand-computed · **Strength:** 🟢 strong

On a circular sequence an ORF starting near the end and stopping after the origin is found once, with end > length; also on the bottom strand.

| # | Input | Expected output |
|---|---|---|
| 1 | 'AAATAA' + 30 C + 'ATGAAA' (42 bp), min 4, linear | empty |
| 2 | same, circular | [Orf{36, 48, fwd}] |
| 3 | revcomp of the ring, circular | [Orf{36, 48, forward=false}] |
| 4 | 'CCATGAAAAAATAACC', circular | exactly 1 ORF (not duplicated by the doubled scan) |

**Review notes:**

- No case for a ring with a start but no stop in a full lap (the 'full lap' guard, bio.rs line 42).
