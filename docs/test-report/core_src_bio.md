# `core/src/bio.rs`

[← all test files](README.md) · source: [`core/src/bio.rs`](../../core/src/bio.rs)

## `finds_orfs_on_both_strands`

**Result:** ✅ passed · [line 86](../../core/src/bio.rs#L86) · **Area:** ORF finder · **Checked against:** hand-computed · **Strength:** 🟢 strong

**Why it matters:** The ORF finder could show ORFs at wrong coordinates, ignore the minimum length, or miss bottom-strand ORFs.

An ATG..stop ORF is found with exact top-strand coordinates, the minimum length is applied, and the same ORF written on the bottom strand is reported reverse.

| # | Input | Expected output |
|---|---|---|
| 1 | 'CCATGAAAAAATAACC', min 4 codons, linear | [Orf{start 2, end 14, forward}] |
| 2 | same, min 5 codons | empty |
| 3 | revcomp of the sequence, min 4 | [Orf{2, 14, forward=false}] |

**Fixed by the audit of November 2026:**

- 13 Nov 2026: gap closed. The proposed checks were added, pass on the current code and fail under the mutation that slipped through before. Was: The bottom-strand claim ('the same ORF written on the bottom strand is reported reverse' at top-strand 2..14) cannot fail for a wrong coordinate mapping: in 'CCATGAAAAAATAACC' (16 bp) the ORF 2..14 is centred, so its mirror n-end = 2 equals its bottom-strand start a = 2.

## `alternative_starts_and_other_tables`

**Result:** ✅ passed · [line 101](../../core/src/bio.rs#L101) · **Area:** ORF finder / genetic codes · **Checked against:** hand-computed · **Strength:** 🟢 strong

**Why it matters:** The ORF finder could use alternative start codons when they are switched off or from the wrong genetic code, or read TGA as a stop in the Mycoplasma code.

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

**Result:** ✅ passed · [line 120](../../core/src/bio.rs#L120) · **Area:** ORF finder · **Checked against:** hand-computed · **Strength:** 🟢 strong

**Why it matters:** An ORF that crosses the origin of a plasmid could be missed or reported twice.

On a circular sequence an ORF starting near the end and stopping after the origin is found once, with end > length; also on the bottom strand.

| # | Input | Expected output |
|---|---|---|
| 1 | 'AAATAA' + 30 C + 'ATGAAA' (42 bp), min 4, linear | empty |
| 2 | same, circular | [Orf{36, 48, fwd}] |
| 3 | revcomp of the ring, circular | [Orf{36, 48, forward=false}] |
| 4 | 'CCATGAAAAAATAACC', circular | exactly 1 ORF (not duplicated by the doubled scan) |
| 5 | ring 'ATGAAATAA' + 'C'*30 + 'ATGAAA' (45 bp), min 3 codons; and its reverse complement | only (39, 54, +1); on the reverse complement only (36, 51, -1) |

**Review notes:**

- No case for a ring with a start but no stop in a full lap (the 'full lap' guard, bio.rs line 42).

**Fixed by the audit of November 2026:**

- 13 Nov 2026: gap closed. The proposed checks were added, pass on the current code and fail under the mutation that slipped through before. Was: The bottom-strand-across-the-origin claim is symmetric: the ORF covers ring positions 36..41 and 0..5 (6 bases each side of the origin), so the mirrored start (42-48) mod 42 = 36 equals the unmirrored one; a wrong mapping passes. \| The 'full lap without a stop' guard (bio.rs line 42) is not exercised (minor, not claimed).
- 13 Nov 2026: app bug fixed in core/src/bio.rs: on a ring an ATG just after the origin, inside an ORF that started before it, was reported as a second, shorter ORF with the same stop ('ATGAAATAA' + 30 C + 'ATGAAA' gave (0,9) and (39,54)). Now one stop gives one ORF, the longest, on both strands; the test checks that ring exactly and random rings for equality with the hand scan.
