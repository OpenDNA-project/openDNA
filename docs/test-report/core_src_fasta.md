# `core/src/fasta.rs`

[← all test files](README.md) · source: [`core/src/fasta.rs`](../../core/src/fasta.rs)

## `a_fasta_rewritten_in_the_same_clock_tick_is_re_indexed`

**Result:** ✅ passed · [line 316](../../core/src/fasta.rs#L316) · **Area:** Indexed FASTA (.fai) · **Checked against:** hand-computed · **Strength:** 🟢 strong

**Why it matters:** After a genome FASTA is rewritten quickly, a stale .fai could be trusted and fetches would return the old bases.

A FASTA rewritten with the same mtime as its saved .fai is still re-indexed (stale index not trusted), so fetch returns the new bases.

| # | Input | Expected output |
|---|---|---|
| 1 | '>r\nACGTACGTACGTACGTACGT\n' opened (writes .fai) | length 20 |
| 2 | rewrite to '>r\nACGTACGTAC\n' and set its mtime equal to the .fai's; reopen | length 10; fetch('r',0,10) == 'ACGTACGTAC' |

**Review notes:**

- Fixed file names in the shared system temp dir (ge_test_same_tick.fa, line 319); concurrent test runs can collide.

## `fetch_across_line_breaks_with_crlf`

**Result:** ✅ passed · [line 334](../../core/src/fasta.rs#L334) · **Area:** Indexed FASTA (.fai) · **Checked against:** hand-computed · **Strength:** 🟢 strong

**Why it matters:** Fetching a region from a Windows (CRLF) or wrapped genome FASTA could return bases shifted by the line endings.

CRLF line endings and wrapped lines are indexed correctly; fetch crosses line breaks, clamps the end, errors on unknown names; the saved .fai reads back identically.

| # | Input | Expected output |
|---|---|---|
| 1 | '>chrA desc\r\nACGT\r\nTTGG\r\nCA\r\n>chrB\r\nGGGG\r\n' | lengths 10 and 4 |
| 2 | fetch chrA 2..9; chrA 8..999; chrB 0..4; 'nope' | 'GTTTGGC'; 'CA'; 'GGGG'; error |
| 3 | read_fai of the written index | == entries |

## `replace_one_sequence_and_search`

**Result:** ✅ passed · [line 349](../../core/src/fasta.rs#L349) · **Area:** Indexed FASTA (.fai) · **Checked against:** hand-computed · **Strength:** 🟢 strong

**Why it matters:** Saving one edited record into a multi-record FASTA could damage the other records, or motif search in a FASTA could miss or duplicate hits.

Exact motif search reports both strands once for palindromes; write_replacing swaps one record's bases and leaves other records and headers untouched; an unknown name errors.

| # | Input | Expected output |
|---|---|---|
| 1 | records a 'ACGTACGTACGT', b 'TTTT', c 'GGGAATTCGG'; search c 'gaattc' | [(2, +1)] only |
| 2 | search a 'CGTA' | [(1,1),(3,-1),(5,1),(7,-1)] |
| 3 | write_replacing(b -> 'CC' x 40) to another file | b == new bases; a and c unchanged; file starts with '>a first\n' |
| 4 | write_replacing name 'zzz' | error |

**Review notes:**

- Wrapping at 60 (80 new bases should be 60+20 lines) and the in-place case dst == src, both promised by write_replacing's doc, are not tested.

## `rejects_uneven_wrapping`

**Result:** ✅ passed · [line 366](../../core/src/fasta.rs#L366) · **Area:** Indexed FASTA (.fai) · **Checked against:** hand-computed · **Strength:** 🟡 adequate

**Why it matters:** An unevenly wrapped FASTA could be indexed anyway, so every fetch from it would silently return wrong bases.

A FASTA whose lines are wrapped unevenly inside a record is rejected, because offset arithmetic would return wrong bases.

| # | Input | Expected output |
|---|---|---|
| 1 | '>x\nACG\nACGT\n' | build_fai error |

**Fixed by the audit of November 2026:**

- 13 Nov 2026: gap closed. The proposed checks were added, pass on the current code and fail under the mutation that slipped through before. Was: Only one kind of unevenness is tested ('>x\nACG\nACGT\n': a line longer than the first). The common case - a short line in the middle of a record ('ACGT/AC/ACGT') - and same base count with a different line width (mixed LF/CRLF) are never fed, so those two rejection branches are unverified.
