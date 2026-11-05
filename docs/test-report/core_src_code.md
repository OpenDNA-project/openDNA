# `core/src/code.rs`

[← all test files](README.md) · source: [`core/src/code.rs`](../../core/src/code.rs)

## `every_table_has_64_codons_and_known_ids`

**Result:** ✅ passed · [line 330](../../core/src/code.rs#L330) · **Area:** Genetic codes · **Checked against:** property/invariant · **Strength:** 🔴 weak

Every built-in NCBI table has 64 amino-acid and 64 start letters; ids 1 and 11 exist and context-dependent table 27 is absent.

| # | Input | Expected output |
|---|---|---|
| 1 | all TABLES | aas.len() == 64 and starts.len() == 64 for each |
| 2 | by_id(1), by_id(11), by_id(27) | Some, Some, None |

**Review notes:**

- Only lengths are checked; the module doc (code.rs line 11) says 'the tests compare every codon' against Biopython: that is done in Python (tests/test_genetic_code.py), not here, so these Rust tests alone would not catch a wrong letter in a table.

## `the_standard_code_and_the_known_differences`

**Result:** ✅ passed · [line 340](../../core/src/code.rs#L340) · **Area:** Genetic codes · **Checked against:** hand-computed · **Strength:** 🟡 adequate

Standard translation works and known table differences hold (TGA=W in table 4, AGA=stop in table 2); lower case and RNA letters are accepted.

| # | Input | Expected output |
|---|---|---|
| 1 | standard.translate('ATGGCCTAAGNA') | 'MA*X' |
| 2 | table 4 aa('TGA'); table 1 aa('TGA') | 'W'; '*' |
| 3 | table 2 aa('AGA'); table 1 aa('AGA') | '*'; 'R' |
| 4 | table 1 aa('aug') | 'M' |

## `a_cds_is_translated_with_its_own_table_and_frame`

**Result:** ✅ passed · [line 353](../../core/src/code.rs#L353) · **Area:** CDS translation · **Checked against:** hand-computed · **Strength:** 🟡 adequate

A CDS is translated with its /transl_table (else the app default) and /codon_start; an internal stop or an incomplete last codon on a non-partial CDS is an error.

| # | Input | Expected output |
|---|---|---|
| 1 | 'ATGTGAGTGTAA', CDS 1..12, default table 1 | error (internal TGA stop) |
| 2 | same with /transl_table=4; plain CDS with default 4; table-4 CDS with default 11 | 'MWV' in all three |
| 3 | /codon_start=2 on 1..12 | error |
| 4 | /codon_start=2 on 1..>12 | 'CEC' |

**Review notes:**

- Error cases use is_err() only (lines 366, 375): the error kind/message is not checked, so an error for a different reason would pass.

## `an_alternative_start_codon_is_read_as_methionine`

**Result:** ✅ passed · [line 372](../../core/src/code.rs#L372) · **Area:** CDS translation · **Checked against:** hand-computed · **Strength:** 🟢 strong

A table start codon other than ATG is translated as M only at the true 5' start of a complete CDS.

| # | Input | Expected output |
|---|---|---|
| 1 | 'TTGAAATAA', CDS 1..9, table 1 | 'MK' |
| 2 | 'ATGTTGTAA', CDS 1..9 | 'ML' (internal TTG = L) |
| 3 | 'TTGAAATAA', CDS <1..9 with codon_start 1 | 'LK' |

## `bottom_strand_cds_and_transl_except`

**Result:** ✅ passed · [line 383](../../core/src/code.rs#L383) · **Area:** CDS translation · **Checked against:** hand-computed · **Strength:** 🟡 adequate

A complement CDS is read from the bottom strand, /transl_except (Sec) is applied, and an unsupported exception is an error rather than a silent translation.

| # | Input | Expected output |
|---|---|---|
| 1 | top = revcomp('ATGTGATAA'); complement(1..9) with transl_except (pos:complement(4..6),aa:Sec) | 'MU' |
| 2 | same without the exception | error (internal stop) |
| 3 | transl_except (pos:complement(4..5),aa:TERM) | error |

**Review notes:**

- Error cases are is_err() only; the reason is not asserted.

## `joined_cds_reads_its_parts_in_order`

**Result:** ✅ passed · [line 396](../../core/src/code.rs#L396) · **Area:** CDS translation · **Checked against:** hand-computed · **Strength:** 🟡 adequate

A join() CDS is translated from its parts concatenated in order.

| # | Input | Expected output |
|---|---|---|
| 1 | 'ATGCCCAAATAA', CDS join(1..3,7..12) | 'MK' |

**Review notes:**

- Parts are in ascending order; a complement(join) or out-of-order join is not tested.

## `pseudo_and_partial_cds_are_skipped`

**Result:** ✅ passed · [line 403](../../core/src/code.rs#L403) · **Area:** CDS translation · **Checked against:** hand-computed · **Strength:** 🟢 strong

/pseudo CDSs and 5'-partial CDSs without /codon_start are skipped with the right reason; /codon_start makes a partial CDS translatable.

| # | Input | Expected output |
|---|---|---|
| 1 | 'ATGAAATAA', 1..9 with /pseudo | Err(Skip::Pseudo) |
| 2 | <1..9, no codon_start | Err(Skip::Partial) |
| 3 | <1..9 with codon_start=1 | 'MK' |
