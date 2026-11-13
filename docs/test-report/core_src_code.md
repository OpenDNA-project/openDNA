# `core/src/code.rs`

[← all test files](README.md) · source: [`core/src/code.rs`](../../core/src/code.rs)

## `a_registered_code_is_found_by_id_and_a_removed_codon_reads_x`

**Result:** ✅ passed · [line 357](../../core/src/code.rs#L357) · **Area:** Genetic codes · **Checked against:** hand-computed · **Strength:** 🟡 adequate

**Why it matters:** A recoded organism's or user-made genetic code could fail to register, be looked up wrongly, or translate a reassigned codon as its old amino acid instead of X.

register() refuses ids below 100 and strings that are not 64 letters, stores a valid code under its id with its name; by_id finds it and translates with it (a codon written 'X' reads as X); an unknown id is None in by_id and the standard code in by_id_or_standard. Every built-in NCBI table has 64 amino-acid and 64 start letters; ids 1 and 11 exist and context-dependent table 27 is absent.

| # | Input | Expected output |
|---|---|---|
| 1 | register(5, 'too low', table-11 aas with TCG and TAG set to X, table-11 starts) | Err (ids below 100 belong to NCBI) |
| 2 | register(100, 'short', 'FFLL', '----') | Err (not 64 letters) |
| 3 | register(100, 'Test recoded', modified table 11); by_id(100) | Ok; name == 'Test recoded' |
| 4 | code 100 translates 'ATGTCGTAGTCATAA' | 'MXXS*' (TCG and TAG read X, TCA still S, TAA stop) |
| 5 | by_id(101); by_id_or_standard(101) | None; standard code (id 1) |
| 6 | all TABLES | aas.len() == 64 and starts.len() == 64 for each |
| 7 | by_id(1), by_id(11), by_id(27) | Some, Some, None |

**Review notes:**

- (was a_registered_code_is_found_by_id_and_a_removed_codon_reads_x) Error cases are is_err() only; the reason is not asserted.
- (was a_registered_code_is_found_by_id_and_a_removed_codon_reads_x) Replacing an existing custom id (documented 'add (or replace)') and the starts string of a registered code are not tested.
- (was a_registered_code_is_found_by_id_and_a_removed_codon_reads_x) The Python binding register_genetic_code (lib.rs) is not called by any Python test; this is the only test of the code registry.
- (was every_table_has_64_codons_and_known_ids) Only lengths are checked; the module doc says 'the tests compare every codon' against Biopython: that is done in Python, not here, so these Rust tests alone would not catch a wrong letter in a table.
- (was every_table_has_64_codons_and_known_ids) Largely superseded by tests/test_genetic_code.py::test_translation_tables_match_biopython, which compares every codon of every table (a wrong length would fail there too); this only adds a fast length/id smoke check.

## `the_standard_code_and_the_known_differences`

**Result:** ✅ passed · [line 380](../../core/src/code.rs#L380) · **Area:** Genetic codes · **Checked against:** hand-computed · **Strength:** 🟡 adequate

**Why it matters:** Translation could get basic codons wrong, ignore table-specific differences (TGA=W in table 4, AGA=stop in table 2), or reject lower-case or RNA input.

Standard translation works and known table differences hold (TGA=W in table 4, AGA=stop in table 2); lower case and RNA letters are accepted.

| # | Input | Expected output |
|---|---|---|
| 1 | standard.translate('ATGGCCTAAGNA') | 'MA*X' |
| 2 | table 4 aa('TGA'); table 1 aa('TGA') | 'W'; '*' |
| 3 | table 2 aa('AGA'); table 1 aa('AGA') | '*'; 'R' |
| 4 | table 1 aa('aug') | 'M' |

## `a_cds_is_translated_with_its_own_table_and_frame`

**Result:** ✅ passed · [line 393](../../core/src/code.rs#L393) · **Area:** CDS translation · **Checked against:** hand-computed · **Strength:** 🟡 adequate

**Why it matters:** A CDS could be translated with the app default instead of its own /transl_table or /codon_start, or a broken CDS (internal stop, incomplete codon) could be translated silently.

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

**Result:** ✅ passed · [line 412](../../core/src/code.rs#L412) · **Area:** CDS translation · **Checked against:** hand-computed · **Strength:** 🟢 strong

**Why it matters:** A CDS starting with an alternative start such as GTG or TTG could show V or L instead of M as its first residue, or internal GTGs could wrongly become M.

A table start codon other than ATG is translated as M only at the true 5' start of a complete CDS.

| # | Input | Expected output |
|---|---|---|
| 1 | 'TTGAAATAA', CDS 1..9, table 1 | 'MK' |
| 2 | 'ATGTTGTAA', CDS 1..9 | 'ML' (internal TTG = L) |
| 3 | 'TTGAAATAA', CDS <1..9 with codon_start 1 | 'LK' |

## `bottom_strand_cds_and_transl_except`

**Result:** ✅ passed · [line 423](../../core/src/code.rs#L423) · **Area:** CDS translation · **Checked against:** hand-computed · **Strength:** 🟡 adequate

**Why it matters:** A complement CDS could be translated from the wrong strand, or a /transl_except (e.g. selenocysteine) could be ignored.

A complement CDS is read from the bottom strand, /transl_except (Sec) is applied, and an unsupported exception is an error rather than a silent translation.

| # | Input | Expected output |
|---|---|---|
| 1 | top = revcomp('ATGTGATAA'); complement(1..9) with transl_except (pos:complement(4..6),aa:Sec) | 'MU' |
| 2 | same without the exception | error (internal stop) |
| 3 | transl_except (pos:complement(4..5),aa:TERM) | error |

**Review notes:**

- Error cases are is_err() only; the reason is not asserted.

**Fixed by the audit of November 2026:**

- 13 Nov 2026: gap closed. The proposed checks were added, pass on the current code and fail under the mutation that slipped through before. Was: 'An unsupported exception is an error rather than a silent translation' is tested with (pos:complement(4..5),aa:TERM), but that input is rejected anyway by the reading-frame check (first base 4 sits at offset 4, not a codon start), so the whole-codon rule itself is never exercised; is_err() hides which check fired.

## `joined_cds_reads_its_parts_in_order`

**Result:** ✅ passed · [line 442](../../core/src/code.rs#L442) · **Area:** CDS translation · **Checked against:** hand-computed · **Strength:** 🟡 adequate

**Why it matters:** A spliced (join) CDS could be translated from its parts in the wrong order.

A join() CDS is translated from its parts concatenated in order.

| # | Input | Expected output |
|---|---|---|
| 1 | 'ATGCCCAAATAA', CDS join(1..3,7..12) | 'MK' |

**Review notes:**

- Parts are in ascending order; a complement(join) or out-of-order join is not tested.

**Fixed by the audit of November 2026:**

- 13 Nov 2026: gap closed. The proposed checks were added, pass on the current code and fail under the mutation that slipped through before. Was: 'Translated from its parts concatenated in order' is tested only with join(1..3,7..12), whose written order equals ascending order, so reading parts sorted (instead of as written) passes; complement(join) reading order is not covered either.

## `pseudo_and_partial_cds_are_skipped`

**Result:** ✅ passed · [line 452](../../core/src/code.rs#L452) · **Area:** CDS translation · **Checked against:** hand-computed · **Strength:** 🟢 strong

**Why it matters:** Pseudogenes and 5'-partial CDSs could get a nonsense /translation instead of being skipped with a reason.

/pseudo CDSs and 5'-partial CDSs without /codon_start are skipped with the right reason; /codon_start makes a partial CDS translatable.

| # | Input | Expected output |
|---|---|---|
| 1 | 'ATGAAATAA', 1..9 with /pseudo | Err(Skip::Pseudo) |
| 2 | <1..9, no codon_start | Err(Skip::Partial) |
| 3 | <1..9 with codon_start=1 | 'MK' |
