# `tests/test_variants.py`

[← all test files](README.md) · source: [`tests/test_variants.py`](../../tests/test_variants.py)

## `test_variants_are_written_into_the_bases_or_annotated_beside_them`

**Result:** ✅ passed in 0.01 s · [line 39](../../tests/test_variants.py#L39) · **Area:** Read variants applied to the sequence · **Checked against:** in-test reference model · **Strength:** 🟢 strong

**Why it matters:** Without it, applying variants found in sequencing reads could corrupt the user's sequence - wrong base, shifted indel, several edits interfering - and the damage would be hard to undo.

Checks the label given to each kind of change and which can be applied at all, that applying one gives exactly the sequence a plain string edit gives (first base, last base, past the end), that it is one undo step with the right follow-up message, that several variants are written last-position-first in any chosen order, that a locked sequence refuses everything, and that annotating writes one correctly placed `variation` feature per variant instead.

| # | Input | Expected output |
|---|---|---|
| 1 | change_text for snp A>G, 5 bp del, 11 bp del, insertion TT, insertion '' | 'A > G'; 'deletion of 5 bp (ACGTA)'; 'deletion of 11 bp' (not spelled out); 'insertion of TT'; 'insertion of ?' |
| 2 | can_apply for snp, del, ins TT, ins '', ins TNT | True, True, True, False, False |
| 3 | each of snp at 5, del of 4 at 8, ins CCC at 8, snp at 0, del of the last base, ins after the last base | returns 1; sequence equals the corresponding plain string edit; exactly one undo step |
| 4 | the BAM entry after each apply | cache cleared; variants None for a SNP (coordinates unchanged), [] for an indel; the 'map the reads again' message appears only for indels |
| 5 | four variants (ins 20, snp 2, del 10, ins 25 with N) applied in both orders | count 3 (the N insertion skipped); the same expected sequence both ways; three undos restore the original and leave three redo steps |
| 6 | apply with dna_locked = True | returns 0, sequence unchanged, no history at all, entry's variants and cache untouched, a 'locked' message emitted |
| 7 | apply with an empty list / only an unsettled insertion | returns 0, nothing changed |
| 8 | variant_feature for snp at 5 (7 of 10 reads), del 8-12, ins at 8, ins at 0, ins at 32, snp with 0 of 0 reads | ('variation','6') with label 'C > T', replace 'T', note '7 of 10 reads (70 %) in sample.bam'; '9..12' with replace ''; '8..9'; '1..1'; '32..32'; note '0 of 0 reads (0 %)' |
| 9 | annotate_variants with a SNP and a deletion, then undo | returns 2, features at '6' and '9..12' covering the reference bases, sequence untouched, one undo removes both; an empty list returns 0 |

**Review notes:**

- One long function covering apply, ordering, locking and annotation; splitting it would make failures easier to read.
