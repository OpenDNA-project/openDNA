# `tests/test_multi_record.py`

[← all test files](README.md) · source: [`tests/test_multi_record.py`](../../tests/test_multi_record.py)

## `test_one_record_is_read_and_saved_in_place`

**Result:** ✅ passed in 0.12 s · [line 45](../../tests/test_multi_record.py#L45) · **Area:** Multi-record files · **Checked against:** property/invariant · **Strength:** 🟢 strong

**Why it matters:** Without it, editing one chromosome or plasmid inside a multi-record GenBank file could corrupt the other records or leave the file half-written when the save fails.

Builds a three-record GenBank from the demo plasmids and checks detection, listing, reading one record, saving an edited record leaving the others byte-identical, and that a failed final rename leaves the file unchanged. Also checks track contig-name matching, the note shown when saving one FASTA record alone, and the real Vibrio two-chromosome file.

| # | Input | Expected output |
|---|---|---|
| 1 | three.gb (pBAD30, pFG051, copy of pBAD30 renamed 'copy'); single-record pBAD30 | is_multi_genbank true for three.gb, false for the plasmid; records list names and lengths in file order |
| 2 | read_record(b); read_record('nope') | sequence equals pFG051; OSError 'no record named' |
| 3 | Insert ACGTACGT at 10 in record a, save_document | record a reads back edited; records b and c bytes identical; order unchanged; no *.saving left |
| 4 | fileio.replace patched to raise PermissionError, save again | PermissionError; file bytes unchanged; no *.saving |
| 5 | match_contig cases | chr2->chr2; ['1','2'],chr2->'2'; ['Chr1','Chr2'],chr1->'Chr1'; ['1','chr1'],chrX->None; ['chr1','1'],chr1->'chr1'; ['chr7'],chr1->None |
| 6 | other_records_note for chr1 of a 3-record FASTA; with count 1 | starts ' (only chr1: the other 2 sequences of g.fa'; '' |
| 7 | vibrio_cholerae_N16961.gb | multi-record; [(NC_002505, 2961149), (NC_002506, 1072315)] |

**Review notes:**

- Skips entirely if the demo plasmids are missing.
- match_contig and other_records_note are main-window helpers unrelated to in-place saving; they would fit better in a tracks or FASTA-save test.

## `test_records_are_added_reordered_and_removed_byte_for_byte`

**Result:** ✅ passed in 0.14 s · [line 101](../../tests/test_multi_record.py#L101) · **Area:** Multi-record files · **Checked against:** property/invariant · **Strength:** 🟢 strong

**Why it matters:** Without it, reordering, removing or adding records in a GenBank or FASTA file could alter untouched records, break CRLF/BOM files, reuse a stale index, or create duplicate record names.

Checks reorder/remove of GenBank and FASTA records keeps each record byte-identical and refuses bad orders; free-name generation; appending a plasmid (and a .dna with primers) to a GenBank file keeps the original bytes; appending FASTA records uses the file's line endings, ends the old last line and rebuilds the index.

| # | Input | Expected output |
|---|---|---|
| 1 | three.gb rewritten as [c,a,b], then [c,b] | names in that order; each record chunk byte-identical to before |
| 2 | Orders [], [c,c], [c,'nope'] | OSError each; file still [c,b]; no *.saving |
| 3 | FASTA with BOM, CRLF, no final newline, existing .fai; reorder to chr3,chr1,chr2 | exact bytes '\xef\xbb\xbf>chr3\r\nTT\r\n>chr1 one\r\nACGT\r\nAC\r\n>chr2\r\nGGGG\r\n'; index lists new order; chr1 = 'ACGTAC' |
| 4 | Rewrite FASTA to [chr2] | {chr2: 4} |
| 5 | record_name(' chr 1 '), ('', ['sequence']), ('a', ['a','a_2']) | 'chr_1', 'sequence_2', 'a_3' |
| 6 | Append pBAD30 to a one-record copy of itself | new name '<name>_2'; file starts with the original bytes; new record's sequence equals the plasmid; source document name unchanged |
| 7 | Append the .dna plasmid as 'from_dna' with its primers as primer_bind | three records; from_dna has primer_bind features |
| 8 | CRLF genome.fa (chrA 120, chrB 40, no final newline, stale .fai) + p1 and chrA from more.fasta (p2 unticked) | {chrA:120, chrB:40, p1:70, chrA_2:8}; only CRLF endings; p1 wrapped 60+10; chrB still 'GGCC'*10 (last line ended, not joined) |

**Review notes:**

- The .dna-with-primers part only checks that some primer_bind exists, not names or positions.
- The test builds the primer_bind list itself, so it does not test how the app's Add Record dialog builds it.

**Fixed by the audit of November 2026:**

- 13 Nov 2026: gap closed. The proposed checks were added, pass on the current code and fail under the mutation that slipped through before. Was: 'a .dna source brings its primers as primer_bind features' is asserted as any(kind == 'primer_bind'): one primer out of 15 satisfies it; names, locations and sequences are not compared.
