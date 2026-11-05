# `tests/test_multi_genbank.py`

[← all test files](README.md) · source: [`tests/test_multi_genbank.py`](../../tests/test_multi_genbank.py)

## `test_records_and_saving`

**Result:** ✅ passed in 1.82 s · [line 44](../../tests/test_multi_genbank.py#L44) · **Area:** Multi-record GenBank I/O · **Checked against:** hand-computed · **Strength:** 🟢 strong

A GenBank file with two records lists them largest first and switches between them; Save writes the edited record back without touching the other; Save As keeps both records; switching keeps tabs and unsaved edits, and one Save writes both.

| # | Input | Expected output |
|---|---|---|
| 1 | two.gb = pBAD30_eilR.gb + pFG051_Chlor.gb; load | is_multi_genbank true (false for the single plasmid); 2 records sorted by size desc |
| 2 | record_requested(second) | editor shows second record, length matches, picker text updated |
| 3 | insert 'ACGTACGT' at 10 in first record, Save | not modified; 2 LOCUS lines; first == old[:10]+'ACGTACGT'+old[10:]; second unchanged; first name once; no *.saving file |
| 4 | Save As copy.gb (dialog patched) | both files are multi-GenBank; editor path == copy.gb |
| 5 | insert 'GGGG' at 5, switch to second and back | tab layout and current tab kept; first record doc still modified; sequence == before[:5]+'GGGG'+before[5:] |
| 6 | insert 'TT' at 7 in second, Save | file has both edits; first record no longer modified |

**Review notes:**

- Lines 50, 63, 79, 89: four MainWindows are created and never closed (two hold unsaved edits); they stay alive for later tests in the session.
- Line 290: Save As only checks that both files are multi-GenBank, not that the copy's records match the originals (content is checked in test_file_history).

## `test_search_primers_and_detect_cover_every_record`

**Result:** ✅ passed in 1.28 s · [line 122](../../tests/test_multi_genbank.py#L122) · **Area:** Multi-record GenBank: search, primers, Detect Features · **Checked against:** hand-computed · **Strength:** 🟡 adequate

Searching a multi-record file finds DNA and feature names in records that are not open (opening them), primers are matched against other records, and Detect Features offers every other record.

| # | Input | Expected output |
|---|---|---|
| 1 | two-record file; search for 24 bp seq[1000:1024] of the second record | hits include the second record; editor switches to it |
| 2 | search for the label of a CDS of the first record (not open) | a 'feature' hit in that record |
| 3 | extra primer 'probe' = other record seq[2000:2022] | ed.elsewhere has record==other, label 'probe', start 2000; primer tree mentions 'other sequences' |
| 4 | search 'probe' | a 'primer' hit in the other record |
| 5 | DetectDialog | last column 'Record'; others == all records except the open one |

**Review notes:**

- Lines 33-41 and 157-161: polling loops with time.sleep(0.02) and 8 s / 6 s timeouts that end silently; on a slow machine the asserts that follow fail without saying the wait timed out.
- Line 133: garbled comment ('only the second (open: no) record is sure to have'); the 24 bp probe is not checked to be absent from the first record.
- Line 138: DNA hit coordinates are not checked (only the record name).

## `test_vibrio_cholerae_two_chromosomes`

**Result:** ⏭️ skipped · [line 177](../../tests/test_multi_genbank.py#L177) · **Area:** Multi-record GenBank on a real genome · **Checked against:** regression value · **Strength:** 🔴 weak

> no V. cholerae file at V:/HM_Softwares/TestData/vibrio_cholerae/vibrio_cholerae_N16961.gb

The two-chromosome V. cholerae N16961 GenBank opens with both records (sizes exact), switches to chromosome II, and a 'toxin' search finds hits on both chromosomes.

| # | Input | Expected output |
|---|---|---|
| 1 | vibrio_cholerae_N16961.gb (OPENDNA_VIBRIO_GB or V:/ path) | fasta_records == [('NC_002505', 2961149), ('NC_002506', 1072315)], circular, > 1000 features |
| 2 | switch to NC_002506; search 'toxin' | length 1072315; hit records == {NC_002505, NC_002506} |

**Review notes:**

- Lines 178-179: runtime pytest.skip when the file is absent; the default is a Windows path (V:/...) outside git, so it is always skipped in CI/Linux (skipped in the lead's run).
- No docstring.
