# `tests/test_input_validation.py`

[← all test files](README.md) · source: [`tests/test_input_validation.py`](../../tests/test_input_validation.py)

## `test_sequence_field_paste`

**Result:** ✅ passed in 0.00 s · [line 45](../../tests/test_input_validation.py#L45) · **Area:** Input validation: DNA / protein fields · **Checked against:** hand-computed · **Strength:** 🟡 adequate

Pasting into a DNA field strips FASTA headers, numbers and line breaks or refuses the whole paste; protein pastes are cleaned the same way, so users can paste raw text safely.

| # | Input | Expected output |
|---|---|---|
| 1 | paste '>seq1 some record\n1 acgt acgt\n9 ACGT\n' into dna_field | text == 'acgtacgtACGT' |
| 2 | paste 'ACGT is not DNA: hello' | text == '' (refused whole) |
| 3 | insert 'ACGU xyz' through the validator | no 'x' in the text |
| 4 | clean_pasted_protein('>p\nMKT AYI\n10 AKQ'); clean_pasted_protein('MKT123?!') | 'MKTAYIAKQ'; raises ValueError |
| 5 | paste '>p\nMKTAYI\nAKQ' into protein_field | text == 'MKTAYIAKQ' |

**Review notes:**

- Line 53: comment promises the refusal comes 'with a tooltip'; the tooltip is not checked.
- Line 55: only checks that 'x' is absent; what the field holds after inserting 'ACGU xyz' (e.g. whether U is kept or the insert refused whole) is not asserted.

## `test_primer_dialogs`

**Result:** ✅ passed in 0.89 s · [line 67](../../tests/test_input_validation.py#L67) · **Area:** Input validation: primer dialogs · **Checked against:** hand-computed · **Strength:** 🟡 adequate

Bad primer names/sequences give a named problem and save nothing; fields are length-limited; PCR conditions set to 0 are clamped so Tm stays positive; the options dialog refuses a missing genome file.

| # | Input | Expected output |
|---|---|---|
| 1 | primer_problem: ('', 'ACGTACGTAC'), ('a\tb', ...), ('ok_short','ACGTAC'), ('ok_long','A'*201), ('ok_bad','ACGTACGTXX') | message contains 'name', 'control character', 'at least 8', 'at most', 'not DNA' |
| 2 | PrimerDialog on selection 10-40; name set to 'x'*500 | maxLength == MAX_NAME / MAX_PRIMER; name truncated to MAX_NAME |
| 3 | accept with name '   ' | warning mentions 'name', dialog not Accepted |
| 4 | name 'good_primer', sequence 'ACGT' | warning 'at least 8', primer not in db |
| 5 | sequence 'ACGT ACGT\nACGT' | db has good_primer with sequence 'ACGTACGTACGT' |
| 6 | PrimerDialog with selection 0-400 | sequence box empty |
| 7 | settings pcr/na_mm and pcr/primer_nm set to 0 | conditions()[key] >= 1.0 and tm('ATGCATGCATGCATGCATGC', **conditions) > 0; option spin boxes minimum >= 1 |
| 8 | OptionsDialog, genome check on, path tmp/nope.fa, accept | warning contains 'no file', not Accepted |

**Review notes:**

- Line 104: `tm(...) > 0` only rules out a negative/zero Tm; a nonsense but positive value (e.g. 500 °C) would pass.
- Line 22: fixtures/helpers (open_doc, close, dna) are imported from another test module (test_cas9), coupling the two files.

## `test_feature_and_new_sequence_dialogs`

**Result:** ✅ passed in 0.08 s · [line 130](../../tests/test_input_validation.py#L130) · **Area:** Input validation: feature and new-sequence dialogs · **Checked against:** hand-computed · **Strength:** 🟡 adequate

The feature dialog refuses empty, unreadable or out-of-range locations and malformed qualifiers with a message and accepts a valid one; the New Sequence dialog stays open on non-DNA text.

| # | Input | Expected output |
|---|---|---|
| 1 | 500 bp doc; location '  ' | 0 features added, message mentions 'location' |
| 2 | location 'banana' | 0 added, message has 'Cannot read the location' and 'complement' |
| 3 | location '5..99999' | 0 added, message has 'outside' |
| 4 | qualifiers '=orphan value'; then 'two words=x' | 0 added with 'key=value' message; 0 added |
| 5 | location ' 10..20 ', qualifiers '/note="fine"\n\n  gene = g1 ', name 'é ü' | 1 added; name maxLength >= 200; note qualifier == 'fine' |
| 6 | New Sequence dialog: text 'ACGT hello world', name 'n'*500, OK clicked | dialog not accepted, a 'Not DNA' label shown, name maxLength <= 100; new_file() returns None |

**Review notes:**

- Line 143: 'two words=x' is only checked for 0 added; no message is checked (docstring says bad qualifiers are refused 'with a message').
- Lines 144-147: the accepted feature's location (10..20), its 'gene=g1' qualifier and the non-ASCII name 'é ü' are not checked, only the note.

## `test_locus_box`

**Result:** ✅ passed in 0.07 s · [line 173](../../tests/test_input_validation.py#L173) · **Area:** Input validation: go-to / locus box · **Checked against:** hand-computed · **Strength:** 🟡 adequate

Typing bad positions in the Alignment view's locus box gives an error message instead of an exception; valid positions and ranges are accepted and selected.

| # | Input | Expected output |
|---|---|---|
| 1 | 500 bp doc; locus ',', '-', ':', '0', '99999999999999999999', '5000000', '–' then go() | editor.error emitted with a non-empty message for each |
| 2 | locus '12' and '  100 - 200 ' | no error; final selection (start, end) == (99, 200) |

**Review notes:**

- Line 186: only the last input's selection is checked; where '12' moved the cursor is not asserted.

## `test_small_limits_and_links`

**Result:** ✅ passed in 0.71 s · [line 189](../../tests/test_input_validation.py#L189) · **Area:** Input validation: limits and linked spin boxes · **Checked against:** hand-computed · **Strength:** 🟡 adequate

The search box has a length cap, linked min/max spin boxes stay ordered, tab names become one short line, the image export refuses a giant raster, and a primer-pair name is validated.

| # | Input | Expected output |
|---|---|---|
| 1 | search box maxLength; setText('>seq\nACGT' + 'N'*6000) and run(1) | maxLength == 5000; no exception |
| 2 | DesignDialog: min_len=30, then max_len=15 | max_len >= 30; then min_len <= 15 |
| 3 | AssemblyPanel: min_ov=50, then max_ov=20 | max_ov >= 50; then min_ov <= 20 |
| 4 | set_tab_name(0, 'a\nvery   long ' + 'x'*300) | no newline, length <= 60 |
| 5 | ImageExportDialog dpi 2400, scale 1000; then SVG format; then first format, dpi 300 | OK disabled and 'Too large' warning; OK enabled; OK enabled |
| 6 | DesignDialog save('bad\tname') and save('x'*100) with one pair row | None, 'control character' message; None, 'too long' message |

**Review notes:**

- Line 194-195: smoke only (no assertion on the result), and setText on a field with maxLength 5000 truncates the 6000+ char text, so an over-long query never actually reaches run().
