# `tests/test_input_parsing.py`

[← all test files](README.md) · source: [`tests/test_input_parsing.py`](../../tests/test_input_parsing.py)

## `test_typed_or_imported_text_is_parsed_or_refused_with_a_reason`

**Result:** ✅ passed in 0.01 s · [line 22](../../tests/test_input_parsing.py#L22) · **Area:** Input checking: typed and imported text · **Checked against:** hand-computed · **Strength:** 🟢 strong

**Why it matters:** Without it the range box could select the wrong bases, pasted protein text could keep junk characters, and bad names, NCBI keys or primer lists could end up stored in the primer database.

Checks the range box parser on common notations and refusals, cleaning of pasted protein text, the name and API-key validators, single-primer problems, a deliberately bad TSV primer list and a FASTA list (each refused line numbered with its reason), a binary file, and that the primer database itself refuses invalid add / add_many / update.

| # | Input | Expected output |
|---|---|---|
| 1 | selection_range for 100..200, a cursor, 3699..3830 on 3780 bp, and the whole sequence, typed back through parse_range as 'first..last' | (101, 200); (None, None); (3700, 50); each round-trips to the same (start, end) |
| 2 | grouped('') / '7' / '1234567' / '0042', and '1,234..2,500' (the boxes' own text) through parse_range on 3780 bp | '' / '7' / '1,234,567' / '42'; (1233, 2500) |
| 3 | parse_range '101..200' on a 3780 bp circle | (100, 200) |
| 4 | ' 101-200 ', '101 to 200', '200..101' (linear) | All equal to the same range |
| 5 | '1,000..2,000'; '3700..50' circular / linear; '42' | (999, 2000); (3699, 3830) across the origin / (49, 3700); (42, 42) cursor |
| 6 | '', 'abc', '1..2..3', '0..10', '100..4000' | ValueError each |
| 7 | clean_pasted_protein('>p\nMKT AYI\n10 AKQ') / ('MKT123?!') | 'MKTAYIAKQ' / ValueError |
| 8 | text_problem 'fine name' / 'a\tb' / 30 chars (max 20); api_key_problem '0123abcdef' / 'abc-def!' | '' / 'control character' / 'too long'; '' / 'letters and digits' |
| 9 | primer_problem: empty name, tab in name, 6 bases, 201 bases, 'ACGTACGTXX' | Mentions 'name', 'control character', 'at least 8', 'at most', 'not DNA' |
| 10 | TSV list: good1, words, short, one-column line, no name, uracil, blank + comment, comma-separated good2 | good1 and good2 imported; rejected lines [3,4,5,6,7] with reasons not DNA / at least 8 / No tab, comma or semicolon / name / U; report says '5 lines not imported' and 'line 3'; rejected names absent from db |
| 11 | FASTA list: fa1 over two lines, fa2 'ACGT-XYZ-ACGT', an unnamed record | [('fa1','ACGTACGTACGTAC')]; rejected lines [4, 6], first 'not DNA' |
| 12 | TSV containing a NUL byte | ValueError 'not a text file' |
| 13 | db().add with tab in name / 'ACGTXX'; add_many with one good, one tab name, one 500-base; update to 'NOT DNA AT ALL' | 'invalid' twice; added 1, invalid ['bad\tone','toolong']; ValueError |

**Review notes:**

- Bundles five unrelated input areas (range box, protein paste, name/key validation, primer lists, primer database) in one test; a failure in the range box hides all the primer checks. Splitting into 2-3 tests would help.
- Writes good1, good2 and fine_one into the primer database of the module's temp data dir; harmless, but state is shared with anything else in the same process using that dir.
