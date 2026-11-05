# `tests/test_bugfixes.py`

[← all test files](README.md) · source: [`tests/test_bugfixes.py`](../../tests/test_bugfixes.py)

## `test_tracks_enzymes_and_files`

**Result:** ✅ passed in 0.00 s · [line 18](../../tests/test_bugfixes.py#L18) · **Area:** Regressions: tracks, enzymes, FASTA, save · **Checked against:** hand-computed · **Strength:** 🟢 strong

Guards earlier bug fixes: track limit keeps the first items, empty GFF Name falls back, no cut past a linear end, FASTA with BOM or CRLF, failed saves leave no .tmp, fuzzy locations kept.

| # | Input | Expected output |
|---|---|---|
| 1 | 100 BED items every 10 bp; query(0,1000,10), (500,1000,3), limit 100, (503,520,5) | first ten f0000..f0090 more=True; f0500,f0510,f0520; more=False; f0500 first |
| 2 | GFF 'Name=;gene=lacZ' and 'ID=g2;Name=' | names ['lacZ','g2'] |
| 3 | BsaI on linear 'GGTCTCA', 'GGTCTCAAAAAAAAAA', circular 'GGTCTCA' | {}; [(7,11)]; a cut |
| 4 | FASTA with UTF-8 BOM; replace chr1 by 'GGGG' | sequences [('chr1',14),('chr2',12)]; fetch correct; output bytes start with BOM and keep chr2 |
| 5 | CRLF FASTA, replace chr2 with 70 G | every line CRLF; chr2 wrapped 60+10 |
| 6 | save_genbank / write_fasta onto a directory | OSError; no *.tmp left |

**Fixed by the audit of November 2026:**

- Removed in the second pass of the audit (fewer tests to review, no biological check lost): the fuzzy-location block was removed: test_core_functions.py::test_locations_in_genbank_files checks the same 10..(20.30) line, a join with fuzzy ends, and that both survive a save.

## `test_closing_a_pane_does_not_leave_the_editor_calling_a_deleted_view`

**Result:** ✅ passed in 0.09 s · [line 84](../../tests/test_bugfixes.py#L84) · **Area:** GUI: workspace panes · **Checked against:** smoke (only 'no crash'/'exists') · **Strength:** 🟡 adequate

Guards against a crash ('Internal C++ object already deleted') when editing after the Protein or Map pane was closed.

| # | Input | Expected output |
|---|---|---|
| 1 | MainWindow with a 200 bp circular doc; open ProteinView and MapView, close both panes, flush deferred deletes; type 'A' at 3 | no exception reaches sys.excepthook |
