# `tests/test_genes_and_history.py`

[← all test files](README.md) · source: [`tests/test_genes_and_history.py`](../../tests/test_genes_and_history.py)

## `test_cds_only_hides_a_gene_at_exactly_a_cds_and_edits_keep_the_bases`

**Result:** ✅ passed in 0.00 s · [line 56](../../tests/test_genes_and_history.py#L56) · **Area:** CDS only / editing · **Checked against:** hand-computed · **Strength:** 🟢 strong

**Why it matters:** Without it, a gene feature that duplicates its CDS could stay drawn (or a gene that is not a duplicate could vanish), and after editing the sequence a feature could end up covering bases it never covered.

Checks that "CDS only" hides a gene feature at exactly the bases and strand of a CDS and keeps the others, and that after an insert before everything, a delete inside the gene, an insert inside the gene but outside the CDS, and an undo, every feature still covers exactly the bases it did.

| # | Input | Expected output |
|---|---|---|
| 1 | 7 features: gene 101..400 + CDS 151..390, a promoter, gene 900..1000 + CDS 900..1000, complement gene 1100..1200 + CDS 1100..1200 | shown = {gA, cA, pA, abcCDS, g2, c2}: only gene abc (same bases and strand as its CDS) is hidden, also from visible_features |
| 2 | set_display(cds_only=False) | every feature is shown |
| 3 | insert of 40 A's at position 0, then delete of bases 300-310 | cA loses its 110th-120th bases; gA loses its 160th-170th; pA, abcCDS, abc, g2 and c2 cover the same bases as before |
| 4 | insert of CCCCC inside the gene but outside the CDS | the CDS covers the same bases; the gene is 5 bp longer |
| 5 | undo of the three edits | every feature covers the bases it covered at the start |

## `test_history_texts`

**Result:** ✅ passed in 0.01 s · [line 90](../../tests/test_genes_and_history.py#L90) · **Area:** Edit history wording · **Checked against:** regression value · **Strength:** 🟡 adequate

**Why it matters:** Without it, the History panel - the user's only record of what an edit actually did - could describe an edit wrongly or lump several edits into one line, making a mistake impossible to trace back.

Checks the exact text written for each kind of edit (insert, delete, short and long replace, long-sequence truncation with an ellipsis, paste with features, feature added/changed/renamed/deleted, many features, reverse complement), that an undone row keeps its text in the redo list, and the texts for circular operations and edits across the origin. Also that Set Origin to Selection Start rotates a circle (and refuses at base 1 or on a line).

| # | Input | Expected output |
|---|---|---|
| 1 | insert GGG at cursor 10 / delete that selection | 'Inserted 3 bp at 11: GGG'; 'Deleted 3 bp at 11-13: GGG' |
| 2 | replace 3 bases at 20 / replace 30 bases with 5 at 30 | '<old> > <new> at 21' (both sides spelled out); 'Replaced 30 bp with 5 bp at 31' (counts only above 12 bp) |
| 3 | insert 40 C's at 0 / delete 40 bases at 100 | text shows the first 30 bases then an ellipsis, in both directions |
| 4 | paste 6 bp with one feature / with two features | 'Pasted 6 bp at 201: GGGCCC; Added feature pf'; 'Pasted 6 bp at 251: TTTAAA; Added features q1, q2' |
| 5 | add a feature, change its location, rename it, delete it | 'Added feature mine'; 'Changed feature mine'; 'Added feature renamed; Deleted feature mine'; 'Deleted feature renamed' |
| 6 | add five features at once / reverse complement | 'Added features n0, n1, n2, n3...' (at most four names); 'Reverse-complemented the whole sequence' |
| 7 | undo the reverse complement | the redo list holds the same text |
| 8 | circular sequence: rotate(7), linearize(0), set_circular, linearize(50), set_circular, then replace 15 bases across the origin | 'Moved base 1 (broke the circle elsewhere)', 'Made the sequence linear', 'Made the sequence circular', 'Made the sequence linear (broke the circle elsewhere)', 'Made the sequence circular', 'Replaced 15 bp across the origin from 191 with GG' |
| 9 | delete 15 bases across the origin | text starts 'Replaced 15 bp across the origin from ' and ends ' (deleted)' |
| 10 | set_origin(10) on the circle; set_origin(0); linear sequence set_origin(5) | sequence becomes seq[10:]+seq[:10], still circular, cursor at 0 (True); then False; then False |

**Review notes:**

- Entirely string-literal assertions on UI wording: it pins the current text but a wrong-but-plausible message (e.g. an off-by-one position) is only caught where the position is spelled out.
- One assertion builds its expected text from the history entry it is checking, so that part cannot fail.

**Fixed by the audit of November 2026:**

- 13 Nov 2026: gap closed. The proposed checks were added, pass on the current code and fail under the mutation that slipped through before. Was: 'the text shows the first 30 of the 40 deleted bases' is not tested: the expected text at the long-delete step is 'Deleted 40 bp at 101-140: ' + ed.doc.history()[0][-1][1] + ellipsis, so the stored bases are compared with themselves. If the undo record kept the last 30 removed bases, bases from the wrong offset, or the inserted bases instead, the row would be wrong on screen and the test would still pass. (The insert side, 'C'*30, is a literal and is pinned.) \| The last circular case is matched only with startswith and endswith, so the base it names ('across the origin from <base>') is never checked - an off-by-one in the position of a delete across the origin slips through, although the same position is spelled out and checked for the replace one line above.
