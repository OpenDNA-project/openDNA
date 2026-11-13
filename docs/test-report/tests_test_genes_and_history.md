# `tests/test_genes_and_history.py`

[← all test files](README.md) · source: [`tests/test_genes_and_history.py`](../../tests/test_genes_and_history.py)

## `test_gene_groups_fold_their_parts_and_survive_edits`

**Result:** ✅ passed in 0.01 s · [line 59](../../tests/test_genes_and_history.py#L59) · **Area:** Gene grouping · **Checked against:** hand-computed · **Strength:** 🟢 strong

**Why it matters:** Without it, the feature list and map could fold the wrong features under a gene (or fail to fold them), and after editing the sequence a feature could end up covering bases it never covered.

Checks that features sharing a locus_tag or /gene fold under their gene while untagged or orphaned ones stay alone, which features each view is given under the cds_only and collapse_genes switches, and that after an insert before everything, a delete inside the gene, an insert inside the gene but outside the CDS, and an undo, every feature still covers exactly the bases it did.

| # | Input | Expected output |
|---|---|---|
| 1 | 9 features: gene/mRNA/CDS sharing L1, an untagged promoter, a solo misc_feature, a CDS with an unused tag L9, a gene+CDS sharing /gene abc, and an empty gene L2 | groups {L1, abc, L2}; L1 master 'gA' with parts ['mA','cA'] (ranked, then by position); abc master 'abc' with part 'abcCDS'; L2 has no parts |
| 2 | part_of mapping | exactly mA, cA and abcCDS are parts; part_of[cA] == 'L1'; pA, solo and orphan are neither parts nor masters |
| 3 | map_features / shown with cds_only on, then off, then collapse_genes off | {gA, pA, solo, orphan, abcCDS, g2} (the gene 'abc' hidden in favour of its single CDS); then {gA, pA, solo, orphan, abc, g2}; then every feature |
| 4 | insert of 40 A's at position 0, then delete of bases 300-310 | groups unchanged; cA loses its 110th-120th bases; gA loses its 160th-170th; mA 10 bp shorter; pA, solo, orphan, abcCDS, abc and g2 cover the same bases as before |
| 5 | genes.span_of on the L1 master and its parts | every part's span lies inside the gene's |
| 6 | insert of CCCCC inside the gene but outside the CDS | the CDS covers the same bases; the gene is 5 bp longer; master still 'gA' |
| 7 | undo_steps(3) | returns 3; every feature covers its original bases and the groups are {L1, abc, L2} again |

## `test_history_texts`

**Result:** ✅ passed in 0.01 s · [line 114](../../tests/test_genes_and_history.py#L114) · **Area:** Edit history wording · **Checked against:** regression value · **Strength:** 🟡 adequate

**Why it matters:** Without it, the History panel - the user's only record of what an edit actually did - could describe an edit wrongly or lump several edits into one line, making a mistake impossible to trace back.

Checks the exact text written for each kind of edit (insert, delete, short and long replace, long-sequence truncation with an ellipsis, paste with features, feature added/changed/renamed/deleted, many features, reverse complement), that an undone row keeps its text in the redo list, and the texts for circular operations and edits across the origin.

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

**Review notes:**

- Entirely string-literal assertions on UI wording: it pins the current text but a wrong-but-plausible message (e.g. an off-by-one position) is only caught where the position is spelled out.
- One assertion builds its expected text from the history entry it is checking, so that part cannot fail.
