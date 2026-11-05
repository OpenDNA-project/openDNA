# `core/src/document.rs`

[← all test files](README.md) · source: [`core/src/document.rs`](../../core/src/document.rs)

## `edit_moves_features_and_undo_restores_exactly`

**Result:** ✅ passed · [line 1205](../../core/src/document.rs#L1205) · **Area:** Document editing / undo · **Checked against:** hand-computed · **Strength:** 🟢 strong

An insertion inside a feature grows it, marks the document modified, and undo/redo restore the exact record.

| # | Input | Expected output |
|---|---|---|
| 1 | sample .dna (20 bp, lacZ 2..10); replace(4,4,'TTT') | lacZ '2..13', len 23, modified |
| 2 | undo; redo | record == original, not modified; then lacZ '2..13' again |

## `rejects_bad_edits`

**Result:** ✅ passed · [line 1220](../../core/src/document.rs#L1220) · **Area:** Document editing / validation · **Checked against:** hand-computed · **Strength:** 🟡 adequate

Inverted ranges, ranges past the end and invalid letters are refused and leave no history.

| # | Input | Expected output |
|---|---|---|
| 1 | replace(5,4,''); replace(0,999,''); replace(0,0,'AXG') | all errors |
| 2 | undo | false (no history) |

**Review notes:**

- Record content after the refused edits is not compared (only the undo stack).

## `wrapped_qualifiers_are_stable_across_saves`

**Result:** ✅ passed · [line 1229](../../core/src/document.rs#L1229) · **Area:** GenBank I/O · **Checked against:** round-trip/self-consistency · **Strength:** 🟢 strong

Qualifiers wrapped over several lines are joined correctly on read (spaces in notes, none in translations) and do not drift over repeated saves.

| # | Input | Expected output |
|---|---|---|
| 1 | GenBank text, 70 bp, CDS 1..9 with a /note wrapped over 3 lines and a /translation wrapped over 2 lines | note 'a long note that the original writer wrapped over two lines.'; translation 'MKVLAAGIVALLLAAGCSSSKEETPAPAAAPAPAAEAPAAPAAEAPAAPAAE' |
| 2 | save, reopen, save, reopen | features equal after each reopen |

**Review notes:**

- Fixed file name ge_test_wrapped.gb in the shared temp dir (line 1280), never removed.

## `quotes_survive_at_every_wrap_position`

**Result:** ✅ passed · [line 1262](../../core/src/document.rs#L1262) · **Area:** GenBank I/O · **Checked against:** round-trip/self-consistency · **Strength:** 🟢 strong

Regression for the patched gb-io writer: a note containing a quoted word survives save/reopen whatever its offset relative to the wrap column.

| # | Input | Expected output |
|---|---|---|
| 1 | for k in 0..120: note 'w' + k*'x' + ' "quoted" end' added to feature 0; save, reopen | reopened feature 0 has a note equal to the original (message names offset k) |

## `feature_edits_are_undoable`

**Result:** ✅ passed · [line 1278](../../core/src/document.rs#L1278) · **Area:** Feature editing / undo · **Checked against:** hand-computed · **Strength:** 🟢 strong

Adding and deleting features are undo steps; multi-line qualifier text is normalised; invalid feature edits are refused without history.

| # | Input | Expected output |
|---|---|---|
| 1 | set_feature(None, misc_feature complement(3..8) note 'a\n  b') | index = old count; note 'a b' |
| 2 | delete_feature(0) | count back to original |
| 3 | set_feature with location 0..999; set_feature(Some(99), ..) | both errors |
| 4 | undo, undo, undo; redo | two undos succeed, third fails; features == original; after redo count = original+1 |

## `several_feature_changes_are_one_undo_step`

**Result:** ✅ passed · [line 1300](../../core/src/document.rs#L1300) · **Area:** Feature editing / undo · **Checked against:** hand-computed · **Strength:** 🟢 strong

edit_features applies deletes, in-place replacements and additions atomically as one undo step; an invalid batch applies nothing.

| # | Input | Expected output |
|---|---|---|
| 1 | [(0, delete), (1, replace with kind CDS), (None, add gene 1..10)] | features == [renamed #1] + original[2..] + [gene] |
| 2 | undo | features == original |
| 3 | [(0, delete), (99, delete)] | error; features unchanged |

## `refreshing_translations_is_one_undo_step_and_only_touches_stale_cds`

**Result:** ✅ passed · [line 1318](../../core/src/document.rs#L1318) · **Area:** CDS /translation refresh · **Checked against:** hand-computed · **Strength:** 🟢 strong

Refresh translations updates only CDSs that already carry a stale /translation, skips pseudogenes, is one undo step, and is a no-op the second time.

| # | Input | Expected output |
|---|---|---|
| 1 | 'ATGAAATAA' x3; CDS 1..9 /translation='MQ', CDS 10..18 none, CDS 19..27 /pseudo /translation='MQ'; refresh(None, table 1) | report [(0, Updated), (2, Skipped(Pseudo))]; translations 'MK', none, 'MQ'; history +1 |
| 2 | refresh again | report[0] = (0, Unchanged); history unchanged |
| 3 | undo | features == before |

## `copy_paste_carries_features`

**Result:** ✅ passed · [line 1351](../../core/src/document.rs#L1351) · **Area:** Copy/paste with features · **Checked against:** hand-computed · **Strength:** 🟡 adequate

Copying a range brings the overlapping part of features with relative coordinates; pasting adds them; undo removes them; features that do not fit the pasted bases are refused.

| # | Input | Expected output |
|---|---|---|
| 1 | sample; features_between(4, 12) | lacZ copy at '1..6' |
| 2 | paste(0,0, seq[4..12], copied) | len 28, count n + copied, feature n at '1..6' |
| 3 | undo; paste(0,0,'AC', copied) | count n; error |

**Review notes:**

- The shift of the existing features by 8 after pasting at 0 is not asserted.

## `reverse_complement_twice_is_identity`

**Result:** ✅ passed · [line 1369](../../core/src/document.rs#L1369) · **Area:** Reverse complement · **Checked against:** hand-computed · **Strength:** 🟢 strong

Reverse-complementing flips bases and feature coordinates; doing it twice, or undoing both, restores the record.

| # | Input | Expected output |
|---|---|---|
| 1 | sample 'AAAACCCCGGGGTTTTACGT'; reverse_complement | seq 'ACGTAAAACCCCGGGGTTTT'; lacZ 2..10 -> 'complement(11..19)' |
| 2 | reverse_complement again; then undo x2 | record == original both times |

## `pcr_product_with_tails_and_across_the_origin`

**Result:** ✅ passed · [line 1383](../../core/src/document.rs#L1383) · **Area:** PCR product · **Checked against:** hand-computed · **Strength:** 🟢 strong

A PCR product is a new linear, unsaved document with primer tails added and features shifted; it can span the origin; invalid ranges and tails are refused; overwrite keeps features.

| # | Input | Expected output |
|---|---|---|
| 1 | pcr_product(1, 10, 'GG', 'CC') | seq 'GGAAACCCCGGCC'; lacZ '3..11'; modified, linear |
| 2 | pcr_product(16, 24) on 20 bp circle | seq 'ACGTAAAA' |
| 3 | pcr_product(5,5); pcr_product(0,4,'XX','') | errors |
| 4 | product(1,10) then overwrite(0,'T') | seq[0]='T'; feature 0 location unchanged |

## `saving_replaces_primer_sites_without_touching_the_document`

**Result:** ✅ passed · [line 1404](../../core/src/document.rs#L1404) · **Area:** GenBank I/O (primers) · **Checked against:** round-trip/self-consistency · **Strength:** 🟢 strong

Saving with a primer list writes those primer_bind features instead of the document's, without changing the in-memory record; an out-of-range site is refused.

| # | Input | Expected output |
|---|---|---|
| 1 | sample (1 primer_bind p1); save with [primer_bind 2..12 'new'] | record unchanged, not modified; reopened file has 1 primer_bind labelled 'new', same total feature count |
| 2 | save with site 0..999 | error |

**Review notes:**

- Comment says 'nothing written' for the bad site but the file content after the failed save is not checked.

## `gc_bins_exact_and_from_blocks`

**Result:** ✅ passed · [line 1427](../../core/src/document.rs#L1427) · **Area:** GC content plot · **Checked against:** hand-computed · **Strength:** 🟡 adequate

GC fractions per bin are exact for short ranges and correct through the cached-block path for long sequences, including cache invalidation after edits and undo.

| # | Input | Expected output |
|---|---|---|
| 1 | 'AAAACCCCGGGGTTTTACGT', 0..20 in 2 bins; 0..8 in 1 bin | [0.6, 0.4]; [0.5] |
| 2 | replace all with 'GC' x 20000; 0..40000 in 2 bins | [1.0, 1.0] |
| 3 | replace all with 'AT' x 20000; then undo | [0.0, 0.0]; then [1.0] |

**Review notes:**

- The block path (lines 1468-1474) is only fed uniform sequences (all GC or all AT), so a wrong block index, off-by-one block edge or bins swapped would still give 1.0/0.0; a mixed sequence is needed to test it.

## `genbank_round_trip`

**Result:** ✅ passed · [line 1441](../../core/src/document.rs#L1441) · **Area:** GenBank I/O · **Checked against:** round-trip/self-consistency · **Strength:** 🟡 adequate

An edited .dna-derived document saved as GenBank reopens with identical bases and features, and saving clears modified.

| # | Input | Expected output |
|---|---|---|
| 1 | sample; replace(0,2,'GG'); save_genbank; open | not modified; seq and features equal |

**Review notes:**

- Topology and definition are not compared after the round trip.

## `copy_and_edit_across_the_origin`

**Result:** ✅ passed · [line 1453](../../core/src/document.rs#L1453) · **Area:** Circular editing across the origin · **Checked against:** hand-computed · **Strength:** 🟢 strong

Copying across the origin brings an origin-spanning feature as one contiguous range; replacing across the origin is one exact undo step; a range not crossing the origin is refused.

| # | Input | Expected output |
|---|---|---|
| 1 | sample (ori complement(join(18..20,1..3))); features_across(15, 25) | one reverse feature, spans [(2,8)], strand -1 |
| 2 | edit_across_origin(17, 23, 'GG') | returns 14, len 16, seq[14..16]='GG' |
| 3 | undo; redo; edit_across_origin(5,10,'') | record == original, no further undo, redo gives len 16; error |

## `rotate_moves_bases_and_features`

**Result:** ✅ passed · [line 1474](../../core/src/document.rs#L1474) · **Area:** Rotate (set origin) · **Checked against:** hand-computed · **Strength:** 🟢 strong

Rotating a circular sequence moves bases and features; an origin-spanning feature becomes one piece; undo is exact; rotate(0) is refused.

| # | Input | Expected output |
|---|---|---|
| 1 | sample; rotate(5) | seq = old[5..] + old[..5]; ori spans [(12,18)] |
| 2 | undo; rotate(0) | record == original; error |

## `history_says_what_each_step_did`

**Result:** ✅ passed · [line 1489](../../core/src/document.rs#L1489) · **Area:** History panel · **Checked against:** hand-computed · **Strength:** 🟢 strong

The history summary records position, removed length, inserted text and feature additions/deletions for each step, and moves entries to redo on undo.

| # | Input | Expected output |
|---|---|---|
| 1 | replace(4,7,'GAG'); replace(0,0,'TT'); delete_feature(0) | 3 undo entries: (pos 4, removed 3, 'GAG'), (0, 'TT', 0), 1 deleted / 0 added; redo empty |
| 2 | undo | 2 undo, 1 redo, redo[0] has 1 deleted feature |

## `empty_edits_change_nothing`

**Result:** ✅ passed · [line 1602](../../core/src/document.rs#L1602) · **Area:** Document editing / undo · **Checked against:** hand-computed · **Strength:** 🟢 strong

Empty replace/overwrite calls at start, middle and end change nothing and create no undo step.

| # | Input | Expected output |
|---|---|---|
| 1 | 50 random bp (Lcg 1), circular, gene 11..20; replace(i,i,'') and overwrite(i,'') for i in 0, 7, 50 | not modified; undo() false; record unchanged |

## `random_edits_undo_redo_follow_a_snapshot_and_string_model`

**Result:** ✅ passed · [line 1618](../../core/src/document.rs#L1618) · **Area:** Document editing / undo (model-based) · **Checked against:** in-test reference model · **Strength:** 🟢 strong

600 random operations per topology (replace, overwrite, paste, rotate, reverse complement, add feature, edit across origin, undo, redo, save) keep the bases equal to a plain string model and the record equal to recorded snapshots, with correct undo/redo availability and modified flag.

| # | Input | Expected output |
|---|---|---|
| 1 | 120 random bp (Lcg 11 circular / 12 linear) with gene 10..40, CDS complement(join(50..60,70..80)), misc <90..>110 (+ origin-spanning join(115..120,1..5) if circular); 600 random ops | after every step: seq == string model; record == snapshot[pos]; record.len == seq.len; modified == (saved != pos); undo/redo availability as modelled; no redo after a new edit; edit_across_origin returns start-(end-n) |
| 2 | end of run | more than 100 edits applied; redo all / undo all walk the snapshots exactly; original seq restored; redo all again matches the last snapshot |

**Review notes:**

- Feature locations are compared only with snapshots taken from the document itself (self-consistency); no independent model checks where features land after each edit (other drift tests cover that).

## `circular_edit_across_the_origin_equals_the_rotated_string_edit`

**Result:** ✅ passed · [line 1766](../../core/src/document.rs#L1766) · **Area:** Circular editing across the origin · **Checked against:** in-test reference model · **Strength:** 🟡 adequate

Replacing a range that wraps the origin gives a sequence that is a rotation of the equivalent string edit, returns where the new bases start, keeps the circle, and is one exact undo step.

| # | Input | Expected output |
|---|---|---|
| 1 | 60 random bp (Lcg 5), circular, genes 11..30 and 36..50; (start,end,new) in (55,63,'GGG'), (50,61,''), (59,61,'TTTTT'), (30,90,''), (30,90,'A'), (59,61,''), (10,70,'CC') | length == model length; result is a rotation of new + rot[end-start..]; got[at..at+len(new)] == new; circular; 1 undo step; undo restores record |

**Review notes:**

- Only 'is some rotation' of the model is checked (line 1832), not the exact rotation, and the two genes' positions after the edit are never asserted.

## `circular_inserts_at_zero_and_at_the_end_and_wrapping_deletes`

**Result:** ✅ passed · [line 1794](../../core/src/document.rs#L1794) · **Area:** Circular editing / feature drift · **Checked against:** hand-computed · **Strength:** 🟢 strong

Inserting at position 0 or n of a circle, and deleting across the origin, keep features (including one spanning the origin) on their own bases.

| # | Input | Expected output |
|---|---|---|
| 1 | 60 random bp (Lcg 6), circular, gene 11..20, rep_origin join(56..60,1..5); replace(at,at,'TTT') for at in 0, 60 | len 63; circular; 'o' covers seq[55..]+seq[..5]; 'g' covers seq[10..20]; len field consistent; undo exact |
| 2 | edit_across_origin(57, 63, '') | len 54; seq == seq[3..57]; 'g' unchanged; 'o' covers seq[55..57]+seq[3..5]; undo exact |

## `paste_at_every_boundary_matches_the_string_model_and_is_one_step`

**Result:** ✅ passed · [line 1824](../../core/src/document.rs#L1824) · **Area:** Copy/paste with features · **Checked against:** in-test reference model · **Strength:** 🟢 strong

Pasting at the start, end, over a selection, the whole sequence and inside a feature produces the string-model bases, a pasted feature exactly on the pasted bases, and correct behaviour of the existing feature.

| # | Input | Expected output |
|---|---|---|
| 1 | 100 random bp (Lcg 7), gene 11..40; paste seq[20..35] with its feature at (s,e) in (0,0),(0,10),(100,100),(90,100),(0,100),(40,50),(40,40),(45,45),(12,14) | seq == seq[..s]+pasted+seq[e..]; last feature spans [(s, s+15)] and covers the pasted bases; 1 undo step; undo exact |
| 2 | existing gene g per case | (0,100): only 1 feature left; (12,14): g covers 30-2+15 bases; others: g still covers seq[10..40] |

**Review notes:**

- For (12,14) only the length of g is checked (line 1897), not its bases.

## `overwrite_keeps_length_moves_nothing_and_is_one_step_per_call`

**Result:** ✅ passed · [line 1856](../../core/src/document.rs#L1856) · **Area:** Overwrite mode · **Checked against:** hand-computed · **Strength:** 🟡 adequate

Overwrite keeps length and never moves or trims features, each call is one undo step, and overwriting past the end is refused.

| # | Input | Expected output |
|---|---|---|
| 1 | 80 random bp (Lcg 8), circular, gene 11..40, complement(join(50..60,70..80)); overwrite(10,'TTTT'), (35,'CCCCCC'), (76,'GGGG') | len 80; features unchanged; 3 history steps |
| 2 | overwrite(78,'ACG') | error, still 3 steps |
| 3 | undo x3 | 3 succeed, 4th fails; record == before |

**Review notes:**

- The overwritten bases themselves are never checked (e.g. seq[10..14]=='TTTT').

## `letters_are_kept_as_typed_and_invalid_ones_are_refused_not_dropped`

**Result:** ✅ passed · [line 1873](../../core/src/document.rs#L1873) · **Area:** Document editing / validation · **Checked against:** hand-computed · **Strength:** 🟢 strong

Spaces, digits, '*', non-ASCII, newlines and '.' are refused by replace, overwrite, paste and new_sequence; lower case, U, IUPAC codes, N and '-' are kept exactly as typed.

| # | Input | Expected output |
|---|---|---|
| 1 | 'A C', 'A1', 'A*', UTF-8 'e-acute', 'A\n', 'A.G' via replace/overwrite/paste/new_sequence | every call errors; record unchanged; no undo step |
| 2 | replace(4,4,'acguRYSWKMBDHVNn-') | len 8+17; inserted bytes identical; undo exact |

## `editing_an_empty_sequence`

**Result:** ✅ passed · [line 1895](../../core/src/document.rs#L1895) · **Area:** Document editing / edge cases · **Checked against:** hand-computed · **Strength:** 🟡 adequate

An empty sequence (linear and circular) accepts inserts, refuses out-of-range edits and origin operations, and does not panic in reverse complement, GC plot or copy.

| # | Input | Expected output |
|---|---|---|
| 1 | new_sequence('e', '', circular in false/true); replace(0,0,'') | len 0; no undo step |
| 2 | replace(0,1,''); replace(1,1,'A') | errors |
| 3 | replace(0,0,'ACG'); undo; redo; replace(0,3,'') | 'ACG'; len 0 and len field 0; 'ACG'; len 0 |
| 4 | rotate(1); edit_across_origin(0,1,''); reverse_complement; gc_bins(0,0,3); features_between(0,0) | errors; no panic; 3 bins; empty |

**Review notes:**

- gc_bins on an empty range only checks there are 3 values (expected NaN), not their value.

## `modified_flag_follows_edits_undo_and_save`

**Result:** ✅ passed · [line 1918](../../core/src/document.rs#L1918) · **Area:** Modified flag · **Checked against:** hand-computed · **Strength:** 🟢 strong

The modified flag tracks the saved state through edits, undo, redo and save, including when the saved state becomes unreachable; a new sequence starts modified.

| # | Input | Expected output |
|---|---|---|
| 1 | sample; edit; undo; redo | modified / not / modified |
| 2 | save; edit; undo; undo; redo | not modified / modified / not / modified / not |
| 3 | undo; new edit; undo | modified both times (saved state gone) |
| 4 | new_sequence('n','ACGT') | modified |

## `external_steps_sit_in_the_history_in_order_and_never_change_the_modified_flag`

**Result:** ✅ passed · [line 1944](../../core/src/document.rs#L1944) · **Area:** History (external primer steps) · **Checked against:** hand-computed · **Strength:** 🟢 strong

App-level 'external' steps (primer database changes) sit in the undo history in order, are announced before undo/redo, never change the modified flag, and are dropped with the redo stack.

| # | Input | Expected output |
|---|---|---|
| 1 | external('Added primer p1','t1'); replace(0,1,'G'); external('Deleted primer p2','t2') | not modified after first; then modified; 3 entries, entry 0 = ('external','Added primer p1','t1'), entry 1 note '' |
| 2 | next_external / undo sequence | 't2' then undo (modified), None (base edit next), redo-side 't2', undo (not modified), 't1', undo (not modified, seq == before), no more undo |
| 3 | redo x3; undo x2; new edit | modified, no 4th redo; next_external(true) None |
| 4 | save; external('Edited primer p3','t3'); undo | not modified throughout |

## `insert_at_the_exact_start_or_end_of_a_feature_does_not_grow_it`

**Result:** ✅ passed · [line 1990](../../core/src/document.rs#L1990) · **Area:** Feature drift on insert · **Checked against:** hand-computed · **Strength:** 🟢 strong

Inserting exactly at a feature's first base shifts it; inserting right after its last base leaves it alone; in both cases it covers the same bases.

| # | Input | Expected output |
|---|---|---|
| 1 | 200 random bp (Lcg 99), gene 51..100; replace(50,50,'TTT') | spans [(53,103)], covers s[50..100] |
| 2 | undo; replace(100,100,'TTT') | spans [(50,100)], covers s[50..100]; undo restores record |

## `delete_overlapping_the_start_the_end_or_the_middle_keeps_the_remaining_bases`

**Result:** ✅ passed · [line 2006](../../core/src/document.rs#L2006) · **Area:** Feature drift on delete · **Checked against:** in-test reference model · **Strength:** 🟢 strong

A deletion overlapping a feature's start, end, middle or edges leaves the feature on exactly its surviving bases.

| # | Input | Expected output |
|---|---|---|
| 1 | gene 51..100 on 200 random bp; delete (40,60), (90,110), (70,80), (50,60), (90,100), (20,70) | covered == s[i] for i in 50..100 not in deleted range; one span; undo exact |

## `deleting_exactly_the_feature_or_more_removes_it_and_undo_restores_it_whole`

**Result:** ✅ passed · [line 2021](../../core/src/document.rs#L2021) · **Area:** Feature drift on delete · **Checked against:** hand-computed · **Strength:** 🟡 adequate

Deleting all bases of a (complement join) feature removes it; a partial delete keeps it; undo restores location and qualifiers; redo/undo is stable.

| # | Input | Expected output |
|---|---|---|
| 1 | CDS complement(join(51..70,81..100)) with note and codon_start; delete (50,100), (49,101), (0,200) | no features left |
| 2 | delete (30,100), (50,130) | also no features left (both ranges cover 0-based 50..100, so the 'keep what is left' branch is not reached) |
| 3 | undo; redo; undo | record == before; features == before |

**Review notes:**

- Cases (30,100) and (50,130) delete the whole feature (0-based 50..100) too, yet line 2083 condition a<=50 && b>=100 sends them to the 'removed' branch; the else branch 'must keep what is left' is never executed by any case, so partial deletion is not tested here.

## `replace_over_part_of_a_feature_follows_the_documented_rules`

**Result:** ✅ passed · [line 2042](../../core/src/document.rs#L2042) · **Area:** Feature drift on replace · **Checked against:** hand-computed · **Strength:** 🟢 strong

Replace = delete then insert at start: new bases join the feature only when the replaced stretch is strictly inside; at the start/end they stay outside; replacing all of it removes it.

| # | Input | Expected output |
|---|---|---|
| 1 | gene 51..100; replace(60,70, 13 T) | covers s[50..60]+13T+s[70..100] |
| 2 | replace(60,70, 10 T); replace(60,70,'GG') | spans [(50,100)] and same-length cover; shrunk cover with GG |
| 3 | replace(40,60, 13 T); replace(90,110, 13 T) | covers s[60..100]; covers s[50..90] |
| 4 | replace(50,100, 13 T) | no features |

## `complement_feature_keeps_its_bases_and_strand_through_edits`

**Result:** ✅ passed · [line 2075](../../core/src/document.rs#L2075) · **Area:** Feature drift (complement) · **Checked against:** hand-computed · **Strength:** 🟢 strong

A bottom-strand feature keeps its strand and bases through edits before/after it and loses exactly the deleted bases inside.

| # | Input | Expected output |
|---|---|---|
| 1 | gene complement(51..100); insert 'ACGTA' at 10, 'CC' at 190 | strand -1; covers s[50..100] |
| 2 | delete 70..75 | covers own[..15]+own[20..]; bio == rc(covered); location 'complement(56..100)' |

**Review notes:**

- bio == rc(covered) (line 2137) is true by construction of the helper for any strand -1 feature, so it adds nothing.

## `join_edit_between_parts_grows_the_gap_inside_a_part_grows_that_part`

**Result:** ✅ passed · [line 2089](../../core/src/document.rs#L2089) · **Area:** Feature drift (joins) · **Checked against:** hand-computed · **Strength:** 🟢 strong

For a two-part join, insertions in the gap or at part boundaries grow the gap, insertions inside a part grow that part, and a delete across the gap trims both parts.

| # | Input | Expected output |
|---|---|---|
| 1 | join(21..40,81..100); insert 'ACGT' at 60, 30, 90, 40, 80, 20, 100, 0 | spans resp. [(20,40),(84,104)], [(20,44),(84,104)], [(20,40),(80,104)], [(20,40),(84,104)], [(20,40),(84,104)], [(24,44),(84,104)], [(20,40),(80,100)], [(24,44),(84,104)]; covered length matches; bases unchanged when size unchanged |
| 2 | delete 30..90 | covers s[20..30]+s[90..100] |

## `complement_join_order_survives_edits_and_reverse_complement`

**Result:** ✅ passed · [line 2120](../../core/src/document.rs#L2120) · **Area:** Feature drift (joins) / reverse complement · **Checked against:** in-test reference model · **Strength:** 🟢 strong

A join and a complement(join) keep their own sequence and part order through outside edits, and reverse-complementing mirrors their parts and flips the strand.

| # | Input | Expected output |
|---|---|---|
| 1 | complement(join(21..40,81..100)) and join(21..40,81..100); check initial bio | rc(s[20..40]+s[80..100]) / s[20..40]+s[80..100] |
| 2 | insert 'GGGG' at 60, delete 0..5 | bio unchanged; 2 parts in ascending order |
| 3 | reverse_complement; reverse_complement | strand flipped, bio unchanged, sorted parts == mirrored (n-b, n-a); twice: parts == old parts |

## `feature_crossing_the_origin_keeps_its_bases`

**Result:** ✅ passed · [line 2150](../../core/src/document.rs#L2150) · **Area:** Feature drift (origin-spanning) · **Checked against:** hand-computed · **Strength:** 🟢 strong

A feature spanning the origin of a circle keeps its bases through edits outside it, grows the right part on inserts inside, trims on deletes, and keeps its bases in order under any rotation.

| # | Input | Expected output |
|---|---|---|
| 1 | 200 random bp (Lcg 98), circular, o = join(191..200,1..10), g = 101..150; insert 'ACGT' / replace at (0,0),(200,200),(100,110),(160,170),(10,10),(190,190) | o covers s[190..]+s[..10] |
| 2 | insert at 195; insert at 5 | o covers s[190..195]+ACGT+s[195..]+s[..10]; s[190..]+s[..5]+ACGT+s[5..10] |
| 3 | delete 0..5; delete 195..200 | o covers s[190..]+s[5..10]; s[190..195]+s[..10] |
| 4 | rotate(k) for k in 1,5,9,10,11,60,100,150,190,191,195,199 | bio(o) == own; bio(g) == s[100..150]; undo exact |

## `edit_across_the_origin_keeps_features_on_both_sides`

**Result:** ✅ passed · [line 2192](../../core/src/document.rs#L2192) · **Area:** Circular editing across the origin · **Checked against:** hand-computed · **Strength:** 🟢 strong

Replacing a stretch across the origin trims only the overlapped bases of features on either side and places pasted features on the new bases.

| # | Input | Expected output |
|---|---|---|
| 1 | 200 random bp (Lcg 97), circular, a 5..30, b 100..150, c 170..195; edit_across_origin(190, 205, 'TTT') | a covers s[5..30]; b covers s[99..150]; c covers s[169..190]; undo exact |
| 2 | same with pasted misc_feature 1..3 'p' | p covers 'TTT' at [(at, at+3)]; b unchanged |

## `fuzzy_locations_survive_edits_and_saves_with_their_markers`

**Result:** ✅ passed · [line 2218](../../core/src/document.rs#L2218) · **Area:** Feature drift (partial < > markers) / GenBank I/O · **Checked against:** hand-computed · **Strength:** 🟢 strong

Partial-end markers < and > are kept and shifted correctly by edits, survive a GenBank save, and come back on undo.

| # | Input | Expected output |
|---|---|---|
| 1 | misc <21..>60, CDS complement(<101..>150) on 200 bp; insert 3 at 0 | '<24..>63' |
| 2 | insert 2 at 40 | '<24..>65' |
| 3 | delete 150..160 | 'complement(<106..>150)' |
| 4 | save, reopen; undo x3 | features equal; '<21..>60' |

**Review notes:**

- Fixed temp file name ge_test_fuzzy.gb (line 2279).

## `reverse_complement_flips_every_feature_onto_the_same_bases`

**Result:** ✅ passed · [line 2236](../../core/src/document.rs#L2236) · **Area:** Reverse complement · **Checked against:** in-test reference model · **Strength:** 🟢 strong

Reverse complement flips the strand of plain, complement, join and complement(join) features while they keep their own sequence, on linear and circular records.

| # | Input | Expected output |
|---|---|---|
| 1 | 200 random bp (Lcg 96), features 21..60, complement(81..120), join(131..140,151..160), complement(join(171..180,191..200)); linear and circular; reverse_complement | seq == rc(s); per feature: own sequence unchanged, strand negated, sorted spans == (200-e, 200-a); 1 qualifier; undo exact |

## `copying_the_reverse_complement_keeps_every_feature_on_the_same_bases`

**Result:** ✅ passed · [line 2268](../../core/src/document.rs#L2268) · **Area:** Copy reverse complement · **Checked against:** in-test reference model · **Strength:** 🟢 strong

Features copied with features_between_reversed sit on the reverse-complemented bases with the same own sequence and opposite strand as a plain copy.

| # | Input | Expected output |
|---|---|---|
| 1 | 200 random bp (Lcg 11), top 61..100, bottom complement(101..140), split join(51..60,111..120); copy 40..160 plain and reversed | same feature count; for each: bio on rc(bases) with flipped features == bio of the plain copy; strand opposite |

## `pasted_features_cover_the_pasted_bases_and_existing_ones_shift_correctly`

**Result:** ✅ passed · [line 2293](../../core/src/document.rs#L2293) · **Area:** Copy/paste with features · **Checked against:** hand-computed · **Strength:** 🟡 adequate

Pasting at the start shifts existing features right while keeping their bases, and the pasted feature copies sit on the pasted bases.

| # | Input | Expected output |
|---|---|---|
| 1 | 200 random bp, g 51..100, h 121..150, k complement(join(11..20,31..40)); features_between(60,130); paste at 0 | 2 features copied; old g spans [(120,170)] on s[50..100]; pasted g spans [(0,40)] = s[60..100]; pasted h [(60,70)] = s[120..130] |

**Review notes:**

- Old h and the complement join k are set up but never checked after the paste, although the name says existing ones shift correctly.

## `qualifiers_are_untouched_by_edits`

**Result:** ✅ passed · [line 2313](../../core/src/document.rs#L2313) · **Area:** Feature drift / qualifiers · **Checked against:** property/invariant · **Strength:** 🟡 adequate

Random edits never alter a feature's qualifiers (note, translation, pseudo).

| # | Input | Expected output |
|---|---|---|
| 1 | CDS join(21..40,81..100) with label, note, translation 'MKVLAAGIVALLL', pseudo; 60 random replaces (Lcg 5, <=5 bases) | qualifiers equal to the original after every edit (stops if the feature disappears) |

## `making_linear_is_one_undo_step_that_rotates_and_undo_gives_the_circle_back`

**Result:** ✅ passed · [line 2335](../../core/src/document.rs#L2335) · **Area:** Topology (make linear) · **Checked against:** hand-computed · **Strength:** 🟡 adequate

Linearize at a base rotates and switches to linear in one undo step; undo restores the circle and old origin; linearize(0) keeps the start; set_topology toggles are single steps and no-ops add none.

| # | Input | Expected output |
|---|---|---|
| 1 | 100 random bp (Lcg 3), circular, g 21..40, o join(91..100,1..10), cut 31..70; linearize(50) | linear; seq = seq[50..]+seq[..50]; g and o keep their bases; 1 step; undo restores record; redo linear |
| 2 | undo; linearize(0); linearize(5) | seq unchanged; error (already linear); undo -> circular |
| 3 | set_topology(false) twice; undo | 1 history step; circular again; not modified |

**Review notes:**

- Line 2396 says linearize(50) 'cuts the feature cut (bases 31..70)', but nothing asserts what happens to 'cut' (two pieces on the right bases), the most interesting case here.

## `edited_ranges_follow_every_edit_undo_and_redo`

**Result:** ✅ passed · [line 2373](../../core/src/document.rs#L2373) · **Area:** Edited-bases highlight · **Checked against:** in-test reference model · **Strength:** 🟢 strong

The yellow 'edited bases' ranges agree with a one-flag-per-base model through 300 random replaces, overwrites, rotations/reverse complements, undos and redos.

| # | Input | Expected output |
|---|---|---|
| 1 | 40 random bp (Lcg 77 linear / 78 circular), no features; 300 random ops | after each op, edited_ranges() == runs of the flag model (deleted bases leave no flag, inserted ones are flagged; rotate/reverse permute flags); same after undo and redo |

**Review notes:**

- Line 2477 'if pos > 0 && doc.undo()' silently skips when undo() unexpectedly returns false, so a missing undo step would not fail the test (it would desync later, or not at all).

## `edited_ranges_basics_origin_and_save`

**Result:** ✅ passed · [line 2448](../../core/src/document.rs#L2448) · **Area:** Edited-bases highlight · **Checked against:** hand-computed · **Strength:** 🟡 adequate

Typed bases merge into one highlighted run that moves and shrinks with deletions, undo removes highlights, save clears them, and an edit across the origin marks the new bases.

| # | Input | Expected output |
|---|---|---|
| 1 | 'ACGT' x5 linear; insert G at 5, 6, 7 | [(5,8)] |
| 2 | delete 0..2; delete 4..5 | [(3,6)]; [(3,5)] |
| 3 | undo x2; undo x3 | [(5,8)]; empty |
| 4 | replace(3,4,'TT'); save; undo | empty after save; still empty after undo |
| 5 | circular, edit_across_origin(18,22,'GGG') | 1 range of width 3; len 19; undo empties, redo gives 1 range |

**Review notes:**

- Docstring (lines 2496-2497) says the origin edit 'marks the new bases on both sides', but lines 2524-2528 only check one range of width 3; its position is never compared with where 'GGG' actually is.
