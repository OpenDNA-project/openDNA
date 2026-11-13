# `tests/test_documents.py`

[← all test files](README.md) · source: [`tests/test_documents.py`](../../tests/test_documents.py)

## `test_saving_never_damages_the_original`

**Result:** ✅ passed in 0.54 s · [line 61](../../tests/test_documents.py#L61) · **Area:** Safe saving · **Checked against:** property/invariant · **Strength:** 🟢 strong

**Why it matters:** Without it, a save that fails (disk path missing, file locked, bad data) could truncate or corrupt the user's only copy of a file or wrongly mark unsaved work as saved.

Forces GenBank and FASTA saves to fail in several ways and checks the original file, the modified flag and the directory are untouched with no temp files left, and that a later save works. Also checks the Windows-retry rename helper and that autosave writes beside the file without counting as a save.

| # | Input | Expected output |
|---|---|---|
| 1 | save_genbank with primer site 1..99999 on a 407 bp doc | OSError; file bytes unchanged; no *.tmp; doc still modified |
| 2 | save into a missing folder | OSError; folder not created; still modified |
| 3 | save onto a path that is an existing folder with a file inside | OSError; folder and its file untouched; no *.tmp |
| 4 | Normal save afterwards | file changes, doc not modified, reload equals doc sequence |
| 5 | write_fasta of a record name not in the file; to a missing folder; onto a folder | OSError ('no sequence named' for the first); source FASTA text unchanged; no *.tmp |
| 6 | Editor: paste GATTACA, save to missing folder | OSError; editor path still ed.gb; still modified; real save then gives 'GATTACA'+'ACGT'*10 |
| 7 | fileio.replace while another handle holds the target (closed after 0.5 s) | target becomes 'new', tmp gone; on Windows it waited >= 0.4 s |
| 8 | fileio.replace onto a folder | OSError in under 2 s (no pointless retry) |
| 9 | autosave_path('plasmid.dna') and ('genome.fa','chr2'); write_genbank to autosave path after an edit | names 'plasmid.dna.autosave.gb' and 'genome.fa.chr2.autosave.gb'; user's file unchanged, doc still modified; autosave file holds the edited sequence |

**Review notes:**

- The rename-wait check only proves the wait on Windows; elsewhere it only checks the result.
- Timing-based assertions (>= 0.4 s, < 2 s) could be flaky on a loaded machine.

## `test_corrupt_and_truncated_files_are_a_clear_error`

**Result:** ✅ passed in 1.21 s · [line 144](../../tests/test_documents.py#L144) · **Area:** Corrupt files · **Checked against:** property/invariant · **Strength:** 🟢 strong

**Why it matters:** Without it, opening a truncated or damaged file could crash the app or, worse, open a partial sequence that the user then saves over the original.

Truncates a real .dna and a GenBank file at many points and checks each open either raises OSError with a message or yields the full record; feeds corrupt/foreign files and checks they are refused; checks a two-record GenBank file is refused as one document with guidance.

| # | Input | Expected output |
|---|---|---|
| 1 | pFG051_Chlor.dna and demo_insert.gb cut by 1..39 bytes and ~120 points across the file | each open is either OSError with non-empty message or a document with the full sequence and feature count; more than 100 errors per file |
| 2 | empty .gb, .dna of one byte, plain text, 1024 binary bytes, LOCUS line only, FASTA named .gb, .dna with a bad packet | OSError with a non-empty message |
| 3 | missing.gb | FileNotFoundError |
| 4 | Real .dna with every 53rd byte from 100 set to 0xFF | OSError |
| 5 | demo GenBank text concatenated twice | OSError matching '2 GenBank records' and mentioning 'reference genome' |

**Review notes:**

- 'errors > 100' is a loose count; it does not say which cut points must fail.
- Parts depend on demo files being present.

## `test_random_edits_keep_bases_and_features_then_undo_redo_exactly`

**Result:** ✅ passed in 16.91 s · [line 206](../../tests/test_documents.py#L206) · **Area:** Editing and undo / redo · **Checked against:** in-test reference model · **Strength:** 🟢 strong

**Why it matters:** Without it, some combination of edits (especially across the origin or with joins) could shift features onto the wrong bases or make undo/redo restore a different document, corrupting annotations silently.

Runs app/selftest.run_random_edits for 500 random edits per run, 3 seeds, linear and circular: after each edit the bases equal a plain-string model and each feature covers exactly the bases a per-base identity model predicts; then undoes and redoes every step comparing full snapshots.

| # | Input | Expected output |
|---|---|---|
| 1 | Seeds 0,1,2 x linear/circular, 500 random edits each (insert, delete, replace, overwrite, paste with features, feature add/delete, reverse complement; rotate and across-origin edits on circular) | after every edit: bases = string model, topology unchanged, exactly one new undo step, features match the identity model |
| 2 | Undo from the top down to the start | each snapshot equals the one recorded at that depth |
| 3 | Redo back to the top | each snapshot equals the recorded one; no further redo; returns 500 |

**Review notes:**

- The models and assertions live in app/selftest.py, not in the test, so the report reader cannot see them here; the test is a one-line wrapper.
- The app's own self-check runs shorter versions (seeds 1 and 2, 150 edits) of the same machinery: overlapping by design.

## `test_edited_base_ranges_follow_edits_undo_redo_and_saves`

**Result:** ✅ passed in 0.01 s · [line 216](../../tests/test_documents.py#L216) · **Area:** Edited-base highlighting · **Checked against:** hand-computed · **Strength:** 🟢 strong

**Why it matters:** Without it, the highlight showing which bases changed since the last save could point at the wrong bases after further edits, undo or rotation, misleading the user about what they changed.

Checks the edited_ranges list after typing, replacing, pasting and deleting; that ranges shift or shrink with later edits; that undo/redo restore them; that a save clears them but an autosave does not; and behaviour across the circular origin.

| # | Input | Expected output |
|---|---|---|
| 1 | Fresh 200 bp doc | [] |
| 2 | Type A,C,G at 10; select 50..54 type 'TT'; paste GGGGG at 100 | [(10,13)]; then [(10,13),(50,52)]; last range (100,105) |
| 3 | Type AAAA at 100, delete 0..10, type CC at 0, delete 93..95 | [(90,94)]; [(0,2),(92,96)]; [(0,2),(92,94)] |
| 4 | Type GG at 20, delete 5..8, undo, undo, redo, redo | [(17,19)] -> [(20,22)] -> [] -> [(20,22)] -> [(17,19)] |
| 5 | Type TTT at 30, write_genbank (autosave), Editor.save | still [(30,33)] after autosave; [] after save |
| 6 | Type A at 0, undo twice (past the save), redo; then modified=False | [(0,1)]; []; [] (back at saved state); [] |
| 7 | Circular: select 190..205 across origin, type 6 G; undo; redo | ranges total 6 bases; []; same ranges again |
| 8 | Type A at 0, rotate(50), then three undos | total 7 bases after rotation; [] after undoing |

**Review notes:**

- Across-origin case only checks the total length of ranges, not their positions.
- Its older sequence is built with a new random.Random(N) per base, so it is one letter repeated; the sections added on 13 Nov 2026 build varied sequences, the older ones still use it.

**Fixed by the audit of November 2026:**

- 13 Nov 2026: gap closed. The proposed checks were added, pass on the current code and fail under the mutation that slipped through before. Was: 'an edit across the circular origin survives ... rotation': after rotate(50) only the total length (7) is checked. A rotate that leaves the ranges at their pre-rotation offsets (not shifted by -50 mod length), or splits them wrongly at the origin, keeps the same total and passes. Same for the across-origin replace: the 6 G could be marked at any 6 positions.

## `test_undo_history_steps`

**Result:** ✅ passed in 0.06 s · [line 334](../../tests/test_documents.py#L334) · **Area:** Undo history · **Checked against:** hand-computed · **Strength:** 🟢 strong

**Why it matters:** Without it, the History panel could show wrong step names, primer database changes could fall out of sync with base edits on undo/redo, the DNA lock could be bypassed, or undo across records of a multi-record file could go in the wrong order.

Checks the names of annotation steps, that primer database add/rename/edit/delete/import are undo steps interleaved with base edits without setting the modified flag, that conflicting external primer changes are left alone with a message, that the DNA lock blocks base edits and their undo but not feature edits, and the FileHistory log ordering across records.

| # | Input | Expected output |
|---|---|---|
| 1 | Add feature alpha, extend it, relabel to beta, delete | History: 'Added feature alpha', 'Changed feature alpha', 'Added feature beta; Deleted feature alpha', 'Deleted feature beta' |
| 2 | On a saved doc: add primer pFwd, type A at 10, add primer pRev | pFwd matched on sequence, doc not modified; history ['Added primer pFwd','Inserted 1 bp at 11: A','Added primer pRev'] |
| 3 | undo x3, extra undo, redo x3 | pRev removed (still modified); A removed (sequence original, not modified); pFwd removed and no primer sites; 4th undo False; redo restores pFwd, the edit, then both primers |
| 4 | undo x2 then type C | redo is False; database has only pFwd |
| 5 | Rename p1->p1b, edit note, try renaming p1b to existing p2, delete p1b,p2 | steps 'Renamed primer p1 to p1b', 'Edited primer p1b'; clash refused and not recorded; 'Deleted primers p1b, p2' |
| 6 | undo x3, redo x3, undo_steps(3), undo_steps(-3) | notes restore 'edited' then 'renamed', then the original p1 row including date; redo empties; multi-step moves 3 each way |
| 7 | Import a1,a2,a3; import a1 again | added 3, step 'Imported 3 primers: a1, a2, a3'; second import exists=1 and no new step |
| 8 | a2 changed by another document, then undo import; redo | a1, a3 removed, a2 kept with its new note; error mentions a2 and 'changed since'; redo restores a1, a3 |
| 9 | Delete a1, someone adds a different a1, undo | the other a1 kept; error mentions a1 and 'used now' |
| 10 | DNA lock on after typing ACGT | type, backspace, paste refused with 'locked' message; feature add allowed and undoable; undo of typed bases refused; after unlock undo restores the sequence |
| 11 | FileHistory observe calls for records A and B (edits, undo in A, undo in B, redo in B, new edit in A), reset, steps() | undo_log/redo_log exactly as listed at each step (e.g. ['A','B','A'] then ['A','B'],['A']); reset empties; steps gives undo [(A,a1),(B,b1),(A,a2)] and redo [(A,a3),(B,b2)] |

**Review notes:**

- Five fairly independent features (annotation step names, primer undo, primer conflicts, DNA lock, FileHistory) in one test; the DNA lock and FileHistory parts could be separate tests.

## `test_copy_reverse_complement_keeps_annotations_flipped`

**Result:** ✅ passed in 0.01 s · [line 455](../../tests/test_documents.py#L455) · **Area:** Copy / paste · **Checked against:** hand-computed · **Strength:** 🟡 adequate

**Why it matters:** Without it, Copy Reverse Complement could paste the right bases but with features on the wrong strand or wrong positions, or break across the plasmid origin.

Copies a selection as reverse complement and pastes it into another document, checking bases, flipped strands and exact flipped location; checks a plain copy keeps strands; checks across-origin selection and that reverse-complementing twice gives the plain copy.

| # | Input | Expected output |
|---|---|---|
| 1 | 300 bp with gene 61..100 (+), gene complement(111..140), CDS join(51..60,121..130); select 40..160, copy RC | text = reverse_complement(seq[40:160]) |
| 2 | Paste at 10 into 50 A | bases 10..130 = RC; features top/bottom/split present with strands -1, +1, -1; top at complement(70..109); split location starts 'complement(join(' |
| 3 | Plain copy of the same selection pasted at 0 | top strand +1, bottom strand -1 |
| 4 | Circular, select 280..320 (across origin), copy RC | text = RC(seq[280:] + seq[:20]) |
| 5 | Paste that, select it, copy RC again | text = seq[280:] + seq[:20] |

**Review notes:**

- Only 'top' gets an exact location check; 'bottom' and the join's parts are checked by strand or prefix only.
- Across-origin case checks the bases but not that features are carried across the origin.
- The RC helper is the app's own reverse_complement, so a bug shared by both would not be caught (it is simple enough that this is minor).
- Its older sequence is built with a new random.Random(N) per base, so it is one letter repeated; the sections added on 13 Nov 2026 build varied sequences, the older ones still use it.

**Fixed by the audit of November 2026:**

- 13 Nov 2026: gap closed. The proposed checks were added, pass on the current code and fail under the mutation that slipped through before. Was: 'a join keeps both parts, in the new reading order' is checked only with startswith('complement(join('): part coordinates and order are free. \| 'bottom' location is unchecked (strand only). \| 'works across the origin' checks bases only; features of an across-origin selection are never inspected.

## `test_text_exports_and_the_project_file`

**Result:** ✅ passed in 0.02 s · [line 518](../../tests/test_documents.py#L518) · **Area:** Exports and project file · **Checked against:** hand-computed · **Strength:** 🟡 adequate

**Why it matters:** Without it, GFF3/BED/feature-table/protein exports could be malformed (e.g. wrong strand convention for NCBI submission) and Open DNA Project files could lose tracks, display settings or the DNA lock.

Checks the text exports of the demo vector and a built two-CDS sequence (GFF3 header and strand column, BED line count, NCBI table header and minus-strand start>stop, translations, FASTA width), the submission file pair, and that a .odp project round-trips sequence, features, tracks, display and lock, finds a moved track, and rejects a non-zip file.

| # | Input | Expected output |
|---|---|---|
| 1 | gff3_text(demo_vector) | starts '##gff-version 3'; at least features+2 lines; first CDS line strand column matches the feature |
| 2 | bed_text | one line per feature span |
| 3 | tbl_text | first line '>Feature <seq_id>'; first CDS start<stop iff plus strand; no ApEinfo or label qualifiers |
| 4 | proteins_fasta(demo) | starts with '>' |
| 5 | Built TT+ATGAAACCCGGGTAA+TT+RC(gene)+TT with CDS 3..17 and complement(20..34) | proteins ['MKPG','MKPG']; table has '3\t17\tCDS', '34\t20\tCDS' and a 'note rev' line |
| 6 | translate('ATGAAATAA'); fasta_text('x', 75 A) | 'MK*'; lines of 70 and 5 |
| 7 | write_submission | .fsa starts with '>seq_id'; .tbl exists |
| 8 | save_project with a GFF track, enzymes on, misc_feature hidden, DNA locked | zip holds exactly sequence.gb and project.json; version, dna_locked, track relative 'extra.gff3' and absolute path, display settings recorded; no *.part left |
| 9 | open_project | same sequence and feature count; metadata equals what was written |
| 10 | project_file for a track whose absolute path is gone, project moved | resolves to moved/extra.gff3 |
| 11 | open_project on 'not a zip' bytes | OSError matching 'not an Open DNA Project' |

**Review notes:**

- GFF3, BED and proteins_fasta checks on the demo vector are shallow (line counts, prefix, one strand column); a wrong coordinate would pass.
- Project round trip compares only sequence and feature count, not feature contents.
- Exports and the project file are two topics merged into one test.

**Fixed by the audit of November 2026:**

- 13 Nov 2026: gap closed. The proposed checks were added, pass on the current code and fail under the mutation that slipped through before. Was: GFF3 and BED coordinates are never checked (only header, line count and one strand column), so off-by-one starts, swapped start/end or a wrong seqid pass; a GFF3 with 0-based starts is malformed. \| Project round trip compares sequence and feature COUNT only: lost or altered qualifiers/locations, and primer sites written via file_primer_features(), are not checked.

## `test_gene_groups_and_display_switches`

**Result:** ✅ passed in 0.01 s · [line 617](../../tests/test_documents.py#L617) · **Area:** Feature display / gene groups · **Checked against:** hand-computed · **Strength:** 🟡 adequate

**Why it matters:** Without it, gene + CDS pairs might not group by locus_tag, toggling display options could re-read all features from the core (slow on genomes) or wipe Local Align hits, and type counts could go stale after edits.

Builds a 1200 bp GenBank with six gene/CDS pairs sharing locus_tags plus two misc_features, and checks grouping, per-type counts, that display switches reuse the cached features (a spy counts core reads) and keep query hits, and that counts update after adding a feature.

| # | Input | Expected output |
|---|---|---|
| 1 | Editor on the generated file | 6 gene groups; group T_F has 2 parts (CDS and tagged misc_feature); kind_counts {gene:6, CDS:6, misc_feature:2} |
| 2 | Spy on doc.features; set_display(hidden_kinds={misc_feature}) | no features() call; query_hits kept; no misc_feature in shown |
| 3 | Toggle collapse_genes, then hidden_kinds=() | still no features() call; misc_feature shown again |
| 4 | Add repeat_region 1001..1010 | kind_counts repeat_region 1, gene still 6 |

**Review notes:**

- Docstring says the switch 'also dropped the Local Align hits' while the test asserts they are kept; it refers to an old bug and reads backwards.
- Group T_F expects 2 parts although gene F has a gene, a CDS and a misc_feature; the intent (the gene feature is the group head, not a part) is not explained in the test.
- Does not check collapse_genes actually changes what is shown.
