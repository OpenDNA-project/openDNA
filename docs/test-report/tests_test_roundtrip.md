# `tests/test_roundtrip.py`

[← all test files](README.md) · source: [`tests/test_roundtrip.py`](../../tests/test_roundtrip.py)

## `test_real_file_round_trips`

**Result:** ✅ passed in 4.26 s · [line 119](../../tests/test_roundtrip.py#L119) · **Area:** GenBank/.dna file I/O round trip · **Checked against:** independent implementation · **Strength:** 🟢 strong

A real plasmid/.dna/genome file loaded, saved as GenBank and reloaded must be identical (bases, feature locations+order, qualifiers, topology, name, definition); re-saving must be byte-stable; and Biopython must read our file the same, so users never lose annotation on save.

| # | Input | Expected output |
|---|---|---|
| 1 | every present file in genbank_sources()+genome_dna (pBAD30_eilR.gb, pFG051_Chlor.gb/.dna, demos, BW25113 genome); skip if none | snap(reloaded)==snap(original) for every key; identical f.spans and f.strand; back.modified False; also asserts original has features |
| 2 | save->reload->save->reload->save of each source to 1.gb/2.gb/3.gb | first==second==third bytes (no drift across saves) |
| 3 | each saved file re-read by Biopython via assert_biopython_agrees | seq upper-cased equal, topology, id/name, equal feature count, bio_parts==our_parts, bio_qualifiers==our_qualifiers |
| 4 | demo/BW25113_2026.gb | length > 4 500 000, circular, > 4000 features |
| 5 | save_genbank, reopen | same sequence and same (kind, location, qualifiers) for every feature |
| 6 | save the reopened document again | byte-identical to the first save |

**Review notes:**

- pytest.skip at line 124 when no real/demo files present: silently passes in a CI checkout lacking demo/ (though files are committed here).
- Only self-consistency: the saved file is never compared with the original or read by an independent parser, so a bug that loses data the same way on load and save would pass. Runtime skip at line 121 if the genome is missing.

**Fixed by the audit of November 2026:**

- Removed in the second pass of the audit (fewer tests to review, no biological check lost): tests/test_real_files.py::test_genome_opens_and_round_trips was merged into this test; the genome is now also in the byte-identical save loop, and its size, topology and feature count are checked here.

## `test_rich_document_round_trip`

**Result:** ✅ passed in 0.00 s · [line 158](../../tests/test_roundtrip.py#L158) · **Area:** GenBank I/O of hard location forms · **Checked against:** independent implementation · **Strength:** 🟢 strong

A synthetic record with joins, complements, fuzzy ends, order() and awkward qualifiers survives save/reload byte-identically and is read the same by Biopython, both circular and linear.

| # | Input | Expected output |
|---|---|---|
| 1 | rich_document(seed=3, 400 bp, 8 features incl join/complement/fuzzy/order/pseudo) | snap(reload)==snap(doc); re-save byte-identical to first |
| 2 | same doc, circular in (True,False), via assert_biopython_agrees | Biopython agrees; 'fuzzy both ends' start is BeforePosition, end is AfterPosition; 'complement of a join' location.strand == -1 |

## `test_topology_survives_genbank_save_and_load`

**Result:** ✅ passed in 0.00 s · [line 184](../../tests/test_roundtrip.py#L184) · **Area:** GenBank topology (circular/linear) · **Checked against:** round-trip/self-consistency · **Strength:** 🟢 strong

Topology is written to the LOCUS line, reloaded correctly, and a post-load toggle is saved and is a single undo step, so circular vs linear is never silently lost.

| # | Input | Expected output |
|---|---|---|
| 1 | 'ACGTACGTTTGACCA'*7 (105 bp), circular in (True,False) | ' circular ' in text iff circular, ' linear ' iff not; reloaded circular matches and sequence preserved |
| 2 | set_topology(not circular), save, reload; then undo | reloaded circular == not circular; undo() True and back.circular restored to original |

## `test_location_forms_survive_a_genbank_round_trip`

**Result:** ✅ passed in 0.01 s · [line 209](../../tests/test_roundtrip.py#L209) · **Area:** GenBank feature location syntax · **Checked against:** independent implementation · **Strength:** 🟢 strong

Every GenBank location form (ranges, complement, join, order, fuzzy <>, cross-origin, and the 20^21 between-base gap) survives a round trip written in the form given and matches Biopython's parsing of the span.

| # | Input | Expected output |
|---|---|---|
| 1 | 13 LOCATIONS on a circular 'ACGT'*100 doc, each as a lone misc_feature | (location,spans,strand) unchanged; after.location==location (written verbatim); Biopython part spans == after.spans |
| 2 | '20^21' on linear 'ACGT'*20 | location text '20^21' survives; Biopython start==end==20; f.spans==[(20,20)] |

## `test_qualifier_text_survives_and_is_written_as_documented`

**Result:** ✅ passed in 0.03 s · [line 243](../../tests/test_roundtrip.py#L243) · **Area:** GenBank qualifier text wrapping/escaping · **Checked against:** independent implementation · **Strength:** 🟢 strong

Awkward qualifier text (doubled quotes, whitespace runs collapsed, long wrapping at spaces, unicode, flag with no value, translation whitespace stripped) is written canonically within 80 bytes/line and read back identically by us and Biopython; very long words are never cut; a value filling the last line has no trailing space.

| # | Input | Expected output |
|---|---|---|
| 1 | 11 qualifiers incl 'say "hi" and ""twice""', newlines/tabs, words(70)~500ch, unicode*12, pseudo=None on misc_feature 10..90 | back.qualifiers == expected (whitespace collapsed to single space; translation whitespace removed); max line <=80 bytes; '日本語' and '/pseudo' present; re-save byte-identical; Biopython note/translation/pseudo/product match |
| 2 | 'x'*300 word in note (before/after and alone) | both notes read back exactly (word not broken even past column width) |
| 3 | note 'sequence: '+dna(n,n) for n in 40..89 on primer_bind | our reader and Biopython read note exactly (no trailing space); every line <=79 unless an unbreakable word |

## `test_primers_survive_save`

**Result:** ✅ passed in 0.35 s · [line 315](../../tests/test_roundtrip.py#L315) · **Area:** Primer features in GenBank save · **Checked against:** round-trip/self-consistency · **Strength:** 🟢 strong

primer_bind features keep name and sequence through a save; .dna primers survive export; and Editor.save writes a primer_bind only for primers imported into the file (not DB-only ones), adding no extra undo step.

| # | Input | Expected output |
|---|---|---|
| 1 | doc.save_genbank(out, [two primer_bind feats with note 'sequence: ...']) on a 124 bp seq | primer_features==[('fwd','TTTTT...','5..24'),('rev','cgatccgg...','complement(60..80)')]; document itself keeps 0 features |
| 2 | pFG051_Chlor.dna (if present): save as GenBank | >=5 named primers; back==primer_features(doc) and non-empty; {name:seq} equals doc.primers() |
| 3 | Editor on 130 bp doc, DB primer 'myprimer'=seq[5:27]; save; import_primers; save | DB-only primer not written; after import exactly 1 imported; reopened sites all 'myprimer' with full seq, one at '6..27'; doc not modified, history()[0] length 1 (no extra undo) |

**Review notes:**

- .dna section (330-338) guarded by if PLASMID_DNA.exists(): silently skipped when the .dna sample is absent.

## `test_tiny_and_empty_sequences_save_and_reload`

**Result:** ✅ passed in 0.01 s · [line 361](../../tests/test_roundtrip.py#L361) · **Area:** GenBank edge-size sequences · **Checked against:** round-trip/self-consistency · **Strength:** 🟢 strong

One-base, mixed-IUPAC, empty, and edited-to-nothing sequences save and reload without error in both topologies, protecting against crashes on degenerate inputs.

| # | Input | Expected output |
|---|---|---|
| 1 | bases in ('A','g','N','AC','acgtNRYKMSWBDHVU') x circular in (False,True), one full-length misc_feature | reloaded (sequence,length,circular,name) match; snap equal; Biopython seq upper-cased equal, 1 feature |
| 2 | empty sequence '' x both topologies | reloaded (length,sequence,circular,name)==(0,'',circular,'empty') |
| 3 | 'ACGT' deleted 0..4 then save | file exists, doc not modified, length 0 |

## `test_gap_letter_and_letter_case_survive`

**Result:** ✅ passed in 0.00 s · [line 395](../../tests/test_roundtrip.py#L395) · **Area:** GenBank base case / gap character · **Checked against:** round-trip/self-consistency · **Strength:** 🟢 strong

A '-' gap letter survives, and lower/mixed-case bases are preserved as typed, written in 10-base blocks, stable across saves, and an edit keeps typed case, so users' case annotations are not clobbered.

| # | Input | Expected output |
|---|---|---|
| 1 | 'ACGT-ACGT' save/reload | sequence()=='ACGT-ACGT' |
| 2 | hand-written GenBank with mixed-case ORIGIN (60 bp) | sequence preserves exact case; slice (8,14)=='acGTAC'; re-save byte-identical; ORIGIN blocks ['acgtacgtac','GTACGTACGT','nnnnRYKMac',...]; Biopython seq == upper() |
| 3 | insert(0,'TTT') | sequence[:6]=='TTTacg' (typed case kept) |

## `test_fasta_reading_and_export`

**Result:** ✅ passed in 0.00 s · [line 455](../../tests/test_roundtrip.py#L455) · **Area:** FASTA I/O and multi-record index · **Checked against:** independent implementation · **Strength:** 🟢 strong

FASTA export is read by Biopython with the same id/seq; the multi-record Fasta index fetches each record and arbitrary substrings matching Biopython; and line-ending/blank-line/case variants all yield the same bases.

| # | Input | Expected output |
|---|---|---|
| 1 | RECORDS=[chr1(250),chr2(130),chr3(61),empty_one(0)]; from_fasta chr2; write_fasta | doc=(chr2,130,linear,RECORDS[1]); exported ids [chr1..empty_one]; seqs match; description 'chr1 first one'; edit then export still matches; doc.modified True |
| 2 | Fasta(multi.fa).sequences() and fetch of full and substrings (0,1),(59,61),(55,130),(249,250),(119,121),(0,250) incl end past length | sequences()==[(chr1,250),(chr2,130),(chr3,61),(empty_one,0)]; every fetch equals Biopython slice; clipped at end; .fai written and reused; from_fasta per record matches |
| 3 | 8 variants (crlf, no final newline, trailing blanks, lower/mixed case, one long line) each its own file | sequences() stable; fetch upper==bases and ==Biopython; case kept from file; from_fasta matches; cross-line fetch (148,152) correct |

## `test_fasta_writing_after_an_edit`

**Result:** ✅ passed in 0.00 s · [line 519](../../tests/test_roundtrip.py#L519) · **Area:** FASTA in-place edit write · **Checked against:** independent implementation · **Strength:** 🟢 strong

Writing FASTA after edits equals the same edits applied to the file string (to a new file and in place), leaving other records byte-for-byte untouched even with CRLF, so saving one chromosome never corrupts its neighbours.

| # | Input | Expected output |
|---|---|---|
| 1 | 3 non-empty RECORDS; from_fasta chr2; insert(5,'GGGGGGG'), delete(100,110), replace(0,3,'tt') | doc.sequence()==same edits on plain string; write to new file equals expected and leaves src unchanged; in-place write matches and replaces stale index; chr1/chr3 intact; no *.tmp left |
| 2 | CRLF file RECORDS[:3]; from_fasta chr2; insert(0,'AAAA'); write in place | chr1 and chr3 copied byte-for-byte (CRLF), chr2 rewritten with \n; Fasta fetch chr2=='AAAA'+orig; Biopython records == expected dict |

## `test_a_failed_save_leaves_the_original_alone`

**Result:** ✅ passed in 0.00 s · [line 567](../../tests/test_roundtrip.py#L567) · **Area:** Atomic save / failure safety · **Checked against:** property/invariant · **Strength:** 🟢 strong

A failed GenBank save (bad primer site, missing folder, folder-in-the-way) leaves the original file and document's modified flag untouched and a later real save works; likewise for FASTA and for the Editor keeping its path and modified flag.

| # | Input | Expected output |
|---|---|---|
| 1 | save rich_document, edit, then save_genbank with primer 'primer_bind 1..99999' | raises OSError; file unchanged, no *.tmp, doc still modified |
| 2 | save to nonexistent subfolder; save over a directory 'taken.gb' containing a file | raises OSError each; folder not created / directory+contents untouched; doc still modified; subsequent real save succeeds, doc not modified, reload equal |
| 3 | write_fasta to source for a doc not in it; to missing folder; Editor.save to missing folder after paste | OSError match 'no sequence named'; src bytes unchanged; Editor keeps path and modified; later save OK and sequence 'GATTACA'+... |

## `test_autosave_and_recovery`

**Result:** ✅ passed in 1.67 s · [line 649](../../tests/test_roundtrip.py#L649) · **Area:** Autosave and crash recovery · **Checked against:** GUI state · **Strength:** 🟢 strong

Autosave writes a .autosave.gb recovery copy without ever touching the user's original (GenBank and .dna), drops declined/stale copies, and offers+recovers unsaved edits after a simulated crash, protecting against data loss.

| # | Input | Expected output |
|---|---|---|
| 1 | open demo_vector.gb copy in MainWindow; autosave with no change; then paste 'GATTACA' + autosave; second edit 'CCC' | no recovery file when unchanged; recovery at vector.gb.autosave.gb with edited seq (len+7), original bytes intact, doc still modified; copy follows second edit; close removes recovery |
| 2 | edit+autosave, simulate crash; reopen with Discard (QMessageBox.No); reopen with stale recovery older than file | declined: asked once, recovery removed, original opens unmodified; stale: asked==[] (dropped without prompt), recovery removed |
| 3 | (if .dna present) plasmid.dna: edit+autosave | recovery named plasmid.dna.autosave.gb, .dna untouched, recovery seq matches editor |
| 4 | edit 'TTTTTTTT'+autosave, crash, reopen with Yes, then save() | asked once with 'unsaved changes'; recovered seq==edited, modified, feature count preserved, original still untouched; save() makes it permanent and removes recovery |
| 5 | load plasmid.gb; autosave timer | active, interval AUTOSAVE_MS == 60000 |
| 6 | autosave with no change | recovery 'plasmid.gb.autosave.gb' not written |
| 7 | type 'GATTACA', autosave; autosave again; type 'C', autosave | recovery exists, original unchanged; not rewritten when unchanged (mtime equal); recovery holds GATTACAC once |
| 8 | save() | recovery removed; original changed; autosave_file None |
| 9 | crash: edit + feature, autosave, close without deleting recovery; make recovery newer; load with _ask_recover->True | asked once; recovered sequence and features match edited; doc modified, title starts '*'; original still unchanged; normal close removes recovery |
| 10 | decline recovery (_ask_recover->False) | original opens unmodified, != edited; recovery deleted |
| 11 | recovery older than file (_ask_recover must not be called) | not offered; recovery deleted; original unchanged |

**Review notes:**

- work_file fixture pytest.skip when no demo files (whole test skipped); .dna branch (697) silently skipped when sample absent. Relies on monkeypatched _ask_recover (no real dialog), which is appropriate.

**Fixed by the audit of November 2026:**

- Removed in the second pass of the audit (fewer tests to review, no biological check lost): tests/test_app_logic.py::test_autosave_and_recovery was merged into this test; its checks that are not repeated here (the 60 s timer, no rewrite without a change, a feature added by the user recovered, the * title) were added.

## `test_a_truncated_file_is_a_clear_error_never_a_crash_or_a_silent_partial_document`

**Result:** ✅ passed in 0.07 s · [line 749](../../tests/test_roundtrip.py#L749) · **Area:** Truncated file robustness · **Checked against:** property/invariant · **Strength:** 🟢 strong

Truncating a .dna or GenBank file at ~150 positions must either raise a clear OSError (with a message) or return the complete record, never crash or silently produce a partial document.

| # | Input | Expected output |
|---|---|---|
| 1 | PLASMID_DNA and a GenBank source, truncated at cut_points(len) (~first 40 + spread), each reopened | either OSError with non-empty str, or doc equals full sequence and full feature count; errors>100 per label (confirms most truncations are rejected, not partially read) |

**Review notes:**

- Lines 742-744: both sources guarded by if not source.exists(): continue with no assertion that at least one ran; if neither sample exists the test passes with zero checks (vacuous). errors>100 guard only applies to files that do exist.

## `test_corrupt_and_foreign_files_give_an_error_message`

**Result:** ✅ passed in 0.17 s · [line 769](../../tests/test_roundtrip.py#L769) · **Area:** Corrupt/foreign file errors · **Checked against:** property/invariant · **Strength:** 🟢 strong

Corrupt, empty, foreign, and multi-record files raise OSError with a message (and a missing file FileNotFoundError), a 2-record GenBank is refused rather than dropping a record, and opening a corrupt file in the window shows an error while keeping the current document.

| # | Input | Expected output |
|---|---|---|
| 1 | 7 cases: empty.gb b'', empty.dna, non-sequence text, binary, LOCUS-only, FASTA-as-.gb, bad .dna packet | each raises OSError with non-blank message; missing.gb raises FileNotFoundError; scrambled real .dna (if present) raises OSError |
| 2 | demo file text duplicated into two.gb | raises OSError matching '2 GenBank records' |
| 3 | truncated demo bytes[:300] loaded into a window that already has a demo open (QMessageBox.critical monkeypatched) | win.load returns False; one critical dialog with 'Cannot Open File' and 'bad.gb'; open document's name/length unchanged |

**Review notes:**

- Line 792 'if not DEMOS: return' silently skips the two-record and window-keeps-document checks when demo files are absent (no skip marker).
