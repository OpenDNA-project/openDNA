# `tests/test_import_export.py`

[← all test files](README.md) · source: [`tests/test_import_export.py`](../../tests/test_import_export.py)

## `test_exports`

**Result:** ✅ passed in 0.32 s · [line 38](../../tests/test_import_export.py#L38) · **Area:** Text/file exports (GFF3/BED/TBL/FASTA/GenBank) · **Checked against:** hand-computed · **Strength:** 🟢 strong

All export formats (GFF3, BED, feature table, proteins, FASTA, submission pair, GenBank/FASTA copies) produce correct content with correct strand signs and the Export menu greys unavailable items, so users get valid downstream files.

| # | Input | Expected output |
|---|---|---|
| 1 | demo_vector editor: gff3_text, bed_text, tbl_text, proteins_fasta | gff starts '##gff-version 3', lines>=features+2; CDS line strand +/-; bed line count==sum spans; tbl header '>Feature <id>', minus-strand parts start>stop, no ApEinfo/label; proteins_fasta starts '>' |
| 2 | built seq with fwd CDS 3..17 and complement(20..34) of 'ATGAAACCCGGGTAA'; translate('ATGAAATAA'); fasta_text('x','A'*75) | protein_of both == 'MKPG'; tbl has '3\t17\tCDS','34\t20\tCDS','\t\t\tnote\trev'; translate=='MK*'; fasta wraps to ['A'*70,'A'*5] |
| 3 | write_submission; export_fasta/export_genbank via monkeypatched dialog; _update_export_menu in several states | fsa starts '>'+seq_id, tbl exists; fasta ends with last 20 bases, genbank starts 'LOCUS', path unchanged (export is a copy); menu enables whole/features, disables alignments (no BAM) and selection until a region selected; empty window all disabled |

## `test_dna_lock_and_project_file`

**Result:** ✅ passed in 0.23 s · [line 95](../../tests/test_import_export.py#L95) · **Area:** DNA lock + .odp project round trip · **Checked against:** round-trip/self-consistency · **Strength:** 🟢 strong

The DNA lock refuses base edits and linearising (but allows feature edits) and unlocks again, and an .odp project zip round-trips the sequence, tracks, display settings and lock state.

| # | Input | Expected output |
|---|---|---|
| 1 | lock editor; try type_bases/backspace/paste; add+undo a feature; undo typed bases; set_circular(False) | base edits return False with 'locked' message; feature add works and undoes; locked base-undo refused (seq != original); linearising refused (circular unchanged) |
| 2 | unlock then undo; add gff3 track, set_display(enzymes,hidden_kinds), lock, write .odp | unlocked undo restores seq; _write True, path==project, not modified; zip holds sequence.gb+project.json; meta dna_locked True, tracks[0].relative=='extra.gff3' |
| 3 | reopen project in new MainWindow | sequence+feature count match; dna_locked True, lock_action checked iconText 'DNA'; tracks ['extra.gff3']; show_enzymes True; hidden_kinds {'misc_feature'} |

## `test_reference_genomes`

**Result:** ✅ passed in 10.11 s · [line 156](../../tests/test_import_export.py#L156) · **Area:** NCBI reference genome download (network faked) · **Checked against:** independent implementation · **Strength:** 🟡 adequate

Searching NCBI Datasets, downloading (and caching) a genome, and opening it with its GFF track and setting it as the off-target reference all work with faked network responses.

| # | Input | Expected output |
|---|---|---|
| 1 | monkeypatched ncbi.fetch returns REPORT; search_genomes('Madeup coli') and accession form | taxon URL with reference_only; assembly (accession,length 60,sequences 1,level 'Complete Genome'); accession query uses accession URL; blank query raises NcbiError |
| 2 | monkeypatched open_url returns fake_zip(); download_genome(a) | not cached, sizes recorded, fasta starts '>chr1', gff has 'abcA', URLs have GENOME_FASTA+GENOME_GFF, no .part left; second call cached, only 1 download |
| 3 | GenomeDialog: search, wait for assemblies, open, use_as_reference (network faked, QMessageBox.Yes) | reference_button enabled (1 seq); opened is window, length 60, track 'GCF_000000001.1.gff'; genome_path ends '.fna' and exists; primers/genome setting saved/restored in finally |

**Review notes:**

- Three busy-wait loops (191,196,201) with 10s timeouts and QTest.qWait could flake under heavy load; they fall through to a normal assert on timeout, which is correct but time-dependent. Shared QSettings 'primers/genome' is mutated but restored in finally.

## `test_example_files_and_feature_direction`

**Result:** ✅ passed in 0.45 s · [line 212](../../tests/test_import_export.py#L212) · **Area:** Example files + feature direction/arrows · **Checked against:** hand-computed · **Strength:** 🟡 adequate

The bundled example plasmids and genome open with features, .dna CDS features keep direction while misc_features come in arrowless, and the feature dialog defaults/edits direction correctly (writing the qualifier only when non-default).

| # | Input | Expected output |
|---|---|---|
| 1 | open pBAD30_eilR.gb(5936), pFG051_Chlor.gb(3780), BW25113_2026.gb(>4e6) | each has features and expected length |
| 2 | (if .dna present) read direction; is_directional checks | CDS keep arrow!=0, some misc_feature arrow==0; is_directional default misc_feature False / CDS True; overridden by 'directional' qualifier |
| 3 | FeatureDialog new misc_feature 'plain', then edit direction; new CDS 'arrowless' | default direction 'None...', arrow 0, no 'directional' written; edit to index1 -> arrow -1, ('directional','yes'), location complement(...); new CDS default 'Forward', set index2 -> arrow 0; views grab without error |

**Review notes:**

- .dna direction section (224) guarded by if DNA_FILE.exists(): silently skipped when sample absent. Final grab() calls (257-259) are smoke-only (no assertion on rendering).

## `test_primer_display_and_hiding`

**Result:** ✅ passed in 3.93 s · [line 263](../../tests/test_import_export.py#L263) · **Area:** GUI: primer hide + display source filter · **Checked against:** GUI state · **Strength:** 🟡 adequate

'Hide in This File' removes a primer from drawing but keeps it found/listed/remembered per file, and the Primers display filters ('from this sequence' vs 'from the database') switch independently, with a primer in both counting as the file's.

| # | Input | Expected output |
|---|---|---|
| 1 | add DB primer 'hideme'=seq[100:122]; match in background; set_primer_hidden True/False | found and drawn; hidden -> not drawn but still found+in DB; PrimerPanel item 'hideme (hidden)'; remembered in a fresh Editor's hidden_primers; unhidden -> drawn again |
| 2 | primers fromdb/fromfile/both; toggle db_primers/file_primers; DisplayPanel | all three drawn; db off -> fromdb gone; file off (db on) -> only fromdb (shared counts as file's); both off -> none, show_primers False; DisplayPanel file_primers restores; a file-only primer's 'Hide in This File' action enabled |

**Review notes:**

- Two background-match busy-wait loops (276,303) with 10s timeouts (QTest.qWait); potential flakiness under load. Cleanup of shared primer DB done in finally (326-328), good.

## `test_tab_names_and_standard_layout`

**Result:** ✅ passed in 0.60 s · [line 334](../../tests/test_import_export.py#L334) · **Area:** GUI: tab naming + standard layout · **Checked against:** GUI state · **Strength:** 🟡 adequate

A workspace tab is named after its visible panes (hiding/showing Alignment renames it), and every file opens with the standard 70/30 one-tab layout while only an .odp project keeps a custom layout.

| # | Input | Expected output |
|---|---|---|
| 1 | demo open; set_genome_visible True/False | tab 'Sequence · Map · Features'; after show starts 'Alignment' with tooltip; after hide pane kept but not visible, 'Alignment' gone from title+tooltip |
| 2 | fresh window loads demo at 1400x850 | root sizes ~[70,30]; inner [50,50]; panes [Sequence,Map,Features]; exactly 1 tab |
| 3 | save rearranged as .odp and reopen; reopen plain demo | project opens with 'Alignment' tab as saved; plain file is 'Sequence · Map · Features' |

## `test_tab_image_export`

**Result:** ✅ passed in 4.97 s · [line 380](../../tests/test_import_export.py#L380) · **Area:** GUI: tab image export (PNG/TIFF/SVG) · **Checked against:** GUI state · **Strength:** 🟡 adequate

The image-export dialog keeps scale and DPI in step, and the current tab exports to PNG/TIFF/SVG at the chosen size/DPI with views drawn as vector paths (not cached pixmaps) in SVG.

| # | Input | Expected output |
|---|---|---|
| 1 | ImageExportDialog dpi 600 then scale 100 | scale == round(600/96*100); options()[1:3]==(1.0,600); 'cm' in info |
| 2 | export_tab_image for png/tiff/svg at scale 2.0, 300 DPI (monkeypatched options/dialog) | returns path, size>1000, ACTIVE reset; png/tiff: image size == 2x root, dotsPerMeter~300, >3 distinct pixels |
| 3 | SVG content inspection | has <path>, image widths<=64 (icons only); font sizes < 2x screen size; svg width_mm == root.width*2/300*25.4; <text> tags have fill-rule='nonzero' |
