# `tests/test_gui.py`

[← all test files](README.md) · source: [`tests/test_gui.py`](../../tests/test_gui.py)

## `test_plasmid_edit_features_find_and_autosave`

**Result:** ✅ passed in 1.46 s · [line 71](../../tests/test_gui.py#L71) · **Area:** GUI: editing a plasmid end to end · **Checked against:** GUI state · **Strength:** 🟢 strong

Protects the core editing loop on the real plasmid pBAD30_eilR.gb (5,936 bp): selecting, typing, undo, saving a .dna, adding/deleting features, reverse-complement, Find and autosave recovery, so a user never loses edits or corrupts a file.

| # | Input | Expected output |
|---|---|---|
| 1 | Select row 3 of the feature table | ed.current_feature set and (start,end) == feature_bounds of that feature |
| 2 | Click base 5 on the top row; type 'gattaca'; Backspace; type 'qzx!1'; undo all | cursor at 5; length n+7 and sequence(5,12)=='gattaca'; title starts with '*'; n+6 after Backspace; non-DNA keys ignored (n+6); after undo length n and not modified |
| 3 | Cursor 10, Shift+Right x4 | selection (10,14) |
| 4 | Open pFG051_Chlor.dna, type 'A', save (save dialog stubbed to saved.gb) | save() True, saved.gb exists, original .dna bytes unchanged, not modified, editor.path == saved.gb |
| 5 | Rows of the sequence view | row_at(offset(r)+1)==r for all rows; more than one distinct row height; ed.unique_sites non-empty |
| 6 | FeatureDialog on selection (9,30): name 'my site', Reverse, '/note="added in test"' | (label,location,strand) == ('my site','complement(10..30)',-1); ('note','added in test') in qualifiers |
| 7 | Edit that feature with location '5..999999' (warning stubbed) | feature count n+1 unchanged, location still 'complement(10..30)' |
| 8 | Select its table row, press Delete, then undo | feature count back to n, then n+1 |
| 9 | reverse_complement() twice | first 20 bases == revcomp of old last 20; after second, same sequence and same (kind, bases covered, strand) set |
| 10 | Find revcomp of sequence(500,515) from cursor 0; Find 'my sit'; translation toggle; autosave then reload with recover=True; close | selection (500,515); current_feature == new feature; total row height grows; recovery file exists and work file unchanged; reloaded doc modified with n+1 features; recovery file removed on close |

**Review notes:**

- L113-115: 'Clicking the map selects something' (named in the docstring) has no assertion at all: click + grab only, smoke.
- L128-129: redundant pytest.skip (open_window at L74 already skipped if the plasmid is missing).
- L167-171: the Delete shortcut depends on the window becoming active; qWaitForWindowActive is 'best effort' (comment L169), so this step can be flaky when other tests leave windows open.

## `test_primers_clipboard_design_pcr_orfs_and_digest`

**Result:** ✅ passed in 7.93 s · [line 218](../../tests/test_gui.py#L218) · **Area:** GUI: clipboard, primer database, primer design, PCR, ORFs, digest · **Checked against:** GUI state · **Strength:** 🟢 strong

Checks copy/paste with annotations, DNA-only typing, the primer database (binding sites, mismatches, tails, checks, import/export, saving as primer_bind), primer design, PCR products, ORF detection, digests and primer secondary structure on pBAD30_eilR.

| # | Input | Expected output |
|---|---|---|
| 1 | Copy the first CDS (Ctrl+C) and paste at 0; paste into another window with Editor.clip=None | length n+size, 2 copies of the label, one with spans [(0,size)]; undo restores n and feature count; other window gets 2 copies via system clipboard |
| 2 | paste('>x\n1 acgt acgt\n'); paste('hello world'); type 'acgtnr' | sequence(0,8)=='acgtacgt'; 'hello world' refused, length n; only 5 bases typed (R refused) |
| 3 | text_color_on('#1a237e'), text_color_on('#ffd54f') | white, black |
| 4 | PrimerDialog 'testF' = 'gaattc' + sequence(2000,2022) | report has '2,001–2,022 → exact'; DB stores it; site spans [(2000,2022)], strand 1, 0 mismatches, tail 6; not a primer_bind feature; file primer oAC1345 in sites; 5 bases typed at 100 move site to (2005,2027) |
| 5 | testR = revcomp of sequence(2300,2325) with base 8 changed | spans [(2300,2325)], strand -1, 1 mismatch, mismatch_at [2308] |
| 6 | Primers panel group testF; click its row; click its primer lane; panel.reload() | collapsed, 1 exact row, Tm == melting_temp within 0.05; current_primer == testF and selection (2000,2022); 'Tm' in status; stays expanded after reload |
| 7 | primers.csv with testR2 (revcomp 2600..2625) and junk; export FASTA; Primer Database tab | '2 new'; testR2 spans [(2600,2625)]; junk in DB but no site; export has '>testF\ngaattc'; last parsed name 'junk'; tab lists testF,testR,junk,oAC1345 and all DB rows |
| 8 | binding_report / primer_checks on sequence(3000,3020); after duplicating those 20 bases; tm_low=Tm+1, margin 0.5 | 'Binds here: 3,001–3,020 → exact' (+ genome line if BW25113 present); 'ok' for Non-specific and Binds here; 'bad' when duplicated; Tm 'bad' |
| 9 | Save with oAC1346 redefined in DB as 'ACGTTTGACCAGTTTACGGA' | primer_bind testF at '2001..2022' with note 'sequence: gaattc…'; testR (1 mismatch) not written; oAC1346 kept once with original sequence containing 'GTAGAGGATCTGCTCATG' |
| 10 | Design primers for selection 1000..1200; design_for_selection; PCR des_F/des_R; mutF with mismatch at 12 (1500..1525) + mutR (1700..1725); OrfDialog; digest BamHI+PvuI; fragments(100,False,{30:['X']}); structure() | fwd end<=1000, rev start>=1200, Tm within 5 °C of 60, seqs match template; temp sites pair1_F/R then pair2; save -> des_F/des_R; product == seq[f.start:r.end], untitled+modified, has features; mutant product starts with mutant, length 225; a >100-codon CDS found as ORF; added ORF labelled 'ORF…'; 2 fragments summing to length, 2 gel lanes; [(30,0,30,[],['X']),(70,30,100,['X'],[])]; '3′ self-dimer' warning for ACGTTTAGCGAATTCGC, none for AGCGGATAACAATTTCACACAGGA |

**Review notes:**

- L383-398: the shared primer DB is altered (oAC1346 replaced) and only restored at L397-398 without try/finally; a failure in between leaves a wrong oAC1346 for every later test. The DB (primer_db.db()) is a process-wide singleton, and testF/testR/junk/mutF/des_F stay in it for later tests.
- L375-379: primer-check settings are reset to hard-coded 52.0/7.0 (the defaults) without try/finally.
- L444-448: the PCR dialog's combo boxes are set, then PcrDialog.pair/exec/_pick_default are monkeypatched, so the dialog's own pair selection is not tested, only pcr_product/new window.
- L352: OptionsDialog(win) is constructed only (smoke).
- L474: only requires that any one long CDS is found as an ORF; L416 Tm tolerance ±5 °C is loose.

## `test_fasta_file_tracks_genome_view_search_and_chromosome_list`

**Result:** ✅ passed in 2.58 s · [line 498](../../tests/test_gui.py#L498) · **Area:** GUI: multi-FASTA, genome view, tracks, whole-file search · **Checked against:** GUI state · **Strength:** 🟢 strong

Protects work on multi-sequence FASTA files: editing one record without touching others, GFF/BED tracks in the genome (alignment) view, track display controls, file-wide search, and the chromosome list ordering.

| # | Input | Expected output |
|---|---|---|
| 1 | genome.fa: chrA 150,000 random bases (seed 7), chrB 20,000 random + 'GGTACCTTGGTACC' | record chrA, length 150,000; picker visible with 2 entries; genome view shown (>100 kb) |
| 2 | Type 'NNNN' at 10 and save | Fasta.fetch('chrA',0,20) == chr_a[:10]+'NNNN'+chr_a[10:16]; chrB unchanged |
| 3 | Add gene 101..200 and save (dialogs stubbed) | info shown, chrA.gb written, fasta_record None |
| 4 | GFF3 with 140 genes at 1000*i+1..1000*i+800 on chrA; locus 'gene42'; '1,000-2,000'; zoom 0.001 | track.seqid chrA; selection (42000,42800) inside view; (999,2000); a 'Bases' row drawn |
| 5 | Click track item gene10 | selection (10000,10800) |
| 6 | Click track name twice; condensed; drag edge +40 px; log10; load GFF+BED via Load Tracks | mode 'hidden' with HIDDEN_H, then original height; condensed smaller; height small+40 and minimumHeight >= layout height; 2 tracks loaded, 3 total |
| 7 | perf.STATS after a selection change | tracks not redrawn, only overlay stat recorded; timings window visible with >=1 row |
| 8 | Legend button / close box | window visible with 'mate unmapped' and '80 %'; closing unchecks button; toggles show/hide |
| 9 | Search whole file for 'GGTACC'; go to hit chrB:20000 | hit count == Python str.count over both records, 2 records searched; chrB opened with selection (20000,20006); counter 'of N'; search spans highlighted |
| 10 | two.fa (chrA 2000 bp, chrB 1500 bp) record_requested('chrB'); unsorted.fa small/big/mid | chrB loaded, length 1500, picker 'chrB'; records and picker ordered ['big','mid','small'], 'big' opened |

**Review notes:**

- L626-627: waits up to 10 s for streamed search hits (timing-dependent on slow machines).
- L516, L547, L556, L593: canvas grabs are paint smoke checks only.

## `test_workspace_tabs_and_view_picker`

**Result:** ✅ passed in 1.46 s · [line 661](../../tests/test_gui.py#L661) · **Area:** GUI: workspace tabs, panes and view picker · **Checked against:** GUI state · **Strength:** 🟢 strong

Checks the tabbed workspace (tab naming, lazy panes, split/close, '+' button, middle-click copy, view switching without crashes, default tabs for new windows) and the grouped view-picker popup and theme colours.

| # | Input | Expected output |
|---|---|---|
| 1 | Fresh window on pBAD30_eilR; add tab Sequence/Detect ORFs/Digest | tabs ['Sequence · Map · Features']; new tab 'Sequence · Detect ORFs · Digest'; its panes have view None until shown |
| 2 | Type 'ATG'+'GCT'*150+'TAA' at 0, wait 600 ms | (0,456,1) in ORF pane list and more ORFs than before |
| 3 | Split Sequence pane, set to Map, select (100,200), close pane | MapView visible, both views share the editor; no MapView left |
| 4 | Click '+'; middle-click tab 0 | '+' right of last tab; 3 tabs, current 2, Sequence pane; copy tab right after with same text and pane views |
| 5 | Switch a pane Alignment/Sequence/Alignment/Map/Alignment with deferred deletes | no exception (errors fixture) |
| 6 | Add tab 'Mine' (Sequence+PCR), close, open new window | PCR view created; new window tabs == ['Sequence · Map · Features'], no 'Mine'; reset_layout gives 1 tab |
| 7 | workspace.grouped_views() | groups Views/Annotations/Tools; Views == ['Alignment','Map','Protein','Sequence']; Annotations == Features, Primers, SEPARATOR, Detect Features, Detect ORFs; Tools sorted and contains Assembly, Digest, Local Align, PCR, Spacer Finder; no Display Filter/Quick Align/Primer Design; count == VIEWS+TAB_TOOLS |
| 8 | Open picker popup; click Map; click Assembly | 3 columns with headings, buttons match groups, 1 divider, padding '5px 16px 5px 16px', theme colours in style sheet, current marked; pane becomes Map and popup closes; Assembly opens its tab, pane stays Map |
| 9 | theme.ROLES light/dark; theme.is_dark stubbed | header bg darker than item bg in both; hex('group_item_bg') picks the dark/light entry |
| 10 | Pane set to Protein; set_tab_name 'My edit' then ''; old saved names Edit/Tab 3/Mine | tab name starts 'Protein'; manual name sticks and is dumped; empty name reverts; old layouts named ['Sequence','Features','Mine'] |

**Review notes:**

- L685: fixed QTest.qWait(600) for a 400 ms refresh timer (app/tools.py:32) leaves 200 ms margin; can flake on a loaded machine.
- L791-796: theme.is_dark is replaced by hand without try/finally or monkeypatch; a failing assert at L793/L795 leaks the stub into later tests.
- L719-723: the view-switching regression is a crash-only (smoke) check.

## `test_primer_pair_selection_and_circular_origin`

**Result:** ✅ passed in 0.49 s · [line 818](../../tests/test_gui.py#L818) · **Area:** GUI: primer-pair selection and circular sequences · **Checked against:** GUI state · **Strength:** 🟢 strong

Shift-clicking a primer pair selects the PCR product (also across the origin); the Circular button linearizes at a chosen base with confirmation and undo; the circular sequence view scrolls endlessly and treats selections across the origin as one stretch.

| # | Input | Expected output |
|---|---|---|
| 1 | Select 'Verif pBAD30_F' (4750..4769) and shift-select 'oAC1346' (complement(5670..5699)), both orders | selection (4749, 5699) |
| 2 | Shift-select another forward primer | current_primer == that primer |
| 3 | oAC736 (join(5931..5936,1..23)) + Verif pBAD30_R (complement(75..94)) | selection (5930, length+94), wraps |
| 4 | Circular toggle: cancel at 1001; confirm No (cuts Ampicilin 521..1381); Yes | stays circular; features_cut_at(1000) == ['Ampicilin']; sequence unchanged; then linear, sequence == seq[1000:]+seq[:1000], modified, Ampicilin has 2 spans |
| 5 | undo / redo / press again / undo | one undo restores circle and sequence; redo linear; re-press circular; undo linear |
| 6 | Click red origin bar at end of last row | loop on, copies >= 3, rows == copies*rpc; origin_at false 3 bases left |
| 7 | Scroll to 2*one+30, one-30, 3 turns line by line | value wraps to one+30, 2*one-30, stays in [one, 2*one) |
| 8 | Cursor n-2, Shift+Right x5 | (n-2, n+3), selected == seq[-2:]+seq[:3], status 'across the origin' |
| 9 | pos_at in second copy; select(10,10+3n); set_cursor(n); feature join(n-9..n,1..10) copy/paste at 1000; delete it | row rpc, pos n+3; clamped to (10,10+n); cursor n; feature selected (n-10,n+10); one pasted feature of 20 bp, sequence(1000,1020)==seq[-10:]+seq[:10]; delete gives seq[10:-10], one undo restores |
| 10 | set_circular(False); select(n-5,n+5) | 1 copy, rows == rpc, selection clamped to (n-5,n) |

**Review notes:**

- L861-865: comment says 'Same strand, or no shift: just the clicked primer' but the no-shift case is not exercised and the same-strand case checks only current_primer, not that the selection is that primer's span.

## `test_display_filters_translation_menus_history_and_new_file`

**Result:** ✅ passed in 1.11 s · [line 951](../../tests/test_gui.py#L951) · **Area:** GUI: display filters, translation, zoom, history panel, New file · **Checked against:** GUI state · **Strength:** 🟢 strong

Covers display filters and quick filters, type colours, the range box, CDS translation, letter zoom, selection-bracket dragging, name-vs-DNA Find, ORF links and motif search, the shortcuts dialog, the history panel, File > New and selection context-menu helpers.

| # | Input | Expected output |
|---|---|---|
| 1 | Filters panel: CDS off; primers off; enzymes toggled via menu | CDS hidden but still in features; drawn_primers == [] and no primer rows; menu and panel stay in sync |
| 2 | Quick filter buttons enzymes/primers/annotations | width <= sizeHint+4 < quarter of row; switches mirror panel/menu; annotations off/on restores hidden_kinds == {'CDS'} |
| 3 | FeatureDialog kind 'terminator'; recolor_by_type() | colour == type_color('terminator'); every feature coloured by type |
| 4 | Range box '101..200'; select_range('{len-9} - 10') | (100,200); (len-10, len+10) |
| 5 | eilR CDS (5119..5697) codons, hover text, context menu 'Show Translation' | first codon at 5119 'M', codon count == length//3, bases (5118,5121); 'eilR' and 'aa' in hover; menu has Edit…, Delete, Show Translation; toggling turns translation off |
| 6 | Ctrl+= twice, Ctrl+0; drag ']' bracket of (100,110) by 2 bases | size 13, top row within 2 lines; size 11; selection (100,112) |
| 7 | Add CDS 'cat' at 200..260; lookup('cat'), lookup('CAT'), find('cat',0) | (199,260,index); index; hit is the feature or bases CAT/ATG |
| 8 | OrfDialog: eilR row; motif = aa 2-5 of eilR | 'eilR (annotated)'; select links feature; hit 'eilR:' selects (start+3, start+15) |
| 9 | ShortcutsDialog; history after typing GAG at 100, NNN over 200..203, 2 backspaces at 300, feature 'mine'; go/drag marker | 'Ctrl+Shift+F' listed, >=5 tabs; rows 'Opened pBAD30_eilR.gb', 'Typed 3 bp at 101: GAG', '<seq[197:200]> > NNN at 201', 'Deleted 2 bp at 299…', 'Added feature mine', NOW; going back/forward restores exact sequences; button texts |
| 10 | File > New with '1 acgtac gtNN\n7 ACGT'; New again; reverse primer on (944,1004); group_as_gene; key filters 'gaATCG1x-Z?T', 'MK1L*j' | 'acgtacgtNNACGT', circular, untitled; second New opens new empty window that accepts typing; primer == revcomp of selection; gene name prefilled 'eilR', 'no CDS' for RBS only; 'gaATCGT'; 'MKL*' |

**Review notes:**

- L1103-1104: 'Find also looks for the bases' passes when find() returns the feature (hit[2]==named.index); Editor.find adds feature hits too (app/editor.py:1245-1247), so a broken DNA search would not be caught.
- L1163-1167: the marker drag is simulated by calling private panel._moved with fabricated args, not a real drag.
- L1123: dialog.tabs.count() >= 5 is loose.

## `test_start_page_sponsors_and_updates`

**Result:** ✅ passed in 0.16 s · [line 1208](../../tests/test_gui.py#L1208) · **Area:** GUI: start page (welcome, recent files, updates, sponsors) · **Checked against:** GUI state · **Strength:** 🟢 strong

The welcome page shows once a day, can be turned off and forced, opens recent files, reports update status in place, pins/removes recent files, shows release notes once, cycles tips, offers the demo, and shows sponsor cards/contact only when configured.

| # | Input | Expected output |
|---|---|---|
| 1 | SponsorsDialog with shipped sponsors.json | cards == []; welcome text has 'Antoine Champie', 'official download', 'github.com/OpenDNA-project/openDNA' |
| 2 | due/show_at_startup; last_shown='2000-01-01'; snooze on/off; never=True; force=True | due True, then False after showing; None when shown; True on a new day; False when snoozed; True when unsnoozed; None when never; forced dialog shown with snooze checked and hidden() still True |
| 3 | recent_files stubbed to [pX.gb]; click item | 1 item 'pX.gb…'; open_path called with the path; dialog accepted |
| 4 | No recent files; check_updates with UPDATES_ENABLED False, then stubbed answers (no update, 9.9.9) | list disabled; 'Not available yet'; 'Up to date'; '9.9.9' and 'install:' |
| 5 | add_recent a then b; pin a; remove a | order [b,a]; a first and bold; a gone from list and pins, file still exists |
| 6 | seen_version removed; then '0.0.1' with WHATS_NEW[VERSION]=['Pinned files'] | whats_new() None on first install; news shown with 'Pinned files', hidden on next open |
| 7 | next-tip link; quick actions with/without window; click demo | tip text changes, TIPS unique and >1; all quick buttons disabled without window; demo enabled, new_file disabled; open_untitled called with DEMO_NAME, accepted |
| 8 | load stubbed with Acme + 2 empty slots; CONTACT_EMAIL '' / 'opendna@example.org' | 3 cards, first has 'Acme', 'SPONSORED' label visible; no 'mailto' / 'mailto:opendna@example.org' in welcome_html |

**Review notes:**

- L1245: depends on the shipped app/sponsors.json having empty slots; adding a real sponsor breaks the test.
- L1248: assumes fresh welcome settings in the shared data dir (no 'never', not shown today); the cleanup at L1344 only runs if every earlier assert passed, so one failure cascades into later runs in the same session.

## `test_gene_groups_and_detect_common_features`

**Result:** ✅ passed in 0.94 s · [line 1328](../../tests/test_gui.py#L1328) · **Area:** GUI: gene groups and Detect Common Features · **Checked against:** GUI state · **Strength:** 🟢 strong

Gene groups (promoter+RBS+CDS) behave as one unit across table, map, edit, copy-paste, save/reopen, ungroup and delete; Detect Common Features finds library parts (also a codon-recoded sfGFP via its protein), adds them grouped, honours settings and user-imported parts.

| # | Input | Expected output |
|---|---|---|
| 1 | Promoter and RBS typed 60..20 and 15..8 bp upstream of the longest CDS; make_gene 'myGene' | parts ['promoter','RBS','CDS']; master span == min/max of parts, strand of CDS; features n+3; undo -> no groups, n+2 |
| 2 | Feature table and map; toggle collapse | only master row until opened, then parts right below; map shows master not parts; parts on map when collapse off |
| 3 | Move promoter 100 bp further out | gene span grows to the new end; undo restores |
| 4 | Copy-paste gene at 0 | 2 groups of 3 parts each |
| 5 | Save grouped.gb, reopen | 4 '/locus_tag="myGene"' lines; groups restored with same part kinds |
| 6 | ungroup; undo; delete_gene | no groups, n+2 features; then no groups, n-1 features |
| 7 | Synthetic FASTA (seed 11): 400 rnd + Lac promoter + 20 rnd + recoded sfGFP (>5% bases changed) + 30 rnd + T7 terminator + 300 rnd + revcomp(AmpR) + 400 rnd | hits ('AmpR',-1,DNA), ('Lac promoter',1), ('T7 terminator',1), ('sfGFP',1,protein) at 400+len(promoter)+20 with identity 1.0; none known |
| 8 | dialog.add() with grouping on; undo; grouping off | features added; gene 'sfGFP' with parts CDS, promoter, terminator; all hits known; one undo step; no groups when off |
| 9 | Settings dialog min_length = 31 | tab label matches tool; value saved 31; Lac promoter (30 bp) no longer a hit |
| 10 | Add misc_feature 'myTag' 101..160, save, import_files, delete feature, rerun detection | import count >= 1; 'myTag' detected and not known |

**Review notes:**

- L1493: import_files adds 'myTag' to the user parts library in the shared data dir and never removes it, so later tests in the session see an extra part.
- L1466: len(features) > n only; the number of added features is not pinned.

## `test_read_colours_local_align_and_map_reads`

**Result:** ✅ passed in 1.34 s · [line 1478](../../tests/test_gui.py#L1478) · **Area:** Read alignment colouring, Local Align, read mapping (bwa) · **Checked against:** hand-computed · **Strength:** 🟢 strong

Checks IGV-style read-pair colour categories, Local Align of a mutated reverse-complemented query, and the full read-mapping workflow (bwa mem, BAM track, variant calling, coverage, applying variants, cancel).

| # | Input | Expected output |
|---|---|---|
| 1 | classify() on hand-made pairs, limits (50,1000), mode 'Pair orientation & insert size' | →← insert 300 'normal' (both mates); 4000 'larger'; 40 'smaller'; →→ 'LL'; ←← 'RR'; ←→ 'RL'; mate on other contig 'elsewhere'; mate unmapped 'normal'; Strand mode reverse read 'reverse' |
| 2 | Query = revcomp of seq[1000:1300] with base 50 changed, 'GGG' inserted at 150, 5 bases deleted at 250 (via clipboard) | prefilled; top hit (1000,1300,-1), 1 mismatch, 3 inserted, 5 deleted; row '1,001–1,300' and INVERTED_TEXT; click selects (1000,1300) |
| 3 | Alignment view after hits; clear(); Pane('Genome') | query_hits count == dlg.hits, 'Query hits' track and a qhit (1000,1300); removed after clear; old name maps to 'Alignment' |
| 4 | 50 kb random ref (seed 3), sample with SNP at 20,000 and 5-bp deletion at 30,000; paired reads every 15 bp + one 4.1 kb pair; map_reads | BAM + .bai exist; alignment entry seqid 'refA'; run folder has ref.fa and 'bwa mem' command |
| 5 | View SNP ±150 | >=10 covering reads all carrying (snp,alt); only 'normal' kind; pileup count == covering; flag_level 'consensus'; condensed rows > 20 |
| 6 | Variant scan | ('snp',20000,'consensus') found; one del within 5 bp of 30000, ref length 5, consensus; exactly 2 variants; red mark at '20,001 ' |
| 7 | Click a read; view 4900..5300; color_by 'Strand' | selection == read span, info has 'MAPQ' and 'insert'; 1 'larger' read; kinds {'forward','reverse'} |
| 8 | annotate_variants / apply_variants (locked then unlocked) | 2 variation features; locked returns 0 and sequence unchanged; unlocked 2, length 49,995, base 20000 == alt, sequence == sample; 2 undos restore ref |
| 9 | Whole-sequence coverage; RunDialog('sleep 30') cancelled | mean of 10 bins in (12,14) (≈13.3 = 200 bases per 15 bp); run not ok, stopped within 20 s |

**Review notes:**

- L1571-1572: runtime pytest.skip when bwa/samtools/minimap2 are missing. CI (windows-latest, .github/workflows/ci.yml) and this machine have none, so the whole mapping half (L1573-1708: mapping, variants, coverage, apply, cancel) never runs and the test is reported SKIPPED even though the classify/Local Align sections ran (confirmed: 'SKIPPED tests/test_gui.py:1572').
- L1603, L1639, L1668, L1691, L1701: polling loops with 20-120 s timeouts (timing-dependent when they do run).

## `test_genome_is_responsive`

**Result:** ✅ passed in 13.74 s · [line 1679](../../tests/test_gui.py#L1679) · **Area:** GUI: performance on a 4.6 Mb genome · **Checked against:** GUI state · **Strength:** 🔴 weak

Opening, scrolling, typing in and matching primers on the E. coli BW25113 genome must work; the name claims it stays responsive.

| # | Input | Expected output |
|---|---|---|
| 1 | Open BW25113_2026.gb (4,631,469 bp); jump scroll bar to 1/3, 1/2, max with repaints | no exception (timings printed only) |
| 2 | Type 'A' at 1000 | length == 4,631,470; undo |
| 3 | db().add('genomeP', sequence(500000,500022)); wait up to 60 s | a genomeP site with spans [(500000,500022)] |

**Review notes:**

- Name promises responsiveness, but no time limit is asserted: load, scroll (L1722) and typing (L1729) timings are only printed.
- No docstring.
- L1736: polling loop up to 60 s for background primer matching.
- L1735: adds 'genomeP' to the shared primer DB and never removes it.
