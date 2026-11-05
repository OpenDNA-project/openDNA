# `tests/test_assembly.py`

[← all test files](README.md) · source: [`tests/test_assembly.py`](../../tests/test_assembly.py)

## `test_fragment_and_gibson_engine`

**Result:** ✅ passed in 0.02 s · [line 92](../../tests/test_assembly.py#L92) · **Area:** Assembly: Fragment + Gibson engine · **Checked against:** independent implementation · **Strength:** 🟢 strong

A Frag is the exact PCR product with correct tails and flips correctly; Gibson reports ends that don't fit then fixes them with new primers; detects existing overlap, non-unique homology and short homology; and every junction sits exactly where it claims in the product.

| # | Input | Expected output |
|---|---|---|
| 1 | Frag('x',doc,100,322, fwd='GGGGCC'+t[100:122], rev='TTTAAA'+rc(t[300:322])) | seq=='GGGGCC'+t[100:322]+rc('TTTAAA'), tails==(6,6), core_seq==t[100:322]; flip() reverse-complements, tails preserved, F/R swap roles; flip back restores |
| 2 | frags() = A(300),B(400),C(250); analyse(f,0,True) then auto_primers then analyse | initially all joints 'none', not ok, 'Circle closes'=='bad', 'no homology' problem; auto_primers==3; then ok, joints [(A,B),(B,C),(C,A)], each 20<=len<=40 and joint_tm>=52, 6 rows, product len 950; rotation_of design_gibson product; 'Circle closes'=='ok' |
| 3 | existing shared dna(30) overlap; duplicated shared; short shared[:16] | existing overlap kind, product 530; duplicated -> 'Homology unique (a → b)'=='bad', not ok; short[:16] -> existing overlap but 'Homology (a → b)'=='warn' |
| 4 | four designs (gibson circ/lin, sidewinder, golden gate) junction positions | joined_bases(d,j)==j.seq for every non-none joint; golden-gate j.a/j.b == (pos,pos+4); sidewinder j.a[0]==j.b[0]==pos+TOEHOLD |

**Review notes:**

- Loop at 143-145 guards the position check with 'if j.kind != "none"'; if a design had only 'none' joints the check would be vacuous, but preceding asserts establish the gibson designs are ok so low risk.

## `test_golden_gate_engine_and_overhang_checks`

**Result:** ✅ passed in 0.03 s · [line 166](../../tests/test_assembly.py#L166) · **Area:** Assembly: Golden Gate + overhangs · **Checked against:** hand-computed · **Strength:** 🟢 strong

Golden Gate turns a vector+inserts into one circle, reads the real overhang each end carries, rejects/flags bad setups (wrong overhang pairs, missing/inner sites, palindromes, GC-biased, cross-ligating), and the overhangs module matches the Potapov 2018 validated set.

| # | Input | Expected output |
|---|---|---|
| 1 | vector_frag + ins(700) + term(220), analyse(f,1,True,'BsaI'); auto_primers==2; flip insert; bad vector setups | initially all 'none'/not ok; auto_primers==2 (vector cut not amplified); then ok, 4 rows, all joints len 4, 3 distinct, rotation of design_golden_gate, unique=='ok', 'BsaI sites left in product'=='ok'; flipped insert not ok until re-primed; ValueError 'needs exactly 2' / 'point the same way' |
| 2 | ring(['GGAG','TCGC','CAGT']) and broken variants | joints seq==ov[1:]+ov[:1], product 3*(4+250), 'validated'=='ok', no bad/warn; mismatched overhang -> joint 'none' naming TCGC/TCGA; missing site -> 'no BsaI site'; inner site -> 'BsaI sites in g0'=='bad' |
| 3 | chips() over poor overhang sets (duplicate, palindrome GATC, GC-only, reverse-complement pair) | ('unique','bad'); ('not palindromic','bad') for GATC not GGTT; ('orthogonal','warn')+('validated','info') for near set; ('GC mix','warn') for AAAA; ('orthogonal','bad') for CAGT/ACTG; clean set has no bad/warn |
| 4 | OH.VALIDATED, OH.assess, OH.pair_risk, OH.junction_checks, A.pick_overhangs | VALIDATED len 25 unique all len4; assess(VALIDATED) no bad/warn; POOR_EXAMPLE flagged warn; pair_risk('GGAG','CTCC')=='bad' (rc), ('GGAG','GGAA')=='warn'; junction_checks all 'ok'; pick_overhangs prefers validated, never repeats |

## `test_sidewinder_product_views_and_citations`

**Result:** ✅ passed in 0.03 s · [line 242](../../tests/test_assembly.py#L242) · **Area:** Assembly: Sidewinder, product, views, citations · **Checked against:** independent implementation · **Strength:** 🟢 strong

Sidewinder builds the simulated product keeping each method's features, enforces short-fragment limits, the junction-view HTML shows sketch/reason/tails, and the citations list names all required DOIs/sources.

| # | Input | Expected output |
|---|---|---|
| 1 | frags a(100,promoter)/b(120)/c(90), analyse(f,2,False); simulate; long frags(); f[:1] | 2 joints, problems [None,None], ok; product doc len 310 with 3 misc_features; long not ok and 'Oligos at most 165 nt'=='bad'; ValueError 'at least two' |
| 2 | Gibson and Golden Gate products keep features; piece_doc | for method 0,1: ok, product keeps 'pro' label, circular, len==product; piece_doc length==len(seqs[0]) |
| 3 | junction_html for plain/gibson/golden gate | plain shows 'no homology' and >=3 ✗; primed html shows 'reverse-complemented','as given', note prefix, one <table> per joint, '[A]','bottom', each j.seq and its complement, tail background count >= 60; golden-gate shows 'overhang' and cut ends |
| 4 | citations_html() | contains 14 needles incl DOIs 10.1021/acssynbio.8b00333, 10.1038/nmeth.1318, 'Sidewinder','SantaLucia','PySide6','PyO3', etc. |

## `test_panel_boxes_check_design_perform_and_golden_gate_vector`

**Result:** ✅ passed in 0.12 s · [line 334](../../tests/test_assembly.py#L334) · **Area:** GUI: Assembly panel flow · **Checked against:** GUI state · **Strength:** 🟢 strong

The Assembly panel starts with gray boxes tracking the fragment count, clicking a gray box asks for a source, and it auto-checks/designs/performs, including Golden Gate vector and enzyme handling.

| # | Input | Expected output |
|---|---|---|
| 1 | empty AssemblyPanel; count 5 then 1; render strip | 3 methods, 3 slots all None, 'still empty', perform disabled; slots resize to 5 then 1; count min 1 max MAX_FRAGMENTS; gray boxes paint |
| 2 | click gray box / filled box / its ▾ | empty box click -> source_asked [1], current 1; filled box click selects+shows file (asked []); click ▾ rect -> asked [2] |
| 3 | filled_panel(0): add_homology, perform, flip/reorder | 'do not fit' + perform disabled + 3 'bad' cards; after add_homology design.ok, 'checks out', perform enabled, 6 rows; perform -> opened len 1800 circular; flip breaks (not ok); rotation (_move_to 2,0) still ok; real reorder (_move_to 1,0) not ok |
| 4 | filled_panel(1) Golden Gate vector/enzyme | toggling vector on plain vector -> 'needs exactly 2'; new 2-site source keeps vector role; add_homology ok+perform enabled; _show_product shows 3 overhang features, product seq matches; switch to BsmBI -> design None, '0 BsmBI site(s)', perform disabled |

## `test_panel_sequence_view_selection_pcr_and_saving`

**Result:** ✅ passed in 0.19 s · [line 409](../../tests/test_assembly.py#L409) · **Area:** GUI: panel sequence view, selection, PCR, saving · **Checked against:** GUI state · **Strength:** 🟢 strong

The panel's embedded sequence view is read-only, a selected stretch becomes a fragment, two primers make a PCR fragment with the correct tails, and fragments/primers are saved correctly.

| # | Input | Expected output |
|---|---|---|
| 1 | filled_panel(0); feed key events A/Backspace/Delete/Ctrl-X/Ctrl-V and a right-click to the view | doc.sequence and version unchanged (nothing typed/cut/pasted/deleted); arrow key moves cursor to 9 |
| 2 | select(10,30) then use_selection; re-select same; select(5,5) | 'Selected 11–30 (20 bp)', by_region checked; slot (10,30) len20 name 'insert 11-30'; file shown; same stretch keeps name; coming back shows selection; empty selection -> 'Select a stretch' |
| 3 | DB primers F1='GGTCTCAAATG'+t[200:224], R1='TTT'+rc(t[600:624]); select both primers; use_selection | fwd strand 1 rev -1, fwd.tail>0; forward/reverse count 1, by_primers enabled; 'F1 + R1' chosen; slot (fw.start,rv.end), fwd[0]=='F1'; seq starts with tail+core, tails[0]==len(tail); ends with rc of rev anneal; 'PCR with F1 + R1'; file-only primer 'mine' also found |
| 4 | filled_panel(0) add_homology; save_fragments to dir; save_primers | save_fragments==3; files insert/term/vector.gb; reloaded insert length/seq match slot; 'added' in status |

## `test_assembly_is_a_tab_opened_from_the_menu`

**Result:** ✅ passed in 0.23 s · [line 497](../../tests/test_assembly.py#L497) · **Area:** GUI: Assembly tab lifecycle · **Checked against:** GUI state · **Strength:** 🟢 strong

Assembly is a standalone tab (not a pane dropdown) opened from the menu, survives closing/reopening and workspace changes, and keeps its fragment state and count across file opens.

| # | Input | Expected output |
|---|---|---|
| 1 | 'Assembly' not in workspace.VIEWS; MainWindow.show_assembly() with no file | assembly starts None; after show, tool_host tab 0 'Assembly' is current; close_tab keeps panel (parent None); count set 2 survives reopen |
| 2 | open a saved p.gb file | tab moves into workspace next to others, 2 tabs, 'Assembly' at idx1 current; layout flags [True,False]; dump() length 1 (tool tab not saved); count still 2; switching tabs keeps sequence_view; duplicate_tab keeps the only pane tab; recent_files includes path |
| 3 | _open_path(0,path); count 1; add_homology; perform; open another file; close+reopen assembly | slot0 name 'p' length 500; one fragment design.ok; perform works; tab+content survive new workspace; close then show_assembly restores with slot0 'p' |

## `test_every_method_puts_its_junction_info_on_a_card_and_the_map_shows_them`

**Result:** ✅ passed in 0.09 s · [line 548](../../tests/test_assembly.py#L548) · **Area:** GUI: junction cards and map · **Checked against:** GUI state · **Strength:** 🟢 strong

Each assembly method renders its junction info on a card (right chips, closing flags, both-strand sequence), a palindromic vector overhang is flagged red, and the map draws cards alternating with fragment boxes and opens detail/report on click.

| # | Input | Expected output |
|---|---|---|
| 1 | junction_cards for Gibson/Golden Gate/Sidewinder and a failing Gibson | Gibson 3 cards closing [F,F,T], both-strand bases, 'bp homology'/'°C'/'GC' in title, chips [length,Tm,GC,unique,hairpin,distinct]; GG chips [ends fit,unique,not palindromic,orthogonal,GC mix,validated]; Sidewinder 2 cards 'toehold'/'barcode', chips [GC,unique,distinct]; failing -> all 'bad' 'do not fit' |
| 2 | vector with palindromic overhang AATT | not ok, card level 'bad', ('not palindromic','bad') in jchecks[0] |
| 3 | filled_panel(0) map rendering and clicks | no tabs/junction_view/table widgets, no QTabWidget/QTableWidget, no 'Gibson et al'; after add_homology 3 card_rects width>=CARD_MIN, boxes width<=130, boxes/cards alternate left-to-right; click card -> 'Junction 2' dialog; show_report -> 'Assembly Report' |

## `test_designed_primers_are_temporary_until_saved_and_explain_themselves_on_hover`

**Result:** ✅ passed in 0.19 s · [line 602](../../tests/test_assembly.py#L602) · **Area:** GUI: temporary designed primers + hover · **Checked against:** GUI state · **Strength:** 🟢 strong

Primers the designer creates are temporary (not in DB) and shown on the correct strand with a tail, their hover tooltip explains them, saving makes them permanent at the same place, and flipping re-creates temporary oligos on the opposite strands.

| # | Input | Expected output |
|---|---|---|
| 1 | filled_panel(0) _show_source(1); add_homology; inspect sites | no 'insert-' sites before; after add_homology fwd/rev sites both temporary, strands (1,-1), fwd.tail==len(fwd[1])-fwd[2]>0; primer_info hidden |
| 2 | hover the fwd primer in the view; describe() | hover_target returns site.key; tooltip has primer name, 'temporary primer','Tm','hairpin','Save primers', and tail lowercase |
| 3 | save_primers; then _flip_at(1)+add_homology | after save same sites non-temporary, tooltip no longer 'temporary' but still 'primer'; flipped: fwd/rev strands (-1,1) |

## `test_actions_menu_pointer_product_box_and_sidewinder_panel`

**Result:** ✅ passed in 0.17 s · [line 644](../../tests/test_assembly.py#L644) · **Area:** GUI: panel actions, menus, product box, Sidewinder · **Checked against:** GUI state · **Strength:** 🟢 strong

Action buttons are a right-hand column with Perform last/highlighted, the fragment right-click menu offers move/flip/vector/clear, the pointer is a hand over boxes/cards, the product box previews the assembled product with colored bars, and the Sidewinder panel lists oligos.

| # | Input | Expected output |
|---|---|---|
| 1 | filled_panel(0) layout of buttons/status/fragment controls | buttons {homology,save,primers,perform}, no move/rc/clear buttons; homology text 'Design primers'; status+buttons right of map, ordered top-to-bottom perform last; perform+use objectName 'primary', highlight styleSheet; fragment controls right of seqstack |
| 2 | filled_panel(1) right-click menu actions and triggers | menu_asked [1] on right-click; texts [Source,Move left,Move right,Reverse complement,Destination vector...,Clear]; Source submenu 3 entries; Move left reorders; first frag can't move left; Reverse complement flips (seq==rc); vector checkable toggles; Clear empties; gray box only 'Source'; Gibson has no vector entry |
| 3 | hover over box/card/empty; product box click before/after add_homology | hand cursor over box and card, arrow over empty; product_text None until fits; click before add_homology -> not shown; after -> design.ok, text ends 'circular', click shows product, view is product_view, seq==product, 3 colored fragment bars + 3 homology features matching joint lengths; _flip_at redraws/gone; _pick(0) returns to source with controls enabled |
| 4 | Sidewinder panel over seeds 0..99 until design.ok | design.ok, 2 cards, not circular, homology button disabled; primer_info shows 's1-coding' and 's1-barcode'; _show_product shows toehold features; perform -> opened len 300 |

**Review notes:**

- Lines 755-763: loop over 100 seeds until a valid Sidewinder design is found; if none succeeded the final assert would fail (not vacuous) but the test is mildly seed-dependent/fragile.
