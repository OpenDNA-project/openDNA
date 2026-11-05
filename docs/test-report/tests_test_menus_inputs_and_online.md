# `tests/test_menus_inputs_and_online.py`

[← all test files](README.md) · source: [`tests/test_menus_inputs_and_online.py`](../../tests/test_menus_inputs_and_online.py)

## `test_right_click_on_several_features_or_primers`

**Result:** ✅ passed in 2.92 s · [line 43](../../tests/test_menus_inputs_and_online.py#L43) · **Area:** GUI: multi-row feature/primer menu · **Checked against:** hand-computed · **Strength:** 🟢 strong

Right-click menu actions apply to all selected rows: hide/show/delete features (one undo step) and hide/show/delete primers, leaving unselected rows untouched.

| # | Input | Expected output |
|---|---|---|
| 1 | 3 features a,b,c; select a,b; features_actions menu | 'Delete 2 Features' and 'Hide All' present, 'Show All on Sequence' absent |
| 2 | trigger 'Hide All' | hidden == {a,b}; c untouched |
| 3 | 'Show All on Sequence' then 'Delete 2 Features' | no hidden; features == [c]; one undo restores a,b,c |
| 4 | primers multi1,multi2: Hide All / Show All / Delete 2 | hidden set then cleared; db().find returns None for both |

**Review notes:**

- Uses the shared global db() (lines 70,82) with no private-db fixture; the module sets a fresh GENOMEEDIT_DATA_DIR at import so cross-file leakage is limited, but names could collide with other tests in this module.

## `test_about_dialog_and_view_theme`

**Result:** ✅ passed in 2.60 s · [line 87](../../tests/test_menus_inputs_and_online.py#L87) · **Area:** GUI: theme switching · **Checked against:** GUI state · **Strength:** 🟢 strong

Help > About shows the licence notices; View > Theme (System, Light, Dark) applies and remembers the choice, starts at System, and every view redraws in the new colours.

| # | Input | Expected output |
|---|---|---|
| 1 | theme menu actions | texts == {System, Light, Dark}; all have non-null icons |
| 2 | ensure a MapView, grab each (caches image), act 'Dark' | styleHints().setColorScheme got Dark; saved()=='dark'; dark action checked; every map.cache is None |
| 3 | act 'System' | saved()=='system'; system action checked |
| 4 | win.about() with QMessageBox.about stubbed | text contains 'GNU Affero General Public License', 'NO WARRANTY', 'Copyright' |
| 5 | Theme actions System/Light/Dark (styleHints stubbed) | exactly these 3, all checkable, System checked by default |
| 6 | Trigger Dark, Light, System | setColorScheme called with Dark/Light/Unknown; theme.saved() == 'dark'/'light'/'system' |
| 7 | theme.apply('nonsense') | returns 'system' |

**Review notes:**

- L833: styleHints is stubbed, so only the call to setColorScheme is checked, not that colours/palette actually change; L841 repaint is smoke.
- L837: 'System checked by default' depends on no earlier test having left another theme in the shared data dir's theme.ini.

**Fixed by the audit of November 2026:**

- Removed in the second pass of the audit (fewer tests to review, no biological check lost): tests/test_gui.py::test_about_dialog_and_light_dark_themes was merged into this test.

## `test_database_status_in_the_databases_menu`

**Result:** ✅ passed in 0.07 s · [line 127](../../tests/test_menus_inputs_and_online.py#L127) · **Area:** GUI: Databases > Status · **Checked against:** hand-computed · **Strength:** 🟢 strong

The Databases>Status submenu lists each known database, each status check returns rows and a valid level, the bundled features library reads OK, and _parts_status classifies usable/empty/non-ACGTN parts and a missing file.

| # | Input | Expected output |
|---|---|---|
| 1 | status submenu actions | texts == [label for label,_ in db_status.DATABASES] |
| 2 | each check(): st and db_status.show dialog | st.rows non-empty; st.level in {ok,warn,error}; dialog windowTitle==st.title |
| 3 | features_library() | level 'ok'; rows[0][0]=='Sources' |
| 4 | _parts_status('Test', tsv with 1 ok + empty + ACGXX, 'error','missing') | Parts starts '1 usable of 3'; 'No sequence' level warn, mentions 'empty'; 'Not only A,C,G,T,N' mentions 'weird'; st.level 'warn' |
| 5 | _parts_status('Gone', missing file,...) | level 'error' |

**Review notes:**

- The menu-labels assertion (line 116) compares the menu against the same db_status.DATABASES constant used to build it — it verifies wiring but is a mild self-reference, not an independent oracle.

## `test_open_recent_lists_the_start_page_files`

**Result:** ✅ passed in 0.06 s · [line 152](../../tests/test_menus_inputs_and_online.py#L152) · **Area:** GUI: File > Open Recent · **Checked against:** hand-computed · **Strength:** 🟢 strong

The Open Recent menu shows 'No recent files' when empty, lists pinned files first with bold font and '&&'-escaped names, sets tooltips to resolved paths, and is attached to the File menu.

| # | Input | Expected output |
|---|---|---|
| 1 | clear recent/pin keys, fill menu | actions == ['No recent files'] |
| 2 | add a.gb and b&c.gb, pin a.gb, emit aboutToShow | texts == ['&1 a.gb','&2 b&&c.gb']; action[0] bold; action[1].toolTip == resolved b path |
| 3 | recent_menu action | present in the &File menu's actions |

## `test_tracks_that_do_not_match_are_loaded_with_a_warning`

**Result:** ✅ passed in 0.78 s · [line 173](../../tests/test_menus_inputs_and_online.py#L173) · **Area:** GUI: track mismatch detection · **Checked against:** hand-computed · **Strength:** 🟢 strong

track_check.mismatches flags a wrong sequence name, items past the end, and VCF REF bases that disagree with the sequence; a mismatching track is still loaded via the menu but with a warning, a matching one loads silently.

| # | Input | Expected output |
|---|---|---|
| 1 | matching GFF3 (name 't') | mismatches == [] |
| 2 | GFF3 named chrZ with a feature at 480..900 on a shorter seq | a problem naming 'chrZ' and 't'; a 'beyond the end' problem mentioning 900 |
| 3 | VCF two variants at pos 10: one matching REF=seq[9], one wrong | exactly 1 problem; '1 of 2 variant(s)' and 'the file says {wrong}' |
| 4 | mismatches(ed, None, 't', 999, 'x.bam') | first problem starts "Its reference 't' is 999 bp" |
| 5 | load_track(other) via menu; then load_track(good) | entry loaded into tracks; warning mentions 'beyond the end' and 'loaded anyway'; good track -> no warning |

## `test_a_control_that_sends_data_online_carries_the_globe`

**Result:** ✅ passed in 0.00 s · [line 202](../../tests/test_menus_inputs_and_online.py#L202) · **Area:** GUI: online-indicator (privacy) · **Checked against:** hand-computed · **Strength:** 🟡 adequate

Controls that send data to the internet (the NCBI Search button) carry the globe icon and a tooltip naming the endpoint, and the dialog no longer has an email/contact field.

| # | Input | Expected output |
|---|---|---|
| 1 | NcbiDialog; online.marked(d) | d.find in marked; icon not null; tooltip contains online.NOTE and 'eutils.ncbi.nlm.nih.gov' |
| 2 | check for email field | no 'email' attr; no QLineEdit with 'contact' placeholder |

## `test_a_bogus_primer_list_imports_nothing_wrong`

**Result:** ✅ passed in 1.90 s · [line 216](../../tests/test_menus_inputs_and_online.py#L216) · **Area:** Input validation: primer import · **Checked against:** hand-computed · **Strength:** 🟢 strong

Only valid primers are imported from TSV/FASTA; each rejected line is reported with its line number and reason; invalid primers never reach the database via import, add, or add_many; binary files are refused.

| # | Input | Expected output |
|---|---|---|
| 1 | TSV with good1, 'words'(not dna), short(ACGT), one-column, empty-name, uracil, comment, good2(comma) | good1 & good2 in db; rejected line numbers == [3,4,5,6,7]; reasons mention 'not DNA','at least 8','No tab/comma/semicolon','name','U'; report says '5 lines not imported' and 'line 3'; words/short/uracil not in db |
| 2 | FASTA '>fa1 desc\nACGTACGT\nACGTAC\n>fa2\nACGT-XYZ-ACGT\n>\nACGTACGTACGT' | out == [('fa1','ACGTACGTACGTAC')]; rejected lines [4,6]; 'not DNA' in reason |
| 3 | binary tsv with NUL byte | raises ValueError match 'not a text file' |
| 4 | db().add('tab\there',...), add('ok_name','ACGTXX'); add_many([fine, bad tab name, 500-A toolong]); update to 'NOT DNA' | add returns 'invalid' both; add_many added==1, invalid names ['bad\tone','toolong']; update raises ValueError |

## `test_parts_and_text_fields_refuse_bad_characters`

**Result:** ✅ passed in 0.07 s · [line 258](../../tests/test_menus_inputs_and_online.py#L258) · **Area:** Input validation: parts/names/keys/qualifiers · **Checked against:** hand-computed · **Strength:** 🟢 strong

Parts library, feature names/qualifiers, gene names, tab names and the NCBI API key all reject tabs/control chars/non-ACGT or over-long text with a named reason, and a trimmed valid key is sent.

| # | Input | Expected output |
|---|---|---|
| 1 | add_parts([p1 ok, 'bad\tname', iupac 'ATGRYAAATAA']) | added==1; rejected names ['bad\tname','iupac']; reasons 'control character' and 'R Y' |
| 2 | validate.text_problem variants; api_key_problem | 'fine name' ok; tab->'control character'; 30 chars->'too long'; '0123abcdef' ok; 'abc-def!'->'letters and digits' |
| 3 | NcbiDialog key 'abc def;' then search; then '  abc123  ' | status 'letters and digits', no request; valid key trimmed to 'abc123' and sent |
| 4 | FeatureDialog name 'my\tfeature'; qualifiers 'note=a\tb'; then valid | messages mention 'control character' / 'tab or control character'; no feature until valid -> label 'my feature' |
| 5 | group_as_gene with name 'gene\x07x'; rename_tab with bad input | returns None; tab text unchanged; 'control character' in messages |

## `test_project_files_are_open_dna_projects`

**Result:** ✅ passed in 0.00 s · [line 309](../../tests/test_menus_inputs_and_online.py#L309) · **Area:** Project file type constants · **Checked against:** hand-computed · **Strength:** 🟡 adequate

The project file extension is .odp (Open DNA Project) and appears in the open-file filter, with no stale 'oge' naming.

| # | Input | Expected output |
|---|---|---|
| 1 | exports.PROJECT_EXT and main.OPEN_FILTER | PROJECT_EXT=='.odp'; '*.odp' in filter; 'oge' not in filter (lowercased); 'Open DNA Project' in filter |
