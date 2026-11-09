"""Genetic codes: every embedded NCBI table against Biopython (used here only, never by the app), alternative start codons
in the ORF finder, a CDS read with its own /transl_table and /codon_start, the Tools > Settings choice and its persistence,
and the /translation refresh (Feature menu): what it updates, what it leaves alone, one undo step, never on autosave."""
import os
import sys
import tempfile
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ["OPENDNA_DATA_DIR"] = tempfile.mkdtemp(prefix="opendna_test_")
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "app"))

import pytest
from Bio.Data import CodonTable
from PySide6.QtWidgets import QApplication

import genetic_code
import tool_settings
from editor import Editor, reverse_complement
from opendna_core import Document, find_protein, translate
from sequence_view import cds_codons, codons
from test_cas9 import close, open_doc

app = QApplication.instance() or QApplication([])

# Tables left out on purpose: the meaning of a codon there depends on its place in the gene (see core/src/code.rs).
LEFT_OUT = {27, 28, 31}
IDS = sorted(CodonTable.unambiguous_dna_by_id.keys() - LEFT_OUT)


@pytest.fixture
def code_setting():
    """Leave the genetic code settings as found (the settings folder is shared by all test modules)."""
    saved = {k: tool_settings.value(genetic_code.TOOL, k) for k in ("table", "alt_starts")}
    yield
    for k, v in saved.items():
        tool_settings.set_value(genetic_code.TOOL, k, v)


@pytest.fixture
def keep_setting():
    saved = genetic_code.keep_translation()
    yield
    genetic_code.set_keep_translation(saved)


def choose(table, alt=False):
    label = next(x for x in genetic_code.LABELS if x.split()[0] == str(table))
    tool_settings.set_value(genetic_code.TOOL, "table", label)
    tool_settings.set_value(genetic_code.TOOL, "alt_starts", alt)


def cds_doc(seq, *features):
    doc = Document.new_sequence("t", seq, False)
    doc.add_features(list(features))
    return doc


def cds(location, label="g", **quals):
    return ("CDS", location, [("label", label)] + [(k, v) for k, v in quals.items()])


# ---- the tables --------------------------------------------------------------------------------------------------------------
def test_the_embedded_tables_match_biopython_and_differ_where_they_should():
    """The embedded tables are the expected ones; every codon of every table on both strands and the start / stop
    sets match Biopython; tables differ where they should; Rust and Python agree on the data; the citation."""
    # --- the embedded tables are the expected ones
    assert [i for i in sorted(genetic_code.TABLES) if not genetic_code.is_custom(i)] == IDS   # (the recoded codes, 100+, are not NCBI's)
    assert IDS[0] == 1 and 11 in IDS and 4 in IDS and 2 in IDS
    # --- every codon of every table matches Biopython on both strands, and the start and stop sets
    for table in IDS:
        ref = CodonTable.unambiguous_dna_by_id[table]
        bases = "ACGT"
        all_codons = [a + b + c for a in bases for b in bases for c in bases]
        expected = {c: "*" if c in ref.stop_codons else ref.forward_table[c] for c in all_codons}
        assert len(expected) == 64
        for c in all_codons:
            assert translate(c, table) == expected[c], (table, c)
            assert genetic_code.codon_map(table)[c] == expected[c], (table, c)
        # The reverse strand: a CDS on the bottom strand reads the reverse complement of the top strand.
        top = "".join(reverse_complement(c) for c in all_codons)
        doc = Document.new_sequence("rev", top, False)
        got = codons(doc, [(0, len(top))], -1, 0, table)
        # Reading backwards, the codons come out in the reverse order of how they were laid down.
        assert [aa for _, aa, _, _ in got] == [expected[c] for c in reversed(all_codons)], table
        # And the forward strand in one go.
        assert translate("".join(all_codons), table) == "".join(expected[c] for c in all_codons), table
        t = genetic_code.TABLES[table]
        assert sorted(t.stop_codons) == sorted(ref.stop_codons), table
        assert sorted(t.start_codons) == sorted(ref.start_codons), table
    # --- tables differ where they should
    assert translate("TGA", 1) == "*" and translate("TGA", 4) == "W"        # mycoplasma: TGA is Trp
    assert translate("AGA", 1) == "R" and translate("AGA", 2) == "*"        # vertebrate mitochondria: AGA/AGG stop
    assert translate("ATA", 1) == "I" and translate("ATA", 2) == "M"
    assert translate("CTG", 12) == "S" and translate("CTG", 1) == "L"       # alternative yeast nuclear
    assert translate("TGA", 25) == "G"                                       # Gracilibacteria
    # Table 11 has the standard amino acids and stops, but more start codons.
    assert genetic_code.TABLES[11].codons == genetic_code.TABLES[1].codons
    assert set(genetic_code.TABLES[1].start_codons) < set(genetic_code.TABLES[11].start_codons)
    # An unknown table falls back to the standard code rather than failing.
    assert translate("TGA", 99) == "*"
    # --- Rust and Python agree on the data
    from opendna_core import genetic_codes
    assert [t[0] for t in genetic_codes()] == IDS
    assert all(len(t[2]) == len(t[3]) == 64 for t in genetic_codes())
    # --- the citation for the genetic codes
    import citations
    text = " ".join(" ".join(str(x) for x in ref) for refs in citations.PANELS.values() for ref in refs)
    assert "Elzanowski" in text and "wprintgc" in text


# ---- ORFs --------------------------------------------------------------------------------------------------------------------
def test_orfs_start_and_stop_follow_the_table():
    """ORFs start at alternative codons only when allowed; stops follow the table; protein motif search uses the table."""
    # --- alternative start codons
    gtg = Document.new_sequence("g", "CCGTGAAAAAATAACC", False)         # GTG AAA AAA TAA
    assert gtg.orfs(4) == []
    assert gtg.orfs(4, 11) == [], "alternative starts are off by default"
    assert gtg.orfs(4, 11, True) == [(2, 14, 1)]
    assert gtg.orfs(4, 1, True) == [], "table 1 has no GTG start"
    ttg = Document.new_sequence("t", "CCTTGAAAAAATAACC", False)
    assert ttg.orfs(4, 1, True) == [(2, 14, 1)] and ttg.orfs(4, 1) == []
    # The reverse strand too.
    assert Document.new_sequence("r", reverse_complement("CCGTGAAAAAATAACC"), False).orfs(4, 11, True) == [(2, 14, -1)]
    # --- stops follow the table
    seq = Document.new_sequence("m", "ATGTGAAAATAA", False)             # ATG TGA AAA TAA
    assert seq.orfs(4) == [] and seq.orfs(4, 4) == [(0, 12, 1)]          # table 4: TGA is not a stop
    assert Document.new_sequence("v", "ATGAGAAAATAA", False).orfs(1, 2) == [(0, 6, 1)]   # table 2: AGA stops
    # --- protein motif search uses the table
    seq = "AAATGAAAAAA"   # ... TGA AAA: W K in table 4, stop K in table 1 (frame +3 starts at base 2)
    assert find_protein(seq, "WK", 0, 100, 4)
    assert not find_protein(seq, "WK", 0, 100, 1)


# ---- a CDS: /transl_table, /codon_start, /transl_except -------------------------------------------------------------------
def test_a_cds_is_read_with_its_own_table_frame_and_exceptions(code_setting):
    """/transl_table beats the setting, /codon_start, start codon shown as M, /transl_except (and the unreproducible
    one flagged), and the tooltip says so."""
    # --- /transl_table regardless of the setting
    doc = cds_doc("ATGTGAGTGTAA", cds("1..12", transl_table="4"), cds("1..12", label="plain"))
    ed = Editor(doc, Path("t.gb"))
    own, plain = ed.features
    for setting in (1, 11):
        choose(setting)
        assert "".join(aa for _, aa, _, _ in cds_codons(doc, own)) == "MWV*"
    choose(4)
    assert "".join(aa for _, aa, _, _ in cds_codons(doc, plain)) == "MWV*", "no /transl_table: the setting applies"
    choose(1)
    assert "".join(aa for _, aa, _, _ in cds_codons(doc, plain)) == "M*V*"
    # --- /codon_start and the start codon as methionine
    doc = cds_doc("ATGTGAGTGTAA" + "TTGAAATAA",
                  cds("1..12", codon_start="2"), cds("13..21", label="ttg"), cds("<13..21", label="partial", codon_start="1"))
    two, ttg, partial = Editor(doc, Path("t.gb")).features
    assert "".join(aa for _, aa, _, _ in cds_codons(doc, two)) == "CEC"
    # A start codon of the table at the start of a complete CDS is shown (and exported) as M, like GenBank's /translation.
    assert "".join(aa for _, aa, _, _ in cds_codons(doc, ttg)) == "MK*"
    assert "".join(aa for _, aa, _, _ in cds_codons(doc, partial)) == "LK*"
    # --- /transl_except is applied and an unreproducible one is flagged
    seq = "ATGTGATAA"
    sec = cds_doc(seq, cds("1..9", transl_except="(pos:4..6,aa:Sec)"))
    assert "".join(aa for _, aa, _, _ in cds_codons(sec, Editor(sec, Path("t.gb")).features[0])) == "MU*"
    odd = cds_doc(seq, cds("1..9", transl_except="(pos:4..5,aa:TERM)"))
    f = Editor(odd, Path("t.gb")).features[0]
    assert "not reproducible" in genetic_code.cds_note(odd, f)
    assert "".join(aa for _, aa, _, _ in cds_codons(odd, f)) == "M**", "drawn plainly, and flagged"
    # --- the tooltip names the table and the unreproducible translation
    from hover import describe
    doc = cds_doc("ATGTGAGTGTAA", cds("1..12", transl_table="4"), cds("1..12", label="x", transl_except="(pos:4..5,aa:TERM)"))
    ed = Editor(doc, Path("t.gb"))
    assert "Table 4" in describe(ed, ed.features[0].index)
    assert "translation not reproducible" in describe(ed, ed.features[1].index)


# ---- the setting ----------------------------------------------------------------------------------------------------------
def test_the_setting_persists_shows_in_the_dialog_and_drives_views_and_exports(code_setting):
    """Persists and defaults to the standard code; the settings dialog shows stops and starts of the selected table;
    the protein view, exports and the sequence view's CDS lane use the chosen table."""
    # --- persists and defaults to the standard code
    assert genetic_code.LABELS[0] == "1 – Standard"
    tool_settings.set_value(genetic_code.TOOL, "table", "x")
    assert genetic_code.table_id() == 1, "an unknown value falls back to the default"
    choose(11, alt=True)
    tool_settings._cache.clear()   # read it back from settings.ini, as a new session would
    assert genetic_code.table_id() == 11 and genetic_code.alt_starts() is True
    assert tool_settings.value(genetic_code.TOOL, "table").startswith("11 – Bacterial")
    # --- the settings dialog shows stops and starts of the selected table
    choose(1)   # back to the standard code, alternative starts off
    dialog = tool_settings.SettingsDialog(None, genetic_code.TOOL)
    p, combo = dialog.widgets[(genetic_code.TOOL, "table")]
    assert combo.count() == len(genetic_code.TABLES) >= len(IDS) + 2 and combo.itemText(0) == "1 – Standard"
    page = dialog.tabs.currentWidget()
    from PySide6.QtWidgets import QLabel
    lines = [w.text() for w in page.findChildren(QLabel) if w.text().startswith("Stop codons")]
    assert lines == ["Stop codons: TAA, TAG, TGA  ·  Start codons: TTG, CTG, ATG"]
    combo.setCurrentText(next(x for x in genetic_code.LABELS if x.startswith("2 ")))
    lines = [w.text() for w in page.findChildren(QLabel) if w.text().startswith("Stop codons")]
    assert lines == ["Stop codons: TAA, TAG, AGA, AGG  ·  Start codons: ATT, ATC, ATA, ATG, GTG"]
    assert dialog.widgets[(genetic_code.TOOL, "alt_starts")][1].isChecked() is False   # off by default
    # --- views and exports use the chosen table
    import exports
    from protein import ProteinView
    doc = cds_doc("CCATGTGAAAATAACC", cds("3..14"))
    ed = Editor(doc, Path("t.gb"))
    f = ed.features[0]
    choose(1)
    assert exports.protein_of(ed, f) == "M*K", "an internal stop is kept: only the final one is stripped"
    assert exports.translate("ATGTGA") == "M*"
    pane = ProteinView(ed)
    ed.select(2, 14)
    pane.refresh(True)
    assert "Table" not in pane.title.text()
    choose(4)
    tool_settings.notifier.changed.emit(genetic_code.TOOL)
    assert exports.protein_of(ed, f) == "MWK"
    assert exports.translate("ATGTGA") == "MW"
    assert "Table 4" in pane.title.text()
    # A selection that is not a CDS is read with the setting, and says so.
    ed.select(2, 14)
    ed.current_feature = None
    pane.refresh(True)
    assert pane.src["kind"] == "selection" and "Table 4" in pane.title.text()
    assert "".join(aa for _, aa, _, _ in pane.src["cods"]) == "MWK*"
    # The sequence view's CDS lane reads the same way.
    win = open_doc(Document.new_sequence("w", "CCATGTGAAAATAACC", False))
    try:
        win.editor.add_features([cds("3..14")])
        assert "".join(aa for _, aa, _, _ in win.sequence_view._codons(win.editor.features[0])) == "MWK*"
        choose(1)
        assert "".join(aa for _, aa, _, _ in win.sequence_view._codons(win.editor.features[0])) == "M*K*"
    finally:
        close(win)


# ---- /translation refresh ------------------------------------------------------------------------------------------------
STALE = "MQ"


def refresh_doc():
    """#0 stale, #1 codon_start 2 stale, #2 pseudo, #3 transl_except, #4 no /translation, #5 table 4, #6 partial without frame."""
    seq = "ATGAAATAA" + "GATGAAATAA" + "ATGAAATAA" + "ATGTGATAA" + "ATGAAATAA" + "ATGTGAAAATAA" + "ATGAAATAA"
    return cds_doc(seq,
                   cds("1..9", translation=STALE),
                   cds("10..19", label="cs2", codon_start="2", translation=STALE),
                   cds("20..28", label="pseudo", pseudo="", translation=STALE),
                   cds("29..37", label="exc", transl_except="(pos:32..34,aa:Sec)", translation="MQ"),
                   cds("38..46", label="none"),
                   cds("47..58", label="t4", transl_table="4", translation=STALE),
                   cds("<59..67", label="partial", translation=STALE))


def trans(ed, label):
    return dict(next(f for f in ed.features if f.label == label).qualifiers).get("translation")


def test_refresh_updates_stale_translations_and_leaves_the_rest_alone(code_setting):
    """/translation refresh: what it updates, what it leaves alone, one undo step; the setting's table for a CDS without
    its own; it matches what the views draw; the report text."""
    # --- updates stale translations and leaves the rest alone
    choose(1)
    ed = Editor(refresh_doc(), Path("t.gb"))
    before = {f.label: f.qualifiers for f in ed.features}
    depth = len(ed.doc.history()[0])
    report = ed.refresh_translations()
    status = {ed.features[i].label: s for i, s, _ in report}
    assert status == {"g": "updated", "cs2": "updated", "pseudo": "pseudo", "exc": "updated", "t4": "updated", "partial": "partial"}
    assert trans(ed, "g") == "MK"
    assert trans(ed, "cs2") == "MK", "/codon_start=2 is honoured"
    assert trans(ed, "t4") == "MWK", "/transl_table=4 wins over the setting"
    assert trans(ed, "exc") == "MU", "the selenocysteine /transl_except is applied"
    assert trans(ed, "pseudo") == STALE and trans(ed, "partial") == STALE
    assert trans(ed, "none") is None, "a CDS without /translation never gets one"
    # Only /translation values changed.
    for f in ed.features:
        assert [kv for kv in f.qualifiers if kv[0] != "translation"] == [kv for kv in before[f.label] if kv[0] != "translation"]
    # One undo step; a second run does nothing and adds no step.
    assert len(ed.doc.history()[0]) == depth + 1
    assert all(s == "unchanged" or s in ("pseudo", "partial") for _, s, _ in ed.refresh_translations())
    assert len(ed.doc.history()[0]) == depth + 1
    ed.undo()
    assert {f.label: f.qualifiers for f in ed.features} == before
    # --- uses the settings table for a CDS without its own
    choose(4)
    doc = cds_doc("ATGTGAAAATAA", cds("1..12", translation="M"))
    ed = Editor(doc, Path("t.gb"))
    assert [s for _, s, _ in ed.refresh_translations()] == ["updated"]
    assert trans(ed, "g") == "MWK"
    choose(1)
    assert [s for _, s, why in ed.refresh_translations()] == ["other"], "the standard code: a stop inside the CDS, left alone"
    assert trans(ed, "g") == "MWK"
    # --- matches the translation drawn in the views: Rust (refresh) and Python (drawing, FASTA export) follow the same
    # plan: selenocysteine, codon_start, tables, strands
    import exports
    choose(11)
    seq = "ATGTGATAA"
    top = reverse_complement(seq)
    doc = cds_doc(top + "GATGAAATAA" + "TTGAAATAA",
                  cds("complement(1..9)", label="rev", transl_except="(pos:complement(4..6),aa:Sec)", translation="X"),
                  cds("10..19", label="cs", codon_start="2", translation="X"),
                  cds("20..28", label="ttg", translation="X"))
    ed = Editor(doc, Path("t.gb"))
    ed.refresh_translations()
    for f in ed.features:
        assert trans(ed, f.label) == exports.protein_of(ed, f), f.label
    assert trans(ed, "rev") == "MU" and trans(ed, "ttg") == "MK"
    # --- the sort of report text
    report = [(0, "updated", ""), (1, "updated", ""), (2, "updated", ""), (3, "pseudo", ""), (4, "unchanged", "")]
    assert genetic_code.summarize(report) == "Updated /translation of 3 CDSs, skipped 1 pseudo"
    assert genetic_code.summarize([(0, "updated", "")]) == "Updated /translation of 1 CDS"


# ---- Feature menu, save and autosave -----------------------------------------------------------------------------------
def test_features_menu_toggle_save_autosave_and_update_now(tmp_path, keep_setting, code_setting):
    """The Feature menu toggle is persistent and on by default; save refreshes translations but autosave never does; save
    leaves them alone when the option is off; Update now acts on the selected CDS, else all, as one undo step."""
    # --- the Feature menu toggle is persistent and on by default
    genetic_code.set_keep_translation(True)
    win = open_doc(Document.new_sequence("w", "ACGT" * 10, False))
    try:
        assert win.keep_translation_action.isChecked() and win.keep_translation_action.isCheckable()
        assert "Keep /translation" in win.keep_translation_action.text().replace("&", "")
        win.keep_translation_action.trigger()
        assert genetic_code.keep_translation() is False
    finally:
        close(win)
    again = open_doc(Document.new_sequence("w", "ACGT" * 10, False))   # a new window reads the saved choice
    try:
        assert not again.keep_translation_action.isChecked()
    finally:
        close(again)
    # --- save refreshes translations but autosave never does
    choose(1)
    genetic_code.set_keep_translation(True)
    win = open_doc(refresh_doc())
    try:
        ed = win.editor
        path = tmp_path / "t.gb"
        ed.path, ed.untitled = path, False
        ed.doc.modified = True
        ed.version += 1   # (the document was filled directly, so autosave would see "nothing new" otherwise)
        version, depth = ed.version, len(ed.doc.history()[0])
        win.autosave()
        # Autosave is a pure write: the document, its version and its undo history are untouched, the copy has the old text.
        assert (ed.version, len(ed.doc.history()[0])) == (version, depth)
        assert trans(ed, "g") == STALE
        recovery = Path(str(path) + ".autosave.gb")
        assert recovery.exists() and "MQ" in recovery.read_text()
        # An explicit save refreshes, writes the new text, and leaves the file saved.
        assert win._write(path)
        assert trans(ed, "g") == "MK"
        assert not ed.doc.modified
        reread = Editor(Document(str(path)), path)   # what is on disk now
        assert trans(reread, "g") == "MK" and trans(reread, "pseudo") == STALE
        assert "Updated /translation of 4 CDSs, skipped 1 partial, 1 pseudo" in win.statusBar().currentMessage()
        # Saving again: nothing to refresh, no new undo step, still not modified (no loop).
        depth_after = len(ed.doc.history()[0])
        assert win._write(path)
        assert len(ed.doc.history()[0]) == depth_after and not ed.doc.modified
        # Undo takes the refresh back, which makes the document differ from the file again.
        ed.undo()
        assert trans(ed, "g") == STALE and ed.doc.modified
    finally:
        close(win)
    # --- save leaves translations alone when the option is off
    genetic_code.set_keep_translation(False)
    win = open_doc(refresh_doc())
    try:
        assert win._write(tmp_path / "off.gb")
        assert trans(win.editor, "g") == STALE
    finally:
        close(win)
    # --- Update now acts on the selected CDS, else all, as one undo step
    genetic_code.set_keep_translation(True)
    choose(1)
    win = open_doc(refresh_doc())
    try:
        ed = win.editor
        depth = len(ed.doc.history()[0])
        ed.select_feature(0)
        win.update_translations()
        assert trans(ed, "g") == "MK" and trans(ed, "cs2") == STALE, "only the selected CDS"
        ed.current_feature = None
        win.update_translations()
        assert trans(ed, "cs2") == "MK" and trans(ed, "t4") == "MWK"
        assert len(ed.doc.history()[0]) == depth + 2
        assert "Updated /translation of 3 CDSs, skipped 1 partial, 1 pseudo" in win.statusBar().currentMessage()
        ed.undo()
        assert trans(ed, "cs2") == STALE and trans(ed, "g") == "MK"
    finally:
        close(win)


def test_the_toolbar_quick_switch_follows_the_setting_and_can_be_hidden(code_setting):
    """The quick switch shows the table in use, changes the same setting as Tools > Settings, and hides with its checkbox."""
    import code_switch
    saved_show = tool_settings.value(genetic_code.TOOL, "show_switch")
    try:
        choose(1)
        heard = []
        tool_settings.notifier.changed.connect(heard.append)
        switch = code_switch.CodeSwitch()
        assert switch.text() == "Code 1" and len(switch.actions_by_label) == len(genetic_code.LABELS)
        label = next(x for x in genetic_code.LABELS if x.split()[0] == "11")
        switch.choose(label)
        assert genetic_code.table_id() == 11 and switch.text() == "Code 11" and genetic_code.TOOL in heard
        assert switch.actions_by_label[label].isChecked()
        choose(2)   # changed in the Settings: the switch follows
        tool_settings.notifier.changed.emit(genetic_code.TOOL)
        assert switch.text() == "Code 2"
        assert code_switch.shown() is True   # shown by default
        tool_settings.set_value(genetic_code.TOOL, "show_switch", False)
        assert code_switch.shown() is False
        tool_settings.notifier.changed.disconnect(heard.append)
    finally:
        tool_settings.set_value(genetic_code.TOOL, "show_switch", saved_show)


def test_settings_pages_follow_the_new_tab_picker(code_setting):
    """Tools > Settings lists one page per view or tool in the picker's groups, then General; a page without options says so,
    only the opened page's group is unfolded, and Primers holds the primer options and the check limits."""
    import primer_panel
    import workspace
    from PySide6.QtWidgets import QLabel
    layout = tool_settings.page_layout()
    assert [g for g, _ in layout] == ["Views", "Annotations", "Tools", "General"]
    picker = [(g, [n for n in names if n]) for g, names in workspace.grouped_views()]
    assert layout[:3] == picker
    assert layout[3][1][:4] == ["Genetic code", "Display", "History", "Updates"]
    dialog = tool_settings.SettingsDialog(None, "Digest")
    tree = dialog.tabs.tree
    assert dialog.tabs.currentPage() == "Digest"
    assert [tree.topLevelItem(i).isExpanded() for i in range(tree.topLevelItemCount())] == [False, False, True, False]
    assert any(w.text() == tool_settings.NO_OPTIONS for w in dialog.tabs.currentWidget().findChildren(QLabel))
    dialog.tabs.show_page("Genetic code")   # General unfolds, Tools folds
    assert [tree.topLevelItem(i).isExpanded() for i in range(tree.topLevelItemCount())] == [False, False, False, True]
    primers = tool_settings.SettingsDialog(None, "Primers")
    page = primers.tabs.currentWidget()
    assert isinstance(primers.custom["Primers"], primer_panel.OptionsForm) and primers.page_for(primer_panel.CHECKS) is page
    assert (primer_panel.CHECKS, "tm_low") in primers.widgets
    assert tool_settings.SettingsDialog(None, "Alignments").tabs.currentPage() == "Alignment"   # a tool key opens its view's page
    assert tool_settings.SettingsDialog(None, "Spacer Finder").tabs.currentPage() == "Spacer Finder"
