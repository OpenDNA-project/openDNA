"""Right-click on several features or primers, the theme button (and every view redrawing at once), Databases > Database Status,
File > Open Recent, tracks that do not match the sequence, the globe on controls that send data over the internet, and the checks
on what is imported or typed (primer lists, parts, names, the NCBI key)."""
import os
import sys
import tempfile
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ["OPENDNA_DATA_DIR"] = tempfile.mkdtemp(prefix="opendna_test_")
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "app"))

import pytest
from PySide6.QtCore import Qt
from PySide6.QtWidgets import QApplication, QInputDialog, QLineEdit, QMenu, QMessageBox

from opendna_core import Document, Track

import validate
from editor import feature_key
from primer_db import db
from test_cas9 import close, dna, open_doc

app = QApplication.instance() or QApplication([])


@pytest.fixture
def win(monkeypatch):
    w = open_doc(Document.new_sequence("t", dna(500), False))
    w.messages = []
    for kind in ("warning", "information", "critical"):
        monkeypatch.setattr(QMessageBox, kind, staticmethod(lambda parent, title, text, *a, **k: w.messages.append(text)))
    yield w
    close(w)


def act(menu, text):
    next(a for a in menu.actions() if a.text() == text).trigger()


# ---- several selected rows ----------------------------------------------------------------------------------------------------

def test_right_click_on_several_features_or_primers(win):
    """Every single-row action applies to all the selected rows: hide / show / delete (one undo step) for features, hide / show /
    delete for primers."""
    from feature_table import features_actions, primers_actions
    ed = win.editor
    ed.add_features([("misc_feature", "10..40", [("label", "a")]), ("CDS", "50..109", [("label", "b")]),
                     ("misc_feature", "200..240", [("label", "c")])])
    chosen = [f.index for f in ed.features if f.label in ("a", "b")]
    menu = QMenu()
    features_actions(menu, None, ed, chosen)
    texts = [a.text() for a in menu.actions()]
    assert "Delete 2 Features" in texts and "Hide All" in texts and "Show All on Sequence" not in texts
    act(menu, "Hide All")
    assert sorted(f.label for f in ed.features if feature_key(f) in ed.hidden_features) == ["a", "b"], "c is untouched"
    menu = QMenu()
    features_actions(menu, None, ed, chosen)
    act(menu, "Show All on Sequence")
    assert not ed.hidden_features
    menu = QMenu()
    features_actions(menu, None, ed, chosen)
    act(menu, "Delete 2 Features")
    assert [f.label for f in ed.features] == ["c"]
    ed.undo()
    assert sorted(f.label for f in ed.features) == ["a", "b", "c"], "one undo brings both back"

    # --- primers ---
    for name in ("multi1", "multi2"):
        db().add(name, "ACGTACGTACGTACGTAC")
    menu = QMenu()
    primers_actions(menu, ed, {"multi1", "multi2"})
    act(menu, "Hide All in This File")
    assert {"multi1", "multi2"} <= ed.hidden_primers
    menu = QMenu()
    primers_actions(menu, ed, {"multi1", "multi2"})
    act(menu, "Show All in This File")
    assert not {"multi1", "multi2"} & ed.hidden_primers
    menu = QMenu()
    primers_actions(menu, ed, {"multi1", "multi2"})
    act(menu, "Delete 2 from Primer Database")
    assert db().find("multi1") is None and db().find("multi2") is None


# ---- theme -------------------------------------------------------------------------------------------------------------------

def test_about_dialog_and_view_theme(win, monkeypatch):
    """Help > About shows the licence notices; View > Theme (System, Light, Dark) applies and remembers the choice, starts at
    System, and every view redraws in the new colours."""
    import theme
    from map_view import MapView
    # --- the About dialog shows the licence notices ---
    shown = []
    monkeypatch.setattr(QMessageBox, "about", lambda parent, title, text: shown.append(text))
    win.about()
    assert "GNU Affero General Public License" in shown[0] and "NO WARRANTY" in shown[0] and "Copyright" in shown[0]

    # --- the theme menu ---
    theme.apply(theme.saved(), remember=False)   # links the real colour-scheme signal first (the stub below has none)
    applied = []
    stub = type("Hints", (), {"setColorScheme": lambda self, s: applied.append(s)})()   # the offscreen platform has no scheme
    monkeypatch.setattr(theme.QGuiApplication, "styleHints", staticmethod(lambda: stub))
    menu = win.theme_menu   # View > Theme (no toolbar button: the theme is chosen once)
    assert {a.text() for a in menu.actions()} == {"System", "Light", "Dark"}
    assert all(a.isCheckable() for a in menu.actions()) and win.theme_actions["system"].isChecked(), "System unless chosen otherwise"
    assert all(not a.icon().isNull() for a in menu.actions()), "a sun, a dark sun, a sun with an A"
    win.workspace.ensure(MapView)
    maps = win.findChildren(MapView)
    assert maps
    for m in maps:
        m.grab()   # draws, so the map keeps a cached image
        assert m.cache is not None
    act(menu, "Dark")
    assert applied[-1] == Qt.ColorScheme.Dark and theme.saved() == "dark" and win.theme_actions["dark"].isChecked()
    assert all(m.cache is None for m in maps), "the map drops its image drawn in the old colours at once"
    win.sequence_view.viewport().repaint()
    act(menu, "Light")
    assert applied[-1] == Qt.ColorScheme.Light and theme.saved() == "light" and win.theme_actions["light"].isChecked()
    act(menu, "System")
    assert applied[-1] == Qt.ColorScheme.Unknown and theme.saved() == "system" and win.theme_actions["system"].isChecked(), \
        "System follows the operating system"
    assert theme.apply("nonsense") == "system"


# ---- database status ---------------------------------------------------------------------------------------------------------

def test_database_status_in_the_databases_menu(win, tmp_path):
    import db_status
    databases = next(a.menu() for a in win.menuBar().actions() if a.text() == "Data&bases")
    status = next(a.menu() for a in databases.actions() if a.menu() and "Status" in a.text())
    assert [a.text() for a in status.actions()] == [label for label, _ in db_status.DATABASES]
    for _, check in db_status.DATABASES:
        st = check()
        assert st.rows and st.level in ("ok", "warn", "error"), st.title
        dialog = db_status.show(win, check)
        assert dialog.windowTitle() == st.title
        dialog.close()
    library = db_status.features_library()
    assert library.level == "ok" and library.rows[0][0] == "Sources", "the bundled library is readable"
    bad = tmp_path / "parts.tsv"
    bad.write_text("#title: test\nname\ttype\tsequence\nok\tpromoter\tACGTACGT\nempty\tRBS\t\nweird\tCDS\tACGXX\n", encoding="utf-8")
    st = db_status._parts_status("Test", bad, "error", "missing")
    rows = {label: (value, level) for label, value, level in st.rows}
    assert rows["Parts"][0].startswith("1 usable of 3")
    assert rows["No sequence"][1] == "warn" and "empty" in rows["No sequence"][0]
    assert "weird" in rows["Not only A, C, G, T, N"][0] and st.level == "warn"
    assert db_status._parts_status("Gone", tmp_path / "none.tsv", "error", "missing").level == "error"


# ---- File > Open Recent ------------------------------------------------------------------------------------------------------

def test_open_recent_lists_the_start_page_files(win, tmp_path):
    import sourcefiles
    sourcefiles._settings().setValue(sourcefiles.KEY, [])
    sourcefiles._settings().setValue(sourcefiles.PIN_KEY, [])
    win._fill_recent_menu()
    assert [a.text() for a in win.recent_menu.actions()] == ["No recent files"]
    a, b = tmp_path / "a.gb", tmp_path / "b&c.gb"
    for p in (a, b):
        Document.new_sequence(p.stem, dna(100), False).save_genbank(str(p))
        sourcefiles.add_recent(p)
    sourcefiles.set_pinned(str(a.resolve()), True)
    win.recent_menu.aboutToShow.emit()
    actions = win.recent_menu.actions()
    assert [x.text() for x in actions] == ["&1 a.gb", "&2 b&&c.gb"], "pinned first, & shown as itself"
    assert actions[0].font().bold() and actions[1].toolTip() == str(b.resolve())
    file_menu = next(x.menu() for x in win.menuBar().actions() if x.text() == "&File")
    assert win.recent_menu.menuAction() in file_menu.actions()


# ---- tracks that do not match ------------------------------------------------------------------------------------------------

def test_tracks_that_do_not_match_are_loaded_with_a_warning(win, tmp_path):
    import track_check
    ed = win.editor
    seq = ed.doc.sequence(0, ed.length)
    good = tmp_path / "good.gff3"
    good.write_text("##gff-version 3\nt\tx\tgene\t5\t50\t.\t+\t.\tName=a\n")
    assert track_check.mismatches(ed, Track(str(good)), "t", path=good) == []
    other = tmp_path / "other.gff3"
    other.write_text("##gff-version 3\nchrZ\tx\tgene\t5\t50\t.\t+\t.\tName=a\nchrZ\tx\tgene\t480\t900\t.\t+\t.\tName=b\n")
    problems = track_check.mismatches(ed, Track(str(other)), "chrZ", path=other)
    assert any("'chrZ'" in p and "'t'" in p for p in problems), "another sequence name"
    assert any("beyond the end" in p and "900" in p for p in problems), "items past the end"
    vcf = tmp_path / "v.vcf"
    wrong = "A" if seq[9] != "A" else "C"
    vcf.write_text(f"##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\nt\t10\t.\t{seq[9]}\tG\nt\t10\t.\t{wrong}\tG\n")
    problems = track_check.mismatches(ed, Track(str(vcf)), "t", path=vcf)
    assert len(problems) == 1 and "1 of 2 variant(s)" in problems[0] and f"the file says {wrong}" in problems[0]
    assert track_check.mismatches(ed, None, "t", 999, "x.bam")[0].startswith("Its reference 't' is 999 bp")
    # --- through the menu: a warning, and the track is loaded anyway ---
    win.messages.clear()
    entry = win.load_track(str(other))
    assert entry is not None and entry in ed.tracks, "loaded: the user decides"
    assert win.messages and "beyond the end" in win.messages[-1] and "loaded anyway" in win.messages[-1]
    win.messages.clear()
    assert win.load_track(str(good)) is not None and not win.messages, "a matching track: no warning"


# ---- the internet symbol -----------------------------------------------------------------------------------------------------

def test_a_control_that_sends_data_online_carries_the_globe():
    import online
    from ncbi_dialog import NcbiDialog
    d = NcbiDialog(None)
    marked = online.marked(d)
    assert d.find in marked and not d.find.icon().isNull(), "Search sends the words to NCBI: globe"
    assert online.NOTE in d.find.toolTip() and "eutils.ncbi.nlm.nih.gov" in d.find.toolTip()
    assert not hasattr(d, "email") and not [e for e in d.findChildren(QLineEdit) if "contact" in e.placeholderText()], \
        "no e-mail field any more"
    d.close()


# ---- what is typed or imported -----------------------------------------------------------------------------------------------

def test_a_bogus_primer_list_imports_nothing_wrong(tmp_path):
    """Only real primers get in; every refused line is named with its number and the reason; nothing invalid reaches the database
    whichever way it comes (a list, add, add_many)."""
    from primer_panel import import_primers, parse_primer_list
    bogus = tmp_path / "bogus.tsv"
    bogus.write_text("name\tsequence\n"
                     "good1\tACGTACGTACGTACGTACGT\n"
                     "words\thello world, this is not dna\n"
                     "short\tACGT\n"
                     "just one column with no separator\n"
                     "\tACGTACGTACGTACGT\n"
                     "uracil\tACGUACGUACGUACGU\n"
                     "\n# a comment\n"
                     "good2,acgtacgtacgtACGTAC\n", encoding="utf-8")
    rejected = []
    report = import_primers(str(bogus), None, rejected)
    assert db().find("good1") and db().find("good2")
    assert [n for n, _, _ in rejected] == [3, 4, 5, 6, 7], rejected
    reasons = dict((n, why) for n, _, why in rejected)
    assert "not DNA" in reasons[3] and "at least 8" in reasons[4] and "No tab, comma or semicolon" in reasons[5]
    assert "name" in reasons[6].lower() and "U" in reasons[7]
    assert "5 lines not imported" in report and "line 3" in report
    for name in ("words", "short", "uracil"):
        assert db().find(name) is None
    # --- FASTA ---
    rejected = []
    out = parse_primer_list(">fa1 desc\nACGTACGT\nACGTAC\n>fa2\nACGT-XYZ-ACGT\n>\nACGTACGTACGT\n", rejected)
    assert out == [("fa1", "ACGTACGTACGTAC")]
    assert [n for n, _, _ in rejected] == [4, 6] and "not DNA" in rejected[0][2]
    # --- binary ---
    binary = tmp_path / "x.tsv"
    binary.write_bytes(b"name\tsequence\nab\x00cd\tACGTACGTACGT\n")
    with pytest.raises(ValueError, match="not a text file"):
        import_primers(str(binary))
    # --- the database itself ---
    assert db().add("tab\there", "ACGTACGTACGT") == "invalid" and db().add("ok_name", "ACGTXX") == "invalid"
    report = db().add_many([("fine_one", "ACGTACGTACGT"), ("bad\tone", "ACGTACGTACGT"), ("toolong", "A" * 500)])
    assert report["added"] == 1 and [n for n, _ in report["invalid"]] == ["bad\tone", "toolong"]
    with pytest.raises(ValueError):
        db().update("fine_one", "fine_one", "NOT DNA AT ALL", "")


def test_parts_and_text_fields_refuse_bad_characters(win, monkeypatch):
    from parts_library import add_parts
    from feature_dialog import FeatureDialog
    from feature_table import group_as_gene
    # --- parts library: a tab would break the file, other letters cannot be searched ---
    rejected = []
    added = add_parts([("p1", "promoter", "ACGTACGTTTGACA", "f.gb"), ("bad\tname", "CDS", "ATGAAATAA", "f.gb"),
                       ("iupac", "CDS", "ATGRYAAATAA", "f.gb")], rejected)
    assert added == 1 and [n for n, _ in rejected] == ["bad\tname", "iupac"]
    assert "control character" in rejected[0][1] and "R Y" in rejected[1][1]
    # --- shared checks ---
    assert validate.text_problem("fine name", "name", 20) == ""
    assert "control character" in validate.text_problem("a\tb", "name", 20) and "too long" in validate.text_problem("x" * 30, "name", 20)
    assert validate.api_key_problem("0123abcdef") == "" and "letters and digits" in validate.api_key_problem("abc-def!")
    # --- the NCBI key: named at once, nothing sent ---
    from ncbi_dialog import NcbiDialog
    d = NcbiDialog(None)
    sent = []
    monkeypatch.setattr(d, "_run", lambda *a: sent.append(a))
    d.key.setText("abc def;")
    assert "letters and digits" in d.status.text()
    d.query.setText("pUC19")
    d.search()
    assert not sent, "a bad key: no request"
    d.key.setText("  abc123  ")
    d.search()
    assert sent and d.key.text() == "abc123", "trimmed, then sent"
    d.close()
    # --- feature dialog: a pasted tab in the name or a qualifier value ---
    ed = win.editor
    ed.select(10, 40)
    f = FeatureDialog(win, ed)
    f.name.setText("my\tfeature")
    f.accept()
    assert "control character" in win.messages[-1] and not ed.features
    f.name.setText("my feature")
    f.qualifiers.setPlainText("note=a\tb")
    f.accept()
    assert "tab or control character" in win.messages[-1] and not ed.features
    f.qualifiers.setPlainText("note=fine")
    f.accept()
    assert [x.label for x in ed.features] == ["my feature"]
    # --- gene name and tab name ---
    monkeypatch.setattr(QInputDialog, "getText", staticmethod(lambda *a, **k: ("gene\x07x", True)))
    assert group_as_gene(win, ed, [ed.features[0].index]) is None and "control character" in win.messages[-1]
    tabs = win.workspace
    before = tabs.tabText(0)
    tabs.rename_tab(0)
    assert tabs.tabText(0) == before and "control character" in win.messages[-1]


def test_project_files_are_open_dna_projects():
    import exports
    import main as M
    assert exports.PROJECT_EXT == ".odp"
    assert "*.odp" in M.OPEN_FILTER and "oge" not in M.OPEN_FILTER.lower() and "Open DNA Project" in M.OPEN_FILTER
