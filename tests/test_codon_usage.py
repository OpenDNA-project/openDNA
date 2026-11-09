"""Codon usage tables (app/codon_usage.py) and the genetic codes of id 100+ (app/genetic_code.py): CSV in and out, your own tables on
disk, the codon adaptation index with its rare codons, the editor, the recoded codes (removed codons read as X) and the toolbar
switch. openDNA ships no codon usage table (no licence allows it), so every table here is made up."""
import os
import sys
import tempfile
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ["OPENDNA_DATA_DIR"] = tempfile.mkdtemp(prefix="opendna_test_")
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "app"))

import pytest
from PySide6.QtWidgets import QApplication

import codon_usage as cu
import genetic_code

app = QApplication.instance() or QApplication([])


def flat(name="flat", **changes):
    """A table where every codon is used 15.625 per thousand, except the codons in `changes`."""
    freq = {c: 15.625 for c in cu.CODONS}
    freq.update(changes)
    return cu._table(name, freq)


@pytest.fixture
def clean():
    """No table of yours; the table in use left as found."""
    for p in cu.folder().glob("*.json"):
        p.unlink()
    saved = cu.active_name()
    yield
    cu.set_active(saved)
    for p in cu.folder().glob("*.json"):
        p.unlink()


def test_no_codon_usage_table_ships_with_the_app(clean):
    assert cu.all_tables() == {} and cu.active() is None
    import checkup
    assert checkup.codon_notes("ATG" + "CTA" * 20 + "TAA") == []   # no table: no note, no error


def test_csv_round_trip_and_the_files_that_are_refused(tmp_path):
    t = flat(AAA=40.0)
    t.name, t.source = "My organism", "a paper"
    cu.write_csv(t, tmp_path / "t.csv")
    back = cu.read_csv(tmp_path / "t.csv")
    assert back.name == "My organism" and back.source == "a paper" and back.per_thousand["AAA"] == 40.0 and len(back.per_thousand) == 64
    pairs = tmp_path / "pairs.txt"   # a plain list of codon, frequency (U accepted for T, no header)
    pairs.write_text("\n".join(f"{c.replace('T', 'U')}\t{i}" for i, c in enumerate(cu.CODONS)), encoding="utf-8")
    assert cu.read_csv(pairs).per_thousand["TTT"] == 0 and cu.read_csv(pairs).per_thousand["GGG"] == 63
    short = tmp_path / "short.csv"
    short.write_text("codon,aa,per_thousand\nTTT,F,10\n", encoding="utf-8")
    with pytest.raises(ValueError, match="64 codons"):
        cu.read_csv(short)
    bad = tmp_path / "bad.csv"
    bad.write_text("codon,aa,per_thousand\nTTT,F,10\nXYZ,F,10\n", encoding="utf-8")
    with pytest.raises(ValueError, match="not a codon"):
        cu.read_csv(bad)


def test_your_tables_are_kept_in_the_data_folder(clean):
    mine = cu.duplicate(flat("Original"), "Recoded host")
    cu.save_user(mine)
    assert (cu.folder() / "Recoded_host.json").is_file() and cu.folder().parent == Path(os.environ["OPENDNA_DATA_DIR"])
    assert list(cu.all_tables()) == ["Recoded host"] and cu.all_tables()["Recoded host"].note == "Copy of Original."
    with pytest.raises(ValueError, match="needs a name"):
        cu.save_user(cu.duplicate(mine, "  "))
    (cu.folder() / "broken.json").write_text("{not json", encoding="utf-8")
    assert list(cu.user_tables()) == ["Recoded host"]   # a stray file is ignored
    cu.set_active("Recoded host")
    assert cu.active().name == "Recoded host"
    cu.delete_user("Recoded host")
    assert cu.all_tables() == {} and cu.active_name() == "" and cu.active() is None   # the table in use went


def test_a_table_that_cannot_be_written_is_reported_not_lost(clean, monkeypatch):
    import fileio

    def locked(tmp, target, *a, **k):
        raise PermissionError(13, "being used by another process")
    monkeypatch.setattr(fileio, "replace", locked)
    with pytest.raises(ValueError, match="could not be saved.*being used by another process"):
        cu.save_user(flat("Mine"))
    assert cu.all_tables() == {}
    aas, starts = genetic_code.recoded(11, ["AGA"])
    with pytest.raises(ValueError, match="could not be saved"):
        genetic_code.save_code("A code of mine", aas, starts)
    assert not any(t.name == "A code of mine" for t in genetic_code.TABLES.values())


def test_cai_and_rare_codons_follow_the_table():
    code = genetic_code.codon_map(1)
    table = flat(CTG=50.0, CTA=1.0)   # Leu: CTG is the best codon, CTA is used 1/50 as much
    assert cu.cai("CTGCTGCTG", table, code) == pytest.approx(1.0)
    assert cu.cai("CTACTA", table, code) == pytest.approx(1 / 50)   # w = 1 / 50 for each
    both = cu.cai("CTGCTA", table, code)
    assert 1 / 50 < both < 1 and both == pytest.approx((1 / 50) ** 0.5)   # the geometric mean
    assert cu.cai("ATGTGGTAA", table, code) is None   # Met, Trp and the stop have no synonym: nothing to score
    assert cu.rare_codons("ATGCTGCTAGGT", table, code) == [(3, "CTA")]
    assert cu.cai("CTA", flat(CTA=0.0, CTG=10.0), code) == pytest.approx(0.01)   # a codon never used counts as 0.01, not 0
    removed_free = cu.adaptiveness(flat(CTA=0.0), code)
    assert removed_free["CTA"] == 0 and "TAA" not in removed_free


def test_the_editor_starts_empty_makes_edits_and_saves_a_table_of_yours(clean, monkeypatch, tmp_path):
    from PySide6.QtWidgets import QInputDialog
    import codon_dialog
    import tool_settings
    holder = codon_dialog.CodonDialog()
    dialog = holder.page
    assert dialog.current is None and "No codon usage table yet" in dialog.title.text()
    monkeypatch.setattr(QInputDialog, "getText", lambda *a, **k: ("Recoded host", True))
    dialog.new()
    assert dialog.current.name == "Recoded host" and dialog.boxes[0][2].isEnabled()
    dialog.boxes[0][2].setText("99.5")   # TTT
    dialog.boxes[1][2].setText("not a number")   # refused: the old number comes back
    assert dialog.boxes[1][2].text() in ("15.62", "15.63")
    dialog.use()
    holder.done(0)   # leaving saves the edit
    assert cu.all_tables()["Recoded host"].per_thousand["TTT"] == 99.5 and cu.active_name() == "Recoded host"
    # the Genetic code page of the Settings holds the choice, and OK stores it
    tool_settings.register_page(genetic_code.TOOL, codon_dialog.UsageBox)
    page = tool_settings.SettingsDialog(None, genetic_code.TOOL)
    box = page.custom[genetic_code.TOOL]
    assert [box.choice.itemText(i) for i in range(box.choice.count())] == ["Recoded host"] and box.choice.currentText() == "Recoded host"
    page.accept()
    assert cu.active_name() == "Recoded host"
    # export, then import under the same name
    holder = codon_dialog.CodonDialog()
    dialog = holder.page
    cu.write_csv(cu.all_tables()["Recoded host"], tmp_path / "out.csv")
    monkeypatch.setattr(codon_dialog.QFileDialog, "getOpenFileName", lambda *a, **k: (str(tmp_path / "out.csv"), ""))
    cu.delete_user("Recoded host")
    dialog.tables = cu.all_tables()
    dialog.import_csv()
    assert dialog.current.name == "Recoded host" and dialog.current.per_thousand["TTT"] == 99.5


def test_one_csv_carries_the_usage_table_and_the_genetic_code(clean, monkeypatch, tmp_path):
    """Export writes the frequencies and the code shown to one file; importing it back brings both, and a code already here that says
    the same thing is used instead of a second copy of it. A file with no letters (an older one) brings the table alone."""
    import codon_dialog
    holder = codon_dialog.CodonDialog()
    page = holder.page
    aas, starts = genetic_code.recoded(11, ["AGA", "AGG"])
    id = genetic_code.save_code("Export test code", aas, starts)
    path = tmp_path / "both.csv"
    try:
        cu.save_user(flat("Export test host", CTG=50.0, CTA=0.5))
        page.tables = cu.all_tables()
        page._fill_usage("Export test host")
        page._fill_codes(id)
        monkeypatch.setattr(codon_dialog.QFileDialog, "getSaveFileName", lambda *a, **k: (str(path), ""))
        page.export_csv()
        text = path.read_text(encoding="utf-8")
        assert "# code: Export test code" in text and "codon,aa,per_thousand,start" in text
        assert "AGA,X," in text and "ATG,M," in text   # the removed codon and the start codon travel too
        # read back: the frequencies and the letters both come through
        back = cu.read_csv(path)
        assert back.aas == aas and back.starts == starts and back.code_name == "Export test code"
        assert back.per_thousand["CTG"] == 50.0 and back.per_thousand["CTA"] == 0.5
        # importing it finds the code that already says the same thing instead of adding another
        cu.delete_user("Export test host")
        page.tables = cu.all_tables()
        monkeypatch.setattr(codon_dialog.QFileDialog, "getOpenFileName", lambda *a, **k: (str(path), ""))
        before = set(genetic_code.TABLES)
        page.import_csv()
        assert set(genetic_code.TABLES) == before and page.current_code.id == id
        assert page.current.per_thousand["CTG"] == 50.0
    finally:
        genetic_code.delete_code(id)
    # the same file once that code is gone: the code comes back as one of yours, named as the file says
    page.import_csv()
    new = page.current_code
    try:
        assert new.id >= genetic_code.FIRST_USER and new.name == "Export test code" and new.aas == aas
        assert "AGA" in genetic_code.removed_codons(new.id)
    finally:
        genetic_code.delete_code(new.id)
    # a file with no amino-acid column is still just a usage table
    plain = tmp_path / "plain.csv"
    plain.write_text("codon,per_thousand\n" + "\n".join(f"{c},{15.625}" for c in cu.CODONS), encoding="utf-8")
    assert cu.read_csv(plain).aas is None


def test_a_csv_with_a_bad_amino_acid_is_refused(tmp_path):
    rows = "\n".join(f"{c},{'B' if c == 'TTT' else 'F'},15.6" for c in cu.CODONS)
    path = tmp_path / "bad.csv"
    path.write_text("codon,aa,per_thousand\n" + rows, encoding="utf-8")
    with pytest.raises(ValueError, match="not an amino acid"):
        cu.read_csv(path)


def test_the_editor_says_so_when_a_table_cannot_be_saved(clean, monkeypatch):
    from PySide6.QtWidgets import QInputDialog, QMessageBox
    import codon_dialog
    import fileio
    shown = []
    monkeypatch.setattr(QMessageBox, "warning", staticmethod(lambda parent, title, text, *a, **k: shown.append(text)))
    monkeypatch.setattr(QInputDialog, "getText", lambda *a, **k: ("Locked one", True))
    holder = codon_dialog.CodonDialog()
    monkeypatch.setattr(fileio, "replace", lambda *a, **k: (_ for _ in ()).throw(PermissionError(13, "being used by another process")))
    holder.page.new()
    assert shown and "could not be saved" in shown[0] and "being used by another process" in shown[0]
    assert cu.all_tables() == {}


def test_plasmid_checkup_notes_a_poor_codon_adaptation_only(clean):
    import checkup
    cu.save_user(flat("Host", CTG=50.0, CTA=0.5))
    cu.set_active("Host")
    assert checkup.codon_notes("ATG" + "CTG" * 20 + "TAA") == []   # all best codons: nothing to say
    notes = checkup.codon_notes("ATG" + "CTA" * 20 + "TAA")
    assert len(notes) == 1 and "adaptation index 0.01 for Host" in notes[0] and "20 rare codons (CTA)" in notes[0]
    cu.set_active("Missing table")   # a table that is gone: another one is used, never an error
    assert isinstance(checkup.codon_notes("ATG" + "CTA" * 20 + "TAA"), list)


def test_recoded_genetic_codes_read_removed_codons_as_x_and_yours_are_kept_in_the_data_folder():
    import opendna_core
    assert genetic_code.removed_codons(101) == ["TCA", "TCG", "TAG"] and len(genetic_code.removed_codons(102)) == 7
    assert opendna_core.translate("ATGTCGTAGTCATAA", 101) == "MXXX*"   # TCG, TAG, TCA are not used; TAA still stops
    assert genetic_code.table_id() in genetic_code.TABLES
    aas, starts = genetic_code.recoded(11, ["AGA"])
    id = genetic_code.save_code("Test recoded host", aas, starts, note="made in a test")
    try:
        assert id >= genetic_code.FIRST_USER and f"{id} – Test recoded host" in genetic_code.LABELS
        assert (genetic_code.folder() / f"{id}.json").is_file()
        assert opendna_core.translate("AGAAGG", id) == "XR"   # AGA removed, AGG still arginine in table 11
        with pytest.raises(ValueError, match="already a genetic code"):
            genetic_code.save_code("test RECODED host", aas, starts)
        with pytest.raises(ValueError, match="built-in"):
            genetic_code.save_code("Edit a built-in", aas, starts, id=101)
        with pytest.raises(ValueError, match="64"):
            genetic_code.save_code("Short", "FFLL", "----")
    finally:
        genetic_code.delete_code(id)
    assert id not in genetic_code.TABLES and not (genetic_code.folder() / f"{id}.json").exists()


def test_the_quick_switch_lists_recoded_codes_after_a_divider_and_has_a_manage_entry():
    import code_switch
    switch = code_switch.CodeSwitch()
    texts = [a.text() for a in switch.menu.actions()]
    assert texts[-1] == "Manage codon tables…" and any(t.startswith("101 –") for t in texts)
    first_custom = next(i for i, a in enumerate(switch.menu.actions()) if a.text().startswith("101 –"))
    assert switch.menu.actions()[first_custom - 1].isSeparator()
    heard = []
    switch.manage.connect(lambda: heard.append(1))
    switch.menu.actions()[-1].trigger()
    assert heard == [1]


def test_the_chart_makes_a_recoded_code_of_yours(monkeypatch):
    from PySide6.QtWidgets import QInputDialog
    import codon_dialog
    holder = codon_dialog.CodonDialog()
    page = holder.page
    assert page.current_code.id == genetic_code.table_id()
    # The 64 codons are the chart, not a column: sixteen rows of four, each cell the codon with its letter, start box and frequency.
    assert (page.grid.rowCount(), page.grid.columnCount()) == (16, 4) and len(page.boxes) == 64
    assert codon_dialog.row_label(0) == "TT" and cu.CODONS[5] == "TCC" and page.grid.cellWidget(1, 1) is not None
    page._fill_codes(101)
    assert "built in" in page.title.text() and not page.boxes[0][0].isEnabled() and not page.save_button.isEnabled()
    assert page.note.text().startswith("61 codons used; not used: TCA, TCG, TAG") and "Fredens" in page.note.text()
    monkeypatch.setattr(QInputDialog, "getText", lambda *a, **k: ("Editor test code", True))
    page.duplicate()
    id = page.current_code.id
    try:
        assert id >= genetic_code.FIRST_USER and page.boxes[0][0].isEnabled() and page.name.text() == "Editor test code"
        page.boxes[cu.CODONS.index("TTT")][0].setCurrentText("X")   # TTT is not used either
        page.save()
        assert genetic_code.TABLES[id].codons["TTT"] == "X" and "TTT" in genetic_code.removed_codons(id)
        monkeypatch.setattr(codon_dialog.QMessageBox, "question", lambda *a, **k: codon_dialog.QMessageBox.Yes)
        page.delete()
        assert id not in genetic_code.TABLES
    finally:
        genetic_code.delete_code(id)


def test_checkup_flags_codons_a_recoded_code_does_not_use():
    import checkup
    import tool_settings
    saved = tool_settings.value(genetic_code.TOOL, "table")
    try:
        label = next(x for x in genetic_code.LABELS if x.startswith("101 "))
        tool_settings.set_value(genetic_code.TOOL, "table", label)
        problems = checkup.protein_problems("ATGTCGAAATCATAA")   # TCG and TCA are removed in Syn61
        assert "uses codons this genetic code does not use: TCA, TCG" in problems
        tool_settings.set_value(genetic_code.TOOL, "table", genetic_code.DEFAULT_LABEL)
        assert checkup.protein_problems("ATGTCGAAATCATAA") == []
    finally:
        tool_settings.set_value(genetic_code.TOOL, "table", saved)
