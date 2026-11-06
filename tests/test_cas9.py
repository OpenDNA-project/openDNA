"""Cas9 guide finder: the engine (candidates, score, other sites here and in a genome) and the tool pane."""
import os
import random
import sys
import tempfile
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ["OPENDNA_DATA_DIR"] = tempfile.mkdtemp(prefix="opendna_test_")
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "app"))

from PySide6.QtTest import QTest
from PySide6.QtWidgets import QApplication

import opendna_core as core
from opendna_core import Document

import cas9
from editor import reverse_complement as rc

app = QApplication.instance() or QApplication([])
rng = random.Random(7)
SP = cas9.NUCLEASES["SpCas9 (NGG, 20 nt)"]
SA = cas9.NUCLEASES["SaCas9 (NNGRRT, 21 nt)"]
GUIDE = "GACTCGTACGTAGCATCGAT"   # 20 nt, 50 % GC, no GG or CC in it: no PAM of its own on either strand


def dna(n):
    """Filler of A and T only: no G or C, so no PAM on either strand (the sites are the ones a test puts in)."""
    return "".join(rng.choice("AT") for _ in range(n))


def bend(guide, positions):
    """The guide with the bases at `positions` changed to another base."""
    out = list(guide)
    for i in positions:
        out[i] = "ACGT"["ACGT".index(out[i]) - 1]
    return "".join(out)


def test_candidates_circles_other_nucleases_score_and_seed_mismatches():
    """Guides on both strands with their cut; an offset region, circles and SaCas9; the score is a plain checklist (best guides first); seed mismatches."""
    # --- guides on both strands with their cut ---
    top = dna(40) + GUIDE + "TGG" + dna(40)
    found = cas9.candidates(top, False, SP)
    assert len(found) == 1
    g = found[0]
    assert (g.start, g.end, g.strand, g.guide, g.pam) == (40, 63, 1, GUIDE, "TGG")
    assert g.cut == 40 + 20 - 3, "Cas9 cuts 3 bp from the PAM, inside the protospacer"
    # The same site on the bottom strand: the top strand holds its reverse complement, PAM on the left.
    bottom = dna(40) + rc(GUIDE + "AGG") + dna(40)
    (r,) = cas9.candidates(bottom, False, SP)
    assert (r.start, r.end, r.strand, r.guide, r.pam) == (40, 63, -1, GUIDE, "AGG")
    assert r.cut == 40 + 3 + 3, "3 bp inside the protospacer, counted from the PAM on the left"
    assert bottom[r.start:r.end] == rc(r.guide + r.pam)

    # --- offset region, circles and other nucleases ---
    seq = dna(60) + GUIDE + "CGG" + dna(30)
    (g,) = cas9.candidates(seq[50:], False, SP, offset=50)
    assert g.start == 60 and g.guide == GUIDE, "a stretch is read at its place on the whole sequence"
    # A site across the origin of a circle.
    site = GUIDE + "TGG"
    rotated = site[12:] + dna(80) + site[:12]
    (c,) = cas9.candidates(rotated, True, SP, wrap=True)
    assert c.guide == GUIDE and c.start == len(rotated) - 12 and c.end > len(rotated)
    assert not cas9.candidates(rotated, False, SP), "a line has no site across its ends"
    # SaCas9: 21 nt and an NNGRRT PAM.
    sa = GUIDE + "A" + "ATGAAT"
    (s,) = [c for c in cas9.candidates(dna(30) + sa + dna(30), False, SA) if c.strand > 0]
    assert len(s.guide) == 21 and s.pam == "ATGAAT"
    assert not [c for c in cas9.candidates(dna(30) + GUIDE + "A" + "ATGCCT" + dna(30), False, SA) if c.strand > 0], "ATGCCT is no NNGRRT"

    # --- the score is a plain checklist ---
    assert cas9.score_guide(GUIDE) == (100, [])
    s, notes = cas9.score_guide("ACGTACGTTTTTACGTACGA")
    assert s < 70 and any("TTTT" in n for n in notes)
    s, notes = cas9.score_guide("GGGGGGGGGGGGGGGGGGGG")
    assert s <= 45 and any("GC" in n for n in notes) and any("run" in n for n in notes)
    s, notes = cas9.score_guide("ATATATATATTATAATATAT")
    assert s <= 60 and any("GC" in n for n in notes)
    assert cas9.score_guide("ACGATCGTAGCTAGCTAGCA")[0] == 100
    # The best guides come first.
    seq = dna(30) + GUIDE + "TGG" + dna(20) + "ACGTACGTTTTTACGTACGA" + "TGG" + dna(30)
    found = cas9.candidates(seq, False, SP)
    assert [g.score for g in found] == sorted((g.score for g in found), reverse=True) and found[0].guide == GUIDE

    # --- seed mismatches ---
    assert cas9.seed_mismatches(1, 20) == 0, "base 0 is at the far end"
    assert cas9.seed_mismatches(1 << 19, 20) == 1, "base 19 touches the PAM"
    assert cas9.seed_mismatches((1 << 19) | (1 << 7) | 1, 20) == 1 and cas9.seed_mismatches((1 << 8) | (1 << 12), 20) == 2


def test_other_sites_here_and_in_a_genome_and_the_rust_search():
    """Counts the other sites here and in a genome (with the verdicts); the Rust search matches a plain Python one."""
    # --- the other sites here and in a genome ---
    near = bend(GUIDE, [2])               # one mismatch, far from the PAM
    far = bend(GUIDE, [0, 1, 5, 19])      # four mismatches, one next to the PAM
    seq = dna(30) + GUIDE + "TGG" + dna(30) + near + "AGG" + dna(30) + far + "CGG" + dna(30)
    (g,) = [x for x in cas9.candidates(seq, False, SP) if x.guide == GUIDE]
    (spec,) = cas9.scan(seq, False, "plasmid", [g], SP, 3, None)
    assert spec.here == [0, 1, 0, 0] and spec.genome is None, "one near site; the four-mismatch one is beyond 3"
    assert cas9.verdict(spec, g)[0] == "warn"
    # A genome holding the guide twice (a perfect copy and one with two mismatches).
    two = bend(GUIDE, [3, 9])
    ref = Document.new_sequence("ref", dna(50) + GUIDE + "AGG" + dna(40) + two + "TGG" + dna(50), True)
    (spec,) = cas9.scan(seq, False, "plasmid", [g], SP, 3, ref)
    assert spec.genome == [1, 0, 1, 0] and spec.genome_name == "ref"
    assert cas9.verdict(spec, g)[0] == "warn"
    # Two perfect copies in the genome, or another perfect site here, is bad.
    ref2 = Document.new_sequence("ref2", GUIDE + "AGG" + dna(30) + GUIDE + "CGG" + dna(30), False)
    (spec,) = cas9.scan(seq, False, "plasmid", [g], SP, 3, ref2)
    assert spec.genome[0] == 2 and cas9.verdict(spec, g)[0] == "bad"
    twice = dna(20) + GUIDE + "TGG" + dna(20) + GUIDE + "AGG" + dna(20)
    (g1, g2) = sorted(cas9.candidates(twice, False, SP), key=lambda x: x.start)
    (s1, s2) = cas9.scan(twice, False, "p", [g1, g2], SP, 3, None)
    assert s1.here[0] == 1 and s2.here[0] == 1 and cas9.verdict(s1, g1)[0] == "bad"
    # The genome is the sequence itself: it is not scanned twice.
    same = Document.new_sequence("plasmid", seq, False)
    (spec,) = cas9.scan(seq, False, "plasmid", [g], SP, 3, same)
    assert spec.genome is None
    assert cas9.verdict(None, g)[0] == "info"

    # --- the Rust search matches a plain Python one ---
    seq = "".join(rng.choice("ACGT") for _ in range(4000))
    guides = [seq[i:i + 20] for i in (100, 1000, 2500)]
    hits = core.find_guide_sites(seq, False, guides, "NGG", 2)
    for guide, found in zip(guides, hits):
        expect = set()
        for strand, s in ((1, seq), (-1, rc(seq))):
            for i in range(len(s) - 22):
                if s[i + 21:i + 23] == "GG" and sum(a != b for a, b in zip(s[i:i + 20], guide)) <= 2:
                    expect.add((i if strand > 0 else len(seq) - i - 23, strand))
        assert {(h[0], h[2]) for h in found} == expect


# ---- the tool -------------------------------------------------------------------------------------------------------


def open_doc(doc):
    import main as M
    win = M.MainWindow()
    win.show()
    win.show_document(doc, Path(tempfile.mkdtemp()) / "t.gb", untitled=True)
    return win


def close(win):
    """Close a test window without the "save changes?" question a new document would ask."""
    win.editor.doc.modified = False
    win.reset_layout()   # closing saves the tabs: leave the default ones for the next window (tests share a data folder)
    win.close()


def wait_for(cond, ms=5000):
    QTest.qWait(50)
    waited = 50
    while not cond() and waited < ms:
        QTest.qWait(50)
        waited += 50
    assert cond()


def test_the_tool_in_a_pane_with_cas9_cas12a_and_the_reference_genome(tmp_path, monkeypatch):
    """The tool lists, selects, checks and adds a guide; it is a pane view in the Tools group; it draws Cas12a's sticky cut; its reference genome is the one of the primers tool."""
    monkeypatch.setattr(cas9, "genome", lambda: None)   # whatever genome this machine has set
    # --- the tool lists, selects, checks and adds a guide ---
    seq = dna(80) + GUIDE + "TGG" + dna(60) + GUIDE + "AGG" + dna(60)
    win = open_doc(Document.new_sequence("t", seq, False))
    ed = win.editor
    tool = cas9.Cas9Dialog(win, ed)
    assert tool.table.rowCount() == 2 and "2 guide" in tool.status.text()
    tool.table.selectRow(0)
    tool.table.setCurrentCell(0, 0)
    assert (ed.start, ed.end) == (tool.guides[0].start, tool.guides[0].end), "selecting a row selects the site in the sequence"
    (mark,) = ed.cut_marks["cas9"]
    assert mark.top == mark.bottom == tool.guides[0].cut == 80 + 20 - 3, "Cas9 cuts both strands 3 bp from the PAM: drawn there"
    # No genome set: only this sequence is checked, and the two copies find each other.
    tool.start_check()
    wait_for(lambda: tool.check.isEnabled() and "–" not in tool.table.item(0, 6).text())
    assert tool.table.item(0, 6).text() == "1/0/0/0" and tool.table.item(0, 7).text() == "–"
    assert "✗" in tool.table.item(0, 8).text() and "other perfect site" in tool.table.item(0, 8).text()
    # Add it as a feature: an arrow on the sequence, with the guide in its note.
    before = len(ed.doc.features())
    tool.table.setCurrentCell(0, 0)
    tool.add()
    feats = ed.doc.features()
    assert len(feats) == before + 1
    f = feats[-1]
    assert f.kind == "misc_feature" and f.label.startswith("gRNA ") and "SpCas9 guide " + GUIDE in dict(f.qualifiers)["note"]
    assert dict(f.qualifiers)["note"].count("cuts 3 bp from the PAM")
    # A selection narrows the search.
    ed.select(0, 60)
    tool.reload()
    assert tool.table.rowCount() == 0
    assert not ed.cut_marks, "no guide, no cut drawn"
    close(win)

    # --- the tool is a pane view in the tools group ---
    import workspace
    assert "Spacer Finder" in dict(workspace.grouped_views())["Tools"]
    win = open_doc(Document.new_sequence("t", dna(100) + GUIDE + "TGG" + dna(50), False))
    pane = win.workspace.ordered_panes(win.workspace.widget(0))[0]
    pane.set_view("Spacer Finder")
    assert isinstance(pane.view or pane.build(), cas9.Cas9Dialog) and "Spacer Finder" in win.workspace.tabText(0)
    assert win.workspace.views(cas9.Cas9Dialog)
    close(win)

    # --- Cas12a in the tool ---
    seq = FILL + "TTTA" + GUIDE23 + FILL
    win = open_doc(Document.new_sequence("t12", seq, False))
    ed = win.editor
    tool = cas9.Cas9Dialog(win, ed)
    tool.nuclease.setCurrentText("Cas12a (TTTV, 23 nt)")
    assert tool.table.rowCount() == 1 and tool.table.item(0, 2).text() == "TTTA " + GUIDE23
    tool.table.setCurrentCell(0, 0)
    (mark,) = ed.cut_marks["cas9"]
    assert (mark.top, mark.bottom, mark.label) == (50 + 4 + 18, 50 + 4 + 23, "Cas12a"), "drawn as a sticky cut"
    tool.add()
    note = dict(ed.doc.features()[-1].qualifiers)["note"]
    assert note.startswith("Cas12a guide " + GUIDE23) and "5′ overhang" in note
    close(win)

    monkeypatch.undo()   # the genome section below needs the real genome lookup, not the "no genome" stub
    # --- the reference genome is the one of the primers tool ---
    import primer_panel
    from PySide6.QtCore import QSettings
    # The primer options live in the user's own settings: use a throw-away file instead of touching them.
    s = QSettings(str(tmp_path / "options.ini"), QSettings.IniFormat)
    monkeypatch.setattr(primer_panel, "settings", lambda: s)
    ref_path = tmp_path / "ref.gb"
    seq = dna(80) + GUIDE + "TGG" + dna(80)
    Document.new_sequence("myref", dna(30) + GUIDE + "CGG" + dna(30), False).save_genbank(ref_path)
    win = open_doc(Document.new_sequence("plasmid", seq, False))
    s.setValue("primers/genome_check", "true")
    s.setValue("primers/genome", str(ref_path))
    pane = win.workspace.ordered_panes(win.workspace.widget(0))[0]
    pane.set_view("Spacer Finder")   # the tool in a pane of the window, where the Primers options reach it
    tool = pane.view or pane.build()
    win.primer_options_changed()
    assert "myref" in tool.genome_label.text() and cas9.genome() is primer_panel.genome(), "one genome for both tools"
    tool.start_check()
    wait_for(lambda: tool.check.isEnabled() and tool.table.item(0, 7).text() != "–")
    assert tool.table.item(0, 7).text() == "1/0/0/0", "the genome holds the guide once"
    assert tool.table.item(0, 8).text().startswith("✓")
    # Turning the genome check off in the Primers options reaches the tool.
    s.setValue("primers/genome_check", "false")
    win.primer_options_changed()
    assert "none set" in tool.genome_label.text()
    close(win)


# ---- Cas12a: the PAM (TTTV) is on the 5' side, and the two strands are cut in different places ---------------------------------------


A12 = cas9.NUCLEASES["Cas12a (TTTV, 23 nt)"]
GUIDE23 = GUIDE + "GCA"      # 23 nt; no TTT in it or in its reverse complement's PAM sense
FILL = "CA" * 25             # no T run, no AAA: no PAM of its own on either strand


def test_cas12a_sites_staggered_cuts_and_other_sites_counted_with_the_pam_first():
    """Cas12a: the PAM (TTTV) is on the 5' side, the two strands are cut in different places, other sites are counted with the PAM first."""
    # --- sites and staggered cuts ---
    top = FILL + "TTTA" + GUIDE23 + FILL
    (g,) = cas9.candidates(top, False, A12)
    assert (g.start, g.end, g.strand, g.guide, g.pam) == (50, 77, 1, GUIDE23, "TTTA")
    assert (g.cut, g.bottom) == (50 + 4 + 18, 50 + 4 + 23), "the PAM's strand is cut after 18 bases, the other after 23: a 5 nt 5' overhang"
    # On the bottom strand the PAM is on the right of the site on the top strand.
    bottom = FILL + rc("TTTC" + GUIDE23) + FILL
    (r,) = cas9.candidates(bottom, False, A12)
    assert (r.start, r.end, r.strand, r.guide, r.pam) == (50, 77, -1, GUIDE23, "TTTC")
    assert (r.cut, r.bottom) == (50, 55) and bottom[r.start:r.end] == rc(r.pam + r.guide)
    assert [x.start for x in cas9.candidates(FILL + "TTTT" + GUIDE23 + FILL, False, A12)] == [51], "TTTT is no PAM; the TTTG after its first T is"
    assert "5′ overhang" in cas9.cut_text(A12, g, len(top)) and "3 bp" in cas9.cut_text(SP, cas9.candidates(GUIDE + "TGG" + FILL, False, SP)[0], 99)

    # --- other sites are counted with the PAM first ---
    text = FILL + "TTTA" + GUIDE23 + FILL + "TTTC" + GUIDE23 + FILL
    guides = cas9.candidates(text, False, A12)
    assert len(guides) == 2
    spec = cas9.scan(text, False, "t", guides[:1], A12, 1)[0]
    assert spec.here == [1, 0], "the other copy, found with the PAM before the protospacer; the guide itself is left out"
