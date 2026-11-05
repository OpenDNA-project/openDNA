"""Cut sites: the hover text of restriction cuts and Cas9 / Cas12a guide sites, the cut marker (not the circle's origin mark),
the guide-site display and the way a saved guide site is stored in a GenBank file (a plain misc_feature other readers accept)."""
import os
import re
import sys
import tempfile
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ["OPENDNA_DATA_DIR"] = tempfile.mkdtemp(prefix="opendna_test_")
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "app"))

import pytest
from Bio import SeqIO
from PySide6.QtGui import QColor
from PySide6.QtWidgets import QApplication

from opendna_core import Document

import cas9
import cut_sites as cs
import enzymes
import hover
import theme
from assembly_panel import FragmentEditor
from sequence_view import SequenceView

app = QApplication.instance() or QApplication([])
pytestmark = pytest.mark.filterwarnings("ignore::Bio.BiopythonParserWarning")
ENZ = {e.name: e for e in enzymes.COMMON}
SP = cs.NUCLEASES["SpCas9 (NGG, 20 nt)"]
A12 = cs.NUCLEASES["Cas12a (TTTV, 23 nt)"]
GUIDE = "GACTCGTACGTAGCATCGAT"
FILL = "ATATTATAATTATATATTAA"


def plain(html):
    """The tooltip as the reader sees it."""
    return re.sub(r"<[^>]+>", "", html.replace("<br>", "\n"))


def tip(name, seq, circular=False, which=0):
    top, bottom = enzymes.find_cuts(seq, circular, [ENZ[name]])[name][which]
    return plain(cs.enzyme_tip(ENZ[name], seq, circular, top, bottom))


# ---- restriction enzyme tooltips -----------------------------------------------------------------------------------------


def test_enzyme_tooltips_by_kind_of_end():
    """The tooltip of a 5' overhang, a 3' overhang and a blunt enzyme, and of a Type IIS enzyme (how far from its site)."""
    # --- 5' overhang ---
    t = tip("EcoRI", "TTGAATTCAA")
    assert "EcoRI" in t and "G^AATTC" in t
    assert "Cuts after base 3 on the top strand and after base 7 on the bottom strand." in t
    assert "5′ overhang of 4 nt: AATT" in t
    assert "bases 3–8" in t and "Type IIS" not in t

    # --- 3' overhang ---
    t = tip("KpnI", "TTGGTACCAA")
    assert "GGTAC^C" in t
    assert "Cuts after base 7 on the top strand and after base 3 on the bottom strand." in t
    assert "3′ overhang of 4 nt: GTAC" in t

    # --- blunt ---
    t = tip("EcoRV", "TTGATATCAA")
    assert "GAT^ATC" in t and "Cuts after base 5 on the top strand and after base 5 on the bottom strand." in t
    assert "Blunt ends" in t and "overhang of" not in t

    # --- a Type IIS enzyme says how far from its site it cuts ---
    t = tip("BsaI", "TTGGTCTCAAGGAATTTT")   # GGTCTC at bases 3-8; cuts 1 nt after on the top strand, 5 after on the bottom
    assert "GGTCTC(1/5)" in t
    assert "Cuts after base 9 on the top strand and after base 13 on the bottom strand." in t
    assert "5′ overhang of 4 nt: AGGA" in t
    assert "Type IIS: it cuts outside its site, 1 nt after the site on the top strand and 5 nt after the site on the bottom strand." in t


def test_enzyme_tooltips_by_site_position():
    """The tooltip of a site on the reverse strand, a site across the origin of a circle, and every common enzyme."""
    # --- a site on the reverse strand ---
    t = tip("BsaI", "TTTTAAAAGGAAGAGACCTT")   # GAGACC = GGTCTC on the bottom strand, at bases 13-18; the cuts are before it
    assert "bases 13–18, on the bottom strand, where it reads GGTCTC; GAGACC on the top strand" in t
    assert "Type IIS" in t and "before the site" in t
    top, bottom = enzymes.find_cuts("TTTTAAAAGGAAGAGACCTT", False, [ENZ["BsaI"]])["BsaI"][0]
    assert (top, bottom) == (7, 11) and f"after base {top} on the top strand and after base {bottom} on the bottom" in t
    # a palindromic enzyme has one site, listed as the top strand
    assert "top strand)" in tip("EcoRI", "TTGAATTCAA")

    # --- a site across the origin of a circle ---
    seq = "ATTCAAAAAAGA"   # EcoRI GAATTC = GA at the end + ATTC at the start
    t = tip("EcoRI", seq, circular=True)
    assert "bases 11–4" in t
    assert "Cuts after base 11 on the top strand and after base 3 on the bottom strand." in t
    assert "5′ overhang of 4 nt: AATT" in t   # AA at the end, TT at the start

    # --- every common enzyme has a tooltip at each cut ---
    seq = "".join(e.site.replace("N", "A") for e in enzymes.COMMON) * 2
    for name, cuts in enzymes.find_cuts(seq, True, enzymes.COMMON).items():
        for top, bottom in cuts:
            t = plain(cs.enzyme_tip(ENZ[name], seq, True, top, bottom))
            assert name in t and "Cuts after base" in t and "Recognition site" in t


# ---- Cas9 / Cas12a tooltips and layout ---------------------------------------------------------------------------------


def guide_doc(strand, nuclease=SP, circular=False, lead=40, text=None):
    """A document with one guide site, as the Spacer Finder saves it."""
    site = GUIDE + "TGG"
    if text is not None:
        site = text
    if strand < 0:
        site = site[::-1].translate(str.maketrans("ACGT", "TGCA"))
    seq = FILL * 3 + site + FILL * 3
    lead = len(FILL) * 3
    n = len(seq)
    doc = Document.new_sequence("g", seq, circular)
    loc = f"{lead + 1}..{lead + len(site)}"
    doc.add_features([("misc_feature", f"complement({loc})" if strand < 0 else loc,
                       [("label", "gRNA GACTCGTA"), (cs.CAS9_KEY, cs.cas9_qualifier(nuclease))])])
    return doc, lead, n


def test_cas9_and_cas12a_tooltips():
    """The Cas9 tooltip and layout on both strands (protospacer, PAM and the blunt cut), and the Cas12a one with two cuts."""
    # --- Cas9 on both strands ---
    for strand in (1, -1):
        doc, lead, n = guide_doc(strand)
        (f,) = [f for f in doc.features() if f.kind == "misc_feature"]
        layout = cs.feature_layout(f, n)
        t = plain(cs.cas9_tip(layout, doc.sequence, n, "gRNA GACTCGTA"))
        assert f"5′-{GUIDE}-3′" in t and "PAM TGG (NGG)" in t and "SpCas9" in t, f"strand {strand}"
        assert "blunt cut" in t and "3 bp from the PAM" in t, f"strand {strand}"
        cut = lead + 20 - 3 if strand > 0 else lead + 3 + 3   # the gap, as the Spacer Finder computes it
        assert f"after base {cut} on both strands" in t, f"strand {strand}"
        assert cs.cut_places(SP, lead, strand) == (cut, cut), f"strand {strand}"
        assert layout.gaps == [cut], f"strand {strand}"
        # the protospacer and the PAM are where the strand says: PAM on the right for the top strand, on the left for the bottom
        if strand > 0:
            assert layout.proto == [(lead, lead + 20)] and layout.pam == [(lead + 20, lead + 23)]
        else:
            assert layout.pam == [(lead, lead + 3)] and layout.proto == [(lead + 3, lead + 23)]

    # --- Cas12a has two cuts ---
    site = "TTTA" + "GACTCGTACGTAGCATCGATAGC"   # PAM TTTV then 23 nt
    doc, lead, n = guide_doc(1, A12, text=site)
    (f,) = [f for f in doc.features() if f.kind == "misc_feature"]
    t = plain(cs.cas9_tip(cs.feature_layout(f, n), doc.sequence, n, "g"))
    assert "Cas12a" in t and "PAM TTTA (TTTV)" in t
    assert f"after base {lead + 4 + 18} on the top strand and after base {lead + 4 + 23} on the bottom strand" in t
    assert "5 nt 5′ overhang" in t


def test_guide_site_origin_qualifier_and_spacer_finder():
    """A guide site across the origin of a circle, the qualifier value read back (a damaged one is ignored), and the
    Spacer Finder and the view agreeing on the cut."""
    # --- a guide site across the origin of a circle ---
    seq = "ATATTATAATTATATATTAA" * 2
    site = GUIDE + "TGG"
    seq = site[10:] + seq + site[:10]   # the site starts 13 bases before the end and goes on at the start
    n = len(seq)
    doc = Document.new_sequence("c", seq, True)
    start = n - 10
    doc.add_features([("misc_feature", f"join({start + 1}..{n},1..13)", [("label", "x"), (cs.CAS9_KEY, cs.cas9_qualifier(SP))])])
    (f,) = [f for f in doc.features() if f.kind == "misc_feature"]
    layout = cs.feature_layout(f, n)
    assert layout.proto_cells[0] == start and layout.proto_cells[-1] == 9, "protospacer runs from the end across the origin"
    assert plain(cs.cas9_tip(layout, doc.sequence, n)).count(GUIDE) == 1

    # --- the qualifier value is read back and a damaged one is ignored ---
    for nuc in cs.NUCLEASES.values():
        assert cs.parse_cas9(cs.cas9_qualifier(nuc)) == nuc
    assert cs.parse_cas9("nonsense") is None and cs.parse_cas9("name=X;pam=NGG") is None and cs.parse_cas9("") is None

    # --- the Spacer Finder and the view agree on the cut ---
    seq = FILL * 3 + GUIDE + "TGG" + FILL
    (g,) = [g for g in cas9.candidates(seq, False, SP) if g.strand > 0]
    cells = [(g.start + i) % len(seq) for i in range(23)]
    assert cs.cas9_layout(SP, cells, 1, len(seq)).gaps == [g.cut]
    (r,) = [g for g in cas9.candidates(seq.translate(str.maketrans("ACGT", "TGCA"))[::-1], False, SP) if g.strand < 0]
    cells = [(r.start + i) for i in range(23)]
    assert cs.cas9_layout(SP, cells, -1, len(seq)).gaps == [r.cut]


# ---- the hover mechanism ------------------------------------------------------------------------------------------------


def shown(ed, size=(700, 300)):
    view = SequenceView(ed)
    view.resize(*size)
    view.show()
    app.processEvents()
    return view


def test_hover_reaches_features_cut_marks_and_unique_cutters():
    """hover.describe reaches features, cut marks and unique cutters; the view finds the cut and the enzyme name under
    the pointer."""
    # --- describe ---
    doc, lead, n = guide_doc(1)
    ed = FragmentEditor(doc, Path("t.gb"))
    (f,) = ed.features
    assert "SpCas9" in plain(hover.describe(ed, f.index))
    ed.set_cut_marks("digest", [(3, 7, "EcoRI", lambda: "tip text")])
    assert hover.describe(ed, ("cut", "digest", 0)) == "tip text"
    ed.set_cut_marks("x", [(3, 7, "EcoRI", "plain tip")])
    assert hover.describe(ed, ("cut", "x", 0)) == "plain tip"
    assert hover.describe(ed, ("cut", "x", 5)) == ""   # a mark that is gone is no error
    ed2 = FragmentEditor(Document.new_sequence("e", "TTGAATTCAA", False), Path("e.gb"))
    ed2.set_display(enzymes=True)
    (pos, name), = [s for s in ed2.shown_sites() if s[1] == "EcoRI"]
    assert "G^AATTC" in plain(hover.describe(ed2, ("site", name, pos)))

    # --- the view finds the cut and the enzyme name under the pointer ---
    ed = FragmentEditor(Document.new_sequence("e", "TTGAATTCAA" * 3 + "GGATCC" + "T" * 70, False), Path("e.gb"))
    ed.set_display(enzymes=True)
    ed.set_cut_marks("digest", [(3, 7, "EcoRI", "tip")])
    view = shown(ed)
    _, seq_top, _, _ = view._parts(0)
    top = view.offset(0)
    xt, xb = view.left + 3 * view.cw, view.left + 7 * view.cw
    assert view.cut_at(xt, top + seq_top + 2) == ("cut", "digest", 0)
    assert view.cut_at(xb, top + seq_top + 2 * view.lh - 2) == ("cut", "digest", 0)
    assert view.cut_at((xt + xb) / 2, top + seq_top + view.lh) == ("cut", "digest", 0)   # on the joining line
    assert view.cut_at(view.left + 15 * view.cw, top + seq_top + 2) is None
    assert view.hover_target(xt, top + seq_top + 2) == ("cut", "digest", 0)
    # the unique cutter's name on the line above the bases (BamHI is the only one that cuts once)
    (pos, name), = [s for s in ed.shown_sites() if s[1] == "BamHI"]
    x = view.left + pos * view.cw + 3
    assert view.site_label_at(x, top + 2) == ("site", "BamHI", pos)
    assert view.site_label_at(x, top + seq_top + 2) is None
    view.close()


# ---- painting -------------------------------------------------------------------------------------------------------------


def count(view, colour, rect=None):
    img = view.viewport().grab().toImage()
    rect = rect or (0, 0, img.width(), img.height())
    x0, y0, x1, y1 = (int(v) for v in rect)
    return sum(1 for x in range(max(x0, 0), min(x1, img.width())) for y in range(max(y0, 0), min(y1, img.height()))
               if QColor(img.pixel(x, y)).rgb() == colour.rgb())


def test_the_cut_marker_is_not_the_origin_mark(monkeypatch):
    """In the light and the dark theme the cut ends in a cross, while the origin mark is a bare bar of another colour."""
    for dark in (False, True):
        monkeypatch.setattr(theme, "is_dark", lambda: dark)
        ed = FragmentEditor(Document.new_sequence("c", "ACGT" * 200, True), Path("c.gb"))
        view = shown(ed)
        bpl = view.bases_per_line()
        assert view.rpc >= 3
        ed.set_cut_marks("digest", [(bpl + 10, bpl + 10, "EcoRV")])
        app.processEvents()
        scroll = view.verticalScrollBar().value()
        colour = theme.color("cut")
        cx = view.left + 10 * view.cw
        top = view.offset(1) + view._parts(1)[1] - scroll   # y of the bases of the row with the cut
        # Around the cut, the crosses reach either side of the stroke (above the top strand and below the bottom one) ...
        for y0, y1 in ((top - 8, top - 1), (top + 2 * view.lh + 1, top + 2 * view.lh + 8)):
            assert count(view, colour, (cx - 4.5, y0, cx - 1.5, y1)) > 0 and count(view, colour, (cx + 1.5, y0, cx + 4.5, y1)) > 0, \
                f"the cut ends in a cross (dark={dark})"
        # ... while the origin mark (the last row, after its last base) is a bare vertical bar of a different colour
        view.verticalScrollBar().setValue(view.verticalScrollBar().maximum())
        app.processEvents()
        scroll = view.verticalScrollBar().value()
        from sequence_view import ORIGIN_COLOR
        origin = QColor(ORIGIN_COLOR)
        last = view.rpc - 1
        ox = view.left + (ed.length - last * bpl) * view.cw + 1
        otop = view.offset(last) + view._parts(last)[1] - scroll
        assert count(view, origin, (ox - 4.5, otop - 8, ox - 1.5, otop + 2 * view.lh + 8)) == 0, f"the origin mark has no cross (dark={dark})"
        assert count(view, origin, (ox - 1, otop - 1, ox + 2, otop + 2 * view.lh)) > 0, f"(the origin bar is there) (dark={dark})"
        assert theme.color("cut").rgb() != origin.rgb()
        view.close()


def test_a_saved_guide_site_is_drawn(monkeypatch):
    """A saved guide site is drawn with protospacer, PAM, name and cut arrow on both strands in both themes (PAM on the
    right of the protospacer for the top strand); edited to another size it is drawn as a plain feature."""
    # --- protospacer, PAM and name, both strands, both themes ---
    for dark in (False, True):
        for strand in (1, -1):
            monkeypatch.setattr(theme, "is_dark", lambda: dark)
            doc, lead, n = guide_doc(strand)
            ed = FragmentEditor(doc, Path("g.gb"))
            view = shown(ed, (1800, 300))
            for role in ("cas9_proto", "cas9_pam", "cas9_label"):
                assert count(view, theme.color(role)) > 50, f"{role} is not painted (dark={dark}, strand={strand})"
            # the boxes: the protospacer is 20 bases wide, the PAM 3
            bpl = view.bases_per_line()
            assert lead + 23 < bpl, "the whole site is in the first row for this test"
            lane_y = int(view._parts(0)[3]) + 2
            img = view.viewport().grab().toImage()
            proto, pam = theme.color("cas9_proto").rgb(), theme.color("cas9_pam").rgb()
            # Box widths are measured column by column over the lane's height, not along one pixel row: the "PAM" name is
            # written inside its box, and on one row through the letters the glyphs of a wide font hide most of the box.
            band = range(max(lane_y - 2, 0), min(lane_y + 8, img.height()))
            colours = [{QColor(img.pixel(x, y)).rgb() for y in band} for x in range(img.width())]
            row = [proto if proto in c else pam if pam in c else None for c in colours]
            assert row.count(proto) > 15 * view.cw and row.count(pam) > 2 * view.cw, f"dark={dark}, strand={strand}"
            assert (row.index(pam) > row.index(proto)) == (strand > 0), \
                f"the PAM is on the right of the protospacer on the top strand (dark={dark}, strand={strand})"
            view.close()

    # --- a guide site edited to another size is drawn as a plain feature ---
    monkeypatch.undo()
    doc, lead, n = guide_doc(1)
    ed = FragmentEditor(doc, Path("g.gb"))
    (f,) = ed.features
    ed.edit_features([(f.index, ("misc_feature", f"{lead + 1}..{lead + 10}", list(f.qualifiers)))])
    assert cs.feature_layout(ed.features[0], ed.length) is None
    view = shown(ed)
    assert count(view, theme.color("cas9_pam")) == 0
    view.close()


# ---- storage ----------------------------------------------------------------------------------------------------------------


def saved_guide(tmp_path, circular=False):
    """(editor, guide) after the Spacer Finder's own "Add selected as feature" saved a guide."""
    import main as M
    win = M.MainWindow()
    win.show()
    win.show_document(Document.new_sequence("g", FILL * 3 + GUIDE + "TGG" + FILL, circular), Path(tmp_path) / "g.gb", untitled=True)
    tool = cas9.Cas9Dialog(win, win.editor)
    (row,) = [i for i, g in enumerate(tool.guides) if g.strand > 0]
    tool.table.setCurrentCell(row, 0)
    tool.add()
    ed, g = win.editor, tool.guides[row]
    assert ed.cut_marks["cas9"][0].tip() and "SpCas9" in plain(ed.cut_marks["cas9"][0].tip()), "the tool's cut hovers too"
    tool.close()
    win.editor.doc.modified = False
    win.reset_layout()
    return ed, g, win


def test_a_saved_guide_in_the_file(tmp_path):
    """A saved guide round-trips byte for byte and Biopython reads it as a plain misc_feature; a reverse-strand guide
    keeps its strand in the file; exports leave the program's qualifier out."""
    # --- round trip, byte for byte, and a plain reader ---
    ed, g, win = saved_guide(tmp_path)
    p1, p2 = tmp_path / "a.gb", tmp_path / "b.gb"
    ed.doc.save_genbank(str(p1))
    back = Document(str(p1))
    back.save_genbank(str(p2))
    assert p1.read_bytes() == p2.read_bytes(), "load and save again must not change the file"
    (f,) = [f for f in back.features() if f.kind == "misc_feature"]
    layout = cs.feature_layout(f, back.length)
    assert layout is not None and layout.gaps == [g.cut]
    assert "".join(back.sequence(a, b) for a, b in layout.proto) == g.guide
    assert "".join(back.sequence(a, b) for a, b in layout.pam) == g.pam
    # a plain reader: Biopython parses it without errors, as a misc_feature with a label, a readable note and the program's qualifier
    with open(p1, encoding="utf-8") as fh:
        rec = SeqIO.read(fh, "genbank")
    (feat,) = [x for x in rec.features if x.type == "misc_feature"]
    assert feat.qualifiers["label"] == [f"gRNA {g.guide[:8]}"]
    assert g.guide in feat.qualifiers["note"][0] and f"PAM {g.pam}" in feat.qualifiers["note"][0]
    assert feat.qualifiers[cs.CAS9_KEY] == [cs.cas9_qualifier(SP)]
    assert int(feat.location.start) == g.start and int(feat.location.end) == g.end and feat.location.strand == 1
    text = p1.read_text(encoding="utf-8")
    assert f'/{cs.CAS9_KEY}="{cs.cas9_qualifier(SP)}"' in text
    assert all(len(line) <= 80 for line in text.splitlines())

    # --- a reverse-strand guide keeps its strand in the file ---
    doc, lead, n = guide_doc(-1)
    path = tmp_path / "r.gb"
    doc.save_genbank(str(path))
    with open(path, encoding="utf-8") as fh:
        (feat,) = [x for x in SeqIO.read(fh, "genbank").features if x.type == "misc_feature"]
    assert feat.location.strand == -1 and int(feat.location.start) == lead
    back = Document(str(path))
    layout = cs.feature_layout([f for f in back.features() if f.kind == "misc_feature"][0], back.length)
    assert layout.pam == [(lead, lead + 3)]

    # --- exports leave the program qualifier out ---
    import exports
    assert cs.CAS9_KEY not in exports.tbl_text(ed) and cs.CAS9_KEY not in exports.gff3_text(ed)
    assert g.guide in exports.tbl_text(ed), "the plain note stays"
