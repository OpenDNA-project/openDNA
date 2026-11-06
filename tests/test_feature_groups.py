"""The Features view's gene groups (a gene and its CDS sharing a /locus_tag): opening or closing one changes only that gene's row
and the rows of its parts under it, the same rows a full rebuild gives, whatever the sort and the filter, and without measuring
the columns again (it used to rebuild every row and reset the table: most of a second on a genome)."""
import os
import random
import sys
import tempfile
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ.setdefault("OPENDNA_DATA_DIR", tempfile.mkdtemp(prefix="opendna_test_"))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "app"))

from PySide6.QtCore import Qt
from PySide6.QtWidgets import QApplication

from opendna_core import Document

from editor import Editor
from feature_table import FeatureModel, FeatureTable

app = QApplication.instance() or QApplication([])


def genbank(path):
    """1,200 bp with six genes (each a gene + CDS sharing a /locus_tag; gene F has two parts) and a feature in no gene."""
    feats = []
    for i, name in enumerate("ABCDEF"):
        s = 1 + i * 180
        loc = f"{s}..{s + 149}" if i % 2 == 0 else f"complement({s}..{s + 149})"
        feats.append(f"     gene            {loc}\n                     /locus_tag=\"T_{name}\"\n                     /gene=\"gen{name}\"")
        feats.append(f"     CDS             {loc}\n                     /locus_tag=\"T_{name}\"\n                     /label=cds{name}")
    feats.append("     misc_feature    901..920\n                     /locus_tag=\"T_F\"\n                     /label=tagF")
    feats.append("     misc_feature    1101..1150\n                     /label=lonely")
    seq = ("ATGAAACGCATTAGCACCACC" * 60)[:1200]
    origin = "\n".join(f"{i + 1:>9} " + " ".join(seq[j:j + 10] for j in range(i, min(i + 60, 1200), 10))
                       for i in range(0, 1200, 60))
    path.write_text("LOCUS       groups      1200 bp    DNA     linear   UNK 01-JAN-1980\nFEATURES             Location/Qualifiers\n"
                    + "\n".join(feats) + "\nORIGIN\n" + origin + "\n//\n")
    return path


def shown(model):
    return [(r[0], r[3]) for r in model.shown]


def rebuilt(model):
    """What a full rebuild shows with the model's expanded genes, filter and sort."""
    fresh = FeatureModel(model.editor)
    fresh.expanded = set(model.expanded)
    fresh.filter_text, fresh.sort_column, fresh.sort_order = model.filter_text, model.sort_column, model.sort_order
    fresh.reload()
    return shown(fresh)


def test_opening_a_gene_changes_only_its_rows_and_matches_a_rebuild(tmp_path):
    ed = Editor(Document(str(genbank(tmp_path / "g.gb"))), tmp_path / "g.gb")
    assert len(ed.groups) == 6
    table = FeatureTable(ed)
    table.resize(800, 500)
    table.show()
    app.processEvents()
    model = table.model
    events = []
    for signal in ("modelReset", "rowsInserted", "rowsRemoved"):
        getattr(model, signal).connect(lambda *a, s=signal: events.append(s))
    measured = []
    fitter = table.view._fitter
    fitter.timer.timeout.connect(lambda: measured.append(1))
    # one gene opened, then closed
    row = next(i for i, r in enumerate(model.shown) if r[4] == "T_F")
    before = len(model.shown)
    model.toggle(row)
    app.processEvents()
    assert len(model.shown) == before + 2 and model.shown[row][0][0].startswith("▾")
    assert [model.shown[row + k][3] for k in (1, 2)] == [p.index for p in ed.groups["T_F"].parts]
    assert shown(model) == rebuilt(model)
    model.toggle(row)
    app.processEvents()
    assert len(model.shown) == before and model.shown[row][0][0].startswith("▸") and shown(model) == rebuilt(model)
    assert "modelReset" not in events and events.count("rowsInserted") == 1 and events.count("rowsRemoved") == 1
    assert not measured and not fitter.timer.isActive(), "the columns are not measured again for a gene opened or closed"
    # random opening and closing, under every sort and with a filter that matches only some parts
    rng = random.Random(4)
    for column in range(len(model.columns)):
        for order in (Qt.AscendingOrder, Qt.DescendingOrder):
            model.sort(column, order)
            for text in ("", "cds", "tag"):
                model.set_filter(text)
                for _ in range(8):
                    genes_on_view = [i for i, r in enumerate(model.shown) if r[4]]
                    if genes_on_view:
                        model.toggle(rng.choice(genes_on_view))
                    assert shown(model) == rebuilt(model), (column, order, text)
    # a part selected elsewhere (the map) opens its gene
    model.set_filter("")
    model.expanded.clear()
    model.reload()
    cds = ed.groups["T_C"].parts[0].index
    at = model.row_of(cds)
    assert model.shown[at][3] == cds and "T_C" in model.expanded and shown(model) == rebuilt(model)
    # an edit still rebuilds and refits
    measured.clear()
    model.reload()
    app.processEvents()
    assert measured, "a real change of the rows is measured again"
