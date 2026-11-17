"""Makes the README pictures (docs/images/*.png) from the real program and the real sequences in demo/, so they can be made
again after any change to the interface:

    python scripts/make_screenshots.py            every picture
    python scripts/make_screenshots.py map cas9   only these

Run it with real fonts (QT_QPA_PLATFORM unset: windows open briefly on screen). The settings go to a temporary data folder,
never yours. The "tracks" picture maps simulated reads with bwa and samtools (in WSL on Windows, as Tools > Map Reads does).
"""
import os
import sys
import tempfile
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DEMO, OUT = ROOT / "demo", ROOT / "docs" / "images"
sys.path.insert(0, str(ROOT / "app"))
os.environ.setdefault("OPENDNA_DATA_DIR", tempfile.mkdtemp(prefix="opendna_screens_"))   # before the app reads any setting

from PySide6.QtWidgets import QApplication  # noqa: E402

app = QApplication.instance() or QApplication(sys.argv[:1])
app.setApplicationName("openDNA")

import main as M  # noqa: E402
import theme  # noqa: E402
from workspace import pane, split  # noqa: E402

SIZE = (1600, 940)


def settle(seconds=0.6):
    """Let the window lay itself out and paint (and background work finish)."""
    end = time.time() + seconds
    while time.time() < end:
        app.processEvents()
        time.sleep(0.02)


def window(path, layout=None, tool=None):
    """`tool`: a view that gets a view tab of its own, shown full size (the first tab keeps `layout`, or the defaults)."""
    win = M.new_window()
    win.window().resize(*SIZE)
    assert win.load(str(path)), path
    if tool is not None:
        tabs = [{"layout": layout if layout is not None else M.DEFAULT_LAYOUT[0]["layout"]}, {"layout": pane(tool)}]
        win._build_workspace(tabs, current=1)
    elif layout is not None:
        win._build_workspace([{"layout": layout}])
    settle()
    return win


def save(win, name):
    settle()
    path = OUT / f"{name}.png"
    win.window().grab().save(str(path))
    win.window().close()
    print(f"{path.relative_to(ROOT)}")


# ------------------------------------------------------------------------------------------------------ the pictures

def editor():
    """The default layout on a real plasmid: sequence, map and feature table."""
    save(window(DEMO / "pBAD30_eilR.gb"), "editor")


def multi_genome():
    """Vibrio cholerae: two chromosomes in one GenBank file; the Map lists them, a click switches."""
    win = window(DEMO / "vibrio_cholerae_N16961.gb", tool="Map")
    save(win, "multi_genome")


def gibson():
    """The Assembly tab, Gibson: pBAD30 rebuilt from its own backbone and its eilR insert, a card per junction."""
    import assembly as A
    win = window(DEMO / "pBAD30_eilR.gb")
    win.show_assembly()
    panel = win.assembly
    panel.method.setCurrentIndex(0)          # Gibson / HiFi
    panel.count.setValue(2)
    doc = win.editor.doc
    for i, (start, end) in enumerate([(0, 5118), (5118, doc.length)]):   # the pBAD30 backbone, then the eilR insert
        panel._set_source(i, A.snapshot(doc), "pBAD30_eilR.gb", None)
        panel._use(i, start, end)
    panel.circular.setChecked(True)
    panel.add_homology()                     # design the overlaps and their primers
    settle(1.2)
    save(win, "assembly_gibson")


def primers():
    """The PCR tab on pBAD30, full size, with the plasmid's own Verif pBAD30_F / _R primers."""
    win = window(DEMO / "pBAD30_eilR.gb", tool="PCR")
    pcr = [p.view for p in win.workspace.panes() if p.view_name == "PCR"][0]
    fwd = next(i for i, s in enumerate(pcr.fwd_sites) if "Verif" in s.label)
    rev = next(i for i, s in enumerate(pcr.rev_sites) if "Verif" in s.label)
    pcr.fwd.setCurrentIndex(fwd)
    pcr.rev.setCurrentIndex(rev)
    settle(1.0)
    save(win, "primers_pcr")


def cas9():
    """The Spacer Finder on pFG051: every SpCas9 guide with its cut, score and other matching sites."""
    win = window(DEMO / "pFG051_Chlor.gb", tool="Spacer Finder")
    gene = next(f for f in win.editor.features if "Transposase" in (f.label or ""))
    start, end = min(a for a, _ in gene.spans), max(b for _, b in gene.spans)
    finder = [p.view for p in win.workspace.panes() if p.view_name == "Spacer Finder"][0]

    def guides_in_the_gene():
        win.editor.select(start, end)   # the finder searches the selection
        finder.reload()
        settle(12.0)                    # every guide is also checked against the 4.6 Mb background genome

    guides_in_the_gene()
    finder.table.setCurrentCell(0, 0)
    finder.add()                        # the first guide becomes a feature, drawn on the sequence above
    guides_in_the_gene()                # adding it selected the guide: search the whole gene again
    save(win, "cas9_guides")


def local_align():
    """Local Align: the araC gene of pBAD30 found on the E. coli chromosome, the hits drawn as a track."""
    from opendna_core import Document
    plasmid = Document(str(DEMO / "pBAD30_eilR.gb"))
    araC = next(f for f in plasmid.features() if f.label == "araC")
    query = "".join(plasmid.sequence(a, b) for a, b in araC.spans)
    win = window(DEMO / "BW25113_2026.gb", tool="Local Align")
    dialog = [p.view for p in win.workspace.panes() if p.view_name == "Local Align"][0]
    dialog.query.setPlainText(f">araC (pBAD30_eilR)\n{query}")
    dialog.run()
    settle(4.0)
    save(win, "local_align")


def reads_bam():
    """A fake BAM of pFG051 written like the unit tests write theirs (no bwa, no WSL): 150-bp reads at about 60x, 40 % of them
    from a clone that lacks 120 bases near 1,900, and a base changed in a third of the reads at one position."""
    import random
    sys.path.insert(0, str(ROOT / "tests"))
    import test_core_formats as T
    import runner
    bam = Path(os.environ["OPENDNA_DATA_DIR"]) / "clone.bam"
    from opendna_core import Document
    doc = Document(str(DEMO / "pFG051_Chlor.gb"))
    ref, rng, n = doc.sequence().upper(), random.Random(55), doc.length
    gap0, gap1, snp, reads = 1850, 1970, 700, []
    for i in range(n * 60 // 150):
        pos, left = rng.randrange(0, n - 150), 150
        deleted = rng.random() < 0.4
        if deleted and gap0 <= pos < gap1:
            pos = gap1
        if deleted and pos < gap0 < pos + left:      # the read runs into the missing stretch: M, D, M
            ops = [(gap0 - pos, "M"), (gap1 - gap0, "D"), (left - (gap0 - pos), "M")]
        else:
            ops = [(left, "M")]
        subs = {snp: "A" if ref[snp] != "A" else "G"} if rng.random() < 0.33 else {}
        reads.append(T.build_read(f"r{i}", pos, ops, ref, subs=subs, qual=[35] * 150, flag=16 * rng.randrange(2)))
    T.write_bam(bam, [(runner.safe_name(doc.name), n)], reads)
    return bam


def genome_tracks():
    """The Alignment view: simulated Illumina reads mapped on pFG051, their coverage and the mismatches they carry."""
    bam = reads_bam()
    win = window(DEMO / "pFG051_Chlor.gb", split("v", [6500, 3500], pane("Alignment"), pane("Sequence")))
    win.load_track(str(bam), quiet=True)
    settle(4.0)
    save(win, "reads_coverage")


PICTURES = {"editor": editor, "map": multi_genome, "gibson": gibson, "primers": primers, "cas9": cas9,
            "align": local_align, "tracks": genome_tracks}


if __name__ == "__main__":
    theme.apply("light", remember=False)
    OUT.mkdir(parents=True, exist_ok=True)
    for name in sys.argv[1:] or PICTURES:
        PICTURES[name]()
