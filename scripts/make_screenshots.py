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


def window(path, layout=None):
    win = M.MainWindow()
    win.resize(*SIZE)
    win.show()
    assert win.load(str(path)), path
    if layout is not None:
        win._build_workspace([{"layout": layout}])
    settle()
    return win


def save(win, name):
    settle()
    path = OUT / f"{name}.png"
    win.grab().save(str(path))
    win.close()
    print(f"{path.relative_to(ROOT)}")


# ------------------------------------------------------------------------------------------------------ the pictures

def editor():
    """The default layout on a real plasmid: sequence, map and feature table."""
    save(window(DEMO / "pBAD30_eilR.gb"), "editor")


def multi_genome():
    """Vibrio cholerae: two chromosomes in one GenBank file; the Map lists them, a click switches."""
    win = window(DEMO / "vibrio_cholerae_N16961.gb",
                 split("h", [6000, 4000], pane("Map"), split("v", [5000, 5000], pane("Sequence"), pane("Features"))))
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
    """The Primers pane and the PCR pane on pBAD30, with the plasmid's own Verif pBAD30_F / _R primers."""
    win = window(DEMO / "pBAD30_eilR.gb",
                 split("h", [4800, 5200], pane("Sequence"), split("v", [5200, 4800], pane("Primers"), pane("PCR"))))
    pcr = [p.view for p in win.workspace.panes() if p.view_name == "PCR"][0]
    fwd = next(i for i, s in enumerate(pcr.fwd_sites) if "Verif" in s.label)
    rev = next(i for i, s in enumerate(pcr.rev_sites) if "Verif" in s.label)
    pcr.fwd.setCurrentIndex(fwd)
    pcr.rev.setCurrentIndex(rev)
    settle(1.0)
    save(win, "primers_pcr")


def cas9():
    """The Spacer Finder on pFG051: every SpCas9 guide with its cut, score and other matching sites, one added to the sequence."""
    win = window(DEMO / "pFG051_Chlor.gb",
                 split("v", [4500, 5500], pane("Sequence"), pane("Spacer Finder")))
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
    win = window(DEMO / "BW25113_2026.gb",
                 split("v", [1800, 4200, 4000], pane("Sequence"), pane("Alignment"), pane("Local Align")))
    dialog = [p.view for p in win.workspace.panes() if p.view_name == "Local Align"][0]
    dialog.query.setPlainText(f">araC (pBAD30_eilR)\n{query}")
    dialog.run()
    settle(4.0)
    save(win, "local_align")


def reads_bam():
    """A BAM of simulated Illumina reads of pFG051, mapped with the program's own recipe (bwa in WSL): made once, kept in
    the scratch folder. Needs bwa and samtools, as Tools > Map Reads does."""
    import random
    sys.path.insert(0, str(ROOT / "tests"))
    import fake_reads as F
    import runner
    folder = Path(os.environ["OPENDNA_DATA_DIR"]) / "reads"
    bam = folder / "clone.bam"
    if bam.exists():
        return bam
    folder.mkdir(parents=True, exist_ok=True)
    from opendna_core import Document
    doc = Document(str(DEMO / "pFG051_Chlor.gb"))
    ref, rng = doc.sequence().upper(), random.Random(55)
    fasta = folder / "ref.fa"
    fasta.write_text(f">{runner.safe_name(doc.name)}\n{ref}\n")
    clone = F.clones(ref, rng)["cloneB_mixed_population"]        # a mixed colony: some reads carry the change, some do not
    haplotypes = [(share, F.apply(ref, events)) for share, events in clone]
    r1, r2 = folder / "R1.fastq", folder / "R2.fastq"
    with open(r1, "w") as f1, open(r2, "w") as f2:
        for n, (a, b, qa, qb) in enumerate(F.element(haplotypes, len(ref), rng)):
            f1.write(f"@r{n}/1\n{a}\n+\n{qa}\n")
            f2.write(f"@r{n}/2\n{b}\n+\n{qb}\n")
    F.map_reads("bwa_mem", fasta, r1, r2, folder / "clone")
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
