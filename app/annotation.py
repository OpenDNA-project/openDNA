"""Detect common features: find known parts (promoters, RBS, origins, markers, fluorescent
proteins, terminators…) in the sequence and add the ones you tick, optionally grouped
into genes (promoter + RBS + CDS + terminator on the same strand).

Parts come from the Features Database (parts_library.py: app/parts/*.tsv, bundled one file per
source, and your own features, the source "Local"); Sources there says which sources are used.
A part may carry the strength its source measured, in that source's units and host.
The search itself is in Rust (core/src/annotate.rs).

pLannotate (McGuffie & Barrick 2021, plannotate_tool.py), when it is installed, searches the same sequence with its own databases,
and its features are listed beside the Features Database's: one list, each feature saying who found it (the column "Found by":
openDNA, pLannotate, or both when they found the same thing in the same place) and in which database.
"""
from pathlib import Path
from types import SimpleNamespace

from PySide6.QtCore import QProcess, Qt
from PySide6.QtWidgets import (
    QAbstractItemView, QDialog, QHBoxLayout, QLabel, QMessageBox, QPlainTextEdit, QPushButton, QVBoxLayout,
)

from opendna_core import Document, find_parts

import genes
import genetic_code
import multigb
import plannotate_tool
import runner
import tool_settings
from editor import HIDDEN_KINDS, reverse_complement, site_location
import parts_library
from parts_library import add_parts, describe, load_parts
from tools import ToolDialog
from data_table import DataTable
from uikit import fit_columns, no_auto_default, plural, run_busy

TOOL = "Detect Features"
tool_settings.register(TOOL, [
    tool_settings.Param("min_dna", "Minimum DNA identity", 95, 70, 100, suffix=" %",
                        tip="A part is found where at least this share of its bases match (no gaps)."),
    tool_settings.Param("min_protein", "Minimum protein identity (CDS)", 90, 50, 100, suffix=" %",
                        tip="Coding parts are also searched as protein, to find codon-optimized versions."),
    tool_settings.Param("min_length", "Ignore parts shorter than", 15, 6, 100, suffix=" bp",
                        tip="Short parts (operators, RBS) match by chance in large sequences."),
    tool_settings.Param("group_genes", "Group promoter + RBS + CDS + terminator into genes", True),
    tool_settings.Param("gene_reach", "…when within", 300, 0, 5000, step=50, suffix=" bp of the CDS"),
    tool_settings.Param("use_plannotate", "Also search with pLannotate (when installed)", True,
                        tip="pLannotate, the Barrick lab's plasmid annotator, is a separate tool that you install once (Install pLannotate). "
                            "Its features are listed beside the Features Database's, each marked with who found it."),
])
KINDS_UPSTREAM = {"promoter": 0, "RBS": 1, "5'UTR": 1}
OWN, THEIRS, BOTH = "openDNA", "pLannotate", "openDNA + pLannotate"   # the "Found by" values
PROTEIN_DATABASES = {"fpbase", "swissprot"}   # pLannotate's searches that compare proteins (DIAMOND), not DNA (BLAST)


def feature_part(doc, f, source):
    """(name, kind, sequence in the feature's orientation, source) of a single-span feature."""
    s, e = f.spans[0]
    seq = doc.sequence(s, e)
    return (f.label, f.kind, reverse_complement(seq) if f.strand < 0 else seq, source)


def import_files(paths, rejected=None):
    """Every named single-span feature (not source/primer/gene) of these files into the
    user library. Returns how many parts were added; `rejected` receives (name, why) of those that cannot be stored."""
    new = []
    for path in paths:
        doc = Document(str(path))
        for f in doc.features():
            if f.kind in HIDDEN_KINDS or f.kind == genes.GENE or len(f.spans) != 1 or not f.label:
                continue
            if 6 <= f.spans[0][1] - f.spans[0][0] <= 10_000:
                new.append(feature_part(doc, f, Path(path).name))
    return add_parts(new, rejected)


def detect(editor, kinds=None, run=None):
    """Hits of library parts on the sequence: SimpleNamespace(name, kind, start, end, strand, identity,
    protein, source, collection, organism, strength, units, measured, known, covered), sorted by start. `known`: a
    feature of the same type already covers exactly the same bases (on the same strand for a CDS); `covered`: one covers most of them. `measured`: the
    strength its source reports, as text ("TS 98 (E. coli)"), or "". `kinds`: only parts of these types. `run(work)`: how the
    search itself (Rust, it lets go of Python's lock) is run, e.g. in a thread behind a busy dialog; by default, here."""
    min_len = tool_settings.value(TOOL, "min_length")
    parts = [p for p in load_parts() if len(p.sequence) >= min_len and (kinds is None or p.kind in kinds)]
    args = (editor.doc.sequence(), editor.doc.circular, [p.sequence for p in parts], [p.kind == "CDS" for p in parts],
            tool_settings.value(TOOL, "min_dna") / 100, tool_settings.value(TOOL, "min_protein") / 100, genetic_code.table_id())
    found = (run or (lambda work: work()))(lambda: find_parts(*args))
    hits = []
    for i, s, e, strand, identity, protein in found:
        p = parts[i]
        hits.append(SimpleNamespace(name=p.name, kind=p.kind, start=s, end=e, strand=strand, identity=identity,
                                    protein=protein, source=p.source, collection=p.collection, organism=p.organism,
                                    strength=p.strength, units=p.units, measured=describe(p), origin=OWN, also=None,
                                    known=already_there(editor, p.kind, s, e, strand),
                                    covered=already_there(editor, p.kind, s, e, strand, exact=False)))
    # Library parts on the same bases (two KanR variants, one promoter annotated on opposite
    # strands by two records, "rrnB" vs "5S RNA, T1, T2, rrnB"): keep one, preferring a
    # specific type over misc_feature, then the best identity, then the longest.
    hits.sort(key=lambda h: (h.kind == "misc_feature", -h.identity, -(h.end - h.start)))
    kept = []
    for h in hits:
        if not any(same_place(k, h) for k in kept):
            kept.append(h)
    kept.sort(key=lambda h: h.start)
    return kept


def plannotate_hits(editor, rows):
    """pLannotate's rows (plannotate_tool.annotate) as hits shaped like detect()'s, `origin` THEIRS, `source` "pLannotate: <its database>"
    (snapgene, fpbase, swissprot, Rfam). A fragment keeps pLannotate's own mark, "(fragment)", after the name."""
    return [SimpleNamespace(name=r["name"] + (" (fragment)" if r["fragment"] else ""), kind=r["kind"], start=r["start"], end=r["end"],
                            strand=r["strand"], identity=r["identity"], protein=r["database"] in PROTEIN_DATABASES,
                            source=f"{plannotate_tool.NAME}: {r['database']}", collection="", organism="", strength="", units="", measured="",
                            origin=THEIRS, also=None, database=r["database"], note=r["note"],
                            known=already_there(editor, r["kind"], r["start"], r["end"], r["strand"]),
                            covered=already_there(editor, r["kind"], r["start"], r["end"], r["strand"], exact=False)) for r in rows]


def same_thing(own, theirs):
    """The two tools found the same feature: mostly the same bases, the same type (or one of them only says misc_feature), and the
    same strand for a coding sequence (other parts are annotated on either strand)."""
    same_kind = own.kind == theirs.kind or "misc_feature" in (own.kind, theirs.kind)
    return same_kind and same_place(own, theirs) and (own.strand == theirs.strand or own.kind != "CDS")


def merge(own, theirs):
    """One list, by position, of what openDNA's Features Database found (`own`) and what pLannotate found (`theirs`). A feature both found in
    the same place is one entry, `origin` BOTH, with the Features Database's name and the other tool's hit in `also`; the rest keep
    their own `origin`. Nothing is dropped: two features of one tool that overlap stay two."""
    out, alone = list(own), []
    for t in theirs:
        mate = next((o for o in out if o.origin == OWN and same_thing(o, t)), None)
        if mate is None:
            alone.append(t)
        else:
            mate.origin, mate.also, mate.source = BOTH, t, f"{mate.source} + {t.source}"
    return sorted(out + alone, key=lambda h: (h.start, h.end))


def detect_all(editor, kinds=None, run=None):
    """(hits, state): detect()'s hits merged with pLannotate's (merge) when the setting is on and pLannotate is installed. `state.text` says
    what happened to pLannotate: not asked for (`state.name` "off"), not installed ("missing"), failed ("error", with the reason), or
    ran ("ran"). `kinds` restricts only the Features Database's search (pLannotate has no such option). `run`: as for detect()."""
    own = detect(editor, kinds, run)
    if not tool_settings.value(TOOL, "use_plannotate"):
        return own, SimpleNamespace(name="off", text="")
    try:
        rows = plannotate_tool.annotate(editor.doc.sequence(), editor.doc.circular, run)
    except plannotate_tool.NotInstalled as e:
        return own, SimpleNamespace(name="missing", text=str(e))
    except plannotate_tool.PlannotateError as e:
        return own, SimpleNamespace(name="error", text=str(e))
    found = plannotate_hits(editor, rows)
    if kinds is not None:
        found = [h for h in found if h.kind in kinds]
    return merge(own, found), SimpleNamespace(name="ran", text="")


def note_for(h):
    """The /note of a feature added from a hit: who found it, in which database, with what identity."""
    if h.origin == THEIRS:
        return f"detected by pLannotate (database {h.database}): {h.identity * 100:.0f} % identity" + (f"; {h.note}" if h.note else "")
    source = h.source[:-len(f" + {h.also.source}")] if h.origin == BOTH else h.source   # (the Features Database's part, before merge added the other)
    note = (f"detected: {h.identity * 100:.0f} % identity{' (protein)' if h.protein else ''} to {source or 'library part'}"
            f"{'; ' + h.measured if h.measured else ''}")
    if h.origin == BOTH:
        note += f"; also found by pLannotate as {h.also.name} (database {h.also.database}, {h.also.identity * 100:.0f} % identity)"
    return note


def overlap(a, b):
    return max(0, min(a.end, b.end) - max(a.start, b.start))


def same_place(a, b, share=0.9):
    """Mostly the same bases: the overlap covers `share` of each of them."""
    return overlap(a, b) >= share * max(a.end - a.start, b.end - b.start)


def already_there(editor, kind, s, e, strand, exact=True):
    """An existing feature covers exactly the bases of this hit (a part that only overlaps one is not annotated: it may be another
    part); with `exact` False, it covers most of them or lies mostly inside them (a shorter annotation of the same thing).
    Strand only matters for coding sequences: origins, terminators and the like are often annotated on either strand."""
    probe = SimpleNamespace(start=s, end=e)
    for f in editor.visible_features(s, min(e, editor.length), overlays=False):   # a tool's temporary ones are not annotations
        if f.kind == genes.GENE or (kind == "CDS" and f.strand != strand):
            continue
        fs, fe = genes.span_of([f])
        if (fs, fe) == (s, e) if exact else overlap(probe, SimpleNamespace(start=fs, end=fe)) >= 0.8 * min(e - s, fe - fs):
            return True
    return False


def gene_groups(hits, reach):
    """[(CDS hit, [part hits])]: for each CDS, the promoter / RBS upstream and terminator
    downstream on its strand within `reach` bp."""
    out = []
    for cds in (h for h in hits if h.kind == "CDS"):
        fwd = cds.strand > 0
        parts = []
        for h in hits:
            if h is cds or h.strand != cds.strand or h.kind == "CDS":
                continue
            upstream = (cds.start - h.end) if fwd else (h.start - cds.end)
            downstream = (h.start - cds.end) if fwd else (cds.start - h.end)
            if h.kind in KINDS_UPSTREAM and -10 <= upstream <= reach:
                parts.append(h)
            elif h.kind == "terminator" and -10 <= downstream <= reach:
                parts.append(h)
        if parts:
            out.append((cds, parts))
    return out


def features_for(editor, chosen, group):
    """Feature tuples for the ticked hits; with `group`, CDS + nearby parts become genes
    (shared locus_tag + a gene feature), all added in one undo step."""
    n = editor.length
    tag_of = {}
    extra = []
    if group:
        taken = set(editor.groups)
        for cds, parts in gene_groups(chosen, tool_settings.value(TOOL, "gene_reach")):
            if any(id(p) in tag_of for p in parts):
                continue  # a promoter shared by two CDS stays ungrouped (ambiguous)
            key = genes.unique_key(taken, cds.name)
            taken.add(key)
            members = [cds] + parts
            for m in members:
                tag_of[id(m)] = (key, cds.name)
            lo, hi = min(m.start for m in members), max(m.end for m in members)
            extra.append((genes.GENE, site_location(lo, hi, cds.strand, n),
                          [("gene", cds.name), ("locus_tag", key), ("label", cds.name)]))
    out = []
    for h in chosen:
        quals = [("label", h.name), ("note", note_for(h))]
        if id(h) in tag_of:
            key, gene = tag_of[id(h)]
            quals += [("gene", gene), ("locus_tag", key)]
        out.append((h.kind, site_location(h.start, h.end, h.strand, n), quals))
    return out + extra


class InstallDialog(QDialog):
    """Runs the install of pLannotate and/or the Terminal's tools in WSL / bash, its text shown as it comes. `ok`: it finished without error."""

    def __init__(self, parent, plannotate=True, terminal=False):
        super().__init__(parent)
        self.ok = False
        self.setWindowTitle("Install pLannotate" if not terminal else "Install the optional tools")
        self.resize(720, 360)
        self.log = QPlainTextEdit(readOnly=True)
        self.button = QPushButton("Stop", clicked=self._stop)
        layout = QVBoxLayout(self)
        layout.addWidget(self.log)
        layout.addWidget(self.button)
        self.process = QProcess(self)
        self.process.setProcessChannelMode(QProcess.MergedChannels)
        self.process.readyRead.connect(lambda: self.log.appendPlainText(bytes(self.process.readAll()).decode("utf-8", "replace").rstrip()))
        self.process.finished.connect(self._finished)
        self.process.errorOccurred.connect(lambda _: self._finished(-1) if self.process.state() == QProcess.NotRunning else None)
        argv = runner.command_line(["-lc", plannotate_tool.install_script(plannotate, terminal)])
        self.process.start(argv[0], argv[1:])

    def _stop(self):
        if self.process.state() == QProcess.NotRunning:
            self.accept()
        else:
            self.process.kill()

    def _finished(self, code, *_):
        self.log.appendPlainText("\nDone." if code == 0 else f"\nStopped or failed (exit code {code}). " + plannotate_tool.WSL_HELP)
        self.button.setText("Close")
        self.ok = code == 0

    def reject(self):   # Escape or the close box must not leave the install running unseen
        if self.process.state() != QProcess.NotRunning:
            self.process.kill()
            self.process.waitForFinished(3000)
        super().reject()


INSTALL_QUESTION = (
    "Install pLannotate, the plasmid annotator of the Barrick lab (McGuffie & Barrick, Nucleic Acids Res 2021, doi:10.1093/nar/gkab374), beside openDNA?\n\n"
    "It is a separate tool: it goes into its own conda environment ('opendna-plannotate') {where}, with its databases and the programs it "
    "uses (BLAST+, DIAMOND, Infernal). That downloads about 250 MB of databases and several hundred MB of programs, and takes a few "
    "minutes. openDNA then lists what it finds beside its own Features Database, each marked with who found it.")


class DetectDialog(ToolDialog):
    COLUMNS = ["Name", "Type", "Location", "Identity", "Strength", "Source", "Found by"]

    def __init__(self, parent, editor):
        super().__init__(parent)
        self.editor, self.hits = editor, []
        self.others = {}   # multi-record file: record name -> its Editor (the other chromosomes are searched too)
        self.columns = self.COLUMNS + (["Record"] if editor.multi_genbank else [])
        self.setWindowTitle("Detect Common Features")   # its name in the Tools menu
        self.resize(820, 480)
        self.table = DataTable(self.columns, filtered=["Type", "Source", "Found by"])
        self.table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        fit_columns(self.table, stretch=0)
        self.table.currentCellChanged.connect(lambda row, *_: self._select(self.table.data_row(row)))
        self.info = QLabel()
        self.credit = QLabel(openExternalLinks=True)   # who made pLannotate, with the paper's DOI, while its features are listed
        self.credit.setStyleSheet('color: gray')
        self.credit.hide()
        add = QPushButton("Add ticked features", clicked=self.add)
        settings = QPushButton("Settings…", clicked=lambda: tool_settings.SettingsDialog(self, TOOL).exec())
        import db_tabs
        library = QPushButton("Sources")   # a push button with a menu: its arrow has room of its own (a tool button drew it over the text)
        library.setToolTip("Which sources of the Features Database are looked for (Databases → Features Database lists them)")
        library.setMenu(db_tabs.sources_menu(library))
        self.install = QPushButton("Install pLannotate…", clicked=self.install_plannotate)   # shown only while it is not installed
        self.install.setToolTip("Also search with pLannotate (Barrick lab), a separate tool installed once in WSL; its features are listed here too")
        self.install.hide()
        row = QHBoxLayout()
        row.addWidget(self.info, 1)
        row.addWidget(self.install)
        row.addWidget(library)
        row.addWidget(settings)
        row.addWidget(add)
        layout = QVBoxLayout(self)
        layout.addWidget(self.table)
        layout.addWidget(self.credit)
        layout.addLayout(row)
        no_auto_default(self)   # Enter does not add every ticked feature
        # A bound method (not a lambda): Qt drops the connection when this dialog is deleted.
        tool_settings.notifier.changed.connect(self._settings_changed)
        parts_library.notifier.changed.connect(self._library_changed)
        self.refresh()

    def _library_changed(self):
        """A source was switched on / off, or features were added: what is found changes."""
        if self.isVisible():
            self.refresh()
        else:
            self.stale = True

    def _settings_changed(self, tool):
        if tool == TOOL:
            self.refresh()

    def refresh(self):
        # (the search, every part of every source against the sequence, in a thread behind a busy dialog if it takes a while)
        busy = lambda work: run_busy(self, "Detecting features…", work)
        self.hits, plannotate = detect_all(self.editor, run=busy)   # (this sequence only: pLannotate is not run on the file's other records)
        self.install.setVisible(plannotate.name == "missing")
        for h in self.hits:
            h.record = None
        if self.editor.multi_genbank:   # the same tool, run on each other record in turn
            self.others = dict(multigb.sibling_editors(self.editor))
            for name, other in self.others.items():
                found = detect(other, run=busy)
                for h in found:
                    h.record = name
                self.hits += found
        with self.table.refilling(len(self.hits)):
            for i, h in enumerate(self.hits):
                length = (self.others[h.record] if h.record else self.editor).length
                loc = f"{h.start + 1:,}–{(h.end - 1) % length + 1:,} {'→' if h.strand > 0 else '←'}"
                ident = f"{h.identity * 100:.0f} %{' (protein)' if h.protein else ''}"
                where = [h.record or self.editor.fasta_record] if self.editor.multi_genbank else []
                keys = {2: (h.record or "", h.start, h.end), 3: h.identity}   # Location and Identity sort by number; the rest by their text
                for c, v in enumerate([h.name, h.kind, loc, ident, h.measured, h.source, h.origin] + where):
                    item = self.table.put(i, c, v, key=keys.get(c), right=c == 3)
                    if c == 0:
                        item.setFlags(item.flags() | Qt.ItemIsUserCheckable)
                        item.setCheckState(Qt.Unchecked if h.known else Qt.Checked)
                        if h.known:
                            item.setToolTip("Already annotated")
                        elif h.origin == BOTH:
                            item.setToolTip(f"Found by both. pLannotate calls it {h.also.name} ({h.also.source.split(': ')[-1]}, {h.also.identity * 100:.0f} % identity).")
        new = sum(not h.known for h in self.hits)
        by = {o: sum(h.origin == o for h in self.hits) for o in (OWN, THEIRS, BOTH)}
        who = f" ({by[OWN]} by openDNA, {by[THEIRS]} by pLannotate, {by[BOTH]} by both)" if plannotate.name == "ran" else ""
        problem = {"missing": " pLannotate is not installed.", "error": f" pLannotate: {plannotate.text}"}.get(plannotate.name, "")
        self.info.setText(f"{plural(len(self.hits), 'feature')} found{who}, {new} not annotated yet (ticked).{problem}")
        self.credit.setVisible(plannotate.name == "ran")
        self.credit.setText(f"pLannotate features: {plannotate_tool.CITATION} <a href='https://doi.org/{plannotate_tool.DOI}'>doi:{plannotate_tool.DOI}</a>")

    def install_plannotate(self):
        where = "inside WSL" if runner.use_wsl() else "(conda)"
        if QMessageBox.question(self, "Install pLannotate", INSTALL_QUESTION.format(where=where), QMessageBox.Yes | QMessageBox.No,
                                QMessageBox.No) != QMessageBox.Yes:
            return
        dialog = InstallDialog(self)
        dialog.exec()
        plannotate_tool.forget()
        self.refresh()

    def _select(self, i):
        if 0 <= i < len(self.hits) and not self.hits[i].record:   # (a hit on another record is not on this sequence)
            h = self.hits[i]
            self.editor.select(h.start, min(h.end, self.editor.length))
            self.editor.reveal.emit(h.start)

    def ticked(self):
        rows = range(self.table.rowCount())
        return [self.hits[self.table.data_row(r)] for r in rows if self.table.item(r, 0).checkState() == Qt.Checked]

    def add(self):
        chosen = self.ticked()
        here = [h for h in chosen if not h.record]
        away = {}
        for h in chosen:
            if h.record:
                away.setdefault(h.record, []).append(h)
        count = sum(map(len, away.values()))
        if away and QMessageBox.question(   # No by default: what is written there cannot be undone
                self, "Other sequences of the file",
                f"{plural(count, 'ticked feature')} {'is' if count == 1 else 'are'} on other sequences of {self.editor.path.name} "
                f"({', '.join(away)}). They are written straight into the file (no undo). Continue?",
                QMessageBox.Yes | QMessageBox.No, QMessageBox.No) != QMessageBox.Yes:
            return
        group = tool_settings.value(TOOL, "group_genes")
        if here:
            self.editor.add_features(features_for(self.editor, here, group))
        for name, hits in away.items():
            other = self.others[name]
            other.add_features(features_for(other, hits, group))
            multigb.save_record_of(other)
        if away:
            multigb.forget_siblings(self.editor)
        if chosen:
            self.refresh()
