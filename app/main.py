"""openDNA main window: one file per window, shown in tabs of split panes where any
view or tool can go (see workspace.py).

Run:  python app/main.py [file.gb | file.dna | genome.fa]
"""
import os
import re
import sys
from pathlib import Path
from types import SimpleNamespace

from PySide6 import __version__ as pyside_version
from PySide6.QtCore import QEvent, QSettings, Qt, QTimer, qVersion
from PySide6.QtGui import QAction, QActionGroup, QGuiApplication, QIcon, QKeySequence
from PySide6.QtWidgets import (
    QApplication, QCheckBox, QComboBox, QDialogButtonBox, QDockWidget, QFileDialog, QFormLayout, QInputDialog, QLabel,
    QLineEdit, QMainWindow, QMenu, QMessageBox, QPlainTextEdit, QDialog, QVBoxLayout, QWidget,
)

from opendna_core import Bam, Document, Fasta, Track

from display_panel import DisplayPanel
from history_panel import HistoryPanel, side_panel_width
from editor import Editor, clean_pasted
from shortcuts import ShortcutsDialog
from sponsors import SponsorsDialog, show_at_startup
from feature_dialog import edit_feature
from feature_table import FeatureTable, delete_feature, features_inside, group_as_gene
from map_view import MapView
from editor import primer_sequence
from primer_db import copy_old_data, data_dir, db
from perf import PerfWindow
from primer_panel import PrimerDialog, PrimerPanel, annealed, melting_temp
from assembly_panel import AssemblyPanel
from citations import CitationsDialog
from editor import remember_hidden_primers
from sourcefiles import add_recent, pinned_files, recent_files
from local_align import LocalAlignDialog
from codon_tools import CompareDialog, RecodeDialog
from folder_search import FolderSearchDialog
from tools import DigestDialog, OrfDialog, PcrDialog
import updater
import exports
from enzymes import EnzymeFilter
import icons
from checkup import TOOL as CHECKUP
import image_export
import ncbi
from ncbi_dialog import GenomeDialog, NcbiDialog
from primer_panel import export_for_ordering, import_primers
from primer_design import design_for_selection
from structure_dialog import StructureDialog
from search_bar import SearchBar
from genome_view import GenomeView
import feedback
import file_assoc
import project
import runner
import wip
import db_status
import db_tabs
from reference_picker import ReferencePicker
import theme
import track_check
from runner import RunDialog
from sequence_view import SequenceView
from workspace import DEFAULT as DEFAULT_LAYOUT, Workspace, pane, split
from annotation import DetectDialog, import_files
from parts_library import user_library
from cas9 import Cas9Dialog
from repeats import RepeatsDialog
import tool_settings
import codon_dialog
import codon_usage
from code_switch import CodeSwitch, shown as code_switch_shown
from tool_settings import SettingsDialog

import genetic_code
import fileio
import shutil
from file_history import FileHistory
import add_record
from add_record import AddRecordsDialog
import multigb
from uikit import plural, run_busy

tool_settings.register_page(genetic_code.TOOL, codon_dialog.UsageBox)   # the codon usage line of the Genetic code page

FASTA_EXT = {".fa", ".fasta", ".fna", ".fas", ".ffn", ".frn", ".mfa"}
OPEN_FILTER = ("Sequence files (*.gb *.gbk *.gbff *.genbank *.dna *.fa *.fasta *.fna *.fas *.odp);;"
               "Open DNA Project (*.odp);;All files (*)")
TRACK_FILTER = ("Tracks (*.gff *.gff3 *.gtf *.bed *.vcf *.vcf.gz *.bam);;Alignments (*.bam);;"
                "Variants (*.vcf *.vcf.gz);;All files (*)")
ANNOTATION_FILTER = "Annotation tracks (*.gff *.gff3 *.gtf *.bed);;All files (*)"
VCF_FILTER = "Variants (*.vcf *.vcf.gz);;All files (*)"
BAM_FILTER = "Aligned reads (*.bam);;All files (*)"
AUTOSAVE_MS = 60_000


def autosave_path(path, record=None):
    """Recovery copy next to the file: plasmid.dna -> plasmid.dna.autosave.gb
    (FASTA: one per sequence, genome.fa.chr2.autosave.gb)."""
    path = Path(path)
    return path.with_name(path.name + (f".{record}" if record else "") + ".autosave.gb")


def is_fasta(path):
    path = Path(path)
    if path.suffix.lower() in FASTA_EXT:
        return True
    try:
        with open(path, "rb") as fh:
            return fh.read(64).lstrip()[:1] == b">"
    except OSError:
        return False


def match_contig(seqs, *names):
    """The sequence of a track file (`seqs`) that is the open one, or None: the same name as one of `names` (the document's
    name, the FASTA record), else the same name once "chr" and the case are ignored ("chr1" is "1") if only one fits."""
    for n in names:
        if n and n in seqs:
            return n
    plain = lambda x: x.lower().removeprefix("chr")
    wanted = {plain(n) for n in names if n}
    close = [q for q in seqs if plain(q) in wanted]
    return close[0] if len(close) == 1 else None


def parse_range(text, length, circular):
    """'101..200', '101-200', '101 to 200', '1,000..2,000' (1-based, both ends included) -> (start, end) for Editor.select;
    on a circular sequence an end below the start runs through the origin; one number N -> the cursor after base N.
    Raises ValueError with a sentence to show."""
    words = re.sub(r"(?<=\d),(?=\d{3}\b)", "", text.strip())   # thousands separators
    numbers = [int(n) for n in re.findall(r"\d+", words)]
    if not numbers or len(numbers) > 2 or re.search(r"[^\d\s.\-–—to]", words, re.I):
        raise ValueError(f"Type a stretch as 101..200 (or one base number), from 1 to {length:,}.")
    if any(n < 1 or n > length for n in numbers) and not (len(numbers) == 1 and numbers[0] == 0):
        raise ValueError(f"The sequence has {length:,} bp: use numbers from 1 to {length:,}.")
    if len(numbers) == 1:
        return numbers[0], numbers[0]
    a, b = numbers
    if b < a and circular:
        return a - 1, b + length   # through the origin
    a, b = sorted((a, b))
    return a - 1, b


def other_records_note(editor, picker_count):
    """' (only chr3: the other 23 sequences of genome.fa are not included)' when a one-sequence copy leaves records behind."""
    if editor.fasta_record and picker_count > 1:
        return (f" (only {editor.fasta_record}: the other {plural(picker_count - 1, 'sequence')} of {editor.path.name} "
                "are not included)")
    return ""


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.editor = None
        self.workspace = None
        self.file_history = FileHistory()   # the order of the edits across the records of a multi-record file
        self.record_docs = {}        # records of the open multi-record file that are not shown, kept in memory with their edits
        self.autosave_file = None    # recovery copy we wrote (deleted after save/close)
        self.autosaved_version = None
        self.hint = QLabel("Open a GenBank or .dna file (Ctrl+O)", alignment=Qt.AlignCenter)
        self.setCentralWidget(self.hint)
        # The Assembly panel is a tab (opened from Tools > Assembly only); `tool_host` holds it while no file is open.
        self.assembly, self.tool_host = None, None
        self.tool_tabs = {}   # title -> widget of the tabs that are tools (Assembly, the databases): kept when their tab closes
        self.status = QLabel()
        # permanent: the selection readout stays while a message (saved, searched...) is shown, instead of being hidden by it
        self.statusBar().addPermanentWidget(self.status)
        # Display filter in a side panel (its content is made per document, in show_document).
        self.display_dock = QDockWidget("Display", self, objectName="display")
        self.display_dock.setWidget(QLabel())
        self._fix_to_the_right(self.display_dock)
        self.addDockWidget(Qt.RightDockWidgetArea, self.display_dock)
        self.display_dock.hide()  # folded away until the Filters button opens it
        # History of changes: same place as the filters, one or the other.
        self.history_dock = QDockWidget("History", self, objectName="history")
        self.history_dock.setWidget(QLabel())
        self._fix_to_the_right(self.history_dock)
        self.addDockWidget(Qt.RightDockWidgetArea, self.history_dock)
        self.history_dock.hide()
        self.display_dock.visibilityChanged.connect(lambda on: on and self.history_dock.hide())
        self.history_dock.visibilityChanged.connect(lambda on: on and self.display_dock.hide())
        self._menus()
        self._order_menus()
        for menu in self.menuBar().findChildren(QMenu):   # the tooltips written on menu items (what an item does) are shown
            menu.setToolTipsVisible(True)
        self.setAcceptDrops(True)
        self._update_title()
        self._restore_window()
        self.autosave_timer = QTimer(self, interval=AUTOSAVE_MS, timeout=self.autosave)
        self.autosave_timer.start()

    # ---- what the window remembers between sessions: its size and place, the View switches. Not the tabs (a file always opens
    # with the standard ones, unless it is a project), not which side panel was open, and no side panel width (fixed).
    VIEW_SWITCHES = ("enzymes_action", "translation_action", "complement_action", "collapse_action", "cds_only_action")

    @staticmethod
    def _window_settings():
        return QSettings(str(data_dir() / "window.ini"), QSettings.IniFormat)

    def _restore_window(self):
        s = self._window_settings()
        geometry = s.value("geometry")
        if geometry is None or not self.restoreGeometry(geometry):
            # first start: most of the screen, centred (a fixed 1400 x 850 was taller than a 1366 x 768 laptop screen)
            screen = (self.screen() or QGuiApplication.primaryScreen()).availableGeometry()
            w, h = min(1400, int(screen.width() * 0.85)), min(850, int(screen.height() * 0.85))
            self.resize(w, h)
            self.move(screen.center().x() - w // 2, screen.center().y() - h // 2)
        for name in self.VIEW_SWITCHES:
            a = getattr(self, name)
            saved = s.value(f"view/{name}")
            if saved is not None:
                a.setChecked(str(saved).lower() == "true")
            a.toggled.connect(lambda on, name=name: self._window_settings().setValue(f"view/{name}", on))

    def _save_window(self):
        s = self._window_settings()
        s.setValue("geometry", self.saveGeometry())   # size, place and maximised (restoreGeometry keeps it on a screen)
        s.remove("panel_width")   # (the side panels were once as wide as dragged: they all have one fixed width now)

    @staticmethod
    def _fix_to_the_right(dock):
        """A side panel that is always on the right: no title bar (no heading, no float / close buttons) and it cannot be moved; the
        toolbar's button (or its shortcut) opens and closes it. (Closable stays: with no feature at all Qt's toggle action cannot
        show the panel.) Every side panel has the same fixed width, the History panel's (side_panel_width): it cannot be dragged
        wider, and no content (the enzyme criteria unfolded, a long feature type...) can push it to half the window."""
        dock.setFeatures(QDockWidget.DockWidgetClosable)
        dock.setTitleBarWidget(QWidget())
        dock.setFixedWidth(side_panel_width())

    MENU_ORDER = ("&File", "&Edit", "&View", "&Import", "E&xport", "Fea&ture", "&Primers", "Data&bases", "T&ools", "&Help")

    def _order_menus(self):
        """The menu bar in the usual order: File, Edit, View first (where every program has them), then this program's own
        menus, Help last. The menus are built in the order their code needs; they are only moved here."""
        bar = self.menuBar()
        by_text = {a.text(): a for a in bar.actions()}
        for text in self.MENU_ORDER:
            a = by_text.pop(text, None)
            if a is not None:
                bar.removeAction(a)
                bar.addAction(a)
        help_action = bar.actions()[-1]
        for a in by_text.values():   # (a menu not named above goes just before Help)
            bar.removeAction(a)
            bar.insertAction(help_action, a)

    def _menus(self):
        self.doc_actions = []  # enabled once a document is open

        def action(menu, text, shortcut, slot, checkable=False, checked=False):
            a = QAction(text, self, shortcut=shortcut or QKeySequence(), triggered=slot,
                        checkable=checkable, checked=checked)
            menu.addAction(a)
            self.doc_actions.append(a)
            return a

        file = self.menuBar().addMenu("&File")
        file.addAction(QAction("&New…", self, shortcut=QKeySequence.New, triggered=self.new_file))
        file.addAction(QAction("&Open…", self, shortcut=QKeySequence.Open, triggered=self.choose_file))
        file.aboutToShow.connect(lambda: self.remove_tracks_action.setEnabled(
            self.editor is not None and bool(self.editor.tracks or self.editor.alignments)))   # greyed with nothing loaded
        self.recent_menu = file.addMenu("Open &Recent")   # the start page's list: pinned first, then newest; filled when opened
        self.recent_menu.aboutToShow.connect(self._fill_recent_menu)
        self._fill_recent_menu()
        action(file, "&Save", QKeySequence.Save, self.save)
        save_as = action(file, "Save &As…", None, self.save_as)
        # QKeySequence.SaveAs has no key on Windows: Ctrl+Shift+S everywhere (once, where the platform's own is the same)
        save_as.setShortcuts(list(dict.fromkeys([*QKeySequence.keyBindings(QKeySequence.SaveAs), QKeySequence("Ctrl+Shift+S")])))
        file.addSeparator()
        self.remove_tracks_action = action(file, "Remove &Tracks…", None, self.remove_tracks)
        file.addSeparator()
        file.addAction(QAction("&Close Window", self, shortcut=QKeySequence.Close, triggered=self.close))
        file.addAction(QAction("&Quit", self, shortcut=QKeySequence.Quit, triggered=QApplication.closeAllWindows))   # every window, each asking about its own unsaved changes

        # Import: by what the file holds. Tracks and reads need an open sequence (greyed until then).
        imp = self.menuBar().addMenu("&Import")
        dna = imp.addMenu("&DNA")
        dna.addAction(QAction("&Sequence File (GenBank, .dna, FASTA, Project)…", self, triggered=self.choose_file))
        action(dna, "Sequence File as New &Chromosome…", None, lambda: self.import_as_chromosome())
        dna.addAction(QAction("From &NCBI (Nucleotide)…", self, triggered=self.import_ncbi))
        dna.addAction(QAction("&Background Genome from NCBI…", self, triggered=self.import_genome))
        features = imp.addMenu("&Features and Tracks")
        action(features, "&Annotation Track (GFF3, GTF, BED)…", None, lambda: self.load_tracks(ANNOTATION_FILTER))
        action(features, "&Variants (VCF)…", None, lambda: self.load_tracks(VCF_FILTER))
        features.addAction(QAction("Import &Features to Features Database…", self, triggered=self.import_parts))
        proteins = imp.addMenu("&Proteins")
        proteins.addAction(wip.mark_action(QAction("Protein &FASTA (3D Structure)…", self, triggered=self.import_protein), "protein"))
        proteins.addAction(wip.mark_action(QAction("From &NCBI (Protein)…", self, triggered=lambda: self.import_ncbi("Protein")),
                                           "protein"))
        alignments = imp.addMenu("&Alignments")
        action(alignments, "Aligned &Reads (BAM)…", None, lambda: self.load_tracks(BAM_FILTER))
        wip.mark_action(action(alignments, "&Map Reads (FASTQ) to This Sequence…", None, self.map_reads), "mapping")
        primer_import = imp.addMenu("P&rimers")
        primer_import.addAction(QAction("Primer &List (TSV, CSV, FASTA)…", self, triggered=self.import_primer_list))

        # Export: each kind greyed while there is nothing of it (see _update_export_menu).
        exp = self.menuBar().addMenu("E&xport")
        self.export_whole = exp.addMenu("&Whole File")
        self.export_whole.addAction(QAction("&GenBank (.gb)…", self, triggered=self.export_genbank))
        self.export_whole.addAction(QAction("Open DNA &Project (.odp: tracks, alignments, display)…", self,
                                            triggered=self.export_project))
        self.export_whole.addAction(QAction("NCBI &Submission (.fsa + .tbl for table2asn)…", self,
                                            triggered=self.export_submission))
        self.export_dna = exp.addMenu("&DNA")
        self.export_dna.addAction(QAction("&FASTA…", self, triggered=lambda: self.export_fasta(False)))
        self.export_selection = QAction("FASTA of the &Selection…", self, triggered=lambda: self.export_fasta(True))
        self.export_dna.addAction(self.export_selection)
        self.export_dna.addAction(QAction("&Plain Text (bases only)…", self, triggered=self.export_plain))
        self.export_features = exp.addMenu("&Features")
        self.export_features.addAction(QAction("&GFF3…", self, triggered=lambda: self.export_text(
            "Export features", "gff3", "GFF3 (*.gff3 *.gff)", exports.gff3_text)))
        self.export_features.addAction(QAction("&BED…", self, triggered=lambda: self.export_text(
            "Export features", "bed", "BED (*.bed)", exports.bed_text)))
        self.export_features.addAction(QAction("NCBI Feature &Table (.tbl)…", self, triggered=lambda: self.export_text(
            "Export feature table", "tbl", "NCBI feature table (*.tbl)", exports.tbl_text)))
        self.export_protein = exp.addMenu("&Protein")
        self.export_cds = QAction("&CDS Translations (FASTA)…", self, triggered=lambda: self.export_text(
            "Export proteins", "faa", "Protein FASTA (*.faa *.fasta)", exports.proteins_fasta))
        self.export_protein.addAction(self.export_cds)
        self.export_translation = QAction("&Translation of the Selection (FASTA)…", self, triggered=self.export_translated)
        self.export_protein.addAction(self.export_translation)
        self.export_alignments = exp.addMenu("&Alignments (SAM)")   # one entry per loaded BAM, made when shown
        self.export_primers = exp.addMenu("P&rimers")
        self.export_primers.addAction(QAction("Primer &Database (TSV, CSV, FASTA)…", self, triggered=self.export_primer_db))
        self.export_primers.addAction(QAction("For &Ordering (IDT plate, CSV)…", self, triggered=self.export_primers_for_ordering))
        exp.addSeparator()
        self.export_image = QAction("&Image of the Current Tab (SVG, PNG, TIFF)…", self, triggered=self.export_tab_image)
        exp.addAction(self.export_image)
        exp.aboutToShow.connect(self._update_export_menu)

        edit = self.menuBar().addMenu("&Edit")
        self.undo_action = action(edit, "&Undo", QKeySequence.Undo, self.undo)
        self.redo_action = action(edit, "&Redo", QKeySequence.Redo, self.redo)
        edit.addSeparator()
        # Cut / Copy / Paste / Select All act on whatever has the keyboard focus (a text field, a table, the sequence). Their keys
        # are only shown here (menu-wide shortcut context): the focused widget keeps handling them itself, as before.
        for text, key, verb in (("Cu&t", QKeySequence.Cut, "cut"), ("&Copy", QKeySequence.Copy, "copy"),
                                ("&Paste", QKeySequence.Paste, "paste"), ("Select &All", QKeySequence.SelectAll, "select_all")):
            a = QAction(text, self, shortcut=key, triggered=lambda _=False, verb=verb: self.clipboard_action(verb))
            a.setShortcutContext(Qt.WidgetShortcut)
            edit.addAction(a)
        action(edit, "Copy Re&verse Complement", "Ctrl+Shift+C", self.copy_reverse_complement)
        action(edit, "Reverse-Complement &Whole Sequence", None, lambda: self.editor.reverse_complement())
        self.circular_action = action(edit, "C&ircular Sequence", None, self.set_circular, True)
        self.lock_action = action(edit, "&Lock DNA", "Ctrl+Shift+L", self.set_dna_locked, True)
        edit.addSeparator()
        action(edit, "&Find…", QKeySequence.Find, self.focus_find)
        # F3 / Shift+F3 on every platform: the standard Find Next of macOS and GNOME is Ctrl+G, which is Go to Position here
        # (and Find Previous would be Ctrl+Shift+G, the Alignment view): two actions on one key, neither would fire
        action(edit, "Find &Next", "F3", self.find_next)
        action(edit, "Find Previou&s", "Shift+F3", self.find_previous)
        action(edit, "&Go to Position…", "Ctrl+G", self.go_to)

        feature = self.menuBar().addMenu("Fea&ture")
        action(feature, "&New Feature from Selection…", "Ctrl+T", lambda: edit_feature(self, self.editor))
        action(feature, "&Edit Feature…", "Ctrl+E", self.edit_current_feature)
        action(feature, "&Delete Feature", None, self.delete_current_feature)
        feature.addSeparator()
        action(feature, "&Group as Gene…", "Ctrl+Shift+E", self.group_as_gene)
        feature.addSeparator()
        action(feature, "&Recolor Annotations by Type", None, lambda: self.statusBar().showMessage(
            f"{plural(self.editor.recolor_by_type(), 'annotation')} recolored (Ctrl+Z undoes it).", 5000))
        feature.addSeparator()
        # The /translation qualifier of CDSs: refreshed when a GenBank file is saved (never by autosave), or on demand.
        # The switch is not a document action: it is a persistent preference and can be changed with no file open.
        self.keep_translation_action = QAction("&Keep /translation up to Date", self, checkable=True,
                                               checked=genetic_code.keep_translation(), triggered=genetic_code.set_keep_translation)
        self.keep_translation_action.setToolTip("When a GenBank file is saved, recompute the /translation of every CDS that has one "
                                                "from its own bases, /codon_start and /transl_table.")
        feature.addAction(self.keep_translation_action)
        action(feature, "Update /&translation Now", None, self.update_translations)

        primers = self.menuBar().addMenu("&Primers")
        action(primers, "&New Primer…", "Ctrl+Shift+P", lambda: PrimerDialog(self, self.editor).exec())
        primers.addAction(QAction("Primer &Database…", self, triggered=lambda: db_tabs.open_tab(self, db_tabs.PRIMERS)))
        action(primers, "&Options…", None, self.primer_options)
        primers.addSeparator()
        action(primers, "D&esign Primers for Selection…", None, lambda: design_for_selection(self, self.editor))
        action(primers, "&PCR Product…", None, self.pcr_product)

        # Each database opens as a tab of this window (db_tabs.py); none needs an open file.
        databases = self.menuBar().addMenu("Data&bases")
        for label, key in db_tabs.MENU:
            databases.addAction(QAction(label, self, shortcut="Ctrl+Shift+D" if key == db_tabs.PRIMERS else QKeySequence(),
                                        triggered=lambda _=False, key=key: db_tabs.open_tab(self, key)))
        databases.addSeparator()
        status = databases.addMenu("Database &Status")   # where each database is, how many entries, what is wrong in it
        for label, check in db_status.DATABASES:
            status.addAction(QAction(label, self, triggered=lambda _=False, check=check: db_status.show(self, check).exec()))

        tools = self.menuBar().addMenu("T&ools")
        action(tools, "S&pacer Finder (Cas9 / Cas12a)…", None, lambda: Cas9Dialog(self, self.editor).exec())
        action(tools, "&Repeats and Dot Plot…", None, lambda: RepeatsDialog(self, self.editor).exec())
        tools.addAction(QAction("&Assembly", self, triggered=self.show_assembly))  # needs no open file
        action(tools, "&Local Align…", "Ctrl+Shift+B", lambda: LocalAlignDialog(self, self.editor).exec())
        action(tools, "Align &Clipboard…", None, lambda: LocalAlignDialog.align_clipboard(self, self.editor))
        action(tools, "Folder S&earch…", None, lambda: FolderSearchDialog(self, self.editor).exec())
        action(tools, "Restriction &Digest…", None, lambda: DigestDialog(self, self.editor).exec())
        tools.addSeparator()
        action(tools, "Detect Common &Features…", "Ctrl+Shift+A", lambda: DetectDialog(self, self.editor).exec())
        action(tools, "Detect &ORFs…", None, lambda: OrfDialog(self, self.editor).exec())
        action(tools, "&Recode CDS for a Codon Usage Table…", None, lambda: RecodeDialog(self, self.editor).exec())
        action(tools, "Compare CDS Under Another &Genetic Code…", None, lambda: CompareDialog(self, self.editor).exec())
        action(tools, "Plasmid C&heckup (expression of each CDS)…", None, self.open_checkup)
        tools.addAction(QAction("&Import Features to Features Database…", self, triggered=self.import_parts))
        tools.addSeparator()
        wip.mark_action(action(tools, "&Map Reads to This Sequence…", None, self.map_reads), "mapping")
        tools.addAction(wip.mark_action(QAction("Check Mapping &Tools…", self, triggered=lambda: RunDialog(
            self, "Mapping tools", [runner.CHECK]).exec()), "mapping"))
        tools.addAction(wip.mark_action(QAction("I&nstall Mapping Tools (conda)…", self, triggered=self.install_tools), "mapping"))
        tools.addSeparator()
        tools.addAction(QAction("&Settings…", self, triggered=lambda: self.open_settings()))

        view = self.menuBar().addMenu("&View")
        self.enzymes_action = action(view, "&Restriction Sites", None,
                                     lambda on: self.editor.set_display(enzymes=on), True, False)  # off by default
        self.translation_action = action(view, "3-Frame &Translation", "Ctrl+Shift+T",
                                         lambda on: self.editor.set_display(translation=on), True)
        self.complement_action = action(view, "&Complementary Strand", None,
                                        lambda on: self.editor.set_display(complement=on), True, True)   # off: the top strand alone
        self.genome_action = action(view, "&Alignment View", "Ctrl+Shift+G",
                                    self.set_genome_visible, True)
        self.collapse_action = action(view, "Genes as One Arrow on &Map", None,
                                      lambda on: self.editor.set_display(collapse_genes=on), True, True)
        self.cds_only_action = action(view, "Hide Gene Arrow when its CDS is A&lone", None,
                                      lambda on: self.editor.set_display(cds_only=on), True, True)
        self.cds_only_action.setToolTip("A gene whose only part is one CDS (the usual NCBI gene + CDS pair) is shown as that CDS only, in every view")
        filter_action = self.display_dock.toggleViewAction()  # checked while the panel is shown
        filter_action.setText("&Display Filter (Side Panel)")
        filter_action.setShortcut("Ctrl+Shift+F")
        view.addAction(filter_action)
        view.addSeparator()
        themes = QActionGroup(self)
        self.theme_menu = theme_menu = view.addMenu("Th&eme")
        self.theme_actions = {}
        for key, (label, _) in theme.THEMES.items():
            a = QAction(label, self, checkable=True, checked=key == theme.saved())
            a.triggered.connect(lambda _checked=False, key=key: self.set_theme(key))
            themes.addAction(a)
            theme_menu.addAction(a)
            self.theme_actions[key] = a
        action(view, "Re&set Tabs", None, self.reset_layout)

        help_menu = self.menuBar().addMenu("&Help")
        help_menu.addAction(QAction("&Keyboard Shortcuts…", self, shortcut="F1",
                                    triggered=lambda: ShortcutsDialog(self).exec()))
        help_menu.addAction(QAction("&Citations…", self, triggered=lambda: CitationsDialog(self).exec()))
        help_menu.addAction(QAction("&About and Licence…", self, triggered=self.about))
        help_menu.addAction(QAction("&Welcome…", self, triggered=lambda: SponsorsDialog(self).exec()))
        help_menu.addSeparator()
        help_menu.addAction(QAction("&Report a Problem…", self, triggered=lambda: feedback.open_report("bug")))
        help_menu.addAction(QAction("Request a &Feature…", self, triggered=lambda: feedback.open_feature()))
        # what helps find a problem, together and out of the way: the log, the diagnostics, the self-test, the timings
        trouble = help_menu.addMenu("&Troubleshooting")
        trouble.addAction(QAction("Open &Log Folder", self, triggered=lambda: feedback.open_log_folder()))
        trouble.addAction(QAction("Copy &Diagnostics", self, triggered=self.copy_diagnostics))
        trouble.addSeparator()
        trouble.addAction(QAction("&Self-Test…", self, triggered=self.show_self_test))
        trouble.addAction(QAction("&Performance Timings…", self, triggered=self.show_timings))
        help_menu.addSeparator()
        help_menu.addAction(QAction("Check for &Updates…", self, triggered=lambda: updater.check_now(self, **self.update_hooks())))
        help_menu.addAction(QAction("File Ass&ociations…", self, triggered=lambda: file_assoc.show_manage_dialog(self)))

        self.search_bar = SearchBar(self)
        self.find_box = self.search_bar.box
        toolbar = self.addToolBar("Find")
        toolbar.setMovable(False)
        # First, so it never falls into the overflow (») of a narrow window. Pressed: the bases cannot change (typing, paste,
        # delete, undo of a base edit...); features and the rest still can.
        self.lock_action.setToolTip("Lock the DNA: the bases cannot be changed (features, primers and the rest still can)")
        self.lock_action.setIconText("DNA")
        self._show_lock(False)
        toolbar.addAction(self.lock_action)
        toolbar.widgetForAction(self.lock_action).setToolButtonStyle(Qt.ToolButtonTextBesideIcon)   # a padlock, then "DNA"
        # The sequence's own two switches side by side, first: locked or not, circular or linear. (Anything that has to fall
        # into the overflow (») of a narrow window comes from the end: the side panels, which also have their shortcuts.)
        self.circular_action.setIconText("Circular")  # pressed = circular; unpress to make linear
        self.circular_action.setToolTip("Circular sequence (unpress to make it linear at a chosen base)")
        toolbar.addAction(self.circular_action)
        toolbar.addSeparator()
        toolbar.addWidget(self.search_bar)
        self.doc_actions.append(self.search_bar)
        # One box selects a stretch: "101..200" (1-based, both included; on a circular sequence "3700..50" runs through the
        # origin), or one number for the cursor after that base. It takes the room two number boxes and a button took.
        self.range_box = QLineEdit(placeholderText="Select 101..200", clearButtonEnabled=True, maxLength=40)
        self.range_box.setFixedWidth(self.range_box.fontMetrics().horizontalAdvance("100,000..100,100") + 30)   # (longer: scrolls)
        self.range_box.setToolTip("Select a stretch of bases: 101..200 (also 101-200 or 101 to 200), both ends included.\n"
                                  "On a circular sequence, 3700..50 selects through the origin. One number puts the cursor after "
                                  "that base. Enter selects.")
        self.range_box.returnPressed.connect(self.select_range)
        pipe = toolbar.addSeparator  # tools apart
        pipe()
        # The genetic code in use, a click away (hidden by Tools > Settings > General > Genetic code).
        self.code_switch = CodeSwitch(self)
        self.code_switch.manage.connect(lambda: codon_dialog.CodonDialog(self).exec())
        self.code_switch_actions = [toolbar.addWidget(self.code_switch), pipe()]
        self._show_code_switch()
        tool_settings.notifier.changed.connect(lambda tool: tool == genetic_code.TOOL and self._show_code_switch())
        toolbar.addWidget(self.range_box)
        self.doc_actions.append(self.range_box)
        pipe()
        # The background genome the primers (and spacers) are checked against: the last ones used, and "Add genome file…".
        self.reference_picker = ReferencePicker(self)
        self.reference_picker.changed.connect(lambda: self.editor is not None and self.primer_options_changed())
        toolbar.addWidget(self.reference_picker)
        pipe()
        # The side panel folds away to the window's edge; this button (or Ctrl+Shift+F) opens it.
        # Flat one-colour icons (icons.py): a funnel for the display filter, an hourglass for the history.
        filters = self.display_dock.toggleViewAction()
        filters.setIconText("Filters")
        filters.setToolTip("Display filter: what the views show (Ctrl+Shift+F)")
        toolbar.addAction(filters)   # (the two side panels side by side, one group)
        history = self.history_dock.toggleViewAction()
        history.setText("&History (Side Panel)")
        history.setIconText("History")
        history.setToolTip("History of the changes (Ctrl+H)")
        history.setShortcut("Ctrl+H")
        toolbar.addAction(history)
        view.addAction(history)
        # (Light / Dark / System is View > Theme: chosen once, it does not need room in the toolbar)
        self.panel_actions = (filters, history)
        self._paint_panel_icons()
        self.picker_separator = pipe()   # shown with the picker only (no line left dangling at the end of the toolbar)
        # Multi-sequence FASTA: which sequence is open.
        self.picker = QComboBox(minimumContentsLength=16)
        self.picker.setToolTip("Sequence (chromosome) of this FASTA file")
        self.picker.activated.connect(lambda i: self.switch_sequence(self.picker.itemText(i)))
        self.picker_action = toolbar.addWidget(self.picker)
        self.picker_action.visibleChanged.connect(lambda: self.picker_separator.setVisible(self.picker_action.isVisible()))
        self.picker_action.setVisible(False)
        self.picker_separator.setVisible(False)
        for a in self.doc_actions:
            a.setEnabled(False)

    # ---- files ---------------------------------------------------------------------
    def new_file(self):
        """File → New: a name, optionally pasted bases, circular or not. Opens here if this
        window is empty, else in a new window. Save asks where to write it."""
        dialog = QDialog(self)
        dialog.setWindowTitle("New Sequence")
        name = QLineEdit("New sequence", maxLength=100)
        bases = QPlainTextEdit(placeholderText="Paste bases here (optional; numbers, spaces and FASTA headers are removed)")
        circular = QCheckBox("Circular (plasmid)", checked=True)
        problem = QLabel(styleSheet=f"color: {theme.hex('status_bad')}", wordWrap=True)   # what is wrong with the bases, shown without closing the dialog
        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)

        def build():
            text = clean_pasted(bases.toPlainText()) if bases.toPlainText().strip() else ""
            return Document.new_sequence(name.text().strip() or "New sequence", text, circular.isChecked())

        def accept():
            """OK: check the bases here, so a mistake keeps the dialog (and what was pasted) open."""
            try:
                build()
            except (OSError, ValueError) as e:
                problem.setText(str(e))
                bases.setFocus()
                return
            dialog.accept()

        buttons.accepted.connect(accept)
        buttons.rejected.connect(dialog.reject)
        form = QFormLayout(dialog)
        form.addRow("Name", name)
        form.addRow("Sequence", bases)
        form.addRow("", problem)
        form.addRow("", circular)
        form.addRow(buttons)
        if dialog.exec() != QDialog.Accepted:
            return None
        try:
            doc = build()
        except (OSError, ValueError) as e:   # not reached through OK (which checks first), but never a traceback
            QMessageBox.warning(self, "New Sequence", str(e))
            return None
        return self.open_untitled(doc, Path(runner.safe_name(doc.name) + ".gb"))

    def open_untitled(self, doc, path):
        """Show a document made elsewhere (a new sequence, an NCBI record): in this window if it is empty, else in a new one, so
        nothing open is touched. Save asks where to write it."""
        if self.editor is None:
            self.show_document(doc, path, untitled=True)
            return self
        return open_window(doc, path, untitled=True)

    def import_ncbi(self, database=None):
        dialog = NcbiDialog(self)
        if database:
            dialog.database.setCurrentText(next(k for k in ncbi.DATABASES if k.startswith(database)))
        dialog.exec()

    def import_genome(self):
        GenomeDialog(self).exec()

    def import_as_chromosome(self, path=None):
        """Import > DNA > Sequence File as New Chromosome: the sequences of another file become new records of the open file. A
        GenBank file with one record becomes a multi-record file (the map then lists the records, one is open at a time); a FASTA
        file gets more sequences. The file is written now (see add_record), then the first added record is shown. True if
        something was added."""
        ed = self.editor
        if ed is None:
            return False
        title = "Import as New Chromosome"
        # A .dna file is never changed, and an unsaved sequence or a project has no GenBank file to add to: save one first.
        if ed.untitled or ed.path.suffix.lower() in (".dna", exports.PROJECT_EXT):
            why = ("is not saved yet" if ed.untitled else
                   "is a .dna file, which openDNA never changes" if ed.path.suffix.lower() == ".dna" else "is a project")
            if QMessageBox.question(self, title, f"{ed.path.name} {why}. Save it as a GenBank file first? The new chromosome is "
                                    "then added to that file.", QMessageBox.Save | QMessageBox.Cancel) != QMessageBox.Save:
                return False
            target, _ = QFileDialog.getSaveFileName(self, "Save As GenBank", str(ed.path.with_suffix(".gb")), "GenBank (*.gb *.gbk)")
            if not target or Path(target).suffix.lower() in (".dna", exports.PROJECT_EXT) or not self._write(Path(target)):
                return False
        fasta = bool(ed.fasta_record) and not ed.multi_genbank
        if path is None:
            path, _ = QFileDialog.getOpenFileName(self, "Import Sequence File as New Chromosome", "", OPEN_FILTER)
            if not path:
                return False
        try:
            records = run_busy(self, f"Reading {Path(path).name}…", lambda: add_record.source_records(path))
        except (OSError, ValueError) as e:
            QMessageBox.critical(self, "Cannot Open File", f"{path}\n\n{e}")
            return False
        dialog = AddRecordsDialog(self, ed.path, fasta, records, add_record.current_names(ed))
        if dialog.exec() != QDialog.Accepted:
            return False
        chosen = dialog.chosen()
        current = ed.fasta_record or add_record.current_names(ed)[0]   # a one-record GenBank file: the record open now

        try:
            docs = run_busy(self, f"Reading {Path(path).name}…", lambda: [(r.load(), name) for r, name in chosen])
            # primers of a .dna file (its primer list) are written as primer_bind features, as Save As GenBank does; read here,
            # in the window's thread (an Editor is a Qt object)
            primers = None if fasta else [Editor(doc, Path(path)).file_primer_features() for doc, _ in docs]

            def work():
                if fasta:
                    multigb.append_fasta(ed.path, [(name, doc) for doc, name in docs])
                else:
                    multigb.append_genbank(ed.path, [multigb.renamed_genbank(doc, name, p) for (doc, name), p in zip(docs, primers)])
            run_busy(self, f"Adding to {ed.path.name}…", work)
        except (OSError, ValueError) as e:
            QMessageBox.critical(self, title, f"{ed.path}\n\n{e}")
            return False
        names = [name for _, name in chosen]
        # The open record stays in memory with its edits (load keeps it, like switching records); the file is now read as one of
        # several records.
        ed.fasta_record, ed.multi_genbank = current, not fasta
        multigb.forget_siblings(ed)
        if self.load(str(ed.path), names[0], reread=True):
            self.statusBar().showMessage(
                f"Added {', '.join(names)} to {self.editor.path.name}: {plural(len(self.editor.fasta_records), 'sequence')} "
                "in the file.", 8000)
        return True

    # Files dropped on the window open like File > Open: here if the window is empty, else each in a window of its own.
    def dragEnterEvent(self, event):
        if any(u.isLocalFile() for u in event.mimeData().urls()):
            event.acceptProposedAction()

    def dropEvent(self, event):
        paths = [u.toLocalFile() for u in event.mimeData().urls() if u.isLocalFile() and Path(u.toLocalFile()).is_file()]
        if paths:
            event.acceptProposedAction()
            for path in paths:
                self.open_path(path)

    def open_path(self, path):
        """Open a file in this window if it is empty, else in a new one. The window, or None if it could not be read."""
        if self.editor is None:
            return self if self.load(path) else None
        win = MainWindow()
        WINDOWS.append(win)
        win.show()
        if win.load(path):
            return win
        win.close()
        WINDOWS.remove(win)
        return None

    def import_protein(self):
        """Import > Proteins > Protein FASTA: its first protein in the 3D structure window."""
        path, _ = QFileDialog.getOpenFileName(self, "Import Protein FASTA", "", "Protein FASTA (*.fa *.faa *.fasta *.txt);;All files (*)")
        if not path:
            return None
        try:
            name, protein = ncbi.protein_of(Path(path))
        except (OSError, UnicodeDecodeError, IndexError) as e:
            QMessageBox.warning(self, "Protein FASTA", f"{path}\n\n{e or 'empty file'}")
            return None
        self.structure_window = StructureDialog(self, name, protein, "")
        self.structure_window.show()
        return self.structure_window

    def import_primer_list(self):
        path, _ = QFileDialog.getOpenFileName(self, "Import Primers", "", "Primer lists (*.tsv *.csv *.txt *.fa *.fasta);;All files (*)")
        if path:
            try:
                rejected = []
                report = import_primers(path, self.editor, rejected)
                (QMessageBox.warning if rejected else QMessageBox.information)(self, "Import primers", report)
            except (OSError, UnicodeDecodeError, ValueError) as e:
                QMessageBox.warning(self, "Import Primers", str(e))

    def _fill_recent_menu(self):
        """File > Open Recent: the files of the start page's list (pinned ones first), each opened here or in a new window."""
        menu = self.recent_menu
        menu.clear()
        pinned = pinned_files()
        paths = pinned + [p for p in recent_files() if p not in pinned]
        for i, path in enumerate(paths):
            a = menu.addAction(f"{'&' + str(i + 1) + ' ' if i < 9 else ''}{Path(path).name.replace('&', '&&')}",
                               lambda _=False, path=path: self.open_path(path))
            a.setToolTip(path)
            a.setStatusTip(path)
            if path in pinned:
                f = a.font()
                f.setBold(True)
                a.setFont(f)
        if not paths:
            menu.addAction("No recent files").setEnabled(False)
        menu.setToolTipsVisible(True)

    def choose_file(self):
        path, _ = QFileDialog.getOpenFileName(self, "Open Sequence", "", OPEN_FILTER)
        if path:
            self.load(path)

    def load(self, path, record=None, reread=False):
        """Open a GenBank/.dna file, or one sequence (`record`, default the first) of
        a FASTA file. Big files load in the background behind a busy dialog. Another record of the file that is open keeps the
        list of records the window has (the file is not read through again for it: seconds on a genome of gigabytes), unless
        `reread` (records were just added to the file)."""
        old = self.editor
        # Another record of the file that is open: the tabs stay as they are, the sequence in them is switched, and the
        # record that was shown stays in memory (with its unsaved edits) until the file is closed or another one opened.
        same_file = bool(record) and old is not None and bool(old.fasta_record) and Path(old.path) == Path(path)
        if same_file and record == old.fasta_record:
            return True
        if not same_file:
            if not self._confirm_discard():
                return False
            self.record_docs.clear()
            self.file_history.reset()
        self._drop_autosave()
        cached = self.record_docs.get(record) if same_file else None
        names, lengths, meta, gb_record = [], [], None, False
        project = Path(path).suffix.lower() == exports.PROJECT_EXT
        try:
            if same_file and not reread and record in dict(old.fasta_records):
                lengths = list(old.fasta_records)
                names, gb_record = [n for n, _ in lengths], old.multi_genbank
            elif is_fasta(path):
                fasta = run_busy(self, f"Indexing {Path(path).name}…", lambda: Fasta(str(path)))
                lengths = [(n, int(l)) for n, l in fasta.sequences()]   # in file order (the map's list can reorder the file)
                names = [n for n, _ in lengths]
                if not names:
                    raise OSError("no sequence in this FASTA file")
                record = record if record in names else names[0]
            elif not project and multigb.is_multi_genbank(path):   # treated like a multi-FASTA: one record open, the map lists them all
                lengths = [(n, int(l)) for n, l, _ in multigb.records(path)]
                names, gb_record = [n for n, _ in lengths], True
                record = record if record in names else names[0]
            else:
                record = None
            recovery = autosave_path(path, record)
            recovered = False
            if cached is None and recovery.exists() and recovery.stat().st_mtime > Path(path).stat().st_mtime and \
                    self._ask_recover(path, recovery):
                doc = Document(str(recovery))
                doc.modified = True  # not saved to the real file yet
                recovered = True
                if project:
                    meta = exports.open_project(path)[1]
            elif cached is not None:
                doc = cached
            else:
                recovery.unlink(missing_ok=True)  # declined or stale
                label = f"Loading {record or Path(path).name}…"
                if project:
                    doc, meta = run_busy(self, label, lambda: exports.open_project(path))
                else:
                    doc = run_busy(self, label, (lambda: multigb.read_record(path, record)) if gb_record
                                   else (lambda: Document.from_fasta(str(path), record)) if record
                                   else (lambda: Document(str(path))))
        except OSError as e:
            QMessageBox.critical(self, "Cannot Open File", f"{path}\n\n{e}")
            return False
        if recovered:
            self.autosave_file = recovery
        add_recent(path)
        # Another sequence of the FASTA file that is open: its tracks (GFF, BAM...) follow, each on the sequence of the same name.
        carry = ([e.path for e in old.tracks + old.alignments if getattr(e, "path", None)]
                 if old is not None and record and old.fasta_record and Path(old.path) == Path(path) else [])
        keep_tabs = same_file and meta is None
        letters = self.sequence_view.size if same_file else None
        if same_file:
            self.record_docs[old.fasta_record] = old.doc
            self.record_docs.pop(record, None)
        # Another record of the same file: nothing is painted until the new record is ready, so the old one stays on screen
        # and then the new one replaces it in one go (the tabs are built again: painted as they were made, a window showed
        # half-built or black panes for a moment on a big record).
        frozen = same_file and self.updatesEnabled()
        if frozen:
            self.setUpdatesEnabled(False)
        try:
            self.show_document(doc, path, record=record, sequences=names, lengths=lengths, genbank=gb_record,
                               layout=self.workspace.dump() if keep_tabs else (meta or {}).get("layout"),
                               current=self.workspace.currentIndex() if keep_tabs else 0,   # a project opens with the tabs it was saved with
                               display=not same_file)   # (the record switched to takes the file's display settings: below)
            if same_file:
                self._carry_settings(old, letters)
            if frozen:
                self._settle_layout()
        finally:
            if frozen:
                self.setUpdatesEnabled(True)
        if carry:
            skipped = [Path(p).name for p in carry if self.load_track(p, quiet=True) is None]
            if skipped:
                self.statusBar().showMessage(f"No data for {record} in: {', '.join(skipped)}", 8000)
        if meta is not None:
            self._apply_project(meta, path)
        return True

    def _settle_layout(self):
        """Show and lay the new tabs out now, while nothing is painted: Qt shows a new central widget from a posted call, then come
        the layouts, the sizes they give (the sequence view lays its bases out for its width), the panes built when first shown,
        all of which otherwise ran after painting was back on. On Windows the window was then painted at once with the new tabs
        not shown or sized yet, showing only its dark background for the second that work took. (Fitting table columns is left
        for after: it only adjusts widths, and took half a second on a genome.)"""
        for _ in range(3):   # the new tabs are shown by a posted call (Qt's own), then laid out; a size change can ask for one more pass
            QApplication.sendPostedEvents()

    def _carry_settings(self, old, letters):
        """Switching to another record of the file: the display settings, the lock and the letter size are the file's, not the
        record's, so they stay as they were."""
        ed = self.editor
        ed.untranslated, ed.hidden_features = set(old.untranslated), set(old.hidden_features)
        ed.set_display(enzymes=old.show_enzymes, enzyme_filter=old.enzyme_filter, translation=old.show_translation,
                       collapse_genes=old.collapse_genes, hidden_kinds=old.hidden_kinds, primers=old.show_primers,
                       cds_translation=old.show_cds_translation, file_primers=old.show_file_primers,
                       db_primers=old.show_db_primers, units=old.show_units, cds_only=old.cds_only, complement=old.show_complement)
        temporary =[(n, s) for n, s in old.extra if (n, s.upper()) in old.temporary_primers]
        ed.set_extra([p for p in old.extra if p not in temporary], temporary)
        self.set_dna_locked(old.dna_locked)
        if letters:
            self.sequence_view.set_size(letters)

    def _apply_project(self, meta, path):
        """What an Open DNA Project holds besides the GenBank part: display, tracks and BAM files, the DNA lock."""
        ed = self.editor
        d = meta.get("display", {})
        ed.untranslated = {tuple(x) for x in d.get("untranslated", [])}
        ed.hidden_features = {tuple(x) for x in d.get("hidden_features", [])}
        ed.set_display(enzymes=d.get("enzymes"), enzyme_filter=EnzymeFilter.from_dict(d.get("enzyme_filter")), translation=d.get("translation"), collapse_genes=d.get("collapse_genes"), cds_only=d.get("cds_only"),
                       hidden_kinds=d.get("hidden_kinds"), primers=d.get("primers"), cds_translation=d.get("cds_translation"),
                       file_primers=d.get("file_primers"), db_primers=d.get("db_primers"), units=d.get("units"),
                       complement=d.get("complement"))
        missing = []
        for entry in meta.get("tracks", []) + meta.get("alignments", []):
            found = exports.project_file(entry, path)
            if found is None:
                missing.append(entry["path"])
            else:
                self.load_track(found, entry.get("seqid"))
        self.set_dna_locked(bool(meta.get("dna_locked")))
        said = [m for m in (genetic_code.adopt(meta.get("genetic_code")), codon_usage.adopt(meta.get("codon_usage"))) if m]
        if said:
            self.statusBar().showMessage(" ".join(said), 8000)
        for name in meta.get("hidden_primers", []):
            ed.set_primer_hidden(name, True)
        if missing:
            QMessageBox.warning(self, "Project", "These files of the project were not found (moved or deleted):\n\n" + "\n".join(missing))

    def _observe_history(self):
        ed = self.editor
        undo, redo = ed.doc.history()
        self.file_history.observe(ed.fasta_record, len(undo), len(redo))

    def _ask_recover(self, path, recovery):
        """Recover the autosaved changes? Buttons that say what happens: declining deletes the recovery copy for good."""
        box = QMessageBox(QMessageBox.Question, "Recover Unsaved Changes",
                          f"{Path(path).name} has unsaved changes from an earlier session (autosaved {recovery.name}).", parent=self)
        box.setInformativeText("Discarding deletes the recovery copy; the file itself is not changed either way.")
        recover = box.addButton("&Recover Changes", QMessageBox.AcceptRole)
        box.addButton("&Discard Changes", QMessageBox.DestructiveRole)
        box.setDefaultButton(recover)
        box.exec()
        return box.clickedButton() is recover

    def _file_histories(self):
        """{record: Document.history()} of every record held in memory, the open one included."""
        out = {name: doc.history() for name, doc in self.record_docs.items()}
        out[self.editor.fasta_record] = self.editor.doc.history()
        return out

    def undo(self):
        """Edit > Undo: the last edit of the file, in whichever record it was (that record is then shown)."""
        self._step_file(self.file_history.next_undo, lambda: self.editor.undo())

    def redo(self):
        self._step_file(self.file_history.next_redo, lambda: self.editor.redo())

    def _step_file(self, which, step):
        ed = self.editor
        record = which() if ed.fasta_record else None
        if record and record != ed.fasta_record and not self.load(ed.path, record):
            return
        step()

    def undo_file_steps(self, count):
        """History panel: `count` steps back through the whole file (forward when negative)."""
        for _ in range(abs(count)):
            (self.undo if count > 0 else self.redo)()

    def switch_sequence(self, name):
        """Picker: open another sequence of the same FASTA file."""
        ed = self.editor
        if name != ed.fasta_record and not self.load(ed.path, name):
            self.picker.setCurrentText(ed.fasta_record)  # cancelled: show what is open

    # ---- the records of a file of several (the map's chromosome list) ------------------------------------------------------------
    def _set_record_list(self, lengths):
        """The file's records are now `lengths` [(name, length)] (reordered or one removed): the map's list and the picker follow."""
        ed = self.editor
        ed.fasta_records = list(lengths)
        multigb.forget_siblings(ed)
        self.picker.clear()
        self.picker.addItems([n for n, _ in lengths])
        self.picker.setCurrentText(ed.fasta_record)
        self.picker_action.setVisible(len(lengths) > 1)
        ed.match_elsewhere()
        ed.display_changed.emit()

    def reorder_records(self, names):
        """Map's list, a line dragged (or Move Up / Down): the records are written into the file in this order, each unchanged.
        The edits not saved yet stay in memory (records are found by name)."""
        ed = self.editor
        lengths = dict(ed.fasta_records)
        if sorted(names) != sorted(lengths) or list(names) == list(lengths):
            return False
        try:
            run_busy(self, f"Reordering {ed.path.name}…", lambda: multigb.rewrite_records(ed.path, list(names)))
        except OSError as e:
            QMessageBox.critical(self, "Reorder Sequences", f"{ed.path}\n\n{e}")
            return False
        self._set_record_list([(n, lengths[n]) for n in names])
        self.statusBar().showMessage(f"Sequences of {ed.path.name} reordered: {', '.join(names[:6])}{'…' if len(names) > 6 else ''}", 5000)
        return True

    def remove_record(self, name):
        """Map's list, right-click > Remove from File: the record is deleted from the file (after a question; there is no undo).
        The open record is first switched to the next one."""
        ed = self.editor
        names = [n for n, _ in ed.fasta_records]
        if name not in names or len(names) < 2:
            return False
        edited = (ed.doc.modified if name == ed.fasta_record else
                  bool(self.record_docs.get(name) is not None and self.record_docs[name].modified))
        if not self._confirm_remove(name, edited):
            return False
        if name == ed.fasta_record:   # show another record first (the open one stays in memory, then is dropped below)
            i = names.index(name)
            if not self.load(str(ed.path), names[i + 1] if i + 1 < len(names) else names[i - 1]):
                return False
            ed = self.editor
        keep = [n for n in names if n != name]
        try:
            run_busy(self, f"Removing {name}…", lambda: multigb.rewrite_records(ed.path, keep))
        except OSError as e:
            QMessageBox.critical(self, "Remove from File", f"{ed.path}\n\n{e}")
            return False
        self.record_docs.pop(name, None)
        self.file_history.forget(name)
        autosave_path(ed.path, name).unlink(missing_ok=True)
        lengths = dict(ed.fasta_records)
        self._set_record_list([(n, lengths[n]) for n in keep])
        self._update_undo_actions()
        self.statusBar().showMessage(f"Removed {name} from {ed.path.name}: {plural(len(keep), 'sequence')} left.", 5000)
        return True

    def _confirm_remove(self, name, edited):
        """Remove the record? Cancel is the default: nothing written can be undone."""
        box = QMessageBox(QMessageBox.Warning, "Remove from File", f"Remove {name} from {self.editor.path.name}?", parent=self)
        box.setInformativeText("The file is written now and this cannot be undone"
                               + (f": the unsaved edits of {name} are lost too." if edited else "."))
        remove = box.addButton("&Remove", QMessageBox.DestructiveRole)
        box.addButton(QMessageBox.Cancel)
        box.setDefaultButton(QMessageBox.Cancel)
        box.exec()
        return box.clickedButton() is remove

    def export_record(self, name, target=None):
        """Map's list, right-click > Export as GenBank: one record of the file as a GenBank file of its own, with its unsaved
        edits (the file itself is not changed)."""
        ed = self.editor
        if target is None:
            suggested = ed.path.with_name(runner.safe_name(name) + ".gb")
            target, _ = QFileDialog.getSaveFileName(self, f"Export {name} as GenBank", str(suggested), "GenBank (*.gb *.gbk)")
            if not target:
                return False
        target = Path(target)
        if target.exists() and target.resolve() == ed.path.resolve():
            QMessageBox.warning(self, "Export as GenBank", "Choose another file: this one holds every record of the file.")
            return False
        try:
            if name == ed.fasta_record:
                doc, primers = ed.doc, ed.file_primer_features()
            elif self.record_docs.get(name) is not None:
                doc, primers = self.record_docs[name], None
            else:
                doc = run_busy(self, f"Reading {name}…", (lambda: multigb.read_record(ed.path, name)) if ed.multi_genbank
                               else (lambda: Document.from_fasta(str(ed.path), name)))
                primers = None
            run_busy(self, f"Exporting {name}…", lambda: doc.write_genbank(str(target), primers))
        except OSError as e:
            QMessageBox.critical(self, "Export as GenBank", f"{target}\n\n{e}")
            return False
        self.statusBar().showMessage(f"Exported {name} to {target}", 5000)
        return True

    def show_document(self, doc, path, untitled=False, record=None, sequences=(), layout=None, lengths=(), genbank=False, current=0,
                      display=True):
        """Put a Document in this window. `untitled`: made in memory (e.g. a PCR product);
        `path` is only the suggested file name, Save asks where to save. `record`: the
        sequence of the FASTA file `path` that `doc` holds. `layout`: the tabs to open with
        (an Open DNA Project's own); everything else opens with the standard ones. `display` False: the caller sets the display
        options itself (another record of the file keeps the file's: laying a genome out twice took most of the switch)."""
        self.editor = ed = Editor(doc, path)
        ed.untitled = untitled
        ed.fasta_record = record
        ed.fasta_records = list(lengths)
        ed.multi_genbank = genbank
        ed.record_requested.connect(self.switch_sequence)
        ed.records_order_requested.connect(self.reorder_records)
        ed.record_remove_requested.connect(self.remove_record)
        ed.record_export_requested.connect(self.export_record)
        ed.record_hit_requested.connect(self.go_to_hit)
        self.picker.clear()
        self.picker.addItems(list(sequences))
        if record:
            self.picker.setCurrentText(record)
        self.picker_action.setVisible(len(sequences) > 1)
        self.autosaved_version = ed.version
        ed.match_elsewhere()   # a multi-record file: where the primers bind on the other records
        ed.changed.connect(self._update_title)
        if record:   # a file of several records has one history: Undo goes back through all of them
            ed.file_history, ed.file_histories, ed.file_undo_steps = self.file_history, self._file_histories, self.undo_file_steps
            ed.changed.connect(lambda: ed is self.editor and self._observe_history())
            self._observe_history()
        ed.selection_changed.connect(self._update_status)
        ed.error.connect(lambda msg: self.statusBar().showMessage(msg, 5000))
        ed.edit_feature_requested.connect(lambda i: edit_feature(self, ed, i))
        ed.edit_primer_requested.connect(lambda name: PrimerDialog(self, ed, name).exec())
        if display:
            ed.set_display(enzymes=self.enzymes_action.isChecked(), translation=self.translation_action.isChecked(),
                           complement=self.complement_action.isChecked(),
                           collapse_genes=self.collapse_action.isChecked(), cds_only=self.cds_only_action.isChecked())
        ed.display_changed.connect(self._sync_display_actions)
        self._show_lock(ed.dna_locked)
        self.search_bar.attach(ed)
        self.display_dock.setWidget(DisplayPanel(ed))
        self.history_dock.setWidget(HistoryPanel(ed))
        self.range_box.clear()

        # Always the standard tabs (see workspace.DEFAULT); only a project carries its own layout.
        self._build_workspace(layout or DEFAULT_LAYOUT, current)
        for a in self.doc_actions:
            a.setEnabled(True)
        ed.changed.connect(lambda: ed is self.editor and self._update_undo_actions())
        self._update_undo_actions()
        self.sequence_view.setFocus()
        self._update_title()
        self._update_status()
        self._collect_primers()

    # The first view of each kind (current tab first), for menus and tests. Built on
    # demand if its pane was never shown.
    sequence_view = property(lambda self: self.workspace.ensure(SequenceView))
    map_view = property(lambda self: self.workspace.ensure(MapView))
    genome_view = property(lambda self: self.workspace.ensure(GenomeView))
    feature_table = property(lambda self: self.workspace.ensure(FeatureTable))
    primer_panel = property(lambda self: self.workspace.ensure(PrimerPanel))

    def set_genome_visible(self, on):
        """View → Alignment View: show/hide every genome pane (adding one if there is none). Switching it off never
        builds a view: an alignment view that was never shown stays unbuilt."""
        self.genome_action.setChecked(on)
        views = self.workspace.views(GenomeView, build="all" if on else False)
        if on and not views:
            views = [self.workspace.ensure(GenomeView)]
        for g in views:
            g.setVisible(on)
            g.parentWidget().setVisible(on)  # its pane, header included
        if views:
            self.workspace.refresh_names()   # the tab is named after the panes you can see

    def open_product(self, doc):
        """Open an assembled sequence in a new unsaved window."""
        folder = self.editor.path.parent if self.editor else Path.cwd()
        return open_window(doc, folder / f"{doc.name}.gb", untitled=True)

    def about(self):
        """Help > About: the notices the AGPL asks an interactive program to show (copyright, no warranty, the licence)."""
        QMessageBox.about(
            self, f"About {project.NAME}",
            f"<b>{project.NAME}</b> {project.VERSION}<br>A local, offline DNA editor.<br><br>"
            f"Developed by {project.AUTHOR}. Copyright © 2026 {project.AUTHOR}.<br>"
            f"{project.NAME} is free software. The official download is "
            f"<a href='{project.REPO_URL}'>github.com/{project.REPO}</a>. Forks are welcome under a name of their own "
            "(TRADEMARKS.md).<br><br>"
            "This program is free software: you can redistribute it and/or modify it under the terms of the "
            "<b>GNU Affero General Public License</b> version 3 (see the LICENSE file).<br><br>"
            "It comes with <b>NO WARRANTY</b>, to the extent permitted by law. Always keep a copy of your original files.<br><br>"
            f"Built with Qt {qVersion()} (PySide6 {pyside_version}), under the LGPL-3.0; see THIRD_PARTY_NOTICES.md. "
            "The folder version of the download lets you replace Qt.")

    def update_hooks(self):
        """What the updater needs from the window: how to quit, and a chance to refuse while there are unsaved changes."""
        return {"quit": self.close_for_update, "confirm": self.confirm_all_saved}

    def confirm_all_saved(self):
        """False if the user cancels while any open window has unsaved changes (each one asks)."""
        return all(w._confirm_discard() for w in [self] + list(WINDOWS) if w.editor is not None and w.editor.doc.modified)

    def close_for_update(self):
        """The new version is in place and started: close this window and the others without asking again."""
        for w in [self] + list(WINDOWS):
            if w.editor is not None:
                w.editor.doc.modified = False
        QApplication.quit()

    def open_pcr_product(self, dialog):
        try:
            doc = dialog.product()
        except (OSError, ValueError) as e:
            QMessageBox.warning(self, "PCR product", str(e))
            return None
        return open_window(doc, self.editor.path.with_name(f"{doc.name}.gb"), untitled=True)

    def open_documents(self):
        """[(file name, Document)] of every file open in a window of this program (for the Assembly panel)."""
        return [(w.editor.path.name, w.editor.doc) for w in [self] + [w for w in WINDOWS if w is not self]
                if w.editor is not None]

    def show_assembly(self):
        """Tools > Assembly: the Assembly panel as a tab."""
        self.assembly = self.show_tool_tab("Assembly", lambda: AssemblyPanel(self))

    def show_tool_tab(self, title, make):
        """A tool as a tab (of the workspace, or on its own while no file is open): made by `make()` the first time, the same
        widget (as it was left) after. Returns the widget."""
        widget = self.tool_tabs.get(title)
        if widget is None:
            widget = self.tool_tabs[title] = make()
        host = self.workspace
        if host is None:
            if self.tool_host is None:
                self.tool_host = Workspace(self, tool_only=True)
                self.tool_host.tool_closed.connect(self._tool_host_emptied)
                self.setCentralWidget(self.tool_host)
            host = self.tool_host
        host.show_tool(widget, title)
        return widget

    def _tool_host_emptied(self):
        if self.tool_host is not None and self.tool_host.count() == 0:
            self.tool_host = None
            self.setCentralWidget(QLabel("Open a GenBank or .dna file (Ctrl+O)", alignment=Qt.AlignCenter))

    def _show_code_switch(self):
        for action in self.code_switch_actions:   # the button and the line after it
            action.setVisible(code_switch_shown())

    def primer_options(self):
        self.open_settings("Primers")

    def open_settings(self, page=None):
        """Tools > Settings, on `page` (a view or tool name) when given; the primer lists redo their matching after OK."""
        if SettingsDialog(self, page).exec() == QDialog.Accepted and self.editor is not None:
            self.primer_options_changed()

    def primer_options_changed(self):
        """The primer options (and the background genome they hold, shared with the Cas9 tool) changed."""
        for panel in self.workspace.views(PrimerPanel):
            panel.reload()
        for tool in self.workspace.views(Cas9Dialog):
            tool.update_genome_label()
        self.editor.selection_changed.emit()  # status bar Tm

    def _build_workspace(self, layout=None, current=0):
        # The tool tabs that are open (Assembly, the databases) move to the new workspace, as they were.
        tools, shown = [], None
        for host in (self.workspace, self.tool_host):
            for title, widget in self.tool_tabs.items():
                if host is not None and host.indexOf(widget) >= 0:
                    tools.append((title, widget))
                    if host.currentWidget() is widget:
                        shown = widget
                    host.take_tool(widget)
        self.workspace = Workspace(self, layout)
        self.workspace.setCurrentIndex(min(current, self.workspace.count() - 1))
        self.setCentralWidget(self.workspace)  # the old one (if any) is deleted by Qt
        self.tool_host = None
        layouts = self.workspace.count()
        for title, widget in tools:
            self.workspace.addTab(widget, title)
        if shown is not None:
            self.workspace.setCurrentWidget(shown)
        else:
            self.workspace.setCurrentIndex(min(current, layouts - 1))
        self.workspace.ensure(SequenceView)   # a document always has one, whatever the layout says
        self.genome_action.setChecked(bool(self.workspace.views(GenomeView)))

    def reset_layout(self):
        """View → Reset Tabs: the default tabs, same document, undo history kept."""
        self._build_workspace(DEFAULT_LAYOUT)

    def _collect_primers(self):
        """Every primer in an opened file goes into the primer database (by name)."""
        ed = self.editor
        found = dict(ed.doc.primers())  # the .dna file's primer list, including unbound ones
        for f in ed.file_primers:
            found.setdefault(f.label, primer_sequence(f) or annealed(ed, f))
        report = db().add_many(found.items())
        if report["added"] or report["conflict"] or report["invalid"]:
            msg = f"{plural(report['added'], 'primer')} from this file added to the primer database."
            if report["conflict"]:
                msg += f" Not added, name used for another sequence: {', '.join(report['conflict'])}"
            if report["invalid"]:
                msg += " Not added, not a valid primer: " + "; ".join(f"{n} ({why})" for n, why in report["invalid"][:5])
            self.statusBar().showMessage(msg, 8000)

    def save(self):
        ed = self.editor
        # .dna files are read-only: saving one always goes through "Save As GenBank".
        if ed.path.suffix.lower() == ".dna" or ed.untitled:
            return self.save_as()
        if ed.fasta_record and ed.multi_genbank:
            return self._write_record()
        if ed.fasta_record:
            if ed.features:
                QMessageBox.information(self, "Save as GenBank",
                                        "FASTA files cannot store features: save as GenBank to keep them."
                                        + other_records_note(ed, self.picker.count()))
                return self.save_as()
            return self._write_fasta()
        return self._write(ed.path)

    def _write_record(self):
        """Save the open record back into its multi-record GenBank file, the other records unchanged."""
        ed = self.editor
        refreshed = self._refresh_before_save()
        try:
            run_busy(self, f"Saving {ed.path.name}…", lambda: multigb.save_record_of(ed))
        except OSError as e:
            QMessageBox.critical(self, "Save Failed", f"{ed.path}\n\n{e}")
            return False
        multigb.forget_siblings(ed)
        for name, doc in self.record_docs.items():   # records edited, then switched away from: saved with it
            if doc.modified:
                multigb.save_document(ed.path, name, doc)
        ed.changed.emit()
        self._drop_autosave()
        self.statusBar().showMessage(f"Saved {ed.fasta_record} in {ed.path}" + (f" · {refreshed}" if refreshed else ""), 3000)
        return True

    def _write_fasta(self):
        """Write the edited sequence back into its FASTA file (other sequences streamed
        through unchanged)."""
        ed = self.editor
        try:
            def work():
                ed.doc.write_fasta(str(ed.path), str(ed.path))
                for doc in self.record_docs.values():   # records edited, then switched away from: saved with it
                    if doc.modified:
                        doc.write_fasta(str(ed.path), str(ed.path))
                        doc.modified = False
            run_busy(self, f"Saving {ed.path.name}…", work)
        except OSError as e:
            QMessageBox.critical(self, "Save Failed", f"{ed.path}\n\n{e}")
            return False
        ed.doc.modified = False
        ed.changed.emit()  # title: no longer modified
        self._drop_autosave()
        self.statusBar().showMessage(f"Saved {ed.fasta_record} in {ed.path}", 3000)
        return True

    def _save_file_as(self):
        """Save As on a file of several records: the whole file is copied, with the edits of every record in it."""
        ed = self.editor
        path, _ = QFileDialog.getSaveFileName(self, "Save As", str(ed.path), "GenBank (*.gb *.gbk)")
        if not path:
            return False
        new = Path(path)
        if new.resolve() == ed.path.resolve():
            return self.save()
        refreshed = self._refresh_before_save()

        def work():
            shutil.copyfile(ed.path, new)
            multigb.save_document(new, ed.fasta_record, ed.doc, ed.file_primer_features())
            for name, doc in self.record_docs.items():
                if doc.modified:
                    multigb.save_document(new, name, doc)
        try:
            run_busy(self, f"Saving {new.name}…", work)
        except OSError as e:
            QMessageBox.critical(self, "Save Failed", f"{new}\n\n{e}")
            return False
        ed.path = new
        if ed.hidden_primers:
            remember_hidden_primers(new, ed.hidden_primers)
        multigb.forget_siblings(ed)
        add_recent(new)
        self._drop_autosave()
        ed.changed.emit()
        self.statusBar().showMessage(f"Saved {new}" + (f" · {refreshed}" if refreshed else ""), 3000)
        return True

    def save_as(self):
        ed = self.editor
        if ed.multi_genbank:
            return self._save_file_as()
        project = ed.path.suffix.lower() == exports.PROJECT_EXT
        suggested = str(ed.path.with_suffix(exports.PROJECT_EXT if project else ".gb"))
        filters = ["GenBank (*.gb *.gbk)", "Open DNA Project, GenBank + tracks and display (*.odp)"]
        path, _ = QFileDialog.getSaveFileName(self, "Save As", suggested, ";;".join(filters[::-1] if project else filters))
        return self._write(Path(path)) if path else False

    def _write(self, path):
        ed = self.editor
        refreshed = self._refresh_before_save()
        try:
            if Path(path).suffix.lower() == exports.PROJECT_EXT:
                run_busy(self, f"Saving {Path(path).name}…", lambda: exports.save_project(ed, path, self.workspace.dump()))
                ed.doc.modified = False
                ed.untitled, ed.path = False, Path(path)
                ed.changed.emit()
            else:
                ed.save(path)
        except OSError as e:
            QMessageBox.critical(self, "Save Failed", f"{path}\n\n{e}")
            return False
        self._drop_autosave()
        note = other_records_note(self.editor, self.picker.count())
        if self.editor.fasta_record:  # now a GenBank document
            self.record_docs.clear()
            self.editor.fasta_record, self.editor.multi_genbank, self.editor.fasta_records = None, False, []
            self.picker_action.setVisible(False)
            self.editor.display_changed.emit()
        self.statusBar().showMessage(f"Saved {path}{note}" + (f" · {refreshed}" if refreshed else ""), 8000 if note or refreshed else 3000)
        return True

    def _refresh_before_save(self):
        """Feature > Keep /translation up to Date: before a GenBank file is written, recompute the /translation of its CDSs.
        Only an explicit save comes here: autosave (`autosave`) writes the document as it is and never edits it. The refresh is
        one undo step, so Ctrl+Z after saving undoes it (and the file then counts as modified); when every /translation is
        already right nothing is pushed, so saving twice changes nothing and there is no loop. Returns the status text, or ''
        when nothing was updated or skipped."""
        if not genetic_code.keep_translation():
            return ""
        report = self.editor.refresh_translations()
        if all(status == "unchanged" for _, status, _ in report):
            return ""
        return genetic_code.summarize(report)

    def update_translations(self):
        """Feature > Update /translation Now: the selected CDS, or every CDS when none is selected (one undo step)."""
        ed = self.editor
        f = ed.by_index.get(ed.current_feature) if ed.current_feature is not None else None
        report = ed.refresh_translations([f.index] if f is not None and f.kind == "CDS" else None)
        if not report:
            self.statusBar().showMessage("No CDS with a /translation qualifier to update.", 5000)
        else:
            self.statusBar().showMessage(genetic_code.summarize(report) + (" (Ctrl+Z undoes it)." if any(s == "updated" for _, s, _ in report) else "."), 8000)

    def open_checkup(self):
        """Tools > Plasmid Checkup: a tab of panes, the checkup beside the document's own map and sequence (which then
        draw what it found). The tab already holding one is shown instead of a second."""
        ws = self.workspace
        for i in range(ws.count()):
            if ws.is_layout(i) and any(p.view_name == CHECKUP for p in ws.panes(ws.widget(i))):
                ws.setCurrentIndex(i)
                return
        ws.add_tab(None, split("h", [560, 900], pane(CHECKUP), split("v", [420, 480], pane("Map"), pane("Sequence"))),
                   select=True)

    def _paint_panel_icons(self):
        filters, history = self.panel_actions
        filters.setIcon(icons.filter_icon(self.palette()))
        self.lock_action.setIcon(icons.lock_icon(self.palette()))
        history.setIcon(icons.history_icon(self.palette()))
        for key, a in self.theme_actions.items():
            a.setIcon(icons.theme_icon(key, self.palette()))

    def set_theme(self, key):
        """View > Theme or the sun button: every window of the program follows (theme.apply repaints all the views)."""
        theme.apply(key)
        for w in QApplication.topLevelWidgets():   # the other windows' menus and sun buttons show the new choice too
            if isinstance(w, MainWindow):
                w.theme_actions[key].setChecked(True)
                w._paint_panel_icons()

    def changeEvent(self, event):
        # A theme switch changes the palette: the one-colour icons follow (black on light, white on dark).
        if event.type() == QEvent.PaletteChange and hasattr(self, "panel_actions"):
            self._paint_panel_icons()
        super().changeEvent(event)

    # ---- export ----------------------------------------------------------------------------------------------------------
    def _update_export_menu(self):
        """Grey each kind of export while there is nothing of it; one SAM entry per loaded BAM."""
        ed = self.editor
        has = ed is not None
        selected = has and ed.end > ed.start
        for menu, on in ((self.export_whole, has), (self.export_dna, has), (self.export_features, has and bool(ed.features)),
                         (self.export_protein, has and (bool(exports.cds_features(ed)) or selected)),
                         (self.export_alignments, has and bool(ed.alignments)), (self.export_primers, bool(db().rows))):
            menu.menuAction().setEnabled(on)
        self.export_selection.setEnabled(selected)
        self.export_image.setEnabled(self.workspace is not None and self.workspace.currentWidget() is not None)
        self.export_cds.setEnabled(has and bool(exports.cds_features(ed)))
        self.export_translation.setEnabled(selected and ed.end - ed.start >= 3)
        self.export_alignments.clear()
        for a in ed.alignments if has else []:
            self.export_alignments.addAction(QAction(f"{a.name} as SAM…", self, triggered=lambda _=False, a=a: self.export_sam(a)))

    def _export_path(self, title, suffix, filters, name=None):
        """Ask where to export: next to the open file, named after the sequence (or `name`)."""
        ed = self.editor
        folder = ed.path.parent if ed is not None and not ed.untitled else Path.cwd()
        name = name or (exports.seq_id(ed) if ed is not None else "export")
        path, _ = QFileDialog.getSaveFileName(self, title, str(folder / f"{name}.{suffix}"), filters)
        return Path(path) if path else None

    def _export(self, path, work):
        try:
            run_busy(self, f"Writing {path.name}…", work)
        except OSError as e:
            QMessageBox.critical(self, "Export Failed", f"{path}\n\n{e}")
            return None
        self.statusBar().showMessage(f"Exported {path}", 5000)
        return path

    def export_text(self, title, suffix, filters, make):
        """Export > a text format made by `make(editor)`."""
        path = self._export_path(title, suffix, filters)
        return path and self._export(path, lambda: path.write_text(make(self.editor), encoding="utf-8"))

    def export_tab_image(self):
        """Export > Image of the Current Tab: every pane of the tab shown, laid out as on screen."""
        root = self.workspace.currentWidget() if self.workspace is not None else None
        if root is None:
            return None
        dialog = image_export.ImageExportDialog(self, root)
        if not dialog.exec():
            return None
        fmt, scale, dpi, headers = dialog.options()
        name = exports.seq_id(self.editor) if self.editor is not None else "export"
        tab = re.sub(r"[^A-Za-z0-9_.-]+", "_", self.workspace.tabText(self.workspace.currentIndex())).strip("_")
        path = self._export_path("Export image", fmt, f"{fmt.upper()} image (*.{fmt}{' *.tif' if fmt == 'tiff' else ''})",
                                 name=f"{name}_{tab}" if tab else name)
        if not path:
            return None
        views = None if headers else self.workspace.tab_views(root)
        try:  # painting must happen on the GUI thread: no run_busy here
            image_export.write(root, path, fmt, scale, dpi, views)
        except OSError as e:
            QMessageBox.critical(self, "Export Failed", f"{path}\n\n{e}")
            return None
        self.statusBar().showMessage(f"Exported {path}", 5000)
        return path

    def export_genbank(self):
        """A GenBank copy (the open file keeps its name; File > Save As renames it)."""
        ed = self.editor
        path = self._export_path("Export GenBank", "gb", "GenBank (*.gb *.gbk)")
        return path and self._export(path, lambda: ed.doc.write_genbank(str(path), ed.file_primer_features()))

    def export_project(self):
        path = self._export_path("Export Open DNA Project", "odp", "Open DNA Project (*.odp)")
        return path and self._export(path, lambda: exports.save_project(self.editor, path, self.workspace.dump()))

    def export_submission(self):
        path = self._export_path("Export for NCBI submission (table2asn)", "tbl", "NCBI feature table + FASTA (*.tbl)")
        return path and self._export(path, lambda: exports.write_submission(self.editor, path))

    def export_fasta(self, selection):
        ed = self.editor
        path = self._export_path("Export FASTA", "fasta", "FASTA (*.fasta *.fa *.fna)")
        if not path:
            return None
        name = exports.seq_id(ed)
        if selection:
            end = ed.end - ed.length if ed.wraps else ed.end
            name, seq = f"{name}:{ed.start + 1}-{end}", ed.selected_sequence()
        else:
            seq = ed.doc.sequence()
        return self._export(path, lambda: path.write_text(exports.fasta_text(name, seq), encoding="utf-8"))

    def export_plain(self):
        path = self._export_path("Export bases", "txt", "Text (*.txt *.seq)")
        return path and self._export(path, lambda: path.write_text(self.editor.doc.sequence() + "\n", encoding="utf-8"))

    def export_translated(self):
        ed = self.editor
        path = self._export_path("Export translation", "faa", "Protein FASTA (*.faa *.fasta)")
        name = f"{exports.seq_id(ed)}:{ed.start + 1}-{ed.end}"
        return path and self._export(path, lambda: path.write_text(
            exports.fasta_text(name, exports.translate(ed.selected_sequence())), encoding="utf-8"))

    def export_sam(self, a):
        path = self._export_path(f"Export {a.name} as SAM", "sam", "SAM (*.sam)", Path(a.name).stem)
        return path and self._export(path, lambda: a.bam.write_sam(str(path)))

    def export_primer_db(self):
        path, _ = QFileDialog.getSaveFileName(self, "Export Primer Database", "primers.tsv",
                                              "Tab-separated (*.tsv);;CSV (*.csv);;FASTA (*.fasta)")
        return path and self._export(Path(path), lambda: db().export(path))

    def export_primers_for_ordering(self):
        """The primers selected in the Primer Database tab when it is on view, else all of them."""
        tab = self.tool_tabs.get(db_tabs.TITLES[db_tabs.PRIMERS])
        export_for_ordering(self, tab.table.chosen() if tab is not None and tab.isVisible() and tab.table is not None else None)

    # ---- DNA lock --------------------------------------------------------------------------------------------------------
    def set_dna_locked(self, on):
        """Lock button / Edit > Lock DNA."""
        if self.editor is not None:
            self.editor.dna_locked = on
        self._show_lock(on)

    def _show_lock(self, on):
        self.lock_action.setChecked(on)   # the padlock icon shows the state (shut when pressed), the text stays "DNA"

    def pcr_product(self):
        dialog = PcrDialog(self, self.editor)
        if dialog.exec() != QDialog.Accepted:
            return None
        return self.open_pcr_product(dialog)

    def set_circular(self, circular):
        """Circular button / Edit → Circular Sequence. Making a plasmid linear asks where to
        break it (that base becomes base 1) and warns about features the break cuts."""
        ed = self.editor
        if circular or not ed.circular:
            ed.set_circular(circular)
            return True
        if not ed.dna_editable():   # making it linear moves base 1: a change of the bases
            self.circular_action.setChecked(True)
            return False
        n = ed.length
        base, ok = QInputDialog.getInt(
            self, "Make linear", f"Break the circle just before base (1–{n:,}):\n"
            "this base becomes base 1 (1 = keep the numbering).", ed.start + 1 if ed.start < n else 1, 1, n)
        cut = ed.features_cut_at(base - 1) if ok else []
        if ok and cut:
            names = ", ".join(f.label for f in cut[:8]) + ("…" if len(cut) > 8 else "")
            ok = QMessageBox.question(self, "Make linear", f"Breaking before base {base:,} cuts through: {names}.\n\n"
                                      "They stay as two-part features. Break there anyway?") == QMessageBox.Yes
        if not ok:
            self.circular_action.setChecked(True)  # still circular
            return False
        ed.linearize(base - 1)
        return True

    def _sync_display_actions(self):
        """The View menu ticks follow changes made in the display panel (and back)."""
        ed = self.editor
        self.enzymes_action.setChecked(ed.show_enzymes)
        self.translation_action.setChecked(ed.show_translation)
        self.complement_action.setChecked(ed.show_complement)
        self.cds_only_action.setChecked(ed.cds_only)

    def select_range(self, text=None):
        """The toolbar's range box (or `text`): select that stretch, or put the cursor after one base. A typing mistake is said
        in the status bar and the box keeps the text, to correct."""
        ed = self.editor
        if ed is None:
            return
        try:
            start, end = parse_range(self.range_box.text() if text is None else text, ed.length, ed.circular)
        except ValueError as e:
            self.statusBar().showMessage(str(e), 6000)
            return
        ed.select(start, end)
        ed.reveal.emit(ed.start)

    def copy_diagnostics(self):
        feedback.copy_diagnostics()
        self.statusBar().showMessage("Diagnostics copied to the clipboard (version, system, last log lines; no sequences).", 5000)

    def show_self_test(self):
        """Help > Troubleshooting > Self-Test: the checks that the editing, annotation and primer code is right."""
        from selftest import SelfTestPanel   # (imported here: it imports the tools)
        dialog = QDialog(self)
        dialog.setWindowTitle("Self-Test")
        dialog.resize(760, 520)
        layout = QVBoxLayout(dialog)
        layout.addWidget(SelfTestPanel())
        buttons = QDialogButtonBox(QDialogButtonBox.Close, rejected=dialog.reject)
        layout.addWidget(buttons)
        dialog.exec()

    def show_timings(self):
        """Help > Troubleshooting > Performance Timings: how long each step of drawing / editing takes."""
        if getattr(self, "timings", None) is None:
            self.timings = PerfWindow(self)
        self.timings.show()
        self.timings.raise_()

    def load_tracks(self, filters=TRACK_FILTER):
        """Import → tracks / reads: one or several files at once (Shift/Ctrl-click in the dialog)."""
        paths, _ = QFileDialog.getOpenFileNames(self, "Load tracks", "", filters)
        return [entry for entry in map(self.load_track, paths) if entry is not None]

    def load_track(self, path, seqid=None, quiet=False):
        """Add a GFF3/GTF/BED/VCF file as a track of the genome view (and to Find), or a
        BAM file as coverage + reads. `quiet`: no questions: None when the track has nothing for this sequence."""
        ed = self.editor
        bam = str(path).lower().endswith(".bam")
        try:
            track = run_busy(self, f"Reading {Path(path).name}…", lambda: (Bam if bam else Track)(str(path)))
        except OSError as e:
            QMessageBox.critical(self, "Cannot Read Track", f"{path}\n\n{e}")
            return None
        lengths = dict(track.sequences())
        seqs = list(lengths)
        if not seqs:
            QMessageBox.warning(self, "Empty track", f"No annotations in {path}.")
            return None
        # Which sequence of the track is this document? Given (a project), same name ("chr1" = "1"), or the only one, or ask.
        seqid = seqid if seqid in seqs else match_contig(seqs, ed.doc.name, ed.fasta_record)
        if seqid is None and len(seqs) == 1 and not quiet:
            seqid = seqs[0]   # the only one: shown here even if its name differs (warned below)
        if seqid is None and quiet:
            return None
        if seqid is None:
            seqid, ok = QInputDialog.getItem(self, "Track sequence",
                                             f"Which sequence of {Path(path).name} is {ed.doc.name}?", seqs, 0, False)
            if not ok:
                return None
        # Does it belong to this sequence? A warning only: the track is loaded anyway (a renamed or edited sequence is common).
        problems = track_check.mismatches(ed, None if bam else track, seqid, lengths[seqid] if bam else None, path)
        if problems and quiet:
            ed.error.emit(f"{Path(path).name} may not match this sequence: " + " ".join(problems))
        elif problems:
            QMessageBox.warning(self, "Track may not match this sequence",
                                "\n\n".join(problems) + "\n\nIt is loaded anyway; check that it is the right file.")
        if bam:
            entry = SimpleNamespace(name=Path(path).name, path=str(path), bam=track, seqid=seqid, length=lengths[seqid],
                                    cache=None, depth_pending=False, error=None, variants=None)
            ed.alignments.append(entry)
        else:
            entry = SimpleNamespace(name=Path(path).name, path=str(path), track=track, seqid=seqid)
            ed.tracks.append(entry)
            ed._search_text = None
        self.set_genome_visible(True)
        ed.display_changed.emit()
        return entry

    def map_reads(self):
        """Map reads to the current sequence with an external mapper, then show the BAM."""
        ed = self.editor
        name = runner.safe_name(ed.fasta_record or ed.doc.name)
        dialog = runner.MapReadsDialog(self, name)
        if dialog.exec() != QDialog.Accepted:
            return None
        folder = runner.run_folder(dialog.out.text())
        reference = folder / "ref.fa"
        try:
            folder.mkdir(parents=True, exist_ok=True)
            seq = ed.doc.sequence()  # current (edited) sequence, so reads line up with it
            run_busy(self, "Writing the reference…", lambda: reference.write_text(f">{name}\n{seq}\n"))
            steps, result = dialog.plan(reference)
        except (OSError, ValueError) as e:
            QMessageBox.critical(self, "Map Reads", str(e))
            return None
        if self.picker.count() > 1:   # a multi-sequence FASTA: the mapper only gets the open sequence
            self.statusBar().showMessage(f"Reads are mapped to {name} only: reads from the other "
                                         f"{plural(self.picker.count() - 1, 'sequence')} of the file stay unmapped.", 10000)
        run = RunDialog(self, f"Mapping reads to {name}", steps, folder)
        if run.exec() == QDialog.Accepted and run.ok:
            return self.load_track(str(result))
        return None

    def import_parts(self, paths=None):
        """Every named feature of your own files goes in the Features Database, as a Local feature."""
        if not paths:
            paths, _ = QFileDialog.getOpenFileNames(self, "Import features to the Features Database", "",
                                                    "Annotated sequences (*.gb *.gbk *.gbff *.dna);;All files (*)")
        if not paths:
            return 0
        rejected = []
        try:
            added = run_busy(self, "Reading features…", lambda: import_files(paths, rejected))
        except OSError as e:
            QMessageBox.critical(self, "Import features to the Features Database", str(e))
            return 0
        text = f"{plural(added, 'new feature')} added to the Features Database, source Local ({user_library()})."
        if rejected:   # named, with the reason, never dropped silently
            text += f"\n\n{plural(len(rejected), 'feature')} not added:\n" + "\n".join(f"• {n[:40]!r}: {why}" for n, why in rejected[:15])
            text += f"\n… and {len(rejected) - 15} more." if len(rejected) > 15 else ""
        (QMessageBox.warning if rejected else QMessageBox.information)(self, "Import features to the Features Database", text)
        return added

    def install_tools(self):
        if QMessageBox.question(
                self, "Install Mapping Tools",
                f"Create the conda environment '{runner.ENV}' with {', '.join(runner.TOOLS)} from bioconda "
                f"({'inside WSL, ' if runner.use_wsl() else ''}downloads a few hundred MB)?") == QMessageBox.Yes:
            RunDialog(self, "Installing mapping tools", [runner.INSTALL]).exec()

    def remove_tracks(self):
        """File > Remove Tracks: every annotation / variant track and read alignment loaded on this sequence (the files are not
        touched; it cannot be undone, so it asks first)."""
        ed = self.editor
        n_tracks, n_reads = len(ed.tracks), len(ed.alignments)
        if not n_tracks and not n_reads:
            return
        what = " and ".join(t for t in (plural(n_tracks, "track") if n_tracks else "", plural(n_reads, "read alignment")
                                        if n_reads else "") if t)
        if QMessageBox.question(self, "Remove Tracks", f"Remove {what} from this window? The files themselves are not changed.",
                                QMessageBox.Yes | QMessageBox.Cancel, QMessageBox.Cancel) != QMessageBox.Yes:
            return
        ed.tracks.clear()
        ed.alignments.clear()
        ed.display_changed.emit()

    def go_to_hit(self, record, start, end):
        """Genome-wide search result: open its sequence if needed, then select it."""
        ed = self.editor
        if record and record != ed.fasta_record and not self.load(ed.path, record):
            return
        self.editor.select(start, end)
        self.editor.reveal.emit(start)

    def autosave(self):
        """Write a recovery copy if there are changes since the last one. Never touches
        the user's file: that only happens on Save."""
        ed = self.editor
        if ed is None or not ed.doc.modified or ed.version == self.autosaved_version:
            return
        target = autosave_path(ed.path, ed.fasta_record)
        try:
            ed.doc.write_genbank(str(target), ed.file_primer_features())
        except OSError as e:
            self.statusBar().showMessage(f"Autosave failed: {e}", 5000)
            return
        if self.autosave_file not in (None, target):
            self._drop_autosave()  # file was renamed by Save As
        self.autosave_file, self.autosaved_version = target, ed.version

    def _drop_autosave(self):
        if self.autosave_file is not None:
            self.autosave_file.unlink(missing_ok=True)
            self.autosave_file = None

    def _confirm_discard(self):
        """True if it is OK to drop the current document (saved, or user agreed)."""
        away = [n for n, d in self.record_docs.items() if d.modified]   # other records of the file, edited earlier
        if self.editor is None or not (self.editor.doc.modified or away):
            return True
        answer = QMessageBox.question(
            self, "Unsaved changes", f"Save changes to {self.editor.path.name}?"
            + (f"\n\nAlso edited: {', '.join(away)}." if away else ""),
            QMessageBox.Save | QMessageBox.Discard | QMessageBox.Cancel)
        if answer == QMessageBox.Save:
            return self.save()
        return answer == QMessageBox.Discard

    def closeEvent(self, event):
        if self._confirm_discard():
            self._drop_autosave()
            self._save_window()
            event.accept()
        else:
            event.ignore()

    # ---- editing commands ----------------------------------------------------------
    def _update_undo_actions(self):
        """Edit > Undo / Redo say what they will undo or redo ("Undo Inserted 5 bp at 301"), greyed when there is nothing."""
        from history_panel import describe
        ed = self.editor
        if getattr(ed, "file_history", None) is not None:
            undo, redo = ed.file_history.steps(ed.file_histories())
            undo, redo = [step for _, step in undo], [step for _, step in redo]
        else:
            undo, redo = ed.doc.history()
        for a, verb, steps, step in ((self.undo_action, "&Undo", undo, undo[-1] if undo else None),
                                     (self.redo_action, "&Redo", redo, redo[0] if redo else None)):
            words = describe(step) if step is not None else ""
            words = words if len(words) <= 40 else words[:39] + "…"
            a.setText(f"{verb} {words.replace('&', '&&')}" if words else verb)   # (a "&" of the words is not a mnemonic)
            a.setEnabled(bool(steps))

    def clipboard_action(self, verb):
        """Edit > Cut / Copy / Paste / Select All on what has the keyboard focus: a text field or table its own way, the sequence
        (or any view of the document) through the editor, as its keys do."""
        w = self.focusWidget()
        own = {"cut": "cut", "copy": "copy", "paste": "paste", "select_all": "selectAll"}[verb]
        if w is not None and not isinstance(w, SequenceView):
            target = w.lineEdit() if hasattr(w, "lineEdit") and callable(w.lineEdit) and w.lineEdit() is not None else w
            if callable(getattr(target, own, None)):
                getattr(target, own)()
                return
        ed = self.editor
        if ed is None:
            return
        clipboard = QGuiApplication.clipboard()
        if verb == "copy":
            clipboard.setMimeData(ed.copy_mime())
        elif verb == "cut" and ed.end > ed.start:
            clipboard.setMimeData(ed.copy_mime())
            ed.backspace()
        elif verb == "paste":
            ed.paste_mime(clipboard.mimeData())
        elif verb == "select_all":
            ed.select(0, ed.length)

    def copy_reverse_complement(self):
        QGuiApplication.clipboard().setMimeData(self.editor.copy_reverse_complement_mime())

    def edit_current_feature(self):
        if self.editor.current_feature is None:
            self.statusBar().showMessage("Click a feature first.", 3000)
        else:
            edit_feature(self, self.editor, self.editor.current_feature)

    def delete_current_feature(self):
        if self.editor.current_feature is None:
            self.statusBar().showMessage("Click a feature first.", 3000)
        else:
            delete_feature(self, self.editor, self.editor.current_feature)  # a gene asks about its parts

    def group_as_gene(self):
        """Features ticked in the feature table, else those inside the selected range."""
        ed = self.editor
        tables = self.workspace.views(FeatureTable)
        chosen = tables[0].selected_indices() if tables else []
        if len(chosen) < 2 and ed.end > ed.start:
            chosen = features_inside(ed)
        if not chosen:
            self.statusBar().showMessage("Select the parts first (feature table rows, or a range covering them).", 5000)
            return None
        return group_as_gene(self, ed, chosen)

    def focus_find(self):
        self.search_bar.focus()

    def find_next(self):
        self.search_bar.step(1)

    def find_previous(self):
        self.search_bar.step(-1)

    def go_to(self):
        ed = self.editor
        pos, ok = QInputDialog.getInt(self, "Go to Position", f"Base (1–{ed.length:,}):",
                                      ed.start + 1, 1, max(1, ed.length))
        if ok:
            ed.set_cursor(pos - 1)
            ed.reveal.emit(pos - 1)
            self.sequence_view.setFocus()

    # ---- status --------------------------------------------------------------------
    def _update_title(self):
        if self.editor is None:
            self.setWindowTitle(project.NAME)
            return
        star = "*" if self.editor.doc.modified else ""
        self.circular_action.setChecked(self.editor.doc.circular)
        self.setWindowTitle(f"{star}{self.editor.path.name} ({self.editor.length:,} bp) — {project.NAME}")
        self._update_status()

    def _update_status(self):
        ed = self.editor
        if ed.end > ed.start:
            gc = ed.gc_percent()
            end = ed.end - ed.length if ed.wraps else ed.end
            text = f"{ed.start + 1:,}–{end:,}  ({ed.end - ed.start:,} bp{' across the origin' if ed.wraps else ''}, GC {gc:.1f}%"
            t = melting_temp(ed.selected_sequence()) if ed.end - ed.start <= 60 else None
            self.status.setText(text + (f", Tm {t:.1f} °C)" if t is not None else ")"))
        else:
            where = "Cursor at start" if ed.start == 0 else f"Cursor after base {ed.start:,}"
            self.status.setText(f"{where}  ·  {ed.length:,} bp")


WINDOWS = []  # extra windows (PCR products...): kept referenced so they stay open


def open_window(doc, path, untitled=False):
    win = MainWindow()
    WINDOWS.append(win)
    win.show()
    win.show_document(doc, path, untitled)
    return win


WELCOME_FLAG = "--welcome"


def start_up_questions(win, force_welcome=False):
    """Over the open window, one after the other (never three dialogs at once): the welcome page (always with force_welcome, the
    --welcome option), then the question about opening sequence files with this program, then (when updates are switched on) the
    silent update check."""
    def then_associations(*_):
        question = file_assoc.ask_at_startup(win)
        if question is not None:
            question.finished.connect(then_updates)
        else:
            then_updates()

    def then_updates(*_):
        updater.startup_check(win, **win.update_hooks())

    welcome = show_at_startup(win, force=force_welcome)
    if welcome is not None:
        welcome.finished.connect(then_associations)
    else:
        then_associations()


# Files the program reads from its own folder at run time (the DATA of packaging/openDNA.spec): a build that misses one fails on that
# feature only, so the smoke test asks the program itself whether it finds them.
BUNDLED_FILES = ("parts/ncbi_vectors.tsv", "parts/engineering_proteins.tsv", "tools.toml", "sponsors.json", "vendor/3Dmol-min.js",
                 "vendor/3Dmol.LICENSE", "logo/icon_256.png", "logo/wordmark.png", "logo/wordmark_dark.png")


def write_smoke_report(win, path):
    """For the smoke test of a build: the window's title, the file and sequence it shows, and the bundled files it cannot find, as
    JSON in `path`. (Only when OPENDNA_SMOKE_REPORT is set: nothing is written otherwise.)"""
    import json
    here = Path(__file__).parent
    ed = win.editor
    report = {"title": win.windowTitle(), "file": str(ed.path) if ed is not None else None,
              "sequence": ed.doc.name if ed is not None else None, "length": ed.length if ed is not None else 0,
              "missing": [f for f in BUNDLED_FILES if not (here / f).exists()], "frozen": bool(getattr(sys, "frozen", False)),
              "platform": sys.platform, "qt": qVersion()}
    tmp = Path(path).with_suffix(".tmp")
    tmp.write_text(json.dumps(report, indent=1), encoding="utf-8")
    fileio.replace(tmp, path)


def app_icon():
    """The window and taskbar icon: every size drawn for itself (packaging/make_icon.py makes app/logo/icon_*.png)."""
    icon = QIcon()
    for png in Path(__file__).with_name("logo").glob("icon_*.png"):
        icon.addFile(str(png))
    return icon


if __name__ == "__main__":
    if sys.platform == "win32":   # own taskbar button and icon, also when run from source (otherwise it groups under python.exe)
        import ctypes
        ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(project.NAME)
    app = QApplication(sys.argv)
    app.setWindowIcon(app_icon())
    app.setApplicationName(project.NAME)   # the data folder is then %APPDATA%\openDNA whether run from source or as the exe
    copy_old_data()   # settings and primers kept under the old name come along, before anything reads the data folder
    feedback.setup_logging()
    theme.apply(theme.saved(), remember=False)
    win = MainWindow()
    win.show()
    # Options: --welcome opens the start page even if it was shown today or turned off. Anything else is a file to open.
    files = [a for a in sys.argv[1:] if a != WELCOME_FLAG]
    QTimer.singleShot(0, lambda: start_up_questions(win, force_welcome=WELCOME_FLAG in sys.argv[1:]))
    if files:
        win.load(files[0])
    if os.environ.get("OPENDNA_SMOKE_REPORT"):   # the smoke test of a build (scripts/smoke_unix.py) reads what the program shows
        QTimer.singleShot(1500, lambda: write_smoke_report(win, os.environ["OPENDNA_SMOKE_REPORT"]))
    code = app.exec()
    feedback.teardown_logging()
    sys.exit(code)
