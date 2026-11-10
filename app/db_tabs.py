"""Databases → Features Database / Primer Database: a tab of the window (like the Assembly tab: it works with no file open, and
closing it only takes it out, so reopening it finds it as it was).

Each tab is the program's own views: on the left the Features view (feature_table.py) or the Primers view (primer_panel.py)
listing the whole database, on the right the sequence view (read-only: select, scroll, copy) of the entry clicked, alone: a
feature on its own bases with the feature over all of them, a primer's own bases with the primer drawn on them.

The list needs a document to read its rows from: every feature (every primer) of the database one after the other, 20 N apart. It
is never shown. Its editor carries the database (`editor.database`): the columns the table shows (Source, with the record or
file, instead of Location for the features; Note and Added for the primers), and the right-click menu of the database instead of that of a file. It
is made again when the database changes (a feature added, a primer edited).

What is filtered and how (right-click a header) is in value_filter.py.
"""
from pathlib import Path

from PySide6.QtCore import Qt, QTimer
from PySide6.QtGui import QKeySequence, QShortcut
from PySide6.QtWidgets import (
    QDialog, QDialogButtonBox, QFrame, QHBoxLayout, QLabel, QMenu, QMessageBox, QPushButton, QSplitter, QVBoxLayout, QWidget,
)

from opendna_core import Document

import db_status
import parts_library
import primer_db
import runner
from editor import Editor
from uikit import fit_columns, plural

SPACER = 20   # N between two entries on the database's (unseen) document
FEATURES, PRIMERS = "features", "primers"
TITLES = {FEATURES: "Features Database", PRIMERS: "Primer Database"}
# Databases menu: (label, key)
MENU = (("&Features Database…", FEATURES), ("&Primer Database…", PRIMERS))


def concatenated(name, sequences):
    """(Document of the sequences one after the other, SPACER N apart, [(start, end)] of each). An empty database: one N."""
    spans, at = [], 0
    for seq in sequences:
        spans.append((at, at + len(seq)))
        at += len(seq) + SPACER
    doc = Document.new_sequence(name, ("N" * SPACER).join(sequences) or "N", False)
    return doc, spans


def part_document(part):
    """A feature of the database on its own: its bases, with the feature over all of them."""
    doc = Document.new_sequence(part.name, part.sequence, False)
    notes = [("label", part.name), ("note", f"Features Database, source: {part.source_name}")]
    if part.source:
        notes.append(("note", f"from {part.source}"))
    doc.add_features([(part.kind, f"1..{len(part.sequence)}", notes)])
    doc.modified = False
    return doc


def primer_document(row):
    """A primer of the database on its own bases (the sequence view draws the primer on them)."""
    doc = Document.new_sequence(row["name"], "".join(row["sequence"].split()), False)
    doc.modified = False
    return doc


# ---- the Features Database ---------------------------------------------------------------------------------------------------

class FeaturesDatabase:
    """`editor.database` of the Features Database's document: which part each feature is, the table's columns, the menu."""
    key, title = FEATURES, TITLES[FEATURES]

    def __init__(self):
        self.parts = parts_library.all_parts()
        self.measured = any(p.strength is not None for p in self.parts)

    def document(self):
        doc, spans = concatenated(self.title, [p.sequence for p in self.parts])
        first = doc.add_features([(p.kind, f"{s + 1}..{e}", [("label", p.name)]) for p, (s, e) in zip(self.parts, spans)]) \
            if self.parts else 0
        doc.modified = False
        self.by_index = {first + i: p for i, p in enumerate(self.parts)}
        return doc

    def columns(self):
        """The Features view's columns here: where a feature comes from instead of where it is on the (unseen) document."""
        return ["Name", "Type", "Source", "Length"] + (["Measured"] if self.measured else [])

    def extra(self):
        """Source: the source's name and the record or file the feature was taken from, in one field."""
        return {i: {"Source": parts_library.source_label(p), "Measured": parts_library.describe(p)} for i, p in self.by_index.items()}

    def menu(self, menu, table, indices):
        local = [i for i in indices if self.by_index.get(i) is not None and self.by_index[i].collection == parts_library.YOURS]
        delete = menu.addAction(f"Delete {plural(len(local), 'Local Feature')} from the Features Database" if local
                                else "Delete from the Features Database", lambda: self.delete(table, indices))
        delete.setEnabled(bool(local))
        delete.setToolTip("Only your own (Local) features can be deleted; a published source is switched off whole (Sources)")
        menu.setToolTipsVisible(True)

    def delete(self, parent, indices):
        parts = [self.by_index[i] for i in indices if i in self.by_index and self.by_index[i].collection == parts_library.YOURS]
        if not parts:
            return
        names = ", ".join(p.name for p in parts[:5]) + ("…" if len(parts) > 5 else "")
        if QMessageBox.question(parent, "Delete features", f"Delete {plural(len(parts), 'Local feature')} from the Features Database?"
                                f"\n\n{names}", QMessageBox.Yes | QMessageBox.No, QMessageBox.No) == QMessageBox.Yes:   # No by default
            parts_library.remove_parts([(p.name, p.sequence) for p in parts])


def sources_menu(parent):
    """A menu of the sources, each ticked when Detect Common Features and the Plasmid Checkup use it."""
    menu = QMenu(parent)

    def fill():
        menu.clear()
        menu.addSection("Sources used by Detect Features and the Checkup")
        for c in parts_library.collections():
            a = menu.addAction(f"{c.title}  ({sum(c.counts.values())})")
            a.setCheckable(True)
            a.setChecked(c.on)
            a.setToolTip("\n".join(x for x in (c.reference, c.doi and f"doi:{c.doi}", c.note) if x) or c.title)
            a.toggled.connect(lambda on, cid=c.id: parts_library.set_on(cid, on))
        if menu.isEmpty():
            menu.addAction("No source").setEnabled(False)
    menu.aboutToShow.connect(fill)
    menu.setToolTipsVisible(True)
    return menu


# ---- the Primer Database -----------------------------------------------------------------------------------------------------

class PrimersDatabase:
    """`editor.database` of the Primer Database's document."""
    key, title = PRIMERS, TITLES[PRIMERS]

    def document(self):
        doc, _ = concatenated(self.title, ["".join(r["sequence"].split()) for r in primer_db.db().rows])
        doc.modified = False
        return doc

    tab = None   # the DatabaseTab showing it

    def edit(self, name):
        self.tab.edit_primer(name)

    def menu(self, menu, panel, names):
        known = [n for n in names if primer_db.db().find(n) is not None]
        edit = menu.addAction("Edit…", lambda: self.edit(known[0]))
        edit.setEnabled(len(known) == 1)
        delete = menu.addAction(f"Delete {len(known)} from Primer Database" if len(known) > 1 else "Delete from Primer Database",
                                lambda: self.delete(panel, known))
        delete.setEnabled(bool(known))

    def delete(self, parent, names):
        if names and QMessageBox.question(parent, "Delete primers", f"Delete {plural(len(names), 'primer')} from the database?\n\n"
                                          + ", ".join(names[:5]) + ("…" if len(names) > 5 else ""),
                                          QMessageBox.Yes | QMessageBox.No, QMessageBox.No) == QMessageBox.Yes:   # No by default
            primer_db.remove_primers(None, names)


# ---- the tab ---------------------------------------------------------------------------------------------------------------

HELP = {
    FEATURES: """<p>Every known feature the program can look for: the published <b>sources</b> that ship with it, and your own
features, the source <b>Local</b>.</p>
<ul>
<li><b>Click a feature</b> to see its sequence on the right.</li>
<li><b>Sort</b>: click a column header. <b>Filter</b>: type in the box above the list (name, type, source...), or <b>right-click a
column header</b> (Type, Strand, Source) to tick the values shown; a filtered column ends with ▼.</li>
<li><b>Source</b> is where a feature comes from (a publication or database, or Local), then the record or file it was taken
from. <b>Measured</b>, when a source measured something, is its strength in that source's units and host.</li>
<li><b>Sources</b> ticks the sources that Detect Common Features and the Plasmid Checkup look for.</li>
<li><b>Add</b> your own: right-click features in a file → Add to Features Database, or Import → Features and Tracks → Import
Features to Features Database. <b>Delete</b> them: right-click here (only Local features can be deleted).</li>
</ul>""",
    PRIMERS: """<p>Every primer you have made, designed or imported, the same in every file: a primer that binds the open sequence is
drawn there.</p>
<ul>
<li><b>Click a primer</b> to see it on the right.</li>
<li><b>Sort</b>: click a column header. <b>Filter</b>: <b>right-click a column header</b> (Match, Structure, Added) to tick the
values shown; a filtered column ends with ▼.</li>
<li><b>Open a primer</b> (its ▸) to see where it binds in the background genome (the toolbar's genome box): the whole primer, or
only its 3′ end. In each site the primer's bases that do not pair with the genome are <b style="color:#d32f2f">red</b>.</li>
<li><b>Structure</b> is the hairpin and self-dimer check.</li>
<li><b>New / Edit / Delete</b> change the database (double-click a primer to edit it). <b>Import / Export</b> read and write TSV,
CSV or FASTA lists. <b>Options</b> are the PCR conditions and the checks.</li>
</ul>""",
}


def help_dialog(parent, key):
    """The tab's Help: how to use it, then where the database is kept and what is wrong in it (db_status)."""
    check = db_status.features_library if key == FEATURES else db_status.primers
    st = check()
    dialog = QDialog(parent)
    dialog.setWindowTitle(f"{TITLES[key]} Help")
    layout = QVBoxLayout(dialog)
    layout.addWidget(QLabel("<h3>How to use it</h3>", textFormat=Qt.RichText))
    layout.addWidget(QLabel(HELP[key], wordWrap=True, textFormat=Qt.RichText))
    layout.addWidget(QFrame(frameShape=QFrame.HLine, frameShadow=QFrame.Sunken))
    layout.addWidget(QLabel("<h3>Where it is kept, and its state</h3>", textFormat=Qt.RichText))
    head, facts = db_status.status_widgets(st)
    layout.addWidget(head)
    layout.addWidget(facts)
    layout.addStretch(1)
    buttons = QDialogButtonBox(QDialogButtonBox.Close)
    buttons.rejected.connect(dialog.reject)
    if st.path is not None:
        buttons.addButton(QPushButton("Open Folder", clicked=lambda: db_status.open_folder(st)), QDialogButtonBox.ActionRole)
    layout.addWidget(buttons)
    dialog.resize(680, dialog.sizeHint().height())
    dialog.status = st   # for the tests
    return dialog


class DatabaseTab(QWidget):
    """A database's tab: a header (its name, how many entries, its buttons, Help), its list (the Features or Primers view) and
    the sequence of the entry clicked."""

    def __init__(self, key):
        super().__init__()
        self.key, self.title = key, TITLES[key]
        self.editor = self.table = self.sequence = self.entry_editor = None
        self.entry = None
        title = QLabel(self.title)
        font = title.font()
        font.setBold(True)
        font.setPointSizeF(font.pointSizeF() * 1.25)
        title.setFont(font)
        self.count = QLabel()
        self.count.setEnabled(False)   # greyed: a detail beside the title
        head = QHBoxLayout()
        head.setContentsMargins(8, 6, 8, 6)
        head.setSpacing(6)
        head.addWidget(title)
        head.addWidget(self.count)
        head.addStretch()
        self.edit_button = self.delete_button = None
        if key == FEATURES:
            # a push button with a menu: its arrow has room of its own (a tool button drew it over the text)
            sources = QPushButton("Sources", toolTip="Which sources Detect Common Features and the Plasmid Checkup look for")
            sources.setMenu(sources_menu(sources))
            head.addWidget(sources)
        else:
            for text, tip, slot in (("New…", "A new primer", self.new_primer),
                                    ("Edit…", "Edit the primer clicked", lambda: self.edit_primer(self._primer_name())),
                                    ("Delete", "Delete the primers selected (Delete key)", self.delete_primers),
                                    ("Import…", "Add the primers of a TSV, CSV or FASTA list", self.import_primers),
                                    ("Export…", "Write the whole database as TSV, CSV or FASTA", self.export_primers),
                                    ("Options…", "PCR conditions and the primer checks", self.primer_options)):
                button = QPushButton(text, clicked=slot, toolTip=tip)
                head.addWidget(button)
                if text == "Edit…":
                    self.edit_button = button
                elif text == "Delete":
                    self.delete_button = button
        separator = QFrame(frameShape=QFrame.VLine, frameShadow=QFrame.Sunken)
        head.addWidget(separator)
        head.addWidget(QPushButton("Help", clicked=lambda: help_dialog(self, key).exec(),
                                   toolTip="How to use this tab, and where the database is kept"))
        self.caption = QLabel(wordWrap=True)
        self.caption.setContentsMargins(6, 4, 6, 4)
        self.preview = QWidget()   # the sequence of the entry clicked (or a hint)
        self.preview_box = QVBoxLayout(self.preview)
        self.preview_box.setContentsMargins(0, 0, 0, 0)
        self.preview_box.setSpacing(0)
        self.preview_box.addWidget(self.caption)
        self.hint = QLabel(f"Click a {'feature' if key == FEATURES else 'primer'} to see its sequence.", alignment=Qt.AlignCenter)
        self.hint.setEnabled(False)
        self.preview_box.addWidget(self.hint, 1)
        self.split = QSplitter(Qt.Horizontal)
        self.split.setChildrenCollapsible(False)
        self.split.addWidget(self.preview)   # the list goes before it (rebuild)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)
        layout.addLayout(head)
        layout.addWidget(QFrame(frameShape=QFrame.HLine, frameShadow=QFrame.Sunken))
        layout.addWidget(self.split, 1)
        # The database changed (often several times in a row: an import): the list is made again once, after.
        self._later = QTimer(self, singleShot=True, interval=0, timeout=self.rebuild)
        (parts_library.notifier if key == FEATURES else primer_db.notifier).changed.connect(self._later.start)
        self.rebuild()
        self.split.setSizes([6500, 3500])   # the list first: the right side is empty until an entry is clicked

    # ---- the list ----------------------------------------------------------------------------------------------------------
    def rebuild(self):
        """Make the list again from the database, keeping the filters and the entry clicked."""
        from feature_table import FeatureTable
        from primer_panel import PrimerPanel
        keep = self._state()
        database = FeaturesDatabase() if self.key == FEATURES else PrimersDatabase()
        database.tab = self
        doc = database.document()
        ed = Editor(doc, Path(f"{self.title}.gb"))
        ed.untitled, ed.dna_locked, ed.database = True, True, database
        if self.key == FEATURES:
            ed.feature_columns, ed.feature_extra = database.columns(), database.extra()
            ed.set_display(primers=False)
        old = self.table
        sizes = self.split.sizes()
        self.editor = ed
        self.table = FeatureTable(ed) if self.key == FEATURES else PrimerPanel(ed)
        if old is not None:
            old.setParent(None)
            old.deleteLater()
        self.split.insertWidget(0, self.table)
        if old is not None:
            self.split.setSizes(sizes)
        view = self.view()
        if self.key == FEATURES:
            # Columns fitted to what they hold; Source (a publication and its record) takes the room left, up to its content.
            fit_columns(view, stretch=database.columns().index("Source"))
            view.selectionModel().currentRowChanged.connect(lambda current, _=None: self._feature_row(current.row()))
            n, sources = len(database.parts), len({p.collection for p in database.parts})
            self.count.setText(f"{plural(n, 'feature')} from {plural(sources, 'source')}")
        else:
            view.currentItemChanged.connect(lambda item, _=None: self._primer_item(item))
            view.itemSelectionChanged.connect(self._can_change)
            # The Delete key deletes the primers selected, as the button does (only while the list has the focus: not while
            # the sequence on the right is being read).
            QShortcut(QKeySequence.Delete, view, self.delete_primers, context=Qt.WidgetShortcut)
            n = len(primer_db.db().rows)
            self.count.setText(plural(n, "primer"))
        self._restore(keep)
        self._can_change()

    def view(self):
        """The table widget of the list (for the filters)."""
        if self.table is None:
            return None
        return self.table.view if self.key == FEATURES else self.table.tree

    def _state(self):
        view = self.view()
        if view is None:
            return None
        text = self.table.filter.text() if self.key == FEATURES else ""
        sort = (self.table.model.sort_column, self.table.model.sort_order) if self.key == FEATURES else self.table.sort
        return dict(view.value_filter.hidden), text, self.entry, sort

    def _restore(self, keep):
        if keep is None:
            return
        hidden, text, entry, sort = keep
        view = self.view()
        if self.key == FEATURES and text:
            self.table.filter.setText(text)
        if self.key == FEATURES:   # the column sorted by stays
            view.horizontalHeader().setSortIndicator(*sort)
        elif sort is not None:
            self.table.sort = sort
            self.table._apply_sort()
        view.value_filter.hidden = hidden
        view.value_filter.apply()
        if entry is None:
            return
        if self.key == FEATURES:
            by_index, wanted = self.editor.database.by_index, (entry.collection, entry.name, entry.sequence)
            row = next((r for r, shown in enumerate(self.table.model.shown) if shown[3] in by_index and
                        (by_index[shown[3]].collection, by_index[shown[3]].name, by_index[shown[3]].sequence) == wanted), None)
            if row is not None:
                view.selectRow(row)
            else:
                self.show_entry(None)   # deleted
        else:
            root = view.invisibleRootItem()
            item = next((root.child(i) for i in range(root.childCount()) if root.child(i).data(0, Qt.UserRole) == entry["name"]), None)
            if item is not None:
                view.setCurrentItem(item)
                self.show_entry(primer_db.db().find(entry["name"]) or None)
            else:
                self.show_entry(None)

    # ---- the entry clicked -------------------------------------------------------------------------------------------------
    def _feature_row(self, row):
        model = self.table.model
        if 0 <= row < len(model.shown):
            part = self.editor.database.by_index.get(model.shown[row][3])
            if part is not None:
                self.show_entry(part)

    def _primer_item(self, item):
        from primer_panel import primer_item
        if item is not None:
            row = primer_db.db().find(primer_item(item).data(0, Qt.UserRole))
            if row is not None:
                self.show_entry(row)

    def show_entry(self, entry):
        """The sequence of `entry` (a part or a primer row) alone on the right; None: the hint."""
        from assembly_panel import ReadOnlyView
        same = entry is not None and self.entry is not None and self.sequence is not None and (
            (entry.name, entry.sequence) == (self.entry.name, self.entry.sequence) if self.key == FEATURES
            else (entry["name"], entry["sequence"]) == (self.entry["name"], self.entry["sequence"]))
        self.entry = entry
        if same:
            return
        if self.sequence is not None:
            self.sequence.setParent(None)
            self.sequence.deleteLater()
            self.sequence = self.entry_editor = None
        self.hint.setVisible(entry is None)
        if entry is None:
            self.caption.clear()
            return
        doc = part_document(entry) if self.key == FEATURES else primer_document(entry)
        name = entry.name if self.key == FEATURES else entry["name"]
        ed = Editor(doc, Path(runner.safe_name(name) + ".gb"))
        ed.untitled, ed.dna_locked = True, True
        if self.key == FEATURES:
            ed.set_display(primers=False)
            measured = parts_library.describe(entry)
            self.caption.setText(f"<b>{_html(entry.name)}</b> · {_html(entry.kind)} · {len(entry.sequence):,} bp · "
                                 f"{_html(entry.source_name)}" + (f" ({_html(entry.source)})" if entry.source else "")
                                 + (f" · measured {_html(measured)}" if measured else ""))
        else:
            seq = "".join(entry["sequence"].split())
            self.caption.setText(f"<b>{_html(entry['name'])}</b> · {len(seq)} nt" + (f" · {_html(entry['note'])}" if entry["note"] else "")
                                 + (f" · added {_html(entry['added'])}" if entry["added"] else ""))
        self.entry_editor = ed
        self.sequence = ReadOnlyView(ed)
        self.preview_box.addWidget(self.sequence, 1)

    # ---- the Primer Database's buttons -------------------------------------------------------------------------------------
    def _can_change(self):
        """Edit needs one primer clicked, Delete at least one selected: otherwise they are greyed."""
        if self.key != PRIMERS or self.table is None:
            return
        chosen = self.table.chosen()
        self.edit_button.setEnabled(self._primer_name() is not None and len(chosen) <= 1 and bool(chosen))
        self.delete_button.setEnabled(bool(chosen))

    def _primer_name(self):
        from primer_panel import primer_item
        item = self.table.tree.currentItem() if self.key == PRIMERS else None
        return primer_item(item).data(0, Qt.UserRole) if item is not None else None

    def _dialog_editor(self):
        """The editor a primer dialog checks against: the primer on view, else the list's."""
        return self.entry_editor or self.editor

    def new_primer(self):
        from primer_panel import PrimerDialog
        ed = self._dialog_editor()
        ed.set_cursor(0)   # no selection: an empty new primer, not a copy of the one on view
        PrimerDialog(self, ed).exec()

    def edit_primer(self, name):
        from primer_panel import PrimerDialog
        if name and primer_db.db().find(name) is not None:
            PrimerDialog(self, self._dialog_editor(), name).exec()

    def delete_primers(self):
        chosen = self.table.chosen()
        if chosen:
            self.editor.database.delete(self, chosen)

    def import_primers(self):
        from primer_panel import import_file
        import_file(self, None)

    def export_primers(self):
        from primer_panel import export_file
        export_file(self)

    def primer_options(self):
        from primer_panel import open_options
        if open_options(self):
            self.table.reload()


def _html(text):
    from html import escape
    return escape(str(text))


def open_tab(window, key):
    """Databases → Features / Primer Database: its tab in `window` (made the first time, shown again after)."""
    return window.show_tool_tab(TITLES[key], lambda: DatabaseTab(key))
