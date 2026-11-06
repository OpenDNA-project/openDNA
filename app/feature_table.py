"""Sortable, filterable feature list. Selecting a row selects the feature everywhere.

Gene groups (see genes.py) show as one row per gene; click its ▸ to list the parts below
it (collapsed by default). Right-click: group selected features as a gene, ungroup, delete.

Rows are computed once per document change and sorted/filtered as plain Python lists:
Qt's QSortFilterProxyModel would call data() ~70 times per row while sorting, which took
seconds on a 10k-feature genome.

Right-click a header (Type, Strand, and the database's Source) to show only some of its values (value_filter.py).
The Features Database tab (db_tabs.py) uses this same table on a document made of the library's features: its editor names
the columns (`feature_columns`, e.g. Source and From instead of Location) and their texts (`feature_extra`), and
`database` makes the right-click menu the database's own.
"""
from PySide6.QtCore import QAbstractTableModel, QEvent, QModelIndex, Qt
from PySide6.QtGui import QFont, QGuiApplication, QKeySequence
from PySide6.QtWidgets import (
    QAbstractItemView, QInputDialog, QLineEdit, QMenu, QMessageBox, QTableView, QVBoxLayout, QWidget,
)

import genes
import validate
import value_filter
from uikit import ask, fit_columns, plural
from annotation import add_parts, feature_part
from parts_library import source_of
from editor import arrow, feature_key, reverse_complement
from primer_db import db, remove_primers

COLUMNS = ["Name", "Type", "Location", "Length", "Strand"]
FILTERED = ["Type", "Strand", "Source"]   # the columns with a filter by values (categories; the others are names / numbers)
TOGGLE_PX = 18  # clicks this close to the left of a gene's name expand/collapse it


class FeatureModel(QAbstractTableModel):
    fit_quiet = False   # True during a change the table's column fitting may skip (see uikit._Fitter)

    def __init__(self, editor):
        super().__init__()
        self.editor = editor
        self.columns = list(getattr(editor, "feature_columns", None) or COLUMNS)
        self.marked = set()   # columns filtered by values: their header shows the mark
        self.sort_column = self.columns.index("Location") if "Location" in self.columns else 0
        self.sort_order, self.filter_text = Qt.AscendingOrder, ""
        self.expanded = set()  # keys of open gene groups
        self.bold = QFont()
        self.bold.setBold(True)
        self.rows = []      # top level (genes, ungrouped features): (texts, sort keys, search, index, group key)
        self.children = {}  # group key -> rows of its parts, in gene order
        self.shown = []     # rows after filter + sort, parts under their open gene
        self.reload()

    def _row(self, f, group=None, part=False):
        length = sum(e - s for s, e in f.spans)
        name = f.label
        if group is not None:
            g = self.editor.groups[group]
            name = f"{'▾' if group in self.expanded else '▸'} {f.label}  ({len(g.parts)} parts)"
        elif part:
            name = f"      {f.label}"
        way = arrow(f)
        base = {"Name": (name, f.label.lower()), "Type": (f.kind, f.kind), "Location": (f.location, f.spans[0][0]),
                "Length": (f"{length:,}", length), "Strand": ({1: "→", -1: "←", 0: "—"}[way], way)}
        extra = getattr(self.editor, "feature_extra", {}).get(f.index, {})
        texts = tuple(base[c][0] if c in base else extra.get(c, "") for c in self.columns)
        keys = tuple(base[c][1] if c in base else extra.get(c, "").lower() for c in self.columns)
        search = " ".join((f.label, f.kind, f.location) + tuple(extra.values())).lower()
        return (texts, keys, search, f.index, group)

    def reload(self):
        self.beginResetModel()
        ed = self.editor
        self.rows = [self._row(f, master_key(ed, f)) for f in ed.features if f.index not in ed.part_of]
        self.children = {k: [self._row(p, part=True) for p in g.parts] for k, g in ed.groups.items()}
        self._arrange()
        self.endResetModel()

    def _arrange(self):
        needle = self.filter_text
        shown = []
        top = sorted(self.rows, key=lambda r: r[1][self.sort_column], reverse=self.sort_order == Qt.DescendingOrder)
        for r in top:
            kids = self.children.get(r[4], []) if r[4] else []
            hits = [k for k in kids if needle in k[2]] if needle else []
            if needle and needle not in r[2] and not hits:
                continue
            shown.append(r)
            if r[4] in self.expanded:
                shown += kids
            elif hits:  # a filter matching a part shows that part under its gene
                shown += hits
        self.shown = shown

    def toggle(self, row):
        key = self.shown[row][4]
        if key:
            self.set_expanded(key, key not in self.expanded)

    def set_expanded(self, key, open_):
        """Open or close one gene: only its row (the arrow) and the rows of its parts under it change. (It rebuilt and sorted the
        rows of every feature and reset the table, then the columns were measured again: most of a second on a genome.)"""
        if (key in self.expanded) == open_:
            return
        self.expanded ^= {key}
        at = next((i for i, r in enumerate(self.shown) if r[4] == key), None)
        if at is None:   # the gene is not on view (filtered out): nothing to redraw
            return
        old = self.shown[at]
        new = self._row(self.editor.by_index[old[3]], key)   # its name now shows ▾ / ▸
        self.rows[next(i for i, r in enumerate(self.rows) if r is old)] = new
        kids = self.children.get(key, [])
        ids = {id(k) for k in kids}
        end = at + 1
        while end < len(self.shown) and id(self.shown[end]) in ids:   # the parts shown under it now
            end += 1
        needle = self.filter_text
        under = kids if open_ else [k for k in kids if needle in k[2]] if needle else []
        self.fit_quiet = True   # (the parts hold what their gene holds: the columns need not be measured again)
        try:
            self.shown[at] = new
            self.dataChanged.emit(self.index(at, 0), self.index(at, len(self.columns) - 1))
            if end > at + 1:
                self.beginRemoveRows(QModelIndex(), at + 1, end - 1)
                del self.shown[at + 1:end]
                self.endRemoveRows()
            if under:
                self.beginInsertRows(QModelIndex(), at + 1, at + len(under))
                self.shown[at + 1:at + 1] = under
                self.endInsertRows()
        finally:
            self.fit_quiet = False

    def set_filter(self, text):
        self.beginResetModel()
        self.filter_text = text.lower()
        self._arrange()
        self.endResetModel()

    def sort(self, column, order=Qt.AscendingOrder):
        self.layoutAboutToBeChanged.emit()
        self.sort_column, self.sort_order = column, order
        self._arrange()
        self.layoutChanged.emit()

    def row_of(self, feature_index):
        """Row of a feature; a part of a closed gene opens it first."""
        key = self.editor.part_of.get(feature_index)
        if key and key not in self.expanded and not any(r[3] == feature_index for r in self.shown):
            self.set_expanded(key, True)
        return next((i for i, r in enumerate(self.shown) if r[3] == feature_index), None)

    def rowCount(self, parent=QModelIndex()):
        return 0 if parent.isValid() else len(self.shown)

    def columnCount(self, parent=QModelIndex()):
        return len(self.columns)

    def headerData(self, section, orientation, role=Qt.DisplayRole):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            name = self.columns[section]
            return name + value_filter.MARK if name in self.marked else name
        return None

    def data(self, index, role=Qt.DisplayRole):
        texts, _, _, feature_index, group = self.shown[index.row()]
        if role == Qt.DisplayRole:
            return texts[index.column()]
        if role == Qt.ToolTipRole and len(str(texts[index.column()])) > 24:   # what an elided cell holds, in full
            return str(texts[index.column()])
        if role == Qt.TextAlignmentRole and self.columns[index.column()] == "Length":   # numbers line up on their last digit
            return int(Qt.AlignRight | Qt.AlignVCenter)
        if role == Qt.DecorationRole and index.column() == 0:
            return self.editor.colors[feature_index]
        if role == Qt.FontRole and group:
            return self.bold
        if role == Qt.ForegroundRole and self._hidden(feature_index):   # hidden on the sequence: greyed (right-click > Show)
            return QGuiApplication.palette().placeholderText()
        return None

    def _hidden(self, feature_index):
        f = self.editor.by_index.get(feature_index)
        return f is not None and feature_key(f) in self.editor.hidden_features

    def repaint_rows(self):
        """A feature was hidden or shown: only the colour of the rows changes."""
        if self.shown:
            self.dataChanged.emit(self.index(0, 0), self.index(len(self.shown) - 1, len(self.columns) - 1), [Qt.ForegroundRole])


def master_key(editor, f):
    """Group key if `f` is the gene feature of a group, else None."""
    key = genes.group_key(f)
    return key if f.kind == genes.GENE and key in editor.groups and editor.groups[key].master is f else None


class GroupTableView(QTableView):
    """A click on the ▸ of a gene row opens/closes it."""

    def mousePressEvent(self, event):
        index = self.indexAt(event.position().toPoint())
        model = self.model()
        if index.isValid() and index.column() == 0 and model.shown[index.row()][4] and \
                event.position().x() - self.visualRect(index).left() < TOGGLE_PX + 16:  # 16 = color icon
            model.toggle(index.row())
            return
        super().mousePressEvent(event)


class FeatureTable(QWidget):
    def __init__(self, editor):
        super().__init__()
        self.editor = editor
        self.model = FeatureModel(editor)

        self.filter = QLineEdit(placeholderText="Filter features…", clearButtonEnabled=True)
        self.filter.textChanged.connect(self.model.set_filter)
        self.view = GroupTableView()
        self.view.setModel(self.model)
        self.view.horizontalHeader().setSortIndicator(self.model.sort_column, Qt.AscendingOrder)
        self.view.setSortingEnabled(True)
        self.view.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.view.setSelectionMode(QAbstractItemView.ExtendedSelection)  # several, to group them
        self.view.setEditTriggers(QAbstractItemView.NoEditTriggers)
        # columns as wide as what they hold, the name taking the room left (never a Strand column spread over the whole tab);
        # rows of one fixed height from the font (no per-row measuring)
        fit_columns(self.view, stretch=0)
        self.view.setContextMenuPolicy(Qt.CustomContextMenu)
        self.view.customContextMenuRequested.connect(self._menu)
        # (the database's Source holds the record too: it is filtered by the source's name)
        value_filter.install(self.view, [c for c in FILTERED if c in self.model.columns], keys={"Source": source_of})

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self.filter)
        layout.addWidget(self.view)

        self.view.selectionModel().currentRowChanged.connect(self._row_changed)
        self.view.doubleClicked.connect(self._double_clicked)
        # The Delete key deletes the selected features: read from the table's own key presses (a window-wide shortcut fired
        # only in the active window, and the key belongs to the table only while it has the keyboard anyway)
        self.view.installEventFilter(self)
        editor.changed.connect(self.model.reload)
        editor.display_changed.connect(self.model.repaint_rows)
        editor.selection_changed.connect(self._follow_selection)
        self._syncing = False

    def _double_clicked(self, index):
        if not getattr(self.editor, "database", None):   # the database's features are not edited here
            self.editor.edit_feature_requested.emit(self.model.shown[index.row()][3])

    def selected_indices(self):
        rows = sorted({i.row() for i in self.view.selectionModel().selectedRows()})
        return [self.model.shown[r][3] for r in rows if r < len(self.model.shown)]

    def _row_changed(self, current, _previous):
        if self._syncing or not current.isValid():
            return
        self.editor.select_feature(self.model.shown[current.row()][3])

    def eventFilter(self, obj, event):
        if obj is self.view and event.type() == QEvent.KeyPress and event.matches(QKeySequence.Delete):
            self._delete()
            return True
        return super().eventFilter(obj, event)

    def _delete(self):
        chosen = self.selected_indices()
        if getattr(self.editor, "database", None):
            self.editor.database.delete(self, chosen)
            return
        if len(chosen) > 1:
            delete_features(self, self.editor, chosen)
            return
        row = self.view.currentIndex().row()
        if 0 <= row < len(self.model.shown):
            delete_feature(self, self.editor, self.model.shown[row][3])

    def _menu(self, pos):
        ed, chosen = self.editor, self.selected_indices()
        menu = QMenu(self)
        if getattr(ed, "database", None):   # the Features Database: its own actions, not the editing ones of a file
            ed.database.menu(menu, self, chosen)
            if not menu.isEmpty():
                menu.exec(self.view.viewport().mapToGlobal(pos))
            return
        if len(chosen) == 1:
            feature_actions(menu, self, ed, chosen[0])
            menu.addSeparator()
        elif len(chosen) > 1:   # the same actions, applied to every selected feature
            features_actions(menu, self, ed, chosen)
            menu.addSeparator()
        menu.addAction("Group as Gene…", lambda: group_as_gene(self, ed, chosen)).setEnabled(bool(chosen))
        keys = {ed.part_of.get(i) or genes.group_key(ed.by_index[i]) for i in chosen if i in ed.by_index}
        keys = [k for k in keys if k in ed.groups]
        if len(keys) == 1:
            menu.addAction("Ungroup Gene", lambda: genes.ungroup(ed, keys[0]))
            menu.addAction("Delete Gene and Its Parts", lambda: genes.delete_gene(ed, keys[0]))
        menu.addSeparator()
        menu.addAction("Add to Features Database", lambda: self.add_to_library(chosen)).setEnabled(bool(chosen))
        menu.exec(self.view.viewport().mapToGlobal(pos))

    def add_to_library(self, indices):
        """Selected features go in the Features Database as Local features (found by Detect Features)."""
        ed = self.editor
        feats = [ed.by_index[i] for i in indices if i in ed.by_index and len(ed.by_index[i].spans) == 1]
        added = add_parts([feature_part(ed.doc, f, ed.path.name) for f in feats if f.kind != genes.GENE])
        ed.error.emit(f"{plural(added, 'feature')} added to the Features Database (Local).")  # status bar message
        return added

    def _follow_selection(self):
        """Highlight the row of a feature selected from another view (map/sequence)."""
        if self.editor.current_feature is None:
            return
        row = self.model.row_of(self.editor.current_feature)
        if row is None or row == self.view.currentIndex().row():
            return
        self._syncing = True
        index = self.model.index(row, 0)
        self.view.setCurrentIndex(index)
        self.view.scrollTo(index)
        self._syncing = False


def feature_actions(menu, parent, editor, index):
    """Edit / Delete (and, for a CDS, its translation on/off) for feature `index`."""
    f = editor.by_index.get(index)
    if f is None:
        return
    menu.addSection(f.label or f.kind)
    menu.addAction("Edit…", lambda: editor.edit_feature_requested.emit(index))
    menu.addAction("Delete", lambda: delete_feature(parent, editor, index))
    hidden = feature_key(f) in editor.hidden_features
    menu.addAction("Show on Sequence" if hidden else "Hide", lambda: editor.set_features_hidden([feature_key(f)], not hidden)).setToolTip(
        "Stop drawing this feature in the views (it stays in the file); Show in the feature table or the display filter brings it back")
    if f.kind == "CDS":
        a = menu.addAction("Show Translation", lambda: editor.toggle_translation(index))
        a.setCheckable(True)
        a.setChecked(editor.translated(f))
        a.setEnabled(editor.show_cds_translation)  # all off in the display panel


def features_actions(menu, parent, editor, indices):
    """Delete / Hide / Show (and, for CDSs, their translation on/off) for several features at once."""
    feats = [editor.by_index[i] for i in indices if i in editor.by_index]
    if not feats:
        return
    menu.addSection(f"{len(feats)} features selected")
    menu.addAction("Edit…", lambda: editor.edit_feature_requested.emit(feats[0].index)).setEnabled(False)
    menu.addAction(f"Delete {len(feats)} Features", lambda: delete_features(parent, editor, [f.index for f in feats]))
    keys = [feature_key(f) for f in feats]
    hidden = [k in editor.hidden_features for k in keys]
    if not all(hidden):
        menu.addAction("Hide All", lambda: editor.set_features_hidden(keys, True)).setToolTip(
            "Stop drawing these features in the views (they stay in the file)")
    if any(hidden):
        menu.addAction("Show All on Sequence", lambda: editor.set_features_hidden(keys, False))
    cds = [f for f in feats if f.kind == "CDS"]
    if cds:
        on = [f for f in cds if editor.translated(f)]
        a = menu.addAction(f"Show Translation ({len(cds)} CDS)", lambda: set_translations(editor, cds, len(on) < len(cds)))
        a.setCheckable(True)
        a.setChecked(len(on) == len(cds))
        a.setEnabled(editor.show_cds_translation)  # all off in the display panel


def set_translations(editor, cds, on):
    """Draw (on) or not the amino acids of these CDSs."""
    changed = False
    for f in cds:
        if editor.translated(f) != on:
            editor.untranslated ^= {(f.location, f.label)}
            changed = True
    if changed:
        editor.display_changed.emit()


def delete_features(parent, editor, indices):
    """Delete several features as one undo step; a gene among them asks once whether its parts go too."""
    chosen = {i for i in indices if i in editor.by_index}
    owned = set()   # parts of the selected genes
    for i in chosen:
        f = editor.by_index[i]
        key = genes.group_key(f) if f.kind == genes.GENE else None
        if key in editor.groups and editor.groups[key].master is f:
            owned |= {p.index for p in editor.groups[key].parts} - chosen
    if owned:
        answer = ask(parent, "Delete Features", f"The selection holds genes with {plural(len(owned), 'other part')}.",
                     [("all", "Delete the Parts Too"), ("selected", "Delete the Selection Only")], default="selected",
                     detail="Deleting the selection only leaves the other parts, ungrouped.")
        if answer is None:
            return
        if answer == "all":
            chosen |= owned
    if chosen:
        editor.edit_features([(i, None) for i in sorted(chosen)])


def feature_menu(parent, editor, index, global_pos):
    """Right-click on a feature in any view."""
    editor.select_feature(index)
    menu = QMenu(parent)
    feature_actions(menu, parent, editor, index)
    menu.exec(global_pos)


def primer_actions(menu, editor, name):
    """Edit / Delete for the primer called `name`: it is a primer of the database, so deleting it takes it out of every sequence
    (a step in this document's history: Undo brings it back). Import to This Sequence makes a primer that is only in the database
    (or a temporary one) part of this file. A primer that is not in the database (a temporary one) cannot be edited or deleted."""
    menu.addSection(name)
    known = db().find(name) is not None
    edit = menu.addAction("Edit…", lambda: editor.edit_primer_requested.emit(name))
    delete = menu.addAction("Delete from Primer Database", lambda: delete_primer(editor, name))
    for a in (edit, delete):
        a.setEnabled(known)
        if not known:
            a.setToolTip("Not in the primer database yet (a temporary primer): save it first")
    menu.addSeparator()
    todo = primers_to_import(editor, [name])
    imp = menu.addAction("Import to This Sequence", lambda: editor.import_primers(todo))
    imp.setEnabled(bool(todo))
    imp.setToolTip("Make this primer part of this file: drawn as a regular primer annotation and saved with the file" if todo
                   else "Already part of this file" if any(s.label == name for s in editor.primer_sites)
                   else "This primer does not bind here")
    menu.addSeparator()
    hidden = name in editor.hidden_primers
    toggle = menu.addAction("Show in This File" if hidden else "Hide in This File", lambda: editor.set_primer_hidden(name, not hidden))
    toggle.setToolTip("Stop drawing this primer on this file; it stays in the primer database (the Primers view shows it again)")


def primers_actions(menu, editor, names):
    """primer_actions for several primers at once: delete, import, hide/show them all as one step."""
    names = sorted(names)
    menu.addSection(f"{len(names)} primers selected")
    menu.addAction("Edit…").setEnabled(False)
    known = [n for n in names if db().find(n) is not None]
    delete = menu.addAction(f"Delete {len(known)} from Primer Database", lambda: delete_primers(editor, known))
    delete.setEnabled(bool(known))
    if len(known) < len(names):
        delete.setToolTip("Temporary primers (not in the primer database) are left out")
    menu.addSeparator()
    todo = primers_to_import(editor, names)
    imp = menu.addAction(f"Import {len(todo)} to This Sequence" if todo else "Import to This Sequence",
                         lambda: editor.import_primers(todo))
    imp.setEnabled(bool(todo))
    imp.setToolTip("Make these primers part of this file" if todo else "Already part of this file, or they do not bind here")
    menu.addSeparator()
    hidden = [n in editor.hidden_primers for n in names]
    if not all(hidden):
        menu.addAction("Hide All in This File", lambda: editor.set_primers_hidden(names, True)).setToolTip(
            "Stop drawing these primers on this file; they stay in the primer database")
    if any(hidden):
        menu.addAction("Show All in This File", lambda: editor.set_primers_hidden(names, False))


def delete_primers(editor, names):
    if not names:
        return
    remove_primers(editor, names)
    editor.error.emit(f"{plural(len(names), 'primer')} deleted from the primer database. Edit → Undo brings them back.")   # status bar


def primers_to_import(editor, names):
    """[(name, sequence)] of the primers called `names` that bind here and are not part of this file yet."""
    out = {}
    for site in editor.primer_sites:
        if site.label in names and not editor.in_file(site.label, site.sequence):
            out[(site.label, site.sequence.upper())] = (site.label, site.sequence)
    return list(out.values())


def delete_primer(editor, name):
    remove_primers(editor, [name])
    editor.error.emit(f"Primer {name} deleted from the primer database. Edit → Undo brings it back.")   # status bar


def primer_menu(parent, editor, key, global_pos):
    """Right-click on a primer site in any view."""
    site = editor.site_by_key.get(key)
    if site is None:
        return
    editor.select_primer(key)
    menu = QMenu(parent)
    primer_actions(menu, editor, site.label)
    menu.exec(global_pos)


def selection_menu(parent, editor, global_pos):
    """Right-click on selected bases: primers, feature, gene from them; copy."""
    from feature_dialog import edit_feature  # here: it imports the tools
    from primer_design import design_for_selection
    from primer_panel import PrimerDialog
    window = parent.window()
    menu = QMenu(parent)
    menu.addSection(f"{editor.end - editor.start:,} bp selected")
    menu.addAction("Add Primer (F)", lambda: PrimerDialog(window, editor).exec())
    menu.addAction("Add Primer (R)", lambda: reverse_primer_dialog(window, editor).exec())
    menu.addAction("Design Primers for This Sequence…", lambda: design_for_selection(window, editor))
    menu.addSeparator()
    menu.addAction("Add Feature…", lambda: edit_feature(window, editor))
    menu.addAction("Add Gene…", lambda: group_as_gene(window, editor, features_inside(editor)))
    menu.addSeparator()
    from sequence_view import add_edit_actions
    add_edit_actions(menu, editor, getattr(parent, "editable", True))   # Cut / Copy / Paste / Select All, as on empty space
    menu.addAction("Copy Reverse Complement\t" + QKeySequence("Ctrl+Shift+C").toString(QKeySequence.NativeText),
                   lambda: QGuiApplication.clipboard().setMimeData(editor.copy_reverse_complement_mime()))
    menu.exec(global_pos)


def reverse_primer_dialog(parent, editor):
    """New primer that anneals to the top strand of the selection: its reverse complement."""
    from primer_panel import PrimerDialog
    dialog = PrimerDialog(parent, editor)
    dialog.sequence.setText(reverse_complement(editor.selected_sequence()))
    return dialog


def features_inside(editor):
    """Indices of the features lying inside the selection (genes excluded)."""
    inside = [(editor.start, editor.end)] if editor.end <= editor.length else editor.selected_spans()
    return [f.index for f in editor.features
            if f.kind != genes.GENE and any(a <= f.spans[0][0] and f.spans[-1][1] <= b for a, b in inside)]


def group_as_gene(parent, editor, indices):
    """Ask a name (default: the largest CDS's; a warning when there is no CDS) and group
    `indices` under a new gene."""
    feats = [editor.by_index[i] for i in indices if i in editor.by_index]
    if not feats:
        QMessageBox.information(parent, "Add gene", "No feature inside the selection to make a gene of.")
        return None
    cds = sorted((f for f in feats if f.kind == "CDS"), key=lambda f: -sum(e - s for s, e in f.spans))
    prompt = f"Gene made of {plural(len(feats), 'feature')}: " + ", ".join(f.label for f in feats[:6]) + (
        "…" if len(feats) > 6 else "") + "\n"
    if not cds:
        prompt += "\nNote: there is no CDS among them.\n"
    name, ok = QInputDialog.getText(parent, "Group as Gene", prompt + "\nGene name:", text=(cds or feats)[0].label)
    if not ok:
        return None
    if not name.strip():
        QMessageBox.warning(parent, "Group as gene", "Give the gene a name.")
        return None
    problem = validate.text_problem(name.strip(), "gene name", 200)
    if problem:
        QMessageBox.warning(parent, "Group as gene", problem)
        return None
    try:
        return genes.make_gene(editor, indices, name.strip())
    except (OSError, ValueError) as e:
        QMessageBox.warning(parent, "Group as gene", str(e))
        return None


def delete_feature(parent, editor, index):
    """Delete a feature; for a gene, ask whether its parts go too."""
    f = editor.by_index.get(index)
    key = genes.group_key(f) if f is not None and f.kind == genes.GENE else None
    if key in editor.groups and editor.groups[key].master is f and editor.groups[key].parts:
        g = editor.groups[key]
        answer = ask(parent, "Delete Gene", f"{f.label} is made of {plural(len(g.parts), 'part')}. Delete them too?",
                     [("all", "Delete Gene and Parts"), ("gene", "Delete Gene Only")], default="gene",
                     detail="Deleting the gene only leaves its parts, ungrouped.")
        if answer == "all":
            genes.delete_gene(editor, key)
        elif answer == "gene":
            genes.ungroup(editor, key)
        return
    editor.delete_feature(index)
