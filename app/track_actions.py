"""The Alignment view's actions on regions and on what tracks hold (the menus of genome_view and track_panel call these):

- a selection: Highlight Selection, Recap Tracks Here, Export This Region as Image; Recap Highlighted Regions;
- an annotation track: Add as Feature (one item), Copy to Features, Next / Previous Item, Export Items as BED;
- a signal track: Regions Above Threshold (as highlights or as features).

The work is done by modules without widgets (track_recap, track_features, track_layout): these are the dialogs around it.
"""
from pathlib import Path

from PySide6.QtCore import Qt
from PySide6.QtGui import QColor, QGuiApplication, QIcon, QPixmap
from PySide6.QtWidgets import (
    QAbstractItemView, QCheckBox, QColorDialog, QComboBox, QDialog, QDialogButtonBox, QDoubleSpinBox, QFileDialog, QFormLayout,
    QHBoxLayout, QHeaderView, QLabel, QLineEdit, QMessageBox, QPushButton, QSpinBox, QTableWidget, QTableWidgetItem, QTreeWidget,
    QTreeWidgetItem, QVBoxLayout, QWidget,
)

import track_features
import track_recap
from editor import split_at_origin
from track_layout import REGION_COLOR, Region, read_regions_bed, regions_bed
from uikit import plural, run_busy

_SWATCHES = {}   # colour -> its icon


def selection_spans(editor):
    """The selected bases as [(start, end)] (two parts when the selection runs through the origin), [] when none."""
    return list(editor.selected_spans()) if editor.end > editor.start else []


def swatch(color):
    """A small square of `color` (an icon: the tracks' and regions' colours in the panel), made once per colour."""
    if color not in _SWATCHES:
        pix = QPixmap(10, 10)
        pix.fill(QColor(color))
        _SWATCHES[color] = QIcon(pix)
    return _SWATCHES[color]


# ---- the menus ---------------------------------------------------------------------------------------------------------------
def region_actions(menu, view):
    """The actions on the selection and on the highlighted regions."""
    ed = view.editor
    spans = selection_spans(ed)
    menu.addSection("Selection" if spans else "Regions")
    menu.addAction("Highlight Selection…", lambda: highlight_selection(view)).setEnabled(bool(spans))
    menu.addAction("Recap Tracks Here…", lambda: show_recap(view, [("selection", s, e) for s, e in spans])).setEnabled(bool(spans))
    menu.addAction("Export This Region as Image…", lambda: export_region_image(view, *spans[0])).setEnabled(len(spans) == 1)
    regions = ed.track_layout.regions
    menu.addAction(f"Recap Highlighted Regions ({len(regions)})…",
                   lambda: show_recap(view, [(r.name, r.start, r.end) for r in regions])).setEnabled(bool(regions))


def track_actions(menu, view, keys):
    """The actions on what tracks `keys` hold: their items into features, a signal's high stretches, moving item to item."""
    ed = view.editor
    entries = [e for k in keys if (e := ed.track_of(k)) is not None and e.name == k]
    annotations = [e for e in entries if e.kind == "track"]
    signals = [e for e in entries if e.kind == "signal"]
    if not (annotations or len(signals) == 1):
        return
    menu.addSeparator()
    if len(annotations) == 1:
        t = annotations[0]
        menu.addAction("Next Item", lambda: go_to_item(view, t, True))
        menu.addAction("Previous Item", lambda: go_to_item(view, t, False))
        menu.addAction("Export Items as BED…", lambda: export_items_bed(view, t))
    if annotations:
        menu.addAction("Copy to Features…" if len(annotations) == 1 else f"Copy {len(annotations)} Tracks to Features…",
                       lambda: copy_to_features(view, annotations))
    if len(signals) == 1:
        menu.addAction("Regions Above Threshold…", lambda: threshold_regions(view, signals[0]))


def item_actions(menu, view, entry, item):
    """The actions on one item of an annotation track (right-click on it)."""
    s, e, _, name, *_ = item
    menu.addSection(f"{name or 'item'}  {s + 1:,}–{e:,}")
    menu.addAction("Add as Feature", lambda: add_item_feature(view, entry, item))
    parts = split_at_origin(s, e, view.editor.length)
    menu.addAction("Highlight This Item", lambda: view.editor.track_layout.add_regions([Region(a, b, name) for a, b in parts]))


# ---- regions -----------------------------------------------------------------------------------------------------------------
class RegionDialog(QDialog):
    """A highlight's name and colour."""

    def __init__(self, parent, title, name="", color=REGION_COLOR):
        super().__init__(parent)
        self.setWindowTitle(title)
        self.name = QLineEdit(name, placeholderText="e.g. promoter, deleted region")
        self.color = color
        self.color_button = QPushButton("Choose…")
        self.color_button.setIcon(swatch(color))
        self.color_button.clicked.connect(self._choose)
        form = QFormLayout(self)
        form.addRow("Name:", self.name)
        form.addRow("Colour:", self.color_button)
        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        form.addRow(buttons)

    def _choose(self):
        c = QColorDialog.getColor(QColor(self.color), self, "Colour of the highlight")
        if c.isValid():
            self.color = c.name()
            self.color_button.setIcon(swatch(self.color))


def highlight_selection(view):
    ed = view.editor
    spans = selection_spans(ed)
    if not spans:
        return
    dialog = RegionDialog(view, "Highlight region")
    if dialog.exec() == QDialog.Accepted:
        name = dialog.name.text().strip()
        ed.track_layout.add_regions([Region(s, e, name, dialog.color) for s, e in spans])


def jump_to(view, start, end):
    """Select start..end and show it with a margin."""
    view.editor.select(start, end)
    view.editor.reveal.emit(start)
    view.show_selection()


def import_regions(view):
    ed = view.editor
    path, _ = QFileDialog.getOpenFileName(view, "Import highlighted regions", "", "BED (*.bed *.txt);;All files (*)")
    if not path:
        return
    try:
        regions, unread = read_regions_bed(Path(path).read_text(encoding="utf-8", errors="replace"), ed.length)
    except OSError as e:
        QMessageBox.critical(view, "Import regions", f"{path}\n\n{e}")
        return
    ed.track_layout.add_regions(regions)
    if unread:   # (named, never dropped silently)
        QMessageBox.warning(view, "Import regions", f"{plural(len(regions), 'region')} imported; {plural(len(unread), 'line')} not read:\n\n"
                            + "\n".join(unread[:20]) + ("\n…" if len(unread) > 20 else ""))


def export_regions(view, regions):
    path, _ = QFileDialog.getSaveFileName(view, "Export highlighted regions", "regions.bed", "BED (*.bed)")
    if path:
        try:
            Path(path).write_text(regions_bed(regions, view.editor.doc.name), encoding="utf-8")
        except OSError as e:
            QMessageBox.critical(view, "Export regions", f"{path}\n\n{e}")


def export_region_image(view, start, end):
    """The view's tracks over start..end, as they are shown (order, looks, highlights), as an SVG or PNG file."""
    from genome_view import write_region
    path, chosen = QFileDialog.getSaveFileName(view, "Export this region as image", f"region_{start + 1}-{end}.svg",
                                               "SVG, vector (*.svg);;PNG (*.png)")
    if not path:
        return
    try:
        write_region(view.editor, start, end, path)
    except OSError as e:
        QMessageBox.critical(view, "Export region", f"{path}\n\n{e}")
        return
    view.info.setText(f"Written: {path}")


# ---- recap -------------------------------------------------------------------------------------------------------------------
def show_recap(view, regions):
    if regions:
        RecapDialog(view, view.editor, regions).show()


class RecapDialog(QDialog):
    """Every track over one or more regions, side by side (track_recap): one row per track (and region), a row per group with
    the mean of its tracks; sortable, copied or exported as CSV."""

    def __init__(self, parent, editor, regions):
        super().__init__(parent)
        self.editor, self.regions = editor, regions
        self.setWindowTitle("Recap of the tracks")
        self.setAttribute(Qt.WA_DeleteOnClose)
        self.resize(1100, 600)
        self.summary = QLabel()
        self.per_base = QCheckBox(f"Value at every base (signals, regions up to {track_recap.PER_BASE:,} bp)")
        self.per_base.toggled.connect(self.fill)
        self.table = QTableWidget(0, len(track_recap.COLUMNS))
        self.table.setHorizontalHeaderLabels([c.capitalize() for c in track_recap.COLUMNS])
        self.table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.table.horizontalHeader().setSectionResizeMode(len(track_recap.COLUMNS) - 1, QHeaderView.Stretch)
        copy, export, close = QPushButton("Copy"), QPushButton("Export CSV…"), QPushButton("Close")
        copy.clicked.connect(lambda: QGuiApplication.clipboard().setText(track_recap.tsv(self.rows)))
        export.clicked.connect(self.export)
        close.clicked.connect(self.close)
        buttons = QHBoxLayout()
        buttons.addWidget(self.per_base)
        buttons.addStretch()
        for b in (copy, export, close):
            buttons.addWidget(b)
        layout = QVBoxLayout(self)
        layout.addWidget(self.summary)
        layout.addWidget(self.table, 1)
        layout.addLayout(buttons)
        self.rows = []
        self.fill()

    def fill(self):
        per_base = self.per_base.isChecked()
        try:
            recaps = run_busy(self.parentWidget(), "Reading the tracks…", lambda: track_recap.recap(self.editor, self.regions, per_base))
        except OSError as e:
            QMessageBox.critical(self, "Recap", str(e))
            return
        self.rows = track_recap.rows(recaps)
        bp = sum(r["length"] for r in recaps)
        self.summary.setText(f"{plural(len(recaps), 'region')} ({bp:,} bp), {plural(len(recaps[0]['tracks']) if recaps else 0, 'track')}. "
                             "Click a column's title to sort; a group's row holds the mean of its tracks.")
        self.table.setSortingEnabled(False)
        self.table.setRowCount(len(self.rows))
        for i, row in enumerate(self.rows):
            for j, column in enumerate(track_recap.COLUMNS):
                value = row[column]
                item = QTableWidgetItem()
                if isinstance(value, (int, float)):
                    item.setData(Qt.DisplayRole, value)   # (sorted as numbers)
                else:
                    item.setText("" if value is None else str(value))
                item.setToolTip(item.text() if column == "details" else "")
                if row["kind"] == "group":
                    font = item.font()
                    font.setBold(True)
                    item.setFont(font)
                self.table.setItem(i, j, item)
        self.table.setSortingEnabled(True)
        self.table.resizeColumnsToContents()

    def export(self):
        path, _ = QFileDialog.getSaveFileName(self, "Export the recap", "recap.csv", "CSV (*.csv)")
        if path:
            try:
                Path(path).write_text(track_recap.tsv(self.rows, ","), encoding="utf-8")
            except OSError as e:
                QMessageBox.critical(self, "Export the recap", f"{path}\n\n{e}")


# ---- tracks into features ----------------------------------------------------------------------------------------------------
def add_item_feature(view, entry, item):
    ed = view.editor
    label = ed.track_layout.label(entry.name)
    feats, skipped = track_features.to_features([item], label, track_features.default_type(label), "label",
                                                track_features.feature_keys(ed.features), ed.length)
    if skipped:
        view.info.setText("That feature is already in the sequence (same type, place and strand): not added again")
    else:
        ed.add_features(feats)


SCOPES = ("Whole sequence", "Shown in the view", "Selection", "Highlighted regions")


class CopyDialog(QDialog):
    """Copy to Features: which items (where), as what type, the name in which qualifier."""

    def __init__(self, parent, labels, has_selection, has_regions):
        super().__init__(parent)
        self.setWindowTitle("Copy to features")
        self.scope = QComboBox()
        self.scope.addItems(SCOPES)
        model = self.scope.model()
        model.item(2).setEnabled(has_selection)
        model.item(3).setEnabled(has_regions)
        self.kind = QComboBox(editable=True)
        self.kind.addItems(["misc_feature", "CDS", "gene", "promoter", "terminator", "regulatory", "ncRNA", "repeat_region"])
        self.kind.setCurrentText(track_features.default_type(labels[0]))
        self.name_key = QComboBox()
        self.name_key.addItems(track_features.NAME_FIELDS)
        form = QFormLayout(self)
        form.addRow(QLabel(f"Items of {', '.join(labels[:3])}{' …' if len(labels) > 3 else ''} become features of this sequence, "
                           "each with a note \"from track …\" and the file's extra columns as notes. An item already there as a "
                           "feature (same type, place and strand) is skipped. One undo step.", wordWrap=True))
        form.addRow("Items:", self.scope)
        form.addRow("Feature type:", self.kind)
        form.addRow("Their name goes in:", self.name_key)
        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        form.addRow(buttons)


IN_VIEW = SCOPES.index("Shown in the view")


def scope_spans(view, scope):
    """[(start, end)] of a scope, SCOPES[scope]."""
    ed = view.editor
    return (lambda: [(0, ed.length)], lambda: [(int(view.vstart), int(view.vend + 0.999))], lambda: selection_spans(ed),
            lambda: [(r.start, r.end) for r in ed.track_layout.regions])[scope]()


def copy_to_features(view, entries):
    ed = view.editor
    labels = [ed.track_layout.label(e.name) for e in entries]
    dialog = CopyDialog(view, labels, bool(selection_spans(ed)), bool(ed.track_layout.regions))
    if dialog.exec() != QDialog.Accepted:
        return
    spans = scope_spans(view, dialog.scope.currentIndex())
    kind, name_key = dialog.kind.currentText().strip() or "misc_feature", dialog.name_key.currentText()
    added, skipped = [], 0
    have = track_features.feature_keys(ed.features)
    for entry, label in zip(entries, labels):
        feats, n = track_features.to_features(track_features.items_in(entry, spans), label, kind, name_key, have, ed.length)
        added += feats
        skipped += n
    if added and QMessageBox.question(view, "Copy to features", f"Add {plural(len(added), 'feature')} ({kind}) to this sequence?"
                                      + (f" {plural(skipped, 'item')} already there skipped." if skipped else ""),
                                      QMessageBox.Yes | QMessageBox.Cancel, QMessageBox.Yes) == QMessageBox.Yes:
        ed.add_features(added)
        view.info.setText(f"{plural(len(added), 'feature')} added" + (f", {skipped} already there skipped" if skipped else ""))
    elif not added:
        QMessageBox.information(view, "Copy to features", "No item to add there" + (f" ({plural(skipped, 'item')} already features)." if skipped else "."))


# ---- moving item to item, BED ------------------------------------------------------------------------------------------------
def go_to_item(view, entry, forward):
    ed = view.editor
    pos = ed.start   # (the next item by start: one that starts inside the selected one is not skipped)
    item = track_features.next_item(entry, pos, forward)
    if item is None:
        view.info.setText(f"No {'next' if forward else 'previous'} item in {ed.track_layout.label(entry.name)}")
        return
    jump_to(view, item[0], min(item[1], ed.length))
    view.info.setText(f"{item[3] or 'item'}  {item[0] + 1:,}–{item[1]:,}")


def export_items_bed(view, entry):
    choice = QMessageBox(QMessageBox.Question, "Export items as BED", "Which items of the track?", parent=view)
    in_view = choice.addButton("Those in the View", QMessageBox.AcceptRole)
    every = choice.addButton("All", QMessageBox.AcceptRole)
    choice.addButton(QMessageBox.Cancel)
    choice.exec()
    if choice.clickedButton() not in (in_view, every):
        return
    spans = [(0, 1 << 62)] if choice.clickedButton() is every else scope_spans(view, IN_VIEW)
    path, _ = QFileDialog.getSaveFileName(view, "Export items as BED", Path(entry.name).stem + ".bed", "BED (*.bed)")
    if path:
        try:
            Path(path).write_text(track_features.items_bed(track_features.items_in(entry, spans), entry.seqid), encoding="utf-8")
        except OSError as e:
            QMessageBox.critical(view, "Export items as BED", f"{path}\n\n{e}")


# ---- a signal's high stretches -----------------------------------------------------------------------------------------------
class ThresholdDialog(QDialog):
    def __init__(self, parent, label, low, high, has_selection):
        super().__init__(parent)
        self.setWindowTitle("Regions above threshold")
        self.threshold = QDoubleSpinBox(decimals=4)
        self.threshold.setRange(-1e12, 1e12)   # (before the value: a value outside the default range would be clamped)
        self.threshold.setValue(low + (high - low) / 2)
        self.min_length = QSpinBox(minimum=1, maximum=10_000_000, value=1, suffix=" bp")
        self.gap = QSpinBox(minimum=0, maximum=10_000_000, value=0, suffix=" bp")
        self.scope = QComboBox()
        self.scope.addItems(SCOPES[:3])
        self.scope.model().item(2).setEnabled(has_selection)
        self.into = QComboBox()
        self.into.addItems(["Highlighted regions", "Features (misc_feature)"])
        form = QFormLayout(self)
        form.addRow(QLabel(f"Where {label} is at least the threshold (its values here: {low:g} to {high:g}).", wordWrap=True))
        form.addRow("Threshold:", self.threshold)
        form.addRow("Join stretches closer than:", self.gap)
        form.addRow("Keep stretches of at least:", self.min_length)
        form.addRow("Where:", self.scope)
        form.addRow("Make them:", self.into)
        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        form.addRow(buttons)


def threshold_regions(view, entry):
    ed = view.editor
    label = ed.track_layout.label(entry.name)
    try:
        _, _, low, high = run_busy(view, "Reading the values…", lambda: entry.signal.stats(entry.seqid, 0, ed.length))
    except OSError as e:
        QMessageBox.critical(view, "Regions above threshold", str(e))
        return
    if low != low:
        QMessageBox.information(view, "Regions above threshold", f"{label} has no value on this sequence.")
        return
    dialog = ThresholdDialog(view, label, low, high, bool(selection_spans(ed)))
    if dialog.exec() != QDialog.Accepted:
        return
    threshold, gap, least = dialog.threshold.value(), dialog.gap.value(), dialog.min_length.value()
    spans = scope_spans(view, dialog.scope.currentIndex())
    found = run_busy(view, "Finding the regions…", lambda: [r for s, e in spans for r in
                                                           track_features.signal_regions(entry.signal, entry.seqid, s, e, threshold, least, gap)])
    if not found:
        QMessageBox.information(view, "Regions above threshold", f"{label} is never at least {threshold:g} over {least} bp there.")
        return
    name = f"{label} ≥ {threshold:g}"
    if dialog.into.currentIndex() == 0:
        ed.track_layout.add_regions([Region(s, e, name) for s, e, _, _ in found])
    else:
        ed.add_features([("misc_feature", f"{s + 1}..{e}", [("label", name), ("note", f"from track {label}: mean {mean:.4g}, highest {top:.4g}")])
                         for s, e, mean, top in found])
    view.info.setText(f"{plural(len(found), 'region')} where {label} ≥ {threshold:g}")


# ---- the Regions tab of the Tracks panel -------------------------------------------------------------------------------------
class RegionsTab(QWidget):
    """The highlighted regions of the file: double-click one to go there, F2 renames it; recap, recolour, delete, BED."""

    def __init__(self, view):
        super().__init__()
        self.view, self.editor = view, view.editor
        self.tree = QTreeWidget()
        self.tree.setHeaderLabels(["Name", "Start", "End", "Length"])
        self.tree.setRootIsDecorated(False)
        self.tree.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.tree.setEditTriggers(QAbstractItemView.EditKeyPressed)
        self.tree.itemDoubleClicked.connect(lambda item, _: jump_to(self.view, *self._region(item)))
        self.tree.itemChanged.connect(self._renamed)
        self.empty = QLabel("No highlighted region yet: select bases in the view, right-click, Highlight Selection.", wordWrap=True)
        self.empty.setEnabled(False)
        bar = QHBoxLayout()
        for text, act, tip in (("Recap", self._recap, "Every track over the selected regions (all when none is selected)"),
                               ("Colour…", self._recolor, "Colour of the selected regions"),
                               ("Delete", self._delete, "Remove the selected highlights"),
                               ("Import…", lambda: import_regions(self.view), "Highlight the regions of a BED file"),
                               ("Export…", lambda: export_regions(self.view, self._chosen() or self.editor.track_layout.regions),
                                "The selected regions (all when none is selected) as a BED file")):
            b = QPushButton(text, toolTip=tip)
            b.clicked.connect(act)
            bar.addWidget(b)
        bar.addStretch()
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 4, 0, 0)
        layout.addLayout(bar)
        layout.addWidget(self.empty)
        layout.addWidget(self.tree, 1)
        self.shown = None   # the regions the list was built for (a change to tracks alone leaves it as it is)
        self.editor.tracks_changed.connect(self.refresh)
        self.refresh()

    def refresh(self):
        try:
            regions = self.editor.track_layout.regions
            shown = [(r.start, r.end, r.name, r.color) for r in regions]
            if shown == self.shown:
                return
            self.shown = shown
            self.tree.blockSignals(True)
            self.tree.clear()
            for i, r in enumerate(regions):
                item = QTreeWidgetItem([r.name, f"{r.start + 1:,}", f"{r.end:,}", f"{r.end - r.start:,}"])
                item.setData(0, Qt.UserRole, i)
                item.setIcon(0, swatch(r.color))
                item.setFlags(item.flags() | Qt.ItemIsEditable)
                self.tree.addTopLevelItem(item)
            self.tree.blockSignals(False)
            self.empty.setVisible(not regions)
        except RuntimeError:
            pass   # the pane is gone (closed) while the editor still signals

    def _region(self, item):
        r = self.editor.track_layout.regions[item.data(0, Qt.UserRole)]
        return r.start, r.end

    def _chosen(self):
        regions = self.editor.track_layout.regions
        return [regions[item.data(0, Qt.UserRole)] for item in self.tree.selectedItems()]

    def _renamed(self, item, column):
        if column == 0:
            region = self.editor.track_layout.regions[item.data(0, Qt.UserRole)]
            self.editor.track_layout.change_regions([region], name=item.text(0).strip())

    def _recap(self):
        regions = self._chosen() or self.editor.track_layout.regions
        show_recap(self.view, [(r.name, r.start, r.end) for r in regions])

    def _recolor(self):
        if chosen := self._chosen():
            c = QColorDialog.getColor(QColor(chosen[0].color), self, "Colour of the highlights")
            if c.isValid():
                self.editor.track_layout.change_regions(chosen, color=c.name())

    def _delete(self):
        if chosen := self._chosen():
            self.editor.track_layout.remove_regions(chosen)
