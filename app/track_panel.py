"""The Tracks panel of the Alignment view (genome_view): every track of the open file in its groups, to show, hide, reorder,
group and restyle many at once. It edits the file's track layout (track_layout, `Editor.track_layout`), and every view of
the file redraws on `Editor.tracks_changed`.

- Tick a box to show a track, untick to hide it; ticking one of several selected tracks ticks them all; a group's box ticks
  its tracks.
- Ctrl- and Shift-click select several; a group selected stands for all its tracks. The buttons and the right-click menu
  (`fill_menu`, also the menu of a track or group in the view) act on the selection.
- Drag tracks to reorder them, or onto a group to move them into it; groups are dragged too (never into another group).
- The filter keeps the tracks whose name has the text (and whole groups whose name has it); Select All picks those shown.
- F2 renames the selected track or group (the name shown: the file is never renamed).

A second tab lists the highlighted regions (track_actions.RegionsTab).
"""
from collections import Counter

from PySide6.QtCore import QItemSelection, QItemSelectionModel, Qt, QTimer
from PySide6.QtGui import QColor, QIcon
from PySide6.QtWidgets import (
    QAbstractItemView, QColorDialog, QDialog, QDialogButtonBox, QFormLayout, QHBoxLayout, QHeaderView, QInputDialog, QLabel,
    QLineEdit, QListWidget, QMenu, QMessageBox, QPushButton, QSpinBox, QTabWidget, QTreeWidget, QTreeWidgetItem, QVBoxLayout,
    QWidget,
)

import track_actions
import track_layout
from track_layout import FIXED, MODES, RENDERS, SCALED, Group
from uikit import plural

ROW, GROUP = "row", "group"
KEY = Qt.UserRole   # (ROW, row key) or (GROUP, group name) of each item of the tree
ROW_FLAGS = Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsDragEnabled | Qt.ItemIsUserCheckable | Qt.ItemIsEditable
GROUP_FLAGS = ROW_FLAGS | Qt.ItemIsDropEnabled   # (its box is set from its tracks here, not by Qt: see _set_group)


def show_menu(menu, global_pos):
    """Run a context menu and free it (its actions hold the tracks they act on)."""
    menu.exec(global_pos)
    menu.deleteLater()


def colored(view, keys):
    """The rows among `keys` drawn in one colour that can be chosen (not features, primers, reads, bases)."""
    return [k for k in keys if view.canvas.row_color(k)]


def confirm_remove(parent, editor, entries):
    """Take loaded tracks out of the window after asking (the files are not touched; it cannot be undone)."""
    reads = sum(e.kind == "bam" for e in entries)
    tracks = len(entries) - reads
    what = entries[0].name if len(entries) == 1 else " and ".join(
        t for t in (plural(tracks, "track") if tracks else "", plural(reads, "read alignment") if reads else "") if t)
    if QMessageBox.question(parent, "Remove Tracks", f"Remove {what} from this window? The files themselves are not changed.",
                            QMessageBox.Yes | QMessageBox.Cancel, QMessageBox.Cancel) == QMessageBox.Yes:
        editor.remove_tracks([e.name for e in entries])


# ---- the menu: the panel's More button, its right-click, and the right-click on a track or group in the view ------------------
def fill_menu(menu, view, keys, group=None):
    """Actions on tracks `keys` (rows of the file's layout) of GenomeView `view`; `group`: they are this group (its own actions)."""
    ed = view.editor
    layout = ed.track_layout
    kinds = {k: ed.row_table[k][0] for k in keys}
    n = len(keys)
    menu.addSection(group.name if group else layout.label(keys[0]) if n == 1 else f"{n} tracks")
    if group is not None:
        menu.addAction("Expand in the View" if group.collapsed else "Collapse in the View",
                       lambda: layout.collapse(group.name, not group.collapsed))
    if not keys:
        return
    menu.addAction("Show", lambda: layout.restyle(keys, shown=True))
    menu.addAction("Hide", lambda: layout.restyle(keys, shown=False))
    styled = [k for k in keys if kinds[k] not in FIXED]
    if styled:
        display = menu.addMenu("Display")
        modes = {layout.style(k).mode for k in styled}
        for mode in MODES:
            a = display.addAction(mode.capitalize(), lambda m=mode: layout.restyle(styled, mode=m))
            a.setCheckable(True)
            a.setChecked(modes == {mode})
        display.addSeparator()
        display.addAction("Height…", lambda: _height(view, styled))
        display.addAction("Reset Height", lambda: layout.restyle(styled, height=None)).setEnabled(
            any(layout.style(k).height for k in styled))
        scaled = [k for k in styled if kinds[k] in SCALED]
        if scaled:
            scale = display.addMenu("Scale")
            scales = {layout.style(k).scale for k in scaled}
            for key in track_layout.SCALES:
                a = scale.addAction(key.capitalize(), lambda s=key: layout.restyle(scaled, scale=s))
                a.setCheckable(True)
                a.setChecked(scales == {key})
    track_actions.track_actions(menu, view, keys)
    signals = [k for k in keys if kinds[k] == "signal"]
    if signals:
        _signal_menu(menu.addMenu("Signal"), view, signals)
    if paint := colored(view, keys):
        color = menu.addMenu("Colour")
        color.addAction("Choose…", lambda: _choose_color(view, paint))
        if len(paint) > 1:
            color.addAction("Light to Dark over These Tracks…", lambda: _gradient(view, paint))
        color.addAction("Back to Their Own Colours", lambda: layout.restyle(paint, color=None)).setEnabled(
            any(layout.style(k).color for k in paint))
    menu.addSeparator()
    move = menu.addMenu("Move to Group")
    for g in layout.groups():
        if g is not group:
            move.addAction(g.name, lambda name=g.name: layout.move(keys, name))
    move.addSeparator()
    move.addAction("New Group…", lambda: _new_group(view, keys))
    move.addAction("No Group", lambda: layout.move(keys, None))
    if n > 1:
        menu.addAction("Group by Name…", lambda: _group_by_name(view, keys))
    if group is not None:
        menu.addAction("Rename Group…", lambda: _rename_group(view, group))
        menu.addAction("Ungroup", lambda: layout.ungroup(group.name))
    elif n == 1:
        menu.addAction("Rename…", lambda: _rename(view, keys[0]))
    entries = list({id(e): e for k in keys if (e := ed.track_of(k)) is not None}.values())
    if entries:
        menu.addSeparator()
        menu.addAction("Remove from This Window…" if len(entries) == 1 else f"Remove {len(entries)} Tracks from This Window…",
                       lambda: confirm_remove(view, ed, entries))


def _signal_menu(menu, view, keys):
    """How signal tracks `keys` are drawn: bars, a line or a heatmap; mean or highest value per pixel; their range; overlaid."""
    layout = view.editor.track_layout
    styles = [layout.style(k) for k in keys]

    def choice(title, field, options):
        sub = menu.addMenu(title)
        now = {getattr(st, field) for st in styles}
        for value, label in options:
            a = sub.addAction(label, lambda v=value: layout.restyle(keys, **{field: v}))
            a.setCheckable(True)
            a.setChecked(now == {value})
        return sub

    choice("Draw As", "render", [(r, r.capitalize()) for r in RENDERS])
    choice("Value per Pixel", "summary", [("mean", "Mean"), ("max", "Highest")])
    ranges = choice("Range", "range", [("auto", "Automatic (Values in View)"), ("group", "Shared with Its Group")])
    fixed = ranges.addAction("Fixed…", lambda: _fixed_range(view, keys))
    fixed.setCheckable(True)
    fixed.setChecked({st.range for st in styles} == {"fixed"})
    menu.addSeparator()
    if len(keys) > 1:
        menu.addAction(f"Overlay in One Row ({len(keys)})", lambda: layout.overlay(keys))
    together = {k for k, st in layout.styles.items() if st.overlay in keys or (st.overlay and k in keys)}
    menu.addAction("Separate the Overlaid Tracks", lambda: layout.separate(keys)).setEnabled(bool(together))


def _fixed_range(view, keys):
    st = view.editor.track_layout.style(keys[0])
    low, ok = QInputDialog.getDouble(view, "Fixed range", "Lowest value (bottom of the row):", st.low, -1e12, 1e12, 3)
    if not ok:
        return
    high, ok = QInputDialog.getDouble(view, "Fixed range", "Highest value (top of the row):", max(st.high, low + 1), low, 1e12, 3)
    if ok and high > low:
        view.editor.track_layout.restyle(keys, range="fixed", low=low, high=high)


def _height(view, keys):
    st = view.editor.track_layout.style(keys[0])
    px, ok = QInputDialog.getInt(view, "Track height", "Height (px):", st.height or 30, 6, 4000)
    if ok:
        view.editor.track_layout.restyle(keys, height=px)


def _choose_color(view, keys):
    color = QColorDialog.getColor(QColor(view.canvas.row_color(keys[0])), view, "Track colour")
    if color.isValid():
        view.editor.track_layout.restyle(keys, color=color.name())


def _gradient(view, keys):
    color = QColorDialog.getColor(QColor(view.canvas.row_color(keys[0])), view, "Colour of the tracks: from light to dark")
    if color.isValid():
        view.editor.track_layout.restyle_each(
            {k: {"color": c} for k, c in zip(keys, track_layout.gradient(color.name(), len(keys)))})


def _new_group(view, keys):
    name, ok = QInputDialog.getText(view, "New group", "Name of the group:")
    if ok and name.strip():
        view.editor.track_layout.move(keys, name.strip())


def _rename(view, key):
    layout = view.editor.track_layout
    text, ok = QInputDialog.getText(view, "Rename track", "Name shown (the file keeps its name):", text=layout.label(key))
    if ok:
        layout.restyle([key], label=None if text.strip() in ("", key) else text.strip())


def _rename_group(view, group):
    text, ok = QInputDialog.getText(view, "Rename group", "Name of the group:", text=group.name)
    old, new = group.name, text.strip()
    if not ok or not new or new == old:
        return
    if view.editor.track_layout.rename_group(old, new):
        view.tracks_panel.renamed(old, new)
    else:
        QMessageBox.warning(view, "Rename group", f"There is already a group named {new!r}.")


class GroupByName(QDialog):
    """Group by Name: the first parts of each name (split on _ - . and spaces) name its group; the list shows the groups it
    makes before anything moves."""

    def __init__(self, parent, labels):
        super().__init__(parent)
        self.setWindowTitle("Group by name")
        self.labels = labels
        self.parts = QSpinBox(minimum=1, maximum=6, value=1)
        self.parts.valueChanged.connect(self._preview)
        self.preview = QListWidget()
        form = QFormLayout()
        form.addRow("Parts of the name that make the group:", self.parts)
        layout = QVBoxLayout(self)
        layout.addLayout(form)
        layout.addWidget(QLabel("Groups it makes (a name alone in its group is not moved):"))
        layout.addWidget(self.preview)
        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)
        self._preview()

    def _preview(self):
        self.preview.clear()
        groups = track_layout.name_groups(self.labels, self.parts.value())
        self.preview.addItems([f"{g} ({plural(n, 'track')})" for g, n in Counter(groups.values()).items()])
        if len(groups) < len(self.labels):
            self.preview.addItem(f"Not grouped: {plural(len(self.labels) - len(groups), 'track')}")


def _group_by_name(view, keys):
    layout = view.editor.track_layout
    dialog = GroupByName(view, list(dict.fromkeys(layout.label(k) for k in keys)))
    if dialog.exec() == QDialog.Accepted:
        layout.group_by_name(keys, dialog.parts.value())


# ---- the panel -------------------------------------------------------------------------------------------------------------------
class TrackTree(QTreeWidget):
    """The tree of groups and tracks. Dropping keeps one level: tracks go into a group or between items, a group never goes
    into another group."""

    def __init__(self, panel):
        super().__init__()
        self.panel = panel
        self.setHeaderHidden(True)
        self.setColumnCount(2)
        self.header().setStretchLastSection(False)
        self.header().setSectionResizeMode(0, QHeaderView.Stretch)
        self.header().setSectionResizeMode(1, QHeaderView.ResizeToContents)
        self.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.setDragDropMode(QAbstractItemView.InternalMove)
        self.setDefaultDropAction(Qt.MoveAction)
        self.setEditTriggers(QAbstractItemView.EditKeyPressed)   # F2: rename (a click selects, a drag moves)
        self.setUniformRowHeights(True)
        self.setIndentation(14)
        self.setContextMenuPolicy(Qt.CustomContextMenu)

    def dropEvent(self, event):
        target = self.itemAt(event.position().toPoint())
        on = self.dropIndicatorPosition() == QAbstractItemView.OnItem
        parent = None if target is None else target if on else target.parent()
        groups_moved = any(i.data(0, KEY)[0] == GROUP for i in self.selectedItems())
        if parent is not None and (parent.data(0, KEY)[0] != GROUP or groups_moved):
            event.ignore()
            return
        super().dropEvent(event)
        QTimer.singleShot(0, self.panel.tree_moved)   # (Qt finishes the move after this returns)


class TracksPanel(QWidget):
    def __init__(self, view):
        super().__init__()
        self.view, self.editor = view, view.editor
        self.setMinimumWidth(170)
        self.filter = QLineEdit(placeholderText="Filter tracks", clearButtonEnabled=True)
        self.filter.textChanged.connect(self._apply_filter)
        select = QPushButton("Select All", toolTip="Select every track the filter shows")
        select.clicked.connect(self.select_shown)
        top = QHBoxLayout()
        top.addWidget(self.filter, 1)
        top.addWidget(select)
        self.count = QLabel()
        self.count.setEnabled(False)   # (greyed: a note, not a control)
        bar = QHBoxLayout()
        layout_ = self.editor.track_layout
        for text, act, tip in (("Show", lambda keys: layout_.restyle(keys, shown=True), "Show the selected tracks"),
                               ("Hide", lambda keys: layout_.restyle(keys, shown=False), "Hide the selected tracks"),
                               ("Heatmap", self._heatmap, "The selected signal tracks as one-line heatmaps"),
                               ("Colour…", self._color, "One colour for the selected tracks (More: light to dark)")):
            b = QPushButton(text, toolTip=tip)
            b.clicked.connect(lambda _=False, act=act: self._on_selection(act))
            bar.addWidget(b)
        more = QPushButton("More", toolTip="Display, height, groups, rename, remove: for the selected tracks")
        more.clicked.connect(lambda: self._menu_at(more.mapToGlobal(more.rect().bottomLeft())))
        bar.addWidget(more)
        bar.addStretch()
        self.tree = TrackTree(self)
        self.tree.itemChanged.connect(self._item_changed)
        self.tree.itemSelectionChanged.connect(self._count)
        self.tree.customContextMenuRequested.connect(lambda pos: self._menu_at(self.tree.viewport().mapToGlobal(pos), pos))
        self.tree.itemExpanded.connect(lambda it: self.folded.discard(it.text(0)))
        self.tree.itemCollapsed.connect(lambda it: self.folded.add(it.text(0)))
        tracks = QWidget()
        layout = QVBoxLayout(tracks)
        layout.setContentsMargins(0, 4, 0, 0)
        layout.addLayout(top)
        layout.addWidget(self.count)
        layout.addLayout(bar)
        layout.addWidget(self.tree, 1)
        self.tabs = QTabWidget(documentMode=True)
        self.tabs.addTab(tracks, "Tracks")
        self.regions = track_actions.RegionsTab(view)
        self.tabs.addTab(self.regions, "Regions")
        outer = QVBoxLayout(self)
        outer.setContentsMargins(4, 4, 0, 4)
        outer.addWidget(self.tabs)
        self.folded = set()      # groups folded in this tree (only here: folding a group in the view is Group.collapsed)
        self.shape = None        # the order and groups the tree was built for (a style change updates the rows in place)
        self.dirty = True        # the tree is built when shown (a hidden panel does no work)
        # Many changes in a row (a group's box ticks all its tracks, one change each) are shown once.
        self.timer = QTimer(self, singleShot=True, interval=0, timeout=self._refresh_now)
        self.editor.tracks_changed.connect(self.refresh)

    # ---- building ------------------------------------------------------------------------------------------------------------
    def refresh(self):
        try:
            if self.isVisible():
                self.timer.start()
            else:
                self.dirty = True
        except RuntimeError:
            pass   # the pane is gone (closed) while the editor still signals

    def showEvent(self, event):
        super().showEvent(event)
        if self.dirty:
            self._refresh_now()

    def _refresh_now(self):
        self.dirty = False
        layout = self.editor.track_layout
        shape = tuple((n.name, tuple(n.keys), n.collapsed) if isinstance(n, Group) else n for n in layout.nodes)
        if shape != self.shape:
            self.shape = shape
            self.rebuild()
            return
        self.tree.blockSignals(True)   # only the look of rows changed: set it on the items there
        for item in self._row_items():
            self._set_row(item, item.data(0, KEY)[1])
        self._set_groups()
        self.tree.blockSignals(False)
        self._apply_filter()   # (a track renamed may match the filter or no longer)

    def rebuild(self):
        layout = self.editor.track_layout
        kept = set(self.selected_keys())
        scroll = self.tree.verticalScrollBar().value()
        self.tree.blockSignals(True)   # (no change, expansion or selection signal for the items built)
        self.tree.clear()
        for node in layout.nodes:
            if isinstance(node, Group):
                item = QTreeWidgetItem([node.name, "collapsed" if node.collapsed else ""])
                item.setData(0, KEY, (GROUP, node.name))
                item.setFlags(GROUP_FLAGS)
                item.addChildren([self._row(k) for k in node.keys])
                self.tree.addTopLevelItem(item)
                item.setExpanded(node.name not in self.folded)
            else:
                self.tree.addTopLevelItem(self._row(node))
        for item in self._row_items():
            item.setSelected(item.data(0, KEY)[1] in kept)
        self._set_groups()
        self.tree.blockSignals(False)
        self.tree.verticalScrollBar().setValue(scroll)
        self._apply_filter()

    def _row(self, key):
        item = QTreeWidgetItem()
        item.setData(0, KEY, (ROW, key))
        item.setFlags(ROW_FLAGS)
        self._set_row(item, key)
        return item

    def _set_row(self, item, key):
        """A row's look in the tree: its name, box, display mode and colour swatch."""
        st = self.editor.track_layout.style(key)
        item.setText(0, self.editor.track_layout.label(key))
        kind = self.editor.row_table[key][0]
        item.setText(1, "hidden" if not st.shown else "" if kind in FIXED else f"in {self.editor.track_layout.label(st.overlay)}"
                     if st.overlay else st.render if kind == "signal" and st.render == "heatmap" else st.mode)
        item.setCheckState(0, Qt.Checked if st.shown else Qt.Unchecked)
        color = self.view.canvas.row_color(key)
        item.setIcon(0, track_actions.swatch(color) if color else QIcon())
        item.setToolTip(0, key if item.text(0) != key else "")

    def _set_groups(self):
        """Each group's box: ticked when all its tracks are shown, empty when none is, half otherwise."""
        for i in range(self.tree.topLevelItemCount()):
            top = self.tree.topLevelItem(i)
            if top.data(0, KEY)[0] == GROUP:
                states = {top.child(j).checkState(0) for j in range(top.childCount())}
                top.setCheckState(0, states.pop() if len(states) == 1 else Qt.PartiallyChecked if states else Qt.Unchecked)

    def renamed(self, old, new):
        """A group was renamed: it stays folded here if it was."""
        if old in self.folded:
            self.folded.discard(old)
            self.folded.add(new)

    def _row_items(self):
        for i in range(self.tree.topLevelItemCount()):
            top = self.tree.topLevelItem(i)
            if top.data(0, KEY)[0] == ROW:
                yield top
            else:
                yield from (top.child(j) for j in range(top.childCount()))

    # ---- reading the selection -----------------------------------------------------------------------------------------------
    def selected_keys(self):
        """The rows selected (a group selected: all its rows), in display order."""
        layout = self.editor.track_layout
        chosen = set()
        for item in self.tree.selectedItems():
            kind, name = item.data(0, KEY)
            if kind == ROW:
                chosen.add(name)
            elif (g := layout.group(name)) is not None:
                chosen.update(g.keys)
        return [k for k in layout.keys() if k in chosen]

    def _selected_group(self):
        items = self.tree.selectedItems()
        if len(items) == 1 and items[0].data(0, KEY)[0] == GROUP:
            return self.editor.track_layout.group(items[0].data(0, KEY)[1])
        return None

    def _on_selection(self, act):
        if keys := self.selected_keys():
            act(keys)
        else:
            self.view.info.setText("Select tracks in the list first (Ctrl- and Shift-click pick several)")

    def _heatmap(self, keys):
        if signals := [k for k in keys if self.editor.row_table[k][0] == "signal"]:
            self.editor.track_layout.restyle(signals, render="heatmap", mode="dense")   # (one line each: many in one screen)
        else:
            self.view.info.setText("Heatmap is for signal tracks (bigWig, bedGraph, wig)")

    def _color(self, keys):
        if paint := colored(self.view, keys):
            _choose_color(self.view, paint)
        else:
            self.view.info.setText("These tracks are drawn in colours of their own (features, primers, reads)")

    def _menu_at(self, global_pos, pos=None):
        if pos is not None and (item := self.tree.itemAt(pos)) is not None and not item.isSelected():
            self.tree.clearSelection()
            item.setSelected(True)
        menu = QMenu(self)
        keys, group = self.selected_keys(), self._selected_group()
        if keys or group is not None:
            fill_menu(menu, self.view, keys, group)
        else:
            menu.addAction("Select tracks first (Ctrl- and Shift-click pick several)").setEnabled(False)
        show_menu(menu, global_pos)

    # ---- changes made in the tree --------------------------------------------------------------------------------------------
    def _item_changed(self, item, column):
        if column != 0:
            return
        layout = self.editor.track_layout
        kind, name = item.data(0, KEY)
        if kind == GROUP:
            new = item.text(0).strip()
            if new != name:
                if layout.rename_group(name, new):
                    self.renamed(name, new)
                else:   # empty, or the name of another group: put back
                    self.tree.blockSignals(True)
                    item.setText(0, name)
                    self.tree.blockSignals(False)
            elif (g := layout.group(name)) is not None:   # its box: all its tracks shown, or all hidden
                layout.restyle(g.keys, shown=item.checkState(0) != Qt.Unchecked)
            return
        shown = item.checkState(0) == Qt.Checked
        label = item.text(0).strip()
        if label != layout.label(name):
            layout.restyle([name], label=None if label in ("", name) else label)
        elif shown != layout.style(name).shown:
            # one box of a selection: all of them (their boxes are set when the panel shows the change)
            keys = self.selected_keys() if item.isSelected() else [name]
            layout.restyle(keys, shown=shown)

    def tree_moved(self):
        """After a drag and drop: the tree's order is the layout's."""
        tree = []
        for i in range(self.tree.topLevelItemCount()):
            top = self.tree.topLevelItem(i)
            kind, name = top.data(0, KEY)
            if kind == GROUP:
                tree.append((name, [top.child(j).data(0, KEY)[1] for j in range(top.childCount())]))
            else:
                tree.append((None, [name]))
        self.editor.track_layout.set_tree(tree)

    # ---- filter --------------------------------------------------------------------------------------------------------------
    def _apply_filter(self):
        text = self.filter.text().strip().lower()

        def match(item):
            return not text or text in item.text(0).lower() or text in item.data(0, KEY)[1].lower()

        for i in range(self.tree.topLevelItemCount()):
            top = self.tree.topLevelItem(i)
            if top.data(0, KEY)[0] == ROW:
                top.setHidden(not match(top))
                continue
            whole = match(top)   # a group whose name has the text: all its tracks
            seen = False
            for j in range(top.childCount()):
                child = top.child(j)
                child.setHidden(not (whole or match(child)))
                seen |= not child.isHidden()
            top.setHidden(bool(text) and not (whole or seen))
            if text and seen:
                top.setExpanded(True)
        self._count()

    def select_shown(self):
        model, chosen = self.tree.model(), QItemSelection()
        for item in self._row_items():
            if not item.isHidden():
                index = self.tree.indexFromItem(item)
                chosen.select(index, model.index(index.row(), 1, index.parent()))
        self.tree.selectionModel().select(chosen, QItemSelectionModel.ClearAndSelect)   # (one selection signal, not one each)

    def _count(self):
        rows = list(self._row_items())
        shown = sum(1 for it in rows if not it.isHidden())
        selected = len(self.selected_keys())
        text = f"{shown} of {plural(len(rows), 'track')}" if self.filter.text().strip() else plural(len(rows), "track")
        self.count.setText(text + (f" · {selected} selected" if selected else ""))
