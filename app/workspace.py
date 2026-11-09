"""Workspace: the tabs of a window. Each tab is a tree of split panes; each pane shows one
view or tool, picked from its header (any tool in any tab, the same tool several times).
All panes of a window share its Editor, so they stay in sync without knowing each other.

Layouts are plain data ({"view": name} or {"split": "h"/"v", "sizes": [...], "children":
[...]}, per tab). Every document opens with DEFAULT below, so the window you get is the same
every time; only an Open DNA Project (.odp) carries a layout of its own, and View > Reset
Tabs brings DEFAULT back after you have rearranged things.

A pane builds its view the first time it is shown, and drops it when it shows something else
or its tab closes: tabs never visited, and panes hidden behind another tab, cost nothing.
"""
import re

import icons
import validate

from PySide6.QtCore import QEvent, Qt, QTimer, Signal
from PySide6.QtWidgets import (
    QHBoxLayout, QInputDialog, QMenu, QMessageBox, QSplitter, QTabBar, QTabWidget, QToolButton, QVBoxLayout,
    QWidget,
)

from annotation import DetectDialog
from cas9 import Cas9Dialog
from checkup import CheckupDialog
from display_panel import QuickFilters
from feature_table import FeatureTable
from genome_view import GenomeView
from map_view import MapView
from local_align import LocalAlignDialog
from primer_panel import PrimerPanel
from protein import ProteinView
import wip
from sequence_view import SequenceView
from view_help import ViewHelpDialog, open_options, options_tip, OPTIONS
from view_picker import ViewPicker
from tools import DigestDialog, OrfDialog, PcrDialog

# Pane choices: name -> factory(window) -> widget.
VIEWS = {
    "Sequence": lambda w: SequenceView(w.editor),
    "Map": lambda w: MapView(w.editor),
    "Alignment": lambda w: GenomeView(w.editor),
    "Features": lambda w: FeatureTable(w.editor),
    "Primers": lambda w: PrimerPanel(w.editor),
    "Protein": lambda w: ProteinView(w.editor),
    "PCR": lambda w: PcrDialog(w, w.editor).embed(w.open_pcr_product),
    "Spacer Finder": lambda w: Cas9Dialog(w, w.editor).embed(),
    "Local Align": lambda w: LocalAlignDialog(w, w.editor).embed(),
    "Detect ORFs": lambda w: OrfDialog(w, w.editor).embed(),
    "Digest": lambda w: DigestDialog(w, w.editor).embed(),
    "Detect Features": lambda w: DetectDialog(w, w.editor).embed(),
    "Plasmid Checkup": lambda w: CheckupDialog(w, w.editor).embed(),
}


# The picker lists the views by group, in the order written here (None = a thin dividing line); every view not named here
# is a tool (listed alphabetically).
SEPARATOR = None
GROUPS = [("Views", ["Alignment", "Map", "Protein", "Sequence"]),
          ("Annotations", ["Features", "Primers", SEPARATOR, "Detect Features", "Detect ORFs"])]
# Tools that are tabs, not panes: choosing one in a picker opens (or shows) its tab and leaves the pane as it was.
TAB_TOOLS = {"Assembly": lambda window: window.show_assembly()}
# Picked in a pane, these tools bring the views they are read beside (added to the tab when it lacks them).
COMPANIONS = {"Plasmid Checkup": ("Map", "Sequence")}
TOOLS = "Tools"


def grouped_views():
    """[(group, [view names, or SEPARATOR])]: Views, Annotations as written in GROUPS, then every other view as a Tool."""
    named = {n for _, names in GROUPS for n in names}
    return [(g, list(names)) for g, names in GROUPS] + [(TOOLS, sorted([n for n in VIEWS if n not in named] + list(TAB_TOOLS)))]


def split(direction, sizes, *children):
    return {"split": direction, "sizes": list(sizes), "children": list(children)}


def pane(view):
    return {"view": view}


# Every document opens with this, whatever was on screen before: the sequence across 70 % of the window, the map above
# the features in the 30 % on the right, in one tab ("+" adds more). Nothing is added to it (no alignment pane for a long
# sequence), and a pane builds its view only when it is shown. The one exception is an
# Open DNA Project (.odp), which carries the layout it was saved with.
DEFAULT = [
    # Sizes are weights, not pixels: Qt shares the window out in this ratio (70 / 30, then half and half).
    {"layout": split("h", [7000, 3000], pane("Sequence"), split("v", [5000, 5000], pane("Map"), pane("Features")))},
]
# Views renamed since layouts were saved (Variants is now a track of the Alignment view).
# Primer Design is now a window opened from a selection (right-click > Design Primers for This Sequence).
OLD_NAMES = {"Genome": "Alignment", "Quick Align": "Local Align", "Cas9 Guides": "Spacer Finder", "ORFs": "Detect ORFs",
             "Variants": "Alignment", "Primer Design": "Primers"}
OLD_TAB_NAMES = {"Edit", "Primers", "Tools"}  # names the old default tabs were saved with: now automatic
MAX_TAB_NAME = 60
MAX_NAMED = 3  # views named in an automatic tab name; the others are counted


class Pane(QWidget):
    """Header (view picker, split right / split below / close) + the view itself."""

    def __init__(self, workspace, view):
        super().__init__()
        view = OLD_NAMES.get(view, view)
        self.workspace, self.view_name, self.view = workspace, view if view in VIEWS else "Sequence", None
        self.picker = ViewPicker(grouped_views)
        self.picker.set_current(self.view_name)
        self.picker.chosen.connect(self.set_view)
        header = QHBoxLayout()
        header.setContentsMargins(2, 1, 2, 1)
        header.addWidget(self.picker)
        self.wip = wip.Badge(wip.VIEWS.get(self.view_name))   # an unfinished view: its "WIP" tag beside its name
        header.addWidget(self.wip)
        # Next to the picker: the options of this view (a cog) and its quick guide (?). Both follow the picker's choice.
        self.cog = QToolButton(autoRaise=True, focusPolicy=Qt.TabFocus)
        self.cog.clicked.connect(lambda: open_options(workspace.window, self.view_name))
        self.help = QToolButton(text="?", toolTip="Quick guide to this view", autoRaise=True, focusPolicy=Qt.TabFocus)
        self.help.clicked.connect(lambda: ViewHelpDialog(self, self.view_name).exec())
        header.addWidget(self.cog)
        header.addWidget(self.help)
        self._update_helpers()
        header.addStretch()
        self.header_buttons = {}
        for name, tip, slot in (("right", "Split: new pane on the right", lambda: workspace.split(self, Qt.Horizontal)),
                                ("below", "Split: new pane below", lambda: workspace.split(self, Qt.Vertical)),
                                ("close", "Close this pane", lambda: workspace.close_pane(self))):
            b = QToolButton(toolTip=tip, autoRaise=True)
            b.clicked.connect(slot)
            header.addWidget(b)
            self.header_buttons[name] = b
        self.close_button = self.header_buttons["close"]
        self._paint_icons()
        self.box = QVBoxLayout(self)
        self.box.setContentsMargins(0, 0, 0, 0)
        self.box.setSpacing(0)
        self.box.addLayout(header)

    def _paint_icons(self):
        """Drawn icons in the theme's text colour (icons.pane_icon), made again on a theme switch."""
        self.cog.setIcon(icons.pane_icon("options", self.palette()))
        for name, b in self.header_buttons.items():
            b.setIcon(icons.pane_icon(name, self.palette()))

    def changeEvent(self, event):
        if event.type() == QEvent.PaletteChange:
            self._paint_icons()
        super().changeEvent(event)

    def _update_helpers(self):
        """The cog is enabled only when the view has options; its tooltip says where they are."""
        self.cog.setEnabled(OPTIONS.get(self.view_name) is not None)
        self.cog.setToolTip(options_tip(self.view_name))

    def showEvent(self, event):
        self.build()
        super().showEvent(event)

    def build(self):
        if self.view is None:
            self.view = VIEWS[self.view_name](self.workspace.window)
            self.box.addWidget(self.view, 1)
            self.view.show()  # a child added to an already shown widget starts hidden
            self.workspace.rewire()
        return self.view

    def set_view(self, name):
        if name in TAB_TOOLS:   # a tab of its own: open it, this pane keeps what it shows
            TAB_TOOLS[name](self.workspace.window)
            self.picker.set_current(self.view_name)
            return
        if name == self.view_name and self.view is not None:
            return
        if self.view is not None:
            self.view.setParent(None)
            self.view.deleteLater()  # its editor connections go with it
            self.view = None
        self.view_name = name
        self.picker.set_current(name)
        self.wip.set_key(wip.VIEWS.get(name))
        self._update_helpers()
        self.workspace.refresh_names()
        if self.isVisible():
            self.build()
        self._add_companions(name)

    def _add_companions(self, name):
        """The views a tool is read beside (COMPANIONS) that its tab does not show yet: the first in a new pane on the
        right, the next ones below it. Closed later, they stay closed (nothing adds them back until the tool is picked)."""
        root = self
        while root is not None and self.workspace.indexOf(root) < 0:
            root = root.parentWidget()
        shown = {p.view_name for p in self.workspace.panes(root)} if root is not None else set()
        beside, orientation = self, Qt.Horizontal
        for view in COMPANIONS.get(name, ()):
            if view not in shown:
                beside = self.workspace.split(beside, orientation)
                beside.set_view(view)
                orientation = Qt.Vertical


class TabBar(QTabBar):
    """Tab bar with a "+" button right after the last tab, like a web browser. Right-click a
    tab: rename / duplicate / close; middle-click: duplicate."""

    def __init__(self, workspace):
        super().__init__()
        self.workspace = workspace
        self.add = QToolButton(self, text="+", toolTip="New tab", autoRaise=True)
        font = self.add.font()
        font.setPointSizeF(font.pointSizeF() * 1.4)
        self.add.setFont(font)
        # "+" asks what the new tab should hold: this picker is never shown, only its popup (under the "+"). Nothing is
        # opened until a view is chosen, and closing the popup opens nothing.
        self.picker = ViewPicker(grouped_views)
        self.picker.chosen.connect(workspace.new_tab)
        self.add.clicked.connect(lambda: self.picker.open_popup(self.add))

    def mouseReleaseEvent(self, event):
        i = self.tabAt(event.position().toPoint())
        if event.button() == Qt.MiddleButton and i >= 0:
            self.workspace.duplicate_tab(i)
        else:
            super().mouseReleaseEvent(event)

    def contextMenuEvent(self, event):
        i = self.tabAt(event.pos())
        if i < 0:
            return
        menu = QMenu(self)
        menu.addAction("Rename…", lambda: self.workspace.rename_tab(i))
        menu.addAction("Duplicate", lambda: self.workspace.duplicate_tab(i)).setEnabled(self.workspace.is_layout(i))
        menu.addAction("Close", lambda: self.workspace.close_tab(i)).setEnabled(self.workspace.can_close(i))
        menu.exec(event.globalPos())

    def _side(self):
        r = self.tabRect(0) if self.count() else self.rect()
        return max(16, r.height() - 6)

    def sizeHint(self):
        """Room for the "+" after the tabs, so the tab widget lays the bar out with it (it is not a tab Qt knows about)."""
        hint = super().sizeHint()
        hint.setWidth(hint.width() + self._side() + 8)
        return hint

    def minimumSizeHint(self):
        hint = super().minimumSizeHint()
        hint.setWidth(hint.width() + self._side() + 8)
        return hint

    def tabInserted(self, index):   # a tab added or taken out: which close buttons can work changes (see refresh_names)
        super().tabInserted(index)
        QTimer.singleShot(0, self.workspace.refresh_names)

    def tabRemoved(self, index):
        super().tabRemoved(index)
        QTimer.singleShot(0, self.workspace.refresh_names)

    def tabLayoutChange(self):  # Qt calls this whenever tabs are added, removed, moved or resized
        super().tabLayoutChange()
        self._place_add()

    def resizeEvent(self, event):
        super().resizeEvent(event)
        self._place_add()

    def _place_add(self):
        """The "+" right after the last tab; when the tabs no longer fit and scroll, pinned just left of the scroll arrows (it
        used to be left after the last tab, out of sight past the arrows, as soon as there were more tabs than room)."""
        side = self._side()
        r = self.tabRect(self.count() - 1) if self.count() else self.rect()
        x = (r.right() if self.count() else 0) + 4
        arrows = [b for b in self.findChildren(QToolButton) if b is not self.add and b.isVisible()]
        limit = min([b.x() for b in arrows], default=self.width()) - side - 2
        self.add.setAutoFillBackground(x > limit)   # over the end of a tab: on the bar's own colour, not on the tab's text
        self.add.setGeometry(min(x, limit), (r.center().y() if self.count() else self.height() // 2) - side // 2, side, side)
        self.add.raise_()


class Workspace(QTabWidget):
    """The tabs of a window. Most tabs are a tree of split panes; a tool tab (the Assembly panel) holds one widget
    that is not a layout: it is not saved with the layouts, not duplicated, and closing it only takes it out
    (its owner keeps it, so reopening it finds it as it was). `tool_only`: a workspace with no panes at all, for a
    window with no file open yet."""
    tool_closed = Signal()

    def __init__(self, window, tabs=None, tool_only=False):
        super().__init__()
        self.window = window
        self.setTabBar(TabBar(self))
        if tool_only:
            self.tabBar().add.hide()
        self.links = []  # (signal, slot) sequence view -> genome view, redone by rewire()
        self.setTabsClosable(True)
        self.setMovable(True)
        self.setDocumentMode(True)
        self.tabCloseRequested.connect(self.close_tab)
        self.tabBarDoubleClicked.connect(self.rename_tab)
        # The quick display filters at the right end of the tab row, as wide as their buttons (the tabs keep the rest).
        self.quick_filters = None
        if not tool_only and getattr(window, "editor", None) is not None:
            self.quick_filters = QuickFilters(window, window.editor)
            self.setCornerWidget(self.quick_filters, Qt.TopRightCorner)
            self.currentChanged.connect(self._quick_filters_for)
        for tab in [] if tool_only else tabs or DEFAULT:
            name = tab.get("name")
            if name in OLD_TAB_NAMES or re.fullmatch(r"Tab \d+", name or ""):
                name = None   # saved by an older version as a fixed name: automatic now
            self.add_tab(name, tab.get("layout", pane("Sequence")))

    def _quick_filters_for(self, i):
        """The quick filters act on the document's views: greyed on a tab that shows none of them (Assembly, a database). Never
        hidden, which would give the tabs their room and take it back at every switch, moving them."""
        if self.quick_filters is not None:
            self.quick_filters.setEnabled(i < 0 or self.is_layout(i))

    # ---- building / saving -------------------------------------------------------------
    def build(self, spec):
        if "split" in spec:
            s = QSplitter(Qt.Horizontal if spec["split"] == "h" else Qt.Vertical)
            s.setChildrenCollapsible(False)
            for child in spec.get("children", []):
                s.addWidget(self.build(child))
            if spec.get("sizes"):
                s.setSizes(spec["sizes"])
            return s
        return Pane(self, spec.get("view", "Sequence"))

    def dump_widget(self, w):
        if isinstance(w, QSplitter):
            return split("h" if w.orientation() == Qt.Horizontal else "v", w.sizes(),
                         *[self.dump_widget(w.widget(i)) for i in range(w.count())])
        return pane(w.view_name)

    def dump(self):
        """The layouts of the tabs of panes; a tab only keeps its name when it was given one by hand."""
        out = []
        for i in range(self.count()):
            if self.is_layout(i):
                custom = getattr(self.widget(i), "custom_name", None)
                out.append({**({"name": custom} if custom else {}), "layout": self.dump_widget(self.widget(i))})
        return out

    # ---- tabs ----------------------------------------------------------------------------
    def is_layout(self, i):
        """True for a tab of split panes, False for a tool tab."""
        return isinstance(self.widget(i), QSplitter)

    def can_close(self, i):
        return not self.is_layout(i) or sum(map(self.is_layout, range(self.count()))) > 1

    def take_tool(self, widget):
        """Take a tool tab out without deleting its widget (its owner keeps it)."""
        i = self.indexOf(widget)
        if i == self.currentIndex() and self.count() > 1:
            # Removing the shown tab directly crashes inside Qt (access violation, PySide 6 on Windows):
            # show a neighbour first.
            self.setCurrentIndex(i - 1 if i else 1)
        self.removeTab(self.indexOf(widget))
        widget.setParent(None)

    def show_tool(self, widget, title):
        """Add `widget` as a tab (unless it already is one) and show it."""
        if self.indexOf(widget) < 0:
            self.addTab(widget, title)
        self.setCurrentWidget(widget)

    def add_tab(self, name, spec, select=False):
        """`name`: a name given by hand, or None to name the tab after the views in it."""
        root = self.build(spec)
        if not isinstance(root, QSplitter):  # the root is always a splitter (see close_pane)
            s = QSplitter()
            s.addWidget(root)
            root = s
        root.custom_name = name or None
        i = self.addTab(root, name or "")
        self.refresh_names()
        if select:
            self.setCurrentIndex(i)
        return root

    def ordered_panes(self, w):
        """The panes under `w` in layout order (left to right, top to bottom)."""
        if isinstance(w, QSplitter):
            return [p for i in range(w.count()) for p in self.ordered_panes(w.widget(i))]
        return [w] if isinstance(w, Pane) else []

    def tab_views(self, root):
        """What a tab shows, without the pane headers: the views of its visible panes, or the tool of a tool tab."""
        if not isinstance(root, QSplitter):
            return [root]
        return [p.view for p in self.ordered_panes(root) if p.view is not None and p.isVisibleTo(root)]

    def named_panes(self, root):
        """The panes a tab's name describes: those you can actually see. The Alignment pane of a short sequence is
        hidden (View → Alignment View brings it back), and a name promising a view that is not on screen is a puzzle.
        `isVisibleTo` asks "would this show if the window were open?", so naming works before the window is shown.
        If every pane is hidden, fall back to all of them rather than calling the tab "Empty"."""
        panes = self.ordered_panes(root)
        return [p for p in panes if p.isVisibleTo(root)] or panes

    def auto_name(self, root):
        """The views of a tab, in layout order: "Sequence · Map · Features", with a count after the first few."""
        names = []
        for p in self.named_panes(root):
            if p.view_name not in names:
                names.append(p.view_name)
        shown = " · ".join(names[:MAX_NAMED]) or "Empty"
        return shown + (f" +{len(names) - MAX_NAMED}" if len(names) > MAX_NAMED else "")

    def refresh_names(self):
        """Tabs without a name of their own are named after the views they hold (the full list is the tooltip). Also greys the
        close buttons that would do nothing: a pane's ✕ when it is the last pane of its tab, a tab's when it is the last layout
        tab (greyed, not hidden: hiding would change the tab's width at every switch)."""
        for i in range(self.count()):
            closer = self.tabBar().tabButton(i, QTabBar.RightSide)
            if closer is not None:
                closer.setEnabled(self.can_close(i))
                closer.setToolTip("Close this tab" if self.can_close(i) else "The last tab cannot be closed")
            if self.is_layout(i):
                panes = self.panes(self.widget(i))
                for p in panes:
                    p.close_button.setEnabled(len(panes) > 1)
                    p.close_button.setToolTip("Close this pane" if len(panes) > 1 else
                                              "The only pane of this tab (close the tab instead)")
                root = self.widget(i)
                auto = self.auto_name(root)
                self.setTabText(i, getattr(root, "custom_name", None) or auto)
                self.setTabToolTip(i, " · ".join(dict.fromkeys(p.view_name for p in self.named_panes(root))))

    def new_tab(self, view):
        """A new tab with just `view` in it (an Assembly choice shows its own tab instead, as in a pane's picker)."""
        if view in TAB_TOOLS:
            TAB_TOOLS[view](self.window)
        else:
            self.add_tab(None, pane(view), select=True)

    def duplicate_tab(self, i):
        """Same layout (views and sizes) in a new tab right after it."""
        if not self.is_layout(i):
            return
        custom = getattr(self.widget(i), "custom_name", None)
        root = self.add_tab(custom + " copy" if custom else None, self.dump_widget(self.widget(i)))
        self.tabBar().moveTab(self.indexOf(root), i + 1)
        self.setCurrentIndex(i + 1)

    def close_tab(self, i):
        if not self.is_layout(i):   # a tool: take it out, its owner keeps it
            self.take_tool(self.widget(i))
            self.tool_closed.emit()
        elif self.can_close(i):
            w = self.widget(i)
            if i == self.currentIndex():
                # Removing the shown tab directly crashes inside Qt (access violation, PySide 6
                # on Windows); showing a neighbour first avoids it.
                self.setCurrentIndex(i - 1 if i else 1)
            self.removeTab(self.indexOf(w))
            w.deleteLater()

    def rename_tab(self, i):
        """Give a tab a name of its own; an empty name (or the automatic one) goes back to naming it after its views."""
        if i >= 0:
            name, ok = QInputDialog.getText(self, "Rename tab", "Name (empty = named after its views):", text=self.tabText(i))
            problem = validate.text_problem(name.strip(), "tab name", MAX_TAB_NAME) if ok else ""
            if problem:
                QMessageBox.warning(self, "Rename tab", problem)
            elif ok:
                self.set_tab_name(i, name)

    def set_tab_name(self, i, name):
        name = " ".join(name.split())[:MAX_TAB_NAME]   # one line, and short enough not to push the other tabs away
        if self.is_layout(i):
            root = self.widget(i)
            root.custom_name = None if not name or name == self.auto_name(root) else name
            self.refresh_names()
        elif name:
            self.setTabText(i, name)

    # ---- panes ---------------------------------------------------------------------------
    def split(self, p, orientation):
        """Add a pane next to `p` (same view as `p`, change it with the picker)."""
        parent = p.parentWidget()
        new = Pane(self, p.view_name)
        i = parent.indexOf(p)
        if parent.orientation() == orientation or parent.count() == 1:
            parent.setOrientation(orientation)
            sizes = parent.sizes()
            parent.insertWidget(i + 1, new)
            half = sizes[i] // 2
            parent.setSizes(sizes[:i] + [half, sizes[i] - half] + sizes[i + 1:])
        else:  # wrap p in a splitter of the other direction
            sizes = parent.sizes()
            s = QSplitter(orientation)
            s.setChildrenCollapsible(False)
            parent.insertWidget(i, s)
            s.addWidget(p)
            s.addWidget(new)
            parent.setSizes(sizes)
        self.refresh_names()
        return new

    def close_pane(self, p):
        parent = p.parentWidget()
        if len(self.panes(self.currentWidget())) <= 1 and parent is self.currentWidget():
            return  # the last pane of a tab: close the tab instead
        p.setParent(None)
        p.deleteLater()
        # A splitter left with one child (not a tab root) is replaced by that child.
        if parent.count() == 1 and isinstance(parent.parentWidget(), QSplitter):
            grand = parent.parentWidget()
            sizes = grand.sizes()
            grand.insertWidget(grand.indexOf(parent), parent.widget(0))
            parent.setParent(None)
            parent.deleteLater()
            grand.setSizes(sizes)
        self.refresh_names()
        self.rewire()

    def panes(self, root=None):
        roots = [root] if root is not None else [self.widget(i) for i in range(self.count())]
        return [p for r in roots for p in r.findChildren(Pane)]

    def views(self, cls, build=False):
        """Views of type `cls`, current tab first. `build`: if there is none yet, create the
        first pane's view that was never shown ("all" = create every such view)."""
        current = self.panes(self.currentWidget())
        ordered = current + [p for p in self.panes() if p not in current]
        found = []
        for p in ordered:
            if p.view is None and build and p.view_name == NAME_OF.get(cls) and (build == "all" or not found):
                p.build()
            if isinstance(p.view, cls):
                found.append(p.view)
        return found

    def ensure(self, cls):
        """The first view of `cls`, adding a pane for it at the top of the current tab if
        there is none."""
        found = self.views(cls, build=True)
        if found:
            return found[0]
        root = self.currentWidget()
        if not isinstance(root, QSplitter):   # a tool tab is shown: add to the first tab of panes
            root = next(self.widget(i) for i in range(self.count()) if self.is_layout(i))
        p = Pane(self, NAME_OF[cls])
        if root.orientation() == Qt.Vertical:
            root.insertWidget(0, p)
        else:
            s = QSplitter(Qt.Horizontal)
            s.setChildrenCollapsible(False)
            while root.count():
                s.addWidget(root.widget(0))
            root.setOrientation(Qt.Vertical)
            root.addWidget(p)
            root.addWidget(s)
        self.refresh_names()
        return p.build()

    def rewire(self):
        """Each genome view frames what the first sequence view of its tab shows (or of the
        window, if its tab has none)."""
        for signal, slot in self.links:
            try:
                signal.disconnect(slot)
            except (RuntimeError, TypeError, SystemError):
                pass  # already gone with its view (PySide reports a deleted one as SystemError)
        self.links = []
        everywhere = self.views(SequenceView)
        for i in range(self.count()):
            tab = [p.view for p in self.panes(self.widget(i)) if p.view is not None]
            seqs = [v for v in tab if isinstance(v, SequenceView)] or everywhere
            for g in (v for v in tab if isinstance(v, GenomeView)):
                if seqs:
                    seqs[0].window_changed.connect(g.set_sequence_window)
                    self.links.append((seqs[0].window_changed, g.set_sequence_window))


NAME_OF = {SequenceView: "Sequence", MapView: "Map", GenomeView: "Alignment", FeatureTable: "Features",
           PrimerPanel: "Primers", ProteinView: "Protein"}
