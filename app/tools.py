"""Tools: PCR product, Detect ORFs, restriction digest with a virtual gel."""
import math

from PySide6.QtCore import QEvent, QObject, QPointF, QRectF, Qt, QTimer, Signal
from PySide6.QtGui import QColor, QPainter
from PySide6.QtWidgets import (
    QAbstractItemView, QComboBox, QDialog, QDialogButtonBox, QFormLayout, QHBoxLayout, QLabel,
    QLineEdit, QListWidget, QMenu, QPushButton, QSizePolicy, QSpinBox, QSplitter, QStyle,
    QTabBar, QTableWidget, QTableWidgetItem, QToolButton, QTreeWidget, QVBoxLayout, QWidget,
)

from editor import protein_field, reverse_complement, site_location
from cut_sites import enzyme_mark_tip
import inventory
from enzymes import COMMON, EnzymeFilter, describe, find_cuts, wanted
import genetic_code
from sequence_view import cds_codons, codons
from uikit import fit_columns, no_auto_default, plural, right_aligned, run_busy
import value_filter


class ToolDialog(QDialog):
    """A tool that opens as a dialog, or lives in a workspace pane after `embed()`. In a
    pane, OK runs `on_accept(self)` instead of closing, and the tool follows the document:
    `refresh()` runs after edits / primer changes (once typing pauses, and only while the
    pane is visible; a hidden one refreshes when shown again)."""
    embedded, stale = False, False

    def embed(self, on_accept=None):
        self.embedded, self.on_accept = True, on_accept
        self.setWindowFlags(Qt.Widget)
        timer = QTimer(self, singleShot=True, interval=400, timeout=self._maybe_refresh)
        self.editor.changed.connect(timer.start)
        self.editor.primers_changed.connect(timer.start)
        return self

    def done(self, result):  # accept()/reject() end here (C++ virtual, so this override works)
        if not self.embedded:
            return super().done(result)
        if result == QDialog.Accepted and self.on_accept:
            self.on_accept(self)

    def _maybe_refresh(self):
        if self.isVisible():
            self.refresh()
        else:
            self.stale = True

    def showEvent(self, event):
        if self.stale:
            self.stale = False
            self.refresh()
        super().showEvent(event)

    def refresh(self):
        pass


class _EmptyHint(QObject):
    """Keeps an empty table's hint centred over its rows' area, and shown only while there is no row."""

    def __init__(self, view, label):
        super().__init__(view)
        self.view, self.label = view, label
        view.viewport().installEventFilter(self)
        model = view.model()
        for signal in (model.rowsInserted, model.rowsRemoved, model.modelReset, model.layoutChanged):
            signal.connect(self.update)
        self.update()

    def update(self, *_):
        try:
            self.label.setGeometry(self.view.viewport().rect().adjusted(12, 12, -12, -12))
            self.label.setVisible(self.view.model().rowCount() == 0)
        except RuntimeError:
            pass   # the table is being deleted (its model says its rows go)

    def eventFilter(self, obj, event):
        if event.type() == QEvent.Resize:
            self.update()
        return False


def empty_hint(view, text):
    """A grey line in the middle of a table while it has no row ("Paste a query and press Align"): what to do, where the
    result will be, instead of a blank grid or a permanent line of instructions. Returns the label (its text can change)."""
    label = QLabel(text, view.viewport(), alignment=Qt.AlignCenter, wordWrap=True)
    label.setEnabled(False)   # greyed, like the database tabs' "Click a feature to see its sequence."
    label.setAttribute(Qt.WA_TransparentForMouseEvents)
    view.empty_hint = _EmptyHint(view, label)
    return label

# ---- PCR ------------------------------------------------------------------------------


def product_range(editor, f, r):
    """(start, end) of the product on the template, end past the origin if it wraps; None
    if the reverse site is not downstream of the forward one."""
    n = editor.length
    start, end = f.start, r.end
    if end <= start:
        if not editor.doc.circular:
            return None
        end += n
    return (start, end) if end - start <= n else None


def pcr_product(editor, f, r):
    """New unsaved Document: forward primer + template between the primers + reverse
    complement of the reverse primer. Primer bases replace template bases where they
    differ (mismatched sites mutate the product, as in real PCR)."""
    start, end = product_range(editor, f, r)
    name = f"{editor.doc.name}_{f.label}_{r.label}"
    prefix = f.sequence[:f.tail]
    suffix = reverse_complement(r.sequence[:r.tail])
    doc = editor.doc.pcr_product(name, start, end, prefix, suffix)
    if f.mismatches:
        doc.overwrite(0, f.sequence)
    if r.mismatches:
        rc = reverse_complement(r.sequence)
        doc.overwrite(doc.length - len(rc), rc)
    return doc


class PcrDialog(ToolDialog):
    def __init__(self, parent, editor):
        super().__init__(parent)
        self.editor = editor
        self.setWindowTitle("PCR Product")
        self.setMinimumWidth(480)
        self.fwd, self.rev = QComboBox(), QComboBox()
        for combo in (self.fwd, self.rev):
            combo.currentIndexChanged.connect(self._update)
            combo.setMaximumWidth(400)   # in a wide pane the form stays a form, on the left, not 1400 px of combo box
        self.info = QLabel(wordWrap=True, maximumWidth=400)
        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel, maximumWidth=400)
        buttons.button(QDialogButtonBox.Ok).setText("Open product")
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        self.ok = buttons.button(QDialogButtonBox.Ok)
        self.cancel = buttons.button(QDialogButtonBox.Cancel)
        form = QFormLayout(self)
        form.setFormAlignment(Qt.AlignLeft | Qt.AlignTop)
        form.addRow("Forward primer", self.fwd)
        form.addRow("Reverse primer", self.rev)
        form.addRow("Product", self.info)
        form.addRow("", buttons)
        self.refresh()

    def embed(self, on_accept=None):
        self.cancel.hide()   # in a pane there is nothing to cancel: the pane stays
        return super().embed(on_accept)

    def refresh(self):
        """(Re)list the primer sites, keeping the chosen pair if it still binds."""
        before = (self.fwd.currentText(), self.rev.currentText())
        self.fwd_sites = [s for s in self.editor.primer_sites if s.strand > 0]
        self.rev_sites = [s for s in self.editor.primer_sites if s.strand < 0]
        for combo, sites in ((self.fwd, self.fwd_sites), (self.rev, self.rev_sites)):
            combo.blockSignals(True)
            combo.clear()
            for s in sites:
                mm = f", {plural(s.mismatches, 'mismatch')}" if s.mismatches else ""
                combo.addItem(f"{s.label}  ({s.start + 1:,}{mm})")
            combo.blockSignals(False)
        if before[0] and self.fwd.findText(before[0]) >= 0 and self.rev.findText(before[1]) >= 0:
            self.fwd.setCurrentText(before[0])
            self.rev.setCurrentText(before[1])
        else:
            self._pick_default()
        self._update()

    def _pick_default(self):
        """The pair giving the shortest product (around the selection if there is one)."""
        best = None
        for i, f in enumerate(self.fwd_sites):
            for j, r in enumerate(self.rev_sites):
                span = product_range(self.editor, f, r)
                if span and (best is None or span[1] - span[0] < best[0]):
                    best = (span[1] - span[0], i, j)
        if best:
            self.fwd.setCurrentIndex(best[1])
            self.rev.setCurrentIndex(best[2])

    def pair(self):
        i, j = self.fwd.currentIndex(), self.rev.currentIndex()
        if i < 0 or j < 0:
            return None
        return self.fwd_sites[i], self.rev_sites[j]

    def _update(self):
        pair = self.pair()
        span = product_range(self.editor, *pair) if pair else None
        if not span:
            self.info.setText("Needs a forward and a reverse primer site, reverse downstream." if not pair
                              else "No product: the reverse primer is not downstream.")
        else:
            f, r = pair
            size = span[1] - span[0] + f.tail + r.tail
            self.info.setText(f"{size:,} bp (template {span[0] + 1:,}–{(span[1] - 1) % self.editor.length + 1:,}"
                              f" + {f.tail + r.tail} bp of tails)")
        self.ok.setEnabled(bool(span))

    def product(self):
        return pcr_product(self.editor, *self.pair())


# ---- ORFs -----------------------------------------------------------------------------


class OrfDialog(ToolDialog):
    def __init__(self, parent, editor):
        super().__init__(parent)
        self.editor, self.orfs = editor, []
        self.setWindowTitle("Detect ORFs")
        self.resize(520, 480)
        # The ORFs are found again once the number stops changing (not at every step of the arrows, a genome takes a moment).
        self.min_aa = QSpinBox(minimum=10, maximum=5000, value=100, suffix=" aa", keyboardTracking=False)
        self._reload_later = QTimer(self, singleShot=True, interval=300, timeout=self.reload)
        self.min_aa.valueChanged.connect(self._reload_later.start)
        self.table = QTableWidget(0, 4)
        self.table.setHorizontalHeaderLabels(["Location", "Strand", "Length (aa)", "Annotation"])
        value_filter.install(self.table, ["Strand", "Annotation"])   # right-click a header: filter by values
        self.table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        fit_columns(self.table, stretch=3)
        self.table.currentCellChanged.connect(lambda row, *_: self._select(row))
        self.add_button = add = QPushButton("Add selected as CDS features", clicked=self.add, enabled=False)
        self.table.itemSelectionChanged.connect(self._can_add)
        top = QHBoxLayout()
        top.addWidget(QLabel("Minimum length"))
        top.addWidget(self.min_aa)
        self.count = QLabel()   # how many found, how many not annotated yet (greyed: a detail)
        self.count.setEnabled(False)
        top.addWidget(self.count)
        top.addStretch()
        # Protein motif search in the translations of the CDSs and of the ORFs listed.
        self.motif = protein_field(QLineEdit(placeholderText="Find in translations (protein motif, e.g. MKKL)",
                                             clearButtonEnabled=True))
        self.motif.textChanged.connect(self.search)
        self.hits = QListWidget()
        self.hits.setMaximumHeight(130)
        self.hits.currentRowChanged.connect(self._go_hit)
        self.found = []  # (start, end) on the sequence of each hit row
        layout = QVBoxLayout(self)
        layout.addLayout(top)
        layout.addWidget(self.table)
        layout.addWidget(add)
        layout.addWidget(self.motif)
        layout.addWidget(self.hits)
        no_auto_default(self)   # Enter in the length box does not add CDS features
        self.reload()

    def refresh(self):
        self.reload()

    def _can_add(self):
        """The button can add something only when a selected ORF is not annotated yet."""
        rows = {i.row() for i in self.table.selectionModel().selectedRows()}
        self.add_button.setEnabled(any(0 <= r < len(self.links) and self.links[r] is None for r in rows))

    def annotated(self, s, e, strand):
        """The CDS feature already annotating ORF [s, e) (same strand, same stop), or None."""
        for f in self.editor.features:
            if f.kind == "CDS" and f.strand == strand and (
                    max(b for _, b in f.spans) == e if strand > 0 else min(a for a, _ in f.spans) == s):
                return f
        return None

    def reload(self):
        # +1 codon: the stop codon is counted by the finder.
        # The genetic code and "alternative start codons" come from Tools > Settings > Genetic code.
        doc, args = self.editor.doc, (self.min_aa.value() + 1, genetic_code.table_id(), genetic_code.alt_starts())
        self.orfs = run_busy(self, "Finding ORFs…", lambda: doc.orfs(*args))   # (a genome takes a moment: the window stays alive)
        self.links = [self.annotated(*orf) for orf in self.orfs]  # annotation of each ORF, or None
        self.table.setRowCount(len(self.orfs))
        link = QColor(self.palette().color(self.palette().ColorRole.Link))
        for i, ((s, e, strand), f) in enumerate(zip(self.orfs, self.links)):
            values = [f"{s + 1:,}–{e:,}", "→" if strand > 0 else "←", f"{(e - s) // 3 - 1:,}",
                      f"→ {f.label} (annotated)" if f else "predicted"]
            for c, v in enumerate(values):
                item = QTableWidgetItem(v)
                if c == 2:
                    right_aligned(item)
                if f:
                    item.setForeground(link)  # a link: selecting the row selects that feature
                self.table.setItem(i, c, item)
        new = sum(f is None for f in self.links)
        self.setWindowTitle("Detect ORFs")   # (the counts beside the length box, not in the title)
        self.count.setText(f"{plural(len(self.orfs), 'ORF')}, {new} not annotated")
        self._can_add()
        self.search()

    def _select(self, row):
        if 0 <= row < len(self.orfs):
            if self.links[row] is not None:
                self.editor.select_feature(self.links[row].index)  # the real annotation
                return
            s, e, _ = self.orfs[row]
            self.editor.select(s, e)
            self.editor.reveal.emit(s)

    def add(self):
        rows = sorted({i.row() for i in self.table.selectionModel().selectedRows()})
        features = [("CDS", site_location(s, e, strand, self.editor.length),
                     [("label", f"ORF {(e - s) // 3 - 1} aa"), ("codon_start", "1")])
                    for i in rows if self.links[i] is None  # annotated ones are there already
                    for s, e, strand in [self.orfs[i]]]
        if features:
            self.editor.add_features(features)

    def search(self):
        """Hits of the protein motif in each CDS's translation, then in each ORF not
        annotated as a CDS (shortcut: exact letters, X matches only X)."""
        motif = self.motif.text().strip().upper()
        self.hits.clear()
        self.found = []
        if not motif:
            return
        ed = self.editor
        table = genetic_code.table_id()
        # A CDS is read as annotated (its /transl_table, /codon_start, /transl_except); an ORF with the chosen table.
        frames = [(f.label, cds_codons(ed.doc, f)) for f in ed.features if f.kind == "CDS"]
        frames += [(f"ORF {s + 1:,}–{e:,}", codons(ed.doc, [(s, e)], strand, 0, table))
                   for (s, e, strand), f in zip(self.orfs, self.links) if f is None]
        for name, cods in frames:
            protein = "".join(aa for _, aa, _, _ in cods)
            i = protein.find(motif)
            while i >= 0 and len(self.found) < 500:
                part = cods[i:i + len(motif)]
                self.found.append((min(c[2] for c in part), max(c[3] for c in part)))
                self.hits.addItem(f"{name}: {motif} at aa {i + 1}")
                i = protein.find(motif, i + 1)
        if not self.found:
            self.hits.addItem("No match")

    def _go_hit(self, row):
        if 0 <= row < len(self.found):
            s, e = self.found[row]
            self.editor.select(s, e)
            self.editor.reveal.emit(s)


# ---- digest + gel ----------------------------------------------------------------------

# Ladders as marketed in general terms: the bands of the usual sizes, not the exact list of one vendor. name -> (bands in bp, the bright band)
LADDERS = {
    "1 kb": ([10000, 8000, 6000, 5000, 4000, 3000, 2500, 2000, 1500, 1000, 750, 500, 250], 3000),
    "1 kb Plus": ([25000, 20000, 15000, 10000, 8000, 6000, 5000, 4000, 3000, 2000, 1500, 1000, 750, 500, 250], 5000),
    "100 bp": ([1500, 1200, 1000, 900, 800, 700, 600, 500, 400, 300, 200, 100], 500),
    "100 bp Plus": ([3000, 2000, 1500, 1200, 1000, 900, 800, 700, 600, 500, 400, 300, 200, 100], 500),
    "50 bp": ([500, 450, 400, 350, 300, 250, 200, 150, 100, 50], 250),
}
DEFAULT_LADDER = "1 kb"


def fragments(length, circular, cuts):
    """[(size, start, end, left enzymes, right enzymes)] for cut positions {pos: [names]}.
    Circular: as many fragments as cuts; linear: one more."""
    pos = sorted(cuts)
    if not pos:
        return [(length, 0, length, [], [])]
    if circular:
        out = []
        for i, a in enumerate(pos):
            b = pos[(i + 1) % len(pos)]
            size = (b - a) % length or length
            out.append((size, a, b, cuts[a], cuts[b]))
        return out
    edges = [0] + pos + [length]
    return [(b - a, a, b, cuts.get(a, []), cuts.get(b, [])) for a, b in zip(edges, edges[1:]) if b > a]


class GelView(QWidget):
    """Agarose gel sketch: band distance ~ log(size). The ladder is the first lane, its name in brackets: click it to choose
    another; the scale follows the ladder and the digests."""

    ladder_changed = Signal(str)

    def __init__(self):
        super().__init__()
        self.ladder, self.digests = DEFAULT_LADDER, []
        self.setMinimumSize(220, 360)
        self.setMouseTracking(True)

    @property
    def lanes(self):
        """[(name, band sizes)]: the ladder, then each digest."""
        return [(f"Ladder ({self.ladder})", LADDERS[self.ladder][0])] + self.digests

    def set_lanes(self, lanes):
        self.digests = list(lanes)
        self.setMinimumWidth(max(220, 50 + 56 * (len(lanes) + 1)))
        self.update()

    def set_ladder(self, name):
        if name in LADDERS and name != self.ladder:
            self.ladder = name
            self.update()
            self.ladder_changed.emit(name)

    def lane_width(self):
        return (self.width() - 40) / max(1, len(self.digests) + 1)

    def over_ladder(self, x):
        return 30 <= x < 30 + self.lane_width()

    def menu(self):
        """The ladders to choose from; the current one is ticked."""
        menu = QMenu(self)
        for name in LADDERS:
            lo, hi = min(LADDERS[name][0]), max(LADDERS[name][0])
            action = menu.addAction(f"{name}   ({lo:,}–{hi:,} bp)", lambda n=name: self.set_ladder(n))
            action.setCheckable(True)
            action.setChecked(name == self.ladder)
        return menu

    def mouseMoveEvent(self, event):
        over = self.over_ladder(event.position().x())
        self.setCursor(Qt.PointingHandCursor if over else Qt.ArrowCursor)
        self.setToolTip("Click to choose another ladder" if over else "")

    def mousePressEvent(self, event):
        if event.button() == Qt.LeftButton and self.over_ladder(event.position().x()):
            self.menu().exec(event.globalPosition().toPoint())

    LINES = 3   # lines of enzyme names over a lane (more: "+n more")

    def lane_lines(self, name):
        """The text over a lane, one line each: a digest's enzymes one per line (the last line counts the rest), the ladder's name
        over its size range."""
        if name.startswith("Ladder ("):
            head, _, rest = name.partition(" (")
            return [head, "(" + rest]
        names = name.split("+")
        if len(names) > self.LINES:
            names = names[:self.LINES - 1] + [f"+{len(names) - self.LINES + 1} more"]
        return names

    def paintEvent(self, event):
        p = QPainter(self)
        p.fillRect(self.rect(), QColor("#1b1b1b"))
        lines = max(len(self.lane_lines(name)) for name, _ in self.lanes)
        top, bottom = 20 + 14 * lines, self.height() - 12
        bands, bright_band = LADDERS[self.ladder]
        sizes = bands + [s for _, lane in self.digests for s in lane]
        hi, lo = math.log10(max(sizes) * 1.15), math.log10(min(sizes) / 1.15)

        def y(size):
            return top + (hi - math.log10(size)) / (hi - lo) * (bottom - top)

        width = self.lane_width()
        small = p.font()
        small.setPointSizeF(max(7.0, small.pointSizeF() - 1))
        p.setFont(small)
        for i, (name, lane) in enumerate(self.lanes):
            x0 = 30 + i * width + 6
            p.setPen(QColor("#bbbbbb"))
            for k, line in enumerate(self.lane_lines(name)):   # elided to the lane when too long
                p.drawText(QRectF(x0 - 6, 3 + 14 * k, width, 15), Qt.AlignCenter, p.fontMetrics().elidedText(line, Qt.ElideRight, int(width)))
            p.fillRect(QRectF(x0, top - 8, width - 12, 5), QColor("#555555"))  # well
            for size in lane:
                if i == 0:
                    bright, thick = (255, 4) if size == bright_band else (170, 3)
                else:
                    bright, thick = 120 + min(135, int(135 * size / 3000)), 3  # bigger fragments carry more DNA
                p.fillRect(QRectF(x0, y(size) - thick / 2, width - 12, thick), QColor(bright, bright, bright))
        p.setPen(QColor("#888888"))
        last = -99.0
        for size in bands:   # the ladder's sizes down the left, as many as fit
            if y(size) - last >= 12:
                p.drawText(QPointF(2, y(size) + 4), f"{size / 1000:g}k" if size >= 1000 else str(size))
                last = y(size)
        p.end()


def _roomy(criteria):
    """The Digest's filter bar is a column of its own (the Display panel's copy of these criteria is squeezed into a side
    panel): roomier number boxes, and the presets as wide as their text."""
    for box in criteria.findChildren(QSpinBox):
        box.setFixedWidth(max(box.width(), box.fontMetrics().horizontalAdvance("8888") + 40))
    for combo in criteria.findChildren(QComboBox):
        combo.setFixedWidth(criteria.min_sites.maximumWidth())
    for button in criteria.findChildren(QPushButton):
        button.setStyleSheet("")
        button.setMinimumWidth(0)
        button.setMaximumWidth(16777215)   # (QWIDGETSIZE_MAX: as wide as its text in the normal font)
        button.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Fixed)


class LaneBar(QTabBar):
    """The Digest's lanes, one tab each, with a "+" right after the last tab (as the window's own tabs have)."""

    def __init__(self):
        super().__init__()
        self.add = QToolButton(self, text="+", toolTip="Add another digest (lane)", autoRaise=True)
        font = self.add.font()
        font.setPointSizeF(font.pointSizeF() * 1.4)
        self.add.setFont(font)
        self.setExpanding(False)
        self.setDrawBase(True)

    def tabLayoutChange(self):  # Qt calls this whenever tabs are added, removed, moved or resized
        super().tabLayoutChange()
        self._place()

    def resizeEvent(self, event):
        super().resizeEvent(event)
        self._place()

    def _place(self):
        r = self.tabRect(self.count() - 1) if self.count() else self.rect()
        side = max(16, r.height() - 6)
        self.add.setGeometry((r.right() if self.count() else 0) + 4, r.center().y() - side // 2, side, side)

    def minimumSizeHint(self):
        hint = super().minimumSizeHint()
        return hint.expandedTo(self.add.sizeHint())


DIGEST_FILTER = EnzymeFilter(1, None, 1)   # the Digest lists every enzyme that cuts, whatever its site, to start with


def ends_of(enzyme):
    """The ends an enzyme leaves, in a word: "5′ 4 nt", "3′ 4 nt" or "blunt" (its first cut, as enzymes.wanted reads it)."""
    top, bottom = enzyme.cuts[0]
    return "blunt" if top == bottom else f"{'5′' if top < bottom else '3′'} {abs(bottom - top)} nt"


class EnzymeTable(QTreeWidget):
    """The Digest's enzymes: a tick (in the selected digest), how many sites here, the site it recognises, the ends it leaves,
    and Owned (in the lab's inventory, inventory.py: the same for every sequence). Click a header to sort."""
    COLUMNS = ["Enzyme", "Sites", "Recognition site", "Ends", "Owned"]
    NAME, SITES, SITE, ENDS, OWNED = range(5)
    owned_changed = Signal(str, bool)
    ticked = Signal(str, bool)

    def __init__(self):
        super().__init__()
        self.setHeaderLabels(self.COLUMNS)
        self.setRootIsDecorated(False)
        self.setUniformRowHeights(True)
        self.setAllColumnsShowFocus(True)
        self.setSortingEnabled(True)
        self.sortByColumn(self.NAME, Qt.AscendingOrder)
        self.header().setStretchLastSection(False)
        self.header().setMinimumSectionSize(16)   # (the style's minimum would widen the narrow columns: Sites, Owned)
        self.headerItem().setToolTip(self.OWNED, "Owned in the lab (your enzyme inventory): tick to tag, the same for every sequence")
        self.headerItem().setToolTip(self.ENDS, "The ends the enzyme leaves: a 5′ or 3′ overhang (its length) or blunt")
        self.itemChanged.connect(self._changed)

    def fill(self, enzymes, cuts, ticked, owned):
        """One row per enzyme of `enzymes` (Enzyme tuples), its sites from `cuts` {name: [cuts]}."""
        self.blockSignals(True)
        self.setSortingEnabled(False)   # (refilling a sorted tree row by row re-sorts it every time)
        self.clear()
        for e in enzymes:
            n = len(cuts.get(e.name, []))
            item = value_filter.SortItem([e.name, str(n), e.site, ends_of(e), ""])
            item.keys = {self.SITES: n}
            item.setData(self.NAME, Qt.UserRole, e.name)
            item.setTextAlignment(self.SITES, Qt.AlignRight | Qt.AlignVCenter)
            item.setFlags(item.flags() | Qt.ItemIsUserCheckable)
            item.setCheckState(self.NAME, Qt.Checked if e.name in ticked else Qt.Unchecked)
            item.setCheckState(self.OWNED, Qt.Checked if e.name in owned else Qt.Unchecked)
            tip = describe(e)
            for c in range(len(self.COLUMNS) - 1):
                item.setToolTip(c, tip)
            self.addTopLevelItem(item)
        self.setSortingEnabled(True)
        self.blockSignals(False)

    ARROW = 14   # room for the sort arrow next to a header's text

    def column_need(self, c):
        """How wide column `c` must be: its widest cell, or its header's text with the style's margins and a sort arrow.
        (Qt's own measure of a header gives every column the room of the sort arrow at the style's size, which the Windows
        style makes large: the list came out about 100 px wider than what it holds.)"""
        header = self.header()
        margin = self.style().pixelMetric(QStyle.PM_HeaderMargin, None, header)
        title = header.fontMetrics().horizontalAdvance(self.headerItem().text(c)) + 2 * margin + self.ARROW
        return max(self.sizeHintForColumn(c) + 4, title)

    def fit(self):
        """Every column as wide as what it shows; the list no wider than its columns."""
        for c in range(len(self.COLUMNS)):
            self.setColumnWidth(c, self.column_need(c))
        width = sum(self.columnWidth(c) for c in range(len(self.COLUMNS))) + self.verticalScrollBar().sizeHint().width() + 2 * self.frameWidth() + 4
        self.setMaximumWidth(width)
        self.setMinimumWidth(min(width, 260))
        return width

    def names(self):
        """The enzymes listed, in the order shown."""
        return [self.topLevelItem(i).data(self.NAME, Qt.UserRole) for i in range(self.topLevelItemCount())]

    def item_of(self, name):
        return next((self.topLevelItem(i) for i in range(self.topLevelItemCount())
                     if self.topLevelItem(i).data(self.NAME, Qt.UserRole) == name), None)

    def _changed(self, item, column):
        name = item.data(self.NAME, Qt.UserRole)
        if column == self.NAME:
            self.ticked.emit(name, item.checkState(self.NAME) == Qt.Checked)
        elif column == self.OWNED:
            self.owned_changed.emit(name, item.checkState(self.OWNED) == Qt.Checked)


class DigestDialog(ToolDialog):
    """Restriction digest: pick enzymes, see the fragments, a virtual gel, and the cuts drawn on the sequence. Each tab is a
    digest with its own enzymes and table and its own lane on the shared gel (+ lane adds one, ✕ removes it)."""

    def __init__(self, parent, editor):
        super().__init__(parent)
        self.editor, self.cuts, self.info = editor, {}, {}
        self.lanes, self.numbers, self.frags = [set()], [1], []   # per digest: its ticked enzymes and its number (its tab's default name)
        self.setWindowTitle("Restriction Digest")
        self.resize(900, 520)
        self.tabs = LaneBar()
        self.tabs.addTab("Digest 1")
        self.tabs.currentChanged.connect(self.show_list)
        self.tabs.tabCloseRequested.connect(self.remove_lane)
        self.add = self.tabs.add   # "+" right after the last lane's tab, where the next one will be
        self.add.clicked.connect(self.add_lane)
        # The filter bar on the left: a name or site, then the same criteria as the Display filter's (but the Digest's own: what
        # is listed here, not what is drawn), every enzyme that cuts to start with.
        from display_panel import EnzymeCriteria
        from enzymes import PRESETS
        self.filter = QLineEdit(placeholderText="Name or site", clearButtonEnabled=True)
        self.criteria = EnzymeCriteria(PRESETS)
        _roomy(self.criteria)
        self.criteria.set_filter(DIGEST_FILTER)
        self.criteria.changed.connect(lambda _f: self.show_list())
        self.note = QLabel(wordWrap=True)
        self.note.setEnabled(False)   # greyed: a detail
        self.enzymes = EnzymeTable()
        self.enzymes.setToolTip("Tick enzymes for the selected digest; click a header to sort")
        self.enzymes.ticked.connect(self._ticked)
        self.enzymes.owned_changed.connect(lambda name, on: inventory.set_owned(name, on))
        # The inventory changed (an Owned box here or in another digest): the list is refilled after that click, not inside it
        # (the row clicked would be deleted under Qt's feet).
        self._list_later = QTimer(self, singleShot=True, interval=0, timeout=self.show_list)
        inventory.notifier.changed.connect(self._list_later.start)   # (a QObject slot: the link goes with the dialog)
        self.filter.textChanged.connect(self.show_list)
        self.table = QTableWidget(0, 4)
        self.table.setHorizontalHeaderLabels(["Size (bp)", "From", "To", "Ends"])
        value_filter.install(self.table, ["Ends"])   # right-click a header: filter by values
        self.table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.table.setSelectionMode(QAbstractItemView.SingleSelection)
        fit_columns(self.table, stretch=3)
        self.table.setToolTip("Click a fragment to select it on the sequence, map and genome view and show the cuts at its two ends")
        empty_hint(self.table, "Tick an enzyme to see its fragments")
        self.table.currentCellChanged.connect(lambda row, *_: self.show_fragment(row))
        self.fragment_span = None   # (start, end) of the fragment selected in the table, while it is
        editor.selection_changed.connect(self._selection_moved)
        self.gel = GelView()
        self.gel.setToolTip("Click the ladder's lane to choose another ladder")
        bar = QWidget()
        left = QVBoxLayout(bar)
        left.setContentsMargins(0, 0, 4, 0)
        left.addWidget(QLabel("<b>Filter</b>"))
        left.addWidget(self.filter)
        left.addWidget(self.criteria)
        left.addWidget(self.note)
        left.addStretch()
        bar.setFixedWidth(max(self.criteria.minimumSizeHint().width(), self.criteria.sizeHint().width()) + 24)   # the criteria, uncramped
        split = QSplitter()
        split.setChildrenCollapsible(False)
        split.addWidget(bar)
        split.addWidget(self.enzymes)
        split.addWidget(self.table)
        split.addWidget(self.gel)
        self.split = split
        layout = QVBoxLayout(self)
        layout.addWidget(self.tabs)
        layout.addWidget(split, 1)  # the lists and gel take all the height
        self.destroyed.connect(lambda _=None, ed=editor: (ed.set_cut_marks("digest", []), ed.set_cut_marks("digest_ends", [])))   # the cuts go with the tool
        no_auto_default(self)   # Enter in the name box does not add a lane
        self.refresh()
        width = self.enzymes.fit()   # the enzyme list as narrow as its columns
        split.setSizes([bar.width(), width, 380, 260])

    @property
    def ticked(self):
        """The enzymes ticked in the selected digest."""
        return self.lanes[self.tabs.currentIndex()]

    def add_lane(self):
        self.lanes.append(set())
        self.numbers.append(max(self.numbers) + 1)
        self.tabs.blockSignals(True)
        self.tabs.setCurrentIndex(self.tabs.addTab(f"Digest {self.numbers[-1]}"))
        self.tabs.blockSignals(False)
        self.tabs.setTabsClosable(True)
        self.show_list()

    def remove_lane(self, i):
        if len(self.lanes) < 2:
            return
        del self.lanes[i], self.numbers[i]
        self.tabs.removeTab(i)   # selects a neighbour, which shows its list
        self.tabs.setTabsClosable(len(self.lanes) > 1)
        self.show_list()

    def refresh(self):
        """Find the sites of the enzymes on the current sequence, keeping the ticked ones."""
        self.info = {e.name: e for e in COMMON}
        self.cuts = find_cuts(self.editor.doc.sequence(), self.editor.doc.circular, COMMON)
        self.sites = {name: sorted({t for t, _ in cuts}) for name, cuts in self.cuts.items()}
        self.show_list()

    def show_list(self):
        """The enzymes that cut, narrowed by the name / site box and the criteria of the filter bar."""
        text, f, mine = self.filter.text().strip().upper(), self.criteria.filter(), inventory.owned()
        # Every common enzyme is looked at: with 0 sites allowed, those that do not cut are listed too.
        keep = lambda e: (not text or text in e.name.upper() or text in e.site) and wanted(f, e, len(self.cuts.get(e.name, [])), e.name in mine)
        shown = [e for e in COMMON if keep(e)]
        self.enzymes.fill(shown, self.cuts, self.ticked, mine)
        self.note.setText(f"{len(shown)} of the {len(COMMON)} common enzymes shown; {len(self.cuts)} cut this sequence.")
        self.reload()

    def _ticked(self, name, on):
        (self.ticked.add if on else self.ticked.discard)(name)
        self.reload()

    def chosen(self, ticked=None):
        """The ticked enzymes (of the selected digest) that cut this sequence (a tick stays when the filter hides the enzyme)."""
        return [name for name in sorted(self.ticked if ticked is None else ticked, key=str.lower) if name in self.cuts]

    def digest(self, ticked):
        """(enzymes, fragments largest first) of one digest."""
        names, cuts = self.chosen(ticked), {}
        for name in names:
            for pos in self.sites[name]:
                cuts.setdefault(pos, []).append(name)
        return names, sorted(fragments(self.editor.length, self.editor.doc.circular, cuts), reverse=True)

    def reload(self):
        digests = [self.digest(t) for t in self.lanes]
        here = self.tabs.currentIndex()
        names, frags = digests[here]
        self.table.blockSignals(True)   # (the old rows go: no fragment is selected any more)
        self.table.setRowCount(0)
        self.table.blockSignals(False)
        self._clear_fragment()
        shown = frags if names else []   # nothing ticked: the hint, not one "fragment" that is the whole sequence uncut
        self.table.setRowCount(len(shown))
        for i, (size, a, b, left, right) in enumerate(shown):
            ends = f"{'/'.join(left) or 'end'} – {'/'.join(right) or 'end'}"
            for c, v in enumerate([f"{size:,}", f"{a + 1:,}", f"{b:,}", ends]):
                item = QTableWidgetItem(v)
                if c < 3:
                    right_aligned(item)
                self.table.setItem(i, c, item)
        # The selected digest's cuts on the sequence, on both strands.
        ed = self.editor
        tip = lambda name, top, bottom: (lambda: enzyme_mark_tip(ed, name, top, bottom))   # made on hover, not for every cut
        self.editor.set_cut_marks("digest", [(top, bottom, name, tip(name, top, bottom)) for name in names for top, bottom in self.cuts[name]])
        lanes = []
        for k, (n, fr) in enumerate(digests):
            label = "+".join(n) or f"Digest {self.numbers[k]}"
            self.tabs.setTabText(k, label if len(label) < 24 else f"Digest {self.numbers[k]} ({plural(len(n), 'enzyme')})")
            lanes.append((label, [f[0] for f in fr]))   # the gel writes the enzymes one per line
        self.gel.set_lanes(lanes)
        self.frags = frags

    # -- a fragment clicked in the table: selected on the shared views, with the cuts at its two ends pinned
    def show_fragment(self, row):
        """Select fragment `row` of the table on the sequence (so the sequence view, map and genome view all highlight it; a
        fragment across the origin of a circle is one selection that passes it) and pin the cut of each enzyme that makes
        one of its two ends: heavier strokes with the enzyme's name above, whose hover is the enzyme's pill tooltip."""
        if not 0 <= row < len(self.frags):
            return self._clear_fragment()
        size, a, b, left, right = self.frags[row]
        ed = self.editor
        marks = []
        for names, pos in ((left, a), (right, b)):
            for name in names:
                pair = next(((t, bt) for t, bt in self.cuts[name] if t == pos), None)
                if pair:
                    marks.append((pair[0], pair[1], name, self._tip(name, *pair), True))
        end = a + size   # (past the length when the fragment runs through the origin)
        self.fragment_span = None
        ed.select(a, end)
        self.fragment_span = (ed.start, ed.end)
        ed.set_cut_marks("digest_ends", marks)
        ed.reveal.emit(a)

    def _tip(self, name, top, bottom):
        ed = self.editor
        return lambda: enzyme_mark_tip(ed, name, top, bottom)   # made on hover

    def _clear_fragment(self):
        self.fragment_span = None
        self.editor.set_cut_marks("digest_ends", [])

    def _selection_moved(self):
        """Something else was selected (a click on the sequence, another tool): the fragment is no longer the selection, so its
        end cuts go and the table row is deselected."""
        if self.fragment_span is not None and (self.editor.start, self.editor.end) != self.fragment_span:
            self._clear_fragment()
            self.table.blockSignals(True)
            self.table.clearSelection()
            self.table.setCurrentCell(-1, -1)
            self.table.blockSignals(False)

