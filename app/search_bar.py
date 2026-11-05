"""The search bar of the main toolbar: type what to find, the views highlight every hit, the counter says "3 of 17", Enter / Shift+Enter
(or the arrows) go to the next / previous one. The Options menu says what is looked for (see search.py for what each option does).

A search runs by itself a moment after typing when the sequence is small (LIVE_MAX); on a bigger one, or across a whole FASTA file, it
runs on Enter (the file in a worker thread, the hits streaming in per sequence).
"""
import threading
from dataclasses import replace

from PySide6.QtCore import QEvent, QObject, QSettings, Qt, QTimer, Signal
from PySide6.QtGui import QPalette
from PySide6.QtWidgets import (
    QButtonGroup, QCheckBox, QHBoxLayout, QLabel, QLineEdit, QMenu, QRadioButton, QSpinBox, QToolButton, QVBoxLayout, QWidget,
    QWidgetAction,
)

import icons
import search
from primer_db import data_dir
from search import Options

MAX_QUERY = 5000   # characters: a longer text is not something to look for (and would stall the live search)
LIVE_MAX = 2_000_000   # bp: below this the hits follow the typing; above it the search waits for Enter
LIVE_MS = 350
SCOPE_TEXT = {"sequence": "This sequence", "selection": "Selection only", "file": "All sequences in this file"}


class _Relay(QObject):
    """Brings the hits of the file search from its thread to the GUI thread (lives forever: a result arriving late is harmless)."""
    found = Signal(object, object, object)   # token, sequence name, [hit]
    done = Signal(object)


_relay = _Relay()


class _Counter(QLabel):
    """The hit counter, inside the right end of the search box ("3 of 12", like a browser's find bar): a fixed width, so its
    text changing never moves anything; a longer text ("Searching 3 / 24…") is cut short there, whole in its tooltip."""

    def __init__(self, box):
        super().__init__(box)
        self.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
        self.setForegroundRole(QPalette.PlaceholderText)   # quieter than what is typed
        self.setFixedWidth(self.fontMetrics().horizontalAdvance("9,999 of 9,999") + 6)
        self.full = ""

    def setText(self, text):
        self.full = text
        box = self.parentWidget()   # the box keeps the counter's room only while it says something: the placeholder is whole
        box.setTextMargins(0, 0, self.width() + 6 if text else 0, 0)
        self.setToolTip(text if self.fontMetrics().horizontalAdvance(text) > self.width() else "")
        super().setText(self.fontMetrics().elidedText(text, Qt.ElideRight, self.width()))

    def text(self):
        return self.full


class _Box(QLineEdit):
    backward = Signal()   # Shift+Enter
    escaped = Signal()
    counter = None

    def resizeEvent(self, e):
        super().resizeEvent(e)
        if self.counter is not None:   # just left of the clear button
            room = self.height()
            self.counter.setGeometry(self.width() - room - 6 - self.counter.width(), 0, self.counter.width(), self.height())

    def keyPressEvent(self, e):
        if e.key() in (Qt.Key_Return, Qt.Key_Enter) and e.modifiers() & Qt.ShiftModifier:
            self.backward.emit()
        elif e.key() == Qt.Key_Escape:
            self.escaped.emit()
        else:
            super().keyPressEvent(e)


class SearchBar(QWidget):
    def __init__(self, window):
        super().__init__()
        self.window, self.editor = window, None
        self.hits, self.index, self.key = [], -1, None   # key: (query, options) the hits are for; None = no valid search
        self.rank, self.file_path = {}, None             # file search: the order of the sequences; the file the hits are from
        self.token, self.stop_flag, self.want, self.total, self.seen = None, threading.Event(), 0, 0, 0
        self.settings = QSettings(str(data_dir() / "search.ini"), QSettings.IniFormat)
        self.timer = QTimer(self, singleShot=True, interval=LIVE_MS, timeout=self._live)

        self.box = _Box(placeholderText="Find sequence or name…", clearButtonEnabled=True)   # (the rest: tooltip)
        self.box.setMaxLength(MAX_QUERY)
        self.box.setMinimumWidth(190)
        self.box.setMaximumWidth(280)
        self.box.setToolTip("DNA (both strands, IUPAC codes), a protein motif, or the name of a feature, primer or enzyme.\n"
                            "Enter: next hit · Shift+Enter: previous · Esc: clear the highlights. What is looked for: Options.")
        # Inside the box (see _Counter): its own fixed width, so its text going from "" to "12 hits · 3 / 24…" while typing
        # never moves anything, and the toolbar keeps the room a counter beside the box took.
        self.counter = self.box.counter = _Counter(self.box)
        self.box.setTextMargins(0, 0, 0, 0)   # room for the counter only while it says something (see _Counter.setText)
        self.back = QToolButton(text="◀", toolTip="Previous hit (Shift+Enter)")
        self.forward = QToolButton(text="▶", toolTip="Next hit (Enter)")
        self.menu = QMenu(self)
        action = QWidgetAction(self.menu)
        action.setDefaultWidget(self._panel())
        self.menu.addAction(action)
        # a cog (what the search looks for: DNA, protein, names...), as compact as the arrows beside it
        self.options_button = QToolButton(popupMode=QToolButton.InstantPopup, autoRaise=True,
                                          toolTip="Search options: what the search looks for")
        self.options_button.setIcon(icons.pane_icon("options", self.palette()))
        self.options_button.setMenu(self.menu)
        self.options_button.setStyleSheet("QToolButton::menu-indicator { image: none; width: 0; }")   # (the cog alone)

        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(3)
        for w in (self.box, self.options_button, self.back, self.forward):
            layout.addWidget(w)
        layout.addStretch()   # room the toolbar gives beyond what these need goes after them, not between them
        self.box.returnPressed.connect(lambda: self.step(1))
        self.box.backward.connect(lambda: self.step(-1))
        self.box.escaped.connect(self.escape)
        self.box.textChanged.connect(self._typed)
        self.forward.clicked.connect(lambda: self.step(1))
        self.back.clicked.connect(lambda: self.step(-1))
        _relay.found.connect(self._found)
        _relay.done.connect(self._done)

    # ---- the options ------------------------------------------------------------------------------------------------------------
    def _panel(self):
        panel = QWidget()
        lay = QVBoxLayout(panel)
        lay.setContentsMargins(10, 8, 10, 8)

        def check(text, tip, key, default):
            box = QCheckBox(text, toolTip=tip)
            box.setChecked(self.settings.value(key, default, type=bool))
            box.toggled.connect(self._options_changed)
            lay.addWidget(box)
            return box
        self.both = check("Both strands", "Also find the reverse complement of DNA.", "both", True)
        self.iupac = check("IUPAC codes", "R, Y, S, W, K, M, B, D, H, V and N in the query stand for several bases.", "iupac", True)
        row = QHBoxLayout()
        self.fuzzy_on = check("Fuzzy: up to", "Allow mismatches (never more than a third of the query).", "fuzzy_on", False)
        lay.removeWidget(self.fuzzy_on)
        self.fuzzy_n = QSpinBox(minimum=1, maximum=10, value=self.settings.value("fuzzy_n", 1, type=int), suffix=" mismatches")
        self.fuzzy_n.valueChanged.connect(self._options_changed)
        row.addWidget(self.fuzzy_on)
        row.addWidget(self.fuzzy_n)
        lay.addLayout(row)
        self.translations = check("Translations", "A protein motif (one-letter amino acids, X = any) in the 6 reading frames.",
                                  "translations", False)
        self.details = check("Feature details", "Look in every qualifier of the features (notes, references...), not only their names.",
                             "details", False)
        self.enzymes = check("Restriction sites", "An enzyme's name (EcoRI) finds its recognition sites.", "enzymes", False)
        self.primers = check("Primers", "A primer's name finds where it binds.", "primers", True)
        lay.addWidget(QLabel("<b>Where</b>"))
        self.scopes, group = {}, QButtonGroup(panel)
        saved = self.settings.value("scope", "sequence")
        for scope, text in SCOPE_TEXT.items():
            radio = QRadioButton(text, checked=scope == saved)
            radio.toggled.connect(self._options_changed)
            group.addButton(radio)
            lay.addWidget(radio)
            self.scopes[scope] = radio
        if not any(r.isChecked() for r in self.scopes.values()):
            self.scopes["sequence"].setChecked(True)
        return panel

    def options(self):
        return Options(both=self.both.isChecked(), iupac=self.iupac.isChecked(),
                       fuzzy=self.fuzzy_n.value() if self.fuzzy_on.isChecked() else 0, translations=self.translations.isChecked(),
                       details=self.details.isChecked(), enzymes=self.enzymes.isChecked(), primers=self.primers.isChecked(),
                       scope=next(s for s, r in self.scopes.items() if r.isChecked()))

    def _options_changed(self, *_):
        o = self.options()
        for key, value in (("both", o.both), ("iupac", o.iupac), ("fuzzy_on", o.fuzzy > 0), ("fuzzy_n", self.fuzzy_n.value()),
                           ("translations", o.translations), ("details", o.details), ("enzymes", o.enzymes), ("primers", o.primers),
                           ("scope", o.scope)):
            self.settings.setValue(key, value)
        self._typed()

    # ---- the editor ---------------------------------------------------------------------------------------------------------------
    def changeEvent(self, event):
        if event.type() == QEvent.PaletteChange:   # the cog in the theme's text colour
            self.options_button.setIcon(icons.pane_icon("options", self.palette()))
        super().changeEvent(event)

    def attach(self, ed):
        """The window shows a new document. Hits of a file search stay when the new one is another sequence of the same file."""
        same_file = self.file_path is not None and ed.fasta_record and str(ed.path) == self.file_path
        self.editor = ed
        ed.changed.connect(lambda: ed is self.editor and self.clear())   # the hits would be stale
        ed.selection_changed.connect(lambda: ed is self.editor and self._show_count())
        if same_file:
            self._mark()
            self._show_count()
        else:
            self.clear()

    def focus(self):
        self.box.setFocus()
        self.box.selectAll()

    def escape(self):
        self.clear()
        self.window.sequence_view.setFocus()

    # ---- searching ----------------------------------------------------------------------------------------------------------------
    def _typed(self, *_):
        if self.editor is None:
            return
        if self.box.text().strip():
            self.timer.start()
        else:
            self.clear()

    def _live(self):
        ed, opts = self.editor, self.options()
        if ed is None:
            return
        if opts.scope != "file" and ed.length <= LIVE_MAX:
            self.run(0)
        else:
            self.clear()
            self.counter.setText("Press Enter")

    def step(self, direction):
        """Go to the next (1) or previous (-1) hit; searches first if the text or the options changed since the last search."""
        if self.editor is None or not self.box.text().strip():
            return
        self.timer.stop()
        if self.key != (self.box.text().strip(), self.options().key()):
            self.run(direction)
        elif self.hits:
            self._go(self._pick(direction))
        else:
            self.window.statusBar().showMessage(f"'{self.box.text().strip()}' not found", 3000)

    def run(self, direction=1):
        """Search now; then go to the next (1) / previous (-1) hit, or stay where we are (0)."""
        ed, query, opts = self.editor, self.box.text().strip(), self.options()
        self.stop_flag.set()
        self.token = None   # results of a file search still on their way are ignored
        self.timer.stop()
        if not query:
            self.clear()
            return
        note = ""
        if ed.multi_genbank and opts.scope != "selection":   # one GenBank file, several chromosomes: the search covers them all
            self._run_file(query, replace(opts, scope="file"), direction)
            return
        if opts.scope == "file" and ed.fasta_record and search.dna_query(query, opts):
            self._run_file(query, opts, direction)
            return
        if opts.scope == "file":   # names and proteins belong to the open sequence
            opts, note = replace(opts, scope="sequence"), ("Searched the open sequence only (a file search looks for DNA)." if ed.fasta_record else "")
        hits, more = search.search_sequence(ed, query, opts)
        self.file_path, self.rank = None, {}
        self.hits, self.index = hits, -1
        self.key = (query, self.options().key())
        self._mark()
        self._show_count()
        if note or more:
            self.window.statusBar().showMessage(note or more, 4000)
        if direction and hits:
            self._go(self._pick(direction))
        elif direction:
            self.window.statusBar().showMessage(f"'{query}' not found", 3000)

    def _run_file(self, query, opts, direction):
        ed = self.editor
        picker = self.window.picker
        names = [picker.itemText(i) for i in range(picker.count())] or [ed.fasta_record]
        self.rank = {n: i for i, n in enumerate(names)}
        self.file_path, self.hits, self.index = str(ed.path), [], -1
        self.key = (query, self.options().key())
        flag, token = threading.Event(), object()
        self.stop_flag, self.token = flag, token
        self.total, self.seen, self.want = len(names), 0, direction
        path, current, bases = ed.path, ed.fasta_record, ed.doc.sequence()   # a copy: the document may change meanwhile
        others = names
        if ed.multi_genbank:   # the open record is searched in full here; the other records in the thread
            own, _ = search.search_sequence(ed, query, replace(opts, scope="sequence"))
            if opts.primers:   # primer binding sites on the other records (found when the file was opened)
                needle = query.strip().lower()
                own += [search.hit(s.start, s.end, "primer", s.label, s.strand, s.mismatches, record=s.record)
                        for s in ed.elsewhere if needle in s.label.lower()]
            for h in own:
                h.record = h.record or current
            others = [n for n in names if n != current]
            self._found(token, current, own)
        self._mark()
        self.counter.setText(f"Searching {self.seen} / {self.total}…")

        def work():
            try:
                search.search_file(path, others, current, bases, query, opts, flag, lambda n, h: _relay.found.emit(token, n, h))
                _relay.done.emit(token)
            except RuntimeError:
                pass   # the program is closing: the relay is already gone
        threading.Thread(target=work, daemon=True).start()

    def _found(self, token, name, hits):
        if token is not self.token:
            return
        self.seen += 1
        self.hits = (self.hits + sorted(hits, key=lambda h: h.start))[:search.MAX_HITS]
        if name == self.editor.fasta_record:
            self._mark()
        self.counter.setText(f"{len(self.hits):,} hits · {self.seen} / {self.total}…")

    def _done(self, token):
        if token is not self.token:
            return
        self._show_count()
        if self.want and self.hits:
            self._go(self._pick(self.want))
        elif self.want:
            self.window.statusBar().showMessage(f"'{self.box.text().strip()}' not found in this file", 3000)
        self.want = 0

    # ---- going to a hit -----------------------------------------------------------------------------------------------------------
    def _pick(self, direction):
        """Index of the next (1) or previous (-1) hit from the cursor, wrapping round."""
        ed = self.editor
        here, strict = (self.rank.get(ed.fasta_record, 0), ed.start), ed.end > ed.start   # a hit under a bare cursor counts
        keys = [(self.rank.get(h.record, 0), h.start) for h in self.hits]
        if direction > 0:
            return next((i for i, k in enumerate(keys) if (k > here if strict else k >= here)), 0)
        return next((i for i in range(len(keys) - 1, -1, -1) if keys[i] < here), len(keys) - 1)

    def _go(self, i):
        h = self.hits[i]
        self.index = i
        ed = self.editor
        if h.record and h.record != ed.fasta_record:
            self.window.go_to_hit(h.record, h.start, h.end)   # opens that sequence: attach() runs for the new editor
            ed = self.editor
        elif h.feature is not None:
            ed.select_feature(h.feature)
        else:
            ed.select(h.start, h.end)
            ed.reveal.emit(h.start)
        self._show_count()
        self.window.statusBar().showMessage(search.describe(h), 5000)

    # ---- showing the result -------------------------------------------------------------------------------------------------------
    def _mark(self):
        """Tell the open editor which hits to highlight (those of its own sequence)."""
        ed = self.editor
        # a hit found on one strand is highlighted on that strand only (features and tracks cover both)
        own = [(h.start, h.end, h.strand if h.kind in ("dna", "protein", "primer", "enzyme") else 0)
               for h in self.hits if not self.file_path or h.record == ed.fasta_record]
        ed.set_search_hits(own)

    def _selected(self, i):
        ed, h = self.editor, self.hits[i] if 0 <= i < len(self.hits) else None
        return h is not None and (h.start, h.end) == (ed.start, ed.end) and (not h.record or h.record == ed.fasta_record)

    def _show_count(self):
        n = len(self.hits)
        if self.key is None:
            self.counter.setText("")
            return
        if not n:
            self.counter.setText("No hits")
            return
        i = self.index if self._selected(self.index) else next((i for i in range(n) if self._selected(i)), -1)
        total = f"{n:,}+" if n >= search.MAX_HITS else f"{n:,}"
        self.counter.setText(f"{i + 1:,} of {total}" if i >= 0 else f"{total} hit{'s' if n > 1 else ''}")

    def clear(self):
        self.stop_flag.set()
        self.token = None
        self.hits, self.index, self.key, self.file_path, self.rank = [], -1, None, None, {}
        self.counter.setText("")
        if self.editor is not None:
            self.editor.set_search_hits([])
