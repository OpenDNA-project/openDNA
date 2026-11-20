"""The Synthesis tab (a pane's picker > Tools > Synthesis): will a company make this DNA? The sequence is in a sequence view on top (a
scratch sequence of its own, not the open file: type or paste in it); below it, one section for each way of finding out. Hovering a
row of the calculator's table highlights the stretches that measure counted; a click keeps them highlighted and scrolls to the first.
The Synthesis Success Calculator (Halper, Hossain & Salis 2020) runs on this computer and sends nothing anywhere. Each synthesis company's
own check is a section of its own, started by its own button, which says what it sends (and the user is asked once more before it does).

All the figuring is in synthesis_salis.py, idt.py, vienna_install.py and secret_store.py; this file only lays it out and runs it in a worker
thread, so the rest of the plasmid can be worked on meanwhile."""
import threading
from html import escape
from pathlib import Path

from PySide6.QtCore import QEvent, QObject, Qt, Signal
from PySide6.QtGui import QGuiApplication
from PySide6.QtWidgets import (
    QCheckBox, QComboBox, QFormLayout, QGroupBox, QHBoxLayout, QHeaderView, QLabel, QLineEdit, QMessageBox, QProgressBar,
    QPushButton, QSplitter, QTableWidget, QTableWidgetItem, QVBoxLayout, QWidget,
)

from opendna_core import Document

import idt
import online
import secret_store
import synthesis_salis as salis
import theme
import vienna_install
from uikit import no_auto_default, plural
from user_error import UserError

SECONDS_PER_KB = 10   # measured: folding is the slow part (about 10 s per 1,000 bases on a laptop)


class _Relay(QObject):
    """Brings the answer of a worker thread to the GUI thread (lives forever: an answer arriving after the tab closed is harmless)."""
    progress = Signal(object, int, int)   # token, done, total
    done = Signal(object, object, object)   # token, result, error text


_relay = _Relay()


def colored(text, role):
    return f"<span style='color:{theme.hex(role)}'>{escape(text)}</span>"


def duration(seconds):
    return f"{seconds:.0f} s" if seconds < 90 else f"{seconds / 60:.0f} min"


class Section(QGroupBox):
    """A box of the tab with a status line, a progress bar and one background job at a time (a newer job makes an older answer stale)."""

    def __init__(self, title, panel):
        super().__init__(title)
        self.panel, self.token, self.stop, self.buttons = panel, None, threading.Event(), []
        self.status = QLabel("Not run yet.", wordWrap=True, textFormat=Qt.RichText)
        self.verdict = QLabel(textFormat=Qt.RichText, wordWrap=True)
        self.progress = QProgressBar()
        self.progress.hide()
        self.cancel_button = QPushButton("Cancel", clicked=self.cancel)
        self.cancel_button.hide()
        _relay.progress.connect(self._progress)
        _relay.done.connect(self._done)

    def button(self, text, slot):
        b = QPushButton(text, clicked=slot)
        self.buttons.append(b)
        return b

    def busy(self, what, cancellable):
        self.token = token = object()
        self.stop = threading.Event()
        for b in self.buttons:
            b.setEnabled(False)
        self.status.setText(what)
        self.verdict.clear()
        self.progress.setRange(0, 100 if cancellable else 0)
        self.progress.setValue(0)
        self.progress.show()
        self.cancel_button.setVisible(cancellable)
        self.cancel_button.setEnabled(True)
        return token, self.stop

    def spawn(self, token, stop, work):
        def target():
            try:
                result, error = work(token, stop), None
            except salis.Cancelled:
                result, error = "cancelled", None
            except UserError as e:
                result, error = None, str(e)
            except Exception as e:   # a bug must not leave the buttons disabled for good
                result, error = None, f"{type(e).__name__}: {e}"
            try:
                _relay.done.emit(token, result, error)
            except RuntimeError:
                pass   # the program is closing: the relay is gone
        threading.Thread(target=target, daemon=True).start()

    def cancel(self):
        self.stop.set()
        self.cancel_button.setEnabled(False)
        self.status.setText("Cancelling…")

    def _progress(self, token, done, total):
        if token is self.token and total:
            self.progress.setRange(0, total)
            self.progress.setValue(done)

    def _done(self, token, result, error):
        if token is not self.token:
            return   # an answer for an older job, or for the other section
        self.progress.hide()
        self.cancel_button.hide()
        for b in self.buttons:
            b.setEnabled(True)
        if error:
            self.status.setText(colored(error, "status_bad"))
            self.failed()
        elif result == "cancelled":
            self.status.setText("Cancelled.")
        else:
            self.finished(result)

    def failed(self):
        pass

    def finished(self, result):
        raise NotImplementedError


class SalisSection(Section):
    def __init__(self, panel):
        super().__init__(f"Synthesis Success Calculator · {salis.CITATION_SHORT} · runs on this computer", panel)
        self.setToolTip(salis.CITATION + "\n\nNothing is sent anywhere: the model runs here. It reads the sequence as a straight piece.")
        self.run_button = self.button("Run", self.run)
        self.install_button = self.button("Install ViennaRNA…", self.install)
        online.mark(self.install_button, "nothing of your sequence: it downloads the ViennaRNA package (about 2 to 14 MB) from pypi.org")
        self.install_button.setToolTip(self.install_button.toolTip() + "\n\nThe hairpin folding of the calculator needs ViennaRNA (Lorenz et al., University of Vienna). "
                                       "openDNA does not ship it; this fetches the official package into openDNA's data folder.")
        self.install_button.hide()
        self.places, self.pinned = {}, None   # where each measure was counted (from the last run); the row clicked
        # The nine measures side by side, ROWS to a column of (Measured, Value) pairs, so all show at once; what each counts is its tooltip.
        self.table = QTableWidget(self.ROWS, 2 * -(-len(salis.FEATURES) // self.ROWS))
        self.table.setHorizontalHeaderLabels(["Measured", "Value"] * (self.table.columnCount() // 2))
        self.table.verticalHeader().hide()
        self.table.setEditTriggers(QTableWidget.NoEditTriggers)
        self.table.setSelectionMode(QTableWidget.SingleSelection)
        self.table.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        self.table.setMouseTracking(True)
        self.table.cellEntered.connect(lambda row, column: self.panel.highlight(self.spans(self.index(row, column))))
        self.table.cellClicked.connect(lambda row, column: self.pin(self.index(row, column)))
        self.table.viewport().installEventFilter(self)   # leaving the table: back to the measure clicked
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.hide()
        row = QHBoxLayout()
        for w in (self.run_button, self.cancel_button, self.install_button):
            row.addWidget(w)
        row.addWidget(self.progress, 1)
        layout = QVBoxLayout(self)
        layout.addLayout(row)
        layout.addWidget(self.status)
        layout.addWidget(self.verdict)
        layout.addWidget(self.table)
        if not salis.have_vienna():
            self.need_vienna()

    def need_vienna(self):
        self.status.setText(colored("ViennaRNA is not installed: the calculator needs it to fold hairpins.", "status_warn"))
        self.install_button.show()
        self.run_button.setEnabled(False)

    def busy(self, what, cancellable):
        self.table.hide()
        self.forget_places()
        return super().busy(what, cancellable)

    ROWS = 3   # measures per column of the table

    def index(self, row, column):
        """The measure (its place in salis.FEATURES) shown in that cell, or -1 for an empty cell."""
        i = column // 2 * self.ROWS + row
        return i if i < len(salis.FEATURES) else -1

    def spans(self, i):
        """The stretches measure `i` counted, [] when it points at nothing (or the sequence changed since)."""
        return self.places.get(salis.FEATURES[i], []) if 0 <= i < len(salis.FEATURES) else []

    def pin(self, i):
        if i < 0:
            return
        self.pinned = i
        spans = self.spans(i)
        self.panel.highlight(spans, reveal=True)
        label = salis.FEATURE_HELP[salis.FEATURES[i]][0]
        self.panel.info.setText(f"{label}: {plural(len(spans), 'stretch', 'stretches')} highlighted" if spans else f"{label}: concerns the whole sequence")

    def forget_places(self):
        """The sequence changed (or a new run starts): the stretches no longer fit it."""
        self.places, self.pinned = {}, None
        self.table.clearSelection()
        self.panel.highlight([])

    def eventFilter(self, watched, event):
        if event.type() == QEvent.Leave:
            self.panel.highlight(self.spans(self.pinned) if self.pinned is not None else [])
        return False

    def run(self):
        try:
            seq = self.panel.sequence()
        except UserError as e:
            self.status.setText(colored(str(e), "status_bad"))
            return
        estimate = duration(len(seq) / 1000 * SECONDS_PER_KB)
        token, stop = self.busy(f"Folding and counting {plural(len(seq), 'base')} (about {estimate}; you can keep working meanwhile)…", True)
        self.spawn(token, stop, lambda t, s: {"mode": "run", **salis.evaluate(
            seq, progress=lambda done, total: _relay.progress.emit(t, done, total), cancelled=s.is_set)})

    def install(self):
        token, stop = self.busy("Downloading ViennaRNA from pypi.org…", False)
        self.spawn(token, stop, lambda t, s: {"mode": "install", "version": vienna_install.install()})

    def failed(self):
        if not salis.have_vienna():
            self.need_vienna()

    def finished(self, result):
        if result["mode"] == "install":
            self.install_button.hide()
            self.run_button.setEnabled(True)
            self.status.setText(colored(f"ViennaRNA {result['version']} installed.", "status_ok"))
            return
        self.places = result["places"]
        quick = result["outcome"]
        self.status.setText(f"{plural(result['features']['size'], 'base')} tested. The model's prediction, not a promise: a company's own check decides.")
        role = "status_ok" if quick else "status_bad"
        words = f"{'Likely' if quick else 'Unlikely'} to be synthesised quickly: probability {result['probability']:.0%}"
        self.verdict.setText(f"<span style='font-size:15px;font-weight:600;color:{theme.hex(role)}'>{escape(words)}</span>")
        for i, name in enumerate(salis.FEATURES):
            label, meaning = salis.FEATURE_HELP[name]
            value = result["features"][name]
            text = f"{value:,}" if isinstance(value, int) else f"{value:.3f}"
            tip = f"{label}: {meaning}. Point at it to highlight where it was counted in the sequence above; click to keep the highlight."
            for k, cell in enumerate((label, text)):
                item = QTableWidgetItem(cell)
                item.setToolTip(tip)
                if k == 1:
                    item.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
                self.table.setItem(i % self.ROWS, 2 * (i // self.ROWS) + k, item)
        self.table.resizeColumnsToContents()
        self.table.resizeRowsToContents()
        # Exactly as tall as its rows: nothing to scroll.
        self.table.setFixedHeight(self.table.horizontalHeader().height() + sum(self.table.rowHeight(r) for r in range(self.ROWS))
                                  + 2 * self.table.frameWidth())
        self.table.show()


class IdtSection(Section):
    """IDT's own complexity screen: needs the user's IDT API account; sends the sequence and the login to idtdna.com."""
    KEY = "idt"   # the name the login is kept under in secret_store

    def __init__(self, panel):
        super().__init__("Integrated DNA Technologies (IDT) · complexity screen · sends your sequence to idtdna.com", panel)
        self.raw = ""
        self.setToolTip("IDT's own check of a gBlocks / gene sequence. It needs your own IDT API account (a username and password, and the client id and "
                        f"client secret IDT gives for API use: {idt.SIGN_UP_URL}). Your sequence and these details go to idtdna.com, and only "
                        "when you press the button and confirm.")
        self.fields = {}
        form = QFormLayout()
        for key, label, hide in (("username", "IDT username", False), ("password", "IDT password", True),
                                 ("client_id", "API client id", False), ("client_secret", "API client secret", True)):
            edit = self.fields[key] = QLineEdit(echoMode=QLineEdit.Password if hide else QLineEdit.Normal)
            form.addRow(label, edit)
        self.product = QComboBox()
        for name, cutoff in idt.PRODUCTS:
            self.product.addItem(f"{name} (score under {cutoff:g})", cutoff)
        form.addRow("Judge the score for", self.product)
        store = secret_store.backend()
        self.remember = QCheckBox(f"Remember these in {store}" if store else "Keep these until openDNA closes (this system has no secure storage)")
        saved = secret_store.load(self.KEY)
        self.remember.setChecked(bool(store) and bool(saved))
        self.remember.setToolTip("Passwords are never written to a file: they go to the system's own store (Windows Credential Manager, macOS Keychain, "
                                 "Linux Secret Service), or stay in memory.")
        self.run_button = self.button("Check with IDT", self.run)
        online.mark(self.run_button, "your sequence (as DNA text) and your IDT username, password, client id and client secret, to idtdna.com")
        self.forget_button = QPushButton("Forget my login", clicked=self.forget)
        self.copy_button = QPushButton("Copy IDT's answer as JSON", clicked=lambda: QGuiApplication.clipboard().setText(self.raw))
        self.copy_button.hide()
        self.table = QTableWidget()
        self.table.setEditTriggers(QTableWidget.NoEditTriggers)
        self.table.setMaximumHeight(200)
        self.table.hide()
        note = QLabel(f"You need your own IDT API account (<a href='{idt.SIGN_UP_URL}'>how to get one</a>). Nothing is sent until you press the button "
                      "and confirm.", openExternalLinks=True, wordWrap=True)
        row = QHBoxLayout()
        for w in (self.run_button, self.cancel_button, self.forget_button, self.copy_button):
            row.addWidget(w)
        row.addWidget(self.progress, 1)
        layout = QVBoxLayout(self)
        layout.addWidget(note)
        layout.addLayout(form)
        layout.addWidget(self.remember)
        layout.addLayout(row)
        layout.addWidget(self.status)
        layout.addWidget(self.verdict)
        layout.addWidget(self.table)
        for key, edit in self.fields.items():
            edit.setText(saved.get(key, ""))

    def credentials(self):
        values = {k: e.text().strip() for k, e in self.fields.items()}
        missing = [label for k, label in (("username", "the username"), ("password", "the password"), ("client_id", "the client id"),
                                          ("client_secret", "the client secret")) if not values[k]]
        if missing:
            raise UserError("Fill in " + ", ".join(missing) + ".")
        return values

    def forget(self):
        secret_store.forget(self.KEY)
        for edit in self.fields.values():
            edit.clear()
        self.remember.setChecked(False)
        self.status.setText("Your IDT login is forgotten.")

    def run(self):
        try:
            seq = self.panel.sequence()
            idt.check_length(seq)
            values = self.credentials()
        except UserError as e:
            self.status.setText(colored(str(e), "status_bad"))
            return
        answer = QMessageBox.question(self, "Send to IDT?", f"This sends the {len(seq):,} bases in the box, and your IDT login, to idtdna.com. Send?",
                                      QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
        if answer != QMessageBox.Yes:
            self.status.setText("Nothing was sent.")
            return
        if self.remember.isChecked():
            if not secret_store.save(self.KEY, values) and secret_store.secure():
                self.status.setText(colored("Could not store the login securely; it is kept until openDNA closes.", "status_warn"))
        else:
            secret_store.forget(self.KEY)
        token, stop = self.busy("Asking IDT…", False)
        self.copy_button.hide()
        self.table.hide()
        self.spawn(token, stop, lambda t, s: idt.screen(seq, idt.Credentials(**values)))

    def finished(self, result):
        self.raw, rows = result["raw"], result["rows"]
        score, (name, cutoff) = result["score"], idt.PRODUCTS[self.product.currentIndex()]
        if score is None:
            self.status.setText(colored("IDT's answer has no \"Score\" to judge. Its full answer is in the table and the JSON.", "status_warn"))
        else:
            ok = score < cutoff
            words = f"IDT complexity score {score:g}: {'below' if ok else 'at or above'} {cutoff:g}, the limit quoted for {name}"
            self.status.setText("IDT answered. Its own ordering page decides: this is the number it returned, set against the limit quoted for the product.")
            self.verdict.setText(f"<span style='font-size:15px;font-weight:600;color:{theme.hex('status_ok' if ok else 'status_bad')}'>{escape(words)}</span>")
        columns = sorted({k for r in rows for k in r}, key=lambda k: (k != "Name", k))
        self.table.clear()
        self.table.setColumnCount(len(columns))
        self.table.setRowCount(len(rows))
        self.table.setHorizontalHeaderLabels(columns)
        for i, r in enumerate(rows):
            for j, k in enumerate(columns):
                self.table.setItem(i, j, QTableWidgetItem("" if k not in r else str(r[k])))
        self.table.resizeColumnsToContents()
        self.table.setVisible(bool(rows))
        self.copy_button.show()


class SynthesisPanel(QWidget):
    def __init__(self, window):
        super().__init__()
        self.window = window
        from editor import Editor
        from sequence_view import SequenceView
        intro = QLabel("Type or paste a DNA sequence (bases, or one FASTA record) in the sequence view to see how likely a synthesis company is to "
                       "make it quickly. Select bases in a file and choose Test Synthesis… on right-click to fill it in. It is a scratch sequence: "
                       "your files are not changed.", wordWrap=True)
        # (a document of its own, linear, never saved: the calculator reads a straight piece)
        self.editor = Editor(Document.new_sequence("Synthesis", "", False), Path("Synthesis.gb"))
        self.view = SequenceView(self.editor)
        self.editor.changed.connect(self.edited)
        self.info = QLabel(textFormat=Qt.RichText)
        self.paste_button = QPushButton("Paste", clicked=lambda: self.set_text(QGuiApplication.clipboard().text()),
                                        toolTip="Replace the sequence with the clipboard's bases (one FASTA record, or bases)")
        self.selection_button = QPushButton("Use selection", clicked=lambda: self.use(True))
        self.whole_button = QPushButton("Use whole sequence", clicked=lambda: self.use(False))
        self.clear_button = QPushButton("Clear", clicked=lambda: self.set_sequence(""))
        buttons = QHBoxLayout()
        for b in (self.paste_button, self.selection_button, self.whole_button, self.clear_button):
            buttons.addWidget(b)
        buttons.addWidget(self.info, 1)
        self.salis = SalisSection(self)
        self.idt = IdtSection(self)
        top = QWidget()
        top_box = QVBoxLayout(top)
        top_box.setContentsMargins(0, 0, 0, 0)
        top_box.addWidget(intro)
        top_box.addWidget(self.view, 1)
        top_box.addLayout(buttons)
        sections = QWidget()
        below = QVBoxLayout(sections)
        below.setContentsMargins(0, 0, 0, 0)
        below.addWidget(self.salis)
        below.addWidget(self.idt)   # (one section for each company, below the local one)
        below.addStretch(1)
        split = QSplitter(Qt.Vertical)
        split.setChildrenCollapsible(False)
        split.addWidget(top)
        split.addWidget(sections)
        split.setSizes([10000, 15000])
        layout = QVBoxLayout(self)
        layout.addWidget(split)
        no_auto_default(self)
        self.count()

    def sequence(self):
        """The bases of the sequence view, checked (UserError names the first base the calculator cannot read)."""
        return salis.sequence_from_text(self.editor.doc.sequence())

    def edited(self):
        self.salis.forget_places()   # (the stretches of the last run no longer fit)
        self.count()

    def count(self):
        try:
            self.info.setText(plural(len(self.sequence()), "base"))
        except UserError as e:
            self.info.setText(colored(str(e), "status_bad") if self.editor.length else "")

    def highlight(self, spans, reveal=False):
        """Highlight `spans` [(start, end)] in the sequence view (the search marks of this scratch sequence); `reveal`: scroll to the first."""
        self.editor.set_search_hits(spans)
        if reveal and spans:
            self.editor.reveal.emit(spans[0][0])

    def set_text(self, text):
        """Replace the sequence with the bases of a pasted text (one FASTA record, or bases); a message when it has none."""
        try:
            self.set_sequence(salis.sequence_from_text(text))
        except UserError as e:
            self.info.setText(colored(str(e), "status_bad"))

    def use(self, selection):
        editor = getattr(self.window, "editor", None)
        if editor is None:
            self.info.setText(colored("No sequence is open.", "status_warn"))
        elif selection and editor.end <= editor.start:
            self.info.setText(colored("Nothing is selected.", "status_warn"))
        else:
            self.set_sequence(editor.selected_sequence() if selection else editor.doc.sequence(0, None))

    def set_sequence(self, seq):
        """Make `seq` the sequence of the view (nothing runs: the user decides)."""
        ed = self.editor
        ed.select(0, ed.length)
        if seq:
            ed.paste(seq)
        elif ed.length:
            ed.backspace()
