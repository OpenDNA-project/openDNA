"""The Terminal view (picked in a pane like Sequence or Map): pipelines on the left (forms built from tools.toml), a persistent bash session on
the right (in WSL on Windows) where every external program runs. A pipeline's Run types its command
into the session. `$REF` is the document as FASTA (rewritten before each command), `$DOC` its file;
`open <file>` brings a result back (runner.open_kind decides where). One session per Terminal pane.
See docs/adr/0001-terminal-runs-every-external-tool.md.
"""
import html
from functools import partial
import os
import re
import shlex
import tempfile
from pathlib import Path
from types import SimpleNamespace

from PySide6.QtCore import QPointF, QProcess, QRectF, QSettings, Qt, QUrl
from PySide6.QtGui import QDesktopServices, QPainter, QPen, QPolygonF
from PySide6.QtWidgets import (
    QCheckBox, QComboBox, QFileDialog, QFormLayout, QGroupBox, QHBoxLayout, QLabel, QLineEdit, QListWidget, QMessageBox,
    QPushButton, QSpinBox, QSplitter, QStackedWidget, QTextBrowser, QVBoxLayout, QWidget,
)

import msa
import multigb
import runner
from add_record import document_records, files_records
from editor import mono_font
from msa_view import MsaView
from primer_db import data_dir
from uikit import run_busy

HISTORY_KEY, HISTORY_MAX = "terminal/history", 500
CUSTOM = "Custom"
READS = {k: runner.INPUT_NAMES[k] + (" (optional)" if k == "reads2" else "") for k in ("reads1", "reads2", "long")}
PATH_RE = re.compile(r"(/mnt/[a-z]/[^\s'\"]+|[A-Za-z]:\\[^\s'\"]+)")


def write_fasta(path, seqs):
    """FASTA with names safe for the tools (and unique: two records of one name are told apart by a number)."""
    names = []
    for name, _ in seqs:
        names.append(multigb.record_name(runner.safe_name(name), names))
    Path(path).write_text(msa.to_fasta(names, [s for _, s in seqs]), encoding="utf-8")


class PipelineDiagram(QWidget):
    """The current pipeline as a flow chart (runner.diagram): its inputs side by side on top, then one box per
    step, down to the result. Redrawn whenever the preset, an option or the inputs shown change."""
    BOX, GAP, PAD, WIDTH = 22, 14, 6, 280

    def __init__(self):
        super().__init__()
        self.nodes = []
        self.setToolTip("The steps the selected pipeline runs, from its inputs to the result opened in openDNA. "
                        "A dashed step is switched off by the current options.")

    def set_nodes(self, nodes):
        self.nodes = nodes
        rows = any(k == "input" for k, _ in nodes) + sum(k != "input" for k, _ in nodes)
        self.setFixedHeight(2 * self.PAD + rows * self.BOX + max(rows - 1, 0) * self.GAP)
        self.update()

    def paintEvent(self, _event):
        p, pal = QPainter(self), self.palette()
        p.setRenderHint(QPainter.Antialiasing)
        fm, w = p.fontMetrics(), self.width()
        inputs = [t for k, t in self.nodes if k == "input"]
        above, y = [], self.PAD
        if inputs:
            slot = (w - 2 * self.PAD) / len(inputs)
            for i, text in enumerate(inputs):
                r = QRectF(self.PAD + i * slot + 2, y, slot - 4, self.BOX)
                self._box(p, r, text, pal.alternateBase().color(), QPen(pal.mid().color(), 1), pal.text().color(), fm, round_=True)
                above.append(r.center().x())
            y += self.BOX + self.GAP
        bw = min(w - 2 * self.PAD, self.WIDTH)
        x = (w - bw) / 2
        for kind, text in self.nodes:
            if kind == "input":
                continue
            r = QRectF(x, y, bw, self.BOX)
            for cx in above:
                self._arrow(p, QPointF(cx, y - self.GAP), QPointF(r.center().x(), y), pal.placeholderText().color())
            if kind == "output":
                self._box(p, r, text, pal.highlight().color(), QPen(pal.highlight().color(), 1), pal.highlightedText().color(), fm, round_=True)
            elif kind == "skipped":
                self._box(p, r, text + " (off)", pal.window().color(), QPen(pal.mid().color(), 1, Qt.DashLine),
                          pal.placeholderText().color(), fm)
            else:
                self._box(p, r, text, pal.base().color(), QPen(pal.highlight().color(), 1.5), pal.text().color(), fm)
            above, y = [r.center().x()], y + self.BOX + self.GAP

    @staticmethod
    def _box(p, r, text, fill, pen, color, fm, round_=False):
        p.setPen(pen)
        p.setBrush(fill)
        radius = r.height() / 2 if round_ else 4
        p.drawRoundedRect(r, radius, radius)
        p.setPen(color)
        p.drawText(r.adjusted(6, 0, -6, 0), Qt.AlignCenter, fm.elidedText(text, Qt.ElideRight, int(r.width()) - 12))

    @staticmethod
    def _arrow(p, a, b, color):
        p.setPen(QPen(color, 1.2))
        p.drawLine(a, QPointF(b.x(), b.y() - 4) if a.x() == b.x() else b)
        p.setBrush(color)
        p.drawPolygon(QPolygonF([b, QPointF(b.x() - 4, b.y() - 6), QPointF(b.x() + 4, b.y() - 6)]))


class HistoryLine(QLineEdit):
    """The command line: Up / Down walk through the commands typed before (kept across restarts)."""

    def __init__(self):
        super().__init__(placeholderText="bash command ($REF = this sequence as FASTA, $DOC = its file, open <file> to view a result)")
        self.settings = QSettings(str(data_dir() / "window.ini"), QSettings.IniFormat)
        saved = self.settings.value(HISTORY_KEY, []) or []
        self.history = [saved] if isinstance(saved, str) else list(saved)   # an .ini gives a one-item list back as a str
        self.at = len(self.history)

    def remember(self, command):
        if not self.history or self.history[-1] != command:
            self.history = (self.history + [command])[-HISTORY_MAX:]
            self.settings.setValue(HISTORY_KEY, self.history)
        self.at = len(self.history)

    def keyPressEvent(self, event):
        step = {Qt.Key_Up: -1, Qt.Key_Down: 1}.get(event.key())
        if step is None or not self.history:
            return super().keyPressEvent(event)
        self.at = max(0, min(len(self.history), self.at + step))
        self.setText(self.history[self.at] if self.at < len(self.history) else "")


class PipelineForm(QWidget):
    """A pipeline's form: preset, inputs, output, options folded under Advanced, Run and Install."""

    def __init__(self, terminal, key, pipeline):
        super().__init__()
        self.terminal, self.key, self.pipeline = terminal, key, pipeline
        form = QFormLayout(self)
        self.preset = QComboBox()
        self.preset.addItems(list(pipeline.get("presets", {})))
        self.preset.activated.connect(lambda: self.apply_preset(self.preset.currentText()))
        form.addRow("Preset", self.preset)
        self.rows, self.edits = {}, {}
        for name in pipeline["inputs"]:
            if name in READS:
                self.edits[name] = edit = QLineEdit()
                self.rows[name] = self._browse_row(edit, lambda _=False, e=edit, n=name: self._pick_reads(e, n))
                form.addRow(READS[name], self.rows[name])
            elif name == "reference":
                self.reference = QComboBox()
                self.reference.activated.connect(self._reference_picked)
                self.rows[name] = self.reference
                form.addRow("Reference", self.reference)
            elif name == "sequences":
                self.own = QCheckBox("This document's records", checked=True)
                self.extra = QLineEdit(placeholderText="more FASTA / GenBank files (optional)")
                self.rows[name] = box = QWidget()
                v = QVBoxLayout(box)
                v.setContentsMargins(0, 0, 0, 0)
                v.addWidget(self.own)
                v.addWidget(self._browse_row(self.extra, self._pick_extra))
                form.addRow("Sequences", box)
        self.out = QLineEdit(placeholderText="chosen from the inputs when left empty")
        form.addRow("Output", self._browse_row(self.out, self._pick_out))
        self.advanced = QGroupBox("Advanced", checkable=True, checked=False)
        inner = QFormLayout(self.advanced)
        self.options = {}
        for name, opt in pipeline.get("options", {}).items():
            if "choices" in opt:
                w = QComboBox()
                for choice in opt["choices"]:
                    w.addItem(choice, choice)
                w.currentIndexChanged.connect(self._changed)
            else:
                w = QSpinBox(minimum=0, maximum=1_000_000)
                w.valueChanged.connect(self._changed)
            self.options[name] = w
            inner.addRow(opt["label"], w)
        self.threads = QSpinBox(minimum=1, maximum=256, value=max(1, (os.cpu_count() or 2) - 1))
        inner.addRow("Threads", self.threads)
        self.advanced.toggled.connect(lambda on: [inner.itemAt(i).widget().setVisible(on) for i in range(inner.count())
                                                  if inner.itemAt(i).widget()])
        self.advanced.toggled.emit(False)
        form.addRow(self.advanced)
        buttons = QHBoxLayout()
        buttons.addWidget(QPushButton("Run", clicked=lambda: terminal.run_pipeline(self),
                                      toolTip="Types the pipeline's command into the session. Each run leaves a folder beside "
                                              "its output (<output>_run) with the inputs, command, log and reports."))
        buttons.addWidget(QPushButton("Install tools", clicked=lambda: terminal.install(self)))
        buttons.addStretch()
        form.addRow(buttons)
        form.addRow(QLabel(f"Tools: {', '.join(pipeline.get('needs', []))}", wordWrap=True))
        self.form = form
        self.apply_preset(self.preset.currentText())

    def _browse_row(self, edit, slot):
        row = QWidget()
        layout = QHBoxLayout(row)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(edit, 1)
        layout.addWidget(QPushButton("…", clicked=slot, maximumWidth=30))
        return row

    def refresh_references(self):
        """The reference list: this document, the other open ones, then any file picked before."""
        if not hasattr(self, "reference"):
            return
        files = [self.reference.itemData(i) for i in range(self.reference.count())
                 if isinstance(self.reference.itemData(i), str) and self.reference.itemData(i)]   # "" is File…
        self.reference.clear()
        self.reference.addItem("This document", None)
        for win in self.terminal.window.open_windows()[1:]:
            self.reference.addItem(f"Open: {win.editor.path.name}", win)
        for path in files:
            self.reference.addItem(Path(path).name, path)
        self.reference.addItem("File…", "")

    def _reference_picked(self):
        if self.reference.currentData() != "":
            return
        path, _ = QFileDialog.getOpenFileName(self, "Reference", "", runner.SEQ_FILTER)
        if not path:
            self.reference.setCurrentIndex(0)
            return
        self.reference.insertItem(self.reference.count() - 1, Path(path).name, path)
        self.reference.setCurrentIndex(self.reference.count() - 2)

    def apply_preset(self, name):
        values = runner.options_for(self.pipeline, name)
        for key, w in self.options.items():
            w.blockSignals(True)
            if isinstance(w, QComboBox):
                w.setCurrentIndex(w.findData(values[key]))
            else:
                w.setValue(int(values[key]))
            w.blockSignals(False)
        self.shown = runner.inputs_for(self.pipeline, name)
        for key, row in self.rows.items():
            self.form.setRowVisible(row, key in self.shown)
        if self.preset.findText(CUSTOM) >= 0:
            self.preset.removeItem(self.preset.findText(CUSTOM))
        self.preset.setCurrentText(name)
        self.preset_name = name
        self.terminal.show_diagram(self)

    def _changed(self):
        if self.values() != runner.options_for(self.pipeline, self.preset_name) and self.preset.findText(CUSTOM) < 0:
            self.preset.addItem(CUSTOM)
            self.preset.setCurrentText(CUSTOM)
        self.terminal.show_diagram(self)

    def values(self):
        return {k: w.currentData() if isinstance(w, QComboBox) else w.value() for k, w in self.options.items()}

    def _pick_reads(self, edit, name):
        path, _ = QFileDialog.getOpenFileName(self, READS[name], edit.text(), runner.FASTQ_FILTER)
        if path:
            edit.setText(path)
            if not self.out.text():
                self.out.setText(self.default_output())

    def _pick_extra(self):
        paths, _ = QFileDialog.getOpenFileNames(self, "Sequence files to align", "", runner.SEQ_FILTER)
        if paths:
            self.extra.setText("; ".join(paths))

    def _pick_out(self):
        path = (QFileDialog.getExistingDirectory(self, "Output folder") if not Path(self.pipeline["output"]).suffix
                else QFileDialog.getSaveFileName(self, "Output", self.out.text() or self.default_output())[0])
        if path:
            self.out.setText(path)

    def reads(self):
        return {k: e.text().strip() if k in self.shown else "" for k, e in self.edits.items()}

    def default_output(self):
        ed = self.terminal.window.editor
        first = next((p for p in self.reads().values() if p), "")
        return runner.default_output(self.pipeline, ed.fasta_record or ed.doc.name, first, ed.path)


def _stop_quietly(process):
    process.blockSignals(True)
    process.kill()
    process.waitForFinished(1000)


class Terminal(QWidget):
    def __init__(self, window):
        super().__init__()
        self.window, self.pid, self.busy, self.run, self.buffer, self.closing = window, None, False, None, "", False
        self.ref_key = None   # what $REF was last written from (see send)
        self.scratch = Path(tempfile.mkdtemp(prefix="opendna_terminal_"))
        self.diagram = PipelineDiagram()
        self.list = QListWidget()
        self.forms = QStackedWidget()
        self.keys = []
        for key, pipeline in runner.pipelines().items():
            self.keys.append(key)
            self.list.addItem(pipeline["label"])
            self.forms.addWidget(PipelineForm(self, key, pipeline))
        self.list.currentRowChanged.connect(self._pipeline_changed)
        left = QWidget()
        lv = QVBoxLayout(left)
        lv.setContentsMargins(0, 0, 0, 0)
        lv.addWidget(QLabel("Pipelines"))
        lv.addWidget(self.list)
        lv.addWidget(self.forms, 1)
        lv.addWidget(QLabel("Pipeline steps"))
        lv.addWidget(self.diagram)
        self.out = QTextBrowser(openLinks=False)
        self.out.setFont(mono_font(9))
        self.out.document().setMaximumBlockCount(20000)
        self.out.anchorClicked.connect(lambda url: self.open_file(url.toLocalFile()))
        self.line = HistoryLine()
        self.line.returnPressed.connect(self._typed)
        self.stop_button = QPushButton("Stop", clicked=self.stop, enabled=False)
        bottom = QHBoxLayout()
        bottom.addWidget(QLabel("$"))
        bottom.addWidget(self.line, 1)
        bottom.addWidget(QPushButton("+ file", clicked=self._insert_files, toolTip="Insert the path of files, quoted"))
        bottom.addWidget(self.stop_button)
        right = QWidget()
        rv = QVBoxLayout(right)
        rv.setContentsMargins(0, 0, 0, 0)
        rv.addWidget(self.out, 1)
        rv.addLayout(bottom)
        split = QSplitter()
        split.addWidget(left)
        split.addWidget(right)
        split.setStretchFactor(1, 1)
        QVBoxLayout(self).addWidget(split)
        self.process = QProcess(self)
        self.process.setProcessChannelMode(QProcess.MergedChannels)
        self.process.readyReadStandardOutput.connect(self._read)
        self.process.finished.connect(self._ended)
        self.process.errorOccurred.connect(self._error)
        # The pane can be closed, switched to another view or rebuilt (Reset Tabs): stop the shell with it, quietly
        # (its `finished` would otherwise reach _ended on a deleted widget and start a new shell).
        self.destroyed.connect(partial(_stop_quietly, self.process))
        self.list.setFixedHeight(self.list.sizeHintForRow(0) * self.list.count() + 2 * self.list.frameWidth() + 4)
        self.list.setCurrentRow(0)
        self.start()

    # ---- the session ---------------------------------------------------------------------------------------------------------------
    def start(self):
        args = runner.command_line(["-ls"])
        self.process.start(args[0], args[1:])
        ed = self.window.editor
        folder = Path(ed.path).parent if ed is not None and ed.path and Path(ed.path).parent.exists() else Path.home()
        try:
            runner.shell_path(folder)
        except ValueError:   # a folder WSL cannot reach (a network share): start in the home folder
            folder = Path.home()
        self.process.write(runner.session_setup(folder).encode())
        self._say(f"bash session{' in WSL' if runner.use_wsl() else ''}. Pick a pipeline on the left, or type a command.")

    def _error(self, error):
        if error == QProcess.FailedToStart:
            self._say(f"Cannot start {self.process.program()}." + (" Install WSL first: wsl --install" if runner.use_wsl() else ""))

    def _ended(self, *_):
        if self.closing:
            return
        started, self.pid, self.buffer = self.pid, None, ""
        self._finish(1)
        if not started:   # it never got going (no WSL distribution, no bash): do not loop
            self._say("The shell did not start." + (" Is a WSL distribution installed? wsl --install" if runner.use_wsl() else ""))
            return
        self._say("The shell ended; starting a new one.")
        self.start()

    def send(self, command, echo=None):
        if self.busy:
            return False
        if self.process.state() == QProcess.NotRunning:   # bash never started (no WSL): nothing would ever answer
            self._say("The shell is not running." + (" Is a WSL distribution installed? wsl --install" if runner.use_wsl() else ""))
            self.run = None
            return False
        ed = self.window.editor
        exports = ""
        if ed is not None:
            ref = self.scratch / f"{runner.safe_name(ed.fasta_record or ed.doc.name)}.fa"
            key = (id(ed), ed.version, ed.fasta_record, str(ed.path))   # rewritten only once the document changed
            try:
                if key != self.ref_key or not ref.exists():
                    write_fasta(ref, document_records(ed))
                    self.ref_key = key
                exports = f"export REF={shlex.quote(runner.shell_path(ref))}"
                if ed.path and Path(ed.path).exists():
                    exports += f" DOC={shlex.quote(runner.shell_path(ed.path))}"
            except (OSError, ValueError) as e:
                self._say(f"$REF not updated: {e}")
        self._append(f"$ {echo or command}", bold=True)
        self.busy = True
        self.line.setEnabled(False)
        self.stop_button.setEnabled(True)
        self.process.write(((exports + "\n") if exports else "").encode() + runner.wrap(command).encode())
        return True

    def _typed(self):
        command = self.line.text().strip()
        if command and self.send(command):
            self.line.remember(command)
            self.line.clear()

    def _read(self):
        self.buffer += bytes(self.process.readAllStandardOutput()).decode("utf-8", "replace")
        # shortcut: output shows line by line (a progress bar drawn with \r waits for its newline)
        *lines, self.buffer = self.buffer.split("\n")
        plain = []   # output lines between markers, shown in one go (one append and one log write per chunk)
        for line in lines:
            line = line.rstrip("\r")
            if not line.startswith(("__PID__ ", "__DONE__ ", "__OPEN__ ", "__OPENMSA__ ")):
                plain.append(line)
                continue
            self._append(plain)
            plain = []
            if line.startswith("__PID__ "):
                self.pid, self.keep = line.split()[1:3]
            elif line.startswith("__DONE__ "):
                self._finish(int(line.split()[1]))
            else:
                tag, path = line.split(" ", 1)
                self.open_file(path, msa=tag == "__OPENMSA__")
        self._append(plain)

    def _finish(self, code):
        if not self.busy:
            return
        self.busy = False
        self.line.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.line.setFocus()
        run = self.run
        if run is None:
            return
        # still the current run while these are written, so they also go into its log.txt
        self._append(runner.summary(run.folder), bold=True)
        self._append([f"Report: {path}" for path in runner.reports(run.folder)])
        self._append("Done." if code == 0 else f"Failed (exit code {code}): see the log above, or {run.folder / 'log.txt'}",
                     bold=True)
        self.run = None

    def stop(self):
        if self.busy and self.pid:
            args = runner.command_line(["-c", runner.stop_children(self.pid, self.keep)])
            QProcess.execute(args[0], args[1:])

    def close_session(self):
        self.closing = True
        self.process.kill()
        self.process.waitForFinished(2000)

    # ---- output --------------------------------------------------------------------------------------------------------------------
    def _append(self, lines, bold=False):
        """Show output (a line, or a list of lines), paths made clickable, and add it to the run's log."""
        lines = [lines] if isinstance(lines, str) else lines
        if not lines:
            return
        if self.run is not None:
            with open(self.run.folder / "log.txt", "a", encoding="utf-8") as fh:
                fh.write("\n".join(lines) + "\n")
        shown = []
        for text in lines:
            parts = []
            for i, piece in enumerate(PATH_RE.split(text)):
                if i % 2:
                    href = QUrl.fromLocalFile(runner.windows_path(piece)).toString()
                    parts.append(f'<a href="{html.escape(href)}">{html.escape(piece)}</a>')
                else:
                    parts.append(html.escape(piece))
            shown.append("".join(parts))
        body = "<br>".join(shown)
        self.out.append(f'<span style="white-space:pre-wrap">{"<b>" + body + "</b>" if bold else body}</span>')

    def _say(self, text):
        self.out.append(f'<span style="color:gray">{html.escape(text)}</span>')

    def _insert_files(self):
        paths, _ = QFileDialog.getOpenFileNames(self, "Insert file paths", "", "All files (*)")
        try:
            self.line.insert(" ".join(shlex.quote(runner.shell_path(p)) for p in paths))
        except ValueError as e:
            self._say(str(e))
        self.line.setFocus()

    # ---- bringing results back -----------------------------------------------------------------------------------------------------
    def open_file(self, path, msa=False):
        """What `open` (or a click on a path) does: see runner.open_kind."""
        path = runner.windows_path(path)
        try:
            kind = runner.open_kind(path, msa)
            if kind == "browser":
                QDesktopServices.openUrl(QUrl.fromLocalFile(path))
            elif kind == "document":
                self.window.open_path(path)
            elif kind == "track":
                self._target().load_track(path)
            else:
                self._show_alignment(path)
        except (OSError, ValueError) as e:
            self._say(f"Cannot open {path}: {e}")

    def _target(self):
        """The window a mapping's BAM belongs to: the one of its reference (opened now if it was a file)."""
        target = self.run.target if self.run is not None else None
        if isinstance(target, str):
            target = self.window.open_path(target) or self.window
            self.run.target = target
        return target or self.window

    def _show_alignment(self, path):
        names, rows = msa.parse_aligned_fasta(Path(path).read_text(encoding="utf-8"))
        view = self.window.workspace.ensure(MsaView)
        view.set_alignment(msa.Alignment(names, rows))
        self._say(f"{len(names)} sequences shown in the Multiple Alignment view.")

    # ---- pipelines -----------------------------------------------------------------------------------------------------------------
    def _pipeline_changed(self, row):
        self.forms.setCurrentIndex(row)
        self.forms.currentWidget().refresh_references()
        self.show_diagram(self.forms.currentWidget())

    def show_diagram(self, form):
        if form is self.forms.currentWidget():
            self.diagram.set_nodes(runner.diagram(form.pipeline, form.values(), form.shown))

    def select(self, key):
        self.list.setCurrentRow(self.keys.index(key))

    def install(self, form):
        tools = form.pipeline.get("needs", [])
        if QMessageBox.question(self, "Install tools", f"Install {', '.join(tools)} from bioconda into the conda environment "
                                f"'{runner.ENV}'{' inside WSL' if runner.use_wsl() else ''}? This can download several hundred MB."
                                ) == QMessageBox.Yes:
            self.send(runner.install(tools), echo=f"install {' '.join(tools)} into {runner.ENV}")

    def run_pipeline(self, form):
        if self.busy:
            self._say("A command is running: wait for it or press Stop.")
            return
        try:
            command, folder, target = self._prepare(form)
        except (OSError, ValueError) as e:
            QMessageBox.warning(self, form.pipeline["label"], str(e))
            return
        (folder / "command.txt").write_text(command + "\n", encoding="utf-8")
        self.run = SimpleNamespace(folder=folder, target=target)
        self.send(command)

    def _prepare(self, form):
        """(command, run folder, target window) for a run, its inputs checked and written to the run folder."""
        reads = form.reads()
        for key, path in reads.items():
            if path and not Path(path).is_file():
                raise ValueError(f"{READS[key]}: '{path}' was not found.")
        if any(k in form.shown for k in READS) and not any(reads.values()):
            raise ValueError("Pick a reads file first.")
        if "reads2" in form.shown and reads.get("reads2") and not reads.get("reads1"):
            raise ValueError("Mate reads need the first reads file too.")
        out = form.out.text().strip() or form.default_output()
        if not Path(out).parent.is_dir():
            raise ValueError(f"The folder of '{out}' does not exist.")
        form.out.setText(out)
        folder = runner.run_folder(out)
        files, target, to_write = dict(reads, out=out), None, []
        if "reference" in form.shown:
            target = form.reference.currentData()
            seqs = (document_records(self.window.editor) if target is None else document_records(target.editor)
                    if not isinstance(target, str) else files_records([target]))
            files["ref"] = str(folder / "ref.fa")
            to_write.append((files["ref"], seqs))
        if "sequences" in form.shown:
            extra = [p.strip() for p in form.extra.text().split(";") if p.strip()]
            seqs = (document_records(self.window.editor) if form.own.isChecked() else []) + files_records(extra)
            if len(seqs) < 2:
                raise ValueError("Pick at least two sequences to align.")
            files["seqs"] = str(folder / "sequences.fa")
            to_write.append((files["seqs"], seqs))
        folder.mkdir(parents=True, exist_ok=True)
        for path, seqs in to_write:
            run_busy(self, "Writing the sequences…", lambda: write_fasta(path, seqs))
        command, _ = runner.build(form.pipeline, form.values(), files, form.threads.value(), folder)
        return command, folder, target
