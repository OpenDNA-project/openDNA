"""Runs external command-line tools (bwa, minimap2, samtools…) and opens their output.

Recipes live in tools.toml (one entry per tool, no code). On Windows the commands run in
WSL through `wsl -e bash -lc`, with Windows paths turned into /mnt/<drive>/… paths;
elsewhere through bash directly. Tools are looked up first in a dedicated conda env
(opendna-tools, see INSTALL), then on the normal PATH.

Each run gets a folder next to its output (<out>_run/) with the reference, the exact
commands (command.txt) and the full log (log.txt), so it can be re-run by hand.
"""
import os
import shlex
import sys
import tomllib
from datetime import datetime
from pathlib import Path

from PySide6.QtCore import QProcess
from PySide6.QtWidgets import (
    QComboBox, QDialog, QDialogButtonBox, QFileDialog, QFormLayout, QHBoxLayout, QLabel, QLineEdit,
    QMessageBox, QPlainTextEdit, QProgressBar, QPushButton, QSpinBox, QVBoxLayout, QWidget,
)

import wip

RECIPES = Path(__file__).with_name("tools.toml")
ENV = "opendna-tools"
# Where the tools env usually is (Miniforge, Miniconda, Mambaforge); searched before PATH.
TOOL_DIRS = [f"$HOME/{base}/envs/{ENV}/bin" for base in ("miniforge3", "miniconda3", "mambaforge")]
TOOLS = ["bwa", "samtools", "minimap2"]
CHECK = "; ".join(f'printf "%-9s " {t}; command -v {t} || echo "MISSING"' for t in TOOLS)
INSTALL = (
    'for c in "$HOME"/miniforge3/bin/conda "$HOME"/miniconda3/bin/conda "$HOME"/mambaforge/bin/conda; do '
    f'if [ -x "$c" ]; then "$c" create -y -n {ENV} -c conda-forge -c bioconda {" ".join(TOOLS)}; exit $?; fi; done; '
    'echo "conda not found: install Miniforge first (https://conda-forge.org/download/)"; exit 1'
)
FASTQ_FILTER = "Reads (*.fastq *.fq *.fastq.gz *.fq.gz *.fasta *.fa *.fa.gz);;All files (*)"


def recipes():
    with open(RECIPES, "rb") as fh:
        return tomllib.load(fh)


def use_wsl():
    return sys.platform == "win32"


def shell_path(path):
    """A path as the shell sees it: on Windows V:\\data\\x.fq -> /mnt/v/data/x.fq (WSL)."""
    path = str(Path(path).resolve())
    if not use_wsl():
        return path
    if len(path) < 2 or path[1] != ":":
        raise ValueError(f"{path}: only paths on a drive letter can be reached from WSL")
    return f"/mnt/{path[0].lower()}{path[2:].replace(chr(92), '/')}"


def fill(template, files, threads):
    """Recipe text -> shell text. File values are quoted (never pasted raw); an empty one
    (e.g. no mate reads) disappears."""
    values = {k: shlex.quote(shell_path(v)) if v else "" for k, v in files.items()}
    return template.format(threads=int(threads), **values)


def script(command):
    """Full bash script: tools env first on PATH, a pipeline fails if any step fails, and
    the shell's PID is printed first so Cancel can stop it (see RunDialog.cancel)."""
    return f'echo "__PID__ $$"; export PATH="{":".join(TOOL_DIRS)}:$PATH"; set -o pipefail; {command}'


def command_line(text):
    return (["wsl", "-e", "bash", "-lc", text] if use_wsl() else ["bash", "-lc", text])


def steps_for(recipe, files, threads):
    """[shell command] to run for a recipe, `prepare` first unless already done."""
    steps = []
    done = recipe.get("prepared")
    if recipe.get("prepare") and not (done and Path(done.format(**files)).exists()):
        steps.append(fill(recipe["prepare"], files, threads))
    steps.append(fill(recipe["run"], files, threads))
    return steps


class RunDialog(QDialog):
    """Runs shell commands one after another with a live log and a Cancel button. Also
    writes command.txt / log.txt into `folder` (if given)."""

    def __init__(self, parent, title, steps, folder=None):
        super().__init__(parent)
        self.setWindowTitle(title)
        self.resize(820, 480)
        self.steps, self.folder = list(steps), folder
        self.ok, self.pid, self.cancelled = False, None, False
        self.log = QPlainTextEdit(readOnly=True)
        self.log.setMaximumBlockCount(20000)
        self.log.setLineWrapMode(QPlainTextEdit.NoWrap)
        self.status = QLabel("Running…")
        self.busy = QProgressBar(minimum=0, maximum=0, textVisible=False, maximumWidth=160)   # moving while a step runs
        self.button = QPushButton("Cancel", clicked=self.ask_cancel)
        bottom = QHBoxLayout()
        bottom.addWidget(self.status, 1)
        bottom.addWidget(self.busy)
        bottom.addWidget(self.button)
        layout = QVBoxLayout(self)
        layout.addWidget(self.log)
        layout.addLayout(bottom)
        self.process = QProcess(self)
        self.process.setProcessChannelMode(QProcess.MergedChannels)
        self.process.readyReadStandardOutput.connect(self._read)
        self.process.finished.connect(self._finished)
        self.process.errorOccurred.connect(self._error)
        if folder:
            Path(folder).mkdir(parents=True, exist_ok=True)
            (Path(folder) / "command.txt").write_text("\n".join(self.steps) + "\n", encoding="utf-8")
        self._next()

    def _append(self, text):
        self.log.appendPlainText(text.rstrip("\n"))
        if self.folder:
            with open(Path(self.folder) / "log.txt", "a", encoding="utf-8") as fh:
                fh.write(text if text.endswith("\n") else text + "\n")

    def _next(self):
        command = self.steps.pop(0)
        self._append(f"$ {command}")
        args = command_line(script(command))
        self.process.start(args[0], args[1:])

    def _read(self):
        text = bytes(self.process.readAllStandardOutput()).decode("utf-8", "replace")
        lines = []
        for line in text.splitlines():
            if line.startswith("__PID__ "):
                self.pid = line.split()[1]
            else:
                lines.append(line)
        if lines:
            self._append("\n".join(lines))

    def _error(self, error):
        if error == QProcess.FailedToStart:
            hint = "Install WSL (wsl --install) and the tools (Tools → Install Mapping Tools)." if use_wsl() else ""
            self._append(f"Cannot start {self.process.program()}. {hint}")
            self._done(False)

    def _finished(self, code, _status):
        self._read()
        if self.cancelled or code != 0:
            self._append("Cancelled." if self.cancelled else f"Failed (exit code {code}).")
            self._done(False, "" if self.cancelled else f"Failed (exit code {code}).")
        elif self.steps:
            self._next()
        else:
            self._done(True)

    def _done(self, ok, why=""):
        self.ok = ok
        self.status.setText("Done." if ok else why or "Stopped.")
        self.busy.hide()
        self.button.setText("Close")
        self.button.clicked.disconnect()
        self.button.clicked.connect(self.accept if ok else self.reject)

    def running(self):
        return self.process.state() != QProcess.NotRunning

    def ask_cancel(self):
        """Cancel, Escape or the window's close button while a step runs: asked first (No by default: a long mapping is lost),
        then the run is stopped. True when nothing runs any more."""
        if not self.running():
            return True
        if QMessageBox.question(self, self.windowTitle(), "Stop the run? What it has written so far is kept in its run folder.",
                                QMessageBox.Yes | QMessageBox.No, QMessageBox.No) != QMessageBox.Yes:
            return False
        self.cancel()
        self.process.waitForFinished(3000)
        return True

    def cancel(self):
        """Stop the running command and the programs it started. Killing wsl.exe alone
        would leave bwa running inside WSL, so the shell's children are stopped first."""
        self.cancelled = True
        if self.pid:
            stop = f"pkill -TERM -P {int(self.pid)}; kill -TERM {int(self.pid)}"
            args = command_line(stop)
            QProcess.execute(args[0], args[1:])
        self.process.kill()

    def reject(self):
        # Escape comes here, not through closeEvent: without this it hid the window while bwa / minimap2 went on running
        if self.ask_cancel():
            super().reject()

    def closeEvent(self, event):
        if not self.ask_cancel():
            event.ignore()
            return
        super().closeEvent(event)


class MapReadsDialog(QDialog):
    """Pick a mapping recipe, reads and the output BAM."""

    def __init__(self, parent, seq_name):
        super().__init__(parent)
        self.setWindowTitle(f"Map Reads to {seq_name}")
        self.recipes = recipes()
        self.recipe = QComboBox()
        for key, r in self.recipes.items():
            self.recipe.addItem(r.get("label", key), key)
        self.reads1, self.reads2 = QLineEdit(), QLineEdit(placeholderText="optional: mate reads (paired-end)")
        self.out = QLineEdit()
        self.threads = QSpinBox(minimum=1, maximum=256, value=max(1, (os.cpu_count() or 2) - 1))
        form = QFormLayout(self)
        mapper = QHBoxLayout()   # the "WIP" tag while mapping is unfinished (wip.FEATURES)
        mapper.addWidget(self.recipe, 1)
        mapper.addWidget(wip.Badge("mapping"))
        form.addRow("Mapper", mapper)
        form.addRow("Reads", self._with_browse(self.reads1, self._pick_reads1))
        self.reads2_row = self._with_browse(self.reads2, lambda: self._pick(self.reads2))
        form.addRow("Mate reads", self.reads2_row)
        form.addRow("Output BAM", self._with_browse(self.out, self._pick_out))
        form.addRow("Threads", self.threads)
        form.addRow(QLabel("Runs in WSL on Windows. The run folder (<output>_run) keeps the reference, "
                           "commands and log.", wordWrap=True))
        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.button(QDialogButtonBox.Ok).setText("Run")
        buttons.accepted.connect(self._accept)
        buttons.rejected.connect(self.reject)
        form.addRow(buttons)
        self.form, self.seq_name = form, seq_name
        self.recipe.currentIndexChanged.connect(self._recipe_changed)
        self._recipe_changed()

    def _with_browse(self, edit, slot):
        row = QWidget()
        layout = QHBoxLayout(row)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(edit, 1)
        layout.addWidget(QPushButton("…", clicked=slot, maximumWidth=30))
        return row

    def _recipe_changed(self):
        paired = "{reads2}" in self.current()["run"]
        self.form.setRowVisible(self.reads2_row, paired)

    def current(self):
        return self.recipes[self.recipe.currentData()]

    def _pick(self, edit):
        path, _ = QFileDialog.getOpenFileName(self, "Reads", edit.text(), FASTQ_FILTER)
        if path:
            edit.setText(path)
        return path

    def _pick_reads1(self):
        if self._pick(self.reads1) and not self.out.text():
            stem = Path(self.reads1.text()).name.split(".")[0]
            for tag in ("_R1", "_1"):  # sample_R1.fastq.gz -> sample
                stem = stem[: -len(tag)] if stem.endswith(tag) else stem
            self.out.setText(str(Path(self.reads1.text()).with_name(f"{stem}_{safe_name(self.seq_name)}.bam")))

    def _pick_out(self):
        path, _ = QFileDialog.getSaveFileName(self, "Output BAM", self.out.text(), "BAM (*.bam)")
        if path:
            self.out.setText(path)

    def _accept(self):
        for edit, what, needed in ((self.reads1, "reads", True), (self.reads2, "mate reads", "{reads2}" in self.current()["run"])):
            path = edit.text().strip()
            if needed and (path or edit is self.reads1) and not Path(path).is_file():
                QMessageBox.warning(self, "Map reads", f"The {what} file '{path}' was not found: use … to pick it." if path
                                    else "Pick the reads file (FASTQ) first.")
                edit.setFocus()
                return
        out = self.out.text().strip()
        if not out or Path(out).name.lower() in (".bam", ""):
            QMessageBox.warning(self, "Map reads", "Choose where to write the output BAM file.")
            self.out.setFocus()
            return
        if not Path(out).parent.is_dir():
            QMessageBox.warning(self, "Map reads", f"The folder of '{out}' does not exist.")
            self.out.setFocus()
            return
        self.out.setText(out if out.lower().endswith(".bam") else out + ".bam")
        self.accept()

    def plan(self, reference):
        """(run folder, steps, result BAM) once the reference is written to `reference`."""
        out = Path(self.out.text())
        files = {"ref": str(reference), "reads1": self.reads1.text(),
                 # The recipe decides (not the row's visibility: the dialog is closed by now).
                 "reads2": self.reads2.text() if "{reads2}" in self.current()["run"] else "",
                 "out": str(out.with_suffix(""))}
        recipe = self.current()
        steps = steps_for(recipe, files, self.threads.value())
        return steps, Path(recipe["result"].format(**files))


def safe_name(name):
    """A sequence name usable in a FASTA header and a file name."""
    return "".join(c if c.isalnum() or c in "._-" else "_" for c in name) or "sequence"


def run_folder(out_bam):
    folder = Path(out_bam).with_name(Path(out_bam).stem + "_run")
    if folder.exists() and any(folder.iterdir()):
        folder = folder.with_name(folder.name + datetime.now().strftime("_%Y%m%d_%H%M%S"))
    return folder
