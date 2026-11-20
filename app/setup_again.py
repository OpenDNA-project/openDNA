"""Help > Run Setup Again: the installer's choices (packaging/openDNA.iss) for a user who installed an older version or skipped them,
without reinstalling. Each row shows the current state; nothing changes until OK. The work is done by what the installer's command
lines run: file_assoc.register (--register-file-types), the Performance setting `cores` (--set-cores), and the install of
pLannotate and the Terminal's tools (oDNA-cli wsltools install, offered only when WSL works, as the installer's WslReady). One choice the
installer does not make: ViennaRNA, which the Synthesis tab's calculator needs and which is downloaded on request (vienna_install)."""
import sys
import threading

from PySide6.QtCore import QTimer
from PySide6.QtWidgets import QCheckBox, QDialog, QDialogButtonBox, QFormLayout, QLabel, QMessageBox

import file_assoc
import plannotate_tool
import runner
import online
import tool_settings
from annotation import InstallDialog
from uikit import run_busy
from user_error import UserError


def run(parent):
    dialog = QDialog(parent)
    dialog.setWindowTitle("Run Setup Again")
    form = QFormLayout(dialog)
    form.addRow(QLabel("The installer's choices. Nothing changes until you press OK."))

    assoc = None
    if sys.platform == "win32":   # elsewhere there is no "Open with" registry to write
        found = set(file_assoc.status().values()) if file_assoc.supported() else set()
        offered = found == {file_assoc.THIS}
        note = ("Only the installed openDNA.exe can offer itself, not openDNA run from source." if not file_assoc.supported()
                else "Now: offered." if offered
                else "Now: offered for another copy of openDNA (OK makes it this one)." if file_assoc.OTHER in found
                else "Now: offered for some file types only (OK adds the others)." if file_assoc.THIS in found
                else "Now: not offered.")
        assoc = QCheckBox('Offer openDNA in "Open with" for sequence files', checked=offered, enabled=file_assoc.supported())
        form.addRow(assoc)
        form.addRow("", QLabel(note))

    param = next(p for p in tool_settings.TOOLS[tool_settings.PERFORMANCE] if p.key == "cores")
    cores = tool_settings.SettingsDialog.param_widget(param, tool_settings.cores())   # its maximum is this PC's number of cores
    shown_cores = cores.value()   # (a saved value above the maximum shows as the maximum: OK alone must not save it)
    form.addRow(param.label + ":", cores)
    form.addRow("", QLabel(f"Default {param.default}; also in Settings > General > Performance.", toolTip=param.tip))

    # The tools are asked in a thread (it starts WSL: seconds, or up to the timeout if WSL hangs): the dialog is usable meanwhile.
    install = QCheckBox("Set up pLannotate and the Terminal's tools" + (" in WSL2" if runner.use_wsl() else "")
                        + " now (several hundred MB, a few minutes)", enabled=False)
    tools_note = QLabel("Now: checking…", wordWrap=True)
    form.addRow(install)
    form.addRow("", tools_note)
    answer, missing = [], {}
    threading.Thread(target=lambda: answer.append(plannotate_tool.tools_state()), daemon=True).start()

    def show_tools():
        if not answer:
            return
        poll.stop()
        state, tools = answer[0]
        if state.state == "no_wsl":
            install.setVisible(False)
            tools_note.setText(f"pLannotate and the Terminal's tools: {plannotate_tool.WSL_HELP}")
            return
        missing.update(plannotate=state.state != "ready", terminal=not all(tools.values()))
        install.setEnabled(any(missing.values()))
        tools_note.setText(f"Now: {state.text}; Terminal tools: {sum(tools.values())} of {len(tools)} installed.")
    poll = QTimer(dialog, interval=100, timeout=show_tools)   # the thread only appends; the widgets are touched here, in the window's thread
    poll.start()

    import synthesis_salis
    have_vienna = synthesis_salis.have_vienna()
    vienna = QCheckBox("Install ViennaRNA for the Synthesis tab's calculator (2 to 14 MB from pypi.org)", checked=False, enabled=not have_vienna)
    online.mark(vienna, "nothing of your sequences: it downloads the ViennaRNA package from pypi.org")
    form.addRow(vienna)
    form.addRow("", QLabel("Now: installed." if have_vienna else "Now: not installed (the calculator cannot run without it)."))

    form.addRow(QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel, accepted=dialog.accept, rejected=dialog.reject))
    accepted = dialog.exec() == QDialog.Accepted
    n, register, setup = cores.value(), assoc is not None and assoc.isChecked() != offered, install.isEnabled() and install.isChecked()
    get_vienna = vienna.isChecked()
    poll.stop()
    dialog.deleteLater()   # (a child of the window: it would stay until the window closes)
    if not accepted:
        return
    if n != shown_cores:
        tool_settings.set_value(tool_settings.PERFORMANCE, "cores", n)
    if get_vienna:
        import vienna_install
        try:
            version = run_busy(parent, "Downloading ViennaRNA from pypi.org…", vienna_install.install)
            QMessageBox.information(parent, "ViennaRNA", f"ViennaRNA {version} is installed: the Synthesis tab's calculator can run.")
        except UserError as e:
            QMessageBox.warning(parent, "ViennaRNA", str(e))
    if setup:   # before the "Open with" notice, which does not block: the two would be on screen together
        InstallDialog(parent, **missing).exec()
        plannotate_tool.forget()
    if register:
        (file_assoc.unregister_now if offered else file_assoc.register_now)(parent)
