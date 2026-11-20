"""Help > Run Setup Again: the installer's choices (packaging/openDNA.iss) for a user who installed an older version or skipped them,
without reinstalling. Each row shows the current state; nothing changes until OK. The work is done by what the installer's command
lines run: file_assoc.register (--register-file-types), the Performance setting `cores` (--set-cores), and the install of
pLannotate and the Terminal's tools (oDNA-cli wsltools install, offered only when WSL works, as the installer's WslReady)."""
import sys

from PySide6.QtWidgets import QCheckBox, QDialog, QDialogButtonBox, QFormLayout, QLabel, QSpinBox

import file_assoc
import plannotate_tool
import runner
import tool_settings
from annotation import InstallDialog
from uikit import run_busy


def run(parent):
    tools = runner.terminal_tools()
    state, have = run_busy(parent, "Checking pLannotate and the Terminal's tools…",   # (it starts WSL: a few seconds)
                           lambda: (plannotate_tool.status(), set(runner.installed_tools(tools))))
    dialog = QDialog(parent)
    dialog.setWindowTitle("Run Setup Again")
    form = QFormLayout(dialog)
    form.addRow(QLabel("The installer's choices. Nothing changes until you press OK."))

    assoc = None
    if sys.platform == "win32":   # elsewhere there is no "Open with" registry to write
        offered = file_assoc.is_registered_for_current_exe()
        assoc = QCheckBox('Offer openDNA in "Open with" for sequence files', checked=offered)
        note = "Now: offered." if offered else "Now: not offered."
        if not file_assoc.supported():
            assoc.setEnabled(False)
            note = "Only the installed openDNA.exe can offer itself, not openDNA run from source."
        form.addRow(assoc)
        form.addRow("", QLabel(note))

    cores = QSpinBox(minimum=1, maximum=tool_settings.CPUS, value=tool_settings.cores(), suffix=f" of {tool_settings.CPUS}")
    form.addRow("CPU cores openDNA may use:", cores)
    form.addRow("", QLabel(f"For Local Align and Folder Search. Default {tool_settings.DEFAULT_CORES}; also in Settings > General > Performance."))

    install = None
    if state.state == "no_wsl":
        form.addRow(QLabel(f"pLannotate and the Terminal's tools: {plannotate_tool.WSL_HELP}", wordWrap=True))
    else:
        missing_plannotate, missing_tools = state.state != "ready", len(have) < len(tools)
        install = QCheckBox("Set up pLannotate and the Terminal's tools" + (" in WSL2" if runner.use_wsl() else "")
                            + " now (several hundred MB, a few minutes)", enabled=missing_plannotate or missing_tools)
        form.addRow(install)
        form.addRow("", QLabel(f"Now: {state.text}; Terminal tools: {len(have)} of {len(tools)} installed."))

    buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel, accepted=dialog.accept, rejected=dialog.reject)
    form.addRow(buttons)
    if dialog.exec() != QDialog.Accepted:
        return
    if cores.value() != tool_settings.cores():
        tool_settings.set_value(tool_settings.PERFORMANCE, "cores", cores.value())
    if assoc is not None and assoc.isEnabled() and assoc.isChecked() != offered:
        (file_assoc.register_now if assoc.isChecked() else file_assoc.unregister_now)(parent)
    if install is not None and install.isChecked():
        InstallDialog(parent, plannotate=missing_plannotate, terminal=missing_tools).exec()
        plannotate_tool.forget()
