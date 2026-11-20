"""Help > Run Setup Again: the installer's choices (packaging/openDNA.iss) for a user who installed an older version or skipped them,
without reinstalling. Each row shows the current state; nothing changes until OK. The work is done by what the installer's command
lines run: file_assoc.register (--register-file-types), the Performance setting `cores` (--set-cores), and the install of
pLannotate and the Terminal's tools (oDNA-cli wsltools install, offered only when WSL works, as the installer's WslReady)."""
import sys

from PySide6.QtWidgets import QCheckBox, QDialog, QDialogButtonBox, QFormLayout, QLabel

import file_assoc
import plannotate_tool
import runner
import tool_settings
from annotation import InstallDialog
from uikit import run_busy


def run(parent):
    state, tools = run_busy(parent, "Checking pLannotate and the Terminal's tools…", plannotate_tool.tools_state)   # (it starts WSL: a few seconds)
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
                else "Now: not offered.")
        assoc = QCheckBox('Offer openDNA in "Open with" for sequence files', checked=offered, enabled=file_assoc.supported())
        form.addRow(assoc)
        form.addRow("", QLabel(note))

    param = next(p for p in tool_settings.TOOLS[tool_settings.PERFORMANCE] if p.key == "cores")
    cores = tool_settings.SettingsDialog._widget(param, tool_settings.cores())   # the Settings page's own spinbox: same limits and suffix
    form.addRow(param.label + ":", cores)
    form.addRow("", QLabel(f"Default {param.default}; also in Settings > General > Performance.", toolTip=param.tip))

    install = None
    if state.state == "no_wsl":
        form.addRow(QLabel(f"pLannotate and the Terminal's tools: {plannotate_tool.WSL_HELP}", wordWrap=True))
    else:
        missing_plannotate, missing_tools = state.state != "ready", not all(tools.values())
        install = QCheckBox("Set up pLannotate and the Terminal's tools" + (" in WSL2" if runner.use_wsl() else "")
                            + " now (several hundred MB, a few minutes)", enabled=missing_plannotate or missing_tools)
        form.addRow(install)
        form.addRow("", QLabel(f"Now: {state.text}; Terminal tools: {sum(tools.values())} of {len(tools)} installed."))

    form.addRow(QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel, accepted=dialog.accept, rejected=dialog.reject))
    if dialog.exec() != QDialog.Accepted:
        return
    if cores.value() != tool_settings.cores():
        tool_settings.set_value(tool_settings.PERFORMANCE, "cores", cores.value())
    if assoc is not None and assoc.isChecked() != offered:
        (file_assoc.register_now if assoc.isChecked() else file_assoc.unregister_now)(parent)
    if install is not None and install.isChecked():
        InstallDialog(parent, plannotate=missing_plannotate, terminal=missing_tools).exec()
        plannotate_tool.forget()
