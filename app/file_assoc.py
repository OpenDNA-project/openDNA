"""Let the .exe offer to open .gb / .dna / FASTA files (double-click, "Open with", Settings > Default apps) without an installer.

What is written (all under HKEY_CURRENT_USER, so no administrator rights and no UAC prompt; nothing outside the user's own hive):

    Software\\Classes\\openDNA.SequenceFile                 the "ProgID": name, icon, and the command that opens a file
        DefaultIcon                                           the exe's own icon
        shell\\open\\command                                   "<exe>" "%1"
    Software\\Classes\\<ext>\\OpenWithProgids                    one value per extension: "openDNA can open this" (what fills the
                                                              "Open with" list and Settings > Default apps > by file type)
    Software\\Classes\\Applications\\openDNA.exe              SupportedTypes, FriendlyAppName, command: shows the program in "Open with >
                                                              Choose another app" when it is picked by file name
    Software\\openDNA\\Capabilities                         name, description, FileAssociations (ext -> ProgID)
    Software\\RegisteredApplications  value "openDNA"       points at the Capabilities key: this is what makes the program appear in
                                                              Settings > Default apps as an app of its own

What the user sees in Windows
    * Right-click a .gb file > Open with: openDNA is in the list (also under "Choose another app").
    * Double-click a .gb file that has no default program yet: Windows asks "How do you want to open this file?" and openDNA is
      one of the choices.
    * Settings > Apps > Default apps: type "openDNA" in the search box (or use the extension box): the page lists every extension
      the program offers, each with a "set default" button.

Caveats (they are Windows' rules, not ours; the program never tries to get around them)
    * Since Windows 8 a program CANNOT make itself the default for an extension. The choice is stored in
      HKCU\\...\\FileExts\\<ext>\\UserChoice together with a hash computed from the user's SID, the extension, the time and the ProgID, and
      Windows ignores a UserChoice written without a valid hash. Writing that hash by hand is exactly what browser hijackers do, and
      it breaks with every Windows update, so it is NOT done here. Registering only makes openDNA a candidate; the user clicks
      "openDNA" in Settings (we open the right page and say what to click: open_default_apps_settings()).
    * `ms-settings:defaultapps?registeredAppUser=openDNA` opens the page of this very app (the value is the name under
      HKCU\\Software\\RegisteredApplications, written by register()). Microsoft documents it for Windows 11 21H2/22H2 with the April 2023
      update and 23H2 or later ("Launch the Default Apps settings page", learn.microsoft.com). On Windows 10 and older Windows 11 builds the
      query string is ignored and the general Default apps page opens: the user then searches for openDNA.
    * The exe is not code-signed. SmartScreen ("Windows protected your PC") judges the file when it is launched, whichever way it is
      launched, and a downloaded file keeps its "downloaded from the Internet" mark: registering does not remove it nor bypass it. If it
      keeps appearing, Properties > "Unblock" on the .exe removes the mark. Nothing here changes any security setting.
    * A one-file PyInstaller exe unpacks itself into a temporary folder at each start: opening a file by double-click takes a few
      seconds longer than for an installed program.
    * The command holds the exe's absolute path. If the exe is moved, renamed or replaced at another place, the association keeps
      pointing at the old path: is_registered_for_current_exe() sees that (status "other") and the startup question is asked again,
      but only for a user who said Yes before (see should_ask).
    * Software\\openDNA may also hold values that older versions kept there with QSettings (settings are .ini files in the data
      folder now), so unregister() removes only its Capabilities subkey, never Software\\openDNA itself.

The registry is reached through a small object (WinRegistry for the real one; the tests use a dictionary-backed fake), so all of the
logic runs, and is tested, on any OS.
"""
import os
import sys
from pathlib import Path

from PySide6.QtCore import QSettings, Qt
from PySide6.QtWidgets import QDialog, QHBoxLayout, QLabel, QMessageBox, QPushButton, QVBoxLayout

from primer_db import data_dir
from project import FILE_TYPES, NAME

PROGID = "openDNA.SequenceFile"
# key under Software\Classes\Applications: by Windows' rules the exe's own file name (the released one is renamed: project.EXE_NAME)
APP_KEY = Path(sys.executable).name if getattr(sys, "frozen", False) else "openDNA.exe"
CLASSES = r"Software\Classes"
CAPABILITIES = r"Software\openDNA\Capabilities"
REGISTERED_APPS = r"Software\RegisteredApplications"
USER_CHOICE = r"Software\Microsoft\Windows\CurrentVersion\Explorer\FileExts\{ext}\UserChoice"   # read only: which program the user chose
SETTINGS_URI = "ms-settings:defaultapps"
SETTINGS_APP_URI = f"{SETTINGS_URI}?registeredAppUser={NAME}"
DESCRIPTION = "openDNA sequence file"

NONE, THIS, OTHER = "none", "this", "other"   # status of one extension: not registered / for this exe / for another path

# What the Windows Settings page asks the user to do (shown after registering).
CLICK_TEXT = ("openDNA is now offered to Windows. Windows does not let a program choose itself, so one last step is yours:\n\n"
              "In Windows Settings (Default apps), click openDNA, then, for each file type you want (.gb, .gbk, .dna, .fa, .fasta ...), "
              "click its button and choose openDNA.\n\n"
              "Or, without Settings: right-click a file > Open with > Choose another app > openDNA > \"Always\".")


# -------------------------------------------------------------------------------------------------- the registry

class WinRegistry:
    """The real registry (HKEY_CURRENT_USER only). Paths are relative to HKCU. A value of None is written as REG_NONE (an empty marker,
    which is what Windows itself writes in OpenWithProgids). Everything raises OSError on a real failure, never for a missing key."""

    def __init__(self):
        import winreg   # only exists on Windows: imported here so the module loads anywhere
        self.w = winreg

    def _open(self, path, write=False):
        access = self.w.KEY_READ | (self.w.KEY_WRITE if write else 0)
        try:
            return self.w.OpenKey(self.w.HKEY_CURRENT_USER, path, 0, access)
        except FileNotFoundError:
            return None

    def get(self, path, name=""):
        key = self._open(path)
        if key is None:
            return None
        with key:
            try:
                value, kind = self.w.QueryValueEx(key, name)
            except FileNotFoundError:
                return None
        return value if isinstance(value, str) else ""   # REG_NONE comes back as bytes

    def set(self, path, name, value):
        with self.w.CreateKeyEx(self.w.HKEY_CURRENT_USER, path, 0, self.w.KEY_WRITE) as key:   # creates the missing parents too
            if value is None:
                self.w.SetValueEx(key, name, 0, self.w.REG_NONE, b"")
            else:
                self.w.SetValueEx(key, name, 0, self.w.REG_SZ, value)

    def delete_value(self, path, name):
        key = self._open(path, write=True)
        if key is None:
            return
        with key:
            try:
                self.w.DeleteValue(key, name)
            except FileNotFoundError:
                pass

    def exists(self, path):
        key = self._open(path)
        if key is None:
            return False
        key.Close()
        return True

    def _enum(self, path, enum):
        key = self._open(path)
        if key is None:
            return []
        with key:
            names, i = [], 0
            while True:
                try:
                    item = enum(key, i)
                except OSError:
                    return names
                names.append(item if isinstance(item, str) else item[0])   # EnumValue gives (name, data, type)
                i += 1

    def subkeys(self, path):
        return self._enum(path, self.w.EnumKey)

    def values(self, path):
        return self._enum(path, self.w.EnumValue)

    def is_empty(self, path):
        """True for a key with no subkey and no value at all (a missing key counts as not empty: there is nothing to remove)."""
        key = self._open(path)
        if key is None:
            return False
        with key:
            subkeys, values, _ = self.w.QueryInfoKey(key)
        return subkeys == 0 and values == 0

    def delete_tree(self, path):
        for sub in self.subkeys(path):
            self.delete_tree(path + "\\" + sub)
        if self.exists(path):
            self.w.DeleteKey(self.w.HKEY_CURRENT_USER, path)

    def delete_key(self, path):
        """Deletes a key that has no subkey (the caller checked is_empty)."""
        if self.exists(path):
            self.w.DeleteKey(self.w.HKEY_CURRENT_USER, path)

    def changed(self):
        """Tells Explorer the associations changed (SHCNE_ASSOCCHANGED), so icons and "Open with" refresh without a log-off."""
        import ctypes
        ctypes.windll.shell32.SHChangeNotify(0x08000000, 0x0000, None, None)   # SHCNE_ASSOCCHANGED, SHCNF_IDLIST


def _registry(registry):
    return registry if registry is not None else WinRegistry()


# -------------------------------------------------------------------------------------------------- what to register

def is_frozen():
    return bool(getattr(sys, "frozen", False))


def supported():
    """Only the shipped Windows exe registers itself: from source there is no exe, and other systems have no such registry."""
    return sys.platform == "win32" and is_frozen()


def launcher(dev=False):
    """The program to start, as a list: [the exe] when frozen (sys.executable is the exe itself even for a one-file PyInstaller build,
    not the temporary unpack folder); for developers only (dev=True) [pythonw.exe, main.py]; None when there is nothing to register."""
    if is_frozen():
        return [str(Path(sys.executable).resolve())]
    if dev:
        python = Path(sys.executable)
        pythonw = python.with_name("pythonw.exe")
        return [str(pythonw if pythonw.exists() else python), str(Path(__file__).resolve().with_name("main.py"))]
    return None


def command_for(launch):
    return " ".join(f'"{part}"' for part in launch) + ' "%1"'


def _norm(command):
    return (command or "").strip().casefold()   # Windows paths do not care about case


def _exts(exts):
    chosen = list(FILE_TYPES) if exts is None else [e if e.startswith(".") else "." + e for e in exts]
    chosen = [e.lower() for e in chosen]
    unknown = [e for e in chosen if e not in FILE_TYPES]
    if unknown:
        raise ValueError(f"not a file type of {NAME}: {', '.join(unknown)}")   # nothing is written for an extension the program does not open
    return chosen


LEGACY_NAMES = ("OpenGenEdit",)   # what this program was called before openDNA: its entries show the old name and the old exe path


def has_legacy(registry=None):
    return any(_registry(registry).exists(f"{CLASSES}\\{n}.SequenceFile") for n in LEGACY_NAMES)


def remove_legacy(registry=None):
    """Deletes what the program registered under its old names. A UserChoice that named the old ProgID is left alone (Windows owns it):
    it then falls back to asking, and the user picks openDNA once."""
    reg = _registry(registry)
    for n in LEGACY_NAMES:
        progid = f"{n}.SequenceFile"
        for ext in set(FILE_TYPES) | {".gbk", ".fa"}:
            reg.delete_value(f"{CLASSES}\\{ext}\\OpenWithProgids", progid)
        for tree in (f"{CLASSES}\\{progid}", f"{CLASSES}\\Applications\\{n}.exe", rf"Software\{n}\Capabilities"):
            if reg.exists(tree):
                reg.delete_tree(tree)
        if reg.get(REGISTERED_APPS, n) is not None:
            reg.delete_value(REGISTERED_APPS, n)


def _progid_command(reg):
    return reg.get(f"{CLASSES}\\{PROGID}\\shell\\open\\command")


# -------------------------------------------------------------------------------------------------- the public functions

def status(exts=None, registry=None, launch=None):
    """{extension: NONE | THIS | OTHER}. THIS: the extension offers openDNA and the command is the one of the running exe;
    OTHER: it offers openDNA but the stored command is another path (the exe moved, or an older copy registered itself);
    NONE: not offered. `launch`: the program list to compare with (default: launcher(), the running exe)."""
    reg = _registry(registry)
    launch = launch if launch is not None else launcher()
    ours = _norm(command_for(launch)) if launch else None
    stored = _norm(_progid_command(reg))
    result = {}
    for ext in _exts(exts):
        if reg.get(f"{CLASSES}\\{ext}\\OpenWithProgids", PROGID) is None:
            result[ext] = NONE
        else:
            result[ext] = THIS if ours is not None and stored == ours else OTHER
    return result


def is_registered_for_current_exe(registry=None, launch=None):
    """True when every extension the program opens is registered for this very exe path."""
    return all(s == THIS for s in status(None, registry, launch).values())


def is_stale(registry=None, launch=None):
    """True when something is registered for another path than the running exe (the exe moved or was replaced)."""
    return OTHER in status(None, registry, launch).values()


def is_default_for(ext, registry=None):
    """True when the user has chosen openDNA for this extension in Windows (reads UserChoice; never writes it)."""
    return _registry(registry).get(USER_CHOICE.format(ext=ext), "ProgId") == PROGID


def register(exts=None, registry=None, launch=None, dev=False):
    """Writes the registration for `exts` (default: every extension of FILE_TYPES). Calling it again just refreshes the paths, so it is
    also the "the exe moved" repair. Raises RuntimeError when there is no exe to register (running from source without dev=True)
    and OSError when the registry refuses. Returns the list of extensions registered."""
    chosen = _exts(exts)
    launch = launch if launch is not None else launcher(dev)
    if not launch:
        raise RuntimeError(f"{NAME} can only register itself from the program file (.exe), not when it runs from source.")
    reg = _registry(registry)
    remove_legacy(reg)
    command = command_for(launch)
    icon = f'"{launch[0]}",0'
    progid = f"{CLASSES}\\{PROGID}"
    reg.set(progid, "", DESCRIPTION)
    reg.set(progid, "FriendlyTypeName", DESCRIPTION)
    reg.set(progid, "FriendlyAppName", NAME)
    reg.set(progid + r"\DefaultIcon", "", icon)
    reg.set(progid + r"\shell\open", "FriendlyAppName", NAME)
    reg.set(progid + r"\shell\open\command", "", command)
    app = f"{CLASSES}\\Applications\\{APP_KEY}"
    reg.set(app, "FriendlyAppName", NAME)
    reg.set(app + r"\DefaultIcon", "", icon)
    reg.set(app + r"\shell\open\command", "", command)
    reg.set(CAPABILITIES, "ApplicationName", NAME)
    reg.set(CAPABILITIES, "ApplicationDescription", f"Offline DNA sequence editor: opens {', '.join(sorted(set(FILE_TYPES.values())))}s.")
    for ext in chosen:
        reg.set(f"{CLASSES}\\{ext}\\OpenWithProgids", PROGID, None)
        reg.set(f"{app}\\SupportedTypes", ext, "")
        reg.set(CAPABILITIES + r"\FileAssociations", ext, PROGID)
    reg.set(REGISTERED_APPS, NAME, CAPABILITIES)
    reg.changed()
    return chosen


def unregister(registry=None):
    """Removes what register() created, and only that: the ProgID tree, the Applications key, the Capabilities key, our value in
    RegisteredApplications, and our value under each extension's OpenWithProgids (the key is deleted too only if it is then empty, so a
    .gb key that other programs use keeps everything of theirs). Never touches UserChoice (the user's own choice stays: Windows falls
    back to asking). Returns True if something was removed."""
    reg = _registry(registry)
    app = f"{CLASSES}\\Applications\\{APP_KEY}"
    # Extensions: those of FILE_TYPES plus any an older version registered (named in our Capabilities key).
    exts = set(FILE_TYPES) | {e for e in reg.values(CAPABILITIES + r"\FileAssociations") if e.startswith(".")}
    found = False
    for ext in sorted(exts):
        path = f"{CLASSES}\\{ext}\\OpenWithProgids"
        if reg.get(path, PROGID) is not None:
            reg.delete_value(path, PROGID)
            found = True
            for key in (path, f"{CLASSES}\\{ext}"):   # only keys that are now empty (ours alone): a shared .gb key keeps its content
                if reg.is_empty(key):
                    reg.delete_key(key)
    for tree in (f"{CLASSES}\\{PROGID}", app, CAPABILITIES):
        if reg.exists(tree):
            reg.delete_tree(tree)
            found = True
    if reg.get(REGISTERED_APPS, NAME) == CAPABILITIES:   # a shared Windows key: only our own value, and only if it is still ours
        reg.delete_value(REGISTERED_APPS, NAME)
        found = True
    if found:
        reg.changed()
    return found


def open_default_apps_settings(opener=None):
    """Opens Windows Settings at this app's Default apps page (Windows 11: its own page; older: the general page, where the user
    searches for openDNA). `opener(uri)` is injectable for tests. Returns True if something was opened."""
    opener = opener or getattr(os, "startfile", None)
    if opener is None:
        return False
    for uri in (SETTINGS_APP_URI, SETTINGS_URI):
        try:
            opener(uri)
            return True
        except OSError:
            continue
    return False


# -------------------------------------------------------------------------------------------------- memory of the user's choice

def settings_path():
    return str(data_dir() / "file_assoc.ini")


def open_settings(path=None):
    return QSettings(path or settings_path(), QSettings.IniFormat)


def _flag(settings, key):
    return str(settings.value(key, "false")).lower() == "true"   # an ini file hands booleans back as text


def should_ask(settings=None, registry=None, launch=None, windows=None, frozen=None):
    """Whether to put the question to the user now: None (no), "first" (never asked, nothing registered) or "stale" (it was registered,
    the user said Yes, and the exe has since moved). The table:

        not Windows, or not the frozen exe        -> None   (nothing to register, nothing to ask)
        "Don't ask again" was chosen              -> None
        already registered for this exe           -> None
        some extension registered for another path
            user said Yes before, not yet asked for this path   -> "stale"       (any other case: None)
        nothing registered
            never asked                           -> "first"    (a Not now / Yes / Don't ask again is an answer: asked once only)
    """
    windows = (sys.platform == "win32") if windows is None else windows
    frozen = is_frozen() if frozen is None else frozen
    if not (windows and frozen):
        return None
    settings = open_settings() if settings is None else settings
    if _flag(settings, "dont_ask"):
        return None
    launch = launch if launch is not None else (launcher() or [str(Path(sys.executable).resolve())])
    states = set(status(None, registry, launch).values())
    if has_legacy(registry):   # an older name is still registered: refresh once per path, whatever the answer was before
        return "stale" if _norm(str(settings.value("stale_asked_for", ""))) != _norm(command_for(launch)) else None
    if states == {THIS}:
        return None
    if OTHER in states:
        stale_for = str(settings.value("stale_asked_for", ""))
        return "stale" if _flag(settings, "accepted") and _norm(stale_for) != _norm(command_for(launch)) else None
    return None if _flag(settings, "asked") else "first"


# -------------------------------------------------------------------------------------------------- the screens

_open_boxes = []   # message boxes shown without exec(): kept referenced until closed (tests read them from here)


def _notice(parent, title, text, icon=QMessageBox.Information, settings_button=False):
    """A message box that does not block (open(), not exec()): the program keeps running behind it and a test can click its buttons.
    settings_button: adds "Open Windows Settings" (which calls open_default_apps_settings())."""
    box = QMessageBox(icon, title, text, QMessageBox.NoButton, parent)
    box.setWindowModality(Qt.WindowModal if parent is not None else Qt.NonModal)
    box.setTextInteractionFlags(Qt.TextSelectableByMouse)
    box.open_button = box.addButton("Open Windows Settings", QMessageBox.AcceptRole) if settings_button else None
    box.close_button = box.addButton("Close" if settings_button else "OK", QMessageBox.RejectRole)
    if box.open_button is not None:
        box.open_button.clicked.connect(lambda: open_default_apps_settings())
    _open_boxes.append(box)
    box.finished.connect(lambda _=0, b=box: _open_boxes.remove(b) if b in _open_boxes else None)
    box.open()
    return box


def register_now(parent=None, registry=None, settings=None, launch=None, dev=False):
    """Registers every extension, remembers that the user wants it (so a moved exe is offered a repair later), and shows what to click
    in Windows Settings. Returns True on success; on failure shows why and returns False."""
    try:
        register(None, registry, launch, dev)
    except (OSError, RuntimeError) as e:
        _notice(parent, "File associations", f"Could not register {NAME}:\n\n{e}", QMessageBox.Warning)
        return False
    settings = open_settings() if settings is None else settings
    settings.setValue("accepted", "true")
    settings.setValue("asked", "true")
    settings.sync()
    _notice(parent, "File associations", CLICK_TEXT, settings_button=True)
    return True


def unregister_now(parent=None, registry=None, settings=None):
    try:
        removed = unregister(registry)
    except OSError as e:
        _notice(parent, "File associations", f"Could not remove the registration:\n\n{e}", QMessageBox.Warning)
        return False
    settings = open_settings() if settings is None else settings
    settings.setValue("accepted", "false")   # removed on purpose: never offer a "repair" for it
    settings.sync()
    _notice(parent, "File associations",
            f"{NAME} is no longer offered for these files." if removed else f"{NAME} was not registered.")
    return True


def ask_at_startup(parent=None, settings=None, registry=None, launch=None, windows=None, frozen=None):
    """The one-time question. Returns the (open, non-blocking) box, or None when nothing is to be asked (see should_ask). Buttons are
    attributes (yes_button, later_button, never_button) so a test clicks them; the answer is remembered in file_assoc.ini."""
    settings = open_settings() if settings is None else settings
    kind = should_ask(settings, registry, launch, windows, frozen)
    if kind is None:
        return None
    if kind == "stale":   # the user already said Yes once: follow the exe to its new place without asking again
        settings.setValue("stale_asked_for", command_for(launch if launch is not None else launcher() or [sys.executable]))
        try:
            register(None, registry, launch)
        except (OSError, RuntimeError):
            pass   # shortcut: a refusal is silent here; Help > File Associations shows the error when the user retries
        settings.sync()
        return None
    if kind == "first":
        text = (f"Open .gb, .dna and FASTA files with {NAME} by default?\n\n"
                f"This adds {NAME} to the \"Open with\" list of these file types (current user only, no administrator rights). "
                "Windows will then ask you to confirm in its Settings. You can change your mind in Help > File Associations.")
    else:
        text = (f"{NAME} was set up to open sequence files, but the program file has moved since.\n\n"
                "Update the registration so that double-clicking a file starts it from the new place?")
    box = QMessageBox(QMessageBox.Question, "File associations", text, QMessageBox.NoButton, parent)
    box.yes_button = box.addButton("Yes", QMessageBox.YesRole)
    box.later_button = box.addButton("Not now", QMessageBox.RejectRole)
    box.never_button = box.addButton("Don't ask again", QMessageBox.DestructiveRole)
    box.setDefaultButton(box.yes_button)
    box.setEscapeButton(box.later_button)

    def answered(button):
        if button is box.yes_button:
            register_now(parent, registry, settings, launch)
        elif button is box.never_button:
            settings.setValue("dont_ask", "true")
        settings.setValue("asked", "true")
        if kind == "stale":   # whatever was answered, do not repeat the question for this same path
            settings.setValue("stale_asked_for", command_for(launch if launch is not None else launcher() or [sys.executable]))
        settings.sync()
        if box in _open_boxes:
            _open_boxes.remove(box)

    box.buttonClicked.connect(answered)
    box.setWindowModality(Qt.WindowModal if parent is not None else Qt.NonModal)
    _open_boxes.append(box)
    box.open()
    return box


class FileAssociationsDialog(QDialog):
    """Help > File Associations…: says which extensions are registered and registers / updates / removes them later.
    Meant to be shown with open() (or show()); no exec() needed, and the buttons are attributes for tests."""

    def __init__(self, parent=None, registry=None, settings=None, launch=None, dev=False):
        super().__init__(parent)
        self.registry, self.settings, self.launch, self.dev = registry, settings, launch, dev
        self.setWindowTitle("File associations")
        self.text = QLabel(wordWrap=True, textInteractionFlags=Qt.TextSelectableByMouse)
        self.text.setMinimumWidth(420)
        self.register_button = QPushButton("Register")
        self.remove_button = QPushButton("Remove registration")
        self.settings_button = QPushButton("Open Windows Settings")
        self.close_button = QPushButton("Close")
        row = QHBoxLayout()
        for b in (self.register_button, self.remove_button, self.settings_button):
            row.addWidget(b)
        row.addStretch(1)
        row.addWidget(self.close_button)
        layout = QVBoxLayout(self)
        layout.addWidget(self.text)
        layout.addLayout(row)
        self.register_button.clicked.connect(self._register)
        self.remove_button.clicked.connect(self._remove)
        self.settings_button.clicked.connect(lambda: open_default_apps_settings())
        self.close_button.clicked.connect(self.accept)
        self.refresh()

    def _launch(self):
        return self.launch if self.launch is not None else launcher(self.dev)

    def _states(self):
        return status(None, self.registry, self._launch()) if self._launch() else {}

    def refresh(self):
        states = self._states()
        values = set(states.values())
        if not self._launch():
            msg = (f"{NAME} is running from source, so there is no program file to register. "
                   "The program file (.exe) can register itself from here.")
        elif values == {THIS}:
            msg = f"{NAME} is registered for this program file and offered for: {', '.join(states)}."
        elif OTHER in values:
            msg = (f"{NAME} is registered, but for another location of the program file (it was moved or replaced). "
                   "Click Update so that double-clicking a file starts this one.")
        else:
            msg = f"{NAME} is not yet offered to Windows for {', '.join(FILE_TYPES)} files."
        if self._launch() and THIS in values and OTHER not in values:
            chosen = [e for e in states if is_default_for(e, self.registry)]
            msg += (f"\n\nDefault for: {', '.join(chosen)}." if chosen else
                    "\n\nIt is not the default yet: choose it in Windows Settings (button below).")
        self.text.setText(msg)
        self.register_button.setText("Update" if OTHER in values else "Register")
        self.register_button.setEnabled(bool(self._launch()) and values != {THIS})
        self.remove_button.setEnabled(bool(values - {NONE}))

    def _register(self):
        register_now(self, self.registry, self.settings, self._launch(), self.dev)
        self.refresh()

    def _remove(self):
        unregister_now(self, self.registry, self.settings)
        self.refresh()


def show_manage_dialog(parent=None, registry=None, settings=None, launch=None, dev=False):
    """For the Help menu action: opens the dialog (non-blocking) and returns it."""
    dlg = FileAssociationsDialog(parent, registry, settings, launch, dev)
    dlg.setAttribute(Qt.WA_DeleteOnClose, parent is not None)
    dlg.open()
    return dlg
