"""File associations (app/file_assoc.py): the registration, status, stale-path detection, unregister, and the start-up question.
Never touches the real registry (a dictionary-backed fake stands in for it) nor the user's settings (private ini files)."""
import os
import sys
import tempfile
from pathlib import Path

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ.setdefault("GENOMEEDIT_DATA_DIR", tempfile.mkdtemp(prefix="genomeedit_test_"))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "app"))

from PySide6.QtWidgets import QApplication

import file_assoc as F
from project import FILE_TYPES

app = QApplication.instance() or QApplication([])

EXE = [r"C:\Apps\OpenGenEdit\OpenGenEdit.exe"]
MOVED = [r"D:\Elsewhere\OpenGenEdit.exe"]


class FakeRegistry:
    """The small surface of WinRegistry on a dict: key path (case-insensitive, like Windows) -> {value name: data}."""

    def __init__(self):
        self.keys = {}
        self.changes = 0
        self.writes = []   # every (path, name) written: to check nothing outside our keys was touched

    @staticmethod
    def _k(path):
        return path.casefold()

    def _mk(self, path):   # creates the parents too, like CreateKeyEx
        parts = path.split("\\")
        for i in range(1, len(parts) + 1):
            self.keys.setdefault(self._k("\\".join(parts[:i])), {})

    def exists(self, path):
        return self._k(path) in self.keys

    def get(self, path, name=""):
        return self.keys.get(self._k(path), {}).get(name.casefold())

    def set(self, path, name, value):
        self._mk(path)
        self.keys[self._k(path)][name.casefold()] = "" if value is None else value
        self.writes.append((path, name))

    def delete_value(self, path, name):
        self.keys.get(self._k(path), {}).pop(name.casefold(), None)

    def subkeys(self, path):
        prefix = self._k(path) + "\\"
        return sorted({k[len(prefix):].split("\\")[0] for k in self.keys if k.startswith(prefix)})

    def values(self, path):
        return list(self.keys.get(self._k(path), {}))

    def is_empty(self, path):
        return self.exists(path) and not self.keys[self._k(path)] and not self.subkeys(path)

    def delete_tree(self, path):
        base = self._k(path)
        for k in [k for k in self.keys if k == base or k.startswith(base + "\\")]:
            del self.keys[k]

    def delete_key(self, path):
        assert not self.subkeys(path), "delete_key on a key with subkeys"
        self.keys.pop(self._k(path), None)

    def changed(self):
        self.changes += 1

    def snapshot(self):
        return {k: dict(v) for k, v in self.keys.items()}


@pytest.fixture
def reg():
    return FakeRegistry()


@pytest.fixture
def ini(tmp_path):
    return F.open_settings(str(tmp_path / "file_assoc.ini"))


def close_open_boxes():
    for box in list(F._open_boxes):   # close what a test left open
        try:
            box.close()
        except RuntimeError:   # already deleted together with its parent dialog
            pass
    F._open_boxes.clear()


@pytest.fixture(autouse=True)
def no_boxes():
    yield
    close_open_boxes()


def test_register(reg, monkeypatch):
    """register: everything it writes, a subset and an unknown extension, only our keys touched, refusing from source, the
    frozen launcher."""
    # --- writes everything
    done = F.register(None, reg, EXE)
    assert done == list(FILE_TYPES)
    pid = r"Software\Classes\OpenGenEdit.SequenceFile"
    assert reg.get(pid + r"\shell\open\command") == '"C:\\Apps\\OpenGenEdit\\OpenGenEdit.exe" "%1"'
    assert reg.get(pid + r"\DefaultIcon") == '"C:\\Apps\\OpenGenEdit\\OpenGenEdit.exe",0'
    assert reg.get(pid, "FriendlyAppName") == "OpenGenEdit"
    for ext in FILE_TYPES:
        assert reg.get(rf"Software\Classes\{ext}\OpenWithProgids", F.PROGID) is not None
        assert reg.get(r"Software\Classes\Applications\OpenGenEdit.exe\SupportedTypes", ext) is not None
        assert reg.get(F.CAPABILITIES + r"\FileAssociations", ext) == F.PROGID
    assert reg.get(F.REGISTERED_APPS, "OpenGenEdit") == F.CAPABILITIES
    assert reg.changes == 1   # the shell was told once

    # --- only writes our keys
    ours = (r"Software\Classes\OpenGenEdit.SequenceFile", r"Software\Classes\Applications\OpenGenEdit.exe", F.CAPABILITIES, F.REGISTERED_APPS)
    for path, name in reg.writes:
        assert path.startswith(ours) or path.endswith("OpenWithProgids"), path
        if path == F.REGISTERED_APPS:
            assert name == "OpenGenEdit"
        if path.endswith("OpenWithProgids"):
            assert name == F.PROGID
    assert not any("userchoice" in p.casefold() for p, _ in reg.writes)

    # --- a subset and an unknown extension
    reg = FakeRegistry()
    assert F.register([".gb", "fa"], reg, EXE) == [".gb", ".fa"]
    st = F.status(None, reg, EXE)
    assert st[".gb"] == F.THIS and st[".dna"] == F.NONE
    assert not F.is_registered_for_current_exe(reg, EXE)
    before = reg.snapshot()
    with pytest.raises(ValueError):
        F.register([".txt"], reg, EXE)
    assert reg.snapshot() == before   # nothing written for a type the program does not open

    # --- from source it refuses politely
    reg = FakeRegistry()
    with monkeypatch.context() as source:
        source.setattr(F, "is_frozen", lambda: False)
        assert F.launcher() is None
        with pytest.raises(RuntimeError, match="program file"):
            F.register(None, reg)
        assert reg.keys == {}
        dev = F.launcher(dev=True)   # developers only: python(w) + main.py
        assert dev[1].endswith("main.py") and F.command_for(dev).endswith('"%1"')

    # --- the frozen launcher is the running exe
    with monkeypatch.context() as frozen:
        frozen.setattr(sys, "frozen", True, raising=False)
        assert F.launcher() == [str(Path(sys.executable).resolve())]


def test_status_default_and_settings_page(reg):
    """status round trip and stale detection (path case ignored), the default-app check reading UserChoice, and opening the
    Windows default-apps page with its fallback."""
    # --- status round trip and stale
    assert set(F.status(None, reg, EXE).values()) == {F.NONE}
    assert not F.is_registered_for_current_exe(reg, EXE) and not F.is_stale(reg, EXE)
    F.register(None, reg, EXE)
    assert set(F.status(None, reg, EXE).values()) == {F.THIS}
    assert F.is_registered_for_current_exe(reg, EXE) and not F.is_stale(reg, EXE)
    # the exe moved: same registry, another running path
    assert set(F.status(None, reg, MOVED).values()) == {F.OTHER}
    assert F.is_stale(reg, MOVED) and not F.is_registered_for_current_exe(reg, MOVED)
    F.register(None, reg, MOVED)   # the repair is just registering again
    assert F.is_registered_for_current_exe(reg, MOVED)
    assert reg.get(r"Software\Classes\OpenGenEdit.SequenceFile\shell\open\command").startswith('"D:')

    # --- path case is ignored
    reg = FakeRegistry()
    F.register(None, reg, EXE)
    assert F.is_registered_for_current_exe(reg, [EXE[0].upper()])

    # --- is_default reads UserChoice
    reg = FakeRegistry()
    assert not F.is_default_for(".gb", reg)
    reg.set(F.USER_CHOICE.format(ext=".gb"), "ProgId", F.PROGID)
    assert F.is_default_for(".gb", reg) and not F.is_default_for(".fa", reg)

    # --- the default-apps settings page
    opened = []
    assert F.open_default_apps_settings(opened.append) is True
    assert opened == ["ms-settings:defaultapps?registeredAppUser=OpenGenEdit"]

    def refuses_deep_link(uri):
        if "?" in uri:
            raise OSError("no handler")
        opened.append(uri)
    opened.clear()
    assert F.open_default_apps_settings(refuses_deep_link) is True
    assert opened == ["ms-settings:defaultapps"]   # falls back to the general page

    def never(uri):
        raise OSError
    assert F.open_default_apps_settings(never) is False


def empty_skeleton(reg):
    """Nothing but empty parent keys is left (Software, Software/Classes...: shared and harmless, never removed)."""
    return not any(reg.keys.values()) and not any(F.PROGID.casefold() in k or "capabilities" in k or "openwithprogids" in k
                                                 or k.startswith(r"software\classes\.") for k in reg.keys)


def test_unregister(reg):
    """unregister removes only our keys (exactly the state before), leaves nothing of ours behind even on an empty registry, and
    removes the extensions of an older version."""
    # --- removes only ours
    # Other programs' keys that must survive: a shared .gb key with their own ProgID, an unrelated app, our QSettings area.
    reg.set(r"Software\Classes\.gb\OpenWithProgids", "Other.Program", None)
    reg.set(r"Software\Classes\.gb", "", "Other.Program")
    reg.set(r"Software\Classes\Applications\Other.exe\shell\open\command", "", '"other.exe" "%1"')
    reg.set(r"Software\OpenGenEdit\OpenGenEdit\view", "zoom", "3")
    reg.set(F.REGISTERED_APPS, "OtherApp", r"Software\Other\Capabilities")
    reg.set(F.USER_CHOICE.format(ext=".gb"), "ProgId", "Other.Program")
    others = reg.snapshot()
    F.register(None, reg, EXE)
    assert F.unregister(reg) is True
    assert reg.snapshot() == others   # exactly the state before: not a key more, not a key less
    assert F.status(None, reg, EXE) == {e: F.NONE for e in FILE_TYPES}
    assert F.unregister(reg) is False   # nothing left to remove: no error, no notification
    assert reg.changes == 2   # register + the first unregister

    # --- removes the keys it created even if nothing else existed
    reg = FakeRegistry()
    F.register(None, reg, EXE)
    F.unregister(reg)
    assert empty_skeleton(reg), reg.snapshot()

    # --- removes extensions of an older version
    reg = FakeRegistry()
    F.register(None, reg, EXE)
    reg.set(F.CAPABILITIES + r"\FileAssociations", ".old", F.PROGID)
    reg.set(r"Software\Classes\.old\OpenWithProgids", F.PROGID, None)
    F.unregister(reg)
    assert not reg.exists(r"Software\Classes\.old")


# ---------------------------------------------------------------------------------------- the start-up question: decision table

def ask(reg, ini, *, windows=True, frozen=True, launch=EXE):
    return F.should_ask(ini, reg, launch, windows, frozen)


def test_should_ask_decision_table(reg, ini, tmp_path):
    """Never asks from source or off Windows, nothing when already registered, "don't ask again" wins, asks once, a stale
    registration asks only after a yes (and again for a new move), and the ini booleans survive a reload."""
    # --- never from source or off Windows
    assert ask(reg, ini, frozen=False) is None
    assert ask(reg, ini, windows=False) is None
    assert ask(reg, ini, windows=False, frozen=False) is None
    assert ask(reg, ini) == "first"

    # --- nothing when already registered
    reg = FakeRegistry()
    ini.clear()
    F.register(None, reg, EXE)
    assert ask(reg, ini) is None

    # --- "don't ask again" wins over everything
    reg = FakeRegistry()
    ini.clear()
    ini.setValue("dont_ask", "true")
    assert ask(reg, ini) is None
    ini.setValue("accepted", "true")
    F.register(None, reg, MOVED)   # stale and accepted, but "don't ask again"
    assert ask(reg, ini) is None

    # --- asks once
    reg = FakeRegistry()
    ini.clear()
    assert ask(reg, ini) == "first"
    ini.setValue("asked", "true")
    assert ask(reg, ini) is None

    # --- stale asks only after a yes
    reg = FakeRegistry()
    ini.clear()
    F.register(None, reg, MOVED)   # registered for another path
    assert ask(reg, ini) is None   # the user never said Yes (the registration came from elsewhere): no nagging
    ini.setValue("accepted", "true")
    assert ask(reg, ini) == "stale"
    ini.setValue("stale_asked_for", F.command_for(EXE))   # answered for this path already
    assert ask(reg, ini) is None
    assert ask(reg, ini, launch=[r"E:\Third\OpenGenEdit.exe"]) == "stale"   # but a new move asks again

    # --- ini booleans survive a reload
    reg = FakeRegistry()
    path = str(tmp_path / "a.ini")
    s = F.open_settings(path)
    s.setValue("asked", "true")
    s.sync()
    assert ask(reg, F.open_settings(path)) is None


# ---------------------------------------------------------------------------------------- the screens (driven without exec)

def quiet_settings_button(monkeypatch):
    opened = []
    monkeypatch.setattr(F, "open_default_apps_settings", lambda opener=None: opened.append(1) or True)
    return opened


def test_startup_question_screens(reg, ini, monkeypatch):
    """The start-up question driven without exec: nothing to ask, Yes registers/remembers/shows what to click, Not now and
    Don't ask again, the stale question and its repair, and a failing registration shown instead of raised."""
    # --- nothing to ask
    assert F.ask_at_startup(None, ini, reg, EXE, windows=True, frozen=False) is None
    assert F.ask_at_startup(None, ini, reg, EXE, windows=False, frozen=True) is None
    assert F._open_boxes == []

    # --- Yes registers, remembers and shows what to click
    opened = quiet_settings_button(monkeypatch)
    box = F.ask_at_startup(None, ini, reg, EXE, windows=True, frozen=True)
    assert box is not None and "FASTA" in box.text()
    assert [b.text() for b in (box.yes_button, box.later_button, box.never_button)] == ["Yes", "Not now", "Don't ask again"]
    box.yes_button.click()
    assert F.is_registered_for_current_exe(reg, EXE)
    assert ask(reg, ini) is None and ini.value("accepted") == "true"
    info = [b for b in F._open_boxes if "Default apps" in b.text()]
    assert len(info) == 1 and "click OpenGenEdit" in info[0].text()
    info[0].open_button.click()   # "Open Windows Settings"
    assert opened == [1]

    # --- Not now and Don't ask again
    close_open_boxes()
    reg = FakeRegistry()
    ini.clear()
    F.ask_at_startup(None, ini, reg, EXE, windows=True, frozen=True).later_button.click()
    assert F.status(None, reg, EXE)[".gb"] == F.NONE   # nothing written
    assert ask(reg, ini) is None and ini.value("dont_ask") is None   # asked once; not a permanent refusal
    ini.clear()
    F.ask_at_startup(None, ini, reg, EXE, windows=True, frozen=True).never_button.click()
    assert ask(reg, ini) is None and ini.value("dont_ask") == "true" and reg.keys == {}

    # --- the stale question and repair
    close_open_boxes()
    reg = FakeRegistry()
    ini.clear()
    ini.setValue("accepted", "true")
    F.register(None, reg, MOVED)
    box = F.ask_at_startup(None, ini, reg, EXE, windows=True, frozen=True)
    assert "moved" in box.text()
    box.yes_button.click()
    assert F.is_registered_for_current_exe(reg, EXE)
    # declining a repair is remembered for that path
    third = [r"E:\Third\OpenGenEdit.exe"]
    F.register(None, reg, MOVED)
    F.ask_at_startup(None, ini, reg, third, windows=True, frozen=True).later_button.click()
    assert ask(reg, ini, launch=third) is None

    # --- a failed registration is shown, not raised
    close_open_boxes()
    ini.clear()

    class Refusing(FakeRegistry):
        def set(self, *a):
            raise PermissionError("denied by policy")
    assert F.register_now(None, Refusing(), ini, EXE) is False
    assert any("denied by policy" in b.text() for b in F._open_boxes)
    assert ini.value("accepted") is None   # a failed registration is not remembered as accepted


def test_manage_dialog(reg, ini, monkeypatch):
    """The management dialog (register, update after a move, remove; from source it is disabled) and the default settings file
    in the data folder."""
    # --- register, update, remove
    quiet_settings_button(monkeypatch)
    dlg = F.FileAssociationsDialog(None, reg, ini, EXE)
    assert "not yet offered" in dlg.text.text()
    assert dlg.register_button.isEnabled() and not dlg.remove_button.isEnabled()
    dlg.register_button.click()
    assert F.is_registered_for_current_exe(reg, EXE)
    assert "registered for this program file" in dlg.text.text()
    assert not dlg.register_button.isEnabled() and dlg.remove_button.isEnabled()
    dlg.launch = MOVED   # the exe moved
    dlg.refresh()
    assert dlg.register_button.text() == "Update" and dlg.register_button.isEnabled()
    dlg.remove_button.click()
    assert empty_skeleton(reg) and ini.value("accepted") == "false"
    assert not dlg.remove_button.isEnabled()

    # --- from source
    monkeypatch.setattr(F, "is_frozen", lambda: False)
    dlg = F.FileAssociationsDialog(None, FakeRegistry(), ini)
    assert "from source" in dlg.text.text()
    assert not dlg.register_button.isEnabled() and not dlg.remove_button.isEnabled()

    # --- the default settings file is in the data dir
    assert Path(F.settings_path()).parent == F.data_dir()   # (only the path: nothing is created here)
