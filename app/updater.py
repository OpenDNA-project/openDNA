"""Self-update of the single openDNA.exe: look for a newer GitHub release, download it, check it, swap it in and restart.

The steps, each usable without Qt (the dialog at the end only wires them together):
  check()            ask GitHub for the latest release (never raises: the answer says what went wrong).
  download_update()  fetch openDNA.exe next to the running exe, verify its SHA-256 against openDNA.exe.sha256, refuse otherwise.
  install_update()   rename the running exe to openDNA.exe.old (Windows allows that, not an overwrite), move the new one in, relaunch.
  cleanup_old()      at the next start-up, delete the leftover .old.

What the checksum does and does not give (be honest with the user about it): the .sha256 file comes from the same release as the exe,
so it catches a corrupted or truncated download and a bad mirror or proxy, but NOT a release an attacker has taken over (they would
replace both files). Only code signing (Authenticode) fixes that. What is enforced here: HTTPS only, and only GitHub hosts (every
redirect hop is checked too); nothing is ever executed except the exe whose hash matched.

Nothing here runs unless the caller asks. The automatic start-up check is also off while project.UPDATES_ENABLED is False (the repository
is private, so its releases cannot be read yet); the manual "Check for Updates..." says so instead.
"""
import gc
import html
import json
import os
import re
import subprocess
import sys
import threading
import time
import urllib.error
import urllib.request
from dataclasses import dataclass
from hmac import compare_digest
from pathlib import Path
from urllib.parse import urlparse

from PySide6.QtCore import QObject, QSettings, Qt, Signal
from PySide6.QtWidgets import (
    QApplication, QDialog, QHBoxLayout, QLabel, QMessageBox, QPlainTextEdit, QProgressBar, QPushButton, QVBoxLayout,
)

import project
import tool_settings
from primer_db import data_dir

ASSET_NAME = f"{project.EXE_NAME}.exe"        # what the release must contain (the exe under its released name, project.EXE_NAME) ...
CHECKSUM_NAME = f"{ASSET_NAME}.sha256"         # ... and its checksum ("<hex>  oDNA_beta.exe", like sha256sum, or just the hex)
# github.com serves the release assets and redirects to the other two (GitHub has used both names for its asset storage).
ALLOWED_HOSTS = ("api.github.com", "github.com", "objects.githubusercontent.com", "release-assets.githubusercontent.com")
USER_AGENT = f"{project.NAME}/{project.VERSION}"   # the GitHub API refuses requests without one
TIMEOUT = 30                       # seconds per network operation (not for the whole download)
MAX_TEXT = 2_000_000               # bytes: the release JSON and the checksum are small; anything bigger is not what we asked for
MAX_EXE = 1_000_000_000            # bytes: stops a wrong URL from filling the disk
CHECK_EVERY = 24 * 3600            # seconds between two automatic checks
TOOL = "Updates"                   # its tab in Tools > Settings
tool_settings.register(TOOL, [
    tool_settings.Param("check", "Look for a new version when openDNA starts (once a day)", True,
                        tip="Asks GitHub for the latest release, nothing else is sent. Tools > Check for Updates always works."),
])
RETRIES = 5                        # antivirus scanners often hold a freshly written exe for a second or two
RETRY_DELAY = 0.6                  # seconds (tests set it to 0)
CHUNK = 256 * 1024


class UpdateError(Exception):
    """Something the user can read: why the update did not happen."""


@dataclass(frozen=True)
class Update:
    version: str        # "0.2.0" (no "v")
    tag: str            # "v0.2.0" as published
    asset_url: str      # where openDNA.exe is
    sha256_url: str     # where openDNA.exe.sha256 is
    notes: str          # the release description (plain text)
    page_url: str       # the release page, for people who prefer to download by hand


@dataclass
class Check:
    """The outcome of a check: `update` when there is a newer installable version; else `error` says why not (None = up to date)."""
    update: Update = None
    error: str = None
    latest: str = ""    # the tag found, when one was
    page_url: str = ""


# ---- versions ------------------------------------------------------------------------------------------------------------------------

def parse_version(text):
    """(0, 10, 0) for "v0.10.0" / "0.10.0"; None for anything else, notably "0.2.0-beta1" (a prerelease) or "latest"."""
    m = re.fullmatch(r"[vV]?(\d+(?:\.\d+){0,3})", str(text or "").strip())
    return tuple(int(n) for n in m.group(1).split(".")) if m else None


def is_newer(candidate, current):
    """True when `candidate` is a higher version than `current` (numbers compared one by one: 0.10.0 > 0.9.0; 1.0 == 1.0.0)."""
    a, b = parse_version(candidate), parse_version(current)
    if a is None or b is None:
        return False
    size = max(len(a), len(b))
    return a + (0,) * (size - len(a)) > b + (0,) * (size - len(b))


# ---- the network ---------------------------------------------------------------------------------------------------------------------

def check_url(url):
    """Raise UpdateError unless `url` is HTTPS on one of ALLOWED_HOSTS (no login part, no odd port)."""
    p = urlparse(str(url or ""))
    try:
        port = p.port
    except ValueError:
        port = -1
    if p.scheme != "https" or (p.hostname or "").lower() not in ALLOWED_HOSTS or p.username or p.password or port not in (None, 443):
        raise UpdateError(f"Refused: {url or 'an empty address'} is not an https address on github.com.")
    return url


class _SafeRedirects(urllib.request.HTTPRedirectHandler):
    """A release asset is answered with a redirect to GitHub's storage: every hop must pass check_url too."""

    def redirect_request(self, req, fp, code, msg, headers, newurl):
        check_url(newurl)
        return super().redirect_request(req, fp, code, msg, headers, newurl)


def fetch(url, timeout=TIMEOUT, sink=None, progress=None, cancel=None):
    """The bytes at `url`; or, when `sink` (a binary file object) is given, write them there in chunks and return None.
    progress(done, total) is called after each chunk (total 0 = unknown); `cancel` is a threading.Event that stops a download.
    UpdateError when the address is not allowed, cannot be reached or says no. Tests replace this function by a fake."""
    check_url(url)
    headers = {"User-Agent": USER_AGENT}
    if urlparse(url).hostname == "api.github.com":
        headers["Accept"] = "application/vnd.github+json"
    host = urlparse(url).hostname
    try:
        with urllib.request.build_opener(_SafeRedirects).open(urllib.request.Request(url, headers=headers), timeout=timeout) as reply:
            total = int(reply.headers.get("Content-Length") or 0)
            limit = MAX_EXE if sink is not None else MAX_TEXT
            if total > limit:
                raise UpdateError(f"{host} sent more data than expected ({total:,} bytes): refused.")
            done, parts = 0, []
            while True:
                if cancel is not None and cancel.is_set():
                    raise UpdateError("Download cancelled.")
                chunk = reply.read(CHUNK)
                if not chunk:
                    break
                done += len(chunk)
                if done > limit:
                    raise UpdateError(f"{host} sent more data than expected: refused.")
                if sink is not None:
                    sink.write(chunk)
                else:
                    parts.append(chunk)
                if progress is not None:
                    progress(done, total)
            if total and done != total:
                raise UpdateError(f"The download from {host} stopped early ({done:,} of {total:,} bytes).")
    except urllib.error.HTTPError as e:
        raise UpdateError(f"{host} answered {e.code} {e.reason}.") from None
    except (urllib.error.URLError, TimeoutError, OSError) as e:
        raise UpdateError(f"Could not reach {host}: {getattr(e, 'reason', e)}") from None
    return None if sink is not None else b"".join(parts)


def _pick(f):
    """The injected fetch, else the module's fetch (looked up now, so a test that replaces updater.fetch is obeyed)."""
    return f if f is not None else fetch


# ---- looking for a release -----------------------------------------------------------------------------------------------------------

def check(current_version=None, fetch=None):
    """Ask GitHub for the latest release and compare. Never raises: see Check."""
    current = current_version or project.VERSION
    try:
        data = json.loads(_pick(fetch)(project.LATEST_RELEASE_API, timeout=TIMEOUT))
        if not isinstance(data, dict):
            raise ValueError
    except UpdateError as e:
        return Check(error=str(e))
    except ValueError:      # includes json.JSONDecodeError and UnicodeDecodeError
        return Check(error="GitHub's answer could not be understood.")
    except Exception as e:  # a fake or an odd library error: still no crash at start-up
        return Check(error=f"The update check failed: {e}")
    tag = str(data.get("tag_name") or "")
    page = data.get("html_url")
    try:
        page = check_url(page)
    except UpdateError:
        page = project.RELEASES_URL
    if data.get("draft") or data.get("prerelease"):
        return Check(latest=tag, page_url=page)               # not for everyone: as if it did not exist
    version = parse_version(tag)
    if version is None:
        return Check(error=f"The latest release has a version name ({tag!r}) that is not understood.", latest=tag, page_url=page)
    if not is_newer(tag, current):
        return Check(latest=tag, page_url=page)
    assets = {a.get("name"): a.get("browser_download_url") for a in data.get("assets") or [] if isinstance(a, dict)}
    missing = [n for n in (ASSET_NAME, CHECKSUM_NAME) if not assets.get(n)]
    if missing:
        return Check(error=f"Version {tag.lstrip('vV')} exists but its release has no {' / '.join(missing)}, so it cannot be installed "
                           "automatically.", latest=tag, page_url=page)
    try:
        check_url(assets[ASSET_NAME]), check_url(assets[CHECKSUM_NAME])
        # The exe must come from this project's own releases, not just from any github.com address.
        if not urlparse(assets[ASSET_NAME]).path.startswith(f"/{project.REPO}/releases/download/"):
            raise UpdateError(f"Refused: {assets[ASSET_NAME]} is not in the releases of {project.REPO}.")
    except UpdateError as e:
        return Check(error=str(e), latest=tag, page_url=page)
    notes = data.get("body") if isinstance(data.get("body"), str) else ""
    update = Update(".".join(map(str, version)), tag, assets[ASSET_NAME], assets[CHECKSUM_NAME], notes[:20000], page)
    return Check(update=update, latest=tag, page_url=page)


def check_for_update(current_version=None, fetch=None):
    """The Update when a newer release can be installed, else None (up to date, or no network: use check() to know which)."""
    return check(current_version, fetch).update


# ---- downloading ---------------------------------------------------------------------------------------------------------------------

def parse_checksum(text, filename=ASSET_NAME):
    """The 64-hex digest for `filename` from a .sha256 file (bytes or str), or None. Accepts "<hex>  name", "<hex> *name" and a bare
    "<hex>". PowerShell's `>` writes UTF-16, so a byte order mark is honoured."""
    if isinstance(text, bytes):
        text = text.decode("utf-16" if text[:2] in (b"\xff\xfe", b"\xfe\xff") else "utf-8-sig", errors="replace")
    for line in text.splitlines():
        parts = line.split()
        if parts and re.fullmatch(r"[0-9a-fA-F]{64}", parts[0]) and (len(parts) == 1 or parts[1].lstrip("*") == filename):
            return parts[0].lower()
    return None


def sha256_of(path):
    import hashlib
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def download_update(update, folder, progress=None, fetch=None, cancel=None):
    """Download the new exe into `folder` (put it next to the running exe: same drive, so the final move is a rename), verify it, and
    return its path (openDNA.exe.new). UpdateError when the checksum is missing or does not match (the file is deleted), when an
    address is not allowed, on a network error or when `cancel` (a threading.Event) is set. Meant to run in a worker thread."""
    get = _pick(fetch)
    check_url(update.asset_url)
    if not update.sha256_url:
        raise UpdateError("This release has no checksum file, so the download cannot be verified: not installed.")
    check_url(update.sha256_url)
    expected = parse_checksum(get(update.sha256_url, timeout=TIMEOUT), ASSET_NAME)
    if expected is None:
        raise UpdateError(f"The checksum file ({CHECKSUM_NAME}) is not in the expected form: not installed.")
    part, final = Path(folder) / (ASSET_NAME + ".part"), Path(folder) / (ASSET_NAME + ".new")
    try:
        try:
            with open(part, "wb") as out:
                get(update.asset_url, timeout=TIMEOUT, sink=out, progress=progress, cancel=cancel)
            actual = sha256_of(part)
            if not compare_digest(actual.encode(), expected.encode()):
                raise UpdateError("The downloaded file does not match its published checksum (corrupted or altered): deleted, not "
                                  "installed. Try again later.")
            os.replace(part, final)
        except OSError as e:
            raise UpdateError(f"Could not write the update in {folder}: {e}") from None
    except BaseException:
        for leftover in (part, final):
            try:
                leftover.unlink(missing_ok=True)
            except OSError:
                pass
        raise
    return final


# ---- installing ----------------------------------------------------------------------------------------------------------------------

def is_frozen():
    return bool(getattr(sys, "frozen", False))


def is_folder_build():
    """True for the one-folder build (openDNA-folder.zip), where the program sits next to its files: PyInstaller's bundle folder is
    then the exe's folder (version 5) or its `_internal` subfolder (6); a one-file exe unpacks into a temp folder instead."""
    meipass = getattr(sys, "_MEIPASS", None)
    if not (is_frozen() and meipass):
        return False
    here = Path(sys.executable).resolve().parent
    return Path(meipass).resolve() in (here, here / "_internal")


def install_blocker(current_exe=None):
    """Why this run cannot replace itself (a sentence), or None. Passing `current_exe` means "act on that file" (tests, and
    nothing else: the program itself passes nothing and then only a PyInstaller one-file exe qualifies)."""
    if current_exe is None and not is_frozen():
        return ("openDNA is running from its source code, which cannot replace itself. Update with git, or download "
                f"{ASSET_NAME} from {project.RELEASES_URL}.")
    if current_exe is None and sys.platform != "win32":
        # The release's one-file program is the Windows exe: a Mac or a Linux build must never swap that in for itself.
        import platform
        mac = f"{project.EXE_NAME}-macos-{'arm64' if platform.machine() == 'arm64' else 'x86_64'}.zip"
        what = f"the new openDNA for macOS ({mac})" if sys.platform == "darwin" else \
            f"the new openDNA for Linux ({project.EXE_NAME}-linux-x86_64.tar.gz)"
        return f"This version of openDNA does not update itself. Download {what} from {project.RELEASES_URL}."
    if current_exe is None and is_folder_build():
        # Swapping one exe would leave the Qt files the user may have put in this folder untouched but mismatched with the new version.
        return ("This is the folder version of openDNA (the one whose Qt files you can replace), so it does not update itself. "
                f"Download the new openDNA-folder.zip from {project.RELEASES_URL}.")
    return None


def _retry(action, what):
    """Run action(); on PermissionError / OSError try again a few times (antivirus), then raise UpdateError."""
    for attempt in range(RETRIES):
        try:
            return action()
        except OSError as e:
            if attempt == RETRIES - 1:
                folder_hint = ""
                if isinstance(e, PermissionError):
                    folder_hint = (" The folder may need administrator rights (for example Program Files) or an antivirus is holding "
                                   f"the file: download the new version by hand from {project.RELEASES_URL} instead.")
                raise UpdateError(f"Could not {what}: {e}.{folder_hint}") from None
            time.sleep(RETRY_DELAY)


def relaunch(exe):
    """Start `exe` as a separate program that outlives us. A PyInstaller --onefile exe that starts itself hands its unpack folder to the
    child through environment variables; that folder disappears when we quit, so they are removed (PYINSTALLER_RESET_ENVIRONMENT is
    the official switch of PyInstaller 6.9 and later)."""
    env = {k: v for k, v in os.environ.items() if not k.startswith("_PYI") and k != "_MEIPASS2"}
    env["PYINSTALLER_RESET_ENVIRONMENT"] = "1"
    flags = getattr(subprocess, "DETACHED_PROCESS", 0x8) | getattr(subprocess, "CREATE_NEW_PROCESS_GROUP", 0x200) if os.name == "nt" else 0
    subprocess.Popen([str(exe)], env=env, cwd=str(Path(exe).parent), creationflags=flags, close_fds=True)


def install_update(new_exe, current_exe=None, relaunch=None):
    """Swap the verified `new_exe` in for the running exe and start it. The caller quits right after this returns.
    UpdateError (with everything put back as it was) when it cannot. `relaunch(path)` is injectable for tests."""
    blocker = install_blocker(current_exe)
    if blocker:
        raise UpdateError(blocker)
    start = relaunch if relaunch is not None else globals()["relaunch"]
    current = Path(current_exe or sys.executable)
    new_exe = Path(new_exe)
    old = current.with_name(current.name + ".old")
    if not new_exe.is_file() or not current.is_file():
        raise UpdateError("The new version or the running program could not be found on disk.")
    if new_exe.parent != current.parent:
        raise UpdateError("The new version must be in the same folder as the program (it is moved by renaming).")
    if old.exists():     # left by an earlier update that could not clean up: the new .old must be the version we are leaving
        _retry(lambda: old.unlink(), f"remove the old file {old.name}")
    # Windows cannot overwrite a running exe but can rename it, even while it runs.
    _retry(lambda: os.replace(current, old), f"rename {current.name}")
    try:
        _retry(lambda: os.replace(new_exe, current), f"put the new {current.name} in place")
    except UpdateError:
        _retry(lambda: os.replace(old, current), f"put {current.name} back")
        raise
    try:
        start(current)
    except Exception as e:
        # The new exe will not start: go back to the version that is running, and keep the new file for another try.
        _retry(lambda: os.replace(current, new_exe), "undo the update")
        _retry(lambda: os.replace(old, current), f"put {current.name} back")
        raise UpdateError(f"The new version could not be started ({e}); the old one was kept.") from None


def cleanup_old(current_exe=None, wait=0.0):
    """Delete what an update leaves behind (openDNA.exe.old and any .part / .new). Call at start-up. The old program may still be
    quitting when the new one starts, so with `wait` seconds it keeps trying. True when nothing is left. Does nothing from source."""
    if current_exe is None and not is_frozen():
        return True
    current = Path(current_exe or sys.executable)
    leftovers = [current.with_name(current.name + s) for s in (".old", ".part", ".new")]
    deadline = time.monotonic() + wait
    while True:
        for f in leftovers:
            try:
                f.unlink(missing_ok=True)
            except OSError:
                pass
        left = any(f.exists() for f in leftovers)
        if not left or time.monotonic() >= deadline:
            return not left
        time.sleep(0.5)


def cleanup_old_in_background():
    """cleanup_old() for the start-up, without delaying the window (it may wait for the old process to exit)."""
    threading.Thread(target=cleanup_old, kwargs={"wait": 15.0}, daemon=True).start()


# ---- what is remembered --------------------------------------------------------------------------------------------------------------

def _state():
    return QSettings(str(data_dir() / "updates.ini"), QSettings.IniFormat)


def check_enabled():
    """The setting "look for a new version when openDNA starts" (Tools > Settings > Updates), on by default."""
    return bool(tool_settings.value(TOOL, "check"))


def skipped_version():
    return str(_state().value("skipped", "") or "")


def skip_version(version):
    s = _state()
    s.setValue("skipped", str(version))
    s.sync()


def mark_checked(now=None):
    s = _state()
    s.setValue("last_check", float(time.time() if now is None else now))
    s.sync()


def auto_check_due(now=None):
    """True when the start-up check should run: updates available, the setting on, and none in the last 24 hours."""
    if not project.UPDATES_ENABLED or not check_enabled():
        return False
    now = time.time() if now is None else now
    try:
        last = float(_state().value("last_check", 0) or 0)
    except (TypeError, ValueError):
        last = 0.0
    return last > now or now - last >= CHECK_EVERY      # last in the future: the clock was changed, check anyway


# ---- threads and the window ----------------------------------------------------------------------------------------------------------

_alive = set()   # running jobs and open message boxes must not be garbage collected

# Python's cycle collector is paused while a worker thread runs. A full-suite run on CI crashed (access violation) in a worker thread that
# was opening a file while the GUI thread was inside a garbage collection; PySide gives up the GIL while it deletes Qt objects, so
# collecting dead widgets can let a worker thread run in the middle of the collection. Reference counting still frees almost everything
# meanwhile; the collector runs again once the last job is done. (Not reproduced on demand: this removes the overlap, it does not
# prove the cause.)
_gc_lock = threading.Lock()
_gc_jobs = 0
_gc_was_on = False


def _gc_pause():
    global _gc_jobs, _gc_was_on
    with _gc_lock:
        if _gc_jobs == 0:
            _gc_was_on = gc.isenabled()
            gc.disable()
        _gc_jobs += 1


def _gc_resume():
    global _gc_jobs
    with _gc_lock:
        _gc_jobs -= 1
        if _gc_jobs == 0 and _gc_was_on:
            gc.enable()


class _Job(QObject):
    """Runs work(progress) in a thread and calls callback(result) in the GUI thread (an exception from work is the result)."""
    done = Signal(object)
    progress = Signal(object, object)

    def __init__(self, work, callback, on_progress=None):
        super().__init__()
        self.callback = callback
        self.done.connect(self._finish)          # a method of an object living in the GUI thread: delivered there, not in the worker
        if on_progress is not None:
            self.progress.connect(on_progress)   # keep to methods of QObjects here, same reason
        _alive.add(self)
        _gc_pause()
        threading.Thread(target=self._run, args=(work,), daemon=True).start()

    def _run(self, work):
        try:
            result = work(self.progress.emit)
        except Exception as e:
            result = e
        self.done.emit(result)

    def _finish(self, result):
        _alive.discard(self)
        _gc_resume()
        try:
            self.callback(result)
        except RuntimeError:     # the window the answer was for has been deleted meanwhile
            pass


def check_in_background(callback, current_version=None, fetch=None):
    """check() in a thread; callback(Check) runs in the GUI thread. Returns the job (keeps itself alive)."""
    return _Job(lambda progress: check(current_version, fetch), callback)


def show_message(parent, title, text):
    """A non-blocking message box (rich text allowed). Tests replace this function."""
    box = QMessageBox(QMessageBox.Information, title, text, QMessageBox.Ok, parent)
    box.setTextFormat(Qt.RichText)
    box.setTextInteractionFlags(Qt.TextBrowserInteraction)
    _alive.add(box)
    box.finished.connect(lambda *_: _alive.discard(box))
    box.open()
    return box


def _link(url, label=None):
    return f'<a href="{html.escape(url)}">{html.escape(label or url)}</a>'


def open_dialog(window, update, **kwargs):
    """Show the UpdateDialog for `update` without blocking (window-modal); kept on window._update_dialog."""
    try:
        dialog = UpdateDialog(update, window, **kwargs)
    except RuntimeError:     # the window is gone
        return None
    window._update_dialog = dialog
    dialog.open()
    return dialog


def startup_check(window, fetch=None, now=None, **dialog_kwargs):
    """Call once when the main window is shown. Silent: does nothing unless an automatic check is due (updates available, setting on,
    not in the last 24 hours); then looks in a thread and opens the dialog only for a new version not skipped. True if it looked."""
    cleanup_old_in_background()
    if not auto_check_due(now):
        return False
    mark_checked(now)

    def got(result):
        if result.update is not None and result.update.version != skipped_version():
            open_dialog(window, result.update, **dialog_kwargs)
    check_in_background(got, fetch=fetch)
    return True


def check_now(window, fetch=None, **dialog_kwargs):
    """The menu action "Check for Updates...". Says "not available yet" while project.UPDATES_ENABLED is False (no network access then);
    else looks, and reports "up to date", an error, or opens the dialog (even for a version that was skipped). True if it looked."""
    if not project.UPDATES_ENABLED:
        show_message(window, "Check for Updates",
                     f"Updates are not available yet.<br>New versions will be published at {_link(project.RELEASES_URL)}.")
        return False
    mark_checked()

    def got(result):
        if result.update is not None:
            open_dialog(window, result.update, **dialog_kwargs)
        elif result.error:
            page = result.page_url or project.RELEASES_URL
            show_message(window, "Check for Updates", f"{html.escape(result.error)}<br><br>Releases: {_link(page)}")
        else:
            show_message(window, "Check for Updates", f"You are up to date: {project.NAME} {project.VERSION} is the latest version.")
    check_in_background(got, fetch=fetch)
    return True


# ---- the dialog ----------------------------------------------------------------------------------------------------------------------

class UpdateDialog(QDialog):
    """Release notes, Install and restart (download with a progress bar, verify, swap, relaunch), Later, Skip this version.

    Use open() or show(), never exec(). The optional arguments are for tests and for the main window:
      confirm()  called after the download, before the swap: return False (after asking about unsaved work) to stop;
      quit()     called after a successful install: the main window should close; default QApplication.quit;
      folder, fetch, install, relaunch  where to download, and replacements of the network / install steps."""

    def __init__(self, update, parent=None, folder=None, fetch=None, install=None, relaunch=None, quit=None, confirm=None,
                 current_version=None):
        super().__init__(parent)
        self.update, self.cancel = update, threading.Event()
        self.folder, self._fetch = Path(folder) if folder else Path(sys.executable).parent, fetch
        self._install = install or (lambda path: install_update(path, relaunch=relaunch))
        self._quit, self._confirm = quit or QApplication.quit, confirm
        self.blocker = None if install else install_blocker()
        self.busy = self.installed = False
        self.setWindowTitle("Update Available")
        self.setMinimumWidth(480)

        current = current_version or project.VERSION
        head = QLabel(f"<b>{html.escape(project.NAME)} {html.escape(update.version)}</b> is available (you have {html.escape(current)}).")
        self.notes = QPlainTextEdit(update.notes or "No release notes.", readOnly=True)   # plain text: a release description is not trusted
        self.notes.setMinimumHeight(140)
        self.link = QLabel(_link(update.page_url, "Open the release page"), openExternalLinks=True)
        self.status = QLabel(wordWrap=True, openExternalLinks=True, textFormat=Qt.RichText)
        self.bar = QProgressBar()
        policy = self.bar.sizePolicy()
        policy.setRetainSizeWhenHidden(True)   # its room is kept while hidden: the dialog does not grow when a download starts
        self.bar.setSizePolicy(policy)
        self.bar.hide()
        self.install_button = QPushButton("Install and Restart", default=True, clicked=self._start)
        self.later = QPushButton("Later", clicked=self.reject)
        self.skip = QPushButton("Skip This Version", clicked=self._skip)
        # each button as wide as the longest text it takes ("Later" / "Cancel", "Install and Restart" / "Try Again"): the
        # buttons do not move when a download starts or fails
        for button, texts in ((self.later, ("Later", "Cancel")), (self.install_button, ("Install and Restart", "Try Again"))):
            button.setMinimumWidth(max(button.fontMetrics().horizontalAdvance(t) for t in texts) + 32)
        if self.blocker:
            self.install_button.setEnabled(False)
            self.status.setText(html.escape(self.blocker).replace(project.RELEASES_URL, _link(project.RELEASES_URL)))
        buttons = QHBoxLayout()
        buttons.addWidget(self.skip)
        buttons.addStretch()
        buttons.addWidget(self.later)
        buttons.addWidget(self.install_button)
        layout = QVBoxLayout(self)
        for w in (head, self.notes, self.link, self.status, self.bar):
            layout.addWidget(w)
        layout.addLayout(buttons)

    def _skip(self):
        skip_version(self.update.version)
        self.reject()

    def reject(self):
        if self.busy:
            self.cancel.set()     # a download in progress stops and deletes its partial file
        super().reject()

    def _set_busy(self, busy):
        self.busy = busy
        self.install_button.setEnabled(not busy and not self.blocker)
        self.skip.setEnabled(not busy)
        self.later.setText("Cancel" if busy else "Later")

    def _start(self):
        self.cancel.clear()
        self._set_busy(True)
        self.bar.setRange(0, 0)
        self.bar.show()
        self.status.setText("Downloading…")
        _Job(lambda progress: download_update(self.update, self.folder, progress, self._fetch, self.cancel), self._downloaded,
             self._on_progress)

    def _on_progress(self, done, total):
        if total:
            self.bar.setRange(0, total)
            self.bar.setValue(min(done, total))
        self.status.setText(f"Downloading… {done / 1e6:.1f} of {total / 1e6:.1f} MB" if total else f"Downloading… {done / 1e6:.1f} MB")

    def _failed(self, message):
        self._set_busy(False)
        self.bar.hide()
        self.install_button.setText("Try Again")
        self.status.setText(html.escape(message))

    def _downloaded(self, result):
        if self.cancel.is_set():
            self._set_busy(False)
            self.bar.hide()
            return
        if isinstance(result, Exception):
            self._failed(str(result) if isinstance(result, UpdateError) else f"The download failed: {result}")
            return
        if self._confirm is not None and not self._confirm():
            Path(result).unlink(missing_ok=True)
            self._failed("The update was not installed.")
            return
        self.status.setText("Installing…")
        try:
            self._install(result)
        except UpdateError as e:
            self._failed(str(e))
            return
        self.installed = True
        self.status.setText("Restarting…")
        self.accept()
        self._quit()
