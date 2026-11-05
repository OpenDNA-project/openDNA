"""Testing feedback: a log file, a friendly "report it?" dialog for unexpected errors, and prefilled GitHub issue links.

What this module does (Help menu items map to the functions in brackets):
  * setup_logging()      a rotating log file <data folder>/logs/genomeedit.log (3 files x 500 KB), a start-up line with the
                         versions, Qt warnings copied into it, and hooks so an error nobody caught is logged and shown once in
                         a small dialog (Report on GitHub / Copy details / Close) instead of killing the program.
  * open_report("bug")   Help > Report a Problem...   opens GitHub's new-issue form, already filled with versions and log lines.
  * open_feature()       Help > Request a Feature...
  * open_log_folder()    Help > Open Log Folder
  * copy_diagnostics()   Help > Troubleshooting > Copy Diagnostics   (versions + the last log lines, ready to paste anywhere)

Privacy rule: nothing that came from the user's files leaves the machine unless the user presses a button, and even then every
text passes through scrub() first: DNA/protein sequences, file names, file paths and the Windows user name are replaced by
placeholders. The log file itself is written through the same scrubbing, so attaching it to an issue is safe.

Tests (tests/test_feedback.py) never open a browser or a dialog: the opener and the dialog are injectable.
"""
import getpass
import logging
import logging.handlers
import os
import platform
import re
import sys
import threading
import traceback
from pathlib import Path
from urllib.parse import quote, urlencode

from PySide6.QtCore import QObject, QUrl, Signal, qInstallMessageHandler, qVersion
from PySide6.QtCore import QtMsgType

import project
from primer_db import data_dir

LOG_NAME = "genomeedit.log"
LOG_MAX_BYTES = 500 * 1024      # one file; with 2 older copies kept that is 3 x 500 KB at most
LOG_BACKUPS = 2
LOG_TAIL_LINES = 30             # lines of the log put into a report
URL_CAP = 7000                  # GitHub refuses very long addresses (about 8 KB); stay well below
SUMMARY_MAX = 120

log = logging.getLogger("genomeedit")


# --------------------------------------------------------------------------------------------------- scrubbing

# Extensions of files a user may have opened or exported: their names can reveal what the user works on.
_FILE_EXTENSIONS = ("gb", "gbk", "gbff", "genbank", "dna", "fa", "fasta", "fna", "fas", "fastq", "fq", "gff", "gff3", "bed",
                    "bam", "bai", "sam", "vcf", "tsv", "csv", "txt", "ab1", "gz", "json")
_RE_FRAME = re.compile(r'File "([^"\n]*)"')
_RE_QUOTED_PATH = re.compile(r"""(['"])((?:[A-Za-z]:[\\/]|\\\\|/)[^'"\n]*[\\/][^'"\n]*)\1""")
_RE_WIN_PATH = re.compile(r"(?<![A-Za-z0-9])[A-Za-z]:[\\/][^\s'\"<>|*?]*")
_RE_UNC_PATH = re.compile(r"\\\\[^\s'\"<>|*?\\]+\\[^\s'\"<>|*?]*")
_RE_POSIX_PATH = re.compile(r"(?<![\w:/.~])/(?:[^\s/'\"<>]+/)+[^\s/'\"<>]*")
_RE_FILE_NAME = re.compile(r"(?<![\w.-])[\w\-.]+\.(?:%s)(?![\w])" % "|".join(_FILE_EXTENSIONS), re.IGNORECASE)
# 12+ nucleotide letters (IUPAC) in a row, or 25+ capitals in a row (a protein), or GenBank-style blocks of ten letters.
_RE_DNA = re.compile(r"(?<![A-Za-z])[ACGTUNRYKMSWBDHVacgtunrykmswbdhv]{12,}(?![A-Za-z])")
_RE_PROTEIN = re.compile(r"(?<![A-Za-z])[A-Z]{25,}(?![A-Za-z])")
_RE_BLOCKS = re.compile(r"(?<![A-Za-z])(?:[ACGTNacgtn]{10}[ \t]+){2,}[ACGTNacgtn]{1,10}(?![A-Za-z])")


def _user_names():
    """Names that identify the person on this machine (account name, home folder name)."""
    names = set()
    for getter in (getpass.getuser, lambda: os.environ.get("USERNAME"), lambda: os.environ.get("USER"),
                   lambda: Path.home().name):
        try:
            value = getter()
        except Exception:  # getpass raises when no account name can be found
            continue
        if value and len(value) >= 3:
            names.add(value)
    return names


def scrub(text):
    """Return `text` with everything private replaced: stack-frame paths keep only the file name, other absolute paths become
    <path>, file names with a sequence-file extension become <file>, long runs of bases or residues become <sequence>, and the
    user's account name becomes <user>. Safe to call on any string (and idempotent)."""
    if not text:
        return ""
    text = str(text)
    text = _RE_FRAME.sub(lambda m: 'File "%s"' % re.split(r"[\\/]", m.group(1))[-1], text)
    text = _RE_QUOTED_PATH.sub(lambda m: "%s<path>%s" % (m.group(1), m.group(1)), text)
    text = _RE_UNC_PATH.sub("<path>", text)
    text = _RE_WIN_PATH.sub("<path>", text)
    text = _RE_POSIX_PATH.sub("<path>", text)
    text = _RE_FILE_NAME.sub("<file>", text)
    text = _RE_BLOCKS.sub("<sequence>", text)
    text = _RE_DNA.sub("<sequence>", text)
    text = _RE_PROTEIN.sub("<sequence>", text)
    for name in sorted(_user_names(), key=len, reverse=True):
        text = re.sub(r"(?<![A-Za-z0-9])%s(?![A-Za-z0-9])" % re.escape(name), "<user>", text, flags=re.IGNORECASE)
    return text


def format_exception(exc_type, exc, tb):
    """The scrubbed traceback of an exception, as one string."""
    return scrub("".join(traceback.format_exception(exc_type, exc, tb)))


# --------------------------------------------------------------------------------------------------- system info and log

def log_folder():
    """The folder holding the log file (created when missing)."""
    path = data_dir() / "logs"
    path.mkdir(parents=True, exist_ok=True)
    return path


def log_path():
    return log_folder() / LOG_NAME


def _qt_versions():
    try:
        import PySide6
        pyside = PySide6.__version__
    except Exception:
        pyside = "?"
    return pyside, qVersion()


def system_info():
    """Lines describing this installation (no account name, no folder names)."""
    pyside, qt = _qt_versions()
    return [
        f"{project.NAME} {project.VERSION}",
        f"OS: {platform.platform()}",
        f"Python: {platform.python_version()} ({platform.machine()})",
        f"PySide6 {pyside} / Qt {qt}",
    ]


def last_log_lines(count=LOG_TAIL_LINES):
    """The last `count` lines of the log file, scrubbed (the file is already written scrubbed; this is a second net)."""
    try:
        text = log_path().read_text(encoding="utf-8", errors="replace")
    except OSError:
        return []
    return scrub(text).splitlines()[-count:]


def diagnostics_text(details=None):
    """Everything useful for a bug report, as plain text: versions, the traceback if any, the last log lines."""
    parts = system_info()
    if details:
        parts += ["", "--- Error ---", scrub(details).rstrip()]
    lines = last_log_lines()
    parts += ["", f"--- Last {LOG_TAIL_LINES} log lines ---"] + (lines or ["(the log is empty)"])
    return "\n".join(parts) + "\n"


def copy_text(text, clipboard=None):
    """Put `text` on the clipboard (`clipboard` is injectable for tests: any object with setText)."""
    if clipboard is None:
        from PySide6.QtGui import QGuiApplication
        clipboard = QGuiApplication.clipboard()
    if clipboard is not None:
        clipboard.setText(text)


def copy_diagnostics(clipboard=None, details=None):
    """Help > Copy Diagnostics. Returns the text copied."""
    text = diagnostics_text(details)
    copy_text(text, clipboard)
    return text


# --------------------------------------------------------------------------------------------------- URLs

def _clean_summary(summary):
    return re.sub(r"\s+", " ", scrub(summary)).strip()[:SUMMARY_MAX]


def _url(template, params):
    return f"{project.ISSUES_URL}/new?" + urlencode([("template", template)] + params, quote_via=quote, safe="")


def report_url(kind, summary="", details=None):
    """A GitHub new-issue link for the bug form, prefilled with the version, OS, Python/Qt versions and either `details` (a
    traceback) or the last log lines. `kind` is "bug" or "crash" (a crash is a bug with an automatic summary); "feature" gives
    feature_url(). The address is at most URL_CAP characters: the oldest lines are cut first. Everything is scrubbed."""
    if kind == "feature":
        return feature_url(summary)
    if kind not in ("bug", "crash"):
        raise ValueError(f"unknown report kind: {kind!r}")
    pyside, qt = _qt_versions()
    title = "[Bug]: " + _clean_summary(summary)
    body = scrub(details).rstrip() if details else "\n".join(last_log_lines())
    fixed = [("title", title), ("version", f"{project.VERSION}"),
             ("os", f"{platform.platform()}; Python {platform.python_version()}; PySide6 {pyside}; Qt {qt}")]

    def build(n):
        logs = body
        if n < len(body):
            logs = "[... earlier lines left out ...]\n" + body[len(body) - n:] if n > 0 else ""
        return _url("bug_report.yml", fixed + ([("logs", logs)] if logs else []))

    url = build(len(body))
    if len(url) > URL_CAP:
        low, high = 0, len(body)           # largest tail of the text that still fits
        while low < high:
            mid = (low + high + 1) // 2
            if len(build(mid)) <= URL_CAP:
                low = mid
            else:
                high = mid - 1
        url = build(low)
    return url


def feature_url(summary=""):
    """A GitHub new-issue link for the feature-request form."""
    return _url("feature_request.yml", [("title", "[Feature]: " + _clean_summary(summary))])


def _default_opener(url):
    from PySide6.QtGui import QDesktopServices
    return QDesktopServices.openUrl(QUrl(url))


def open_report(kind="bug", summary="", details=None, opener=None):
    """Help > Report a Problem... Opens the prefilled form in the browser and returns the address."""
    url = report_url(kind, summary, details)
    (opener or _default_opener)(url)
    return url


def open_feature(summary="", opener=None):
    """Help > Request a Feature..."""
    url = feature_url(summary)
    (opener or _default_opener)(url)
    return url


def open_log_folder(opener=None):
    """Help > Open Log Folder. Returns the folder."""
    folder = log_folder()
    (opener or _default_opener)(QUrl.fromLocalFile(str(folder)).toString())
    return folder


# --------------------------------------------------------------------------------------------------- logging and hooks

class _ScrubFormatter(logging.Formatter):
    """Every line written to the log file passes through scrub()."""

    def format(self, record):
        return scrub(super().format(record))


class _Bridge(QObject):
    """Lets a thread other than the main one ask for the dialog: the signal is delivered in the thread that made the bridge."""
    requested = Signal(str, str)


_state = {"handler": None, "hooks": False, "old_excepthook": None, "old_threading_hook": None, "qt_handler": None,
          "qt_installed": False, "old_qt_handler": None, "old_level": None, "dialog": None, "opener": None, "clipboard": None,
          "seen": set(), "showing": False, "bridge": None, "path": None}

_QT_LEVELS = {QtMsgType.QtDebugMsg: logging.DEBUG, QtMsgType.QtInfoMsg: logging.INFO, QtMsgType.QtWarningMsg: logging.WARNING,
              QtMsgType.QtCriticalMsg: logging.ERROR, QtMsgType.QtFatalMsg: logging.CRITICAL}


def _qt_message(mode, context, message):
    try:
        logging.getLogger("qt").log(_QT_LEVELS.get(mode, logging.WARNING), "%s", message)
        if sys.stderr is not None:
            print(scrub(message), file=sys.stderr)   # Qt's own handler printed it; keep doing so for developers
    except Exception:
        pass  # a logging problem must never come back into Qt


def _default_dialog(summary, details):
    """Ask the user what to do with an error. Returns "report", "copy" or "close" (None when no window can be shown)."""
    from PySide6.QtWidgets import QApplication, QMessageBox
    if QApplication.instance() is None:
        return None
    box = QMessageBox(QMessageBox.Warning, "Something went wrong", "Something went wrong. Report it?")
    box.setInformativeText(f"{summary}\n\nThe program keeps running; your work is not lost. The report is prefilled and contains "
                           "no sequences or file names; you see it before sending anything.")
    report = box.addButton("Report on GitHub", QMessageBox.AcceptRole)
    copy = box.addButton("Copy details", QMessageBox.ActionRole)
    close = box.addButton("Close", QMessageBox.RejectRole)
    box.setDefaultButton(close)
    box.exec()
    clicked = box.clickedButton()
    return "report" if clicked is report else "copy" if clicked is copy else "close"


def _present(summary, details):
    """Show the dialog (main thread) and act on the answer. Never raises."""
    if _state["showing"]:
        return
    _state["showing"] = True
    try:
        choice = (_state["dialog"] or _default_dialog)(summary, details)
        if choice == "report":
            open_report("crash", summary, details, opener=_state["opener"])
        elif choice == "copy":
            copy_diagnostics(clipboard=_state["clipboard"], details=details)
    except Exception:
        log.exception("The error dialog itself failed")
    finally:
        _state["showing"] = False


def _handle_exception(exc_type, exc, tb, where):
    if issubclass(exc_type, KeyboardInterrupt):
        sys.__excepthook__(exc_type, exc, tb)
        return
    try:
        log.error("Unhandled exception (%s)", where, exc_info=(exc_type, exc, tb))
        frames = traceback.extract_tb(tb)
        last = frames[-1] if frames else None
        key = (exc_type.__qualname__, re.split(r"[\\/]", last.filename)[-1] if last else "", last.lineno if last else 0)
        if key in _state["seen"]:
            return                      # the same error again: logged, but the user is not asked twice
        _state["seen"].add(key)
        details = format_exception(exc_type, exc, tb)
        summary = scrub(f"{exc_type.__name__}: {exc}").splitlines()[0][:200]
        if _state["dialog"] is not None or threading.current_thread() is threading.main_thread():
            _present(summary, details)
        elif _state["bridge"] is not None:
            _state["bridge"].requested.emit(summary, details)
        if sys.stderr is not None:
            sys.__excepthook__(exc_type, exc, tb)   # the usual traceback for a developer running from a terminal
    except Exception:
        pass  # whatever goes wrong here, the program must go on


def _excepthook(exc_type, exc, tb):
    _handle_exception(exc_type, exc, tb, "main thread")


def _threading_hook(args):
    if args.exc_type is SystemExit:
        return
    _handle_exception(args.exc_type, args.exc_value, args.exc_traceback, f"thread {getattr(args.thread, 'name', '?')}")


def setup_logging(hooks=True, qt_handler=True, show_dialog=None, opener=None, clipboard=None, level=logging.INFO):
    """Start logging to <data folder>/logs/genomeedit.log and (optionally) install the error hooks. Call once, right after the
    QApplication is created. `show_dialog(summary, details)`, `opener(url)` and `clipboard` replace the real dialog, browser and clipboard (tests).
    Calling it again first undoes the previous call. Returns the log file path."""
    teardown_logging()
    path = log_path()
    handler = logging.handlers.RotatingFileHandler(path, maxBytes=LOG_MAX_BYTES, backupCount=LOG_BACKUPS, encoding="utf-8")
    handler.setFormatter(_ScrubFormatter("%(asctime)s %(levelname)-7s %(name)s: %(message)s"))
    root = logging.getLogger()
    _state["old_level"] = root.level
    root.addHandler(handler)
    root.setLevel(level)
    _state.update(handler=handler, path=path, dialog=show_dialog, opener=opener, clipboard=clipboard, seen=set(), showing=False)
    pyside, qt = _qt_versions()
    log.info("Start: %s %s | Python %s | PySide6 %s / Qt %s | %s", project.NAME, project.VERSION, platform.python_version(),
             pyside, qt, platform.platform())
    if hooks:
        _state["bridge"] = _Bridge()
        _state["bridge"].requested.connect(_present)
        _state["old_excepthook"], sys.excepthook = sys.excepthook, _excepthook
        _state["old_threading_hook"], threading.excepthook = threading.excepthook, _threading_hook
        _state["hooks"] = True
    if qt_handler:
        _state["old_qt_handler"] = qInstallMessageHandler(_qt_message)
        _state["qt_installed"] = True
    return path


def teardown_logging():
    """Undo setup_logging(): remove the file handler (closing the file), the hooks and the Qt message handler."""
    if _state["qt_installed"]:
        qInstallMessageHandler(_state["old_qt_handler"])
        _state["qt_installed"] = False
    if _state["hooks"]:
        if sys.excepthook is _excepthook:
            sys.excepthook = _state["old_excepthook"]
        if threading.excepthook is _threading_hook:
            threading.excepthook = _state["old_threading_hook"]
        _state["hooks"] = False
        _state["bridge"] = None
    handler = _state["handler"]
    if handler is not None:
        root = logging.getLogger()
        root.removeHandler(handler)
        handler.close()
        if _state["old_level"] is not None:
            root.setLevel(_state["old_level"])
        _state["handler"] = None
    _state.update(dialog=None, opener=None, clipboard=None)
