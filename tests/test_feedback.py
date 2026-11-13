"""Feedback (app/feedback.py): scrubbing of private text, prefilled GitHub links (encoding, length cap, nothing private), the
rotating log file and the error hooks (one dialog per distinct error). Offline: no browser, no dialog, no clipboard is ever touched; the
opener, the dialog and the clipboard are replaced by fakes. Each test installs and removes the hooks itself (pytest swaps
threading.excepthook around every phase, so a fixture could not keep it)."""
import contextlib
import logging
import os
import sys
import tempfile
import threading
from pathlib import Path
from urllib.parse import parse_qs, urlsplit

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ["OPENDNA_DATA_DIR"] = tempfile.mkdtemp(prefix="opendna_test_")
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "app"))

import feedback
import project

DNA = "ATGGCTAGCTAGGATCCGATCGATCGTAGCTAGCTAGC"
SECRET_PATH = r"C:\Users\bobsmith\Documents\Lab Notes\my_secret_plasmid.gb"


@pytest.fixture(autouse=True)
def private_data_dir(tmp_path, monkeypatch):
    monkeypatch.setenv("OPENDNA_DATA_DIR", str(tmp_path))
    monkeypatch.setenv("USERNAME", "bobsmith")
    monkeypatch.setenv("USER", "bobsmith")
    monkeypatch.setattr("getpass.getuser", lambda: "bobsmith")
    monkeypatch.setattr(Path, "home", classmethod(lambda cls: Path(r"C:\Users\bobsmith")))
    return tmp_path


@contextlib.contextmanager
def active(**kwargs):
    """setup_logging with fakes; everything is put back afterwards (hooks, handlers, root level)."""
    root = logging.getLogger()
    before = (sys.excepthook, threading.excepthook, list(root.handlers), root.level)
    calls = {"dialog": [], "urls": [], "copied": []}

    class Clip:
        def setText(self, text):
            calls["copied"].append(text)

    choice = kwargs.pop("choice", "close")

    def dialog(summary, details):
        calls["dialog"].append((summary, details))
        if choice == "boom":
            raise RuntimeError("dialog broke")
        return choice

    feedback.setup_logging(show_dialog=dialog, opener=calls["urls"].append, clipboard=Clip(), **kwargs)
    try:
        yield calls
    finally:
        feedback.teardown_logging()
        after = (sys.excepthook, threading.excepthook, list(root.handlers), root.level)
        assert after == before, "setup_logging/teardown_logging did not restore hooks, handlers or level"


def boom(message="x"):
    raise ValueError(message)


def boom_other():
    raise KeyError("k")


def raised(fn, *args):
    try:
        fn(*args)
    except BaseException:
        return sys.exc_info()


def query(url):
    return parse_qs(urlsplit(url).query, keep_blank_values=True)


def test_scrub_removes_private_text():
    """scrub removes the username and home path, quoted and UNC paths, traceback folders (keeping the file name), file names and
    sequences; ordinary text and URLs are left alone; it is idempotent."""
    # --- username and home path
    out = feedback.scrub(r'Cannot write C:\Users\bobsmith\AppData\Roaming\openDNA\x.tsv for bobsmith')
    assert "bobsmith" not in out and "Users" not in out and "AppData" not in out
    assert "/home/bobsmith/data" not in feedback.scrub("failed: /home/bobsmith/data/run1")
    assert "bobsmith" not in feedback.scrub("failed: /home/bobsmith/data/run1")
    assert "bobsmith" not in feedback.scrub("owner BobSmith")          # case does not matter

    # --- paths quoted with spaces and UNC
    out = feedback.scrub(f"FileNotFoundError: [Errno 2] No such file: {SECRET_PATH!r}")
    assert "Lab Notes" not in out and "secret" not in out and "<path>" in out
    assert "Lab" not in feedback.scrub(r"\\fileserver\lab\project\plasmid.fa not found")

    # --- traceback frames keep only the file name
    out = feedback.scrub('  File "C:\\Users\\bobsmith\\miniconda3\\lib\\site-packages\\x\\y.py", line 12, in run')
    assert 'File "y.py", line 12, in run' in out and "bobsmith" not in out

    # --- file names and sequences
    out = feedback.scrub(f"Cannot open my_secret_plasmid.gb and reads.FASTA; sequence {DNA}; protein " + "ACDEFGHIKLMNPQRSTVWY" * 2)
    assert "secret" not in out and "reads" not in out and DNA not in out and "ACDEFGHIK" not in out
    assert out.count("<sequence>") == 2 and out.count("<file>") == 2
    assert "acgtacgtac gtacgtacgt acgtacgtac" not in feedback.scrub("1 acgtacgtac gtacgtacgt acgtacgtac gt")
    assert "acgt" not in feedback.scrub("1 acgtacgtac gtacgtacgt acgtacgtac gt")

    # --- ordinary text and URLs are left alone
    text = "Unhandled exception in main thread: initialization configuration failed (see https://github.com/OpenDNA-project/openDNA/issues)"
    assert feedback.scrub(text) == text
    assert feedback.scrub("") == "" and feedback.scrub(None) == ""

    # --- idempotent
    once = feedback.scrub(f"{SECRET_PATH} {DNA} bobsmith")
    assert feedback.scrub(once) == once


def test_report_links_are_encoded_capped_and_carry_nothing_private():
    """Prefilled GitHub issue links: form fields and system info, special characters encoded, the length cap that keeps the end
    of the log, feature requests and the unknown kind, injectable openers; without details the link carries the last 30 log
    lines; sequences, paths and names never reach any link."""
    # --- form fields and system info
    url = feedback.report_url("bug", "Crash when opening a file")
    assert url.startswith(project.ISSUES_URL + "/new?template=bug_report.yml&")
    q = query(url)
    assert q["template"] == ["bug_report.yml"]
    assert q["title"] == ["[Bug]: Crash when opening a file"]
    assert q["version"] == [project.VERSION]
    assert "Python" in q["os"][0] and "Qt" in q["os"][0]
    assert " " not in url and "\n" not in url

    # --- special characters are encoded
    details = "line one & two = three #4\nnaive caf\u00e9 100% <b>\nlast line"
    q = query(feedback.report_url("bug", "a&b #c?", details))
    assert q["title"] == ["[Bug]: a&b #c?"]
    assert q["logs"] == [details]

    # --- capped, keeps the end
    details = "\n".join(f"line {i:05d} " + "x" * 60 for i in range(2000)) + "\nValueError: the real problem"
    url = feedback.report_url("crash", "s" * 500, details)
    assert len(url) <= feedback.URL_CAP
    logs = query(url)["logs"][0]
    assert logs.endswith("ValueError: the real problem") and "earlier lines left out" in logs
    assert len(query(url)["title"][0]) <= len("[Bug]: ") + feedback.SUMMARY_MAX

    # --- feature link and unknown kind
    url = feedback.feature_url("Dark mode for the map")
    assert url.startswith(project.ISSUES_URL + "/new?template=feature_request.yml&")
    assert query(url)["title"] == ["[Feature]: Dark mode for the map"]
    assert feedback.report_url("feature", "x") == feedback.feature_url("x")
    with pytest.raises(ValueError):
        feedback.report_url("nonsense")

    # --- openers are injectable
    urls = []
    assert feedback.open_report("bug", "hello", opener=urls.append) == urls[0]
    assert feedback.open_feature("idea", opener=urls.append) == urls[1]
    assert "template=bug_report.yml" in urls[0] and "template=feature_request.yml" in urls[1]
    opened = []
    folder = feedback.open_log_folder(opener=opened.append)
    assert folder == feedback.log_folder() and folder.name == "logs" and folder.is_dir()
    assert opened and opened[0].startswith("file:")

    # --- log lines used when there are no details
    with active():
        logging.getLogger("opendna.test").warning("marker line %d", 42)
        q = query(feedback.report_url("bug"))
    assert "marker line 42" in q["logs"][0]
    assert "Start: openDNA" in q["logs"][0]

    # --- only the last thirty log lines
    with active():
        for i in range(100):
            logging.getLogger("opendna.test").warning("entry-%03d", i)
        logs = query(feedback.report_url("bug"))["logs"][0]
    assert "entry-099" in logs and "entry-069" not in logs and "entry-070" in logs

    # --- never any sequences, paths or names
    details = f"Traceback...\nFileNotFoundError: {SECRET_PATH}\nValueError: bad base in {DNA} ({SECRET_PATH!r})"
    with active():
        logging.getLogger("opendna.test").error("opened %s with %s as bobsmith", SECRET_PATH, DNA)
        for url in (feedback.report_url("bug", f"cannot open {SECRET_PATH}", details), feedback.report_url("bug"),
                    feedback.feature_url(f"want {DNA} from my_secret_plasmid.gb")):
            decoded = " ".join(sum(query(url).values(), []))
            for private in (DNA, "bobsmith", "secret", "Lab Notes", "Users", "Documents"):
                assert private not in decoded, private


def test_log_file_and_error_hooks(private_data_dir, monkeypatch):
    """The scrubbed, bounded log file (Qt warnings included, handlers not stacked, diagnostics text) and the error hooks: every
    unhandled error (main thread or worker) logged, one dialog per distinct error with scrubbed text, its report and copy
    buttons, a failing dialog never raising, KeyboardInterrupt not reported, the hooks optional and teardown idempotent."""
    # --- a scrubbed file in the data folder
    with active() as calls:
        path = feedback.log_path()
        assert path == private_data_dir / "logs" / "opendna.log"
        logging.getLogger("opendna.test").info("hello from the test, file %s", SECRET_PATH)
        feedback.teardown_logging()                     # closes (flushes) the file
        text = path.read_text(encoding="utf-8")
    assert f"Start: {project.NAME} {project.VERSION}" in text and "Python" in text and "Qt" in text
    assert "hello from the test" in text and "INFO" in text and "opendna.test" in text
    assert "bobsmith" not in text and "secret" not in text
    assert calls["dialog"] == []

    # --- Qt warnings go to the log
    from PySide6.QtCore import qWarning
    with active():
        qWarning("QWidget: something odd happened in plasmid.gb")
        feedback.teardown_logging()
    text = (private_data_dir / "logs" / "opendna.log").read_text(encoding="utf-8")
    assert "WARNING" in text and "qt:" in text and "something odd happened" in text and "plasmid.gb" not in text

    # --- setting up twice does not stack handlers
    root = logging.getLogger()
    n = len(root.handlers)
    with active():
        feedback.setup_logging(hooks=False, qt_handler=False)
        assert len(root.handlers) == n + 1
    assert len(root.handlers) == n

    # --- rotation keeps a bounded number of files
    with monkeypatch.context() as small:
        small.setattr(feedback, "LOG_MAX_BYTES", 2000)
        with active():
            for i in range(400):
                logging.getLogger("opendna.test").info("filler line %04d %s", i, "y" * 40)
    names = sorted(p.name for p in (private_data_dir / "logs").iterdir())
    assert names == ["opendna.log", "opendna.log.1", "opendna.log.2"]
    assert all(p.stat().st_size <= 2000 + 200 for p in (private_data_dir / "logs").iterdir())

    # --- diagnostics text and copy
    with active() as calls:
        logging.getLogger("opendna.test").info("a log line")
        text = feedback.copy_diagnostics(clipboard=type("C", (), {"setText": lambda self, t: calls["copied"].append(t)})())
    assert text == calls["copied"][0]
    assert f"{project.NAME} {project.VERSION}" in text and "OS:" in text and "a log line" in text
    assert "bobsmith" not in text
    monkeypatch.setenv("OPENDNA_DATA_DIR", str(private_data_dir / "fresh"))
    assert "(the log is empty)" in feedback.diagnostics_text()   # no log file yet in a fresh data folder

    # --- one dialog per distinct error, every occurrence logged (in a fresh data folder: only these errors in its log)
    hooks_dir = private_data_dir / "hooks"
    monkeypatch.setenv("OPENDNA_DATA_DIR", str(hooks_dir))
    with active() as calls:
        assert sys.excepthook is not sys.__excepthook__
        for _ in range(3):
            sys.excepthook(*raised(boom))               # same error three times: one dialog
        assert len(calls["dialog"]) == 1
        sys.excepthook(*raised(boom_other))             # a different error: a second dialog
        assert len(calls["dialog"]) == 2
        sys.excepthook(*raised(boom, "other message"))  # same type and place, other text: still the same error
        assert len(calls["dialog"]) == 2
        feedback.teardown_logging()
    text = (hooks_dir / "logs" / "opendna.log").read_text(encoding="utf-8")
    assert text.count("Unhandled exception") == 5       # every occurrence is logged, only the dialog is once
    assert "Traceback (most recent call last)" in text and "ValueError: x" in text and "KeyError" in text

    # --- the dialog receives only scrubbed text
    with active() as calls:
        sys.excepthook(*raised(boom, f"bad base in {DNA} from {SECRET_PATH}"))
    summary, details = calls["dialog"][0]
    assert summary.startswith("ValueError:")
    for private in (DNA, "bobsmith", "secret", "Lab Notes"):
        assert private not in summary and private not in details

    # --- report button opens the prefilled link, copy button copies
    with active(choice="report") as calls:
        sys.excepthook(*raised(boom, "needs a report"))
        assert len(calls["urls"]) == 1
        url = calls["urls"][0]
        assert "template=bug_report.yml" in url and "needs+a+report" not in url
        assert "ValueError: needs a report" in query(url)["logs"][0]
        assert "ValueError: needs a report" in query(url)["title"][0]
    with active(choice="copy") as calls:
        sys.excepthook(*raised(boom, "copy me"))
        assert len(calls["copied"]) == 1 and "ValueError: copy me" in calls["copied"][0]
        assert calls["urls"] == []

    # --- a failing dialog never raises, KeyboardInterrupt is not reported
    with active(choice="boom") as calls:
        sys.excepthook(*raised(boom, "dialog will fail"))     # must not raise
        assert len(calls["dialog"]) == 1
    with active() as calls:
        sys.excepthook(*raised(lambda: (_ for _ in ()).throw(KeyboardInterrupt())))
        assert calls["dialog"] == []

    # --- errors from other threads
    with active() as calls:
        def work():
            boom("from a worker")
        t = threading.Thread(target=work, name="worker-1")
        t.start()
        t.join()
        assert len(calls["dialog"]) == 1
        assert calls["dialog"][0][0] == "ValueError: from a worker"
        t = threading.Thread(target=work)
        t.start()
        t.join()
        assert len(calls["dialog"]) == 1                    # same error again: not shown again

    # --- hooks are optional, teardown is idempotent
    before = (sys.excepthook, threading.excepthook)
    feedback.setup_logging(hooks=False, qt_handler=False)
    try:
        assert (sys.excepthook, threading.excepthook) == before
    finally:
        feedback.teardown_logging()
        feedback.teardown_logging()
    assert (sys.excepthook, threading.excepthook) == before
