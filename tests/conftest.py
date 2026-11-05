"""Shared test hooks."""
import os
import sys

import pytest

_status = 0


@pytest.fixture(autouse=True)
def _fresh_window_settings():
    """Each test starts with no remembered window size, View switches or column widths (window.ini): the windows remember
    those between sessions, and one test's choices must not become the next test's defaults."""
    try:
        from primer_db import data_dir
    except ImportError:   # a test that does not use the app
        yield
        return
    (data_dir() / "window.ini").unlink(missing_ok=True)
    yield


def pytest_sessionfinish(session, exitstatus):
    """Let every background job finish before Python shuts down.

    Several tests leave work running in Qt's thread pool (primer matching on the 4.6 Mb genome, off-target searches, ...).
    A thread still inside the Rust core, or releasing Qt objects, while the interpreter tears down can crash the process.
    Waiting here (windows are left open on purpose: closing one can ask about unsaved changes) avoids that. (Neither gc.collect()
    nor processEvents() here: both run the destructors of the abandoned windows, which is what crashes, see pytest_unconfigure.)
    """
    global _status
    _status = int(exitstatus)
    try:
        from PySide6.QtCore import QThreadPool
    except ImportError:   # a machine without Qt: nothing to wait for
        return
    pool = QThreadPool.globalInstance()
    pool.waitForDone(60_000)   # at most a minute: a hung job must not hang the run for ever


def pytest_unconfigure(config):
    """End the process by hand once pytest has printed its summary, keeping its exit status.

    The tests build many windows and never delete them. When the process then shuts down, Qt's teardown (in Python's
    finalisation, or in the DLL detach code at ExitProcess) ends in an access violation (exit code 0xC0000005 /
    0xC0000409) AFTER every test has passed, which fails a CI job with no traceback. Nothing the tests check
    happens during that teardown, so skipping it changes no result; it only stops a crash that is not a test failure from
    being reported as one. Only done when a Qt application was created.
    """
    try:
        from PySide6.QtWidgets import QApplication
    except ImportError:
        return
    if QApplication.instance() is None:
        return
    sys.stdout.flush()
    sys.stderr.flush()
    if sys.platform == "win32":
        # os._exit() still runs the DLLs' detach code (ExitProcess), which is where Qt crashes; TerminateProcess does not.
        import ctypes
        kernel32 = ctypes.windll.kernel32
        kernel32.GetCurrentProcess.restype = ctypes.c_void_p
        kernel32.TerminateProcess.argtypes = [ctypes.c_void_p, ctypes.c_uint]
        kernel32.TerminateProcess(kernel32.GetCurrentProcess(), _status)
    os._exit(_status)
