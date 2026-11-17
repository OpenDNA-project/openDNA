"""One openDNA process per user (per data folder): a second launch (a double-click on a file in the file manager) hands its files to
the one that is running over a local socket and ends; the running one opens them as View > Open Files in Tabs says."""
import hashlib
import json

from PySide6.QtCore import QObject, Signal
from PySide6.QtNetwork import QLocalServer, QLocalSocket

from primer_db import data_dir

MAX_MESSAGE = 1_000_000   # bytes: a list of paths is far smaller


def _name():
    # (the data folder is in the name: a run with a data folder of its own, a test or a smoke run, is a program of its own)
    return "openDNA-" + hashlib.sha1(str(data_dir()).encode("utf-8")).hexdigest()[:12]


def send_to_running(files, timeout=1000, answer_timeout=10_000):
    """True if a running openDNA took `files` (full paths; it opens them and shows its window), False if none is running (or it did
    not answer: then this launch starts a program of its own, as before). Waits for the answer: a launch that ended first could
    take the message with it."""
    socket = QLocalSocket()
    socket.connectToServer(_name())
    if not socket.waitForConnected(timeout):
        return False
    socket.write(json.dumps({"files": list(files)}).encode("utf-8") + b"\n")
    socket.waitForBytesWritten(timeout)
    taken = socket.waitForReadyRead(answer_timeout) and bytes(socket.readAll()).startswith(b"ok")
    socket.disconnectFromServer()
    return taken


class Listener(QObject):
    """Receives what other launches hand over: `files` is emitted with their list of paths (maybe empty: just show the window)."""
    files = Signal(list)

    def __init__(self, parent=None):
        super().__init__(parent)
        QLocalServer.removeServer(_name())   # the leftover socket of a run that crashed (a live one was found by send_to_running)
        self.server = QLocalServer(self)
        self.server.setSocketOptions(QLocalServer.UserAccessOption)   # only this user
        self.listening = self.server.listen(_name())
        self.server.newConnection.connect(self._accept)

    def _accept(self):
        while self.server.hasPendingConnections():
            socket = self.server.nextPendingConnection()
            socket.setProperty("buffer", b"")
            socket.readyRead.connect(lambda socket=socket: self._read(socket))
            socket.disconnected.connect(socket.deleteLater)
            if socket.bytesAvailable():   # (what was sent before we got to accept it: readyRead will not say)
                self._read(socket)

    def _read(self, socket):
        buffer = bytes(socket.property("buffer")) + bytes(socket.readAll())
        socket.setProperty("buffer", buffer)
        if b"\n" not in buffer:
            if len(buffer) > MAX_MESSAGE:
                socket.abort()
            return
        socket.setProperty("buffer", b"")
        try:
            paths = json.loads(buffer.split(b"\n", 1)[0].decode("utf-8"))["files"]
        except (ValueError, KeyError, TypeError):
            return
        if isinstance(paths, list) and all(isinstance(p, str) for p in paths):
            socket.write(b"ok\n")   # (before opening: the sender is only waiting for this)
            socket.flush()
            self.files.emit(paths)
