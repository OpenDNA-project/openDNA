"""Where the time goes: named steps timed with `with timed("genome: reads"):`, the last 50
durations of each kept, shown live by Help → Troubleshooting → Performance Timings (PerfWindow).

Costs ~1 µs per timed step, so it stays on all the time.
"""
import time
from collections import deque
from contextlib import contextmanager

from PySide6.QtCore import Qt, QTimer
from PySide6.QtWidgets import (
    QHBoxLayout, QHeaderView, QPushButton, QTableWidget, QTableWidgetItem, QVBoxLayout, QWidget,
)

STATS = {}  # step name -> deque of the last durations, ms


def record(name, seconds):
    STATS.setdefault(name, deque(maxlen=50)).append(seconds * 1000)


@contextmanager
def timed(name):
    t = time.perf_counter()
    try:
        yield
    finally:
        record(name, time.perf_counter() - t)


class Laps:
    """Consecutive steps of one job: `lap("GC")` records the time since the previous lap
    (or since this was made) as "<prefix>: GC"; `done()` records the whole job."""

    def __init__(self, prefix):
        self.prefix = prefix
        self.t0 = self.t = time.perf_counter()

    def __call__(self, name):
        now = time.perf_counter()
        record(f"{self.prefix}: {name}", now - self.t)
        self.t = now

    def done(self):
        record(f"{self.prefix} (total)", time.perf_counter() - self.t0)


class PerfWindow(QWidget):
    """Table of every timed step: last / mean / max ms over its last 50 runs, and how many
    times it ran; refreshed twice a second while open."""
    COLUMNS = ["Step", "Last ms", "Mean ms", "Max ms", "Runs (≤ 50)"]

    def __init__(self, parent=None):
        super().__init__(parent, Qt.Window)
        self.setWindowTitle("Performance timings")
        self.resize(560, 420)
        self.table = QTableWidget(0, len(self.COLUMNS))
        self.table.setHorizontalHeaderLabels(self.COLUMNS)
        self.table.verticalHeader().hide()
        self.table.horizontalHeader().setSectionResizeMode(0, QHeaderView.Stretch)
        self.table.setEditTriggers(QTableWidget.NoEditTriggers)
        reset = QPushButton("Reset", clicked=self.reset)
        buttons = QHBoxLayout()
        buttons.addStretch()
        buttons.addWidget(reset)
        layout = QVBoxLayout(self)
        layout.addWidget(self.table)
        layout.addLayout(buttons)
        self.timer = QTimer(self, interval=500, timeout=self.refresh)

    def showEvent(self, event):
        self.refresh()
        self.timer.start()
        super().showEvent(event)

    def hideEvent(self, event):
        self.timer.stop()
        super().hideEvent(event)

    def reset(self):
        STATS.clear()
        self.refresh()

    def refresh(self):
        # Slowest (by mean) first: that is where to look.
        rows = sorted(STATS.items(), key=lambda kv: -sum(kv[1]) / max(1, len(kv[1])))
        self.table.setRowCount(len(rows))
        for i, (name, times) in enumerate(rows):
            values = [name, f"{times[-1]:.1f}", f"{sum(times) / len(times):.1f}", f"{max(times):.1f}", str(len(times))]
            for j, v in enumerate(values):
                item = QTableWidgetItem(v)
                if j:
                    item.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
                self.table.setItem(i, j, item)
