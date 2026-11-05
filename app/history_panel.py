"""History side panel (toolbar → History): the changes made to the sequence, its features and
the primer database, oldest first, in plain words ("TTT > GAG at 1,204", "Pasted 120 bp + 2
features at 3,001", "Added feature GFP", "Deleted primer oAC1"). Bases typed
or deleted one after the other are one row. Select a row to go back to just before it
(undoing it and everything after), or to redo undone rows (greyed) up to it.

It reads Document.history() (summaries kept by Rust with each undo step), so it is exactly
the Edit → Undo / Redo history, not a separate log.
"""
from PySide6.QtCore import Qt, QTimer
from PySide6.QtGui import QColor, QFont
from PySide6.QtWidgets import QLabel, QListWidget, QListWidgetItem, QPushButton, QVBoxLayout, QWidget

SHOWN = 300  # rows listed (the oldest beyond are still in the undo history)


def _bases(text, count):
    return text + ("…" if count > len(text) else "")


def describe(step):
    """One history step (a Document.history() tuple) in words."""
    pos, removed, n_removed, inserted, n_inserted, added, deleted, note = step
    at = f"{pos + 1:,}"
    if note == "external":   # not in the sequence: a primer database change, its text is in `inserted`
        return inserted
    if note == "reverse complement":
        return "Reverse-complemented the whole sequence"
    if note == "made linear":
        return "Made the sequence linear" + (" (broke the circle elsewhere)" if n_removed else "")
    if note == "made circular":
        return "Made the sequence circular"
    if note == "rotated":
        return "Moved base 1 (broke the circle elsewhere)"
    if note == "across the origin":
        return f"Replaced {n_removed:,} bp across the origin from {at}" + (
            f" with {_bases(inserted, n_inserted)}" if n_inserted else " (deleted)")
    changed = sorted(set(added) & set(deleted))
    words = []
    if n_removed and n_inserted:
        words.append(f"{_bases(removed, n_removed)} > {_bases(inserted, n_inserted)} at {at}" if max(n_removed, n_inserted) <= 12
                     else f"Replaced {n_removed:,} bp with {n_inserted:,} bp at {at}")
    elif n_inserted:
        verb = "Pasted" if added else "Inserted"
        words.append(f"{verb} {n_inserted:,} bp at {at}: {_bases(inserted, n_inserted)}")
    elif n_removed:
        words.append(f"Deleted {n_removed:,} bp at {at}–{pos + n_removed:,}: {_bases(removed, n_removed)}")
    for verb, names in (("Changed", changed), ("Added", [a for a in added if a not in changed]),
                        ("Deleted", [d for d in deleted if d not in changed])):
        if names:
            words.append(f"{verb} feature{'s' if len(names) > 1 else ''} " + ", ".join(names[:4])
                         + ("…" if len(names) > 4 else ""))
    return "; ".join(words) or "No change"


def _typing(step):
    pos, removed, n_removed, inserted, n_inserted, added, deleted, note = step
    if note or added or deleted:
        return None
    if n_removed == 0 and n_inserted == 1:
        return "type"
    if n_inserted == 0 and n_removed == 1:
        return "delete"
    return None


def group(steps):
    """Steps -> rows [(first step, last step, text)]: bases typed one after the other, or
    deleted one after the other (Backspace or Delete), become one row."""
    rows = []  # [first, last, kind, pos, text, start, end]
    for i, step in enumerate(steps):
        kind, pos = _typing(step), step[0]
        last = rows[-1] if rows else None
        if kind == "type" and last and last[2] == "type" and pos == last[6]:
            last[1], last[4], last[6] = i, last[4] + step[3], pos + 1
        elif kind == "delete" and last and last[2] == "delete" and pos in (last[5] - 1, last[5]):
            # Backspace eats leftwards (pos one before), Delete keeps the same pos.
            last[1] = i
            last[4] = step[1] + last[4] if pos < last[5] else last[4] + step[1]
            last[5] = pos
        else:
            rows.append([i, i, kind, pos, step[3] if kind == "type" else step[1], pos, pos + 1])
    out = []
    for first, last, kind, pos, text, start, _ in rows:
        if kind and last > first:
            what = "Typed" if kind == "type" else "Deleted"
            out.append((first, last, f"{what} {len(text):,} bp at {start + 1:,}: {_bases(text[:30], len(text))}"))
        else:
            out.append((first, last, describe(steps[first])))
    return out


def file_group(steps):
    """Like `group`, for the steps of a file of several records, [(record, step)]: rows [(first, last, text)] over the whole list,
    each text starting with its record's name; only steps of the same record, one after the other, are put in one row."""
    out, i = [], 0
    while i < len(steps):
        j = i
        while j + 1 < len(steps) and steps[j + 1][0] == steps[i][0]:
            j += 1
        out += [(i + a, i + b, f"{steps[i][0]}: {text}") for a, b, text in group([s for _, s in steps[i:j + 1]])]
        i = j + 1
    return out


NOW = "▶ Current sequence"


class HistoryPanel(QWidget):
    def __init__(self, editor):
        super().__init__()
        self.editor = editor
        self.list = QListWidget()
        self.list.setWordWrap(True)
        self.list.currentRowChanged.connect(self._update_button)
        self.list.itemDoubleClicked.connect(lambda _: self.go())
        # The "current sequence" row can be dragged up (go back) or down (redo).
        self.list.setDragDropMode(QListWidget.InternalMove)
        self.list.model().rowsMoved.connect(self._moved)
        self.button = QPushButton("Go back to before this change", clicked=self.go)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(4, 4, 4, 4)
        layout.addWidget(QLabel("<b>History</b> (oldest first; greyed = undone).<br>"
                                "Drag ▶ Current sequence up or down to go back or forward."))
        layout.addWidget(self.list, 1)
        layout.addWidget(self.button)
        editor.changed.connect(self.refresh)
        self.refresh()

    def theme_changed(self):
        """Light / dark switch (theme.refresh_views): the greyed rows and the current-sequence box take the new colours."""
        self.refresh()

    def refresh(self):
        ed = self.editor
        if getattr(ed, "file_history", None) is not None:   # a file of several records: one history for all of them
            undo, redo = ed.file_history.steps(ed.file_histories())
            rows_undo, rows_redo = file_group(undo), file_group(redo)
        else:
            undo, redo = ed.doc.history()
            rows_undo, rows_redo = group(undo), group(redo)
        self.rows = []  # per list row: ("undo" | "redo" | "now", step count to undo / redo)
        self.list.clear()
        grey = QColor(self.palette().color(self.palette().ColorRole.PlaceholderText))
        still = Qt.ItemIsEnabled | Qt.ItemIsSelectable  # not draggable
        opened = QListWidgetItem(f"Opened {self.editor.path.name}")
        italic = QFont(self.list.font())
        italic.setItalic(True)
        opened.setFont(italic)
        opened.setFlags(still)
        self.list.addItem(opened)
        self.rows.append(("undo", len(undo)) if undo else ("start", 0))  # back to the file as opened
        for first, _, text in rows_undo[-SHOWN:]:
            item = QListWidgetItem(text)
            item.setFlags(still)
            self.list.addItem(item)
            self.rows.append(("undo", len(undo) - first))  # undo this row and all after it
        now = QListWidgetItem(NOW)
        bold = QFont(self.list.font())
        bold.setBold(True)
        now.setFont(bold)
        # A tinted box (selection colour, see-through: fits light and dark themes), so the
        # marker stands out and looks grabbable.
        box = QColor(self.palette().color(self.palette().ColorRole.Highlight))
        box.setAlpha(70)
        now.setBackground(box)
        now.setFlags(Qt.ItemIsEnabled | Qt.ItemIsDragEnabled)
        now.setToolTip("Drag up to go back to before a change, down to redo")
        self.list.addItem(now)
        self.now_row = self.list.count() - 1
        self.rows.append(("now", 0))
        for _, last, text in rows_redo:
            item = QListWidgetItem(text)
            item.setForeground(grey)
            item.setFlags(still)
            self.list.addItem(item)
            self.rows.append(("redo", last + 1))  # redo up to and including this row
        self.list.scrollToItem(now)
        self._update_button(self.list.currentRow())

    def _moved(self, _parent, start, _end, _dest, row):
        """The current-sequence row was dropped before old row `row`: undo the rows it
        passed going up, or redo those it passed going down."""
        if start != self.now_row:
            return
        if row < self.now_row:
            kind, count = self.rows[row]
            steps = count if kind == "undo" else 0
        elif row > self.now_row + 1:
            kind, count = self.rows[row - 1]
            steps = -count if kind == "redo" else 0
        else:
            steps = 0
        # After Qt finishes the move: then rebuild the list from the real history.
        QTimer.singleShot(0, lambda: self._apply(steps) if steps else self.refresh())

    def _apply(self, count):
        """Undo `count` steps (redo when negative): of the file when it has several records (that may show another one: the
        panel is replaced, so nothing here touches it afterwards)."""
        ed = self.editor
        if getattr(ed, "file_history", None) is not None:
            ed.file_undo_steps(count)
        else:
            ed.undo_steps(count)

    def _update_button(self, row):
        kind = self.rows[row][0] if 0 <= row < len(self.rows) else "now"
        self.button.setEnabled(kind in ("undo", "redo"))
        self.button.setText("Redo up to this change" if kind == "redo" else
                            "Back to the file as opened" if row == 0 else "Go back to before this change")

    def go(self):
        row = self.list.currentRow()
        if not 0 <= row < len(self.rows) or self.rows[row][0] not in ("undo", "redo"):
            return
        kind, count = self.rows[row]
        count = count if kind == "undo" else -count
        if getattr(self.editor, "file_history", None) is not None:
            QTimer.singleShot(0, lambda: self._apply(count))
        else:
            self.editor.undo_steps(count)
