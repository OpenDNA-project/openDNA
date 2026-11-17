"""The toolbar's two number boxes, From and To: the selected stretch (1-based, both ends included). They show what is selected and
changing either one changes the selection at once (main.MainWindow wires them)."""
from PySide6.QtCore import Signal
from PySide6.QtWidgets import QHBoxLayout, QLabel, QLineEdit, QWidget


class RangeBoxes(QWidget):
    edited = Signal()   # the user typed in a box (not when `show_range` sets them)
    submitted = Signal()   # Enter

    def __init__(self):
        super().__init__()
        self.first = QLineEdit(placeholderText="From", clearButtonEnabled=True, maxLength=14)
        self.last = QLineEdit(placeholderText="To", clearButtonEnabled=True, maxLength=14)
        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(4)
        for box, sep in ((self.first, QLabel("–")), (self.last, None)):
            box.setFixedWidth(box.fontMetrics().horizontalAdvance("100,000,000") + 36)   # (longer: scrolls)
            box.textEdited.connect(self.edited)
            box.returnPressed.connect(self.submitted)
            layout.addWidget(box)
            if sep is not None:
                layout.addWidget(sep)
        self.setToolTip("The selected stretch of bases, first and last (1-based, both included); it follows the selection, and typing "
                        "in a box changes the selection. On a circular sequence a last base before the first runs through the "
                        "origin (3700 to 50). With only From filled, the cursor goes after that base.")

    def text(self):
        """What the boxes say as main.parse_range reads it: "101..200", or "101" when only From is filled."""
        a, b = self.first.text().strip(), self.last.text().strip()
        return f"{a}..{b}" if a and b else a or b

    def show_range(self, first, last):
        """Show the selection (numbers, or None for a cursor with nothing selected)."""
        for box, n in ((self.first, first), (self.last, last)):
            box.setText("" if n is None else str(n))

    def clear(self):
        self.show_range(None, None)
