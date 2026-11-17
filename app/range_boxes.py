"""The toolbar's two number boxes, From and To: the selected stretch (1-based, both ends included). They show what is selected and
changing either one changes the selection at once (main.MainWindow wires them). Numbers are shown with a comma every three digits;
the mouse wheel over a box changes its number (Shift: by 10, Ctrl: by 100)."""
from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import QHBoxLayout, QLabel, QLineEdit, QSizePolicy, QWidget


def grouped(digits):
    """'1234567' (or 1234567) -> '1,234,567' ('' stays '')."""
    return f"{int(digits):,}" if str(digits) else ""


class _Box(QLineEdit):
    def __init__(self, owner, **kw):
        super().__init__(**kw)
        self.owner = owner
        self.turned = 0   # wheel movement not yet worth a step
        self.textEdited.connect(self._typed)

    def value(self):
        digits = "".join(c for c in self.text() if c.isdigit())
        return int(digits) if digits else None

    def _typed(self, text):
        """Only digits, grouped as they are typed, the cursor staying after the same digit."""
        digits = "".join(c for c in text if c.isdigit())
        left = sum(c.isdigit() for c in text[:self.cursorPosition()])   # digits before the cursor stay before it
        self.setText(grouped(digits))
        self.setCursorPosition(len(grouped(digits[:left])))
        self.owner.edited.emit()

    def wheelEvent(self, event):
        # A notch is 120, but a smooth wheel or a touchpad sends smaller steps: they add up (and // would round the upward ones to 0)
        self.turned += event.angleDelta().y()
        steps = int(self.turned / 120)
        self.turned -= steps * 120
        if steps:
            mods = event.modifiers()
            self.owner.scroll(self, steps * (100 if mods & Qt.ControlModifier else 10 if mods & Qt.ShiftModifier else 1))
        event.accept()


class RangeBoxes(QWidget):
    edited = Signal()   # the user typed in a box or turned the wheel (not when `show_range` sets them)
    submitted = Signal()   # Enter

    def __init__(self):
        super().__init__()
        self.length = 0   # of the sequence: the wheel stays between 1 and this
        self.first = _Box(self, placeholderText="From", clearButtonEnabled=True, maxLength=14)
        self.last = _Box(self, placeholderText="To", clearButtonEnabled=True, maxLength=14)
        self.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Preferred)   # (no room added between the boxes)
        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(3)
        for box in (self.first, self.last):
            box.setFixedWidth(box.fontMetrics().horizontalAdvance("100,000,000") + 36)   # (longer: scrolls)
            box.returnPressed.connect(self.submitted)
        for w in (self.first, QLabel("–"), self.last):
            layout.addWidget(w)
        self.setToolTip("The selected stretch of bases, first and last (1-based, both included); it follows the selection, and typing "
                        "in a box changes the selection. Turn the mouse wheel over a box to change its number (Shift: by 10, "
                        "Ctrl: by 100). On a circular sequence a last base before the first runs through the origin (3700 to 50). "
                        "With only From filled, the cursor goes after that base.")

    def scroll(self, box, by):
        """The wheel over `box`: its number by `by`, kept between 1 and the length and (with both boxes filled) on the right side of
        the other one, so turning it never jumps through the origin. An empty box starts from the other one's number, or 1."""
        other = self.last if box is self.first else self.first
        value = box.value() or other.value() or 1
        low, high = 1, max(1, self.length)
        if other.value():
            low, high = (low, other.value()) if box is self.first else (other.value(), high)
        box.setText(grouped(str(min(max(value + by, low), high))))
        self.edited.emit()

    def text(self):
        """What the boxes say as main.parse_range reads it: "101..200", or "101" when only From is filled."""
        a, b = self.first.text().strip(), self.last.text().strip()
        return f"{a}..{b}" if a and b else a or b

    def show_range(self, first, last):
        """Show the selection (numbers, or None for a cursor with nothing selected)."""
        for box, n in ((self.first, first), (self.last, last)):
            box.setText("" if n is None else grouped(n))
