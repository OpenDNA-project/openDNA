"""The PCR tab: what two primers amplify, in three columns like KLD. Left: the template (the open sequence), editable, where a region
can be selected. Middle: the primers: two opposing primers that bind the template, or primer pairs designed for the selected
region. Right: the product, read only. The last choice made (a primer, or a designed pair) is the one shown on the right."""
from types import SimpleNamespace

from PySide6.QtCore import Qt, QTimer
from PySide6.QtGui import QColor
from PySide6.QtWidgets import (
    QAbstractItemView, QComboBox, QDoubleSpinBox, QFormLayout, QHBoxLayout, QInputDialog, QLabel, QMessageBox, QPushButton,
    QSpinBox, QSplitter, QTableWidget, QTableWidgetItem, QVBoxLayout, QWidget,
)

from html import escape

from opendna_core import dimer

import theme
from primer_db import add_primers
from primer_design import DEFAULT_MAX_LEN, DEFAULT_MIN_LEN, DEFAULT_TM, DEFAULT_WINDOW, design_region, salt
from primer_panel import (MARKS, MAX_NAME, ReportView, conditions, is_degenerate, limit, primer_checks, primer_problem, warm_genome,
                          worst)
from tools import ToolDialog, pcr_product, product_range
from uikit import fit_columns, plural, right_aligned, run_busy


class PcrPage(ToolDialog):
    COLUMNS = ["Forward 5'→3'", "Reverse 5'→3'", "Tm F / R °C", "Product bp", "Penalty", "Notes"]

    def __init__(self, parent, editor):
        super().__init__(parent)
        from sequence_view import SequenceView
        self.editor, self.pairs, self.product_doc = editor, [], None
        self.version = None   # the editor's version the product on the right was made from
        self.temp = None   # the designed pair clicked in the table: {"fwd", "rev" (sites), "checks", "row"}
        self.fwd_sites = self.rev_sites = []
        self.picked = {}   # strand -> the primer clicked on the template (any primer of the sequence, listed beside the file's own)
        self.right = None
        self.setWindowTitle("PCR")
        # -- the middle: two primers that bind, or pairs designed for the selection
        self.fwd, self.rev = QComboBox(), QComboBox()
        for combo in (self.fwd, self.rev):
            combo.currentIndexChanged.connect(self._pair_chosen)
        self.tm = QDoubleSpinBox(minimum=45, maximum=75, value=DEFAULT_TM, decimals=1, suffix=" °C")
        self.min_len = QSpinBox(minimum=12, maximum=35, value=DEFAULT_MIN_LEN)
        self.max_len = QSpinBox(minimum=12, maximum=40, value=DEFAULT_MAX_LEN)
        self.min_len.valueChanged.connect(lambda v: self.max_len.setValue(max(v, self.max_len.value())))   # shortest <= longest
        self.max_len.valueChanged.connect(lambda v: self.min_len.setValue(min(v, self.min_len.value())))
        self.window = QSpinBox(minimum=0, maximum=2000, value=DEFAULT_WINDOW, suffix=" bp")
        self.window.setToolTip("How far outside the selected region the primers may sit, on each side. 0: the product is exactly "
                               "the selected region.")
        self.target = QLabel()
        lengths = QHBoxLayout()
        lengths.addWidget(self.min_len)
        lengths.addWidget(QLabel("to"))
        lengths.addWidget(self.max_len)
        form = QFormLayout()
        form.addRow("Forward primer", self.fwd)
        form.addRow("Reverse primer", self.rev)
        design_form = QFormLayout()
        design_form.addRow("Region", self.target)
        design_form.addRow("Tm", self.tm)
        design_form.addRow("Length", lengths)
        design_form.addRow("Search window", self.window)
        self.table = QTableWidget(0, len(self.COLUMNS))
        self.table.setHorizontalHeaderLabels(self.COLUMNS)
        self.table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.table.setSelectionMode(QAbstractItemView.SingleSelection)
        self.table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        fit_columns(self.table, stretch=5)
        self.table.currentCellChanged.connect(lambda row, *_: self._design_chosen(row))
        self.destroyed.connect(lambda *_, ed=editor: self._forget(ed))   # the tab is closed: its primers go with it
        self.details = ReportView(16)
        self.details.anchorClicked.connect(self._toggle_check)
        self.open_checks = set()   # the collapsed rows that were opened (by their link)
        self.shown_checks = None   # what the details show: (checks of the forward and reverse primer, heterodimer), or None
        self.info = QLabel(wordWrap=True)
        self.open_button = QPushButton("Open Product", clicked=self.accept, enabled=False,
                                       toolTip="The product as a new sequence")
        self.save_button = QPushButton("Save Pair to Primer Database…", clicked=self.save, enabled=False)
        buttons = QHBoxLayout()
        buttons.addWidget(self.open_button)
        buttons.addWidget(self.save_button)
        buttons.addStretch()
        middle = QVBoxLayout()
        middle.setContentsMargins(0, 0, 0, 0)
        middle.addWidget(QLabel("<b>Primers</b>"))
        middle.addLayout(form)
        middle.addWidget(QLabel("<b>or design for the selected region</b>"))
        middle.addLayout(design_form)
        middle.addWidget(theme.mark_primary(QPushButton("Design", clicked=self.design)))
        middle.addWidget(self.table, 1)
        middle.addWidget(self.details)
        middle.addWidget(self.info)
        middle.addLayout(buttons)
        # -- left and right
        self.left_view = SequenceView(editor)
        self.right_box = QVBoxLayout()
        self.right_box.setContentsMargins(0, 0, 0, 0)
        self.right_box.addWidget(QLabel("<b>Product</b> (read only)"))
        left = QWidget()
        lbox = QVBoxLayout(left)
        lbox.setContentsMargins(0, 0, 0, 0)
        lbox.addWidget(QLabel("<b>Template</b> (editable: select a region to design for)"))
        lbox.addWidget(self.left_view, 1)
        mid, right = QWidget(), QWidget()
        mid.setLayout(middle)
        right.setLayout(self.right_box)
        split = QSplitter(Qt.Horizontal)
        for w in (left, mid, right):
            split.addWidget(w)
        split.setSizes([420, 420, 420])
        QVBoxLayout(self).addWidget(split)
        editor.selection_changed.connect(self._set_target)
        editor.changed.connect(self._template_edited)
        editor.primer_clicked.connect(self._primer_clicked)
        self._set_target()
        self.refresh()
        if parent.__dict__.pop("design_on_open", False) and editor.end > editor.start:   # right-click > Design Primers
            QTimer.singleShot(0, self.design)

    # -- the right: the product of whatever was chosen last
    def _show(self, doc, text, error=False):
        """Show `doc` (None: `text` says why there is none) on the right."""
        from assembly_panel import FragmentEditor, ReadOnlyView
        self.product_doc = doc
        self.open_button.setEnabled(doc is not None)
        self.info.setText(text)
        if self.right is not None:
            self.right.hide()   # (reparented while shown, a widget flashes as a window of its own)
            self.right_box.removeWidget(self.right)
            self.right.deleteLater()
        if doc is None:
            self.right = QLabel(text, wordWrap=True, alignment=Qt.AlignCenter)
            self.right.setEnabled(not error)
        else:
            self.right = ReadOnlyView(FragmentEditor(doc, "PCR product"))
        self.right_box.addWidget(self.right, 1)

    def _from_pair(self):
        """The product of the two primers chosen in the lists: the lists are the one source of what is shown."""
        pair = self.pair()
        if not pair:
            self.save_button.setEnabled(False)
            self._show_details(None)
            return self._show(None, "Pick a forward and a reverse primer, or design a pair for the selected region.")
        f, r = pair
        self.version = self.editor.version
        span = product_range(self.editor, f, r)
        designed = self.temp is not None and f is self.temp["fwd"] and r is self.temp["rev"]   # (the pair clicked in the table)
        self.save_button.setEnabled(designed)
        self._show_details(self.temp["checks"] if designed else None, heterodimer(f.sequence, r.sequence))
        if not span:
            return self._show(None, "No product: the reverse primer is not downstream.", True)
        try:
            doc = pcr_product(self.editor, f, r)
        except (OSError, ValueError) as e:
            return self._show(None, str(e), True)
        if f.temporary or r.temporary:   # a designed primer is in no file or database: put it on the product to be drawn
            from assembly import transfer_primers
            transfer_primers(doc, [SimpleNamespace(fwd=(f.label, f.sequence), rev=(r.label, r.sequence))])
        size = span[1] - span[0] + f.tail + r.tail
        self._show(doc, f"{size:,} bp (template {span[0] + 1:,}–{(span[1] - 1) % self.editor.length + 1:,}"
                        f" + {f.tail + r.tail} bp of tails)")

    # -- the lists of primers: the one source of what is chosen
    @staticmethod
    def _key(site):
        return (site.label, site.start, site.strand)

    def _selected(self):
        """(key of the forward primer chosen, key of the reverse primer chosen): None for a list with nothing chosen."""
        return tuple(self._key(sites[combo.currentIndex()]) if combo.currentIndex() >= 0 else None
                     for combo, sites in ((self.fwd, self.fwd_sites), (self.rev, self.rev_sites)))

    def _populate(self, keep=None):
        """Fill the lists: the file's own primers (not the database's), then the pair clicked in the table while there is one.
        The choice (`keep`: the keys of a pair, else what is chosen now) stays when it is still listed; signals are blocked."""
        keep = keep or self._selected()
        sites = [s for s in self.editor.primer_sites if s.from_file and not s.temporary]
        listed = {self._key(s) for s in sites}
        extra = [s for s in self.picked.values() if self._key(s) not in listed]
        extra += [self.temp["fwd"], self.temp["rev"]] if self.temp else []
        self.fwd_sites = [s for s in sites + extra if s.strand > 0]
        self.rev_sites = [s for s in sites + extra if s.strand < 0]
        for combo, found, i in ((self.fwd, self.fwd_sites, 0), (self.rev, self.rev_sites, 1)):
            combo.blockSignals(True)
            combo.clear()
            for s in found:
                mm = f", {plural(s.mismatches, 'mismatch')}" if s.mismatches else ""
                combo.addItem(f"{s.label}  ({s.start + 1:,}{mm})" + ("  · designed" if s.temporary else
                              "  · clicked" if s in self.picked.values() else ""))
            keys = [self._key(s) for s in found]
            combo.setCurrentIndex(keys.index(keep[i]) if keep and keep[i] in keys else -1)
            combo.blockSignals(False)

    def refresh(self):
        """The primers of the template came or changed: list them again, keeping the choice. Only a pair that is not the one
        shown is drawn: the primers come back a moment after each design (they are searched in the background), and that must
        not make a new PCR of what was chosen."""
        before = self._selected()
        self._populate()
        if self._selected() == (None, None):
            self._pick_default()
        if self.product_doc is None or self._selected() != before or self.version != self.editor.version:
            self._from_pair()

    def _template_edited(self):
        """The template changed: the designed pairs and the primers clicked were for the sequence as it was."""
        if self.pairs or self.picked:
            self.picked = {}
            self._clear_designs()
            self.refresh()

    def _primer_clicked(self, key):
        """A primer clicked (or shift-clicked) on the template: it is chosen in the list of its strand (and added to it when the
        list does not have it: a primer of the database only). Click one on each strand to make a pair."""
        site = self.editor.site_by_key.get(key)
        if site is None:
            return
        combo = self.fwd if site.strand > 0 else self.rev
        if self._key(site) not in self._listed(site):
            self.picked[site.strand] = site
            self._populate()
        combo.setCurrentIndex(self._listed(site).index(self._key(site)))   # (draws through _pair_chosen)

    def _listed(self, site):
        """The keys listed in the list of this site's strand."""
        return [self._key(s) for s in (self.fwd_sites if site.strand > 0 else self.rev_sites)]

    def _pick_default(self):
        """The pair of the file's own primers giving the shortest product (a designed pair is never picked for you)."""
        best = None
        for i, f in enumerate(self.fwd_sites):
            for j, r in enumerate(self.rev_sites):
                span = None if f.temporary or r.temporary else product_range(self.editor, f, r)
                if span and (best is None or span[1] - span[0] < best[0]):
                    best = (span[1] - span[0], i, j)
        for combo, i in ((self.fwd, 1), (self.rev, 2)):
            combo.blockSignals(True)
            combo.setCurrentIndex(best[i] if best else -1)   # (none: nothing chosen, not the first of the list)
            combo.blockSignals(False)

    def pair(self):
        i, j = self.fwd.currentIndex(), self.rev.currentIndex()
        if i < 0 or j < 0:
            return None
        return self.fwd_sites[i], self.rev_sites[j]

    def _pair_chosen(self, *_):
        """A primer chosen in a list. The designed pair clicked in the table stays in the lists, chosen, until another pair is
        clicked or both lists are moved off it: then it goes (and its row is let go)."""
        pair = self.pair()
        if self.temp is not None and not (pair and any(s.temporary for s in pair)):
            self.temp = None
            self.table.blockSignals(True)
            self.table.clearSelection()
            self.table.setCurrentCell(-1, -1)
            self.table.blockSignals(False)
            self._forget(self.editor)
            self._populate()
        self._from_pair()

    # -- the pairs designed for the selected region
    def _set_target(self):
        ed = self.editor
        self.target.setText(f"{ed.start + 1:,}–{ed.end:,} ({ed.end - ed.start:,} bp, the selection)" if ed.end > ed.start
                            else "none: select a region on the template")

    def design(self):
        ed = self.editor
        if ed.end <= ed.start:
            QMessageBox.warning(self, "No region", "Select the region to amplify on the template first.")
            return
        self._clear_designs()
        self._populate()
        args = dict(start=ed.start, end=ed.end, tm_target=self.tm.value(), min_len=self.min_len.value(), max_len=self.max_len.value(),
                    window=self.window.value())
        # (thousands of candidates: in a thread behind a busy dialog, so the window is drawn and stays alive meanwhile)
        self.pairs = run_busy(self, "Designing primers…", lambda: design_region(ed.doc, **args))
        self.table.setRowCount(len(self.pairs))
        # The genome search (every chromosome) first, with a progress dialog; the checks below then read its results.
        warm_genome(self, ed, [p[k]["seq"] for p in self.pairs for k in ("fwd", "rev")])
        for i, p in enumerate(self.pairs):
            values = [p["fwd"]["seq"], p["rev"]["seq"], f"{p['fwd']['tm']:.1f} / {p['rev']['tm']:.1f}",
                      f"{p['size']:,}", f"{p['penalty']:.1f}", "; ".join(p["notes"])]
            # The pair's mark: the worst of all checks of both primers (details below the table).
            level = worst(primer_checks(ed, p["fwd"]["seq"], "cached") + primer_checks(ed, p["rev"]["seq"], "cached"))
            mark, color = MARKS["warn" if level == "ok" and p["notes"] else level]
            values[5] = f"{mark} {values[5]}".strip()
            for c, v in enumerate(values):
                item = QTableWidgetItem(v)
                if c in (2, 3, 4):
                    right_aligned(item)
                if c == 5:
                    item.setForeground(QColor(color))
                self.table.setItem(i, c, item)
        if self.pairs:
            self.table.selectRow(0)
        else:
            self._from_pair()   # (back to what the lists say)
            self.details.setText("No pair found: widen the search window, the length range or the Tm.")

    def _show_details(self, checks, hetero=None):
        """The bottom section: the checks of each primer side by side, then the pair's heterodimer."""
        self.shown_checks = (checks, hetero)
        self._draw_details()

    def _draw_details(self):
        """One table, a row per check, so a check is on the same line for both primers: name, forward, reverse; the pair's
        heterodimer under them, over both columns."""
        checks, hetero = self.shown_checks
        rows = []
        if checks:
            rows.append("<tr><td></td><td><b>Forward</b></td><td><b>Reverse</b></td></tr>")
            for i, (f, r) in enumerate(zip(*checks)):
                rows.append(f"<tr><td><b>{escape(f[0])}</b>&nbsp;&nbsp;</td>"
                            f"<td>{self._cell(f'F{i}', *f)}</td><td>{self._cell(f'R{i}', *r)}</td></tr>")
        if hetero:
            rows.append(f"<tr><td><b>{escape(hetero[0])}</b>&nbsp;&nbsp;</td><td colspan='2'>{self._cell('H', *hetero)}</td></tr>")
        self.details.setText("<table cellspacing='2' cellpadding='1'>" + "".join(rows) + "</table>")

    def _cell(self, key, name, level, text):
        """A check's mark and text. The lists of places a primer also matches are collapsed to what comes before the list; the
        link at their end (the cell's `key`) opens the list and closes it again."""
        mark, color = MARKS[level]
        body = escape(text)
        parts = text.split("; ")
        exact = [x for x in parts if "mismatch" not in x]
        if name == "Binds here" and len(exact) < len(parts):   # the sites with mismatches, after the exact one(s)
            shorter = escape("; ".join(exact) or "no exact site") + " · and places with mismatches"
        elif name.startswith("Non-specific") and level != "ok" and ": " in text:
            shorter = escape(text.split(": ")[0])
        else:
            shorter = None
        if shorter is not None:
            body = f"{body} <a href='{key}'>◂ hide</a>" if key in self.open_checks else f"{shorter} <a href='{key}'>▸ show</a>"
        return f"<b><font color='{color}'>{mark}</font></b>&nbsp;{body}"

    def _toggle_check(self, url):
        key = url.toString().lstrip("#")
        self.open_checks ^= {key}
        if self.shown_checks is not None:
            self._draw_details()

    def _clear_designs(self):
        self.pairs, self.temp = [], None
        self.table.setRowCount(0)
        self._forget(self.editor)
        self.save_button.setEnabled(False)

    @staticmethod
    def _forget(editor):
        """Stop drawing the designed primers (they are not in the database)."""
        try:
            editor.set_extra([])
        except RuntimeError:
            pass   # the editor is already being deleted with its window

    def temporary(self, row):
        """The pair on `row` as temporary primers: [(name, sequence)], named pair 1 F / pair 1 R..."""
        p = self.pairs[row]
        return [(f"pair{row + 1}_F", p["fwd"]["seq"]), (f"pair{row + 1}_R", p["rev"]["seq"])]

    def _design_chosen(self, row):
        """A designed pair clicked: its two primers are added to the lists and chosen there (the product, the checks and the
        drawing on the template follow from the lists)."""
        if not 0 <= row < len(self.pairs):
            return
        p = self.pairs[row]
        names = self.temporary(row)
        sites = []
        for (name, seq), key, strand in zip(names, ("fwd", "rev"), (1, -1)):
            c = p[key]
            sites.append(SimpleNamespace(key=name, label=name, sequence=seq, strand=strand, start=c["start"], end=c["end"],
                                         mismatches=0, tail=0, temporary=True))
        checks = [primer_checks(self.editor, p[k]["seq"], "cached") for k in ("fwd", "rev")]
        self.temp = {"fwd": sites[0], "rev": sites[1], "checks": checks, "row": row}
        self.editor.set_extra([], names)
        self._populate(keep=(self._key(sites[0]), self._key(sites[1])))
        self._from_pair()

    def save(self, base=None):
        """Save Pair to Primer Database: the designed pair on the chosen row, as <name>_F and <name>_R."""
        pair = self.pair()
        if self.temp is None or not pair or pair[0] is not self.temp["fwd"] or pair[1] is not self.temp["rev"]:
            return None   # (the designed pair is not the one chosen in the lists)
        row = self.temp["row"]
        if base is None:
            base, ok = QInputDialog.getText(self, "Save pair", "Name (saved as <name>_F and <name>_R):",
                                            text=f"{self.editor.doc.name}_{self.pairs[row]['fwd']['start'] + 1}")
            if not ok or not base.strip():
                return None
        p = self.pairs[row]
        names = [f"{base.strip()}_F", f"{base.strip()}_R"]
        problem = primer_problem(names[0], p["fwd"]["seq"]) or primer_problem(names[1], p["rev"]["seq"])
        if problem or len(names[0]) > MAX_NAME:   # a tab in the name, or one too long for the database's name box
            QMessageBox.warning(self, "Save pair", problem or f"The name is too long (at most {MAX_NAME - 2} characters).")
            return None
        if add_primers(self.editor, [(names[0], p["fwd"]["seq"]), (names[1], p["rev"]["seq"])], "Saved", in_file=True)["conflict"]:
            QMessageBox.warning(self, "Name in use", "A primer with that name and another sequence exists.")
            return None
        self._clear_designs()   # now regular primers, under their new names
        self.shown_checks = None
        self.details.setText(f"Saved {names[0]} and {names[1]} to the primer database.")
        return names

    def product(self):
        """The product on the right as a Document (Open Product)."""
        if self.product_doc is None:
            raise ValueError("There is no product to open.")
        return self.product_doc


def heterodimer(forward, reverse):
    """The ("Heterodimer (F × R)", level, text) check of a pair: how well the two primers pair with each other (ΔG, kcal/mol, anywhere
    and at the 3′ ends: a 3′ pairing is what makes primer dimers), judged by the self-dimer limits of Tools > Settings."""
    name = "Heterodimer (F × R)"
    if is_degenerate(forward) or is_degenerate(reverse):
        return (name, "info", "not calculated for a degenerate primer")
    anywhere, end = dimer(forward.upper(), reverse.upper(), **salt(conditions()))
    level = "bad" if end < limit("dimer3_bad") or anywhere < limit("dimer_bad") else "warn" if end < -5 or anywhere < -9 else "ok"
    return (name, level, f"ΔG {anywhere:.1f} (3′ end {end:.1f}) kcal/mol")
