"""Databases > Database Status (and the Help of each database's tab): one entry per database the program keeps, each showing where its file is, whether it can be
read, how many entries it holds and what is wrong in it (rows with no sequence, letters that are not bases, names used twice).

Every check returns a Status (rows of (label, value, level), level "ok", "warn" or "error"), so the tests read the same
facts the dialog shows; nothing here changes a database."""
import csv
import datetime
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

from PySide6.QtCore import Qt, QUrl
from PySide6.QtGui import QDesktopServices
from PySide6.QtWidgets import QDialog, QDialogButtonBox, QFormLayout, QLabel, QPushButton, QVBoxLayout, QWidget

from uikit import plural
import theme

DNA = set("ACGTN")
IUPAC = set("ACGTRYSWKMBDHVN")
LEVEL_COLORS = theme.Live({"ok": "status_ok", "warn": "status_warn", "error": "status_bad"})
LEVEL_MARKS = {"ok": "✓", "warn": "⚠", "error": "✗"}


@dataclass
class Status:
    title: str
    path: Path | None
    rows: list = field(default_factory=list)   # (label, value, level)

    def add(self, label, value, level="ok"):
        self.rows.append((label, str(value), level))

    @property
    def level(self):
        """The worst level of the rows: what the database is, at a glance."""
        levels = {lvl for _, _, lvl in self.rows}
        return "error" if "error" in levels else "warn" if "warn" in levels else "ok"


def _file_rows(st, path, missing="error", missing_text="No file"):
    """File, size, last change; False when there is no file to read."""
    st.add("File", str(path))
    if not path.exists():
        st.add("State", missing_text, missing)
        return False
    info = path.stat()
    st.add("Size", f"{info.st_size:,} bytes")
    st.add("Last changed", datetime.datetime.fromtimestamp(info.st_mtime).strftime("%Y-%m-%d %H:%M"))
    return True


def _few(items, n=5):
    items = list(items)
    return ", ".join(items[:n]) + ("…" if len(items) > n else "")


def _check_rows(st, rows, columns, prefix=""):
    """Rows of a parts collection: usable parts, types, problems. Returns how many parts are usable."""
    lacking = [c for c in ("name", "type", "sequence") if c not in columns]
    if lacking:
        st.add(f"{prefix}Columns", "missing: " + ", ".join(lacking), "error")
    empty = [r.get("name") or "?" for r in rows if not (r.get("sequence") or "").strip()]
    bad = [r.get("name") or "?" for r in rows if (r.get("sequence") or "").strip() and set((r.get("sequence") or "").strip().upper()) - DNA]
    usable = len(rows) - len(empty) - len(bad)
    st.add(f"{prefix}Parts", f"{usable:,} usable" + (f" of {len(rows):,} rows" if usable != len(rows) else ""), "ok" if usable else "warn")
    if empty:
        st.add(f"{prefix}No sequence", f"{len(empty)} (ignored): {_few(empty)}", "warn")
    if bad:
        st.add(f"{prefix}Not only A, C, G, T, N", f"{len(bad)} (ignored): {_few(bad)}", "warn")
    return usable


def _read_collection(st, path, prefix=""):
    """(rows, columns) of a collection file, or None after adding why it cannot be read."""
    import parts_library
    try:
        _, rows = parts_library._meta_and_rows(path)
    except (OSError, UnicodeDecodeError, csv.Error) as e:
        st.add(f"{prefix}State", f"Cannot be read: {e}", "error")
        return None
    columns = list(rows[0].keys()) if rows else list(parts_library.FIELDS)
    return rows, columns


def _parts_status(title, path, missing, missing_text):
    """One source file on its own (used by the tests on a file of their own)."""
    st = Status(title, path)
    if not _file_rows(st, path, missing, missing_text):
        return st
    read = _read_collection(st, path)
    if read is None:
        return st
    rows, columns = read
    st.add("State", "Readable")
    usable = _check_rows(st, rows, columns)
    kinds = Counter(r.get("type") or "misc_feature" for r in rows if (r.get("sequence") or "").strip())
    if kinds and usable:
        st.add("Types", ", ".join(f"{k} {n}" for k, n in kinds.most_common(8)) + ("…" if len(kinds) > 8 else ""))
    return st


def features_library():
    """The Features Database: every source, the bundled ones (app/parts/*.tsv) and your own (Local), one block of rows each."""
    import parts_library
    folder = parts_library.BUNDLED_DIR
    st = Status("Features Database", folder)
    st.add("Built-in folder", str(folder))
    st.add("Your own features", str(parts_library.user_library()))
    paths = sorted(folder.glob("*.tsv")) if folder.is_dir() else []
    if not paths:
        st.add("Built-in sources", "Missing: Detect Common Features and the Checkup find nothing built in", "error")
    names = {c.id: c.title for c in parts_library.collections()}
    total, kinds, counted = 0, Counter(), 0
    for cid, path in [(p.stem, p) for p in paths] + [(parts_library.YOURS, parts_library.user_library())]:
        prefix = f"{names.get(cid) or parts_library.source_name(cid)}: "
        if cid == parts_library.YOURS and not path.exists():
            st.add(f"{prefix}Features", "none yet: right-click a feature → Add to Features Database, "
                                        "or Import → Features and Tracks → Import Features to Features Database")
            continue
        read = _read_collection(st, path, prefix)
        if read is None:
            continue
        counted += 1
        total += _check_rows(st, *read, prefix=prefix)
        kinds.update(r.get("type") or "misc_feature" for r in read[0] if (r.get("sequence") or "").strip())
        if not parts_library.is_on(cid):
            st.add(f"{prefix}Used", "no: switched off in Databases → Features Database", "warn")
    st.rows.insert(0, ("Sources", f"{counted}, {total:,} usable features in all", "ok" if total else "error"))
    if kinds:
        st.rows.insert(1, ("Types", ", ".join(f"{k} {n}" for k, n in kinds.most_common(8)) + ("…" if len(kinds) > 8 else ""), "ok"))
    return st


def primers():
    from primer_db import data_dir, db
    path = data_dir() / "primers.tsv"
    st = Status("Primer Database", path)
    if not _file_rows(st, path, "ok", "Empty (no file yet): Primers → New Primer… adds the first one"):
        return st
    try:
        with open(path, newline="", encoding="utf-8") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            on_disk = list(reader)
            columns = reader.fieldnames or []
    except (OSError, UnicodeDecodeError, csv.Error) as e:
        st.add("State", f"Cannot be read: {e}", "error")
        return st
    lacking = [c for c in ("name", "sequence") if c not in columns]
    st.add("State", "Readable" if not lacking else "Columns missing: " + ", ".join(lacking), "error" if lacking else "ok")
    loaded = len(db().rows)
    st.add("Primers", f"{len(on_disk):,}", "ok" if on_disk else "warn")
    if loaded != len(on_disk):
        st.add("Loaded in the program", f"{loaded:,} (the file was changed outside the program: restart to read it again)", "warn")
    names = Counter((r.get("name") or "").strip() for r in on_disk)
    twice = [n for n, k in names.items() if n and k > 1]
    if twice:
        st.add("Same name twice", f"{len(twice)}: {_few(twice)}", "warn")
    unnamed = sum(1 for r in on_disk if not (r.get("name") or "").strip())
    if unnamed:
        st.add("No name", plural(unnamed, "row"), "warn")
    empty = [r.get("name") or "?" for r in on_disk if not "".join((r.get("sequence") or "").split())]
    if empty:
        st.add("No sequence", f"{len(empty)}: {_few(empty)}", "warn")
    bad = [r.get("name") or "?" for r in on_disk
           if "".join((r.get("sequence") or "").split()) and set("".join((r.get("sequence") or "").split()).upper()) - IUPAC]
    if bad:
        st.add("Letters that are not bases", f"{len(bad)}: {_few(bad)}", "warn")
    lengths = [len("".join((r.get("sequence") or "").split())) for r in on_disk if (r.get("sequence") or "").strip()]
    if lengths:
        st.add("Lengths", f"{min(lengths)}–{max(lengths)} nt")
    return st


DATABASES = (("&Features Database…", features_library), ("&Primer Database…", primers))


def folder_of(st):
    """The folder a database is kept in (None when it has no file at all)."""
    if st.path is None:
        return None
    return st.path if st.path.is_dir() else st.path.parent


def open_folder(st):
    folder = folder_of(st)
    if folder is not None:
        QDesktopServices.openUrl(QUrl.fromLocalFile(str(folder)))


def status_widgets(st):
    """(heading, form) showing a Status: the overall verdict, then one line per fact, a problem in its colour. A form, not a
    table: nothing to select or edit, no grid, the labels as wide as the longest one, a long value wrapped."""
    head = QLabel(f"<b>{st.title}</b>: <span style='color:{LEVEL_COLORS[st.level]}'>{LEVEL_MARKS[st.level]} "
                  f"{ {'ok': 'OK', 'warn': 'usable, with warnings', 'error': 'problem'}[st.level] }</span>")
    body = QWidget()
    form = QFormLayout(body)
    form.setContentsMargins(0, 0, 0, 0)
    form.setFieldGrowthPolicy(QFormLayout.AllNonFixedFieldsGrow)
    form.setLabelAlignment(Qt.AlignLeft)
    for label, value, level in st.rows:
        name = QLabel(label)
        name.setEnabled(False)   # greyed: the label, the value is what is read
        text = QLabel(value if level == "ok" else f"{LEVEL_MARKS[level]} {value}", wordWrap=True, textFormat=Qt.PlainText)
        text.setTextInteractionFlags(Qt.TextSelectableByMouse)   # a path can be copied
        if level != "ok":
            text.setStyleSheet(f"color: {LEVEL_COLORS[level]}")
        form.addRow(name, text)
    body.rows = form   # for the tests
    return head, body


def show(parent, check):
    """A dialog with the rows of `check()` and a button that opens the database's folder."""
    st = check()
    dialog = QDialog(parent)
    dialog.setWindowTitle(st.title)
    layout = QVBoxLayout(dialog)
    head, facts = status_widgets(st)
    layout.addWidget(head)
    layout.addWidget(facts)
    buttons = QDialogButtonBox(QDialogButtonBox.Close)
    buttons.rejected.connect(dialog.reject)
    if st.path is not None:
        buttons.addButton(QPushButton("Open Folder", clicked=lambda: open_folder(st)), QDialogButtonBox.ActionRole)
    layout.addWidget(buttons)
    dialog.resize(640, dialog.sizeHint().height())
    dialog.status = st   # for the tests
    return dialog
