"""Startup page: the welcome message, what is new, recent files, quick actions and an update check on the left; sponsor slots,
a tip and links on the right. Shown at the first start of each day (unless turned off) and from Help → Welcome.

Slots come from sponsors.json next to this file (bundled in the exe): name, text, url,
image. An empty slot shows "Your ad here". The page is labelled as sponsored content
(Canada's Competition Act expects advertising to be identifiable as such) and fetches
nothing: no user data leaves the computer.
"""
import html
import json
from datetime import date
from pathlib import Path

from PySide6.QtCore import QSettings, Qt, QUrl
from PySide6.QtGui import QDesktopServices, QFont, QPixmap
from PySide6.QtWidgets import (
    QCheckBox, QDialog, QDialogButtonBox, QFrame, QHBoxLayout, QLabel, QListWidget, QListWidgetItem, QMenu, QPushButton,
    QVBoxLayout,
)

import feedback
import project
import updater
from citations import CitationsDialog
from primer_db import data_dir
from sourcefiles import pinned_files, recent_files, remove_recent, set_pinned

FILE = Path(__file__).with_name("sponsors.json")
# The demo plasmid: next to the modules in the exe (the spec bundles it), in the repository's demo folder when run from source.
DEMO_NAME = "pBAD30_eilR.gb"
DEMO = next((p for p in (Path(__file__).with_name("demo") / DEMO_NAME, Path(__file__).parent.parent / "demo" / DEMO_NAME)
             if p.exists()), None)
USER_GUIDE_URL = f"{project.REPO_URL}/blob/main/docs/USER_GUIDE.md"
TIPS = [
    "Ctrl+G goes to a position; Ctrl+T makes a feature from the selected bases.",
    "Right-click selected bases to add a forward or reverse primer there.",
    "Shift+click a primer, with one on the other strand selected, to select the PCR product.",
    "The lock button at the left of the toolbar keeps the bases from being changed; features and primers still can be.",
    "Tools > Align Clipboard aligns whatever sequence you copied against the open file.",
    "Ctrl+Shift+A finds common features (promoters, tags, resistance genes...) and annotates them.",
    "Ctrl+Shift+T shows the 3-frame translation under the sequence.",
    "Ctrl+Shift+F opens the Display Filter: choose which feature types are drawn.",
    "Ctrl+wheel makes the letters bigger or smaller; Ctrl+0 brings back the default size.",
    "On a circular sequence, scrolling past the end continues at base 1.",
    "In the Alignment view, Shift+drag selects bases and a click on a track name hides it.",
    "Tools > Plasmid Checkup tells how well each CDS of a plasmid should be expressed.",
    "F1 lists every keyboard shortcut, for the menus and for each panel.",
]


def load():
    """(slots, contact) from sponsors.json; three empty slots if it is missing or broken."""
    try:
        data = json.loads(FILE.read_text(encoding="utf-8"))
        return list(data.get("slots") or [{}, {}, {}])[:3], data.get("contact", "")
    except (OSError, ValueError, AttributeError):
        return [{}, {}, {}], ""


def settings():
    return QSettings(str(data_dir() / "workspace.ini"), QSettings.IniFormat)


def hidden():
    """The user ticked "Don't show at startup" (Help > Welcome still opens the page, and unticking it there undoes this)."""
    return settings().value("sponsors/never", False, type=bool)


def due():
    """Show the page at startup? Once a day, at the first start of the day, unless the user turned it off."""
    return not hidden() and settings().value("sponsors/last_shown", "") != date.today().isoformat()


def welcome_html():
    """The opening message: who made it, that it is free, and where the official download is."""
    by = f"{project.AUTHOR} (<a href='mailto:{project.CONTACT_EMAIL}'>{project.CONTACT_EMAIL}</a>)" if project.CONTACT_EMAIL else project.AUTHOR
    return (f"Version {project.VERSION}. Developed by {by}.<br><br>"
            f"{project.NAME} is <b>free software</b> (AGPL-3.0). The official download is "
            f"<a href='{project.REPO_URL}'>github.com/{project.REPO}</a>.<br><br>"
            f"Bugs and feature requests: <a href='{project.ISSUES_URL}'>GitHub Issues</a> "
            "(Help &gt; Report a Problem).")


def whats_new():
    """The notes of this version, the first time the page opens after an update (not on a first install: nothing was "before");
    None otherwise. Remembers the version, so the notes show once."""
    seen = settings().value("sponsors/seen_version", "")
    settings().setValue("sponsors/seen_version", project.VERSION)
    if seen and seen != project.VERSION:
        return project.WHATS_NEW.get(project.VERSION)
    return None


def tip_of_the_day(offset=0):
    """A different tip each day; `offset` steps through the next ones."""
    return TIPS[(date.today().toordinal() + offset) % len(TIPS)]


def wordmark(dark=False):
    """The logo with the name (app/logo/wordmark*.png, drawn at twice the size shown so it stays sharp on high-DPI screens)."""
    pix = QPixmap(str(Path(__file__).with_name("logo") / ("wordmark_dark.png" if dark else "wordmark.png")))
    pix.setDevicePixelRatio(2.0)
    return pix


def _heading(text):
    label = QLabel(text)
    label.setStyleSheet("color: gray; letter-spacing: 2px; font-size: 9pt")
    return label


class SponsorsDialog(QDialog):
    """The page shown at startup (and from Help > Welcome). The logo on top; on the left the welcome message, what is new, the
    recently opened files (a click opens one; right-click pins or removes one), quick actions and an update check; on the right
    the sponsors when there are any, a tip and links."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle(f"Welcome to {project.NAME}")
        self.setMinimumWidth(820)
        slots, contact = load()
        logo = QLabel(alignment=Qt.AlignCenter, accessibleName=project.NAME)
        logo.setPixmap(wordmark(self.palette().window().color().lightness() < 128))
        title = QLabel(welcome_html(), wordWrap=True, openExternalLinks=True)

        notes = whats_new()
        self.news = QLabel(f"<b>New in {html.escape(project.VERSION)}</b><ul style='margin: 0'>"
                           + "".join(f"<li>{html.escape(n)}</li>" for n in notes) + "</ul>" if notes else "",
                           wordWrap=True, frameShape=QFrame.StyledPanel)
        self.news.setVisible(bool(notes))

        # Recent files: pinned ones first (bold), then the others; the whole path as tooltip (two files of the same name differ by
        # their folder).
        self.recent = QListWidget(contextMenuPolicy=Qt.CustomContextMenu)
        self.recent.setMinimumHeight(150)
        self.recent.itemClicked.connect(self._open_recent)
        self.recent.customContextMenuRequested.connect(self._recent_menu)
        self.fill_recent()

        # Quick actions: each needs the main window behind the page (absent in tests), and leaves the page.
        actions = QHBoxLayout()
        self.quick = {}
        for label, method in (("Open File…", "choose_file"), ("New Sequence…", "new_file"), ("Import from NCBI…", "import_ncbi"),
                              ("Assembly", "show_assembly")):
            button = QPushButton(label, clicked=lambda _=False, m=method: self._run(m))
            button.setEnabled(hasattr(parent, method))
            actions.addWidget(button)
            self.quick[method] = button
        demo = QPushButton("Open the Demo Plasmid", clicked=self.open_demo, toolTip=DEMO_NAME)
        demo.setEnabled(DEMO is not None and hasattr(parent, "open_untitled"))
        actions.addWidget(demo)
        self.quick["demo"] = demo
        actions.addStretch()

        # Update check: the answer shows here, on the page, rather than in a box over it.
        self.update_button = QPushButton("Check for Updates", clicked=self.check_updates)
        self.update_status = QLabel(wordWrap=True, openExternalLinks=False, textFormat=Qt.RichText)
        self.update_status.linkActivated.connect(self._link)
        self.update_found = None
        update_row = QHBoxLayout()
        update_row.addWidget(self.update_button)
        update_row.addWidget(self.update_status, 1)

        left = QVBoxLayout()
        left.addWidget(title)
        left.addWidget(self.news)
        left.addSpacing(8)
        left.addWidget(_heading("RECENT FILES"))
        left.addWidget(self.recent, 1)
        left.addLayout(actions)
        left.addSpacing(8)
        left.addLayout(update_row)

        # Right: sponsors one above the other (none yet: no empty "Your ad here" boxes in a test version), a tip, links.
        sponsored = any(s.get("name") for s in slots)
        self.cards = [self._card(s, contact) for s in slots] if sponsored else []
        right = QVBoxLayout()
        if sponsored:
            right.addWidget(_heading("SPONSORED"))  # advertising, identified as such
            for card in self.cards:
                right.addWidget(card)
            right.addSpacing(8)
        self.tip_offset = 0
        self.tip = QLabel(wordWrap=True, textFormat=Qt.RichText)
        self.tip.linkActivated.connect(self._link)
        self._show_tip()
        right.addWidget(_heading("TIP"))
        right.addWidget(self.tip)
        right.addSpacing(8)
        right.addWidget(_heading("LINKS"))
        links = QLabel(f"<a href='{USER_GUIDE_URL}'>User Guide</a><br><a href='app:shortcuts'>Keyboard Shortcuts</a><br>"
                       "<a href='app:citations'>Citations</a><br><a href='app:report'>Report a Problem</a><br>"
                       f"<a href='{project.RELEASES_URL}'>All Releases</a>", textFormat=Qt.RichText)
        links.linkActivated.connect(self._link)
        right.addWidget(links)
        right.addStretch()

        columns = QHBoxLayout()
        columns.addLayout(left, 3)
        columns.addSpacing(16)
        columns.addLayout(right, 2)

        self.snooze = QCheckBox("Don't show at startup (Help > Welcome opens this page)", checked=hidden())
        buttons = QDialogButtonBox(QDialogButtonBox.Ok)
        buttons.button(QDialogButtonBox.Ok).setText("Continue")
        buttons.accepted.connect(self.accept)
        layout = QVBoxLayout(self)
        layout.addWidget(logo)
        layout.addLayout(columns, 1)
        layout.addWidget(self.snooze)
        layout.addWidget(buttons)

    # ---- recent files ----

    def fill_recent(self):
        self.recent.clear()
        pinned = pinned_files()
        for path in pinned + [p for p in recent_files() if p not in pinned]:
            item = QListWidgetItem(f"{Path(path).name}    {Path(path).parent}")
            item.setData(Qt.UserRole, path)
            item.setToolTip(path + ("\nPinned" if path in pinned else ""))
            if path in pinned:
                font = QFont(item.font())
                font.setBold(True)
                item.setFont(font)
            self.recent.addItem(item)
        empty = not self.recent.count()
        if empty:
            self.recent.addItem(QListWidgetItem("No recent files"))
        self.recent.setEnabled(not empty)

    def _recent_menu(self, pos):
        item = self.recent.itemAt(pos)
        path = item.data(Qt.UserRole) if item else None
        if not path:
            return
        menu = QMenu(self)
        pinned = path in pinned_files()
        menu.addAction("Unpin" if pinned else "Pin to the Top", lambda: self.pin(path, not pinned))
        menu.addAction("Remove from the List", lambda: self.remove(path))
        menu.exec(self.recent.viewport().mapToGlobal(pos))

    def pin(self, path, pinned):
        set_pinned(path, pinned)
        self.fill_recent()

    def remove(self, path):
        remove_recent(path)
        self.fill_recent()

    def _open_recent(self, item):
        """Open the file in the window behind (or a new one), then leave the page."""
        path = item.data(Qt.UserRole)
        if path and hasattr(self.parent(), "open_path"):
            self.accept()
            self.parent().open_path(path)

    # ---- quick actions ----

    def _run(self, method):
        self.accept()
        getattr(self.parent(), method)()

    def open_demo(self):
        """The demo plasmid as a new unsaved file: Save asks where, so the bundled copy is never overwritten."""
        from opendna_core import Document
        self.accept()
        self.parent().open_untitled(Document(str(DEMO)), Path(DEMO_NAME))

    # ---- update check, tip, links ----

    def check_updates(self, fetch=None):
        """Ask GitHub for the latest release (in a thread) and say what came back next to the button."""
        if not project.UPDATES_ENABLED:
            self.update_status.setText(f"Not available yet. New versions: <a href='{project.RELEASES_URL}'>Releases</a>")
            return
        updater.mark_checked()
        self.update_button.setEnabled(False)
        self.update_status.setText("Checking…")

        def got(result):
            try:
                self.update_button.setEnabled(True)
                if result.update is not None:
                    self.update_found = result.update
                    self.update_status.setText(f"<b>{html.escape(result.update.version)}</b> is available. "
                                               "<a href='install:'>Install…</a>")
                elif result.error:
                    self.update_status.setText(f"{html.escape(result.error)} "
                                               f"<a href='{result.page_url or project.RELEASES_URL}'>Releases</a>")
                else:
                    self.update_status.setText(f"Up to date ({project.VERSION} is the latest).")
            except RuntimeError:   # the page was closed while it looked
                pass
        updater.check_in_background(got, fetch=fetch)

    def _show_tip(self):
        self.tip.setText(f"{html.escape(tip_of_the_day(self.tip_offset))}<br><a href='app:next-tip'>Next tip</a>")

    def _link(self, url):
        """Links of the page: app:… acts in the program, install: leaves the page for the update dialog (over the main window);
        the others open in the browser."""
        if url == "app:next-tip":
            self.tip_offset += 1
            self._show_tip()
        elif url == "app:citations":
            CitationsDialog(self).exec()
        elif url == "app:shortcuts" and hasattr(self.parent(), "menuBar"):   # the list is read from the main window's menus
            from shortcuts import ShortcutsDialog
            ShortcutsDialog(self.parent()).exec()
        elif url == "app:report":
            feedback.open_report("bug")
        elif url == "install:":
            window = self.parent() or self
            self.accept()
            hooks = window.update_hooks() if hasattr(window, "update_hooks") else {}
            updater.open_dialog(window, self.update_found, **hooks)
        else:
            QDesktopServices.openUrl(QUrl(url))

    @staticmethod
    def _card(slot, contact):
        """One slot: logo or name, a line of text, a link; or the placeholder."""
        card = QFrame(frameShape=QFrame.StyledPanel)
        card.setMinimumSize(210, 110)
        box = QVBoxLayout(card)
        if not slot.get("name"):
            text = "<b>Your ad here</b>" + (f"<br><br>{contact}" if contact else "")
            card.setStyleSheet("QFrame { border: 2px dashed #9e9e9e; border-radius: 6px; }")
        else:
            text = f"<b>{slot['name']}</b><br>{slot.get('text', '')}"
            if slot.get("url"):
                text += f"<br><a href='{slot['url']}'>{slot['url']}</a>"
            image = FILE.with_name(slot["image"]) if slot.get("image") else None
            if image and image.exists():
                logo = QLabel(alignment=Qt.AlignCenter)
                logo.setPixmap(QPixmap(str(image)).scaled(180, 70, Qt.KeepAspectRatio, Qt.SmoothTransformation))
                box.addWidget(logo)
        label = QLabel(text, alignment=Qt.AlignCenter, wordWrap=True, openExternalLinks=True)
        label.setStyleSheet("border: none")
        box.addWidget(label, 1)
        return card

    def done(self, result):
        settings().setValue("sponsors/never", self.snooze.isChecked())
        super().done(result)


def show_at_startup(parent, force=False):
    """Over the main window without blocking it (a file opened from the command line keeps
    loading behind). Counts as today's showing even if closed at once. force (the --welcome
    option): shown even if already shown today or turned off."""
    if force or due():
        settings().setValue("sponsors/last_shown", date.today().isoformat())
        dialog = SponsorsDialog(parent)
        dialog.open()
        return dialog
    return None
