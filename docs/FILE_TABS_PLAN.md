# File tabs: several files in one window

Status: plan agreed, no code written yet.

## Goal

Today one window holds one file; opening a second file opens a second window. Add an optional mode where
every file is a **file tab** of one window. Default for everyone: Tabs. Switchable in one click for anyone who
prefers the old window-per-file behaviour.

## Words (go in `GLOSSARY.md`)

- **File tab**: the outer tab, one per open file.
- **View tab**: the inner tab (Map, Sequence...), today's `Workspace` tabs. "Re&set Tabs" becomes "Reset View Tabs".

## Decisions

| # | Decision |
|---|----------|
| Architecture | `MainWindow` becomes a **host**: menu bar, toolbar, status bar, file-tab strip. Each file is a **page** (its own editor, `Workspace`, Filters and History docks, search bar and text, range box, lock and circular state, autosave timer, background search tokens). Windows mode = a host with exactly one page, so there is one code path. |
| Isolation | Nothing stateful is shared between pages and nothing is rebound: the host flips a stacked widget. Menu actions and shortcuts go through `host.current_page`; `MainWindow.editor` becomes a property returning that page's editor (keeps the ~100 call sites). Every async callback (search `_found`, `run_busy`, timers, threads) captures its own page, never `self.editor`. Status bar is shared; a page's messages show only while it is active. |
| Setting | One global setting "Open files in Tabs" (checkable, **View** menu; also a line in Settings > General). Default on. Only affects future opens; never moves files already open. |
| Layout (same in both modes) | Row 1: Lock, Circular, then the file tabs (stretching, scrolling), then search and range on the right. Row 2: view tabs on the left; right corner: code switch, reference picker, Filters, History, FASTA sequence picker, quick filters. Overflow `»` as today; Lock and Circular stay in row 1. In Windows mode the strip is hidden. |
| What opens in a tab | Every in-app open: File > Open, drag-drop, recent, New, NCBI import, Assembly and PCR products. A file already open: its tab is selected, no second copy. An empty window opens the file in place. |
| Closing, keys | Close tab or window asks per modified file (reuse `_confirm_discard`; Cancel stops everything). Closing the last tab leaves the empty screen. Ctrl+W, Ctrl+Tab / Ctrl+Shift+Tab, Ctrl+1..9. Modified marker on the tab. Right-click: Close, Close Others, Copy Path, Reveal in Folder, Move to New Window. |
| Tab dragging | Reorder inside the strip only. No tabs shared between windows. |
| Move to New Window | Reopens the file by path in a new window, then closes the tab. Modified: asks "Save first?". Untitled (NCBI, products): Save As first. Cancel anywhere leaves everything untouched. Old tab closes only after the new window has loaded. |
| Autosave | Moves from beside the file to `<data dir>/autosave/<name>-<8 hex of full-path hash>[.<record>].autosave.gb`. The prompt still names the original file. An old `*.autosave.gb` found beside a file is offered once, then deleted. Copies of files that no longer exist and are older than 30 days are swept. |
| Single instance | Always one process. A second launch hands all its `argv` files to the running instance (local socket) and exits. Setting = Tabs: tabs in the last-used window. Setting = Windows: new windows of the same process. |
| Assembly and other cross-file features | List every open file in every host window and tab by file name, read-only, as today. |
| Restart | Only the setting is remembered. No restore of open files or tabs in v1 (autosave recovery covers unsaved work). |

## Delivery: four steps, one PR, committed together at the end

Each step leaves Windows mode fully working and updates its own docs, `tests.json` entries and the rebuilt test report. Nothing is committed until all four are done and tested; then
`/simplify` and `/code-review` run once, before the single PR (CLAUDE.md). I report what changed and how to test after each step.

### Step 1: Page / host refactor (no visible change)
- Page object per file; host with one page; `editor` property; async-callback audit.
- Pure-logic class for tab bookkeeping (open files, active one, dedupe by path, close order, unsaved files).
- Tests: bookkeeping logic; two `Editor`/`Document` objects, edit one, the other is untouched. No widget tests.
- `GLOSSARY.md`; "Reset View Tabs" rename.
- **You test (everything must behave exactly as before):** open a file; edit, undo, redo; save, Save As; search (typing, next, previous, results while editing); range box; lock; circular; Filters and History docks open/close/size; view tabs add, rename, duplicate, split, Reset View Tabs; FASTA sequence picker; primers panel; Assembly panel with two windows; New sequence; NCBI import; close with unsaved changes; drop a file on a window; open a second file (still a new window in this PR).

### Step 2: File-tab strip, two-row layout, View setting, Tabs by default
- Strip, scrolling, modified marker, tooltips with full path, context menu, keys, reorder, dedupe.
- Toolbar rearranged as in the layout row. View > "Open Files in Tabs".
- Window title = active file.
- **You test:** open 3 files; switch tabs; type different search text in each, each keeps its own; edit in one, the others stay unmodified; lock one, others stay editable; docks differ per file; close modified tab (Save / Discard / Cancel); close window with 2 modified tabs; close last tab; open an already open file; Ctrl+W, Ctrl+Tab, Ctrl+1..9; Move to New Window (clean, modified, untitled); toggle View setting then open (tab vs window); Assembly lists all tabs; narrow window overflow; many tabs scroll; start a search and switch tab before it finishes.

### Step 3: Autosave to the data folder
- New path scheme, one-time migration of old copies, sweep, prompt text; `.gitignore`, USER_GUIDE, README, PRIVACY, FEEDBACK, TESTING docs; `test_documents.py` entries.
- **You test:** edit and wait 60 s: copy appears in `%APPDATA%\openDNA\autosave`, none beside the file; kill the program, reopen: recovery offered; save: copy gone; two same-named files in different folders do not mix; an old `.autosave.gb` beside a file is offered once then removed.

### Step 4: Single instance
- Local server and client hand-off, all `argv` files, honours the setting.
- **You test:** program running, double-click a file in Explorer: opens as a tab (Tabs) or new window (Windows), no second process in Task Manager; double-click several files at once; program not running: starts normally; after a crash a stale socket does not block startup.

## Risks

- PR 1 touches most of `main.py` and nothing exercises the UI automatically: the manual checklist above is the safety net.
- Tabs as default changes behaviour for every existing user: say so in the release notes.
- Qt 6.8 vs 6.12 differ in object lifetime (see CLAUDE.md): check the pages' deletion on tab close on the lower version too.
