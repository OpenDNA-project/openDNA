# Every external program runs through the Terminal

openDNA runs external command-line programs (mappers, assemblers, multiple aligners) only through
the Terminal tab: one persistent bash session per document (in WSL on Windows), fed through a pipe,
with pipelines being forms that type their commands into it. We chose this over the earlier
one-dialog-per-tool approach (`RunDialog`, Map Reads) so the user always sees, can copy and can edit
exactly what ran, and so adding a pipeline never means a new dialog. Everything that comes back into
the app goes through the session's `open` command.

## Consequences

- The session is a pipe, not a terminal emulator: full-screen programs (`vim`, `top`, `less`) do not
  work, and one command runs at a time per document.
- HTML reports (fastp, QUAST) open in the system browser, not inside the app: the app ships
  `PySide6-Essentials` without QtWebEngine (~200 MB) and a browser already shows them well. The
  Terminal itself stays plain text.
