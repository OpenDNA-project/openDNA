# `tests/test_feedback.py`

[← all test files](README.md) · source: [`tests/test_feedback.py`](../../tests/test_feedback.py)

## `test_scrub_removes_private_text`

**Result:** ✅ passed in 0.00 s · [line 84](../../tests/test_feedback.py#L84) · **Area:** Error reports: scrubbing private text · **Checked against:** hand-computed · **Strength:** 🟡 adequate

**Why it matters:** Without it a crash report or log could leak the user's name, folder paths, file names or their DNA/protein sequences into a public GitHub issue.

Checks feedback.scrub on strings containing the username, home paths, quoted and UNC paths, traceback frames, file names and DNA/protein sequences; ordinary text and URLs must pass through unchanged, and scrubbing twice gives the same result.

| # | Input | Expected output |
|---|---|---|
| 1 | 'Cannot write C:\Users\bobsmith\AppData\...\x.tsv for bobsmith' | No 'bobsmith', 'Users' or 'AppData' |
| 2 | 'failed: /home/bobsmith/data/run1'; 'owner BobSmith' | Home path and the name removed, case-insensitively |
| 3 | FileNotFoundError with a quoted Windows path containing 'Lab Notes' | No 'Lab Notes' or 'secret'; '<path>' present |
| 4 | '\\fileserver\lab\project\plasmid.fa not found' | No 'Lab' |
| 5 | Traceback line 'File "C:\Users\bobsmith\...\y.py", line 12, in run' | 'File "y.py", line 12, in run' kept; name removed |
| 6 | Text with my_secret_plasmid.gb, reads.FASTA, a 38-base DNA and a 40-letter protein | Exactly 2 '<sequence>' and 2 '<file>'; none of the originals left |
| 7 | GenBank-style '1 acgtacgtac gtacgtacgt acgtacgtac gt' | No 'acgt' left |
| 8 | Plain English sentence with a GitHub URL; '' and None | Unchanged; '' and '' |
| 9 | scrub(scrub(path + DNA + name)) | Equal to a single scrub |

**Review notes:**

- Mostly 'not in' checks: a scrub that removed too much would pass most rows (only the plain-sentence row, the traceback row and the placeholder counts guard against over-scrubbing).
- The /home/bobsmith input is scrubbed and checked twice in a row with overlapping assertions (lines 90-91).

## `test_report_links_are_encoded_capped_and_carry_nothing_private`

**Result:** ✅ passed in 0.12 s · [line 120](../../tests/test_feedback.py#L120) · **Area:** Error reports: GitHub issue links · **Checked against:** hand-computed · **Strength:** 🟢 strong

**Why it matters:** Without it the 'Report a bug' / 'Request a feature' links could be malformed, too long for the browser, lose the end of the log, or carry private data to GitHub.

Builds bug and feature issue URLs and decodes their query: template, title prefix, version and system info, encoding of special characters, the length cap keeping the end of the log, the last 30 log lines used when no details are given, injectable openers, and that no sequence, path or name ever appears.

| # | Input | Expected output |
|---|---|---|
| 1 | report_url('bug', 'Crash when opening a file') | Starts with ISSUES_URL/new?template=bug_report.yml&; title '[Bug]: Crash when opening a file'; version = project VERSION; os field mentions Python and Qt; no spaces or newlines in URL |
| 2 | Title 'a&b #c?' and details with &, =, #, newline, e-acute, %, <b> | Both decode back exactly |
| 3 | 2000-line details ending 'ValueError: the real problem', 500-char summary | URL <= URL_CAP; logs end with the last line and say 'earlier lines left out'; title <= prefix + SUMMARY_MAX |
| 4 | feature_url('Dark mode for the map'); report_url('feature','x'); report_url('nonsense') | feature_request.yml template, '[Feature]: ...' title; same as feature_url('x'); ValueError |
| 5 | open_report / open_feature / open_log_folder with a recording opener | Return the opened URL; log folder is data_dir/logs, exists, opened as a file: URL |
| 6 | Logging active, one warning 'marker line 42', report_url('bug') with no details | Logs field contains the marker and 'Start: openDNA' |
| 7 | 100 warnings entry-000..entry-099 | entry-070..099 present, entry-069 absent (last 30 lines) |
| 8 | Summary, details, logged error and feature text containing a Windows path, DNA and the username | None of DNA, 'bobsmith', 'secret', 'Lab Notes', 'Users', 'Documents' in any decoded field |

**Review notes:**

- The privacy section partly re-tests scrub (previous test) through the link builder; reasonable as an end-to-end guard.

## `test_log_file_and_error_hooks`

**Result:** ✅ passed in 0.10 s · [line 191](../../tests/test_feedback.py#L191) · **Area:** Error reports: log file and crash dialog · **Checked against:** hand-computed · **Strength:** 🟢 strong

**Why it matters:** Without it unhandled errors could go unlogged or unreported, the same error could pop a dialog repeatedly, the log could grow without bound or contain private data, or the error dialog itself could crash the app.

Installs the logging and exception hooks with fake dialog, opener and clipboard and checks the scrubbed log file, Qt warning capture, no handler stacking, size rotation, diagnostics text, one dialog per distinct error (main and worker threads) with scrubbed text, the report and copy buttons, a failing dialog and KeyboardInterrupt handled, and hooks restored on teardown.

| # | Input | Expected output |
|---|---|---|
| 1 | Log one INFO line containing a private path | data_dir/logs/opendna.log has 'Start: openDNA <version>', Python, Qt, the message, INFO and the logger name; no 'bobsmith' or 'secret'; no dialog |
| 2 | qWarning('QWidget: something odd happened in plasmid.gb') | Log has WARNING, 'qt:' and the message, but not 'plasmid.gb' |
| 3 | setup_logging called a second time inside an active session | Only one extra handler; back to the original count after teardown |
| 4 | LOG_MAX_BYTES=2000, 400 lines logged | Exactly opendna.log, .log.1, .log.2; each <= 2200 bytes |
| 5 | copy_diagnostics with a fake clipboard; diagnostics_text in an empty data folder | Copied text = returned text, contains name, version, 'OS:', the log line, no username; '(the log is empty)' |
| 6 | Same ValueError 3x, a KeyError, the ValueError with another message | 2 dialogs total; log has 5 'Unhandled exception' and full tracebacks |
| 7 | ValueError message containing DNA and a private path | Dialog summary starts 'ValueError:'; no private text in summary or details |
| 8 | Dialog answers 'report' / 'copy' | One URL opened with the bug template and the error in title and logs / one clipboard copy, no URL |
| 9 | Dialog raising RuntimeError; a KeyboardInterrupt | Hook does not raise / no dialog |
| 10 | Same ValueError raised in two worker threads | One dialog, summary 'ValueError: from a worker' |
| 11 | setup_logging(hooks=False); teardown twice | sys/threading excepthooks unchanged; double teardown harmless |

**Review notes:**

- Very long test bundling about ten independent features; splitting it would make a failure point at one feature.
- The assertion that 'needs+a+report' is not in the URL is an indirect check of the encoding style; its intent is unclear.
