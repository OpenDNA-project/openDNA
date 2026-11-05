# `tests/test_feedback.py`

[← all test files](README.md) · source: [`tests/test_feedback.py`](../../tests/test_feedback.py)

## `test_scrub_private_text`

**Result:** ✅ passed in 0.00 s · [line 87](../../tests/test_feedback.py#L87) · **Area:** Feedback: scrubbing · **Checked against:** hand-computed · **Strength:** 🟢 strong

feedback.scrub removes private data from text before it leaves the machine: usernames, home paths, quoted/UNC paths, traceback frame paths (file name kept), file names and DNA/protein sequences, while leaving ordinary text and URLs; it is idempotent.

| # | Input | Expected output |
|---|---|---|
| 1 | '...C:\Users\bobsmith\AppData\...\x.tsv for bobsmith'; '/home/bobsmith/data/run1'; 'owner BobSmith' | no 'bobsmith'/'Users'/'AppData'; username removed case-insensitively |
| 2 | FileNotFoundError with SECRET_PATH quoted; UNC '\\fileserver\lab\...' | 'Lab Notes'/'secret' gone, '<path>' present; 'Lab' gone |
| 3 | traceback frame 'File "C:\...\y.py", line 12, in run' | 'File "y.py", line 12, in run' kept, 'bobsmith' gone |
| 4 | file names my_secret_plasmid.gb, reads.FASTA; DNA; protein 'ACDEF...'x2; a spaced ACGT run | no 'secret'/'reads'/DNA/'ACDEFGHIK'; '<sequence>' x2, '<file>' x2; acgt run replaced |
| 5 | ordinary text with a github URL; '' ; None | returned unchanged; '' ; '' (None -> '') |
| 6 | scrub(scrub(SECRET_PATH+DNA+'bobsmith')) | equals the once-scrubbed text (idempotent) |

## `test_report_and_feature_urls`

**Result:** ✅ passed in 0.04 s · [line 125](../../tests/test_feedback.py#L125) · **Area:** Feedback: prefilled GitHub URLs · **Checked against:** hand-computed · **Strength:** 🟢 strong

report_url/feature_url build prefilled GitHub issue links: template, title, version and OS fields, URL-encoding of special characters, a length cap that keeps the end, feature requests and an unknown kind error, and injectable openers.

| # | Input | Expected output |
|---|---|---|
| 1 | report_url('bug','Crash when opening a file') | starts ISSUES_URL+'/new?template=bug_report.yml&'; title '[Bug]: ...'; version==VERSION; os has 'Python' and 'Qt'; no space/newline in url |
| 2 | details with & = # newline cafe 100% <b> | title '[Bug]: a&b #c?'; logs field == details (decoded) |
| 3 | 2000 long lines + 'ValueError: the real problem', summary 's'*500 | len(url) <= URL_CAP; logs end with the ValueError and contain 'earlier lines left out'; title <= 7 + SUMMARY_MAX |
| 4 | feature_url('Dark mode...'); report_url('feature','x'); report_url('nonsense') | feature template, '[Feature]:' title; report('feature',x)==feature_url(x); 'nonsense' raises ValueError |
| 5 | open_report/open_feature with opener list; open_log_folder | return value == recorded url; bug vs feature templates; folder == log_folder(), name 'logs', is_dir; opener got a 'file:' url |

## `test_report_url_from_the_log`

**Result:** ✅ passed in 0.01 s · [line 171](../../tests/test_feedback.py#L171) · **Area:** Feedback: URL from log · **Checked against:** property/invariant · **Strength:** 🟢 strong

With no explicit details the report link carries the last 30 log lines, and private data (sequences, paths, names) never reaches any link.

| # | Input | Expected output |
|---|---|---|
| 1 | log 'marker line 42', then report_url('bug') | logs field has 'marker line 42' and 'Start: openDNA' |
| 2 | 100 log entries, report_url('bug') | logs has entry-099 and entry-070, not entry-069 (last 30 only) |
| 3 | log + report_url/feature_url built with SECRET_PATH and DNA | decoded query has none of DNA,'bobsmith','secret','Lab Notes','Users','Documents' |

## `test_log_file`

**Result:** ✅ passed in 0.05 s · [line 200](../../tests/test_feedback.py#L200) · **Area:** Feedback: log file · **Checked against:** hand-computed · **Strength:** 🟢 strong

The log file is a scrubbed file in the data folder, captures Qt warnings, does not stack handlers when set up twice, rotates within a bounded size/count, and feeds the diagnostics text/copy.

| # | Input | Expected output |
|---|---|---|
| 1 | active(); log an INFO line mentioning SECRET_PATH; teardown; read file | path data/logs/opendna.log; has 'Start: openDNA <ver>', Python, Qt, the message, INFO, logger name; no 'bobsmith'/'secret'; no dialog shown |
| 2 | qWarning('... in plasmid.gb') | log has WARNING, 'qt:', the message; 'plasmid.gb' scrubbed out |
| 3 | setup_logging twice (hooks/qt off) | handler count n -> n+1 inside, back to n after teardown |
| 4 | LOG_MAX_BYTES=2000, write 400 lines | files == opendna.log, .1, .2; each <= 2200 bytes |
| 5 | copy_diagnostics with a fake clipboard; then fresh empty data dir | returned text == copied text, has name/version, 'OS:', the log line, no 'bobsmith'; empty folder -> '(the log is empty)' |

**Review notes:**

- L227-228 depends on setup_logging adding exactly one handler with hooks=False,qt_handler=False; it checks the count, not which handler.

## `test_excepthook_dialog_and_buttons`

**Result:** ✅ passed in 0.01 s · [line 254](../../tests/test_feedback.py#L254) · **Area:** Feedback: sys.excepthook · **Checked against:** hand-computed · **Strength:** 🟢 strong

The installed sys.excepthook shows one dialog per distinct error (every occurrence still logged), passes the dialog only scrubbed text, wires the report and copy buttons, never raises when the dialog itself fails, and ignores KeyboardInterrupt.

| # | Input | Expected output |
|---|---|---|
| 1 | boom (ValueError) x3, boom_other (KeyError), boom('other message') | dialogs 1, then 2, then still 2 (same type+place = same error); log has 5 'Unhandled exception', traceback, 'ValueError: x', 'KeyError' |
| 2 | boom with DNA+SECRET_PATH in the message | summary starts 'ValueError:'; neither summary nor details contain DNA,'bobsmith','secret','Lab Notes' |
| 3 | choice='report'; boom('needs a report') | one url opened, bug template, not raw 'needs+a+report' text in url but 'ValueError: needs a report' in logs and title |
| 4 | choice='copy'; boom('copy me') | one clipboard copy with 'ValueError: copy me'; no url |
| 5 | choice='boom' dialog raises; then KeyboardInterrupt | excepthook does not raise, one dialog; KeyboardInterrupt -> no dialog |

## `test_threading_hook_and_optional_hooks`

**Result:** ✅ passed in 0.01 s · [line 302](../../tests/test_feedback.py#L302) · **Area:** Feedback: threading hook · **Checked against:** hand-computed · **Strength:** 🟢 strong

threading.excepthook reports an error from another thread once (same error not re-shown), and the hooks are optional with idempotent teardown.

| # | Input | Expected output |
|---|---|---|
| 1 | worker thread raises ValueError('from a worker'), twice | one dialog, summary 'ValueError: from a worker'; second time still one dialog |
| 2 | setup_logging(hooks=False, qt_handler=False); teardown twice | sys/threading excepthook unchanged throughout |

**Review notes:**

- L308-316 the two worker threads have different names (worker-1 then default), so 'same error not re-shown' relies on the error key ignoring thread name; if it included the name the second dialog would appear - worth confirming the key is thread-independent by design.
