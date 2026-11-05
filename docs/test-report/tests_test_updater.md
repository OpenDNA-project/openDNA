# `tests/test_updater.py`

[← all test files](README.md) · source: [`tests/test_updater.py`](../../tests/test_updater.py)

## `test_versions`

**Result:** ✅ passed in 0.00 s · [line 124](../../tests/test_updater.py#L124) · **Area:** Self-update: version comparison · **Checked against:** hand-computed · **Strength:** 🟢 strong

Makes sure the updater only offers a release whose version is really higher than the running one: numbers are compared one by one, and pre-releases or junk tags are never offered.

| # | Input | Expected output |
|---|---|---|
| 1 | is_newer pairs: ('v0.10.0','0.9.0'), ('0.2.0','0.1.0'), ('v1.0.0','0.99.99'), ('0.1.1','0.1.0'), ('1.0.0.1','1.0') | True for each |
| 2 | is_newer pairs: ('1.0','1.0.0'), ('0.1.0','0.1.0'), ('v0.1.0','0.1.0'), ('0.9.0','0.10.0'), ('0.0.9','0.1.0') | False for each (equal or older; 1.0 == 1.0.0) |
| 3 | is_newer with '0.2.0-beta1', 'latest', '', None as candidate, or current 'dev' | False for each |
| 4 | parse_version('v0.10.0'), parse_version(' V2.1 ') | (0, 10, 0) and (2, 1) |
| 5 | parse_version('1.2.3-rc1'), parse_version('1.2.3+build') | None for both |

## `test_release_json`

**Result:** ✅ passed in 0.00 s · [line 142](../../tests/test_updater.py#L142) · **Area:** Self-update: reading the GitHub release · **Checked against:** hand-computed · **Strength:** 🟢 strong

Checks that a newer GitHub release is turned into an Update with the right download, checksum and page addresses, and that up-to-date, older, draft, prerelease and badly named releases offer nothing.

| # | Input | Expected output |
|---|---|---|
| 1 | check('0.1.0') on a fake release JSON with tag 'v0.2.0', body 'Fixes things.', exe + .sha256 + other.zip assets | error None; version '0.2.0', tag 'v0.2.0', notes 'Fixes things.', asset_url .../v0.2.0/OpenGenEdit.exe, sha256_url ...exe.sha256, page_url .../releases/tag/v0.2.0; check_for_update returns the same Update |
| 2 | check(fetch=fake) with no version, release tag = v{project.VERSION} | update None; the only request is project.LATEST_RELEASE_API |
| 3 | updater.fetch replaced by the fake, check_for_update('0.1.0') with no fetch argument | version '0.2.0' (the module-level fetch is used) |
| 4 | current 0.1.0; release v0.1.0, v0.0.9, v0.3.0 draft, v0.3.0 prerelease, v0.3.0-beta1 | update None each; error None or containing 'not understood' |
| 5 | release tag 'nightly' | update None and 'nightly' in the error |

**Fixed by the audit of November 2026:**

- The draft/prerelease/older loop now asserts that only the 'v0.3.0-beta1' tag gives the 'not understood' message and the others give no error.

## `test_release_assets_and_errors`

**Result:** ✅ passed in 0.00 s · [line 174](../../tests/test_updater.py#L174) · **Area:** Self-update: release assets and failures · **Checked against:** hand-computed · **Strength:** 🟢 strong

Protects the user from installing an incomplete or foreign release: a release missing the exe or its checksum, or pointing to another host or repository, is refused, and network errors or garbage answers never crash the check.

| # | Input | Expected output |
|---|---|---|
| 1 | newer release with no assets / only the exe / only the .sha256 | update None; error names the missing OpenGenEdit.exe.sha256 or OpenGenEdit.exe; page_url set |
| 2 | exe asset URL 'https://evil.example.com/OpenGenEdit.exe', 'http://github.com/x/OpenGenEdit.exe', 'https://github.com/someone/else/releases/download/v9/OpenGenEdit.exe' | update None and 'Refused' in error |
| 3 | html_url='https://evil.example.com/x' | update offered with page_url == project.RELEASES_URL |
| 4 | fetch raising UpdateError('Could not reach api.github.com: offline') / RuntimeError('boom') | error ends with 'offline' / contains 'failed'; no exception |
| 5 | answers 'not json', '[1, 2]', '', b'\xff\xfe\x00' | update None and 'understood' in error |
| 6 | check_for_update with a fetch raising UpdateError('x') | None |

## `test_allowed_addresses`

**Result:** ✅ passed in 0.00 s · [line 210](../../tests/test_updater.py#L210) · **Area:** Self-update: allowed download hosts · **Checked against:** hand-computed · **Strength:** 🟢 strong

Ensures the updater only ever talks to HTTPS GitHub hosts, refuses look-alike or tricked addresses before any network access, and refuses a redirect to another host.

| # | Input | Expected output |
|---|---|---|
| 1 | check_url on https api.github.com, github.com, objects.githubusercontent.com, release-assets.githubusercontent.com, 'https://GitHub.com/a' | returns the URL unchanged |
| 2 | check_url on http://github.com, ftp://, https://evil.example.com/github.com, github.com.evil.example.com, evilgithub.com, user:pw@github.com, github.com@evil.example.com, github.com:8443, gist.github.com, '', None, file:///C:/x.exe | raises UpdateError |
| 3 | real updater.fetch('http://github.com/a/b') and fetch('https://evil.example.com/a') with build_opener patched to fail if called | raises UpdateError (no connection attempted) |
| 4 | _SafeRedirects.redirect_request 302 from a github.com asset to evil.example.com / to objects.githubusercontent.com | raises UpdateError / returns a request (not None) |

## `test_checksum_and_download`

**Result:** ✅ passed in 0.00 s · [line 248](../../tests/test_updater.py#L248) · **Area:** Self-update: checksum and download · **Checked against:** hand-computed · **Strength:** 🟢 strong

Checks that the downloaded program is verified against its published SHA-256 before use, and that any failure (mismatch, bad checksum file, foreign address, cancel, network error, unwritable folder) leaves no partial file behind.

| # | Input | Expected output |
|---|---|---|
| 1 | parse_checksum on '<digest>  OpenGenEdit.exe', '<digest> *OpenGenEdit.exe', bare digest, upper-case digest, a two-line file (other.zip first), UTF-16 with BOM, UTF-8 BOM | the lower-case digest each time |
| 2 | parse_checksum on '', 'hello', 'abc  OpenGenEdit.exe', '<digest>  other.zip', 63-char and 65-char hex | None each |
| 3 | download_update with the fake serving EXE_BYTES ('MZ' + 'fake exe '*1000) and the matching digest | returns folder/OpenGenEdit.exe.new with EXE_BYTES; no .part left; last progress = (len, len) |
| 4 | exe bytes + b'tampered' | UpdateError 'does not match'; folder empty |
| 5 | sha256_url='' / checksum text '<html>404</html>' | UpdateError 'no checksum' with no .exe request / 'expected form'; folder empty |
| 6 | asset_url on evil.example.com / sha256_url on http://github.com | UpdateError 'Refused'; no request made |
| 7 | cancel Event already set / a fetch that writes b'half' then raises 'Could not reach' | UpdateError 'cancelled' / 'Could not reach'; folder empty both times |
| 8 | folder tmp_path/'nowhere' (does not exist) | UpdateError 'Could not write' |

**Review notes:**

- Lines 294-301: the cancel is honoured by the test's own Fake (line 65-66), not by updater.fetch; the real fetch's cancel check, size caps (MAX_TEXT/MAX_EXE) and 'stopped early' Content-Length check are not tested anywhere in the suite.

## `test_install_swaps_and_refuses`

**Result:** ✅ passed in 0.01 s · [line 341](../../tests/test_updater.py#L341) · **Area:** Self-update: installing the new exe · **Checked against:** hand-computed · **Strength:** 🟢 strong

Checks the file swap that installs an update (running exe to .old, new exe in place, relaunch) and that it refuses to act when running from source, as the folder build, or with missing / misplaced files.

| # | Input | Expected output |
|---|---|---|
| 1 | folder with OpenGenEdit.exe=b'OLD' and .exe.new=b'NEW'; install_update(new, current, relaunch=recorder) | exe = b'NEW', .exe.old = b'OLD', .new gone, relaunch called with [current] |
| 2 | same with a stale .exe.old = b'ANCIENT' | .exe.old now b'OLD' |
| 3 | not frozen, install_update(new) without current_exe | UpdateError 'source code', relaunch never called, .new kept; install_blocker() mentions 'source code', install_blocker(Path('x.exe')) is None |
| 4 | sys.frozen=True, _MEIPASS a temp _MEI1234 folder, then _MEIPASS = exe folder/_internal and exe folder | one-file: not folder build, blocker None; folder layouts: is_folder_build True, blocker mentions 'OpenGenEdit-folder.zip', install raises 'folder version'; files unchanged |
| 5 | sys.frozen=True, sys.executable=current; install_update(new) | current now b'NEW' |
| 6 | new exe path that does not exist / new exe in another folder | UpdateError 'could not be found' / 'same folder'; current still b'OLD' |

## `test_install_rolls_back`

**Result:** ✅ passed in 0.00 s · [line 400](../../tests/test_updater.py#L400) · **Area:** Self-update: install failure and rollback · **Checked against:** hand-computed · **Strength:** 🟢 strong

Ensures a failed install never leaves the user without a working program: each failure puts the old exe back, a read-only folder points to the releases page, and transient antivirus locks are retried.

| # | Input | Expected output |
|---|---|---|
| 1 | os.replace fails (PermissionError) when moving the new file in | UpdateError 'administrator'; exe b'OLD', .new b'NEW', no .old |
| 2 | os.replace fails for every rename of the running exe | UpdateError containing project.RELEASES_URL and 'administrator'; exactly RETRIES (5) attempts; files unchanged |
| 3 | rename of the running exe fails twice then succeeds | install succeeds after 2 refusals; exe b'NEW' |
| 4 | relaunch raises OSError('not a valid Win32 application') | UpdateError 'old one was kept'; exe b'OLD', .new b'NEW', no .old |

## `test_relaunch_and_cleanup`

**Result:** ✅ passed in 0.00 s · [line 445](../../tests/test_updater.py#L445) · **Area:** Self-update: relaunch and cleanup · **Checked against:** hand-computed · **Strength:** 🟢 strong

Checks that the restarted program gets a clean PyInstaller environment and that leftovers of an update are deleted at start-up (waiting for a file still held), without touching anything when run from source.

| # | Input | Expected output |
|---|---|---|
| 1 | relaunch(tmp/OpenGenEdit.exe) with Popen faked, env _MEIPASS2 and _PYI_APPLICATION_HOME_DIR set | args [exe path]; env PYINSTALLER_RESET_ENVIRONMENT='1'; no _MEIPASS2 and no _PYI* keys |
| 2 | folder with OpenGenEdit.exe, .exe.old, .exe.part; cleanup_old(current) twice | True; only OpenGenEdit.exe remains (still b'OLD'); second call True |
| 3 | cleanup_old() with no exe, not frozen | True and .old untouched |
| 4 | Path.unlink raising PermissionError for the first 2 calls, sleep faked; cleanup_old(current, wait=60) | True and .old deleted |

**Review notes:**

- Line 460: only .old and .part leftovers are created; the .new leftover that cleanup_old also deletes is never tested. Popen's cwd/creationflags are not checked.

## `test_remembered_state`

**Result:** ✅ passed in 0.01 s · [line 488](../../tests/test_updater.py#L488) · **Area:** Self-update: settings and throttle · **Checked against:** hand-computed · **Strength:** 🟢 strong

Checks that the automatic update check obeys the user setting and the project switch, runs at most once every 24 hours (even if the clock goes back), and that a skipped version is remembered.

| # | Input | Expected output |
|---|---|---|
| 1 | default state; setting 'check' False; back True; project.UPDATES_ENABLED False | due/enabled True; both False; due True; due False |
| 2 | mark_checked(1e6); auto_check_due at +60 s, +23 h, +24 h, -1 h | False, False, True, True |
| 3 | skipped_version() before and after skip_version('0.2.0') | '' then '0.2.0' |

## `test_startup_check`

**Result:** ✅ passed in 0.10 s · [line 519](../../tests/test_updater.py#L519) · **Area:** Self-update: start-up check · **Checked against:** GUI state · **Strength:** 🟡 adequate

Checks that at start-up a new version opens the update dialog once a day, while up-to-date, offline, skipped and disabled cases stay silent.

| # | Input | Expected output |
|---|---|---|
| 1 | VERSION '0.1.0'; startup_check(window, fake v0.2.0, now=1e6) | True; within 5 s window._update_dialog exists, shows version '0.2.0', visible |
| 2 | startup_check again at now=1e6+100 | False; fake saw only one request (the API) |
| 3 | release v0.1.0 at now=1e6 / fetch error 'offline' at now=2e6 | True (looked); after jobs finish no _update_dialog |
| 4 | skip_version('0.2.0') then startup_check with v0.2.0 at now=1e7 | no dialog, no message |
| 5 | UPDATES_ENABLED False | False and no request |

**Review notes:**

- Lines 525, 536, 541: relies on a real background thread and wait_for polling up to 5 s; could be flaky on a loaded CI machine (failure would show as a missing dialog).

## `test_menu_action`

**Result:** ✅ passed in 0.06 s · [line 553](../../tests/test_updater.py#L553) · **Area:** Self-update: Check for Updates menu · **Checked against:** GUI state · **Strength:** 🟡 adequate

Checks the Help-menu action: it says 'not available yet' without network while updates are disabled, reports 'up to date' and errors with the releases link, and opens the dialog even for a skipped version.

| # | Input | Expected output |
|---|---|---|
| 1 | UPDATES_ENABLED False; check_now(window) | False; no request; message contains 'not available yet' and RELEASES_URL |
| 2 | VERSION 0.1.0; release v0.1.0 | True; message contains 'up to date' |
| 3 | fetch error 'Could not reach api.github.com: offline' | second message contains 'offline' and RELEASES_URL |
| 4 | skip_version('0.2.0'); check_now with v0.2.0 release | window._update_dialog appears |

**Review notes:**

- Line 572: only checks a dialog appears, not that it is for version 0.2.0. Uses 5 s thread waits (lines 567, 569, 572).

## `test_job_gc_handling`

**Result:** ✅ passed in 0.09 s · [line 588](../../tests/test_updater.py#L588) · **Area:** Self-update: background job / garbage collector · **Checked against:** property/invariant · **Strength:** 🟢 strong

Checks the crash workaround that pauses Python's cycle collector while a worker thread runs and restores it afterwards, without turning it on if it was off.

| # | Input | Expected output |
|---|---|---|
| 1 | two _Job workers blocked on one gate | gc disabled while they run; after the gate opens both callbacks arrive, gc enabled again, _gc_jobs == 0 |
| 2 | gc.disable() then one _Job returning 1 | callback arrives and gc is still disabled |

**Review notes:**

- Line 594 asserts gc.isenabled() at start: depends on no earlier test leaving gc disabled.

**Fixed by the audit of November 2026:**

- The two jobs now have their own gates: after the first ends the test asserts the collector is still off (one job left), and on only after the second.

## `test_dialog_install_flow`

**Result:** ✅ passed in 0.39 s · [line 620](../../tests/test_updater.py#L620) · **Area:** Self-update: update dialog · **Checked against:** GUI state · **Strength:** 🟢 strong

Checks the update dialog end to end with fakes: it shows the notes, installs only a verified download and then quits, offers 'Try Again' after a checksum failure, lets the unsaved-work confirm stop the install, and remembers 'Skip' but not 'Later'.

| # | Input | Expected output |
|---|---|---|
| 1 | dialog for version 0.2.0 with notes 'Notes'; click Install | header contains '0.2.0', notes text 'Notes'; installed == [folder/OpenGenEdit.exe.new], quit called once, file == EXE_BYTES |
| 2 | fake exe b'MZ-tampered'; click Install | nothing installed, quit 0; status 'does not match'; button enabled, text 'Try Again'; folder empty |
| 3 | confirm=lambda: False; click Install | nothing installed; folder empty |
| 4 | click 'Skip This Version' | skipped_version() == '0.2.0', dialog hidden, nothing installed |
| 5 | click 'Later' | skipped_version() == '' and dialog hidden |

**Review notes:**

- Lines 623, 631, 640: background download waited with 5 s polling (timing-dependent).

## `test_dialog_cancel_and_from_source`

**Result:** ✅ passed in 0.07 s · [line 660](../../tests/test_updater.py#L660) · **Area:** Self-update: update dialog · **Checked against:** GUI state · **Strength:** 🟡 adequate

Checks that closing the dialog during a download cancels it and deletes the partial file, and that when running from source the dialog cannot install and points to the releases page.

| # | Input | Expected output |
|---|---|---|
| 1 | slow fake that writes b'part' then waits; click Install then Cancel | busy and Later reads 'Cancel' during download; after cancel nothing installed, quit 0, folder empty |
| 2 | UpdateDialog without injected install (running from source) | Install disabled; status contains 'source code' and RELEASES_URL; Skip and Later enabled |

**Review notes:**

- Line 666: the cancellation is raised by the Fake's own check (line 65), so this proves the Event is passed through, not that the real fetch stops streaming.
