# `tests/test_updater.py`

[← all test files](README.md) · source: [`tests/test_updater.py`](../../tests/test_updater.py)

## `test_release_check_offers_only_a_newer_complete_release`

**Result:** ✅ passed in 0.00 s · [line 125](../../tests/test_updater.py#L125) · **Area:** Updater: which release is offered · **Checked against:** hand-computed · **Strength:** 🟢 strong

**Why it matters:** Without it the app could offer users an older, draft, prerelease or incomplete release, one without the file their system installs, or one whose files are hosted off the project's GitHub repository, as an update.

Checks version parsing and numeric comparison, then feeds updater.check fake GitHub 'releases/latest' JSON and checks that only a newer, final release with both the file this system installs (installer, .app zip or AppImage) and its .sha256 on the project's own repo is offered; every other case returns no update with an error text and never raises.

| # | Input | Expected output |
|---|---|---|
| 1 | is_newer on 15 pairs, e.g. v0.10.0 vs 0.9.0, 1.0 vs 1.0.0, 1.0.0.1 vs 1.0, 0.2.0-beta1 vs 0.1.0, 'latest', '', None, current 'dev' | True only for genuinely higher numeric versions (0.10.0 > 0.9.0; 1.0 == 1.0.0); prerelease, non-version or None -> False |
| 2 | parse_version('v0.10.0'), ' V2.1 ', '1.2.3-rc1', '1.2.3+build' | (0,10,0), (2,1), None, None |
| 3 | release_asset('v0.7.0', 'win32' / 'linux' / 'darwin' with machine arm64 and x86_64); a Linux copy checking a release that only has the Windows installer | oDNA_beta_0_7_0-setup.exe / -linux-x86_64.AppImage / -macos-arm64.zip and -macos-x86_64.zip; no update for the Linux copy |
| 4 | the names release_asset gives for win32, darwin (arm64) and linux, against the text of .github/workflows/release.yml | -setup.exe, -macos-arm64.zip and -linux-x86_64.AppImage; each suffix (its stem and its extension) appears in release.yml |
| 5 | check('0.1.0') with a v0.2.0 release carrying exe + .sha256 + other.zip | Update with version 0.2.0, tag v0.2.0, notes 'Fixes things.', exact asset, sha256 and page URLs; check_for_update returns the same object |
| 6 | check() with no version given, release tag = current project version | No update; exactly one request, to LATEST_RELEASE_API |
| 7 | check_for_update('0.1.0') with module-level fetch monkeypatched | Version 0.2.0 (default fetch is used) |
| 8 | Releases tagged v0.1.0 (same), v0.0.9 (older), v0.3.0 draft, v0.3.0 prerelease, v0.3.0-beta1 | No update; only the -beta1 tag produces a 'not understood' error |
| 9 | Release tagged 'nightly' | No update; error mentions 'nightly' |
| 10 | Release with no assets / only the exe / only the checksum | No update; error names the missing asset; page_url still set |
| 11 | Exe URL on evil.example.com, http://github.com, or another GitHub repo | No update; error contains 'Refused' |
| 12 | html_url pointing to evil.example.com | Update offered but page_url replaced by project RELEASES_URL |
| 13 | fetch raising UpdateError('...offline') / RuntimeError('boom') | Error ends with 'offline' / contains 'failed'; no exception |
| 14 | Bodies 'not json', '[1, 2]', '', b'\xff\xfe\x00' | No update; error contains 'understood' |
| 15 | check_for_update with fetch raising UpdateError | Returns None |

**Review notes:**

- Large multi-section test; a failure in an early section hides the later ones.
- Error texts are checked by substring ('Refused', 'understood', 'failed'), not exact message.
- The release.yml check is by substring (stem and extension separately), so it catches a renamed suffix but not every possible mismatch.

**Fixed by the audit of November 2026:**

- 13 Nov 2026: gap closed. The proposed checks were added, pass on the current code and fail under the mutation that slipped through before. Was: Claim 'both the exe and its .sha256 on the project's own repo': only the EXE URL is ever moved off-host in the test; the checksum URL is never tested, so check() may offer a release whose checksum comes from any host/scheme (code also does not enforce the own-repo path for the .sha256, so the json claim overstates the code).

## `test_download_comes_only_from_github_and_is_kept_only_when_it_matches_its_checksum`

**Result:** ✅ passed in 0.03 s · [line 223](../../tests/test_updater.py#L223) · **Area:** Updater: download safety · **Checked against:** hand-computed · **Strength:** 🟢 strong

**Why it matters:** Without it the updater could download a program from a look-alike or hijacked host, follow a redirect off GitHub, or install a tampered or partial file.

Checks the HTTPS GitHub host allow-list (direct and on redirects), the real fetch refusing bad URLs before any network call, the checksum-file parser, and download_update: verified file returned, mismatch / missing / unreadable checksum and off-GitHub URLs refused, and no partial file left after cancel or failure.

| # | Input | Expected output |
|---|---|---|
| 1 | check_url on api.github.com, github.com, objects/release-assets.githubusercontent.com, 'GitHub.com' | Accepted, URL returned unchanged |
| 2 | check_url on http://, ftp://, evil.example.com/github.com, github.com.evil.example.com, evilgithub.com, user:pw@github.com, github.com@evil.example.com, github.com:8443, gist.github.com, '', None, file:/// | UpdateError for each |
| 3 | Real updater.fetch('http://github.com/a/b') and ('https://evil.example.com/a') with build_opener patched to fail | UpdateError raised before any connection attempt |
| 4 | _SafeRedirects.redirect_request to evil.example.com / to objects.githubusercontent.com | UpdateError / a new request object |
| 5 | parse_checksum on '<hex>  openDNA.exe', '<hex> *openDNA.exe', bare hex, upper-case hex, multi-line list, UTF-16 with BOM, UTF-8 BOM | The lower-case digest each time |
| 6 | parse_checksum on '', 'hello', 'abc  openDNA.exe', digest for other.zip, 63 or 65 hex chars | None |
| 7 | download_update with matching exe bytes and a progress callback | Returns folder/.openDNA-update-openDNA.exe.new with the exact bytes; no .part left; last progress = (total, total) |
| 8 | Exe bytes + b'tampered' | UpdateError 'does not match'; folder empty |
| 9 | Update with sha256_url='' / checksum body '<html>404</html>' | UpdateError 'no checksum' (exe never requested) / 'expected form'; folder empty |
| 10 | asset_url on evil.example.com / sha256_url on http://github.com | UpdateError 'Refused'; no request made |
| 11 | cancel event already set | UpdateError 'cancelled'; folder empty |
| 12 | fetch that writes 'half' then raises 'Could not reach' | That UpdateError propagates; folder empty |
| 13 | Target folder that does not exist | UpdateError 'Could not write' |

**Review notes:**

- The real network path is not exercised (correct for a unit test); the redirect check calls the handler method directly rather than going through a real redirect.

**Fixed by the audit of November 2026:**

- 13 Nov 2026: gap closed. The proposed checks were added, pass on the current code and fail under the mutation that slipped through before. Was: 'enforced by the real fetch ... on redirects': the test calls _SafeRedirects.redirect_request directly; nothing shows fetch() installs that handler, so a fetch that follows a redirect to any host would pass.

## `test_install_swaps_relaunches_and_cleans_up`

**Result:** ✅ passed in 0.03 s · [line 337](../../tests/test_updater.py#L337) · **Area:** Updater: install and cleanup · **Checked against:** hand-computed · **Strength:** 🟢 strong

**Why it matters:** Without it a successful update could fail to replace the installed program (installer not started, .app or AppImage not swapped), fail to restart it, start it with a broken PyInstaller environment, or leave downloads and .old files behind (or delete the user's own files).

On fake files in temp folders, checks that _swap (the step install_update ends with on macOS and Linux) swaps new over current keeping the old as .old (a file, an AppImage named by $APPIMAGE, a macOS .app folder, a zip unpacked with ditto where ditto exists), that on Windows it runs the downloaded installer silently with /RELAUNCH=1 instead of swapping, that relaunch clears PyInstaller environment variables, and that cleanup_old removes only leftovers and downloads, does nothing from source and retries a held file.

| # | Input | Expected output |
|---|---|---|
| 1 | current='OLD', new='NEW', install_update with a recording relaunch | current now 'NEW', openDNA.exe.old = 'OLD', .new gone, relaunch called once with current |
| 2 | Same with an existing .old = 'ANCIENT' | .old now holds 'OLD' |
| 3 | frozen, platform linux, $APPIMAGE = current, no target given | current replaced with 'NEW' |
| 4 | macOS .app folders: openDNA.app 'OLD', openDNA.app.new 'NEW', stale openDNA.app.old; install_update, then cleanup_old(app) | app holds 'NEW', .new gone, relaunch called with the app, .old holds 'OLD'; cleanup removes the .old folder and keeps the app |
| 5 | (only where ditto exists, i.e. on a Mac) a zip made with ditto holding an openDNA.app 'ZIPPED', frozen darwin from that .app | the app now holds 'ZIPPED' and the downloaded zip is gone |
| 6 | frozen win32 beside unins000.exe, installer downloaded as .openDNA-update-oDNA_beta_0_2_0-setup.exe.new in the data folder (spawn faked) | the file is renamed to .exe and started with /VERYSILENT /SUPPRESSMSGBOXES /NORESTART /CLOSEAPPLICATIONS /RELAUNCH=1; the installed exe is untouched; relaunch is not used |
| 7 | then an old .openDNA-update-...-setup.exe.part in the data folder, the user's notes.txt and a browser's oDNA_beta_0_9_0-setup.exe.part beside the exe; cleanup_old() | True; the data folder is empty (installer and .part removed); notes.txt and the browser download are still there |
| 8 | relaunch() with _MEIPASS2 and _PYI_APPLICATION_HOME_DIR set (Popen faked) | args=[exe path]; env has PYINSTALLER_RESET_ENVIRONMENT=1 and no _MEIPASS2 or _PYI* keys |
| 9 | cleanup_old(current) with .old and .part beside the exe; then again | True; only openDNA.exe remains with 'OLD'; second call True |
| 10 | cleanup_old() with no argument, not frozen | True and the .old file is untouched |
| 11 | Path.unlink raising PermissionError twice then succeeding, wait=60 | True and the .old file is deleted |

## `test_install_refuses_or_rolls_back_leaving_the_old_program`

**Result:** ✅ passed in 0.02 s · [line 449](../../tests/test_updater.py#L449) · **Area:** Updater: refusal and rollback · **Checked against:** hand-computed · **Strength:** 🟢 strong

**Why it matters:** Without it a failed update (no write access, antivirus lock, broken new program, a portable copy, a copy that is not installed) could leave the user with no working program, or a portable copy could replace itself.

Checks that install refuses when run from source or from a portable copy (not set up by the installer: pointing to the releases page), that installed_copy recognises an installed copy on each system (Windows: unins000.exe beside the exe; macOS: a .app; Linux: an existing $APPIMAGE) and that a folder it may not write to blocks the update, that a missing file or a file in another folder is refused, and that a failing swap or a new program that does not start rolls back to the old program; a briefly held file is retried until it succeeds.

| # | Input | Expected output |
|---|---|---|
| 1 | Not frozen: install_update(new) | UpdateError 'source code'; nothing started; new untouched; install_blocker() mentions 'source code', install_blocker(Path('x.exe')) is None |
| 2 | Frozen, no unins000.exe / $APPIMAGE / .app folder, on win32, darwin and linux; install_update(new) | installed_copy() None; install_blocker says 'portable' and names RELEASES_URL; UpdateError 'portable'; files unchanged |
| 3 | Frozen win32 with unins000.exe beside the exe; darwin with the exe inside openDNA.app/Contents/MacOS; linux with $APPIMAGE = an existing file | installed_copy is the exe / the .app folder / the AppImage file; install_blocker None |
| 4 | Linux with os.access saying the folder is not writable; $APPIMAGE pointing to a file that is gone | blocker says 'no write access'; installed_copy None |
| 5 | download_folder() when run from source | the data folder's 'update' subfolder |
| 6 | os.replace failing when moving the new file in | UpdateError 'administrator'; current 'OLD', new 'NEW', no .old left |
| 7 | os.replace always failing when moving current out | UpdateError naming RELEASES_URL and 'administrator'; exactly RETRIES attempts; files unchanged |
| 8 | os.replace failing twice then succeeding | Install succeeds after 2 refusals; current 'NEW' |
| 9 | relaunch raising OSError('not a valid Win32 application') | UpdateError 'old one was kept'; current 'OLD', new 'NEW', no .old |

## `test_automatic_check_obeys_the_switches_runs_once_a_day_and_remembers_a_skip`

**Result:** ✅ passed in 0.17 s · [line 548](../../tests/test_updater.py#L548) · **Area:** Updater: automatic check schedule · **Checked against:** hand-computed · **Strength:** 🟡 adequate

**Why it matters:** Without it the app could check for updates when the user turned it off, check on every start, stop checking after the clock goes back, or forget a version the user chose to skip.

Checks that the automatic check follows the user setting and the project-wide UPDATES_ENABLED flag, is due at most once per 24 h (and again if the clock went back), and that skip_version is stored.

| # | Input | Expected output |
|---|---|---|
| 1 | Fresh state | check_enabled() and auto_check_due() True |
| 2 | Setting Updates/check=False, then True | Not due / due again |
| 3 | project.UPDATES_ENABLED=False | Not due |
| 4 | mark_checked(t), then ask at t+60 s, t+23 h, t+24 h, t-1 h | Not due, not due, due, due |
| 5 | skipped_version() fresh, then skip_version('0.2.0') | '' then '0.2.0' |

**Review notes:**

- Only checks that a skipped version is stored, not that the automatic check then stays quiet about it (that logic, if in main.py, is untested here).

**Fixed by the audit of November 2026:**

- 13 Nov 2026: gap closed. The proposed checks were added, pass on the current code and fail under the mutation that slipped through before. Was: 'remembers a skip' (role: 'forget a version the user chose to skip'): only the store round trip is tested; the one place the skip has an effect, startup_check's `result.update.version != skipped_version()`, is untested, so a skipped version would be offered again every day. \| startup_check's own obedience to auto_check_due and its mark_checked call are untested too.
