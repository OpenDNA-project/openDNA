# `tests/test_file_assoc.py`

[← all test files](README.md) · source: [`tests/test_file_assoc.py`](../../tests/test_file_assoc.py)

## `test_register`

**Result:** ✅ passed in 0.00 s · [line 107](../../tests/test_file_assoc.py#L107) · **Area:** File associations: register() · **Checked against:** hand-computed · **Strength:** 🟢 strong

register() writes the full Windows file-association registration (ProgID, DefaultIcon, OpenWithProgids, Capabilities, RegisteredApplications) touching only the program's own keys, supports a subset/unknown extensions, refuses when run from source, and uses the frozen exe as launcher.

| # | Input | Expected output |
|---|---|---|
| 1 | register(None, fakeReg, EXE) | returns list(FILE_TYPES); shell/open/command == '"...OpenGenEdit.exe" "%1"'; DefaultIcon ==',0'; FriendlyAppName=='OpenGenEdit'; every ext has OpenWithProgids/SupportedTypes/FileAssociations; RegisteredApps set; reg.changes==1 |
| 2 | inspect reg.writes | every write path under our 4 roots or endsWith OpenWithProgids; no 'userchoice' written |
| 3 | register(['.gb','fa'], reg, EXE) | returns ['.gb','.fa']; status .gb==THIS, .dna==NONE; not registered-for-current-exe |
| 4 | register(['.txt']) | raises ValueError; snapshot unchanged (nothing written) |
| 5 | is_frozen->False: launcher() / register(None) | launcher None; RuntimeError match 'program file'; reg.keys=={}; launcher(dev=True) ends main.py; command ends '"%1"' |
| 6 | sys.frozen=True | launcher()==[resolved sys.executable] |

## `test_status_default_and_settings_page`

**Result:** ✅ passed in 0.00 s · [line 162](../../tests/test_file_assoc.py#L162) · **Area:** File associations: status/stale/default · **Checked against:** hand-computed · **Strength:** 🟢 strong

status() round-trips NONE/THIS/OTHER, detects a stale (moved-exe) registration ignoring path case, is_default_for reads UserChoice, and the default-apps settings page opens with a deep-link fallback.

| # | Input | Expected output |
|---|---|---|
| 1 | unregistered status | all NONE; not registered, not stale |
| 2 | register EXE; status | all THIS; registered; not stale |
| 3 | status with MOVED launch path | all OTHER; is_stale True; register(MOVED) repairs; command starts '"D:' |
| 4 | register EXE, query with EXE[0].upper() | is_registered_for_current_exe True (case ignored) |
| 5 | set UserChoice ProgId==PROGID for .gb | is_default_for('.gb') True, '.fa' False |
| 6 | open_default_apps_settings(opener) | opens 'ms-settings:defaultapps?registeredAppUser=OpenGenEdit'; opener that rejects '?' falls back to 'ms-settings:defaultapps'; opener that always raises -> returns False |

## `test_unregister`

**Result:** ✅ passed in 0.00 s · [line 213](../../tests/test_file_assoc.py#L213) · **Area:** File associations: unregister() · **Checked against:** round-trip/self-consistency · **Strength:** 🟢 strong

unregister() removes exactly the keys register() created (restoring other programs' shared keys byte-for-byte), leaves no trace on an otherwise-empty registry, and removes extensions an older version registered.

| # | Input | Expected output |
|---|---|---|
| 1 | pre-seed other programs' .gb/OpenWithProgids, Other.exe, QSettings, RegisteredApps, UserChoice; snapshot; register EXE; unregister | unregister True; snapshot == pre-seeded others (exact); status all NONE; second unregister False; reg.changes==2 |
| 2 | register then unregister on empty registry | empty_skeleton(reg) True (no PROGID/capabilities/openwithprogids/.ext keys left) |
| 3 | register, add .old to Capabilities FileAssociations + .old OpenWithProgids, unregister | Software\Classes\.old no longer exists |

## `test_should_ask_decision_table`

**Result:** ✅ passed in 0.01 s · [line 253](../../tests/test_file_assoc.py#L253) · **Area:** File associations: should_ask() decision table · **Checked against:** hand-computed · **Strength:** 🟢 strong

should_ask() returns the right prompt ('first'/'stale'/None) across Windows/source/frozen, already-registered, 'don't ask again', asked-once, and stale-only-after-accept cases, and persisted ini booleans survive a reload.

| # | Input | Expected output |
|---|---|---|
| 1 | frozen=False / windows=False combos | None; plain Windows+frozen -> 'first' |
| 2 | already registered for EXE | None |
| 3 | dont_ask='true' (even with accepted+stale MOVED) | None |
| 4 | asked='true' | 'first' before set, None after |
| 5 | register MOVED (stale); accepted='true' | None before accept, 'stale' after; stale_asked_for==command_for(EXE) -> None; new path E:\Third -> 'stale' again |
| 6 | write asked='true', sync, reopen settings | ask -> None (booleans survive reload) |

## `test_startup_question_screens`

**Result:** ✅ passed in 0.05 s · [line 312](../../tests/test_file_assoc.py#L312) · **Area:** File associations: startup dialog flow · **Checked against:** GUI state · **Strength:** 🟢 strong

The start-up association dialog (driven without exec) does nothing when there's nothing to ask, and Yes/Not now/Don't ask again/stale-repair/failed-registration each register, remember, and message correctly.

| # | Input | Expected output |
|---|---|---|
| 1 | ask_at_startup with frozen=False then windows=False | None; _open_boxes empty |
| 2 | Yes button click | registered; ask->None; ini accepted=='true'; a 'Default apps' info box with 'click OpenGenEdit'; its open button calls open_default_apps_settings (opened==[1]) |
| 3 | Not now (later) click | nothing written (status .gb NONE); ask None; dont_ask unset |
| 4 | Don't ask again (never) click | ask None; dont_ask=='true'; reg.keys=={} |
| 5 | accepted + register MOVED, stale box | box text has 'moved'; Yes registers EXE; later-click on third path remembered (ask None) |
| 6 | Refusing registry (set raises PermissionError), register_now | returns False; a box contains 'denied by policy'; ini accepted stays None |

## `test_manage_dialog`

**Result:** ✅ passed in 0.01 s · [line 372](../../tests/test_file_assoc.py#L372) · **Area:** File associations: manage dialog · **Checked against:** GUI state · **Strength:** 🟢 strong

The FileAssociationsDialog enables register/update/remove buttons correctly as state changes (register, moved-exe->Update, remove), disables all when run from source, and places the default settings file in the data dir.

| # | Input | Expected output |
|---|---|---|
| 1 | open dialog unregistered | 'not yet offered'; register enabled, remove disabled |
| 2 | click register | registered; 'registered for this program file'; register disabled, remove enabled |
| 3 | set launch=MOVED, refresh | register button text 'Update', enabled |
| 4 | click remove | empty_skeleton(reg); ini accepted=='false'; remove disabled |
| 5 | is_frozen->False dialog | 'from source'; both buttons disabled |
| 6 | settings_path() | parent == data_dir() (nothing created) |
