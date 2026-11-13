# `tests/test_file_assoc.py`

[← all test files](README.md) · source: [`tests/test_file_assoc.py`](../../tests/test_file_assoc.py)

## `test_register_writes_only_our_keys_and_status_reads_them_back`

**Result:** ✅ passed in 0.00 s · [line 99](../../tests/test_file_assoc.py#L99) · **Area:** Windows file associations: register · **Checked against:** in-test reference model · **Strength:** 🟢 strong

**Why it matters:** Without it double-clicking a .gb/.fa file might not open openDNA, or registering could overwrite other programs' entries or the user's default-app choice in the registry.

Against a dictionary-backed fake registry, checks the keys register writes (ProgID command and icon, OpenWithProgids, SupportedTypes, Capabilities, RegisteredApplications), that nothing outside them or UserChoice is touched, refusals (unknown extension, running from source), status and stale detection for a moved exe, case-insensitive path matching, UserChoice reading, the default-apps settings link with fallback, and removal of the old OpenGenEdit registration.

| # | Input | Expected output |
|---|---|---|
| 1 | register(None, reg, ['C:\Apps\openDNA\openDNA.exe']) | Returns all FILE_TYPES; command '"...openDNA.exe" "%1"', icon '...,0', FriendlyAppName 'openDNA'; per-extension entries; RegisteredApplications -> CAPABILITIES; shell notified once |
| 2 | Every recorded write | Only under our ProgID, Applications\openDNA.exe, CAPABILITIES, RegisteredApplications (name 'openDNA') or OpenWithProgids (our ProgID); never UserChoice |
| 3 | register(['.gb','fa']); then register(['.txt']) | ['.gb','.fa']; status .gb THIS, .dna NONE; not fully registered; .txt -> ValueError and registry unchanged |
| 4 | Not frozen: launcher(), register(None, reg) | None; RuntimeError 'program file'; registry empty; launcher(dev=True) ends with main.py and its command ends with "%1" |
| 5 | sys.frozen=True: launcher() | [resolved sys.executable] |
| 6 | Status before/after register, then with the exe moved to D:\Elsewhere | All NONE -> all THIS -> all OTHER with is_stale True; registering again for the moved path fixes it, command starts with "D: |
| 7 | is_registered_for_current_exe with the path upper-cased | True |
| 8 | is_default_for('.gb') before/after UserChoice ProgId set | False / True; .fa still False |
| 9 | open_default_apps_settings with a working opener / an opener refusing '?' URIs / one always failing | 'ms-settings:defaultapps?registeredAppUser=openDNA', True / falls back to 'ms-settings:defaultapps', True / False |
| 10 | Registry holding OpenGenEdit ProgID, OpenWithProgids and RegisteredApplications entries, then register | has_legacy True before, False after; old entries removed |

**Review notes:**

- The fake registry is the test's own model of winreg; the real WinRegistry wrapper is not exercised here, so a bug in that thin layer would go unnoticed.

## `test_unregister_removes_exactly_what_we_added`

**Result:** ✅ passed in 0.00 s · [line 208](../../tests/test_file_assoc.py#L208) · **Area:** Windows file associations: unregister · **Checked against:** round-trip/self-consistency · **Strength:** 🟢 strong

**Why it matters:** Without it removing the file associations could delete other programs' registry keys or the user's settings, or leave openDNA entries behind.

Seeds the fake registry with other programs' keys, registers then unregisters, and checks the registry returns exactly to the earlier snapshot; also on an empty registry and for an extension added by an older version.

| # | Input | Expected output |
|---|---|---|
| 1 | Other programs' .gb ProgID and default, another app, openDNA QSettings key, another RegisteredApplications entry, UserChoice; register then unregister | unregister True; snapshot identical to before; status all NONE |
| 2 | unregister again | False; shell notified only twice in total |
| 3 | register + unregister on an empty registry | Only empty shared parent keys remain |
| 4 | Extra '.old' extension under our Capabilities and OpenWithProgids, then unregister | Software\Classes\.old removed |

## `test_startup_question_decision_table`

**Result:** ✅ passed in 0.04 s · [line 242](../../tests/test_file_assoc.py#L242) · **Area:** Windows file associations: start-up question · **Checked against:** hand-computed · **Strength:** 🟢 strong

**Why it matters:** Without it the app could nag users at every start about file associations, ask when running from source or off Windows, ignore 'don't ask again', or fail to repair the association after the exe was moved.

Walks the should_ask / ask_at_startup decision table with the fake registry and a private ini file: when to return 'first', 'stale' or nothing, and the silent repair cases.

| # | Input | Expected output |
|---|---|---|
| 1 | frozen=False and/or windows=False | None |
| 2 | Frozen on Windows, nothing registered, fresh ini | 'first' |
| 3 | Already registered for this exe | None |
| 4 | dont_ask=true (also with accepted=true and a stale registration) | None |
| 5 | asked=true | None |
| 6 | Registered for another path: never accepted / accepted / stale_asked_for = this path / launched from a third path | None / 'stale' / None / 'stale' |
| 7 | asked=true written to an ini file and reloaded | None (the boolean survives the reload) |
| 8 | accepted=true, registered for another path, ask_at_startup | None, and the registration now points to this exe (silent repair) |
| 9 | asked=true but an OpenGenEdit registration present | should_ask 'stale'; ask_at_startup repairs silently and removes the legacy keys; then None |
