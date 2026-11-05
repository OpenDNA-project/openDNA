# `tests/test_view_help.py`

[← all test files](README.md) · source: [`tests/test_view_help.py`](../../tests/test_view_help.py)

## `test_every_picker_view_has_a_guide_and_an_options_entry`

**Result:** ✅ passed in 0.00 s · [line 27](../../tests/test_view_help.py#L27) · **Area:** GUI: per-view help and options (data) · **Checked against:** property/invariant · **Strength:** 🟡 adequate

Every view in the picker has a short guide and an options entry, and every options target names a real settings tab, so the '?' and cog buttons never dead-end.

| # | Input | Expected output |
|---|---|---|
| 1 | Each picker view name in view_help.HELP | 4 to 8 lines, <=130 words total, present in OPTIONS |
| 2 | Key sets | set(HELP) == set(OPTIONS) == set(picker views) |
| 3 | OPTIONS targets of kind 'settings' | target[1] in tool_settings.TOOLS |

**Review notes:**

- L37-39: only 'settings' targets are validated; other target kinds (e.g. the primer-options one) are not checked here.

## `test_cog_and_help_buttons_in_every_pane`

**Result:** ✅ passed in 0.16 s · [line 42](../../tests/test_view_help.py#L42) · **Area:** GUI: pane cog and '?' buttons · **Checked against:** GUI state · **Strength:** 🟡 adequate

In every view, the pane's cog opens the matching options page (or is disabled with 'No options for this view'), and '?' opens the guide popup.

| # | Input | Expected output |
|---|---|---|
| 1 | 2,000-bp A/T document; pane set to each workspace.VIEWS entry; SettingsDialog.exec and primer_options stubbed | cog enabled iff OPTIONS[name] not None; tooltip == options_tip(name) ('No options for this view' when none); click opens tab labelled LABELS[target] or 'Primer options' |
| 2 | Digest pane, click '?' (ViewHelpDialog.exec stubbed) | one dialog shown, text contains 'restriction enzymes'; help tooltip non-empty; cog and help focusPolicy != NoFocus |

**Review notes:**

- L57, L63: expected tooltip and tab label are computed with view_help.options_tip and tool_settings.LABELS, the same functions/tables the code uses, so a wrong mapping in those tables would not be caught (only the wiring is).
- L17: imports helpers close/dna/open_doc from tests/test_cas9.py, coupling the two modules.
