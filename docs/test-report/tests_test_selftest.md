# `tests/test_selftest.py`

[← all test files](README.md) · source: [`tests/test_selftest.py`](../../tests/test_selftest.py)

## `test_selftest_checks`

**Result:** ✅ passed in 0.84 s · [line 18](../../tests/test_selftest.py#L18) · **Area:** In-app self-test · **Checked against:** in-test reference model · **Strength:** 🟡 adequate

Runs every in-app self-test check (Settings > Self-test) under pytest and checks the panel reports them all passed.

| # | Input | Expected output |
|---|---|---|
| 1 | each of the 9 selftest.CHECKS (editing vs string, 200 edits + undo, annotation drift, 150-step random edits x2 seeds, linearize undo, primer mismatches, 5' tail / 3' mismatch) | all return ok; failures listed together |

**Fixed by the audit of November 2026:**

- Removed in the second pass of the audit (fewer tests to review, no biological check lost): the panel re-run was removed (the panel's counting and messages are tested with fake checks in test_app_logic.py); the test now asserts that the 9 checks are registered, so it cannot pass with none.
