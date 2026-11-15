# `tests/test_selftest.py`

[← all test files](README.md) · source: [`tests/test_selftest.py`](../../tests/test_selftest.py)

## `test_every_selftest_check_passes`

**Result:** ✅ passed in 0.39 s · [line 14](../../tests/test_selftest.py#L14) · **Area:** In-app self-test · **Checked against:** property/invariant · **Strength:** 🟡 adequate

**Why it matters:** Without it the Help > Troubleshooting > Self-Test checks could start failing in a release without anyone noticing before a user runs them.

Runs every check registered in app/selftest.py (editing vs plain string edits, random edits with undo/redo, annotation drift, linearising a plasmid, primer binding mismatches) and fails listing every failing check; also requires at least 9 checks to be registered.

| # | Input | Expected output |
|---|---|---|
| 1 | selftest.CHECKS | At least 9 checks registered |
| 2 | run_check(fn) for each check | All pass; otherwise every failure listed as 'area: what: message' |

**Review notes:**

- The real checks live in app/selftest.py; this test's strength is theirs and cannot be judged from the test file.
- Overlap: tests/test_documents.py::test_random_edits_keep_bases_and_features_then_undo_redo_exactly runs the same selftest.run_random_edits machinery with 500 edits; this test re-runs the short in-app versions. Still worth keeping as the guard that the in-app self-test itself passes.
