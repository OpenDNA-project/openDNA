# `tests/test_editing.py`

[← all test files](README.md) · source: [`tests/test_editing.py`](../../tests/test_editing.py)

## `test_500_random_edits_keep_bases_and_features_then_undo_redo_exactly`

**Result:** ✅ passed in 24.18 s · [line 15](../../tests/test_editing.py#L15) · **Area:** Editing / annotation drift / undo · **Checked against:** in-test reference model · **Strength:** 🟢 strong

Property test: random edits of every kind must keep bases equal to a plain string and every feature on exactly its own bases, and undo/redo must restore every snapshot exactly. This protects users' annotations during editing.

| # | Input | Expected output |
|---|---|---|
| 1 | selftest.run_random_edits(seed 0..2, linear and circular, 500 steps) on a 300 bp doc with many feature shapes; edits: replace, insert, delete, overwrite, paste with features, add/delete feature, reverse complement, and on circles edits across the origin and rotate | after each edit: bases == string model, topology unchanged, exactly one undo step, features match the per-base IdModel; full undo and redo give identical snapshots; returns 500 |
