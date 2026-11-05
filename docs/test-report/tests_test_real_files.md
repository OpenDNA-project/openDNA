# `tests/test_real_files.py`

[← all test files](README.md) · source: [`tests/test_real_files.py`](../../tests/test_real_files.py)

## `test_plasmid_files`

**Result:** ✅ passed in 0.01 s · [line 28](../../tests/test_real_files.py#L28) · **Area:** File I/O: SnapGene .dna and GenBank · **Checked against:** independent implementation · **Strength:** 🟢 strong

Checks that a real SnapGene .dna file and its GenBank version hold the same sequence, and that an edit inside a feature moves it correctly through a GenBank save that Biopython reads back.

| # | Input | Expected output |
|---|---|---|
| 1 | demo/pFG051_Chlor.dna and pFG051_Chlor.gb | both 3780 bp, circular, same sequence (case-insensitive); .dna has 15 primers |
| 2 | demo/pBAD30_eilR.gb | circular, 5936 bp |
| 3 | insert 'GATTACA' 10 bases into the first single-span non-source feature > 20 bp, save_genbank | feature spans (start, end+7); modified True then False; Biopython reads same sequence, same feature count, that feature at (start, end+7) |
| 4 | undo | feature back to (start, end) |

**Review notes:**

- Line 83: runtime pytest.skip via need() if demo files are missing (they are tracked in git). Module docstring line 1 says the files are 'kept out of git', which is stale.
