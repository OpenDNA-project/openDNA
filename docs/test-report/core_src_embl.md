# `core/src/embl.rs`

[← all test files](README.md) · source: [`core/src/embl.rs`](../../core/src/embl.rs)

## `a_last_record_without_the_end_line_is_read_and_said`

**Result:** ✅ passed · [line 79](../../core/src/embl.rs#L79) · **Area:** EMBL reader · **Checked against:** hand-computed · **Strength:** 🟡 adequate

**Why it matters:** A cut-short EMBL file could lose its last record without a word.

A whole file gives no note; a last record without // is still written (when it holds bases) with a note; a cut second record is read next to the whole first one; a cut record with no bases is not read and is said.

| # | Input | Expected output |
|---|---|---|
| 1 | one record with // | 1 LOCUS, no notes |
| 2 | one record without // | 1 LOCUS, the bases present, note 'The last record has no `//` line: the file may be cut short.' |
| 3 | record, //, then a second record without // | 2 LOCUS, 1 note |
| 4 | an ID line and a DE line only | no GenBank text, note says it holds no bases |

**Review notes:**

- Only the conversion to GenBank text is checked here; that gb-io parses it is covered by the Python read of a cut file.
