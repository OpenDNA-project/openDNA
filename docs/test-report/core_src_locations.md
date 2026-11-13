# `core/src/locations.rs`

[← all test files](README.md) · source: [`core/src/locations.rs`](../../core/src/locations.rs)

## `delete_rules`

**Result:** ✅ passed · [line 248](../../core/src/locations.rs#L248) · **Area:** Feature locations (delete) · **Checked against:** hand-computed · **Strength:** 🟡 adequate

**Why it matters:** Deleting bases could leave features shifted, wrongly trimmed, or not removed when fully deleted, especially multi-part complement features.

Delete rules: shift, trim at start, shrink inside, remove when covered; a complement(join) losing one part collapses to complement(range).

| # | Input | Expected output |
|---|---|---|
| 1 | Range 10..20; delete 0..5, 5..15, 12..14, 10..20 | Some(5..15), Some(5..10), Some(10..18), None |
| 2 | complement(join(0..5, 10..20)); delete 0..5 | Some(complement(5..15)) |

**Review notes:**

- No case for a deletion after the range (unchanged) or overlapping only the end.

## `spans_and_strand`

**Result:** ✅ passed · [line 262](../../core/src/locations.rs#L262) · **Area:** Feature locations (flatten) · **Checked against:** hand-computed · **Strength:** 🟡 adequate

**Why it matters:** Feature locations could be flattened into wrong ranges or the wrong strand, which every drawing and export relies on.

spans() flattens a location into its ranges and reports the strand.

| # | Input | Expected output |
|---|---|---|
| 1 | complement(join(0..5, 10..20)) | ([(0,5),(10,20)], -1) |
| 2 | range 3..4 | ([(3,4)], 1) |

**Review notes:**

- Between (zero-width site) handling is not tested.

## `clip_window_keeps_parts_in_reading_order`

**Result:** ✅ passed · [line 269](../../core/src/locations.rs#L269) · **Area:** Feature locations (clip window / rotate) · **Checked against:** hand-computed · **Strength:** 🟢 strong

**Why it matters:** Copying or rotating across the origin could scramble the order of a feature's parts, lose its strand, or keep parts outside the copied window.

Clipping a location to a window on a circle (copy across origin, rotate) keeps parts in reading order, merges touching parts, shifts, drops or clips ranges, and keeps complement strand.

| # | Input | Expected output |
|---|---|---|
| 1 | join(190..200, 0..10), window start 5, len 200, n 200 | spans [(185,200),(0,5)] |
| 2 | 10..20 window 5+50; 100..120 window 5+50; 0..20 window 10+50 (n 200) | Some(5..15); None; Some(0..10) |
| 3 | complement(join(17..20, 0..3)), window 15+10, n 20 | ([(2,8)], -1) |
| 4 | 8..14 window 10+20 n 20; window 12+20 | [(18,20),(0,4)]; [(16,20),(0,2)] |

**Review notes:**

- Fuzzy-marker dropping on cut ends (locations.rs line 159) is not tested.
