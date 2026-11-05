# `tests/test_copy_and_search_strand.py`

[← all test files](README.md) · source: [`tests/test_copy_and_search_strand.py`](../../tests/test_copy_and_search_strand.py)

## `test_copy_and_copy_reverse_complement_keep_annotations`

**Result:** ✅ passed in 0.06 s · [line 50](../../tests/test_copy_and_search_strand.py#L50) · **Area:** Copy / paste with features · **Checked against:** hand-computed · **Strength:** 🟡 adequate

Copy Reverse Complement pastes the reverse-complemented bases with features flipped; plain Copy keeps them unflipped; works across the origin and twice gives the original; the selection menu offers both.

| # | Input | Expected output |
|---|---|---|
| 1 | random 300 bp (seed 9), linear; genes 61..100 'top', complement(111..140) 'bottom', CDS join(51..60,121..130) 'split'; select(40,160), copy RC, paste at 10 in 50 A | clipboard text and pasted bases == revcomp(SEQ[40:160]); labels {top,bottom,split} |
| 2 | same paste | strands top -1, bottom +1, split -1; top at complement(71..110); split starts 'complement(join(' |
| 3 | plain copy of 40..160 pasted at 0 | top strand +1, bottom -1 |
| 4 | circular, select(280, 320) (280..299 + 0..19), copy RC, paste, copy RC again | text == revcomp(SEQ[280:]+SEQ[:20]); second RC == SEQ[280:]+SEQ[:20] |
| 5 | selection menu on select(40, 80) | actions 'Copy\tCtrl+C' and 'Copy Reverse Complement\tCtrl+Shift+C' |

**Review notes:**

- comment line 64 says each feature's bases read along its strand are unchanged, but only the strand is checked (line 65-66); helper own_bases (line 35) is unused and misnamed
- 'bottom' location and the join parts of 'split' are not checked (line 70 only startswith)
- across-origin: 'features included' (line 87) but the selection holds no features and paste_mime result is not asserted (line 90)
- menu text read via QTimer.singleShot(50) during a modal popup (line 106): if no popup is active the callback errors and the menu may never close

## `test_search_hit_strand`

**Result:** ✅ passed in 0.55 s · [line 122](../../tests/test_copy_and_search_strand.py#L122) · **Area:** GUI: search highlighting · **Checked against:** GUI state · **Strength:** 🟢 strong

A search hit is painted only on the strand it matched (both if none), and the editor records each hit's strand (0 for none or a palindrome).

| # | Input | Expected output |
|---|---|---|
| 1 | SequenceView 900x400; set_search_hits [(10,20,1)], [(10,20,-1)], [(10,20)]; compare pixel rows with a no-hit grab | all three change rows; top rows all above bottom rows; no-strand rows == top \| bottom |
| 2 | set_search_hits [(10,20,-1), (40,50)] | search_strand {(10,20): -1, (40,50): 0} |
| 3 | [(5,15,1), (5,15,-1)] (palindrome) | search_strand[(5,15)] == 0 |
| 4 | set_search_hits([]) | search_strand == {} and search_spans == [] |

**Review notes:**

- section 2 header says the search bar passes the strand, but the search bar is never used; set_search_hits is called directly (line 148-150)
