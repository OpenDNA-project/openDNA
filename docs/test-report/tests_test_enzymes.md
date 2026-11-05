# `tests/test_enzymes.py`

[← all test files](README.md) · source: [`tests/test_enzymes.py`](../../tests/test_enzymes.py)

## `test_cut_finder_and_the_common_list`

**Result:** ✅ passed in 0.00 s · [line 34](../../tests/test_enzymes.py#L34) · **Area:** Restriction enzymes: cut finder · **Checked against:** hand-computed · **Strength:** 🟢 strong

The cut finder handles degenerate (IUPAC) sites, 3' overhangs, enzymes cutting on both sides, reverse-strand sites with mirrored cuts and sites across a circular origin; the common enzyme list keeps its cut positions and descriptions.

| # | Input | Expected output |
|---|---|---|
| 1 | AccI GTMKAC (2,4) on 'AAAGTAGACTTTGTCTACAAAGTGGAC' | [(5,7), (14,16)] (GTGGAC does not match) |
| 2 | BglI GCCNNNNNGGC (7,4) on 'TTGCCAAGTTGGCTT' | [(9, 6)] (3' overhang) |
| 3 | BcgI CGANNNNNNTGC ((-10,-12),(24,22)) at position 40 in random flanks (seed 11) | [(30,28), (64,62)] |
| 4 | BsaI GGTCTC (7,11): 'AA'+'GGTCTC'+A*20 and A*20+'GAGACC'+'AA' | [(9,13)] and [(15,19)] (4 nt apart) |
| 5 | EcoRI on circular 'ATTCAAAAAAGA' | [(11, 3)] |
| 6 | COMMON list entries | EcoRI (1,5), KpnI (5,1), EcoRV blunt, BsaI (7,11), SapI (8,11) |
| 7 | find_sites('TTGAATTCAAGAGACCTTTTTTCAA', linear) | {'EcoRI': [3], 'BsaI': [5]} |
| 8 | describe() for EcoRI, KpnI, EcoRV, BcgI, BsaI, BbsI, SapI | 'GAATTC: G^AATTC, 5′ overhang 4 nt'; ends '3′ overhang 4 nt'; 'blunt'; 'more place'; 'GGTCTC: GGTCTC(1/5), 5′ overhang 4 nt'; 'GAAGAC: GAAGAC(2/6),'; 'GCTCTTC: GCTCTTC(1/4), 5′ overhang 3 nt' |

## `test_the_digest_dialog_filter_gel_ladders_and_lanes`

**Result:** ✅ passed in 0.06 s · [line 81](../../tests/test_enzymes.py#L81) · **Area:** GUI: Digest dialog and gel · **Checked against:** hand-computed · **Strength:** 🟢 strong

The Digest tool lists cutting enzymes with correct columns, filters by name/site/criteria while keeping ticks, offers five gel ladders chosen by clicking, and keeps several side-by-side digests (lanes) each with its own enzymes, fragment table and gel lane.

| # | Input | Expected output |
|---|---|---|
| 1 | random 3006 bp with one GAATTC at 1500 (seed 11, re-drawn until unique) | EcoRI listed; default filter EnzymeFilter(1, None, 1); list == cutters, sorted by name |
| 2 | EcoRI row; EcoRV row if present | ['1', 'GAATTC', '5′ 4 nt'], Owned unchecked; EcoRV ends 'blunt' |
| 3 | column widths / sort by sites descending | widths == column_need(c); max width == sum + scrollbar + frame + 4; counts sorted numerically |
| 4 | filter 'ecor', 'GAATTC'; max_sites=1; ends='5' | EcoRI kept and every listed enzyme matches the filter |
| 5 | tick EcoRI, filter 'zzzz', then clear | list empty, chosen ['EcoRI'], cut marks [(1501,1505)]; still ticked; tooltip has 'G^AATTC' and '5′ overhang 4 nt' |
| 6 | LADDERS | ranges 1 kb 250-10000, 1 kb Plus 250-25000, 100 bp 100-1500, 100 bp Plus 100-3000, 50 bp 50-500; bright band in each |
| 7 | gel 300x400, menu action '100 bp', mouse move at x=40 | over_ladder(40) true; ladder '100 bp', signal emitted, lane 0 renamed; pointing-hand cursor |
| 8 | seq with one EcoRI and one BamHI; tick EcoRI; add_lane; tick BamHI+EcoRI | 2 frags, tab 'EcoRI'; new lane uncut; 3 frags summing to len(seq), tab 'BamHI+EcoRI'; gel lanes named; lane_lines('A+B+C+D+E')==['A','B','+3 more'] |
| 9 | back to tab 0, close it, remove_lane(0) | EcoRI-only state restored; one lane left (not closable); last lane cannot be removed |

**Review notes:**

- Line 97: `item_of('EcoRV') is None or ...` passes vacuously if EcoRV does not cut; it does cut today only because the seeded rng is consumed by the previous test first (run alone, the sequence differs).
- Line 102: widths are compared with table.column_need(c), the same function fit() uses to set them (self-consistency); only line 103's loose upper bound is independent.

**Fixed by the audit of November 2026:**

- The docstring now comes first, so the test has a __doc__.
