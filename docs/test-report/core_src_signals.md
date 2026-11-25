# `core/src/signals.rs`

[← all test files](README.md) · source: [`core/src/signals.rs`](../../core/src/signals.rs)

## `bedgraph_and_wig_bins_by_hand`

**Result:** ✅ passed · [line 326](../../core/src/signals.rs#L326) · **Area:** Alignment view tracks · **Checked against:** hand-computed · **Strength:** 🟢 strong

**Why it matters:** Without it, a bedGraph or wig signal track could show wrong values in the Alignment view (shifted by one base, averaged over the wrong bases, or wrong when zoomed out) with nothing to show it.

Reads a bedGraph and a wig written in the test and checks the bins the view draws: the mean over the bases that have a value and the highest value, NaN where there is none, through the raw intervals and through the zoom levels built for wide views; lines that are not valid are counted in the notes, a file with none valid is an error.

| # | Input | Expected output |
|---|---|---|
| 1 | bedGraph c: [0,100)=1, [100,200)=3, [500,600)=2, plus a line 'c x 1 1' and a track line | sequences [(c, 600)]; notes: 1 line skipped (line 5, start, end or value is not a number) |
| 2 | bins(c, 0, 600, 6) | (1,1), (3,3), NaN x3, (2,2) |
| 3 | bins(c, 50, 150, 1); bins(c, 150, 550, 1) | (2, 3) (50 bases at 1, 50 at 3); (2.5, 3) (50 at 3, 50 at 2) |
| 4 | bins(c, 0, 600, 1) (through the 64 bp zoom blocks); bins(other, 0, 10, 2) | (2, 3) = (100 + 300 + 200) / 300 bases; NaN |
| 5 | wig: fixedStep start=1 step=10 span=5 values 4, 8; variableStep span=10 at 101 value 6 | bins(c, 0, 20, 4) = (4,4), NaN, (8,8), NaN; bins(c, 100, 110, 1) = (6, 6) |
| 6 | a file with only 'c 1'; is_bigwig on the bedGraph | error; false |

## `bigwig_and_bigbed_of_the_demo_hub`

**Result:** ✅ passed · [line 355](../../core/src/signals.rs#L355) · **Area:** Alignment view tracks · **Checked against:** hand-computed · **Strength:** 🟢 strong

**Why it matters:** Without it, a bigWig signal or a bigBed annotation track (what UCSC track hubs hold) could be read wrong, or not at all, with nothing to show it.

Reads the demo hub's plus.bw (written by pybigtools: 10 over [0,2000), 30 over [2000,4000), 20 to 5,936) through the bigtools library and checks bins by hand; reads its sites.bb as an annotation track (the BED lines its entries stand for).

| # | Input | Expected output |
|---|---|---|
| 1 | Signal(plus.bw) | is a bigWig, no notes, sequences [(pBAD30_eilR, 5936)] |
| 2 | bins(1000, 3000, 2); bins(0, 5936, 1); bins(1990, 2010, 2); bins(other, 0, 10, 1) | (10,10), (30,30); (20, 30) = (2000x10 + 2000x30 + 1936x20) / 5936; (10,10), (30,30); NaN |
| 3 | Track(sites.bb).query(0, 5936) | (100, 400, +, site_a, 255,0,0), (2500, 2600, -, site_b, 0,0,255), (5000, 5100, +, site_c, 0,128,0) |

**Review notes:**

- The demo bigWig is too small to have zoom levels: the zoom path of a bigWig is checked against the bedGraph of the real Meso_v7 hub only by hand, not here.
