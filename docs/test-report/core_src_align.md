# `core/src/align.rs`

[← all test files](README.md) · source: [`core/src/align.rs`](../../core/src/align.rs)

## `exact_hit_forward_and_reverse`

**Result:** ✅ passed · [line 499](../../core/src/align.rs#L499) · **Area:** Quick alignment (BLAST-like) · **Checked against:** hand-computed · **Strength:** 🟢 strong

A query copied exactly from the template is found at the right place with 100% identity, on the forward strand and as its reverse complement.

| # | Input | Expected output |
|---|---|---|
| 1 | template = random_dna(2000, seed 1); query = t[500..800]; linear, max_hits=10 | best hit start=500, end=800, forward=true, identity=1.0, qstart..qend = 0..300 |
| 2 | query = revcomp(t[500..800]) | best hit 500..800, forward=false; exactly 1 hit |

**Review notes:**

- Forward case (line 502-505) does not assert hits.len()==1, so spurious extra hits would pass; the reverse case does (line 509).

## `mismatch_insertion_and_deletion`

**Result:** ✅ passed · [line 513](../../core/src/align.rs#L513) · **Area:** Quick alignment (BLAST-like) · **Checked against:** hand-computed · **Strength:** 🟡 adequate

A query with a SNP, a 3-base insertion and a 5-base deletion still aligns over the whole template region and reports each difference, on both strands.

| # | Input | Expected output |
|---|---|---|
| 1 | template = random_dna(3000, seed 2); query = t[1000..1300] with bases 250..255 deleted, 'GGG' inserted at 150, base 50 mutated (A->C, else ->A) | best hit forward, 1000..1300; mismatches == [(1050, new base)]; sum of insertion lengths = 3; sum of deletion lengths = 5; identity > 0.95 |
| 2 | revcomp of the edited query | best hit forward=false, 1000..1300 |

**Review notes:**

- Positions of the insertion (template 1150) and deletion (1250..1255) are never asserted, only their total sizes (lines 525-526); a misplaced indel would pass.
- Reverse-strand case (line 530) checks only strand and span, not the mismatch/indel reports.

## `across_the_origin`

**Result:** ✅ passed · [line 534](../../core/src/align.rs#L534) · **Area:** Quick alignment (BLAST-like) · **Checked against:** hand-computed · **Strength:** 🟡 adequate

On a circular template a query spanning the origin is reported as one hit whose end passes the sequence length; a linear template does not invent such a hit.

| # | Input | Expected output |
|---|---|---|
| 1 | template = random_dna(1000, seed 3), circular; query = t[900..] + t[..100] (200 bp) | best hit start=900, end=1100 |
| 2 | same query, linear template | best hit end <= 1000 |

**Review notes:**

- Comment line 540 says 'only the longer half is found' but both halves are 100 bp; the linear check (line 542) is only end <= 1000 and does not assert the number of hits (expected 2) or their coordinates.
- Circular case does not check forward or qstart/qend.

## `no_hit_for_random_query`

**Result:** ✅ passed · [line 546](../../core/src/align.rs#L546) · **Area:** Quick alignment (BLAST-like) · **Checked against:** property/invariant · **Strength:** 🟡 adequate

An unrelated random query gives no hits (no false positives above the score threshold).

| # | Input | Expected output |
|---|---|---|
| 1 | template = random_dna(5000, seed 4); query = random_dna(500, seed 99); linear | empty hit list |

**Review notes:**

- Single seed pair only; a threshold regression that lets chance hits through on other seeds could pass.

## `query_hanging_off_the_template_end_aligns_only_the_overlap`

**Result:** ✅ passed · [line 556](../../core/src/align.rs#L556) · **Area:** Quick alignment (BLAST-like) · **Checked against:** hand-computed · **Strength:** 🟢 strong

A query that runs past either end of a linear template is aligned only over the overlapping bases, with exact template and query coordinates on both strands.

| # | Input | Expected output |
|---|---|---|
| 1 | template = random_dna(1000, seed 11); query = t[960..] + random_dna(60, seed 12) | exactly 1 hit: 960..1000, forward, qstart..qend = 0..40 |
| 2 | revcomp of that query | best hit 960..1000, forward=false, qstart..qend = 60..100 |
| 3 | query = random 60 bp + t[..40] | best hit 0..40, qstart..qend = 60..100 |

## `overhanging_query_wraps_on_a_circle_and_a_tiny_overlap_is_no_hit`

**Result:** ✅ passed · [line 579](../../core/src/align.rs#L579) · **Area:** Quick alignment (BLAST-like) · **Checked against:** hand-computed · **Strength:** 🟢 strong

End+start of the template is one hit across the origin on a circle and two pieces on a line; a 15-base overlap is below the reporting threshold.

| # | Input | Expected output |
|---|---|---|
| 1 | template = random_dna(1000, seed 11); query = t[960..] + t[..60]; circular | exactly 1 hit: start 960, end 1060, qstart..qend 0..100 |
| 2 | same query, linear | 2 hits |
| 3 | query = t[985..] (15 bp) + random_dna(60, seed 12); linear and circular | no hit (min score min(1.2*75, 40) = 40 points, i.e. 20 matching bases) |

**Review notes:**

- Linear case (line 587) checks only the count 2, not the coordinates of the two halves.

## `every_piece_of_a_patchwork_query_is_reported_including_the_inverted_one`

**Result:** ✅ passed · [line 636](../../core/src/align.rs#L636) · **Area:** Quick alignment (structural variants) · **Checked against:** hand-computed · **Strength:** 🟢 strong

A query built from 4 template pieces (one inverted) plus a 200-bp random insert yields one hit per piece with correct strand, query range and template range.

| # | Input | Expected output |
|---|---|---|
| 1 | template random_dna(2000, seed 21); query = A(100@100) + 200 random + B(600) + C(1100, inverted) + D(1600); linear and circular | 4 pieces (q 0..100 fwd t100..200; q300..400 fwd 600..700; q400..500 rev 1100..1200; q500..600 fwd 1600..1700), each coordinate within +-2 |
| 2 | pieces in another order: A=1400, B=900, C=500 (inverted), D=100; seed 6 | same 4-piece pattern at those template places, +-2 |
| 3 | pieces close together: 500, 650, 800, 950; seed 7 | same 4-piece pattern, +-2 |

## `pieces_on_one_diagonal_with_a_long_unmatched_middle_are_not_merged`

**Result:** ✅ passed · [line 651](../../core/src/align.rs#L651) · **Area:** Quick alignment (structural variants) · **Checked against:** hand-computed · **Strength:** 🟢 strong

Pieces lying on the same diagonal but separated by long unmatched stretches are reported separately instead of as one long low-identity hit.

| # | Input | Expected output |
|---|---|---|
| 1 | template random_dna(2000, seed 22); sv_query(A=500, B=800, C=1500 inverted, D=1000, seed 8): A, B, D share diagonal 500 | 4 pieces as listed (+-2); every hit spans <= 102 template bases and has identity > 0.97 |

## `max_hits_still_limits_the_number_of_pieces_and_keeps_the_best`

**Result:** ✅ passed · [line 663](../../core/src/align.rs#L663) · **Area:** Quick alignment (structural variants) · **Checked against:** property/invariant · **Strength:** 🔴 weak

max_hits caps the number of reported pieces and hits are sorted best first.

| # | Input | Expected output |
|---|---|---|
| 1 | 4-piece query (seed 9) on random_dna(2000, seed 23); max_hits = 2, 1, 0 | 2, 1 and 0 hits |
| 2 | max_hits = 10 | scores non-increasing |

**Review notes:**

- Name promises 'keeps the best' but lines 666-667 only check counts; nothing compares align(..,2) with the top 2 of align(..,10), so returning the 2 worst pieces would pass.

## `a_piece_repeated_in_the_query_is_found_twice_on_the_same_template_place`

**Result:** ✅ passed · [line 674](../../core/src/align.rs#L674) · **Area:** Quick alignment (structural variants) · **Checked against:** hand-computed · **Strength:** 🟢 strong

A tandem duplication in the query (same 100 bp twice) gives two hits on the same template place, not one deduplicated hit.

| # | Input | Expected output |
|---|---|---|
| 1 | template random_dna(2000, seed 24); query = 60 random + t[700..800] + 60 random + t[700..800] + 60 random | 2 pieces: q 60..160 -> t 700..800 fwd, q 220..320 -> t 700..800 fwd (+-2) |

## `an_unmatched_stretch_is_cut_out_but_a_short_one_is_still_bridged`

**Result:** ✅ passed · [line 686](../../core/src/align.rs#L686) · **Area:** Quick alignment (structural variants) · **Checked against:** hand-computed · **Strength:** 🟢 strong

A short unmatched stretch (8 bp) is bridged into one hit, while 40 bp or 120 bp unmatched stretches split the alignment into two pieces (MAX_DRAWDOWN rule).

| # | Input | Expected output |
|---|---|---|
| 1 | template random_dna(2000, seed 25); query = t[300..360] + 8 random + t[368..428] | exactly 1 hit, 300..428 |
| 2 | t[300..360] + 40 random + t[400..460] | 2 pieces: q0..60 -> 300..360, q100..160 -> 400..460 (+-2) |
| 3 | t[300..360] + 120 random + t[480..540] | 2 pieces: q0..60 -> 300..360, q180..240 -> 480..540 (+-2) |

**Review notes:**

- Comment line 689 says the 8 changed bases are 'shown as mismatches' but the mismatch list is not asserted.
