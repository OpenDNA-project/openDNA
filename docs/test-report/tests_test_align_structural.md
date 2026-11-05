# `tests/test_align_structural.py`

[← all test files](README.md) · source: [`tests/test_align_structural.py`](../../tests/test_align_structural.py)

## `test_every_piece_of_a_structural_variant_query_is_found`

**Result:** ✅ passed in 0.09 s · [line 91](../../tests/test_align_structural.py#L91) · **Area:** Core quick aligner (structural variant) · **Checked against:** hand-computed · **Strength:** 🟢 strong

core.align reports every matching piece of a structural-variant query (A-gap-B-C(inverted)-D) separately with exact template/query coordinates and strand, across many template layouts, circular and linear, and symmetrically for the reverse-complemented query.

| # | Input | Expected output |
|---|---|---|
| 1 | 6 LAYOUTS x circular in (False,True); scenario() builds 2000 bp template, query = A+200spacer+B+rc(C)+D, draws retried so no chance base extends a piece | assert_pieces: exactly 4 hits covering each expected (qstart,qend,strand,start,end) within TOL=2; each hit identity>=0.95 and score>=190 |
| 2 | same query typed as rc(query) | same 4 pieces with every strand flipped and query coords mirrored (m-qe,m-qs) |

## `test_alignment_edge_cases`

**Result:** ✅ passed in 0.02 s · [line 105](../../tests/test_align_structural.py#L105) · **Area:** Core quick aligner edge cases · **Checked against:** hand-computed · **Strength:** 🟢 strong

The aligner finds an origin-wrapping piece on a circle (two halves on a line), honours max_hits, never merges collinear pieces over a mismatching gap, reports a tandem-duplicated query piece twice at one template place, and keeps a near-match piece separate from an exact one with its diff fields populated.

| # | Input | Expected output |
|---|---|---|
| 1 | wrap_d scenario, circular=True vs linear | D template span (1950,2050); 4 pieces found with wrap end>=2050; on a line D is two halves (500,550,1,1950,2000)+(550,600,1,0,50) |
| 2 | collinear A/B/D with 200 bp gaps; max_hits 2/1/0 | len 2/1/0 respectively; full align gives 4 hits none spanning the gap (max qlen<=PIECE+TOL); sorted best-first by score |
| 3 | tandem dup: 60 random + t[700:800] + 60 + t[700:800] + 60 | two hits (60,160,1,700,800) and (220,320,1,700,800) at same template place |
| 4 | exact t[900:1000] + 150 random + messy t[300:420] (1 sub, 3 bp deletion) | 2 hits; first identity 1.0 at (900,1000); second ~(300,420) with mismatches (h[7]) and deletions (h[9]) non-empty, identity<1.0 |
