# `core/src/alignments.rs`

[← all test files](README.md) · source: [`core/src/alignments.rs`](../../core/src/alignments.rs)

## `reads_mismatches_and_coverage`

**Result:** ✅ passed · [line 505](../../core/src/alignments.rs#L505) · **Area:** BAM reads, coverage, pileup and variants · **Checked against:** hand-computed · **Strength:** 🟢 strong

**Why it matters:** Imported BAM reads could be drawn with wrong positions, strands, mismatches, indels or soft clips, and the coverage, pileup and variant tracks could show wrong numbers.

A small synthetic BAM is indexed, exported as SAM, and decoded into reads with correct spans, strand, mismatches, insertions, deletions, soft clips, coverage, pileup counts and called variants.

| # | Input | Expected output |
|---|---|---|
| 1 | BAM on 1000 bp 'ref' with 5 reads: 10M@10, 10M@15 reverse with G at 17, 5M2I5M@30, 5M3D5M@50, 2S8M@900; Bam::open | .bai index created next to the BAM |
| 2 | write_sam | text has '@SQ SN:ref LN:1000', 5 non-header lines, a line 'r2 0 ref 31 60 5M2I5M'; sequences() = [('ref', 1000)] |
| 3 | reads('ref', 0, 100, lowercase 'a' reference, limit 100) | 4 reads, not truncated; r0 = (10, 20, fwd, mapq 60); r1 reverse with mismatches [(17,'G')]; r2 end 40, insertions [(35,2)], no mismatches; r3 end 63, deletions [(55,58)]; limit 2 -> truncated=true |
| 4 | reads 890..1000 | read at 900..908, 0 mismatches (soft-clipped CC ignored) |
| 5 | coverage 10..30 in 4 bins; coverage 0..1024 in 2 bins (block path) | [1.0, 2.0, 1.0, 0.0]; has_depth false before, true after; [40/512, 8/512]; unknown 'chrX' -> error |
| 6 | pileup 15..20; pileup 50..60; pileup with min_mapq 61 | pos 17 = [A1,C0,G1,T0,N0,del0,ins0]; DEL count at 55 = 1; all zeros with mapq filter |
| 7 | variants(all-A ref, min_depth 1, min_frac 0.5); min_depth 3 | [(17,snp,A,G,1,2), (35,ins,'',TT,1,1), (55,del,AAA,'',1,1)]; empty with min_depth 3 |

**Review notes:**

- Writes fixed file names ge_test_reads.bam/.sam into the shared system temp dir and never deletes them; two concurrent cargo test runs could clobber each other.
- Overlaps tests/test_core_formats.py::test_bam_reads_pileup_coverage_and_variants_match_a_model, which covers the same BAM decoding from Python with a model check on random reads; this is the faster unit-level copy.
