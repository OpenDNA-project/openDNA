# `core/src/tracks.rs`

[← all test files](README.md) · source: [`core/src/tracks.rs`](../../core/src/tracks.rs)

## `gff_bed_and_queries`

**Result:** ✅ passed · [line 245](../../core/src/tracks.rs#L245) · **Area:** Annotation tracks (GFF/GTF/BED/VCF) · **Checked against:** hand-computed · **Strength:** 🟡 adequate

**Why it matters:** Loaded GFF/GTF/BED/VCF annotation tracks could show features off by one base, lose escaped names, or miss features in a region.

GFF3/GTF, BED and VCF files load into an index: region lines skipped, ##FASTA ends parsing, URL-decoding, 1-based to 0-based, overlap queries via max_end, limit, density, name search, and bad lines rejected.

| # | Input | Expected output |
|---|---|---|
| 1 | GFF3: region chr1 1-1000 (skipped), genes lacZ 10-100 +, 'big%2Cgene' 50-900 -, CDS araC 200-300, GTF-style chr2 tiny 5-6, then ##FASTA | sequences [(chr1,3),(chr2,1)] |
| 2 | query chr1 150..160; 0..1000 | ['big,gene']; ['lacZ','big,gene','araC'] |
| 3 | query chr1 0..1000 limit 2; chr2 4..5; chr1 9..10; density 0..1000 10 bins; search 'ARA' | truncated=true; 'tiny'; start 9; bin 0 = 2; 'araC' |
| 4 | BED 'track name=x', 'chr1 0 10 first 0 -', 'chr1 20 30' (spaces); query 5..25; then 'chr1 not-a-number 10' | 2 hits, first strand -1; open error |
| 5 | VCF: chr1 100 A>G; chr1 200 rs7 ACGT>A | (99, 100, 'A>G'); end 203, 'rs7 ACGT>A', kind 'variant' |

**Review notes:**

- Fixed file names ge_test.gff3/.bed/.vcf in the shared temp dir, never removed.
- Overlaps tests/test_core_formats.py::test_tracks_gff_bed_parsing_query_density_and_search (GFF/GTF/BED with a brute-force query check); only the VCF part is covered here alone.

**Fixed by the audit of November 2026:**

- The limit check compared a query with itself; it now asserts the two items kept (lacZ, big,gene), the cut flag, and no cut at a limit of 3.
