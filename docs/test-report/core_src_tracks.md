# `core/src/tracks.rs`

[← all test files](README.md) · source: [`core/src/tracks.rs`](../../core/src/tracks.rs)

## `gff_bed_and_queries`

**Result:** ✅ passed · [line 391](../../core/src/tracks.rs#L391) · **Area:** Annotation tracks (GFF/GTF/BED/VCF) · **Checked against:** hand-computed · **Strength:** 🟡 adequate

**Why it matters:** Loaded GFF/GTF/BED/VCF annotation tracks could show features off by one base, lose escaped names, or miss features in a region.

GFF3/GTF, BED and VCF files load into an index: every type counted (a whole-sequence region flagged), the kept types indexed, ##FASTA ends parsing, URL-decoding, 1-based to 0-based, overlap queries via max_end, limit, density, name search, and bad lines rejected.

| # | Input | Expected output |
|---|---|---|
| 1 | GFF3: region chr1 1-1000, genes lacZ 10-100 +, 'big%2Cgene' 50-900 -, CDS araC 200-300, GTF-style chr2 tiny 5-6, then ##FASTA; keep gene and CDS | before: sequences [(chr1,4),(chr2,1)], types [(CDS,1,not whole),(gene,3,not whole),(region,1,whole)]; keep -> 0 dropped, [(chr1,3),(chr2,1)] |
| 2 | query chr1 150..160; 0..1000 | ['big,gene']; ['lacZ','big,gene','araC'] |
| 3 | query chr1 0..1000 limit 2; chr2 4..5; chr1 9..10; density 0..1000 10 bins; search 'ARA' | truncated=true; 'tiny'; start 9; bin 0 = 2; 'araC' |
| 4 | BED 'track name=x', 'chr1 0 10 first 0 -', 'chr1 20 30' (spaces); query 5..25; then 'chr1 not-a-number 10' | 2 hits, first strand -1; open error |
| 5 | VCF: chr1 100 A>G; chr1 200 rs7 ACGT>A | (99, 100, 'A>G'); end 203, 'rs7 ACGT>A', kind 'variant' |

**Review notes:**

- Fixed file names ge_test.gff3/.bed/.vcf in the shared temp dir, never removed.
- Overlaps tests/test_core_formats.py::test_tracks_gff_bed_parsing_query_density_and_search (GFF/GTF/BED with a brute-force query check); only the VCF part is covered here alone.

**Fixed by the audit of November 2026:**

- The limit check compared a query with itself; it now asserts the two items kept (lacZ, big,gene), the cut flag, and no cut at a limit of 3.
- 13 Nov 2026: gap closed. The proposed checks were added, pass on the current code and fail under the mutation that slipped through before. Was: density is checked only at bin 0 (= 2), and both items in bin 0 start there, so a density that ignores item ends (or miscomputes the end bin) passes; the 'density' claim is effectively untested. \| Other parts (region skip, ##FASTA, URL decoding, 0-based conversion, max_end query, limit flag, BED, VCF, bad line) are pinned exactly.

## `partial_regions_colours_notes_and_features_across_the_origin`

**Result:** ✅ passed · [line 442](../../core/src/tracks.rs#L442) · **Area:** Annotation tracks (GFF/GTF/BED/VCF) · **Checked against:** hand-computed · **Strength:** 🟢 strong

**Why it matters:** A GFF3 file of partial region lines could be read as empty again, a feature across the origin could be refused or split in two, or a colour or note in the file could be lost.

A region that does not cover its sequence is a type worth showing; a line whose end comes before its start is one item on a circle (its end past the length, found from both sides of the origin) and is left out and counted on a linear sequence; GFF colour/colour and Note/product and the BED itemRgb column are kept as written.

| # | Input | Expected output |
|---|---|---|
| 1 | ##sequence-region c 1 1000; region 100-500 color=#257215; region 950-50 colour=185,21,21 Note=across%20the%20origin; gene 600-700 product=lacZ | types [(gene,1,not whole),(region,2,not whole)]; field counts (2 colours, 2 notes); colours [('#257215',1),('185,21,21',1)]; 3 items on c |
| 2 | keep region+gene on a circle of 1000; query 990..1000, 1000..1010, 600..601 | 0 dropped; item 949..1050 note 'across the origin'; 1 item; note 'lacZ' |
| 3 | keep region only on a linear sequence | 1 dropped; 1 item left |
| 4 | BED itemRgb '255,0,0' and '0' | colours [('255,0,0',1)] ('0' is no colour) |
| 5 | one gene 1..5000, no ##sequence-region | types [(gene,1,not whole)]: only a type naming a sequence, or a stated length, makes a line 'whole' |

**Review notes:**

- Fixed file names in the shared temp dir, never removed (as the other track test).
