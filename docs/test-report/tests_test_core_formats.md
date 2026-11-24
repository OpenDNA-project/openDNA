# `tests/test_core_formats.py`

[← all test files](README.md) · source: [`tests/test_core_formats.py`](../../tests/test_core_formats.py)

## `test_bam_reads_pileup_coverage_and_variants_match_a_model`

**Result:** ✅ passed in 0.07 s · [line 261](../../tests/test_core_formats.py#L261) · **Area:** BAM reads, pileup, coverage, variants · **Checked against:** in-test reference model · **Strength:** 🟢 strong

**Why it matters:** Without it, the read alignment track could show reads at wrong places, miscount depth, ignore duplicate/QC filters or call wrong SNPs and indels from a user's BAM file.

Writes BAM files with the standard library (hand-designed reads, 300 random reads, a SNP set) and compares core.Bam reads, pileup, coverage and variants with a hand count and a Python model of the CIGAR walk and filters.

| # | Input | Expected output |
|---|---|---|
| 1 | hand BAM chrA 2000 / chrB 500 | sequences() lists both; .bai is written; unknown chrZ raises OSError |
| 2 | reads() over chrA | unmapped/duplicate/QC-fail excluded; perfect (10,60,+1); reverse with mismatch at 45; soft clip starts 100 with no mismatch; deletion (220,225); insertion [(320,3)]; hard clip 700..740; mate position/other-ref flag/insert size for pairs |
| 3 | windows and limits | 100..130 gives only 'softclip'; limit 2 gives 2 reads and more=True; limit 4 gives more=False; lower-case reference same mismatches |
| 4 | pileup 100..150, 215..230, 315..325 | 45 single counts then 0; deleted column for 220..225 counted in depth; insertion flag before 320 |
| 5 | pileup whole chrA; past end; min_baseq 20 / min_mapq 20 / 41 / 61 at 500..550 | equals model; cut to 10 columns; 99 / 50 / 0 / 0 bases |
| 6 | coverage narrow bins (0..200/20, 95..155/6, 200..260/60) and wide bins via blocks | equal mean model depth within 1e-6; has_depth turns True after first wide query; empty window [0,0,0]; 0 bins [] |
| 7 | 300 random reads, 5 filter combinations | pileup and coverage equal the model |
| 8 | SNP BAM: 10 reads, SNPs at 120 (5/10), 125 (1/10), 130 (3/10) | default 0.2: 120 and 130 called with counts (5,10), (3,10); fraction and depth thresholds include/exclude exactly as computed |
| 9 | reference with alt base at 120; lower case ref; N at 120; chrZ | alt/ref swapped; same calls; no call; OSError |
| 10 | hand BAM variants min_depth 1, 0.5 | del (220, ref[220:225]); ins GGG at 320; SNP at 45 with (1,2); nothing from soft clip |
| 11 | random BAM variants min_depth 5, 0.3, baseq 10, mapq 20 | SNP list equals model; alt counts match |

**Review notes:**

- Very large test covering four functions (reads, pileup, coverage, variants); a failure early hides later sections. Could be split.

## `test_tracks_gff_bed_parsing_query_density_and_search`

**Result:** ✅ passed in 0.10 s · [line 427](../../tests/test_core_formats.py#L427) · **Area:** Annotation tracks (GFF3/GTF/BED) · **Checked against:** in-test reference model · **Strength:** 🟢 strong

**Why it matters:** Without it, imported annotation tracks could show features with off-by-one coordinates, wrong names or strands, crash on malformed lines, drop features across the origin or a whole type, lose the file's colours, or miss features when zooming and searching.

Parses GFF3/GTF/BED files written in the test and checks names, coordinates, strands, error handling; compares query and density against a brute-force overlap on 400 random BED features; checks text search. It also checks that a bad line among good ones is skipped and counted in Track.notes (the first one named), while a file with no valid line is still an error. Every type is counted (a whole-sequence region flagged) before the kept ones are chosen; a feature across the origin is one item on a circle and left out (counted) on a linear sequence; the colour of a line is read in every way files write it.

| # | Input | Expected output |
|---|---|---|
| 1 | GFF3 with escaped Name, GTF-style attributes, region line, color=red, ##FASTA; keep the types that do not cover the whole sequence | types include (region,1,whole) and (gene,4,not whole); colours [('red',1)], field counts (1,0); keep drops 0; names 'big,gene', 'a;b=c', '100% é', 'araC', 'bad%ZZ', 'gtfname', 'tRNA', 'onebase' in start order; 1-based 10..100 -> (9,100,+1); '.' strand 0; no region item |
| 2 | name attribute preference: Name, gene_name, gene, locus_tag, gene_id, ID, transcript_id | highest present key wins |
| 3 | start 0, non-number, too few columns, plain text, missing file | OSError |
| 4 | c: region 950..50 (Name=wraps, color=%23B91515) and region 100..200; keep region on a circle of 1000; and on a linear sequence | sequences [(c,2)]; circle: 0 dropped, query 990..1000 -> (949,1050), query 1000..1010 -> 'wraps'; linear: 1 dropped, only 'in' |
| 5 | read_color of #257215, #2a1, B91515, 0xB91515, 185,21,21, '185 21 21', rgb(185, 21, 21), darkgreen, 2, 17, 18, 300,0,0, green-ish | #257215, #22aa11, #b91515 (x5), #006400, #ff0000, #ffbfbf, None, None, None; colors_of the wrapping track -> {'#B91515': '#b91515'}, no warning |
| 6 | BED with 3..12 columns, tab and space separated | 0-based coordinates as written; name from column 4; strand from column 6; same with spaces |
| 7 | BED with track/browser/comment/blank lines and zero-length feature | 3 items; zero-length 10..10 found in 5..15 and 0..11, not in 10..20; bad lines raise OSError |
| 8 | 400 random BED features, 5 windows | query equals brute-force overlap, in start order; limit 399 -> more=True, 400 -> False; unknown sequence ([], False) |
| 9 | GFF with empty Name= | falls back to gene, then ID |
| 10 | density 100 / 7 / 40 bins | equals per-bin overlap count; 100 tidy features in 100 bins -> [1]*100; empty cases zeros/[] |
| 11 | search 'lac', 'LAC', 'LaC', 'écoli', 'ecoli', limit 3 | 12 matches case-insensitive; unicode case folding; no accent folding; limit respected; row tuple ('chr1',0,10,0,'LacZ','bed') |
| 12 | a GFF with one good line and two bad ones (end before start, start 0); a BED with one bad and one good line; a clean BED | the good line is loaded; notes == ['2 lines were skipped as not valid GFF/GTF (first: line 3, start is not a number).'], the BED's '1 line was skipped as not valid BED (first: line 1, end before start).', and [] for the clean one |
| 13 | t.search('lac', 3, 'chr3') and t.search('lac', 100, 'chr1') / 'nope' on a track of three sequences | three items of chr3 in start order; LacZ and lacY; [] (the limit counts only the sequence asked for) |

## `test_feature_locations_round_trip`

**Result:** ✅ passed in 0.01 s · [line 621](../../tests/test_core_formats.py#L621) · **Area:** GenBank feature locations · **Checked against:** round-trip/self-consistency · **Strength:** 🟢 strong

**Why it matters:** Without it, opening and saving a GenBank file could silently change or drop complex feature locations (joins, complement, fuzzy ends, between-base sites, remote references).

Adds each of 20 GenBank location forms to a document and checks the written text, spans and strand, that the written form re-parses to itself, and that a GenBank file with these locations survives save and reload. It also checks that an EMBL file cut short is still read and says so.

| # | Input | Expected output |
|---|---|---|
| 1 | 20 forms: 10..20, 5, 5..5, complement, join, order, complement(join), join(complement...), <, >, 10^11, 200^1, one-of, bond, remote J00194.1, gap(10) | written form, 0-based spans and strand as listed (e.g. 5..5 -> '5'; 10^11 -> (10,10); join with remote part keeps only local span) |
| 2 | written form added again; set_feature with the text | same location/spans/strand |
| 3 | 10..300, 0..5, 20..10, 'bad', '' | refused (OSError/ValueError), nothing added; type 'two words' refused |
| 4 | GenBank file with 10..(20.30) and join(1..5,(10.20)..(30.40)) | kept word for word after load and save; typing it in raises 'outside' |
| 5 | GenBank file with every non-gap form, saved and reloaded | same (location, spans, strand) list |
| 6 | complement(join(190..200,1..10)) on 200 bp circle | spans [(189,200),(0,10)], strand -1 |
| 7 | an EMBL file whose last record has no // line, then the same with // | the sequence is read and Document.notes says the file may be cut short; with // the notes are [] |

## `test_genome_fasta_fetch_and_search_match_slicing_and_brute_force`

**Result:** ✅ passed in 0.02 s · [line 700](../../tests/test_core_formats.py#L700) · **Area:** Indexed genome FASTA · **Checked against:** in-test reference model · **Strength:** 🟢 strong

**Why it matters:** Without it, a reference genome opened from a FASTA file could return wrong bases at line breaks or ends, and motif/sequence search on it could miss hits or report wrong strands.

Writes a multi-record 60-column FASTA and checks core.Fasta fetch against Python slicing, search against a brute-force scan on both strands, and search_motif against core.find_motif on the same string. Also the topology of a FASTA record: its header tag, else one record = circular, and the search and the document follow it.

| # | Input | Expected output |
|---|---|---|
| 1 | 3 records chr1 250, chr2 120, chr3 7 | sequences() with lengths; .fai written |
| 2 | fetch first base, last base, whole, end past length | equal slicing, end cut to length |
| 3 | ranges across line breaks (58..62, 59..61, 118..122, 10..245, ...) | equal slicing |
| 4 | chr2 last base at exact line multiple; empty/inverted/past-end ranges; chrZ | correct base; ''; OSError |
| 5 | 300 random ranges | equal slicing |
| 6 | 6 bases straddling the first line break | found at 57 forward |
| 7 | search_motif GARTTC/0, ACGT/1, NNGG/0 | equals core.find_motif on the string (a file of several records: each linear) |
| 8 | one-record FASTA '>p / TTCAAAAAGA', then the same with '[topology=linear]' in the header; search_motif GATTC (one strand) | no tag: Fasta.circular and Document.from_fasta circular, hit (8,13) across the origin; tag linear: not circular, no hit |
| 9 | A/C filler with GATTACA fwd at 30 and rc at 77 | [(30,1),(77,-1)]; rc motif swaps strands; lower case same |
| 10 | AAAA in 9 A; palindrome GAATTC | six overlapping positions; palindrome once on top strand |
| 11 | motifs AC, ACA, CACACACACA, GAT, TACAC | equal brute-force scan both strands |
| 12 | empty motif; 500 G | [] |

**Review notes:**

- Overlaps a little with the FASTA genome part of test_primer_binding.py::test_binding_agrees_with_the_documented_rules_written_in_python and the file-scope part of test_search.py.

## `test_bam_variants_leave_nothing_out_near_the_ends_and_read_the_origin_of_a_circle`

**Result:** ✅ passed in 0.04 s · [line 781](../../tests/test_core_formats.py#L781) · **Area:** BAM reads · **Checked against:** hand-computed · **Strength:** 🟢 strong

**Why it matters:** A break, a gap or a deletion near an end of the reference, or through the origin of a plasmid, could be hidden or mistaken for something else.

On a 2000 bp reference: linear, uncovered ends are gaps and a clip at the very end is a break; circular, reads over the origin (also with a missing base in the clip) are not breaks, a clip that is something else is, hard-clipped reads are reported with a note, reads that jump over the origin are one deletion through it, an uncovered start and end are one gap through it.

| # | Input | Expected output |
|---|---|---|
| 1 | linear: reads over 100..199 and three clipped at the very end | gaps (0,100), (200,1750) and a clip at 1999 |
| 2 | circular: three reads over the origin (primary + supplementary), the clip lacking a base | no break at all; read as linear, breaks at 0 and 1999 |
| 3 | circular: the same reads with a clip that is not the circle | a break at 1999 |
| 4 | circular: supplementary pieces hard-clipped | a clip at 0 saying it could not be checked |
| 5 | circular: 80 bases (1961..2000 and 1..40) missing | one 'del' of 80 bases at 1960, alt 'across the origin'; no clip and no gap |
| 6 | circular: nothing over the first and the last 100 bases | one gap of 200 bases at 1900 through the origin |

**Review notes:**

- The reads are error-free; the edit-distance tolerance is checked only with one missing base.
