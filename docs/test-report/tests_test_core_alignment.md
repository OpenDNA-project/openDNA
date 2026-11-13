# `tests/test_core_alignment.py`

[← all test files](README.md) · source: [`tests/test_core_alignment.py`](../../tests/test_core_alignment.py)

## `test_align_single_queries_coordinates_indels_and_overhangs`

**Result:** ✅ passed in 0.01 s · [line 52](../../tests/test_core_alignment.py#L52) · **Area:** Local alignment (core.align) · **Checked against:** hand-computed · **Strength:** 🟢 strong

**Why it matters:** Without it, Local Align could report wrong coordinates, strands, indel positions or invented alignments when a user maps a read or construct onto a sequence.

Checks core.align with planted queries in random templates and exact expected hit tuples: strands, query offsets, score, indels, mismatches, max_hits, overhangs at template ends (linear and circular), split-end queries and the reporting threshold.

| # | Input | Expected output |
|---|---|---|
| 1 | 50 bp flank + template[1000:1300], and its rc | (1000,1300,+1), query part 50..350, identity 1.0; rc: strand -1, query 0..300 |
| 2 | template[1000:1300] alone | score 600 (+2 per base) |
| 3 | query with 5-base unslidable deletion, both strands | deletion list [(a,a+5)], no mismatches/insertions |
| 4 | query with GGG inserted before p | insertions [(p,3)] |
| 5 | one substituted base at 1150 | mismatch [(1150, base)], identity 399/400 |
| 6 | 80 bp unit placed fwd, fwd, rc | hits at (300,1),(680,1),(1060,-1); max_hits 2/1/0 respected; scores non-increasing |
| 7 | random 300 bp query; 'ACG'; empty template | [] each |
| 8 | query = 60 + 150 bp template + 60 | (0,150,1), query 60..210; rc strand -1 |
| 9 | query over origin of 1000 bp circle | circular (900,1100,±1); linear: one 100-base half |
| 10 | last 40 template bases + 60 foreign bases (and mirror, and rc) | one hit (960,1000) or (0,40), never past the ends; query offsets as typed |
| 11 | end40+start60 on line vs circle | line: two hits (0,60) and (960,1000); circle: one hit (960,1060) |
| 12 | ends match, middle 150 random / 10 random / 60 random | two hits / one bridged hit identity < 0.97 / two hits |
| 13 | 5, 10, 15 base end overlaps; 20 base overlap | no hit; exactly 20 bases (score 40) is reported |

**Review notes:**

- Leftover 'for _ in range(100): ... break' loop around pick_unslidable does nothing.
- The linear over-origin check only asserts the hit is 100 bases and ends <= 1000, not which half.

## `test_align_finds_every_piece_of_a_structural_variant_query`

**Result:** ✅ passed in 0.06 s · [line 299](../../tests/test_core_alignment.py#L299) · **Area:** Local alignment, structural variants · **Checked against:** hand-computed · **Strength:** 🟢 strong

**Why it matters:** Without it, Local Align could merge, drop or misplace pieces of a rearranged construct (insertion, inversion, shuffled or duplicated segments), which is what users rely on to diagnose structural variants.

Builds A-gap-B-C(inverted)-D queries for six template layouts, on lines and circles and as reverse complements, and checks that each piece is reported once with its strand and coordinates within 2 bases; plus origin-spanning, max_hits, duplicated and mismatched pieces.

| # | Input | Expected output |
|---|---|---|
| 1 | 6 layouts x linear/circular, 100 bp pieces, 200 bp spacer | one hit per piece, strands (+,+,-,+), coordinates within 2 bases, identity >= 0.95, score >= 190 |
| 2 | same query reverse-complemented | same pieces, strands flipped, query coordinates mirrored |
| 3 | D piece = template[1950:2050] on a 2000 bp circle | circular: one hit ending >= 2050; linear: D split into (1950,2000) and (0,50) |
| 4 | A,B,D collinear layout, max_hits 2/1/0 | 2 / 1 / 0 hits; no hit spans the 200 bp gap; 4 hits best first |
| 5 | same 100 bp twice in the query | two hits, both at template 700..800 |
| 6 | exact 100 bp piece + piece with a substitution and a 3-base deletion | two hits; second covers 300..420, has mismatch and deletion, identity < 1 |

**Review notes:**

- Overlaps somewhat with the 'ends match, middle does not' section of test_align_single_queries_coordinates_indels_and_overhangs.

## `test_local_align_query_cleaning_and_junction_texts`

**Result:** ✅ passed in 0.03 s · [line 356](../../tests/test_core_alignment.py#L356) · **Area:** Local Align query input and report text · **Checked against:** hand-computed · **Strength:** 🟡 adequate

**Why it matters:** Without it, pasting a FASTA, GenBank ORIGIN or RNA query could be rejected or mangled, and the junction descriptions (insertion, deletion, inversion) in the Local Align report could be wrong.

Checks local_align.clean_query on accepted formats and on errors with their messages, and junction_text / summarize on hand-made hit pairs.

| # | Input | Expected output |
|---|---|---|
| 1 | >seq1 description + lower case with spaces | ACGTACGTGGCC |
| 2 | GenBank ORIGIN-style lines with numbers; ACGU UGCA; AC-GT--AC; IUPAC | ACGTACGTGGCCTT; ACGTTGCA; ACGTAC; kept as is |
| 3 | two FASTA records; CRLF and leading blank lines | only the first record; ACGTACGT |
| 4 | '', blanks, header only, '123 456' | QueryError mentioning empty / only a FASTA header / No DNA |
| 5 | protein MKVLAAGIVGLLLAQPSEEF; 'ACGT hello 12' | QueryError 'protein'; 'Not DNA' |
| 6 | MAX_QUERY+1 bases; exactly MAX_QUERY | error 'too long ... 200,000'; accepted |
| 7 | adjacent hits; 200 bp query gap; 50 bp template skip; template jumps back; strand change | 'contiguous'; 'insertion in the query: 200 bp unmatched'; 'deletion from the query: 50 bp...'; 'template jumps back'; 'inversion breakpoint' |
| 8 | two reverse pieces running backwards; circle 4900..5000 then 50..150 | 'contiguous'; 'deletion from the query: 50 bp' |
| 9 | summarize 1 hit; 3 hits with gap and inversion; [] | '1 piece: 1 forward'; '3 pieces: 2 forward, 1 inverted ... 200 bp unmatched (#1→#2)' with both junction labels; '' |

**Review notes:**

- Junction texts are checked with startswith only, so the rest of each message (sizes, positions) is not verified.
- Duplication wording is mentioned in the docstring but no duplication junction is tested.

## `test_folder_search_finds_a_query_in_its_file_record_place_and_strand`

**Result:** ✅ passed in 1.30 s · [line 399](../../tests/test_core_alignment.py#L399) · **Area:** Folder Search · **Checked against:** hand-computed · **Strength:** 🟡 adequate

**Why it matters:** Without it, Folder Search could fail to find a sequence in the user's files, report the wrong file, record, position or strand, or crash on an unreadable file.

Searches the demo folder and a temporary folder with search_folder and checks the hit location, strand, the skipped-file count, the multi-record name, the tab-separated hits text and cancellation.

| # | Input | Expected output |
|---|---|---|
| 1 | pFG051_Chlor.gb bases 1000..1060 searched in demo/ | hit in that file at (1000,1060,+1), identity 1.0; nothing skipped; >= 5 files |
| 2 | its reverse complement | (1000,1060,-1) |
| 3 | random 60 bp (seed 4171) | no hits |
| 4 | many.fa with records first/second + broken.gb, query second[500:560] | files 2, skipped 1; one hit ('second',500,560,+1); text line 'many.fa\tsecond\t501-560\tforward\t100.0%\t60' |
| 5 | cancelled=lambda: True | no hits |

**Review notes:**

- The cancellation check is weak: the query ACGTACGTACGT would almost certainly not be found in that folder anyway, so a search that ignores cancelled() would still pass.
- 'files >= 5' depends on the contents of demo/.
