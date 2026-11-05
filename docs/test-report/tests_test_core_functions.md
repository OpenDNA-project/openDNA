# `tests/test_core_functions.py`

[← all test files](README.md) · source: [`tests/test_core_functions.py`](../../tests/test_core_functions.py)

## `test_melting_temperature`

**Result:** ✅ passed in 0.00 s · [line 72](../../tests/test_core_functions.py#L72) · **Area:** Primer Tm · **Checked against:** independent implementation · **Strength:** 🟢 strong

Checks that core.tm gives the textbook nearest-neighbour melting temperature and responds correctly to length, GC, case, U and the buffer conditions a user types in. Primer design depends on this number.

| # | Input | Expected output |
|---|---|---|
| 1 | M13 fwd 'GTAAAACGACGGCCAGT', lac 'AGCGGATAACAATTTCACACAGGA', T7 'TAATACGACTCACTATAGGG'; na_mm=1000, mg=0, dntp=0, primer_nm=250 | equals a hand-coded SantaLucia 1998 formula within 1e-6 C, and Biopython Tm_NN(DNA_NN3, saltcorr=0) within 0.5 C |
| 2 | same primers at default PCR conditions | 45 < Tm < 70 |
| 3 | prefixes of length 10..40 of a random 40-mer (seed 3) | Tm strictly increases with each added base |
| 4 | AT-only, mixed ATGC, GC-only 20-mers | Tm(AT) < Tm(mixed) < Tm(GC) |
| 5 | 'acgtacgtacgtacgt' vs upper case; 'ACGUACGUACGUACGU' vs T | identical Tm |
| 6 | lac primer, base 50 mM Na, 1.5 mM Mg, 0.6 mM dNTP, 50 nM; vary Na 100/20, Mg 3.0, dNTP 1.0, primer 500 nM | more Na/Mg/primer raises Tm, less Na or more dNTP lowers it |
| 7 | mg_mm=0.5, dntp_mm=0.6 vs mg=0, dntp=0 | equal Tm (no free Mg2+) |

**Review notes:**

- Line 79 (abs=0.5) is redundant with line 80 (abs=1e-6); cosmetic only.

## `test_hairpins_and_dimers`

**Result:** ✅ passed in 0.01 s · [line 112](../../tests/test_core_functions.py#L112) · **Area:** Primer hairpins and dimers · **Checked against:** property/invariant · **Strength:** 🟡 adequate

Checks that hairpin and dimer free energies find palindromic stems, return 0 when nothing can pair, are symmetric, and rank stability sensibly. These values drive primer warnings.

| # | Input | Expected output |
|---|---|---|
| 1 | hairpin('GGATCCTTTTGGATCC') | any dG < -1 and 3' dG == any dG |
| 2 | dimer('GAATTCGAATTC', same) | any dG < -5 and 3' dG == any dG |
| 3 | poly-A 20-mer hairpin/dimer; hairpin('AC'), hairpin(''), dimer('', 'ACGT') | (0.0, 0.0) |
| 4 | 300 random pairs, lengths 8-30 (seed 5) | both dGs <= 0, any <= 3', dimer(a,b) == dimer(b,a) |
| 5 | lac primer vs its rc at lengths 24/12/8; GCGCGCGC vs ATATATAT duplexes | longer and GC-richer duplexes are more negative |
| 6 | GCGCGC stem with loop 4 vs 12 T; GC vs AT stem | short loop < long loop < 0; GC stem < AT stem |
| 7 | hairpin('GCGCGCTTTTGCGCGC' + 10 A) | any < -3 and 3' dG == 0.0 |
| 8 | na_mm 200 vs 20 for hairpin and dimer | more salt gives more negative dG |

**Review notes:**

- No absolute dG value is compared to a reference (e.g. Primer3 or a hand calculation); only signs and orderings, so a consistently mis-scaled energy model would pass.

## `test_gc_bins_uniform_mixed_and_block_counts`

**Result:** ✅ passed in 0.03 s · [line 155](../../tests/test_core_functions.py#L155) · **Area:** GC content plot · **Checked against:** in-test reference model · **Strength:** 🟢 strong

Checks the GC-fraction bins behind the GC plot, both when reading bases directly and when zoomed out on cached 1024-base blocks.

| # | Input | Expected output |
|---|---|---|
| 1 | 'G'*100 / 'A'*100 in 4 bins; 'gcgcgcgc' in 2; 'NNGC' in 1 | [1.0]*4, [0.0]*4, [1.0,1.0], [0.5] |
| 2 | 'G'*10000 in 1 bin; 'A'*10001 in 2 bins (block path) | [1.0]; [0.0, 0.0] |
| 3 | 'A'*5120 + 'G'*5120 in 2 bins | [0.0, 1.0] |
| 4 | random 100 kb (seed 9), 4 bins of 25 000 (block path) | exact per-bin GC within 0.02 |
| 5 | same sequence, 100 bins of 1000 (base path) | exact per-bin GC within 1e-6 |

## `test_gc_bins_edges_end_past_the_sequence_empty_bins_and_more_bins_than_bases`

**Result:** ✅ passed · [line 180](../../tests/test_core_functions.py#L180) · **Area:** GC content plot · **Checked against:** hand-computed · **Strength:** 🟢 strong

Checks the GC bins at the edges: ends past the sequence, more bins than bases, empty or inverted windows and zero bins.

| # | Input | Expected output |
|---|---|---|
| 1 | circular 'G'*10+'A'*10, gc_bins(0,1000,2) and (0,20,2) | both [1.0, 0.0] (no wrap) |
| 2 | same, gc_bins(5,15,1) | [0.5] |
| 3 | 'GA', 6 bins | 6 values: 4 NaN, the others [1.0, 0.0] |
| 4 | gc_bins(7,7,3); gc_bins(15,5,2) | all NaN (3 values for the first) |
| 5 | gc_bins(0,20,0) | [] |

## `test_orfs_lengths_strands_and_nesting`

**Result:** ✅ passed · [line 218](../../tests/test_core_functions.py#L218) · **Area:** ORF finder · **Checked against:** hand-computed · **Strength:** 🟢 strong

Checks ORF coordinates, minimum length (counting the stop), both strands, all three stops and the nested / back-to-back / no-start / no-stop cases.

| # | Input | Expected output |
|---|---|---|
| 1 | 'CCATGAAAAAATAACC', min 4 / 5 / 1 codons | [(2,14,1)] / [] / [(2,14,1)] |
| 2 | reverse complement of the same | [(2,14,-1)] |
| 3 | 'ccatgaaa'+stop+'cc' for TAA/TAG/TGA, lower case | [(2,11,1)] each |
| 4 | 'ATGATGAAATAA' | [(0,12,1)] only (inner ATG not reported) |
| 5 | 'ATGAAATAA'+'ATGCCCTGA' | [(0,9,1),(9,18,1)] |
| 6 | 'AAA'*30+'TAA'; 'ATG'+'AAA'*30; 'TAATAGTGA' | [] each |

## `test_orfs_agree_with_a_brute_force_scan_and_wrap_the_origin`

**Result:** ✅ passed in 0.01 s · [line 250](../../tests/test_core_functions.py#L250) · **Area:** ORF finder · **Checked against:** in-test reference model · **Strength:** 🟡 adequate

Checks the ORF finder against a brute-force six-frame scan on random sequences, and that an ORF across the origin of a plasmid is found.

| # | Input | Expected output |
|---|---|---|
| 1 | random 500 and 3001 bp (seed 21), min_codons 1/10/30 | equals orfs_by_hand (six-frame scan written in the test) |
| 2 | random 6000 bp, min 40 codons | non-empty; each ORF starts ATG, ends in a stop, no in-frame stop inside, length % 3 == 0 |
| 3 | circular 'AAATAA'+'C'*30+'ATGAAA', min 4 | exactly 1 ORF, strand +1 |

**Review notes:**

- Line 271: the origin-crossing ORF is only checked for count and strand; its coordinates (36, 48, 1) are not asserted, so a wrong start/end across the origin would pass. The brute-force model (line 200) has no circular mode, so wrapping is never compared with it.

## `test_cuts_overhangs_strands_origin_and_degenerate_sites`

**Result:** ✅ passed in 0.00 s · [line 277](../../tests/test_core_functions.py#L277) · **Area:** Restriction sites · **Checked against:** hand-computed · **Strength:** 🟢 strong

Checks restriction cut positions and overhangs for 5', 3', blunt and type IIS enzymes, on both strands, across the origin and with IUPAC-degenerate sites.

| # | Input | Expected output |
|---|---|---|
| 1 | A-padded EcoRI/KpnI/EcoRV sites at 10, 26, 42 | EcoRI [(11,15)] overhang AATT; KpnI [(31,27)] GTAC; EcoRV [(45,45)] blunt |
| 2 | each COMMON enzyme planted at 30 in 'A'*30 padding, both orientations | (30+top, 30+bottom); reverse site mirrored (30+L-bottom, 30+L-top) |
| 3 | BsaI GGTCTC in 'A'*20 padding, forward / reverse | [(27,31)] / [(15,19)] |
| 4 | circular 'ATTC'+'A'*20+'GA' (EcoRI over origin) | [(25,3)] on circle, absent on line; site_spans [(24,30,1)] |
| 5 | 'GGTCTC' linear / 'GGTCTCA' circular with BsaI | no cut / [(0,4)] |
| 6 | BsmBI 'CGTCTC' in T padding | [(17,21)] |
| 7 | AccI GTMKAC with GTAGAC/GTCGAC/GTATAC/GTCTAC vs GTGGAC/GTTGAC/GTAGCC | {'AccI': [(5,7)]} vs {} |
| 8 | BglI GCCNNNNNGGC cuts (7,4) on 'TTGCCATCGTGGCTT' | {'BglI': [(9,6)]} |

## `test_digest_on_the_real_plasmid_equals_biopython_restriction`

**Result:** ✅ passed in 0.10 s · [line 357](../../tests/test_core_functions.py#L357) · **Area:** Restriction digest · **Checked against:** independent implementation · **Strength:** 🟢 strong

Checks that the cut positions and fragment sizes on a real plasmid match Biopython.Restriction for every common enzyme, circular and linear.

| # | Input | Expected output |
|---|---|---|
| 1 | demo/pFG051_Chlor.dna (circular, 3780 bp), each of the 53 COMMON enzymes | top-strand cut positions equal Biopython's (1-based minus 1); fragments equal and sum to length; at least 10 enzymes cut |
| 2 | same sequence as linear, first 20 enzymes | cut positions equal Biopython linear=True |

**Review notes:**

- Line 358: runtime pytest.skip if the demo .dna file is missing (it is tracked in git, so not expected in CI). Only the top-strand cut is compared; bottom-strand positions/overhangs are not.

## `test_digest_with_every_enzyme_planted_in_a_random_circle_equals_biopython`

**Result:** ✅ passed in 0.07 s · [line 379](../../tests/test_core_functions.py#L379) · **Area:** Restriction digest · **Checked against:** independent implementation · **Strength:** 🟢 strong

Checks every common enzyme planted in random circles, including sites across the origin, against Biopython on circular and linear topology.

| # | Input | Expected output |
|---|---|---|
| 1 | 5 random 4000 bp circles (seed 33), every COMMON site planted in random orientation; EcoRI (k=3) and BsaI (k=2) planted across the origin | for each enzyme, circular and linear: cut positions equal Biopython; fragments sum to 4000 |
| 2 | same | >= 40 enzymes cut per circle; EcoRI and BsaI both cut |

**Fixed by the audit of November 2026:**

- EcoRI and BsaI now take turns crossing the origin (circles 0, 2, 4: EcoRI; 1, 3: BsaI) instead of BsaI overwriting EcoRI, and the cut of the site across the origin is asserted at its exact position.

## `test_find_parts_identity_thresholds_short_parts_and_the_origin`

**Result:** ✅ passed in 0.00 s · [line 418](../../tests/test_core_functions.py#L418) · **Area:** Feature detection (find_parts) · **Checked against:** hand-computed · **Strength:** 🟢 strong

Checks that known parts are found on both strands with the right identity and thresholds, that short parts need an exact match, and that a part across the origin is found on a circle.

| # | Input | Expected output |
|---|---|---|
| 1 | random 100 bp part between 300 bp flanks (seed 41), as is and rc | [(0,300,400,strand,1.0,False)] |
| 2 | part with 4 spaced mismatches | one hit 300..400, identity 0.96 |
| 3 | 6 mismatches; then min_dna=0.9 | []; then one hit identity 0.94 |
| 4 | 1 mismatch, min_dna=0.995 | [] |
| 5 | 8-mer 'GAATTCAC' in A/C host; then 1 mismatch with min_dna=0.5 | span (200,208,1); [] (short parts exact only) |
| 6 | 60 bp part split over the origin of a 560 bp circle | [(0, n-20, n+40, 1, 1.0, False)] on circle, [] on line |

## `test_find_parts_coding_and_repeated_parts`

**Result:** ✅ passed in 0.00 s · [line 452](../../tests/test_core_functions.py#L452) · **Area:** Feature detection (find_parts) · **Checked against:** hand-computed · **Strength:** 🟢 strong

Checks protein-level matching of coding parts (synonymous recoding, min_protein) and how overlapping, repeated and tandem parts are reported.

| # | Input | Expected output |
|---|---|---|
| 1 | 60 aa protein (seed 43), recoded with synonymous codons (DNA identity < 0.9), in AT host; coding=True | one hit (400, 400+183, 1, 1.0, True) |
| 2 | same with coding=False | [] |
| 3 | rc of host, coding=True | strand -1, identity 1.0, start = len(host)-400-183 |
| 4 | 3 of 60 amino acids changed | one hit, identity 1-3/60; [] at min_protein=0.97; found at 0.9 |
| 5 | parts a (120 bp) and b[:60]=a[60:] in one host | both reported: (0,100,220), (1,160,220) |
| 6 | a, rc(a), a separated by 50 bp | [(50,1),(220,-1),(390,1)] |
| 7 | 40 bp unit x3 tandem | spans [(30,70),(70,110),(110,150)] |

## `test_translation_of_all_codons_incomplete_codons_and_n`

**Result:** ✅ passed in 0.00 s · [line 512](../../tests/test_core_functions.py#L512) · **Area:** Translation / protein search · **Checked against:** independent implementation · **Strength:** 🟢 strong

Checks the standard genetic code for all 64 codons on both strands via find_protein, plus incomplete last codons and N-containing codons.

| # | Input | Expected output |
|---|---|---|
| 1 | each of 64 codons; its rc | (0,3,1,0,1) and (0,3,-1,0,1) in find_protein(codon, Biopython aa); a different aa not found |
| 2 | all 64 codons concatenated (and lower case) | find_protein == [(0,192,1,0,1)] |
| 3 | 'ATGGCCAAAGG': 'MAK', 'MAKG', 'WPK' | [(0,9,1,0,1)]; [] forward; (1,10,1,0,2) present |
| 4 | 'ATGNNNGCC': 'MXA', 'MAA', 'MAA' max_mismatches=1 | [(0,9,1,0,1)]; []; [(0,9,1,1,1)] |

## `test_align_strands_indels_hits_and_query_longer_than_the_template`

**Result:** ✅ passed in 0.01 s · [line 563](../../tests/test_core_functions.py#L563) · **Area:** Local alignment · **Checked against:** hand-computed · **Strength:** 🟢 strong

Checks local alignment coordinates, strand, score, exact indel and mismatch positions, max_hits, no-hit cases, queries longer than the template and hits across the origin.

| # | Input | Expected output |
|---|---|---|
| 1 | 3000 bp template (seed 51), query = 50 bp flank + t[1000:1300]; and its rc | (1000,1300,+1, q 50..350, id 1.0); rc: strand -1, q 0..300 |
| 2 | t[1000:1300] alone | score 600 (+2 per base) |
| 3 | window 1000..1400 with an unslidable 5 bp deletion (both strands); 3 bp 'GGG' insertion; one mismatch at 1150 | deletions [(a,a+5)]; insertions [(p,3)]; mismatches [(1150, base)], identity 399/400 |
| 4 | 80 bp unit at 300, 680 and rc at 1060 | hits (300,1),(680,1),(1060,-1); max_hits 2/1/0 give 2/1/0; scores non-increasing |
| 5 | random 300 bp query; 'ACG'; empty template | [] each |
| 6 | query = 60 + 150 bp template + 60 | (0,150,1), q 60..210; rc strand -1 |
| 7 | circular 1000 bp, q = t[900:]+t[:100] and rc | (900,1100,+/-1), id 1.0; on a line first hit is 100 bp ending <= 1000 |
| 8 | circular, q = t[:80] | exactly one hit |

**Review notes:**

- Lines 594-596: a 'for _ in range(100): ... break' loop that runs once; dead code. Line 638 comment says 'only the longer half is found' but both halves are 100 bp; the assertion is loose (any 100 bp hit ending <= 1000).

## `test_align_overhangs_and_partial_matches`

**Result:** ✅ passed in 0.00 s · [line 646](../../tests/test_core_functions.py#L646) · **Area:** Local alignment · **Checked against:** hand-computed · **Strength:** 🟢 strong

Checks partial alignments: queries hanging off a linear end, wrapping on a circle, matching only at both ends, and overlaps just below/at the reporting threshold.

| # | Input | Expected output |
|---|---|---|
| 1 | 1000 bp t (seed 61); q = t[-40:] + 60 bp random; and rc | one hit (960,1000,+1,1.0) q 0..40; rc: strand -1, q 60..100 |
| 2 | q = 60 bp random + t[:40]; and rc | (0,40,+1) q 60..100; rc q 0..40 |
| 3 | linear, q = t[-40:]+t[:60] | two hits [(0,60,40,100),(960,1000,0,40)] |
| 4 | circular, same q and rc | one hit (960,1060,strand,1.0), q 0..100 |
| 5 | circular, t[-40:]+random overhang | start 960, q from 0, 40 <= q_end < 60, end - 1000 <= 12 |
| 6 | t[300:340] + 150 random + t[490:550]; and rc; t[300:340]+150+t[600:660] | 2 hits: (490,550) q 190..250 first, then 300..340-342 |
| 7 | t[300:340]+10 random+t[350:410]; with 60 random + t[400:460] | one bridged hit (300,410) q 0..110, identity < 0.97 with mismatches; two pieces |
| 8 | overlaps of 5/10/15 bases at either end (linear, rc, circular) | [] each |
| 9 | 20 base overlap at start / end | [(0,20,40,60,80)] / [(980,1000,40,0,20)] |

## `test_bam_hand_made_reads_pileup_and_coverage`

**Result:** ✅ passed in 0.00 s · [line 913](../../tests/test_core_functions.py#L913) · **Area:** BAM reads / pileup / coverage · **Checked against:** in-test reference model · **Strength:** 🟢 strong

Checks BAM reading on one read of every kind (soft/hard clips, indels, reverse, ignored flags, pairs), plus pileup and coverage against a hand count. Users rely on this to inspect sequencing data.

| # | Input | Expected output |
|---|---|---|
| 1 | BAM written in-test (chrA 2000, chrB 500), 14 reads | sequences [('chrA',2000),('chrB',500)]; .bai created; reads('chrZ') raises OSError |
| 2 | reads over chrA 0..2000 | unmapped/duplicate/qcfail excluded; perfect (10,60,+1, no diffs); reverse (30,80,-1) mismatch at 45; softclip (100,145); deletion (200,250) del [(220,225)]; insertion (300,347) ins [(320,3)]; hardclip (700,740) |
| 3 | pair1/pair2/perfect mate fields | (100,True,0), (1150,False,200), (None,False,0) |
| 4 | windows 100..130, 0..200 with limit 2/4, chrB; lower-case reference | ['softclip']; ['perfect','reverse'] more=True; more=False at 4; ['chrB_read']; same mismatches |
| 5 | pileup 100..150, 215..230, 315..325, 0..2000, 1990..3000 | 45 x depth 1 then 0; deleted at 220..225; insertion before 320; whole equals model_pileup; length 10 |
| 6 | pileup 500..550 with min_baseq 20 / min_mapq 20 / 41 / 61 | 100 -> 99 / 50 bases, equal to model; 0 and 0 |
| 7 | coverage narrow bins (0,200,20),(95,155,6),(200,260,60); wide (0,1024,2),(0,2048,4); empty / 0 bins | mean depth per bin within 1e-6; has_depth false then true; [0,0,0] and [] |

**Review notes:**

- Line 988 asserts has_depth is False on the module-scoped hand_bam fixture: it holds only because no earlier test asks that BAM for wide coverage (test order dependence, currently harmless).

## `test_bam_random_reads_and_variants`

**Result:** ✅ passed in 0.06 s · [line 1035](../../tests/test_core_functions.py#L1035) · **Area:** BAM pileup / variant calling · **Checked against:** in-test reference model · **Strength:** 🟢 strong

Checks pileup and coverage on 300 random reads against a Python model, and variant calling (SNP thresholds, reference base, deletions, insertions) on designed and random reads.

| # | Input | Expected output |
|---|---|---|
| 1 | random_bam: ~300 reads (seed 62) with clips/indels/errors/flags; filters (0,0),(20,0),(0,30),(30,60),(41,0) | pileup 0..2000 and 700..913 equal model_pileup; coverage 100 bins and 2 wide bins equal model depth |
| 2 | snp_bam: 10 reads, SNP at 120 in 5, at 125 in 1, at 130 in 3; min_depth 5, min_fraction 0.2 | [(120,'snp',ref,alt,5,10),(130,'snp',ref,alt,3,10)]; nothing at 125 |
| 3 | min_fraction 0.5/0.6/0.4/0.05; min_depth 10/11 | {120}/{}/{120}/{120,125,130}; {120}/{} |
| 4 | min_baseq=41; min_mapq=61 | [] |
| 5 | min_depth 1, fraction 0.1 | SNP list equals model_snps |
| 6 | reference with alt at 120; lower-case ref; N at 120; chrZ | (120,'snp',alt,ref,5,10); identical; no call; OSError |
| 7 | hand_bam, min_depth 1, fraction 0.5 | del (220,'del',ref[220:225],'',1,1); ins (320,'ins','','GGG',1,1); one SNP at 45 with (1,2); no call 95..99 |
| 8 | random_bam, variants(5, 0.3, 10, 20) | SNPs equal model_snps; alt count equals model pileup |

**Review notes:**

- Line 1047: 'any(c[5] ...)' only checks some deleted base exists; the comment's claim that deleted bases are in pileup depth but not in coverage is not compared.

## `test_gff3_attributes_and_malformed_lines`

**Result:** ✅ passed in 0.00 s · [line 1141](../../tests/test_core_functions.py#L1141) · **Area:** Track import: GFF3/GTF · **Checked against:** hand-computed · **Strength:** 🟢 strong

Checks GFF3 parsing: percent-escaped attributes, GTF quoting, name preference order, 1-based to 0-based coordinates, and that malformed lines raise errors.

| # | Input | Expected output |
|---|---|---|
| 1 | GFF3 with region, genes with Name=big%2Cgene, a%3Bb%3Dc, 100%25 %C3%A9, CDS gene=araC, bad%ZZ, GTF-quoted mRNA, tRNA with no name, one-base gene, ##FASTA | names in start order ['big,gene','a;b=c','100% e-acute','araC','bad%ZZ','gtfname','tRNA','onebase']; region skipped |
| 2 | coordinates / strands | big,gene (9,100,1); a;b=c (119,180,-1); '.' strand 0; araC (309,320,1) CDS; gtfname (349,360,1); onebase (399,400) |
| 3 | attributes Name, gene_name, gene, locus_tag, gene_id, ID, transcript_id, dropping the first k | name is v_<first remaining key> |
| 4 | end < start; start 0; 'ten'; 8 columns; 'just some text'; missing file | OSError each |

## `test_bed_columns_headers_and_errors`

**Result:** ✅ passed in 0.01 s · [line 1191](../../tests/test_core_functions.py#L1191) · **Area:** Track import: BED · **Checked against:** hand-computed · **Strength:** 🟢 strong

Checks BED files with 3 to 12 columns, tab or space separated, header/comment lines, zero-length features and bad lines.

| # | Input | Expected output |
|---|---|---|
| 1 | chr1 100 200 geneX 500 - ... with 3..12 columns; chr2 line space-separated | (100,200,'bed'); name 'geneX' if >= 4 columns else ''; strand -1 if >= 6 else 0; chr2 (5,15) same; sequences [('chr1',1),('chr2',1)] |
| 2 | track/browser/#/blank lines, then 'chr1 10 10 empty', '20 30 n .', '40 50 p +' | [(10,10,0,'empty'),(20,30,0,'n'),(40,50,1,'p')] |
| 3 | query windows 5..15, 10..20, 0..11 | ['empty'], [], ['empty'] |
| 4 | 'not-a-number', end < start, 2 columns | OSError |

## `test_track_query_density_and_search`

**Result:** ✅ passed in 0.01 s · [line 1234](../../tests/test_core_functions.py#L1234) · **Area:** Track query / density / search · **Checked against:** in-test reference model · **Strength:** 🟢 strong

Checks track window queries against a brute-force overlap, the limit/more flag, density bins and case-insensitive name search.

| # | Input | Expected output |
|---|---|---|
| 1 | 400 random BED items (seed 71) on chr1; windows (0,10000),(1000,1100),(4321,4322),(8900,9700),(20000,30000) | names equal brute-force overlap, in start order, more=False |
| 2 | limit 399 / 400 / 1; chr2; default limit | 399 more=True; 400 more=False; more=True; ([],False); 400 items |
| 3 | density (0,10000,100),(0,10000,7),(2000,6000,40) | per-bin overlap count equals brute force |
| 4 | 100 items of 5 bp every 10 bp, 100 bins; empty windows; unknown seq | [1]*100; [0,0,0,0]; []; [0,0,0]; [0,0] |
| 5 | search 'lac'/'LAC'/'LaC', 'ecoli' with accent, limits | 12 names; 'ECOLI-gene' for e-acute only; 3 with limit 3; row ('chr1',0,10,0,'LacZ','bed'); 'zzz' [] |

## `test_locations_round_trip_and_refusals`

**Result:** ✅ passed in 0.00 s · [line 1317](../../tests/test_core_functions.py#L1317) · **Area:** Feature locations · **Checked against:** hand-computed · **Strength:** 🟢 strong

Checks that every GenBank location form parses to the right spans/strand, is written back unchanged, and that invalid locations are refused without leaving a feature.

| # | Input | Expected output |
|---|---|---|
| 1 | 20 forms on a 200 bp circle: '10..20', '5', 'complement(...)', join/order, '<10..20', '10^11', '200^1', one-of, bond, remote 'J00194.1:100..202', gap(10) | location text as listed, spans e.g. [(9,20)], [(0,5),(9,20)], [(10,10)], [(0,0)]; strand +/-1; re-parsing gives the same; set_feature accepts it |
| 2 | '10..300', '0..5', '20..10', 'bad', '' | OSError or ValueError, no feature left |
| 3 | feature kind 'two words' | OSError |

## `test_locations_in_genbank_files`

**Result:** ✅ passed in 0.00 s · [line 1342](../../tests/test_core_functions.py#L1342) · **Area:** Feature locations / GenBank I/O · **Checked against:** round-trip/self-consistency · **Strength:** 🟡 adequate

Checks that unusual locations in GenBank files survive load and save word for word, and that an origin-spanning complement join gets the right spans.

| # | Input | Expected output |
|---|---|---|
| 1 | GenBank with '10..(20.30)' and 'join(1..5,(10.20)..(30.40))' | kept verbatim after load and after write_genbank; typing '10..(20.30)' raises OSError matching 'outside' |
| 2 | file with all LOCATION_FORMS except the gap one | 19 features; (location, spans, strand) identical after save and reload |
| 3 | add 'complement(join(190..200,1..10))' | spans [(189,200),(0,10)], strand -1 |

**Review notes:**

- Line 1360 silently drops the 'gap(10)' form from the file round trip with no comment on why; that form is only covered by the in-memory test.

## `test_genome_fasta_fetch_and_motif_search`

**Result:** ✅ passed in 0.00 s · [line 1394](../../tests/test_core_functions.py#L1394) · **Area:** Indexed genome FASTA · **Checked against:** in-test reference model · **Strength:** 🟡 adequate

Checks that an indexed multi-record FASTA returns the right bases across line breaks and at the ends, and that motif search on the file works.

| # | Input | Expected output |
|---|---|---|
| 1 | FASTA (seed 81): chr1 250 bp, chr2 120 bp, chr3 7 bp, 60 columns | sequences listed with lengths; .fai written |
| 2 | fetch first/last base, whole, end past length, last 3 | equal Python slices |
| 3 | chr1 ranges across line breaks (58,62),(59,61),(60,61),(118,122),(179,181),(10,245),...; chr2 119..120 | equal Python slices |
| 4 | (5,5), (9,5), (250,260); 'chrZ' | ''; OSError |
| 5 | 300 random ranges (seed 82) | equal Python slices |
| 6 | search s[57:63]; search_motif 'GARTTC'/0, 'ACGT'/1, 'NNGG'/0 | hit at 57 strand +1; search_motif equals core.find_motif on the string |

**Review notes:**

- Lines 1429-1432: Fasta.search_motif is implemented as fetch + search::find_motif (core/src/fasta.rs:113-115), so comparing it to core.find_motif reuses the code under test; it only checks the plumbing, and nothing asserts the hit lists are non-empty. Fetch checks are strong.

## `test_genome_fasta_search_both_strands_and_overlapping_hits`

**Result:** ✅ passed in 0.00 s · [line 1438](../../tests/test_core_functions.py#L1438) · **Area:** Indexed genome FASTA search · **Checked against:** in-test reference model · **Strength:** 🟢 strong

Checks exact search on an indexed FASTA: both strands, case, overlapping hits, palindromes reported once, and agreement with a brute-force scan.

| # | Input | Expected output |
|---|---|---|
| 1 | A/C filler with 'GATTACA' at 30 and its rc at 77 | search == [(30,1),(77,-1)]; rc motif gives [(30,-1),(77,1)]; lower case same |
| 2 | 'AAAA' in a run of nine A | six consecutive positions |
| 3 | 'GAATTC' palindrome | [(e,1)] once |
| 4 | motifs 'AC','ACA','CACACACACA','GAT','TACAC' | equal a brute-force scan of both strands |
| 5 | '' and 'G'*500 | [] |
