# `tests/test_core_sequence.py`

[← all test files](README.md) · source: [`tests/test_core_sequence.py`](../../tests/test_core_sequence.py)

## `test_melting_temperature_matches_santalucia_and_biopython`

**Result:** ✅ passed in 0.02 s · [line 74](../../tests/test_core_sequence.py#L74) · **Area:** Primer Tm · **Checked against:** independent implementation · **Strength:** 🟢 strong

**Why it matters:** Without it, the primer melting temperatures shown everywhere in the editor (primer list, PCR design) could drift from the published SantaLucia/Primer3 values and nobody would notice.

Checks core.tm against a hand-coded SantaLucia 1998 calculation and against Biopython Tm_NN (salt correction 5) on many random primers and salt conditions, including self-complementary primers, the direction of each parameter, and tm_range for degenerate primers.

| # | Input | Expected output |
|---|---|---|
| 1 | M13 fwd, lac and T7 primers at 1 M Na, no Mg/dNTP, 250 nM | equal to the hand SantaLucia 1998 formula within 1e-6 and to Biopython (saltcorr 0) within 0.5 C; at default PCR conditions 45-70 C |
| 2 | prefixes of a random 40-mer, lengths 10..40 | Tm strictly increases with each added base |
| 3 | 20-mers AT-only / ATGC repeat / GC-only | Tm(AT) < Tm(mixed) < Tm(GC) |
| 4 | lower case, and U instead of T | same Tm as the upper-case DNA |
| 5 | lac primer at Na 50, Mg 1.5, dNTP 0.6, 50 nM, varying one parameter | more Na, more Mg, more primer raise Tm; more dNTP or less Na lowers it; dNTP >= Mg equals no Mg at all |
| 6 | 200 random non-palindromic 12-40-mers x 3 conditions (50/1.5/0.6/50, 50/0/0/250, 100/2/0.8/500) | equal to Biopython Tm_NN(DNA_NN3, saltcorr 5, dnac1=dnac2=C/2) within 0.01 C |
| 7 | ACGTNACGT | None |
| 8 | 6 palindromes (GAATTCGAATTC, ...) x 3 conditions | equal to Biopython selfcomp=True within 0.01 C, and >1 C away from the non-self-complementary formula |
| 9 | degenerate AGCGGRTAACAATYTCACACNGGA (16 expansions) | tm_range = (min, max) of the 16 expanded oligos; tm() of it is None; plain primer gives (tm, tm); ACGTX gives None |

**Review notes:**

- The abs=0.5 comparison to the hand formula is made redundant by the abs=1e-6 one on the next line.

## `test_hairpins_and_dimers`

**Result:** ✅ passed in 0.01 s · [line 147](../../tests/test_core_sequence.py#L147) · **Area:** Primer hairpins and dimers · **Checked against:** property/invariant · **Strength:** 🟡 adequate

**Why it matters:** Without it, the hairpin and primer-dimer warnings in primer checks could report nonsense (positive energies, asymmetric pairs, wrong 3'-end flag) without being caught.

Checks core.hairpin and core.dimer on designed palindromes and stems, and general sign, symmetry and ordering properties on random primers. There is no external reference for the actual dG values.

| # | Input | Expected output |
|---|---|---|
| 1 | hairpin GGATCC TTTT GGATCC | any dG < -1 and 3' dG equals it (stem ends at 3' end) |
| 2 | dimer of palindrome GAATTCGAATTC with itself | dG < -5, 3' dG equals the any-dG |
| 3 | poly-A 20-mer; 'AC'; '' | hairpin and dimer = (0.0, 0.0) |
| 4 | 300 random pairs, 8-30 bases | all dG <= 0; any-dG <= 3'-dG; dimer(a,b) == dimer(b,a) |
| 5 | lac primer with its full rc, first 12, first 8 | longer perfect duplex more stable: full < 12 < 8 < 0 |
| 6 | GCGCGCGC vs ATATATAT duplex | GC more stable |
| 7 | GCGCGC loop TTTT vs loop of 12 T; AT stem | short loop < long loop < 0; GC stem < AT stem |
| 8 | hairpin followed by 10 A | any dG < -3 but 3' dG = 0.0 |
| 9 | Na 200 vs 20 mM | hairpin and dimer more negative at higher salt |

**Review notes:**

- No absolute dG value is checked against Primer3 or a hand calculation: a scaled or shifted energy model would still pass.

**Fixed by the audit of November 2026:**

- 13 Nov 2026: gap closed. The proposed checks were added, pass on the current code and fail under the mutation that slipped through before. Was: No absolute dG is asserted, so a wrong NN parameter, a wrong loop_dg table entry, a scaled energy or a broken salt correction (0.368*(N-1)*ln(Na_eq)) pass: only signs, symmetry and orderings are checked (proven with a 0.7x scale and with salt ignored). \| The 3'-end flag is only checked on 2 hand cases plus any<=3' on random pairs; a 3' flag that is set for every hairpin ending within the last base or two is not distinguished.

## `test_gc_bins_match_an_exact_count`

**Result:** ✅ passed in 0.03 s · [line 213](../../tests/test_core_sequence.py#L213) · **Area:** GC content plot · **Checked against:** in-test reference model · **Strength:** 🟢 strong

**Why it matters:** Without it, the GC-content track could show wrong percentages, especially zoomed out where cached 1024-base blocks are used, or misbehave at sequence ends.

Checks Document.gc_bins on uniform and mixed sequences, compares block-based and base-by-base bins to an exact count on a random 100 kb sequence, and covers edge cases (end past length, empty window, more bins than bases).

| # | Input | Expected output |
|---|---|---|
| 1 | 100 G / 100 A, 4 bins; 'gcgcgcgc'; 'NNGC' | [1.0]*4, [0.0]*4, [1.0,1.0] (lower case counts), [0.5] (N in denominator) |
| 2 | 10 000 G in 1 bin; 10 001 A in 2 bins (block path) | [1.0]; [0.0, 0.0] |
| 3 | 5120 A + 5120 G, 2 bins | [0.0, 1.0] |
| 4 | random 100 kb, 4 bins of 25 000 (block sums) | within 0.02 of the exact count |
| 5 | same, 100 bins of 1000 (base by base) | exact count within 1e-6 |
| 6 | circular G10A10: window 0..1000 in 2 bins; 5..15 in 1 bin | end cut to length: [1.0, 0.0]; [0.5] |
| 7 | 'GA' in 6 bins | 6 values: 1.0, 0.0 and 4 NaN |
| 8 | empty (7,7), inverted (15,5), 0 bins | all NaN of the requested length; [] |

**Review notes:**

- The block-based path is only checked to 0.02 tolerance, so a small rounding bug at block edges would pass.

## `test_orfs_match_a_brute_force_scan`

**Result:** ✅ passed in 0.03 s · [line 299](../../tests/test_core_sequence.py#L299) · **Area:** ORF finder · **Checked against:** in-test reference model · **Strength:** 🟢 strong

**Why it matters:** Without it, the ORF track could miss open reading frames, report nested duplicates, or place reverse-strand ORFs at wrong coordinates.

Checks Document.orfs on hand-made cases (minimum length, all three stops, strands, nested/back-to-back/open ORFs) and against a brute-force six-frame scan written in the test on random sequences.

| # | Input | Expected output |
|---|---|---|
| 1 | CC ATG AAA AAA TAA CC | orfs(4) = [(2,14,1)]; orfs(5) = [] (stop counts as a codon); orfs(1) same as 4 |
| 2 | reverse complement of it | [(2,14,-1)] in top-strand coordinates |
| 3 | ccatgaaa + each stop (lower case) + cc | [(2,11,1)] for TAA, TAG, TGA |
| 4 | ATGATGAAATAA | one ORF [(0,12,1)], inner ATG not reported |
| 5 | ATGAAATAA ATGCCCTGA | [(0,9,1), (9,18,1)] |
| 6 | no ATG; ATG without stop; stops only | [] in each case |
| 7 | random 500 and 3001 bp, min codons 1/10/30 | equal to brute-force orfs_by_hand |
| 8 | random 6000 bp, min 40 codons | each ORF starts with ATG, ends with the first in-frame stop, length multiple of 3 |
| 9 | circular AAATAA + C*30 + ATGAAA (ORF over origin) | exactly one forward ORF |
| 10 | ring 'ATGAAATAA' + 'C'*30 + 'ATGAAA' (45 bp), min 3 codons; and its reverse complement | only (39, 54, +1); on the reverse complement only (36, 51, -1) |

**Review notes:**

- The ORF across the origin is only checked for count and strand, not for its start/end coordinates.
- Partly overlaps the ORF part of test_genetic_code.py::test_translation_tables_match_biopython (that one covers alternative tables and starts).

**Fixed by the audit of November 2026:**

- 13 Nov 2026: gap closed. The proposed checks were added, pass on the current code and fail under the mutation that slipped through before. Was: ORF across the origin of a circle: only `len(orfs)==1 and strand==1` is asserted, so wrong start/end (or an end not past the length) passes (proven: shifted by 3). No reverse-strand origin case and no random circular comparison: circular ORF coordinates are untested.
- 13 Nov 2026: app bug fixed in core/src/bio.rs: on a ring an ATG just after the origin, inside an ORF that started before it, was reported as a second, shorter ORF with the same stop ('ATGAAATAA' + 30 C + 'ATGAAA' gave (0,9) and (39,54)). Now one stop gives one ORF, the longest, on both strands; the test checks that ring exactly and random rings for equality with the hand scan.

## `test_restriction_cuts_and_digests_match_biopython`

**Result:** ✅ passed in 0.22 s · [line 392](../../tests/test_core_sequence.py#L392) · **Area:** Restriction digest · **Checked against:** independent implementation · **Strength:** 🟢 strong

**Why it matters:** Without it, enzyme cut sites, overhangs and digest fragment sizes shown on the map could be wrong (strand, Type IIS offsets, origin-spanning sites, degenerate sites) and cloning plans built on them would fail.

Checks enzymes.find_cuts / site_spans on hand-made sequences (overhang types, Type IIS, origin, IUPAC sites) and compares cut positions and digest fragments with Biopython's Restriction module on a real plasmid and on random circles with every common enzyme planted.

| # | Input | Expected output |
|---|---|---|
| 1 | A10 GAATTC A10 GGTACC A10 GATATC A10, linear | EcoRI (11,15) overhang AATT; KpnI (31,27) 3' overhang GTAC; EcoRV (45,45) blunt |
| 2 | every COMMON enzyme planted between 30 A on each side, both orientations | cut = declared offsets + 30, mirrored on the reverse site |
| 3 | BsaI GGTCTC fwd / GAGACC rev with 20 A pads | (27,31) and (15,19); 4-base 5' overhang |
| 4 | ATTC + A20 + GA (EcoRI split over origin) | circular: [(25,3)]; linear: not found; site_spans [(24,30,1)] |
| 5 | BsaI site_spans on fwd+rev | [(20,6,1), (66,6,-1)] |
| 6 | GGTCTC linear; GGTCTCA circular | no BsaI cut; circular wraps to (0,4) |
| 7 | BsmBI CGTCTC; AccI GTMKAC 4 matching / 3 non-matching sites; BglI GCCNNNNNGGC | (17,21); AccI (5,7) only for the allowed sites and on rc; BglI (9,6) |
| 8 | sample plasmid (circular and linear), every COMMON enzyme | cut positions equal Biopython; fragments equal and sum to length; at least 10 enzymes cut |
| 9 | 5 random 4000 bp circles with every enzyme planted, EcoRI or BsaI over the origin | origin cut found; every enzyme matches Biopython circular and linear; at least 40 enzymes cut |

**Review notes:**

- The real-plasmid part silently does nothing if the sample file is missing (no skip reported).
- Picks EcoRI/BsaI by index COMMON[14]/COMMON[9]; guarded by a name assertion, but will fail if the enzyme list is reordered.

**Fixed by the audit of November 2026:**

- 13 Nov 2026: gap closed. The proposed checks were added, pass on the current code and fail under the mutation that slipped through before. Was: Claim 'digest fragments equal Biopython / sum to length': the fragment sizes the app shows (tools.fragments) are never exercised; a bug there (off-by-one on the wrap-around fragment, a missing last linear fragment, the single-cut circle size) cannot fail this test.

## `test_find_parts_identity_coding_and_repeated_parts`

**Result:** ✅ passed in 0.00 s · [line 535](../../tests/test_core_sequence.py#L535) · **Area:** Feature detection (find_parts) · **Checked against:** hand-computed · **Strength:** 🟢 strong

**Why it matters:** Without it, automatic annotation of known parts could miss mutated or recoded copies, ignore identity thresholds, or report duplicated/overlapping hits.

Checks core.find_parts on planted 100-bp parts with spaced mutations (identity thresholds), short exact-only parts, an origin-spanning part, recoded coding parts found at protein level, and repeated/overlapping/tandem parts.

| # | Input | Expected output |
|---|---|---|
| 1 | 100 bp part in 300+300 random flanks, both strands | [(0,300,400,strand,1.0,False)] |
| 2 | 4 mutations (identity 0.96) / 6 mutations (0.94) | found at 0.96; not found by default (0.95); found at 0.94 with min_dna=0.9 |
| 3 | 1 mutation with min_dna=0.995 | not found |
| 4 | 8-mer GAATTCAC in A/C host; same with 1 mismatch at min_dna 0.5 | found at (200,208,1); mismatched short part not found (exact only) |
| 5 | 60 bp part split over the origin of a 560 bp circle | circular: (0, n-20, n+40, 1, 1.0, False); linear: [] |
| 6 | 60-codon protein recoded with synonymous codons (DNA identity < 0.9), coding=True | found at (400, 583, +1), identity 1.0, coding flag True; not found with coding=False; reverse strand found at mirrored start |
| 7 | 3 of 60 amino acids changed | identity 1 - 3/60; dropped at min_protein 0.97, kept at 0.9 |
| 8 | two overlapping different parts | both reported: (0,100,220), (1,160,220) |
| 9 | same part twice forward and once reversed | [(50,1),(220,-1),(390,1)] |
| 10 | 40 bp unit x 3 tandem | [(30,70),(70,110),(110,150)] |
