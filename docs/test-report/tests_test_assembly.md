# `tests/test_assembly.py`

[← all test files](README.md) · source: [`tests/test_assembly.py`](../../tests/test_assembly.py)

## `test_fragment_and_gibson_products`

**Result:** ✅ passed in 0.06 s · [line 100](../../tests/test_assembly.py#L100) · **Area:** Gibson assembly · **Checked against:** hand-computed · **Strength:** 🟢 strong

**Why it matters:** Without it, Gibson designs could hand the user primers or a predicted product whose junctions do not actually exist in the assembled DNA, and silently accept ends that share no homology.

Checks that a fragment behaves as a PCR product (tails, core, flipping swaps the primers), that Gibson analysis reports non-fitting ends and then fixes them with auto-designed primers, that homology the fragments already share is recognised and judged for length and uniqueness, and that every junction really sits at the position it claims in the product.

| # | Input | Expected output |
|---|---|---|
| 1 | 600 bp template, fwd primer 'GGGGCC'+bases 100-122, rev primer 'TTTAAA'+rc(300-322), Frag(100,322) | seq = GGGGCC + template[100:322] + rc(TTTAAA); tails (6,6); core_seq = template[100:322] |
| 2 | the same fragment flipped | seq becomes its reverse complement, tails unchanged, flipped True, fwd/rev primer names swap (F<->R); flipping again restores |
| 3 | three random fragments A(300) B(400) C(250), analyse circular Gibson | all three joints kind 'none', design not ok, check 'Circle closes' = bad, first problem mentions 'no homology' |
| 4 | auto_primers on the same three fragments | returns 3; re-analysis ok; joints (A,B),(B,C),(C,A); every overlap 20-40 nt with Tm >= 52; 6 primer rows; product 950 bp; same product (up to rotation) as design_gibson on plain copies; 'Circle closes' = ok |
| 5 | a(200 bp + 30 shared) and b(same 30 shared + 300) | ok, joint kind 'existing overlap', overlap >= 30 nt, product 530 bp |
| 6 | the shared 30 bases present twice in b | check 'Homology unique (a -> b)' = bad, design not ok |
| 7 | only 16 shared bases | joint still 'existing overlap' but check 'Homology (a -> b)' = warn |
| 8 | Gibson circular/linear, Sidewinder, and Golden Gate designs | for every joint with kind != 'none', the bases at joint.pos in the product equal joint.seq (read round the circle) |
| 9 | Golden Gate joints / Sidewinder joints | Golden Gate: a and b spans both (pos, pos+4); Sidewinder: a[0] == b[0] == pos + TOEHOLD |

**Review notes:**

- Long test covering four separate behaviours; a failure in the last block does not say which engine is wrong without reading the loop.
- Tm threshold 52 is hard-coded here and also in the primer tests.

## `test_golden_gate_products_and_overhang_checks`

**Result:** ✅ passed in 0.02 s · [line 165](../../tests/test_assembly.py#L165) · **Area:** Golden Gate assembly · **Checked against:** hand-computed · **Strength:** 🟢 strong

**Why it matters:** Without it, a Golden Gate design could produce a circle with the wrong overhangs, accept overhangs that cross-ligate or repeat, or fail to warn about leftover enzyme sites, so the cloning would not work at the bench.

Checks a vector-plus-inserts assembly end to end (primer design skips the vector, product is one circle, no BsaI site left), the overhang actually carried by each end, the errors for a mismatched end, a missing site, an internal site or a malformed destination vector, the per-junction quality checks, and the published validated overhang set.

| # | Input | Expected output |
|---|---|---|
| 1 | vector (two outward BsaI sites, end bases GGAG and TCGC) + inserts 700 bp and 220 bp, BsaI | before primers: joints all 'none', not ok; auto_primers returns 2 (vector is cut, not amplified); then ok, 4 rows, every junction 4 nt, 3 distinct overhangs, product a rotation of design_golden_gate's, all 'unique' checks ok, 'BsaI sites left in product' ok |
| 2 | insert reverse-complemented | design not ok until auto_primers is run again, then ok |
| 3 | vector with no / with two same-direction BsaI sites | ValueError matching 'needs exactly 2' and 'point the same way' |
| 4 | ring of three fragments with overhangs GGAG, TCGC, CAGT | ok; joint overhangs read back as TCGC, CAGT, GGAG; product 3*(4+250) bp; all 'validated' checks ok; no bad or warn check |
| 5 | one fragment's left overhang changed to TCGA | joint kind 'none', problem names both TCGC and TCGA, not ok |
| 6 | one fragment given no BsaI site | problem text contains 'no BsaI site', not ok |
| 7 | a BsaI site planted inside a fragment core | check 'BsaI sites in g0' = bad, not ok |
| 8 | overhang sets: [GGAG,GGAG,TCGC] / [GATC,...] / [GGTT,GGTA,TCGC] / [AAAA,...] / [CAGT,ACTG,TCGC] | ('unique','bad'); ('not palindromic','bad') (absent for GGTT); ('orthogonal','warn') and ('validated','info'); ('GC mix','warn'); ('orthogonal','bad') for a reverse-complement pair |
| 9 | overhangs.VALIDATED and POOR_EXAMPLE | 25 distinct 4 nt overhangs, no bad/warn from assess(); the poor set yields warns; pair_risk(GGAG,CTCC)='bad', pair_risk(GGAG,GGAA)='warn'; junction_checks(GGAG,[TCGC,CAGT]) all ok |
| 10 | pick_overhangs candidates | prefers the validated one (GGAG over ACCA) and never repeats an overhang already used |

**Review notes:**

- Also re-tests the overhangs module directly (OH.assess, pair_risk, junction_checks), which could live in its own test.

## `test_sidewinder_oligos_toeholds_and_barcodes`

**Result:** ✅ passed in 0.05 s · [line 262](../../tests/test_assembly.py#L262) · **Area:** Sidewinder assembly · **Checked against:** property/invariant · **Strength:** 🟢 strong

**Why it matters:** Without it, the ordered oligos could fail to tile the intended product, repeat a toehold so the assembly joins in the wrong place, or reuse barcodes, and the user would only discover it after paying for synthesis.

Checks that Sidewinder's toeholds are taken from the product itself (scarless), straddle the fragment boundaries, slide when they would be ambiguous, that the coding oligos tile the product exactly, that barcode arms are the right length, composition and mutually orthogonal, that the product is always linear and carries the fragments' features, and that the same input gives the same oligos.

| # | Input | Expected output |
|---|---|---|
| 1 | three fragments 100/120/90 bp (first seed of 60 that gives a clean design) | d.ok; product is the plain concatenation, 310 bp; two joints of kind 'toehold' (a,b) and (b,c) |
| 2 | each toehold | length TOEHOLD, equals product[pos:pos+TOEHOLD], centre within MAX_SHIFT of the boundary (100, 220), occurs exactly once in the product |
| 3 | the -coding oligo rows | a/b/c coding strands equal product[:e0], product[e0:e1], product[e1:]; a and b each end with their toehold |
| 4 | the -barcode oligo rows | each starts with the reverse complement of its coding strand; arms of BARCODE, BARCODE, 0 nt; the two arms differ; every oligo <= MAX_OLIGO nt |
| 5 | simulate(parts, d) | feature 'pro' at (0,50), 'cds1' at (109,160); no fragment-name features; sequence equals the product. With preview=True, markers a/b/c at (0,100),(100,220),(220,310) |
| 6 | fragment c made to contain the exact 10 bases the boundary toehold would use | the planted repeat occurs twice in the product; the chosen toehold differs from it and is unique; still within MAX_SHIFT of base 100; check text 'occurs once in the product'; the coding oligos still concatenate to the product |
| 7 | the two barcode arms, 8-mers of both strands | no shared 8-mer between the arms; GC count 8-10; no 4-base homopolymer; no 8-mer shared with any toehold |
| 8 | re-analysing the same fragments / two different seeds | identical oligo rows; equal product length but different barcode arms |
| 9 | analyse with circular=True | design still linear, two joints, 'Linear product' check present, simulated document not circular, same product as the linear run; Gibson in contrast closes with 2 joints and circular=True |
| 10 | one fragment only / three ~300 bp fragments | ValueError 'at least two'; and for long fragments not ok with check 'Oligos at most <MAX_OLIGO> nt' = bad |

**Review notes:**

- The opening loop over 60 seeds hides how often a random design is rejected; if the toehold picker regressed to accept almost nothing, the test would still pass on a later seed.

## `test_products_keep_every_orf_on_its_own_bases`

**Result:** ✅ passed in 0.24 s · [line 402](../../tests/test_assembly.py#L402) · **Area:** Feature geometry in assembly products · **Checked against:** independent implementation · **Strength:** 🟢 strong

**Why it matters:** Without it, an assembled product could annotate a CDS a few bases off, so the user reads a shifted frame and the wrong protein without any error being shown.

Builds fragments that are each a single ORF with the CDS annotated over all of it, assembles them with every method, and verifies in the product that each CDS still covers exactly its original bases, keeps codon_start=1 and translates as M...* with no internal stop (Biopython translation).

| # | Input | Expected output |
|---|---|---|
| 1 | Gibson, linear, overlaps forced to 20,21,...,25 bp | both joints 'existing overlap', product length = sum of fragments - 2*overlap, every ORF intact and translatable |
| 2 | the same, circular | three 'existing overlap' joints, product shorter by 3*overlap, product circular, every ORF intact (including over the origin) |
| 3 | Golden Gate, vector + three ORF fragments, seeds 0-5 | four joints of kind 'overhang', circular product, every ORF still annotated over all of its bases despite the overhang eating 0-4 bases |
| 4 | Sidewinder on three short ORFs (30/40/35 codons) | two joints; ORFs intact |
| 5 | middle ORF on the minus strand, with Gibson linear, Gibson circular and Golden Gate | no joint of kind 'none'; ORFs intact on both strands |
| 6 | fragments whose right end carries 30 unannotated bases of the next ORF (shared overlap) | two 'existing overlap' joints of 30 nt; all three ORFs annotated over exactly their designed bases |

**Review notes:**

- Slowest test in the file (many seeds x methods); the Golden Gate loop over 6 seeds adds little once one seed passes.

## `test_the_product_carries_the_features_and_the_primers_of_its_fragments`

**Result:** ✅ passed in 0.05 s · [line 460](../../tests/test_assembly.py#L460) · **Area:** Assembly product annotation · **Checked against:** hand-computed · **Strength:** 🟢 strong

**Why it matters:** Without it, the simulated product could lose or duplicate the fragments' features and primer_bind annotations, or credit bases to a primer tail that the fragments already shared, misleading the user about what the design actually adds.

Checks that features inside an overlap are annotated once and at the right place, that primers become primer_bind features covering exactly where they bind (origin-crossing included) and are recognised by the editor, that preview markers and tail/joint colouring are kept out of the plain product, and that junction_facts/junction_summary attribute each junction's bases to the right primer or to shared bases.

| # | Input | Expected output |
|---|---|---|
| 1 | 1000 bp plasmid cut into 0-500 and 470-1000, feature 476..495 inside the overlap | one 'existing overlap' joint of 30 nt, 1000 bp product; features exactly [('in overlap',[(475,495)]), ('left',[(50,120)]), ('right',[(800,900)])] |
| 2 | Gibson and Golden Gate circular designs with an annotated insert | design ok; feature 'pro' appears exactly once; product circular with length = product(); piece_doc length equals the fragment's |
| 3 | three ORF fragments, Gibson circular, 20 bp overlaps | 6 distinct primers; in both plain and preview products, primer_bind features match name+sequence exactly, each covers the full primer length, F on +1 and R on -1; one runs across the origin in two spans |
| 4 | Editor opened on the product | file_primer_set and primer_sites contain all 6 primers (found by their bases) |
| 5 | colour qualifiers of the preview product | TAIL_COLOR and JOINT_COLOR distinct from fragment colours; 3 'homology' features; tail features total 60 bp (3 x 20 bp made of two 10 nt tails) |
| 6 | simulate plain / preview with tails=False | plain product carries no fragment, tail or homology marker; tails=False keeps fragment names and 3 joint features but no tails |
| 7 | junction_facts / junction_summary for each of the 3 junctions | length 20, shared 0, two added tails of 10 nt, one on 'start' and one on 'end'; text contains '20 bp homology', 'adds 10 nt' twice, 'start of Fragment' and 'end of Fragment' |
| 8 | Golden Gate product / Sidewinder product | 2*(n fragments)-2 primer_bind features and 'overhang' in the junction summary; Sidewinder product has no primer_bind feature |
| 9 | two fragments overlapping 30 bp with no primers | facts (30, no added, shared 30); summary says the homology 'comes from the fragments' own bases' and 'no primer adds anything'; no primer_bind or tail features |
| 10 | overlap of 10 shared bases plus a 10 nt primer tail | facts (20, [('right-F',1,'start',10)], shared 10); summary ends with the sentence crediting 10 nt to right-F and 10 bp to the fragments; the tail feature sits at (480,490) |

**Review notes:**

- Very long; mixes feature transfer, primer annotation, preview colouring and junction wording in one function, so a failure needs the line number to interpret.
- One assertion checks a full English sentence verbatim, which will break on any rewording.

## `test_restriction_ligation_products`

**Result:** ✅ passed in 0.01 s · [line 583](../../tests/test_assembly.py#L583) · **Area:** Restriction / ligation cloning · **Checked against:** hand-computed · **Strength:** 🟢 strong

**Why it matters:** Without it, the predicted ligation product could be the wrong arc of the vector or have the insert at the wrong place, so the user would clone a construct different from the one shown.

Checks the exact base sequence of a directional EcoRI/BamHI ligation, of a compatible-but-different BamHI/BglII pair (hybrid sites destroyed), and of a blunt EcoRV/SmaI ligation; then that vector and insert features travel into the product correctly, and that a real cut-and-ligate of the demo plasmids reproduces the expected sequence.

| # | Input | Expected output |
|---|---|---|
| 1 | vector clean(300)+GAATTC+clean(20)+GGATCC+clean(400), insert cut EcoRI..BamHI, enzymes (EcoRI,BamHI) | ok; joints 'overhang' with sequences AATT and GATC; product = vec[327:] + vec[:301] + 'AATTC' + mid + 'G'; circular; overhang bases sit at the joint positions; both sites regenerated; simulate gives the same bases; junction cards titled '4 bp overhang' |
| 2 | vector with one BamHI site, insert flanked by two BglII sites | ok; both junctions GATC; 'ends fit' text names both enzymes and "5' GATC overhang"; 'sites' text says hybrid junction / both sites destroyed; product contains GGATCT and AGATCC but neither GGATCC nor AGATCT |
| 3 | EcoRV vector + SmaI-flanked insert | ok; joints 'blunt end'; junction card title 'blunt end'; product = blunt_vec[303:] + blunt_vec[:303] + 'GGG' + mid2 + 'CCC'; 'blunt end' appears in the junction HTML |
| 4 | vector features: CDS over the origin, misc_feature in the discarded stretch, CDS 40..100; insert CDS on the bottom strand | 'dropped' absent from the product; 'plain' and 'gene' cover the same bases as in their sources; 'gene' on the minus strand, 'plain' on the plus; the origin-crossing CDS still 51 bp |
| 5 | demo/pBAD30_RC.gb as vector, eilR cassette cut BamHI..XhoI from demo/pBAD30_eilR.gb | ok; product = seq[xho+1:] + seq[:bam+1] + seq[bam+1:xho+1], same length as the original plasmid; eilR feature present and starting with ATG; 'Vector self-ligation (pBAD30)' and 'Orientation' checks ok |

**Review notes:**

- Depends on two demo GenBank files; if those are edited the expected arcs change silently.

## `test_restriction_ligation_end_rules_and_warnings`

**Result:** ✅ passed in 0.01 s · [line 646](../../tests/test_assembly.py#L646) · **Area:** Restriction / ligation end compatibility · **Checked against:** hand-computed · **Strength:** 🟢 strong

**Why it matters:** Without it, the program could pair ends that cannot actually ligate (or refuse ones that can) and stay silent about an insert that can go in backwards or a vector that will re-close on itself.

Checks the end-compatibility rule base by base (same sticky bases, same overhang kind, blunt with blunt, 3' with 3'), the Orientation and Vector self-ligation warnings and where they are shown, the message for ends that do not match, the errors for a vector or insert not cut the right number of times, and the pieces a double digest produces.

| # | Input | Expected output |
|---|---|---|
| 1 | cut_end for EcoRI/BamHI/XhoI/EcoRV/SmaI/PstI/NsiI | EcoRI fits EcoRI; EcoRI not BamHI; sticky never fits blunt; EcoRV fits SmaI; PstI fits NsiI and is kind "3'"; same bases declared 5' do not pair with a 3' end |
| 2 | vector with one EcoRI site, insert flanked by two EcoRI sites | ok but Orientation = warn, text names 'ins' and '2 different products', and the check is in global_checks |
| 3 | the directional EcoRI/BamHI pair / a blunt EcoRV-EcoRV pair | Orientation ok for directional; warn for the blunt one |
| 4 | the same EcoRI/EcoRI vector | 'Vector self-ligation (vec)' = warn, text contains 'close on itself' and 'Dephosphorylate', shown globally |
| 5 | vector with BamHI and BglII sites (compatible ends) / directional vector | self-ligation warn; ok for the directional one |
| 6 | directional insert reverse-complemented | not ok, first joint 'none', problem says 'ends do not match' and names EcoRI and BamHI; junction card level 'bad', brief ('ends do not fit',''); flipping the insert makes it ok |
| 7 | blunt vector with a sticky EcoRI insert | design not ok |
| 8 | uncut vector / insert cut only once | ValueError 'needs 1 or 2' and ValueError 'cut once' |
| 9 | sequence with EcoRI/BamHI/EcoRI/BamHI sites, digest linear and circular | pieces with end enzymes (EcoRI,BamHI),(BamHI,EcoRI),(EcoRI,BamHI) plus one extra (BamHI,EcoRI) when circular; re-cutting each piece gives the same two ends, no extra site, same kept length |
| 10 | the same sequence digested with EcoRI only | one piece of 102 bp |

## `test_assembly_reports_bad_input_instead_of_failing`

**Result:** ✅ passed in 0.26 s · [line 707](../../tests/test_assembly.py#L707) · **Area:** Assembly input validation · **Checked against:** property/invariant · **Strength:** 🟡 adequate

**Why it matters:** Without it, a zero-length or very short fragment could crash the assembly panel with a traceback instead of showing the user a red check saying what is wrong.

Runs every method (Gibson, Golden Gate, Sidewinder), circular and linear, with a fragment of 0, 5 and 18 bp and checks that analysis and automatic primer design report bad or warning checks rather than raising; also that an empty stretch of a source gives an empty fragment, and that a repeated Golden Gate overhang blocks the design.

| # | Input | Expected output |
|---|---|---|
| 1 | fragments 300 / {0,5,18} / 250 bp, methods 0-2, circular and linear, analysed as given | never ok, at least one bad check, every bad check carries explanatory text |
| 2 | the same after auto_primers | no exception; Sidewinder or size 0 still not ok with a bad check; otherwise 'Fragment length' = warn naming fragment B |
| 3 | the 0 bp fragment | doc.length == 0 and seq == '' |
| 4 | Frag over ranges (5,5), (0,0) and (50,20) of a 100 bp source | seq and core_seq empty, tails (0,0) |
| 5 | Golden Gate ring with overhangs GGAG, GGAG, TCGC | a bad '...unique' check whose text contains 'same overhang'; not ok |
| 6 | the same ring with GGAG, TCGC, CAGT | design ok |

**Review notes:**

- Mostly asserts 'not ok and some bad check', so a wrong but non-crashing message would pass.
- The repeated-overhang case overlaps with the ('unique','bad') case in test_golden_gate_products_and_overhang_checks.
