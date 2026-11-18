# `tests/test_agent_cli.py`

[← all test files](README.md) · source: [`tests/test_agent_cli.py`](../../tests/test_agent_cli.py)

## `test_agent_reads_the_plasmid`

**Result:** ✅ passed in 0.02 s · [line 47](../../tests/test_agent_cli.py#L47) · **Area:** Command line for agents · **Checked against:** independent implementation · **Strength:** 🟡 adequate

**Why it matters:** An agent would plan 200 edits from wrong positions, strands or amino acids, and every one of them would be wrong.

inspect, feature and region give what the GenBank file says: length, topology, feature count, a reverse-strand CDS read from its own start codon, its protein, and 1-based inclusive regions.

| # | Input | Expected output |
|---|---|---|
| 1 | inspect demo/pBAD30_eilR.gb | length and topology of the Document, one entry per feature |
| 2 | feature araC (3655..4533, strand -1) | sequence = reverse complement of bases 3655..4533, starts ATG; protein starts M, ends *, 292 aa |
| 3 | region 1 12 | the first 12 bases of the Document |

**Review notes:**

- One demo plasmid only: a join or a feature across the origin is not exercised.

## `test_agent_mistakes_are_one_error_line`

**Result:** ✅ passed in 0.01 s · [line 64](../../tests/test_agent_cli.py#L64) · **Area:** Command line for agents · **Checked against:** hand-computed · **Strength:** 🟡 adequate

**Why it matters:** A mistake in the agent's input could end as a traceback or silent output that the agent cannot read or pass to the user.

Each user mistake gives exit code 1, nothing on stdout and exactly one `ERROR:` line on stderr that says what is wrong.

| # | Input | Expected output |
|---|---|---|
| 1 | missing file; unknown feature; region 0..12; region past the end; unknown recipe | exit 1, empty stdout, one ERROR line naming the problem |

## `test_every_recipe_example_still_runs`

**Result:** ✅ passed in 2.16 s · [line 72](../../tests/test_agent_cli.py#L72) · **Area:** Command line for agents · **Checked against:** smoke · **Strength:** 🟡 adequate

**Why it matters:** A recipe could keep telling agents to run a command that no longer exists or whose arguments changed.

Every command in the Example block of every recipe (and the index) runs and succeeds, in a temporary folder, with the one network call answered by canned replies, so the recipes cannot drift from the tool. The block is found even when a sentence stands between the heading and the code (kld, pcr and recombineering examples were not run before this was fixed).

| # | Input | Expected output |
|---|---|---|
| 1 | skills list, then each recipe's Example block, run in a temporary folder | the 14 function recipes and the index are listed; every example command exits 0 |

**Review notes:**

- Checks that examples run, not that their output is what the recipe says: each command has its own test for that.

## `test_hgvs_mutations_become_the_change_they_name`

**Result:** ✅ passed in 0.00 s · [line 111](../../tests/test_agent_cli.py#L111) · **Area:** Command line for agents · **Checked against:** hand-computed · **Strength:** 🟢 strong

**Why it matters:** A mutation written in HGVS could be turned into the wrong bases (off by one, wrong strand, wrong codon) and 200 wrong plasmids and primer pairs would be made from it.

Each HGVS form (substitution, deletion, insertion, deletion-insertion, duplication; c., o. and protein; one or three letters) becomes the exact plasmid change, with its HGVS text, predicted protein and effect, on a forward and on a reverse-strand copy of the same gene.

| # | Input | Expected output |
|---|---|---|
| 1 | gx (CDS 101..121 forward) c.4G>A, p.Ala2Thr, A2T, o.104G>A | (103, 104, 'A'), o.104G>A; the first also p.(Ala2Thr) missense |
| 2 | gx p.Glu3del; c.9_10insCCC; c.9_10insC; c.7_9dup | (106,109,''); (109,109,'CCC') p.(Glu3_Lys4insPro); frameshift; (109,109,'GAA') p.(Glu3dup) |
| 3 | gx p.Glu3_Lys4delinsGln; c.7_12delinsCAA | (106, 110, 'C') for both (only GAAA -> C differs) |
| 4 | gy (complement(401..421)) c.4G>A, p.Ala2Thr, p.Glu3del, c.9_10insCCC | (417,418,'T') o.418C>T; (417,418,'T'); (412,415,''); (412,412,'GGG') |
| 5 | gx:p.Ter7Gln ; gx:p.Met1Val | effect 'stop lost', protein p.(Ter7Glnext*?), no warning ; no warning for the changed start codon |

**Review notes:**

- Features made of several parts (joins) and changes beyond the ends of a CDS are not exercised.

## `test_hgvs_mistakes_say_what_to_fix_and_ambiguities_only_warn`

**Result:** ✅ passed in 0.00 s · [line 134](../../tests/test_agent_cli.py#L134) · **Area:** Command line for agents · **Checked against:** hand-computed · **Strength:** 🟡 adequate

**Why it matters:** A wrong reference base or amino acid in a mutation table could be applied silently, and an ambiguous codon choice could be made without telling the user.

Every mistake (wrong reference, wrong residue, nothing to change, outside the plasmid, past the CDS end, unknown feature, non-adjacent insertion, position 0, missing feature name, unreadable text, c.*N) raises one readable error; a unique nearest codon gives no warning, tied codons give a default and a warning naming the alternatives.

| # | Input | Expected output |
|---|---|---|
| 1 | gx:c.4C>A, p.Lys2Thr, p.Ala2Ala, o.600A>T, c.50G>A, nope:c.1A>T, c.9_11insA, c.0A>T, c.4G>A, 'lacZ del 5', c.*3A>T | an error each, with the words 'not C', 'is Ala (GCT), not Lys', 'already', 'outside the plasmid', 'past the end', 'No feature called', 'neighbouring', 'no position 0', 'name the feature', 'Cannot read', 'stop codon' |
| 2 | gx:p.Ala2Gly; p.Ala2Ser; p.Glu3Asp | G with no warning; T; C with a warning naming GAC and GAT |

**Review notes:**

- The warning for a change that can slide along a repeat is not tested.

## `test_a_change_through_the_origin_is_made_and_the_file_keeps_its_origin`

**Result:** ✅ passed in 0.10 s · [line 149](../../tests/test_agent_cli.py#L149) · **Area:** Command line for agents · **Checked against:** hand-computed · **Strength:** 🟢 strong

**Why it matters:** A mutation that spans the origin of the plasmid could be refused or give a shifted, wrong plasmid.

On a 521 bp circle, changes written through the origin (first position after the last) are made, and each output file equals the hand-edited string with its origin kept at base 1 unless the change replaced it.

| # | Input | Expected output |
|---|---|---|
| 1 | o.520_3del | hgvs o.520_3del; file = bases 4..519 (new origin after the change) |
| 2 | o.521_1insATG | file = plasmid + ATG (origin unchanged) |
| 3 | o.520_2delinsGG | file = bases 3..519 + GG |
| 4 | o.520_2dup | file = bases 1..2, then the duplicated 520..2, then 3.. (origin kept) |
| 5 | o.3_1del (nearly the whole plasmid) | error row, the others are made |

## `test_kld_batch_gives_the_files_and_primers_the_panel_gives`

**Result:** ✅ passed in 0.08 s · [line 164](../../tests/test_agent_cli.py#L164) · **Area:** Command line for agents · **Checked against:** independent implementation · **Strength:** 🟢 strong

**Why it matters:** The batch could give files or primers different from what the user gets by hand in the KLD panel, or overwrite earlier results, or accept a malformed table.

A 4-row table (substitution, reverse-strand substitution, deletion, wrong reference) writes the 3 good edited plasmids with both primers drawn on them, a primers table, a report and the settings; the plasmids are the string edits of the original; the primers are those of kld.design_selected; the failed row is reported and gives exit code 1; a second run into the same folder and a table with extra columns are refused.

| # | Input | Expected output |
|---|---|---|
| 1 | kld run t.gb m.csv --out out (rows sub, rev, del, bad) | exit 1; summary 3 ok 1 error; out holds del.gb primers.csv report.jsonl rev.gb settings.json sub.gb |
| 2 | each written .gb | bases = original with the change; primer_bind features KLD-F and KLD-R; primers equal kld.design_selected's |
| 3 | same run again; table with a note column | 'already exists'; 'exactly id,mutation' |

**Review notes:**

- Settings are the defaults of the test machine: a saved setting that changes the primers is not exercised.

## `test_pcr_batch_gives_the_pairs_the_panel_gives`

**Result:** ✅ passed in 1.61 s · [line 198](../../tests/test_agent_cli.py#L198) · **Area:** Command line for agents · **Checked against:** independent implementation · **Strength:** 🟢 strong

**Why it matters:** The batch could give primers different from what the PCR panel gives for the same region, or hide that a target has no pair, or overwrite earlier results.

pcr check and run give, for 150..400 and 200..450, the pairs of primer_design.design_region (same sequences and sizes), sitting on the two ends of the target; a 12 bp target and an unknown feature are error rows with the reason (exit 1); check writes nothing; run writes primers.csv, report.jsonl and settings.json and shows an overridden Tm as an override; a second run into the same folder is refused.

| # | Input | Expected output |
|---|---|---|
| 1 | pcr check t.gb with rows a 150..400, b 200..450, short 150..161, bad nope | exit 1; a and b ok, short 'No pair found', bad 'No feature called'; no output folder |
| 2 | pairs of a and b | equal to design_region(doc, 149, 400) and (199, 450); product 251 bp with fwd.start 150 and rev.end 400; fwd primer = the plasmid's own bases |
| 3 | pcr run --out out --tm 62, then again | tm override reported; three files; primers.csv lines = pairs found; the second run refused ('already exists') |

**Review notes:**

- The genome checks run against the demo genome of the repository, whatever the machine's own setting.

## `test_recombineering_batch_gives_the_designs_the_panel_gives`

**Result:** ✅ passed in 0.04 s · [line 231](../../tests/test_agent_cli.py#L231) · **Area:** Command line for agents · **Checked against:** hand-computed · **Strength:** 🟢 strong

**Why it matters:** The batch could place the cassette in the wrong orientation or place, give primers whose tails are not the genome's homologies, or overwrite earlier results.

For a 521 bp circular genome, a 120 bp cassette and rows replacing a CDS, a 11 bp region in reverse orientation and a CDS deleted with no cassette, the written genomes are the hand-built strings, the primers carry the 50 bp homologies as tails and equal recombineering.design's rows, the removed feature is warned about, a bad feature name and a bad orientation are error rows, nothing is overwritten and a homology longer than the genome is refused for the row.

| # | Input | Expected output |
|---|---|---|
| 1 | rows ko (gx, cas.gb), flip (150..160, reverse), del (gy, none), bad (nope), ori (orientation up) | exit 1; ok ok ok error error; files del.gb flip.gb ko.gb primers.csv report.jsonl settings.json |
| 2 | ko.gb and flip.gb, del.gb | seq[:100]+cas+seq[121:]; seq[:149]+rc(cas)+seq[160:]; seq[:400]+seq[421:] (no primers) |
| 3 | ko primers | forward = seq[50:100] + the first bases of the cassette; reverse = rc(seq[121:171]) + the first bases of rc(cassette); equal to design(); Rec-F and Rec-R drawn on ko.gb; warning 'removes the CDS gx' |
| 4 | second run into the same folder; --homology 600 with --primers-only | 'already exists'; 'no room' for the row and no .gb written |

**Review notes:**

- A genome with a linear topology and a region near its end is not exercised.

## `test_batches_keep_going_past_bad_rows_and_never_collide_on_file_names`

**Result:** ✅ passed in 0.09 s · [line 262](../../tests/test_agent_cli.py#L262) · **Area:** Command line for agents · **Checked against:** hand-computed · **Strength:** 🟢 strong

**Why it matters:** One bad row could end a 200-row run, two rows could overwrite each other's file, or a value the windows would refuse could pass without the user knowing.

A table (saved as cp1252, as Excel does) with ids differing only by case, unsafe characters, a Windows device name, a ragged row, an empty cell, a mutation that changes nothing and a row that fails in an unforeseen way: every other row is processed, every ok row has its own file (unique in any letter case), renames are told in the row's warnings, the errors are loud on stderr (with the traceback of the unforeseen one), a no-op mutation is accepted with a warning and no primers, and a Tm beyond the window's range is used with a responsibility warning. The same promises are then checked on pcr (a ragged row and an empty cell are that row's error, the others are processed) and on recombineering (ids that differ by case get two files, an empty cassette is that row's error).

| # | Input | Expected output |
|---|---|---|
| 1 | ids A, a, 'bad id!', CON, café; rows ragged (3 cells), empty (no mutation), noop (o.104G>G), boom (resolve raises RuntimeError) | files A.gb _CON.gb a_2.gb bad_id_.gb caf_.gb noop.gb; ragged '3 cells but the header has 2', empty 'mutation is empty', boom 'Unexpected RuntimeError' with a traceback; exit 1 |
| 2 | noop | effect 'no change', no primers, 'o.104=', warning 'no change'; its file is the plasmid unchanged |
| 3 | kld check --tm 80 | exit 0; settings max 72.0; top-level warning 'responsibility'; WARNING on stderr |
| 4 | pcr check with rows ok, ragged (3 cells), empty target; recombineering run with ids A, a, a ragged row, an empty cassette | pcr: ok, error '3 cells but the header has 2', error 'target is empty', exit 1; recombineering: files A.gb and a_2.gb, 'cassette is empty' |

**Review notes:**

- Codon and gibson are not in this table-level test: codon has its own bad-row case in its test, gibson stops on a bad row by design.

## `test_sites_command_gives_the_cuts_and_fragments_of_the_windows`

**Result:** ✅ passed in 0.02 s · [line 326](../../tests/test_agent_cli.py#L326) · **Area:** Command line for agents · **Checked against:** hand-computed · **Strength:** 🟢 strong

**Why it matters:** An agent could be told an enzyme cuts at the wrong base, or get fragment sizes that differ from the Digest tool's.

On a 400 bp circle with GAATTC planted at 201..206 and GGATCC at 301..306, `sites` reports the sites, the notation, where each strand is cut and the digest fragments; the default list is the window's (unique cutters of 6 or more bases).

| # | Input | Expected output |
|---|---|---|
| 1 | --enzymes EcoRI,BamHI,NotI --digest EcoRI,BamHI | EcoRI G^AATTC site 201..206, top cut after base 201 and bottom after 205; BamHI top cut after 301; NotI count 0; fragments 300 (302..201, BamHI to EcoRI) and 100 (202..301) |
| 2 | sites with no option | exactly the enzymes find_cuts says cut once and recognise 6 bases or more |
| 3 | --enzymes Nope | ERROR listing the enzymes openDNA knows |

**Review notes:**

- The enzyme list is openDNA's 52 common enzymes: an enzyme outside it is an error, not a search.

## `test_guides_command_gives_the_guides_and_marks_of_the_spacer_finder`

**Result:** ✅ passed in 0.02 s · [line 342](../../tests/test_agent_cli.py#L342) · **Area:** Command line for agents · **Checked against:** independent implementation · **Strength:** 🟢 strong

**Why it matters:** An agent could report a guide at the wrong place, or a different list or mark from the Spacer Finder window's.

A planted SpCas9 site (20 nt guide + AGG at 101..123) is found with its PAM, strand and cut position; the whole list and the marks equal cas9.candidates / scan / verdict (what the window calls); a region keeps to its bases; a mismatch limit beyond the window's range is used with a responsibility warning.

| # | Input | Expected output |
|---|---|---|
| 1 | --no-genome --top 400 | guide GCGCATGCAAGTCCGATGAC at 101..123 strand +, PAM AGG, cut before base 118; guides and marks equal the module's |
| 2 | --region 90..130 --max-mm 9 | every guide inside 90..130; warning 'responsibility' |
| 3 | --nuclease Cas99 | ERROR 'No nuclease called' |
| 4 | --top -5; --top 2 | ERROR '--top must be at least 1'; exactly 2 guides |

**Review notes:**

- The background-genome counts are not compared (the test runs --no-genome); they come from the same cas9.scan.

## `test_repeats_command_gives_the_stretches_the_compare_view_gives`

**Result:** ✅ passed in 0.01 s · [line 364](../../tests/test_agent_cli.py#L364) · **Area:** Command line for agents · **Checked against:** hand-computed · **Strength:** 🟢 strong

**Why it matters:** An agent could miss a repeat that can recombine in the cell, or place it wrongly.

A 30 bp stretch planted twice (21..50 and 221..250) is reported as one direct repeat of length 30, equal to repeats.pair_repeats; with a second file the shared stretch is found at the planted places in both.

| # | Input | Expected output |
|---|---|---|
| 1 | repeats p.gb --min-len 20 | direct, start_in_a 21, start_in_b 221, length 30, 0 differences; same as pair_repeats |
| 2 | repeats p.gb q.gb --min-len 25 | direct stretches at (21, 101) and (221, 101) |

**Review notes:**

- Inverted repeats are not planted here; they are covered by the repeats module's own tests.

## `test_search_command_finds_dna_on_both_strands_and_names`

**Result:** ✅ passed in 0.02 s · [line 375](../../tests/test_agent_cli.py#L375) · **Area:** Command line for agents · **Checked against:** hand-computed · **Strength:** 🟢 strong

**Why it matters:** An agent could say a sequence is absent when it is on the other strand, or give a position that is off by one.

A 19 nt site planted at 51..69: its reverse complement is found as a minus-strand hit at 51..69, not with --single-strand, a one-base change is found with --fuzzy 1, a name finds the feature, and the hits equal search.search_sequence.

| # | Input | Expected output |
|---|---|---|
| 1 | search p.gb <reverse complement of the site> | one dna hit 51..69 strand -1 |
| 2 | same with --single-strand | no hit |
| 3 | site with its last base changed, --fuzzy 1 | hit at 51 on the + strand |
| 4 | search demo araC | first hit is the feature araC, same as search_sequence |

## `test_every_command_has_a_recipe_example_and_a_test_that_compares_it_with_the_window`

**Result:** ✅ passed in 0.01 s · [line 389](../../tests/test_agent_cli.py#L389) · **Area:** Command line for agents · **Checked against:** property/invariant · **Strength:** 🟡 adequate

**Why it matters:** A command could be added with no recipe for agents, an example that never runs, or no test that holds it to the window's function.

Every command of the parser except the four that read a file has a runnable example in some recipe and a test named test_<command>_... in this file.

| # | Input | Expected output |
|---|---|---|
| 1 | every subcommand of oDNA-cli | an example line starting with the command in a recipe, and a test function test_<command>_... |

**Review notes:**

- Checks that a test exists by name, not that it compares the command with the window: that is for the review of each test.

## `test_align_command_finds_pieces_and_structural_changes_like_local_align`

**Result:** ✅ passed in 0.02 s · [line 399](../../tests/test_agent_cli.py#L399) · **Area:** Command line for agents · **Checked against:** hand-computed · **Strength:** 🟢 strong

**Why it matters:** An agent could report a read as matching where it does not, or miss that a query has a deletion or is inverted.

A 60 bp piece of a 400 bp plasmid is found at 101..160 with identity 1.0; its reverse complement as a minus-strand piece; a query with 30 bp left out is two pieces, the second 'deletion from the query'; the pieces equal local_align.find_hits (what the window calls); a protein is refused.

| # | Input | Expected output |
|---|---|---|
| 1 | seq[100:160] | one piece, query 1..60, plasmid 101..160, strand +, identity 1.0 |
| 2 | reverse complement of it | 101..160, strand -1 |
| 3 | seq[100:130] + seq[160:190] | two pieces from 101..130 to ..190, 28 to 30 bases apart, 'deletion from the query' |
| 4 | MKLVPQ | ERROR |

**Review notes:**

- The second piece may start a base or two early (bases next to a junction can match by chance): the test allows 28 to 30.

## `test_msa_command_aligns_the_files_and_never_overwrites`

**Result:** ✅ passed in 0.04 s · [line 415](../../tests/test_agent_cli.py#L415) · **Area:** Command line for agents · **Checked against:** hand-computed · **Strength:** 🟢 strong

**Why it matters:** An agent could report a wrong identity between two variants, or overwrite an alignment file.

A 200 bp sequence and a copy with 10 bases deleted are aligned: removing the gaps gives back each sequence, the copy has 10 gaps, identity to the first is 100 and 95 percent, the FASTA written parses back to the same rows, an existing output file is refused, one file is refused.

| # | Input | Expected output |
|---|---|---|
| 1 | a.gb b.gb --out aln.fa | width 200; rows without gaps equal a and b; identities 100.0 and 95.0; aln.fa parses to the same rows |
| 2 | same --out again | ERROR 'already exists' |
| 3 | a.gb alone | ERROR 'at least two' |

**Review notes:**

- Only one deletion is tested: large rearrangements are beyond the built-in aligner (the recipe says so).

## `test_codon_batch_recodes_in_place_keeping_the_protein`

**Result:** ✅ passed in 0.04 s · [line 431](../../tests/test_agent_cli.py#L431) · **Area:** Command line for agents · **Checked against:** hand-computed · **Strength:** 🟢 strong

**Why it matters:** An agent could write a recoded plasmid whose protein changed, whose reverse-strand gene was written the wrong way round, or whose bases differ from what the Recode window gives.

With the Highest frequency method each codon of the planted CDS becomes the most used synonym of the table (hand-read from the CSV: ATG GCG GAA AAA TTT GGC TAA); the plasmid has only that stretch changed, for a reverse-strand CDS the reverse complement is written; with a seed the frequency-matched result equals codon_tools.recode and keeps the protein; a bad name is that row's error; GC limits that contradict are refused.

| # | Input | Expected output |
|---|---|---|
| 1 | gx and gy with --method 'Highest frequency' | both rows ATGGCGGAAAAATTTGGCTAA, 2 codons changed; fwd.gb = seq[:100] + best + seq[121:]; rev.gb = seq[:400] + rc(best) + seq[421:] |
| 2 | a row naming a missing feature | that row 'No feature called', exit 1 |
| 3 | --gc-low 70 | ERROR '--gc-low is above' |
| 4 | codon check ... --no-gc | settings show gc_limits 'none' and no gc_low / gc_high / gc_window: they say what was applied |

**Review notes:**

- The usage table is an approximate E. coli table in the examples: the test checks the logic against the table given, not that the table is accurate.

## `test_export_command_writes_the_formats_and_never_overwrites`

**Result:** ✅ passed in 0.04 s · [line 460](../../tests/test_agent_cli.py#L460) · **Area:** Command line for agents · **Checked against:** independent implementation · **Strength:** 🟡 adequate

**Why it matters:** An agent could hand the user a file in the wrong format, or overwrite one.

Each format is written and equals the text the Export menu's functions give (fasta, gff3, bed, tbl, proteins); the proteins file holds MAEKFG for gx; the GenBank copy keeps the sequence; an existing target is refused.

| # | Input | Expected output |
|---|---|---|
| 1 | export t.gb <format> out | file text equals exports.<format> of the same document |
| 2 | proteins | >gx ... MAEKFG |
| 3 | genbank | same sequence as the source |
| 4 | an existing file | ERROR 'already exists' |

**Review notes:**

- The expected text comes from the same exports functions: it proves the command adds nothing and loses nothing, not that the formats are right (see the export tests).

## `test_order_command_writes_the_supplier_file_or_nothing`

**Result:** ✅ passed in 0.03 s · [line 476](../../tests/test_agent_cli.py#L476) · **Area:** Command line for agents · **Checked against:** hand-computed · **Strength:** 🟢 strong

**Why it matters:** An agent could order a primer with a typo, or write half an order.

Two primers become an IDT plate file (wells A1, B1, scale and purification as asked, sequences upper-cased) or a plain name,sequence file; a table with a primer that is too short or has a bad base writes nothing and names both rows; an existing file is refused.

| # | Input | Expected output |
|---|---|---|
| 1 | --scale 100nm --purification PAGE | Well Position,Name,Sequence,Scale,Purification then A1 F1 ... B1 R1 ... |
| 2 | --format plain | F1,ATGC... / R1,TTAG... |
| 3 | a table with 'short' (4 bases) and 'bad' (an X) | ERROR naming both, no file |

## `test_variants_command_finds_what_the_simulated_reads_were_given`

**Result:** ✅ passed in 0.04 s · [line 493](../../tests/test_agent_cli.py#L493) · **Area:** Command line for agents · **Checked against:** independent implementation · **Strength:** 🟢 strong

**Why it matters:** An agent could tell a user that a mutation is in their plasmid when it is not, or miss one that is.

On reads simulated from the demo plasmid with known changes (truth.tsv): the 100 percent substitutions at 4000 (C>A) and 5300 are consensus, the 30 percent one at 4200 is a variant (a mixed colony), the 1 bp insertion at 4700 and the 3 bp deletion at 5400 are found, the 100 bp deletion shows as a gap; the list equals alignments.scan_variants; --apply writes the consensus changes into a copy (length 5936 + 1 - 3 - 1).

| # | Input | Expected output |
|---|---|---|
| 1 | variants demo reads.bam | snp 4000 C>A consensus; snp 4200 A>T variant, fraction 0.2..0.8; ins 4700 A; del 5400 GGA; gap at 2000 |
| 2 | --apply fixed.gb | length 5933, base 4000 is A (C in the source), changes = the consensus ones that can be applied |
| 3 | --min-depth 0 | warning 'responsibility' |
| 4 | a GenBank file as the reads | ERROR 'Cannot read' |

**Review notes:**

- The truth is the simulation's own list: the extra changes the reads show (5564, 5680) are sequencing errors the simulation put in and are not asserted.

## `test_tool_settings_override_restores_the_outer_value_and_is_never_saved`

**Result:** ✅ passed in 0.00 s · [line 518](../../tests/test_agent_cli.py#L518) · **Area:** Command line for agents · **Checked against:** hand-computed · **Strength:** 🟡 adequate

**Why it matters:** A value given on the command line could stay in force after its command, or a nested override could lose the one around it.

Two nested overrides of one setting: the inner value is read inside, the outer one is back after it, and once both are left the saved value is read again.

| # | Input | Expected output |
|---|---|---|
| 1 | override 41, inside it override 42 | value 42 inside, 41 after the inner block |
| 2 | both blocks left | no override left; the saved value is read (not 41) |

**Review notes:**

- Single-threaded only: the overrides are process-wide (the command line runs one command at a time).

## `test_checkup_command_gives_the_verdicts_of_the_checkup_window`

**Result:** ✅ passed in 0.03 s · [line 529](../../tests/test_agent_cli.py#L529) · **Area:** Command line for agents · **Checked against:** independent implementation · **Strength:** 🟡 adequate

**Why it matters:** An agent could say a gene will be expressed when nothing drives it, or the reverse.

On a plasmid with two CDS and nothing else both are 'Not expressed' (no bacterial promoter upstream); with an annotated promoter, RBS and terminator around the forward one it is 'Expressed' with those elements named and the other stays 'Not expressed'; the verdicts equal checkup.run on the same document; a threshold beyond the window's range warns.

| # | Input | Expected output |
|---|---|---|
| 1 | checkup t.gb --no-predict --no-library | gx and gy Not expressed, reason 'no bacterial promoter upstream on this strand' |
| 2 | with promoter 40..70, RBS 88..99, terminator 130..160 annotated | gx Expressed, promoter Pgx, terminator Tgx; elements all 'annotated'; gy Not expressed |
| 3 | --min-promoter 99 | warning 'responsibility' |

**Review notes:**

- The verdict logic belongs to the checkup module and its own tests; this test only shows the command reports it and leaves the user's saved settings alone.

## `test_annotate_command_finds_library_parts_and_adds_only_new_ones`

**Result:** ✅ passed in 0.04 s · [line 552](../../tests/test_agent_cli.py#L552) · **Area:** Command line for agents · **Checked against:** independent implementation · **Strength:** 🟡 adequate

**Why it matters:** An agent could annotate a plasmid twice, or with parts that are not in it.

The longest bundled terminator planted into a random 600 bp sequence is found at the planted place with identity 1.0 and not yet annotated; the hits equal annotation.detect; --add writes a copy with the feature, and a second pass over that copy adds nothing because the part is now annotated.

| # | Input | Expected output |
|---|---|---|
| 1 | annotate p.gb --kinds terminator | a hit starting at base 101, strand +, identity 1.0, already_annotated false |
| 2 | --add annotated.gb | the copy has a terminator feature at the same place |
| 3 | annotate annotated.gb --add again.gb | 0 features added |
| 4 | --min-dna 10 | warning 'responsibility' |

**Review notes:**

- The part used is whichever bundled terminator is longest: if the library changes, the test follows it.

## `test_gibson_command_adds_overlaps_and_builds_the_product_the_panel_builds`

**Result:** ✅ passed in 0.07 s · [line 570](../../tests/test_agent_cli.py#L570) · **Area:** Command line for agents · **Checked against:** independent implementation · **Strength:** 🟡 adequate

**Why it matters:** An agent could give primers whose overlaps do not join the fragments, or a product that is not what the Assembly panel makes.

Two unrelated 300 bp fragments (the second cut to 11..290) become a 580 bp circle: the product, turned to start at the first fragment, is fragment 1 + the cut fragment 2; the primer of fragment 2 carries the end of fragment 1 as its 5' tail and the reverse primer of fragment 1 the start of fragment 2, the two tails adding to 20 to 40 bases; primers, check levels and product equal auto_primers / analyse / simulate; a bad row, a duplicate id and an existing folder are refused, an overlap beyond the window's range warns.

| # | Input | Expected output |
|---|---|---|
| 1 | gibson run f.csv | 580 bp circular product containing fa + fb[10:290]; fb-F tail is the end of fa; fa-R tail is the reverse complement of the start of fb; 2 junctions |
| 2 | region 0..5 | ERROR 'positions must satisfy' |
| 3 | two rows with the same id | ERROR 'must be different' |
| 4 | --min-overlap 5 | warning 'responsibility' |

**Review notes:**

- Only Gibson is covered: Golden Gate, Sidewinder and restriction/ligation have no command yet.
- The product's origin is the window's (a few bases before fragment 1): the test compares it as a circle.

## `test_ncbi_command_searches_and_fetches_through_the_import_dialogs_functions`

**Result:** ✅ passed in 0.03 s · [line 606](../../tests/test_agent_cli.py#L606) · **Area:** Command line for agents · **Checked against:** hand-computed · **Strength:** 🟡 adequate

**Why it matters:** An agent could send more than the user's words to NCBI, fetch a record other than the one named, or overwrite a file.

With the one network call replaced by canned replies: a search lists the record with its accession, length and topology, only the search words are sent (no key in the URL), `get` of an accession without its version writes the GenBank text of the current version, the results equal ncbi.search; an existing file, a missing --out, an accession that is not in the results and an empty search are each one ERROR line.

| # | Input | Expected output |
|---|---|---|
| 1 | ncbi search pUC19 --limit 3 | L09137.2, 60 bases, circular; URLs only to eutils.ncbi.nlm.nih.gov with retmax=3 |
| 2 | ncbi get L09137 --out puc.gb | a 60 bp GenBank file, record L09137.2 |
| 3 | get again to the same file / without --out / XYZ / search '  ' | ERROR already exists / needs --out / is not the accession / Type an accession |
| 4 | search --limit 0; get onto an existing --out | ERROR '--limit must be at least 1' / 'already exists', and no request is made to the (fake) NCBI |

**Review notes:**

- The real service is never called by the tests: an NCBI change in its answers is not caught here.
