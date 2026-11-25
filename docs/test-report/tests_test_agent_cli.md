# `tests/test_agent_cli.py`

[← all test files](README.md) · source: [`tests/test_agent_cli.py`](../../tests/test_agent_cli.py)

## `test_agent_reads_the_plasmid`

**Result:** ✅ passed in 0.02 s · [line 49](../../tests/test_agent_cli.py#L49) · **Area:** Command line for agents · **Checked against:** independent implementation · **Strength:** 🟡 adequate

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

**Result:** ✅ passed in 0.02 s · [line 66](../../tests/test_agent_cli.py#L66) · **Area:** Command line for agents · **Checked against:** hand-computed · **Strength:** 🟡 adequate

**Why it matters:** A mistake in the agent's input could end as a traceback or silent output that the agent cannot read or pass to the user.

Each user mistake gives exit code 1, nothing on stdout and exactly one `ERROR:` line on stderr that says what is wrong.

| # | Input | Expected output |
|---|---|---|
| 1 | missing file; unknown feature; region 0..12; region past the end; unknown recipe | exit 1, empty stdout, one ERROR line naming the problem |

## `test_every_recipe_example_still_runs`

**Result:** ✅ passed in 3.28 s · [line 74](../../tests/test_agent_cli.py#L74) · **Area:** Command line for agents · **Checked against:** smoke · **Strength:** 🟡 adequate

**Why it matters:** A recipe could keep telling agents to run a command that no longer exists or whose arguments changed.

Every command in the Example block of every recipe (and the index) runs and succeeds, in a temporary folder, with the one network call answered by canned replies, so the recipes cannot drift from the tool. The block is found even when a sentence stands between the heading and the code (kld, pcr and recombineering examples were not run before this was fixed).

| # | Input | Expected output |
|---|---|---|
| 1 | skills list, then each recipe's Example block, run in a temporary folder | the 14 function recipes and the index are listed; every example command exits 0 |

**Review notes:**

- Checks that examples run, not that their output is what the recipe says: each command has its own test for that.

## `test_hgvs_mutations_become_the_change_they_name`

**Result:** ✅ passed in 0.00 s · [line 115](../../tests/test_agent_cli.py#L115) · **Area:** Command line for agents · **Checked against:** hand-computed · **Strength:** 🟢 strong

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

**Result:** ✅ passed in 0.00 s · [line 138](../../tests/test_agent_cli.py#L138) · **Area:** Command line for agents · **Checked against:** hand-computed · **Strength:** 🟡 adequate

**Why it matters:** A wrong reference base or amino acid in a mutation table could be applied silently, and an ambiguous codon choice could be made without telling the user.

Every mistake (wrong reference, wrong residue, nothing to change, outside the plasmid, past the CDS end, unknown feature, non-adjacent insertion, position 0, missing feature name, unreadable text, c.*N) raises one readable error; a unique nearest codon gives no warning, tied codons give a default and a warning naming the alternatives.

| # | Input | Expected output |
|---|---|---|
| 1 | gx:c.4C>A, p.Lys2Thr, p.Ala2Ala, o.600A>T, c.50G>A, nope:c.1A>T, c.9_11insA, c.0A>T, c.4G>A, 'lacZ del 5', c.*3A>T | an error each, with the words 'not C', 'is Ala (GCT), not Lys', 'already', 'outside the plasmid', 'past the end', 'No feature called', 'neighbouring', 'no position 0', 'name the feature', 'Cannot read', 'stop codon' |
| 2 | gx:p.Ala2Gly; p.Ala2Ser; p.Glu3Asp | G with no warning; T; C with a warning naming GAC and GAT |

**Review notes:**

- The warning for a change that can slide along a repeat is not tested.

## `test_a_change_through_the_origin_is_made_and_the_file_keeps_its_origin`

**Result:** ✅ passed in 0.03 s · [line 153](../../tests/test_agent_cli.py#L153) · **Area:** Command line for agents · **Checked against:** hand-computed · **Strength:** 🟢 strong

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

**Result:** ✅ passed in 0.02 s · [line 168](../../tests/test_agent_cli.py#L168) · **Area:** Command line for agents · **Checked against:** independent implementation · **Strength:** 🟢 strong

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

**Result:** ✅ passed in 0.07 s · [line 202](../../tests/test_agent_cli.py#L202) · **Area:** Command line for agents · **Checked against:** independent implementation · **Strength:** 🟢 strong

**Why it matters:** The batch could give primers different from what the PCR panel gives for the same region, or hide that a target has no pair, or overwrite earlier results.

pcr check and run give, for 150..400 and 200..450, the pairs of primer_design.design_region (same sequences and sizes), sitting on the two ends of the target; a 12 bp target and an unknown feature are error rows with the reason (exit 1); check writes nothing; run writes primers.csv, report.jsonl and settings.json and shows an overridden Tm as an override; a second run into the same folder is refused. When every pair is marked bad the row warns and offers a wider window, and the window stays 0. When more pairs were found than the best 5 listed, the row warns of it. A window wider than the circle allows is narrowed and said.

| # | Input | Expected output |
|---|---|---|
| 1 | pcr check t.gb with rows a 150..400, b 200..450, short 150..161, bad nope | exit 1; a and b ok, short 'No pair found', bad 'No feature called'; no output folder |
| 2 | pairs of a and b | equal to design_region(doc, 149, 400) and (199, 450); product 251 bp with fwd.start 150 and rev.end 400; fwd primer = the plasmid's own bases |
| 3 | pcr run --out out --tm 62, then again | tm override reported; three files; primers.csv lines = pairs found; the second run refused ('already exists') |
| 4 | the pairs printed by default, then pcr check --full, and report.jsonl of the run | default: no 'checks', 'checks_passed' + 'issues' = 8 and issues are warn/bad only; --full and report.jsonl: all 8 checks |
| 5 | pcr check 150..400 with the pair marks forced to bad | the row's warning names '--window 0' and says to offer a wider one; settings window value 0 |
| 6 | the same check with the pair marks forced to bad | a warning 'Only the best 5 pairs are listed' as well as the wider-window one |
| 7 | pcr check 150..400 --window 150 on the 521 bp circle | warning 'narrowed to 135 bases' |

**Review notes:**

- The genome checks run against the demo genome of the repository, whatever the machine's own setting.

## `test_recombineering_batch_gives_the_designs_the_panel_gives`

**Result:** ✅ passed in 0.04 s · [line 252](../../tests/test_agent_cli.py#L252) · **Area:** Command line for agents · **Checked against:** hand-computed · **Strength:** 🟢 strong

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

**Result:** ✅ passed in 0.05 s · [line 283](../../tests/test_agent_cli.py#L283) · **Area:** Command line for agents · **Checked against:** hand-computed · **Strength:** 🟢 strong

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

**Result:** ✅ passed in 0.02 s · [line 347](../../tests/test_agent_cli.py#L347) · **Area:** Command line for agents · **Checked against:** hand-computed · **Strength:** 🟢 strong

**Why it matters:** An agent could be told an enzyme cuts at the wrong base, or get fragment sizes that differ from the Digest tool's.

On a 400 bp circle with GAATTC planted at 201..206 and GGATCC at 301..306, `sites` reports the sites, the notation, where each strand is cut and the digest fragments; the default list is the window's (unique cutters of 6 or more bases).

| # | Input | Expected output |
|---|---|---|
| 1 | --enzymes EcoRI,BamHI,NotI --digest EcoRI,BamHI | EcoRI G^AATTC site 201..206, top cut after base 201 and bottom after 205; BamHI top cut after 301; NotI count 0; fragments 300 (302..201, BamHI to EcoRI) and 100 (202..301) |
| 2 | sites with no option | exactly the enzymes find_cuts says cut once and recognise 6 bases or more |
| 3 | --enzymes Nope | ERROR listing the enzymes openDNA knows |
| 4 | sites --digest EcoRI,BamHI (no --enzymes) | the list stays the default single cutters and the digest gives fragments 100 and 300 |

**Review notes:**

- The enzyme list is openDNA's 52 common enzymes: an enzyme outside it is an error, not a search.

## `test_guides_command_gives_the_guides_and_marks_of_the_spacer_finder`

**Result:** ✅ passed in 0.03 s · [line 365](../../tests/test_agent_cli.py#L365) · **Area:** Command line for agents · **Checked against:** independent implementation · **Strength:** 🟢 strong

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

**Result:** ✅ passed in 0.02 s · [line 387](../../tests/test_agent_cli.py#L387) · **Area:** Command line for agents · **Checked against:** hand-computed · **Strength:** 🟢 strong

**Why it matters:** An agent could miss a repeat that can recombine in the cell, or place it wrongly.

A 30 bp stretch planted twice (21..50 and 221..250) is reported as one direct repeat of length 30, equal to repeats.pair_repeats; with a second file the shared stretch is found at the planted places in both.

| # | Input | Expected output |
|---|---|---|
| 1 | repeats p.gb --min-len 20 | direct, start_in_a 21, start_in_b 221, length 30, 0 differences; same as pair_repeats |
| 2 | repeats p.gb q.gb --min-len 25 | direct stretches at (21, 101) and (221, 101) |
| 3 | repeats a.gb twin.gb (another file, the same bases) | identical_sequences true; the whole length is the first repeat (direct, 1..1); one file with itself still says false |

**Review notes:**

- Inverted repeats are not planted here; they are covered by the repeats module's own tests.

## `test_synthesis_command_gives_what_the_synthesis_module_gives`

**Result:** ✅ passed in 3.37 s · [line 405](../../tests/test_agent_cli.py#L405) · **Area:** Command line for agents · **Checked against:** round-trip/self-consistency · **Strength:** 🟡 adequate

**Why it matters:** The command line and the Synthesis tab could drift apart if the command had its own copy of the logic.

The synthesis command returns the same probability, features and verdict as synthesis_salis.evaluate for the region, warns that a circular file is read straight, and refuses a backwards region. Its places equal evaluate's, 1-based on the file.

| # | Input | Expected output |
|---|---|---|
| 1 | random 260 bp circle, --region 1..200 | equal to evaluate(seq[:200]); region {1, 200}; places == evaluate's places + 1 (start) |
| 2 | whole circular file | warning 'straight piece' |
| 3 | --region 200..100 | ERROR line |

**Review notes:**

- Compares with the module, not with an external value: the module is checked in test_synthesis_salis.py.

## `test_search_command_finds_dna_on_both_strands_and_names`

**Result:** ✅ passed in 0.02 s · [line 420](../../tests/test_agent_cli.py#L420) · **Area:** Command line for agents · **Checked against:** hand-computed · **Strength:** 🟢 strong

**Why it matters:** An agent could say a sequence is absent when it is on the other strand, or give a position that is off by one.

A 19 nt site planted at 51..69: its reverse complement is found as a minus-strand hit at 51..69, not with --single-strand, a one-base change is found with --fuzzy 1, a name finds the feature, and the hits equal search.search_sequence.

| # | Input | Expected output |
|---|---|---|
| 1 | search p.gb <reverse complement of the site> | one dna hit 51..69 strand -1 |
| 2 | same with --single-strand | no hit |
| 3 | site with its last base changed, --fuzzy 1 | hit at 51 on the + strand |
| 4 | search demo araC | first hit is the feature araC, same as search_sequence |

## `test_every_command_has_a_recipe_example_and_a_test_that_compares_it_with_the_window`

**Result:** ✅ passed in 0.01 s · [line 434](../../tests/test_agent_cli.py#L434) · **Area:** Command line for agents · **Checked against:** property/invariant · **Strength:** 🟡 adequate

**Why it matters:** A command could be added with no recipe for agents, an example that never runs, or no test that holds it to the window's function.

Every command of the parser except the four that read a file has a runnable example in some recipe and a test named test_<command>_... in this file.

| # | Input | Expected output |
|---|---|---|
| 1 | every subcommand of oDNA-cli | an example line starting with the command in a recipe, and a test function test_<command>_... |

**Review notes:**

- Checks that a test exists by name, not that it compares the command with the window: that is for the review of each test.

## `test_align_command_finds_pieces_and_structural_changes_like_local_align`

**Result:** ✅ passed in 0.05 s · [line 445](../../tests/test_agent_cli.py#L445) · **Area:** Command line for agents · **Checked against:** hand-computed · **Strength:** 🟢 strong

**Why it matters:** An agent could report a read as matching where it does not, or miss that a query has a deletion or is inverted.

A 60 bp piece of a 400 bp plasmid is found at 101..160 with identity 1.0; its reverse complement as a minus-strand piece; a query with 30 bp left out is two pieces, the second 'deletion from the query'; the pieces equal local_align.find_hits (what the window calls); a protein is refused. Also: a mismatch just past the origin is named at its own place; limits reached are warnings (query longer than the circle, time limit) with --time-limit and --threads shown in settings.

| # | Input | Expected output |
|---|---|---|
| 1 | seq[100:160] | one piece, query 1..60, plasmid 101..160, strand +, identity 1.0; aligned_bp 60 |
| 2 | reverse complement of it | 101..160, strand -1 |
| 3 | seq[100:130] + seq[160:190] | two pieces from 101..130 to ..190, 28 to 30 bases apart, 'deletion from the query' |
| 4 | MKLVPQ | ERROR |
| 5 | one query base changed, forward and reverse-complemented | mismatch_list names position 131, the plasmid base and the query base on the top strand for both strands; 5 missing bases give a deletion 131..135 or two pieces |
| 6 | seq[-30:] + seq[0] + changed base + seq[2:30] (across the origin) | mismatch_list == [(position 2, changed base)] (it used to crash with IndexError) |
| 7 | seq + seq[:100] (longer than the 400 bp circle) | a warning 'longer than the circular template' |
| 8 | seq[100:160] --time-limit 0.000000001 --threads 2 | no piece; warnings 'The time limit...' and '... --time-limit 0 ...'; settings.align time_limit 1e-09, cores 2 |
| 9 | the gap query: notes of find_hits and the command's warnings | both empty |

**Review notes:**

- The second piece may start a base or two early (bases next to a junction can match by chance): the test allows 28 to 30.

## `test_msa_command_aligns_the_files_and_never_overwrites`

**Result:** ✅ passed in 0.09 s · [line 480](../../tests/test_agent_cli.py#L480) · **Area:** Command line for agents · **Checked against:** hand-computed · **Strength:** 🟢 strong

**Why it matters:** An agent could report a wrong identity between two variants, or overwrite an alignment file.

A 200 bp sequence and a copy with 10 bases deleted are aligned: removing the gaps gives back each sequence, the copy has 10 gaps, identity to the first is 100 and 95 percent, the FASTA written parses back to the same rows, an existing output file is refused, one file is refused. A multi-FASTA of the same two sequences is aligned whole, and FILE#NAME picks a record.

| # | Input | Expected output |
|---|---|---|
| 1 | a.gb b.gb --out aln.fa | width 200; rows without gaps equal a and b; identities 100.0 and 95.0; aln.fa parses to the same rows |
| 2 | same --out again | ERROR 'already exists' |
| 3 | a.gb alone | ERROR 'at least two' |
| 4 | the two files | each sequence names its file |
| 5 | both.fasta (two records) alone; then both.fasta#a both.fasta#b | same rows as the two files; each sequence names FILE#NAME |

**Review notes:**

- Only one deletion is tested: large rearrangements are beyond the built-in aligner (the recipe says so).

## `test_codon_batch_recodes_in_place_keeping_the_protein`

**Result:** ✅ passed in 0.04 s · [line 503](../../tests/test_agent_cli.py#L503) · **Area:** Command line for agents · **Checked against:** hand-computed · **Strength:** 🟢 strong

**Why it matters:** An agent could write a recoded plasmid whose protein changed, whose reverse-strand gene was written the wrong way round, or whose bases differ from what the Recode window gives.

With the Highest frequency method each codon of the planted CDS becomes the most used synonym of the table (hand-read from the CSV: ATG GCG GAA AAA TTT GGC TAA); the plasmid has only that stretch changed, for a reverse-strand CDS the reverse complement is written; with a seed the frequency-matched result equals codon_tools.recode and keeps the protein; a bad name is that row's error; GC limits that contradict are refused.

| # | Input | Expected output |
|---|---|---|
| 1 | gx and gy with --method 'Highest frequency' | both rows ATGGCGGAAAAATTTGGCTAA, 2 codons changed; fwd.gb = seq[:100] + best + seq[121:]; rev.gb = seq[:400] + rc(best) + seq[421:] |
| 2 | a row naming a missing feature | that row 'No feature called', exit 1 |
| 3 | --gc-low 70 | ERROR '--gc-low is above' |
| 4 | codon check ... --no-gc | settings show gc_limits 'none' and no gc_low / gc_high / gc_window: they say what was applied |
| 5 | --gc-high 20 --tries 3; --tries 99 | seed_used in 1..3 and the warning says 'above 20 %'; settings show one usage table; recode_tries keeps the fewest unmet; ERROR for --tries 99 |

**Review notes:**

- The usage table is an approximate E. coli table in the examples: the test checks the logic against the table given, not that the table is accurate.

## `test_export_command_writes_the_formats_and_never_overwrites`

**Result:** ✅ passed in 0.04 s · [line 539](../../tests/test_agent_cli.py#L539) · **Area:** Command line for agents · **Checked against:** independent implementation · **Strength:** 🟡 adequate

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

**Result:** ✅ passed in 0.03 s · [line 555](../../tests/test_agent_cli.py#L555) · **Area:** Command line for agents · **Checked against:** hand-computed · **Strength:** 🟢 strong

**Why it matters:** An agent could order a primer with a typo, or write half an order.

Two primers become an IDT plate file (wells A1, B1, scale and purification as asked, sequences upper-cased) or a plain name,sequence file; a table with a primer that is too short or has a bad base writes nothing and names both rows; an existing file is refused.

| # | Input | Expected output |
|---|---|---|
| 1 | --scale 100nm --purification PAGE | Well Position,Name,Sequence,Scale,Purification then A1 F1 ... B1 R1 ... |
| 2 | --format plain | F1,ATGC... / R1,TTAG... |
| 3 | a table with 'short' (4 bases) and 'bad' (an X) | ERROR naming both, no file |

## `test_variants_command_finds_what_the_simulated_reads_were_given`

**Result:** ✅ passed in 0.04 s · [line 572](../../tests/test_agent_cli.py#L572) · **Area:** Command line for agents · **Checked against:** independent implementation · **Strength:** 🟢 strong

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

## `test_tool_settings_override_limits_and_checks`

**Result:** ✅ passed in 0.03 s · [line 597](../../tests/test_agent_cli.py#L597) · **Area:** Settings · **Checked against:** hand-computed · **Strength:** 🟢 strong

**Why it matters:** A value given on the command line could stay in force after its command; a hand-edited settings file could put a setting beyond its limits or contradict another (more exact 3' bases than bases compared); or a missing genome file could be accepted.

Nested overrides; values beyond their limits read within them; the Primer binding check brings settings that limit each other in range and says so, and refuses a missing genome file.

| # | Input | Expected output |
|---|---|---|
| 1 | override 41, inside it override 42 | value 42 inside, 41 after the inner block |
| 2 | both blocks left | no override left; the saved value is read (not 41) |
| 3 | settings.ini: seed 99, exact 3' 40, mismatches 7 | binding() == (30, 30, 0) |
| 4 | check {seed 15, exact 20, mismatches 5}; check of consistent values; genome 'no_such_genome.fa' switched on | exact 15, mismatches 0, note 'exact 3′ bases lowered to 15'; ''; UserError 'no file' |

**Review notes:**

- Single-threaded only: the overrides are process-wide (the command line runs one command at a time).
- Writes the test data folder's settings.ini and removes the two tools' keys afterwards.

## `test_checkup_command_gives_the_verdicts_of_the_checkup_window`

**Result:** ✅ passed in 0.03 s · [line 629](../../tests/test_agent_cli.py#L629) · **Area:** Command line for agents · **Checked against:** independent implementation · **Strength:** 🟡 adequate

**Why it matters:** An agent could say a gene will be expressed when nothing drives it, or the reverse.

On a plasmid with two CDS and nothing else both are 'Not expressed' (no bacterial promoter upstream); with an annotated promoter, RBS and terminator around the forward one it is 'Expressed' with those elements named and the other stays 'Not expressed'; the verdicts equal checkup.run on the same document; a threshold beyond the window's range warns. An annotated p15A origin is reported as fields (label, positions, name, group), not as a Python object's text.

| # | Input | Expected output |
|---|---|---|
| 1 | checkup t.gb --no-predict --no-library | gx and gy Not expressed, reason 'no bacterial promoter upstream on this strand' |
| 2 | with promoter 40..70, RBS 88..99, terminator 130..160 annotated | gx Expressed, promoter Pgx, terminator Tgx; elements all 'annotated'; gy Not expressed |
| 3 | --min-promoter 99 | warning 'responsibility' |
| 4 | a rep_origin 200..400 labelled p15A | origins[0] has label p15A, start 200, end 400, name p15A, incompatibility_group p15A |

**Review notes:**

- The verdict logic belongs to the checkup module and its own tests; this test only shows the command reports it and leaves the user's saved settings alone.

## `test_annotate_command_finds_library_parts_and_adds_only_new_ones`

**Result:** ✅ passed in 0.09 s · [line 656](../../tests/test_agent_cli.py#L656) · **Area:** Command line for agents · **Checked against:** independent implementation · **Strength:** 🟡 adequate

**Why it matters:** An agent could annotate a plasmid twice, or with parts that are not in it.

The longest bundled terminator planted into a random 600 bp sequence is found at the planted place with identity 1.0 and not yet annotated; the hits equal annotation.detect; --add writes a copy with the feature, and a second pass over that copy adds nothing because the part is now annotated. A promoter annotated the other way is reported (annotated_other_way) and --add turns it round instead of adding a second one.

| # | Input | Expected output |
|---|---|---|
| 1 | annotate p.gb --kinds terminator | a hit starting at base 101, strand +, identity 1.0, already_annotated false |
| 2 | --add annotated.gb | the copy has a terminator feature at the same place |
| 3 | annotate annotated.gb --add again.gb | 0 features added |
| 4 | --min-dna 10 | warning 'responsibility' |
| 5 | the longest library promoter planted at 101, annotated complement(...) as 'mine' | hit at 101: annotated_other_way, not already_annotated, strand 1 |
| 6 | --add fixed.gb | turned_round >= 1; one 'mine' feature, strand +1, same bases; a second pass: already annotated, nothing to turn |

**Review notes:**

- The part used is whichever bundled terminator is longest: if the library changes, the test follows it.

## `test_annotate_command_lists_plannotate_hits_beside_the_library_and_says_who_found_each`

**Result:** ✅ passed in 0.04 s · [line 702](../../tests/test_agent_cli.py#L702) · **Area:** Command line for agents · **Checked against:** in-test reference model · **Strength:** 🟡 adequate

**Why it matters:** A user could be shown pLannotate's features as openDNA's own, lose one of the two tools' hits in the merge, or get a feature added with no trace of who found it.

Canned pLannotate rows stand in for the search: one on exactly the planted terminator, one elsewhere (a protein-database fragment on the other strand). The merged list keeps both (nothing dropped or doubled), marks each found_by, names the Features Database's part for the shared one, and says the database of each.

| # | Input | Expected output |
|---|---|---|
| 1 | annotate p.gb --kinds terminator,CDS (pLannotate rows: rrnB T1 on the planted bases, AmpR fragment at 301..420 on the minus strand) | 2 hits: 'openDNA + pLannotate' (library name, plannotate = rrnB T1 / snapgene / 0.99) and 'pLannotate' (AmpR (fragment), source 'pLannotate: swissprot', protein_match true); plannotate_citation holds the paper's DOI |
| 2 | --no-plannotate | plannotate 'off', every hit found_by openDNA, the setting reported False |
| 3 | --add annotated.gb | the fragment's /note starts 'detected by pLannotate (database swissprot): 90 % identity' and carries pLannotate's description; the shared feature's note says it was also found by pLannotate as rrnB T1 |
| 4 | pLannotate not installed | plannotate 'missing', a warning that names `wsltools install`, only openDNA hits |

**Review notes:**

- pLannotate itself is not run (its rows are canned): the real search was tried by hand on demo/pBAD30_eilR.gb, 21 rows, 3 of the 7 library hits found by both. The rows' format is checked in test_plannotate_rows_become_positions_the_way_the_windows_draw_them.

## `test_detect_features_merges_only_the_same_thing_in_the_same_place`

**Result:** ✅ passed in 0.00 s · [line 733](../../tests/test_agent_cli.py#L733) · **Area:** Detect features · **Checked against:** hand-computed · **Strength:** 🟢 strong

**Why it matters:** Two different features could be fused into one row (a gene and its antisense neighbour), or one thing listed twice, in the merged list of the two tools.

annotation.merge on hand-made hits: a terminator on the other strand with the same bases merges; a CDS on the other strand does not; a misc_feature on a promoter's bases merges; a feature elsewhere stays alone; the list is in order of position.

| # | Input | Expected output |
|---|---|---|
| 1 | terminator 100..200 (+) and 105..200 (-) | one entry, found by both, `also` the other tool's hit, source 'openDNA + pLannotate' |
| 2 | CDS 300..600 on + and on - | two entries |
| 3 | promoter 700..760 and misc_feature 700..760 | one entry, found by both |
| 4 | terminator 1000..1050 from pLannotate only | kept, found by pLannotate |

## `test_plannotate_rows_become_positions_the_way_the_windows_draw_them`

**Result:** ✅ passed in 0.00 s · [line 747](../../tests/test_agent_cli.py#L747) · **Area:** Detect features · **Checked against:** hand-computed · **Strength:** 🟡 adequate

**Why it matters:** A pLannotate feature across the origin, or on the minus strand, could be drawn in the wrong place.

plannotate_run.rows_from on two rows in pLannotate's columns: a plain CDS, and an origin of replication that starts at 5900 and ends at 40 on a 5936 bp circle on the minus strand, 99.34 % identity, a fragment.

| # | Input | Expected output |
|---|---|---|
| 1 | AmpR qstart 520, qend 1381, sframe 1, pident 100 | start 520, end 1381, strand 1, kind CDS, identity 1.0 |
| 2 | origin of replication qstart 5900, qend 40, sframe -1, pident 99.34, fragment | start 5900, end 5976 (past the length, like a library hit across the origin), strand -1, kind rep_origin, identity 0.9934, fragment true |

**Review notes:**

- Only the conversion is tested; the column names are those of pLannotate 2.0.0 (checked by running it), and a later pLannotate could rename them.

## `test_wsltools_command_says_what_is_missing_and_never_installs_wsl`

**Result:** ✅ passed in 0.02 s · [line 756](../../tests/test_agent_cli.py#L756) · **Area:** Command line for agents · **Checked against:** property/invariant · **Strength:** 🔴 weak

**Why it matters:** The installer's optional step, or an agent, could try to install WSL itself, or fail with a traceback where WSL is missing; or the install could leave out pLannotate's databases or a Terminal tool; or the status (also shown by Help > Run Setup Again) could report the wrong tools as installed.

With every call into WSL failing, `wsltools status` reports wsl false and the Terminal's tools, and `wsltools install` is one ERROR line that tells the user to run `wsl --install` themselves. With WSL working, one WSL start answers both the pLannotate status and the installed tools (the list after the mark), each tool true only when listed. The install script names the pLannotate environment, runs `plannotate setupdb`, lists every tool the pipelines need, ends with the combined exit status; pLannotate alone leaves the Terminal's environment out.

| # | Input | Expected output |
|---|---|---|
| 1 | wsltools status (no WSL) | wsl false, plannotate no_wsl, one entry per tool of tools.toml |
| 2 | wsltools install (no WSL) | exit 1, 'ERROR: ... wsl --install' |
| 3 | wsltools status, WSL answering 'ready 2.0.0', the mark, 'fastp', 'samtools' | one WSL call; wsl true, plannotate ready 2.0.0; fastp and samtools true, every other tool false |
| 4 | install_script() | has opendna-plannotate, setupdb, every needed tool, ends 'exit $r' |
| 5 | install_script(terminal=False) | no fastp, no opendna-tools |

**Review notes:**

- The script text is checked, not run: running it downloads hundreds of MB. It was run for real by hand (pLannotate 2.0.0 and its databases installed, `wsltools status` ready) but nothing in the suite does.

## `test_detect_features_calls_a_part_known_only_when_a_feature_covers_exactly_its_bases`

**Result:** ✅ passed in 0.00 s · [line 775](../../tests/test_agent_cli.py#L775) · **Area:** Command line for agents · **Checked against:** hand-computed · **Strength:** 🟡 adequate

**Why it matters:** A part could be reported as already annotated just because another feature overlaps it, hiding a real missing part.

A stub editor with one terminator feature: the same bases (either strand for a terminator) count as known; a feature shifted by 5 bp or larger than the part does not, but still counts as covered for the Plasmid Checkup; a CDS on the other strand does not count. Also: the same type on the same bases in the other direction is found (reversed_twin), and turned() turns a location round.

| # | Input | Expected output |
|---|---|---|
| 1 | feature = the hit's bases, strand +1 or -1, a terminator | known |
| 2 | feature shifted by 5 bp, or 20 bp larger each side | not known; covered (exact=False) for the shifted one |
| 3 | a CDS hit, feature on the other strand | not known |
| 4 | a promoter -1 on the hit's bases, hit +1 | reversed_twin finds it; none for the same direction or another type |
| 5 | 101..160 / complement(101..160) / join(complement(5..9),complement(1..3)) | complement(101..160) / 101..160 / join(1..3,5..9) |

**Review notes:**

- Uses a stand-in editor, not a real one; Detect Features' window itself has no widget test.

## `test_gibson_command_adds_overlaps_and_builds_the_product_the_panel_builds`

**Result:** ✅ passed in 0.04 s · [line 791](../../tests/test_agent_cli.py#L791) · **Area:** Command line for agents · **Checked against:** independent implementation · **Strength:** 🟡 adequate

**Why it matters:** An agent could give primers whose overlaps do not join the fragments, or a product that is not what the Assembly panel makes.

Two unrelated 300 bp fragments (the second cut to 11..290) become a 580 bp circle: the product, turned to start at the first fragment, is fragment 1 + the cut fragment 2; the primer of fragment 2 carries the end of fragment 1 as its 5' tail and the reverse primer of fragment 1 the start of fragment 2, the two tails adding to 20 to 40 bases; primers, check levels and product equal auto_primers / analyse / simulate; a bad row, a duplicate id and an existing folder are refused, an overlap beyond the window's range warns.

| # | Input | Expected output |
|---|---|---|
| 1 | gibson run f.csv | 580 bp circular product containing fa + fb[10:290]; fb-F tail is the end of fa; fa-R tail is the reverse complement of the start of fb; 2 junctions |
| 2 | region 0..5 | ERROR 'positions must satisfy' |
| 3 | two rows with the same id | ERROR 'must be different' |
| 4 | --min-overlap 5 | warning 'responsibility' |
| 5 | primers | tail_nt equals the non-annealing 5' part of each primer |

**Review notes:**

- Only Gibson is covered: Golden Gate, Sidewinder and restriction/ligation have no command yet.
- The product's origin is the window's (a few bases before fragment 1): the test compares it as a circle.

## `test_ncbi_command_searches_and_fetches_through_the_import_dialogs_functions`

**Result:** ✅ passed in 0.03 s · [line 828](../../tests/test_agent_cli.py#L828) · **Area:** Command line for agents · **Checked against:** hand-computed · **Strength:** 🟡 adequate

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

## `test_hub_command_lists_and_copies_like_the_import`

**Result:** ✅ passed in 0.10 s · [line 849](../../tests/test_agent_cli.py#L849) · **Area:** Agent command line · **Checked against:** hand-computed · **Strength:** 🟢 strong

**Why it matters:** Without it, `hub list` / `hub get` could disagree with Import > Track Hub, or a hub's groups, labels, colours and hidden tracks could be read wrong (or a track silently dropped), with nothing to show it.

Lists and copies the demo hub and compares with track_hub (the module the window calls); then a hub built in the test with two genomes, an included trackDb file, a composite with a view and a subtrack switched off, a colour given by the composite, a super track not shown, a bam track and a track hidden by its own visibility; then the same hub on a server (track_hub.open_url replaced: no network).

| # | Input | Expected output |
|---|---|---|
| 1 | hub list demo/hub | plus (Coverage, shown, #c82828), minus (Coverage, hidden, #2828c8), sites (Demo hub, shown, no colour); the same as track_hub.hub_tracks |
| 2 | hub get demo/hub --out copy; again into copy; get with no --out | plus.bw, minus.bw, sites.bb with the same bytes; 'not empty'; 'needs --out' |
| 3 | the two-genome hub without --genome; with --genome g1 | 'several genomes'; a1 (Comp, shown, #010203), a2 (hidden: parent off), s1 (Super, hidden: super track not shown), reads (bam, not read, said in warnings), quiet (hidden) |
| 4 | hub get https://example.org/h/hub.txt --genome g1 (open_url serves the files) | a1.bw, a2.bw, s1.bb, q.bw with their bytes; ../a2.bw from g1/ asked at the hub's folder; r.bam not copied; settings say it sends requests to its server |

**Review notes:**

- The window's dialog and its copy with a progress bar (hub_dialog, main.import_hub) are not tested: no widget tests; checked by running the app offscreen on the demo hub.

## `test_recap_command_gives_the_numbers_of_the_window_recap`

**Result:** ✅ passed in 0.06 s · [line 908](../../tests/test_agent_cli.py#L908) · **Area:** Agent command line · **Checked against:** hand-computed · **Strength:** 🟢 strong

**Why it matters:** Without it, `recap` could disagree with the window's Recap Tracks Here, or a recap could give a signal's mean, sum or items wrong, mix up a group's mean row, or miss a project's highlighted regions, with nothing to show it.

Recaps the demo plasmid with the demo hub's two bigWigs and bigBed and compares the JSON with track_recap.recap on the same tracks loaded the same way (cli_tools.track_editor); checks the numbers by hand, the per-base values across a step, and a project saved with a group and a highlight recapped with no place given, as a table.

| # | Input | Expected output |
|---|---|---|
| 1 | recap pBAD30_eilR.gb 1001..2600 --track plus.bw minus.bw sites.bb | = track_recap.recap; plus: 1600 bases, mean 17.5, min 10, max 30, sum 28,000; minus mean 5, sum 8,000; sites: site_b 2501-2600 '-' covers 0.0625; the sequence's features overlapping |
| 2 | recap 1995..2004 --track plus.bw --per-base | values [10 x 6, 30 x 4] |
| 3 | recap 1..500 --track w.bed (one item 5800..6236 on the 5,936 bp circle) | item wrap 5801-6236 covers 0.6 (300 of 500 bases, past the origin) |
| 4 | a project with plus and minus in group Coverage and a highlight 'here' 1001..2600; recap p.oDNA --tsv | header = COLUMNS; rows sites.bb (1 item), plus.bw (1600, 17.5, 28000), minus.bw (1600, 5.0, 8000), 'Coverage (mean of 2)' (1600, 11.25, 18000) |
| 5 | recap with tracks but no place, no --regions, not a project | error 'no region' |

**Review notes:**

- A BAM file's depth and variants in a recap are not checked here (no small BAM in the demo); the window's dialog around the same rows is checked by running the app only.

## `test_print_command_draws_the_project_as_the_window_export`

**Result:** ✅ passed in 0.16 s · [line 945](../../tests/test_agent_cli.py#L945) · **Area:** Agent command line · **Checked against:** smoke · **Strength:** 🟡 adequate

**Why it matters:** Without it, `print` could draw a project without its tracks' names as shown or its highlights, at the wrong width, or write over a file, with nothing to show it.

Saves a project of the demo plasmid with a renamed signal track and a highlight, prints a region to SVG and checks the text drawn; prints a sequence with --track to PNG at a chosen width; refuses an existing file and an unknown extension.

| # | Input | Expected output |
|---|---|---|
| 1 | print p.oDNA 1001..4000 r.svg (plus.bw shown as 'Coverage plus', sites.bb, highlight 'Highlighted here') | tracks_shown 2, regions 1, start 1001, end 4000; the SVG holds 'Coverage' (the name wraps in the narrow names column), 'Highlighted here', 'sites' |
| 2 | print pBAD30_eilR.gb 1001..4000 r.png --track plus.bw --width 900 | a PNG 900 px wide |
| 3 | print again to r.svg; print to r.jpg | 'already exists'; 'must end in' |

**Review notes:**

- Checks the text and size of the image, not the drawing itself (the bars, heatmaps and bands are checked by looking at images from the app).

## `test_track_features_command_adds_what_copy_to_features_adds`

**Result:** ✅ passed in 0.04 s · [line 966](../../tests/test_agent_cli.py#L966) · **Area:** Agent command line · **Checked against:** hand-computed · **Strength:** 🟢 strong

**Why it matters:** Without it, copying a track's items into the sequence could lose an item's name, its 'from track' note or its file's extra columns, add the same item twice, or give the wrong type, with nothing to show it.

Copies a BED file of two items (one with extra columns, one on the minus strand) into the demo plasmid with track_features (what Copy to Features adds), reads the written GenBank file back, runs it again on the result, and copies only a part of the demo bigBed.

| # | Input | Expected output |
|---|---|---|
| 1 | track-features pBAD30_eilR.gb Sites_CDS.bed --name-in gene | added 2, skipped 0, type CDS (from the name); CDS 101..400 gene=site_a, notes 'from track Sites_CDS.bed', 'height=81', 'RSPM=12'; CDS complement(2501..2600) gene=site_b |
| 2 | default_type('RefSeq_2022_CDS.bed'), default_type('sites.bb') | CDS, misc_feature |
| 3 | the same on the written file | added 0, skipped 2 |
| 4 | sites.bb over 1..1000; a bigWig as the track | added 1; 'not an annotation track' |

## `test_signal_regions_command_finds_the_stretches_of_regions_above_threshold`

**Result:** ✅ passed in 0.02 s · [line 988](../../tests/test_agent_cli.py#L988) · **Area:** Agent command line · **Checked against:** hand-computed · **Strength:** 🟢 strong

**Why it matters:** Without it, Regions Above Threshold and `signal-regions` could miss, split or merge stretches wrongly, give a wrong mean, or write a wrong BED file, with nothing to show it.

Finds the stretches of the demo bigWig over 20 (with its BED file), then of a bedGraph written in the test with joining, a least length and a region, all by hand.

| # | Input | Expected output |
|---|---|---|
| 1 | signal-regions plus.bw --threshold 20 --out r.bed | one stretch 2001-5936 (3936 bp), mean 98,720 / 3936, max 30; BED 'pBAD30_eilR 2000 5936 plus_1 <mean> .' |
| 2 | bedGraph 5 over [0,10), 1 over [10,15), 5 over [15,20), 9 over [100,103); threshold 5 | 1-10, 16-20, 101-103 |
| 3 | threshold 5, gap 5, min length 4 | 1-20 only, mean 4 ((10 x 5 + 5 x 1 + 5 x 5) / 20) |
| 4 | --region 12..103; --region 0..5 | 16-20, 101-103; error 'START..END' |

## `test_split_command_writes_one_file_per_record_and_file_hash_name_reads_one`

**Result:** ✅ passed in 0.06 s · [line 1026](../../tests/test_agent_cli.py#L1026) · **Area:** Command line for agents · **Checked against:** round-trip/self-consistency · **Strength:** 🟡 adequate

**Why it matters:** An agent could be unable to use one chromosome of a multi-record file, or split a file into files that lose bases or have clashing names.

A two-record GenBank file is refused by a bare name with the record names and the way out; FILE#NAME reads one record; split writes one file per record with the right lengths and topology and refuses a used folder or a one-record file; a multi-FASTA with two names differing only by case gives two FASTA files and a warning.

| # | Input | Expected output |
|---|---|---|
| 1 | inspect multi.gb | ERROR naming p_one, p_two and `split` |
| 2 | region multi.gb#p_two 1 10; #nope | the first 10 bases of that record; ERROR 'no record called' |
| 3 | split multi.gb --out parts | p_one.gb (120, circular), p_two.gb (90, linear); p_two.gb has the original bases |
| 4 | split again; split a one-record file | ERROR 'not empty'; ERROR 'one record' |
| 5 | split m.fasta (>x, >X) | x.fasta and X_2.fasta, a warning, X_2 holds the second sequence |

**Review notes:**

- Sequences are random and short; a real GenBank file with features is only covered by the demo recipes' examples.

## `test_merge_command_joins_files_and_records_and_keeps_what_the_format_keeps`

**Result:** ✅ passed in 0.08 s · [line 1049](../../tests/test_agent_cli.py#L1049) · **Area:** Command line for agents · **Checked against:** round-trip/self-consistency · **Strength:** 🟡 adequate

**Why it matters:** An agent could merge files into a multi-record file that loses a sequence, overwrites a file, or silently drops features or topology.

Three GenBank inputs (one repeated) are joined in order into a multi-record file with the repeated name numbered and a warning; every record reads back with its bases; an existing output and a single sequence are refused; a record picked by FILE#NAME joins a .fasta output, which warns that only the bases are kept. A record whose name is taken is named after its file; --names sets every name, and a wrong count is refused.

| # | Input | Expected output |
|---|---|---|
| 1 | merge a.gb b.gb a.gb --out all.gb | records pa (120, circular), pb (90), a (120, circular: the repeated file is named after itself); a warning |
| 2 | region all.gb#pb 1 90 and #pa_2 1 120 | the original bases |
| 3 | same --out again; one input | ERROR 'already exists'; ERROR 'two or more' |
| 4 | merge all.gb#pb a.gb --out two.fasta | records pb, pa; warning 'FASTA keeps only the bases'; two.fasta#pa starts with a's bases |
| 5 | merge a.gb (record pa) other.gb (record pa) | records pa and other; a warning |
| 6 | --names x,y; --names x | x and y, no warning; ERROR 'one name each' |

**Review notes:**

- Sequences are random and short; a real GenBank file with features is only covered by the demo recipes' examples.

## `test_circular_tools_give_the_same_answers_on_a_plasmid_turned_round_the_origin`

**Result:** ✅ passed in 0.32 s · [line 1072](../../tests/test_agent_cli.py#L1072) · **Area:** Command line for agents · **Checked against:** property/invariant · **Strength:** 🟢 strong

**Why it matters:** A restriction site, guide, repeat or gene that straddles the origin of a plasmid could be missed or misplaced by a tool, while the same plasmid cut elsewhere is fine.

A 900 bp circle holds EcoRI, a guide with its PAM, a 30 bp repeat in two copies and a CDS. It is turned so that each piece crosses the origin in turn, and sites, guides, repeats, DNA and protein search, feature, pcr and kld are run on every turned copy: they must give what they give on the original, with positions mapped back. codon also recodes the CDS across the origin: the same recoded bases on every turned copy; the written plasmid keeps the protein, the gene's place and every base outside the gene. A residue split by the origin is an error, not a change over the whole plasmid.

| # | Input | Expected output |
|---|---|---|
| 1 | sites, search GAATTC on the original and on copies turned by 302 and 110 | the same site positions (mapped back), the same guide with its strand |
| 2 | repeats --min-len 25 on copies turned by 64 and 614 (each copy of the repeat crosses the origin) | the same one direct repeat of about 30 bp |
| 3 | search the protein MAKDAHHHF, feature wrapCDS, pcr target wrapCDS and kld p.Ala2Gly on the copy turned by 720 | same protein hit, same 60 bp feature and protein, product 60 bp with the same first primer, p.(Ala2Gly) missense |
| 4 | region START END with START > END on the turned copy; inspect | the CDS bases; the feature is marked crosses_origin |
| 5 | codon check on every turned copy; codon run on the copy turned by 720 | same codons changed and same recoded sequence; only bases of the CDS (880..900, 1..39) change; protein unchanged; feature still join(880..900,1..39) |
| 6 | kld wrapCDS:p.His8Ala on the copy turned by 721 (codon 8 is split by the origin) | row error 'cross the origin' |

**Review notes:**

- Variants have their own test (they need a BAM).
- Random sequence: a stretch of 30 bases that repeats by chance next to the planted one would change a length (the check allows 30 to 32).

## `test_a_limit_that_is_reached_is_always_a_warning`

**Result:** ✅ passed in 0.05 s · [line 1132](../../tests/test_agent_cli.py#L1132) · **Area:** Command line for agents · **Checked against:** in-test reference model · **Strength:** 🟡 adequate

**Why it matters:** A list cut short or a search stopped at its limit could be taken for the complete answer.

Each limit is lowered to 1 (guides --top 1, search hits, repeats listed, align pieces, NCBI --limit with a canned total of 40): the answer of every command carries a warning that says what was left out. Bases that are N are counted in the guides warning.

| # | Input | Expected output |
|---|---|---|
| 1 | guides --top 1 | 'Only the best 1 of the ... guides' |
| 2 | search with the hit limit 1; repeats with the list limit 1; align with 1 piece | 'stopped at'; 'only the longest 1 of'; 'Only 1 pieces are listed' |
| 3 | ncbi search --limit 1 with a canned count of 40 | '40 matches; only the first 1 are shown' |
| 4 | guides on a sequence with NNNN | warning '4 base(s) are N ...' |

**Review notes:**

- The limits are monkeypatched down, so it shows each command reports a reached limit, not that the real limits (500, 20,000, 30) are the right ones.

## `test_variants_command_sees_changes_on_both_sides_of_the_origin_and_apply_writes_them`

**Result:** ✅ passed in 0.03 s · [line 1154](../../tests/test_agent_cli.py#L1154) · **Area:** Command line for agents · **Checked against:** hand-computed · **Strength:** 🟢 strong

**Why it matters:** A change at the first or last base of a plasmid, or a deletion across the origin, could be missed or mistaken for a break, because reads over the origin are mapped as two clipped pieces.

A 600 bp circle with SNPs at bases 1, 5 and 600 is read by 100-base reads every 10 bases, those over the origin mapped as a primary piece and a supplementary piece (written with the test suite's own BAM writer): the three SNPs are consensus with 10 of 10 reads, nothing else is reported, and --apply writes a circular copy equal to the sample. Reads of a circle lacking 80 bases across the origin show breaks at 41 and 561 only. --apply writes the deletion through the origin; reads hard-clipped at the origin are reported with a note.

| # | Input | Expected output |
|---|---|---|
| 1 | variants on the SNP sample, --min-depth 3 --apply fixed.gb | snps at 1, 5, 600 (alt as planted), 10 of 10 reads, consensus; fixed.gb equals the sample and is circular |
| 2 | variants on reads of the circle without reference bases 561..600 and 1..40, --apply | one consensus deletion of 80 bases at 561 (through the origin); the written plasmid equals the shortened sequence |
| 3 | reads of the unchanged circle whose supplementary pieces are hard-clipped at the origin | one clip at 1 whose text says it could not be checked |

**Review notes:**

- The reads are error-free and evenly spaced; real reads would add noise (the thresholds are tested in test_core_formats).
