# `tests/test_agent_cli.py`

[← all test files](README.md) · source: [`tests/test_agent_cli.py`](../../tests/test_agent_cli.py)

## `test_agent_reads_the_plasmid`

**Result:** ✅ passed in 0.01 s · [line 37](../../tests/test_agent_cli.py#L37) · **Area:** Command line for agents · **Checked against:** independent implementation · **Strength:** 🟡 adequate

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

**Result:** ✅ passed in 0.01 s · [line 54](../../tests/test_agent_cli.py#L54) · **Area:** Command line for agents · **Checked against:** hand-computed · **Strength:** 🟡 adequate

**Why it matters:** A mistake in the agent's input could end as a traceback or silent output that the agent cannot read or pass to the user.

Each user mistake gives exit code 1, nothing on stdout and exactly one `ERROR:` line on stderr that says what is wrong.

| # | Input | Expected output |
|---|---|---|
| 1 | missing file; unknown feature; region 0..12; region past the end; unknown recipe | exit 1, empty stdout, one ERROR line naming the problem |

## `test_every_recipe_example_still_runs`

**Result:** ✅ passed in 0.01 s · [line 62](../../tests/test_agent_cli.py#L62) · **Area:** Command line for agents · **Checked against:** smoke · **Strength:** 🟡 adequate

**Why it matters:** A recipe could keep telling agents to run a command that no longer exists or whose arguments changed.

Every command in the Example block of every recipe (and the index) runs and succeeds, so the recipes cannot drift from the tool.

| # | Input | Expected output |
|---|---|---|
| 1 | skills list, then each recipe's Example block | index and inspect are listed; every example command exits 0 |

**Review notes:**

- Checks that examples run, not that their output is what the recipe says.

## `test_hgvs_mutations_become_the_change_they_name`

**Result:** ✅ passed in 0.00 s · [line 99](../../tests/test_agent_cli.py#L99) · **Area:** Command line for agents · **Checked against:** hand-computed · **Strength:** 🟢 strong

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

**Result:** ✅ passed in 0.00 s · [line 122](../../tests/test_agent_cli.py#L122) · **Area:** Command line for agents · **Checked against:** hand-computed · **Strength:** 🟡 adequate

**Why it matters:** A wrong reference base or amino acid in a mutation table could be applied silently, and an ambiguous codon choice could be made without telling the user.

Every mistake (wrong reference, wrong residue, nothing to change, outside the plasmid, past the CDS end, unknown feature, non-adjacent insertion, position 0, missing feature name, unreadable text, c.*N) raises one readable error; a unique nearest codon gives no warning, tied codons give a default and a warning naming the alternatives.

| # | Input | Expected output |
|---|---|---|
| 1 | gx:c.4C>A, p.Lys2Thr, p.Ala2Ala, o.600A>T, c.50G>A, nope:c.1A>T, c.9_11insA, c.0A>T, c.4G>A, 'lacZ del 5', c.*3A>T | an error each, with the words 'not C', 'is Ala (GCT), not Lys', 'already', 'outside the plasmid', 'past the end', 'No feature called', 'neighbouring', 'no position 0', 'name the feature', 'Cannot read', 'stop codon' |
| 2 | gx:p.Ala2Gly; p.Ala2Ser; p.Glu3Asp | G with no warning; T; C with a warning naming GAC and GAT |

**Review notes:**

- The warning for a change that can slide along a repeat is not tested.

## `test_a_change_through_the_origin_is_made_and_the_file_keeps_its_origin`

**Result:** ✅ passed in 0.07 s · [line 137](../../tests/test_agent_cli.py#L137) · **Area:** Command line for agents · **Checked against:** hand-computed · **Strength:** 🟢 strong

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

**Result:** ✅ passed in 0.06 s · [line 152](../../tests/test_agent_cli.py#L152) · **Area:** Command line for agents · **Checked against:** independent implementation · **Strength:** 🟢 strong

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

**Result:** ✅ passed in 1.75 s · [line 186](../../tests/test_agent_cli.py#L186) · **Area:** Command line for agents · **Checked against:** independent implementation · **Strength:** 🟢 strong

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

**Result:** ✅ passed in 0.03 s · [line 219](../../tests/test_agent_cli.py#L219) · **Area:** Command line for agents · **Checked against:** hand-computed · **Strength:** 🟢 strong

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

**Result:** ✅ passed in 0.05 s · [line 250](../../tests/test_agent_cli.py#L250) · **Area:** Command line for agents · **Checked against:** hand-computed · **Strength:** 🟢 strong

**Why it matters:** One bad row could end a 200-row run, two rows could overwrite each other's file, or a value the windows would refuse could pass without the user knowing.

A table (saved as cp1252, as Excel does) with ids differing only by case, unsafe characters, a Windows device name, a ragged row, an empty cell, a mutation that changes nothing and a row that fails in an unforeseen way: every other row is processed, every ok row has its own file (unique in any letter case), renames are told in the row's warnings, the errors are loud on stderr (with the traceback of the unforeseen one), a no-op mutation is accepted with a warning and no primers, and a Tm beyond the window's range is used with a responsibility warning.

| # | Input | Expected output |
|---|---|---|
| 1 | ids A, a, 'bad id!', CON, café; rows ragged (3 cells), empty (no mutation), noop (o.104G>G), boom (resolve raises RuntimeError) | files A.gb _CON.gb a_2.gb bad_id_.gb caf_.gb noop.gb; ragged '3 cells but the header has 2', empty 'mutation is empty', boom 'Unexpected RuntimeError' with a traceback; exit 1 |
| 2 | noop | effect 'no change', no primers, 'o.104=', warning 'no change'; its file is the plasmid unchanged |
| 3 | kld check --tm 80 | exit 0; settings max 72.0; top-level warning 'responsibility'; WARNING on stderr |

**Review notes:**

- Only the KLD batch is exercised; pcr and recombineering share the same helpers (run_rows, finish) but are not run with bad rows.
