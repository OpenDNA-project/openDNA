# `tests/test_agent_cli.py`

[← all test files](README.md) · source: [`tests/test_agent_cli.py`](../../tests/test_agent_cli.py)

## `test_agent_reads_the_plasmid`

**Result:** ✅ passed in 0.01 s · [line 29](../../tests/test_agent_cli.py#L29) · **Area:** Command line for agents · **Checked against:** independent implementation · **Strength:** 🟡 adequate

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

**Result:** ✅ passed in 0.01 s · [line 46](../../tests/test_agent_cli.py#L46) · **Area:** Command line for agents · **Checked against:** hand-computed · **Strength:** 🟡 adequate

**Why it matters:** A mistake in the agent's input could end as a traceback or silent output that the agent cannot read or pass to the user.

Each user mistake gives exit code 1, nothing on stdout and exactly one `ERROR:` line on stderr that says what is wrong.

| # | Input | Expected output |
|---|---|---|
| 1 | missing file; unknown feature; region 0..12; region past the end; unknown recipe | exit 1, empty stdout, one ERROR line naming the problem |

## `test_every_recipe_example_still_runs`

**Result:** ✅ passed in 0.01 s · [line 54](../../tests/test_agent_cli.py#L54) · **Area:** Command line for agents · **Checked against:** smoke · **Strength:** 🟡 adequate

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

**Review notes:**

- Features made of several parts (joins) and changes beyond the ends of a CDS are not exercised.

## `test_hgvs_mistakes_say_what_to_fix_and_ambiguities_only_warn`

**Result:** ✅ passed in 0.00 s · [line 118](../../tests/test_agent_cli.py#L118) · **Area:** Command line for agents · **Checked against:** hand-computed · **Strength:** 🟡 adequate

**Why it matters:** A wrong reference base or amino acid in a mutation table could be applied silently, and an ambiguous codon choice could be made without telling the user.

Every mistake (wrong reference, wrong residue, nothing to change, outside the plasmid, past the CDS end, unknown feature, non-adjacent insertion, position 0, missing feature name, unreadable text, c.*N) raises one readable error; a unique nearest codon gives no warning, tied codons give a default and a warning naming the alternatives.

| # | Input | Expected output |
|---|---|---|
| 1 | gx:c.4C>A, p.Lys2Thr, p.Ala2Ala, o.600A>T, c.50G>A, nope:c.1A>T, c.9_11insA, c.0A>T, c.4G>A, 'lacZ del 5', c.*3A>T | an error each, with the words 'not C', 'is Ala (GCT), not Lys', 'already', 'outside the plasmid', 'past the end', 'No feature called', 'neighbouring', 'no position 0', 'name the feature', 'Cannot read', 'stop codon' |
| 2 | gx:p.Ala2Gly; p.Ala2Ser; p.Glu3Asp | G with no warning; T; C with a warning naming GAC and GAT |

**Review notes:**

- The warning for a change that can slide along a repeat is not tested.

## `test_kld_batch_gives_the_files_and_primers_the_panel_gives`

**Result:** ✅ passed in 0.06 s · [line 133](../../tests/test_agent_cli.py#L133) · **Area:** Command line for agents · **Checked against:** independent implementation · **Strength:** 🟢 strong

**Why it matters:** The batch could give files or primers different from what the user gets by hand in the KLD panel, or overwrite earlier results, or accept a malformed table.

A 4-row table (substitution, reverse-strand substitution, deletion, wrong reference) writes the 3 good edited plasmids with both primers drawn on them, a primers table, a report and the settings; the plasmids are the string edits of the original; the primers are those of kld.design_selected; the failed row is reported and gives exit code 1; a second run into the same folder, extra columns, repeated or bad ids are refused.

| # | Input | Expected output |
|---|---|---|
| 1 | kld run t.gb m.csv --out out (rows sub, rev, del, bad) | exit 1; summary 3 ok 1 error; out holds del.gb primers.csv report.jsonl rev.gb settings.json sub.gb |
| 2 | each written .gb | bases = original with the change; primer_bind features KLD-F and KLD-R; primers equal kld.design_selected's |
| 3 | same run again; table with a note column; ids 'a b' and a repeated id | 'already exists'; 'exactly id,mutation'; 'used twice' |

**Review notes:**

- Settings are the defaults of the test machine: a saved setting that changes the primers is not exercised.

## `test_pcr_batch_gives_the_pairs_the_panel_gives`

**Result:** ✅ passed in 1.74 s · [line 170](../../tests/test_agent_cli.py#L170) · **Area:** Command line for agents · **Checked against:** independent implementation · **Strength:** 🟢 strong

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

**Result:** ✅ passed in 0.04 s · [line 204](../../tests/test_agent_cli.py#L204) · **Area:** Command line for agents · **Checked against:** hand-computed · **Strength:** 🟢 strong

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
