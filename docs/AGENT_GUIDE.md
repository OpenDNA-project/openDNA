# Using openDNA through an AI agent

*Draft for the GitHub documentation.*

openDNA has a command line with no window, `oDNA-cli`. It runs the same functions as the windows and gives the same results, so an AI agent (Claude Code, or any assistant that can run commands) can design primers, edit plasmids and run analyses for you while you describe what you want in plain words.

You do not need to learn the commands. The agent reads the instructions that ship with openDNA and does the work. This page tells you how to set that up, what the agent can do, and how to get good results.

## 1. Set it up (two minutes)

1. Install openDNA. The Windows folder build contains `oDNA-cli.exe` next to `openDNA.exe`. From a source checkout, the same program is `python -m app.cli` (run from the repository root, in its virtual environment).
2. Give your agent a shell in the folder that holds your files, and start with one sentence:

   > You can run openDNA from the command line: `C:\path\to\oDNA-cli.exe`. Run `skills show index` first, then follow the recipes it points to.

   `skills list` names every recipe and `skills show <name>` prints one, so the agent finds everything it needs from the program itself, whichever version you have.
3. Ask for what you want.

> **Tip.** Put that sentence in your agent's project instructions (for Claude Code, a `CLAUDE.md` file in your working folder), so every session starts knowing how to reach openDNA.

## 2. What you can ask for

| You want to | The agent uses | What you get |
|---|---|---|
| See what is in a GenBank file | `inspect`, `feature`, `region` | length, topology, every feature, the protein of each CDS, exact bases |
| Make point mutations, deletions, insertions or duplications by inverse PCR + KLD | `kld` | one primer pair and one edited plasmid per mutation (written as HGVS, e.g. `lacZ:p.Asp120Ala`) |
| Design PCR primers for several regions | `pcr` | the best primer pairs per region, with the Primers panel's checks |
| Replace a region of a genome with a cassette (lambda Red) | `recombineering` | the two cassette primers with homology tails and the recombinant genome |
| Recode a CDS for a codon usage table | `codon` | recoded plasmids (the protein is kept), GC and codon statistics |
| Assemble fragments by Gibson / HiFi | `gibson` | primers with overlaps, junction checks, the assembled plasmid |
| Find restriction sites, unique cutters, digest fragments | `sites` | enzymes, cut positions, fragment sizes |
| Find Cas9 / Cas12a guides | `guides` | ranked guides with off-target counts and a verdict, as in the Spacer Finder |
| Find repeats shared by two sequences (or one with itself) | `repeats` | direct and inverted repeats |
| Search a sequence, a protein motif or a feature name | `search` | hits with position and strand |
| Locate a sequence (a Sanger read, an insert) in a plasmid | `align` | pieces, strands, mismatches, structural changes |
| Align several sequences | `msa` | alignment, consensus, identity (optional FASTA file) |
| Check a plasmid for expression (**beta**: treat the verdicts as hints) | `checkup` | for each CDS: expressed or not, with promoter, RBS and terminator |
| Annotate known parts | `annotate` | promoters, terminators, origins, resistance genes (optional annotated copy) |
| Compare sequencing reads with the plasmid | `variants` | SNPs, insertions, deletions with read counts (optional corrected copy) |
| Convert to GenBank, FASTA, GFF3, BED, NCBI feature table or protein FASTA | `export` | one new file |
| Turn a primer list into an order file | `order` | an IDT plate file or a plain name,sequence file |
| Get a sequence from NCBI | `ncbi` | search results, or a GenBank / FASTA file |

Example requests:

- "In `pBAD30_eilR.gb`, make T50A in eilR and the deletion of residues 12 to 14 of eilR. Show me the primers before you make the files." (Name the protein: a bare `D120A` could be in any of the plasmid's CDS.)
- "Which restriction enzymes cut `pUC19.gb` exactly once? What are the fragment sizes of an EcoRI + HindIII digest?"
- "Here are three Sanger reads (`r1.ab1` converted to FASTA). Where do they map on the plasmid, and are there mismatches?"
- "Find the best Cas9 guides in the 300 bp around position 1200 and tell me which ones have no off-targets in my background genome."
- "Download pUC19 from NCBI and annotate its parts."

## 3. How a session goes

**Batch functions** (`kld`, `pcr`, `recombineering`, `codon`, `gibson`) work from a table of rows, one row per mutation, region or fragment. The agent:

1. Looks at your plasmid (`inspect`).
2. Turns your list into the exact table the tool reads. If your list is clear it translates it itself; if not, it asks.
3. Runs `check`, which validates every row and writes nothing, and fixes the problems it finds.
4. Shows you the settings (Tm, codon table, PCR conditions, limits) and waits for your OK.
5. Runs `run`, which writes a new folder of files (edited plasmids, a primer table, a report).

**Analyses** (`sites`, `guides`, `repeats`, `search`, `align`, `msa`, `checkup`, `annotate`, `variants`) are read-only: one command, one answer. They write a file only when you ask for one (`--out`, `--add`, `--apply`).

Positions are 1-based and inclusive, as in GenBank files.

## 4. What it guarantees

- **Same results as the windows.** The command line calls the functions the windows call, and the test suite compares each command with its window function.
- **Your files are never overwritten.** Results go in a new folder, or a new file; a name that already exists is refused. If two rows would give the same file name, one is renamed and you are told.
- **One bad row never stops the others.** The failed rows are listed with the reason, and the rest are done.
- **Nothing silent.** When something is ambiguous (several codons for one amino acid, a mutation that changes nothing), the tool picks a default and adds a warning to that row. The agent must pass the warnings on to you.
- **Errors are plain.** A problem you can fix is one line starting `ERROR:`.
- **Your settings are used and shown.** Thresholds you saved in openDNA (Settings) apply, and each answer says where each value came from. A value the agent passes on the command line applies to that run only; it never changes your saved settings.
- **Limits are the windows' limits.** Every adjustable setting shows the range the windows accept. A value outside it is allowed, with a warning that the result is your responsibility.

## 5. Good practices

**Describe the goal and the guard.** Give the mutation with what it should find ("T50A", "360G>C") and **say which protein or gene it is in** ("T50A in eilR"): a bare number could be in any of the plasmid's genes. openDNA checks that the reference base or amino acid matches, which catches a wrong coordinate before any primer is made.

**Let the agent look first.** Ask it to inspect the file and tell you what it sees (the CDS it will use, its strand) before it builds anything. Most mistakes are a wrong feature or a wrong numbering.

**Read the settings it shows you.** The values that change the result most are the annealing Tm and the codon usage table (primers and the codon chosen), the overlap length for Gibson, and the background genome for guides. If a value is not what you want, say so before the run.

**Read the warnings.** They are where the tool tells you what it decided for you. A good agent lists them; if yours does not, ask "were there any warnings?".

**Name things yourself.** The `id` of a row names its output file. Give your own names, or ask for a pattern such as `<gene>_<change>`.

**Check the answer against the file, not the summary.** Ask the agent to open the edited plasmid with `inspect` or `region` and show you the changed bases, or open the file in openDNA. The agent's summary is not a substitute for looking at the result of an experiment you are about to order.

**Treat the output as a design, not a result.** Primers and assemblies are predictions. Check them the way you would one made in the window (the Primers panel's checks, the assembly report) before ordering.

**Keep your data private.** Only `ncbi` uses the network, and it sends only the search words or accession, never your sequence. The agent must tell you the words before the first call. Do not put an unpublished sequence or a project name in a query.

**Start small.** Try the recipe's example first (each recipe in `skills show <name>` ends with an `## Example` that runs on files shipped with openDNA), then your own data.

**One folder per job.** Keep inputs, tables and results of one job together. The agent writes `PLASMID_kld`, `PLASMID_codon` and similar folders next to your file by default; use `--out` to choose another place.

**Ask for a re-run of the failures only.** If three rows out of fifty failed, have the agent fix and re-run those rows rather than the whole table.

## 6. If something goes wrong

| What you see | What to do |
|---|---|
| `ERROR: File not found` | Give the agent the full path, or start it in the folder that holds the file. |
| `ERROR: ... already exists` | Choose a new file or folder name; openDNA never overwrites. |
| `check` reports rows in error | Read the reasons (a reference base that does not match is the usual one) and correct the table; the other rows are fine. |
| A warning that a value is "outside what openDNA's windows accept" | The agent was given a value beyond the windows' range. Confirm it is what you meant. |
| The agent invents a command or a column | Ask it to run `skills show <name>` and follow the recipe. |
| The results differ from the window | This should not happen. Report it with the command, the file and both results (see `docs/FEEDBACK.md`). |

## 7. For other tools and scripts

The same commands work without an agent: every command prints JSON (the `skills` command prints Markdown), errors go to the error stream with exit code 1, and tables are plain CSV. A script can call `oDNA-cli` the way an agent does. `oDNA-cli --help` lists the commands and `oDNA-cli <command> --help` its options.

## 8. Try it on the demo files

Every file below is in the repository's `demo/` folder. Paste a request to your agent (start the agent in the repository folder, or give it the full paths) and compare what it does with section 3. The main plasmid is `pBAD30_eilR.gb`: 5,936 bp, circular, with the genes eilR (5119..5697, 192 amino acids), araC (3655..4533, reverse strand) and Ampicilin (521..1381).

| Function | A request to try |
|---|---|
| `inspect` | What is in `demo/pBAD30_eilR.gb`? List its genes (CDS) with their positions and strand, and give me the first 60 amino acids of eilR and tell me how long that protein is. |
| `kld` | In `demo/pBAD30_eilR.gb`, make two mutations with KLD: T50A in eilR (protein numbering), and the deletion of residues 12 to 14 of eilR. Show me the primers before you make the files. |
| `pcr` | Design PCR primers for two regions of `demo/pBAD30_eilR.gb`: the whole eilR gene, and the p15A origin (positions 2839 to 3384). Show me the best pair for each before you write files. |
| `recombineering` | In the E. coli genome `demo/BW25113_2026.gb`, replace the gene araC with the cassette `demo/demo_insert.gb` (use the whole file, forward orientation) by lambda Red. I only want the primers, not the recombinant genome. Show me the primers before you write files. |
| `codon` | Recode the eilR gene of `demo/pBAD30_eilR.gb` using the codon usage table `demo/codon_usage_ecoli.csv`: no BamHI or SalI sites in the result, GC between 40 and 60 % in every window. Keep the protein identical. If some rule cannot be met, try more random draws and tell me which one you kept. Show me the settings and the GC before and after, then write the plasmid. |
| `gibson` | Using Gibson assembly, build a circular plasmid from three fragments in this order: all of `demo/demo_insert.gb` (1 to 700), then `demo/demo_vector.gb` positions 428 to 1834, then `demo/demo_vector.gb` positions 1 to 307 (this is the vector with its 'dropout' feature replaced by the insert). Show me the overlaps and any problem before you write files. |
| `sites` | In `demo/pBAD30_eilR.gb`: which restriction enzymes cut exactly once, and where? What fragment sizes do I get from a BamHI + SalI digest? |
| `guides` | Find SpCas9 guides in the first 300 bp of the eilR gene (positions 5119 to 5418) of `demo/pBAD30_eilR.gb`. Give me the best 5 with their position, strand, off-target counts and verdict, and tell me which one you would use. |
| `repeats` | In `demo/pBAD30_eilR.gb`, are there direct or inverted repeats of at least 15 bp within the plasmid? Give me the longest three with their positions. Then compare `demo/pBAD30_eilR.gb` with `demo/pBAD30_RC.gb`: what do they share? |
| `search` | In `demo/pBAD30_eilR.gb`, find every BamHI site (GGATCC) on both strands, and every feature whose name contains 'Verif'. |
| `align` | Where does this sequence match `demo/pBAD30_eilR.gb`: `atgggctatctgaatcgcgtagaacgtcgggaaaccatcatgcaggccgc` ? Tell me the positions, the strand, and each mismatch. |
| `msa` | Align `demo/pBAD30_eilR.gb` and `demo/pBAD30_RC.gb`. What is the identity of the second to the first, and what does that tell you about how the two files are related? |
| `annotate` | Which known parts (promoters, terminators, origins, resistance genes) are present in `demo/pBAD30_eilR.gb` but not annotated yet? Write an annotated copy with the missing ones added. |
| `variants` | I sequenced a clone: the reads are in `demo/pBAD30_eilR_clone_reads.bam`, mapped to `demo/pBAD30_eilR.gb`. List the positions where the reads differ from the plasmid, say which ones most reads agree on, and write a corrected copy of the plasmid with those applied. |
| `export` | From `demo/pBAD30_eilR.gb`, give me three files: the sequence as FASTA, the features as GFF3, and the protein sequences of its genes as FASTA. |
| `order` | Turn the primers in `demo/primers_to_order.csv` into an IDT order file (25 nmol scale, standard desalting). |
| `ncbi` | Find the NCBI record with accession L09137 (the pUC19 cloning vector) and download it as a GenBank file. |

The `ncbi` request uses the network and sends only the accession. The others work offline. `checkup` is still a beta feature and has no request here.
