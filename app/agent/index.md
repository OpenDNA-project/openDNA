# openDNA for agents: start here

openDNA has a command line with no window. Every recipe in this folder is for one function; read `skills show <name>` before using it.

## How to run it

- Installed: `oDNA-cli.exe <command>` (in the install folder, next to `openDNA.exe`).
- From the repository (run from its root): `.venv\Scripts\python.exe -m app.cli <command>`.
- `skills list` names every recipe. `skills show <name>` prints one.

## Rules for every recipe

1. **You turn the user's data into the tool's strict input.** Tables become CSV with the exact column names of the recipe (unknown columns are an error). If the user's input is clear, translate it yourself; if it is not, ask.
2. **Look before you write.** Run `inspect <file.gb>` first: it gives the length, topology, every feature with its position and strand, and the protein of each CDS. Use `feature` and `region` for exact bases.
3. **`check`, then `run`.** `check` validates the whole input and computes nothing; it reports every problem at once. Fix them, then `run`.
4. **Show the settings before `run`.** `check` prints the effective settings (value and where it came from: your saved openDNA settings, a default, or an override). Show the ones that matter to the user and wait for their OK.
5. **Relay every warning to the user, and every limit reached is one.** A list cut short (`--top`, the 500 repeats, the pieces of `align`, the matches of `ncbi search`, the best 5 PCR pairs), a search that stopped, a sequence read only in part: each says so in `warnings`, never only in a side field. A cut answer is not the whole answer: tell the user what is missing and how to get more.
   Relay every warning (each row's `warnings`, and the top-level `warnings` about settings). The tool never stops on an ambiguity: it picks a default and adds a warning to that row (for example "3 codons give Ala: used GCG, alternatives GCC, GCA, GCT"). The user must see them.
6. **Errors are loud.** A message starting `ERROR:` (exit code 1) is for you to fix or to pass to the user in plain words. Never retry blindly.
7. **A bad row never stops the others.** A malformed row (wrong number of cells, an empty cell) or one that fails is reported as an error record, loudly (in the JSON and on stderr), and every other row is processed. A wrong header (the columns) stops the run, since the table cannot be read. Tell the user which rows failed and why, and offer to fix and re-run only those.
8. **Names are yours to choose, files are never lost.** The `id` of a row names its output file. Choose a naming convention that suits the user (their own names, or a clear pattern such as `<gene>_<change>`). If two ids would give the same file (also differing only by letter case, which a Windows or macOS disk cannot tell apart), or an id has characters unsafe in a file name, openDNA renames the file (`a_2.gb`, `bad_id_.gb`) and says so in that row's `warnings`; each record's `file` is the real name, and the run checks that every file it names is there. Relay such renames.
9. **Limits.** Each adjustable setting is shown in `settings` with its `min` and `max`, the range openDNA's windows accept. Tell the user these values. A value outside them is used as given, with a `warning` ("the result is your responsibility"): only pass one if the user chose it, and repeat the warning to them.
10. **Inputs are never overwritten.** Outputs go to a new folder (or a new file: a command that writes one refuses a name that exists).
11. **Two kinds of command.** The batches (`kld`, `pcr`, `recombineering`, `codon`, `gibson`) take a table and follow rules 3 to 8. The analyses (`sites`, `guides`, `repeats`, `search`, `align`, `msa`, `checkup`, `annotate`, `variants`, `synthesis`) are read-only: one command, one JSON answer with `settings` and `warnings` (rules 4, 5, 6 and 9 still apply), nothing is written unless you pass an output option (`--out`, `--add`, `--apply`). `export` and `order` write one file. `ncbi` is the only analysis that uses the network: tell the user which words it sends before the first call. `wsltools install` also downloads (the optional pLannotate and Terminal tools, several hundred MB): only when the user asks.
12. **Settings.** Thresholds the user saved in openDNA's Settings are used and shown with their source. A value you pass on the command line is for that run only and never changes their saved settings.

13. **Files of several records.** A multi-FASTA or a GenBank file of several LOCUS records is read one record at a time: write `FILE#NAME` wherever a command takes a file (`inspect "genome.gb#NC_002506"`); a bare name is refused with the list of names when the file holds more than one. `msa` takes such a file whole. `split FILE` writes one file per record, `merge A B ... --out NEW` joins files (or `FILE#NAME` records) into one multi-record file (see `skills show split`). Never cut a file by hand with sed.
14. **Do not invent what the user did not say.** When a request leaves out something the command needs (which feature when there are several CDS, which mutation, which region, which cassette, a scale or purification) ask, or, if you must proceed to show something, say plainly what you chose and why and offer to change it. The settings the user did not choose (for example the PCR window) stay at their defaults: when a result is poor, tell the user and propose a change; do not make it yourself.

15. **Circular plasmids.** The tools read a circular plasmid as a circle: sites, guides, repeats, searches (also for a protein), alignments and parts are found through the origin. A region through it is written START > END (`region f.gb 880 20`, a `pcr` target `880..39`, a `--region`, a recombineering or gibson region); a feature that crosses it is used as it is. `inspect` marks it `crosses_origin: true` with `end` below `start`. `codon` also recodes a CDS that crosses the origin.

## Positions and mutations

Positions are 1-based and inclusive (GenBank and HGVS style). Mutations are written in HGVS, with the feature name where HGVS writes a reference ID:

| Type | In a feature (`c.`) | On the plasmid (`o.`, `g.` also accepted) | Protein (`p.`) |
|---|---|---|---|
| Substitution | `lacZ:c.360G>C` | `o.1234G>C` | `lacZ:p.Asp120Ala` or `lacZ:D120A` |
| Deletion | `lacZ:c.360_365del` | `o.1234_1240del` | `lacZ:p.Lys12_Gly14del` |
| Insertion | `lacZ:c.360_361insATG` | `o.1234_1235insATG` | `lacZ:p.Lys12_Gly13insAla` |
| Deletion-insertion | `lacZ:c.360_365delinsATG` | `o.1234_1240delinsATG` | `lacZ:p.Lys12_Gly14delinsAla` |
| Duplication | `lacZ:c.360_365dup` | `o.1234_1240dup` | `lacZ:p.Gly13dup` |

`c.1` is the A of the start codon; `c.-N` is upstream of it. A change through the origin of the plasmid is written with `o.` positions, first to last (`o.3779_2del`, `o.3780_1insATG`); the file keeps its origin.
