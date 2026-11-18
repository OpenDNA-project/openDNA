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
5. **Relay every warning to the user.** The tool never stops on an ambiguity: it picks a default and adds a warning to that row (for example "3 codons give Ala: used GCG, alternatives GCC, GCA, GCT"). The user must see them.
6. **Errors are loud.** A message starting `ERROR:` (exit code 1) is for you to fix or to pass to the user in plain words. Never retry blindly.
7. **Inputs are never overwritten.** Outputs go to a new folder.

## Positions and mutations

Positions are 1-based and inclusive (GenBank and HGVS style). Mutations are written in HGVS, with the feature name where HGVS writes a reference ID:

| Type | In a feature (`c.`) | On the plasmid (`o.`, `g.` also accepted) | Protein (`p.`) |
|---|---|---|---|
| Substitution | `lacZ:c.360G>C` | `o.1234G>C` | `lacZ:p.Asp120Ala` or `lacZ:D120A` |
| Deletion | `lacZ:c.360_365del` | `o.1234_1240del` | `lacZ:p.Lys12_Gly14del` |
| Insertion | `lacZ:c.360_361insATG` | `o.1234_1235insATG` | `lacZ:p.Lys12_Gly13insAla` |
| Deletion-insertion | `lacZ:c.360_365delinsATG` | `o.1234_1240delinsATG` | `lacZ:p.Lys12_Gly14delinsAla` |
| Duplication | `lacZ:c.360_365dup` | `o.1234_1240dup` | `lacZ:p.Gly13dup` |

`c.1` is the A of the start codon; `c.-N` is upstream of it. A change across the origin of the plasmid is not accepted: it is an error.
