# To do

What is planned, in the order we intend to do it. Tick an item when it is done and tested (add the test name in
[TESTING.md](TESTING.md)). Sizes: S = under an hour, M = a few hours, L = a day or more.

## Next

1. ~~**Temporary primers**~~ — done for the Assembly tab: designed primers are drawn in teal ("temporary") with their details on hover;
   *Save primers to database* makes them regular. Primer Design too.
2. ~~**Assembly preview**~~ — done: a *Product* box ends the Assembly map; clicking it shows the exact sequence in the bottom view,
   fragments coloured, homologies / overhangs / toeholds amber.
3. ~~**Digest: ladders and lanes**~~ — done: five ladders chosen by clicking the ladder lane (name in brackets), and a tab per digest
   (+ lane / ✕) with its own enzymes, table and gel lane. The bands are typical sets, not one vendor's exact list.

## User interface

- [x] Assembly is in the view drop-down again (under Tools): choosing it opens the Assembly tab.
- [x] *Detect Features* is in the Annotations group of the view picker.
- [x] The *Display Filter* view is gone from the picker (the Filters button and View menu still open the same panel).
- [x] The pointer becomes a hand over anything clickable in the sequence, map and alignment views and the Assembly map (a feature, a primer, a read, a box, a card).
- [x] History records adding, removing and editing annotations (they already were steps of the document) and now also adding, editing,
  renaming, deleting and importing primers; those can be undone and redone with Edit → Undo / Redo, in order with the base edits. They are
  steps of the document the change was made from (the database is shared by every open file), they never make the file "modified", and
  an undo that finds the primer changed or its name taken since leaves it alone and says so in the status bar. *Not yet*: the Assembly
  tab's *Save primers to database* is not a step (it belongs to no one document).
- [x] Primers get a menu like features: right-click a primer in the sequence view, the map (its thin mark), the alignment view or the
  primer panel for Edit… and Delete from Primer Database (greyed for a temporary primer that is not in the database yet). A click on a
  primer mark in the map selects it.
- [x] Primer Design uses the temporary primers, and opens by default as a half tab with a sequence view beside it showing where
  each primer binds.

## Search

- [x] **One upgraded search.** The main search bar searches DNA and feature names today. Next to it, checkboxes (each narrows or widens
  what is looked for):
  - *Both strands* (the reverse complement as well; on by default for DNA),
  - *Translations*: a protein motif in the translations of the CDSs and of the three frames,
  - *Fuzzy*: up to N mismatches (a number box next to it),
  - *IUPAC codes* in the query (R, Y, N, ...: on by default),
  - *Feature details*: also the qualifiers of features (gene, product, note), not only their names,
  - *Restriction sites*: an enzyme name finds its sites; *Primers*: a primer name finds where it binds,
  - *Selection only* / *All sequences in this file* (what the *Search Sequences* tool does today),
  - a hit counter with next / previous, and every hit highlighted in the views.
- [x] Then *Search Sequences* is removed from the view picker (the new search replaces it).
- [x] *Quick Align* is now *Local Align* (not "BLAST": it is a seed-and-extend search with a banded Smith-Waterman alignment, and BLAST is a
  registered trademark of the National Library of Medicine).

## Protein

- [x] **Placeholder protein viewer (S).** A view that shows the amino acid sequence of whatever is selected (a CDS, or the
  translation of a selection); to be upgraded later.
- [x] **Protein pane (M).** Translated CDS with numbering, features and basic properties (length, mass, pI): the *Protein* view (both
  items are this one view).
- [x] **3D structure in a new window (L).** A provider layer (`structure.get(sequence, source, accession, command)`): look-up in RCSB PDB and
  AlphaFold DB, prediction with ESMFold (the ESM Atlas API) or a user-supplied local command; a source selector (Auto, PDB, AlphaFold DB,
  ESMFold, Local command); cache by sequence hash; pLDDT colouring; a notice that the sequence leaves the machine; a citation for each
  provider. The viewer is a self-contained 3Dmol.js page opened in the browser (BSD-3). *Not done*: the "ColabFold server" provider (the
  public ColabFold service only builds alignments; the structure is predicted by software on your computer, which the local command
  covers); a native viewer.

## Data

- [x] **Import from NCBI (L).** File > Import from NCBI: nucleotide, genome and protein records by accession or search, through the
  NCBI E-utilities, with an optional API key, a cache and a rate limit. A protein record opens in the 3D structure window (and
  can be saved as FASTA). *Not done*: a protein *document* type (editing / viewing a protein file in a window of its own); the other
  databases (ENA, UniProt, Addgene); NCBI Datasets for genomes over 100 million bases.
- [ ] **NEB enzyme catalogue (blocked: permission).** The restriction enzymes New England Biolabs sells, in the Digest tool. REBASE (the
  source) marks its files "All rights reserved" and publishes no licence, so we ask Dr Roberts / REBASE for written permission before
  bundling anything from it. The Digest code already reads degenerate sites and
  two-sided cutters, so a catalogue file is all that is missing. An alternative with no permission needed: a list written from NEB's
  own published catalogue, if NEB allows it.

## Quality

- [ ] The remaining items of [TESTING.md](TESTING.md) (2: a comparison with the primer list stored in a .dna file with primers, and locations with fuzzy range ends).
- [x] Continuous integration: `.github/workflows/ci.yml` (Windows; Rust tests, pyflakes, pytest), first run green.
- [x] Cas12a in the Cas9 tool (TTTV PAM on the 5' side, 23 nt, staggered cut; the Rust search takes a `pam_first` flag).

## Known problems (found by the tests; the fixes are tested in `tests/test_bugfixes.py`)

- [x] A location with a fuzzy end range such as `10..(20.30)` was cut to its first base, silently. It is now kept word for word (the
  feature has no extent, so it is not drawn; typing one in is refused with a message).
- [x] The GenBank writer put the closing quote of a qualifier alone on the last line (Biopython then read a trailing space).
- [x] `Track.query(limit=N)` dropped the first items of the window instead of the last; a GFF line with an empty `Name=` gave an empty name;
  `enzymes.find_cuts` clamped a bottom-strand cut that falls past the end of a linear sequence; a FASTA with a UTF-8 BOM indexed as zero
  sequences; saving one record of a CRLF FASTA wrote it with LF endings; a failed rename left a `*.tmp` file.
- [ ] Still open: `(10.20)..30` and a lone remote location (`J00194.1:100..202`) are kept in the file but have no extent, so they are not
  drawn; the Features table and the sequence view skip them.
