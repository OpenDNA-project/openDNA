# openDNA user guide

Reference for every feature. New here? Start with [TUTORIAL.md](TUTORIAL.md). Developers: [../README.md](../README.md#build-from-source), [TESTING.md](TESTING.md), [RUST_PRIMER.md](RUST_PRIMER.md).

Python API (positions 0-based, end-exclusive):
```python
from opendna_core import Document
doc = Document("demo/pFG051_Chlor.dna")   # or .gb; the format is detected from the content
doc.features_in(100, 500); doc.primers()
doc.insert(120, "GATTACA"); doc.delete(10, 20); doc.undo(); doc.redo()
i = doc.set_feature(None, "misc_feature", "complement(10..30)", [("label", "site")])  # None = new
doc.delete_feature(i); doc.reverse_complement()
doc.binding_sites("gaattcAGCGGATAACAATTTCACACAGGA", 1)  # [(start, end, strand, mismatches)]
from opendna_core import tm; tm("AGCGGATAACAATTTCACACAGGA", na_mm=50, mg_mm=1.5)
doc.save_genbank("pFG051_edited.gb")      # .dna files are read-only
```

The primer database is `primers.tsv` in the app data folder (`%APPDATA%\openDNA\` on
Windows; set `OPENDNA_DATA_DIR` to use another folder). Export/import it from Primers → Primer Database.

FASTA references (`.fa`, `.fasta`, `.fna`, multi-sequence) open like any other file: pick the sequence in
the toolbar; Save writes the edited sequence back into the same FASTA (others are copied through
unchanged). Only one sequence is open at a time: tracks (GFF3/BED/VCF/BAM) follow you to another sequence of the file (each on the
contig of the same name; "chr1" matches "1"), Save as GenBank / a project writes only the open sequence (the status bar says how many
others stay behind), and mapping reads sends the mapper only the open sequence. A GenBank file of several records (a genome and its plasmids) opens the
same way: the map lists the records, one is open at a time, and Save puts it back among the others; as a *background genome* it is
read in full.

**Import → DNA → Sequence File as New Chromosome** adds the sequences of another file (GenBank, `.dna`, FASTA, project; a
multi-record file lists its records, tick the ones you want and rename them if you like) to the open file as new records,
after the ones it has: a one-record GenBank file becomes a multi-record GenBank file, a FASTA file gets more sequences (bases
only: FASTA keeps no features). The new records are written into the file at once (the records already in it are copied byte for
byte, and Undo does not remove the new ones); the sequence that was open keeps its unsaved edits until you save it. An open `.dna`
file, project or unsaved sequence is first saved as GenBank, and the records are added to that file. The alignment view (View → Alignment View, Ctrl+Shift+G; off until you ask for it, on any sequence) shows
tracks loaded with File → Load Track: GFF3/GTF/BED annotations, VCF variants and BAM alignments
(coverage + reads; the `.bai` is made if missing; several files at once). Click a track's name to hide
it, right-click it for expanded/condensed/hidden and a log scale, drag the line under it to resize; the
Legend button in the genome view's header opens a window explaining every colour and mark. Edit → Circular
Sequence (or the Circular toolbar button) switches topology, asking where to break when making it
linear (that base becomes base 1); circular sequences loop in the sequence view (shown twice, a red bar at the
origin), so selections, copy/paste and deletions can run through the origin. Help → Troubleshooting → Performance Timings
shows how long each drawing/editing step takes (`app/perf.py`).

The Display side panel (the funnel in the toolbar, View → Display Filter, Ctrl+Shift+F) switches feature types, primers,
cut sites, CDS translations (amino acids inside CDS arrows; also per CDS from its right-click menu), the 3-frame
translation, the complementary strand (off, the sequence view shows the top strand alone: one line of bases; also View → Complementary Strand) and transcription units on and off: one filter per document, for every view of it (what a tool such as the
Plasmid Checkup draws in them too). **Primers: in this file or only in the database.** A primer that is part of the file (stored in it, made from its sequence, or imported to it) is drawn in the regular
purple. A primer that only the primer database knows, and a designed one not saved yet, is drawn in the *temporary* teal, in the sequence view, the map and
the Alignment view; hover one to read "from the primer database (not in this file)" or "primer in this file". A primer in both is the file's: regular. In the Display
panel, Primers has two switches, *From this sequence* and *From the database*, each followed by a sample of the line it switches. Right-click a teal primer (or a
row of the Primers panel; select several rows to import them together) and choose *Import to This Sequence*: it becomes a regular primer annotation of the
file (one step in the History, the file is marked modified) and Save writes it as a primer_bind feature. Save writes only the file's own primers, so the
database's hundreds of primers no longer end up in every file. Edit…, Delete from Primer Database and Hide in This File stay as they were.
Right-click a feature in any view to edit or delete it. In the feature table, select several rows (click, then Shift- or
Ctrl-click) and right-click: Delete, Hide All / Show All and Show Translation apply to all of them at once (deleting is one
undo step), next to Group as Gene. The toolbar's range box selects a stretch: type 101..200 (or 101-200, 101 to 200)
and press Enter; on a circular sequence 3700..50 runs through the origin, and one number puts the cursor after that base. Ctrl+wheel (or Ctrl+= / Ctrl+- / Ctrl+0) sizes the sequence view's letters.
History (the hourglass in the toolbar, Ctrl+H; shares the side panel with the filter) lists the changes in words, typing grouped,
and goes back to before any of them (or forward again; or drag ▶ Current sequence). Adding, editing and deleting features, and
adding, editing, renaming, deleting and importing primers, are steps too: Edit → Undo brings a deleted primer back. A primer change is a
step of the file you made it from (the primer database is shared by every open file) and does not mark the file as modified.
**Search** (the bar in the toolbar, Ctrl+F; Enter or F3 = next hit, Shift+Enter or Shift+F3 = previous, Esc clears): type DNA, a protein
motif or a name. Every hit is highlighted in the sequence, map and alignment views, and the counter says "3 of 17". It searches by itself
as you type on sequences up to 2 Mb; on a bigger one, or across a file, press Enter. *Options ▾* chooses what is looked for (remembered):
*Both strands* (the reverse complement too), *IUPAC codes* (R, Y, N... in the query; a query shorter than 5 letters of codes such as "cat"
is a name), *Fuzzy: up to N mismatches* (never more than a third of the query), *Translations* (a protein motif, X = any, in the six
frames), *Feature details* (every qualifier of the features, not only the names), *Restriction sites* (an enzyme name finds its sites),
*Primers* (a primer name finds where it binds), and *Where*: this sequence, the selection only, or all sequences of the FASTA file (DNA only;
it opens the sequence of the hit). It stops at 20,000 hits.
**Import and Export menus.** *Import* is grouped by what a file holds: DNA (GenBank, .dna, FASTA or an Open DNA Project; a
record from NCBI; a background genome from NCBI), Features and Tracks (GFF3/GTF/BED, VCF, *Import Features to Features Database*), Proteins (a protein FASTA or an
NCBI protein, shown in the 3D structure window), Alignments (BAM, or map FASTQ reads) and Primers (a primer list). A track, VCF or BAM file is checked against the open sequence: another sequence name, items past its end, VCF REF bases that are not the bases here or a BAM reference of another length give a warning naming what differs, and the file is loaded anyway.

**Primer lists** (Import → Primers): either FASTA (`>name` then the sequence, on one or more lines) or a table with one primer per line, the name then the sequence, separated by a tab, a comma or a semicolon (more columns are ignored). A first line naming the columns (`name`, `sequence`...), blank lines and lines starting with `#` are skipped. A primer is 8 to 200 letters of A, C, G, T or IUPAC codes (either case) with a one-line name of at most 60 characters. Every other line is refused, never guessed: the report names each refused line with its number and the reason, and nothing invalid reaches the primer database. *Export* is grouped the
same way, each group greyed while there is nothing of it: Whole File (a GenBank copy; an **Open DNA Project** `.odp`: the
GenBank file plus the loaded tracks and BAM files (as paths, absolute and relative to the project, since BAM files can be huge), the DNA lock
and the display options; the **NCBI submission** pair `.fsa` + `.tbl` for table2asn), DNA (FASTA, FASTA of the selection, plain bases),
Features (GFF3, BED, NCBI feature table), Protein (the translations of the CDS features, or of the selection), Alignments (each loaded BAM as
SAM) and Primers (the primer database). `.dna` is read but not written (a closed format). File → Save As writes GenBank or a project.
Export → **Image of the Current Tab** draws every pane of the tab you are looking at, laid out as on screen, into an SVG (vector: lines
and text stay sharp at any size, for figures), PNG or TIFF (LZW-compressed). Choose the resolution (DPI, written in the file) and the
size (100 % = one output pixel per screen pixel; changing the DPI keeps the printed size of what you see, so 300 DPI gives 313 %); the
dialog shows the size in pixels and in centimetres. The pane headers (view pickers, buttons) are left out unless you tick them. The
image uses the colours on screen, so switch to the light theme first for a white background.

**Edited bases**: bases you typed, pasted, replaced or inserted since the file was opened or last saved are drawn as blue letters with a thin blue underline in the sequence view (in every tab that shows it; there is no background, so they are never mistaken for the yellow search hits, and a hit over an edited base shows both). Deleting leaves nothing to mark. The marks follow later edits, and Undo / Redo keep them right (undoing an insertion removes its mark). Save and Save As clear them; the automatic backup does not. The blue follows the light or dark theme.

**DNA lock** (🔓 DNA / 🔒 DNA button at the start of the toolbar, Edit → Lock DNA, Ctrl+Shift+L): while it is pressed the bases cannot change (typing, paste,
delete, reverse-complement, making a plasmid linear, applying variants, and Undo/Redo of a base change are refused, with a message in the
status bar). Features, primers, the display and everything else can still be edited.

**Background genome from NCBI** (Import → DNA → Background Genome from NCBI…): type an organism (Escherichia coli, or its taxonomy number) or
an assembly accession (GCF_000005845.2); *NCBI reference genomes only* keeps NCBI's reference assembly of each species. *Download and open* fetches
the sequences (FASTA) and annotation (GFF3) from NCBI Datasets into the `genomes` folder of the app data (once), and opens the FASTA with the
annotation as a track. *Use for off-target checks* makes the downloaded FASTA itself the reference of the Primers and Spacer Finder checks, with every one of its sequences (all chromosomes).

**Import from NCBI** (Import → DNA → From NCBI…): choose *Nucleotide / genome* or *Protein*, type an accession (NC_000913.3, P69905) or
words (pUC19 cloning vector), and pick a record from the list. *Open* shows a nucleotide record with all its annotation in a new window (an
empty window takes it; nothing that is open is touched; Save asks where to write it), and a protein record in the 3D structure window; *Save as
file…* keeps the GenBank or FASTA file. Records are kept in the `ncbi` folder of the app data, so asking again does not download again. Only
your search words (or accession) and, if you type them, an API key (10 requests a second instead of 3; free in an NCBI account, kept in plain
text in the app data) go to NCBI. Records over 20 million bases ask first;
over 100 million are not fetched (use NCBI Datasets).
**Reporting a problem or asking for a feature.** The Help menu has four items. *Report a Problem…* opens a GitHub form already filled with your version, system and the last lines of the log, and *Request a Feature…* opens the idea form. *Open Log Folder* shows where `opendna.log` is kept, and *Copy Diagnostics* copies the versions and the last log lines so you can paste them into an email. If the program meets an unexpected error it does not close: it records it in the log and shows one small window, "Something went wrong. Report it?", with *Report on GitHub*, *Copy details* and *Close*. Nothing is sent by the program itself, and before anything goes into a report or the log, sequences, file names, folder paths and your user name are replaced by placeholders. You always see the report on GitHub before you send it. Each issue carries a status label (new, under review, planned, in progress, done, released, declined), so you can follow what happens to it.
**Opening files by double-click (File Associations).** The first time the program starts, it asks whether you want to open .gb, .dna and FASTA files with openDNA by default. Answer Yes and openDNA is added to Windows' "Open with" list for these file types. This is for your user account only and needs no administrator rights and no installer. Windows does not allow a program to make itself the default, so you finish in Windows Settings (Default apps): the program opens that page; click openDNA, then choose it for each file type you want (or right-click a file, then Open with > Choose another app > openDNA > Always). You can register or remove the registration at any time in Help > File Associations…. If you move the program file, the registration still points at the old place; the program notices at the next start and offers to update it. Because the program file is not code-signed, Windows SmartScreen may warn the first time you run a downloaded copy; right-click the .exe, choose Properties and tick Unblock. Each double-click starts a new window with its own copy of the program, and it takes a few seconds to open because the single-file program unpacks itself first.
**Updates.** openDNA can look for a newer version on GitHub once a day when it starts (Tools > Settings > Updates turns this off), and Help > Check for Updates… looks right away. When a new version exists, a window shows what changed. "Install and restart" downloads the new `oDNA_beta.exe`, checks it against the published SHA-256 checksum (a download that does not match is deleted, never installed), replaces the running program and starts the new one. "Later" asks again next time, and "Skip this version" stays quiet until the next release. Only github.com is contacted, over HTTPS, and nothing but the version request is sent. The checksum protects against a damaged or tampered download, not against a hijacked release, so download from the official releases page only. If openDNA is in a protected folder (such as Program Files) or an antivirus holds the file, it says so and links to the releases page so you can replace the file by hand. Run from source, it cannot replace itself. Until the project is public, updates are switched off and the menu item says so.
**Start page** (Help > Welcome): opens at the first start of each day. On the left: what changed (once, after an update), the recent files (also in File → Open Recent; click one to open it; right-click to pin it to the top or remove it from the list, the file itself is not touched), buttons for Open File, New Sequence, Import from NCBI, Assembly and the demo plasmid (opened as a new unsaved file), and Check for Updates, whose answer shows next to the button. On the right: the sponsors, a tip of the day and links. Tick *Don't show at startup* to stop it opening by itself; untick it from Help > Welcome to bring it back. To open it at startup anyway, start the program with `--welcome` (`oDNA_beta.exe --welcome`, or `--welcome` after a file name); it ignores "already shown today" and the tick, and leaves the tick as it was.
**Primer design** (select the target, then right-click it → *Design Primers for This Sequence…*, or Primers → Design Primers for
Selection…): opens the designer in its own window over the main one, already run on the selection; the selected pair is drawn in the sequence view as **temporary primers** (teal, "pair1_F" /
"pair1_R": they are not in the primer database and have no menu yet) so you see where each one binds, and each row you click shows its pair.
*Save pair to database…* names them and makes them regular primers; closing the tab forgets the unsaved ones.
**Protein view** (pick *Protein* in any pane's view picker): click a CDS, or select a stretch of DNA, and the view shows its amino acids numbered in
blocks of ten, with the length, mass (kDa), isoelectric point and charge at pH 7 (up to the first stop codon; estimates for the free chain),
and the features inside the CDS (click one to select it). For a stretch of DNA, choose the reading frame (+1 to −3) beside the title. Select
letters of the protein to select their codons in the sequence; the residues of whatever you select in the sequence are highlighted. While a
CDS is shown, selecting inside it keeps it shown.
**3D structure** (the *3D structure…* button of the Protein view): a window to choose the source: *Auto* (the cache, then AlphaFold DB when the
CDS has a UniProt accession in its db_xref or you type one, then the PDB, then ESMFold), *PDB* (RCSB's sequence search, 90 % identity or
more, 25 amino acids or more), *AlphaFold DB* (by UniProt accession; only the accession is sent), *ESMFold* (a prediction from the sequence
alone, up to 400 amino acids, on Meta's ESM Atlas server) or *Local command* (your own AlphaFold / ColabFold / Boltz: a command with
{fasta} and {pdb} in it). The window says which service the sequence goes to before you press *Get structure*; nothing is sent until then,
and the sequence only leaves this computer for the PDB, ESMFold and Auto sources. Structures are kept in the `structures` folder of the
app data (the same sequence is not asked for twice; delete the folder to forget them). *Open 3D view in the browser* shows a self-contained
page (3Dmol.js): cartoon, sticks or surface, coloured by confidence (pLDDT, as AlphaFold DB does), rainbow or chain. *Save PDB…* keeps the
file. A prediction is a model, not a measurement: low-confidence (orange) parts are not reliable.
**Right-click a primer** (sequence view, map, alignment view or Primers panel) for *Edit…* and *Delete from Primer Database*. *Hide in This File* stops drawing that primer on this file only (it stays in the primer database, and the
choice is remembered for the file and kept in projects); the Primers view still lists it, greyed and marked "(hidden)": right-click it there →
*Show in This File*. Select several primers in the Primers panel and right-click: delete them from the database, import them, hide or
show them all at once.

**Feature direction.** A feature is drawn with an arrow (its strand) unless it has no direction: *misc_feature* has none by default, and
the feature dialog's *Direction* (Forward → / Reverse ← / None) changes it for any feature. GenBank has no "no direction" location, so a
choice that differs from the type's default is saved as a `/directional="yes"` or `"no"` qualifier; `.dna` files keep their own
direction setting. The feature table shows → / ← / —. The selection has
draggable blue [ ] ends. Tools → Detect ORFs (also in the *Annotations* group of a pane's view list; on a circular sequence also those that cross the origin) links ORFs already annotated as CDSs and searches protein
motifs in all translations. Stop codons show as red blocks in translations.

**Genetic code.** Tools → Settings → *Genetic code* picks the NCBI translation table (1 standard, 2 vertebrate mitochondrial, 4 mould /
protozoan / mycoplasma, 11 bacterial / archaeal / plant plastid, and the other tables of the NCBI list except 27, 28 and 31, where a
codon is a stop or an amino acid depending on its place in the gene). The default is table 1. The line under the choice lists the
stop and start codons of the selected table. The choice drives the 3-frame translation, the Protein pane for a selected stretch,
Detect ORFs (stops), protein-motif search in Find and in Detect ORFs, Detect Common Features at the protein level, the Plasmid Checkup's
stop-codon check and the FASTA export of a translated selection. When it is not table 1 the headers say so ("Table 11": the Protein pane, and
the Display panel's 3-frame translation switch). *Allow alternative start codons* (off by default) lets Detect ORFs start at the table's
other start codons (GTG, TTG, CTG, ATT... in table 11) as well as ATG. A CDS is different: one with `/transl_table=N` is always read with
table N, whatever the setting; `/codon_start` (1-3) is honoured in the sequence view, the Protein pane and the exports; a first codon that
is a start codon of the table is shown as Met (as GenBank writes `/translation`) unless the 5' end is marked partial; and a simple
`/transl_except` (`(pos:100..102,aa:Sec)`, also `complement(...)`) is applied. A `/transl_except` that cannot be followed (for example a stop
completed by poly-A) leaves the CDS drawn plainly and its tooltip says "translation not reproducible".

**Keeping `/translation` right.** Feature → *Keep /translation up to Date* (on by default, remembered between sessions): whenever a
GenBank file is saved (File → Save, Save As), the `/translation` of every CDS that already has one is recomputed from its own bases, its
`/codon_start` and its `/transl_table` (else the Genetic code setting), without the final `*`, like GenBank writes it. A CDS with no
`/translation` never gets one. Left alone: `/pseudo` or `/pseudogene`, a partial 5' end with no `/codon_start` (the frame is undefined),
a `/transl_except` that cannot be reproduced, a stop codon inside the reading frame, or a CDS whose length is not a multiple of 3 (unless
its 3' end is partial). The status bar reports it ("Saved ... · Updated /translation of 3 CDSs, skipped 1 pseudo"). The refresh is one
undo step, added just before the file is written (Ctrl+Z after saving undoes it, and the file then counts as modified); when nothing
differs no step is added, so saving twice changes nothing. Autosave never does this: it only writes a recovery copy. Feature → *Update
/translation Now* does the same on demand for the selected CDS, or for every CDS when none is selected (one undo step).
Help → Keyboard Shortcuts (F1) lists every shortcut (menus read live; per-panel ones in `app/shortcuts.py`).

Tools → Local Align (Ctrl+Shift+B) aligns a pasted sequence (FASTA or plain bases; the clipboard is
offered) to the open one, both strands, with mismatches and gaps (seed-and-extend with a banded
Smith-Waterman, `align` in `core/src/align.rs`; ~0.1 s on a 4.6 Mb genome). Every matching piece of the query is
reported on its own (up to 30), so a query that is a patchwork (a structural variant) reads like a story: one row per piece, **in query
order**, with its query range, its template range, the strand in words (*forward*, or *↺ reverse complement (inverted)*, shaded purple),
the identity, and a *Before this piece* column that says how it follows the previous piece: *contiguous*, *insertion in the query*
(N bp unmatched in the query), *deletion from the query* (N bp of the template skipped), *inversion breakpoint (strand changes)* or
*template jumps back* (a duplication, or pieces in another order). A line above the table sums it up, for example
"4 pieces: 3 forward, 1 inverted (reverse complement); gaps in the query between pieces: 200 bp unmatched (#1→#2); inversion breakpoint (#2→#3)".
A piece is never joined to the next across a long stretch that matches nothing (about 40 unrelated bases or more); SNPs and small
indels (up to about 20 bases) stay inside one piece. Click a row to select its template span; the hits are drawn in the alignment view as a
"Query hits" track with **one row per alignment**, numbered like the table's rows (overlapping or repeated hits never sit on top of each other; inverted ones
hatched and purple, with a ↺ in the label). Each row has its label (#, → or ↺, identity) in the names column at the left, which stays put when you pan;
click a label or a bar to pick that alignment: its bases are selected, its row is tinted, and the same row is selected in the table (and the other way round).
Clear hits removes the track. Loaded BAM/SAM files keep their usual coverage, variants and packed reads tracks (thousands of reads are not given a row each).
The query box takes plain bases or FASTA (the header line, spaces, line breaks, digits and `-` are ignored, lowercase and RNA's U are
accepted, only the first FASTA record is used, up to 200,000 bp); anything else, such as a protein, gets a message instead of a result.
*Paste from clipboard* replaces the query with the clipboard's DNA (a message says when the clipboard is empty or has no DNA), and
Tools → Align Clipboard… opens Local Align with the clipboard already aligned.

Tools → Assembly opens the Assembly tab (a tool tab next to the others, not a pane: not in the pane drop-down, never split; closing it keeps
its content). Choose the method (Gibson/HiFi, Golden Gate with BsaI, BsmBI, BbsI or SapI, or Sidewinder) and the number of fragments
(1-12); each is a gray box. Right-click a fragment for everything that can be done to it (source, move left or right, reverse complement, destination vector, clear); the actions are a column on the right of the map (the verdict, Design primers, Save fragments, Save primers to database and, last and highlighted, Perform assembly), and the controls for what the fragment is made of are a column on the right of the sequence view. Click a box to take its DNA from a currently open file, a recent file or a file you browse to; the file is
shown in the bottom half with the editor's own sequence view (read-only), where you select a stretch or click a primer and shift-click the one on the
other strand (primers of the file or of the primer database) to use their PCR product, tails included. The junctions are checked as soon as every box is filled, on the real ends of the fragments: Gibson needs the end of one
to be the start of the next (15+ bp), Golden Gate needs a site at each end and the same 4 nt overhang on both sides of a junction, and a
destination vector (right-click it, *Destination vector*) is cut at its two sites. Each junction is a compact pill between the fragment boxes with what joins the two (its length), Tm and GC, in the colour of its verdict; a failed check is named on the pill (a
red error if the homology is not unique in the product), hover it for every reason. Click it for the **junction inspector**: the product around the joint in the program's own sequence view, on the top strand alone, with the joint selected, each fragment a bar in its colour, the joint a slate bar, what a primer's tail added an orange bar, the primers drawn as primers (hover one, a tail or a feature as in any view), a sentence saying how much of the joint each primer adds and how much the two fragments already shared, and the checks of the junction listed under it; *Full report…* lists every junction, check and oligo. A circular assembly has one pill more, the last, and an arrow from it back into Fragment 1; Golden Gate overhangs are checked for repeats, palindromes, reverse complements, ones one
base apart, and against the 25 high-fidelity overhangs of Potapov et al. 2018. *Design primers* designs new primers carrying what the method needs, drawn on each fragment's sequence below (arrow on the binding site, tail dashed);
the product beside the map (a ring for a circle, a bar for a line, each fragment in its colour with a dot on every junction; clickable once every junction fits) shows the exact sequence below like a fragment, each fragment a coloured bar and each homology, overhang or toehold slate, with the primers that made the fragments (typed, picked or designed) drawn on it like any primer (what each primer's tail added is shown in the junction inspector, not here); the product *Perform assembly* opens has no such bars: it keeps the fragments' own annotations and those primers, so a file saved from it stores them;
*Perform assembly* (enabled when every junction checks out) opens that product in a new window; *Save fragments…* writes each PCR product to
a GenBank file; *Save primers to database* keeps the oligos. Sidewinder uses coding/barcode oligos with 10 nt toeholds (after Robinson et al.,
Nature 2026; simplified). Help → Citations lists the references (engine `app/assembly.py`, panel `app/assembly_panel.py`, overhang data
`app/overhangs.py`).

**Spacer Finder** (Tools → Spacer Finder…, or the *Spacer Finder* view in any pane): Cas9 / Cas12a guides.: lists every guide of SpCas9 (NGG, 20 nt), SaCas9
(NNGRRT, 21 nt) or Cas12a (TTTV, 23 nt) in the selection or the whole sequence, on both strands, best first, with where the nuclease cuts (Cas9: 3 bp from the PAM, both strands at once; Cas12a: the PAM is on the 5′ side, shown first, and the two strands are cut in different places, 18 and 23 bases after the PAM, leaving a 5 nt 5′ overhang; both drawn on the sequence), GC
and a checklist score (GC 40-70 %, no TTTT, no long runs: a rule of thumb, not a trained model). *Check specificity* counts, for each
guide, the other sites it matches up to 1-4 mismatches (as 0/1/2/… mismatches) in the open sequence and in the background genome,
and marks it ✓ / ? / ✗ (another perfect site is a ✗). The background genome is the one of the Primers tool (Primers → Options): both
tools read the same setting; every sequence of a multi-record FASTA or GenBank reference is searched (counts are added up). Selecting a row selects the site in the sequence; *Add selected as feature* annotates it: the guide is drawn as an orange box (the protospacer), a yellow one (the PAM), a blue bar under both with its name (all text black, in both themes) and a hollow red arrow pointing up at the cut, and is saved as an ordinary `misc_feature` (a `/label`, a plain-words `/note` with the guide, PAM and cut, and a `/OpenDNA_cas9` qualifier that only this program reads), so other programs show a normal annotation; hover it to see where it cuts. **Automatic check:** when a background genome is set, the first 100 guides of the list are checked by themselves in the background (a progress bar and *Cancel check* appear; the table stays usable and fills as results come in; columns sort by a click on their header). Guides beyond the first 100 keep *Check specificity*, and the status line says "100 of N checked". With no background genome set nothing is checked automatically; the *Genome* column header says so. **Right-click** a guide (or several selected ones) → *Add as feature*, the same as the button.
The search is
done in Rust (`core/src/cas9.rs`): about 0.6 s for 100 guides against a 5 Mb genome.

**Restriction digest** (Tools → Restriction Digest…): tick enzymes to see the fragments, the virtual gel and the cuts on the sequence. The
**filter bar** on the left narrows the list: a box for a name or a recognition site, then the same criteria as the Display filter (number of
sites here, size of the site, palindromic, Type IIS, the ends left, owned) and its one-click presets; it is the Digest's own (it starts with
every enzyme that cuts) and does not change what the views draw. Sites can be 0: *0 to 0* lists the enzymes that do not cut
this sequence (to clone into a site it lacks), *0 to any* every common enzyme. The enzyme list shows, for each enzyme, the number of **Sites** in this
sequence, the **Recognition site**, the **Ends** it leaves (5′ or 3′ overhang and its length, or blunt) and **Owned**: tick it for the
enzymes your lab has (the enzyme inventory, the same for every sequence; the Owned criterion of the filters reads it). Click a header to
sort; hover an enzyme for where it cuts. Each lane of the gel is named by its enzymes, one per line. **Click a fragment** in the table: it is selected on the sequence, map and genome views (a fragment across the origin of a circle too) and the cuts at its two ends are shown heavier with the enzyme names above them; hover them for the pill. Selecting anything else clears them. A tick stays when the filter hides the enzyme. **+ lane** adds another digest next to the first: each tab has its own
enzymes and fragment table and its own lane on the shared gel (✕ on a tab removes it; the cuts drawn on the sequence are the selected tab's). Click the ladder
lane (its name is in brackets) to choose another: 1 kb, 1 kb Plus, 100 bp, 100 bp Plus or 50 bp (typical band sets of those ladder types, not one vendor's exact list). The finder reads degenerate (IUPAC) sites and enzymes that cut on
both sides of their site, so a longer enzyme list only needs more entries (the list is the 47 common enzymes of `app/enzymes.py`).

**Cuts on the sequence.** The restriction digest (ticked enzymes) and the Cas9 tool (the selected guide) draw their cuts in the sequence view,
on each strand where it is cut: a red stroke through the top strand and one through the bottom strand, each ending in a small cross (not the
plain red bar of a circular sequence's origin), joined by a dotted line when the ends are sticky (EcoRI cuts the top strand after base 1 of its site and the bottom strand after base 5; Cas9 cuts both in the same place,
3 bp from the PAM). Hover a cut (or an enzyme name shown by *View → Restriction Sites*) to see a *pill*: the enzyme's name and the recognition site drawn as double-stranded DNA that has really been cut (both strands, base-pair ticks, the two pieces pulled apart, the strands ending where they are cut, so a 5′ or 3′ overhang or a blunt end is visible; a Type IIS cut far from its site is marked with …), and one line of words: the site with the cut marked (G^AATTC), the overhang and its bases, the base each strand is cut after and, for a Type IIS enzyme, the distance from its
site; a Cas9 / Cas12a cut says the guide, PAM, nuclease, strand and cut position. They go away with the tool. A primer whose 5' tail hangs off the start of a sequence gets room: the bases move right
until the whole primer shows. Clicking a base (or a feature, read or primer) in the alignment view moves the sequence view to it.

**Cog and ? of a pane.** Next to each pane's view picker, the cog (⚙) opens the options of that view: the matching tab of
Tools → Settings (Sequence and Map: Display; Alignment: Alignments; Protein and Detect ORFs: Genetic code; Detect Features; Plasmid Checkup), or the Primers options
(Primers, Spacer Finder: the background genome). A view with no options has a greyed cog whose tooltip says so. The ? opens a
short guide to the view (what it shows and its main clicks and keys); the texts are in `app/view_help.py`. Both buttons can be
reached with Tab and pressed with Space.

**Databases.** The Databases menu opens each database the program keeps as a tab of the window (like the Assembly tab: it works
with no file open, follows the window when another file is opened, and closing it only takes it out, so opening it again finds it as
it was): the **Features Database** and the **Primer Database** (Ctrl+Shift+D, also Primers → Primer Database). Each tab has a header
(its name, how many entries, its buttons, **Help**), the program's own list on the left (the Features view; the Primers view) and, on
the right, the sequence of the entry you click, alone and read-only: a feature on its own bases with the feature over them, a primer
on its own bases. Click a column header to sort (again: the other way); right-click it to filter by values.
The Features Database lists each feature with its **Source** instead of a location: the publication or database it comes from (or
**Local** for your own), then the record or file it was taken from, e.g. "Cloning vectors (NCBI records) · pBR322 (J01749)" (filtering
the column by values goes by the source's name); **Measured** appears only when a source measured a
strength. *Sources* in the header ticks the sources that Detect Common Features and the Plasmid Checkup look for. Right-click a
feature to delete it, which is only possible for Local features: a published source is switched off whole. The Primer Database adds
the **Note** and the day each primer was **Added**; open a primer (its ▸) to see where it binds in the background genome; its header has New, Edit, Delete, Import, Export and Options (double-click a
primer to edit it). **Help** explains the tab and shows where the database is kept, how many entries it holds and what is wrong in it
(rows with no sequence, letters that are not bases, names used twice), with *Open Folder*; Databases → Database Status shows the same.
The restriction enzymes are listed in Tools → Restriction Digest, with the Owned box of the enzyme inventory.

**Binding sites of a primer.** In the Primers view, open a primer (its ▸) to see its sites, grouped by where they are: *On this
sequence*, *On the other sequences of this file* (a file of several records), and *On the background genome* (the toolbar's genome
box), the last searched when you open the primer (it says *searching…* meanwhile; a large genome takes a few seconds). Each site
shows the primer's annealing bases with those that do not pair with the template in **red**, its length and Tm, and how it matches
(in the genome: the whole primer, or only its 3′ end). Click a column header to sort.

**Filter by values.** Right-click a column header to tick which of its values are shown, like a spreadsheet's filter (a search box
narrows long lists; *All*, *None*, *Clear filter*). A filtered column's header ends with ▼. It is on the columns that hold categories:
the Features view (Type, Strand; Source in the Features Database), the Primers view (Match, Structure, Added), Detect
Common Features (Type, Source), Detect ORFs (Strand, Annotation), the Digest's fragments (Ends), the Spacer Finder (Strand, Here,
Genome, Verdict), Local Align (Strand), the Plasmid Checkup (Type, From; Replicates in, Copies) and the NCBI searches (Organism;
Assembly, Level). The filter stays when the table is refilled; a value it has not seen before is shown.

**Background genome.** The genome box of the toolbar (a ring icon and the file name), right of the range box, is the genome the primers and the Spacer Finder's guides
are checked against: it lists the last three used, *No background genome* switches the check off, and *Add genome file…* adds another file (FASTA
or GenBank). It is the same setting as Primers → Options → Genome file, and a genome downloaded with Import → DNA → Background Genome
from NCBI goes in the list too.

**Views, tabs and themes.** Each pane has a view picker: a popup with a column per category (Views, Annotations, Tools), each in
alphabetical order. A tab is named after the views in it (the first three, then a count) unless you rename it (double-click it; an empty
name goes back to automatic). The sun button of the toolbar (or View → Theme) switches between Light (a sun), Dark (a dark sun) and System (a sun with an A: follows the computer), remembered; every view, the map and the alignment view included, redraws in the new colours at once; the colours the picker uses live in
`app/theme.py` (`ROLES`) so other widgets can reuse them. Help → About shows the licence notices.

**Every file opens the same way:** one tab, the sequence across 70 % of the window, the map above the feature table in the 30 %
on the right ("+" adds more tabs). **A crowded map keeps its largest features:** when there are more than it can show legibly at its size (about one
per 20 px of backbone, in lanes at least 5 px thick), the smallest are left out and a note in its corner says how many.
Enlarging the pane brings them back; the sequence view and the feature table always have them all, and a map that is not
crowded draws everything. **A position scale** runs inside the features (under them on a linear map): ticks every 1, 2 or
5 x 10^k bp, the step picked from the map's size so the labels never crowd (every 10 bp on a 100 bp oligo, every 500 bp on a
6 kb plasmid, every 0.5-1 Mb on a 10 Mb genome), with smaller unlabelled ticks between when there is room, in bp, kb or Mb to
suit the step. Rearranging the tabs lasts as long as the window; the next file you open starts from the standard layout again, so what
you see never depends on what you left behind. The one exception is an Open DNA Project (`.odp`), which opens with the tabs
it was saved with. Windows show tabs of split panes: pick any view or tool in a pane's header ("+" after the last tab adds
one; right-click a tab to rename/duplicate it, middle-click duplicates; View → Reset Tabs puts the standard layout back).
A pane builds its view the first time you look at it, and drops it when it shows something else, so tabs you never visit and
views you have switched away from cost nothing. Click a primer, then shift-click one on the other strand to select what they amplify.
Feature → Group as Gene makes a gene own its promoter/RBS/CDS/terminator (saved as a shared
`/locus_tag`). Tools → Detect Common Features finds known features from the **Features Database**
(Databases → Features Database). It is one database made of **sources**, one file each (`app/parts/*.tsv`, named by their
publication or database, plus your own features in `parts.tsv` in the data folder, the source **Local**: right-click a feature → Add
to Features Database, or Import → Features and Tracks → Import Features to Features Database; each keeps the file it came from).
*Sources* (in the tool and in the database's tab) ticks the sources that are looked for; a source switched off is used by nothing, the Plasmid Checkup included. A part may carry the **strength**
its source measured, with the units and the host it was measured in; units differ per source, so a number is only
comparable inside its own source, and a part nothing measured shows none. What ships: cloning-vector features from
NCBI records, and engineering parts (markers, counter-selection, recombinases, CRISPR nucleases, polymerases, reporters,
regulators) whose protein comes from UniProt and whose DNA comes from the GenBank record UniProt cites.
Tools → **Plasmid Checkup** opens a tab of panes, the checkup beside the document's own map and sequence (it is also a pane in
the *Tools* group, to put anywhere: picked in a pane, it adds a map and a sequence pane beside it when the tab has none, which you can close). It traces the transcription units of the construct and says for each
CDS whether it should be expressed: *Expressed*, *Likely* (it relies on a predicted promoter or RBS, or a weak RBS / Kozak),
*Unlikely* (no RBS, the transcription start after the RBS, no polyA signal, no start/stop codon or a stop codon inside) or *Not
expressed* (no promoter upstream on its strand, a terminator before it or inside it, a Pol III promoter). The verdict is given per
host: *Bacteria* always, and *Eukaryotic cells* when the plasmid has a eukaryotic promoter (CMV, EF-1α, CAG, SV40, PGK, GAL1, TEF...).
Each kind of unit ends at its own kind of terminator: a bacterial unit reads through polyA signals and yeast terminators, a
eukaryotic one through bacterial terminators and ends at a polyA signal (or a yeast terminator); a eukaryotic CDS needs a Kozak
context (purine at -3, G at +4) instead of an RBS. The *Replication origins* table names each origin (pUC, pMB1/ColE1, p15A,
pSC101, R6K, RK2, RSF1010, pBBR1, F, f1, SV40, 2µ, ARS/CEN...), where it replicates, its usual copy number and what it needs (R6K: a
pir+ host), and warns when there is no bacterial origin or two of the same incompatibility group. One tree lists the
transcription units, each with the CDSs it carries and their verdict per host (folded units stay folded), then the CDSs no unit
reaches ("Not transcribed"). While the panel is shown, the units and the library / predicted elements it found are drawn in the
document's views as temporary features, under the document's display filter: what was found and how far each unit runs can be
checked by eye. They can be selected and hovered, not edited; they go when the panel is hidden (another tab), and the document
is not changed until you add them. It uses the
annotated promoters, terminators and RBS (and `regulatory` features with that /regulatory_class),
the promoters and terminators of the parts library, and, in sequences up to 300 kb, motifs it finds itself: sigma-70 promoters
(-35 TTGACA, -10 TATAAT), intrinsic terminators (G-C rich hairpin + T tract) and Shine-Dalgarno sites. A unit runs from its promoter to
the first terminator downstream on the same strand (all round a circular plasmid if there is none). The notes name what a promoter
needs (IPTG for lac/tac/trc, arabinose for araBAD, aTc for tet, T7 RNA polymerase for T7...) and whether the plasmid carries the
regulator. Select a row to select it in the document's views; *Add units to the document* adds the units as mRNA
features, and the ticked library or predicted elements can be added too. A promoter with no terminator downstream makes a unit that
runs all round a circular plasmid (to the end of a linear one): every CDS on that strand is reached, and those more than 3 kb away
get a note. **Transcription units are drawn as units**, in the sequence view and on the map: a thin line over the promoter, a
thicker one from the transcription start (+1) on, and a T bar where the terminator stops it (an arrow head when nothing does: it
runs on). A "+1" is written at the start only when its base is known: an annotated TSS (a `TSS` feature, or one named "+1",
"TSS" or "transcription start", inside the promoter or up to 50 bp after it) or a T7 / T3 / SP6 promoter, whose +1 is the base
after its conserved box (…CACTATA|G for T7). Elsewhere the start is estimated (the promoter's 3' end, ~7 bp after a predicted
-10 box): the line still widens there, without a "+1", and the unit's notes say "+1 estimated". Overlapping units each start at their own promoter, in lanes of their own. Their name says how they end, for bacteria
the kind of termination: *intrinsic* (a hairpin + T tract was found in the terminator, or it was predicted from one) or
*rho-dependent?* (an annotated terminator with no such hairpin: rho-dependent, or not a terminator at all); eukaryotic units end
at a polyA signal. Units added to the document are mRNA features with a program-specific `/OpenDNA_tu` qualifier (like
the Cas9 guides' `/OpenDNA_cas9`): other programs see an ordinary mRNA; GFF3 and the NCBI table leave the qualifier out. It is a check of the design, not of expression levels: terminators leak,
and predicted elements are motifs only. Tools → Settings tunes each tool. In the Alignment view, each BAM has a
*variants* track under its coverage: every position where the reads disagree with the sequence (solid: most reads; faded: a minority).
Right-click a variant to *Change the Sequence to Match* or *Annotate This Variant* (a `variation` feature); right-click the track for both on
all consensus changes at once.

Tools → Map Reads to This Sequence runs bwa mem or minimap2 (recipes in `app/tools.toml`) and opens
the resulting BAM. `python tests/fake_reads.py [plasmid] [folder]` simulates ONT and Element reads of
a plasmid with typical cloning errors (mutations, indels, inversion, duplication, mixed populations) and
maps them the same way; `truth.tsv` lists what was planted. On Windows they run in WSL: Tools → Install Mapping Tools creates the conda env
`opendna-tools` (bwa, samtools, minimap2 from bioconda) inside WSL's conda.

Unsaved work is autosaved every minute to `<file>.autosave.gb` next to the file; reopening the
file after a crash offers to recover it. The copy is deleted on save or a clean close.

`core/vendor/gb-io` is a patched copy of the gb-io crate (GenBank writer wrapping fix);
see the note in `core/Cargo.toml`.

Use the dedicated `opendna` env, **not** conda `base`: base ships ICU 73 in `Library\bin`, which
conda's Python puts on the DLL path; PySide6's `Qt6Core.dll` then loads it instead of Windows' ICU and
fails with "DLL load failed ... The specified procedure could not be found".

**Background genome with several sequences (off-target / mis-priming check).** Primers → Options → *Genome file* may be a GenBank, `.dna` or
FASTA file with one sequence or many (a 24-chromosome genome, a multi-record GenBank). The Primers dialog and the primer designer search
every sequence, one at a time (a FASTA through its `.fai` index in windows of 32 Mb, so a genome never sits in memory whole), on a worker
thread: the dialog shows "searching chr3 (3 of 24), 40 %", the designer a progress dialog with Cancel. Hits name their sequence
(`chr3:1,200–1,221 →`); red = the whole primer matches somewhere, yellow = only its 3′ end, green = nowhere. The open sequence is left out of
the comparison (when it is one record of the reference). A search that passes the time limit (10 minutes) or is cancelled is yellow ("only 5 of
24 sequences were searched: not conclusive"), never green. The reference label reads "24 sequences, 3.1 Gb". A FASTA record is treated as
linear (the format cannot say it is circular); a GenBank record keeps its own topology. The cost is per primer: a human-size genome takes minutes
for each.

**Input checks.** Every field that takes a value was audited (`docs/INPUT_AUDIT.md`): sequence fields take pasted FASTA or
numbered blocks and say why when the text is not DNA, numbers outside the possible range are refused with a message, long text
has a maximum length, and a dialog that finds a mistake stays open with what you typed.
