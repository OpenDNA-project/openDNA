<h1 align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/images/wordmark-dark.svg">
    <img src="docs/images/wordmark.svg" alt="oDNA: open Design, Navigation and Assembly" width="380">
  </picture>
</h1>

<p align="center"><strong>openDNA: a free, open-source, offline DNA sequence editor.</strong></p>

<p align="center">Plasmids and whole genomes in the same window. Edit a sequence and every annotation stays on the right bases.
Design primers, Cas9 guides and assemblies, on your own computer, with nothing sent anywhere.</p>

## Downloads

- **Windows** (10 or 11, 64-bit) — **[download openDNA](https://github.com/OpenDNA-project/openDNA/releases/latest/download/oDNA_beta.exe)** (one file, no installer, no Python needed).
  Also [as a folder](https://github.com/OpenDNA-project/openDNA/releases/latest/download/oDNA_beta-folder.zip) if you want to replace the Qt libraries yourself.
- **macOS** — [not available yet](https://github.com/OpenDNA-project/openDNA/releases); [run from source](CONTRIBUTING.md#run-from-source) meanwhile.
- **Linux** — [not available yet](https://github.com/OpenDNA-project/openDNA/releases); [run from source](CONTRIBUTING.md#run-from-source) meanwhile.

> **Beta (0.5): keep copies of your files.** The exe is unsigned, so Windows SmartScreen warns you: check its SHA-256 against the
> release notes, then *More info* → *Run anyway* ([SECURITY.md](SECURITY.md)). New here? The [tutorial](docs/TUTORIAL.md) takes 15 minutes.

## Edit a sequence without losing your annotations

![The sequence, the circular map and the feature table of pBAD30_eilR, side by side](docs/images/editor.png)

Insert, delete, replace, paste, rotate or reverse-complement, and every feature — including `join`, `complement`, fuzzy and
origin-crossing ones — stays on the bases it belongs to. Tests make hundreds of random edits and compare every base against an
independent model.

Also here: [unlimited undo with a plain-English History](docs/USER_GUIDE.md) ("Inserted 5 bp at 301: …"), a
[search bar](docs/USER_GUIDE.md) for DNA on both strands, IUPAC codes, fuzzy matches, protein motifs in six frames, feature names and
primers, a [protein view and 3D structures](docs/USER_GUIDE.md) from RCSB PDB, AlphaFold DB or ESMFold, the
[construct check-up with transcription units](docs/USER_GUIDE.md), [ORF detection](docs/USER_GUIDE.md), light and dark themes,
autosave and crash recovery. Reads and writes GenBank, reads `.dna` and FASTA.

## Design an assembly and see every junction

![The Assembly tab rebuilding pBAD30 from two fragments: a card per junction, each checked, and the designed primers on the sequence](docs/images/assembly_gibson.png)

Pick the fragments, and openDNA designs the primers and shows you each junction: the exact joining bases, their length, Tm and GC,
and whether the circle closes. Above, pBAD30 is rebuilt from its own backbone and its eilR insert — two junctions, four oligos.

Also here: [Golden Gate](docs/USER_GUIDE.md) with overhang sets checked against published high-fidelity data, and
[Sidewinder](docs/USER_GUIDE.md).

## Primers, PCR products and mis-priming

![The primer list of pBAD30 with Tm, matches and secondary structure, and the PCR product of its two verification primers](docs/images/primers_pcr.png)

Every primer in the file and in your own database, with its Tm, hairpins, self-dimers, binding sites and mismatches. Choose two and
see the product. Above: the verification primers of pBAD30 give a 1,281 bp product, and one primer is flagged for a self-dimer.

Also here: [primer design for a selection](docs/USER_GUIDE.md), [a check for other binding sites](docs/USER_GUIDE.md) in the
sequence or in a background genome, and a [restriction digest with a virtual gel](docs/USER_GUIDE.md).

## Find Cas9 guides and check their specificity

![The Spacer Finder listing 183 SpCas9 guides inside one gene of pFG051, with scores and off-target counts](docs/images/cas9_guides.png)

Every guide in a gene or a whole plasmid, with its cut site, GC, a checklist score (a rule of thumb, not a trained model) and the other
sites it matches — in this sequence and in a background genome. Add one as a feature and it is drawn with its protospacer, PAM and cut.

Also here: [SaCas9 and Cas12a](docs/USER_GUIDE.md), and the [off-target list](docs/USER_GUIDE.md) of any guide.

## Map reads and see what your clone really is

![The alignment view of pFG051: coverage, variants and the mapped reads, with a clear drop where a deletion is](docs/images/reads_coverage.png)

An IGV-like view with coverage, variants, mismatches and the reads themselves, beside the features and primers of the same sequence.
Above are simulated Illumina reads of a mixed colony; the dip near 1,900 is a deletion carried by part of the population.

Also here: [read mapping with bwa or minimap2](docs/USER_GUIDE.md) (on Windows through WSL), and
[GFF / GTF / BED / VCF / BAM tracks](docs/USER_GUIDE.md).

## Search one sequence inside another

![Local Align finding the araC gene of pBAD30 on the E. coli chromosome at 99 % identity, drawn as a track](docs/images/local_align.png)

A seed-and-extend search with gaps, on both strands, with the hits drawn as a track. It is not BLAST, and it needs no server: above,
the araC gene of pBAD30 is found on the *E. coli* chromosome at 99.0 % identity, inverted, with its 9 mismatches listed.

Also here: a [background genome](docs/USER_GUIDE.md) shared by the primer and guide tools, and
[import from NCBI](docs/USER_GUIDE.md) by accession or search.

## Whole genomes, as easily as a plasmid

![The map of Vibrio cholerae N16961: both chromosomes listed, the open one drawn with its genes](docs/images/multi_genome.png)

A multi-record file lists its sequences and a click opens one: above, both chromosomes of *Vibrio cholerae* N16961 (2.96 Mb and
1.07 Mb) in one GenBank file. A 4.6 Mb genome opens in the same window as a 5 kb plasmid.

Also here: [multi-record GenBank and FASTA](docs/USER_GUIDE.md), searching every record at once, and
[whole genomes from NCBI](docs/USER_GUIDE.md).

## Privacy

openDNA is offline. It reaches the internet only for the features you choose to use — NCBI import, 3D structure lookup, the update
check — and tells you before it does. [PRIVACY.md](PRIVACY.md) says exactly what is sent, and where settings and caches are kept.

## Documentation

| | |
|---|---|
| [docs/TUTORIAL.md](docs/TUTORIAL.md) | Hands-on walk-through: edit, annotate, design primers, assemble |
| [docs/USER_GUIDE.md](docs/USER_GUIDE.md) | Every feature, the Python API, the primer database, read mapping |
| [PRIVACY.md](PRIVACY.md) | What is sent over the network, what is stored where, how to delete it |
| [SECURITY.md](SECURITY.md) | Reporting a vulnerability; the unsigned exe and how to verify it |
| [CONTRIBUTING.md](CONTRIBUTING.md) | Bugs, feature requests, running from source, building, the tests |
| [docs/TESTING.md](docs/TESTING.md) | The test plan and its progress |
| [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md) | The components openDNA is built with and their licences |
| Help > Citations | What each panel is based on, and the packages used |

## Report a bug, ask for a feature

**[Open an issue](https://github.com/OpenDNA-project/openDNA/issues/new/choose)**, or use *Help > Report a Problem* in the program: it
fills the form with the version and the last log lines (sequences, file names and paths are blanked out first, and nothing is sent until
you press Submit). Security problems: [SECURITY.md](SECURITY.md), not a public issue.

Want to help with the code? [CONTRIBUTING.md](CONTRIBUTING.md) has how to run it from source, build the exe and run the tests.

## Licence

Copyright (C) 2026 Antoine Champie. openDNA is free software under the **GNU Affero General Public License v3.0**
([LICENSE](LICENSE)): you may use, study and change it, but anything you distribute, or run as a network service, built on it must be
offered under the same licence with its source. Third-party components keep their own licences, including Qt under the LGPL:
see [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).

**Forks are welcome, under another name**: the names "openDNA" and "oDNA" and the logo stay with this project
([TRADEMARKS.md](TRADEMARKS.md)).

openDNA builds on PySide6 / Qt, PyO3, `gb-io`, `noodles`, `quick-xml` and 3Dmol.js. The methods it implements are cited in
**Help > Citations**: Gibson et al. 2009, Engler et al. 2008, Potapov et al. 2018, SantaLucia 1998, the Cas9 / Cas12a papers and the
Sidewinder method (Robinson et al., Nature 2026).
