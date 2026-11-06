# The parts library: how it was built, what went wrong, and the curation rules to apply

Hand-off note for whoever works on `app/parts/` next. It covers the data model, the two build
methods used so far, an audit of what they actually produced (the bad parts are listed by name),
and the rules that should govern what gets in from now on.

Status at the time of writing: `app/parts/ncbi_vectors.tsv` (94 rows) and
`app/parts/engineering_proteins.tsv` (47 rows). Neither has been curated against the rules below.
**The first job is to apply them**, not to add more sources.


## 1. What the library is for

The library is the set of sequences that **Tools → Detect Common Features** looks for in an open
plasmid, and that the **Plasmid Checkup** uses to find promoters, terminators, RBS and origins it
can reason about. A part earns its place only if finding it in someone's plasmid tells them
something they did not already know.

That is a higher bar than "it was annotated in a GenBank record". It excludes anything specific to
one construct, anything the program already handles another way, and anything so short that it
matches by chance.

### The matcher it feeds

`annotation.detect()` → `find_parts()` (Rust, `core/src/annotate.rs`). Defaults, in
Tools → Settings → Detect Features:

| Setting | Default | Meaning |
|---|---|---|
| `min_dna` | 95 % | share of bases that must match, no gaps |
| `min_protein` | 90 % | coding parts are also searched as protein, to catch codon-optimised copies |
| `min_length` | 15 bp | parts shorter than this are skipped entirely |

Two consequences the curation rules below depend on:

- **Protein matching makes near-variants redundant.** A part only needs the parent sequence; a
  codon-optimised or point-mutated version of a CDS still matches it. Storing sfGFP *and* EGFP *and*
  GFP as separate rows buys very little and costs a duplicate hit.
- **`min_length` of 15 bp is far too permissive.** A 17 bp part matches at random roughly once per
  4^17 ≈ 1.7 × 10^10 bases in theory, but real plasmid sequence is not random — it is full of shared
  vector backbone, and 95 % identity over 17 bp means *one* mismatch. Short parts are the main
  source of noise.


## 2. The data model

One file per source in `app/parts/`, plus the user's own `parts.tsv` in the app data folder.
`app/parts_library.py` loads them; the Library… dialog switches any of them off everywhere
(Detect Common Features and the Plasmid Checkup read the same loader).

```
#title: Engineering parts (UniProt + GenBank)
#reference: UniProt Consortium. ... Nucleic Acids Res 53, D609-D617 (2025).
#doi: 10.1093/nar/gkae1010
#note: Swiss-Prot entries, CC BY 4.0; DNA from the GenBank record each entry cites.
name<TAB>type<TAB>sequence<TAB>organism<TAB>strength<TAB>units<TAB>source
```

- `strength` is the number the source measured, `units` says what it means (TS, TE, RPU…).
  **Units are never comparable across sources**, which is why the unit travels with the number.
  A part nobody measured leaves both empty — never estimate one.
- `organism` is the host the measurement was made in.
- `source` must identify the exact record, e.g. `AE004092.2 (AAK33936.1), UniProt Q99ZW2`.
- Parts are de-duplicated **by sequence** at load time, last file wins (user's own is last).
  Files sort alphabetically, so `engineering_proteins` currently loses to `ncbi_vectors` for an
  identical sequence — see defect D5.


## 3. Method A — annotated features of NCBI vector records

`scripts/build_parts.py`. Fetches ~26 cloning-vector GenBank records by accession and turns every
annotated feature of a useful type into a part: name from `/label`, `/gene`, `/product` or a short
`/note`; sequence in the feature's own orientation; source = vector name + accession.

**Why it was chosen:** NCBI content is public domain, the records are real, and attribution is a
single accession.

**Why it produced junk:** it trusts the record's own annotation completely. GenBank submitters
annotate *their construct*, not a parts catalogue. They label primer binding sites, polylinkers,
and whichever 34 bases of an FRT site happened to survive their cloning. The script has a `JUNK`
substring list (`"derived from"`, `"fragment"`, `"region"`, …) and an `ALIASES` table, but these
are far too weak — see the audit in §5.

**It also does not deduplicate across records.** Six vectors carrying KanR give six KanR rows.


## 4. Method B — UniProt entry → the GenBank record it cites

`scripts/build_engineering_parts.py`. For each wanted part, a UniProtKB query (reviewed/Swiss-Prot
preferred, TrEMBL as fallback — both CC BY 4.0) returns the canonical protein and its EMBL/GenBank
cross-references. The coding sequence is then pulled from that record through the protein
accession's `/coded_by`, and kept only if it passes two gates:

1. **Length gate** — the DNA must be exactly `3 × protein_length + 3` (or `+ 0`). A cross-reference
   pointing at a partial CDS or the wrong feature comes back the wrong length and is refused.
2. **Translation check** (`scratchpad/verify_parts.py`) — multiple of 3, a start codon, exactly one
   stop and it is the last codon, no ambiguous bases.

**This worked.** All 47 rows pass both gates. The gates caught three real errors that would
otherwise have shipped silently:

| Part | Fetched | Correct | Cause |
|---|---|---|---|
| SacB | 204 bp | 1,422 bp | cross-reference to a partial CDS |
| MBP (malE) | 144 bp | 1,191 bp | same |
| "Tn5 transposase" | 3,048 bp | — | **P03008 is the Tn3 transposase.** The accession was wrong, not the fetch |

The third is the important lesson: **the length gate cannot catch a confidently wrong accession.**
P03008 returned a perfectly valid CDS of the wrong gene. Only reading the returned protein name
caught it. Any accession written from memory must be checked against the name the API returns.

**Where it falls down:** Swiss-Prot holds natural proteins. Engineered variants (sfGFP, mScarlet,
dCas9, SB100X, piggyBac) are simply absent, and entries that cross-reference only mRNA or a
polyprotein yield no clean CDS. 11 wanted parts were not obtained: Flp, Bxb1 integrase, Tn5 and Tn7
transposases, T3 and SP6 RNA polymerases, StCas9, DsRed, Renilla luciferase, TEV protease.


## 5. Audit — what is actually wrong with the current files

This is the user's complaint, made specific. **These are the things to fix.**

### D1. Primer binding sites — must be removed outright

```
pUC/M13 forward primer binding site    17 bp   pGEM-luc (X65316)
pUC/M13 reverse primer binding site    17 bp   pGEM-luc (X65316)
```

Primers are not parts. The program has a whole primer subsystem — a primer database, binding-site
matching with mismatch tolerance, and its own drawing in every view. Annotating "M13 forward" as a
*feature* duplicates that badly and clutters the map. **No primer, primer site or sequencing site
belongs in the parts library, ever.**

### D2. The same biological part stored many times

| Part | Rows | Lengths | Comment |
|---|---|---|---|
| KanR | 6 | 795 ×5, 816 | one canonical entry is enough |
| AmpR | 3 (+1 in the other file) | 861 | identical sequence, three records |
| lacZα | 3 | 231, 285, 366 | three different *truncations* of the same gene |
| R6K gamma origin | 3 | 387, 399, 400 | same origin, different boundaries |
| p15A origin | 3 | 545, 546, 913 | also three different names |
| ColE1/pMB1 origin | 3 | 281, 595, 698 | named "pBR322 origin" ×2 and "pMB1 replication origin" |
| FRT site | 5 | 34, 35, 42, 46, 48 | see D3 |
| multiple cloning site | 4 | 34, 39, 51, 63 | see D4 |
| rrnB terminator | 3 | 260, 425, 426 | one is named "5S RNA, T1, T2, rrnB" |
| sacB | 2 | 1,422 | identical |
| CmR | 2 | 639, 660 | |
| lac promoter | 3 | 30, 31, 55 | "Lac promoter", "E. coli lac operon promoter", "PlacZ; lacZ promoter" |
| Pc promoter | 2 | 29 | "Pc" and "Pc promoter" |

Identical sequences collapse at load time, so they waste space rather than produce double hits. The
*different-length* ones are worse: they are genuinely separate entries that will all match the same
stretch of a plasmid, and `detect()` then has to pick between them.

### D3. Fragments of a part stored as if they were parts

```
natural FRT site                        48 bp
FRT site                                46 bp
distal 42-nt of natural FRT site        42 bp
distal 35-nt of natural FRT site        35 bp
distal 34-nt of of natural FRT site     34 bp     <- note the "of of" typo, copied straight from the record
AmpR  (sig_peptide)                     69 bp     <- the signal peptide of a CDS already present in full
```

An FRT site is 34 bp of core plus spacer. Storing five overlapping slices of it means a plasmid with
one FRT site lights up five times. Keep **one** canonical FRT. The `sig_peptide` row is a fragment
of the AmpR CDS stored two rows above it.

The `of of` typo also shows that record `/note` text is being copied into user-visible names with no
cleaning at all.

### D4. Cloning scaffolding and construct-specific oddities

```
multiple cloning site / multiple cloning sites ×2 / pUC18 MCS; polylinker
H-strand Y effector site        157 bp    pBR322 curiosity
L-strand Y effector site         64 bp    pBR322 curiosity
SITS  (regulatory)              114 bp    unexplained acronym
prescission cleavage site        24 bp    protease site, too short to match safely
operator O1 / O2 / I2 + I1       22 / 18 / 39 bp
```

An MCS is defined by the enzymes in it, which the program already computes and draws (unique
cutters, the digest view). It is different in every vector, so it cannot be a library part. The
pBR322 "Y effector sites" are annotations almost nobody wants. The araBAD operators are real
biology but 18–39 bp, i.e. in the chance-match zone, and they are already implied by the araBAD
promoter next to them.

### D5. Cross-collection duplicates

Eight sequences are byte-identical in both files:

```
rop == Rop / Rom          exonuclease == Lambda Red exo      AmpR == AmpR (bla, TEM-1)
GmR == GmR (aacC1)        araC == AraC                       gam protein == Lambda Red gam
KanR == KanR (nptII, Tn5) bet protein == Lambda Red beta
```

The loader keeps one, and because files sort alphabetically `ncbi_vectors` wins — so the *worse*
row wins. The `engineering_proteins` row has the organism, a UniProt accession and a proper name;
the `ncbi_vectors` one says `bet protein` and cites only a plasmid. Either resolve the overlap by
curating the duplicate out of `ncbi_vectors`, or give the format an explicit precedence key so the
richer collection wins.

### D6. Ambiguous names

`T1` appears twice with different sequences (44 bp from pEX18Gm, 105 bp from pSEVA221): two
different terminators sharing a name. `T0`, `T2`, `tL3`, `Pc` have the same problem — these are
*positions in a particular vector*, not part names. They need qualifying (`T1 (rrnB)`,
`T0 (λ tL3)`) or dropping.


## 6. Curation rules — apply these

### Inclusion

A row belongs in the library only if **all** of these hold.

1. **It is a reusable biological element.** Someone could deliberately put it in a new construct and
   name it. Test: would it appear in a parts catalogue, or only in one lab's plasmid map?
2. **Finding it is informative.** It tells the user something the program does not already show.
3. **It is canonical and complete.** The whole promoter, the whole CDS, the whole terminator — not a
   slice that one record happened to annotate.
4. **It is long enough to be specific** (see the table below).
5. **It is traceable.** `source` names a record and accession a reader can check.
6. **Its name is the name people use.** `AmpR`, not `beta-lactamase TEM-1 precursor, partial cds`.

### Exclusion — never include

- Primer, sequencing-primer or probe binding sites (D1).
- Multiple cloning sites, polylinkers, spacers, linkers, stuffer fragments.
- Anything whose name contains `fragment`, `partial`, `region`, `derived from`, `distal`, `-nt of`,
  `site of`, or a bare coordinate.
- Signal peptides, domains or other sub-parts of a CDS the library already holds in full.
- Restriction sites (computed), ORFs (computed), transcription units (computed by the Checkup).
- Anything whose function you cannot state in one line (`SITS`).

### Minimum lengths — replace the single 15 bp setting

The global `min_length` is the wrong tool because the right threshold depends on the type. Suggested
per-type floors, enforced at build time rather than at search time:

These are `FLOORS` in `scripts/build_parts.py`, enforced at build time by `curate()`:

| Type | Floor | Why |
|---|---|---|
| CDS | 150 bp | shorter coding parts are peptide tags; handle those separately if wanted |
| promoter | 19 bp | **not 30.** A phage promoter really is short — T7 is 19 bp, SP6 20 — and an exact 19-mer is specific enough. The chance matches came from 17–18 bp primer sites and generic operators, which the name rules remove anyway |
| terminator | 30 bp | stem-loop plus T-tract |
| rep_origin | 200 bp | origins are long; a short one is a fragment |
| oriT | 50 bp | |
| protein_bind | 25 bp | lacO, tetO, araO — named, never generic |
| misc_recomb | 25 bp | loxP 34, FRT 34+spacer, attB/attP |
| polyA_signal | 50 bp | |
| centromere | 100 bp | |
| misc_feature | — | **banned.** Give it a real type or drop it |

`misc_feature` is the tell: 22 of the 94 NCBI rows are `misc_feature`, and nearly every defect above
is one of them. A part whose type you cannot name is a part you have not understood.

### One part, one row

For each biological part keep exactly one row. **Which row matters more than it looks**, because the
matcher needs 95 % identity across the whole part with no gaps, so *a longer part is the harder one
to find*. A 426 bp rrnB taken from pBAD18 does not match a plasmid that carries 260 bp of it.

- **CDS** — keep the **longest**: that is the complete gene, and a partial one is a fragment.
- **Everything else** — keep the **shortest** surviving row. Different records annotate different
  extents of one element, and after the junk rules have run, the shortest is its conserved core,
  which is found in more constructs.

This was learnt the hard way: keeping the longest row of each origin and terminator silently stopped
the Plasmid Checkup finding rrnB on a demo plasmid it had always found.

Record the alternatives nowhere — the matcher's protein mode and 95 % threshold cover the variants.


## 7. Worklist

**Done** (defects D1–D6 are all closed; §9 records the result):

1. ~~Curate `ncbi_vectors.tsv` down~~ — 94 → 31 rows, as `curate()` in `scripts/build_parts.py`, so a
   rebuild applies the same rules. `--recurate` re-runs them over the file without fetching.
2. ~~Resolve the D5 overlap~~ — the 16 CDS that `engineering_proteins` holds with a UniProt accession
   and an organism are now dropped from `ncbi_vectors` (`SUPERSEDED`). No sequence is in both files.
3. ~~Enforce the per-type floors~~ — `FLOORS`, applied in both build scripts.

**Still to do:**

4. **Add the parts §4 could not fetch**, individually: Flp, Bxb1 integrase, Tn5 and Tn7 transposases,
   T3 and SP6 RNA polymerases, StCas9, DsRed, Renilla luciferase, TEV protease. Each needs a named
   record and the same two gates. Check the returned protein name against what was asked for —
   remember P03008.
5. **Build the short-DNA collection** (promoters, terminators, loxP/FRT/attB, origins). These have no
   protein to check against, so each one needs a cited record and an explicit length and boundary
   decision. This is the collection where `strength` can finally be filled in, and where the rules
   above matter most, because short parts are exactly what goes wrong.
6. **Restore what curation removed, properly.** Dropping a row was the right call when its identity
   could not be established, but several are real parts that belong in the collection of §5 once
   someone decides their boundaries against a cited record: the T-DNA right border (25–26 bp, below
   the `misc_recomb` floor), the araBAD operators O1/O2/I1/I2, and the pEX18Gm/pSEVA terminators
   labelled only `T0`/`T1`/`T2`.

### How to verify any change

```bash
# every CDS really is one
python scratchpad/verify_parts.py

# nothing is duplicated within or across collections, and no row breaks the floors
python -m pytest tests/test_parts_library.py -q
```

`tests/test_parts_library.py` checks the format and the per-source switch, and now also holds the
curation line against the files that actually ship: `test_no_part_is_junk` (the `NEVER` patterns, and
no `misc_feature`/`regulatory`/`sig_peptide`), `test_every_part_clears_the_floor_for_its_type`,
`test_one_part_one_row` (no repeated name, no sequence stored twice in or across collections) and
`test_every_part_says_where_it_came_from`. They import the rules from `scripts/build_parts.py`, so
there is one source of truth. **Adding a part that breaks a rule now fails the suite.**


## 9. What the curation actually removed

`ncbi_vectors.tsv`: **94 → 31 rows.** 21 parts kept their sequence but took a canonical name
(`p15A Origin`/`p15A origin of replication`/`p15A_ori` → `p15A ori`; `rrnB`/`rrnB transcription
terminator`/`5S RNA, T1, T2, rrnB` → `rrnB T1T2 terminator`; `natural FRT site` → `FRT`;
`yeast centromere` → `CEN (yeast centromere)`, retyped from `misc_feature` to `centromere`).
63 rows went: both primer binding sites, all four cloning sites, the three `distal …-nt of` FRT
slices, the `AmpR` signal peptide, the two pBR322 "Y effector sites", `SITS`, the prescission site,
the three generic araBAD operators, the 16 CDS superseded by `engineering_proteins`, the bare
`T0`/`T1`/`T2` labels, and every duplicate row of a part kept elsewhere.

`engineering_proteins.tsv`: **47 → 46 rows.** The two chloramphenicol acetyltransferases
(`CmR (cat)` from *Morganella morganii*, `CAT reporter (cat, Tn9)` from *E. coli*) were one part
under two names; the Tn9 *E. coli* allele that vectors actually carry is kept, as `CmR (cat, Tn9)`.

`app/editor.py` gained colours for `centromere` and `telomere`, which the retyping made reachable.

The library is now **77 parts**: 31 vector features and 46 verified coding parts, no part twice.

Two tests had to change, and both were asserting something the curation made untrue:

- `test_gui.py::test_gene_groups_and_detect_common_features` looked parts up as `"Lac promoter"` and
  `"AmpR"`; those are now `"lac promoter"` and `"AmpR (bla, TEM-1)"`.
- `test_checkup.py::test_library_parts_on_a_real_plasmid` relied on a 27 bp row called `T2` — a
  fragment of rrnB T1T2, under the terminator floor and ambiguously named. With it gone, **the
  library contributes nothing at all to pBAD30**, because everything it would find is already
  annotated in the file. That is the right behaviour, not a regression: the library exists for what a
  file leaves out. The test now asserts exactly that, and checks the library still supplies the
  terminator when the same sequence is stripped of its annotations.


## 8. Sources: what is usable, what is not

Established while choosing sources; worth not re-deriving.

| Source | Verdict |
|---|---|
| **UniProt REST** (`rest.uniprot.org`) | **Use.** No auth, CC BY 4.0, cross-references to the nucleotide record |
| **NCBI E-utilities** | **Use.** Public domain. ≤3 requests/s anonymously, so space them ~350 ms |
| **NCBI AMRFinderPlus** (`AMR_CDS.fa`) | Usable, public domain, but clinical AMR — ~40 of 11,554 rows are relevant |
| **ResFinder** | Usable, Apache-2.0, acquired resistance only |
| **pLannotate** | **Avoid.** 1,367 features covering everything wanted, but the data is SnapGene-derived with no data licence, inside a GPL-3.0 repo. Shipping a commercial editor's curated feature list inside this AGPL program is the exact risk `docs/PUBLISH_AUDIT.md` was written about. **No REST API either** — it is a CLI/library plus a Streamlit page |
| **CARD** | **Cannot use.** README bars use by any commercial organisation without a McMaster licence |
| **Cello UCF** (Voigt/CIDAR) | **Cannot use.** Sequences + response-function parameters in one JSON, but the repo carries no LICENSE at all |
| **Cambray 2013** (NAR, doi 10.1093/nar/gkt163) | **Cannot bundle.** ~76 *named* terminators with measured TE — exactly the right shape — but CC BY-**NC** 3.0, incompatible with this program's AGPL-3.0 |
| **La Fleur 2022** (Nat Commun, PMC9440211) | Open (CC BY 4.0) and numeric, but the rows are unnamed random 40-mers for model training, not reusable parts |
| **Chen 2013** (582 terminators, doi 10.1038/nmeth.2515) | Not in PMC, not open access. No lawful bulk route found |
| **iGEM Registry** | CloudFront returns 403 to scripted requests |
| **SynBioHub** | REST returns "Login required" |
| **Addgene** | Developer API needs an approved token and an accepted data licence (~5 business days) |
| **NCBI free-text search** | **Does not work for curation.** "plasmid RP4 complete sequence" returns the *Streptomyces coelicolor* genome; "transposon Tn7" returns an *E. coli* chromosome. Use accessions you have verified, never a name search |

The licence constraint that governs all of this: **openDNA is AGPL-3.0**, so any bundled data
must permit commercial use and redistribution. That rules out NC and non-commercial-research
licences no matter how good the data is.
