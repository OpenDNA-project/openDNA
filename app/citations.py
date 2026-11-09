"""Help > Citations: what each panel is based on, and the packages openDNA is built with.

Whenever a panel uses a published method or dataset, add it to PANELS below (and say so where it is
used): the user asked for every such source to be cited by panel.
"""
from html import escape

from PySide6.QtWidgets import QDialog, QDialogButtonBox, QTextBrowser, QVBoxLayout

# panel -> [(reference, doi or None, what it is used for)]
PANELS = {
    "Assembly: Gibson / HiFi": [
        ("Gibson DG, Young L, Chuang RY, Venter JC, Hutchison CA 3rd, Smith HO. Enzymatic assembly of DNA molecules up to "
         "several hundred kilobases. Nat Methods 6, 343-345 (2009).", "10.1038/nmeth.1318",
         "The method: neighbouring fragments share end homology that the exonuclease exposes and the polymerase and ligase seal."),
    ],
    "Assembly: Golden Gate": [
        ("Engler C, Kandzia R, Marillonnet S. A one pot, one step, precision cloning method with high throughput capability. "
         "PLoS ONE 3(11), e3647 (2008).", "10.1371/journal.pone.0003647",
         "The method: Type IIS sites outside the fragment leave 4 nt overhangs that are ligated in one pot."),
        ("Potapov V, Ong JL, Kucera RB, et al. Comprehensive profiling of four base overhang ligation fidelity by T4 DNA ligase "
         "and application to DNA assembly. ACS Synth Biol 7(11), 2665-2674 (2018).", "10.1021/acssynbio.8b00333",
         "The 25 overhangs openDNA prefers and checks against (supplementary Tables S3 and S4 of the bioRxiv version, "
         "doi:10.1101/322297: the junctions of the 12- and 24-fragment lac cassettes). The paper's full ligation matrix is not "
         "included, so the one-base-apart warning is a rule of thumb, not a measured fidelity."),
    ],
    "Assembly: Sidewinder": [
        ("Robinson NE, Zhang W, Ghosh R, Gerber B, Zhang H, Sanfiorenzo C, Wang S, Di Carlo D, Wang K. Construction of complex "
         "and diverse DNA sequences using DNA three-way junctions. Nature 651, 491-500 (2026).", "10.1038/s41586-025-10006-0",
         "The Sidewinder method (three-way junctions from toeholds and barcodes), in a simplified form: check the strand "
         "layout against the paper before ordering oligos."),
    ],
    "Spacer Finder (Cas9 / Cas12a guides)": [
        ("Jinek M, Chylinski K, Fonfara I, Hauer M, Doudna JA, Charpentier E. A programmable dual-RNA-guided DNA endonuclease in "
         "adaptive bacterial immunity. Science 337, 816-821 (2012).", "10.1126/science.1225829",
         "Cas9 cuts a site next to an NGG PAM, 3 bp from the PAM, inside the 20 nt the guide pairs with."),
        ("Hsu PD, Scott DA, Weinstein JA, et al. DNA targeting specificity of RNA-guided Cas9 nucleases. Nat Biotechnol 31, "
         "827-832 (2013).", "10.1038/nbt.2647",
         "Why other sites with a few mismatches are counted, and why those next to the PAM (the seed) matter most."),
        ("Ran FA, Cong L, Yan WX, et al. In vivo genome editing using Staphylococcus aureus Cas9. Nature 520, 186-191 (2015).",
         "10.1038/nature14299", "SaCas9: an NNGRRT PAM and a 21 nt guide."),
        ("Zetsche B, Gootenberg JS, Abudayyeh OO, et al. Cpf1 is a single RNA-guided endonuclease of a class 2 CRISPR-Cas system. Cell 163, "
         "759-771 (2015).", "10.1016/j.cell.2015.09.038",
         "Cas12a (Cpf1): a TTTV PAM on the 5' side of the target, and a staggered cut (after base 18 of the PAM's strand and base 23 of the "
         "other) leaving a 5 nt 5' overhang. openDNA uses a 23 nt protospacer; check the cut positions for the Cas12a you use "
         "(AsCas12a, LbCas12a and others differ a little)."),
        ("(no paper) The score is a plain checklist of common design advice (GC content, no TTTT, no long runs), not a trained "
         "model such as Doench's Rule Set 2.", None, "Use it to discard poor guides, not to rank good ones against each other."),
    ],
    "Primers and melting temperature (Primers, Primer Design, PCR, Assembly)": [
        ("SantaLucia J Jr. A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor thermodynamics. "
         "Proc Natl Acad Sci USA 95, 1460-1465 (1998).", "10.1073/pnas.95.4.1460",
         "Nearest-neighbour parameters for the Tm (the same model as Primer3 and Biopython's Tm_NN)."),
        ("von Ahsen N, Wittwer CT, Schutz E. Oligonucleotide melting temperatures under PCR conditions: nearest-neighbor "
         "corrections for Mg2+, deoxynucleotide triphosphate, and dimethyl sulfoxide concentrations. Clin Chem 47, 1956-1961 (2001).",
         "10.1093/clinchem/47.11.1956", "Converting Mg2+ to an equivalent Na+ in the salt correction."),
    ],
    "Protein view": [
        ("Rice P, Longden I, Bleasby A. EMBOSS: the European Molecular Biology Open Software Suite. Trends Genet 16, 276-277 (2000).",
         "10.1016/S0168-9525(00)02024-2",
         "The pKa values of the isoelectric point and charge (those of EMBOSS's iep). The mass is the sum of average residue masses "
         "(as on ExPASy) plus water; both are estimates for the free chain."),
        ("Kyte J, Doolittle RF. A simple method for displaying the hydropathic character of a protein. J Mol Biol 157, 105-132 (1982).",
         "10.1016/0022-2836(82)90515-0", "The hydrophobicity colours of the letters and of the 3D view (the Kyte-Doolittle scale)."),
        ("Zamyatnin AA. Protein volume in solution. Prog Biophys Mol Biol 24, 107-123 (1972).", "10.1016/0079-6107(72)90005-3",
         "The size colours (the volumes of the amino acids)."),
        ("Waterhouse AM, Procter JB, Martin DMA, Clamp M, Barton GJ. Jalview Version 2 - a multiple sequence alignment editor and "
         "analysis workbench. Bioinformatics 25, 1189-1191 (2009).", "10.1093/bioinformatics/btp033",
         "The property colours (Jalview's Zappo scheme: aliphatic, aromatic, positive, negative, polar, P and G, C)."),
    ],
    "Import from NCBI": [
        ("Entrez Programming Utilities Help [Internet]. Bethesda (MD): National Center for Biotechnology Information (US); 2010-.", None,
         "The E-utilities (esearch, esummary, efetch) that find and download GenBank / RefSeq records and proteins "
         "(https://www.ncbi.nlm.nih.gov/books/NBK25501/). Records are public NCBI data; check a record's own terms before reusing it."),
    ],
    "Background genomes from NCBI (Import > DNA > Background Genome from NCBI)": [
        ("O'Leary NA, Cox E, Holmes JB, et al. Exploring and retrieving sequence and metadata for species across the tree of life "
         "with NCBI Datasets. Sci Data 11, 732 (2024).", "10.1038/s41597-024-03571-y",
         "The NCBI Datasets v2 API (https://www.ncbi.nlm.nih.gov/datasets/docs/v2/api/) that lists assemblies and downloads their "
         "sequences (FASTA) and annotation (GFF3)."),
    ],
    "Features Database (Databases > Features Database)": [
        ("Built from the annotated features of public NCBI GenBank records of well-known cloning vectors, one entry per part with "
         "its source record.", None,
         "NCBI says most database content is free to use, but a submitter may claim rights in a record; the library keeps names, short feature "
         "sequences and their record numbers, not whole records. Your own features (Local: Add to Features Database) are yours."),
        ("UniProt Consortium. UniProt: the Universal Protein Knowledgebase in 2025. Nucleic Acids Res 53, D609-D617 (2025).",
         "10.1093/nar/gkae1010",
         "The engineering parts collection: markers, counter-selection genes, recombinases, CRISPR "
         "nucleases, polymerases, reporters and regulators. UniProtKB (CC BY 4.0) names each protein and cites the record holding its coding "
         "sequence; the DNA is then taken from that GenBank record (public domain), and kept only when its length matches the protein exactly."),
    ],
    "Plasmid Checkup (Tools > Plasmid Checkup)": [
        ("Harley CB, Reynolds RP. Analysis of E. coli promoter sequences. Nucleic Acids Res 15(5), 2343-2361 (1987).",
         "10.1093/nar/15.5.2343",
         "The sigma-70 consensus the predicted promoters are scored against: -35 TTGACA, a 17 bp spacer (15-19 accepted), "
         "-10 TATAAT."),
        ("d'Aubenton Carafa Y, Brody E, Thermes C. Prediction of rho-independent Escherichia coli transcription terminators. "
         "A statistical analysis of their RNA stem-loop structures. J Mol Biol 216(4), 835-858 (1990).",
         "10.1016/S0022-2836(05)80005-0",
         "What a predicted intrinsic terminator looks like: a G-C rich stem-loop followed by a T tract. openDNA uses a "
         "simplified version (perfect stem with G-T pairs allowed, a pair score instead of a free energy)."),
        ("Shine J, Dalgarno L. The 3'-terminal sequence of Escherichia coli 16S ribosomal RNA: complementarity to nonsense "
         "triplets and ribosome binding sites. Proc Natl Acad Sci USA 71(4), 1342-1346 (1974).", "10.1073/pnas.71.4.1342",
         "The Shine-Dalgarno sequence (AGGAGG) looked for before each start codon when no RBS is annotated."),
    ],
    "3D structure": [
        ("Jumper J, Evans R, Pritzel A, et al. Highly accurate protein structure prediction with AlphaFold. Nature 596, 583-589 (2021).",
         "10.1038/s41586-021-03819-2", "AlphaFold models, with the pLDDT confidence in the B-factor column."),
        ("Varadi M, Anyango S, Deshpande M, et al. AlphaFold Protein Structure Database: massively expanding the structural coverage of "
         "protein-sequence space with high-accuracy models. Nucleic Acids Res 50, D439-D444 (2022).", "10.1093/nar/gkab1061",
         "AlphaFold DB (alphafold.ebi.ac.uk): models looked up by UniProt accession; the confidence colours are its."),
        ("Lin Z, Akin H, Rao R, et al. Evolutionary-scale prediction of atomic-level protein structure with a language model. Science 379, "
         "1123-1130 (2023).", "10.1126/science.ade2574", "ESMFold, run by the ESM Metagenomic Atlas API: a prediction from the sequence alone."),
        ("Berman HM, Westbrook J, Feng Z, et al. The Protein Data Bank. Nucleic Acids Res 28, 235-242 (2000).", "10.1093/nar/28.1.235",
         "RCSB PDB (rcsb.org): experimental structures found by sequence search (90 % identity or more) and downloaded."),
        ("Rego N, Koes D. 3Dmol.js: molecular visualization with WebGL. Bioinformatics 31, 1322-1324 (2015).",
         "10.1093/bioinformatics/btu829", "The 3D page in the browser (3Dmol.js, BSD-3-Clause, bundled in app/vendor)."),
        ("Labesse G, Colloc'h N, Pothier J, Mornon JP. P-SEA: a new efficient assignment of secondary structure from C alpha trace of "
         "proteins. Comput Appl Biosci 13, 291-295 (1997).", "10.1093/bioinformatics/13.3.291",
         "Helices and strands of the 3D view of the Protein panel, from the alpha carbons alone (its distance criteria)."),
    ],
    "Genetic code (translation, ORFs, protein search)": [
        ("Elzanowski A, Ostell J. The Genetic Codes. NCBI Taxonomy, https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi "
         "(a continuously updated web page: no year or DOI given).", None,
         "The translation tables 1-6, 9-16, 21-26, 29, 30, 32 and 33 (names, codons, start codons), embedded in core/src/code.rs "
         "and checked codon by codon against Biopython's copy of the NCBI tables. Tables 27, 28 and 31, where a codon is a stop or an "
         "amino acid depending on its place in the gene, are not offered."),
    ],
}

PACKAGES = [
    ("PySide6 / Qt for Python", "the user interface", "LGPL"),
    ("PyO3", "Rust-to-Python bindings for the core", "MIT / Apache-2.0"),
    ("maturin", "builds the Rust core as a Python module", None),
    ("gb-io", "GenBank reading and writing (a patched copy is in core/vendor/gb-io)", "MIT"),
    ("noodles", "BAM / VCF / index reading", "MIT"),
    ("quick-xml", "reading .dna files", "MIT"),
    ("3Dmol.js", "the 3D structure page (bundled in app/vendor)", "BSD-3-Clause"),
    ("PyInstaller", "builds the standalone .exe", None),
    ("pytest, Biopython, pyflakes", "tests and checks only (not shipped)", None),
]


def citations_html():
    out = ["<h2>Citations</h2><p>What each part of openDNA is based on.</p>"]
    for panel, refs in PANELS.items():
        out.append(f"<h3>{escape(panel)}</h3><ul>")
        for text, doi, use in refs:
            link = f" <a href='https://doi.org/{doi}'>doi:{doi}</a>" if doi else ""
            out.append(f"<li>{escape(text)}{link}<br><i>{escape(use)}</i></li>")
        out.append("</ul>")
    out.append("<h3>Packages</h3><ul>")
    for name, what, licence in PACKAGES:
        out.append(f"<li><b>{escape(name)}</b>: {escape(what)}" + (f" ({escape(licence)})" if licence else "") + "</li>")
    out.append("</ul>")
    return "".join(out)


class CitationsDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Citations")
        self.resize(760, 640)
        self.view = QTextBrowser(openExternalLinks=True)
        self.view.setHtml(citations_html())
        buttons = QDialogButtonBox(QDialogButtonBox.Close)
        buttons.rejected.connect(self.reject)
        box = QVBoxLayout(self)
        box.addWidget(self.view)
        box.addWidget(buttons)
