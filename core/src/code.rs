//! Genetic codes (NCBI translation tables) and the translation of annotated CDS features.
//!
//! Rust primer for readers new to the language: docs/RUST_PRIMER.md. Short version of what is used here:
//! - `&'static str` is text that lives for the whole program (string literals in the source).
//! - `Option<T>` is "a T or nothing" (`Some(x)` / `None`); `Result<T, E>` is "a T or an error" (`Ok` / `Err`).
//! - `?` after a Result/Option returns early from the function when it is an Err/None.
//! - `match` picks code by the shape of a value; `if let` is a one-arm match.
//!
//! Data source: the tables below are the NCBI genetic codes (Elzanowski & Ostell, NCBI Taxonomy "Genetic Codes",
//! https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi, file gc.prt), copied through Biopython's
//! `Bio.Data.CodonTable.unambiguous_dna_by_id`, which is also what the tests compare every codon against.
//! Each table is written the way NCBI writes it: 64 letters, codons in the order TCAG x TCAG x TCAG (TTT, TTC, TTA, TTG,
//! TCT, ... GGG), the first letter = first codon base.
//!   - `aas`:    the amino acid of each codon, '*' = stop (NCBI "ncbieaa").
//!   - `starts`: 'M' where the codon can start a protein, '-' elsewhere (NCBI "sncbieaa").
//! Left out on purpose: tables 27, 28 and 31 (Karyorelict, Condylostoma, Blastocrithidia nuclear codes), where the
//! meaning of a codon (stop or amino acid) depends on its position in the gene: one letter per codon cannot say that.

use gb_io::seq::{Feature, Location};

use crate::primers::revcomp;

/// One genetic code.
pub struct GeneticCode {
    pub id: u8,
    pub name: &'static str,
    pub aas: &'static str,
    pub starts: &'static str,
}

/// `const` item: fixed at compile time. A slice (`&[...]`) is a read-only view of an array.
pub const TABLES: &[GeneticCode] = &[
    GeneticCode { id: 1, name: "Standard", aas: "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", starts: "---M---------------M---------------M----------------------------" },
    GeneticCode { id: 2, name: "Vertebrate Mitochondrial", aas: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG", starts: "--------------------------------MMMM---------------M------------" },
    GeneticCode { id: 3, name: "Yeast Mitochondrial", aas: "FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG", starts: "----------------------------------MM---------------M------------" },
    GeneticCode { id: 4, name: "Mold, Protozoan, Coelenterate Mitochondrial and Mycoplasma/Spiroplasma", aas: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", starts: "--MM---------------M------------MMMM---------------M------------" },
    GeneticCode { id: 5, name: "Invertebrate Mitochondrial", aas: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG", starts: "---M----------------------------MMMM---------------M------------" },
    GeneticCode { id: 6, name: "Ciliate, Dasycladacean and Hexamita Nuclear", aas: "FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", starts: "-----------------------------------M----------------------------" },
    GeneticCode { id: 9, name: "Echinoderm and Flatworm Mitochondrial", aas: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG", starts: "-----------------------------------M---------------M------------" },
    GeneticCode { id: 10, name: "Euplotid Nuclear", aas: "FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", starts: "-----------------------------------M----------------------------" },
    GeneticCode { id: 11, name: "Bacterial, Archaeal and Plant Plastid", aas: "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", starts: "---M---------------M------------MMMM---------------M------------" },
    GeneticCode { id: 12, name: "Alternative Yeast Nuclear", aas: "FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", starts: "-------------------M---------------M----------------------------" },
    GeneticCode { id: 13, name: "Ascidian Mitochondrial", aas: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG", starts: "---M------------------------------MM---------------M------------" },
    GeneticCode { id: 14, name: "Alternative Flatworm Mitochondrial", aas: "FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG", starts: "-----------------------------------M----------------------------" },
    GeneticCode { id: 15, name: "Blepharisma Macronuclear", aas: "FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", starts: "-----------------------------------M----------------------------" },
    GeneticCode { id: 16, name: "Chlorophycean Mitochondrial", aas: "FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", starts: "-----------------------------------M----------------------------" },
    GeneticCode { id: 21, name: "Trematode Mitochondrial", aas: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG", starts: "-----------------------------------M---------------M------------" },
    GeneticCode { id: 22, name: "Scenedesmus obliquus Mitochondrial", aas: "FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", starts: "-----------------------------------M----------------------------" },
    GeneticCode { id: 23, name: "Thraustochytrium Mitochondrial", aas: "FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", starts: "--------------------------------M--M---------------M------------" },
    GeneticCode { id: 24, name: "Rhabdopleuridae (Pterobranchia) Mitochondrial", aas: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG", starts: "---M---------------M---------------M---------------M------------" },
    GeneticCode { id: 25, name: "Candidate Division SR1 and Gracilibacteria", aas: "FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", starts: "---M-------------------------------M---------------M------------" },
    GeneticCode { id: 26, name: "Pachysolen tannophilus Nuclear", aas: "FFLLSSSSYY**CC*WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", starts: "-------------------M---------------M----------------------------" },
    GeneticCode { id: 29, name: "Mesodinium Nuclear", aas: "FFLLSSSSYYYYCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", starts: "-----------------------------------M----------------------------" },
    GeneticCode { id: 30, name: "Peritrich Nuclear", aas: "FFLLSSSSYYEECC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", starts: "-----------------------------------M----------------------------" },
    GeneticCode { id: 32, name: "Balanophoraceae Plastid", aas: "FFLLSSSSYY*WCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", starts: "---M---------------M------------MMMM---------------M------------" },
    GeneticCode { id: 33, name: "Cephalodiscidae Mitochondrial (UAA-Tyr)", aas: "FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG", starts: "---M---------------M---------------M---------------M------------" },
];

/// The standard code (table 1): the default everywhere.
pub fn standard() -> &'static GeneticCode {
    &TABLES[0]
}

/// The table with this NCBI id, if we have it.
pub fn by_id(id: u8) -> Option<&'static GeneticCode> {
    TABLES.iter().find(|t| t.id == id)
}

/// Like `by_id`, but an unknown id falls back to the standard code (used for the app's setting, which must never fail).
pub fn by_id_or_standard(id: u8) -> &'static GeneticCode {
    by_id(id).unwrap_or_else(standard)
}

/// Position of a codon in the 64-letter strings (TCAG order). `None` if it holds anything but A, C, G, T/U
/// (an N, say): such a codon has no amino acid. Lower case and U are accepted.
fn codon_index(codon: &[u8]) -> Option<usize> {
    if codon.len() != 3 {
        return None;
    }
    let code = |b: u8| match b.to_ascii_uppercase() {
        b'T' | b'U' => Some(0),
        b'C' => Some(1),
        b'A' => Some(2),
        b'G' => Some(3),
        _ => None,
    };
    Some(code(codon[0])? * 16 + code(codon[1])? * 4 + code(codon[2])?)
}

impl GeneticCode {
    /// Amino acid of a 3-base codon: a letter, '*' for a stop, 'X' when the codon is not plain ACGT.
    pub fn aa(&self, codon: &[u8]) -> u8 {
        match codon_index(codon) {
            Some(i) => self.aas.as_bytes()[i],
            None => b'X',
        }
    }

    pub fn is_stop(&self, codon: &[u8]) -> bool {
        self.aa(codon) == b'*'
    }

    /// Can this codon start a protein? With `alternative` false only ATG counts (the app's default, and what the ORF
    /// finder always did); with it true every start codon of the table counts (e.g. GTG and TTG in table 11).
    pub fn is_start(&self, codon: &[u8], alternative: bool) -> bool {
        if !alternative {
            return codon_index(codon).is_some() && codon_index(codon) == codon_index(b"ATG");
        }
        codon_index(codon).is_some_and(|i| self.starts.as_bytes()[i] == b'M')
    }

    /// Amino acids of `seq` read from its first base ('X' for codons with N or other letters).
    pub fn translate(&self, seq: &[u8]) -> Vec<u8> {
        seq.chunks_exact(3).map(|c| self.aa(c)).collect()
    }

    /// The DNA codons whose letter in `letters` (`aas` or `starts`) is `wanted`, in TCAG order (for display).
    pub fn codons_where(&self, letters: &str, wanted: u8) -> Vec<String> {
        const B: &[u8; 4] = b"TCAG";
        (0..64)
            .filter(|&i| letters.as_bytes()[i] == wanted)
            .map(|i| String::from_utf8(vec![B[i / 16], B[i / 4 % 4], B[i % 4]]).unwrap())
            .collect()
    }
}

// ---- translating an annotated CDS ----------------------------------------------------------------------------------------

/// Everything needed to translate one CDS feature the way its annotation says.
pub struct CdsPlan {
    /// NCBI table: the feature's own /transl_table, else the default passed in.
    pub table: &'static GeneticCode,
    /// Bases skipped before the first codon: /codon_start - 1 (0, 1 or 2).
    pub frame: usize,
    /// True when /codon_start was written in the feature.
    pub has_codon_start: bool,
    /// The 5' end (where translation starts) / 3' end is marked partial, "<" or ">" in the location.
    pub partial5: bool,
    pub partial3: bool,
    /// Genomic positions (0-based) of the CDS bases in reading order, before `frame` is skipped.
    pub positions: Vec<usize>,
    /// True when the CDS is on the bottom strand (bases read complemented, backwards).
    pub reverse: bool,
    /// /transl_except: (codon number counted from the first codon after `frame`, amino acid letter or '*' for TERM).
    pub exceptions: Vec<(usize, u8)>,
}

/// GenBank amino acid abbreviations used in /transl_except (aa:Sec). "TERM" is a stop codon.
fn aa_letter(name: &str) -> Option<u8> {
    Some(match name {
        "Ala" => b'A', "Arg" => b'R', "Asn" => b'N', "Asp" => b'D', "Cys" => b'C', "Gln" => b'Q', "Glu" => b'E',
        "Gly" => b'G', "His" => b'H', "Ile" => b'I', "Leu" => b'L', "Lys" => b'K', "Met" => b'M', "Phe" => b'F',
        "Pro" => b'P', "Ser" => b'S', "Thr" => b'T', "Trp" => b'W', "Tyr" => b'Y', "Val" => b'V',
        "Sec" => b'U', "Pyl" => b'O', "TERM" => b'*',
        _ => return None,
    })
}

/// A plain CDS location: its ranges in the order written, whether it is complemented, and the partial flags
/// (5' and 3' in reading direction).
type PlainLocation = (Vec<(i64, i64)>, bool, bool, bool);

/// Read a CDS location of the usual kinds: `a..b`, `join(...)`, `order(...)` and `complement` of those. `None` for
/// what we do not handle: complement inside a join, "between" positions, external references, one-of positions.
fn plain_location(loc: &Location) -> Option<PlainLocation> {
    // (start, start is "<", end, end is ">") of a Range; None for anything else.
    fn range(l: &Location) -> Option<(i64, bool, i64, bool)> {
        match l {
            Location::Range((s, before), (e, after)) => Some((*s, before.0, *e, after.0)),
            _ => None,
        }
    }
    // The ranges of a Range, or of a Join/Order made only of Ranges.
    fn parts(l: &Location) -> Option<Vec<(i64, bool, i64, bool)>> {
        match l {
            Location::Join(v) | Location::Order(v) => v.iter().map(range).collect(),
            other => range(other).map(|r| vec![r]),
        }
    }
    let (reverse, inner) = match loc {
        Location::Complement(inner) => (true, inner.as_ref()),
        other => (false, other),
    };
    let ps = parts(inner)?;
    let first = ps.first()?;
    let last = ps.last()?;
    // On the top strand reading starts at the first range's start; on the bottom strand at the last range's end.
    let (partial5, partial3) = if reverse { (last.3, first.1) } else { (first.1, last.3) };
    Some((ps.iter().map(|p| (p.0, p.2)).collect(), reverse, partial5, partial3))
}

/// Parse one /transl_except value, "(pos:1000..1002,aa:Sec)" or "(pos:complement(1000..1002),aa:TERM)", into
/// (0-based first base of the codon in reading order, complement?, amino acid letter).
/// Only a full 3-base codon is accepted: a stop completed by poly-A (pos spanning 1 or 2 bases) cannot be
/// reproduced from the bases, so it is an error, and so is any form we do not understand.
fn parse_exception(text: &str) -> Result<(i64, bool, u8), String> {
    let bad = || format!("/transl_except \"{text}\" not understood");
    let t = text.trim().trim_start_matches('(').trim_end_matches(')');
    let (pos, aa) = t.split_once(",aa:").ok_or_else(bad)?;
    let pos = pos.strip_prefix("pos:").ok_or_else(bad)?;
    let (complement, pos) = match pos.strip_prefix("complement(") {
        Some(rest) => (true, rest.strip_suffix(')').ok_or_else(bad)?),
        None => (false, pos),
    };
    let (a, b) = pos.split_once("..").ok_or_else(bad)?;
    let a: i64 = a.parse().map_err(|_| bad())?;
    let b: i64 = b.parse().map_err(|_| bad())?;
    if b - a != 2 {
        return Err(format!("/transl_except \"{text}\" is not a whole codon"));
    }
    let aa = aa_letter(aa.trim().trim_end_matches(')')).ok_or_else(bad)?;
    // The codon's first base in reading order: the lowest base on the top strand, the highest on the bottom strand.
    Ok((if complement { b - 1 } else { a - 1 }, complement, aa))
}

/// Work out how to translate CDS `f` of a record of `n` bases. `default_table` is used when the feature has no
/// /transl_table. `Err` explains why the feature cannot be translated from its own bases (bad table number,
/// bad /codon_start, a location we do not handle, a /transl_except we cannot reproduce).
pub fn plan_cds(f: &Feature, n: usize, default_table: u8) -> Result<CdsPlan, String> {
    // /transl_table=N overrides the setting for this feature.
    let table_id = match f.qualifier_values("transl_table").next() {
        Some(v) => v.trim().parse::<u8>().map_err(|_| format!("bad /transl_table \"{v}\""))?,
        None => default_table,
    };
    let table = by_id(table_id).ok_or_else(|| format!("genetic code table {table_id} is not available"))?;
    let (has_codon_start, frame) = match f.qualifier_values("codon_start").next() {
        Some(v) => match v.trim() {
            "1" => (true, 0),
            "2" => (true, 1),
            "3" => (true, 2),
            other => return Err(format!("bad /codon_start \"{other}\"")),
        },
        None => (false, 0),
    };
    let (ranges, reverse, partial5, partial3) = plain_location(&f.location).ok_or("unusual location")?;
    if ranges.iter().any(|&(s, e)| s < 0 || e < s || e as usize > n) {
        return Err("location outside the sequence".into());
    }
    // Base positions in reading order: the ranges as written, backwards as a whole on the bottom strand.
    let mut positions: Vec<usize> = ranges.iter().flat_map(|&(s, e)| s as usize..e as usize).collect();
    if reverse {
        positions.reverse();
    }
    let mut exceptions = Vec::new();
    for text in f.qualifier_values("transl_except") {
        let (first, complement, aa) = parse_exception(text)?;
        if complement != reverse {
            return Err(format!("/transl_except \"{text}\" is on the other strand"));
        }
        // The codon must sit in the reading frame, on three consecutive bases of the CDS.
        let offset = positions.iter().position(|&p| p as i64 == first).ok_or("/transl_except outside the CDS")?;
        if offset < frame || (offset - frame) % 3 != 0 || offset + 3 > positions.len() {
            return Err(format!("/transl_except \"{text}\" is not on a codon of the reading frame"));
        }
        exceptions.push(((offset - frame) / 3, aa));
    }
    Ok(CdsPlan { table, frame, has_codon_start, partial5, partial3, positions, reverse, exceptions })
}

/// Why a CDS keeps its old /translation (see `cds_protein`).
#[derive(Debug, PartialEq, Clone, Copy)]
pub enum Skip {
    Pseudo,
    Partial,
    Other,
}

/// The protein of CDS `f`, recomputed from its bases the way GenBank writes /translation: its own /transl_table and
/// /codon_start, a start codon of the table read as methionine, no final stop. `Err((kind, why))` when it should be
/// left alone: /pseudo or /pseudogene, a 5' end marked partial with no /codon_start (the frame is unknown), a
/// /transl_except we cannot reproduce, a stop codon inside the protein, a length that is not a multiple of 3, a CDS
/// that is not even one whole codon.
pub fn cds_protein(seq: &[u8], f: &Feature, default_table: u8) -> Result<String, (Skip, String)> {
    if f.qualifiers.iter().any(|(k, _)| k == "pseudo" || k == "pseudogene") {
        return Err((Skip::Pseudo, "pseudogene".into()));
    }
    let plan = plan_cds(f, seq.len(), default_table).map_err(|e| (Skip::Other, e))?;
    if plan.partial5 && !plan.has_codon_start {
        return Err((Skip::Partial, "partial 5' end without /codon_start: the reading frame is undefined".into()));
    }
    // The bases in reading order (complemented on the bottom strand), minus the first `frame` bases.
    let mut bases: Vec<u8> = plan.positions.iter().map(|&p| seq[p].to_ascii_uppercase()).collect();
    if plan.reverse {
        // `positions` is already in reading order (reversed), so each base only needs complementing:
        // the reverse complement of a single base is its complement.
        bases = bases.iter().map(|&b| revcomp(&[b])[0]).collect();
    }
    if bases.len() < plan.frame + 3 {
        return Err((Skip::Other, "shorter than one codon".into()));
    }
    let bases = &bases[plan.frame..];
    if bases.len() % 3 != 0 && !plan.partial3 {
        return Err((Skip::Other, "length is not a multiple of 3".into()));
    }
    let mut protein = plan.table.translate(bases);
    // A /transl_except codon replaces whatever the table says (selenocysteine, pyrrolysine, a stop).
    for &(codon, aa) in &plan.exceptions {
        if codon < protein.len() {
            protein[codon] = aa;
        }
    }
    // GenBank translates a start codon of the table as methionine, but only when the CDS really begins there.
    let exception_at_start = plan.exceptions.iter().any(|&(c, _)| c == 0);
    if !plan.partial5 && plan.frame == 0 && !exception_at_start && plan.table.is_start(&bases[..3], true) {
        protein[0] = b'M';
    }
    if protein.last() == Some(&b'*') {
        protein.pop(); // the final stop is not part of /translation
    }
    if protein.contains(&b'*') {
        return Err((Skip::Other, "stop codon inside the reading frame".into()));
    }
    Ok(String::from_utf8(protein).unwrap_or_default())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A CDS feature from GenBank location text and (key, value) qualifiers.
    fn cds(location: &str, quals: &[(&str, &str)]) -> Feature {
        Feature {
            kind: "CDS".into(),
            location: Location::from_gb_format(location).unwrap(),
            qualifiers: quals.iter().map(|(k, v)| ((*k).to_string().into(), Some((*v).to_string()))).collect(),
        }
    }

    #[test]
    fn every_table_has_64_codons_and_known_ids() {
        for t in TABLES {
            assert_eq!(t.aas.len(), 64, "table {}", t.id);
            assert_eq!(t.starts.len(), 64, "table {}", t.id);
        }
        assert!(by_id(1).is_some() && by_id(11).is_some());
        assert!(by_id(27).is_none(), "context-dependent tables are left out");
    }

    #[test]
    fn the_standard_code_and_the_known_differences() {
        let std = standard();
        assert_eq!(std.translate(b"ATGGCCTAAGNA"), b"MA*X");
        // Table 4 (mycoplasma): TGA is tryptophan, not a stop. Table 2: AGA and AGG are stops.
        assert_eq!(by_id(4).unwrap().aa(b"TGA"), b'W');
        assert_eq!(std.aa(b"TGA"), b'*');
        assert_eq!(by_id(2).unwrap().aa(b"AGA"), b'*');
        assert_eq!(std.aa(b"AGA"), b'R');
        // Lower case and RNA letters are fine.
        assert_eq!(std.aa(b"aug"), b'M');
    }

    #[test]
    fn a_cds_is_translated_with_its_own_table_and_frame() {
        let seq = b"ATGTGAGTGTAA"; // ATG TGA GTG TAA
        let plain = cds("1..12", &[]);
        // Standard code: TGA at codon 2 is a stop inside the CDS, so it is left alone.
        assert!(cds_protein(seq, &plain, 1).is_err());
        // The same bases with /transl_table=4: M W V (GTG is a start codon, but not the first codon, so V).
        let t4 = cds("1..12", &[("transl_table", "4")]);
        assert_eq!(cds_protein(seq, &t4, 1).unwrap(), "MWV");
        // The setting's table is used when the feature has none; the feature wins when it has one.
        assert_eq!(cds_protein(seq, &plain, 4).unwrap(), "MWV");
        assert_eq!(cds_protein(seq, &t4, 11).unwrap(), "MWV");
        // /codon_start=2 skips the A: TGT GAG TGT AA. The 2 last bases are an incomplete codon: fine only if partial.
        let cs2 = cds("1..12", &[("codon_start", "2")]);
        assert!(cds_protein(seq, &cs2, 1).is_err());
        let cs2p = cds("1..>12", &[("codon_start", "2")]);
        assert_eq!(cds_protein(seq, &cs2p, 1).unwrap(), "CEC");
    }

    #[test]
    fn an_alternative_start_codon_is_read_as_methionine() {
        // TTG is a start codon of table 1 (and 11): GenBank writes M for it at the start of a CDS.
        let seq = b"TTGAAATAA";
        assert_eq!(cds_protein(seq, &cds("1..9", &[]), 1).unwrap(), "MK");
        // Not at the start: a plain leucine. And not when the 5' end is partial.
        let inner = b"ATGTTGTAA";
        assert_eq!(cds_protein(inner, &cds("1..9", &[]), 1).unwrap(), "ML");
        assert_eq!(cds_protein(seq, &cds("<1..9", &[("codon_start", "1")]), 1).unwrap(), "LK");
    }

    #[test]
    fn bottom_strand_cds_and_transl_except() {
        // CDS on the bottom strand: ATG TGA TAA, with the TGA (bases 4..6 of the CDS) read as selenocysteine.
        let top = revcomp(b"ATGTGATAA"); // 9 bases; the CDS is complement(1..9)
        let f = cds("complement(1..9)", &[("transl_except", "(pos:complement(4..6),aa:Sec)")]);
        assert_eq!(cds_protein(&top, &f, 1).unwrap(), "MU");
        // The same without the exception has a stop inside: left alone.
        assert!(cds_protein(&top, &cds("complement(1..9)", &[]), 1).is_err());
        // A form we cannot reproduce is an error, never translated silently.
        let odd = cds("complement(1..9)", &[("transl_except", "(pos:complement(4..5),aa:TERM)")]);
        assert!(cds_protein(&top, &odd, 1).is_err());
    }

    #[test]
    fn joined_cds_reads_its_parts_in_order() {
        // join(1..3,7..12): ATG + AAATAA -> M K stop.
        let seq = b"ATGCCCAAATAA";
        assert_eq!(cds_protein(seq, &cds("join(1..3,7..12)", &[]), 1).unwrap(), "MK");
    }

    #[test]
    fn pseudo_and_partial_cds_are_skipped() {
        let seq = b"ATGAAATAA";
        assert_eq!(cds_protein(seq, &cds("1..9", &[("pseudo", "")]), 1).map_err(|e| e.0), Err(Skip::Pseudo));
        assert_eq!(cds_protein(seq, &cds("<1..9", &[]), 1).map_err(|e| e.0), Err(Skip::Partial));
        // ... unless /codon_start says where the frame is.
        assert_eq!(cds_protein(seq, &cds("<1..9", &[("codon_start", "1")]), 1).unwrap(), "MK");
    }
}
