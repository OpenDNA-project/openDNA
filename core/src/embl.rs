//! EMBL flat files, read by turning each record into GenBank text that gb-io parses. The feature table
//! (FT lines) has the same columns as GenBank's, so only the header and the sequence block are rewritten.
//! Reading only: a document is always saved as GenBank.

/// True for text that starts like an EMBL record ("ID   ...").
pub fn looks_like_embl(head: &[u8]) -> bool {
    head.starts_with(b"ID   ")
}

/// The GenBank text of an EMBL file (one GenBank record per EMBL record).
pub fn to_genbank(text: &str) -> String {
    let mut out = String::new();
    let (mut locus, mut definition, mut features, mut bases) = (String::new(), String::new(), String::new(), String::new());
    let mut in_sequence = false;
    for line in text.lines() {
        if in_sequence && !line.starts_with("//") {
            bases.extend(line.chars().filter(|c| c.is_ascii_alphabetic())); // drops the spaces and the position numbers
        } else if let Some(id) = line.strip_prefix("ID   ") {
            locus = locus_line(id);
        } else if let Some(de) = line.strip_prefix("DE   ") {
            if !definition.is_empty() {
                definition.push(' ');
            }
            definition.push_str(de.trim());
        } else if let Some(ft) = line.strip_prefix("FT") {
            features.push_str("  ");
            features.push_str(ft);
            features.push('\n');
        } else if line.starts_with("SQ   ") {
            in_sequence = true;
        } else if line.starts_with("//") {
            out.push_str(&format!("{locus}\nDEFINITION  {}\nFEATURES             Location/Qualifiers\n{features}ORIGIN\n", if definition.is_empty() { "." } else { &definition }));
            for (i, row) in bases.as_bytes().chunks(60).enumerate() {
                out.push_str(&format!("{:>9}", i * 60 + 1));
                row.chunks(10).for_each(|g| out.push_str(&format!(" {}", String::from_utf8_lossy(g))));
                out.push('\n');
            }
            out.push_str("//\n");
            (definition, features, bases, in_sequence) = (String::new(), String::new(), String::new(), false);
        }
    }
    out
}

/// "HQ324066; SV 1; linear; genomic DNA; STD; PLN; 84 BP." (or the older "X  standard; circular DNA; CON; 84 BP.") as a LOCUS line.
fn locus_line(id: &str) -> String {
    let name = id.split(|c: char| c == ';' || c.is_whitespace()).next().unwrap_or("");
    let lower = id.to_lowercase();
    let topology = if lower.contains("circular") { "circular" } else { "linear" };
    let molecule = if lower.contains("rna") { "RNA" } else { "DNA" };
    let len = id.split_whitespace().rev().skip_while(|w| !w.trim_end_matches('.').eq_ignore_ascii_case("bp")).nth(1).unwrap_or("0");
    format!("LOCUS       {name}                {len} bp    {molecule}     {topology}   UNK 01-JAN-1980")
}
