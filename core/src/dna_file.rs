//! Read-only `.dna` reader. Output is a gb-io `Seq`, i.e. the same in-memory
//! record a GenBank file produces, so the rest of the app never knows the difference.
//!
//! File layout: a stream of packets `[type: 1 byte][length: 4 bytes big-endian][payload]`.
//! Packets we use:
//!   0x09  cookie, payload starts with "SnapGene"   (file signature)
//!   0x00  DNA: 1 flag byte (bit 0 = circular) followed by the bases
//!   0x0A  features, XML
//!   0x05  primers, XML
//!   0x06  notes, XML (we take <Description> as the GenBank DEFINITION)
//! Everything else (undo history, cached enzyme sites, alignments, display settings)
//! is skipped on purpose.
//!
//! Primers become `primer_bind` features, exactly like a GenBank export of the same file,
//! so they survive "save as GenBank".

use std::borrow::Cow;
use std::io;
use std::path::Path;

use gb_io::seq::{Feature, Location, Seq, Topology};
use quick_xml::events::{BytesStart, Event};
use quick_xml::{Reader, XmlVersion};

/// A primer as stored in the file (name + full sequence including any 5' tail).
#[derive(Debug, Clone, PartialEq)]
pub struct Primer {
    pub name: String,
    pub sequence: String,
}

/// Parse a `.dna` file from disk.
pub fn read_file(path: &Path) -> io::Result<(Seq, Vec<Primer>)> {
    let data = std::fs::read(path)?;
    // `file_stem` = file name without extension; it becomes the record name.
    let name = path.file_stem().map(|s| s.to_string_lossy().into_owned());
    let (mut seq, primers) = parse(&data)?;
    seq.name = name;
    Ok((seq, primers))
}

/// Parse `.dna` bytes already in memory.
pub fn parse(data: &[u8]) -> io::Result<(Seq, Vec<Primer>)> {
    let mut seq = Seq::empty();
    let mut primers = Vec::new();
    let mut have_dna = false;
    // Features/primers XML need the sequence length (for origin-spanning ranges), and the
    // DNA packet comes first in every file seen so far; still, collect XML and parse at the end
    // so packet order never matters.
    let mut features_xml: Option<&[u8]> = None;
    let mut primers_xml: Option<&[u8]> = None;

    let mut pos = 0usize;
    let mut first = true;
    while pos < data.len() {
        // `.get(a..b)` returns None instead of crashing if the file is truncated.
        let header = data.get(pos..pos + 5).ok_or_else(|| bad("truncated packet header"))?;
        let kind = header[0];
        let len = u32::from_be_bytes([header[1], header[2], header[3], header[4]]) as usize;
        let body = data.get(pos + 5..pos + 5 + len).ok_or_else(|| bad("truncated packet"))?;
        pos += 5 + len;

        if first {
            if kind != 0x09 || !body.starts_with(b"SnapGene") {
                return Err(bad("not a .dna file"));
            }
            first = false;
            continue;
        }
        match kind {
            0x00 => {
                // `split_first` -> (first byte, rest). An empty packet is malformed.
                let (flags, bases) = body.split_first().ok_or_else(|| bad("empty DNA packet"))?;
                seq.topology = if flags & 1 == 1 { Topology::Circular } else { Topology::Linear };
                seq.seq = bases.to_vec();
                have_dna = true;
            }
            0x0A => features_xml = Some(body),
            0x05 => primers_xml = Some(body),
            0x06 => seq.definition = xml_text(body, "Description"),
            _ => {} // history, caches, display settings: not needed
        }
    }
    if !have_dna {
        return Err(bad("no DNA sequence in file (protein or empty .dna file?)"));
    }
    seq.molecule_type = Some("DNA".into());
    seq.len = Some(seq.seq.len());

    let seq_len = seq.seq.len() as i64;
    if let Some(xml) = features_xml {
        seq.features = parse_features(xml, seq_len)?;
    }
    if let Some(xml) = primers_xml {
        let (list, bind_features) = parse_primers(xml, seq_len)?;
        primers = list;
        seq.features.extend(bind_features);
    }
    Ok((seq, primers))
}

/// Turn the features XML into gb-io features.
///
/// Shape: <Feature name=".." type=".." directionality="1|2">   (absent or 0: no direction, 3: both)
///          <Segment range="10-200" color="#ff0000" type="standard"/>   (1-based, inclusive)
///          <Q name="note"><V text="..."/></Q>
///        </Feature>
fn parse_features(xml: &[u8], seq_len: i64) -> io::Result<Vec<Feature>> {
    let mut reader = Reader::from_reader(xml);
    let mut features = Vec::new();

    // State of the feature currently being read (reset at each <Feature>).
    let mut kind = String::new();
    let mut name = String::new();
    let mut reverse = false;
    // The `directionality` attribute as written (None when absent). GenBank has no way to say "no direction" in a
    // location, so a feature drawn without an arrow gets a /directional="no" qualifier (see DIRECTIONLESS_KINDS in
    // app/editor.py: misc_feature has no arrow unless /directional="yes").
    let mut direction: Option<String> = None;
    let mut segments: Vec<Location> = Vec::new();
    let mut quals: Vec<(Cow<'static, str>, Option<String>)> = Vec::new();
    let mut color: Option<String> = None;
    let mut q_name: Option<String> = None;

    loop {
        // Start = <Tag ...>, Empty = <Tag .../>, End = </Tag>. Both Start and Empty carry
        // attributes, so we handle them together.
        match reader.read_event().map_err(xml_err)? {
            Event::Start(e) | Event::Empty(e) => match e.name().as_ref() {
                "Feature" => {
                    kind = attr(&e, "type")?.unwrap_or_else(|| "misc_feature".into());
                    name = attr(&e, "name")?.unwrap_or_default();
                    direction = attr(&e, "directionality")?;
                    // `as_deref()` turns Option<String> into Option<&str>, so it can be compared with a literal.
                    reverse = direction.as_deref() == Some("2");
                    segments.clear();
                    quals.clear();
                    color = None;
                }
                "Segment" => {
                    // "gap" segments are spacers between parts, not part of the feature.
                    if attr(&e, "type")?.as_deref() == Some("gap") {
                        continue;
                    }
                    if let Some(range) = attr(&e, "range")? {
                        segments.extend(parse_range(&range, seq_len, 1)?);
                    }
                    if color.is_none() {
                        color = attr(&e, "color")?;
                    }
                }
                "Q" => q_name = attr(&e, "name")?,
                "V" => {
                    // A value is stored as text="..", int=".." or predef=".." (a preset word).
                    let value = match attr(&e, "text")? {
                        Some(v) => Some(v),
                        None => match attr(&e, "int")? {
                            Some(v) => Some(v),
                            None => attr(&e, "predef")?,
                        },
                    };
                    if let Some(q) = &q_name {
                        quals.push((Cow::Owned(q.clone()), value));
                    }
                }
                _ => {}
            },
            Event::End(e) if e.name().as_ref() == "Feature" => {
                if segments.is_empty() {
                    continue; // malformed feature without a position: nothing to keep
                }
                // `std::mem::take` moves the Vec out and leaves an empty one behind,
                // ready for the next feature (no copy).
                let mut loc = if segments.len() == 1 {
                    segments.pop().unwrap()
                } else {
                    Location::Join(std::mem::take(&mut segments))
                };
                if reverse {
                    loc = Location::Complement(Box::new(loc));
                }
                let mut qualifiers = std::mem::take(&mut quals);
                // The file keeps the display name as an attribute; GenBank uses /label.
                if !name.is_empty() && !qualifiers.iter().any(|(k, _)| k == "label") {
                    qualifiers.insert(0, (Cow::Borrowed("label"), Some(name.clone())));
                }
                // Direction: the file's flag against this program's default for the type (misc_feature: none).
                // `matches!` is a yes/no test of a value against patterns ("1" | "2" | "3" = any of them).
                let directional = matches!(direction.as_deref(), Some("1") | Some("2") | Some("3"));
                if directional != (kind != "misc_feature") {
                    let value = if directional { "yes" } else { "no" };
                    qualifiers.push((Cow::Borrowed("directional"), Some(value.to_string())));
                }
                // Colors as ApE/Benchling-style qualifiers, which other tools understand.
                if let Some(c) = color.take() {
                    qualifiers.push((Cow::Borrowed("ApEinfo_fwdcolor"), Some(c.clone())));
                    qualifiers.push((Cow::Borrowed("ApEinfo_revcolor"), Some(c)));
                }
                features.push(Feature { kind: Cow::Owned(kind.clone()), location: loc, qualifiers });
            }
            Event::Eof => break,
            _ => {}
        }
    }
    Ok(features)
}

/// Primers XML: <Primer name=".." sequence=".."><BindingSite location="a-b" boundStrand="0|1">
/// Returns the primer list plus one `primer_bind` feature per binding site.
fn parse_primers(xml: &[u8], seq_len: i64) -> io::Result<(Vec<Primer>, Vec<Feature>)> {
    let mut reader = Reader::from_reader(xml);
    let mut primers: Vec<Primer> = Vec::new();
    let mut features = Vec::new();
    loop {
        match reader.read_event().map_err(xml_err)? {
            Event::Start(e) | Event::Empty(e) => match e.name().as_ref() {
                "Primer" => primers.push(Primer {
                    name: attr(&e, "name")?.unwrap_or_default(),
                    sequence: attr(&e, "sequence")?.unwrap_or_default(),
                }),
                // The file lists each site twice; the copy marked simplified="1" is a duplicate.
                "BindingSite" if attr(&e, "simplified")?.is_none() => {
                    let Some(p) = primers.last() else { continue };
                    let Some(range) = attr(&e, "location")? else { continue };
                    let mut parts = parse_range(&range, seq_len, 0)?;
                    let mut loc = if parts.len() == 1 { parts.pop().unwrap() } else { Location::Join(parts) };
                    if attr(&e, "boundStrand")?.as_deref() == Some("1") {
                        loc = Location::Complement(Box::new(loc));
                    }
                    features.push(Feature {
                        kind: Cow::Borrowed("primer_bind"),
                        location: loc,
                        qualifiers: vec![
                            (Cow::Borrowed("label"), Some(p.name.clone())),
                            (Cow::Borrowed("note"), Some(format!("sequence: {}", p.sequence))),
                        ],
                    });
                }
                _ => {}
            },
            Event::Eof => break,
            _ => {}
        }
    }
    Ok((primers, features))
}

/// "10-200" (inclusive) -> [Range(9, 200)]. On a circular sequence the file writes
/// origin-spanning ranges as "4000-120": that becomes two ranges, [3999, len) and [0, 120),
/// the same `join(4000..len,1..120)` GenBank uses.
/// `first` is the number of the first base: 1 for feature segments, but 0 for primer
/// binding sites (checked against the bases in real files: real files mix both).
fn parse_range(text: &str, seq_len: i64, first: i64) -> io::Result<Vec<Location>> {
    let (a, b) = text.split_once('-').ok_or_else(|| bad(&format!("bad range '{text}'")))?;
    let num = |s: &str| s.trim().parse::<i64>().map_err(|_| bad(&format!("bad range '{text}'")));
    // Convert to 1-based numbering, then the same rules for both.
    let (a, b) = (num(a)? + 1 - first, num(b)? + 1 - first);
    if a < 1 || b < 1 || a > seq_len || b > seq_len {
        return Err(bad(&format!("range '{text}' outside sequence of {seq_len} bp")));
    }
    Ok(if a <= b {
        vec![Location::simple_range(a - 1, b)]
    } else {
        vec![Location::simple_range(a - 1, seq_len), Location::simple_range(0, b)]
    })
}

/// Value of attribute `key` on an XML tag, with XML escapes (&lt; etc.) decoded.
fn attr(e: &BytesStart, key: &str) -> io::Result<Option<String>> {
    for a in e.attributes() {
        let a = a.map_err(xml_err)?;
        if a.key.as_ref() == key {
            return Ok(Some(a.normalized_value(XmlVersion::Implicit1_0).map_err(xml_err)?.into_owned()));
        }
    }
    Ok(None)
}

/// Text of the first `<tag>...</tag>` in a small XML blob (used for the Notes packet).
fn xml_text(xml: &[u8], tag: &str) -> Option<String> {
    let s = String::from_utf8_lossy(xml);
    let open = format!("<{tag}>");
    let start = s.find(&open)? + open.len();
    let end = start + s[start..].find(&format!("</{tag}>"))?;
    let text = s[start..end].trim();
    if text.is_empty() { None } else { Some(text.to_string()) }
}

// `impl std::fmt::Display` = "anything printable"; lets one helper accept every XML error type.
fn xml_err(e: impl std::fmt::Display) -> io::Error {
    bad(&format!(".dna XML: {e}"))
}

fn bad(msg: &str) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, msg.to_string())
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;

    /// Build a packet: type byte + big-endian length + payload.
    fn packet(kind: u8, body: &[u8]) -> Vec<u8> {
        let mut p = vec![kind];
        p.extend((body.len() as u32).to_be_bytes());
        p.extend(body);
        p
    }

    /// A tiny synthetic circular .dna: 20 bp, one forward feature, one origin-spanning
    /// reverse feature, a feature without direction and a misc_feature with one, one primer
    /// on the reverse strand, plus a junk packet to skip.
    pub(crate) fn sample_dna() -> Vec<u8> {
        let mut f = packet(0x09, b"SnapGene\x00\x01\x00\x0f\x00\x14");
        f.extend(packet(0x00, b"\x01AAAACCCCGGGGTTTTACGT"));
        f.extend(packet(0x07, b"history blob, ignored"));
        f.extend(packet(0x0A, br##"<?xml version="1.0"?><Features>
            <Feature name="lacZ" type="CDS" directionality="1"><Segment range="2-10" color="#ff0000" type="standard"/>
              <Q name="note"><V text="a &lt;b&gt; c"/></Q><Q name="codon_start"><V int="1"/></Q></Feature>
            <Feature name="ori" type="rep_origin" directionality="2"><Segment range="18-3" color="#ffff00" type="standard"/></Feature>
            <Feature name="box" type="protein_bind"><Segment range="5-8" type="standard"/></Feature>
            <Feature name="tag" type="misc_feature" directionality="1"><Segment range="11-14" type="standard"/></Feature>
            </Features>"##));
        f.extend(packet(0x05, br#"<?xml version="1.0"?><Primers><Primer name="p1" sequence="ttttGGGG">
            <BindingSite location="9-12" boundStrand="1"/><BindingSite simplified="1" location="9-12" boundStrand="1"/>
            </Primer></Primers>"#));
        f.extend(packet(0x06, b"<Notes><Description>test plasmid</Description></Notes>"));
        f
    }

    #[test]
    fn reads_synthetic_file() {
        let (seq, primers) = parse(&sample_dna()).unwrap();
        assert_eq!(seq.seq, b"AAAACCCCGGGGTTTTACGT");
        assert_eq!(seq.topology, Topology::Circular);
        assert_eq!(seq.definition.as_deref(), Some("test plasmid"));
        assert_eq!(seq.features.len(), 5); // 4 features + 1 primer_bind (duplicate site skipped)

        let lacz = &seq.features[0];
        assert_eq!(lacz.location, Location::simple_range(1, 10));
        assert_eq!(lacz.qualifier_values("label").next(), Some("lacZ"));
        assert_eq!(lacz.qualifier_values("note").next(), Some("a <b> c")); // XML unescaped
        assert_eq!(lacz.qualifier_values("codon_start").next(), Some("1"));

        // "18-3" on a 20 bp circle = complement(join(18..20,1..3))
        let ori = &seq.features[1].location;
        assert_eq!(ori.to_gb_format(), "complement(join(18..20,1..3))");

        // Direction: only what differs from the type's default is written (lacZ: a CDS with an arrow, nothing added).
        assert_eq!(lacz.qualifier_values("directional").next(), None);
        assert_eq!(seq.features[2].qualifier_values("directional").next(), Some("no")); // no directionality attribute
        assert_eq!(seq.features[3].qualifier_values("directional").next(), Some("yes")); // a misc_feature with an arrow

        assert_eq!(primers, vec![Primer { name: "p1".into(), sequence: "ttttGGGG".into() }]);
        assert_eq!(seq.features[4].kind, "primer_bind");
        // Binding sites are 0-based: "9-12" = bases 10..13 in GenBank numbering.
        assert_eq!(seq.features[4].location.to_gb_format(), "complement(10..13)");
    }

    #[test]
    fn rejects_other_files() {
        // A GenBank text: its "length" bytes are huge, so it fails as a truncated packet.
        assert!(parse(b"LOCUS       x").unwrap_err().to_string().contains("truncated"));
        let mut no_dna = packet(0x09, b"SnapGene\x00");
        no_dna.extend(packet(0x06, b"<Notes/>"));
        assert!(parse(&no_dna).is_err());
        // --- well-formed packets with a DNA packet but the wrong cookie: only the signature check rejects it ---
        let mut f = packet(0x09, b"NotSnapGene");
        f.extend(packet(0x00, b"\x00ACGT"));
        let e = parse(&f).unwrap_err();
        assert!(e.to_string().contains("not a .dna file"), "{e}");
    }
}
