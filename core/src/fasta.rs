//! Random access to FASTA files through a `.fai` index (same format as `samtools faidx`).
//!
//! Why: a human genome FASTA is ~3 GB. Instead of loading it, we record for each sequence
//! where its bases start in the file and how lines are wrapped. Then fetching
//! `chr1:1,000,000-1,000,100` is one `seek` + one small read, whatever the file size.
//!
//! `.fai` line format (tab-separated):  name  length  offset  line_bases  line_width
//!   offset     = byte position of the first base
//!   line_bases = bases per full line (e.g. 60)
//!   line_width = bytes per full line including the newline (61 for "\n", 62 for "\r\n")

// `use` = Python's `import`. `std` is Rust's standard library.
use std::ffi::OsString;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};

/// One line of a `.fai` file. (`///` comments are documentation, like Python docstrings.)
///
/// `#[derive(...)]` auto-generates boilerplate: `Debug` lets us print it with `{:?}`,
/// `Clone` adds a `.clone()` method, `PartialEq` allows `==` (used in tests).
#[derive(Debug, Clone, PartialEq)]
pub struct FaiEntry {
    // `pub` = visible outside this file. Rust fields are private by default.
    pub name: String,
    pub length: u64, // u64 = unsigned 64-bit integer; genomes exceed the 4 GB limit of u32
    pub offset: u64,
    pub line_bases: u64,
    pub line_width: u64,
}

/// An open FASTA file plus its index. This is what Python sees as `opendna_core.Fasta`.
pub struct IndexedFasta {
    file: File, // private: only code in this file touches the file handle
    pub entries: Vec<FaiEntry>, // Vec<T> = growable list of T, like a Python list with one type
}

// `impl` block = the methods of a struct (like the body of a Python class).
impl IndexedFasta {
    /// Open `path`, using `path.fai` if it exists, otherwise building (and saving) it.
    ///
    /// Return type `io::Result<Self>` means "either Ok(an IndexedFasta) or Err(an I/O error)".
    /// Rust has no exceptions: errors are ordinary return values that the caller must handle.
    pub fn open(path: &Path) -> io::Result<Self> {
        let fai = fai_path_for(path);
        // `if` is an expression in Rust: its value is assigned to `entries`.
        // An index that is not NEWER than the FASTA may describe another file version: rebuild it. (Strictly newer, not
        // "same age or newer": a FASTA rewritten right after its index was saved can carry the very same timestamp, because
        // file times are coarse, and then the old index - with the old record lengths - would be trusted and reads would
        // run past the end of the shorter file. Rebuilding costs one scan, and the saved index is newer from then on.)
        let fresh = |f: &Path| -> Option<bool> {
            Some(std::fs::metadata(f).ok()?.modified().ok()? > std::fs::metadata(path).ok()?.modified().ok()?)
        };
        let entries = if fai.exists() && fresh(&fai) == Some(true) {
            read_fai(&fai)? // `?` = "if this is an Err, return it from open() right now"
        } else {
            let built = build_fai(path)?;
            // Saving the index is only a speed-up for next time; ignore failure
            // (e.g. the FASTA sits in a read-only folder). `let _ =` discards the result.
            let _ = write_fai(&fai, &built);
            built
        };
        // Last expression without `;` is the return value (no `return` keyword needed).
        Ok(IndexedFasta { file: File::open(path)?, entries })
    }

    /// Bases of sequence `name` in the 0-based half-open range [start, end), like Python slicing.
    /// `end` is clamped to the sequence length.
    ///
    /// `&mut self` = this method may modify the struct (seeking moves the file cursor).
    /// `&str` = a borrowed, read-only string (we only look at the name, never keep it).
    pub fn fetch(&mut self, name: &str, start: u64, end: u64) -> io::Result<Vec<u8>> {
        // `.iter().find(...)` scans the list; `|e| ...` is a closure (Python: lambda e: ...).
        // `find` returns an Option: Some(entry) or None. `ok_or_else` turns None into an error.
        let e = self
            .entries
            .iter()
            .find(|e| e.name == name)
            .ok_or_else(|| {
                io::Error::new(io::ErrorKind::NotFound, format!("no sequence named '{name}'"))
            })?;

        let end = end.min(e.length); // shadowing: a new `end` replaces the old one
        if start >= end {
            return Ok(Vec::new());
        }
        // Byte position of base number `b`: skip whole lines, then step into the current line.
        let byte_of = |b: u64| e.offset + (b / e.line_bases) * e.line_width + b % e.line_bases;
        let first = byte_of(start);
        let last = byte_of(end - 1) + 1;

        // `vec![0u8; n]` = a list of n zero bytes (u8 = one byte). We read straight into it.
        let mut raw = vec![0u8; (last - first) as usize];
        self.file.seek(SeekFrom::Start(first))?;
        self.file.read_exact(&mut raw)?;
        // The range may span several lines: drop the newline bytes, keep only bases.
        raw.retain(|&b| b != b'\n' && b != b'\r');
        Ok(raw)
    }
}

impl IndexedFasta {
    /// `search::find_motif` on sequence `name` (read once, whole): IUPAC letters, mismatches, either strand.
    pub fn search_motif(
        &mut self,
        name: &str,
        motif: &[u8],
        max_mismatches: usize,
        iupac: bool,
        both: bool,
        limit: usize,
    ) -> io::Result<Vec<crate::search::MotifHit>> {
        let length = self.entries.iter().find(|e| e.name == name).map(|e| e.length).unwrap_or(0);
        let seq = self.fetch(name, 0, length)?;
        Ok(crate::search::find_motif(&seq, motif, max_mismatches, iupac, both, false, limit))
    }

    /// Positions of `motif` (exact, case-insensitive) in sequence `name`, on both strands:
    /// (start, strand). The whole sequence is read once (one chromosome at a time).
    pub fn search(&mut self, name: &str, motif: &[u8]) -> io::Result<Vec<(u64, i8)>> {
        let length = self.entries.iter().find(|e| e.name == name).map(|e| e.length).unwrap_or(0);
        let mut seq = self.fetch(name, 0, length)?;
        seq.make_ascii_uppercase();
        let fwd = motif.to_ascii_uppercase();
        let rev = crate::primers::revcomp(&fwd);
        let mut hits = Vec::new();
        if fwd.is_empty() || fwd.len() > seq.len() {
            return Ok(hits);
        }
        for (i, w) in seq.windows(fwd.len()).enumerate() {
            if w == fwd.as_slice() {
                hits.push((i as u64, 1));
            } else if w == rev.as_slice() {
                hits.push((i as u64, -1));
            }
        }
        Ok(hits)
    }
}

/// The three bytes some programs write at the very start of a UTF-8 text file.
const UTF8_BOM: &[u8] = &[0xEF, 0xBB, 0xBF];

/// Copy FASTA `src` to `dst`, replacing the bases of sequence `name` with `seq` (wrapped at
/// 60). Streams line by line, so a 3 GB genome is never in memory. `dst` may be `src`: we
/// write a temporary file and rename it over the target at the end. The old index is
/// removed (it will be rebuilt on the next open).
pub fn write_replacing(src: &Path, dst: &Path, name: &str, seq: &[u8]) -> io::Result<()> {
    let mut tmp = dst.as_os_str().to_owned();
    tmp.push(".tmp");
    let mut found = false;
    {
        let mut out = io::BufWriter::new(File::create(&tmp)?);
        let mut reader = BufReader::new(File::open(src)?);
        let mut line = Vec::new();
        let mut skipping = false;
        let mut first_line = true;
        loop {
            line.clear();
            if reader.read_until(b'\n', &mut line)? == 0 {
                break;
            }
            // A file saved by some Windows programs starts with a byte order mark (EF BB BF): it is kept in the
            // copy, but it is not part of the first header.
            let bom = if first_line && line.starts_with(UTF8_BOM) { UTF8_BOM.len() } else { 0 };
            first_line = false;
            if line[bom] == b'>' {
                let header = String::from_utf8_lossy(&line[bom + 1..]);
                skipping = header.split_whitespace().next() == Some(name);
                out.write_all(&line)?;
                // The new bases end their lines like this header does ("\r\n" in a Windows file), so the record
                // does not come out with mixed line endings.
                let eol: &[u8] = if line.ends_with(b"\r\n") { b"\r\n" } else { b"\n" };
                if !line.ends_with(b"\n") {
                    out.write_all(b"\n")?; // header was the last line, without a newline
                }
                if skipping {
                    found = true;
                    for chunk in seq.chunks(60) {
                        out.write_all(chunk)?;
                        out.write_all(eol)?;
                    }
                }
            } else if !skipping {
                out.write_all(&line)?;
            }
        }
        out.flush()?;
    }
    if !found {
        let _ = std::fs::remove_file(&tmp);
        return Err(io::Error::new(io::ErrorKind::NotFound, format!("no sequence named '{name}' in {}", src.display())));
    }
    if let Err(e) = std::fs::rename(&tmp, dst) {
        let _ = std::fs::remove_file(&tmp); // leave no stray .tmp file behind
        return Err(e);
    }
    let _ = std::fs::remove_file(fai_path_for(dst));
    Ok(())
}

/// "genome.fa" -> "genome.fa.fai"
fn fai_path_for(path: &Path) -> PathBuf {
    // OsString instead of String: Windows paths are not guaranteed to be valid UTF-8.
    let mut s: OsString = path.as_os_str().to_owned();
    s.push(".fai");
    PathBuf::from(s)
}

/// Scan the whole FASTA once and compute the index entries.
///
/// Rejects files whose lines are wrapped unevenly inside a sequence: with those, the
/// offset arithmetic in `fetch` would silently return the wrong bases.
pub fn build_fai(path: &Path) -> io::Result<Vec<FaiEntry>> {
    // BufReader reads the file in large chunks, so line-by-line reading stays fast.
    let mut reader = BufReader::new(File::open(path)?);
    let mut entries: Vec<FaiEntry> = Vec::new();
    let mut line: Vec<u8> = Vec::new(); // reused buffer: no allocation per line
    let mut pos: u64 = 0; // byte offset of the start of the next line
    let mut ended = false; // saw the short last line of the current sequence?

    loop {
        line.clear();
        // read_until keeps the '\n' (and any '\r' before it) in `line`; returns bytes read.
        let n = reader.read_until(b'\n', &mut line)? as u64;
        if n == 0 {
            break; // end of file
        }
        let line_start = pos;
        pos += n;
        // A byte order mark at the start of the file is skipped when looking at the first line (the byte
        // positions recorded below still count it, so they stay right).
        let bom = if line_start == 0 && line.starts_with(UTF8_BOM) { UTF8_BOM.len() } else { 0 };

        if line[bom] == b'>' {
            // Header: the name is the first word after '>' (same rule as samtools).
            let header = String::from_utf8_lossy(&line[bom + 1..]);
            let name = header.split_whitespace().next().unwrap_or("").to_string();
            entries.push(FaiEntry { name, length: 0, offset: pos, line_bases: 0, line_width: 0 });
            ended = false;
            continue;
        }

        // `let Some(x) = ... else { ... }` = unpack the Option or run the else branch.
        // Text before the first '>' has no sequence to belong to: skip it.
        let Some(cur) = entries.last_mut() else { continue };
        let bases = line.iter().filter(|b| !b.is_ascii_whitespace()).count() as u64;
        if bases == 0 {
            ended = true; // blank line: only allowed at the end of a sequence
            continue;
        }
        if cur.line_width == 0 {
            // First sequence line defines the wrapping for this record.
            cur.line_bases = bases;
            cur.line_width = n;
        } else if ended || bases > cur.line_bases || (bases == cur.line_bases && n != cur.line_width) {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "'{}': uneven line lengths at byte {line_start}; re-wrap the file (e.g. `seqkit seq -w 60`)",
                    cur.name
                ),
            ));
        }
        if bases < cur.line_bases {
            ended = true; // a short line must be the last one of the record
        }
        cur.length += bases;
    }
    Ok(entries)
}

fn read_fai(path: &Path) -> io::Result<Vec<FaiEntry>> {
    let text = std::fs::read_to_string(path)?;
    let bad = || io::Error::new(io::ErrorKind::InvalidData, format!("malformed index {}", path.display()));
    let mut entries = Vec::new();
    for l in text.lines().filter(|l| !l.is_empty()) {
        let f: Vec<&str> = l.split('\t').collect();
        if f.len() < 5 {
            return Err(bad());
        }
        // `.parse::<u64>()` converts text to a number; `.map_err(|_| bad())` swaps in our error.
        let num = |s: &str| s.parse::<u64>().map_err(|_| bad());
        entries.push(FaiEntry {
            name: f[0].to_string(),
            length: num(f[1])?,
            offset: num(f[2])?,
            line_bases: num(f[3])?,
            line_width: num(f[4])?,
        });
    }
    Ok(entries)
}

fn write_fai(path: &Path, entries: &[FaiEntry]) -> io::Result<()> {
    let mut out = File::create(path)?;
    for e in entries {
        writeln!(out, "{}\t{}\t{}\t{}\t{}", e.name, e.length, e.offset, e.line_bases, e.line_width)?;
    }
    Ok(()) // `()` = "nothing", like Python's None
}

// Unit tests. `#[cfg(test)]` = compiled only by `cargo test`, never shipped.
#[cfg(test)]
mod tests {
    use super::*; // import everything from the parent module (this file)

    fn write_tmp(name: &str, content: &str) -> PathBuf {
        let p = std::env::temp_dir().join(name);
        std::fs::write(&p, content).unwrap(); // unwrap = "crash the test if this fails"
        let _ = std::fs::remove_file(fai_path_for(&p)); // stale index from a previous run
        p
    }

    #[test]
    fn a_fasta_rewritten_in_the_same_clock_tick_is_re_indexed() {
        // Open a file (this saves its .fai), then rewrite it SHORTER and give it exactly the index's timestamp, as happens
        // on a coarse file clock when the rewrite follows at once. The old index (record of 20 bases) must not be trusted.
        let p = write_tmp("ge_test_same_tick.fa", ">r
ACGTACGTACGTACGTACGT
");
        assert_eq!(IndexedFasta::open(&p).unwrap().entries[0].length, 20);
        std::fs::write(&p, ">r
ACGTACGTAC
").unwrap();
        let stamp = std::fs::metadata(fai_path_for(&p)).unwrap().modified().unwrap();
        std::fs::File::options().write(true).open(&p).unwrap().set_modified(stamp).unwrap();
        let mut fa = IndexedFasta::open(&p).unwrap();
        assert_eq!(fa.entries[0].length, 10, "the stale index was used");
        assert_eq!(fa.fetch("r", 0, 10).unwrap(), b"ACGTACGTAC");
    }

    #[test]
    fn fetch_across_line_breaks_with_crlf() {
        // Windows line endings + wrapping at 4 bases, two records.
        let p = write_tmp("ge_test_crlf.fa", ">chrA desc\r\nACGT\r\nTTGG\r\nCA\r\n>chrB\r\nGGGG\r\n");
        let mut fa = IndexedFasta::open(&p).unwrap();
        assert_eq!(fa.entries[0].length, 10);
        assert_eq!(fa.entries[1].length, 4);
        assert_eq!(fa.fetch("chrA", 2, 9).unwrap(), b"GTTTGGC");
        assert_eq!(fa.fetch("chrA", 8, 999).unwrap(), b"CA"); // end clamped
        assert_eq!(fa.fetch("chrB", 0, 4).unwrap(), b"GGGG");
        assert!(fa.fetch("nope", 0, 1).is_err());
        // The saved .fai must read back identically.
        assert_eq!(read_fai(&fai_path_for(&p)).unwrap(), fa.entries);
    }

    #[test]
    fn replace_one_sequence_and_search() {
        let p = write_tmp("ge_test_replace.fa", ">a first\nACGTACGTAC\nGT\n>b\nTTTT\n>c\nGGGAATTCGG\n");
        let mut fa = IndexedFasta::open(&p).unwrap();
        assert_eq!(fa.search("c", b"gaattc").unwrap(), vec![(2, 1)]); // palindrome: + only
        assert_eq!(fa.search("a", b"CGTA").unwrap(), vec![(1, 1), (3, -1), (5, 1), (7, -1)]); // rc = TACG
        let out = std::env::temp_dir().join("ge_test_replace_out.fa");
        let new_b: Vec<u8> = b"CC".repeat(40);
        write_replacing(&p, &out, "b", &new_b).unwrap();
        let mut back = IndexedFasta::open(&out).unwrap();
        assert_eq!(back.fetch("b", 0, 999).unwrap(), new_b);
        assert_eq!(back.fetch("a", 0, 99).unwrap(), b"ACGTACGTACGT"); // others untouched
        assert_eq!(back.fetch("c", 0, 99).unwrap(), b"GGGAATTCGG");
        assert!(std::fs::read_to_string(&out).unwrap().starts_with(">a first\n"));
        assert!(write_replacing(&p, &out, "zzz", b"A").is_err());
    }

    #[test]
    fn rejects_uneven_wrapping() {
        let p = write_tmp("ge_test_uneven.fa", ">x\nACG\nACGT\n");
        assert!(build_fai(&p).is_err());
        // --- a short line that is not the last of its record, and the same bases on a wider (CRLF) line ---
        let p = write_tmp("ge_test_uneven.fa", ">x\nACGT\nAC\nACGT\n");
        assert!(build_fai(&p).is_err());
        let p = write_tmp("ge_test_uneven.fa", ">x\nACGT\nACGT\r\nAC\n");
        assert!(build_fai(&p).is_err());
        // --- control: a short LAST line, a blank line and another record are fine ---
        let p = write_tmp("ge_test_uneven.fa", ">x\nACGT\nACGT\nAC\n\n>y\nA\n");
        let fai = build_fai(&p).unwrap();
        assert_eq!(fai.iter().map(|e| (e.name.as_str(), e.length)).collect::<Vec<_>>(), vec![("x", 10), ("y", 1)]);
    }
}
