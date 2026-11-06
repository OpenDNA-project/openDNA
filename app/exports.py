"""File writers behind the Export menu (no Qt here, so tests call them directly), and the
Open DNA Project file.

- FASTA / plain text of the sequence or the selection;
- features as GFF3 (https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md)
  or BED (one line per part of a feature);
- an NCBI submission pair: the sequence as .fsa and the features as a 5-column feature table
  (.tbl), the input of NCBI's table2asn (https://www.ncbi.nlm.nih.gov/genbank/feature_table/);
- the proteins of the CDS features (standard genetic code) as FASTA.

Open DNA Project file (.odp): a zip holding `sequence.gb` (the document, as Save writes it)
and `project.json`, what a GenBank file cannot hold and the program otherwise keeps only in
memory: the loaded tracks and BAM files (their paths, absolute and relative to the project,
since BAM files can be gigabytes), the DNA lock and the display options. Gene groups and
primers are already in the GenBank part.
"""
import json
import re
import shutil
import tempfile
import zipfile
from pathlib import Path

from cut_sites import CAS9_KEYS
from units import TU_KEYS
from editor import COLOR_KEYS
import genetic_code
from sequence_view import cds_codons

PROJECT_EXT = ".odp"
PROJECT_VERSION = 1


def seq_id(editor):
    """The sequence's name with no spaces (FASTA / GFF / table2asn ids)."""
    return re.sub(r"[^A-Za-z0-9_.:|-]", "_", editor.fasta_record or editor.doc.name or "sequence") or "sequence"


def fasta_text(name, seq, header_extra="", width=70):
    body = "\n".join(seq[i:i + width] for i in range(0, len(seq), width))
    return f">{name}{(' ' + header_extra) if header_extra else ''}\n{body}\n"


# ---- features ---------------------------------------------------------------------------------------------------------
def _gff_escape(text):
    return "".join(f"%{ord(c):02X}" if c in "%;=&,\t\n\r" else c for c in str(text))


def gff3_text(editor):
    """Every feature (sources and primer sites excepted) as GFF3; a feature in several parts is one line per part,
    sharing its ID."""
    sid = seq_id(editor)
    lines = ["##gff-version 3", f"##sequence-region {sid} 1 {editor.length}"]
    for f in editor.features:
        attrs = [f"ID=f{f.index}", f"Name={_gff_escape(f.label)}"]
        attrs += [f"{_gff_escape(k)}={_gff_escape(v)}" for k, v in f.qualifiers
                  if v is not None and k not in COLOR_KEYS and k not in ("label", "translation", *CAS9_KEYS, *TU_KEYS)]
        strand = "+" if f.strand > 0 else "-" if f.strand < 0 else "."
        for s, e in f.spans:
            lines.append("\t".join([sid, "openDNA", f.kind, str(s + 1), str(e), ".", strand, ".", ";".join(attrs)]))
    return "\n".join(lines) + "\n"


def bed_text(editor):
    """BED6: one line per part of each feature (0-based starts)."""
    sid = seq_id(editor)
    out = []
    for f in editor.features:
        strand = "+" if f.strand > 0 else "-" if f.strand < 0 else "."
        for s, e in f.spans:
            out.append("\t".join([sid, str(s), str(e), f.label.replace("\t", " ").replace(" ", "_"), "0", strand]))
    return "\n".join(out) + ("\n" if out else "")


TBL_SKIP = set(COLOR_KEYS) | {"label", "translation", *CAS9_KEYS, *TU_KEYS}   # program-specific, or computed by table2asn


def tbl_text(editor):
    """NCBI 5-column feature table. Parts are listed 5'→3' (a minus-strand part reads start > stop); a partial
    end is marked < (5') or > (3'). The label becomes a note when nothing else names the feature."""
    lines = [f">Feature {seq_id(editor)}"]
    for f in editor.features:
        minus = f.strand < 0
        parts = [(e, s + 1) for s, e in reversed(f.spans)] if minus else [(s + 1, e) for s, e in f.spans]
        low_partial, high_partial = "<" in f.location, ">" in f.location
        # GenBank marks the low / high coordinate; the table marks the 5' (<) and 3' (>) ends.
        five, three = (high_partial, low_partial) if minus else (low_partial, high_partial)
        for i, (a, b) in enumerate(parts):
            a_txt = ("<" if i == 0 and five else "") + str(a)
            b_txt = (">" if i == len(parts) - 1 and three else "") + str(b)
            lines.append(f"{a_txt}\t{b_txt}" + (f"\t{f.kind}" if i == 0 else ""))
        quals = [(k, v) for k, v in f.qualifiers if k not in TBL_SKIP]
        if not any(k in ("gene", "product", "note", "locus_tag") for k, _ in quals) and f.label != f.kind:
            quals.append(("note", f.label))
        for k, v in quals:
            lines.append(f"\t\t\t{k}" + (f"\t{v}" if v is not None else ""))
    return "\n".join(lines) + "\n"


def write_submission(editor, path):
    """table2asn input: `path` with .fsa (sequence) and .tbl (features) suffixes; returns both paths."""
    fsa, tbl = Path(path).with_suffix(".fsa"), Path(path).with_suffix(".tbl")
    extra = "[topology=circular]" if editor.circular else ""
    fsa.write_text(fasta_text(seq_id(editor), editor.doc.sequence(), extra), encoding="utf-8")
    tbl.write_text(tbl_text(editor), encoding="utf-8")
    return fsa, tbl


# ---- proteins ---------------------------------------------------------------------------------------------------------
def cds_features(editor):
    return [f for f in editor.features if f.kind == "CDS"]


def protein_of(editor, f):
    """Amino acids of CDS `f` without the final stop: its /transl_table (else the table chosen in Tools > Settings),
    /codon_start and /transl_except are followed, as in the sequence view."""
    return "".join(aa for _, aa, _, _ in cds_codons(editor.doc, f)).rstrip("*")


def proteins_fasta(editor):
    out = []
    for f in cds_features(editor):
        name = re.sub(r"\s+", "_", f.label)
        out.append(fasta_text(name, protein_of(editor, f), f"[{f.location}]"))
    return "".join(out)


def translate(seq, table=None):
    """Amino acids of `seq` (frame +1), with the genetic code chosen in Tools > Settings (`table` overrides it)."""
    seq = seq.upper()
    code = genetic_code.codon_map(table)
    return "".join(code.get(seq[i:i + 3], "X") for i in range(0, len(seq) - 2, 3))


# ---- project file -----------------------------------------------------------------------------------------------------
def _relative(path, base):
    try:
        return str(Path(path).resolve().relative_to(Path(base).resolve().parent))
    except ValueError:
        return None   # another drive: absolute path only


def project_meta(editor, project_path, layout=None):
    ed = editor
    return {
        "version": PROJECT_VERSION,
        "layout": layout,   # the tabs this project opens with (every other file opens with the standard ones)
        "dna_locked": ed.dna_locked,
        "hidden_primers": sorted(ed.hidden_primers),
        "tracks": [{"path": t.path, "relative": _relative(t.path, project_path), "seqid": t.seqid}
                   for t in ed.tracks if getattr(t, "path", None)],
        "alignments": [{"path": a.path, "relative": _relative(a.path, project_path), "seqid": a.seqid}
                       for a in ed.alignments if getattr(a, "path", None)],
        "display": {"enzymes": ed.show_enzymes, "enzyme_filter": ed.enzyme_filter._asdict(),
                    "hidden_features": sorted(list(k) for k in ed.hidden_features), "translation": ed.show_translation, "complement": ed.show_complement, "collapse_genes": ed.collapse_genes, "cds_only": ed.cds_only,
                    "hidden_kinds": sorted(ed.hidden_kinds), "primers": ed.show_primers, "file_primers": ed.show_file_primers,
                    "db_primers": ed.show_db_primers, "units": ed.show_units,
                    "cds_translation": ed.show_cds_translation, "untranslated": sorted(ed.untranslated)},
    }


def save_project(editor, path, layout=None):
    """Write the document and its in-memory extras as an Open DNA Project (.odp). `layout`: the window's tabs."""
    path = Path(path)
    with tempfile.TemporaryDirectory() as tmp:
        gb = Path(tmp) / "sequence.gb"
        editor.doc.write_genbank(str(gb), editor.file_primer_features())
        part = path.with_name(path.name + ".part")   # a crash mid-write never leaves a broken project
        with zipfile.ZipFile(part, "w", zipfile.ZIP_DEFLATED) as z:
            z.write(gb, "sequence.gb")
            z.writestr("project.json", json.dumps(project_meta(editor, path, layout), indent=1))
        part.replace(path)


def read_project(path, folder):
    """Unpack project `path` into `folder`: (GenBank file path, metadata dict). OSError when it is not a project."""
    try:
        with zipfile.ZipFile(path) as z:
            meta = json.loads(z.read("project.json").decode("utf-8"))
            gb = Path(folder) / "sequence.gb"
            with z.open("sequence.gb") as src, open(gb, "wb") as dst:
                shutil.copyfileobj(src, dst)
    except (zipfile.BadZipFile, KeyError, ValueError) as e:
        raise OSError(f"not an Open DNA Project ({e})") from None
    return gb, meta


def open_project(path):
    """(Document, metadata) of project `path` (the GenBank part is read from a temporary folder)."""
    from opendna_core import Document
    with tempfile.TemporaryDirectory() as tmp:
        gb, meta = read_project(path, tmp)
        return Document(str(gb)), meta


def project_file(entry, project_path):
    """Where a track listed in a project is now: its absolute path, else beside the project (moved together)."""
    if Path(entry["path"]).exists():
        return entry["path"]
    if entry.get("relative") and (p := Path(project_path).parent / entry["relative"]).exists():
        return str(p)
    return None

