"""Reading a track file into a loaded track (the entry the Editor keeps: `Editor.add_track`), with no window: the window's
Import (main.load_track, which adds the questions) and the commands that draw or recap tracks (`print`, `recap`: `load`)
both go through here, so they read a file the same way.

A loaded track is a SimpleNamespace with `name`, `kind` ("track": GFF, GTF, BED, bigBed, VCF; "signal": bigWig, bedGraph,
wig; "bam"), `path`, `seqid` (the file's name for this sequence) and what its kind needs (see `make_entry`).
"""
from pathlib import Path
from types import SimpleNamespace

from PySide6.QtGui import QColor

from opendna_core import Bam, Signal, Track, is_bigwig
import track_check
import track_import

SIGNAL_EXT = (".bw", ".bigwig", ".bedgraph", ".bg", ".bdg", ".wig")


def track_kind(path):
    """How a track file is read: "bam", "signal" (bigWig by its contents, bedGraph, wig) or "track" (GFF, GTF, BED, bigBed, VCF)."""
    lower = str(path).lower()
    if lower.endswith(".bam"):
        return "bam"
    return "signal" if lower.endswith(SIGNAL_EXT) or is_bigwig(str(path)) else "track"


def read(path, kind):
    """The file read as its kind says (OSError when it cannot be)."""
    return {"bam": Bam, "signal": Signal}.get(kind, Track)(str(path))


def is_gff(found):
    """Do these types (Track.types) come from a GFF or GTF file, which has types to choose from (BED and VCF items have one each)?"""
    return not {kind for kind, _, _ in found} <= {"bed", "variant"}


def choose(track, types, fields):
    """(types, fields) of an annotation track to keep: `types` given (a project's list, or track_import.ALL), or None: every type
    but those whose items only name their sequence (a GFF `region` line); BED and VCF keep everything. `fields`: given, or every
    extra field the file has (None or ALL)."""
    found = track.types()
    if types == track_import.ALL or not is_gff(found):
        types = [kind for kind, _, _ in found]
    elif types is None:
        types = [kind for kind, _, whole in found if not whole]
    if fields in (None, track_import.ALL):
        fields = list(track_import.detected_fields(track))
    return types, fields


def make_entry(editor, kind, reader, path, seqid, name=None, types=None, fields=()):
    """The loaded track of a file read by `read`, on its sequence `seqid`: (entry, warning about its colours or "")."""
    name = name or Path(path).name
    if kind == "bam":
        return SimpleNamespace(name=name, kind="bam", path=str(path), bam=reader, seqid=seqid, length=dict(reader.sequences())[seqid],
                               cache=None, depth_pending=False, error=None, variants=None), ""
    if kind == "signal":
        return SimpleNamespace(name=name, kind="signal", path=str(path), signal=reader, seqid=seqid, cache={}), ""
    colors, unread = track_import.colors_of(reader) if "color" in fields else ({}, "")
    # wraps: items end past the length (they cross the origin of this circle), so drawing also looks one turn on
    return SimpleNamespace(name=name, kind="track", path=str(path), track=reader, seqid=seqid, types=types, fields=fields,
                           colors={raw: QColor(c) for raw, c in colors.items() if c}, show_notes="note" in fields,
                           wraps=editor.circular and bool(reader.query(seqid, editor.length, 1 << 62, 1)[0])), unread


def load(editor, path, name=None, seqid=None, types=None, fields=None, group=None):
    """Load a track file onto `editor` with no questions (the command line's --track, a project's tracks): (entry or None,
    [warnings]). Its sequence: `seqid` when the file has it, else the one named like this sequence, else its only one; types and
    fields as `choose` says; items across the origin of a linear sequence are left out, and said."""
    file_name = Path(path).name
    kind = track_kind(path)
    try:
        reader = read(path, kind)
    except OSError as e:
        return None, [f"{file_name}: not read ({e})"]
    warnings = [f"{file_name}: {note}" for note in getattr(reader, "notes", [])]
    seqs = [s for s, _ in reader.sequences()]
    seqid = seqid if seqid in seqs else track_check.match_contig(seqs, editor.doc.name, editor.fasta_record) or (seqs[0] if len(seqs) == 1 else None)
    if seqid is None:
        return None, warnings + [f"{file_name}: none of its sequences ({', '.join(seqs[:5]) or 'none'}) is this one ({editor.doc.name}): not loaded"]
    if kind == "track":
        types, fields = choose(reader, types, fields)
        if dropped := reader.keep(types, seqid, editor.length if editor.circular else None):
            warnings.append(f"{file_name}: {dropped} item(s) cross the origin and this sequence is linear: left out")
    lengths = dict(reader.sequences())
    problems = track_check.mismatches(editor, reader if kind == "track" else None, seqid, None if kind == "track" else lengths[seqid],
                                      path, signal=reader if kind == "signal" else None)
    if problems:
        warnings.append(f"{file_name}: may not match this sequence: " + " ".join(problems))
    entry, unread = make_entry(editor, kind, reader, path, seqid, name, types, fields or ())
    if unread:
        warnings.append(f"{file_name}: {unread}")
    editor.add_track(entry, group)
    return entry, warnings

