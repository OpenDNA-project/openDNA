"""A GenBank file with several records (a bacterium's two chromosomes, a genome and its plasmids), treated like a multi-FASTA
file: the map lists the records, one is open at a time, and Save puts it back in place among the others.

The Rust document holds one record, so a record is cut out of the file into a temporary one-record file and read from there
(the file is streamed: nothing but the open record is held in memory). Records are told apart by their LOCUS name.
"""
import os
import tempfile
from pathlib import Path

from genomeedit_core import Document

from reference import genbank_records

GENBANK_EXT = {".gb", ".gbk", ".gbff", ".genbank", ".gen"}


def is_multi_genbank(path):
    """True for a text GenBank file with more than one LOCUS line (the file is streamed, and stops at the second)."""
    path = Path(path)
    if path.suffix.lower() not in GENBANK_EXT:
        return False
    try:
        with open(path, encoding="utf-8", errors="replace") as fh:
            n = 0
            for line in fh:
                if line.startswith("LOCUS"):
                    n += 1
                    if n > 1:
                        return True
    except OSError:
        pass
    return False


def records(path):
    """[(name, length, circular)] of the records, in file order."""
    return genbank_records(path)


def _chunks(path):
    """(LOCUS name, lines of the record) for each record of the file, one at a time."""
    name, chunk = None, []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if name is None:
                if not line.startswith("LOCUS"):
                    continue   # blank lines between records
                parts = line.split()
                name = parts[1] if len(parts) > 1 else ""
                chunk = [line]
                continue
            chunk.append(line)
            if line.startswith("//"):
                yield name, chunk
                name, chunk = None, []
        if name is not None:
            yield name, chunk


def read_record(path, name):
    """The Document of record `name` of the file."""
    for found, lines in _chunks(path):
        if found == name:
            handle, tmp = tempfile.mkstemp(suffix=".gb", prefix="genomeedit_record_")
            with os.fdopen(handle, "w", encoding="utf-8", newline="\n") as out:
                out.writelines(lines)
            try:
                return Document(tmp)
            finally:
                Path(tmp).unlink(missing_ok=True)
    raise OSError(f"{Path(path).name}: no record named {name}")


def sibling_editors(editor):
    """[(name, Editor)] for the other records of the open editor's file: each is an ordinary Editor of its record, so every tool
    that works on "the open sequence" can be run on them one after the other (Detect Features, ...). Kept on `editor`
    (reading a chromosome takes a moment); `forget_siblings` drops them when the file changes."""
    from editor import Editor   # here: editor imports this module
    cache = editor.__dict__.setdefault("_siblings", {})
    out = []
    for name, _ in editor.fasta_records:
        if name == editor.fasta_record:
            continue
        if name not in cache:
            other = Editor(read_record(editor.path, name), editor.path)
            other.fasta_record, other.fasta_records, other.multi_genbank = name, editor.fasta_records, True
            cache[name] = other
        out.append((name, cache[name]))
    return out


def forget_siblings(editor):
    editor.__dict__.pop("_siblings", None)


def save_document(path, name, doc, primers=None):
    """Write `doc` (the record `name` of the multi-record file `path`) into the file, the other records unchanged."""
    handle, tmp_name = tempfile.mkstemp(suffix=".gb", prefix="genomeedit_save_")
    os.close(handle)
    tmp = Path(tmp_name)
    try:
        doc.save_genbank(str(tmp), primers)
        write_record(path, name, tmp.read_text(encoding="utf-8"))
    finally:
        tmp.unlink(missing_ok=True)


def save_record_of(editor):
    """Write the Document of `editor` (one record of a multi-record file) back into its file, the others unchanged."""
    save_document(editor.path, editor.fasta_record, editor.doc, editor.file_primer_features())


def write_record(path, name, text):
    """Replace record `name` of the file by `text` (one GenBank record); the others pass through unchanged. Written to a
    temporary file first and renamed over the original, so a crash never leaves half a file."""
    path = Path(path)
    tmp = path.with_name(path.name + ".saving")
    replaced = False
    try:
        with open(tmp, "w", encoding="utf-8", newline="\n") as out:
            for found, lines in _chunks(path):
                if found == name and not replaced:
                    out.write(text if text.endswith("\n") else text + "\n")
                    replaced = True
                else:
                    out.writelines(lines)
        if not replaced:
            raise OSError(f"{path.name}: no record named {name}")
        os.replace(tmp, path)
    finally:
        tmp.unlink(missing_ok=True)
