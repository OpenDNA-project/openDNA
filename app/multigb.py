"""A GenBank file with several records (a bacterium's two chromosomes, a genome and its plasmids), treated like a multi-FASTA
file: the map lists the records, one is open at a time, and Save puts it back in place among the others.

The Rust document holds one record, so a record is cut out of the file into a temporary one-record file and read from there
(the file is streamed: nothing but the open record is held in memory). Records are told apart by their LOCUS name.
"""
import os
import re
import tempfile
from pathlib import Path

from opendna_core import Document

import fileio
from reference import genbank_records, is_fasta_file


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
            handle, tmp = tempfile.mkstemp(suffix=".gb", prefix="opendna_record_")
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
    handle, tmp_name = tempfile.mkstemp(suffix=".gb", prefix="opendna_save_")
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


def record_name(name, taken=()):
    """A record name for a file: no spaces (a LOCUS / FASTA name is one word), not one of `taken` (then "_2", "_3"... is added)."""
    base = re.sub(r"\s+", "_", (name or "").strip()) or "sequence"
    taken = set(taken)
    out, n = base, 1
    while out in taken:
        n += 1
        out = f"{base}_{n}"
    return out


def renamed_genbank(doc, name, primers=None):
    """The GenBank text of `doc` (one record) with `name` as its LOCUS name; the document itself is not changed."""
    handle, tmp_name = tempfile.mkstemp(suffix=".gb", prefix="opendna_add_")
    os.close(handle)
    tmp = Path(tmp_name)
    try:
        doc.write_genbank(str(tmp), primers)
        text = tmp.read_text(encoding="utf-8")
    finally:
        tmp.unlink(missing_ok=True)
    text, count = re.subn(r"^LOCUS([ \t]+)\S+", lambda m: "LOCUS" + m.group(1) + name, text, count=1, flags=re.M)
    if not count:
        raise OSError(f"{name}: the GenBank text has no LOCUS line")
    return text if text.endswith("\n") else text + "\n"


def _write_appended(path, write_more):
    """Copy `path` to a temporary file, let `write_more(out, eol)` add to its end, then rename it over the original (a crash
    never leaves half a file). `eol` is the file's own line ending ("\\r\\n" for a Windows file), so nothing comes out mixed."""
    path = Path(path)
    tmp = path.with_name(path.name + ".saving")
    try:
        with open(path, "rb") as src, open(tmp, "wb") as out:
            first = src.readline()
            eol = b"\r\n" if first.endswith(b"\r\n") else b"\n"
            out.write(first)
            last = first
            while block := src.read(1 << 20):
                out.write(block)
                last = block
            if last and not last.endswith(b"\n"):
                out.write(eol)   # the old last line had no line ending
            write_more(out, eol)
        fileio.replace(tmp, path)
    finally:
        tmp.unlink(missing_ok=True)


def append_genbank(path, texts):
    """Add GenBank records (their texts, see `renamed_genbank`) at the end of the GenBank file `path`, which then holds several
    records. The records already in it are copied byte for byte."""
    def more(out, eol):
        for text in texts:
            out.write(text.replace("\r\n", "\n").replace("\n", eol.decode()).encode("utf-8"))
    _write_appended(path, more)


FASTA_WIDTH = 60
_BLOCK = FASTA_WIDTH * 100_000   # bases read from a document at a time (a chromosome is never copied whole into Python)


def append_fasta(path, records):
    """Add sequences [(name, Document)] at the end of the FASTA file `path` (only the bases: FASTA has no features). Its `.fai`
    index is deleted, so it is built again, with the new sequences, the next time the file is read."""
    def more(out, eol):
        for name, doc in records:
            header = f">{name}" + (f" {doc.definition.strip()}" if (doc.definition or "").strip() else "")
            out.write(header.encode("utf-8") + eol)
            for start in range(0, doc.length, _BLOCK):
                bases = doc.sequence(start, min(start + _BLOCK, doc.length))
                for i in range(0, len(bases), FASTA_WIDTH):
                    out.write(bases[i:i + FASTA_WIDTH].encode("ascii", "replace") + eol)
    _write_appended(path, more)
    Path(str(path) + ".fai").unlink(missing_ok=True)


def record_spans(path):
    """(bytes before the first record, [(name, start, end)]) of the records of a GenBank or FASTA file, in file order: each
    record runs from its LOCUS / '>' line to the next one (blank lines between records go with the one before). The file is
    scanned in blocks for the start of a line ("\\nLOCUS", "\\n>"), not read line by line, so a genome of gigabytes is quick."""
    marker = b">" if is_fasta_file(path) else b"LOCUS"
    needle = b"\n" + marker
    starts = []
    with open(path, "rb") as fh:
        head = fh.read(len(marker) + 3)
        bom = 3 if head.startswith(b"\xef\xbb\xbf") else 0
        if head[bom:bom + len(marker)] == marker:
            starts.append(bom)
        fh.seek(0)
        offset, tail = 0, b""
        while block := fh.read(1 << 22):
            data, base = tail + block, offset - len(tail)
            i = data.find(needle)
            while i != -1:
                starts.append(base + i + 1)
                i = data.find(needle, i + 1)
            tail = data[-(len(needle) - 1):]   # a line start cut in two by the block boundary is found next time
            offset += len(block)
        spans = []
        for k, start in enumerate(starts):
            fh.seek(start)
            line = fh.readline(65536).decode("utf-8", "replace")
            words = line[1:].split() if marker == b">" else line.split()[1:]   # "> chr1 x" names chr1, as the FASTA index does
            name = words[0] if words else ""
            spans.append((name, start, starts[k + 1] if k + 1 < len(starts) else offset))
    return (starts[0] if starts else offset), spans


def rewrite_records(path, names):
    """Write the file with its records in the order of `names` (record names it has; a record left out is removed). Each record
    is copied byte for byte, through a temporary file renamed over the original; a FASTA file's `.fai` is deleted (rebuilt the
    next time the file is read)."""
    path = Path(path)
    prefix, spans = record_spans(path)
    where = {n: (s, e) for n, s, e in spans}
    if len(where) != len(spans):
        raise OSError(f"{path.name}: two records have the same name")
    if not names:
        raise OSError(f"{path.name}: a file keeps at least one sequence")
    if len(set(names)) != len(names) or any(n not in where for n in names):
        raise OSError(f"{path.name}: the records do not match the file (was it changed by another program?)")
    tmp = path.with_name(path.name + ".saving")
    try:
        with open(path, "rb") as src, open(tmp, "wb") as out:
            eol = b"\r\n" if src.readline().endswith(b"\r\n") else b"\n"
            src.seek(0)
            out.write(src.read(prefix))
            for name in names:
                start, end = where[name]
                src.seek(start)
                left, last = end - start, b""
                while left > 0:
                    block = src.read(min(left, 1 << 22))
                    if not block:
                        break
                    out.write(block)
                    left -= len(block)
                    last = block
                if not last.endswith(b"\n"):
                    out.write(eol)   # the old last record had no line ending: it may not be last any more
        fileio.replace(tmp, path)
    finally:
        tmp.unlink(missing_ok=True)
    if is_fasta_file(path):
        Path(str(path) + ".fai").unlink(missing_ok=True)


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
        fileio.replace(tmp, path)
    finally:
        tmp.unlink(missing_ok=True)
