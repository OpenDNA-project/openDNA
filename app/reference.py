"""The reference genome of the off-target checks (Primers, Spacer Finder): every sequence of the file, read one at a time.

A reference is a file holding one sequence (GenBank, .dna) or several (a multi-record FASTA such as a 24-chromosome genome,
or a multi-record GenBank file). `Reference.search` looks for the binding sites of a primer in each record in turn and
names the record in every hit ("chr3:1,200–1,221"), so a copy on the third chromosome is never missed. Nothing here needs Qt:
the search runs on a worker thread, reports progress and can be cancelled.

Memory: a FASTA is read through its .fai index, a record at a time, and a long linear record in windows of WINDOW bases
(that overlap, so a site on a window edge is found exactly once). A multi-record GenBank file is split into temporary
one-record files, again one at a time. Circular or linear: a GenBank / .dna record keeps its own topology, a FASTA record is
linear (the format cannot say; a plasmid's last bases then do not pair with its first).
"""

import tempfile
import time
from pathlib import Path

from genomeedit_core import Document, Fasta

FASTA_EXT = {".fa", ".fasta", ".fna", ".fas", ".ffn", ".frn", ".mfa"}
WINDOW = 32_000_000        # bases of a linear FASTA record searched in one go (memory bound: ~3 copies of this)
MAX_SECONDS = 600          # a search that takes longer stops and says how far it got
MAX_HITS = 200             # sites kept per list (the verdict only needs to know there are some)


def is_fasta_file(path):
    path = Path(path)
    if path.suffix.lower() in FASTA_EXT:
        return True
    try:
        with open(path, "rb") as fh:
            return fh.read(64).lstrip()[:1] == b">"
    except OSError:
        return False


def genbank_records(path):
    """[(name, length, circular)] of every record of a GenBank file, from the LOCUS lines (the file is streamed)."""
    out = []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("LOCUS"):
                parts = line.split()
                try:
                    out.append((parts[1], int(parts[2]), "circular" in line.lower()))
                except (IndexError, ValueError):
                    out.append((parts[1] if len(parts) > 1 else f"record{len(out) + 1}", 0, False))
    return out


def split_genbank(path, folder):
    """Yield (name, Path of a temporary one-record file) for each record of a multi-record GenBank file, one at a time."""
    chunk, n = [], 0
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            chunk.append(line)
            if line.startswith("//"):
                n += 1
                part = Path(folder) / f"record{n}.gb"
                part.write_text("".join(chunk), encoding="utf-8")
                chunk = []
                yield part
                part.unlink(missing_ok=True)


class GenomeHits:
    """What a search found. `whole`: places where the whole primer matches exactly; `hits`: where its 3' end matches
    (within the mismatches allowed). Each is (record name, start, end, strand, mismatches, record length, circular); `bases` holds
    the template's bases of each listed site, so the Primers view can show which of the primer's bases do not pair."""

    def __init__(self, total):
        self.whole, self.hits = [], []
        self.total, self.searched = total, 0   # sequences in the reference / looked at so far
        self.skipped = []                      # names of the records left out (the open sequence itself)
        self.stopped = ""                      # "" = complete, else "time" or "cancelled"
        self.n_whole = self.n_hits = 0         # real counts (the lists above stop at MAX_HITS)
        self.bases = {}                        # (name, start, end, strand) -> the template's bases there (top strand)

    @property
    def complete(self):
        return not self.stopped


class Reference:
    """A reference file: `records` = [(name, length, circular)], read lazily with `search` / `sequences`."""

    def __init__(self, path):
        self.path = Path(path)
        if not self.path.exists():
            raise OSError(f"{self.path}: no such file")
        self.fasta = is_fasta_file(self.path)
        self._doc = None
        self._cache = {}
        if self.fasta:
            self.records = [(n, int(l), False) for n, l in Fasta(str(self.path)).sequences()]
            if not self.records:
                raise OSError(f"{self.path.name}: no FASTA sequence found")
        else:
            try:
                self._doc = Document(str(self.path))
                self.records = [(self._doc.name, self._doc.length, self._doc.circular)]
            except OSError as e:
                if "GenBank records" not in str(e):
                    raise
                self.records = genbank_records(self.path)   # several records: read one at a time
        self.single = len(self.records) == 1

    # -- what the file is
    @property
    def count(self):
        return len(self.records)

    @property
    def length(self):
        return sum(l for _, l, _ in self.records)

    @property
    def name(self):
        return self.records[0][0] if self.single else f"{self.count} sequences"

    @property
    def circular(self):
        return self.single and self.records[0][2]

    def label(self):
        """'BW25113 (4.6 Mb)' for one sequence, '24 sequences, 3.1 Gb' for several."""
        size = self.length
        text = (f"{size / 1e9:,.1f} Gb" if size >= 1e9 else f"{size / 1e6:,.1f} Mb" if size >= 1e6 else f"{size:,} bp")
        return f"{self.records[0][0]} ({text})" if self.single else f"{self.count} sequences, {text}"

    # -- reading
    def sequences(self):
        """Yield (name, text, circular), one record at a time (the caller must not keep them: that is the memory bound)."""
        if self.fasta:
            fa = Fasta(str(self.path))
            for name, length, circular in self.records:
                yield name, fa.fetch(name, 0, length), circular
        elif self._doc is not None:
            yield self._doc.name, self._doc.sequence(0, None), self._doc.circular
        else:
            with tempfile.TemporaryDirectory() as folder:
                for part in split_genbank(self.path, folder):
                    d = Document(str(part))
                    yield d.name, d.sequence(0, None), d.circular

    def _windows(self, fa, name, length, overlap):
        """Yield (first base owned, bases) pieces of a linear FASTA record; a piece owns [first, first + WINDOW)."""
        for first in range(0, max(length, 1), WINDOW):
            lo, hi = max(0, first - overlap), min(length, first + WINDOW + overlap)
            yield first, lo, fa.fetch(name, lo, hi)

    def search(self, sequence, core_part, mismatches, seed, skip=None, progress=None, cancel=None, seconds=MAX_SECONDS):
        """Binding sites of the primer `sequence` (whole annealing part `core_part`, no 5' tail) in every record.
        `skip` = (name, length) of a record to leave out (the open sequence: never compared with itself).
        `progress(done bases, total bases, record name, record number, records)`; `cancel()` -> True stops the search."""
        key = (sequence.upper(), core_part.upper(), mismatches, seed, skip)
        if key in self._cache:
            return self._cache[key]
        result = GenomeHits(self.count)
        total, done, started = self.length, 0, time.monotonic()
        overlap = max(1000, len(sequence) + 50)
        fa = Fasta(str(self.path)) if self.fasta else None
        todo = [(n, l, c) for n, l, c in self.records if (n, l) != skip]
        result.skipped = [n for n, l, _ in self.records if (n, l) == skip]
        result.total = len(todo)

        def add(name, length, circular, sites, whole, text, offset):
            for s, e, strand, mm in sites:
                result.n_whole += whole
                result.n_hits += not whole
                lst = result.whole if whole else result.hits
                if len(lst) < MAX_HITS:
                    lst.append((name, s + offset, e + offset, strand, mm, length, circular))
                    # (a site of a circular record may run through its origin)
                    result.bases[(name, s + offset, e + offset, strand)] = text[s:e] + (text[:e - len(text)] if e > len(text) else "")

        def search_text(name, length, circular, text, offset=0, keep=None):
            """Sites in `text` (a piece of the record starting at `offset`); `keep(start)` filters to the piece's own."""
            d = Document.new_sequence(name, text, circular)
            for sites, whole in ((d.binding_sites(core_part, 0, len(core_part)), True),
                                 (d.binding_sites(sequence, mismatches, seed), False)):
                add(name, length, circular, [(s, e, st, mm) for s, e, st, mm in sites if keep is None or keep(s + offset)],
                    whole, text, offset)

        def stop():
            if cancel is not None and cancel():
                result.stopped = "cancelled"
            elif time.monotonic() - started > seconds:
                result.stopped = "time"
            return bool(result.stopped)

        def run(name, length, circular, number):
            for first, lo, text in self._windows(fa, name, length, overlap):
                if stop():
                    return
                search_text(name, length, circular, text, lo, lambda p, f=first: f <= p < f + WINDOW)
                if progress:
                    progress(done + min(first + WINDOW, length), total, name, number, len(todo))

        if self.fasta:
            for number, (name, length, circular) in enumerate(todo, 1):
                if stop():
                    break
                if progress:
                    progress(done, total, name, number, len(todo))
                run(name, length, circular, number)
                if result.stopped:
                    break
                done += length
                result.searched += 1
        else:
            wanted = {(n, l) for n, l, _ in todo}
            for number, (name, text, circular) in enumerate(self.sequences(), 1):
                if (name, len(text)) not in wanted:
                    continue
                if stop():
                    break
                if progress:
                    progress(done, total, name, number, len(todo))
                search_text(name, len(text), circular, text)
                done += len(text)
                result.searched += 1
        if result.stopped != "cancelled":   # a time-limited search is kept: running it again would stop at the same place
            self._cache[key] = result
        return result


def as_reference(ref):
    """`ref` as a Reference: a Document (an open or built sequence) becomes a one-sequence reference."""
    if ref is None or isinstance(ref, Reference):
        return ref
    return DocumentReference(ref)


class DocumentReference(Reference):
    """A Document used as a reference (tests, tools): the same interface, nothing read from a file."""

    def __init__(self, doc):
        self.path, self.fasta, self._doc, self._cache = None, False, doc, {}
        self.records = [(doc.name, doc.length, doc.circular)]
        self.single = True
