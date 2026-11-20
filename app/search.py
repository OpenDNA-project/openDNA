"""The search bar's engine: what a query finds in a sequence, as a list of hits in order of position.

A query can be
- DNA (A C G T N, and with *IUPAC codes* R Y S W K M B D H V too): found as written and, with *both strands*, as its reverse
  complement, allowing up to N *mismatches* (a "fuzzy" search);
- a protein motif (one-letter amino acids, X = any): found in the translations of the three frames of both strands;
- a name: of a feature (its label, gene, product...; with *details* any value of any qualifier), of a primer, of a restriction
  enzyme (its sites), or of an item of a loaded annotation track.

The search runs in Rust (search.rs) for DNA and proteins, so a chromosome takes about a second. `search_sequence` searches the open
sequence; `search_file` the other sequences of a FASTA file (DNA only, the names belong to the open sequence).
"""
from dataclasses import dataclass
from types import SimpleNamespace

from opendna_core import Fasta, find_motif, find_protein

import enzymes
import genetic_code
import multigb
from editor import HIDDEN_KINDS, NAME_KEYS, feature_bounds
from reference import is_fasta_file

DNA = set("ACGTN")
IUPAC = set("ACGTRYSWKMBDHVN")
AMINO = set("ACDEFGHIKLMNPQRSTVWYX*")
MAX_HITS = 20_000        # a search that finds more stops here (the counter says "20,000+")
MIN_CODE_QUERY = 5       # a query of IUPAC letters shorter than this ("cat", "SsrA") is looked up as a name, not as a DNA pattern
MIN_PROTEIN = 3          # amino acids
SCOPES = ("sequence", "selection", "file")


@dataclass(frozen=True)
class Options:
    """What to look for; see the module docstring."""
    both: bool = True            # DNA also as its reverse complement
    iupac: bool = True           # R, Y, N... in the query stand for several bases
    fuzzy: int = 0               # mismatches allowed (0 = exact)
    translations: bool = False   # a protein motif in the translations
    details: bool = False        # feature qualifiers, not only names
    enzymes: bool = False        # an enzyme name finds its sites
    primers: bool = True         # a primer name finds where it binds
    scope: str = "sequence"      # "sequence" | "selection" | "file" (every sequence of the FASTA file)

    def key(self):
        return tuple(self.__dict__.values())


def hit(start, end, kind, label, strand=0, mismatches=0, feature=None, record=None):
    """kind: dna, protein, feature, primer, enzyme or track. strand: 1 / -1 / 0 (not about a strand). feature: index of the feature."""
    return SimpleNamespace(start=start, end=end, kind=kind, label=label, strand=strand, mismatches=mismatches, feature=feature, record=record)


def compact(query):
    return "".join(query.split()).upper()


def dna_query(query, opts):
    """The query as a DNA pattern to search for, or None when it is not one (then only names, and maybe proteins, are looked up)."""
    q = compact(query)
    if not q:
        return None
    if set(q) <= DNA or (opts.iupac and set(q) <= IUPAC and len(q) >= MIN_CODE_QUERY):
        return q
    return None


def protein_query(query, opts):
    q = compact(query)
    return q if opts.translations and len(q) >= MIN_PROTEIN and set(q) <= AMINO else None


def allowed_mismatches(opts, length):
    """The mismatches that make sense for a query of this length: never more than a third of it."""
    return max(0, min(opts.fuzzy, length // 3))


def describe(h):
    """One line on a hit, for the status bar."""
    where = f"{h.start + 1:,}–{h.end:,}"
    strand = {1: " →", -1: " ←"}.get(h.strand, "")
    extra = f", {h.mismatches} mismatch{'es' if h.mismatches > 1 else ''}" if h.mismatches else ""
    return f"{h.label}  {where}{strand}{extra}" + (f"  in {h.record}" if h.record else "")


# ---- the open sequence ---------------------------------------------------------------------------------------------------------------


def enzymes_named(needle):
    """The enzymes a lower-case query names: the one called exactly that (EcoRI, not also EcoRII), else those whose name starts
    with it (3 letters or more, so typing "Bsa" finds BsaI, BsaAI...)."""
    names = lambda e: [e.name.lower(), *(a.lower() for a in e.aliases)]   # (an isoschizomer's name finds the enzyme the views show)
    exact = [e for e in enzymes.COMMON if needle in names(e)]
    return exact or ([e for e in enzymes.COMMON if any(n.startswith(needle) for n in names(e))] if len(needle) >= 3 else [])


def _feature_text(f, details):
    if details:   # every value of every qualifier (the translation is a protein, not a name), the type too
        return " ".join([f.label, f.kind] + [v for k, v in f.qualifiers if v and k != "translation"]).lower()
    return (" ".join(v for k, v in f.qualifiers if k in NAME_KEYS and v) + " " + f.label).lower()


CAP_NOTE = f"The search stopped at {MAX_HITS:,} hits: there may be more (look for something longer or more specific)."


def protein_hits(seq, protein, mm, circular):
    """`find_protein` of `seq`, and on a circular sequence also the matches that run through the origin (read in the sequence turned by
    half a lap: a match there that spans the old origin). Frame: the codon phase counted from the top strand's first base (reverse: from its last)."""
    table, n = genetic_code.table_id(), len(seq)
    hits = list(find_protein(seq, protein, mm, MAX_HITS, table))
    if circular and n > 3 * len(protein):
        half = n // 2
        for s, e, strand, m, _ in find_protein(seq[half:] + seq[:half], protein, mm, MAX_HITS, table):
            if s < n - half < e:
                s0 = (s + half) % n
                e0 = s0 + e - s
                hits.append((s0, e0, strand, m, (s0 % 3 if strand > 0 else (n - e0) % 3) + 1))
    return hits


def search_sequence(ed, query, opts):
    """Every hit of `query` in the open sequence of editor `ed`: [hit] in order of position, and a note (what was not searched,
    why nothing was looked for...), as (hits, note). Hits past MAX_HITS are dropped."""
    q, note = query.strip(), ""
    hits = []
    if not q:
        return hits, note
    seq, circular = ed.doc.sequence(), ed.doc.circular
    motif, protein = dna_query(q, opts), protein_query(q, opts)
    if motif:
        mm = allowed_mismatches(opts, len(motif))
        for s, e, strand, m in find_motif(seq, motif, mm, opts.iupac, opts.both, circular, MAX_HITS):
            hits.append(hit(s, e, "dna", motif, strand, m))
    if protein:
        mm = allowed_mismatches(opts, len(protein))
        for s, e, strand, m, frame in protein_hits(seq, protein, mm, circular):
            hits.append(hit(s, e, "protein", f"{protein} (frame {'+' if strand > 0 else '−'}{frame})", strand, m))
    needle = q.lower()
    for f in ed.features:
        if needle in _feature_text(f, opts.details):
            hits.append(hit(*feature_bounds(f), "feature", f.label or f.kind, f.strand, feature=f.index))
    for t in ed.tracks:
        for name, s, e, *_ in t.track.search(q, 500):
            if name == t.seqid:
                hits.append(hit(s, e, "track", q))
    if opts.primers:
        for site in ed.primer_sites:
            if needle in site.label.lower():
                hits.append(hit(site.start, site.end, "primer", site.label, site.strand, site.mismatches))
    if opts.enzymes:
        for e in enzymes_named(needle):
            hits += [hit(s, en, "enzyme", f"{e.name} site", strand) for s, en, strand in enzymes.site_spans(seq, circular, e)]
    if opts.scope == "selection":
        if ed.end > ed.start:
            laps = (0, ed.length) if getattr(ed, "circular", False) else (0,)   # (a selection through the origin runs past the end: a hit may be one lap on)
            hits = [h for h in hits if any(ed.start <= h.start + lap and h.end + lap <= ed.end for lap in laps)]
        else:
            note = "Nothing is selected: searched the whole sequence."
    seen, out = set(), []
    for h in sorted(hits, key=lambda h: (h.start, h.end, h.kind, h.label)):
        key = (h.start, h.end, h.kind, h.label)   # the two strands of a palindrome or an ambiguous site are one hit
        if key not in seen:
            seen.add(key)
            out.append(h)
    if len(out) >= MAX_HITS:
        note = (note + " " + CAP_NOTE).strip()
    return out[:MAX_HITS], note


# ---- the other sequences of a FASTA file ------------------------------------------------------------------------------------------------


def search_file(path, records, current, current_seq, query, opts, stop, on_record):
    """Search every sequence named in `records` of the FASTA file `path` for the DNA `query` (nothing else is searched in a
    file: names belong to the open sequence). `current` / `current_seq`: the open record and its bases, which may have unsaved
    edits and so are not read from the file. Calls `on_record(name, [hit])` as each one is done; stops when `stop` is set.
    Runs in a worker thread: it opens its own file handle. Returns the DNA pattern searched, or None if the query is not one.
    A multi-record GenBank file is searched for names too (`current` is then not searched: the caller does it, with everything the
    open record has: primers, tracks, enzymes, translations)."""
    genbank = not is_fasta_file(path)   # a multi-record GenBank file: DNA and the names of the features of each record
    motif = dna_query(query, opts)
    if motif is None and not genbank:
        return None
    mm = allowed_mismatches(opts, len(motif or ""))
    fasta = None
    for name in records:
        if stop.is_set():
            break
        if genbank:
            doc = multigb.read_record(path, name)
            hits = []
            if motif:
                hits += [hit(s, e, "dna", motif, strand, m, record=name)
                         for s, e, strand, m in find_motif(doc.sequence(), motif, mm, opts.iupac, opts.both, doc.circular, MAX_HITS)]
            if prot := protein_query(query, opts):
                for s, e, strand, m, frame in protein_hits(doc.sequence(), prot, allowed_mismatches(opts, len(prot)), doc.circular):
                    hits.append(hit(s, e, "protein", f"{prot} (frame {'+' if strand > 0 else '−'}{frame})", strand, m, record=name))
            needle = query.strip().lower()
            if opts.enzymes:
                for en in enzymes_named(needle):
                    hits += [hit(s, e, "enzyme", f"{en.name} site", strand, record=name)
                             for s, e, strand in enzymes.site_spans(doc.sequence(), doc.circular, en)]
            for f in doc.features():
                if f.kind not in HIDDEN_KINDS and needle in _feature_text(f, opts.details):
                    hits.append(hit(*feature_bounds(f), "feature", f.label or f.kind, f.strand, record=name))
            on_record(name, sorted(hits, key=lambda h: (h.start, h.end))[:MAX_HITS])
            continue
        if name == current:
            found = find_motif(current_seq, motif, mm, opts.iupac, opts.both, False, MAX_HITS)
        else:
            fasta = fasta or Fasta(str(path))
            found = fasta.search_motif(name, motif, mm, opts.iupac, opts.both, MAX_HITS)
        seen, hits = set(), []
        for s, e, strand, m in found:
            if (s, e) not in seen:
                seen.add((s, e))
                hits.append(hit(s, e, "dna", motif, strand, m, record=name))
        on_record(name, hits)
    return motif
