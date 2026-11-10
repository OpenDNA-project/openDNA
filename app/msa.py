"""Multiple sequence alignment without Qt: the Alignment class the MSA view shows, a built-in progressive aligner for the common small
case (2 to a few tens of sequences of a few kb) so the view works on a machine with nothing installed, and a reader for the
aligned FASTA that MAFFT and Clustal Omega write.

The built-in aligner is deliberately simple: sequences are compared by shared k-mers, the two most alike are aligned first,
then the sequence closest to what is already aligned is merged in, and so on. Each merge is a banded Needleman-Wunsch of the sequence
against the profile (the column counts) of the alignment so far, with a linear gap penalty. Anything bigger or more divergent
than that is what MAFFT / Clustal Omega are for.
"""
from collections import Counter

GAP = "-"
MATCH, MISMATCH, GAP_COST = 2, -1, -3
K = 6                 # k-mer length of the similarity used for the order
MAX_BUILTIN = 20_000  # bp: a longer sequence is left to MAFFT / Clustal Omega (the built-in one is pure Python)
MIN_BAND = 64         # columns either side of the diagonal the aligner looks at (more for sequences of different lengths)
BASES = "ACGT"


class Alignment:
    """Aligned rows (equal length, '-' for a gap) with the names they are shown under; row `reference` is what identity is measured
    against, or None for the consensus."""

    def __init__(self, names, rows, reference=0):
        if not rows or len(set(map(len, rows))) != 1 or len(names) != len(rows):
            raise ValueError("an alignment needs one name per row and rows of one length")
        self.names, self.rows, self.reference = list(names), [r.upper() for r in rows], reference

    @property
    def width(self):
        return len(self.rows[0])

    def column(self, c):
        return [r[c] for r in self.rows]

    def consensus(self):
        """The most frequent residue of each column; a gap where gaps are more than half the rows (ties go to the first seen)."""
        out = []
        for c in range(self.width):
            counts = Counter(self.column(c))
            gaps = counts.pop(GAP, 0)
            out.append(GAP if not counts or gaps * 2 > len(self.rows) else counts.most_common(1)[0][0])
        return "".join(out)

    def conservation(self):
        """Per column, the share of the rows that carry its most frequent residue (gaps count against), 0..1."""
        n = len(self.rows)
        return [max((v for k, v in Counter(self.column(c)).items() if k != GAP), default=0) / n for c in range(self.width)]

    def identity(self, row, against=None):
        """Percentage identity of `row` to `against` (a row, the consensus string, or by default the reference / consensus): the
        columns where both hold the same residue, over the columns where either holds a residue (so a gap counts as a difference)."""
        if against is None:
            against = self.rows[self.reference] if self.reference is not None else self.consensus()
        elif isinstance(against, int):
            against = self.rows[against]
        mine = self.rows[row]
        same = sum(a == b != GAP for a, b in zip(mine, against))
        used = sum(a != GAP or b != GAP for a, b in zip(mine, against))
        return 100.0 * same / used if used else 0.0

    def identities(self):
        return [self.identity(i) for i in range(len(self.rows))]

    def position(self, row, col):
        """The 0-based place in sequence `row` of the residue at column `col` (in a gap: of the next residue, so a click lands
        where an insertion would go)."""
        return sum(ch != GAP for ch in self.rows[row][:col])

    def column_of(self, row, pos):
        """The column holding residue `pos` of sequence `row`."""
        seen = 0
        for c, ch in enumerate(self.rows[row]):
            if ch != GAP:
                if seen == pos:
                    return c
                seen += 1
        return self.width - 1


def parse_aligned_fasta(text):
    """(names, rows) of an aligned FASTA text, as MAFFT (--auto) and Clustal Omega (--outfmt=fa) write it."""
    names, rows = [], []
    for line in text.splitlines():
        line = line.strip()
        if line.startswith(">"):
            names.append(line[1:].strip())
            rows.append([])
        elif line and rows:
            rows[-1].append(line)
    return names, ["".join(r) for r in rows]


def to_fasta(names, seqs):
    return "".join(f">{n}\n{s}\n" for n, s in zip(names, seqs))


def kmers(seq):
    return {seq[i:i + K] for i in range(len(seq) - K + 1)}


def similarity(a, b):
    """Share of the smaller k-mer set found in the other: 1 for the same sequence, near 0 for unrelated ones."""
    if not a or not b:
        return 0.0
    return len(a & b) / min(len(a), len(b))


def _profile(rows):
    """Per column: the shares of A, C, G, T among the rows that are not gaps in it (gap-only column: all zero)."""
    out = []
    for c in range(len(rows[0])):
        counts = Counter(r[c] for r in rows)
        n = sum(counts[b] for b in BASES) + sum(v for k, v in counts.items() if k not in BASES and k != GAP)
        out.append(tuple(counts[b] / n for b in BASES) if n else (0, 0, 0, 0))
    return out


def _merge(rows, seq):
    """Align `seq` against the alignment `rows`: (rows with new gap columns put in, seq with gaps). Banded, so the cost is about
    length x band: an insertion or deletion longer than the band cannot be found (shortcut: the band widens with the length
    difference only; use MAFFT for large rearrangements)."""
    prof, n, m = _profile(rows), len(rows[0]), len(seq)
    low, high = min(0, m - n) - MIN_BAND, max(0, m - n) + MIN_BAND   # allowed j - i
    NEG = -10 ** 9
    score = [[NEG] * (high - low + 1) for _ in range(n + 1)]   # score[i][j - i - low]
    back = [[0] * (high - low + 1) for _ in range(n + 1)]      # 0 diagonal, 1 gap in seq (column kept), 2 gap in the rows
    for i in range(n + 1):
        for j in range(max(0, i + low), min(m, i + high) + 1):
            k = j - i - low
            if i == 0 and j == 0:
                score[0][k] = 0
                continue
            best, how = NEG, 0
            if i and j:   # the previous cell (i-1, j-1) has the same offset k
                s = prof[i - 1]
                b = BASES.find(seq[j - 1])
                step = (s[b] * MATCH + (1 - s[b]) * MISMATCH) if b >= 0 and sum(s) else 0
                if score[i - 1][k] > NEG:
                    best, how = score[i - 1][k] + step, 0
            if i and k + 1 < len(score[i]) and score[i - 1][k + 1] > NEG and score[i - 1][k + 1] + GAP_COST > best:
                best, how = score[i - 1][k + 1] + GAP_COST, 1
            if j and k - 1 >= 0 and score[i][k - 1] > NEG and score[i][k - 1] + GAP_COST > best:
                best, how = score[i][k - 1] + GAP_COST, 2
            score[i][k], back[i][k] = best, how
    new_rows, new_seq, i, j = [[] for _ in rows], [], n, m
    while i or j:
        how = back[i][j - i - low]
        if i and how in (0, 1) and (how == 1 or j):
            for r, row in zip(new_rows, rows):
                r.append(row[i - 1])
            new_seq.append(seq[j - 1] if how == 0 else GAP)
            i, j = i - 1, j - (how == 0)
        else:
            for r in new_rows:
                r.append(GAP)
            new_seq.append(seq[j - 1])
            j -= 1
    return ["".join(reversed(r)) for r in new_rows], "".join(reversed(new_seq))


def align(names, seqs, reference=0):
    """The Alignment of `seqs` (DNA strings) by the built-in progressive aligner. Rows come back in the order given."""
    seqs = [s.upper().replace("U", "T") for s in seqs]
    if len(seqs) < 2:
        raise ValueError("Pick at least two sequences to align.")
    if any(len(s) > MAX_BUILTIN for s in seqs):
        raise ValueError(f"The built-in aligner takes sequences up to {MAX_BUILTIN:,} bp: use MAFFT or Clustal Omega for longer ones.")
    sets = [kmers(s) for s in seqs]
    sim = {(i, j): similarity(sets[i], sets[j]) for i in range(len(seqs)) for j in range(i)}
    first, second = max(sim, key=sim.get)
    order = [second, first]
    rows, seq2 = _merge([seqs[second]], seqs[first])
    rows = rows + [seq2]
    todo = set(range(len(seqs))) - set(order)
    while todo:   # the sequence most like any already aligned one is next
        nxt = max(todo, key=lambda t: max(sim[(max(t, d), min(t, d))] for d in order))
        rows, new = _merge(rows, seqs[nxt])
        rows.append(new)
        order.append(nxt)
        todo.discard(nxt)
    by_input = {src: row for src, row in zip(order, rows)}
    return Alignment(names, [by_input[i] for i in range(len(seqs))], reference)
