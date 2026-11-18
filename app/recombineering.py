"""Recombineering engine (λ Red, dsDNA cassette): design the primers and the recombinant genome for a cassette that
replaces a region of the genome, or check a finished fragment (or a primer pair) against a genome. No widgets here;
see docs/specs/recombineering.md."""
from types import SimpleNamespace

from opendna_core import find_motif, find_primers

from assembly import anneal_len, snapshot
from editor import reverse_complement as rc

MAX_HOMOLOGY = 50_000   # longest homology the panel offers
HOMOLOGY = 50      # nt of homology on each side of the cassette
DEFAULT_TM = 60.0   # annealing Tm of the cassette primers
MAX_MISMATCHES = 5  # check mode: mismatches tolerated in an arm, so non-specific sites show up


def _bases(doc):
    return doc.sequence(0, None).upper()


def apply(genome, start, end, bases, features):
    """A new Document: `genome` with [start, end) replaced by `bases` (and their `features`, located inside `bases`). On a
    circular genome `end` may pass the length: the region then runs through the origin, which stays where it was."""
    n = genome.length
    if not 0 <= start <= n or end < start or end - start > n or (end > n and not genome.circular):
        raise ValueError("the region to replace does not fit the genome")
    doc = snapshot(genome, f"{genome.name} (recombinant)")
    if end > n:
        doc.paste(start, n, bases, features)
        doc.paste(0, end - n, "", [])
    else:
        doc.paste(start, end, bases, features)
    return doc


def _arm(seq, i, j, circular):
    """seq[i:j], looping round the origin of a circular genome."""
    n = len(seq)
    if i < 0 or j > n:
        if not circular or j - i > n:
            raise ValueError("no room for the homology next to the region")
        return (seq + seq)[i % n:i % n + j - i]
    return seq[i:j]


def primer_rows(left, right, c, tm=DEFAULT_TM, stem="Rec"):
    """The two primers that amplify cassette `c` (bases, in the orientation it goes in) with the homologies as 5' tails."""
    if not c:
        return []
    nf, nr = anneal_len(c, tm), anneal_len(rc(c), tm)
    return [(f"{stem}-F", left + c[:nf], nf, "recombineering"), (f"{stem}-R", rc(right) + rc(c)[:nr], nr, "recombineering")]


def design(genome, cassette, start, end, left=HOMOLOGY, right=HOMOLOGY, flip=False, tm=DEFAULT_TM):
    """The cassette (a Document) replaces genome[start:end]. `flip`: it goes in reverse-complemented. An empty cassette is a
    deletion (no primers: nothing to amplify), start == end an insertion. Returns SimpleNamespace(left, right: the
    homologies; rows: [(name, 5'->3' sequence, annealing length, "recombineering")]; product: the recombinant Document)."""
    g = _bases(genome)
    arms = (_arm(g, start - left, start, genome.circular), _arm(g, end, end + right, genome.circular))
    cas = snapshot(cassette)
    if flip:
        cas.reverse_complement()
    c = _bases(cas)
    rows = primer_rows(arms[0], arms[1], c, tm)
    feats = [(f.kind, f.location, f.qualifiers) for f in cas.copy_features(0, cas.length)]
    return SimpleNamespace(left=arms[0], right=arms[1], rows=rows, product=apply(genome, start, end, c, feats))


def amplicon(template, circular, fwd, rev):
    """The PCR product of primers `fwd` and `rev` (5'->3', tails allowed) on a template string, tails included."""
    f, r = (found[:1] for found in find_primers(template, circular, [fwd, rev], 1))
    if not f or not r or f[0][2] != 1 or r[0][2] != -1:
        raise ValueError("the primers do not face each other on the template")
    (_, fe, _, _), (rs, re, _, _) = f[0], r[0]
    if re < fe and not circular:
        raise ValueError("the primers do not face each other on the template")
    mid = template[fe:re] if re >= fe else template[fe:] + template[:re]
    return fwd + mid + rc(rev[:len(rev) - (re - rs)])


def check(genome, fragment, arm=HOMOLOGY, max_mismatches=MAX_MISMATCHES):
    """Where would `fragment` (a string, tails included; its first and last `arm` bases are the homologies) integrate? Each arm
    is searched on both strands with up to `max_mismatches` mismatches. The intended site is the best pair of hits that
    puts the fragment in order on the genome (either orientation); every other hit is non-specific. Returns SimpleNamespace(
    region: (start, end) or None, forward: orientation of the fragment, hits: {"left"|"right": [(start, end, strand,
    mismatches)]}, checks: [(name, level, text)], product: the recombinant Document or None)."""
    g, n = _bases(genome), genome.length
    frag = fragment.upper()
    if len(frag) < 2 * arm:
        raise ValueError("the fragment is shorter than its two homologies")
    best = None
    for forward in (True, False):
        f = frag if forward else rc(frag)
        hits = [find_motif(g, a, max_mismatches, False, True, genome.circular) for a in (f[:arm], f[-arm:])]
        for lh in (h for h in hits[0] if h[2] == 1):
            for rh in (h for h in hits[1] if h[2] == 1):
                e = rh[0] + (n if rh[0] < lh[1] else 0)
                if e - lh[1] <= n and (e == rh[0] or genome.circular):
                    cost = (lh[3] + rh[3], e - lh[1])
                    if best is None or cost < best[0]:
                        best = (cost, forward, f, hits, lh, rh, e)
    checks = []
    if best is None:
        hits = [find_motif(g, a, max_mismatches, False, True, genome.circular) for a in (frag[:arm], frag[-arm:])]
        checks.append(("Integration", "bad", "the two homologies do not land in order on the genome"))
        return SimpleNamespace(region=None, forward=True, hits=dict(zip(("left", "right"), hits)), checks=checks, product=None)
    _, forward, f, hits, lh, rh, e = best
    for side, own, found in (("left", lh, hits[0]), ("right", rh, hits[1])):
        if own[3]:
            checks.append((f"{side} homology", "bad", f"{own[3]} mismatch(es) against the genome"))
        extra = [h for h in found if h != own]
        if extra:
            checks.append((f"{side} homology", "bad", f"also matches {len(extra)} other place(s): " + ", ".join(
                f"{h[0] + 1}..{h[1]} ({h[3]} mm)" for h in extra[:5])))
    # Exact arms are the genome's own bases: replace only what lies between them. Otherwise the fragment's bases win.
    exact = not (lh[3] or rh[3])
    inner, s, t = (f[arm:-arm], lh[1], e) if exact else (f, lh[0], e + arm)
    product = apply(genome, s, t, inner, [])
    return SimpleNamespace(region=(lh[1], e), forward=forward, hits=dict(zip(("left", "right"), hits)), checks=checks, product=product)
