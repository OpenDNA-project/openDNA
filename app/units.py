"""Transcription units as features: an mRNA feature from a promoter to the end of its terminator, with a program-specific
qualifier (like the Cas9 guides' OpenDNA_cas9) that tells the views to draw it as a unit rather than an arrow:

    ───────━━━━━━━━━━━━━━━━━━━━━━━━━━┫      (minus strand: mirrored)
    promoter  +1 (transcription start)   terminator

a thin line over the promoter, a thicker one from the transcription start (+1) on, and a T bar where the terminator stops
it (an arrow head when nothing stops it: it runs on). Other programs see an ordinary mRNA feature.

The qualifier value is "tss=<bases from the 5' end to +1>;end=<how it ends>;stopped=<1 or 0>;exact=<1 or 0>", exact
when +1 is known to the base (an annotated TSS, a T7 / T3 / SP6 promoter): only then do the views write "+1" there. No Qt
here.
"""
from types import SimpleNamespace

TU_KEY = "OpenDNA_tu"
TU_KEYS = (TU_KEY, "OpenGenEdit_tu")   # read both: files saved before the rename to openDNA carry the old key
KIND = "mRNA"


def tu_qualifier(tss, ending, stopped, exact=False):
    return TU_KEY, (f"tss={max(0, int(tss))};end={ending.replace(';', ',')};stopped={1 if stopped else 0};"
                    f"exact={1 if exact else 0}")


def tu_of(feature):
    """SimpleNamespace(tss, ending, stopped) of a feature saved as a unit, or None."""
    value = next((v for k, v in feature.qualifiers if k in TU_KEYS), None)
    if value is None:
        return None
    fields = dict(part.split("=", 1) for part in value.split(";") if "=" in part)
    try:
        tss = int(fields.get("tss", 0))
    except ValueError:
        tss = 0
    return SimpleNamespace(tss=tss, ending=fields.get("end", ""), stopped=fields.get("stopped", "1") != "0",
                           exact=fields.get("exact", "0") == "1")


def is_tu(feature):
    return any(k in TU_KEYS for k, _ in feature.qualifiers)


def extent(f, n, circular):
    """(start, end) of a feature on the sequence, end past the origin when it wraps (join(x..n,1..y))."""
    s, e = min(a for a, _ in f.spans), max(b for _, b in f.spans)
    if circular and len(f.spans) > 1 and s == 0 and e == n:
        first = next((a for a, b in f.spans if b == n), s)
        last = next((b for a, b in f.spans if a == 0), e)
        return first, n + last
    return s, e


def geometry(f, n, circular):
    """How to draw unit `f`: SimpleNamespace(pieces [(a, b, thick)] on the sequence (none crossing the origin),
    tss (the gap where +1 starts), end (the gap at the 3' end), stopped, strand). None if `f` is not a unit."""
    tu = tu_of(f)
    if tu is None or not n:
        return None
    strand = f.strand or 1
    start, end = extent(f, n, circular)
    tss = min(tu.tss, end - start)
    q = (start + tss) % n if strand > 0 else (end - tss) % n

    def downstream(x):   # base x is at or after +1 (bases counted from the 5' end)
        d = (x - start) % n if strand > 0 else (end - 1 - x) % n
        return d >= tss

    pieces = []
    for a, b in f.spans:
        if a < q < b:
            pieces += [(a, q, strand < 0), (q, b, strand > 0)]   # thick on the side +1 reads into
        elif b > a:
            pieces.append((a, b, downstream(a)))
    three = (end % n or n) if strand > 0 else start % n
    return SimpleNamespace(pieces=pieces, tss=q, end=three, stopped=tu.stopped, strand=strand, ending=tu.ending,
                           exact=tu.exact)
