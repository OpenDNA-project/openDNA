"""Primer pair design around a selected target (our own scorer; no Primer3, GPL).

The forward primer starts (5') at most `window` bases before the target, the reverse primer
at most `window` bases after it, so the product always contains the whole target (with
window 0 it is exactly the target). Each candidate gets a penalty (lower is better):
- distance from the target Tm (1 per °C),
- GC outside 40-60 % (1 per 5 %),
- no G/C at the 3' end (+1), more than 3 G/C in the last 5 bases (+1),
- runs of 5+ identical bases (+2),
- hairpin / 3' self-dimer ΔG beyond the usual limits (1 per kcal/mol),
- 3' end (12 nt) matching elsewhere in the sequence (+5 per site).
A pair adds the Tm difference and their 3' cross-dimer ΔG beyond -5 kcal/mol.
"""
import re
from types import SimpleNamespace

from opendna_core import dimer, hairpin, tm

from editor import reverse_complement
from primer_panel import conditions, homology_settings
from uikit import plural

KEEP = 25  # best candidates per side that get the (slower) uniqueness check and pairing
DEFAULT_TM, DEFAULT_MIN_LEN, DEFAULT_MAX_LEN, DEFAULT_WINDOW = 60.0, 18, 25, 0   # the Primers panel's starting values and the command line's defaults


def salt(cond):
    return {k: cond[k] for k in ("na_mm", "mg_mm", "dntp_mm")}


def candidate(seq, tm_target, cond, tolerance=5):
    """Penalty and properties of one primer sequence (uppercase), or None if hopeless (its Tm more than `tolerance` °C from
    the target). `cond`: PCR conditions, read once per design run (settings reads are slow)."""
    t = tm(seq, **cond)
    if t is None or abs(t - tm_target) > tolerance:
        return None
    gc = 100 * sum(c in "GC" for c in seq) / len(seq)
    penalty = abs(t - tm_target) + max(0.0, abs(gc - 50) - 10) / 5
    notes = []
    if seq[-1] not in "GC":
        penalty += 1
        notes.append("no 3′ G/C")
    if sum(c in "GC" for c in seq[-5:]) > 3:
        penalty += 1
        notes.append("GC-rich 3′")
    if re.search(r"(A{5,}|C{5,}|G{5,}|T{5,})", seq):
        penalty += 2
        notes.append("run")
    hp, hp3 = hairpin(seq, **salt(cond))
    sd, sd3 = dimer(seq, seq, **salt(cond))
    penalty += max(0.0, -2 - hp3) + max(0.0, -3 - hp) + max(0.0, -5 - sd3) + max(0.0, -9 - sd)
    if hp3 < -2 or hp < -3:
        notes.append(f"hairpin {min(hp, hp3):.1f}")
    if sd3 < -5:
        notes.append(f"3′ dimer {sd3:.1f}")
    return {"seq": seq, "tm": t, "gc": gc, "penalty": penalty, "notes": notes}


def design(editor, start, end, tm_target=60.0, min_len=18, max_len=25, window=0, count=5, tm_tolerance=5):
    """Best primer pairs whose product contains [start, end). Returns a list of dicts:
    fwd/rev candidates, their positions, product size, pair penalty. The forward primer starts (5') at most `window` bases before
    `start`, the reverse primer starts (5') at most `window` bases after `end`: with the default 0 the product is exactly the
    target, the primers sitting on its two ends. A primer whose Tm is more than `tm_tolerance` °C from `tm_target` is not
    considered (None: any Tm, ranked by its penalty)."""
    n = editor.length
    seq = editor.doc.sequence().upper()
    cond = conditions()
    fwds, revs = [], []
    off, total = 0, n
    tolerance = float("inf") if tm_tolerance is None else tm_tolerance
    if editor.doc.circular:
        # Three laps with the target in the middle one, so primers (and the target) may pass through the origin.
        # The two search windows must not meet on the other side of the circle.
        off, total, seq = n, 3 * n, seq * 3
        window = max(0, min(window, (n - (end - start)) // 2))
        start, end = start + off, end + off
    # Forward: bases [a, a + L), a in [start - window, start].
    for a in range(max(0, start - window), start + 1):
        for length in range(min_len, max_len + 1):
            if a + length > total:
                break
            c = candidate(seq[a:a + length], tm_target, cond, tolerance)
            if c:
                c.update(start=a % n, end=a % n + length, raw=a)   # as binding sites: end may pass n
                fwds.append(c)
    # Reverse: anneals [b - L, b), b in [end, end + window] (its 5' end at b).
    for b in range(end, min(total, end + window) + 1):
        for length in range(min_len, max_len + 1):
            if b - length < 0:
                break
            c = candidate(reverse_complement(seq[b - length:b]), tm_target, cond, tolerance)
            if c:
                c.update(start=(b - length) % n, end=(b - length) % n + length, raw=b)
                revs.append(c)
    seed, _ = homology_settings()
    for side in (fwds, revs):
        side.sort(key=lambda c: c["penalty"])
        del side[KEEP:]
        for c in side:  # other places the 3' end could prime, in this sequence
            extra = len(editor.doc.binding_sites(c["seq"], 0, seed)) - 1
            if extra > 0:
                c["penalty"] += 5 * extra
                c["notes"].append(f"3′ matches {plural(extra, 'other site')}")
    pairs = []
    for f in fwds:
        for r in revs:
            if r["raw"] - f["raw"] < max(len(f["seq"]), len(r["seq"])):
                continue   # the product would be shorter than a primer
            _, cross3 = dimer(f["seq"], r["seq"], **salt(cond))
            penalty = f["penalty"] + r["penalty"] + abs(f["tm"] - r["tm"]) + max(0.0, -5 - cross3)
            notes = [f"F: {x}" for x in f["notes"]] + [f"R: {x}" for x in r["notes"]]
            if cross3 < -5:
                notes.append(f"cross-dimer {cross3:.1f}")
            pairs.append({"fwd": f, "rev": r, "size": r["raw"] - f["raw"], "penalty": penalty, "notes": notes})
    pairs.sort(key=lambda p: p["penalty"])
    return pairs[:count]


def editor_like(doc):
    """What `design` and the primer checks read of an editor: its length and its document."""
    return SimpleNamespace(length=doc.length, doc=doc)


def design_region(doc, start, end, tm_target=DEFAULT_TM, min_len=DEFAULT_MIN_LEN, max_len=DEFAULT_MAX_LEN, window=DEFAULT_WINDOW):
    """`design` as the PCR panel calls it, for the region [start, end) of `doc`: the longest allowed length is never below the shortest, and
    on a fixed place (window 0) the best pair whatever its Tm, else only primers within 5 degrees of the target."""
    return design(editor_like(doc), start, end, tm_target=tm_target, min_len=min_len, max_len=max(min_len, max_len),
                  window=window, tm_tolerance=None if window == 0 else 5)
