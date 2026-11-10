"""Plasmid Checkup: trace the transcription units of a construct and say, for each CDS, whether it
should be expressed.

Elements come from three places, in this order of trust:
- annotated: promoter / terminator / RBS features (or `regulatory` ones with a matching
  /regulatory_class), and polyA_signal for eukaryotic promoters;
- library: known parts found by Detect Features (app/parts.tsv and your own parts.tsv);
- predicted: motifs searched here, in bacterial terms:
  * sigma-70 promoters: -35 TTGACA, 15-19 bp spacer, -10 TATAAT (consensus of Harley & Reynolds
    1987), the conserved T-A...T of the -10 box required, an extended -10 (TG) counting extra;
  * intrinsic (rho-independent) terminators: a GC-rich hairpin followed by a T tract
    (d'Aubenton Carafa et al. 1990, simplified: perfect stem, G-T pairs allowed, no bulge);
  * ribosome binding sites: the Shine-Dalgarno AGGAGG (Shine & Dalgarno 1974) 4-20 bp before
    the start codon.

A transcription unit starts at a promoter (its 3' end, about where transcription starts) and runs
on its strand to the end of the first terminator downstream (around the plasmid if none, to the
end of a linear sequence). A CDS is expressed when a unit covers it whole, it has an RBS in
the transcript (bacterial promoters), and it reads from a start codon to a stop with no stop in
between. Promoters known to need an inducer or a special polymerase (lac / tac / trc, araBAD,
tet, T7...) are named in the notes. This is a check of the design, not a model of expression
levels: real terminators leak, and predicted elements are only motifs.

Everything is worked out on "oriented" coordinates: positions along the strand read 5'->3', so
the minus strand is handled like the plus strand of the reverse complement (x = length - pos).
"""
import re
from types import SimpleNamespace

from PySide6.QtCore import Qt
from PySide6.QtGui import QColor, QFont, QFontMetrics
from PySide6.QtWidgets import (
    QAbstractItemView, QHBoxLayout, QLabel, QPushButton, QSplitter, QTableWidget, QTableWidgetItem, QTreeWidget,
    QTreeWidgetItem, QVBoxLayout, QWidget,
)

import annotation
import codon_usage
import genetic_code
import tool_settings
from editor import COLOR_KEYS, reverse_complement, site_location
import units
from tools import ToolDialog
from uikit import fit_columns
import value_filter
import theme

TOOL = "Plasmid Checkup"
tool_settings.register(TOOL, [
    tool_settings.Param("predict", "Predict promoters, terminators and RBS from motifs", True,
                        tip="Off: only annotated features and library parts (Detect Features) are used."),
    tool_settings.Param("library", "Use library parts (Detect Features)", True),
    tool_settings.Param("min_promoter", "Predicted promoter: minimum score", 10, 8, 13,
                        tip="Bases matching the -35 TTGACA and -10 TATAAT boxes (12 = both perfect), +1 for an "
                            "extended -10 (TG), -0.5 per base the spacer is away from 17."),
    tool_settings.Param("min_stem", "Predicted terminator: minimum stem score", 18, 12, 40,
                        tip="Hairpin stem: 3 per G-C pair, 2 per A-T, 1 per G-T."),
    tool_settings.Param("rbs_reach", "RBS within", 20, 5, 50, suffix=" bp of the start codon"),
    tool_settings.Param("max_predict", "Predict only in sequences up to", 300, 10, 10_000, step=50, suffix=" kb",
                        tip="Motifs alone find many false hits in a genome; annotated elements are always used."),
])

PROMOTER, TERMINATOR, RBS, POLYA, ORIGIN, TSS = "promoter", "terminator", "RBS", "polyA_signal", "rep_origin", "TSS"
# An annotated transcription start: a feature of type TSS, or one named so (a misc_feature "+1", "TSS", "transcription start").
TSS_NAME = r"^\+1$|\bTSS\b|transcription start"
# Phage promoters whose +1 is fixed by their sequence: the first base after these (the G of ...TATA|GGG).
PHAGE_PLUS1 = {"T7": "TAATACGACTCACTATA", "T3": "AATTAACCCTCACTAAA", "SP6": "ATTTAGGTGACACTATA"}
TSS_REACH = 50   # bp after a promoter's end where an annotated TSS still counts as its +1
BACTERIA, EUKARYOTES = "Bacteria", "Eukaryotic cells"
REGULATORY_CLASS = {"promoter": PROMOTER, "terminator": TERMINATOR, "ribosome_binding_site": RBS,
                    "polyA_signal_sequence": POLYA}
STARTS = {"ATG", "GTG", "TTG"}
STOPS = {"TAA", "TAG", "TGA"}
BOX35, BOX10, SD = "TTGACA", "TATAAT", "AGGAGG"
SPACERS = range(15, 20)
TSS_AFTER_10 = 7                # transcription starts ~7 bp after the -10 box
STEMS, LOOPS, TAIL = range(5, 15), range(3, 10), 8
PAIR_SCORE = {("G", "C"): 3, ("C", "G"): 3, ("A", "T"): 2, ("T", "A"): 2, ("G", "T"): 1, ("T", "G"): 1}

# Verdicts, best first.
YES, LIKELY, UNLIKELY, NO = "Expressed", "Likely", "Unlikely", "Not expressed"
VERDICT_COLORS = theme.Live({YES: "status_ok", LIKELY: "status_likely", UNLIKELY: "status_warn", NO: "status_bad"})

# Promoters needing something: (pattern on the label, note, repressor/activator CDS that makes it inducible).
REGULATED = [
    (r"\bT7\b|T7 ?promoter|pT7", "needs T7 RNA polymerase (a DE3 host such as BL21(DE3))", None),
    (r"\bT3\b|\bSP6\b", "needs its phage RNA polymerase", None),
    (r"\blac|\btac\b|\btrc\b|lacUV5|Ptac|Ptrc", "IPTG-inducible", r"\blacI"),
    (r"araBAD|\bpBAD|\bPBAD|\bara\b", "arabinose-inducible, glucose-repressed", r"\baraC"),
    (r"\btet|\bPtet|TRE", "anhydrotetracycline-inducible", r"\btetR"),
    (r"\brha|PrhaB", "rhamnose-inducible", r"\brha[RS]"),
    (r"lambda P[LR]|λ ?P[LR]|\bP[LR]\b", "repressed by lambda cI (heat-inducible with cI857)", r"\bcI"),
]
# Eukaryotic promoters: they do nothing in bacteria; their units end at a polyA signal. U6 / H1 / 7SK are Pol III
# (guide RNAs, shRNAs), so a CDS behind them is not translated.
EUKARYOTIC = r"\bCMV\b|SV40 promoter|\bSV40\b|EF-?1|\bCAG\b|\bPGK\b|\bUbC\b|\bRSV\b|\bGAL1|\bTEF1?\b|\bADH1\b|" \
             r"\bCYC1\b|\bGPD\b|\bTDH3\b|\bU6\b|\bH1\b|\b7SK\b|\bCaMKII|\bsynapsin|\bhSyn|\bGFAP\b|\bpolyhedrin|\bAOX1\b"
POL3 = r"\bU6\b|\bH1\b|\b7SK\b"
FAR = 3000   # bp: a CDS this far from its promoter, with no terminator in between, gets a note

# Replication origins, by label: (pattern, name, host, copies per cell, incompatibility group, note). Copy numbers are the
# usual textbook values in E. coli; the group says which origins cannot be kept together in one cell. The first
# pattern that matches wins: specific names first, the generic ColE1-type one ("ori") last.
ORIGINS = [
    (r"\bpUC\b", "pUC", BACTERIA, "500-700", "ColE1", "high copy (a pMB1 origin with the rop gene removed and a point "
     "mutation; temperature-dependent)"),
    (r"p15A", "p15A", BACTERIA, "10-12", "p15A", "compatible with ColE1-type origins"),
    (r"oriR101|pSC101.?ts|repA101", "pSC101 ts", BACTERIA, "~5", "pSC101", "temperature-sensitive: lost at 37-42 °C"),
    (r"pSC101", "pSC101", BACTERIA, "~5", "pSC101", "low copy; needs repA (a ts variant stops replicating at 42 °C)"),
    (r"R6K", "R6K", BACTERIA, "15-20", "R6K", "replicates only in a pir+ host (pir / pir116 strains) unless pir is on the plasmid"),
    (r"oriV|RK2|\bIncP", "RK2 (IncP)", BACTERIA, "4-7", "IncP", "broad host range; needs trfA"),
    (r"RSF1010|IncQ", "RSF1010 (IncQ)", BACTERIA, "10-12", "IncQ", "broad host range; needs repA, repB, repC"),
    (r"pBBR1", "pBBR1", BACTERIA, "~30", "pBBR1", "broad host range"),
    (r"\bF ori|oriS|\bBAC\b|\bF plasmid", "F (BAC)", BACTERIA, "1", "IncF", "single copy"),
    (r"\bf1\b|f1.?ori|M13", "f1 / M13", "Phage", "—", "", "phage origin: single-stranded DNA with a helper phage, "
     "not plasmid replication"),
    (r"SV40.?ori|SV40 origin", "SV40 ori", EUKARYOTES, "—", "", "replicates only in cells expressing SV40 large T "
     "(COS-7, HEK293T)"),
    (r"EBV|oriP", "EBV oriP", EUKARYOTES, "—", "", "episomal in cells expressing EBNA1"),
    (r"2.?(µ|u|micron)|2mu", "2µ", "Yeast", "~40", "", "yeast high copy"),
    (r"\bARS|CEN", "ARS / CEN", "Yeast", "1-2", "", "yeast low copy"),
    (r"\boriC\b", "oriC", BACTERIA, "1 (chromosome)", "", "the chromosome's own origin"),
    (r"pMB1|ColE1|pBR322|\bori\b|high.copy", "pMB1 / ColE1", BACTERIA, "15-20 (500+ without rop)", "ColE1", ""),
]


# ---- coordinates --------------------------------------------------------------------------------------------------
def to_oriented(start, end, strand, n):
    """[start, end) on the sequence (end may pass the origin) -> [a, b) along `strand` read 5'->3'."""
    if strand > 0:
        return start, end
    a = (n - end) % n if n else 0
    return a, a + (end - start)


def to_genomic(a, b, strand, n):
    if strand > 0:
        s = a % n
    else:
        s = (n - b) % n
    return s, s + (b - a)


class Strand:
    """One strand read 5'->3' (the reverse complement for the minus strand), with slicing across the origin."""

    def __init__(self, seq, strand, circular):
        seq = seq.upper()
        self.seq = seq if strand > 0 else reverse_complement(seq)
        self.n, self.strand, self.circular = len(seq), strand, circular

    def slice(self, a, b):
        n = self.n
        if not self.circular:
            return self.seq[max(0, a):max(0, min(b, n))]
        if b - a >= n:
            return self.seq
        length = b - a
        a %= n
        b = a + length
        return self.seq[a:b] if b <= n else self.seq[a:] + self.seq[:b - n]

    def extended(self, pad):
        """(sequence padded with `pad` bases from the other end when circular, offset of position 0 in it)."""
        if not self.circular or self.n < 2 * pad:
            return self.seq, 0
        return self.seq[-pad:] + self.seq + self.seq[:pad], pad

    def distance(self, frm, to):
        """Bases from oriented `frm` to `to` going downstream; None if `to` is upstream on a linear sequence."""
        d = to - frm
        if self.circular:
            return d % self.n
        return d if d >= 0 else None


# ---- elements -----------------------------------------------------------------------------------------------------
def element(kind, start, end, strand, label, origin, detail="", feature=None, score=None):
    return SimpleNamespace(kind=kind, start=start, end=end, strand=strand, label=label, origin=origin, detail=detail,
                           feature=feature, score=score)


def feature_kind(f):
    """PROMOTER / TERMINATOR / RBS / POLYA for a feature that is one, else None."""
    if f.kind in (PROMOTER, TERMINATOR, RBS, POLYA):
        return f.kind
    if f.kind == "regulatory":
        cls = next((v for k, v in f.qualifiers if k == "regulatory_class"), None)
        return REGULATORY_CLASS.get((cls or "").strip())
    return None


def feature_extent(f, n, circular):
    """(start, end) of a feature on the sequence, end past the origin when it wraps."""
    s, e = min(a for a, _ in f.spans), max(b for _, b in f.spans)
    if circular and len(f.spans) > 1 and s == 0 and e == n:   # join(x..n,1..y)
        first = next((a for a, b in f.spans if b == n), s)
        last = next((b for a, b in f.spans if a == 0), e)
        return first, n + last
    return s, e


def annotated_elements(editor):
    n, circ = editor.length, editor.doc.circular
    out = []
    for f in editor.features:
        kind = feature_kind(f)
        if not kind and (f.kind == TSS or re.search(TSS_NAME, f.label or "", re.IGNORECASE)) \
                and f.kind in (TSS, "misc_feature", "regulatory", "misc_signal"):
            kind = TSS
        if kind or f.kind == ORIGIN:
            kind = kind or ORIGIN
            s, e = feature_extent(f, n, circ)
            for strand in ([f.strand] if f.strand else [1, -1]):   # no direction: acts on both strands
                out.append(element(kind, s, e, strand, f.label, "annotated", feature=f.index))
    return out


def library_elements(editor):
    out = []
    for h in annotation.detect(editor, {PROMOTER, TERMINATOR, RBS, ORIGIN}):
        if not h.known:
            out.append(element(h.kind, h.start, h.end, h.strand, h.name, "library",
                               f"{h.identity * 100:.0f} % identical to {h.source or 'a library part'}"))
    return out


def matches(a, b):
    return sum(x == y for x, y in zip(a, b))


def predict_promoters(st, min_score):
    """Sigma-70 promoters on one strand: SimpleNamespace(a, b, score, detail) in oriented coordinates, the best one
    kept where several overlap."""
    seq, off = st.extended(60)
    found = []
    for m in re.finditer(r"(?=TA...T)", seq):
        i10 = m.start()
        if not off <= i10 < off + st.n:
            continue
        m10 = matches(seq[i10:i10 + 6], BOX10)
        extended = seq[i10 - 3:i10 - 1] == "TG" if i10 >= 3 else False
        best = None
        for spacer in SPACERS:
            i35 = i10 - spacer - 6
            if i35 < 0:
                continue
            box = seq[i35:i35 + 6]
            m35 = matches(box, BOX35)
            if m35 < 3:
                continue
            score = m10 + m35 + extended - 0.5 * abs(spacer - 17)
            if best is None or score > best[0]:
                best = (score, i35, spacer, box, m35)
        if best and best[0] >= min_score:
            score, i35, spacer, box, m35 = best
            detail = (f"-35 {box} ({m35}/6), spacer {spacer} bp, -10 {seq[i10:i10 + 6]} ({m10}/6)"
                      f"{', extended -10' if extended else ''}")
            found.append(SimpleNamespace(a=i35 - off, b=i10 + 6 + TSS_AFTER_10 - off, score=score, detail=detail))
    return best_of(found)


def stem_at(seq, e, k, loop):
    """Score, G-C pairs and G-T pairs of a stem of `k` pairs whose right arm ends at `e` (exclusive) after a loop of
    `loop` bases; None if a base does not pair."""
    left = e - 2 * k - loop
    if left < 0:
        return None
    score = gc = gt = 0
    for t in range(k):
        pair = (seq[left + t], seq[e - 1 - t])
        s = PAIR_SCORE.get(pair)
        if s is None:
            return None
        score += s
        gc += s == 3
        gt += s == 1
    return score, gc, gt


def predict_terminators(st, min_stem):
    """Intrinsic terminators on one strand: a hairpin (stem >= 5 pairs, >= 3 G-C, at most one G-T) ending 0-2 bases
    before a T tract (>= 5 T in 8 bases, starting with T)."""
    seq, off = st.extended(60)
    found = []
    for m in re.finditer(r"(?<!T)(?=T)", seq):
        j = m.start()
        if not off <= j < off + st.n:
            continue
        tail = seq[j:j + TAIL].count("T")
        if tail < 5:
            continue
        best = None
        for e in range(j - 2, j + 1):
            for k in STEMS:
                for loop in LOOPS:
                    got = stem_at(seq, e, k, loop)
                    if got and got[1] >= 3 and got[2] <= 1 and got[0] >= min_stem and (best is None or got[0] > best[0]):
                        best = (got[0], e, k, loop, got[1])
        if best:
            score, e, k, loop, gc = best
            found.append(SimpleNamespace(a=e - 2 * k - loop - off, b=j + TAIL - off, score=score + tail,
                                         detail=f"hairpin {k} bp ({gc} G-C), loop {loop}, T tract {tail}/{TAIL}"))
    return best_of(found)


def best_of(found):
    """The best-scoring of overlapping predictions."""
    kept = []
    for h in sorted(found, key=lambda h: -h.score):
        if not any(h.a < k.b and k.a < h.b for k in kept):
            kept.append(h)
    return sorted(kept, key=lambda h: h.a)


def sd_match(st, cds_start, reach):
    """Best Shine-Dalgarno match in the `reach` bases before the start codon (ending >= 4 bases before it):
    (matching bases of AGGAGG, oriented start, end) or None."""
    window_a = cds_start - reach if st.circular else max(0, cds_start - reach)   # a line's slice starts at base 0
    window = st.slice(window_a, cds_start - 4)
    best = None
    for i in range(len(window) - len(SD) + 1):
        m = matches(window[i:i + 6], SD)
        if best is None or m > best[0]:
            best = (m, window_a + i, window_a + i + 6)
    return best


def overlaps_any(el, others, n):
    """`el` covers mostly the same bases as one of `others` (same kind and strand)."""
    for o in others:
        if o.kind == el.kind and o.strand == el.strand:
            ov = min(el.end, o.end) - max(el.start, o.start)
            if n:   # also compare across the origin
                ov = max(ov, min(el.end, o.end + n) - max(el.start, o.start + n),
                         min(el.end + n, o.end) - max(el.start + n, o.start))
            if ov > 0.5 * min(el.end - el.start, o.end - o.start):
                return True
    return False


def gather_elements(editor, strands):
    """Annotated + library + predicted promoters, terminators and RBS (predicted ones only where nothing known is)."""
    n = editor.length
    known = annotated_elements(editor)
    if tool_settings.value(TOOL, "library"):
        known += [e for e in library_elements(editor) if not overlaps_any(e, known, n)]
    predicted = []
    if tool_settings.value(TOOL, "predict") and n <= tool_settings.value(TOOL, "max_predict") * 1000:
        for st in strands.values():
            for kind, hits in ((PROMOTER, predict_promoters(st, tool_settings.value(TOOL, "min_promoter"))),
                               (TERMINATOR, predict_terminators(st, tool_settings.value(TOOL, "min_stem")))):
                for h in hits:
                    s, e = to_genomic(h.a, h.b, st.strand, n)
                    el = element(kind, s, e, st.strand, f"predicted {kind}", "predicted", h.detail, score=h.score)
                    if not overlaps_any(el, known, n):
                        predicted.append(el)
    return known + predicted


# ---- transcription units and verdicts -----------------------------------------------------------------------------
def is_eukaryotic(label):
    return bool(re.search(EUKARYOTIC, label or "", re.IGNORECASE))


def regulation_notes(label, cds_labels):
    """Notes on what a promoter needs, naming the repressor / activator when the plasmid carries it."""
    out = []
    for pattern, note, partner in REGULATED:
        if re.search(pattern, label or "", re.IGNORECASE):
            if partner:
                have = next((c for c in cds_labels if re.search(partner, c, re.IGNORECASE)), None)
                note += f" ({have} on the plasmid)" if have else " (the host must provide the regulator)"
            out.append(note)
            break
    return out


def ends_unit(t, eukaryotic):
    """Whether `t` ends a unit of a eukaryotic / bacterial promoter. Bacterial RNA polymerase reads through polyA signals
    and eukaryotic terminators (yeast CYC1, ADH1...); the eukaryotic polymerases ignore bacterial terminators."""
    if t.kind == POLYA:
        return eukaryotic
    return t.kind == TERMINATOR and is_eukaryotic(t.label) == eukaryotic


def oriented_of(el, n):
    return to_oriented(el.start, el.end, el.strand, n)


def transcription_units(elements, strands, cds_labels):
    """One unit per promoter: SimpleNamespace(promoter, strand, a, b (oriented: promoter start to terminator end),
    terminator or None, end_note, eukaryotic, pol3, notes). The first terminator downstream of the transcription start
    (the promoter's 3' end) ends it."""
    units = []
    for el in elements:
        if el.kind != PROMOTER:
            continue
        st = strands[el.strand]
        n = st.n
        start, tss = oriented_of(el, n)
        exact = exact_tss(el, st, elements)
        if exact is not None:
            tss = exact
        euk = is_eukaryotic(el.label)
        best = None
        for t in elements:
            if t.strand != el.strand or not ends_unit(t, euk):
                continue
            ta, tb = oriented_of(t, n)
            d = st.distance(tss, ta)
            if d is not None and (best is None or d < best[0]) and d < n:
                best = (d, t, tb - ta)
        if best:
            d, term, tlen = best
            b, end_note = tss + d + tlen, f"ends at {term.label}"
        elif st.circular:
            b, end_note = start + n, "no terminator: runs around the whole plasmid"
        else:
            b, end_note = n, "no terminator: runs to the end of the sequence"
        notes = regulation_notes(el.label, cds_labels)
        term = best[1] if best else None
        ending, ending_detail = termination(term, euk, st)
        # The unit starts with its promoter: where transcription starts inside an annotated promoter is not known
        # (many annotations run up to the RBS), so a CDS or RBS inside the promoter counts as transcribed. `tss`: +1.
        units.append(SimpleNamespace(promoter=el, strand=el.strand, a=start, b=b, tss=tss, tss_exact=exact is not None,
                                     terminator=term,
                                     end_note=end_note, eukaryotic=euk, pol3=bool(re.search(POL3, el.label or "", re.I)),
                                     notes=notes, ending=ending, ending_detail=ending_detail))
    return units


def exact_tss(promoter, st, elements):
    """The oriented position of +1 when it is known to the base, else None: an annotated TSS on the promoter's strand,
    inside it or up to TSS_REACH bp after it; or the fixed +1 of a T7 / T3 / SP6 promoter found in its sequence. Where
    +1 is only estimated (a promoter's 3' end, ~7 bp after a predicted -10 box) there is none."""
    n = st.n
    a, b = oriented_of(promoter, n)
    for t in elements:
        if t.kind == TSS and t.strand == promoter.strand:
            ta, _ = oriented_of(t, n)
            d = st.distance(a, ta)
            if d is not None and d <= b - a + TSS_REACH:
                return a + d
    lo = a - 5 if st.circular else max(0, a - 5)   # a line's slice starts at its first base
    window = st.slice(lo, b + 5)
    for box in PHAGE_PLUS1.values():
        i = window.find(box)
        if i >= 0:
            return lo + i + len(box)
    return None


def termination(term, eukaryotic, st):
    """(short word, detail) on how a unit ends. For bacteria, an annotated or library terminator is searched for the
    hairpin + T tract of an intrinsic terminator (with the predicted terminators' settings): found, it is "intrinsic";
    not found, it may be rho-dependent (or not a terminator at all), which this check cannot tell."""
    if term is None:
        return "runs on", "no terminator: transcription is not stopped"
    if term.kind == POLYA:
        return "polyA", f"polyA signal ({term.label}): cleavage and polyadenylation"
    if eukaryotic:
        return "terminator", f"eukaryotic terminator ({term.label})"
    if term.origin == "predicted":
        return "intrinsic", f"intrinsic (rho-independent): {term.detail}"
    a, b = oriented_of(term, st.n)
    piece = Strand(st.slice(a - 10, b + 10), 1, False)
    hits = predict_terminators(piece, tool_settings.value(TOOL, "min_stem"))
    if hits:
        return "intrinsic", f"intrinsic (rho-independent): {term.label}, {hits[0].detail}"
    return "rho-dependent?", (f"{term.label}: no hairpin + T tract found in it, so not intrinsic as far as this check "
                              f"can tell; rho-dependent, or not a terminator")


def cds_oriented(f, n, circular):
    s, e = feature_extent(f, n, circular)
    return to_oriented(s, e, f.strand or 1, n)


def protein_problems(dna):
    """Problems of a coding sequence read from its first base: no start / no stop / a stop before the end. The stop codons
    are those of the genetic code chosen in Tools > Settings (the start codons stay the bacterial ATG / GTG / TTG)."""
    code = genetic_code.codon_map()
    out = []
    if dna[:3] not in STARTS:
        out.append(f"starts with {dna[:3] or 'nothing'}, not a start codon")
    if len(dna) % 3:
        out.append(f"length {len(dna)} is not a multiple of 3")
    codons = [dna[i:i + 3] for i in range(0, len(dna) - 2, 3)]
    early = next((i for i, c in enumerate(codons[:-1]) if code.get(c) == "*"), None)
    if early is not None:
        out.append(f"stop codon at codon {early + 1} of {len(codons)}: truncated protein")
    if codons and code.get(codons[-1]) != "*":
        out.append("does not end with a stop codon")
    forbidden = sorted({c for c in codons if code.get(c) == genetic_code.REMOVED})   # a recoded organism does not use them
    if forbidden:
        out.append(f"uses codons this genetic code does not use: {', '.join(forbidden)}")
    return out


LOW_CAI = 0.5      # a coding sequence below this codon adaptation index gets a note
MANY_RARE = 0.05   # ... and one with this share of rare codons or more


def codon_notes(dna):
    """Notes on how well the codons of a coding sequence suit the codon usage table chosen in Tools > Settings (empty when it
    suits it, or when no table is available): the adaptation index, and the rare codons."""
    table = codon_usage.active()
    index = codon_usage.cai(dna, table) if table is not None else None
    if index is None:
        return []
    rare = codon_usage.rare_codons(dna, table)
    scored = max(1, len(dna) // 3)
    if index >= LOW_CAI and len(rare) / scored < MANY_RARE:
        return []
    shown = ", ".join(dict.fromkeys(c for _, c in rare[:50]))
    return [f"codon adaptation index {index:.2f} for {table.name}" + (f", {len(rare)} rare codons ({shown})" if rare else "")]


def kozak(st, a):
    """Kozak context of a start codon at oriented `a`: (score 0-2, text). The purine at -3 and the G at +4 matter most
    (gccRccATGG)."""
    m3, p4 = st.slice(a - 3, a - 2), st.slice(a + 3, a + 4)
    score = (m3 in "AG" and m3 != "") + (p4 == "G")
    return score, f"{st.slice(a - 6, a).lower()}{st.slice(a, a + 4)} ({['weak', 'adequate', 'strong'][score]} Kozak)"


def check_host(f, host, units, elements, strands, reach):
    """SimpleNamespace(verdict, unit, rbs, reasons) of one CDS in one host: only the units of that host's promoters."""
    strand = f.strand or 1
    st = strands[strand]
    n = st.n
    a, b = cds_oriented(f, n, st.circular)
    length = b - a
    euk = host == EUKARYOTES
    problems = protein_problems(st.slice(a, b))
    reasons, predicted = [], False
    # The closest unit reaching the start codon; one that covers the whole CDS wins over one that stops inside.
    covering, partial = None, None
    for u in units:
        if u.strand != strand or u.eukaryotic != euk:
            continue
        d = st.distance(u.a, a)
        if d is None or d >= u.b - u.a:
            continue
        if d + length <= u.b - u.a:
            if covering is None or d < covering[0]:
                covering = (d, u)
        elif partial is None or d < partial[0]:
            partial = (d, u)
    rbs = None
    if covering is None:
        if partial:
            u = partial[1]
            reasons.append(f"{u.terminator.label if u.terminator else 'the end of the sequence'} stops the transcript "
                           f"from {u.promoter.label} inside the CDS")
        else:
            blocked = blocked_by(a, strand, elements, strands, euk)
            what = "eukaryotic promoter" if euk else "bacterial promoter"
            reasons.append(f"no {what} upstream on this strand" + (f" ({blocked})" if blocked else ""))
        return SimpleNamespace(verdict=NO, unit=None, rbs=None, reasons=reasons + problems)
    d, unit = covering
    predicted |= unit.promoter.origin == "predicted"
    if unit.pol3:
        return SimpleNamespace(verdict=NO, unit=unit, rbs=None, reasons=[
            f"{unit.promoter.label} is a Pol III promoter (guide RNA / shRNA): its transcripts are not translated"])
    if euk:
        score, text = kozak(st, a)
        rbs = SimpleNamespace(label=text, origin="sequence", detail=text, weak=score == 0, oriented=(a - 6, a + 4),
                              element=None)
        verdict = YES if score else LIKELY
        if not score:
            reasons.append(f"weak Kozak context: {text}")
        if not unit.terminator:
            reasons.append("no polyA signal downstream: the mRNA is not polyadenylated")
            verdict = UNLIKELY
    else:
        rbs = find_rbs(a, strand, elements, st, reach)
        if rbs is None:
            reasons.append("no RBS before the start codon")
            verdict = UNLIKELY
        else:
            ra, _ = rbs.oriented
            if st.distance(unit.a, ra) is None or st.distance(unit.a, ra) > d:
                reasons.append(f"transcription starts after the RBS ({rbs.label}): it is not in the mRNA")
                verdict = UNLIKELY
            else:
                verdict = YES
                if rbs.weak:
                    reasons.append(f"weak RBS ({rbs.detail})")
                    verdict = LIKELY
                predicted |= rbs.origin == "predicted"
    if verdict == YES and predicted:
        verdict = LIKELY
    if predicted:
        reasons.append("relies on predicted elements")
    reasons += unit.notes
    if not unit.terminator:
        if not euk:
            reasons.append(unit.end_note)
        if d > FAR:
            reasons.append(f"{d / 1000:.1f} kb downstream of {unit.promoter.label}, with no terminator in between")
    if problems:
        reasons += problems
        verdict = {YES: UNLIKELY, LIKELY: UNLIKELY}.get(verdict, verdict)
    return SimpleNamespace(verdict=verdict, unit=unit, rbs=rbs, reasons=reasons + codon_notes(st.slice(a, b)))


RANK = {YES: 0, LIKELY: 1, UNLIKELY: 2, NO: 3}


def check_cds(f, units, elements, strands, reach, hosts=(BACTERIA,)):
    """SimpleNamespace(feature, hosts {host: check_host result}, host (the best one), verdict, unit, rbs, reasons) for
    one CDS; verdict / unit / rbs / reasons are those of the host where it does best (bacteria on a tie)."""
    by_host = {h: check_host(f, h, units, elements, strands, reach) for h in hosts}
    best = min(hosts, key=lambda h: RANK[by_host[h].verdict])
    r = by_host[best]
    return SimpleNamespace(feature=f, hosts=by_host, host=best, verdict=r.verdict, unit=r.unit, rbs=r.rbs,
                           reasons=r.reasons)


def find_rbs(cds_a, strand, elements, st, reach):
    """The RBS of a CDS starting at oriented `cds_a`: an annotated / library one ending within `reach` bases before it,
    else a Shine-Dalgarno match (>= 4 of 6 bases; 4 is weak). SimpleNamespace(label, origin, detail, weak, oriented,
    element)."""
    n = st.n
    for el in elements:
        if el.kind == RBS and el.strand == strand:
            ra, rb = oriented_of(el, n)
            gap = st.distance(rb, cds_a)
            if gap is not None and (gap <= reach or gap > n - 3):   # a few bases of overlap are fine
                return SimpleNamespace(label=el.label, origin=el.origin, detail=el.detail, weak=False, oriented=(ra, rb),
                                       element=el)
    m = sd_match(st, cds_a, reach)
    if m and m[0] >= 4:
        box = st.slice(m[1], m[2])
        detail = f"{box} ({m[0]}/6 of AGGAGG), {cds_a - m[2]} bp before the start"
        s, e = to_genomic(m[1], m[2], strand, n)
        el = element(RBS, s, e, strand, "predicted RBS", "predicted", detail, score=m[0])
        return SimpleNamespace(label=el.label, origin="predicted", detail=detail, weak=m[0] < 5, oriented=(m[1], m[2]),
                               element=el)
    return None


def blocked_by(cds_a, strand, elements, strands, eukaryotic=False):
    """A promoter upstream whose transcript a terminator stops before the CDS: "P blocked by T", else ""."""
    st = strands[strand]
    n = st.n
    best = None
    for p in elements:
        if p.kind == PROMOTER and p.strand == strand and is_eukaryotic(p.label) == eukaryotic:
            d = st.distance(oriented_of(p, n)[1], cds_a)
            if d is not None and d < n and (best is None or d < best[0]):
                best = (d, p)
    if not best:
        return ""
    d, p = best
    tss = oriented_of(p, n)[1]
    for t in elements:
        if t.strand == strand and ends_unit(t, eukaryotic):
            td = st.distance(tss, oriented_of(t, n)[0])
            if td is not None and td < d:
                return f"{p.label} is stopped by {t.label} before it"
    return ""


def run(editor):
    """SimpleNamespace(elements, units, cds, strands, note) for the open document."""
    seq, circ = editor.doc.sequence(), editor.doc.circular
    strands = {1: Strand(seq, 1, circ), -1: Strand(seq, -1, circ)}
    n = editor.length
    note = ""
    if tool_settings.value(TOOL, "predict") and n > tool_settings.value(TOOL, "max_predict") * 1000:
        note = f"Over {tool_settings.value(TOOL, 'max_predict')} kb: elements are not predicted, only annotated / library ones are used."
    elements = gather_elements(editor, strands) if n else []
    cds_features = [f for f in editor.features if f.kind == "CDS"]
    labels = [f.label for f in cds_features]
    units = transcription_units(elements, strands, labels)
    reach = tool_settings.value(TOOL, "rbs_reach")
    # Eukaryotic cells are a host of their own only when the plasmid has a eukaryotic promoter.
    hosts = (BACTERIA, EUKARYOTES) if any(u.eukaryotic for u in units) else (BACTERIA,)
    cds = [check_cds(f, units, elements, strands, reach, hosts) for f in cds_features]
    # RBS predicted for a CDS become elements too (so they can be added as features).
    for c in cds:
        rbs = c.hosts[BACTERIA].rbs
        if rbs and rbs.origin == "predicted" and not any(e is rbs.element for e in elements):
            elements.append(rbs.element)
    for u in units:
        u.cds = [c.feature.label for c in cds if any(h.unit is u for h in c.hosts.values())]
    elements.sort(key=lambda e: (e.start, e.kind))
    origins, origin_notes = check_origins(elements, labels, circ)
    return SimpleNamespace(elements=elements, units=units, cds=cds, strands=strands, note=note, hosts=hosts,
                           origins=origins, origin_notes=origin_notes)


def origin_kind(label):
    """The ORIGINS row whose pattern matches the label first (pUC before pMB1), or None."""
    return next((row for row in ORIGINS if re.search(row[0], label or "", re.IGNORECASE)), None)


def check_origins(elements, cds_labels, circular):
    """([SimpleNamespace(element, name, host, copies, group, note)], [notes on the plasmid as a whole])."""
    seen, origins = set(), []
    for el in elements:
        if el.kind != ORIGIN or (el.start, el.end) in seen:
            continue   # an origin without direction is listed once, not per strand
        seen.add((el.start, el.end))
        row = origin_kind(el.label)
        _, name, host, copies, group, note = row or (None, el.label, "?", "?", "", "origin not recognised by name")
        if name == "R6K" and any(re.search(r"\bpir\b", c, re.IGNORECASE) for c in cds_labels):
            note = "pir is on the plasmid"
        if name.startswith("RK2") and any(re.search(r"trfA", c, re.IGNORECASE) for c in cds_labels):
            note = "broad host range; trfA on the plasmid"
        origins.append(SimpleNamespace(element=el, name=name, host=host, copies=copies, group=group, note=note))
    notes = []
    if not circular:
        notes.append("linear sequence: a plasmid replicates as a circle")
    if not any(o.host == BACTERIA for o in origins):
        notes.append("no bacterial replication origin found: it will not replicate in E. coli (fine for a suicide or "
                     "integration vector, or a fragment)")
    groups = {}
    for o in origins:
        if o.group:
            groups.setdefault(o.group, []).append(o.name)
    for group, names in groups.items():
        if len(names) > 1:
            notes.append(f"{len(names)} origins of the {group} group ({', '.join(names)}): one is enough")
    if len(groups) > 1:
        notes.append("origins of different groups: the plasmid may replicate from either (a shuttle or dual-origin vector)")
    return origins, notes


def unit_extent(u, n):
    """(start, end) of a unit on the sequence (end past the origin when it wraps), capped at one turn."""
    return to_genomic(u.a, min(u.b, u.a + n), u.strand, n)


UNIT_COLOR = "#9575cd"


def unit_label(u):
    """"Ptac → rrnB T1 (intrinsic)": the promoter, the terminator and, for bacteria, how it terminates."""
    if u.terminator is None:
        return f"{u.promoter.label} → runs on"
    how = "" if u.eukaryotic else f" ({u.ending})"
    return f"{u.promoter.label} → {u.terminator.label}{how}"


def unit_feature(u, n):
    """The unit as an mRNA feature the views draw as a unit (units.py): thin over the promoter, thick from +1, a T at
    its terminator."""
    s, e = unit_extent(u, n)
    tss = min(u.tss - u.a, e - s)
    return (units.KIND, site_location(s, e, u.strand, n),
            [("label", unit_label(u)),
             ("note", f"transcription unit (Plasmid Checkup): {u.ending_detail}; CDS: {', '.join(u.cds) or 'none'}"),
             units.tu_qualifier(tss, u.ending_detail, u.terminator is not None, u.tss_exact)]
            + [(k, UNIT_COLOR) for k in COLOR_KEYS])


def element_feature(el, n):
    quals = [("label", el.label)]
    if el.detail:
        quals.append(("note", f"{el.origin}: {el.detail}"))
    kind = el.kind if el.kind != POLYA else "polyA_signal"
    return (kind, site_location(el.start % n, el.start % n + (el.end - el.start), el.strand, n), quals)


# ---- the panel ----------------------------------------------------------------------------------------------------
SMALL = 0.88   # the panel's text, relative to the application font: the tables are dense


def small_font(widget, scale=SMALL):
    font = QFont(widget.font())
    font.setPointSizeF(max(6.5, font.pointSizeF() * scale))
    return font


def _compact(view):
    """Dense rows and columns: a smaller font, one line per row, columns as wide as their content (the last one takes
    the rest), small headers."""
    view.setFont(small_font(view))
    view.setWordWrap(False)
    view.setTextElideMode(Qt.ElideRight)
    view.setEditTriggers(QAbstractItemView.NoEditTriggers)
    view.setSelectionBehavior(QAbstractItemView.SelectRows)
    view.setStyleSheet("QTableView::item, QTreeView::item { padding: 0px 3px; }")
    header = view.header() if isinstance(view, QTreeWidget) else view.horizontalHeader()
    header.setFont(small_font(view, SMALL * 0.95))
    header.setMinimumSectionSize(24)
    header.setStretchLastSection(True)
    if isinstance(view, QTableWidget):
        view.verticalHeader().hide()
        view.verticalHeader().setDefaultSectionSize(QFontMetrics(view.font()).height() + 3)
    return view


def _fit(view, stretch):
    """Columns as wide as what they hold (fit_columns: "Replicates in" and "Unit / CDS" were cut off while Notes took 1000 px
    for 300 px of text), the descriptive column `stretch` taking the room left; rows stay dense."""
    fit_columns(view, stretch=stretch, cap=260)
    if isinstance(view, QTableWidget):
        view.verticalHeader().setDefaultSectionSize(QFontMetrics(view.font()).height() + 3)


def _table(columns):
    t = QTableWidget(0, len(columns))
    t.setHorizontalHeaderLabels(columns)
    return _compact(t)


def _key(view):
    """What the current row of a table or tree is, to find it again after the view is filled anew: its first column's text
    (and its parent's, in a tree), else None."""
    if isinstance(view, QTreeWidget):
        item = view.currentItem()
        return (item.parent().text(0) if item.parent() else None, item.text(0)) if item is not None else None
    item = view.item(view.currentRow(), 0) if view.currentRow() >= 0 else None
    return (view.currentRow(), item.text()) if item is not None else None


def _restore(view, key, scroll):
    """Make the row `key` (see _key) current again, quietly (not a click: nothing is selected in the sequence), and put the
    scroll bar back where it was."""
    view.blockSignals(True)
    try:
        if isinstance(view, QTreeWidget) and key is not None:
            parent, text = key
            tops = [view.topLevelItem(i) for i in range(view.topLevelItemCount())]
            items = tops if parent is None else [t.child(j) for t in tops if t.text(0) == parent for j in range(t.childCount())]
            found = next((i for i in items if i.text(0) == text), None)
            if found is not None:
                view.setCurrentItem(found)
        elif key is not None:
            row, text = key
            rows = [r for r in range(view.rowCount()) if view.item(r, 0) is not None and view.item(r, 0).text() == text]
            row = min(rows, key=lambda r: abs(r - row)) if rows else min(row, view.rowCount() - 1)
            if row >= 0:
                view.setCurrentCell(row, 0)
    finally:
        view.blockSignals(False)
    view.verticalScrollBar().setValue(scroll)


def _box(title, widget, *buttons):
    w = QWidget()
    lay = QVBoxLayout(w)
    lay.setContentsMargins(0, 1, 0, 0)
    lay.setSpacing(1)
    head = QHBoxLayout()
    name = QLabel(title)   # a plain title, smaller than the text: the tables are the content
    name.setFont(small_font(name))
    head.addWidget(name)
    head.addStretch()
    for b in buttons:
        b.setFont(small_font(b))
        head.addWidget(b)
    lay.addLayout(head)
    lay.addWidget(widget)
    return w


def overlay_of(result, n, added=()):
    """The checkup's results as temporary features of the document (Editor.set_overlay): the units (but those already
    added to the document: `added`, their locations), and the library / predicted elements.
    [(kind, spans, strand, qualifiers)]."""
    out, seen = [], set()

    def spans(s, e):
        s %= n
        return [(s, s + (e - s))] if s + (e - s) <= n else [(s, n), (0, s + (e - s) - n)]

    for el in result.elements:
        key = (el.kind, el.start, el.end)
        if el.origin != "annotated" and key not in seen:   # an element without direction is listed per strand
            seen.add(key)
            kind, _, quals = element_feature(el, n)
            out.append((kind, spans(el.start, el.end), el.strand, quals + [("note", "temporary (Plasmid Checkup)")]))
    for u in result.units:
        kind, location, quals = unit_feature(u, n)
        if location not in added:
            out.append((kind, spans(*unit_extent(u, n)), u.strand, quals))
    return out


def clear_overlay(editor):
    try:
        editor.set_overlay(TOOL, [])
    except RuntimeError:   # the window is closing and the editor's Qt side is already gone: nothing left to clear
        pass


class CheckupDialog(ToolDialog):
    """A pane like any other (put it beside the document's map and sequence; Tools > Plasmid Checkup opens such a tab).
    While it is shown, the units and the library / predicted elements it found are drawn in the document's views as
    temporary features, under the document's display filter (Editor.set_overlay); they go when it is hidden.

    Tables: the transcription units with the CDSs each one carries (verdict per host) and those no unit carries, the
    replication origins, and the elements used (library / predicted ones ticked, to add to the document)."""
    ORIGIN_COLUMNS = ["Origin", "Location", "Replicates in", "Copies", "Notes"]
    ELEMENT_COLUMNS = ["Element", "Type", "Location", "From", "Details"]
    NONE = "Not transcribed"

    def __init__(self, parent, editor):
        super().__init__(parent)
        self.editor, self.result = editor, None
        self.setWindowTitle(TOOL)
        self.resize(760, 640)
        self.units_tree = _compact(QTreeWidget())
        self.units_tree.setUniformRowHeights(True)
        self.units_tree.setIndentation(12)
        self.units_tree.setRootIsDecorated(True)
        self.origin_table = _table(self.ORIGIN_COLUMNS)
        self.element_table = _table(self.ELEMENT_COLUMNS)
        _fit(self.origin_table, self.ORIGIN_COLUMNS.index("Notes"))
        _fit(self.element_table, self.ELEMENT_COLUMNS.index("Details"))
        value_filter.install(self.origin_table, ["Replicates in", "Copies"])   # right-click a header: filter by values
        value_filter.install(self.element_table, ["Type", "From"])
        self.units_tree.currentItemChanged.connect(self._select_tree)
        self.origin_table.currentCellChanged.connect(lambda row, *_: self._select_origin(row))
        self.element_table.currentCellChanged.connect(lambda row, *_: self._select_element(row))
        self.info = QLabel(wordWrap=True)
        self.info.setFont(small_font(self.info))
        self.origin_info = QLabel(wordWrap=True)
        self.origin_info.setFont(small_font(self.origin_info))
        settings = QPushButton("Settings…", clicked=lambda: tool_settings.SettingsDialog(self, TOOL).exec())
        settings.setFont(small_font(settings))
        add_units = QPushButton("Add units to the document", clicked=self.add_units,
                                toolTip="Add each transcription unit to the document as an mRNA feature")
        add_elements = QPushButton("Add ticked to the document", clicked=self.add_elements)
        origins = QWidget()
        lay = QVBoxLayout(origins)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.setSpacing(1)
        lay.addWidget(self.origin_table)
        lay.addWidget(self.origin_info)
        tables = QSplitter(Qt.Vertical)
        tables.addWidget(_box("Transcription units and the CDSs they carry", self.units_tree, add_units))
        tables.addWidget(_box("Replication origins", origins))
        tables.addWidget(_box("Promoters, terminators, RBS, origins", self.element_table, add_elements))
        tables.setSizes([380, 90, 200])
        top = QHBoxLayout()
        top.setContentsMargins(0, 0, 0, 0)
        top.addWidget(self.info, 1)
        top.addWidget(settings)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(4, 4, 4, 4)
        layout.setSpacing(2)
        layout.addLayout(top)
        layout.addWidget(tables, 1)
        tool_settings.notifier.changed.connect(self._settings_changed)  # bound method: dropped with the dialog
        # The temporary features go with the panel (closed pane, closed window). A lambda holding the editor only:
        # by the time `destroyed` fires the panel's Python side may be gone.
        editor_ = editor
        self.destroyed.connect(lambda *_: clear_overlay(editor_))
        self.refresh()

    def _settings_changed(self, tool):
        if tool in (TOOL, annotation.TOOL):
            self.refresh()

    def _overlay(self):
        added = {f.location for f in self.editor.features if units.is_tu(f)}
        self.editor.set_overlay(TOOL, overlay_of(self.result, self.editor.length, added))

    def showEvent(self, event):
        super().showEvent(event)
        if self.result is not None:
            self._overlay()

    def hideEvent(self, event):
        clear_overlay(self.editor)   # another tab is shown: the views there show the document as it is
        super().hideEvent(event)

    def _loc(self, s, e, strand):
        n = self.editor.length
        return f"{s % n + 1:,}–{(e - 1) % n + 1:,} {'→' if strand > 0 else '←'}"

    def refresh(self):
        self.result = r = run(self.editor)
        n = self.editor.length
        # After an edit the tables are filled anew: the row you were on and the scroll position are kept.
        views = (self.units_tree, self.origin_table, self.element_table)
        kept = [(_key(v), v.verticalScrollBar().value()) for v in views]
        self._fill_units(r, n)
        self.origin_table.setRowCount(len(r.origins))
        for i, o in enumerate(r.origins):
            el = o.element
            name = el.label if o.name == el.label else f"{el.label} ({o.name})"
            for col, v in enumerate([name, self._loc(el.start, el.end, 1), o.host, o.copies, o.note]):
                item = QTableWidgetItem(v)
                item.setToolTip(v)
                self.origin_table.setItem(i, col, item)
        self.origin_info.setText("\n".join(f"⚠ {x}" for x in r.origin_notes))
        self.origin_info.setVisible(bool(r.origin_notes))
        self.element_table.setRowCount(len(r.elements))
        for i, el in enumerate(r.elements):
            for col, v in enumerate([el.label, el.kind, self._loc(el.start, el.end, el.strand), el.origin, el.detail]):
                item = QTableWidgetItem(v)
                item.setToolTip(v)
                if col == 0 and el.origin != "annotated":
                    item.setFlags(item.flags() | Qt.ItemIsUserCheckable)
                    item.setCheckState(Qt.Checked)
                self.element_table.setItem(i, col, item)
        _fit(self.units_tree, self.units_tree.columnCount() - 1)   # (its columns change with the hosts: fitted again)
        for view, (key, scroll) in zip(views, kept):
            _restore(view, key, scroll)
        counts = {v: sum(c.verdict == v for c in r.cds) for v in VERDICT_COLORS}
        summary = ", ".join(f"{k} {v.lower()}" for v, k in counts.items() if k) or "no CDS annotated"
        self.info.setText(f"{len(r.units)} transcription units · {len(r.cds)} CDS: {summary} (best host). "
                          f"Shown in the document's views while this panel is. {r.note}")
        if self.isVisible():
            self._overlay()

    def _fill_units(self, r, n):
        """One row per unit (sorted by position), its CDSs under it; the CDSs no unit carries under "Not transcribed"."""
        tree = self.units_tree
        collapsed = {tree.topLevelItem(i).text(0) for i in range(tree.topLevelItemCount())
                     if not tree.topLevelItem(i).isExpanded()}   # folded by hand: they stay folded
        tree.clear()
        tree.setColumnCount(3 + len(r.hosts))
        tree.setHeaderLabels(["Unit / CDS", "Location"] + list(r.hosts) + ["Notes"])
        bold = QFont(tree.font())
        bold.setBold(True)

        def cds_item(parent, c):
            f = c.feature
            s, e = feature_extent(f, n, self.editor.doc.circular)
            verdicts = [c.hosts[h].verdict for h in r.hosts]
            notes = "; ".join(c.reasons)
            item = QTreeWidgetItem(parent, [f.label, self._loc(s, e, f.strand or 1)] + verdicts + [notes])
            item.setData(0, Qt.UserRole, ("cds", f.index))
            tip = "\n\n".join(f"{h}: {c.hosts[h].verdict}\n" + "\n".join(c.hosts[h].reasons) for h in r.hosts)
            for col in range(2, 3 + len(r.hosts)):
                item.setToolTip(col, tip)
            for col, v in enumerate(verdicts, start=2):
                item.setForeground(col, QColor(VERDICT_COLORS[v]))
            return item

        placed = set()
        for k, u in sorted(enumerate(r.units), key=lambda ku: unit_extent(ku[1], n)[0]):
            s, e = unit_extent(u, n)
            notes = ("eukaryotic · " if u.eukaryotic else "") + u.ending_detail + (
                f" · {'; '.join(u.notes)}" if u.notes else "") + ("" if u.tss_exact else " · +1 estimated")
            item = QTreeWidgetItem(tree, [unit_label(u), f"{self._loc(s, e, u.strand)}  {e - s:,} bp"]
                                   + [""] * len(r.hosts) + [notes])
            item.setData(0, Qt.UserRole, ("unit", k))
            item.setFont(0, bold)
            item.setToolTip(0, item.text(0))
            item.setToolTip(len(r.hosts) + 2, notes)
            for c in r.cds:
                if c.unit is u:
                    cds_item(item, c)
                    placed.add(c.feature.index)
            item.setExpanded(item.text(0) not in collapsed)
        rest = [c for c in r.cds if c.feature.index not in placed]
        if rest:
            item = QTreeWidgetItem(tree, [self.NONE, "", *([""] * len(r.hosts)), "no unit reaches them on their strand"])
            item.setFont(0, bold)
            for c in rest:
                cds_item(item, c)
            item.setExpanded(self.NONE not in collapsed)

    def _show(self, s, e, feature=None):
        """Select [s, e) in the document's views (the feature itself when it has one)."""
        n = self.editor.length
        if feature is not None:
            self.editor.select_feature(feature)
            return
        s, e = s % n, min(s % n + (e - s), n)
        self.editor.select(s, e)
        self.editor.reveal.emit(s)

    def _select_tree(self, item, _previous=None):
        what = item.data(0, Qt.UserRole) if item is not None else None
        if not what or self.result is None:
            return
        kind, key = what
        if kind == "cds":
            self.editor.select_feature(key)
        elif 0 <= key < len(self.result.units):
            self._show(*unit_extent(self.result.units[key], self.editor.length))

    def _select_origin(self, row):
        if self.result and 0 <= row < len(self.result.origins):
            el = self.result.origins[row].element
            self._show(el.start, el.end, feature=el.feature)

    def _select_element(self, row):
        if self.result and 0 <= row < len(self.result.elements):
            el = self.result.elements[row]
            self._show(el.start, el.end, feature=el.feature)

    def ticked(self):
        return [el for i, el in enumerate(self.result.elements)
                if el.origin != "annotated" and self.element_table.item(i, 0).checkState() == Qt.Checked]

    def add_elements(self):
        chosen = self.ticked()
        if chosen:
            self.editor.add_features([element_feature(el, self.editor.length) for el in chosen])

    def add_units(self):
        if self.result and self.result.units:
            self.editor.add_features([unit_feature(u, self.editor.length) for u in self.result.units])
