"""Synthesis Success Calculator: will this DNA fragment be made quickly by a synthesis company?

A port of the model of

    Halper SM, Hossain A, Salis HM. Synthesis success calculator: predicting the rapid synthesis of DNA
    fragments with machine learning. ACS Synth Biol 2020;9(7):1563-1571. doi:10.1021/acssynbio.9b00460

Original code: https://github.com/hsalis/SalisLabCode (Synthesis_Success_Calculator), GPL v3,
(c) Sean Halper, Ayaan Hossain and Howard Salis. This file follows `seq_eval.py` (features), `RepeatFinder.py`
(repeats, in core/src/salis.rs), `PyVRNA.py` (reading the folds) and `SSC_forest.pkl` (the trained random forest,
exported to app/salis_forest.json.gz by scripts/salis_export_forest.py). Melting temperatures follow Biopython's
`Tm_NN(seq, dnac1=250.0, saltcorr=7)` as the original calls it.

The forest looks at 9 numbers: size, total_GC, longest_repeat, repeat_9_metric, richest_hairpin, GC_short_l,
GC_term_l, Tm_low and Tm_high. The results are the original's, to the last digit (tests/test_synthesis_salis.py
compares them with the original's own output on the sequences of the paper). What is the original's choice, kept
as published, is listed here so nobody has to rediscover it:

- the model reads the sequence as a straight piece of DNA, never a circle;
- a 20-base window is counted only if it starts before size - 20, so the very last window is never looked at;
  the hairpin folds likewise stop one window early (size - 101);
- hairpins are folded with ViennaRNA's DNA parameters (Mathews 2004, no GU pairs, G-quadruplexes on). ViennaRNA is
  NOT bundled (its licence does not allow redistribution for a fee and asks commercial products to ask the authors):
  the user installs it (pip install ViennaRNA), and without it the score cannot be computed.
"""
import gzip
import importlib.util
import json
import math
import struct
from bisect import bisect
from collections import deque
from pathlib import Path

import opendna_core
from editor import complement, reverse_complement
from user_error import UserError

CITATION_SHORT = "Salis, ACS Synth Biol 2020"
CITATION = ("Halper SM, Hossain A, Salis HM. Synthesis success calculator: predicting the rapid synthesis of DNA fragments "
            "with machine learning. ACS Synth Biol 2020;9(7):1563-1571. doi:10.1021/acssynbio.9b00460")
FEATURES = ("longest_repeat", "repeat_9_metric", "richest_hairpin", "GC_short_l", "GC_term_l", "Tm_low", "Tm_high", "total_GC", "size")
VIENNA_HELP = "ViennaRNA is not installed. Install it with:  pip install ViennaRNA   (openDNA does not bundle it)"

_MODEL = None


class Cancelled(Exception):
    pass


def _model():
    global _MODEL
    if _MODEL is None:
        with gzip.open(Path(__file__).with_name("salis_forest.json.gz"), "rt", encoding="ascii") as f:
            _MODEL = json.load(f)
    return _MODEL


def have_vienna():
    """True when ViennaRNA can be imported: from the environment, or from the copy the Synthesis tab installed in the data folder."""
    import vienna_install
    vienna_install.activate()
    return importlib.util.find_spec("RNA") is not None


# ---- sequence and small pieces of the original --------------------------------------------------------------



def clean(seq):
    """Upper-case A/C/G/T text; anything else is a UserError that names the first bad base (nothing is dropped)."""
    s = "".join(str(seq).split()).upper()
    if not s:
        raise UserError("There is no sequence to test.")
    for i, b in enumerate(s):
        if b not in "ACGT":
            raise UserError(f"Base {i + 1} is '{b}': the calculator accepts A, C, G and T only.")
    return s


def sequence_from_text(text):
    """The bases of a pasted text: bare bases, or one FASTA record (its '>' line is dropped). UserError for several records or a base the
    calculator cannot read (see clean)."""
    lines = str(text).splitlines()
    if sum(line.lstrip().startswith(">") for line in lines) > 1:
        raise UserError("The box holds more than one FASTA record: test one sequence at a time.")
    return clean("".join(line for line in lines if not line.lstrip().startswith(">")))


# what each of the nine numbers counts, for the tab and the recipe
FEATURE_HELP = {
    "size": ("Size", "bases"),
    "total_GC": ("GC content", "fraction of G and C over the whole fragment"),
    "longest_repeat": ("Longest repeat", "bases of the longest repeat over 11 bases (0: none)"),
    "repeat_9_metric": ("Repeat load", "9 for every base beyond the 8th of each copy of each repeat of 9+ bases, per base of the fragment"),
    "richest_hairpin": ("Richest hairpin", "GC fraction of the most GC-rich hairpin stem found"),
    "GC_short_l": ("Low-GC windows", "20-base windows with 20% GC or less"),
    "GC_term_l": ("Low-GC ends", "20-base windows with under 30% GC among the first and last 50 bases"),
    "Tm_low": ("Cold windows", "20-base windows melting at 40 °C or below"),
    "Tm_high": ("Hot windows", "20-base windows melting at 70 °C or above"),
}


def gc_fraction(seq, window):
    """seq_eval.GC_count: (G + C) / window, 0 when there is no G or C."""
    g, c = seq.count("G"), seq.count("C")
    return float(g + c) / window if g or c else 0


# DNA_NN3 of Biopython (Allawi & SantaLucia 1997): (dH kcal/mol, dS cal/mol/K)
_NN = {
    "AA/TT": (-7.9, -22.2), "AT/TA": (-7.2, -20.4), "TA/AT": (-7.2, -21.3), "CA/GT": (-8.5, -22.7), "GT/CA": (-8.4, -22.4),
    "CT/GA": (-7.8, -21.0), "GA/CT": (-8.2, -22.2), "CG/GC": (-10.6, -27.2), "GC/CG": (-9.8, -24.4), "GG/CC": (-8.0, -19.9),
}
_INIT_AT, _INIT_GC = (2.3, 4.1), (0.1, -2.8)
_LOG_MON = math.log(50 * 1e-3)   # Na = 50 mM, no Mg: Biopython's salt correction 7 takes its Mg-free branch


def tm_nn(seq):
    """Biopython Tm_NN(seq, dnac1=250.0, saltcorr=7) for a perfectly paired A/C/G/T oligo."""
    ends = seq[0] + seq[-1]
    at = ends.count("A") + ends.count("T")
    gc = ends.count("G") + ends.count("C")
    dh = _INIT_AT[0] * at + _INIT_GC[0] * gc
    ds = _INIT_AT[1] * at + _INIT_GC[1] * gc
    comp = complement(seq)
    for i in range(len(seq) - 1):
        a, b = seq[i:i + 2], comp[i:i + 2]
        key = a + "/" + b
        h, s = _NN[key] if key in _NN else _NN[key[::-1]]
        dh += h
        ds += s
    k = (250.0 - (25 / 2.0)) * 1e-9
    tm = (1000 * dh) / (ds + (1.987 * math.log(k))) - 273.15
    gcf = (seq.count("G") + seq.count("C")) / len(seq)
    corr = (4.29 * gcf - 3.95) * 1e-5 * _LOG_MON + 9.40e-6 * _LOG_MON ** 2
    return 1 / (1 / (tm + 273.15) + corr) - 273.15


# ---- the nine features ----------------------------------------------------------------------------------

def _repeat_features(seq, places=None):
    """longest_repeat (longest repeat over 11 bases, else 0) and repeat_9_metric (9 per base beyond the 8th,
    for every copy of every repeat of 9+ bases, over the size). `places` (a dict), when given, gets where each is."""
    longest, total = 0, 0
    repeats = opendna_core.salis_repeats(seq, 9)
    for length, _inverted, starts in repeats:
        total += 9 * (length - 8) * len(starts)
        if length > 11 and length >= longest:
            longest = length
    if places is not None:
        places["repeat_9_metric"] = merged((s, s + length) for length, _i, starts in repeats for s in starts)
        places["longest_repeat"] = merged((s, s + length) for length, _i, starts in repeats for s in starts if length == longest and longest)
    return longest, total / float(len(seq))


def _window_features(seq, places=None):
    """GC_short_l, GC_term_l, Tm_low, Tm_high, total_GC. `places` (a dict), when given, gets the windows each counts."""
    size = len(seq)
    short_l = tm_low = tm_high = term_l = 0
    at = {"GC_short_l": [], "GC_term_l": [], "Tm_low": [], "Tm_high": []}   # the 20-base windows counted, by feature
    for i in range(size - 20):
        w = seq[i:i + 20]
        if gc_fraction(w, 20) <= 0.20:
            short_l += 1
            at["GC_short_l"].append(i)
        tm = tm_nn(w)
        if tm <= 40:
            tm_low += 1
            at["Tm_low"].append(i)
        if tm >= 70:
            tm_high += 1
            at["Tm_high"].append(i)
    if size > 50:
        for i in range(30):
            if gc_fraction(seq[i:i + 20], 20) < 0.30:
                term_l += 1
                at["GC_term_l"].append(i)
            if gc_fraction(seq[size - 50 + i:size - 50 + i + 20], 20) < 0.30:
                term_l += 1
                at["GC_term_l"].append(size - 50 + i)
    else:
        for i in range(size - 20):
            f = gc_fraction(seq[i:i + 20], 20)
            term_l += 2 * (f < 0.30)   # five- and three-prime windows are the same window here
            if f < 0.30:
                at["GC_term_l"].append(i)
    if places is not None:
        places.update({name: merged((i, i + 20) for i in starts) for name, starts in at.items()})
        places["total_GC"] = places["size"] = []   # (the whole sequence: nothing to point at)
    return short_l, term_l, tm_low, tm_high, gc_fraction(seq, size)


def _fold_details(structure, energy, seq):
    """PyVRNA.vienna2bp + seq_eval.get_hairpin_lengths: (helix sizes, (hairpin text, energy) list, G-quadruplex texts)."""
    stack, bp = [], {}
    has_gquad = False
    for index, ch in enumerate(structure):
        if ch == "(":
            stack.append(index + 1)
            bp[index + 1] = -1
        elif ch == ")":
            bp[stack.pop()] = index + 1
        elif ch == "+":
            has_gquad = True
    bp_x, bp_y = list(bp), list(bp.values())
    lengths, details = [], []
    while bp_x:
        n = bisect(bp_x, bp_y[0])
        lengths.append(n)
        details.append((seq[bp_x[0] - 1:bp_y[0]], energy, bp_x[0] - 1))   # (the third: where the hairpin starts in `seq`)
        bp_x, bp_y = bp_x[n:], bp_y[n:]
    gquads = []
    if has_gquad:
        rows, layer, i = [], 0, 0
        while i < len(structure):
            if structure[i] == "+" and layer in (0, 2):
                if layer == 0:
                    rows.append([deque(), deque(), deque(), deque()])
                while i < len(structure) and structure[i] == "+":
                    rows[-1][layer].append(i + 1)
                    i += 1
                layer = (layer + 1) % 4
            elif structure[i] == "+":
                while i < len(structure) and structure[i] == "+":
                    rows[-1][layer].appendleft(i + 1)
                    i += 1
                layer = (layer + 1) % 4
            else:
                i += 1
        tuples = []
        for row in rows:
            tuples.extend(zip(*row))
        for t in range(len(tuples) // 3):
            first, last = tuples[3 * t][0], tuples[3 * t][3]
            gquads.append(seq[first - 1:last])
    return lengths, details, gquads


def _vienna_richest(seq, progress=None, cancelled=None, places=None):
    """The highest GC fraction of a folded 'strong hairpin', as seq_eval.find_hairpins_windowed records it.
    The folds of every 100-base window (both strands, spans 50, 100 and 50) also fill the table of hairpin texts
    already seen, and a text seen once is not counted again, so all three spans are walked as the original does.
    `places` (a dict), when given, gets "richest_hairpin": [(start, end)] of the first hairpin with that GC, on the top strand."""
    import RNA
    md = RNA.md()
    md.temperature, md.dangles, md.gquad, md.noGU, md.noLP = 37.0, 0, 1, 1, 0
    for name, value in (("temperature", 37.0), ("dangles", 0), ("gquad", 1), ("noGU", 1), ("noLP", 0)):
        setattr(RNA.cvar, name, value)
    RNA.params_load_DNA_Mathews2004()
    size = len(seq)
    reverse = reverse_complement(seq)
    seen, richest, where = set(), 0, []
    spans, thresholds, mfes = (50, 100, 50), (11, 17, 21), (-10, -15, -20)

    def fold(window, span):
        md.max_bp_span = span
        structure, energy = RNA.fold_compound(window, md).mfe()
        return _fold_details(structure, energy, window)

    total = max(size - 100, 0)
    for j in range(total):
        if cancelled and cancelled():
            raise Cancelled()
        if progress:
            progress(j, total)
        windows = (seq[j:j + 100], reverse[size - (j + 100):size - j])
        cache = {}
        for span, lt, mfe in zip(spans, thresholds, mfes):
            for window in windows:
                if (window, span) not in cache:
                    cache[(window, span)] = fold(window, span)
                lengths, details, gquads = cache[(window, span)]
                for n, (hairpin, energy, at) in zip(lengths, details):
                    if not (n > lt and energy < mfe):
                        continue
                    if lt < 16:
                        gc = gc_fraction(hairpin, len(hairpin))
                        if gc > 0.8 and hairpin not in seen:
                            seen.add(hairpin)
                            if gc > richest:   # (the top strand: a hairpin of the reverse window lies mirrored in seq[j:j + 100])
                                where = [(j + at, j + at + len(hairpin)) if window is windows[0] else (j + 100 - at - len(hairpin), j + 100 - at)]
                            richest = max(richest, gc)
                    else:
                        seen.add(hairpin)
                seen.update(gquads)
    if progress:
        progress(total, total)
    if places is not None:
        places["richest_hairpin"] = where
    return richest


def merged(spans):
    """(start, end) spans sorted, overlapping or touching ones joined."""
    out = []
    for a, b in sorted(spans):
        if out and a <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], b))
        else:
            out.append((a, b))
    return out


def features(seq, progress=None, cancelled=None, places=None):
    """The nine numbers of the model for `seq`. Needs ViennaRNA (see have_vienna). `places` (a dict), when given, gets for each
    number the (start, end) stretches of `seq` it counted (merged; none for total_GC and size), for the window to point at."""
    seq = clean(seq)
    if not have_vienna():
        raise UserError(VIENNA_HELP)
    longest, repeat9 = _repeat_features(seq, places)
    short_l, term_l, tm_low, tm_high, total_gc = _window_features(seq, places)
    folded = {}
    vienna = _vienna_richest(seq, progress, cancelled, folded)
    rule, rule_at = opendna_core.salis_hairpin_gc(seq)
    richest = max(vienna, rule)
    if places is not None:   # (where the richest one is: the fold's, unless the rule-based one is richer)
        places["richest_hairpin"] = folded["richest_hairpin"] if vienna >= rule else ([rule_at] if rule_at else [])
    return {"longest_repeat": longest, "repeat_9_metric": repeat9, "richest_hairpin": richest, "GC_short_l": short_l,
            "GC_term_l": term_l, "Tm_low": tm_low, "Tm_high": tm_high, "total_GC": total_gc, "size": len(seq)}


def predict(values):
    """(probability of 'Synthesized', outcome) from the nine numbers: the mean of the forest's 1512 trees, as
    scikit-learn 0.20.4 computes it (inputs rounded to 32-bit floats first, a split goes left when <=)."""
    model = _model()
    x = [struct.unpack("f", struct.pack("f", float(values[name])))[0] for name in model["features"]]
    s0 = s1 = 0.0
    for t in model["trees"]:
        left, right, feature, threshold = t["left"], t["right"], t["feature"], t["threshold"]
        n = 0
        while left[n] != -1:
            n = left[n] if x[feature[n]] <= threshold[n] else right[n]
        s0 += t["p0"][n]
        s1 += t["p1"][n]
    count = len(model["trees"])
    p0, p1 = s0 / count, s1 / count
    return p1, p1 > p0


def evaluate(seq, progress=None, cancelled=None):
    """{'features': ..., 'probability': ..., 'outcome': True if the fragment is predicted to be synthesised quickly}."""
    places = {}
    values = features(seq, progress, cancelled, places)
    probability, outcome = predict(values)
    return {"features": values, "probability": probability, "outcome": outcome, "places": places}
