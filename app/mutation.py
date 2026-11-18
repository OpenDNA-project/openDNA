"""Mutations written in HGVS (https://hgvs-nomenclature.org), turned into changes of a plasmid. No windows.

    lacZ:c.360G>C  lacZ:c.360_365del  lacZ:c.360_361insATG  lacZ:c.360_365delinsATG  lacZ:c.360_365dup    (DNA, numbered from the A of the start codon)
    o.1234G>C  o.1234_1240del  ...                                                                      (DNA, numbered along the plasmid; `g.` is accepted)
    lacZ:p.Asp120Ala  lacZ:D120A  lacZ:p.Lys12_Gly14del  lacZ:p.Lys12_Gly13insAla  lacZ:p.Gly13dup      (protein; one or three letters)

`resolve` gives the change as a replacement of plasmid[start:end] (0-based, forward strand), with the HGVS text of what it does and
its effect on the proteins. A mistake the sender can fix raises UserError (its text is shown as it is); anything that could be read in
more than one way (several codons for an amino acid, a repeat the change could slide along) gives a warning and a default, never a stop.
"""
import re
from types import SimpleNamespace

import genetic_code
from editor import complement, reverse_complement
from user_error import UserError

THREE = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
         "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
         "Ter": "*"}
LETTER3 = {v: k for k, v in THREE.items()}
AA = "(?:" + "|".join(THREE) + "|[A-Z*])"
POS = r"-?\d+"
DNA_FORMS = [   # (kind, pattern); delins before del
    ("sub", rf"(?P<a>{POS})(?P<ref>[ACGT])>(?P<alt>[ACGT])"),
    ("delins", rf"(?P<a>{POS})(?:_(?P<b>{POS}))?delins(?P<alt>[ACGT]+)"),
    ("del", rf"(?P<a>{POS})(?:_(?P<b>{POS}))?del(?P<ref>[ACGT]*)"),
    ("ins", rf"(?P<a>{POS})_(?P<b>{POS})ins(?P<alt>[ACGT]+)"),
    ("dup", rf"(?P<a>{POS})(?:_(?P<b>{POS}))?dup(?P<ref>[ACGT]*)"),
]
PROTEIN_FORMS = [
    ("sub", rf"(?P<ra>{AA})(?P<a>\d+)(?P<alt>{AA})"),
    ("delins", rf"(?P<ra>{AA})(?P<a>\d+)(?:_(?P<rb>{AA})(?P<b>\d+))?delins(?P<alt>(?:{AA})+)"),
    ("del", rf"(?P<ra>{AA})(?P<a>\d+)(?:_(?P<rb>{AA})(?P<b>\d+))?del"),
    ("ins", rf"(?P<ra>{AA})(?P<a>\d+)_(?P<rb>{AA})(?P<b>\d+)ins(?P<alt>(?:{AA})+)"),
    ("dup", rf"(?P<ra>{AA})(?P<a>\d+)(?:_(?P<rb>{AA})(?P<b>\d+))?dup"),
]
ACCEPTED = ("Accepted forms (1-based): lacZ:c.360G>C, lacZ:c.360_365del, lacZ:c.360_361insATG, lacZ:c.360_365delinsATG, lacZ:c.360_365dup; "
            "o.1234G>C and the same with o. positions along the plasmid; lacZ:p.Asp120Ala (or lacZ:D120A), lacZ:p.Lys12_Gly14del, "
            "lacZ:p.Lys12_Gly13insAla, lacZ:p.Lys12_Gly14delinsAla, lacZ:p.Gly13dup.")

def hull(f):
    """(first, last) base of a feature, 1-based inclusive, over all its parts."""
    return min(s for s, _ in f.spans) + 1, max(e for _, e in f.spans)


def find_feature(doc, name):
    """The one feature called `name` (its label, case ignored, or #INDEX). A gene and its CDS share a label: the CDS is the one meant."""
    if name.startswith("#") and name[1:].isdigit():
        found = [f for f in doc.features() if f.index == int(name[1:])]
    else:
        found = [f for f in doc.features() if f.label.lower() == name.lower() and f.kind != "source"]
    if not found:
        raise UserError(f"No feature called {name!r}. Run `inspect` for the names.")
    if len(found) > 1:
        cds = [f for f in found if f.kind == "CDS"]
        if len(cds) == 1:
            return cds[0]
        raise UserError(f"{len(found)} features are called {name!r}: " + ", ".join(f"#{f.index} {f.kind} {hull(f)[0]}..{hull(f)[1]}" for f in found)
                        + ". Use #INDEX.")
    return found[0]


def coding_positions(f):
    """Plasmid index (0-based) of every base of a feature, in the direction of the feature."""
    idx = [i for s, e in f.spans for i in range(s, e)]
    return idx[::-1] if f.strand < 0 else idx


def parse(spec):
    """(feature name or None, scope 'c' | 'o' | 'p', kind, groups) of one mutation text."""
    text = spec.strip()
    name, _, body = text.partition(":") if ":" in text else (None, "", text)
    name = name or None
    scope = "c"
    m = re.match(r"([cogp])\.(.*)$", body)
    if m:
        scope, body = ("o" if m.group(1) == "g" else m.group(1)), m.group(2)
    elif re.fullmatch(rf"{AA}\d+{AA}", body):
        scope = "p"
    else:
        raise UserError(f"Cannot read the mutation {spec!r}. {ACCEPTED}")
    if scope in "co" and "*" in body:
        raise UserError(f"{spec!r}: positions after the stop codon (c.*N) are not supported; use o. positions along the plasmid.")
    for kind, pattern in (PROTEIN_FORMS if scope == "p" else DNA_FORMS):
        m = re.fullmatch(pattern, body)
        if m:
            return name, scope, kind, m.groupdict()
    raise UserError(f"Cannot read the mutation {spec!r}. {ACCEPTED}")


def one(a):
    """One-letter amino acid from a one- or three-letter one."""
    return THREE.get(a, a)


def letters(text):
    return [one(a) for a in re.findall(AA, text)]


def resolve(doc, spec, usage=None):
    """The mutation `spec` on `doc`: SimpleNamespace(spec, feature, start, end, text, warnings, hgvs_dna, hgvs_protein, effect),
    the change being doc[start:end] -> text (0-based, forward strand). `usage` is a codon usage table (`.per_thousand`) or None."""
    name, scope, kind, g = parse(spec)
    n = doc.length
    warnings = []
    seq = doc.sequence(0, n).upper()
    f = None
    if scope != "o":   # (o. numbers run along the plasmid: a name before it is not needed)
        if not name:
            raise UserError(f"{spec!r}: name the feature, as in lacZ:{spec}")
        f = find_feature(doc, name)
    if scope == "p":
        if f.kind != "CDS":
            raise UserError(f"{spec!r}: {name!r} is a {f.kind}, not a CDS: protein positions need a CDS.")
        start, end, text = _protein_change(doc, f, seq, kind, g, usage, warnings, spec)
    else:
        start, end, text = _dna_change(f if scope == "c" else None, seq, kind, g, warnings, spec)
    asked = (start, end)
    while start < end and text and seq[start] == text[0]:   # a replacement keeps only the bases that differ (a codon change of one base is a substitution)
        start, text = start + 1, text[1:]
    while start < end and text and seq[end - 1] == text[-1]:
        end, text = end - 1, text[:-1]
    if start == end and not text:   # the plasmid already has what was asked: accepted, with a warning, and nothing to make
        warnings.append(f"{spec}: no change, the plasmid already has these bases ({seq[asked[0]:asked[1]] or 'nothing to add'}): no primers are designed.")
        r = SimpleNamespace(spec=spec, feature=None if f is None else f.label, start=start, end=end, text=text, warnings=warnings)
        r.hgvs_dna = f"o.{asked[0] + 1}{'' if asked[1] - asked[0] < 2 else f'_{asked[1]}'}="
        r.hgvs_protein, r.effect = "p.(=)", "no change"
        return r
    if start == end and start in (0, n):
        raise UserError(f"{spec!r}: an insertion at the very start or end touches the origin of the plasmid, which is not supported.")
    slide = sliding(seq, start, end, text)
    if slide:
        warnings.append(f"{spec}: in a repeat, the same plasmid results if the change is placed up to {slide} base{'s' if slide > 1 else ''} further on "
                        "(HGVS writes it at the most 3' place; the primers differ slightly with the place).")
    r = SimpleNamespace(spec=spec, feature=None if f is None else f.label, start=start, end=end, text=text, warnings=warnings)
    r.hgvs_dna = describe(seq, start, end, text)
    r.hgvs_protein, r.effect = effect(doc, seq, start, end, text, warnings)
    return r


def _map(f, pos, spec):
    """Plasmid index of c.`pos` (1 = first base of the feature, -1 = the one before it)."""
    idx = coding_positions(f)
    if pos == 0:
        raise UserError(f"{spec!r}: there is no position 0 (c.1 is the first base, c.-1 the one before it).")
    if pos > 0:
        if pos > len(idx):
            raise UserError(f"{spec!r}: position {pos} is past the end of {f.label!r} ({len(idx)} bases).")
        return idx[pos - 1]
    return idx[0] + pos if f.strand > 0 else idx[0] - pos


def _dna_change(f, seq, kind, g, warnings, spec):
    n = len(seq)
    a = int(g["a"])
    b = int(g["b"]) if g.get("b") else a
    if b < a:
        raise UserError(f"{spec!r}: the first position is after the last one.")
    if f is None:   # o.: along the plasmid
        for p in (a, b):
            if not 1 <= p <= n:
                raise UserError(f"{spec!r}: position {p} is outside the plasmid (1..{n}).")
        ia, ib, strand = a - 1, b - 1, 1
    else:
        ia, ib, strand = _map(f, a, spec), _map(f, b, spec), f.strand
        for p in (ia, ib):
            if not 0 <= p < n:
                raise UserError(f"{spec!r}: that position is off the end of the plasmid (the origin is not crossed).")
    out = reverse_complement if strand < 0 else str   # the bases as written -> forward strand
    if kind == "ins":
        if abs(ia - ib) != 1:
            raise UserError(f"{spec!r}: an insertion goes between two neighbouring bases (c.360_361insATG).")
        at = max(ia, ib)
        return at, at, out(g["alt"])
    lo, hi = min(ia, ib), max(ia, ib) + 1
    if a * b > 0 and hi - lo != b - a + 1:
        raise UserError(f"{spec!r}: that range crosses the origin or a gap of the feature: write it with o. positions.")
    old = seq[lo:hi]
    expected = g.get("ref") or ""
    if expected and out(expected) != old:   # the bases the sender expects must be the ones there
        raise UserError(f"{spec!r}: the plasmid has {old if strand > 0 else reverse_complement(old)} at that place, not {expected}.")
    if kind in ("sub", "delins"):
        return lo, hi, out(g["alt"])
    if kind == "del":
        return lo, hi, ""
    return hi, hi, old   # dup: a second copy straight after


def _codons(doc, f, spec):
    """(coding-order plasmid indices, codon table, the protein) of a CDS read as annotated, or UserError."""
    plan = genetic_code.cds_plan(doc, f)
    if plan.why or plan.exceptions:
        raise UserError(f"{spec!r}: {f.label!r} has an annotation openDNA cannot follow ({plan.why or '/transl_except'}): use c. or o. positions.")
    if plan.frame:
        raise UserError(f"{spec!r}: {f.label!r} has /codon_start={plan.frame + 1}: use c. or o. positions.")
    idx = coding_positions(f)
    code = genetic_code.codon_map(plan.table)
    seq = doc.sequence(0, doc.length).upper()
    cod = [("".join(seq[i] for i in idx[k:k + 3])) for k in range(0, len(idx) - 2, 3)]
    if f.strand < 0:   # (idx runs against the plasmid: the bases are read complemented)
        cod = [complement(c) for c in cod]
    return idx, code, cod


def _best(code, usage, aa, current=None):
    """The codon for `aa`: the fewest base changes from `current` (if given), then the most used; and the others equally near."""
    options = sorted(c for c, a in code.items() if a == aa)
    if not options:
        raise UserError(f"No codon for {aa!r} in the genetic code in use.")
    dist = (lambda c: sum(x != y for x, y in zip(c, current))) if current else (lambda c: 0)
    use = (lambda c: -usage.per_thousand[c]) if usage is not None else (lambda c: 0)
    options.sort(key=lambda c: (dist(c), use(c), c))
    best = options[0]
    return best, [c for c in options[1:] if dist(c) == dist(best)]


def _protein_change(doc, f, seq, kind, g, usage, warnings, spec):
    idx, code, cod = _codons(doc, f, spec)
    protein = [code.get(c, "X") for c in cod]
    a = int(g["a"])
    b = int(g["b"]) if g.get("b") else a
    for p, letter in ((a, g["ra"]), (b, g.get("rb") or g["ra"])):
        if not 1 <= p <= len(protein):
            raise UserError(f"{spec!r}: residue {p} is outside {f.label!r} (1..{len(protein)}).")
        if protein[p - 1] != one(letter):
            raise UserError(f"{spec!r}: residue {p} of {f.label!r} is {LETTER3.get(protein[p - 1], protein[p - 1])} ({cod[p - 1]}), not {letter}.")
    if b < a:
        raise UserError(f"{spec!r}: the first residue is after the last one.")
    span = idx[3 * (a - 1):3 * b]
    lo, hi = min(span), max(span) + 1
    out = reverse_complement if f.strand < 0 else str
    if kind == "sub":
        new = one(g["alt"])
        if protein[a - 1] == new:
            raise UserError(f"{spec!r}: residue {a} is already {g['alt']}: nothing to change.")
        codon, others = _best(code, usage, new, cod[a - 1])
        if others:
            warnings.append(f"{spec}: {codon} chosen for {LETTER3.get(new, new)}; as few base changes: {', '.join(others)}"
                            + (" (ties broken by the most used codon)" if usage is not None else " (no codon usage table: alphabetical)"))
        return lo, hi, out(codon)
    if kind == "del":
        return lo, hi, ""
    if kind == "dup":
        return hi, hi, seq[lo:hi]
    new = letters(g["alt"])
    chosen = [_best(code, usage, aa) for aa in new]
    note = ", ".join(f"{LETTER3.get(aa, aa)} {c}" for aa, (c, _) in zip(new, chosen))
    warnings.append(f"{spec}: codons chosen for the new residues (the most used in the codon usage table, else alphabetical): {note}")
    text = out("".join(c for c, _ in chosen))
    if kind == "delins":
        return lo, hi, text
    if b != a + 1:
        raise UserError(f"{spec!r}: an insertion goes between two neighbouring residues (p.Lys12_Gly13insAla).")
    # (kind == "ins": a and b are the neighbours; `span` covers both residues)
    mid = idx[3 * a - 1]   # last base of residue a in coding order
    at = mid + 1 if f.strand > 0 else mid
    return at, at, text


def sliding(seq, start, end, text):
    """How many bases a deletion or insertion can move right and give the same plasmid (0 when it cannot)."""
    k = 0
    if start < end and not text:
        while end + k < len(seq) and seq[start + k] == seq[end + k]:
            k += 1
    elif start == end and text:
        while start + k < len(seq) and seq[start + k] == text[k % len(text)]:
            k += 1
    return k


def describe(seq, start, end, text):
    """The change as HGVS on the plasmid: o.1234G>C, o.1234_1240del, o.1234_1235insATG, o.1234_1240delinsATG."""
    old = seq[start:end]
    if start == end:
        return f"o.{start}_{start + 1}ins{text}"
    if not text:
        return f"o.{start + 1}del" if end - start == 1 else f"o.{start + 1}_{end}del"
    if len(old) == len(text) == 1:
        return f"o.{start + 1}{old}>{text}"
    return f"o.{start + 1}{'' if end - start == 1 else f'_{end}'}delins{text}"


def translate(dna, code):
    return "".join(code.get(dna[i:i + 3], "X") for i in range(0, len(dna) - len(dna) % 3, 3))


def t3(aa):
    return LETTER3.get(aa, aa)


def to_stop(protein):
    """`protein` up to and including its first stop."""
    return protein[:protein.index("*") + 1] if "*" in protein else protein


def residues(aas, first):
    """HGVS text of the amino acids `aas` that start at residue number `first`: Lys12, or Lys12_Gly14."""
    return f"{t3(aas[0])}{first}" + ("" if len(aas) == 1 else f"_{t3(aas[-1])}{first + len(aas) - 1}")


def effect(doc, seq, start, end, text, warnings):
    """(HGVS protein text of the CDS the change falls in, effect word); ("", "non-coding") outside every CDS."""
    delta = len(text) - (end - start)
    cds = [f for f in doc.features() if f.kind == "CDS"]
    hits = [f for f in cds if len(f.spans) == 1 and f.spans[0][0] <= start and end <= f.spans[0][1] and not (start == end and start in f.spans[0])]
    if not hits:
        touched = [f for f in cds if any(s < end and start < e for s, e in f.spans)]
        if touched:
            warnings.append(f"The change overlaps the edge of {touched[0].label!r}: its protein effect is not worked out.")
            return "", "overlaps CDS edge"
        return "", "non-coding"
    f = hits[0]
    s, e = f.spans[0]
    plan = genetic_code.cds_plan(doc, f)
    code = genetic_code.codon_map(plan.table)
    old_cds, new_cds = seq[s:e], seq[s:start] + text + seq[end:e]   # (only the CDS: the rest of the plasmid does not matter)
    if f.strand < 0:
        old_cds, new_cds = reverse_complement(old_cds), reverse_complement(new_cds)
    old_p, new_p = to_stop(translate(old_cds[plan.frame:], code)), to_stop(translate(new_cds[plan.frame:], code))
    p = 0
    while p < min(len(old_p), len(new_p)) and old_p[p] == new_p[p]:
        p += 1
    if old_p == new_p:
        return "p.(=)", "silent"
    if delta % 3:
        warnings.append(f"Frameshift in {f.label!r} from residue {p + 1}.")
        return f"p.({t3(old_p[p]) if p < len(old_p) else 'Ter'}{p + 1}fs)", "frameshift"
    if p < len(new_p) and new_p[p] == "*" and p < len(old_p) and old_p[p] != "*":   # a stop where there was an amino acid
        warnings.append(f"A stop codon appears in {f.label!r} at residue {p + 1}: the protein ends there.")
        return f"p.({t3(old_p[p])}{p + 1}Ter)", "nonsense"
    o_mid_end, n_mid_end = len(old_p), len(new_p)
    while o_mid_end > p and n_mid_end > p and old_p[o_mid_end - 1] == new_p[n_mid_end - 1]:
        o_mid_end, n_mid_end = o_mid_end - 1, n_mid_end - 1
    old_mid, new_mid = old_p[p:o_mid_end], new_p[p:n_mid_end]
    if p == 0 and old_p[0] == "M":
        warnings.append(f"The start codon of {f.label!r} is changed.")
    if len(old_mid) == len(new_mid) == 1:
        return f"p.({t3(old_mid)}{p + 1}{t3(new_mid)})", "missense"
    if "*" in new_mid[:-1] or (new_mid.endswith("*") and not old_mid.endswith("*")):
        warnings.append(f"A premature stop codon appears in {f.label!r}.")
    if not new_mid:
        return f"p.({residues(old_mid, p + 1)}del)", "in-frame deletion"
    if not old_mid and p >= len(new_mid) and old_p[p - len(new_mid):p] == new_mid:
        return f"p.({residues(new_mid, p - len(new_mid) + 1)}dup)", "in-frame duplication"
    if not old_mid:
        return f"p.({t3(old_p[p - 1])}{p}_{t3(old_p[p])}{p + 1}ins{''.join(map(t3, new_mid))})", "in-frame insertion"
    return f"p.({residues(old_mid, p + 1)}delins{''.join(map(t3, new_mid))})", "in-frame deletion-insertion"
