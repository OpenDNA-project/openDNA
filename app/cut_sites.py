"""Cut sites in words and in numbers: what a restriction enzyme or a Cas nuclease does to a sequence, shared by the sequence
view (drawing), the hover tooltips and the Spacer Finder. No GUI code here.

Two parts:
* Restriction enzymes: `enzyme_tip` says where an enzyme cuts each strand (1-based, "after base N"), the ends it leaves and,
  for a Type IIS enzyme, how far from its site.
* Cas9 / Cas12a guide sites: `NUCLEASES` (the nucleases the Spacer Finder knows), `cut_places`, the way a saved guide is
  written in a GenBank feature (`cas9_qualifier`, read back by `parse_cas9`), `cas9_layout` (where the protospacer, the PAM and
  the cut are, which the view draws) and `cas9_tip`.

How a saved guide site is stored (so any GenBank reader accepts it): a `misc_feature` over protospacer + PAM, on the strand
the guide reads along (complement(...) for the bottom strand), with
    /label="gRNA GACTCGTA"
    /note="SpCas9 guide GACTCGTACGTAGCATCGAT, PAM TGG; blunt cut 3 bp from the PAM, after base 57"   (plain words)
    /OpenDNA_cas9="name=SpCas9;pam=NGG;length=20;cut=3"                                          (for this program)
Other programs show an ordinary annotated feature (and ignore the unknown qualifier, like the /ApEinfo_ colour ones); this
program redraws the protospacer, PAM and cut from the qualifier. The guide and the PAM are not stored a second time: they are
the bases under the feature (read along the feature's strand), so editing the sequence keeps them right.
"""
from collections import namedtuple
from types import SimpleNamespace

from enzymes import notation, revcomp_site, site_spans

# ---- Cas nucleases ----------------------------------------------------------------------------------------------------

# name -> (PAM as IUPAC letters, protospacer length, then for Cas9: the bases between the cut and the PAM, and for Cas12a, whose
# PAM is on the 5' side and which cuts the two strands in different places: the bases from the PAM to the cut of the PAM's strand and
# of the other strand)
Nuclease = namedtuple("Nuclease", "name pam length cut_from_pam pam_first cuts")
NUCLEASES = {
    "SpCas9 (NGG, 20 nt)": Nuclease("SpCas9", "NGG", 20, 3, False, None),
    "SaCas9 (NNGRRT, 21 nt)": Nuclease("SaCas9", "NNGRRT", 21, 3, False, None),
    "Cas12a (TTTV, 23 nt)": Nuclease("Cas12a", "TTTV", 23, None, True, (18, 23)),
}
CAS9_KEY = "OpenDNA_cas9"   # the program-specific qualifier (like ApEinfo_fwdcolor): ignored by other readers
CAS9_KEYS = (CAS9_KEY, "OpenGenEdit_cas9")   # read both: files saved before the rename to openDNA carry the old key


def cut_places(nuclease, start, strand):
    """(top-strand cut, bottom-strand cut) of the site that starts at `start` (on the top strand): each the gap between two bases,
    counted as the base after it. Cas9 cuts both strands in the same place, `cut_from_pam` bases inside the protospacer; Cas12a
    cuts the strand with the PAM `cuts[0]` bases after the PAM and the other strand `cuts[1]` bases after it (a 5' overhang)."""
    L, P = nuclease.length, len(nuclease.pam)
    if nuclease.pam_first:
        near, far = nuclease.cuts
        # PAM strand = the strand the guide reads along. Forward: the PAM is on the left of the site, the PAM strand is the top.
        return (start + P + near, start + P + far) if strand > 0 else (start + L - far, start + L - near)
    cut = start + L - nuclease.cut_from_pam if strand > 0 else start + P + nuclease.cut_from_pam
    return cut, cut


def cas9_qualifier(nuclease):
    """The value of the /OpenDNA_cas9 qualifier: what is needed to know the PAM, protospacer and cut from the feature's
    location. `key=value` pairs, no quotes or line breaks, so it is valid GenBank text as it is."""
    parts = [f"name={nuclease.name}", f"pam={nuclease.pam}", f"length={nuclease.length}"]
    if nuclease.pam_first:
        parts += ["pam_first=yes", "cuts=%d,%d" % tuple(nuclease.cuts)]
    else:
        parts.append(f"cut={nuclease.cut_from_pam}")
    return ";".join(parts)


def parse_cas9(value):
    """The Nuclease a /OpenDNA_cas9 value describes, or None if it is not a usable one (damaged by hand)."""
    try:
        d = dict(part.split("=", 1) for part in (value or "").split(";"))
        pam_first = d.get("pam_first") == "yes"
        cuts = tuple(int(x) for x in d["cuts"].split(",")) if pam_first else None
        if pam_first and len(cuts) != 2:
            return None
        return Nuclease(d["name"], d["pam"], int(d["length"]), None if pam_first else int(d["cut"]), pam_first, cuts)
    except (KeyError, ValueError):
        return None


def cas9_of(feature):
    """The Nuclease of a feature saved as a guide site, or None for any other feature."""
    for key, value in feature.qualifiers:
        if key in CAS9_KEYS:
            return parse_cas9(value)
    return None


def runs(cells):
    """[(start, end)] of the consecutive stretches of a list of positions."""
    out = []
    for c in cells:
        if out and out[-1][1] == c:
            out[-1][1] = c + 1
        else:
            out.append([c, c + 1])
    return [tuple(r) for r in out]


def site_cells(feature, length):
    """The bases of a feature in top-strand order, left to right (a site across the origin of a circle: the end, then the start)."""
    spans = list(feature.spans)
    if [s for s, _ in spans] != sorted(s for s, _ in spans) and len(spans) == 2:
        spans.sort(key=lambda span: -span[0])
    return [c for s, e in spans for c in range(s, e)]


def cas9_layout(nuclease, cells, strand, length):
    """Where the parts of a guide site are. `cells`: the positions of its bases, left to right on the top strand (protospacer and
    PAM together); `strand`: the one the guide reads along. Returns a namespace: `proto` and `pam` (lists of (start, end) runs),
    `proto_cells`, `pam_cells`, `gaps` (the places the nuclease cuts, as the gap before that base: one for Cas9, two for Cas12a,
    top strand's first), `top_gap`, `bottom_gap`. None if `cells` is not the length of a site."""
    L, P = nuclease.length, len(nuclease.pam)
    if len(cells) != L + P:
        return None
    # The PAM is 3' of the protospacer on the guide's strand (5' of it for Cas12a); on the bottom strand that is the left.
    pam_left = nuclease.pam_first == (strand > 0)
    pam, proto = (cells[:P], cells[P:]) if pam_left else (cells[L:], cells[:L])

    def gap(offset):   # the gap before cells[offset]; one past the last base is the end of the site
        return cells[offset - 1] + 1 if offset > 0 else cells[0]   # (== length at the end of a circle: the origin)

    top, bottom = cut_places(nuclease, 0, strand)
    gaps = [gap(top)] if top == bottom else [gap(top), gap(bottom)]
    return SimpleNamespace(nuclease=nuclease, strand=strand, cells=cells, proto_cells=proto, pam_cells=pam,
                           proto=runs(proto), pam=runs(pam), gaps=gaps, top_gap=gap(top), bottom_gap=gap(bottom))


def feature_layout(feature, length):
    """cas9_layout of a feature saved as a guide site (None for another feature, or one whose size no longer fits)."""
    nuclease = cas9_of(feature)
    if nuclease is None or length == 0:
        return None
    return cas9_layout(nuclease, site_cells(feature, length), feature.strand if feature.strand else 1, length)


def _word(strand):
    return "top (forward)" if strand > 0 else "bottom (reverse)"


def _after(q, n, circular):
    """'after base N' for a cut at gap `q` (the cut falls before the base at index q), bases numbered from 1."""
    base = (q - 1) % n + 1 if circular or q > 0 else 0
    return f"after base {base:,}"


def cas9_tip(layout, fetch, length, label=None):
    """Tooltip (rich text) for a guide site. `fetch(a, b)` returns the bases [a, b) of the sequence."""
    nuc = layout.nuclease
    read = lambda cells: "".join(fetch(a, b) for a, b in runs(cells))
    guide, pam = read(layout.proto_cells), read(layout.pam_cells)
    if layout.strand < 0:   # the guide reads along the bottom strand
        comp = str.maketrans("ACGTURYSWKMBDHVNacgturyswkmbdhvn", "TGCAAYRSWMKVHDBNtgcaayrswmkvhdbn")
        guide, pam = guide.translate(comp)[::-1], pam.translate(comp)[::-1]
    first, last = layout.cells[0] % length + 1, layout.cells[-1] % length + 1
    lines = [f"<b>{label}</b> &nbsp;<i>{nuc.name} guide site</i>" if label else f"<b>{nuc.name} guide site</b>",
             f"Guide (protospacer) 5′-{guide}-3′, {len(guide)} nt; PAM {pam} ({nuc.pam})",
             f"Bases {first:,}–{last:,}; the guide pairs with the {'bottom' if layout.strand > 0 else 'top'} strand and "
             f"reads along the {_word(layout.strand)} strand."]
    circ = True   # gaps are already reduced to 0..length
    if layout.top_gap == layout.bottom_gap:
        lines.append(f"{nuc.name} makes a blunt cut {_after(layout.top_gap, length, circ)} on both strands, "
                     f"{nuc.cut_from_pam} bp from the PAM, inside the protospacer.")
    else:
        lines.append(f"{nuc.name} cuts {_after(layout.top_gap, length, circ)} on the top strand and "
                     f"{_after(layout.bottom_gap, length, circ)} on the bottom strand, leaving a "
                     f"{abs(layout.bottom_gap - layout.top_gap)} nt 5′ overhang, away from the PAM.")
    return "<br>".join(lines)


def cas9_mark_tip(nuclease, start, strand, length, fetch, label=None):
    """The tooltip of a guide site found by the Spacer Finder (`start`: top-strand position of its first base)."""
    cells = [(start + i) % length for i in range(nuclease.length + len(nuclease.pam))]
    layout = cas9_layout(nuclease, cells, strand, length)
    return cas9_tip(layout, fetch, length, label) if layout else ""


# ---- restriction enzymes ------------------------------------------------------------------------------------------------


def locate_site(enzyme, seq, circular, top, bottom):
    """(start, end, strand, top, bottom) of the recognition site that makes the cut (`top`, `bottom`) (gaps, as find_cuts gives
    them), with the cuts counted on from the site's start (not reduced, so a cut across the origin keeps its distance to the
    site); None if no site of `enzyme` gives that cut."""
    n, L = len(seq), len(enzyme.site)
    for s, e, strand in site_spans(seq, circular, enzyme):
        pairs = enzyme.cuts if strand > 0 else tuple((L - b, L - t) for t, b in enzyme.cuts)
        for t, b in pairs:
            ct, cb = s + t, s + b
            if (ct, cb) == (top, bottom) or (circular and n and (ct - top) % n == 0 and (cb - bottom) % n == 0):
                return s, e, strand, ct, cb
    return None


def _distance(q, s, e):
    """Where a cut is relative to a site [s, e): inside it, or so many bases before / after."""
    if q >= e:
        return f"{q - e} nt after the site"
    if q <= s:
        return f"{s - q} nt before the site"
    return "inside the site"


def enzyme_tip(enzyme, seq, circular, top, bottom):
    """Tooltip (rich text) for the cut of `enzyme` (an enzymes.Enzyme) at (`top`, `bottom`) on `seq`: the recognition site with
    its cut marked (G^AATTC), the bases each strand is cut after (numbered from 1), the ends it leaves, and for an enzyme that
    cuts outside its site how far away."""
    n = len(seq)
    found = locate_site(enzyme, seq, circular, top, bottom)
    lines = [f"<b>{enzyme.name}</b> &nbsp;<i>restriction enzyme</i>"]
    if found:
        s, e, strand, ut, ub = found
        where = f"bases {s % n + 1:,}–{(e - 1) % n + 1:,}, " + ("top strand" if strand > 0 else f"on the bottom strand, where it reads {enzyme.site}; "
                                                               f"{revcomp_site(enzyme.site)} on the top strand")
        lines.append(f"Recognition site {notation(enzyme)} ({where})")
    else:
        ut, ub, s, e = top, bottom, None, None
        lines.append(f"Recognition site {notation(enzyme)}")
    lines.append(f"Cuts {_after(top, n, circular)} on the top strand and {_after(bottom, n, circular)} on the bottom strand.")
    if ut == ub:
        lines.append("Blunt ends (no overhang).")
    else:
        lo, hi = min(ut, ub), max(ut, ub)
        bases = "".join(seq[i % n] for i in range(lo, hi)) if n else ""
        kind = "5′" if ut < ub else "3′"
        lines.append(f"{kind} overhang of {hi - lo} nt: {bases} (read along the top strand).")
    if s is not None and (min(ut, ub) >= e or max(ut, ub) <= s):
        lines.append(f"Type IIS: it cuts outside its site, {_distance(ut, s, e)} on the top strand and "
                     f"{_distance(ub, s, e)} on the bottom strand.")
    return "<br>".join(lines)


def _rich(enzyme, seq, circular, top, bottom):
    """The hover tooltip of a restriction cut: the pill picture of the cut site and a line of words (cut_pill.py). Falls back to
    the plain description if the picture cannot be made (no GUI)."""
    try:
        from cut_pill import enzyme_rich_tip
        return enzyme_rich_tip(enzyme, seq, circular, top, bottom)
    except Exception:   # shortcut: a tooltip must never raise; the text version carries the same facts
        return enzyme_tip(enzyme, seq, circular, top, bottom)


def hover_site(editor, key):
    """What the sequence view draws while the pointer rests on a restriction cut: (start, end, top cut, bottom cut) of its
    recognition site, the cuts counted on from the site (not reduced, see locate_site); None if `key` (see hover.describe_cut)
    is not the cut of a known enzyme (a Cas9 guide has no recognition site)."""
    from enzymes import COMMON, find_cuts
    seq, circular = editor.doc.sequence(), editor.doc.circular
    if key[0] == "site":
        name, top = key[1], key[2]
        pair = None
    else:
        marks = editor.cut_marks.get(key[1], [])
        if key[2] >= len(marks):
            return None
        name, top = marks[key[2]].label, marks[key[2]].top
        pair = (top, marks[key[2]].bottom)
    enzyme = next((e for e in COMMON if e.name == name), None)
    if enzyme is None:
        return None
    if pair is None:
        pair = next(((t, b) for t, b in find_cuts(seq, circular, [enzyme]).get(name, []) if t == top), None)
        if pair is None:
            return None
    found = locate_site(enzyme, seq, circular, *pair)
    return (found[0], found[1], found[3], found[4]) if found else None


def unique_site_tip(editor, name, pos):
    """Tooltip for an enzyme name drawn on the sequence (View > Restriction Sites): the enzyme `name` cutting the top strand at
    gap `pos`."""
    from enzymes import COMMON, find_cuts
    enzyme = next((e for e in COMMON if e.name == name), None)
    if enzyme is None:
        return ""
    seq = editor.doc.sequence()
    pair = next(((t, b) for t, b in find_cuts(seq, editor.doc.circular, [enzyme]).get(name, []) if t == pos), None)
    return _rich(enzyme, seq, editor.doc.circular, *pair) if pair else ""


def enzyme_mark_tip(editor, name, top, bottom):
    """Tooltip for a digest's cut mark (enzyme `name`, top and bottom cut gaps)."""
    from enzymes import COMMON
    enzyme = next((e for e in COMMON if e.name == name), None)
    if enzyme is None:
        return ""
    return _rich(enzyme, editor.doc.sequence(), editor.doc.circular, top, bottom)
