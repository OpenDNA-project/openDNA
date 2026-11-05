"""Restriction enzymes: a built-in list of common ones (ENZYMES) and a site finder. The finder reads degenerate (IUPAC) sites and
enzymes that cut on both sides of their site, so a longer list only needs more entries.

Cut positions are offsets from the start of the recognition site, on the top strand
(`top`) and on the bottom strand (`bottom`), both counted on top-strand coordinates.
EcoRI G^AATTC cuts after base 1 on top and after base 5 on the bottom: (1, 5).
Type IIS enzymes cut outside their site: BsaI GGTCTC(1/5) -> (7, 11).
"""
import re
from collections import namedtuple

# name; site (IUPAC letters); cuts: ((top strand offset, bottom strand offset), ...) as offsets of the gap from the start of the
# site, one pair for the usual enzyme and two for one that cuts on both sides of its site; blunt: leaves flush ends.
Enzyme = namedtuple("Enzyme", "name site cuts blunt")

ENZYMES = [
    ("AatII", "GACGTC", 5, 1), ("AflII", "CTTAAG", 1, 5), ("AgeI", "ACCGGT", 1, 5),
    ("ApaI", "GGGCCC", 5, 1), ("AscI", "GGCGCGCC", 2, 6), ("AvrII", "CCTAGG", 1, 5),
    ("BamHI", "GGATCC", 1, 5), ("BbsI", "GAAGAC", 8, 12), ("BglII", "AGATCT", 1, 5),
    ("BsaI", "GGTCTC", 7, 11), ("BsmBI", "CGTCTC", 7, 11), ("BspHI", "TCATGA", 1, 5),
    ("BsrGI", "TGTACA", 1, 5), ("ClaI", "ATCGAT", 2, 4), ("EcoRI", "GAATTC", 1, 5),
    ("EcoRV", "GATATC", 3, 3), ("FseI", "GGCCGGCC", 6, 2), ("HindIII", "AAGCTT", 1, 5),
    ("HpaI", "GTTAAC", 3, 3), ("KpnI", "GGTACC", 5, 1), ("MfeI", "CAATTG", 1, 5),
    ("MluI", "ACGCGT", 1, 5), ("NcoI", "CCATGG", 1, 5), ("NdeI", "CATATG", 2, 4),
    ("NheI", "GCTAGC", 1, 5), ("NotI", "GCGGCCGC", 2, 6), ("NruI", "TCGCGA", 3, 3),
    ("NsiI", "ATGCAT", 5, 1), ("PacI", "TTAATTAA", 5, 3), ("PmeI", "GTTTAAAC", 4, 4),
    ("PstI", "CTGCAG", 5, 1), ("PvuI", "CGATCG", 4, 2), ("PvuII", "CAGCTG", 3, 3),
    ("SacI", "GAGCTC", 5, 1), ("SacII", "CCGCGG", 4, 2), ("SalI", "GTCGAC", 1, 5),
    ("SapI", "GCTCTTC", 8, 11), ("SbfI", "CCTGCAGG", 6, 2), ("ScaI", "AGTACT", 3, 3),
    ("SmaI", "CCCGGG", 3, 3), ("SpeI", "ACTAGT", 1, 5), ("SphI", "GCATGC", 5, 1),
    ("StuI", "AGGCCT", 3, 3), ("SwaI", "ATTTAAAT", 4, 4), ("XbaI", "TCTAGA", 1, 5),
    ("XhoI", "CTCGAG", 1, 5), ("XmaI", "CCCGGG", 1, 5),
    # short cutters (4 bp sites): they cut often, so the "6+" filter of the display leaves them out
    ("AluI", "AGCT", 2, 2), ("HaeIII", "GGCC", 2, 2), ("MspI", "CCGG", 1, 3), ("RsaI", "GTAC", 2, 2),
    ("Sau3AI", "GATC", 0, 4), ("TaqI", "TCGA", 1, 3),
]
_ACGT = str.maketrans("ACGT", "TGCA")
_IUPAC_COMPLEMENT = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")
_IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC", "B": "CGT",
          "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT"}
COMMON = [Enzyme(name, site, ((top, bottom),), top == bottom) for name, site, top, bottom in ENZYMES]


class EnzymeFilter(namedtuple("EnzymeFilter", "min_sites max_sites min_size max_size palindromic type_iis ends inventory",
                              defaults=(1, 1, 6, None, None, None, None, None))):
    """Which restriction sites the views draw (Display filter). Independent criteria, all of which must hold: how many times
    the enzyme cuts here (min_sites..max_sites), how long its recognition site is (min_size..max_size bases; a max of None is
    no limit), `palindromic` (the site reads the same on both strands), `type_iis` (cuts outside its site) and `inventory` (the
    enzyme is in the user's inventory, see inventory.py): True, False or None for either; `ends`: "5" / "3" overhang, "blunt"
    or None for any. The default is the unique cutters of 6 or more bases."""
    __slots__ = ()

    @classmethod
    def from_dict(cls, d):
        """The filter saved in a project (None if there is none)."""
        return cls(**{k: v for k, v in d.items() if k in cls._fields}) if isinstance(d, dict) else None

    def summary(self):
        """The filter in words, for the notes of the panels."""
        def span(lo, hi, unit):
            return f"{unit}: " + (f"{lo}" if lo == hi else f"{lo}–{hi}" if hi is not None else f"{lo}+")
        parts = [span(self.min_sites, self.max_sites, "sites"), span(self.min_size, self.max_size, "size")]
        for flag, name in ((self.palindromic, "palindromic"), (self.type_iis, "Type IIS")):
            if flag is not None:
                parts.append(name if flag else f"not {name}")
        if self.inventory is not None:
            parts.append("in my inventory" if self.inventory else "not in my inventory")
        if self.ends:
            parts.append({"5": "5′ overhang", "3": "3′ overhang", "blunt": "blunt"}[self.ends])
        return ", ".join(parts)


def is_palindromic(enzyme):
    return enzyme.site == revcomp_site(enzyme.site)


def is_type_iis(enzyme):
    """Cuts outside its recognition site (BsaI, BbsI...)."""
    n = len(enzyme.site)
    cuts = [c for pair in enzyme.cuts for c in pair]
    return min(cuts) >= n or max(cuts) <= 0


# The buttons under the criteria (Display filter): (label, filter, what it shows). Only a few: the criteria do the rest.
PRESETS = (
    ("Unique 6+", EnzymeFilter(), "Enzymes that cut once and recognise 6 or more bases (the default)"),
    ("Type IIS", EnzymeFilter(1, None, 1, None, None, True), "Every site of the enzymes that cut outside their site (BsaI, BbsI...)"),
    ("All sites", EnzymeFilter(1, None, 1), "Every site of every enzyme"),
)


def wanted(f, enzyme, count, owned=False):
    """Does the display filter `f` (an EnzymeFilter) draw the sites of `enzyme`, which cuts `count` times here? `owned`: the
    enzyme is in the user's inventory."""
    if f.inventory is not None and owned != f.inventory:
        return False
    size = len(enzyme.site)
    if not (f.min_sites <= count and (f.max_sites is None or count <= f.max_sites)):
        return False
    if not (f.min_size <= size and (f.max_size is None or size <= f.max_size)):
        return False
    if f.palindromic is not None and is_palindromic(enzyme) != f.palindromic:
        return False
    if f.type_iis is not None and is_type_iis(enzyme) != f.type_iis:
        return False
    if f.ends:
        top, bottom = enzyme.cuts[0]
        if f.ends != ("blunt" if top == bottom else "5" if top < bottom else "3"):
            return False
    return True


def revcomp_site(site):
    """The recognition site as it reads on the other strand (degenerate letters complemented too)."""
    return site.translate(_IUPAC_COMPLEMENT)[::-1]


def notation(e):
    """The site with its cuts marked the way suppliers write them: G^AATTC for EcoRI, GGTCTC(1/5) for a Type IIS enzyme (the two
    numbers are the bases between the site and the cut on each strand). An enzyme that also cuts upstream (BcgI) stays plain."""
    top, bottom = e.cuts[0]
    n = len(e.site)
    if top > n:                              # cuts downstream of its site
        return f"{e.site}({top - n}/{bottom - n})"
    if top < 0:                              # cuts upstream too (BcgI): keep the site plain
        return e.site
    return e.site[:top] + "^" + e.site[top:]


def describe(e):
    """One line on an enzyme: where it cuts and the ends it leaves ("GAATTC: G^AATTC, 5' overhang 4 nt"). An enzyme that
    cuts outside its site (Type IIS: BsaI, BbsI...) is written the way suppliers do, GGTCTC(1/5), so the overhang's sequence
    is whatever you put there."""
    top, bottom = e.cuts[0]
    ends = "blunt ends" if top == bottom else f"{'5′' if top < bottom else '3′'} overhang {abs(bottom - top)} nt"
    return f"{e.site}: {notation(e)}, {ends}" + (f", also cuts {len(e.cuts) - 1} more place(s)" if len(e.cuts) > 1 else "")


def _starts(text, site, limit):
    """Where `site` (IUPAC letters) starts in `text`, below `limit`. Plain A/C/G/T sites use str.find, which is much faster."""
    if all(c in "ACGT" for c in site):
        i = text.find(site)
        while i != -1 and i < limit:
            yield i
            i = text.find(site, i + 1)
        return
    pattern = re.compile("(?=(" + "".join(c if len(_IUPAC[c]) == 1 else f"[{_IUPAC[c]}]" for c in site) + "))")
    for m in pattern.finditer(text):
        if m.start() >= limit:
            break
        yield m.start()


def find_cuts(seq, circular, enzymes=None):
    """{enzyme name: sorted [(top cut, bottom cut)]} for every enzyme of `enzymes` (default: COMMON) that cuts: where it cuts each
    strand, as the position of the gap between two bases on the top strand's numbering (the cut falls before the base at that
    index). The two differ for sticky ends: EcoRI G^AATTC cuts the top strand after base 1 of its site and the bottom strand
    after base 5, so (site + 1, site + 5). A blunt cutter has the same number twice; an enzyme that cuts on both sides of its
    site has two pairs. Searches both strands (an enzyme that is not its own reverse complement is looked for in both
    orientations, the cuts mirrored); on circular sequences also finds sites spanning the origin."""
    seq = seq.upper()
    n = len(seq)
    out = {}
    for e in enzymes or COMMON:
        L = len(e.site)
        rc = revcomp_site(e.site)
        # Read the other way, the cut that was `t` from the start of the site is `L - t` from its end, on the other strand.
        orientations = [(e.site, e.cuts)] + ([] if rc == e.site else [(rc, tuple((L - b, L - t) for t, b in e.cuts))])
        reach = L + max(abs(x) for pair in e.cuts for x in pair) + 1   # bases past the end that a site across the origin needs
        text = seq + seq[:reach] if circular else seq
        cuts = set()
        for pattern, pairs in orientations:
            for i in _starts(text, pattern, n):
                for top, bottom in pairs:
                    t, b = i + top, i + bottom
                    if circular:
                        cuts.add((t % n, b % n))
                    elif 0 <= t <= n and 0 <= b <= n:   # a cut that falls off the end of a linear sequence is no cut
                        cuts.add((t, b))
        if cuts:
            out[e.name] = sorted(cuts)
    return out


def site_spans(seq, circular, enzyme):
    """[(start, end, strand)] of every recognition site of `enzyme` (an Enzyme), on either strand, in order of position; an enzyme
    that is its own reverse complement is listed once (strand 1). On a circular sequence a site can run across the origin (end
    passes the length)."""
    seq, L = seq.upper(), len(enzyme.site)
    rc = revcomp_site(enzyme.site)
    text = seq + seq[:L - 1] if circular else seq
    spans = [(i, i + L, 1) for i in _starts(text, enzyme.site, len(seq))]
    if rc != enzyme.site:
        spans += [(i, i + L, -1) for i in _starts(text, rc, len(seq))]
    return sorted(spans)


def find_sites(seq, circular):
    """{enzyme name: sorted top-strand cut positions} for every common enzyme that cuts (see find_cuts)."""
    return {name: sorted({t for t, _ in cuts}) for name, cuts in find_cuts(seq, circular).items()}


if __name__ == "__main__":
    # Self-check: EcoRI forward, BsaI on the reverse strand, a site across the origin.
    s = "TTGAATTCAA" + "GAGACCTTTTT" + "TCAA"
    sites = find_sites(s, circular=False)
    assert sites["EcoRI"] == [3], sites
    assert sites["BsaI"] == [5], sites  # GAGACC at 10: cuts 5 bases upstream
    assert find_sites("ATTCAAAAAAGA", circular=True)["EcoRI"] == [11]
    assert "EcoRI" not in find_sites("ATTCAAAAAAGA", circular=False)
    print("ok")
