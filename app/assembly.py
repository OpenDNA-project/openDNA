"""Assembly engine: what a set of fragments (PCR products, or plain stretches) joins into with
Gibson / HiFi, Golden Gate or Sidewinder: the junctions found at their real ends, the checks, the
primers that would add what is missing, and the exact product. The widgets are in assembly_panel.py."""
import random
from bisect import bisect_right
from collections import namedtuple
from html import escape
from types import SimpleNamespace

from PySide6.QtGui import QColor

from opendna_core import Document, find_primers

import overhangs as OH
from editor import MAX_SITES, reverse_complement as rc, site_location
from enzymes import find_sites
from primer_panel import MARKS, band, limit, melting_temp, structure, worst
import theme

# name -> (enzyme in enzymes.py, recognition site, spacer between site and overhang, overhang length)
GG_ENZYMES = {"BsaI": ("BsaI", "GGTCTC", 1, 4), "BsmBI / Esp3I": ("BsmBI", "CGTCTC", 1, 4),
              "BbsI": ("BbsI", "GAAGAC", 2, 4), "SapI": ("SapI", "GCTCTTC", 1, 3)}
PAD = "ttt"                # extra bases 5' of a Golden Gate site so the enzyme can bind near an end
MIN_EXISTING = 15          # bp: adjacent fragments overlapping at least this much need no tail
TOEHOLD, BARCODE, MAX_OLIGO = 10, 18, 165   # Sidewinder: toehold nt, barcode nt, longest oligo
MAX_SHIFT = 20             # nt a Sidewinder toehold may slide from the fragment boundary to dodge a repeat
# Orange is kept for what a primer adds (TAIL_COLOR, like the tails written in orange in the sequence view), so no fragment is orange;
# a junction (the homology, overhang or toehold) is a dark slate no fragment has.
FRAG_COLORS = ["#4fc3f7", "#9ccc65", "#fdd835", "#ba68c8", "#4db6ac", "#f06292", "#a1887f", "#7986cb"]
TAIL_COLOR, JOINT_COLOR = "#ef6c00", "#546e7a"


def snapshot(doc, name=None):
    """An independent copy of a Document (bases, features, topology): fragments never share
    state with a file that is open in an editor."""
    out = Document.new_sequence(name or doc.name, doc.sequence(0, None), doc.circular)
    feats = doc.copy_features(0, doc.length)
    if feats:
        out.add_features([(f.kind, f.location, f.qualifiers) for f in feats])
    return out


class Frag:
    """One fragment of the assembly. `src` is a copy of the file it comes from; the fragment is
    the stretch [start, end) of it (end may pass the origin of a circular source), made by PCR with
    `fwd` / `rev` = (name, full 5'->3' sequence, annealing length) in the fragment's own orientation,
    or just cut out when there are none. The 5' part of a primer that does not anneal is its tail.
    `vector`: a Golden Gate destination vector (kept whole, cut at its two sites by the design).
    `doc` is the finished fragment (tails included), `core` the stretch without the tails."""
    def __init__(self, name, src, start=0, end=None, fwd=None, rev=None, vector=False, source=""):
        self.name, self.src, self.start = name, src, start
        self.end = src.length if end is None else end
        self.fwd, self.rev, self.vector, self.source = fwd, rev, vector, source
        self.flipped = False
        self.rebuild()

    def rebuild(self):
        """Make `doc` and `core` from the source, the primers and the orientation."""
        if self.vector:
            self.doc = snapshot(self.src, self.name)
            self.doc.circular = True    # a vector is a plasmid: its backbone may span the origin
            self.core, self.tails = self.doc, (0, 0)
            if self.flipped:
                self.doc.reverse_complement()
            return
        if self.end <= self.start:   # an empty stretch: a fragment with no bases, which the checks then report (the core cannot cut it)
            self.doc = self.core = Document.new_sequence(self.name, "", False)
            self.tails = (0, 0)
            return
        tail = lambda p: p[1][:len(p[1]) - p[2]].upper() if p else ""
        tf, tr = tail(self.fwd), tail(self.rev)
        # `fwd`/`rev` are in the final orientation: the final fragment is tf + core + rc(tr). Cut from
        # the source it is that, or (flipped) its reverse complement, which is then flipped back.
        pre, suf = (tr, rc(tf)) if self.flipped else (tf, rc(tr))
        self.doc = self.src.pcr_product(self.name, self.start, self.end, pre, suf)
        self.core = self.src.pcr_product(self.name, self.start, self.end, "", "")
        if self.flipped:
            self.doc.reverse_complement()
            self.core.reverse_complement()
        self.tails = (len(tf), len(tr))

    def flip(self):
        """Reverse-complement the fragment: the same molecule read from the other strand."""
        self.fwd, self.rev = self.rev, self.fwd
        self.flipped = not self.flipped
        self.rebuild()

    @property
    def seq(self):
        return self.doc.sequence(0, None).upper()

    @property
    def core_seq(self):
        return self.core.sequence(0, None).upper()


# One junction between fragment `left` and the next one, `right`. `kind`: what joins them
# ("existing overlap", "overhang", "toehold") or "none" when the two ends do not fit. `seq`: the
# bases that join them, starting at `pos` in the product (counted on from 0, so a junction that
# closes a circle sits at len(product) and wraps round). `a` / `b`: how the left and the right
# fragment's two strands end there, as (top strand edge, bottom strand edge, unused). The edges
# are product positions: the left fragment's strands end there, the right one's start there.
Joint = namedtuple("Joint", "left right kind seq pos a b")


class Design:
    """The result of one analysis: oligos, checks, junction texts and how the product is made."""
    def __init__(self, frags, circular):
        self.names = [f.name for f in frags]
        self.flipped = [f.flipped for f in frags]   # reverse-complemented from what was added
        self.note = ""       # one sentence on how this method joins the fragments
        self.seqs = [f.seq for f in frags]   # as used in the product (a vector: its backbone; Golden Gate: digested)
        self.offsets = [0] * len(frags)      # where self.seqs[i] starts inside the fragment's own bases
        self.ranges = [None] * len(frags)    # per fragment: (start, end) of self.seqs[i] on its document
        self.pieces = [(i, 0, 0) for i in range(len(frags))]  # (fragment, bases cut at start, at end)
        # Per fragment: how many bases before / after its piece are still its own bases in the product, because the
        # neighbour's piece carries the very same bases there (a Gibson overlap, a Golden Gate overhang). The product keeps
        # the fragment's features over them, or an ORF that starts in the overlap would start late and lose its frame.
        self.shared = [(0, 0)] * len(frags)
        self.tails = [[] for _ in frags]     # per fragment: (start, end) in the product of bases a primer tail added
        self.tail_primers = [[] for _ in frags]   # ... and, in the same order, whose they are: (primer name, "start" | "end" of the fragment)
        self.circular = circular
        self.rows = []       # (name, sequence, length of annealing part (0 = all), kind)
        self.checks = []     # (check, level, text)
        self.junctions = []
        self.joints = []     # one Joint per junction, in order
        self.problems = []   # per junction: why the two ends do not fit, or None
        self.jlevels = []    # per junction: ok / warn / bad
        self.cuts = []       # per junction: lines of text that show how it is made (Golden Gate: the two ends with site and cut)
        self.jchecks = {}    # junction -> [(short name, level, text)]: what its card shows
        self.fchecks = {}    # fragment -> [(short name, level, text)]
        self.local = set()   # indexes in self.checks of the ones shown on a junction card or a fragment box

    def product(self):
        return "".join(self.seqs[i][a:len(self.seqs[i]) - b] for i, a, b in self.pieces)

    def check(self, name, level, text, junction=None, fragment=None, short=None):
        """Add a check. `junction` / `fragment`: also show it on that junction's card / fragment's box."""
        self.checks.append((name, level, text))
        if junction is not None or fragment is not None:
            self.local.add(len(self.checks) - 1)
        if junction is not None:
            self.jchecks.setdefault(junction, []).append((short or name, level, text))
        if fragment is not None:
            self.fchecks.setdefault(fragment, []).append((short or name, level, text))

    @property
    def global_checks(self):
        """The checks that are not about one junction or one fragment (the product, the circle, leftover sites...)."""
        return [c for i, c in enumerate(self.checks) if i not in self.local]

    @property
    def ok(self):
        """Every junction fits and nothing is red: the assembly can be performed."""
        return bool(self.joints or len(self.names) == 1) and all(j.kind != "none" for j in self.joints) \
            and not any(level == "bad" for _, level, _ in self.checks)


# ---- helpers ----------------------------------------------------------------------------------


def anneal_len(seq, target):
    """Length (18-35) of the start of `seq` whose Tm reaches `target`; as much as there is if short."""
    for n in range(18, min(35, len(seq)) + 1):
        t = melting_temp(seq[:n])
        if t is not None and t >= target:
            return n
    return min(35, len(seq))


def add_primers(d, name, seq, left, right, target, kind, left_ext="", right_ext=""):
    """Forward/reverse primer for `seq`. `left`/`right` = neighbour bases added 5' of the
    annealing part, `*_ext` = bases in front of those (site, spacer). Extensions are lowercase."""
    nf, nr = anneal_len(seq, target), anneal_len(rc(seq), target)
    d.rows.append((f"{name}-F", (left_ext + left).lower() + seq[:nf], nf, kind))
    d.rows.append((f"{name}-R", (right_ext + rc(right)).lower() + rc(seq)[:nr], nr, kind))


def oligo_checks(row):
    """Tm and GC of the annealing part, hairpin and self-dimer of the whole oligo."""
    _, seq, n, _ = row
    core = seq[len(seq) - n:] if n else seq
    out, t = [], melting_temp(core)
    if t is not None:
        out.append(("Tm (annealing part)", band(t, limit("tm_low"), limit("tm_high"), limit("tm_margin")),
                    f"{t:.1f} °C ({len(core)} nt)"))
    gc = 100 * sum(c in "GC" for c in core.upper()) / max(len(core), 1)
    out.append(("GC", band(gc, limit("gc_low"), limit("gc_high"), limit("gc_margin")), f"{gc:.0f} %"))
    if len(seq) >= 8:
        s = structure(seq)
        out.append(("Hairpin", hairpin_level(s), f"ΔG {min(s['hairpin'], 0):.1f} kcal/mol"))
        out.append(("Self-dimer", "bad" if s["self3"] < limit("dimer3_bad") or s["self"] < limit("dimer_bad") else
                    "warn" if any("dimer" in w for w in s["warnings"]) else "ok",
                    f"ΔG {s['self']:.1f} (3′ end {s['self3']:.1f}) kcal/mol"))
    return out


def hairpin_level(s):
    return "bad" if s["hairpin"] < limit("hairpin_bad") else "warn" if any(w.startswith("hairpin") for w in s["warnings"]) else "ok"


def occurrences(text, pattern, circular):
    """How often `pattern` or its reverse complement occurs in `text` (looping if circular)."""
    if circular:
        text += text[:len(pattern) - 1]
    count = lambda p: sum(text.startswith(p, i) for i in range(len(text) - len(p) + 1))
    return count(pattern) + (0 if rc(pattern) == pattern else count(rc(pattern)))


def shared_kmer(a, b, k=10):
    """True if a and b (either strand) share a k-mer."""
    kmers = {a[i:i + k] for i in range(len(a) - k + 1)}
    return any(b[i:i + k] in kmers or rc(b[i:i + k]) in kmers for i in range(len(b) - k + 1))


def joint_tm(seq):
    """Tm of the bases joining two fragments, or None when there are too few for it to mean
    anything: a 4 nt Golden Gate overhang is held by the ligase, not by its melting temperature."""
    return melting_temp(seq) if len(seq) >= 8 else None


def joint_label(joint):
    """Short text for the fragment map: how many bases join the two fragments, and how warm."""
    kind, seq = joint.kind, joint.seq
    if kind == "none":
        return "✗ no match"
    t = joint_tm(seq)
    return f"{len(seq)} bp " + (f"· {t:.0f} °C" if t is not None else kind)


def joints_html(d):
    """The junction table: what joins each pair of fragments, how long it is and how warm."""
    head = ("<tr><th align='left'>From</th><th align='left'>To</th><th align='left'>Joined by</th>"
            "<th align='left'>Sequence 5′→3′</th><th align='right'>bp</th>"
            "<th align='right'>Tm</th><th align='right'>GC</th></tr>")
    rows = []
    for left, right, kind, seq, *_ in d.joints:
        if kind == "none":
            continue
        t = joint_tm(seq)
        gc = 100 * sum(c in "GC" for c in seq.upper()) / max(len(seq), 1)
        rows.append(f"<tr><td>{escape(left)}&nbsp;&nbsp;</td><td>{escape(right)}&nbsp;&nbsp;</td>"
                    f"<td>{escape(kind)}&nbsp;&nbsp;</td><td><tt>{escape(seq)}</tt>&nbsp;&nbsp;</td>"
                    f"<td align='right'>{len(seq)}&nbsp;&nbsp;</td>"
                    f"<td align='right'>{f'{t:.1f} °C' if t is not None else '—'}&nbsp;&nbsp;</td>"
                    f"<td align='right'>{gc:.0f} %</td></tr>")
    return f"<table cellspacing='1' cellpadding='1'>{head}{''.join(rows)}</table>"


JOINED_BY = {"existing overlap": "homology", "added overlap": "homology", "overhang": "overhang", "toehold": "toehold"}
TAIL_BG, TAIL_FG = "#fff3e0", "#bf360c"   # bases a primer tail adds: lowercase, orange on cream


COMP = str.maketrans("ACGTacgt", "TGCAtgca")


def sketch_html(d, j):
    """A junction the way it is drawn in the manual: the left fragment's top strand running into the
    joining bases, and under them their complement running into the right fragment."""
    jt, n = d.joints[j], len(d.names)
    left, right = f"[{d.names[j % n][:20]}]", f"[{d.names[(j + 1) % n][:20]}]"
    pad = "&nbsp;" * len(left)
    top = f"<b>{escape(left)}</b><span style='font-weight:600;text-decoration:underline'>{escape(jt.seq)}</span>"
    pipes = pad + "|" * len(jt.seq)
    bottom = pad + escape(jt.seq.translate(COMP)) + f"<b>{escape(right)}</b>"
    return f"<tt>{top}<br>{pipes}<br>{bottom}</tt>"


def junction_html(d, context=16, only=None):
    """Every junction as the bases around it on both strands, the way the method uses them:
    the product (each base coloured by the fragment it comes from, the joining bases bold and
    underlined) and, for the fragment on each side, its top and bottom strand where they stop.
    Lowercase on orange = bases a primer tail adds. Selectable, so a sequence can be copied."""
    if not d.joints:
        return "<p>No junctions.</p>"
    product, n = d.product(), len(d.names)
    N = len(product)
    ends, total = [], 0
    for i, a, b in d.pieces:
        total += len(d.seqs[i]) - a - b
        ends.append(total)
    base = lambda q: product[q % N] if d.circular else (product[q] if 0 <= q < N else "")
    color = lambda i: FRAG_COLORS[i % len(FRAG_COLORS)]
    light = lambda c: QColor(c).lighter(150).name()
    is_tail = lambda idx, q: any(a <= q % N < b for a, b in d.tails[idx])

    def cell(ch, bg=None, fg="#1b1b1b", bold=False):
        if not ch:
            return "&nbsp;"
        style = (f"background-color:{bg};color:{fg};" if bg else "") + ("font-weight:600;text-decoration:underline;" if bold else "")
        return f"<span style='{style}'>{escape(ch)}</span>"

    out = [f"<p>{escape(d.note)}</p><p>Lowercase on orange = bases added by a primer tail. "
           "Bold and underlined = what joins the two fragments.</p>"]
    for j, jt in enumerate(d.joints):
        if only is not None and j != only:
            continue
        ia, ib = j % n, (j + 1) % n
        orient = lambda i: f"{escape(d.names[i])} {'◀ reverse-complemented' if d.flipped[i] else '▶ as given'}"
        if jt.kind == "none":
            reason = d.problems[j] if j < len(d.problems) and d.problems[j] else "the two ends do not fit"
            out.append(f"<p><b>{j + 1}. {escape(jt.left)} → {escape(jt.right)}</b> — "
                       f"<span style='color:{theme.hex('status_bad')};font-weight:600'>✗ {escape(reason)}</span><br>"
                       f"{orient(ia)} &nbsp;·&nbsp; {orient(ib)}</p>"
                       f"<p><tt>{escape(d.names[ia][:20])} ends&nbsp;&nbsp;…{escape(d.seqs[ia][-30:])}<br>"
                       f"{escape(d.names[ib][:20])} starts {escape(d.seqs[ib][:30])}…</tt></p>")
            continue
        L = len(jt.seq)
        lo, hi = jt.pos - context, jt.pos + L + context
        word = JOINED_BY[jt.kind]
        inside = lambda q: jt.pos <= q < jt.pos + L
        t = joint_tm(jt.seq)
        gc = 100 * sum(c in "GC" for c in jt.seq.upper()) / max(L, 1)
        level = d.jlevels[j] if j < len(d.jlevels) else "ok"
        mark = {"ok": "✓", "warn": "⚠", "bad": "✗"}.get(level, "")
        colour = theme.hex({"ok": "status_ok", "warn": "status_warn", "bad": "status_bad"}.get(level, "status_info"))
        out.append(f"<p><b>{j + 1}. {escape(jt.left)} → {escape(jt.right)}</b> — <span style='color:{colour};font-weight:600'>{mark}</span> "
                   f"{jt.kind}, {L} bp" + (f", Tm {t:.1f} °C" if t is not None else "") + f", GC {gc:.0f} %<br>"
                   f"{orient(ia)} &nbsp;·&nbsp; {orient(ib)}</p>")
        out.append(sketch_html(d, j))
        if j < len(d.cuts) and d.cuts[j]:
            out.append("<p><tt>" + "<br>".join(escape(line).replace(" ", "&nbsp;") for line in d.cuts[j]) + "</tt></p>")
        rows = []

        def owner(q):
            return min(bisect_right(ends, q % N), n - 1)

        qs = range(lo, hi)
        for text, bottom in (("product 5′→3′", False), ("3′←5′", True)):
            cells = "".join(cell(rc(base(q)) if bottom and base(q) else base(q),
                                 (light if bottom else str)(color(owner(q))) if base(q) else None, bold=inside(q))
                            for q in qs)
            rows.append((text, None, cells))
        marks = [" "] * (hi - lo)
        for q in range(jt.pos, jt.pos + L):
            marks[q - lo] = "─"
        marks[jt.pos - lo], marks[jt.pos + L - 1 - lo] = ("│", "│") if L == 1 else ("└", "┘")
        label = f" {L} bp {word} "
        if len(label) <= L - 2:
            start = jt.pos - lo + 1 + (L - 2 - len(label)) // 2
            marks[start:start + len(label)] = label
        else:   # too short to hold it: write it just after the bracket
            start = jt.pos + L - lo + 1
            marks[start:] = label.strip()[:len(marks) - start]
        rows.append(("", None, escape("".join(marks)).replace(" ", "&nbsp;")))
        for idx, side, left in ((ia, jt.a, True), (ib, jt.b, False)):
            for text, edge in ((f"{d.names[idx][:14]} top", side[0]), (f"{d.names[idx][:14]} bottom", side[1])):
                bottom = text.endswith("bottom")
                cells = ""
                for q in qs:
                    ch = base(q) if (q < edge if left else q >= edge) else ""
                    if bottom and ch:
                        ch = rc(ch)
                    if ch and is_tail(idx, q):
                        cells += cell(ch.lower(), TAIL_BG, TAIL_FG)
                    else:
                        cells += cell(ch, (light if bottom else str)(color(idx)) if ch else None)
                rows.append((text, color(idx), cells))

        def label_cell(text, chip):
            if not chip:
                return escape(text)
            return f"<span style='background-color:{chip};color:#1b1b1b;'>&nbsp;{escape(text)}&nbsp;</span>"

        out.append("<table cellspacing='0' cellpadding='0'>" + "".join(
            f"<tr><td>{label_cell(text, chip)}&nbsp;&nbsp;</td><td><tt>{cells}</tt></td></tr>"
            for text, chip, cells in rows) + "</table><br>")
    return "".join(out)


def check_product(d):
    """The product itself: size and topology, and — when it should be circular — that the bases
    joining the last fragment back to the first really are there, once, in the finished loop."""
    product = d.product()
    missing = sum(j.kind == "none" for j in d.joints)
    d.check("Product", "info", f"{len(product):,} bp, {'circular' if d.circular else 'linear'}, "
                               f"{len(d.names)} fragment(s), {len(d.joints)} junction(s)"
                               + (f", {missing} not joined" if missing else ""))
    if not d.circular:
        d.check("Ends", "info", "a linear product: its two ends are not joined to each other")
        return
    if len(d.joints) != len(d.names) or d.joints[-1].kind == "none":
        d.check("Circle closes", "bad", f"{len(d.names) - missing} of {len(d.names)} junctions fit"
                                        + (": the last fragment is not joined back to the first"
                                           if len(d.joints) == len(d.names) and d.joints[-1].kind == "none" else ""))
        return
    left, right, _, seq, *_ = d.joints[-1]
    # A long join must also be the only place it could happen. A 4 nt Golden Gate overhang
    # occurs by chance all over any product, and its specificity comes from where the enzyme
    # cuts, not from the bases being unique, so it is checked against the other overhangs
    # ("Overhangs unique") instead of against the product.
    hits = occurrences(product, seq, True) if len(seq) >= 8 else 1
    d.check("Circle closes", "ok" if hits == 1 else "bad",
            f"every fragment is joined on both sides; {left} → {right} closes the loop with "
            f"{len(seq)} bp ({seq})"
            + ("" if hits == 1 else f", but those bases occur {hits}× in the product: the ends could join elsewhere"))


def check_lengths(d, minimum=200):
    small = [(i, n) for i, (n, s, v) in enumerate(zip(d.names, d.seqs, d.vectors)) if len(s) < minimum and not v]
    d.check("Fragment length", "warn" if small else "ok",
            f"shorter than {minimum} bp (hard to handle): " + ", ".join(n for _, n in small) if small else f"all at least {minimum} bp")
    for i, _ in small:
        d.fchecks.setdefault(i, []).append(("short", "warn", f"shorter than {minimum} bp: hard to handle"))


# ---- Gibson / HiFi ------------------------------------------------------------------------------


def design_gibson(frags, circular, min_ov=20, max_ov=40, min_tm=52.0, target=60.0):
    d = Design(frags, circular)
    d.vectors = [False] * len(frags)
    n, seqs = len(frags), d.seqs
    left, right, skip, trim, overlaps, pending = [""] * n, [""] * n, [0] * n, [0] * n, [], []
    for j in range(n if circular else n - 1):
        a, b = j, (j + 1) % n
        sa, sb = seqs[a], seqs[b]
        # Longest existing overlap: end of A = start of B (no tail needed then).
        ex = 0 if a == b else next((k for k in range(min(60, len(sa), len(sb)), MIN_EXISTING - 1, -1) if sa[-k:] == sb[:k]), 0)
        ha = hb = 0
        if ex:
            ov = sa[-ex:]
            if b == 0:
                trim[a] = ex     # closing the circle: drop the duplicate from the last fragment
            else:
                skip[b] = ex
            d.junctions.append(f"{d.names[a]} → {d.names[b]}: existing {ex} bp overlap, no tail")
        else:
            for length in range(min_ov, max_ov + 1):   # shortest overlap that reaches the Tm
                ha, hb = length // 2, length - length // 2
                t = melting_temp(sa[-ha:] + sb[:hb])
                if t is not None and t >= min_tm:
                    break
            right[a], left[b] = sb[:hb], sa[-ha:]
            ov = sa[-ha:] + sb[:hb]
            d.junctions.append(f"{d.names[a]} → {d.names[b]}: {len(ov)} bp overlap added "
                               f"({hb} nt on {d.names[a]}-R, {ha} nt on {d.names[b]}-F)")
        overlaps.append((f"{d.names[a]} → {d.names[b]}", ov))
        pending.append((ex, ha, hb, b == 0, ov))
    d.pieces = [(i, skip[i], trim[i]) for i in range(n)]
    d.note = ("Each fragment is amplified with the start (or end) of its neighbour as a 5′ tail (lowercase). "
              "The exonuclease chews back the 5′ ends and the homologous 3′ overhangs anneal.")
    ends, total = [], 0
    for i in range(n):
        total += len(seqs[i]) - skip[i] - trim[i]
        ends.append(total)
    for j, (ex, ha, hb, closing, ov) in enumerate(pending):
        edge = ends[j]       # where fragment j ends in the product
        pos = edge if ex and closing else edge - (ex or ha)
        top = pos + len(ov)  # both fragments carry the whole overlap: the left one up to here, the right one from `pos`
        d.joints.append(Joint(d.names[j], d.names[(j + 1) % n], "existing overlap" if ex else "added overlap", ov, pos,
                              (top, top, (edge, top) if hb else None), (pos, pos, (pos, edge) if ha else None)))
    for i in range(n):
        add_primers(d, d.names[i], seqs[i], left[i], right[i], target, "Gibson")
    product = d.product()
    for label, ov in overlaps:
        t = melting_temp(ov)
        gc = 100 * sum(c in "GC" for c in ov) / len(ov)
        d.check(f"Overlap Tm ({label})", band(t, min_tm, 80, 4) if t else "warn", f"{t:.1f} °C, {len(ov)} bp" if t else "?")
        d.check(f"Overlap GC ({label})", band(gc, 35, 65, 10), f"{gc:.0f} %")
        hits = occurrences(product, ov, circular)
        d.check(f"Overlap unique ({label})", "ok" if hits == 1 else "bad",
                "occurs once in the product" if hits == 1 else f"occurs {hits} times in the product (wrong joins)")
        d.check(f"Overlap hairpin ({label})", hairpin_level(structure(ov)), f"ΔG {min(structure(ov)['hairpin'], 0):.1f}")
    clash = [f"{overlaps[i][0]} / {overlaps[j][0]}" for i in range(len(overlaps)) for j in range(i)
             if shared_kmer(overlaps[i][1], overlaps[j][1])]
    d.check("Overlaps distinct", "warn" if clash else "ok", "share a 10-mer: " + "; ".join(clash) if clash else "no two share a 10-mer")
    check_product(d)
    check_lengths(d)
    return d


# ---- Golden Gate --------------------------------------------------------------------------------


def backbone(seq, enzyme, doc_length):
    """Destination vector between its two sites: (bases incl. both overhangs, start on the
    vector, error). The backbone is the arc that holds no site (the vector without the dropout)."""
    name, _, _, k = GG_ENZYMES[enzyme]
    cuts = find_sites(seq, True).get(name, [])
    if len(cuts) != 2:
        return None, 0, f"destination vector has {len(cuts)} {name} cut(s), needs exactly 2"
    c1, c2 = cuts
    arcs = [(c1, c2 - c1 + k), (c2, len(seq) - c2 + c1 + k)]
    arcs = [((seq + seq)[s:s + n], s) for s, n in arcs]
    arcs.sort(key=lambda a: (len(find_sites(a[0], False).get(name, [])), -len(a[0])))
    return arcs[0][0], arcs[0][1], None


def pick_overhangs(cands, fixed, k):
    """One (position, overhang) per junction, greedily: no palindromes, none equal to or one
    base off another overhang (or its reverse complement), not all G/C or all A/T, and the
    overhangs of the published high-fidelity set (see overhangs.py) first."""
    chosen = []
    for cand, fix in zip(cands, fixed):
        best = None
        for p, o in [fix] if fix else cand:
            pen = OH.preference(o, [q for _, q in chosen]) + 0.01 * p
            if best is None or pen < best[0]:
                best = (pen, p, o)
        chosen.append(best[1:])
    return chosen


def design_golden_gate(frags, circular, enzyme, target=60.0):
    ename, site, spacer, k = GG_ENZYMES[enzyme]
    d = Design(frags, circular)
    d.vectors = [f.vector for f in frags]
    for i, f in enumerate(frags):
        if f.vector:
            d.seqs[i], start, err = backbone(f.seq, enzyme, f.doc.length)
            if err:
                raise ValueError(f"{f.name}: {err}")
            d.ranges[i] = (start, start + len(d.seqs[i]))
    n, seqs = len(frags), d.seqs
    product = d.product()
    N = len(product)
    starts = [sum(map(len, seqs[:i])) for i in range(n)]
    ends = [s + len(q) for s, q in zip(starts, seqs)]
    njunc = n if circular else n - 1
    at = lambda i: product[i % N]                    # the product as a loop
    cands, fixed = [], []
    for j in range(njunc):
        b = ends[j]                                  # junction j sits after fragment j
        cands.append([(p, "".join(at(p + x) for x in range(k))) for p in range(b - k, b + 1)])
        # Next to a vector the overhang is the vector's own end bases.
        fixed.append(cands[j][-1] if d.vectors[(j + 1) % n] else cands[j][0] if d.vectors[j] else None)
    chosen = pick_overhangs(cands, fixed, k)
    overhangs = [o for _, o in chosen]
    d.note = ("Each fragment is digested next to its sites. Its top strand ends where the 4 nt overhang starts and "
              "its bottom strand carries the overhang (a 5′ overhang on the bottom strand of the left fragment, on the "
              "top strand of the right one); matching overhangs are ligated.")
    ext = (PAD + site + "a" * spacer).lower()
    for i in range(n):
        if d.vectors[i]:
            continue
        has_l, has_r = i > 0 or circular, i < njunc
        lt = starts[i] - chosen[i - 1][0] if has_l else 0   # neighbour bases before this fragment
        rt = chosen[i][0] + k - ends[i] if has_r else 0     # ... and after it
        left = "".join(at(starts[i] - lt + x) for x in range(lt))
        right = "".join(at(ends[i] + x) for x in range(rt))
        add_primers(d, d.names[i], seqs[i], left, right, target, f"Golden Gate {ename}",
                    ext if has_l else "", ext if has_r else "")
    for j, o in enumerate(overhangs):
        d.junctions.append(f"{d.names[j]} → {d.names[(j + 1) % n]}: overhang {o}")
        p = chosen[j][0]   # the left fragment's top strand ends here, its bottom strand k bases on
        d.joints.append(Joint(d.names[j], d.names[(j + 1) % n], "overhang", o, p, (p, p + k, None), (p, p + k, None)))
    for f, s, v in zip(frags, seqs, d.vectors):
        if not v:
            cuts = find_sites(s, False).get(ename, [])
            d.check(f"{ename} sites in {f.name}", "bad" if cuts else "ok",
                    "internal site(s) at " + ", ".join(f"{c + 1:,}" for c in cuts) + " (domesticate: mutate them)" if cuts else "none")
    left_over = find_sites(product, circular).get(ename, [])
    d.check(f"{ename} sites left in product", "warn" if left_over else "ok",
            "at " + ", ".join(f"{c + 1:,}" for c in left_over) + " (the product can be re-cut)" if left_over else "none")
    dup = sorted({o for o in overhangs if overhangs.count(o) > 1})
    d.check("Overhangs unique", "bad" if dup else "ok", "repeated: " + ", ".join(dup) if dup else "all different")
    pal = [o for o in overhangs if o == rc(o)]
    d.check("Overhangs not palindromic", "bad" if pal else "ok", "palindromic: " + ", ".join(pal) if pal else "none")
    rcs = [f"{o}/{q}" for i, o in enumerate(overhangs) for q in overhangs[:i] if o == rc(q)]
    d.check("No overhang is another's reverse complement", "bad" if rcs else "ok", ", ".join(rcs) or "none")
    near = [f"{o}~{q}" for i, o in enumerate(overhangs) for q in overhangs[:i]
            if any(o != r and sum(x == y for x, y in zip(o, r)) >= k - 1 for r in (q, rc(q)))]
    d.check("Overhangs differ in 2+ positions", "warn" if near else "ok",
            f"share {k - 1} of {k} bases: " + ", ".join(near) if near else "all differ by 2 or more")
    weak = [o for o in overhangs if set(o) <= set("GC") or set(o) <= set("AT")]
    d.check("Overhang GC mix", "warn" if weak else "ok", "all G/C or all A/T: " + ", ".join(weak) if weak else "mixed")
    check_product(d)
    check_lengths(d)
    return d


# ---- Sidewinder ---------------------------------------------------------------------------------


def toehold_cuts(product, bounds, lengths):
    """Where each coding strand ends: half a toehold past the fragment boundary, slid by up to
    MAX_SHIFT nt when the toehold there occurs more than once in the product or sits within 2 bases
    of one already chosen. The toehold is taken from the product wherever it sits, so sliding stays
    scarless: it only moves which oligo carries those bases. Fragments too short to give the slide
    room keep the boundary itself, so a tiny or empty one still fails the toehold checks."""
    # shortcut: greedy left to right, so a clash only a later junction could resolve is left for the checks to report;
    # upgrade to a search over all junctions together if real designs hit it.
    half, cuts, chosen = TOEHOLD // 2, [], []
    for j, b in enumerate(bounds):
        room = max(0, min(MAX_SHIFT, lengths[j] // 2, lengths[j + 1] // 2))
        low = (cuts[-1] if cuts else 0) + TOEHOLD      # every coding strand must hold a whole toehold
        wanted = [b + half + s for s in sorted(range(-room, room + 1), key=abs) if low <= b + half + s <= len(product)]
        fits = lambda t: occurrences(product, t, False) == 1 and all(
            min(sum(x != y for x, y in zip(t, r)) for r in (u, rc(u))) >= 3 for u in chosen)
        cut = next((c for c in wanted if fits(product[c - TOEHOLD:c])), wanted[0] if wanted else b + half)
        cuts.append(cut)
        chosen.append(product[max(0, cut - TOEHOLD):cut])
    return cuts


def design_sidewinder(frags, target=60.0):
    """Sidewinder (Robinson et al., Nature 2026; bioRxiv 2026 oligo-pool follow-up): each fragment is
    two annealed oligos, a coding strand (the final sequence) and a barcode strand (its reverse
    complement plus an 18 nt barcode arm that builds the third helix of a three-way junction with
    the neighbour). Neighbours meet at 10 nt toeholds taken from the product itself, so the result
    is scarless; Taq ligase seals the junctions and PCR removes the barcode arms. Simplified: coding
    strand i runs from the end of toehold i-1 to the end of toehold i; barcode strand i is its
    reverse complement plus barcode i (none on the last fragment). Check the strand layout
    against the paper before ordering."""
    d = Design(frags, False)
    d.vectors = [False] * len(frags)
    n = len(frags)
    product = d.product()
    bounds = [sum(map(len, d.seqs[:i + 1])) for i in range(n - 1)]
    cuts = toehold_cuts(product, bounds, [len(s) for s in d.seqs])
    toeholds = [product[max(0, c - TOEHOLD):c] for c in cuts]   # shorter than TOEHOLD only next to an empty fragment
    coding = [product[a:b] for a, b in zip([0] + cuts, cuts + [len(product)])]
    # Seeded on the product itself, not its length: two designs of the same length must not share barcodes,
    # or their arms could cross-hybridise when the oligos are pooled.
    rng, barcodes = random.Random(product), []
    taboo = {x[i:i + 8] for t in toeholds for x in (t, rc(t)) for i in range(len(x) - 7)}
    while len(barcodes) < n - 1:   # random 18-mers: balanced GC, no repeats, no shared 8-mer
        for _ in range(5000):
            b = "".join(rng.choice("ACGT") for _ in range(BARCODE))
            kmers = {x[i:i + 8] for x in (b, rc(b)) for i in range(BARCODE - 7)}
            if 8 <= sum(c in "GC" for c in b) <= 10 and not any(c * 4 in b for c in "ACGT") \
                    and not kmers & taboo and structure(b)["hairpin"] > -2:
                barcodes.append(b)
                taboo |= kmers
                break
        else:
            break               # too few: a bad check below, not an exception, so the panel still shows the design
    for i in range(n):
        d.rows.append((f"{d.names[i]}-coding", coding[i], 0, "Sidewinder oligo"))
        d.rows.append((f"{d.names[i]}-barcode", rc(coding[i]) + (barcodes[i] if i < len(barcodes) else ""), 0, "Sidewinder oligo"))
    add_primers(d, "Sidewinder product", product, "", "", target, "Sidewinder PCR")   # removes the barcode arms
    d.note = ("Each fragment is a coding strand (top) annealed to its barcode strand (bottom). A coding strand ends "
              "with the toehold it shares with the next one, and the next coding strand starts right after it.")
    for j, t in enumerate(toeholds):
        arm = barcodes[j] if j < len(barcodes) else "none found"
        d.junctions.append(f"{d.names[j]} → {d.names[j + 1]}: toehold {t}, barcode {arm}")
        top = cuts[j]            # the coding strand of fragment j ends with the toehold
        d.joints.append(Joint(d.names[j], d.names[j + 1], "toehold", t, top - TOEHOLD, (top, top, None), (top, top, None)))
        gc = 100 * sum(c in "GC" for c in t) / max(1, len(t))
        d.check(f"Toehold {j + 1} GC", band(gc, 30, 70, 10), f"{gc:.0f} %", junction=j, short="GC")
        if len(t) < TOEHOLD:
            d.check(f"Toehold {j + 1} length", "bad", f"only {len(t)} of {TOEHOLD} nt: a fragment next to it is empty or too short",
                    junction=j, short="length")
        hits = occurrences(product, t, False) if t else 0
        d.check(f"Toehold {j + 1} unique", "ok" if hits == 1 else "bad",
                "occurs once in the product" if hits == 1 else f"occurs {hits} times (could join the wrong fragments)",
                junction=j, short="unique")
        near = [k + 1 for k, u in enumerate(toeholds) if k != j and
                min(sum(x != y for x, y in zip(t, r)) for r in (u, rc(u))) < 3]
        d.check(f"Toehold {j + 1} distinct", "bad" if near else "ok",
                f"fewer than 3 bases from toehold {', '.join(map(str, near))}" if near else "3 or more bases from every other toehold",
                junction=j, short="distinct")
        d.cuts.append([f"barcode {arm}"])
    long = [r[0] for r in d.rows if r[3] == "Sidewinder oligo" and len(r[1]) > MAX_OLIGO]
    d.check(f"Oligos at most {MAX_OLIGO} nt", "bad" if long else "ok",
            "too long (split the fragment): " + ", ".join(long) if long else "every strand fits in one oligo")
    # The oligo that is ordered is the whole strand, barcode arm included: the arm can fold back on the
    # coding strand it is joined to, which screening the barcode on its own does not see.
    risky = []
    for name, seq, _, kind in d.rows:
        if kind == "Sidewinder oligo":
            warnings = structure(seq)["warnings"]
            if warnings:
                risky.append(f"{name} ({', '.join(warnings)})")
    d.check("Oligo structure", "warn" if risky else "ok",
            "folds or self-dimerises: " + "; ".join(risky) if risky else "no hairpin or self-dimer worth reporting")
    tiny = [name for name, s in zip(d.names[1:-1], d.seqs[1:-1]) if len(s) < 2 * TOEHOLD]
    d.check("Fragments hold two toeholds", "warn" if tiny else "ok", "short: " + ", ".join(tiny) if tiny else "ok")
    d.check("Barcodes orthogonal", "bad" if len(barcodes) < n - 1 else "ok",
            f"only {len(barcodes)} of {n - 1} found: use fewer fragments" if len(barcodes) < n - 1 else
            f"{len(barcodes)} × {BARCODE} nt, 8-10 G/C, no 8-mer shared with another barcode or toehold")
    d.check("Linear product", "info", "Sidewinder gives a linear product: clone it into a vector afterwards")
    check_product(d)
    return d


# ---- analysis: do the fragments, as they are, fit? ------------------------------------------------


def frag_rows(frags, kind):
    """The primers the fragments carry, as oligo rows."""
    return [(p[0], p[1], p[2], kind) for f in frags for p in (f.fwd, f.rev) if p]


def tail_intervals(d, frags):
    """d.tails: for each fragment, where in the product the bases its primers' tails added sit; d.tail_primers: whose they are. A tail
    counts where a neighbour carries the same bases (d.shared): they are still the fragment's own bases there, as when each of two
    primers adds half of a Gibson overlap. Round the origin of a circle an interval is cut in two."""
    n = sum(len(d.seqs[i]) - a - b for i, a, b in d.pieces)
    pos = 0
    for i, a, b in d.pieces:
        start = d.offsets[i] + a                                            # where its piece starts in the fragment's own bases
        end = d.offsets[i] + len(d.seqs[i]) - b                             # ... and ends
        lo, hi = start - d.shared[i][0], end + d.shared[i][1]               # the bases that are its own in the product
        full = len(frags[i].seq)
        # (the fragment is tail + core + rc(tail): the forward primer's tail is in front)
        ends = [(frags[i].fwd, "start", 0, frags[i].tails[0]), (frags[i].rev, "end", full - frags[i].tails[1], full)]
        d.tails[i], d.tail_primers[i] = [], []
        for primer, side, s, e in ends:
            s, e = max(s, lo), min(e, hi)
            if e <= s or not n:
                continue
            first, size = pos + s - start, e - s                            # in the product (before the origin wraps)
            if d.circular:
                first %= n
                parts = [(first, first + size)] if first + size <= n else [(first, n), (0, first + size - n)]
            else:
                parts = [(max(first, 0), min(first + size, n))]
            for part in parts:
                if part[1] > part[0]:
                    d.tails[i].append(part)
                    d.tail_primers[i].append((primer[0] if primer else "", side))
        pos += end - start


def analyse_gibson(frags, circular, min_ov=20, min_tm=52.0):
    """Gibson / HiFi: neighbours fit when the end of one is the start of the next (>= MIN_EXISTING bp)."""
    d = Design(frags, circular)
    d.vectors = [False] * len(frags)
    n, seqs = len(frags), d.seqs
    skip, trim, found = [0] * n, [0] * n, []
    for j in range(n if circular else n - 1):
        a, b = j, (j + 1) % n
        sa, sb = seqs[a], seqs[b]
        top = len(sa) // 2 if a == b else min(100, len(sa) - 1, len(sb) - 1)   # not the whole fragment
        ex = next((k for k in range(top, MIN_EXISTING - 1, -1) if sa[-k:] == sb[:k]), 0)
        if ex:
            if b == 0:
                trim[a] = ex     # closing the circle: drop the duplicate from the last fragment
            else:
                skip[b] = ex
        found.append(ex)
    d.pieces = [(i, skip[i], trim[i]) for i in range(n)]
    d.shared = [(skip[i], trim[i]) for i in range(n)]   # a fragment is cut where a neighbour has the same bases: they stay its own
    d.note = ("Neighbouring fragments share the same bases at their ends (one fragment's 3′ end is the next one's 5′ start). "
              "The exonuclease chews back the 5′ ends and the complementary 3′ overhangs anneal. "
              "Lowercase = a tail a primer added.")
    ends, total = [], 0
    for i in range(n):
        total += len(seqs[i]) - skip[i] - trim[i]
        ends.append(total)
    for j, ex in enumerate(found):
        a, b = j, (j + 1) % n
        edge = ends[j]                      # where fragment j's part of the product ends
        if ex:
            pos = edge if b == 0 else edge - ex
            ov = seqs[a][-ex:]
            top = pos + ex
            d.joints.append(Joint(d.names[a], d.names[b], "existing overlap", ov, pos, (top, top, None), (pos, pos, None)))
            d.problems.append(None)
            d.junctions.append(f"{d.names[a]} → {d.names[b]}: {ex} bp homology")
        else:
            d.joints.append(Joint(d.names[a], d.names[b], "none", "", edge, (edge, edge, None), (edge, edge, None)))
            d.problems.append(f"no homology: the last 15+ bases of {d.names[a]} are not the first bases of {d.names[b]}")
            d.junctions.append(f"{d.names[a]} → {d.names[b]}: no homology")
    tail_intervals(d, frags)
    product = d.product()
    ovs = [(j, jt.seq) for j, jt in enumerate(d.joints) if jt.kind != "none"]
    for j, jt in enumerate(d.joints):
        label = f"{jt.left} → {jt.right}"
        if jt.kind == "none":
            d.check(f"Homology ({label})", "bad", d.problems[j], junction=j, short="homology")
            continue
        ov, t = jt.seq, melting_temp(jt.seq)
        gc = 100 * sum(c in "GC" for c in ov) / len(ov)
        hits = occurrences(product, ov, circular)
        d.check(f"Homology ({label})", "ok" if len(ov) >= min_ov else "warn", f"{len(ov)} bp: {ov}"
                + ("" if len(ov) >= min_ov else f" (shorter than {min_ov} bp: lower efficiency)"), junction=j, short="length")
        d.check(f"Homology Tm ({label})", band(t, min_tm, 80, 4) if t else "warn", f"{t:.1f} °C" if t else "?", junction=j, short="Tm")
        d.check(f"Homology GC ({label})", band(gc, 35, 65, 10), f"{gc:.0f} %", junction=j, short="GC")
        d.check(f"Homology unique ({label})", "ok" if hits == 1 else "bad",
                "occurs once in the product" if hits == 1 else f"occurs {hits} times in the product (wrong joins)",
                junction=j, short="unique")
        d.check(f"Homology hairpin ({label})", hairpin_level(structure(ov)), f"ΔG {min(structure(ov)['hairpin'], 0):.1f}",
                junction=j, short="hairpin")
        clash = [f"{d.joints[k].left} → {d.joints[k].right}" for k, other in ovs if k != j and shared_kmer(ov, other)]
        d.check(f"Homology distinct ({label})", "warn" if clash else "ok",
                "shares a 10-mer with " + ", ".join(clash) if clash else "shares no 10-mer with another homology",
                junction=j, short="distinct")
    check_product(d)
    check_lengths(d)
    d.rows = frag_rows(frags, "Gibson")
    return d


def site_events(seq, enzyme, circular):
    """Where a Type IIS enzyme cuts `seq`: sorted [(top-strand cut, orientation)]. Orientation +1: the
    site is on the top strand and the enzyme cuts to its right; -1: on the bottom strand, cut to its
    left. The overhang is seq[cut:cut + k]; the bottom strand is cut k bases further right."""
    _, site, spacer, k = GG_ENZYMES[enzyme]
    n, L = len(seq), len(site)
    text = seq + seq[:L + spacer + k] if circular else seq
    out = set()
    for orient, pattern in ((1, site), (-1, rc(site))):
        i = text.find(pattern)
        while i != -1 and i < n:
            t = i + L + spacer if orient > 0 else i - spacer - k
            if circular:
                out.add((t % n, orient))
            elif t >= 0 and t + k <= n:
                out.add((t, orient))
            i = text.find(pattern, i + 1)
    return sorted(out)


def analyse_golden_gate(frags, circular, enzyme):
    """Golden Gate: each fragment is digested at the sites it carries (a site near its start facing in, one
    near its end facing in); neighbours fit when the overhang one ends with is the one the next starts with."""
    ename, site, spacer, k = GG_ENZYMES[enzyme]
    d = Design(frags, circular)
    d.vectors = [f.vector for f in frags]
    n = len(frags)
    njunc = n if circular else n - 1
    info = []   # per fragment: (left overhang or None, right overhang or None, internal sites)
    around = []   # per fragment: the bases (before, overhang, after) at its left end and at its right end
    for i, f in enumerate(frags):
        seq = f.seq
        if f.vector:
            ev = site_events(seq, enzyme, True)
            if len(ev) != 2:
                raise ValueError(f"{f.name}: destination vector has {len(ev)} {ename} site(s), needs exactly 2")
            (t1, o1), (t2, o2) = ev
            if (o1, o2) == (1, -1):      # the arc between the cuts holds neither site
                start, length = t1, t2 - t1
            elif (o1, o2) == (-1, 1):    # the other arc, round the origin
                start, length = t2, t1 + len(seq) - t2
            else:
                raise ValueError(f"{f.name}: the two {ename} sites point the same way; in a destination vector they face "
                                 "away from the backbone, so the enzyme removes them with the dropout")
            ring = seq * 3
            d.seqs[i], d.ranges[i] = ring[start:start + length], (start, start + length)
            info.append((ring[start:start + k], ring[start + length:start + length + k], []))
            s0, e0 = start + len(seq), start + length + len(seq)
            around.append(((ring[s0 - 10:s0].lower(), ring[s0:s0 + k], ring[s0 + k:s0 + k + 6]),
                           (ring[e0 - 6:e0], ring[e0:e0 + k], ring[e0 + k:e0 + k + 10].lower())))
            continue
        ev = site_events(seq, enzyme, False)
        tops = [e for e in ev if e[1] > 0]
        bottoms = [e for e in ev if e[1] < 0]
        left = tops[0] if tops else None
        right = bottoms[-1] if bottoms else None
        cl = left[0] if left else 0
        cr = right[0] if right else len(seq)
        if left and right and cr < cl + k:
            left = right = None             # the two sites overlap or face the wrong way: use neither
            cl, cr = 0, len(seq)
        extra = [e for e in ev if e not in (left, right)]
        d.seqs[i], d.offsets[i], d.ranges[i] = seq[cl:cr], cl, (cl, cr)
        info.append((seq[cl:cl + k] if left else None, seq[cr:cr + k] if right else None, extra))
        around.append(((seq[max(0, cl - 10):cl].lower(), seq[cl:cl + k], seq[cl + k:cl + k + 6]),
                       (seq[max(0, cr - 6):cr], seq[cr:cr + k], seq[cr + k:cr + k + 10].lower())))
    d.pieces = [(i, 0, 0) for i in range(n)]
    d.note = (f"Each fragment is digested next to its {ename} sites. Its top strand ends where the {k} nt overhang starts "
              "and its bottom strand carries the overhang (a 5′ overhang on the bottom strand of the left fragment, on the "
              "top strand of the right one); matching overhangs are ligated. Lowercase = a tail a primer added.")
    d.cuts = []
    ends, total = [], 0
    for i in range(n):
        total += len(d.seqs[i])
        ends.append(total)
    good = []
    for j in range(njunc):
        a, b = j, (j + 1) % n
        ro, lo = info[a][1], info[b][0]
        p = ends[j]
        if ro and lo and ro == lo:
            d.joints.append(Joint(d.names[a], d.names[b], "overhang", ro, p, (p, p + k, None), (p, p + k, None)))
            d.problems.append(None)
            d.junctions.append(f"{d.names[a]} → {d.names[b]}: overhang {ro}")
            good.append(ro)
        else:
            reason = (f"overhangs differ: {d.names[a]} ends with {ro}, {d.names[b]} starts with {lo}" if ro and lo else
                      f"{d.names[a]} has no {ename} site at its right end" if lo else
                      f"{d.names[b]} has no {ename} site at its left end" if ro else
                      f"no {ename} site at {d.names[a]}'s right end or at {d.names[b]}'s left end")
            d.joints.append(Joint(d.names[a], d.names[b], "none", "", p, (p, p, None), (p, p, None)))
            d.problems.append(reason)
            d.junctions.append(f"{d.names[a]} → {d.names[b]}: {reason}")
        (rp, rov, rq), (lp, lov, lq) = around[a][1], around[b][0]
        d.cuts.append([f"…{rp}|{rov}|{rq}…", f"…{lp}|{lov}|{lq}…"] if ro and lo else None)
    # A fragment's top strand stops where the overhang starts, but the overhang is still its own bases: the next piece
    # starts with them (an ORF's last codon can lie in it).
    d.shared = [(0, k if i < njunc and d.joints[i].kind != "none" else 0) for i in range(n)]
    tail_intervals(d, frags)
    for j, jt in enumerate(d.joints):
        if jt.kind == "none":
            d.check(f"Junction {j + 1}: {jt.left} → {jt.right}", "bad", d.problems[j], junction=j, short="ends fit")
            continue
        d.check(f"Junction {j + 1} ends", "ok", f"both ends carry {jt.seq}", junction=j, short="ends fit")
        others = [q.seq for i, q in enumerate(d.joints) if i != j and q.kind != "none"]
        for short, level, text in OH.junction_checks(jt.seq, others):
            d.check(f"Overhang {jt.seq}: {short}", level, text, junction=j, short=short)
    for i, (f, (lo, ro, extra)) in enumerate(zip(frags, info)):
        if not f.vector:
            d.check(f"{ename} sites in {f.name}", "bad" if extra else "ok",
                    "extra or misplaced site(s) at " + ", ".join(f"{c + 1:,}" for c, _ in extra)
                    + " (domesticate: mutate them)" if extra else "only the two end sites", fragment=i, short="extra site")
    product = d.product()
    left_over = find_sites(product, circular).get(ename, [])
    d.check(f"{ename} sites left in product", "warn" if left_over else "ok",
            "at " + ", ".join(f"{c + 1:,}" for c in left_over) + " (the product can be re-cut)" if left_over else "none")
    check_product(d)
    check_lengths(d)
    d.rows = frag_rows(frags, f"Golden Gate {ename}")
    return d


def junction_cards(d):
    """What the map shows in place of each junction: SimpleNamespace(level, title, brief, lines, chips, closing, text).
    `brief`: (what joins, "Tm · GC"), the two short lines of the compact form; the full card is `title`, `lines` and `chips`.
    `lines`: a few lines of sequence (the joining bases on both strands, or the two ends with their sites);
    `chips`: [(short name, level, text)], the checks of that junction."""
    n, cards = len(d.names), []
    for j, jt in enumerate(d.joints):
        a, b = j % n, (j + 1) % n
        chips = d.jchecks.get(j, [])
        if jt.kind == "none":
            title, brief = "✗ ends do not fit", ("ends do not fit", "")
            lines = [f"…{d.seqs[a][-14:]}", f"{d.seqs[b][:14]}…"]
        else:
            L, t = len(jt.seq), joint_tm(jt.seq)
            gc = 100 * sum(c in "GC" for c in jt.seq.upper()) / max(L, 1)
            unit, word = "nt" if jt.kind == "toehold" else "bp", JOINED_BY[jt.kind]
            title = f"{L} {unit} {word}" + (f" · {t:.0f} °C" if t is not None else "") + f" · GC {gc:.0f} %"
            # Gibson joins by homology only, so the word adds nothing there; the other methods say what they join by.
            brief = (f"{L} {unit}" + ("" if word == "homology" else f" {word}"),
                     f"{t:.0f} °C · {gc:.0f} %" if t is not None else f"GC {gc:.0f} %")
            lines = [jt.seq, "|" * L, jt.seq.translate(COMP)]
            if j < len(d.cuts) and d.cuts[j]:
                lines = list(d.cuts[j]) if jt.kind == "overhang" else lines + list(d.cuts[j])
        tip = "\n".join(f"{MARKS[level][0]} {name}: {text}" for name, level, text in chips)
        cards.append(SimpleNamespace(level=d.jlevels[j] if j < len(d.jlevels) else "ok", title=title, brief=brief, lines=lines, chips=chips,
                                     closing=d.circular and j == n - 1 and n > 1 or (n == 1 and d.circular),
                                     text=(d.problems[j] if jt.kind == "none" and j < len(d.problems) and d.problems[j] else tip)))
    return cards


def analyse(frags, method, circular, enzyme="BsaI", min_ov=20, min_tm=52.0, tm=60.0):
    """Check `frags` as they are (0 Gibson, 1 Golden Gate, 2 Sidewinder). Raises ValueError for a setup that
    cannot be analysed at all (e.g. a destination vector without two sites)."""
    if method == 2:
        if len(frags) < 2:
            raise ValueError("Sidewinder needs at least two fragments")
        d = design_sidewinder(frags, tm)
        d.problems = [None] * len(d.joints)
        tail_intervals(d, frags)
    else:
        d = analyse_golden_gate(frags, circular, enzyme) if method == 1 else analyse_gibson(frags, circular, min_ov, min_tm)
    # One verdict per junction: the worst of its checks (an end that does not fit is always red).
    d.jlevels = ["bad" if jt.kind == "none" else {"info": "ok"}.get(worst(d.jchecks.get(j, [])), worst(d.jchecks.get(j, [])))
                 for j, jt in enumerate(d.joints)]
    return d


def auto_primers(frags, method, circular, enzyme="BsaI", min_ov=20, max_ov=40, min_tm=52.0, tm=60.0):
    """Give every fragment new primers carrying what the method needs at each junction (a homology
    tail for Gibson; site, spacer and overhang for Golden Gate), designed on the fragments' cores, so
    primers chosen earlier are replaced. Returns how many fragments got primers."""
    views = [SimpleNamespace(name=f.name, seq=f.core_seq, vector=f.vector, flipped=f.flipped, doc=f.core) for f in frags]
    d = design_gibson(views, circular, min_ov, max_ov, min_tm, tm) if method == 0 else design_golden_gate(views, circular, enzyme, tm)
    rows = {r[0]: r for r in d.rows}
    count = 0
    for f in frags:
        fw, rv = rows.get(f"{f.name}-F"), rows.get(f"{f.name}-R")
        if not f.vector and fw and rv:
            f.fwd, f.rev = (fw[0], fw[1].upper(), fw[2]), (rv[0], rv[1].upper(), rv[2])
            f.rebuild()
            count += 1
    return count


# ---- product ------------------------------------------------------------------------------------


def wrapped(feature, length):
    """The location of `feature` when the part of it past base `length` carries on at the start of a circle that long.
    Written again from its base ranges (a fuzzy `<` or `>` mark is lost, which only matters for a feature that wraps)."""
    parts = []
    for s, e in feature.spans:
        if e <= length:
            parts.append((s, e))
        elif s >= length:
            parts.append((s - length, e - length))
        else:
            parts += [(s, length), (0, e - length)]
    where = ",".join(f"{s + 1}..{e}" for s, e in parts)
    where = f"join({where})" if len(parts) > 1 else where
    return f"complement({where})" if feature.strand < 0 else where


def placed_features(doc, window, at, length, circular):
    """The features of `doc` over its bases [window[0], window[1]) as (kind, location, qualifiers), clipped to the window,
    for a product of `length` bases where window[0] sits at product position `at`. On a circular product what runs past
    its end carries on at its start (the bases that close the circle are the first bases of the first fragment). The core
    can only move features to the right (a paste does it), so they go through a scratch PCR product with `at` filler
    bases in front of the window."""
    lo, hi = window
    if at < 0:
        lo, at = lo - at, 0           # the window starts before the product does: nothing of it can be placed
    hi = min(hi, lo + (length if circular else length - at), lo + doc.length if doc.circular else doc.length)
    if hi <= lo or lo >= doc.length:
        return []
    scratch = doc.pcr_product("scratch", lo, hi, "N" * at, "")
    return [(x.kind, wrapped(x, length) if circular and any(e > length for _, e in x.spans) else x.location, x.qualifiers)
            for x in scratch.features()]


def transfer_primers(doc, frags):
    """The primers that made the fragments (typed, picked in a file or designed) as primer_bind features of the product, at the
    sites where they bind it: the product shows them like any primer of the program, and a file saved from it stores them. A
    primer is found by its bases, the way the sequence view finds every primer: a Gibson tail is the neighbour's own end, so
    the whole primer binds; a Golden Gate site and spacer are cut off, so they hang off the product as a tail."""
    primers = []
    for f in frags:
        for p in (f.fwd, f.rev):
            if p and (p[0], p[1].upper()) not in primers:
                primers.append((p[0], p[1].upper()))
    if not primers or not doc.length:
        return
    have = {(f.label, f.location) for f in doc.features() if f.kind == "primer_bind"}   # (one a fragment's file already had)
    for (name, seq), found in zip(primers, find_primers(doc.sequence(), doc.circular, [s for _, s in primers], 0)):
        for s, e, strand, _ in found[:MAX_SITES]:
            where = site_location(s, e, strand, doc.length)
            if (name, where) not in have:
                doc.set_feature(None, "primer_bind", where, [("label", name), ("note", f"sequence: {seq}")])


def simulate(frags, d, preview=False, tails=True):
    """The assembled product as an unsaved Document: each fragment's features moved into place and the primers that made the
    fragments (transfer_primers), nothing more, so the file Perform assembly opens carries only what the sequences had. `preview`:
    also what the panel draws to tell the parts apart, one misc_feature per fragment named after it in the fragment's colour and
    one per junction (the homology, overhang or toehold) in slate; `tails`: and, in a preview, one per primer tail in orange (the
    junction inspector shows them, the product view of the panel does not)."""
    doc = Document.new_sequence("Assembly product", "", False)
    spans = []
    for i, a, b in d.pieces:
        seq = d.seqs[i]
        start = doc.length
        doc.paste(start, start, seq[a:len(seq) - b], [])
        spans.append((d.names[i], start, doc.length))
    doc.circular = d.circular
    seen = set()    # both fragments of an overlap can carry the same annotation of it: the product shows it once
    for (i, a, b), (_, start, _) in zip(d.pieces, spans):
        rng, (before, after) = d.ranges[i] or (0, len(d.seqs[i])), d.shared[i]
        feats = placed_features(frags[i].doc, (rng[0] + a - before, rng[1] - b + after), start - before, doc.length, d.circular)
        new = [x for x in feats if (x[0], x[1], tuple(x[2])) not in seen]
        seen.update((x[0], x[1], tuple(x[2])) for x in new)
        if new:
            doc.add_features(new)
    if preview:
        for (name, s, e), (i, _, _) in zip(spans, d.pieces):
            doc.set_feature(None, "misc_feature", f"{s + 1}..{e}", [("label", name), ("ApEinfo_fwdcolor", FRAG_COLORS[i % len(FRAG_COLORS)])])
    transfer_primers(doc, frags)
    if preview:
        n = doc.length
        for j in d.joints:
            if j.kind == "none" or not j.seq or n == 0:
                continue
            s, e = j.pos % n, j.pos % n + len(j.seq)
            what = {"toehold": "toehold", "overhang": "overhang"}.get(j.kind, "homology")
            where = f"{s + 1}..{e}" if e <= n else f"join({s + 1}..{n},1..{e - n})"   # the closing junction wraps round the origin
            doc.set_feature(None, "misc_feature", where, [("label", f"{what} {j.seq if len(j.seq) <= 8 else str(len(j.seq)) + ' bp'}"),
                                                          ("note", f"joins {j.left} to {j.right}"), ("ApEinfo_fwdcolor", JOINT_COLOR)])
        for i, (spans, whose) in enumerate(zip(d.tails, d.tail_primers) if tails else ()):   # what the primers added to a fragment (orange)
            for (s, e), (name, side) in zip(spans, whose):
                if 0 <= s < e <= n:
                    doc.set_feature(None, "misc_feature", f"{s + 1}..{e}",
                                    [("label", f"tail {e - s} nt"),
                                     ("note", f"{e - s} nt that primer {name or '?'} adds to the {side} of {d.names[i]}: not in its own sequence"),
                                     ("ApEinfo_fwdcolor", TAIL_COLOR)])
    return doc


def junction_facts(d, j):
    """What the bases that join junction `j` are made of: SimpleNamespace(length, added, shared). `added`: [(primer, fragment index,
    "start" | "end" of that fragment, bases)], the bases of the joint a primer's tail puts there; `shared`: the rest, which both
    neighbours already had."""
    jt = d.joints[j]
    n, length = len(d.product()), len(jt.seq)
    spans = []
    if n and length:
        s = jt.pos % n
        spans = [(s, s + length)] if s + length <= n else [(s, n), (0, s + length - n)]   # the closing junction wraps round the origin
    added = []
    for i, (tails, whose) in enumerate(zip(d.tails, d.tail_primers)):
        for (s, e), (name, side) in zip(tails, whose):
            count = sum(max(0, min(e, b) - max(s, a)) for a, b in spans)
            if count:
                added.append((name, i, side, count))
    return SimpleNamespace(length=length, added=added, shared=max(0, length - sum(c for *_, c in added)))


def junction_summary(d, j):
    """One or two sentences on junction `j` for the inspector: what joins the two fragments, how much of it a primer adds and how much
    they already shared (the part that was easy to lose track of in the old text)."""
    jt, n = d.joints[j], len(d.names)
    a, b = j % n, (j + 1) % n
    left, right = f"Fragment {a + 1}", f"Fragment {b + 1}"
    if jt.kind == "none":
        return d.problems[j] if j < len(d.problems) and d.problems[j] else f"The ends of {left} and {right} do not fit."
    facts = junction_facts(d, j)
    unit, word = ("nt", "toehold") if jt.kind == "toehold" else ("bp", JOINED_BY[jt.kind])
    what = f"The {facts.length} {unit} {word} between {left} and {right}"
    if not facts.added:
        return f"{what} comes from the fragments' own bases: no primer adds anything here."
    parts = [f"primer {name or '?'} adds {count} nt to the {side} of Fragment {i + 1}" for name, i, side, count in facts.added]
    rest = f"; the other {facts.shared} {unit} were already shared by the two fragments" if facts.shared else ""
    return f"{what}: " + " and ".join(parts) + rest + "."


def junction_window(doc, d, j, context=50):
    """(a Document of the product around junction `j`, where the joint starts in it, where the document starts in the product): the
    joint and `context` bases each side, round the origin for the junction that closes a circle, with the features and primers there
    (a primer reaches about 60 bases from it)."""
    jt, n = d.joints[j], doc.length
    length = max(len(jt.seq), 1)
    start = jt.pos % n
    if doc.circular:
        span = length + 2 * context
        if span >= n:   # a product no longer than that: all of it
            return snapshot(doc, f"Junction {j + 1}"), start, 0
        first = (start - context) % n
        return doc.pcr_product(f"Junction {j + 1}", first, first + span, "", ""), context, first
    lo, hi = max(0, start - context), min(n, start + length + context)
    return doc.pcr_product(f"Junction {j + 1}", lo, hi, "", ""), start - lo, lo


def piece_doc(frags, d, i):
    """Fragment i as the product uses it (a vector: its backbone, a Golden Gate insert: digested) as an unsaved Document."""
    seq, rng, f = d.seqs[i], d.ranges[i] or (0, len(d.seqs[i])), frags[i]
    doc = Document.new_sequence(f"{f.name} piece", "", False)
    doc.paste(0, 0, seq, [(x.kind, x.location, x.qualifiers) for x in f.doc.copy_features(*rng)])
    return doc
