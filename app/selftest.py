"""Self-test: small checks that the editing, annotation and primer code is right, run from
Help → Troubleshooting → Self-Test (and by pytest). Each check builds its own in-memory Document,
so it never touches the user's files. Add one with @check("Area", "what it verifies"); it
passes when it returns, fails when it raises AssertionError. The full plan is docs/TESTING.md.
"""
import random
import time

from PySide6.QtCore import Qt
from PySide6.QtGui import QColor
from PySide6.QtWidgets import QHBoxLayout, QLabel, QPushButton, QTableWidget, QTableWidgetItem, QVBoxLayout, QWidget

from opendna_core import Document
from primer_panel import MARKS
from uikit import fit_columns, right_aligned, run_busy

CHECKS = []  # (area, what, function), in registration order


def check(area, what):
    def register(fn):
        CHECKS.append((area, what, fn))
        return fn
    return register


def run_check(fn):
    """(ok, milliseconds, message)."""
    t = time.perf_counter()
    try:
        fn()
        ok, message = True, ""
    except AssertionError as e:
        ok, message = False, str(e) or "assertion failed"
    except Exception as e:  # a crash is a failure too
        ok, message = False, f"{type(e).__name__}: {e}"
    return ok, (time.perf_counter() - t) * 1000, message


def dna(n, seed=1):
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def revcomp(s):
    return s[::-1].translate(str.maketrans("ACGT", "TGCA"))


def covered(doc, label):
    """The bases a feature covers, in order."""
    f = next(f for f in doc.features() if f.label == label)
    return "".join(doc.sequence(s, e) for s, e in f.spans)


@check("Editing", "insert, delete and replace at the start, end and whole sequence match plain string edits")
def _():
    seq = dna(60)
    for op, args in [("insert", (0, "AC")), ("insert", (60, "GT")), ("delete", (0, 5)), ("delete", (55, 60)),
                     ("delete", (0, 60)), ("replace", (0, 60, "TTT")), ("replace", (10, 10, "G")), ("insert", (30, ""))]:
        doc = Document.new_sequence("t", seq, False)
        getattr(doc, op)(*args)
        want = {"insert": lambda p, b: seq[:p] + b + seq[p:], "delete": lambda s, e: seq[:s] + seq[e:],
                "replace": lambda s, e, b: seq[:s] + b + seq[e:]}[op](*args)
        assert doc.sequence() == want, f"{op}{args}: got {doc.sequence()!r}"


@check("Editing", "200 random edits then undo all restores the original bases and features exactly")
def _():
    rng = random.Random(7)
    doc = Document.new_sequence("t", dna(300), True)
    doc.add_features([("gene", "20..80", [("label", "a")]), ("CDS", "complement(join(100..120,150..170))", [("label", "b")])])
    before = (doc.sequence(), [(f.location, f.qualifiers) for f in doc.features()])
    steps = len(doc.history()[0])  # adding the features was a step too: undo only the edits
    for _ in range(200):
        s = rng.randint(0, doc.length)
        doc.replace(s, min(doc.length, s + rng.randint(0, 4)), "".join(rng.choice("ACGT") for _ in range(rng.randint(0, 4))))
    while len(doc.history()[0]) > steps:
        assert doc.undo()
    assert (doc.sequence(), [(f.location, f.qualifiers) for f in doc.features()]) == before, "undo did not restore the document"


def feature_doc():
    """40 bases, a 60-base marked feature (1-based 41..100), 40 bases."""
    seq = dna(40, 2) + dna(60, 3) + dna(40, 4)
    doc = Document.new_sequence("t", seq, False)
    doc.add_features([("misc_feature", "41..100", [("label", "m")])])
    return doc, seq[40:100]


@check("Annotation drift", "an insertion before, after or inside a feature leaves it covering its own bases (inside: grown)")
def _():
    doc, own = feature_doc()
    doc.insert(10, "AAAAA")
    assert covered(doc, "m") == own, "insert before moved the feature off its bases"
    doc.insert(doc.length - 10, "CCC")
    assert covered(doc, "m") == own, "insert after changed the feature"
    doc.insert(45 + 20, "GGG")  # inside: 20 bases into the feature (it now starts at 45)
    assert covered(doc, "m") == own[:20] + "GGG" + own[20:], "insert inside did not grow the feature"


@check("Annotation drift", "100 random edits outside a feature never change the bases it covers")
def _():
    rng = random.Random(3)
    doc, own = feature_doc()
    for _ in range(100):
        s, e = doc.features()[0].spans[0]
        lo, hi = rng.choice([(0, s), (e, doc.length)])
        pos = rng.randint(lo, hi)
        if rng.random() < 0.5:
            doc.insert(pos, dna(rng.randint(1, 5), rng.randint(0, 99)))
        else:
            doc.delete(pos, min(hi, pos + 3))
        assert covered(doc, "m") == own, "feature drifted"


# ---- the two property tests' machinery (also used by tests/test_editing.py) -------------------

def snapshot(doc):
    """Everything a user could see of a document: bases, topology, every feature exactly."""
    return (doc.sequence(), doc.circular, [(f.kind, f.location, f.qualifiers) for f in doc.features()])


class IdModel:
    """The 'every base has an identity' model of what an edit must do to annotations.

    Each base carries an id that never changes (bases the user adds get new ids). A feature
    is the list of ids it covers, part by part. An edit is mirrored on the id list:
      * ids in a deleted stretch leave every feature (a feature left with none is gone);
      * ids of inserted bases join a feature only when the insertion falls strictly inside
        one of its parts (between two bases of the same part): at a feature's very start
        or end, or between two parts, they stay outside;
      * rotate / reverse complement keep every feature's own sequence (same ids in reading
        order) and flip the strand for reverse complement;
      * a pasted feature covers exactly the pasted copies of the bases it covered.
    `check` then asserts the real Document agrees, feature by feature."""

    def __init__(self, doc):
        self.ids = list(range(doc.length))
        self.next_id = doc.length
        self.parts = [[self.ids[s:e] for s, e in f.spans] for f in doc.features()]
        self.meta = [(f.kind, f.qualifiers, f.strand) for f in doc.features()]

    def fresh(self, count):
        ids = list(range(self.next_id, self.next_id + count))
        self.next_id += count
        return ids

    @staticmethod
    def bio(parts, strand):
        flat = [i for p in parts for i in p]
        return flat if strand > 0 else flat[::-1]

    def replace(self, start, end, new_ids):
        """Bases [start, end) replaced by `new_ids` (delete, then insert at `start`)."""
        deleted = set(self.ids[start:end])
        kept = self.ids[:start] + self.ids[end:]
        left = kept[start - 1] if start > 0 else None
        right = kept[start] if start < len(kept) else None
        parts_out, meta_out = [], []
        for parts, meta in zip(self.parts, self.meta):
            parts = [[i for i in p if i not in deleted] for p in parts]
            parts = [p for p in parts if p]
            if not parts:
                continue  # every base it covered is gone: so is the feature
            if new_ids and left is not None and right is not None:
                for p in parts:
                    if left in p and right in p:  # strictly inside this part
                        k = p.index(right)
                        p[k:k] = new_ids
                        break
            parts_out.append(parts)
            meta_out.append(meta)
        self.ids = kept[:start] + list(new_ids) + kept[start:]
        self.parts, self.meta = parts_out, meta_out

    def copies(self, a, b):
        """Fresh ids for copies of bases [a, b), and the features a copy of them carries:
        [(set of copied ids it covers, its meta)]."""
        copied = self.ids[a:b]
        new_ids = self.fresh(len(copied))
        twin = dict(zip(copied, new_ids))
        carried = []
        for parts, m in zip(self.parts, self.meta):
            covered = {twin[i] for p in parts for i in p if i in twin}
            if covered:
                carried.append((covered, m))
        return new_ids, carried

    def add_pasted(self, carried):
        for covered, m in carried:
            self.parts.append(("set", covered))
            self.meta.append(m)

    def rotate(self, k):
        self.parts = [("bio", self.bio(p, m[2]), m[2]) for p, m in zip(self.parts, self.meta)]
        self.ids = self.ids[k:] + self.ids[:k]

    def reverse_complement(self):
        self.parts = [("bio", self.bio(p, m[2]), -m[2]) for p, m in zip(self.parts, self.meta)]
        self.meta = [(k, q, -s) for k, q, s in self.meta]
        self.ids = self.ids[::-1]

    def add_feature(self, kind, qualifiers, strand, start, end):
        self.parts.append([self.ids[start:end]])
        self.meta.append((kind, qualifiers, strand))

    def delete_feature(self, index):
        del self.parts[index], self.meta[index]

    def check(self, doc, where):
        """Assert `doc` matches the model; `where` says which edit just happened."""
        assert doc.length == len(self.ids), f"{where}: length {doc.length}, model {len(self.ids)}"
        features = doc.features()
        assert len(features) == len(self.parts), (
            f"{where}: {len(features)} features, model has {len(self.parts)} "
            f"(document has {[f.label for f in features]})")
        for i, (f, want, meta) in enumerate(zip(features, self.parts, self.meta)):
            got = [self.ids[s:e] for s, e in f.spans]
            if isinstance(want, tuple):  # pending: pasted, rotated or reverse-complemented
                if want[0] == "set":
                    assert {x for p in got for x in p} == want[1], f"{where}: pasted feature #{i} covers the wrong bases"
                else:
                    assert f.strand == want[2], f"{where}: feature #{i} {f.label}: strand {f.strand}, want {want[2]}"
                    assert self.bio(got, f.strand) == want[1], (
                        f"{where}: feature #{i} {f.label} {f.location}: its own sequence changed "
                        f"(ids now {self.bio(got, f.strand)[:12]}..., want {want[1][:12]}...)")
                self.parts[i] = got  # from here on, follow the real structure
            else:
                assert got == want, (
                    f"{where}: feature #{i} {f.label} {f.location} covers ids "
                    f"{[(p[0], p[-1], len(p)) for p in got]}, want {[(p[0], p[-1], len(p)) for p in want]}")
            kind, qualifiers, strand = meta
            assert (f.kind, f.qualifiers, f.strand) == (kind, qualifiers, strand), (
                f"{where}: feature #{i} {f.label}: kind/qualifiers/strand changed")


# Feature shapes on a 300 bp sequence: simple, complement, join with a gap, complement(join),
# fuzzy ends, one base, a whole-sequence feature, overlapping features, touching parts.
SHAPES = [("gene", "21..60", "a"), ("gene", "complement(81..120)", "b"), ("CDS", "join(131..150,171..190)", "c"),
          ("CDS", "complement(join(201..220,231..250))", "d"), ("misc_feature", "<261..>280", "e"),
          ("misc_feature", "290..290", "one"), ("misc_feature", "1..300", "whole"), ("misc_feature", "41..100", "overlap"),
          ("CDS", "join(141..160,161..170)", "touch"), ("gene", "complement(<5..>15)", "fuzzy_c")]
ORIGIN_SHAPES = [("rep_origin", "join(291..300,1..10)", "ori"), ("CDS", "complement(join(286..300,1..5))", "cori")]


def shape_doc(seed, circular):
    """A 300 bp document with the feature shapes above (the origin-crossing ones too if circular)."""
    doc = Document.new_sequence("t", dna(300, seed), circular)
    shapes = SHAPES + (ORIGIN_SHAPES if circular else [])
    doc.add_features([(k, loc, [("label", label), ("note", f"note of {label}")]) for k, loc, label in shapes])
    return doc


def rand_bases(rng, count):
    return "".join(rng.choice("ACGT") for _ in range(count))


def random_edit(rng, doc, model, string):
    """Make one random edit on `doc`, mirror it on `model` (an IdModel) and on `string` (a
    one-item list holding the plain bases). Returns a text saying what was done."""
    n = doc.length
    edges = [0, n] + [x for f in doc.features() for sp in f.spans for x in sp]

    def pos():  # half of the time at or next to a feature boundary: that is where bugs live
        if rng.random() < 0.5:
            return min(n, max(0, rng.choice(edges) + rng.choice((-1, 0, 0, 1))))
        return rng.randint(0, n)

    kinds = ["replace", "insert", "delete", "overwrite", "paste", "add feature", "delete feature", "reverse complement"]
    weights = [3, 4, 4, 1, 3, 1, 1, 0.3]
    if doc.circular and n >= 12:
        kinds += ["across", "across with features", "rotate"]
        weights += [1.5, 1, 1]
    if n < 60:
        kinds, weights = ["insert", "paste"], [1, 1]
    kind = rng.choices(kinds, weights)[0]
    if kind == "delete feature" and not doc.features():
        kind = "add feature"
    seq = string[0]

    if kind in ("replace", "insert", "delete"):
        s = pos() if kind == "insert" else min(pos(), n - 1)
        e = min(n, s + (0 if kind == "insert" else rng.randint(0 if kind == "replace" else 1, 6)))
        new = "" if kind == "delete" else rand_bases(rng, rng.randint(1 if kind == "insert" else 0, 6))
        if s == e and not new:
            new = "A"  # an empty edit makes no undo step; keep every step a real one
        what = f"{kind}({s}, {e}, {new!r})"
        doc.replace(s, e, new)
        model.replace(s, e, model.fresh(len(new)))
        string[0] = seq[:s] + new + seq[e:]
    elif kind == "overwrite":
        s = min(pos(), n - 1)
        new = rand_bases(rng, min(rng.randint(1, 5), n - s))
        what = f"overwrite({s}, {new!r})"
        doc.overwrite(s, new)
        string[0] = seq[:s] + new + seq[s + len(new):]
    elif kind == "paste":
        a = rng.randint(0, n - 1)
        b = min(n, a + rng.randint(1, 40))
        s = pos()
        e = min(n, s + rng.choice((0, 0, 1, 5)))
        copied = [(f.kind, f.location, f.qualifiers) for f in doc.copy_features(a, b)]
        what = f"paste({s}, {e}) of bases {a}..{b} with {len(copied)} features"
        text = doc.sequence(a, b)
        doc.paste(s, e, text, copied)
        new_ids, carried = model.copies(a, b)
        model.replace(s, e, new_ids)
        model.add_pasted(carried)
        string[0] = seq[:s] + text + seq[e:]
    elif kind == "add feature":
        s = rng.randint(0, n - 1)
        e = rng.randint(s + 1, min(n, s + 50))
        strand = rng.choice((1, -1))
        loc = f"{s + 1}..{e}" if strand > 0 else f"complement({s + 1}..{e})"
        quals = [("label", f"new{rng.randint(0, 999)}")]
        what = f"add feature {loc}"
        doc.add_features([("misc_feature", loc, quals)])
        model.add_feature("misc_feature", quals, strand, s, e)
    elif kind == "delete feature":
        feats = doc.features()
        if not feats:
            return "nothing to delete"
        i = rng.randrange(len(feats))
        what = f"delete feature #{i} {feats[i].label}"
        doc.delete_feature(i)
        model.delete_feature(i)
    elif kind == "reverse complement":
        what = "reverse complement"
        doc.reverse_complement()
        model.reverse_complement()
        string[0] = revcomp(seq)
    elif kind == "rotate":
        k = rng.randint(1, n - 1)
        what = f"rotate({k})"
        doc.rotate(k)
        model.rotate(k)
        string[0] = seq[k:] + seq[:k]
    else:  # across the origin (circular only): bases start..n and 0..head are replaced
        start = n - rng.randint(1, 8)
        end = n + rng.randint(1, 8)
        head = end - n
        a = rng.randint(0, n - 1)
        b = min(n, a + rng.randint(0, 5))
        if kind == "across":
            b = a  # plain bases, no features
        text = doc.sequence(a, b) if kind == "across with features" else rand_bases(rng, rng.randint(0, 5))
        copied = [(f.kind, f.location, f.qualifiers) for f in doc.copy_features(a, b)] if kind == "across with features" else None
        what = f"edit_across_origin({start}, {end}, {len(text)} bases, {len(copied) if copied is not None else 'no'} features)"
        doc.edit_across_origin(start, end, text, copied)
        if kind == "across with features":
            new_ids, carried = model.copies(a, b)
        else:
            new_ids, carried = model.fresh(len(text)), []
        model.replace(0, head, [])
        model.replace(start - head, n - head, new_ids)
        model.add_pasted(carried)
        string[0] = seq[head:start] + text
    return what


def run_random_edits(seed, circular, steps):
    """The two property tests in one run: `steps` random edits on a document with many feature
    shapes. After every edit the bases must equal a plain string model and every feature must
    cover exactly the bases (ids) the per-base identity model says. Then undo walks back
    through every snapshot taken on the way, and redo forward again; each must be identical.
    Returns the number of undo steps made."""
    rng = random.Random(seed)
    doc = shape_doc(seed, circular)
    base_depth = len(doc.history()[0])
    model, string = IdModel(doc), [doc.sequence()]
    snaps = {base_depth: snapshot(doc)}
    for step in range(steps):
        before = len(doc.history()[0])
        what = random_edit(rng, doc, model, string)
        where = f"seed {seed} {'circular' if circular else 'linear'} edit {step + 1}: {what}"
        assert doc.sequence() == string[0], f"{where}: bases differ from the plain string model"
        assert doc.circular == circular, f"{where}: topology changed"
        assert len(doc.history()[0]) == before + 1, f"{where}: not exactly one undo step"
        model.check(doc, where)
        snaps[len(doc.history()[0])] = snapshot(doc)
    top = len(doc.history()[0])
    for depth in range(top - 1, base_depth - 1, -1):
        assert doc.undo(), f"seed {seed}: undo failed at depth {depth + 1}"
        assert snapshot(doc) == snaps[depth], f"seed {seed}: undoing to step {depth} did not restore the document exactly"
    for depth in range(base_depth + 1, top + 1):
        assert doc.redo(), f"seed {seed}: redo failed at step {depth}"
        assert snapshot(doc) == snaps[depth], f"seed {seed}: redo to step {depth} differs from the edited document"
    assert not doc.redo()
    return top - base_depth


@check("Editing", "random edits of every kind on a plasmid and a linear sequence match a plain string; undo / redo walk back and forth exactly")
def _():
    for circular in (False, True):
        run_random_edits(1, circular, 150)


@check("Annotation drift", "after each random edit every feature (join, complement, fuzzy, across the origin) covers exactly its own bases")
def _():
    for circular in (False, True):
        run_random_edits(2, circular, 150)


@check("Editing", "making a plasmid linear at a base is one undo step; undo gives back the circle and its origin")
def _():
    doc = shape_doc(3, True)
    before, steps = snapshot(doc), len(doc.history()[0])
    seq = doc.sequence()
    doc.linearize(100)
    assert doc.sequence() == seq[100:] + seq[:100] and not doc.circular, "wrong bases or still circular"
    assert len(doc.history()[0]) == steps + 1, "not one undo step"
    assert doc.undo() and snapshot(doc) == before, "undo did not give back the circle"
    doc.set_topology(False)
    assert not doc.circular and doc.undo() and doc.circular, "the topology toggle is not undoable"


PRIMER = dna(24, 5)


def primer_doc():
    return Document.new_sequence("t", dna(50, 6) + PRIMER + dna(50, 7), False)


@check("Primer binding", "a mismatch is counted once and the site is found on both strands at the same place")
def _():
    doc = primer_doc()
    odd = PRIMER[:10] + ("A" if PRIMER[10] != "A" else "C") + PRIMER[11:]
    for seq, strand in [(odd, 1), (revcomp(odd), -1)]:
        assert doc.binding_sites(seq, 0) == [], "mismatch primer found with 0 mismatches allowed"
        assert doc.binding_sites(seq, 1) == [(50, 74, strand, 1)], doc.binding_sites(seq, 1)
    assert doc.binding_sites(PRIMER, 1) == [(50, 74, 1, 0)]


@check("Primer binding", "a 5' tail is not counted as mismatches, and a mismatch at the 3' end never binds")
def _():
    doc = primer_doc()
    (start, end, strand, mm), = doc.binding_sites("GAATTC" + PRIMER, 0)  # a tail base may match by chance: the start can move left
    assert (end, strand, mm) == (74, 1, 0) and start <= 50, "tail changed the site or the count"
    bad = PRIMER[:-1] + ("A" if PRIMER[-1] != "A" else "C")
    assert doc.binding_sites(bad, 2) == [], "primer with a 3' mismatch bound"


class SelfTestPanel(QWidget):
    """A table of every check: Run all fills in ✓/✗, time and the failure message."""

    def __init__(self):
        super().__init__()
        self.table = QTableWidget(len(CHECKS), 5)
        self.table.setHorizontalHeaderLabels(["Area", "What it verifies", "", "ms", "Message"])
        self.table.setEditTriggers(QTableWidget.NoEditTriggers)
        fit_columns(self.table, stretch=4, cap=480)   # (the message takes the room left; "what it verifies" read whole, within reason)
        for row, (area, what, _) in enumerate(CHECKS):
            self.table.setItem(row, 0, QTableWidgetItem(area))
            self.table.setItem(row, 1, QTableWidgetItem(what))
        self.status = QLabel("Not run")
        run = QPushButton("Run all")
        run.clicked.connect(self.run_all)
        bar = QHBoxLayout()
        bar.addWidget(run)
        bar.addWidget(self.status, 1)
        layout = QVBoxLayout(self)
        layout.addWidget(self.table)
        layout.addLayout(bar)

    def run_all(self):
        # Every check in a thread, behind a busy dialog if they take a while (the window stays alive); then the table is filled.
        self.status.setText(f"Running {len(CHECKS)} checks…")
        results = run_busy(self, "Running the self-test…", lambda: [run_check(fn) for _, _, fn in CHECKS])
        failed = 0
        for row, (ok, ms, message) in enumerate(results):
            failed += not ok
            mark, color = MARKS["ok" if ok else "bad"]
            cell = QTableWidgetItem(mark)
            cell.setForeground(QColor(color))
            cell.setTextAlignment(Qt.AlignCenter)
            self.table.setItem(row, 2, cell)
            self.table.setItem(row, 3, right_aligned(QTableWidgetItem(f"{ms:.0f}")))
            self.table.setItem(row, 4, QTableWidgetItem(message))
        self.status.setText(f"{len(CHECKS) - failed} passed, {failed} failed")
