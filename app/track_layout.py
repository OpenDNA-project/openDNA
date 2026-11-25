"""The order, groups and look of the Alignment view's tracks for one open file (`Editor.track_layout`), kept apart from Qt so
it can be tested and saved in a project (.oDNA).

A row is a line of the view, named by a key: a built-in track ("GC %", "Bases", "Features", "Primers", "Query hits"), a
loaded track (its name), or one of a BAM file's three rows (its name, "<name> variants", "<name> reads"): the Editor keeps the
table of rows (`Editor.row_table`). The layout is a list of nodes, each a row key or a `Group` of row keys (one level: a group
holds rows, never groups). Each row has a `Style`. Every change made through a method here calls `on_change(edited)` (the
Editor's: the views and the Tracks panel follow it, and a change made by hand, `edited`, is something new for a project to save;
adding and removing rows and loading a saved layout are not).

The layout also holds the highlighted regions (`regions`, a list of `Region`): bands drawn across every track, listed in the
panel's Regions tab, recapped and exported as BED, saved with the rest.
"""
import re
from collections import Counter
from dataclasses import asdict, dataclass, field

BUILTIN = {"Bases": "bases", "Features": "features", "Primers": "primers", "Query hits": "hits", "GC %": "gc"}   # key -> kind (bases at the top)
FROZEN = ("Bases", "Features")             # the rows kept in view when the tracks scroll, until the user unfreezes them
HIDDEN = ("GC %",)                         # the rows hidden until the user shows them
SEQUENCE = "Sequence"     # the group the built-in tracks start in
MODES = ("expanded", "condensed", "dense")   # dense: one line (UCSC); hiding is `Style.shown`
SCALES = ("linear", "log2", "log10")
SCALED = ("coverage", "signal")   # kinds with a scale choice (bars: an annotation track has none to choose)
COLORED = ("gc", "track", "coverage", "signal")   # kinds drawn in one colour that can be chosen (features, primers, reads have their own)
RENDERS = ("bars", "line", "heatmap")      # how a signal track is drawn
SUMMARIES = ("mean", "max")                # what a signal shows where a pixel holds many bases: their mean, or the highest value
RANGES = ("auto", "group", "fixed")        # a signal's range: its own values in view, its group's, or low..high
FIXED = ("bases",)                         # kinds that are only shown or hidden (no display mode, height or colour)
REGION_COLOR = "#ffb300"                   # a highlight's colour unless one is chosen
SPLIT = re.compile(r"[_\-. ]+")   # what separates the parts of a track's name ("Group by name")


@dataclass
class Style:
    shown: bool = True
    mode: str = "expanded"
    height: int = None        # px, None: the mode's own
    scale: str = "linear"     # bar tracks: one of SCALES
    color: str = None         # "#rrggbb" chosen by hand (over the file's colours), None: the track's own
    label: str = None         # the name shown, None: the key
    render: str = "bars"      # signal tracks: one of RENDERS
    summary: str = "mean"     # signal tracks: one of SUMMARIES
    range: str = "auto"       # signal tracks: one of RANGES
    low: float = 0.0          # signal tracks: the fixed range
    high: float = 1.0
    overlay: str = None       # signal tracks: the key of the signal row it is drawn in (an overlay), None: a row of its own
    frozen: bool = False      # kept in view at the top when the tracks scroll


DEFAULT = asdict(Style())


def default_style(key):
    """A row's style before anything is chosen: the sequence and the annotations frozen, the GC % hidden."""
    return Style(shown=key not in HIDDEN, frozen=key in FROZEN)


@dataclass
class Region:
    """A highlighted region: bases start..end (0-based, end excluded), a name and a colour ("#rrggbb")."""
    start: int
    end: int
    name: str = ""
    color: str = REGION_COLOR


@dataclass
class Group:
    name: str
    keys: list = field(default_factory=list)
    collapsed: bool = False


class Layout:
    def __init__(self, on_change=None):
        self.nodes = [Group(SEQUENCE, list(BUILTIN))]
        self.styles = {}
        self.regions = []
        self.on_change = on_change or (lambda edited=True: None)

    # ---- reading ---------------------------------------------------------------------------------------------------------------
    def style(self, key):
        st = self.styles.get(key)
        if st is None:
            st = self.styles[key] = default_style(key)
        return st

    def label(self, key):
        st = self.styles.get(key)
        return st.label or key if st else key

    def keys(self):
        """Every row key in display order."""
        return [k for node in self.nodes for k in (node.keys if isinstance(node, Group) else [node])]

    def groups(self):
        return [n for n in self.nodes if isinstance(n, Group)]

    def group(self, name):
        return next((g for g in self.groups() if g.name == name), None)

    # ---- changing (each calls on_change) --------------------------------------------------------------------------------------
    def add(self, keys, group=None):
        """New rows at the end: of group `group` (made at the end of the list if new), else of the list."""
        have = set(self.keys())
        self._put([k for k in keys if k not in have], group)
        self.on_change(False)

    def remove(self, keys):
        """Rows taken out (their styles too); a group left empty goes."""
        self._take(keys)
        for k in keys:
            self.styles.pop(k, None)
        self.nodes = [n for n in self.nodes if not (isinstance(n, Group) and not n.keys)]
        self.on_change(False)

    def move(self, keys, group=None):
        """Rows `keys` (in their display order) to the end of group `group` (made if new), or of the top level (None)."""
        self._move(keys, group)
        self.on_change(True)

    def restyle(self, keys, **changes):
        """Set style fields of rows `keys`; a new display mode drops a height set by hand (each mode has its own)."""
        self.restyle_each({k: changes for k in keys})

    def restyle_each(self, changes):
        """{key: {field: value}}: each row its own values (a colour range, the heights Fit gives), in one change."""
        for k, fields in changes.items():
            st = self.style(k)
            if "mode" in fields and fields["mode"] != st.mode:
                st.height = None
            for name, value in fields.items():
                setattr(st, name, value)
        self.on_change(True)

    def collapse(self, name, collapsed):
        """Fold a group in the view to its header line, or unfold it."""
        if (g := self.group(name)) is not None:
            g.collapsed = collapsed
            self.on_change(True)

    def overlay(self, keys):
        """Draw signal rows `keys` in one row: the first one's (the others keep their place, to come back to it)."""
        self.restyle_each({k: {"overlay": keys[0]} for k in keys[1:]} | {keys[0]: {"overlay": None}})

    def separate(self, keys):
        """Rows `keys` and those drawn in them get rows of their own again."""
        hosts = set(keys)
        self.restyle_each({k: {"overlay": None} for k, st in self.styles.items() if k in hosts or st.overlay in hosts})

    def add_regions(self, regions):
        """Highlight `regions` (Region), kept in order of start."""
        self.regions = sorted(self.regions + list(regions), key=lambda r: (r.start, r.end))
        self.on_change(True)

    def change_regions(self, regions, **fields):
        """Rename or recolour highlights (`name`, `color`)."""
        for r in regions:
            for name, value in fields.items():
                setattr(r, name, value)
        self.on_change(True)

    def remove_regions(self, regions):
        gone = {id(r) for r in regions}
        self.regions = [r for r in self.regions if id(r) not in gone]
        self.on_change(True)

    def ungroup(self, name):
        """The group's rows go to the top level, where the group was."""
        g = self.group(name)
        if g is not None:
            i = self.nodes.index(g)
            self.nodes[i:i + 1] = g.keys
            self.on_change(True)

    def rename_group(self, old, new):
        g = self.group(old)
        if g is None or not new or (new != old and self.group(new) is not None):
            return False
        g.name = new
        self.on_change(True)
        return True

    def group_by_name(self, keys, parts):
        """Group by Name: rows `keys` go into groups named by the first `parts` parts of their names (see name_groups); a group
        emptied by it goes. Returns {label: group name} of what moved."""
        labels = {k: self.label(k) for k in keys}
        groups = name_groups(list(dict.fromkeys(labels.values())), parts)
        for name in dict.fromkeys(groups.values()):
            self._move([k for k in keys if groups.get(labels[k]) == name], name)
        self.nodes = [n for n in self.nodes if not (isinstance(n, Group) and not n.keys)]
        self.on_change(True)
        return groups

    def set_tree(self, tree):
        """The whole order at once (the panel after a drag and drop): [(group name or None, [keys])], a None entry being rows
        at the top level. Each group keeps its collapsed state; a key not given is kept at the end (never lost)."""
        old = {g.name: g.collapsed for g in self.groups()}
        nodes, seen = _nodes([(name, keys, old.get(name, False)) for name, keys in tree])
        left = [k for k in self.keys() if k not in seen]
        self.nodes = nodes
        self._put(left, None)
        self.on_change(True)

    def _move(self, keys, group):
        wanted = set(keys)
        order = [k for k in self.keys() if k in wanted]
        self._take(order)
        self._put(order, group)

    def _take(self, keys):
        gone = set(keys)
        self.nodes = [n for n in self.nodes if not (isinstance(n, str) and n in gone)]
        for g in self.groups():
            g.keys = [k for k in g.keys if k not in gone]

    def _put(self, keys, group):
        if group is None:
            self.nodes.extend(keys)
        else:
            g = self.group(group)
            if g is None:
                g = Group(group)
                self.nodes.append(g)
            g.keys.extend(keys)

    # ---- project --------------------------------------------------------------------------------------------------------------
    def to_json(self):
        nodes = [{"group": n.name, "keys": list(n.keys), "collapsed": n.collapsed} if isinstance(n, Group) else n for n in self.nodes]
        styles = {}
        for k, s in self.styles.items():
            base = asdict(default_style(k))
            styles[k] = {f: v for f, v in asdict(s).items() if v != base[f]}
        return {"nodes": nodes, "styles": {k: s for k, s in styles.items() if s}, **({"regions": [asdict(r) for r in self.regions]} if self.regions else {})}

    def load_json(self, data, present, length=None):
        """A saved layout, for the rows this file has now (`present`): rows it no longer has are dropped, rows it did not know
        are added at the end, and a style value that cannot be read is left at its default (an older or hand-edited project); a
        highlight that cannot be read, or that ends past `length` (the sequence's), is left out."""
        if not data:
            return
        present = list(present)
        known = set(present)
        tree = [(str(n["group"]), [k for k in n.get("keys", []) if isinstance(k, str)], bool(n.get("collapsed")))
                if isinstance(n, dict) else (None, [n], False)
                for n in data.get("nodes", []) if isinstance(n, str) or (isinstance(n, dict) and n.get("group"))]
        self.nodes, seen = _nodes(tree, known)
        self._put([k for k in present if k not in seen], None)
        self.styles = {}
        for k, values in (data.get("styles") or {}).items():
            if k in known and isinstance(values, dict):
                for f, v in values.items():
                    if f in DEFAULT and _valid(f, v):
                        setattr(self.style(k), f, v)
        self.regions = [Region(r["start"], r["end"], str(r.get("name", "")), r.get("color") or REGION_COLOR)
                        for r in data.get("regions") or [] if isinstance(r, dict) and _region(r) and (length is None or r["end"] <= length)]
        self.on_change(False)


def _nodes(tree, known=None):
    """(nodes, keys placed) from [(group name or None, keys, collapsed)]: each key once (and only `known` ones, when given)."""
    nodes, seen = [], set()
    for name, keys, collapsed in tree:
        keys = [k for k in keys if k not in seen and (known is None or k in known)]
        seen.update(keys)
        if name is None:
            nodes.extend(keys)
        else:
            nodes.append(Group(name, keys, collapsed))
    return nodes, seen


def _valid(name, value):
    return {"shown": lambda v: isinstance(v, bool), "mode": lambda v: v in MODES,
            "height": lambda v: v is None or (type(v) is int and v > 0), "scale": lambda v: v in SCALES,
            "color": lambda v: v is None or bool(re.fullmatch(r"#[0-9a-fA-F]{6}", str(v))),
            "label": lambda v: v is None or isinstance(v, str), "render": lambda v: v in RENDERS,
            "summary": lambda v: v in SUMMARIES, "range": lambda v: v in RANGES,
            "low": lambda v: type(v) in (int, float), "high": lambda v: type(v) in (int, float),
            "overlay": lambda v: v is None or isinstance(v, str), "frozen": lambda v: isinstance(v, bool)}[name](value)


def _region(r):
    """Can this saved highlight be read (whole positions, start before end, a colour or none)?"""
    return type(r.get("start")) is int and type(r.get("end")) is int and 0 <= r["start"] < r["end"] and _valid("color", r.get("color"))


def regions_bed(regions, seqid):
    """Highlights as BED lines (sequence, start, end, name, score 0, strand ., thickStart, thickEnd, itemRgb), so a genome
    browser shows them in their colours."""
    def rgb(color):
        return ",".join(str(int(color[i:i + 2], 16)) for i in (1, 3, 5))
    return "".join(f"{seqid}\t{r.start}\t{r.end}\t{r.name or '.'}\t0\t.\t{r.start}\t{r.end}\t{rgb(r.color)}\n" for r in regions)


def read_regions_bed(text, length):
    """([Region], [what was not read]) from BED text: name and itemRgb kept when there; a line that is not a region of a
    sequence of `length` bases is named, not dropped silently. The sequence column is not checked (a BED of another name
    for the same sequence is common)."""
    regions, unread = [], []
    for n, line in enumerate(text.splitlines(), 1):
        f = line.split("\t") if "\t" in line else line.split()
        if not line.strip() or f[0].startswith(("#", "track", "browser")):
            continue
        try:
            start, end = int(f[1]), int(f[2])
        except (IndexError, ValueError):
            unread.append(f"line {n}: not sequence, start, end")
            continue
        if not 0 <= start < end <= length:
            unread.append(f"line {n}: {start:,}-{end:,} is not within 0-{length:,}")
            continue
        color = REGION_COLOR
        if len(f) > 8 and re.fullmatch(r"\d+,\d+,\d+", f[8]) and all(int(c) < 256 for c in f[8].split(",")):   # (as tracks read it)
            color = "#%02x%02x%02x" % tuple(int(c) for c in f[8].split(","))
        regions.append(Region(start, end, f[3] if len(f) > 3 and f[3] != "." else "", color))
    return regions, unread


def name_groups(labels, parts=1):
    """{label: group name} for "Group by name": the first `parts` parts of each name (split on _ - . and spaces). A label with
    no more parts than that, or whose group would hold it alone, gets no group (it stays where it is)."""
    prefix = {}
    for label in labels:
        bits = [b for b in SPLIT.split(label) if b]
        if len(bits) > parts:
            prefix[label] = "_".join(bits[:parts])
    counts = Counter(prefix.values())
    return {label: g for label, g in prefix.items() if counts[g] > 1}


def gradient(color, n):
    """n colours from a light tint of `color` ("#rrggbb") to a dark shade of it, in order: P1 to P6 shade light to dark."""
    rgb = [int(color[i:i + 2], 16) for i in (1, 3, 5)]
    light = [c + (255 - c) * 0.65 for c in rgb]
    dark = [c * 0.6 for c in rgb]
    if n == 1:
        return [color.lower()]
    return ["#%02x%02x%02x" % tuple(round(a + (b - a) * i / (n - 1)) for a, b in zip(light, dark)) for i in range(n)]


def fit_heights(bounds, available, gap):
    """{key: px} sharing `available` px between rows (each with its own gap under it): the Fit button. `bounds`: {key: (least,
    natural)}; a row is never made taller than its natural height (room to spare stays spare) nor under its least (when even
    that does not fit, the view scrolls)."""
    if not bounds:
        return {}
    share = int((available - gap * len(bounds)) // len(bounds))
    return {k: max(least, min(natural, share)) for k, (least, natural) in bounds.items()}
