"""How BAM reads and variants are shown: the "Alignments" settings, IGV-style read colours
and the variant / consensus thresholds. Used by the genome view and the Variants pane.

Read colours (mode "Pair orientation & insert size", like IGV's default):
- grey: a normal pair (left read forward, right read reverse, insert within limits), or
  an unpaired read;
- red / blue: insert larger / smaller than expected (deletion / insertion in the sample,
  or a chimeric fragment);
- teal (both forward), indigo (both reverse), green (left reverse, right forward): pairs
  in the wrong orientation (inversion, duplication, circular junction…);
- orange: mate on another sequence; purple outline: mate unmapped.
Insert limits come from the reads themselves (0.5th / 99.5th percentile of proper pairs,
IGV's default) when enough pairs are in view, else from the settings.
"""
from types import SimpleNamespace

from PySide6.QtCore import QObject, QThreadPool, Signal

import tool_settings
from tool_settings import Param

TOOL = "Alignments"
COLOR_MODES = ["Pair orientation & insert size", "Strand", "None"]
tool_settings.register(TOOL, [
    Param("color_by", "Colour reads by", COLOR_MODES[0], choices=COLOR_MODES),
    Param("insert_auto", "Expected insert size from the reads", True,
          tip="0.5th–99.5th percentile of properly paired reads in view (IGV's default)."),
    Param("min_insert", "…otherwise smaller than expected below", 50, 0, 100_000, step=10, suffix=" bp"),
    Param("max_insert", "…and larger than expected above", 1000, 0, 1_000_000, step=50, suffix=" bp"),
    Param("reads_span", "Show reads when zoomed under", 30_000, 1000, 1_000_000, step=1000, suffix=" bp",
          tip="IGV's visibility window: wider views show coverage only."),
    Param("variant_pct", "Colour coverage where reads differ from the reference above", 20, 1, 100, suffix=" %",
          tip="IGV's allele-frequency threshold: the coverage bar is split by base there."),
    Param("consensus_pct", "Flag as consensus change (SNP / indel) above", 80, 1, 100, suffix=" %",
          tip="Most reads disagree with the sequence: flagged in the genome view and ticked in Variants."),
    Param("min_depth", "Minimum depth to flag", 5, 1, 10_000, suffix=" reads"),
    Param("min_baseq", "Minimum base quality counted", 13, 0, 60,
          tip="Phred score; lower-quality bases are left out of variant counts."),
    Param("min_mapq", "Minimum mapping quality counted", 0, 0, 60,
          tip="Reads below it are drawn but left out of variant counts."),
])

COLORS = {
    "normal": "#bdbdbd", "larger": "#e53935", "smaller": "#1e88e5", "LL": "#26a69a", "RR": "#3949ab",
    "RL": "#43a047", "elsewhere": "#fb8c00", "forward": "#d7a9a5", "reverse": "#a5bcd7",
}
LABELS = {
    "normal": "", "larger": "insert larger than expected", "smaller": "insert smaller than expected",
    "LL": "pair both forward", "RR": "pair both reverse", "RL": "pair outward (RL)",
    "elsewhere": "mate on another sequence", "forward": "", "reverse": "",
}
MATE_UNMAPPED_COLOR = "#8e24aa"
# Local Align: a piece of the query that matches the REVERSE complement (an inverted piece) is drawn in this colour
# (forward pieces use COLORS["forward"]), with a hatch and a "↺" in its label.
INVERTED_COLOR = "#7e57c2"
PAIRED, PROPER, MATE_UNMAPPED, MATE_REVERSE = 0x1, 0x2, 0x8, 0x20


def setting(key):
    return tool_settings.value(TOOL, key)


def insert_limits(reads):
    """(smaller-than, larger-than) bp for the reads in view: percentiles of proper pairs
    when there are at least 100 of them and the setting asks for it, else the settings."""
    if setting("insert_auto"):
        sizes = sorted(abs(r[11]) for r in reads if r[8] & PROPER and r[11])
        if len(sizes) >= 100:
            return sizes[int(len(sizes) * 0.005)], sizes[min(len(sizes) - 1, int(len(sizes) * 0.995))]
    return setting("min_insert"), setting("max_insert")


def classify(read, limits, mode=None):
    """Colour category of a read tuple (see Bam.reads) under the current colour mode."""
    mode = mode or setting("color_by")
    start, strand, flags, mate_start, elsewhere, insert = read[0], read[2], read[8], read[9], read[10], read[11]
    if mode == "Strand":
        return "forward" if strand > 0 else "reverse"
    if mode == "None" or not flags & PAIRED or flags & MATE_UNMAPPED or mate_start is None:
        return "normal"
    if elsewhere:
        return "elsewhere"
    forward, mate_forward = strand > 0, not flags & MATE_REVERSE
    left_fwd, right_fwd = (forward, mate_forward) if start <= mate_start else (mate_forward, forward)
    orientation = {(True, False): "LR", (True, True): "LL", (False, False): "RR", (False, True): "RL"}[(left_fwd, right_fwd)]
    if orientation != "LR":
        return orientation
    lo, hi = limits
    size = abs(insert)
    if size and size > hi:
        return "larger"
    if size and size < lo:
        return "smaller"
    return "normal"


def flag_level(count, depth, limits=None):
    """"consensus", "variant" or None for `count` of `depth` reads differing. `limits`:
    (min depth, variant %, consensus %), read from the settings if not given."""
    min_depth, variant_pct, consensus_pct = limits or (setting("min_depth"), setting("variant_pct"),
                                                       setting("consensus_pct"))
    if depth < min_depth or not depth:
        return None
    frac = 100 * count / depth
    if frac >= consensus_pct:
        return "consensus"
    if frac >= variant_pct:
        return "variant"
    return None


class _Relay(QObject):
    """Tells views (GUI thread) that a background variant scan finished. Lives forever."""
    variants_ready = Signal()


relay = _Relay()


def describe(pos, kind, ref, alt, count, depth):
    pct = f"{100 * count / max(1, depth):.0f} % of {depth}"
    if kind == "snp":
        return f"{pos + 1:,} {ref}>{alt}  {pct}"
    if kind == "del":
        return f"{pos + 1:,} deletion of {len(ref)} bp ({ref if len(ref) <= 12 else ref[:12] + '…'})  {pct}"
    return f"{pos + 1:,} insertion of {alt or '?'}  {pct}"


def start_scan(entry, editor):
    """Scan `entry`'s reads for positions differing from the current sequence, in a
    background thread; sets entry.variants = [SimpleNamespace(pos, kind, ref, alt, count,
    depth, level, text)] and emits relay.variants_ready. Settings changes start a new scan
    (the token makes a late result of an older scan be dropped)."""
    if getattr(entry, "scanning", False) and entry.variants is None and entry.scan_settings == _scan_settings():
        return
    entry.scan_token = getattr(entry, "scan_token", 0) + 1
    entry.scanning, entry.scan_settings = True, _scan_settings()
    token, reference = entry.scan_token, editor.doc.sequence()
    min_depth, variant_pct, consensus_pct, min_baseq, min_mapq = entry.scan_settings
    bam, seqid = entry.bam, entry.seqid

    def work():  # pool thread; Rust releases Python's lock
        try:
            found = bam.variants(seqid, reference, min_depth, variant_pct / 100, min_baseq, min_mapq)
        except OSError as e:
            found, entry.error = [], str(e)
        result = [SimpleNamespace(pos=pos, kind=kind, ref=ref, alt=alt, count=n, depth=d,
                                  level=flag_level(n, d, (min_depth, variant_pct, consensus_pct)) or "variant",
                                  text=describe(pos, kind, ref, alt, n, d))
                  for pos, kind, ref, alt, n, d in found]
        if token == entry.scan_token:
            entry.variants, entry.scanning = result, False
            relay.variants_ready.emit()

    QThreadPool.globalInstance().start(work)


def _scan_settings():
    return (setting("min_depth"), setting("variant_pct"), setting("consensus_pct"), setting("min_baseq"),
            setting("min_mapq"))
