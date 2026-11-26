"""UCSC track hubs: read a hub (its hub.txt, genomes file and trackDb), list its tracks in the hub's own groups, and copy their
files into a folder. No Qt: the window (main.import_hub) and the command line (`hub`) call the same functions.

A hub is a hub.txt (a URL or a local file, or a folder holding one): `genomesFile` names the genomes file, whose stanzas each give
a `genome` and the `trackDb` of its tracks. A trackDb is stanzas of "key value" lines, blank lines between them (indentation
means nothing; `include other.txt` reads another file in place). A track with a `bigDataUrl` has data; one without is a container
(a composite, a super track, a view) that others name as their `parent`. Settings a track does not give are its parent's (type,
color, itemRgb), as in the UCSC browser.

openDNA reads the bigWig and bigBed tracks; any other type (bam, vcfTabix, hic...) is listed with the reason it is not read,
never dropped silently. The tracks are always copied (never read from the server: an edit must land on the user's own copy).
"""
import posixpath
import queue
import re
import threading
from collections import Counter
import urllib.parse
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path

import twobit
from user_error import UserError

READ = ("bigWig", "bigBed")   # the track types openDNA reads
UCSC_2BIT = "https://hgdownload.soe.ucsc.edu/goldenPath/{0}/bigZips/{0}.2bit"   # a genome the hub names without a twoBitPath (hg38, mm10...)


def is_url(source):
    return bool(re.match(r"(?i)https?://", str(source)))


def open_url(url, timeout=60):
    """An open reply for `url` (the one network call here: tests replace it), with the program's name and version as the user
    agent (some servers refuse a bare name: see structure.USER_AGENT)."""
    import structure
    return urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": structure.USER_AGENT}), timeout=timeout)


def fetch(url, timeout=60):
    """The bytes at `url` (a hub's text files: small)."""
    with open_url(url, timeout) as reply:
        return reply.read()


def read_text(source):
    try:
        data = fetch(source) if is_url(source) else Path(source).read_bytes()
    except OSError as e:   # (urllib's errors are OSErrors)
        raise UserError(f"cannot read {source}: {e}") from None
    return data.decode("utf-8", errors="replace")


def resolve(base, rel):
    """`rel` (a URL, an absolute path or a path relative to the file `base`) as a URL or a path."""
    if is_url(rel) or Path(rel).is_absolute():
        return rel
    if is_url(base):
        return urllib.parse.urljoin(base, rel)
    return str(Path(base).parent / rel)


def stanzas(text):
    """[{key: value}] of a hub file: blank lines between stanzas, `#` comments, indentation ignored."""
    out, current = [], {}
    for line in text.splitlines() + [""]:
        line = line.strip()
        if not line:
            if current:
                out.append(current)
            current = {}
        elif not line.startswith("#"):
            key, *value = line.split(None, 1)   # (a space or a tab between them)
            current[key] = value[0].strip() if value else ""
    return out


@dataclass
class Hub:
    name: str
    label: str
    source: str           # where its hub.txt is
    genomes: dict         # genome -> where its trackDb is
    twobit: dict = field(default_factory=dict)   # genome -> where its .2bit sequence is (only those the hub gives)


@dataclass
class HubTrack:
    name: str             # the trackDb's `track`
    label: str            # its shortLabel (what the view shows)
    url: str              # its bigDataUrl, made whole (a URL or a path)
    kind: str             # the first word of its type (bigWig, bigBed, bam...)
    group: str            # the label of the container it is in (a composite, a super track), else the hub's
    hidden: bool          # off when the hub opens (visibility hide, an `off` subtrack, a super track not shown)
    color: str = None     # "#rrggbb" from `color`, None: none given
    item_rgb: bool = False   # bigBed items carry their own colours
    why: str = ""         # why openDNA does not read it ("" when it does)

    @property
    def file_name(self):
        return posixpath.basename(urllib.parse.urlparse(self.url).path) if is_url(self.url) else Path(self.url).name


def open_hub(source):
    """The hub at `source`: a hub.txt as a URL or a path, or a folder holding hub.txt."""
    source = str(source)
    if not is_url(source) and Path(source).is_dir():
        source = str(Path(source) / "hub.txt")
    first = (stanzas(read_text(source)) or [{}])[0]
    if "hub" not in first or "genomesFile" not in first:
        raise UserError(f"{source} is not a track hub (its first lines must name `hub` and `genomesFile`)")
    genomes_file = resolve(source, first["genomesFile"])
    genomes = {s["genome"]: resolve(genomes_file, s["trackDb"]) for s in stanzas(read_text(genomes_file)) if "genome" in s and "trackDb" in s}
    if not genomes:
        raise UserError(f"{genomes_file} names no genome with a trackDb")
    twobit = {s["genome"]: resolve(genomes_file, s["twoBitPath"]) for s in stanzas(read_text(genomes_file)) if "genome" in s and "twoBitPath" in s}
    return Hub(first["hub"], first.get("shortLabel") or first["hub"], source, genomes, twobit)


def _trackdb(location, seen=None):
    """The stanzas of a trackDb with its `include` files read in place (each once)."""
    seen = set() if seen is None else seen
    seen.add(location)
    out = []
    for s in stanzas(read_text(location)):
        if set(s) == {"include"}:
            inner = resolve(location, s["include"].split()[0])
            if inner not in seen:
                out += _trackdb(inner, seen)
        else:
            out.append(dict(s, _base=location))
    return out


def hub_tracks(hub, genome):
    """Every track of `genome` that has data, in the order of the trackDb, with its group, label, colour and whether it starts
    hidden; the types openDNA does not read are listed with `why`."""
    if genome not in hub.genomes:
        raise UserError(f"the hub has no genome {genome!r} (it has: {', '.join(hub.genomes)})")
    all_ = [s for s in _trackdb(hub.genomes[genome]) if "track" in s]
    by_name = {s["track"].split()[0]: s for s in all_}

    def parents(s):
        seen = set()
        while (p := s.get("parent")) and (name := p.split()[0]) in by_name and name not in seen:
            seen.add(name)
            s = by_name[name]
            yield s

    def setting(s, key):
        return next((x[key] for x in (s, *parents(s)) if key in x), None)

    def off(s):
        """Hidden when the hub opens: visibility hide (here or above), a subtrack switched off, a super track not shown."""
        for x in (s, *parents(s)):
            if x.get("visibility") == "hide" or x.get("parent", "").split()[1:2] == ["off"]:
                return True
            if x is not s and x.get("superTrack", "").split()[:1] == ["on"] and "show" not in x["superTrack"].split():
                return True
        return False

    out = []
    for s in all_:
        if "bigDataUrl" not in s:
            continue
        kind = (setting(s, "type") or "").split()[:1] or [""]
        container = next((p for p in parents(s) if "view" not in p), None)   # (a composite's views are not groups of their own)
        color = setting(s, "color")
        rgb = re.fullmatch(r"\s*(\d{1,3}),(\d{1,3}),(\d{1,3})\s*", color or "")
        out.append(HubTrack(
            name=s["track"].split()[0], label=s.get("shortLabel") or s["track"].split()[0],
            url=resolve(s["_base"], s["bigDataUrl"]), kind=kind[0],
            group=(container.get("shortLabel") or container["track"].split()[0]) if container else hub.label,
            hidden=off(s), color="#%02x%02x%02x" % tuple(min(255, int(v)) for v in rgb.groups()) if rgb else None,
            item_rgb=(setting(s, "itemRgb") or "").lower() == "on",
            why="" if kind[0] in READ else f"type {kind[0] or 'not given'}: openDNA reads bigWig and bigBed tracks"))
    return out


def local_names(tracks):
    """{track name: file name} for copying `tracks` into one folder: the file's own name, the track's name before it when two
    files share a name."""
    counts = Counter(t.file_name for t in tracks)
    return {t.name: t.file_name if counts[t.file_name] == 1 else f"{t.name}_{t.file_name}" for t in tracks}


STALL = 15   # seconds a server may stay silent (connecting, or between two chunks) before the file is given up
WORKERS = 4   # files copied at once
CHUNK = 1 << 16


def copy_file(location, path, progress=None):
    """Copy (or download) the file at `location` (a URL or a path) to `path`, through a temporary file (a stopped or failed copy
    leaves no half file). `progress(done, total)` (total None when a server does not say) is told after each chunk and may return
    False to stop: copy_file then returns False. A server silent for STALL seconds raises OSError."""
    part = Path(str(path) + ".part")
    stopped, remote = False, is_url(location)
    try:
        with (open_url(location, timeout=STALL) if remote else open(location, "rb")) as reply, open(part, "wb") as out:
            length = getattr(reply, "headers", {}).get("Content-Length") if remote else str(Path(location).stat().st_size)
            total = int(length) if length and length.isdigit() else None
            done = 0
            while chunk := reply.read(CHUNK):   # (streamed: a bigWig can be hundreds of MB)
                out.write(chunk)
                done += len(chunk)
                if progress is not None and progress(done, total) is False:
                    stopped = True
                    break
        if not stopped:
            part.replace(path)
    except OSError as e:
        if isinstance(e, TimeoutError) or isinstance(getattr(e, "reason", None), TimeoutError):   # (a connect timeout comes as a URLError)
            raise OSError(f"the server did not answer for {STALL} seconds") from None
        raise
    finally:
        part.unlink(missing_ok=True)
    return not stopped


class Transfer:
    """One file being copied by a worker thread, for the window to show: `state` is waiting, copying, done, failed (`error` says
    why) or cancelled; `done` and `total` count bytes (total 0 when unknown), `text` adds a word ("converting chr1"). `cancel()`
    takes effect at once for the window; the worker gives up at its next chunk (or when the server's silence times out), unseen."""

    def __init__(self, label, job, track=None, path=None):
        self.label, self.job, self.track, self.path = label, job, track, path
        self.state, self.done, self.total, self.text, self.error, self.result = "waiting", 0, 0, "", "", None
        self.stop = threading.Event()

    def update(self, done, total, text=""):
        """The job's progress call: False once cancelled."""
        self.done, self.total, self.text = done, total or 0, text
        return not self.stop.is_set()

    def cancel(self):
        self.stop.set()
        if self.state in ("waiting", "copying"):
            self.state = "cancelled"

    def again(self):
        """A new Transfer of the same file for Retry (the old one may still have a thread stuck on a silent server: it is left alone)."""
        return Transfer(self.label, self.job, self.track, self.path)

    def finished(self):
        return self.state in ("done", "failed", "cancelled")

    def run(self):
        if self.stop.is_set():
            return
        self.state = "copying"
        if self.stop.is_set():   # (cancelled between the two lines above)
            self.state = "cancelled"
            return
        try:
            self.result = self.job(self)
        except (OSError, UserError) as e:
            self.error = str(e)
        self.state = "cancelled" if self.stop.is_set() else "failed" if self.error else "done"
        if self.state == "cancelled":
            self.error = ""   # (a cancelled file's late error is unseen)


def start(transfers, workers=WORKERS):
    """Run `transfers` on `workers` threads, the rest waiting their turn. Returns the threads (daemon: one stuck on a silent
    server never keeps the program from closing)."""
    todo = queue.Queue()
    for t in transfers:
        todo.put(t)

    def work():
        while True:
            try:
                todo.get_nowait().run()
            except queue.Empty:
                return
    threads = [threading.Thread(target=work, daemon=True) for _ in range(min(workers, len(transfers)))]
    for t in threads:
        t.start()
    return threads


def track_transfers(tracks, folder):
    """A Transfer for each of `tracks` into `folder` (made if needed). A file already there is kept (copying the hub again
    reuses it)."""
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    names = local_names(tracks)

    def one(path, location):
        return lambda t: path.exists() or copy_file(location, path, t.update)
    return [Transfer(t.label, one(folder / names[t.name], t.url), t, folder / names[t.name]) for t in tracks]


def copied(transfers):
    """[(track, path or None, error)], one per track transfer: one that failed or was stopped is said, never dropped."""
    return [(t.track, t.path if t.state == "done" and t.result is not False else None, t.error or "stopped") for t in transfers]


def copy_tracks(tracks, folder):
    """Copy the files of `tracks` into `folder`, several at once (the command line; the window shows the Transfers).
    Returns what `copied` says."""
    transfers = track_transfers(tracks, folder)
    for thread in start(transfers):
        thread.join()
    return copied(transfers)


def pick_genome(hub, names):
    """The hub's genome for a sequence known by `names` (its name, its record...): the one of the same name (any case), else the
    only one, else None (the caller asks)."""
    wanted = {str(n).lower() for n in names if n}
    return next((g for g in hub.genomes if g.lower() in wanted), None) or (next(iter(hub.genomes)) if len(hub.genomes) == 1 else None)


def genome_source(hub, genome):
    """Where the sequence of `genome` is: the hub's own `twoBitPath`, else UCSC's file of that name (a genome the hub only maps
    its tracks on, such as hg38)."""
    return hub.twobit.get(genome) or UCSC_2BIT.format(genome)


def genome_size(location):
    """The size in bytes of the .2bit file at `location` (a URL or a path), None when a server does not say."""
    if not is_url(location):
        return Path(location).stat().st_size
    with open_url(location) as reply:
        length = reply.headers.get("Content-Length")
    return int(length) if length and length.isdigit() else None


def copy_genome(hub, genome, folder, progress=None):
    """The sequence of `genome` as the FASTA file `folder`/<genome>.fa (the .2bit file downloaded or read, then written as FASTA;
    the downloaded copy is not kept). Returns (path, [(sequence, length)]), None when `progress` stopped it. `progress(done,
    total, name)` is told after each chunk of the download (name "": bytes) and after each sequence (name: the sequence). A FASTA file already there is kept. UserError when the sequence cannot be read."""
    source = genome_source(hub, genome)
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    fasta = folder / f"{Path(genome).name}.fa"
    local = folder / f"{Path(genome).name}.2bit"
    downloaded = is_url(source)
    try:
        if downloaded:
            if not local.exists() and not copy_file(source, local, None if progress is None else lambda d, n: progress(d, n, "")):
                return None
        else:
            local = Path(source)
        lengths = twobit.sequences(local)
        if not fasta.exists() and twobit.to_fasta(local, fasta, progress) is None:
            return None
    except OSError as e:   # (urllib's errors are OSErrors)
        raise UserError(f"cannot get the sequence of {genome} from {source}: {e}") from None
    finally:
        if downloaded:
            local.unlink(missing_ok=True)
    return fasta, lengths
