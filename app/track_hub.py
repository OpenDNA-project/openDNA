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
import re
import shutil
from collections import Counter
import urllib.parse
import urllib.request
from dataclasses import dataclass
from pathlib import Path

from user_error import UserError

READ = ("bigWig", "bigBed")   # the track types openDNA reads


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
    return Hub(first["hub"], first.get("shortLabel") or first["hub"], source, genomes)


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


def copy_track(track, path):
    """Copy (or download) one track's file to `path`, through a temporary file (a stopped download leaves no half file)."""
    part = Path(str(path) + ".part")
    if is_url(track.url):
        with open_url(track.url, timeout=300) as reply, open(part, "wb") as out:   # (streamed: a bigWig can be hundreds of MB)
            shutil.copyfileobj(reply, out)
    else:
        shutil.copyfile(track.url, part)
    part.replace(path)


def copy_tracks(tracks, folder, progress=None):
    """Copy the files of `tracks` into `folder` (made if needed). A file already there is kept (copying the hub again reuses
    it). `progress(done, total, track)` is told after each file and may return False to stop. Returns [(track, path or None,
    error or "")], one per track: a track that failed is said, never dropped."""
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    names = local_names(tracks)
    out = []
    for i, t in enumerate(tracks):
        path = folder / names[t.name]
        try:
            if not path.exists():
                copy_track(t, path)
            out.append((t, path, ""))
        except OSError as e:
            out.append((t, None, str(e)))
        if progress is not None and progress(i + 1, len(tracks), t) is False:
            out += [(rest, None, "stopped") for rest in tracks[i + 1:]]
            break
    return out


def pick_genome(hub, names):
    """The hub's genome for a sequence known by `names` (its name, its record...): the one of the same name (any case), else the
    only one, else None (the caller asks)."""
    wanted = {str(n).lower() for n in names if n}
    return next((g for g in hub.genomes if g.lower() in wanted), None) or (next(iter(hub.genomes)) if len(hub.genomes) == 1 else None)
