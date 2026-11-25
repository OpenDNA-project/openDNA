"""The commands of the command line for agents that are not table batches (see cli.py, which registers them and holds what they share).
Each one calls the function the matching window calls and answers with one JSON: the settings that shaped it (`settings`) and `warnings` to relay.
Read-only analyses (no `check`/`run` pair, nothing written): `sites`, `guides`, `repeats`, `search`, `align`, `msa` (`--out` writes one FASTA),
`checkup`, `annotate` (`--add` writes a copy), `variants` (`--apply` writes a copy). Files: `export`, `order`. Network: `ncbi`.
Two are batches with `check` then `run`, a table in and files out, like the ones in cli.py: `codon` and `gibson`."""
import json
import re
import sys
from contextlib import contextmanager
from pathlib import Path
from types import SimpleNamespace

C = None   # the cli module (set by register): its helpers (load, given, ...) are shared, not copied
UserError = None


def csv_list(text):
    return [t.strip() for t in (text or "").split(",") if t.strip()]


def answer(settings, caps=(), **fields):
    """The JSON of an analysis. `caps`: a text for every limit the run reached (a list cut short, a search stopped, a sequence not read in
    full): each is a warning, so a cut answer is never taken for the whole one."""
    return {"settings": settings, "warnings": C.setting_warnings(settings) + [c for c in caps if c], **fields}


def pick_enzymes(names):
    """The enzymes called `names` (each once, in order); an unknown name is an error that lists the known ones."""
    import enzymes
    by_name = {n.lower(): e for e in enzymes.COMMON for n in (e.name, *e.aliases)}   # an isoschizomer's name picks the enzyme listed
    picked = {}
    for name in names:
        if name.lower() not in by_name:
            raise UserError(f"No enzyme called {name!r}. The enzymes openDNA knows: {', '.join(e.name for e in enzymes.COMMON)}.")
        picked[by_name[name.lower()].name] = by_name[name.lower()]
    return list(picked.values())


def end1(end, n):
    """The 1-based last base of a 0-based end (exclusive) on a circle of n, so a piece that runs over the origin ends where it should."""
    return (end - 1) % n + 1


def region_bounds(doc, text):
    """(start, end) 0-based, end exclusive, of `START..END` (1-based, inclusive) or of the whole sequence when `text` is empty."""
    if not text:
        return 0, doc.length
    first, last = C.target_range(doc, text)
    return first - 1, last


# ---- sites: restriction enzymes and digests ---------------------------------------------------------------------------------------

def cmd_sites(args):
    import enzymes
    import tools
    doc = C.load(args.file)
    chosen = pick_enzymes(csv_list(args.enzymes))
    cutting = pick_enzymes(csv_list(args.digest))
    digest = [e.name for e in cutting]
    lo = 1 if args.min_sites is None else args.min_sites
    flt = enzymes.EnzymeFilter(min_sites=lo, max_sites=1 if args.min_sites is None and args.max_sites is None else args.max_sites)   # (the window's default: unique cutters)
    seq = doc.sequence(0, None).upper()
    cuts = enzymes.find_cuts(seq, doc.circular, list({e.name: e for e in chosen + cutting}.values()) if chosen else None)   # (None: all)
    listed = []
    for e in chosen or enzymes.COMMON:
        found = cuts.get(e.name, [])
        if chosen or enzymes.wanted(flt, e, len(found)):
            listed.append({"enzyme": e.name, "site": e.site, "notation": enzymes.notation(e), "ends": enzymes.describe(e).split(", ", 1)[1], "count": len(found),
                           "cuts": [{"top_after_base": t, "bottom_after_base": b} for t, b in found],   # (0 = before the first base)
                           "sites": [{"start": s + 1, "end": end1(en, doc.length), "strand": st} for s, en, st in enzymes.site_spans(seq, doc.circular, e)]})
    out = answer({"filter": {"min_sites": flt.min_sites, "max_sites": flt.max_sites, "min_site_length": flt.min_size,
                             "applies": "no: the enzymes were named" if chosen else "yes"}}, length=doc.length, circular=doc.circular, enzymes=listed)
    if digest:
        positions = {}
        for name in digest:
            for t, _ in cuts.get(name, []):
                positions.setdefault(t, []).append(name)
        out["digest"] = {"enzymes": digest, "fragments": [
            {"size": size, "start": a + 1, "end": b, "left": left, "right": right}
            for size, a, b, left, right in sorted(tools.fragments(doc.length, doc.circular, positions), reverse=True)]}
    return out


# ---- guides: Cas9 / Cas12a ---------------------------------------------------------------------------------------------------------

def cmd_guides(args):
    import cas9
    from cut_sites import NUCLEASES
    doc = C.load(args.file)
    by_name = {n.name.lower(): n for n in NUCLEASES.values()}
    if args.nuclease.lower() not in by_name:
        raise UserError(f"No nuclease called {args.nuclease!r}. Choices: {', '.join(n.name for n in by_name.values())}.")
    nuc = by_name[args.nuclease.lower()]
    a, b = region_bounds(doc, args.region)
    offset, seq, wrap, error = cas9.stretch(doc, (a, b) if args.region else None)
    if error:
        raise UserError(error)
    found = cas9.candidates(seq, doc.circular, nuc, offset, wrap)
    if args.top is not None and args.top < 1:
        raise UserError("--top must be at least 1.")
    shown = found[:args.top or cas9.MAX_CANDIDATES]
    max_mm = C.given(args.max_mm, cas9.DEFAULT_MAX_MM, cas9.MAX_MM_LIMITS, "--max-mm")
    ref, specs = None, []
    if shown:   # (the background genome is only opened when there is a guide to check)
        import primer_panel
        ref = None if args.no_genome else primer_panel.genome()
        specs = cas9.scan(doc.sequence(0, None), doc.circular, doc.name, shown, nuc, max_mm["value"], ref)
    n = doc.length
    guides = []
    for rank, (g, spec) in enumerate(zip(shown, specs), 1):
        guides.append({"rank": rank, **cas9.guide_row(nuc, g, spec, n)})
    genome_name = ref.name if ref is not None else None
    settings = {"cas9": {"nuclease": {"value": nuc.name, "source": "command line"}, "max_mismatches": max_mm,
                         "region": {"value": args.region or "whole sequence", "source": "command line" if args.region else "default"}},
                "background_genome": {"value": genome_name, "source": "none: --no-genome" if args.no_genome else "openDNA settings (Primers > Options)" if ref else "none set: no genome check"}}
    cap = f"Only the best {len(shown)} of the {len(found)} guides found are listed and checked (--top N lists more)." if len(found) > len(shown) else ""
    return answer(settings, caps=[cap, cas9.ambiguity_note(seq)], guides_found=len(found), guides=guides)


# ---- repeats: what two sequences share (the Compare view of Local Align) -----------------------------------------------------------

def cmd_repeats(args):
    import repeats
    a = C.load(args.file)
    b = C.load(args.other) if args.other else a
    lim = repeats.LIMITS
    settings = {"repeats": {"min_length": C.given(args.min_len, repeats.DEFAULT_MIN_LEN, lim["min_len"], "--min-len"),
                            "mismatches": C.given(args.mismatches, repeats.DEFAULT_MISMATCHES, lim["mismatches"], "--mismatches")}}
    sa = a.sequence(0, None)
    same = b is a
    sb = sa if same else b.sequence(0, None)
    circular = a.circular and b.circular   # (a stretch may run through the origin when both are circles)
    res = repeats.pair_repeats(sa, sb, settings["repeats"]["min_length"]["value"], settings["repeats"]["mismatches"]["value"], same=same, circular=circular)
    identical = not same and sa.upper() == sb.upper()
    mixed = "One sequence is circular and the other linear: stretches through an origin are not looked for." if a.circular != b.circular else ""
    note = "the two sequences are identical: they share their whole length (listed first)" if identical else ""
    return answer(settings, caps=[res.capped, mixed], compared={"a": Path(args.file).name, "b": Path(args.other or args.file).name, "same_file": same, "identical_sequences": identical},
                  repeats=[{"kind": k, "start_in_a": i + 1, "end_in_a": i + n, "start_in_b": j + 1, "end_in_b": j + n, "length": n, "differences": d}
                           for k, i, j, n, d in res.repeats], note=note)


# ---- synthesis: will a synthesis company make this fragment quickly? (the Synthesis tab's local row) --------------------------------

def cmd_synthesis(args):
    import synthesis_salis
    from codon_tools import circular_bases
    doc = C.load(args.file)
    a, b = region_bounds(doc, args.region)
    seq = circular_bases(doc, a, b)   # (a region through the origin is joined into one straight piece)
    result = synthesis_salis.evaluate(seq)
    settings = {"synthesis": {"model": {"value": "Synthesis Success Calculator (random forest, 9 features)", "source": "fixed"},
                              "region": {"value": args.region or "whole sequence", "source": "command line" if args.region else "default"}}}
    circle = "The file is a circular sequence, but the model reads it as a straight piece, as published." if doc.circular and not args.region else ""
    return answer(settings, caps=[circle], citation=synthesis_salis.CITATION, bases=len(seq), region=C.span_fields(doc, a + 1, b),   # (a region through the origin: end below start, crosses_origin)
                  predicted_quick_synthesis=result["outcome"], probability=result["probability"], features=result["features"],
                  # where each measure was counted, on the file (1-based; a stretch through the origin ends below its start)
                  places={name: [{"start": (a + x) % doc.length + 1, "end": end1(a + y, doc.length)} for x, y in spans] for name, spans in result["places"].items()})


# ---- search: the search bar's engine -----------------------------------------------------------------------------------------------

def cmd_search(args):
    import search
    doc = C.load(args.file)
    opts = search.Options(both=not args.single_strand, fuzzy=args.fuzzy or 0, translations=args.protein, details=args.details, enzymes=args.enzymes, primers=False)
    ed = doc_editor(doc)
    hits, note = search.search_sequence(ed, args.query, opts)
    n = doc.length
    return answer({"search": {"both_strands": {"value": opts.both}, "fuzzy_mismatches": {"value": opts.fuzzy}, "protein_motif_in_translations": {"value": opts.translations},
                              "feature_details": {"value": opts.details}, "enzyme_names": {"value": opts.enzymes}}},
                  caps=[note], query=args.query, hits=[{"kind": h.kind, "label": h.label, "start": h.start + 1, "end": end1(h.end, n), "strand": h.strand, "mismatches": h.mismatches,
                                                      **({"feature_index": h.feature} if h.feature is not None else {})} for h in hits])


# ---- align: where a query matches the plasmid (the Local Align window) ----------------------------------------------------------------

def query_bases(text):
    """The DNA of a QUERY argument: a file (GenBank or FASTA: its first record) when one is named, else the bases written there."""
    import local_align
    if Path(text).is_file():
        return C.load(text).sequence(0, None).upper()
    try:
        return local_align.clean_query(text)
    except local_align.QueryError as e:
        raise UserError(str(e))


def cmd_align(args):
    import local_align
    doc = C.load(args.file)
    query = query_bases(args.query)
    import tool_settings
    with tool_settings_block(local_align.TOOL, {"time_limit": "--time-limit"}, args) as limit,             tool_settings_block(tool_settings.PERFORMANCE, {"cores": "--threads"}, args) as cores:
        hits, notes = local_align.find_hits(doc.sequence(0, None), query, doc.circular)
    n = doc.length
    top = doc.sequence(0, None)
    rows = [{"piece": i + 1, "query_start": h.qstart + 1, "query_end": h.qend, "aligned_bp": h.qend - h.qstart, "start": h.start + 1, "end": end1(h.end, n), "strand": h.strand, "identity": round(h.identity, 4),
             "score": h.score, "mismatches": len(h.mismatches), "insertions": len(h.insertions), "deletions": len(h.deletions),
             # (bases as on the plasmid's top strand; a piece across the origin has positions past n: wrapped)
             "mismatch_list": [{"position": p % n + 1, "plasmid_base": top[p % n], "query_base": q} for p, q in h.mismatches],
             "insertion_list": [{"after_plasmid_base": end1(p, n), "extra_bases_in_query": k} for p, k in h.insertions],
             "deletion_list": [{"start": a % n + 1, "end": end1(b, n)} for a, b in h.deletions],   # (plasmid bases the query lacks)
             "before": local_align.junction_text(hits[i - 1], h, n, doc.circular) if i else "start of the query"} for i, h in enumerate(hits)]
    cap = f"Only {local_align.MAX_HITS} pieces are listed (the limit): weaker ones are left out." if len(hits) >= local_align.MAX_HITS else ""
    if local_align.time_limit_note(notes):
        notes.append(f"To search every candidate place, run again with --time-limit 0 (or a larger value); the window's setting is {local_align.WHERE_TIME_LIMIT}.")
    return answer({"align": {"query_bp": {"value": len(query)}, "max_pieces": {"value": local_align.MAX_HITS}, **limit, **cores}}, caps=[cap, *notes],
                  summary=local_align.summarize(hits, n, doc.circular) or "No hit", pieces=rows)


# ---- msa: several sequences aligned (the Multiple alignment tab) --------------------------------------------------------------------

def new_path(path):
    """`path` as a Path, refused when it exists: nothing is overwritten."""
    if Path(path).exists():
        raise UserError(f"{path} already exists: choose a new file name (nothing is overwritten).")
    return Path(path)


def records_of(files):
    """[(label, Document)] for the files given: a file of several records gives all of them (label FILE#NAME), `FILE#NAME` gives that one."""
    out = []
    for f in files:
        names = C.file_records(f) if Path(f).is_file() else None
        if names is None or len(names) == 1:
            out.append((f, C.load(f)))
        else:
            out += [(f"{f}#{n}", C.load(f"{f}#{n}")) for n in names]
    return out


def cmd_msa(args):
    import msa
    labels, docs = zip(*records_of(args.files)) if args.files else ((), ())
    names = msa.unique_names([d.name for d in docs])
    out = new_path(args.out) if args.out else None
    try:
        aln = msa.align(names, [d.sequence(0, None) for d in docs])
    except ValueError as e:
        raise UserError(str(e))
    if out:
        out.write_text(msa.to_fasta(aln.names, aln.rows), encoding="utf-8")
    return answer({"msa": {"aligner": {"value": f"built-in progressive aligner (up to {msa.MAX_BUILTIN} bp each)"}, "identity_against": {"value": names[0]}}},
                  width=aln.width, consensus=aln.consensus(),
                  sequences=[{"name": nm, "file": str(f), "length": d.length, "identity_to_first": round(aln.identity(i), 2), "aligned": row}
                             for i, (nm, f, d, row) in enumerate(zip(aln.names, labels, docs, aln.rows))], **({"file": str(out)} if out else {}))


# ---- split / merge: files of several records (a multi-FASTA, a GenBank file of several LOCUS records) ------------------------------------

def file_stem(name, taken):
    """A file name (no extension) for a record name: safe characters, and not one of `taken` (letter case ignored, as a Windows disk does)."""
    base = re.sub(r"[^\w.-]", "_", name) or "sequence"
    out, n = base, 1
    while out.lower() in taken:
        n += 1
        out = f"{base}_{n}"
    taken.add(out.lower())
    return out


def write_records(path, records, fasta):
    """A new file `path` with [(name, Document)] as records: GenBank (features and topology kept) or FASTA (only the bases)."""
    import multigb
    Path(path).write_bytes(b"")
    if fasta:
        multigb.append_fasta(path, records)
    else:
        multigb.append_genbank(path, [multigb.renamed_genbank(d, n) for n, d in records])


def cmd_split(args):
    from reference import is_fasta_file
    names = C.file_records(args.file) if Path(args.file).is_file() else None
    if names is None:
        raise UserError(f"{args.file} holds one record (or is not found): nothing to split.")
    fasta = is_fasta_file(args.file)
    out = Path(args.out) if args.out else Path(args.file).with_name(Path(args.file).stem + "_split")
    if out.exists() and any(out.iterdir()):
        raise UserError(f"The output folder {out} already exists and is not empty: choose a new one with --out (nothing is overwritten).")
    out.mkdir(parents=True, exist_ok=True)
    taken, files, warnings = set(), [], []
    for name in names:
        doc = C.load(f"{args.file}#{name}")
        stem = file_stem(name, taken)
        if stem != name:
            warnings.append(f"record {name!r} was written as {stem}{'.fasta' if fasta else '.gb'} (a file name cannot hold some characters, or two records would share it)")
        write_records(out / (stem + (".fasta" if fasta else ".gb")), [(name, doc)], fasta)
        files.append({"record": name, "file": f"{stem}{'.fasta' if fasta else '.gb'}", "length": doc.length, "circular": doc.circular})
    result = answer({"split": {"format": {"value": "fasta" if fasta else "genbank"}}}, out=str(out), files=files)
    result["warnings"] += warnings
    return result


def cmd_merge(args):
    from sourcefiles import FASTA_EXT
    out = new_path(args.out)
    fasta = out.suffix.lower() in FASTA_EXT
    taken, records, warnings = set(), [], []
    sources, given = records_of(args.files), csv_list(args.names)
    if given and len(given) != len(sources):
        raise UserError(f"--names gives {len(given)} names for {len(sources)} sequences: one name each, in the order of the files.")
    for k, (label, doc) in enumerate(sources):
        name = given[k] if given else re.sub(r"\s+", "_", doc.name or Path(label).stem) or "sequence"
        stem = re.sub(r"\s+", "_", Path(label.split("#")[0]).stem)
        if not given and name.lower() in taken and stem and stem.lower() not in taken:
            name = stem   # (two records with one name: the file's own name tells them apart)
        base, n = name, 1
        while name.lower() in taken:   # (a name two records share gets _2, _3...)
            n += 1
            name = f"{base}_{n}"
        taken.add(name.lower())
        if name != doc.name and not (given and name == given[k]):
            warnings.append(f"{label}: the record is called {name} in the merged file (its name was empty, had spaces, or another record has it)")
        records.append((name, doc, label))
    if len(records) < 2:
        raise UserError("merge needs two or more sequences: give several files (a file of several records counts for all of them).")
    write_records(out, [(n, d) for n, d, _ in records], fasta)
    if fasta:
        warnings.append("FASTA keeps only the bases: features and circular topology are not written (use a .gb file name to keep them).")
    result = answer({"merge": {"format": {"value": "fasta" if fasta else "genbank"}}}, file=str(out), bytes=out.stat().st_size,
                    records=[{"name": n, "source": lb, "length": d.length, "circular": d.circular} for n, d, lb in records])
    result["warnings"] += warnings
    return result


# ---- codon: recode a CDS for a codon usage table (the Recode window) --------------------------------------------------------------

def gc_range(args):
    import codon_tools as CT
    return None if args.no_gc else (CT.DEFAULT_GC[0] if args.gc_low is None else args.gc_low, CT.DEFAULT_GC[1] if args.gc_high is None else args.gc_high)


def codon_row(doc, row, args, usage, avoid, out_dir):
    import codon_tools as CT
    from assembly import snapshot
    from mutation import find_feature
    rec = {"id": row["id"], "feature": row["feature"], "status": "error"}
    try:
        f = find_feature(doc, row["feature"])
        if f.kind != "CDS":
            raise UserError(f"{row['feature']!r} is a {f.kind}, not a CDS.")
        s = CT.feature_source(doc, f)
        if isinstance(s, str):
            raise UserError(s)
        r = CT.recode_tries(s, usage, args.tries, args.seed, method=args.method, avoid=avoid, max_run=CT.MAX_RUN if args.max_run is None else args.max_run,
                             rbs=not args.no_rbs, gc_range=gc_range(args), window=CT.GC_WINDOW if args.window is None else args.window)
    except (UserError, ValueError, RuntimeError) as e:
        rec["error"] = str(e)
        return rec
    rec.update(status="ok", **C.span_fields(doc, s.start + 1, s.end), strand=s.strand, genetic_code=s.table, codons=r.codons, codons_changed=r.changed, seed_used=r.seed, gc_before=round(r.gc_before, 1),
               gc_after=round(r.gc_after, 1), warnings=[f"Could not be met: {u}" for u in r.unmet], sequence=r.seq)
    if out_dir is not None:
        made = snapshot(doc, doc.name)
        CT.write_recoded(made, s, r)
        made.write_genbank(str(out_dir / f"{row['file_stem']}.gb"), None)
        rec["file"] = f"{row['file_stem']}.gb"
    return rec


def cmd_codon(args):
    import codon_tools as CT
    import codon_usage
    doc, rows = C.load(args.plasmid), C.read_table(args.targets, "feature")
    C.assign_files(rows)
    tables = codon_usage.all_tables()
    if args.usage and Path(args.usage).is_file():   # a table in a CSV file (the window's Import reads the same file)
        try:
            usage = codon_usage.read_csv(args.usage)
        except ValueError as e:
            raise UserError(f"{args.usage}: {e}")
    else:
        usage = tables.get(args.usage) if args.usage else codon_usage.active()
    if args.usage and usage is None:
        raise UserError(f"No codon usage table called {args.usage!r}. Choices: {', '.join(tables) or 'none (import one in openDNA: Tools > Settings > Codon usage)'}.")
    if usage is None:
        raise UserError("No codon usage table is set: name one with --usage, or choose one in openDNA (Tools > Settings > Codon usage).")
    if args.method not in (CT.MATCHED, CT.HIGHEST):
        raise UserError(f"--method must be {CT.MATCHED!r} or {CT.HIGHEST!r}.")
    avoid = pick_enzymes(csv_list(args.avoid))
    gc = gc_range(args)
    if gc and gc[0] > gc[1]:
        raise UserError("--gc-low is above --gc-high.")
    if not 1 <= args.tries <= 50:
        raise UserError("--tries is between 1 and 50.")
    out_dir = C.make_out_dir(args, args.plasmid, "_codon")
    lim, g = CT.LIMITS, CT.DEFAULT_GC
    settings = C.effective_settings({"codon": {
        "usage_table": {"value": usage.name, "source": "command line" if args.usage else "openDNA settings"}, "method": C.given(args.method if args.method != CT.MATCHED else None, CT.MATCHED),
        "avoid_sites": {"value": [e.name for e in avoid]},
        **({"gc_limits": {"value": "none", "source": "command line: --no-gc"}} if args.no_gc else
           {"gc_low": C.given(args.gc_low, g[0], lim["gc"], "--gc-low"), "gc_high": C.given(args.gc_high, g[1], lim["gc"], "--gc-high"),
            "gc_window": C.given(args.window, CT.GC_WINDOW, lim["window"], "--window")}), "max_run": C.given(args.max_run, CT.MAX_RUN, lim["max_run"], "--max-run"),
        "no_internal_rbs": {"value": not args.no_rbs}, "seed": {"value": args.seed, "source": "command line" if args.seed != 1 else "default (the window's first preview)"},
        "tries": {"value": args.tries, "source": "command line" if args.tries != 1 else "default (one seed, as the window's first preview)"}}})
    settings["codon_usage_table"] = {"value": usage.name, "source": "command line" if args.usage else "openDNA settings"}   # (the table used, not the saved one)
    records = C.run_rows(rows, lambda row: codon_row(doc, row, args, usage, avoid, out_dir), "feature")
    return C.finish(out_dir, settings, records, ["id", "feature", "file", "codons", "codons_changed", "gc_before", "gc_after", "warnings"],
                    lambda r: [[r["id"], r["feature"], r.get("file", ""), r["codons"], r["codons_changed"], r["gc_before"], r["gc_after"], " | ".join(r["warnings"])]],
                    lambda r: bool(r.get("warnings")))


# ---- export: the other file formats (Export menu) ---------------------------------------------------------------------------------

EXPORTS = ("genbank", "fasta", "gff3", "bed", "tbl", "proteins")


def cmd_export(args):
    import exports
    doc = C.load(args.file)
    out = new_path(args.out)
    ed = doc_editor(doc)
    if args.format == "genbank":
        doc.write_genbank(str(out), None)   # (also how a Benchling, APE or EMBL file becomes a plain GenBank file)
    elif args.format == "fasta":
        out.write_text(exports.sequence_fasta(ed), encoding="utf-8")
    else:
        out.write_text({"gff3": exports.gff3_text, "bed": exports.bed_text, "tbl": exports.tbl_text, "proteins": exports.proteins_fasta}[args.format](ed), encoding="utf-8")
    return answer({"export": {"format": {"value": args.format}}}, file=str(out), bytes=out.stat().st_size, features=len(ed.features))


# ---- order: a primer table as a supplier's order file (Primers > Export order) --------------------------------------------------------

def cmd_order(args):
    import primer_db
    rows = C.read_table(args.primers, "sequence")
    fmt = {"idt": primer_db.ORDER_FORMATS[0], "plain": primer_db.ORDER_FORMATS[1]}[args.format]
    problems = [f"line {r['_line']} ({r['id']}): " + why for r in rows if (why := r.get("_error") or primer_db.primer_problem(r["id"], r["sequence"].upper()))]
    if problems:
        raise UserError("nothing was written, the table has problems: " + "; ".join(problems))
    if args.scale not in primer_db.SCALES or args.purification not in primer_db.PURIFICATIONS:
        raise UserError(f"--scale is one of {', '.join(primer_db.SCALES)}; --purification one of {', '.join(primer_db.PURIFICATIONS)}.")
    out = new_path(args.out)
    written = primer_db.export_order([{"name": r["id"], "sequence": r["sequence"]} for r in rows], out, fmt, args.scale, args.purification)
    return answer({"order": {"format": {"value": fmt}, "scale": {"value": args.scale}, "purification": {"value": args.purification}}}, primers=len(rows), files=[str(w) for w in written])


# ---- variants: where sequencing reads (a BAM) differ from the plasmid (the Variants pane) -----------------------------------------

VARIANT_FLAGS = {"min_depth": "--min-depth", "variant_pct": "--variant-pct", "consensus_pct": "--consensus-pct", "min_baseq": "--min-baseq", "min_mapq": "--min-mapq"}


def cmd_variants(args):
    import alignments
    import track_check
    import variants as V
    from opendna_core import Bam
    doc = C.load(args.file)
    if not Path(args.reads).is_file():
        raise UserError(f"File not found: {args.reads}")
    try:
        bam = Bam(str(args.reads))
    except Exception as e:   # (the core raises OSError with its own words)
        raise UserError(f"Cannot read {args.reads} as a BAM file: {e}")
    lengths = dict(bam.sequences())
    if not lengths:
        raise UserError(f"{args.reads} has no sequence in its header.")
    seqid = args.seqid
    if seqid is None:
        seqid = track_check.match_contig(list(lengths), doc.name) or (next(iter(lengths)) if len(lengths) == 1 else None)   # (as the window picks it)
        if seqid is None:
            raise UserError(f"{args.reads} holds reads for several sequences ({', '.join(lengths)}): name the one that is this plasmid with --seqid.")
    elif seqid not in lengths:
        raise UserError(f"{args.reads} has no sequence called {seqid!r} (it has: {', '.join(lengths)}).")
    warnings = track_check.mismatches(SimpleNamespace(fasta_record=None, doc=doc, length=doc.length), None, seqid, lengths[seqid], args.reads)
    with tool_settings_block(alignments.TOOL, VARIANT_FLAGS, args) as settings:
        values = tuple(settings[k]["value"] for k in ("min_depth", "variant_pct", "consensus_pct", "min_baseq", "min_mapq"))
        try:
            found = alignments.scan_variants(bam, seqid, doc.sequence(0, None), values, doc.circular)
        except OSError as e:
            raise UserError(f"Cannot read the reads of {args.reads}: {e}")
    out = answer({"variants": settings}, reads=str(args.reads), sequence=seqid, reference_length=lengths[seqid], variants=[
        {"position": v.pos + 1, "kind": v.kind, "reference": v.ref, "alternative": v.alt, "reads": v.count, "depth": v.depth, "fraction": round(v.count / max(1, v.depth), 3),
         "level": v.level, "change": V.change_text(v), "can_apply": V.can_apply(v)} for v in found])
    out["warnings"] += warnings
    if args.apply:
        chosen = [v for v in found if v.level == "consensus" and V.can_apply(v)]
        target = new_path(args.apply)
        V.apply_to_document(doc, chosen)   # (the document was read for this command alone)
        doc.write_genbank(str(target), None)
        out["applied"] = {"file": str(target), "changes": len(chosen),
                          "note": "only the consensus changes (most reads disagree with the plasmid) that can be written in; after an insertion or deletion the coordinates move: map the reads again to check" if any(v.kind != "snp" for v in chosen) else "substitutions only: the coordinates did not move"}
    return out


# ---- annotate (Detect Features) and checkup (Plasmid Checkup): the tools that read their thresholds from Settings --------------------

def doc_editor(doc):
    """What the modules that take an editor (annotation, checkup, search, exports) read of one, for a document with no window: no selection,
    tracks or primer sites. A module that reads one more attribute needs it here (tests/test_agent_cli.py runs every command)."""
    from editor import HIDDEN_KINDS
    return SimpleNamespace(doc=doc, length=doc.length, circular=doc.circular, features=doc.features(), groups={}, fasta_record=None, tracks=[], primer_sites=[], start=0, end=0,
                           visible_features=lambda start, end, overlays=True: [f for f in doc.features_in(start, end) if f.kind not in HIDDEN_KINDS])


@contextmanager
def tool_settings_block(tool, flags, args):
    """{key: {value, source, min, max}} for the settings of `tool` that have a command-line flag (`flags`: key -> flag; the value is `args.<key>`). A value
    given on the command line is used for the length of the block only (tool_settings.override: it is never saved); otherwise the user's saved value."""
    import tool_settings
    saved = tool_settings.store()
    out = {}
    for p in tool_settings.TOOLS[tool]:
        if p.key not in flags:
            continue
        given = getattr(args, p.key)
        source = "saved in openDNA (Settings)" if saved.contains(f"{tool}/{p.key}") else "default"
        if p.minimum is None:   # a switch (True/False)
            out[p.key] = {"value": tool_settings.value(tool, p.key) if given is None else given, "source": "command line" if given is not None else source}
        else:
            out[p.key] = C.given(given, tool_settings.value(tool, p.key), (p.minimum, p.maximum), flags[p.key])
            if given is None:
                out[p.key]["source"] = source
    with tool_settings.override({(tool, k): v["value"] for k, v in out.items()}):
        yield out


ANNOTATE_FLAGS = {"min_dna": "--min-dna", "min_protein": "--min-protein", "min_length": "--min-length", "use_plannotate": "--no-plannotate"}


def hit_fields(h, n):
    """One hit of the Detect Features list as the command line prints it. `found_by`: openDNA, pLannotate or both (then `plannotate` is its own hit)."""
    out = {"name": h.name, "kind": h.kind, "start": h.start + 1, "end": end1(h.end, n), "strand": h.strand, "identity": round(h.identity, 4), "protein_match": bool(h.protein),
           "already_annotated": h.known, "annotated_other_way": h.reversed is not None, "found_by": h.origin, "source": h.source, "strength": h.measured}
    if h.also:
        out["plannotate"] = {"name": h.also.name, "database": h.also.database, "identity": round(h.also.identity, 4)}
    return out


def cmd_annotate(args):
    import annotation
    import plannotate_tool
    doc = C.load(args.file)
    ed = doc_editor(doc)
    args.use_plannotate = False if args.no_plannotate else None   # (a switch left alone keeps the saved value)
    with tool_settings_block(annotation.TOOL, ANNOTATE_FLAGS, args) as settings:
        hits, plannotate = annotation.detect_all(ed, set(csv_list(args.kinds)) or None)
    n = doc.length
    out = answer({"detect_features": settings}, parts_in_library=len(annotation.load_parts()), plannotate=plannotate.name, hits=[hit_fields(h, n) for h in hits])
    if plannotate.name == "ran":
        out["plannotate_citation"] = f"{plannotate_tool.CITATION} doi:{plannotate_tool.DOI}"   # the user who publishes a result cites it
    if plannotate.text:
        out["warnings"].append(f"pLannotate was not used: {plannotate.text}" + (" Install it with: oDNA-cli wsltools install --only plannotate" if plannotate.name == "missing" else ""))
    if args.add:
        target = new_path(args.add)
        fresh = [h for h in hits if not h.known and h.reversed is None]
        turn = [h for h in hits if h.reversed is not None and h.kind in annotation.TURN_KINDS]   # (as the window ticks them)
        annotation.turn_features(doc, turn)   # (before the new ones are added: the feature indexes are the ones read)
        if fresh:
            doc.add_features(annotation.features_for(ed, fresh))
        doc.write_genbank(str(target), None)
        out["added"] = {"file": str(target), "features": len(fresh), "turned_round": len(turn), "skipped_already_annotated": len(hits) - len(fresh) - len(turn)}
    return out


CHECKUP_FLAGS = {"min_promoter": "--min-promoter", "min_stem": "--min-stem", "rbs_reach": "--rbs-reach"}


def cmd_checkup(args):
    import checkup
    doc = C.load(args.file)
    args.predict, args.library = (False if args.no_predict else None), (False if args.no_library else None)   # (a switch left alone keeps the saved value)
    with tool_settings_block(checkup.TOOL, {**CHECKUP_FLAGS, "predict": "--no-predict", "library": "--no-library"}, args) as settings:
        res = checkup.run(doc_editor(doc))
    from mutation import arc
    n = doc.length
    span = lambda el: {"start": el.start + 1, "end": end1(el.end, n), "strand": el.strand}
    cds = []
    for c in res.cds:
        f = c.feature
        cds.append({"cds": f.label, **C.span_fields(doc, *arc(f, n, doc.circular)), "strand": f.strand, "verdict": c.verdict, "best_host": c.host, "reasons": c.reasons,
                    "hosts": {h: {"verdict": r.verdict, "promoter": r.unit.promoter.label if r.unit else None, "terminator": r.unit.terminator.label if r.unit and r.unit.terminator else None,
                                  "rbs": r.rbs.element.label if r.rbs else None, "reasons": r.reasons} for h, r in c.hosts.items()}})
    return answer({"checkup": settings}, note=res.note, hosts=list(res.hosts), cds=cds,
                  elements=[{"kind": e.kind, "label": e.label, "origin": e.origin, "detail": e.detail, **span(e)} for e in res.elements],
                  origins=[{"label": o.element.label, **span(o.element), "name": o.name, "host": o.host, "copies_per_cell": o.copies, "incompatibility_group": o.group, "note": o.note} for o in res.origins], origin_notes=list(res.origin_notes))


# ---- gibson: Gibson / HiFi assembly of fragments (the Assembly panel's Gibson method) -------------------------------------------------

def cmd_gibson(args):
    import csv
    import assembly as A
    import recombineering as R
    from assembly import Frag, snapshot
    rows = C.read_table(args.fragments, "source", optional=("region", "orientation"))
    bad = [r["_error"] for r in rows if "_error" in r]
    if bad:   # (unlike the batch commands, an assembly needs every fragment: one bad row stops it)
        raise UserError("the fragments cannot be assembled with a row missing: " + " ".join(bad))
    names = [r["id"] for r in rows]
    if len(set(names)) != len(names):
        raise UserError(f"The ids must be different (they name the fragments and their primers): {', '.join(sorted({n for n in names if names.count(n) > 1}))}.")
    if len(rows) < 2:
        raise UserError("An assembly needs at least two fragments.")
    if not re.fullmatch(r"[A-Za-z0-9._-]+", args.name):
        raise UserError(f"--name {args.name!r} is not usable as a file name: letters, digits, . _ - only.")
    folder = Path(args.fragments).resolve().parent
    docs, frags = {}, []
    for r in rows:
        path = Path(r["source"])
        path = str(path if path.is_absolute() else folder / path)
        if path not in docs:
            docs[path] = C.load(path)
        src = docs[path]
        first, last = (1, src.length) if not r.get("region") else C.target_range(src, r["region"])
        orientation = (r.get("orientation") or "forward").lower()
        if orientation not in ("forward", "reverse"):
            raise UserError(f"Row {r['id']!r}: orientation must be forward or reverse, not {orientation!r}.")
        frag = Frag(r["id"], snapshot(src, src.name), first - 1, last, source=Path(path).name)
        if orientation == "reverse":
            frag.flip()
        frags.append(frag)
    lim, dflt = A.GIBSON_LIMITS, A.GIBSON_DEFAULTS
    settings = C.effective_settings({"gibson": {
        "circular": {"value": not args.linear, "source": "command line" if args.linear else "default"},
        "min_overlap": C.given(args.min_overlap, dflt["min_overlap"], lim["min_overlap"], "--min-overlap"), "max_overlap": C.given(args.max_overlap, dflt["max_overlap"], lim["max_overlap"], "--max-overlap"),
        "overlap_tm": C.given(args.overlap_tm, dflt["overlap_tm"], lim["overlap_tm"], "--overlap-tm"), "annealing_tm": C.given(args.tm, R.DEFAULT_TM, R.LIMITS["annealing_tm"], "--tm")}})
    g = {k: v["value"] for k, v in settings["gibson"].items()}
    if g["min_overlap"] > g["max_overlap"]:
        raise UserError("--min-overlap is above --max-overlap.")
    out_dir = C.make_out_dir(args, args.fragments, "_gibson")
    circular = g["circular"]
    A.auto_primers(frags, 0, circular, min_ov=g["min_overlap"], max_ov=g["max_overlap"], min_tm=g["overlap_tm"], tm=g["annealing_tm"])
    try:
        d = A.analyse(frags, 0, circular, min_ov=g["min_overlap"], min_tm=g["overlap_tm"], tm=g["annealing_tm"], max_ov=100)
    except ValueError as e:
        raise UserError(str(e))
    product = A.simulate(frags, d)
    import primer_panel
    primers = [{"fragment": f.name, "name": p[0], "sequence": p[1], "annealing_length": p[2], "tail_nt": len(p[1]) - p[2], "annealing_tm": primer_panel.melting_temp(p[1][-p[2]:]), "role": role}
               for f in frags for role, p in (("forward", f.fwd), ("reverse", f.rev)) if p]
    levels = [lv for _, lv, _ in d.checks]
    verdict = "bad" if not d.ok else "warn" if "warn" in levels else "ok"   # (d.ok: every junction fits and no check is red, as the panel decides)
    result = {"settings": settings, "warnings": C.setting_warnings(settings), "verdict": verdict,
              "fragments": [{"id": f.name, "source": f.source, "start": f.start + 1, "end": f.end, "orientation": "reverse" if f.flipped else "forward", "length": f.core.length} for f in frags],
              "junctions": [{"text": t, "level": lv} for t, lv in zip(d.junctions, d.jlevels)], "checks": [{"check": c, "level": lv, "text": t} for c, lv, t in d.checks],
              "primers": primers, "product": {"length": product.length, "circular": product.circular}}
    if out_dir is not None:
        A.snapshot(product, args.name).write_genbank(str(out_dir / f"{args.name}.gb"), None)
        result["product"]["file"] = f"{args.name}.gb"
        with open(out_dir / "primers.csv", "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow(["fragment", "role", "name", "sequence", "annealing_length", "tail_nt", "annealing_tm"])
            w.writerows([p["fragment"], p["role"], p["name"], p["sequence"], p["annealing_length"], p["tail_nt"], C.tm_of(p)] for p in primers)
        (out_dir / "report.json").write_text(json.dumps(result, indent=1), encoding="utf-8")
        result["out"] = str(out_dir)
    for w in result["warnings"]:
        print(f"WARNING: {w}", file=sys.stderr)
    return result


# ---- ncbi: search and fetch records (Import from NCBI) -------------------------------------------------------------------------------

def cmd_ncbi(args):
    import ncbi
    import shutil
    db = args.db
    if args.action == "get" and not args.out:
        raise UserError("`ncbi get` needs --out NEW_FILE.")
    if args.limit is not None and args.limit < 1:
        raise UserError("--limit must be at least 1.")
    out = new_path(args.out) if args.action == "get" else None   # (before the network is used)
    key = ncbi.saved_key().strip()
    settings = {"ncbi": {"database": {"value": args.db, "source": "command line" if args.db != "nuccore" else "default"},
                         "api_key": {"value": "set" if key else "none", "source": "saved in openDNA (Import from NCBI)" if key else "default: 3 requests a second"},
                         "sends": {"value": f"only {args.query!r} to eutils.ncbi.nlm.nih.gov"}}}
    try:
        found = ncbi.search(db, args.query, key, ncbi.MAX_RESULTS if args.limit is None else args.limit)
        records = [{"accession": r.accession, "title": r.title, "length": r.length, "organism": r.organism, "updated": r.updated, "topology": r.topology, "source": r.source}
                   for r in found]
        if args.action == "search":
            return answer(settings, caps=[ncbi.more_note(found, ncbi.MAX_RESULTS if args.limit is None else args.limit)], records=records)
        q = args.query.strip().lower()
        wanted = [r for r in found if r.accession.lower() == q or ("." not in q and r.accession.lower().split(".")[0] == q)]
        if not wanted:
            raise UserError(f"{args.query!r} is not the accession of a {args.db} record (`ncbi search` lists what matches).")
        record = wanted[0]
        if record.db == "nuccore" and record.length > ncbi.WARN_BASES and not args.allow_large:
            raise UserError(f"{record.accession} has {record.length:,} bases: a download of that size takes minutes; add --allow-large if the user wants it.")
        path, cached = ncbi.get(record, key)
        shutil.copyfile(path, out)
    except ncbi.NcbiError as e:
        raise UserError(str(e))
    return answer(settings, record=records[found.index(record)], file=str(out), bytes=out.stat().st_size, from_cache=cached)


def cmd_hub(args):
    """A UCSC track hub (a URL, a hub.txt or a folder holding one): `list` its tracks (in the hub's groups, with labels, types and
    whether the hub hides them at first), `get` copies the files of the tracks openDNA reads into a new folder. The same reading
    and copying as Import > Track Hub; the window then loads them."""
    import track_hub
    out = None
    if args.action == "get":
        if not args.out:
            raise UserError("`hub get` needs --out NEW_FOLDER.")
        out = Path(args.out)
        if out.exists() and any(out.iterdir()):
            raise UserError(f"The output folder {out} already exists and is not empty: choose a new one with --out (nothing is overwritten).")
    hub = track_hub.open_hub(args.source)
    genome = args.genome or track_hub.pick_genome(hub, [])
    if genome is None:
        raise UserError(f"the hub has several genomes ({', '.join(hub.genomes)}): name one with --genome")
    tracks = track_hub.hub_tracks(hub, genome)
    remote = track_hub.is_url(hub.source)
    settings = {"hub": {"source": {"value": hub.source}, "genome": {"value": genome, "source": "command line" if args.genome else "the hub's only genome"},
                         "sends": {"value": "requests for the hub's own files, to its server" if remote else "nothing: a hub on this computer"}}}
    unread = [t for t in tracks if t.why]
    caps = [f"{len(unread)} track(s) of a type openDNA does not read, listed but not copied: {', '.join(t.label for t in unread[:10])}"
            + (" …" if len(unread) > 10 else "")] if unread else []
    listed = [{"name": t.name, "label": t.label, "group": t.group, "type": t.kind, "hidden": t.hidden, "color": t.color, "url": t.url,
               "read": not t.why, **({"why": t.why} if t.why else {})} for t in tracks]
    if args.action == "list":
        return answer(settings, caps=caps, hub={"name": hub.name, "label": hub.label, "genomes": list(hub.genomes)}, tracks=listed)
    try:
        copied = track_hub.copy_tracks([t for t in tracks if not t.why], out)
    except OSError as e:   # (the folder cannot be made)
        raise UserError(f"cannot write in {out}: {e}") from None
    failed = [{"name": t.name, "label": t.label, "error": why} for t, path, why in copied if path is None]
    caps += [f"{f['label']}: not copied ({f['error']})" for f in failed]
    files = [{"name": t.name, "label": t.label, "group": t.group, "hidden": t.hidden, "file": str(path), "bytes": path.stat().st_size}
             for t, path, _ in copied if path is not None]
    return answer(settings, caps=caps, folder=str(out), files=files, failed=failed)


# ---- tracks: recap, print, track-features, signal-regions (the Alignment view's actions on regions and tracks) -----------------

def track_editor(source, tracks=()):
    """An Editor with no window for `source` (a sequence file, or an Open DNA Project with its display, tracks, their order,
    groups, look and highlighted regions, read as the window reads them), with the track files `tracks` added at the end:
    (editor, [warnings]). What `print` and `recap` draw and count is what the window shows."""
    import exports
    import track_files
    import track_import
    from editor import Editor
    warnings, meta = [], None
    if exports.is_project(source):
        try:
            doc, meta = exports.open_project(source)
        except OSError as e:
            raise UserError(f"{source}: {e}") from None
    else:
        doc = C.load(source)
    ed = Editor(doc, source)
    if meta is not None:
        ed.apply_project_display(meta)
        for saved, found in exports.project_tracks(meta, source):
            if found is None:
                warnings.append(f"{saved['path']}: a track of the project, not found (moved or deleted)")
            else:
                warnings += track_files.load(ed, found, saved.get("name"), saved.get("seqid"), saved.get("types", track_import.ALL),
                                             saved.get("fields", track_import.ALL))[1]
        ed.track_layout.load_json(meta.get("track_layout"), ed.row_keys(), ed.length)
    for path in tracks or ():
        if not Path(path).is_file():
            raise UserError(f"Track file not found: {path}")
        warnings += track_files.load(ed, path)[1]
    return ed, warnings


def spans_of(editor, text):
    """[(start, end)] (0-based, end excluded) of START..END or a feature's name (cli.target_range): two parts through the origin."""
    from editor import split_at_origin
    first, last = C.target_range(editor.doc, text)
    return split_at_origin(first - 1, last, editor.length)


def named(regions):
    """[(name, start, end)] of highlights (track_layout.Region), an unnamed one named by its place."""
    return [(r.name or f"{r.start + 1}..{r.end}", r.start, r.end) for r in regions]


def cmd_recap(args):
    """Every loaded track over regions (track_recap: what Recap Tracks Here shows): the place given, the regions of a BED file,
    else the project's highlighted regions."""
    import track_layout
    import track_recap
    ed, warnings = track_editor(args.source, args.track)
    if not (args.where or args.regions):
        regions = named(ed.track_layout.regions)
    else:
        regions = [(args.where, s, e) for s, e in spans_of(ed, args.where)] if args.where else []
    if args.regions:
        try:
            found, unread = track_layout.read_regions_bed(Path(args.regions).read_text(encoding="utf-8", errors="replace"), ed.length)
        except OSError as e:
            raise UserError(f"cannot read {args.regions}: {e}") from None
        regions += named(found)
        warnings += [f"{args.regions}: {u}" for u in unread]
    if not regions:
        raise UserError("no region: give START..END (or a feature's name), --regions FILE.bed, or a project with highlighted regions.")
    if not ed.loaded():
        warnings.append("no track loaded: add --track FILE, or give a project that has tracks")
    long = [n for n, s, e in regions if e - s > track_recap.PER_BASE]
    if args.per_base and long:
        warnings.append(f"--per-base: values base by base only for regions up to {track_recap.PER_BASE:,} bp; not for {', '.join(long[:5])}")
    recaps = track_recap.recap(ed, regions, args.per_base)
    if args.tsv:
        for w in warnings:
            print(f"WARNING: {w}", file=sys.stderr)
        args.text = True
        return track_recap.tsv(track_recap.rows(recaps)).rstrip("\n")
    return answer({"recap": {"source": {"value": str(args.source)}, "per_base": {"value": bool(args.per_base)}}}, caps=warnings, regions=recaps)


def cmd_print(args):
    """The Alignment view over START..END as an image file (genome_view.write_region: what Export This Region as Image writes)."""
    from genome_view import write_region
    ed, warnings = track_editor(args.source, args.track)
    spans = spans_of(ed, args.where)
    if len(spans) > 1:
        raise UserError(f"{args.where} runs through the origin: print the two parts one at a time.")
    if not 200 <= args.width <= 20_000:
        raise UserError("--width is in pixels, from 200 to 20,000.")
    out = new_path(args.out)
    try:
        write_region(ed, *spans[0], out, args.width)
    except OSError as e:
        raise UserError(str(e)) from None
    shown = [k for k in ed.track_layout.keys() if ed.track_of(k) is not None and ed.track_layout.style(k).shown]
    return answer({"print": {"width": {"value": args.width}}}, caps=warnings, file=str(out), bytes=out.stat().st_size,
                  start=spans[0][0] + 1, end=spans[0][1], tracks_shown=len(shown), regions=len(ed.track_layout.regions))


def cmd_track_features(args):
    """A track's items as features of the sequence, written to a new GenBank file (track_features: what Copy to Features adds)."""
    import track_features
    ed, warnings = track_editor(args.genome, [args.track_file])
    entry = next(iter(ed.loaded()), None)
    if entry is None or entry.kind != "track":
        raise UserError(f"{args.track_file} is not an annotation track of this sequence (GFF, GTF, BED, bigBed, VCF): " + "; ".join(warnings))
    label = Path(args.track_file).name
    kind = args.type or track_features.default_type(label)
    if args.name_in not in track_features.NAME_FIELDS:
        raise UserError(f"--name-in is one of {', '.join(track_features.NAME_FIELDS)}.")
    spans = spans_of(ed, args.where) if args.where else [(0, ed.length)]
    out = new_path(args.out)
    feats, skipped = track_features.to_features(track_features.items_in(entry, spans), label, kind, args.name_in,
                                                track_features.feature_keys(ed.features), ed.length)
    if feats:
        ed.doc.add_features(feats)
    ed.doc.write_genbank(str(out), None)
    return answer({"track-features": {"type": {"value": kind, "source": "command line" if args.type else "from the track's name"},
                                      "name_in": {"value": args.name_in}}},
                  caps=warnings + ([f"{skipped} item(s) already features of the sequence (same type, place and strand): not added again"] if skipped else []),
                  file=str(out), added=len(feats), skipped=skipped, features=len(ed.doc.features()))


def cmd_signal_regions(args):
    """The stretches where a signal (bigWig, bedGraph, wig) is at least a threshold, as JSON and optionally BED (track_features:
    what Regions Above Threshold finds)."""
    import track_features
    import track_files
    if not Path(args.signal).is_file():
        raise UserError(f"File not found: {args.signal}")
    if track_files.track_kind(args.signal) != "signal":
        raise UserError(f"{args.signal} is not a signal track: give a bigWig, bedGraph or wig file.")
    try:
        signal = track_files.read(args.signal, "signal")
    except OSError as e:
        raise UserError(f"cannot read {args.signal}: {e}") from None
    lengths = dict(signal.sequences())
    seq = args.seq or (next(iter(lengths)) if len(lengths) == 1 else None)
    if seq not in lengths:
        raise UserError(f"name the sequence with --seq: one of {', '.join(list(lengths)[:10])}.")
    start, end = 0, lengths[seq]
    if args.where:
        m = re.fullmatch(r"(\d+)\.\.(\d+)", args.where)
        if not m or not 1 <= int(m.group(1)) <= int(m.group(2)) <= lengths[seq]:
            raise UserError(f"{args.where!r}: give START..END with 1 <= START <= END <= {lengths[seq]:,}.")
        start, end = int(m.group(1)) - 1, int(m.group(2))
    if args.min_length < 1 or args.gap < 0:
        raise UserError("--min-length is at least 1 and --gap at least 0 (bp).")
    out = new_path(args.out) if args.out else None
    found = track_features.signal_regions(signal, seq, start, end, args.threshold, args.min_length, args.gap)
    if out:
        out.write_text(track_features.regions_bed(found, seq, Path(args.signal).stem), encoding="utf-8")
    return answer({"signal-regions": {"threshold": {"value": args.threshold}, "min_length": {"value": args.min_length}, "gap": {"value": args.gap}}},
                  caps=signal.notes, sequence=seq, regions=[{"start": s + 1, "end": e, "length": e - s, "mean": float(f"{mean:.6g}"),
                                                             "max": float(f"{top:.6g}")} for s, e, mean, top in found],
                  **({"file": str(out)} if out else {}))


def cmd_wsltools(args):
    """The tools that run beside openDNA (inside WSL on Windows): pLannotate and the Terminal's. `install` is what the Windows installer's
    optional task runs; it installs WSL itself never (that needs an administrator and a restart: the user does it)."""
    import plannotate_tool as P
    state, tools = P.tools_state()
    if state.state == "no_wsl":
        if args.action == "install":
            raise UserError(P.WSL_HELP)
        return answer({}, wsl=False, plannotate={"state": "no_wsl"}, terminal_tools=tools, note=P.WSL_HELP)
    if args.action == "install":
        code = P.install(plannotate=args.only != "terminal", terminal=args.only != "plannotate")
        if code != 0:
            raise UserError(f"The install did not finish (exit code {code}); the text above says why. It is safe to run it again.")
        state, tools = P.tools_state()
    return answer({}, wsl=True, plannotate={"state": state.state, "version": state.version}, terminal_tools=tools)


def register(sub, cli):
    """Add the commands to the argument parser `sub`; `cli` is the cli module."""
    global C, UserError
    C, UserError = cli, cli.UserError
    s = sub.add_parser("sites", help="restriction sites (unique 6+ cutters by default) and the fragments of a digest")
    s.add_argument("file")
    s.add_argument("--enzymes", help="comma-separated names: only these (counted whatever their number of sites)")
    s.add_argument("--min-sites", type=int, help="list enzymes that cut at least this many times (default 1)")
    s.add_argument("--max-sites", type=int, help="list enzymes that cut at most this many times (default 1: unique cutters)")
    s.add_argument("--digest", help="comma-separated enzymes cut together: the fragments")
    s.set_defaults(run=cmd_sites)
    s = sub.add_parser("guides", help="Cas9/Cas12a guides ranked, with their off-target counts and a mark (Spacer Finder)")
    s.add_argument("file")
    s.add_argument("--nuclease", default="SpCas9", help="SpCas9 (default), SaCas9 or Cas12a")
    s.add_argument("--region", help="START..END (1-based) to look in (default: the whole sequence)")
    s.add_argument("--max-mm", type=int, help="count other sites up to this many mismatches (default 3)")
    s.add_argument("--top", type=int, help="how many of the best guides to report and check (default: 300)")
    s.add_argument("--no-genome", action="store_true", help="do not count sites in the background genome")
    s.set_defaults(run=cmd_guides)
    s = sub.add_parser("repeats", help="stretches shared by two files (or by one with itself): direct and inverted repeats")
    s.add_argument("file")
    s.add_argument("other", nargs="?", help="a second file (default: the first one again)")
    s.add_argument("--min-len", type=int, help="shortest shared stretch in bp (default 12)")
    s.add_argument("--mismatches", type=int, help="differences allowed inside a stretch (default 0)")
    s.set_defaults(run=cmd_repeats)
    s = sub.add_parser("synthesis", help="Synthesis Success Calculator (Halper, Hossain & Salis 2020): will a synthesis company make this fragment quickly? Needs ViennaRNA")
    s.add_argument("file")
    s.add_argument("--region", help="START..END (1-based) to test (default: the whole sequence)")
    s.set_defaults(run=cmd_synthesis)
    s = sub.add_parser("search", help="a DNA sequence (with IUPAC codes), a protein motif or a name, in the sequence and its features")
    s.add_argument("file")
    s.add_argument("query")
    s.add_argument("--fuzzy", type=int, help="mismatches allowed (default 0)")
    s.add_argument("--single-strand", action="store_true", help="do not look for the reverse complement")
    s.add_argument("--protein", action="store_true", help="also look for a protein motif in the six-frame translations")
    s.add_argument("--details", action="store_true", help="names match any qualifier of a feature, not just its name")
    s.add_argument("--enzymes", action="store_true", help="an enzyme name finds its sites")
    s.set_defaults(run=cmd_search)
    s = sub.add_parser("align", help="where a query (bases, or a file) matches the plasmid: pieces, strands, mismatches, structural changes (Local Align)")
    s.add_argument("file")
    s.add_argument("query", help="DNA text, or a GenBank/FASTA file")
    s.add_argument("--time-limit", type=float, help="seconds for the search, 0 = none (default: the window's setting, Local Align > Time limit)")
    s.add_argument("--threads", type=int, dest="cores", help="CPU cores to use (default: the window's setting, General > Performance)")
    s.set_defaults(run=cmd_align)
    s = sub.add_parser("split", help="write one file per record of a multi-FASTA or multi-record GenBank file (new folder)")
    s.add_argument("file")
    s.add_argument("--out", help="the new folder (default: next to the file, FILE_split)")
    s.set_defaults(run=cmd_split)
    s = sub.add_parser("merge", help="join sequences (files, or FILE#NAME records) into one multi-FASTA (.fasta) or multi-record GenBank (.gb) file")
    s.add_argument("files", nargs="+")
    s.add_argument("--names", help="the record names, comma separated, one per sequence in order (default: the record's own name; on a clash, the file's name)")
    s.add_argument("--out", required=True, help="the new file; its extension picks the format (.fasta/.fa keeps only the bases, anything else is GenBank)")
    s.set_defaults(run=cmd_merge)
    s = sub.add_parser("msa", help="multiple alignment of two or more sequences (files, or one file of several records), with identity and consensus")
    s.add_argument("files", nargs="+")
    s.add_argument("--out", help="also write the alignment as FASTA to this new file")
    s.set_defaults(run=cmd_msa)
    s = sub.add_parser("codon", help="recode CDS features for a codon usage table keeping the protein (check: report only; run: also write plasmids)")
    s.add_argument("action", choices=["check", "run"])
    s.add_argument("plasmid")
    s.add_argument("targets", help="CSV with the columns id,feature")
    s.add_argument("--out", help="run: a new output folder (default: PLASMID_codon next to the plasmid)")
    s.add_argument("--usage", help="a codon usage table: its name in openDNA, or a CSV file (default: the one chosen in openDNA)")
    s.add_argument("--method", default="Frequency matched", help='"Frequency matched" (default) or "Highest frequency"')
    s.add_argument("--avoid", help="comma-separated enzymes whose sites must not appear")
    s.add_argument("--gc-low", type=int, help="lowest GC percent in every window (default 35)")
    s.add_argument("--gc-high", type=int, help="highest GC percent (default 65)")
    s.add_argument("--no-gc", action="store_true", help="no GC limits")
    s.add_argument("--window", type=int, help="bases per GC window (default 50)")
    s.add_argument("--max-run", type=int, help="longest run of one base (default 5)")
    s.add_argument("--no-rbs", action="store_true", help="allow an internal ribosome binding site")
    s.add_argument("--seed", type=int, default=1, help="the random draw of a frequency-matched recoding (default 1; the same seed gives the same sequence)")
    s.add_argument("--tries", type=int, default=1, help="try this many seeds (seed, seed+1, ...) and keep the recoding with the fewest rules not met (default 1; at most 50)")
    s.set_defaults(run=cmd_codon, returns_code=True)
    s = sub.add_parser("export", help="write the sequence as genbank, fasta, gff3, bed, tbl (NCBI feature table) or proteins (FASTA of the CDS)")
    s.add_argument("file")
    s.add_argument("format", choices=EXPORTS)
    s.add_argument("out", help="a new file (never overwritten)")
    s.set_defaults(run=cmd_export)
    s = sub.add_parser("order", help="a primer table (columns id,sequence) as an IDT plate or plain name,sequence order file")
    s.add_argument("primers")
    s.add_argument("out", help="a new file (never overwritten)")
    s.add_argument("--format", choices=["idt", "plain"], default="idt")
    s.add_argument("--scale", default="25nm")
    s.add_argument("--purification", default="STD")
    s.set_defaults(run=cmd_order)
    s = sub.add_parser("variants", help="where the reads of a BAM file differ from the plasmid (SNPs, insertions, deletions, breaks); --apply writes a corrected copy")
    s.add_argument("file")
    s.add_argument("reads", help="a sorted, indexed BAM file mapped to this plasmid")
    s.add_argument("--seqid", help="the BAM's name for this plasmid (default: the one that matches its name, or the only one)")
    s.add_argument("--min-depth", dest="min_depth", type=int, help="fewest reads to count a position (default: the Alignments setting, 5)")
    s.add_argument("--variant-pct", dest="variant_pct", type=int, help="percent of reads that must differ (setting default 20)")
    s.add_argument("--consensus-pct", dest="consensus_pct", type=int, help="percent of reads that make a change a consensus change (setting default 80)")
    s.add_argument("--min-baseq", dest="min_baseq", type=int, help="lowest base quality counted (setting default 13)")
    s.add_argument("--min-mapq", dest="min_mapq", type=int, help="lowest mapping quality counted (setting default 0)")
    s.add_argument("--apply", help="write a copy of the plasmid with the consensus changes applied to this new GenBank file")
    s.set_defaults(run=cmd_variants)
    s = sub.add_parser("annotate", help="find known parts (promoters, terminators, origins, resistance genes...) in the plasmid (Detect Features); --add writes a copy with them annotated")
    s.add_argument("file")
    s.add_argument("--kinds", help="comma-separated feature types to look for (promoter, terminator, CDS, rep_origin...)")
    s.add_argument("--min-dna", dest="min_dna", type=int, help="minimum DNA identity in percent (the Detect Features setting, 95)")
    s.add_argument("--min-protein", dest="min_protein", type=int, help="minimum protein identity of a CDS in percent (90)")
    s.add_argument("--min-length", dest="min_length", type=int, help="ignore parts shorter than this many bp (15)")
    s.add_argument("--no-plannotate", action="store_true", help="search only openDNA's Features Database, not pLannotate (which is used when it is installed)")
    s.add_argument("--add", help="write a copy of the plasmid with the parts not annotated yet added as features, to this new GenBank file")
    s.set_defaults(run=cmd_annotate)
    s = sub.add_parser("wsltools", help="status of pLannotate and of the Terminal's tools (WSL on Windows); install: set them up (several hundred MB, a few minutes)")
    s.add_argument("action", choices=["status", "install"])
    s.add_argument("--only", choices=["plannotate", "terminal"], help="install: only pLannotate (with its databases), or only the Terminal's tools (default: both)")
    s.set_defaults(run=cmd_wsltools)
    s = sub.add_parser("checkup", help="for each CDS: should it be expressed? promoters, RBS, terminators, transcription units (Plasmid Checkup)")
    s.add_argument("file")
    s.add_argument("--min-promoter", dest="min_promoter", type=int, help="predicted promoter: minimum score (the Plasmid Checkup setting, 10)")
    s.add_argument("--min-stem", dest="min_stem", type=int, help="predicted terminator: minimum stem score (18)")
    s.add_argument("--rbs-reach", dest="rbs_reach", type=int, help="an RBS counts within this many bp of the start codon (20)")
    s.add_argument("--no-predict", action="store_true", help="use only annotated features and library parts, no motif prediction")
    s.add_argument("--no-library", action="store_true", help="do not use library parts")
    s.set_defaults(run=cmd_checkup)
    s = sub.add_parser("gibson", help="Gibson / HiFi assembly of fragments in a table (source file, region, orientation): primers with overlaps, checks and the assembled plasmid")
    s.add_argument("action", choices=["check", "run"])
    s.add_argument("fragments", help="CSV with the columns id,source (and optionally region,orientation), one row per fragment, in order")
    s.add_argument("--out", help="run: a new output folder (default: next to the table, FRAGMENTS_gibson)")
    s.add_argument("--name", default="assembly", help="run: name of the assembled plasmid's file (default assembly)")
    s.add_argument("--linear", action="store_true", help="the product is linear (default: circular, the last fragment joins the first)")
    s.add_argument("--min-overlap", dest="min_overlap", type=int, help="shortest overlap in bp (default 20)")
    s.add_argument("--max-overlap", dest="max_overlap", type=int, help="longest overlap in bp (default 40)")
    s.add_argument("--overlap-tm", dest="overlap_tm", type=float, help="the overlap's lowest Tm in deg C (default 52)")
    s.add_argument("--tm", type=float, help="annealing Tm of the primers in deg C (default 60)")
    s.set_defaults(run=cmd_gibson)
    s = sub.add_parser("ncbi", help="search NCBI Nucleotide or Protein, or fetch one record as a GenBank (nucleotide) or FASTA (protein) file; sends only the query to NCBI")
    s.add_argument("action", choices=["search", "get"])
    s.add_argument("query", help="search words, or (get) an accession such as NC_000913.3 or J01749")
    s.add_argument("--db", choices=["nuccore", "protein"], default="nuccore")
    s.add_argument("--limit", type=int, help="search: how many records (default 25)")
    s.add_argument("--out", help="get: the new file to write")
    s.add_argument("--allow-large", action="store_true", help="get: a record of more than 20 million bases")
    s.set_defaults(run=cmd_ncbi)
    s = sub.add_parser("hub", help="a UCSC track hub: list its tracks (groups, labels, types, hidden at first), or copy their files into a new folder")
    s.add_argument("action", choices=["list", "get"])
    s.add_argument("source", help="the hub: a URL of its hub.txt, a hub.txt, or a folder holding one")
    s.add_argument("--genome", help="which genome of the hub (needed when it has several)")
    s.add_argument("--out", help="get: a new folder for the files")
    s.set_defaults(run=cmd_hub)
    tracks_help = "a track file to add (annotation, signal or BAM); several --track for several"
    s = sub.add_parser("recap", help="what every track holds over a region (items, signal mean / min / max / sum, depth), as Recap Tracks Here")
    s.add_argument("source", help="a sequence file, or an Open DNA Project (.oDNA) with its tracks and highlighted regions")
    s.add_argument("where", nargs="?", help="START..END (1-based) or a feature's name; none: --regions, else the project's highlighted regions")
    s.add_argument("--track", action="append", help=tracks_help)
    s.add_argument("--regions", help="a BED file of regions to recap, each on its own")
    s.add_argument("--per-base", dest="per_base", action="store_true", help="signals: the value at every base too (regions up to 1,000 bp)")
    s.add_argument("--tsv", action="store_true", help="one tab-separated row per region and track (and per group: the mean of its tracks) instead of JSON")
    s.set_defaults(run=cmd_recap)
    s = sub.add_parser("print", help="the Alignment view over START..END as an SVG or PNG file, as the project shows its tracks")
    s.add_argument("source", help="an Open DNA Project (.oDNA: its tracks, order, looks, highlights) or a sequence file with --track")
    s.add_argument("where", help="START..END (1-based) or a feature's name")
    s.add_argument("out", help="the new image file: .svg (vector), .png or .tif")
    s.add_argument("--track", action="append", help=tracks_help)
    s.add_argument("--width", type=int, default=1600, help="the image's width in pixels (default 1600)")
    s.set_defaults(run=cmd_print)
    s = sub.add_parser("track-features", help="a track's items as features of the sequence, in a new GenBank file (as Copy to Features)")
    s.add_argument("genome", help="the sequence (GenBank)")
    s.add_argument("track_file", metavar="track", help="the annotation track (GFF, GTF, BED, bigBed, VCF)")
    s.add_argument("where", nargs="?", help="only the items overlapping START..END (1-based) or a feature (default: the whole sequence)")
    s.add_argument("--out", required=True, help="the new GenBank file")
    s.add_argument("--type", help="the features' type (default: CDS for a track named *_CDS, else misc_feature)")
    s.add_argument("--name-in", dest="name_in", default="label", help="the qualifier the item's name goes in (label, gene, locus_tag, product, note)")
    s.set_defaults(run=cmd_track_features)
    s = sub.add_parser("signal-regions", help="the stretches where a signal track (bigWig, bedGraph, wig) is at least a threshold, optionally as BED")
    s.add_argument("signal")
    s.add_argument("--threshold", type=float, required=True, help="the lowest value of a stretch")
    s.add_argument("--min-length", dest="min_length", type=int, default=1, help="leave out stretches shorter than this (bp, default 1)")
    s.add_argument("--gap", type=int, default=0, help="join stretches this close or closer (bp, default 0)")
    s.add_argument("--seq", help="the sequence of the file (needed when it has several)")
    s.add_argument("--region", dest="where", help="only START..END (1-based)")
    s.add_argument("--out", help="a new BED file of the stretches")
    s.set_defaults(run=cmd_signal_regions)
