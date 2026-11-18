"""Command line for agents (and scripts): the same functions as the windows, with no window.

    oDNA-cli.exe <command> ...                      (installed: next to openDNA.exe)
    .venv\\Scripts\\python.exe -m app.cli <command>  (from the repository, run from its root)

Commands: `skills list|show NAME` (the recipes in app/agent/), `inspect FILE`, `feature FILE NAME`, `region FILE START END`,
`kld check|run PLASMID MUTATIONS.csv` (primers and edited plasmids for a table of HGVS mutations).
Everything prints JSON (skills print Markdown). An error is one line `ERROR: ...` on stderr and exit code 1, never a traceback.
Positions are 1-based and inclusive, as in the GenBank file and HGVS (`region f.gb 10 20` is bases 10 to 20).
"""
import argparse
import json
import os
import sys
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")   # some modules import Qt; nothing is ever shown
sys.path.insert(0, str(Path(__file__).resolve().parent))   # the modules import each other by flat name

AGENT_DIR = Path(__file__).with_name("agent")   # the recipes; the PyInstaller spec bundles the folder next to the modules


from user_error import UserError as CliError   # a mistake the agent can fix: its message is printed as it is


def qt_setup():
    """What the modules that read settings need: a Qt application with the program's name (so the data folder is the app's)."""
    from PySide6.QtWidgets import QApplication
    import project
    app = QApplication.instance() or QApplication(sys.argv[:1])
    app.setApplicationName(project.NAME)
    return app


def load(path):
    from opendna_core import Document
    if not Path(path).is_file():
        raise CliError(f"File not found: {path}")
    try:
        return Document(str(path))
    except Exception as e:   # the core raises ValueError / OSError with its own words
        raise CliError(f"Cannot read {path} as GenBank: {e}")


def feature_dna(doc, f):
    """The bases of a feature read in its own direction (a join: the parts in order), or None for a part that wraps the origin."""
    from mutation import reverse_complement
    if any(s < 0 or e > doc.length for s, e in f.spans):
        return None
    seq = "".join(doc.sequence(s, e) for s, e in f.spans)
    return reverse_complement(seq) if f.strand < 0 else seq


def protein(doc, f):
    """{frame, table, length, sequence} of a CDS, or {"note": why none}."""
    import genetic_code
    plan = genetic_code.cds_plan(doc, f)
    dna = feature_dna(doc, f)
    if dna is None or plan.why:
        return {"note": plan.why or "crosses the origin: not translated here"}
    dna = dna[plan.frame:].upper()
    code = genetic_code.codon_map(plan.table)
    aa = "".join(code.get(dna[i:i + 3], "X") for i in range(0, len(dna) - 2, 3))
    return {"frame": plan.frame + 1, "table": plan.table, "length": len(aa.rstrip("*")), "sequence": aa}


def cmd_inspect(args):
    from mutation import hull
    doc = load(args.file)
    out = {"file": str(args.file), "name": doc.name, "length": doc.length, "circular": doc.circular, "features": []}
    for f in doc.features():
        start, end = hull(f)
        item = {"index": f.index, "kind": f.kind, "label": f.label, "start": start, "end": end, "strand": f.strand, "location": f.location}
        if f.kind == "CDS":
            item["protein"] = protein(doc, f)
        out["features"].append(item)
    return out


def cmd_feature(args):
    from mutation import find_feature, hull
    doc = load(args.file)
    f = find_feature(doc, args.name)
    start, end = hull(f)
    dna = feature_dna(doc, f)
    if dna is None:
        raise CliError(f"{args.name!r} crosses the origin of the plasmid: read it with `region` in two parts.")
    out = {"index": f.index, "kind": f.kind, "label": f.label, "start": start, "end": end, "strand": f.strand, "length": len(dna), "sequence": dna}
    if f.kind == "CDS":
        out["protein"] = protein(doc, f)
    return out


def cmd_region(args):
    doc = load(args.file)
    if not 1 <= args.start <= args.end <= doc.length:
        raise CliError(f"Positions must satisfy 1 <= START <= END <= {doc.length} (a region across the origin: ask for it in two parts).")
    return {"start": args.start, "end": args.end, "length": args.end - args.start + 1, "sequence": doc.sequence(args.start - 1, args.end)}


def given(value, default):
    """A setting for the agent to show: the value in use and where it came from."""
    return {"value": default if value is None else value, "source": "default" if value is None else "override (command line)"}


def effective_settings(extra):
    """Every setting that changes the answer, with its value and where it came from: for the agent to show before `run`. `extra`: the
    settings of the command itself."""
    import codon_usage
    import genetic_code
    import primer_panel
    from primer_db import settings
    saved = settings()
    usage = codon_usage.active()
    seed, mismatches = primer_panel.homology_settings()
    genome = primer_panel.default_genome()
    return {
        **extra,
        "pcr_conditions": {k: {"value": v, "source": "saved in openDNA" if saved.contains(f"pcr/{k}") else "default"}
                           for k, v in primer_panel.conditions().items()},
        "genetic_code_table": {"value": genetic_code.table_id(), "source": "openDNA settings"},
        "codon_usage_table": {"value": usage.name if usage else None,
                              "source": "openDNA settings" if usage else "none: ties between codons are broken alphabetically"},
        "primer_checks": {"homology_seed": {"value": seed, "source": "openDNA settings"}, "homology_mismatches": {"value": mismatches, "source": "openDNA settings"},
                          "background_genome": {"value": genome, "source": "openDNA settings (Primers > Options)" if genome else "none set: no genome check"}},
    }


def read_table(path, *columns, optional=()):
    """Rows of a CSV with the columns id and `columns` (and the `optional` ones, if present; nothing else): each row a dict."""
    import csv
    import re
    if not Path(path).is_file():
        raise CliError(f"File not found: {path}")
    with open(path, newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh)
        found, need = list(reader.fieldnames or []), ["id", *columns]
        if not set(need) <= set(found) or not set(found) <= set(need) | set(optional):
            raise CliError(f"{path}: the columns must be exactly {','.join(need)}" + (f" (and optionally {','.join(optional)})" if optional else "")
                           + f"; found: {', '.join(found) or 'none'}.")
        rows = [{k: (v or "").strip() for k, v in r.items()} for r in reader]
    if not rows:
        raise CliError(f"{path}: no rows.")
    seen, problems = set(), []
    for n, r in enumerate(rows, 2):   # (line 1 is the header)
        if not re.fullmatch(r"[A-Za-z0-9._-]+", r["id"]):
            problems.append(f"line {n}: id {r['id']!r} must be letters, digits, . _ - (it names the output file)")
        elif r["id"] in seen:
            problems.append(f"line {n}: id {r['id']!r} is used twice")
        seen.add(r["id"])
    if problems:
        raise CliError(f"{path}: " + "; ".join(problems))
    return rows


def kld_row(source, row, args, usage, out_dir):
    """One mutation: the result record, and the edited plasmid with its primers written to `out_dir` (when there is one)."""
    import kld
    import mutation
    import primer_panel
    from assembly import snapshot
    rec = {"id": row["id"], "mutation": row["mutation"], "status": "error"}
    try:
        r = mutation.resolve(source, row["mutation"], usage)
        d = kld.design_selected(source.sequence(0, None).upper(), [(r.start, r.end, r.text)], source.circular,
                                kld.DEFAULT_TM if args.tm is None else args.tm, kld.DEFAULT_GROUP if args.group is None else args.group,
                                kld.DEFAULT_WARN if args.warn is None else args.warn)
    except (CliError, ValueError) as e:   # (design raises ValueError with a sentence: too short, too many new bases, ...)
        rec["error"] = str(e)
        return rec
    rec.update(status="ok", feature=r.feature, hgvs_dna=r.hgvs_dna, hgvs_protein=r.hgvs_protein, effect=r.effect,
               warnings=r.warnings + [t for level, t in d.notes if level == "warn"], notes=[t for level, t in d.notes if level != "warn"],
               primers=[{"name": n, "sequence": s, "annealing_length": a, "annealing_tm": primer_panel.melting_temp(s[-a:])} for n, s, a in d.rows])
    if out_dir is not None:
        doc = snapshot(source, source.name)
        doc.replace(r.start, r.end, r.text)
        kld.with_primers(doc, d).write_genbank(str(out_dir / f"{row['id']}.gb"), None)
        rec["file"] = f"{row['id']}.gb"
    return rec


def cmd_kld(args):
    import codon_usage
    import kld
    source, rows = load(args.plasmid), read_table(args.mutations, "mutation")
    out_dir = None
    if args.action == "run":
        out_dir = Path(args.out) if args.out else Path(args.plasmid).with_name(Path(args.plasmid).stem + "_kld")
        if out_dir.exists() and any(out_dir.iterdir()):
            raise CliError(f"The output folder {out_dir} already exists and is not empty: choose a new one with --out (nothing is overwritten).")
        out_dir.mkdir(parents=True, exist_ok=True)
    settings = effective_settings({"kld": {"tm": given(args.tm, kld.DEFAULT_TM), "group": given(args.group, kld.DEFAULT_GROUP),
                                           "warn_insert": given(args.warn, kld.DEFAULT_WARN)}})
    usage = codon_usage.active()
    records = [kld_row(source, row, args, usage, out_dir) for row in rows]
    failed = [r for r in records if r["status"] == "error"]
    result = {"settings": settings, "summary": {"rows": len(rows), "ok": len(rows) - len(failed), "errors": len(failed),
                                                 "with_warnings": sum(bool(r.get("warnings")) for r in records)}, "rows": records}
    if out_dir is not None:
        import csv
        with open(out_dir / "primers.csv", "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow(["id", "mutation", "hgvs_dna", "hgvs_protein", "forward_primer", "forward_tm", "reverse_primer", "reverse_tm", "warnings"])
            for r in records:
                if r["status"] == "ok":
                    f, rv = r["primers"]
                    tm = lambda p: round(p["annealing_tm"], 1) if p["annealing_tm"] else ""
                    w.writerow([r["id"], r["mutation"], r["hgvs_dna"], r["hgvs_protein"], f["sequence"], tm(f), rv["sequence"], tm(rv),
                                " | ".join(r["warnings"])])
        (out_dir / "report.jsonl").write_text("".join(json.dumps(r) + "\n" for r in records), encoding="utf-8")
        (out_dir / "settings.json").write_text(json.dumps(settings, indent=1), encoding="utf-8")
        result["out"] = str(out_dir)
    return result, 1 if failed else 0


def target_range(doc, text):
    """(first, last), 1-based inclusive, of a feature (name or #INDEX) or of START..END."""
    import re
    from mutation import find_feature, hull
    m = re.fullmatch(r"(\d+)\.\.(\d+)", text)
    if not m:
        return hull(find_feature(doc, text))
    first, last = int(m.group(1)), int(m.group(2))
    if not 1 <= first < last <= doc.length:
        raise CliError(f"{text!r}: positions must satisfy 1 <= START < END <= {doc.length}.")
    return first, last


def pcr_row(doc, row, args):
    """The best primer pairs for one target (a feature name or START..END, 1-based inclusive): the record of the row."""
    from types import SimpleNamespace
    import primer_design
    import primer_panel
    rec = {"id": row["id"], "target": row["target"], "status": "error"}
    try:
        first, last = target_range(doc, row["target"])
        opts = {k: v for k, v in dict(tm_target=args.tm, min_len=args.min_len, max_len=args.max_len, window=args.window).items() if v is not None}
        pairs = primer_design.design_region(doc, first - 1, last, **opts)
    except CliError as e:
        rec["error"] = str(e)
        return rec
    if not pairs:
        rec["error"] = "No pair found: widen the search window (--window), the length range (--min-len, --max-len) or the Tm (--tm)."
        return rec
    ed = SimpleNamespace(length=doc.length, doc=doc)   # (what primer_checks reads of an editor)
    out = []
    for rank, p in enumerate(pairs, 1):
        checks = {k: primer_panel.primer_checks(ed, p[k]["seq"], "wait") for k in ("fwd", "rev")}
        out.append({"rank": rank, "product_size": p["size"], "penalty": round(p["penalty"], 2), "notes": p["notes"],
                    "mark": primer_panel.worst(checks["fwd"] + checks["rev"]),
                    **{k: {"sequence": p[k]["seq"], "tm": round(p[k]["tm"], 1), "gc": round(p[k]["gc"]), "start": p[k]["start"] + 1, "end": p[k]["end"],
                           "checks": [{"check": n, "level": lv, "text": t} for n, lv, t in checks[k]]} for k in ("fwd", "rev")}})
    rec.update(status="ok", start=first, end=last, pairs=out)
    return rec


def cmd_pcr(args):
    import primer_design
    doc, rows = load(args.plasmid), read_table(args.targets, "target")
    out_dir = None
    if args.action == "run":
        out_dir = Path(args.out) if args.out else Path(args.plasmid).with_name(Path(args.plasmid).stem + "_pcr")
        if out_dir.exists() and any(out_dir.iterdir()):
            raise CliError(f"The output folder {out_dir} already exists and is not empty: choose a new one with --out (nothing is overwritten).")
        out_dir.mkdir(parents=True, exist_ok=True)
    settings = effective_settings({"pcr_design": {"tm": given(args.tm, primer_design.DEFAULT_TM), "min_len": given(args.min_len, primer_design.DEFAULT_MIN_LEN),
                                                  "max_len": given(args.max_len, primer_design.DEFAULT_MAX_LEN),
                                                  "window": given(args.window, primer_design.DEFAULT_WINDOW)}})
    records = [pcr_row(doc, row, args) for row in rows]
    failed = [r for r in records if r["status"] == "error"]
    result = {"settings": settings, "summary": {"rows": len(rows), "ok": len(rows) - len(failed), "errors": len(failed),
                                                 "with_warnings": sum(any(p["mark"] in ("warn", "bad") for p in r.get("pairs", [])) for r in records)},
              "rows": records}
    if out_dir is not None:
        import csv
        with open(out_dir / "primers.csv", "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow(["id", "target", "rank", "forward_primer", "forward_tm", "reverse_primer", "reverse_tm", "product_size", "penalty", "mark", "notes"])
            for r in records:
                for p in r.get("pairs", []):
                    w.writerow([r["id"], r["target"], p["rank"], p["fwd"]["sequence"], p["fwd"]["tm"], p["rev"]["sequence"], p["rev"]["tm"],
                                p["product_size"], p["penalty"], p["mark"], " | ".join(p["notes"])])
        (out_dir / "report.jsonl").write_text("".join(json.dumps(r) + "\n" for r in records), encoding="utf-8")
        (out_dir / "settings.json").write_text(json.dumps(settings, indent=1), encoding="utf-8")
        result["out"] = str(out_dir)
    return result, 1 if failed else 0


def rec_row(genome, row, args, folder, cassettes, out_dir):
    """One recombineering design: the cassette (a GenBank file, or `none` for a deletion) replaces a region of the genome."""
    import primer_design
    import primer_panel
    import recombineering as R
    from assembly import transfer_primers
    from types import SimpleNamespace
    rec = {"id": row["id"], "region": row["region"], "cassette": row["cassette"], "status": "error"}
    try:
        first, last = target_range(genome, row["region"])
        orientation = row.get("orientation", "").lower() or "forward"
        if orientation not in ("forward", "reverse"):
            raise CliError(f"orientation must be forward or reverse, not {orientation!r}.")
        if row["cassette"].lower() in ("none", "-"):   # a deletion
            cassette = Document_new("none", "")
        else:
            path = Path(row["cassette"])
            path = path if path.is_absolute() else folder / path
            if str(path) not in cassettes:
                cassettes[str(path)] = load(path)
            cassette = cassettes[str(path)]
        left = R.HOMOLOGY if args.homology is None else args.homology
        d = R.design(genome, cassette, first - 1, last, left, left, orientation == "reverse",
                     primer_design.DEFAULT_TM if args.tm is None else args.tm)
    except (CliError, ValueError) as e:
        rec["error"] = str(e)
        return rec
    warnings = []
    for f in genome.features():   # (by the parts of a feature: one that joins the two ends of the plasmid is not "all of it")
        if f.kind == "source" or not any(s < last and first - 1 < e for s, e in f.spans):
            continue
        if all(first - 1 <= s and e <= last for s, e in f.spans):
            warnings.append(f"The region removes the {f.kind} {f.label!r} ({f.location}).")
        else:
            warnings.append(f"The region cuts the {f.kind} {f.label!r} ({f.location}): only part of it is replaced.")
    rec.update(status="ok", start=first, end=last, replaced_bp=last - first + 1, cassette_bp=sum(len(s) for s in [cassette.sequence(0, None)]),
               orientation=orientation, left_homology=d.left, right_homology=d.right, warnings=warnings,
               primers=[{"name": n, "sequence": s, "annealing_length": a, "annealing_tm": primer_panel.melting_temp(s[-a:])} for n, s, a, _ in d.rows])
    if out_dir is not None:
        if d.rows:
            transfer_primers(d.product, [SimpleNamespace(fwd=d.rows[0][:2], rev=d.rows[1][:2])])
        if not args.primers_only:
            d.product.write_genbank(str(out_dir / f"{row['id']}.gb"), None)
            rec["file"] = f"{row['id']}.gb"
    return rec


def Document_new(name, bases):
    from opendna_core import Document
    return Document.new_sequence(name, bases, False)


def cmd_recombineering(args):
    import recombineering as R
    genome = load(args.genome)
    rows = read_table(args.table, "region", "cassette", optional=("orientation",))
    out_dir = None
    if args.action == "run":
        out_dir = Path(args.out) if args.out else Path(args.genome).with_name(Path(args.genome).stem + "_rec")
        if out_dir.exists() and any(out_dir.iterdir()):
            raise CliError(f"The output folder {out_dir} already exists and is not empty: choose a new one with --out (nothing is overwritten).")
        out_dir.mkdir(parents=True, exist_ok=True)
    import primer_design
    settings = effective_settings({"recombineering": {"homology": given(args.homology, R.HOMOLOGY), "annealing_tm": given(args.tm, primer_design.DEFAULT_TM),
                                                      "write_genomes": {"value": not args.primers_only, "source": "default" if not args.primers_only else "override (command line)"}}})
    cassettes = {}
    records = [rec_row(genome, row, args, Path(args.table).resolve().parent, cassettes, out_dir) for row in rows]
    failed = [r for r in records if r["status"] == "error"]
    result = {"settings": settings, "summary": {"rows": len(rows), "ok": len(rows) - len(failed), "errors": len(failed),
                                                 "with_warnings": sum(bool(r.get("warnings")) for r in records)}, "rows": records}
    if out_dir is not None:
        import csv
        with open(out_dir / "primers.csv", "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow(["id", "region", "cassette", "orientation", "forward_primer", "forward_tm", "reverse_primer", "reverse_tm", "warnings"])
            for r in records:
                if r["status"] == "ok":
                    f, rv = (r["primers"] + [{"sequence": "", "annealing_tm": None}] * 2)[:2]   # (a deletion has no primers: nothing to amplify)
                    tm = lambda p: round(p["annealing_tm"], 1) if p["annealing_tm"] else ""
                    w.writerow([r["id"], r["region"], r["cassette"], r["orientation"], f["sequence"], tm(f), rv["sequence"], tm(rv), " | ".join(r["warnings"])])
        (out_dir / "report.jsonl").write_text("".join(json.dumps(r) + "\n" for r in records), encoding="utf-8")
        (out_dir / "settings.json").write_text(json.dumps(settings, indent=1), encoding="utf-8")
        result["out"] = str(out_dir)
    return result, 1 if failed else 0


def recipes():
    """{name: path} of the recipes (the index is `index`)."""
    found = {p.stem: p for p in sorted((AGENT_DIR / "skills").glob("*.md"))}
    return {"index": AGENT_DIR / "index.md", **found}


def cmd_skills(args):
    names = recipes()
    if args.action == "list":
        return "\n".join(f"{n}: " + next((ln.lstrip("# ").strip() for ln in p.read_text(encoding="utf-8").splitlines() if ln.strip()), "") for n, p in names.items())
    if args.name not in names:
        raise CliError(f"No recipe called {args.name!r}. Choices: {', '.join(names)}.")
    return names[args.name].read_text(encoding="utf-8")


def parser():
    p = argparse.ArgumentParser(prog="oDNA-cli", description=__doc__.split("\n\n")[0], formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="command", required=True)
    s = sub.add_parser("skills", help="the recipes for agents: `skills list`, `skills show index` first")
    s.add_argument("action", choices=["list", "show"])
    s.add_argument("name", nargs="?", default="index")
    s.set_defaults(run=cmd_skills, text=True)
    s = sub.add_parser("inspect", help="features of a GenBank file with positions, strand and the protein of each CDS")
    s.add_argument("file")
    s.set_defaults(run=cmd_inspect)
    s = sub.add_parser("feature", help="the sequence (and protein) of one feature, by label or #index")
    s.add_argument("file")
    s.add_argument("name")
    s.set_defaults(run=cmd_feature)
    s = sub.add_parser("kld", help="KLD primers and edited plasmids for a table of HGVS mutations (check: report only; run: also write files)")
    s.add_argument("action", choices=["check", "run"])
    s.add_argument("plasmid")
    s.add_argument("mutations", help="CSV with the columns id,mutation")
    s.add_argument("--out", help="run: a new output folder (default: PLASMID_kld next to the plasmid)")
    s.add_argument("--tm", type=float, help="annealing Tm in deg C (default 60)")
    s.add_argument("--group", type=int, help="changes closer than this many bp share one primer pair (default 30)")
    s.add_argument("--warn", type=int, help="warn for insertions above this many bp (default 120)")
    s.set_defaults(run=cmd_kld, returns_code=True)
    s = sub.add_parser("pcr", help="PCR primer pairs for a table of targets (a feature name or START..END), with the Primers panel's checks")
    s.add_argument("action", choices=["check", "run"])
    s.add_argument("plasmid")
    s.add_argument("targets", help="CSV with the columns id,target")
    s.add_argument("--out", help="run: a new output folder (default: PLASMID_pcr next to the plasmid)")
    s.add_argument("--tm", type=float, help="target Tm in deg C (default 60)")
    s.add_argument("--min-len", type=int, help="shortest primer (default 18)")
    s.add_argument("--max-len", type=int, help="longest primer (default 25)")
    s.add_argument("--window", type=int, help="how far from the target ends a primer may start, in bp (default 0: on the ends)")
    s.set_defaults(run=cmd_pcr, returns_code=True)
    s = sub.add_parser("recombineering", help="lambda Red: primers (homology tails) and the recombinant genome for cassettes that replace regions")
    s.add_argument("action", choices=["check", "run"])
    s.add_argument("genome")
    s.add_argument("table", help="CSV with the columns id,region,cassette (and optionally orientation)")
    s.add_argument("--out", help="run: a new output folder (default: GENOME_rec next to the genome)")
    s.add_argument("--homology", type=int, help="bp of homology on each side of the cassette (default 50)")
    s.add_argument("--tm", type=float, help="annealing Tm in deg C (default 60)")
    s.add_argument("--primers-only", action="store_true", help="run: do not write the recombinant genomes (a bacterial genome is megabytes each)")
    s.set_defaults(run=cmd_recombineering, returns_code=True)
    s = sub.add_parser("region", help="the bases START..END (1-based, inclusive)")
    s.add_argument("file")
    s.add_argument("start", type=int)
    s.add_argument("end", type=int)
    s.set_defaults(run=cmd_region)
    return p


def main(argv=None):
    args = parser().parse_args(argv)
    for stream in (sys.stdout, sys.stderr):   # a Windows console is not UTF-8: primer notes and recipes have symbols
        if hasattr(stream, "reconfigure"):
            stream.reconfigure(encoding="utf-8")
    try:
        if args.command != "skills":
            qt_setup()   # (genetic_code reads the settings file: the same table choice as the windows)
        result = args.run(args)
    except CliError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1
    code = 0
    if getattr(args, "returns_code", False):   # a batch: failed rows make the exit code 1, the rest is still reported
        result, code = result
    print(result if getattr(args, "text", False) else json.dumps(result, indent=1))
    return code


if __name__ == "__main__":
    sys.exit(main())
