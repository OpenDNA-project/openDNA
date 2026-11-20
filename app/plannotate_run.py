"""The half of the pLannotate link that runs inside its own conda environment (plannotate_tool.py starts it there):
`python plannotate_run.py IN.json OUT.json`.

It calls pLannotate's own function, `plannotate.annotate.annotate`, and hands back the rows it returns, with positions as
openDNA draws them. Nothing of pLannotate is copied or changed here. It imports nothing of openDNA and nothing outside the
standard library (but pLannotate, inside the function that runs it), so the environment needs no other package.
"""
import json
import os
import sys


def rows_from(records, length):
    """openDNA's view of pLannotate's rows (dicts with its column names): [{name, kind, start, end, strand, identity, match,
    database, fragment, note}]. start 0-based, end exclusive and past `length` for a feature across the origin (pLannotate gives
    qstart > qend there), as openDNA's own hits do; strand 1 or -1; identity and match in 0..1."""
    out = []
    for r in records:
        start, end = int(r["qstart"]), int(r["qend"])
        if end <= start:
            end += length
        kind = str(r.get("type") or "misc_feature")
        out.append({"name": str(r["name"]), "kind": "rep_origin" if kind == "origin of replication" else kind, "start": start, "end": end,
                    "strand": -1 if r.get("sframe") == -1 else 1, "identity": round(float(r["pident"]) / 100, 4),
                    "match": round(float(r["abs percmatch"]) / 100, 4), "database": str(r["db"]), "fragment": bool(r["fragment"]),
                    "note": str(r.get("blurb") or "")})
    out.sort(key=lambda x: (x["start"], x["end"]))
    return out


def run(sequence, circular, cores):
    import importlib.metadata
    from plannotate.annotate import annotate
    # is_detailed: the type of each hit (CDS, promoter, terminator...) instead of one anonymous class; cores: every core of the
    # machine for the searches (pLannotate's default is 1; the results do not depend on it)
    table = annotate(sequence, linear=not circular, is_detailed=True, cores=cores)
    return {"version": importlib.metadata.version("plannotate"), "rows": rows_from(table.to_dict("records"), len(sequence))}


def status():
    """One line for the caller: `ready VERSION`, or `nodb VERSION` when `plannotate setupdb` has not been run."""
    import importlib.metadata
    from plannotate import _package_data
    return ("ready " if _package_data.databases_exist() else "nodb ") + importlib.metadata.version("plannotate")


def main(argv):
    if argv[1] == "--status":
        print(status())
        return
    with open(argv[1], encoding="utf-8") as fh:
        job = json.load(fh)
    try:
        answer = run(job["sequence"], job["circular"], job.get("cores") or os.cpu_count() or 1)
    except Exception as e:   # handed back as text: the caller shows it
        answer = {"error": f"{type(e).__name__}: {e}"}
    with open(argv[2], "w", encoding="utf-8") as fh:
        json.dump(answer, fh)


if __name__ == "__main__":
    main(sys.argv)
