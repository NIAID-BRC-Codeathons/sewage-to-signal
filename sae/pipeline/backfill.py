"""Load work directories written before the store existed into the lake.

There are two older shapes on disk and this reads both.

**Fragment runs** wrote one parquet per stage plus a ``<output>.manifest.json``
naming the level and role. Those are read from the manifest, which is exact.

**Pre-fragment runs** wrote FASTAs and TSVs - ``proteins.faa``, ``genes.tsv``,
``classification.tsv``, ``nr.faa`` - and are reconstructed by reading them. Less
exact, because those files record less, but it is the difference between a
cohort that starts with eleven samples and one that starts empty.

Idempotent: a sample is replaced wholesale, so re-running is safe.

    python backfill.py --work ../../work --lake ../../data/sae.ducklake
    python backfill.py --work ../../work --dry-run
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
import time
from pathlib import Path

import lake
from entities import LEVELS, gene_row, gene_schema
from stage import registry

# Directories that sit in a work root without being samples.
NOT_SAMPLES = {"uploads", "logs", ".logs"}


def manifests(sample_dir: Path) -> list[tuple[Path, dict]]:
    out = []
    for stage_dir in sorted(p for p in sample_dir.iterdir() if p.is_dir()):
        for mf in sorted(stage_dir.glob("*.manifest.json")):
            try:
                out.append((stage_dir, json.loads(mf.read_text())))
            except (json.JSONDecodeError, OSError):
                continue
    return out


def resolve(stage_dir: Path, recorded: str) -> Path | None:
    """A fragment's parquet, allowing for paths recorded inside a container.

    A run executed in the image records ``/work/META1/s03_genes/...`` because
    run.sh binds the work root to /work. On the host that path does not exist,
    so fall back to the same filename in the directory the manifest came from.
    """
    p = Path(recorded)
    if p.is_file():
        return p
    alt = stage_dir / p.name
    return alt if alt.is_file() else None


# --------------------------------------------------------------------------
# pre-fragment readers
# --------------------------------------------------------------------------
def legacy_genes(stage_dir: Path):
    """gene rows from proteins.faa plus the coordinate TSV beside it."""
    from common import read_fasta

    faa = next(iter(sorted(stage_dir.glob("*.proteins.faa"))), None)
    if not faa:
        return []
    coords = {}
    tsv = next(iter(sorted(stage_dir.glob("*.genes.tsv"))), None)
    if tsv:
        with open(tsv) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                coords[r["gene_id"]] = r
    rows = []
    for header, seq in read_fasta(faa):
        gid = header.split()[0]
        c = coords.get(gid, {})
        rows.append(gene_row(
            gid, seq.rstrip("*"), contig=c.get("contig"),
            begin=int(c["begin"]) if c.get("begin") else None,
            end=int(c["end"]) if c.get("end") else None,
            strand=c.get("strand"),
            partial=bool(int(c["partial"])) if c.get("partial") else False))
    return rows


def legacy_classification(stage_dir: Path) -> list[dict]:
    tsv = next(iter(sorted(stage_dir.glob("*.classification.tsv"))), None)
    if not tsv:
        return []
    out = []
    with open(tsv) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            out.append({
                "gene_id": r["gene_id"], "category": r.get("category") or None,
                "family": r.get("family") or None,
                "family_acc": r.get("family_acc") or None,
                "evalue": float(r["evalue"]) if r.get("evalue") else None,
                "coverage": float(r["coverage"]) if r.get("coverage") else None,
                "n_domains": int(r["n_domains"]) if r.get("n_domains") else None,
            })
    return out


def legacy_derep(stage_dir: Path) -> list[dict]:
    """rep_id per gene, from the old representative->member mapping."""
    tsv = next(iter(sorted(stage_dir.glob("*.derep.tsv"))), None)
    if not tsv:
        return []
    out, seen = [], set()
    with open(tsv) as fh:
        first = fh.readline()
        if not first.lower().startswith("representative"):
            fh.seek(0)
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            rep, mem = parts[0], parts[1]
            if mem in seen:
                continue
            seen.add(mem)
            out.append({"gene_id": mem, "rep_id": rep,
                        "is_representative": rep == mem})
    return out


# --------------------------------------------------------------------------
# loading one sample
# --------------------------------------------------------------------------
def load_sample(con, sample: str, sample_dir: Path, reg, verbose=False) -> dict:
    import pyarrow as pa
    import pyarrow.parquet as pq

    counts: dict[str, int] = {}
    docs = manifests(sample_dir)
    fragment_runs = [(d, doc) for d, doc in docs if doc.get("tables")]

    def note(what, n):
        counts[what] = counts.get(what, 0) + n
        if verbose:
            print(f"      {what}: {n}")

    if fragment_runs:
        # Base fragments first, so annotations have rows to land on.
        ordered = sorted(
            ((sd, doc, t) for sd, doc in fragment_runs for t in doc["tables"]),
            key=lambda x: (x[2].get("role") != "base", x[1].get("written") or ""))
        for stage_dir, doc, t in ordered:
            path = resolve(stage_dir, t.get("path", ""))
            if path is None:
                continue
            table = pq.read_table(path)
            level, role = t["level"], t.get("role", "annotation")
            _write(con, level, sample, table, role)
            note(f"{level} ({role})", table.num_rows)
    else:
        rows = []
        for stage_dir, _ in docs:
            rows = rows or legacy_genes(stage_dir)
        if not rows:
            for stage_dir in sorted(p for p in sample_dir.iterdir() if p.is_dir()):
                rows = rows or legacy_genes(stage_dir)
        if rows:
            _write(con, "gene", sample,
                   pa.Table.from_pylist(rows, schema=gene_schema()), "base")
            note("gene (legacy faa)", len(rows))
        for stage_dir in sorted(p for p in sample_dir.iterdir() if p.is_dir()):
            for reader, cols in ((legacy_derep, ["gene_id", "rep_id", "is_representative"]),
                                 (legacy_classification, None)):
                got = reader(stage_dir)
                if got:
                    _write(con, "gene", sample, pa.Table.from_pylist(got),
                           "annotation")
                    note(f"gene (legacy {stage_dir.name})", len(got))
            for fp in sorted(stage_dir.glob("*.sae_features.parquet")):
                t = pq.read_table(fp)
                _write(con, "feature_hit", sample, t, "base")
                note("feature_hit (legacy)", t.num_rows)

    # Provenance, preserved verbatim: params and the recorded predicate are what
    # is_current compares, so rewriting them here would mark every stage stale.
    for stage_dir, doc in docs:
        stage = doc.get("stage") or stage_dir.name
        tables = doc.get("tables") or []
        lake.record_run(
            con, sample, stage, doc.get("params") or {}, doc.get("inputs") or [],
            next((t.get("where") for t in tables if t.get("where")), None),
            doc.get("stats") or {}, doc.get("tools") or {},
            seconds=doc.get("seconds"))
    lake.register_sample(con, sample, "backfill", str(sample_dir))
    return counts


def _write(con, level: str, sample: str, table, role: str) -> None:
    """Insert or update without a Stage to ask - role comes from the manifest."""
    from entities import _annotate, _insert

    lv = LEVELS[level]
    con.register("_incoming", table)
    try:
        if role == "base":
            _insert(con, lv, sample, table)
        else:
            cols = [c for c in table.column_names
                    if c not in lv.key and c != "sample"]
            _annotate(con, lv, sample, cols)
    finally:
        con.unregister("_incoming")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--work", type=Path, action="append", required=True,
                   help="work root to read; repeatable")
    p.add_argument("--lake", default=None, help="target store")
    p.add_argument("--sample", action="append",
                   help="load only these samples; repeatable")
    p.add_argument("--dry-run", action="store_true")
    p.add_argument("--rename", action="append", default=[], metavar="OLD=NEW",
                   help="resolve a name collision between two work roots")
    p.add_argument("--verbose", action="store_true")
    a = p.parse_args()

    repo = Path(__file__).resolve().parent.parent.parent
    target = a.lake or lake.default_target(repo)
    renames = dict(r.split("=", 1) for r in a.rename)
    reg = registry()

    # Gather first, so a collision is reported before anything is written.
    found: dict[str, Path] = {}
    clashes: list[str] = []
    for root in a.work:
        root = root.resolve()
        if not root.is_dir():
            sys.exit(f"no such work root: {root}")
        for d in sorted(x for x in root.iterdir() if x.is_dir()):
            if d.name.startswith(".") or d.name in NOT_SAMPLES:
                continue
            if not any(x.is_dir() for x in d.iterdir()):
                continue
            name = renames.get(d.name, d.name)
            if a.sample and name not in a.sample:
                continue
            if name in found:
                clashes.append(f"{name}: {found[name]} and {d}")
                continue
            found[name] = d

    if clashes:
        print("Name collisions across work roots - `sample` is the key, so two "
              "different samples cannot share one:", file=sys.stderr)
        for c in clashes:
            print(f"  {c}", file=sys.stderr)
        sys.exit("Pass --rename OLD=NEW to disambiguate.")

    print(f"lake: {target}")
    print(f"{len(found)} sample(s) to load\n")
    if a.dry_run:
        for name, d in found.items():
            stages = [x.name for x in sorted(d.iterdir()) if x.is_dir()]
            print(f"  {name:<24} {', '.join(stages)}")
        return

    t0 = time.time()
    for name, d in found.items():
        with lake.open(target) as con:
            lake.ensure_schema(con, reg, LEVELS)
            try:
                counts = load_sample(con, name, d, reg, verbose=a.verbose)
            except Exception as exc:
                print(f"  {name:<24} FAILED: {type(exc).__name__}: {exc}")
                continue
        if counts:
            summary = ", ".join(f"{k}={v}" for k, v in counts.items())
        else:
            # Not a failure: a sample whose stages all produce files (QC,
            # assembly) has nothing tabular yet, and so does one whose stage
            # directory exists because a run was interrupted before it wrote.
            stages = [x.name for x in sorted(d.iterdir()) if x.is_dir()]
            summary = f"no rows yet ({', '.join(stages) or 'empty'})"
        print(f"  {name:<24} {summary}")

    with lake.read(target) as con:
        print(f"\nloaded in {time.time() - t0:.1f}s")
        for level in LEVELS:
            try:
                n = con.execute(f"SELECT count(*) FROM {level}").fetchone()[0]
            except Exception:
                n = 0
            if n:
                print(f"  {level:<14}{n:>9} rows")


if __name__ == "__main__":
    main()
