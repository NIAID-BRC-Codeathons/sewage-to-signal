"""Stage 04 - dereplication.

Wastewater surveillance resequences the same sewershed over and over, so the
protein set across samples is heavily redundant. Collapsing it before the GPU
stage is the cheapest large saving in the whole pipeline.

This stage does not remove anything. It *labels*: every gene gets the id of its
cluster representative and a flag saying whether it is that representative.
Dereplication then stops being a position in the pipeline and becomes a
predicate - ``is_representative`` - that any later stage can opt into or
ignore. That matters because the right answer differs by question: the GPU
stage wants one protein per cluster, while counting how much of a sample a
family covers wants all of them. The old ``nr.faa`` silently made that choice
once, for everybody.

Exact dedup (by sequence hash) is pure Python and always available. Clustering
at sub-100% identity needs MMseqs2; when it is missing the stage still runs and
reports how much exact dedup alone achieved.
"""

from __future__ import annotations

import argparse
import hashlib
import subprocess
import time
from pathlib import Path

import lake
from common import StageResult, which, workdir, write_fasta
from entities import fasta_records, table_from_fasta, write_rows
from stage import Column, Param, Stage, Tool

COLUMN_HELP = {
    "rep_id": "gene_id of this gene's cluster representative",
    "is_representative": "true when the gene represents its own cluster",
    "cluster_size": "how many genes collapsed into this one's cluster",
}


def run(
    rows,
    out_dir: Path,
    sample: str,
    con=None,
    source: Path | None = None,
    where: str | None = None,
    identity: float = 0.95,
    coverage: float = 0.8,
    threads: int = 4,
    use_mmseqs: bool = True,
    force: bool = False,
) -> StageResult:
    out_dir = Path(out_dir)
    engine = "mmseqs" if (use_mmseqs and which("mmseqs")) else "exact"
    params = {"identity": identity, "coverage": coverage, "engine": engine,
              "where": where}
    deps = [lake.fingerprint(source)] if source else []
    if not force and lake.is_current(con, sample, "s04_derep", params, deps):
        return StageResult("s04_derep", out_dir, {"engine": engine}, skipped=True,
                           produced={"gene": None})

    import pyarrow as pa

    t0 = time.time()
    records = fasta_records(rows)
    # Exact dedup first: it is free and shrinks the input to any clusterer.
    first: dict[str, str] = {}          # sha1 -> first gene_id with it
    rep_of: dict[str, str] = {}
    for gid, seq in records:
        h = hashlib.sha1(seq.encode()).hexdigest()
        rep_of[gid] = first.setdefault(h, gid)
    n_in = len(records)
    exact_reps = [(gid, seq) for gid, seq in records if rep_of[gid] == gid]
    stats = {"genes_in": n_in, "after_exact": len(exact_reps)}

    if engine == "mmseqs":
        tmp = out_dir / "_mm"
        tmp.mkdir(exist_ok=True)
        exact_fa = tmp / "exact.faa"
        write_fasta(exact_fa, exact_reps)
        pref = tmp / "clu"
        subprocess.run(
            ["mmseqs", "easy-linclust", str(exact_fa), str(pref), str(tmp / "tmp"),
             "--min-seq-id", str(identity), "-c", str(coverage),
             "--threads", str(threads)],
            check=True, capture_output=True,
        )
        clu_tsv = Path(str(pref) + "_cluster.tsv")
        if clu_tsv.exists():
            # mmseqs clusters the exact representatives; fold its mapping back
            # onto every original gene so the column is complete.
            second = {}
            for line in clu_tsv.read_text().splitlines():
                rep, mem = line.split("\t")[:2]
                second[mem] = rep
            rep_of = {gid: second.get(r, r) for gid, r in rep_of.items()}
        stats["after_cluster"] = len(set(rep_of.values()))

    sizes: dict[str, int] = {}
    for r in rep_of.values():
        sizes[r] = sizes.get(r, 0) + 1
    out_rows = [{"gene_id": gid, "rep_id": rep_of[gid],
                 "is_representative": rep_of[gid] == gid,
                 "cluster_size": sizes[rep_of[gid]]}
                for gid, _ in records]
    table = pa.Table.from_pylist(out_rows, schema=pa.schema([
        ("gene_id", pa.string()), ("rep_id", pa.string()),
        ("is_representative", pa.bool_()), ("cluster_size", pa.int32()),
    ]))
    write_rows(con, "gene", sample, table, STAGE)

    n_reps = sum(1 for r in out_rows if r["is_representative"])
    stats["representatives"] = n_reps
    stats["reduction"] = round(1 - n_reps / n_in, 4) if n_in else 0.0
    stats["engine"] = engine
    el = time.time() - t0
    lake.record_run(con, sample, "s04_derep", params, deps, where, stats,
                    seconds=el, rows=len(out_rows))
    return StageResult("s04_derep", out_dir, stats, seconds=el,
                       produced={"gene": None})


STAGE = Stage(
    name="s04_derep",
    title="Dereplicate",
    summary="Label each gene with its cluster representative. Removes nothing; "
            "downstream stages select on is_representative when they want one "
            "protein per cluster.",
    run=run,
    consumes="gene",
    produces="gene",
    order=40,
    selectable=True,
    roles=("dereplication",),
    adds=(
        Column("rep_id", "string", COLUMN_HELP["rep_id"]),
        Column("is_representative", "bool", COLUMN_HELP["is_representative"]),
        Column("cluster_size", "int32", COLUMN_HELP["cluster_size"]),
    ),
    params=(
        Param("identity", float, 0.95, group="derep",
              help="MMseqs2 clustering identity; ignored without mmseqs"),
        Param("coverage", float, 0.8, group="derep",
              help="MMseqs2 alignment coverage"),
        Param("use_mmseqs", bool, True, group="derep",
              help="cluster below 100% identity when mmseqs is on PATH"),
    ),
    requires=(Tool("mmseqs", optional=True,
                   hint="without it only exact duplicates collapse"),),
)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("proteins", type=Path)
    p.add_argument("--sample", required=True)
    p.add_argument("--work", type=Path, default=Path("work"))
    p.add_argument("--identity", type=float, default=0.95)
    p.add_argument("--no-mmseqs", action="store_true")
    p.add_argument("--force", action="store_true")
    a = p.parse_args()
    r = run(table_from_fasta(a.proteins), workdir(a.work, a.sample, "s04_derep"), a.sample,
            source=a.proteins, identity=a.identity,
            use_mmseqs=not a.no_mmseqs, force=a.force)
    print(r.describe())


if __name__ == "__main__":
    main()
