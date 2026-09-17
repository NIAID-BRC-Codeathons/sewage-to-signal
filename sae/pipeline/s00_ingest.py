"""Stage 00 - import proteins somebody else called.

Entering the pipeline with a protein FASTA used to be a special case in the
driver: a flag that quietly set the start stage and skipped the first three.
That made "where does this pipeline begin" a property of ``run.py`` rather than
of the stages, and it meant the web UI had to know the same rule.

As a stage it is just another edge in the graph - ``proteins`` (a file) to
``gene`` (a level) - so the planner discovers it, and a caller that has
proteins and wants features gets this stage at the front of the plan without
anybody encoding that.

It is the only other writer of the gene level's base rows, so it takes the
schema from ``entities`` rather than defining its own. Coordinates are null:
an imported protein has no contig to have come from.
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path

import lake
from common import StageResult, read_fasta, workdir
from entities import gene_columns, gene_row, gene_schema, write_rows
from stage import Param, Stage


def run(
    proteins: Path,
    out_dir: Path,
    sample: str,
    con=None,
    min_aa: int = 0,
    force: bool = False,
) -> StageResult:
    proteins = Path(proteins)
    params = {"min_aa": min_aa, "source": "imported"}
    deps = [lake.fingerprint(proteins)]
    if not force and lake.is_current(con, sample, "s00_ingest", params, deps):
        return StageResult("s00_ingest", proteins, {}, skipped=True,
                           produced={"gene": None})

    import pyarrow as pa

    t0 = time.time()
    rows, seen, n_in, n_short, n_dup = [], set(), 0, 0, 0
    for header, seq in read_fasta(proteins):
        n_in += 1
        seq = seq.rstrip("*")
        if min_aa and len(seq) < min_aa:
            n_short += 1
            continue
        gid = header.split()[0]
        if gid in seen:
            # Duplicate ids would break every join downstream, so they are a
            # hard error rather than a silently-kept last-one-wins.
            n_dup += 1
            raise SystemExit(
                f"{proteins}: duplicate sequence id {gid!r}. Gene ids have to "
                f"be unique within a sample - they are the level's key.")
        seen.add(gid)
        rows.append(gene_row(gid, seq))

    if not rows:
        raise SystemExit(f"no sequences kept from {proteins}")

    table = pa.Table.from_pylist(rows, schema=gene_schema())
    n = write_rows(con, "gene", sample, table, STAGE)

    stats = {"proteins_in": n_in, "genes": len(rows), "too_short": n_short}
    el = time.time() - t0
    lake.record_run(con, sample, "s00_ingest", params, deps, None, stats,
                    seconds=el, rows=n)
    return StageResult("s00_ingest", proteins, stats, seconds=el,
                       produced={"gene": None})


STAGE = Stage(
    name="s00_ingest",
    title="Import proteins",
    summary="Load an externally called protein FASTA as the gene level, so a "
            "run can start from proteins without the driver special-casing it.",
    run=run,
    consumes="proteins",
    produces="gene",
    order=0,
    roles=("gene_source",),
    adds=gene_columns(),
    params=(
        Param("min_aa", int, 0, group="genes",
              help="drop imported sequences shorter than this (0 keeps all)"),
    ),
    requires=(),
)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("proteins", type=Path)
    p.add_argument("--sample", required=True)
    p.add_argument("--work", type=Path, default=Path("work"))
    p.add_argument("--min-aa", type=int, default=0)
    p.add_argument("--force", action="store_true")
    a = p.parse_args()
    r = run(a.proteins, workdir(a.work, a.sample, "s00_ingest"), a.sample,
            min_aa=a.min_aa, force=a.force)
    print(r.describe())


if __name__ == "__main__":
    main()
