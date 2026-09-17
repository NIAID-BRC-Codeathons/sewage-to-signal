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

from common import (StageResult, is_current, read_fasta, workdir,
                    write_fragment, write_manifest)
from entities import GENE_COLUMN_HELP, gene_row, gene_schema
from stage import Column, Param, Stage


def run(
    proteins: Path,
    out_dir: Path,
    sample: str,
    min_aa: int = 0,
    force: bool = False,
) -> StageResult:
    proteins = Path(proteins)
    out = Path(out_dir) / f"{sample}.genes.parquet"
    params = {"min_aa": min_aa, "source": "imported"}
    if not force and is_current(out, [proteins], params):
        return StageResult("s00_ingest", out, {}, skipped=True,
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
    frag = write_fragment(out, table, "gene", role="base", help=GENE_COLUMN_HELP)

    stats = {"proteins_in": n_in, "genes": len(rows), "too_short": n_short}
    el = time.time() - t0
    write_manifest(out, [proteins], params, stats, seconds=el, tables=[frag],
                   stage="s00_ingest")
    return StageResult("s00_ingest", out, stats, seconds=el,
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
    adds=tuple(Column(n, "string", h) for n, h in GENE_COLUMN_HELP.items()),
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
