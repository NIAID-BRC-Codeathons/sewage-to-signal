"""Stage 03 - gene calling.

pyrodigal in metagenomic mode: a pure wheel, no external binary, and the right
mode for mixed-organism contigs.

This is where the pipeline stops passing files and starts passing rows. The
stage emits the **gene** level: one row per called gene, carrying coordinates,
the amino acid sequence and a content hash. Every later stage adds columns to
these rows rather than writing a filtered FASTA, and every later stage picks
its input with a predicate over them.

The sequence is a column because it makes any subset reconstructible: a
predicate is enough to rebuild the exact FASTA any stage was given, which is
strictly more auditable than the fixed handful of subsets the pipeline used to
materialise. A protein FASTA is still written to the work directory, for tools
outside this pipeline.

The hash is ``sha1`` of the sequence, and it is what makes cross-sample
predicates possible - ``gene_id`` is per-assembly, so without it "also dark in
CHI-A" has nothing to join on.
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path

import lake
from common import StageResult, read_fasta, workdir, write_fasta
from entities import gene_columns, gene_row, gene_schema, write_rows
from stage import Param, Stage


def run(
    contigs: Path,
    out_dir: Path,
    sample: str,
    con=None,
    min_aa: int = 60,
    keep_partial: bool = True,
    force: bool = False,
) -> StageResult:
    contigs = Path(contigs)
    faa = Path(out_dir) / f"{sample}.proteins.faa"
    params = {"min_aa": min_aa, "keep_partial": keep_partial}
    deps = [lake.fingerprint(contigs)]
    if not force and lake.is_current(con, sample, "s03_genes", params, deps):
        return StageResult("s03_genes", faa, {}, skipped=True,
                           produced={"gene": None})

    import pyarrow as pa
    import pyrodigal

    t0 = time.time()
    finder = pyrodigal.GeneFinder(meta=True)
    records, rows = [], []
    n_contigs = n_called = n_partial = 0
    for header, seq in read_fasta(contigs):
        n_contigs += 1
        contig_id = header.split()[0]
        for i, gene in enumerate(finder.find_genes(seq.encode()), start=1):
            n_called += 1
            aa = gene.translate().rstrip("*")
            partial = bool(gene.partial_begin or gene.partial_end)
            if partial:
                n_partial += 1
            if len(aa) < min_aa or (partial and not keep_partial):
                continue
            gid = f"{contig_id}_{i}"
            strand = "+" if gene.strand > 0 else "-"
            records.append((f"{gid} {gene.begin}-{gene.end}({strand}) partial={int(partial)}", aa))
            rows.append(gene_row(gid, aa, contig=contig_id, begin=gene.begin,
                                 end=gene.end, strand=strand, partial=partial))

    write_fasta(faa, records)
    table = pa.Table.from_pylist(rows, schema=gene_schema())
    n = write_rows(con, "gene", sample, table, STAGE)

    stats = {
        "contigs": n_contigs, "genes_called": n_called, "genes_kept": len(rows),
        "partial": n_partial,
        "mean_aa": round(sum(r["aa_len"] for r in rows) / len(rows), 1) if rows else 0,
    }
    el = time.time() - t0
    lake.record_run(con, sample, "s03_genes", params, deps, None, stats,
                    tools={"pyrodigal": pyrodigal.__version__}, seconds=el, rows=n)
    return StageResult("s03_genes", faa, stats, seconds=el, produced={"gene": None})


STAGE = Stage(
    name="s03_genes",
    title="Call genes",
    summary="pyrodigal in metagenomic mode. Emits the gene level: one row per "
            "called gene with coordinates, sequence and content hash.",
    run=run,
    consumes="contigs",
    produces="gene",
    order=30,
    roles=("gene_source",),
    adds=gene_columns(),
    params=(
        Param("min_aa", int, 60, group="genes",
              help="drop genes shorter than this many amino acids"),
        Param("keep_partial", bool, True, group="genes",
              help="keep genes that run off the end of their contig"),
    ),
    requires=(),
)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("contigs", type=Path)
    p.add_argument("--sample", required=True)
    p.add_argument("--work", type=Path, default=Path("work"))
    p.add_argument("--min-aa", type=int, default=60)
    p.add_argument("--drop-partial", action="store_true")
    p.add_argument("--force", action="store_true")
    a = p.parse_args()
    r = run(a.contigs, workdir(a.work, a.sample, "s03_genes"), a.sample,
            min_aa=a.min_aa, keep_partial=not a.drop_partial, force=a.force)
    print(r.describe())


if __name__ == "__main__":
    main()
