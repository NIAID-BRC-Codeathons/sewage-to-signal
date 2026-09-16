"""Stage 03 - gene calling.

pyrodigal in metagenomic mode: a pure wheel, no external binary, and the right
mode for mixed-organism contigs. Emits a protein FASTA plus a TSV of
coordinates, strand and partial-gene flags.
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path

from common import StageResult, is_current, read_fasta, workdir, write_fasta, write_manifest


def run(
    contigs: Path,
    out_dir: Path,
    sample: str,
    min_aa: int = 60,
    keep_partial: bool = True,
    force: bool = False,
) -> StageResult:
    contigs = Path(contigs)
    out = Path(out_dir) / f"{sample}.proteins.faa"
    coords = Path(out_dir) / f"{sample}.genes.tsv"
    params = {"min_aa": min_aa, "keep_partial": keep_partial}
    if not force and is_current(out, [contigs], params):
        return StageResult("s03_genes", out, {}, skipped=True)

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
            rows.append((gid, contig_id, gene.begin, gene.end, strand, int(partial), len(aa)))

    write_fasta(out, records)
    with open(coords, "w") as fh:
        fh.write("gene_id\tcontig\tbegin\tend\tstrand\tpartial\taa_len\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")

    stats = {
        "contigs": n_contigs, "genes_called": n_called, "genes_kept": len(records),
        "partial": n_partial,
        "mean_aa": round(sum(r[6] for r in rows) / len(rows), 1) if rows else 0,
    }
    el = time.time() - t0
    write_manifest(out, [contigs], params, stats,
                   tools={"pyrodigal": pyrodigal.__version__}, seconds=el)
    return StageResult("s03_genes", out, stats, seconds=el)


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
