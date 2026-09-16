"""Pipeline driver.

Runs stages in order, starting from whatever input you have. Each stage is
idempotent, so re-running resumes rather than recomputing.

    python run.py --fastq ../../data/fastq/SRR40033132.fastq.gz --sample S1
    python run.py --contigs contigs.fa --sample S1 --from s03_genes
    python run.py --proteins proteins.faa --sample S1 --from s06_embed
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import s01_qc, s02_assemble, s03_genes, s04_derep, s05_prefilter, s06_embed, s07_match
from common import MissingTool, workdir

ORDER = ["s01_qc", "s02_assemble", "s03_genes", "s04_derep",
         "s05_prefilter", "s06_embed", "s07_match"]


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("--fastq", type=Path)
    p.add_argument("--fastq2", type=Path, help="second mate (with --fastq)")
    src.add_argument("--contigs", type=Path)
    src.add_argument("--proteins", type=Path)
    p.add_argument("--sample", required=True)
    p.add_argument("--work", type=Path, default=Path("work"))
    p.add_argument("--from", dest="start", choices=ORDER)
    p.add_argument("--to", dest="stop", choices=ORDER, default="s07_match")
    p.add_argument("--max-reads", type=int)
    p.add_argument("--min-aa", type=int, default=60)
    p.add_argument("--hmm", type=Path)
    p.add_argument("--ref", type=Path)
    p.add_argument("--top-k", type=int, default=16)
    p.add_argument("--batch-size", type=int, default=8)
    p.add_argument("--limit", type=int, help="cap proteins sent to the GPU stage")
    p.add_argument("--device", default="auto")
    p.add_argument("--reps", type=Path,
                   default=Path(__file__).resolve().parent.parent / "representative_proteins.parquet")
    p.add_argument("--force", action="store_true")
    a = p.parse_args()

    if a.fastq:
        start = a.start or "s01_qc"
        current = a.fastq
    elif a.contigs:
        start = a.start or "s03_genes"
        current = a.contigs
    else:
        start = a.start or "s04_derep"
        current = a.proteins

    begin, end = ORDER.index(start), ORDER.index(a.stop)
    if begin > end:
        sys.exit(f"--from {start} comes after --to {a.stop}")

    mate = a.fastq2
    wd = lambda st: workdir(a.work, a.sample, st)
    results = []
    for stage in ORDER[begin : end + 1]:
        try:
            if stage == "s01_qc":
                r = s01_qc.run(current, wd(stage), a.sample, fastq2=mate,
                               max_reads=a.max_reads, force=a.force)
                mate = r.mate
            elif stage == "s02_assemble":
                r = s02_assemble.run(current, wd(stage), a.sample,
                                     fastq2=mate, force=a.force)
            elif stage == "s03_genes":
                r = s03_genes.run(current, wd(stage), a.sample,
                                  min_aa=a.min_aa, force=a.force)
            elif stage == "s04_derep":
                r = s04_derep.run(current, wd(stage), a.sample, force=a.force)
            elif stage == "s05_prefilter":
                r = s05_prefilter.run(current, wd(stage), a.sample,
                                      hmm=a.hmm, ref=a.ref, force=a.force)
            elif stage == "s06_embed":
                r = s06_embed.run(current, wd(stage), a.sample, top_k=a.top_k,
                                  batch_size=a.batch_size, limit=a.limit,
                                  device=a.device, force=a.force)
            else:
                r = s07_match.run(current, wd(stage), a.sample, reps=a.reps)
        except MissingTool as exc:
            print(f"\n[{stage}] BLOCKED: {exc}\n", file=sys.stderr)
            print("Completed stages:", file=sys.stderr)
            for done in results:
                print("  " + done.describe(), file=sys.stderr)
            sys.exit(2)
        print(r.describe(), flush=True)
        results.append(r)
        current = r.output
    print(f"\nDone. Final output: {current}")


if __name__ == "__main__":
    main()
