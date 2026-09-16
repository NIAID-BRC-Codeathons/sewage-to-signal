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

# Backbone, SAE repo and layer have to agree; picking them separately is an
# easy way to get a silently wrong answer, so offer them as one choice.
# ESMC-6B needs ~12 GB for weights alone — on anything smaller, or on CPU,
# use 300m. Only the 6B layer-60 SAE has a published feature description
# table, so 300m gives retrieval but no human-readable summaries in s07.
MODELS = {
    "6b": ("biohub/ESMC-6B",
           "biohub/ESMC-6B-sae-layer60-k64-codebook16384", 60),
    "300m": ("biohub/ESMC-300M",
             "biohub/ESMC-300M-sae-layer23-k64-codebook16384", 23),
}


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
    p.add_argument("--evalue", type=float, default=1e-5)
    p.add_argument("--confident-evalue", type=float, default=1e-20)
    p.add_argument("--min-coverage", type=float, default=0.80)
    p.add_argument("--bit-cutoffs", choices=["gathering", "noise", "trusted"])
    p.add_argument("--model", choices=sorted(MODELS), default="6b",
                   help="6b needs ~12 GB and a GPU to be practical; 300m runs "
                        "on CPU but has no feature description table (default: 6b)")
    p.add_argument("--backbone", help="override the --model backbone")
    p.add_argument("--sae-repo", help="override the --model SAE repo")
    p.add_argument("--layer", type=int, help="override the --model SAE layer")
    p.add_argument("--max-len", type=int, default=1022)
    p.add_argument("--top-k", type=int, default=16)
    p.add_argument("--batch-size", type=int, default=8)
    p.add_argument("--limit", type=int, help="cap proteins sent to the GPU stage")
    p.add_argument("--device", default="auto")
    p.add_argument("--reps", type=Path,
                   default=Path(__file__).resolve().parent.parent / "representative_proteins.parquet")
    p.add_argument("--force", action="store_true")
    a = p.parse_args()

    backbone, sae_repo, layer = MODELS[a.model]
    backbone = a.backbone or backbone
    sae_repo = a.sae_repo or sae_repo
    layer = layer if a.layer is None else a.layer

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
                                      hmm=a.hmm, ref=a.ref, evalue=a.evalue,
                                      confident_evalue=a.confident_evalue,
                                      min_coverage=a.min_coverage,
                                      bit_cutoffs=a.bit_cutoffs, force=a.force)
            elif stage == "s06_embed":
                r = s06_embed.run(current, wd(stage), a.sample, top_k=a.top_k,
                                  backbone=backbone, sae_repo=sae_repo,
                                  layer=layer, max_len=a.max_len,
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
