"""Stage 02 - assembly.

Reads are 151 bp; a 50 aa translated fragment is too short to carry a domain
and one indel error garbles the frame. Assembling first is what makes the
downstream protein language model meaningful, so this stage is not optional
for read input - but it does require an external assembler.

No pure-Python metagenome assembler exists, so this stage shells out to MEGAHIT
(fast, low memory) or metaSPAdes (slower, better contiguity).
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import time
from pathlib import Path

from stage import Param, Stage, Tool
from common import MissingTool, StageResult, is_current, read_fasta, which, workdir, write_manifest

INSTALL_HINT = (
    "Install an assembler, e.g.\n"
    "    conda install -c bioconda megahit\n"
    "  (this machine needs `conda tos accept` first - it is interactive),\n"
    "  or `brew install megahit`, or run this stage on a Linux host.\n"
    "  Alternatively skip to s03_genes with contigs you assembled elsewhere."
)


def run(
    fastq: Path,
    out_dir: Path,
    sample: str,
    fastq2: Path | None = None,
    min_contig: int = 500,
    threads: int = 4,
    memory_frac: float = 0.5,
    assembler: str = "auto",
    force: bool = False,
) -> StageResult:
    fastq = Path(fastq)
    out = Path(out_dir) / f"{sample}.contigs.fa"
    paired = fastq2 is not None
    deps = [fastq] + ([Path(fastq2)] if paired else [])
    params = {"min_contig": min_contig, "assembler": assembler,
              "threads": threads, "paired": paired}
    if not force and is_current(out, deps, params):
        return StageResult("s02_assemble", out, {}, skipped=True)

    tool = assembler
    if tool == "auto":
        tool = "megahit" if which("megahit") else ("spades" if which("spades.py") else "")
    if not tool or not which("megahit" if tool == "megahit" else "spades.py"):
        raise MissingTool(tool or "megahit/spades.py", INSTALL_HINT)

    t0 = time.time()
    tmp = Path(out_dir) / "_asm"
    if tmp.exists():
        shutil.rmtree(tmp)
    if tool == "megahit":
        # Paired input (-1/-2) gives far better contiguity. The NCBI fastq
        # endpoint returns concatenated reads rather than split mates, so those
        # downloads are single-end (--read); the ENA-streamed CASPER runs are
        # properly paired.
        cmd = ["megahit", "-o", str(tmp), "--min-contig-len", str(min_contig),
               "-t", str(threads), "--memory", str(memory_frac)]
        cmd += (["-1", str(fastq), "-2", str(fastq2)] if paired
                else ["--read", str(fastq)])
        subprocess.run(cmd, check=True, capture_output=True)
        produced = tmp / "final.contigs.fa"
    else:
        cmd = ["spades.py", "--meta", "-o", str(tmp), "-t", str(threads)]
        cmd += (["-1", str(fastq), "-2", str(fastq2)] if paired
                else ["-s", str(fastq)])
        subprocess.run(cmd, check=True, capture_output=True)
        produced = tmp / "contigs.fasta"
    shutil.copy(produced, out)

    lengths = [len(s) for _, s in read_fasta(out)]
    lengths.sort(reverse=True)
    total = sum(lengths)
    n50, acc = 0, 0
    for L in lengths:
        acc += L
        if acc >= total / 2:
            n50 = L
            break
    stats = {
        "assembler": tool, "contigs": len(lengths), "total_bp": total,
        "n50": n50, "longest": lengths[0] if lengths else 0,
    }
    el = time.time() - t0
    write_manifest(out, deps, params, stats, tools={tool: "external"}, seconds=el)
    return StageResult("s02_assemble", out, stats, seconds=el)


STAGE = Stage(
    name="s02_assemble",
    title="Assemble",
    summary="MEGAHIT or metaSPAdes. Not optional for read input: a 50 aa "
            "translated fragment is too short to carry a domain, and one indel "
            "garbles the frame.",
    run=run,
    consumes="reads",
    produces="contigs",
    order=20,
    input_arg="fastq",
    params=(
        Param("min_contig", int, 500, group="assembly",
              help="drop contigs shorter than this"),
        Param("assembler", str, "auto", choices=("auto", "megahit", "spades"),
              group="assembly",
              help="megahit is fast and low memory; spades is slower with "
                   "better contiguity"),
        Param("memory_frac", float, 0.5, group="assembly",
              help="fraction of system memory the assembler may use"),
    ),
    requires=(Tool("megahit", optional=True, hint="or metaspades.py"),),
)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("fastq", type=Path)
    p.add_argument("--fastq2", type=Path, help="second mate for paired input")
    p.add_argument("--sample", required=True)
    p.add_argument("--work", type=Path, default=Path("work"))
    p.add_argument("--min-contig", type=int, default=500)
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--assembler", default="auto", choices=["auto", "megahit", "spades"])
    p.add_argument("--force", action="store_true")
    a = p.parse_args()
    r = run(a.fastq, workdir(a.work, a.sample, "s02_assemble"), a.sample,
            fastq2=a.fastq2, min_contig=a.min_contig, threads=a.threads, assembler=a.assembler,
            force=a.force)
    print(r.describe())


if __name__ == "__main__":
    main()
