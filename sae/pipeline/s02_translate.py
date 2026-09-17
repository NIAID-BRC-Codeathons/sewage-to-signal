"""Stage 02 (alternative) - six-frame translation, skipping the assembler.

The assembly route exists because a 151 bp read is a 50 aa peptide: too short
to carry a domain, and one indel garbles the frame. That argument has not
changed, and this stage does not answer it. What it answers is a different
question - *what is in this sample, roughly, right now* - for which waiting on
MEGAHIT, and on having MEGAHIT installed at all, is the wrong trade.

So this is the short road: reads to peptides to embeddings, with no assembly
and no external binary. Each read is translated in all six frames and each
frame is cut at its stop codons, which is the standard translated-search shape
(the same one DIAMOND blastx uses) and is what makes the length filter do real
work: a frame that is not coding hits a stop every ~21 codons on average, so
its pieces fall under the floor and vanish, while the one frame that is coding
yields a single long piece. The filter is a frame-selector that costs nothing.

**Read this before believing a result from it.** Every peptide here is at most
``read_length / 3`` residues - 50 aa for a 151 bp read - so nothing it produces
can contain a complete domain. s05 will call almost all of it ``dark``, not
because the proteins are novel but because they are fragments. Treat the map it
produces as a triage picture of what the sample is made of, and assemble before
claiming anything about a protein.

Measured on SRR38294894 (CASPER influent RNA-seq, 151 bp):

| min_aa | peptides per read | unique after exact dedup |
|---|---|---|
| 30 | 2.97 | 1.80 |
| 40 | 1.51 | 1.05 |
| 50 | 0.51 | 0.41 |

Translation runs at ~38,000 reads/s, so it is never the cost. The GPU is: at
the default floor a million reads is ~1.5 M peptides, which is days of ESMC
even at 300M throughput. Cap the input (``--set s01_qc.max_reads=``) or the
batch (``--set s06_embed.limit=``) rather than discovering this on the GPU.

Duplicates are left in. s04_derep already collapses them on the gene level and
records how many collapsed, and doing it here as well would throw away the
count while re-implementing the stage that owns it.
"""

from __future__ import annotations

import argparse
import itertools
import time
from pathlib import Path

import lake
from common import (StageResult, open_maybe_gzip, workdir,
                    write_fasta)
from stage import Param, Stage

# The standard genetic code, laid out in TCAG order so the 64 codons fall out
# of one product rather than a hand-typed table nobody can proofread.
_BASES = "TCAG"
_AAS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODONS = {"".join(c): _AAS[i]
          for i, c in enumerate(itertools.product(_BASES, repeat=3))}
COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def peptides(seq: str, min_aa: int):
    """Stop-free peptides of at least `min_aa` from all six frames.

    Yields (strand, frame, peptide). A codon with an ambiguous base has no
    translation, so it becomes X and the piece carrying it is dropped by the
    caller - after QC those are rare, and a peptide with an unknown residue in
    it is not worth a forward pass.
    """
    seq = seq.upper()
    rc = seq.translate(COMPLEMENT)[::-1]
    for strand, s in ((1, seq), (-1, rc)):
        for frame in range(3):
            aa = "".join([CODONS.get(s[i:i + 3], "X")
                          for i in range(frame, len(s) - 2, 3)])
            for piece in aa.split("*"):
                if len(piece) >= min_aa:
                    yield strand, frame, piece


def _fastq_seqs(path: Path, limit: int | None = None):
    """Sequence lines of a FASTQ, streamed. Four lines per record."""
    n = 0
    with open_maybe_gzip(path) as fh:
        for i, line in enumerate(fh):
            if i % 4 != 1:
                continue
            yield line.strip()
            n += 1
            if limit and n >= limit:
                return


def run(
    fastq: Path,
    out_dir: Path,
    sample: str,
    con=None,
    fastq2: Path | None = None,
    min_aa: int = 40,
    max_reads: int | None = None,
    force: bool = False,
) -> StageResult:
    fastq = Path(fastq)
    fastq2 = Path(fastq2) if fastq2 else None
    out = Path(out_dir) / f"{sample}.translated.faa"
    params = {"min_aa": min_aa, "max_reads": max_reads,
              "paired": fastq2 is not None}
    deps = [lake.fingerprint(f) for f in [fastq] + ([fastq2] if fastq2 else [])]
    if not force and lake.is_current(con, sample, "s02_translate", params, deps):
        return StageResult("s02_translate", out, {}, skipped=True,
                           produced={"proteins": out})

    t0 = time.time()
    n_reads = n_kept = n_ambiguous = 0
    longest = 0

    def records():
        """(id, peptide) for every kept frame of every read, both mates.

        The id is a counter rather than the read name: a mate pair shares one
        name, a FASTQ may repeat one, and a duplicate id would break every
        join downstream - s00_ingest rejects them outright. The read it came
        from stays in the description, which nothing keys on but a person
        reading the FASTA wants.
        """
        nonlocal n_reads, n_kept, n_ambiguous, longest
        mates = [(fastq, 1)] + ([(fastq2, 2)] if fastq2 else [])
        for path, mate in mates:
            for seq in _fastq_seqs(path, max_reads):
                n_reads += 1
                for strand, frame, pep in peptides(seq, min_aa):
                    if "X" in pep:
                        n_ambiguous += 1
                        continue
                    n_kept += 1
                    longest = max(longest, len(pep))
                    yield (f"{sample}_t{n_kept:08d}",
                           f"mate={mate} strand={strand:+d} frame={frame}"), pep

    written = write_fasta(out, ((f"{i} {d}", p) for (i, d), p in records()))
    if not n_kept:
        raise SystemExit(
            f"no peptides of at least {min_aa} aa from {n_reads} reads; "
            f"lower min_aa or check that the input is nucleotide reads")

    el = time.time() - t0
    stats = {
        "reads_in": n_reads, "peptides": written,
        "per_read": round(written / n_reads, 2) if n_reads else 0,
        "dropped_ambiguous": n_ambiguous, "longest_aa": longest,
        "reads_per_s": round(n_reads / el) if el else None,
    }
    lake.record_run(con, sample, "s02_translate", params, deps, None,
                    stats, seconds=el)
    return StageResult("s02_translate", out, stats, seconds=el,
                       produced={"proteins": out})


STAGE = Stage(
    name="s02_translate",
    title="Translate reads",
    summary="Six-frame translation of reads into peptides, cut at stop codons. "
            "The no-assembly route: fast and binary-free, but every peptide is "
            "a read-length fragment and cannot carry a complete domain.",
    run=run,
    consumes="reads",
    produces="proteins",
    # Above s02_assemble, so that when both routes reach the gene level in the
    # same number of steps - reads->contigs->gene against reads->proteins->gene
    # - the planner still picks assembly. This one is asked for, never chosen.
    order=25,
    roles=("translation",),
    input_arg="fastq",
    params=(
        Param("min_aa", int, 40, group="translate",
              help="drop peptides shorter than this; the floor is what selects "
                   "the coding frame, so lowering it mostly adds junk"),
        Param("max_reads", int, None, group="translate",
              help="translate only the first N reads of each mate"),
    ),
    requires=(),
)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("fastq", type=Path)
    p.add_argument("--fastq2", type=Path)
    p.add_argument("--sample", required=True)
    p.add_argument("--work", type=Path, default=Path("work"))
    p.add_argument("--min-aa", type=int, default=40)
    p.add_argument("--max-reads", type=int)
    p.add_argument("--force", action="store_true")
    a = p.parse_args()
    r = run(a.fastq, workdir(a.work, a.sample, "s02_translate"), a.sample,
            fastq2=a.fastq2, min_aa=a.min_aa, max_reads=a.max_reads,
            force=a.force)
    print(r.describe())


if __name__ == "__main__":
    main()
