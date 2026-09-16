"""Stage 05 - homology prefilter.

The biggest lever in the pipeline. Homology search is orders of magnitude
cheaper than an ESMC-6B forward pass, and a protein with a confident hit to a
known family does not need an SAE to identify it. The value of the SAE is in
the *unannotated* remainder, so this stage splits proteins into:

* ``<sample>.known.faa`` - confident hit, annotated conventionally
* ``<sample>.dark.faa``  - no hit; these go to the GPU

Backends, in preference order:
  --hmm  Pfam-A.hmm      -> pyhmmer  (domain-level, most informative)
  --ref  reference.faa   -> pyswrd   (fast heuristic + Smith-Waterman)
Neither given -> everything is treated as dark, with a warning. That is a
pass-through, not a filter, and it will cost you GPU time.
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path

from common import StageResult, is_current, read_fasta, workdir, write_fasta, write_manifest


def _hmm_hits(proteins, hmm_path: Path, threads: int, evalue: float) -> set[str]:
    import pyhmmer

    alphabet = pyhmmer.easel.Alphabet.amino()
    seqs = [
        pyhmmer.easel.TextSequence(name=g.encode(), sequence=s).digitize(alphabet)
        for g, s in proteins
    ]
    hit_ids: set[str] = set()
    with pyhmmer.plan7.HMMFile(hmm_path) as hf:
        for hits in pyhmmer.hmmsearch(hf, seqs, cpus=threads, E=evalue):
            for h in hits:
                if h.evalue <= evalue:
                    hit_ids.add(h.name.decode())
    return hit_ids


def _swrd_hits(proteins, ref_path: Path, evalue: float) -> set[str]:
    import pyswrd

    queries = [s for _, s in proteins]
    names = [g for g, _ in proteins]
    targets = [s for _, s in read_fasta(ref_path)]
    hit_ids: set[str] = set()
    for hit in pyswrd.search(queries, targets, threads=0):
        if hit.evalue <= evalue:
            hit_ids.add(names[hit.query_index])
    return hit_ids


def run(
    proteins: Path,
    out_dir: Path,
    sample: str,
    hmm: Path | None = None,
    ref: Path | None = None,
    evalue: float = 1e-5,
    threads: int = 4,
    force: bool = False,
) -> StageResult:
    proteins = Path(proteins)
    dark = Path(out_dir) / f"{sample}.dark.faa"
    known = Path(out_dir) / f"{sample}.known.faa"
    backend = "pyhmmer" if hmm else ("pyswrd" if ref else "none")
    params = {
        "backend": backend, "evalue": evalue,
        "reference": str(hmm or ref) if (hmm or ref) else None,
    }
    deps = [proteins] + ([Path(hmm)] if hmm else []) + ([Path(ref)] if ref else [])
    if not force and is_current(dark, deps, params):
        return StageResult("s05_prefilter", dark, {"backend": backend}, skipped=True)

    t0 = time.time()
    records = [(h.split()[0], s) for h, s in read_fasta(proteins)]
    if backend == "pyhmmer":
        hits = _hmm_hits(records, Path(hmm), threads, evalue)
    elif backend == "pyswrd":
        hits = _swrd_hits(records, Path(ref), evalue)
    else:
        print(
            "  WARNING: no --hmm or --ref given; no prefilter applied. Every "
            "protein is passed to the GPU stage, which is the expensive path.",
        )
        hits = set()

    dark_recs = [(g, s) for g, s in records if g not in hits]
    known_recs = [(g, s) for g, s in records if g in hits]
    write_fasta(dark, dark_recs)
    write_fasta(known, known_recs)

    stats = {
        "backend": backend, "proteins_in": len(records),
        "known": len(known_recs), "dark": len(dark_recs),
        "dark_frac": round(len(dark_recs) / len(records), 4) if records else 0.0,
    }
    el = time.time() - t0
    write_manifest(dark, deps, params, stats, seconds=el)
    return StageResult("s05_prefilter", dark, stats, seconds=el)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("proteins", type=Path)
    p.add_argument("--sample", required=True)
    p.add_argument("--work", type=Path, default=Path("work"))
    p.add_argument("--hmm", type=Path, help="Pfam-A.hmm for pyhmmer")
    p.add_argument("--ref", type=Path, help="reference protein FASTA for pyswrd")
    p.add_argument("--evalue", type=float, default=1e-5)
    p.add_argument("--force", action="store_true")
    a = p.parse_args()
    r = run(a.proteins, workdir(a.work, a.sample, "s05_prefilter"), a.sample,
            hmm=a.hmm, ref=a.ref, evalue=a.evalue, force=a.force)
    print(r.describe())


if __name__ == "__main__":
    main()
