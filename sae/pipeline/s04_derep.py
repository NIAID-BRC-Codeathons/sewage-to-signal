"""Stage 04 - dereplication.

Wastewater surveillance resequences the same sewershed over and over, so the
protein set across samples is heavily redundant. Collapsing it before the GPU
stage is the cheapest large saving in the whole pipeline.

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

from common import StageResult, is_current, read_fasta, which, workdir, write_fasta, write_manifest


def run(
    proteins: Path,
    out_dir: Path,
    sample: str,
    identity: float = 0.95,
    coverage: float = 0.8,
    threads: int = 4,
    use_mmseqs: bool = True,
    force: bool = False,
) -> StageResult:
    proteins = Path(proteins)
    out = Path(out_dir) / f"{sample}.nr.faa"
    mapping = Path(out_dir) / f"{sample}.derep.tsv"
    engine = "mmseqs" if (use_mmseqs and which("mmseqs")) else "exact"
    params = {"identity": identity, "coverage": coverage, "engine": engine}
    if not force and is_current(out, [proteins], params):
        return StageResult("s04_derep", out, {"engine": engine}, skipped=True)

    t0 = time.time()
    # Exact dedup first: it is free and shrinks the input to any clusterer.
    seen: dict[str, str] = {}
    members: dict[str, list[str]] = {}
    n_in = 0
    for header, seq in read_fasta(proteins):
        n_in += 1
        gid = header.split()[0]
        h = hashlib.sha1(seq.encode()).hexdigest()
        if h in seen:
            members[seen[h]].append(gid)
            continue
        seen[h] = gid
        members[gid] = [gid]
    uniq = [(g, s) for (g, s) in ((seen[h], None) for h in seen)]
    # re-read to pull the representative sequences in one pass
    rep_ids = set(seen.values())
    reps = [(hd.split()[0], sq) for hd, sq in read_fasta(proteins) if hd.split()[0] in rep_ids]
    seen_once, deduped = set(), []
    for gid, sq in reps:
        if gid not in seen_once:
            seen_once.add(gid)
            deduped.append((gid, sq))

    stats = {"proteins_in": n_in, "after_exact": len(deduped)}

    if engine == "mmseqs":
        tmp = Path(out_dir) / "_mm"
        tmp.mkdir(exist_ok=True)
        exact_fa = tmp / "exact.faa"
        write_fasta(exact_fa, deduped)
        pref = tmp / "clu"
        subprocess.run(
            ["mmseqs", "easy-linclust", str(exact_fa), str(pref), str(tmp / "tmp"),
             "--min-seq-id", str(identity), "-c", str(coverage),
             "--threads", str(threads)],
            check=True, capture_output=True,
        )
        rep_fa = Path(str(pref) + "_rep_seq.fasta")
        final = list(read_fasta(rep_fa))
        final = [(h.split()[0], s) for h, s in final]
        clu_tsv = Path(str(pref) + "_cluster.tsv")
        if clu_tsv.exists():
            mapping.write_text(clu_tsv.read_text())
        stats["after_cluster"] = len(final)
    else:
        final = deduped
        with open(mapping, "w") as fh:
            fh.write("representative\tmember\n")
            for rep, mem in members.items():
                for m in mem:
                    fh.write(f"{rep}\t{m}\n")

    write_fasta(out, final)
    stats["representatives"] = len(final)
    stats["reduction"] = round(1 - len(final) / n_in, 4) if n_in else 0.0
    stats["engine"] = engine
    el = time.time() - t0
    write_manifest(out, [proteins], params, stats, seconds=el)
    return StageResult("s04_derep", out, stats, seconds=el)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("proteins", type=Path)
    p.add_argument("--sample", required=True)
    p.add_argument("--work", type=Path, default=Path("work"))
    p.add_argument("--identity", type=float, default=0.95)
    p.add_argument("--no-mmseqs", action="store_true")
    p.add_argument("--force", action="store_true")
    a = p.parse_args()
    r = run(a.proteins, workdir(a.work, a.sample, "s04_derep"), a.sample,
            identity=a.identity, use_mmseqs=not a.no_mmseqs, force=a.force)
    print(r.describe())


if __name__ == "__main__":
    main()
