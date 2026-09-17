"""Stage 05 - homology triage.

The biggest lever in the pipeline. Homology search is orders of magnitude
cheaper than an ESMC-6B forward pass, so this stage decides what is worth
spending the GPU on. Proteins fall into three classes:

* ``known``   - a confident hit that explains *most of the protein*. Fully
                accounted for by conventional annotation, so nothing downstream
                asks for it.
* ``partial`` - a hit exists but is weak, or it covers only a fraction of the
                sequence. Either way the protein is not explained, so it is
                analysed - carrying its family label, which is the conditioning
                key for reference-relative (within-family) comparison.
* ``dark``    - no significant hit at all. Analysed with no prior.

Coverage is what separates ``known`` from ``partial``, and it is the reason a
single E-value threshold is not enough: a 1273 aa protein matching one 150 aa
domain at E=1e-40 has a superb E-value and still leaves 88% of its sequence
unexplained. Coverage is the union of all significant domain envelopes across
*all* families, so a genuine multi-domain protein is correctly called complete.

Backends, in preference order:
  --hmm  Pfam-A.hmm      -> pyhmmer  (domain-level, gives real coverage)
  --ref  reference.faa   -> pyswrd   (heuristic + Smith-Waterman)
Neither given -> everything is dark, with a warning. That is a pass-through,
not a filter, and it will cost you GPU time.

With Pfam-A, prefer ``--bit-cutoffs gathering``: Pfam ships curated per-family
thresholds, which are a better significance test than any flat E-value.

The stage writes columns, not FASTAs. It used to emit four - one per class plus
``analyze.faa`` - which were only ever a way of handing the next stage a
subset; the class was the real output and the files were a transport. Now
``category`` is a column and the GPU stage selects ``category <> 'known'``.
Anyone who wants a different cut (say, dark proteins over 200 aa that are also
dark in another sample) writes that predicate instead of asking for a fifth
file.
"""

from __future__ import annotations

import argparse
import time
from dataclasses import dataclass, field
from pathlib import Path

from common import (StageResult, is_current, read_fasta, workdir,
                    write_fragment, write_manifest)
from entities import fasta_records, table_from_fasta
from stage import Column, Param, Stage, Tool

CATEGORIES = ("known", "partial", "dark")

COLUMN_HELP = {
    "category": "known | partial | dark - see the module docstring",
    "family": "name of the best-scoring family, if any",
    "family_acc": "accession of that family",
    "evalue": "E-value of the best hit; null when there was none",
    "coverage": "fraction of the protein covered by all significant envelopes",
    "n_domains": "number of significant domain hits across all families",
}


def _txt(v) -> str | None:
    """pyhmmer returns names as str in 0.12.x and bytes in older builds."""
    if v is None:
        return None
    return v.decode() if isinstance(v, bytes) else str(v)


@dataclass
class Assignment:
    """Per-protein homology verdict."""

    gene_id: str
    aa_len: int
    family: str | None = None
    family_acc: str | None = None
    evalue: float | None = None
    coverage: float | None = None
    n_domains: int = 0
    envelopes: list[tuple[int, int]] = field(default_factory=list)
    category: str = "dark"

    def classify(self, confident_evalue: float, min_coverage: float) -> str:
        if self.n_domains == 0:
            self.category = "dark"
        elif self.coverage is None:
            # No coverage available (heuristic backend). Completeness cannot be
            # established, so never discard - fall through to analysis.
            self.category = "partial"
        elif self.evalue <= confident_evalue and self.coverage >= min_coverage:
            self.category = "known"
        else:
            self.category = "partial"
        return self.category


def _covered(intervals: list[tuple[int, int]]) -> int:
    """Total residues spanned by a set of 1-based inclusive intervals."""
    if not intervals:
        return 0
    total = 0
    cur_s, cur_e = (ivs := sorted(intervals))[0]
    for s, e in ivs[1:]:
        if s <= cur_e + 1:            # overlapping or abutting
            cur_e = max(cur_e, e)
        else:
            total += cur_e - cur_s + 1
            cur_s, cur_e = s, e
    return total + cur_e - cur_s + 1


def _hmm_assign(proteins, hmm_path: Path, threads: int, evalue: float,
                bit_cutoffs: str | None) -> dict[str, Assignment]:
    """Search every HMM, accumulating per-protein envelopes and the best hit.

    ``hmmsearch`` iterates over HMMs and yields one TopHits per model, so a
    protein's domains arrive spread across many iterations and have to be
    accumulated rather than read off a single result.
    """
    import pyhmmer

    alphabet = pyhmmer.easel.Alphabet.amino()
    seqs = [
        pyhmmer.easel.TextSequence(name=g.encode(), sequence=s).digitize(alphabet)
        for g, s in proteins
    ]
    out = {g: Assignment(gene_id=g, aa_len=len(s)) for g, s in proteins}

    opts = {"bit_cutoffs": bit_cutoffs} if bit_cutoffs else {"E": evalue}
    with pyhmmer.plan7.HMMFile(hmm_path) as hf:
        for hits in pyhmmer.hmmsearch(hf, seqs, cpus=threads, **opts):
            fam = _txt(hits.query.name)
            acc = _txt(hits.query.accession)
            for hit in hits:
                if not (hit.included if bit_cutoffs else hit.evalue <= evalue):
                    continue
                a = out[_txt(hit.name)]
                for d in hit.domains:
                    if not (d.included if bit_cutoffs else d.i_evalue <= evalue):
                        continue
                    a.envelopes.append((d.env_from, d.env_to))
                    a.n_domains += 1
                if a.evalue is None or hit.evalue < a.evalue:
                    a.evalue, a.family, a.family_acc = hit.evalue, fam, acc

    for a in out.values():
        if a.n_domains:
            a.coverage = round(_covered(a.envelopes) / a.aa_len, 4)
    return out


def _swrd_assign(proteins, ref_path: Path, evalue: float) -> dict[str, Assignment]:
    """Heuristic backend. Coverage only when the result carries alignment bounds.

    pyswrd returns one of several result types depending on the algorithm; only
    the full alignment carries query bounds. When they are absent, coverage is
    left as None and the protein can never be discarded.
    """
    import pyswrd

    names = [g for g, _ in proteins]
    out = {g: Assignment(gene_id=g, aa_len=len(s)) for g, s in proteins}
    targets = [s for _, s in read_fasta(ref_path)]
    for hit in pyswrd.search([s for _, s in proteins], targets, threads=0):
        if hit.evalue > evalue:
            continue
        a = out[names[hit.query_index]]
        a.n_domains += 1
        if a.evalue is None or hit.evalue < a.evalue:
            a.evalue = hit.evalue
            a.family = f"target_{hit.target_index}"
        r = getattr(hit, "result", None)
        qs, qe = getattr(r, "query_start", None), getattr(r, "query_end", None)
        if qs is not None and qe is not None:
            a.envelopes.append((qs + 1, qe + 1))    # pyopal bounds are 0-based

    for a in out.values():
        if a.envelopes:
            a.coverage = round(_covered(a.envelopes) / a.aa_len, 4)
    return out


def run(
    rows,
    out_dir: Path,
    sample: str,
    source: Path | None = None,
    where: str | None = None,
    hmm: Path | None = None,
    ref: Path | None = None,
    evalue: float = 1e-5,
    confident_evalue: float = 1e-20,
    min_coverage: float = 0.80,
    bit_cutoffs: str | None = None,
    threads: int = 4,
    force: bool = False,
) -> StageResult:
    out_dir = Path(out_dir)
    out = out_dir / f"{sample}.classification.parquet"

    backend = "pyhmmer" if hmm else ("pyswrd" if ref else "none")
    params = {
        "backend": backend, "evalue": evalue,
        "confident_evalue": confident_evalue, "min_coverage": min_coverage,
        "bit_cutoffs": bit_cutoffs, "where": where,
        "reference": str(hmm or ref) if (hmm or ref) else None,
    }
    deps = ([Path(source)] if source else []) \
        + ([Path(hmm)] if hmm else []) + ([Path(ref)] if ref else [])
    if not force and is_current(out, deps, params):
        return StageResult("s05_prefilter", out, {"backend": backend},
                           skipped=True, produced={"gene": None})

    import pyarrow as pa

    t0 = time.time()
    records = fasta_records(rows)
    if backend == "pyhmmer":
        assigned = _hmm_assign(records, Path(hmm), threads, evalue, bit_cutoffs)
    elif backend == "pyswrd":
        assigned = _swrd_assign(records, Path(ref), evalue)
    else:
        print(
            "  WARNING: no --hmm or --ref given; no triage applied. Every "
            "protein is passed to the GPU stage, which is the expensive path.",
        )
        assigned = {g: Assignment(gene_id=g, aa_len=len(s)) for g, s in records}

    for a in assigned.values():
        a.classify(confident_evalue, min_coverage)

    counts = {c: 0 for c in CATEGORIES}
    rows = []
    for gid, _ in records:
        a = assigned[gid]
        counts[a.category] += 1
        rows.append({
            "gene_id": gid, "category": a.category, "family": a.family,
            "family_acc": a.family_acc, "evalue": a.evalue,
            "coverage": a.coverage, "n_domains": a.n_domains,
        })
    table = pa.Table.from_pylist(rows, schema=pa.schema([
        ("gene_id", pa.string()), ("category", pa.string()),
        ("family", pa.string()), ("family_acc", pa.string()),
        ("evalue", pa.float64()), ("coverage", pa.float64()),
        ("n_domains", pa.int32()),
    ]))
    frag = write_fragment(out, table, "gene", where=where, help=COLUMN_HELP)

    n = len(records)
    n_analyze = n - counts["known"]
    stats = {
        "backend": backend, "proteins_in": n,
        "known": counts["known"], "partial": counts["partial"],
        "dark": counts["dark"], "analyzed": n_analyze,
        "discarded_frac": round(counts["known"] / n, 4) if n else 0.0,
        "analyzed_frac": round(n_analyze / n, 4) if n else 0.0,
    }
    el = time.time() - t0
    write_manifest(out, deps, params, stats, seconds=el, tables=[frag],
                   stage="s05_prefilter")
    return StageResult("s05_prefilter", out, stats, seconds=el,
                       produced={"gene": None})


STAGE = Stage(
    name="s05_prefilter",
    title="Homology triage",
    summary="Label each gene known / partial / dark by homology, with the "
            "coverage that separates a fully explained protein from one that "
            "merely has a hit. The pipeline's main cost lever.",
    run=run,
    consumes="gene",
    produces="gene",
    order=50,
    selectable=True,
    roles=("classification", "triage"),
    adds=(
        Column("category", "string", COLUMN_HELP["category"]),
        Column("family", "string", COLUMN_HELP["family"]),
        Column("family_acc", "string", COLUMN_HELP["family_acc"]),
        Column("evalue", "double", COLUMN_HELP["evalue"]),
        Column("coverage", "double", COLUMN_HELP["coverage"]),
        Column("n_domains", "int32", COLUMN_HELP["n_domains"]),
    ),
    params=(
        Param("hmm", str, None, group="prefilter", path=True,
              suffixes=(".hmm",), prefer_available=True,
              help="Pfam-A.hmm (or any HMM database) to search with pyhmmer"),
        Param("ref", str, None, group="prefilter", path=True,
              suffixes=(".faa", ".fasta", ".fa"),
              help="reference protein FASTA to search with pyswrd instead"),
        Param("evalue", float, 1e-5, group="prefilter",
              help="a hit below this is significant"),
        Param("confident_evalue", float, 1e-20, group="prefilter",
              help="a hit below this counts as confident"),
        Param("min_coverage", float, 0.80, group="prefilter",
              help="fraction a confident hit must span to count as explained"),
        Param("bit_cutoffs", str, None,
              choices=("gathering", "noise", "trusted"), group="prefilter",
              help="use the HMM's curated per-family thresholds instead of "
                   "--evalue; recommended with Pfam-A"),
    ),
    requires=(),
)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("proteins", type=Path)
    p.add_argument("--sample", required=True)
    p.add_argument("--work", type=Path, default=Path("work"))
    p.add_argument("--hmm", type=Path, help="Pfam-A.hmm for pyhmmer")
    p.add_argument("--ref", type=Path, help="reference protein FASTA for pyswrd")
    p.add_argument("--evalue", type=float, default=1e-5,
                   help="a hit below this is significant (default 1e-5)")
    p.add_argument("--confident-evalue", type=float, default=1e-20,
                   help="a hit below this counts as confident (default 1e-20)")
    p.add_argument("--min-coverage", type=float, default=0.80,
                   help="fraction of the protein a confident hit must span to "
                        "count as fully explained (default 0.80)")
    p.add_argument("--bit-cutoffs", choices=["gathering", "noise", "trusted"],
                   help="use the HMM's curated per-family thresholds instead "
                        "of --evalue for significance (recommended with Pfam-A)")
    p.add_argument("--force", action="store_true")
    a = p.parse_args()
    r = run(table_from_fasta(a.proteins), workdir(a.work, a.sample, "s05_prefilter"), a.sample,
            source=a.proteins, hmm=a.hmm, ref=a.ref, evalue=a.evalue,
            confident_evalue=a.confident_evalue, min_coverage=a.min_coverage,
            bit_cutoffs=a.bit_cutoffs, force=a.force)
    print(r.describe())


if __name__ == "__main__":
    main()
