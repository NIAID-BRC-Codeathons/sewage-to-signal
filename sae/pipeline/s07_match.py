"""Stage 07 - annotate features and match clusters.

Joins the parquet from s06 against the published SAE feature table (adding
human-readable descriptions) and then against the local ESM Atlas
representatives.

The cluster join bridges through UniRef accessions because the local tables
carry no SAE-feature column. That bridge covers only ~2.3% of the feature
table's accessions, so treat recovered clusters as leads, not as a survey.
The durable fix is to run s06 over the atlas representatives themselves and
build a real feature -> cluster index; see the pipeline README.
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path

from common import StageResult, workdir, write_manifest

FEATURE_REPO = "biohub/ESMC-SAE-Features"
FEATURE_FILE = "uniref90_feature_table.parquet"


def load_features(feature_ids):
    import pyarrow as pa
    import pyarrow.compute as pc
    import pyarrow.parquet as pq
    from huggingface_hub import hf_hub_download

    path = hf_hub_download(FEATURE_REPO, FEATURE_FILE, repo_type="dataset")
    t = pq.read_table(path, columns=[
        "feature_id", "summary", "category", "exemplar_protein_families",
        "threshold", "uniref90_frequency", "top_100_uniref_ids",
    ])
    keep = pc.is_in(t.column("feature_id"), value_set=pa.array(sorted(feature_ids)))
    return {r["feature_id"]: r for r in t.filter(keep).to_pylist()}


def run(
    features_parquet: Path,
    out_dir: Path,
    sample: str,
    reps: Path,
    uniref_per_feature: int = 100,
    skip_clusters: bool = False,
    force: bool = False,
) -> StageResult:
    import pyarrow as pa
    import pyarrow.compute as pc
    import pyarrow.dataset as ds
    import pyarrow.parquet as pq

    features_parquet = Path(features_parquet)
    out = Path(out_dir) / f"{sample}.annotated.parquet"
    clusters_out = Path(out_dir) / f"{sample}.clusters.parquet"
    t0 = time.time()

    feats = pq.read_table(features_parquet)
    fids = set(feats.column("feature_id").to_pylist())
    table = load_features(fids)

    ann = feats.append_column(
        "summary",
        pa.array([(table.get(f) or {}).get("summary") for f in feats.column("feature_id").to_pylist()]),
    ).append_column(
        "category",
        pa.array([(table.get(f) or {}).get("category") for f in feats.column("feature_id").to_pylist()]),
    )
    pq.write_table(ann, out, compression="zstd")
    stats = {"rows": ann.num_rows, "distinct_features": len(fids)}

    if not skip_clusters and Path(reps).exists():
        nominations: dict[int, dict[str, float]] = {}
        for fid in fids:
            row = table.get(fid)
            if not row:
                continue
            nominations[fid] = {
                e["uniref_id"]: e["activation"]
                for e in (row.get("top_100_uniref_ids") or [])[:uniref_per_feature]
            }
        wanted = {a for m in nominations.values() for a in m}
        if wanted:
            probes = pa.array(sorted(
                {f"{p}_{a}" for a in wanted for p in ("UniRef90", "UniRef100")}
            ))
            matched = ds.dataset(reps, format="parquet").to_table(
                columns=["protein_hash", "uniref_match_accession", "lca_taxonomy",
                         "product_name", "cluster_top_pfam_names"],
                filter=ds.field("uniref_match_accession").isin(probes),
            )
            acc_to_feats: dict[str, list[int]] = {}
            for fid, accs in nominations.items():
                for a in accs:
                    acc_to_feats.setdefault(a, []).append(fid)
            rows = []
            for r in matched.to_pylist():
                bare = (r["uniref_match_accession"] or "").split("_", 1)[-1]
                for fid in acc_to_feats.get(bare, []):
                    rows.append({
                        "feature_id": fid,
                        "cluster_rep_protein_hash": r["protein_hash"],
                        "uniref_match_accession": r["uniref_match_accession"],
                        "lca_taxonomy": r["lca_taxonomy"],
                        "product_name": r["product_name"],
                        "pfam": "; ".join(f"{a} ({b})" for a, b in (r["cluster_top_pfam_names"] or [])[:3]),
                    })
            if rows:
                pq.write_table(pa.Table.from_pylist(rows), clusters_out, compression="zstd")
            stats["cluster_rows"] = len(rows)
            stats["uniref_probed"] = len(wanted)

    el = time.time() - t0
    write_manifest(out, [features_parquet], {"uniref_per_feature": uniref_per_feature},
                   stats, seconds=el)
    return StageResult("s07_match", out, stats, seconds=el)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("features", type=Path)
    p.add_argument("--sample", required=True)
    p.add_argument("--work", type=Path, default=Path("work"))
    p.add_argument("--reps", type=Path,
                   default=Path(__file__).resolve().parent.parent / "representative_proteins.parquet")
    p.add_argument("--skip-clusters", action="store_true")
    a = p.parse_args()
    r = run(a.features, workdir(a.work, a.sample, "s07_match"), a.sample,
            reps=a.reps, skip_clusters=a.skip_clusters)
    print(r.describe())


if __name__ == "__main__":
    main()
