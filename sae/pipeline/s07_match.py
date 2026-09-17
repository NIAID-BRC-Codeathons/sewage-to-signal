"""Stage 07 - annotate features and match clusters.

Joins the parquet from s06 against the published SAE feature table (adding
human-readable descriptions) and then against the local ESM Atlas
representatives.

It annotates the **feature** level, not the hit table: a summary is a property
of a codebook feature, so it is stored once per feature rather than once per
gene that happened to activate it. Anything wanting both joins them - which is
what the level key is for. Cluster matches are many-per-feature, so they are
their own level.

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

import lake
from common import StageResult, workdir
from entities import write_rows
from stage import Column, Param, Stage

FEATURE_HELP = {
    "summary": "human-readable description of what the feature detects",
    "feature_category": "the feature table's own category for it",
    "threshold": "activation the published table considers meaningful",
    "uniref90_frequency": "how common the feature is across UniRef90",
}
CLUSTER_HELP = {
    "cluster_rep_protein_hash": "ESM Atlas cluster representative",
    "uniref_match_accession": "accession that bridged feature to cluster",
    "lca_taxonomy": "lowest common ancestor of the cluster",
    "product_name": "product name recorded for the cluster",
    "pfam": "top Pfam names for the cluster",
}

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
    rows,
    out_dir: Path,
    sample: str,
    con=None,
    source: Path | None = None,
    where: str | None = None,
    reps: Path | None = None,
    uniref_per_feature: int = 100,
    skip_clusters: bool = False,
    force: bool = False,
) -> StageResult:
    import pyarrow as pa
    import pyarrow.dataset as ds

    out_dir = Path(out_dir)
    params = {"uniref_per_feature": uniref_per_feature,
              "skip_clusters": skip_clusters, "where": where}
    deps = [lake.fingerprint(reps)] if reps and Path(reps).exists() else []
    if not force and lake.is_current(con, sample, "s07_match", params, deps):
        return StageResult("s07_match", out_dir, {}, skipped=True,
                           produced={"feature": None, "cluster_hit": None})

    t0 = time.time()
    fids = sorted(set(rows.column("feature_id").to_pylist()))
    if not fids:
        raise SystemExit("the selection is empty - no features to annotate"
                         + (f" (predicate: {where})" if where else ""))
    table = load_features(fids)

    # One row per feature, not per hit: `category` is renamed because s05
    # already owns that name on the gene level, and two columns called
    # `category` meaning different things is exactly the confusion the
    # fragment model is meant to prevent.
    meta = pa.table({
        "feature_id": pa.array(fids, pa.int32()),
        "summary": pa.array([(table.get(f) or {}).get("summary") for f in fids]),
        "feature_category": pa.array([(table.get(f) or {}).get("category") for f in fids]),
        "threshold": pa.array([(table.get(f) or {}).get("threshold") for f in fids],
                              pa.float64()),
        "uniref90_frequency": pa.array(
            [(table.get(f) or {}).get("uniref90_frequency") for f in fids], pa.float64()),
    })
    write_rows(con, "feature", sample, meta, STAGE)
    stats = {"features": len(fids),
             "described": sum(1 for f in fids if f in table)}

    if not skip_clusters and reps and Path(reps).exists():
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
                write_rows(con, "cluster_hit", sample,
                           pa.Table.from_pylist(rows), STAGE)
            stats["cluster_rows"] = len(rows)
            stats["uniref_probed"] = len(wanted)

    el = time.time() - t0
    lake.record_run(con, sample, "s07_match", params, deps, where, stats,
                    seconds=el, rows=len(fids))
    return StageResult("s07_match", out_dir, stats, seconds=el,
                       produced={"feature": None, "cluster_hit": None})


STAGE = Stage(
    name="s07_match",
    title="Annotate features",
    summary="Join SAE features to the published description table, and bridge "
            "through UniRef to candidate ESM Atlas clusters.",
    run=run,
    consumes="feature",
    produces="feature",
    also_produces=("cluster_hit",),
    order=70,
    roles=("descriptions",),
    adds=(
        Column("summary", "string", FEATURE_HELP["summary"]),
        Column("feature_category", "string", FEATURE_HELP["feature_category"]),
        Column("threshold", "double", FEATURE_HELP["threshold"]),
        Column("uniref90_frequency", "double", FEATURE_HELP["uniref90_frequency"]),
    ),
    also_adds=(
        ("cluster_hit", (
            Column("cluster_rep_protein_hash", "string",
                   CLUSTER_HELP["cluster_rep_protein_hash"]),
            Column("uniref_match_accession", "string",
                   CLUSTER_HELP["uniref_match_accession"]),
            Column("lca_taxonomy", "string", CLUSTER_HELP["lca_taxonomy"]),
            Column("product_name", "string", CLUSTER_HELP["product_name"]),
            Column("pfam", "string", CLUSTER_HELP["pfam"]),
        )),
    ),
    params=(
        Param("reps", str, None, group="match", path=True,
              suffixes=(".parquet",),
              help="ESM Atlas representative_proteins.parquet; without it the "
                   "cluster bridge is skipped"),
        Param("uniref_per_feature", int, 100, group="match",
              help="UniRef accessions probed per feature"),
        Param("skip_clusters", bool, False, group="match",
              help="descriptions only; do not touch the atlas tables"),
    ),
    requires=(),
)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("features", type=Path)
    p.add_argument("--sample", required=True)
    p.add_argument("--work", type=Path, default=Path("work"))
    p.add_argument("--reps", type=Path,
                   default=Path(__file__).resolve().parent.parent / "representative_proteins.parquet")
    p.add_argument("--skip-clusters", action="store_true")
    p.add_argument("--force", action="store_true")
    a = p.parse_args()
    import pyarrow.parquet as pq

    feats = pq.read_table(a.features, columns=["feature_id"])
    r = run(feats, workdir(a.work, a.sample, "s07_match"), a.sample,
            source=a.features, reps=a.reps, skip_clusters=a.skip_clusters,
            force=a.force)
    print(r.describe())


if __name__ == "__main__":
    main()
