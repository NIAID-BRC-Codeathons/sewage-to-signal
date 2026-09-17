"""Stage 06 - SAE feature extraction.

The scalable form of the prototype in ``sae_testing_script.py``. Two changes
matter at volume:

* **Length-bucketed batching.** Sorting by length before batching keeps padding
  waste low. Note that on Apple Silicon the 6B is memory-bandwidth-bound, so
  batching buys ~20% at best (measured: 1.03 -> 1.24 seq/s); on a CUDA GPU the
  gain is the usual large one.
* **No densification.** The prototype called ``.to_dense()``, materialising an
  L x 16,384 matrix per sequence (and a second copy for normalisation). Here
  the top-K is taken straight off the sparse COO values via a scatter-reduce
  into a (batch x codebook) buffer - 2 MB at batch=32 regardless of length.

Output is long-format parquet: one row per (gene_id, feature_id, activation),
which joins directly against the SAE feature table.
"""

from __future__ import annotations

import argparse
import logging
import os
import time
from pathlib import Path

from common import StageResult, is_current, workdir, write_fragment, write_manifest
from entities import fasta_records, table_from_fasta
from stage import Column, Param, Stage

# Backbone, SAE repo and layer have to agree; picking them separately is an
# easy way to get a silently wrong answer, so they are offered as one choice.
# ESMC-6B needs ~12 GB for weights alone - on anything smaller, or on CPU, use
# 300m. Only the 6B layer-60 SAE has a published feature description table, so
# 300m gives retrieval but no human-readable summaries in s07.
#
# This lives with the stage rather than with the driver: it is a fact about
# what this stage can run, and the web UI reads it off the parameter's choices
# rather than being told separately.
MODELS: dict[str, tuple[str, str, int]] = {
    "6b": ("biohub/ESMC-6B",
           "biohub/ESMC-6B-sae-layer60-k64-codebook16384", 60),
    "300m": ("biohub/ESMC-300M",
             "biohub/ESMC-300M-sae-layer23-k64-codebook16384", 23),
}
DEFAULT_MODEL = "6b"
DEFAULT_BACKBONE, DEFAULT_SAE, DEFAULT_LAYER = MODELS[DEFAULT_MODEL]
HIT_HELP = {
    "gene_id": "gene this activation belongs to",
    "feature_id": "SAE codebook feature",
    "activation": "idf-weighted, max-normalised activation",
    "raw_activation": "mean-pooled activation before weighting",
}
FEATURE_HELP = {
    "n_genes": "genes in this sample that activated the feature",
    "max_activation": "strongest activation seen in this sample",
    "mean_activation": "mean activation across the genes that have it",
}
GENE_HELP = {
    "embedded": "the gene was sent through the model",
    "n_features": "how many features survived the top-K cut",
    "top_feature": "strongest-activating feature for this gene",
    "top_activation": "that feature's activation",
}

CODEBOOK = 16384


def load_models(backbone, sae_repo, layer, device, dtype):
    from esm.models.esmc import EsmcForMaskedLM, EsmcSaeModel, EsmcTokenizer

    tok = EsmcTokenizer.from_pretrained(backbone)
    model = EsmcForMaskedLM.from_pretrained(backbone, dtype=dtype).eval().to(device)
    sae = EsmcSaeModel.from_pretrained(sae_repo, device=device, dtype=dtype)
    sae.initialize_layers([layer], device=device, dtype=dtype)
    # Single-layer repos auto-load on CPU and initialize_layers skips anything
    # already present, so neither call reliably honours device/dtype.
    sae_layer = sae.layers[str(layer)].to(device=device, dtype=dtype)
    model.add_sae_models([sae_layer])
    return tok, model, sae_layer


def feature_stats(codebook_dim: int):
    """UniRef90 idf/max from the published table.

    The layer-60 SAE ships these buffers as all ones, making the library's own
    ``normalize_sae=True`` a no-op; without the real statistics, ranking is
    dominated by features that fire on nearly every protein.
    """
    import pyarrow.parquet as pq
    import torch
    from huggingface_hub import hf_hub_download

    path = hf_hub_download(
        "biohub/ESMC-SAE-Features", "uniref90_feature_table.parquet", repo_type="dataset"
    )
    t = pq.read_table(path, columns=["feature_id", "uniref90_idf", "uniref90_max_activation"])
    if t.num_rows != codebook_dim:
        return None
    idf, mx = torch.ones(codebook_dim), torch.ones(codebook_dim)
    ids = t.column("feature_id").to_pylist()
    idf[ids] = torch.tensor(t.column("uniref90_idf").to_pylist())
    mx[ids] = torch.tensor(t.column("uniref90_max_activation").to_pylist())
    return idf, mx.clamp(min=1e-6)


def _pool_sparse(fm, counts, batch, codebook, device):
    """Max-pool sparse SAE output per sequence, excluding <cls>/<eos>.

    `fm` rows are the batch's non-pad tokens concatenated in order, so each
    sequence owns rows [offset, offset+count).
    """
    import torch

    idx, val = fm.indices(), fm.values()
    # .contiguous(): indices() returns a non-contiguous view and bucketize
    # would otherwise copy it on every batch.
    rows, cols = idx[0].contiguous(), idx[1].contiguous()
    offsets = torch.cat([torch.zeros(1, dtype=counts.dtype, device=device),
                         counts.cumsum(0)[:-1]])
    seq_of_row = torch.bucketize(rows, offsets, right=True) - 1
    starts = offsets[seq_of_row]
    ends = starts + counts[seq_of_row] - 1
    keep = (rows != starts) & (rows != ends)          # drop <cls> and <eos>
    seq_of_row, cols, val = seq_of_row[keep], cols[keep], val[keep]

    pooled = torch.zeros(batch * codebook, device=device, dtype=torch.float32)
    pooled.scatter_reduce_(0, seq_of_row * codebook + cols, val.float(),
                           reduce="amax", include_self=True)
    return pooled.view(batch, codebook)


def run(
    rows,
    out_dir: Path,
    sample: str,
    source: Path | None = None,
    where: str | None = None,
    model: str = DEFAULT_MODEL,
    backbone: str | None = None,
    sae_repo: str | None = None,
    layer: int | None = None,
    top_k: int = 16,
    batch_size: int = 8,
    max_len: int = 1022,
    device: str = "auto",
    dtype: str = "auto",
    limit: int | None = None,
    force: bool = False,
) -> StageResult:
    # An explicit override wins over the model preset, one field at a time.
    if model not in MODELS:
        raise SystemExit(f"unknown model {model!r}; have {', '.join(MODELS)}")
    m_backbone, m_sae, m_layer = MODELS[model]
    backbone = backbone or m_backbone
    sae_repo = sae_repo or m_sae
    layer = m_layer if layer is None else layer

    out_dir = Path(out_dir)
    out = out_dir / f"{sample}.sae_features.parquet"
    feat_out = out_dir / f"{sample}.features.parquet"
    gene_out = out_dir / f"{sample}.embedded.parquet"
    params = {
        "model": model,
        "backbone": backbone, "sae_repo": sae_repo, "layer": layer,
        "top_k": top_k, "max_len": max_len, "limit": limit, "where": where,
    }
    deps = [Path(source)] if source else []
    if not force and is_current(out, deps, params):
        return StageResult("s06_embed", out, {}, skipped=True,
                           produced={"feature_hit": None, "feature": None,
                                     "gene": None})

    os.environ.setdefault("HF_HUB_DISABLE_PROGRESS_BARS", "1")
    logging.getLogger("esm").setLevel(logging.ERROR)
    import pyarrow as pa
    import pyarrow.parquet as pq
    import torch

    if device == "auto":
        device = "cuda" if torch.cuda.is_available() else (
            "mps" if torch.backends.mps.is_available() else "cpu")
    dev = torch.device(device)
    dt = getattr(torch, dtype) if dtype != "auto" else (
        torch.bfloat16 if "6B" in backbone else torch.float32)

    records = [(gid, seq[:max_len]) for gid, seq in fasta_records(rows)]
    if limit:
        records = records[:limit]
    if not records:
        raise SystemExit(
            "the selection is empty - no sequences to embed"
            + (f" (predicate: {where})" if where else ""))
    # Length-bucketed batching: sort so each batch pads to a similar length.
    records.sort(key=lambda r: len(r[1]))

    t0 = time.time()
    tok, model, sae_layer = load_models(backbone, sae_repo, layer, dev, dt)
    stats_tensors = feature_stats(CODEBOOK)
    idf, mx = (stats_tensors if stats_tensors else
               (sae_layer.idf.float().cpu(), sae_layer.max.float().cpu()))
    idf, mx = idf.to(dev), mx.to(dev)
    load_s = time.time() - t0

    gene_ids, feat_ids, acts, raws = [], [], [], []
    t1 = time.time()
    done = 0
    with torch.no_grad():
        for i in range(0, len(records), batch_size):
            chunk = records[i : i + batch_size]
            enc = tok([s for _, s in chunk], padding=True, return_tensors="pt").to(dev)
            out_b = model.esmc(**enc, compute_sae=True, normalize_sae=False)
            fm = out_b.sae_outputs[f"layer{layer}"]
            counts = enc["attention_mask"].sum(1)
            pooled_raw = _pool_sparse(fm, counts, len(chunk), CODEBOOK, dev)
            pooled_norm = pooled_raw / mx * idf
            k = min(top_k, CODEBOOK)
            top = torch.topk(pooled_norm, k, dim=1)
            ti, tv = top.indices.cpu(), top.values.cpu()
            tr = torch.gather(pooled_raw, 1, top.indices).cpu()
            for j, (gid, _) in enumerate(chunk):
                for r in range(k):
                    if tv[j, r] <= 0:
                        continue
                    gene_ids.append(gid)
                    feat_ids.append(int(ti[j, r]))
                    acts.append(float(tv[j, r]))
                    raws.append(float(tr[j, r]))
            done += len(chunk)
            if done % (batch_size * 20) == 0 or done == len(records):
                rate = done / (time.time() - t1)
                print(f"    {done}/{len(records)} proteins  {rate:.2f} seq/s", flush=True)

    table = pa.table({
        "gene_id": pa.array(gene_ids),
        "feature_id": pa.array(feat_ids, pa.int32()),
        "activation": pa.array(acts, pa.float32()),
        "raw_activation": pa.array(raws, pa.float32()),
    })
    frags = [write_fragment(out, table, "feature_hit", role="base",
                            where=where, help=HIT_HELP)]

    # The codebook feature is its own entity: s07 annotates a feature once,
    # not once per gene that activated it. Emitting the level here keeps that
    # join well-defined even when only some features were ever seen.
    per_feat: dict[int, list[float]] = {}
    per_gene: dict[str, list[tuple[int, float]]] = {g: [] for g, _ in records}
    for g, f, a in zip(gene_ids, feat_ids, acts):
        per_feat.setdefault(f, []).append(a)
        per_gene[g].append((f, a))
    fids = sorted(per_feat)
    frags.append(write_fragment(feat_out, pa.table({
        "feature_id": pa.array(fids, pa.int32()),
        "n_genes": pa.array([len(per_feat[f]) for f in fids], pa.int32()),
        "max_activation": pa.array([max(per_feat[f]) for f in fids], pa.float32()),
        "mean_activation": pa.array(
            [sum(per_feat[f]) / len(per_feat[f]) for f in fids], pa.float32()),
    }), "feature", role="base", where=where, help=FEATURE_HELP))

    # ...and the genes get told what happened to them, so a later predicate can
    # ask for what has not been embedded yet without reading the hit table.
    gids = sorted(per_gene)
    tops = [max(per_gene[g], key=lambda t: t[1], default=(None, None)) for g in gids]
    frags.append(write_fragment(gene_out, pa.table({
        "gene_id": pa.array(gids),
        "embedded": pa.array([True] * len(gids), pa.bool_()),
        "n_features": pa.array([len(per_gene[g]) for g in gids], pa.int32()),
        "top_feature": pa.array([t[0] for t in tops], pa.int32()),
        "top_activation": pa.array([t[1] for t in tops], pa.float32()),
    }), "gene", where=where, help=GENE_HELP))

    el = time.time() - t1
    stats = {
        "proteins": len(records), "rows": table.num_rows,
        "distinct_features": len(fids),
        "seq_per_s": round(len(records) / el, 2) if el else None,
        "load_s": round(load_s, 1), "device": str(dev), "dtype": str(dt),
        "idf_source": "feature_table" if stats_tensors else "sae_buffers",
    }
    write_manifest(out, deps, params, stats, seconds=el + load_s, tables=frags,
                   stage="s06_embed")
    return StageResult("s06_embed", out, stats, seconds=el,
                       produced={"feature_hit": None, "feature": None,
                                 "gene": None})


STAGE = Stage(
    name="s06_embed",
    title="SAE features",
    summary="ESMC forward pass plus the SAE. Emits feature hits, the features "
            "themselves, and columns telling each gene it was embedded.",
    run=run,
    consumes="gene",
    produces="feature_hit",
    also_produces=("feature", "gene"),
    order=60,
    selectable=True,
    roles=("projection", "embedding"),
    # Today's default: representatives that homology could not fully explain.
    # `category <> 'known'` is null-safe in the useful direction - a gene s05
    # never looked at has a null category and is not selected - and the driver
    # drops the clause entirely when no stage has written those columns.
    default_where="is_representative AND category <> 'known'",
    adds=(
        Column("feature_id", "int32", HIT_HELP["feature_id"]),
        Column("activation", "float", HIT_HELP["activation"]),
        Column("raw_activation", "float", HIT_HELP["raw_activation"]),
    ),
    params=(
        Param("model", str, "6b", choices=("6b", "300m"), group="model",
              help="6b needs ~12 GB and a GPU; 300m runs on CPU but has no "
                   "published feature description table"),
        Param("backbone", str, None, group="model",
              help="override the backbone the model choice implies"),
        Param("sae_repo", str, None, group="model",
              help="override the SAE repo the model choice implies"),
        Param("layer", int, None, group="model",
              help="override the SAE layer the model choice implies"),
        Param("top_k", int, 16, group="model",
              help="features kept per protein"),
        Param("batch_size", int, 8, group="model"),
        Param("max_len", int, 1022, group="model",
              help="sequences are truncated to this many residues"),
        Param("limit", int, None, group="model",
              help="cap how many proteins reach the GPU"),
        Param("device", str, "auto", group="model",
              choices=("auto", "cpu", "cuda", "mps")),
    ),
    requires=(),
)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("proteins", type=Path)
    p.add_argument("--sample", required=True)
    p.add_argument("--work", type=Path, default=Path("work"))
    p.add_argument("--model", choices=sorted(MODELS), default=DEFAULT_MODEL)
    p.add_argument("--backbone")
    p.add_argument("--sae-repo")
    p.add_argument("--layer", type=int)
    p.add_argument("--top-k", type=int, default=16)
    p.add_argument("--batch-size", type=int, default=8)
    p.add_argument("--device", default="auto")
    p.add_argument("--dtype", default="auto")
    p.add_argument("--limit", type=int)
    p.add_argument("--force", action="store_true")
    a = p.parse_args()
    r = run(table_from_fasta(a.proteins), workdir(a.work, a.sample, "s06_embed"), a.sample,
            source=a.proteins, model=a.model,
            backbone=a.backbone, sae_repo=a.sae_repo, layer=a.layer,
            top_k=a.top_k, batch_size=a.batch_size, device=a.device,
            dtype=a.dtype, limit=a.limit, force=a.force)
    print(r.describe())


if __name__ == "__main__":
    main()
