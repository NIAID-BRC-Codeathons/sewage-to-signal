#!/usr/bin/env python3
"""Map a query protein sequence onto ESM Atlas clusters via ESMC SAE features.

Pipeline
--------
1. Run the query through an ESMC backbone with a sparse-autoencoder (SAE) head
   attached to one layer, yielding a sparse per-residue activation over the SAE
   codebook.
2. Pool over residues and keep the top-K features.
3. Annotate those features from the published feature table
   (``biohub/ESMC-SAE-Features``), which describes what each of the 16,384
   features of the ESMC-6B layer-60 k64 SAE responds to.
4. Vote for candidate clusters: each top feature nominates the UniRef90
   proteins it fires hardest on (``top_100_uniref_ids``); those accessions are
   looked up in the local ``representative_proteins.parquet`` through its
   ``uniref_match_accession`` column.

On the feature -> cluster join
-----------------------------
The local tables carry no SAE-feature column, so step 4 cannot match features
against clusters directly. It bridges through UniRef accessions instead. That
bridge is sparse: only ~2.3% of the accessions named in the feature table are
present in this local atlas subset, so a query typically recovers a handful of
clusters, not an exhaustive list. Absence of a cluster here is not evidence
that it lacks the feature.

Feature ids are only meaningful for the SAE they came from. The descriptions in
``biohub/ESMC-SAE-Features`` are for the ESMC-6B layer-60 k64 codebook-16384
SAE, which is why that is the default. Pointing --backbone/--sae-repo at a
different variant leaves the ids valid but the annotations meaningless, so the
script drops annotations in that case.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent

# The (backbone, sae repo, layer) the public feature table describes.
DEFAULT_BACKBONE = "biohub/ESMC-6B"
DEFAULT_SAE_REPO = "biohub/ESMC-6B-sae-layer60-k64-codebook16384"
DEFAULT_LAYER = 60
FEATURE_TABLE_REPO = "biohub/ESMC-SAE-Features"
FEATURE_TABLE_FILE = "uniref90_feature_table.parquet"

DEMO_SEQUENCE = (
    "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQF"
)

REP_COLUMNS = [
    "protein_hash",
    "uniref_match_accession",
    "uniref_match_identity",
    "lca_taxonomy",
    "product_name",
    "cluster_top_pfam_names",
    "cluster_pct_characterized",
]


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    src = p.add_mutually_exclusive_group()
    src.add_argument("--sequence", help="query amino-acid sequence")
    src.add_argument("--fasta", type=Path, help="FASTA file; first record is used")
    p.add_argument(
        "--input-type",
        default="auto",
        choices=["auto", "protein", "nucleotide"],
        help="auto-detects nucleotide input and translates it (default auto)",
    )
    p.add_argument(
        "--orfs",
        default="longest",
        choices=["longest", "all"],
        help="for nucleotide input, analyse only the longest ORF or every ORF",
    )
    p.add_argument(
        "--min-orf-aa", type=int, default=30, help="minimum ORF length in aa (default 30)"
    )
    p.add_argument(
        "--no-gene-caller",
        action="store_true",
        help="skip pyrodigal and use plain six-frame translation",
    )
    p.add_argument("--top-k", type=int, default=5, help="features to keep (default 5)")
    p.add_argument("--backbone", default=DEFAULT_BACKBONE)
    p.add_argument("--sae-repo", default=DEFAULT_SAE_REPO)
    p.add_argument("--layer", type=int, default=DEFAULT_LAYER)
    p.add_argument(
        "--reps",
        type=Path,
        default=HERE / "representative_proteins.parquet",
        help="local representative_proteins.parquet",
    )
    p.add_argument(
        "--clusters-dir",
        type=Path,
        default=HERE / "sae_clusters",
        help="directory of cluster_members_*.parquet (used only with --member-counts)",
    )
    p.add_argument(
        "--member-counts",
        action="store_true",
        help="count cluster members; scans ~27 GB of parquet, slow",
    )
    p.add_argument(
        "--uniref-per-feature",
        type=int,
        default=100,
        help="how many of each feature's top UniRef proteins to use (default 100)",
    )
    p.add_argument("--max-clusters", type=int, default=15, help="clusters to print")
    p.add_argument("--device", default="auto", choices=["auto", "cpu", "mps", "cuda"])
    p.add_argument(
        "--dtype",
        default="auto",
        choices=["auto", "float32", "bfloat16", "float16"],
        help="auto = bfloat16 for the 6B backbone, float32 otherwise",
    )
    p.add_argument(
        "--no-idf",
        action="store_true",
        help="rank by raw activation instead of IDF-weighted (see --help notes)",
    )
    p.add_argument(
        "--offline",
        action="store_true",
        help="use only the local HF cache; skips per-run ETag revalidation",
    )
    p.add_argument(
        "--verbose",
        action="store_true",
        help="show HF progress bars and the esm fused-kernel warnings",
    )
    p.add_argument("--json-out", type=Path, help="also write results as JSON")
    p.add_argument(
        "--skip-clusters",
        action="store_true",
        help="stop after feature extraction/annotation",
    )
    return p.parse_args(argv)


def configure_environment(args) -> None:
    """Quiet the per-run cache-validation bars and unfixable kernel warnings.

    Nothing is re-downloaded between runs - the "Fetching N files" bar is
    huggingface_hub revalidating an already-populated cache (~0.25 s of ETag
    checks). --offline skips even that.

    The esm fused-kernel warnings fire once at import of esm.models.esmc, so
    the logger has to be silenced before that import happens (it lives inside
    compute_top_features). Transformer Engine and flash-attn are CUDA-only, so
    on Apple Silicon these warnings cannot be resolved by installing anything.
    """
    import logging
    import os

    if args.offline:
        # Read by huggingface_hub at import time, so set it before that import.
        os.environ["HF_HUB_OFFLINE"] = "1"
    if not args.verbose:
        os.environ.setdefault("HF_HUB_DISABLE_PROGRESS_BARS", "1")
        logging.getLogger("esm").setLevel(logging.ERROR)


def read_query(args) -> tuple[str, str]:
    """Return (name, sequence)."""
    if args.sequence:
        return "query", args.sequence.strip().upper()
    if args.fasta:
        name, chunks = None, []
        for line in args.fasta.read_text().splitlines():
            if line.startswith(">"):
                if name is not None:
                    break
                name = line[1:].strip() or "query"
            elif name is not None:
                chunks.append(line.strip())
        if not chunks:
            sys.exit(f"no sequence found in {args.fasta}")
        return name, "".join(chunks).upper()
    return "demo", DEMO_SEQUENCE


# Standard genetic code. Translation tables 1 and 11 share this codon->AA
# mapping; table 11 differs only in which codons may act as starts, which is
# the gene caller's concern, not the translator's.
CODON_TABLE = {}
_BASES = "TCAG"
_AAS = (
    "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
)
for _i, _aa in enumerate(_AAS):
    CODON_TABLE[_BASES[_i // 16] + _BASES[(_i // 4) % 4] + _BASES[_i % 4]] = _aa

_COMPLEMENT = str.maketrans("ACGTUNRYSWKMBDHVacgtunryswkmbdhv",
                            "TGCAANYRSWMKVHDBtgcaanyrswmkvhdb")
NUCLEOTIDE_CHARS = set("ACGTUN")


def looks_like_nucleotide(sequence: str) -> bool:
    """True when the sequence is overwhelmingly ACGT/U/N.

    A real protein can in principle be spelled entirely from {A,C,G,T} (Ala,
    Cys, Gly, Thr), so this is a heuristic. --input-type overrides it.
    """
    if len(sequence) < 20:
        return False
    hits = sum(1 for c in sequence if c in NUCLEOTIDE_CHARS)
    return hits / len(sequence) >= 0.9


def reverse_complement(nt: str) -> str:
    return nt.translate(_COMPLEMENT)[::-1]


def translate(nt: str) -> str:
    """Translate a nucleotide string in frame 0; unknown codons become 'X'."""
    nt = nt.upper().replace("U", "T")
    return "".join(
        CODON_TABLE.get(nt[i : i + 3], "X") for i in range(0, len(nt) - 2, 3)
    )


def six_frame_orfs(nt: str, min_aa: int):
    """ORFs from all six frames, split on stop codons.

    Start codons are not required: metagenomic contigs are routinely fragments,
    so demanding a Met start would discard genuine partial genes.
    """
    orfs = []
    for strand, seq in (("+", nt), ("-", reverse_complement(nt))):
        for frame in range(3):
            protein = translate(seq[frame:])
            offset = 0
            for segment in protein.split("*"):
                if len(segment) >= min_aa:
                    start_nt = frame + offset * 3
                    orfs.append(
                        {
                            "aa": segment,
                            "strand": strand,
                            "frame": frame + 1,
                            "nt_start": start_nt + 1,
                            "nt_end": start_nt + len(segment) * 3,
                            "caller": "six-frame",
                        }
                    )
                offset += len(segment) + 1
    return orfs


def call_genes(nt: str, args):
    """Prefer pyrodigal's metagenomic gene caller; fall back to six frames."""
    if not args.no_gene_caller:
        try:
            import pyrodigal
        except ImportError:
            print(
                "pyrodigal not installed; falling back to six-frame translation "
                "(pip install pyrodigal for proper gene calling).",
                file=sys.stderr,
            )
        else:
            # meta mode uses pre-trained profiles, which is what you want for
            # mixed-organism contigs where per-genome training is impossible.
            finder = pyrodigal.GeneFinder(meta=True)
            genes = finder.find_genes(nt.encode())
            orfs = []
            for gene in genes:
                aa = gene.translate().rstrip("*")
                if len(aa) >= args.min_orf_aa:
                    orfs.append(
                        {
                            "aa": aa,
                            "strand": "+" if gene.strand > 0 else "-",
                            "frame": None,
                            "nt_start": gene.begin,
                            "nt_end": gene.end,
                            "partial": bool(gene.partial_begin or gene.partial_end),
                            "caller": "pyrodigal",
                        }
                    )
            if orfs:
                return orfs
            print(
                "pyrodigal called no genes above the length cutoff; "
                "falling back to six-frame translation.",
                file=sys.stderr,
            )
    return six_frame_orfs(nt.upper().replace("U", "T"), args.min_orf_aa)


def prepare_queries(args, name: str, sequence: str):
    """Return [(name, protein, meta)] for the query, translating if needed."""
    if args.input_type == "protein":
        is_nt = False
    elif args.input_type == "nucleotide":
        is_nt = True
    else:
        is_nt = looks_like_nucleotide(sequence)
        if is_nt:
            print(
                f"Detected nucleotide input ({len(sequence)} nt); translating. "
                "Use --input-type protein to override.",
                file=sys.stderr,
            )

    if not is_nt:
        return [(name, sequence, {"input_type": "protein"})]

    orfs = call_genes(sequence, args)
    if not orfs:
        sys.exit(
            f"No ORF of at least {args.min_orf_aa} aa found in {len(sequence)} nt. "
            "Lower --min-orf-aa, or pass a protein sequence."
        )
    orfs.sort(key=lambda o: len(o["aa"]), reverse=True)
    called = len(orfs)
    if args.orfs == "longest":
        orfs = orfs[:1]
    print(
        f"Translated {len(sequence)} nt -> {called} ORF(s) >= {args.min_orf_aa} aa "
        f"via {orfs[0]['caller']} (longest {len(orfs[0]['aa'])} aa)"
    )
    if called > len(orfs):
        # Report the discard explicitly: on a multi-gene contig, --orfs longest
        # silently analysing 1 of N genes is easy to mistake for N being 1.
        print(
            f"  Analysing the longest only; {called - len(orfs)} other ORF(s) "
            "ignored. Use --orfs all to analyse every gene."
        )
    out = []
    for i, orf in enumerate(orfs, start=1):
        label = f"{name}|orf{i}:{orf['nt_start']}-{orf['nt_end']}({orf['strand']})"
        meta = {"input_type": "nucleotide", **{k: v for k, v in orf.items() if k != "aa"}}
        out.append((label, orf["aa"], meta))
    return out


STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")


def validate_sequence(sequence: str) -> None:
    odd = sorted(set(sequence) - STANDARD_AA)
    if odd:
        print(
            f"Warning: sequence contains non-standard characters {odd}. "
            "ESMC tokenizes '.' as a gap and '-' as an insertion, so a literal "
            "'...' is read as three gaps rather than 'and so on'.",
            file=sys.stderr,
        )


def resolve_device(choice: str):
    import torch

    if choice != "auto":
        return torch.device(choice)
    if torch.cuda.is_available():
        return torch.device("cuda")
    if torch.backends.mps.is_available():
        return torch.device("mps")
    return torch.device("cpu")


def resolve_dtype(choice: str, backbone: str):
    import torch

    if choice != "auto":
        return getattr(torch, choice)
    # The 6B weights are 25 GB in fp32; bf16 keeps it inside a 32 GB machine.
    return torch.bfloat16 if "6B" in backbone else torch.float32


def load_feature_stats(codebook_dim: int):
    """Per-feature (idf, max) from the published table, as dense tensors.

    The layer-60 SAE ships its ``idf``/``max`` buffers as all-ones placeholders,
    so the library's own ``normalize_sae=True`` is a no-op. The real UniRef90
    statistics live in the feature table instead; without them, ranking is
    dominated by high-magnitude features that fire on almost everything.
    """
    import pyarrow.parquet as pq
    import torch
    from huggingface_hub import hf_hub_download

    path = hf_hub_download(FEATURE_TABLE_REPO, FEATURE_TABLE_FILE, repo_type="dataset")
    table = pq.read_table(
        path, columns=["feature_id", "uniref90_idf", "uniref90_max_activation"]
    )
    if table.num_rows != codebook_dim:
        print(
            f"Warning: feature table has {table.num_rows} rows but the SAE "
            f"codebook is {codebook_dim}; skipping IDF weighting.",
            file=sys.stderr,
        )
        return None
    idf = torch.ones(codebook_dim)
    mx = torch.ones(codebook_dim)
    ids = table.column("feature_id").to_pylist()
    idf[ids] = torch.tensor(table.column("uniref90_idf").to_pylist())
    mx[ids] = torch.tensor(table.column("uniref90_max_activation").to_pylist())
    return idf, mx.clamp(min=1e-6)


def load_backbone(args):
    """Load backbone + SAE once, so multiple ORFs reuse one 24 GB load."""
    import torch
    from esm.models.esmc import EsmcForMaskedLM, EsmcSaeModel, EsmcTokenizer

    device = resolve_device(args.device)
    dtype = resolve_dtype(args.dtype, args.backbone)
    print(f"Loading {args.backbone} on {device} ({dtype})...", flush=True)

    tokenizer = EsmcTokenizer.from_pretrained(args.backbone)
    model = EsmcForMaskedLM.from_pretrained(args.backbone, dtype=dtype)
    model.eval().to(device)

    print(f"Loading SAE {args.sae_repo} (layer {args.layer})...", flush=True)
    sae = EsmcSaeModel.from_pretrained(args.sae_repo, device=device, dtype=dtype)
    sae.initialize_layers([args.layer], device=device, dtype=dtype)
    # A single-layer repo is auto-loaded by from_pretrained, and
    # initialize_layers skips layers that are already present - so neither call
    # is guaranteed to have honoured device/dtype. Move it explicitly.
    sae_layer = sae.layers[str(args.layer)].to(device=device, dtype=dtype)
    model.add_sae_models([sae_layer])
    return {"tokenizer": tokenizer, "model": model, "sae_layer": sae_layer,
            "device": device, "torch": torch}


def compute_top_features(args, sequence: str, bundle, stats=None):
    """Run backbone + SAE and return the top-k features for the sequence."""
    torch = bundle["torch"]
    tokenizer, model = bundle["tokenizer"], bundle["model"]
    sae_layer, device = bundle["sae_layer"], bundle["device"]

    enc = tokenizer(sequence, return_tensors="pt").to(device)
    with torch.no_grad():
        # normalize_sae=False keeps raw magnitudes, which are what the feature
        # table's `threshold` column is expressed in; the idf/max scaling used
        # for ranking is applied below from the same buffers.
        out = model.esmc(**enc, compute_sae=True, normalize_sae=False)

    raw = out.sae_outputs[f"layer{args.layer}"].to_dense().float()
    n_tokens = int(enc["input_ids"].shape[1])
    if raw.shape[0] != n_tokens:
        raise RuntimeError(f"expected {n_tokens} SAE rows, got {raw.shape[0]}")
    # Drop <cls>/<eos> so pooling sees residues only.
    raw = raw[1:-1]

    raw = raw.cpu()
    if stats is not None and not args.no_idf:
        idf, mx = stats
    else:
        # Fallback: whatever the SAE itself shipped (often all ones).
        idf = sae_layer.idf.float().cpu()
        mx = sae_layer.max.float().cpu()
    normalized = raw / mx * idf

    pooled_norm = normalized.max(dim=0).values
    pooled_raw = raw.max(dim=0).values
    top = torch.topk(pooled_norm, args.top_k)

    features = []
    for rank, (fid, score) in enumerate(
        zip(top.indices.tolist(), top.values.tolist()), start=1
    ):
        residue = int(normalized[:, fid].argmax())
        features.append(
            {
                "rank": rank,
                "feature_id": int(fid),
                "normalized_activation": float(score),
                "raw_activation": float(pooled_raw[fid]),
                # 1-based position in the input sequence
                "peak_residue": residue + 1,
                "peak_aa": sequence[residue] if residue < len(sequence) else "?",
            }
        )
    return features


def annotates_with_table(args) -> bool:
    """The published descriptions only apply to the 6B layer-60 k64 16384 SAE."""
    return args.sae_repo == DEFAULT_SAE_REPO and args.layer == DEFAULT_LAYER


def load_feature_table(feature_ids):
    import pyarrow as pa
    import pyarrow.compute as pc
    import pyarrow.parquet as pq
    from huggingface_hub import hf_hub_download

    path = hf_hub_download(
        FEATURE_TABLE_REPO, FEATURE_TABLE_FILE, repo_type="dataset"
    )
    table = pq.read_table(
        path,
        columns=[
            "feature_id",
            "summary",
            "category",
            "exemplar_protein_families",
            "threshold",
            "uniref90_frequency",
            "top_100_uniref_ids",
        ],
    )
    keep = pc.is_in(table.column("feature_id"), value_set=pa.array(sorted(feature_ids)))
    return {row["feature_id"]: row for row in table.filter(keep).to_pylist()}


def strip_uniref_prefix(accession: str) -> str:
    return accession.split("_", 1)[1] if "_" in accession else accession


def find_candidate_clusters(args, features, annotations):
    """Vote for local cluster representatives from the top features' UniRef hits."""
    import pyarrow as pa
    import pyarrow.compute as pc
    import pyarrow.dataset as ds

    # feature -> {bare accession: activation}
    nominations: dict[int, dict[str, float]] = {}
    for feat in features:
        row = annotations.get(feat["feature_id"])
        if not row:
            continue
        hits = (row.get("top_100_uniref_ids") or [])[: args.uniref_per_feature]
        nominations[feat["feature_id"]] = {
            e["uniref_id"]: e["activation"] for e in hits
        }

    wanted = {acc for m in nominations.values() for acc in m}
    if not wanted:
        return []

    # The local column keeps the cluster prefix (UniRef90_/UniRef100_); match both.
    probes = pa.array(sorted({f"{p}_{a}" for a in wanted for p in ("UniRef90", "UniRef100")}))
    print(
        f"Scanning {args.reps.name} for {len(wanted)} UniRef accessions...", flush=True
    )
    dataset = ds.dataset(args.reps, format="parquet")
    matched = dataset.to_table(
        columns=REP_COLUMNS,
        filter=ds.field("uniref_match_accession").isin(probes),
    )
    print(f"  {matched.num_rows} representative rows matched.", flush=True)

    weight = {f["feature_id"]: f["normalized_activation"] for f in features}
    clusters: dict[str, dict] = {}
    for row in matched.to_pylist():
        bare = strip_uniref_prefix(row["uniref_match_accession"] or "")
        entry = clusters.setdefault(
            row["protein_hash"],
            {
                "cluster_rep_protein_hash": row["protein_hash"],
                "uniref_match_accession": row["uniref_match_accession"],
                "uniref_match_identity": row["uniref_match_identity"],
                "lca_taxonomy": row["lca_taxonomy"],
                "product_name": row["product_name"],
                "pfam": [list(t) for t in (row["cluster_top_pfam_names"] or [])],
                "pct_characterized": row["cluster_pct_characterized"],
                "supporting_features": [],
                "score": 0.0,
            },
        )
        for fid, accs in nominations.items():
            if bare in accs:
                entry["supporting_features"].append(
                    {"feature_id": fid, "feature_activation_on_hit": accs[bare]}
                )
                entry["score"] += weight.get(fid, 0.0)

    ranked = sorted(
        clusters.values(),
        key=lambda c: (len(c["supporting_features"]), c["score"]),
        reverse=True,
    )
    return ranked


def add_member_counts(args, clusters):
    import pyarrow as pa
    import pyarrow.dataset as ds

    files = sorted(args.clusters_dir.glob("cluster_members_*.parquet"))
    if not files:
        print(f"  no cluster_members_*.parquet in {args.clusters_dir}; skipping.")
        return
    probes = pa.array([c["cluster_rep_protein_hash"] for c in clusters])
    print(f"Counting members across {len(files)} shards (slow)...", flush=True)
    dataset = ds.dataset(files, format="parquet")
    hits = dataset.to_table(
        columns=["cluster_rep_protein_hash"],
        filter=ds.field("cluster_rep_protein_hash").isin(probes),
    )
    counts = (
        hits.group_by("cluster_rep_protein_hash").aggregate([([], "count_all")])
    ).to_pylist()
    lookup = {r["cluster_rep_protein_hash"]: r["count_all"] for r in counts}
    for c in clusters:
        c["member_count"] = lookup.get(c["cluster_rep_protein_hash"], 0)


def shorten(text, limit=200):
    if not text:
        return "-"
    text = re.sub(r"\s+", " ", str(text)).strip()
    return text if len(text) <= limit else text[: limit - 1] + "…"


def main(argv=None):
    args = parse_args(argv)
    configure_environment(args)
    name, raw_sequence = read_query(args)
    if not raw_sequence:
        sys.exit("empty query sequence")
    queries = prepare_queries(args, name, raw_sequence)

    use_table = annotates_with_table(args)
    stats = None
    if use_table and not args.no_idf:
        # Codebook width is fixed for this SAE variant; assert it matches.
        stats = load_feature_stats(16384)

    bundle = load_backbone(args)
    all_results = []
    for name, sequence, meta in queries:
        print(f"\nQuery '{name}': {len(sequence)} residues")
        validate_sequence(sequence)
        all_results.append(
            analyse_one(args, name, sequence, meta, bundle, stats, use_table)
        )

    result = {"queries": all_results}
    if args.json_out:
        args.json_out.write_text(json.dumps(result, indent=2))
        print(f"\nWrote {args.json_out}")
    return result


def analyse_one(args, name, sequence, meta, bundle, stats, use_table):
    features = compute_top_features(args, sequence, bundle, stats)

    annotations = {}
    if use_table:
        annotations = load_feature_table({f["feature_id"] for f in features})
    else:
        print(
            "\nNote: --sae-repo/--layer differ from the variant the public feature\n"
            "table describes, so feature descriptions are omitted.",
            file=sys.stderr,
        )

    ranking = "raw activation" if (stats is None or args.no_idf) else "IDF-weighted"
    print(f"\n=== Top {len(features)} SAE features for '{name}' (ranked by {ranking}) ===")
    for f in features:
        row = annotations.get(f["feature_id"], {})
        thr = row.get("threshold")
        f["above_threshold"] = (
            None if thr is None else bool(f["raw_activation"] >= thr)
        )
        f["summary"] = row.get("summary")
        f["category"] = row.get("category")
        f["exemplar_protein_families"] = row.get("exemplar_protein_families")
        f["uniref90_frequency"] = row.get("uniref90_frequency")

        flag = (
            ""
            if f["above_threshold"] is None
            else ("" if f["above_threshold"] else "  [below description threshold]")
        )
        print(
            f"\n[{f['rank']}] feature {f['feature_id']}  "
            f"norm={f['normalized_activation']:.3f}  raw={f['raw_activation']:.3f}  "
            f"peak={f['peak_aa']}{f['peak_residue']}{flag}"
        )
        if row:
            print(f"     category: {row.get('category')}")
            print(f"     summary : {shorten(row.get('summary'), 300)}")

    result = {"query": name, "sequence": sequence, "source": meta,
              "features": features}

    if not args.skip_clusters and annotations:
        if not args.reps.exists():
            print(f"\n{args.reps} not found; skipping cluster lookup.", file=sys.stderr)
        else:
            print()
            clusters = find_candidate_clusters(args, features, annotations)
            if args.member_counts and clusters:
                add_member_counts(args, clusters[: args.max_clusters])
            result["candidate_clusters"] = clusters

            print(f"\n=== Candidate ESM Atlas clusters ({len(clusters)} found) ===")
            if not clusters:
                print(
                    "None. The feature->cluster bridge runs through UniRef "
                    "accessions and only ~2.3% of the feature table's proteins "
                    "are in this local subset, so this is common."
                )
            for c in clusters[: args.max_clusters]:
                fids = sorted({s["feature_id"] for s in c["supporting_features"]})
                pfam = ", ".join(f"{a} ({b})" for a, b in c["pfam"][:3]) or "-"
                members = (
                    f"  members={c['member_count']}" if "member_count" in c else ""
                )
                print(
                    f"\n  {c['cluster_rep_protein_hash']}  features={fids}  "
                    f"score={c['score']:.2f}{members}"
                )
                print(f"    product : {shorten(c['product_name'], 110)}")
                print(f"    taxonomy: {c['lca_taxonomy']}")
                print(f"    pfam    : {shorten(pfam, 110)}")
                print(
                    f"    uniref  : {c['uniref_match_accession']} "
                    f"({c['uniref_match_identity']:.1f}% id)"
                    if c["uniref_match_identity"] is not None
                    else f"    uniref  : {c['uniref_match_accession']}"
                )

    return result


if __name__ == "__main__":
    main()
