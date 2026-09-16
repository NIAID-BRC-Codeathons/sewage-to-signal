"""A fixed 2D projection, fitted once and reused for every run.

Projecting each run on its own makes a pretty picture and a useless one: the
layout is arbitrary, so two samples cannot be compared and re-running moves
every point. Fitting UMAP once over a reference corpus and calling
``transform`` thereafter gives every protein a stable address, which is what
the ESM Atlas does when it serves precomputed ``umap_1``/``umap_2`` columns.

This is a UMAP-specific capability. t-SNE has no ``transform`` — there is no
way to place a new point in an existing t-SNE layout without refitting.

The corpus should span what you expect to see, because ``transform`` places a
point relative to the fitted manifold and has nothing useful to say about a
region the fit never covered. So it wants both ends: proteins you care about,
and enough empirical wastewater to cover the unannotated bulk.

    python sae/web/reference_map.py build \\
        --out data/reference_map.joblib \\
        work/*/s06_embed/*.sae_features.parquet

Rebuild it when the corpus changes. Everyone sharing a map must share the file:
two people who each fit their own have incomparable coordinates, which is the
problem this exists to solve.
"""

from __future__ import annotations

import argparse
import glob
import time
from pathlib import Path

CODEBOOK = 16384
FORMAT_VERSION = 1


def sparse_features(path: Path):
    """Long-format parquet -> (gene_ids, CSR matrix over the codebook)."""
    import numpy as np
    import pyarrow.parquet as pq
    from scipy.sparse import csr_matrix

    t = pq.read_table(path, columns=["gene_id", "feature_id", "activation"])
    genes = t.column("gene_id").to_pylist()
    feats = np.asarray(t.column("feature_id").to_pylist(), dtype=np.int32)
    vals = np.asarray(t.column("activation").to_pylist(), dtype=np.float32)

    order: dict[str, int] = {}
    rows = np.empty(len(genes), dtype=np.int32)
    for i, g in enumerate(genes):
        rows[i] = order.setdefault(g, len(order))
    ids = [g for g, _ in sorted(order.items(), key=lambda kv: kv[1])]
    return ids, csr_matrix((vals, (rows, feats)), shape=(len(ids), CODEBOOK))


def normalise(m):
    """L2 by row: which features fire is the signal, not how hard."""
    from sklearn.preprocessing import normalize
    return normalize(m, norm="l2", copy=True)


def _versions() -> dict:
    import numpy, scipy, sklearn, umap
    return {"umap": umap.__version__, "numpy": numpy.__version__,
            "scipy": scipy.__version__, "sklearn": sklearn.__version__}


def build(parquets: list[Path], out: Path, seed: int = 0,
          n_neighbors: int = 15, min_dist: float = 0.1) -> dict:
    import joblib
    import numpy as np
    import umap
    from scipy.sparse import vstack

    blocks, ids, source = [], [], []
    for p in parquets:
        i, m = sparse_features(Path(p))
        blocks.append(m)
        ids.extend(i)
        source.extend([Path(p).name] * len(i))
    if not blocks:
        raise SystemExit("no input parquets")
    m = normalise(vstack(blocks))
    if m.shape[0] < 4:
        raise SystemExit(f"need at least 4 proteins, got {m.shape[0]}")

    t0 = time.time()
    reducer = umap.UMAP(n_components=2, metric="cosine", random_state=seed,
                        n_neighbors=max(2, min(n_neighbors, m.shape[0] - 1)),
                        min_dist=min_dist)
    xy = reducer.fit_transform(m)
    seconds = time.time() - t0

    counts = np.bincount(m.tocoo().col, minlength=CODEBOOK)
    payload = {
        "format_version": FORMAT_VERSION,
        "reducer": reducer,
        "ids": ids,
        "source": source,
        "xy": np.asarray(xy, dtype=np.float32),
        "inputs": [str(p) for p in parquets],
        "params": {"seed": seed, "n_neighbors": n_neighbors,
                   "min_dist": min_dist, "metric": "cosine"},
        # Pickled estimators are version-sensitive; record what built this so a
        # mismatch is reported rather than silently producing a wrong layout.
        "versions": _versions(),
        "built": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "n": m.shape[0],
        "distinct_features": int((counts > 0).sum()),
        "shared_features": int((counts > 1).sum()),
        "seconds": round(seconds, 1),
    }
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    joblib.dump(payload, out, compress=3)
    return payload


def load(path: Path) -> dict | None:
    """Load a reference map, or None if absent. Raises on a bad one."""
    path = Path(path)
    if not path.is_file():
        return None
    import joblib
    ref = joblib.load(path)
    if ref.get("format_version") != FORMAT_VERSION:
        raise ValueError(f"{path} has format_version "
                         f"{ref.get('format_version')}, expected {FORMAT_VERSION}")
    now, was = _versions(), ref.get("versions", {})
    drift = {k: (was.get(k), v) for k, v in now.items() if was.get(k) != v}
    ref["version_drift"] = drift or None
    return ref


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    b = sub.add_parser("build")
    b.add_argument("parquets", nargs="+")
    b.add_argument("--out", type=Path, default=Path("data/reference_map.joblib"))
    b.add_argument("--seed", type=int, default=0)
    b.add_argument("--n-neighbors", type=int, default=15)
    b.add_argument("--min-dist", type=float, default=0.1)
    i = sub.add_parser("info")
    i.add_argument("--out", type=Path, default=Path("data/reference_map.joblib"))
    a = p.parse_args()

    if a.cmd == "build":
        files = [Path(f) for pat in a.parquets for f in sorted(glob.glob(pat))] \
            or [Path(f) for f in a.parquets]
        print(f"  corpus: {len(files)} file(s)")
        r = build(files, a.out, a.seed, a.n_neighbors, a.min_dist)
        print(f"  fitted {r['n']} proteins in {r['seconds']}s")
        print(f"  features: {r['distinct_features']} distinct, "
              f"{r['shared_features']} shared")
        print(f"  wrote {a.out} "
              f"({Path(a.out).stat().st_size / 1048576:.1f} MB)")
    else:
        r = load(a.out)
        if r is None:
            raise SystemExit(f"no reference map at {a.out}")
        for k in ("built", "n", "distinct_features", "shared_features",
                  "params", "versions", "inputs", "version_drift"):
            print(f"  {k}: {r.get(k)}")


if __name__ == "__main__":
    main()
