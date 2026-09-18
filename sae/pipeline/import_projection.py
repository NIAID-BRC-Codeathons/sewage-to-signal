"""Load a fitted UMAP layout into the store as coordinates.

A projection used to be a `.joblib` beside the store: a pickled reducer plus
its coordinates, 109 MB, and invisible to anyone you handed the lake to. The
dashboard would fall back to whatever single map they happened to have.

Coordinates are a table, so they belong in the lake. Once they are there a dump
of the store is the whole dashboard - same projections in the picker, same
points in the same places - and drawing becomes a join rather than a UMAP
transform, which is the difference between a minute and no time at all.

The fitted reducer is deliberately *not* stored. It is only needed to place a
protein that was not in the fit, and a pickled estimator is version-fragile in
a way a table of numbers is not: this same file format already failed to
unpickle once, above 4096 points, on a pynndescent/numba mismatch. A sample
added after a projection was fitted simply has no coordinates in it, which the
dashboard reports rather than papers over. Refit to include it.

    python import_projection.py ../../data/reference_map_pathogens_full.joblib
    python import_projection.py ../../data/*.joblib --lake ../../data/sae.ducklake
    python import_projection.py --list
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import lake

REPO = Path(__file__).resolve().parent.parent.parent

# Points in a fit come from one of two places: the cohort being studied, or a
# reference panel anchored alongside it. The dashboard draws them differently,
# so the distinction is stored rather than re-derived later.
#
# Decided by whether the store actually has that sample, not by the name of the
# parquet it was exported from. Filenames are a convention and drift; a sample
# either is in the lake or is not.


def split_id(gene_key: str) -> tuple[str, str]:
    """`<sample>|<gene_id>` back into its parts.

    Anchor ids carry more structure (`Mpox|YP_010377176.1|OPG210`), so only the
    first separator is significant - the rest is the gene's own name.
    """
    sample, sep, rest = gene_key.partition("|")
    return (sample, rest) if sep else ("", gene_key)


def to_table(ref: dict, cohort_samples: set[str]):
    import pyarrow as pa

    ids, xy = ref["ids"], ref["xy"]
    samples, genes, roles = [], [], []
    for key in ids:
        s, g = split_id(key)
        samples.append(s)
        genes.append(g)
        roles.append("cohort" if s in cohort_samples else "anchor")
    return pa.table({
        "sample": pa.array(samples), "gene_id": pa.array(genes),
        "x": pa.array([float(p[0]) for p in xy], pa.float32()),
        "y": pa.array([float(p[1]) for p in xy], pa.float32()),
        "role": pa.array(roles),
    })


def main():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("joblib", nargs="*", type=Path,
                   help="fitted layouts to import; the file stem is the name")
    p.add_argument("--lake", default=None)
    p.add_argument("--name", help="override the name (one file only)")
    p.add_argument("--list", action="store_true",
                   help="show what the store already holds")
    p.add_argument("--drop", help="remove a projection from the store")
    a = p.parse_args()

    target = a.lake or lake.default_target(REPO)

    if a.list or a.drop:
        with lake.open(target) as con:
            from stage import registry

            from entities import LEVELS
            lake.ensure_schema(con, registry(), LEVELS)
            if a.drop:
                con.execute("DELETE FROM projection WHERE name = ?", [a.drop])
                con.execute("DELETE FROM projection_meta WHERE name = ?", [a.drop])
                print(f"dropped {a.drop}")
            for pr in lake.projections(con):
                roles = ", ".join(f"{k}={v}" for k, v in sorted(pr["by_role"].items()))
                print(f"  {pr['name']:<32}{pr['stored']:>8} points  {roles}")
                print(f"  {'':<32}built {pr['built']}")
        return

    if not a.joblib:
        p.error("give at least one .joblib, or --list")
    if a.name and len(a.joblib) > 1:
        p.error("--name takes one file")

    import joblib as jl

    for path in a.joblib:
        if not path.is_file():
            print(f"  {path}: no such file", file=sys.stderr)
            continue
        ref = jl.load(path)
        with lake.read(target) as con:
            cohort_samples = {r[0] for r in con.execute(
                "SELECT DISTINCT sample FROM gene").fetchall()}
        rows = to_table(ref, cohort_samples)
        name = a.name or path.stem
        matched = sum(1 for r in rows.column("role").to_pylist() if r == "cohort")
        if matched == 0:
            print(f"  {path.name}: none of its points name a sample this store "
                  f"has - the layout was fitted against a different cohort. "
                  f"Importing anyway; the dashboard will fall back to "
                  f"transforming.", file=sys.stderr)
        with lake.open(target) as con:
            from stage import registry

            from entities import LEVELS
            lake.ensure_schema(con, registry(), LEVELS)
            n = lake.store_projection(
                con, name, rows, ref.get("n") or rows.num_rows,
                ref.get("params") or {},
                sorted({Path(i).name for i in (ref.get("inputs") or [])}),
                ref.get("note") or "")
        anchors = sum(1 for r in rows.column("role").to_pylist() if r == "anchor")
        print(f"  {name:<32}{n:>8} points  ({n - anchors} cohort, {anchors} anchors)")
    print(f"\n-> {target}")


if __name__ == "__main__":
    main()
