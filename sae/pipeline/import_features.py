"""Load an SAE feature parquet produced somewhere else into the store.

s06 writes its features straight to the lake, but a run done on a cluster
arrives as a file. This takes that file and makes it a sample: the feature hits
themselves, the feature level they imply, and the gene rows they hang off.

Those gene rows carry no sequence. The feature parquet records which features
fired, not what the protein was, so `seq` and `aa_len` are null and anything
needing them - dereplication, homology triage, re-embedding, the atlas's FASTA
download - cannot run for this sample until the peptides are imported too.
Everything else works: the points are on the map, and any annotation keyed on
gene_id still lands.

    python import_features.py SRR35987665_1M.sae_features.parquet \\
        --sample SRR35987665_1M

The read id is taken to be the peptide id minus its frame suffix, which is what
a six-frame translation leaves behind, and is stored as `contig` - the same
place the in-pipeline translate route puts it, so a read-keyed annotation
(s08_taxonomy) joins the same way for both.
"""

from __future__ import annotations

import argparse
import re
import sys
import time
from pathlib import Path

import lake
from entities import LEVELS, gene_schema, write_rows
from stage import registry

REPO = Path(__file__).resolve().parent.parent.parent
FRAME = re.compile(r"_[0-9]+$")


def main():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("parquet", type=Path,
                   help="long-format (gene_id, feature_id, activation, ...)")
    p.add_argument("--sample", required=True)
    p.add_argument("--lake", default=None)
    p.add_argument("--read-from-contig", action="store_true", default=True,
                   help="store the frame-stripped peptide id as `contig`")
    a = p.parse_args()

    if not a.parquet.is_file():
        sys.exit(f"no such file: {a.parquet}")
    target = a.lake or lake.default_target(REPO)
    reg = registry()

    import duckdb
    import pyarrow as pa

    t0 = time.time()
    scan = duckdb.connect(":memory:")
    scan.execute("SET preserve_insertion_order = false")
    src = "'" + str(a.parquet).replace("'", "''") + "'"

    hits = scan.execute(f"""
        SELECT gene_id, CAST(feature_id AS INTEGER) AS feature_id,
               CAST(activation AS FLOAT) AS activation,
               CAST(raw_activation AS FLOAT) AS raw_activation
        FROM read_parquet({src})""").arrow().read_all()
    print(f"  {hits.num_rows:,} feature hits", flush=True)

    genes = scan.execute(f"""
        SELECT DISTINCT gene_id FROM read_parquet({src}) ORDER BY gene_id
    """).arrow().read_all().column("gene_id").to_pylist()
    print(f"  {len(genes):,} peptides", flush=True)

    # The feature level: one row per codebook feature seen, as s06 would write.
    feats = scan.execute(f"""
        SELECT CAST(feature_id AS INTEGER)       AS feature_id,
               CAST(count(*) AS INTEGER)         AS n_genes,
               CAST(max(activation) AS FLOAT)    AS max_activation,
               CAST(avg(activation) AS FLOAT)    AS mean_activation
        FROM read_parquet({src}) GROUP BY feature_id ORDER BY feature_id
    """).arrow().read_all()

    # And what each gene got out of it, the columns s06 writes back onto genes.
    per_gene = scan.execute(f"""
        SELECT gene_id,
               CAST(count(*) AS INTEGER) AS n_features,
               CAST(arg_max(feature_id, activation) AS INTEGER) AS top_feature,
               CAST(max(activation) AS FLOAT) AS top_activation
        FROM read_parquet({src}) GROUP BY gene_id ORDER BY gene_id
    """).arrow().read_all()
    scan.close()

    rows = [{**{k: None for k in gene_schema().names},
             "gene_id": g, "contig": FRAME.sub("", g), "partial": False}
            for g in genes]
    base = pa.Table.from_pylist(rows, schema=gene_schema())

    with lake.open(target) as con:
        lake.ensure_schema(con, reg, LEVELS)
        write_rows(con, "gene", a.sample, base, reg["s03_genes"])
        write_rows(con, "feature_hit", a.sample, hits, reg["s06_embed"])
        write_rows(con, "feature", a.sample, feats, reg["s06_embed"])
        embedded = per_gene.append_column(
            "embedded", pa.array([True] * per_gene.num_rows, pa.bool_()))
        write_rows(con, "gene", a.sample, embedded, reg["s06_embed"])
        lake.register_sample(con, a.sample, "imported", str(a.parquet))
        lake.record_run(
            con, a.sample, "s06_embed",
            {"imported_from": str(a.parquet), "note": "run elsewhere"},
            [lake.fingerprint(a.parquet)], None,
            {"proteins": len(genes), "rows": hits.num_rows,
             "distinct_features": feats.num_rows},
            seconds=time.time() - t0, rows=hits.num_rows)

    print(f"\n  {a.sample}: {len(genes):,} genes, {hits.num_rows:,} hits, "
          f"{feats.num_rows:,} features -> {target}")
    print("  note: no sequences - s04/s05/s06 and the FASTA download cannot "
          "run for this sample until the peptides are imported")


if __name__ == "__main__":
    main()
