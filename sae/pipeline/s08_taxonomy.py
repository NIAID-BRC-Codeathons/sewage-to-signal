"""Stage 08 - taxonomic assignment from a Kraken 2 run.

Homology triage (s05) asks what a protein *is*; this asks where it came from.
They are independent: a read can be confidently assigned to Pseudomonas and
still carry a protein no family explains, and that combination - a known
organism with an unexplained protein - is the interesting one.

Kraken classifies reads, not proteins, so the join is through the read. For
genes called on translated reads the ``contig`` column *is* the read id, which
is what this joins on; where it is null it falls back to stripping the frame
suffix a six-frame translation leaves on the peptide id.

Kraken's per-read output is enormous - 8.3 GB for 90 million reads on the run
this was written for - and all but a fraction of it describes reads that were
never embedded. It is streamed and semi-joined against the genes being
annotated, so what is read is the file and what is held is the key set. The
fifth column, the per-k-mer LCA trace, is not read: it is most of the size and
it answers how a read was classified rather than what it was classified as.

    --kraken         the per-read output (C/U, read, taxid, length, LCA)
    --kraken-report  the summary, which is where taxid gets a name and a rank
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path

import lake
from common import StageResult, workdir
from entities import table_from_fasta, write_rows
from stage import Column, Param, Stage

COLUMN_HELP = {
    "taxid": "NCBI taxon id Kraken assigned to the read (0 when unclassified)",
    "taxon": "that taxon's name, from the Kraken report",
    "taxon_rank": "its rank code - S species, G genus, F family, D domain, U unclassified",
    "taxon_classified": "Kraken classified the read at all",
}

# Kraken's report is six columns and headerless: percent, clade reads, direct
# reads, rank, taxid, name. The name is indented by tree depth.
REPORT_COLUMNS = {"pct": "DOUBLE", "clade": "BIGINT", "direct": "BIGINT",
                  "rank": "VARCHAR", "taxid": "BIGINT", "name": "VARCHAR"}
# The per-read output is five, of which the fifth is not worth parsing.
KRAKEN_COLUMNS = {"state": "VARCHAR", "read_id": "VARCHAR", "taxid": "VARCHAR",
                  "length": "VARCHAR", "lca": "VARCHAR"}


def _sql_str(v) -> str:
    return "'" + str(v).replace("'", "''") + "'"


def _cols(spec: dict) -> str:
    return "{" + ", ".join(f"'{k}':'{v}'" for k, v in spec.items()) + "}"


def run(
    rows,
    out_dir: Path,
    sample: str,
    con=None,
    source: Path | None = None,
    where: str | None = None,
    kraken: Path | None = None,
    kraken_report: Path | None = None,
    force: bool = False,
) -> StageResult:
    out_dir = Path(out_dir)
    params = {"kraken": str(kraken) if kraken else None,
              "kraken_report": str(kraken_report) if kraken_report else None,
              "where": where}
    deps = [lake.fingerprint(p) for p in (kraken, kraken_report) if p]
    if not force and lake.is_current(con, sample, "s08_taxonomy", params, deps):
        return StageResult("s08_taxonomy", out_dir, {}, skipped=True,
                           produced={"gene": None})
    if not kraken:
        raise SystemExit("s08_taxonomy needs --kraken (the per-read output)")

    import pyarrow as pa

    t0 = time.time()
    # The read behind each gene. `contig` holds it for translated reads; the
    # fallback covers a peptide id that carries its frame and nothing else.
    gene_ids = rows.column("gene_id").to_pylist()
    contigs = (rows.column("contig").to_pylist()
               if "contig" in rows.column_names else [None] * len(gene_ids))
    import re

    reads = [c if c else re.sub(r"_[0-9]+$", "", g)
             for g, c in zip(gene_ids, contigs)]
    want = pa.table({"gene_id": pa.array(gene_ids), "read_id": pa.array(reads)})

    # A throwaway connection: this reads flat files and must not hold the store
    # while it streams gigabytes.
    import duckdb

    scan = duckdb.connect(":memory:")
    scan.execute("SET preserve_insertion_order = false")
    scan.register("want", want)
    if kraken_report:
        scan.execute(f"""
            CREATE TABLE taxon AS
            SELECT taxid, trim(name) AS taxon, rank
            FROM read_csv({_sql_str(kraken_report)}, delim='\\t', header=false,
                          columns={_cols(REPORT_COLUMNS)})""")
    else:
        scan.execute("CREATE TABLE taxon (taxid BIGINT, taxon VARCHAR, rank VARCHAR)")

    joined = scan.execute(f"""
        SELECT w.gene_id,
               CAST(k.taxid AS BIGINT)        AS taxid,
               t.taxon                        AS taxon,
               t.rank                         AS taxon_rank,
               k.state = 'C'                  AS taxon_classified
        FROM want w
        JOIN (
          SELECT read_id, state, taxid FROM read_csv(
              {_sql_str(kraken)}, delim='\\t', header=false,
              columns={_cols(KRAKEN_COLUMNS)})
        ) k ON k.read_id = w.read_id
        LEFT JOIN taxon t ON t.taxid = CAST(k.taxid AS BIGINT)
    """).arrow().read_all()
    scan.close()

    if joined.num_rows == 0:
        raise SystemExit(
            f"no read in {kraken} matched any of the {len(gene_ids)} genes "
            f"selected - check that this Kraken run is for this sample")

    n = write_rows(con, "gene", sample, joined, STAGE)

    classified = sum(1 for v in joined.column("taxon_classified").to_pylist() if v)
    ranks: dict[str, int] = {}
    for r in joined.column("taxon_rank").to_pylist():
        ranks[r or "-"] = ranks.get(r or "-", 0) + 1
    stats = {
        "genes_in": len(gene_ids), "matched": joined.num_rows,
        "classified": classified,
        "classified_frac": round(classified / joined.num_rows, 4),
        "distinct_taxa": len(set(joined.column("taxid").to_pylist())),
        "species_level": ranks.get("S", 0),
    }
    el = time.time() - t0
    lake.record_run(con, sample, "s08_taxonomy", params, deps, where, stats,
                    tools={"kraken2": "external"}, seconds=el, rows=n)
    return StageResult("s08_taxonomy", out_dir, stats, seconds=el,
                       produced={"gene": None})


STAGE = Stage(
    name="s08_taxonomy",
    title="Taxonomy",
    summary="Join Kraken 2 read classifications onto genes. Independent of "
            "homology triage: a read can be confidently Pseudomonas and still "
            "carry a protein no family explains.",
    run=run,
    consumes="gene",
    produces="gene",
    order=80,
    selectable=True,
    roles=("taxonomy",),
    adds=(
        Column("taxid", "int64", COLUMN_HELP["taxid"]),
        Column("taxon", "string", COLUMN_HELP["taxon"]),
        Column("taxon_rank", "string", COLUMN_HELP["taxon_rank"]),
        Column("taxon_classified", "bool", COLUMN_HELP["taxon_classified"]),
    ),
    params=(
        Param("kraken", str, None, group="taxonomy", path=True,
              suffixes=(".tsv", ".txt", ".out"),
              help="Kraken 2 per-read output for this sample"),
        Param("kraken_report", str, None, group="taxonomy", path=True,
              suffixes=(".tsv", ".txt", ".report"),
              help="Kraken 2 report; without it taxa have ids but no names"),
    ),
    requires=(),
)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("proteins", type=Path)
    p.add_argument("--sample", required=True)
    p.add_argument("--work", type=Path, default=Path("work"))
    p.add_argument("--kraken", type=Path, required=True)
    p.add_argument("--kraken-report", type=Path)
    p.add_argument("--force", action="store_true")
    a = p.parse_args()
    with lake.open(lake.default_target(
            Path(__file__).resolve().parent.parent.parent)) as con:
        r = run(table_from_fasta(a.proteins),
                workdir(a.work, a.sample, "s08_taxonomy"), a.sample, con=con,
                kraken=a.kraken, kraken_report=a.kraken_report, force=a.force)
    print(r.describe())


if __name__ == "__main__":
    main()
