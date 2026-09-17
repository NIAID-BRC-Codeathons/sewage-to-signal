"""Entity levels: what a row is, and how a stage writes columns onto one.

An **entity level** is a kind of row - a gene, an SAE feature, a feature hit.
A stage either

* **creates** rows in a level (it consumes something else and emits these), or
* **annotates** a level - same rows, new columns.

Both are tables in the lake (see ``lake.py``); creating is an INSERT and
annotating is an UPDATE of just the columns that stage declares. Nothing here
reconstructs a table from files any more: the store holds it.

Levels are declared here rather than by the stages, because a level is a shared
vocabulary - two stages annotating ``gene`` have to agree on what a gene is.
Stages declare which level they read and write and what columns they add; the
DDL is generated from those declarations.

Selection is SQL over the level, and because every table carries ``sample``,
a cohort question is an ordinary query:

    # within a sample
    category <> 'known' AND is_representative AND aa_len > 200

    # genes here that another sample also found dark
    seq_sha1 IN (SELECT seq_sha1 FROM gene
                 WHERE sample = 'CHI-A' AND category = 'dark')

    # the question a 381-run cohort exists to ask
    seq_sha1 IN (SELECT seq_sha1 FROM gene WHERE category = 'dark'
                 GROUP BY seq_sha1 HAVING count(DISTINCT sample) >= 3)
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

from common import read_fasta
from lake import quote

# Reads and contigs stay files. A per-read table is the one place this model
# does not pay: CASPER is ~400 billion reads, and nothing downstream wants to
# predicate over them individually - s01 and s02 report summary stats instead.
# Everything from gene calling on is tabular, which is where predicates are
# actually interesting.


@dataclass(frozen=True)
class Level:
    """A kind of row. ``key`` is what makes a row unique *within a sample*."""

    name: str
    key: tuple[str, ...]
    parent: str | None = None          # level this one hangs off, if any
    title: str = ""
    summary: str = ""

    @property
    def full_key(self) -> tuple[str, ...]:
        """The key including ``sample`` - what is unique across the whole store.

        ``gene_id`` is per-assembly, so it means nothing on its own once every
        sample lives in one table. Every join and every update keys on this.
        """
        return ("sample", *self.key)


LEVELS: dict[str, Level] = {
    lv.name: lv
    for lv in (
        Level("gene", ("gene_id",), parent=None, title="Genes",
              summary="One row per called protein-coding gene."),
        Level("feature", ("feature_id",), parent=None, title="SAE features",
              summary="One row per SAE codebook feature seen in this sample."),
        Level("feature_hit", ("gene_id", "feature_id"), parent="gene",
              title="Feature hits",
              summary="One row per (gene, feature) activation."),
        Level("cluster_hit", ("feature_id", "cluster_rep_protein_hash"),
              parent="feature", title="Cluster matches",
              summary="Candidate ESM Atlas clusters nominated by a feature."),
    )
}

# File artifacts that are not tables. Stages name these as inputs/outputs the
# same way they name levels; the registry tells them apart by lookup here.
FILE_KINDS: dict[str, str] = {
    "reads": "sequencing reads (FASTQ)",
    "contigs": "assembled contigs (FASTA)",
    "proteins": "protein sequences (FASTA)",
}


def is_level(port: str) -> bool:
    return port in LEVELS


class NoSuchLevel(KeyError):
    pass


# --------------------------------------------------------------------------
# writing
# --------------------------------------------------------------------------
def write_rows(con, level: str, sample: str, table, stage) -> int:
    """Write one stage's contribution to one level, for one sample.

    A stage that *creates* rows replaces this sample's rows wholesale, so
    re-running it is idempotent rather than cumulative. A stage that
    *annotates* updates only the columns it declares, leaving every other
    stage's columns untouched - which is what lets s04, s05 and s06 all write
    to ``gene`` without coordinating.

    The caller holds the lake; see ``lake.open``. Keep the window short.
    """
    lv = LEVELS[level]
    declared = [c.name for c in stage.columns_for(level)]
    have = set(table.column_names)
    missing = [c for c in declared if c not in have]
    if missing:
        raise ValueError(
            f"{stage.name} declares {missing} on {level} but did not write them")

    con.register("_incoming", table)
    try:
        if level == stage.consumes:
            return _annotate(con, lv, sample, declared)
        return _insert(con, lv, sample, table)
    finally:
        con.unregister("_incoming")


def _insert(con, lv: Level, sample: str, table) -> int:
    cols = [c for c in table.column_names if c != "sample"]
    con.execute(f"DELETE FROM {quote(lv.name)} WHERE sample = ?", [sample])
    con.execute(
        f"INSERT INTO {quote(lv.name)} (sample, {', '.join(quote(c) for c in cols)}) "
        f"SELECT ?, {', '.join(quote(c) for c in cols)} FROM _incoming", [sample])
    return table.num_rows


def _annotate(con, lv: Level, sample: str, declared: list[str]) -> int:
    """UPDATE just this stage's columns, matched on the level's full key."""
    if not declared:
        return 0
    sets = ", ".join(f"{quote(c)} = _incoming.{quote(c)}" for c in declared)
    on = " AND ".join(f"{quote(lv.name)}.{quote(k)} = _incoming.{quote(k)}"
                      for k in lv.key)
    con.execute(
        f"UPDATE {quote(lv.name)} SET {sets} FROM _incoming "
        f"WHERE {quote(lv.name)}.sample = ? AND {on}", [sample])
    return con.execute("SELECT count(*) FROM _incoming").fetchone()[0]


# --------------------------------------------------------------------------
# reading
# --------------------------------------------------------------------------
def select(con, level: str, where: str | None = None,
           columns: Iterable[str] | None = None, limit: int | None = None,
           sample: str | None = None):
    """Rows of ``level`` as a pyarrow table.

    ``sample=None`` means the whole cohort. ``where`` is a SQL boolean
    expression over the level's columns; it is user-supplied and evaluated, so
    callers taking it from a network request must pass it through
    ``guard_predicate`` first.
    """
    cols = ", ".join(quote(c) for c in columns) if columns else "*"
    sql = f"SELECT {cols} FROM {quote(level)}"
    clauses, params = [], []
    if sample is not None:
        clauses.append("sample = ?")
        params.append(sample)
    if where:
        clauses.append(f"({where})")
    if clauses:
        sql += " WHERE " + " AND ".join(clauses)
    if limit:
        sql += f" LIMIT {int(limit)}"
    # .arrow() hands back a lazy reader over a connection the caller is about
    # to detach; read it now.
    return con.execute(sql, params).arrow().read_all()


def count(con, level: str, where: str | None = None,
          sample: str | None = None) -> int:
    sql = f"SELECT count(*) FROM {quote(level)}"
    clauses, params = [], []
    if sample is not None:
        clauses.append("sample = ?")
        params.append(sample)
    if where:
        clauses.append(f"({where})")
    if clauses:
        sql += " WHERE " + " AND ".join(clauses)
    return int(con.execute(sql, params).fetchone()[0])


def columns_of(con, level: str) -> list[dict]:
    """The level's columns, with the stage that owns each one."""
    from stage import registry

    owner = {}
    for st in registry():
        for c in st.columns_for(level):
            owner.setdefault(c.name, (st.name, c.help))
    rows = con.execute(
        "SELECT column_name, data_type FROM information_schema.columns "
        "WHERE table_name = ? ORDER BY ordinal_position", [level]).fetchall()
    out = []
    for name, typ in rows:
        stage_name, help_text = owner.get(name, ("", ""))
        out.append({"name": name, "type": typ, "stage": stage_name,
                    "help": help_text})
    return out


def levels_present(con, sample: str | None = None) -> list[str]:
    """Levels that actually hold rows, for a sample or for the whole store."""
    out = []
    for name in LEVELS:
        try:
            if count(con, name, sample=sample) > 0:
                out.append(name)
        except Exception:
            continue
    return out


def samples(con) -> list[str]:
    return [r[0] for r in con.execute(
        "SELECT DISTINCT sample FROM gene ORDER BY sample").fetchall()]


# --------------------------------------------------------------------------
# predicates
# --------------------------------------------------------------------------
# Statements that would write, attach or shell out. A predicate is an
# expression, so any of these appearing in one means it is not a predicate.
_FORBIDDEN = (
    "attach", "copy", "delete", "drop", "export", "insert",
    "install", "load", "pragma", "set ", "update", "call", "system",
    "create",
)


def guard_predicate(where: str) -> str:
    """Reject anything that is not a read-only boolean expression.

    DuckDB will happily run ``COPY ... TO`` or ``INSTALL`` from inside a WHERE
    clause via a subquery, and the web server writes files. So predicates
    arriving over HTTP are checked here and evaluated on a read-only attach.
    """
    if not where or not where.strip():
        return ""
    s = where.strip().rstrip(";")
    if ";" in s:
        raise ValueError("a predicate is one expression; ';' is not allowed")
    low = " " + s.lower().replace("(", " ( ").replace(",", " , ") + " "
    for word in _FORBIDDEN:
        if f" {word.strip()} " in low:
            raise ValueError(f"{word.strip()!r} is not allowed in a predicate")
    return s


def columns_named(where: str, known: Sequence[str]) -> set[str]:
    """Which known column names a predicate mentions.

    Used to decide whether a stage's default selection can apply at all: a
    default naming a column no stage has filled yet is dropped rather than
    silently matching nothing.
    """
    import re

    if not where:
        return set()
    words = set(re.findall(r"[A-Za-z_][A-Za-z0-9_]*", where))
    return {c for c in known if c in words}


# --------------------------------------------------------------------------
# the gene level's shared schema
# --------------------------------------------------------------------------
# Two stages can put genes on the table - calling them from contigs, or
# importing a protein FASTA somebody else produced - and everything downstream
# has to see the same columns either way. So the base schema lives with the
# level rather than with whichever stage happened to write it.
GENE_COLUMN_HELP = {
    "gene_id": "unique within the sample; from the FASTA header",
    "contig": "contig the gene was called on, when it came from an assembly",
    "begin": "1-based start on the contig",
    "end": "1-based end on the contig",
    "strand": "+ or -",
    "partial": "gene runs off the end of its contig",
    "aa_len": "length in amino acids",
    "seq_sha1": "sha1 of the sequence; the cross-sample join key",
    "seq": "amino acid sequence",
}


def gene_schema():
    import pyarrow as pa

    return pa.schema([
        ("gene_id", pa.string()), ("contig", pa.string()),
        ("begin", pa.int64()), ("end", pa.int64()), ("strand", pa.string()),
        ("partial", pa.bool_()), ("aa_len", pa.int32()),
        ("seq_sha1", pa.string()), ("seq", pa.string()),
    ])


def gene_columns():
    """The gene level's base columns, as stage declarations.

    Both writers of the gene base - calling genes from contigs, and importing a
    protein FASTA - take these from here rather than writing their own list.
    The store's DDL is generated from the declarations, so a stage that declared
    `aa_len` as a string would make it a string for everybody: the first
    declaration wins, and there is no reason for two writers of the same rows to
    have two opinions about their types.
    """
    from stage import Column

    return tuple(Column(f.name, str(f.type), GENE_COLUMN_HELP.get(f.name, ""))
                 for f in gene_schema())


def gene_row(gene_id: str, seq: str, **extra) -> dict:
    """One gene-level row, with the hash filled in."""
    import hashlib

    row = {k: None for k in GENE_COLUMN_HELP}
    row.update({
        "gene_id": gene_id, "seq": seq, "aa_len": len(seq),
        "seq_sha1": hashlib.sha1(seq.encode()).hexdigest(),
        "partial": False,
    })
    row.update(extra)
    return row


def fasta_records(table, id_col: str = "gene_id", seq_col: str = "seq"):
    """(id, sequence) pairs from a selected table, for stages that need FASTA.

    Every level-consuming stage is handed the selected rows as an Arrow table;
    the ones that actually want sequences call this. Keeping the handover
    uniform is what lets the driver treat all of them the same.
    """
    if seq_col not in table.column_names:
        raise ValueError(
            f"{seq_col!r} is not a column of this table; the level has to carry "
            f"sequences for this stage to run")
    ids = table.column(id_col).to_pylist()
    seqs = table.column(seq_col).to_pylist()
    return [(i, s) for i, s in zip(ids, seqs) if s]


def table_from_fasta(path: Path):
    """Gene-level rows straight from a FASTA.

    Only for running a stage standalone from the command line, where there is
    no store to select from. The driver never uses this.
    """
    import pyarrow as pa

    rows = [gene_row(h.split()[0], s.rstrip("*")) for h, s in read_fasta(path)]
    return pa.Table.from_pylist(rows, schema=gene_schema())
