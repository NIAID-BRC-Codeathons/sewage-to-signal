"""Entity levels, column fragments, and predicates over them.

The pipeline used to pass files: each stage read a FASTA and wrote a smaller
FASTA. That works, but it hides what most stages are actually doing. ``s05``
does not transform proteins - it *labels* them, and then writes four FASTAs
purely to hand the next stage a subset. ``s04`` is the same: its real output is
a ``gene -> representative`` mapping, and ``nr.faa`` exists only to be somebody
else's input.

So the interface is a table, not a file. An **entity level** is a kind of row
(a gene, an SAE feature, a feature hit). A stage either

* **annotates** a level - same rows in, new *columns* out; or
* **emits** a level - new rows, carrying a parent id back to the level above.

Each stage writes its columns as its own parquet **fragment**, never rewriting
anybody else's. The logical table for a level is the base fragment LEFT JOINed
to every annotation fragment, assembled on demand. That keeps the existing
one-directory-per-stage layout, keeps each stage independently re-runnable, and
means a stage can be added without touching any reader.

Reading is SQL, over DuckDB. Every sample in every work root is registered as a
schema, so a predicate can reach across samples:

    # within a sample
    category <> 'known' AND is_representative AND aa_len > 200

    # across samples - genes here that were also dark in CHI-A
    seq_sha1 IN (SELECT seq_sha1 FROM "CHI-A".gene WHERE category = 'dark')

Levels are declared here rather than by the stages, because a level is a shared
vocabulary: two stages annotating ``gene`` have to agree on what a gene is.
Stages declare which level they read and write, and what columns they add.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Sequence

from common import read_fasta

# Reads and contigs stay files. A per-read table is the one place this model
# does not pay: CASPER is ~400 billion reads, and nothing downstream wants to
# predicate over them individually - s01 and s02 report summary stats instead.
# Everything from gene calling on is tabular, which is where predicates are
# actually interesting.


@dataclass(frozen=True)
class Level:
    """A kind of row. ``key`` is what makes a row unique within a sample."""

    name: str
    key: tuple[str, ...]
    parent: str | None = None          # level this one hangs off, if any
    title: str = ""
    summary: str = ""

    @property
    def key_sql(self) -> str:
        return ", ".join(f'"{k}"' for k in self.key)


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


# --------------------------------------------------------------------------
# fragments
# --------------------------------------------------------------------------
@dataclass
class Fragment:
    """One stage's contribution of columns to one level."""

    path: Path
    level: str
    key: tuple[str, ...]
    role: str                      # "base" (defines rows) or "annotation"
    columns: list[dict]            # [{"name", "type", "help"}]
    stage: str = ""
    written: str | None = None
    where: str | None = None       # predicate the producing stage selected on

    @property
    def added(self) -> list[str]:
        """Columns this fragment contributes, excluding the join key."""
        return [c["name"] for c in self.columns if c["name"] not in self.key]

    def to_json(self) -> dict:
        return {
            "path": str(self.path), "level": self.level, "key": list(self.key),
            "role": self.role, "columns": self.columns, "stage": self.stage,
            "written": self.written, "where": self.where,
        }


def describe_table(table, level: str, role: str, *, where: str | None = None,
                   help: dict[str, str] | None = None) -> dict:
    """Manifest entry for a pyarrow table a stage just wrote.

    Types come from the table itself rather than a declaration, so a stage
    cannot drift from what it actually produced.
    """
    lv = LEVELS[level]
    help = help or {}
    return {
        "level": level,
        "key": list(lv.key),
        "role": role,
        "where": where,
        "columns": [{"name": f.name, "type": str(f.type), "help": help.get(f.name, "")}
                    for f in table.schema],
    }


def fragments(sample_dir: Path, level: str | None = None) -> list[Fragment]:
    """Every fragment under a sample, discovered from stage manifests.

    Nothing here knows the stage list: it reads whatever manifests are on disk,
    so a work directory written by a pipeline this process has never heard of
    still describes itself completely.
    """
    out: list[Fragment] = []
    if not Path(sample_dir).is_dir():
        return out
    for stage_dir in sorted(p for p in Path(sample_dir).iterdir() if p.is_dir()):
        for mf in sorted(stage_dir.glob("*.manifest.json")):
            try:
                doc = json.loads(mf.read_text())
            except (json.JSONDecodeError, OSError):
                continue
            for t in doc.get("tables") or []:
                if level and t.get("level") != level:
                    continue
                p = Path(t["path"]) if t.get("path") else \
                    stage_dir / str(mf.name).replace(".manifest.json", "")
                if not p.is_absolute():
                    p = stage_dir / p
                if not p.is_file():
                    continue
                out.append(Fragment(
                    path=p, level=t["level"], key=tuple(t.get("key") or ()),
                    role=t.get("role", "annotation"), columns=t.get("columns") or [],
                    stage=doc.get("stage") or stage_dir.name,
                    written=doc.get("written"), where=t.get("where"),
                ))
    out.sort(key=lambda f: (f.role != "base", f.written or "", f.stage))
    return out


def levels_present(sample_dir: Path) -> list[str]:
    seen = []
    for f in fragments(sample_dir):
        if f.level not in seen:
            seen.append(f.level)
    return seen


# --------------------------------------------------------------------------
# assembling a level into one queryable view
# --------------------------------------------------------------------------
@dataclass
class Assembled:
    """A level's fragments resolved into a single SELECT."""

    level: str
    sql: str
    columns: list[dict]            # name, type, stage, help - provenance per column
    collisions: list[str] = field(default_factory=list)
    base: Fragment | None = None
    parts: list[Fragment] = field(default_factory=list)


def _q(p: Path) -> str:
    return "read_parquet('" + str(p).replace("'", "''") + "')"


def assemble(sample_dir: Path, level: str) -> Assembled | None:
    """Base LEFT JOIN every annotation fragment, in write order.

    A column claimed by two stages is kept from the first writer and reported
    as a collision rather than silently shadowed - the alternative is a
    predicate that quietly means something other than it says.
    """
    parts = fragments(sample_dir, level)
    if not parts:
        return None
    lv = LEVELS.get(level) or Level(level, tuple(parts[0].key))
    base = next((f for f in parts if f.role == "base"), parts[0])
    rest = [f for f in parts if f is not base]

    taken = {k: base.stage for k in lv.key}
    cols: list[dict] = []
    for c in base.columns:
        taken.setdefault(c["name"], base.stage)
        cols.append({**c, "stage": base.stage})

    select = ["b.*"]
    collisions: list[str] = []
    joins: list[str] = []
    for i, f in enumerate(rest):
        alias = f"a{i}"
        joins.append(
            f'LEFT JOIN {_q(f.path)} AS {alias} ON '
            + " AND ".join(f'b."{k}" = {alias}."{k}"' for k in lv.key)
        )
        for c in f.columns:
            n = c["name"]
            if n in lv.key:
                continue
            if n in taken:
                collisions.append(f"{n} (kept from {taken[n]}, also written by {f.stage})")
                continue
            taken[n] = f.stage
            select.append(f'{alias}."{n}"')
            cols.append({**c, "stage": f.stage})

    sql = (f"SELECT {', '.join(select)} FROM {_q(base.path)} AS b "
           + " ".join(joins)).strip()
    return Assembled(level=level, sql=sql, columns=cols, collisions=collisions,
                     base=base, parts=parts)


# --------------------------------------------------------------------------
# querying
# --------------------------------------------------------------------------
class NoSuchLevel(KeyError):
    pass


def connect(sample_dir: Path, roots: Sequence[Path] = ()) -> Any:
    """A DuckDB connection with this sample's levels as plain views.

    Sibling samples are registered as schemas, so a predicate can reference
    ``"CHI-A".gene`` and compare one sample against another. Only samples that
    actually have fragments are registered, and only as views over parquet -
    nothing is copied.
    """
    import duckdb

    # In-memory and thrown away per query. Every view is a read over parquet,
    # so nothing a predicate does can outlive the call - which is most of why
    # accepting one over HTTP is tolerable. The rest is ``guard_predicate``.
    con = duckdb.connect(":memory:")
    sample_dir = Path(sample_dir)

    for lvl in levels_present(sample_dir):
        a = assemble(sample_dir, lvl)
        if a:
            con.execute(f'CREATE OR REPLACE VIEW "{lvl}" AS {a.sql}')

    seen: set[str] = set()
    for root in roots:
        if not Path(root).is_dir():
            continue
        for other in sorted(p for p in Path(root).iterdir() if p.is_dir()):
            if other.name.startswith(".") or other.name in seen:
                continue
            present = levels_present(other)
            if not present:
                continue
            seen.add(other.name)
            con.execute(f'CREATE SCHEMA IF NOT EXISTS "{other.name}"')
            for lvl in present:
                a = assemble(other, lvl)
                if a:
                    con.execute(
                        f'CREATE OR REPLACE VIEW "{other.name}"."{lvl}" AS {a.sql}')
    return con


def select(sample_dir: Path, level: str, where: str | None = None,
           columns: Iterable[str] | None = None, limit: int | None = None,
           roots: Sequence[Path] = ()):
    """Rows of ``level`` matching ``where``, as a pyarrow table.

    ``where`` is a SQL boolean expression over the level's columns. It is
    user-supplied and evaluated, so callers that accept it from a network
    request must be read-only - see ``guard_predicate``.
    """
    con = connect(sample_dir, roots)
    try:
        if level not in [r[0] for r in con.execute(
                "SELECT table_name FROM information_schema.tables "
                "WHERE table_schema = 'main'").fetchall()]:
            raise NoSuchLevel(f"sample has no {level!r} table yet")
        cols = ", ".join(f'"{c}"' for c in columns) if columns else "*"
        sql = f'SELECT {cols} FROM "{level}"'
        if where:
            sql += f" WHERE ({where})"
        if limit:
            sql += f" LIMIT {int(limit)}"
        # .arrow() hands back a RecordBatchReader, which is lazy over a
        # connection this function is about to close. Read it now.
        return con.execute(sql).arrow().read_all()
    finally:
        con.close()


def count(sample_dir: Path, level: str, where: str | None = None,
          roots: Sequence[Path] = ()) -> int:
    con = connect(sample_dir, roots)
    try:
        sql = f'SELECT count(*) FROM "{level}"'
        if where:
            sql += f" WHERE ({where})"
        return int(con.execute(sql).fetchone()[0])
    finally:
        con.close()


# Statements that would write, attach or shell out. A predicate is an
# expression, so any of these appearing in one means it is not a predicate.
_FORBIDDEN = (
    "attach", "copy", "create", "delete", "drop", "export", "insert",
    "install", "load", "pragma", "set ", "update", "call", "system",
)


def guard_predicate(where: str) -> str:
    """Reject anything that is not a read-only boolean expression.

    DuckDB will happily run ``COPY ... TO`` or ``INSTALL`` from inside a WHERE
    clause via a subquery, and this server writes files. So predicates arriving
    over HTTP are checked here and the connection that runs them is a
    throwaway in-memory one over read-only parquet views.
    """
    if not where or not where.strip():
        return ""
    s = where.strip().rstrip(";")
    if ";" in s:
        raise ValueError("a predicate is one expression; ';' is not allowed")
    low = " " + s.lower().replace("(", " ( ").replace(",", " , ") + " "
    for word in _FORBIDDEN:
        if f" {word.strip()} " in low or low.startswith(f" {word.strip()} "):
            raise ValueError(f"{word.strip()!r} is not allowed in a predicate")
    return s


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
            f"{seq_col!r} is not a column of this table; the level's base "
            f"fragment has to carry sequences for this stage to run")
    ids = table.column(id_col).to_pylist()
    seqs = table.column(seq_col).to_pylist()
    return [(i, s) for i, s in zip(ids, seqs) if s]


def table_from_fasta(path: Path):
    """Gene-level rows straight from a FASTA.

    Only for running a stage standalone from the command line, where there is
    no work directory to select from. The driver never uses this.
    """
    import pyarrow as pa

    rows = [gene_row(h.split()[0], s.rstrip("*")) for h, s in read_fasta(path)]
    return pa.Table.from_pylist(rows, schema=gene_schema())
