"""The store. One DuckLake database holds every sample's rows.

The pipeline used to write a parquet file per stage per sample and reconstruct
tables from them by parsing manifest sidecars. That was multi-writer safe, which
the batch path needs, but it was not a database: nothing to point at, no schema,
no integrity between levels, and a cohort question like "dark in three or more
samples" could not be written down at all.

Here there is one lake. ``gene``, ``feature_hit``, ``feature`` and
``cluster_hit`` are real tables with real types and a ``sample`` column, so a
cohort query is ordinary SQL. ``stage_run`` records what produced what, which is
the provenance the manifests used to carry.

The schema is generated from the stage registry rather than written out here:
``Stage.adds`` already declares which columns a stage owns and ``LEVELS``
declares each level's key, so a new stage's columns appear by ``ADD COLUMN``
with nothing to migrate.

Holding the lake
----------------
DuckLake keeps its metadata in a DuckDB file, and DuckDB is one writer *or* many
readers - with the lock held for as long as the ``ATTACH`` lives, not just the
transaction. Measured: eight concurrent writers, seven fail at ``ATTACH`` before
reaching a commit. Retrying the *commit* never helps, because the connection
never opened.

So the rule is that nobody holds the lake across compute. Attach, write, detach;
retry the attach. Measured that way, eight writers doing an insert and two
updates each finish in 0.7 s with 54 attach retries and no failures, and a
reader polling throughout gets in every time (median 3 ms). Real stages spend
minutes on a GPU between writes, so contention is far lower in practice.

That rule is the one thing that breaks everything if broken, so this module only
hands out a context manager. There is no ``connect()`` to keep in a variable.

    with lake.open(path) as con:          # writers: brief
        con.execute("INSERT INTO gene ...")

    with lake.read(path) as con:          # readers: brief, read-only
        con.execute("SELECT ...").arrow()

Where the lake lives
--------------------
The catalog and the parquet are chosen separately, so the same pipeline runs
against a file on a laptop or a shared lake on object storage:

    SAE_LAKE=data/sae.ducklake                       # default: local, local
    SAE_LAKE=postgres:dbname=sae host=db.internal \
    SAE_LAKE_DATA=s3://sae-cohort/lake/              # shared catalog + S3

A local catalog is a DuckDB file, so it needs working POSIX locks: keep it off
Lustre/NFS and let only the parquet live on scratch. A Postgres catalog removes
both that constraint and the reader/writer exclusion above, which is the right
answer once the cohort is written from more than one machine.

Credentials are DuckDB's own - ``CREATE SECRET`` or the usual AWS_*/GOOGLE_*/
AZURE_* environment variables. Nothing here handles a secret, so nothing here
can leak one.
"""

from __future__ import annotations

import os
import random
import time
from contextlib import contextmanager
from pathlib import Path
from typing import Iterator

# Long enough to outlast a burst of writers, short enough that a genuinely stuck
# lock surfaces as an error rather than a hang.
ATTACH_BUDGET = float(os.environ.get("SAE_LAKE_ATTACH_BUDGET", "60"))
EXTENSION = "ducklake"

# Where DuckDB keeps its extensions. Unset means DuckDB's own default,
# ~/.duckdb, which is right on a workstation and wrong everywhere else: the
# container image is read-only and a compute node usually has neither a
# writable HOME nor a route to the extension repository. The image bakes its
# extensions in at build time and points this at them.
EXTENSION_DIR = os.environ.get("SAE_DUCKDB_EXTENSIONS")

# DuckLake separates the *catalog* (metadata: which files hold which snapshot of
# which table) from the *data path* (the parquet itself). Both can be remote,
# and they are chosen independently:
#
#   catalog     a local file, or postgres:/sqlite:/mysql: for a shared one
#   data path   a local directory, or s3:// gs:// az:// r2:// for object storage
#
# A local catalog is a DuckDB file, so it is one writer or many readers and the
# lock lives as long as the ATTACH - which is why everything here is a context
# manager and why attaching retries. A Postgres catalog has real MVCC and lifts
# that restriction entirely; the retry then costs nothing and never fires.
#
# Credentials are DuckDB's, not ours: run CREATE SECRET, or set the usual
# AWS_*/GOOGLE_*/AZURE_* environment variables. Nothing here handles a secret,
# so nothing here can leak one.
CATALOG_SCHEMES = ("postgres:", "sqlite:", "mysql:", "md:", "motherduck:")
REMOTE_DATA_SCHEMES = ("s3://", "gs://", "gcs://", "az://", "azure://", "r2://",
                       "http://", "https://")


def default_target(repo: Path) -> str:
    """The lake this checkout uses unless told otherwise."""
    return os.environ.get("SAE_LAKE") or str(repo / "data" / "sae.ducklake")


def is_remote_catalog(target: str) -> bool:
    return str(target).lower().startswith(CATALOG_SCHEMES)


def is_remote_data(path: str) -> bool:
    return str(path).lower().startswith(REMOTE_DATA_SCHEMES)


def attach_target(target) -> str:
    """Normalise a catalog into the string DuckLake attaches.

    A bare path means a local DuckDB catalog. Anything with a scheme is passed
    through, so `postgres:dbname=sae host=db.internal` and `sqlite:/srv/sae.db`
    work without this module knowing what they are.
    """
    t = str(target)
    if t.startswith("ducklake:"):
        return t
    return "ducklake:" + t


def data_path(target) -> str:
    """Where the lake keeps its parquet. Beside a local catalog by default.

    A remote catalog has no "beside", so pointing one at object storage means
    setting SAE_LAKE_DATA as well.
    """
    env = os.environ.get("SAE_LAKE_DATA")
    if env:
        return env if is_remote_data(env) else str(Path(env))
    t = str(target)
    if is_remote_catalog(t) or t.startswith("ducklake:"):
        raise ValueError(
            "SAE_LAKE_DATA must be set when the catalog is not a local file: "
            "a remote catalog has no directory to put parquet beside")
    return str(Path(t)) + ".files"


class LakeBusy(RuntimeError):
    """Someone else held the lake for longer than we were willing to wait."""


def _base(target: str, files: str):
    import duckdb

    con = duckdb.connect(":memory:")
    if EXTENSION_DIR:
        con.execute(f"SET extension_directory = {_sql_str(EXTENSION_DIR)}")
    for ext in _needed_extensions(target, files):
        _ensure_extension(con, ext)
    return con


def _ensure_extension(con, name: str) -> None:
    """Load an extension, installing it only if it is not already present.

    An unconditional INSTALL writes to the extension directory and reaches the
    network, and fails on both counts inside a read-only image or on an offline
    node - with an error about creating ~/.duckdb that says nothing about
    either. So try the load first; install only where installing can work.
    """
    try:
        con.execute(f"LOAD {name}")
        return
    except Exception as exc:
        first = exc
    try:
        con.execute(f"INSTALL {name}")
        con.execute(f"LOAD {name}")
    except Exception as exc:
        raise RuntimeError(
            f"the {name!r} DuckDB extension is neither installed nor "
            f"installable here"
            + (f" (looked in {EXTENSION_DIR})" if EXTENSION_DIR else "")
            + f". Load said: {first}. Install said: {exc}") from exc


def _needed_extensions(target: str, files: str) -> list[str]:
    """DuckLake plus whatever the catalog and the data path happen to need."""
    ext = [EXTENSION]
    t = str(target).lower()
    if "postgres:" in t:
        ext.append("postgres")
    elif "sqlite:" in t:
        ext.append("sqlite")
    elif "mysql:" in t:
        ext.append("mysql")
    if is_remote_data(files):
        ext.append("httpfs")
    return ext


def _sql_str(v) -> str:
    return "'" + str(v).replace("'", "''") + "'"


@contextmanager
def _attached(target, read_only: bool, budget: float) -> Iterator:
    """Attach with retry, yield, always detach.

    The retry is on the ATTACH, not on the commit: with a local catalog the
    lock is taken when the connection opens, so a writer that is already in has
    to finish and leave before anyone else can start. With a Postgres catalog
    there is nothing to wait for and the loop succeeds first time.
    """
    files = data_path(target)
    if not is_remote_catalog(str(target)) and not str(target).startswith("ducklake:"):
        Path(target).parent.mkdir(parents=True, exist_ok=True)
    if not is_remote_data(files):
        Path(files).mkdir(parents=True, exist_ok=True)
        files = str(files) + os.sep
    elif not files.endswith("/"):
        files += "/"

    con = _base(target, files)
    opts = f"DATA_PATH {_sql_str(files)}"
    if read_only:
        opts += ", READ_ONLY"
    stmt = f"ATTACH {_sql_str(attach_target(target))} AS lake ({opts})"

    deadline, attempt, last = time.monotonic() + budget, 0, None
    while True:
        try:
            con.execute(stmt)
            break
        except Exception as exc:                  # the lock is held; wait and retry
            last = exc
            if time.monotonic() >= deadline:
                con.close()
                raise LakeBusy(
                    f"could not attach {target} within {budget:g}s - another "
                    f"process is holding it. Last error: {exc}") from exc
            # Full jitter, so a batch of writers does not retry in lockstep.
            time.sleep(random.uniform(0.005, min(0.25, 0.005 * 2 ** attempt)))
            attempt += 1
    try:
        con.execute("USE lake")
        yield con
    finally:
        try:
            con.execute("DETACH lake")
        except Exception:
            pass
        con.close()


@contextmanager
def open(target, budget: float = ATTACH_BUDGET) -> Iterator:
    """Writable lake, held for as short a time as the caller can manage."""
    with _attached(target, read_only=False, budget=budget) as con:
        yield con


@contextmanager
def read(target, budget: float = ATTACH_BUDGET) -> Iterator:
    """Read-only lake. Exclusive against a writer when the catalog is a local
    file, concurrent when it is Postgres. Keep it brief either way."""
    with _attached(target, read_only=True, budget=budget) as con:
        yield con


def exists(target) -> bool:
    """Whether the lake is already there. Only answerable for a local catalog."""
    if is_remote_catalog(str(target)):
        return True            # ask the server, not the filesystem
    return Path(str(target)).is_file()


# --------------------------------------------------------------------------
# schema, derived from the stage registry
# --------------------------------------------------------------------------
# Types a stage declares in `Column.type` mapped to SQL. Stages describe their
# columns with arrow-ish names because that is what they write.
_SQL_TYPE = {
    "string": "VARCHAR", "str": "VARCHAR", "utf8": "VARCHAR", "large_string": "VARCHAR",
    "bool": "BOOLEAN", "boolean": "BOOLEAN",
    "int8": "TINYINT", "int16": "SMALLINT", "int32": "INTEGER", "int64": "BIGINT",
    "uint8": "UTINYINT", "uint16": "USMALLINT", "uint32": "UINTEGER", "uint64": "UBIGINT",
    "float": "FLOAT", "float32": "FLOAT", "double": "DOUBLE", "float64": "DOUBLE",
}

# Column names the store owns. A stage may not declare these.
RESERVED = {"sample"}


def sql_type(declared: str) -> str:
    return _SQL_TYPE.get(str(declared).lower(), "VARCHAR")


def quote(name: str) -> str:
    return '"' + str(name).replace('"', '""') + '"'


def ensure_schema(con, registry, levels) -> dict[str, list[str]]:
    """Create or extend every table the registry implies. Idempotent.

    Returns level -> column names, so a caller can see what the store actually
    holds rather than what it hoped for.
    """
    con.execute("""
        CREATE TABLE IF NOT EXISTS sample (
            sample VARCHAR NOT NULL,
            created TIMESTAMP,
            source_kind VARCHAR,
            source_path VARCHAR
        )""")
    # Provenance. Replaces <output>.manifest.json: same fields, one table.
    con.execute("""
        CREATE TABLE IF NOT EXISTS stage_run (
            sample VARCHAR NOT NULL,
            stage VARCHAR NOT NULL,
            params_hash VARCHAR,
            params VARCHAR,
            inputs VARCHAR,
            where_clause VARCHAR,
            stats VARCHAR,
            tools VARCHAR,
            seconds DOUBLE,
            rows_written BIGINT,
            snapshot_id BIGINT,
            written TIMESTAMP
        )""")

    out: dict[str, list[str]] = {}
    for name, lv in levels.items():
        cols: dict[str, str] = {"sample": "VARCHAR"}
        for k in lv.key:
            cols[k] = "VARCHAR"
        for st in registry:
            # `columns_for` is per level, because a stage that fills several
            # levels writes different columns to each.
            for c in st.columns_for(name):
                if c.name in RESERVED:
                    raise ValueError(
                        f"{st.name} declares a column named {c.name!r}, which the "
                        f"store owns")
                cols.setdefault(c.name, sql_type(c.type))
        # Keys are typed by the level, not by whoever declared them first.
        for k in lv.key:
            cols[k] = _key_type(registry, name, k)

        exists_already = con.execute(
            "SELECT count(*) FROM information_schema.tables "
            "WHERE table_name = ?", [name]).fetchone()[0]
        if not exists_already:
            body = ", ".join(f"{quote(c)} {t}" for c, t in cols.items())
            con.execute(f"CREATE TABLE {quote(name)} ({body})")
            # Partitioning by sample is what keeps one sample's writes from
            # rewriting another's files.
            con.execute(f"ALTER TABLE {quote(name)} SET PARTITIONED BY (sample)")
        else:
            have = {r[0] for r in con.execute(
                "SELECT column_name FROM information_schema.columns "
                "WHERE table_name = ?", [name]).fetchall()}
            for c, t in cols.items():
                if c not in have:
                    con.execute(f"ALTER TABLE {quote(name)} ADD COLUMN {quote(c)} {t}")
        out[name] = list(cols)
    return out


def _key_type(registry, level: str, key: str) -> str:
    """A key column's type, taken from whichever stage declares it."""
    for st in registry:
        for c in st.adds:
            if c.name == key:
                return sql_type(c.type)
    return "VARCHAR"


def snapshot_id(con) -> int | None:
    try:
        return con.execute(
            "SELECT max(snapshot_id) FROM lake.snapshots()").fetchone()[0]
    except Exception:
        return None


# --------------------------------------------------------------------------
# provenance — what stage_run replaces
# --------------------------------------------------------------------------
# Every stage used to drop a <output>.manifest.json recording inputs, params,
# stats, timing and tools, and idempotency was a comparison against that file.
# Same fields, one table, so "which stages have run, with what, and how long did
# they take" is a query instead of a directory walk.


def params_hash(params: dict, inputs: list) -> str:
    """Stable identity for a stage invocation.

    Inputs are fingerprinted by size and mtime, as the manifests did - content
    hashing a 30 GB FASTQ to decide whether to skip a stage is not a trade worth
    making.
    """
    import hashlib
    import json

    payload = json.dumps({"params": params, "inputs": inputs},
                         sort_keys=True, default=str)
    return hashlib.sha1(payload.encode()).hexdigest()


def fingerprint(path) -> dict:
    p = Path(path)
    if not p.exists():
        return {"path": str(p), "exists": False}
    st = p.stat()
    return {"path": str(p), "exists": True, "bytes": st.st_size,
            "mtime": round(st.st_mtime, 3)}


def record_run(con, sample: str, stage: str, params: dict, inputs: list,
               where: str | None, stats: dict, tools: dict | None = None,
               seconds: float | None = None, rows: int | None = None) -> None:
    import json

    con.execute("DELETE FROM stage_run WHERE sample = ? AND stage = ?",
                [sample, stage])
    con.execute(
        "INSERT INTO stage_run (sample, stage, params_hash, params, inputs, "
        "where_clause, stats, tools, seconds, rows_written, snapshot_id, written) "
        "VALUES (?,?,?,?,?,?,?,?,?,?,?, now()::TIMESTAMP)",
        [sample, stage, params_hash(params, inputs), json.dumps(params, default=str),
         json.dumps(inputs, default=str), where, json.dumps(stats, default=str),
         json.dumps(tools or {}, default=str),
         None if seconds is None else round(seconds, 2), rows, snapshot_id(con)])


def is_current(con, sample: str, stage: str, params: dict, inputs: list) -> bool:
    """True when this stage already ran for this sample with exactly these inputs."""
    row = con.execute(
        "SELECT params_hash FROM stage_run WHERE sample = ? AND stage = ?",
        [sample, stage]).fetchone()
    return bool(row) and row[0] == params_hash(params, inputs)


def completed_stages(con, sample: str) -> set[str]:
    return {r[0] for r in con.execute(
        "SELECT stage FROM stage_run WHERE sample = ?", [sample]).fetchall()}


def runs(con, sample: str | None = None) -> list[dict]:
    """Every stage run, newest first. This is what the dashboard reads."""
    import json

    sql = ("SELECT sample, stage, params, inputs, where_clause, stats, tools, "
           "seconds, rows_written, snapshot_id, written::VARCHAR FROM stage_run")
    args = []
    if sample is not None:
        sql += " WHERE sample = ?"
        args.append(sample)
    sql += " ORDER BY written DESC"
    out = []
    for r in con.execute(sql, args).fetchall():
        def _j(v):
            try:
                return json.loads(v) if v else {}
            except (TypeError, ValueError):
                return {}
        out.append({
            "sample": r[0], "stage": r[1], "params": _j(r[2]), "inputs": _j(r[3]),
            "where": r[4], "stats": _j(r[5]), "tools": _j(r[6]),
            "seconds": r[7], "rows": r[8], "snapshot": r[9], "written": r[10],
        })
    return out


def register_sample(con, sample: str, source_kind: str = "",
                    source_path: str = "") -> None:
    con.execute("DELETE FROM sample WHERE sample = ?", [sample])
    con.execute("INSERT INTO sample VALUES (?, now()::TIMESTAMP, ?, ?)",
                [sample, source_kind, source_path])
