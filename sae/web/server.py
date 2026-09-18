"""Local web UI for the pipeline: progress, uploads, and launching runs.

The pipeline already writes everything a dashboard needs. Every stage drops a
``<output>.manifest.json`` next to its output recording inputs, params, stats,
timing and tool versions, so progress monitoring is a read over ``work/``
rather than new instrumentation. Runs started outside this UI show up too.

Two sources of truth, deliberately:

* **manifests** - durable state. A stage is done when its manifest exists.
* **the job log** - live state. ``run.py`` prints one line per stage as it
  finishes, so a running job's position comes from tailing its log. Manifests
  cannot supply this: they appear only on completion.

Stdlib only, to keep the project's "no dependency you did not already need"
property. That costs a little hand-rolled routing and buys no new pins.

    python sae/web/server.py                  # http://127.0.0.1:8765
    python sae/web/server.py --read-only      # progress only, no writes
    container/run.sh web                      # inside the image

It lives under ``sae/`` so the container's ``COPY sae /opt/sae/sae`` carries it
into the image with no extra recipe step, and ``SAE_CODE`` shadowing covers it
too. Paths, the interpreter and the bind address all default differently inside
a container, where work and data are bind mounts rather than repo subdirs.

Binding: Apptainer shares the host network namespace, so localhost inside is
localhost outside and the default is right. Docker isolates it, so the default
there is 0.0.0.0 *inside the container* and ``run.sh`` publishes it to the
host's loopback only. This launches subprocesses and writes files; it is a
development tool, not something to expose.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
import threading
import time
import uuid
from dataclasses import dataclass, field
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, urlparse

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent                     # <repo> on a host, /opt/sae inside
PIPELINE = REPO / "sae" / "pipeline"

IN_DOCKER = Path("/.dockerenv").exists()
IN_APPTAINER = bool(os.environ.get("APPTAINER_CONTAINER")
                    or os.environ.get("SINGULARITY_CONTAINER"))
IN_CONTAINER = IN_DOCKER or IN_APPTAINER

# The pipeline describes itself. Stages, their parameters, the entity levels
# and which stage fills which role all come from the registry, so this server
# has no list of stages to keep in step and no schema to mirror - a pipeline it
# has never seen renders the same way. REGISTRY is None when the pipeline
# cannot be imported at all; every reader below degrades to what the work
# directory says rather than failing.
sys.path.insert(0, str(PIPELINE))
try:
    import entities                          # type: ignore
    import lake                              # type: ignore
    from stage import Registry, registry     # type: ignore

    REGISTRY = registry()
except Exception as exc:                     # pragma: no cover - defensive
    print(f"  warning: cannot load the pipeline registry: {exc}", file=sys.stderr)
    entities = lake = REGISTRY = None


def stage_names() -> list[str]:
    return REGISTRY.names if REGISTRY else []


def stage_dirs(run_dir: Path) -> list[str]:
    """Stage directories actually present, registry order first.

    Reading the directory rather than the registry is what lets the UI display
    a run produced by a different pipeline, or by a version of this one with
    stages that have since been removed.
    """
    try:
        found = [p.name for p in sorted(run_dir.iterdir())
                 if p.is_dir() and not p.name.startswith(".")]
    except OSError:
        return []
    known = [n for n in stage_names() if n in found]
    return known + [n for n in found if n not in known]

# Sample names and uploaded filenames become path components, so they are
# restricted rather than escaped - no separators, no leading dot, no traversal.
SAFE_NAME = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,63}$")
MAX_UPLOAD = 16 * 1024**3                    # 16 GiB; FASTQ gets large
LOG_TAIL = 64 * 1024
MAX_LISTING = 500
SEQ_SUFFIXES = (".gz", ".fa", ".fasta", ".fna", ".faa", ".fastq")

# Artifact preview. The frontend knows only three shapes - text, table, json -
# and never a stage's schema, so a stage may change its columns, or a new stage
# may appear, without touching any JavaScript. A stage that wants a curated
# summary just drops a .md or .tsv in its work directory and it shows up.
TEXT_SUFFIXES = (".faa", ".fa", ".fasta", ".fna", ".fastq", ".txt", ".md",
                 ".log", ".err", ".out", ".sam", ".gfa", ".bed")
TABLE_SUFFIXES = (".tsv", ".csv", ".parquet")
PREVIEW_LINES = 200
PREVIEW_ROWS = 100
PREVIEW_BYTES = 2 * 1024**2

# Feature map. s06 writes long-format (gene_id, feature_id, activation) with
# top-K per gene, so a run is a sparse matrix over the 16,384-wide codebook.
# Reduced to 2D it gives the same kind of picture as the ESM Atlas map, at a
# scale that needs no tiling: a sample is thousands of points, not millions.
# UMAP only. It is the one method here with a `transform`, which is what makes
# a fixed reference layout possible; t-SNE cannot place a new point in an
# existing layout at all, and a linear projection was never good enough.
MAX_CONTEXT_POINTS = 8000
MAX_PROJECTION_POINTS = 20000
# The cohort view draws every sample at once. The ceiling is the browser's,
# not the store's: these are plain SVG circles with no handlers, and past
# roughly this many the first paint starts to drag.
MAX_ATLAS_POINTS = 30000
# Selecting on the atlas. A drag is resolved in the browser and arrives here
# as indices into the layout the browser was handed - not as gene ids and not
# as a predicate. The ids behind those indices are derived from the store
# deterministically at a snapshot, so index `i` means the same protein on both
# sides for as long as that snapshot holds; the request carries the snapshot so
# a stale one is refused rather than answered with the wrong proteins. It also
# keeps a 30,000-point selection to a few kilobytes of integers.
MAX_SELECTION_SHOWN = 200        # rows sent back for the panel to display
MAX_SELECTION_FASTA = 20000      # sequences one download may carry
SELECT_CHUNK = 500               # gene ids per IN list, so the SQL stays sane
# How many projections to keep placed at once, so switching between them is
# free after the first visit to each.
# Enough for every layout on disk plus "fit on this cohort", so a demo can
# move between them without paying for a re-transform. Coordinates only.
ATLAS_CACHE_KEEP = 6

# UMAP is numba, and numba's default `workqueue` threading layer is not
# threadsafe: called from two Python threads at once it does not raise, it
# aborts the process. This is a ThreadingHTTPServer, so two overlapping
# projections - trivially reachable now that a run can be drawn against a
# second one - took the whole dashboard down. Serialising them costs nothing
# real, because a projection is CPU-bound and gains nothing from running
# beside another; the worst case is that the second request waits.
_PROJECT_LOCK = threading.Lock()

# What a launch request may contain is derived from the registry rather than
# listed here: a parameter exists if some stage declares it, and it is valid if
# that stage's own Param accepts it. The hand-kept mirror of run.py's argparse
# that used to live here could only ever drift.


def defaults() -> dict:
    """Container and host disagree about where everything lives."""
    if IN_CONTAINER:
        return {
            "work": [Path("/work")],
            "uploads": Path("/work/uploads"),
            "data": [Path("/data")],
            "python": "/opt/venv/bin/python",
            "reference": Path("/data/reference_map.joblib"),
            # Apptainer shares the host network namespace, so localhost is
            # already right. Nested in Docker it is not, and run.sh passes an
            # explicit --host 0.0.0.0 for that case.
            "host": "0.0.0.0" if IN_DOCKER else "127.0.0.1",
        }
    return {
        "work": [REPO / "work", PIPELINE / "work"],
        "uploads": REPO / "uploads",
        "data": [REPO / "data"],
        "python": str(REPO / "sae" / ".venv" / "bin" / "python"),
        "reference": REPO / "data" / "reference_map.joblib",
        "host": "127.0.0.1",
    }


# --------------------------------------------------------------------------
# reading pipeline state
# --------------------------------------------------------------------------
def _inner_suffix(p: Path) -> str:
    """Suffix ignoring a trailing .gz, so foo.tsv.gz reads as a table."""
    return (p.suffixes[-2] if p.suffix == ".gz" and len(p.suffixes) > 1
            else p.suffix).lower()


def artifact_kind(p: Path) -> str:
    suf = _inner_suffix(p)
    if suf == ".json":
        return "json"
    if suf in TABLE_SUFFIXES:
        return "table"
    if suf in TEXT_SUFFIXES:
        return "text"
    return "binary"


def _open_text(p: Path):
    import gzip
    return gzip.open(p, "rt", errors="replace") if p.suffix == ".gz" \
        else open(p, "r", errors="replace")


def preview(p: Path, limit: int) -> dict:
    """Normalise any artifact into text, table or json. Never stage-specific."""
    kind = artifact_kind(p)
    if kind == "json":
        raw = p.read_bytes()[:PREVIEW_BYTES]
        try:
            return {"kind": "json", "json": json.loads(raw)}
        except json.JSONDecodeError:
            return {"kind": "text", "text": raw.decode(errors="replace")}

    if kind == "table":
        suf = _inner_suffix(p)
        if suf == ".parquet":
            # pyarrow is a pipeline dependency, not a web one. Import it only
            # when a parquet is actually asked for, so the server still runs
            # anywhere without it.
            try:
                import pyarrow.parquet as pq
            except ImportError:
                return {"kind": "text",
                        "text": "pyarrow is not installed, so this parquet "
                                "cannot be previewed here."}
            f = pq.ParquetFile(p)
            cols = [c.name for c in f.schema_arrow]
            rows: list[list] = []
            for batch in f.iter_batches(batch_size=min(limit, 1000)):
                for row in batch.to_pylist():
                    rows.append([row.get(c) for c in cols])
                    if len(rows) >= limit:
                        break
                if len(rows) >= limit:
                    break
            return {"kind": "table", "columns": cols, "rows": rows,
                    "total_rows": f.metadata.num_rows}
        import csv
        delim = "\t" if suf == ".tsv" else ","
        with _open_text(p) as fh:
            r = csv.reader(fh, delimiter=delim)
            try:
                cols = next(r)
            except StopIteration:
                return {"kind": "table", "columns": [], "rows": []}
            rows = [row for _, row in zip(range(limit), r)]
        return {"kind": "table", "columns": cols, "rows": rows}

    if kind == "text":
        with _open_text(p) as fh:
            lines = [ln.rstrip("\n") for _, ln in zip(range(limit), fh)]
        return {"kind": "text", "text": "\n".join(lines)}

    return {"kind": "binary", "text": f"{p.name} is not a previewable format."}


# One definition of "read s06 output" and "normalise it", shared with the
# reference-map builder so a run is projected exactly as the corpus was fitted.
from launcher import (Job, failure_hint,  # noqa: E402
                      load_jobs, make_launcher, save_job, script_text)
import reference_map                                                        # noqa: E402
from reference_map import load as load_reference                            # noqa: E402
from reference_map import normalise, sparse_from_table                      # noqa: E402
import metadata                                           # noqa: E402


# UMAP runs on numba, whose default workqueue threading layer is **not**
# threadsafe: two projections at once do not merely contend, they terminate the
# process ("Concurrent access has been detected"). This is a ThreadingHTTPServer,
# so two browser tabs - or one tab and one curl - are enough. Every layout goes
# through this lock.
_LAYOUT_LOCK = threading.Lock()


def _fit_layout(m):
    """The fallback fit, for when there is no reference map to borrow.

    Same parameters as ``reference_map.build``, so a picture fitted here and
    one transformed into the corpus differ by what they were fitted on and by
    nothing else.
    """
    import umap

    with _LAYOUT_LOCK:
        return umap.UMAP(n_components=2, metric="cosine", random_state=0,
                         n_neighbors=max(2, min(15, m.shape[0] - 1)),
                         min_dist=0.1).fit_transform(m)


# Densifying the whole matrix at once is 2.1 GB at cohort scale, which is fine
# on a workstation and not fine everywhere. Transform in slices instead: the
# layout is per-row, so the result is identical and the peak is bounded.
DENSE_CHUNK = 4000


def _transform(ref, m):
    """Place rows in the reference layout. Same lock, same reason.

    A map fitted on a dense matrix will not accept a sparse one - and above
    ~4096 points fitting dense is the only way to get a layout that survives
    being saved and loaded again, so this is the normal case for any map built
    over a real cohort.
    """
    import numpy as np

    dense = ref.get("input_form") == "dense"
    with _LAYOUT_LOCK:
        if not dense:
            return ref["reducer"].transform(m)
        out = []
        for i in range(0, m.shape[0], DENSE_CHUNK):
            block = m[i:i + DENSE_CHUNK]
            out.append(ref["reducer"].transform(
                np.asarray(block.todense(), dtype=np.float32)))
        return np.vstack(out) if len(out) > 1 else out[0]


def _cap_per_sample(ids, m, limit):
    """Thin to `limit` rows while keeping every sample represented.

    A flat stride over the cohort would thin each sample in proportion to its
    size, which is fine until the smallest sample rounds to nothing. Here the
    budget is shared out proportionally with a floor, so a nine-gene sample
    still appears.
    """
    import numpy as np

    if len(ids) <= limit:
        return ids, m, False
    groups: dict[str, list[int]] = {}
    for i, key in enumerate(ids):
        groups.setdefault(key.split("\x1f", 1)[0], []).append(i)
    floor = min(200, limit // max(1, len(groups)))
    keep: list[int] = []
    for members in groups.values():
        share = min(len(members), max(floor, round(limit * len(members) / len(ids))))
        if share >= len(members):
            keep.extend(members)
            continue
        # Evenly spaced picks rather than a stride: a stride can only halve,
        # so asking for 93% of a sample would hand back 50%.
        pick = np.unique(np.linspace(0, len(members) - 1, share).round().astype(int))
        keep.extend(members[i] for i in pick)
    keep.sort()
    return [ids[i] for i in keep], m[np.asarray(keep)], True


def _cap(ids, m, limit):
    """Thin a sample to `limit` rows by a deterministic stride."""
    if len(ids) <= limit:
        return ids, m, False
    step = -(-len(ids) // limit)                  # ceil, so the result fits
    return ids[::step], m[::step], True


def _points(ids, xy, rows: dict, column: str | None, domain: dict | None,
            keep: set | None, side: str | None) -> list[dict]:
    """Drawable points, carrying the value the plot is coloured by.

    `keep` is the filter: None means no filter was applied, which is not the
    same as a filter that matched nothing. The slot is resolved here rather
    than in the browser so that binning a numeric column - and deciding what
    counts as the folded tail - happens in exactly one place.
    """
    out = []
    for i, g in enumerate(ids):
        if keep is not None and g not in keep:
            continue
        p = {"id": g, "x": round(float(xy[i][0]), 3),
             "y": round(float(xy[i][1]), 3)}
        if column:
            v = (rows.get(g) or {}).get(column)
            p["v"] = v
            p["s"] = metadata.slot_of(v, domain) if domain else None
        if side:
            p["side"] = side
        out.append(p)
    return out


def _fasta_record(row: dict, names: list[str], width: int = 60) -> str:
    """One FASTA record for a selected protein, wrapped like the pipeline's.

    The header is ``sample|gene_id`` followed by every non-empty column the
    gene level holds, as ``key=value``. Listing them all rather than a chosen
    few is the same rule the rest of this server follows: nothing here knows
    what a stage's columns mean, so a download from a store with a taxon call
    in it carries the taxon call without anyone adding it.

    Whitespace is stripped from values because it would split the header into
    fields that were never there.
    """
    fields = []
    for n in names:
        if n == "gene_id" or row.get(n) is None:
            continue
        v = "_".join(str(row[n]).split())
        if v:
            fields.append(f"{n}={v}")
    seq = row["seq"]
    head = f">{row['sample']}|{row['gene_id']}"
    lines = [head + (" " + " ".join(fields) if fields else "")]
    lines += [seq[i:i + width] for i in range(0, len(seq), width)]
    return "\n".join(lines) + "\n"


def _summarise(mf: Path) -> dict | None:
    try:
        d = json.loads(mf.read_text())
    except (json.JSONDecodeError, OSError):
        return None
    out = d.get("output") or {}
    return {
        "output": Path(out.get("path", "")).name,
        "bytes": out.get("bytes"),
        "stats": d.get("stats") or {},
        "params": d.get("params") or {},
        "tools": d.get("tools") or {},
        "seconds": d.get("seconds"),
        "written": d.get("written"),
        "inputs": [Path(i.get("path", "")).name for i in (d.get("inputs") or [])],
        # Columns this stage contributed, and to which level. The frontend
        # renders these without knowing what any of them mean.
        "tables": [{"level": t.get("level"), "role": t.get("role"),
                    "where": t.get("where"),
                    "columns": [c.get("name") for c in (t.get("columns") or [])
                                if c.get("name") not in (t.get("key") or [])]}
                   for t in (d.get("tables") or [])],
    }


def _contributed(stage_name: str, where: str | None) -> list[dict]:
    """The columns a stage writes, per level, for the run detail view."""
    if REGISTRY is None or stage_name not in REGISTRY:
        return [{"level": "", "role": "", "where": where, "columns": []}] \
            if where else []
    st = REGISTRY[stage_name]
    out = []
    for level in st.outputs:
        cols = [c.name for c in st.columns_for(level)]
        if not cols:
            continue
        key = entities.LEVELS[level].key if level in entities.LEVELS else ()
        out.append({
            "level": level,
            "role": "annotation" if level == st.consumes else "base",
            "where": where,
            "columns": [c for c in cols if c not in key],
        })
    return out


def scan_lake(target) -> list[dict]:
    """Every sample in the store, with what has run against it.

    This replaced a walk that re-read every manifest under every sample on a
    two-second poll - roughly 5,300 JSON parses per tick at cohort scale. It is
    now two queries, and the read attach is brief so it does not hold the store
    against a running job.
    """
    if lake is None:
        return []
    try:
        with lake.read(target, budget=3) as con:
            runs = lake.runs(con)
            levels: dict[str, list[str]] = {}
            for lvl in entities.LEVELS:
                try:
                    for (smp,) in con.execute(
                            f"SELECT DISTINCT sample FROM {lvl}").fetchall():
                        levels.setdefault(smp, []).append(lvl)
                except Exception:
                    continue
    except Exception:
        return []

    by_sample: dict[str, dict] = {}
    for r in runs:
        d = by_sample.setdefault(r["sample"], {"stages": {}, "order": [],
                                               "updated": None})
        d["stages"][r["stage"]] = {
            "state": "done", "stats": r["stats"], "params": r["params"],
            "seconds": r["seconds"], "written": r["written"],
            "output": f"{r['rows'] or 0} rows" if r["rows"] else "",
            "bytes": None, "inputs": [Path(i.get("path", "")).name
                                      for i in (r["inputs"] or []) if isinstance(i, dict)],
            # Which columns this stage contributed, and to which level. Taken
            # from the registry rather than from the run record: the stage
            # declares them, the store's DDL is generated from that same
            # declaration, so it is the truth rather than a copy of it.
            "tables": _contributed(r["stage"], r["where"]),
        }
        if r["written"] and (d["updated"] is None or r["written"] > d["updated"]):
            d["updated"] = r["written"]

    # A sample can hold rows without a stage_run record - anything backfilled
    # from a work directory whose manifest was missing. It still has data, so
    # it still belongs in the list.
    for sample in levels:
        by_sample.setdefault(sample, {"stages": {}, "order": [], "updated": None})

    known = REGISTRY.names if REGISTRY else []
    out = []
    for sample, d in by_sample.items():
        have = list(d["stages"])
        order = [n for n in known if n in have] + [n for n in have if n not in known]
        out.append({
            "id": f"0:{sample}", "sample": sample, "root": str(target),
            "stages": d["stages"], "order": order,
            "levels": levels.get(sample, []),
            "updated": d["updated"],
            "done": len(d["stages"]),
        })
    out.sort(key=lambda r: (r["updated"] or ""), reverse=True)
    return out


def scan(roots: list[Path]) -> list[dict]:
    """Every sample across every work root, with per-stage manifest summaries.

    Sample names collide across roots (a CHI-A in each is normal), so the id
    carries the root index rather than the bare name.
    """
    runs = []
    for idx, root in enumerate(roots):
        if not root.is_dir():
            continue
        for sample_dir in sorted(p for p in root.iterdir() if p.is_dir()):
            if sample_dir.name.startswith("."):
                continue
            stages, latest = {}, None
            for stage in stage_dirs(sample_dir):
                sd = sample_dir / stage
                mfs = sorted(sd.glob("*.manifest.json"))
                if not mfs:
                    stages[stage] = {"state": "started"}   # dir but no manifest
                    continue
                s = _summarise(mfs[0])
                if s is None:
                    stages[stage] = {"state": "unreadable"}
                    continue
                stages[stage] = {"state": "done", **s}
                if s["written"] and (latest is None or s["written"] > latest):
                    latest = s["written"]
            if stages:
                runs.append({
                    "id": f"{idx}:{sample_dir.name}",
                    "sample": sample_dir.name,
                    "root": str(root),
                    "stages": stages,
                    # The order these ran in, as this run actually has them -
                    # not the registry's, so a run from another pipeline still
                    # draws a sensible strip.
                    "order": list(stages),
                    "levels": entities.levels_present(sample_dir) if entities else [],
                    "updated": latest,
                    "done": sum(1 for v in stages.values() if v["state"] == "done"),
                })
    runs.sort(key=lambda r: (r["updated"] or ""), reverse=True)
    return runs


# --------------------------------------------------------------------------
# launching runs
# --------------------------------------------------------------------------
def stage_from_log(log: Path, live: bool) -> str | None:
    """run.py prints '[stage] ...' as each stage completes."""
    try:
        text = log.read_text(errors="replace")
    except OSError:
        return None
    seen = re.findall(r"^\[(\w+)\]", text, re.MULTILINE)
    if not seen:
        return None
    last = seen[-1]
    if not live:
        return last
    # The line is printed on completion, so the next stage is in flight. Which
    # one that is comes from the plan run.py printed, not from a fixed order -
    # a run with --only or --skip has a different sequence, and a run from
    # another pipeline has different stages entirely.
    planned = re.findall(r"^\s{2}(\w+)\s+\S+ -> ", text, re.MULTILINE) \
        or stage_names()
    if last in planned:
        i = planned.index(last)
        return planned[i + 1] if i + 1 < len(planned) else last
    return last


class Jobs:
    """Registry over a launcher. Submission and polling are the launcher's."""

    def __init__(self, launcher, log_dir: Path):
        self.launcher = launcher
        self.log_dir = log_dir
        self._lock = threading.Lock()
        # Adopt anything from a previous session. A job outlives this process,
        # so its state is recovered from the pid and the exit-code file rather
        # than from anything we remembered.
        adopted = load_jobs(log_dir)
        self._jobs: dict[str, Job] = {j.id: j for j in adopted}
        # The launcher has to know too: one still running must be waited for
        # before anything else starts, and one that never started has to be
        # queued again or it sits pending forever.
        launcher.adopt(adopted)

    def submit(self, argv: list[str], sample: str, cwd: Path) -> Job:
        job = self.launcher.submit(argv, sample, cwd)
        save_job(job, self.log_dir)
        with self._lock:
            self._jobs[job.id] = job
        return job

    def _as_dict(self, job: Job) -> dict:
        self.launcher.refresh(job)
        live = job._state in ("running", "pending")
        return {
            "id": job.id, "backend": job.backend,
            "backend_id": job.backend_id, "sample": job.sample,
            "state": job._state, "returncode": job._rc,
            "started": job.started, "elapsed": job.elapsed(),
            "argv": job.argv, "stage": stage_from_log(job.log, live),
            "hint": failure_hint(job),
            "log_path": str(job.log),
            "script_path": str(job.script) if job.script else None,
        }

    def all(self) -> list[dict]:
        with self._lock:
            jobs = list(self._jobs.values())
        return sorted((self._as_dict(j) for j in jobs),
                      key=lambda d: d["started"], reverse=True)

    def get(self, jid: str) -> Job | None:
        with self._lock:
            return self._jobs.get(jid)

    def cancel(self, jid: str) -> bool:
        job = self.get(jid)
        return bool(job) and self.launcher.cancel(job)


def build_argv(body: dict, allowed: list[Path]) -> list[str]:
    """Translate a request body into run.py flags.

    Validation is the registry's: a parameter is accepted because some stage
    declares it and that stage's own ``Param`` coerces the value, so this
    function has nothing to keep in step with the pipeline. Predicates go
    through ``guard_predicate`` for the same reason the pipeline does - the
    server writes files, and a WHERE clause can reach a COPY.

        {"sample": "CHI-A", "input_kind": "contigs", "input_path": "...",
         "target": "feature",
         "params": {"s05_prefilter": {"hmm": "/data/pfam/Pfam-A.hmm"}},
         "where":  {"s06_embed": "category = 'dark' AND aa_len > 200"}}
    """
    if REGISTRY is None:
        raise ValueError("the pipeline registry is unavailable on this server")

    argv: list[str] = []
    kind = body.get("input_kind")
    if kind:
        if entities and kind not in entities.FILE_KINDS:
            raise ValueError(f"input_kind must be one of "
                             f"{', '.join(entities.FILE_KINDS)}")
        flag = {"reads": "fastq"}.get(kind, kind)
        argv += [f"--{flag}", checked_path(body.get("input_path"), allowed)]
        if body.get("input_path2"):
            argv += ["--fastq2", checked_path(body["input_path2"], allowed)]
    elif body.get("start"):
        argv += ["--from", str(body["start"])]
    else:
        raise ValueError("give input_kind with input_path, or a start port")

    if body.get("target"):
        argv += ["--to", str(body["target"])]
    if body.get("force"):
        argv.append("--force")
    for name in body.get("only") or []:
        if name not in REGISTRY:
            raise ValueError(f"no stage named {name!r}")
        argv += ["--only", name]
    for name in body.get("skip") or []:
        if name not in REGISTRY:
            raise ValueError(f"no stage named {name!r}")
        argv += ["--skip", name]

    for stage_name, values in (body.get("params") or {}).items():
        if stage_name not in REGISTRY:
            raise ValueError(f"no stage named {stage_name!r}")
        st = REGISTRY[stage_name]
        for key, value in (values or {}).items():
            param = st.param(key)
            if param is None:
                raise ValueError(
                    f"{stage_name} has no parameter {key!r}; have "
                    + ", ".join(p.name for p in st.params))
            if value in (None, ""):
                continue
            if param.path:
                # A path parameter is still a path, whoever declared it.
                argv += ["--set", f"{stage_name}.{key}={checked_path(value, allowed)}"]
                continue
            coerced = param.coerce(value)          # raises ValueError on a bad one
            if param.type is bool and not coerced:
                continue
            argv += ["--set", f"{stage_name}.{key}={coerced}"]

    for stage_name, expr in (body.get("where") or {}).items():
        if stage_name not in REGISTRY:
            raise ValueError(f"no stage named {stage_name!r}")
        if not REGISTRY[stage_name].selectable:
            raise ValueError(f"{stage_name} does not take a selection")
        if not expr or not str(expr).strip():
            continue
        argv += ["--where", f"{stage_name}={entities.guard_predicate(str(expr))}"]
    return argv


def checked_path(raw, allowed: list[Path]) -> str:
    """Resolve an input path and require it to sit under an allowed root."""
    if not raw:
        raise ValueError("input_path is required")
    p = Path(str(raw)).expanduser()
    p = (p if p.is_absolute() else (REPO / p)).resolve()
    if not p.is_file():
        raise ValueError(f"no such file: {p}")
    if not any(p == root or root in p.parents for root in allowed):
        raise ValueError("input is outside the permitted directories")
    return str(p)


# --------------------------------------------------------------------------
# HTTP
# --------------------------------------------------------------------------
class Handler(BaseHTTPRequestHandler):
    server_version = "sae-web"
    cfg: dict = {}

    def log_message(self, fmt, *args):
        if self.cfg.get("verbose"):
            super().log_message(fmt, *args)

    # -- helpers
    def _send(self, code: int, body: bytes, ctype: str,
              filename: str | None = None):
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        if filename:
            # The one response meant to become a file on disk rather than a
            # value in the page. The name is composed here, never echoed from
            # the request, so there is nothing in it to escape.
            self.send_header("Content-Disposition",
                             f'attachment; filename="{filename}"')
        self.end_headers()
        if self.command != "HEAD":
            self.wfile.write(body)

    def _json(self, payload, code: int = 200):
        self._send(code, json.dumps(payload).encode(), "application/json")

    def _err(self, code: int, msg: str):
        self._json({"error": msg}, code)

    def _guard_writes(self) -> bool:
        if self.cfg["read_only"]:
            self._err(403, "server is running in --read-only mode")
            return False
        return True

    # -- routes
    def do_GET(self):
        u = urlparse(self.path)
        q = parse_qs(u.query)
        if u.path in ("/", "/index.html"):
            page = HERE / "index.html"
            if not page.is_file():
                return self._err(500, "index.html missing")
            return self._send(200, page.read_bytes(), "text/html; charset=utf-8")
        if u.path == "/api/state":
            return self._json({
                "runs": scan_lake(self.cfg["lake"]),
                "jobs": self.cfg["jobs"].all(),
                "stages": stage_names(),
                "roots": [str(r) for r in self.cfg["roots"]],
                "lake": str(self.cfg["lake"]),
                # Moves exactly when something is written, so a client can tell
                # a cached view is stale without polling the view itself.
                "snapshot": self._snapshot(),
                "read_only": self.cfg["read_only"],
                "uploads": str(self.cfg["uploads"]),
                "data_files": self._candidate_data(),
                "container": "docker" if IN_DOCKER else
                             ("apptainer" if IN_APPTAINER else None),
                "launcher": self.cfg["jobs"].launcher.describe(),
            })
        if u.path == "/api/log":
            job = self.cfg["jobs"].get((q.get("id") or [""])[0])
            if job is None:
                return self._err(404, "no such job")
            self.cfg["jobs"].launcher.refresh(job)
            try:
                data = job.log.read_bytes()[-LOG_TAIL:]
            except OSError:
                data = b""
            return self._json({"id": job.id, "state": job._state,
                               "log": data.decode(errors="replace")})
        if u.path == "/api/script":
            # Static once submitted, so it is its own route rather than a
            # rider on /api/log: the log is polled every couple of seconds and
            # the script would be re-sent with every poll for nothing.
            job = self.cfg["jobs"].get((q.get("id") or [""])[0])
            if job is None:
                return self._err(404, "no such job")
            return self._json({
                "id": job.id,
                "path": str(job.script) if job.script else None,
                "script": script_text(job),
            })
        if u.path == "/api/pipeline":
            return self._pipeline()
        if u.path == "/api/plan":
            return self._plan((q.get("have") or [""])[0],
                              (q.get("want") or [""])[0],
                              q.get("skip") or [])
        if u.path == "/api/columns":
            return self._columns((q.get("run") or [""])[0],
                                 (q.get("level") or [""])[0])
        if u.path == "/api/query":
            try:
                limit = int((q.get("limit") or ["20"])[0])
            except ValueError:
                limit = 20
            return self._query((q.get("run") or [""])[0],
                               (q.get("level") or ["gene"])[0],
                               (q.get("where") or [""])[0], limit,
                               (q.get("scope") or ["sample"])[0])
        if u.path == "/api/inputs":
            return self._json({"files": self._candidate_inputs()})
        if u.path == "/api/artifacts":
            return self._artifacts((q.get("run") or [""])[0])
        if u.path == "/api/projection":
            # `b` overlays a second run in the same layout. Optional, so the
            # single-run URL is unchanged.
            return self._projection((q.get("run") or [""])[0],
                                    (q.get("mode") or ["reference"])[0],
                                    (q.get("b") or [""])[0],
                                    (q.get("color") or [""])[0],
                                    (q.get("filter") or [""])[0],
                                    (q.get("map") or [""])[0])
        if u.path == "/api/maps":
            return self._json({
                "maps": self._maps(),
                "default": Path(self.cfg["reference_path"]).stem,
            })
        if u.path == "/api/atlas":
            try:
                limit = int((q.get("limit") or ["0"])[0]) or MAX_ATLAS_POINTS
            except ValueError:
                limit = MAX_ATLAS_POINTS
            return self._atlas((q.get("color") or [""])[0],
                               max(500, min(limit, MAX_ATLAS_POINTS)),
                               (q.get("map") or [""])[0])
        if u.path == "/api/preview":
            try:
                limit = int((q.get("limit") or ["0"])[0]) or PREVIEW_ROWS
            except ValueError:
                limit = PREVIEW_ROWS
            return self._preview((q.get("run") or [""])[0],
                                 (q.get("stage") or [""])[0],
                                 (q.get("file") or [""])[0], limit)
        return self._err(404, "not found")

    def do_POST(self):
        u = urlparse(self.path)
        q = parse_qs(u.query)
        if u.path == "/api/upload":
            return self._upload((q.get("name") or [""])[0])
        if u.path == "/api/run":
            return self._launch()
        if u.path == "/api/selection":
            # A read, so no write guard: it resolves picked points to the rows
            # behind them and hands back their sequences.
            return self._selection((q.get("format") or ["json"])[0])
        if u.path == "/api/cancel":
            if not self._guard_writes():
                return
            jid = (q.get("id") or [""])[0]
            if self.cfg["jobs"].get(jid) is None:
                return self._err(404, "no such job")
            return self._json({"id": jid,
                               "cancelled": self.cfg["jobs"].cancel(jid)})
        return self._err(404, "not found")

    # -- handlers
    def _candidate_data(self) -> list[dict]:
        """Files under data/ that some path parameter could be pointed at.

        The suffixes come from the parameters themselves, so a new stage that
        wants a ``.dmnd`` gets a working file picker without this server
        learning what Diamond is. Without a Pfam database, for instance, s05 is
        a pass-through and every protein reaches the GPU - the UI should not
        make finding it a matter of knowing a path.
        """
        wanted = {sfx.lower()
                  for st in (REGISTRY or []) for prm in st.params
                  if prm.path for sfx in prm.suffixes}
        if not wanted:
            return []
        out = []
        for root in self.cfg["data"]:
            if not root.is_dir():
                continue
            for p in sorted(root.rglob("*")):
                if len(out) >= MAX_LISTING:
                    break
                if not p.is_file() or p.name.startswith("."):
                    continue
                if p.suffix.lower() not in wanted:
                    continue
                out.append({"path": str(p), "name": p.name,
                            "suffix": p.suffix.lower(),
                            "bytes": p.stat().st_size})
        return out

    def _candidate_inputs(self) -> list[dict]:
        """Files that can start a run: uploads plus anything under data/."""
        seen, out = set(), []
        for root in [self.cfg["uploads"], *self.cfg["data"]]:
            if not root.is_dir():
                continue
            for p in sorted(root.rglob("*")):
                if len(out) >= MAX_LISTING:
                    break
                if not p.is_file() or p.name.startswith("."):
                    continue
                if p.suffix not in SEQ_SUFFIXES:
                    continue
                sp = str(p)
                if sp in seen:
                    continue
                seen.add(sp)
                out.append({"path": sp, "name": p.name, "bytes": p.stat().st_size})
        return out

    def _sample(self, run_id: str) -> str:
        """`0:<sample>` -> sample. The index is vestigial: there is one store."""
        _, _, sample = str(run_id).partition(":")
        sample = sample or str(run_id)
        if not SAFE_NAME.match(sample):
            raise ValueError("bad run id")
        return sample

    def _run_dir(self, run_id: str) -> Path:
        """Resolve '<root index>:<sample>' to a sample directory, or raise."""
        idx, _, sample = run_id.partition(":")
        roots = self.cfg["roots"]
        if not idx.isdigit() or not SAFE_NAME.match(sample):
            raise ValueError("bad run id")
        i = int(idx)
        if i >= len(roots):
            raise ValueError("bad run id")
        d = (roots[i] / sample).resolve()
        if roots[i] not in d.parents or not d.is_dir():
            raise ValueError("no such run")
        return d

    def _pipeline(self):
        """The whole pipeline, described well enough to render a UI from.

        Stages, their parameters with types / defaults / choices / help, the
        entity levels and their keys, which external tools are actually
        installed. The frontend builds its form out of this and holds no
        knowledge of any stage.
        """
        if REGISTRY is None:
            return self._err(503, "the pipeline registry is unavailable")
        body = REGISTRY.to_json()
        body["roles"] = sorted({r for s in REGISTRY for r in s.roles})
        return self._json(body)

    def _plan(self, have: str, want: str, skip: list[str]):
        """The stages that would run, resolved by the driver's own planner.

        The page used to work this out itself, walking ports greedily. That
        held only while one stage consumed each port: as soon as two did -
        assemble and translate both take reads - the greedy walk returned the
        union of both routes and could not tell that skipping the assembler
        leaves no route to contigs at all. So the preview asks the planner
        instead of imitating it, and what it shows is what will run.
        """
        if REGISTRY is None:
            return self._err(503, "the pipeline registry is unavailable")
        names = [n for n in skip if n in REGISTRY]
        reg = (Registry([s for s in REGISTRY if s.name not in names])
               if names else REGISTRY)
        body = {"have": have, "want": want, "skip": names}
        try:
            body["stages"] = [s.name for s in reg.plan(have, want)]
        except Exception as exc:
            # No route is an answer, not a failure: it is what tells the page
            # an option is not on offer.
            body["stages"] = []
            body["unreachable"] = str(exc)
        return self._json(body)

    def _columns(self, run_id: str, level: str):
        """What can be predicated on, and which stage wrote each column.

        Read from the store's own schema plus the registry, so a column is
        offered because it exists and is explained because a stage declared it.
        """
        if entities is None:
            return self._err(503, "the pipeline registry is unavailable")
        try:
            sample = self._sample(run_id)
        except ValueError as exc:
            return self._err(404, str(exc))
        try:
            with lake.read(self.cfg["lake"], budget=5) as con:
                levels = entities.levels_present(con, sample)
                lvl = level or (levels[0] if levels else "")
                if not lvl:
                    return self._json({"run": run_id, "sample": sample,
                                       "levels": [], "level": None, "columns": []})
                if lvl not in entities.LEVELS:
                    return self._err(404, f"no level named {lvl!r}")
                cols = entities.columns_of(con, lvl)
                n = entities.count(con, lvl, sample=sample)
                cohort = entities.count(con, lvl)
                everyone = entities.samples(con)
        except Exception as exc:
            return self._err(503, f"cannot read the store: {exc}")
        return self._json({
            "run": run_id, "sample": sample, "levels": levels, "level": lvl,
            "rows": n, "cohort_rows": cohort,
            "key": list(entities.LEVELS[lvl].key),
            "columns": cols, "collisions": [],
            # Every sample is queryable by name in a predicate, because
            # `sample` is a column rather than a separate schema.
            "samples": everyone,
        })

    def _query(self, run_id: str, level: str, where: str, limit: int,
               scope: str = "sample"):
        """How many rows a predicate selects, and a look at them.

        The point is to see what a selection takes before spending a GPU on it.
        `scope=cohort` drops the sample filter, so the same box answers a
        cohort question - which is the thing the old per-sample schemas could
        not express at all.
        """
        if entities is None:
            return self._err(503, "the pipeline registry is unavailable")
        try:
            sample = self._sample(run_id)
        except ValueError as exc:
            return self._err(404, str(exc))
        try:
            expr = entities.guard_predicate(where) or None
        except ValueError as exc:
            return self._err(400, str(exc))
        only = None if scope == "cohort" else sample
        try:
            with lake.read(self.cfg["lake"], budget=5) as con:
                n = entities.count(con, level, expr, sample=only)
                # `seq` is megabytes of amino acids nobody reads in a preview.
                names = [c["name"] for c in entities.columns_of(con, level)
                         if c["name"] != "seq"]
                rows = entities.select(con, level, expr, columns=names or None,
                                       limit=max(1, min(limit, 500)), sample=only)
        except Exception as exc:
            return self._err(422, f"{type(exc).__name__}: {exc}".strip()[:400])
        return self._json({"run": run_id, "sample": sample, "level": level,
                           "where": expr, "scope": scope, "matched": n,
                           "columns": rows.column_names, "rows": rows.to_pylist()})

    def _artifacts(self, run_id: str):
        try:
            run_dir = self._run_dir(run_id)
        except ValueError as exc:
            return self._err(404, str(exc))
        out: dict[str, list] = {}
        for stage in stage_dirs(run_dir):
            sd = run_dir / stage
            files = []
            for f in sorted(sd.iterdir()):
                if not f.is_file() or f.name.startswith("."):
                    continue
                files.append({"name": f.name, "bytes": f.stat().st_size,
                              "kind": artifact_kind(f)})
            if files:
                out[stage] = files
        return self._json({"run": run_id, "stages": out})

    def _maps(self) -> list[dict]:
        """Every projection that can be offered, the store's first.

        A projection kept in the store needs nothing else to be usable - which
        is the point: hand someone a dump of the lake and they get the same
        picker with the same layouts, no 100 MB pickles to pass alongside.
        Listing them is a query, so this stays cheap enough to call on load.

        A `.joblib` on disk that the store does not have is still offered,
        because it can place points by transforming even though it cannot do it
        by lookup. One that the store *does* have is not listed twice.
        """
        out, seen = [], set()
        if lake is not None:
            try:
                with lake.read(self.cfg["lake"], budget=5) as con:
                    for pr in lake.projections(con):
                        seen.add(pr["name"])
                        out.append({
                            "name": pr["name"], "n": pr["n"] or pr["stored"],
                            "built": pr["built"], "note": pr["note"],
                            "params": pr["params"],
                            # What it was fitted on, by name - that is the
                            # difference a person is choosing between. Role
                            # counts are carried separately rather than
                            # standing in for it.
                            "sources": [{"name": Path(x).stem, "n": None}
                                        for x in (pr["sources"] or [])]
                                       or [{"name": k, "n": v} for k, v
                                           in sorted(pr["by_role"].items())],
                            "by_role": pr["by_role"],
                            "where": "store", "stored": pr["stored"],
                        })
            except Exception as exc:
                self.cfg["maps_error"] = str(exc)[:200]

        cache = self.cfg.setdefault("map_summaries", {})
        for path in reference_map.discover(self.cfg["reference_dir"]):
            if path.stem in seen:
                continue
            key = (str(path), path.stat().st_mtime)
            if key not in cache:
                try:
                    ref = self._reference(path.stem)
                    cache[key] = (reference_map.summarise(path, ref)
                                  if ref else None)
                except Exception as exc:
                    cache[key] = {"name": path.stem, "path": str(path),
                                  "error": str(exc)[:200]}
            if cache[key]:
                out.append({**cache[key], "where": "file"})
        return out

    def _reference(self, name: str | None = None):
        """One layout by name, loaded once and kept. None when there is none.

        `name` is a file stem, so it is whatever the chooser listed. Falling
        back to the default rather than erroring keeps a stale bookmark or a
        deleted map from breaking the page.
        """
        loaded = self.cfg.setdefault("references", {})
        chosen = self._map_path(name)
        if chosen is None:
            return None
        key = str(chosen)
        if key not in loaded:
            try:
                loaded[key] = load_reference(chosen)
            except Exception as exc:
                loaded[key] = None
                self.cfg.setdefault("reference_errors", {})[key] = str(exc)
        return loaded[key]

    def _map_path(self, name: str | None) -> Path | None:
        """Resolve a layout name to a file, or the default when it does not."""
        if name:
            for p in reference_map.discover(self.cfg["reference_dir"]):
                if p.stem == name:
                    return p
        default = Path(self.cfg["reference_path"])
        return default if default.is_file() else None

    def _features(self, run_id: str):
        """A sample's feature hits as (sample, gene ids, L2-normalised matrix).

        Straight out of the store now - the hits are a table, so there is no
        file to find and no stage name to know. ValueError means "nothing to
        plot" and is the caller's 404; ImportError is let through so a missing
        dependency is reported as one rather than as a bad request.
        """
        sample = self._sample(run_id)
        with lake.read(self.cfg["lake"], budget=5) as con:
            t = entities.select(con, "feature_hit",
                                columns=["gene_id", "feature_id", "activation"],
                                sample=sample)
        if t.num_rows == 0:
            raise ValueError(f"{sample} has nothing to project")
        ids, m = sparse_from_table(t)
        return sample, ids, normalise(m)

    def _meta(self, sample: str):
        """(column descriptors, gene_id -> row) for one sample, cached.

        Colouring and filtering the map read every gene column, so this is the
        whole gene level for that sample. Invalidated on the store's latest
        snapshot rather than on file mtimes: one number, and it moves exactly
        when something was written.
        """
        cache = self.cfg.setdefault("meta_cache", {})
        try:
            with lake.read(self.cfg["lake"], budget=5) as con:
                stamp = lake.snapshot_id(con)
                hit = cache.get(sample)
                if hit and hit[0] == stamp:
                    return hit[1], hit[2]
                descs, rows = metadata.load(con, sample)
        except Exception:
            return [], {}
        cache[sample] = (stamp, descs, rows)
        return descs, rows

    def _color_choice(self, requested: str, usable: list, comparing: bool):
        """Which column the plot is coloured by, and its fixed domain.

        An overlay defaults to colouring by sample, because telling the two
        runs apart is the question it was opened to answer. A single run
        defaults to the triage class, which is what the map has always shown.
        A request naming a column this run does not have falls back rather
        than failing - switching runs should not 400.
        """
        # Only a column worth colouring by. `columns` in the response is the
        # whole list, because filtering on an identifier is useful even where
        # colouring by it is not - one colour per protein says nothing.
        by_name = {d["name"]: d for d in usable if d.get("usable")}
        if comparing and requested in ("", "sample"):
            return None, None, None
        want = requested or ("category" if "category" in by_name else "")
        d = by_name.get(want)
        if d is None:
            return None, None, None
        return d["name"], metadata.color_domain(d), d

    def _layout(self, mode, ids, m, b_ids, b_m, comparing):
        """Place one or two samples in 2D. Called under _PROJECT_LOCK.

        Returns (xy, b_xy, mode_used, context, ref), or a (status, message)
        pair for the caller to report - the lock is released either way, and a
        2-tuple is the one shape the success return can never be mistaken for.
        """
        from scipy.sparse import vstack

        ref = self._reference(map_name) if mode == "reference" else None
        context: list = []
        try:
            if ref is not None:
                # Fixed layout: place these proteins in the corpus's space, so
                # coordinates mean the same thing across runs and re-runs.
                xy = _transform(ref, m)
                b_xy = _transform(ref, b_m) if comparing else None
                used = "reference"
                step = max(1, len(ref["ids"]) // MAX_CONTEXT_POINTS)
                context = [{"x": round(float(x), 3), "y": round(float(y), 3)}
                           for x, y in ref["xy"][::step]]
            elif comparing:
                # No shared layout to borrow, so fit one over exactly these two
                # samples. That is internally comparable - the two sit in one
                # space - and comparable with nothing else.
                if len(ids) + len(b_ids) < 4:
                    return (422, "need at least 4 proteins across the two runs "
                                 "to fit a layout")
                both = _fit_layout(vstack([m, b_m]))
                xy, b_xy = both[:len(ids)], both[len(ids):]
                used = "pair"
            else:
                if len(ids) < 4:
                    return (422, "need at least 4 proteins to fit a layout; "
                                 "build a reference map instead")
                xy, b_xy = _fit_layout(m), None
                used = "run"
        except ImportError:
            return (501, "umap-learn is not installed in this environment; "
                         "rebuild the image or "
                         "`uv pip install -r requirements.txt`")
        except Exception as exc:
            return (422, f"projection failed: {exc}")
        return xy, b_xy, used, context, ref

    def _snapshot(self):
        """The store's latest snapshot id, or None if it cannot be read."""
        if lake is None:
            return None
        try:
            with lake.read(self.cfg["lake"], budget=3) as con:
                return lake.snapshot_id(con)
        except Exception:
            return None

    def _atlas_data(self, limit: int):
        """The cohort's activations and metadata, thinned. Independent of map.

        Reading half a million activations and building the sparse matrix costs
        the same whichever layout places the result, and the gene metadata is
        identical across all of them - so it is cached once per snapshot rather
        than once per layout. Before this split, holding three layouts meant
        holding three copies of the same 14,000-gene metadata.
        """
        cache = self.cfg.setdefault("atlas_data", {})
        with lake.read(self.cfg["lake"], budget=10) as con:
            stamp = lake.snapshot_id(con)
            hit = cache.get(limit)
            if hit and hit[0] == stamp:
                return (stamp, *hit[1])
            # gene_id is unique only within a sample, so the matrix is keyed on
            # both - otherwise two samples' genes would collapse into one row.
            t = entities.select(
                con, "feature_hit",
                columns=["sample", "gene_id", "feature_id", "activation"])
            descs, rows = metadata.load_cohort(con)

        if t.num_rows == 0:
            raise ValueError("no embeddings in the store yet")

        import pyarrow as pa

        keyed = t.append_column("key", pa.array(
            [metadata.cohort_key(s_, g) for s_, g in
             zip(t.column("sample").to_pylist(), t.column("gene_id").to_pylist())]))
        ids, m = sparse_from_table(
            keyed.select(["key", "feature_id", "activation"])
                 .rename_columns(["gene_id", "feature_id", "activation"]))
        m = normalise(m)
        # Thin per sample rather than over the whole cohort, so a small sample
        # is not rounded away by a large one - the point of the view is that
        # every sample is on it.
        ids, m, thinned = _cap_per_sample(ids, m, limit)

        cache.clear()                       # one snapshot's worth is enough
        cache[limit] = (stamp, (ids, m, thinned, descs, rows))
        return stamp, ids, m, thinned, descs, rows

    def _stored_xy(self, map_name: str, ids):
        """Coordinates for these points from the store, if it has them.

        A projection kept in the lake is a table, so placing a point is a
        lookup rather than a UMAP transform - no minute of waiting, and no
        pickled estimator needed at all. Returns None when the store has no
        such projection or when it names none of these points, which is what
        happens to a layout fitted against a cohort that has since been renamed.
        """
        if lake is None or not map_name or map_name == "fit":
            return None
        try:
            with lake.read(self.cfg["lake"], budget=5) as con:
                known = {r[0] for r in con.execute(
                    "SELECT name FROM projection_meta").fetchall()}
                if map_name not in known:
                    return None
                found = lake.projection_xy(con, map_name)
        except Exception:
            return None
        if not found:
            return None
        import numpy as np

        xy, miss = np.zeros((len(ids), 2), dtype=np.float32), 0
        for i, key in enumerate(ids):
            sample, _, gene = key.partition(metadata.KEY_SEP)
            hit = found.get((sample, gene))
            if hit is None:
                miss += 1
            else:
                xy[i] = (hit[0], hit[1])
        if miss == len(ids):
            return None                      # names nothing we are drawing
        return xy, miss

    def _atlas_layout(self, limit: int, map_name: str | None = None):
        """Coordinates for one layout. Only this part differs between maps."""
        stamp, ids, m, thinned, descs, rows = self._atlas_data(limit)
        cache = self.cfg.setdefault("atlas_xy", {})
        ckey = (limit, map_name or "")
        hit = cache.get(ckey)
        missing = 0
        if hit and hit[0] == stamp:
            xy, mode, missing = hit[1]
        else:
            stored = self._stored_xy(map_name or "", ids)
            if stored is not None:
                xy, missing = stored
                mode = "stored"
            else:
                ref = None if map_name == "fit" else self._reference(map_name)
                try:
                    if ref is not None:
                        xy, mode = _transform(ref, m), "reference"
                    else:
                        xy, mode = _fit_layout(m), "cohort"
                except Exception:
                    xy, mode = _fit_layout(m), "cohort"
            # Keep several. Comparing layouts is the reason a chooser exists,
            # and re-transforming 30,000 points into a dense map costs a
            # minute - paying that on every toggle would make the comparison
            # not worth making. Coordinates only, so each entry is small.
            cache[ckey] = (stamp, (xy, mode, missing))
            for old in list(cache)[:-ATLAS_CACHE_KEEP]:
                cache.pop(old, None)

        return {"ids": ids, "xy": xy, "mode": mode, "descs": descs,
                "rows": rows, "thinned": thinned, "total": len(ids),
                "stamp": stamp, "map": map_name or "", "missing": missing}

    def _atlas(self, color: str, limit: int, map_name: str = ""):
        """The cohort as a backdrop, with per-sample membership carried along.

        Everything is sent once and hovering is done in the browser: the
        highlight is a restyle, not a request, so it is instant and the layout
        never moves under the cursor.
        """
        if entities is None:
            return self._err(503, "the pipeline registry is unavailable")
        try:
            built = self._atlas_layout(limit, map_name or None)
        except ValueError as exc:
            return self._err(404, str(exc))
        except ImportError as exc:
            return self._err(501, f"the atlas needs scipy and umap-learn: {exc}")
        except Exception as exc:
            return self._err(422, f"cannot build the atlas: {exc}")

        usable = built["descs"]
        column, domain, _ = self._color_choice(color, usable, False)
        order: dict[str, int] = {}
        points = []
        for i, key in enumerate(built["ids"]):
            sample, _, gene = key.partition(metadata.KEY_SEP)
            si = order.setdefault(sample, len(order))
            p = {"x": round(float(built["xy"][i][0]), 2),
                 "y": round(float(built["xy"][i][1]), 2), "i": si}
            if column:
                v = (built["rows"].get(key) or {}).get(column)
                p["v"] = v
                p["s"] = metadata.slot_of(v, domain) if domain else None
            points.append(p)

        samples = [s for s, _ in sorted(order.items(), key=lambda kv: kv[1])]
        counts = [0] * len(samples)
        for p in points:
            counts[p["i"]] += 1
        return self._json({
            "snapshot": built["stamp"], "map": built.get("map") or "",
            "mode": built["mode"], "n": len(points), "points": points,
            "samples": samples, "counts": counts,
            "color": column, "domain": domain, "columns": usable,
            "thinned": built["thinned"],
            # Points the chosen layout has no stored coordinates for - a sample
            # embedded after it was fitted. Reported rather than hidden: they
            # would otherwise be drawn at the origin and read as a real cluster.
            "missing": built.get("missing", 0),
            # Echoed so a selection can name the layout it was made on. The
            # cap is what decides which proteins are in `points` at all, so a
            # selection resolved against a different one would resolve to
            # different proteins.
            "limit": limit,
        })

    def _gene_seqs(self, keys: list[str]) -> dict[str, str]:
        """key -> amino acid sequence, read from the store a chunk at a time.

        ``seq`` is the one gene column the atlas deliberately does not cache:
        it is megabytes of amino acids that nothing on the plot can show, and
        it is wanted only for the handful of proteins a selection covers. So
        it is fetched here, per selection, rather than carried by every point.

        Keyed by sample because ``gene_id`` is unique only within one - two
        samples can both hold a ``gene_1``, and a query that forgot the sample
        would hand back the wrong protein's sequence.
        """
        by_sample: dict[str, list[str]] = {}
        for k in keys:
            sample, _, gene = k.partition(metadata.KEY_SEP)
            by_sample.setdefault(sample, []).append(gene)
        out: dict[str, str] = {}
        with lake.read(self.cfg["lake"], budget=10) as con:
            for sample, gids in by_sample.items():
                for i in range(0, len(gids), SELECT_CHUNK):
                    lits = ", ".join("'" + g.replace("'", "''") + "'"
                                     for g in gids[i:i + SELECT_CHUNK])
                    t = entities.select(con, "gene", where=f"gene_id IN ({lits})",
                                        columns=["gene_id", "seq"], sample=sample)
                    for gid, seq in zip(t.column("gene_id").to_pylist(),
                                        t.column("seq").to_pylist()):
                        if seq:
                            out[metadata.cohort_key(sample, gid)] = seq
        return out

    def _selection(self, fmt: str):
        """Picked atlas points, resolved to the rows and sequences behind them.

        Two shapes from one route, because they answer the same question at
        two sizes: ``json`` is what the panel shows - a couple of hundred rows
        with every column the gene level carries - and ``fasta`` is the whole
        selection as a file. Splitting them into two routes would have meant
        two copies of the index-to-protein resolution, which is the only part
        that can go subtly wrong.
        """
        if entities is None:
            return self._err(503, "the pipeline registry is unavailable")
        try:
            length = int(self.headers.get("Content-Length", "0"))
            body = json.loads(self.rfile.read(length) or b"{}")
        except (ValueError, json.JSONDecodeError):
            return self._err(400, "body must be JSON")
        if not isinstance(body, dict):
            return self._err(400, "body must be a JSON object")
        raw = body.get("indices")
        if not isinstance(raw, list) or not raw:
            return self._err(400, "indices must be a non-empty array")
        if len(raw) > MAX_ATLAS_POINTS:
            return self._err(413, "more indices than the atlas has points")
        try:
            limit = int(body.get("limit") or MAX_ATLAS_POINTS)
        except (TypeError, ValueError):
            limit = MAX_ATLAS_POINTS
        limit = max(500, min(limit, MAX_ATLAS_POINTS))

        try:
            stamp, ids, _m, thinned, descs, rows = self._atlas_data(limit)
        except ValueError as exc:
            return self._err(404, str(exc))
        except ImportError as exc:
            return self._err(501, f"the atlas needs scipy and umap-learn: {exc}")
        except Exception as exc:
            return self._err(422, f"cannot read the atlas: {exc}")

        # A selection is only meaningful against the point set it was drawn
        # on. Refusing a stale one is the whole reason the snapshot travels
        # with it: silently answering would hand back other proteins.
        want = body.get("snapshot")
        if want is not None and want != stamp:
            return self._err(409, "the store has changed since this atlas was "
                                  "drawn - rebuild it and select again")

        keys, seen = [], set()
        for i in raw:
            try:
                j = int(i)
            except (TypeError, ValueError):
                return self._err(400, "indices must be integers")
            if not 0 <= j < len(ids):
                return self._err(400, "an index is outside this atlas")
            if ids[j] not in seen:                # a box may be added twice
                seen.add(ids[j])
                keys.append(ids[j])

        # Over every picked point, not just the ones sent back: a box that
        # covers 4,000 proteins should say which samples they came from even
        # though the panel shows 200 of them.
        tally: dict[str, int] = {}
        for k in keys:
            name = k.partition(metadata.KEY_SEP)[0]
            tally[name] = tally.get(name, 0) + 1
        samples = [{"name": s, "n": n} for s, n in
                   sorted(tally.items(), key=lambda kv: (-kv[1], kv[0]))]

        cap = MAX_SELECTION_FASTA if fmt == "fasta" else MAX_SELECTION_SHOWN
        shown = keys[:cap]
        try:
            seqs = self._gene_seqs(shown)
        except Exception as exc:
            return self._err(503, f"cannot read the store: {exc}")

        # Every column the gene level carries, in the store's own order, plus
        # the two the atlas drops from its cache - `sample`, which the point
        # carried instead, and `seq`, which is what all of this is for. No
        # column is named here, so a stage added tomorrow is another column in
        # the table and another field in the FASTA header with no change.
        names = [d["name"] for d in descs]
        out = []
        for k in shown:
            sample, _, gene = k.partition(metadata.KEY_SEP)
            row = dict(rows.get(k) or {})
            row.pop("sample", None)
            r = {"sample": sample, **{n: row.get(n) for n in names},
                 "seq": seqs.get(k)}
            if r.get("gene_id") is None:       # no gene row behind the point
                r["gene_id"] = gene
            out.append(r)
        missing = sum(1 for r in out if not r["seq"])

        text = "".join(_fasta_record(r, names) for r in out if r["seq"])
        if fmt == "fasta":
            return self._send(200, text.encode(), "text/x-fasta; charset=utf-8",
                              filename=f"atlas-selection-{len(out) - missing}.faa")
        return self._json({
            "snapshot": stamp, "limit": limit, "n": len(keys),
            "shown": len(out), "truncated": len(keys) > len(out),
            "cap": cap, "fasta_cap": MAX_SELECTION_FASTA,
            "thinned": thinned, "missing": missing, "samples": samples,
            "columns": ["sample", *names, "seq"], "rows": out,
            # The same records the download would carry, for the rows sent.
            # Rendered once, here, so what the page shows and what lands on
            # disk cannot drift apart - the alternative was a second FASTA
            # writer in JavaScript that had to agree with this one.
            "fasta": text,
        })

    def _projection(self, run_id: str, mode: str, b_id: str = "",
                    color: str = "", filt: str = "", map_name: str = ""):
        """One run's proteins in 2D, or two runs drawn in the same 2D.

        The overlay is the reason the layout is fixed rather than fitted per
        run: coordinates that mean the same thing across runs are what make two
        samples drawn together readable rather than decorative.
        """
        if mode not in ("reference", "run"):
            return self._err(400, "mode must be reference or run")
        if b_id and b_id == run_id:
            return self._err(400, "cannot compare a run with itself")
        try:
            sample_a, ids, m = self._features(run_id)
            other = self._features(b_id) if b_id else None
        except ImportError as exc:
            return self._err(501, f"projection needs scipy and pyarrow: {exc}")
        except ValueError as exc:
            return self._err(404, str(exc))
        except Exception as exc:
            return self._err(422, f"cannot read features: {exc}")
        if len(ids) < 2 or (other is not None and len(other[1]) < 2):
            return self._err(422, "need at least 2 proteins in each run")

        # Two samples can double the point count, so cap what is drawn - by a
        # stride rather than a random draw, so the same pair always yields the
        # same picture instead of reshuffling on every request.
        limit = MAX_PROJECTION_POINTS // (2 if other else 1)
        ids, m, cut_a = _cap(ids, m, limit)
        cut_b = False
        if other is not None:
            sample_b, b_ids, b_m = other
            b_ids, b_m, cut_b = _cap(b_ids, b_m, limit)

        from scipy.sparse import vstack

        # Serialised: numba aborts the process if two threads call into it at
        # once, and this server hands every request its own thread.
        with _PROJECT_LOCK:
            laid = self._layout(mode, ids, m,
                                b_ids if other is not None else None,
                                b_m if other is not None else None,
                                other is not None)
        if len(laid) == 2:                            # (status, message)
            return self._err(*laid)
        xy, b_xy, used, context, ref = laid

        # Colour and filter by any column the gene level carries. Which columns
        # exist is read from the run, so a stage added later is another thing
        # to colour by with no change here or in the page.
        descs_a, rows_a = self._meta(sample_a)
        descs_b, rows_b = self._meta(sample_b) if other is not None else ([], {})
        # An overlay may only offer columns both sides have: a scale shown over
        # two runs has to mean the same thing on both.
        usable = metadata.merge(descs_a, descs_b) if other is not None else descs_a
        column, domain, _ = self._color_choice(color, usable, other is not None)

        keep_a = keep_b = None
        where = ""
        if filt:
            try:
                terms = json.loads(filt)
            except json.JSONDecodeError:
                return self._err(400, "filter must be a JSON array of terms")
            try:
                where = metadata.build_where(terms, usable)
            except metadata.BadFilter as exc:
                return self._err(400, str(exc))
            if where:
                keep_a = metadata.matching(self.cfg["lake"], sample_a, where)
                if other is not None:
                    keep_b = metadata.matching(self.cfg["lake"], sample_b, where)

        points = _points(ids, xy, rows_a, column, domain, keep_a,
                         "a" if other is not None else None)
        shown_a = len(points)
        if other is not None:
            points += _points(b_ids, b_xy, rows_b, column, domain, keep_b, "b")
        shown_b = len(points) - shown_a
        if not points:
            return self._err(422, "the filter matched nothing in this run")

        # With top-K over a 16,384-wide codebook two proteins may share no
        # features at all, and then the layout is noise. Say so rather than let
        # it be read as biology.
        import numpy as np
        allm = vstack([m, b_m]) if other is not None else m
        counts = np.bincount(allm.tocoo().col, minlength=16384)
        distinct = int((counts > 0).sum())
        shared = int((counts > 1).sum())
        body = {
            "run": run_id, "mode": used, "n": len(points), "points": points,
            # `n` is what is drawn and `total` what was projected: a filter is
            # unreadable without both.
            "total": len(ids) + (len(b_ids) if other is not None else 0),
            "context": context,
            # What the plot is coloured by, the columns it could be coloured or
            # filtered by instead, and the filter actually applied.
            "color": column, "domain": domain, "columns": usable,
            "filter": where or None,
            "distinct_features": distinct, "shared_features": shared,
            "shared_frac": round(shared / distinct, 3) if distinct else 0.0,
        }
        if other is not None:
            body["samples"] = [
                {"side": "a", "run": run_id, "sample": sample_a,
                 "n": len(ids), "shown": shown_a, "subsampled": cut_a},
                {"side": "b", "run": b_id, "sample": sample_b,
                 "n": len(b_ids), "shown": shown_b, "subsampled": cut_b},
            ]
        elif cut_a:
            body["subsampled"] = True
        if ref is not None:
            body["reference"] = {
                "n": ref["n"], "built": ref["built"],
                "inputs": [Path(i).name for i in ref["inputs"]],
                "version_drift": ref.get("version_drift"),
            }
        elif self.cfg.get("reference_error"):
            body["reference_error"] = self.cfg["reference_error"]
        return self._json(body)

    def _preview(self, run_id: str, stage: str, name: str, limit: int):
        try:
            run_dir = self._run_dir(run_id)
        except ValueError as exc:
            return self._err(404, str(exc))
        if stage not in stage_dirs(run_dir):
            return self._err(400, "unknown stage")
        if "/" in name or "\\" in name or not SAFE_NAME.match(name):
            return self._err(400, "bad filename")
        p = (run_dir / stage / name).resolve()
        if run_dir not in p.parents or not p.is_file():
            return self._err(404, "no such artifact")
        limit = max(1, min(limit, 5000))
        try:
            body = preview(p, limit)
        except Exception as exc:                 # a corrupt artifact is data,
            return self._err(422, f"cannot preview: {exc}")   # not a crash
        body.update({"name": name, "stage": stage,
                     "bytes": p.stat().st_size, "limit": limit})
        return self._json(body)

    def _upload(self, name: str):
        if not self._guard_writes():
            return
        if "/" in name or "\\" in name or not SAFE_NAME.match(name):
            return self._err(400, "filename must be simple: letters, digits, . _ -")
        try:
            length = int(self.headers.get("Content-Length", "0"))
        except ValueError:
            return self._err(400, "bad Content-Length")
        if length <= 0:
            return self._err(400, "empty upload")
        if length > MAX_UPLOAD:
            return self._err(413, f"too large (max {MAX_UPLOAD // 1024**3} GiB)")

        uploads = self.cfg["uploads"]
        try:
            uploads.mkdir(parents=True, exist_ok=True)
        except OSError as exc:
            return self._err(500, f"cannot create {uploads}: {exc}")
        dest = uploads / name
        tmp = dest.with_name(dest.name + ".part")
        remaining = length
        try:
            with open(tmp, "wb") as fh:
                while remaining > 0:
                    chunk = self.rfile.read(min(1 << 20, remaining))
                    if not chunk:
                        break
                    fh.write(chunk)
                    remaining -= len(chunk)
            if remaining:
                tmp.unlink(missing_ok=True)
                return self._err(400, "upload truncated")
            tmp.replace(dest)
        except OSError as exc:
            tmp.unlink(missing_ok=True)
            return self._err(500, f"write failed: {exc}")
        self._json({"path": str(dest), "name": name, "bytes": length})

    def _container_path(self, host: str) -> str | None:
        """Host path -> the path run.sh binds it to inside the image."""
        p = Path(host).resolve()
        for root, inside in self.cfg["binds"]:
            try:
                rel = p.relative_to(Path(root).resolve())
            except ValueError:
                continue
            return str(Path(inside) / rel) if str(rel) != "." else inside
        return None

    def _job_argv(self, sample: str, argv_tail: list[str]) -> tuple[list[str], str]:
        """What the scheduler should actually run.

        A submitted job lands on a compute node, which has the image but not
        this server's interpreter — so it runs container/run.sh, and any path
        in the arguments is rewritten to the path the image sees.
        """
        if self.cfg["runner"] == "python":
            return ([self.cfg["python"], str(PIPELINE / "run.py"),
                     "--work", str(self.cfg["roots"][0]),
                     "--sample", sample, *argv_tail], "python")

        out: list[str] = []
        skip = False
        for i, tok in enumerate(argv_tail):
            if skip:
                skip = False
                continue
            if tok in ("--fastq", "--fastq2", "--contigs", "--proteins",
                       "--hmm", "--ref"):
                inside = self._container_path(argv_tail[i + 1])
                if inside is None:
                    raise ValueError(
                        f"{argv_tail[i + 1]} is not under a directory the "
                        f"container can see ({', '.join(b for _, b in self.cfg['binds'])})")
                out += [tok, inside]
                skip = True
            else:
                out.append(tok)
        # Pin run.sh's binds to what this server is actually showing.
        # Its defaults are relative to $PWD, so a job would otherwise write
        # into whatever directory it happened to start in rather than the work
        # root the UI lists.
        env = [f"{var}={path}" for var, path in (
            ("SAE_WORK", self.cfg["roots"][0]),
            ("SAE_DATA", self.cfg["data"][0]),
            ("SAE_ATLAS", REPO / "sae"),
        )]
        # run.sh's pipeline app supplies --work /work itself.
        return (["env", *env, str(REPO / "container" / "run.sh"), "pipeline",
                 "--sample", sample, *out], "container")

    def _launch(self):
        if not self._guard_writes():
            return
        try:
            length = int(self.headers.get("Content-Length", "0"))
            body = json.loads(self.rfile.read(length) or b"{}")
        except (ValueError, json.JSONDecodeError):
            return self._err(400, "body must be JSON")

        sample = str(body.get("sample", "")).strip()
        if not SAFE_NAME.match(sample):
            return self._err(400, "sample must be letters, digits, . _ - (max 64)")
        try:
            argv_tail = build_argv(body, self.cfg["allowed"])
        except ValueError as exc:
            return self._err(400, str(exc))

        try:
            argv, runner = self._job_argv(sample, argv_tail)
        except ValueError as exc:
            return self._err(400, str(exc))
        if runner == "python" and not Path(self.cfg["python"]).is_file():
            return self._err(500, f"interpreter not found: {self.cfg['python']}")
        try:
            job = self.cfg["jobs"].submit(
                argv, sample, REPO if runner == "container" else PIPELINE)
        except (OSError, RuntimeError) as exc:
            return self._err(500, f"submission failed: {exc}")
        self._json(self.cfg["jobs"]._as_dict(job), 201)


def main():
    d = defaults()
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--work", type=Path, action="append",
                   help=f"work root to display; repeatable. Runs launch into "
                        f"the first. Default: {', '.join(str(x) for x in d['work'])}")
    p.add_argument("--data", type=Path, action="append",
                   help=f"directory offered as run input. Default: "
                        f"{', '.join(str(x) for x in d['data'])}")
    p.add_argument("--uploads", type=Path, default=d["uploads"])
    p.add_argument("--lake", default=None,
                   help="the store every run writes to. A path, or "
                        "postgres:/sqlite: for a shared catalog; SAE_LAKE_DATA "
                        "points the parquet elsewhere, including s3://")
    p.add_argument("--reference", type=Path, default=d["reference"],
                   help="fixed UMAP layout from reference_map.py; without one, "
                        "each run is projected on its own and coordinates are "
                        "not comparable between runs")
    p.add_argument("--python", default=d["python"])
    p.add_argument("--host", default=d["host"])
    p.add_argument("--port", type=int, default=8765)
    p.add_argument("--job-cpus", type=int, default=4,
                   help="exported as OMP_NUM_THREADS to each job")
    p.add_argument("--job-time", default="08:00:00",
                   help="wall clock per job; over it the job is terminated")
    p.add_argument("--job-runner", choices=["auto", "python", "container"],
                   default="auto",
                   help="what a job runs. 'container' runs container/run.sh "
                        "so the job needs only the image; 'python' runs this "
                        "server's interpreter. auto picks container when a "
                        ".sif is present")
    p.add_argument("--read-only", action="store_true",
                   help="serve progress only; reject upload and launch")
    p.add_argument("--published", action="store_true",
                   help="something in front of this process controls exposure "
                        "(run.sh passes it when publishing a container port), "
                        "so a non-loopback bind is intentional")
    p.add_argument("--verbose", action="store_true")
    a = p.parse_args()

    sif = Path(os.environ.get("SAE_SIF", REPO / "container" / "sae.sif"))
    runner = a.job_runner
    if runner == "auto":
        runner = "container" if sif.is_file() else "python"

    roots = [r.resolve() for r in (a.work or d["work"])]
    data = [r.resolve() for r in (a.data or d["data"])]
    uploads = a.uploads.resolve()
    Handler.cfg = {
        "roots": roots, "data": data, "uploads": uploads, "python": a.python,
        "reference_path": a.reference,
        # Layouts are discovered here rather than listed, so building one makes
        # it selectable.
        "reference_dir": Path(a.reference).parent,
        "lake": a.lake or lake.default_target(REPO) if lake else None,
        "runner": runner,
        # Mirrors run.sh's bind table, so a host path can be rewritten to the
        # path a job sees inside the image.
        "binds": [(roots[0], "/work"), (data[0], "/data"), (REPO / "sae", "/atlas")],
        "read_only": a.read_only, "verbose": a.verbose,
        # A run may only read from these; see checked_path.
        "allowed": [uploads, *data, *roots],
        "jobs": Jobs(make_launcher(uploads / ".logs", cpus=a.job_cpus,
                                   time_limit=a.job_time),
                     uploads / ".logs"),
    }

    where = f"{'docker' if IN_DOCKER else 'apptainer'} container" \
        if IN_CONTAINER else "host"
    print(f"  environment: {where}")
    print(f"  work roots : {', '.join(str(r) for r in roots)}")
    print(f"  data       : {', '.join(str(r) for r in data)}")
    print(f"  uploads    : {uploads}")
    print(f"  lake       : {Handler.cfg['lake']}")
    print(f"  interpreter: {a.python}")
    print(f"  job runner : {runner}"
          + ("" if runner == "container" else
             "  — needs this interpreter visible wherever the job runs"))
    _l = Handler.cfg["jobs"].launcher.describe()
    print(f"  launcher   : {_l['backend']} — one run at a time "
          f"({_l['cpus']} cpus, {_l['time']})")
    if a.read_only:
        print("  mode       : read-only (upload and launch disabled)")
    # A non-loopback bind is only alarming when this process is what decides
    # reachability. Nested in a container whose port the runtime publishes to
    # host loopback, 0.0.0.0 is the container's own namespace and is correct.
    if a.host not in ("127.0.0.1", "localhost", "::1") and not a.published:
        print(f"\n  WARNING: bound to {a.host}, not localhost. This server "
              f"launches\n           subprocesses and accepts uploads. Do not "
              f"expose it.\n")
    if a.host == "0.0.0.0":
        # Bound to every interface, which means the reachable address belongs
        # to whatever published the port (run.sh), not to us. Naming a URL here
        # would name the wrong one.
        via = " (published by the runtime)" if a.published else ""
        print(f"\n  listening on {a.host}:{a.port}{via}\n", flush=True)
    else:
        print(f"\n  http://{a.host}:{a.port}\n", flush=True)
    try:
        ThreadingHTTPServer((a.host, a.port), Handler).serve_forever()
    except KeyboardInterrupt:
        print("\nstopped")


if __name__ == "__main__":
    main()
