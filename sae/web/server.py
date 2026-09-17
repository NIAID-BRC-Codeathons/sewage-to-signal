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

# Single source of truth for the stage list; fall back if run.py cannot import.
try:
    sys.path.insert(0, str(PIPELINE))
    from run import ORDER as STAGES          # type: ignore
except Exception:                            # pragma: no cover - defensive
    STAGES = ["s01_qc", "s02_assemble", "s03_genes", "s04_derep",
              "s05_prefilter", "s06_embed", "s07_match"]

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

# Options forwarded to run.py. Anything not listed here is rejected, so the
# request body can never introduce a new flag.
OPTIONS: dict[str, type] = {
    "from": str, "to": str, "max_reads": int, "min_aa": int,
    "hmm": str, "ref": str, "evalue": float, "confident_evalue": float,
    "min_coverage": float, "bit_cutoffs": str, "top_k": int,
    "batch_size": int, "limit": int, "device": str, "force": bool,
    "model": str, "max_len": int,
}
CHOICES = {"from": STAGES, "to": STAGES,
           "bit_cutoffs": ["gathering", "noise", "trusted"],
           "model": ["6b", "300m"]}


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
from reference_map import load as load_reference          # noqa: E402
from reference_map import normalise, sparse_features      # noqa: E402


def _groups_from_classification(tsv: Path) -> dict[str, str]:
    """gene_id -> s05 category, when the run has one. Optional by design."""
    out: dict[str, str] = {}
    try:
        with open(tsv) as fh:
            header = fh.readline().rstrip("\n").split("\t")
            gi, ci = header.index("gene_id"), header.index("category")
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) > max(gi, ci):
                    out[parts[gi]] = parts[ci]
    except (OSError, ValueError):
        pass
    return out


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
    }


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
            for stage in STAGES:
                sd = sample_dir / stage
                if not sd.is_dir():
                    continue
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
    if last not in STAGES:
        return last
    # The line is printed on completion, so the next stage is in flight.
    i = STAGES.index(last)
    return STAGES[i + 1] if (live and i + 1 < len(STAGES)) else last


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
    """Translate a request body into run.py flags, rejecting anything unlisted."""
    kind = body.get("input_kind")
    if kind not in ("fastq", "contigs", "proteins"):
        raise ValueError("input_kind must be fastq, contigs or proteins")
    argv = [f"--{kind}", checked_path(body.get("input_path"), allowed)]
    if kind == "fastq" and body.get("input_path2"):
        argv += ["--fastq2", checked_path(body["input_path2"], allowed)]

    for key, value in (body.get("options") or {}).items():
        if key not in OPTIONS:
            raise ValueError(f"unknown option: {key}")
        if value in (None, ""):
            continue
        typ = OPTIONS[key]
        flag = "--" + key.replace("_", "-")
        if typ is bool:
            if value:
                argv.append(flag)
            continue
        if key in ("hmm", "ref"):          # these are paths too
            argv += [flag, checked_path(value, allowed)]
            continue
        try:
            coerced = typ(value)
        except (TypeError, ValueError):
            raise ValueError(f"{key} must be {typ.__name__}")
        if key in CHOICES and str(coerced) not in CHOICES[key]:
            raise ValueError(f"{key} must be one of {', '.join(CHOICES[key])}")
        argv += [flag, str(coerced)]
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
    def _send(self, code: int, body: bytes, ctype: str):
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
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
                "runs": scan(self.cfg["roots"]),
                "jobs": self.cfg["jobs"].all(),
                "stages": STAGES,
                "roots": [str(r) for r in self.cfg["roots"]],
                "read_only": self.cfg["read_only"],
                "uploads": str(self.cfg["uploads"]),
                "hmms": self._candidate_hmms(),
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
        if u.path == "/api/inputs":
            return self._json({"files": self._candidate_inputs()})
        if u.path == "/api/artifacts":
            return self._artifacts((q.get("run") or [""])[0])
        if u.path == "/api/projection":
            return self._projection((q.get("run") or [""])[0],
                                    (q.get("mode") or ["reference"])[0])
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
    def _candidate_hmms(self) -> list[dict]:
        """Profile databases for s05. Without one the stage is a pass-through
        and every protein goes to the GPU, so the UI should not make finding
        it a matter of knowing a path."""
        out = []
        for root in self.cfg["data"]:
            if not root.is_dir():
                continue
            for p in sorted(root.rglob("*.hmm")):
                if p.is_file() and not p.name.startswith("."):
                    out.append({"path": str(p), "name": p.name,
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

    def _artifacts(self, run_id: str):
        try:
            run_dir = self._run_dir(run_id)
        except ValueError as exc:
            return self._err(404, str(exc))
        out: dict[str, list] = {}
        for stage in STAGES:
            sd = run_dir / stage
            if not sd.is_dir():
                continue
            files = []
            for f in sorted(sd.iterdir()):
                if not f.is_file() or f.name.startswith("."):
                    continue
                files.append({"name": f.name, "bytes": f.stat().st_size,
                              "kind": artifact_kind(f)})
            if files:
                out[stage] = files
        return self._json({"run": run_id, "stages": out})

    def _reference(self):
        """The shared layout, loaded once. None when there is no map yet."""
        if "reference" not in self.cfg:
            try:
                self.cfg["reference"] = load_reference(self.cfg["reference_path"])
            except Exception as exc:
                self.cfg["reference"] = None
                self.cfg["reference_error"] = str(exc)
        return self.cfg["reference"]

    def _projection(self, run_id: str, mode: str):
        if mode not in ("reference", "run"):
            return self._err(400, "mode must be reference or run")
        try:
            run_dir = self._run_dir(run_id)
        except ValueError as exc:
            return self._err(404, str(exc))

        found = sorted((run_dir / "s06_embed").glob("*.sae_features.parquet"))
        if not found:
            return self._err(404, "this run has no s06_embed output")
        try:
            ids, m = sparse_features(found[0])
        except ImportError as exc:
            return self._err(501, f"projection needs scipy and pyarrow: {exc}")
        except Exception as exc:
            return self._err(422, f"cannot read features: {exc}")
        if len(ids) < 2:
            return self._err(422, f"need at least 2 proteins, got {len(ids)}")

        m = normalise(m)
        ref = self._reference() if mode == "reference" else None
        context: list = []
        try:
            if ref is not None:
                # Fixed layout: place these proteins in the corpus's space, so
                # coordinates mean the same thing across runs and re-runs.
                xy = ref["reducer"].transform(m)
                used = "reference"
                step = max(1, len(ref["ids"]) // MAX_CONTEXT_POINTS)
                context = [{"x": round(float(x), 3), "y": round(float(y), 3)}
                           for x, y in ref["xy"][::step]]
            else:
                import umap
                if len(ids) < 4:
                    return self._err(422, "need at least 4 proteins to fit a "
                                          "layout; build a reference map instead")
                xy = umap.UMAP(n_components=2, metric="cosine", random_state=0,
                               n_neighbors=max(2, min(15, len(ids) - 1)),
                               min_dist=0.1).fit_transform(m)
                used = "run"
        except ImportError:
            return self._err(501, "umap-learn is not installed in this "
                                  "environment; rebuild the image or "
                                  "`uv pip install -r requirements.txt`")
        except Exception as exc:
            return self._err(422, f"projection failed: {exc}")

        groups = _groups_from_classification(
            next(iter(sorted((run_dir / "s05_prefilter").glob("*.classification.tsv"))),
                 run_dir / "missing"))
        points = [{"id": g, "x": round(float(xy[i][0]), 3),
                   "y": round(float(xy[i][1]), 3), "group": groups.get(g)}
                  for i, g in enumerate(ids)]

        # With top-K over a 16,384-wide codebook two proteins may share no
        # features at all, and then the layout is noise. Say so rather than let
        # it be read as biology.
        import numpy as np
        counts = np.bincount(m.tocoo().col, minlength=16384)
        distinct = int((counts > 0).sum())
        shared = int((counts > 1).sum())
        body = {
            "run": run_id, "mode": used, "n": len(points), "points": points,
            "groups": sorted({p["group"] for p in points if p["group"]}),
            "context": context,
            "distinct_features": distinct, "shared_features": shared,
            "shared_frac": round(shared / distinct, 3) if distinct else 0.0,
        }
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
        if stage not in STAGES:
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
