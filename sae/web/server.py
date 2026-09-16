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

# Options forwarded to run.py. Anything not listed here is rejected, so the
# request body can never introduce a new flag.
OPTIONS: dict[str, type] = {
    "from": str, "to": str, "max_reads": int, "min_aa": int,
    "hmm": str, "ref": str, "evalue": float, "confident_evalue": float,
    "min_coverage": float, "bit_cutoffs": str, "top_k": int,
    "batch_size": int, "limit": int, "device": str, "force": bool,
}
CHOICES = {"from": STAGES, "to": STAGES,
           "bit_cutoffs": ["gathering", "noise", "trusted"]}


def defaults() -> dict:
    """Container and host disagree about where everything lives."""
    if IN_CONTAINER:
        return {
            "work": [Path("/work")],
            "uploads": Path("/work/uploads"),
            "data": [Path("/data")],
            "python": "/opt/venv/bin/python",
            # Docker isolates the network namespace, so localhost inside is
            # unreachable from the host; run.sh publishes to host loopback.
            "host": "0.0.0.0" if IN_DOCKER else "127.0.0.1",
        }
    return {
        "work": [REPO / "work", PIPELINE / "work"],
        "uploads": REPO / "uploads",
        "data": [REPO / "data"],
        "python": str(REPO / "sae" / ".venv" / "bin" / "python"),
        "host": "127.0.0.1",
    }


# --------------------------------------------------------------------------
# reading pipeline state
# --------------------------------------------------------------------------
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
@dataclass
class Job:
    id: str
    sample: str
    argv: list[str]
    log: Path
    started: float
    proc: subprocess.Popen = field(repr=False)

    def state(self) -> str:
        rc = self.proc.poll()
        if rc is None:
            return "running"
        return "finished" if rc == 0 else "failed"

    def as_dict(self) -> dict:
        return {
            "id": self.id, "sample": self.sample, "state": self.state(),
            "returncode": self.proc.poll(), "started": self.started,
            "elapsed": round(time.time() - self.started, 1),
            "argv": self.argv, "stage": self._stage_from_log(),
        }

    def _stage_from_log(self) -> str | None:
        """run.py prints '[stage] ...' as each stage completes."""
        try:
            text = self.log.read_text(errors="replace")
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
        return STAGES[i + 1] if (self.state() == "running"
                                 and i + 1 < len(STAGES)) else last


class Jobs:
    def __init__(self, log_dir: Path):
        self.log_dir = log_dir
        self._jobs: dict[str, Job] = {}
        self._lock = threading.Lock()

    def launch(self, python: str, work: Path, sample: str,
               argv_tail: list[str]) -> Job:
        self.log_dir.mkdir(parents=True, exist_ok=True)
        jid = uuid.uuid4().hex[:12]
        log = self.log_dir / f"{jid}.log"
        argv = [python, str(PIPELINE / "run.py"),
                "--work", str(work), "--sample", sample, *argv_tail]
        fh = open(log, "w", buffering=1)
        fh.write(f"$ {' '.join(argv)}\n\n")
        proc = subprocess.Popen(argv, cwd=str(PIPELINE), stdout=fh,
                                stderr=subprocess.STDOUT, text=True)
        job = Job(jid, sample, argv, log, time.time(), proc)
        with self._lock:
            self._jobs[jid] = job
        return job

    def all(self) -> list[dict]:
        with self._lock:
            jobs = list(self._jobs.values())
        return sorted((j.as_dict() for j in jobs),
                      key=lambda d: d["started"], reverse=True)

    def get(self, jid: str) -> Job | None:
        with self._lock:
            return self._jobs.get(jid)


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
                "container": "docker" if IN_DOCKER else
                             ("apptainer" if IN_APPTAINER else None),
            })
        if u.path == "/api/log":
            job = self.cfg["jobs"].get((q.get("id") or [""])[0])
            if job is None:
                return self._err(404, "no such job")
            try:
                data = job.log.read_bytes()[-LOG_TAIL:]
            except OSError:
                data = b""
            return self._json({"id": job.id, "state": job.state(),
                               "log": data.decode(errors="replace")})
        if u.path == "/api/inputs":
            return self._json({"files": self._candidate_inputs()})
        return self._err(404, "not found")

    def do_POST(self):
        u = urlparse(self.path)
        q = parse_qs(u.query)
        if u.path == "/api/upload":
            return self._upload((q.get("name") or [""])[0])
        if u.path == "/api/run":
            return self._launch()
        return self._err(404, "not found")

    # -- handlers
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

        python = self.cfg["python"]
        if not Path(python).is_file():
            return self._err(500, f"interpreter not found: {python}")
        try:
            job = self.cfg["jobs"].launch(python, self.cfg["roots"][0],
                                          sample, argv_tail)
        except OSError as exc:
            return self._err(500, f"launch failed: {exc}")
        self._json(job.as_dict(), 201)


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
    p.add_argument("--python", default=d["python"])
    p.add_argument("--host", default=d["host"])
    p.add_argument("--port", type=int, default=8765)
    p.add_argument("--read-only", action="store_true",
                   help="serve progress only; reject upload and launch")
    p.add_argument("--verbose", action="store_true")
    a = p.parse_args()

    roots = [r.resolve() for r in (a.work or d["work"])]
    data = [r.resolve() for r in (a.data or d["data"])]
    uploads = a.uploads.resolve()
    Handler.cfg = {
        "roots": roots, "data": data, "uploads": uploads, "python": a.python,
        "read_only": a.read_only, "verbose": a.verbose,
        # A run may only read from these; see checked_path.
        "allowed": [uploads, *data, *roots],
        "jobs": Jobs(uploads / ".logs"),
    }

    where = f"{'docker' if IN_DOCKER else 'apptainer'} container" \
        if IN_CONTAINER else "host"
    print(f"  environment: {where}")
    print(f"  work roots : {', '.join(str(r) for r in roots)}")
    print(f"  data       : {', '.join(str(r) for r in data)}")
    print(f"  uploads    : {uploads}")
    print(f"  interpreter: {a.python}")
    if a.read_only:
        print("  mode       : read-only (upload and launch disabled)")
    # Inside Docker 0.0.0.0 is the container's own namespace, not the host's.
    if a.host not in ("127.0.0.1", "localhost", "::1") and not IN_DOCKER:
        print(f"\n  WARNING: bound to {a.host}, not localhost. This server "
              f"launches\n           subprocesses and accepts uploads. Do not "
              f"expose it.\n")
    if IN_DOCKER:
        # The published host port is run.sh's business, not ours; printing a
        # localhost URL here would name the container-internal port.
        print(f"\n  listening on {a.host}:{a.port} inside the container\n",
              flush=True)
    else:
        print(f"\n  http://{a.host}:{a.port}\n", flush=True)
    try:
        ThreadingHTTPServer((a.host, a.port), Handler).serve_forever()
    except KeyboardInterrupt:
        print("\nstopped")


if __name__ == "__main__":
    main()
