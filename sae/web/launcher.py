"""How the web UI starts a pipeline run.

The server used to fork the pipeline as its own child. On a cluster that is
wrong twice over: the run dies with the server, and it is confined to whatever
allocation the *UI* was given, so a dashboard sized for browsing cannot start
real work. Runs are submitted to a scheduler instead.

Two backends behind one interface, chosen by whether ``sbatch`` is on PATH:

* ``SlurmLauncher``  - submits with ``sbatch``, polls ``squeue`` then ``sacct``.
* ``LocalLauncher``  - forks, for a machine with no scheduler at all.

The point of keeping both is that the deployment should not change shape
between a laptop and a cluster. Running SLURM locally (``container/slurm-local``)
gets you the scheduler path everywhere, so what you test is what you deploy;
the local fork remains for when even that is unavailable.
"""

from __future__ import annotations

import os
import re
import shlex
import shutil
import signal
import subprocess
import time
import uuid
from dataclasses import dataclass, field
from pathlib import Path

STATE_RUNNING = "running"
STATE_FINISHED = "finished"
STATE_FAILED = "failed"
STATE_PENDING = "pending"

# SLURM states that are not terminal. Anything else means the job has stopped.
_LIVE = {"PENDING", "CONFIGURING", "RUNNING", "COMPLETING", "RESIZING",
         "SUSPENDED", "REQUEUED", "SIGNALING"}


@dataclass
class Job:
    id: str                       # ours, stable, names the log file
    sample: str
    argv: list[str]
    log: Path
    started: float
    backend: str
    backend_id: str | None = None   # SLURM job id, when there is one
    ended: float | None = None
    _state: str = STATE_RUNNING
    _rc: int | None = None
    proc: subprocess.Popen | None = field(default=None, repr=False)

    def elapsed(self) -> float:
        return round((self.ended or time.time()) - self.started, 1)


class LocalLauncher:
    """Fork the run as a child of the server. No queue, no wall clock."""

    name = "local"

    def __init__(self, log_dir: Path):
        self.log_dir = log_dir

    def submit(self, argv: list[str], sample: str, cwd: Path) -> Job:
        self.log_dir.mkdir(parents=True, exist_ok=True)
        jid = uuid.uuid4().hex[:12]
        log = self.log_dir / f"{jid}.log"
        fh = open(log, "w", buffering=1)
        fh.write(f"$ {shlex.join(argv)}\n\n")
        proc = subprocess.Popen(argv, cwd=str(cwd), stdout=fh,
                                stderr=subprocess.STDOUT, text=True)
        return Job(jid, sample, argv, log, time.time(), self.name, proc=proc)

    def refresh(self, job: Job) -> None:
        rc = job.proc.poll() if job.proc else None
        if rc is None:
            job._state = STATE_RUNNING
            return
        if job.ended is None:
            job.ended = time.time()
        job._rc = rc
        job._state = STATE_FINISHED if rc == 0 else STATE_FAILED

    def cancel(self, job: Job) -> bool:
        if job.proc and job.proc.poll() is None:
            job.proc.terminate()
            return True
        return False

    def describe(self) -> dict:
        return {"backend": self.name, "queue": None}


class SlurmLauncher:
    """Submit with sbatch; the run outlives the server and gets its own quota."""

    name = "slurm"

    def __init__(self, log_dir: Path, *, partition: str | None = None,
                 account: str | None = None, cpus: int = 4,
                 mem: str = "16G", time_limit: str = "08:00:00",
                 gres: str | None = None, extra: list[str] | None = None):
        self.log_dir = log_dir
        self.partition, self.account = partition, account
        self.cpus, self.mem, self.time_limit = cpus, mem, time_limit
        self.gres, self.extra = gres, list(extra or [])

    def _directives(self, jid: str, sample: str, log: Path) -> list[str]:
        d = ["--parsable", f"--job-name=sae-{sample}", f"--output={log}",
             f"--cpus-per-task={self.cpus}", f"--mem={self.mem}",
             f"--time={self.time_limit}"]
        # Site-specific options are only sent when set: an undefined gres or a
        # partition that does not exist is rejected at submission, not later.
        if self.partition:
            d.append(f"--partition={self.partition}")
        if self.account:
            d.append(f"--account={self.account}")
        if self.gres:
            d.append(f"--gres={self.gres}")
        return d + self.extra

    def submit(self, argv: list[str], sample: str, cwd: Path) -> Job:
        self.log_dir.mkdir(parents=True, exist_ok=True)
        jid = uuid.uuid4().hex[:12]
        log = self.log_dir / f"{jid}.log"
        script = self.log_dir / f"{jid}.sbatch"
        # A real script rather than --wrap: it survives as a record of exactly
        # what was submitted, and avoids a second layer of shell quoting.
        script.write_text(
            "#!/bin/bash\n"
            "set -euo pipefail\n"
            f"cd {shlex.quote(str(cwd))}\n"
            f"exec {shlex.join(argv)}\n"
        )
        script.chmod(0o755)
        out = subprocess.run(
            ["sbatch", *self._directives(jid, sample, log), str(script)],
            capture_output=True, text=True, cwd=str(cwd),
            env={**os.environ, "SBATCH_EXPORT": "ALL"},
        )
        if out.returncode != 0:
            raise RuntimeError(f"sbatch rejected the job: "
                               f"{(out.stderr or out.stdout).strip()}")
        # --parsable gives "jobid" or "jobid;cluster".
        backend_id = out.stdout.strip().split(";")[0]
        job = Job(jid, sample, argv, log, time.time(), self.name,
                  backend_id=backend_id)
        job._state = STATE_PENDING
        return job

    def refresh(self, job: Job) -> None:
        if job.backend_id is None:
            return
        q = subprocess.run(["squeue", "-h", "-j", job.backend_id, "-o", "%T"],
                           capture_output=True, text=True)
        state = q.stdout.strip().split("\n")[0].strip() if q.stdout.strip() else ""
        if state in _LIVE:
            job._state = STATE_PENDING if state in ("PENDING", "CONFIGURING") \
                else STATE_RUNNING
            return
        # Gone from the queue: ask accounting how it ended. sacct needs a
        # storage plugin, so treat its absence as "finished, code unknown"
        # rather than failing the whole poll.
        if job.ended is None:
            job.ended = time.time()
        a = subprocess.run(
            ["sacct", "-n", "-X", "-P", "-j", job.backend_id,
             "-o", "State,ExitCode"], capture_output=True, text=True)
        line = a.stdout.strip().split("\n")[0] if a.stdout.strip() else ""
        if not line:
            job._state = STATE_FINISHED
            return
        st, _, code = line.partition("|")
        st = st.strip().split()[0] if st.strip() else ""
        m = re.match(r"(\d+):(\d+)", code.strip())
        rc, sig = (int(m.group(1)), int(m.group(2))) if m else (None, 0)
        job._rc = -sig if sig else rc
        job._state = STATE_FINISHED if st == "COMPLETED" and not sig else STATE_FAILED

    def cancel(self, job: Job) -> bool:
        if job.backend_id and job._state in (STATE_RUNNING, STATE_PENDING):
            subprocess.run(["scancel", job.backend_id], capture_output=True)
            return True
        return False

    def describe(self) -> dict:
        return {"backend": self.name, "partition": self.partition,
                "cpus": self.cpus, "mem": self.mem, "time": self.time_limit,
                "gres": self.gres,
                "queue": shutil.which("sbatch")}


def failure_hint(job: Job) -> str | None:
    """Name a signalled death; its log is usually empty."""
    rc = job._rc
    if rc is None or rc >= 0:
        return None
    try:
        name = signal.Signals(-rc).name
    except ValueError:
        return f"killed by signal {-rc}"
    if -rc in (signal.SIGKILL, signal.SIGABRT):
        return (f"killed by {name} — usually out of memory. ESMC-6B needs "
                f"~12 GB for weights alone; try model=300m, or raise the "
                f"job's --mem.")
    return f"killed by {name}"


def make_launcher(log_dir: Path, prefer: str = "auto", **slurm_kw):
    """slurm when sbatch is on PATH, else a local fork. 'prefer' forces one."""
    if prefer not in ("auto", "slurm", "local"):
        raise ValueError("launcher must be auto, slurm or local")
    has_sbatch = shutil.which("sbatch") is not None
    if prefer == "slurm" and not has_sbatch:
        raise RuntimeError("--launcher slurm but sbatch is not on PATH; "
                           "see container/slurm-local to run one locally")
    if prefer == "local" or (prefer == "auto" and not has_sbatch):
        return LocalLauncher(log_dir)
    return SlurmLauncher(log_dir, **slurm_kw)
