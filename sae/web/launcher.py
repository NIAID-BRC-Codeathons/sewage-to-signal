"""How the web UI starts a pipeline run.

The server used to fork the pipeline as its own child. On a cluster that is
wrong twice over: the run dies with the server, and it is confined to whatever
allocation the *UI* was given, so a dashboard sized for browsing cannot start
real work. Runs are submitted to a scheduler instead.

There is exactly one way a run executes: ``sbatch``. There used to be a local
fork as a fallback, which meant two execution paths to keep working and a
default that quietly chose the weaker one — the deployment changed shape
depending on where it ran, so what you tested was not what you deployed.

Where there is no cluster, run one: ``container/slurm-local`` starts SLURM and
Apptainer in Docker and puts ``sbatch`` on PATH. Without a scheduler the server
still serves the dashboard — reading manifests needs nothing — and refuses to
launch, which is a missing capability rather than a second code path.
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
from dataclasses import dataclass
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

    def elapsed(self) -> float:
        return round((self.ended or time.time()) - self.started, 1)


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


NO_SCHEDULER = (
    "no scheduler: sbatch is not on PATH, and runs are only ever submitted, "
    "never forked from this server. Start one locally with "
    "`container/slurm-local/up.sh` then "
    "`eval \"$(container/slurm-local/up.sh env)\"`, or run this server on a "
    "login node."
)


def make_launcher(log_dir: Path, **slurm_kw):
    """A SlurmLauncher, or None when there is no sbatch to submit to.

    None is deliberate: the dashboard still works without a scheduler, so the
    server starts and only launching is unavailable.
    """
    if shutil.which("sbatch") is None:
        return None
    return SlurmLauncher(log_dir, **slurm_kw)
