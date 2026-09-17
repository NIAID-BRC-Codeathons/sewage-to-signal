"""How the web UI starts a pipeline run.

The server used to fork the pipeline as its own child. On a cluster that was
wrong twice over: the run died with the server, and it was confined to whatever
allocation the *UI* was given, so a dashboard sized for browsing could not start
real work. Runs were submitted to a scheduler instead.

The deployment is one server now, not a cluster. The allocation half of that
argument is gone with it — there is no allocation to be confined to — but the
first half stands: a run must not die because the dashboard was restarted. So
runs execute here, detached, and two things stand in for what the scheduler
used to provide.

``start_new_session`` gives each job its own session, so it outlives this
process and can be signalled as a group. And the job records its own exit code
in ``<id>.rc``, which is the local stand-in for ``sacct``: without it, a job
reaped by nobody after a restart could only be reported as "gone", which looks
exactly like a failure.

Jobs run **one at a time**, each seeing the whole machine. There is no slot
pool and nothing sets CUDA_VISIBLE_DEVICES: on a single host the useful
question is not which GPU a run gets but whether two runs are competing, and
a queue of one answers it.

There is still exactly one way a run executes. There used to be two backends
and a default that quietly chose the weaker one, so what you tested was not
what you deployed; that is the mistake this keeps avoiding, not the specific
choice of scheduler.
"""

from __future__ import annotations

import json
import os
import shlex
import signal
import subprocess
import threading
import time
import uuid
from dataclasses import dataclass
from pathlib import Path

STATE_RUNNING = "running"
STATE_FINISHED = "finished"
STATE_FAILED = "failed"
STATE_PENDING = "pending"

SCRIPT_MAX = 64 * 1024          # a submission is a few hundred bytes


@dataclass
class Job:
    id: str                       # ours, stable, names the log file
    sample: str
    argv: list[str]
    log: Path
    started: float
    backend: str
    backend_id: str | None = None   # the pid, once it has started
    script: Path | None = None      # the script, kept as the record
    ended: float | None = None
    _state: str = STATE_RUNNING
    _rc: int | None = None

    def elapsed(self) -> float:
        return round((self.ended or time.time()) - self.started, 1)


def rc_path(job: Job) -> Path:
    """Where a job records its own exit code."""
    return job.log.with_suffix(".rc")


def _read_rc(job: Job) -> tuple[int, float] | None:
    """The recorded exit code and when it was recorded, or None if not yet."""
    p = rc_path(job)
    try:
        raw = p.read_text().strip()
        when = p.stat().st_mtime
    except (OSError, ValueError):
        return None
    if not raw:
        return None
    try:
        rc = int(raw)
    except ValueError:
        return None
    # The wrapper can only report what the shell tells it, and the shell
    # reports a signalled child as 128+N. The worker overwrites this with the
    # signed value when it is still around to reap, so this branch is the
    # degraded path after a restart: 137 could in principle be a real exit
    # code, but a pipeline stage that exits 137 on purpose does not exist.
    if 128 < rc < 192:
        rc = -(rc - 128)
    return rc, when


def _pid_alive(pid: int) -> bool:
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True             # exists, just not ours to signal
    except OSError:
        return False
    return True


def parse_limit(limit: str) -> float | None:
    """Seconds from a SLURM-shaped time limit, or None if it cannot be read.

    The format is kept because the flag is unchanged and people write it from
    memory: [D-]HH:MM:SS, MM:SS, or a bare count of minutes.
    """
    if not limit:
        return None
    days = 0
    text = limit.strip()
    if "-" in text:
        d, _, text = text.partition("-")
        try:
            days = int(d)
        except ValueError:
            return None
    parts = text.split(":")
    try:
        nums = [int(p) for p in parts]
    except ValueError:
        return None
    if len(parts) == 3:
        h, m, s = nums
    elif len(parts) == 2:
        h, m, s = 0, nums[0], nums[1]
    elif len(parts) == 1:
        h, m, s = 0, nums[0], 0
    else:
        return None
    return days * 86400 + h * 3600 + m * 60 + s


class LocalLauncher:
    """Run jobs here, one at a time, detached so they outlive the server."""

    name = "local"

    def __init__(self, log_dir: Path, *, cpus: int = 4,
                 time_limit: str = "08:00:00"):
        self.log_dir = log_dir
        self.cpus = cpus
        self.time_limit = time_limit
        self._limit_s = parse_limit(time_limit)
        self._lock = threading.Lock()
        self._queue: list[Job] = []
        self._known: list[Job] = []     # everything we might have to wait on
        self._wake = threading.Event()
        self._worker: threading.Thread | None = None

    # --- what the server calls ---------------------------------------------

    def submit(self, argv: list[str], sample: str, cwd: Path) -> Job:
        self.log_dir.mkdir(parents=True, exist_ok=True)
        jid = uuid.uuid4().hex[:12]
        log = self.log_dir / f"{jid}.log"
        script = self.log_dir / f"{jid}.sh"
        # A real script rather than a bare argv: it survives as a record of
        # exactly what was run, and avoids a second layer of shell quoting.
        # The limits are in the file for the same reason the SLURM directives
        # were — split between the file and the caller, the record on disk
        # would be missing half of what decided how the job behaved.
        script.write_text(
            "#!/bin/bash\n"
            f"# sae job {jid}  sample={sample}\n"
            f"# cpus={self.cpus}  time={self.time_limit}\n"
            f"export OMP_NUM_THREADS={self.cpus}\n"
            f"cd {shlex.quote(str(cwd))}\n"
            f"{shlex.join(argv)}\n"
            # Recorded by the job itself, because after a restart there is no
            # parent left to reap it. The worker overwrites this with the
            # signed value when it is still here; this is the fallback.
            "rc=$?\n"
            f"echo \"$rc\" > {shlex.quote(str(self.log_dir / f'{jid}.rc'))}\n"
            "exit \"$rc\"\n"
        )
        script.chmod(0o755)
        job = Job(jid, sample, argv, log, time.time(), self.name, script=script)
        job._state = STATE_PENDING
        with self._lock:
            self._queue.append(job)
            self._known.append(job)
        self._ensure_worker()
        self._wake.set()
        return job

    def refresh(self, job: Job) -> None:
        if job.backend and job.backend != self.name:
            # A record from a different backend. Its id means something there,
            # not here — reading a SLURM job id as a pid would probe an
            # unrelated process, and on Linux the low pids are always-present
            # kernel threads, so every old record would report as running
            # forever. Report what is on disk and nothing more.
            self._settle_without_pid(job)
            return
        if job.backend_id is None:
            # Queued, or adopted before it ever started. Either way it has not
            # run; leaving the dataclass default would report it as running.
            if job._state != STATE_PENDING and job.ended is None:
                job._state = STATE_PENDING
            return
        rec = _read_rc(job)
        if rec is None:
            if _pid_alive(int(job.backend_id)):
                job._state = STATE_RUNNING
                return
            # Dead with nothing recorded: the job stopped but no exit code
            # reached disk. Finished, code unknown — the same call the SLURM
            # backend made when a site had no accounting storage for sacct.
            self._settle_without_pid(job)
            return
        rc, when = rec
        if job.ended is None:
            job.ended = when        # not now: a restart would inflate elapsed
        job._rc = rc
        job._state = STATE_FINISHED if rc == 0 else STATE_FAILED

    def _settle_without_pid(self, job: Job) -> None:
        """Report a job we cannot probe, from whatever reached disk."""
        rec = _read_rc(job)
        if rec is not None:
            rc, when = rec
            if job.ended is None:
                job.ended = when
            job._rc = rc
            job._state = STATE_FINISHED if rc == 0 else STATE_FAILED
            return
        if job.ended is None:
            # The last thing it wrote is the best evidence of when it stopped.
            # time.time() would count every hour since as runtime and report a
            # job from last week as having taken a week.
            try:
                job.ended = job.log.stat().st_mtime
            except OSError:
                job.ended = time.time()
        job._state = STATE_FINISHED

    def cancel(self, job: Job) -> bool:
        if job._state not in (STATE_RUNNING, STATE_PENDING):
            return False
        if job.backend_id is None:
            # Never started. Drop it from the queue and record the outcome the
            # same way a running job would, so it does not sit pending forever.
            with self._lock:
                self._queue = [j for j in self._queue if j.id != job.id]
            self._record_rc(job, -signal.SIGTERM)
            job._rc = -signal.SIGTERM
            job.ended = time.time()
            job._state = STATE_FAILED
            return True
        return self._signal(job, signal.SIGTERM)

    def describe(self) -> dict:
        with self._lock:
            queued = len(self._queue)
        return {"backend": self.name, "cpus": self.cpus,
                "time": self.time_limit, "concurrency": 1, "queued": queued}

    def adopt(self, jobs: list[Job]) -> None:
        """Take responsibility for jobs from a previous session.

        Two reasons this is not just bookkeeping. One still running must be
        waited for, or a restart would start a second run beside it and break
        the one-at-a-time promise. One that never started has no pid and no
        exit code, so without re-queueing it would sit pending forever.
        """
        requeue = []
        with self._lock:
            for j in jobs:
                self._known.append(j)
                if j.backend_id is None and _read_rc(j) is None:
                    j._state = STATE_PENDING
                    requeue.append(j)
            # Oldest first: a restart should not reorder what was waiting.
            self._queue = sorted(requeue, key=lambda j: j.started) + self._queue
        if requeue:
            self._ensure_worker()
            self._wake.set()

    # --- the worker ---------------------------------------------------------

    def _ensure_worker(self) -> None:
        with self._lock:
            if self._worker and self._worker.is_alive():
                return
            self._worker = threading.Thread(target=self._run_queue,
                                            name="sae-jobs", daemon=True)
            self._worker.start()

    def _busy(self) -> bool:
        """Is anything we know about still running?"""
        with self._lock:
            known = list(self._known)
        for j in known:
            # Same reason refresh() checks this: another backend's id is not a
            # pid. A stale SLURM record whose id happens to match a live kernel
            # thread would otherwise read as "busy" forever and the queue would
            # never start anything again.
            if j.backend != self.name or not j.backend_id:
                continue
            if _read_rc(j) is None and _pid_alive(int(j.backend_id)):
                return True
        return False

    def _run_queue(self) -> None:
        while True:
            self._wake.wait(timeout=2.0)
            self._wake.clear()
            while True:
                with self._lock:
                    if not self._queue:
                        break
                if self._busy():
                    break           # come back on the next tick
                with self._lock:
                    if not self._queue:
                        break
                    job = self._queue.pop(0)
                try:
                    self._execute(job)
                except Exception as exc:                # noqa: BLE001
                    # A worker that dies takes the queue with it, so no failure
                    # here may escape. Report it on the job instead.
                    job.ended = time.time()
                    job._state = STATE_FAILED
                    try:
                        with open(job.log, "a") as fh:
                            fh.write(f"\n[launcher] could not start: {exc}\n")
                    except OSError:
                        pass

    def _execute(self, job: Job) -> None:
        job.started = time.time()   # queued time is not run time
        with open(job.log, "ab", buffering=0) as out:
            proc = subprocess.Popen(
                ["/bin/bash", str(job.script)],
                stdout=out, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL,
                cwd=str(job.script.parent), start_new_session=True,
            )
        job.backend_id = str(proc.pid)
        job._state = STATE_RUNNING
        # Persist the pid: without it a restart cannot tell a running job from
        # one that never started.
        save_job(job, self.log_dir)
        try:
            proc.wait(timeout=self._limit_s)
        except subprocess.TimeoutExpired:
            try:
                with open(job.log, "a") as fh:
                    fh.write(f"\n[launcher] over the {self.time_limit} limit "
                             f"— terminating\n")
            except OSError:
                pass
            self._signal(job, signal.SIGTERM)
            try:
                proc.wait(timeout=30)
            except subprocess.TimeoutExpired:
                self._signal(job, signal.SIGKILL)
                proc.wait()
        # Popen gives the signed form directly, which is better than the 128+N
        # the wrapper could record, so overwrite what it wrote.
        self._record_rc(job, proc.returncode)
        job._rc = proc.returncode
        job.ended = time.time()
        job._state = STATE_FINISHED if proc.returncode == 0 else STATE_FAILED
        self._wake.set()

    def _signal(self, job: Job, sig: int) -> bool:
        if not job.backend_id:
            return False
        pid = int(job.backend_id)
        try:
            # The group, not the pid: the wrapper's child is the actual work,
            # and killing only the shell would leave the pipeline running.
            os.killpg(os.getpgid(pid), sig)
        except (ProcessLookupError, PermissionError, OSError):
            return False
        return True

    def _record_rc(self, job: Job, rc: int | None) -> None:
        if rc is None:
            return
        try:
            rc_path(job).write_text(f"{rc}\n")
        except OSError:
            pass


def save_job(job: Job, log_dir: Path) -> None:
    """Record enough to adopt this job after a restart.

    A job outlives the server, so one held only in memory is orphaned by a
    restart: still running, but invisible — which looks exactly like a failure.
    """
    try:
        (log_dir / f"{job.id}.json").write_text(json.dumps({
            "id": job.id, "sample": job.sample, "argv": job.argv,
            "log": str(job.log), "started": job.started,
            "backend": job.backend, "backend_id": job.backend_id,
            "script": str(job.script) if job.script else None,
        }))
    except OSError:
        pass                      # losing the record must not fail the submit


def load_jobs(log_dir: Path, limit: int = 200) -> list[Job]:
    """Jobs from previous sessions, newest first. State is recovered from the
    exit-code file and the pid afterwards, so what is stored here is only
    identity."""
    out: list[Job] = []
    if not log_dir.is_dir():
        return out
    for f in sorted(log_dir.glob("*.json"), key=lambda p: p.stat().st_mtime,
                    reverse=True)[:limit]:
        try:
            d = json.loads(f.read_text())
            # Records written before the path was stored still have their
            # script where submit() put it. Two extensions, because jobs
            # submitted to a scheduler were .sbatch.
            script = d.get("script")
            if not script:
                for ext in (".sh", ".sbatch"):
                    cand = log_dir / f"{d['id']}{ext}"
                    if cand.exists():
                        script = cand
                        break
                else:
                    script = log_dir / f"{d['id']}.sh"
            out.append(Job(d["id"], d["sample"], d["argv"], Path(d["log"]),
                           d["started"], d.get("backend", "local"),
                           backend_id=d.get("backend_id"),
                           script=Path(script)))
        except (json.JSONDecodeError, OSError, KeyError, TypeError):
            continue
    return out


def script_text(job: Job) -> str | None:
    """The script this job ran as, or None when it is not on disk.

    Reading it back rather than reconstructing it is the point: a rebuilt
    command shows what the server *would* run today, which is not evidence
    about a job that ran last week under different settings.
    """
    if job.script is None:
        return None
    try:
        # Bounded like the log tail, from the front: a submission is a few
        # hundred bytes, so anything past this is not one.
        with open(job.script, errors="replace") as fh:
            return fh.read(SCRIPT_MAX)
    except OSError:
        return None


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
                f"~12 GB for weights alone; try model=300m, or run it "
                f"somewhere with more memory.")
    if -rc == signal.SIGTERM:
        return f"killed by {name} — cancelled, or over the job time limit"
    return f"killed by {name}"


def make_launcher(log_dir: Path, **kw) -> LocalLauncher:
    """The one backend. There is no second path and no None case: running a
    job needs nothing that could be absent."""
    return LocalLauncher(log_dir, **kw)
