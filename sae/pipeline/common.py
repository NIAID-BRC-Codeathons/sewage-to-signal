"""Shared plumbing for the SAE metagenomics pipeline.

Every stage follows the same contract:

* it reads declared input files and writes declared output files under
  ``work/<sample>/``;
* it writes a sidecar ``<output>.manifest.json`` recording inputs (with size +
  mtime), parameters, counts, timing and tool versions;
* it is idempotent - re-running with unchanged inputs and parameters is a
  no-op unless ``force=True``.

Stages are independently runnable, so you can enter the pipeline at whatever
point your data already reaches (reads, contigs, or proteins).

Most stages also *annotate* rather than transform: they add columns to an
entity level instead of writing a filtered copy of their input. Those columns
go in a parquet fragment and are recorded in the manifest's ``tables`` list,
which is how a reader discovers them without being told the stage list. See
``entities.py``.
"""

from __future__ import annotations

import gzip
import json
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterator

SCHEMA_VERSION = 2          # 2 added the `tables` fragment record


def open_maybe_gzip(path: Path, mode: str = "rt"):
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, mode)
    return open(path, mode)


def read_fasta(path: Path) -> Iterator[tuple[str, str]]:
    """Yield (header, sequence). Small dependency-free reader."""
    name, chunks = None, []
    with open_maybe_gzip(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name, chunks = line[1:], []
            elif name is not None:
                chunks.append(line)
    if name is not None:
        yield name, "".join(chunks)


def write_fasta(path: Path, records, width: int = 60) -> int:
    n = 0
    with open_maybe_gzip(path, "wt") as fh:
        for name, seq in records:
            fh.write(f">{name}\n")
            for i in range(0, len(seq), width):
                fh.write(seq[i : i + width] + "\n")
            n += 1
    return n


def fingerprint(path: Path) -> dict:
    p = Path(path)
    if not p.exists():
        return {"path": str(p), "exists": False}
    st = p.stat()
    return {
        "path": str(p),
        "exists": True,
        "bytes": st.st_size,
        "mtime": round(st.st_mtime, 3),
    }


def manifest_path(output: Path) -> Path:
    return Path(str(output) + ".manifest.json")


def is_current(output: Path, inputs: list[Path], params: dict) -> bool:
    """True when `output` was built from exactly these inputs and params."""
    mp = manifest_path(output)
    if not Path(output).exists() or not mp.exists():
        return False
    try:
        old = json.loads(mp.read_text())
    except (json.JSONDecodeError, OSError):
        return False
    if old.get("schema_version") != SCHEMA_VERSION:
        return False
    if old.get("params") != params:
        return False
    return old.get("inputs") == [fingerprint(p) for p in inputs]


def write_manifest(output: Path, inputs, params, stats, tools=None, seconds=None,
                   tables=None, stage=None):
    payload = {
        "schema_version": SCHEMA_VERSION,
        "stage": stage or Path(output).parent.name,
        "output": fingerprint(output),
        "inputs": [fingerprint(p) for p in inputs],
        "params": params,
        "stats": stats,
        "tools": tools or {},
        "tables": tables or [],
        "seconds": None if seconds is None else round(seconds, 2),
        "written": time.strftime("%Y-%m-%dT%H:%M:%S"),
    }
    manifest_path(output).write_text(json.dumps(payload, indent=2))
    return payload


def write_fragment(path: Path, table, level: str, role: str = "annotation",
                   where: str | None = None, help: dict | None = None) -> dict:
    """Write a stage's columns for one entity level, and describe them.

    The return value goes straight into ``write_manifest(tables=[...])``; that
    record is the only thing a reader needs in order to find these columns and
    know what they mean, which is what keeps readers from having to know the
    stage list.
    """
    import pyarrow.parquet as pq

    from entities import describe_table

    path = Path(path)
    pq.write_table(table, path, compression="zstd")
    return {"path": str(path),
            **describe_table(table, level, role, where=where, help=help)}


@dataclass
class StageResult:
    name: str
    output: Path
    stats: dict = field(default_factory=dict)
    skipped: bool = False
    seconds: float = 0.0
    mate: Path | None = None          # second mate, for paired stages
    # Ports this stage filled, for the driver to carry forward: a file kind
    # maps to a path, a level maps to None because a level lives in the work
    # directory rather than in any one file.
    produced: dict = field(default_factory=dict)

    def describe(self) -> str:
        tag = "cached" if self.skipped else f"{self.seconds:.1f}s"
        bits = " ".join(f"{k}={v}" for k, v in self.stats.items())
        return f"[{self.name}] {tag}  {bits}"


class MissingTool(RuntimeError):
    """Raised when a stage needs an external binary that is not installed."""

    def __init__(self, tool: str, hint: str):
        super().__init__(f"{tool!r} not found on PATH.\n  {hint}")
        self.tool = tool


def which(tool: str) -> str | None:
    import shutil

    return shutil.which(tool)


def workdir(root: Path, sample: str, stage: str) -> Path:
    d = Path(root) / sample / stage
    d.mkdir(parents=True, exist_ok=True)
    return d
