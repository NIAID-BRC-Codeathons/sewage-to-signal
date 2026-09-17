"""Shared plumbing for the SAE metagenomics pipeline.

Every stage follows the same contract:

* it reads its declared inputs and writes its declared columns;
* it records what it did in ``stage_run``, with inputs (size + mtime),
  parameters, counts, timing and tool versions;
* it is idempotent - re-running with unchanged inputs and parameters is a
  no-op unless ``force=True``.

Stages are independently runnable, so you can enter the pipeline at whatever
point your data already reaches (reads, contigs, or proteins).

Most stages also *annotate* rather than transform: they add columns to an
entity level instead of writing a filtered copy of their input. Those columns
go into the lake (``lake.py``), and what ran is recorded in ``stage_run`` -
which is what the ``<output>.manifest.json`` sidecars used to carry.

What is left here is the file-shaped plumbing: FASTA/FASTQ reading and writing
for the stages that still deal in files, and the work directory those files
live in.
"""

from __future__ import annotations

import gzip
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterator


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


@dataclass
class StageResult:
    name: str
    output: Path
    stats: dict = field(default_factory=dict)
    skipped: bool = False
    seconds: float = 0.0
    mate: Path | None = None          # second mate, for paired stages
    # Ports this stage filled, for the driver to carry forward: a file kind
    # maps to a path, a level maps to None because a level lives in the store
    # rather than in any one file.
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
    """Scratch for the stages that still produce files (reads, contigs, logs).

    Tabular output no longer lands here - it goes to the lake - so this is a
    working directory rather than the record of what happened.
    """
    d = Path(root) / sample / stage
    d.mkdir(parents=True, exist_ok=True)
    return d
