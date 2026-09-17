"""Stage 01 - read QC.

Sliding-window quality trimming, length and N filtering, with optional
subsampling. Implemented with pyfastx so it needs no external binary; if fastp
is on PATH it is preferred, being far faster and adaptor-aware.
"""

from __future__ import annotations

import argparse
import subprocess
import time
from pathlib import Path

from stage import Param, Stage, Tool
import lake
from common import (
    StageResult,
    open_maybe_gzip,
    which,
    workdir,
)


def _trim(seq: str, qual: str, min_q: int, window: int) -> tuple[str, str]:
    """Trim from the 3' end until a `window`-wide mean quality clears `min_q`."""
    scores = [ord(c) - 33 for c in qual]
    end = len(scores)
    while end >= window:
        if sum(scores[end - window : end]) / window >= min_q:
            break
        end -= 1
    start = 0
    while start < end and scores[start] < min_q:
        start += 1
    return seq[start:end], qual[start:end]


def _passes(seq, qual, min_q, window, min_len, max_n_frac):
    s, q = _trim(seq, qual, min_q, window)
    if len(s) < min_len or (s.count("N") / max(len(s), 1)) > max_n_frac:
        return None
    return s, q


def run(
    fastq: Path,
    out_dir: Path,
    sample: str,
    con=None,
    fastq2: Path | None = None,
    min_q: int = 20,
    window: int = 4,
    min_len: int = 50,
    max_n_frac: float = 0.1,
    max_reads: int | None = None,
    use_fastp: bool = True,
    force: bool = False,
) -> StageResult:
    """QC one FASTQ, or a mate pair when `fastq2` is given.

    Paired mode keeps a pair only when *both* mates pass, so the two output
    files stay index-synchronised - assemblers reject mismatched mates.
    """
    fastq = Path(fastq)
    paired = fastq2 is not None
    fastq2 = Path(fastq2) if paired else None
    engine = "fastp" if (use_fastp and which("fastp")) else "pyfastx"
    out = Path(out_dir) / (f"{sample}.qc_1.fastq.gz" if paired else f"{sample}.qc.fastq.gz")
    out2 = Path(out_dir) / f"{sample}.qc_2.fastq.gz" if paired else None
    params = {
        "min_q": min_q, "window": window, "min_len": min_len,
        "max_n_frac": max_n_frac, "max_reads": max_reads,
        "engine": engine, "paired": paired,
    }
    deps = [lake.fingerprint(f) for f in [fastq] + ([fastq2] if paired else [])]
    if not force and lake.is_current(con, sample, "s01_qc", params, deps):
        r = StageResult("s01_qc", out, {"engine": engine, "paired": paired}, skipped=True)
        r.mate = out2
        return r

    t0 = time.time()
    if engine == "fastp":
        json_report = Path(out_dir) / f"{sample}.fastp.json"
        cmd = ["fastp", "-i", str(fastq), "-o", str(out)]
        if paired:
            cmd += ["-I", str(fastq2), "-O", str(out2)]
        cmd += [
            "--cut_tail", "--cut_tail_mean_quality", str(min_q),
            "--cut_window_size", str(window), "--length_required", str(min_len),
            "--json", str(json_report), "--html", "/dev/null", "--thread", "4",
        ]
        if max_reads:
            cmd += ["--reads_to_process", str(max_reads)]
        subprocess.run(cmd, check=True, capture_output=True)
        stats = {"engine": "fastp", "paired": paired, "report": json_report.name}
    else:
        import pyfastx

        kept = dropped = seen = 0
        if paired:
            it1 = pyfastx.Fastq(str(fastq), build_index=False)
            it2 = pyfastx.Fastq(str(fastq2), build_index=False)
            with open_maybe_gzip(out, "wt") as f1, open_maybe_gzip(out2, "wt") as f2:
                for (n1, s1, q1), (n2, s2, q2) in zip(it1, it2):
                    seen += 1
                    if max_reads and seen > max_reads:
                        break
                    a = _passes(s1, q1, min_q, window, min_len, max_n_frac)
                    b = _passes(s2, q2, min_q, window, min_len, max_n_frac)
                    if a is None or b is None:
                        dropped += 1          # drop the pair, not one mate
                        continue
                    f1.write(f"@{n1}\n{a[0]}\n+\n{a[1]}\n")
                    f2.write(f"@{n2}\n{b[0]}\n+\n{b[1]}\n")
                    kept += 1
            unit = "pairs"
        else:
            with open_maybe_gzip(out, "wt") as fh:
                # build_index=False streams instead of writing a .fxi sidecar,
                # which matters for multi-GB inputs.
                for name, seq, qual in pyfastx.Fastq(str(fastq), build_index=False):
                    seen += 1
                    if max_reads and seen > max_reads:
                        break
                    a = _passes(seq, qual, min_q, window, min_len, max_n_frac)
                    if a is None:
                        dropped += 1
                        continue
                    fh.write(f"@{name}\n{a[0]}\n+\n{a[1]}\n")
                    kept += 1
            unit = "reads"
        stats = {
            "engine": "pyfastx", "paired": paired,
            f"{unit}_in": seen, f"{unit}_kept": kept, f"{unit}_dropped": dropped,
            "kept_frac": round(kept / seen, 4) if seen else 0.0,
        }

    el = time.time() - t0
    lake.record_run(con, sample, "s01_qc", params, deps, None, stats,
                    seconds=el)
    r = StageResult("s01_qc", out, stats, seconds=el)
    r.mate = out2
    return r


STAGE = Stage(
    name="s01_qc",
    title="Read QC",
    summary="Sliding-window quality trimming, length and N filtering, with "
            "optional subsampling. Reads stay files: a per-read table is the "
            "one place the annotation model does not pay.",
    run=run,
    consumes="reads",
    produces="reads",
    order=10,
    input_arg="fastq",
    params=(
        Param("min_q", int, 20, group="qc",
              help="mean quality a window must reach"),
        Param("window", int, 4, group="qc", help="sliding window size"),
        Param("min_len", int, 50, group="qc",
              help="drop reads shorter than this after trimming"),
        Param("max_n_frac", float, 0.1, group="qc",
              help="drop reads with more than this fraction of Ns"),
        Param("max_reads", int, None, group="qc",
              help="subsample to this many reads"),
    ),
    requires=(Tool("fastp", optional=True,
                   hint="far faster and adaptor-aware; pyfastx is used without it"),),
)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("fastq", type=Path)
    p.add_argument("--fastq2", type=Path, help="second mate for paired input")
    p.add_argument("--sample", required=True)
    p.add_argument("--work", type=Path, default=Path("work"))
    p.add_argument("--min-q", type=int, default=20)
    p.add_argument("--min-len", type=int, default=50)
    p.add_argument("--max-reads", type=int)
    p.add_argument("--no-fastp", action="store_true")
    p.add_argument("--force", action="store_true")
    a = p.parse_args()
    r = run(
        a.fastq, workdir(a.work, a.sample, "s01_qc"), a.sample, fastq2=a.fastq2,
        min_q=a.min_q, min_len=a.min_len, max_reads=a.max_reads,
        use_fastp=not a.no_fastp, force=a.force,
    )
    print(r.describe())


if __name__ == "__main__":
    main()
