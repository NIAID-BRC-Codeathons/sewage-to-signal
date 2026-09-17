"""Pipeline driver.

Resolves a path through the stage graph and runs it. The driver knows how to
call *a* stage - resolve its input, bind its declared parameters, hand it the
rows a predicate selected - and nothing about which stages exist. Adding one is
dropping a module with a ``STAGE`` in it into this directory.

    python run.py --fastq reads.fastq.gz --sample S1
    python run.py --contigs contigs.fa --sample S1
    python run.py --proteins prot.faa --sample S1 --to feature_hit

Stages are idempotent, so re-running resumes rather than recomputes.

Selection replaces the old chain of filtered FASTAs. Each stage that reads an
entity level takes a SQL predicate over that level's columns, defaulting to
whatever reproduces the linear pipeline:

    # only the long dark proteins
    python run.py --proteins p.faa --sample S1 \\
        --where s06_embed="category = 'dark' AND aa_len > 200"

    # only what another sample also found dark
    python run.py --proteins p.faa --sample S1 --where s06_embed=\\
        "seq_sha1 IN (SELECT seq_sha1 FROM \\"CHI-A\\".gene WHERE category='dark')"

``--list`` prints the graph, ``--plan`` prints what would run and stops, and
``--next`` asks what could run against a sample as it stands - the question an
interactive driver asks instead of "what is stage 5".
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import entities
from common import MissingTool, workdir
from entities import LEVELS, is_level
from stage import Request, Stage, registry

DEFAULT_TARGET = "feature"


# --------------------------------------------------------------------------
# running one stage
# --------------------------------------------------------------------------
def resolve_where(stage: Stage, sample_dir: Path, requested: str | None,
                  roots: list[Path], quiet: bool = False) -> str | None:
    """Pick the predicate for this stage and check it can actually run.

    A stage's ``default_where`` describes the pipeline's usual shape, not a
    requirement - ``s06`` prefers representatives that homology could not
    explain, but if nothing has written those columns the clause names nothing
    and would be an error. So a default that does not bind is dropped with a
    note, while a predicate somebody actually asked for is always an error if
    it does not bind. Guessing on the user's behalf is how you silently embed
    the wrong 40,000 proteins.
    """
    if requested is not None:
        where = entities.guard_predicate(requested) or None
        _check(sample_dir, stage.consumes, where, roots)
        return where
    where = stage.default_where
    if not where:
        return None
    try:
        _check(sample_dir, stage.consumes, where, roots)
    except Exception as exc:
        if not quiet:
            print(f"    note: default selection ({where}) does not apply here "
                  f"- {_brief(exc)}; taking all rows", flush=True)
        return None
    return where


def _brief(exc: Exception) -> str:
    return str(exc).strip().splitlines()[0][:160]


def _check(sample_dir: Path, level: str, where: str | None, roots: list[Path]):
    con = entities.connect(sample_dir, roots)
    try:
        sql = f'SELECT 1 FROM "{level}"'
        if where:
            sql += f" WHERE ({where})"
        con.execute(sql + " LIMIT 0")
    finally:
        con.close()


def run_stage(stage: Stage, sample_dir: Path, sample: str, values: dict,
              where: str | None, current: dict, roots: list[Path],
              force: bool = False):
    """Invoke one stage: resolve its input, bind its params, call it."""
    out_dir = workdir(sample_dir.parent, sample, stage.name)
    inputs: dict = {}
    selected = None

    if is_level(stage.consumes):
        where = resolve_where(stage, sample_dir, where, roots)
        selected = entities.select(sample_dir, stage.consumes, where, roots=roots)
        if selected.num_rows == 0:
            raise SystemExit(
                f"[{stage.name}] the selection over {stage.consumes} is empty"
                + (f" (predicate: {where})" if where else "")
                + ". Nothing to do.")
        base = next((f for f in entities.fragments(sample_dir, stage.consumes)
                     if f.role == "base"), None)
        inputs = {"rows": selected, "source": base.path if base else None}
        print(f"    selecting {selected.num_rows} {stage.consumes} rows"
              + (f" where {where}" if where else " (all)"), flush=True)
    else:
        path = current.get(stage.consumes)
        if path is None:
            raise SystemExit(f"[{stage.name}] nothing produced {stage.consumes!r}")
        inputs = {stage.input_arg or stage.consumes: path}
        mate = current.get(stage.consumes + "2")
        if mate and stage.mate_arg:
            inputs[stage.mate_arg] = mate

    req = Request(stage=stage, out_dir=out_dir, sample=sample, values=values,
                  where=where, force=force, inputs=inputs)
    result = stage.run(**req.kwargs())

    # Carry forward what this stage filled. Levels live in the work directory
    # rather than in a path, so they map to None and are found by reading it.
    produced = {stage.produces: result.output if not is_level(stage.produces) else None}
    produced.update(result.produced or {})
    if result.mate:
        produced[stage.produces + "2"] = result.mate
    return result, produced


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------
def parse_assignments(items: list[str] | None, reg) -> dict[str, str]:
    """``--where stage=expr`` / ``--set stage.param=value`` into a dict."""
    out: dict[str, str] = {}
    for item in items or []:
        key, sep, value = item.partition("=")
        if not sep:
            sys.exit(f"expected name=value, got {item!r}")
        out[key.strip()] = value
    return out


def describe_graph(reg) -> str:
    lines = ["Stages (a stage is a module here exporting STAGE):", ""]
    for s in reg:
        tools = ", ".join(
            f"{t.name}{'' if t.available() else ' (missing)'}" for t in s.requires)
        lines.append(f"  {s.name:<16} {s.consumes:>12} -> {s.produces:<12} "
                     f"{s.kind:<9} {s.title}")
        if s.summary:
            lines.append(f"  {'':<16} {s.summary}")
        if s.default_where:
            lines.append(f"  {'':<16} selects: {s.default_where}")
        if tools:
            lines.append(f"  {'':<16} tools: {tools}")
        lines.append("")
    lines.append("Levels:")
    for lv in LEVELS.values():
        lines.append(f"  {lv.name:<14} key ({', '.join(lv.key)})  {lv.summary}")
    return "\n".join(lines)


def main():
    reg = registry()
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    src = p.add_mutually_exclusive_group()
    src.add_argument("--fastq", type=Path)
    src.add_argument("--contigs", type=Path)
    src.add_argument("--proteins", type=Path)
    p.add_argument("--fastq2", type=Path, help="second mate (with --fastq)")
    p.add_argument("--sample")
    p.add_argument("--work", type=Path, default=Path("work"))
    p.add_argument("--to", dest="target", default=DEFAULT_TARGET,
                   help=f"port to produce - a level or a file kind "
                        f"(default: {DEFAULT_TARGET})")
    p.add_argument("--from", dest="start",
                   help="port you already have, if not implied by the input")
    p.add_argument("--only", action="append",
                   help="run just this stage; repeatable")
    p.add_argument("--skip", action="append", default=[],
                   help="drop this stage from the plan; repeatable")
    p.add_argument("--where", action="append", metavar="STAGE=EXPR",
                   help="SQL predicate selecting what a stage reads")
    p.add_argument("--set", action="append", metavar="STAGE.PARAM=VALUE",
                   dest="settings", help="override a stage parameter")
    p.add_argument("--force", action="store_true")
    p.add_argument("--plan", action="store_true", help="print the plan and stop")
    p.add_argument("--list", action="store_true", help="print the stage graph")
    p.add_argument("--next", action="store_true", dest="next_",
                   help="print what could run against this sample right now")
    # Kept because they are the flags people have in their notes; each one is
    # just a --set against the stage that declares it.
    for legacy, target in (("--max-reads", "s01_qc.max_reads"),
                           ("--min-aa", "s03_genes.min_aa"),
                           ("--hmm", "s05_prefilter.hmm"),
                           ("--ref", "s05_prefilter.ref"),
                           ("--evalue", "s05_prefilter.evalue"),
                           ("--confident-evalue", "s05_prefilter.confident_evalue"),
                           ("--min-coverage", "s05_prefilter.min_coverage"),
                           ("--bit-cutoffs", "s05_prefilter.bit_cutoffs"),
                           ("--model", "s06_embed.model"),
                           ("--backbone", "s06_embed.backbone"),
                           ("--sae-repo", "s06_embed.sae_repo"),
                           ("--layer", "s06_embed.layer"),
                           ("--top-k", "s06_embed.top_k"),
                           ("--batch-size", "s06_embed.batch_size"),
                           ("--max-len", "s06_embed.max_len"),
                           ("--limit", "s06_embed.limit"),
                           ("--device", "s06_embed.device"),
                           ("--reps", "s07_match.reps")):
        p.add_argument(legacy, dest=f"legacy::{target}", default=None)
    a = p.parse_args()

    if a.list:
        print(describe_graph(reg))
        return
    if not a.sample:
        p.error("--sample is required")

    sample_dir = Path(a.work).resolve() / a.sample
    roots = [Path(a.work).resolve()]

    if a.next_:
        present = entities.levels_present(sample_dir)
        if not present:
            print(f"{a.sample}: no entity levels yet; start from a file input.")
            return
        print(f"{a.sample} has: {', '.join(present)}\n")
        for s in reg.next_steps(present):
            print(f"  {s.name:<16} reads {s.consumes:<12} {s.title}")
        return

    # -- parameter overrides: --set wins over a legacy flag for the same param
    settings: dict[str, dict] = {}
    for key, value in vars(a).items():
        if key.startswith("legacy::") and value is not None:
            stage_name, _, param = key[len("legacy::"):].partition(".")
            settings.setdefault(stage_name, {})[param] = value
    for key, value in parse_assignments(a.settings, reg).items():
        stage_name, _, param = key.rpartition(".")
        if stage_name not in reg:
            sys.exit(f"--set: no stage named {stage_name!r}")
        if reg[stage_name].param(param) is None:
            sys.exit(f"--set: {stage_name} has no parameter {param!r}; have "
                     + ", ".join(pp.name for pp in reg[stage_name].params))
        settings.setdefault(stage_name, {})[param] = value

    wheres = parse_assignments(a.where, reg)
    for name in wheres:
        if name not in reg:
            sys.exit(f"--where: no stage named {name!r}")

    # -- what we have, and the plan to get where we are going
    current: dict[str, Path] = {}
    if a.fastq:
        current["reads"] = a.fastq.resolve()
        if a.fastq2:
            current["reads2"] = a.fastq2.resolve()
        have = "reads"
    elif a.contigs:
        current["contigs"] = a.contigs.resolve()
        have = "contigs"
    elif a.proteins:
        current["proteins"] = a.proteins.resolve()
        have = "proteins"
    else:
        have = a.start or (entities.levels_present(sample_dir) or [None])[-1]
        if not have:
            p.error("give an input (--fastq/--contigs/--proteins) or --from")
    have = a.start or have

    if a.only:
        unknown = [n for n in a.only if n not in reg]
        if unknown:
            sys.exit(f"--only: no stage named {', '.join(unknown)}")
        plan = [reg[n] for n in a.only]
    else:
        try:
            plan = reg.plan(have, a.target)
        except ValueError as exc:
            sys.exit(str(exc))
        plan = [s for s in plan if s.name not in a.skip]

    if a.plan:
        print(f"{have} -> {a.target}:")
        for s in plan:
            sel = wheres.get(s.name, s.default_where)
            print(f"  {s.name:<16} {s.consumes:>12} -> {s.produces:<12}"
                  + (f"  where {sel}" if sel and is_level(s.consumes) else ""))
        return

    results = []
    for s in plan:
        try:
            r, produced = run_stage(
                s, sample_dir, a.sample, settings.get(s.name, {}),
                wheres.get(s.name), current, roots, force=a.force)
        except MissingTool as exc:
            print(f"\n[{s.name}] BLOCKED: {exc}\n", file=sys.stderr)
            print("Completed stages:", file=sys.stderr)
            for done in results:
                print("  " + done.describe(), file=sys.stderr)
            sys.exit(2)
        print(r.describe(), flush=True)
        results.append(r)
        current.update({k: v for k, v in produced.items() if v is not None})
        for port in produced:
            current.setdefault(port, None)

    levels = entities.levels_present(sample_dir)
    print(f"\nDone. {a.sample} has levels: {', '.join(levels) or 'none'}")


if __name__ == "__main__":
    main()
