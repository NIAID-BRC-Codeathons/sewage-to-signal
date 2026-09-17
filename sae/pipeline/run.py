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
        "seq_sha1 IN (SELECT seq_sha1 FROM gene
                      WHERE sample='CHI-A' AND category='dark')"

    # the question a cohort exists to ask
    python run.py --proteins p.faa --sample S1 --where s06_embed=\\
        "seq_sha1 IN (SELECT seq_sha1 FROM gene WHERE category='dark'
                      GROUP BY seq_sha1 HAVING count(DISTINCT sample) >= 3)"

``--list`` prints the graph, ``--plan`` prints what would run and stops, and
``--next`` asks what could run against a sample as it stands - the question an
interactive driver asks instead of "what is stage 5".
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import entities
import lake
from common import MissingTool, workdir
from entities import LEVELS, is_level
from stage import Registry, Request, Stage, registry

DEFAULT_TARGET = "feature"
REPO = Path(__file__).resolve().parent.parent.parent


def _print_table(t) -> None:
    """A query result, wide enough to read and narrow enough to fit."""
    cols = t.column_names
    rows = [[("" if v is None else str(v))[:40] for v in r.values()]
            for r in t.to_pylist()]
    width = [max(len(c), *(len(r[i]) for r in rows)) if rows else len(c)
             for i, c in enumerate(cols)]
    print("  ".join(c.ljust(w) for c, w in zip(cols, width)))
    print("  ".join("-" * w for w in width))
    for r in rows:
        print("  ".join(v.ljust(w) for v, w in zip(r, width)))
    print(f"\n{len(rows)} row{'' if len(rows) == 1 else 's'}")


# --------------------------------------------------------------------------
# running one stage
# --------------------------------------------------------------------------
def resolve_where(con, reg: Registry, stage: Stage, sample: str,
                  requested: str | None, quiet: bool = False) -> str | None:
    """Pick the predicate for this stage, and say so when a default cannot apply.

    A stage's ``default_where`` describes the pipeline's usual shape, not a
    requirement: ``s06`` prefers representatives that homology could not
    explain, but on a sample where ``s05`` never ran there is nothing to mean
    by ``category``.

    That used to be detected by letting the query fail to bind, which worked
    only because each sample had its own columns. Now every sample shares one
    table, so the column exists and is merely NULL - the predicate would bind,
    match nothing, and the stage would exit claiming an empty selection. So ask
    the question directly instead: a default applies only if every stage that
    owns a column it names has actually run for this sample. The registry knows
    the ownership and ``stage_run`` knows what ran.

    A predicate somebody actually asked for is never softened - guessing on the
    user's behalf is how you silently embed the wrong 40,000 proteins.
    """
    if requested is not None:
        return entities.guard_predicate(requested) or None
    where = stage.default_where
    if not where:
        return None
    done = lake.completed_stages(con, sample)
    needed = _owners(reg, stage.consumes, where)
    absent = sorted(n for n in needed if n not in done)
    if absent:
        if not quiet:
            print(f"    note: default selection ({where}) needs "
                  f"{', '.join(absent)}, which has not run for {sample}; "
                  f"taking all rows", flush=True)
        return None
    return where


def _owners(reg: Registry, level: str, where: str) -> set[str]:
    """Stages owning any column the predicate names."""
    owner = {c.name: st.name for st in reg for c in st.columns_for(level)}
    return {owner[c] for c in entities.columns_named(where, list(owner))}


def run_stage(reg: Registry, stage: Stage, target, work: Path, sample: str,
              values: dict, where: str | None, current: dict,
              force: bool = False):
    """Invoke one stage: resolve its input, bind its params, call it.

    The lake is opened once for the whole stage and handed to it. That is a
    compromise: a stage holds it across its own compute, which for s06 is a GPU
    pass. It is the right trade anyway, because a stage that wrote a partial
    result and then failed to re-attach would leave the store inconsistent -
    one attach per stage means one transaction boundary per stage. Parallel
    batches should pin one sample per process, which is what run_batch.sh does.
    """
    out_dir = workdir(work, sample, stage.name)
    inputs: dict = {}

    with lake.open(target) as con:
        lake.ensure_schema(con, reg, LEVELS)
        if is_level(stage.consumes):
            where = resolve_where(con, reg, stage, sample, where)
            selected = entities.select(con, stage.consumes, where, sample=sample)
            if selected.num_rows == 0:
                raise SystemExit(
                    f"[{stage.name}] the selection over {stage.consumes} is empty"
                    + (f" (predicate: {where})" if where else "")
                    + ". Nothing to do.")
            inputs = {"rows": selected, "source": current.get(stage.consumes)}
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

        inputs["con"] = con
        req = Request(stage=stage, out_dir=out_dir, sample=sample, values=values,
                      where=where, force=force, inputs=inputs)
        result = stage.run(**req.kwargs())

    # Carry forward what this stage filled. Levels live in the store rather
    # than in a path, so they map to None and are found by querying it.
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
    p.add_argument("--work", type=Path, default=Path("work"),
                   help="scratch for the stages that still produce files "
                        "(reads, contigs, logs). Tabular output goes to the lake")
    p.add_argument("--lake", default=lake.default_target(REPO),
                   help="the store: a path, or postgres:/sqlite: for a shared "
                        "catalog. SAE_LAKE_DATA points the parquet elsewhere, "
                        "including s3://")
    p.add_argument("--to", dest="target", default=DEFAULT_TARGET,
                   help=f"port to produce - a level or a file kind "
                        f"(default: {DEFAULT_TARGET})")
    p.add_argument("--from", dest="start",
                   help="port you already have, if not implied by the input")
    p.add_argument("--only", action="append",
                   help="run just this stage; repeatable")
    p.add_argument("--fast", action="store_true",
                   help="route around assembly: translate reads directly into "
                        "peptides. Fast and needs no assembler, but every "
                        "peptide is a read-length fragment")
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
    p.add_argument("--sql", help="run one read-only query against the store "
                                 "and print it - the cohort, not one sample")
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
    if a.sql:
        # A whole query, not a predicate, and it asks about the cohort rather
        # than a sample - so no --sample, and no guard: this is a local CLI and
        # the attach is read-only. The web server's /api/query is the guarded
        # path, because that one takes input over a socket.
        with lake.read(a.lake) as con:
            _print_table(con.execute(a.sql).arrow().read_all())
        return
    if not a.sample:
        p.error("--sample is required")

    work = Path(a.work).resolve()

    if a.next_:
        with lake.read(a.lake) as con:
            present = entities.levels_present(con, a.sample)
            done = sorted(lake.completed_stages(con, a.sample))
        if not present:
            print(f"{a.sample}: nothing in the store yet; start from a file input.")
            return
        print(f"{a.sample} has: {', '.join(present)}")
        print(f"already run: {', '.join(done) or 'nothing'}\n")
        for s in reg.next_steps(present):
            mark = "  (again)" if s.name in done else ""
            print(f"  {s.name:<16} reads {s.consumes:<12} {s.title}{mark}")
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
        with lake.read(a.lake) as con:
            present = entities.levels_present(con, a.sample)
        have = a.start or (present or [None])[-1]
        if not have:
            p.error("give an input (--fastq/--contigs/--proteins) or --from")
    have = a.start or have

    # --fast names a capability, not a stage: anything that fills the assembly
    # role is taken out and the planner finds whatever other route exists.
    skip = set(a.skip) | ({s.name for s in reg.by_role("assembly")}
                          if a.fast else set())
    unknown = [n for n in skip if n not in reg]
    if unknown:
        sys.exit(f"--skip: no stage named {', '.join(sorted(unknown))}")

    if a.only:
        unknown = [n for n in a.only if n not in reg]
        if unknown:
            sys.exit(f"--only: no stage named {', '.join(unknown)}")
        plan = [reg[n] for n in a.only]
    else:
        # Routing happens on a registry without the skipped stages, so a skip
        # re-plans around the gap instead of leaving a hole. Dropping
        # s02_assemble from a finished plan used to strand s03_genes with no
        # contigs; now it finds the translation route on its own.
        routable = Registry([s for s in reg if s.name not in skip]) if skip else reg
        try:
            plan = routable.plan(have, a.target)
        except ValueError as exc:
            sys.exit(str(exc) + (
                f"\n(skipping {', '.join(sorted(skip))} removed the only route)"
                if skip else ""))

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
                reg, s, a.lake, work, a.sample, settings.get(s.name, {}),
                wheres.get(s.name), current, force=a.force)
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

    with lake.read(a.lake) as con:
        levels = entities.levels_present(con, a.sample)
    print(f"\nDone. {a.sample} holds: {', '.join(levels) or 'nothing'}")
    print(f"      lake: {a.lake}")


if __name__ == "__main__":
    main()
