"""Stage descriptors, the registry, and planning.

A stage used to be a name in a list in ``run.py`` plus a branch in an if/elif
chain that knew its keyword arguments, and the web UI mirrored that list and
those arguments by hand. Adding a stage meant editing four files, and the UI
could only ever show the one pipeline somebody had typed out.

Here a stage describes itself:

    STAGE = Stage(
        name="s04_derep",
        title="Dereplicate",
        consumes="gene", produces="gene",
        adds=(Column("rep_id", "string", "cluster representative"), ...),
        params=(Param("identity", float, 0.95, help="clustering identity"),),
        requires=(Tool("mmseqs", optional=True, hint="..."),),
        run=run,
    )

and everything else is derived. ``consumes``/``produces`` name an entity level
or a file kind, which makes the pipeline a graph rather than a line: ``plan()``
answers "I have contigs and I want features" without anybody writing that
sequence down, which is also what lets a caller ask "what could run next on
this sample?" - the question an interactive or adaptive driver has to answer.

Discovery is a directory scan for modules exporting ``STAGE``, so installing a
stage is dropping in a file.
"""

from __future__ import annotations

import importlib
import sys
from collections import deque
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable

from entities import FILE_KINDS, LEVELS, is_level

HERE = Path(__file__).resolve().parent

# Modules that are infrastructure rather than stages.
_NOT_STAGES = {"common", "entities", "lake", "stage", "run", "backfill",
               "import_projection"}


@dataclass(frozen=True)
class Param:
    """One tunable, described well enough for a UI to render it unaided."""

    name: str
    type: type = str
    default: Any = None
    choices: tuple[str, ...] | None = None
    help: str = ""
    group: str = "general"
    path: bool = False             # value names a file; callers must check it
    # Extensions a path parameter accepts. A UI offering a file picker reads
    # this instead of knowing that "the hmm field wants .hmm files".
    suffixes: tuple[str, ...] = ()
    # Point this at the first matching file when there is one. True only where
    # leaving it empty is the worse default - s05 without a profile database
    # is a pass-through that sends every protein to the GPU. Alternatives are
    # left empty: picking one of several mutually exclusive backends on the
    # strength of a file extension is how you silently run the wrong search.
    prefer_available: bool = False
    cli: str | None = None         # flag name, when it is not --<name>

    @property
    def flag(self) -> str:
        return self.cli or "--" + self.name.replace("_", "-")

    def to_json(self) -> dict:
        return {
            "name": self.name, "type": self.type.__name__,
            "default": self.default,
            "choices": list(self.choices) if self.choices else None,
            "help": self.help, "group": self.group, "path": self.path,
            "suffixes": list(self.suffixes),
            "prefer_available": self.prefer_available, "flag": self.flag,
        }

    def coerce(self, value):
        if value in (None, ""):
            return None
        if self.type is bool:
            return bool(value) if not isinstance(value, str) else \
                value.lower() not in ("", "0", "false", "no")
        try:
            v = self.type(value)
        except (TypeError, ValueError):
            raise ValueError(f"{self.name} must be {self.type.__name__}")
        if self.choices and str(v) not in self.choices:
            raise ValueError(f"{self.name} must be one of {', '.join(self.choices)}")
        return v


@dataclass(frozen=True)
class Column:
    """A column a stage adds to its output level."""

    name: str
    type: str = "string"
    help: str = ""
    # Names a thing rather than describing it. Colouring a plot by an
    # identifier gives every point its own colour and filtering by one selects
    # a single row, so a UI offers neither. Declared rather than measured: an
    # id is unique *within a sample*, so across a cohort it looks merely
    # high-cardinality and a count-based rule lets it through.
    identifier: bool = False

    def to_json(self) -> dict:
        return {"name": self.name, "type": self.type, "help": self.help,
                "identifier": self.identifier}


@dataclass(frozen=True)
class Tool:
    """An external binary. Optional ones change behaviour; required ones block."""

    name: str
    optional: bool = True
    hint: str = ""

    def available(self) -> bool:
        import shutil
        return shutil.which(self.name) is not None

    def to_json(self) -> dict:
        return {"name": self.name, "optional": self.optional,
                "hint": self.hint, "available": self.available()}


@dataclass(frozen=True)
class Stage:
    name: str
    run: Callable
    consumes: str                      # level name or file kind
    produces: str                      # level name or file kind - the main one
    # Levels a stage fills as a side effect. s06 emits feature hits, but it
    # also defines the feature level and tells genes they were embedded, and a
    # planner that only knew the main output could not route through it.
    also_produces: tuple[str, ...] = ()
    title: str = ""
    summary: str = ""
    # Columns this stage writes onto `produces`. Load-bearing: the store's DDL
    # is generated from these, so a column that is written but not declared will
    # not exist.
    adds: tuple[Column, ...] = ()
    # Columns written onto a level in `also_produces`, as (level, columns).
    # A stage that fills several levels writes different columns to each, so
    # `adds` alone cannot describe it - s06 writes activations to feature_hit
    # but tells genes they were embedded.
    also_adds: tuple[tuple[str, tuple[Column, ...]], ...] = ()
    params: tuple[Param, ...] = ()
    requires: tuple[Tool, ...] = ()
    roles: tuple[str, ...] = ()        # what this stage's output is good for
    selectable: bool = False           # accepts a predicate over `consumes`
    default_where: str | None = None   # selection reproducing the linear default
    optional: bool = False             # a plan may skip it when nothing asks for it
    order: int = 0                     # tie-break when several paths are equal
    # A stage consuming a level is always called run(rows, out_dir, sample,
    # source=, where=, ...). A stage consuming a *file* names the keyword its
    # path arrives on, plus a second one when it can take paired mates.
    input_arg: str = ""
    mate_arg: str | None = "fastq2"

    # -- derived
    @property
    def kind(self) -> str:
        if not is_level(self.produces):
            return "file"
        return "annotate" if self.consumes == self.produces else "emit"

    @property
    def outputs(self) -> tuple[str, ...]:
        return (self.produces, *self.also_produces)

    @property
    def level(self) -> str | None:
        """The level this stage writes columns to, if any."""
        return self.produces if is_level(self.produces) else None

    def columns_for(self, level: str) -> tuple[Column, ...]:
        """The columns this stage writes onto one level, or () if it writes none."""
        if level == self.produces:
            return self.adds
        for name, cols in self.also_adds:
            if name == level:
                return cols
        return ()

    def param(self, name: str) -> Param | None:
        return next((p for p in self.params if p.name == name), None)

    def missing_tools(self) -> list[Tool]:
        return [t for t in self.requires if not t.optional and not t.available()]

    def to_json(self) -> dict:
        return {
            "name": self.name, "title": self.title or self.name,
            "summary": self.summary, "kind": self.kind,
            "consumes": self.consumes, "produces": self.produces,
            "also_produces": list(self.also_produces),
            "consumes_level": is_level(self.consumes),
            "produces_level": is_level(self.produces),
            "adds": [c.to_json() for c in self.adds],
            "also_adds": {lvl: [c.to_json() for c in cols]
                          for lvl, cols in self.also_adds},
            "params": [p.to_json() for p in self.params],
            "requires": [t.to_json() for t in self.requires],
            "roles": list(self.roles), "selectable": self.selectable,
            "default_where": self.default_where, "optional": self.optional,
        }


# --------------------------------------------------------------------------
# registry
# --------------------------------------------------------------------------
class Registry:
    """Every stage this installation can run, plus the graph they form."""

    def __init__(self, stages: list[Stage]):
        self.stages: dict[str, Stage] = {}
        for s in sorted(stages, key=lambda s: (s.order, s.name)):
            if s.name in self.stages:
                raise ValueError(f"two stages named {s.name!r}")
            self.stages[s.name] = s

    def __iter__(self):
        return iter(self.stages.values())

    def __contains__(self, name):
        return name in self.stages

    def __len__(self):
        return len(self.stages)

    def __getitem__(self, name: str) -> Stage:
        try:
            return self.stages[name]
        except KeyError:
            raise KeyError(f"no stage named {name!r}; have "
                           f"{', '.join(self.stages)}") from None

    @property
    def names(self) -> list[str]:
        return list(self.stages)

    def by_role(self, role: str) -> list[Stage]:
        """Stages whose output fills a role, e.g. 'projection' or 'triage'.

        Callers ask for a capability instead of naming a stage, so a pipeline
        that fills the role differently still works.
        """
        return [s for s in self if role in s.roles]

    def producing(self, port: str) -> list[Stage]:
        return [s for s in self if port in s.outputs]

    def consuming(self, port: str) -> list[Stage]:
        return [s for s in self if s.consumes == port]

    # -- planning
    def plan(self, have: str, want: str) -> list[Stage]:
        """Shortest stage path from a port you have to one you want.

        Both ends are ports - a level name or a file kind - so "I have contigs,
        I want feature_hit" is answerable without a hardcoded order. Annotators
        consume and produce the same level, so a straight BFS over edges would
        never schedule them; they are inserted afterwards, in declared order,
        wherever their level is live.
        """
        # A port that only ever appears as an *annotator's* side output - s07's
        # cluster hits - is never reached by walking main outputs, because the
        # walk skips stages that return to the level they read. Aim at the
        # level that annotator consumes instead and let the weave pick it up.
        if not any(s.produces == want for s in self.producing(want)):
            side = self.producing(want)
            if side:
                want = side[0].consumes

        if have == want:
            return self._weave(want, [])
        seen = {have}
        q: deque[tuple[str, list[Stage]]] = deque([(have, [])])
        best: list[Stage] | None = None
        while q:
            port, path = q.popleft()
            for s in sorted(self.consuming(port), key=lambda s: (s.order, s.name)):
                if s.produces == port:            # annotator; added below
                    continue
                nxt = path + [s]
                if want in s.outputs:
                    best = nxt
                    q.clear()
                    break
                for out in s.outputs:
                    if out not in seen:
                        seen.add(out)
                        q.append((out, nxt))
        if best is None:
            raise ValueError(
                f"no path from {have!r} to {want!r}. Stages produce: "
                + ", ".join(sorted({s.produces for s in self})))
        return self._weave(want, best)

    def _weave(self, want: str, path: list[Stage]) -> list[Stage]:
        """Insert same-port stages into a path of port-to-port ones.

        The BFS only walks edges that change port, so it never schedules a
        stage that returns to the port it read - every annotator, and ``s01``,
        which rewrites reads as reads. Those belong *before* whatever consumes
        that port, because their whole purpose is to have already run when it
        does: ``s06`` selecting on ``category`` is only meaningful once ``s05``
        has written it.
        """
        out: list[Stage] = []
        added: set[str] = set()

        def same_port(port: str) -> list[Stage]:
            return [a for a in sorted(self.consuming(port),
                                      key=lambda a: (a.order, a.name))
                    if a.produces == port and a.name not in added]

        for s in path:
            for a in same_port(s.consumes):
                out.append(a)
                added.add(a.name)
            out.append(s)
            added.add(s.name)
        for a in same_port(want):
            out.append(a)
            added.add(a.name)
        return out

    def next_steps(self, present: list[str]) -> list[Stage]:
        """Stages runnable right now, given the levels/kinds a sample has.

        This is the question an interactive driver asks: not "what is stage 5"
        but "what could I do to this sample". Ordering stays out of it.
        """
        return [s for s in self if s.consumes in present]

    def to_json(self) -> dict:
        return {
            "stages": [s.to_json() for s in self],
            "levels": [{"name": lv.name, "key": list(lv.key), "parent": lv.parent,
                        "title": lv.title, "summary": lv.summary}
                       for lv in LEVELS.values()],
            "file_kinds": [{"name": k, "summary": v} for k, v in FILE_KINDS.items()],
        }


def discover(package_dir: Path = HERE) -> Registry:
    """Import every module here that exports ``STAGE``.

    A stage is installed by being present. Nothing lists them, so nothing has
    to be edited to add one - which is the whole point, since the web UI builds
    itself out of whatever this returns.
    """
    if str(package_dir) not in sys.path:
        sys.path.insert(0, str(package_dir))
    found: list[Stage] = []
    for path in sorted(package_dir.glob("*.py")):
        if path.stem.startswith("_") or path.stem in _NOT_STAGES:
            continue
        try:
            mod = importlib.import_module(path.stem)
        except Exception as exc:                  # a broken stage is not fatal
            print(f"  warning: cannot load stage module {path.name}: {exc}",
                  file=sys.stderr)
            continue
        st = getattr(mod, "STAGE", None)
        if isinstance(st, Stage):
            found.append(st)
    return Registry(found)


_REGISTRY: Registry | None = None


def registry() -> Registry:
    global _REGISTRY
    if _REGISTRY is None:
        _REGISTRY = discover()
    return _REGISTRY


# --------------------------------------------------------------------------
# calling a stage
# --------------------------------------------------------------------------
@dataclass
class Request:
    """Everything a driver needs to invoke one stage once."""

    stage: Stage
    out_dir: Path
    sample: str
    values: dict[str, Any] = field(default_factory=dict)
    where: str | None = None
    force: bool = False
    inputs: dict[str, Any] = field(default_factory=dict)

    def kwargs(self) -> dict:
        """Declared params, coerced, plus whatever the driver resolved."""
        out = {"out_dir": self.out_dir, "sample": self.sample, **self.inputs}
        for p in self.stage.params:
            if p.name in self.values:
                v = p.coerce(self.values[p.name])
            else:
                v = p.default
            if v is not None or p.name in self.values:
                out[p.name] = v
        if self.stage.selectable:
            out["where"] = self.where
        out["force"] = self.force
        return out
