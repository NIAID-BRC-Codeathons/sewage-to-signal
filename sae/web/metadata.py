"""Gene-level metadata for the map: what can be coloured by, and filtered on.

The store holds one ``gene`` table, so a gene already carries every column any
stage has written about it - coordinates and length from s03, cluster membership from s04, the homology
triage from s05, embedding status from s06. That is the metadata. Nothing here
adds a column; it reads what the pipeline already recorded and describes it
well enough for a UI to build its own controls.

Which is the point: **the frontend never learns a stage's schema.** A stage
added tomorrow that annotates genes with a taxon call shows up as another
column to colour by, with no change to the server or the page - the same
property the artifact browser already has.

Two things this owes the caller.

* **A stable colour domain.** The value -> slot mapping is computed over the
  *unfiltered* column, so filtering the plot never repaints the survivors. A
  colour that changes meaning when you narrow the selection is worse than no
  colour.
* **A safe filter.** The browser sends structured terms - column, operator,
  values - and the SQL is composed here against the column list the sample
  actually has. No predicate string crosses the wire, so there is nothing to
  escape at the boundary; literals are coerced by the column's own type and
  the result is still put through ``entities.guard_predicate``.

Old work directories predate the fragment manifests and have no gene level at
all. They fall back to whatever s05 wrote beside its output, which yields the
one ``category`` column and nothing else - enough to keep the map coloured as
it always was, rather than going blank on a run nobody has re-run yet.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path
from typing import Any, Sequence

# Colours are assigned in fixed slot order and a scatter validates only three
# of them all-pairs, so a categorical column shows its three commonest values
# and folds the rest into one neutral "other". Raising this silently produces
# a plot whose colours cannot be told apart - see the palette note in
# index.html.
COLOR_SLOTS = 3
# Distinct values offered in a filter dropdown. Far more than can be coloured,
# because picking one value out of 195 families is a perfectly good filter.
MAX_FILTER_VALUES = 300
# Ramp steps. Five, because an eight-step blue ramp fails the ordinal
# gate in both themes - adjacent lightness gaps come out at 0.047 against
# a 0.06 floor, and on the light surface the lightest step drops to
# 1.85:1 and sinks into the background. Five clears both.
NUMERIC_BINS = 5
# A column with one value per row is an identifier, not metadata: colouring by
# it gives every point its own colour and filtering by it selects one protein.
# This is the backstop. The primary signal is the stage's own `identifier`
# declaration, because a gene id is unique only *within* a sample - pooled
# across a cohort it looks merely high-cardinality and slips under any ratio.
IDENTITY_FRAC = 0.98

# Columns that are never useful to colour or filter by, whatever they contain.
HIDDEN = {"seq", "sample"}

# Where a level's values have a meaning-order, it wins over frequency. The
# triage classes are the case that matters: known/partial/dark is a documented
# gradient and the UI has shown them in fixed colours since the map existed, so
# ordering them by count would repaint established meaning.
CANONICAL_ORDER: dict[str, tuple[str, ...]] = {
    "category": ("known", "partial", "dark"),
}

_NUMERIC_PREFIXES = ("int", "float", "double", "decimal")


def _is_numeric(arrow_type: str) -> bool:
    return str(arrow_type).lower().startswith(_NUMERIC_PREFIXES)


# --------------------------------------------------------------------------
# loading
# --------------------------------------------------------------------------
def _entities():
    """`entities` from the pipeline, or None where duckdb is absent.

    The web server is meant to run in places the pipeline does not, so a
    missing dependency degrades to the old behaviour instead of a 500.
    """
    try:
        import entities
        import duckdb                                   # noqa: F401
        return entities
    except Exception:
        return None


def load(con, sample: str) -> tuple[list[dict], dict[str, dict]]:
    """(column descriptors, gene_id -> row) for one sample.

    One query for the whole level rather than a round trip per column: the
    descriptors need counts over every value anyway, and a gene table is small
    next to the feature hits it explains.

    The caller owns the connection, because the store is held briefly and
    deliberately - see ``lake.py``.
    """
    E = _entities()
    if E is None:
        return [], {}
    try:
        cols = [c for c in E.columns_of(con, "gene")
                if c["name"] not in HIDDEN]
        t = E.select(con, "gene", columns=[c["name"] for c in cols],
                     sample=sample)
    except Exception:
        return [], {}
    if t.num_rows == 0:
        return [], {}
    help_by = {c["name"]: c for c in cols}
    rows = {r["gene_id"]: r for r in t.to_pylist()}
    return _describe(t, help_by), rows


# --------------------------------------------------------------------------
# describing
# --------------------------------------------------------------------------
# A gene is unique per (sample, gene_id), so anything holding the whole cohort
# in one dict needs both. Unit separator: it cannot occur in a sample name
# (SAFE_NAME) or in a FASTA header id.
KEY_SEP = "\x1f"


def cohort_key(sample: str, gene_id: str) -> str:
    return f"{sample}{KEY_SEP}{gene_id}"


def load_cohort(con) -> tuple[list[dict], dict[str, dict]]:
    """(column descriptors, (sample, gene_id) -> row) for every sample.

    The per-sample `load` describes one run; this describes the store. The
    descriptors are computed over the whole cohort on purpose, so a colour
    means the same thing whichever sample you are looking at - which is the
    point of drawing them in one layout.
    """
    E = _entities()
    if E is None:
        return [], {}
    try:
        cols = [c for c in E.columns_of(con, "gene") if c["name"] not in HIDDEN]
        names = [c["name"] for c in cols]
        t = E.select(con, "gene", columns=["sample", *names])
    except Exception:
        return [], {}
    if t.num_rows == 0:
        return [], {}
    rows = {}
    for r in t.to_pylist():
        rows[cohort_key(r["sample"], r["gene_id"])] = r
    # Describe without `sample`, which is constant per run and already carried
    # by the point itself.
    return _describe(t.drop_columns(["sample"]),
                     {c["name"]: c for c in cols}), rows


def _describe(table, help_by: dict) -> list[dict]:
    from collections import Counter

    n = table.num_rows
    out: list[dict] = []
    for field in table.schema:
        name = field.name
        if name in HIDDEN:
            continue
        meta = help_by.get(name) or {}
        vals = table.column(name).to_pylist()
        present = [v for v in vals if v is not None]
        n_null = n - len(present)
        numeric = _is_numeric(field.type) and not isinstance(
            next(iter(present), None), bool)

        d: dict[str, Any] = {
            "name": name, "type": str(field.type),
            "stage": meta.get("stage", ""), "help": meta.get("help", ""),
            "n_null": n_null, "n": n,
        }
        if not present:
            d.update(kind="empty", usable=False, n_distinct=0)
            out.append(d)
            continue

        if numeric:
            lo, hi = min(present), max(present)
            # A column spanning many orders of magnitude - an E-value runs from
            # 1e-29 to 1e-5 - puts every point in the first bin on a linear
            # ramp. Log it when the values allow, and say so in the legend.
            log = bool(lo > 0 and hi > 0 and hi / lo >= 1e3)
            d.update(kind="numeric", n_distinct=len(set(present)),
                     min=lo, max=hi, log=log,
                     usable=lo != hi and not meta.get("identifier"),
                     bins=NUMERIC_BINS)
            out.append(d)
            continue

        counts = Counter(str(v) for v in present)
        order = CANONICAL_ORDER.get(name)
        if order:
            # Declared order first, then anything the declaration did not
            # anticipate, so an unexpected value still appears.
            ranked = [v for v in order if v in counts] + \
                     sorted((v for v in counts if v not in order),
                            key=lambda v: (-counts[v], v))
        else:
            ranked = sorted(counts, key=lambda v: (-counts[v], v))
        d.update(
            kind="categorical", n_distinct=len(counts),
            values=[{"value": v, "n": counts[v]} for v in ranked[:MAX_FILTER_VALUES]],
            truncated=len(counts) > MAX_FILTER_VALUES,
            # More distinct values than slots means the tail is one colour; the
            # legend has to say so rather than imply the plot shows them all.
            folded=max(0, len(counts) - COLOR_SLOTS),
            usable=(len(counts) >= 2 and len(counts) < n * IDENTITY_FRAC
                    and not meta.get("identifier")),
            slots=COLOR_SLOTS,
        )
        out.append(d)
    return out


def color_domain(desc: dict) -> dict:
    """The fixed value -> slot mapping, or the numeric domain, for a column.

    Computed from the whole column so it does not move when the plot is
    filtered.
    """
    if desc["kind"] == "categorical":
        top = [v["value"] for v in desc["values"][:COLOR_SLOTS]]
        return {"kind": "categorical", "slots": top,
                "folded": desc.get("folded", 0), "n_null": desc["n_null"]}
    if desc["kind"] == "numeric":
        lo, hi = desc["min"], desc["max"]
        log = desc.get("log", False)
        edges = _bin_edges(lo, hi, NUMERIC_BINS, log)
        return {"kind": "numeric", "min": lo, "max": hi, "log": log,
                "edges": edges, "n_null": desc["n_null"]}
    return {"kind": "empty"}


def _bin_edges(lo: float, hi: float, bins: int, log: bool) -> list[float]:
    """Equal-width edges, on the log scale when the column is logged.

    Equal width rather than quantiles so the legend reads as even ranges; the
    log transform is what handles the skew that would otherwise pile every
    point into the first bin.
    """
    if log:
        a, b = math.log10(lo), math.log10(hi)
        return [10 ** (a + (b - a) * i / bins) for i in range(bins + 1)]
    return [lo + (hi - lo) * i / bins for i in range(bins + 1)]


def bin_of(value, domain: dict) -> int | None:
    """Which ramp step a numeric value falls in. None for a missing value."""
    if value is None:
        return None
    edges = domain["edges"]
    if domain["log"] and value <= 0:
        return 0
    v = math.log10(value) if domain["log"] else value
    lo = math.log10(edges[0]) if domain["log"] else edges[0]
    hi = math.log10(edges[-1]) if domain["log"] else edges[-1]
    if hi <= lo:
        return 0
    frac = (v - lo) / (hi - lo)
    return max(0, min(NUMERIC_BINS - 1, int(frac * NUMERIC_BINS)))


# --------------------------------------------------------------------------
# filtering
# --------------------------------------------------------------------------
class BadFilter(ValueError):
    pass


def _literal(value, desc: dict) -> str:
    """One SQL literal, typed by the column rather than by what arrived."""
    if desc["kind"] == "numeric":
        try:
            f = float(value)
        except (TypeError, ValueError):
            raise BadFilter(f"{desc['name']} takes a number, not {value!r}")
        if math.isnan(f) or math.isinf(f):
            raise BadFilter(f"{desc['name']} takes a finite number")
        return repr(f)
    if str(desc["type"]).startswith("bool"):
        return "TRUE" if str(value).lower() in ("1", "true", "yes") else "FALSE"
    return "'" + str(value).replace("'", "''") + "'"


def build_where(terms: Sequence[dict], descs: Sequence[dict]) -> str:
    """Structured filter terms -> one SQL predicate.

    Terms are ANDed. A column name is only ever emitted after matching one the
    sample actually has, and every value becomes a literal typed by that
    column, so the browser contributes no SQL text of its own.
    """
    by_name = {d["name"]: d for d in descs}
    parts: list[str] = []
    for t in terms or []:
        if not isinstance(t, dict):
            raise BadFilter("each filter term must be an object")
        name = t.get("col")
        desc = by_name.get(name)
        if desc is None:
            raise BadFilter(f"no column named {name!r} in this sample")
        col = '"' + str(name).replace('"', '""') + '"'
        op = t.get("op") or "in"

        if op == "in":
            vals = t.get("values") or []
            if not vals:
                continue                       # an empty pick is no constraint
            # A null is not "in" anything in SQL, so an explicit null pick is
            # ORed on rather than listed - and it must be kept out of the IN
            # list, or it would be rendered as the literal string 'None' and
            # match nothing while looking like it worked.
            want_null = any(v is None for v in vals)
            present = [v for v in vals if v is not None]
            clauses = []
            if present:
                clauses.append(f"{col} IN ({', '.join(_literal(v, desc) for v in present)})")
            if want_null:
                clauses.append(f"{col} IS NULL")
            parts.append(clauses[0] if len(clauses) == 1
                         else "(" + " OR ".join(clauses) + ")")
        elif op == "range":
            lo, hi = t.get("min"), t.get("max")
            if lo is None and hi is None:
                continue
            if lo is not None:
                parts.append(f"{col} >= {_literal(lo, desc)}")
            if hi is not None:
                parts.append(f"{col} <= {_literal(hi, desc)}")
        elif op == "notnull":
            parts.append(f"{col} IS NOT NULL")
        elif op == "isnull":
            parts.append(f"{col} IS NULL")
        else:
            raise BadFilter(f"unknown filter operator {op!r}")
    return " AND ".join(parts)


def matching(target, sample: str, where: str) -> set[str] | None:
    """gene_ids passing a predicate, or None when it cannot be evaluated.

    None and the empty set mean different things - "cannot filter here" versus
    "nothing matched" - so the caller can tell a sample with no gene rows from
    a selection that excluded everything.
    """
    if not where:
        return None
    E = _entities()
    if E is None:
        return None
    try:
        import lake

        safe = E.guard_predicate(where)          # belt and braces; we built it
        with lake.read(target, budget=5) as con:
            t = E.select(con, "gene", where=safe, columns=["gene_id"],
                         sample=sample)
    except Exception:
        return None
    return set(t.column("gene_id").to_pylist())


def merge(a: Sequence[dict], b: Sequence[dict]) -> list[dict]:
    """Descriptors for two samples reduced to one shared set.

    Only columns both samples have: a colour scale or a filter shown over an
    overlay has to mean the same thing on both sides, and a column one run
    lacks cannot. Domains are unioned rather than taken from either side, so
    the scale covers everything on screen.
    """
    from collections import Counter

    by_b = {d["name"]: d for d in b}
    out: list[dict] = []
    for da in a:
        db = by_b.get(da["name"])
        if db is None or da["kind"] != db["kind"]:
            continue
        d = dict(da)
        d["n"] = da["n"] + db["n"]
        d["n_null"] = da["n_null"] + db["n_null"]
        if da["kind"] == "numeric":
            d["min"] = min(da["min"], db["min"])
            d["max"] = max(da["max"], db["max"])
            d["log"] = bool(d["min"] > 0 and d["max"] > 0
                            and d["max"] / d["min"] >= 1e3)
            d["usable"] = d["min"] != d["max"]
        elif da["kind"] == "categorical":
            counts = Counter()
            for src in (da, db):
                for v in src["values"]:
                    counts[v["value"]] += v["n"]
            order = CANONICAL_ORDER.get(da["name"])
            if order:
                ranked = [v for v in order if v in counts] + \
                         sorted((v for v in counts if v not in order),
                                key=lambda v: (-counts[v], v))
            else:
                ranked = sorted(counts, key=lambda v: (-counts[v], v))
            d["values"] = [{"value": v, "n": counts[v]}
                           for v in ranked[:MAX_FILTER_VALUES]]
            d["n_distinct"] = len(counts)
            d["truncated"] = len(counts) > MAX_FILTER_VALUES
            d["folded"] = max(0, len(counts) - COLOR_SLOTS)
            d["usable"] = len(counts) >= 2 and len(counts) < d["n"] * IDENTITY_FRAC
        out.append(d)
    return out


def slot_of(value, domain: dict) -> int | None:
    """Which colour slot a value takes, given the fixed domain.

    ``None`` is a missing value and ``-1`` is the folded tail; both are drawn
    in neutrals, and both are named in the legend rather than left to be
    guessed at.
    """
    if domain["kind"] == "numeric":
        return bin_of(value, domain)
    if value is None:
        return None
    slots = domain["slots"]
    s = str(value)
    return slots.index(s) if s in slots else -1
