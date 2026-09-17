# Web UI

Progress dashboard for `sae/pipeline`, plus uploading an input and starting a
run. Stdlib only — it adds no entry to `requirements.txt`.

**It has no list of stages.** The launch form, the parameter fields, the stage
strip, the plan preview and the predicate boxes are all generated from
`GET /api/pipeline`, which is the pipeline describing itself. A stage that is
added, renamed or removed shows up with no change here, and a work directory
written by a different pipeline renders from its own manifests.

```bash
python sae/web/server.py            # host    -> http://127.0.0.1:8765
container/run.sh web                # container -> http://127.0.0.1:8765
python sae/web/server.py --read-only
```

## Where the data is

`data/sae.ducklake`. Progress, columns, predicates and the feature map are all
queries against it, so the server holds no stage list and reads no manifests.
It attaches **read-only and briefly** for each request: a local DuckLake catalog
is a DuckDB file, so a held connection would lock out a running job. Work
directories still exist for the file-shaped stages (FASTQ, contigs, logs) and
the artifact browser still lists them.

## Where progress comes from

Nothing was added to the pipeline to support this. Every stage already writes
`<output>.manifest.json` recording inputs, params, stats, timing and tool
versions, so the dashboard is a read over the work roots. Runs started from the
CLI, from a different machine, or inside the container all appear.

Two sources, because neither alone is enough:

| | supplies | limitation |
|---|---|---|
| manifests | durable state; a stage is done when its manifest exists | only written on completion |
| job log | live state — `run.py` prints `[stage]` per stage | only for runs started here |

So a run started from the CLI shows accurate completed stages but no live
cursor; a run started from the UI shows both.

## Building the form from the pipeline

`GET /api/pipeline` returns every stage with its parameters — type, default,
choices, help text, whether the value is a path and which suffixes it takes —
along with the entity levels and the roles each stage fills. The frontend turns
that into controls:

* a **plan strip** showing which stages a given input and target imply;
* a **fieldset per stage**, with its summary, its missing tools, and one
  control per declared parameter;
* a **selection box** for every stage that says it takes a predicate, showing
  the stage's own default as the placeholder.

Typing a predicate against a sample that already has rows queries
`GET /api/query`, so the match count updates as you type — you see what a
selection takes before spending a GPU on it. Cross-sample predicates work here
too, because sibling samples are registered as SQL schemas:
`seq_sha1 IN (SELECT seq_sha1 FROM "CHI-A".gene WHERE category='dark')`.

Capabilities are looked up by **role**, never by stage name. The feature map is
drawn when some stage declaring `projection` has run; the point colours come
from whatever wrote the `category` column. A pipeline that fills those roles
differently still gets a map.

## Viewing stage outputs

Expanding a run lists the files in each stage directory; clicking one previews
it. The point of the design is that **the frontend never learns a stage's
schema.** The server normalises every artifact into one of three shapes:

| shape | from | rendered as |
|---|---|---|
| `table` | `.parquet`, `.tsv`, `.csv` | columns + rows, whatever they are |
| `json` | `.json` | pretty-printed |
| `text` | `.faa`, `.fastq`, `.txt`, `.md`, `.log`, … | head of the file |

So `s07_match` can add, drop or rename a column in `clusters.parquet` and the
UI keeps working — it renders whatever columns come back. A new stage needs no
frontend change either. And a stage that wants a curated summary rather than
its raw output just drops a `.md` or `.tsv` into its work directory; it appears
automatically, with arbitrary content.

`.gz` is transparent, and the suffix underneath decides the shape, so
`foo.tsv.gz` is still a table. Anything unrecognised is reported as `binary`
with its size and no preview. Previews are bounded — 200 lines or 100 rows by
default, `&limit=` to raise it, capped at 5000 — so opening a 30 GB FASTQ is
cheap.

Parquet needs `pyarrow`, which is a *pipeline* dependency, not a web one. It is
imported only when a parquet is actually requested, so the server still runs
where it is absent; you get a message in place of the table.

## The atlas

Every embedding in the store in one layout, drawn colourless, with one sample
lit at a time. It answers a question the per-run map cannot: where a sample sits
in *everything that has been embedded*, rather than where its own proteins sit
relative to each other.

`GET /api/atlas?color=<column>&limit=` returns the whole cohort once — each
point carrying its sample index, its value for the chosen column, and the colour
slot the server resolved. **Hovering is done in the browser**: the backdrop is
drawn once and a hover rewrites only the highlight layer, so lighting a sample
costs no request and the layout cannot shift under the cursor.

The expensive half — reading every activation and running UMAP over it — is
cached on the store's snapshot id, so changing the colour column is ~100 ms
against ~10 s for the first build. When the snapshot moves the view says it is
stale and offers a rebuild rather than taking one: an automatic rebuild would
fire on every stage a running batch finishes, and throw away the highlight you
were reading.

Thinning is per sample, not over the cohort, so a nine-gene sample still appears
next to a twenty-thousand-gene one. Points are picked at even spacing rather
than by a stride, because a stride can only halve: asking for 93% of a sample
would otherwise hand back 50%.

A sample with embeddings but no gene rows — anything backfilled from a work
directory that only ran s06 — draws hollow, because absence is not a colour.

**Concurrency note.** UMAP runs on numba, whose default threading layer is not
threadsafe: two projections at once terminate the process rather than merely
contending. This is a `ThreadingHTTPServer`, so two tabs are enough. Every
layout goes through one lock.

## The feature map

`s06` writes long-format `(gene_id, feature_id, activation)` with top-K per
gene, so a run is a sparse matrix over the 16,384-wide codebook. Expanding a
run plots it in 2D, coloured by the `s05` class — the same kind of picture as
the ESM Atlas map, at a scale that needs no tiling.

```
GET /api/projection?run=<id>&mode=reference|run[&b=<id>]
```

**UMAP only.** It is the one method here with a `transform`, and that is what
makes a fixed layout possible — t-SNE cannot place a new point in an existing
layout without refitting, and a linear projection was not good enough.

### The shared layout

Projecting each run on its own gives an arbitrary layout: two samples cannot be
compared, and re-running moves every point. So UMAP is fitted **once** over a
reference corpus and every run is `transform`ed into that space, exactly as the
ESM Atlas serves precomputed `umap_1`/`umap_2` columns.

```bash
python sae/web/reference_map.py build \
    --out data/reference_map.joblib \
    "work/*/s06_embed/*.sae_features.parquet"

python sae/web/reference_map.py info        # corpus, params, library versions
```

The run's proteins are drawn over the corpus, which appears as faint context.
Without a map the server falls back to fitting each run alone and labels it
plainly, because those coordinates mean something different.

The corpus has to span what you expect to see: `transform` places a point
relative to the fitted manifold and has nothing useful to say about a region
the fit never covered. That is why it wants both ends — proteins you care about
*and* enough empirical wastewater to cover the unannotated bulk.

Two honest limits. `transform` is an approximation, so a protein in the corpus
does not land exactly where the fit put it (median drift 0.73 on a span of 18.4
in our build, ~4%). And the map is a pickle, so it is version-sensitive; the
build records `umap`, `numpy`, `scipy` and `sklearn` versions and the response
flags any drift rather than quietly returning a wrong layout.

**Everyone sharing a map must share the file.** Two people who each fit their
own have incomparable coordinates, which is the problem this exists to solve.
It is gitignored on purpose: our corpus today is "whatever is in `work/`",
which is not reproducible, so committing one would enshrine an arbitrary
sample. A corpus worth sharing should be defined first.

### Two samples at once

`&b=<run id>` draws a second run's proteins in the same layout — the **compare
with…** picker beside the mode buttons. This is what the fixed layout is *for*:
two samples only overlay meaningfully if their coordinates already mean the
same thing, which is exactly the property a per-run fit destroys.

The overlay colours by **sample** rather than by s05 class. Two samples times
three classes is six series, which a scatter of several thousand points cannot
carry; the class stays in the hover text, so it is demoted rather than lost.
The two sample colours are a separate pair from the class colours for the same
reason — reusing "known" blue for sample A would make one swatch mean two
different things across two modes of the same plot.

Points are merged in proportion rather than drawn sample-by-sample. Drawing all
of A and then all of B buries A under whichever sample is larger, and the
picture then reads as "B is everywhere" regardless of what is there.

With `mode=run` and a `b`, the fallback fits UMAP over **both** samples
together rather than over one. That is a real comparison — the two sit in one
space — but a private one: those coordinates match no other run and not the
shared map either, and the panel says so. It is the honest option when there is
no reference map to borrow.

Two runs can double the point count, so the total drawn is capped
(`MAX_PROJECTION_POINTS`) by a deterministic stride, and the panel says when it
thinned. A stride rather than a random draw, so re-opening a comparison shows
the same picture instead of reshuffling.

A run has to have reached `s06_embed` to appear in the picker, and a run cannot
be compared with itself.

### Colour and filter by anything a stage wrote

The gene level is the join of every stage's column fragments, so a gene already
carries coordinates and length from s03, cluster membership from s04, the
homology triage from s05 and embedding status from s06. That is the metadata,
and the map colours and filters by any of it:

```
GET /api/projection?run=<id>&color=<column>&filter=<json terms>
```

**The page has no list of columns.** The response carries the columns the run
actually has, with their types, value counts and ranges, and the controls are
built from that — so a stage added tomorrow that annotates genes with, say, a
taxon call becomes another thing to colour by with no change to the server or
to `index.html`. It is the same property the artifact browser has, applied to
the map.

A column is offered for **colouring** only if it discriminates: an identifier
gives every protein its own colour and a constant gives them all one, so
`gene_id`, `seq_sha1` and `rep_id` are filterable but not colourable.

#### Three colours, and what happens past three

A scatter puts every pair of series on screen at once, and only **three** hues
clear the colour-blind separation floors under that condition — a fourth cannot
(`scripts/validate_palette.js --pairs all`). So a categorical column shows its
three commonest values and folds the rest into one neutral, which the legend
names as `other (n)` rather than implying the plot shows them all. Colouring
1474 proteins by `family` — 195 distinct values — is still useful; it just
answers "where are the three big families" and says so.

Two things a colour cannot express are drawn rather than left out. A **missing
value** is a hollow ring, at reduced strength so that absence recedes instead
of outdrawing the findings; 71% of these proteins have no `family`, and at full
strength the rings were the loudest thing on the plot. A **numeric** column
gets a five-step single-hue ramp, log-scaled when the values span three orders
of magnitude or more — an E-value runs from 1e-158 to 1e-5 here, and on a
linear ramp every point lands in the first bin.

The value → colour mapping is computed over the **unfiltered** column and held
fixed, so narrowing the plot never repaints the points that survive. The legend
counts, by contrast, are counted over the points actually drawn: the gene level
holds every gene but the plot holds only the embedded ones, which is why
`category` legends read `known (0)` — s06 does not embed a protein a family
already explains.

#### Filtering

Filter controls are built from the same descriptors: a value picker for a
categorical column, a min/max pair for a numeric one, one term per column,
ANDed. **No predicate text crosses the wire.** The browser sends structured
terms — column, operator, values — and the SQL is composed server-side against
the column list the sample actually has, with every literal typed by its own
column, so there is nothing to escape at the boundary; the result still goes
through `entities.guard_predicate`. The same selection written by hand is what
`run.py --where` takes, and it means the same thing.

Filtering is a display operation, not a re-projection: everything is projected
and then filtered, so the layout does not move when the filter changes. That
matters most in the fallback modes, where filtering first would refit the
layout around whatever survived.

When two runs are overlaid, the colour and filter columns are the ones **both**
have — a scale shown over two samples has to mean the same thing on both — and
the domains are unioned so the scale covers everything on screen. Colour still
defaults to sample there; choose a column instead and run B keeps a dark ring,
which reads at a 2.6px mark where a different shape does not.

A work directory written before the column fragments existed has no gene level.
It falls back to whatever s05 wrote beside its output, which yields `category`
and nothing else — enough to keep the map coloured as it always was, rather
than going blank on a run nobody has re-run.

### Reading it

Rows are L2-normalised before fitting: otherwise a protein's activation
*magnitude* dominates the leading components and everything else collapses
toward the origin. The response reports how many codebook features are *shared*
between proteins, which is the health check — with top-16 over 16,384, two
proteins may share none, and then the layout is noise rather than biology.

Measured on 1142 CASPER proteins plus the SARS-CoV-2 reference: 1033 distinct
features, 476 shared, fitted in 5.2 s, 0.2 MB on disk. The classes separate in
the shared layout — 10-nearest-neighbour same-class rate 0.612 against 0.461
expected by chance — along a known → partial → dark gradient.

### Two things the map had to learn

**Projections are serialised.** UMAP is numba, and numba's default
`workqueue` threading layer is not threadsafe — called from two Python threads
at once it does not raise, it aborts the process. This is a
`ThreadingHTTPServer`, so two overlapping projection requests took the whole
dashboard down, which became easy to hit once a run could be drawn against a
second one. They now run under one lock. That costs nothing real: a projection
is CPU-bound and gains nothing from running beside another, and the worst case
is a request that waits instead of a server that dies.

**The plan preview asks rather than imitates.** The page used to resolve the
stage list itself by walking ports. That held only while one stage consumed
each port; the moment two did — `s02_assemble` and `s02_translate` both take
reads — the greedy walk returned the *union* of both routes and could not tell
that skipping the assembler leaves no route to contigs at all. `/api/plan` runs
the driver's planner, so the preview is what will run.

## Where runs execute

The server runs the pipeline here, on this machine, **one run at a time**. Each
job is detached into its own session, so it outlives the server rather than
dying with it.

```
GET  /api/state      -> .launcher tells you which backend is in use
GET  /api/script?id= -> the script exactly as it was run
POST /api/cancel?id= -> SIGTERM to the job's process group
```

Every job is written as a shell script next to its log, limits included, and
that file — read back from disk, not reconstructed — is what the Jobs panel
shows beside the job's output. So a job is answerable after the fact: what it
asked for, what it ran, and what it printed, rather than what the server would
run today under whatever settings it now has.

The **Runs** panel stays output-only. It reads manifests, which is how runs
started outside this UI appear at all, and those have no script to show.

### Serial, and why that is the whole scheduler

One run at a time, and it sees every GPU on the box. On a single host the
useful question is not which GPU a run gets but whether two runs are competing
for the same one, and a queue of one answers it without a slot table to keep
correct. A second submission waits, visibly, as `pending`.

For a cohort rather than a single run, `container/run_batch.sh` is the other
trade: it fans samples out across all the GPUs at once, one sample pinned per
GPU. Use the UI to watch one run; use the batch script to process many.

### Surviving a restart

Jobs survive a restart of the server, which is the one property the scheduler
used to provide for free. Two mechanisms replace it.

`start_new_session` puts each job in its own session, so killing the server
leaves it running and `killpg` can still take down the whole tree — the
container child included, not just the wrapper shell.

The job records its own exit code in `<id>.rc`. That file is the local
stand-in for `sacct`: after a restart there is no parent left to reap the
process, so without it a finished job could only be reported as "gone", which
looks exactly like a failure. On start-up the records next to the logs are
adopted and their state recovered — from the `.rc` file if it is there, from
the pid if it is not.

| SLURM used to | Now |
|---|---|
| `squeue` says live | no `.rc` yet **and** the pid answers |
| `sacct` gives State/ExitCode | `.rc` holds the exit code |
| no accounting plugin → finished, code unknown | no `.rc`, pid gone → the same |

That last row is a real gap, not a tidy one: a job whose server died *and*
which never wrote its code is reported finished with `returncode: null`. It is
the same imprecision the SLURM backend accepted on sites without accounting
storage, and it is preferable to the alternative of calling a job failed
because nobody was watching when it ended.

A job the shell reaped is recorded as `128+N` when it dies by signal; one the
server reaped is recorded as `-N`, which is more precise. Both decode to the
same hint. The ambiguity only bites a job that exits 137 on purpose, which no
stage does.

### One execution path, still

The server used to submit every run with `sbatch`, and before that it forked
them as its own children. Forking was wrong on a cluster twice over: the run
died with the server, and it was confined to the *UI's* allocation, so a
dashboard sized for browsing could never start real work.

The deployment is one server now. The allocation half of that argument has no
target — there is no allocation — but the first half stood, and is what the
detaching and the `.rc` file are for.

What has not changed is that there is exactly **one** backend. A local fork
once existed as a *fallback* beside `sbatch`, and that was the actual mistake:
two paths meant the deployment changed shape depending on where it ran, and the
default quietly chose the weaker one. Replacing the backend keeps that
property; adding a second one would not.

### Sizing

`--job-cpus` is exported to each job as `OMP_NUM_THREADS` — without a scheduler
nothing enforces a core count, so the flag sets the one knob that actually
reaches the work. `--job-time` is enforced by the server: over it, the job is
terminated.

There is no `--job-mem`. Memory cannot be capped without cgroups, and a flag
that reports a limit nothing applies is worse than no flag — an OOM kill would
look like the limit working.

By default a job runs `container/run.sh`, and any path in its arguments is
rewritten to the path the image sees — `/data`, `/work`, `/atlas`, mirroring
run.sh's bind table. `--job-runner python` runs this server's interpreter
instead.

## Host vs container

Paths, interpreter and bind address all differ, and the server detects which it
is in (`/.dockerenv`, `APPTAINER_CONTAINER`, `SINGULARITY_CONTAINER`).

| | host | container |
|---|---|---|
| work roots | `<repo>/work`, `sae/pipeline/work` | `/work` |
| data | `<repo>/data` | `/data` |
| uploads | `<repo>/uploads` | `/work/uploads` |
| interpreter | `sae/.venv/bin/python` | `/opt/venv/bin/python` |
| bind | `127.0.0.1` | `127.0.0.1`; run.sh passes `0.0.0.0` when nested in Docker |

The bind difference is the subtle one. Apptainer shares the host network
namespace, so localhost inside is localhost outside and the default is right.
Nested in Docker it shares the *Docker container's* namespace instead, which
`-p` cannot reach, so `run.sh` passes `--host 0.0.0.0` and publishes that to
the host's loopback only — never `0.0.0.0` on the host. `SAE_PORT` sets the
host port.

It lives under `sae/` so the image's existing `%files sae /opt/sae/sae` entry
carries it in with no extra recipe step, and `SAE_CODE` shadowing covers it.

## Security

This accepts uploads and launches subprocesses, so it is a development tool,
bound to localhost, and should not be exposed. Beyond that:

* Sample names and uploaded filenames must match `[A-Za-z0-9][A-Za-z0-9._-]{0,63}`
  — they become path components, so they are restricted rather than escaped.
* Any parameter a stage declared as a path is resolved and must land under an
  allowed root: uploads, a `--data` directory, or a work root.
* A parameter reaches `run.py` only if some stage declares it, and only after
  that stage's own `Param` has coerced it; an unknown stage or key is a 400, so
  a request body cannot introduce a new flag.
* Predicates are checked by `entities.guard_predicate` before being forwarded
  or evaluated: one expression, no `;`, and none of the statements that would
  write, attach or install. They are evaluated against a throwaway in-memory
  DuckDB connection whose only tables are read-only views over parquet. This
  matters because a `WHERE` clause can otherwise reach `COPY ... TO`, and this
  server writes files.
* `subprocess` is called with an argument list and never a shell.
* Uploads are capped at 16 GiB and streamed to a `.part` file, renamed only on
  a complete transfer.
* `--read-only` serves progress and rejects every write.

## API

| | |
|---|---|
| `GET /api/state` | runs, jobs, stage list, roots, the lake, environment; `data_files` lists files matching any path parameter's declared suffixes |
| `GET /api/pipeline` | the whole pipeline: stages, parameters, levels, roles, tool availability |
| `GET /api/columns?run=&level=` | a level's columns, with which stage wrote each |
| `GET /api/query?run=&level=&where=&limit=&scope=` | how many rows a predicate selects, plus a look at them; `scope=cohort` drops the sample filter |
| `GET /api/inputs` | files eligible to start a run |
| `GET /api/plan?have=&want=&skip=` | the stages that would run, from the driver's own planner |
| `GET /api/artifacts?run=` | files in each stage directory, with shape and size |
| `GET /api/projection?run=&mode=&b=&color=&filter=` | the feature map: layout, colour domain, filterable columns |
| `GET /api/preview?run=&stage=&file=&limit=` | one artifact as text, table or json |
| `GET /api/log?id=` | tail of a job's log |
| `GET /api/script?id=` | the script the job was run as |
| `POST /api/upload?name=` | raw body is the file; no multipart, so no `cgi` |
| `POST /api/run` | JSON `{sample, input_kind, input_path, input_path2?, target?, params?, where?, only?, skip?, force?}` — `params` and `where` are keyed by stage name |
