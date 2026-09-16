# Web UI

Progress dashboard for `sae/pipeline`, plus uploading an input and starting a
run. Stdlib only — it adds no entry to `requirements.txt`.

```bash
python sae/web/server.py            # host    -> http://127.0.0.1:8765
container/run.sh web                # container -> http://127.0.0.1:8765
python sae/web/server.py --read-only
```

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

## The feature map

`s06` writes long-format `(gene_id, feature_id, activation)` with top-K per
gene, so a run is a sparse matrix over the 16,384-wide codebook. Expanding a
run plots it in 2D, coloured by the `s05` class — the same kind of picture as
the ESM Atlas map, at a scale that needs no tiling.

```
GET /api/projection?run=<id>&mode=reference|run
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

## Where runs execute

The server does not run the pipeline. It **submits** it, and polls the
scheduler. There is exactly one execution path: `sbatch`. Nothing is ever
forked from the server.

```
GET  /api/state      -> .launcher tells you which backend is in use
POST /api/cancel?id= -> scancel, or terminate for a local fork
```

Jobs survive a restart of the server. Each submission records its identity
next to its log, and on start-up those are adopted and their state refreshed
from the scheduler — otherwise restarting orphans a running job, which looks
exactly like a failure even though the scheduler is still running it.

Forking was wrong on a cluster twice over: the run died with the server, and it
was confined to the *UI's* allocation, so a dashboard sized for browsing could
never start real work. Submitting removes both problems and the sizing question
with them.

There was briefly a local fork as a fallback. It is gone — two execution paths
meant the deployment changed shape depending on where it ran, and the default
quietly chose the weaker one. **Without a scheduler the server still serves the
dashboard** (reading manifests needs nothing) and refuses to launch with a
message saying how to get one. That is a missing capability, not a second code
path.

**Run the server on the login node.** It only reads manifests and submits, so
it needs no allocation — one long-lived lightweight process, with every
expensive thing in its own job.

A submitted job lands on a compute node that has the image but not this
server's interpreter, so by default it runs `container/run.sh` and any path in
its arguments is rewritten to the path the image sees — `/data`, `/work`,
`/atlas`, mirroring run.sh's bind table. `--job-runner python` submits the
interpreter instead, which only works where that path is visible on the node.

`--job-cpus`, `--job-mem`, `--job-time`, `--partition`, `--account` and
`--job-gres` size the jobs. Site-specific options are only sent when set: an
undefined gres or a missing partition is rejected at submission, not later.

### Running SLURM locally

So the deployment does not change shape between a laptop and a cluster:

```bash
container/slurm-local/up.sh                 # SLURM + Apptainer in Docker
eval "$(container/slurm-local/up.sh env)"   # put the shims on PATH
./sae/.venv/bin/python sae/web/server.py    # now submits instead of forking
container/slurm-local/up.sh down
```

The node carries Apptainer as well as SLURM, so a job there executes the same
`sae.sif` a cluster node would — otherwise the local setup would test
scheduling and never execution. The repo is mounted at its own absolute path,
so a path that resolves on the host resolves identically inside a job.

Verified end to end: the server submits, SLURM schedules, the job runs
Apptainer against `sae.sif`, and `s05` writes its output back to the host's
`work/`. One caveat — the toy cluster has no accounting storage, so `sacct`
returns nothing and a finished job reports `returncode: null` rather than 0.
Real clusters have it.

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
* Input paths (including `--hmm` and `--ref`) are resolved and must land under
  an allowed root: uploads, a `--data` directory, or a work root.
* Only the options in `OPTIONS` reach `run.py`; an unlisted key is a 400, so a
  request body cannot introduce a new flag.
* `subprocess` is called with an argument list and never a shell.
* Uploads are capped at 16 GiB and streamed to a `.part` file, renamed only on
  a complete transfer.
* `--read-only` serves progress and rejects every write.

## API

| | |
|---|---|
| `GET /api/state` | runs, jobs, stage list, roots, environment |
| `GET /api/inputs` | files eligible to start a run |
| `GET /api/state` | also lists `hmms` — profile databases found under the data dirs |
| `GET /api/artifacts?run=` | files in each stage directory, with shape and size |
| `GET /api/preview?run=&stage=&file=&limit=` | one artifact as text, table or json |
| `GET /api/log?id=` | tail of a job's log |
| `POST /api/upload?name=` | raw body is the file; no multipart, so no `cgi` |
| `POST /api/run` | JSON `{sample, input_kind, input_path, input_path2?, options}` |
