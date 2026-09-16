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
| `GET /api/artifacts?run=` | files in each stage directory, with shape and size |
| `GET /api/preview?run=&stage=&file=&limit=` | one artifact as text, table or json |
| `GET /api/log?id=` | tail of a job's log |
| `POST /api/upload?name=` | raw body is the file; no multipart, so no `cgi` |
| `POST /api/run` | JSON `{sample, input_kind, input_path, input_path2?, options}` |
