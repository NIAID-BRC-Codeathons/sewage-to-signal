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

## Host vs container

Paths, interpreter and bind address all differ, and the server detects which it
is in (`/.dockerenv`, `APPTAINER_CONTAINER`, `SINGULARITY_CONTAINER`).

| | host | container |
|---|---|---|
| work roots | `<repo>/work`, `sae/pipeline/work` | `/work` |
| data | `<repo>/data` | `/data` |
| uploads | `<repo>/uploads` | `/work/uploads` |
| interpreter | `sae/.venv/bin/python` | `/opt/venv/bin/python` |
| bind | `127.0.0.1` | `0.0.0.0` under Docker, `127.0.0.1` under Apptainer |

The bind difference is the subtle one. Apptainer shares the host network
namespace, so localhost inside is localhost outside. Docker isolates it, so the
server binds `0.0.0.0` *within its own namespace* and `run.sh` publishes that
to the host's loopback only (`-p 127.0.0.1:8765:8765`) — never `0.0.0.0` on the
host. Override the host port with `SAE_PORT`.

It lives under `sae/` so the image's existing `COPY sae /opt/sae/sae` carries
it in with no extra recipe step, and `SAE_CODE` shadowing covers it too.

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
| `GET /api/log?id=` | tail of a job's log |
| `POST /api/upload?name=` | raw body is the file; no multipart, so no `cgi` |
| `POST /api/run` | JSON `{sample, input_kind, input_path, input_path2?, options}` |
