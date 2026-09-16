# Container image

Apptainer/Singularity image for cluster use. It solves the dependency problem
that blocks the pipeline on a bare machine: `megahit`, `mmseqs2` and `fastp`
have no pip equivalent, and installing them locally needed an interactive
`conda tos accept`. In the image they are pinned and present.

```
sae.def               Apptainer definition — the primary artifact
Dockerfile            same recipe, for building without root and local testing
build.sh              build via apptainer / docker+convert / remote builder
run.sh                runtime wrapper that sets the bind mounts
slurm_example.sbatch  example job submission
```

## Build

```bash
container/build.sh              # picks the best available route
container/build.sh apptainer    # native; needs root or --fakeroot, Linux only
container/build.sh docker       # build with Docker, convert to SIF
container/build.sh remote       # Sylabs remote builder
```

Apptainer cannot build on macOS. From a Mac, use the `docker` route and convert
on a host that has Apptainer, or build on the cluster directly.

## Run

One command runs all seven stages; `run.sh` sets the bind mounts and passes
everything else through to `pipeline/run.py`.

```bash
container/run.sh pipeline --fastq /data/fastq_rnaseq/SRR38294894_1.fastq.gz \
                          --fastq2 /data/fastq_rnaseq/SRR38294894_2.fastq.gz \
                          --sample CHI-A
container/run.sh web          # progress UI, published to host loopback
container/run.sh query --fasta /data/contig.fna --top-k 8
container/run.sh manifest     # exact package versions baked into this image
container/run.sh test         # self-check
container/run.sh shell        # interactive
```

### The web UI

`container/run.sh web` serves the dashboard from inside the image
(`sae/web/`, carried in by the existing `COPY sae`). Under Docker the server
binds `0.0.0.0` in its own network namespace and run.sh publishes it to the
host's loopback only, `-p 127.0.0.1:8765:8765`; `SAE_PORT` changes the host
port. Under Apptainer the network namespace is shared, so it simply binds
localhost. Uploads land in `/work/uploads`, which is your `$SAE_WORK` bind.

### Runtimes

`run.sh` picks the first of **apptainer**, **singularity**, **docker** on PATH;
`SAE_RUNTIME=docker` overrides. The subcommands are identical across runtimes.

| | image | selected by |
|---|---|---|
| apptainer / singularity | `$SAE_SIF` (default `container/sae.sif`) | the cluster path |
| docker | `$SAE_IMAGE` (default `wastewater-sae:latest`) | local, where Apptainer is unavailable |

The docker branch mirrors the `%apprun` entrypoints in `sae.def` rather than
relying on the Dockerfile's `CMD`, so `pipeline`, `query`, `manifest`, `test`
and `shell` behave the same either way. It adds `--platform linux/amd64`
(`$SAE_PLATFORM`), `--gpus all` when a driver is present, and on Linux
`--user $(id -u):$(id -g)` so runs do not leave root-owned files in `/work`.
Docker Desktop maps ownership itself, so that flag is skipped on macOS.

**On Apple Silicon this runs under emulation.** The image is `linux/amd64`, so
`s06` executes on emulated CPU with no MPS — much slower than the native venv.
The container's value on a Mac is narrow but real: it is the only way past
`s02_assemble` without a local `megahit`. Assembling in the container and
returning to the native venv for `s06` is the sensible split.

**The image is built from `git archive HEAD`**, so uncommitted work is not in
it. Either commit first, or shadow the baked-in code with
`SAE_CODE=$PWD/sae` — that works on both runtimes.

## What is in the image, and what is not

**In:** Ubuntu 24.04, Python 3.12 (provisioned by uv, not the base image's
python3), all of `requirements.txt`, and `megahit=1.2.9`, `mmseqs2=18.8cc5c`,
`fastp=1.3.7` from bioconda. Also the pipeline code, so the image is a complete
artifact — `SAE_CODE=/path/to/sae` shadows it with a working copy when you want
to iterate without rebuilding.

**Out:** model weights (~27 GB) and data (~29 GB). Baking them would produce an
unusable image and they are shared between runs anyway. They are bind mounts:

| Inside | Variable | Default |
|---|---|---|
| `/data` | `SAE_DATA` | `<repo>/data` |
| `/hf` | `SAE_HF` | `${HF_HOME:-$HOME/.cache/huggingface}` |
| `/work` | `SAE_WORK` | `$PWD/work` |
| `/atlas` | `SAE_ATLAS` | `<repo>/sae` |

Put `SAE_HF` on shared scratch rather than `$HOME` — it is ~27 GB and every job
reuses it.

## Two design decisions

**No CUDA base image.** PyTorch's Linux wheels bundle their own CUDA runtime
(the `nvidia-*-cu12` dependencies), so the image needs only the host driver,
which `apptainer run --nv` provides. That keeps the image smaller and avoids
pinning a CUDA base tag that has to match whatever driver the cluster runs.
`run.sh` probes for `nvidia-smi` and adds `--nv` only when a driver is present,
so the same image works on CPU and GPU nodes.

**The image is the lockfile.** `/opt/image-manifest.txt` records the resolved
bioconda and pip versions at build time, readable with
`container/run.sh manifest`. `requirements.txt` pins Python packages, but it is
not cross-platform — `torch==2.11.0` resolves to a CUDA build on Linux and a
CPU/MPS build on macOS. The image is what makes a run reproducible, not the
requirements file alone.

`setup.sh` detects the container (via `APPTAINER_CONTAINER`,
`SINGULARITY_CONTAINER` or `/.dockerenv`) and skips venv creation, reporting the
baked-in environment instead. Data targets still work, so you can fetch inside
or outside.

## Verification status — read this before trusting the recipe

Verified on this machine (macOS/arm64, Docker only):

* bioconda resolves and installs the three tools on `linux/amd64`, and all
  three binaries execute: `MEGAHIT v1.2.9`, `mmseqs 18.8cc5c`, `fastp 1.3.7`.
  The pins in `sae.def` and `Dockerfile` are those verified versions.
* `setup.sh` container detection works — it reports the baked-in environment
  and stops asking for uv.
* Shell syntax of every script here.

**Not verified here**, because this machine has no Apptainer and no NVIDIA GPU:

* `apptainer build` of `sae.def` end to end, and the `%test` section.
* The Docker→SIF conversion.
* GPU execution under `--nv`, and CUDA-enabled `torch` at all.
* The SLURM script, whose partition/account lines are placeholders.

So treat the first build on the cluster as the real test. Run
`apptainer test sae.sif` immediately after — the `%test` section checks the
Python imports, the three binaries, and that all seven pipeline stages import.

## Worth adding on a GPU cluster

`flash-attn` and `transformer-engine` are CUDA-only, so they could not be
installed locally — the `esm` library warns about their absence on every run.
On a CUDA cluster they are installable and would address the throughput ceiling
measured locally (1.03 seq/s at batch 1, rising only to 1.24 at batch 8, because
the 6B model is memory-bandwidth-bound without fused kernels). They are
deliberately not in the image: both are slow and fragile to build, and pinning
them against the cluster's exact CUDA version is better done once you know it.
