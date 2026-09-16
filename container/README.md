# Container image

Apptainer/Singularity image for cluster use. It solves the dependency problem
that blocks the pipeline on a bare machine: `megahit`, `mmseqs2` and `fastp`
have no pip equivalent, and installing them locally needed an interactive
`conda tos accept`. In the image they are pinned and present.

```
sae.def               the only recipe
build.sh              build via apptainer / nested / remote builder
run.sh                runtime wrapper that sets the bind mounts
slurm_example.sbatch  example job submission
```

## One recipe, one artifact

There was a `Dockerfile` mirroring `sae.def`, so a Mac without Apptainer could
still build something. It is gone, because the duplication never paid for
itself and quietly cost correctness.

Every commit that ever touched the recipes touched *both* of them — `d00bc12`,
`59e26ac`, `be0936a`, `83780bf`, `fc18e5f`, five for five. There was no case
where they legitimately differed; the README simply asked you to keep them in
step.

Worse, a Docker build cannot exercise Apptainer's semantics, and that is where
the bugs were. `83780bf` is the example: uv installed the interpreter under
`$HOME`, Apptainer bind-mounts the host home over the container's, and the
symlink target vanished. The Docker build passed. The `%test` block passed
*during* the build. Only `apptainer test sae.sif` caught it, and every cluster
run would have failed. A green Docker build was a false signal for the runtime
actually being deployed to.

So the Mac now runs Apptainer *inside* Docker and produces the same `.sif` the
cluster runs. Docker is a host for Apptainer, never a second image format.

## Getting the image

On a cluster, **pull rather than build**. `apptainer build` needs root or
`--fakeroot`, and HPC sites commonly disable fakeroot — so pulling may be the
only route that works there at all.

```bash
container/build.sh pull         # fetch the prebuilt image from GHCR
```

`.github/workflows/container.yml` builds on a native amd64 runner and publishes
to `ghcr.io/niaid-brc-codeathons/sewage-to-signal/sae`, tagged with the commit
SHA and `latest`. `SAE_ORAS` points the pull somewhere else.

The image is 3.21 GB, which rules out the usual GitHub routes — the 100 MB file
limit for the repo itself, and the 2 GB per-file ceiling on both Git LFS and
release assets. GHCR has no such limit and Apptainer speaks ORAS natively, so a
`.sif` is a first-class registry artifact.

Almost all of that size is the CUDA stack, and almost none of it is ours:

| | uncompressed |
|---|---|
| `nvidia/*` (16 CUDA wheels) | 2.7 GB |
| `torch` | 1.1 GB |
| `triton` | 639 MB |
| bioconda (megahit, mmseqs2, fastp) | 595 MB |
| **this project's code** | **188 KB** |

That ratio sets the rebuild policy. CI triggers only on `requirements.txt` and
`sae.def` — the files that change the *environment* — because rebuilding and
pushing 3.2 GB for a 188 KB code change is waste. Code changes ride along at
run time via `SAE_CODE`, which shadows the baked-in copy.

## Build

```bash
container/build.sh              # picks the best available route
container/build.sh apptainer    # native; needs root or --fakeroot, Linux only
container/build.sh nested       # Apptainer inside Docker, for macOS
container/build.sh remote       # Sylabs remote builder
```

Building a SIF creates user namespaces, so the nested route needs
`--privileged`. Without it the build fails at the `%post` scriptlet with
`Failed to create user namespace`. On Docker Desktop that privilege is confined
to its Linux VM, not macOS — but some managed Docker installations forbid
`--privileged` entirely, and there the fallback is `remote`, or building on the
cluster.

`SAE_APPTAINER_IMAGE` overrides the image carrying Apptainer (default
`quay.io/singularity/singularity:v4.1.0`). It holds Apptainer and nothing of
this project, so it is a pinned tool reference, not a recipe to maintain.

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
(`sae/web/`, carried in by the existing `%files sae` entry). Under native
Apptainer the network namespace is shared with the host, so the server binds
localhost and is reachable directly. Nested in Docker it is not, so run.sh
passes `--host 0.0.0.0` and publishes to the host's loopback only,
`-p 127.0.0.1:8765:8765`; `SAE_PORT` changes the host port. Uploads land in
`/work/uploads`, which is your `$SAE_WORK` bind.

### The UI on an HPC node

`container/slurm_web.sbatch` runs the dashboard on a compute node:

```bash
sbatch container/slurm_web.sbatch
tail -f logs/sae-web-<jobid>.out      # prints the exact ssh command to use
```

Compute nodes are not reachable from outside the cluster, so the server binds
loopback on the node and you tunnel to it, terminating the tunnel *on the
node*:

```bash
ssh -N -L PORT:127.0.0.1:PORT -J you@login you@node
```

Binding loopback is not just convention. This server accepts uploads and
launches subprocesses, so on a shared cluster a wider bind hands every other
user on that network a way to run commands as you. If your site refuses ssh
straight to compute nodes, the fallback is forwarding through the login node
against the node's hostname — which does require a wider bind, so pair it with
`--read-only`.

The job picks a free high port rather than assuming 8765, since several people
may do this on one node.

**Sizing is the decision to make.** The server launches pipeline runs as child
processes *inside its own allocation*, so a small allocation means a small
pipeline. Two sensible shapes:

* *Watching* — a modest allocation, `--read-only`, and real work submitted
  separately with `slurm_example.sbatch`. Best for a shared dashboard.
* *Working* — a real allocation (`--gres=gpu:1` for `s06_embed`) and launch
  from the UI. Bounded by the job's wall time.

Submitting a SLURM job per pipeline run from the UI would be the better model
and is not implemented; the server shells out with `subprocess` directly.

### Runtimes

`run.sh` uses **apptainer** or **singularity** when either is on PATH, and
otherwise nests Apptainer in **docker**. `SAE_RUNTIME` forces
`apptainer|singularity|nested`. There is one image either way, `$SAE_SIF`
(default `container/sae.sif`), so every subcommand behaves identically.

The nested path translates binds twice: Docker puts the host paths under
`/mnt`, then Apptainer binds those onto `/data`, `/hf`, `/work`, `/atlas`. It
needs `--privileged` to run a SIF for the same user-namespace reason the build
does, adds `--platform linux/amd64` (`$SAE_PLATFORM`) and `--gpus all` with
`--nv` when a driver is present.

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

## Verification status

Built and tested on macOS/arm64. Apptainer runs inside privileged Docker, so
the artifact tested here is the same `.sif` a cluster would run.

Verified:

* **`container/build.sh nested` builds `sae.def` in a single pass** — one
  `apptainer build`, exit 0, 3.2 GB squashfs. This closes what was previously
  the largest gap: `%post` and `%files` had only ever been exercised through a
  Docker build plus a conversion, never as one command.
* **`apptainer test sae.sif` passes**: interpreter resolves to
  `/opt/uv-python/...` and not through a home mount, `torch 2.11.0+cu130`,
  `esm 3.4.1`, `pyrodigal 3.7.1`, all three binaries (`MEGAHIT v1.2.9`,
  `mmseqs 18.8cc5c`, `fastp 1.3.7`), all seven stages import.
* **Apptainer sections**: `%environment`, `%runscript`, `%test`, `%labels`,
  `%help`, `%files`, and the `%apprun` SCIF apps.
* **`run.sh` nested**: `manifest`, `test`, `pipeline` and `web` all work, with
  binds translated across both layers. A pipeline run wrote
  `work/NESTED/s05_prefilter/` on the host from inside two containers, and the
  web UI was reachable on the host's loopback.
* **`SAE_CODE` shadowing works nested**, so a code change can be tested against
  a built image without a rebuild.
* `torch` carries its own CUDA (`+cu130`), confirming the plain-Ubuntu base
  plus `--nv` is sound.

Not verified, and unavoidable here:

* **GPU execution under `--nv`** — no NVIDIA device on this machine.
* **Native `apptainer build` on Linux.** The nested route is what was
  exercised; a cluster build runs the same definition without the Docker layer.
* **The rendered web UI in a browser.** Every endpoint behind it is tested, but
  the layout and JavaScript are not.
* The SLURM script's partition/account lines are placeholders.

### Bugs this testing caught

Worth recording, because none were visible from inspection:

1. `micro.mamba.pm` unreachable — switched to the GitHub releases static
   binary, which also removes the tar step.
2. `curl` without `-f` wrote an error page that `tar` reported as bzip2
   corruption, hiding a network failure as a format error.
3. `uv venv` installs no `pip`, so `python -m pip freeze` exited 1 and failed
   the build; `2>/dev/null` hid the message but not the status.
4. **`uv` installs its managed interpreter under `$HOME`.** Apptainer
   bind-mounts the host home over the container's, hiding the target of
   `/opt/venv/bin/python`. The `%test` block passed at build time but
   `apptainer test` failed with `python: not found`, and every `run.sh` call on
   a cluster would have failed the same way. Docker does not mount over the
   home, so **only the Apptainer test could find this.** Fixed with
   `UV_PYTHON_INSTALL_DIR=/opt/uv-python`; `%test` now resolves the symlink and
   fails loudly if it regresses.
5. **Nested, Apptainer shares the *Docker container's* network namespace**, not
   the host's, so the web UI binding localhost was unreachable through `-p`.
   `run.sh` passes `--host 0.0.0.0 --published` for that case only; the server
   warns about a non-loopback bind unless something in front of it is known to
   control exposure.

### Two Pythons in the image

`/opt/venv` holds ours (3.12.14). The bioconda environment brings its own
(3.14.7) as a dependency of the three tools. `PATH` puts `/opt/venv/bin` first,
so `python` resolves to 3.12.14 — verified. Both appear in the image manifest,
which is expected rather than a packaging error.

## Worth adding on a GPU cluster

`flash-attn` and `transformer-engine` are CUDA-only, so they could not be
installed locally — the `esm` library warns about their absence on every run.
On a CUDA cluster they are installable and would address the throughput ceiling
measured locally (1.03 seq/s at batch 1, rising only to 1.24 at batch 8, because
the 6B model is memory-bandwidth-bound without fused kernels). They are
deliberately not in the image: both are slow and fragile to build, and pinning
them against the cluster's exact CUDA version is better done once you know it.
