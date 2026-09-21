# ab_initio_pipeline

A generalised, configure-and-run template for **ab-initio complex modelling
with [ABCfold](https://github.com/rigdenlab/ABCFold)**, using Snakemake for
all three stages and a pluggable Snakemake executor for stage 2's cluster
jobs (SLURM / PBS / HTCondor). Distilled from three one-off sibling
pipelines (NPF transporters; DRB2/DCL4/dsRNA; AGO1/FBW2/miR) that all ran
the same 3-stage protocol and diverged only in hand-edited scripts.

Stage 2 was built and validated against **SLURM on the IFB Core Cluster**
(France Bioinformatique) — that's the only scheduler this pipeline has
actually been run against. PBS and HTCondor are supported the same way
config-wise, but untested — see "Non-SLURM schedulers" below before relying
on either.

Adding a new dataset should require editing **only**:

1. one line in [`config.yaml`](config.yaml)'s `systems:` list, and
2. a new `configs/<system>.yaml` chain spec (schema: [`configs/_schema.md`](configs/_schema.md)).

No script edits.

---

## The three stages

```
┌─ STAGE 1  preprocessing (local, needs internet) ──────────────────────────┐
│  workflows/preprocessing/Snakefile          snakemake --use-conda --cores 2 │
│                                                                            │
│  fold_input.json ─► MMseqs2 MSA+templates ─► fold_input.resolved.json       │
│  sequence annotation:  domains (ScanProsite + InterProScan REST)            │
│                        disorder (IUPred3 long/short + ANCHOR2 + AIUPred)     │
│                        MoRF (MoRFchibi2, tools/MC2)                          │
│                     ─► data/annotation/<sys>/annotation.yaml  [HUMAN REVIEW] │
│  ─► reports/preprocessing.zip  (unzip → preprocessing/report.html)          │
└────────────────────────────────┬───────────────────────────────────────────┘
                     rsync data/fold_inputs/  ►  IFB
┌─ STAGE 2  processing (IFB login node) ─────────────────────────────────────┐
│  workflows/processing/Snakefile   snakemake --profile workflows/processing/profiles/ifb │
│                                                                            │
│  one SLURM job / system:  abcfold -abcopr  (6 backends together)            │
│  ─► metadata compression (arrays.h5 + model_metadata.parquet)              │
│  ─► reports/processing.zip                                                  │
└────────────────────────────────┬───────────────────────────────────────────┘
              rsync results/abcfold/ results/metadata/  ◄  IFB
┌─ STAGE 3  postprocessing (local) ─────────────────────────────────────────┐
│  workflows/postprocessing/Snakefile         snakemake --use-conda --cores 4 │
│                                                                            │
│  anchor rigid-core Kabsch alignment (one per anchor_chains entry)          │
│  interface metrics over the FULL ensemble:  ipSAE · iLIS · Pinc             │
│  pose clustering (HDBSCAN + Optuna/DBCV, or GMM/BIC)                        │
│  top-N/cluster ─► ChimeraX minimize ─► fix_pdb ─► PLIP ×passes ─► aggregate │
│  ─► reports/postprocessing.zip   (dim-reduction panels, PLIP heatmaps, …)   │
└────────────────────────────────────────────────────────────────────────────┘
```

Each stage's **default target builds a self-contained report bundle**
(`reports/<stage>.zip`) via an `onsuccess:` hook — unzip it and open
`<stage>/report.html`. It ships as a `.zip` rather than a bare `.html` so the
interactive datavzrd tables can be embedded (browsers block the `data:` URLs
Snakemake would otherwise use for embedded HTML). A dark theme is the default,
from `report/custom.css` (injected via `--report-stylesheet`); figures are
SVG/PDF only, never raster.

---

## Running on IFB vs. another SLURM cluster

This pipeline was developed and validated end-to-end on the **IFB Core
Cluster** (`core.cluster.france-bioinformatique.fr`, France
Bioinformatique). If that's you, `config.yaml`'s `hpc.slurm:` block and
`abcfold.model_params` already carry the tested values — you only need your
own `account` (prompted at generation time, or edit `hpc.account` directly)
and the `workflows/processing/profiles/ifb/` executor profile works as
shipped.

If you're on a **different SLURM cluster**, nothing about the pipeline
logic assumes IFB — only these values do, and they're plain `config.yaml`
keys, not build-time template variables, so just edit them directly for
your site before running stage 2:

| Key | Where | IFB value (yours will differ) |
|---|---|---|
| `hpc.slurm.partition` | `config.yaml` | `"gpu"` |
| `hpc.slurm.exclude_nodes` | `config.yaml` | `"gpu-node-7,gpu-node-9"` (known-bad IFB nodes — set `""` elsewhere) |
| `hpc.cuda_module_version` | `config.yaml` | `"12.9.1"` (IFB's `cuda-toolkit` module) |
| `hpc.slurm.size_classes.*.gres` | `config.yaml` | `"gpu:l40s:1"` / `"gpu:h200:1"` (IFB GPU names) |
| `abcfold.model_params` | `config.yaml` | `/shared/bank/alphafold3/current` (IFB's AF3 weights) |
| `workflows/processing/profiles/ifb/` | executor profile | rename/copy for your own site if its defaults don't fit |

## Non-SLURM schedulers

`scheduler` is a cookiecutter prompt (`slurm` / `pbs` / `htcondor`) that
sets `config.yaml`'s `hpc.scheduler`, which in turn picks which block of
`hpc:` the stage-2 Snakefile (`workflows/processing/Snakefile`) reads to
build each `run_abcfold` job's resource request, and which
`workflows/processing/profiles/<scheduler>/` you point `--profile` at.
Switching later (not just at generation time) is just editing
`hpc.scheduler` and using the matching `--profile` — no code changes.

**Only `slurm` is validated** — the whole pipeline was developed and run
end-to-end on it. `pbs` and `htcondor` are wired up against each Snakemake
executor plugin's own documented resource keys (see the comments in
`config.yaml`'s `hpc:` block and in each profile), but neither has been run
against a live PBS/Torque or HTCondor pool. Known gaps if you try one:

- **PBS** — via the community `snakemake-executor-plugin-pbs`
  (not yet on PyPI as of 2026-09; install command is in `envs/controller.yaml`'s
  comments). This plugin takes `queue`/`account`/extra qsub args as
  profile-level settings only (`workflows/processing/profiles/pbs/config.yaml`),
  not per-rule — so unlike SLURM, they can't vary by a system's `size_class`.
  It also has no native GPU resource; you'll likely need a
  `pbs-extra-qsub-args` GPU request in the profile, which is very
  site-specific and left as a commented placeholder.
- **HTCondor** — via the official `snakemake-executor-plugin-htcondor`
  (`envs/controller.yaml` installs it unconditionally — safe, official
  package). Assumes a shared filesystem between submit and execute nodes by
  default; opportunistic pools without one (OSG/CHTC-style) need
  `shared-fs-usage` tuned and file transfer configured — see the comments
  in `workflows/processing/profiles/htcondor/config.yaml`. `hpc.account`
  isn't wired in here (most HTCondor pools don't use SLURM/PBS-style
  accounting).

If you hit something wrong in either path, it's very likely a resource-key
or profile-key naming mismatch against a newer/older version of that
plugin — check the plugin's own README against what's in the relevant
`config.yaml` comment block and `workflows/processing/profiles/<scheduler>/`.

## Quick start

```bash
# 0. controller env (once)
conda env create -f envs/controller.yaml && conda activate ab-initio-pipeline

# 1. edit config.yaml `systems:` + add configs/<system>.yaml   (start from configs/example_toy.yaml)

# 2. STAGE 1 — local
snakemake -s workflows/preprocessing/Snakefile --use-conda --cores 2
#    review data/annotation/<system>/annotation.yaml, curate the `domains:` block
#    into configs/<system>.yaml, set `annotation_reviewed: true`

# 3. STAGE 2 — on an IFB login node, after rsync-ing data/fold_inputs/
#    (different SLURM cluster, or PBS/HTCondor? see "Running on IFB vs.
#    another SLURM cluster" / "Non-SLURM schedulers" above first)
#    `module load snakemake/9.4.0` also exists but its slurm executor plugin
#    (2.6.0) hangs silently at job submission — use an existing project conda
#    env's snakemake (>=9.24, slurm plugin >=2.7) instead:
export PATH=/shared/projects/<your_project>/conda/envs/<env_with_snakemake>/bin:$PATH
export CONDA_PKGS_DIRS="$PWD/.conda_pkgs"   # isolate rule-env builds from a shared cache
snakemake -s workflows/processing/Snakefile --unlock 2>/dev/null || true
snakemake -s workflows/processing/Snakefile --until prime_backends -c1 -p   # once, no GPU
snakemake -s workflows/processing/Snakefile --profile workflows/processing/profiles/ifb

# 4. STAGE 3 — local, after rsync-ing results/abcfold/ + results/metadata/
snakemake -s workflows/postprocessing/Snakefile --use-conda --cores 4
```

---

## Layout

| Path | What |
|---|---|
| `config.yaml` | all shared defaults, every stage |
| `config.local.yaml` | per-machine / per-user overrides (gitignored; copy `.example`) |
| `configs/<system>.yaml` | per-system chain spec + anchors + curated domains + PLIP passes |
| `envs/` | one conda env per rule group, built by `--use-conda` |
| `tools/` | vendored: `MC2` (MoRFchibi2), `AFF`, `ipsae`, `afm_lis`, `pinc` |
| `scripts/` | carried from the DRB2 pipeline + new (`lib_pipeline.py`, `annotate_*`, `compute_interface_metrics.py`, `cluster_poses.py`, plotters) |
| `workflows/<stage>/Snakefile` | the three stage workflows |
| `workflows/processing/profiles/<scheduler>/` | executor-plugin profile — `ifb` (slurm, validated), `pbs`, `htcondor` (both untested) |
| `report/` | per-stage `.rst` captions, `custom.css` (dark theme), `datavzrd/*.datavzrd.yaml` table specs |
| `data/ results/ logs/ reports/` | runtime outputs (gitignored) |

## Worked example

`configs/example_toy.yaml` (a ubiquitin homodimer smoke test) is kept as a
worked reference for the schema in `configs/_schema.md` — use it as the
template when adding a new `configs/<system>.yaml`.

## Status

**All three stages are implemented and validated end-to-end**, including
the interactive datavzrd table bundles (dark-skinned, embedded in each
`reports/<stage>.zip`): stage 1 locally, stage 2 on IFB SLURM (GPU), stage 3
locally (ChimeraX + Docker/PLIP).
