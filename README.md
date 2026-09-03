# ab_initio_pipeline

A generalised, configure-and-run template for **ab-initio complex modelling
with [ABCfold](https://github.com/rigdenlab/ABCFold)** on the IFB Core
Cluster. Distilled from three one-off sibling pipelines (NPF transporters;
DRB2/DCL4/dsRNA; AGO1/FBW2/miR) that all ran the same 3-stage protocol and
diverged only in hand-edited scripts.

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
snakemake -s workflows/processing/Snakefile \
    --profile workflows/processing/profiles/ifb --until prime_backends   # once
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
| `workflows/processing/profiles/ifb/` | SLURM executor-plugin profile |
| `report/` | per-stage `.rst` captions, `custom.css` (dark theme), `datavzrd/*.datavzrd.yaml` table specs |
| `data/ results/ logs/ reports/` | runtime outputs (gitignored) |

## Carried real examples

`configs/drb2_drb4.yaml`, `rna_ds_drb2_drb4.yaml`, `rna_ds_dcl4_drb2_drb4.yaml`
and `ago1_fbw2*.yaml` are the real systems from the source pipelines, kept as
worked references. `lib_pipeline.load_system()` upgrades their pre-schema
shape on the fly (`anchor_chain` → `anchor_chains: [..]`, the legacy
`plip:` / `plip_rna_ligands:` / `plip_drb2_drb4:` triple → `plip_passes:`).

## Status

**Stage 1 is implemented and DAG-validated**, including the interactive
datavzrd table bundles (`report/datavzrd/preprocessing.datavzrd.yaml`,
embedded dark-skinned in `reports/preprocessing.zip`). Stages 2–3 and their
own datavzrd bundles are scaffolded per [`../ab_initio_pipeline` plan]; see
`git log` / the plan file for what remains.
