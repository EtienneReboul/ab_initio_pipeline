# `configs/<system>.yaml` — per-system schema

One file per entry in `config.yaml`'s `systems:` list. This is the **only**
file (besides one line in `config.yaml`) you edit to add a new dataset.

```yaml
name: my_system            # must equal the filename stem and the systems: entry
size_class: small          # key into config.yaml hpc.<scheduler>.size_classes {small,large}

# ── Chains ────────────────────────────────────────────────────────────────
# Order defines the fold_input.json chain order. ids must be unique A,B,C…
sequences:
  - {type: protein, id: A, name: MyReceptor, sequence: "MABC…"}
  - {type: protein, id: B, name: MyPartner,  sequence: "MDEF…"}
  - {type: rna,     id: C, name: guide,      sequence: "ACGU…"}   # optional
  - {type: dna,     id: D, name: target,     sequence: "ACGT…"}   # optional
  # type: ligand with a `ccd` or `smiles` key is also accepted (passed through
  # to fold_input.json); ligands are ignored by annotation and pose clustering.

# ── Rigid structural anchors for pose clustering ─────────────────────────
# A LIST. Each anchor gets its own iterative rigid-core refinement and its
# own rmsf_profile figure. Partner chains are what the pose feature vector
# is built from (protein Ca only — RNA/DNA ride along but don't drive
# cluster assignment).
anchor_chains:  [A]
partner_chains: [B]

# ── Annotation review gate ──────────────────────────────────────────────
# Stage 1 writes data/annotation/<system>/annotation.yaml for you to merge
# into the `domains:` block below. Flip this to true once you've reviewed
# it — stage 2 refuses to run for a system still set false.
annotation_reviewed: false

# ── Domain / region map (hand-curated from stage-1 annotation) ──────────
# Per chain id, an ordered list of non-overlapping segments. `kind` is one
# of: domain | linker | disordered | morf. 1-based, inclusive. Used by the
# stage-3 PLIP domain×domain heatmaps and the energy/over-fold filter.
domains:
  A:
    - {name: dom1,  start: 1,   end: 76,  kind: domain}
  B:
    - {name: dom1,  start: 1,   end: 60,  kind: domain}
    - {name: tail,  start: 61,  end: 95,  kind: disordered}
    - {name: morf1, start: 78,  end: 88,  kind: morf}

# ── PLIP passes ────────────────────────────────────────────────────────
# Generalises the DRB2 pipeline's hard-coded plip / plip_rna_ligands /
# plip_drb2_drb4 triple. One `run_plip` + aggregate per entry.
#   chains: PLIP --chains string "[['A'],['B','C']]" (receptor groups first),
#           or null to omit --chains (auto-detect every non-receptor as its
#           own ligand — needed for per-nucleotide RNA contacts).
#   dnareceptor: pass PLIP --dnareceptor (folds NA into the receptor).
#   fix_pdb: run pdb4amber + PDBFixer first (needed whenever NA is present,
#            and for ABCfold RNA atom naming even in the RNA-ligand pass).
plip_passes:
  - {name: main, dnareceptor: false, fix_pdb: false, chains: "[['A'],['B']]"}

plip_host:
  docker_memory: "8g"
  docker_platform: "linux/amd64"   # "" for native linux
```

## Worked example

`example_toy.yaml` (a ubiquitin homodimer smoke test) is included as a
worked reference for this schema — use it as the template when adding a
new `configs/<system>.yaml`.
