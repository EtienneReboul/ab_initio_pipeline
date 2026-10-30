# `configs/<system>.yaml` — per-system schema

One file per entry in `config.yaml`'s `systems:` list. This is the **only**
file (besides one line in `config.yaml`) you edit to add a new dataset.

```yaml
name: my_system            # must equal the filename stem and the systems: entry
size_class: small          # key into config.yaml hpc.<scheduler>.size_classes {small,large}

# ── Chains ────────────────────────────────────────────────────────────────
# Order defines the fold_input.json chain order. ids must be unique A,B,C…
sequences:
  - {type: protein, id: A, name: MyReceptor, uniprot: P12345, sequence: "MABC…"}
  - {type: protein, id: B, name: MyPartner,  uniprot: Q67890, region: [1, 120], sequence: "MDEF…"}
  - {type: rna,     id: C, name: guide,      sequence: "ACGU…"}   # optional
  - {type: dna,     id: D, name: target,     sequence: "ACGT…"}   # optional
  - {type: ligand,  id: E, name: bridge,     ccd: GLY}            # optional
  # uniprot (optional, protein): the UniProt entry this chain IS — stage 1 then
  # reads InterPro's precomputed domains for it (fast, no email, carries the
  # representative-domain flags). `region: [start, end]` (1-based, inclusive)
  # when the chain is a slice of that entry; hits are projected onto it. Left
  # out, it is found by locating the chain's sequence in the entry (when it
  # occurs once — polyubiquitin repeats, for example, need it spelled out).
  # Without it (or if the sequence doesn't match) InterProScan runs on the
  # sequence itself (needs annotate.interproscan.email).
  # type: ligand takes exactly one of `ccd` (code or list of codes) or `smiles`
  # (passed through to fold_input.json); ligands are ignored by annotation and
  # pose clustering.

# ── Covalent bonds (optional) ─────────────────────────────────────────────
# -> fold_input.json `bondedAtomPairs`. [[chain, residue, atom], [chain,
# residue, atom]], 1-based residue index within the chain, CCD atom names.
# AF3 rejects polymer-polymer bonds, so each bond must touch a ligand chain
# (bridge e.g. a Ub-Ub isopeptide through the donor's C-terminal GLY as a
# ligand). Boltz/Chai-1/Protenix get them via ABCfold's converters;
# OpenFold3 and RosettaFold3 silently drop them.
bonds:
  - [[A, 75, C], [E, 1, N]]
  - [[E, 1, C], [B, 63, NZ]]

# ── Rigid structural anchors for pose clustering ─────────────────────────
# A LIST. Each anchor gets its own iterative rigid-core refinement and its
# own rmsf_profile figure. Partner chains are what the pose feature vector
# is built from (protein Ca only — RNA/DNA ride along but don't drive
# cluster assignment).
anchor_chains:  [A]
partner_chains: [B]

# ── Domain / region map (OPTIONAL override) ────────────────────────────
# Stage 1 writes the domain map stage 3 uses to data/annotation/<system>/
# annotation.yaml: InterPro's representative domains, named after their
# InterPro entries, the gaps between them split into disordered / linker by
# AIUPred, plus AIUPred MoRFs. No review step. A chain listed here replaces
# that chain's map (others keep theirs). Per chain, an ordered list of
# non-overlapping segments; `kind` is one of: domain | linker | disordered |
# morf. 1-based, inclusive. Used by the stage-3 PLIP domain×domain heatmaps.
domains:
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
