Stage 1 — Preprocessing
=======================

Inputs prepared for ABCfold, plus sequence annotation, for every system in
``config.yaml``'s ``systems:`` list.

* **MSA coverage** — per protein chain, the ColabFold-webserver MSA embedded
  in ``fold_input.resolved.json``: depth and its positional unevenness.
* **Template hits** — PDB templates MMseqs2 attached to each chain.
* **Domain map** — curated (``annotation_reviewed: true``) or auto-proposed
  segment map per chain: domain / linker / disordered / MoRF.
* **Disorder / MoRF** — superposed IUPred3 (long/short), ANCHOR2, AIUPred and
  MoRFchibi2 tracks, consensus-disordered regions shaded.
* **Input stats** — chains, residues, tokens, MSA depth, and per-backend
  token-cap risk flags.

The auto-annotation at ``data/annotation/<system>/annotation.yaml`` is a
**proposal for human review** — curate it into ``configs/<system>.yaml`` and
set ``annotation_reviewed: true`` before running Stage 2.
