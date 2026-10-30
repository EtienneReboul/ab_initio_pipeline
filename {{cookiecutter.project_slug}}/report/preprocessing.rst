Stage 1 — Preprocessing
=======================

Inputs prepared for ABCfold, plus sequence annotation, for every system in
``config.yaml``'s ``systems:`` list.

* **MSA coverage** — per protein chain, the ColabFold-webserver MSA embedded
  in ``fold_input.resolved.json``, drawn like ColabFold's own "Sequence
  coverage" plot: one row per aligned sequence, coloured by identity to the
  query; the line is the per-position depth.
* **Template hits** — PDB templates MMseqs2 attached to each chain.
* **Domain map + disorder** — per protein chain, laid out like InterPro's
  protein viewer: AIUPred's disorder / binding (MoRF) / linker tracks; the
  protein's family; the **domain map** stage 3 uses; one row per InterPro
  entry with the member databases that agree on it; unintegrated signatures;
  MobiDB-lite, coiled coils, signal peptides, TM helices and AIUPred MoRFs.
  The domain map is InterPro's representative domain set (plus InterPro
  domain entries it leaves uncovered), each region named after its InterPro
  entry, with the gaps between domains called disordered or linker from the
  AIUPred disorder track. A chain's ``domains:`` in ``configs/<system>.yaml``
  replaces it (labelled "curated").
* **Input stats** — chains, residues, tokens, MSA depth, and per-backend
  token-cap risk flags.

The domain map at ``data/annotation/<system>/annotation.yaml`` is used as is
by Stage 3 — there is no review step before Stage 2. Give a chain a
``uniprot:`` accession in ``configs/<system>.yaml`` to annotate it from
InterPro's precomputed matches (otherwise InterProScan runs on the sequence).
