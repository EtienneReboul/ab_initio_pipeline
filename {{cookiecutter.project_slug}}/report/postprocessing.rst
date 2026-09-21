Stage 3 — Postprocessing
========================

The ABCfold ensemble (all backends × seeds × samples), pooled per system,
aligned on its rigid anchor and analysed.

* **Rigid core / RMSF** — per-residue RMSF of the anchor chain(s) after
  iterative rigid-core superposition, and the partner-chain mobility trim
  that defines the pose feature vector.
* **Dimensionality reduction** — PCA / UMAP / t-SNE / MDS of the anchor-aligned
  partner-Cα pose vectors, coloured by ABCfold backend (marker size ∝ a
  rescoring metric) and by pose cluster.
* **Interface metrics** — ipSAE, iLIS and Pinc computed for every model in the
  ensemble from its PAE matrix (``interface_metrics.parquet``).
* **Minimization energy** — ChimeraX energy traces per selected model and a
  per-backend failure-rate table (diverged / NaN).
* **PLIP heatmaps** — domain × domain (or domain × nucleotide) interaction
  counts, one panel set per PLIP pass, stratified total / per backend /
  per pose cluster.
