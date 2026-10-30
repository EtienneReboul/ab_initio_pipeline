Stage 3 — Postprocessing
========================

The ABCfold ensemble (all backends × seeds × samples), pooled per system,
aligned on its rigid anchor and analysed.

* **Rigid core / RMSF** — per-residue RMSF of the anchor chain(s) after
  iterative rigid-core superposition, and the partner-chain mobility trim
  that defines the pose feature vector.
* **Dimensionality reduction** — PCA / UMAP / t-SNE / MDS of the anchor-aligned
  partner-Cα pose vectors: marker shape = ABCfold backend, marker color = a
  rescoring metric (shared colorbar), and a second panel coloured by pose
  cluster.
* **Interface metrics** — ipSAE, iLIS and Pinc computed for every model in the
  ensemble from its PAE matrix (``interface_metrics.parquet``).
* **Minimization energy** — ChimeraX energy traces per selected model and a
  per-backend failure-rate table (diverged / NaN).
* **PLIP heatmaps** — domain × domain (or domain × nucleotide) interaction
  counts, one panel set per PLIP pass, stratified total / per backend /
  per pose cluster.

Score definitions and references
--------------------------------

*Global fold confidence* — the predicted complex as a whole (high values mean
a well-predicted structure, not necessarily an interaction):

* **pTM** — predicted TM-score of the whole complex (0-1). Jumper et al. 2021,
  https://doi.org/10.1038/s41586-021-03819-2; Zhang & Skolnick 2004,
  https://doi.org/10.1002/prot.20264
* **mean pLDDT** — mean per-residue predicted lDDT (0-1). Jumper et al. 2021;
  Mariani et al. 2013, https://doi.org/10.1093/bioinformatics/btt473
* **ranking score** — each backend's own ranking score (AF3: 0.8 ipTM +
  0.2 pTM + disorder/clash terms). Abramson et al. 2024,
  https://doi.org/10.1038/s41586-024-07487-w

*Interface* — the anchor-partner interface only:

* **ipTM** — interface pTM over inter-chain residue pairs (0-1). Evans et al.
  2021, https://doi.org/10.1101/2021.10.04.463034
* **ipSAE** — ipTM-like score restricted to residue pairs with inter-chain PAE
  under the cutoff, max of the two directions (0-1). Dunbrack 2025,
  https://doi.org/10.1101/2025.02.10.637595
* **LIS** — Local Interaction Score: mean (cutoff - PAE) / cutoff over
  inter-chain pairs with PAE <= 12 Å (0-1). Kim et al. 2024,
  https://doi.org/10.1101/2024.02.19.580970
* **cLIS** — LIS over the pairs also within 8 Å (Cβ-Cβ); **iLIS** =
  sqrt(LIS × cLIS) (0-1). Kim et al. 2026, https://doi.org/10.64898/2026.04.14.718529
* **Pinc** — mean PAE-derived contact probability over inter-chain pairs within
  12 Å (0-1). Badonyi & Toth-Petroczy 2026, https://doi.org/10.1002/pro.70760
* **pDockQ** — predicted DockQ from interface pLDDT and contact count (0-1).
  Bryant et al. 2022, https://doi.org/10.1038/s41467-022-28865-w
* **pDockQ2** — predicted DockQ from interface pLDDT and PAE (0-1). Zhu et al.
  2023, https://doi.org/10.1093/bioinformatics/btad424
* **PLIP contacts** — non-covalent interactions detected in the minimized
  model. Salentin et al. 2015, https://doi.org/10.1093/nar/gkv315; Adasme et
  al. 2021, https://doi.org/10.1093/nar/gkab294

(References as compiled, and DOIs resolved against Crossref, in NanoPulldown's
``nanopulldown.core.metric_refs``.)
