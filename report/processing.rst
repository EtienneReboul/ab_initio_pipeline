Stage 2 — Processing (ABCfold on IFB)
====================================

One SLURM job per system launched all selected ABCfold backends
(AlphaFold3 / Boltz-2 / Chai-1 / OpenFold3 / Protenix / RosettaFold3)
against the same ``fold_input.resolved.json``, then compressed each run's
raw per-sample confidence sprawl to ``arrays.h5`` + ``model_metadata.parquet``.

* **Run timing** — wall-clock per system (from ``sacct`` when available),
  SLURM exit state, and — where a backend produced no usable models —
  the failure reason parsed from its log. ``abcfold`` returns success even
  when individual backends fail, so this table is the real completeness check.
* **Metadata compression** — bytes before/after and the ratio, per system.
* **Per-backend model yield** — how many models each backend actually
  contributed. A backend at zero (token cap, OOM, crash) shows up here.
