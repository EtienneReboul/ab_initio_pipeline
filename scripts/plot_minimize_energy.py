#!/usr/bin/env python3
"""
scripts/plot_minimize_energy.py — Stage 3 report figure + table
=============================================================
Reads every ChimeraX minimize energy trace under
results/<system>/minimized/**/<stem>_energy.csv (step,energy_kJ_mol, written
by scripts/minimize_cif.py) and produces:

  * <out>.svg           seaborn line plot of energy vs step, one line per
                        model, coloured by backend (backend parsed from path)
  * <out_table>.csv     per-backend failure-rate table:
                        backend, n_models, n_diverged, n_nan, frac_failed
    diverged = |final energy| > --max-abs-energy  OR  final energy is NaN
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

try:
    import seaborn as sns
    _HAS_SNS = True
except Exception:                                          # noqa: BLE001
    _HAS_SNS = False

BACKEND_RE = re.compile(r"(alphafold3|boltz|chai1|openfold3|protenix|rosettafold3)")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--system", required=True)
    ap.add_argument("--results-root", default="results")
    ap.add_argument("--max-abs-energy", type=float, default=1e15)
    ap.add_argument("--out", required=True)
    ap.add_argument("--out-table", required=True)
    ap.add_argument("--formats", default="svg")
    a = ap.parse_args()

    mroot = Path(a.results_root) / a.system / "minimized"
    traces = []
    finals = []
    for csv in mroot.rglob("*_energy.csv"):
        bm = BACKEND_RE.search(str(csv))
        backend = bm.group(1) if bm else "unknown"
        try:
            d = pd.read_csv(csv)
        except Exception:                                  # noqa: BLE001
            continue
        ecol = next((c for c in d.columns if "energy" in c.lower()), None)
        scol = next((c for c in d.columns if "step" in c.lower()), None)
        if ecol is None:
            continue
        step = d[scol] if scol else np.arange(len(d))
        e = pd.to_numeric(d[ecol], errors="coerce")
        traces.append(pd.DataFrame({"step": step, "energy": e,
                                    "backend": backend, "model": csv.stem}))
        fe = e.dropna()
        final = fe.iloc[-1] if len(fe) else np.nan
        finals.append({"backend": backend, "model": csv.stem,
                       "final_energy": final,
                       "nan": bool(e.isna().any() or np.isnan(final)),
                       "diverged": bool(np.isfinite(final) and abs(final) > a.max_abs_energy)})

    out = Path(a.out); out.parent.mkdir(parents=True, exist_ok=True)
    tbl = Path(a.out_table); tbl.parent.mkdir(parents=True, exist_ok=True)

    fdf = pd.DataFrame(finals)
    if fdf.empty:
        pd.DataFrame(columns=["backend", "n_models", "n_diverged", "n_nan", "frac_failed"]).to_csv(tbl, index=False)
        fig, ax = plt.subplots(figsize=(5, 3))
        ax.text(0.5, 0.5, f"{a.system}: no *_energy.csv found", ha="center")
        for fmt in a.formats.split(","):
            fig.savefig(f"{out.with_suffix('')}.{fmt.strip()}")
        return 0

    summ = (fdf.assign(failed=fdf["nan"] | fdf["diverged"])
               .groupby("backend")
               .agg(n_models=("model", "count"),
                    n_diverged=("diverged", "sum"),
                    n_nan=("nan", "sum"),
                    frac_failed=("failed", "mean")).reset_index())
    summ["frac_failed"] = summ["frac_failed"].round(3)
    summ.to_csv(tbl, index=False)

    tr = pd.concat(traces, ignore_index=True)
    fig, ax = plt.subplots(figsize=(8, 4.5), constrained_layout=True)
    if _HAS_SNS:
        sns.lineplot(tr, x="step", y="energy", units="model", hue="backend",
                     estimator=None, lw=0.6, alpha=0.5, ax=ax)
    else:
        for (b, mdl), g in tr.groupby(["backend", "model"]):
            ax.plot(g["step"], g["energy"], lw=0.5, alpha=0.5)
    ax.set_yscale("symlog")
    ax.set_xlabel("minimization step"); ax.set_ylabel("energy (kJ/mol, symlog)")
    ax.set_title(f"{a.system}: ChimeraX minimization energy traces")
    for fmt in a.formats.split(","):
        fig.savefig(f"{out.with_suffix('')}.{fmt.strip()}")
    plt.close(fig)
    print(f"[plot_minimize_energy] -> {out} ; table -> {tbl}")
    print(summ.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
