#!/usr/bin/env python3
"""
scripts/plot_dimreduction.py — Stage 3 report figure
====================================================
Dimensionality-reduction panels of the anchor-aligned partner-Ca pose
vectors (results/<system>/aligned_partner_ca.npy). A configurable set of
methods {pca, umap, tsne, mds}, each rendered twice:

  row 1 — coloured by ABCfold backend, marker size by a rescoring metric
          (iLIS by default, from interface_metrics.parquet)
  row 2 — coloured by pose cluster (pose_clusters.csv)

One SVG grid.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402


def embed(X, method, seed=42):
    if method == "pca":
        from sklearn.decomposition import PCA
        return PCA(n_components=2, random_state=seed).fit_transform(X)
    if method == "mds":
        from sklearn.manifold import MDS
        return MDS(n_components=2, random_state=seed, normalized_stress="auto").fit_transform(X)
    if method == "tsne":
        from sklearn.manifold import TSNE
        p = max(5, min(30, X.shape[0] // 4))
        return TSNE(n_components=2, random_state=seed, perplexity=p).fit_transform(X)
    if method == "umap":
        import umap
        return umap.UMAP(n_components=2, random_state=seed).fit_transform(X)
    raise ValueError(method)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--system", required=True)
    ap.add_argument("--results-root", default="results")
    ap.add_argument("--methods", default="pca,umap,tsne,mds")
    ap.add_argument("--metric", default="ilis")
    ap.add_argument("--out", required=True)
    ap.add_argument("--formats", default="svg")
    a = ap.parse_args()

    base = Path(a.results_root) / a.system
    X = np.load(base / "aligned_partner_ca.npy").reshape(
        np.load(base / "aligned_partner_ca.npy").shape[0], -1)
    pc = pd.read_csv(base / "pose_clusters.csv")
    n = len(pc)
    X = X[:n]

    metric_vals = np.full(n, np.nan)
    mpath = base / "interface_metrics.parquet"
    if mpath.exists() and {"backend", "seed", "sample_index"}.issubset(pc.columns):
        im = pd.read_parquet(mpath)
        agg = (im.groupby(["backend", "seed", "sample_index"])[a.metric]
                 .mean().reset_index())
        key = pc.merge(agg, on=["backend", "seed", "sample_index"], how="left")
        metric_vals = key[a.metric].to_numpy()
    sizes = 20 + 120 * np.nan_to_num(
        (metric_vals - np.nanmin(metric_vals)) /
        (np.nanmax(metric_vals) - np.nanmin(metric_vals) + 1e-9))

    methods = [m.strip() for m in a.methods.split(",") if m.strip()]
    fig, axes = plt.subplots(2, len(methods),
                             figsize=(3.6 * len(methods), 7), squeeze=False,
                             constrained_layout=True)
    backends = sorted(pc["backend"].unique())
    bcol = {b: plt.get_cmap("tab10")(i) for i, b in enumerate(backends)}
    clusters = sorted(pc["cluster"].unique())
    ccol = {c: plt.get_cmap("tab20")(i) for i, c in enumerate(clusters)}

    for j, m in enumerate(methods):
        try:
            E = embed(X, m)
        except Exception as e:                             # noqa: BLE001
            for r in (0, 1):
                axes[r][j].text(0.5, 0.5, f"{m}: {e}", fontsize=7,
                                ha="center", transform=axes[r][j].transAxes)
            continue
        for b in backends:
            s = (pc["backend"] == b).to_numpy()
            axes[0][j].scatter(E[s, 0], E[s, 1], s=sizes[s], color=bcol[b],
                               alpha=.8, edgecolor="white", lw=.3, label=b)
        axes[0][j].set_title(f"{m} — backend / size∝{a.metric}", fontsize=9)
        for c in clusters:
            s = (pc["cluster"] == c).to_numpy()
            axes[1][j].scatter(E[s, 0], E[s, 1], s=26, color=ccol[c],
                               alpha=.85, edgecolor="white", lw=.3, label=f"cl{c}")
        axes[1][j].set_title(f"{m} — pose cluster", fontsize=9)
        for r in (0, 1):
            axes[r][j].set_xticks([]); axes[r][j].set_yticks([])
    axes[0][0].legend(fontsize=6, frameon=False, loc="best")
    axes[1][0].legend(fontsize=6, frameon=False, loc="best")
    fig.suptitle(f"{a.system}: pose-vector dimensionality reduction")

    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    stem = out.with_suffix("")
    for fmt in a.formats.split(","):
        fig.savefig(f"{stem}.{fmt.strip()}")
    plt.close(fig)
    print(f"[plot_dimreduction] -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
