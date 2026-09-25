#!/usr/bin/env python3
"""
scripts/plot_dimreduction.py — Stage 3 report figure
====================================================
Dimensionality-reduction panels of the anchor-aligned partner-Ca pose
vectors (results/<system>/aligned_partner_ca.npy). A configurable set of
methods {pca, umap, tsne, mds}, each rendered twice:

  row 1 — marker shape = ABCfold backend, marker color = --metric (one of
          DIMR_METRICS in workflows/postprocessing/Snakefile), with a
          shared colorbar bounded to that metric's known theoretical range
          where one exists (falls back to this run's data range otherwise)
  row 2 — coloured by pose cluster (pose_clusters.csv)

One SVG per metric (workflows/postprocessing/Snakefile's fig_dimreduction
expands over DIMR_METRICS), so the report browses them like the
msa_coverage_{chain} tabs — same embeddings (fixed random_state), different
coloring, side by side for comparing what each metric actually shows.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

# chain-pair rows in interface_metrics.parquet, averaged per model below
INTERFACE_METRICS = {"ipsae", "ilis", "pinc", "pdockq", "pdockq2", "lis", "clis"}
# one row per model in model_metadata.parquet, no averaging needed
MODEL_METRICS = {"ptm", "iptm", "ranking_score", "mean_plddt"}

# Theoretical (not just this run's observed) min/max, so the colorbar reads
# as an absolute confidence scale rather than stretching to whatever range
# this one system happened to produce. Omit a metric here to fall back to
# its data-driven min/max (labeled as such on the colorbar).
#   ipsae, ilis, pinc, ptm, iptm: each defined as a probability / PAE-derived
#     score in [0, 1] by construction (tools/ipsae/ipsae.py, tools/afm_lis,
#     tools/pinc/README.md).
#   pdockq: 0.724/(1+exp(-0.052*(x-152.611)))+0.018, x=mean_plddt*log10(npairs)
#     >= 0; ceiling is the sigmoid's asymptote as x->inf (Bryant et al. 2022).
#   pdockq2: 1.31/(1+exp(-0.075*(x-84.733)))+0.005, x=mean_plddt*mean_ptm.
#     The sigmoid's own asymptote is 1.315, but x is bounded above by
#     100*1=100 (mean_plddt<=100, mean_ptm<=1), and f(100)=0.999 — so 1.0 is
#     the practically-reachable ceiling given those bounded inputs.
#   mean_plddt: standard 0-100 pLDDT scale.
#   ranking_score: deliberately absent — backend-defined (AF3-style:
#     weighted ptm/iptm minus a clash penalty), observed as low as -99.7 in
#     this project's own data, so there is no fixed theoretical range to
#     bound it to.
METRIC_BOUNDS = {
    "ipsae": (0.0, 1.0),
    "ilis": (0.0, 1.0),
    "pinc": (0.0, 1.0),
    "pdockq": (0.0, 0.742),
    "pdockq2": (0.0, 1.0),
    "ptm": (0.0, 1.0),
    "iptm": (0.0, 1.0),
    "mean_plddt": (0.0, 100.0),
}


def embed(X, method, seed=42):
    if method == "pca":
        from sklearn.decomposition import PCA
        p = PCA(n_components=2, random_state=seed)
        E = p.fit_transform(X)
        return E, p.explained_variance_ratio_
    if method == "mds":
        from sklearn.manifold import MDS
        return MDS(n_components=2, random_state=seed, normalized_stress="auto").fit_transform(X), None
    if method == "tsne":
        from sklearn.manifold import TSNE
        p = max(5, min(30, X.shape[0] // 4))
        return TSNE(n_components=2, random_state=seed, perplexity=p).fit_transform(X), None
    if method == "umap":
        import umap
        return umap.UMAP(n_components=2, random_state=seed).fit_transform(X), None
    raise ValueError(method)


def load_metric(metric, base, metadata_root, system, pc):
    """Per-model metric values aligned to pose_clusters.csv row order."""
    n = len(pc)
    if metric in pc.columns:
        # pose_clusters.csv already carries ptm/iptm/ranking_score straight
        # from align_anchor's own model_metadata read — use it directly
        # rather than re-merging model_metadata.parquet, which would collide
        # on the column name and get pandas-suffixed (ptm_x/ptm_y) instead.
        return pc[metric].to_numpy()
    if metric in INTERFACE_METRICS:
        mpath = base / "interface_metrics.parquet"
        if not mpath.exists() or not {"backend", "seed", "sample_index"}.issubset(pc.columns):
            return np.full(n, np.nan)
        im = pd.read_parquet(mpath)
        agg = (im.groupby(["backend", "seed", "sample_index"])[metric]
                 .mean().reset_index())
        key = pc.merge(agg, on=["backend", "seed", "sample_index"], how="left")
        return key[metric].to_numpy()
    if metric in MODEL_METRICS:
        mpath = Path(metadata_root) / system / "model_metadata.parquet"
        if not mpath.exists() or not {"backend", "seed", "sample_index"}.issubset(pc.columns):
            return np.full(n, np.nan)
        mm = pd.read_parquet(mpath)
        # (backend, seed, sample_index) isn't a unique key for every backend
        # — e.g. OpenFold3 stores sample_index as NaN for all 5 diffusion
        # samples of a seed, so a plain merge fans out (crashes downstream
        # on the row-count mismatch). groupby+mean collapses those instead
        # of fanning out, same as the INTERFACE_METRICS path above; coarser
        # (seed-level-averaged) for whichever backend doesn't have a real
        # per-sample index, rather than wrong-shaped.
        agg = (mm.groupby(["backend", "seed", "sample_index"])[metric]
                 .mean().reset_index())
        key = pc.merge(agg, on=["backend", "seed", "sample_index"], how="left")
        return key[metric].to_numpy()
    raise ValueError(f"unknown metric {metric!r} (not in INTERFACE_METRICS or MODEL_METRICS)")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--system", required=True)
    ap.add_argument("--results-root", default="results")
    ap.add_argument("--metadata-root", default="results/metadata")
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

    metric_vals = load_metric(a.metric, base, a.metadata_root, a.system, pc)

    theoretical = a.metric in METRIC_BOUNDS
    if theoretical:
        vmin, vmax = METRIC_BOUNDS[a.metric]
    elif np.isfinite(metric_vals).any():
        vmin = float(np.nanmin(metric_vals))
        vmax = float(np.nanmax(metric_vals))
        if vmin == vmax:
            vmax = vmin + 1e-9
    else:
        vmin, vmax = 0.0, 1.0

    MARKERS = ["o", "s", "^", "D", "v", "P", "X", "*"]

    methods = [m.strip() for m in a.methods.split(",") if m.strip()]
    fig, axes = plt.subplots(2, len(methods),
                             figsize=(3.6 * len(methods), 7), squeeze=False,
                             constrained_layout=True)
    backends = sorted(pc["backend"].unique())
    bmark = {b: MARKERS[i % len(MARKERS)] for i, b in enumerate(backends)}
    clusters = sorted(pc["cluster"].unique())
    ccol = {c: plt.get_cmap("tab20")(i) for i, c in enumerate(clusters)}
    cmap = plt.get_cmap("viridis")

    last_sc = None
    for j, m in enumerate(methods):
        try:
            E, var_ratio = embed(X, m)
        except Exception as e:                             # noqa: BLE001
            for r in (0, 1):
                axes[r][j].text(0.5, 0.5, f"{m}: {e}", fontsize=7,
                                ha="center", transform=axes[r][j].transAxes)
            continue
        for b in backends:
            s = (pc["backend"] == b).to_numpy()
            last_sc = axes[0][j].scatter(
                E[s, 0], E[s, 1], marker=bmark[b], c=metric_vals[s],
                cmap=cmap, vmin=vmin, vmax=vmax, s=42, alpha=.9,
                edgecolor="white", lw=.3, label=b)
        axes[0][j].set_title(f"{m} — shape=backend, color={a.metric}", fontsize=9)
        for c in clusters:
            s = (pc["cluster"] == c).to_numpy()
            axes[1][j].scatter(E[s, 0], E[s, 1], s=26, color=ccol[c],
                               alpha=.85, edgecolor="white", lw=.3, label=f"cl{c}")
        axes[1][j].set_title(f"{m} — pose cluster", fontsize=9)
        for r in (0, 1):
            axes[r][j].set_xticks([]); axes[r][j].set_yticks([])
        if var_ratio is not None:
            for r in (0, 1):
                axes[r][j].set_xlabel(f"PC1 ({var_ratio[0]*100:.1f}%)", fontsize=7)
                axes[r][j].set_ylabel(f"PC2 ({var_ratio[1]*100:.1f}%)", fontsize=7)

    from matplotlib.lines import Line2D
    shape_handles = [Line2D([0], [0], marker=bmark[b], linestyle="none",
                            markerfacecolor="0.6", markeredgecolor="white",
                            markersize=7, label=b) for b in backends]
    axes[0][0].legend(handles=shape_handles, fontsize=6, frameon=False, loc="best")
    axes[1][0].legend(fontsize=6, frameon=False, loc="best")
    if last_sc is not None:
        cb = fig.colorbar(last_sc, ax=axes[0, :].tolist(), fraction=0.025,
                          pad=0.02)
        cb.ax.tick_params(labelsize=6)
        bound_note = "theoretical range" if theoretical else "data range"
        cb.set_label(f"{a.metric} ({bound_note})", fontsize=7)
    fig.suptitle(f"{a.system}: pose-vector dimensionality reduction — {a.metric}")

    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    stem = out.with_suffix("")
    for fmt in a.formats.split(","):
        fig.savefig(f"{stem}.{fmt.strip()}")
    plt.close(fig)
    print(f"[plot_dimreduction] {a.metric} -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
