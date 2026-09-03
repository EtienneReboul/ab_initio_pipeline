#!/usr/bin/env python3
"""
scripts/plot_plip_heatmaps.py — Stage 3 report figure
====================================================
Contact-count heatmaps from a PLIP aggregate CSV
(all_selected_summary_<pass>.csv: columns replica, model, resnr, restype,
reschain, resnr_lig, restype_lig, reschain_lig, dist, interaction_type).

Axes = domain names from configs/<system>.yaml `domains:` for protein
chains, or `nt<resnr>` buckets for RNA/DNA chains. Cell = number of PLIP
interactions between that receptor-domain and ligand-domain.

Strata (configurable): total | backend | cluster  — one heatmap each.
`backend`/`cluster` are parsed from the `model` string (rank_NN_<backend>_
seed<S>_sample<M>) and, for `cluster`, joined to pose_clusters.csv.
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
import yaml  # noqa: E402

MODEL_RE = re.compile(r"_(alphafold3|boltz|chai1|openfold3|protenix|rosettafold3)_"
                      r"seed(\d+)_sample([0-9.]+)")


def domain_map(spec, chain, nt_bucket):
    segs = (spec.get("domains") or {}).get(chain)
    if segs:
        def _lookup(resnr):
            for d in segs:
                if int(d["start"]) <= resnr <= int(d["end"]) and d["kind"] != "morf":
                    return f"{chain}:{d['name']}"
            return f"{chain}:?"
        return _lookup
    seqtype = next((s["type"] for s in spec["sequences"] if s["id"] == chain), "protein")
    if seqtype in ("rna", "dna"):
        return lambda r: f"{chain}:nt{((r - 1) // nt_bucket) * nt_bucket + 1}"
    return lambda r: f"{chain}:{((r - 1) // 50) * 50 + 1}"


def one_heatmap(ax, df, rmap, lmap, title):
    if df.empty:
        ax.text(0.5, 0.5, "no contacts", ha="center", transform=ax.transAxes)
        ax.set_title(title, fontsize=8)
        return
    df = df.assign(rdom=df.apply(lambda x: rmap.get(x["reschain"], lambda r: x["reschain"])(int(x["resnr"])), axis=1),
                   ldom=df.apply(lambda x: lmap.get(x["reschain_lig"], lambda r: x["reschain_lig"])(int(x["resnr_lig"])), axis=1))
    piv = df.pivot_table(index="rdom", columns="ldom", values="dist",
                         aggfunc="count", fill_value=0)
    im = ax.imshow(piv.values, cmap="magma", aspect="auto")
    ax.set_xticks(range(len(piv.columns))); ax.set_xticklabels(piv.columns, rotation=60, ha="right", fontsize=6)
    ax.set_yticks(range(len(piv.index))); ax.set_yticklabels(piv.index, fontsize=6)
    ax.set_title(title, fontsize=8)
    for (i, j), v in np.ndenumerate(piv.values):
        if v:
            ax.text(j, i, int(v), ha="center", va="center", fontsize=5,
                    color="white" if v < piv.values.max() * 0.6 else "black")
    return im


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--system", required=True)
    ap.add_argument("--spec", required=True)
    ap.add_argument("--summary", required=True, help="all_selected_summary_<pass>.csv")
    ap.add_argument("--pose-clusters", default="")
    ap.add_argument("--strata", default="total,backend,cluster")
    ap.add_argument("--nt-bucket", type=int, default=1)
    ap.add_argument("--out", required=True)
    ap.add_argument("--formats", default="svg")
    a = ap.parse_args()

    spec = yaml.safe_load(Path(a.spec).read_text())
    df = pd.read_csv(a.summary)
    if df.empty:
        Path(a.out).parent.mkdir(parents=True, exist_ok=True)
        fig, ax = plt.subplots(figsize=(4, 2))
        ax.text(0.5, 0.5, f"{a.system}: PLIP summary empty", ha="center")
        for fmt in a.formats.split(","):
            fig.savefig(f"{Path(a.out).with_suffix('')}.{fmt.strip()}")
        return 0

    m = df["model"].str.extract(MODEL_RE)
    df["backend"] = m[0]
    df["seed"] = pd.to_numeric(m[1], errors="coerce")
    df["sample_index"] = pd.to_numeric(m[2], errors="coerce")
    if a.pose_clusters and Path(a.pose_clusters).exists():
        pc = pd.read_csv(a.pose_clusters)[["backend", "seed", "sample_index", "cluster"]]
        df = df.merge(pc, on=["backend", "seed", "sample_index"], how="left")
    else:
        df["cluster"] = np.nan

    chains = sorted(set(df["reschain"]) | set(df["reschain_lig"]))
    rmap = {c: domain_map(spec, c, a.nt_bucket) for c in chains}
    lmap = rmap

    strata = [s.strip() for s in a.strata.split(",") if s.strip()]
    panels = []
    for st in strata:
        if st == "total":
            panels.append(("total", df))
        elif st == "backend":
            panels += [(f"backend={b}", g) for b, g in df.groupby("backend")]
        elif st == "cluster":
            panels += [(f"cluster={int(c)}", g) for c, g in df.dropna(subset=["cluster"]).groupby("cluster")]
    if not panels:
        panels = [("total", df)]

    ncol = min(4, len(panels))
    nrow = int(np.ceil(len(panels) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.2 * ncol, 3.6 * nrow),
                             squeeze=False, constrained_layout=True)
    for k, (title, g) in enumerate(panels):
        one_heatmap(axes[k // ncol][k % ncol], g, rmap, lmap, title)
    for k in range(len(panels), nrow * ncol):
        axes[k // ncol][k % ncol].axis("off")
    fig.suptitle(f"{a.system}: PLIP domain×domain contact counts ({Path(a.summary).stem})")

    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    for fmt in a.formats.split(","):
        fig.savefig(f"{out.with_suffix('')}.{fmt.strip()}")
    plt.close(fig)
    print(f"[plot_plip_heatmaps] -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
