#!/usr/bin/env python3
"""
scripts/plot_backend_yield.py — Stage 2 report figure
====================================================
Grouped bar: how many models each ABCfold backend contributed, per system.
A backend at zero (token cap / OOM / crash) is immediately visible — the
key check that `abcfold`'s overall exit-0 didn't hide a backend failure.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

BACKENDS = ["alphafold3", "boltz", "chai1", "openfold3", "protenix", "rosettafold3"]
COLORS = ["#4c72b0", "#dd8452", "#55a868", "#c44e52", "#8172b3", "#937860"]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--systems", nargs="+", required=True)
    ap.add_argument("--metadata-root", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--formats", default="svg")
    a = ap.parse_args()

    root = Path(a.metadata_root)
    mat = np.zeros((len(a.systems), len(BACKENDS)), dtype=int)
    for i, s in enumerate(a.systems):
        pq = root / s / "model_metadata.parquet"
        if pq.exists():
            vc = pd.read_parquet(pq, columns=["backend"])["backend"].value_counts()
            for j, b in enumerate(BACKENDS):
                mat[i, j] = int(vc.get(b, 0))

    fig, ax = plt.subplots(figsize=(1.6 + 1.3 * len(a.systems), 4),
                           constrained_layout=True)
    x = np.arange(len(a.systems))
    w = 0.8 / len(BACKENDS)
    for j, b in enumerate(BACKENDS):
        ax.bar(x + (j - len(BACKENDS) / 2) * w + w / 2, mat[:, j], w,
               label=b, color=COLORS[j])
    ax.set_xticks(x)
    ax.set_xticklabels(a.systems, rotation=20, ha="right")
    ax.set_ylabel("models produced")
    ax.set_title("ABCfold per-backend model yield")
    ax.legend(fontsize=7, ncol=3, frameon=False)
    for i in range(len(a.systems)):
        for j in range(len(BACKENDS)):
            if mat[i, j] == 0:
                ax.text(x[i] + (j - len(BACKENDS) / 2) * w + w / 2, 0.5, "0",
                        ha="center", va="bottom", fontsize=6, color="#c44e52")

    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    stem = out.with_suffix("")
    for fmt in a.formats.split(","):
        fig.savefig(f"{stem}.{fmt.strip()}")
    plt.close(fig)
    print(f"[plot_backend_yield] -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
