#!/usr/bin/env python3
"""
scripts/plot_disorder.py — Stage 1 report figure
================================================
Superposed per-residue disorder / MoRF curves for one protein chain
(IUPred3 long/short, ANCHOR2, AIUPred, MoRFchibi2), with residues where a
majority of tracks exceed --cutoff shaded.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

TRACK_STYLE = {
    "metapredict":   ("#1f77b4", "-"),
    "iupred3_long":  ("#4c72b0", "-"),
    "iupred3_short": ("#4c72b0", "--"),
    "anchor2":       ("#8172b3", "-"),
    "aiupred":       ("#c44e52", "-"),
    "morfchibi":     ("#55a868", "-"),
}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--disorder", required=True)
    ap.add_argument("--morf", required=True)
    ap.add_argument("--chain", required=True)
    ap.add_argument("--cutoff", type=float, default=0.5)
    ap.add_argument("--out", required=True)
    ap.add_argument("--formats", default="svg")
    args = ap.parse_args()

    dis = pd.read_csv(args.disorder, sep="\t")
    dis = dis[dis["chain"] == args.chain]
    morf = pd.read_csv(args.morf, sep="\t")
    morf = morf[morf["chain"] == args.chain].assign(track="morfchibi").rename(
        columns={"score": "score"})[["chain", "resi", "track", "score"]]
    allrows = pd.concat([dis[["chain", "resi", "track", "score"]], morf], ignore_index=True)

    fig, ax = plt.subplots(figsize=(9, 3.2), constrained_layout=True)
    if allrows.empty:
        ax.text(0.5, 0.5, f"chain {args.chain}: no disorder tracks",
                ha="center", va="center", transform=ax.transAxes)
    else:
        wide = allrows.pivot_table(index="resi", columns="track", values="score")
        wide = wide.sort_index()
        for t in wide.columns:
            c, ls = TRACK_STYLE.get(t, ("#333333", "-"))
            ax.plot(wide.index, wide[t], color=c, ls=ls, lw=1.2, label=t)
        dcols = [c for c in wide.columns if c != "morfchibi"]
        if dcols:
            frac = (wide[dcols] >= args.cutoff).sum(axis=1) / len(dcols)
            strong = frac >= 0.5
            ax.fill_between(wide.index, 0, 1, where=strong, color="#dd8452",
                            alpha=0.15, step="mid", label="consensus disordered")
        ax.axhline(args.cutoff, color="grey", lw=.7, ls=":")
        ax.set_ylim(0, 1)
        ax.legend(loc="upper right", fontsize=7, ncol=2, frameon=False)
    ax.set_xlabel("residue"); ax.set_ylabel("score")
    ax.set_title(f"chain {args.chain} — disorder / MoRF", fontsize=9)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    stem = out.with_suffix("")
    for fmt in args.formats.split(","):
        fig.savefig(f"{stem}.{fmt.strip()}")
    plt.close(fig)
    print(f"[plot_disorder] -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
