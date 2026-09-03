#!/usr/bin/env python3
"""
scripts/plot_msa_coverage.py — Stage 1 report figure
====================================================
ColabFold-style MSA coverage plot for one protein chain, read straight from
the embedded `unpairedMsa` a3m in fold_input.resolved.json.

Top: per-position coverage (number of non-gap sequences). Background image:
the alignment sorted by identity to the query, coloured non-gap / gap, so
depth and its positional unevenness are both visible. SVG (+ any extra
--formats).
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib_pipeline import norm_id  # noqa: E402


def parse_a3m(a3m: str) -> tuple[list[str], np.ndarray]:
    names, seqs, cur = [], [], []
    for ln in a3m.splitlines():
        if ln.startswith(">"):
            if cur:
                seqs.append("".join(cur)); cur = []
            names.append(ln[1:].split()[0])
        elif ln.strip():
            cur.append(ln.strip())
    if cur:
        seqs.append("".join(cur))
    # drop lowercase (insertions relative to query) -> match-state columns
    clean = ["".join(c for c in s if not c.islower()) for s in seqs]
    L = len(clean[0])
    clean = [s for s in clean if len(s) == L]
    arr = np.array([[0 if c in "-." else 1 for c in s] for s in clean], dtype=np.int8)
    return names, arr


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--resolved", required=True)
    ap.add_argument("--chain", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--formats", default="svg")
    args = ap.parse_args()

    doc = json.loads(Path(args.resolved).read_text())
    a3m = None
    for entry in doc.get("sequences", []):
        p = entry.get("protein")
        if p and norm_id(p.get("id")) == args.chain:
            a3m = p.get("unpairedMsa") or ""
            break
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    fig, (ax0, ax1) = plt.subplots(
        2, 1, figsize=(9, 4.5), height_ratios=[1, 2.4], sharex=True,
        constrained_layout=True)

    if not a3m:
        for ax in (ax0, ax1):
            ax.text(0.5, 0.5, f"chain {args.chain}: no MSA embedded",
                    ha="center", va="center", transform=ax.transAxes)
            ax.set_xticks([]); ax.set_yticks([])
    else:
        names, arr = parse_a3m(a3m)
        cov = arr.sum(axis=0)
        L = arr.shape[1]
        order = np.argsort(-arr.mean(axis=1))          # densest rows first
        ax0.fill_between(np.arange(1, L + 1), cov, step="mid", color="#4c72b0", alpha=.9)
        ax0.set_ylabel("seqs"); ax0.set_title(
            f"chain {args.chain} — {arr.shape[0]} sequences, "
            f"median coverage {int(np.median(cov))}")
        ax1.imshow(arr[order], aspect="auto", interpolation="nearest",
                   cmap="Greys", vmin=0, vmax=1,
                   extent=[0.5, L + 0.5, arr.shape[0], 0])
        ax1.set_ylabel("sequence (by identity)")
        ax1.set_xlabel("query residue")

    stem = out.with_suffix("")
    for fmt in args.formats.split(","):
        fig.savefig(f"{stem}.{fmt.strip()}", bbox_inches="tight")
    plt.close(fig)
    print(f"[plot_msa_coverage] -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
