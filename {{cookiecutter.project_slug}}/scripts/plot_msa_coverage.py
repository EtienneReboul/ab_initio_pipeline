#!/usr/bin/env python3
"""
scripts/plot_msa_coverage.py — Stage 1 report figure
====================================================
ColabFold-style MSA coverage plot for one protein chain, read straight from
the embedded `unpairedMsa` a3m in fold_input.resolved.json.

Reproduces the AlphaFold2 Colab notebook's "Sequence coverage" figure: each
row is one aligned hit (match-state columns only, sorted by identity to the
query), coloured by that hit's fractional identity to the query wherever it
has a residue, left blank where it has a gap. A black line overlays the
per-position depth (non-gap sequence count) on the same axis. SVG (+ any
extra --formats).
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


def parse_a3m(a3m: str) -> list[str]:
    """Match-state sequences (insertions relative to the query dropped)."""
    seqs, cur = [], []
    for ln in a3m.splitlines():
        if ln.startswith(">"):
            if cur:
                seqs.append("".join(cur)); cur = []
        elif ln.strip():
            cur.append(ln.strip())
    if cur:
        seqs.append("".join(cur))
    clean = ["".join(c for c in s if not c.islower()) for s in seqs]
    L = len(clean[0])
    return [s for s in clean if len(s) == L]


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

    fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    clean = parse_a3m(a3m) if a3m else []
    if len(clean) < 2:
        ax.text(0.5, 0.5, f"chain {args.chain}: no MSA embedded",
                 ha="center", va="center", transform=ax.transAxes)
        ax.set_xticks([]); ax.set_yticks([])
    else:
        query = np.array(list(clean[0]))
        msa = np.array([list(s) for s in clean])          # (N, L)
        L = msa.shape[1]
        gap = msa != "-"
        # Identity to query over the *full* length (as colabfold.plot.plot_msa_v2
        # does), not just the row's own non-gap span — so a short, otherwise
        # perfect partial hit still reads as low identity, matching the
        # ColabFold notebook's colouring.
        match = msa == query[None, :]
        seqid = match.sum(axis=1) / L

        lines = np.where(gap, seqid[:, None], np.nan)
        lines = lines[np.argsort(seqid)]                    # low->high identity
        coverage = gap.sum(axis=0)                          # per-position depth

        cmap = matplotlib.colormaps["rainbow_r"].copy()
        cmap.set_bad(alpha=0)
        im = ax.imshow(np.ma.masked_invalid(lines), interpolation="nearest",
                        aspect="auto", cmap=cmap, vmin=0, vmax=1, origin="lower",
                        extent=(0.5, L + 0.5, 0, lines.shape[0]))
        ax.plot(np.arange(1, L + 1), coverage, color="black", linewidth=1)

        ax.set_xlim(0.5, L + 0.5)
        ax.set_ylim(0, lines.shape[0])
        ax.set_xlabel("positions")
        ax.set_ylabel("sequences")
        ax.set_title(f"chain {args.chain} — {msa.shape[0]} sequences, "
                     f"median depth {int(np.median(coverage))}")
        cbar = fig.colorbar(im, ax=ax, fraction=0.05, pad=0.02)
        cbar.set_label("sequence identity to query")

    stem = out.with_suffix("")
    for fmt in args.formats.split(","):
        fig.savefig(f"{stem}.{fmt.strip()}", facecolor="white",
                    bbox_inches="tight")
    plt.close(fig)
    print(f"[plot_msa_coverage] -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
