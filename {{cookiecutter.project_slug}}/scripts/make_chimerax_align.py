#!/usr/bin/env python3
"""
scripts/make_chimerax_align.py
==============================
Write ChimeraX command scripts (.cxc) that open the selected models and
superpose them on the same anchor rigid core pose_cluster_anchor.py used
for clustering (rigid_core.json), so what you see in ChimeraX is the frame
the pose clusters were computed in.

Outputs, in --out-dir:
  - align_selected.cxc       every selected model, all clusters
  - align_cluster_<k>.cxc    one per pose cluster
  - align_<backend>.cxc      one per backend (e.g. alphafold3, openfold3)

Each script opens the models by paths relative to itself (ChimeraX resolves
`open` paths in a .cxc against the script's folder), names the core
`rigidcore`, aligns every model's core Cα onto the first model's (the
top-ranked one), and colours by chain with the core tinted.

    chimerax results/<system>/chimerax/align_cluster_1.cxc
"""
from __future__ import annotations

import argparse
import csv
import json
import os
from collections import defaultdict
from itertools import groupby
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--rigid-core", required=True)
    p.add_argument("--selected-csv", required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--system", required=True)
    return p.parse_args()


def ranges(resids: list[int]) -> str:
    """[1,2,3,7,8,10] -> '1-3,7-8,10'"""
    out = []
    for _, grp in groupby(enumerate(sorted(resids)), key=lambda x: x[1] - x[0]):
        g = [r for _, r in grp]
        out.append(f"{g[0]}-{g[-1]}" if len(g) > 1 else f"{g[0]}")
    return ",".join(out)


def core_spec(core: list[list]) -> str:
    """[[A,1],[A,2],[B,5]] -> '/A:1-2/B:5' (ChimeraX atomspec, no model)"""
    by_chain = defaultdict(list)
    for chain, resid in core:
        by_chain[chain].append(int(resid))
    return "".join(f"/{c}:{ranges(r)}" for c, r in by_chain.items())


def write_cxc(path: Path, rows: list[dict], core: dict, system: str, label: str):
    spec = core_spec(core["core"])
    n_core = len(core["core"])
    lines = [
        f"# {system}: {label} - {len(rows)} models superposed on the anchor rigid core",
        f"# core = {n_core} Ca of chain(s) {'+'.join(core['anchor_chains'])} "
        f"(pose_cluster_anchor.py, median RMSF target {core['core_rmsf_target']} A)",
        "# #1 is the reference (top-ranked model); core RMSDs print in the log",
        "close session",
    ]
    for i, r in enumerate(rows, 1):
        rel = os.path.relpath(r["staged_cif"], path.parent)
        lines.append(f"open {rel}")
        lines.append(f"rename #{i} c{r['cluster']}_{Path(r['staged_cif']).stem}")
    lines.append(f"name frozen rigidcore {spec}")
    for i in range(2, len(rows) + 1):
        lines.append(f"align #{i} & rigidcore & @CA toAtoms #1 & rigidcore & @CA move structures")
    lines += [
        "hide atoms",
        "cartoon",
        "color bychain",
        "color rigidcore & protein lightgray target c",
        "view #1",
    ]
    path.write_text("\n".join(lines) + "\n")


def main():
    args = parse_args()
    core = json.loads(Path(args.rigid_core).read_text())
    with open(args.selected_csv) as fh:
        rows = sorted(csv.DictReader(fh), key=lambda r: (int(r["cluster"]), int(r["rank"])))
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    write_cxc(out / "align_selected.cxc", rows, core, args.system, "all selected models")
    for k, grp in groupby(rows, key=lambda r: r["cluster"]):
        write_cxc(out / f"align_cluster_{k}.cxc", list(grp), core, args.system, f"pose cluster {k}")
    for be in sorted({r["backend"] for r in rows}):
        write_cxc(out / f"align_{be}.cxc", [r for r in rows if r["backend"] == be],
                  core, args.system, f"{be} models")
    print(f"[chimerax_align] {args.system}: {len(rows)} models, "
          f"core {len(core['core'])} Ca -> {out}")


if __name__ == "__main__":
    main()
