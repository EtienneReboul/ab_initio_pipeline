#!/usr/bin/env python3
"""
scripts/abcfold_input_stats.py — Stage 1d
=========================================
Two modes, both reading data/fold_inputs/<system>/fold_input.resolved.json:

  (default)          -> input_stats.tsv   one row per chain + a TOTAL row:
                        chain type length tokens msa_depth n_templates
                        plus TOTAL with n_seeds, json_bytes, models_enabled,
                        and per-backend token-cap risk flags.
  --templates-table  -> templates.csv     one row per embedded template hit:
                        chain idx pdb_id query_start query_end n_res mmcif_bytes
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib_pipeline import BACKEND_TOKEN_CAPS, load_config, norm_id  # noqa: E402

CFG = load_config()
_NUC = {"rna", "dna"}


def _iter_chains(doc: dict):
    for entry in doc.get("sequences", []):
        for ctype in ("protein", "rna", "dna", "ligand"):
            if ctype in entry:
                yield ctype, entry[ctype]


def _msa_depth(a3m: str | None) -> int:
    if not a3m:
        return 0
    return sum(1 for ln in a3m.splitlines() if ln.startswith(">"))


def _pdb_id_from_mmcif(txt: str) -> str:
    m = re.search(r"^_entry\.id\s+(\S+)", txt, re.MULTILINE)
    if m:
        return m.group(1).strip("'\"")
    m = re.search(r"^data_(\S+)", txt, re.MULTILINE)
    return m.group(1) if m else ""


def templates_table(doc: dict, out: Path) -> None:
    rows = []
    for ctype, ch in _iter_chains(doc):
        if ctype != "protein":
            continue
        for i, t in enumerate(ch.get("templates", []) or []):
            mmcif = t.get("mmcif") or t.get("mmcifPath") or ""
            qi = t.get("queryIndices") or []
            rows.append({
                "chain": norm_id(ch.get("id", "?")),
                "idx": i,
                "pdb_id": _pdb_id_from_mmcif(mmcif) if "\n" in str(mmcif) else Path(str(mmcif)).stem,
                "query_start": (min(qi) + 1) if qi else "",
                "query_end": (max(qi) + 1) if qi else "",
                "n_res": len(t.get("templateIndices") or qi),
                "mmcif_bytes": len(mmcif) if isinstance(mmcif, str) else "",
            })
    cols = ["chain", "idx", "pdb_id", "query_start", "query_end", "n_res", "mmcif_bytes"]
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w") as fh:
        fh.write(",".join(cols) + "\n")
        for r in rows:
            fh.write(",".join(str(r[c]) for c in cols) + "\n")
    print(f"[input_stats] {len(rows)} template row(s) -> {out}")


def stats_table(doc: dict, resolved: Path, out: Path) -> None:
    rows, total_tokens = [], 0
    for ctype, ch in _iter_chains(doc):
        seq = ch.get("sequence", "") or ""
        is_poly = ctype in {"protein"} | _NUC
        tokens = len(seq) if is_poly else 1
        total_tokens += tokens
        rows.append({
            "chain": norm_id(ch.get("id", "?")), "type": ctype, "length": len(seq),
            "tokens": tokens,
            "msa_depth": _msa_depth(ch.get("unpairedMsa")) if ctype == "protein" else "",
            "n_templates": len(ch.get("templates", []) or []) if ctype == "protein" else "",
        })

    models = [m for m, on in CFG["abcfold"]["models"].items() if on]
    risk = [f"{b}(>{cap})" for b, cap in BACKEND_TOKEN_CAPS.items()
            if b in models and total_tokens > cap]

    cols = ["chain", "type", "length", "tokens", "msa_depth", "n_templates"]
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")
        fh.write("\t".join([
            "TOTAL", f"n_seeds={len(doc.get('modelSeeds', []))}",
            f"json_bytes={resolved.stat().st_size}",
            str(total_tokens),
            f"models={'+'.join(models)}",
            f"token_cap_risk={';'.join(risk) if risk else 'none'}",
        ]) + "\n")
    print(f"[input_stats] {len(rows)} chains, {total_tokens} tokens, "
          f"risk={risk or 'none'} -> {out}")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--system", required=True)
    ap.add_argument("--spec")
    ap.add_argument("--resolved", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--templates-table", action="store_true")
    args = ap.parse_args()

    resolved = Path(args.resolved)
    doc = json.loads(resolved.read_text())
    if args.templates_table:
        templates_table(doc, Path(args.out))
    else:
        stats_table(doc, resolved, Path(args.out))
    return 0


if __name__ == "__main__":
    sys.exit(main())
