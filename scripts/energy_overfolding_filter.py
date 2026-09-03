#!/usr/bin/env python3
"""
scripts/energy_overfolding_filter.py — Stage 3 (optional analysis)
================================================================
Port of the DRB2 pipeline's
notebooks/rna_complexes_energy_and_overfolding_filter.ipynb.

Two successive filters over the minimized ensemble:

  1. Garbage-energy filter — MAD-z on each model's final ChimeraX
     minimization energy; a model is an energy outlier if |z| > --energy-z.
     A whole backend is dropped if > --backend-drop-frac of its assessed
     models are outliers (RosettaFold3 numeric blow-up in the DRB2 runs).
     Models with no parseable energy are kept.

  2. Over-folding filter — for the NON-MoRF residues of every `kind:
     disordered` segment in configs/<system>.yaml `domains:` (i.e. the
     disordered span minus any overlapping `kind: morf` sub-segment), drop a
     model if those residues show > --helix-frac helix OR a contiguous helix
     run >= --helix-run. Needs dssp_summary.csv (fname,cluster,chain,resnum,
     restype,ss_type) from scripts/dssp_summary.py over the minimized PDBs.

Output:
  results/<system>/filtered_models.csv           (per-model keep/drop + reasons)
  results/<system>/figures/domain_analysis/energy_overfolding_filter/*.svg
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

BACKEND_RE = re.compile(r"(alphafold3|boltz|chai1|openfold3|protenix|rosettafold3)")


def non_morf_disordered(spec) -> dict[str, list[tuple[int, int]]]:
    out: dict[str, list[tuple[int, int]]] = {}
    for chain, segs in (spec.get("domains") or {}).items():
        morfs = [(int(d["start"]), int(d["end"])) for d in segs if d["kind"] == "morf"]
        spans = []
        for d in segs:
            if d["kind"] != "disordered":
                continue
            res = set(range(int(d["start"]), int(d["end"]) + 1))
            for ms, me in morfs:
                res -= set(range(ms, me + 1))
            if res:
                r = sorted(res)
                # collapse to contiguous ranges
                start = prev = r[0]
                for x in r[1:]:
                    if x != prev + 1:
                        spans.append((start, prev)); start = x
                    prev = x
                spans.append((start, prev))
        if spans:
            out[chain] = spans
    return out


def final_energies(results_root: Path, system: str) -> pd.DataFrame:
    rows = []
    for csv in (results_root / system / "minimized").rglob("*_energy.csv"):
        bm = BACKEND_RE.search(str(csv))
        try:
            d = pd.read_csv(csv)
        except Exception:                                  # noqa: BLE001
            continue
        ecol = next((c for c in d.columns if "energy" in c.lower()), None)
        e = pd.to_numeric(d[ecol], errors="coerce").dropna() if ecol else pd.Series(dtype=float)
        rows.append({"model": csv.stem.replace("_energy", ""),
                     "backend": bm.group(1) if bm else "unknown",
                     "final_energy": float(e.iloc[-1]) if len(e) else np.nan})
    return pd.DataFrame(rows)


def mad_z(x):
    x = np.asarray(x, float)
    med = np.nanmedian(x)
    mad = np.nanmedian(np.abs(x - med)) or 1e-9
    return 0.6745 * (x - med) / mad


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--system", required=True)
    ap.add_argument("--spec", required=True)
    ap.add_argument("--results-root", default="results")
    ap.add_argument("--dssp", default="", help="dssp_summary.csv (optional; skips filter 2 if absent)")
    ap.add_argument("--energy-z", type=float, default=3.5)
    ap.add_argument("--backend-drop-frac", type=float, default=0.5)
    ap.add_argument("--helix-frac", type=float, default=0.20)
    ap.add_argument("--helix-run", type=int, default=12)
    ap.add_argument("--out", required=True)
    ap.add_argument("--formats", default="svg")
    a = ap.parse_args()

    spec = yaml.safe_load(Path(a.spec).read_text())
    rroot = Path(a.results_root)
    en = final_energies(rroot, a.system)
    if en.empty:
        raise SystemExit(f"[overfold] no *_energy.csv under {rroot / a.system / 'minimized'}")

    # ── filter 1: garbage energy ───────────────────────────────────────────
    assessable = en["final_energy"].notna()
    en["energy_z"] = np.nan
    en.loc[assessable, "energy_z"] = mad_z(en.loc[assessable, "final_energy"])
    en["energy_outlier"] = en["energy_z"].abs() > a.energy_z
    drop_backends = set()
    for b, g in en[assessable].groupby("backend"):
        if len(g) and g["energy_outlier"].mean() > a.backend_drop_frac:
            drop_backends.add(b)
    en["drop_energy"] = en["backend"].isin(drop_backends) | en["energy_outlier"].fillna(False)
    print(f"[overfold] filter 1: dropped backends={sorted(drop_backends) or 'none'}; "
          f"{int(en['drop_energy'].sum())}/{len(en)} models flagged")

    # ── filter 2: over-folding of non-MoRF disordered residues ─────────────
    en["drop_overfold"] = False
    en["overfold_detail"] = ""
    spans = non_morf_disordered(spec)
    if a.dssp and Path(a.dssp).exists() and spans:
        dssp = pd.read_csv(a.dssp)
        for model, md in dssp.groupby("fname"):
            hits = []
            for chain, rngs in spans.items():
                sub = md[(md["chain"] == chain)]
                for s, e in rngs:
                    seg = sub[(sub["resnum"] >= s) & (sub["resnum"] <= e)].sort_values("resnum")
                    if seg.empty:
                        continue
                    ishelix = (seg["ss_type"] == "helix").to_numpy()
                    frac = ishelix.mean()
                    run = 0; best = 0
                    for v in ishelix:
                        run = run + 1 if v else 0
                        best = max(best, run)
                    if frac > a.helix_frac or best >= a.helix_run:
                        hits.append(f"{chain}:{s}-{e} helix{frac:.0%}/run{best}")
            if hits:
                en.loc[en["model"] == model, "drop_overfold"] = True
                en.loc[en["model"] == model, "overfold_detail"] = "; ".join(hits)
        print(f"[overfold] filter 2: {int(en['drop_overfold'].sum())}/{len(en)} models over-fold")
    else:
        print("[overfold] filter 2 skipped (no dssp_summary.csv or no disordered segments)")

    en["keep"] = ~(en["drop_energy"] | en["drop_overfold"])
    out = Path(a.out); out.parent.mkdir(parents=True, exist_ok=True)
    en.to_csv(out, index=False)
    print(f"[overfold] {int(en['keep'].sum())}/{len(en)} models survive -> {out}")

    figdir = rroot / a.system / "figures" / "domain_analysis" / "energy_overfolding_filter"
    figdir.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4), constrained_layout=True)
    for b, g in en.groupby("backend"):
        axes[0].scatter(g["backend"], g["final_energy"], s=12, alpha=.6)
    axes[0].set_yscale("symlog"); axes[0].set_ylabel("final energy (kJ/mol)")
    axes[0].set_title("filter 1: minimization energy by backend")
    axes[0].tick_params(axis="x", rotation=30)
    surv = en.groupby("backend")["keep"].sum()
    tot = en.groupby("backend")["keep"].count()
    axes[1].bar(surv.index, tot, color="#dddddd", label="assessed")
    axes[1].bar(surv.index, surv, color="#55a868", label="survive")
    axes[1].set_title("models surviving both filters"); axes[1].legend(frameon=False)
    axes[1].tick_params(axis="x", rotation=30)
    for fmt in a.formats.split(","):
        fig.savefig(figdir / f"summary.{fmt.strip()}")
    plt.close(fig)
    print(f"[overfold] figures -> {figdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
