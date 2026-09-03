#!/usr/bin/env python3
"""
scripts/annotate_disorder.py — Stage 1c
=======================================
Per-residue disorder / binding-propensity tracks for every protein chain.

  * metapredict   — metapredict v3 (pip, no gate). ALWAYS ON. Primary track.
  * iupred3_long / iupred3_short / anchor2
                  — IUPred3 + ANCHOR2. OPTIONAL: academic-download gated, so
                    not on PyPI. Drop the downloaded package at tools/iupred3/
                    (or `pip install` it into envs/annotate.yaml) and these
                    tracks appear automatically.
  * aiupred       — AIUPred. OPTIONAL: gated, tools/aiupred/. No web fallback
                    (the public REST endpoint returns 403).

Output: long-format TSV  chain  resi  restype  track  score
Missing optional tools are simply skipped — the stage still completes.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib_pipeline import REPO_ROOT, load_config, load_system, protein_chains  # noqa: E402

CFG = load_config()
for _p in (REPO_ROOT / "tools" / "iupred3", REPO_ROOT / "tools" / "aiupred"):
    if _p.is_dir():
        sys.path.insert(0, str(_p))


def metapredict_track(seq: str) -> list[float] | None:
    try:
        import metapredict as meta
    except Exception as e:                               # noqa: BLE001
        print(f"[disorder] metapredict unavailable: {e}", file=sys.stderr)
        return None
    try:
        out = meta.predict_disorder(seq)
        # v2 returns list; v3 may return an object with .disorder
        return list(map(float, getattr(out, "disorder", out)))
    except Exception as e:                               # noqa: BLE001
        print(f"[disorder] metapredict failed: {e}", file=sys.stderr)
        return None


def iupred3_tracks(seq: str, types, anchor2: bool) -> dict[str, list[float]]:
    try:
        try:
            from iupred3 import iupred3_lib as lib
        except Exception:
            import iupred3_lib as lib                    # flat layout
    except Exception as e:                               # noqa: BLE001
        print(f"[disorder] IUPred3 not vendored (tools/iupred3/): {e}", file=sys.stderr)
        return {}
    out: dict[str, list[float]] = {}
    for t in types:
        try:
            out[f"iupred3_{t}"] = list(map(float, lib.iupred(seq, t)[0]))
        except Exception as e:                           # noqa: BLE001
            print(f"[disorder] IUPred3 {t} failed: {e}", file=sys.stderr)
    if anchor2:
        try:
            rl = lib.iupred(seq, "long")
            anc = lib.anchor2(seq, rl[0], rl[1] if len(rl) > 1 else None)
            out["anchor2"] = list(map(float, anc))
        except Exception as e:                           # noqa: BLE001
            print(f"[disorder] ANCHOR2 failed: {e}", file=sys.stderr)
    return out


def aiupred_track(seq: str) -> list[float] | None:
    try:
        import aiupred
    except Exception as e:                               # noqa: BLE001
        print(f"[disorder] AIUPred not vendored (tools/aiupred/): {e}", file=sys.stderr)
        return None
    try:
        if hasattr(aiupred, "predict"):
            return list(map(float, aiupred.predict(seq)))
        from aiupred import aiupred_lib
        return list(map(float, aiupred_lib.predict_disorder(seq)))
    except Exception as e:                               # noqa: BLE001
        print(f"[disorder] AIUPred failed: {e}", file=sys.stderr)
        return None


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--system", required=True)
    ap.add_argument("--spec", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    spec = load_system(args.system, Path(args.spec).resolve().parent.parent)
    acfg = CFG["annotate"]
    rows: list[tuple] = []

    for ch in protein_chains(spec):
        cid, seq = ch["id"], ch["sequence"]
        print(f"[disorder] chain {cid} ({len(seq)} aa)", flush=True)
        tracks: dict[str, list[float]] = {}
        if acfg.get("metapredict", {}).get("enabled", True):
            mp = metapredict_track(seq)
            if mp and len(mp) == len(seq):
                tracks["metapredict"] = mp
        tracks.update(iupred3_tracks(seq, acfg["iupred3"]["types"],
                                     acfg["anchor2"]["enabled"]))
        if acfg.get("aiupred", {}).get("enabled", False):
            ai = aiupred_track(seq)
            if ai and len(ai) == len(seq):
                tracks["aiupred"] = ai
        for i, aa in enumerate(seq, start=1):
            for tname, vals in tracks.items():
                if i - 1 < len(vals):
                    rows.append((cid, i, aa, tname, round(float(vals[i - 1]), 5)))

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w") as fh:
        fh.write("chain\tresi\trestype\ttrack\tscore\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")
    seen = sorted({r[3] for r in rows})
    print(f"[disorder] wrote {len(rows)} rows, tracks={seen} -> {out}")
    if not rows:
        print("[disorder] WARNING: no tracks produced — is metapredict installed "
              "(envs/annotate.yaml)?", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
