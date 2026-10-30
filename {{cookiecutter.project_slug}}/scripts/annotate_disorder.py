#!/usr/bin/env python3
"""
scripts/annotate_disorder.py — Stage 1c
=======================================
Per-residue disorder, MoRF/binding-propensity and linker tracks for every
protein chain, from AIUPred's docker image (ghcr.io/doszilab/aiupred) — one
container run over every chain at once (multi-FASTA; the networks load once).
Ported from NanoPulldown (nanopulldown.annotate.disorder). Replaces the old
metapredict / IUPred3+ANCHOR2 / MoRFchibi2 split: AIUPred's binding mode is
what merge_annotation.py calls MoRFs from.

  * disorder  — AIUPred disorder score
  * binding   — AIUPred binding (ANCHOR2-like) score -> MoRF segments
  * linker    — AIUPred linker score (drawn, never a segment call by itself)

Fallback: if docker / the image is unavailable, metapredict v3 (pip, in
envs/annotate.yaml) supplies a `disorder` track only — no MoRFs then. The
stage still completes; the log says which predictor ran.

Container contract (ghcr.io/doszilab/aiupred:cpu v2.1.0; the image has no
ENTRYPOINT, so it runs as `python3 -m aiupred.cli`):
  python3 -m aiupred.cli -i in.fasta -o out.tsv -b -l [--force-cpu]
Output: per record a `#>id` line, `#` comment lines, then
  resi  restype  disorder  binding  linker

Output: long-format TSV  chain  resi  restype  track  score
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
import tempfile
from pathlib import Path

# metapredict's pip/PyPI torch wheel bundles its own libomp.dylib, which
# conflicts with conda-forge numpy's libomp on macOS and aborts the process
# (OMP Error #15) the moment `import metapredict` triggers torch's load.
# Must be set before that import, not just before torch's own import.
os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib_pipeline import load_config, load_system, protein_chains  # noqa: E402

CFG = load_config()
TRACKS = ("disorder", "binding", "linker")     # AIUPred output columns 3..5


def aiupred_command(acfg: dict, workdir: Path) -> list[str]:
    gpu = bool(acfg.get("use_gpu"))
    cmd = ["docker", "run", "--rm", "-v", f"{workdir}:/work"]
    if acfg.get("docker_platform"):
        cmd += ["--platform", acfg["docker_platform"]]
    if gpu:
        cmd += ["--gpus", "all"]
    cmd += ["--entrypoint", "python3", acfg["image_gpu"] if gpu else acfg["image_cpu"],
            "-m", "aiupred.cli", "-i", "/work/input.fasta", "-o", "/work/aiupred.tsv", "-b", "-l"]
    return cmd + ([] if gpu else ["--force-cpu"])


def parse_aiupred(text: str) -> dict[str, dict[str, list[float]]]:
    """{record id: {track: [score per residue]}} from AIUPred's multi-record TSV."""
    out: dict[str, dict[str, list[float]]] = {}
    cur = None
    for line in text.splitlines():
        if line.startswith("#>"):
            cur = line[2:].strip()
            out[cur] = {t: [] for t in TRACKS}
            continue
        line = line.strip()
        if not line or line.startswith("#") or cur is None:
            continue
        parts = line.split()
        for t, v in zip(TRACKS, parts[2:5]):
            try:
                out[cur][t].append(float(v))
            except ValueError:
                pass
    return {k: {t: v for t, v in tr.items() if v} for k, tr in out.items()}


def aiupred_tracks(chains: list[dict], acfg: dict) -> dict[str, dict[str, list[float]]] | None:
    """{chain id: tracks} from ONE container run, or None if AIUPred can't run."""
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        (td / "input.fasta").write_text("".join(f">{c['id']}\n{c['sequence']}\n" for c in chains))
        cmd = aiupred_command(acfg, td)
        print(f"[disorder] {' '.join(cmd)}", flush=True)
        try:
            proc = subprocess.run(cmd, capture_output=True, text=True,
                                  timeout=int(acfg.get("timeout_s", 3600)))
        except (FileNotFoundError, subprocess.TimeoutExpired) as e:
            print(f"[disorder] AIUPred docker run failed: {e}", file=sys.stderr)
            return None
        out = td / "aiupred.tsv"
        if proc.returncode != 0 or not out.exists():
            print(f"[disorder] AIUPred exited {proc.returncode}: {proc.stderr[-500:]}", file=sys.stderr)
            return None
        parsed = parse_aiupred(out.read_text())
    for c in chains:
        n = len(parsed.get(c["id"], {}).get("disorder", []))
        if n != len(c["sequence"]):
            print(f"[disorder] AIUPred chain {c['id']}: {n} scores for {len(c['sequence'])} residues",
                  file=sys.stderr)
            return None
    return parsed


def metapredict_tracks(chains: list[dict]) -> dict[str, dict[str, list[float]]]:
    try:
        import metapredict as meta
    except Exception as e:                                  # noqa: BLE001
        print(f"[disorder] metapredict unavailable: {e}", file=sys.stderr)
        return {}
    out = {}
    for c in chains:
        try:
            out[c["id"]] = {"disorder": [float(x) for x in meta.predict_disorder(c["sequence"])]}
        except Exception as e:                              # noqa: BLE001
            print(f"[disorder] metapredict failed on chain {c['id']}: {e}", file=sys.stderr)
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--system", required=True)
    ap.add_argument("--spec", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    spec = load_system(args.system, Path(args.spec).resolve().parent.parent)
    chains = protein_chains(spec)
    acfg = CFG["annotate"]["aiupred"]

    tracks = aiupred_tracks(chains, acfg) if acfg.get("enabled", True) else None
    if tracks is None:
        print("[disorder] falling back to metapredict (disorder track only, no MoRFs)", file=sys.stderr)
        tracks = metapredict_tracks(chains)
    else:
        print("[disorder] AIUPred: disorder + binding + linker tracks", flush=True)

    rows: list[tuple] = []
    for c in chains:
        for i, aa in enumerate(c["sequence"], start=1):
            for tname, vals in tracks.get(c["id"], {}).items():
                if i - 1 < len(vals):
                    rows.append((c["id"], i, aa, tname, round(vals[i - 1], 5)))

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w") as fh:
        fh.write("chain\tresi\trestype\ttrack\tscore\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")
    print(f"[disorder] wrote {len(rows)} rows, tracks={sorted({r[3] for r in rows})} -> {out}")
    if not rows:
        print("[disorder] WARNING: no tracks produced — is docker running (AIUPred image "
              f"{acfg['image_cpu']}) or metapredict installed (envs/annotate.yaml)?", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
