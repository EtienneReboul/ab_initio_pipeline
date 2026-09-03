#!/usr/bin/env python3
"""
scripts/annotate_morf.py — Stage 1c
===================================
Run MoRFchibi 2.0 (vendored tools/MC2) on every protein chain and emit a
per-residue MoRF-probability track plus called segments.

MC2 is a run-in-place script (hardcoded param.py / input.fasta / relative
Models,Stuff paths, writes output/<AC>.caid). We therefore stage a scratch
copy of tools/MC2 under data/annotation/<system>/_mc2/, patch its param.py
`aff_path` to the vendored tools/AFF, drop our FASTA in, run it there, and
read back output/<chain>.caid (format: `resi<TAB>restype<TAB>score`).

Output:
  morf.tsv       long-format  chain  resi  restype  score
  morf.segments.tsv (sidecar) chain  start  end  length  mean_score
"""
from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib_pipeline import REPO_ROOT, load_system, protein_chains  # noqa: E402


def _segments(scores: list[float], thr: float, min_len: int) -> list[tuple[int, int, float]]:
    segs, run = [], []
    for i, s in enumerate(scores, start=1):
        if s >= thr:
            run.append((i, s))
        elif run:
            if len(run) >= min_len:
                segs.append((run[0][0], run[-1][0], sum(x[1] for x in run) / len(run)))
            run = []
    if run and len(run) >= min_len:
        segs.append((run[0][0], run[-1][0], sum(x[1] for x in run) / len(run)))
    return segs


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--system", required=True)
    ap.add_argument("--spec", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--score-threshold", type=float, default=0.725)
    ap.add_argument("--min-len", type=int, default=5)
    args = ap.parse_args()

    spec = load_system(args.system, Path(args.spec).resolve().parent.parent)
    prots = protein_chains(spec)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    work = out.parent / "_mc2"
    if work.exists():
        shutil.rmtree(work)
    shutil.copytree(REPO_ROOT / "tools" / "MC2", work,
                    ignore=shutil.ignore_patterns("__pycache__", "output", ".git"))
    aff = (REPO_ROOT / "tools" / "AFF").resolve()
    param = work / "param.py"
    lines = param.read_text().splitlines()
    lines = [f"aff_path = '{aff}/'" if ln.startswith("aff_path") else ln
             for ln in lines]
    param.write_text("\n".join(lines) + "\n")
    (work / "input.fasta").write_text(
        "".join(f">{c['id']}\n{c['sequence']}\n" for c in prots)
    )

    print(f"[morf] running MoRFchibi2 on {len(prots)} chain(s) in {work}", flush=True)
    r = subprocess.run([sys.executable, "MC2.py"], cwd=work,
                       capture_output=True, text=True)
    print(r.stdout)
    if r.returncode != 0:
        print(r.stderr, file=sys.stderr)
        raise SystemExit(f"MC2.py exited {r.returncode}")

    rows, segrows = [], []
    for c in prots:
        caid = work / "output" / f"{c['id']}.caid"
        if not caid.exists():
            print(f"[morf] WARNING: no output for chain {c['id']}", file=sys.stderr)
            continue
        scores = []
        for line in caid.read_text().splitlines():
            if line.startswith(">") or not line.strip():
                continue
            p = line.split("\t")
            resi, aa, sc = int(p[0]), p[1], float(p[2])
            scores.append(sc)
            rows.append((c["id"], resi, aa, round(sc, 5)))
        for s, e, m in _segments(scores, args.score_threshold, args.min_len):
            segrows.append((c["id"], s, e, e - s + 1, round(m, 5)))

    with out.open("w") as fh:
        fh.write("chain\tresi\trestype\tscore\n")
        for x in rows:
            fh.write("\t".join(map(str, x)) + "\n")
    seg = out.with_suffix(".segments.tsv")
    with seg.open("w") as fh:
        fh.write("chain\tstart\tend\tlength\tmean_score\n")
        for x in segrows:
            fh.write("\t".join(map(str, x)) + "\n")
    print(f"[morf] wrote {len(rows)} residues, {len(segrows)} MoRF segment(s) -> {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
