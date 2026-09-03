#!/usr/bin/env python3
"""
scripts/collect_run_stats.py — Stage 2 report tables
===================================================
Emits two CSVs from a finished Stage 2 run:

  --timings      system,backend,n_models,status,wall_clock,note
                 n_models    : rows in model_metadata.parquet for that backend
                 status      : ok | NO MODELS
                 wall_clock  : from `sacct` if a job id is discoverable, else ""
                 note        : failure signature grepped from the run log when
                               a backend produced nothing (token cap / OOM /
                               segfault / traceback)
  --compression  system,raw_mb,stored_mb,ratio,freed_mb
                 parsed from logs/processing/compress/<system>.log (the
                 compress_abcfold_metadata.py stdout), falling back to a live
                 `du` of results/metadata/<system>/.

Best-effort: anything unknown is left blank, never fatal.
"""
from __future__ import annotations

import argparse
import re
import subprocess
from pathlib import Path

import pandas as pd

FAIL_SIGNS = [
    (r"Too many tokens in input.*?(\d+) > (\d+)", "token cap: {0} > {1}"),
    (r"does not support n_token > (\d+)", "token cap: n_token > {0}"),
    (r"[Oo]ut of memory|CUDA out of memory|ran out of memory", "GPU OOM"),
    (r"SIGSEGV|Segmentation fault", "segfault"),
    (r"Traceback \(most recent call last\)", "python traceback"),
]


def _run_log_note(logdir: Path, system: str, backend: str) -> str:
    for cand in (logdir / "run_abcfold" / f"{system}.log",
                 logdir.parent / "processing" / "run_abcfold" / f"{system}.log"):
        if not cand.exists():
            continue
        txt = cand.read_text(errors="ignore")
        # narrow to the backend's section if the log delimits them
        seg = txt
        m = re.search(rf"{backend}.*?(?=\n\S+ (?:START|DONE|backend)|\Z)", txt,
                      re.IGNORECASE | re.DOTALL)
        if m:
            seg = m.group(0)
        for pat, tmpl in FAIL_SIGNS:
            mm = re.search(pat, seg)
            if mm:
                return tmpl.format(*mm.groups())
    return ""


def _sacct_wall(system: str) -> str:
    try:
        out = subprocess.run(
            ["sacct", "--name", f"run_abcfold_{system}", "--format", "Elapsed",
             "--noheader", "--parsable2"],
            capture_output=True, text=True, timeout=20).stdout.strip()
        return out.splitlines()[0] if out else ""
    except Exception:                                   # noqa: BLE001
        return ""


def timings(systems, meta_root: Path, log_root: Path, out: Path) -> None:
    rows = []
    for s in systems:
        pq = meta_root / s / "model_metadata.parquet"
        counts = {}
        if pq.exists():
            df = pd.read_parquet(pq, columns=["backend"])
            counts = df["backend"].value_counts().to_dict()
        wall = _sacct_wall(s)
        backends = counts or {b: 0 for b in
                              ["alphafold3", "boltz", "chai1", "openfold3",
                               "protenix", "rosettafold3"]}
        for b, n in sorted(backends.items()):
            rows.append({
                "system": s, "backend": b, "n_models": int(n),
                "status": "ok" if n else "NO MODELS",
                "wall_clock": wall if b == sorted(backends)[0] else "",
                "note": "" if n else _run_log_note(log_root, s, b),
            })
    out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out, index=False)
    print(f"[run_stats] timings -> {out} ({len(rows)} rows)")


_RAW_RE = re.compile(
    r"([\d.]+)MB raw -> ([\d.]+)MB .*?ratio\s*([\d.]+)?", re.IGNORECASE)
_FREED_RE = re.compile(r"freed:\s*([\d.]+)GB", re.IGNORECASE)


def _du_mb(p: Path) -> float:
    if not p.exists():
        return 0.0
    total = sum(f.stat().st_size for f in p.rglob("*") if f.is_file())
    return round(total / 1e6, 1)


def compression(systems, abc_root: Path, meta_root: Path, log_root: Path, out: Path) -> None:
    rows = []
    for s in systems:
        log = None
        for cand in (log_root / "compress" / f"{s}.log",):
            if cand.exists():
                log = cand.read_text(errors="ignore")
        raw = stored = ratio = freed = ""
        if log:
            m = _RAW_RE.search(log)
            if m:
                raw, stored = m.group(1), m.group(2)
                ratio = m.group(3) or (
                    round(float(raw) / float(stored), 1) if float(stored) else "")
            fm = _FREED_RE.search(log)
            if fm:
                freed = round(float(fm.group(1)) * 1000, 1)
        stored_live = _du_mb(meta_root / s)
        rows.append({
            "system": s,
            "raw_mb": raw,
            "stored_mb": stored or stored_live,
            "ratio": ratio,
            "freed_mb": freed,
            "abcfold_remaining_mb": _du_mb(abc_root / s),
        })
    out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out, index=False)
    print(f"[run_stats] compression -> {out} ({len(rows)} rows)")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--systems", nargs="+", required=True)
    ap.add_argument("--abcfold-root", required=True)
    ap.add_argument("--metadata-root", required=True)
    ap.add_argument("--timings", required=True)
    ap.add_argument("--compression", required=True)
    ap.add_argument("--log-root", default="logs/processing")
    a = ap.parse_args()
    timings(a.systems, Path(a.metadata_root), Path(a.log_root), Path(a.timings))
    compression(a.systems, Path(a.abcfold_root), Path(a.metadata_root),
                Path(a.log_root), Path(a.compression))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
