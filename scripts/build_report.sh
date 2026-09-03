#!/usr/bin/env bash
# scripts/build_report.sh <stage> [repo_root]
# Rebuilds reports/<stage>.zip from the just-completed Snakemake run.
# Invoked from each stage Snakefile's `onsuccess:` hook, so it must never
# fail the pipeline — a broken report build is logged and swallowed.
#
# The report is a .zip (not a bare .html) because `snakemake --report` only
# embeds interactive HTML items — the datavzrd table bundles — for the zip
# form. Unzip it and open <stage>/report.html; the dark theme lives in
# report/custom.css, injected via --report-stylesheet.
set -uo pipefail

STAGE="${1:?usage: build_report.sh <preprocessing|processing|postprocessing> [repo_root]}"
REPO="${2:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
SNAKEFILE="$REPO/workflows/$STAGE/Snakefile"
OUT="$REPO/reports/${STAGE}.zip"
CSS="$REPO/report/custom.css"

mkdir -p "$REPO/reports"
echo "[build_report] $STAGE -> $OUT"

args=(-s "$SNAKEFILE" --report "$OUT")
[[ -f "$CSS" ]] && args+=(--report-stylesheet "$CSS")

if command -v snakemake >/dev/null 2>&1; then
    ( cd "$REPO" && snakemake "${args[@]}" ) \
        && echo "[build_report] ok: $OUT" \
        || echo "[build_report] WARNING: report build failed for $STAGE (non-fatal)" >&2
else
    echo "[build_report] WARNING: snakemake not on PATH — skipping report build" >&2
fi
exit 0
