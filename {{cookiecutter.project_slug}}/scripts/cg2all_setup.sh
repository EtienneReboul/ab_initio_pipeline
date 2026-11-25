#!/usr/bin/env bash
# Fetch cg2all's checkpoint (~49 MB, Zenodo, into the active env) by rebuilding a
# 12-ALA helix, then write that env's bin directory to $1.
set -euo pipefail
out=$1
bin=$(dirname "$(command -v convert_cg2all)")
tmp=$(mktemp -d "${TMPDIR:-/tmp}/cg2all.XXXXXX")
trap 'rm -rf "$tmp"' EXIT
python - "$tmp/ca.pdb" <<'PY'
import math, sys
with open(sys.argv[1], "w") as fh:
    for i in range(12):
        a = math.radians(100 * i)
        fh.write("ATOM  %5d  CA  ALA A%4d    %8.3f%8.3f%8.3f  1.00  0.00           C\n"
                 % (i + 1, i + 1, 2.3 * math.cos(a), 2.3 * math.sin(a), 1.5 * i))
    fh.write("END\n")
PY
DGLBACKEND=pytorch convert_cg2all -p "$tmp/ca.pdb" -o "$tmp/all.pdb" --cg CalphaBasedModel --fix --device cpu
[ -s "$tmp/all.pdb" ]
mkdir -p "$(dirname "$out")"
echo "$bin" > "$out"
