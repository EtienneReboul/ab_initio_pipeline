#!/usr/bin/env python3
"""
scripts/merge_annotation.py — Stage 1c (the human-review checkpoint)
==================================================================
Fuse domains.raw.tsv + disorder.tsv + morf.tsv into ONE proposed
`domains:` block per protein chain, written to
data/annotation/<system>/annotation.yaml for a human to curate into
configs/<system>.yaml. Also emits a plain-text diff against whatever
`domains:` the config already has.

Segment logic, per chain, left to right:
  * folded-domain intervals  = merged InterProScan/PROSITE hits (profiles
    preferred over patterns; overlapping same-DB hits merged)
  * consensus disorder mask  = residues where >= `disorder_min_tracks`
    disorder tracks exceed `disorder_cutoff`
  * gaps between domains      -> kind: disordered if mostly in the mask,
    else kind: linker
  * MoRF segments (morf.segments.tsv) added as separate kind: morf entries
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import sys
from collections import defaultdict
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib_pipeline import load_config, load_system, protein_chains  # noqa: E402

CFG = load_config()
DOMAIN_DBS_PREFERRED = ("PfamA", "Pfam", "PROSITEProfiles", "PROSITE_profiles", "SMART", "Gene3D")


def _read_tsv(path: Path) -> list[dict]:
    if not path.exists():
        return []
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _merge_intervals(ivals: list[tuple[int, int, str]], gap: int = 10) -> list[tuple[int, int, str]]:
    if not ivals:
        return []
    ivals = sorted(ivals)
    out = [list(ivals[0])]
    for s, e, name in ivals[1:]:
        if s <= out[-1][1] + gap:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e, name])
    return [(s, e, n) for s, e, n in out]


def _chain_segments(cid, seqlen, dom_hits, dis_rows, morf_segs, acfg):
    # folded domains from hits
    pref = [(int(h["start"]), int(h["end"]), h["name"] or h["accession"])
            for h in dom_hits if h["db"] in DOMAIN_DBS_PREFERRED]
    allh = [(int(h["start"]), int(h["end"]), h["name"] or h["accession"]) for h in dom_hits]
    domains = _merge_intervals(pref or allh)

    # consensus disorder mask
    by_res = defaultdict(dict)
    for r in dis_rows:
        by_res[int(r["resi"])][r["track"]] = float(r["score"])
    cutoff = acfg["consensus"]["disorder_cutoff"]
    n_avail = len({r["track"] for r in dis_rows})
    # Cap at however many tracks actually ran: IUPred3/ANCHOR2/AIUPred are
    # academic-gated (not on PyPI) so a fresh checkout usually only has
    # metapredict. disorder_min_tracks (default 2) would then be
    # unsatisfiable forever, silently pinning disordered_fraction at 0.0 and
    # every gap at kind: linker regardless of actual scores. Mirrors
    # plot_disorder.py's own shading, which already uses a majority of
    # whatever tracks are present rather than a fixed count.
    need = min(acfg["consensus"]["disorder_min_tracks"], max(n_avail, 1))
    mask = {i for i in range(1, seqlen + 1)
            if sum(1 for v in by_res.get(i, {}).values() if v >= cutoff) >= need}

    segs: list[dict] = []
    cursor = 1
    di = 1
    for s, e, name in domains:
        s, e = max(s, 1), min(e, seqlen)
        if s > cursor:                       # gap before this domain
            gap_res = range(cursor, s)
            frac = sum(1 for i in gap_res if i in mask) / max(1, len(gap_res))
            segs.append({"name": f"region{di}", "start": cursor, "end": s - 1,
                         "kind": "disordered" if frac >= 0.5 else "linker"})
            di += 1
        segs.append({"name": name[:24] or f"dom{di}", "start": s, "end": e, "kind": "domain"})
        cursor = e + 1
    if cursor <= seqlen:
        gap_res = range(cursor, seqlen + 1)
        frac = sum(1 for i in gap_res if i in mask) / max(1, len(gap_res))
        segs.append({"name": f"region{di}", "start": cursor, "end": seqlen,
                     "kind": "disordered" if frac >= 0.5 else "linker"})

    for s in morf_segs:
        segs.append({"name": f"morf_{s['start']}_{s['end']}", "start": int(s["start"]),
                     "end": int(s["end"]), "kind": "morf"})

    disfrac = len(mask) / max(1, seqlen)
    return segs, {"disordered_fraction": round(disfrac, 3),
                  "tracks": sorted({r["track"] for r in dis_rows})}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--system", required=True)
    ap.add_argument("--spec", required=True)
    ap.add_argument("--domains", required=True)
    ap.add_argument("--disorder", required=True)
    ap.add_argument("--morf", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    repo = Path(args.spec).resolve().parent.parent
    spec = load_system(args.system, repo)
    acfg = CFG["annotate"]

    dom_all = _read_tsv(Path(args.domains))
    dis_all = _read_tsv(Path(args.disorder))
    morf_seg_all = _read_tsv(Path(args.morf).with_suffix(".segments.tsv"))

    proposed, summary = {}, {}
    for ch in protein_chains(spec):
        cid = ch["id"]
        segs, summ = _chain_segments(
            cid, len(ch["sequence"]),
            [h for h in dom_all if h["chain"] == cid],
            [r for r in dis_all if r["chain"] == cid],
            [s for s in morf_seg_all if s["chain"] == cid],
            acfg,
        )
        proposed[cid] = segs
        summary[cid] = summ

    # diff vs current config
    current = spec.get("domains") or {}
    diff_lines = []
    for cid in proposed:
        cur = current.get(cid)
        if not cur:
            diff_lines.append(f"  {cid}: NEW ({len(proposed[cid])} segments proposed)")
        else:
            cur_set = {(d["start"], d["end"], d["kind"]) for d in cur}
            new_set = {(d["start"], d["end"], d["kind"]) for d in proposed[cid]}
            if cur_set != new_set:
                diff_lines.append(f"  {cid}: CHANGED  (config has {len(cur)}, proposal {len(proposed[cid])})")
            else:
                diff_lines.append(f"  {cid}: unchanged")

    doc = {
        "system": args.system,
        "generated": dt.datetime.now().isoformat(timespec="seconds"),
        "domains": proposed,
        "disorder_summary": summary,
        "review": (
            "Curate these segments, paste the `domains:` block into "
            f"configs/{args.system}.yaml, then set `annotation_reviewed: true`.\n"
            "Diff vs current config:\n" + "\n".join(diff_lines)
        ),
    }
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(yaml.safe_dump(doc, sort_keys=False, width=100))
    print(f"[merge_annotation] -> {out}")
    print("\n".join(diff_lines))
    if not spec.get("annotation_reviewed"):
        print(f"\n[merge_annotation] REVIEW REQUIRED: configs/{args.system}.yaml still has "
              "annotation_reviewed: false — stage 2 will refuse to run this system.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
