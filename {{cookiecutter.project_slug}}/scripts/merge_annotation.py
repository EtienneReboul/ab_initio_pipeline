#!/usr/bin/env python3
"""
scripts/merge_annotation.py — Stage 1c
======================================
Fuse domains.raw.tsv + features.json (annotate_domains.py) and disorder.tsv
(annotate_disorder.py) into ONE `domains:` map per protein chain, written to
data/annotation/<system>/annotation.yaml. Stage 3 uses it directly
(lib_pipeline.system_domains) — no human review step. A hand-written
`domains:` block in configs/<system>.yaml still takes precedence, chain by
chain; the diff against it is printed and stored under `review:`.

Ported from NanoPulldown (nanopulldown.annotate.merge). Segment logic, per
chain, left to right (select_domains):
  * folded domains = InterPro's representative domain set — the regions its
    protein viewer summarises (`representative` flag from the InterPro API;
    for InterProScan-fallback hits, which carry no flag, the same selection
    is recomputed here: overlapping candidates grouped, the combination with
    no pair overlapping >= 30% of the shorter that covers most residues wins,
    ties to more Pfam — InterProScan 6's ProcessReprLocations) ...
  * ... plus InterPro domain / repeat ENTRIES no selected region covers
  * each region is NAMED after the InterPro domain/repeat entry that best
    matches it (overlap >= half the shorter; best Jaccard x member-database
    consensus), so a Gene3D superfamily region reads "PAZ domain", not
    "PAZ domain superfamily"; else its own entry or signature name
  * no InterPro data at all -> ScanProsite hits, merged
  * disorder mask           = residues whose "disorder" track >= disorder_cutoff
  * gaps between domains    -> split into disordered / linker runs along the
                               mask (runs < consensus.min_region_len merged into
                               a neighbour, each merged run called by majority)
  * AIUPred binding-track runs (>= morf.score_threshold, >= morf.min_len
    residues) -> separate kind: morf entries (none on the metapredict fallback)

annotation.yaml also gets an `interpro:` block per chain (family, entries,
unintegrated signatures, feature tracks) for plot_domain_bar.py.
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib_pipeline import load_config, load_system, protein_chains  # noqa: E402

CFG = load_config()
# member DBs whose domain/repeat/homologous-superfamily signatures are
# representative-domain candidates (InterPro: Representative Domain and Family
# selection), as the InterPro API / InterProScan TSV label them
DOMAIN_DBS = {"Pfam", "CDD", "PROSITE profiles", "ProSiteProfiles", "SMART", "NCBIfam", "NCBIFAM",
              "CATH-Gene3D", "Gene3D", "SUPERFAMILY"}
FAMILY_DBS = {"PANTHER", "PIRSF", "HAMAP", "NCBIfam", "NCBIFAM", "SFLD", "Pfam"}
DOMAIN_TYPES = {"domain", "repeat", "homologous_superfamily"}
OVERLAP = 0.3                 # InterProScan 6: max overlap (fraction of the shorter) inside a set
MAX_PER_GROUP = 20
_ACCESSION_RE = re.compile(r"^(IPR\d+|PF\d+|SM\d+|PS\d+|PTHR\d+|G3DSA:.*|cd\d+)$")


def _read_tsv(path: Path) -> list[dict]:
    if not path.exists():
        return []
    with path.open() as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    for r in rows:
        r["start"], r["end"] = int(r["start"]), int(r["end"])
        for k in ("type", "ipr", "ipr_name"):
            r[k] = r.get(k) or None
        rep = r.get("representative")
        r["representative"] = None if rep in (None, "") else rep == "True"
    return rows


def _label(h: dict) -> str:
    """Human-readable hit name: `name` unless it is just an accession."""
    name = (h.get("name") or "").strip()
    if name and not _ACCESSION_RE.match(name):
        return name
    desc = (h.get("description") or "").strip()
    return (desc if desc not in ("", "-") else "") or name or h["accession"]


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


def _runs(by_res: dict[int, float], threshold: float, min_len: int) -> list[tuple[int, int]]:
    """Contiguous residue runs with score >= threshold, at least min_len long."""
    segs, start, prev = [], None, None
    for i in sorted(by_res):
        hit = by_res[i] >= threshold
        if hit and start is None:
            start = i
        if not hit and start is not None:
            if prev - start + 1 >= min_len:
                segs.append((start, prev))
            start = None
        prev = i
    if start is not None and prev - start + 1 >= min_len:
        segs.append((start, prev))
    return segs


def _ov(a, b) -> int:
    return max(0, min(a[1], b[1]) - max(a[0], b[0]) + 1)


def _frac_shorter(a, b) -> float:
    return _ov(a, b) / max(1, min(a[1] - a[0] + 1, b[1] - b[0] + 1))


def _jaccard(a, b) -> float:
    o = _ov(a, b)
    return o / max(1, (a[1] - a[0] + 1) + (b[1] - b[0] + 1) - o)


def representative_set(cands: list[dict]) -> list[dict]:
    """InterProScan 6's representative-location selection (ProcessReprLocations):
    group overlapping candidates; per group keep the 20 longest, then the
    combination whose members pairwise overlap < 30% of the shorter and that
    covers the most residues (ties: more Pfam)."""
    cands = sorted(cands, key=lambda h: (h["start"], h["end"]))
    groups, cur, stop = [], [], -1
    for h in cands:
        if cur and h["start"] > stop:
            groups.append(cur)
            cur, stop = [], -1
        cur.append(h)
        stop = max(stop, h["end"])
    if cur:
        groups.append(cur)
    chosen = []
    for g in groups:
        g = sorted(g, key=lambda h: (-(h["end"] - h["start"]), h["db"] != "Pfam"))[:MAX_PER_GROUP]
        best, best_key = [], (-1, -1)

        def rec(i, picked):
            nonlocal best, best_key
            if i == len(g):
                cov = len({r for h in picked for r in range(h["start"], h["end"] + 1)})
                key = (cov, sum(h["db"] == "Pfam" for h in picked))
                if key > best_key:
                    best, best_key = list(picked), key
                return
            h = g[i]
            if all(_frac_shorter((h["start"], h["end"]), (o["start"], o["end"])) < OVERLAP for o in picked):
                rec(i + 1, picked + [h])
            rec(i + 1, picked)
        rec(0, [])
        chosen += best
    return chosen


def _entry_members(hits: list[dict]) -> dict[str, set[str]]:
    m: dict[str, set[str]] = defaultdict(set)
    for h in hits:
        if h.get("ipr") and h["db"] != "InterPro":
            m[h["ipr"]].add(h["db"])
    return m


def _region_name(region, sig: dict | None, entries: list[dict], members: dict[str, set[str]]) -> str:
    """Best-matching InterPro domain/repeat entry; else the region's own entry
    (sites keep the signature's name, e.g. 'Jas motif'), else the signature."""
    best, best_score = None, 0.0
    for e in entries:
        if e.get("type") not in ("domain", "repeat"):
            continue
        span = (e["start"], e["end"])
        if _frac_shorter(region, span) < 0.5:
            continue
        score = _jaccard(region, span) * (1 + 0.1 * len(members.get(e["accession"], ())))
        if score > best_score:
            best, best_score = e, score
    if best:
        return best["name"]
    if sig:
        own = next((e for e in entries if e["accession"] == sig.get("ipr")), None)
        if own and own.get("type") != "site":
            return own["name"]
        if sig.get("ipr_name") and not own and sig.get("type") != "site":
            return sig["ipr_name"]
        return _label(sig)
    return "domain"


def select_domains(dom_hits: list[dict]) -> list[tuple[int, int, str]]:
    """(start, end, name) folded-domain regions of one chain (module docstring)."""
    ip = [h for h in dom_hits if h.get("source") in ("InterPro", "InterProScan")]
    if not ip:                                   # ScanProsite only
        return _merge_intervals([(h["start"], h["end"], _label(h)) for h in dom_hits])
    entries = [h for h in ip if h["db"] == "InterPro"]
    sigs = [h for h in ip if h["db"] != "InterPro"]
    members = _entry_members(sigs)
    if any(h.get("representative") is not None for h in sigs):
        rep = [h for h in sigs if h.get("representative") and h.get("type") != "family"]
    else:                                        # InterProScan TSV: no flags -> recompute
        rep = representative_set([h for h in sigs if h["db"] in DOMAIN_DBS
                                  and (h.get("type") in DOMAIN_TYPES or h.get("type") is None)])
    regions = [((h["start"], h["end"]), h) for h in rep]
    for e in sorted(entries, key=lambda e: -(e["end"] - e["start"])):
        span = (e["start"], e["end"])
        if e.get("type") in ("domain", "repeat") and all(
                _frac_shorter(span, r) < OVERLAP for r, _ in regions):
            regions.append((span, None))
    return sorted((s, e, _region_name((s, e), sig, entries, members)) for (s, e), sig in regions)


def interpro_tracks(dom_hits: list[dict], features: dict) -> dict:
    """The InterPro-viewer view of one chain, for plot_domain_bar.py."""
    ip = [h for h in dom_hits if h.get("source") in ("InterPro", "InterProScan")]
    entries = [h for h in ip if h["db"] == "InterPro"]
    sigs = [h for h in ip if h["db"] != "InterPro"]
    members = _entry_members(sigs)

    def span(h, **kw):
        return {"start": h["start"], "end": h["end"], **kw}
    fam_entries = [span(e, name=e["name"], accession=e["accession"]) for e in entries if e.get("type") == "family"]
    fam_sigs = [h for h in sigs if h.get("type") == "family" and h["db"] in FAMILY_DBS]
    rep_fam = [h for h in fam_sigs if h.get("representative")] or fam_sigs
    by_entry: dict[str, dict] = {}
    for e in entries:
        if e.get("type") == "family":
            continue
        d = by_entry.setdefault(e["accession"], {"accession": e["accession"], "name": e["name"],
                                                 "type": e.get("type"),
                                                 "members": sorted(members.get(e["accession"], ())),
                                                 "locations": []})
        d["locations"].append([e["start"], e["end"]])
    unint = [span(h, name=_label(h), accession=h["accession"], db=h["db"], type=h.get("type"))
             for h in sigs if not h.get("ipr") and h.get("type") != "family"]
    return {
        "family": fam_entries or [span(h, name=_label(h), accession=h["accession"]) for h in rep_fam],
        "family_sources": sorted({f"{h['db']} {h['accession']}" for h in rep_fam}),
        "entries": sorted(by_entry.values(), key=lambda d: min(l[0] for l in d["locations"])),
        "unintegrated": sorted(unint, key=lambda d: d["start"]),
        "features": {k: [list(map(int, l)) for l in v] for k, v in (features or {}).items()},
    }


def chain_segments(seqlen: int, dom_hits: list[dict], dis_rows: list[dict], acfg: dict):
    by_res: dict[str, dict[int, float]] = defaultdict(dict)
    for r in dis_rows:
        by_res[r["track"]][int(r["resi"])] = float(r["score"])
    cutoff = acfg["consensus"]["disorder_cutoff"]
    mask = {i for i, v in by_res.get("disorder", {}).items() if v >= cutoff}

    min_run = int(acfg["consensus"].get("min_region_len", 10))

    def gap(s, e):
        """A gap between domains, split into its disordered / ordered runs of
        the AIUPred mask: a run shorter than min_run is absorbed into the run
        before it (the first one into the next), then each run is called by
        majority and same-kind neighbours are joined."""
        runs: list[list[int]] = []
        for i in range(s, e + 1):
            if runs and (i in mask) == (runs[-1][1] in mask):
                runs[-1][1] = i
            else:
                runs.append([i, i])
        merged: list[list[int]] = []
        for a, b in runs:
            if merged and b - a + 1 < min_run:
                merged[-1][1] = b
            elif merged and merged[-1][1] - merged[-1][0] + 1 < min_run:
                merged[-1][1] = b              # a short first run joins the next
            else:
                merged.append([a, b])
        out: list[dict] = []
        for a, b in merged:
            frac = sum(1 for i in range(a, b + 1) if i in mask) / (b - a + 1)
            kind = "disordered" if frac >= 0.5 else "linker"
            if out and out[-1]["kind"] == kind:
                out[-1]["end"] = b
            else:
                out.append({"start": a, "end": b, "kind": kind})
        return out

    segs: list[dict] = []
    cursor = 1
    for s, e, name in select_domains(dom_hits):
        s, e = max(s, 1, cursor), min(e, seqlen)    # trim the small overlaps a set may keep
        if e < s:
            continue
        if s > cursor:
            segs += gap(cursor, s - 1)
        segs.append({"name": name[:60] or "domain", "start": s, "end": e, "kind": "domain"})
        cursor = e + 1
    if cursor <= seqlen:
        segs += gap(cursor, seqlen)
    for g in segs:                             # name gaps by kind and span
        g.setdefault("name", f"{g['kind']}_{g['start']}_{g['end']}")
    segs = [{"name": g["name"], "start": g["start"], "end": g["end"], "kind": g["kind"]} for g in segs]

    mcfg = acfg["morf"]
    for s, e in _runs(by_res.get("binding", {}), mcfg["score_threshold"], mcfg["min_len"]):
        segs.append({"name": f"morf_{s}_{e}", "start": s, "end": e, "kind": "morf"})
    return segs, {"disordered_fraction": round(len(mask) / max(1, seqlen), 3),
                  "tracks": sorted(by_res)}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--system", required=True)
    ap.add_argument("--spec", required=True)
    ap.add_argument("--domains", required=True)
    ap.add_argument("--features", required=True)
    ap.add_argument("--disorder", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    spec = load_system(args.system, Path(args.spec).resolve().parent.parent)
    acfg = CFG["annotate"]
    dom_all = _read_tsv(Path(args.domains))
    dis_all = list(csv.DictReader(Path(args.disorder).open(), delimiter="\t"))
    feats = json.loads(Path(args.features).read_text())

    proposed, summary, tracks = {}, {}, {}
    for ch in protein_chains(spec):
        cid = ch["id"]
        hits = [h for h in dom_all if h["chain"] == cid]
        tracks[cid] = interpro_tracks(hits, feats.get(cid))
        proposed[cid], summary[cid] = chain_segments(
            len(ch["sequence"]), hits, [r for r in dis_all if r["chain"] == cid], acfg)

    current = spec.get("domains") or {}
    diff_lines = []
    for cid, segs in proposed.items():
        cur = current.get(cid)
        if not cur:
            diff_lines.append(f"  {cid}: annotated ({len(segs)} segments, used as is)")
        elif {(d["start"], d["end"], d["kind"]) for d in cur} != {(d["start"], d["end"], d["kind"]) for d in segs}:
            diff_lines.append(f"  {cid}: OVERRIDDEN by configs/{args.system}.yaml "
                              f"({len(cur)} curated vs {len(segs)} annotated segments)")
        else:
            diff_lines.append(f"  {cid}: curated == annotated")

    doc = {
        "system": args.system,
        "generated": dt.datetime.now().isoformat(timespec="seconds"),
        "domains": proposed,
        "disorder_summary": summary,
        "interpro": tracks,
        "review": ("Used directly by stage 3. To override a chain, put its segments under "
                   f"`domains:` in configs/{args.system}.yaml.\n" + "\n".join(diff_lines)),
    }
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(yaml.safe_dump(doc, sort_keys=False, width=100))
    print(f"[merge_annotation] -> {out}")
    print("\n".join(diff_lines))
    return 0


if __name__ == "__main__":
    sys.exit(main())
