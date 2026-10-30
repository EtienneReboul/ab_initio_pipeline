#!/usr/bin/env python3
"""
scripts/annotate_domains.py — Stage 1c
======================================
Per protein chain of a system, dump every domain / family / motif hit to a
flat TSV, plus InterPro's feature tracks to a JSON sidecar:

  * ScanProsite  (PROSITE patterns + profiles)  — prosite.expasy.org REST
  * InterPro     (Pfam / SMART / PROSITE profiles / CATH-Gene3D / CDD /
                  PANTHER ... AND the InterPro entries themselves) —
      - fast path: the precomputed matches of a UniProt entry from the
        InterPro API (one request, no email, carries InterPro's
        `representative` domain flag and each signature's InterPro entry),
        used for every chain that declares `uniprot:` in configs/<system>.yaml
        and whose sequence IS that entry (or, with `region: [start, end]`,
        that slice of it: hits are then clipped/shifted onto the chain);
      - fallback: InterProScan 5 on the sequence itself (EBI REST, needs
        annotate.interproscan.email); no representative flags — merge_
        annotation.py recomputes InterProScan's own representative set.

Ported from NanoPulldown (nanopulldown.annotate.domains), adapted to
configs/<system>.yaml chains.

Outputs:
  domains.raw.tsv   chain source db accession name start end score description
                    type ipr ipr_name representative
  features.json     {chain: {disorder|coils|signal_peptide|tm: [[s, e], ...]}}
                    (InterPro extra features: MobiDB-lite, COILS, SignalP/
                    Phobius, TMbed/Phobius — fast path only)

merge_annotation.py turns these into the `domains:` map stage 3 uses (no
human review needed; a hand-written `domains:` in configs/<system>.yaml
still overrides it).
"""
from __future__ import annotations

import argparse
import json
import sys
import time
import xml.etree.ElementTree as ET
from pathlib import Path

import requests

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib_pipeline import load_config, load_system, protein_chains  # noqa: E402

CFG = load_config()
UA = {"User-Agent": "ab-initio-pipeline/annotate_domains"}
COLUMNS = ["chain", "source", "db", "accession", "name", "start", "end", "score", "description",
           "type", "ipr", "ipr_name", "representative"]


def _get(url: str, tries: int = 4, **kw) -> requests.Response:
    """GET with a linear backoff (InterPro / ExPASy answer 5xx / 408 under load)."""
    headers = {**UA, **kw.pop("headers", {})}
    for i in range(tries):
        try:
            r = requests.get(url, headers=headers, timeout=180, **kw)
            if r.status_code in (408, 429, 500, 502, 503, 504):
                raise requests.HTTPError(f"{r.status_code} from {url}")
            return r
        except requests.RequestException as e:
            if i == tries - 1:
                raise
            print(f"  retry {i + 1}/{tries - 1} in {15 * (i + 1)} s: {e}", file=sys.stderr)
            time.sleep(15 * (i + 1))
    raise AssertionError


# ── ScanProsite ────────────────────────────────────────────────────────────

def scanprosite(seq: str, url: str, skip_profiles: bool) -> list[dict]:
    """ScanProsite REST — GET with output=xml, parse <match> elements
    (namespace-agnostic). Layout confirmed live 2026-08-31:
      <match><sequence_ac>..</sequence_ac><start>..</start><stop>..</stop>
             <signature_ac>PS#####</signature_ac><score>..</score>...</match>
    """
    params = {"seq": seq, "output": "xml"}
    if skip_profiles:
        params["skip"] = "on"          # patterns only
    r = _get(url, params=params)
    r.raise_for_status()
    root = ET.fromstring(r.text)

    def _txt(el, tag):
        for c in el.iter():
            if c.tag.rsplit("}", 1)[-1] == tag and c.text:
                return c.text.strip()
        return None

    hits = []
    for m in root.iter():
        if m.tag.rsplit("}", 1)[-1] != "match":
            continue
        start, stop = _txt(m, "start"), _txt(m, "stop")
        sig = _txt(m, "signature_ac")
        if not (start and stop and sig):
            continue
        hits.append({
            "source": "ScanProsite", "db": "PROSITE",
            "accession": sig, "name": sig,
            "start": int(start), "end": int(stop),
            "score": _txt(m, "score") or "",
            "description": "",
        })
    return hits


# ── InterPro API by UniProt accession (fast path) ──────────────────────────

# InterPro API source_database -> the db labels InterProScan's TSV uses, so
# merge_annotation.py sees the same values either way.
_IPR_DB = {"interpro": "InterPro", "pfam": "Pfam", "smart": "SMART", "profile": "PROSITE profiles",
           "prosite": "PROSITE patterns", "cathgene3d": "CATH-Gene3D", "cath-gene3d": "CATH-Gene3D",
           "ssf": "SUPERFAMILY", "cdd": "CDD", "panther": "PANTHER", "prints": "PRINTS",
           "pirsf": "PIRSF", "hamap": "HAMAP", "ncbifam": "NCBIfam", "sfld": "SFLD",
           "mobidblt": "MobiDB-lite", "mobidb-lite": "MobiDB-lite", "antifam": "AntiFam",
           "funfam": "CATH-FunFam", "pirsr": "PIRSR"}
_SITE_TYPES = {"conserved_site", "active_site", "binding_site", "ptm", "site"}


def _norm_type(t: str | None) -> str | None:
    t = (t or "").lower() or None
    return "site" if t in _SITE_TYPES else t


def uniprot_sequence(acc: str) -> str:
    r = _get(f"https://rest.uniprot.org/uniprotkb/{acc}.fasta")
    r.raise_for_status()
    return "".join(r.text.splitlines()[1:])


def interpro_by_accession(acc: str, seq_len: int, api: str) -> list[dict] | None:
    """Every member-signature hit AND InterPro entry (db "InterPro") of a
    UniProt entry, with InterPro's `representative` flag. None when the API
    has nothing for it (204/404) or computed on a different length."""
    url = f"{api}/entry/all/protein/uniprot/{acc}?page_size=200"
    hits: list[dict] = []
    while url:
        r = _get(url, headers={"Accept": "application/json"})
        if r.status_code in (204, 404):
            return None
        r.raise_for_status()
        data = r.json()
        for e in data.get("results", []):
            meta = e["metadata"]
            for prot in e.get("proteins", []):
                if int(prot.get("protein_length", seq_len)) != seq_len:
                    return None
                for loc in prot.get("entry_protein_locations", []):
                    frs = loc.get("fragments", [])
                    if not frs:
                        continue
                    hits.append({
                        "source": "InterPro",
                        "db": _IPR_DB.get(meta["source_database"].lower(), meta["source_database"]),
                        "accession": meta["accession"], "name": meta.get("name") or meta["accession"],
                        # discontinuous domains: one span, first to last fragment
                        "start": int(frs[0]["start"]), "end": int(frs[-1]["end"]),
                        "score": "" if loc.get("score") is None else loc["score"],
                        "description": meta.get("name") or "",
                        "type": _norm_type(meta.get("type")),
                        "ipr": meta.get("integrated"), "ipr_name": None,
                        "representative": bool(loc.get("representative")),
                    })
        url = data.get("next")
    names = {h["accession"]: h["name"] for h in hits if h["db"] == "InterPro"}
    for h in hits:
        if h["ipr"]:
            h["ipr_name"] = names.get(h["ipr"])
    return hits


# InterPro extra-feature source -> track; first source listed wins per track
_FEATURE_SOURCES = {
    "disorder": [("mobidblt", None)],
    "coils": [("coils", None)],
    "signal_peptide": [("signalp_e", None), ("signalp_g+", None), ("signalp_g-", None),
                       ("phobius", "SIGNAL_PEPTIDE")],
    "tm": [("tmbed", "TMhelix"), ("phobius", "TRANSMEMBRANE")],
}


def interpro_features(acc: str, api: str) -> dict[str, list[list[int]]]:
    r = _get(f"{api}/protein/uniprot/{acc}?extra_features", headers={"Accept": "application/json"})
    if r.status_code in (204, 404):
        return {}
    r.raise_for_status()
    found: dict[tuple[str, str], list[list[int]]] = {}
    for key, v in (r.json() or {}).items():
        db = (v.get("source_database") or "").lower()
        for loc in v.get("locations", []):
            for fr in loc.get("fragments", []):
                found.setdefault((db, key), []).append([int(fr["start"]), int(fr["end"])])
    out = {}
    for track, sources in _FEATURE_SOURCES.items():
        for db, prefix in sources:
            locs = [l for (d, k), ls in found.items() if d == db and (prefix is None or k.startswith(prefix))
                    for l in ls]
            if locs:
                out[track] = sorted(map(list, {tuple(l) for l in locs}))
                break
    return out


def clip_span(start: int, end: int, region: tuple[int, int]) -> tuple[int, int] | None:
    """Full-length span -> chain numbering of `region` (1-based, inclusive)."""
    a, b = region
    s, e = max(start, a), min(end, b)
    return (s - a + 1, e - a + 1) if s <= e else None


def interpro_record(acc: str, api: str, cache_dir: Path) -> dict | None:
    """{"sequence", "hits", "features"} of a UniProt entry (full-length
    numbering), cached on disk per accession, so a rerun never depends on the
    InterPro API being up. The feature tracks are best-effort: when only that
    request fails the hits are still cached, and the features retried next run."""
    cache = cache_dir / f"{acc}.json"
    rec = json.loads(cache.read_text()) if cache.exists() else None
    if rec is None:
        full = uniprot_sequence(acc)
        hits = interpro_by_accession(acc, len(full), api)
        if hits is None:
            return None
        rec = {"sequence": full, "hits": hits, "features": None}
    if rec["features"] is None:
        try:
            rec["features"] = interpro_features(acc, api)
        except Exception as e:                      # noqa: BLE001
            print(f"  InterPro features for {acc} FAILED ({e}) — retried next run", file=sys.stderr)
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache.write_text(json.dumps(rec))
    return rec


def interpro_for_chain(ch: dict, api: str, cache_dir: Path) -> tuple[list[dict], dict] | None:
    """(hits, features) of a chain from its `uniprot:` entry, in chain numbering,
    or None when the chain declares none / isn't that entry (or slice of it)."""
    acc = ch.get("uniprot")
    if not acc:
        return None
    rec = interpro_record(acc, api, cache_dir)
    if rec is None:
        return None
    full = rec["sequence"]
    region = ch.get("region")
    if not region:              # no `region:` -> the slice of the entry the chain is, if unambiguous
        i = full.find(ch["sequence"])
        unique = i >= 0 and full.find(ch["sequence"], i + 1) < 0
        region = (i + 1, i + len(ch["sequence"])) if unique else (1, len(full))
    region = tuple(region)
    if full[region[0] - 1:region[1]] != ch["sequence"]:
        print(f"  chain {ch['id']}: sequence is not {acc}[{region[0]}-{region[1]}] — "
              "InterProScan fallback", file=sys.stderr)
        return None
    hits, feats = rec["hits"], rec["features"] or {}
    hits = [{**h, "start": s, "end": e} for h in hits
            for s, e in [clip_span(h["start"], h["end"], region) or (0, 0)] if s]
    feats = {k: [[s, e] for a, b in v for s, e in [clip_span(a, b, region) or (0, 0)] if s]
             for k, v in feats.items()}
    return hits, {k: v for k, v in feats.items() if v}


# ── InterProScan 5 on the sequence (EBI REST, fallback) ────────────────────

def _ebi_run(seq: str, cfg: dict) -> str:
    base = cfg["base_url"]
    payload = {"email": cfg["email"], "sequence": seq, "stype": "p"}
    appl = ",".join(cfg.get("applications", []))
    if appl:
        payload["appl"] = appl
    r = requests.post(f"{base}/run", data=payload, headers=UA, timeout=120)
    if r.status_code == 400 and "appl" in r.text and appl:
        # bad application name(s) — fall back to InterProScan's default set
        print(f"  [interproscan] server rejected appl='{appl}', retrying with defaults",
              file=sys.stderr)
        payload.pop("appl")
        r = requests.post(f"{base}/run", data=payload, headers=UA, timeout=120)
    r.raise_for_status()
    job_id = r.text.strip()
    poll = int(cfg.get("poll_s", 30))
    while True:
        s = requests.get(f"{base}/status/{job_id}", headers=UA, timeout=60).text.strip()
        if s == "FINISHED":
            break
        if s in ("ERROR", "FAILURE", "NOT_FOUND"):
            raise RuntimeError(f"InterProScan job {job_id} status={s}")
        time.sleep(poll)
    out = requests.get(f"{base}/result/{job_id}/tsv", headers=UA, timeout=120)
    out.raise_for_status()
    return out.text


# InterProScan analyses that are feature tracks, not domain signatures
_IPS_FEATURES = {"mobidblite": "disorder", "mobidb-lite": "disorder", "coils": "coils",
                 "signalp_euk": "signal_peptide", "signalp": "signal_peptide",
                 "tmhmm": "tm", "phobius": "tm"}


def interproscan(seq: str, cfg: dict) -> tuple[list[dict], dict]:
    hits, feats = [], {}
    for line in _ebi_run(seq, cfg).splitlines():
        p = line.rstrip("\n").split("\t")
        # IPRScan5 tsv: 0 md5 1 md5hex 2 len 3 analysis 4 sig_acc 5 sig_desc
        # 6 start 7 stop 8 score 9 status 10 date 11 ipr_acc 12 ipr_desc ...
        if len(p) < 8:
            continue
        try:
            start, stop = int(p[6]), int(p[7])
        except ValueError:
            continue
        feat = _IPS_FEATURES.get(p[3].lower())
        if feat:
            feats.setdefault(feat, []).append([start, stop])
            continue
        ipr = p[11] if len(p) > 11 and p[11] not in ("", "-") else None
        hits.append({
            "source": "InterProScan", "db": p[3], "accession": p[4],
            "name": next((x for x in (p[5] if len(p) > 5 else "", p[4]) if x not in ("", "-")), p[4]),
            "start": start, "end": stop,
            "score": p[8] if len(p) > 8 else "",
            "description": next((x for x in (p[12] if len(p) > 12 else "", p[5] if len(p) > 5 else "")
                                 if x not in ("", "-")), ""),
            "type": None, "ipr": ipr,
            "ipr_name": p[12] if ipr and len(p) > 12 and p[12] not in ("", "-") else None,
            "representative": None,
        })
    return hits, feats


# ── main ────────────────────────────────────────────────────────────────

def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--system", required=True)
    ap.add_argument("--spec", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--features", required=True, help="features.json sidecar")
    args = ap.parse_args()

    spec = load_system(args.system, Path(args.spec).resolve().parent.parent)
    acfg = CFG["annotate"]
    api = acfg["interproscan"].get("api_url", "https://www.ebi.ac.uk/interpro/api")
    email = (acfg["interproscan"].get("email") or "").strip()
    cache_dir = Path(CFG["dirs"]["annotation"]) / "_cache" / "interpro"
    rows: list[dict] = []
    features: dict[str, dict] = {}

    # identical chains (homodimers etc.) — query once per (sequence, entry, region)
    groups: dict[tuple, list[dict]] = {}
    for ch in protein_chains(spec):
        key = (ch["sequence"], ch.get("uniprot"), tuple(ch.get("region") or ()))
        groups.setdefault(key, []).append(ch)

    for (seq, acc, _), chs in groups.items():
        cids = [c["id"] for c in chs]
        print(f"[annotate_domains] chains {cids} ({len(seq)} aa, uniprot={acc or '-'})", flush=True)
        try:
            ps = scanprosite(seq, acfg["prosite"]["scan_url"], acfg["prosite"].get("skip_profiles", False))
            print(f"  ScanProsite: {len(ps)} hit(s)")
        except Exception as e:                       # noqa: BLE001
            print(f"  ScanProsite FAILED: {e}", file=sys.stderr)
            ps = []
        ip, feats = [], {}
        try:
            got = interpro_for_chain(chs[0], api, cache_dir)
        except Exception as e:                       # noqa: BLE001
            print(f"  InterPro API FAILED: {e}", file=sys.stderr)
            got = None
        if got is not None:
            ip, feats = got
            print(f"  InterPro API ({acc}): {len(ip)} hit(s), features {sorted(feats)}")
        elif email:
            try:
                ip, feats = interproscan(seq, acfg["interproscan"])
                print(f"  InterProScan: {len(ip)} hit(s)")
            except Exception as e:                   # noqa: BLE001
                print(f"  InterProScan FAILED: {e}", file=sys.stderr)
        else:
            print("  no `uniprot:` match and no annotate.interproscan.email — ScanProsite only",
                  file=sys.stderr)
        for cid in cids:
            features[cid] = feats
            rows += [{**h, "chain": cid} for h in ps + ip]

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w") as fh:
        fh.write("\t".join(COLUMNS) + "\n")
        for r in sorted(rows, key=lambda x: (x["chain"], x["start"], x["end"])):
            fh.write("\t".join("" if r.get(c) is None else str(r[c]).replace("\t", " ")
                               for c in COLUMNS) + "\n")
    Path(args.features).write_text(json.dumps(features, indent=1) + "\n")
    print(f"[annotate_domains] wrote {len(rows)} row(s) -> {out}, features -> {args.features}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
