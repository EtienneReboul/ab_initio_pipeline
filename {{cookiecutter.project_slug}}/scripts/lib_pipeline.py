#!/usr/bin/env python3
"""
scripts/lib_pipeline.py
========================
Shared helpers for all three stage Snakefiles and their scripts: config
loading (config.yaml + optional config.local.yaml overlay), per-system spec
loading (configs/<system>.yaml), and small chain/sequence utilities.

Kept dependency-free (stdlib + pyyaml) so it can be imported from any of the
per-rule conda envs.
"""
from __future__ import annotations

import copy
from pathlib import Path
from typing import Any

import yaml

REPO_ROOT = Path(__file__).resolve().parent.parent


# ── Config ──────────────────────────────────────────────────────────────────

def _deep_merge(base: dict, over: dict) -> dict:
    out = copy.deepcopy(base)
    for k, v in (over or {}).items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = _deep_merge(out[k], v)
        else:
            out[k] = copy.deepcopy(v)
    return out


def load_config(repo_root: Path | str = REPO_ROOT) -> dict:
    """config.yaml with config.local.yaml deep-merged on top when present."""
    root = Path(repo_root)
    cfg = yaml.safe_load((root / "config.yaml").read_text())
    local = root / "config.local.yaml"
    if local.exists():
        cfg = _deep_merge(cfg, yaml.safe_load(local.read_text()) or {})
    return cfg


# ── Per-system spec ─────────────────────────────────────────────────────────

def load_system(name: str, repo_root: Path | str = REPO_ROOT) -> dict:
    """configs/<name>.yaml, with a few normalisations:
      - anchor_chains always a list (old configs used scalar anchor_chain)
      - partner_chains always a list
      - annotation_reviewed defaults False
      - plip_passes falls back to a single 'main' pass built from a legacy
        top-level `plip:` block if present
    """
    root = Path(repo_root)
    spec = yaml.safe_load((root / "configs" / f"{name}.yaml").read_text())

    if "anchor_chains" not in spec and "anchor_chain" in spec:      # legacy
        spec["anchor_chains"] = [spec["anchor_chain"]]
    spec["anchor_chains"] = list(spec.get("anchor_chains", []))
    spec["partner_chains"] = list(spec.get("partner_chains", []))
    spec.setdefault("annotation_reviewed", False)
    spec.setdefault("size_class", "small")

    if "plip_passes" not in spec:
        passes = []
        if isinstance(spec.get("plip"), dict):
            p = spec["plip"]
            passes.append({"name": "main", "dnareceptor": p.get("dnareceptor", False),
                           "fix_pdb": p.get("fix_pdb", False), "chains": p.get("chains")})
        for legacy_key, nm in (("plip_rna_ligands", "rna_ligs"),
                               ("plip_drb2_drb4", "drb2_drb4")):
            if isinstance(spec.get(legacy_key), dict):
                p = spec[legacy_key]
                passes.append({"name": nm, "dnareceptor": p.get("dnareceptor", False),
                               "fix_pdb": p.get("fix_pdb", False), "chains": p.get("chains")})
        spec["plip_passes"] = passes or [
            {"name": "main", "dnareceptor": False, "fix_pdb": False, "chains": None}
        ]

    if isinstance(spec.get("plip"), dict) and "plip_host" not in spec:
        spec["plip_host"] = {
            "docker_memory": spec["plip"].get("docker_memory", "8g"),
            "docker_platform": spec["plip"].get("docker_platform", ""),
        }
    spec.setdefault("plip_host", {"docker_memory": "8g", "docker_platform": ""})
    return spec


# ── Chain helpers ──────────────────────────────────────────────────────────

_PROTEIN_TYPES = {"protein"}
_NUCLEIC_TYPES = {"rna", "dna"}


def chains_of_type(spec: dict, types: set[str]) -> list[dict]:
    return [s for s in spec.get("sequences", []) if s.get("type") in types]


def protein_chains(spec: dict) -> list[dict]:
    return chains_of_type(spec, _PROTEIN_TYPES)


def nucleic_chains(spec: dict) -> list[dict]:
    return chains_of_type(spec, _NUCLEIC_TYPES)


def norm_id(x) -> str:
    """fold_input.json chain ids may be a list (['A']) for multi-copy chains;
    configs/<system>.yaml uses a bare string. Normalise to the first string."""
    if isinstance(x, (list, tuple)):
        return str(x[0]) if x else ""
    return str(x)


def chain_by_id(spec: dict, cid: str) -> dict | None:
    for s in spec.get("sequences", []):
        if s.get("id") == cid:
            return s
    return None


def token_count(spec: dict) -> int:
    """Rough ABCfold token count: 1 per protein residue + 1 per nucleotide."""
    n = 0
    for s in spec.get("sequences", []):
        if s.get("type") in _PROTEIN_TYPES | _NUCLEIC_TYPES:
            n += len(s.get("sequence", ""))
    return n


# Documented hard token caps observed on real IFB runs (see the DRB2 pipeline
# project notes). Used by abcfold_input_stats.py to flag risk up front.
BACKEND_TOKEN_CAPS = {
    "chai1": 2048,
    "protenix": 2560,
    "boltz": 2600,        # de-facto: OOM on H200 above ~this for a dense complex
}
