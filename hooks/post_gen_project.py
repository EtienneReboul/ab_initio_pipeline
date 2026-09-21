#!/usr/bin/env python3
"""
Post-generation hook: cleans up sibling-pipeline leftovers that don't belong
in a fresh project, writes config.local.yaml (gitignored — must exist before
stage 1 will run, since EBI's InterProScan REST needs a real contact
address), optionally builds configs/<first_system_name>.yaml from a FASTA,
and git-inits the generated project.

Runs from the generated project's root (cookiecutter chdirs here first).
See TEMPLATING.md in the template repo for the design behind this.
"""

import shutil
import string
import subprocess
import sys
from pathlib import Path

FIRST_SYSTEM_NAME = "{{ cookiecutter.first_system_name }}"
CONTACT_EMAIL = "{{ cookiecutter.contact_email }}"
SETUP_FIRST_SYSTEM_NOW = "{{ cookiecutter.setup_first_system_now }}" == "yes"

ROOT = Path.cwd()

# cookiecutter renders hook scripts as Jinja too, so writing two adjacent
# brace characters anywhere below (even inside a Python string or comment)
# breaks the render. These two single-brace variables let the YAML-building
# code produce them without ever having two adjacent in this file's source.
B_OPEN, B_CLOSE = "{", "}"


# ── 1. Drop sibling-pipeline leftovers that aren't generic ──────────────────
def drop_leftovers():
    for path in (
        "scripts/symlink_overfolded_samples.py",
        "scripts/energy_overfolding_filter.py",
        "workflows/processing/_reference",
        "report/captions",
    ):
        p = ROOT / path
        if p.is_dir():
            shutil.rmtree(p)
            print(f"  dropped {path}/ (not generic — see TEMPLATING.md)")
        elif p.is_file():
            p.unlink()
            print(f"  dropped {path} (not generic — see TEMPLATING.md)")


# ── 2. config.local.yaml — required, gitignored, never in config.yaml ───────
def write_config_local():
    content = f'''# Local overrides — gitignored. Loaded on top of config.yaml by every stage.
# Written by cookiecutter's post-gen hook from the contact_email prompt.
annotate:
  interproscan:
    email: "{CONTACT_EMAIL}"

# ── Quick-test knobs (uncomment for a fast smoke test / GPU-scarce run) ─────
# af3:
#   n_replicas: 3
# abcfold:
#   number_of_models: 1
#   models:
#     protenix: false

# ── Host tool paths, SLURM account, etc. ────────────────────────────────────
# See config.local.yaml.example for the full set of overridable keys.
'''
    (ROOT / "config.local.yaml").write_text(content)
    print("  wrote config.local.yaml (annotate.interproscan.email set)")


# ── 3. Optional: build configs/<first_system_name>.yaml from a FASTA ────────
def parse_fasta(path: Path):
    records = []
    header, chunks = None, []
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(chunks)))
            header, chunks = line[1:].strip(), []
        else:
            chunks.append(line)
    if header is not None:
        records.append((header, "".join(chunks)))
    return records


def prompt(msg, default=""):
    suffix = f" [{default}]" if default else ""
    val = input(f"{msg}{suffix}: ").strip()
    return val or default


def build_first_system():
    target = ROOT / "configs" / f"{FIRST_SYSTEM_NAME}.yaml"
    if target.exists():
        ans = prompt(f"configs/{FIRST_SYSTEM_NAME}.yaml already exists — overwrite? (y/N)", "n")
        if ans.lower() not in ("y", "yes"):
            print("  kept existing configs/%s.yaml untouched" % FIRST_SYSTEM_NAME)
            return

    fasta_path = prompt("Path to a FASTA file with this system's chains (blank to skip)")
    if not fasta_path:
        print("  no FASTA given — skipping configs/%s.yaml generation" % FIRST_SYSTEM_NAME)
        print("  (hand-write it later from configs/_schema.md / configs/example_toy.yaml)")
        return
    fasta = Path(fasta_path).expanduser()
    if not fasta.is_file():
        print(f"  '{fasta}' not found — skipping configs/{FIRST_SYSTEM_NAME}.yaml generation")
        return

    records = parse_fasta(fasta)
    if not records:
        print(f"  no records parsed from '{fasta}' — skipping")
        return

    chain_ids = list(string.ascii_uppercase)
    sequences = []
    for i, (header, seq) in enumerate(records):
        if i >= len(chain_ids):
            print("  more than 26 chains in the FASTA — stopping at Z")
            break
        chain_id = chain_ids[i]
        default_name = header.split()[0] if header else f"chain_{chain_id}"
        name = prompt(f"  chain {chain_id} ({len(seq)} aa) name", default_name)
        chain_type = prompt(f"  chain {chain_id} type [protein/rna/dna]", "protein")
        sequences.append({"type": chain_type, "id": chain_id, "name": name, "sequence": seq})

    all_ids = [s["id"] for s in sequences]
    anchor_default = all_ids[0] if all_ids else ""
    partner_default = ",".join(all_ids[1:])
    anchor_chains = [c.strip() for c in prompt("anchor_chains (comma-separated)", anchor_default).split(",") if c.strip()]
    partner_chains = [c.strip() for c in prompt("partner_chains (comma-separated)", partner_default).split(",") if c.strip()]
    size_class = prompt("size_class [small/large]", "small")

    seq_lines = "\n".join(
        '  - {ob}type: {type}, id: {id}, name: "{name}", sequence: "{sequence}"{cb}'.format(
            ob=B_OPEN, cb=B_CLOSE, **s
        )
        for s in sequences
    )
    groups = [anchor_chains, partner_chains] if partner_chains else [anchor_chains]
    chains_str = "[" + ",".join("[" + ",".join(f"'{c}'" for c in g) + "]" for g in groups) + "]"

    content = f'''# {FIRST_SYSTEM_NAME} — generated by cookiecutter's post-gen hook from {fasta.name}.
# Review anchor_chains/partner_chains below, run stage 1, curate the
# `domains:` block from data/annotation/{FIRST_SYSTEM_NAME}/annotation.yaml,
# then set annotation_reviewed: true before stage 2 (see README.md).

name: {FIRST_SYSTEM_NAME}
size_class: {size_class}

sequences:
{seq_lines}

anchor_chains:  [{", ".join(anchor_chains)}]
partner_chains: [{", ".join(partner_chains)}]

annotation_reviewed: false

domains: {B_OPEN}{B_CLOSE}

plip_passes:
  - {B_OPEN}name: main, dnareceptor: false, fix_pdb: false, chains: "{chains_str}"{B_CLOSE}

plip_host:
  docker_memory: "8g"
  docker_platform: "linux/amd64"
'''
    target.write_text(content)
    print(f"  wrote configs/{FIRST_SYSTEM_NAME}.yaml ({len(sequences)} chain(s))")


# ── 4. git init ───────────────────────────────────────────────────────────
def git_init():
    try:
        subprocess.run(["git", "init", "-q"], cwd=ROOT, check=True)
        subprocess.run(["git", "add", "-A"], cwd=ROOT, check=True)
        subprocess.run(
            [
                "git", "-c", "user.name=abcfold-pipeline-cookiecutter",
                "-c", f"user.email={CONTACT_EMAIL}",
                "commit", "-q", "-m", "Initial commit generated by the ab_initio_pipeline cookiecutter template",
            ],
            cwd=ROOT, check=True,
        )
        print("  git-initialized the project with an initial commit")
    except Exception as exc:  # not fatal — the project is usable either way
        print(f"  NOTE: git init/commit did not complete cleanly ({exc}); run it yourself.")


def main():
    print("Post-generation setup:")
    drop_leftovers()
    write_config_local()
    if SETUP_FIRST_SYSTEM_NOW:
        build_first_system()
    else:
        print(
            f"  systems: still points at '{FIRST_SYSTEM_NAME}' — make sure "
            f"configs/{FIRST_SYSTEM_NAME}.yaml exists (configs/example_toy.yaml "
            "ships as a worked example/smoke test) before running stage 1."
        )
    git_init()
    print("\nDone. See README.md 'Quick start' for the next steps.")


if __name__ == "__main__":
    sys.exit(main())
