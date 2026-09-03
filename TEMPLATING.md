# Turning this repo into a cookiecutter template (later)

This repo is deliberately a **plain configure-and-run repo** right now: the
generalisation lives entirely in the config layer (`config.yaml` +
`configs/<system>.yaml` + `config.local.yaml`), and no source file hard-codes
a system name, chain id, organism, or cluster account.

That means a `cookiecutter/` wrapper can be added on top with no code
changes. When we do, `cookiecutter.json` should expose:

| Variable | Feeds | Notes |
|---|---|---|
| `project_slug` | repo dir name | e.g. `abcfold_ifb_<family>_complexes` |
| `first_system_name` | `configs/<name>.yaml` + `config.yaml systems:` | replaces `example_toy` |
| `chains` | `configs/<name>.yaml sequences:` | list of `{type,id,name,sequence}` — prompt via a post-gen hook that reads a FASTA |
| `anchor_chains` / `partner_chains` | same file | |
| `size_class` | same file | `small` \| `large` |
| `ifb_account` | `config.yaml slurm.account` | |
| `ifb_partition` / `exclude_nodes` | `config.yaml slurm.*` | |
| `contact_email` | `config.yaml annotate.interproscan.email` | EBI REST requirement |
| `af3_model_params` | `config.yaml abcfold.model_params` | |
| `chimerax_bin` | `config.yaml chimerax_bin` | host-specific → really belongs in `config.local.yaml` |
| `backends` | `config.yaml abcfold.models` | which of a/b/c/o/p/r |

Post-generation hook responsibilities:
1. drop the carried example configs (`drb2_drb4`, `ago1_*`, …) unless
   `keep_examples=yes`;
2. `git init`;
3. optionally build `envs/controller.yaml`.

Everything else (all of `scripts/`, `workflows/`, `report/`, `tools/`) is
copied verbatim into the generated project.
