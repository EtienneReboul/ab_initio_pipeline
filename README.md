# ab_initio_pipeline — cookiecutter template

This repository **is** a [cookiecutter](https://cookiecutter.readthedocs.io/)
template, not a runnable pipeline by itself. It generates new, ready-to-run
copies of a 3-stage ABCfold ab-initio complex-modelling pipeline
(preprocessing → processing → postprocessing), each pre-configured for a
system you specify at generation time.

`cookiecutter.json` and `hooks/` (this repo's template machinery) live at
the root. Everything the generated pipeline actually needs — `scripts/`,
`workflows/`, `config.yaml`, `configs/`, `report/`, `tools/`, `envs/`, its
own `README.md`, etc. — lives under `{{cookiecutter.project_slug}}/` and is
copied (and partly filled in) into a new project directory each time you
run `cookiecutter` against this repo.

If you just want to **use** an already-generated project, this file isn't
what you want — open that project's own `README.md` instead (its "Quick
start" section covers running the three stages). This file is about
generating that project in the first place.

---

## 1. Prerequisites

You only need this to *generate* a project — the generated project has its
own, separate dependencies (a conda env, ChimeraX, Docker, etc.; see its
`README.md` and `envs/controller.yaml`).

- **Python 3.9+** and **pip**
- **[cookiecutter](https://cookiecutter.readthedocs.io/)** itself
- **git** (the generated project gets `git init`-ed automatically; also
  needed if you generate from the GitHub URL instead of a local clone)

Install cookiecutter — a virtualenv or `pipx` keeps it out of your other
Python environments, but a plain `pip install --user` works fine too:

```bash
python3 -m pip install --user cookiecutter
# or, isolated:
pipx install cookiecutter
```

Verify it's on your `PATH`:

```bash
cookiecutter --version
```

---

## 2. Generating a project

From a local clone of this repo:

```bash
cookiecutter /path/to/ab_initio_pipeline
```

Directly from GitHub, without cloning this repo yourself first:

```bash
cookiecutter https://github.com/EtienneReboul/ab_initio_pipeline.git
# or, with the gh: shorthand:
cookiecutter gh:EtienneReboul/ab_initio_pipeline
```

By default cookiecutter creates the new project as a sibling directory of
wherever you run it from (named after your `project_slug` answer below).
Use `-o <dir>` to put it somewhere else:

```bash
cookiecutter -o ~/projects gh:EtienneReboul/ab_initio_pipeline
```

You'll then be asked a series of questions (next section), and once you
answer the last one, cookiecutter renders the project and a post-generation
script finishes the setup automatically (also below).

---

## 3. The prompts

Asked in this order — later defaults can depend on earlier answers (e.g.
`project_slug`'s default embeds whatever you typed for `first_system_name`).

| Prompt | What it's for | Default |
|---|---|---|
| `first_system_name` | Name of your first system: becomes `config.yaml`'s `systems:` entry and (if you opt into the FASTA flow below) the filename of `configs/<name>.yaml`. Must start with a letter, then letters/digits/underscores only — it's used both as a YAML key and a filename. | `example_toy` |
| `project_slug` | The generated directory's name. Must start with a letter, then letters/digits/`_`/`-`. | `abcfold_ifb_<first_system_name>_complex` |
| `contact_email` | **Required, checked before anything is generated.** EBI's InterProScan REST service (used by stage 1's domain annotation) rejects requests without a real contact address. Never leave this as the placeholder — generation refuses to proceed if you do, or if it doesn't contain an `@`. | `you@example.org` (rejected) |
| `scheduler` | Which cluster scheduler stage 2 targets: `slurm`, `pbs`, or `htcondor`. **Only `slurm` has actually been run and validated** (on the IFB Core Cluster) — `pbs`/`htcondor` follow their Snakemake executor plugins' documented interfaces but are untested against a live cluster. See the generated project's `README.md` → "Non-SLURM schedulers" before picking either. | `slurm` |
| `cluster_account` | Your account/project name on that scheduler (SLURM account, PBS account — not used for htcondor). Every site's value is different; there's no sane default. | `your_cluster_account` |
| `chimerax_bin` | Path to the ChimeraX executable used by stage 3's minimization step (runs locally, not on the cluster). Override again later per-machine in `config.local.yaml` if needed. | `/Applications/ChimeraX.app/Contents/MacOS/ChimeraX` |
| `backends` | Comma-separated subset of the ABCfold backends to run: `alphafold3,boltz,chai1,openfold3,protenix,rosettafold3`. Drop any you don't have set up (e.g. `protenix` is currently broken on IFB — see the generated `README.md`). | all six |
| `setup_first_system_now` | `yes` walks you through building `configs/<first_system_name>.yaml` interactively from a FASTA file right after generation (see below). `no` leaves it for you to write by hand later (or does nothing if you kept `first_system_name: example_toy`, since that config ships pre-written as a worked example). | `no` |

---

## 4. What happens automatically after you answer

Once the last prompt is answered, generation happens in two steps you don't
need to do anything for:

1. **The project is rendered** — `config.yaml`, `config.local.yaml.example`,
   `README.md`, and `workflows/processing/profiles/pbs/config.yaml` get your
   answers substituted in; everything else (all scripts, both Snakefiles,
   `tools/`, `envs/`, `report/`, the other executor profiles) is copied
   byte-for-byte unchanged.
2. **A post-generation script runs** inside the new project and:
   - Drops a few sibling-pipeline leftovers that aren't generic
     (`scripts/symlink_overfolded_samples.py`,
     `scripts/energy_overfolding_filter.py`,
     `workflows/processing/_reference/`, an empty `report/captions/`).
   - Writes `config.local.yaml` (gitignored, never committed) with your
     `contact_email` under `annotate.interproscan.email` — this file has
     never existed until now; every stage loads it on top of `config.yaml`.
   - **If you answered `setup_first_system_now: yes`**, prompts you
     interactively, right there in your terminal, for:
     - a path to a FASTA file with your system's chains,
     - a display name and type (`protein`/`rna`/`dna`) for each chain it
       finds (defaults to the FASTA header's first token / `protein`),
     - which chains are the rigid anchor(s) vs. the partner(s) for pose
       clustering,
     - a `size_class` (`small`/`large`, picks which cluster resource
       profile stage 2 uses),

     then writes `configs/<first_system_name>.yaml` from your answers
     (leaving `domains: {}` and `annotation_reviewed: false` — you curate
     those after stage 1's annotation step, same as any system). Answering
     nothing at the FASTA-path prompt skips this and leaves the file for
     you to write later.
   - Runs `git init` and creates an initial commit, so the generated
     project starts life as its own git repo.

You'll see a short log of each of these as it runs. If git init/commit
fails for some reason (rare — e.g. no git installed), it's reported but
non-fatal; the generated project is still fully usable, you'd just run
`git init` yourself.

---

## 5. Non-interactive generation

Useful for scripting, CI, or just answering everything on one line instead
of being prompted. Pass `--no-input` plus any `key=value` overrides — unset
keys fall back to their `cookiecutter.json` defaults (which means
`contact_email` **must** be passed explicitly, since its default is the
rejected placeholder):

```bash
cookiecutter --no-input gh:EtienneReboul/ab_initio_pipeline \
  contact_email=you@yourinstitute.org \
  first_system_name=my_complex \
  scheduler=slurm \
  cluster_account=my_slurm_account
```

Note `--no-input` also skips the interactive FASTA prompts inside the
post-gen hook even if you pass `setup_first_system_now=yes` — those use
Python's `input()`, not a cookiecutter prompt, so there's nothing to read
from in a non-interactive run. Leave `setup_first_system_now` at its `no`
default for non-interactive generation, or pipe answers to it manually:

```bash
printf '%s\n' /path/to/chains.fasta '' '' '' '' '' '' '' \
  | cookiecutter --no-input gh:EtienneReboul/ab_initio_pipeline \
      contact_email=you@yourinstitute.org setup_first_system_now=yes
```

(blank lines accept each prompt's own default — chain names/types, then
anchor/partner chains, then `size_class`; add more/fewer blanks to match
however many chains your FASTA has).

---

## 6. After generation

`cd` into the new project directory and follow **its own `README.md`**
("Quick start") to actually run the pipeline — install its controller conda
env, review the stage-1 annotation output, and run the three stages in
order. That README also has the full "Running on IFB vs. another SLURM
cluster" / "Non-SLURM schedulers" guidance if you picked `pbs`/`htcondor`
or a non-IFB SLURM site.

To regenerate — e.g. to add a second system, or after changing your answers
— just run `cookiecutter` again; it won't overwrite an existing directory
of the same `project_slug` without asking. For adding *more* systems to an
already-generated project, you don't need to regenerate at all: just add
another `configs/<name>.yaml` (using `configs/_schema.md` / the shipped
`configs/example_toy.yaml` as a template) and add it to `config.yaml`'s
`systems:` list, exactly as that project's own README describes.

---

## 7. Troubleshooting generation

- **"contact_email is required..."** — you left it at the placeholder or
  it has no `@`. Pass a real address.
- **"project_slug '...' is not valid..." / "first_system_name '...' is not
  valid..."** — these become a directory name and a YAML key/filename
  respectively; stick to letters, digits, `_`/`-`, starting with a letter.
- **Generation succeeds but stage 1 fails immediately** — check
  `config.local.yaml` got written (it's gitignored, so `git status` won't
  show it — `cat config.local.yaml` in the generated project) and that
  `annotate.interproscan.email` is your real address.
- **PBS/HTCondor**: expect to need to tune queue names, GPU request syntax,
  and (for PBS) install `snakemake-executor-plugin-pbs` yourself — see the
  generated project's `README.md` → "Non-SLURM schedulers" and the comments
  in its `envs/controller.yaml` / `config.yaml` `hpc:` block for exactly
  what's untested and why.
