#!/usr/bin/env python3
"""
scripts/minimize_openmm.py — Stage 3d
=====================================
OpenMM minimization of one selected model (replaces the ChimeraX minimize
script). Ported from NanoPulldown's nanopulldown.structure.minimize, extended
to RNA/DNA and ligand chains.

  proteins, RNA, DNA   amber14 (ff14SB / OL3 / OL15) + GBn2 implicit solvent
  monatomic ions       amber14 ion templates
  free ligands         GAFF2 + AM1-BCC from scripts/parametrize_ligands.py
  covalent / multi-CCD ligands   dropped (no template), as before

Sanitize (5'-phosphates and waters off) -> PDBFixer (missing heavy atoms,
hydrogens at pH 7) -> vacuum, charge-free pre-relaxation that pulls overlaps
apart -> GBn2 minimization in chunks with the divergence guards -> <out>.pdb +
<out stem>_energy.csv (step, energy_kJ_mol).

Heavy clashes (minimize.cg2all_fallback): a model whose first evaluation is
infinite or NaN (atoms ~0.01 A apart) gets its PROTEIN chains rebuilt by cg2all
(Heo & Feig; CA trace in, all-atom out, CA kept) and is minimized again.
RNA/DNA/ligand chains keep their model coordinates. Nothing else is retried.

There is no implicit CPU fallback: the first platform of
minimize.platform_preference that starts is used, else the job fails.
Chain ids and residue numbers are kept (PLIP passes rely on them).
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import subprocess
import sys
import tempfile
import time
from contextlib import contextmanager
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib_pipeline import ligand_kind, load_config, load_system, norm_id  # noqa: E402

NUCLEIC = {"A", "C", "G", "U", "DA", "DC", "DG", "DT"}
WATER = {"HOH", "WAT", "H2O", "DOD"}
PHOSPHATE_5P = {"P", "OP1", "OP2", "OP3", "O1P", "O2P", "O3P"}


def log(msg: str) -> None:
    print(f"[minimize] {msg}", flush=True)


class NonFiniteStart(RuntimeError):
    """Infinite/NaN energy or force at the start: the cg2all-fallback case."""


@contextmanager
def finite_guard():
    """Turn OpenMM's non-finite failures into NonFiniteStart (PDBFixer and
    addHydrogens minimize internally, so they can hit them too)."""
    try:
        yield
    except ZeroDivisionError as e:                         # exactly coincident atoms
        raise NonFiniteStart("coincident atoms") from e
    except Exception as e:                                 # noqa: BLE001
        if "infinite or NaN" in str(e) or "Particle coordinate is NaN" in str(e):
            raise NonFiniteStart(str(e).splitlines()[0]) from e
        raise


# ── Platform ────────────────────────────────────────────────────────────────

def _opencl_gpu_platform(platform) -> int | None:
    """Index of the first OpenCL platform of a GPU vendor (never a CPU one)."""
    from openmm import Context, System, VerletIntegrator

    probe = System()
    probe.addParticle(1.0)
    for i in range(8):
        try:
            ctx = Context(probe, VerletIntegrator(0.001), platform, {"OpenCLPlatformIndex": str(i)})
        except Exception:                                  # noqa: BLE001
            break
        vendor = platform.getPropertyValue(ctx, "OpenCLPlatformName")
        del ctx
        if any(v in vendor for v in ("NVIDIA", "AMD Accelerated", "Apple")):
            return i
    return None


def make_simulation(topology, system, preference: list[str]):
    from openmm import LangevinMiddleIntegrator, Platform, unit
    from openmm.app import Simulation

    for name in preference:
        try:
            platform = Platform.getPlatformByName(name)
            props = {}
            if name == "OpenCL":
                idx = _opencl_gpu_platform(platform)
                if idx is None:
                    raise RuntimeError("no GPU OpenCL platform")
                props = {"OpenCLPlatformIndex": str(idx)}
            integrator = LangevinMiddleIntegrator(300 * unit.kelvin, 1 / unit.picosecond,
                                                  0.002 * unit.picoseconds)
            sim = Simulation(topology, system, integrator, platform, props)
        except Exception as e:                             # noqa: BLE001
            log(f"platform {name} unusable ({str(e).splitlines()[0][:120]}) -- trying the next")
            continue
        log(f"platform: {name}")
        return sim
    raise RuntimeError(f"no OpenMM platform could start (tried {', '.join(preference)}); "
                       "add CPU to minimize.platform_preference to allow a (slow) CPU run")


# ── Ligands ─────────────────────────────────────────────────────────────────

def _pose_ligand(chain_id, atoms, sdf: Path):
    """RDKit molecule with hydrogens: template chemistry on the model's heavy
    atom coordinates. `atoms`: [(name, element, (x, y, z) in A)]. Mapped by atom
    name when the CCD names all match, else by perceived connectivity."""
    from rdkit import Chem
    from rdkit.Chem import rdDetermineBonds
    from rdkit.Geometry import Point3D

    tmpl = next(iter(Chem.SDMolSupplier(str(sdf))))
    names = tmpl.GetProp("atom_names").split(",")
    by_name = {n: p for n, _, p in atoms}
    if len(atoms) != tmpl.GetNumAtoms():
        raise RuntimeError(f"ligand {chain_id}: model has {len(atoms)} heavy atoms, "
                           f"its template has {tmpl.GetNumAtoms()}")
    if len(by_name) == len(atoms) and set(names) == set(by_name):
        mol, coords, out_names = Chem.Mol(tmpl), [by_name[n] for n in names], names
    else:
        rw = Chem.RWMol()
        conf = Chem.Conformer(len(atoms))
        for i, (_, el, p) in enumerate(atoms):
            rw.AddAtom(Chem.Atom(el))
            conf.SetAtomPosition(i, Point3D(*p))
        rw.AddConformer(conf)
        raw = rw.GetMol()
        rdDetermineBonds.DetermineConnectivity(raw)
        # match on elements + connectivity only (no bond orders, no charges)
        graph = Chem.RWMol()
        for a in tmpl.GetAtoms():
            graph.AddAtom(Chem.Atom(a.GetAtomicNum()))
        for b in tmpl.GetBonds():
            graph.AddBond(b.GetBeginAtomIdx(), b.GetEndAtomIdx(), Chem.BondType.SINGLE)
        match = raw.GetSubstructMatch(graph.GetMol())
        if not match:
            raise RuntimeError(f"ligand {chain_id}: geometry does not match its template "
                               "(clashed or broken model?)")
        mol = Chem.Mol(tmpl)
        coords = [atoms[k][2] for k in match]
        out_names = [atoms[k][0] for k in match]
    mol.RemoveAllConformers()
    conf = Chem.Conformer(mol.GetNumAtoms())
    for i, p in enumerate(coords):
        conf.SetAtomPosition(i, Point3D(*p))
    mol.AddConformer(conf, assignId=True)
    return Chem.AddHs(mol, addCoords=True), out_names


def _ligand_topology(chain_id, res_id, res_name, mol, names):
    from openmm import Vec3, unit
    from openmm.app import Element, Topology

    top = Topology()
    res = top.addResidue(res_name, top.addChain(chain_id), res_id)
    added, n_h = [], 0
    for a in mol.GetAtoms():
        if a.GetAtomicNum() == 1:
            n_h += 1
            name = f"H{n_h}"
        else:
            name = names[a.GetIdx()]
        added.append(top.addAtom(name, Element.getByAtomicNumber(a.GetAtomicNum()), res))
    for b in mol.GetBonds():
        top.addBond(added[b.GetBeginAtomIdx()], added[b.GetEndAtomIdx()])
    pos = [Vec3(*p) * 0.1 for p in mol.GetConformer().GetPositions().tolist()]
    return top, pos * unit.nanometer


# ── cg2all rebuild ──────────────────────────────────────────────────────────

def cg2all_rebuild(modeller, fb: dict, bindir: Path):
    """Heavy-atom Modeller whose protein chains are rebuilt by cg2all (all-atom
    -> CA trace -> all-atom, CA fixed); other chains are unchanged."""
    from openmm.app import Modeller, PDBFile
    from pdbfixer import PDBFixer

    def is_protein(chain):
        return any(a.name == "CA" and a.element.symbol == "C"
                   for r in chain.residues() for a in r.atoms())

    chains = list(modeller.topology.chains())
    prot = Modeller(modeller.topology, modeller.positions)
    prot.delete([c for c in chains if not is_protein(c)])
    rest = Modeller(modeller.topology, modeller.positions)
    rest.delete([c for c in chains if is_protein(c)])
    if not prot.topology.getNumResidues():
        raise RuntimeError("no protein chain for cg2all to rebuild")
    model = fb["model"]
    with tempfile.TemporaryDirectory() as td:
        src, cg, out, heavy = (Path(td) / f for f in ("in.pdb", "cg.pdb", "all.pdb", "heavy.pdb"))
        with open(src, "w") as fh:
            PDBFile.writeFile(prot.topology, prot.positions, fh, keepIds=True)
        env = {**os.environ, "DGLBACKEND": "pytorch",
               "PATH": f"{bindir}{os.pathsep}{os.environ.get('PATH', '')}"}
        for cmd, made in (
            (["convert_all2cg", "-p", src, "-o", cg, "--cg", model], cg),
            (["convert_cg2all", "-p", cg, "-o", out, "--cg", model, "--fix",
              "--standard-name", "--device", fb["device"]], out),
        ):
            p = subprocess.run([str(bindir / cmd[0]), *map(str, cmd[1:])], env=env,
                               capture_output=True, text=True, timeout=int(fb["timeout_s"]))
            if p.returncode != 0 or not made.exists():
                raise RuntimeError(f"{cmd[0]} failed ({p.returncode}): {p.stderr[-600:]}")
        heavy.write_text("".join(
            ln for ln in out.read_text().splitlines(keepends=True)
            if not (ln.startswith(("ATOM", "HETATM"))
                    and (ln[76:78].strip() == "H" or ln[12:16].strip().lstrip("0123456789")[:1] == "H"))))
        fixer = PDBFixer(filename=str(heavy))
    old, new = list(prot.topology.residues()), list(fixer.topology.residues())
    if [r.name for r in old] != [r.name for r in new]:
        raise RuntimeError("cg2all changed the residue sequence")
    for o, n in zip(old, new):                             # its own numbering is not ours
        n.id, n.insertionCode, n.chain.id = o.id, o.insertionCode, o.chain.id
    if [c.id for c in prot.topology.chains()] != [c.id for c in fixer.topology.chains()]:
        raise RuntimeError("cg2all changed the chain layout")
    done = Modeller(fixer.topology, fixer.positions)
    if rest.topology.getNumAtoms():
        done.add(rest.topology, rest.positions)
    return done


# ── Structure preparation ───────────────────────────────────────────────────

def prepare(cif: Path, spec: dict, lig_dir: Path | None, rebuild=None):
    """(modeller, extra force-field files, residue-template names)"""
    from openmm import unit
    from openmm.app import Modeller, Topology
    from pdbfixer import PDBFixer

    kinds = {norm_id(s["id"]): ligand_kind(s, spec)
             for s in spec["sequences"] if s.get("type") == "ligand"}
    index = json.loads((lig_dir / "ligands.json").read_text()) if lig_dir else {}

    fixer = PDBFixer(filename=str(cif))
    present = [c.id for c in fixer.topology.chains()]
    missing = {norm_id(s["id"]) for s in spec["sequences"]} - set(present)
    if missing:
        raise RuntimeError(f"chain(s) {sorted(missing)} of the config absent from the model "
                           f"(model chains: {present})")

    xyz = fixer.positions.value_in_unit(unit.angstrom)
    ligands, delete = [], []
    for chain in fixer.topology.chains():
        kind = kinds.get(chain.id)
        if kind is not None:
            delete.append(chain)
            if kind == "drop":
                log(f"dropping ligand chain {chain.id} (covalent or multi-residue: no template)")
                continue
            for res in chain.residues():
                heavy = [(a.name, a.element.symbol, tuple(xyz[a.index]))
                         for a in res.atoms() if a.element.symbol != "H"]
                ligands.append((chain.id, kind, res.id, res.name, heavy))
            continue
        for k, res in enumerate(chain.residues()):
            if res.name in WATER:
                delete.append(res)
            elif k == 0 and res.name in NUCLEIC:           # amber's 5' end is a 5'-OH
                delete += [a for a in res.atoms() if a.name in PHOSPHATE_5P]

    mod = Modeller(fixer.topology, fixer.positions)
    mod.delete(delete)
    if rebuild:                                            # before anything that minimizes
        mod = rebuild(mod)
    fixer.topology, fixer.positions = mod.topology, mod.positions
    fixer.missingResidues = {}                             # models are complete
    fixer.findMissingAtoms()
    for res in list(fixer.missingAtoms):                   # keep the 5' end bare
        if res.name in NUCLEIC and res is next(res.chain.residues()):
            fixer.missingAtoms[res] = [a for a in fixer.missingAtoms[res] if a.name not in PHOSPHATE_5P]
    for res, names in list(fixer.missingTerminals.items()):
        if res.name in NUCLEIC:
            fixer.missingTerminals[res] = [n for n in names if n not in PHOSPHATE_5P]
    with finite_guard():
        fixer.addMissingAtoms()
    modeller = Modeller(fixer.topology, fixer.positions)
    with finite_guard():
        modeller.addHydrogens(pH=7.0)                      # termini / 5'-OH included

    extra, ion_chains = [], {}
    for chain_id, kind, res_id, res_name, heavy in ligands:
        if kind == "ion":
            from openmm.app import Element

            top = Topology()
            res = top.addResidue(res_name, top.addChain(chain_id), res_id)
            for name, el, _ in heavy:
                top.addAtom(name, Element.getBySymbol(el), res)
            modeller.add(top, [p for *_, p in heavy] * unit.angstrom)
            ion_chains[chain_id] = res_name
            continue
        if chain_id not in index:
            raise RuntimeError(f"ligand chain {chain_id} has no parameters in {lig_dir}")
        entry = index[chain_id]
        mol, names = _pose_ligand(chain_id, heavy, lig_dir / entry["mol"])
        # 3 columns in a PDB; backends invent names like LIG_E for SMILES ligands
        rname = res_name if len(res_name) <= 3 else "LIG"
        modeller.add(*_ligand_topology(chain_id, res_id, rname, mol, names))
        for f in (lig_dir / "gaff.xml", lig_dir / entry["params"]):
            if str(f) not in extra:
                extra.append(str(f))
        log(f"ligand chain {chain_id} ({rname}): GAFF2, {mol.GetNumAtoms()} atoms")
    templates = {r: ion_chains[r.chain.id] for r in modeller.topology.residues()
                 if r.chain.id in ion_chains}
    return modeller, extra, templates


# ── Minimization ────────────────────────────────────────────────────────────

def ca_volumes(topology, positions) -> dict:
    """Signed CA(N, C, CB) volume per amino acid: > 0 L, < 0 D."""
    from openmm import unit

    xyz, out = positions.value_in_unit(unit.angstrom), {}
    for res in topology.residues():
        at = {a.name: xyz[a.index] for a in res.atoms() if a.name in ("CA", "N", "C", "CB")}
        if len(at) == 4:
            n, c, cb = (at[k] - at["CA"] for k in ("N", "C", "CB"))
            out[f"{res.chain.id}:{res.name}{res.id}"] = (
                n[0] * (c[1] * cb[2] - c[2] * cb[1]) + n[1] * (c[2] * cb[0] - c[0] * cb[2])
                + n[2] * (c[0] * cb[1] - c[1] * cb[0]))
    return out


def minimize(cif: Path, pdb: Path, mz: dict, spec: dict, lig_dir: Path | None, rebuild=None) -> None:
    from openmm import LocalEnergyMinimizer, NonbondedForce, unit
    from openmm.app import CutoffNonPeriodic, ForceField, PDBFile

    log(f"loading {cif.name}")
    modeller, extra, templates = prepare(cif, spec, lig_dir, rebuild)
    top = modeller.topology
    has_ions = ["amber14/tip3p.xml"] if templates else []

    def step(context, n) -> None:
        with finite_guard():
            LocalEnergyMinimizer.minimize(context, maxIterations=n)

    def energy_kj(sim) -> float:
        return sim.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(
            unit.kilojoule_per_mole)

    system = ForceField(*mz["forcefield"], *has_ions, *extra).createSystem(
        top, residueTemplates=templates)
    sim = make_simulation(top, system, mz["platform_preference"])
    sim.context.setPositions(modeller.positions)

    e0 = energy_kj(sim)
    trace = [(0, e0)] if math.isfinite(e0) else []
    log(f"{top.getNumAtoms()} atoms, initial energy: {e0:.1f} kJ/mol")

    interval = max(1, int(mz["trace_interval"]))
    max_steps, probe_steps = int(mz["max_steps"]), int(mz["probe_steps"])
    max_abs, delta = float(mz["max_abs_energy"]), float(mz["convergence_delta"])
    budget_s = 60.0 * float(mz.get("max_minutes") or 0)
    t0 = time.monotonic()

    # Pre-relaxation: raw models start at ~1e12 kJ/mol or worse, where OpenMM's
    # minimizer falls back to a very slow double-precision pass. Charges are off:
    # overlapping opposite charges would collapse onto each other in vacuum.
    pre, pre_done, e_start = mz.get("pre_relax") or {}, 0, e0
    if pre.get("enabled"):
        psys = ForceField(*pre["forcefield"], *has_ions, *extra).createSystem(
            top, nonbondedMethod=CutoffNonPeriodic, residueTemplates=templates,
            nonbondedCutoff=float(pre["cutoff_nm"]) * unit.nanometer)
        for force in psys.getForces():
            if isinstance(force, NonbondedForce):
                for i in range(force.getNumParticles()):
                    _, sigma, eps = force.getParticleParameters(i)
                    force.setParticleParameters(i, 0.0, sigma, eps)
                for k in range(force.getNumExceptions()):
                    p1, p2, _, sigma, eps = force.getExceptionParameters(k)
                    force.setExceptionParameters(k, p1, p2, 0.0, sigma, eps)
        psim = make_simulation(top, psys, mz["platform_preference"])
        psim.context.setPositions(modeller.positions)
        e_prev = energy_kj(psim)
        while pre_done < int(pre["max_steps"]):
            n = min(interval, int(pre["max_steps"]) - pre_done)
            step(psim.context, n)
            pre_done += n
            e = energy_kj(psim)
            if not math.isfinite(e):
                raise NonFiniteStart(f"Pre-relaxation diverged ({e}) at step {pre_done}")
            if abs(e - e_prev) < delta:
                break
            e_prev = e
        sim.context.setPositions(psim.context.getState(getPositions=True).getPositions())
        del psim
        e_start = energy_kj(sim)
        trace.append((pre_done, e_start))
        log(f"pre-relaxed in {pre_done} steps, {time.monotonic() - t0:.0f}s: "
            f"full-system energy {e_start:.1f} kJ/mol")

    done, probed, prev = 0, False, e_start
    while done < max_steps:
        n = min(interval, max_steps - done)
        step(sim.context, n)
        done += n
        e = energy_kj(sim)
        trace.append((pre_done + done, e))
        bad = not math.isfinite(e) or abs(e) > max_abs
        if not probed and done >= probe_steps:
            probed = True
            bad = bad or (e_start < 0 < e and e > float(mz["divergence_factor"]) * abs(e_start))
        if bad:
            raise RuntimeError(f"Minimization diverged at step {done} (initial={e0:.1f}, "
                               f"current={e:.1f} kJ/mol) -- not saving a garbage structure")
        if done < max_steps and abs(e - prev) < delta:
            log(f"converged at step {done}")
            break
        if probed and budget_s and time.monotonic() - t0 > budget_s:
            log(f"time budget reached at step {done} -- saving unconverged")
            break
        prev = e

    # A heavily clashed model can come out with an inverted CA; the energy
    # does not show it, so say so in the log.
    final = sim.context.getState(getPositions=True).getPositions()
    raw_v = ca_volumes(top, modeller.positions)
    d_res = [f"{k} (was {'D' if raw_v[k] < -1 else 'L' if raw_v[k] > 1 else 'flat'})"
             for k, v in ca_volumes(top, final).items() if v < 0]
    if d_res:
        log(f"WARNING: {len(d_res)} D-amino acid(s) after minimization: {', '.join(d_res)}")

    pdb.parent.mkdir(parents=True, exist_ok=True)
    with open(pdb.with_name(pdb.stem + "_energy.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["step", "energy_kJ_mol"])
        w.writerows(trace)
    tmp = pdb.with_suffix(".tmp")
    with open(tmp, "w") as fh:
        PDBFile.writeFile(top, final, fh, keepIds=True)
    tmp.rename(pdb)
    log(f"saved {pdb.name} (final energy {trace[-1][1]:.1f} kJ/mol, "
        f"{time.monotonic() - t0:.0f}s)")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cif", required=True, type=Path)
    ap.add_argument("--pdb", required=True, type=Path)
    ap.add_argument("--system", required=True)
    ap.add_argument("--ligand-params", type=Path, help="parametrize_ligands.py --out-dir")
    ap.add_argument("--cg2all-bin", type=Path, help="file holding the cg2all env's bin directory")
    a = ap.parse_args()
    mz, spec = load_config()["minimize"], load_system(a.system)
    fb = mz.get("cg2all_fallback") or {}
    try:
        minimize(a.cif, a.pdb, mz, spec, a.ligand_params)
    except NonFiniteStart as e:
        if not (fb.get("enabled") and a.cg2all_bin):
            raise
        log(f"{e} -- rebuilding protein chains with cg2all (CA fixed) and minimizing that")
        bindir = Path(a.cg2all_bin.read_text().strip())
        minimize(a.cif, a.pdb, mz, spec, a.ligand_params,
                 rebuild=lambda m: cg2all_rebuild(m, fb, bindir))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
