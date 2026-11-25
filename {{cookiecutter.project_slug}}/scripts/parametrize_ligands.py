#!/usr/bin/env python3
"""
scripts/parametrize_ligands.py — Stage 3d'
==========================================
GAFF2 + AM1-BCC parameters (antechamber / parmchk2 / sqm, the same AmberTools
route ChimeraX's addcharge takes) for every free `type: ligand` chain of
configs/<system>.yaml. Run once per system, not once per model.

Not parametrized here (see lib_pipeline.ligand_kind):
  * monatomic ions      -> OpenMM's own amber14 ion templates
  * covalent ligands    (named in `bonds:`) -> dropped before minimization
  * multi-CCD ligands   -> dropped before minimization

Chemistry comes from the config: `smiles`, or the CCD definition of `ccd`
(downloaded from RCSB, cached in --ccd-cache). Protonation is rule-based for
pH 7 (acids off, aliphatic amines / amidines on); `protonation: as_is` on the
ligand entry keeps the config/CCD state instead.

Outputs in --out-dir: ligands.json (chain -> files), gaff.xml, <chain>.xml
(OpenMM residue template), <chain>.sdf (template with explicit charges).
"""
from __future__ import annotations

import argparse
import json
import shutil
import sys
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

import yaml
from rdkit import Chem
from rdkit.Chem import AllChem

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib_pipeline import ligand_code, ligand_kind, norm_id  # noqa: E402

GAFF = "gaff-2.2.20"
CCD_URL = "https://files.rcsb.org/ligands/download/{}.cif"
_BO = {"SING": Chem.BondType.SINGLE, "DOUB": Chem.BondType.DOUBLE, "TRIP": Chem.BondType.TRIPLE}

_ACIDS = ["[CX3](=O)[OX2H1:1]", "[P,S](=O)[OX2H1:1]"]
_AMIDINE = "[CX3;!$(C=O)](-[NX3])=[NX2;!$(N~[N,O]);!$(N-a):1]"
_AMINE = ("[NX3;+0;!$(N~[!#6;!#1]);!$(N-[#6]=,:[#7,#8,#16]);!$(N-a);"
          "!$(N-C=C);!$(N-C#N);!$(N#*);!$(N=*):1]")


def mol_from_ccd(code: str, cache: Path) -> Chem.Mol:
    import gemmi

    path = cache / f"{code}.cif"
    if not path.exists():
        cache.mkdir(parents=True, exist_ok=True)
        tmp = path.with_suffix(".tmp")
        with urllib.request.urlopen(CCD_URL.format(code), timeout=60) as r:
            tmp.write_bytes(r.read())
        tmp.rename(path)
    block = gemmi.cif.read(str(path)).sole_block()
    rw, idx = Chem.RWMol(), {}
    for name, el, q in block.find("_chem_comp_atom.", ["atom_id", "type_symbol", "charge"]):
        a = Chem.Atom(el.capitalize())
        a.SetFormalCharge(int(q) if q not in "?." else 0)
        a.SetNoImplicit(True)
        a.SetProp("atomName", gemmi.cif.as_string(name))
        idx[gemmi.cif.as_string(name)] = rw.AddAtom(a)
    for n1, n2, order in block.find("_chem_comp_bond.", ["atom_id_1", "atom_id_2", "value_order"]):
        rw.AddBond(idx[gemmi.cif.as_string(n1)], idx[gemmi.cif.as_string(n2)], _BO[order.upper()])
    mol = rw.GetMol()
    Chem.SanitizeMol(mol)
    mol = Chem.RemoveHs(mol)
    for a in mol.GetAtoms():
        a.SetNoImplicit(False)
    Chem.SanitizeMol(mol)
    return mol


def protonate(mol: Chem.Mol) -> Chem.Mol:
    mol = Chem.RWMol(mol)

    def hits(smarts):
        q = Chem.MolFromSmarts(smarts)
        k = next(a.GetIdx() for a in q.GetAtoms() if a.GetAtomMapNum() == 1)
        return sorted({m[k] for m in mol.GetSubstructMatches(q)})

    def shift(i, dq):
        a = mol.GetAtomWithIdx(i)
        n = a.GetTotalNumHs()
        a.SetNoImplicit(True)
        a.SetNumExplicitHs(n + dq)
        a.SetFormalCharge(a.GetFormalCharge() + dq)
        mol.UpdatePropertyCache()

    for s in _ACIDS:
        for i in hits(s):
            shift(i, -1)
    # one charge per amidine / per cluster of amines (piperazine is +1, not +2)
    dist, taken = Chem.GetDistanceMatrix(mol), []
    for s in (_AMIDINE, _AMINE):
        for i in hits(s):
            if all(dist[i][j] > 3 for j in taken):
                shift(i, 1)
                taken.append(i)
    Chem.SanitizeMol(mol)
    return mol.GetMol()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--spec", required=True)
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--ccd-cache", required=True, type=Path)
    a = ap.parse_args()

    from openff.toolkit import AmberToolsToolkitWrapper, Molecule
    from openmmforcefields.generators import GAFFTemplateGenerator
    import openmmforcefields

    spec = yaml.safe_load(Path(a.spec).read_text())
    a.out_dir.mkdir(parents=True, exist_ok=True)
    base = list(Path(openmmforcefields.__file__).parent.rglob(f"{GAFF}.xml"))
    if len(base) != 1:
        raise RuntimeError(f"{GAFF}.xml not found in openmmforcefields: {base}")
    shutil.copy(base[0], a.out_dir / "gaff.xml")

    index, done = {}, {}
    for s in spec["sequences"]:
        if s.get("type") != "ligand" or ligand_kind(s, spec) != "gaff":
            continue
        chain = norm_id(s["id"])
        as_is = s.get("protonation") == "as_is"
        key = (ligand_code(s), s.get("smiles"), as_is)
        if key not in done:
            mol = (Chem.MolFromSmiles(s["smiles"]) if s.get("smiles")
                   else mol_from_ccd(ligand_code(s), a.ccd_cache))
            if mol is None:
                raise ValueError(f"chain {chain}: unreadable SMILES {s['smiles']!r}")
            if not as_is:
                mol = protonate(mol)
            molh = Chem.AddHs(mol)
            if AllChem.EmbedMolecule(molh, randomSeed=1) != 0 and \
               AllChem.EmbedMolecule(molh, randomSeed=1, useRandomCoords=True) != 0:
                raise RuntimeError(f"chain {chain}: no 3D conformer")
            AllChem.UFFOptimizeMolecule(molh)
            print(f"[ligand] chain {chain}: {Chem.MolToSmiles(mol)} "
                  f"(net charge {Chem.GetFormalCharge(mol):+d}) -> {GAFF} + AM1-BCC", flush=True)
            off = Molecule.from_rdkit(molh, allow_undefined_stereo=True)
            off.assign_partial_charges("am1bcc", use_conformers=off.conformers,
                                       toolkit_registry=AmberToolsToolkitWrapper())
            # template names are molecule-derived, so identical ligands share one
            gen = GAFFTemplateGenerator(molecules=[off], forcefield=GAFF)
            # gaff.xml already defines the atom types and their LJ terms, and
            # OpenMM refuses to load them twice
            root = ET.fromstring(gen.generate_residue_template(off))
            for tag in ("AtomTypes", "NonbondedForce"):
                for el in root.findall(tag):
                    root.remove(el)
            (a.out_dir / f"{chain}.xml").write_text(ET.tostring(root, encoding="unicode"))
            # atom names travel as a property list: SDF has no name column
            mol.SetProp("atom_names", ",".join(
                x.GetProp("atomName") if x.HasProp("atomName") else "" for x in mol.GetAtoms()))
            with Chem.SDWriter(str(a.out_dir / f"{chain}.sdf")) as w:
                w.write(Chem.RemoveHs(mol))
            done[key] = chain
        index[chain] = {"params": f"{done[key]}.xml", "mol": f"{done[key]}.sdf"}
    (a.out_dir / "ligands.json").write_text(json.dumps(index, indent=1) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
