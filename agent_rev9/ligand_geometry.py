from pathlib import Path

import numpy as np

from .screening_gates import cif_atoms
from .support.af3_backbone.structure_qc import AA3_TO_1, _kabsch


def protein_atoms(path, chain):
    """Observed residues in file order; labels retain original residue numbering."""
    path = Path(path)
    groups = {}
    if ".cif" in path.suffixes:
        from Bio.PDB.MMCIF2Dict import MMCIF2Dict

        meta = MMCIF2Dict(str(path))
        count = len(meta.get("_atom_site.label_atom_id", []))
        selected = {}
        columns = {}

        def column(key, default):
            if key not in columns:
                columns[key] = meta.get("_atom_site." + key, [default] * count)
            return columns[key]
        first_model = column("pdbx_PDB_model_num", "1")[0] if count else "1"
        for i in range(count):
            if column("pdbx_PDB_model_num", "1")[i] != first_model:
                continue
            cid = column("label_asym_id", "")[i]
            name, atom = column("label_comp_id", "")[i], column("label_atom_id", "")[i]
            number = column("auth_seq_id", "?")[i]
            if number in {".", "?"}:
                number = column("label_seq_id", "?")[i]
            insertion = column("pdbx_PDB_ins_code", "")[i]
            insertion = "" if insertion in {".", "?"} else insertion
            key = (cid, number, insertion, atom)
            occupancy = column("occupancy", "1")[i]
            occupancy = float(occupancy) if occupancy not in {".", "?"} else 0.0
            xyz = tuple(float(column(axis, "nan")[i]) for axis in ("Cartn_x", "Cartn_y", "Cartn_z"))
            row = (cid, (number, insertion), name, atom, column("type_symbol", "")[i], xyz,
                   f"{cid}:{name}{number}{insertion}")
            if key not in selected or occupancy > selected[key][0]:
                selected[key] = (occupancy, row)
        records = [row for _, row in selected.values()]
    else:
        from Bio.PDB import PDBParser

        model = next(PDBParser(QUIET=True).get_structure("receptor", str(path)).get_models())
        records = [
            (c.id, (r.id[1], r.id[2]), r.resname, a.name, a.element, tuple(a.coord),
             f"{c.id}:{r.resname}{r.id[1]}{r.id[2].strip()}")
            for c in model for r in c if r.id[0] == " " for a in r
        ]
    for c, rid, name, atom, element, xyz, label in records:
        if c != chain or name not in AA3_TO_1 or element.upper() in {"H", "D"}:
            continue
        group = groups.setdefault(
            rid, {"aa": AA3_TO_1[name], "label": label, "atoms": [], "ca": None}
        )
        group["atoms"].append(tuple(float(v) for v in xyz))
        if atom == "CA":
            group["ca"] = tuple(float(v) for v in xyz)
    residues = list(groups.values())
    if not residues:
        raise ValueError("Protein chain has no observed residues")
    return residues


def receptor_alignment(reference_path, moving_path, reference_chain="A", moving_chain="A"):
    from Bio.Align import PairwiseAligner

    reference = protein_atoms(reference_path, reference_chain)
    moving = protein_atoms(moving_path, moving_chain)
    a = "".join(r["aa"] for r in reference)
    b = "".join(r["aa"] for r in moving)
    aligner = PairwiseAligner(
        mode="global",
        match_score=2,
        mismatch_score=-100,
        open_gap_score=-10,
        extend_gap_score=-0.5,
    )
    alignment = aligner.align(a, b)[0]
    mapping = {}
    for ra, rb in zip(*alignment.aligned):
        for i, j in zip(range(*ra), range(*rb)):
            if a[i] != b[j]:
                raise ValueError("Receptor sequence mismatch; residue mapping is ambiguous")
            mapping[j] = i
    if len(mapping) != min(len(a), len(b)):
        raise ValueError("Receptor sequences cannot be fully mapped without substitutions")

    pairs = [
        (reference[i]["ca"], moving[j]["ca"])
        for j, i in mapping.items()
        if reference[i]["ca"] is not None and moving[j]["ca"] is not None
    ]
    if len(pairs) < 3:
        raise ValueError("Fewer than three matched receptor CA atoms")
    ref, mob = np.asarray([x[0] for x in pairs]), np.asarray([x[1] for x in pairs])
    if (
        min(
            np.linalg.matrix_rank(ref - ref.mean(0)),
            np.linalg.matrix_rank(mob - mob.mean(0)),
        )
        < 2
    ):
        raise ValueError("Receptor alignment is geometrically underdetermined")
    rotation, center, reference_center = _kabsch(mob, ref, np)
    aligned = (mob - center) @ rotation.T + reference_center
    return (
        rotation,
        center,
        reference_center,
        mapping,
        reference,
        moving,
        float(np.sqrt(np.mean(np.sum((aligned - ref) ** 2, axis=1)))),
    )


def molecule_from_atoms(elements, coordinates, smiles):
    """Recover connectivity from coordinates, then require the supplied chemical graph.

    No positional atom-index correspondence to the input SMILES is assumed.
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdDetermineBonds

    template = Chem.RemoveHs(Chem.MolFromSmiles(smiles))
    heavy = [
        (element, xyz)
        for element, xyz in zip(elements, coordinates)
        if element.upper() not in {"H", "D"}
    ]
    if len(heavy) != template.GetNumAtoms():
        raise ValueError("Ligand atom count mismatch")
    xyz = str(len(heavy)) + "\nligand\n" + "".join(f"{e} {x[0]} {x[1]} {x[2]}\n" for e, x in heavy)
    mol = Chem.MolFromXYZBlock(xyz)
    if mol is None:
        raise ValueError("Ligand coordinates are malformed")
    rdDetermineBonds.DetermineConnectivity(mol)
    try:
        mol = AllChem.AssignBondOrdersFromTemplate(template, mol)
    except (ValueError, RuntimeError) as exc:
        raise ValueError("Ligand atom mapping failed") from exc

    if mol.GetNumBonds() != template.GetNumBonds():
        raise ValueError("Ligand chemical graph mismatch")

    mapping = mol.GetSubstructMatch(template, useChirality=False)
    if len(mapping) != template.GetNumAtoms():
        raise ValueError("Ligand graph mapping failed")
    result = Chem.Mol(template)
    conformer = Chem.Conformer(result.GetNumAtoms())
    for index, source_index in enumerate(mapping):
        conformer.SetAtomPosition(index, mol.GetConformer().GetAtomPosition(source_index))
    result.RemoveAllConformers()
    result.AddConformer(conformer)
    Chem.AssignStereochemistryFrom3D(result, replaceExistingTags=True)
    return result


def ligand_from_cif(path, chain, smiles):
    atoms = [a for a in cif_atoms(str(path)) if a.label_chain_id == chain and not a.is_hydrogen]
    if not atoms:
        raise ValueError("Ligand chain missing from AF3 model")
    return molecule_from_atoms([a.element for a in atoms], [a.coord for a in atoms], smiles)


def heavy_mol(mol):
    from rdkit import Chem

    result = Chem.RemoveHs(Chem.Mol(mol))
    if not result.GetNumConformers():
        raise ValueError("Ligand has no coordinates")
    xyz = result.GetConformer().GetPositions()
    if not np.all(np.isfinite(xyz)):
        raise ValueError("Non-finite ligand coordinates")
    return result


def chirality_valid(mol, smiles):
    from rdkit import Chem

    reference = Chem.RemoveHs(Chem.MolFromSmiles(smiles))
    observed = heavy_mol(mol)
    if Chem.MolToSmiles(observed, isomericSmiles=False) != Chem.MolToSmiles(
        reference, isomericSmiles=False
    ):
        return False
    try:
        Chem.AssignStereochemistryFrom3D(observed, replaceExistingTags=True)

        if any(str(s.specified) == "Unspecified" for s in Chem.FindPotentialStereo(reference)):
            return None
        return Chem.MolToSmiles(observed, isomericSmiles=True) == Chem.MolToSmiles(
            reference, isomericSmiles=True
        )
    except (ValueError, RuntimeError):
        return None


def aligned_ligand(mol, alignment):
    from rdkit import Chem

    result = heavy_mol(mol)
    rotation, center, reference_center = alignment[:3]
    xyz = (result.GetConformer().GetPositions() - center) @ rotation.T + reference_center
    conformer = Chem.Conformer(result.GetNumAtoms())
    for i, point in enumerate(xyz):
        conformer.SetAtomPosition(i, tuple(point))
    result.RemoveAllConformers()
    result.AddConformer(conformer)
    return result


def symmetry_rmsd(reference, moving):
    """Heavy-atom CalcRMS in the already aligned receptor frame, in angstrom.

    CalcRMS does not refit the ligand, unlike GetBestRMS. Enumerate graph
    automorphisms with chirality; never compare indices by their input order.
    """
    from rdkit.Chem import rdMolAlign

    ref, mob = heavy_mol(reference), heavy_mol(moving)
    if ref.GetNumAtoms() != mob.GetNumAtoms() or ref.GetNumBonds() != mob.GetNumBonds():
        raise ValueError("mapping_failed: atom count mismatch")
    matches = ref.GetSubstructMatches(mob, uniquify=False, useChirality=True, maxMatches=10001)
    if not matches or len(matches) > 10000:
        raise ValueError("mapping_failed: no mapping or automorphism limit exceeded")
    mapping = [list(enumerate(match)) for match in matches]
    return float(rdMolAlign.CalcRMS(mob, ref, map=mapping, symmetrizeConjugatedTerminalGroups=True))


def jaccard(a, b):
    if a is None or b is None:
        return None
    a, b = set(a), set(b)
    return len(a & b) / len(a | b) if a | b else None


def pose_geometry(mol, receptor, box, *, contact_distance, clash_distance, labels=None):
    xyz = heavy_mol(mol).GetConformer().GetPositions()
    protein_xyz = np.asarray([point for residue in receptor for point in residue["atoms"]])
    atom_labels = [
        labels[i] if labels else residue["label"]
        for i, residue in enumerate(receptor)
        for _ in residue["atoms"]
    ]
    distances = np.linalg.norm(protein_xyz[:, None, :] - xyz[None, :, :], axis=2)
    contacts = sorted(
        {
            atom_labels[i]
            for i in np.where(np.any(distances <= contact_distance, axis=1))[0]
            if atom_labels[i]
        }
    )
    center, size = np.asarray(box["center"]), np.asarray(box["size"])
    return {
        "pocket_occupancy": float(np.mean(np.all(np.abs(xyz - center) <= size / 2, axis=1))),
        "minimum_receptor_ligand_distance": float(distances.min()),
        "severe_clash_count": int(np.sum(distances < clash_distance)),
        "contact_residue_count": len(contacts),
        "contact_residues": contacts,
        "binding_site_center_distance": float(np.linalg.norm(xyz.mean(0) - center)),
    }
