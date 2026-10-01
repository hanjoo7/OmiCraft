"""Ligand validation and preparation in the supplied chemical state."""

import shutil
from pathlib import Path

from .small_molecule_io import sha256, write_json
from .support.pipeline_common.chem_validation import (
    validate_smiles_rdkit,
)


def validate_candidate(candidate, *, seed=11):
    original = candidate.get("smiles", "")
    result = {
        "candidate_id": candidate.get("candidate_id"),
        "smiles": original,
        "canonical_smiles": None,
        "prepared_smiles": candidate.get("prepared_smiles"),
        "smiles_valid": False,
        "rdkit_parse_status": "invalid_input",
        "molecular_weight": None,
        "formal_charge": None,
        "fragment_count": None,
        "stereochemistry_status": "not_assessed",
        "metal_present": None,
        "conformer_status": "not_run",
        "validation_message": "",
        "docking_eligible": False,
        "af3_eligible": False,
        "warning": [],
        "transformations": [],
    }
    if not isinstance(candidate.get("candidate_id"), str) or not candidate["candidate_id"].strip():
        result["validation_message"] = "candidate_id is required"
        return result
    if not isinstance(original, str):
        result["validation_message"] = "SMILES must be a string"
        return result
    basic = validate_smiles_rdkit(candidate["candidate_id"], original)
    result.update(
        smiles_valid=basic.smiles_valid,
        rdkit_parse_status=basic.rdkit_parse_status,
        canonical_smiles=basic.canonical_smiles or None,
        validation_message=basic.validation_message,
    )
    if not basic.smiles_valid:
        return result
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors

    mol = Chem.MolFromSmiles(basic.canonical_smiles)
    result.update(
        molecular_weight=Descriptors.MolWt(mol),
        formal_charge=Chem.GetFormalCharge(mol),
        fragment_count=len(Chem.GetMolFrags(mol)),
    )
    nonmetals = {
        1,
        2,
        5,
        6,
        7,
        8,
        9,
        10,
        14,
        15,
        16,
        17,
        18,
        32,
        33,
        34,
        35,
        36,
        52,
        53,
        54,
        85,
        86,
    }
    result["metal_present"] = any(a.GetAtomicNum() not in nonmetals for a in mol.GetAtoms())
    stereo = Chem.FindPotentialStereo(mol)
    undefined = any(str(s.specified) == "Unspecified" for s in stereo)
    result["stereochemistry_status"] = (
        "unspecified" if undefined else "specified" if stereo else "none_required"
    )
    if result["fragment_count"] != 1:
        result["warning"].append("MULTI_FRAGMENT_REQUIRES_EXPLICIT_PREPARATION")
    if result["metal_present"]:
        result["warning"].append("METAL_CONTAINING_LIGAND_REQUIRES_REVIEW")
    if undefined:
        result["warning"].append("UNSPECIFIED_STEREOCHEMISTRY_NOT_ENUMERATED")
    if result["prepared_smiles"]:
        prepared = Chem.MolFromSmiles(result["prepared_smiles"])
        if (
            prepared is None
            or Chem.MolToSmiles(prepared, isomericSmiles=True) != basic.canonical_smiles
        ):
            result["warning"].append("CHEMICAL_STATE_CHANGE_REQUIRES_SEPARATE_CANDIDATE")
            result["validation_message"] = (
                "Changed prepared SMILES must be an explicit variant with its own provenance"
            )
            return result
        result["transformations"].append("user_provided_identical_chemical_state")

    if result["fragment_count"] != 1 or result["metal_present"] or undefined:
        return result
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    try:
        conformer = Chem.AddHs(mol)
        result["conformer_status"] = (
            "success" if AllChem.EmbedMolecule(conformer, params) == 0 else "failed"
        )
    except (ValueError, RuntimeError):
        result["conformer_status"] = "failed"
    result["docking_eligible"] = result["af3_eligible"] = result["conformer_status"] == "success"
    if not result["docking_eligible"]:
        result["validation_message"] = "3D conformer generation failed"
    result["warning"] += ["PROTONATION_AS_SUPPLIED", "TAUTOMER_NOT_ENUMERATED"]
    return result


def chemical_identity(mol):
    from rdkit import Chem

    return Chem.MolToSmiles(Chem.RemoveHs(mol), canonical=True, isomericSmiles=True)


def read_ligand(path):
    from rdkit import Chem

    molecules = list(Chem.SDMolSupplier(str(path), removeHs=False))
    if len(molecules) != 1 or molecules[0] is None or not molecules[0].GetNumConformers():
        raise ValueError("Expected one parseable 3D SDF ligand")
    return molecules[0]


def prepare_ligand(candidate, validation, output_dir, *, seed=11):
    try:
        from rdkit import Chem, rdBase
        from rdkit.Chem import AllChem
    except ImportError:
        result = {
            "status": "dependency_missing",
            "prepared_ligand_path": None,
            "error_message": "RDKit is required",
            "warning": [],
        }
        write_json(Path(output_dir) / "ligand_preparation.json", result)
        return result
    directory = Path(output_dir)
    directory.mkdir(parents=True, exist_ok=True)
    result = {
        "status": "invalid_input",
        "raw_ligand_path": candidate.get("raw_ligand_path") or candidate.get("ligand_path"),
        "prepared_ligand_path": None,
        "original_smiles": candidate.get("smiles"),
        "canonical_smiles": validation["canonical_smiles"],
        "prepared_smiles": validation["prepared_smiles"],
        "conformer_generation_method": None,
        "protonation_method": "as_supplied",
        "target_ph": None,
        "charge_assignment": "formal_charge_from_input",
        "stereochemistry_preservation": None,
        "tool": "RDKit",
        "version": rdBase.rdkitVersion,
        "warning": list(validation["warning"]),
        "transformations": [],
        "error_message": "",
    }
    try:
        if not validation["docking_eligible"]:
            raise ValueError(
                validation["validation_message"] or "Explicit supported chemical state required"
            )
        if result["raw_ligand_path"]:
            mol = read_ligand(result["raw_ligand_path"])
            if chemical_identity(mol) != validation["canonical_smiles"]:
                raise ValueError("Supplied ligand differs from validated SMILES")
            result["conformer_generation_method"] = "provided_3d_coordinates"
            result["raw_sha256"] = sha256(result["raw_ligand_path"])
        else:
            mol = Chem.AddHs(Chem.MolFromSmiles(validation["canonical_smiles"]))
            params = AllChem.ETKDGv3()
            params.randomSeed = seed
            if AllChem.EmbedMolecule(mol, params) != 0:
                raise ValueError("RDKit conformer generation failed")
            result["conformer_generation_method"] = "ETKDGv3"
            result["transformations"] = [
                "explicit_hydrogens_from_supplied_valence",
                "ETKDGv3_coordinates",
            ]
        if chemical_identity(mol) != validation["canonical_smiles"]:
            raise ValueError("Ligand chemical identity changed during preparation")
        path = directory / "prepared_ligand.sdf"
        if (
            result["raw_ligand_path"]
            and path.resolve() == Path(result["raw_ligand_path"]).resolve()
        ):
            raise ValueError("Raw ligand cannot be overwritten")
        with Chem.SDWriter(str(path)) as writer:
            writer.write(mol)
        result.update(
            status="success",
            prepared_smiles=chemical_identity(mol),
            prepared_ligand_path=str(path.resolve()),
            prepared_sha256=sha256(path),
            stereochemistry_preservation=True,
        )
    except (OSError, ValueError, RuntimeError) as exc:
        result["error_message"] = str(exc)
    write_json(directory / "ligand_preparation.json", result)
    return result


def prepare_receptor(context, output_dir, *, is_mock=False):
    """Accept an explicitly prepared receptor, retaining an untouched raw source.

    Chemical protonation/repair is not inferred from merely having a PDB file.
    """
    directory = Path(output_dir)
    directory.mkdir(parents=True, exist_ok=True)
    raw = context.get("raw_receptor_path") or context.get("receptor_path")
    provided = context.get("prepared_receptor_path") or (
        raw if context.get("receptor_prepared") is True else None
    )
    metadata = context.get("receptor_preparation") or {}
    if not isinstance(metadata, dict):
        raise ValueError("receptor_preparation must be an object")
    result = {
        "status": "not_run",
        "raw_receptor_path": raw,
        "prepared_receptor_path": None,
        "chain": context.get("target_chain"),
        "selected_biological_assembly": metadata.get("biological_assembly"),
        "water_policy": metadata.get("water_policy", "not_assessed"),
        "cofactor_policy": metadata.get("cofactor_policy", "not_assessed"),
        "missing_atoms": metadata.get("missing_atoms", "not_assessed"),
        "missing_residues": metadata.get("missing_residues", "not_modeled"),
        "added_hydrogens": metadata.get("added_hydrogens"),
        "ph": metadata.get("ph"),
        "protonation_method": metadata.get("protonation_method", "not_assessed"),
        "histidine_tautomers": metadata.get("histidine_tautomers", "not_assessed"),
        "tool": metadata.get("tool", "provided_copy"),
        "version": metadata.get("version"),
        "is_mock": is_mock,
        "warning": [],
        "error_message": "",
    }
    if not raw or not Path(raw).is_file():
        result.update(status="invalid_input", error_message="Raw receptor is missing")
    elif not provided and not is_mock:
        result.update(
            status="dependency_missing",
            error_message="Provide an explicitly prepared receptor and its preparation metadata",
        )
    else:
        source = Path(provided or raw)
        if not source.is_file() or source.suffix.lower() not in {".pdb", ".pdbqt"}:
            result.update(
                status="invalid_input",
                error_message="Prepared docking receptor must be PDB or PDBQT",
            )
        else:
            destination = directory / ("prepared_receptor" + source.suffix)
            if (
                destination.resolve() == Path(raw).resolve()
                or destination.resolve() == source.resolve()
            ):
                raise ValueError("Receptor preparation cannot overwrite an input")
            try:
                import numpy as np

                from .ligand_geometry import protein_atoms

                reference = {
                    r["label"]: r
                    for r in protein_atoms(raw, context["target_chain"])
                    if r["ca"] is not None
                }
                prepared = {
                    r["label"]: r
                    for r in protein_atoms(source, context["target_chain"])
                    if r["ca"] is not None
                }
                shared = reference.keys() & prepared.keys()

                if len(shared) < 3 or any(
                    not np.allclose(
                        reference[label]["ca"],
                        prepared[label]["ca"],
                        atol=0.002,
                        rtol=0,
                    )
                    for label in shared
                ):
                    raise ValueError(
                        "Prepared receptor must retain the raw binding-site coordinate frame and residue identifiers"
                    )
            except (OSError, ValueError, KeyError) as exc:
                result.update(status="invalid_input", error_message=str(exc))
                write_json(directory / "receptor_preparation.json", result)
                return result
            shutil.copyfile(source, destination)
            result.update(
                status="success",
                prepared_receptor_path=str(destination.resolve()),
                raw_sha256=sha256(raw),
                prepared_sha256=sha256(destination),
            )
            if result["protonation_method"] == "not_assessed":
                result["warning"].append("RECEPTOR_PROTONATION_NOT_ASSESSED")
            if is_mock:
                result["warning"].append("MOCK_RECEPTOR_COPY_NO_PROTONATION")
    write_json(directory / "receptor_preparation.json", result)
    return result
