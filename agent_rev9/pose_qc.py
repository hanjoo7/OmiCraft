import importlib.util

from .ligand_geometry import chirality_valid, heavy_mol, pose_geometry, protein_atoms
from .ligand_preparation import read_ligand


class PoseQCTool:
    def __init__(self, config):
        self.config = config

    def run(self, pose, context, smiles):
        row = {
            "candidate_id": pose.get("candidate_id"),
            "pose_rank": pose.get("pose_rank"),
            "pose_qc_status": "not_run",
            "pocket_occupancy": None,
            "minimum_receptor_ligand_distance": None,
            "severe_clash_count": None,
            "chirality_valid": None,
            "geometry_valid": None,
            "contact_residue_count": None,
            "contact_residues": None,
            "interaction_fingerprint": None,
            "posebusters_status": "not_run",
            "prolif_status": "not_run",
            "warning": [],
            "error_message": "",
        }
        if pose.get("docking_status") != "success":
            return row
        try:
            from rdkit import Chem

            mol = heavy_mol(read_ligand(pose["pose_path"]))
            expected = Chem.RemoveHs(Chem.MolFromSmiles(smiles))
            if mol.GetNumAtoms() != expected.GetNumAtoms():
                raise ValueError("LIGAND_ATOM_MISSING_OR_EXTRA")
            if Chem.MolToSmiles(mol, isomericSmiles=False) != Chem.MolToSmiles(
                expected, isomericSmiles=False
            ):
                raise ValueError("LIGAND_CHEMICAL_GRAPH_MISMATCH")

            Chem.SanitizeMol(mol)
            row["geometry_valid"] = True
            row["chirality_valid"] = chirality_valid(mol, smiles)
            receptor = protein_atoms(context["raw_receptor_path"], context["target_chain"])
            row.update(
                pose_geometry(
                    mol,
                    receptor,
                    context["docking_box"],
                    contact_distance=self.config.contact_distance_angstrom,
                    clash_distance=self.config.severe_clash_distance_angstrom,
                )
            )
            failures = []
            if row["pocket_occupancy"] < 1:
                failures.append("POSE_OUTSIDE_POCKET")
            if row["severe_clash_count"]:
                failures.append("SEVERE_ATOMIC_CLASH")
            if row["chirality_valid"] is False:
                failures.append("CHIRALITY_MISMATCH")
            if row["chirality_valid"] is None:
                row["warning"].append("CHIRALITY_NOT_ASSESSED")
            expected_contacts = set(context.get("expected_pocket_residues") or [])
            row["expected_pocket_contact_count"] = (
                len(expected_contacts & set(row["contact_residues"])) if expected_contacts else None
            )
            if expected_contacts and not row["expected_pocket_contact_count"]:
                row["warning"].append("EXPECTED_POCKET_CONTACT_MISSING")
            self._optional(row, pose["pose_path"], context["prepared_receptor_path"], mol)
            if row["posebusters_status"] == "fail":
                failures.append("POSEBUSTERS_FAILED")
            row["pose_qc_status"] = (
                "fail" if failures else "pass_with_warning" if row["warning"] else "pass"
            )
            row["error_message"] = "; ".join(failures)
        except (OSError, ValueError, RuntimeError, KeyError) as exc:
            row.update(pose_qc_status="fail", geometry_valid=False, error_message=str(exc))
        return row

    def _optional(self, row, ligand_path, receptor_path, mol):
        if self.config.run_posebusters:
            if importlib.util.find_spec("posebusters") is None:
                row["posebusters_status"] = "dependency_missing"
                row["warning"].append("POSEBUSTERS_MISSING")
            else:
                try:
                    from posebusters import PoseBusters

                    table = PoseBusters(config="dock", top_n=None).bust(
                        mol_pred=ligand_path, mol_cond=receptor_path
                    )
                    tests = [
                        bool(value)
                        for value in table.iloc[0].values
                        if type(value).__name__ in {"bool", "bool_"}
                    ]
                    row["posebusters_checks"] = {
                        str(k): bool(v)
                        for k, v in table.iloc[0].items()
                        if type(v).__name__ in {"bool", "bool_"}
                    }
                    row["posebusters_status"] = (
                        "pass" if tests and all(tests) else "fail" if tests else "not_run"
                    )
                except Exception as exc:
                    row["posebusters_status"] = "not_run"
                    row["warning"].append("POSEBUSTERS_ERROR: " + str(exc))
        if self.config.run_prolif:
            fingerprint = prolif_fingerprint(mol, receptor_path)
            row["warning"].extend(fingerprint.pop("warning"))
            row.update(fingerprint)
        if row["posebusters_status"] != "pass":
            row["warning"].append("FULL_BOND_STRAIN_QC_NOT_CONFIRMED")


def prolif_fingerprint(mol, receptor_path):
    """Compare poses against one explicitly prepared reference receptor.

    AddHs derives ligand hydrogens from supplied valence only, in a temporary
    molecule; it does not guess a pH or alter raw/model structures.
    """
    result = {
        "interaction_fingerprint": None,
        "prolif_status": "not_run",
        "warning": [],
    }
    if importlib.util.find_spec("prolif") is None:
        result.update(prolif_status="dependency_missing", warning=["PROLIF_MISSING"])
        return result
    try:
        import prolif as plf
        from rdkit import Chem

        protein = Chem.MolFromPDBFile(receptor_path, removeHs=False)
        if protein is None or not any(a.GetAtomicNum() == 1 for a in protein.GetAtoms()):
            raise ValueError(
                "A prepared protein PDB with explicit hydrogens is required for hydrogen-dependent interactions"
            )
        ligand = Chem.AddHs(Chem.Mol(mol), addCoords=True)
        data = plf.Fingerprint().generate(
            plf.Molecule.from_rdkit(ligand),
            plf.Molecule.from_rdkit(protein),
            metadata=True,
        )
        result.update(
            interaction_fingerprint=sorted(
                f"{residues[1]}:{name}"
                for residues, interactions in data.items()
                for name in interactions
            ),
            prolif_status="pass",
            warning=["IFP_COMMON_PREPARED_RECEPTOR;LIGAND_H_FROM_SUPPLIED_VALENCE"],
            fingerprint_method="ProLIF against common prepared reference receptor; ligand aligned to that frame",
        )
    except Exception as exc:
        result["warning"].append("PROLIF_ERROR: " + str(exc))
    return result
