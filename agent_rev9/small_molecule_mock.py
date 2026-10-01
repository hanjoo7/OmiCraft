from pathlib import Path

import numpy as np

from .af3_ligand_tool import AF3LigandTool
from .ligand_geometry import heavy_mol, protein_atoms
from .ligand_preparation import read_ligand
from .small_molecule_io import write_json
from .support.af3_backbone.structure_qc import AA3_TO_1


class SyntheticAF3Tool(AF3LigandTool):
    def run(self, *, candidate, context, output_dir, **kwargs):
        directory = Path(output_dir)
        directory.mkdir(parents=True, exist_ok=True)
        input_path = self.prepare_input(candidate, context, directory)

        mol = heavy_mol(read_ligand(candidate["preparation"]["ligand"]["prepared_ligand_path"]))
        xyz = mol.GetConformer().GetPositions()
        xyz = xyz - xyz.mean(0) + np.asarray(context["docking_box"]["center"])
        residues = protein_atoms(context["raw_receptor_path"], context["target_chain"])
        names = {value: key for key, value in AA3_TO_1.items()}
        for seed in self.config.seeds:
            for sample in range(self.config.num_samples):
                path = directory / f"seed-{seed}_sample-{sample}"
                path.mkdir()
                rows = []
                for i, residue in enumerate(residues, 1):
                    for j, point in enumerate(residue["atoms"]):
                        atom = "CA" if tuple(point) == residue["ca"] else f"C{j + 1}"
                        rows.append(
                            f"ATOM {len(rows) + 1} C {atom} {names[residue['aa']]} {self.config.protein_chain_id} {i} {point[0]} {point[1]} {point[2]}"
                        )
                protein_atom_count = len(rows)
                for i, (atom, point) in enumerate(zip(mol.GetAtoms(), xyz), 1):
                    rows.append(
                        f"HETATM {len(rows) + 1} {atom.GetSymbol()} {atom.GetSymbol()}{i} LIG {self.config.ligand_chain_id} . {point[0]} {point[1]} {point[2]}"
                    )
                header = "data_synthetic\nloop_\n" + "".join(
                    "_atom_site." + name + "\n"
                    for name in (
                        "group_PDB",
                        "id",
                        "type_symbol",
                        "label_atom_id",
                        "label_comp_id",
                        "label_asym_id",
                        "label_seq_id",
                        "Cartn_x",
                        "Cartn_y",
                        "Cartn_z",
                    )
                )
                (path / "model.cif").write_text(header + "\n".join(rows) + "\n#\n")
                write_json(
                    path / "summary_confidences.json",
                    {
                        "ranking_score": 0.8,
                        "ptm": 0.8,
                        "iptm": 0.8,
                        "has_clash": False,
                        "chain_ids": [
                            self.config.protein_chain_id,
                            self.config.ligand_chain_id,
                        ],
                        "chain_iptm": [0.8, 0.8],
                        "chain_pair_iptm": [[0.8, 0.8], [0.8, 0.8]],
                        "chain_pair_pae_min": [[0, 2], [2, 0]],
                        "is_mock": True,
                    },
                )
                tokens = [self.config.protein_chain_id] * len(residues) + [
                    self.config.ligand_chain_id
                ] * mol.GetNumAtoms()
                atom_chains = [self.config.protein_chain_id] * protein_atom_count + [
                    self.config.ligand_chain_id
                ] * mol.GetNumAtoms()
                write_json(
                    path / "confidences.json",
                    {
                        "atom_chain_ids": atom_chains,
                        "atom_plddts": [85.0] * len(atom_chains),
                        "token_chain_ids": tokens,
                        "contact_probs": [[0.8] * len(tokens) for _ in tokens],
                        "is_mock": True,
                    },
                )
        result = self.parse_output(directory, candidate["candidate_id"])
        for row in result["samples"]:
            row["is_mock"] = True
        result.update(
            is_mock=True,
            backend="synthetic_af3",
            input_json=str(input_path),
            warning=["SYNTHETIC_AF3_COORDINATES_AND_CONFIDENCES"],
        )
        return result
