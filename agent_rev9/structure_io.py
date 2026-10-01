"""Read standard protein residues; mmCIF columns are resolved by name."""

import gzip
from pathlib import Path


def protein_residues(structure_path: str) -> list[tuple[str, str, str]]:
    from Bio.Data.PDBData import protein_letters_3to1
    from Bio.PDB import PDBParser
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict

    path = Path(structure_path)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as handle:
        if path.name.lower().endswith((".cif", ".cif.gz")):
            data = MMCIF2Dict(handle)
            names = data.get("_atom_site.label_comp_id", [])
            chains = data.get("_atom_site.label_asym_id", data.get("_atom_site.auth_asym_id", []))
            numbers = data.get("_atom_site.label_seq_id", data.get("_atom_site.auth_seq_id", []))
            models = data.get("_atom_site.pdbx_PDB_model_num", ["1"] * len(names))
            rows = [
                (c, n, r)
                for c, n, r, m in zip(chains, numbers, names, models)
                if m == models[0] and n not in (".", "?")
            ]
        elif path.name.lower().endswith((".pdb", ".pdb.gz")):
            structure = PDBParser(QUIET=True).get_structure("protein", handle)
            model = next(structure.get_models(), None)
            rows = (
                []
                if model is None
                else [(c.id, str(r.id[1]) + r.id[2].strip(), r.resname) for c in model for r in c]
            )
        else:
            raise ValueError("Expected .pdb, .cif, .pdb.gz or .cif.gz")
    result = []
    seen = set()
    for chain, number, residue in rows:
        if residue in protein_letters_3to1 and (chain, number) not in seen:
            seen.add((chain, number))
            result.append((chain, number, protein_letters_3to1[residue]))
    return result
