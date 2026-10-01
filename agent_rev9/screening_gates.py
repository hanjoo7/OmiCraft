from __future__ import annotations

import gzip
import math
import tempfile
from contextlib import contextmanager
from dataclasses import asdict, dataclass, field
from pathlib import Path

from .protein_design_route import extract_protein_sequence
from .structure_io import protein_residues


@dataclass(frozen=True)
class CascadeThresholds:
    secondary_iptm_min: float = 0.6
    secondary_ptm_min: float = 0.5
    binding_site_jaccard_min: float = 0.5
    binder_ca_rmsd_max: float = 5.0
    gnina_energy_max: float = -5.0
    gnina_cnn_min: float = 0.5
    gnina_cnn_fail_below: float = 0.2
    gnina_cnnaffinity_min: float = 5.0
    boltz_ligand_iptm_min: float = 0.5
    boltz_probability_min: float = 0.5

    admet_limits: dict = field(default_factory=dict)

    def __post_init__(self):
        for name, value in asdict(self).items():
            if name == "admet_limits":
                if not isinstance(value, dict):
                    raise ValueError("admet_limits must be an endpoint -> {min,max} object")
                for endpoint, limits in value.items():
                    if (
                        not isinstance(endpoint, str)
                        or not isinstance(limits, dict)
                        or not limits
                        or set(limits) - {"min", "max"}
                    ):
                        raise ValueError("Invalid ADMET endpoint limits")
                    if any(not finite(v) for v in limits.values()) or limits.get(
                        "min", -math.inf
                    ) > limits.get("max", math.inf):
                        raise ValueError("Invalid ADMET numeric limits")
            elif not finite(value):
                raise ValueError(f"Invalid threshold: {name}")
            elif (
                name.endswith(
                    (
                        "_iptm_min",
                        "_ptm_min",
                        "_jaccard_min",
                        "_cnn_min",
                        "_fail_below",
                        "_probability_min",
                    )
                )
                and not 0 <= value <= 1
            ):
                raise ValueError(f"{name} must be in [0, 1]")
        if self.binder_ca_rmsd_max <= 0 or self.gnina_cnn_fail_below > self.gnina_cnn_min:
            raise ValueError("Invalid RMSD or GNINA threshold ordering")


def finite(value):
    return type(value) in (int, float) and math.isfinite(value)


def verdict_result(verdict, reasons=(), **metrics):
    return {"verdict": verdict, "failure_reasons": list(dict.fromkeys(reasons)), **metrics}


def overlap(first, second):
    a, b = set(first), set(second)
    return len(a & b) / len(a | b) if a or b else None


@contextmanager
def decoded_structure(path):
    path = Path(path)
    if path.suffix not in (".gz", ".zst"):
        yield path
        return
    with tempfile.TemporaryDirectory(prefix="screening-cif-") as directory:
        decoded = Path(directory) / path.stem
        if path.suffix == ".gz":
            with gzip.open(path, "rb") as stream:
                decoded.write_bytes(stream.read())
        else:
            import zstandard

            with (
                path.open("rb") as stream,
                zstandard.ZstdDecompressor().stream_reader(stream) as reader,
            ):
                decoded.write_bytes(reader.read())
        yield decoded


def cif_atoms(path):
    from .support.af3_backbone.structure_qc import parse_cif_atoms

    with decoded_structure(path) as decoded:
        return parse_cif_atoms(decoded)


def sequence_for_chain(path, chain):
    with decoded_structure(path) as decoded:
        return extract_protein_sequence(decoded, chain)


def contacts(target_xyz, partner_xyz, residue_ids, cutoff):
    import numpy as np

    a, b = np.asarray(target_xyz, dtype=float), np.asarray(partner_xyz, dtype=float)
    if not len(a) or not len(b) or not np.isfinite(a).all() or not np.isfinite(b).all():
        raise ValueError("Missing or non-finite structure coordinates")
    nearest = np.full(len(a), np.inf)
    count = 0
    for start in range(0, len(a), 256):
        distances = np.linalg.norm(a[start : start + 256, None, :] - b[None, :, :], axis=2)
        nearest[start : start + 256] = distances.min(axis=1)
        count += int((distances <= cutoff).sum())
    by_residue = {}
    for residue, distance in zip(residue_ids, nearest):
        if residue is not None:
            by_residue[residue] = min(by_residue.get(residue, math.inf), float(distance))
    return {
        "residue_min_distance_angstrom": by_residue,
        "min_interchain_distance_angstrom": float(nearest.min()),
        "interface_contact_count": count,
        "interface_target_residue_ids": sorted(
            {r for r, d in zip(residue_ids, nearest) if r is not None and d <= cutoff}
        ),
    }


def receptor_atoms(path, chain, target_sequence):
    """Exact observed-residue mapping only; preserve source residue identities."""
    from Bio.PDB import PDBParser

    residues = [(n, aa) for c, n, aa in protein_residues(path) if c == chain]
    if "".join(aa for _, aa in residues) != target_sequence:
        raise ValueError("TARGET_SEQUENCE_MAPPING_UNAVAILABLE")
    mapping = {number: i + 1 for i, (number, _) in enumerate(residues)}
    model = next(PDBParser(QUIET=True).get_structure("receptor", path).get_models())
    atoms, indices = [], []
    for residue in model[chain]:
        number = str(residue.id[1]) + residue.id[2].strip()
        if number not in mapping:
            continue
        for atom in residue:
            if atom.element not in {"H", "D"} and atom.get_altloc() in {" ", "A"}:
                atoms.append(atom.coord)
                indices.append(mapping[number])
    return atoms, indices


def ligand_pose(path, expected_smiles=None):
    import numpy as np
    from rdkit import Chem

    mols = [m for m in Chem.SDMolSupplier(str(path), removeHs=False) if m is not None]
    if len(mols) != 1:
        raise ValueError("Selected ligand pose must contain exactly one valid molecule")
    mol = mols[0]
    if not mol.GetNumConformers() or not mol.GetConformer().Is3D():
        raise ValueError("Ligand pose has no 3D conformer")
    canonical = Chem.MolToSmiles(Chem.RemoveHs(mol), isomericSmiles=True)
    if expected_smiles:
        expected = Chem.MolFromSmiles(expected_smiles)
        if expected is None or canonical != Chem.MolToSmiles(expected, isomericSmiles=True):
            raise ValueError("Ligand pose chemical identity differs from input SMILES")
    coords = np.array(
        [
            mol.GetConformer().GetAtomPosition(a.GetIdx())
            for a in mol.GetAtoms()
            if a.GetAtomicNum() > 1
        ]
    )
    if not len(coords) or not np.isfinite(coords).all():
        raise ValueError("Invalid ligand coordinates")
    return mol, coords


class BinderCrossModelGate:
    def __init__(self, thresholds, structural_thresholds):
        self.thresholds, self.structural_thresholds = thresholds, structural_thresholds

    def run(self, primary, secondary, context):
        if primary["verdict"] == "FAIL":
            return verdict_result("FAIL", primary["failure_reasons"], status="SKIPPED_UPSTREAM")
        if secondary.get("status") != "success":
            return verdict_result(
                "REVIEW",
                [secondary.get("error_message") or "PROTENIX_NOT_EXECUTED"],
                status="NOT_EVALUATED",
            )
        from .screening_tools import StructuralScreeningTool

        t, m = self.thresholds, secondary.get("metrics", {})
        reasons = []
        for key, minimum in (("iptm", t.secondary_iptm_min), ("ptm", t.secondary_ptm_min)):
            if not finite(m.get(key)) or not minimum <= m[key] <= 1:
                reasons.append("PROTENIX_" + key.upper() + "_LOW_OR_MISSING")
        if m.get("has_clash") is not False:
            reasons.append("PROTENIX_CLASH_OR_CLASH_STATUS_MISSING")
        geometry = StructuralScreeningTool(self.structural_thresholds)._geometry(
            m.get("model_path"), context
        )
        agreement = orientation = None
        hotspot = None
        try:
            for chain, expected in (
                ("A", context["target_sequence"]),
                ("B", context["binder_sequence"]),
            ):
                if sequence_for_chain(m["model_path"], chain) != expected:
                    raise ValueError("PROTENIX_SEQUENCE_MAPPING_UNAVAILABLE")
                if sequence_for_chain(primary["af3_metrics"]["model_path"], chain) != expected:
                    raise ValueError("AF3_SEQUENCE_MAPPING_UNAVAILABLE")
            if geometry["status"] != "COMPLETED":
                raise ValueError("PROTENIX_GEOMETRY_UNAVAILABLE")
            if (
                geometry["interface_contact_count"]
                < self.structural_thresholds.min_interface_contacts
            ):
                reasons.append("PROTENIX_INSUFFICIENT_INTERFACE_CONTACTS")
            if (
                geometry["min_interchain_distance_angstrom"]
                < self.structural_thresholds.severe_clash_distance_angstrom
            ):
                reasons.append("PROTENIX_SEVERE_CLASH")
            agreement = overlap(
                primary["structural_metrics"].get("interface_target_residue_ids", []),
                geometry.get("interface_target_residue_ids", []),
            )
            if agreement is None or agreement < t.binding_site_jaccard_min:
                reasons.append("BINDING_SITE_DISAGREEMENT_OR_MISSING")
            orientation = self._binder_rmsd(primary["af3_metrics"]["model_path"], m["model_path"])
            if orientation is None or orientation > t.binder_ca_rmsd_max:
                reasons.append("BINDER_ORIENTATION_DISAGREEMENT_OR_UNAVAILABLE")
            if context.get("hotspot_residues"):
                a, b = primary.get("hotspot_preserved"), geometry.get("hotspot_preserved")
                hotspot = a == b if a is not None and b is not None else None
                if a is not True or b is not True:
                    reasons.append("HOTSPOT_AGREEMENT_UNCONFIRMED")
        except (OSError, KeyError, ValueError, TypeError) as exc:
            reasons.append(str(exc))
        verdict = "PASS" if primary["verdict"] == "PASS" and not reasons else "REVIEW"
        return verdict_result(
            verdict,
            primary["failure_reasons"] + reasons,
            status="COMPLETED",
            binding_site_agreement=agreement,
            interface_agreement=agreement,
            hotspot_agreement=hotspot,
            binder_ca_rmsd_angstrom=orientation,
            secondary_structural_metrics=geometry,
        )

    @staticmethod
    def _binder_rmsd(reference_path, moving_path):
        import numpy as np

        from .support.af3_backbone.structure_qc import _kabsch, _protein_ca_by_residue

        a, b = cif_atoms(reference_path), cif_atoms(moving_path)
        reference = _protein_ca_by_residue(a, "A")
        moving = _protein_ca_by_residue(b, "A")
        keys = sorted(set(reference) & set(moving))
        if len(keys) < 3:
            return None
        ref, mob = np.array([reference[k] for k in keys]), np.array([moving[k] for k in keys])
        if (
            np.linalg.matrix_rank(ref - ref.mean(axis=0)) < 2
            or np.linalg.matrix_rank(mob - mob.mean(axis=0)) < 2
        ):
            return None
        rotation, center, ref_center = _kabsch(mob, ref, np)
        reference = _protein_ca_by_residue(a, "B")
        moving = _protein_ca_by_residue(b, "B")
        if not reference or set(reference) != set(moving):
            return None
        keys = sorted(reference)

        aligned = (np.array([moving[k] for k in keys]) - center) @ rotation.T + ref_center
        return float(
            np.sqrt(
                np.mean(np.sum((aligned - np.array([reference[k] for k in keys])) ** 2, axis=1))
            )
        )


class PoseScreeningTool:
    def __init__(self, thresholds, structural_thresholds):
        self.thresholds, self.structural_thresholds = thresholds, structural_thresholds

    def run(self, docking, context):
        if docking.get("status") != "success":
            return verdict_result(
                "FAIL" if docking.get("status") == "failed" else "REVIEW",
                [docking.get("error_message") or "GNINA_NOT_EXECUTED"],
                status="NOT_EVALUATED",
            )
        import numpy as np

        fail, review, geometry = [], [], {}
        try:
            _, xyz = ligand_pose(docking["pose_path"], context["smiles"])
            box = context["docking_box"]
            center, size = np.array(box["center"]), np.array(box["size"])
            inside = bool(np.all(np.abs(xyz - center) <= size / 2))
            if not inside:
                fail.append("POSE_OUTSIDE_BINDING_SITE")
            atoms, residues = receptor_atoms(
                context["receptor_path"], context["target_chain"], context["target_sequence"]
            )
            geometry = contacts(
                atoms, xyz, residues, self.structural_thresholds.contact_distance_angstrom
            )
            geometry["inside_binding_site"] = inside
            if (
                geometry["min_interchain_distance_angstrom"]
                < self.structural_thresholds.severe_clash_distance_angstrom
            ):
                fail.append("POSE_SEVERE_CLASH")
            if (
                geometry["interface_contact_count"]
                < self.structural_thresholds.min_interface_contacts
            ):
                fail.append("POSE_NO_TARGET_CONTACT")
        except (OSError, KeyError, ValueError, TypeError) as exc:
            review.append("POSE_GEOMETRY_UNAVAILABLE: " + str(exc))
        scores, t = docking.get("scores", {}), self.thresholds
        cnn = scores.get("CNNscore")
        if finite(cnn) and 0 <= cnn < t.gnina_cnn_fail_below:
            fail.append("GNINA_CNN_BELOW_HARD_LIMIT")
        elif not finite(cnn) or not t.gnina_cnn_min <= cnn <= 1:
            review.append("GNINA_CNN_LOW_OR_MISSING")
        energy = scores.get("minimizedAffinity")
        if not finite(energy) or energy > t.gnina_energy_max:
            review.append("GNINA_ENERGY_LOW_OR_MISSING")
        affinity = scores.get("CNNaffinity")
        if not finite(affinity) or affinity < t.gnina_cnnaffinity_min:
            review.append("GNINA_CNNAFFINITY_LOW_OR_MISSING")
        if context.get("receptor_prepared") is not True:
            review.append("RECEPTOR_PREPARATION_UNCONFIRMED")
        return verdict_result(
            "FAIL" if fail else "REVIEW" if review else "PASS",
            fail + review,
            status="COMPLETED",
            interface_metrics=geometry,
        )

    def boltz_consistency(self, pose_gate, boltz, context):
        if boltz.get("status") != "success":
            return verdict_result(
                "REVIEW",
                [boltz.get("error_message") or "BOLTZ2_NOT_EXECUTED"],
                status="NOT_EVALUATED",
            )
        reasons, geometry, agreement = [], {}, None
        m, t = boltz.get("metrics", {}), self.thresholds
        for key, minimum in (
            ("ligand_iptm", t.boltz_ligand_iptm_min),
            ("affinity_probability_binary", t.boltz_probability_min),
        ):
            if not finite(m.get(key)) or not minimum <= m[key] <= 1:
                reasons.append("BOLTZ2_" + key.upper() + "_LOW_OR_MISSING")
        if not finite(m.get("affinity_pred_value")):
            reasons.append("BOLTZ2_AFFINITY_MISSING")
        try:
            path = m["model_path"]
            if sequence_for_chain(path, "A") != context["target_sequence"]:
                raise ValueError("BOLTZ2_TARGET_MAPPING_UNAVAILABLE")
            atoms = cif_atoms(path)
            target = [a for a in atoms if a.label_chain_id == "A" and not a.is_hydrogen]
            ligand = [a for a in atoms if a.label_chain_id == "B" and not a.is_hydrogen]
            from collections import Counter

            from rdkit import Chem

            molecule = Chem.MolFromSmiles(context["smiles"])
            expected = Counter(
                a.GetSymbol().upper() for a in molecule.GetAtoms() if a.GetAtomicNum() > 1
            )
            observed = Counter(a.element.upper() for a in ligand)
            if expected != observed:
                raise ValueError("BOLTZ2_LIGAND_ELEMENT_COMPOSITION_MISMATCH")
            geometry = contacts(
                [a.coord for a in target],
                [a.coord for a in ligand],
                [a.label_seq_id for a in target],
                self.structural_thresholds.contact_distance_angstrom,
            )
            if (
                geometry["min_interchain_distance_angstrom"]
                < self.structural_thresholds.severe_clash_distance_angstrom
            ):
                return verdict_result(
                    "FAIL", ["BOLTZ2_SEVERE_CLASH"], status="COMPLETED", interface_metrics=geometry
                )
            if (
                geometry["interface_contact_count"]
                < self.structural_thresholds.min_interface_contacts
            ):
                return verdict_result(
                    "FAIL",
                    ["BOLTZ2_NO_TARGET_CONTACT"],
                    status="COMPLETED",
                    interface_metrics=geometry,
                )
            agreement = overlap(
                pose_gate.get("interface_metrics", {}).get("interface_target_residue_ids", []),
                geometry["interface_target_residue_ids"],
            )
            if agreement is None or agreement < t.binding_site_jaccard_min:
                reasons.append("GNINA_BOLTZ2_POCKET_DISAGREEMENT_OR_MISSING")
        except (OSError, KeyError, ValueError, TypeError) as exc:
            reasons.append(str(exc))
        return verdict_result(
            "REVIEW" if reasons else "PASS",
            reasons,
            status="COMPLETED",
            binding_site_agreement=agreement,
            interface_agreement=agreement,
            pose_agreement=None,
            pose_agreement_status="NOT_EVALUATED_LIGAND_ATOM_MAPPING",
            interface_metrics=geometry,
        )

    def admet(self, output):
        if output.get("status") != "success":
            return verdict_result("REVIEW", [output.get("error_message") or "ADMET_NOT_EXECUTED"])
        limits = self.thresholds.admet_limits
        if not limits:
            return verdict_result("REVIEW", ["ADMET_THRESHOLDS_NOT_CONFIGURED"])
        failed, missing = [], []
        values = output.get("metrics", {})
        for name, bound in limits.items():
            value = values.get(name)
            if not finite(value):
                missing.append("ADMET_ENDPOINT_MISSING: " + name)
            elif value < bound.get("min", -math.inf) or value > bound.get("max", math.inf):
                failed.append("ADMET_LIMIT_EXCEEDED: " + name)
        return verdict_result(
            "FAIL" if failed else "REVIEW" if missing else "PASS", failed + missing
        )
