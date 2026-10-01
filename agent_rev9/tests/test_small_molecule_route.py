"""Real CPU chemistry/geometry with synthetic model data; no model/network execution."""

import copy
import csv
import json
from pathlib import Path

import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from ..docking_backend import MockDockingBackend
from ..ligand_geometry import (
    jaccard,
    receptor_alignment,
    symmetry_rmsd,
)
from ..ligand_preparation import prepare_ligand, validate_candidate
from ..screening_tools import ScreeningOptions
from ..small_molecule_config import (
    ConsensusThresholds,
    SmallMoleculeConfig,
)
from ..small_molecule_io import unavailable
from ..small_molecule_mock import SyntheticAF3Tool
from ..small_molecule_workflow import run_small_molecule_workflow
from ..structural_consensus import classify_consensus
from . import test_screening_tools as binder_fixtures

isolated = binder_fixtures.isolated
state = binder_fixtures.state
fake_models = binder_fixtures.fake_models


@pytest.fixture
def ligand_case(tmp_path):
    receptor = tmp_path / "raw.pdb"
    receptor.write_text(
        "".join(
            f"ATOM  {i:5d}  CA  ALA A{i:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           C\n"
            for i, (x, y, z) in enumerate([(0, 0, 0), (4, 0, 0), (4, 4, 0), (0, 4, 0)], 1)
        )
        + "END\n"
    )
    context = {
        "raw_receptor_path": str(receptor),
        "target_chain": "A",
        "target_sequence": "AAAA",
        "docking_box": {"center": [2, 2, 4], "size": [14, 14, 14]},
        "gene_name": "SYNTHETIC",
    }

    # Numbers here are synthetic fixture assertions, not biological thresholds.
    config = SmallMoleculeConfig(
        docking={"backend": "mock"},
        af3={"seeds": [11, 22], "num_samples": 1},
        thresholds={
            "minimum_successful_af3_seeds": 2,
            "maximum_cross_seed_pose_rmsd": 0.1,
            "maximum_docking_af3_pose_rmsd": 0.1,
            "minimum_contact_residue_jaccard": 0.5,
            "minimum_pocket_reproducibility": 1,
        },
    )
    options = ScreeningOptions(
        max_gnina_candidates=20, max_af3_candidates=20, max_admet_candidates=20
    )
    return context, config, options


def run_case(tmp_path, ligand_case, candidates=None, tools=None):
    context, config, options = ligand_case
    return run_small_molecule_workflow(
        candidates or [{"candidate_id": "ok", "smiles": "CCO"}],
        context,
        config,
        tmp_path / "run",
        options=options,
        dry_run=False,
        tools=tools,
    )


def test_synthetic_complete_flow(tmp_path, ligand_case):
    result = run_case(tmp_path, ligand_case)
    c = result["candidates"][0]
    assert c["consensus"]["consensus_status"] == "CONSENSUS_SUPPORTED", c["consensus"]
    assert c["validation"]["eligible_for_admet"], c["validation"]
    assert c["admet"]["status"] == "success", c["admet"]
    assert c["critic"]["verdict"] == "NEEDS_REVIEW"
    assert c["is_mock"] and c["verdict"] == "REVIEW"
    assert c["smiles"] == "CCO"
    assert c["af3_ligand"]["observed_sample_count"] == 2
    assert len(c["consensus"]["comparisons"]) == 2
    for name in ("run", "docking", "af3", "consensus", "admet"):
        manifest = json.loads((tmp_path / "run" / f"{name}_manifest.json").read_text())
        assert manifest["route"] == "small_molecule"
    with Path(result["candidate_summary_path"]).open() as handle:
        row = next(csv.DictReader(handle))
    assert row["is_mock"] == "True"
    assert row["docking_score"] == ""


class ScenarioDocking(MockDockingBackend):
    def run(self, job):
        cid = job.request.candidate_id
        if cid == "dockfail":
            return unavailable(
                "backend_failed",
                "synthetic docking failure",
                backend="mock",
                is_mock=True,
                poses=[],
            )
        result = super().run(job)
        if cid == "posefail":
            path = result["poses"][0]["pose_path"]
            mol = Chem.SDMolSupplier(path, removeHs=False)[0]
            conf = mol.GetConformer()
            for i, p in enumerate(conf.GetPositions()):
                conf.SetAtomPosition(i, tuple(p + [40, 0, 0]))
            with Chem.SDWriter(path) as w:
                w.write(mol)
        return result


class ScenarioAF3(SyntheticAF3Tool):
    def run(self, **kwargs):
        cid = kwargs["candidate"]["candidate_id"]
        if cid == "af3fail":
            return unavailable("backend_failed", "synthetic AF3 failure", is_mock=True, samples=[])
        result = super().run(**kwargs)
        if cid in {"discordant", "seedfail"}:
            for row in result["samples"]:
                if cid == "seedfail" and row["seed"] == 11:
                    continue
                path = Path(row["model_path"])
                lines = []
                for line in path.read_text().splitlines():
                    if line.startswith("HETATM"):
                        cells = line.split()
                        cells[-3] = str(float(cells[-3]) + 1.0)
                        line = " ".join(cells)
                    lines.append(line)
                path.write_text("\n".join(lines) + "\n")
        return result


class ScenarioADMET:
    def __init__(self):
        self.calls = []

    def run(self, candidate, *args, **kwargs):
        from ..small_molecule_admet import SmallMoleculeADMETTool

        self.calls.append(candidate["candidate_id"])
        result = SmallMoleculeADMETTool().run(candidate, *args, **kwargs)
        if candidate["candidate_id"] == "partial":
            result["status"] = "partial"
            result["missing_endpoint_count"] = 1
            result["available_endpoint_count"] -= 1
            endpoint = result["endpoints"][-1]["endpoint_name"]
            result["metrics"].pop(endpoint)
            result["endpoints"][-1].update(
                endpoint_value=None, prediction_status="prediction_failed"
            )
        return result


def test_eight_small_molecule_integration_scenarios(tmp_path, ligand_case):
    _, config, _ = ligand_case
    candidates = [
        {"candidate_id": name, "smiles": "not-a-smiles" if name == "invalid" else "CCO"}
        for name in (
            "ok",
            "discordant",
            "invalid",
            "dockfail",
            "posefail",
            "af3fail",
            "seedfail",
            "partial",
        )
    ]
    admet = ScenarioADMET()
    result = run_case(
        tmp_path,
        ligand_case,
        candidates,
        {"docking": ScenarioDocking(), "af3": ScenarioAF3(config.af3), "admet": admet},
    )
    rows = {c["candidate_id"]: c for c in result["candidates"]}
    assert admet.calls == ["ok", "partial"]
    assert rows["discordant"]["consensus"]["consensus_status"] == "DISCORDANT"
    assert rows["discordant"]["validation"]["screening_call"] == "INDETERMINATE"
    assert rows["invalid"]["workflow_status"] == "invalid_input"
    assert rows["dockfail"]["workflow_status"] == "docking_failed"
    assert rows["posefail"]["pose_qc_status"] == "fail"
    assert rows["posefail"]["af3_ligand"]["status"] == "not_run"
    assert rows["af3fail"]["workflow_status"] == "af3_failed"
    assert rows["seedfail"]["consensus"]["af3_consistency_status"] == "inconsistent"
    assert not rows["seedfail"]["validation"]["eligible_for_admet"]
    assert rows["partial"]["workflow_status"] == "partial"
    assert rows["partial"]["admet"]["endpoints"][-1]["endpoint_value"] is None
    assert len(result["candidate_rankings"]) == 8
    assert all(c["is_mock"] and c["verdict"] != "PASS" for c in rows.values())


@pytest.mark.parametrize(
    "smiles,valid,eligible",
    [
        ("CCO", True, True),
        ("broken", False, False),
        ("C=[N](C)(C)C", False, False),
        ("CC.[Na+]", True, False),
        ("[Cu+2]", True, False),
        ("CC(O)F", True, False),
        ("C[C@H](O)F", True, True),
        ("[NH4+]", True, True),
    ],
)
def test_rdkit_validation(smiles, valid, eligible):
    result = validate_candidate({"candidate_id": "x", "smiles": smiles})
    assert result["smiles"] == smiles
    assert result["smiles_valid"] is valid
    assert result["docking_eligible"] is eligible
    if valid:
        assert result["molecular_weight"] is not None
    else:
        assert result["molecular_weight"] is None


def test_preparation_preserves_identity_and_raw(tmp_path):
    c = {"candidate_id": "x", "smiles": "C[C@H](O)F"}
    validation = validate_candidate(c)
    prepared = prepare_ligand(c, validation, tmp_path / "prep")
    raw = Path(prepared["prepared_ligand_path"])
    before = raw.read_bytes()
    c["raw_ligand_path"] = str(raw)
    copied = prepare_ligand(c, validation, tmp_path / "copy")
    assert copied["status"] == "success"
    assert copied["prepared_ligand_path"] != str(raw)
    assert raw.read_bytes() == before
    changed = validate_candidate({**c, "prepared_smiles": "C[C@@H](O)F"})
    assert changed["smiles_valid"] and not changed["docking_eligible"]
    assert "CHEMICAL_STATE_CHANGE_REQUIRES_SEPARATE_CANDIDATE" in changed["warning"]


def test_bad_id_duplicate_unknown_route_stop_before_backends(tmp_path, ligand_case):
    assert not validate_candidate({"smiles": "CC"})["smiles_valid"]
    for candidates in (
        [{"candidate_id": "x", "smiles": "CC", "route": "wrong"}],
        [{"candidate_id": "x", "smiles": "CC"}] * 2,
        [{"candidate_id": "x", "smiles": "CC", "route": "protein_binder_rfd3"}],
    ):
        with pytest.raises(ValueError):
            run_case(tmp_path, ligand_case, candidates)


def embedded(smiles):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMolecule(mol, randomSeed=7)
    return Chem.RemoveHs(mol)


def test_symmetry_rmsd_does_not_refit_or_compare_indices():
    mol = embedded("c1ccccc1")
    reordered = Chem.RenumberAtoms(mol, [3, 2, 1, 0, 5, 4])
    assert symmetry_rmsd(mol, reordered) == pytest.approx(0, abs=1e-8)
    shifted = Chem.Mol(reordered)
    for i, xyz in enumerate(shifted.GetConformer().GetPositions()):
        shifted.GetConformer().SetAtomPosition(i, tuple(xyz + [3, 0, 0]))
    assert symmetry_rmsd(mol, shifted) == pytest.approx(3)
    with pytest.raises(ValueError, match="mapping_failed"):
        symmetry_rmsd(mol, embedded("CCO"))
    assert jaccard(["A:ALA1", "A:GLY2"], ["A:ALA1"]) == 0.5
    assert jaccard([], []) is None
    assert jaccard(None, ["x"]) is None


def test_receptor_alignment_rigid_transform(tmp_path, ligand_case):
    source = Path(ligand_case[0]["raw_receptor_path"])
    moved = tmp_path / "moved.pdb"
    lines = []
    for line in source.read_text().splitlines():
        if line.startswith("ATOM"):
            x, y, z = (float(line[i : i + 8]) for i in (30, 38, 46))
            line = line[:30] + f"{-y + 7:8.3f}{x + 3:8.3f}{z + 2:8.3f}" + line[54:]
        lines.append(line)
    moved.write_text("\n".join(lines) + "\n")
    assert receptor_alignment(source, moved)[-1] < 1e-8


@pytest.mark.parametrize(
    "changes,status",
    [
        ({}, "CONSENSUS_SUPPORTED"),
        ({"docking_af3_pose_rmsd": 2}, "DISCORDANT"),
        ({"af3_cross_seed_pose_rmsd": 2}, "POCKET_SUPPORTED_POSE_UNCERTAIN"),
        ({"af3_supported": False}, "DOCKING_SUPPORTED_ONLY"),
        ({"docking_supported": False}, "AF3_SUPPORTED_ONLY"),
        (
            {
                "af3_supported": False,
                "docking_supported": False,
                "both_structurally_failed": True,
            },
            "STRUCTURALLY_UNSUPPORTED",
        ),
        ({"mapping_status": "mapping_failed"}, "INDETERMINATE"),
        ({"docking_af3_pose_rmsd": None}, "INDETERMINATE"),
    ],
)
def test_consensus_statuses(changes, status):
    t = ConsensusThresholds(
        minimum_successful_af3_seeds=2,
        maximum_cross_seed_pose_rmsd=0.1,
        maximum_docking_af3_pose_rmsd=0.1,
    )
    metrics = {
        "docking_supported": True,
        "af3_supported": True,
        "mapping_status": "success",
        "af3_successful_seed_count": 2,
        "af3_cross_seed_pose_rmsd": 0,
        "docking_af3_pose_rmsd": 0,
        **changes,
    }
    assert classify_consensus(metrics, t)["consensus_status"] == status
    assert classify_consensus(metrics, ConsensusThresholds())["threshold_checks"] == {}


def test_missing_metrics_disallow_pass_and_audit_mixture(tmp_path, ligand_case):
    from ..small_molecule_critic import audit_candidate
    from ..validation_agent import validate_small_molecule

    c = run_case(tmp_path, ligand_case)["candidates"][0]
    c["af3_ligand"]["samples"][0]["ligand_chain_iptm"] = None
    v = validate_small_molecule(c, ligand_case[1].validation)
    assert v["screening_call"] == "INDETERMINATE" and not v["eligible_for_admet"]
    c["docking"]["is_mock"] = False
    assert "MOCK_REAL_MIXTURE" in audit_candidate(c)["contradictions"]
    c["config_unchanged"] = False
    assert "CONFIG_CHANGED_DURING_RUN" in audit_candidate(c)["contradictions"]


def test_manifest_handoff_integrity_and_relative_paths(tmp_path, ligand_case):
    from ..small_molecule_io import redact, resolve_path, sha256
    from ..support.admet_pipeline.handoff import load_handoff

    result = run_case(tmp_path, ligand_case)
    path = tmp_path / "run" / "run_manifest.json"
    manifest = json.loads(path.read_text())
    assert sha256(resolve_path(manifest["summary_path"], path)) == manifest["summary_sha256"]
    for entry in manifest["inputs"]:
        assert sha256(resolve_path(entry["path"], path)) == entry["sha256"]
    hmanifest = Path(result["candidates"][0]["admet"]["handoff_manifest"])
    hcsv = hmanifest.parent / "pipeline1_handoff.csv"
    if not hcsv.exists():
        hcsv = next(hmanifest.parent.glob("*.csv"))
    assert load_handoff(hcsv, hmanifest).rows[0]["candidate_id"] == "ok"
    hcsv.write_text(hcsv.read_text() + "\n")
    with pytest.raises(ValueError, match="(?i)checksum|sha256"):
        load_handoff(hcsv, hmanifest)
    assert redact(
        {
            "api_key": "secret",
            "model_dir": "/weights",
            "environment": {"TOKEN": "secret"},
        }
    ) == {
        "api_key": "<redacted>",
        "model_dir": "<redacted>",
        "environment": "<redacted>",
    }


def test_native_missing_backend_never_falls_back(tmp_path, ligand_case):
    context, config, options = ligand_case
    config.docking.backend = "gnina"
    context["prepared_receptor_path"] = context["raw_receptor_path"]
    c = run_case(tmp_path, ligand_case)["candidates"][0]
    assert c["docking"]["status"] == "dependency_missing"
    assert not c["is_mock"]
    assert c["af3_ligand"]["status"] == "not_run"


def test_parser_preserves_all_samples_and_missing_confidence(tmp_path, ligand_case):
    from ..af3_ligand_tool import AF3LigandTool

    c = run_case(tmp_path, ligand_case)["candidates"][0]
    a = c["af3_ligand"]
    first = a["samples"][0]
    assert first["ligand_atom_plddt_mean"] == 85
    assert first["ligand_chain_iptm"] == 0.8
    assert first["protein_ligand_pae_min"] == 2
    assert first["contact_probability_summary"]["mean"] == pytest.approx(0.8)
    Path(first["confidences_path"]).unlink()
    result = AF3LigandTool(ligand_case[1].af3).parse_output(
        Path(first["model_path"]).parent.parent, "ok"
    )
    assert result["status"] == "partial" and len(result["samples"]) == 2
    assert result["samples"][0]["ligand_atom_plddt_mean"] is None
    assert result["is_mock"] is True


def test_input_is_independent_and_ccd_identity_checked(tmp_path, ligand_case):
    from ..af3_ligand_tool import AF3LigandTool

    tool = AF3LigandTool(ligand_case[1].af3)
    candidate = {
        "candidate_id": "x",
        "smiles": "CCO",
        "canonical_smiles": "CCO",
        "pose_path": "DO_NOT_USE.sdf",
    }
    path = tool.prepare_input(candidate, ligand_case[0], tmp_path)
    raw = path.read_text()
    payload = json.loads(raw)
    assert "DO_NOT_USE" not in raw and payload["modelSeeds"] == [11, 22]
    assert payload["sequences"][1]["ligand"] == {"id": "B", "smiles": "CCO"}
    ccd = tmp_path / "CCD.cif"
    ccd.write_text(
        "data_X\n_chem_comp.id X\nloop_\n_pdbx_chem_comp_descriptor.type\n_pdbx_chem_comp_descriptor.descriptor\nSMILES 'CCN'\n#\n"
    )
    with pytest.raises(ValueError, match="IDENTITY_MISMATCH"):
        tool.prepare_input(
            {**candidate, "ccd_code": "X", "ccd_path": str(ccd)},
            ligand_case[0],
            tmp_path,
        )


def test_ninth_integration_binder_route_stays_on_existing_cascade(state, fake_models, monkeypatch):
    from .. import screening_agent
    from ..docking_backend import NativeDockingBackend
    from ..small_molecule_admet import SmallMoleculeADMETTool

    def forbidden(*args, **kwargs):
        raise AssertionError("Small-molecule tool received a binder")

    monkeypatch.setattr(SmallMoleculeADMETTool, "run", forbidden)
    monkeypatch.setattr(NativeDockingBackend, "run", forbidden)
    state.update(route="protein_binder_rfd3", use_rfd3=True, dry_run=True)
    result = screening_agent.screening_node(state)["screening_result"]
    assert fake_models[0], result
    assert result["candidates"]
    assert all(c["route"] == "protein_binder_rfd3" for c in result["candidates"])
    assert all("rfd3_structure_path" in c and "af3_validation" in c for c in result["candidates"])
    assert all(c["admet_status"] == "skipped_not_applicable" for c in result["candidates"])


def test_conformer_failure_and_geometry_failures(tmp_path, ligand_case, monkeypatch):
    monkeypatch.setattr(AllChem, "EmbedMolecule", lambda *args, **kwargs: -1)
    r = validate_candidate({"candidate_id": "x", "smiles": "CCO"})
    assert r["conformer_status"] == "failed" and not r["af3_eligible"]


def test_native_docking_parser_keeps_scores_separate_and_all_poses(tmp_path):
    from ..docking_backend import NativeDockingBackend
    from ..small_molecule_config import DockingConfig

    first = embedded("CCO")
    second = Chem.Mol(first)
    first.SetProp("minimizedAffinity", "-7.1")
    first.SetProp("CNNscore", ".7")
    first.SetProp("CNNaffinity", "6.2")
    second.SetProp("minimizedAffinity", "NaN")
    with Chem.SDWriter(str(tmp_path / "poses.sdf")) as w:
        w.write(first)
        w.write(second)
    with (tmp_path / "poses.sdf").open("a") as handle:
        handle.write("broken\n$$$$\n")
    rows = NativeDockingBackend(DockingConfig(backend="gnina")).parse(tmp_path, "CCO")
    assert len(rows) == 3 and rows[2]["docking_status"] == "backend_failed"
    assert (
        rows[0]["docking_score"],
        rows[0]["cnn_score"],
        rows[0]["cnn_affinity"],
    ) == (-7.1, 0.7, 6.2)
    assert rows[1]["docking_score"] is None and rows[1]["cnn_score"] is None
    pdb = Chem.MolToPDBBlock(first)
    block = "\n".join(
        line + " " + first.GetAtomWithIdx(int(line[6:11]) - 1).GetSymbol()
        for line in pdb.splitlines()
        if line.startswith("HETATM")
    )
    (tmp_path / "poses.pdbqt").write_text(
        "MODEL 1\nREMARK VINA RESULT: -6.5 0 0\n" + block + "\nENDMDL\nMODEL 2\nBROKEN\nENDMDL\n"
    )
    rows = NativeDockingBackend(DockingConfig(backend="vina")).parse(tmp_path, "CCO")
    assert len(rows) == 2 and rows[0]["docking_score"] == -6.5, rows
    assert rows[0]["cnn_score"] is None and rows[0]["cnn_affinity"] is None
    assert rows[1]["docking_status"] == "backend_failed"


def test_admet_real_parser_preserves_partial_and_absence(tmp_path):
    from ..small_molecule_admet import SmallMoleculeADMETTool
    from ..small_molecule_io import write_csv, write_json

    write_json(tmp_path / "run_manifest.json", {"backend": "admet_ai", "is_mock": False})
    write_csv(
        tmp_path / "candidate_summary.csv",
        [
            {
                "candidate_id": "x",
                "prediction_status": "partial",
                "available_endpoint_count": 1,
                "missing_endpoint_count": 1,
            }
        ],
    )
    write_csv(
        tmp_path / "admet_predictions_long.csv",
        [
            {
                "candidate_id": "x",
                "endpoint_name": "A",
                "endpoint_value": 0.5,
                "prediction_status": "success",
                "is_mock": False,
            },
            {
                "candidate_id": "x",
                "endpoint_name": "B",
                "endpoint_value": None,
                "prediction_status": "failed",
                "is_mock": False,
            },
        ],
    )
    parsed = SmallMoleculeADMETTool.parse_output(tmp_path, "x")
    assert parsed["status"] == "partial" and parsed["metrics"] == {"A": 0.5}
    assert parsed["missing_endpoint_count"] == 1
    write_json(tmp_path / "run_manifest.json", {"backend": "admet_ai", "is_mock": True})
    with pytest.raises(ValueError, match="provenance"):
        SmallMoleculeADMETTool.parse_output(tmp_path, "x")


def test_cli_preserves_venv_interpreter_and_dry_run(tmp_path, ligand_case):
    from ..small_molecule_cli import local_paths, main
    from ..small_molecule_io import write_json

    runtime = tmp_path / "venv/bin"
    runtime.mkdir(parents=True)
    interpreter = runtime / "python"
    interpreter.symlink_to("/usr/bin/python3")
    assert local_paths({"python_path": "venv/bin/python"}, tmp_path)["python_path"] == str(
        interpreter
    )
    context, config, _ = ligand_case
    config.docking.backend = "gnina"
    payload = {
        "route": "small_molecule",
        "context": context,
        "candidates": [{"candidate_id": "x", "smiles": "CCO"}],
    }
    inp = write_json(tmp_path / "input.json", payload)
    cfg = write_json(tmp_path / "config.json", {"small_molecule": config.model_dump()})
    assert (
        main(
            [
                "--input",
                str(inp),
                "--config",
                str(cfg),
                "--output",
                str(tmp_path / "planned"),
                "--gpu-device",
                "4",
            ]
        )
        == 0
    )
    assert list((tmp_path / "planned").rglob("af3_ligand_input.json"))
    saved = json.loads((tmp_path / "planned/candidates.json").read_text())[0]
    assert saved["af3_ligand"]["status"] == "not_run"
    assert not saved["af3_ligand"].get("command")


def test_all_seeds_including_clashes_prevent_selective_pass(tmp_path, ligand_case):
    from ..validation_agent import validate_small_molecule

    c = run_case(tmp_path, ligand_case)["candidates"][0]
    c["af3_ligand"]["samples"][1]["has_clash"] = True
    v = validate_small_molecule(c, ligand_case[1].validation)
    assert not v["eligible_for_admet"] and v["screening_call"] == "INDETERMINATE"


def test_optional_qc_adapter_contract_and_fingerprint(tmp_path, ligand_case, monkeypatch):
    import sys
    import types

    import pandas as pd

    from .. import pose_qc

    class FakePB:
        def __init__(self, **kwargs):
            pass

        def bust(self, **kwargs):
            return pd.DataFrame([{"bond_lengths": True, "stereochemistry": True}])

    class FakeFingerprint:
        def generate(self, *args, **kwargs):
            return {("LIG1.B", "ALA1.A"): {"Hydrophobic": True}}

    monkeypatch.setitem(sys.modules, "posebusters", types.SimpleNamespace(PoseBusters=FakePB))
    monkeypatch.setitem(
        sys.modules,
        "prolif",
        types.SimpleNamespace(
            Fingerprint=FakeFingerprint,
            Molecule=types.SimpleNamespace(from_rdkit=lambda m: m),
        ),
    )
    original = pose_qc.importlib.util.find_spec
    monkeypatch.setattr(
        pose_qc.importlib.util,
        "find_spec",
        lambda name: object() if name in {"prolif", "posebusters"} else original(name),
    )
    monkeypatch.setattr(
        Chem,
        "MolFromPDBFile",
        lambda *args, **kwargs: Chem.AddHs(Chem.MolFromSmiles("CC")),
    )
    c = run_case(tmp_path, ligand_case)["candidates"][0]
    assert c["pose_qc"][0]["posebusters_status"] == "pass"
    assert c["pose_qc"][0]["prolif_status"] == "pass"
    assert c["consensus"]["interaction_fingerprint_similarity"] == 1


def test_prepared_receptor_cannot_silently_change_coordinate_frame(tmp_path, ligand_case):
    from ..ligand_preparation import prepare_receptor

    context = dict(ligand_case[0])
    source = Path(context["raw_receptor_path"])
    changed = tmp_path / "changed.pdb"
    changed.write_text(
        source.read_text().replace("   0.000   0.000   0.000", "  40.000   0.000   0.000")
    )
    context["prepared_receptor_path"] = str(changed)
    result = prepare_receptor(context, tmp_path / "prep")
    assert result["status"] == "invalid_input"
    assert result["prepared_receptor_path"] is None


def test_ranking_is_stable_and_failed_admet_is_not_normal(tmp_path, ligand_case):
    from ..candidate_ranking import rank_candidates
    from ..small_molecule_critic import audit_candidate

    c = run_case(tmp_path, ligand_case)["candidates"][0]
    other = copy.deepcopy(c)
    other["candidate_id"] = "failed"
    other["admet"].update(
        status="prediction_failed",
        available_endpoint_count=None,
        missing_endpoint_count=None,
    )
    rows = rank_candidates([other, c])
    assert rows[0]["candidate_id"] == "ok" and rows[1]["ranking_status"] == "admet_failed"
    assert rank_candidates([c, other]) == rows
    other["ranking"]["ranking_status"] = "ranked"
    assert "FAILED_ADMET_RANKED_NORMAL" in audit_candidate(other)["contradictions"]


def test_cpu_profile_does_not_launch_af3(tmp_path, ligand_case):
    from unittest.mock import Mock

    from ..af3_ligand_tool import AF3LigandTool
    from ..small_molecule_admet import SmallMoleculeADMETTool

    context, config, options = ligand_case
    from dataclasses import replace

    ligand_case = context, config, replace(options, resource_profile="C")
    af3 = Mock(spec=AF3LigandTool)
    admet = Mock(spec=SmallMoleculeADMETTool)
    result = run_case(tmp_path, ligand_case, tools={"af3": af3, "admet": admet})
    af3.run.assert_not_called()
    admet.run.assert_not_called()
    assert result["counts"]["af3_tool_calls"] == 0


def test_mock_input_cannot_launch_a_real_backend(tmp_path, ligand_case):
    from unittest.mock import Mock

    context, config, _ = ligand_case
    config.docking.backend = "gnina"
    context["prepared_receptor_path"] = context["raw_receptor_path"]
    backend = Mock()
    backend.name = "gnina"
    result = run_case(
        tmp_path,
        ligand_case,
        [{"candidate_id": "synthetic", "smiles": "CCO", "is_mock": True}],
        {"docking": backend},
    )
    backend.prepare.assert_not_called()
    backend.run.assert_not_called()
    assert (
        result["candidates"][0]["docking"]["error_message"]
        == "MOCK_INPUT_REQUIRES_EXPLICIT_MOCK_BACKEND"
    )


def test_independent_admet_preserves_uncertain_structure(tmp_path, ligand_case):
    from ..small_molecule_admet import SmallMoleculeADMETTool
    _, config, _ = ligand_case
    config.admet_policy = 'valid_input'
    result = run_case(tmp_path, ligand_case,
                      [{'candidate_id': 'af3fail', 'smiles': 'CCO'}],
                      {'af3': ScenarioAF3(config.af3), 'admet': SmallMoleculeADMETTool()})
    candidate = result['candidates'][0]
    assert candidate['af3_ligand']['status'] == 'backend_failed'
    assert not candidate['validation']['eligible_for_admet']
    assert candidate['admet']['status'] == 'success'
    assert candidate['admet']['execution_policy'] == 'valid_input'
    assert candidate['final_verdict'] != 'PASS'
    handoff = Path(candidate['admet']['handoff_manifest'])
    metadata = json.loads(handoff.read_text())
    assert metadata['structural_eligible_for_admet'] is False
    assert metadata['admet_execution_policy'] == 'valid_input'
    summary = handoff.parent.parent / 'candidate_summary.csv'
    assert 'structurally uncertain' in summary.read_text()
    assert 'AF3-supported pose' not in summary.read_text()


def test_independent_admet_does_not_bypass_invalid_input_or_veto(tmp_path, ligand_case):
    from ..small_molecule_admet import eligible_for_property_prediction
    _, config, _ = ligand_case
    config.admet_policy = 'valid_input'
    result = run_case(tmp_path, ligand_case, [{'candidate_id': 'bad', 'smiles': 'not-smiles'}])
    assert result['counts']['admet_tool_calls'] == 0
    valid = {'route': 'small_molecule', 'rdkit_validation': {'smiles_valid': True},
             'validation': {'eligible_for_admet': False}, 'admet_execution_policy': 'valid_input',
             'qualification': {'safety_verdict': 'REJECT'}}
    assert not eligible_for_property_prediction(valid)
    valid['qualification'] = {}
    valid['design_workflow'] = {'safety_veto': True}
    assert not eligible_for_property_prediction(valid)


def test_progress_records_backend_failure_and_final_report(tmp_path, ligand_case):
    context, config, options = ligand_case
    events = []
    result = run_small_molecule_workflow(
        [{'candidate_id': 'af3fail', 'smiles': 'CCO'}], context, config,
        tmp_path / 'stages', options=options, dry_run=False,
        tools={'af3': ScenarioAF3(config.af3)}, progress=lambda *row: events.append(row))
    assert ('af3', 'FAILED', 'af3fail') in events
    assert events[-1][:2] == ('report', 'COMPLETED')
    assert result['candidates'][0]['final_verdict'] != 'PASS'
