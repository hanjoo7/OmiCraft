"""Both cascades use stubbed models and real deterministic gates. No inference."""

import json
from pathlib import Path
from unittest.mock import Mock

import pytest
from rdkit import Chem

from .. import screening_agent as screening
from ..critic_agent import critic_node
from ..protein_design_route import AlphaFold3BinderAdapter
from ..screening_gates import BinderCrossModelGate
from ..screening_model_tools import ADMETTool, Boltz2Tool, GNINATool, ProtenixTool
from ..screening_tools import ScreeningOptions
from .test_screening_tools import complex_cif


@pytest.fixture(autouse=True)
def legacy_boltz_contract(request):
    """Keep the original Boltz cascade contract as an explicit compatibility test."""
    from ..configuration import get_config

    request.getfixturevalue("isolated")
    get_config().small_molecule.structural_backend = "boltz2_legacy"


def write_ligand(path, offset=0.0, score=0.9):
    path.parent.mkdir(parents=True, exist_ok=True)
    mol = Chem.MolFromSmiles("CCO")
    conformer = Chem.Conformer(3)
    conformer.Set3D(True)
    for index, xyz in enumerate(((5.0, 3.0, 0.0), (6.4, 3.0, 0.2), (7.3, 3.8, 0.4))):
        conformer.SetAtomPosition(index, (xyz[0] + offset, xyz[1], xyz[2]))
    mol.AddConformer(conformer)
    for name, value in {"CNNscore": score, "CNNaffinity": 7.0, "minimizedAffinity": -8.0}.items():
        mol.SetProp(name, str(value))
    with Chem.SDWriter(str(path)) as writer:
        writer.write(mol)
    return path


def ligand_complex(path, offset=0.0):
    path.parent.mkdir(parents=True, exist_ok=True)
    header = complex_cif().split("ATOM 1")[0]
    rows = [
        line + "\n"
        for line in complex_cif().splitlines()
        if line.startswith("ATOM") and line.split()[5] == "A"
    ]
    for i, (element, name, xyz) in enumerate(
        (("C", "C1", (5.0, 3.0, 0.0)), ("C", "C2", (6.4, 3.0, 0.2)), ("O", "O1", (7.3, 3.8, 0.4))),
        start=9,
    ):
        rows.append(f"HETATM {i} {element} {name} LIG B . {xyz[0] + offset} {xyz[1]} {xyz[2]}\n")
    path.write_text(header + "".join(rows) + "#\n")
    return path


@pytest.fixture
def sm_state(tmp_path):
    receptor = tmp_path / "receptor.pdb"
    receptor.write_text(
        "".join(
            f"ATOM  {n:5d}  CA  ALA A{n:4d}    {n * 3.8:8.3f}{(n % 3) * 0.4:8.3f}{(n % 2) * 0.3:8.3f}  1.00 80.00           C\n"
            for n in range(1, 9)
        )
        + "END\n"
    )
    design = {
        "gene_name": "SM",
        "modality": "SMALL_MOLECULE",
        "target_sequence": "A" * 8,
        "target_chain": "A",
        "receptor_path": str(receptor),
        "receptor_prepared": True,
        "docking_box": {"center": [15, 3, 0], "size": [40, 20, 20]},
        "candidates": [{"candidate_id": "ligand-1", "smiles": "CCO"}],
    }
    path = tmp_path / "sm_design.json"
    path.write_text(json.dumps({"designs": [design]}))
    return {
        "advance_targets": [{"gene_name": "SM", "modality": "SMALL_MOLECULE"}],
        "design_results_path": str(path),
        "dry_run": False,
        "use_rfd3": False,
        "screening_options": {
            "cascade_thresholds": {"admet_limits": {"hERG": {"max": 0.5}, "AMES": {"max": 0.5}}}
        },
    }


@pytest.fixture
def fake_sm_models(monkeypatch):
    calls = []

    def docking(self, *, candidate, context, output_dir, dry_run):
        calls.append(("GNINA", candidate["candidate_id"], dry_run))
        pose = write_ligand(Path(output_dir) / "poses.sdf")
        return {
            "status": "success",
            **GNINATool.parse_output(pose),
            "is_mock": dry_run,
            "provenance": {"inference_started": not dry_run},
        }

    def boltz(self, *, target_sequence, smiles, output_dir, dry_run, msa_path=None):
        calls.append(("Boltz-2", smiles, dry_run))
        model = ligand_complex(Path(output_dir) / "model.cif")
        return {
            "status": "success",
            "metrics": {
                "model_path": str(model),
                "confidence_score": 0.9,
                "ligand_iptm": 0.8,
                "affinity_probability_binary": 0.9,
                "affinity_pred_value": -0.5,
            },
            "is_mock": dry_run,
            "provenance": {"inference_started": not dry_run},
        }

    def admet(self, *, candidate, output_dir, dry_run):
        calls.append(("ADMET", candidate["candidate_id"], dry_run))
        return {
            "status": "success",
            "metrics": {"hERG": 0.1, "AMES": 0.1},
            "is_mock": dry_run,
            "provenance": {"inference_started": not dry_run},
        }

    monkeypatch.setattr(GNINATool, "run", docking)
    monkeypatch.setattr(Boltz2Tool, "run", boltz)
    monkeypatch.setattr(ADMETTool, "run", admet)
    return calls


def test_small_molecule_full_chain_and_evidence(sm_state, fake_sm_models, monkeypatch):
    monkeypatch.setattr(
        ProtenixTool, "run", Mock(side_effect=AssertionError("No Protenix on molecules"))
    )
    update = screening.screening_node(sm_state)
    result = update["screening_result"]
    assert [x[0] for x in fake_sm_models] == ["GNINA", "Boltz-2", "ADMET"]
    assert result["summary"] == {"pass": 1, "review": 0, "fail": 0}
    assert result["status"] == "COMPLETED" and "use_rfd3" not in update
    assert all(
        result[key] == 1 for key in ("gnina_completed", "boltz2_completed", "admet_completed")
    )
    c = result["candidates"][0]
    assert c["cross_model_validation"]["binding_site_agreement"] == 1
    assert c["cross_model_validation"]["pose_agreement"] is None
    assert c["admet_screening"]["verdict"] == "PASS"
    assert {c["source"] for c in update["evidence_cards"]} >= {
        "GNINA",
        "Boltz-2",
        "ADMET",
        "Screening",
    }
    assert result["profiles"] == {} and result["evidence_store"]["status"] == "SAVED"


def test_small_molecule_mock_full_chain_is_never_real_success(sm_state, fake_sm_models):
    sm_state["dry_run"] = True
    result = screening.screening_node(sm_state)["screening_result"]
    assert [x[0] for x in fake_sm_models] == ["GNINA", "Boltz-2", "ADMET"]
    assert result["execution_mode"] == "mock" and result["status"] == "NOT_EXECUTED"
    assert (
        result["summary"] == {"pass": 0, "review": 0, "fail": 0}
        and result["accepted_candidate_ids"] == []
    )
    assert result["mock_summary"]["pass"] == 1
    assert result["gnina_completed"] == result["boltz2_completed"] == result["admet_completed"] == 0


def test_binder_mock_never_calls_secondary(state, fake_models, monkeypatch):
    state["dry_run"] = True
    call = Mock(side_effect=AssertionError("Protenix is excluded"))
    monkeypatch.setattr(ProtenixTool, "run", call)
    result = screening.screening_node(state)["screening_result"]
    call.assert_not_called()
    assert result["summary"]["pass"] == 0
    assert result["protenix_completed"] == 0


@pytest.mark.parametrize("iptm,verdict", [(0.2, "FAIL"), (0.5, "REVIEW"), (0.85, "PASS")])
def test_primary_af3_gate_without_secondary(state, fake_models, monkeypatch, iptm, verdict):
    original = AlphaFold3BinderAdapter.run
    def af3(self, **kwargs):
        output = original(self, **kwargs)
        output.metrics["iptm"] = iptm
        return output
    monkeypatch.setattr(AlphaFold3BinderAdapter, "run", af3)
    monkeypatch.setattr(ProtenixTool, "run", Mock(side_effect=AssertionError("excluded")))
    result = screening.screening_node(state)["screening_result"]
    assert all(c["final_verdict"] == verdict for c in result["candidates"])
    assert result["protenix_tool_calls"] == 0


@pytest.mark.parametrize("budget", [0, 1, 10])
def test_protenix_budget_ignored(state, fake_models, monkeypatch, budget):
    state["screening_options"] = {"max_protenix_candidates": budget}
    call = Mock(side_effect=RuntimeError("unavailable"))
    monkeypatch.setattr(ProtenixTool, "run", call)
    update = screening.screening_node(state)
    result = update["screening_result"]
    assert result["status"] == "COMPLETED" and result["summary"]["pass"] == 4
    assert critic_node({**state, **update})["critic_verdict"] == "APPROVE"
    call.assert_not_called()


def test_secondary_rigid_transform_is_agreement(tmp_path):
    a = tmp_path / "a.cif"
    b = tmp_path / "b.cif"
    a.write_text(complex_cif())
    lines = []
    for line in complex_cif().splitlines():
        if line.startswith("ATOM"):
            fields = line.split()
            x, y, z = map(float, fields[7:10])
            fields[7:10] = map(str, (-y + 21, x - 17, z + 4))
            line = " ".join(fields)
        lines.append(line)
    b.write_text("\n".join(lines) + "\n")
    assert BinderCrossModelGate._binder_rmsd(a, b) < 1e-8


def test_excluded_secondary_has_no_effect(state, fake_models, monkeypatch):
    monkeypatch.setattr(ProtenixTool, "run", Mock(side_effect=AssertionError("excluded")))
    result = screening.screening_node(state)["screening_result"]
    assert result["summary"]["pass"] == 4
    assert all(c["cross_validation_status"] == "not_used" for c in result["candidates"])


@pytest.mark.parametrize(
    "failure_stage,expected_status,expected_calls",
    [
        ("GNINA", "FAILED_GNINA", ["GNINA"]),
        ("Boltz-2", "FAILED_BOLTZ2", ["GNINA", "Boltz-2"]),
        ("ADMET", "FAILED_ADMET", ["GNINA", "Boltz-2", "ADMET"]),
    ],
)
def test_small_molecule_failure_stops_only_downstream(
    sm_state, fake_sm_models, monkeypatch, failure_stage, expected_status, expected_calls
):
    cls = {"GNINA": GNINATool, "Boltz-2": Boltz2Tool, "ADMET": ADMETTool}[failure_stage]

    def failed(self, **kwargs):
        fake_sm_models.append((failure_stage, "failed", False))
        raise RuntimeError(failure_stage + " failed")

    monkeypatch.setattr(cls, "run", failed)
    result = screening.screening_node(sm_state)["screening_result"]
    assert [x[0] for x in fake_sm_models] == expected_calls
    assert result["candidates"][0]["status"] == expected_status
    assert result["summary"]["pass"] == 0


def test_pose_failure_skips_boltz(sm_state, fake_sm_models, monkeypatch):
    original = GNINATool.run

    def outside(self, **kwargs):
        result = original(self, **kwargs)
        write_ligand(Path(result["pose_path"]), offset=100.0)
        return result

    monkeypatch.setattr(GNINATool, "run", outside)
    result = screening.screening_node(sm_state)["screening_result"]
    assert [x[0] for x in fake_sm_models] == ["GNINA"]
    assert result["summary"]["fail"] == 1


def test_missing_binding_site_does_not_launch(sm_state, monkeypatch):
    path = Path(sm_state["design_results_path"])
    design = json.loads(path.read_text())
    del design["designs"][0]["docking_box"]
    path.write_text(json.dumps(design))
    monkeypatch.setattr(GNINATool, "_execute", Mock(side_effect=AssertionError("must not launch")))
    result = screening.screening_node(sm_state)["screening_result"]
    assert result["status"] == "NOT_CONFIGURED"
    assert "MISSING_BINDING_SITE" in result["candidates"][0]["failure_reasons"][0]


def test_admet_warning_is_separate_and_blocks_pass(sm_state, fake_sm_models, monkeypatch):
    monkeypatch.setattr(
        ADMETTool,
        "run",
        lambda self, **kw: {
            "status": "success",
            "metrics": {"hERG": 0.95, "AMES": 0.1},
            "is_mock": False,
        },
    )
    result = screening.screening_node(sm_state)["screening_result"]
    c = result["candidates"][0]
    assert (
        c["cross_model_validation"]["verdict"] == "PASS"
        and c["admet_screening"]["verdict"] == "FAIL"
    )
    assert c["final_verdict"] == "FAIL" and "ADMET_LIMIT_EXCEEDED: hERG" in c["failure_reasons"]


def test_missing_admet_policy_is_review(sm_state, fake_sm_models):
    sm_state["screening_options"] = {}
    result = screening.screening_node(sm_state)["screening_result"]
    assert result["candidates"][0]["failure_reasons"] == ["ADMET_THRESHOLDS_NOT_CONFIGURED"]


def test_small_molecule_top_pose_budget(sm_state, fake_sm_models):
    path = Path(sm_state["design_results_path"])
    payload = json.loads(path.read_text())
    payload["designs"][0]["candidates"] *= 3
    path.write_text(json.dumps(payload))
    sm_state["screening_options"].update(
        max_gnina_candidates=2, max_boltz2_candidates=1, max_admet_candidates=1
    )
    result = screening.screening_node(sm_state)["screening_result"]
    assert [x[0] for x in fake_sm_models] == ["GNINA", "GNINA", "Boltz-2", "ADMET"]
    assert result["gnina_completed"] == 2 and result["summary"] == {
        "pass": 1,
        "review": 2,
        "fail": 0,
    }


def test_mixed_branches_have_one_evidence_store_and_global_target_cap(
    state, sm_state, fake_models, fake_sm_models
):
    path = Path(state["design_results_path"])
    payload = json.loads(path.read_text())
    payload["designs"] += json.loads(Path(sm_state["design_results_path"]).read_text())["designs"]
    path.write_text(json.dumps(payload))
    state["advance_targets"] += sm_state["advance_targets"]
    state["screening_options"] = {**sm_state["screening_options"], "max_targets": 2}
    result = screening.screening_node(state)["screening_result"]
    assert result["modality"] == "MIXED" and len(result["candidates"]) == 5
    assert result["summary"]["pass"] == 5
    assert result["protein_design_result"]["candidate_count"] == 4
    state["budget"] = {"max_candidates": 4}
    fake_sm_models.clear()
    result = screening.screening_node(state)["screening_result"]
    assert fake_sm_models == [] and result["gnina_completed"] == 0


@pytest.mark.parametrize("modality", ["ADC", "RNA_THERAPEUTIC", "BIOMARKER"])
def test_other_modalities_are_not_misrouted(modality, monkeypatch):
    monkeypatch.setattr(GNINATool, "run", Mock(side_effect=AssertionError("wrong modality")))
    with pytest.raises(ValueError, match="Unsupported modality"):
        screening.screening_node({"advance_targets": [{"gene_name": "X", "modality": modality}]})


@pytest.mark.parametrize("cls", [ProtenixTool, GNINATool, Boltz2Tool, ADMETTool])
def test_unconfigured_tools_never_invent_a_command(cls):
    assert cls().is_configured()[0] is False


def test_native_protenix_and_boltz_inputs_are_independent(tmp_path):
    protenix = ProtenixTool({"executable": "/missing/protenix"})
    result = protenix.run(
        target_sequence="A" * 8,
        binder_sequence="G" * 8,
        output_dir=str(tmp_path / "p"),
        dry_run=True,
    )
    data = json.loads(Path(result["input_json"]).read_text())[0]
    assert data["sequences"][0]["proteinChain"]["id"] == ["A"]
    assert "templatesPath" not in str(data) and "pocket" not in data
    assert (
        result["provenance"]["command"][result["provenance"]["command"].index("-n") + 1]
        == "protenix-v2"
    )
    boltz = Boltz2Tool({"executable": "/missing/boltz"})
    result = boltz.run(
        target_sequence="A" * 8, smiles="CCO", output_dir=str(tmp_path / "b"), dry_run=True
    )
    import yaml

    data = yaml.safe_load(Path(result["input_yaml"]).read_text())
    assert data["properties"] == [{"affinity": {"binder": "B"}}] and "templates" not in data
    assert result["is_mock"] and not result["provenance"]["inference_started"]


def test_protenix_and_boltz_native_parsers(tmp_path):
    protenix = tmp_path / "p"
    protenix.mkdir()
    (protenix / "binder_11_sample_0.cif").write_text(complex_cif())
    (protenix / "binder_11_summary_confidence_sample_0.json").write_text(
        json.dumps({"iptm": 0.8, "ranking_score": 0.9, "has_clash": False})
    )
    parsed = ProtenixTool.parse_output(protenix)
    assert parsed["iptm"] == 0.8 and "ptm" not in parsed
    boltz = tmp_path / "b"
    boltz.mkdir()
    ligand_complex(boltz / "input_model_0.cif")
    (boltz / "confidence_input_model_0.json").write_text(
        json.dumps({"ligand_iptm": 0.7, "confidence_score": 0.8})
    )
    (boltz / "affinity_input.json").write_text(
        json.dumps({"affinity_pred_value": -1.2, "affinity_probability_binary": 0.8})
    )
    parsed = Boltz2Tool.parse_output(boltz)
    assert parsed["affinity_pred_value"] == -1.2 and "has_clash" not in parsed


def test_existing_admet_output_parser_rejects_mock(tmp_path):
    manifest = {"backend": "admet_ai", "is_mock": False}
    (tmp_path / "run_manifest.json").write_text(json.dumps(manifest))
    (tmp_path / "admet_predictions_long.csv").write_text(
        "candidate_id,endpoint_name,endpoint_value,prediction_status,is_mock\nligand,AMES,0.1,success,false\nother,hERG,0.9,success,false\n"
    )
    result = ADMETTool.parse_output(tmp_path, "ligand")
    assert result["metrics"] == {"AMES": 0.1}
    manifest["is_mock"] = True
    (tmp_path / "run_manifest.json").write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match="not a real"):
        ADMETTool.parse_output(tmp_path, "ligand")


def test_gnina_preparation_reuses_validation_and_preserves_identity(tmp_path):
    path, provenance = GNINATool.prepare_ligand(
        {"candidate_id": "x", "smiles": "CCO"}, tmp_path, 11
    )
    assert path.is_file() and provenance["canonical_smiles"] == "CCO"
    assert "as_supplied" in provenance["protonation"]
    with pytest.raises(ValueError, match="single-component"):
        GNINATool.prepare_ligand({"candidate_id": "x", "smiles": "CCO.[Na+]"}, tmp_path, 11)


@pytest.mark.parametrize(
    "options",
    [
        {"max_protenix_candidates": -1},
        {"max_gnina_candidates": True},
        {"cascade_thresholds": {"gnina_cnn_min": 2}},
        {"cascade_thresholds": {"admet_limits": {"x": {"min": 1, "max": 0}}}},
        {"boltz2_options": {"timeout_seconds": 0}},
        {"gnina_options": {"unknown_flag": 1}},
    ],
)
def test_invalid_new_configuration(options):
    with pytest.raises((ValueError, TypeError)):
        ScreeningOptions.from_state({"screening_options": options})


def test_explicit_small_molecule_modality_wins_over_stale_binder_flag(
    sm_state, fake_sm_models, monkeypatch
):
    sm_state["use_rfd3"] = True
    monkeypatch.setattr(
        screening, "_screen_binders", Mock(side_effect=AssertionError("stale binder routing"))
    )
    result = screening.screening_node(sm_state)["screening_result"]
    assert result["modality"] == "SMALL_MOLECULE" and result["summary"]["pass"] == 1


def test_compressed_cross_model_structures(tmp_path):
    import gzip

    import zstandard

    from ..screening_gates import sequence_for_chain

    a = tmp_path / "a.cif.gz"
    b = tmp_path / "b.cif.zst"
    a.write_bytes(gzip.compress(complex_cif().encode()))
    b.write_bytes(zstandard.ZstdCompressor().compress(complex_cif().encode()))
    assert sequence_for_chain(a, "A") == "A" * 8 and sequence_for_chain(b, "B") == "G" * 8
    assert BinderCrossModelGate._binder_rmsd(a, b) < 1e-8


def test_boltz_wrong_ligand_cannot_pass(sm_state, fake_sm_models, monkeypatch):
    original = Boltz2Tool.run

    def wrong_ligand(self, **kwargs):
        result = original(self, **kwargs)
        path = Path(result["metrics"]["model_path"])
        path.write_text(path.read_text().replace(" O O1 LIG", " N N1 LIG"))
        return result

    monkeypatch.setattr(Boltz2Tool, "run", wrong_ligand)
    result = screening.screening_node(sm_state)["screening_result"]
    assert result["summary"]["pass"] == 0
    assert (
        "BOLTZ2_LIGAND_ELEMENT_COMPOSITION_MISMATCH" in result["candidates"][0]["failure_reasons"]
    )
