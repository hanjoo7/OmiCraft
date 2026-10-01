"""Screening-only tests: all inference/network calls are forbidden or stubbed."""

import json
import socket
import subprocess
from dataclasses import replace
from pathlib import Path
from unittest.mock import Mock

import pytest

from .. import screening_agent as screening
from ..configuration import DataPaths, OmiCraftConfig, get_config, set_config
from ..evidence_store import EvidenceStore
from ..protein_design_route import AlphaFold3BinderAdapter, AlphaFold3ValidationOutput
from ..protein_mpnn_tool import ProteinMPNNTool
from ..rfdiffusion3_stage import RFdiffusion3Output, RFdiffusion3Runner, _mock_cif
from ..screening_model_tools import ProtenixTool
from ..screening_tools import StructuralScreeningTool, StructuralThresholds


@pytest.fixture(autouse=True)
def isolated(tmp_path, monkeypatch, request):
    old = get_config()
    set_config(OmiCraftConfig(data=DataPaths(base_dir=str(tmp_path / "data"))))

    def forbidden(*args, **kwargs):
        raise AssertionError("Actual inference and network are forbidden")

    if not any(request.node.get_closest_marker(m) for m in ("external_tool", "local_io")):
        monkeypatch.setattr(subprocess, "run", forbidden)
        monkeypatch.setattr(subprocess, "Popen", forbidden)
        monkeypatch.setattr(socket.socket, "connect", forbidden)
    yield
    set_config(old)


def complex_cif(distance=3.0):
    header = _mock_cif(1).split("ATOM 1")[0]
    rows = []
    for chain, y in (("A", 0), ("B", distance)):
        for n in range(1, 9):
            rows.append(
                f"ATOM {len(rows) + 1} C CA {'ALA' if chain == 'A' else 'GLY'} {chain} {n} {n * 3.8} {y + (n % 3) * 0.4} {(n % 2) * 0.3}\n"
            )
    return header + "".join(rows) + "#\n"


def af3_fixture(
    directory,
    *,
    distance=3.0,
    iptm=0.85,
    pair_iptm=0.82,
    pair_pae=4.0,
    plddt=85.0,
    clash=0.0,
    prefix="sample_",
):
    directory.mkdir(parents=True, exist_ok=True)
    (directory / f"{prefix}model.cif").write_text(complex_cif(distance))
    summary = {
        "ptm": 0.8,
        "iptm": iptm,
        "ranking_score": 0.84,
        "has_clash": clash,
        "chain_ids": ["A", "B"],
        "chain_pair_iptm": [[0.8, pair_iptm], [pair_iptm, 0.8]],
        "chain_pair_pae_min": [[0, pair_pae], [pair_pae, 0]],
        "chain_ptm": [0.8, 0.8],
        "chain_iptm": [0.85, 0.85],
    }
    (directory / f"{prefix}summary_confidences.json").write_text(json.dumps(summary))
    (directory / f"{prefix}confidences.json").write_text(
        json.dumps(
            {
                "atom_plddts": [plddt] * 16,
                "atom_chain_ids": ["A"] * 8 + ["B"] * 8,
                "token_chain_ids": ["A"] * 8 + ["B"] * 8,
            }
        )
    )
    return {
        "status": "success",
        "output_dir": str(directory),
        "metrics": AlphaFold3BinderAdapter.parse_output(directory),
    }


@pytest.fixture
def state(tmp_path):
    target = tmp_path / "target.cif"
    target.write_text(_mock_cif(8, "T"))
    handoff = tmp_path / "design_output.json"
    handoff.write_text(
        json.dumps(
            {
                "n_designed": 1,
                "designs": [
                    {
                        "gene_name": "TEST",
                        "modality": "DE_NOVO_BINDER",
                        "target_structure": str(target),
                        "target_sequence": "A" * 8,
                        "target_chain": "T",
                        "hotspots": ["T3"],
                        "binder_length": 8,
                        "num_designs": 4,
                        "design_steps": [
                            {"step": "proteinmpnn_sequence", "tools": ["ProteinMPNN"]}
                        ],
                    }
                ],
            }
        )
    )
    return {
        "advance_targets": [{"gene_name": "TEST", "modality": "DE_NOVO_BINDER"}],
        "design_results_path": str(handoff),
        "dry_run": False,
        "use_rfd3": False,
    }


@pytest.fixture
def fake_models(monkeypatch):
    generated_requests, af3_calls = [], []

    def generate(self, request, *, dry_run=False):
        generated_requests.append(request)
        directory = Path(request.output_dir)
        directory.mkdir(parents=True)
        paths = []
        for index in range(request.num_designs):
            path = directory / f"candidate_{index}.cif"
            path.write_text(_mock_cif(request.design_length, request.binder_chain))
            paths.append(str(path))
        return RFdiffusion3Output(
            status="dry_run" if dry_run else "success",
            generated_structures=tuple(paths),
            metadata={"is_mock": dry_run},
        )

    def sequence_design(self, *, backbone, output_dir, num_sequences, dry_run):
        path = Path(output_dir) / "seqs.fa"
        path.parent.mkdir(parents=True, exist_ok=True)
        sequence = "G" * len(backbone["rfd3_sequence"])
        path.write_text(
            ">backbone, designed_chains=['A'], seed=11\n"
            + backbone["rfd3_sequence"]
            + "\n"
            + "".join(
                f">T=0.1, sample={index + 1}, score=0.8, global_score=0.9\n{sequence}\n"
                for index in range(num_sequences)
            )
        )
        results = self.parse_output(
            path,
            candidate_id=backbone["candidate_id"],
            expected_length=len(sequence),
            binder_chain="A",
            max_sequences=num_sequences,
        )
        for result in results:
            result["is_mock"] = dry_run
        return {
            "status": "dry_run" if dry_run else "success",
            "is_mock": dry_run,
            "sequences": results,
        }

    def predict(self, *, target_sequence, binder_sequence, output_dir, seeds=(11,), dry_run=False):
        af3_calls.append({"sequence": binder_sequence, "seeds": seeds, "dry_run": dry_run})
        path = self.prepare_input(
            target_sequence=target_sequence,
            binder_sequence=binder_sequence,
            output_dir=output_dir,
            seeds=seeds,
        )
        result = af3_fixture(Path(output_dir))
        return AlphaFold3ValidationOutput(
            status="dry_run" if dry_run else "success",
            input_json=str(path),
            output_dir=output_dir,
            metrics=result["metrics"],
        )

    monkeypatch.setattr(ProtenixTool, "run", Mock(side_effect=AssertionError("Protenix excluded")))
    monkeypatch.setattr(ProteinMPNNTool, "run", sequence_design)
    monkeypatch.setattr(RFdiffusion3Runner, "run", generate)
    monkeypatch.setattr(AlphaFold3BinderAdapter, "run", predict)
    return generated_requests, af3_calls


def context_from(state):
    design = json.loads(Path(state["design_results_path"]).read_text())["designs"][0]
    return {
        "target_structure": design["target_structure"],
        "target_sequence": design["target_sequence"],
        "target_chain": "T",
        "hotspot_residues": ("T3",),
    }


@pytest.mark.parametrize(
    "changes,verdict,reason",
    [
        ({}, "PASS", None),
        ({"iptm": 0.5}, "REVIEW", "IPTM_BORDERLINE"),
        ({"iptm": 0.2}, "FAIL", "IPTM_BELOW_HARD_LIMIT"),
        ({"pair_iptm": 0.2}, "FAIL", "PAIR_IPTM_BELOW_HARD_LIMIT"),
        ({"pair_pae": 25.0}, "FAIL", "PAIR_PAE_ABOVE_HARD_LIMIT"),
        ({"pair_pae": 15.0}, "REVIEW", "PAIR_PAE_BORDERLINE"),
        ({"plddt": 40.0}, "FAIL", "BINDER_PLDDT_BELOW_HARD_LIMIT"),
        ({"clash": 1.0}, "FAIL", "AF3_HAS_CLASH"),
        ({"distance": 0.5}, "FAIL", "SEVERE_INTERCHAIN_CLASH"),
        ({"distance": 30.0}, "FAIL", "INSUFFICIENT_INTERFACE_CONTACTS"),
    ],
)
def test_structural_verdicts(tmp_path, state, changes, verdict, reason):
    result = StructuralScreeningTool().run(
        af3_fixture(tmp_path / "af3", **changes), context_from(state)
    )
    assert result["verdict"] == verdict
    if reason:
        assert reason in result["failure_reasons"]


def test_hotspot_contact_true_false_unknown(tmp_path, state):
    tool = StructuralScreeningTool()
    context = context_from(state)
    assert tool.run(af3_fixture(tmp_path / "near"), context)["hotspot_preserved"] is True
    assert (
        tool.run(af3_fixture(tmp_path / "far", distance=30), context)["hotspot_preserved"] is False
    )
    unknown = tool.run(af3_fixture(tmp_path / "unknown"), {**context, "target_sequence": "GGGG"})
    assert unknown["hotspot_preserved"] is None and unknown["verdict"] == "REVIEW"

    # Missing mapping must never erase a known hard geometry failure.
    clash = tool.run(
        af3_fixture(tmp_path / "clash", distance=0.5), {**context, "target_structure": "/missing"}
    )
    assert clash["verdict"] == "FAIL" and "SEVERE_INTERCHAIN_CLASH" in clash["failure_reasons"]


def test_thresholds_are_configurable(tmp_path, state):
    af3 = af3_fixture(tmp_path / "af3", iptm=0.55)
    assert StructuralScreeningTool().run(af3, context_from(state))["verdict"] == "REVIEW"
    thresholds = replace(StructuralThresholds(), iptm_pass_min=0.5)
    assert StructuralScreeningTool(thresholds).run(af3, context_from(state))["verdict"] == "PASS"
    with pytest.raises(ValueError):
        StructuralThresholds(iptm_fail_below=0.9, iptm_pass_min=0.4)


def test_missing_metrics_review_not_pass(tmp_path, state):
    af3 = af3_fixture(tmp_path / "af3")
    af3["metrics"]["chain_pair_iptm"] = None
    af3["metrics"]["ptm"] = float("nan")
    result = StructuralScreeningTool().run(af3, context_from(state))
    assert result["verdict"] == "REVIEW"


def test_actual_parser_numeric_clash_unprefixed_files(tmp_path):
    result = af3_fixture(tmp_path / "af3", prefix="")
    metrics = result["metrics"]
    assert metrics["has_clash"] is False
    assert metrics["chain_pair_pae_min"][0][1] == 4.0
    assert metrics["chain_pair_iptm"][0][1] == 0.82
    assert metrics["plddt_summary"]["by_chain_mean"]["B"] == 85.0
    assert metrics["model_path"].endswith("/model.cif")


def test_chain_pair_order_and_no_invented_mapping(tmp_path, state):
    af3 = af3_fixture(tmp_path / "af3")
    af3["metrics"]["chain_ids"] = ["B", "A"]
    assert StructuralScreeningTool().run(af3, context_from(state))["verdict"] == "PASS"
    af3["metrics"]["chain_ids"] = []
    assert StructuralScreeningTool().run(af3, context_from(state))["verdict"] == "REVIEW"


def test_routing_budget_counts_and_evidence(state, fake_models):
    state["screening_options"] = {"max_rfd3_candidates": 3, "max_af3_candidates": 2}
    update = screening.screening_node(state)
    result = update["screening_result"]
    assert result["pipeline"] == [
        "RFdiffusion3",
        "ProteinMPNN",
        "AlphaFold3",
        "structural_screening",
        "Validation",
        "Critic",
    ]
    assert result["modality"] == "DE_NOVO_BINDER"
    assert result["requested_designs"] == 4 and result["generated_designs"] == 3
    assert result["af3_evaluated"] == 2 and result["summary"] == {"pass": 2, "review": 1, "fail": 0}
    assert result["candidates"][-1]["status"] == "SKIPPED_RESOURCE"
    requests, calls = fake_models
    assert requests[0].num_designs == 3 and requests[0].hotspot_residues == ("T3",)
    assert len(calls) == 2 and all(c["sequence"] == "G" * 8 for c in calls)
    assert update["use_rfd3"] is True  # Existing Critic can consume the result.
    store = EvidenceStore(result["evidence_store"]["path"])
    try:
        rows = store.get_by_run(result["provenance"]["run_id"])
        assert len(rows) == len(update["evidence_cards"])
        payloads = [json.loads(row["content"]) for row in rows]
        assert any(p.get("metric") == "chain_pair_iptm" for p in payloads)
    finally:
        store.close()
    for candidate in result["candidates"][:2]:
        payload = json.loads(Path(candidate["af3_validation"]["input_json"]).read_text())
        assert payload["sequences"][1]["protein"]["sequence"] == candidate["sequence"]
        assert all("templates" not in entity["protein"] for entity in payload["sequences"])


def test_af3_one_failure_does_not_stop_remaining(state, fake_models, monkeypatch):
    original = AlphaFold3BinderAdapter.run
    count = 0

    def one_failure(self, **kwargs):
        nonlocal count
        count += 1
        if count == 2:
            raise RuntimeError("fixture AF3 failure")
        return original(self, **kwargs)

    monkeypatch.setattr(AlphaFold3BinderAdapter, "run", one_failure)
    result = screening.screening_node(state)["screening_result"]
    assert count == 4 and result["af3_evaluated"] == 3
    assert result["candidates"][1]["status"] == "FAILED_AF3"
    assert result["summary"] == {"pass": 3, "review": 0, "fail": 1}


def test_dry_run_counts_are_not_simulated(state, fake_models):
    state["dry_run"] = True
    result = screening.screening_node(state)["screening_result"]
    assert result["execution_mode"] == "mock" and result["dry_run"] is True
    assert (
        result["generated_designs"] == result["af3_evaluated"] == result["final_binder_hits"] == 0
    )
    assert result["summary"] == {"pass": 0, "review": 0, "fail": 0}
    assert all(c["is_mock"] and c["verdict"] == "REVIEW" for c in result["candidates"])


def test_resource_skip_no_tools_or_candidates(state, monkeypatch):
    generate = Mock(side_effect=AssertionError("RFD3 must not run"))
    predict = Mock(side_effect=AssertionError("AF3 must not run"))
    monkeypatch.setattr(
        ProteinMPNNTool, "run", Mock(side_effect=AssertionError("MPNN must not run"))
    )
    monkeypatch.setattr(RFdiffusion3Runner, "run", generate)
    monkeypatch.setattr(AlphaFold3BinderAdapter, "run", predict)
    state["screening_options"] = {"resource_profile": "C"}
    result = screening.screening_node(state)["screening_result"]
    assert result["status"] == "SKIPPED_RESOURCE"
    assert not result["candidates"] and result["summary"]["pass"] == 0
    generate.assert_not_called()
    predict.assert_not_called()


def test_no_rfd3_dependency_not_configured(state, monkeypatch):
    monkeypatch.setattr(RFdiffusion3Runner, "available", lambda self: False)
    result = screening.screening_node(state)["screening_result"]
    assert result["status"] == "NOT_CONFIGURED" and not result["candidates"]
    assert result["summary"]["pass"] == 0


def test_design_placeholder_does_not_become_hotspot(state, fake_models):
    path = Path(state["design_results_path"])
    data = json.loads(path.read_text())
    data["designs"][0]["hotspots"] = "qualification 단계에서 도출된 결합 잔기"
    path.write_text(json.dumps(data))
    result = screening.screening_node(state)["screening_result"]
    assert all(c["verdict"] == "REVIEW" for c in result["candidates"])
    assert fake_models[0][0].hotspot_residues == ()


def test_mixed_modality_preserves_small_branch(state, fake_models):
    state["advance_targets"].append({"gene_name": "SM", "modality": "SMALL_MOLECULE"})
    result = screening.screening_node(state)["screening_result"]
    assert result["profiles"] == {}
    assert result["skipped_targets"][0]["modality"] == "SMALL_MOLECULE"
    assert result["generated_designs"] == 4


def test_invalid_screening_config_does_not_launch(state, monkeypatch):
    monkeypatch.setattr(
        RFdiffusion3Runner, "run", Mock(side_effect=AssertionError("must not launch"))
    )
    state["screening_options"] = {"max_rfd3_candidates": -1}
    result = screening.screening_node(state)["screening_result"]
    assert result["status"] == "NOT_CONFIGURED" and not result["candidates"]


def test_existing_primary_output_reused(state, fake_models, monkeypatch):
    initial = screening.screening_node(state)
    monkeypatch.setattr(RFdiffusion3Runner, "run", Mock(side_effect=AssertionError("RFD3 rerun")))
    monkeypatch.setattr(
        AlphaFold3BinderAdapter, "run", Mock(side_effect=AssertionError("AF3 rerun"))
    )
    result = screening.screening_node(
        {**state, "protein_design_result": initial["protein_design_result"]}
    )["screening_result"]
    assert result["af3_attempted"] == 0 and result["af3_evaluated"] == 4
    assert result["summary"]["pass"] == 4


def test_failed_af3_with_matching_input_is_retried(state, fake_models, monkeypatch):
    state["screening_options"] = {
        "max_rfd3_candidates": 1,
        "max_mpnn_candidates": 1,
        "max_af3_candidates": 1,
    }
    initial = screening.screening_node(state)
    previous = initial["protein_design_result"]
    candidate = previous["candidates"][0]
    candidate["af3_validation"].update(status="failed", metrics={}, error_message="HMMER missing")
    monkeypatch.setattr(RFdiffusion3Runner, "run", Mock(side_effect=AssertionError("RFD3 rerun")))
    monkeypatch.setattr(ProteinMPNNTool, "run", Mock(side_effect=AssertionError("MPNN rerun")))
    af3_calls_before = len(fake_models[1])
    result = screening.screening_node({**state, "protein_design_result": previous})[
        "screening_result"
    ]
    assert len(fake_models[1]) == af3_calls_before + 1
    assert result["af3_attempted"] == 1 and result["af3_evaluated"] == 1
    assert result["mpnn_attempted"] == 0
    assert result["summary"] == {"pass": 1, "review": 0, "fail": 0}


def test_global_budget_across_targets(state, fake_models):
    path = Path(state["design_results_path"])
    data = json.loads(path.read_text())
    data["designs"].append({**data["designs"][0], "gene_name": "SECOND"})
    path.write_text(json.dumps(data))
    state["advance_targets"].append({"gene_name": "SECOND", "modality": "DE_NOVO_BINDER"})
    state["screening_options"] = {
        "max_targets": 2,
        "max_rfd3_candidates": 5,
        "max_af3_candidates": 3,
    }
    result = screening.screening_node(state)["screening_result"]
    assert [r.num_designs for r in fake_models[0]] == [4, 1]
    assert result["requested_designs"] == 8 and result["generated_designs"] == 5
    assert result["af3_evaluated"] == 3 and result["summary"] == {"pass": 3, "review": 2, "fail": 0}


def test_compressed_af3_artifacts(tmp_path, state):
    import zstandard

    directory = tmp_path / "af3"
    af3_fixture(directory)
    for name in ("sample_model.cif", "sample_confidences.json"):
        path = directory / name
        path.with_suffix(path.suffix + ".zst").write_bytes(
            zstandard.ZstdCompressor().compress(path.read_bytes())
        )
        path.unlink()
    metrics = AlphaFold3BinderAdapter.parse_output(directory)
    assert metrics["model_path"].endswith(".zst")
    result = StructuralScreeningTool().run(
        {"status": "success", "metrics": metrics}, context_from(state)
    )
    assert result["verdict"] == "PASS" and result["hotspot_preserved"] is True


def test_invalid_handoff_is_not_configured(state, monkeypatch):
    Path(state["design_results_path"]).write_text("{ broken JSON")
    monkeypatch.setattr(RFdiffusion3Runner, "run", Mock(side_effect=AssertionError("Must not run")))
    result = screening.screening_node(state)["screening_result"]
    assert result["status"] == "NOT_CONFIGURED" and not result["candidates"]


def test_model_artifact_missing_is_failed_af3(state, fake_models, monkeypatch):
    def missing(self, **kwargs):
        return AlphaFold3ValidationOutput(
            status="success", metrics={"model_paths": ["/missing/model.cif"]}
        )

    monkeypatch.setattr(AlphaFold3BinderAdapter, "run", missing)
    result = screening.screening_node(state)["screening_result"]
    assert result["af3_evaluated"] == 0 and result["summary"]["pass"] == 0
    assert all(c["status"] == "FAILED_AF3" for c in result["candidates"])


def test_cli_demo_no_inference(tmp_path, monkeypatch):
    import sys

    monkeypatch.setattr(
        sys, "argv", ["screening_agent", "--demo", "--output-dir", str(tmp_path / "demo")]
    )
    assert screening.main() == 0
    output = json.loads((tmp_path / "demo" / "agent_results" / "screening_output.json").read_text())
    assert output["status"] == "NOT_EXECUTED" and output["generated_designs"] == 0
    assert output["af3_attempted"] == 0 and output["summary"]["pass"] == 0


def test_dossier_binder_fallback_never_enters_small_molecule_simulation(monkeypatch):
    directory = Path(get_config().data.agent_results)
    directory.mkdir(parents=True)
    (directory / "final_dossier.json").write_text(
        json.dumps({"advance_targets": [{"gene_name": "TEST", "modality": "DE_NOVO_BINDER"}]})
    )
    (directory / "design_output.json").write_text(
        json.dumps({"designs": [{"gene_name": "TEST", "modality": "SMALL_MOLECULE"}]})
    )
    monkeypatch.setattr(
        screening, "_screen_small_molecules", Mock(side_effect=AssertionError("Wrong modality"))
    )
    result = screening.screening_node({})["screening_result"]
    assert result["modality"] == "DE_NOVO_BINDER" and result["status"] == "NOT_CONFIGURED"
