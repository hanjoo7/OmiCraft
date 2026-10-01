"""Offline tests: official CLI is stubbed; GPU inference and network are forbidden."""

import json
import subprocess
import sys
from pathlib import Path
from unittest.mock import Mock

import pytest

from .. import screening_agent as screening
from ..configuration import ProteinMPNNConfig, get_config
from ..protein_mpnn_tool import ProteinMPNNTool
from ..rfdiffusion3_stage import RFdiffusion3Output, RFdiffusion3Runner
from ..screening_tools import AF3Tool, RFD3Tool, ScreeningOptions
from .test_screening_tools import af3_fixture

REAL_MPNN_RUN = ProteinMPNNTool.run


def write_fasta(path, sequence="G" * 8, count=1, designed="B"):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        f">backbone, score=0.01, fixed_chains=['A'], designed_chains=['{designed}'], model_name=v_48_020, seed=11\nAAAAAAAA\n"
        + "".join(
            f">T=0.1, sample={i + 1}, score=0.8, global_score=0.9\n{sequence}\n"
            for i in range(count)
        )
    )
    return path


@pytest.fixture
def adapter_input(tmp_path):

    # Synthetic coordinates for parser tests, never a model prediction.
    path = tmp_path / "source.pdb"
    rows = []
    for chain, offset in (("A", 0), ("B", 10)):
        for residue in range(8):
            for name, dx, element in (
                ("N", 0, "N"),
                ("CA", 1.2, "C"),
                ("C", 2.4, "C"),
                ("O", 2.9, "O"),
            ):
                rows.append(
                    f"ATOM  {len(rows) + 1:5d} {name:^4s} ALA {chain}{residue + 101:4d}    "
                    f"{residue * 3.8 + dx:8.3f}{offset:8.3f}{0:8.3f}{1:6.2f}{80:6.2f}          {element:>2s}\n"
                )
    path.write_text("".join(rows) + "END\n")
    script = tmp_path / "protein_mpnn_run.py"
    script.write_text("# Test placeholder; never executed.\n")
    weights = tmp_path / "weights"
    weights.mkdir()
    (weights / "v_48_020.pt").write_text("TEST_ONLY_NOT_WEIGHTS")
    config = ProteinMPNNConfig(
        python_path=sys.executable,
        script_path=str(script),
        model_weights_dir=str(weights),
        timeout_seconds=7,
    )
    backbone = {
        "candidate_id": "backbone_001",
        "rfd3_structure_path": str(path),
        "binder_chain": "B",
        "rfd3_sequence": "A" * 8,
        "is_mock": False,
    }
    return ProteinMPNNTool(config, gpu_device=1), backbone, tmp_path / "output"


def test_official_command_fixed_target_fasta_and_gpu(adapter_input, monkeypatch):
    tool, backbone, directory = adapter_input
    calls = []
    monkeypatch.setenv("CUDA_VISIBLE_DEVICES", "3,5")

    def execute(command, **kwargs):
        calls.append((command, kwargs))
        args = dict(zip(command[2::2], command[3::2]))
        assert args["--pdb_path_chains"] == "B"
        assert args["--batch_size"] == "1" and args["--num_seq_per_target"] == "2"
        assert args["--sampling_temp"] == "0.1" and args["--omit_AAs"] == "CX"
        assert kwargs["timeout"] == 7 and kwargs["env"]["CUDA_VISIBLE_DEVICES"] == "5"
        assert "shell" not in kwargs
        from Bio.PDB import PDBParser

        model = PDBParser(QUIET=True).get_structure("x", args["--pdb_path"])[0]
        assert [c.id for c in model] == ["A", "B"]  # Target retained as fixed context.
        assert len(model["A"]) == len(model["B"]) == 8
        write_fasta(Path(args["--out_folder"]) / "seqs" / "backbone.fa", count=2)
        return subprocess.CompletedProcess(command, 0)

    monkeypatch.setattr(subprocess, "run", execute)
    result = tool.run(backbone=backbone, output_dir=str(directory), num_sequences=2, dry_run=False)
    assert result["status"] == "success" and len(calls) == 1
    assert [s["sequence"] for s in result["sequences"]] == ["G" * 8, "G" * 8]
    assert all(s["score"] == 0.8 and s["status"] == "MPNN_COMPLETED" for s in result["sequences"])
    assert result["provenance"]["mapping"]["residue_mapping"]["B"][0] == {
        "source_residue": "101",
        "pdb_residue": 1,
    }
    assert result["provenance"]["script_sha256"] and result["provenance"]["inference_started"]
    assert json.loads((directory / "protein_mpnn_run.json").read_text())["status"] == "success"


def test_cif_export_preserves_label_chains_and_coordinates(adapter_input, tmp_path):
    from Bio.PDB import MMCIFIO, PDBParser
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict

    tool, backbone, directory = adapter_input
    writer = MMCIFIO()
    writer.set_structure(
        PDBParser(QUIET=True).get_structure("source", backbone["rfd3_structure_path"])
    )
    cif = tmp_path / "source.cif"
    writer.save(str(cif))
    data = MMCIF2Dict(str(cif))
    data["_atom_site.label_asym_id"] = [
        "LONG" if c == "B" else c for c in data["_atom_site.label_asym_id"]
    ]
    writer.set_dict(data)
    writer.save(str(cif))
    directory.mkdir()
    exported, mapping = tool.prepare_backbone(str(cif), "LONG", directory)
    assert mapping["chain_mapping"] == {"A": "A", "LONG": "B"}
    model = PDBParser(QUIET=True).get_structure("exported", exported)[0]
    assert list(model["B"][1]["CA"].coord) == pytest.approx([1.2, 10, 0])


@pytest.mark.parametrize(
    "change", ["missing_script", "missing_weights", "bad_model", "all_omitted"]
)
def test_not_configured_does_not_execute(adapter_input, change):
    tool, backbone, directory = adapter_input
    if change == "missing_script":
        tool.config.script_path = ""
    elif change == "missing_weights":
        tool.config.model_weights_dir = ""
    elif change == "bad_model":
        tool.config.model_name = "../model"
    else:
        tool.config.omit_AAs = "ACDEFGHIKLMNPQRSTVWY"
    result = tool.run(backbone=backbone, output_dir=str(directory), num_sequences=1, dry_run=False)
    assert result["execution_status"] == "NOT_CONFIGURED" and not result["sequences"]
    assert result["provenance"]["inference_started"] is False


def test_timeout_is_recorded(adapter_input, monkeypatch):
    tool, backbone, directory = adapter_input
    monkeypatch.delenv("CUDA_VISIBLE_DEVICES", raising=False)

    def timeout(command, **kwargs):
        raise subprocess.TimeoutExpired(command, kwargs["timeout"])

    monkeypatch.setattr(subprocess, "run", timeout)
    result = tool.run(backbone=backbone, output_dir=str(directory), num_sequences=1, dry_run=False)
    assert result["status"] == "failed" and "TimeoutExpired" in result["error_message"]


def test_stale_output_is_rejected(adapter_input):
    tool, backbone, directory = adapter_input
    write_fasta(directory / "mpnn_predictions" / "seqs" / "backbone.fa")
    result = tool.run(backbone=backbone, output_dir=str(directory), num_sequences=1, dry_run=False)
    assert result["status"] == "failed" and "must be empty" in result["error_message"]
    assert result["provenance"]["inference_started"] is False


def test_missing_atoms_rejected_before_mpnn(adapter_input):
    tool, backbone, directory = adapter_input
    path = Path(backbone["rfd3_structure_path"])
    path.write_text(
        "".join(
            line
            for line in path.read_text().splitlines(keepends=True)
            if line[12:16].strip() != "O"
        )
    )
    result = tool.run(backbone=backbone, output_dir=str(directory), num_sequences=1, dry_run=False)
    assert result["status"] == "failed" and "Incomplete" in result["error_message"]
    assert result["provenance"]["inference_started"] is False


def test_fasta_native_excluded_bad_sample_isolated(tmp_path):
    path = write_fasta(tmp_path / "result.fa", count=2)
    path.write_text(path.read_text().replace("sample=2, score=0.8", "sample=2, score=nan"))
    results = ProteinMPNNTool.parse_output(
        path, candidate_id="c", expected_length=8, binder_chain="B", max_sequences=2
    )
    assert len(results) == 2 and results[0]["score"] == 0.8
    assert results[1]["status"] == "FAILED_MPNN" and results[1]["score"] is None
    assert results[0]["sequence"] == "G" * 8  # Never the native A sequence.


@pytest.mark.parametrize("sequence", ["GGG", "G" * 7 + "X", "AAAA/GGGG"])
def test_fasta_invalid_binder_sequences(tmp_path, sequence):
    path = write_fasta(tmp_path / "result.fa", sequence=sequence)
    result = ProteinMPNNTool.parse_output(
        path, candidate_id="c", expected_length=8, binder_chain="B", max_sequences=1
    )[0]
    assert result["status"] == "FAILED_MPNN"


def test_fasta_wrong_designed_chain_rejected(tmp_path):
    with pytest.raises(ValueError, match="requested designed binder"):
        ProteinMPNNTool.parse_output(
            write_fasta(tmp_path / "x.fa", designed="A"),
            candidate_id="c",
            expected_length=8,
            binder_chain="B",
            max_sequences=1,
        )


def test_stage_order_and_mpnn_sequence_reaches_af3(state, fake_models, monkeypatch):
    events = []
    for cls, name in ((RFD3Tool, "rfd3"), (ProteinMPNNTool, "mpnn"), (AF3Tool, "af3")):
        original = cls.run

        def wrapper(self, *args, method=original, stage=name, **kwargs):
            events.append((stage, kwargs.get("backbone", {}).get("candidate_id")))
            return method(self, *args, **kwargs)

        monkeypatch.setattr(cls, "run", wrapper)
    result = screening.screening_node(state)["screening_result"]
    assert events[0][0] == "rfd3"
    assert [e[0] for e in events].count("mpnn") == 4
    assert next(i for i, e in enumerate(events) if e[0] == "mpnn") < next(
        i for i, e in enumerate(events) if e[0] == "af3"
    )
    assert result["mpnn_completed"] == result["af3_evaluated"] == 4
    assert all(
        c["rfd3_sequence"] == "A" * 8 and c["binder_sequence"] == "G" * 8
        for c in result["candidates"]
    )
    assert all(call["sequence"] == "G" * 8 for call in fake_models[1])


def test_mpnn_failure_isolated_and_af3_skipped(state, fake_models, monkeypatch):
    original = ProteinMPNNTool.run
    count = 0

    def fail_second(self, **kwargs):
        nonlocal count
        count += 1
        if count == 2:
            raise RuntimeError("Test MPNN error")
        return original(self, **kwargs)

    monkeypatch.setattr(ProteinMPNNTool, "run", fail_second)
    update = screening.screening_node(state)
    result = update["screening_result"]
    failed = result["candidates"][1]
    assert count == 4 and result["mpnn_completed"] == result["af3_evaluated"] == 3
    assert failed["status"] == "FAILED_MPNN" and failed["af3_validation"]["status"] == "not_run"
    assert failed["verdict"] == "FAIL" and "Test MPNN error" in failed["failure_reasons"][0]
    assert result["summary"] == {"pass": 3, "review": 0, "fail": 1}
    evidence = [
        json.loads(c["content"]) for c in update["evidence_cards"] if c["source"] == "ProteinMPNN"
    ]
    assert any(
        e["candidate_status"] == "FAILED_MPNN" and e["metric"] == "mpnn_status" for e in evidence
    )


def test_rfd3_candidate_failure_skips_mpnn_only_for_that_candidate(state, fake_models, monkeypatch):
    original = RFdiffusion3Runner.run

    def broken(self, request, **kwargs):
        result = original(self, request, **kwargs)
        Path(result.generated_structures[1]).unlink()
        return result

    monkeypatch.setattr(RFdiffusion3Runner, "run", broken)
    result = screening.screening_node(state)["screening_result"]
    assert result["candidates"][1]["status"] == "FAILED_RFD3"
    assert result["mpnn_tool_calls"] == result["mpnn_completed"] == result["af3_evaluated"] == 3


def test_rfd3_batch_failure_records_target_without_fabricated_candidates(state, monkeypatch):
    monkeypatch.setattr(
        RFdiffusion3Runner,
        "run",
        lambda *a, **k: RFdiffusion3Output(status="failed", error_message="test RFD3 failure"),
    )
    result = screening.screening_node(state)["screening_result"]
    assert result["status"] == "FAILED_EXECUTION" and not result["candidates"]
    assert result["targets"][0]["failure_reasons"] == ["test RFD3 failure"]
    assert result["mpnn_completed"] == result["af3_evaluated"] == 0


def test_mpnn_not_configured_stops_before_af3(state, fake_models, monkeypatch):
    monkeypatch.setattr(ProteinMPNNTool, "run", REAL_MPNN_RUN)
    result = screening.screening_node(state)["screening_result"]
    assert result["status"] == "NOT_CONFIGURED" and result["mpnn_completed"] == 0
    assert not fake_models[1]
    assert all(
        c["status"] == "NOT_CONFIGURED" and c["verdict"] == "REVIEW" for c in result["candidates"]
    )


def test_mpnn_and_planner_budget_caps(state, fake_models):
    state["budget"] = {"max_candidates": 3, "screening": {"max_mpnn_candidates": 2}}
    result = screening.screening_node(state)["screening_result"]
    assert (
        result["generated_rfd3"] == 3 and result["mpnn_completed"] == result["af3_evaluated"] == 2
    )
    assert (
        result["mpnn_requested"] == 2 and result["candidates"][-1]["status"] == "SKIPPED_RESOURCE"
    )
    assert result["requested_rfd3"] == result["requested_designs"] == 4


def test_explicit_multiple_sequences_do_not_multiply_backbone_count(state, fake_models):
    get_config().protein_mpnn.num_sequences_per_backbone = 2
    state["screening_options"] = {
        "max_rfd3_candidates": 2,
        "max_mpnn_candidates": 3,
        "max_af3_candidates": 2,
    }
    result = screening.screening_node(state)["screening_result"]
    assert (
        result["generated_rfd3"] == 2
        and result["mpnn_completed"] == 3
        and result["af3_evaluated"] == 2
    )
    assert (
        len(result["candidates"]) == 3
        and len({c["candidate_id"] for c in result["candidates"]}) == 3
    )


def test_design_mpnn_parameters_and_total_budget_reused(state, fake_models):
    path = Path(state["design_results_path"])
    data = json.loads(path.read_text())
    data["designs"][0]["design_steps"][0]["params"] = {
        "model_version": "v_48_010",
        "sampling_temp": 0.2,
        "omit_AAs": "X",
        "num_seqs": 2,
    }
    path.write_text(json.dumps(data))
    settings, cap = screening._mpnn_plan(data["designs"][0])
    assert (
        settings.model_name == "v_48_010"
        and settings.sampling_temp == 0.2
        and settings.omit_AAs == "X"
    )
    assert cap == 2 and settings.num_sequences_per_backbone == 1
    assert screening.screening_node(state)["screening_result"]["mpnn_completed"] == 2


def test_legacy_af3_cannot_bypass_mpnn(state, fake_models):
    old = screening.screening_node(state)["protein_design_result"]
    for candidate in old["candidates"]:
        candidate.pop("protein_mpnn")
    prior_af3 = len(fake_models[1])
    result = screening.screening_node({**state, "protein_design_result": old})["screening_result"]
    assert len(fake_models[0]) == 1  # Backbone reused.
    assert result["mpnn_tool_calls"] == 4 and len(fake_models[1]) - prior_af3 == 4


def test_current_artifacts_reused_and_mismatched_af3_input_recomputed(
    state, fake_models, monkeypatch
):
    old = screening.screening_node(state)["protein_design_result"]
    path = Path(old["candidates"][1]["af3_validation"]["input_json"])
    data = json.loads(path.read_text())
    data["sequences"][1]["protein"]["sequence"] = "W" * 8
    path.write_text(json.dumps(data))
    monkeypatch.setattr(
        ProteinMPNNTool, "run", Mock(side_effect=AssertionError("MPNN must be reused"))
    )
    previous_calls = len(fake_models[1])
    result = screening.screening_node({**state, "protein_design_result": old})["screening_result"]
    assert result["mpnn_tool_calls"] == 0 and result["mpnn_completed"] == 4
    assert len(fake_models[1]) - previous_calls == 1 and result["af3_evaluated"] == 4


def test_mock_integration_runs_real_structural_logic_without_real_counts(
    state, fake_models, monkeypatch, tmp_path
):
    from ..critic_agent import _apply_deterministic_rules

    original = RFdiffusion3Runner.run

    def mock_generate(self, request, **kwargs):
        result = original(self, request, **kwargs)
        return RFdiffusion3Output(
            status="success",
            generated_structures=result.generated_structures,
            metadata={"is_mock": True},
        )

    monkeypatch.setattr(RFdiffusion3Runner, "run", mock_generate)

    def mock_predict(self, **kwargs):
        return af3_fixture(Path(kwargs["output_dir"]))

    monkeypatch.setattr(AF3Tool, "run", mock_predict)
    update = screening.screening_node(state)
    result = update["screening_result"]
    assert result["execution_mode"] == "mock" and result["is_mock"]
    assert result["generated_rfd3"] == result["mpnn_completed"] == result["af3_evaluated"] == 0
    assert result["summary"] == {"pass": 0, "review": 0, "fail": 0}
    assert result["mock_summary"] == {"pass": 0, "review": 4, "fail": 0}
    assert all(
        c["structural_metrics"]["status"] == "COMPLETED" and c["hotspot_preserved"] is True
        for c in result["candidates"]
    )
    assert _apply_deterministic_rules({**state, **update})["verdict"] == "HOLD"


def test_missing_mpnn_sample_artifact_cannot_pass(state, fake_models, monkeypatch):
    original = ProteinMPNNTool.run

    def missing(self, **kwargs):
        output = original(self, **kwargs)
        Path(output["sequences"][0]["sequence_path"]).unlink()
        return output

    monkeypatch.setattr(ProteinMPNNTool, "run", missing)
    result = screening.screening_node(state)["screening_result"]
    assert result["mpnn_completed"] == result["af3_evaluated"] == 0
    assert all(c["status"] == "FAILED_MPNN" for c in result["candidates"])


def test_invalid_mpnn_config_is_not_configured(state):
    path = Path(state["design_results_path"])
    data = json.loads(path.read_text())
    data["designs"][0]["design_steps"][0]["params"] = {"sampling_temp": -1}
    path.write_text(json.dumps(data))
    result = screening.screening_node(state)["screening_result"]
    assert result["status"] == "NOT_CONFIGURED" and not result["candidates"]


def test_zero_sequence_budget_never_calls_mpnn(state, fake_models, monkeypatch):
    state["screening_options"] = {"max_mpnn_candidates": 0}
    monkeypatch.setattr(ProteinMPNNTool, "run", Mock(side_effect=AssertionError("MPNN called")))
    result = screening.screening_node(state)["screening_result"]
    assert (
        result["status"] == "SKIPPED_RESOURCE"
        and result["mpnn_completed"] == result["af3_evaluated"] == 0
    )


def test_limits_keep_existing_profiles():
    assert ScreeningOptions(resource_profile="A").limits() == (8, 8, 8)
    assert ScreeningOptions(resource_profile="B").limits() == (4, 4, 4)
    assert ScreeningOptions(resource_profile="C").limits() == (0, 0, 0)


def test_global_mpnn_cap_across_targets(state, fake_models):
    path = Path(state["design_results_path"])
    data = json.loads(path.read_text())
    data["designs"].append({**data["designs"][0], "gene_name": "SECOND"})
    path.write_text(json.dumps(data))
    state["advance_targets"].append({"gene_name": "SECOND", "modality": "DE_NOVO_BINDER"})
    state["screening_options"] = {
        "max_targets": 2,
        "max_rfd3_candidates": 5,
        "max_mpnn_candidates": 3,
    }
    result = screening.screening_node(state)["screening_result"]
    assert (
        result["generated_rfd3"] == 5 and result["mpnn_requested"] == result["mpnn_completed"] == 3
    )
    assert result["af3_evaluated"] == 3 and result["status"] == "PARTIAL_EXECUTION"
    assert (
        result["candidates"][-1]["gene_name"] == "SECOND"
        and result["candidates"][-1]["status"] == "SKIPPED_RESOURCE"
    )


def test_existing_rfd3_wrapper_timeout_is_configurable(state, tmp_path, monkeypatch):
    from ..rfdiffusion3_stage import RFdiffusion3Input

    design = json.loads(Path(state["design_results_path"]).read_text())["designs"][0]
    request = RFdiffusion3Input(
        target_structure=design["target_structure"],
        target_chain="T",
        design_type="protein_binder",
        output_dir=str(tmp_path / "rfd3"),
        design_length=8,
        num_designs=1,
    )
    monkeypatch.setattr(RFdiffusion3Runner, "available", lambda self: True)

    def timeout(command, **kwargs):
        assert kwargs["timeout"] == 12
        raise subprocess.TimeoutExpired(command, 12)

    monkeypatch.setattr(subprocess, "run", timeout)
    result = RFdiffusion3Runner(timeout_seconds=12).run(request)
    assert result.status == "failed" and not result.generated_structures


def test_existing_af3_wrapper_timeout_is_configurable(tmp_path, monkeypatch):
    from ..protein_design_route import AlphaFold3BinderAdapter

    script = tmp_path / "run_alphafold.py"
    script.write_text("# Never executed\n")

    def timeout(command, **kwargs):
        assert kwargs["timeout"] == 13
        raise subprocess.TimeoutExpired(command, 13)

    monkeypatch.setattr(subprocess, "run", timeout)
    adapter = AlphaFold3BinderAdapter(
        python_path=sys.executable,
        script_path=str(script),
        model_dir=str(tmp_path),
        database_dir=str(tmp_path),
        timeout_seconds=13,
    )
    result = adapter.run(
        target_sequence="AAAA", binder_sequence="GGGG", output_dir=str(tmp_path / "af3")
    )
    assert result.status == "failed"


def test_stage_counts_preserve_real_upstream_when_later_stage_is_mock(
    state, fake_models, monkeypatch
):
    original = ProteinMPNNTool.run

    def mocked_mpnn(self, **kwargs):
        output = original(self, **kwargs)
        output["is_mock"] = True
        return output

    monkeypatch.setattr(ProteinMPNNTool, "run", mocked_mpnn)
    result = screening.screening_node(state)["screening_result"]
    assert result["execution_mode"] == "mock"
    assert (
        result["generated_rfd3"] == 4
    )  # RFD3's execution record remains separate from mocked MPNN.
    assert result["mpnn_completed"] == result["af3_evaluated"] == result["summary"]["pass"] == 0
