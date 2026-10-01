"""Regressions for graph budgets, shared GPU isolation and reproducible inputs."""

import json
import os
import subprocess
import sys
from dataclasses import replace
from pathlib import Path

import pytest

from ..configuration import RFD3Config, get_config
from ..model_runtime import gpu_environment
from ..protein_design_route import AlphaFold3BinderAdapter, protein_design_node
from ..rfdiffusion3_stage import RFdiffusion3Input, RFdiffusion3Runner, _mock_cif
from ..screening_agent import screening_node
from ..screening_tools import AF3Tool


def test_selected_route_zero_budget_never_launches_models(state, fake_models):
    state.update(
        use_rfd3=True,
        screening_options={
            "max_rfd3_candidates": 0,
            "max_mpnn_candidates": 0,
            "max_af3_candidates": 0,
            "max_protenix_candidates": 0,
            "gpu_device": 3,
        },
    )
    state.update(protein_design_node(state))
    assert state["protein_design_result"]["status"] == "not_run"
    result = screening_node(state)["screening_result"]
    assert fake_models == ([], [])
    assert result["status"] == "SKIPPED_RESOURCE"
    assert result["af3_tool_calls"] == result["mpnn_tool_calls"] == 0


def test_screening_passes_gpu_timeout_and_sampler_options(state, fake_models, monkeypatch):
    config = get_config().rfd3
    config.timeout_seconds = 7
    config.is_non_loopy = True
    config.step_scale = 3.0
    config.gamma_0 = 0.2
    state["screening_options"] = {
        "gpu_device": 3,
        "max_rfd3_candidates": 2,
        "max_af3_candidates": 1,
    }
    previous = RFdiffusion3Runner.run
    observed = []

    def run(self, request, **kwargs):
        observed.append((self.gpu_device, self.timeout_seconds, request))
        return previous(self, request, **kwargs)

    monkeypatch.setattr(RFdiffusion3Runner, "run", run)
    screening_node(state)
    gpu, timeout, request = observed[0]
    assert (gpu, timeout, request.num_designs) == (3, 7, 2)
    assert (request.is_non_loopy, request.step_scale, request.gamma_0) == (True, 3.0, 0.2)
    assert len(fake_models[1]) == 1


@pytest.mark.parametrize("gpu", [4, 5, 6, 7])
def test_runtime_gpu_mask_and_timeout_are_shared(tmp_path, monkeypatch, gpu):
    monkeypatch.setenv("CUDA_VISIBLE_DEVICES", ",".join("GPU-" + str(i) for i in range(4, 8)))
    target = tmp_path / "target.cif"
    target.write_text(_mock_cif(8, "T"))
    request = RFdiffusion3Input(
        str(target), "T", "protein_binder", str(tmp_path / "rfd3"), design_length=8
    )
    script = tmp_path / "af3.py"
    script.write_text("# never executed")
    models = tmp_path / "models"
    models.mkdir()
    database = tmp_path / "db"
    database.mkdir()
    calls = []

    def execute(command, **kwargs):
        assert kwargs["env"]["CUDA_VISIBLE_DEVICES"] == "GPU-" + str(gpu)
        assert kwargs["timeout"] == 7
        calls.append(command)
        return subprocess.CompletedProcess(command, 0)

    monkeypatch.setattr(subprocess, "run", execute)
    monkeypatch.setattr(RFdiffusion3Runner, "available", lambda self: True)
    runner = RFdiffusion3Runner(gpu_device=gpu - 4, timeout_seconds=7)
    runner.run(request)
    adapter = AlphaFold3BinderAdapter(
        python_path=sys.executable,
        script_path=str(script),
        model_dir=str(models),
        database_dir=str(database),
        gpu_device=gpu - 4,
        timeout_seconds=7,
        environment={"LIBCIFPP_DATA_DIR": "/test/cifpp"},
        hmmer_bin_dir="/test/hmmer",
    )
    adapter.run(target_sequence="A" * 8, binder_sequence="G" * 8, output_dir=str(tmp_path / "af3"))
    assert len(calls) == 2
    assert "--gpu_device=0" in calls[1]  # Device zero inside the isolated subprocess.
    assert os.environ["CUDA_VISIBLE_DEVICES"] == "GPU-4,GPU-5,GPU-6,GPU-7"


def test_runtime_rejects_invisible_and_overridden_gpu(monkeypatch):
    monkeypatch.setenv("CUDA_VISIBLE_DEVICES", "GPU-7")
    for gpu in (-1, 1, True):
        with pytest.raises(ValueError):
            gpu_environment(gpu)
    with pytest.raises(ValueError):
        gpu_environment(0, environment={"CUDA_VISIBLE_DEVICES": "0"})
    assert gpu_environment(0)["CUDA_VISIBLE_DEVICES"] == "GPU-7"
    monkeypatch.delenv("CUDA_VISIBLE_DEVICES")
    assert gpu_environment(7)["CUDA_VISIBLE_DEVICES"] == "7"


def test_rfd3_sampler_config_roundtrip_and_validation(tmp_path):
    config = RFD3Config(is_non_loopy=True, step_scale=3, gamma_0=0.2)
    assert RFD3Config(**config.model_dump()) == config
    target = tmp_path / "target.cif"
    target.write_text(_mock_cif(8, "T"))
    request = RFdiffusion3Input(
        str(target),
        "T",
        "protein_binder",
        str(tmp_path / "rfd3"),
        is_non_loopy=True,
        step_scale=config.step_scale,
        gamma_0=config.gamma_0,
    )
    result = RFdiffusion3Runner().run(request, dry_run=True)
    assert result.command_or_config["input_specification"]["protein_binder"]["is_non_loopy"] is True
    assert "inference_sampler.step_scale=3.0" in result.command_or_config["command"]
    assert "inference_sampler.gamma_0=0.2" in result.command_or_config["command"]
    for change in ({"step_scale": 0}, {"gamma_0": float("nan")}, {"is_non_loopy": "false"}):
        assert (
            RFdiffusion3Runner().run(replace(request, **change), dry_run=True).status
            == "invalid_input"
        )


def test_af3_cached_msa_templates_and_seed_cache_compatibility(tmp_path):
    cache = tmp_path / "target_msa.json"
    cache.write_text(
        json.dumps({"sequence": "AAAA", "unpairedMsa": ">query\nAAAA\n", "pairedMsa": ""})
    )
    adapter = AlphaFold3BinderAdapter(
        python_path="/not/executed",
        script_path="/not/executed",
        model_dir="",
        database_dir="",
        target_msa_path=str(cache),
        binder_msa_mode="single_sequence",
        use_templates=False,
    )
    result = adapter.run(
        target_sequence="AAAA",
        binder_sequence="GGGG",
        output_dir=str(tmp_path / "af3"),
        seeds=(11,),
        dry_run=True,
    )
    payload = json.loads(Path(result.input_json).read_text())
    target, binder = [item["protein"] for item in payload["sequences"]]
    assert target["unpairedMsa"] == ">query\nAAAA\n" and target["templates"] == []
    assert binder["unpairedMsa"] == binder["pairedMsa"] == "" and binder["templates"] == []
    output = {"status": "success", "input_json": result.input_json}
    assert AF3Tool.input_matches(output, "AAAA", "GGGG", expected_payload=payload)
    changed = adapter.input_payload(target_sequence="AAAA", binder_sequence="GGGG", seeds=(12,))
    assert not AF3Tool.input_matches(output, "AAAA", "GGGG", expected_payload=changed)
    adapter.binder_msa_mode = "search"
    changed = adapter.input_payload(target_sequence="AAAA", binder_sequence="GGGG", seeds=(11,))
    assert not AF3Tool.input_matches(output, "AAAA", "GGGG", expected_payload=changed)
    mismatch = adapter.run(
        target_sequence="CCCC",
        binder_sequence="GGGG",
        output_dir=str(tmp_path / "mismatch"),
        dry_run=True,
    )
    assert mismatch.status == "failed" and "mismatch" in mismatch.error_message
