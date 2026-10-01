"""Offline tests: subprocess and sockets are forbidden unless explicitly mocked."""

import json
import socket
import subprocess
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock

import pytest

from ..configuration import DataPaths, OmiCraftConfig, get_config, set_config
from ..protein_design_route import (
    AlphaFold3BinderAdapter,
    extract_protein_sequence,
    protein_route_options,
)
from ..rfdiffusion3_stage import RFdiffusion3Input, RFdiffusion3Runner, _mock_cif


@pytest.fixture(autouse=True)
def offline(monkeypatch, tmp_path):
    monkeypatch.setenv("CUDA_VISIBLE_DEVICES", "0")

    def forbidden(*args, **kwargs):
        raise AssertionError("External process/network is forbidden in offline tests")

    monkeypatch.setattr(subprocess, "run", forbidden)
    monkeypatch.setattr(subprocess, "Popen", forbidden)
    monkeypatch.setattr(socket.socket, "connect", forbidden)
    previous = get_config()
    set_config(OmiCraftConfig(data=DataPaths(base_dir=str(tmp_path / "results"))))
    yield
    set_config(previous)


@pytest.fixture
def request_data(tmp_path):
    target = tmp_path / "target.cif"
    target.write_text(_mock_cif(12, "T"))
    return RFdiffusion3Input(
        target_structure=str(target),
        target_chain="T",
        design_type="protein_binder",
        output_dir=str(tmp_path / "designs"),
        design_length=8,
        hotspot_residues=("T3",),
        motif_residues=("T5",),
        binding_site=("T4",),
        seed=17,
    )


@pytest.fixture
def adapter():
    return AlphaFold3BinderAdapter(
        python_path="/missing/python",
        script_path="/missing/af3.py",
        model_dir="/missing/models",
        database_dir="/missing/db",
    )


def test_default_and_explicit_selection():
    assert protein_route_options({}) == {
        "use_rfd3": False,
        "design_mode": "small_molecule",
        "dry_run": True,
    }
    get_config().rfd3.enabled = True
    assert protein_route_options({})["design_mode"] == "protein_binder"
    assert not protein_route_options({"use_rfd3": False})["use_rfd3"]
    with pytest.raises(ValueError):
        protein_route_options({"use_rfd3": True, "design_mode": "small_molecule"})


def test_rfd3_dry_input_and_mock(request_data):
    result = RFdiffusion3Runner("missing-rfd3").run(request_data, dry_run=True)
    assert result.status == "dry_run"
    spec = result.command_or_config["input_specification"]["protein_binder"]
    assert "9999" not in spec["contig"] and "T12" in spec["contig"]
    assert spec["select_hotspots"] == "T3,T4"
    assert spec["select_fixed_atoms"] is True
    assert "seed=17" in result.command_or_config["command"]
    assert extract_protein_sequence(result.generated_structures[0]) == "A" * 8
    assert RFdiffusion3Runner.parse_output(request_data.output_dir) == ()


@pytest.mark.parametrize(
    "changes",
    [
        {"target_structure": "/does/not/exist"},
        {"target_chain": "Z"},
        {"design_type": "small_molecule"},
        {"design_length": 2},
        {"num_designs": 0},
        {"hotspot_residues": ("T99",)},
    ],
)
def test_invalid_inputs_never_launch(request_data, changes):
    result = RFdiffusion3Runner().run(replace(request_data, **changes), dry_run=True)
    assert result.status == "invalid_input" and not result.generated_structures


def test_pdb_and_reordered_cif_parsing(tmp_path):
    pdb = tmp_path / "binder.pdb"
    pdb.write_text(
        "ATOM      1  CA  GLY B   1       0.000   0.000   0.000  1.00 20.00           C  \n"
        "ATOM      2  CA  SER B   2       3.800   0.000   0.000  1.00 20.00           C  \nEND\n"
    )
    assert extract_protein_sequence(str(pdb), "B") == "GS"
    cif = tmp_path / "binder.cif"
    cif.write_text(
        "data_test\nloop_\n_atom_site.label_seq_id\n_atom_site.label_asym_id\n"
        "_atom_site.label_comp_id\n2 B SER\n3 B GLY\n3 B GLY\n#\n"
    )
    assert extract_protein_sequence(str(cif), "B") == "SG"


def test_af3_missing_dependency(adapter, tmp_path):
    result = adapter.run(target_sequence="AAA", binder_sequence="CCC", output_dir=str(tmp_path))
    assert result.status == "dependency_not_found"


def test_af3_metrics_belong_to_same_model(tmp_path):
    for name, data in [
        ("best", {"ranking_score": 0.9, "iptm": 0.7, "ptm": 0.6, "has_clash": False}),
        ("other", {"ranking_score": 0.5, "iptm": 0.95, "ptm": 0.95, "has_clash": True}),
    ]:
        (tmp_path / f"{name}_model.cif").write_text(_mock_cif())
        (tmp_path / f"{name}_summary_confidences.json").write_text(json.dumps(data))
    metrics = AlphaFold3BinderAdapter.parse_output(tmp_path)
    assert metrics["iptm"] == 0.7 and metrics["has_clash"] is False
    assert metrics["model_path"].endswith("best_model.cif")


@pytest.mark.parametrize("outcome", [0, 1, OSError("cannot execute")])
def test_rfd3_process_failure_or_empty_output(request_data, monkeypatch, outcome):
    monkeypatch.setattr(RFdiffusion3Runner, "available", lambda self: True)
    run = Mock(
        side_effect=outcome if isinstance(outcome, Exception) else None,
        return_value=SimpleNamespace(returncode=outcome),
    )
    monkeypatch.setattr(subprocess, "run", run)
    result = RFdiffusion3Runner().run(request_data)
    assert result.status == "failed" and not result.generated_structures
    assert run.call_count == 1


def test_rfd3_stale_outputs_are_not_success(request_data, monkeypatch):
    monkeypatch.setattr(RFdiffusion3Runner, "available", lambda self: True)
    directory = Path(request_data.output_dir)
    directory.mkdir()
    (directory / "old.cif").write_text(_mock_cif())
    assert RFdiffusion3Runner().run(request_data).status == "invalid_input"


def test_af3_launch_error_is_recorded(tmp_path, monkeypatch):
    executable = tmp_path / "python"
    executable.touch()
    script = tmp_path / "run_alphafold.py"
    script.touch()
    adapter = AlphaFold3BinderAdapter(
        python_path=str(executable),
        script_path=str(script),
        model_dir=str(tmp_path),
        database_dir=str(tmp_path),
    )
    monkeypatch.setattr(subprocess, "run", Mock(side_effect=OSError("cannot execute")))
    result = adapter.run(
        target_sequence="AAA", binder_sequence="CCC", output_dir=str(tmp_path / "output")
    )
    assert result.status == "failed" and "cannot execute" in result.error_message
