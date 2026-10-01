import ast
import json
from pathlib import Path

import pytest
from pydantic import ValidationError

from .. import __main__, screening_agent
from ..configuration import OmiCraftConfig
from ..input_paths import local_paths
from ..orchestrator import run_screening
from ..routes import validate_handoff_routes, validate_route


@pytest.mark.parametrize("route", ["antibody_adc", "rna_therapeutic"])
def test_unsupported_routes_are_rejected(route):
    with pytest.raises(ValueError, match="Unknown route"):
        validate_route(route)


@pytest.mark.parametrize("modality", ["RNA_THERAPEUTIC"])
def test_unsupported_handoff_is_rejected(modality):
    with pytest.raises(ValueError, match="Unsupported modality"):
        validate_handoff_routes([{"gene_name": "TARGET", "modality": modality}])


@pytest.mark.parametrize("field", ["therapeutic_modalities", "llm", "cohort", "analysis"])
def test_configuration_rejects_unrelated_sections(field):
    with pytest.raises(ValidationError):
        OmiCraftConfig(**{field: {}})


def test_relative_model_paths_preserve_virtualenv_symlink(tmp_path):
    executable = tmp_path / "python-real"
    executable.touch()
    link = tmp_path / "env/bin/python"
    link.parent.mkdir(parents=True)
    link.symlink_to(executable)
    result = local_paths(
        {
            "python_path": "env/bin/python",
            "target_structure": "target.cif",
            "required_files": ["checkpoint.pt"],
            "protenix_msa_paths": {"A": "target.a3m"},
        },
        tmp_path,
    )
    assert result["python_path"] == str(link)
    assert result["target_structure"] == str(tmp_path / "target.cif")
    assert result["required_files"] == [str(tmp_path / "checkpoint.pt")]
    assert result["protenix_msa_paths"]["A"] == str(tmp_path / "target.a3m")


def test_handoff_paths_are_relative_to_handoff_file(tmp_path):
    path = tmp_path / "inputs/design.json"
    path.parent.mkdir()
    path.write_text(
        json.dumps(
            {
                "designs": [
                    {
                        "gene_name": "TARGET",
                        "modality": "DE_NOVO_BINDER",
                        "target_structure": "target.cif",
                    }
                ]
            }
        )
    )
    designs, errors = screening_agent._design_handoff({"design_results_path": str(path)})
    assert not errors
    assert designs[0]["target_structure"] == str(path.parent / "target.cif")


def test_agent_screening_requires_explicit_selection(state, fake_models):
    result = run_screening(state)
    assert result["execution_status"] == "BLOCKED"
    assert not result["modality_results"]


def test_direct_cli_preserves_audited_dossier(state, fake_models):
    from ..orchestrator import _run_screening_cli
    result = _run_screening_cli(state)
    saved = json.loads(Path(result["dossier_path"]).read_text())
    assert saved["critic"]["verdict"] == "APPROVE"
    assert saved["screening_result"]["summary"]["pass"] == 4


def test_mock_dossier_stays_on_hold(state, fake_models):
    from ..orchestrator import _run_screening_cli
    result = _run_screening_cli({**state, "dry_run": True})
    assert result["critic_verdict"] == "HOLD"


def test_safety_veto_is_preserved(state, fake_models):
    from ..orchestrator import _run_screening_cli
    state["advance_targets"][0]["safety_verdict"] = "REJECT"
    assert _run_screening_cli(state)["critic_verdict"] == "REJECT"


def test_binder_cli_runs_without_model_processes(tmp_path):
    assert __main__.main(["binder", "--demo", "--output-dir", str(tmp_path / "run")]) == 0
    saved = json.loads((tmp_path / "run/agent_results/final_dossier.json").read_text())
    assert saved["critic"]["verdict"] == "HOLD"
    assert saved["screening_result"]["af3_attempted"] == 0


def test_source_imports_do_not_reference_other_projects():
    root = Path(screening_agent.__file__).parent
    forbidden = ("challenge.", "agent_rev2", "agent_rev4", "agent_rev5", "agent_rev6", "omicraft.agent", "therapeutic_routes", "deg_analysis", "preprocessing")
    for path in root.rglob("*.py"):
        if "tests" in path.relative_to(root).parts:
            continue
        for node in ast.walk(ast.parse(path.read_text())):
            if isinstance(node, ast.Import):
                names = [alias.name for alias in node.names]
            elif isinstance(node, ast.ImportFrom):
                names = [node.module or ""] if node.level == 0 else []
            else:
                continue
            assert not any(name.startswith(forbidden) for name in names), path
