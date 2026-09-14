from __future__ import annotations

import pytest

pytest.importorskip("PySide6")

from pathlib import Path

from phyto_reason.agent.conversation_store import ConversationStore, SessionData
from phyto_reason.agent.tool_handlers import handle_run_dam_analysis
from phyto_reason.desktop.client import UploadWorker
from phyto_reason.models.planner_state import PlannerState
from phyto_reason.tools.network.correlation_tool import CorrelationTool
from phyto_reason.workflows.execution_nodes import _infer_sample_groups


def test_desktop_upload_worker_parses_meta_csv_as_697_by_15():
    results = []
    fixture = Path(__file__).parent / "fixtures" / "meta_synthetic.csv"
    worker = UploadWorker(str(fixture), "metabolite")
    worker.result.connect(results.append)
    worker.execute()
    assert results
    payload = results[0]
    assert payload["kind"] == "metabolite"
    assert len(payload["payload"]) == 697
    assert len(next(iter(payload["payload"].values()))) == 15


def test_metadata_groups_take_precedence_and_report_names():
    matrix = {"M": {s: 1.0 for s in ("F-1", "F-2", "S-1", "S-2")}}
    metadata = {
        "samples": [
            {"sample_id": "F-1", "condition": "Flower", "tissue": "flower"},
            {"sample_id": "F-2", "condition": "Flower", "tissue": "flower"},
            {"sample_id": "S-1", "condition": "Stem", "tissue": "stem"},
            {"sample_id": "S-2", "condition": "Stem", "tissue": "stem"},
        ]
    }
    assert _infer_sample_groups(matrix, metadata) == {
        "Flower": ["F-1", "F-2"], "Stem": ["S-1", "S-2"]
    }


def test_session_persists_metadata_and_alignment_report():
    store = ConversationStore(":memory:")
    session = store.get_or_create("metadata-session")
    session.set_data(
        metabolite={"M": {"F-1": 1.0, "S-1": 2.0}},
        sample_metadata={"samples": [{"sample_id": "F-1", "condition": "Flower"}]},
        sample_alignment_report={"n_common": 1, "groups": {"Flower": ["F-1"]}},
    )
    loaded = store.get("metadata-session")
    assert loaded.sample_metadata["samples"][0]["condition"] == "Flower"
    assert loaded.sample_alignment_report["n_common"] == 1


def test_dam_report_has_small_sample_discipline_and_correct_metabolite_name():
    session = SessionData("dam-small")
    session.metabolite_matrix = {
        "Berberine": {"F-1": 1.0, "F-2": 1.1, "S-1": 4.0, "S-2": 4.1}
    }
    text = handle_run_dam_analysis({}, session=session)
    assert "Berberine" in text
    assert "feasibility=exploratory" in text
    assert "不得仅凭 q 值" in text


def test_correlation_default_is_strict_and_explicit_low_threshold_is_exploratory():
    tool = CorrelationTool()
    assert tool.parameters[5].default == 0.8
    matrix_expr = {"TF1": {f"S{i}": float(i) for i in range(1, 6)}}
    matrix_meta = {"M1": {f"S{i}": float(i) for i in range(1, 6)}}
    result = tool.run(
        tf_ids=["TF1"], metabolite_names=["M1"],
        expression_matrix=matrix_expr, metabolite_matrix=matrix_meta,
        corr_threshold=0.5,
    )
    assert result.metadata["corr_threshold"] == 0.5
    assert result.metadata["exploratory"] is True
    assert "探索性相关" in result.warnings[0]
