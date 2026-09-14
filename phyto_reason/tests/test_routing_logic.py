"""
test_routing_logic.py — Tests for routing_logic.py: expand_analysis_selection, gate functions.
"""

from __future__ import annotations

import pytest

from phyto_reason.workflows.routing_logic import (
    expand_analysis_selection,
    STEP_TO_NODES,
    ALWAYS_RUN_NODES,
    ALL_PIPELINE_NODES,
    STEP_OPTIONS,
    gate_after_quality,
    route_after_classification,
)
from phyto_reason.workflows.runtime_state import RuntimeState
from phyto_reason.models.planner_state import PlannerState
from langgraph.graph import END


class TestExpandAnalysisSelection:

    def test_expand_empty_returns_all(self):
        """Empty list returns all pipeline nodes (backward compatible)."""
        result = expand_analysis_selection([])
        assert sorted(result) == sorted(ALL_PIPELINE_NODES)

    def test_expand_all_returns_all(self):
        """'all' keyword returns all pipeline nodes."""
        result = expand_analysis_selection(["all"])
        assert sorted(result) == sorted(ALL_PIPELINE_NODES)

    def test_expand_deg_only(self):
        """DEG expands to deg_analysis + always-run nodes."""
        result = expand_analysis_selection(["DEG"])
        assert "deg_analysis" in result
        for node in ALWAYS_RUN_NODES:
            assert node in result
        # Should NOT include nodes not in DEG
        assert "multiomics_integration" not in result

    def test_expand_dam_only(self):
        """DAM expands to dam_analysis + always-run nodes."""
        result = expand_analysis_selection(["DAM"])
        assert "dam_analysis" in result
        for node in ALWAYS_RUN_NODES:
            assert node in result
        assert "deg_analysis" not in result

    def test_expand_multiomics(self):
        """multiomics expands to deg + dam + multiomics + always-run."""
        result = expand_analysis_selection(["multiomics"])
        assert "deg_analysis" in result
        assert "dam_analysis" in result
        assert "multiomics_integration" in result
        for node in ALWAYS_RUN_NODES:
            assert node in result

    def test_expand_hypothesis(self):
        """hypothesis expands to regulation, contradiction, falsification, competition, synthesis + always-run."""
        result = expand_analysis_selection(["hypothesis"])
        assert "regulation_evidence" in result
        assert "contradiction_check" in result
        assert "falsification" in result
        assert "hypothesis_competition" in result
        assert "hypothesis_synthesis" in result
        for node in ALWAYS_RUN_NODES:
            assert node in result

    def test_expand_multiple_steps_union(self):
        """Multiple steps union their node sets."""
        result = expand_analysis_selection(["DEG", "tf_narrowing"])
        assert "deg_analysis" in result
        assert "tf_narrowing" in result
        for node in ALWAYS_RUN_NODES:
            assert node in result
        # No duplicates
        assert len(result) == len(set(result))

    def test_expand_unknown_step_ignored(self):
        """Unknown step IDs are silently ignored."""
        result = expand_analysis_selection(["NOT_A_STEP", "DEG"])
        assert "deg_analysis" in result

    def test_expand_always_run_nodes_are_first(self):
        """Always-run nodes appear first in the ordered result."""
        result = expand_analysis_selection(["hypothesis"])
        # Check that all always-run nodes come before any selected nodes
        selected_nodes = {"regulation_evidence", "contradiction_check",
                          "falsification", "hypothesis_competition", "hypothesis_synthesis"}
        found_selected = False
        for node in result:
            if node in selected_nodes:
                found_selected = True
            elif found_selected:
                # once we've seen a selected node, remaining should NOT be always-run
                assert node not in ALWAYS_RUN_NODES


class TestSTEP_OPTIONS:
    """Verify STEP_OPTIONS structure for frontend consumption."""

    def test_step_options_has_eleven_categories(self):
        assert len(STEP_OPTIONS) == 11

    def test_step_options_each_has_id_label_description(self):
        for opt in STEP_OPTIONS:
            assert "id" in opt
            assert "label" in opt
            assert "description" in opt

    def test_step_options_ids_match_step_to_nodes(self):
        for opt in STEP_OPTIONS:
            assert opt["id"] in STEP_TO_NODES


class TestGateAfterQuality:

    def test_passed_returns_metabolite_validation(self, runtime_state):
        """When data quality passed, route to metabolite_validation."""
        runtime_state.data_quality_passed = True
        assert gate_after_quality(runtime_state) == "metabolite_validation"

    def test_failed_returns_end(self, runtime_state):
        """When data quality failed, route to END."""
        runtime_state.data_quality_passed = False
        assert gate_after_quality(runtime_state) == END


class TestRouteAfterClassification:

    def test_transport_returns_transport_evidence(self, runtime_state):
        """transport_redistribution mechanism routes to transport_evidence."""
        runtime_state.mechanism_type = "transport_redistribution"
        assert route_after_classification(runtime_state) == "transport_evidence"

    def test_stress_returns_stress_evidence(self, runtime_state):
        """stress_induced_redistribution mechanism routes to stress_evidence."""
        runtime_state.mechanism_type = "stress_induced_redistribution"
        assert route_after_classification(runtime_state) == "stress_evidence"

    def test_transcriptional_returns_tf_narrowing(self, runtime_state):
        """transcriptional_regulation mechanism routes to tf_narrowing."""
        runtime_state.mechanism_type = "transcriptional_regulation"
        assert route_after_classification(runtime_state) == "tf_narrowing"

    def test_unknown_returns_tf_narrowing(self, runtime_state):
        """Unknown mechanism defaults to tf_narrowing."""
        runtime_state.mechanism_type = "unknown"
        assert route_after_classification(runtime_state) == "tf_narrowing"

    def test_empty_returns_tf_narrowing(self, runtime_state):
        """Empty mechanism type defaults to tf_narrowing."""
        runtime_state.mechanism_type = ""
        assert route_after_classification(runtime_state) == "tf_narrowing"
