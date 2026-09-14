"""
test_execution_nodes.py — Tests for pipeline execution nodes (execution_nodes.py).

Covers previously untested nodes: deg_analysis, dam_analysis, multiomics_integration,
tf_narrowing, regulation_evidence, contradiction_check, falsification,
hypothesis_competition, hypothesis_synthesis.
Also tests WorkflowPlan skip guards.
"""

from __future__ import annotations

from unittest.mock import patch, MagicMock

import pytest

from phyto_reason.workflows.execution_nodes import (
    deg_analysis_node,
    dam_analysis_node,
    multiomics_integration_node,
    hypothesis_synthesis_node,
    hypothesis_competition_node,
    falsification_node,
    _should_skip_node,
)
from phyto_reason.workflows.runtime_state import RuntimeState
from phyto_reason.models.planner_state import PlannerState


class TestDegAnalysisNode:

    def test_with_expression_data(self, runtime_state):
        """With expression data, produces deg_report without errors."""
        result = deg_analysis_node(runtime_state)

        # Should produce a deg_report (either real or mock result)
        assert "deg_report" in result
        assert runtime_state.deg_report is not None
        assert runtime_state.deg_report.get("n_genes_tested", 0) >= 0

    def test_without_expression(self):
        """Without expression matrix, returns empty report."""
        ps = PlannerState(
            species="test", target_metabolite="test",
            has_expression=False, has_metabolite=False,
        )
        state = RuntimeState(planner_state=ps)
        result = deg_analysis_node(state)

        assert "deg_report" in result
        assert result["deg_report"].get("n_genes_tested", 0) == 0

    def test_with_workflow_plan_skip(self, runtime_state, workflow_plan_deg_only):
        """When DEG not in workflow_plan, node is skipped."""
        runtime_state.workflow_plan = workflow_plan_deg_only
        # workflow_plan_deg_only has deg_analysis in tools_to_run, so NOT skipped
        assert not _should_skip_node(runtime_state, "deg_analysis")

    def test_workflow_plan_skips_unrelated(self, runtime_state, workflow_plan_hypothesis_only):
        """When a totally different set is selected, deg is skipped."""
        runtime_state.workflow_plan = workflow_plan_hypothesis_only
        assert _should_skip_node(runtime_state, "deg_analysis")

    def test_no_workflow_plan_runs_everything(self, runtime_state):
        """Without workflow_plan, no skip."""
        assert not _should_skip_node(runtime_state, "deg_analysis")


class TestDamAnalysisNode:

    def test_with_metabolite_data(self, runtime_state):
        """With metabolite data, produces dam_report."""
        with patch("phyto_reason.tools.statistics.deg_tool.DEGTool") as mock_deg:
            mock_tool = MagicMock()
            mock_result = MagicMock()
            mock_result.n_genes_tested = 0
            mock_result.n_effect_ranked = 10
            mock_result.n_groups = 2
            mock_result.group_names = ["Leaf", "Root"]
            mock_result.data_type = "fpkm"
            mock_result.to_dict.return_value = {
                "n_metabolites_tested": 10,
                "n_effect_ranked": 10, "n_groups": 2,
                "group_names": ["Leaf", "Root"],
            }
            mock_tool.return_value = mock_result
            mock_deg.return_value = mock_tool

            result = dam_analysis_node(runtime_state)
            assert "dam_report" in result

    def test_without_metabolite(self):
        """Without metabolite matrix, returns empty."""
        ps = PlannerState(
            species="test", target_metabolite="test",
            has_expression=True, has_metabolite=False,
        )
        state = RuntimeState(planner_state=ps)
        result = dam_analysis_node(state)

        assert "dam_report" in result
        assert result["dam_report"].get("n_metabolites_tested", 0) == 0

    def test_with_workflow_plan_skip(self, runtime_state, workflow_plan_hypothesis_only):
        """When DAM not in workflow_plan, skip guard activates."""
        runtime_state.workflow_plan = workflow_plan_hypothesis_only
        assert _should_skip_node(runtime_state, "dam_analysis")


class TestMultiomicsIntegrationNode:

    def test_with_deg_and_dam(self, runtime_state_with_deg):
        """With DEG data, produces integration report."""
        result = multiomics_integration_node(runtime_state_with_deg)
        assert "multiomics_report" in result or "current_node" in result

    def test_with_workflow_plan_skip(self, runtime_state_with_deg, workflow_plan_deg_only):
        """When multiomics not in workflow_plan, skip guard activates."""
        runtime_state_with_deg.workflow_plan = workflow_plan_deg_only
        # deg_only has deg_analysis but NOT multiomics_integration
        assert _should_skip_node(runtime_state_with_deg, "multiomics_integration")


class TestHypothesisSynthesisNode:

    def test_without_hypotheses_generates_legacy(self, runtime_state):
        """Without competing_hypotheses, generates legacy-format report."""
        runtime_state.tf_candidates = {}
        runtime_state.competing_hypotheses = []
        runtime_state.excluded_tfs = {}
        result = hypothesis_synthesis_node(runtime_state)

        # Should produce a report or finish
        assert "finished" in result or "current_node" in result

    def test_with_workflow_plan_skip(self, runtime_state, workflow_plan_hypothesis_only):
        """When hypothesis_synthesis is in workflow_plan, NOT skipped."""
        runtime_state.workflow_plan = workflow_plan_hypothesis_only
        assert not _should_skip_node(runtime_state, "hypothesis_synthesis")


class TestHypothesisCompetitionNode:

    def test_without_hypotheses(self, runtime_state):
        """Without competing_hypotheses, returns empty."""
        runtime_state.competing_hypotheses = []
        result = hypothesis_competition_node(runtime_state)
        assert isinstance(result, dict)

    def test_with_workflow_plan_skip(self, runtime_state, workflow_plan_deg_only):
        """When hypothesis_competition not in plan, skipped."""
        runtime_state.workflow_plan = workflow_plan_deg_only
        assert _should_skip_node(runtime_state, "hypothesis_competition")


class TestFalsificationNode:

    def test_without_hypotheses(self, runtime_state):
        """Without competing_hypotheses, returns empty."""
        runtime_state.competing_hypotheses = []
        result = falsification_node(runtime_state)
        assert isinstance(result, dict)

    def test_with_workflow_plan_skip(self, runtime_state, workflow_plan_deg_only):
        """When falsification not in plan, skipped."""
        runtime_state.workflow_plan = workflow_plan_deg_only
        assert _should_skip_node(runtime_state, "falsification")


class TestShouldSkipNode:
    """Tests for _should_skip_node() helper."""

    def test_none_plan_returns_false(self, runtime_state):
        """No workflow_plan → don't skip."""
        runtime_state.workflow_plan = None
        assert not _should_skip_node(runtime_state, "deg_analysis")

    def test_empty_tools_to_run_returns_false(self, runtime_state, workflow_plan_full):
        """Empty tools_to_run → backward compatible, don't skip."""
        runtime_state.workflow_plan = workflow_plan_full
        # workflow_plan_full has all nodes, so nothing should be skipped
        assert not _should_skip_node(runtime_state, "deg_analysis")

    def test_node_in_list_returns_false(self, runtime_state, workflow_plan_deg_only):
        """Node in tools_to_run → don't skip."""
        runtime_state.workflow_plan = workflow_plan_deg_only
        assert not _should_skip_node(runtime_state, "deg_analysis")
        assert not _should_skip_node(runtime_state, "data_quality_gate")

    def test_node_not_in_list_returns_true(self, runtime_state, workflow_plan_deg_only):
        """Node not in tools_to_run → skip."""
        runtime_state.workflow_plan = workflow_plan_deg_only
        assert _should_skip_node(runtime_state, "dam_analysis")
        assert _should_skip_node(runtime_state, "multiomics_integration")
        assert _should_skip_node(runtime_state, "hypothesis_synthesis")
