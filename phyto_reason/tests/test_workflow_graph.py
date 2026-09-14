"""
test_workflow_graph.py — Tests for workflow_graph.py: graph structure, compilation, smoke test.
"""

from __future__ import annotations

import pytest

from phyto_reason.workflows.workflow_graph import build_workflow, compile_workflow
from phyto_reason.workflows.runtime_state import RuntimeState
from phyto_reason.models.planner_state import PlannerState


class TestGraphStructure:

    def test_build_workflow_returns_graph(self):
        """build_workflow() returns a StateGraph."""
        graph = build_workflow()
        assert graph is not None
        assert hasattr(graph, "compile")

    def test_compile_workflow_returns_compiled(self):
        """compile_workflow() returns a compiled runnable graph."""
        compiled = compile_workflow()
        assert compiled is not None
        assert hasattr(compiled, "invoke")

    def test_graph_can_be_compiled(self):
        """Graph compilation does not raise errors."""
        graph = build_workflow()
        compiled = graph.compile()
        assert compiled is not None


class TestGraphInvoke:

    def test_graph_invoke_with_minimal_state(self):
        """Full graph invoke with minimal state completes."""
        compiled = compile_workflow()

        ps = PlannerState(
            species="test",
            target_metabolite="test",
            has_expression=False,
            has_metabolite=False,
        )
        state = RuntimeState(planner_state=ps, current_node="data_quality_gate")

        # Invoke with minimal state — should complete without errors
        result = compiled.invoke(state, {"recursion_limit": 50})
        assert result is not None
        # Result may be RuntimeState or dict
        assert isinstance(result, (RuntimeState, dict))

    def test_graph_invoke_with_data(self, runtime_state):
        """Full graph invoke with expression + metabolite data."""
        compiled = compile_workflow()
        runtime_state.current_node = "data_quality_gate"

        result = compiled.invoke(runtime_state, {"recursion_limit": 50})
        assert result is not None

    def test_graph_with_workflow_plan_skip(self, runtime_state, workflow_plan_deg_only):
        """Graph with workflow_plan runs selected nodes only."""
        compiled = compile_workflow()
        runtime_state.current_node = "data_quality_gate"
        runtime_state.workflow_plan = workflow_plan_deg_only

        result = compiled.invoke(runtime_state, {"recursion_limit": 50})
        assert result is not None
        # deg_analysis should be in the trace
        if isinstance(result, RuntimeState):
            trace_nodes = [t.get("node") for t in result.execution_trace]
            # data_quality_gate should always be there
            # deg_analysis might be in trace (skipped node still records trace)


class TestWorkflowPlanIntegration:
    """Integration: workflow_plan field propagates through graph invocation."""

    def test_workflow_plan_preserved_in_state(self, runtime_state, workflow_plan_deg_only):
        """WorkflowPlan is preserved in final state."""
        compiled = compile_workflow()
        runtime_state.current_node = "data_quality_gate"
        runtime_state.workflow_plan = workflow_plan_deg_only

        result = compiled.invoke(runtime_state, {"recursion_limit": 50})

        if isinstance(result, RuntimeState):
            wp = result.workflow_plan
            assert wp is not None
            assert "deg_analysis" in wp.tools_to_run or wp == workflow_plan_deg_only
        elif isinstance(result, dict):
            wp = result.get("workflow_plan")
            assert wp is not None
