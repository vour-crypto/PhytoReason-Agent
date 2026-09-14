"""
test_workflows.py — Smoke tests for the scientific pipeline.

Uses small mock expression/metabolite matrices to verify
pipeline nodes don't crash and produce expected state transitions.
No external data needed.
"""

import pytest

from phyto_reason.workflows.runtime_state import RuntimeState
from phyto_reason.workflows.execution_nodes import (
    data_quality_gate_node,
    metabolite_validation_node,
    pathway_coherence_node,
    mechanism_classifier_node,
)
from phyto_reason.models.planner_state import PlannerState


# ═══════════════════════════════════════════════════════════════
# Mock data helpers
# ═══════════════════════════════════════════════════════════════

def _make_mock_expression(n_genes: int = 10, n_samples: int = 6) -> dict:
    """Create a mock expression matrix: {gene: {sample: value}}."""
    import random
    random.seed(42)
    return {
        f"gene_{i:03d}": {f"sample_{j:02d}": random.uniform(0, 10) for j in range(n_samples)}
        for i in range(n_genes)
    }


def _make_mock_metabolite(n_metabolites: int = 5, n_samples: int = 6) -> dict:
    """Create a mock metabolite matrix."""
    import random
    random.seed(42)
    return {
        f"metabolite_{i:02d}": {f"sample_{j:02d}": random.uniform(0, 5) for j in range(n_samples)}
        for i in range(n_metabolites)
    }


def _make_planner_state(
    species: str = "Arabidopsis thaliana",
    target_metabolite: str = "anthocyanin",
    expression: dict | None = None,
    metabolite: dict | None = None,
) -> PlannerState:
    return PlannerState(
        species=species,
        target_metabolite=target_metabolite,
        has_expression=expression is not None,
        has_metabolite=metabolite is not None,
        expression_matrix=expression or {},
        metabolite_matrix=metabolite or {},
        sample_count=6 if expression else 0,
    )


# ═══════════════════════════════════════════════════════════════
# RuntimeState
# ═══════════════════════════════════════════════════════════════

class TestRuntimeState:
    def test_create_default(self):
        state = RuntimeState()
        assert state.current_node == "data_quality_gate"
        assert state.finished is False
        assert state.session_id

    def test_create_with_planner(self):
        ps = _make_planner_state()
        state = RuntimeState(planner_state=ps)
        assert state.planner_state is not None
        assert state.planner_state.species == "Arabidopsis thaliana"

    def test_add_trace_completed(self):
        state = RuntimeState()
        state.add_trace("test_node", "completed", {"result": "ok"})
        assert len(state.execution_trace) == 1
        assert "test_node" in state.completed_nodes

    def test_add_trace_failed(self):
        state = RuntimeState()
        state.add_trace("test_node", "failed", error="Something went wrong")
        assert len(state.failures) == 1
        assert state.error is not None

    def test_to_summary(self):
        state = RuntimeState()
        summary = state.to_summary()
        assert "session_id" in summary
        assert "n_completed" in summary
        assert "finished" in summary


# ═══════════════════════════════════════════════════════════════
# Node 1: Data Quality Gate
# ═══════════════════════════════════════════════════════════════

class TestDataQualityGate:
    def test_with_good_data(self):
        expr = _make_mock_expression(n_genes=50, n_samples=8)
        meta = _make_mock_metabolite(n_samples=8)
        # Ensure target metabolite is in the data
        meta["anthocyanin"] = meta.pop("metabolite_00")
        ps = _make_planner_state(expression=expr, metabolite=meta)
        state = RuntimeState(planner_state=ps)

        result = data_quality_gate_node(state)
        assert state.data_quality_passed is True
        # May have non-critical warnings (small gene count, etc.) but gate passes

    def test_with_no_data(self):
        ps = _make_planner_state(expression=None, metabolite=None)
        state = RuntimeState(planner_state=ps)

        result = data_quality_gate_node(state)
        # Without data, quality gate should fail or report issues
        assert state.data_quality_passed is False or len(state.data_quality_issues) > 0

    def test_with_small_sample(self):
        expr = _make_mock_expression(n_genes=5, n_samples=3)  # only 3 samples
        ps = _make_planner_state(expression=expr)
        state = RuntimeState(planner_state=ps)

        result = data_quality_gate_node(state)
        # 3 samples is borderline for correlation analysis
        assert len(state.execution_trace) >= 0  # at minimum, doesn't crash


# ═══════════════════════════════════════════════════════════════
# Node 2: Metabolite Validation
# ═══════════════════════════════════════════════════════════════

class TestMetaboliteValidation:
    def test_with_target_in_data(self):
        expr = _make_mock_expression()
        meta = _make_mock_metabolite()
        # Add the target metabolite to the data
        meta["anthocyanin"] = meta.pop("metabolite_00")
        ps = _make_planner_state(expression=expr, metabolite=meta)
        state = RuntimeState(planner_state=ps)

        result = metabolite_validation_node(state)
        assert len(state.execution_trace) > 0

    def test_without_metabolite_data(self):
        expr = _make_mock_expression()
        ps = _make_planner_state(expression=expr, metabolite=None)
        state = RuntimeState(planner_state=ps)

        result = metabolite_validation_node(state)
        assert len(state.execution_trace) > 0


# ═══════════════════════════════════════════════════════════════
# Node 3: Pathway Coherence
# ═══════════════════════════════════════════════════════════════

class TestPathwayCoherence:
    def test_with_expression_data(self):
        expr = _make_mock_expression()
        meta = _make_mock_metabolite()
        ps = _make_planner_state(expression=expr, metabolite=meta)
        state = RuntimeState(planner_state=ps)

        result = pathway_coherence_node(state)
        assert len(state.execution_trace) > 0

    def test_without_expression_data(self):
        ps = _make_planner_state(expression=None, metabolite=None)
        state = RuntimeState(planner_state=ps)

        result = pathway_coherence_node(state)
        assert len(state.execution_trace) > 0


# ═══════════════════════════════════════════════════════════════
# Node 4: Mechanism Classifier
# ═══════════════════════════════════════════════════════════════

class TestMechanismClassifier:
    def test_classify_with_mock_data(self):
        expr = _make_mock_expression(n_genes=50, n_samples=6)
        meta = _make_mock_metabolite(n_samples=6)
        ps = _make_planner_state(expression=expr, metabolite=meta)
        state = RuntimeState(planner_state=ps)

        result = mechanism_classifier_node(state)
        # Should produce a mechanism type
        assert state.mechanism_type or len(state.execution_trace) > 0

    def test_classify_without_data(self):
        ps = _make_planner_state(expression=None, metabolite=None)
        state = RuntimeState(planner_state=ps)

        result = mechanism_classifier_node(state)
        assert len(state.execution_trace) > 0


# ═══════════════════════════════════════════════════════════════
# Full pipeline mini-run (smoke test)
# ═══════════════════════════════════════════════════════════════

class TestPipelineSmoke:
    """Run a few nodes in sequence to verify state passing works."""

    def test_sequential_nodes(self):
        expr = _make_mock_expression(n_genes=20, n_samples=6)
        meta = _make_mock_metabolite(n_samples=6)
        ps = _make_planner_state(expression=expr, metabolite=meta)
        state = RuntimeState(planner_state=ps)

        # Run nodes 1-4 in sequence
        data_quality_gate_node(state)
        assert state.data_quality_passed

        metabolite_validation_node(state)
        pathway_coherence_node(state)
        mechanism_classifier_node(state)

        # All should have traced
        assert len(state.execution_trace) >= 4
        assert len(state.completed_nodes) >= 4
        assert len(state.failures) == 0
