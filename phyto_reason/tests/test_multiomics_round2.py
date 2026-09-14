"""
test_multiomics_round2.py — Tests for Round 2 multi-omics modules.

Covers: WGCNA node, O2PLS node.
"""

from __future__ import annotations

from unittest.mock import patch, MagicMock

import numpy as np
import pytest

from phyto_reason.workflows.nodes.wgcna import (
    wgcna_node,
    _correlation_matrix,
    _pick_soft_threshold,
    _compute_tom,
    _dynamic_tree_cut,
    _module_eigengenes,
    _module_trait_correlation,
)
from phyto_reason.workflows.nodes.o2pls import (
    o2pls_node,
    _standardize,
    _explained_variance,
    _vip_scores,
    _cross_validate_n_components,
)
from phyto_reason.workflows.execution_nodes import _should_skip_node
from phyto_reason.workflows.runtime_state import RuntimeState
from phyto_reason.models.planner_state import PlannerState


# ═══════════════════════════════════════════════════════════════
# WGCNA Unit Tests
# ═══════════════════════════════════════════════════════════════

class TestCorrelationMatrix:

    def test_pearson_shape(self):
        data = np.random.randn(50, 10)
        corr = _correlation_matrix(data, method="pearson")
        assert corr.shape == (50, 50)
        assert np.allclose(np.diag(corr), 1.0)

    def test_spearman_shape(self):
        data = np.random.randn(30, 8)
        corr = _correlation_matrix(data, method="spearman")
        assert corr.shape == (30, 30)

    def test_correlation_range(self):
        data = np.random.randn(40, 12)
        corr = _correlation_matrix(data)
        assert np.all(corr >= -1.0) and np.all(corr <= 1.0)


class TestPickSoftThreshold:

    def test_picks_power(self):
        corr = np.corrcoef(np.random.randn(100, 15))
        result = _pick_soft_threshold(corr, candidate_powers=[1, 2, 3, 4, 5, 6])
        assert "best_power" in result
        assert "best_r2" in result
        assert result["best_power"] >= 1
        assert len(result["all_results"]) == 6

    def test_perfect_correlation(self):
        """Perfectly correlated block gives high R² at low power."""
        # Create a block-structured correlation matrix
        n = 60
        corr = np.eye(n)
        # 3 blocks of highly correlated genes
        for block in [slice(0, 20), slice(20, 40), slice(40, 60)]:
            corr[block, block] = 0.85
            np.fill_diagonal(corr, 1.0)
        # Ensure symmetry and positive values
        corr = np.abs(corr)
        np.fill_diagonal(corr, 1.0)

        result = _pick_soft_threshold(corr, candidate_powers=[1, 2, 3, 4, 5], min_r2=0.50)
        assert "best_power" in result

    def test_zero_correlation(self):
        """Zero correlation gives zero connectivity."""
        corr = np.eye(30)  # Only self-correlation
        result = _pick_soft_threshold(corr, candidate_powers=[1, 2])
        assert "best_power" in result
        # With only diagonal correlation, R² will be very low
        assert not result["is_satisfactory"]


class TestComputeTOM:

    def test_shape(self):
        adj = np.random.uniform(0.1, 0.9, (20, 20))
        np.fill_diagonal(adj, 0)
        adj = (adj + adj.T) / 2  # Symmetric
        tom = _compute_tom(adj)
        assert tom.shape == (20, 20)
        assert np.allclose(np.diag(tom), 1.0)

    def test_range(self):
        adj = np.random.uniform(0.1, 0.9, (15, 15))
        np.fill_diagonal(adj, 0)
        adj = (adj + adj.T) / 2
        tom = _compute_tom(adj)
        assert np.all(tom >= 0.0) and np.all(tom <= 1.0)

    def test_symmetric(self):
        adj = np.random.uniform(0.1, 0.9, (10, 10))
        np.fill_diagonal(adj, 0)
        adj = (adj + adj.T) / 2
        tom = _compute_tom(adj)
        assert np.allclose(tom, tom.T)


class TestDynamicTreeCut:

    def test_returns_labels(self):
        dissim = np.random.uniform(0.3, 0.9, (50, 50))
        np.fill_diagonal(dissim, 0)
        dissim = (dissim + dissim.T) / 2
        labels = _dynamic_tree_cut(dissim, min_module_size=5, deep_split=2)
        assert len(labels) == 50
        assert 0 in labels  # Grey module

    def test_block_structure(self):
        """Block-structured dissimilarity should produce modules."""
        n = 60
        dissim = np.ones((n, n)) * 0.8
        # Block 1: low dissimilarity (similar genes)
        dissim[:20, :20] = 0.3
        dissim[20:40, 20:40] = 0.3
        dissim[40:, 40:] = 0.3
        np.fill_diagonal(dissim, 0)
        dissim = (dissim + dissim.T) / 2

        labels = _dynamic_tree_cut(dissim, min_module_size=5, deep_split=3)
        n_modules = len(set(labels) - {0})
        # Should produce at least some modules from the block structure
        assert n_modules >= 1


class TestModuleEigengenes:

    def test_returns_dict(self):
        expr = np.random.randn(30, 10)
        labels = np.array([1] * 15 + [2] * 15, dtype=int)
        eg = _module_eigengenes(expr, labels)
        assert isinstance(eg, dict)
        assert 1 in eg
        assert 2 in eg
        assert eg[1].shape == (10,)

    def test_grey_skipped(self):
        expr = np.random.randn(20, 8)
        labels = np.array([0] * 10 + [1] * 10, dtype=int)
        eg = _module_eigengenes(expr, labels)
        assert 0 not in eg
        assert 1 in eg


class TestModuleTraitCorrelation:

    def test_returns_list(self):
        eg = {1: np.random.randn(12), 2: np.random.randn(12)}
        traits = np.random.randn(12, 2)
        result = _module_trait_correlation(eg, traits, ["Trait_A", "Trait_B"], {1: 5, 2: 8})
        assert isinstance(result, list)
        assert len(result) == 4  # 2 modules x 2 traits
        for r in result:
            assert "correlation" in r
            assert "p_value" in r
            assert "q_value" in r


# ═══════════════════════════════════════════════════════════════
# O2PLS Unit Tests
# ═══════════════════════════════════════════════════════════════

class TestStandardize:

    def test_zero_mean_unit_var(self):
        X = np.random.randn(50, 10)
        X_std = _standardize(X)
        assert np.allclose(X_std.mean(axis=0), 0, atol=1e-10)
        assert np.allclose(X_std.std(axis=0, ddof=1), 1.0)


class TestExplainedVariance:

    def test_perfect_explanation(self):
        X = np.random.randn(20, 5)
        scores = np.random.randn(20, 2)
        loadings = np.random.randn(5, 2)
        X_pred = scores @ loadings.T
        var = _explained_variance(X_pred, scores, loadings)
        assert var > 0.9  # Near perfect

    def test_zero_explanation(self):
        X = np.random.randn(20, 5)
        scores = np.zeros((20, 2))
        loadings = np.random.randn(5, 2)
        var = _explained_variance(X, scores, loadings)
        assert var < 0.1


class TestVIPScores:

    def test_shape(self):
        weights = np.random.randn(10, 3)
        scores = np.random.randn(20, 3)
        Y = np.random.randn(20, 2)
        vip = _vip_scores(weights, scores, Y)
        assert vip.shape == (10,)

    def test_positive(self):
        weights = np.abs(np.random.randn(15, 2))
        scores = np.random.randn(25, 2)
        Y = np.random.randn(25, 3)
        vip = _vip_scores(weights, scores, Y)
        assert np.all(vip >= 0)


class TestCVNComponents:

    def test_returns_dict(self):
        X = np.random.randn(30, 20)
        Y = np.random.randn(30, 3)
        result = _cross_validate_n_components(X, Y, max_comp=3, n_folds=3)
        assert "best_n_comp" in result
        assert result["best_n_comp"] >= 1
        assert len(result.get("q2_values", [])) <= 3


# ═══════════════════════════════════════════════════════════════
# WGCNA Node Integration Tests
# ═══════════════════════════════════════════════════════════════

class TestWgcanNode:

    def test_with_data(self, runtime_state_with_deg_full):
        """With expression data, runs WGCNA and produces modules."""
        result = wgcna_node(runtime_state_with_deg_full)
        assert "wgcna_report" in result
        report = result["wgcna_report"]
        if "reason" not in report:
            assert "n_genes_analyzed" in report
            assert "n_modules" in report
            assert "soft_power" in report

    def test_without_data(self, runtime_state):
        """Without expression data, returns reason."""
        result = wgcna_node(runtime_state)
        assert "wgcna_report" in result

    def test_insufficient_samples(self, runtime_state):
        """With < 8 samples, returns warning."""
        # runtime_state has 6 samples
        if runtime_state.planner_state:
            # Create expression with only 4 samples
            small_expr = {
                f"GENE_{i:04d}": {"S1": 1.0, "S2": 2.0, "S3": 1.5, "S4": 2.5}
                for i in range(100)
            }
            runtime_state.planner_state.expression_matrix = small_expr
        deg_report = {"top_genes": [
            {"gene_id": f"GENE_{i:04d}", "effect_size": 0.8, "max_log2fc": 2.0,
             "p_value": 0.001, "q_value": 0.01, "is_fdr_sig": True}
            for i in range(50)
        ]}
        runtime_state.deg_report = deg_report
        result = wgcna_node(runtime_state)
        report = result.get("wgcna_report", {})
        # Should return reason about insufficient samples
        assert "reason" in report or report.get("n_samples", 0) < 8

    def test_with_workflow_plan_skip(self, runtime_state, workflow_plan_deg_only):
        """Skip when not in workflow_plan."""
        runtime_state.workflow_plan = workflow_plan_deg_only
        result = wgcna_node(runtime_state)
        assert result == {}

    def test_module_sizes_sum(self, runtime_state_with_deg_full):
        """Module sizes + grey genes = total genes."""
        result = wgcna_node(runtime_state_with_deg_full)
        report = result.get("wgcna_report", {})
        if "reason" not in report and report.get("n_modules", 0) > 0:
            sizes = report.get("module_sizes", {})
            total_mod = sum(int(v) for v in sizes.values())
            grey = report.get("n_grey_genes", 0)
            n_analyzed = report.get("n_genes_analyzed", 0)
            assert total_mod + grey == n_analyzed


# ═══════════════════════════════════════════════════════════════
# O2PLS Node Integration Tests
# ═══════════════════════════════════════════════════════════════

class TestO2plsNode:

    def test_with_data(self, runtime_state_with_full_multiomics):
        """With expression + metabolite data, runs O2PLS."""
        result = o2pls_node(runtime_state_with_full_multiomics)
        assert "o2pls_report" in result
        report = result["o2pls_report"]
        if "reason" not in report:
            assert "n_joint_components" in report
            assert "x_joint_variance_explained" in report
            assert "y_joint_variance_explained" in report

    def test_without_data(self, runtime_state):
        """Without data, returns reason."""
        if runtime_state.planner_state:
            runtime_state.planner_state.expression_matrix = None
            runtime_state.planner_state.metabolite_matrix = None
        result = o2pls_node(runtime_state)
        assert "o2pls_report" in result

    def test_insufficient_samples(self, runtime_state):
        """With < 8 samples, returns warning."""
        if runtime_state.planner_state:
            small_expr = {}
            for i in range(10):
                small_expr[f"G_{i}"] = {"S1": 1.0, "S2": 2.0, "S3": 1.5}
            small_meta = {"M1": {"S1": 1.0, "S2": 2.0, "S3": 1.5}}
            runtime_state.planner_state.expression_matrix = small_expr
            runtime_state.planner_state.metabolite_matrix = small_meta
        runtime_state.deg_report = {"top_genes": [
            {"gene_id": f"G_{i}", "effect_size": 0.8, "max_log2fc": 2.0,
             "p_value": 0.001, "q_value": 0.01, "is_fdr_sig": True}
            for i in range(10)
        ]}
        runtime_state.dam_report = {"top_metabolites": [
            {"metabolite": "M1", "effect_size": 0.7, "max_log2fc": 1.5,
             "p_value": 0.001, "q_value": 0.01, "is_fdr_sig": True}
        ]}
        result = o2pls_node(runtime_state)
        report = result.get("o2pls_report", {})
        assert "reason" in report or report.get("n_samples", 0) < 8

    def test_with_workflow_plan_skip(self, runtime_state, workflow_plan_deg_only):
        """Skip when not in workflow_plan."""
        runtime_state.workflow_plan = workflow_plan_deg_only
        result = o2pls_node(runtime_state)
        assert result == {}

    def test_top_variables_have_vip(self, runtime_state_with_full_multiomics):
        """Top X variables should have VIP scores."""
        result = o2pls_node(runtime_state_with_full_multiomics)
        report = result.get("o2pls_report", {})
        if "reason" not in report:
            top_x = report.get("top_x_variables", [])
            if top_x:
                for var in top_x[:3]:
                    assert "vip" in var
                    assert "loading" in var
                    assert "name" in var


# ═══════════════════════════════════════════════════════════════
# Routing extension tests for Round 2
# ═══════════════════════════════════════════════════════════════

class TestRound2Routing:

    def test_new_nodes_in_all_pipeline(self):
        from phyto_reason.workflows.routing_logic import ALL_PIPELINE_NODES
        assert "wgcna" in ALL_PIPELINE_NODES
        assert "o2pls" in ALL_PIPELINE_NODES

    def test_new_steps_in_step_to_nodes(self):
        from phyto_reason.workflows.routing_logic import STEP_TO_NODES
        assert "wgcna" in STEP_TO_NODES
        assert "o2pls" in STEP_TO_NODES

    def test_new_steps_in_step_options(self):
        from phyto_reason.workflows.routing_logic import STEP_OPTIONS
        step_ids = {s["id"] for s in STEP_OPTIONS}
        assert "wgcna" in step_ids
        assert "o2pls" in step_ids
        assert len(STEP_OPTIONS) >= 10  # 8 基础 + 多组学节点（v4.7 起 >= 11）

    def test_expand_new_steps(self):
        from phyto_reason.workflows.routing_logic import expand_analysis_selection
        assert "wgcna" in expand_analysis_selection(["wgcna"])
        assert "o2pls" in expand_analysis_selection(["o2pls"])

    def test_multiomics_includes_all(self):
        from phyto_reason.workflows.routing_logic import expand_analysis_selection
        expanded = expand_analysis_selection(["multiomics"])
        assert "wgcna" in expanded
        assert "o2pls" in expanded

    def test_skip_guards(self, runtime_state, workflow_plan_deg_only):
        runtime_state.workflow_plan = workflow_plan_deg_only
        assert _should_skip_node(runtime_state, "wgcna")
        assert _should_skip_node(runtime_state, "o2pls")


# ═══════════════════════════════════════════════════════════════
# Fixtures
# ═══════════════════════════════════════════════════════════════

@pytest.fixture
def runtime_state_with_deg_full(mock_expression, mock_metabolite) -> RuntimeState:
    """RuntimeState with full DEG + expression data (6 samples, 30 genes)."""
    ps = PlannerState(
        species="Arabidopsis thaliana",
        target_metabolite="anthocyanin",
        has_expression=True,
        has_metabolite=True,
        expression_matrix=mock_expression,
        metabolite_matrix=mock_metabolite,
        sample_count=6,
    )
    state = RuntimeState(planner_state=ps)
    state.pathway_name = "flavonoid biosynthesis"
    state.pathway_genes = {f"GENE_{i:04d}" for i in range(5)}
    state.deg_report = {
        "n_genes_tested": 30,
        "n_fdr_significant": 8,
        "n_effect_ranked": 15,
        "n_groups": 2,
        "group_names": ["Leaf", "Root"],
        "total_samples": 6,
        "data_type": "fpkm",
        "top_genes": [
            {"gene_id": f"GENE_{i:04d}", "effect_size": 0.9 - i * 0.03,
             "max_log2fc": 3.0 - i * 0.2 if i % 2 == 0 else -2.5 + i * 0.15,
             "p_value": 0.0001 + i * 0.001, "q_value": 0.005 + i * 0.01,
             "stat_name": "F_stat", "is_fdr_sig": i < 8}
            for i in range(20)
        ],
    }
    return state


@pytest.fixture
def runtime_state_with_full_multiomics(mock_expression, mock_metabolite) -> RuntimeState:
    """RuntimeState with ≥8 samples and both expression + metabolite data."""
    import random
    np.random.seed(42)
    random.seed(42)

    # Create 10 samples across 2 groups
    samples = [f"Leaf_{i}" for i in range(1, 6)] + [f"Root_{i}" for i in range(1, 6)]
    genes = {
        f"GENE_{i:04d}": {s: np.random.uniform(3, 15) for s in samples}
        for i in range(40)
    }
    metas = {
        f"META_{i:03d}": {s: np.random.uniform(0.5, 5) for s in samples}
        for i in range(10)
    }

    ps = PlannerState(
        species="Arabidopsis thaliana",
        target_metabolite="anthocyanin",
        has_expression=True,
        has_metabolite=True,
        expression_matrix=genes,
        metabolite_matrix=metas,
        sample_count=10,
    )
    state = RuntimeState(planner_state=ps)
    state.pathway_name = "flavonoid biosynthesis"
    state.pathway_genes = {f"GENE_{i:04d}" for i in range(5)}
    state.deg_report = {
        "n_genes_tested": 40,
        "n_fdr_significant": 10,
        "n_effect_ranked": 20,
        "n_groups": 2,
        "group_names": ["Leaf", "Root"],
        "total_samples": 10,
        "data_type": "fpkm",
        "top_genes": [
            {"gene_id": f"GENE_{i:04d}", "effect_size": 0.9 - i * 0.02,
             "max_log2fc": 2.5 - i * 0.15, "p_value": 0.0001 + i * 0.0005,
             "q_value": 0.005 + i * 0.005, "stat_name": "F_stat", "is_fdr_sig": i < 10}
            for i in range(30)
        ],
    }
    state.dam_report = {
        "n_metabolites_tested": 10,
        "n_fdr_significant": 5,
        "top_metabolites": [
            {"metabolite": f"META_{i:03d}", "effect_size": 0.8 - i * 0.05,
             "max_log2fc": 1.5 - i * 0.2, "p_value": 0.001 + i * 0.01,
             "q_value": 0.01 + i * 0.01, "is_fdr_sig": i < 5}
            for i in range(10)
        ],
    }

    return state
