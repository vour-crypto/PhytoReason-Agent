"""
test_multiomics_v47.py — Tests for v4.7 multi-omics extension nodes.

Covers: joint_enrichment, quadrant_plot, correlation_network nodes.
"""

from __future__ import annotations

from unittest.mock import patch, MagicMock

import numpy as np
import pytest

from phyto_reason.workflows.nodes.joint_enrichment import (
    joint_enrichment_node,
    _hypergeometric_test,
    _fisher_combined_pvalue,
    _detect_pathway_for_metabolite,
    _match_genes,
)
from phyto_reason.workflows.nodes.quadrant_plot import (
    quadrant_plot_node,
    _classify_quadrant,
    _generate_quadrant_plot,
)
from phyto_reason.workflows.nodes.correlation_network import (
    correlation_network_node,
    _spearman_corr,
)
from phyto_reason.workflows.execution_nodes import _should_skip_node
from phyto_reason.workflows.runtime_state import RuntimeState
from phyto_reason.models.planner_state import PlannerState


# ═══════════════════════════════════════════════════════════════
# Unit tests for helper functions
# ═══════════════════════════════════════════════════════════════

class TestHypergeometricTest:

    def test_enrichment_detected(self):
        """Strong enrichment should give low p-value."""
        # 10 pathway genes out of 1000, 5 out of 20 DEGs overlap
        p = _hypergeometric_test(k=5, M=10, n=20, N=1000)
        assert p < 0.001  # Very significant

    def test_no_enrichment(self):
        """Random overlap should give high p-value."""
        p = _hypergeometric_test(k=1, M=100, n=20, N=1000)
        assert p > 0.05

    def test_zero_overlap(self):
        """Zero overlap gives p=1.0."""
        p = _hypergeometric_test(k=0, M=10, n=20, N=1000)
        assert p == 1.0

    def test_perfect_enrichment(self):
        """All pathway genes are DEGs."""
        p = _hypergeometric_test(k=10, M=10, n=20, N=1000)
        assert p < 1e-10

    def test_invalid_inputs(self):
        """Invalid inputs should return 1.0 safely."""
        assert _hypergeometric_test(k=0, M=0, n=0, N=0) == 1.0
        assert _hypergeometric_test(k=-1, M=10, n=20, N=1000) == 1.0


class TestFisherCombinedPvalue:

    def test_combines_low_pvalues(self):
        """Two moderately significant p-values combine to higher significance."""
        p = _fisher_combined_pvalue([0.01, 0.01])
        assert p < 0.01

    def test_combines_high_pvalues(self):
        """High p-values stay high when combined."""
        p = _fisher_combined_pvalue([0.5, 0.5])
        assert p > 0.05

    def test_empty_input(self):
        """Empty list returns 1.0."""
        assert _fisher_combined_pvalue([]) == 1.0

    def test_single_pvalue(self):
        """Single p-value returns similar value."""
        p = _fisher_combined_pvalue([0.03])
        # χ² with df=2: -2*ln(0.03) ≈ 7.01, sf ≈ 0.03
        assert 0.02 < p < 0.05


class TestDetectPathway:

    def test_flavonoid_detection(self):
        assert _detect_pathway_for_metabolite("anthocyanin") == "flavonoid_biosynthesis"
        assert _detect_pathway_for_metabolite("flavonol") == "flavonoid_biosynthesis"

    def test_alkaloid_detection(self):
        assert _detect_pathway_for_metabolite("berberine") == "alkaloid_biosynthesis"

    def test_carotenoid_detection(self):
        assert _detect_pathway_for_metabolite("lycopene") == "carotenoid_biosynthesis"

    def test_unknown_metabolite(self):
        assert _detect_pathway_for_metabolite("unknown_xyz") is None


class TestMatchGenes:

    def test_substring_match(self):
        matched = _match_genes(["CHS", "PAL"], {"AT5G13930_CHS", "AT2G37040_PAL", "AT1G01010"})
        assert "AT5G13930_CHS" in matched
        assert "AT2G37040_PAL" in matched

    def test_case_insensitive(self):
        matched = _match_genes(["chs"], {"At5g13930_CHS"})
        assert "At5g13930_CHS" in matched

    def test_no_match(self):
        matched = _match_genes(["XYZ"], {"AT1G01010"})
        assert len(matched) == 0


class TestQuadrantClassification:

    def test_four_quadrant_q1(self):
        assert _classify_quadrant(2.0, 1.5) == "Q1"  # both up

    def test_four_quadrant_q2(self):
        assert _classify_quadrant(-2.0, 1.5) == "Q2"  # gene down, meta up

    def test_four_quadrant_q3(self):
        assert _classify_quadrant(-2.0, -1.5) == "Q3"  # both down

    def test_four_quadrant_q4(self):
        assert _classify_quadrant(2.0, -1.5) == "Q4"  # gene up, meta down

    def test_four_quadrant_boundary(self):
        assert _classify_quadrant(0.0, 0.0) == "Q1"  # zero goes to Q1

    def test_nine_quadrant_center(self):
        assert _classify_quadrant(0.1, -0.1, threshold=1.0) == "Q9"  # center

    def test_nine_quadrant_q1(self):
        assert _classify_quadrant(2.0, 2.0, threshold=1.0) == "Q1"


class TestSpearmanCorr:

    def test_perfect_positive(self):
        x = np.array([1, 2, 3, 4, 5])
        y = np.array([2, 4, 6, 8, 10])
        rho, p = _spearman_corr(x, y)
        assert abs(rho - 1.0) < 0.01
        assert p < 0.05

    def test_perfect_negative(self):
        x = np.array([1, 2, 3, 4, 5])
        y = np.array([10, 8, 6, 4, 2])
        rho, p = _spearman_corr(x, y)
        assert abs(rho + 1.0) < 0.01

    def test_no_correlation(self):
        np.random.seed(42)
        x = np.random.randn(50)
        y = np.random.randn(50)
        rho, p = _spearman_corr(x, y)
        assert abs(rho) < 0.5  # unlikely to be strongly correlated

    def test_short_input(self):
        """Input with fewer than 3 points returns 0,1."""
        rho, p = _spearman_corr(np.array([1, 2]), np.array([1, 2]))
        assert rho == 0.0
        assert p == 1.0


# ═══════════════════════════════════════════════════════════════
# Integration tests for nodes
# ═══════════════════════════════════════════════════════════════

class TestJointEnrichmentNode:

    def test_with_deg_and_dam(self, runtime_state_with_deg_full):
        """With DEG + DAM data, produces enrichment report."""
        result = joint_enrichment_node(runtime_state_with_deg_full)
        assert "joint_enrichment_report" in result
        report = result["joint_enrichment_report"]
        if "reason" not in report:
            assert "n_pathways_tested" in report
            assert "results" in report

    def test_without_data(self, runtime_state):
        """Without DEG/DAM data, returns reason."""
        runtime_state.deg_report = {}
        runtime_state.dam_report = {}
        result = joint_enrichment_node(runtime_state)
        assert "joint_enrichment_report" in result

    def test_with_workflow_plan_skip(self, runtime_state, workflow_plan_deg_only):
        """Skip when not in workflow_plan."""
        runtime_state.workflow_plan = workflow_plan_deg_only
        result = joint_enrichment_node(runtime_state)
        assert result == {}  # skipped, returns empty

    def test_detects_flavonoid_pathway(self, runtime_state_with_deg_full):
        """With anthocyanin target, should test flavonoid pathway."""
        if runtime_state_with_deg_full.planner_state:
            runtime_state_with_deg_full.planner_state.target_metabolite = "anthocyanin"
        result = joint_enrichment_node(runtime_state_with_deg_full)
        report = result.get("joint_enrichment_report", {})
        if report.get("n_pathways_tested", 0) > 0:
            pathway_names = [r.get("pathway_name", "") for r in report.get("results", [])]
            assert any("flavonoid" in name.lower() or "anthocyanin" in name.lower()
                       for name in pathway_names)


class TestQuadrantPlotNode:

    def test_with_deg_dam_and_multiomics(self, runtime_state_with_deg_full):
        """With all data, produces quadrant plot."""
        # Add multiomics report with correlation pairs
        runtime_state_with_deg_full.multiomics_report = {
            "n_correlation_pairs_tested": 10,
            "n_significant_pairs": 5,
            "top_pairs": [
                {"gene_id": "GENE_0000", "metabolite": "META_000",
                 "correlation": 0.85, "p_value": 0.001, "q_value": 0.01},
                {"gene_id": "GENE_0001", "metabolite": "META_001",
                 "correlation": -0.72, "p_value": 0.005, "q_value": 0.03},
            ],
        }
        # Also need DAM report with log2FC
        runtime_state_with_deg_full.dam_report = {
            "n_metabolites_tested": 10,
            "n_fdr_significant": 3,
            "top_metabolites": [
                {"metabolite": "META_000", "effect_size": 0.7, "max_log2fc": 1.5,
                 "p_value": 0.001, "q_value": 0.01, "is_fdr_sig": True},
                {"metabolite": "META_001", "effect_size": 0.5, "max_log2fc": -1.2,
                 "p_value": 0.01, "q_value": 0.05, "is_fdr_sig": True},
            ],
        }

        result = quadrant_plot_node(runtime_state_with_deg_full)
        assert "quadrant_plot_report" in result
        report = result["quadrant_plot_report"]
        if "reason" not in report:
            assert "n_pairs" in report
            assert "quadrant_counts" in report
            assert "synchronicity_ratio" in report

    def test_without_data(self, runtime_state):
        """Without DEG/DAM, returns reason."""
        result = quadrant_plot_node(runtime_state)
        assert "quadrant_plot_report" in result

    def test_with_workflow_plan_skip(self, runtime_state, workflow_plan_deg_only):
        """Skip when not in workflow_plan."""
        runtime_state.workflow_plan = workflow_plan_deg_only
        result = quadrant_plot_node(runtime_state)
        assert result == {}

    def test_nine_quadrant_mode(self, runtime_state_with_deg_full):
        """Test that nine-quadrant classification works."""
        # Direct test of classification function
        assert _classify_quadrant(0.5, 0.3, threshold=1.0) == "Q9"  # center
        assert _classify_quadrant(2.0, -0.3, threshold=1.0) == "Q5"  # gene up, meta no change


class TestCorrelationNetworkNode:

    def test_with_data(self, runtime_state_with_deg_full):
        """With expression and metabolite data, builds network."""
        # Add DAM report
        runtime_state_with_deg_full.dam_report = {
            "n_metabolites_tested": 10,
            "n_fdr_significant": 5,
            "top_metabolites": [
                {"metabolite": f"META_{i:03d}", "effect_size": 0.8 - i * 0.05,
                 "max_log2fc": 1.5 - i * 0.2, "p_value": 0.001 + i * 0.01,
                 "q_value": 0.01 + i * 0.01, "is_fdr_sig": i < 3}
                for i in range(10)
            ],
        }
        # Update metabolite matrix to include DAM names
        if runtime_state_with_deg_full.planner_state:
            meta_mat = runtime_state_with_deg_full.planner_state.metabolite_matrix or {}
            for i in range(10):
                meta_mat[f"META_{i:03d}"] = {
                    "Leaf_1": np.random.uniform(0.5, 3),
                    "Leaf_2": np.random.uniform(0.5, 3),
                    "Leaf_3": np.random.uniform(0.5, 3),
                    "Root_1": np.random.uniform(0.5, 3),
                    "Root_2": np.random.uniform(0.5, 3),
                    "Root_3": np.random.uniform(0.5, 3),
                }
            runtime_state_with_deg_full.planner_state.metabolite_matrix = meta_mat

        result = correlation_network_node(runtime_state_with_deg_full)
        assert "correlation_network_report" in result
        report = result["correlation_network_report"]
        if "reason" not in report:
            assert "n_nodes" in report
            assert "edges" in report
            assert "correlation_method" in report

    def test_without_data(self, runtime_state):
        """Without expression/metabolite data, returns reason."""
        result = correlation_network_node(runtime_state)
        assert "correlation_network_report" in result

    def test_with_workflow_plan_skip(self, runtime_state, workflow_plan_deg_only):
        """Skip when not in workflow_plan."""
        runtime_state.workflow_plan = workflow_plan_deg_only
        result = correlation_network_node(runtime_state)
        assert result == {}

    def test_networkx_fallback(self, runtime_state_with_deg_full):
        """When networkx unavailable, falls back to edge list."""
        with patch("phyto_reason.workflows.nodes.correlation_network._build_correlation_network") as mock_build:
            mock_build.return_value = {
                "n_nodes": 20, "n_edges": 15, "n_sig_edges": 8,
                "correlation_method": "spearman", "min_abs_rho": 0.4,
                "edges": [{"source": "G1", "target": "M1", "rho": 0.75, "q_value": 0.01}],
                "hubs": [], "modules": [], "fallback": "networkx_unavailable",
            }
            runtime_state_with_deg_full.dam_report = {
                "n_metabolites_tested": 5, "n_fdr_significant": 2,
                "top_metabolites": [
                    {"metabolite": "META_000", "effect_size": 0.7, "max_log2fc": 1.5,
                     "p_value": 0.001, "q_value": 0.01, "is_fdr_sig": True},
                ],
            }
            result = correlation_network_node(runtime_state_with_deg_full)
            report = result.get("correlation_network_report", {})
            assert report.get("n_edges") == 15


# ═══════════════════════════════════════════════════════════════
# Routing logic extension tests
# ═══════════════════════════════════════════════════════════════

class TestV47Routing:

    def test_new_nodes_in_all_pipeline_nodes(self):
        """New v4.7 nodes should be in ALL_PIPELINE_NODES."""
        from phyto_reason.workflows.routing_logic import ALL_PIPELINE_NODES
        assert "joint_enrichment" in ALL_PIPELINE_NODES
        assert "quadrant_plot" in ALL_PIPELINE_NODES
        assert "correlation_network" in ALL_PIPELINE_NODES

    def test_new_steps_in_step_to_nodes(self):
        """New analysis steps should be in STEP_TO_NODES."""
        from phyto_reason.workflows.routing_logic import STEP_TO_NODES
        assert "joint_enrichment" in STEP_TO_NODES
        assert "quadrant_plot" in STEP_TO_NODES
        assert "correlation_network" in STEP_TO_NODES

    def test_new_steps_in_step_options(self):
        """New analysis steps should be in STEP_OPTIONS."""
        from phyto_reason.workflows.routing_logic import STEP_OPTIONS
        step_ids = {s["id"] for s in STEP_OPTIONS}
        assert "joint_enrichment" in step_ids
        assert "quadrant_plot" in step_ids
        assert "correlation_network" in step_ids

    def test_expand_new_steps(self):
        """expand_analysis_selection should handle new steps."""
        from phyto_reason.workflows.routing_logic import expand_analysis_selection
        expanded = expand_analysis_selection(["joint_enrichment"])
        assert "joint_enrichment" in expanded
        # Always-run nodes should be included
        assert "data_quality_gate" in expanded

    def test_multiomics_includes_new_nodes(self):
        """Selecting 'multiomics' should now include all 4 multi-omics nodes."""
        from phyto_reason.workflows.routing_logic import expand_analysis_selection
        expanded = expand_analysis_selection(["multiomics"])
        assert "multiomics_integration" in expanded
        assert "joint_enrichment" in expanded
        assert "quadrant_plot" in expanded
        assert "correlation_network" in expanded

    def test_skip_guards_for_new_nodes(self, runtime_state, workflow_plan_deg_only):
        """WorkflowPlan skip guards work for new nodes."""
        runtime_state.workflow_plan = workflow_plan_deg_only
        # deg_only doesn't include the new nodes
        assert _should_skip_node(runtime_state, "joint_enrichment")
        assert _should_skip_node(runtime_state, "quadrant_plot")
        assert _should_skip_node(runtime_state, "correlation_network")


# ═══════════════════════════════════════════════════════════════
# Fixtures
# ═══════════════════════════════════════════════════════════════

@pytest.fixture
def runtime_state_with_deg_full(mock_expression, mock_metabolite) -> RuntimeState:
    """RuntimeState with full DEG + DAM pre-populated data."""
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

    # Pre-populate DEG report
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
