"""
conftest.py — Shared pytest fixtures for PhytoReason-Agent tests.
"""

from __future__ import annotations

import random
from unittest.mock import patch, MagicMock

import numpy as np
import pytest

from phyto_reason.workflows.runtime_state import RuntimeState
from phyto_reason.models.planner_state import PlannerState
from phyto_reason.models.workflow_plan import WorkflowPlan
from phyto_reason.models.candidate_gene import CandidateGene
from phyto_reason.models.evidence import Evidence, EvidenceType
from phyto_reason.fusion.fusion_result import UnifiedFusionResult, ConfidenceLevel


# ── Deterministic random seed for all tests ────────────────

@pytest.fixture(autouse=True)
def set_seed():
    """Set deterministic seed for all tests."""
    random.seed(42)
    np.random.seed(42)
    yield


# ── Mock LLM client ────────────────────────────────────────

@pytest.fixture
def mock_llm_client():
    """Patch LLMClient to return controlled responses."""
    with patch("phyto_reason.llm.llm_client.LLMClient") as mock_cls:
        mock_instance = MagicMock()
        mock_instance.is_available.return_value = True

        # Structured output for _assess()
        mock_instance.chat_structured.return_value = {
            "research_needed": True,
            "analysis_needed": True,
            "species": "Arabidopsis thaliana",
            "metabolite": "berberine",
            "pathway": "",
        }

        # Tool-use chat for agents
        mock_instance.chat_with_tools.return_value = (
            "Analysis complete.",
            [],  # no tool calls
        )

        # Plain chat fallback
        mock_instance.chat.return_value = "Mock response."

        mock_cls.return_value = mock_instance
        yield mock_instance


# ── Mock expression and metabolite matrices ────────────────

@pytest.fixture
def mock_expression() -> dict:
    """Expression matrix: 30 genes x 6 samples, 2 tissue groups."""
    return {
        f"GENE_{i:04d}": {
            "Leaf_1": random.uniform(5, 15),
            "Leaf_2": random.uniform(5, 15),
            "Leaf_3": random.uniform(5, 15),
            "Root_1": random.uniform(5, 15),
            "Root_2": random.uniform(5, 15),
            "Root_3": random.uniform(5, 15),
        }
        for i in range(30)
    }


@pytest.fixture
def mock_metabolite() -> dict:
    """Metabolite matrix: 10 metabolites x 6 samples."""
    return {
        f"META_{i:03d}": {
            "Leaf_1": random.uniform(0.5, 5),
            "Leaf_2": random.uniform(0.5, 5),
            "Leaf_3": random.uniform(0.5, 5),
            "Root_1": random.uniform(0.5, 5),
            "Root_2": random.uniform(0.5, 5),
            "Root_3": random.uniform(0.5, 5),
        }
        for i in range(10)
    }


# ── Mock SessionData ───────────────────────────────────────

@pytest.fixture
def mock_session():
    """Mock SessionData with expression and metabolite matrices."""
    from phyto_reason.agent.conversation_store import store
    session = store.get_or_create("test-session-orch")
    session.expression_matrix = {
        f"GENE_{i:04d}": {
            "Leaf_1": 10.0, "Leaf_2": 12.0, "Leaf_3": 11.0,
            "Root_1": 5.0, "Root_2": 6.0, "Root_3": 5.5,
        }
        for i in range(30)
    }
    session.metabolite_matrix = {
        "berberine": {"Leaf_1": 1.0, "Leaf_2": 1.2, "Leaf_3": 1.1,
                      "Root_1": 4.0, "Root_2": 4.2, "Root_3": 3.8},
    }
    session.species = "Arabidopsis thaliana"
    session.target_metabolite = "berberine"
    session.has_data = True
    return session


# ── RuntimeState prefilled with data ───────────────────────

@pytest.fixture
def runtime_state(mock_expression, mock_metabolite) -> RuntimeState:
    """RuntimeState with planner_state populated."""
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
    return state


@pytest.fixture
def runtime_state_with_deg(runtime_state) -> RuntimeState:
    """RuntimeState with pre-populated DEG report."""
    runtime_state.deg_report = {
        "n_genes_tested": 30,
        "n_fdr_significant": 5,
        "n_effect_ranked": 10,
        "n_groups": 2,
        "group_names": ["Leaf", "Root"],
        "top_genes": [
            {"gene_id": f"GENE_{i:04d}", "effect_size": 0.8 - i * 0.02,
             "max_log2fc": 2.0 - i * 0.1, "p_value": 0.001, "q_value": 0.01 + i * 0.005,
             "stat_name": "F_stat", "is_fdr_sig": i < 5}
            for i in range(10)
        ],
    }
    return runtime_state


# ── WorkflowPlan fixtures ──────────────────────────────────

@pytest.fixture
def workflow_plan_deg_only() -> WorkflowPlan:
    return WorkflowPlan(
        selected_workflow="custom",
        tools_to_run=[
            "data_quality_gate", "metabolite_validation", "pathway_coherence",
            "mechanism_classifier", "deg_analysis",
        ],
        original_selection=["DEG"],
        reasoning="User selected DEG only",
    )


@pytest.fixture
def workflow_plan_full() -> WorkflowPlan:
    from phyto_reason.workflows.routing_logic import ALL_PIPELINE_NODES
    return WorkflowPlan(
        selected_workflow="full",
        tools_to_run=sorted(ALL_PIPELINE_NODES),
        original_selection=["all"],
        reasoning="Full pipeline",
    )


@pytest.fixture
def workflow_plan_hypothesis_only() -> WorkflowPlan:
    return WorkflowPlan(
        selected_workflow="custom",
        tools_to_run=[
            "data_quality_gate", "metabolite_validation", "pathway_coherence",
            "mechanism_classifier", "regulation_evidence", "contradiction_check",
            "falsification", "hypothesis_competition", "hypothesis_synthesis",
        ],
        original_selection=["hypothesis"],
        reasoning="User selected hypothesis only",
    )


# ── Pre-built reports ──────────────────────────────────────

@pytest.fixture
def deg_report() -> dict:
    return {
        "n_genes_tested": 100,
        "n_fdr_significant": 15,
        "n_effect_ranked": 50,
        "n_groups": 3,
        "group_names": ["Leaf", "Root", "Stem"],
        "total_samples": 9,
        "data_type": "fpkm",
        "top_genes": [
            {"gene_id": "GENE_0001", "effect_size": 0.95, "max_log2fc": 3.2,
             "p_value": 0.00001, "q_value": 0.001, "is_fdr_sig": True},
            {"gene_id": "GENE_0002", "effect_size": 0.88, "max_log2fc": 2.8,
             "p_value": 0.00005, "q_value": 0.003, "is_fdr_sig": True},
        ],
    }


@pytest.fixture
def tf_candidates() -> dict:
    return {
        "TF_MYB": CandidateGene(
            gene_id="TF_MYB", is_tf=True, correlation_score=0.85,
            motif_score=0.7, target_enzymes=["ENZYME_A"],
            target_metabolites=["berberine"],
        ),
        "TF_WRKY": CandidateGene(
            gene_id="TF_WRKY", is_tf=True, correlation_score=0.72,
            motif_score=0.4, target_enzymes=["ENZYME_B"],
            target_metabolites=["berberine"],
        ),
    }


@pytest.fixture
def temp_db_path(tmp_path) -> str:
    """Temporary path for SQLite session DB."""
    return str(tmp_path / "test_sessions.db")
