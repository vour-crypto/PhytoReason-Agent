"""
test_models.py — Unit tests for data models.

Validates Pydantic model construction, defaults, and serialization.
No external data needed.
"""

import pytest

from phyto_reason.models.evidence import Evidence, EvidenceType, EvidenceRelationType
from phyto_reason.models.candidate_gene import CandidateGene
from phyto_reason.models.mechanistic_hypothesis import (
    MechanisticHypothesis, MechanismType, EvidenceRecord, CausalChain, FalsificationTest,
)
from phyto_reason.models.workflow_plan import WorkflowPlan
from phyto_reason.models.tool_result import ToolResult
from phyto_reason.models.enums import ConfidenceLevel, BiologicalPlausibility


# ═══════════════════════════════════════════════════════════════
# Evidence
# ═══════════════════════════════════════════════════════════════

class TestEvidence:
    def test_create_minimal(self):
        ev = Evidence(evidence_type=EvidenceType.CORRELATION)
        assert ev.evidence_type == EvidenceType.CORRELATION
        assert ev.score == 0.0
        assert ev.confidence == 1.0

    def test_create_full(self):
        ev = Evidence(
            evidence_type=EvidenceType.MOTIF_BINDING,
            source="FIMO scan",
            score=0.85,
            confidence=0.9,
            description="E-box motif found in BBE promoter",
            relation_type=EvidenceRelationType.SUPPORTS,
        )
        assert ev.score == 0.85
        assert ev.relation_type == EvidenceRelationType.SUPPORTS

    def test_score_range_enforced(self):
        """Score must be 0-1 (Pydantic validation)."""
        with pytest.raises(Exception):
            Evidence(evidence_type=EvidenceType.CORRELATION, score=1.5)
        with pytest.raises(Exception):
            Evidence(evidence_type=EvidenceType.CORRELATION, score=-0.1)

    def test_to_summary(self):
        ev = Evidence(evidence_type=EvidenceType.CORRELATION, source="bicor", score=0.78)
        summary = ev.to_summary()
        assert "correlation" in summary.lower()
        assert "bicor" in summary

    def test_all_evidence_types_exist(self):
        """All expected evidence types should be in the enum."""
        types = {e.value for e in EvidenceType}
        assert "correlation" in types
        assert "motif_binding" in types
        assert "literature" in types
        assert "tf_family_prior" in types


class TestEvidenceRelationType:
    def test_all_relations_exist(self):
        relations = {e.value for e in EvidenceRelationType}
        assert "supports" in relations
        assert "contradicts" in relations

    def test_default_is_supports(self):
        er = EvidenceRelationType.SUPPORTS
        assert er == "supports"


# ═══════════════════════════════════════════════════════════════
# CandidateGene
# ═══════════════════════════════════════════════════════════════

class TestCandidateGene:
    def test_create_minimal(self):
        cg = CandidateGene(gene_id="test_gene_1")
        assert cg.gene_id == "test_gene_1"
        assert cg.confidence_score == 0.0

    def test_add_evidence(self):
        cg = CandidateGene(gene_id="MYB1")
        ev = Evidence(evidence_type=EvidenceType.CORRELATION, score=0.7)
        cg.evidence_list.append(ev)
        assert len(cg.evidence_list) == 1
        assert cg.evidence_list[0].score == 0.7


# ═══════════════════════════════════════════════════════════════
# MechanisticHypothesis (v4.0 core output)
# ═══════════════════════════════════════════════════════════════

class TestMechanisticHypothesis:
    def test_create_minimal(self):
        h = MechanisticHypothesis(
            id="H1",
            title="Test hypothesis",
            mechanism_type=MechanismType.TRANSCRIPTIONAL,
        )
        assert h.id == "H1"
        assert h.mechanism_type == MechanismType.TRANSCRIPTIONAL
        assert h.uncertainty_level == "high"  # default

    def test_add_supporting_evidence(self):
        h = MechanisticHypothesis(id="H1", title="Test", mechanism_type=MechanismType.TRANSCRIPTIONAL)
        ev = EvidenceRecord(
            source="correlation",
            evidence_type="correlation",
            score=0.85,
            relation_type=EvidenceRelationType.SUPPORTS,
            description="bicor r=0.78",
        )
        h.supporting_evidence.append(ev)
        assert len(h.supporting_evidence) == 1
        assert h.supporting_evidence[0].score == 0.85

    def test_add_contradictory_evidence(self):
        h = MechanisticHypothesis(id="H1", title="Test", mechanism_type=MechanismType.TRANSCRIPTIONAL)
        ev = EvidenceRecord(
            source="motif_scan",
            evidence_type="motif_binding",
            score=0.3,
            relation_type=EvidenceRelationType.CONTRADICTS,
            description="No motif found",
            contradiction_severity="high",
        )
        h.contradictory_evidence.append(ev)
        assert len(h.contradictory_evidence) == 1
        assert h.contradictory_evidence[0].contradiction_severity == "high"

    def test_falsification_test(self):
        ft = FalsificationTest(
            test_name="Knockout experiment",
            prediction="If MYB1 is knocked out, metabolite levels should drop",
            contradictory_prediction="If MYB1 is not causal, levels unchanged",
        )
        h = MechanisticHypothesis(
            id="H1", title="Test",
            mechanism_type=MechanismType.TRANSCRIPTIONAL,
            falsification_tests=[ft],
        )
        assert len(h.falsification_tests) == 1
        assert h.falsification_tests[0].test_name == "Knockout experiment"

    def test_all_mechanism_types_exist(self):
        """All mechanism types from the spec should be in the enum."""
        types = {e.value for e in MechanismType}
        expected = [
            "transcriptional_regulation", "enzymatic_regulation",
            "transport_redistribution", "stress_induced_redistribution",
            "technical_artifact", "unknown",
        ]
        for t in expected:
            assert t in types, f"Missing mechanism type: {t}"


class TestCausalChain:
    def test_create_chain(self):
        from phyto_reason.models.mechanistic_hypothesis import CausalStep
        chain = CausalChain(
            chain_description="MYB1 binds DFR promoter → activates anthocyanin synthesis",
            chain_steps=[
                CausalStep(
                    position=1, subject="MYB1", predicate="binds",
                    object="DFR promoter",
                ),
                CausalStep(
                    position=2, subject="DFR", predicate="catalyzes",
                    object="anthocyanin synthesis",
                ),
            ],
            assumption_level="moderate",
        )
        assert len(chain.chain_steps) == 2
        assert chain.chain_steps[0].subject == "MYB1"


# ═══════════════════════════════════════════════════════════════
# WorkflowPlan
# ═══════════════════════════════════════════════════════════════

class TestWorkflowPlan:
    def test_create_default(self):
        wp = WorkflowPlan()
        assert wp.selected_workflow == ""
        assert wp.tools_to_run == []
        assert wp.tools_to_skip == []

    def test_create_with_tools(self):
        wp = WorkflowPlan(
            selected_workflow="full_analysis",
            tools_to_run=["correlation", "motif_scan"],
            tools_to_skip=["wgcna"],
            reasoning="User has expression data",
        )
        assert len(wp.tools_to_run) == 2
        assert "wgcna" in wp.tools_to_skip


# ═══════════════════════════════════════════════════════════════
# ToolResult
# ═══════════════════════════════════════════════════════════════

class TestToolResult:
    def test_create_success(self):
        tr = ToolResult(status="success")
        assert tr.status == "success"

    def test_create_failure(self):
        tr = ToolResult(status="failed", warnings=["FIMO not found"])
        assert tr.status == "failed"
        assert "FIMO" in tr.warnings[0]

    def test_merge(self):
        tr1 = ToolResult(status="success", candidates=[], evidence_list=[])
        tr2 = ToolResult(status="failed", warnings=["error"])
        tr1.merge(tr2)
        assert tr1.status == "partial"


# ═══════════════════════════════════════════════════════════════
# Enums
# ═══════════════════════════════════════════════════════════════

class TestEnums:
    def test_confidence_levels(self):
        levels = {e.value for e in ConfidenceLevel}
        assert "Gold" in levels
        assert "Weak" in levels

    def test_biological_plausibility(self):
        bp = {e.value for e in BiologicalPlausibility}
        assert "plausible" in bp
        assert "weak" in bp


# ── 矩阵加载 GBK 回退（真实数据暴露的问题）───────────────

def test_matrix_loader_gbk_encoding(tmp_path):
    """代谢物中文名 CSV 常为 GBK 编码——必须能加载。"""
    import pandas as pd
    from phyto_reason.ingestion.preprocessing.matrix_loader import MatrixLoader
    df = pd.DataFrame({"S1": [1.0, 2.0], "S2": [3.0, 4.0]}, index=["小檗碱", "黄连素"])
    path = tmp_path / "meta_gbk.csv"
    df.to_csv(path, encoding="gbk")
    matrix = MatrixLoader.load_metabolite_matrix(path)
    assert "小檗碱" in matrix
    assert matrix["小檗碱"]["S1"] == 1.0
