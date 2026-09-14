"""
mechanistic_hypothesis.py — PhytoReason-Agent 核心科学对象 v3.0。

替代旧的 CandidateTF 作为系统最终输出。

代表一个"经证据支持的、可被反驳的、有竞争者的机理假设"。
不是"鉴定出的调控因子"。

Architecture Freeze Spec v3.0 定义。
"""

from __future__ import annotations

from enum import Enum

from pydantic import BaseModel, Field


class MechanismType(str, Enum):
    """代谢变化的可能机制类型。
    v3.0: 默认先排除非转录机制，最后才考虑转录调控。"""
    TRANSCRIPTIONAL = "transcriptional_regulation"
    ENZYMATIC = "enzymatic_regulation"
    TRANSPORT = "transport_redistribution"
    DEGRADATION = "degradation_alteration"
    STRESS_RESPONSE = "stress_induced_redistribution"
    DEVELOPMENTAL = "developmental_program"
    COMPOSITIONAL = "cell_composition_change"
    PRECURSOR_LIMITATION = "precursor_limitation"
    TECHNICAL_ARTIFACT = "technical_artifact"
    INSUFFICIENT_DATA = "insufficient_data"
    MIXED = "mixed_mechanisms"
    UNKNOWN = "unknown"


class EvidenceRelationType(str, Enum):
    """证据与假设之间的语义关系。"""
    SUPPORTS = "supports"
    CONTRADICTS = "contradicts"
    CONDITIONALLY_SUPPORTS = "conditional"
    INDIRECTLY_SUPPORTS = "indirect"
    SPECIES_LIMITED = "species_limited"
    TEMPORAL = "temporal"
    CONTEXT_DEPENDENT = "context_dependent"
    INFERRED_BY_HOMOLOGY = "inferred_homology"
    INSUFFICIENT = "insufficient"


class EvidenceRecord(BaseModel):
    """带语义关系的证据记录（v3.0 SemanticEvidence 兼容层）。"""
    source: str
    evidence_type: str
    score: float = 0.0
    relation_type: EvidenceRelationType = EvidenceRelationType.INSUFFICIENT
    description: str = ""
    context: str = ""
    species_scope: str = ""
    tissue_scope: str = ""
    condition_scope: str = ""
    contradiction_severity: str = ""  # high | medium | low（仅矛盾证据使用）
    metadata: dict = Field(default_factory=dict)


class CausalStep(BaseModel):
    """因果链中的一步。每个断言可独立验证或反驳。"""
    position: int
    subject: str
    predicate: str
    object: str
    evidence_for: list[str] = Field(default_factory=list)
    evidence_against: list[str] = Field(default_factory=list)
    testable: bool = False


class CausalChain(BaseModel):
    """显式的因果链 — 替代隐式 scoring。

    一个假设必须显式声明其因果逻辑。
    """
    chain_description: str = ""
    chain_steps: list[CausalStep] = Field(default_factory=list)
    known_gaps: list[str] = Field(default_factory=list)
    assumption_level: str = "high"  # low | moderate | high


class FalsificationTest(BaseModel):
    """一个可执行的证伪测试。

    每个假设必须提出：如果我是错的，应该观察到什么。
    """
    test_name: str
    prediction: str
    contradictory_prediction: str = ""
    required_data: str = ""
    feasibility: str = "unknown"  # high | medium | low | unknown


class UncertaintyModel(BaseModel):
    """假设的不确定性分解 — 不只是一个标签。"""
    overall_level: str = "high"  # low | moderate | high | insufficient
    statistical_uncertainty: float = 0.0
    evidence_gap_uncertainty: float = 0.0
    contradiction_uncertainty: float = 0.0
    assumption_uncertainty: float = 0.0
    details: list[str] = Field(default_factory=list)


class MechanisticHypothesis(BaseModel):
    """系统的核心科学输出 — 一个可比较、可反驳的机理假设。

    不是"找到了什么"，而是"在当前证据下什么机理模型最合理"。
    """
    id: str = ""
    title: str = ""

    mechanism_type: MechanismType = MechanismType.UNKNOWN

    upstream_signal: str | None = None

    proposed_regulators: list[str] = Field(default_factory=list)

    target_pathway: str = ""

    target_metabolites: list[str] = Field(default_factory=list)

    causal_chain: CausalChain = Field(default_factory=CausalChain)

    supporting_evidence: list[EvidenceRecord] = Field(default_factory=list)

    contradictory_evidence: list[EvidenceRecord] = Field(default_factory=list)

    missing_evidence: list[str] = Field(default_factory=list)

    uncertainty_level: str = "high"  # low | moderate | high | insufficient
    uncertainty_components: dict = Field(default_factory=dict)

    falsification_tests: list[FalsificationTest] = Field(default_factory=list)

    why_wrong: str = ""

    validation_experiments: list[str] = Field(default_factory=list)

    next_best_evidence: list[str] = Field(default_factory=list)

    def add_support(self, record: EvidenceRecord) -> None:
        self.supporting_evidence.append(record)

    def add_contradiction(self, record: EvidenceRecord) -> None:
        self.contradictory_evidence.append(record)

    def to_diagnostic(self) -> str:
        """输出当前假设的诊断信息（非结论性措辞）。"""
        parts = [f"Hypothesis: {self.title}"]
        parts.append(f"  Type: {self.mechanism_type.value}")
        parts.append(f"  Uncertainty: {self.uncertainty_level}")
        parts.append(f"  Supports: {len(self.supporting_evidence)} sources")
        parts.append(f"  Contradictions: {len(self.contradictory_evidence)}")
        parts.append(f"  Missing evidence: {len(self.missing_evidence)} items")
        if self.falsification_tests:
            parts.append(f"  Falsification tests: {len(self.falsification_tests)} proposed")
        return "\n".join(parts)
