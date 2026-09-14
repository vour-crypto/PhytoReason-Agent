"""
analysis_service.py — 分析服务层。

唯一直接调用 WorkflowRunner 的模块。
API routes 只能通过此服务调用 workflow。
"""

from __future__ import annotations

import logging
from datetime import datetime
from typing import Any

from phyto_reason.api.schemas.request_models import AnalysisRequest
from phyto_reason.api.schemas.response_models import AnalysisResponse, CandidateInfo
from phyto_reason.workflows.workflow_runner import WorkflowRunner

logger = logging.getLogger("analysis_service")


class AnalysisService:
    """分析服务 — API 与 Workflow 之间的桥梁。"""

    def __init__(self) -> None:
        self.runner = WorkflowRunner()

    def run_analysis(self, request: AnalysisRequest) -> AnalysisResponse:
        """执行完整分析流程。"""
        try:
            logger.info(f"Starting analysis: species={request.species}, "
                        f"target={request.target_metabolite}")

            state = self.runner.run(
                research_question=request.research_question,
                species=request.species,
                target_metabolite=request.target_metabolite,
                target_pathway=request.target_pathway,
                sample_count=request.sample_count,
                has_expression=request.has_expression,
                has_metabolite=request.has_metabolite,
                has_annotation=request.has_annotation,
                has_promoter=request.has_promoter,
                has_itak=request.has_itak,
                data_path=request.data_path,
            )

            # 构建响应
            candidates = self._build_candidates(state)
            trace = self._build_trace(state)
            validation = self._build_validation(state)

            return AnalysisResponse(
                success=not state.error,
                session_id=state.session_id,
                workflow_summary=state.to_summary(),
                top_candidates=candidates,
                confidence={
                    "label": state.evaluation.label if state.evaluation else "unknown",
                    "score": state.evaluation.overall_score if state.evaluation else 0.0,
                    "explanation": state.evaluation.explanation if state.evaluation else "",
                } if state.evaluation else {},
                hypothesis=state.hypothesis,
                suggested_validation=validation,
                execution_trace=trace,
                error=state.error,
            )

        except Exception as e:
            logger.error(f"Analysis failed: {e}", exc_info=True)
            return AnalysisResponse(
                success=False,
                error=f"分析执行失败: {e}",
                execution_trace=[{"error": str(e)}],
            )

    @staticmethod
    def _build_candidates(state: Any) -> list[CandidateInfo]:
        """从 RuntimeState 提取候选列表。"""
        candidates: list[CandidateInfo] = []

        fusion_results = getattr(state, "fusion_results", None) or {}
        for gene_id, result in fusion_results.items():
            cand = CandidateInfo(
                gene_id=gene_id,
                confidence_score=result.calibrated_score,
                confidence_level=result.confidence_level.value if hasattr(result.confidence_level, 'value') else str(result.confidence_level),
                n_evidences=len(result.evidence_list),
                top_scores={ev.source: ev.normalized_score for ev in result.evidence_list[:5]},
            )
            candidates.append(cand)

        candidates.sort(key=lambda c: c.confidence_score, reverse=True)
        return candidates[:20]

    @staticmethod
    def _build_trace(state: Any) -> list[dict]:
        """从 RuntimeState 提取执行跟踪。"""
        trace = getattr(state, "execution_trace", None) or []
        return [{
            "node": step.node,
            "status": step.status,
            "error": step.error,
        } for step in trace]

    @staticmethod
    def _build_validation(state: Any) -> list[str]:
        """从 RuntimeState 提取验证建议。"""
        suggestions = []

        fusion_results = getattr(state, "fusion_results", None) or {}
        gold_count = sum(
            1 for r in fusion_results.values()
            if r.confidence_level in ("Gold", "gold", "GOLD")
        )

        if gold_count >= 1:
            suggestions.append(f"发现 {gold_count} 个 Gold 候选TF，建议优先进行 CRISPR/Cas9 敲除验证")
        silver = sum(1 for r in fusion_results.values()
                     if r.confidence_level in ("Silver", "silver", "SILVER"))
        if silver >= 1:
            suggestions.append(f"发现 {silver} 个 Silver 候选TF，建议进行 Y1H/Dual-LUC 结合实验")
        if not gold_count and not silver:
            suggestions.append("未发现高置信候选，建议补充更多组学数据")

        return suggestions
