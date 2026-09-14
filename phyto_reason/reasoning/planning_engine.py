"""
planning_engine.py — 研究规划引擎。
"""

from __future__ import annotations

from phyto_reason.models.planner_state import PlannerState

CORRELATION_MIN_SAMPLES = 4
WGCNA_MIN_SAMPLES = 15


class PlanningEngine:
    """研究规划引擎。"""

    def __init__(self) -> None:
        pass

    def create_plan(
        self,
        research_question: str = "",
        species: str = "",
        target_metabolite: str = "",
        target_pathway: str = "",
        sample_count: int = 0,
        has_expression: bool = False,
        has_metabolite: bool = False,
        has_annotation: bool = False,
        has_promoter: bool = False,
        has_itak: bool = False,
    ) -> PlannerState:
        state = PlannerState(
            research_question=research_question,
            species=species,
            target_metabolite=target_metabolite,
            target_pathway=target_pathway,
            sample_count=sample_count,
            has_expression=has_expression,
            has_metabolite=has_metabolite,
            has_annotation=has_annotation,
            has_promoter=has_promoter,
            has_itak=has_itak,
        )
        state.add_decision(
            stage="initialization",
            decision=f"开始研究规划: 物种={species}, 目标={target_metabolite}",
            rationale=f"样本量={sample_count}, 表达矩阵={'有' if has_expression else '无'}, "
                      f"代谢物矩阵={'有' if has_metabolite else '无'}",
        )
        total_checks = 5
        missing = sum([not has_expression, not has_metabolite, not has_annotation, not has_promoter, not has_itak])
        state.missing_data_ratio = missing / total_checks
        if state.missing_data_ratio > 0.6:
            state.add_warning(f"数据缺失严重 ({missing}/{total_checks} 项缺失)，分析结果置信度将受限。")
        # Inline workflow plan (replaces deprecated WorkflowSelector)
        tools = []
        skipped = []
        if has_expression:
            tools.append("correlation_analysis")
            if sample_count >= WGCNA_MIN_SAMPLES:
                tools.append("wgcna")
            else:
                skipped.append("wgcna (样本不足)")
        else:
            skipped.extend(["correlation_analysis", "wgcna", "deg_analysis"])
        if has_promoter:
            tools.append("motif_scan")
        else:
            skipped.append("motif_scan")
        if has_annotation:
            tools.append("tf_annotation")
        else:
            skipped.append("tf_annotation")
        tools.append("literature_search")
        tools.append("kegg_pathway")
        state.selected_tools = tools
        state.skipped_analyses = skipped
        state.add_decision(
            stage="tool_routing",
            decision=f"选定分析: {', '.join(tools[:4])}",
            rationale=f"基于可用数据 ({'expression' if has_expression else 'no expression'}) 自动选择工具",
            alternatives=[f"跳过: {s}" for s in skipped],
        )
        self._assess_feasibility(state)
        return state

    @staticmethod
    def _assess_feasibility(state: PlannerState) -> None:
        if state.sample_count < 3:
            state.overall_feasibility = "infeasible"
            state.confidence_assessment = "数据严重不足，任何分析结果均不可靠。"
        elif state.sample_count < CORRELATION_MIN_SAMPLES:
            state.overall_feasibility = "limited"
            state.confidence_assessment = "样本量有限，仅适合探索性分析。"
        elif not state.has_expression or not state.has_metabolite:
            state.overall_feasibility = "partial"
            state.confidence_assessment = "缺少部分组学数据，分析维度受限。"
        elif state.sample_count >= WGCNA_MIN_SAMPLES and state.has_expression and state.has_metabolite:
            state.overall_feasibility = "feasible"
            state.confidence_assessment = "数据充分，适合完整分析。"
        else:
            state.overall_feasibility = "feasible"
            state.confidence_assessment = "基本可行，部分分析受到数据限制。"
        state.add_decision(
            stage="feasibility_assessment",
            decision=f"研究可行性: {state.overall_feasibility}",
            rationale=state.confidence_assessment,
        )
