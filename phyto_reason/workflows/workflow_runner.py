"""
workflow_runner.py — 工作流执行入口（v3.0 + ingestion）。
"""

from __future__ import annotations

import logging
from datetime import datetime
from pathlib import Path

from phyto_reason.workflows.runtime_state import RuntimeState
from phyto_reason.workflows.workflow_graph import compile_workflow
from phyto_reason.models.planner_state import PlannerState
from phyto_reason.models.workflow_plan import WorkflowPlan

logger = logging.getLogger("workflow_runner")


class WorkflowRunner:
    """科研工作流执行器。"""

    def __init__(self) -> None:
        self.app = compile_workflow()

    def run(
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
        data_path: str = "",
        verbose: bool = False,
        expression_matrix: dict | None = None,
        metabolite_matrix: dict | None = None,
        promoter_sequences: dict | None = None,
        sample_metadata: dict | None = None,
        annotation_index=None,
        workflow_plan: WorkflowPlan | None = None,
    ) -> RuntimeState:
        """执行完整科研工作流。

        Args:
            研究参数 (同 PIAgent.plan)
            data_path: 包含表达/代谢物/FASTA 文件的目录路径
            expression_matrix: 预加载的表达矩阵 {gene: {sample: value}}
            metabolite_matrix: 预加载的代谢物矩阵 {metabolite: {sample: value}}
            promoter_sequences: 预加载的启动子序列 {gene: sequence}

        Returns:
            最终 RuntimeState
        """
        loaded_expr = expression_matrix
        loaded_meta = metabolite_matrix
        loaded_prom = promoter_sequences
        if loaded_expr:
            has_expression = True
            first_key = next(iter(loaded_expr), None)
            if first_key:
                sample_count = max(sample_count, len(loaded_expr[first_key]))
        if loaded_meta:
            has_metabolite = True
        if loaded_prom:
            has_promoter = True

        ps = PlannerState(
            research_question=research_question,
            species=species, target_metabolite=target_metabolite,
            target_pathway=target_pathway, sample_count=sample_count,
            has_expression=has_expression, has_metabolite=has_metabolite,
            has_annotation=has_annotation, has_promoter=has_promoter,
            has_itak=has_itak, data_path=data_path,
            expression_matrix=loaded_expr, metabolite_matrix=loaded_meta,
            promoter_sequences=loaded_prom,
            sample_metadata=sample_metadata,
        )
        initial_state = RuntimeState(
            planner_state=ps, current_node="data_quality_gate",
            annotation_index=annotation_index,
        )
        if workflow_plan is not None:
            initial_state.workflow_plan = workflow_plan

        logger.info(f"Starting workflow: {research_question[:80]}")
        logger.info(f"  species={species}, target={target_metabolite}, samples={sample_count}")
        if data_path:
            logger.info(f"  data_path={data_path}")

        config = {"recursion_limit": 50}
        if verbose:
            logger.info("  verbose mode enabled")

        try:
            # Emit start progress
            initial_state.emit_progress("workflow", "running", f"Starting analysis for {target_metabolite} in {species}")

            result = self.app.invoke(initial_state, config)
            final_state = result if isinstance(result, RuntimeState) else RuntimeState(**result)
            _merge_inplace_state(initial_state, final_state)

            # Build progress from final state (LangGraph uses state copies, so trace
            # from initial_state may not propagate; we reconstruct from final state)
            final_state.progress_messages = _build_progress(final_state)

            logger.info(f"Workflow completed: {final_state.to_summary()}")
            return final_state
        except Exception as e:
            logger.error(f"Workflow failed: {e}")
            initial_state.emit_progress("workflow", "failed", str(e)[:100])
            initial_state.error = str(e)
            initial_state.finished = True
            return initial_state


    def run_from_files(
        self,
        expression_file: str | Path | None = None,
        metabolite_file: str | Path | None = None,
        metadata_file: str | Path | None = None,
        species: str = "",
        target_metabolite: str = "",
        target_pathway: str = "",
        verbose: bool = False,
    ) -> RuntimeState:
        """Parse uploaded matrices, align samples, then execute the workflow."""
        from phyto_reason.ingestion.parsers.expression_parser import parse_expression
        from phyto_reason.ingestion.parsers.metabolite_parser import parse_metabolite
        from phyto_reason.ingestion.parsers.metadata_parser import parse_metadata
        from phyto_reason.ingestion.validators.sample_alignment_validator import validate_alignment
        from phyto_reason.ingestion.validators.matrix_integrity_validator import validate_matrix_integrity
        from phyto_reason.ingestion.validators.missing_value_validator import validate_missing_values
        from phyto_reason.ingestion.mappers.species_mapper import map_species
        from phyto_reason.ingestion.mappers.annotation_mapper import load_annotation
        from phyto_reason.workflows.execution_nodes import _infer_sample_groups

        warnings: list[str] = []
        errors: list[str] = []
        parsed_expr = parsed_meta = parsed_md = None

        if expression_file:
            try:
                parsed_expr = parse_expression(expression_file)
                warnings.extend(parsed_expr.warnings)
            except Exception as exc:
                errors.append(f"Expression parse failed: {exc}")
        if metabolite_file:
            try:
                parsed_meta = parse_metabolite(metabolite_file)
                warnings.extend(parsed_meta.warnings)
            except Exception as exc:
                errors.append(f"Metabolite parse failed: {exc}")
        if metadata_file:
            try:
                parsed_md = parse_metadata(metadata_file)
                warnings.extend(parsed_md.warnings)
            except Exception as exc:
                errors.append(f"Metadata parse failed: {exc}")

        annotation_index = None
        annotation_candidates: list[Path] = []
        if expression_file:
            annotation_candidates.extend(Path(expression_file).parent.glob("*annotation*"))
            annotation_candidates.extend(Path(expression_file).parent.glob("*Annotation*"))
        for annotation_path in annotation_candidates[:1]:
            try:
                annotation_index = load_annotation(annotation_path)
                warnings.extend(annotation_index.warnings)
            except Exception as exc:
                warnings.append(f"Annotation load failed: {exc}")

        expr_dict = parsed_expr.to_workflow_dict() if parsed_expr else None
        meta_dict = parsed_meta.to_workflow_dict() if parsed_meta else None
        if expr_dict is not None:
            expr_qc = validate_matrix_integrity(expr_dict, name="expression")
            if not expr_qc.all_passed and expr_qc.sample_count_warning:
                warnings.append(f"Expression matrix: {expr_qc.malformed_warning}")
            missing = validate_missing_values(expr_dict)
            if not missing.passed:
                warnings.extend(missing.warnings[:3])
        if meta_dict is not None:
            meta_qc = validate_matrix_integrity(meta_dict, name="metabolite")
            if not meta_qc.all_passed and meta_qc.sample_count_warning:
                warnings.append(f"Metabolite matrix: {meta_qc.malformed_warning}")

        alignment = validate_alignment(
            expression_samples=parsed_expr.sample_ids if parsed_expr else None,
            metabolite_samples=parsed_meta.sample_ids if parsed_meta else None,
            metadata_samples=parsed_md.sample_ids if parsed_md else None,
        )
        if alignment.n_common < 3:
            warnings.append(f"Only {alignment.n_common} common samples after alignment (< 3)")

        metadata_dict = None
        if parsed_md:
            metadata_dict = {
                "source_file": parsed_md.source_file,
                "sample_ids": parsed_md.sample_ids,
                "columns_present": parsed_md.columns_present,
                "columns_missing": parsed_md.columns_missing,
                "samples": [sample.model_dump() for sample in parsed_md.samples],
            }
        sample_alignment_report = {
            "n_common": alignment.n_common,
            "common_samples": alignment.common_samples,
            "dropped_from_expression": alignment.dropped_from_expression,
            "dropped_from_metabolite": alignment.dropped_from_metabolite,
            "dropped_from_metadata": alignment.dropped_from_metadata,
            "groups": _infer_sample_groups(expr_dict or meta_dict or {}, metadata_dict or {}),
            "group_source": "metadata" if metadata_dict else "sample_name_fallback",
        }
        species_result = map_species(species) if species else None
        canonical_species = species_result.canonical if species_result else species
        if species_result and species_result.warnings:
            warnings.extend(species_result.warnings)
        sample_count = alignment.n_common or (parsed_expr.n_samples if parsed_expr else 0)
        state = self.run(
            research_question="",
            species=canonical_species,
            target_metabolite=target_metabolite,
            target_pathway=target_pathway,
            sample_count=sample_count,
            has_expression=parsed_expr is not None,
            has_metabolite=parsed_meta is not None,
            expression_matrix=expr_dict,
            metabolite_matrix=meta_dict,
            sample_metadata=metadata_dict,
            annotation_index=annotation_index,
            verbose=verbose,
        )
        state.ingestion_warnings = warnings
        state.parsing_errors = errors
        state.sample_alignment_report = sample_alignment_report
        if parsed_expr:
            state.data_quality_report.setdefault("expression_provenance", parsed_expr.provenance)
        if parsed_meta:
            state.data_quality_report.setdefault("metabolite_provenance", parsed_meta.provenance)
        if errors or alignment.n_common < 3:
            state.error = "; ".join(errors + warnings[:2])
            state.finished = True
        return state


def _merge_inplace_state(initial_state: RuntimeState, final_state: RuntimeState) -> None:
    """LangGraph 只合并节点返回的 dict。

    节点内原地写入的字段（execution_trace / completed_nodes / failures /
    未随返回 dict 传播的中间结果等）会丢失——从 initial_state 补回 final。
    仅当 initial 值非空时覆盖（节点在 initial 对象上原地变更，是最新视图）。
    """
    empty = (None, [], {}, "", 0)
    for field_name in initial_state.model_fields:
        val = getattr(initial_state, field_name, None)
        if val not in empty:
            setattr(final_state, field_name, val)


def _build_progress(state: RuntimeState) -> list[str]:
    """Build a progress log from the pipeline state after execution."""
    msgs = []

    def _add(node: str, status: str, detail: str = ""):
        ts = datetime.now().strftime("%H:%M:%S")
        emoji = {"pass": "[OK]", "fail": "[FAIL]", "warn": "[WARN]", "skip": "[SKIP]"}.get(status, "  ")
        msg = f"{ts} {emoji} {node}"
        if detail:
            msg += f": {detail}"
        msgs.append(msg)

    # Phase 0: DEG analysis
    deg_report = getattr(state, "deg_report", None) or {}
    if deg_report.get("n_genes_tested", 0) > 0:
        n_sig = deg_report.get("n_fdr_significant", 0)
        n_eff = deg_report.get("n_effect_ranked", 0)
        _add("deg_analysis", "pass",
             f"n_tested={deg_report['n_genes_tested']}, FDR-sig={n_sig}, effect-ranked={n_eff}")
    elif deg_report:
        _add("deg_analysis", "skip", deg_report.get("reason", "no expression data"))

    # Phase 0b: DAM analysis
    dam_report = getattr(state, "dam_report", None) or {}
    if dam_report.get("n_metabolites_tested", 0) > 0:
        _add("dam_analysis", "pass",
             f"n_tested={dam_report['n_metabolites_tested']}, FDR-sig={dam_report.get('n_fdr_significant', 0)}")

    # Phase 0c: Multi-omics integration
    multiomics = getattr(state, "multiomics_report", None) or {}
    if multiomics.get("n_correlation_pairs_tested", 0) > 0:
        _add("multiomics_integration", "pass",
             f"pairs={multiomics['n_correlation_pairs_tested']}, "
             f"sig={multiomics.get('n_significant_pairs', 0)}, "
             f"modules={len(multiomics.get('modules', {}))}")

    # Phase 0d: Joint KEGG enrichment
    enrichment = getattr(state, "joint_enrichment_report", None) or {}
    if enrichment.get("n_pathways_tested", 0) > 0:
        _add("joint_enrichment", "pass",
             f"pathways={enrichment['n_pathways_tested']}, "
             f"sig={enrichment.get('n_significant_pathways', 0)}")

    # Phase 0e: Quadrant plot
    quadrant = getattr(state, "quadrant_plot_report", None) or {}
    if quadrant.get("n_pairs", 0) > 0:
        _add("quadrant_plot", "pass",
             f"pairs={quadrant['n_pairs']}, sync={quadrant.get('synchronicity_ratio', 0):.2f}")

    # Phase 0f: Correlation network
    corr_net = getattr(state, "correlation_network_report", None) or {}
    if corr_net.get("n_edges", 0) > 0:
        _add("correlation_network", "pass",
             f"edges={corr_net['n_edges']}, hubs={len(corr_net.get('hubs', []))}, "
             f"modules={len(corr_net.get('modules', []))}")

    # Phase 0g: WGCNA
    wgcna = getattr(state, "wgcna_report", None) or {}
    if wgcna.get("n_modules", 0) > 0:
        _add("wgcna", "pass",
             f"modules={wgcna['n_modules']}, power={wgcna.get('soft_power', 0)}, "
             f"genes={wgcna.get('n_genes_analyzed', 0)}")

    # Phase 0h: O2PLS
    o2pls = getattr(state, "o2pls_report", None) or {}
    if o2pls.get("n_joint_components", 0) > 0:
        _add("o2pls", "pass",
             f"joint_comp={o2pls['n_joint_components']}, "
             f"X_var={o2pls.get('x_joint_variance_explained', 0):.2f}, "
             f"Y_var={o2pls.get('y_joint_variance_explained', 0):.2f}")

    # Phase 1: Data quality
    if state.data_quality_passed:
        _add("data_quality_gate", "pass", f"n_genes={state.data_quality_report.get('n_genes', '?')}, "
             f"n_samples={state.data_quality_report.get('n_samples', '?')}")
    else:
        _add("data_quality_gate", "fail", f"{len(state.data_quality_issues)} issues")

    # Phase 2: Metabolite validation
    if state.metabolite_validated:
        _add("metabolite_validation", "pass")
    elif state.metabolite_validation_note:
        _add("metabolite_validation", "warn", state.metabolite_validation_note[:80])
    else:
        _add("metabolite_validation", "skip", "no metabolite data")

    # Phase 3: Pathway coherence
    if state.pathway_coherence_passed:
        _add("pathway_coherence", "pass", f"pathway={state.pathway_name}, score={state.pathway_coherence_score:.2f}")
    elif state.pathway_name:
        _add("pathway_coherence", "warn", f"score={state.pathway_coherence_score:.2f}")
    else:
        _add("pathway_coherence", "skip", "no pathway identified")

    # Phase 4: Mechanism classification
    mt = state.mechanism_type or "unknown"
    _add("mechanism_classifier", "pass", f"type={mt}")

    # Phase 5: TF narrowing
    n_tf = len(state.tf_candidates)
    if n_tf > 0:
        _add("tf_narrowing", "pass", f"{n_tf} TF candidates identified")
    else:
        _add("tf_narrowing", "skip", "no TF candidates (insufficient data)")

    # Phase 6: Regulation evidence
    n_fusion = len(state.fusion_results)
    if n_fusion > 0:
        _add("regulation_evidence", "pass", f"{n_fusion} candidates fused with evidence")
    else:
        _add("regulation_evidence", "skip", "no fusion results")

    # Phase 7: Contradiction + Falsification
    n_excluded = len(state.excluded_tfs)
    if n_excluded > 0:
        _add("contradiction_check", "pass", f"{n_excluded} TFs excluded by contradiction")
    else:
        _add("contradiction_check", "pass", "no contradictions found")

    # Phase 8: Hypotheses
    n_hyp = len(state.competing_hypotheses)
    if n_hyp > 0:
        _add("hypothesis_synthesis", "pass", f"{n_hyp} competing hypotheses generated")
    else:
        _add("hypothesis_synthesis", "warn", "no hypotheses — data may be insufficient")

    # Degradation
    for tool in state.degraded_tools:
        _add(f"degraded:{tool}", "warn", "evidence reduced")

    # Error
    if state.error:
        _add("pipeline", "fail", state.error[:100])

    return msgs
