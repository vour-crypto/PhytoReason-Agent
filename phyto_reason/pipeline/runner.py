"""
runner.py — Clean wrapper around WorkflowRunner.

Exposes run_scientific_pipeline as a pure function:
    Input: species, target_metabolite, data matrices
    Output: structured PipelineResult

This is what the LLM orchestration layer calls.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Any

from phyto_reason.workflows.workflow_runner import WorkflowRunner
from phyto_reason.models.mechanistic_hypothesis import MechanisticHypothesis

logger = logging.getLogger("pipeline_runner")


@dataclass
class PipelineResult:
    """Structured result from the scientific pipeline."""
    success: bool = False
    error: str | None = None

    # Summary
    species: str = ""
    target_metabolite: str = ""
    target_pathway: str = ""
    n_samples: int = 0

    # Quality
    quality_passed: bool = False
    quality_issues: list[str] = field(default_factory=list)

    # Analysis
    mechanism_type: str = "unknown"
    mechanism_probabilities: dict = field(default_factory=dict)
    mechanism_uncertainty: str = "high"

    # Results
    competing_hypotheses: list[MechanisticHypothesis] = field(default_factory=list)
    hypothesis_text: str = ""
    execution_trace: list[dict] = field(default_factory=list)

    # Context
    ingestion_warnings: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        return {
            "success": self.success,
            "species": self.species,
            "target_metabolite": self.target_metabolite,
            "quality_passed": self.quality_passed,
            "mechanism_type": self.mechanism_type,
            "n_hypotheses": len(self.competing_hypotheses),
            "top_hypothesis": self.competing_hypotheses[0].title if self.competing_hypotheses else "",
        }


def run_scientific_pipeline(
    species: str = "",
    target_metabolite: str = "",
    target_pathway: str = "",
    expression_matrix: dict | None = None,
    metabolite_matrix: dict | None = None,
    promoter_sequences: dict | None = None,
    verbose: bool = False,
) -> PipelineResult:
    """Run the complete scientific analysis pipeline.

    Pure function — no side effects, no internal state.

    Args:
        species: Species name (e.g., 'Coptis chinensis')
        target_metabolite: Target metabolite (e.g., 'berberine')
        target_pathway: Optional target pathway
        expression_matrix: {gene_id: {sample_id: value}}
        metabolite_matrix: {metabolite_id: {sample_id: value}}
        promoter_sequences: {gene_id: sequence_string}
        verbose: Detailed logging

    Returns:
        PipelineResult with ranked hypotheses and evidence
    """
    result = PipelineResult(
        species=species,
        target_metabolite=target_metabolite,
        target_pathway=target_pathway,
    )

    if not expression_matrix:
        result.success = False
        result.error = "No expression data provided."
        result.quality_issues = ["Missing expression matrix"]
        return result

    # Determine sample count from data
    sample_count = 0
    first_key = next(iter(expression_matrix), None)
    if first_key:
        sample_count = len(expression_matrix[first_key])
    result.n_samples = sample_count

    try:
        runner = WorkflowRunner()
        state = runner.run(
            species=species,
            target_metabolite=target_metabolite,
            target_pathway=target_pathway,
            sample_count=sample_count,
            has_expression=bool(expression_matrix),
            has_metabolite=bool(metabolite_matrix),
            has_promoter=bool(promoter_sequences),
            expression_matrix=expression_matrix,
            metabolite_matrix=metabolite_matrix,
            promoter_sequences=promoter_sequences,
            verbose=verbose,
        )

        # Populate result from state
        result.success = not state.error
        result.error = state.error
        result.quality_passed = state.data_quality_passed
        result.quality_issues = state.data_quality_issues

        result.mechanism_type = getattr(state, "mechanism_type", "") or "unknown"
        result.mechanism_probabilities = getattr(state, "mechanism_probabilities", {}) or {}
        result.mechanism_uncertainty = getattr(state, "mechanism_uncertainty", "high") or "high"

        result.competing_hypotheses = getattr(state, "competing_hypotheses", []) or []
        result.hypothesis_text = getattr(state, "mechanism_hypothesis", "") or state.hypothesis

        result.execution_trace = [
            {"node": step.node, "status": step.status, "error": step.error}
            for step in (state.execution_trace or [])
        ]

        result.ingestion_warnings = getattr(state, "ingestion_warnings", []) or []

        logger.info(
            f"Pipeline completed: species={species}, metabolite={target_metabolite}, "
            f"mechanism={result.mechanism_type}, hypotheses={len(result.competing_hypotheses)}"
        )

        return result

    except Exception as e:
        logger.error(f"Pipeline failed: {e}", exc_info=True)
        result.success = False
        result.error = str(e)
        return result


# Convenience: run from file paths
def run_from_files(
    expression_file: str = "",
    metabolite_file: str = "",
    species: str = "",
    target_metabolite: str = "",
    target_pathway: str = "",
) -> PipelineResult:
    """Run pipeline by loading data from CSV/TSV files."""
    from phyto_reason.workflows.workflow_runner import WorkflowRunner

    runner = WorkflowRunner()
    state = runner.run_from_files(
        expression_file=expression_file or None,
        metabolite_file=metabolite_file or None,
        species=species,
        target_metabolite=target_metabolite,
        target_pathway=target_pathway,
    )

    result = PipelineResult(
        success=not state.error,
        error=state.error,
        species=species,
        target_metabolite=target_metabolite,
        target_pathway=target_pathway,
        quality_passed=state.data_quality_passed,
        quality_issues=state.data_quality_issues,
        mechanism_type=getattr(state, "mechanism_type", "unknown"),
        competing_hypotheses=getattr(state, "competing_hypotheses", []),
        hypothesis_text=getattr(state, "mechanism_hypothesis", "") or state.hypothesis,
        ingestion_warnings=getattr(state, "ingestion_warnings", []),
    )
    return result
