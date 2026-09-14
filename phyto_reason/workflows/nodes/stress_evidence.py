"""
stress_evidence.py — Stress mechanism evidence acquisition node.

When mechanism_classifier identifies STRESS as a likely mechanism,
this node checks stress markers and ontology associations.
"""

from __future__ import annotations

import logging
import math

import numpy as np

from phyto_reason.workflows.runtime_state import RuntimeState

logger = logging.getLogger("stress_evidence")


# Known stress marker genes (conserved across plants)
STRESS_MARKERS = [
    "HSF", "HSP", "DREB", "ERF", "WRKY", "NAC", "LEA", "PR1", "PR2",
    "PAL", "CHS", "LOX", "AOS", "CAT", "SOD", "APX", "GPX", "GST",
    "P5CS", "NCED", "ZEP", "ABA", "RD29", "COR", "KIN",
]


def stress_evidence_node(state: RuntimeState) -> dict:
    """Gather evidence for stress-induced metabolic redistribution.

    Checks:
      1. Stress marker gene expression levels
      2. Ontology stress associations for target metabolite
      3. Coordinated stress response signals in the dataset
    """
    try:
        target = _get_meta(state)
        expr = _get_expr_matrix(state)

        if not expr:
            state.add_trace("stress_evidence", "completed", summary={
                "stress_markers_activated": 0, "reason": "no expression data"
            })
            return {"current_node": "regulation_evidence"}

        # 1. Check stress marker expression
        marker_expressions = {}
        for gene_id in expr:
            upper = gene_id.upper()
            for marker in STRESS_MARKERS:
                if marker.upper() in upper:
                    values = [v for v in expr[gene_id].values()
                             if isinstance(v, (int, float)) and not math.isnan(v)]
                    if values:
                        marker_expressions[gene_id] = {
                            "mean": float(np.mean(values)),
                            "cv": float(np.std(values) / (abs(np.mean(values)) + 1e-10)),
                            "marker_family": marker,
                        }
                    break

        # 2. Ontology stress associations
        from phyto_reason.ontology.metabolite_metadata import resolve_ontology

        onto = resolve_ontology(target) if target else None
        stress_associations = onto.stress_associations if onto else []

        # 3. Assess stress signal strength
        high_cv_markers = [
            gid for gid, info in marker_expressions.items()
            if info["cv"] > 0.5
        ]
        stress_signal_strong = len(high_cv_markers) >= 3

        # Build findings
        if stress_signal_strong:
            interpretation = (
                f"Multiple stress markers ({len(high_cv_markers)}) show high variability (CV>0.5). "
                "This suggests stress-induced transcriptional reprogramming may be occurring. "
                "The observed metabolite changes could be part of a general stress response "
                "rather than specific transcriptional regulation of the target pathway."
            )
        elif marker_expressions:
            interpretation = (
                f"{len(marker_expressions)} stress markers detected, but variability is low. "
                "Stress response is not strongly activated in this dataset."
            )
        else:
            interpretation = (
                "No canonical stress markers found in the expression data. "
                "However, this does not rule out stress effects — "
                "stress markers may be species-specific or not annotated."
            )

        if stress_associations:
            interpretation += f" Target metabolite '{target}' is known to be associated with: {', '.join(stress_associations[:5])}."

        state.add_trace("stress_evidence", "completed", summary={
            "n_stress_markers": len(marker_expressions),
            "high_cv_markers": len(high_cv_markers),
            "stress_associations": stress_associations[:3],
            "stress_signal": "strong" if stress_signal_strong else "weak",
        })

        return {
            "current_node": "regulation_evidence",
            "tool_results": state.tool_results,
        }

    except Exception as e:
        logger.error(f"stress_evidence failed: {e}")
        state.add_trace("stress_evidence", "failed", error=str(e))
        return {"error": str(e), "current_node": "__end__", "finished": True}


def _get_meta(state: RuntimeState) -> str:
    return (state.planner_state.target_metabolite or "").lower().strip() if state.planner_state else ""


def _get_expr_matrix(state: RuntimeState) -> dict | None:
    if state.planner_state:
        return state.planner_state.expression_matrix
    return None
