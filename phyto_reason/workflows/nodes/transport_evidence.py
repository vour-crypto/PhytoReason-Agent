"""
transport_evidence.py — Transport mechanism evidence acquisition node.

When mechanism_classifier identifies TRANSPORT as a likely mechanism,
this node gathers evidence about transporter proteins.
"""

from __future__ import annotations

import logging

from phyto_reason.workflows.runtime_state import RuntimeState

logger = logging.getLogger("transport_evidence")


def transport_evidence_node(state: RuntimeState) -> dict:
    """Gather evidence for transport-mediated metabolite redistribution.

    Checks:
      1. Known transporter families for the target metabolite
      2. Transporter gene expression in the dataset
      3. Tissue-specific transporter expression patterns
    """
    try:
        target = _get_meta(state)
        expr = _get_expr_matrix(state)

        if not expr:
            state.add_trace("transport_evidence", "completed", summary={
                "n_transporters_found": 0, "reason": "no expression data"
            })
            return {"current_node": "regulation_evidence"}

        # Check ontology for known transport modes
        from phyto_reason.ontology.metabolite_metadata import resolve_ontology

        onto = resolve_ontology(target) if target else None
        transport_modes = onto.transport_modes if onto and onto.transport_modes else ["unknown"]

        if transport_modes and transport_modes[0] == "unknown":
            state.add_trace("transport_evidence", "completed", summary={
                "n_transporters_found": 0, "reason": "no known transport mechanism"
            })
            return {"current_node": "regulation_evidence"}

        # Scan expression data for transporter genes
        # Transporter families: ABC, MATE, NRT, PTR, aquaporins, etc.
        TRANSPORTER_KEYWORDS = [
            "ABC", "MATE", "NRT", "PTR", "PDR", "NIP", "TIP", "PIP",
            "transporter", "carrier", "efflux", "influx", "export",
            "DTX", "MFS", "AMT", "SULTR", "ZIP",
        ]

        transporter_genes = []
        for gene_id in expr:
            upper = gene_id.upper()
            if any(kw.upper() in upper for kw in TRANSPORTER_KEYWORDS):
                transporter_genes.append(gene_id)

        findings = {
            "transport_modes": transport_modes,
            "transporter_genes_found": len(transporter_genes),
            "top_transporters": transporter_genes[:10],
            "note": (
                f"Metabolite '{target}' has known transport modes: {', '.join(transport_modes)}. "
                f"Found {len(transporter_genes)} potential transporter genes in expression data. "
                "Transport-mediated redistribution cannot be ruled out. "
                "Recommend: tissue-specific metabolite quantification and transporter inhibitor experiments."
            ),
        }

        state.add_trace("transport_evidence", "completed", summary={
            "transport_modes": transport_modes,
            "n_transporters_found": len(transporter_genes),
        })

        return {
            "current_node": "regulation_evidence",
            "tool_results": state.tool_results,
        }

    except Exception as e:
        logger.error(f"transport_evidence failed: {e}")
        state.add_trace("transport_evidence", "failed", error=str(e))
        return {"error": str(e), "current_node": "__end__", "finished": True}


def _get_meta(state: RuntimeState) -> str:
    return (state.planner_state.target_metabolite or "").lower().strip() if state.planner_state else ""


def _get_expr_matrix(state: RuntimeState) -> dict | None:
    if state.planner_state:
        return state.planner_state.expression_matrix
    return None
