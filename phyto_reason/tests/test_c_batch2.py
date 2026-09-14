from types import SimpleNamespace

import numpy as np


def _five_tissue_state():
    from phyto_reason.models.planner_state import PlannerState
    from phyto_reason.workflows.runtime_state import RuntimeState

    tissues = ["Flower", "Leaf", "Root", "Stem", "Seed"]
    metadata = {"samples": [
        {"sample_id": f"{t}-{i}", "condition": "", "tissue": t}
        for t in tissues for i in range(3)
    ]}
    matrix = {
        "quercetin": {
            f"{t}-{i}": float(10 + 20 * tissues.index(t) + i * 0.01)
            for t in tissues for i in range(3)
        },
        "stable_metabolite": {
            f"{t}-{i}": 2.0 + i * 0.01 for t in tissues for i in range(3)
        },
    }
    return RuntimeState(planner_state=PlannerState(
        metabolite_matrix=matrix, sample_metadata=metadata,
    ))


def test_dam_pairwise_runs_all_ten_tissue_pairs_and_reports_feasibility():
    from phyto_reason.workflows.execution_nodes import dam_analysis_node

    state = _five_tissue_state()
    report = dam_analysis_node(state)["dam_report"]
    assert report["n_total"] == 15
    assert len(report["group_sizes"]) == 5
    assert report["feasibility"] == "standard"
    assert report["pairwise_method"].startswith("moderated_t")
    assert report["n_pairwise_tests"] == 10
    assert {r["group_a"] for r in report["pairwise_results"]}
    assert all({"metabolite", "group_a", "group_b", "log2fc", "p_value", "q_value", "direction"} <= set(r)
               for r in report["pairwise_results"])


def test_joint_enrichment_uses_fisher_combined_deg_dam(monkeypatch):
    from phyto_reason.workflows.nodes import joint_enrichment as je
    from phyto_reason.models.planner_state import PlannerState
    from phyto_reason.workflows.runtime_state import RuntimeState

    monkeypatch.setattr(je, "_get_pathway_data", lambda key: {
        "map_id": "map00941", "name": "Flavonoid biosynthesis", "genes": ["GENE1", "GENE2"]
    })
    state = RuntimeState(planner_state=PlannerState(
        expression_matrix={f"GENE{i}": {f"S{j}": 1.0 for j in range(6)} for i in (1, 2, 3)},
        target_metabolite="flavonoid",
    ))
    state.deg_report = {"top_genes": [
        {"gene_id": "GENE1", "effect_size": 0.8, "is_fdr_sig": True, "p_value": 0.01},
        {"gene_id": "GENE2", "effect_size": 0.8, "is_fdr_sig": True, "p_value": 0.02},
    ], "n_genes_tested": 3}
    state.dam_report = {"top_metabolites": [
        {"metabolite": "flavonoid", "effect_size": 0.8, "is_fdr_sig": True, "p_value": 0.01}
    ], "n_metabolites_tested": 1}
    report = je.joint_enrichment_node(state)["joint_enrichment_report"]
    row = report["results"][0]
    assert row["p_dam"] < 1
    assert row["p_combined"] if "p_combined" in row else row["p_value"]
    assert "q_combined" in row
    assert "Fisher" in row["combined_test"]
    assert "Fisher_combined_DEG_DAM" in report["method"]


def test_joint_enrichment_preserves_tiny_dam_pvalue(monkeypatch):
    from phyto_reason.workflows.nodes import joint_enrichment as je
    from phyto_reason.models.planner_state import PlannerState
    from phyto_reason.workflows.runtime_state import RuntimeState

    monkeypatch.setattr(je, "_get_pathway_data", lambda key: {
        "map_id": "map00941", "name": "Flavonoid biosynthesis", "genes": ["GENE1"]
    })
    state = RuntimeState(planner_state=PlannerState(
        expression_matrix={"GENE1": {f"S{j}": 1.0 for j in range(6)}},
        target_metabolite="flavonoid",
    ))
    state.deg_report = {"top_genes": [
        {"gene_id": "GENE1", "effect_size": 0.8, "is_fdr_sig": True, "p_value": 0.01},
    ], "n_genes_tested": 1}
    state.dam_report = {"top_metabolites": [
        {"metabolite": "flavonoid", "effect_size": 0.8, "is_fdr_sig": True,
         "p_value": 0.0, "p_value_raw": 1e-12}
    ], "n_metabolites_tested": 1}
    report = je.joint_enrichment_node(state)["joint_enrichment_report"]
    row = report["results"][0]
    assert row["p_dam_raw"] < 1e-6
    assert row["p_combined_raw"] < 1e-6


def test_ms2_handler_calls_remote_chain_when_query_is_given(monkeypatch):
    from phyto_reason.agent.tool_handlers import handle_annotate_ms2_spectrum

    fake = SimpleNamespace(
        status="partial",
        warnings=["MassBank unavailable; downgraded to MoNA"],
        metadata={"source": "MoNA", "count": 1},
        evidence_list=[SimpleNamespace(metadata={"actual_source": "MoNA"})],
    )
    monkeypatch.setattr("phyto_reason.tools.spectra_tools.remote_annotation_fallback", lambda *a, **k: fake)
    result = handle_annotate_ms2_spectrum({
        "precursor_mz": 336.1230,
        "fragments": [321.0996, 306.0761, 278.0812],
        "species": "黄连",
        "target_metabolite": "berberine",
    })
    assert "Remote spectral comparison" in result
    assert "actual=MoNA" in result
    assert "MassBank unavailable" in result


def test_pubchem_batch_preserves_source_and_missing_warning(monkeypatch):
    from phyto_reason.tools.spectra_tools import PubChemCompoundPropertiesTool

    class Response:
        def __init__(self, payload): self.payload = payload
        def raise_for_status(self): pass
        def json(self): return self.payload

    def fake_get(url, **kwargs):
        if "/synonyms/" in url:
            return Response({"InformationList": {"Information": [{"Synonym": ["123-45-6"]}]}})
        if "missing" in url:
            return Response({"PropertyTable": {"Properties": []}})
        return Response({"PropertyTable": {"Properties": [{"CID": 1, "InChIKey": "TEST"}]}})

    monkeypatch.setattr("phyto_reason.tools.spectra_tools.requests.get", fake_get)
    result = PubChemCompoundPropertiesTool().run(
        names=["known_batch_name", "missing_batch_name"], properties=["InChIKey", "CAS"]
    )
    assert result.status == "partial"
    assert result.metadata["count"] == 2
    assert result.metadata["resolved"] == 1
    assert result.metadata["records"][0]["actual_source"] == "PubChem"
    assert result.metadata["records"][0]["properties"]["CAS"] == "123-45-6"
    assert result.metadata["records"][1]["properties"] == {}


def test_metabolite_parser_records_gb18030_detection(tmp_path):
    from phyto_reason.ingestion.parsers.metabolite_parser import parse_metabolite

    path = tmp_path / "gbk.csv"
    path.write_bytes("代谢物,S-1,S-2\n黄酮,1,2\n".encode("gb18030"))
    parsed = parse_metabolite(path)
    assert parsed.provenance["encoding"] in {"gb18030", "gbk", "cp936"}
    assert parsed.matrix["黄酮"]["S-1"] == 1.0
    assert any("encoding" in warning.lower() or "detected" in warning.lower() for warning in parsed.warnings)


def test_utf8_cjk_metabolite_name_stays_utf8_and_is_not_garbled(tmp_path):
    from phyto_reason.ingestion.parsers.metabolite_parser import parse_metabolite

    path = tmp_path / "utf8.csv"
    path.write_text("代谢物,F-1,F-2\n黄酮,1,2\n", encoding="utf-8")
    parsed = parse_metabolite(path)
    assert parsed.provenance["encoding"] == "utf-8"
    assert "黄酮" in parsed.matrix
    assert "鍠" not in "".join(parsed.matrix)


def test_multiomics_node_returns_report_for_five_common_samples():
    from phyto_reason.models.planner_state import PlannerState
    from phyto_reason.workflows.runtime_state import RuntimeState
    from phyto_reason.workflows.execution_nodes import multiomics_integration_node

    samples = [f"S{i}" for i in range(6)]
    expr = {"G1": {s: float(i + 1) for i, s in enumerate(samples)},
            "G2": {s: float((i + 1) * 2) for i, s in enumerate(samples)}}
    meta = {"M1": {s: float(i + 1) for i, s in enumerate(samples)},
            "M2": {s: float(6 - i) for i, s in enumerate(samples)}}
    state = RuntimeState(planner_state=PlannerState(expression_matrix=expr, metabolite_matrix=meta))
    result = multiomics_integration_node(state)
    report = result["multiomics_report"]
    assert "error" not in report
    assert report["n_correlation_pairs_tested"] > 0
    assert "summary" not in report or isinstance(report["summary"], str)


def test_metabolite_parser_excludes_aggregated_and_bare_numeric_metadata_columns(tmp_path):
    from phyto_reason.ingestion.parsers.metabolite_parser import parse_metabolite

    path = tmp_path / "aggregated.csv"
    path.write_text(
        "Name,F-1,F-2,Average-F,mass,area\nquercetin,1,2,1.5,302.1,88\n",
        encoding="utf-8",
    )
    parsed = parse_metabolite(path)
    assert parsed.sample_ids == ["F-1", "F-2"]
    assert "Average-F" in parsed.classification["aggregated_columns"]
    assert any("aggregated column excluded" in warning for warning in parsed.warnings)
