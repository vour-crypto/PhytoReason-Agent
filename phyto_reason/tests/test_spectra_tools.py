import json
from pathlib import Path

from phyto_reason.tools.spectra_tools import (
    MassBankSpectrumSearchTool,
    MonaSpectrumSearchTool,
    PubChemCompoundPropertiesTool,
    parse_massbank_records,
)


def test_massbank_sample_prefix_is_parsed_truthfully():
    path = Path(__file__).parents[2] / "massbank_sample.json"
    raw = path.read_text(encoding="utf-8")
    records = parse_massbank_records(raw)
    assert records
    first = records[0]
    assert first["accession"]
    assert first["compound"]["link"][0]["database"] == "CAS"
    values = first["peak"]["peak"]["values"]
    assert {"mz", "intensity", "rel"}.issubset(values[0])


def test_massbank_limit_is_client_bounded(monkeypatch):
    payload = [{"accession": f"A{i}", "peak": {"peak": {"values": []}}, "compound": {"link": []}} for i in range(40)]

    class Response:
        encoding = "utf-8"
        def raise_for_status(self): pass
        def iter_content(self, chunk_size=0): yield json.dumps(payload).encode()

    monkeypatch.setattr("phyto_reason.tools.spectra_tools.requests.get", lambda *a, **k: Response())
    result = MassBankSpectrumSearchTool().run(compound_name="quercetin", limit=1)
    assert result.status == "partial"
    assert result.metadata["count"] == 1
    assert "limit" in result.warnings[0]  # 措辞随 6.4 Step 2 复测结论更新


def test_mona_requires_nonempty_query():
    result = MonaSpectrumSearchTool().run(query="")
    assert result.status == "partial"
    assert "至少包含" in result.warnings[0]


def test_pubchem_failure_uses_explicit_local_fallback(monkeypatch):
    def fail(*args, **kwargs):
        raise TimeoutError("timed out")

    monkeypatch.setattr("phyto_reason.tools.spectra_tools.requests.get", fail)
    tool = PubChemCompoundPropertiesTool()
    first = tool.run(name="quercetin")
    second = tool.run(name="quercetin")
    assert first.status == "partial"
    assert second.status == "partial"
    assert first.metadata["count"] == 1
    assert "降级" in first.warnings[0]


def test_pubchem_profile_fallback_uses_shared_compound_profiles(monkeypatch):
    def fail(*args, **kwargs):
        raise TimeoutError("timed out")

    monkeypatch.setattr("phyto_reason.tools.spectra_tools.requests.get", fail)
    tool = PubChemCompoundPropertiesTool()
    result = tool.run(name="berberine")
    assert result.status == "partial"
    assert result.metadata["count"] == 1
    record = result.metadata["records"][0]
    assert record["source"] == "compound_profiles"
    assert record["profile_class"] == "berberine_protoberberine"
    assert record["MolecularFormula"] == "C20H18NO4"
    assert "CAS" not in record
    assert "SMILES" not in record


def test_pubchem_unknown_profile_fallback_is_empty_and_explicit(monkeypatch):
    def fail(*args, **kwargs):
        raise TimeoutError("timed out")

    monkeypatch.setattr("phyto_reason.tools.spectra_tools.requests.get", fail)
    result = PubChemCompoundPropertiesTool().run(name="not_in_compound_profiles")
    assert result.status == "partial"
    assert result.metadata["count"] == 0
    assert result.metadata["records"] == []
    assert "未收录" in result.warnings[-1]
