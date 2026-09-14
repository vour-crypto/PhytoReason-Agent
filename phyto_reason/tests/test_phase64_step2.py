"""Phase 6.4 Step 2: Legacy 收口、TLS 恢复、长表上传与谱图工具挂账测试。"""

from __future__ import annotations

import json
import logging
from pathlib import Path

import pytest
import requests
from fastapi.testclient import TestClient


# ── Legacy 配置收口 ─────────────────────────────────────────

def test_workbench_db_uses_platform_paths(monkeypatch, tmp_path):
    monkeypatch.setenv("PHYTOREASON_DATA_DIR", str(tmp_path / "ud"))
    from phyto_reason.desktop import workbench_db
    from phyto_reason.platform_paths import workbench_db_path

    assert workbench_db.default_workbench_path() == workbench_db_path()
    db = workbench_db.WorkbenchDB()
    assert db.path == workbench_db_path()


def test_legacy_llm_settings_priority_config_over_db_over_env(monkeypatch, tmp_path):
    monkeypatch.setenv("PHYTOREASON_DATA_DIR", str(tmp_path / "ud"))
    from phyto_reason.config import llm_config
    from phyto_reason.desktop.client import _legacy_llm_settings
    from phyto_reason.desktop.workbench_db import WorkbenchDB

    db = WorkbenchDB(":memory:")
    monkeypatch.setenv("OPENAI_API_KEY", "sk-env")
    monkeypatch.setenv("OPENAI_BASE_URL", "https://env.example/v1")
    monkeypatch.delenv("LLM_MODEL", raising=False)

    settings = _legacy_llm_settings(db)
    assert settings["api_key"] == "sk-env"
    assert settings["base_url"] == "https://env.example/v1"
    assert settings["model"] == "deepseek-chat"

    db.set_setting("api_key", "sk-db")
    db.set_setting("model", "db-model")
    settings = _legacy_llm_settings(db)
    assert settings["api_key"] == "sk-db"
    assert settings["model"] == "db-model"
    assert settings["base_url"] == "https://env.example/v1"

    llm_config.save_config({
        "provider": "custom", "base_url": "https://file.example/v1",
        "model": "file-model", "api_key": "sk-file",
    })
    settings = _legacy_llm_settings(db)
    assert settings["api_key"] == "sk-file"
    assert settings["model"] == "file-model"
    assert settings["base_url"] == "https://file.example/v1"


def test_record_probe_never_persists_env_secrets(monkeypatch, tmp_path):
    monkeypatch.setenv("PHYTOREASON_DATA_DIR", str(tmp_path / "ud"))
    monkeypatch.setenv("OPENAI_API_KEY", "sk-env-secret-9999")
    monkeypatch.setenv("LLM_MODEL", "env-model")
    from phyto_reason.config import llm_config

    llm_config.record_probe({"connectivity": False, "error": "probe only"})
    persisted = json.loads(llm_config.config_path().read_text(encoding="utf-8"))
    assert persisted["probe"]["error"] == "probe only"
    dumped = json.dumps(persisted)
    assert "sk-env-secret-9999" not in dumped
    assert "env-model" not in dumped


# ── TLS 校验恢复 ────────────────────────────────────────────

def test_kegg_calls_verify_tls_and_degrade_structurally(monkeypatch, caplog):
    from phyto_reason.reasoning.pathway_activity_engine import PathwayActivityEngine

    calls: list[dict] = []

    def fake_get(*args, **kwargs):
        calls.append(kwargs)
        raise requests.exceptions.SSLError("certificate verify failed")

    monkeypatch.setattr("requests.get", fake_get)
    monkeypatch.setattr("time.sleep", lambda seconds: None)
    with caplog.at_level(logging.WARNING, logger="phyto_reason.reasoning.pathway_activity_engine"):
        result = PathwayActivityEngine._fetch_kegg_enzymes("berberine")

    assert result == []
    assert calls and all("verify" not in kwargs for kwargs in calls)
    assert any("degraded" in record.getMessage() for record in caplog.records)


def test_pathway_engine_source_has_no_verify_false():
    from phyto_reason.reasoning import pathway_activity_engine

    source = Path(pathway_activity_engine.__file__).read_text(encoding="utf-8")
    assert "verify=False" not in source


# ── MassBank 分页语义 / MoNA 语法提示 ───────────────────────

def _massbank_payload(n: int, duplicate_first: bool = False) -> list[dict]:
    records = [
        {
            "accession": f"MSBNK-X-{i:03d}", "title": f"t{i}",
            "compound": {"names": [{"name": "quercetin"}], "link": []},
            "peak": {"peak": {"values": [{"mz": 100 + i, "intensity": 10, "rel": 5}]}},
        }
        for i in range(n)
    ]
    if duplicate_first:
        records.append(dict(records[0]))
    return records


def _mock_stream_response(monkeypatch, payload):
    class Response:
        encoding = "utf-8"

        def raise_for_status(self) -> None:
            pass

        def iter_content(self, chunk_size: int = 0):
            yield json.dumps(payload).encode()

    monkeypatch.setattr(
        "phyto_reason.tools.spectra_tools.requests.get", lambda *a, **k: Response()
    )


def test_massbank_reports_total_has_more_and_dedups(monkeypatch):
    import phyto_reason.tools.spectra_tools as st

    _mock_stream_response(monkeypatch, _massbank_payload(12, duplicate_first=True))
    result = st.MassBankSpectrumSearchTool().run(compound_name="quercetin", limit=10)
    assert result.metadata["total_count"] == 12
    assert result.metadata["duplicates_dropped"] == 1
    assert result.metadata["has_more"] is True
    assert result.metadata["returned"] == 10
    assert result.status == "partial"
    assert any("不执行 limit" in warning for warning in result.warnings)


def test_massbank_has_more_false_when_all_returned(monkeypatch):
    import phyto_reason.tools.spectra_tools as st

    _mock_stream_response(monkeypatch, _massbank_payload(3))
    result = st.MassBankSpectrumSearchTool().run(compound_name="quercetin", limit=10)
    assert result.metadata["has_more"] is False
    assert result.metadata["total_count"] == 3
    assert result.status == "success"
    assert result.warnings == []


def test_mona_surfaces_server_syntax_hint(monkeypatch):
    import phyto_reason.tools.spectra_tools as st

    class Resp:
        status_code = 400

        def json(self):
            return {"error": "invalid query syntax", "hint": "queries use spring-filter syntax, for example exists(metaData.name:'ionization mode')"}

        def raise_for_status(self):
            raise requests.HTTPError("400 Client Error", response=self)

    monkeypatch.setattr("phyto_reason.tools.spectra_tools.requests.get", lambda *a, **k: Resp())
    result = st.MonaSpectrumSearchTool().run(query='metaData.qcString=="Quercetin"')
    assert result.status == "partial"
    assert any("spring-filter" in warning for warning in result.warnings)


# ── 长表 m/z-RT-intensity ───────────────────────────────────

def _write_csv(tmp_path: Path, name: str, text: str) -> Path:
    path = tmp_path / name
    path.write_text(text, encoding="utf-8")
    return path


def test_long_table_pivots_to_wide_matrix(tmp_path):
    from phyto_reason.ingestion.parsers.metabolite_parser import parse_metabolite

    path = _write_csv(tmp_path, "long.csv",
                      "sample,m/z,RT,intensity,formula\n"
                      "F-1,348.1230,3.25,1000,C21H18NO4\n"
                      "F-1,177.05,1.10,2000,\n"
                      "L-1,348.1230,3.25,1500,\n")
    parsed = parse_metabolite(path)
    assert parsed.n_samples == 2
    assert parsed.sample_ids == ["F-1", "L-1"]
    assert "mz_348.1230_rt_3.250" in parsed.matrix
    assert parsed.matrix["mz_348.1230_rt_3.250"]["F-1"] == 1000.0
    assert parsed.matrix["mz_348.1230_rt_3.250"]["L-1"] == 1500.0
    assert any("long table" in warning.lower() for warning in parsed.warnings)


def test_long_table_trio_without_sample_column_is_refused(tmp_path):
    from phyto_reason.ingestion.parsers.metabolite_parser import parse_metabolite

    path = _write_csv(tmp_path, "trio.csv",
                      "mz,RT,intensity\n348.1,3.2,1000\n177.0,1.1,2000\n")
    with pytest.raises(ValueError, match="sample 列"):
        parse_metabolite(path)


def test_wide_peak_table_with_mz_rt_columns_zero_regression(tmp_path):
    from phyto_reason.ingestion.parsers.metabolite_parser import parse_metabolite

    path = _write_csv(tmp_path, "wide.csv",
                      "Name,mz,RT,F-1,L-1\n"
                      "berberine,336.1,2.5,100,200\n"
                      "palmatine,352.1,2.9,300,400\n")
    parsed = parse_metabolite(path)
    assert parsed.sample_ids == ["F-1", "L-1"]
    assert parsed.n_features == 2
    assert "berberine" in parsed.matrix
    assert parsed.matrix["berberine"]["F-1"] == 100.0


def test_upload_route_sniffs_long_table_and_honors_explicit_type(tmp_path):
    from phyto_reason.api.app import app

    client = TestClient(app)
    content = "sample,mz,RT,intensity\nF-1,348.1,3.2,100\nL-1,348.1,3.2,150\n"

    sniffed = client.post("/upload", files={"file": ("run12.csv", content.encode(), "text/csv")})
    assert sniffed.status_code == 200
    summary = sniffed.json()["summary"]
    assert summary["file_type"] == "metabolite_matrix"
    assert summary.get("detection") == "long_table_header_sniff"

    explicit = client.post(
        "/upload",
        files={"file": ("exp1.csv", content.encode(), "text/csv")},
        data={"data_type": "metabolite"},
    )
    assert explicit.json()["summary"]["file_type"] == "metabolite_matrix"

    forced_expression = client.post(
        "/upload",
        files={"file": ("whatever.csv", "S1,S2\n10,20\n".encode(), "text/csv")},
        data={"data_type": "expression"},
    )
    assert forced_expression.json()["summary"]["file_type"] == "expression_matrix"


# ── 前端槽位④与失败提示语义 ────────────────────────────────

def test_frontend_long_slot_enabled_and_failure_semantics_fixed():
    import phyto_reason.api.app as app_module

    static_dir = Path(app_module.__file__).resolve().parent / "static"
    html = (static_dir / "index.html").read_text(encoding="utf-8")
    assert 'data-slot="long"' in html
    assert 'class="slot disabled"' not in html
    assert "sample, m/z, RT, intensity" in html

    js = (static_dir / "assets" / "app.js").read_text(encoding="utf-8")
    assert "SLOT_BY_INTENT" in js
    assert "slots[3]" in js
    assert 'form.append("data_type"' in js
    assert "尚未对齐" in js
