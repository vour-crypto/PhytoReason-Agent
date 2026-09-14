"""Phase 6.4 Step 1: 设置页、数据目录合同与备份导入导出测试。"""

from __future__ import annotations

import io
import json
import zipfile
from pathlib import Path

import pytest
from fastapi.testclient import TestClient


@pytest.fixture()
def client(monkeypatch, tmp_path):
    monkeypatch.setenv("PHYTOREASON_DATA_DIR", str(tmp_path / "user-data"))
    # 隔离全局会话存储单例：指向临时 DB 并清空内存缓存，
    # 防止测试读写用户真实会话库（曾发生测试误删真实空会话）
    from phyto_reason.agent import conversation_store as cs
    monkeypatch.setattr(cs.store, "_db_path", str(tmp_path / "sessions.db"))
    monkeypatch.setattr(cs.store, "_sessions", {})
    cs.store._init_db()  # 在临时 DB 上建表
    from phyto_reason.api.app import app
    return TestClient(app)


# ── 数据目录合同 ────────────────────────────────────────────

def test_memory_db_path_follows_user_dir_and_memory_store_uses_it(monkeypatch, tmp_path):
    monkeypatch.setenv("PHYTOREASON_DATA_DIR", str(tmp_path / "user-data"))
    from phyto_reason.agent.memory_store import MemoryStore
    from phyto_reason.platform_paths import memory_db_path, sessions_db_path, user_data_dir

    assert user_data_dir() == tmp_path / "user-data"
    assert memory_db_path() == sessions_db_path()
    store = MemoryStore()
    assert Path(store._db_path) == memory_db_path()


def test_copy_if_missing_never_overwrites_existing_target(tmp_path):
    from phyto_reason.platform_paths import copy_if_missing

    source = tmp_path / "src" / "x.db"
    source.parent.mkdir(parents=True)
    source.write_text("new", encoding="utf-8")

    occupied = tmp_path / "dst" / "x.db"
    occupied.parent.mkdir(parents=True)
    occupied.write_text("old", encoding="utf-8")
    copy_if_missing(source, occupied)
    assert occupied.read_text(encoding="utf-8") == "old"

    fresh = tmp_path / "dst2" / "y.db"
    copy_if_missing(source, fresh)
    assert fresh.read_text(encoding="utf-8") == "new"


def test_migration_marker_written_once(monkeypatch, tmp_path):
    monkeypatch.setenv("PHYTOREASON_DATA_DIR", str(tmp_path / "user-data"))
    from phyto_reason.platform_paths import user_data_dir

    root = user_data_dir()
    marker = root / ".migration_v1.complete"
    assert marker.exists()
    first = marker.read_text(encoding="utf-8")
    assert user_data_dir() == root
    assert marker.read_text(encoding="utf-8") == first


# ── 设置页接线 ──────────────────────────────────────────────

def test_settings_page_wired_and_model_card_moved():
    import phyto_reason.api.app as app_module

    static_dir = Path(app_module.__file__).resolve().parent / "static"
    html = (static_dir / "index.html").read_text(encoding="utf-8")
    assert 'id="page-settings"' in html
    assert 'data-page="settings"' in html
    assert "本期未启用" not in html
    assert 'id="gear-btn"' in html
    assert 'id="llm-settings-card"' in html
    data_section = html.split('id="page-data"', 1)[1].split('<div class="page"', 1)[0]
    assert "llm-settings-card" not in data_section
    assert "数据与备份" in html and "backup-import" in html

    js = (static_dir / "assets" / "app.js").read_text(encoding="utf-8")
    assert "bindSettings" in js and "/settings/backup/import" in js
    assert '"dashboard", "analysis", "hyp", "data", "settings"' in js


# ── 数据目录 API ────────────────────────────────────────────

def test_data_info_lists_paths_and_never_leaks_api_key(client):
    from phyto_reason.config import llm_config

    llm_config.save_config({
        "provider": "custom", "base_url": "https://example.test/v1",
        "model": "m", "api_key": "sk-verysecret-9876",
    })
    response = client.get("/settings/data")
    assert response.status_code == 200
    body = json.dumps(response.json())
    payload = response.json()
    assert payload["paths"][0]["label"] == "数据目录"
    assert Path(payload["data_dir"]).is_absolute()
    assert "sk-verysecret" not in body
    assert "api_key" not in body


def test_open_data_dir_calls_file_manager(monkeypatch, client):
    import phyto_reason.api.routes.settings as settings_routes
    from phyto_reason.platform_paths import user_data_dir

    calls: list[Path] = []
    monkeypatch.setattr(settings_routes, "_open_in_file_manager", lambda p: calls.append(Path(p)))
    response = client.post("/settings/data/open", json={"target": "data"})
    assert response.status_code == 200
    assert calls and calls[0] == user_data_dir()

    bad = client.post("/settings/data/open", json={"target": "nope"})
    assert bad.status_code == 400


# ── 备份导出/导入 ───────────────────────────────────────────

def _export_backup(client, target_dir: Path) -> Path:
    response = client.post("/settings/backup/export", json={"target_dir": str(target_dir)})
    assert response.status_code == 200, response.text
    return Path(response.json()["path"])


def _import_zip_bytes(client, data: bytes):
    return client.post(
        "/settings/backup/import",
        files={"file": ("backup.zip", data, "application/zip")},
    )


def test_backup_export_import_roundtrip_with_bak_protection(client, tmp_path):
    from phyto_reason.config import llm_config
    from phyto_reason.platform_paths import user_data_dir

    llm_config.save_config({
        "provider": "deepseek", "base_url": "https://api.deepseek.com/v1",
        "model": "roundtrip-model", "api_key": "sk-secret-1234",
    })
    probe = user_data_dir() / "data" / "probe.txt"
    probe.parent.mkdir(parents=True, exist_ok=True)
    probe.write_text("v1", encoding="utf-8")

    out_dir = tmp_path / "out"
    out_dir.mkdir(parents=True, exist_ok=True)
    zip_path = _export_backup(client, out_dir)

    # 破坏现场后再导入
    probe.write_text("v2", encoding="utf-8")
    llm_config.save_config({
        "provider": "deepseek", "base_url": "https://api.deepseek.com/v1",
        "model": "mutated-model", "api_key": "",
    })

    response = _import_zip_bytes(client, zip_path.read_bytes())
    assert response.status_code == 200, response.text
    assert response.json()["restored"] > 0
    assert probe.read_text(encoding="utf-8") == "v1"
    bak = Path(str(probe) + ".bak")
    assert bak.exists() and bak.read_text(encoding="utf-8") == "v2"
    assert llm_config.load_config()["model"] == "roundtrip-model"
    assert llm_config.public_config()["api_key"] == "sk-***1234"


def test_backup_import_rejects_zip_slip(client, tmp_path):
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w") as zf:
        zf.writestr("../evil.txt", "pwned")
    response = _import_zip_bytes(client, buffer.getvalue())
    assert response.status_code == 400
    assert "绝对路径" in response.json()["error"] or ".." in response.json()["error"]
    assert not (tmp_path / "evil.txt").exists()
    assert not (tmp_path / "user-data" / "evil.txt").exists()


def test_backup_import_enforces_decompressed_limit(client, monkeypatch):
    import phyto_reason.api.routes.settings as settings_routes

    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w", zipfile.ZIP_DEFLATED) as zf:
        zf.writestr("big.txt", "A" * 5000)
    monkeypatch.setattr(settings_routes, "MAX_IMPORT_BYTES", 1024)
    response = _import_zip_bytes(client, buffer.getvalue())
    assert response.status_code == 400
    assert "上限" in response.json()["error"]


def test_backup_bak_rotates_instead_of_stacking(client, tmp_path):
    from phyto_reason.platform_paths import user_data_dir

    target = user_data_dir() / "x.txt"
    target.parent.mkdir(parents=True, exist_ok=True)

    def zip_with(content: str) -> bytes:
        buffer = io.BytesIO()
        with zipfile.ZipFile(buffer, "w") as zf:
            zf.writestr("x.txt", content)
        return buffer.getvalue()

    target.write_text("A", encoding="utf-8")
    assert _import_zip_bytes(client, zip_with("B")).status_code == 200
    assert target.read_text(encoding="utf-8") == "B"
    assert Path(str(target) + ".bak").read_text(encoding="utf-8") == "A"

    assert _import_zip_bytes(client, zip_with("C")).status_code == 200
    assert target.read_text(encoding="utf-8") == "C"
    assert Path(str(target) + ".bak").read_text(encoding="utf-8") == "B"
    assert not Path(str(target) + ".bak.bak").exists()


# ── 静态路由优先级 ─────────────────────────────────────────

def test_user_figures_route_serves_mount_directory(client):
    from phyto_reason.platform_paths import figures_dir

    probe = figures_dir() / "phase64_probe.png"
    probe.write_bytes(b"\x89PNG-probe")
    try:
        response = client.get("/static/figures/phase64_probe.png")
        assert response.status_code == 200
        assert response.content == b"\x89PNG-probe"
    finally:
        probe.unlink(missing_ok=True)


# ── 产物目录自定义（用户反馈：422 / 自定义保存位置）──────────

def test_directories_override_roundtrip(monkeypatch, tmp_path):
    from phyto_reason.platform_paths import (
        dir_overrides, figures_dir, report_dir, set_dir_overrides,
    )

    custom_figs = tmp_path / "custom" / "figs"
    custom_reports = tmp_path / "custom" / "reports"
    saved = set_dir_overrides({"figures_dir": str(custom_figs), "reports_dir": str(custom_reports)})
    assert saved["figures_dir"] == str(custom_figs)
    assert figures_dir() == custom_figs
    assert report_dir() == custom_reports
    assert dir_overrides()["reports_dir"] == str(custom_reports)

    with pytest.raises(ValueError, match="绝对路径"):
        set_dir_overrides({"figures_dir": "relative/path"})

    set_dir_overrides({"figures_dir": "", "reports_dir": ""})
    assert figures_dir() != custom_figs
    assert dir_overrides() == {}


def test_directories_api_roundtrip(client, tmp_path):
    response = client.get("/settings/directories")
    assert response.status_code == 200
    assert Path(response.json()["backup_default"]).is_absolute()

    custom = tmp_path / "chosen" / "figs"
    saved = client.post("/settings/directories", json={"figures_dir": str(custom), "reports_dir": ""})
    assert saved.status_code == 200
    assert saved.json()["saved"]["figures_dir"] == str(custom)

    rejected = client.post("/settings/directories", json={"figures_dir": "nope", "reports_dir": ""})
    assert rejected.status_code == 400
    assert "绝对路径" in rejected.json()["error"]


def test_export_backup_to_custom_target_dir(client, tmp_path):
    from phyto_reason.platform_paths import user_data_dir

    (user_data_dir() / "data").mkdir(parents=True, exist_ok=True)
    (user_data_dir() / "data" / "marker.txt").write_text("v1", encoding="utf-8")
    target = tmp_path / "chosen_backup"
    target.mkdir(parents=True, exist_ok=True)

    response = client.post("/settings/backup/export", json={"target_dir": str(target)})
    assert response.status_code == 200
    zip_path = Path(response.json()["path"])
    assert zip_path.parent == target
    assert zip_path.exists()

    restored = client.post("/settings/backup/import", files={"file": ("b.zip", zip_path.read_bytes(), "application/zip")})
    assert restored.status_code == 200
    assert restored.json()["restored"] >= 1


def test_frontend_directory_inputs_and_content_type_wired():
    import phyto_reason.api.app as app_module

    static_dir = Path(app_module.__file__).resolve().parent / "static"
    html = (static_dir / "index.html").read_text(encoding="utf-8")
    assert 'id="figures-dir-input"' in html
    assert 'id="reports-dir-input"' in html
    assert 'id="backup-dir-input"' in html
    assert 'id="dirs-save"' in html

    js = (static_dir / "assets" / "app.js").read_text(encoding="utf-8")
    assert '"Content-Type": "application/json"' in js
    assert "/settings/directories" in js
    assert "target_dir: targetDir" in js


# ── 上传类型菜单 / 目录选择弹窗 / 备份默认位置（用户反馈批次 3）──

def test_default_backup_dir_prefers_downloads():
    from phyto_reason.platform_paths import default_backup_dir

    downloads = Path.home() / "Downloads"
    if downloads.is_dir():
        assert default_backup_dir() == downloads
    else:
        from phyto_reason.platform_paths import user_data_dir
        assert default_backup_dir() == user_data_dir() / "backups"


def test_export_backup_default_uses_configured_dir(client, monkeypatch, tmp_path):
    import phyto_reason.api.routes.settings as settings_routes
    from phyto_reason.platform_paths import user_data_dir

    (user_data_dir() / "data").mkdir(parents=True, exist_ok=True)
    (user_data_dir() / "data" / "marker.txt").write_text("v1", encoding="utf-8")
    fake_downloads = tmp_path / "downloads"
    fake_downloads.mkdir()
    monkeypatch.setattr(settings_routes, "default_backup_dir", lambda: fake_downloads)

    response = client.post("/settings/backup/export", json={})
    assert response.status_code == 200
    assert Path(response.json()["path"]).parent == fake_downloads


def test_directory_picker_endpoint(client, monkeypatch):
    import phyto_reason.api.routes.settings as settings_routes

    monkeypatch.setattr(settings_routes, "_native_directory_picker", lambda: "D:/chosen/dir")
    picked = client.post("/settings/directories/picker")
    assert picked.status_code == 200
    assert picked.json() == {"picked": True, "path": "D:\chosen\dir"}

    monkeypatch.setattr(settings_routes, "_native_directory_picker", lambda: None)
    cancelled = client.post("/settings/directories/picker")
    assert cancelled.status_code == 501


def test_frontend_upload_menu_and_browse_buttons():
    import phyto_reason.api.app as app_module

    static_dir = Path(app_module.__file__).resolve().parent / "static"
    html = (static_dir / "index.html").read_text(encoding="utf-8")
    assert "点击各行槽位上传对应文件" in html
    assert html.count('class="btn btn-secondary btn-sm browse-btn"') == 3

    js = (static_dir / "assets" / "app.js").read_text(encoding="utf-8")
    assert "function chooseUploadType" in js
    assert "/settings/directories/picker" in js
    assert 'browseButton.dataset.browse' in js


def test_frontend_start_analysis_and_thinking_indicator():
    import phyto_reason.api.app as app_module

    static_dir = Path(app_module.__file__).resolve().parent / "static"
    html = (static_dir / "index.html").read_text(encoding="utf-8")
    assert 'id="start-analysis"' in html
    assert 'id="start-analysis-hint"' in html
    assert 'id="analysis-picker"' in html
    assert 'data-analysis="o2pls"' in html
    # 对话驱动流程：物种/目标代谢物输入 + 按钮只生成草稿不自动分析
    assert 'id="upload-species"' in html
    assert 'id="upload-target"' in html
    assert "生成分析指令" in html
    assert "在对话框确认/补充" in html
    assert "自动跳转分析页执行" not in html

    js = (static_dir / "assets" / "app.js").read_text(encoding="utf-8")
    assert "function updateStartButton()" in js
    assert "function buildAnalysisPrompt(selected)" in js
    assert "function refreshAnalysisOptions()" in js
    assert "function selectedAnalyses()" in js
    assert "DEFAULT_ANALYSES" in js
    assert "uploadedKinds" in js
    assert "分析中，请稍候" in js
    assert "🧠 大模型规划分析中" not in js
    # 用户要求：不得默认一键跑全管线
    assert "请运行完整多组学分析管线" not in js
    assert "state.streaming" in js
    assert "startThinkingIndicator" in js and "stopThinkingIndicator" in js
    # 草稿模式：填入后聚焦，不自动发送
    assert "input.value = buildAnalysisPrompt(selected);" in js
    assert "input.focus();" in js
    assert 'form.append("species", speciesInput.value.trim())' in js
    assert "form.append(\"target_metabolite\", targetInput.value.trim())" in js
    assert "请在此说明你的物种" in js


def test_session_and_figure_clear_endpoints(client):
    from phyto_reason.visualization._common import new_figure, save_figure

    created = client.post("/sessions")
    assert created.status_code == 200
    session_id = created.json()["session_id"]
    assert session_id

    listed = client.get("/sessions")
    assert listed.status_code == 200
    assert any(s["session_id"] == session_id for s in listed.json()["sessions"])

    import matplotlib.pyplot as plt

    fig, ax = new_figure("Clear Probe")
    ax.plot([1], [1])
    url = save_figure(fig, "clear_probe")
    plt.close(fig)
    name = url.rsplit("/", 1)[-1]
    from phyto_reason.platform_paths import figures_dir
    assert (figures_dir() / name).exists()

    cleared = client.delete("/figures")
    assert cleared.status_code == 200
    assert cleared.json()["deleted"] >= 1
    assert not (figures_dir() / name).exists()
    assert not (figures_dir() / name.replace(".png", ".svg")).exists()


def test_frontend_session_list_and_figure_clear_wired():
    import phyto_reason.api.app as app_module

    static_dir = Path(app_module.__file__).resolve().parent / "static"
    html = (static_dir / "index.html").read_text(encoding="utf-8")
    assert 'id="new-session"' in html
    assert 'id="session-list"' in html and 'id="sessions-empty"' in html
    assert 'id="figures-clear"' in html

    js = (static_dir / "assets" / "app.js").read_text(encoding="utf-8")
    assert "function switchSession(" in js
    assert "async function createNewSession()" in js
    assert 'fetch("/sessions"' in js
    assert 'fetch("/figures", { method: "DELETE" })' in js
    assert "function resetChatArea()" in js


# ── 论文式图号 + 按实际分析分节的报告（用户反馈批次 8）────────

def _fake_pipeline_state():
    from types import SimpleNamespace
    return SimpleNamespace(
        completed_nodes=["data_qc", "deg_analysis", "dam_analysis"],
        deg_report={"n_genes_tested": 3000, "n_fdr_significant": 2522,
                    "summary": "n=15; method=ANOVA; q<0.05",
                    "normalization": "median_of_ratios_size_factors + log2",
                    "figure_url": "/static/figures/volcano_deg_abc123.png"},
        dam_report={"n_metabolites_tested": 697, "n_fdr_significant": 652,
                    "figure_url": "/static/figures/volcano_dam_def456.png"},
        qc_report={},
    )


def test_analysis_digest_only_covers_completed_sections():
    from phyto_reason.reports.analysis_report import build_analysis_digest, render_markdown

    digest = build_analysis_digest(_fake_pipeline_state(), species="Zanthoxylum nitidum",
                                   target_metabolite="nitidine", figure_registry={})
    section_titles = [s["title_zh"] for s in digest["sections"]]
    assert any("差异表达" in title for title in section_titles)
    assert any("差异代谢物" in title for title in section_titles)
    # 未运行的 WGCNA/O2PLS 绝不出现在分节里
    assert not any("WGCNA" in title or "O2PLS" in title for title in section_titles)
    assert digest["figure_urls"] == ["/static/figures/volcano_deg_abc123.png",
                                     "/static/figures/volcano_dam_def456.png"]

    markdown = render_markdown(digest)
    assert "**Figure 1**" in markdown and "差异表达火山图" in markdown
    assert "**Figure 2**" in markdown
    assert "如何解读" in markdown
    assert "median_of_ratios_size_factors" in markdown
    assert "WGCNA" not in markdown


def test_figure_numbering_is_stable_and_incremental():
    from phyto_reason.reports.analysis_report import build_analysis_digest

    registry = {}
    first = build_analysis_digest(_fake_pipeline_state(), figure_registry=registry)
    assert [f["no"] for f in first["sections"][0]["figures"]] == [1]
    assert [f["no"] for f in first["sections"][1]["figures"]] == [2]
    # 第二次运行：新图编号继续递增，旧图编号不变
    state2 = _fake_pipeline_state()
    state2.deg_report["figure_url"] = "/static/figures/volcano_deg_new000.png"
    second = build_analysis_digest(state2, figure_registry=registry)
    # 旧文件编号在注册表中保持不变；新文件继续递增
    assert registry["volcano_deg_abc123.png"]["no"] == 1
    assert registry["volcano_dam_def456.png"]["no"] == 2
    numbers = {Path(f["url"]).name: f["no"] for s in second["sections"] for f in s["figures"]}
    assert numbers["volcano_deg_new000.png"] == 3


def test_report_export_endpoint(client):
    from phyto_reason.agent.conversation_store import store

    created = client.post("/sessions")
    session_id = created.json()["session_id"]
    empty = client.post("/report/export", json={"session_id": session_id})
    assert empty.status_code == 409

    session = store.get_or_create(session_id)
    session.analysis_digest = {
        "created_at": "2026-09-13T00:00:00", "species": "Test", "target_metabolite": "",
        "sections": [{"node": "deg_analysis", "title_zh": "差异表达分析（DEG）",
                      "title_en": "DEG", "stats": {"n_fdr_significant": 10},
                      "figures": [{"no": 1, "url": "/static/figures/x.png", "file": "x.png",
                                   "caption_en": "DEG volcano", "caption_zh": "火山图",
                                   "how_to_read": "右上为上调"}]}],
        "figure_urls": ["/static/figures/x.png"],
    }
    response = client.post("/report/export", json={"session_id": session_id})
    assert response.status_code == 200
    assert "Figure 1" in response.json()["markdown"]
    assert response.json()["filename"].endswith(".md")


def test_frontend_report_export_and_captions_wired():
    import phyto_reason.api.app as app_module

    static_dir = Path(app_module.__file__).resolve().parent / "static"
    js = (static_dir / "assets" / "app.js").read_text(encoding="utf-8")
    assert 'fetch("/report/export"' in js
    assert "Figure " + '" + (index + 1) + "' in js
    assert "报告要求：回复仅覆盖上述所选步骤" in js
    assert "不要展开转录因子或假设内容" in js


# ── 会话删除/清空 + 假设池（用户反馈批次 10）────────────────

def test_session_delete_and_clear_empty_endpoints(client):
    created = client.post("/sessions")
    sid = created.json()["session_id"]
    other = client.post("/sessions").json()["session_id"]
    assert client.get("/sessions").json()["sessions"]  # 有会话存在

    deleted = client.delete(f"/sessions/{sid}")
    assert deleted.status_code == 200
    remaining = [s["session_id"] for s in client.get("/sessions").json()["sessions"]]
    assert sid not in remaining
    assert client.get(f"/session/{sid}").json().get("error") == "Session not found"

    cleared = client.post("/sessions/clear-empty", json={"keep": other})
    assert cleared.status_code == 200
    after = [s["session_id"] for s in client.get("/sessions").json()["sessions"]]
    assert other in after  # 保留当前会话


def test_clear_empty_keeps_data_sessions(client, tmp_path):
    from phyto_reason.agent.conversation_store import store

    with_data = client.post("/sessions").json()["session_id"]
    data_session = store.get_or_create(with_data)
    data_session.has_data = True
    data_session._persist()  # has_data 落库（与 UI 上传路径一致）
    empty_one = client.post("/sessions").json()["session_id"]
    empty_two = client.post("/sessions").json()["session_id"]

    cleared = client.post("/sessions/clear-empty", json={"keep": ""})
    assert cleared.status_code == 200
    remaining = [s["session_id"] for s in client.get("/sessions").json()["sessions"]]
    assert with_data in remaining
    assert empty_one not in remaining and empty_two not in remaining


def test_hypotheses_endpoint_roundtrip(client):
    from phyto_reason.agent.conversation_store import store

    sid = client.post("/sessions").json()["session_id"]
    session = store.get_or_create(sid)
    session.set_pipeline_result("result text", [
        {"title": "MYB7 协同调控小檗碱积累", "mechanism_type": "enzymatic_regulation",
         "uncertainty_level": "moderate", "proposed_regulators": ["MYB7", "bHLH34"],
         "target_metabolites": ["berberine"]},
    ])
    response = client.get("/hypotheses", params={"session_id": sid})
    assert response.status_code == 200
    hypotheses = response.json()["hypotheses"]
    assert len(hypotheses) == 1
    assert hypotheses[0]["title"] == "MYB7 协同调控小檗碱积累"
    assert hypotheses[0]["uncertainty_level"] == "moderate"
    assert "MYB7" in hypotheses[0]["proposed_regulators"]


def test_frontend_session_cleanup_and_hyp_pool_wired():
    import phyto_reason.api.app as app_module

    static_dir = Path(app_module.__file__).resolve().parent / "static"
    html = (static_dir / "index.html").read_text(encoding="utf-8")
    assert 'id="clear-empty-sessions"' in html
    assert 'id="hyp-rating-filter"' in html and 'id="hyp-search-input"' in html
    assert 'id="hyp-list"' in html
    # 顶栏死组件已拆除
    assert "proj-switch" not in html and 'data-ctx="project"' not in html
    assert 'class="search"' not in html.split("</div>")[0] or "搜索项目" not in html

    js = (static_dir / "assets" / "app.js").read_text(encoding="utf-8")
    assert '"session-del"' in js
    assert 'fetch("/sessions/clear-empty"' in js
    assert 'fetch("/sessions/" + encodeURIComponent(session.session_id), { method: "DELETE" })' in js
    assert "function renderHypotheses()" in js and "async function refreshHypotheses()" in js
    assert 'fetch("/hypotheses?session_id="' in js


def test_tool_definitions_all_dispatched_or_registered():
    """一致性合同：LLM 可见的工具定义必须已接线（分发链或注册表）。

    防止再次出现"菜单可勾、调用必落空（Unknown tool）"的缺陷。
    """
    import re

    import phyto_reason
    from phyto_reason.tools import TOOL_REGISTRY

    defs_path = Path(phyto_reason.__file__).resolve().parent / "agent" / "tool_definitions.py"
    orch_path = Path(phyto_reason.__file__).resolve().parent / "agent" / "orchestrator.py"
    defined = set(re.findall(r'"name":\s*"([a-z_]+)"', defs_path.read_text(encoding="utf-8")))
    dispatched = set(re.findall(r'tool_name == "([a-z_]+)"', orch_path.read_text(encoding="utf-8")))
    registered = set(TOOL_REGISTRY.list_tools())
    missing = defined - dispatched - registered
    assert not missing, f"以下工具有定义但未接线: {sorted(missing)}"


def test_check_data_quality_session_handler_produces_figures(client):
    import numpy as np

    from phyto_reason.agent.conversation_store import store
    from phyto_reason.agent.tool_handlers import handle_check_data_quality

    rng = np.random.default_rng(2)
    samples = ["L-1", "L-2", "L-3", "R-1", "R-2", "R-3"]
    expr = {f"G{i}": {s: float(rng.normal(100 + i, 10)) for s in samples} for i in range(30)}
    session = store.get_or_create("qc-fig-test")
    session.expression_matrix = expr
    session.metabolite_matrix = {"M1": {s: float(rng.normal(50, 5)) for s in samples}}
    session.sample_metadata = {"samples": [
        {"sample_id": s, "condition": s.rsplit("-", 1)[0], "tissue": s.rsplit("-", 1)[0]}
        for s in samples
    ]}

    reply = handle_check_data_quality({}, session=session)
    assert "Unknown" not in reply and "Error" not in reply
    assert reply.count("/static/figures/") >= 1


def test_upload_tf_fasta_runs_prediction_and_injects_annotation(client, monkeypatch, tmp_path):
    """上传序列 FASTA → TF 预测 → 注释注入会话分析链（数据页槽位⑤）。"""
    import phyto_reason
    fixture = Path(phyto_reason.__file__).resolve().parent / "tests" / "fixtures" / "pfam_test_subset.hmm"
    monkeypatch.setenv("PHYTOREASON_PFAM_HMM", str(fixture))
    monkeypatch.setenv("PHYTOREASON_TF_CACHE", str(tmp_path / "tfcache"))

    # 从夹具 HMM 解析共识序列（必定高分命中）
    consensus = ""
    for line in fixture.read_text(encoding="ascii").splitlines():
        parts = line.split()
        if len(parts) >= 22 and parts[0].isdigit():
            letters = [ch for ch in parts if len(ch) == 1 and ch.isalpha() and ch != "-"]
            if letters:
                consensus += letters[-1]
    assert consensus

    content = ">geneW\n" + consensus + "\n"
    response = client.post(
        "/upload",
        files={"file": ("transcripts.fasta", content.encode(), "text/plain")},
        data={"data_type": "tf_fasta"},
    )
    assert response.status_code == 200
    summary = response.json()["summary"]
    assert summary["file_type"] == "tf_prediction"
    assert summary["n_sequences"] == 1
    assert summary["n_predicted_tf"] >= 1

    from phyto_reason.agent.conversation_store import store
    from phyto_reason.agent.tool_handlers import _make_runtime_state

    sid = response.json()["session_id"]
    session = store.get(sid)
    assert Path(session.tf_annotation_csv).exists()
    state = _make_runtime_state(session)
    assert state.annotation_index is not None
    terms = state.annotation_index.get_terms("geneW")
    assert any("TF family" in term for term in terms)


def test_frontend_tf_fasta_slot_wired():
    import phyto_reason.api.app as app_module

    static_dir = Path(app_module.__file__).resolve().parent / "static"
    html = (static_dir / "index.html").read_text(encoding="utf-8")
    assert 'data-slot="tffasta"' in html
    assert "序列 FASTA（TF 预测）" in html

    js = (static_dir / "assets" / "app.js").read_text(encoding="utf-8")
    assert '"tffasta"' in js
    assert "tf_fasta" in js
    assert "正在预测转录因子" in js


def test_upload_ms2_spectra_annotation(client):
    """上传 MS/MS 谱图（MGF）→ 注释（MSI 分级 + 物种锚定）→ 注释表落盘。"""
    mgf = (
        "BEGIN IONS\nTITLE=scan=1\nPEPMASS=336.1230\nCHARGE=1+\n"
        "336.1230 100.0\n321.0996 80.0\n306.0761 50.0\n278.0812 30.0\nEND IONS\n"
        "BEGIN IONS\nTITLE=scan=2\nPEPMASS=303.0499\nCHARGE=1+\n"
        "303.0499 100.0\n153.0182 60.0\n137.0233 40.0\nEND IONS\n"
    )
    response = client.post(
        "/upload",
        files={"file": ("query_spectra.mgf", mgf.encode(), "text/plain")},
        data={"data_type": "ms2", "species": "两面针"},
    )
    assert response.status_code == 200
    summary = response.json()["summary"]
    assert summary["file_type"] == "ms2_spectra"
    assert summary["n_spectra"] == 2
    assert summary["n_annotated"] == 2
    assert summary["species_anchored"] == "两面针"
    assert summary["msi_counts"]
    csv_path = Path(summary["annotation_csv"])
    assert csv_path.exists()
    text = csv_path.read_text(encoding="utf-8-sig")
    assert "precursor_mz" in text and "336.123" in text

    from phyto_reason.agent.conversation_store import store
    sid = response.json()["session_id"]
    assert store.get(sid).ms2_annotation_csv == str(csv_path)


def test_frontend_ms2_slot_wired():
    import phyto_reason.api.app as app_module

    static_dir = Path(app_module.__file__).resolve().parent / "static"
    html = (static_dir / "index.html").read_text(encoding="utf-8")
    assert 'data-slot="ms2"' in html
    assert "MS/MS 谱图（MGF/CSV/MSP）" in html
    js = (static_dir / "assets" / "app.js").read_text(encoding="utf-8")
    assert '"ms2"' in js
    assert "正在注释谱图" in js
    assert ".mgf,.msp" in js


def test_delete_all_sessions_endpoint(client):
    ids = [client.post("/sessions").json()["session_id"] for _ in range(3)]
    from phyto_reason.agent.conversation_store import store
    store.get_or_create(ids[0]).has_data = True
    store.get_or_create(ids[0])._persist()

    wiped = client.delete("/sessions")
    assert wiped.status_code == 200
    assert wiped.json()["deleted_count"] == 3
    assert client.get("/sessions").json()["sessions"] == []
    assert store.get(ids[0]) is None


def test_frontend_clear_all_wired():
    import phyto_reason.api.app as app_module

    static_dir = Path(app_module.__file__).resolve().parent / "static"
    html = (static_dir / "index.html").read_text(encoding="utf-8")
    assert 'id="clear-all-sessions"' in html
    assert "不可恢复" in html

    js = (static_dir / "assets" / "app.js").read_text(encoding="utf-8")
    assert 'fetch("/sessions", { method: "DELETE" })' in js
    assert "此操作不可恢复" in js

    css = (static_dir / "assets" / "app.css").read_text(encoding="utf-8")
    assert "#session-list" in css and "overflow-y:auto" in css
