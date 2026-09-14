from __future__ import annotations

import json
from pathlib import Path

from fastapi.testclient import TestClient


def test_user_paths_follow_environment_override(monkeypatch, tmp_path):
    monkeypatch.setenv("PHYTOREASON_DATA_DIR", str(tmp_path / "user-data"))
    import phyto_reason.platform_paths as paths
    assert paths.user_data_dir() == tmp_path / "user-data"
    assert paths.figures_dir().parent == paths.user_data_dir()
    assert paths.sessions_db_path().parent == paths.user_data_dir() / "data"


def test_settings_masks_key_and_rejects_remote_http(monkeypatch, tmp_path):
    monkeypatch.setenv("PHYTOREASON_DATA_DIR", str(tmp_path / "user-data"))
    from phyto_reason.config import llm_config
    assert llm_config.save_config({"provider": "custom", "base_url": "https://example.test/v1", "model": "x", "api_key": "sk-secret-1234"})["api_key"] == "sk-secret-1234"
    assert llm_config.public_config()["api_key"] == "sk-***1234"
    try:
        llm_config.validate_config({"provider": "custom", "base_url": "http://example.test/v1", "model": "x", "api_key": "x"})
    except ValueError:
        pass
    else:
        raise AssertionError("remote http URL must be rejected")


def test_dynamic_figures_route_precedes_generic_static(monkeypatch, tmp_path):
    """图表走动态路由（目录覆盖立即生效），且优先于通用 /static。"""
    monkeypatch.setenv("PHYTOREASON_DATA_DIR", str(tmp_path / "user-data"))
    import phyto_reason.api.app as app_module
    routes = [getattr(route, "path", "") for route in app_module.app.routes]
    figure_routes = [i for i, path in enumerate(routes) if path.startswith("/static/figures")]
    assert figure_routes, "动态图表路由不存在"
    assert min(figure_routes) < routes.index("/static")


def test_figures_served_from_override_immediately(monkeypatch, tmp_path):
    """目录覆盖保存后无需重启：动态路由按请求解析新目录。"""
    monkeypatch.setenv("PHYTOREASON_DATA_DIR", str(tmp_path / "user-data"))
    from fastapi.testclient import TestClient
    from phyto_reason.api.app import app
    from phyto_reason.platform_paths import set_dir_overrides
    from phyto_reason.visualization._common import new_figure, save_figure

    client = TestClient(app)
    import matplotlib.pyplot as plt
    fig, ax = new_figure("Override Probe")
    ax.plot([1], [1])
    url = save_figure(fig, "override_probe")
    plt.close(fig)
    # 默认目录可访问
    assert client.get(url).status_code == 200

    custom = tmp_path / "custom" / "figs"
    set_dir_overrides({"figures_dir": str(custom)})
    # 立即生效：默认目录的旧图 404，新写入 custom 的图可访问
    assert client.get(url).status_code == 404
    fig2, ax2 = new_figure("Override Probe 2")
    ax2.plot([1], [1])
    url2 = save_figure(fig2, "override_probe2")
    plt.close(fig2)
    assert url2.startswith("/static/figures/")
    assert client.get(url2).status_code == 200
    assert (custom / url2.rsplit("/", 1)[-1]).exists()
    # 路径穿越拒绝
    assert client.get("/static/figures/../config.json").status_code in (404, 400)


def test_probe_without_key_is_structured_and_persisted(monkeypatch, tmp_path):
    monkeypatch.setenv("PHYTOREASON_DATA_DIR", str(tmp_path / "user-data"))
    monkeypatch.delenv("OPENAI_API_KEY", raising=False)
    from phyto_reason.config import llm_config
    result = llm_config.probe_config()
    assert result["connectivity"] is False
    assert "error" in result
    persisted = json.loads(llm_config.config_path().read_text(encoding="utf-8"))
    assert "probe" in persisted
    assert "api_key" not in persisted or persisted["api_key"] == ""
