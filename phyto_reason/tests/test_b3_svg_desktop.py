"""Phase B3: SVG 矢量导出、图表列表 API 与桌面上传链补齐测试。"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
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


# ── SVG 双格式导出 ──────────────────────────────────────────

def test_save_figure_writes_png_and_svg_twin(monkeypatch, tmp_path):
    monkeypatch.setenv("PHYTOREASON_DATA_DIR", str(tmp_path / "ud"))
    from phyto_reason.platform_paths import figures_dir
    from phyto_reason.visualization._common import new_figure, save_figure

    fig, ax = new_figure("Volcano Test")
    ax.plot([1, 2, 3], [1, 4, 9])
    url = save_figure(fig, "b3_probe")
    name = url.rsplit("/", 1)[-1]
    assert url.startswith("/static/figures/") and name.endswith(".png")
    assert (figures_dir() / name).exists()
    svg = figures_dir() / name.replace(".png", ".svg")
    assert svg.exists()
    svg_text = svg.read_text(encoding="utf-8", errors="replace")
    assert "<svg" in svg_text[:200]


def test_figure_exporter_png_gets_svg_twin(tmp_path):
    from phyto_reason.visualization.figure_exporter import FigureExporter

    fig, ax = matplotlib.pyplot.subplots()
    ax.plot([0, 1], [0, 1])
    target = FigureExporter.save(fig, str(tmp_path / "volcano_x"), dpi=90)
    matplotlib.pyplot.close(fig)
    assert Path(target).suffix == ".png"
    assert (tmp_path / "volcano_x.svg").exists()


def test_figures_endpoint_lists_png_with_svg_url(client):
    from phyto_reason.visualization._common import new_figure, save_figure

    fig, ax = new_figure("Fig Endpoint Probe")
    ax.plot([1], [1])
    url = save_figure(fig, "b3_endpoint")
    matplotlib.pyplot.close(fig)

    response = client.get("/figures")
    assert response.status_code == 200
    payload = response.json()
    names = [item["name"] for item in payload["figures"]]
    probe = [item for item in payload["figures"] if item["png"] == url]
    assert probe, f"{url} not listed in {names[:5]}"
    assert probe[0]["svg"] and probe[0]["svg"].endswith(".svg")
    assert Path(payload["directory"]).is_absolute()


# ── 桌面端上传链 ────────────────────────────────────────────

pytest.importorskip("PySide6")


def test_desktop_upload_worker_parses_long_table(tmp_path):
    from phyto_reason.desktop.client import UploadWorker

    path = tmp_path / "long_run.csv"
    path.write_text(
        "sample,m/z,RT,intensity\nF-1,348.1,3.2,100\nL-1,348.1,3.2,150\nS-1,200.5,1.0,80\n",
        encoding="utf-8",
    )
    results = []
    worker = UploadWorker(str(path), "metabolite_long")
    worker.result.connect(results.append)
    worker.execute()
    assert results
    payload = results[0]
    assert payload["kind"] == "metabolite"
    assert "长表" in payload["detection_reason"]
    matrix = payload["payload"]
    assert matrix and all(len(row) == 3 for row in matrix.values())


def test_desktop_upload_worker_parses_grouping_table(tmp_path):
    from phyto_reason.desktop.client import UploadWorker

    path = tmp_path / "groups.csv"
    path.write_text(
        "sample_id,tissue,condition\nF-1,flower,Flower\nL-1,leaf,Leaf\n",
        encoding="utf-8",
    )
    results = []
    worker = UploadWorker(str(path), "metadata")
    worker.result.connect(results.append)
    worker.execute()
    assert results
    payload = results[0]
    assert payload["kind"] == "metadata"
    assert payload["payload"]["sample_ids"] == ["F-1", "L-1"]
    assert "tissue" in payload["payload"]["columns_present"]


def test_metadata_only_upload_preserves_existing_expression():
    from phyto_reason.agent.conversation_store import SessionData

    session = SessionData(session_id="b3")
    session.set_data(expression={"GENE1": {"F-1": 1.0}})
    session.set_data(sample_metadata={"sample_ids": ["F-1"], "samples": []})
    assert session.expression_matrix == {"GENE1": {"F-1": 1.0}}
    assert session.sample_metadata == {"sample_ids": ["F-1"], "samples": []}


# ── 前端接线断言 ────────────────────────────────────────────

def test_frontend_figure_panel_wired():
    import phyto_reason.api.app as app_module

    static_dir = Path(app_module.__file__).resolve().parent / "static"
    html = (static_dir / "index.html").read_text(encoding="utf-8")
    assert 'id="figure-list"' in html
    assert 'id="figures-empty"' in html

    js = (static_dir / "assets" / "app.js").read_text(encoding="utf-8")
    assert "async function refreshFigures()" in js
    assert 'fetch("/figures" + query)' in js
    assert 'link.download = figure.name + ".svg";' in js
    css = (static_dir / "assets" / "app.css").read_text(encoding="utf-8")
    assert ".figure-card" in css
