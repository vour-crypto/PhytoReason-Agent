"""用户实测反馈修复回归测试：count 归一化、图库会话隔离、导出报告、图修复。"""

from __future__ import annotations

import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import numpy as np
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


# ── count 文库大小归一化 ────────────────────────────────────

def test_size_factors_track_library_difference():
    from phyto_reason.utils.stats_utils import compute_size_factors, compute_size_factors_array

    counts = {"g1": {"S1": 100.0, "S2": 200.0}, "g2": {"S1": 50.0, "S2": 100.0},
              "g3": {"S1": 80.0, "S2": 160.0}}
    factors = compute_size_factors(counts)
    assert abs(factors["S2"] / factors["S1"] - 2.0) < 1e-9

    array_factors = compute_size_factors_array(np.array([[100.0, 200.0], [50.0, 100.0], [80.0, 160.0]]))
    assert abs(array_factors[1] / array_factors[0] - 2.0) < 1e-9

    # 含零值基因被排除在估计外，但不崩溃
    counts["zero_gene"] = {"S1": 0.0, "S2": 5.0}
    assert abs(compute_size_factors(counts)["S2"] / compute_size_factors(counts)["S1"] - 2.0) < 1e-9


def test_deg_node_normalizes_counts_and_reports_it():
    from phyto_reason.models.planner_state import PlannerState
    from phyto_reason.workflows.execution_nodes import deg_analysis_node
    from phyto_reason.workflows.runtime_state import RuntimeState

    # 两个文库差 2 倍的组；高变基因只在组间差异
    groups = {"A": ["A1", "A2", "A3"], "B": ["B1", "B2", "B3"]}
    expr = {}
    rng = np.random.default_rng(7)
    for i in range(60):
        row = {}
        for j, sample in enumerate(["A1", "A2", "A3"]):
            base = 100.0 * (2 ** j)  # 文库因子 1/2/4
            row[sample] = base * (3.0 if i == 0 else 1.0) * (1 + 0.01 * rng.standard_normal())
        for j, sample in enumerate(["B1", "B2", "B3"]):
            base = 200.0 * (2 ** j)
            row[sample] = base * (3.0 if i == 0 else 1.0) * (1 + 0.01 * rng.standard_normal())
        expr[f"G{i}"] = row
    ps = PlannerState(species="test", target_metabolite="", has_expression=True,
                      has_metabolite=False, expression_matrix=expr,
                      sample_metadata={"samples": [
                          {"sample_id": s, "condition": g, "tissue": g}
                          for g, samples in groups.items() for s in samples
                      ]})
    state = RuntimeState(planner_state=ps)
    result = deg_analysis_node(state)
    report = result.get("deg_report", {})
    assert report.get("normalization") == "median_of_ratios_size_factors + log2(count/size_factor + 1)"


# ── 图修复 ─────────────────────────────────────────────────

def test_scale_free_plot_uses_r2_key(monkeypatch, tmp_path):
    monkeypatch.setenv("PHYTOREASON_DATA_DIR", str(tmp_path / "ud"))
    from phyto_reason.visualization.multiomics_figures import plot_wgcna_scale_free

    power_results = [{"power": p, "r2": 0.1 * p} for p in (2, 4, 6, 8)]
    url = plot_wgcna_scale_free(power_results, best_power=8)
    assert url.endswith(".png")
    name = url.rsplit("/", 1)[-1]
    from phyto_reason.platform_paths import figures_dir
    assert (figures_dir() / name).exists()
    assert (figures_dir() / name.replace(".png", ".svg")).exists()


def test_qc_heatmap_hides_labels_for_many_genes():
    from phyto_reason.visualization.heatmap import HeatmapPlotter

    rng = np.random.default_rng(3)
    matrix = {f"G{i}": {f"S{j}": float(rng.normal(i, 1)) for j in range(5)} for i in range(80)}
    fig = HeatmapPlotter.clustered_heatmap(matrix, list(matrix), title="many genes")
    assert fig.axes[0].get_yticklabels() == []
    matplotlib.pyplot.close(fig)

    small = {k: matrix[k] for k in list(matrix)[:10]}
    fig2 = HeatmapPlotter.clustered_heatmap(small, list(small), title="few genes")
    assert len(fig2.axes[0].get_yticklabels()) == 10
    matplotlib.pyplot.close(fig2)


# ── 图库会话隔离 ───────────────────────────────────────────

def test_figures_endpoint_since_filter(client):
    from phyto_reason.visualization._common import new_figure, save_figure

    fig, ax = new_figure("Since Probe")
    ax.plot([1], [1])
    url = save_figure(fig, "since_probe")
    matplotlib.pyplot.close(fig)

    past = client.get("/figures", params={"since": time.time() - 3600})
    assert past.status_code == 200
    assert any(item["png"] == url for item in past.json()["figures"])

    future = client.get("/figures", params={"since": time.time() + 3600})
    assert future.status_code == 200
    assert all(item["png"] != url for item in future.json()["figures"])


def test_frontend_gallery_scoping_and_export_wired():
    import phyto_reason.api.app as app_module

    static_dir = Path(app_module.__file__).resolve().parent / "static"
    html = (static_dir / "index.html").read_text(encoding="utf-8")
    assert 'id="gallery-history"' in html and 'id="export-report"' in html

    js = (static_dir / "assets" / "app.js").read_text(encoding="utf-8")
    assert "gallerySince" in js and "galleryIncludeHistory" in js
    assert '"?since=" + state.gallerySince' in js
    assert "async function exportReport()" in js
    assert "phytoreason_report_" in js
