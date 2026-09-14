from __future__ import annotations

from pathlib import Path


APP_JS = (Path(__file__).resolve().parents[1] / "api" / "static" / "assets" / "app.js").read_text(encoding="utf-8")


def test_upload_parses_success_json_after_status_check_and_has_non_json_fallback():
    assert 'if (!response.ok)' in APP_JS
    assert 'const errorPayload = await response.json();' in APP_JS
    assert 'payload = await response.json();' in APP_JS
    assert 'throw new Error("HTTP " + response.status);' in APP_JS


def test_metadata_upload_failure_lands_on_matrix_slot():
    """6.4 Step 2 语义修正：矩阵失败落在矩阵槽位；分组表保持待对齐，不再冒充失败来源。"""
    assert 'const intentSlot = SLOT_BY_INTENT[type] || "expr";' in APP_JS
    assert 'updateSlot(intentSlot, file, null, "上传失败: " + error.message);' in APP_JS
    assert '矩阵上传失败，尚未对齐' in APP_JS


def test_degraded_tools_reset_on_first_event_of_each_stream():
    assert 'if (!state.streamStarted)' in APP_JS
    assert 'state.degradedTools = [];' in APP_JS
    assert 'state.streamStarted = false;' in APP_JS


def test_markdown_renders_escaped_contiguous_unordered_lists():
    assert 'output.push("<ul>"' in APP_JS
    assert 'return "<li>" + escapeHtml(item) + "</li>";' in APP_JS
