"""
test_api_chat.py — API chat 端点测试（mock LLMClient，无网络）。

覆盖:
  1. POST /chat 非流式（基础问答 + species 注入）
  2. POST /chat/stream SSE 事件序列（stream_start → chunk → done）
  3. GET /session/{id}
"""

from __future__ import annotations

import json
from unittest.mock import MagicMock, patch

import pytest

# 必须在导入 api.app 前 patch LLMClient（模块级 _orchestrator 创建于 import 时）
_llm_mock = MagicMock()
_llm_mock.is_available.return_value = True
_llm_mock.chat_with_tools.return_value = ("Mock response.", [])
_llm_mock.chat_stream.return_value = iter(["Mock ", "response."])

with patch("phyto_reason.agent.orchestrator.LLMClient", return_value=_llm_mock):
    from fastapi.testclient import TestClient
    from phyto_reason.api.app import app

client = TestClient(app)


def test_health():
    r = client.get("/health")
    assert r.status_code == 200


def test_chat_non_streaming():
    r = client.post("/chat", json={"message": "测试", "species": "黄连", "target_metabolite": "berberine"})
    assert r.status_code == 200
    data = r.json()
    assert "response" in data
    assert data["session_id"]
    assert data["error"] is None


def test_chat_stream_sse_sequence():
    r = client.post("/chat/stream", json={"message": "测试流式"})
    assert r.status_code == 200
    assert r.headers["content-type"].startswith("text/event-stream")

    events = []
    for line in r.text.splitlines():
        if line.startswith("data: "):
            events.append(json.loads(line[6:]))

    types = [e["type"] for e in events]
    assert types[0] == "stream_start"
    assert "chunk" in types
    assert types[-1] == "done"
    # chunk 内容拼接为 mock 的流式输出
    chunks = "".join(e.get("content", "") for e in events if e["type"] == "chunk")
    assert chunks == "Mock response."


def test_chat_stream_with_species():
    r = client.post("/chat/stream", json={
        "message": "黄芩素", "species": "黄芩", "target_metabolite": "baicalein",
    })
    assert r.status_code == 200
    # 物种注入不报错，事件正常结束
    assert '"type": "done"' in r.text


def test_get_session():
    r = client.post("/chat", json={"message": "创建会话"})
    sid = r.json()["session_id"]
    r2 = client.get(f"/session/{sid}")
    assert r2.status_code == 200


def test_upload_metadata_persists_alignment_and_groups():
    matrix = b"Name,F-1,F-2,S-1,S-2\nM1,1,2,3,4\n"
    design = b"sample_id,condition,tissue\nF-1,Flower,flower\nF-2,Flower,flower\nS-1,Stem,stem\nS-2,Stem,stem\n"
    response = client.post(
        "/upload",
        files={
            "file": ("meta.csv", matrix, "text/csv"),
            "metadata_file": ("design.csv", design, "text/csv"),
        },
        data={"session_id": "api-metadata-test"},
    )
    assert response.status_code == 200
    payload = response.json()
    assert payload["summary"]["file_type"] == "metabolite_matrix"
    report = payload["diagnostics"]["sample_alignment_report"]
    assert report["n_common"] == 4
    assert report["groups"] == {"Flower": ["F-1", "F-2"], "Stem": ["S-1", "S-2"]}
