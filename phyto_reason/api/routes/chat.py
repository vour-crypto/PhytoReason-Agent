"""
chat.py — Chat endpoint.

POST /chat            — conversational (AgentOrchestrator, Layer 0-3)
POST /chat/stream     — SSE streaming（AgentOrchestrator.handle_stream 事件）
GET  /session/{id}    — session info
"""

from __future__ import annotations

import json

from fastapi import APIRouter
from fastapi.responses import StreamingResponse
from pydantic import BaseModel, Field

from phyto_reason.agent.orchestrator import AgentOrchestrator

router = APIRouter(tags=["chat"])

_orchestrator = AgentOrchestrator()


def _sse(event: dict) -> str:
    """SSE 帧。"""
    return f"data: {json.dumps(event, ensure_ascii=False)}\n\n"


class ChatRequest(BaseModel):
    message: str = Field(..., description="User message")
    session_id: str = Field(default="", description="Session ID for multi-turn")
    species: str = Field(default="", description="会话物种（学名/俗名/slug，注入 profile 或降级声明）")
    target_metabolite: str = Field(default="", description="目标代谢物")


class ChatResponse(BaseModel):
    response: str
    session_id: str
    tools_used: list[str] = Field(default_factory=list)
    error: str | None = None


@router.post("/chat")
async def chat(request: ChatRequest):
    """Conversational endpoint (AgentOrchestrator, 物种感知)."""
    sid = request.session_id or ""
    if request.species or request.target_metabolite:
        sid = _orchestrator.upload_data(
            sid,
            species=request.species or "",
            target_metabolite=request.target_metabolite or "",
        )
    result = _orchestrator.handle(
        request.message,
        session_id=sid or None,
    )
    return ChatResponse(
        response=result["response"],
        session_id=result["session_id"],
        tools_used=result["tool_calls_made"],
        error=result["error"],
    )


@router.post("/chat/stream")
async def chat_stream(request: ChatRequest):
    """SSE 流式对话（AgentOrchestrator.handle_stream）。

    Events:
        data: {"type":"stream_start"}
        data: {"type":"tool_start","tool":"..."}
        data: {"type":"tool_end","tool":"...","ok":true/false}
        data: {"type":"chunk","content":"..."}
        data: {"type":"error","message":"..."}     — 异常时
        data: {"type":"done","session_id":"..."}   — 结束
    """

    async def generate():
        sid = request.session_id or ""
        if request.species or request.target_metabolite:
            sid = _orchestrator.upload_data(
                sid,
                species=request.species or "",
                target_metabolite=request.target_metabolite or "",
            )
        done_sent = False
        yield _sse({"type": "stream_start"})
        try:
            for event in _orchestrator.handle_stream(request.message, session_id=sid or None):
                done_sent = done_sent or event.get("type") == "done"
                yield _sse(event)
        except Exception as e:
            yield _sse({"type": "error", "message": str(e)})
        finally:
            if not done_sent:
                yield _sse({"type": "done", "session_id": sid})

    return StreamingResponse(
        generate(),
        media_type="text/event-stream",
        headers={
            "Cache-Control": "no-cache",
            "Connection": "keep-alive",
            "X-Accel-Buffering": "no",
            "Access-Control-Allow-Origin": "*",
        },
    )


@router.get("/session/{session_id}")
async def get_session(session_id: str):
    """Get session info."""
    info = _orchestrator.get_session_info(session_id)
    return info


@router.get("/sessions")
async def list_sessions():
    """会话列表（新→旧），供前端会话栏渲染。"""
    return {"sessions": _orchestrator._store.list_sessions()}


@router.post("/sessions")
async def create_session():
    """显式新建空白会话。"""
    session = _orchestrator._store.get_or_create(None)
    return {"session_id": session.session_id}


@router.delete("/sessions/{session_id}")
async def delete_session(session_id: str):
    """删除单个会话（内存 + 数据库）。"""
    _orchestrator._store.delete(session_id)
    return {"deleted": session_id}


class ClearEmptyRequest(BaseModel):
    keep: str = Field(default="", description="保留的会话 ID（当前会话不被清除）")


@router.post("/sessions/clear-empty")
async def clear_empty_sessions(request: ClearEmptyRequest):
    """批量清除无数据的空会话（历史残留），保留当前会话与所有有数据的会话。"""
    empty_ids = [
        s["session_id"] for s in _orchestrator._store.list_sessions()
        if not s.get("has_data") and s["session_id"] != request.keep
    ]
    for sid in empty_ids:
        _orchestrator._store.delete(sid)
    return {"deleted_count": len(empty_ids), "deleted": empty_ids}


@router.get("/hypotheses")
async def list_hypotheses(session_id: str = ""):
    """会话的候选机制假设（假设池数据源）。"""
    if not session_id:
        return {"hypotheses": []}
    session = _orchestrator._store.get(session_id)
    raw = getattr(session, "pipeline_hypotheses", []) or [] if session else []
    hypotheses = []
    for item in raw:
        if hasattr(item, "model_dump"):
            hypotheses.append(item.model_dump())
        elif isinstance(item, dict):
            hypotheses.append(item)
    return {"session_id": session_id, "hypotheses": hypotheses}


@router.delete("/sessions")
async def delete_all_sessions():
    """清空全部会话（含有数据的会话）。破坏性操作，前端需二次确认。"""
    all_ids = [s["session_id"] for s in _orchestrator._store.list_sessions()]
    for sid in all_ids:
        _orchestrator._store.delete(sid)
    return {"deleted_count": len(all_ids)}
