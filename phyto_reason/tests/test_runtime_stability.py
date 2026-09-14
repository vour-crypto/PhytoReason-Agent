"""Regression tests for runtime contracts at the agent/API boundary."""

from __future__ import annotations

import tempfile
import os
import gc

from phyto_reason.agent.conversation_store import ConversationStore
from phyto_reason.agent.orchestrator import AgentOrchestrator


class _OneShotLLM:
    last_error = None

    def chat_with_tools(self, messages, tools, temperature=0.3):
        return "stable response", []

    def is_available(self):
        return True


def test_memory_store_uses_one_connection_for_in_memory_database():
    store = ConversationStore(":memory:")
    session = store.get_or_create()
    session.add_turn("user", "hello")

    loaded = store.get(session.session_id)
    assert loaded is not None
    assert loaded.history[-1]["content"] == "hello"


def test_new_sessions_are_unique_even_when_created_back_to_back():
    store = ConversationStore(":memory:")
    first = store.get_or_create()
    second = store.get_or_create()
    assert first.session_id != second.session_id


def test_stream_uses_one_final_model_response_and_persists_chain():
    fd, path = tempfile.mkstemp(suffix=".db")
    os.close(fd)
    try:
        store = ConversationStore(path)
        agent = AgentOrchestrator()
        agent._store = store
        agent.llm = _OneShotLLM()

        events = list(agent.handle_stream("hello"))
        assert [event["type"] for event in events] == ["chunk", "done"]
        assert events[0]["content"] == "stable response"

        session = store.get(events[-1]["session_id"])
        assert session is not None
        assert session.message_chain[-1]["content"] == "stable response"
    finally:
        del agent, store
        gc.collect()
        try:
            os.unlink(path)
        except FileNotFoundError:
            pass


def test_session_is_data_backed_by_any_uploaded_matrix():
    store = ConversationStore(":memory:")
    session = store.get_or_create("data-only")
    session.set_data(promoter={"gene1": "ATGC"})
    assert session.has_data is True
