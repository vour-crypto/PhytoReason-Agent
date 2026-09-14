"""
orchestrator.py — LLM-first Agent Orchestrator (v4.0).

Single entry point for all user interactions.
LLM drives the conversation, decides which tools to call, and interprets results.

Architecture:
    User message → LLM (with tool definitions) → tool calls? → execute → LLM → response

Supports 4 operating layers:
    Layer 0: Zero-data — literature + knowledge base + KEGG
    Layer 1: Cross-species — ortholog inference without user data
    Layer 2: Single-omics — DEG/DAM analysis with public knowledge
    Layer 3: Full multi-omics — complete 10-node scientific pipeline
"""

from __future__ import annotations

import json
import logging
import re
from typing import Any

from phyto_reason.agent.system_prompt import build_system_prompt
from phyto_reason.agent.tool_definitions import TOOL_DEFINITIONS
from phyto_reason.agent.tool_handlers import (
    handle_search_literature,
    handle_query_kegg,
    handle_query_plantcyc,
    handle_annotate_ms2_spectrum,
    handle_run_pipeline,
    handle_check_quality,
    handle_check_data_quality,
    handle_run_deg_analysis,
    handle_run_dam_analysis,
    handle_run_multiomics,
    handle_run_wgcna,
    handle_run_tf_analysis,
    handle_run_hypothesis_synthesis,
    handle_search_knowledge_base,
    handle_cross_species_infer,
    handle_search_rag,
    handle_web_search,
    handle_recall_memory,
    handle_query_public_expression,
)
from phyto_reason.agent.conversation_store import ConversationStore, SessionData, store
from phyto_reason.llm.llm_client import LLMClient

logger = logging.getLogger("orchestrator")

MAX_TOOL_ROUNDS = 5  # max tool-call iterations per user message


def _detect_chinese(text: str) -> bool:
    """Detect if text contains Chinese (CJK) characters."""
    return bool(re.search(r'[一-鿿㐀-䶿]', text))


# Language enforcement instruction injected into user messages
_LANG_INSTRUCTION_CN = (
    "[LANGUAGE: The user is writing in Chinese. You MUST respond ENTIRELY in Chinese (中文). "
    "Every header, every sentence, every label MUST be in Chinese. "
    "Only gene names, database IDs, pathway codes, and species Latin names stay in original form. "
    "Do NOT write any section heading or body text in English.]"
)

_LANG_INSTRUCTION_EN = (
    "[LANGUAGE: The user is writing in English. You MUST respond in English. "
    "Do not use Chinese in your response.]"
)


class AgentOrchestrator:
    """LLM-first agent orchestrator.

    Usage:
        orchestrator = AgentOrchestrator()
        result = orchestrator.handle("帮我分析小檗碱的调控", session_id="...")
        # result contains the assistant's response
    """

    def __init__(self) -> None:
        self.llm = LLMClient()
        self._store: ConversationStore = store

    # ═══════════════════════════════════════════════════════════
    # Public API
    # ═══════════════════════════════════════════════════════════

    def handle(self, user_message: str, session_id: str | None = None) -> dict[str, Any]:
        """Process a user message and return the assistant's response.

        Args:
            user_message: The user's natural language message.
            session_id: Optional session ID for multi-turn conversations.

        Returns:
            dict with keys:
                - "response": str (the assistant's text response)
                - "session_id": str
                - "tool_calls_made": list[str] (names of tools called)
                - "error": str | None
        """
        refresh = getattr(self.llm, "refresh_if_changed", None)
        if refresh:
            refresh()
        session = self._store.get_or_create(session_id)

        try:
            # 1. Build messages for LLM
            messages = self._build_messages(user_message, session)

            # 2. LLM + tool loop
            tool_calls_made: list[str] = []
            response_text = ""

            for round_num in range(MAX_TOOL_ROUNDS):
                content, tool_calls = self.llm.chat_with_tools(
                    messages=messages,
                    tools=TOOL_DEFINITIONS,
                    temperature=0.3,
                )

                # No tool calls → LLM is done, this is the final response
                if not tool_calls:
                    llm_error = getattr(self.llm, "last_error", None)
                    if not content and isinstance(llm_error, Exception):
                        raise llm_error
                    response_text = content
                    break

                # Add assistant message (with tool calls) to history
                # Arguments must be JSON strings for the API
                serialized_calls = []
                for tc in tool_calls:
                    sc = dict(tc)
                    fn_args = sc.get("function", {}).get("arguments", {})
                    if isinstance(fn_args, dict):
                        sc["function"]["arguments"] = json.dumps(fn_args, ensure_ascii=False)
                    serialized_calls.append(sc)

                assistant_msg: dict = {"role": "assistant", "content": content or ""}
                assistant_msg["tool_calls"] = serialized_calls
                messages.append(assistant_msg)

                # Execute each tool call
                for tc in tool_calls:
                    fn_name = tc.get("function", {}).get("name", "")
                    fn_args = tc.get("function", {}).get("arguments", {})

                    if isinstance(fn_args, str):
                        try:
                            fn_args = json.loads(fn_args)
                        except json.JSONDecodeError:
                            fn_args = {}

                    tool_result = self._execute_tool(fn_name, fn_args, session)

                    # If tool returned an error, inject degradation note
                    if isinstance(tool_result, str) and tool_result.startswith("Error"):
                        tool_result = f"[DEGRADED] {tool_result}\nThe analysis continues with reduced confidence for this evidence type."

                    tool_calls_made.append(fn_name)

                    # Add tool result to messages
                    messages.append({
                        "role": "tool",
                        "tool_call_id": tc.get("id", ""),
                        "content": tool_result,
                    })

                # Safety: prevent infinite loops
                if round_num >= MAX_TOOL_ROUNDS - 1:
                    # Force LLM to synthesize without tools
                    messages.append({
                        "role": "user",
                        "content": "Please synthesize all the above tool results into a comprehensive answer. "
                                   "Do NOT call any more tools. Give me the final analysis now."
                    })
                    final_content, _ = self.llm.chat_with_tools(
                        messages=messages,
                        tools=[],  # no tools allowed
                        temperature=0.3,
                    )
                    response_text = final_content or "Analysis complete."
                    break

            # 3. If no response yet (e.g., all tool calls, no final message)
            if not response_text:
                messages.append({
                    "role": "user",
                    "content": "Summarize the findings above concisely."
                })
                final_content, _ = self.llm.chat_with_tools(
                    messages=messages,
                    tools=[],
                    temperature=0.3,
                )
                response_text = final_content or "Analysis complete. Please ask if you have questions."

            # 4. Record turns + auto-capture memory
            session.add_turn("user", user_message)
            session.add_turn("assistant", response_text)

            # Save full message chain for context continuity across turns
            # Include all tool calls + results so follow-up questions have full context
            # Strip the system prompt (index 0) to avoid duplication in next turn
            chain = [msg for msg in messages if msg.get("role") != "system"]
            # Add the final assistant response to the chain
            chain.append({"role": "assistant", "content": response_text})
            # Limit chain size (keep last ~30 messages to stay within context limits)
            if len(chain) > 40:
                chain = chain[-40:]
            session.message_chain = chain
            session._persist()

            self._capture_turn_memory(user_message, response_text, session)

            return {
                "response": response_text,
                "session_id": session.session_id,
                "tool_calls_made": tool_calls_made,
                "error": None,
            }

        except Exception as e:
            logger.error(f"Orchestrator error: {e}", exc_info=True)

            # Distinguish LLM failures from other errors
            from phyto_reason.llm.llm_client import LLMUnavailableError
            if isinstance(e, LLMUnavailableError) or (
                hasattr(self, "llm") and not self.llm.is_available()
            ):
                fallback = (
                    "I'm currently unable to reach the language model. "
                    "You can still:\n"
                    "1. Upload data and I'll check its quality offline\n"
                    "2. Browse the built-in knowledge base for TF-metabolite relationships\n"
                    "3. Check KEGG pathway information\n\n"
                    "Please try again in a moment, or check your API configuration."
                )
            else:
                fallback = (
                    "I encountered an error processing your request. "
                    "Please try again or rephrase your question."
                )

            session.add_turn("user", user_message)
            session.add_turn("assistant", fallback)
            return {
                "response": fallback,
                "session_id": session.session_id,
                "tool_calls_made": [],
                "error": str(e),
            }

    def handle_streamlit(self, user_message: str, session_id: str | None = None) -> str:
        """Convenience method for Streamlit: returns just the response text."""
        result = self.handle(user_message, session_id)
        return result["response"]

    def handle_stream(self, user_message: str, session_id: str | None = None):
        """Streaming version of handle() — yields SSE-style events.

        Events:
            {"type":"tool_start", "tool":"search_literature"}
            {"type":"tool_end", "tool":"search_literature", "ok":true}
            {"type":"chunk", "content":"..."}
            {"type":"figure", "url":"...", "caption":"..."}
            {"type":"done", "session_id":"..."}
        """
        refresh = getattr(self.llm, "refresh_if_changed", None)
        if refresh:
            refresh()
        session = self._store.get_or_create(session_id)

        try:
            messages = self._build_messages(user_message, session)
            tool_calls_made: list[str] = []
            response_text = ""

            for round_num in range(MAX_TOOL_ROUNDS):
                content, tool_calls = self.llm.chat_with_tools(
                    messages=messages,
                    tools=TOOL_DEFINITIONS,
                    temperature=0.3,
                )

                if not tool_calls:
                    llm_error = getattr(self.llm, "last_error", None)
                    if not content and isinstance(llm_error, Exception):
                        raise llm_error
                    # The tool-aware call already produced the final answer.
                    response_text = content or ""
                    if response_text:
                        yield {"type": "chunk", "content": response_text}
                    break

                # Serialize tool calls for API
                serialized_calls = []
                for tc in tool_calls:
                    sc = dict(tc)
                    fn_args = sc.get("function", {}).get("arguments", {})
                    if isinstance(fn_args, dict):
                        sc["function"]["arguments"] = json.dumps(fn_args, ensure_ascii=False)
                    serialized_calls.append(sc)

                assistant_msg = {"role": "assistant", "content": content or ""}
                assistant_msg["tool_calls"] = serialized_calls
                messages.append(assistant_msg)

                for tc in tool_calls:
                    fn_name = tc.get("function", {}).get("name", "")
                    fn_args = tc.get("function", {}).get("arguments", {})

                    if isinstance(fn_args, str):
                        try:
                            fn_args = json.loads(fn_args)
                        except json.JSONDecodeError:
                            fn_args = {}

                    # Inject session_id
                    fn_args["session_id"] = session.session_id

                    yield {"type": "tool_start", "tool": fn_name}

                    try:
                        tool_result = self._execute_tool(fn_name, fn_args, session)
                        ok = not (isinstance(tool_result, str) and tool_result.startswith("Error"))
                    except Exception as e:
                        tool_result = f"Error: {e}"
                        ok = False

                    tool_calls_made.append(fn_name)
                    yield {"type": "tool_end", "tool": fn_name, "ok": ok}

                    if isinstance(tool_result, str) and tool_result.startswith("Error"):
                        tool_result = f"[DEGRADED] {tool_result}\nThe analysis continues with reduced confidence."

                    messages.append({
                        "role": "tool",
                        "tool_call_id": tc.get("id", ""),
                        "content": tool_result,
                    })

                if round_num >= MAX_TOOL_ROUNDS - 1:
                    messages.append({
                        "role": "user",
                        "content": "Please synthesize all the above tool results into a comprehensive answer. "
                                   "Do NOT call any more tools. Give me the final analysis now."
                    })
                    final_content, _ = self.llm.chat_with_tools(
                        messages=messages, tools=[], temperature=0.3,
                    )
                    llm_error = getattr(self.llm, "last_error", None)
                    if not final_content and isinstance(llm_error, Exception):
                        raise llm_error
                    response_text = final_content or ""
                    if response_text:
                        yield {"type": "chunk", "content": response_text}
                    break

            # Record turns
            session.add_turn("user", user_message)
            session.add_turn("assistant", response_text or "(empty response)")
            chain = [msg for msg in messages if msg.get("role") != "system"]
            chain.append({"role": "assistant", "content": response_text})
            session.message_chain = chain[-40:]
            session._persist()

            yield {"type": "done", "session_id": session.session_id}

        except Exception as e:
            logger.error(f"Stream error: {e}", exc_info=True)
            yield {"type": "error", "message": str(e)}
            yield {"type": "done", "session_id": session.session_id}

    def upload_data(
        self,
        session_id: str,
        expression: dict | None = None,
        metabolite: dict | None = None,
        promoter: dict | None = None,
        sample_metadata: dict | None = None,
        sample_alignment_report: dict | None = None,
        species: str = "",
        target_metabolite: str = "",
    ) -> str:
        """Upload omics data to a session.

        Returns the session ID (created if necessary).
        """
        session = self._store.get_or_create(session_id)
        effective_metadata = sample_metadata if sample_metadata is not None else session.sample_metadata
        if effective_metadata is not None and sample_alignment_report is None:
            sample_alignment_report = self._build_sample_alignment_report(
                expression or session.expression_matrix,
                metabolite or session.metabolite_matrix,
                effective_metadata,
            )
        session.set_data(
            expression=expression, metabolite=metabolite, promoter=promoter,
            sample_metadata=sample_metadata,
            sample_alignment_report=sample_alignment_report,
        )
        if species:
            session.species = species
        if target_metabolite:
            session.target_metabolite = target_metabolite
        logger.info(f"Data uploaded to session {session.session_id}: "
                    f"expr={bool(expression)}, meta={bool(metabolite)}, "
                    f"species={species}, target={target_metabolite}")
        return session.session_id

    @staticmethod
    def _build_sample_alignment_report(
        expression: dict | None,
        metabolite: dict | None,
        sample_metadata: dict,
    ) -> dict:
        """Build a serializable alignment/design summary for uploaded data."""
        def matrix_samples(matrix: dict | None) -> list[str]:
            if not matrix:
                return []
            first = next(iter(matrix.values()), {})
            return [str(s) for s in first.keys()] if isinstance(first, dict) else []

        def metadata_samples(value: dict) -> list[str]:
            if "sample_ids" in value:
                return [str(s) for s in value.get("sample_ids", [])]
            if "samples" in value:
                return [str(row.get("sample_id")) for row in value.get("samples", []) if row.get("sample_id")]
            return [str(s) for s in value.keys() if s not in {"source_file", "columns_present", "columns_missing"}]

        expr_samples = matrix_samples(expression)
        meta_samples = matrix_samples(metabolite)
        design_samples = metadata_samples(sample_metadata)
        sides = [set(x) for x in (expr_samples, meta_samples, design_samples) if x]
        common = sorted(set.intersection(*sides)) if sides else []
        groups: dict[str, list[str]] = {}
        rows = sample_metadata.get("samples", []) if isinstance(sample_metadata, dict) else []
        for row in rows:
            sid = str(row.get("sample_id", ""))
            group = str(row.get("condition") or row.get("tissue") or row.get("treatment") or "").strip()
            if sid and group:
                groups.setdefault(group, []).append(sid)
        return {
            "n_common": len(common),
            "common_samples": common,
            "dropped_from_expression": sorted(set(expr_samples) - set(common)),
            "dropped_from_metabolite": sorted(set(meta_samples) - set(common)),
            "dropped_from_metadata": sorted(set(design_samples) - set(common)),
            "groups": groups,
            "group_source": "metadata" if groups else "none",
        }

    def get_session_info(self, session_id: str) -> dict:
        """Get session summary for the UI."""
        session = self._store.get(session_id)
        if not session:
            return {"error": "Session not found"}
        return session.to_summary()

    # ═══════════════════════════════════════════════════════════
    # Memory auto-capture
    # ═══════════════════════════════════════════════════════════

    def _capture_turn_memory(self, user_msg: str, response: str, session: SessionData) -> None:
        """Auto-save an observation after each conversation turn.
        Detects user corrections and saves them as high-importance feedback.
        """
        try:
            from phyto_reason.agent.memory_store import memory_store, MemoryEntry

            short = user_msg[:200]
            tags = []
            if session.species:
                tags.append(session.species)
            if session.target_metabolite:
                tags.append(session.target_metabolite)

            # Detect user corrections (negative feedback patterns)
            correction_patterns = [
                r'不对', r'错了', r'应该是', r'其实是', r'不是.*而是',
                r'纠正', r'修正', r'搞错了', r'说错了', r'弄错了',
                r'你错了', r'wrong', r'incorrect', r"that's not",
                r'actually', r'should be', r'correct that',
                r'it is not', r'that is wrong',
            ]
            is_correction = any(
                re.search(pat, user_msg, re.IGNORECASE)
                for pat in correction_patterns
            )

            if is_correction:
                # Extract the factual claim from the correction
                correction_text = short[:300]
                memory_store.add(MemoryEntry(
                    session_id=session.session_id,
                    content=f"USER CORRECTION: {correction_text}",
                    memory_type="observation",
                    importance=9,
                    tags=tags + ["correction", "feedback"],
                    species=session.species,
                    metabolite=session.target_metabolite,
                ))
                logger.info("Correction memory saved: %s", correction_text[:100])
                return

            # Regular observation
            memory_store.add_observation(
                content=f"User asked: {short}",
                session_id=session.session_id,
                importance=3,
                tags=tags,
                species=session.species,
                metabolite=session.target_metabolite,
            )
        except Exception as e:
            logger.warning("Failed to capture turn memory: %s", e)

    def _capture_pipeline_memory(self, result: str, session: SessionData) -> None:
        """Auto-save a result memory after pipeline execution."""
        try:
            from phyto_reason.agent.memory_store import memory_store

            # Extract key findings from the first 500 chars
            summary = result[:500]
            if len(result) > 500:
                summary += "..."

            tags = [session.species, session.target_metabolite]
            if session.target_pathway:
                tags.append(session.target_pathway)
            tags = [t for t in tags if t]

            memory_store.add_result(
                content=summary,
                session_id=session.session_id,
                importance=7,
                tags=tags,
                species=session.species,
                metabolite=session.target_metabolite,
            )
        except Exception as e:
            logger.warning("Failed to capture pipeline memory: %s", e)

    # ═══════════════════════════════════════════════════════════
    # Internal
    # ═══════════════════════════════════════════════════════════

    def _build_messages(self, user_message: str, session: SessionData) -> list[dict]:
        """Build the full message list for the LLM call."""
        messages = []

        # System prompt（物种感知模板：有 profile 直接答，无则降级声明）
        context_summary = session.get_context_summary()
        use_compact = session.turn_count > 10
        system_text = build_system_prompt(
            species=session.species or "",
            target_metabolite=session.target_metabolite or "",
            compact=use_compact,
        )

        messages.append({
            "role": "system",
            "content": system_text + "\n\n## Current Session Context\n" + context_summary,
        })

        # Conversation history (full message chain with tool calls/results)
        # Preserves context across turns so the LLM remembers what tools were called
        if session.message_chain:
            # Use the stored full message chain from the previous turn
            # Skip the system prompt from stored messages (we already have one)
            for msg in session.message_chain:
                if msg.get("role") != "system":
                    messages.append(msg)
        else:
            # Fallback: simple text history
            for turn in session.history[-20:]:
                messages.append(turn)

        # Recall relevant user corrections from memory
        correction_prefix = ""
        try:
            from phyto_reason.agent.memory_store import memory_store

            # Two-tier correction search:
            # 1. Semantic search by query (TF-IDF relevance)
            # 2. Tag-based fallback (always finds corrections regardless of keyword overlap)
            semantic_results = memory_store.search(query=user_message, top_k=2)
            tag_results = memory_store.get_by_tag("correction", n=5)

            # Merge, deduplicate, and keep top matches
            seen_ids: set[int] = set()
            all_corrections: list[str] = []

            for m in semantic_results + tag_results:
                mid = m.get("id")
                if mid is not None and mid in seen_ids:
                    continue
                seen_ids.add(mid)

                content = m.get("content", "")
                if "correction" in m.get("tags", []) or "USER CORRECTION" in content:
                    all_corrections.append(content[:200])

            if all_corrections:
                correction_prefix = (
                    "[PREVIOUS CORRECTIONS — The user has corrected the following. "
                    "Use this information to avoid repeating mistakes:]\n"
                    + "\n".join(f"  - {c}" for c in all_corrections[:3])
                    + "\n\n"
                )
        except Exception as e:
            logger.warning("Failed to recall correction memory: %s", e)

        # Current user message with language enforcement + correction prefix
        if _detect_chinese(user_message):
            messages.append({"role": "user", "content": correction_prefix + f"{_LANG_INSTRUCTION_CN}\n\n{user_message}"})
        elif _detect_chinese("".join(
            t.get("content", "") for t in session.history[-4:]
            if t.get("role") == "user"
        )):
            messages.append({"role": "user", "content": correction_prefix + f"{_LANG_INSTRUCTION_CN}\n\n{user_message}"})
        else:
            messages.append({"role": "user", "content": correction_prefix + user_message})

        return messages

    def _execute_tool(self, tool_name: str, arguments: dict, session: SessionData) -> str:
        """Execute a tool and return the result string."""
        logger.info(f"Executing tool: {tool_name} with args: {arguments}")

        try:
            if tool_name == "search_literature":
                return handle_search_literature(arguments)

            elif tool_name == "query_kegg":
                return handle_query_kegg(arguments)

            elif tool_name == "query_plantcyc":
                return handle_query_plantcyc(arguments)

            elif tool_name == "annotate_ms2_spectrum":
                return handle_annotate_ms2_spectrum(arguments)

            elif tool_name == "run_scientific_pipeline":
                # Merge context from session into arguments
                if not arguments.get("species") and session.species:
                    arguments["species"] = session.species
                if not arguments.get("target_metabolite") and session.target_metabolite:
                    arguments["target_metabolite"] = session.target_metabolite

                result = handle_run_pipeline(
                    arguments,
                    expression_matrix=session.expression_matrix,
                    metabolite_matrix=session.metabolite_matrix,
                    promoter_sequences=session.promoter_sequences,
                    species=session.species,
                    target_metabolite=session.target_metabolite,
                    session=session,
                )

                # Cache result in session for follow-up questions
                session.set_pipeline_result(result)
                session.update_context(
                    species=arguments.get("species", ""),
                    metabolite=arguments.get("target_metabolite", ""),
                    pathway=arguments.get("target_pathway", ""),
                )

                # Auto-capture pipeline result to memory stream
                if not result.startswith("Error") and not result.startswith("Pipeline error"):
                    self._capture_pipeline_memory(result, session)

                return result

            elif tool_name == "check_data_quality":
                return handle_check_data_quality(arguments, session=session)

            elif tool_name == "run_deg_analysis":
                return handle_run_deg_analysis(arguments, session=session)

            elif tool_name == "run_dam_analysis":
                return handle_run_dam_analysis(arguments, session=session)

            elif tool_name == "run_multiomics":
                return handle_run_multiomics(arguments, session=session)

            elif tool_name == "run_wgcna":
                return handle_run_wgcna(arguments, session=session)

            elif tool_name == "run_tf_analysis":
                return handle_run_tf_analysis(arguments, session=session)

            elif tool_name == "run_hypothesis_synthesis":
                return handle_run_hypothesis_synthesis(arguments, session=session)

            elif tool_name == "search_knowledge_base":
                return handle_search_knowledge_base(arguments)

            elif tool_name == "cross_species_infer":
                return handle_cross_species_infer(arguments)

            elif tool_name == "search_rag":
                return handle_search_rag(arguments)

            elif tool_name == "web_search":
                return handle_web_search(arguments)

            elif tool_name == "recall_memory":
                return handle_recall_memory(arguments)
            elif tool_name == "query_public_expression":
                return handle_query_public_expression(arguments)

            elif tool_name in {"massbank_spectrum_search", "mona_spectrum_search", "pubchem_compound_properties"}:
                from phyto_reason.tools import TOOL_REGISTRY
                tool = TOOL_REGISTRY.get(tool_name)
                if tool is None:
                    return f"Unknown tool: {tool_name}"
                result = tool.run(**arguments)
                return json.dumps({
                    "status": result.status,
                    "warnings": result.warnings,
                    "metadata": result.metadata,
                }, ensure_ascii=False, default=str)

            else:
                return f"Unknown tool: {tool_name}"

        except Exception as e:
            logger.error(f"Tool execution failed: {tool_name}: {e}", exc_info=True)
            return f"Error executing {tool_name}: {e}"
