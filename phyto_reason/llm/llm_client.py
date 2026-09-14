"""
llm_client.py — DeepSeek / OpenAI-compatible LLM client (v4.0).

从 .env 读取配置:
  OPENAI_API_KEY
  OPENAI_BASE_URL
  LLM_MODEL

Features:
  - Exponential backoff retry (3 attempts: 1s → 2s → 4s)
  - Graceful degradation when API is unavailable
"""

from __future__ import annotations

import json
import logging
import os
import time
from typing import Any

from openai import OpenAI
from openai import (
    APIError,
    APIConnectionError,
    RateLimitError,
    APITimeoutError,
)
from dotenv import load_dotenv

logger = logging.getLogger("llm_client")

load_dotenv(os.path.join(os.path.dirname(__file__), "..", "..", ".env"))

# Retry config
MAX_RETRIES = 3
BASE_DELAY = 1.0  # seconds, doubles each retry


class LLMUnavailableError(Exception):
    """Raised when the LLM is completely unavailable after all retries."""
    pass


class LLMClient:
    """LLM 客户端 (OpenAI-compatible / DeepSeek)。"""

    def __init__(self, model: str | None = None) -> None:
        self._config_fingerprint = None
        self.api_key = ""
        self.base_url = ""
        self.model = ""
        self._load_current_config(model)
        self._client: OpenAI | None = None
        self.last_error: Exception | None = None

    def _load_current_config(self, model: str | None = None) -> None:
        from phyto_reason.config.llm_config import config_fingerprint, load_config
        value = load_config()
        self.api_key = value.get("api_key", "")
        self.base_url = value.get("base_url", "https://api.deepseek.com/v1")
        self.model = model or value.get("model", "deepseek-chat")
        self._config_fingerprint = config_fingerprint()

    def refresh_if_changed(self) -> None:
        from phyto_reason.config.llm_config import config_fingerprint
        fingerprint = config_fingerprint()
        if fingerprint != self._config_fingerprint:
            self._load_current_config()
            self._client = None

    @property
    def client(self) -> OpenAI:
        if self._client is None:
            headers = {"Authorization": f"Bearer {self.api_key}"}
            self._client = OpenAI(
                api_key=self.api_key,
                base_url=self.base_url,
                default_headers=headers,
            )
        return self._client

    # ═══════════════════════════════════════════════════════════
    # Retry wrapper
    # ═══════════════════════════════════════════════════════════

    def _call_with_retry(self, fn, *args, **kwargs):
        """Call a function with exponential backoff retry.

        Retries on: ConnectionError, Timeout, RateLimitError, APIConnectionError.
        Does NOT retry on: AuthenticationError, BadRequestError (user error).
        """
        last_error = None
        for attempt in range(MAX_RETRIES + 1):
            try:
                return fn(*args, **kwargs)
            except (APIConnectionError, APITimeoutError, ConnectionError, OSError) as e:
                last_error = e
                if attempt < MAX_RETRIES:
                    delay = BASE_DELAY * (2 ** attempt)
                    logger.warning(
                        f"LLM connection error (attempt {attempt + 1}/{MAX_RETRIES + 1}), "
                        f"retrying in {delay:.1f}s: {e}"
                    )
                    time.sleep(delay)
                else:
                    logger.error(f"LLM unavailable after {MAX_RETRIES + 1} attempts: {e}")
            except RateLimitError as e:
                last_error = e
                if attempt < MAX_RETRIES:
                    delay = BASE_DELAY * (3 ** attempt)  # longer backoff for rate limits
                    logger.warning(
                        f"LLM rate limited (attempt {attempt + 1}/{MAX_RETRIES + 1}), "
                        f"waiting {delay:.1f}s: {e}"
                    )
                    time.sleep(delay)
                else:
                    logger.error(f"LLM rate limited after {MAX_RETRIES + 1} attempts: {e}")
            except APIError as e:
                # Server errors (5xx) → retry; client errors (4xx) → don't
                status = getattr(e, 'status_code', 0) or getattr(e, 'http_status', 0)
                if 500 <= status < 600 and attempt < MAX_RETRIES:
                    delay = BASE_DELAY * (2 ** attempt)
                    logger.warning(f"LLM server error {status} (attempt {attempt + 1}), retrying in {delay:.1f}s")
                    time.sleep(delay)
                    last_error = e
                else:
                    raise  # client error or auth error — don't retry

        raise LLMUnavailableError(
            f"LLM unavailable after {MAX_RETRIES + 1} attempts. "
            f"Last error: {last_error}"
        )

    # ═══════════════════════════════════════════════════════════
    # Public API
    # ═══════════════════════════════════════════════════════════

    def chat(self, messages: list[dict], temperature: float = 0.3,
             max_tokens: int = 4096, **kwargs) -> str:
        """普通对话补全（带重试）。"""

        def _call():
            resp = self.client.chat.completions.create(
                model=self.model,
                messages=messages,
                temperature=temperature,
                max_tokens=max_tokens,
                **kwargs,
            )
            return resp.choices[0].message.content or ""

        self.last_error = None
        try:
            return self._call_with_retry(_call)
        except Exception as e:
            self.last_error = e
            logger.error(f"chat() failed: {e}")
            return ""

    def chat_structured(self, messages: list[dict], response_format: dict,
                        temperature: float = 0.1, **kwargs) -> dict:
        """结构化输出 (JSON mode, 带重试)。

        response_format 归一化：调用方可传裸 JSON Schema（含 "properties"），
        会自动归一化为 API 可接受的 ``{"type": "json_object"}``（DeepSeek 不认
        裸 schema；json_object 模式下由提示词约束结构，返回后由调用方校验）。
        """

        def _normalize(fmt: dict) -> dict:
            if fmt.get("type") in ("json_object", "json_schema", "regex", "text"):
                return fmt
            return {"type": "json_object"}  # 裸 schema → json_object 模式

        def _call():
            resp = self.client.chat.completions.create(
                model=self.model,
                messages=messages,
                temperature=temperature,
                response_format=_normalize(response_format),
                **kwargs,
            )
            content = resp.choices[0].message.content or "{}"
            try:
                return json.loads(content)
            except json.JSONDecodeError:
                return {"raw": content, "parse_error": True}

        self.last_error = None
        try:
            return self._call_with_retry(_call)
        except Exception as e:
            self.last_error = e
            logger.error(f"chat_structured() failed: {e}")
            return {"error": str(e), "parse_error": True}

    def chat_with_tools(self, messages: list[dict], tools: list[dict],
                        temperature: float = 0.1, **kwargs) -> tuple[str, list[dict]]:
        """Function calling — 返回 (text_content, tool_calls)。带重试。"""

        def _call():
            resp = self.client.chat.completions.create(
                model=self.model,
                messages=messages,
                tools=tools,
                temperature=temperature,
                **kwargs,
            )
            msg = resp.choices[0].message
            content = msg.content or ""
            calls = []
            if msg.tool_calls:
                for tc in msg.tool_calls:
                    args = {}
                    try:
                        args = json.loads(tc.function.arguments)
                    except json.JSONDecodeError:
                        args = {"_raw": tc.function.arguments}
                    calls.append({
                        "id": tc.id,
                        "type": tc.type,
                        "function": {"name": tc.function.name, "arguments": args},
                    })
            return content, calls

        self.last_error = None
        try:
            return self._call_with_retry(_call)
        except Exception as e:
            self.last_error = e
            logger.error(f"chat_with_tools() failed: {e}")
            return "", []

    def chat_stream(self, messages: list[dict], temperature: float = 0.3,
                    max_tokens: int = 4096, emit_reasoning: bool = False,
                    **kwargs):
        """流式对话补全 — 逐 chunk yield 文本。

        Args:
            messages: 对话消息列表。
            temperature: 采样温度。
            max_tokens: 最大生成长度。
            emit_reasoning: 若 True，思考模式 (deepseek-reasoner) 的
                reasoning_content 也会被 yield，并用标记包裹，便于调用方区分。
                yield 格式:
                  {"type": "reasoning", "content": "..."}   — 思考过程
                  {"type": "content", "content": "..."}     — 最终回答
                若 emit_reasoning=False (默认)，只 yield 最终回答文本
                (与旧行为一致，向后兼容)。
            **kwargs: 透传给 API 的额外参数。

        Yields:
            见 emit_reasoning 说明。
        """
        try:
            stream = self.client.chat.completions.create(
                model=self.model,
                messages=messages,
                temperature=temperature,
                max_tokens=max_tokens,
                stream=True,
                **kwargs,
            )
            for chunk in stream:
                delta = chunk.choices[0].delta
                # DeepSeek reasoner: reasoning_content 携带思考链
                reasoning = getattr(delta, "reasoning_content", None)
                if reasoning:
                    if emit_reasoning:
                        yield {"type": "reasoning", "content": reasoning}
                    # 不 emit 时不输出思考过程（保持旧行为）
                if delta.content:
                    if emit_reasoning:
                        yield {"type": "content", "content": delta.content}
                    else:
                        yield delta.content
        except Exception as e:
            self.last_error = e
            logger.warning(f"chat_stream failed, falling back: {e}")
            # Fallback to non-streaming
            text = self.chat(messages, temperature=temperature, max_tokens=max_tokens)
            if text:
                if emit_reasoning:
                    yield {"type": "content", "content": text}
                else:
                    yield text

    def is_available(self) -> bool:
        """Check if LLM is configured (API key present)."""
        return bool(self.api_key) and len(self.api_key) > 10

    def check_connectivity(self) -> bool:
        """Quick connectivity check (lightweight ping)."""
        if not self.is_available():
            return False
        try:
            self._call_with_retry(
                lambda: self.client.models.list(),
            )
            return True
        except Exception:
            return False
