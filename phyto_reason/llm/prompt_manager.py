"""
prompt_manager.py — Prompt 模板管理器。

管理 system prompts 和 reasoning prompts 的填充与组合。
"""

from __future__ import annotations

from typing import Any

from phyto_reason.llm.system_prompts import SYSTEM_PROMPTS
from phyto_reason.llm.reasoning_prompts import REASONING_PROMPTS


class PromptManager:
    """Prompt 模板管理器。"""

    @staticmethod
    def build_messages(
        system_key: str,
        reasoning_key: str,
        user_message: str,
        template_vars: dict[str, Any] | None = None,
        conversation_history: list[dict] | None = None,
    ) -> list[dict]:
        """构建 LLM messages 列表。

        Args:
            system_key: SYSTEM_PROMPTS 中的 key
            reasoning_key: REASONING_PROMPTS 中的 key
            user_message: 用户消息 (如果 reasoning template 不需要单独的 user msg)
            template_vars: 模板变量
            conversation_history: 对话历史

        Returns:
            messages: 适用于 LLM chat completion 的 message list
        """
        messages = []

        system = SYSTEM_PROMPTS.get(system_key, SYSTEM_PROMPTS["default"])
        messages.append({"role": "system", "content": system})

        if conversation_history:
            for turn in conversation_history[-10:]:
                messages.append({"role": "user", "content": turn.get("user", "")})
                messages.append({"role": "assistant", "content": turn.get("assistant", "")})

        prompt = REASONING_PROMPTS.get(reasoning_key, "")
        if template_vars:
            try:
                prompt = prompt.format(**template_vars)
            except KeyError as e:
                prompt = f"{prompt}\n\n[Template error: missing {e}]"

        if prompt:
            messages.append({"role": "user", "content": prompt})
        else:
            messages.append({"role": "user", "content": user_message})

        return messages

    @staticmethod
    def render(template_key: str, **kwargs) -> str:
        """直接渲染 prompt 模板。"""
        template = REASONING_PROMPTS.get(template_key, "")
        try:
            return template.format(**kwargs)
        except KeyError as e:
            return f"Template error: missing {e}"
