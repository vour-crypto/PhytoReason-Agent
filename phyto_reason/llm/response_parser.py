"""
response_parser.py — LLM 结构化输出解析器。

将 LLM 的文本/JSON 输出解析为结构化科研数据。
"""

from __future__ import annotations

import json
import re
from typing import Any

from pydantic import BaseModel, Field


class StructuredResult(BaseModel):
    raw_text: str = ""
    parsed: bool = False
    confidence: str = ""
    reasoning_chain: list[str] = Field(default_factory=list)
    hypothesis: str = ""
    suggestions: list[str] = Field(default_factory=list)
    evidence_summary: str = ""
    uncertainty_notes: list[str] = Field(default_factory=list)
    raw_json: dict = Field(default_factory=dict)


class ResponseParser:
    """LLM 响应解析器。"""

    @staticmethod
    def parse_structured(text: str) -> StructuredResult:
        """尝试从 LLM 响应中解析结构化数据。"""
        result = StructuredResult(raw_text=text)

        # 尝试 JSON 解析
        json_data = ResponseParser._extract_json(text)
        if json_data:
            result.raw_json = json_data
            result.parsed = True
            if "confidence" in json_data:
                result.confidence = str(json_data["confidence"])
            if "reasoning_chain" in json_data:
                result.reasoning_chain = json_data["reasoning_chain"]
            if "hypothesis" in json_data:
                result.hypothesis = json_data["hypothesis"]
            if "suggestions" in json_data:
                result.suggestions = json_data["suggestions"]
            if "uncertainty_notes" in json_data:
                result.uncertainty_notes = json_data["uncertainty_notes"]
            return result

        # 非 JSON 解析
        result.parsed = False
        lines = text.strip().split("\n")
        current_section = ""
        for line in lines:
            stripped = line.strip()
            if stripped.lower().startswith("confidence"):
                result.confidence = ResponseParser._extract_value(stripped)
                current_section = "confidence"
            elif stripped.lower().startswith("hypothesis"):
                result.hypothesis = ResponseParser._extract_value(stripped)
                current_section = "hypothesis"
            elif stripped.lower().startswith("suggestion") or stripped.lower().startswith("recommend"):
                result.suggestions.append(ResponseParser._extract_value(stripped))
                current_section = "suggestions"
            elif stripped.lower().startswith("uncertainty"):
                result.uncertainty_notes.append(ResponseParser._extract_value(stripped))
                current_section = "uncertainty"
            elif stripped.startswith("- ") and current_section == "suggestions":
                result.suggestions.append(stripped[2:])
            elif stripped.startswith("- ") and current_section == "uncertainty":
                result.uncertainty_notes.append(stripped[2:])
            elif stripped.startswith("- ") and current_section in ("", "confidence", "hypothesis"):
                result.reasoning_chain.append(stripped[2:])

        return result

    @staticmethod
    def parse_tool_calls(text: str) -> list[dict]:
        """从文本中提取 tool call JSON。"""
        calls = []
        for match in re.finditer(r'```(?:json)?\s*({.*?})\s*```', text, re.DOTALL):
            try:
                data = json.loads(match.group(1))
                if "tool_name" in data or "name" in data:
                    calls.append(data)
            except json.JSONDecodeError:
                pass
        return calls

    @staticmethod
    def _extract_json(text: str) -> dict | None:
        """尝试从文本中提取 JSON 块。"""
        for match in re.finditer(r'```(?:json)?\s*({.*?})\s*```', text, re.DOTALL):
            try:
                return json.loads(match.group(1))
            except json.JSONDecodeError:
                pass
        try:
            return json.loads(text.strip())
        except json.JSONDecodeError:
            pass
        return None

    @staticmethod
    def _extract_value(line: str) -> str:
        """从 'key: value' 行提取值。"""
        if ":" in line:
            return line.split(":", 1)[1].strip()
        return line.strip()
