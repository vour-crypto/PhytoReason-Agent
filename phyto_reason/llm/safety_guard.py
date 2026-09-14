"""
safety_guard.py — LLM 输出安全守卫。

防止:
  - 编造生物学结论
  - 伪造文献引用
  - 未经支持的因果声明
  - 编造验证实验
"""

from __future__ import annotations

import re
from typing import Any

from pydantic import BaseModel, Field


class SafetyVerdict(BaseModel):
    passed: bool = True
    warnings: list[str] = Field(default_factory=list)
    blocked_phrases: list[str] = Field(default_factory=list)
    modified: bool = False


class SafetyGuard:
    """LLM 输出安全守卫。"""

    # 禁止的绝对化断言
    ABSOLUTE_CLAIMS = [
        r"\bproves?\b",
        r"\bdefinitely\b",
        r"\bcertainly\b",
        r"\bwithout any doubt\b",
        r"\bconfirmed that\b",
        r"\bit is known that\b",
        r"\bestablished fact\b",
        r"\bunquestionably\b",
    ]

    # 需要标注的 speculative 词汇
    SPECULATIVE_WORDS = [
        "suggests", "may", "might", "could", "possibly",
        "indicates", "is consistent with", "raises the possibility",
    ]

    # 已知可能被编造的 PMID 格式
    FAKE_PMID_PATTERN = re.compile(r"PMID[:\s]*(\d{8,})")

    # 禁止的虚假实验声明
    FABRICATED_EXPERIMENTS = [
        r"our previous (study|work|experiment)s? (showed|demonstrated|proved)",
        r"we have (already|previously) (validated|confirmed|demonstrated)",
        r"in our (lab|laboratory|unpublished) (data|results|finding)",
    ]

    @staticmethod
    def check(response: str, context: dict[str, Any] | None = None) -> SafetyVerdict:
        """检查 LLM 响应安全性。

        Args:
            response: LLM 生成的响应文本
            context: 可选的上下文 (用于验证 claims)

        Returns:
            SafetyVerdict
        """
        verdict = SafetyVerdict()

        # 1. 检查绝对化断言
        for pattern in SafetyGuard.ABSOLUTE_CLAIMS:
            if re.search(pattern, response, re.IGNORECASE):
                verdict.warnings.append(f"绝对化断言被标记: {pattern}")
                verdict.blocked_phrases.append(pattern)

        # 2. 检查虚假文献号
        fake_pmids = SafetyGuard.FAKE_PMID_PATTERN.findall(response)
        for pmid in fake_pmids:
            if len(pmid) > 8:
                verdict.warnings.append(f"可疑的 PMID 编号: {pmid} (位数过长)")
                verdict.blocked_phrases.append(f"PMID:{pmid}")

        # 3. 检查编造实验声明
        for pattern in SafetyGuard.FABRICATED_EXPERIMENTS:
            if re.search(pattern, response, re.IGNORECASE):
                verdict.warnings.append(f"编造实验声明: {pattern}")
                verdict.blocked_phrases.append(pattern)

        # 4. 检查是否使用了 speculative 词汇 (好现象)
        speculative_count = sum(
            1 for w in SafetyGuard.SPECULATIVE_WORDS
            if w.lower() in response.lower()
        )

        # 检查是否表达了不确定性
        has_uncertainty = speculative_count >= 1
        if not has_uncertainty and not verdict.warnings:
            # 没有 speculative 词汇也没有 warning
            pass  # response is aggressive but not necessarily wrong

        verdict.passed = len(verdict.blocked_phrases) == 0
        return verdict

    @staticmethod
    def sanitize(response: str) -> str:
        """清理有问题的内容 (替换而非完全阻止)。"""
        result = response

        # 替换绝对化断言
        replacements = [
            (r"\bproves?\b", "strongly suggests"),
            (r"\bdefinitely\b", "with high confidence"),
            (r"\bcertainly\b", "is highly likely"),
            (r"\bconfirmed that\b", "evidence indicates that"),
        ]
        for pattern, replacement in replacements:
            result = re.sub(pattern, replacement, result, flags=re.IGNORECASE)

        # 移除虚假 PMID
        result = SafetyGuard.FAKE_PMID_PATTERN.sub("[PMID reference]", result)

        return result
