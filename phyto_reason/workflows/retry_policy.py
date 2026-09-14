"""
retry_policy.py — 重试策略。
"""

from __future__ import annotations

from datetime import datetime, timedelta
from enum import Enum


class RetryDecision(str, Enum):
    RETRY = "retry"
    SKIP = "skip"
    FAIL = "fail"


class RetryPolicy:
    """重试策略配置。

    对每个 node 可以独立配置最大重试次数和退避策略。
    """

    DEFAULT_MAX_RETRIES: int = 2
    DEFAULT_BACKOFF_SECONDS: float = 1.0

    NODE_CONFIG: dict[str, dict] = {
        "research_planning": {"max_retries": 1, "backoff": 0.5},
        "workflow_selection": {"max_retries": 1, "backoff": 0.5},
        "tool_execution": {
            "max_retries": 2,
            "backoff": 2.0,
            "non_retryable_errors": [
                "validation_failed",
                "not_found",
                "no_data",
            ],
        },
        "biological_reasoning": {"max_retries": 1, "backoff": 0.5},
        "evidence_fusion": {"max_retries": 1, "backoff": 0.5},
    }

    def __init__(self, node_name: str) -> None:
        self.node_name = node_name
        config = self.NODE_CONFIG.get(node_name, {})
        self.max_retries = config.get("max_retries", self.DEFAULT_MAX_RETRIES)
        self.backoff = config.get("backoff", self.DEFAULT_BACKOFF_SECONDS)
        self.non_retryable = config.get("non_retryable_errors", [])

    def should_retry(self, attempt: int, error_message: str = "") -> RetryDecision:
        if attempt >= self.max_retries:
            return RetryDecision.FAIL
        for err in self.non_retryable:
            if err in error_message.lower():
                return RetryDecision.SKIP
        return RetryDecision.RETRY

    def get_backoff_seconds(self, attempt: int) -> float:
        return self.backoff * (attempt + 1)
