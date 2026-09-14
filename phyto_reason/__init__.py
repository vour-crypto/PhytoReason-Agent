"""
PhytoReason-Agent v5.0
多物种药用植物次生代谢调控科研推理 Agent
"""
from __future__ import annotations

import sys
from pathlib import Path

_root = Path(__file__).resolve().parent.parent
if str(_root) not in sys.path:
    sys.path.insert(0, str(_root))

from phyto_reason.agents.cli import main

__version__ = "5.0.0"
__all__ = ["main"]
