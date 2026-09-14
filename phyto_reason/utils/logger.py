"""
logging configuration for production-grade scientific logging.
"""

from __future__ import annotations

import logging
import os
import sys
from datetime import datetime
from pathlib import Path

from phyto_reason.platform_paths import log_dir

LOG_DIR = log_dir()


def setup_logging(name: str = "plantomics", level: int = logging.INFO,
                   log_file: str | None = None) -> logging.Logger:
    """配置科研日志系统。

    同时输出到:
      - 文件 (logs/plantomics_YYYYMMDD.log)
      - stdout
    """
    logger = logging.getLogger(name)
    logger.setLevel(level)

    if logger.handlers:
        return logger

    formatter = logging.Formatter(
        fmt="%(asctime)s | %(levelname)-8s | %(name)s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    if log_file is None:
        LOG_DIR.mkdir(parents=True, exist_ok=True)
        log_file = str(LOG_DIR / f"{name}_{datetime.now():%Y%m%d}.log")

    fh = logging.FileHandler(log_file, encoding="utf-8")
    fh.setLevel(level)
    fh.setFormatter(formatter)
    logger.addHandler(fh)

    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(level)
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    return logger
