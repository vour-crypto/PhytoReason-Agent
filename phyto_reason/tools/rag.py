"""
RAG — 基于 TF-IDF 的语义检索知识库。

扫描 outputs/ 下所有 CSV/Markdown，建为向量索引，
支持语义搜索。无需 GPU / PyTorch。
"""

from __future__ import annotations

import json
import math
import os
import re
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from phyto_reason.config.settings import Settings


def _tokenize(text: str) -> list[str]:
    return re.findall(r"[a-zA-Z\u4e00-\u9fff]+", text.lower())


def _csv_chunks(path: Path, rel: str) -> list[dict]:
    try:
        df = pd.read_csv(path)
    except Exception:
        return []
    if df.empty or df.shape[1] < 2:
        return []
    chunks = []
    for col in df.columns:
        vals = df[col].dropna().astype(str).str[:200].unique()
        top = [v for v in vals if v not in ("", "--", "nan")]
        if not top:
            continue
        chunks.append({
            "text": f"【{rel}】列「{col}」: {'、'.join(top[:15])}",
            "meta": {"source": str(rel), "column": col},
        })
    return chunks


def _md_chunks(text: str, rel: str) -> list[dict]:
    chunks = []
    section = "概述"
    lines: list[str] = []
    for line in text.split("\n"):
        if line.startswith("## "):
            if lines:
                chunks.append({"text": "\n".join(lines).strip(), "meta": {"source": rel, "section": section}})
            section = line.strip("# ").strip()
            lines = [line]
        else:
            lines.append(line)
    if lines:
        chunks.append({"text": "\n".join(lines).strip(), "meta": {"source": rel, "section": section}})
    return chunks


class RAGEngine:
    """TF-IDF 语义检索引擎。"""

    def __init__(self, settings: Settings | None = None) -> None:
        self.settings = settings or Settings()
        self.chunks: list[dict] = []
        self.vocab: dict[str, int] = {}
        self._vectors: np.ndarray | None = None
        self._ready = False

    def build(self, scan_dir: Path | None = None) -> int:
        """构建索引。"""
        base = scan_dir or self.settings.output_dir
        chunks: list[dict] = []

        for csv_file in sorted(base.rglob("*.csv")):
            if csv_file.stat().st_size == 0:
                continue
            rel = csv_file.relative_to(self.settings.project_root)
            try:
                chunks.extend(_csv_chunks(csv_file, str(rel)))
            except Exception:
                continue

        for md_file in sorted(base.rglob("*.md")):
            try:
                rel = md_file.relative_to(self.settings.project_root)
                chunks.extend(_md_chunks(md_file.read_text(encoding="utf-8"), str(rel)))
            except Exception:
                continue

        return self._build_tfidf(chunks)

    def _build_tfidf(self, chunks: list[dict]) -> int:
        self.chunks = chunks
        if not chunks:
            self._ready = False
            return 0

        stopwords = {"the", "a", "an", "is", "are", "was", "were", "of", "in",
                     "to", "and", "for", "on", "with", "as", "at", "by", "or"}
        texts = [c["text"] for c in chunks]
        doc_freq: Counter = Counter()
        all_tokens: list[list[str]] = []
        for t in texts:
            tokens = [w for w in _tokenize(t) if w not in stopwords and len(w) > 1]
            all_tokens.append(tokens)
            doc_freq.update(set(tokens))

        top_words = [w for w, _ in doc_freq.most_common(2000)]
        self.vocab = {w: i for i, w in enumerate(top_words)}

        n = len(chunks)
        idf = np.zeros(len(self.vocab))
        for w, i in self.vocab.items():
            df = doc_freq[w]
            idf[i] = math.log((n + 1) / (df + 1)) + 1

        vecs = np.zeros((n, len(self.vocab)), dtype=np.float32)
        for i, tokens in enumerate(all_tokens):
            tf = Counter(tokens)
            for w, c in tf.items():
                if w in self.vocab:
                    vecs[i, self.vocab[w]] = c * idf[self.vocab[w]]

        norms = np.linalg.norm(vecs, axis=1, keepdims=True)
        norms[norms == 0] = 1
        self._vectors = vecs / norms
        self._ready = True
        return len(chunks)

    def search(self, query: str, top_k: int = 8) -> list[dict[str, Any]]:
        if not self._ready or self._vectors is None:
            return []

        qt = _tokenize(query)
        qv = np.zeros(len(self.vocab), dtype=np.float32)
        for w in qt:
            if w in self.vocab:
                qv[self.vocab[w]] += 1
        qn = np.linalg.norm(qv)
        if qn > 0:
            qv /= qn

        scores = self._vectors @ qv
        top = np.argsort(scores)[::-1][:top_k]
        results = []
        for idx in top:
            if scores[idx] <= 0:
                continue
            results.append({
                "score": float(scores[idx]),
                "text": self.chunks[idx]["text"],
                "metadata": self.chunks[idx].get("meta", {}),
            })
        return results

    def format_context(self, results: list[dict], max_chars: int = 4000) -> str:
        parts = []
        total = 0
        for r in results:
            meta = r.get("metadata", {})
            src = meta.get("source", "?")
            col = meta.get("column", "")
            ext = f" | {col}" if col else ""
            block = f"[{src}{ext} | score: {r['score']:.3f}]\n{r['text']}\n"
            total += len(block)
            if total > max_chars:
                break
            parts.append(block)
        return "\n".join(parts) if parts else "（无相关结果）"


# ── Singleton ──────────────────────────────────────────────────

_rag_engine: RAGEngine | None = None


def get_rag_engine() -> RAGEngine:
    """Get or create the global RAG engine singleton (lazy init).

    Scans outputs/ for CSV/MD files and builds TF-IDF index on first call.
    Returns an engine with an empty index if no documents are found.
    """
    global _rag_engine
    if _rag_engine is None:
        import logging
        _rag_engine = RAGEngine()
        n = _rag_engine.build()
        if n > 0:
            logging.getLogger("rag").info("RAG index built: %d chunks", n)
        else:
            logging.getLogger("rag").warning(
                "RAG index is empty — no CSV/MD files found in outputs/. "
                "Run a pipeline first to populate the index."
            )
    return _rag_engine


def refresh_rag_engine() -> int:
    """Force rebuild the RAG index (call after new pipeline results are saved)."""
    global _rag_engine
    import logging
    engine = RAGEngine()
    n = engine.build()
    _rag_engine = engine
    logging.getLogger("rag").info("RAG index refreshed: %d chunks", n)
    return n
