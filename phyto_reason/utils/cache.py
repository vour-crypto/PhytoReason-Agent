"""
cache.py — In-memory TTL cache for tool and knowledge layer results.

Avoids recomputing expensive operations within the same session:
  - KEGG API calls (network round-trip)
  - PubMed searches (rate-limited API)
  - Correlation computations (CPU-intensive)
  - Motif scans (FIMO subprocess)

Usage:
    from phyto_reason.utils.cache import result_cache

    @result_cache(ttl=3600, key_fn=lambda **kw: f"kegg:{kw.get('compound')}")
    def query_kegg(compound, ...):
        ...

    # Or programmatically:
    result_cache.get("key")
    result_cache.set("key", value, ttl=600)
"""

from __future__ import annotations

import hashlib
import json
import time
import threading
from functools import wraps
from typing import Any, Callable


class TTLCache:
    """Thread-safe in-memory TTL cache."""

    def __init__(self, name: str = "default") -> None:
        self.name = name
        self._store: dict[str, tuple[Any, float]] = {}  # key → (value, expires_at)
        self._lock = threading.Lock()
        self._hits = 0
        self._misses = 0

    def get(self, key: str) -> Any | None:
        """Get a cached value. Returns None if expired or missing."""
        with self._lock:
            entry = self._store.get(key)
            if entry is None:
                self._misses += 1
                return None
            value, expires_at = entry
            if expires_at > 0 and time.time() > expires_at:
                del self._store[key]
                self._misses += 1
                return None
            self._hits += 1
            return value

    def set(self, key: str, value: Any, ttl: int = 3600) -> None:
        """Store a value with TTL in seconds. ttl=0 means no expiry."""
        with self._lock:
            expires_at = time.time() + ttl if ttl > 0 else 0
            self._store[key] = (value, expires_at)

    def delete(self, key: str) -> None:
        """Remove a key from the cache."""
        with self._lock:
            self._store.pop(key, None)

    def clear(self) -> None:
        """Clear all cached entries."""
        with self._lock:
            self._store.clear()

    def stats(self) -> dict:
        """Return hit/miss statistics."""
        with self._lock:
            total = self._hits + self._misses
            return {
                "name": self.name,
                "entries": len(self._store),
                "hits": self._hits,
                "misses": self._misses,
                "hit_rate": self._hits / total if total > 0 else 0.0,
            }

    def __len__(self) -> int:
        return len(self._store)


# ═══════════════════════════════════════════════════════════════
# Global caches
# ═══════════════════════════════════════════════════════════════

# Tool execution cache (shared by all tools)
tool_cache = TTLCache("tools")

# Pipeline-level cache (for full pipeline results on same input)
pipeline_cache = TTLCache("pipeline")

# Knowledge layer cache
knowledge_cache = TTLCache("knowledge")


# ═══════════════════════════════════════════════════════════════
# Decorator
# ═══════════════════════════════════════════════════════════════

def cached(cache: TTLCache, ttl: int = 3600, key_fn: Callable | None = None):
    """Decorator: cache function results.

    Args:
        cache: TTLCache instance to use
        ttl: Time-to-live in seconds
        key_fn: Optional function(**kwargs) → str to build cache key.
                Default: md5(str(sorted(kwargs.items())))
    """
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            if key_fn:
                cache_key = key_fn(**kwargs)
            else:
                # Default: hash the sorted kwargs + function name
                raw = json.dumps({
                    "func": func.__name__,
                    "args": str(args[1:]) if args else "",
                    "kwargs": sorted(kwargs.items()),
                }, sort_keys=True, default=str)
                cache_key = hashlib.md5(raw.encode()).hexdigest()

            cached_val = cache.get(cache_key)
            if cached_val is not None:
                return cached_val

            result = func(*args, **kwargs)
            cache.set(cache_key, result, ttl=ttl)
            return result
        return wrapper
    return decorator


# ═══════════════════════════════════════════════════════════════
# Convenience: cache key builders for common tool calls
# ═══════════════════════════════════════════════════════════════

def tool_key(tool_name: str, **kwargs) -> str:
    """Build a cache key for a tool execution."""
    raw = json.dumps({"tool": tool_name, "args": sorted(kwargs.items())}, sort_keys=True, default=str)
    return hashlib.md5(raw.encode()).hexdigest()


def pipeline_key(species: str, metabolite: str, expression_hash: str = "", n_samples: int = 0) -> str:
    """Build a cache key for pipeline results."""
    return f"pipeline:{species}:{metabolite}:{expression_hash}:n={n_samples}"


def knowledge_key(domain: str, query: str) -> str:
    """Build a cache key for knowledge layer queries."""
    return f"{domain}:{query.lower().strip()}"
