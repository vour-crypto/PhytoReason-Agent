"""Deterministic text encoding detection for uploaded tabular files."""

from __future__ import annotations

from pathlib import Path
import re


def resolve_text_encoding(path: str | Path, requested: str = "utf-8") -> tuple[str, dict, list[str]]:
    """Return a usable encoding plus auditable provenance and warnings.

    UTF-8 is attempted first when requested.  ``charset-normalizer`` is then
    used for legacy files; common Chinese aliases are normalized to GB18030.
    A latin-1 fallback is retained only as an explicit last resort and is
    reported as such instead of silently changing the input interpretation.
    """
    path = Path(path)
    raw = path.read_bytes()
    warnings: list[str] = []
    requested = requested or "utf-8"
    try:
        raw.decode(requested)
        encoding = requested
        method = "requested_strict"
    except (UnicodeDecodeError, LookupError):
        encoding = ""
        method = "charset_normalizer"
        # This legacy heuristic is deliberately inside the strict-decode
        # failure branch. Valid UTF-8 must never be reinterpreted.
        if not encoding:
            try:
                gb_text = raw.decode("gb18030")
                if re.search(r"[\u3400-\u9fff]", gb_text):
                    encoding = "gb18030"
                    method = "gb18030_cjk_heuristic"
            except UnicodeDecodeError:
                pass
        try:
            from charset_normalizer import from_bytes

            if not encoding:
                match = from_bytes(raw).best()
                if match is not None:
                    encoding = str(match.encoding or "")
        except Exception as exc:  # pragma: no cover - optional dependency path
            warnings.append(f"charset-normalizer detection failed: {exc}")

        aliases = {"gb2312": "gb18030", "gbk": "gb18030", "cp936": "gb18030"}
        encoding = aliases.get(encoding.lower(), encoding) if encoding else ""
        if not encoding:
            for candidate in ("gb18030", "utf-8-sig", "cp1252", "latin-1"):
                try:
                    raw.decode(candidate)
                    encoding = candidate
                    method = "ordered_fallback"
                    break
                except UnicodeDecodeError:
                    continue
        if not encoding:
            raise UnicodeDecodeError("unknown", raw, 0, min(1, len(raw)), "no usable text encoding")
        warnings.append(f"Input encoding {requested!r} was not valid; detected {encoding!r}.")

    provenance = {
        "encoding": encoding,
        "encoding_requested": requested,
        "encoding_detection": method,
    }
    return encoding, provenance, warnings
