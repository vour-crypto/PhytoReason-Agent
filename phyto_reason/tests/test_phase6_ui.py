from pathlib import Path
import re


STATIC = Path(__file__).resolve().parents[1] / "api" / "static"


def test_phase6_assets_are_split_without_inline_code():
    html = (STATIC / "index.html").read_text(encoding="utf-8")
    html_without_comments = re.sub(r"<!--.*?-->", "", html, flags=re.S)
    assert '<link rel="stylesheet" href="/static/assets/app.css">' in html
    assert '<script src="/static/assets/app.js"></script>' in html
    assert not re.search(r"<style(?:\s[^>]*)?>", html_without_comments, re.I)
    assert not re.search(r"<script(?!\s+src=)(?:\s[^>]*)?>", html_without_comments, re.I)


def test_phase6_css_retains_frozen_visual_tokens():
    css = (STATIC / "assets" / "app.css").read_text(encoding="utf-8")
    for token in (
        "#f4f3ef",
        "#e6cdb1",
        "#d4e0ce",
        "#b45309",
        "border-radius:16px",
        "font-family:'Microsoft YaHei UI'",
        "210px",
        "340px",
    ):
        assert token in css


def test_phase6_frontend_wires_are_local_and_sse_enabled():
    js = (STATIC / "assets" / "app.js").read_text(encoding="utf-8")
    assert 'fetch("/upload"' in js
    assert 'fetch("/chat/stream"' in js
    assert 'fetch("/health"' in js
    assert "/static/figures/" in js
    assert "capsule" in js and "degraded" in js and "done" in js
    assert not re.search(r"https?://", js)


def test_phase6_baseline_contains_no_demo_terms():
    source = "\n".join(
        (STATIC / name).read_text(encoding="utf-8")
        for name in ("index.html", "assets/app.css", "assets/app.js")
    )
    for term in ("水稻", "OsNAC", "盐胁迫", "置信度 6"):
        assert term not in source
