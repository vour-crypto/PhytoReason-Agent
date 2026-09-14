from __future__ import annotations

import os
import urllib.request

import pytest

pytest.importorskip("PySide6")
pytest.importorskip("PySide6.QtWebEngineWidgets")

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ.setdefault("QTWEBENGINE_CHROMIUM_FLAGS", "--disable-gpu --no-sandbox")

from PySide6.QtWidgets import QApplication

from phyto_reason.desktop.webview_host import WebViewWindow


def test_webview_shell_starts_health_and_stops_cleanly():
    app = QApplication.instance() or QApplication([])
    window = WebViewWindow()
    try:
        assert window.server.wait_ready(timeout=10)
        with urllib.request.urlopen(window.server.url.rstrip("/") + "/health", timeout=2) as response:
            assert response.status == 200
        assert window.width() == 1380
        assert window.height() == 860
    finally:
        window.close()
        app.processEvents()
    assert window.server.thread is not None
    assert not window.server.thread.is_alive()
