"""PySide6 desktop shell for the local PhytoReason Web UI."""

from __future__ import annotations

import socket
import threading
import time
import urllib.error
import urllib.request
from typing import Any

try:
    import uvicorn
    from PySide6.QtCore import QUrl, QTimer, Qt
    from PySide6.QtGui import QDesktopServices
    from PySide6.QtWidgets import QApplication, QLabel, QMainWindow, QStackedWidget
    from PySide6.QtWebEngineCore import QWebEnginePage
    from PySide6.QtWebEngineWidgets import QWebEngineView
except ImportError:  # pragma: no cover - exercised when the desktop extra is absent
    uvicorn = None
    QApplication = None
    QUrl = QTimer = Qt = QDesktopServices = QLabel = QMainWindow = QStackedWidget = None
    QWebEnginePage = QWebEngineView = None


def find_free_port() -> int:
    """Ask the OS for a currently unused local TCP port."""
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
        sock.bind(("127.0.0.1", 0))
        return int(sock.getsockname()[1])


class LocalApiServer:
    """Run the existing FastAPI app in a stoppable background thread."""

    def __init__(self, port: int | None = None) -> None:
        self.port = port or find_free_port()
        self.url = f"http://127.0.0.1:{self.port}/"
        self.server: Any = None
        self.thread: threading.Thread | None = None
        self.error: Exception | None = None
        self._stop_requested = threading.Event()

    def start(self) -> None:
        if uvicorn is None:
            self.error = RuntimeError("桌面依赖未安装，请运行 pip install 'phyto-reason-agent[desktop]'")
            return

        def serve() -> None:
            try:
                config = uvicorn.Config(
                    "phyto_reason.api.app:app",
                    host="127.0.0.1",
                    port=self.port,
                    log_level="error",
                    access_log=False,
                    log_config=None,
                )
                self.server = uvicorn.Server(config)
                if self._stop_requested.is_set():
                    self.server.should_exit = True
                self.server.run()
            except Exception as exc:  # pragma: no cover - startup failures are environment-specific
                self.error = exc

        self.thread = threading.Thread(target=serve, name="phyto-reason-api", daemon=True)
        self.thread.start()

    def wait_ready(self, timeout: float = 10.0) -> bool:
        """Poll the existing health endpoint until the local app is ready."""
        deadline = time.monotonic() + timeout
        health_url = self.url.rstrip("/") + "/health"
        while time.monotonic() < deadline:
            if self.error is not None:
                return False
            try:
                with urllib.request.urlopen(health_url, timeout=0.5) as response:
                    if response.status == 200:
                        return True
            except (OSError, urllib.error.URLError):
                time.sleep(0.05)
        return False

    def stop(self) -> None:
        self._stop_requested.set()
        if self.server is not None:
            self.server.should_exit = True
        if self.thread is not None and self.thread.is_alive():
            self.thread.join(timeout=5)


if QWebEnginePage is not None:
    class LocalWebPage(QWebEnginePage):
        """Keep the shell local and hand external navigations to the OS browser."""

        def __init__(self, local_port: int, parent=None) -> None:
            super().__init__(parent)
            self.local_port = local_port

        def acceptNavigationRequest(self, url, navigation_type, is_main_frame):  # noqa: N802
            if url.scheme() in {"http", "https"}:
                is_local = url.host() in {"127.0.0.1", "localhost"} and url.port() == self.local_port
                if not is_local:
                    QDesktopServices.openUrl(url)
                    return False
            return super().acceptNavigationRequest(url, navigation_type, is_main_frame)


if QMainWindow is not None:
    class WebViewWindow(QMainWindow):
        """Desktop window that hosts the existing Web UI without changing it."""

        def __init__(self) -> None:
            super().__init__()
            self.setWindowTitle("PhytoReason")
            self.resize(1380, 860)
            self.server = LocalApiServer()
            self.stack = QStackedWidget(self)
            self.starting_label = QLabel("正在启动...", self)
            self.starting_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
            self.error_label = QLabel(self)
            self.error_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
            self.error_label.setWordWrap(True)
            self.webview = QWebEngineView(self)
            self.webview.setPage(LocalWebPage(self.server.port, self.webview))
            self.stack.addWidget(self.starting_label)
            self.stack.addWidget(self.error_label)
            self.stack.addWidget(self.webview)
            self.setCentralWidget(self.stack)
            self._ready_timer = QTimer(self)
            self._ready_timer.timeout.connect(self._check_backend)
            self.server.start()
            self._ready_timer.start(100)

        def _check_backend(self) -> None:
            if self.server.error is not None:
                self._show_error(self.server.error)
                return
            if not self.server.wait_ready(timeout=0.01):
                if self.server.thread is not None and not self.server.thread.is_alive():
                    self._show_error(RuntimeError("本地 API 线程已退出"))
                return
            self._ready_timer.stop()
            self.stack.setCurrentWidget(self.webview)
            self.webview.load(QUrl(self.server.url))

        def _show_error(self, error: Exception) -> None:
            self._ready_timer.stop()
            self.error_label.setText(
                "PhytoReason 启动失败\n\n"
                f"{error}\n\n"
                "请检查桌面依赖、Python 环境和端口占用后重试。"
            )
            self.stack.setCurrentWidget(self.error_label)

        def closeEvent(self, event) -> None:  # noqa: N802
            self._ready_timer.stop()
            self.webview.stop()
            self.server.stop()
            event.accept()
else:
    WebViewWindow = None


def run() -> None:
    if QApplication is None or QWebEngineView is None:
        raise SystemExit("Install the desktop extra: pip install 'phyto-reason-agent[desktop]'")
    app = QApplication.instance() or QApplication([])
    window = WebViewWindow()
    window.show()
    app.exec()


if __name__ == "__main__":
    run()
