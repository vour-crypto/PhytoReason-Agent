"""SQLite application data layer for the PhytoReason desktop workbench."""

from __future__ import annotations

import sqlite3
import uuid
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterator


SCHEMA = """
PRAGMA journal_mode=WAL;
CREATE TABLE IF NOT EXISTS app_meta (key TEXT PRIMARY KEY, value TEXT NOT NULL);
CREATE TABLE IF NOT EXISTS projects (
    project_id TEXT PRIMARY KEY, name TEXT NOT NULL, description TEXT DEFAULT '',
    created_at TEXT NOT NULL, updated_at TEXT NOT NULL
);
CREATE TABLE IF NOT EXISTS sessions (
    session_id TEXT PRIMARY KEY, project_id TEXT, title TEXT DEFAULT '',
    status TEXT DEFAULT 'ready', created_at TEXT NOT NULL, updated_at TEXT NOT NULL,
    FOREIGN KEY(project_id) REFERENCES projects(project_id)
);
CREATE TABLE IF NOT EXISTS datasets (
    dataset_id TEXT PRIMARY KEY, project_id TEXT, session_id TEXT, name TEXT NOT NULL,
    dataset_type TEXT DEFAULT '', file_path TEXT DEFAULT '', created_at TEXT NOT NULL,
    FOREIGN KEY(project_id) REFERENCES projects(project_id),
    FOREIGN KEY(session_id) REFERENCES sessions(session_id)
);
CREATE TABLE IF NOT EXISTS hypotheses (
    hypothesis_id TEXT PRIMARY KEY, session_id TEXT, statement TEXT NOT NULL,
    rating TEXT NOT NULL DEFAULT 'Insufficient', confidence REAL,
    evidence_json TEXT DEFAULT '[]', created_at TEXT NOT NULL, updated_at TEXT NOT NULL,
    FOREIGN KEY(session_id) REFERENCES sessions(session_id)
);
CREATE TABLE IF NOT EXISTS hypothesis_ratings (
    rating_id INTEGER PRIMARY KEY AUTOINCREMENT, hypothesis_id TEXT NOT NULL,
    rating TEXT NOT NULL CHECK(rating IN ('Plausible', 'Weak', 'Insufficient')),
    confidence REAL, timestamp TEXT NOT NULL, note TEXT DEFAULT '',
    FOREIGN KEY(hypothesis_id) REFERENCES hypotheses(hypothesis_id)
);
CREATE TABLE IF NOT EXISTS charts (
    chart_id TEXT PRIMARY KEY, session_id TEXT, title TEXT NOT NULL,
    chart_type TEXT NOT NULL, source_refs TEXT DEFAULT '[]', file_path TEXT NOT NULL,
    created_at TEXT NOT NULL, FOREIGN KEY(session_id) REFERENCES sessions(session_id)
);
CREATE TABLE IF NOT EXISTS tasks (
    task_id TEXT PRIMARY KEY, project_id TEXT, session_id TEXT, task_type TEXT NOT NULL,
    status TEXT NOT NULL DEFAULT 'queued', progress REAL DEFAULT 0,
    message TEXT DEFAULT '', created_at TEXT NOT NULL, updated_at TEXT NOT NULL
);
CREATE TABLE IF NOT EXISTS app_settings (
    key TEXT PRIMARY KEY, value TEXT NOT NULL, updated_at TEXT NOT NULL
);
CREATE INDEX IF NOT EXISTS idx_hypothesis_ratings_time ON hypothesis_ratings(timestamp);
CREATE INDEX IF NOT EXISTS idx_charts_session ON charts(session_id);
CREATE INDEX IF NOT EXISTS idx_sessions_updated ON sessions(updated_at DESC);
"""


def default_workbench_path() -> Path:
    """Return the per-user workbench database path.

    Phase 6.4 Step 2: routed through platform_paths so the desktop shell
    shares the same user-data contract as the API (no stray ~/.phytoreason).
    """
    from phyto_reason.platform_paths import workbench_db_path
    return workbench_db_path()


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


class WorkbenchDB:
    """Small SQLite facade. ``:memory:`` is supported for tests."""

    def __init__(self, path: str | Path | None = None) -> None:
        self._memory = str(path) == ":memory:"
        self.path = Path(path) if not self._memory and path else default_workbench_path()
        self._persistent_conn: sqlite3.Connection | None = None
        if self._memory:
            self._persistent_conn = sqlite3.connect(":memory:", check_same_thread=False)
        else:
            self.path.parent.mkdir(parents=True, exist_ok=True)
        with self.connection() as conn:
            conn.executescript(SCHEMA)
            conn.execute(
                "INSERT OR REPLACE INTO app_meta(key, value) VALUES('schema_version', '1')"
            )

    @contextmanager
    def connection(self) -> Iterator[sqlite3.Connection]:
        conn = self._persistent_conn or sqlite3.connect(self.path)
        conn.row_factory = sqlite3.Row
        conn.execute("PRAGMA foreign_keys=ON")
        conn.execute("PRAGMA busy_timeout=5000")
        try:
            yield conn
            conn.commit()
        finally:
            if self._persistent_conn is None:
                conn.close()

    def counts(self) -> dict[str, int]:
        with self.connection() as conn:
            return {
                key: conn.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]
                for key, table in {
                    "projects": "projects", "datasets": "datasets",
                    "hypotheses": "hypotheses", "charts": "charts",
                }.items()
            }

    def rating_trend(self) -> list[sqlite3.Row]:
        with self.connection() as conn:
            return conn.execute(
                "SELECT timestamp, confidence, rating FROM hypothesis_ratings ORDER BY timestamp"
            ).fetchall()

    def get_setting(self, key: str, default: str = "") -> str:
        with self.connection() as conn:
            row = conn.execute("SELECT value FROM app_settings WHERE key=?", (key,)).fetchone()
            return str(row[0]) if row else default

    def set_setting(self, key: str, value: str) -> None:
        with self.connection() as conn:
            conn.execute(
                "INSERT INTO app_settings(key,value,updated_at) VALUES(?,?,?) "
                "ON CONFLICT(key) DO UPDATE SET value=excluded.value, updated_at=excluded.updated_at",
                (key, value, utc_now()),
            )

    def create_project(self, name: str, description: str = "") -> str:
        project_id = f"PRJ-{uuid.uuid4().hex[:12]}"
        now = utc_now()
        with self.connection() as conn:
            conn.execute("INSERT INTO projects VALUES(?,?,?,?,?)", (project_id, name, description, now, now))
        return project_id

    def list_projects(self) -> list[sqlite3.Row]:
        with self.connection() as conn:
            return conn.execute("SELECT * FROM projects ORDER BY updated_at DESC").fetchall()

    def create_session(self, title: str, project_id: str | None = None) -> str:
        session_id = f"SES-{uuid.uuid4().hex[:12]}"
        now = utc_now()
        with self.connection() as conn:
            conn.execute("INSERT INTO sessions VALUES(?,?,?,?,?,?)", (session_id, project_id, title, "ready", now, now))
        return session_id

    def list_sessions(self, limit: int = 50) -> list[sqlite3.Row]:
        with self.connection() as conn:
            return conn.execute("SELECT * FROM sessions ORDER BY updated_at DESC LIMIT ?", (limit,)).fetchall()

    def list_hypotheses(self) -> list[sqlite3.Row]:
        with self.connection() as conn:
            return conn.execute("SELECT * FROM hypotheses ORDER BY updated_at DESC").fetchall()

    def rating_distribution(self) -> dict[str, int]:
        distribution = {"Plausible": 0, "Weak": 0, "Insufficient": 0}
        with self.connection() as conn:
            rows = conn.execute("SELECT rating, COUNT(*) AS n FROM hypotheses GROUP BY rating").fetchall()
        for row in rows:
            if row["rating"] in distribution:
                distribution[row["rating"]] = int(row["n"])
        return distribution

    def add_dataset(self, name: str, dataset_type: str, file_path: str, session_id: str | None = None, project_id: str | None = None) -> str:
        dataset_id = f"DS-{uuid.uuid4().hex[:12]}"
        with self.connection() as conn:
            conn.execute("INSERT INTO datasets VALUES(?,?,?,?,?,?,?)", (dataset_id, project_id, session_id, name, dataset_type, file_path, utc_now()))
        return dataset_id

    def add_task(self, task_type: str, message: str = "", session_id: str | None = None) -> str:
        task_id = f"TASK-{uuid.uuid4().hex[:12]}"
        now = utc_now()
        with self.connection() as conn:
            conn.execute("INSERT INTO tasks VALUES(?,?,?,?,?,?,?,?,?)", (task_id, None, session_id, task_type, "queued", 0, message, now, now))
        return task_id

    def list_tasks(self, active_only: bool = True) -> list[sqlite3.Row]:
        with self.connection() as conn:
            if active_only:
                return conn.execute("SELECT * FROM tasks WHERE status NOT IN ('completed','failed','cancelled') ORDER BY updated_at DESC").fetchall()
            return conn.execute("SELECT * FROM tasks ORDER BY updated_at DESC").fetchall()

    def update_task_status(self, task_id: str, status: str, progress: float | None = None, message: str | None = None) -> None:
        with self.connection() as conn:
            current = conn.execute("SELECT progress, message FROM tasks WHERE task_id=?", (task_id,)).fetchone()
            if current is None:
                return
            conn.execute("UPDATE tasks SET status=?, progress=?, message=?, updated_at=? WHERE task_id=?", (status, current[0] if progress is None else progress, current[1] if message is None else message, utc_now(), task_id))

    def created_since(self, iso_cutoff: str) -> int:
        """Count hypotheses created after the given ISO-8601 cutoff."""
        with self.connection() as conn:
            row = conn.execute(
                "SELECT COUNT(*) FROM hypotheses WHERE created_at >= ?",
                (iso_cutoff,),
            ).fetchone()
            return int(row[0])
