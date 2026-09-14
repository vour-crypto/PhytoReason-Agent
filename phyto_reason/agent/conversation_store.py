"""
conversation_store.py — Session-level conversation state and data storage (v4.3).

Stores:
- Conversation history (user/assistant turns)
- Uploaded data (expression, metabolite, promoter matrices)
- Pipeline results (cached for follow-up questions)

v4.3: SQLite persistence — sessions survive server restarts.
      Write-through cache: all mutations persist to DB immediately.
      On startup, existing sessions are loaded from DB into memory.
"""

from __future__ import annotations

import json
import logging
import sqlite3
import threading
import uuid
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any

from phyto_reason.models.workflow_plan import WorkflowPlan

logger = logging.getLogger("conversation_store")

# ── DB schema ──────────────────────────────────────────────────

CREATE_SESSIONS_TABLE = """
CREATE TABLE IF NOT EXISTS sessions (
    session_id     TEXT PRIMARY KEY,
    created_at     TEXT NOT NULL,
    species        TEXT DEFAULT '',
    target_metabolite TEXT DEFAULT '',
    target_pathway TEXT DEFAULT '',
    has_data       INTEGER DEFAULT 0,
    turn_count     INTEGER DEFAULT 0,
    history        TEXT DEFAULT '[]',
    message_chain  TEXT DEFAULT '[]',
    expression_matrix    TEXT DEFAULT NULL,
    metabolite_matrix    TEXT DEFAULT NULL,
    promoter_sequences   TEXT DEFAULT NULL,
    sample_metadata      TEXT DEFAULT NULL,
    sample_alignment_report TEXT DEFAULT '{}',
    pipeline_result      TEXT DEFAULT '',
    pipeline_hypotheses  TEXT DEFAULT '[]',
    pipeline_mechanism_type TEXT DEFAULT '',
    pipeline_pathway     TEXT DEFAULT ''
);
"""

CREATE_INDEX = """
CREATE INDEX IF NOT EXISTS idx_sessions_created
ON sessions(created_at DESC);
"""

# ── Default DB path (project_root/data/sessions.db) ────────────

def _default_db_path() -> Path:
    """Resolve the per-user session database path."""
    from phyto_reason.platform_paths import sessions_db_path
    return sessions_db_path()


# ═══════════════════════════════════════════════════════════════
# SessionData
# ═══════════════════════════════════════════════════════════════

@dataclass
class SessionData:
    """Per-session state with optional persistence callback."""

    session_id: str
    created_at: str = field(default_factory=lambda: datetime.now().isoformat())

    # Conversation history: list of {"role": "user"|"assistant", "content": str}
    history: list[dict] = field(default_factory=list)

    # User-provided context (extracted by LLM over the conversation)
    species: str = ""
    target_metabolite: str = ""
    target_pathway: str = ""

    # Uploaded data
    expression_matrix: dict | None = None
    metabolite_matrix: dict | None = None
    promoter_sequences: dict | None = None
    sample_metadata: dict | None = None
    sample_alignment_report: dict = field(default_factory=dict)
    has_data: bool = False

    # Pipeline results (cached for follow-up questions)
    pipeline_result: str = ""
    pipeline_hypotheses: list = field(default_factory=list)
    pipeline_mechanism_type: str = ""

    # 报告支持：最近一次管线运行的分节摘要 + 图表编号注册表（内存级，重启后重算）
    analysis_digest: dict = field(default_factory=dict)
    figure_captions: dict = field(default_factory=dict)

    # TF 预测：上传 FASTA 生成的注释 CSV 路径（会话内注入 TF 分析链）
    tf_annotation_csv: str = ""

    # MS/MS 谱图注释表（数据页槽位⑥生成）
    ms2_annotation_csv: str = ""
    pipeline_pathway: str = ""

    # Workflow plan for step-selective analysis (in-memory, not persisted)
    workflow_plan: WorkflowPlan | None = field(default=None, repr=False, compare=False)

    # Structured candidate scores for visualization (not persisted, derived from pipeline)
    candidate_scores: dict = field(default_factory=dict)
    # Format: {"WRKY1": {"correlation": 0.85, "motif": 0.6, "pathway": 0.7, ...}, ...}

    # DEG/DAM data for volcano plots (not persisted, loaded from pipeline outputs)
    deg_data: dict = field(default_factory=dict)
    # Format: {"GeneID": {"log2fc": 2.5, "padj": 0.001}, ...}
    dam_data: dict = field(default_factory=dict)
    # Format: {"Metabolite": {"log2fc": 1.8, "pvalue": 0.003}, ...}

    # v5.0 ReAct agent: cached pipeline stage results (in-memory, not persisted)
    # Each field is populated by the corresponding tool handler and read by downstream tools
    deg_report: dict = field(default_factory=dict)
    dam_report: dict = field(default_factory=dict)
    multiomics_report: dict = field(default_factory=dict)
    wgcna_report: dict = field(default_factory=dict)
    qc_report: dict = field(default_factory=dict)
    tf_candidates: dict = field(default_factory=dict)

    # Turn counter
    turn_count: int = 0

    # Full message chain from the most recent turn (preserves tool calls + results)
    # Format: list of {"role": ..., "content": ..., "tool_calls": [...], "tool_call_id": ...}
    # Stored so follow-up turns have full context of what the agent discovered
    message_chain: list[dict] = field(default_factory=list)

    # Back-reference to store for persistence (set by ConversationStore)
    _store: ConversationStore | None = field(default=None, repr=False, compare=False)

    # ── Mutation methods (self-persisting) ─────────────────

    def add_turn(self, role: str, content: str) -> None:
        """Record a conversation turn."""
        self.history.append({"role": role, "content": content})
        self.turn_count += 1

        # Keep history manageable (last 20 turns = 40 messages)
        if len(self.history) > 40:
            self.history = self.history[-40:]

        self._persist()

    def set_data(
        self,
        expression: dict | None = None,
        metabolite: dict | None = None,
        promoter: dict | None = None,
        sample_metadata: dict | None = None,
        sample_alignment_report: dict | None = None,
    ) -> None:
        """Store uploaded data matrices."""
        if expression is not None:
            self.expression_matrix = expression
        if metabolite is not None:
            self.metabolite_matrix = metabolite
        if promoter is not None:
            self.promoter_sequences = promoter
        if sample_metadata is not None:
            self.sample_metadata = sample_metadata
        if sample_alignment_report is not None:
            self.sample_alignment_report = sample_alignment_report
        self.has_data = any((
            bool(self.expression_matrix),
            bool(self.metabolite_matrix),
            bool(self.promoter_sequences),
        ))
        self._persist()

    def set_pipeline_result(self, result_text: str, hypotheses: list | None = None) -> None:
        """Cache the last pipeline result for follow-up questions."""
        self.pipeline_result = result_text
        if hypotheses is not None:
            self.pipeline_hypotheses = hypotheses
        self._persist()

    def update_context(self, species: str = "", metabolite: str = "", pathway: str = "") -> None:
        """Update extracted context parameters."""
        if species:
            self.species = species
        if metabolite:
            self.target_metabolite = metabolite
        if pathway:
            self.target_pathway = pathway
        self._persist()

    # ── Read-only helpers ─────────────────────────────────

    def get_context_summary(self) -> str:
        """Summarize current session context for the LLM."""
        parts = []
        if self.species:
            parts.append(f"Species: {self.species}")
        if self.target_metabolite:
            parts.append(f"Target metabolite: {self.target_metabolite}")
        if self.target_pathway:
            parts.append(f"Target pathway: {self.target_pathway}")
        parts.append(f"Data uploaded: {'Yes' if self.has_data else 'No'}")
        if self.pipeline_result:
            parts.append("Pipeline has been run. Results are cached.")
        return "\n".join(parts)

    def to_summary(self) -> dict:
        return {
            "session_id": self.session_id,
            "turns": self.turn_count,
            "species": self.species,
            "target_metabolite": self.target_metabolite,
            "has_data": self.has_data,
            "has_results": bool(self.pipeline_result),
            "sample_metadata": self.sample_metadata,
            "sample_alignment_report": self.sample_alignment_report,
        }

    # ── Internal ──────────────────────────────────────────

    def _persist(self) -> None:
        """Notify the store to persist this session to DB (if attached)."""
        if self._store is not None:
            self._store._save_session(self)


# ═══════════════════════════════════════════════════════════════
# ConversationStore (SQLite-backed, write-through cache)
# ═══════════════════════════════════════════════════════════════

class ConversationStore:
    """Thread-safe session store backed by SQLite.

    Write-through cache: every mutation to a SessionData immediately
    writes to SQLite. On init, existing sessions are loaded from DB.

    Usage:
        store = ConversationStore()                 # default path
        store = ConversationStore(":memory:")       # in-memory only (legacy behavior)
        store = ConversationStore("/path/to/custom.db")
    """

    def __init__(self, db_path: str | Path | None = None) -> None:
        if db_path is None:
            db_path = str(_default_db_path())

        self._db_path = str(db_path)
        self._sessions: dict[str, SessionData] = {}
        self._lock = threading.Lock()
        self._persistent_conn: sqlite3.Connection | None = None
        if self._db_path == ":memory:":
            self._persistent_conn = sqlite3.connect(self._db_path, check_same_thread=False)
            self._persistent_conn.execute("PRAGMA foreign_keys=ON")
            self._persistent_conn.row_factory = sqlite3.Row

        # Initialize DB
        self._init_db()

        # Load existing sessions from DB (skip for :memory:)
        if self._db_path != ":memory:":
            self._load_all_sessions()

        logger.info(
            "ConversationStore ready: db=%s, loaded=%d sessions",
            self._db_path, len(self._sessions),
        )

    # ═══════════════════════════════════════════════════════════
    # Public API (unchanged)
    # ═══════════════════════════════════════════════════════════

    def get_or_create(self, session_id: str | None = None) -> SessionData:
        """Get existing session or create a new one (thread-safe).

        If session_id exists in memory → return it.
        If session_id exists in DB but not memory → load & return.
        Otherwise → create new, persist to DB.
        """
        # Fast path: already in memory
        if session_id and session_id in self._sessions:
            return self._sessions[session_id]

        with self._lock:
            # Double-check under lock
            if session_id and session_id in self._sessions:
                return self._sessions[session_id]

            # Try loading from DB
            if session_id:
                loaded = self._load_session_from_db(session_id)
                if loaded is not None:
                    self._sessions[session_id] = loaded
                    return loaded

            # Create new
            sid = session_id or f"PA-{uuid.uuid4().hex[:16]}"
            session = SessionData(session_id=sid)
            session._store = self
            self._sessions[sid] = session
            self._insert_session_to_db(session)
            logger.info("Session created: %s", sid)
            return session

    def get(self, session_id: str) -> SessionData | None:
        """Get a session by ID. Checks memory first, then DB."""
        if session_id in self._sessions:
            return self._sessions[session_id]

        with self._lock:
            if session_id in self._sessions:
                return self._sessions[session_id]

            loaded = self._load_session_from_db(session_id)
            if loaded is not None:
                self._sessions[session_id] = loaded
                return loaded

        return None

    def delete(self, session_id: str) -> None:
        """Delete a session from memory and DB."""
        with self._lock:
            self._sessions.pop(session_id, None)
            self._delete_session_from_db(session_id)
        logger.info("Session deleted: %s", session_id)

    def list_sessions(self) -> list[dict]:
        """List all sessions (from DB, for consistent ordering)."""
        return self._list_sessions_from_db()

    # ═══════════════════════════════════════════════════════════
    # Internal: DB initialization
    # ═══════════════════════════════════════════════════════════

    def _init_db(self) -> None:
        """Create tables and indexes if they don't exist."""
        try:
            with self._get_conn() as conn:
                conn.execute(CREATE_SESSIONS_TABLE)
                conn.execute(CREATE_INDEX)
                # Migration: add message_chain column if missing (v4.3 → v4.4)
                try:
                    conn.execute(
                        "ALTER TABLE sessions ADD COLUMN message_chain TEXT DEFAULT '[]'"
                    )
                except sqlite3.OperationalError:
                    pass  # column already exists
                for column, definition in (
                    ("sample_metadata", "TEXT DEFAULT NULL"),
                    ("sample_alignment_report", "TEXT DEFAULT '{}'"),
                ):
                    try:
                        conn.execute(f"ALTER TABLE sessions ADD COLUMN {column} {definition}")
                    except sqlite3.OperationalError:
                        pass  # column already exists
                conn.commit()
        except sqlite3.Error as e:
            logger.error("Failed to initialize DB: %s", e)
            raise

    def _get_conn(self) -> sqlite3.Connection:
        """Get a new SQLite connection (thread-safe, WAL mode)."""
        if self._persistent_conn is not None:
            return self._persistent_conn
        conn = sqlite3.connect(self._db_path, check_same_thread=False)
        conn.execute("PRAGMA journal_mode=WAL")
        conn.execute("PRAGMA foreign_keys=ON")
        conn.row_factory = sqlite3.Row
        return conn

    # ═══════════════════════════════════════════════════════════
    # Internal: session persistence
    # ═══════════════════════════════════════════════════════════

    def _save_session(self, session: SessionData) -> None:
        """Persist a session to DB (UPSERT). Called by SessionData._persist()."""
        try:
            with self._get_conn() as conn:
                conn.execute("""
                    INSERT OR REPLACE INTO sessions (
                        session_id, created_at, species, target_metabolite,
                        target_pathway, has_data, turn_count, history,
                        message_chain,
                        expression_matrix, metabolite_matrix, promoter_sequences,
                        sample_metadata, sample_alignment_report,
                        pipeline_result, pipeline_hypotheses,
                        pipeline_mechanism_type, pipeline_pathway
                    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """, (
                    session.session_id,
                    session.created_at,
                    session.species,
                    session.target_metabolite,
                    session.target_pathway,
                    1 if session.has_data else 0,
                    session.turn_count,
                    json.dumps(session.history, ensure_ascii=False),
                    json.dumps(session.message_chain, ensure_ascii=False),
                    json.dumps(session.expression_matrix, ensure_ascii=False) if session.expression_matrix is not None else None,
                    json.dumps(session.metabolite_matrix, ensure_ascii=False) if session.metabolite_matrix is not None else None,
                    json.dumps(session.promoter_sequences, ensure_ascii=False) if session.promoter_sequences is not None else None,
                    json.dumps(session.sample_metadata, ensure_ascii=False) if session.sample_metadata is not None else None,
                    json.dumps(session.sample_alignment_report, ensure_ascii=False),
                    session.pipeline_result,
                    json.dumps(session.pipeline_hypotheses, ensure_ascii=False),
                    session.pipeline_mechanism_type,
                    session.pipeline_pathway,
                ))
                conn.commit()
        except sqlite3.Error as e:
            logger.error("Failed to save session %s: %s", session.session_id, e)

    def _insert_session_to_db(self, session: SessionData) -> None:
        """Insert a brand-new session (only called during creation)."""
        self._save_session(session)

    def _load_session_from_db(self, session_id: str) -> SessionData | None:
        """Load a single session from DB. Returns None if not found."""
        try:
            with self._get_conn() as conn:
                row = conn.execute(
                    "SELECT * FROM sessions WHERE session_id = ?", (session_id,)
                ).fetchone()

            if row is None:
                return None

            return self._row_to_session(row)

        except sqlite3.Error as e:
            logger.error("Failed to load session %s: %s", session_id, e)
            return None

    def _load_all_sessions(self) -> None:
        """Load all sessions from DB into memory on startup."""
        try:
            with self._get_conn() as conn:
                rows = conn.execute(
                    "SELECT * FROM sessions ORDER BY created_at DESC"
                ).fetchall()

            for row in rows:
                session = self._row_to_session(row)
                self._sessions[session.session_id] = session

            if rows:
                logger.info("Loaded %d sessions from %s", len(rows), self._db_path)

        except sqlite3.Error as e:
            logger.error("Failed to load sessions: %s", e)

    def _delete_session_from_db(self, session_id: str) -> None:
        """Remove a session from DB."""
        try:
            with self._get_conn() as conn:
                conn.execute(
                    "DELETE FROM sessions WHERE session_id = ?", (session_id,)
                )
                conn.commit()
        except sqlite3.Error as e:
            logger.error("Failed to delete session %s: %s", session_id, e)

    def _list_sessions_from_db(self) -> list[dict]:
        """List all sessions from DB."""
        try:
            with self._get_conn() as conn:
                rows = conn.execute(
                    "SELECT session_id, created_at, turn_count, has_data, species, target_metabolite "
                    "FROM sessions ORDER BY created_at DESC"
                ).fetchall()

            return [
                {
                    "session_id": r["session_id"],
                    "created_at": r["created_at"],
                    "turns": r["turn_count"],
                    "has_data": bool(r["has_data"]),
                    "species": r["species"],
                    "target_metabolite": r["target_metabolite"],
                }
                for r in rows
            ]
        except sqlite3.Error as e:
            logger.error("Failed to list sessions: %s", e)
            return []

    # ═══════════════════════════════════════════════════════════
    # Internal: row → SessionData
    # ═══════════════════════════════════════════════════════════

    @staticmethod
    def _safe_json_load(text: str | None, default=None):
        """Safely parse JSON, returning default on failure."""
        if not text:
            return default
        try:
            return json.loads(text)
        except (json.JSONDecodeError, TypeError):
            return default

    def _row_to_session(self, row: sqlite3.Row) -> SessionData:
        """Convert a DB row to a SessionData instance."""
        session = SessionData(
            session_id=row["session_id"],
            created_at=row["created_at"],
            species=row["species"] or "",
            target_metabolite=row["target_metabolite"] or "",
            target_pathway=row["target_pathway"] or "",
            has_data=bool(row["has_data"]),
            turn_count=row["turn_count"] or 0,
            history=self._safe_json_load(row["history"], []),
            expression_matrix=self._safe_json_load(row["expression_matrix"]),
            metabolite_matrix=self._safe_json_load(row["metabolite_matrix"]),
            promoter_sequences=self._safe_json_load(row["promoter_sequences"]),
            sample_metadata=self._safe_json_load(row["sample_metadata"]),
            sample_alignment_report=self._safe_json_load(row["sample_alignment_report"], {}),
            pipeline_result=row["pipeline_result"] or "",
            pipeline_hypotheses=self._safe_json_load(row["pipeline_hypotheses"], []),
            pipeline_mechanism_type=row["pipeline_mechanism_type"] or "",
            pipeline_pathway=row["pipeline_pathway"] or "",
            message_chain=self._safe_json_load(row["message_chain"], []),
        )
        # Recompute this flag for sessions written by older versions, where
        # only expression data counted as "uploaded".
        session.has_data = any((
            bool(session.expression_matrix),
            bool(session.metabolite_matrix),
            bool(session.promoter_sequences),
        ))
        # Attach store reference for future persistence
        session._store = self
        return session


# ── Global store instance ──────────────────────────────────────

store = ConversationStore()
