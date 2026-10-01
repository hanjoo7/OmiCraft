from __future__ import annotations
"""
Evidence Store — SQLite 영속화
모든 Evidence Card를 SQLite DB에 저장하여 세션 간 재사용·추적이 가능하다.

제안서: "근거·절차·실패·에피소드 메모리를 분리 저장하여
반복 실행에서 동일 후보의 탈락 사유를 재계산하지 않으며,
각 항목은 사실이 아닌 가설로 취급해 Critic이 무효화할 수 있고
근거 DB 버전이 갱신되면 만료·재평가"
"""

import os
import sqlite3
import time
from typing import Optional

DB_PATH = os.path.join("runs", "evidence_store.db")


def get_connection(db_path: str = None) -> sqlite3.Connection:
    path = db_path or DB_PATH
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    conn = sqlite3.connect(path)
    conn.row_factory = sqlite3.Row
    _init_tables(conn)
    return conn


def _init_tables(conn: sqlite3.Connection):
    conn.executescript("""
    CREATE TABLE IF NOT EXISTS evidence_cards (
        id TEXT PRIMARY KEY,
        source TEXT NOT NULL,
        type TEXT NOT NULL CHECK(type IN ('SUPPORTIVE','CONTRADICTORY','NOT_APPLICABLE','UNKNOWN')),
        content TEXT,
        gene TEXT,
        confidence REAL DEFAULT 0.0,
        run_id TEXT,
        agent TEXT,
        created_at REAL,
        invalidated INTEGER DEFAULT 0,
        invalidated_by TEXT
    );

    CREATE TABLE IF NOT EXISTS run_history (
        run_id TEXT PRIMARY KEY,
        question TEXT,
        started_at REAL,
        finished_at REAL,
        final_verdict TEXT,
        n_advance INTEGER,
        n_reject INTEGER,
        n_hold INTEGER,
        node_path TEXT
    );

    CREATE TABLE IF NOT EXISTS failure_memory (
        gene TEXT NOT NULL,
        reason TEXT,
        agent TEXT,
        run_id TEXT,
        created_at REAL,
        recheck_after TEXT,
        PRIMARY KEY (gene, reason)
    );

    CREATE INDEX IF NOT EXISTS idx_ev_gene ON evidence_cards(gene);
    CREATE INDEX IF NOT EXISTS idx_ev_run ON evidence_cards(run_id);
    CREATE INDEX IF NOT EXISTS idx_fail_gene ON failure_memory(gene);
    """)
    conn.commit()


class EvidenceStore:
    def __init__(self, db_path: str = None):
        self.conn = get_connection(db_path)

    def add_evidence(self, card: dict, run_id: str = "", agent: str = ""):
        self.conn.execute(
            "INSERT OR REPLACE INTO evidence_cards (id,source,type,content,gene,confidence,run_id,agent,created_at) "
            "VALUES (?,?,?,?,?,?,?,?,?)",
            (
                card.get("id", ""),
                card.get("source", ""),
                card.get("type", "UNKNOWN"),
                card.get("content", ""),
                card.get("gene", ""),
                card.get("confidence", 0),
                run_id,
                agent,
                time.time(),
            ),
        )
        self.conn.commit()

    def add_many(self, cards: list, run_id: str = "", agent: str = ""):
        for c in cards:
            self.add_evidence(c, run_id, agent)

    def get_by_gene(self, gene: str) -> list:
        rows = self.conn.execute(
            "SELECT * FROM evidence_cards WHERE gene=? AND invalidated=0 ORDER BY created_at DESC",
            (gene,),
        ).fetchall()
        return [dict(r) for r in rows]

    def get_by_run(self, run_id: str) -> list:
        rows = self.conn.execute(
            "SELECT * FROM evidence_cards WHERE run_id=? ORDER BY created_at", (run_id,)
        ).fetchall()
        return [dict(r) for r in rows]

    def invalidate(self, evidence_id: str, reason: str = ""):
        self.conn.execute(
            "UPDATE evidence_cards SET invalidated=1, invalidated_by=? WHERE id=?",
            (reason, evidence_id),
        )
        self.conn.commit()

    def save_run(
        self,
        run_id: str,
        question: str,
        verdict: str,
        n_advance: int = 0,
        n_reject: int = 0,
        n_hold: int = 0,
        node_path: str = "",
    ):
        self.conn.execute(
            "INSERT OR REPLACE INTO run_history VALUES (?,?,?,?,?,?,?,?,?)",
            (
                run_id,
                question,
                time.time(),
                time.time(),
                verdict,
                n_advance,
                n_reject,
                n_hold,
                node_path,
            ),
        )
        self.conn.commit()

    def add_failure(self, gene: str, reason: str, agent: str = "", run_id: str = ""):
        self.conn.execute(
            "INSERT OR REPLACE INTO failure_memory VALUES (?,?,?,?,?,?)",
            (gene, reason, agent, run_id, time.time(), None),
        )
        self.conn.commit()

    def is_known_failure(self, gene: str) -> Optional[dict]:
        row = self.conn.execute(
            "SELECT * FROM failure_memory WHERE gene=? ORDER BY created_at DESC LIMIT 1", (gene,)
        ).fetchone()
        return dict(row) if row else None

    def get_all_runs(self) -> list:
        rows = self.conn.execute("SELECT * FROM run_history ORDER BY started_at DESC").fetchall()
        return [dict(r) for r in rows]

    def stats(self) -> dict:
        ev_count = self.conn.execute(
            "SELECT COUNT(*) FROM evidence_cards WHERE invalidated=0"
        ).fetchone()[0]
        run_count = self.conn.execute("SELECT COUNT(*) FROM run_history").fetchone()[0]
        fail_count = self.conn.execute("SELECT COUNT(*) FROM failure_memory").fetchone()[0]
        return {"evidence_cards": ev_count, "runs": run_count, "failures": fail_count}

    def close(self):
        self.conn.close()


if __name__ == "__main__":
    store = EvidenceStore()
    print("Evidence Store initialized")
    print(f"Stats: {store.stats()}")

    # 테스트: Evidence Card 저장
    store.add_evidence(
        {
            "id": "EV-TEST-001",
            "source": "test",
            "type": "SUPPORTIVE",
            "content": "Test evidence card",
            "gene": "PARP1",
            "confidence": 0.9,
        },
        run_id="test_run",
        agent="test",
    )

    store.add_failure("RRM2", "Safety veto — common essential", agent="qualification")

    print(f"After insert: {store.stats()}")
    print(f"PARP1 evidence: {store.get_by_gene('PARP1')}")
    print(f"RRM2 failure: {store.is_known_failure('RRM2')}")

    store.close()
    print("Done")
