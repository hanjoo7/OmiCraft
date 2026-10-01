"""
⑩ Artifact Contract — run_id + sha256 검증
워크플로우 Section 4: 실행 계약

각 노드 실행 시:
  - run_id 부여
  - 입력 artifact의 sha256 확인
  - 출력 artifact 기록
  - 상태(PENDING→RUNNING→COMPLETED/FAILED)
"""

import os
import hashlib
import json
import time
from datetime import datetime, timezone


def sha256_file(path: str) -> str:
    """파일 SHA256 해시."""
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


def sha256_dict(d: dict) -> str:
    """Dict → JSON → SHA256."""
    return hashlib.sha256(json.dumps(d, sort_keys=True, default=str).encode()).hexdigest()


class ArtifactContract:
    """실행 계약 관리."""

    def __init__(self, run_id: str = None):
        self.run_id = run_id or datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S")
        self.nodes = []
        self.start_time = time.time()

    def begin_node(self, node_id: str, inputs: list = None, parameters: dict = None,
                   tool: str = None) -> dict:
        """노드 실행 시작."""
        contract = {
            "run_id": self.run_id,
            "node_id": node_id,
            "status": "RUNNING",
            "started_at": datetime.now(timezone.utc).isoformat(),
            "inputs": [],
            "parameters": parameters or {},
            "tool": tool or "",
            "outputs": [],
            "evidence_ids": [],
        }

        # 입력 파일 해시
        if inputs:
            for inp in inputs:
                if os.path.exists(inp):
                    contract["inputs"].append({
                        "path": inp,
                        "sha256": sha256_file(inp),
                        "size": os.path.getsize(inp),
                    })

        self.nodes.append(contract)
        return contract

    def complete_node(self, node_id: str, outputs: list = None, evidence_ids: list = None,
                      status: str = "COMPLETED") -> dict:
        """노드 실행 완료."""
        for node in self.nodes:
            if node["node_id"] == node_id and node["status"] == "RUNNING":
                node["status"] = status
                node["completed_at"] = datetime.now(timezone.utc).isoformat()

                if outputs:
                    for out in outputs:
                        if os.path.exists(out):
                            node["outputs"].append({
                                "path": out,
                                "sha256": sha256_file(out),
                                "size": os.path.getsize(out),
                            })

                if evidence_ids:
                    node["evidence_ids"] = evidence_ids

                return node
        return None

    def verify_input(self, node_id: str, expected_hash: str, actual_path: str) -> bool:
        """입력 파일 해시 검증."""
        if not os.path.exists(actual_path):
            return False
        return sha256_file(actual_path) == expected_hash

    def save(self, output_dir: str):
        """계약 저장."""
        path = os.path.join(output_dir, f"contract_{self.run_id}.json")
        data = {
            "run_id": self.run_id,
            "total_time_sec": round(time.time() - self.start_time, 1),
            "nodes": self.nodes,
        }
        with open(path, "w", encoding="utf-8") as f:
            json.dump(data, f, indent=2, ensure_ascii=False, default=str)
        return path


if __name__ == "__main__":
    import sys
    if sys.stdout.encoding != "utf-8":
        sys.stdout.reconfigure(encoding="utf-8")

    contract = ArtifactContract()
    print(f"Run ID: {contract.run_id}")

    # 테스트: DEG 결과 파일 해시
    deg_path = os.path.join(
        os.path.dirname(__file__), "..", "..", "..", "dataset", "tcga_brca", "processed",
        "deg_results", "deg_tnbc_vs_nontnbc_full.csv"
    )
    if os.path.exists(deg_path):
        contract.begin_node("discovery", inputs=[deg_path], tool="DESeq2")
        h = sha256_file(deg_path)
        print(f"DEG file hash: {h[:16]}...")
        contract.complete_node("discovery", outputs=[deg_path])

    out_dir = os.path.join(
        os.path.dirname(__file__), "..", "..", "..", "dataset", "tcga_brca", "processed", "agent_results"
    )
    os.makedirs(out_dir, exist_ok=True)
    path = contract.save(out_dir)
    print(f"Contract saved: {path}")
    print(f"Nodes: {len(contract.nodes)}")
