"""
OmiCraft ngrok 런처 — agent_rev9
NGROK_DOMAIN으로 고정 도메인을 지정할 수 있습니다.

사용법:
  # NGROK_AUTHTOKEN 환경 변수로 실행
  set NGROK_AUTHTOKEN=<your-token>
  python run_ngrok.py

  # 또는 토큰 직접 지정
  python run_ngrok.py --token <your-token>

환경 변수:
  NGROK_AUTHTOKEN   ngrok 인증 토큰 (https://dashboard.ngrok.com/get-started/your-authtoken)
  OMICRAFT_PORT     로컬 포트 (기본 8000)
  OMICRAFT_CONFIG   설정 파일 경로
"""

from __future__ import annotations

import argparse
import os
import sys
import threading
import time
from pathlib import Path

# ── 패키지 경로 ──────────────────────────────────────────────
ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parent.parent))  # src/ 경로 추가

NGROK_DOMAIN = os.environ.get("NGROK_DOMAIN")
DEFAULT_PORT = 8000


# ── 서버 스레드 ──────────────────────────────────────────────

def _load_package():
    if __package__:
        from .launch import load_package
    else:
        from launch import load_package
    load_package()


def _start_uvicorn(port: int):
    import uvicorn
    _load_package()
    uvicorn.run(
        "agent_rev9.server:app",
        host="127.0.0.1",
        port=port,
        proxy_headers=True,       # ngrok X-Forwarded-* 헤더 수용
        forwarded_allow_ips="*",  # ngrok IP 신뢰
        log_level="info",
    )


# ── 메인 ────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="OmiCraft ngrok 런처")
    parser.add_argument("--token", help="ngrok 인증 토큰 (NGROK_AUTHTOKEN 환경 변수 대체)")
    parser.add_argument("--port", type=int, default=int(os.environ.get("OMICRAFT_PORT", DEFAULT_PORT)),
                        help=f"로컬 포트 (기본 {DEFAULT_PORT})")
    parser.add_argument("--domain", default=NGROK_DOMAIN,
                        help=f"ngrok 정적 도메인 (기본 {NGROK_DOMAIN})")
    args = parser.parse_args()

    port = args.port
    domain = args.domain

    # 설정 파일 기본값
    if (ROOT / "configs/upstream.local.json").is_file():
        os.environ.setdefault("OMICRAFT_CONFIG", str(ROOT / "configs/upstream.local.json"))

    # ngrok 인증 토큰
    auth_token = (args.token or os.environ.get("NGROK_AUTHTOKEN", "")).strip()

    from pyngrok import conf as pyngrok_conf, ngrok

    if auth_token:
        pyngrok_conf.get_default().auth_token = auth_token
    else:
        # 이미 저장된 토큰이 있는지 확인
        try:
            ngrok.get_ngrok_process()
        except Exception:
            pass
        saved = pyngrok_conf.get_default().auth_token
        if not saved:
            print("=" * 60)
            print("❌  NGROK_AUTHTOKEN 이 설정되지 않았습니다.")
            print()
            print("  1. https://dashboard.ngrok.com/get-started/your-authtoken")
            print("     에서 토큰을 복사하세요.")
            print()
            print("  2. 아래 중 하나로 실행하세요:")
            print("     set NGROK_AUTHTOKEN=<token>")
            print("     python run_ngrok.py --token <token>")
            print("=" * 60)
            sys.exit(1)

    # ── 1. FastAPI 서버 시작 ──────────────────────────────
    print("=" * 60)
    print("  OmiCraft Agent Dashboard")
    print("=" * 60)
    print(f"\n[1/3] FastAPI 서버 시작 중 (port={port}) ...")
    server_thread = threading.Thread(target=_start_uvicorn, args=(port,), daemon=True)
    server_thread.start()

    # 서버 준비 대기 (최대 10초)
    import urllib.request
    for _ in range(20):
        try:
            urllib.request.urlopen(f"http://127.0.0.1:{port}/health", timeout=1)
            break
        except Exception:
            time.sleep(0.5)

    # ── 2. ngrok 터널 연결 ───────────────────────────────
    print(f"[2/3] ngrok 터널 연결 중 → {domain} ...")
    try:
        tunnel = ngrok.connect(
            addr=str(port),
            proto="http",
            **({"hostname": domain} if domain else {}),
        )
        public_url = tunnel.public_url
        # http → https 정규화
        if public_url.startswith("http://"):
            public_url = "https://" + public_url[7:]
    except Exception as exc:
        print(f"\n❌  ngrok 터널 연결 실패: {exc}")
        print()
        print("  확인 사항:")
        print("  · 인증 토큰이 올바른지 확인하세요")
        print(f"  · 도메인 '{domain}' 이 내 계정에 등록되어 있는지 확인하세요")
        print("    → https://dashboard.ngrok.com/domains")
        sys.exit(1)

    # ── 3. 완료 ──────────────────────────────────────────
    print(f"[3/3] 터널 연결 완료!\n")
    print("  ┌─────────────────────────────────────────────────┐")
    print(f"  │  공개 URL  : {public_url:<34} │")
    print(f"  │  로컬 URL  : http://127.0.0.1:{port:<26} │")
    print("  └─────────────────────────────────────────────────┘")
    print()
    print("  처음 접속 시 ngrok 경고 페이지가 표시될 수 있습니다.")
    print("  'Visit Site' 를 클릭하면 대시보드로 이동합니다.")
    print()
    print("  Ctrl+C 로 종료\n")

    try:
        while True:
            time.sleep(1)
    except KeyboardInterrupt:
        print("\n종료 중...")
        ngrok.kill()
        print("완료.")


if __name__ == "__main__":
    main()
