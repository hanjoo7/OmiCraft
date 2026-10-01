# OmiCraft rev9

rev8의 분석·설계 파이프라인을 독립 배포용으로 정리한 버전입니다. rev8 원본은 변경하지 않습니다.

## 기능

- 연구질문 → Planner → DESeq2/GSEA → aPEAR 네트워크 시각화 → 세포 맥락·DepMap → 표적 적격성 → Critic → 선택/설계 → Dossier.
- **Small molecule, De novo binder, ADC, Degrader** 4개 모달리티와 각각의 검증·판정 경로.
- 대시보드, Structure Lab, 3D 구조/리간드 표시, 실행 이력, 보고서 내보내기, 결과 재사용, 작업 중단·시간 제한.
- 직접 실행 CLI, 모의 실행, 고정 시나리오, 쌍둥이 환자 분석, 선택적 ngrok/프록시 실행.

## 설치 및 실행

Python 3.11 이상. 이 디렉터리에서 실행합니다.

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e '.[web,graph,test]'
omicraft-serve backend --port 8088
```

브라우저에서 `http://127.0.0.1:8088`을 엽니다. 소스에서 직접 실행할 때는 `python launch.py backend --port 8088`도 사용할 수 있습니다. 외부 접속은 `--host`로 지정합니다.

GPU 모델을 호출하지 않는 CLI 예제:

```bash
omicraft-screen binder --demo --output-dir ./runs/binder-demo
omicraft-screen small-molecule --input examples/small_molecule/mock.input.json --config examples/small_molecule/mock.config.json --mock --output ./runs/small-molecule-demo
```

`omicraft-screen --help`에서 전체 경로를 확인할 수 있습니다. 실제 설계에는 해당 모델·가중치·데이터와 실행 설정이 필요합니다. 모의 결과는 실제 성공 판정으로 승격하지 않습니다.

## 설정

`configs/*.example.json`을 로컬 설정으로 복사하고 데이터·모델 경로를 수정합니다. 소스 실행에서는 같은 이름의 `*.local.json`을 우선 사용합니다. 설치된 패키지를 수정하지 않으려면 아래 환경변수로 외부 JSON을 지정합니다.

| 환경변수 | 용도 |
| --- | --- |
| `OMICRAFT_WORKSPACE` | 데이터·모델·결과 기준 디렉터리, 기본값은 실행 디렉터리 |
| `OMICRAFT_CONFIG` | 파이프라인 설정 (`upstream.example.json` 참조); `--config`로도 지정 가능 |
| `OMICRAFT_MODELS_CONFIG` | 모달리티별 실행 모델 설정 (`workbench_models.example.json`) |
| `OMICRAFT_ASSETS_CONFIG` | Structure Lab 구조 목록 (`structure_assets.example.json`) |
| `OMICRAFT_DESIGN_INPUTS_CONFIG` | 자동 설계 입력 (`automatic_design_inputs.example.json`) |
| `OMICRAFT_LLM_CONFIG` | LLM 연결 설정; 기본값은 `configs/llm.json` |
| `DACON_API_KEY` | 기본 LLM 설정의 인증 키; 다른 제공자는 LLM JSON의 `api_key_env` 수정 |
| `OMICRAFT_RUNS_ROOT` | 보고서 검색 디렉터리; 파이프라인 설정의 `data.base_dir`와 같은 경로 권장 |
| `OMICRAFT_REFERENCE_RUN` | 기존 분석 결과를 재사용할 디렉터리 |

예제 JSON의 `${OMICRAFT_WORKSPACE}`와 `${OMICRAFT_PACKAGE}`는 설정을 읽을 때 확장됩니다. 예제에 포함된 외부 모델/데이터 경로는 사용 환경에 맞게 수정해야 합니다. R 분석 환경은 `configs/environment-omics-r.yml`, Python 분석 의존성은 `requirements-upstream.txt`를 참고합니다. RFdiffusion3, ProteinMPNN, AlphaFold3, 도킹 도구 및 모델 가중치는 별도로 준비합니다. MSA 생성·캐시 재사용 기능은 유지하며 대용량 캐시 자체는 배포하지 않습니다.

프록시가 필요하면 `OMICRAFT_UPSTREAM`, `OMICRAFT_PUBLIC_RUN`을 지정하고 `omicraft-serve proxy --port 8089`를 실행합니다. 공개 실행 허용 설정은 해당 run의 `cloudflare_access/execution_policy.json`에서 관리합니다. ngrok은 `pip install '.[ngrok]'` 후 `python run_ngrok.py`로 실행하며 `NGROK_AUTHTOKEN`, 선택적으로 `NGROK_DOMAIN`을 지정합니다.

## 검증 및 배포

```bash
pytest
node tests/test_dashboard_progress.cjs
pip wheel --no-deps . -w dist
```

기본 테스트는 실제 모델 추론·외부 네트워크를 차단합니다. `local_io` 테스트는 로컬 HTTP/CPU 작업 프로세스를 검증합니다. `external_tool` 테스트는 설치된 도구를 명시적으로 사용하는 선택 검증이며 기본 실행에서 제외됩니다.

과거 설계 문서·중복본·미참조 모듈·ZIP·실행 캐시를 제거하고 서버 화면을 `templates/`로 분리했습니다. 테스트, 설정 예제, R 스크립트, 구조 예제, 의존성 목록과 3Dmol 라이선스는 유지합니다. 개인 키, `*.local.json`, 실행 결과, 모델 가중치는 `.gitignore`로 제외합니다.
