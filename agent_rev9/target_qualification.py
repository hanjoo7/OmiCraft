"""
Target Qualification 파이프라인 — 워크플로우 v5 Section 06–07

UniProt / HPA / Open Targets / ChEMBL 조회 →
Evidence Card 조립 → Biological confidence (B) →
Shortlist 규칙 → hypotheses.tsv

주의사항 (워크플로우 §06–07):
  - database association score / 생존 상관관계를 인과성으로 바꾸지 않는다.
  - 정상 RNA 고발현 하나만으로 모든 modality를 일괄 REJECT하지 않는다.
  - critical safety rule이 충족될 때만 veto 적용; 불명확한 위험은 HOLD.
  - veto는 Agent가 점수로 상쇄할 수 없다.
  - intervention_role은 FUNCTION_MODULATION / CELL_TARGETING / BIOMARKER_ONLY.
"""

from __future__ import annotations

import csv
import hashlib
import json
import tempfile
import time
import urllib.error
import urllib.parse
import urllib.request
from functools import lru_cache
from pathlib import Path

from .tier_system import (
    assign_bio_confidence,
)

# ── Evidence 상태 열거 ─────────────────────────────────────

EvidenceStatus = str  # SUPPORTIVE | CONTRADICTORY | NOT_APPLICABLE | UNKNOWN

INTERVENTION_ROLES = {
    "FUNCTION_MODULATION": "기능 억제/활성화",
    "CELL_TARGETING":      "세포 표적화 (ADC/CAR-T 등)",
    "BIOMARKER_ONLY":      "진단/바이오마커 전용",
}

# ── API 헬퍼 ──────────────────────────────────────────────

_TIMEOUT = 10  # 초


def _http_get_json(url: str, headers: dict | None = None, payload: bytes | None = None) -> dict | list | None:
    from .configuration import get_config
    from .resource_usage import record_usage
    cache_root = get_config().upstream.cache_dir
    cache = Path(cache_root) / "api" / (hashlib.sha256(url.encode() + (payload or b"")).hexdigest() + ".json") if cache_root else None
    if cache and cache.is_file() and time.time() - cache.stat().st_mtime < 7 * 86400:
        try:
            data = json.loads(cache.read_text())["response"]
            record_usage("cache_hit", tool="qualification API")
            return data
        except (OSError, ValueError, KeyError, TypeError):
            pass
    req = urllib.request.Request(url, data=payload, headers=headers or {
        "Accept": "application/json",
        "User-Agent": "OmiCraft/5.0",
    })
    try:
        record_usage("tool_call", tool="qualification API")
        with urllib.request.urlopen(req, timeout=_TIMEOUT) as resp:
            data = json.loads(resp.read().decode("utf-8"))
            if cache:
                cache.parent.mkdir(parents=True, exist_ok=True)
                temporary = None
                try:
                    with tempfile.NamedTemporaryFile(mode="w", dir=cache.parent, delete=False, encoding="utf-8") as handle:
                        temporary = Path(handle.name)
                        json.dump({"url": url, "request_body": payload.decode() if payload else None,
                                   "retrieved_at": time.time(), "response": data}, handle)
                    temporary.replace(cache)
                finally:
                    if temporary:
                        temporary.unlink(missing_ok=True)
            return data
    except (urllib.error.URLError, urllib.error.HTTPError, json.JSONDecodeError,
            OSError):
        return None


# ── UniProt 조회 ──────────────────────────────────────────

def query_uniprot(gene: str, organism: str = "9606") -> dict:
    """UniProt에서 단백질 위치·isoform·기능을 조회한다.

    Returns:
        {
          "status": str,
          "uniprot_id": str | None,
          "protein_name": str,
          "gene_names": list[str],
          "subcellular_location": list[str],
          "function": str,
          "n_isoforms": int,
          "is_secreted": bool,
          "has_transmembrane": bool,
          "source": "UniProt API",
        }
    """
    url = (
        f"https://rest.uniprot.org/uniprotkb/search"
        f"?query=gene_exact:{urllib.parse.quote(gene)}+AND+organism_id:{organism}"
        f"+AND+reviewed:true"
        f"&fields=id,gene_names,protein_name,cc_subcellular_location,ft_transmem,xref_pdb,"
        f"cc_function,cc_alternative_products"
        f"&format=json&size=1"
    )
    data = _http_get_json(url)
    if not data or not data.get("results"):
        return {"status": "UNAVAILABLE" if data is None else "NOT_FOUND", "uniprot_id": None, "gene": gene,
                "subcellular_location": [], "is_secreted": False,
                "has_transmembrane": False, "n_isoforms": 0,
                "function": "", "protein_name": ""}

    entry = data["results"][0]
    uid   = entry.get("primaryAccession", "")

    # 위치 파싱
    loc_list: list[str] = []
    for loc in (entry.get("comments") or []):
        if loc.get("commentType") == "SUBCELLULAR LOCATION":
            for subloc in loc.get("subcellularLocations") or []:
                loc_list.append(subloc.get("location", {}).get("value", ""))

    is_secreted     = any("secreted" in location.lower() for location in loc_list)
    has_transmem    = bool(entry.get("features") and
                           any(f.get("type") == "Transmembrane"
                               for f in entry["features"]))

    # 기능 텍스트
    func_text = ""
    for c in (entry.get("comments") or []):
        if c.get("commentType") == "FUNCTION":
            texts = c.get("texts") or []
            func_text = " ".join(t.get("value", "") for t in texts)[:500]
            break

    # isoform 수
    n_isoforms = 0
    for c in (entry.get("comments") or []):
        if c.get("commentType") == "ALTERNATIVE SEQUENCE":
            n_isoforms += 1

    # gene_names
    gene_obj  = entry.get("genes") or []
    gene_names = []
    for g in gene_obj:
        gn = g.get("geneName", {}).get("value")
        if gn:
            gene_names.append(gn)

    return {
        "status": "FOUND",
        "gene": gene,
        "uniprot_id": uid,
        "protein_name": (entry.get("proteinDescription") or {})
                        .get("recommendedName", {}).get("fullName", {}).get("value", ""),
        "gene_names": gene_names,
        "subcellular_location": [location for location in loc_list if location],
        "is_secreted": is_secreted,
        "has_transmembrane": has_transmem,
        "n_isoforms": n_isoforms,
        "function": func_text,
        "source": "UniProt REST API",
        "pdb_ids": [x["id"] for x in entry.get("uniProtKBCrossReferences", []) if x.get("database") == "PDB"],
    }


def _location_to_target_info(uniprot_result: dict) -> dict:
    """UniProt 결과에서 target_info (tier_system 입력) 를 생성."""
    locs = [location.lower() for location in uniprot_result.get("subcellular_location", [])]
    is_surface = (
        any(kw in " ".join(locs) for kw in
            ("cell membrane", "plasma membrane", "extracellular", "cell surface"))
    )
    is_intracellular = any(kw in " ".join(locs) for kw in
                           ("nucleus", "cytoplasm", "mitochondr", "endoplasmic"))
    return {
        "location": "SURFACE" if is_surface else
                    "INTRACELLULAR" if is_intracellular else
                    "SECRETED" if uniprot_result.get("is_secreted") else "UNKNOWN",
        "surface_accessible": is_surface,
        "has_transmembrane": uniprot_result.get("has_transmembrane", False),
        "is_secreted": uniprot_result.get("is_secreted", False),
        "has_structure": False,  # PDB 조회에서 채움
        "has_pocket": False,     # 구조 분석에서 채움
    }


# ── HPA 조회 ──────────────────────────────────────────────

@lru_cache(maxsize=4)
def _hpa_index(path, mtime):
    data = {}
    with Path(path).open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            for key in {row.get("Gene"), row.get("Gene name")} - {None, ""}:
                data.setdefault(key, []).append(row)
    return data


def query_hpa_normal(gene: str, gene_id: str | None = None) -> dict:
    from .configuration import get_config
    path = get_config().data.hpa_path
    base = {"status": "UNAVAILABLE", "gene": gene, "normal_tissue": [],
            "high_expression_tissues": [], "critical_tissue_expression": False,
            "source": "Human Protein Atlas 25.1 normal IHC", "source_file": path}
    if not path or not Path(path).is_file():
        return {**base, "reason": "Normal IHC snapshot not configured"}
    index = _hpa_index(path, Path(path).stat().st_mtime_ns)
    rows = index.get((gene_id or "").split(".")[0]) or index.get(gene, [])
    if not rows:
        return {**base, "status": "NOT_FOUND"}
    tissues = [{"tissue": r.get("Tissue", ""), "cell_type": r.get("Cell type", ""),
                "level": r.get("Level", ""), "reliability": r.get("Reliability", "")} for r in rows]
    high = sorted({r["tissue"] for r in tissues if r["level"] in {"High", "Medium"}
                   and r["reliability"] in {"Enhanced", "Supported", "Approved"}})
    critical = {"heart muscle", "cerebral cortex", "cerebellum", "kidney", "liver", "bone marrow"}
    return {**base, "status": "FOUND", "normal_tissue": tissues, "high_expression_tissues": high,
            "critical_tissue_expression": any(t.lower() in critical for t in high)}


# ── Open Targets 조회 ────────────────────────────────────

def query_open_targets(gene_ensembl_id: str, disease_id: str = "MONDO_0004989") -> dict:
    """Open Targets Platform에서 유전자-질환 근거 점수를 조회한다.

    disease_id: EFO_0000305 = breast carcinoma

    Returns:
        {
          "status": str,
          "gene": str,
          "disease_id": str,
          "overall_score": float | None,
          "genetic_association_score": float | None,
          "somatic_mutation_score": float | None,
          "literature_score": float | None,
          "source": "OpenTargets GraphQL",
        }
    """
    query = """
    query TargetDiseaseEvidence($ensgId: String!, $efoId: String!) {
      target(ensemblId: $ensgId) {
        id
        approvedSymbol
      }
      disease(efoId: $efoId) {
        associatedTargets(Bs: [$ensgId]) {
          rows {
            score
            datatypeScores { id score }
          }
        }
      }
    }
    """
    url = "https://api.platform.opentargets.org/api/v4/graphql"
    payload = json.dumps({
        "query": query,
        "variables": {"ensgId": gene_ensembl_id, "efoId": disease_id},
    }).encode("utf-8")
    data = _http_get_json(url, headers={"Content-Type": "application/json", "Accept": "application/json"}, payload=payload)
    if not data or data.get("errors"):
        return {"status": "UNAVAILABLE", "gene": gene_ensembl_id,
                "disease_id": disease_id,
                "overall_score": None, "genetic_association_score": None,
                "somatic_mutation_score": None, "literature_score": None,
                "source": "OpenTargets GraphQL"}

    rows = ((data.get("data") or {}).get("disease") or {}).get("associatedTargets", {}).get("rows") or []
    if not rows:
        return {"status": "NOT_FOUND", "gene": gene_ensembl_id,
                "disease_id": disease_id,
                "overall_score": None, "genetic_association_score": None,
                "somatic_mutation_score": None, "literature_score": None,
                "source": "OpenTargets GraphQL"}

    row = rows[0]
    overall = row.get("score")
    scores = {s["id"]: s["score"]
              for s in (row.get("datatypeScores") or [])}

    return {
        "status": "FOUND",
        "gene": gene_ensembl_id,
        "disease_id": disease_id,
        "overall_score": overall,
        "genetic_association_score": scores.get("genetic_association"),
        "somatic_mutation_score": scores.get("somatic_mutation"),
        "literature_score": scores.get("literature"),
        "source": "OpenTargets GraphQL",
        "note": "association score는 인과성이 아님",
    }


# ── ChEMBL 조회 ──────────────────────────────────────────

def query_chembl(gene: str) -> dict:
    """ChEMBL에서 해당 유전자 표적의 임상 자산을 조회한다.

    Returns:
        {
          "status": str,
          "gene": str,
          "chembl_target_id": str | None,
          "max_phase": int,
          "n_compounds": int,
          "approved_drugs": list[str],
          "source": "ChEMBL REST API",
        }
    """
    # 1단계: 타겟 검색
    url = (f"https://www.ebi.ac.uk/chembl/api/data/target/search.json"
           f"?q={urllib.parse.quote(gene)}&limit=20")
    data = _http_get_json(url)
    if not data or not data.get("targets"):
        return {"status": "UNAVAILABLE" if data is None else "NOT_FOUND", "gene": gene,
                "chembl_target_id": None, "max_phase": 0,
                "n_compounds": 0, "approved_drugs": [],
                "source": "ChEMBL REST API"}

    targets = data["targets"]
    target_id = None
    for target in targets:
        if target.get("organism") != "Homo sapiens" or target.get("target_type") != "SINGLE PROTEIN":
            continue
        synonyms = {str(item.get("component_synonym", "")).upper()
                    for component in target.get("target_components", [])
                    for item in component.get("target_component_synonyms", [])}
        if gene.upper() in synonyms:
            target_id = target["target_chembl_id"]
            break
    if target_id is None:
        return {"status": "NOT_FOUND", "gene": gene, "chembl_target_id": None,
                "max_phase": 0, "n_compounds": 0, "approved_drugs": [],
                "reason": "No exact human single-protein gene match", "source": "ChEMBL REST API"}

    # 2단계: 활성 화합물 조회
    act_url = (f"https://www.ebi.ac.uk/chembl/api/data/activity.json"
               f"?target_chembl_id={target_id}&pchembl_value__gte=6&standard_relation=%3D&assay_type=B&limit=1000&format=json")
    act_data = _http_get_json(act_url)
    compounds: set[str] = set()
    max_phase = 0
    approved: list[str] = []

    if act_data:
        for row in (act_data.get("activities") or []):
            mol_id = row.get("molecule_chembl_id")
            if mol_id:
                compounds.add(mol_id)
            phase = row.get("max_phase") or 0
            try:
                phase = int(phase)
            except (ValueError, TypeError):
                phase = 0
            if phase > max_phase:
                max_phase = phase
            if phase == 4:
                mol_name = row.get("molecule_pref_name") or mol_id
                if mol_name and mol_name not in approved:
                    approved.append(mol_name)

    mechanism_url = (f"https://www.ebi.ac.uk/chembl/api/data/mechanism.json"
                     f"?target_chembl_id={target_id}&limit=1000")
    mechanisms = _http_get_json(mechanism_url)
    for mechanism in (mechanisms or {}).get("mechanisms", []):
        if mechanism.get("direct_interaction") != 1:
            continue
        phase = mechanism.get("max_phase")
        if isinstance(phase, (int, float)):
            max_phase = max(max_phase, phase)
            if phase == 4 and mechanism.get("molecule_chembl_id") not in approved:
                approved.append(mechanism["molecule_chembl_id"])

    return {
        "status": "FOUND",
        "gene": gene,
        "chembl_target_id": target_id,
        "max_phase": max_phase,
        "n_compounds": len(compounds),
        "activity_status": "QUERIED" if act_data is not None else "UNAVAILABLE",
        "approved_drugs": approved[:5],
        "activity_filter": "Human single protein; binding assay, exact relation, pChEMBL >= 6; first 1000 records",
        "clinical_phase_source": mechanism_url,
        "clinical_phase_status": "QUERIED" if mechanisms is not None else "UNAVAILABLE",
        "clinical_phase_note": "Molecule development phase with a recorded direct target mechanism; not disease-specific efficacy",
        "source": "ChEMBL REST API",
    }


# ── Evidence Card 조립 ────────────────────────────────────

def build_evidence_card(
    gene: str,
    source: str,
    claim: str,
    status: EvidenceStatus,
    result: dict,
    scope: str = "",
) -> dict:
    """워크플로우 §4 Evidence Card 형식."""
    import datetime
    import hashlib
    content = json.dumps(result, ensure_ascii=False, sort_keys=True)
    return {
        "id": f"EV-{source.upper()}-{gene}-{hashlib.md5(content.encode()).hexdigest()[:8]}",
        "gene": gene,
        "source": source,
        "claim": claim,
        "status": status,
        "scope": scope,
        "result": result,
        "version": result.get("source", source),
        "timestamp": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "content_hash": hashlib.md5(content.encode()).hexdigest(),
    }


# ── safety veto 판단 ──────────────────────────────────────

def apply_safety_veto(
    gene: str,
    hpa_result: dict,
    uniprot_result: dict,
    modality: str,
    rules_version: str = "omicraft-v5",
) -> dict:
    """명확히 사전 정의된 critical safety rule에만 veto를 적용한다.

    워크플로우 §06:
      - 정상 RNA 고발현 하나만으로 모든 modality 일괄 REJECT 하지 않는다.
      - 개입 방식·노출·근거 수준에 맞는 범위를 명시한다.
      - 불명확한 위험은 HOLD + 필요한 검증으로 남긴다.
      - veto를 Agent가 점수로 상쇄할 수 없다.

    Returns:
        {"has_veto": bool, "veto_scope": str, "veto_reason": str,
         "hold_reasons": list[str], "rules_version": str}
    """
    veto = False
    veto_scope = ""
    veto_reason = ""
    hold_reasons: list[str] = []

    critical_tissues = hpa_result.get("high_expression_tissues", [])
    is_critical = hpa_result.get("critical_tissue_expression", False)
    is_surface = _location_to_target_info(uniprot_result)["surface_accessible"]

    # CELL_TARGETING modality (ADC/CAR-T)
    if modality in ("ADC", "CAR_T", "CELL_TARGETING"):
        if is_critical:
            hold_reasons.append(
                f"정상 critical 조직 ({', '.join(critical_tissues[:3])}) 발현 확인 — "
                "ADC/CAR-T on-target/off-tumor 위험 평가 필요 (HOLD, 자동 REJECT 아님)"
            )
        # 표면 접근성 없으면 N/A
        if not is_surface and uniprot_result.get("status") == "FOUND":
            veto = True
            veto_scope = f"{modality}:{gene}"
            veto_reason = "세포 표면 접근성 없음 — 항원 전달 불가"

    # FUNCTION_MODULATION
    elif modality in ("SMALL_MOLECULE", "INHIBITOR", "FUNCTION_MODULATION"):
        if is_critical and not is_surface:
            hold_reasons.append(
                f"정상 critical 조직 발현 ({', '.join(critical_tissues[:3])}) + "
                "세포 내 표적 — 선택성·노출 평가 필요 (HOLD)"
            )

    return {
        "has_veto": veto,
        "veto_scope": veto_scope,
        "veto_reason": veto_reason,
        "hold_reasons": hold_reasons,
        "rules_version": rules_version,
    }


# ── Shortlist 규칙 ────────────────────────────────────────

def apply_shortlist_rules(
    gene: str,
    bio_confidence: str,
    tier: str,
    has_scope_issue: bool = False,
    has_critical_veto: bool = False,
    dual_contrast: str = "partial_evidence",
    has_external_replication: bool = False,
    prevalence: float = 0.0,
    n_unresolved_evidence: int = 0,
) -> dict:
    """워크플로우 §07 Shortlist 규칙 적용.

    규칙 순서:
      1. scope/검정 품질 불충족 → HOLD
      2. B0 또는 critical veto → REJECT
      3. B1 → B2 → B3 층화
      4. 같은 층: 두 비교 지지 → 독립 cohort → prevalence → 미해결 수

    Returns:
        {"judgment": "ADVANCE"|"HOLD"|"REJECT", "tier": str,
         "bio_confidence": str, "reason": str, "sort_key": tuple}
    """
    if has_scope_issue:
        return {"judgment": "HOLD", "tier": tier, "bio_confidence": bio_confidence,
                "reason": "scope/검정 품질 불충족", "sort_key": (99, 0)}

    if has_critical_veto or bio_confidence == "B0":
        reason = "B0 반증" if bio_confidence == "B0" else "critical safety veto"
        return {"judgment": "REJECT", "tier": tier, "bio_confidence": bio_confidence,
                "reason": reason, "sort_key": (100, 0)}

    # HOLD 조건 (필수 항목 UNKNOWN)
    if tier in ("HOLD", "UNASSESSED"):
        return {"judgment": "HOLD", "tier": tier, "bio_confidence": bio_confidence,
                "reason": "필수 항목 UNKNOWN — 추가 근거 필요", "sort_key": (98, 0)}

    b_rank = {"B1": 0, "B2": 1, "B3": 2}.get(bio_confidence, 3)
    dual_rank  = 0 if dual_contrast == "dual_contrast_supported" else 1
    ext_rank   = 0 if has_external_replication else 1
    prev_rank  = -round(prevalence, 3) if prevalence is not None else 1.0  # 높을수록 우선
    unres_rank = n_unresolved_evidence

    sort_key = (b_rank, dual_rank, ext_rank, prev_rank, unres_rank)

    judgment = "ADVANCE" if bio_confidence in ("B1", "B2") and tier not in ("HOLD", "REJECT", "NOT_APPLICABLE") \
               else "HOLD" if bio_confidence == "B3" \
               else "HOLD"

    return {
        "judgment": judgment,
        "tier": tier,
        "bio_confidence": bio_confidence,
        "reason": f"B={bio_confidence}, Tier={tier}, dual={dual_contrast}",
        "sort_key": sort_key,
    }


# ── 통합 유전자 qualification ─────────────────────────────

def qualify_gene(
    gene: str,
    de_evidence: dict,
    cell_context: dict,
    depmap_result: dict,
    modalities: list[str] | None = None,
    offline: bool = False,
) -> dict:
    """단일 유전자의 전체 qualification을 수행한다.

    Args:
        gene: 유전자명.
        de_evidence: {"dual_contrast": str, "log2FC": float, "padj": float, ...}
        cell_context: {"context_status": str, "top_class": str, ...}
        depmap_result: {"dependency_status": str, "crispr_summary": dict, ...}
        modalities: 평가할 modality 목록. None이면 기본 목록.
        offline: True이면 외부 API 조회 건너뜀 (테스트용).

    Returns:
        {
          "gene": str,
          "uniprot": dict,
          "hpa": dict,
          "chembl": dict,
          "target_info": dict,
          "asset_info": dict,
          "bio_confidence": str,
          "tier_results": list[dict],
          "shortlist": dict,
          "evidence_cards": list[dict],
          "intervention_role": str,
          "safety_veto": dict,
          "missing_evidence": list[str],
        }
    """
    if modalities is None:
        modalities = ["SMALL_MOLECULE", "DE_NOVO_BINDER", "ADC", "DEGRADER"]

    evidence_cards: list[dict] = []
    missing_evidence: list[str] = []

    # ── 1. 외부 API 조회 ─────────────────────────────────
    if offline:
        uniprot_result = {"status": "SKIPPED", "gene": gene,
                          "subcellular_location": [], "is_secreted": False,
                          "has_transmembrane": False, "n_isoforms": 0,
                          "function": ""}
        hpa_result     = {"status": "SKIPPED", "gene": gene,
                          "normal_tissue": [], "high_expression_tissues": [],
                          "critical_tissue_expression": False}
        chembl_result  = {"status": "SKIPPED", "gene": gene,
                          "max_phase": 0, "n_compounds": 0}
    else:
        uniprot_result = query_uniprot(gene)
        hpa_result     = query_hpa_normal(gene, de_evidence.get("gene_id"))
        chembl_result  = query_chembl(gene)

        # Evidence Cards
        if uniprot_result["status"] == "FOUND":
            evidence_cards.append(build_evidence_card(
                gene, "UniProt",
                f"{gene} 단백질 위치: {uniprot_result['subcellular_location']}",
                "SUPPORTIVE", uniprot_result,
            ))
        else:
            missing_evidence.append(f"UniProt_{gene}")

        if hpa_result["status"] == "FOUND":
            evidence_cards.append(build_evidence_card(
                gene, "HPA",
                f"{gene} 정상조직 발현 (high: {hpa_result['high_expression_tissues'][:3]})",
                "SUPPORTIVE" if not hpa_result["critical_tissue_expression"] else "CONTRADICTORY",
                hpa_result,
            ))
        else:
            missing_evidence.append(f"HPA_{gene}")

    gene_id = (de_evidence.get("gene_id") or "").split(".")[0]
    association = query_open_targets(gene_id) if gene_id.startswith("ENSG") and not offline else {"status": "SKIPPED"}
    if association.get("status") == "FOUND":
        evidence_cards.append(build_evidence_card(gene, "OpenTargets", "Breast cancer association; not causal or functional validation",
                                                 "SUPPORTIVE", association))
    elif not offline:
        missing_evidence.append("OpenTargets_" + gene)

    # ── 2. target_info / asset_info 생성 ─────────────────
    target_info = _location_to_target_info(uniprot_result)
    target_info["has_structure"] = bool(uniprot_result.get("pdb_ids"))

    asset_info = {
        "chembl_max_phase": chembl_result.get("max_phase", 0),
        "chembl_n_compounds": chembl_result.get("n_compounds", 0),
        "has_structure": target_info["has_structure"],
        "approved_drugs": chembl_result.get("approved_drugs", []),
    }

    # ── 3. intervention_role 판단 ─────────────────────────
    if target_info["surface_accessible"] and "ADC" in modalities:
        intervention_role = "CELL_TARGETING"
    elif target_info["is_secreted"]:
        intervention_role = "FUNCTION_MODULATION"
    else:
        intervention_role = "FUNCTION_MODULATION"

    # ── 4. Biological confidence (B) ─────────────────────
    dep_status = depmap_result.get("dependency_status", "UNAVAILABLE")
    context_status = cell_context.get("context_status", "INSUFFICIENT")
    dependency = depmap_result.get("crispr_summary") or {}
    if dependency.get("cohort_scope") == "BREAST_LINEAGE":
        dep_status = "UNAVAILABLE"
        missing_evidence.append("TNBC-specific dependency: receptor metadata unavailable")
    if dep_status == "CONCORDANT" and dependency.get("median_effect", 0) is not None and dependency.get("median_effect", 0) >= -0.5:
        dep_status = "UNAVAILABLE"
        missing_evidence.append("Concordant non-dependency is not functional support")

    bio_cell_status = context_status
    if context_status == "SINGLE_SUPPORTED" and cell_context.get("broad_top_class") != "malignant":
        bio_cell_status = "UNRESOLVED"
        missing_evidence.append("Malignant origin unconfirmed; cell-type specificity is descriptive")

    bio_evidence = {
        "dual_contrast": de_evidence.get("dual_contrast", "partial_evidence"),
        "cell_origin": bio_cell_status,
        "depmap_status": dep_status,
        "crispr_effect": dependency.get("median_effect", 0) or 0,
        "has_functional_evidence": de_evidence.get("has_functional_evidence", False),
        "has_safety_veto": False,  # safety veto 적용 후 갱신
        "has_external_replication": de_evidence.get("has_external_replication", False),
        "has_causal_evidence": de_evidence.get("has_causal_evidence", False),
        "ot_genetic_association": association.get("genetic_association_score") if association.get("status") == "FOUND" else None,
        "ot_somatic_mutation": association.get("somatic_mutation_score") if association.get("status") == "FOUND" else None,
    }
    b, b_reasons = assign_bio_confidence(bio_evidence)

    # ── 5. Safety veto ────────────────────────────────────
    safety_veto_results: dict[str, dict] = {}
    overall_veto = False
    for mod in modalities:
        sv = apply_safety_veto(gene, hpa_result, uniprot_result, mod)
        safety_veto_results[mod] = sv
        if sv["has_veto"] and sv["veto_scope"] == "GLOBAL":
            overall_veto = True

    if overall_veto:
        bio_evidence["has_safety_veto"] = True
        b = "B0"
        b_reasons = [f"Critical safety veto ({gene})"]

    # ── 6. modality × B/M/D Tier ─────────────────────────
    from .tier_system import evaluate_target_modality
    tier_results = evaluate_target_modality(
        gene=gene,
        bio_evidence=bio_evidence,
        target_info=target_info,
        asset_info=asset_info,
        modalities=modalities,
    )

    for row in tier_results:
        safety = safety_veto_results[row["modality"]]
        row["safety_veto"] = safety
        if row["modality"] == "ADC" and uniprot_result.get("status") != "FOUND":
            row["Tier"] = "HOLD"
            row["M"] = "UNKNOWN"
        if safety["has_veto"]:
            row["Tier"] = "NOT_APPLICABLE" if "표면" in safety["veto_reason"] else "REJECT"
        elif safety["hold_reasons"] and row["Tier"] not in {"REJECT", "NOT_APPLICABLE"}:
            row["Tier"] = "HOLD"

    # 대표 Tier (가장 좋은 Tier)
    tier_order = {"1A": 0, "1B": 1, "2A": 2, "2B": 3, "3": 4,
                  "HOLD": 5, "NOT_APPLICABLE": 6, "REJECT": 7}
    best_tier = min(tier_results, key=lambda r: tier_order.get(r["Tier"], 99),
                    default={"Tier": "HOLD"})["Tier"]

    # ── 7. Shortlist 판정 ─────────────────────────────────
    shortlist = apply_shortlist_rules(
        gene=gene,
        bio_confidence=b,
        tier=best_tier,
        has_scope_issue=de_evidence.get("has_scope_issue", False),
        has_critical_veto=overall_veto,
        dual_contrast=de_evidence.get("dual_contrast", "partial_evidence"),
        has_external_replication=de_evidence.get("has_external_replication", False),
        prevalence=de_evidence.get("prevalence", 0.0),
        n_unresolved_evidence=len(missing_evidence),
    )

    return {
        "gene": gene,
        "uniprot": uniprot_result,
        "hpa": hpa_result,
        "chembl": chembl_result,
        "open_targets": association,
        "target_info": target_info,
        "asset_info": asset_info,
        "bio_confidence": b,
        "bio_reasons": b_reasons,
        "tier_results": tier_results,
        "best_tier": best_tier,
        "shortlist": shortlist,
        "evidence_cards": evidence_cards,
        "intervention_role": intervention_role,
        "safety_veto": safety_veto_results,
        "missing_evidence": missing_evidence,
    }


# ── 복수 유전자 qualification 파이프라인 ───────────────────

def run_qualification_pipeline(
    genes: list[str],
    de_evidences: dict[str, dict],
    cell_contexts: dict[str, dict],
    depmap_results: dict[str, dict],
    output_dir: str | Path,
    modalities: list[str] | None = None,
    offline: bool = False,
) -> dict:
    """복수 유전자의 qualification을 실행하고 결과를 저장.

    Args:
        genes: 평가할 유전자 목록.
        de_evidences: {gene: de_evidence_dict}
        cell_contexts: {gene: cell_context_dict}
        depmap_results: {gene: depmap_evaluation_dict}
        output_dir: 출력 디렉토리.
        modalities: 평가할 modality 목록.
        offline: 외부 API 건너뜀 여부.

    Returns:
        {
          "status": str,
          "n_genes": int,
          "advance": list[str],
          "hold": list[str],
          "reject": list[str],
          "results": {gene: qualification_dict},
          "all_evidence_cards": list[dict],
          "hypotheses_tsv": str | None,
        }
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    results: dict[str, dict] = {}
    all_cards: list[dict] = []
    advance, hold, reject = [], [], []

    from concurrent.futures import ThreadPoolExecutor

    from .configuration import configuration_scope, get_config
    config = get_config()
    def qualify(gene):
        with configuration_scope(config):
            return qualify_gene(gene, de_evidences.get(gene, {}), cell_contexts.get(gene, {}),
                                depmap_results.get(gene, {}), modalities=modalities, offline=offline)
    with ThreadPoolExecutor(max_workers=4) as pool:
        from contextvars import copy_context
        futures = [pool.submit(copy_context().run, qualify, gene) for gene in genes]
        evaluated = [future.result() for future in futures]
    for gene, qr in zip(genes, evaluated):
        results[gene] = qr
        all_cards.extend(qr.get("evidence_cards", []))

        judgment = qr["shortlist"]["judgment"]
        if judgment == "ADVANCE":
            advance.append(gene)
        elif judgment == "REJECT":
            reject.append(gene)
        else:
            hold.append(gene)

    # hypotheses.tsv 저장
    hyp_path = _write_hypotheses_tsv(genes, results, output_dir / "hypotheses.tsv")

    # evidence_cards.jsonl 저장
    cards_path = output_dir / "evidence_cards.jsonl"
    try:
        with cards_path.open("w", encoding="utf-8") as f:
            for card in all_cards:
                f.write(json.dumps(card, ensure_ascii=False) + "\n")
    except OSError:
        pass

    (output_dir / "qualification_results.json").write_text(
        json.dumps(results, ensure_ascii=False, indent=2), encoding="utf-8")
    result = {
        "status": "COMPLETED",
        "n_genes": len(genes),
        "advance": advance,
        "hold": hold,
        "reject": reject,
        "results": results,
        "all_evidence_cards": all_cards,
        "hypotheses_tsv": hyp_path,
        "evidence_cards_jsonl": str(cards_path),
    }

    try:
        output_dir.joinpath("qualification_status.json").write_text(
            json.dumps({k: v for k, v in result.items()
                        if k not in ("results", "all_evidence_cards")},
                       ensure_ascii=False, indent=2),
            encoding="utf-8",
        )
    except OSError:
        pass

    return result


def _write_hypotheses_tsv(
    genes: list[str],
    results: dict[str, dict],
    path: Path,
) -> str:
    """hypotheses.tsv 작성 (워크플로우 §4 형식)."""
    header = [
        "gene", "cell_context", "intervention_role", "therapeutic_direction",
        "causal_support", "safety_scope", "bio_confidence", "judgment",
        "tier", "missing_evidence",
    ]
    rows = []
    for gene in genes:
        qr = results.get(gene, {})
        rows.append({
            "gene": gene,
            "cell_context": qr.get("target_info", {}).get("location", "UNKNOWN"),
            "intervention_role": qr.get("intervention_role", "UNKNOWN"),
            "therapeutic_direction": "UNRESOLVED",
            "causal_support": "functional_evidence"
                if qr.get("bio_confidence") in ("B1", "B2") else "UNKNOWN",
            "safety_scope": "critical_veto"
                if any(sv.get("has_veto")
                       for sv in (qr.get("safety_veto") or {}).values())
                else "no_veto",
            "bio_confidence": qr.get("bio_confidence", "B3"),
            "judgment": qr.get("shortlist", {}).get("judgment", "HOLD"),
            "tier": qr.get("best_tier", "HOLD"),
            "missing_evidence": "; ".join(qr.get("missing_evidence", [])),
        })

    lines = ["\t".join(header)]
    for r in rows:
        lines.append("\t".join(str(r.get(k, "")) for k in header))
    try:
        path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    except OSError:
        pass
    return str(path)
