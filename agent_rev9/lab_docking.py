"""
Lab Docking Pipeline — ADC / SMALL_MOLECULE / PROTAC 모달리티별 구조 분석
서버의 /lab 페이지에서 버튼 클릭 시 호출된다.

지원 기능:
  1. ADC      — Boltz2 단백질 단독 구조 + 표면 에피토프 하이라이트
  2. SM       — PDB 결정 구조 다운로드 → AutoDock Vina 도킹 → 복합체 3D 뷰
  3. PROTAC   — 표적 단백질 + E3 ligase 구조 + 개념 삼자 복합체 뷰

ESR1 기준 기본값:
  - PDB ID (SM):     3ERT (ESR1 + 타목시펜)
  - PDB ID (PROTAC): 3ERT (표적) + 4CI2 (CRBN E3 ligase)
  - 결합 포켓:       ESR1 LBD (타목시펜 결합 위치 기준)
"""

from __future__ import annotations

import json
import logging
import re
import subprocess
import urllib.request
from pathlib import Path
from typing import Optional

from .docking_backend import resolve_executable

log = logging.getLogger(__name__)

VINA_EXE = resolve_executable("vina")["path"] or "vina"

# ESR1 타목시펜 결합 포켓 좌표 (3ERT 기준)
ESR1_POCKET = {
    "center": (38.0, -15.5, -11.0),
    "size":   (22.0,  22.0,  22.0),
}

# 리간드 SMILES
LIGANDS = {
    name: record['response']['PropertyTable']['Properties'][0]['SMILES']
    for name, record in json.loads((Path(__file__).parent / 'configs/lab_ligands.json').read_text()).items()
}

PROTAC_WARHEAD_SMILES = "COc1ccc2c(c1)N(C)C(=O)c1cc(-c3ccc(NC(=O)CCl)cc3)ccc1-2"  # 개념 warhead

# E3 ligase PDB: 4TZ4 = Cereblon(CRBN)+thalidomide (PROTAC 연구 표준)
DEFAULT_E3_PDB = "4TZ4"


# ── PDB 다운로드 ────────────────────────────────────────────

def download_pdb(pdb_id: str, cache_dir: Path) -> Path:
    """RCSB에서 PDB 파일 다운로드 (캐시 사용)."""
    cache_dir.mkdir(parents=True, exist_ok=True)
    out = cache_dir / f"{pdb_id.lower()}.pdb"
    if out.exists():
        return out
    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
    log.info("Downloading PDB %s ...", pdb_id)
    try:
        urllib.request.urlretrieve(url, str(out))
    except Exception as e:
        raise RuntimeError(f"PDB {pdb_id} 다운로드 실패: {e}") from e
    return out


def extract_ligand_center(pdb_path: Path, hetatm_code: str) -> tuple[float, float, float]:
    """PDB에서 리간드 HETATM 좌표의 중심 계산."""
    xs, ys, zs = [], [], []
    with open(pdb_path, encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith("HETATM") and hetatm_code in line[17:20]:
                try:
                    xs.append(float(line[30:38]))
                    ys.append(float(line[38:46]))
                    zs.append(float(line[46:54]))
                except ValueError:
                    pass
    if not xs:
        return ESR1_POCKET["center"]
    return (sum(xs)/len(xs), sum(ys)/len(ys), sum(zs)/len(zs))


def strip_hetatm(pdb_path: Path, output_path: Path) -> Path:
    """PDB에서 HETATM(리간드/물) 제거 — receptor only."""
    lines = []
    with open(pdb_path, encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not (line.startswith("HETATM") or line.startswith("CONECT")):
                lines.append(line)
    output_path.write_text("".join(lines), encoding="utf-8")
    return output_path


# ── PDBQT 변환 ─────────────────────────────────────────────

_AD4_MAP = {
    "C": "C", "N": "NA", "O": "OA", "S": "SA",
    "H": "HD", "P": "P", "F": "F", "CL": "Cl", "BR": "Br", "I": "I",
}


def _pdb_to_pdbqt(pdb_path: Path, output_path: Path) -> Path:
    """간이 PDB → PDBQT 변환 (charge=0, AutoDock atom type 매핑).

    PDBQT 컬럼 규격:
      1-66  : 표준 PDB (record~B-factor)
      67-76 : partial charge (10.3f)
      77    : space
      78-79 : AutoDock atom type (2 chars, 대소문자 구분)
    """
    lines = []
    with open(pdb_path, encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                element = line[76:78].strip().upper() if len(line) > 76 else ""
                if not element:
                    element = line[12:16].strip().lstrip("0123456789")[:1].upper()
                ad4_type = _AD4_MAP.get(element, "C")
                # col 1-66: 표준 PDB 필드 (없으면 공백 패딩)
                pdb_part = line[:66].ljust(66)
                lines.append(f"{pdb_part}{0.0:10.3f} {ad4_type:<2s}\n")
            elif line.startswith(("TER", "REMARK", "MODEL", "ENDMDL", "END")):
                lines.append(line)
    output_path.write_text("".join(lines), encoding="utf-8")
    return output_path


def _smiles_to_pdbqt(smiles: str, output_path: Path, name: str = "LIG") -> Path:
    """RDKit으로 SMILES → 3D SDF → 간이 PDBQT 변환."""
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"SMILES 파싱 실패: {smiles[:40]}")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    AllChem.MMFFOptimizeMolecule(mol)

    conf = mol.GetConformer()
    lines = ["ROOT\n"]
    atom_idx = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 1:
            continue  # 수소 제외 (Vina heavy-atom only)
        atom_idx += 1
        pos = conf.GetAtomPosition(atom.GetIdx())
        symbol = atom.GetSymbol().upper()
        ad4_type = _AD4_MAP.get(symbol, "C")
        lines.append(
            f"HETATM{atom_idx:5d}  {symbol:<3s} {name:>3s}     1    "
            f"{pos.x:8.3f}{pos.y:8.3f}{pos.z:8.3f}"
            f"  1.00  0.00{0.0:10.3f} {ad4_type:<2s}\n"
        )
    lines.append("ENDROOT\n")
    lines.append("TORSDOF 0\n")
    output_path.write_text("".join(lines), encoding="utf-8")
    return output_path


# ── Vina 도킹 ──────────────────────────────────────────────

def run_vina(
    receptor_pdbqt: Path,
    ligand_pdbqt: Path,
    center: tuple[float, float, float],
    box_size: tuple[float, float, float],
    output_dir: Path,
    exhaustiveness: int = 8,
    num_modes: int = 5,
) -> dict:
    """AutoDock Vina 실행 → 결과 파싱."""
    output_dir.mkdir(parents=True, exist_ok=True)
    out_pdbqt = output_dir / "docked.pdbqt"
    _log_file  = output_dir / "vina.log"

    cmd = [
        VINA_EXE,
        "--receptor",      str(receptor_pdbqt),
        "--ligand",        str(ligand_pdbqt),
        "--center_x",      str(center[0]),
        "--center_y",      str(center[1]),
        "--center_z",      str(center[2]),
        "--size_x",        str(box_size[0]),
        "--size_y",        str(box_size[1]),
        "--size_z",        str(box_size[2]),
        "--out",           str(out_pdbqt),
        "--exhaustiveness", str(exhaustiveness),
        "--num_modes",      str(num_modes),
        "--cpu",            "4",
    ]

    log.info("Vina 실행: %s", " ".join(cmd))
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=300)

    if proc.returncode != 0:
        raise RuntimeError(f"Vina 실패:\n{proc.stderr}")

    # stdout에서 결과 파싱 (--log 없이 stdout으로 출력됨)
    poses = []
    log_text = proc.stdout + proc.stderr
    for line in log_text.splitlines():
        m = re.match(r"\s*(\d+)\s+(-?\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)", line)
        if m:
            poses.append({
                "mode":      int(m.group(1)),
                "affinity":  float(m.group(2)),  # kcal/mol
                "rmsd_lb":   float(m.group(3)),
                "rmsd_ub":   float(m.group(4)),
            })

    best = poses[0] if poses else {}
    return {
        "status":       "DONE",
        "poses":        poses,
        "best_affinity": best.get("affinity"),
        "out_pdbqt":    str(out_pdbqt) if out_pdbqt.exists() else None,
        "log":          log_text[:2000],
    }


def pdbqt_to_pdb_str(pdbqt_path: Path, model_idx: int = 0) -> str:
    """PDBQT 파일에서 첫 번째 모델의 PDB 형식 문자열 반환."""
    lines = []
    in_model = False
    model_count = -1
    with open(pdbqt_path, encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith("MODEL"):
                model_count += 1
                in_model = (model_count == model_idx)
                continue
            if line.startswith("ENDMDL"):
                if in_model:
                    break
                in_model = False
                continue
            if in_model or (model_count == -1):
                if line.startswith(("ATOM", "HETATM")):
                    lines.append(line[:80])
    return "".join(lines)


# ── 3Dmol.js HTML 생성 ─────────────────────────────────────

def make_3dmol_html(
    receptor_pdb_str: str,
    ligand_pdb_str: Optional[str] = None,
    ligand2_pdb_str: Optional[str] = None,
    title: str = "구조 뷰어",
    receptor_color: str = "spectrum",
    ligand_color: str = "greenCarbon",
    ligand2_color: str = "yellowCarbon",
    receptor_style: str = "cartoon",   # cartoon | surface | stick
) -> str:
    """3Dmol.js 기반 단백질(+리간드) HTML 뷰어 생성."""

    lig_js = ""
    if ligand_pdb_str:
        lig_escaped = ligand_pdb_str.replace("`", "\\`")
        lig_js = f"""
  var ligData = `{lig_escaped}`;
  viewer.addModel(ligData, 'pdb');
  viewer.setStyle({{model: 1}}, {{stick: {{colorscheme: '{ligand_color}', radius: 0.15}}}});
"""

    lig2_js = ""
    if ligand2_pdb_str:
        lig2_escaped = ligand2_pdb_str.replace("`", "\\`")
        lig2_js = f"""
  var lig2Data = `{lig2_escaped}`;
  viewer.addModel(lig2Data, 'pdb');
  viewer.setStyle({{model: 2}}, {{stick: {{colorscheme: '{ligand2_color}', radius: 0.15}}}});
"""

    rec_escaped = receptor_pdb_str.replace("`", "\\`")
    style_js = (
        f"{{cartoon: {{color: '{receptor_color}'}}}}"
        if receptor_style == "cartoon"
        else "{surface: {opacity: 0.7, color: 'white'}}"
    )

    return f"""<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>{title}</title>
<script src="https://3Dmol.org/build/3Dmol-min.js"></script>
<style>
  body{{margin:0;background:#0d1117}}
  #viewer{{width:100vw;height:100vh}}
  #info{{position:fixed;bottom:12px;left:12px;background:rgba(0,0,0,.7);
         color:#e0e0e0;padding:8px 14px;border-radius:6px;font:12px/1.6 monospace}}
</style>
</head>
<body>
<div id="viewer"></div>
<div id="info">{title}</div>
<script>
var viewer = $3Dmol.createViewer(document.getElementById('viewer'), {{backgroundColor:'#0d1117'}});
var recData = `{rec_escaped}`;
viewer.addModel(recData, 'pdb');
viewer.setStyle({{model: 0}}, {style_js});
{lig_js}
{lig2_js}
viewer.zoomTo();
viewer.render();
</script>
</body>
</html>"""


# ── 모달리티별 메인 함수 ─────────────────────────────────────

def run_adc_view(gene: str, structure_dir: Path, output_dir: Path) -> dict:
    """ADC: 단백질 표면 구조 + 에피토프 하이라이트."""
    gene_l = gene.lower()
    cif_path = structure_dir / gene_l / "output" / f"boltz_results_{gene_l}" / "predictions" / gene_l / f"{gene_l}_model_0.cif"
    pdb_path = structure_dir / gene_l / f"{gene_l}_boltz2.pdb"

    # CIF가 있으면 PDB로 변환 시도 (간이)
    if not pdb_path.exists() and cif_path.exists():
        try:
            from Bio.PDB import PDBIO, MMCIFParser
            parser = MMCIFParser(QUIET=True)
            structure = parser.get_structure(gene_l, str(cif_path))
            io = PDBIO()
            io.set_structure(structure)
            io.save(str(pdb_path))
        except Exception as e:
            log.warning("CIF→PDB 변환 실패: %s", e)

    # Boltz2 HTML이 있으면 그걸 사용
    html_path = structure_dir / gene_l / f"{gene_l}_boltz2_3d.html"
    if html_path.exists():
        return {
            "status": "DONE",
            "mode": "adc",
            "viewer_url": f"/structure/{gene.upper()}/3d",
            "message": f"{gene.upper()} 표면 구조 (Boltz2 예측) — ADC 에피토프 접근성 확인용",
        }

    return {"status": "NO_STRUCTURE", "mode": "adc",
            "message": "Boltz2 구조 파일 없음"}


def run_sm_docking(
    gene: str,
    pdb_id: str,
    ligand_name: str,
    ligand_smiles: str,
    cache_dir: Path,
    output_dir: Path,
) -> dict:
    """SMALL_MOLECULE: PDB 다운로드 → Vina 도킹 → 복합체 뷰어 HTML."""
    output_dir.mkdir(parents=True, exist_ok=True)

    # 1. PDB 다운로드
    pdb_path = download_pdb(pdb_id, cache_dir)

    # 2. 결합 포켓 좌표 자동 추출
    hetatm_code_map = {
        "3ERT": "TAM",   # tamoxifen
        "2Q70": "ICI",   # fulvestrant (ICI 182,780)
        "3UUC": "Z78",
    }
    hetatm_code = hetatm_code_map.get(pdb_id.upper(), "LIG")
    center = extract_ligand_center(pdb_path, hetatm_code)
    if center == ESR1_POCKET["center"]:
        # fallback to known coordinates
        center = ESR1_POCKET["center"]
    box_size = (22.0, 22.0, 22.0)

    # 3. Receptor 준비
    rec_pdb  = output_dir / "receptor.pdb"
    rec_pdbqt = output_dir / "receptor.pdbqt"
    strip_hetatm(pdb_path, rec_pdb)
    _pdb_to_pdbqt(rec_pdb, rec_pdbqt)

    # 4. Ligand 준비
    lig_pdbqt = output_dir / "ligand.pdbqt"
    _smiles_to_pdbqt(ligand_smiles, lig_pdbqt, name="LIG")

    # 5. Vina 실행
    dock_result = run_vina(
        receptor_pdbqt=rec_pdbqt,
        ligand_pdbqt=lig_pdbqt,
        center=center,
        box_size=box_size,
        output_dir=output_dir / "vina_out",
    )

    # 6. 시각화 HTML 생성
    html_path = output_dir / "complex_viewer.html"
    receptor_pdb_str = rec_pdb.read_text(encoding="utf-8", errors="ignore")
    ligand_pdb_str = ""
    if dock_result.get("out_pdbqt"):
        ligand_pdb_str = pdbqt_to_pdb_str(Path(dock_result["out_pdbqt"]), model_idx=0)

    html = make_3dmol_html(
        receptor_pdb_str=receptor_pdb_str,
        ligand_pdb_str=ligand_pdb_str or None,
        title=f"{gene.upper()} + {ligand_name} (Vina 도킹, ΔG={dock_result.get('best_affinity', '?')} kcal/mol)",
        receptor_style="cartoon",
    )
    html_path.write_text(html, encoding="utf-8")

    return {
        "status": dock_result["status"],
        "mode": "small_molecule",
        "gene": gene,
        "pdb_id": pdb_id,
        "ligand": ligand_name,
        "best_affinity_kcal_mol": dock_result.get("best_affinity"),
        "poses": dock_result.get("poses", []),
        "viewer_html_path": str(html_path),
        "center": center,
        "message": (
            f"Vina 도킹 완료 — {ligand_name} → {gene.upper()} "
            f"최적 결합 에너지: {dock_result.get('best_affinity', '?')} kcal/mol"
        ),
    }


def run_protac_view(
    gene: str,
    target_pdb_id: str,
    e3_pdb_id: str,
    cache_dir: Path,
    output_dir: Path,
) -> dict:
    """PROTAC: 표적 + E3 ligase 이중 구조 + 삼자 복합체 개념 뷰어."""
    output_dir.mkdir(parents=True, exist_ok=True)

    # 두 구조 다운로드
    target_pdb = download_pdb(target_pdb_id, cache_dir)
    e3_pdb     = download_pdb(e3_pdb_id, cache_dir)

    # 수용체 정리
    target_rec = output_dir / "target_rec.pdb"
    e3_rec     = output_dir / "e3_rec.pdb"
    strip_hetatm(target_pdb, target_rec)
    strip_hetatm(e3_pdb, e3_rec)

    # E3 구조 평행이동 (겹치지 않도록 x+60)
    e3_lines = []
    with open(e3_rec, encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                try:
                    x = float(line[30:38]) + 60.0
                    line = line[:30] + f"{x:8.3f}" + line[38:]
                except Exception:
                    pass
            e3_lines.append(line)
    e3_shifted = output_dir / "e3_shifted.pdb"
    e3_shifted.write_text("".join(e3_lines), encoding="utf-8")

    # 시각화 HTML (표적=파랑 cartoon, E3=주황 cartoon, 개념 linker=stick)
    target_str = target_rec.read_text(encoding="utf-8", errors="ignore")
    e3_str     = e3_shifted.read_text(encoding="utf-8", errors="ignore")

    html = make_3dmol_html(
        receptor_pdb_str=target_str,
        ligand_pdb_str=e3_str,
        title=f"PROTAC 삼자 복합체 개념 — {gene.upper()}(표적) + {e3_pdb_id}(E3 ligase)",
        receptor_color="blue",
        ligand_color="orangeCarbon",
        receptor_style="cartoon",
    )
    html_path = output_dir / "protac_viewer.html"
    html_path.write_text(html, encoding="utf-8")

    return {
        "status": "DONE",
        "mode": "protac",
        "gene": gene,
        "target_pdb": target_pdb_id,
        "e3_pdb": e3_pdb_id,
        "viewer_html_path": str(html_path),
        "message": (
            f"PROTAC 삼자 복합체 개념 뷰 — "
            f"{gene.upper()}({target_pdb_id}) + E3 ligase({e3_pdb_id})"
        ),
    }
