"""Display titles for report figures, independent of artifact filenames."""

import csv
import re
from pathlib import Path


FIGURE_TITLES = {
    "before_batch_pca_group": "PCA by Tumor Group — Before Batch Correction",
    "before_batch_pca_tss": "PCA by Tissue Source Site — Before Batch Correction",
    "before_batch_sample_distance": "Sample Distance Heatmap — Before Batch Correction",
    "after_batch_pca_group": "PCA by Tumor Group — After Batch Correction",
    "after_batch_pca_tss": "PCA by Tissue Source Site — After Batch Correction",
    "after_batch_sample_distance": "Sample Distance Heatmap — After Batch Correction",
    "deg_heatmap": "Differential Gene Expression Heatmap",
    "volcano": "Differential Expression Volcano Plot",
    "gobp_nes_top": "GO Biological Process Enrichment",
    "hallmark_nes_top": "Hallmark Pathway Enrichment",
    "cell_specificity": "Cell-type Gene Specificity",
    "dependency_effects": "Gene Dependency Effects",
    "quality_map": "Binder Quality Assessment",
    "docking_scores": "Docking Pose Scores",
    "gnina_pose_scores": "GNINA Docking Pose Scores",
    "docking_af3_rmsd": "Docking and AF3 Pose Comparison",
    "conjugation_sasa": "Conjugation Site Surface Accessibility",
    "degrader_2d": "Degrader Chemical Structure",
    "degrader_3d": "Degrader Molecular Conformation",
}
_IMAGE_SUFFIXES = {".png", ".jpg", ".jpeg", ".svg", ".webp"}


def sequence_variation_notice(path):
    """Explain empty bars from the saved source table, never from image pixels."""
    match = re.fullmatch(r"sequence_variation_(\d+)", path.stem.lower())
    if not match:
        return None
    try:
        with (path.parent / 'sequence_variation.csv').open(newline='') as stream:
            rows = [r for r in csv.DictReader(stream) if int(r['group']) == int(match[1])]
        if not rows:
            return '이 backbone의 서열 변화 분석 데이터가 저장되지 않았습니다.'
        counts = {int(r['sequence_count']) for r in rows}
        if counts == {1}:
            return '비교할 서열 부족 · 이 backbone에는 고유 서열이 1개만 있어 서열 간 변화량을 비교할 수 없습니다.'
        if counts and min(counts) >= 2 and all(float(r['mutation_frequency']) == 0 for r in rows):
            return '모든 위치의 변화 빈도가 0입니다. 저장된 서열 간 차이가 없습니다.'
    except (OSError, ValueError, KeyError, TypeError):
        pass
    return None


def label_figure(artifact):
    """Set a display title while retaining paths and download labels."""
    path = Path(artifact.get("path", ""))
    if path.suffix.lower() not in _IMAGE_SUFFIXES:
        return artifact
    notice = sequence_variation_notice(path)
    if notice:
        artifact['display_notice'] = notice
    else:
        artifact.pop('display_notice', None)
    for key in ("title", "label"):
        title = str(artifact.get(key) or "").strip()
        if title and Path(title).suffix.lower() not in _IMAGE_SUFFIXES and title != path.stem:
            artifact["title"] = title
            return artifact
    stem = path.stem.lower()
    title = FIGURE_TITLES.get(stem)
    if stem.startswith("gobp_network"):
        title = "GO Biological Process Network"
    elif match := re.fullmatch(r"sequence_variation_(\d+)", stem):
        title = f"Binder Sequence Variation — Backbone {int(match[1])}"
    elif re.fullmatch(r"binder_\d+_confidence", stem):
        title = "Binder Structure Confidence"
    elif re.fullmatch(r"binder_\d+_pae", stem):
        title = "Predicted Aligned Error (PAE)"
    if not title:
        words = re.sub(r"[_-]+", " ", path.stem).split()
        acronyms = {"af3", "pca", "pae", "rmsd", "deg", "gsea", "gnina", "adc", "sasa"}
        title = " ".join(word.upper() if word.lower() in acronyms else word.capitalize() for word in words)
    artifact["title"] = title or "Analysis Results"
    return artifact
