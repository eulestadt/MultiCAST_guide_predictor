"""Map MultiCAST PAM/guides onto insertion amino acids, including the ~49 bp offset."""

from __future__ import annotations

import re
from dataclasses import dataclass

from Bio.Seq import Seq

from .config import CastConfig
from .gene import GeneRecord
from .structure import DisruptionWindow, StructureAnalysis

PAM_PATTERN = re.compile(r"(?=([ATCG]{3}C[ATCG].{32}))")


@dataclass
class GuideCandidate:
    pam_region: str
    guide_sequence: str
    strand: str  # coding | template
    pam_n_cds: int
    guide_3prime_cds: int
    insertion_nt: int
    insertion_aa: int
    oligo: str


def _build_oligo(guide: str, cfg: CastConfig) -> str:
    return f"{cfg.oligo_left}{guide}{cfg.oligo_right}"


def scan_guides(rec: GeneRecord, cfg: CastConfig) -> list[GuideCandidate]:
    """
    Type I-F3 Tn6677 binds a 32-nt protospacer next to a 5'-CN-3' PAM and
    inserts ~49 bp *downstream of the protospacer 3' end*.

    Coding-strand target: insertion moves toward the C-terminus.
    Template-strand target: insertion moves toward the N-terminus on the CDS.
    """
    seq = rec.cds.upper()
    L = len(seq)
    offset = cfg.insertion_offset_bp
    out: list[GuideCandidate] = []

    for match in PAM_PATTERN.finditer(seq):
        pos = match.start()
        full = match.group(1)
        pam = full[:5]
        guide = full[5:37]
        guide_3 = pos + 36
        insertion = guide_3 + offset
        if insertion < 0 or insertion >= L:
            continue
        out.append(
            GuideCandidate(
                pam_region=pam,
                guide_sequence=guide,
                strand="coding",
                pam_n_cds=pos + 4,
                guide_3prime_cds=guide_3,
                insertion_nt=insertion,
                insertion_aa=insertion // 3,
                oligo=_build_oligo(guide, cfg),
            )
        )

    rc = str(Seq(seq).reverse_complement())
    for match in PAM_PATTERN.finditer(rc):
        pos = match.start()
        full = match.group(1)
        pam = full[:5]
        guide = full[5:37]
        guide_3_rc = pos + 36
        guide_3_cds = L - 1 - guide_3_rc
        insertion = guide_3_cds - offset
        if insertion < 0 or insertion >= L:
            continue
        out.append(
            GuideCandidate(
                pam_region=pam,
                guide_sequence=guide,
                strand="template",
                pam_n_cds=L - pos - 1,
                guide_3prime_cds=guide_3_cds,
                insertion_nt=insertion,
                insertion_aa=insertion // 3,
                oligo=_build_oligo(guide, cfg),
            )
        )
    return out


def jitter_core_fraction(
    candidate: GuideCandidate,
    window: DisruptionWindow,
    cfg: CastConfig,
    protein_len: int,
) -> float:
    hits = 0
    total = 0
    for off in range(cfg.jitter_min_bp, cfg.jitter_max_bp + 1):
        if candidate.strand == "coding":
            nt = candidate.guide_3prime_cds + off
        else:
            nt = candidate.guide_3prime_cds - off
        if nt < 0:
            continue
        aa = nt // 3
        if aa < 0 or aa >= protein_len:
            continue
        total += 1
        if aa in window.core_residues:
            hits += 1
    return hits / total if total else 0.0


def dominant_negative_risk(insertion_aa: int, analysis: StructureAnalysis) -> bool:
    """True if a complete folded domain would remain entirely N-terminal of the cut."""
    for domain in analysis.domains:
        if domain.end <= insertion_aa and domain.mean_plddt >= 70 and domain.length >= 40:
            return True
    return False


def disruption_score(
    candidate: GuideCandidate,
    analysis: StructureAnalysis,
    jitter_frac: float,
) -> tuple[float, str]:
    window = analysis.window
    aa = candidate.insertion_aa
    n = analysis.residue_count
    reasons: list[str] = []
    score = 0.0

    if aa in window.core_residues:
        score += 0.55
        reasons.append("buried core")
    elif aa in window.structured_residues:
        score += 0.30
        reasons.append("structured but not buried")
    elif window.cds_lo_aa <= aa < window.cds_hi_aa:
        score += 0.10
        reasons.append("first-half CDS")
    else:
        reasons.append("outside preferred window")

    score += 0.25 * jitter_frac
    if jitter_frac >= 0.8:
        reasons.append("jitter-robust")

    if window.domain is not None and window.domain.start <= aa < window.domain.end:
        if window.domain.index == 0:
            score += 0.15
            reasons.append("N-terminal domain")
        else:
            score += 0.05
            reasons.append(f"domain {window.domain.index}")

    if dominant_negative_risk(aa, analysis):
        score -= 0.20
        reasons.append("DN-risk: intact N-terminal domain")

    if aa > int(0.7 * n):
        score -= 0.15
        reasons.append("late insertion")

    return float(max(0.0, min(1.0, score))), "; ".join(reasons)
