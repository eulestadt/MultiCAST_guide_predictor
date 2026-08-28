"""Parse AlphaFold pLDDT / PAE / burial and choose a disruption window."""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser

from .alphafold import AlphaFoldEntry
from .config import StructureConfig
from .gene import GeneRecord


@dataclass
class Domain:
    index: int
    start: int  # 0-based inclusive
    end: int  # 0-based exclusive
    mean_plddt: float
    mean_pae: float

    @property
    def length(self) -> int:
        return self.end - self.start


@dataclass
class DisruptionWindow:
    domain: Domain | None
    core_residues: set[int]
    structured_residues: set[int]
    cds_lo_aa: int
    cds_hi_aa: int
    target_aa: int | None
    method: str
    notes: list[str] = field(default_factory=list)


@dataclass
class StructureAnalysis:
    plddt: np.ndarray
    rsa_proxy: np.ndarray
    pae: np.ndarray | None
    domains: list[Domain]
    window: DisruptionWindow
    residue_count: int
    mean_plddt: float


def parse_pdb_ca(pdb_path: Path) -> tuple[np.ndarray, np.ndarray]:
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("af", str(pdb_path))
    coords: list[list[float]] = []
    plddt: list[float] = []
    for atom in structure.get_atoms():
        if atom.get_name() != "CA":
            continue
        coords.append(atom.get_coord().tolist())
        plddt.append(float(atom.get_bfactor()))
    if not coords:
        raise ValueError(f"No CA atoms in {pdb_path}")
    return np.asarray(coords, dtype=float), np.asarray(plddt, dtype=float)


def load_pae(path: Path) -> np.ndarray:
    data = json.loads(path.read_text(encoding="utf-8"))
    if isinstance(data, list) and data and isinstance(data[0], dict):
        data = data[0]
    if isinstance(data, dict):
        if "predicted_aligned_error" in data:
            arr = data["predicted_aligned_error"]
        elif "pae" in data:
            arr = data["pae"]
        else:
            raise KeyError(f"No PAE matrix in {path}")
    else:
        arr = data
    return np.asarray(arr, dtype=float)


def rsa_proxy(coords: np.ndarray, cutoff: float) -> np.ndarray:
    """Invert CA neighbor counts into a 0–1 exposure proxy (1 = exposed)."""
    n = len(coords)
    if n == 0:
        return np.zeros(0)
    diff = coords[:, None, :] - coords[None, :, :]
    dist = np.sqrt((diff ** 2).sum(axis=2))
    counts = (dist < cutoff).sum(axis=1) - 1
    return np.clip(1.0 - counts / 24.0, 0.0, 1.0)


def _prefix2d(pae: np.ndarray) -> np.ndarray:
    s = np.pad(pae, ((1, 0), (1, 0)))
    return np.cumsum(np.cumsum(s, axis=0), axis=1)


def _block_mean(s: np.ndarray, a0: int, a1: int, b0: int, b1: int) -> float:
    area = (a1 - a0) * (b1 - b0)
    if area <= 0:
        return 0.0
    total = s[a1, b1] - s[a0, b1] - s[a1, b0] + s[a0, b0]
    return float(total / area)


def split_domains(
    pae: np.ndarray,
    plddt: np.ndarray,
    start: int,
    end: int,
    cfg: StructureConfig,
    prefix: np.ndarray | None = None,
) -> list[tuple[int, int]]:
    n = end - start
    if n < cfg.min_domain_aa * 2:
        return [(start, end)]
    if prefix is None:
        prefix = _prefix2d(pae)
    best_score = -1e9
    best_s = None
    step = 4 if n > 80 else 1
    for s in range(start + cfg.min_domain_aa, end - cfg.min_domain_aa + 1, step):
        intra_a = _block_mean(prefix, start, s, start, s)
        intra_b = _block_mean(prefix, s, end, s, end)
        inter = 0.5 * (
            _block_mean(prefix, start, s, s, end) + _block_mean(prefix, s, end, start, s)
        )
        linker = float(plddt[max(start, s - 4) : min(end, s + 4)].mean())
        score = inter - 0.5 * (intra_a + intra_b)
        if linker < 60:
            score += 2.0
        if score > best_score:
            best_score = score
            best_s = s
    if best_s is not None and best_score > cfg.pae_split_delta:
        return split_domains(pae, plddt, start, best_s, cfg, prefix) + split_domains(
            pae, plddt, best_s, end, cfg, prefix
        )
    return [(start, end)]


def domains_from_pae(pae: np.ndarray | None, plddt: np.ndarray, cfg: StructureConfig) -> list[Domain]:
    n = len(plddt)
    if pae is None or pae.shape[0] != n:
        # pLDDT-only: merge contiguous high-confidence stretches
        domains: list[Domain] = []
        i = 0
        while i < n:
            if plddt[i] < cfg.plddt_min:
                i += 1
                continue
            j = i + 1
            while j < n and plddt[j] >= cfg.plddt_min - 10:
                j += 1
            if j - i >= cfg.min_domain_aa:
                domains.append(
                    Domain(
                        index=len(domains),
                        start=i,
                        end=j,
                        mean_plddt=float(plddt[i:j].mean()),
                        mean_pae=float("nan"),
                    )
                )
            i = j
        if not domains:
            domains.append(
                Domain(0, 0, n, float(plddt.mean()), float("nan"))
            )
        return domains

    spans = split_domains(pae, plddt, 0, n, cfg)
    out: list[Domain] = []
    for a, b in spans:
        if b - a < max(15, cfg.min_domain_aa // 2):
            continue
        out.append(
            Domain(
                index=len(out),
                start=a,
                end=b,
                mean_plddt=float(plddt[a:b].mean()),
                mean_pae=float(np.mean(np.diag(pae[a:b, a:b]))),
            )
        )
    if not out:
        out.append(Domain(0, 0, n, float(plddt.mean()), float(np.mean(np.diag(pae)))))
    return out


def pick_window(
    rec: GeneRecord,
    plddt: np.ndarray,
    rsa: np.ndarray,
    domains: list[Domain],
    cfg: StructureConfig,
    cds_fraction: tuple[float, float],
) -> DisruptionWindow:
    n = len(plddt)
    lo_frac, hi_frac = cds_fraction
    cds_lo = int(n * lo_frac)
    cds_hi = max(cds_lo + 1, int(n * hi_frac))
    notes: list[str] = []

    structured = {i for i, p in enumerate(plddt) if p >= cfg.plddt_min}
    core = {
        i
        for i in structured
        if rsa[i] <= cfg.rsa_max
    }

    chosen: Domain | None = None
    if cfg.prefer_n_terminal_domain:
        nterm_candidates = [
            d
            for d in domains
            if d.mean_plddt >= cfg.plddt_min - 5 and d.start < cds_hi
        ]
        if nterm_candidates:
            # first well-folded domain, but do not leave it intact
            chosen = nterm_candidates[0]
            notes.append(
                f"Targeting N-terminal domain {chosen.index} "
                f"(residues {chosen.start + 1}-{chosen.end}) so the first fold cannot survive as a fragment."
            )
    if chosen is None and domains:
        # largest domain overlapping the first half
        overlapping = [d for d in domains if d.start < cds_hi and d.end > cds_lo]
        chosen = max(overlapping or domains, key=lambda d: d.length)
        notes.append(
            f"Targeting domain {chosen.index} (residues {chosen.start + 1}-{chosen.end})."
        )

    core_in_domain: set[int] = set()
    if chosen is not None:
        lo = chosen.start + cfg.edge_buffer_aa
        hi = chosen.end - cfg.edge_buffer_aa
        core_in_domain = {i for i in core if lo <= i < hi}
        # If this domain spans past the first-half CDS backstop, keep cuts
        # in the N-terminal half so a C-terminal fragment is not the plan.
        if chosen.end > cds_hi:
            clipped = {i for i in core_in_domain if cds_lo <= i < cds_hi}
            if len(clipped) >= 15:
                core_in_domain = clipped
                notes.append(
                    f"Domain spans past CDS fraction {hi_frac:.0%}; restricting cuts to residues {cds_lo + 1}-{cds_hi}."
                )
        if cfg.skip_disordered:
            structured = {i for i in structured if chosen.start <= i < chosen.end}

    if len(core_in_domain) < 15:
        # Relative burial fallback: lowest-RSA structured residues in the window.
        pool = [
            i
            for i in range(n)
            if i in structured and cds_lo <= i < (chosen.end if chosen is not None else cds_hi)
        ]
        if pool:
            rsa_vals = np.array([rsa[i] for i in pool])
            cutoff = float(np.quantile(rsa_vals, 0.40))
            core_in_domain = {i for i in pool if rsa[i] <= max(cutoff, cfg.rsa_max)}
            notes.append(
                f"Expanded core with relative burial (RSA ≤ {cutoff:.2f} within the window)."
            )

    if not core_in_domain:
        core_in_domain = {i for i in core if cds_lo <= i < cds_hi}
        notes.append("No buried core inside the preferred domain; falling back to first-half buried residues.")
    if not core_in_domain:
        core_in_domain = set(range(cds_lo, cds_hi))
        notes.append("No AlphaFold core available; using CDS-fraction window only (MultiCAST first-half rule).")

    target = int(np.median(sorted(core_in_domain))) if core_in_domain else None
    return DisruptionWindow(
        domain=chosen,
        core_residues=core_in_domain,
        structured_residues=structured,
        cds_lo_aa=cds_lo,
        cds_hi_aa=cds_hi,
        target_aa=target,
        method="alphafold_core" if chosen is not None else "cds_fraction",
        notes=notes,
    )


def analyze_structure(
    rec: GeneRecord,
    entry: AlphaFoldEntry | None,
    cfg: StructureConfig,
    cds_fraction: tuple[float, float],
) -> StructureAnalysis:
    notes_prefix: list[str] = []
    if entry is None:
        n = rec.length_aa
        plddt = np.full(n, 50.0)
        rsa = np.full(n, 0.5)
        domains = [Domain(0, 0, n, 50.0, float("nan"))]
        window = pick_window(rec, plddt, rsa, domains, cfg, cds_fraction)
        window.method = "cds_fraction_fallback"
        window.notes.insert(0, "No AlphaFold DB entry; using first-half CDS targeting from the MultiCAST paper.")
        return StructureAnalysis(plddt, rsa, None, domains, window, n, 50.0)

    coords, plddt = parse_pdb_ca(entry.pdb_path)
    # Align length to protein if AF is shorter/longer
    n = min(len(plddt), rec.length_aa or len(plddt))
    plddt = plddt[:n]
    coords = coords[:n]
    rsa = rsa_proxy(coords, cfg.burial_cutoff_angstrom)
    pae = None
    if entry.pae_path:
        try:
            pae_full = load_pae(entry.pae_path)
            pae = pae_full[:n, :n]
        except Exception as exc:
            notes_prefix.append(f"PAE parse failed ({exc}); using pLDDT domains.")
    domains = domains_from_pae(pae, plddt, cfg)
    window = pick_window(rec, plddt, rsa, domains, cfg, cds_fraction)
    window.notes = notes_prefix + window.notes
    return StructureAnalysis(
        plddt=plddt,
        rsa_proxy=rsa,
        pae=pae,
        domains=domains,
        window=window,
        residue_count=n,
        mean_plddt=float(plddt.mean()),
    )
