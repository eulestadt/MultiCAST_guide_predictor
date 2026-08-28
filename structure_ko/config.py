"""Load and merge the full pipeline configuration."""

from __future__ import annotations

from copy import deepcopy
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

import yaml

REPO_ROOT = Path(__file__).resolve().parents[1]

# Paper: 5'-TACTACTGCAAAGTAGCTGATAAC-[32 nt guide]-CTTTACTGCTGAATAAGTAGATAACTAC-3'
DEFAULT_OLIGO_LEFT = "TACTACTGCAAAGTAGCTGATAAC"
DEFAULT_OLIGO_RIGHT = "CTTTACTGCTGAATAAGTAGATAACTAC"


@dataclass
class OrganismConfig:
    name: str = "Escherichia coli K-12 MG1655"
    preset: str | None = "ecoli_k12"
    taxon_id: int | None = 83333
    assembly: str | None = "GCF_000005845.2"
    uniprot_proteome: str | None = "UP000000625"
    genome: str | None = "example/GCF_000005845.2_ASM584v2_genomic.fna"
    gff3: str | None = "example/GCF_000005845.2_ASM584v2_genomic.gff"


@dataclass
class CastConfig:
    insertion_offset_bp: int = 49
    jitter_min_bp: int = 44
    jitter_max_bp: int = 55
    guide_length: int = 32
    pam: str = "CN"
    cds_fraction: tuple[float, float] = (0.05, 0.50)
    oligo_left: str = DEFAULT_OLIGO_LEFT
    oligo_right: str = DEFAULT_OLIGO_RIGHT


@dataclass
class StructureConfig:
    plddt_min: float = 70.0
    rsa_max: float = 0.30
    pae_split_delta: float = 5.0
    min_domain_aa: int = 40
    edge_buffer_aa: int = 8
    prefer_n_terminal_domain: bool = True
    skip_disordered: bool = True
    burial_cutoff_angstrom: float = 10.0


@dataclass
class ScoringConfig:
    model: str = "model/model.joblib"
    threshold: float = 0.5
    top_n_per_gene: int = 8
    activity_weight: float = 0.55
    disruption_weight: float = 0.45


@dataclass
class OutputConfig:
    dir: str = "results/structure_ko"
    cache_dir: str = "cache"


@dataclass
class PipelineConfig:
    organism: OrganismConfig = field(default_factory=OrganismConfig)
    genes: list[str] = field(default_factory=list)
    genes_file: str | None = None
    cast: CastConfig = field(default_factory=CastConfig)
    structure: StructureConfig = field(default_factory=StructureConfig)
    scoring: ScoringConfig = field(default_factory=ScoringConfig)
    output: OutputConfig = field(default_factory=OutputConfig)

    def resolved_path(self, maybe_rel: str | None) -> Path | None:
        if not maybe_rel:
            return None
        p = Path(maybe_rel)
        if not p.is_absolute():
            p = REPO_ROOT / p
        return p


ORGANISM_PRESETS: dict[str, dict[str, Any]] = {
    "ecoli_k12": {
        "name": "Escherichia coli K-12 MG1655",
        "taxon_id": 83333,
        "assembly": "GCF_000005845.2",
        "uniprot_proteome": "UP000000625",
        "genome": "example/GCF_000005845.2_ASM584v2_genomic.fna",
        "gff3": "example/GCF_000005845.2_ASM584v2_genomic.gff",
    },
    "ecoli_example": {
        "name": "Escherichia coli (bundled GCF_008369605.1)",
        "assembly": "GCF_008369605.1",
        "genome": "example/GCF_008369605.1.fna",
        "gff3": "example/GCF_008369605.1.gff",
    },
}


def _merge(dst: dict[str, Any], src: dict[str, Any]) -> dict[str, Any]:
    out = deepcopy(dst)
    for k, v in src.items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = _merge(out[k], v)
        else:
            out[k] = v
    return out


def _as_tuple_pair(value: Any, default: tuple[float, float]) -> tuple[float, float]:
    if value is None:
        return default
    if isinstance(value, (list, tuple)) and len(value) == 2:
        return (float(value[0]), float(value[1]))
    return default


def default_dict() -> dict[str, Any]:
    cfg = PipelineConfig()
    d = {
        "organism": asdict(cfg.organism),
        "genes": cfg.genes,
        "genes_file": cfg.genes_file,
        "cast": asdict(cfg.cast),
        "structure": asdict(cfg.structure),
        "scoring": asdict(cfg.scoring),
        "output": asdict(cfg.output),
    }
    d["cast"]["cds_fraction"] = list(cfg.cast.cds_fraction)
    return d


def load_config(
    path: str | Path | None = None,
    *,
    genes: list[str] | None = None,
    genes_file: str | None = None,
    organism: str | None = None,
    overrides: dict[str, Any] | None = None,
) -> PipelineConfig:
    data = default_dict()
    if path:
        with open(path, encoding="utf-8") as fh:
            loaded = yaml.safe_load(fh) or {}
        data = _merge(data, loaded)

    if organism:
        preset = ORGANISM_PRESETS.get(organism)
        if preset:
            data["organism"] = _merge(data.get("organism") or {}, {**preset, "preset": organism})
        else:
            data["organism"] = _merge(
                data.get("organism") or {},
                {"name": organism, "preset": None},
            )

    if genes_file:
        data["genes_file"] = genes_file
    if genes:
        data["genes"] = list(genes)
    if overrides:
        data = _merge(data, overrides)

    org = data["organism"]
    preset_name = org.get("preset")
    if preset_name and preset_name in ORGANISM_PRESETS:
        org = _merge(ORGANISM_PRESETS[preset_name], org)
        org["preset"] = preset_name

    cast = data["cast"]
    structure = data["structure"]
    scoring = data["scoring"]
    output = data["output"]

    gene_list = [str(g).strip() for g in (data.get("genes") or []) if str(g).strip()]
    gf = data.get("genes_file")
    if gf:
        gf_path = Path(gf)
        if not gf_path.is_absolute():
            # Prefer CWD, then repo root
            if not gf_path.exists():
                alt = REPO_ROOT / gf_path
                if alt.exists():
                    gf_path = alt
        if gf_path.exists():
            for line in gf_path.read_text(encoding="utf-8").splitlines():
                line = line.strip()
                if line and not line.startswith("#"):
                    gene_list.append(line.split(",")[0].strip())

    # de-dupe, keep order
    seen: set[str] = set()
    unique_genes: list[str] = []
    for g in gene_list:
        if g not in seen:
            seen.add(g)
            unique_genes.append(g)

    cds_frac = _as_tuple_pair(cast.get("cds_fraction"), (0.05, 0.50))

    return PipelineConfig(
        organism=OrganismConfig(
            name=org.get("name") or "unknown",
            preset=org.get("preset"),
            taxon_id=org.get("taxon_id"),
            assembly=org.get("assembly"),
            uniprot_proteome=org.get("uniprot_proteome"),
            genome=org.get("genome"),
            gff3=org.get("gff3"),
        ),
        genes=unique_genes,
        genes_file=data.get("genes_file"),
        cast=CastConfig(
            insertion_offset_bp=int(cast.get("insertion_offset_bp", 49)),
            jitter_min_bp=int(cast.get("jitter_min_bp", 44)),
            jitter_max_bp=int(cast.get("jitter_max_bp", 55)),
            guide_length=int(cast.get("guide_length", 32)),
            pam=str(cast.get("pam", "CN")),
            cds_fraction=cds_frac,
            oligo_left=str(cast.get("oligo_left") or DEFAULT_OLIGO_LEFT),
            oligo_right=str(cast.get("oligo_right") or DEFAULT_OLIGO_RIGHT),
        ),
        structure=StructureConfig(
            plddt_min=float(structure.get("plddt_min", 70)),
            rsa_max=float(structure.get("rsa_max", 0.30)),
            pae_split_delta=float(structure.get("pae_split_delta", 5.0)),
            min_domain_aa=int(structure.get("min_domain_aa", 40)),
            edge_buffer_aa=int(structure.get("edge_buffer_aa", 8)),
            prefer_n_terminal_domain=bool(structure.get("prefer_n_terminal_domain", True)),
            skip_disordered=bool(structure.get("skip_disordered", True)),
            burial_cutoff_angstrom=float(structure.get("burial_cutoff_angstrom", 10.0)),
        ),
        scoring=ScoringConfig(
            model=str(scoring.get("model", "model/model.joblib")),
            threshold=float(scoring.get("threshold", 0.5)),
            top_n_per_gene=int(scoring.get("top_n_per_gene", 8)),
            activity_weight=float(scoring.get("activity_weight", 0.55)),
            disruption_weight=float(scoring.get("disruption_weight", 0.45)),
        ),
        output=OutputConfig(
            dir=str(output.get("dir", "results/structure_ko")),
            cache_dir=str(output.get("cache_dir", "cache")),
        ),
    )


def dump_config(cfg: PipelineConfig) -> dict[str, Any]:
    d = {
        "organism": asdict(cfg.organism),
        "genes": cfg.genes,
        "genes_file": cfg.genes_file,
        "cast": asdict(cfg.cast),
        "structure": asdict(cfg.structure),
        "scoring": asdict(cfg.scoring),
        "output": asdict(cfg.output),
    }
    d["cast"]["cds_fraction"] = list(cfg.cast.cds_fraction)
    return d
