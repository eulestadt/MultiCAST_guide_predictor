"""Load genome + gene index from SnapGene, FASTA+GFF, or NCBI assembly."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

from Bio.Seq import Seq

from .config import PipelineConfig, REPO_ROOT
from .gene import ensure_genome_files, load_genome, parse_gff_index
from .snapgene import SnapGeneBundle, load_snapgene

SNAPGENE_SUFFIXES = {".dna", ".snapgene"}


@dataclass
class GenomeBundle:
    source: str
    seqid: str | None
    genome: dict[str, Seq]
    index: dict[str, dict[str, Any]]
    gene_names: list[str]
    snapgene_path: Path | None = None
    fasta_path: Path | None = None
    gff_path: Path | None = None


def _is_snapgene(path: Path) -> bool:
    return path.suffix.lower() in SNAPGENE_SUFFIXES


def load_genome_bundle(cfg: PipelineConfig) -> GenomeBundle:
    genome_path = cfg.resolved_path(cfg.organism.genome)
    gff_path = cfg.resolved_path(cfg.organism.gff3)

    # Mode 1: SnapGene file (annotations built-in — easiest for EBV BAC)
    if genome_path and genome_path.exists() and _is_snapgene(genome_path):
        sg = load_snapgene(genome_path)
        return GenomeBundle(
            source="snapgene",
            seqid=sg.seqid,
            genome=sg.genome,
            index=sg.index,
            gene_names=sg.gene_names,
            snapgene_path=sg.path,
        )

    # Mode 2: local FASTA + GFF
    if genome_path and gff_path and genome_path.exists() and gff_path.exists():
        genome = load_genome(genome_path)
        index = parse_gff_index(gff_path)
        names = sorted(
            {
                str(v.get("gene") or v.get("Name") or "")
                for v in index.values()
                if v.get("gene") or v.get("Name")
            }
        )
        return GenomeBundle(
            source="fasta_gff",
            seqid=None,
            genome=genome,
            index=index,
            gene_names=[n for n in names if n],
            fasta_path=genome_path,
            gff_path=gff_path,
        )

    # Mode 3: download assembly from NCBI
    fasta, gff = ensure_genome_files(cfg)
    genome = load_genome(fasta)
    index = parse_gff_index(gff)
    names = sorted(
        {
            str(v.get("gene") or v.get("Name") or "")
            for v in index.values()
            if v.get("gene") or v.get("Name")
        }
    )
    return GenomeBundle(
        source="ncbi",
        seqid=None,
        genome=genome,
        index=index,
        gene_names=[n for n in names if n],
        fasta_path=fasta,
        gff_path=gff,
    )
