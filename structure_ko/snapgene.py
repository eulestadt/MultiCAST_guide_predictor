"""Read SnapGene .dna files — sequence + gene names without exporting FASTA/GFF."""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from Bio import SeqIO
from Bio.Seq import Seq

EBV_GENE_RE = re.compile(r"\b([A-Z][A-Z0-9]{1,5}(?:\.[0-9]+)?)\b")


def _first(vals: list[str] | None) -> str | None:
    if not vals:
        return None
    return vals[0]


def _gene_tokens_from_feature(qualifiers: dict[str, list[str]]) -> list[str]:
    """Pull EBV-style gene symbols from SnapGene feature qualifiers."""
    tokens: list[str] = []
    for key in ("name", "label", "gene", "product", "note"):
        raw = _first(qualifiers.get(key))
        if not raw:
            continue
        if key == "name" and raw.endswith(" CDS"):
            tokens.append(raw[:-4].strip())
        if key == "label" and raw.endswith(" mRNA"):
            tokens.append(raw[:-5].strip())
        for m in EBV_GENE_RE.finditer(raw.upper()):
            sym = m.group(1)
            if sym not in {"CDS", "RNA", "DNA", "MRNA"} and len(sym) >= 3:
                tokens.append(sym)
        if key in ("name", "label") and raw and not raw.startswith("CDS_"):
            tokens.append(raw.split()[0])
    # de-dupe, preserve order
    seen: set[str] = set()
    out: list[str] = []
    for t in tokens:
        t = t.strip()
        if not t:
            continue
        for variant in (t, t.upper(), t.lower()):
            if variant not in seen:
                seen.add(variant)
                out.append(variant)
    return out


def _location_bounds(loc) -> tuple[int, int, str]:
    """Return 1-based inclusive start/end and strand from a Biopython location."""
    strand = "+" if loc.strand >= 0 else "-"
    parts = []
    if hasattr(loc, "parts"):
        parts = list(loc.parts)
    else:
        parts = [loc]
    start = min(int(p.start) for p in parts) + 1
    end = max(int(p.end) for p in parts)
    return start, end, strand


@dataclass
class SnapGeneBundle:
    """In-memory genome + gene index from a SnapGene file."""

    path: Path
    seqid: str
    sequence: Seq
    index: dict[str, dict[str, Any]]
    gene_names: list[str]

    @property
    def genome(self) -> dict[str, Seq]:
        return {self.seqid: self.sequence}


def load_snapgene(path: Path) -> SnapGeneBundle:
    path = path.expanduser().resolve()
    record = SeqIO.read(str(path), "snapgene")
    seqid = record.id if record.id and record.id != "-" else path.stem
    index: dict[str, dict[str, Any]] = {}
    display_names: list[str] = []

    for feat in record.features:
        if feat.type != "CDS":
            continue
        start, end, strand = _location_bounds(feat.location)
        q = feat.qualifiers
        gene_name = None
        for candidate in (_first(q.get("name")), _first(q.get("label"))):
            if candidate and candidate.endswith(" CDS"):
                gene_name = candidate[:-4].strip()
                break
            if candidate and not candidate.startswith("CDS_"):
                gene_name = candidate.split()[0]
                break
        if not gene_name:
            note = _first(q.get("note")) or ""
            m = re.search(r"^([A-Z][A-Z0-9]{1,5}(?:\.[0-9]+)?)\b", note)
            gene_name = m.group(1) if m else _first(q.get("label"))

        transl_table = int(_first(q.get("transl_table")) or "1")
        info: dict[str, Any] = {
            "seqid": seqid,
            "start": start,
            "end": end,
            "strand": strand,
            "type": "CDS",
            "ID": _first(q.get("label")),
            "Name": gene_name,
            "locus_tag": _first(q.get("label")),
            "gene": gene_name,
            "product": _first(q.get("note")),
            "protein_id": None,
            "uniprot": None,
            "transl_table": transl_table,
        }
        keys = _gene_tokens_from_feature(q) + ([gene_name] if gene_name else [])
        for key in keys:
            if not key:
                continue
            prev = index.get(key)
            if prev is None:
                index[key] = info
            low = key.lower()
            if low not in index:
                index[low] = info
        if gene_name and gene_name not in display_names:
            display_names.append(gene_name)

    display_names.sort()
    return SnapGeneBundle(
        path=path,
        seqid=seqid,
        sequence=record.seq,
        index=index,
        gene_names=display_names,
    )


def list_genes(path: Path) -> list[str]:
    return load_snapgene(path).gene_names
