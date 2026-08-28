"""Resolve a gene name / locus tag / UniProt ID to a CDS without pasting sequence."""

from __future__ import annotations

import json
import re
import zipfile
from dataclasses import dataclass, asdict
from io import BytesIO
from pathlib import Path
from typing import Any
from urllib.parse import quote
from urllib.request import Request, urlopen

from Bio import SeqIO
from Bio.Seq import Seq

from .config import PipelineConfig, REPO_ROOT

UNIPROT_RE = re.compile(r"UniProtKB(?:/Swiss-Prot|/TrEMBL)?:([A-Z0-9]+)", re.I)
HTTP_HEADERS = {
    "User-Agent": "MultiCAST-structure-ko/0.1 (https://github.com/eulestadt/MultiCAST_guide_predictor)"
}


@dataclass
class GeneRecord:
    query: str
    name: str
    locus_tag: str | None
    seqid: str
    start: int
    end: int
    strand: str
    cds: str
    protein: str
    uniprot: str | None
    protein_id: str | None
    product: str | None
    organism: str

    @property
    def length_nt(self) -> int:
        return len(self.cds)

    @property
    def length_aa(self) -> int:
        # exclude stop codon if present
        aa = len(self.protein.rstrip("*"))
        return aa


def http_json(url: str, timeout: int = 60) -> Any:
    req = Request(url, headers={**HTTP_HEADERS, "Accept": "application/json"})
    with urlopen(req, timeout=timeout) as resp:
        return json.loads(resp.read().decode("utf-8"))


def http_bytes(url: str, timeout: int = 120) -> bytes:
    req = Request(url, headers=HTTP_HEADERS)
    with urlopen(req, timeout=timeout) as resp:
        return resp.read()


def _attr_map(col9: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for part in col9.split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            out[k] = v
    return out


def _first_uniprot(dbxref: str | None) -> str | None:
    if not dbxref:
        return None
    m = UNIPROT_RE.search(dbxref)
    return m.group(1) if m else None


def parse_gff_index(gff3_path: Path) -> dict[str, dict[str, Any]]:
    """Index gene/CDS features by ID, Name, locus_tag, gene, and UniProt."""
    index: dict[str, dict[str, Any]] = {}
    with open(gff3_path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            ftype = fields[2]
            if ftype not in ("gene", "CDS"):
                continue
            attrs = _attr_map(fields[8])
            info: dict[str, Any] = {
                "seqid": fields[0],
                "start": int(fields[3]),
                "end": int(fields[4]),
                "strand": fields[6],
                "type": ftype,
                "ID": attrs.get("ID"),
                "Name": attrs.get("Name"),
                "locus_tag": attrs.get("locus_tag"),
                "gene": attrs.get("gene"),
                "product": attrs.get("product"),
                "protein_id": attrs.get("protein_id"),
                "uniprot": _first_uniprot(attrs.get("Dbxref")),
                "dbxref": attrs.get("Dbxref"),
            }
            keys = [
                attrs.get("ID"),
                attrs.get("Name"),
                attrs.get("locus_tag"),
                attrs.get("gene"),
                info["uniprot"],
                attrs.get("protein_id"),
            ]
            # Prefer CDS over gene when both exist
            for key in keys:
                if not key:
                    continue
                prev = index.get(key)
                if prev is None or (prev.get("type") == "gene" and ftype == "CDS"):
                    index[key] = info
                # case-insensitive alias
                low = key.lower()
                prev_l = index.get(low)
                if prev_l is None or (prev_l.get("type") == "gene" and ftype == "CDS"):
                    index[low] = info
    return index


def _translate_cds(cds: str) -> str:
    seq = Seq(cds)
    # bacterial table 11; incomplete 3' still translates with to_stop=False
    try:
        prot = str(seq.translate(table=11, to_stop=False, cds=False))
    except Exception:
        prot = str(seq.translate(to_stop=False))
    return prot.rstrip("*")


def extract_cds(genome: dict[str, Seq], info: dict[str, Any]) -> str:
    seqid = info["seqid"]
    if seqid not in genome:
        # NCBI FASTA headers are often "NC_000913.3 Escherichia coli ..."
        for k in genome:
            if k == seqid or k.startswith(seqid + " ") or k.split()[0] == seqid:
                seqid = k
                break
        else:
            raise KeyError(f"seqid {info['seqid']} not in genome FASTA")
    frag = genome[seqid][info["start"] - 1 : info["end"]]
    if info["strand"] == "-":
        frag = frag.reverse_complement()
    return str(frag).upper()


def load_genome(fasta_path: Path) -> dict[str, Seq]:
    recs = {}
    for record in SeqIO.parse(str(fasta_path), "fasta"):
        recs[record.id] = record.seq
        recs[record.id.split()[0]] = record.seq
    return recs


def uniprot_search(query: str, cfg: PipelineConfig) -> str | None:
    """Look up a reviewed UniProt accession from a gene name or accession-like token."""
    q = query.strip()
    if re.fullmatch(r"[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2}", q):
        return q
    parts = [f"(gene_exact:{q}) OR (gene:{q}) OR (locus:{q})"]
    if cfg.organism.uniprot_proteome:
        parts.append(f"(proteome:{cfg.organism.uniprot_proteome})")
    elif cfg.organism.taxon_id:
        parts.append(f"(organism_id:{cfg.organism.taxon_id})")
    elif cfg.organism.name:
        parts.append(f"(organism_name:\"{cfg.organism.name}\")")
    uql = " AND ".join(parts)
    url = (
        "https://rest.uniprot.org/uniprotkb/search?"
        f"query={quote(uql)}&fields=accession,reviewed,gene_names&size=5&format=json"
    )
    try:
        payload = http_json(url)
    except Exception:
        return None
    results = payload.get("results") or []
    if not results:
        return None
    reviewed = [r for r in results if r.get("entryType") == "UniProtKB reviewed (Swiss-Prot)"]
    pick = (reviewed or results)[0]
    return pick.get("primaryAccession")


def ensure_genome_files(cfg: PipelineConfig) -> tuple[Path, Path]:
    genome = cfg.resolved_path(cfg.organism.genome)
    gff3 = cfg.resolved_path(cfg.organism.gff3)
    if genome and gff3 and genome.exists() and gff3.exists():
        return genome, gff3

    assembly = cfg.organism.assembly
    if not assembly:
        raise FileNotFoundError(
            "No local genome/GFF found and no assembly accession set. "
            "Add organism.genome / organism.gff3 or organism.assembly to the config."
        )

    cache = cfg.resolved_path(cfg.output.cache_dir) or (REPO_ROOT / "cache")
    dest = cache / "genomes" / assembly
    dest.mkdir(parents=True, exist_ok=True)
    existing_fna = list(dest.glob("*.fna")) + list(dest.glob("*.fa")) + list(dest.glob("*.fasta"))
    existing_gff = list(dest.glob("*.gff")) + list(dest.glob("*.gff3"))
    if existing_fna and existing_gff:
        return existing_fna[0], existing_gff[0]

    url = (
        "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/"
        f"{assembly}/download?include_annotation_type=GENOME_FASTA&include_annotation_type=GENOME_GFF"
    )
    raw = http_bytes(url, timeout=180)
    with zipfile.ZipFile(BytesIO(raw)) as zf:
        for name in zf.namelist():
            lower = name.lower()
            if lower.endswith((".fna", ".fa", ".fasta", ".gff", ".gff3")):
                target = dest / Path(name).name
                target.write_bytes(zf.read(name))
    existing_fna = list(dest.glob("*.fna")) + list(dest.glob("*.fa")) + list(dest.glob("*.fasta"))
    existing_gff = list(dest.glob("*.gff")) + list(dest.glob("*.gff3"))
    if not existing_fna or not existing_gff:
        raise FileNotFoundError(f"NCBI download for {assembly} did not contain FASTA+GFF")
    return existing_fna[0], existing_gff[0]


def resolve_gene(query: str, cfg: PipelineConfig, genome: dict[str, Seq], index: dict[str, dict[str, Any]]) -> GeneRecord:
    q = query.strip()
    info = index.get(q) or index.get(q.lower())
    if info is None:
        raise KeyError(
            f"Gene '{q}' not found in GFF. Use a gene name, locus_tag, protein_id, or UniProt accession."
        )
    cds = extract_cds(genome, info)
    protein = _translate_cds(cds)
    uniprot = info.get("uniprot") or uniprot_search(info.get("gene") or info.get("locus_tag") or q, cfg)
    return GeneRecord(
        query=q,
        name=info.get("gene") or info.get("Name") or q,
        locus_tag=info.get("locus_tag"),
        seqid=str(info["seqid"]),
        start=int(info["start"]),
        end=int(info["end"]),
        strand=str(info["strand"]),
        cds=cds,
        protein=protein,
        uniprot=uniprot,
        protein_id=info.get("protein_id"),
        product=info.get("product"),
        organism=cfg.organism.name,
    )


def gene_to_dict(rec: GeneRecord) -> dict[str, Any]:
    d = asdict(rec)
    d["length_nt"] = rec.length_nt
    d["length_aa"] = rec.length_aa
    d.pop("cds", None)
    d.pop("protein", None)
    d["cds_length"] = rec.length_nt
    d["protein_length"] = rec.length_aa
    return d
