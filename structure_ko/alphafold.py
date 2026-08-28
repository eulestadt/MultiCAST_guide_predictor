"""Fetch AlphaFold DB models (no de novo folding unless the entry is missing)."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any
from urllib.error import HTTPError
from urllib.request import Request, urlopen

from .config import PipelineConfig, REPO_ROOT
from .gene import HTTP_HEADERS, http_json

AFDB_API = "https://alphafold.ebi.ac.uk/api/prediction/{accession}"


@dataclass
class AlphaFoldEntry:
    uniprot: str
    sequence: str
    pdb_path: Path
    pae_path: Path | None
    model_created: str | None
    source: str


def _download(url: str, dest: Path, timeout: int = 120) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    req = Request(url, headers=HTTP_HEADERS)
    with urlopen(req, timeout=timeout) as resp:
        dest.write_bytes(resp.read())


def fetch_alphafold(uniprot: str, cfg: PipelineConfig) -> AlphaFoldEntry | None:
    cache = cfg.resolved_path(cfg.output.cache_dir) or (REPO_ROOT / "cache")
    entry_dir = cache / "alphafold" / uniprot
    meta_path = entry_dir / "prediction.json"
    pdb_path = entry_dir / f"{uniprot}.pdb"
    pae_path = entry_dir / f"{uniprot}_pae.json"

    payload: list[dict[str, Any]] | None = None
    if meta_path.exists() and pdb_path.exists():
        payload = json.loads(meta_path.read_text(encoding="utf-8"))
    else:
        try:
            payload = http_json(AFDB_API.format(accession=uniprot))
        except HTTPError as exc:
            if exc.code == 404:
                return None
            raise
        except Exception:
            return None
        if not payload:
            return None
        entry_dir.mkdir(parents=True, exist_ok=True)
        meta_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    rec = payload[0]
    if not pdb_path.exists():
        pdb_url = rec.get("pdbUrl")
        if not pdb_url:
            return None
        _download(pdb_url, pdb_path)
    if not pae_path.exists():
        pae_url = rec.get("paeDocUrl") or rec.get("bcifUrl")
        # paeDocUrl is the JSON PAE; skip binary cif
        if rec.get("paeDocUrl"):
            try:
                _download(rec["paeDocUrl"], pae_path)
            except Exception:
                pae_path = None  # type: ignore[assignment]
        else:
            pae_path = None  # type: ignore[assignment]

    return AlphaFoldEntry(
        uniprot=uniprot,
        sequence=rec.get("uniprotSequence") or "",
        pdb_path=pdb_path,
        pae_path=pae_path if pae_path and Path(pae_path).exists() else None,
        model_created=rec.get("modelCreatedDate"),
        source="alphafold_db",
    )
