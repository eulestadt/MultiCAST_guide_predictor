"""Run structure-guided MultiCAST guide design for one or many genes."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import yaml
from joblib import load as joblib_load

from MultiCAST_guide_predictor import build_feature_dicts

from .alphafold import fetch_alphafold
from .cast import (
    disruption_score,
    dominant_negative_risk,
    jitter_core_fraction,
    scan_guides,
)
from .config import PipelineConfig, REPO_ROOT, dump_config
from .gene import (
    ensure_genome_files,
    gene_to_dict,
    load_genome,
    parse_gff_index,
    resolve_gene,
)
from .structure import analyze_structure


def _load_model(cfg: PipelineConfig):
    path = cfg.resolved_path(cfg.scoring.model)
    if path is None or not path.exists():
        raise FileNotFoundError(f"Model not found: {cfg.scoring.model}")
    return joblib_load(str(path))


def design_gene(rec, cfg: PipelineConfig, model) -> tuple[pd.DataFrame, dict]:
    entry = fetch_alphafold(rec.uniprot, cfg) if rec.uniprot else None
    analysis = analyze_structure(rec, entry, cfg.structure, cfg.cast.cds_fraction)
    candidates = scan_guides(rec, cfg.cast)
    if not candidates:
        meta = {
            "gene": gene_to_dict(rec),
            "alphafold": None if entry is None else {
                "uniprot": entry.uniprot,
                "source": entry.source,
                "model_created": entry.model_created,
            },
            "structure": _structure_meta(analysis),
            "n_guides": 0,
            "warning": "No CN-PAM + 32-nt windows found in this CDS.",
        }
        return pd.DataFrame(), meta

    feat_df = pd.DataFrame(
        [
            {
                "gene": rec.name,
                "pam_region": c.pam_region,
                "guide_sequence": c.guide_sequence,
                "sequence_number": i,
                "strand": c.strand,
                "center": "Yes"
                if analysis.window.cds_lo_aa <= c.insertion_aa < analysis.window.cds_hi_aa
                else "No",
            }
            for i, c in enumerate(candidates, 1)
        ]
    )
    proba = model.predict_proba(build_feature_dicts(feat_df))[:, 1]

    rows = []
    protein_len = analysis.residue_count
    for c, p in zip(candidates, proba):
        jitter = jitter_core_fraction(c, analysis.window, cfg.cast, protein_len)
        dscore, why = disruption_score(c, analysis, jitter)
        combined = (
            cfg.scoring.activity_weight * float(p)
            + cfg.scoring.disruption_weight * dscore
        )
        plddt = (
            float(analysis.plddt[c.insertion_aa])
            if 0 <= c.insertion_aa < len(analysis.plddt)
            else None
        )
        rsa = (
            float(analysis.rsa_proxy[c.insertion_aa])
            if 0 <= c.insertion_aa < len(analysis.rsa_proxy)
            else None
        )
        rows.append(
            {
                "gene": rec.name,
                "query": rec.query,
                "locus_tag": rec.locus_tag,
                "uniprot": rec.uniprot,
                "product": rec.product,
                "pam_region": c.pam_region,
                "guide_sequence": c.guide_sequence,
                "oligo": c.oligo,
                "strand": c.strand,
                "insertion_nt": c.insertion_nt,
                "insertion_aa": c.insertion_aa,
                "jitter_core_frac": round(jitter, 3),
                "plddt_at_insertion": None if plddt is None else round(plddt, 1),
                "rsa_at_insertion": None if rsa is None else round(rsa, 3),
                "domain_index": None
                if analysis.window.domain is None
                else analysis.window.domain.index,
                "dn_risk": dominant_negative_risk(c.insertion_aa, analysis),
                "in_core": c.insertion_aa in analysis.window.core_residues,
                "in_first_half": analysis.window.cds_lo_aa
                <= c.insertion_aa
                < analysis.window.cds_hi_aa,
                "proba_pos": float(p),
                "pred_label_thr": int(p >= cfg.scoring.threshold),
                "disruption_score": round(dscore, 3),
                "combined_score": round(combined, 3),
                "why": why,
            }
        )

    df = pd.DataFrame(rows)
    df["combined_rank"] = df["combined_score"].rank(method="min", ascending=False).astype(int)
    df["proba_rank"] = df["proba_pos"].rank(method="min", ascending=False).astype(int)
    df = df.sort_values(["combined_score", "proba_pos"], ascending=False).reset_index(drop=True)

    meta = {
        "gene": gene_to_dict(rec),
        "alphafold": None
        if entry is None
        else {
            "uniprot": entry.uniprot,
            "source": entry.source,
            "model_created": entry.model_created,
            "pdb": str(entry.pdb_path),
            "pae": None if entry.pae_path is None else str(entry.pae_path),
        },
        "structure": _structure_meta(analysis),
        "n_guides": int(len(df)),
        "n_core_hits": int(df["in_core"].sum()) if len(df) else 0,
    }
    return df, meta


def _structure_meta(analysis) -> dict:
    window = analysis.window
    return {
        "residue_count": analysis.residue_count,
        "mean_plddt": round(analysis.mean_plddt, 1),
        "n_domains": len(analysis.domains),
        "domains": [
            {
                "index": d.index,
                "start_aa": d.start + 1,
                "end_aa": d.end,
                "length": d.length,
                "mean_plddt": round(d.mean_plddt, 1),
            }
            for d in analysis.domains
        ],
        "window": {
            "method": window.method,
            "target_aa": None if window.target_aa is None else window.target_aa + 1,
            "cds_fraction_aa": [window.cds_lo_aa + 1, window.cds_hi_aa],
            "n_core_residues": len(window.core_residues),
            "domain_index": None if window.domain is None else window.domain.index,
            "notes": window.notes,
        },
    }


def run_pipeline(cfg: PipelineConfig) -> pd.DataFrame:
    if not cfg.genes:
        raise SystemExit("No genes provided. Pass --genes, --genes-file, or config.genes.")

    genome_path, gff_path = ensure_genome_files(cfg)
    genome = load_genome(genome_path)
    index = parse_gff_index(gff_path)
    model = _load_model(cfg)

    outdir = cfg.resolved_path(cfg.output.dir) or (REPO_ROOT / cfg.output.dir)
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "run_config.yaml").write_text(
        yaml.safe_dump(dump_config(cfg), sort_keys=False), encoding="utf-8"
    )

    frames: list[pd.DataFrame] = []
    summaries: list[dict] = []
    missing: list[str] = []

    for query in cfg.genes:
        print(f"[structure_ko] {query}: resolving gene…")
        try:
            rec = resolve_gene(query, cfg, genome, index)
        except Exception as exc:
            print(f"  ! {exc}")
            missing.append(query)
            summaries.append({"query": query, "error": str(exc)})
            continue
        print(
            f"  {rec.name} ({rec.locus_tag or 'no locus'})  "
            f"{rec.length_aa} aa  UniProt={rec.uniprot or 'NA'}"
        )
        print("  fetching AlphaFold DB + scoring MultiCAST guides…")
        df, meta = design_gene(rec, cfg, model)
        summaries.append(meta)
        gene_dir = outdir / rec.name
        gene_dir.mkdir(parents=True, exist_ok=True)
        (gene_dir / "structure.json").write_text(
            json.dumps(meta, indent=2), encoding="utf-8"
        )
        if df.empty:
            print("  no guides found")
            continue
        top = df.head(cfg.scoring.top_n_per_gene)
        df.to_csv(gene_dir / "guides.csv", index=False)
        top.to_csv(gene_dir / "top_guides.csv", index=False)
        frames.append(df)
        print(
            f"  {len(df)} guides, {int(df['in_core'].sum())} in structural core, "
            f"top combined={df.iloc[0]['combined_score']:.3f}"
        )

    (outdir / "summary.json").write_text(json.dumps(summaries, indent=2), encoding="utf-8")
    if missing:
        (outdir / "missing_genes.txt").write_text("\n".join(missing) + "\n", encoding="utf-8")

    if not frames:
        return pd.DataFrame()

    all_df = pd.concat(frames, ignore_index=True)
    all_df.to_csv(outdir / "all_guides.csv", index=False)
    tops = (
        all_df.sort_values(["gene", "combined_score"], ascending=[True, False])
        .groupby("gene", as_index=False)
        .head(cfg.scoring.top_n_per_gene)
    )
    tops.to_csv(outdir / "top_guides.csv", index=False)
    oligo_df = tops[["gene", "guide_sequence", "oligo", "combined_score", "why"]].copy()
    oligo_df.to_csv(outdir / "oligos.csv", index=False)
    print(f"[structure_ko] wrote {outdir}")
    return all_df
