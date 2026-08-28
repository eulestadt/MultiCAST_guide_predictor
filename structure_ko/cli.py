"""CLI: design MultiCAST guides from gene names + a YAML config."""

from __future__ import annotations

import argparse
from pathlib import Path

from .config import ORGANISM_PRESETS, load_config
from .pipeline import run_pipeline


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="python -m structure_ko",
        description=(
            "Design MultiCAST knockout guides from gene names. "
            "Looks up CDS from a genome/GFF (bundled or downloaded), fetches AlphaFold "
            "structures, and ranks guides so the ~49 bp insertion lands in a buried fold."
        ),
    )
    p.add_argument(
        "--config",
        "-c",
        help="YAML config with organism, genes, CAST offset, structure filters, and scoring.",
    )
    p.add_argument(
        "--genes",
        "-g",
        nargs="+",
        help="Gene names, locus tags, or UniProt accessions (no sequences required).",
    )
    p.add_argument(
        "--genes-file",
        "-l",
        help="Text/CSV file with one gene identifier per line.",
    )
    p.add_argument(
        "--organism",
        "-o",
        help=(
            "Preset name or organism string. Presets: "
            + ", ".join(sorted(ORGANISM_PRESETS))
        ),
    )
    p.add_argument("--top-n", type=int, help="Override scoring.top_n_per_gene.")
    p.add_argument("--outdir", help="Override output.dir.")
    p.add_argument(
        "--list-presets",
        action="store_true",
        help="Print organism presets and exit.",
    )
    return p


def main(argv: list[str] | None = None) -> None:
    args = build_parser().parse_args(argv)
    if args.list_presets:
        for name, spec in ORGANISM_PRESETS.items():
            print(f"{name:16}  {spec.get('name')}  assembly={spec.get('assembly')}")
        return

    overrides = {}
    if args.top_n:
        overrides["scoring"] = {"top_n_per_gene": args.top_n}
    if args.outdir:
        overrides["output"] = {"dir": args.outdir}

    cfg = load_config(
        Path(args.config) if args.config else None,
        genes=args.genes,
        genes_file=args.genes_file,
        organism=args.organism,
        overrides=overrides or None,
    )
    run_pipeline(cfg)
