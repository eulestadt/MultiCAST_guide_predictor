"""CLI: design MultiCAST guides from gene names + a YAML config."""

from __future__ import annotations

import argparse
from pathlib import Path

from .config import ORGANISM_PRESETS, load_config
from .snapgene import list_genes


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="python -m structure_ko",
        description=(
            "Design MultiCAST knockout guides from gene names. "
            "Point at a SnapGene .dna file, a preset organism, or FASTA+GFF — "
            "then type gene names (BALF5, lacZ, …) without pasting sequences."
        ),
    )
    sub = p.add_subparsers(dest="command")

    run_p = sub.add_parser("run", help="Design guides (default)")
    _add_run_args(run_p)

    list_p = sub.add_parser(
        "list-genes",
        help="List gene names found in a SnapGene .dna file",
    )
    list_p.add_argument(
        "--dna",
        required=True,
        help="Path to SnapGene .dna / .snapgene file",
    )
    list_p.add_argument(
        "--search",
        help="Optional filter substring (e.g. BALF)",
    )

    # Default command = run (backward compatible flags on root parser)
    _add_run_args(p)
    p.add_argument(
        "--list-presets",
        action="store_true",
        help="Print organism presets and exit.",
    )
    return p


def _add_run_args(p: argparse.ArgumentParser) -> None:
    p.add_argument("--config", "-c", help="YAML config file.")
    p.add_argument(
        "--genes",
        "-g",
        nargs="+",
        help="Gene names (e.g. BALF5 BXLF1 lacZ). No sequences.",
    )
    p.add_argument("--genes-file", "-l", help="One gene name per line.")
    p.add_argument(
        "--dna",
        help="SnapGene .dna file — replaces genome FASTA + GFF3 (e.g. your EBV BAC).",
    )
    p.add_argument(
        "--organism",
        "-o",
        help="Preset: " + ", ".join(sorted(ORGANISM_PRESETS)),
    )
    p.add_argument("--top-n", type=int, help="Top guides per gene.")
    p.add_argument("--outdir", help="Output directory.")


def main(argv: list[str] | None = None) -> None:
    args = build_parser().parse_args(argv)

    if getattr(args, "list_presets", False):
        for name, spec in ORGANISM_PRESETS.items():
            print(f"{name:16}  {spec.get('name')}")
        return

    if args.command == "list-genes":
        path = Path(args.dna).expanduser()
        names = list_genes(path)
        if args.search:
            needle = args.search.upper()
            names = [n for n in names if needle in n.upper()]
        print(f"{len(names)} genes in {path.name}:")
        for n in names:
            print(n)
        return

    # run (explicit subcommand or default root invocation)
    overrides: dict = {}
    if args.top_n:
        overrides["scoring"] = {"top_n_per_gene": args.top_n}
    if args.outdir:
        overrides["output"] = {"dir": args.outdir}
    if args.dna:
        dna = str(Path(args.dna).expanduser())
        overrides["organism"] = {
            "genome": dna,
            "gff3": None,
            "assembly": None,
            "preset": "snapgene",
            "name": Path(dna).stem,
        }
        if args.organism == "ebv_bac" or "EBV" in Path(dna).name.upper():
            overrides["organism"].update(ORGANISM_PRESETS["ebv_bac"])
            overrides["organism"]["genome"] = dna

    cfg = load_config(
        Path(args.config) if args.config else None,
        genes=args.genes,
        genes_file=args.genes_file,
        organism=args.organism,
        overrides=overrides or None,
    )
    from .pipeline import run_pipeline

    run_pipeline(cfg)
