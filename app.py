# app.py — MultiCAST structure-guided guide design portal
from __future__ import annotations

import contextlib
import io
import json
import tempfile
from pathlib import Path

import pandas as pd
import streamlit as st
from joblib import load as joblib_load

from MultiCAST_guide_predictor import build_feature_dicts, extract_guides
from structure_ko.config import ORGANISM_PRESETS, load_config
from structure_ko.genome_bundle import load_genome_bundle
from structure_ko.pipeline import run_pipeline as run_structure_pipeline
from structure_ko.snapgene import load_snapgene

REPO = Path(__file__).parent
MODEL_PATH = REPO / "model" / "model.joblib"

st.set_page_config(
    page_title="MultiCAST Guide Portal",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

if "workdir" not in st.session_state:
    st.session_state.workdir = tempfile.mkdtemp(prefix="multicast_portal_")
WORKDIR = Path(st.session_state.workdir)


def save_upload(uploaded, name: str) -> Path:
    dest = WORKDIR / name
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_bytes(uploaded.getbuffer())
    return dest


@st.cache_resource(show_spinner=False)
def load_model():
    return joblib_load(str(MODEL_PATH))


@st.cache_data(show_spinner="Reading SnapGene annotations…")
def cached_snapgene(path: str) -> tuple[list[str], int]:
    sg = load_snapgene(Path(path))
    return sg.gene_names, len(sg.sequence)


def parse_gene_lines(text: str) -> list[str]:
    genes: list[str] = []
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        genes.append(line.split(",")[0].strip())
    seen: set[str] = set()
    out: list[str] = []
    for g in genes:
        if g and g not in seen:
            seen.add(g)
            out.append(g)
    return out


def build_cfg(
    genes: list[str],
    *,
    preset: str | None,
    genome_path: str | None,
    gff_path: str | None,
    top_n: int,
    threshold: float,
    outdir: Path,
) -> object:
    org: dict = {}
    if genome_path:
        org = {
            "genome": genome_path,
            "gff3": gff_path,
            "assembly": None,
            "name": Path(genome_path).stem,
        }
        if preset == "ebv_bac" or "EBV" in Path(genome_path).name.upper():
            org.update(ORGANISM_PRESETS["ebv_bac"])
            org["genome"] = genome_path
    overrides = {
        "scoring": {"top_n_per_gene": top_n, "threshold": threshold},
        "output": {"dir": str(outdir)},
    }
    if org:
        overrides["organism"] = org
    return load_config(
        genes=genes,
        organism=preset if not genome_path else (preset or None),
        overrides=overrides,
    )


# ── Sidebar ──────────────────────────────────────────────────────────────────
with st.sidebar:
    st.header("Settings")
    top_n = st.number_input("Top guides per gene", 1, 50, 8)
    threshold = st.slider("Activity threshold", 0.0, 1.0, 0.5, 0.01)
    st.divider()
    st.caption("Structure scoring")
    activity_w = st.slider("Activity weight", 0.0, 1.0, 0.55, 0.05)
    st.caption(f"Disruption weight: {1 - activity_w:.2f}")
    st.divider()
    st.markdown(
        "[Paper](https://doi.org/10.1101/2025.10.31.685825) · "
        "[GitHub](https://github.com/eulestadt/MultiCAST_guide_predictor)"
    )

# ── Header ───────────────────────────────────────────────────────────────────
st.title("MultiCAST Guide Portal")
st.markdown(
    "Design **structure-guided** CAST knockout guides from **gene names only**. "
    "Upload a SnapGene `.dna` (e.g. your EBV BAC) or pick a bacterial preset — "
    "no FASTA/GFF hunting required."
)

tab_design, tab_classic, tab_about = st.tabs(
    ["Structure-guided design", "Classic predictor", "How it works"]
)

# ══════════════════════════════════════════════════════════════════════════════
# Structure-guided design
# ══════════════════════════════════════════════════════════════════════════════
with tab_design:
    col_src, col_genes = st.columns([1, 1], gap="large")

    with col_src:
        st.subheader("1 · Genome source")
        source = st.radio(
            "Where are your genes?",
            options=["snapgene", "preset", "fasta_gff"],
            format_func=lambda x: {
                "snapgene": "SnapGene file (.dna) — EBV BAC, plasmid, etc.",
                "preset": "Built-in bacterial genome (E. coli)",
                "fasta_gff": "FASTA + GFF3 upload",
            }[x],
            label_visibility="collapsed",
        )

        preset = None
        genome_path = None
        gff_path = None
        available_genes: list[str] = []

        if source == "snapgene":
            dna_up = st.file_uploader(
                "SnapGene file",
                type=["dna", "snapgene"],
                help="Your annotated .dna file. Gene names are read automatically.",
            )
            if dna_up:
                genome_path = str(save_upload(dna_up, f"upload/{dna_up.name}"))
                available_genes, seq_len = cached_snapgene(genome_path)
                st.success(f"Loaded **{dna_up.name}** — {seq_len:,} bp, **{len(available_genes)}** named CDS")
                preset = st.selectbox(
                    "Organism hint (for UniProt lookup)",
                    options=["ebv_bac", "none"],
                    format_func=lambda x: "EBV / herpesvirus" if x == "ebv_bac" else "Auto / none",
                )
                if preset == "none":
                    preset = None
                else:
                    preset = "ebv_bac"

        elif source == "preset":
            preset = st.selectbox(
                "Preset",
                options=[k for k in ORGANISM_PRESETS if k != "ebv_bac"],
                format_func=lambda k: f"{k} — {ORGANISM_PRESETS[k].get('name', k)}",
            )
            if preset != "ebv_bac":
                try:
                    cfg_probe = load_config(organism=preset)
                    bundle = load_genome_bundle(cfg_probe)
                    available_genes = bundle.gene_names[:500]
                    st.caption(f"{len(bundle.gene_names)} genes in annotation (showing up to 500 in picker).")
                except Exception as exc:
                    st.warning(f"Could not preload gene list: {exc}")

        else:
            c1, c2 = st.columns(2)
            with c1:
                fasta_up = st.file_uploader("Genome FASTA", type=["fna", "fa", "fasta"])
            with c2:
                gff_up = st.file_uploader("GFF3", type=["gff", "gff3"])
            if fasta_up and gff_up:
                genome_path = str(save_upload(fasta_up, "upload/genome.fna"))
                gff_path = str(save_upload(gff_up, "upload/annotation.gff3"))
                st.success("FASTA + GFF3 uploaded.")
                try:
                    cfg_probe = load_config(
                        overrides={"organism": {"genome": genome_path, "gff3": gff_path, "assembly": None}}
                    )
                    bundle = load_genome_bundle(cfg_probe)
                    available_genes = bundle.gene_names[:500]
                except Exception as exc:
                    st.warning(f"Could not parse annotations: {exc}")

    with col_genes:
        st.subheader("2 · Genes to knock out")
        picked: list[str] = []
        if available_genes:
            filter_q = st.text_input("Filter gene list", placeholder="e.g. BALF")
            shown = (
                [g for g in available_genes if filter_q.upper() in g.upper()]
                if filter_q
                else available_genes
            )
            picked = st.multiselect(
                "Pick from file",
                options=shown,
                default=shown[:3] if len(shown) >= 3 and not filter_q else [],
                help="Select one or more genes. You can also type more below.",
            )
        gene_text = st.text_area(
            "Additional genes (one per line)",
            value="",
            height=100,
            placeholder="BALF5\nBXLF1\nBRLF1",
        )
        typed = parse_gene_lines(gene_text)
        genes = list(dict.fromkeys(picked + typed))

        if genes:
            st.info(f"**{len(genes)}** gene(s): {', '.join(genes[:12])}" + (" …" if len(genes) > 12 else ""))
        else:
            st.warning("Select or type at least one gene.")

    st.divider()
    run_design = st.button("Design knockout guides", type="primary", disabled=not genes)

    if run_design and genes:
        outdir = WORKDIR / "structure_ko"
        outdir.mkdir(parents=True, exist_ok=True)
        cfg = build_cfg(
            genes,
            preset=preset,
            genome_path=genome_path,
            gff_path=gff_path,
            top_n=int(top_n),
            threshold=float(threshold),
            outdir=outdir,
        )
        cfg.scoring.activity_weight = float(activity_w)
        cfg.scoring.disruption_weight = float(1.0 - activity_w)

        log = io.StringIO()
        summary: list = []
        try:
            with st.spinner("Resolving genes · AlphaFold · scoring CAST guides…"):
                with contextlib.redirect_stdout(log):
                    df = run_structure_pipeline(cfg)
            if log.getvalue().strip():
                with st.expander("Run log"):
                    st.code(log.getvalue(), language="text")

            summary_path = outdir / "summary.json"
            if summary_path.exists():
                summary = json.loads(summary_path.read_text())
                errors = [s for s in summary if s.get("error")]
                if errors:
                    for e in errors:
                        st.error(f"{e.get('query')}: {e.get('error')}")

            if df is None or df.empty:
                st.warning("No guides returned. Check gene names against your file.")
            else:
                st.success(
                    f"Done — **{df['gene'].nunique()}** gene(s), **{len(df)}** guides scored."
                )

                m1, m2, m3 = st.columns(3)
                m1.metric("Genes", df["gene"].nunique())
                m2.metric("Guides", len(df))
                m3.metric("In structural core", int(df["in_core"].sum()))

                st.subheader("Top guides")
                show_cols = [
                    "gene",
                    "guide_sequence",
                    "strand",
                    "insertion_aa",
                    "in_core",
                    "jitter_core_frac",
                    "proba_pos",
                    "disruption_score",
                    "combined_score",
                    "oligo",
                    "why",
                ]
                st.dataframe(
                    df.sort_values("combined_score", ascending=False)[show_cols].head(100),
                    use_container_width=True,
                    hide_index=True,
                )

                tops = (
                    df.sort_values("combined_score", ascending=False)
                    .groupby("gene", as_index=False)
                    .head(int(top_n))
                )
                oligo_df = tops[["gene", "guide_sequence", "oligo", "combined_score", "why"]]

                d1, d2, d3 = st.columns(3)
                d1.download_button(
                    "all_guides.csv",
                    df.to_csv(index=False).encode(),
                    "all_guides.csv",
                    "text/csv",
                )
                d2.download_button(
                    "oligos.csv",
                    oligo_df.to_csv(index=False).encode(),
                    "oligos.csv",
                    "text/csv",
                )
                if summary_path.exists():
                    d3.download_button(
                        "structure_summary.json",
                        summary_path.read_text().encode(),
                        "summary.json",
                        "application/json",
                    )

                for gene in df["gene"].unique():
                    with st.expander(f"Structure notes — {gene}"):
                        gene_summary = next(
                            (s for s in summary if s.get("gene", {}).get("name") == gene),
                            None,
                        )
                        if gene_summary and "structure" in gene_summary:
                            st.json(gene_summary["structure"])
                        else:
                            st.caption("No structure metadata for this gene.")

        except Exception as exc:
            st.error(str(exc))
            if log.getvalue().strip():
                with st.expander("Run log"):
                    st.code(log.getvalue(), language="text")

# ══════════════════════════════════════════════════════════════════════════════
# Classic predictor
# ══════════════════════════════════════════════════════════════════════════════
with tab_classic:
    st.subheader("Classic MultiCAST activity predictor")
    st.caption("Original tool: sequence features only, no AlphaFold windowing.")

    use_examples = st.toggle("Use bundled example data", value=False)
    example_dir = REPO / "example"

    if use_examples:
        genome_path_c = example_dir / "GCF_008369605.1.fna"
        gff3_path_c = example_dir / "GCF_008369605.1.gff"
        genes_path_c = example_dir / "gene.csv"
    else:
        u1, u2, u3 = st.columns(3)
        with u1:
            g_up = st.file_uploader("FASTA", type=["fna", "fa", "fasta"], key="c_fasta")
        with u2:
            f_up = st.file_uploader("GFF3", type=["gff", "gff3"], key="c_gff")
        with u3:
            l_up = st.file_uploader("Gene CSV", type=["csv"], key="c_csv")
        genome_path_c = save_upload(g_up, "classic/genome.fna") if g_up else None
        gff3_path_c = save_upload(f_up, "classic/annotation.gff3") if f_up else None
        genes_path_c = save_upload(l_up, "classic/genes.csv") if l_up else None

    ready = use_examples or all(p is not None for p in [genome_path_c, gff3_path_c, genes_path_c])
    if st.button("Run classic prediction", disabled=not ready):
        try:
            with st.spinner("Scoring guides…"):
                df_g, missing = extract_guides(
                    str(genome_path_c), str(gff3_path_c), str(genes_path_c)
                )
                if df_g.empty:
                    raise RuntimeError("No guides found.")
                model = load_model()
                proba = model.predict_proba(build_feature_dicts(df_g))[:, 1]
                out = df_g.copy()
                out["proba_pos"] = proba
                out["pred_label_thr"] = (proba >= threshold).astype(int)
            if missing:
                st.warning(f"{len(missing)} gene(s) not found in GFF3.")
            st.dataframe(out.head(80), use_container_width=True)
            st.download_button(
                "predictions.csv",
                out.to_csv(index=False).encode(),
                "predictions.csv",
                "text/csv",
            )
        except Exception as exc:
            st.error(str(exc))

# ══════════════════════════════════════════════════════════════════════════════
# About
# ══════════════════════════════════════════════════════════════════════════════
with tab_about:
    st.markdown(
        """
### What this portal does

1. **You provide gene names** (`BALF5`, `lacZ`, …) — not raw sequences.
2. **Genome + annotations** come from a SnapGene `.dna` file, a bacterial preset, or FASTA+GFF upload.
3. **AlphaFold DB** structures are fetched when available.
4. Guides are ranked so the Tn6677 insertion (~**49 bp** downstream of the protospacer) lands in a **buried fold**, not a linker after a complete domain.
5. The published **MultiCAST XGBoost** model scores transposition activity.

### Gibson oligo format

```
TACTACTGCAAAGTAGCTGATAAC + [32-nt guide] + CTTTACTGCTGAATAAGTAGATAACTAC
```

### Run locally

```bash
pip install -r requirements.txt
./run_portal.sh
# or: streamlit run app.py
```

### EBV BAC example

Upload `EBV-BAC-p2089.dna`, pick `BALF5` / `BXLF1`, run. No FASTA or GFF export needed.

> Activity scores were trained on *E. coli* chromosomal CAST. Treat probabilities as approximate for BAC/plasmid targets until validated in your system.
        """
    )
