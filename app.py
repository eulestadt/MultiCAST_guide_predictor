# app.py
import tempfile
from pathlib import Path

import streamlit as st
from joblib import load as joblib_load

from MultiCAST_guide_predictor import extract_guides, build_feature_dicts
from structure_ko.config import ORGANISM_PRESETS, load_config
from structure_ko.pipeline import run_pipeline as run_structure_pipeline

st.set_page_config(page_title="MultiCAST Guide Predictor", layout="centered")
st.title("MultiCAST Guide Predictor")
st.caption("Fork with AlphaFold-guided knockout windows. Type a gene name — no FASTA paste.")

if "workdir" not in st.session_state:
    st.session_state["workdir"] = tempfile.mkdtemp(prefix="multicast_")
WORKDIR = Path(st.session_state["workdir"])

tab_struct, tab_classic = st.tabs(["Structure-guided KO", "Classic (FASTA + GFF)"])


def save_upload(uploaded_file, dest_name: str) -> Path:
    out = WORKDIR / dest_name
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "wb") as f:
        f.write(uploaded_file.getbuffer())
    return out


@st.cache_resource(show_spinner=False)
def load_model(model_path: Path):
    return joblib_load(str(model_path))


def run_classic(genome_p: Path, gff3_p: Path, genes_p: Path, model_p: Path, threshold: float):
    df_guides, missing_genes = extract_guides(str(genome_p), str(gff3_p), str(genes_p))
    if df_guides.empty:
        raise RuntimeError("No guides found. Check that genes exist in GFF3 and that PAM sites are present.")
    model = load_model(model_p)
    feats = build_feature_dicts(df_guides)
    proba = model.predict_proba(feats)[:, 1]
    yhat = (proba >= threshold).astype(int)
    out_df = df_guides.copy()
    out_df["proba_pos"] = proba
    out_df["pred_label_thr"] = yhat
    out_df["proba_rank"] = out_df["proba_pos"].rank(method="average", ascending=True)
    return out_df, missing_genes


with tab_struct:
    st.markdown(
        "Looks up the CDS from a bundled (or downloaded) genome, fetches AlphaFold, "
        "and ranks MultiCAST guides so the ~49 bp insertion hits a buried N-terminal fold."
    )
    preset = st.selectbox(
        "Organism preset",
        options=list(ORGANISM_PRESETS.keys()),
        format_func=lambda k: f"{k} — {ORGANISM_PRESETS[k].get('name', k)}",
    )
    gene_text = st.text_area(
        "Genes (one per line)",
        value="lacZ\nrecA",
        help="Gene name, locus_tag, protein_id, or UniProt accession. Not a nucleotide sequence.",
        height=120,
    )
    top_n = st.number_input("Top guides per gene", min_value=1, max_value=50, value=8)
    run_s = st.button("Design knockout guides", type="primary")
    if run_s:
        genes = [ln.strip().split(",")[0] for ln in gene_text.splitlines() if ln.strip() and not ln.strip().startswith("#")]
        if not genes:
            st.error("Enter at least one gene name.")
        else:
            try:
                cfg = load_config(
                    genes=genes,
                    organism=preset,
                    overrides={
                        "scoring": {"top_n_per_gene": int(top_n)},
                        "output": {"dir": str(Path(WORKDIR) / "structure_ko")},
                    },
                )
                with st.spinner("Resolving genes, fetching AlphaFold, scoring guides…"):
                    df = run_structure_pipeline(cfg)
                if df is None or df.empty:
                    st.warning("No guides returned. Check gene names against the preset GFF.")
                else:
                    st.success(f"{df['gene'].nunique()} gene(s), {len(df)} scored guides.")
                    cols = [
                        "gene",
                        "uniprot",
                        "guide_sequence",
                        "strand",
                        "insertion_aa",
                        "in_core",
                        "jitter_core_frac",
                        "proba_pos",
                        "disruption_score",
                        "combined_score",
                        "why",
                    ]
                    show = df[cols].head(80)
                    st.dataframe(show, use_container_width=True)
                    st.download_button(
                        "Download all_guides.csv",
                        df.to_csv(index=False).encode(),
                        file_name="all_guides.csv",
                        mime="text/csv",
                    )
                    oligo = df.sort_values("combined_score", ascending=False).groupby("gene").head(int(top_n))
                    st.download_button(
                        "Download oligos.csv",
                        oligo[["gene", "guide_sequence", "oligo", "combined_score", "why"]].to_csv(index=False).encode(),
                        file_name="oligos.csv",
                        mime="text/csv",
                    )
            except Exception as e:
                st.error(f"{e}")

with tab_classic:
    st.markdown(
        "Original predictor: upload a genome FASTA, GFF3, and a gene-ID list. "
        "No AlphaFold windowing."
    )
    with st.sidebar:
        st.header("Classic settings")
        threshold = st.slider("Decision threshold", 0.0, 1.0, 0.5, 0.01)
        use_examples = st.toggle("Use bundled example data", value=False)

    example_dir = Path(__file__).parent / "example"
    model_path = Path(__file__).parent / "model" / "model.joblib"

    if use_examples:
        genome_path = example_dir / "GCF_008369605.1.fna"
        gff3_path = example_dir / "GCF_008369605.1.gff"
        genes_path = example_dir / "gene.csv"
        missing = [p for p in [genome_path, gff3_path, genes_path, model_path] if not p.exists()]
        if missing:
            st.error(f"Example files missing: {', '.join(str(m) for m in missing)}")
    else:
        c1, c2, c3 = st.columns(3)
        with c1:
            genome_up = st.file_uploader("Genome FASTA (.fna/.fa/.fasta)", type=["fna", "fa", "fasta"])
        with c2:
            gff3_up = st.file_uploader("Annotation GFF3 (.gff/.gff3)", type=["gff", "gff3"])
        with c3:
            genes_up = st.file_uploader("Gene list CSV (one ID per line)", type=["csv"])
        genome_path = save_upload(genome_up, "genome.fna") if genome_up else None
        gff3_path = save_upload(gff3_up, "annotation.gff3") if gff3_up else None
        genes_path = save_upload(genes_up, "genes.csv") if genes_up else None

    can_run = (
        use_examples and all(p and p.exists() for p in [genome_path, gff3_path, genes_path, model_path])
    ) or (
        (not use_examples) and all(p is not None for p in [genome_path, gff3_path, genes_path, model_path])
    )
    run_btn = st.button("Run prediction", disabled=not can_run)
    if run_btn and can_run:
        try:
            with st.spinner("Extracting guides and running model..."):
                preds_df, missing_genes = run_classic(genome_path, gff3_path, genes_path, model_path, threshold)
            if missing_genes:
                st.warning(f"{len(missing_genes)} gene(s) from the CSV were not found in the GFF3.")
                with st.expander("Show missing gene IDs"):
                    st.code("\n".join(str(g) for g in missing_genes), language="text")
            st.success("Done! Preview below.")
            st.dataframe(preds_df.head(51), use_container_width=True)
            st.download_button(
                "Download predictions.csv",
                preds_df.to_csv(index=False).encode(),
                file_name="predictions.csv",
                mime="text/csv",
            )
        except Exception as e:
            st.error(f"Error: {e}")

    with st.expander("Notes"):
        st.markdown(
            """
- Gene list: one identifier per line (`ID`, `Name`, `locus_tag`, or `gene` in the GFF3).
- Prefer the Structure-guided KO tab if you only have gene names.
"""
        )
