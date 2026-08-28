# MultiCAST Guide Predictor (structure-guided fork)

Fork of [YiyanYang0728/MultiCAST_guide_predictor](https://github.com/YiyanYang0728/MultiCAST_guide_predictor) with an AlphaFold layer so MultiCAST insertions land in a buried fold instead of a linker that could leave a stable fragment.

Upstream paper: Basta et al., *Adapting CRISPR-associated transposons for rapid and high-throughput reverse genetics* ([bioRxiv 2025.10.31.685825](https://doi.org/10.1101/2025.10.31.685825)).

**Put in a gene name. Do not paste sequences.**

### Your EBV BAC (SnapGene — no FASTA/GFF hunting)

Your `EBV-BAC-p2089.dna` already has ~87 named viral CDS features (`BALF5`, `BXLF1`, …). Point at the file:

```bash
pip install -r requirements.txt

# What genes can I target?
python -m structure_ko list-genes --dna /Users/phoenixwang/Downloads/EBV-BAC-p2089.dna
python -m structure_ko list-genes --dna /Users/phoenixwang/Downloads/EBV-BAC-p2089.dna --search BALF

# Design knockout guides
python -m structure_ko --dna /Users/phoenixwang/Downloads/EBV-BAC-p2089.dna \
  --genes BALF5 BXLF1 BRLF1 --organism ebv_bac
```

Or edit `examples/ebv_bac.yaml` once (genome path) and run `python -m structure_ko --config examples/ebv_bac.yaml`.

### Bacterial presets (auto genome + GFF)

```bash
# one gene
python -m structure_ko --genes lacZ --organism ecoli_k12

# a list
python -m structure_ko --genes-file examples/genes.txt --organism ecoli_k12

# everything in one YAML
python -m structure_ko --config examples/ecoli_k12.yaml
python -m structure_ko --config examples/config.full.yaml
```

Full algorithm write-up: [PIPELINE.md](PIPELINE.md). Every knob: [examples/config.full.yaml](examples/config.full.yaml).

### Web portal

```bash
pip install -r requirements.txt

# macOS only — required once for XGBoost
brew install libomp

./run_portal.sh
# opens http://localhost:8501
```

In the browser:

1. **Structure-guided design** → upload your SnapGene `.dna` (or pick E. coli preset)
2. Pick genes from the dropdown (`BALF5`, …) or type names
3. Download `oligos.csv` with Gibson-ready sequences

Deploy to [Streamlit Community Cloud](https://share.streamlit.io) by pointing at this repo and `app.py` (model file must be in the repo).

## What you get

For each gene, `results/structure_ko/<gene>/`:

| file | contents |
| --- | --- |
| `structure.json` | UniProt, AlphaFold domains, disruption window, notes |
| `guides.csv` | every CN-PAM guide with insertion AA, jitter robustness, activity, disruption, combined score |
| `top_guides.csv` | top N for that gene |

Plus `oligos.csv` with Gibson oligos:

`TACTACTGCAAAGTAGCTGATAAC` + 32-nt guide + `CTTTACTGCTGAATAAGTAGATAACTAC`

## How a gene becomes a guide

1. Resolve `lacZ` / `b0344` / `P00722` against the genome GFF (bundled for *E. coli* K-12).
2. Fetch the AlphaFold DB model (free; no local folding).
3. Split domains from the PAE matrix; take buried, high-pLDDT residues in the N-terminal fold.
4. Back-calculate the CAST guide window: Tn6677 inserts **~49 bp downstream** of the protospacer 3′ end (jitter ~44–55 bp is 1–2 amino acids).
5. Score remaining PAMs with the published MultiCAST XGBoost model.
6. Rank `combined = 0.55 × activity + 0.45 × disruption`.

## Organism presets

```bash
python -m structure_ko --list-presets
```

| preset | genome |
| --- | --- |
| `ecoli_k12` | MG1655 `GCF_000005845.2` (bundled) |
| `ecoli_example` | bundled `GCF_008369605.1` |

For another bacterium, set `organism.assembly` in YAML. FASTA+GFF download once into `cache/`.

---

## Original activity-only predictor

The upstream CLI is unchanged: it still needs a genome FASTA, GFF3, and a gene-ID file.

Online: [multicastguidepredictor-v1.streamlit.app](https://multicastguidepredictor-v1.streamlit.app/)

```bash
python MultiCAST_guide_predictor.py \
  -g example/GCF_008369605.1.fna \
  -f example/GCF_008369605.1.gff \
  -l example/gene.csv \
  -m model/model.joblib \
  -o results/predictions
```

- `-g, --genome`: Genome FASTA
- `-f, --gff3`: Annotation GFF3
- `-l, --genes`: One ID per line (`ID` / `Name` / `locus_tag` / `gene` in the GFF)
- `-m, --model`: `model/model.joblib`
- `-o, --outprefix`: default `results/predictions`
- `--threshold`: default `0.5`

Outputs `proba_pos`, `pred_label_thr`, and `proba_rank`.
