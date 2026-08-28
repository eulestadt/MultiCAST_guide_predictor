# Structure-guided MultiCAST knockout design

This document is the specification for the algorithm: **what it should do**, **why**,
and **how the code implements it**. It extends the MultiCAST guide predictor
(Basta et al., bioRxiv 2025.10.31.685825) with AlphaFold so insertions are
chosen to destroy the fold, not leave a stable N-terminal fragment.

## Why this exists

A MultiCAST mini-transposon in the middle of a bacterial gene almost always
knocks out function. The cargo is kilobases long, the R/L ends contain stop
codons in every frame, and *E. coli* quality control (tmRNA / ClpXP / Lon)
usually degrades the truncated N-terminal peptide.

The remaining failure mode is **not** “the insertion is 44 bp instead of 49 bp”.
A 1–2 amino-acid jitter cannot rescue a domain. The failure mode is **where**
the cut lands:

- If the insertion falls in a **linker after a complete globular domain**, that
  domain can fold on its own, evade proteases, and act as a dominant-negative.
- If the insertion is **late in the protein**, a nearly full-length chain may
  retain catalysis.
- If the insertion is **outside the ORF** (rare long-distance hops, essential
  genes in the paper), the protein is intact.

AlphaFold is used only to pick the **disruption window**. Transposition
activity is still scored by the published MultiCAST XGBoost model.

## What the algorithm should do

```
gene name  →  CDS + UniProt  →  AlphaFold structure
                                    ↓
                         domains · buried core · first half
                                    ↓
                    back-calculate guide window (−49 bp)
                                    ↓
              scan 5′-CN-3′ PAMs  →  XGBoost activity
                                    ↓
              rank by activity × structural disruption
                                    ↓
              oligo = MultiCAST homology arms + 32-nt guide
```

### 1. Take a gene, not a sequence

The user should type `lacZ` or `b0344`, not paste a FASTA.

Resolution order:

1. Match `ID`, `Name`, `locus_tag`, `gene`, `protein_id`, or UniProt in the GFF.
2. Extract the CDS from the genome FASTA (reverse-complement if needed).
3. Read `Dbxref=UniProtKB/Swiss-Prot:…` when present; otherwise search UniProt
   with the gene name and the configured proteome / taxon.

For *E. coli* K-12 the genome and GFF are already in `example/`. For another
species, set `organism.assembly` and the pipeline downloads FASTA+GFF once.

### 2. Fetch a structure (do not fold 80 genes from scratch)

If the protein is in the AlphaFold Database, download the PDB (pLDDT is the
B-factor) and the PAE JSON. That is free and takes seconds. Only missing
entries need ColabFold later; the code records a first-half-CDS fallback
rather than blocking the run.

### 3. Partition domains from PAE, not from “looks like a helix”

- **PAE split:** recursively split the chain where the mean PAE *between*
  two blocks is high relative to the PAE *inside* each block, preferring
  splits in low-pLDDT linkers.
- **pLDDT core:** residues with pLDDT ≥ 70 are ordered.
- **Burial:** CA atoms with many neighbors within 8 Å are treated as buried
  (RSA proxy). Exposed loops are weaker disruption sites than the hydrophobic
  core.
- **N-terminal preference:** cut **inside** the first well-folded domain so
  that domain cannot persist as a fragment. Do not cut in the linker after it.

### 4. Convert a residue into a CAST guide window

Tn6677 does **not** insert in the protospacer. Cascade binds a 32-nt target
next to a 5′-CN-3′ PAM; TnsABC inserts **~49 bp downstream of the protospacer
3′ end**.

| Target strand | Insertion on the CDS |
| --- | --- |
| Coding | `guide_3′ + 49 bp` (toward the C-terminus) |
| Template | `guide_3′ − 49 bp` (toward the N-terminus) |

Jitter of ~44–55 bp is 1–2 amino acids. The algorithm therefore:

- aims at the **center** of a core secondary-structure block, not a domain edge
- keeps an `edge_buffer_aa` (default 8 residues)
- scores `jitter_core_frac`: fraction of landings in 44–55 bp that still hit core

A ±2 aa shift cannot save a 100–300 residue domain. It *can* miss a 4-residue
catalytic motif if you aimed at the motif edge — which is why we aim at cores.

### 5. Score activity, then combine

Every CN-PAM + 32-nt window is featurized and scored with the published
MultiCAST model (`model/model.joblib`; AUROC 0.87 on held-out genes).

```
combined = 0.55 * P(high activity) + 0.45 * disruption_score
```

`disruption_score` rewards buried core, N-terminal domain, and jitter
robustness; it penalizes intact upstream domains (dominant-negative risk)
and C-terminal insertions.

### 6. Emit order-ready oligos

From the MultiCAST methods:

```
5′-TACTACTGCAAAGTAGCTGATAAC-[32 nt guide]-CTTTACTGCTGAATAAGTAGATAACTAC-3′
```

The flanking arms are atypical CRISPR repeats plus homology for Gibson
assembly into PacI/SpeI-linearized mini-Tn donor.

## What this does *not* claim

- It does not simulate tmRNA tagging or proteasome kinetics.
- It does not replace experimental western blots. It only lowers the chance
  of a stable fragment relative to “any PAM in the first half of the gene”.
- It does not model H-NS occlusion. The original paper’s H-NS warning still
  applies for pooled screens.

## Mapping to code

| Step | Module |
| --- | --- |
| Config + gene list | `structure_ko/config.py` |
| Name → CDS + UniProt | `structure_ko/gene.py` |
| AlphaFold DB fetch | `structure_ko/alphafold.py` |
| PAE domains, burial, window | `structure_ko/structure.py` |
| +49 bp offset, jitter, oligos | `structure_ko/cast.py` |
| XGBoost + ranking | `structure_ko/pipeline.py` (uses `MultiCAST_guide_predictor.py`) |
| CLI | `python -m structure_ko` |

## Suggested defaults for ~80 knockouts

Leave the YAML defaults. Cost is $0 if proteins are in AlphaFold DB.
Runtime is dominated by one HTTP fetch per gene plus a few seconds of
NumPy on the PAE matrix. The XGBoost step is milliseconds per guide.
