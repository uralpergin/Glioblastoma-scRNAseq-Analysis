# Glioblastoma Single-Cell RNA-seq Analysis Pipeline

Computational pipeline for analyzing single-cell RNA sequencing data from glioblastoma tumors integrated with normal fetal brain tissue. The pipeline performs batch correction, trajectory inference, pathway enrichment, and neural differential equation modeling to understand tumor cell hierarchies and stem cell differentiation.

---

## Datasets

| Dataset | Source | Cells |
|---|---|---|
| Human GBM | GSE84465 (Darmanis et al. 2017, Smart-seq2) | ~3,567 tumor cells |
| Normal fetal brain | HumanFetalBrainPool | ~99,859 cells (subsampled from full dataset) |
| **Combined** | | **~103,426 cells × 16,388 genes** |

---

## Pipeline Overview

```
download.py                  # Download GBM dataset from NCBI GEO
cache_normal_brain.py        # Subsample 100K normal brain cells from HDF5

scvi_training.ipynb          # Batch correction, scVI embedding
dpt_analysis.ipynb           # Diffusion pseudotime, PAGA trajectory
enrichment_analysis.ipynb    # PROGENy, DoRothEA, Hallmark ORA

step1_vae_without_balanced_sampling.ipynb   # Baseline VAE
step2_stem_markers_only.ipynb               # Stem marker cosine similarity
step3_vae_ultra_optimized.ipynb             # VAE + contrastive learning
step3_pure_contrastive_optimized.ipynb      # Pure contrastive encoder

scdiffeq_tumor_analysis.ipynb              # Neural SDE trajectory simulation
concord_tumor_normal_embeddings.ipynb      # CONCORD batch integration
```

---

## Setup

### 1. Clone the repository
```bash
git clone https://github.com/your-username/your-repo-name.git
cd your-repo-name
```

### 2. Create the conda environment
```bash
conda env create -f env.yml
conda activate bio
```

### 3. Download and prepare data

**Tumor dataset** — downloads automatically from NCBI GEO:
```bash
python download.py
```

**Normal brain dataset** — place `HumanFetalBrainPool.h5` in `data/HumanFetalBrainPool/`, then run:
```bash
python cache_normal_brain.py
```

---

## Running the Pipeline

Each notebook has a `LOAD_EXISTING_MODEL` or `TRAIN_MODE` flag at the top — set these to `False` to train from scratch, or `True` to load previously saved checkpoints.

---

## Key Outputs

| File | Description |
|---|---|
| `models/tumor_normal_integrated.h5ad` | Integrated AnnData (shared across all notebooks) |
| `models/scvi_model_improved/` | Trained scVI model |
| `models/best_vae_model.pt` | Trained VAE weights |
| `models/vae_contrastive_model.pt` | VAE + contrastive model weights |
| `outputs/` | Figures and embeddings from each notebook |

> **Note:** Model files (`.pt`, `.ckpt`, `.h5ad`) and raw data are excluded from the repository. See `.gitignore`.

---

## Methods Summary

- **scVI** — variational autoencoder for batch-corrected joint embedding of tumor and normal cells
- **DPT / PAGA** — diffusion pseudotime and graph abstraction for trajectory inference
- **Enrichment** — PROGENy pathway activity, DoRothEA TF activity, Hallmark ORA via decoupler
- **VAE + Contrastive Learning** — custom VAE guided by stem marker similarity to produce biologically-informed embeddings
- **scDiffEq** — neural stochastic differential equations for trajectory simulation from progenitor cells
- **CONCORD** — contrastive learning-based integration as an alternative to scVI

---

