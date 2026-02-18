#!/usr/bin/env python

# =========================
# CONFIG
# =========================
BASEDIR = "/project/rrg-tperkins/OBCF/active/BHI_single_cell_processing/analysis/integration/brain"
SCRIPTDIR = "/project/rrg-tperkins/OBCF/active/BHI_single_cell_processing/bhi-scrna-integration-pipeline"
PREPROCESSED_H5AD = f"{BASEDIR}/output/anndata/preprocessed.h5ad"
CELLASSIGN_PREDICTIONS = f"{BASEDIR}/output/cellassign/predictions.csv"
SCVI_MODEL_DIR = f"{BASEDIR}./output/models/scvi"
OUTPUT_DIR = f"{BASEDIR}/output"
MODEL_DIR = f"{OUTPUT_DIR}/models/scanvi"
EMBEDDING_NPZ = f"{OUTPUT_DIR}/embeddings/scanvi/embedding.npz"
LABEL_KEY = "celltype_pred"
UNLABELED_CATEGORY = "Unknown"
# =========================

import scanpy as sc
import pandas as pd
import numpy as np
import scvi
import torch
import json
import os
from datetime import datetime
import matplotlib.pyplot as plt

# Create output directories
os.makedirs(f"{OUTPUT_DIR}/models/scanvi", exist_ok=True)
os.makedirs(f"{OUTPUT_DIR}/embeddings/scanvi", exist_ok=True)

# Set precision
torch.set_float32_matmul_precision("high")

print("=" * 60)
print("STEP 03B: scANVI INTEGRATION")
print("=" * 60)

# Load preprocessed data
print(f"\n Loading preprocessed data: {PREPROCESSED_H5AD}")
adata = sc.read_h5ad(PREPROCESSED_H5AD)
print(f"  Cells: {adata.n_obs:,}")
print(f"  Genes (HVGs): {adata.n_vars:,}")

# Load CellAssign predictions
print(f"\n Loading CellAssign predictions: {CELLASSIGN_PREDICTIONS}")
predictions = pd.read_csv(CELLASSIGN_PREDICTIONS)
predictions = predictions.set_index('cell_id')

# Add labels to adata
labels = predictions['celltype_pred'].reindex(adata.obs_names)

# Handle missing labels
if labels.isna().any():
    print(f"   Warning: {labels.isna().sum()} cells without predictions, labeling as '{UNLABELED_CATEGORY}'")
    labels = labels.fillna(UNLABELED_CATEGORY)

# Convert to categorical
if UNLABELED_CATEGORY not in labels.values:
    labels = labels.astype('category')
else:
    labels = labels.astype('category')
    if UNLABELED_CATEGORY not in labels.cat.categories:
        labels = labels.cat.add_categories([UNLABELED_CATEGORY])

adata.obs[LABEL_KEY] = labels

print(f"\n Label distribution:")
print(adata.obs[LABEL_KEY].value_counts())

# Load scVI model
print(f"\n Loading scVI model: {SCVI_MODEL_DIR}")
vae = scvi.model.SCVI.load(SCVI_MODEL_DIR, adata=adata)

# Convert to scANVI
print("\n️  Building scANVI from scVI model...")
scanvi = scvi.model.SCANVI.from_scvi_model(
    vae,
    labels_key=LABEL_KEY,
    unlabeled_category=UNLABELED_CATEGORY
)

print(f"  Using labels from: {LABEL_KEY}")
print(f"  Unlabeled category: {UNLABELED_CATEGORY}")

# Train scANVI
print("\n Training scANVI (fine-tuning with labels)...")
scanvi.train(
    max_epochs=800,
    early_stopping=True,
    early_stopping_monitor="elbo_validation",
    plan_kwargs={"lr": 2e-4, "weight_decay": 1e-4},
    gradient_clip_val=1.0
)

# Get latent representation
print("\n Extracting latent representation...")
X_scanvi = scanvi.get_latent_representation()

# Save model
print(f"\n Saving model: {MODEL_DIR}")
scanvi.save(MODEL_DIR, overwrite=True)

# Save embedding as NPZ
print(f" Saving embedding: {EMBEDDING_NPZ}")
np.savez_compressed(
    EMBEDDING_NPZ,
    embedding=X_scanvi.astype(np.float32),
    obs_names=adata.obs_names.to_numpy()
)

# Save metadata
metadata = {
    "method": "scANVI",
    "date": datetime.now().isoformat(),
    "n_cells": int(adata.n_obs),
    "n_genes": int(adata.n_vars),
    "n_batches": int(adata.obs['sample_id'].nunique()) if 'sample_id' in adata.obs else None,
    "n_latent": scanvi.module.n_latent,
    "labels_key": LABEL_KEY,
    "n_cell_types": int(adata.obs[LABEL_KEY].nunique()),
    "scvi_tools_version": scvi.__version__
}

with open(f"{OUTPUT_DIR}/embeddings/scanvi/metadata.json", "w") as f:
    json.dump(metadata, f, indent=2)

# Quick UMAP preview for QC
print("\n️  Generating UMAP preview...")
adata.obsm["X_scanvi"] = X_scanvi
sc.pp.neighbors(adata, use_rep="X_scanvi")
sc.tl.umap(adata)

fig, axes = plt.subplots(1, 2, figsize=(16, 7))
sc.pl.umap(adata, color=LABEL_KEY, ax=axes[0], show=False, legend_loc='right margin')
axes[0].set_title("scANVI: UMAP colored by cell type")
sc.pl.umap(adata, color='sample_id', ax=axes[1], show=False, legend_loc='on data', legend_fontsize=6)
axes[1].set_title("scANVI: UMAP colored by batch")
plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/embeddings/scanvi/umap_preview.png", dpi=150)
plt.close()

print("\n scANVI integration complete!")
print(f"  Model: {MODEL_DIR}")
print(f"  Embedding: {EMBEDDING_NPZ}")
print(f"  Preview: {OUTPUT_DIR}/embeddings/scanvi/umap_preview.png")
