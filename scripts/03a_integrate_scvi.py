#!/usr/bin/env python

# =========================
# CONFIG
# =========================
PREPROCESSED_H5AD = "/home/dcook/projects/def-dcook/active/rare_ov/output/anndata/preprocessed.h5ad"
OUTPUT_DIR = "/home/dcook/projects/def-dcook/active/rare_ov/output"
MODEL_DIR = f"{OUTPUT_DIR}/models/scvi"
EMBEDDING_NPZ = f"{OUTPUT_DIR}/embeddings/scvi/embedding.npz"
BATCH_KEY = "sample_id"
N_LATENT = 20
N_LAYERS = 2
# =========================

import scanpy as sc
import numpy as np
import scvi
import torch
import json
import os
from datetime import datetime
import matplotlib.pyplot as plt

# Create output directories
os.makedirs(f"{OUTPUT_DIR}/models/scvi", exist_ok=True)
os.makedirs(f"{OUTPUT_DIR}/embeddings/scvi", exist_ok=True)

# Set precision for faster training
torch.set_float32_matmul_precision("high")

print("=" * 60)
print("STEP 03A: scVI INTEGRATION")
print("=" * 60)

# Load preprocessed data
print(f"\n📂 Loading preprocessed data: {PREPROCESSED_H5AD}")
adata = sc.read_h5ad(PREPROCESSED_H5AD)
print(f"   Cells: {adata.n_obs:,}")
print(f"   Genes (HVGs): {adata.n_vars:,}")
print(f"   Batches: {adata.obs[BATCH_KEY].nunique()}")

# Setup scVI
print("\n⚙️  Setting up scVI model...")
scvi.model.SCVI.setup_anndata(
    adata,
    layer="counts",
    batch_key=BATCH_KEY
)

vae = scvi.model.SCVI(
    adata,
    n_layers=N_LAYERS,
    n_latent=N_LATENT,
    gene_likelihood="nb"
)

print(f"   Latent dimensions: {N_LATENT}")
print(f"   Layers: {N_LAYERS}")
print(f"   Gene likelihood: negative binomial")

# Train model
print("\n🚀 Training scVI (this may take a while)...")
vae.train(
    max_epochs=400,
    early_stopping=True,
    early_stopping_monitor="elbo_validation",
    early_stopping_patience=10
)

# Get latent representation
print("\n📊 Extracting latent representation...")
X_scvi = vae.get_latent_representation()

# Save model
print(f"\n💾 Saving model: {MODEL_DIR}")
vae.save(MODEL_DIR, overwrite=True)

# Save embedding as NPZ
print(f"💾 Saving embedding: {EMBEDDING_NPZ}")
np.savez_compressed(
    EMBEDDING_NPZ,
    embedding=X_scvi.astype(np.float32),
    obs_names=adata.obs_names.to_numpy()
)

# Save metadata
metadata = {
    "method": "scVI",
    "date": datetime.now().isoformat(),
    "n_cells": int(adata.n_obs),
    "n_genes": int(adata.n_vars),
    "n_batches": int(adata.obs[BATCH_KEY].nunique()),
    "n_latent": N_LATENT,
    "n_layers": N_LAYERS,
    "gene_likelihood": "nb",
    "scvi_tools_version": scvi.__version__
}

with open(f"{OUTPUT_DIR}/embeddings/scvi/metadata.json", "w") as f:
    json.dump(metadata, f, indent=2)

# Quick UMAP preview for QC
print("\n🗺️  Generating UMAP preview...")
adata.obsm["X_scvi"] = X_scvi
sc.pp.neighbors(adata, use_rep="X_scvi")
sc.tl.umap(adata)

fig, ax = plt.subplots(figsize=(8, 8))
sc.pl.umap(adata, color=BATCH_KEY, ax=ax, show=False, legend_loc='on data', legend_fontsize=6)
ax.set_title("scVI: UMAP colored by batch")
plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/embeddings/scvi/umap_preview.png", dpi=150)
plt.close()

print("\n✅ scVI integration complete!")
print(f"   Model: {MODEL_DIR}")
print(f"   Embedding: {EMBEDDING_NPZ}")
print(f"   Preview: {OUTPUT_DIR}/embeddings/scvi/umap_preview.png")
