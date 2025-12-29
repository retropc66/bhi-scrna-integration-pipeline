#!/usr/bin/env python

# =========================
# CONFIG
# =========================
PREPROCESSED_H5AD = "/home/dcook/projects/def-dcook/active/rare_ov/output/anndata/preprocessed.h5ad"
OUTPUT_DIR = "/home/dcook/projects/def-dcook/active/rare_ov/output"
MODEL_DIR = f"{OUTPUT_DIR}/models/sysvi"
EMBEDDING_NPZ = f"{OUTPUT_DIR}/embeddings/sysvi/embedding.npz"
BATCH_KEY = "sample_id"
N_LATENT = 30
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
os.makedirs(f"{OUTPUT_DIR}/models/sysvi", exist_ok=True)
os.makedirs(f"{OUTPUT_DIR}/embeddings/sysvi", exist_ok=True)

# Set precision
torch.set_float32_matmul_precision("high")

print("=" * 60)
print("STEP 03D: SysVI INTEGRATION")
print("=" * 60)

# Load preprocessed data
print(f"\n📂 Loading preprocessed data: {PREPROCESSED_H5AD}")
adata = sc.read_h5ad(PREPROCESSED_H5AD)
print(f"   Cells: {adata.n_obs:,}")
print(f"   Genes (HVGs): {adata.n_vars:,}")
print(f"   Batches: {adata.obs[BATCH_KEY].nunique()}")

# Setup SysVI
print("\n⚙️  Setting up SysVI model...")
scvi.model.SYSVI.setup_anndata(
    adata,
    layer="counts",
    batch_key=BATCH_KEY
)

sysvi = scvi.model.SYSVI(
    adata,
    n_latent=N_LATENT,
    n_layers=N_LAYERS
)

print(f"   Latent dimensions: {N_LATENT}")
print(f"   Layers: {N_LAYERS}")
print(f"   System-level variational inference")

# Train model
print("\n🚀 Training SysVI (this may take a while)...")
sysvi.train(
    max_epochs=800,
    early_stopping=True,
    early_stopping_monitor="elbo_validation",
    early_stopping_patience=10
)

# Get latent representation
print("\n📊 Extracting latent representation...")
X_sysvi = sysvi.get_latent_representation()

# Save model
print(f"\n💾 Saving model: {MODEL_DIR}")
sysvi.save(MODEL_DIR, overwrite=True)

# Save embedding as NPZ
print(f"💾 Saving embedding: {EMBEDDING_NPZ}")
np.savez_compressed(
    EMBEDDING_NPZ,
    embedding=X_sysvi.astype(np.float32),
    obs_names=adata.obs_names.to_numpy()
)

# Save metadata
metadata = {
    "method": "SysVI",
    "date": datetime.now().isoformat(),
    "n_cells": int(adata.n_obs),
    "n_genes": int(adata.n_vars),
    "n_batches": int(adata.obs[BATCH_KEY].nunique()),
    "n_latent": N_LATENT,
    "n_layers": N_LAYERS,
    "scvi_tools_version": scvi.__version__
}

with open(f"{OUTPUT_DIR}/embeddings/sysvi/metadata.json", "w") as f:
    json.dump(metadata, f, indent=2)

# Quick UMAP preview for QC
print("\n🗺️  Generating UMAP preview...")
adata.obsm["X_sysvi"] = X_sysvi
sc.pp.neighbors(adata, use_rep="X_sysvi")
sc.tl.umap(adata)

fig, ax = plt.subplots(figsize=(8, 8))
sc.pl.umap(adata, color=BATCH_KEY, ax=ax, show=False, legend_loc='on data', legend_fontsize=6)
ax.set_title("SysVI: UMAP colored by batch")
plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/embeddings/sysvi/umap_preview.png", dpi=150)
plt.close()

print("\n✅ SysVI integration complete!")
print(f"   Model: {MODEL_DIR}")
print(f"   Embedding: {EMBEDDING_NPZ}")
print(f"   Preview: {OUTPUT_DIR}/embeddings/sysvi/umap_preview.png")
