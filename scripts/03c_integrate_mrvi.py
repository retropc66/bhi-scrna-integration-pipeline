#!/usr/bin/env python

# =========================
# CONFIG
# =========================
PREPROCESSED_H5AD = "/home/dcook/projects/def-dcook/active/rare_ov/output/anndata/preprocessed.h5ad"
OUTPUT_DIR = "/home/dcook/projects/def-dcook/active/rare_ov/output"
MODEL_DIR = f"{OUTPUT_DIR}/models/mrvi"
EMBEDDING_NPZ = f"{OUTPUT_DIR}/embeddings/mrvi/embedding.npz"
SAMPLE_KEY = "sample_id"     # Target covariate (what you compare)
BATCH_KEY = None             # Set if you have separate technical batches
N_LATENT = 30
# =========================

import scanpy as sc
import numpy as np
import scvi
from scvi.external.mrvi_torch import TorchMRVI as MRVI  # Explicitly use PyTorch backend
import torch
import json
import os
from datetime import datetime
import matplotlib.pyplot as plt

# Create output directories
os.makedirs(f"{OUTPUT_DIR}/models/mrvi", exist_ok=True)
os.makedirs(f"{OUTPUT_DIR}/embeddings/mrvi", exist_ok=True)

# Set precision
torch.set_float32_matmul_precision("high")

print("=" * 60)
print("STEP 03C: MrVI INTEGRATION (PyTorch backend)")
print("=" * 60)

# Load preprocessed data
print(f"\n📂 Loading preprocessed data: {PREPROCESSED_H5AD}")
adata = sc.read_h5ad(PREPROCESSED_H5AD)
print(f"   Cells: {adata.n_obs:,}")
print(f"   Genes (HVGs): {adata.n_vars:,}")
print(f"   Samples: {adata.obs[SAMPLE_KEY].nunique()}")

# Setup MrVI - CRITICAL: use sample_key, not batch_key
print("\n⚙️  Setting up MrVI model...")
MRVI.setup_anndata(
    adata,
    layer="counts",
    sample_key=SAMPLE_KEY,  # Target covariate for sample comparisons
    batch_key=BATCH_KEY     # Optional: technical batch correction
)

mrvi = MRVI(
    adata,
    n_latent=N_LATENT
)

print(f"   Latent dimensions: {N_LATENT}")
print(f"   Samples (target covariate): {adata.obs[SAMPLE_KEY].nunique()}")

# Train model
print("\n🚀 Training MrVI (this may take a while)...")
mrvi.train(
    max_epochs=800,
    early_stopping=True,
    early_stopping_monitor="elbo_validation",
    early_stopping_patience=10
)

# Get BOTH latent representations
print("\n📊 Extracting latent representations...")
u_mrvi = mrvi.get_latent_representation(give_z=False)  # Sample-unaware (integrated)
z_mrvi = mrvi.get_latent_representation(give_z=True)   # Sample-aware

# Save model
print(f"\n💾 Saving model: {MODEL_DIR}")
mrvi.save(MODEL_DIR, overwrite=True)

# Save embeddings as NPZ (both u and z)
print(f"💾 Saving embeddings: {EMBEDDING_NPZ}")
np.savez_compressed(
    EMBEDDING_NPZ,
    u=u_mrvi.astype(np.float32),           # Primary integration embedding
    z=z_mrvi.astype(np.float32),           # Sample-aware embedding
    embedding=u_mrvi.astype(np.float32),   # Alias for consistency with other methods
    obs_names=adata.obs_names.to_numpy()
)

# Save metadata
metadata = {
    "method": "MrVI",
    "backend": "PyTorch",
    "date": datetime.now().isoformat(),
    "n_cells": int(adata.n_obs),
    "n_genes": int(adata.n_vars),
    "n_samples": int(adata.obs[SAMPLE_KEY].nunique()),
    "n_latent": N_LATENT,
    "sample_key": SAMPLE_KEY,
    "batch_key": BATCH_KEY,
    "embeddings_saved": ["u", "z"],
    "primary_embedding": "u",  # For integration benchmarking
    "scvi_tools_version": scvi.__version__
}

with open(f"{OUTPUT_DIR}/embeddings/mrvi/metadata.json", "w") as f:
    json.dump(metadata, f, indent=2)

# Quick UMAP preview for QC (using u embedding)
print("\n🗺️  Generating UMAP preview...")
adata.obsm["X_mrvi_u"] = u_mrvi
adata.obsm["X_mrvi_z"] = z_mrvi
sc.pp.neighbors(adata, use_rep="X_mrvi_u")
sc.tl.umap(adata)

fig, axes = plt.subplots(1, 2, figsize=(14, 6))
sc.pl.umap(adata, color=SAMPLE_KEY, ax=axes[0], show=False, 
           legend_loc='on data', legend_fontsize=6)
axes[0].set_title("MrVI u-space: UMAP colored by sample\n(sample-unaware / integrated)")

# Also show a cell type coloring if available
if "cell_type" in adata.obs.columns:
    sc.pl.umap(adata, color="cell_type", ax=axes[1], show=False)
    axes[1].set_title("MrVI u-space: UMAP colored by cell type")
else:
    sc.pl.umap(adata, color=SAMPLE_KEY, ax=axes[1], show=False)
    axes[1].set_title("MrVI u-space: UMAP")

plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/embeddings/mrvi/umap_preview.png", dpi=150)
plt.close()

print("\n✅ MrVI integration complete!")
print(f"   Model: {MODEL_DIR}")
print(f"   Embedding: {EMBEDDING_NPZ}")
print(f"   - u (sample-unaware): use for integration benchmarking")
print(f"   - z (sample-aware): use for DE/DA analysis")
print(f"   Preview: {OUTPUT_DIR}/embeddings/mrvi/umap_preview.png")
