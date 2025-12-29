#!/usr/bin/env python

# =========================
# CONFIG
# =========================
PREPROCESSED_H5AD = "/home/dcook/projects/def-dcook/active/rare_ov/output/anndata/preprocessed.h5ad"
OUTPUT_DIR = "/home/dcook/projects/def-dcook/active/rare_ov/output"
EMBEDDING_NPZ = f"{OUTPUT_DIR}/embeddings/harmony/embedding.npz"
BATCH_KEY = "sample_id"
N_PCS = 50
# =========================

import scanpy as sc
import numpy as np
import json
import os
from datetime import datetime
import matplotlib.pyplot as plt

try:
    from harmonypy import run_harmony
except ImportError:
    raise ImportError("harmonypy not installed. Install with: pip install harmonypy")

os.makedirs(f"{OUTPUT_DIR}/embeddings/harmony", exist_ok=True)

print("=" * 60)
print("STEP 03E: HARMONY INTEGRATION")
print("=" * 60)

# Load preprocessed data (already log-normalized)
print(f"\n📂 Loading preprocessed data: {PREPROCESSED_H5AD}")
adata = sc.read_h5ad(PREPROCESSED_H5AD)
print(f"   Cells: {adata.n_obs:,}")
print(f"   Genes (HVGs): {adata.n_vars:,}")
print(f"   Batches: {adata.obs[BATCH_KEY].nunique()}")

# Scale for PCA
print("\n📊 Scaling...")
sc.pp.scale(adata, max_value=10)

# Run PCA
print(f"🔬 Running PCA (n_comps={N_PCS})...")
sc.tl.pca(adata, n_comps=N_PCS)

# Run Harmony
print(f"\n🚀 Running Harmony on {N_PCS} PCs...")
harmony_out = run_harmony(
    adata.obsm['X_pca'],
    adata.obs,
    BATCH_KEY,
    max_iter_harmony=20
)

X_harmony = harmony_out.Z_corr.T  # Transpose to cells x dims

print(f"   Output shape: {X_harmony.shape}")

# Save embedding as NPZ
print(f"\n💾 Saving embedding: {EMBEDDING_NPZ}")
np.savez_compressed(
    EMBEDDING_NPZ,
    embedding=X_harmony.astype(np.float32),
    obs_names=adata.obs_names.to_numpy()
)

# Save metadata
metadata = {
    "method": "Harmony",
    "date": datetime.now().isoformat(),
    "n_cells": int(adata.n_obs),
    "n_genes": int(adata.n_vars),
    "n_batches": int(adata.obs[BATCH_KEY].nunique()),
    "n_pcs": N_PCS,
    "max_iter_harmony": 20
}

with open(f"{OUTPUT_DIR}/embeddings/harmony/metadata.json", "w") as f:
    json.dump(metadata, f, indent=2)

# Quick UMAP preview for QC
print("\n🗺️  Generating UMAP preview...")
adata.obsm["X_harmony"] = X_harmony
sc.pp.neighbors(adata, use_rep="X_harmony")
sc.tl.umap(adata)

fig, ax = plt.subplots(figsize=(8, 8))
sc.pl.umap(adata, color=BATCH_KEY, ax=ax, show=False, legend_loc='on data', legend_fontsize=6)
ax.set_title("Harmony: UMAP colored by batch")
plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/embeddings/harmony/umap_preview.png", dpi=150)
plt.close()

print("\n✅ Harmony integration complete!")
print(f"   Embedding: {EMBEDDING_NPZ}")
print(f"   Preview: {OUTPUT_DIR}/embeddings/harmony/umap_preview.png")
