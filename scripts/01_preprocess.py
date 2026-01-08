#!/usr/bin/env python

# =========================
# CONFIG
# =========================
RAW_H5AD = "/home/dcook/projects/def-dcook/active/rare_ov/output/anndata/raw.h5ad"
PREPROCESSED_H5AD = "/home/dcook/projects/def-dcook/active/rare_ov/output/anndata/preprocessed.h5ad"
N_TOP_GENES = 3000
BATCH_KEY = "sample_id"
# =========================

import scanpy as sc

print("=" * 60)
print("STEP 01: PREPROCESSING")
print("=" * 60)

print(f"\n📂 Loading raw data: {RAW_H5AD}")
adata = sc.read_h5ad(RAW_H5AD)
print(f"   Cells: {adata.n_obs:,}")
print(f"   Genes: {adata.n_vars:,}")
print(f"   Batches: {adata.obs[BATCH_KEY].nunique()}")

# Preserve raw counts in a layer
print("\n📋 Creating counts layer...")
adata.layers["counts"] = adata.X.copy()

# Select highly variable genes (batch-aware, uses counts layer)
print(f"\n🔍 Identifying {N_TOP_GENES} highly variable genes (batch-aware)...")
sc.pp.highly_variable_genes(
    adata,
    n_top_genes=N_TOP_GENES,
    batch_key=BATCH_KEY,
    flavor="seurat_v3",
    layer="counts",
    subset=False  # Don't subset yet, just mark them
)

n_hvg = adata.var["highly_variable"].sum()
print(f"   Found {n_hvg} HVGs")

# Subset to HVGs
print("\n✂️  Subsetting to HVGs...")
adata = adata[:, adata.var["highly_variable"]].copy()
print(f"   Dimensions: {adata.n_obs:,} cells × {adata.n_vars:,} genes")

# Normalize (standard scvi-tools convention: X=normalized, layers["counts"]=raw)
print("\n📊 Normalizing...")
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Save preprocessed data
print(f"\n💾 Saving preprocessed data: {PREPROCESSED_H5AD}")
adata.write_h5ad(PREPROCESSED_H5AD, compression="gzip")

print("\n✅ Preprocessing complete!")
print(f"   Output: {PREPROCESSED_H5AD}")
print(f"   adata.X = log-normalized counts")
print(f"   adata.layers['counts'] = raw counts")
