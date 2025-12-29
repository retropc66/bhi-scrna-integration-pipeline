# Quick Start Guide

## Setup

1. **Install dependencies**
```bash
pip install scanpy scvi-tools anndata pandas numpy scipy matplotlib seaborn scikit-learn
pip install scib  # Optional, for advanced benchmarking
```

2. **Set up directory structure**
```bash
mkdir -p output/{anndata,cellassign,embeddings,models,benchmark}
mkdir -p data/markers
```

3. **Update paths in each script**
   - Edit the CONFIG section at the top of each script
   - Update base path from `/home/dcook/projects/def-dcook/active/rare_ov/` to your project directory

## Execution Order

### Step 0: Aggregate (User-provided)
Create `raw.h5ad` with:
- Raw count matrix in `.X`
- `sample_id` column in `.obs`
- See README.md for detailed requirements

### Step 1: Preprocess
```bash
python 01_preprocess.py
```
**Output**: `preprocessed.h5ad` (HVG-subset, immutable)  
**Node**: CPU, ~4x data size memory

### Step 2: CellAssign
```bash
python 02_cellassign.py
```
**Output**: `cellassign/predictions.csv`, `probabilities.csv`, QC plots  
**Node**: CPU or GPU

### Step 3: Integration (GPU)
Run methods independently or in parallel:

```bash
# Option A: Run sequentially
python 03a_integrate_scvi.py
python 03b_integrate_scanvi.py  # Depends on 02 and 03a
python 03c_integrate_mrvi.py
python 03d_integrate_sysvi.py

# Option B: Submit as SLURM jobs (parallel)
sbatch job_03a_scvi.sh
sbatch job_03c_mrvi.sh
sbatch job_03d_sysvi.sh
# scANVI must wait for scVI:
sbatch --dependency=afterok:<scvi_job_id> job_03b_scanvi.sh
```

**Output**: `models/{method}/`, `embeddings/{method}/embedding.npz`  
**Node**: GPU with CUDA

### Step 4: Benchmark
```bash
python 04_benchmark.py
```
**Output**: `benchmark/metrics.csv` and comparison plots  
**Node**: CPU

### Step 5: Finalize
```bash
# Basic usage - use defaults
python 05_finalize.py --method scanvi

# Custom parameters
python 05_finalize.py --method scanvi --leiden-res 0.5 --umap-min-dist 0.5

# Try different methods
python 05_finalize.py --method scvi
python 05_finalize.py --method mrvi
```

**Output**: `anndata/integrated_{method}.h5ad`  
**Node**: CPU, ~4x data size memory

## Common Use Cases

### Quick Test Run (Small Dataset)
```bash
python 01_preprocess.py
python 02_cellassign.py
python 03a_integrate_scvi.py
python 05_finalize.py --method scvi
```

### Full Comparison (All Methods)
```bash
python 01_preprocess.py
python 02_cellassign.py

# Run all integrations
for method in 03a 03b 03c 03d; do
    python ${method}_integrate_*.py
done

python 04_benchmark.py
python 05_finalize.py --method scanvi  # Use best method
```

### Re-cluster with Different Parameters
```bash
# No need to rerun integration, just finalize differently
python 05_finalize.py --method scanvi --leiden-res 0.1  # Fewer clusters
python 05_finalize.py --method scanvi --leiden-res 1.0  # More clusters
```

## Troubleshooting

### "Module not found: scvi"
```bash
pip install scvi-tools
```

### "No marker genes found"
- Check that `cellassign_markers_v2.csv` is in `data/markers/`
- Verify gene names match between marker file and your data (case-sensitive)
- Marker file should be: genes as rows (index), cell types as columns, binary values

### "CUDA out of memory"
- Reduce batch size in training (add to config): `train_size=0.9` 
- Use gradient accumulation
- Run on node with more GPU memory
- Try reducing n_latent or n_layers

### "File not found" during finalize
- Ensure you ran the corresponding integration script (03x)
- Check that embeddings/{method}/embedding.npz exists

### Integration taking too long
- Reduce `max_epochs` (e.g., from 800 to 400)
- Increase `early_stopping_patience` for faster stopping
- For testing, use subset of data

## File Structure After Running

```
project_dir/
├── output/
│   ├── anndata/
│   │   ├── raw.h5ad                    # Step 0 (immutable)
│   │   ├── preprocessed.h5ad           # Step 1 (immutable)
│   │   ├── integrated_scvi.h5ad        # Step 5
│   │   ├── integrated_scanvi.h5ad      # Step 5
│   │   ├── integrated_mrvi.h5ad        # Step 5
│   │   └── integrated_sysvi.h5ad       # Step 5
│   ├── cellassign/                     # Step 2
│   │   ├── predictions.csv
│   │   ├── probabilities.csv
│   │   └── qc_plots/
│   ├── embeddings/                     # Step 3
│   │   ├── scvi/
│   │   ├── scanvi/
│   │   ├── mrvi/
│   │   └── sysvi/
│   ├── models/                         # Step 3
│   │   ├── scvi/
│   │   ├── scanvi/
│   │   ├── mrvi/
│   │   └── sysvi/
│   └── benchmark/                      # Step 4
│       ├── metrics.csv
│       └── *.png
└── data/
    └── markers/
        └── cellassign_markers_v2.csv
```

## Performance Tips

1. **Memory Management**
   - Always use sparse matrices for counts
   - Scripts preserve sparsity throughout
   - Use high-memory nodes for Steps 1 and 5

2. **GPU Efficiency**
   - Set `torch.set_float32_matmul_precision("high")` (already in scripts)
   - Use `gradient_clip_val` to stabilize training (already in scripts)
   - Monitor GPU usage: `nvidia-smi`

3. **Parallelization**
   - Steps 03a, 03c, 03d can run simultaneously
   - Each uses one GPU
   - Can run on separate nodes

4. **Checkpointing**
   - All models auto-saved during training
   - Can resume interrupted training (see scvi-tools docs)
   - All outputs are files - safe to re-run any step

## Next Steps

After creating `integrated_{method}.h5ad`:

```python
import scanpy as sc

# Load integrated object
adata = sc.read_h5ad("output/anndata/integrated_scanvi.h5ad")

# Basic visualization
sc.pl.umap(adata, color=['celltype_pred', 'leiden', 'sample_id'])

# Find marker genes for clusters
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
sc.pl.rank_genes_groups(adata, n_genes=20)

# Differential expression between cell types
sc.tl.rank_genes_groups(adata, 'celltype_pred', method='wilcoxon')

# Export for visualization in other tools
adata.write_csvs("output/for_export/", skip_data=True)
```

## Support

For scvi-tools issues: https://docs.scvi-tools.org/  
For scanpy issues: https://scanpy.readthedocs.io/
