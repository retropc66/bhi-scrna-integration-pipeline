# Environment Setup on DRAC Nibi Cluster

This guide sets up a Python environment for the scRNA-seq integration pipeline on DRAC's Nibi cluster (H100 GPUs).

## Key Principles

1. **Use virtualenv, NOT conda** — Conda causes CUDA compatibility issues on DRAC
2. **Match CUDA versions** — PyTorch CUDA version must match the loaded CUDA module
3. **Use DRAC wheels when available** — Use `--no-index` for optimized packages
4. **Environment location** — Create in `$HOME` or `/project` (not `$SLURM_TMPDIR`)

---

## Step 1: Check Available Modules

First, SSH into Nibi and check available Python and CUDA versions:

```bash
ssh <username>@nibi.alliancecan.ca

# Check Python versions
module spider python

# Check CUDA versions  
module spider cuda

# Check what's recommended
module spider StdEnv
```

## Step 2: Load Modules

Load the standard environment with Python 3.11 (recommended for scvi-tools compatibility):

```bash
module load StdEnv/2023 python/3.11
```

Verify:
```bash
python --version  # Should show 3.11.x
which python      # Should be in /cvmfs/...
```

## Step 3: Create Virtual Environment

Create a persistent virtualenv in your project directory:

```bash
# Navigate to your project
cd /home/dcook/projects/def-dcook/active/rare_ov

# Create environment
virtualenv --no-download ~/envs/scrna_integration

# Activate it
source ~/envs/scrna_integration/bin/activate

# Upgrade pip (use DRAC wheel)
pip install --no-index --upgrade pip
```

## Step 4: Install PyTorch with Correct CUDA Version

**Critical**: Install PyTorch matching the CUDA module you'll use in jobs.

Check which CUDA version Nibi has:
```bash
module spider cuda
```

Install PyTorch (assuming CUDA 12.x is available on Nibi):

```bash
# Option A: Use DRAC wheel if available (preferred)
pip install --no-index torch torchvision torchaudio

# Option B: If DRAC wheel has wrong CUDA version, install from PyTorch directly
# Match this to your CUDA module version (e.g., cuda/12.4 → cu124)
# pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu124
```

Verify PyTorch CUDA version:
```bash
python -c "import torch; print(f'PyTorch CUDA: {torch.version.cuda}')"
```

## Step 5: Install Core Dependencies

Install packages available as DRAC wheels first:

```bash
# Core scientific packages (DRAC optimized)
pip install --no-index numpy scipy pandas scikit-learn

# These may or may not be in DRAC wheels - try --no-index first
pip install --no-index matplotlib seaborn h5py
```

## Step 6: Install scRNA-seq Specific Packages

These are typically NOT in DRAC wheels, install from PyPI:

```bash
# scvi-tools and dependencies
pip install scvi-tools

# scanpy and anndata
pip install scanpy anndata harmonypy

# Optional: scib for advanced benchmarking
pip install scib
```

---

## Step 7: MrVI Backend Note

MrVI in scvi-tools 1.4+ supports both JAX and PyTorch backends. **We use the PyTorch backend (`TorchMRVI`)** which runs in the same `scrna_integration` environment as scVI and scANVI. This avoids the complexity of managing separate JAX/cuDNN dependencies.

The script `03c_integrate_mrvi.py` explicitly imports:
```python
from scvi.external.mrvi_torch import TorchMRVI as MRVI
```

No separate environment is needed.

## Step 8: Verify PyTorch Installation

Create a test script to verify GPU access:

```bash
cat > ~/test_gpu.py << 'EOF'
import torch
print(f"PyTorch version: {torch.__version__}")
print(f"CUDA available: {torch.cuda.is_available()}")
print(f"PyTorch CUDA version: {torch.version.cuda}")
if torch.cuda.is_available():
    print(f"GPU count: {torch.cuda.device_count()}")
    print(f"GPU name: {torch.cuda.get_device_name(0)}")

import scvi
print(f"scvi-tools version: {scvi.__version__}")

import scanpy as sc
print(f"scanpy version: {sc.__version__}")
EOF
```

**Important**: GPU tests must run on a compute node, not login node.

## Step 9: Test on GPU Node

Request an interactive GPU session:

```bash
salloc --account=def-dcook --time=0:30:00 --gres=gpu:h100:1 --mem=32G --cpus-per-task=4

# Once allocated, load modules and activate env
module load StdEnv/2023 python/3.11 cuda/12.9
source ~/envs/scrna_integration/bin/activate

# Run test
python ~/test_gpu.py
```

If `torch.cuda.is_available()` returns `False`, see Troubleshooting below.

---

## Complete requirements.txt

For reproducibility, here's the full requirements:

```
# requirements.txt
# Install with: pip install -r requirements.txt

# Core (try --no-index first)
numpy
scipy
pandas
scikit-learn
matplotlib
seaborn
h5py

# scRNA-seq specific (from PyPI)
scanpy>=1.9.0
anndata>=0.9.0
scvi-tools>=1.0.0

# Optional benchmarking
# scib>=1.1.0
```

---

## Job Submission Template

Example SLURM script for GPU jobs (`job_gpu.sh`):

```bash
#!/bin/bash
#SBATCH --account=def-dcook
#SBATCH --job-name=scvi_train
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --gres=gpu:h100:1
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

# Load modules - CUDA version must match PyTorch
module load StdEnv/2023 python/3.11 cuda/12.9

# Activate environment
source ~/envs/scrna_integration/bin/activate

# Verify GPU (optional)
nvidia-smi
python -c "import torch; print(f'CUDA available: {torch.cuda.is_available()}')"

# Run your script
python 03a_integrate_scvi.py
```

Example for CPU jobs (`job_cpu.sh`):

```bash
#!/bin/bash
#SBATCH --account=def-dcook
#SBATCH --job-name=preprocess
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=4:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

module load StdEnv/2023 python/3.11

source ~/envs/scrna_integration/bin/activate

python 01_preprocess.py
```

---

## Troubleshooting

### `torch.cuda.is_available()` returns `False`

**Cause**: CUDA version mismatch between PyTorch and loaded module.

**Fix**:
```bash
# Check versions
module list  # Note CUDA version
python -c "import torch; print(torch.version.cuda)"  # Note PyTorch CUDA

# If mismatch, reinstall PyTorch with correct CUDA
pip uninstall torch torchvision torchaudio -y

# For CUDA 12.9 (Nibi):
pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu124

# Note: PyTorch doesn't have cu129 wheels yet, but cu124 is forward-compatible
# with CUDA 12.9 drivers. The driver is backward-compatible with older toolkit versions.
```

### `ModuleNotFoundError` for packages

Ensure environment is activated:
```bash
source ~/envs/scrna_integration/bin/activate
which python  # Should point to your env
```

### Out of memory during pip install

Use a compute node for installation:
```bash
salloc --account=def-dcook --time=1:00:00 --mem=16G --cpus-per-task=4
# Then install packages
```

### scvi-tools training hangs

Set PyTorch precision for H100:
```python
import torch
torch.set_float32_matmul_precision("high")  # Already in your scripts
```

### Slow I/O during training

Copy data to `$SLURM_TMPDIR` at job start:
```bash
# In SLURM script, before running Python:
cp /path/to/preprocessed.h5ad $SLURM_TMPDIR/
# Then modify script to read from $SLURM_TMPDIR
```

---

## Quick Setup Summary

```bash
# One-time setup (on login node)
module load StdEnv/2023 python/3.11
virtualenv --no-download ~/envs/scrna_integration
source ~/envs/scrna_integration/bin/activate
pip install --no-index --upgrade pip
pip install --no-index torch torchvision torchaudio numpy scipy pandas scikit-learn matplotlib seaborn h5py
pip install scvi-tools scanpy anndata

# For each session
module load StdEnv/2023 python/3.11 cuda/12.9
source ~/envs/scrna_integration/bin/activate
```

All scripts (scVI, scANVI, MrVI, CellAssign) use the same `scrna_integration` environment.

---

## Notes for Nibi Specifically

- Nibi has H100 GPUs (80GB HBM3 memory)
- Use `--gres=gpu:h100:1` for GPU requests
- Nibi is the successor to Niagara for SciNet/Toronto users
- CUDA version on Nibi is 12.9 (`module load cuda/12.9`)
- Login nodes do NOT have GPUs — always test on compute nodes
- All pipeline scripts use the single `scrna_integration` environment
