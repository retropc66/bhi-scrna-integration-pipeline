#!/bin/bash
#SBATCH --job-name=02_cellassign
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --gres=gpu:h100:1
#SBATCH -o /project/rrg-tperkins/OBCF/active/BHI_single_cell_processing/logs/slurm/02_cellassign.sh.%A_%a.out

# Step 02: CellAssign cell type annotation (GPU)
# Runs on raw.h5ad, produces predictions.csv

module load StdEnv/2023 python/3.11 cuda/12.9

source /project/rrg-tperkins/OBCF/active/BHI_single_cell_processing/envs/scrna_integration/bin/activate

echo "=========================================="
echo "Job: $SLURM_JOB_NAME"
echo "Job ID: $SLURM_JOB_ID"
echo "Started: $(date)"
echo "=========================================="

# Verify GPU
echo ""
echo "GPU Status:"
nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv
python -c "import torch; print(f'CUDA available: {torch.cuda.is_available()}')"
echo ""

python /project/rrg-tperkins/OBCF/active/BHI_single_cell_processing/bhi-scrna-integration-pipeline/scripts/02_cellassign.py

echo ""
echo "Completed: $(date)"
