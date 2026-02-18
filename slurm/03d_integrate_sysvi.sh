#!/bin/bash
#SBATCH --job-name=03d_sysvi
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --gres=gpu:h100:1
#SBATCH -o /project/rrg-tperkins/OBCF/active/BHI_single_cell_processing/logs/slurm/03d_integrate_sysvi.sh.%A_%a.out

# Step 03d: SysVI integration (GPU)
# Requires: preprocessed.h5ad from Step 01

module load StdEnv/2023 python/3.11 cuda/12.9

source /project/rrg-tperkins/OBCF/active/BHI_single_cell_processing/envs/scrna_integration/bin/activate

echo "=========================================="
echo "Job: $SLURM_JOB_NAME"
echo "Job ID: $SLURM_JOB_ID"
echo "Started: $(date)"
echo "=========================================="

echo ""
echo "GPU Status:"
nvidia-smi --query-gpu=name,memory.total,memory.free --format=csv
echo ""

python /project/rrg-tperkins/OBCF/active/BHI_single_cell_processing/bhi-scrna-integration-pipeline/scripts/03d_integrate_sysvi.py

echo ""
echo "Completed: $(date)"
