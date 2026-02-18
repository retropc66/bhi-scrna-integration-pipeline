#!/bin/bash
#SBATCH --job-name=01_preprocess
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=4:00:00
#SBATCH -o /project/rrg-tperkins/OBCF/active/BHI_single_cell_processing/logs/slurm/01_preprocess.sh.%A_%a.out

# Step 01: Preprocessing (CPU)
# Identifies HVGs and creates preprocessed.h5ad

module load StdEnv/2023 python/3.11

source /project/rrg-tperkins/OBCF/active/BHI_single_cell_processing/envs/scrna_integration/bin/activate

echo "=========================================="
echo "Job: $SLURM_JOB_NAME"
echo "Job ID: $SLURM_JOB_ID"
echo "Started: $(date)"
echo "=========================================="

python /project/rrg-tperkins/OBCF/active/BHI_single_cell_processing/bhi-scrna-integration-pipeline/scripts/01_preprocess.py

echo ""
echo "Completed: $(date)"
