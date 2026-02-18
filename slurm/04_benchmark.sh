#!/bin/bash
#SBATCH --job-name=04_benchmark
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH -o /project/rrg-tperkins/OBCF/active/BHI_single_cell_processing/logs/slurm/04_benchmark.sh.%A_%a.out

# Step 04: Benchmark integration methods using scib-metrics
# 
# Requires:
#   - Preprocessed data (from Step 01)
#   - CellAssign predictions (from Step 02)
#   - Embeddings: harmony, scvi, scanvi, mrvi (from Step 03)
#
# Dependencies: scib-metrics, scanpy

module load StdEnv/2023 python/3.11

source /project/rrg-tperkins/OBCF/active/BHI_single_cell_processing/envs/scrna_integration/bin/activate

echo "=========================================="
echo "Job: $SLURM_JOB_NAME"
echo "Job ID: $SLURM_JOB_ID"
echo "Started: $(date)"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "=========================================="

python /project/rrg-tperkins/OBCF/active/BHI_single_cell_processing/bhi-scrna-integration-pipeline/scripts/04_benchmark.py

echo ""
echo "Completed: $(date)"
