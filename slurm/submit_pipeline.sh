#!/bin/bash
# submit_pipeline.sh - Submit the full scRNA-seq integration pipeline
#
# Usage:
#   ./submit_pipeline.sh [options]
#
# Options:
#   --from STEP     Start from a specific step (1-5, default: 1)
#   --to STEP       Stop at a specific step (1-5, default: 5)
#   --method METHOD Method for finalize step (default: scanvi)
#   --dry-run       Show what would be submitted without actually submitting
#   --help          Show this help message
#
# Examples:
#   ./submit_pipeline.sh                    # Run full pipeline
#   ./submit_pipeline.sh --from 3           # Start from integration (skip 01, 02)
#   ./submit_pipeline.sh --from 3 --to 4    # Only run integration + benchmark
#   ./submit_pipeline.sh --dry-run          # Preview submission plan

set -e

# Defaults
FROM_STEP=1
TO_STEP=5
METHOD="scanvi"
DRY_RUN=false

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --from)
            FROM_STEP="$2"
            shift 2
            ;;
        --to)
            TO_STEP="$2"
            shift 2
            ;;
        --method)
            METHOD="$2"
            shift 2
            ;;
        --dry-run)
            DRY_RUN=true
            shift
            ;;
        --help)
            head -20 "$0" | tail -18
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Create logs directory
# mkdir -p logs

echo "=================================================="
echo "scRNA-seq Integration Pipeline Submission"
echo "=================================================="
echo "Steps: $FROM_STEP to $TO_STEP"
echo "Finalize method: $METHOD"
echo "Dry run: $DRY_RUN"
echo ""

# Function to submit job
submit_job() {
    local script=$1
    local deps=$2
    local extra_args=$3
    
    if [[ -n "$deps" ]]; then
        dep_flag="--dependency=afterok:$deps"
    else
        dep_flag=""
    fi
    
    if [[ "$DRY_RUN" == true ]]; then
        echo "[DRY RUN] sbatch $dep_flag $script $extra_args"
        echo "FAKE_JOB_ID_$$"
    else
        result=$(sbatch $dep_flag $script $extra_args)
        job_id=$(echo $result | awk '{print $4}')
        echo "$result"
        echo "$job_id"
    fi
}

# Track job IDs for dependencies
JOB_01=""
JOB_02=""
JOB_03A=""
JOB_03B=""
JOB_03C=""
JOB_03D=""
JOB_04=""

# ==========================================
# Step 01: Preprocessing
# ==========================================
if [[ $FROM_STEP -le 1 && $TO_STEP -ge 1 ]]; then
    echo "Submitting Step 01: Preprocessing..."
    JOB_01=$(submit_job "01_preprocess.sh" "" "" | tail -1)
    echo "  Job ID: $JOB_01"
    echo ""
fi

# ==========================================
# Step 02: CellAssign
# ==========================================
if [[ $FROM_STEP -le 2 && $TO_STEP -ge 2 ]]; then
    echo "Submitting Step 02: CellAssign..."
    # CellAssign uses raw.h5ad, so no dependency on Step 01
    JOB_02=$(submit_job "02_cellassign.sh" "" "" | tail -1)
    echo "  Job ID: $JOB_02"
    echo ""
fi

# ==========================================
# Step 03: Integration (parallel except scANVI)
# ==========================================
if [[ $FROM_STEP -le 3 && $TO_STEP -ge 3 ]]; then
    # Build dependency string for Step 03
    STEP3_DEPS=""
    if [[ -n "$JOB_01" ]]; then
        STEP3_DEPS="$JOB_01"
    fi
    if [[ -n "$JOB_02" ]]; then
        if [[ -n "$STEP3_DEPS" ]]; then
            STEP3_DEPS="$STEP3_DEPS:$JOB_02"
        else
            STEP3_DEPS="$JOB_02"
        fi
    fi
    
    echo "Submitting Step 03: Integration methods..."
    echo "  Dependencies: ${STEP3_DEPS:-none}"
    
    # 03a: scVI (can run in parallel with 03c, 03d)
    echo "  03a: scVI..."
    JOB_03A=$(submit_job "03a_integrate_scvi.sh" "$STEP3_DEPS" "" | tail -1)
    echo "    Job ID: $JOB_03A"
    
    # 03c: MrVI (parallel with scVI)
    echo "  03c: MrVI..."
    JOB_03C=$(submit_job "03c_integrate_mrvi.sh" "$STEP3_DEPS" "" | tail -1)
    echo "    Job ID: $JOB_03C"
    
    # 03d: SysVI (parallel with scVI)
    echo "  03d: SysVI..."
    JOB_03D=$(submit_job "03d_integrate_sysvi.sh" "$STEP3_DEPS" "" | tail -1)
    echo "    Job ID: $JOB_03D"
    
    # 03b: scANVI (must wait for scVI to complete)
    SCANVI_DEPS="$JOB_03A"
    if [[ -n "$JOB_02" ]]; then
        SCANVI_DEPS="$SCANVI_DEPS:$JOB_02"
    fi
    echo "  03b: scANVI (depends on scVI)..."
    JOB_03B=$(submit_job "03b_integrate_scanvi.sh" "$SCANVI_DEPS" "" | tail -1)
    echo "    Job ID: $JOB_03B"
    
    echo ""
fi

# ==========================================
# Step 04: Benchmark
# ==========================================
if [[ $FROM_STEP -le 4 && $TO_STEP -ge 4 ]]; then
    # Depends on all integration jobs completing
    STEP4_DEPS=""
    for job in $JOB_03A $JOB_03B $JOB_03C $JOB_03D; do
        if [[ -n "$job" ]]; then
            if [[ -n "$STEP4_DEPS" ]]; then
                STEP4_DEPS="$STEP4_DEPS:$job"
            else
                STEP4_DEPS="$job"
            fi
        fi
    done
    
    echo "Submitting Step 04: Benchmark..."
    echo "  Dependencies: ${STEP4_DEPS:-none}"
    JOB_04=$(submit_job "04_benchmark.sh" "$STEP4_DEPS" "" | tail -1)
    echo "  Job ID: $JOB_04"
    echo ""
fi

# ==========================================
# Step 05: Finalize
# ==========================================
if [[ $FROM_STEP -le 5 && $TO_STEP -ge 5 ]]; then
    # Depends on benchmark (or specific integration job if running alone)
    STEP5_DEPS=""
    if [[ -n "$JOB_04" ]]; then
        STEP5_DEPS="$JOB_04"
    else
        # If not running benchmark, depend on relevant integration job
        case $METHOD in
            scvi)   STEP5_DEPS="$JOB_03A" ;;
            scanvi) STEP5_DEPS="$JOB_03B" ;;
            mrvi)   STEP5_DEPS="$JOB_03C" ;;
            sysvi)  STEP5_DEPS="$JOB_03D" ;;
        esac
    fi
    
    echo "Submitting Step 05: Finalize (method=$METHOD)..."
    echo "  Dependencies: ${STEP5_DEPS:-none}"
    JOB_05=$(submit_job "05_finalize.sh" "$STEP5_DEPS" "$METHOD" | tail -1)
    echo "  Job ID: $JOB_05"
    echo ""
fi

# ==========================================
# Summary
# ==========================================
echo "=================================================="
echo "Pipeline Submitted!"
echo "=================================================="
echo ""
echo "Job Summary:"
[[ -n "$JOB_01" ]] && echo "  01_preprocess: $JOB_01"
[[ -n "$JOB_02" ]] && echo "  02_cellassign: $JOB_02"
[[ -n "$JOB_03A" ]] && echo "  03a_scvi:      $JOB_03A"
[[ -n "$JOB_03B" ]] && echo "  03b_scanvi:    $JOB_03B (waits for 03a)"
[[ -n "$JOB_03C" ]] && echo "  03c_mrvi:      $JOB_03C"
[[ -n "$JOB_03D" ]] && echo "  03d_sysvi:     $JOB_03D"
[[ -n "$JOB_04" ]] && echo "  04_benchmark:  $JOB_04"
[[ -n "$JOB_05" ]] && echo "  05_finalize:   $JOB_05"
echo ""
echo "Monitor with: squeue -u \$USER"
echo "Cancel all:   scancel -u \$USER"
echo "View logs:    tail -f logs/<jobname>_<jobid>.out"
