#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
BENCH_ROOT=${BENCH_ROOT:-"$(cd "$PROJECT_ROOT/.." && pwd)/PolyploidHLA_simulation"}
mode=${1:-smoke}

mkdir -p "$BENCH_ROOT/logs/slurm"

case "$mode" in
    smoke)
        rescue_job=$(sbatch --parsable --array=1-3%3 --export=ALL,VALIDATION_MODE=smoke \
            "$SCRIPT_DIR/run_class_i_rescue_validation.sbatch")
        baseline_job=$(sbatch --parsable --dependency="afterok:$rescue_job" \
            --array=1-3%3 --export=ALL,VALIDATION_MODE=smoke \
            "$SCRIPT_DIR/run_class_i_rescue_baseline.sbatch")
        ;;
    formal)
        rescue_job=$(sbatch --parsable --array=1-60%20 --export=ALL,VALIDATION_MODE=formal \
            "$SCRIPT_DIR/run_class_i_rescue_validation.sbatch")
        baseline_job=$(sbatch --parsable --dependency="afterok:$rescue_job" \
            --array=1-60%20 --export=ALL,VALIDATION_MODE=formal \
            "$SCRIPT_DIR/run_class_i_rescue_baseline.sbatch")
        ;;
    *)
        echo "Usage: $0 [smoke|formal]" >&2
        exit 2
        ;;
esac

    echo "Submitted rescue job $rescue_job"
    echo "Submitted paired baseline job $baseline_job (afterok:$rescue_job)"