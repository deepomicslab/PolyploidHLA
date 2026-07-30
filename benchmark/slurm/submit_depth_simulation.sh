#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
BENCH_ROOT=${BENCH_ROOT:-"$(cd "$PROJECT_ROOT/.." && pwd)/PolyploidHLA_simulation"}

mkdir -p "$BENCH_ROOT/logs/slurm" "$BENCH_ROOT/results" "$BENCH_ROOT/metrics"
array_job=$(sbatch --parsable --array=1-60%12 "$SCRIPT_DIR/run_depth_simulation.sbatch")
score_job=$(sbatch --parsable --dependency="afterok:${array_job}" "$SCRIPT_DIR/score_depth_simulation.sbatch")

printf 'depth_array_job=%s\n' "$array_job"
printf 'depth_score_job=%s\n' "$score_job"
printf 'batch_result=%s\n' "$BENCH_ROOT/results/accuracy_depth_v1.all_samples.tsv"
printf 'score_summary=%s\n' "$BENCH_ROOT/metrics/accuracy_depth_v1/allele_set_summary.tsv"