#!/usr/bin/env bash
#SBATCH --job-name=hla_norm_summary
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:20:00
#SBATCH --output=/data2/wangxuedong/polyploid_benchmark_v2/PolyploidHLA_simulation/logs/slurm/class_i_norm_summary_%j.out
#SBATCH --error=/data2/wangxuedong/polyploid_benchmark_v2/PolyploidHLA_simulation/logs/slurm/class_i_norm_summary_%j.err

set -euo pipefail

PROJECT_ROOT=/data2/wangxuedong/polyploid_benchmark_v2/PolyploidHLA
BENCH_ROOT=/data2/wangxuedong/polyploid_benchmark_v2/PolyploidHLA_simulation
PYTHON=/data3/wangxuedong/app/miniconda3/envs/polyploid-hla/bin/python
RESULT_ROOT="$BENCH_ROOT/results/class_i_normalized_v1"

"$PYTHON" "$PROJECT_ROOT/diagnostics/summarize_class_i_normalized_depths.py" \
    --input-glob "$BENCH_ROOT/results/hybrid_50x_bf180/shards/class1_hybrid_bf180_shard??.tsv" \
    --input-glob "$RESULT_ROOT/100x/shards/class1_normalized_v1_shard??.tsv" \
    --input-glob "$RESULT_ROOT/300x/shards/class1_normalized_v1_shard??.tsv" \
    --detail-output "$RESULT_ROOT/class_i_normalized_v1.per_locus.tsv" \
    --summary-output "$RESULT_ROOT/class_i_normalized_v1.summary.tsv"