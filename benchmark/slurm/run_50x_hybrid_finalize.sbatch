#!/usr/bin/env bash
#SBATCH --job-name=hla_50x_summary
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:20:00
#SBATCH --output=/data2/wangxuedong/polyploid_benchmark_v2/PolyploidHLA_simulation/logs/slurm/hybrid_50x_summary_%j.out
#SBATCH --error=/data2/wangxuedong/polyploid_benchmark_v2/PolyploidHLA_simulation/logs/slurm/hybrid_50x_summary_%j.err

set -euo pipefail

PROJECT_ROOT=/data2/wangxuedong/polyploid_benchmark_v2/PolyploidHLA
BENCH_ROOT=/data2/wangxuedong/polyploid_benchmark_v2/PolyploidHLA_simulation
PYTHON=/data3/wangxuedong/app/miniconda3/envs/polyploid-hla/bin/python
RESULT_ROOT="$BENCH_ROOT/results/hybrid_50x_bf180"

"$PYTHON" "$PROJECT_ROOT/diagnostics/summarize_50x_hybrid_strategy.py" \
    --class-i-glob "$RESULT_ROOT/shards/class1_hybrid_bf180_shard??.tsv" \
    --class-ii "$BENCH_ROOT/metrics/joint_quartet_50x_class2_v2.tsv" \
    --log-bayes-factor-threshold 180 \
    --detail-output "$RESULT_ROOT/accuracy_depth_v1.hybrid_bf180_class2_joint_v2.per_locus.tsv" \
    --summary-output "$RESULT_ROOT/accuracy_depth_v1.hybrid_bf180_class2_joint_v2.summary.tsv"