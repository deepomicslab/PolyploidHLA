#!/usr/bin/env bash
#SBATCH --job-name=hla_depth_score
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=/data2/wangxuedong/polyploid_benchmark_v2/PolyploidHLA_simulation/logs/slurm/depth_score_%j.out
#SBATCH --error=/data2/wangxuedong/polyploid_benchmark_v2/PolyploidHLA_simulation/logs/slurm/depth_score_%j.err

set -euo pipefail

PROJECT_ROOT=/data2/wangxuedong/polyploid_benchmark_v2/PolyploidHLA
BENCH_ROOT=/data2/wangxuedong/polyploid_benchmark_v2/PolyploidHLA_simulation
PYTHON=/data3/wangxuedong/app/miniconda3/envs/polyploid-hla/bin/python

"$PYTHON" "$PROJECT_ROOT/benchmark/scripts/score_simulation_allele_sets.py" \
    --bench-root "$BENCH_ROOT" \
    --experiment-glob 'accuracy_depth_v1_cov*_shard*' \
    --out-dir "$BENCH_ROOT/metrics/accuracy_depth_v1"