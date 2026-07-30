#!/usr/bin/env bash
#SBATCH --job-name=hla_50x_hybrid
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --array=1-20%12
#SBATCH --output=/data2/wangxuedong/polyploid_benchmark_v2/PolyploidHLA_simulation/logs/slurm/hybrid_50x_%A_%a.out
#SBATCH --error=/data2/wangxuedong/polyploid_benchmark_v2/PolyploidHLA_simulation/logs/slurm/hybrid_50x_%A_%a.err

set -euo pipefail

PROJECT_ROOT=/data2/wangxuedong/polyploid_benchmark_v2/PolyploidHLA
BENCH_ROOT=/data2/wangxuedong/polyploid_benchmark_v2/PolyploidHLA_simulation
PYTHON=/data3/wangxuedong/app/miniconda3/envs/polyploid-hla/bin/python
IMGT="$PROJECT_ROOT/resources/spechla/db/ref/hla_gen.format.filter.extend.DRB.no26789.v2.fasta"
SHARD_ROOT="$BENCH_ROOT/results/hybrid_50x_bf180/shards"

task_id=${SLURM_ARRAY_TASK_ID:?This script must run as a Slurm array task}
shard=$(printf '%02d' "$task_id")
sample_start=$((1 + (task_id - 1) * 4))
experiment="accuracy_depth_v1_cov0050_shard${shard}"
work_root="${TMPDIR:-/tmp}/polyploid_hybrid_50x_${SLURM_ARRAY_JOB_ID}_${shard}"
proposal="$SHARD_ROOT/source_agnostic_class1_shard${shard}.tsv"
hybrid="$SHARD_ROOT/class1_hybrid_bf180_shard${shard}.tsv"

mkdir -p "$SHARD_ROOT"
rm -rf "$work_root"
trap 'rm -rf "$work_root"' EXIT

echo "[$(date '+%Y-%m-%d %H:%M:%S')] START shard=$shard samples=$sample_start-$((sample_start + 3))"

"$PYTHON" "$PROJECT_ROOT/diagnostics/offline_joint_quartet_posterior.py" \
    --bench-root "$BENCH_ROOT" \
    --experiment-glob "$experiment" \
    --genes HLA-A HLA-B HLA-C \
    --out "$proposal" \
    --summary "$SHARD_ROOT/source_agnostic_class1_shard${shard}.summary.tsv"

"$PYTHON" "$PROJECT_ROOT/benchmark/scripts/simulate_wgsim_benchmark.py" \
    --project-root "$PROJECT_ROOT" \
    --bench-root "$work_root" \
    --imgt-fasta "$IMGT" \
    --caller-alleles "$IMGT" \
    --experiment "$experiment" \
    --scenario distinct4 \
    --sample-start "$sample_start" \
    --individuals 4 \
    --graft-fractions 0.10 0.20 0.30 \
    --coverages 50 \
    --seed 20260728

"$PYTHON" "$PROJECT_ROOT/diagnostics/offline_class_i_hybrid_quartet.py" \
    --proposal-tsv "$proposal" \
    --fastq-bench-root "$work_root" \
    --min-log-bayes-factor 180 \
    --out "$hybrid" \
    --summary "$SHARD_ROOT/class1_hybrid_bf180_shard${shard}.summary.tsv"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] DONE shard=$shard"