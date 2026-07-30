#!/usr/bin/env bash
#SBATCH --job-name=hla_norm_depth
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --time=03:00:00
#SBATCH --array=1-40%12
#SBATCH --output=/data2/wangxuedong/polyploid_benchmark_v2/PolyploidHLA_simulation/logs/slurm/class_i_norm_depth_%A_%a.out
#SBATCH --error=/data2/wangxuedong/polyploid_benchmark_v2/PolyploidHLA_simulation/logs/slurm/class_i_norm_depth_%A_%a.err

set -euo pipefail

PROJECT_ROOT=/data2/wangxuedong/polyploid_benchmark_v2/PolyploidHLA
BENCH_ROOT=/data2/wangxuedong/polyploid_benchmark_v2/PolyploidHLA_simulation
PYTHON=/data3/wangxuedong/app/miniconda3/envs/polyploid-hla/bin/python
IMGT="$PROJECT_ROOT/resources/spechla/db/ref/hla_gen.format.filter.extend.DRB.no26789.v2.fasta"
RESULT_ROOT="$BENCH_ROOT/results/class_i_normalized_v1"

task_id=${SLURM_ARRAY_TASK_ID:?This script must run as a Slurm array task}
if (( task_id <= 20 )); then
    depth=100
    shard_number=$task_id
    experiment_prefix=accuracy_depth_v1_cov0100
else
    depth=300
    shard_number=$((task_id - 20))
    experiment_prefix=accuracy_main_v4
fi

shard=$(printf '%02d' "$shard_number")
sample_start=$((1 + (shard_number - 1) * 4))
experiment="${experiment_prefix}_shard${shard}"
shard_root="$RESULT_ROOT/${depth}x/shards"
work_root="${TMPDIR:-/tmp}/polyploid_class_i_norm_${SLURM_ARRAY_JOB_ID}_${depth}x_${shard}"
proposal="$shard_root/source_agnostic_class1_shard${shard}.tsv"
hybrid="$shard_root/class1_normalized_v1_shard${shard}.tsv"

mkdir -p "$shard_root"
rm -rf "$work_root"
trap 'rm -rf "$work_root"' EXIT

echo "[$(date '+%Y-%m-%d %H:%M:%S')] START depth=${depth}x shard=$shard samples=$sample_start-$((sample_start + 3))"

if [[ -s "$hybrid" ]] && [[ $(wc -l < "$hybrid") -eq 37 ]]; then
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] REUSE completed hybrid=$hybrid"
    exit 0
fi

if [[ -s "$proposal" ]] && [[ $(wc -l < "$proposal") -eq 37 ]]; then
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] REUSE proposal=$proposal"
else
    "$PYTHON" "$PROJECT_ROOT/diagnostics/offline_joint_quartet_posterior.py" \
        --bench-root "$BENCH_ROOT" \
        --experiment-glob "$experiment" \
        --genes HLA-A HLA-B HLA-C \
        --out "$proposal" \
        --summary "$shard_root/source_agnostic_class1_shard${shard}.summary.tsv"
fi

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
    --coverages "$depth" \
    --seed 20260728

"$PYTHON" "$PROJECT_ROOT/diagnostics/offline_class_i_hybrid_quartet.py" \
    --proposal-tsv "$proposal" \
    --fastq-bench-root "$work_root" \
    --conditions "graft10_cov${depth}x" "graft20_cov${depth}x" "graft30_cov${depth}x" \
    --min-log-bayes-factor=-inf \
    --min-log-bayes-factor-per-informative-pair 0.10 \
    --min-proposal-private-pairs 10 \
    --min-private-pair-ratio 5.0 \
    --out "$hybrid" \
    --summary "$shard_root/class1_normalized_v1_shard${shard}.summary.tsv"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] DONE depth=${depth}x shard=$shard"