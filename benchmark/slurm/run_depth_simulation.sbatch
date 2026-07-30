#!/usr/bin/env bash
#SBATCH --job-name=hla_depth_v1
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=3-00:00:00
#SBATCH --array=1-60%12
#SBATCH --output=/data2/wangxuedong/polyploid_benchmark_v2/PolyploidHLA_simulation/logs/slurm/depth_v1_%A_%a.out
#SBATCH --error=/data2/wangxuedong/polyploid_benchmark_v2/PolyploidHLA_simulation/logs/slurm/depth_v1_%A_%a.err

set -euo pipefail

PROJECT_ROOT=/data2/wangxuedong/polyploid_benchmark_v2/PolyploidHLA
BENCH_ROOT=/data2/wangxuedong/polyploid_benchmark_v2/PolyploidHLA_simulation
CONDA_PREFIX=/data3/wangxuedong/app/miniconda3/envs/polyploid-hla
export PATH="$CONDA_PREFIX/bin:$PATH"
PYTHON="$CONDA_PREFIX/bin/python"
export THREADS=${SLURM_CPUS_PER_TASK:-4}
export SAMTOOLS_THREADS=$THREADS

task_id=${SLURM_ARRAY_TASK_ID:?This script must run as a Slurm array task}
if (( task_id < 1 || task_id > 60 )); then
    echo "Invalid array task: $task_id" >&2
    exit 2
fi

depths=(50 100 1000)
depth_index=$(((task_id - 1) / 20))
shard=$(((task_id - 1) % 20 + 1))
coverage=${depths[$depth_index]}
sample_start=$((1 + (shard - 1) * 4))
individuals=4
CALLER_ALLELES="$PROJECT_ROOT/resources/spechla/db/ref/hla_gen.format.filter.extend.DRB.no26789.v2.fasta"
experiment=$(printf "accuracy_depth_v1_cov%04d_shard%02d" "$coverage" "$shard")

echo "[$(date '+%Y-%m-%d %H:%M:%S')] START job=$SLURM_JOB_ID task=$task_id depth=$coverage shard=$shard samples=$sample_start-$((sample_start + individuals - 1))"

"$PYTHON" -u "$PROJECT_ROOT/benchmark/scripts/simulate_wgsim_benchmark.py" \
    --project-root "$PROJECT_ROOT" \
    --bench-root "$BENCH_ROOT" \
    --imgt-fasta "$CALLER_ALLELES" \
    --caller-alleles "$CALLER_ALLELES" \
    --experiment "$experiment" \
    --scenario distinct4 \
    --sample-start "$sample_start" \
    --individuals "$individuals" \
    --graft-fractions 0.10 0.20 0.30 \
    --coverages "$coverage" \
    --seed 20260728 \
    --batch-results "$BENCH_ROOT/results/accuracy_depth_v1.all_samples.tsv" \
    --cleanup-caller-intermediates \
    --run-caller

echo "[$(date '+%Y-%m-%d %H:%M:%S')] DONE job=$SLURM_JOB_ID task=$task_id depth=$coverage shard=$shard"