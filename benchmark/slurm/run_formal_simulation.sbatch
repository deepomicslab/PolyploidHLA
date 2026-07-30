#!/usr/bin/env bash
#SBATCH --job-name=hla_sim_v4
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=1-12:00:00
#SBATCH --array=1-20%20
#SBATCH --output=/data2/wangxuedong/polyploid_benchmark_v2/PolyploidHLA_simulation/logs/slurm/formal_main_v4_%A_%a.out
#SBATCH --error=/data2/wangxuedong/polyploid_benchmark_v2/PolyploidHLA_simulation/logs/slurm/formal_main_v4_%A_%a.err

set -euo pipefail

PROJECT_ROOT=/data2/wangxuedong/polyploid_benchmark_v2/PolyploidHLA
BENCH_ROOT=/data2/wangxuedong/polyploid_benchmark_v2/PolyploidHLA_simulation
CONDA_PREFIX=/data3/wangxuedong/app/miniconda3/envs/polyploid-hla
export PATH="$CONDA_PREFIX/bin:$PATH"
PYTHON="$CONDA_PREFIX/bin/python"
export THREADS=${SLURM_CPUS_PER_TASK:-8}
export SAMTOOLS_THREADS=$THREADS

task_id=${SLURM_ARRAY_TASK_ID:?This script must run as a Slurm array task}
if (( task_id < 1 || task_id > 20 )); then
    echo "Invalid array task: $task_id" >&2
    exit 2
fi
sample_start=$((1 + (task_id - 1) * 4))
individuals=4
CALLER_ALLELES="$PROJECT_ROOT/resources/spechla/db/ref/hla_gen.format.filter.extend.DRB.no26789.v2.fasta"

experiment=$(printf "accuracy_main_v4_shard%02d" "$task_id")
echo "[$(date '+%Y-%m-%d %H:%M:%S')] START job=$SLURM_JOB_ID task=$task_id samples=$sample_start-$((sample_start + individuals - 1))"

"$PYTHON" -u "$PROJECT_ROOT/benchmark/scripts/simulate_wgsim_benchmark.py" \
    --project-root "$PROJECT_ROOT" \
    --bench-root "$BENCH_ROOT" \
    --imgt-fasta "$CALLER_ALLELES" \
    --caller-alleles "$CALLER_ALLELES" \
    --experiment "$experiment" \
    --scenario distinct4 \
    --sample-start "$sample_start" \
    --individuals "$individuals" \
    --graft-fractions 0.10 0.20 0.30 0.40 0.50 \
    --coverages 300 \
    --seed 20260728 \
    --batch-results "$BENCH_ROOT/results/accuracy_main_v4.all_samples.tsv" \
    --run-caller

echo "[$(date '+%Y-%m-%d %H:%M:%S')] DONE job=$SLURM_JOB_ID task=$task_id"