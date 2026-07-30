#!/usr/bin/env bash
#SBATCH --job-name=hla_shared1_joint
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=1-00:00:00
#SBATCH --output=/data2/wangxuedong/polyploid_benchmark_v2/PolyploidHLA_simulation/logs/slurm/shared1_joint_50x_%j.out
#SBATCH --error=/data2/wangxuedong/polyploid_benchmark_v2/PolyploidHLA_simulation/logs/slurm/shared1_joint_50x_%j.err

set -euo pipefail

PROJECT_ROOT=/data2/wangxuedong/polyploid_benchmark_v2/PolyploidHLA
BENCH_ROOT=/data2/wangxuedong/polyploid_benchmark_v2/PolyploidHLA_simulation
CONDA_PREFIX=/data3/wangxuedong/app/miniconda3/envs/polyploid-hla
export PATH="$CONDA_PREFIX/bin:$PATH"
PYTHON="$CONDA_PREFIX/bin/python"
export THREADS=${SLURM_CPUS_PER_TASK:-4}
export SAMTOOLS_THREADS=$THREADS

EXPERIMENT=shared1_joint_50x_v1
CALLER_ALLELES="$PROJECT_ROOT/resources/spechla/db/ref/hla_gen.format.filter.extend.DRB.no26789.v2.fasta"
METRICS_ROOT="$BENCH_ROOT/metrics/$EXPERIMENT"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] START job=$SLURM_JOB_ID experiment=$EXPERIMENT"

"$PYTHON" -u "$PROJECT_ROOT/benchmark/scripts/simulate_wgsim_benchmark.py" \
    --project-root "$PROJECT_ROOT" \
    --bench-root "$BENCH_ROOT" \
    --imgt-fasta "$CALLER_ALLELES" \
    --caller-alleles "$CALLER_ALLELES" \
    --experiment "$EXPERIMENT" \
    --scenario shared1 \
    --individuals 4 \
    --genes DRB1 DQB1 DPB1 \
    --graft-fractions 0.10 0.20 0.30 \
    --coverages 50 \
    --seed 20260728 \
    --batch-results "$BENCH_ROOT/results/$EXPERIMENT.all_samples.tsv" \
    --cleanup-caller-intermediates \
    --run-caller

"$PYTHON" "$PROJECT_ROOT/benchmark/scripts/score_simulation_allele_sets.py" \
    --bench-root "$BENCH_ROOT" \
    --experiment-glob "$EXPERIMENT" \
    --genes HLA-DRB1 HLA-DQB1 HLA-DPB1 \
    --out-dir "$METRICS_ROOT/baseline"

"$PYTHON" "$PROJECT_ROOT/diagnostics/offline_joint_quartet_posterior.py" \
    --bench-root "$BENCH_ROOT" \
    --experiment-glob "$EXPERIMENT" \
    --genes HLA-DRB1 HLA-DQB1 HLA-DPB1 \
    --out "$METRICS_ROOT/joint_quartet.tsv" \
    --summary "$METRICS_ROOT/joint_quartet.summary.tsv"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] DONE job=$SLURM_JOB_ID experiment=$EXPERIMENT"