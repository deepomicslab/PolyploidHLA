#!/usr/bin/env bash
#SBATCH --job-name=hla_ci_rescue
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=1-12:00:00
#SBATCH --output=/data2/wangxuedong/polyploid_benchmark_v2/PolyploidHLA_simulation/logs/slurm/class_i_rescue_%A_%a.out
#SBATCH --error=/data2/wangxuedong/polyploid_benchmark_v2/PolyploidHLA_simulation/logs/slurm/class_i_rescue_%A_%a.err

set -euo pipefail

PROJECT_ROOT=/data2/wangxuedong/polyploid_benchmark_v2/PolyploidHLA
BENCH_ROOT=/data2/wangxuedong/polyploid_benchmark_v2/PolyploidHLA_simulation
CONDA_PREFIX=/data3/wangxuedong/app/miniconda3/envs/polyploid-hla
export PATH="$CONDA_PREFIX/bin:$PATH"
PYTHON="$CONDA_PREFIX/bin/python"
export THREADS=${SLURM_CPUS_PER_TASK:-4}
export SAMTOOLS_THREADS=$THREADS

mode=${VALIDATION_MODE:-smoke}
task_id=${SLURM_ARRAY_TASK_ID:?This script must run as a Slurm array task}
scenarios=(distinct4 shared1 homozygous_graft)

case "$mode" in
    smoke)
        if (( task_id < 1 || task_id > 3 )); then
            echo "Smoke mode requires array task 1-3, got $task_id" >&2
            exit 2
        fi
        scenario=${scenarios[$((task_id - 1))]}
        sample_start=1
        individuals=1
        experiment="class_i_rescue_smoke_v1_${scenario}"
        fractions=(0.10)
        ;;
    formal)
        if (( task_id < 1 || task_id > 60 )); then
            echo "Formal mode requires array task 1-60, got $task_id" >&2
            exit 2
        fi
        scenario_index=$(((task_id - 1) / 20))
        shard=$(((task_id - 1) % 20 + 1))
        scenario=${scenarios[$scenario_index]}
        sample_start=$((1 + (shard - 1) * 4))
        individuals=4
        experiment=$(printf "class_i_rescue_validation_v1_%s_shard%02d" "$scenario" "$shard")
        fractions=(0.10 0.20 0.30)
        ;;
    *)
        echo "Unsupported VALIDATION_MODE=$mode (expected smoke or formal)" >&2
        exit 2
        ;;
esac

CALLER_ALLELES="$PROJECT_ROOT/resources/spechla/db/ref/hla_gen.format.filter.extend.DRB.no26789.v2.fasta"
export EM_REFINE_CLASS_I_DISTINCT_RESCUE=1
export CLASS2_DPB1_RARE_COLLAPSE=0

echo "[$(date '+%Y-%m-%d %H:%M:%S')] START mode=$mode job=$SLURM_JOB_ID task=$task_id scenario=$scenario samples=$sample_start-$((sample_start + individuals - 1))"

"$PYTHON" -u "$PROJECT_ROOT/benchmark/scripts/simulate_wgsim_benchmark.py" \
    --project-root "$PROJECT_ROOT" \
    --bench-root "$BENCH_ROOT" \
    --imgt-fasta "$CALLER_ALLELES" \
    --caller-alleles "$CALLER_ALLELES" \
    --experiment "$experiment" \
    --scenario "$scenario" \
    --sample-start "$sample_start" \
    --individuals "$individuals" \
    --graft-fractions "${fractions[@]}" \
    --coverages 300 \
    --seed 20260729 \
    --run-caller

echo "[$(date '+%Y-%m-%d %H:%M:%S')] DONE mode=$mode job=$SLURM_JOB_ID task=$task_id"