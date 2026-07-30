#!/usr/bin/env bash
#SBATCH --job-name=hla_ci_base
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=1-12:00:00
#SBATCH --output=/data2/wangxuedong/polyploid_benchmark_v2/PolyploidHLA_simulation/logs/slurm/class_i_baseline_%A_%a.out
#SBATCH --error=/data2/wangxuedong/polyploid_benchmark_v2/PolyploidHLA_simulation/logs/slurm/class_i_baseline_%A_%a.err

set -euo pipefail

PROJECT_ROOT=/data2/wangxuedong/polyploid_benchmark_v2/PolyploidHLA
BENCH_ROOT=/data2/wangxuedong/polyploid_benchmark_v2/PolyploidHLA_simulation
CONDA_PREFIX=/data3/wangxuedong/app/miniconda3/envs/polyploid-hla
export PATH="$CONDA_PREFIX/bin:$PATH"
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
        experiment="class_i_rescue_smoke_v1_${scenario}"
        sample_start=1
        individuals=1
        conditions=(graft10_cov300x)
        ;;
    formal)
        if (( task_id < 1 || task_id > 60 )); then
            echo "Formal mode requires array task 1-60, got $task_id" >&2
            exit 2
        fi
        scenario_index=$(((task_id - 1) / 20))
        shard=$(((task_id - 1) % 20 + 1))
        scenario=${scenarios[$scenario_index]}
        experiment=$(printf "class_i_rescue_validation_v1_%s_shard%02d" "$scenario" "$shard")
        sample_start=$((1 + (shard - 1) * 4))
        individuals=4
        conditions=(graft10_cov300x graft20_cov300x graft30_cov300x)
        ;;
    *)
        echo "Unsupported VALIDATION_MODE=$mode (expected smoke or formal)" >&2
        exit 2
        ;;
esac

export EM_REFINE_CLASS_I_DISTINCT_RESCUE=0
export CLASS2_DPB1_RARE_COLLAPSE=0
baseline_experiment="${experiment}_baseline"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] START mode=$mode job=$SLURM_JOB_ID task=$task_id scenario=$scenario baseline=$baseline_experiment"

for condition in "${conditions[@]}"; do
    for ((offset = 0; offset < individuals; offset++)); do
        sample_id=$(printf "SIM%04d" "$((sample_start + offset))")
        fq1="$BENCH_ROOT/reads/$experiment/$condition/$sample_id.R1.fastq.gz"
        fq2="$BENCH_ROOT/reads/$experiment/$condition/$sample_id.R2.fastq.gz"
        if [[ ! -s "$fq1" || ! -s "$fq2" ]]; then
            echo "Missing paired validation FASTQs: $fq1 $fq2" >&2
            exit 2
        fi
        work_dir="$BENCH_ROOT/runs/$baseline_experiment/$condition/$sample_id"
        log_dir="$BENCH_ROOT/logs/$baseline_experiment/$condition"
        mkdir -p "$work_dir" "$log_dir"
        FQ1="$fq1" FQ2="$fq2" SAMPLE="$sample_id" RECIPIENT_MAJOR=1 \
            FB_MIN_AF=0.03 GT_DROP_FP_AF=0.05 WORK_DIR="$work_dir" \
            bash "$PROJECT_ROOT/polyphase_v2.sh" \
            >"$log_dir/$sample_id.polyphase.log" 2>&1
    done
done

echo "[$(date '+%Y-%m-%d %H:%M:%S')] DONE mode=$mode job=$SLURM_JOB_ID task=$task_id"