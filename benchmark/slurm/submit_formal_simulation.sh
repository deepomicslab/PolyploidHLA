#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
BENCH_ROOT=${BENCH_ROOT:-"$(cd "$PROJECT_ROOT/.." && pwd)/PolyploidHLA_simulation"}

mkdir -p "$BENCH_ROOT/logs/slurm"
sbatch --array=1-20%20 "$SCRIPT_DIR/run_formal_simulation.sbatch"