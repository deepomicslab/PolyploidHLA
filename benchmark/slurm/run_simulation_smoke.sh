#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
BENCH_ROOT=${BENCH_ROOT:-"$(cd "$PROJECT_ROOT/.." && pwd)/PolyploidHLA_simulation"}
CONDA_ENV=${CONDA_ENV:-polyploid-hla}
RUN_CALLER=${RUN_CALLER:-1}

if [[ "${CONDA_DEFAULT_ENV:-}" != "$CONDA_ENV" ]]; then
    ENV_PREFIX=$(conda run -n "$CONDA_ENV" python -c 'import sys; print(sys.prefix)')
    exec conda run --no-capture-output -n "$CONDA_ENV" \
        env PATH="$ENV_PREFIX/bin:$PATH" bash "$0" "$@"
fi

export PATH="$CONDA_PREFIX/bin:$PATH"

args=(
    "$PROJECT_ROOT/benchmark/scripts/simulate_wgsim_benchmark.py"
    --bench-root "$BENCH_ROOT"
    --experiment smoke_closedset
    --scenario distinct4
    --individuals 1
    --graft-fractions 0.10
    --coverages 50
    --overwrite
)

if [[ "$RUN_CALLER" == "1" ]]; then
    args+=(--run-caller)
fi

python "${args[@]}"