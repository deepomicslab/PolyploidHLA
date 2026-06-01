#!/bin/bash
set -euo pipefail

ROOT=${ROOT:-/data6/wangxuedong/polyploid_hla}
REALSETS=${REALSETS:-/data2/wangxuedong/polyploid-hla-realsets}

export GATK_CALLER=${GATK_CALLER:-Mutect2}
export CALLER_TAG=${CALLER_TAG:-mutect2}
export SCORE_PREFIX=${SCORE_PREFIX:-mutect2_chi_count_strict}
export GT_DROP_UNASSIGNABLE=${GT_DROP_UNASSIGNABLE:-1}
export VALID_SPECHLA_ROOT=${VALID_SPECHLA_ROOT:-${REALSETS}/spechla_out_gendx_amp_abc_mutect2_chi_count_strict_20260526}
export VALID_ASM_ROOT=${VALID_ASM_ROOT:-${REALSETS}/asm_v2_gendx_amp_abc_mutect2_chi_count_strict_20260526}
export DIAG_DIR=${DIAG_DIR:-${ROOT}/diagnostics/gendx_mutect2_chi_count_strict_20260526}

exec bash "${ROOT}/scripts/run_gendx_gatk_ploidy4_validation.sh"
