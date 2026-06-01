#!/bin/bash
set -euo pipefail

# Read-only diagnostic: starting from remaining missing alleles after the
# guarded combined manifest, scan raw GenDx FASTQs and current gene-bin FASTQs
# for allele-specific exon k-mer support.

ROOT=${ROOT:-/data6/wangxuedong/polyploid_hla}
FQ_ROOT=${FQ_ROOT:-${ROOT}/fqs/amp/GenDx}
SPECHLA_ROOT=${SPECHLA_ROOT:-/data2/wangxuedong/polyploid-hla-realsets/spechla_out_gendx_amp_abc_20260520}
QUARTET_SUMMARY=${QUARTET_SUMMARY:-${ROOT}/diagnostics/gendx_quartet_summary_noA4.tsv}
COMBINED_MANIFEST=${COMBINED_MANIFEST:-${ROOT}/diagnostics/gendx_combined_rescue_noA4.manifest.tsv}
OUT=${OUT:-${ROOT}/diagnostics/gendx_remaining_raw_support_after_combined_noA4.tsv}
SUMMARY=${SUMMARY:-${ROOT}/diagnostics/gendx_remaining_raw_support_after_combined_noA4.summary.tsv}
PYBIN=${PYBIN:-python}

cd "$ROOT"
"$PYBIN" -m py_compile scripts/diagnostics/diagnose_gendx_remaining_raw_support.py
"$PYBIN" scripts/diagnostics/diagnose_gendx_remaining_raw_support.py \
    --quartet-summary "$QUARTET_SUMMARY" \
    --combined-manifest "$COMBINED_MANIFEST" \
    --fq-root "$FQ_ROOT" \
    --spechla-root "$SPECHLA_ROOT" \
    --out "$OUT" \
    --summary "$SUMMARY" \
    --k 31 \
    --k 51

printf 'summary\t%s\n' "$SUMMARY"
printf 'detail\t%s\n' "$OUT"
