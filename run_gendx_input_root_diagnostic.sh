#!/bin/bash
set -euo pipefail

# Truth-free input-layer diagnostic for GenDx amplicon samples.
# It measures how many gene-unique read pairs are present in the full FASTQs,
# how many were retained by strict SpecHLA gene binning, and how many would be
# eligible for conservative read-bin rescue.

FQ_ROOT=${FQ_ROOT:-/data6/wangxuedong/polyploid_hla/fqs/amp/GenDx}
SPECHLA_ROOT=${SPECHLA_ROOT:-/data2/wangxuedong/polyploid-hla-realsets/spechla_out_gendx_amp_abc_20260520}
OUT_DIR=${OUT_DIR:-/data6/wangxuedong/polyploid_hla/diagnostics/gendx_input_root_diagnostic_noA4}
PYBIN=${PYBIN:-python}
SCRIPTS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EXCLUDE_SAMPLES=${EXCLUDE_SAMPLES:-}
GENES=${GENES:-"HLA-A HLA-B HLA-C HLA-DRB1 HLA-DPB1 HLA-DQB1"}

mkdir -p "$OUT_DIR"

is_excluded() {
    local sample="$1"
    for excluded in $EXCLUDE_SAMPLES; do
        [[ "$sample" == "$excluded" ]] && return 0
    done
    return 1
}

for fq1 in "$FQ_ROOT"/*_R1_001.fastq.gz; do
    sample=$(basename "$fq1" _R1_001.fastq.gz)
    is_excluded "$sample" && continue
    fq2="$FQ_ROOT/${sample}_R2_001.fastq.gz"
    fq_dir="$SPECHLA_ROOT/$sample"
    manifest="$OUT_DIR/${sample}.read_bin_rescue.tsv"
    log="$OUT_DIR/${sample}.log"
    gene_args=()
    for gene in $GENES; do
        gene_args+=(--gene "$gene")
    done
    "$PYBIN" "$SCRIPTS_DIR/rescue_gene_binned_reads.py" \
        --fq1 "$fq1" \
        --fq2 "$fq2" \
        --fq-dir "$fq_dir" \
        --require-both-mates \
        --retention-gate \
        --dry-run \
        --manifest "$manifest" \
        "${gene_args[@]}" > "$log" 2>&1
    echo "$sample"
done

summary="$OUT_DIR/summary.tsv"
{
    printf 'sample\tgene\toriginal_pairs\tfull_gene_pairs\tretained_gene_pairs\tretained_fraction\tmissing_fraction\trescued_pairs\tstatus\tretention_reason\n'
    for manifest in "$OUT_DIR"/*.read_bin_rescue.tsv; do
        sample=$(basename "$manifest" .read_bin_rescue.tsv)
        awk -F '\t' -v sample="$sample" 'NR > 1 && $1 !~ /^#/ {print sample"\t"$1"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$14"\t"$13}' "$manifest"
    done
} > "$summary"

echo "wrote $summary"
