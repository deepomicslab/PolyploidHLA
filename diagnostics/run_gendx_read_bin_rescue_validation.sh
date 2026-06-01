#!/bin/bash
set -euo pipefail

# Independent GenDx read-bin rescue validation rerun.
#
# This does not modify the baseline SpecHLA/ASM roots and does not change the
# main pipeline. For each sample, it copies only the reusable pre-BWA baseline
# intermediates into a validation root, appends conservatively rescued gene-bin
# read pairs there, copies baseline ASM calls as the unchanged-gene baseline,
# and then lets polyphase_v2.sh rerun downstream steps only for genes whose
# rescue manifest status is written.

ROOT=${ROOT:-/data6/wangxuedong/polyploid_hla}
REALSETS=${REALSETS:-/data2/wangxuedong/polyploid-hla-realsets}
FQ_ROOT=${FQ_ROOT:-${ROOT}/fqs/amp/GenDx}
BASE_SPECHLA_ROOT=${BASE_SPECHLA_ROOT:-${REALSETS}/spechla_out_gendx_amp_abc_20260520}
BASE_ASM_ROOT=${BASE_ASM_ROOT:-${REALSETS}/asm_v2_gendx_amp_abc_20260520}
VALID_SPECHLA_ROOT=${VALID_SPECHLA_ROOT:-${REALSETS}/spechla_out_gendx_amp_abc_readbin_rescue_20260522}
VALID_ASM_ROOT=${VALID_ASM_ROOT:-${REALSETS}/asm_v2_gendx_amp_abc_readbin_rescue_20260522}
LOG_ROOT=${LOG_ROOT:-${REALSETS}/logs_gendx_amp_abc_readbin_rescue_20260522}
DIAG_DIR=${DIAG_DIR:-${ROOT}/diagnostics/gendx_read_bin_rescue_validation_20260522}
QUARTET_SUMMARY=${QUARTET_SUMMARY:-${ROOT}/diagnostics/gendx_quartet_summary_noA4.tsv}
PYBIN=${PYBIN:-python}
THREADS=${THREADS:-8}
GENES=${GENES:-"HLA-A HLA-B HLA-C HLA-DRB1 HLA-DPB1 HLA-DQB1"}
EXCLUDE_SAMPLES=${EXCLUDE_SAMPLES:-}
SAMPLES=${SAMPLES:-}
RESET_SAMPLE=${RESET_SAMPLE:-0}
STREAM_LOG=${STREAM_LOG:-1}
REUSE_RESCUE_MANIFEST=${REUSE_RESCUE_MANIFEST:-0}

# Truth-free rescue gates. The defaults are intentionally conservative, but the
# minimum rescued-pair threshold is lower than the generic helper default because
# GenDx amplicons can have very uneven gene-unique support by locus.
READ_BIN_RESCUE_K=${READ_BIN_RESCUE_K:-31}
READ_BIN_RESCUE_MIN_HITS=${READ_BIN_RESCUE_MIN_HITS:-1}
READ_BIN_RESCUE_MIN_MARGIN=${READ_BIN_RESCUE_MIN_MARGIN:-1}
READ_BIN_RESCUE_MAX_FRACTION=${READ_BIN_RESCUE_MAX_FRACTION:-0.25}
READ_BIN_RESCUE_MAX_PAIRS=${READ_BIN_RESCUE_MAX_PAIRS:-100000}
READ_BIN_RESCUE_RETENTION_MIN_FULL_PAIRS=${READ_BIN_RESCUE_RETENTION_MIN_FULL_PAIRS:-50}
READ_BIN_RESCUE_RETENTION_MAX_RETAINED_FRACTION=${READ_BIN_RESCUE_RETENTION_MAX_RETAINED_FRACTION:-0.10}
READ_BIN_RESCUE_RETENTION_MIN_MISSING_FRACTION=${READ_BIN_RESCUE_RETENTION_MIN_MISSING_FRACTION:-0.30}
READ_BIN_RESCUE_RETENTION_MIN_RESCUE_PAIRS=${READ_BIN_RESCUE_RETENTION_MIN_RESCUE_PAIRS:-5}
READ_BIN_RESCUE_MAX_RERUN_RESCUED_PAIRS=${READ_BIN_RESCUE_MAX_RERUN_RESCUED_PAIRS:-3000}
READ_BIN_RESCUE_MAX_RERUN_ORIGINAL_PAIRS=${READ_BIN_RESCUE_MAX_RERUN_ORIGINAL_PAIRS:-10000}

# Bounded DPB1 validation branch. These safeguards are truth-free and only
# affect this independent read-bin rescue validation runner unless callers pass
# the same environment knobs to the main pipeline explicitly.
READ_BIN_RESCUE_DPB1_MAX_RESCUED_PAIRS=${READ_BIN_RESCUE_DPB1_MAX_RESCUED_PAIRS:-1000}
DPB1_ASSEMBLE_PREFILTER_TOP=${DPB1_ASSEMBLE_PREFILTER_TOP:-80}
DPB1_ASSEMBLE_TOP_N_PER_BLOCK=${DPB1_ASSEMBLE_TOP_N_PER_BLOCK:-6}
DPB1_ASSEMBLE_GLOBAL_POOL_CAP=${DPB1_ASSEMBLE_GLOBAL_POOL_CAP:-20}
DPB1_ASSEMBLE_MAX_BLOCKS=${DPB1_ASSEMBLE_MAX_BLOCKS:-40}
DPB1_ASSEMBLE_MAX_BLOCK_HAPS=${DPB1_ASSEMBLE_MAX_BLOCK_HAPS:-160}
DPB1_ASSEMBLE_MAX_PAIRED_COMBOS=${DPB1_ASSEMBLE_MAX_PAIRED_COMBOS:-50000}
DPB1_ASSEMBLE_MAX_WALL_SECONDS=${DPB1_ASSEMBLE_MAX_WALL_SECONDS:-600}

SCRIPTS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
mkdir -p "$VALID_SPECHLA_ROOT" "$VALID_ASM_ROOT" "$LOG_ROOT" "$DIAG_DIR"

is_excluded() {
    local sample="$1"
    for excluded in $EXCLUDE_SAMPLES; do
        [[ "$sample" == "$excluded" ]] && return 0
    done
    return 1
}

discover_samples() {
    if [[ -n "$SAMPLES" ]]; then
        for sample in $SAMPLES; do
            is_excluded "$sample" || printf '%s\n' "$sample"
        done
        return
    fi
    for fq1 in "$FQ_ROOT"/*_R1_001.fastq.gz; do
        [[ -e "$fq1" ]] || continue
        sample=$(basename "$fq1" _R1_001.fastq.gz)
        is_excluded "$sample" && continue
        printf '%s\n' "$sample"
    done | sort
}

copy_if_present() {
    local src="$1"
    local dst="$2"
    if [[ -e "$src" && ! -e "$dst" ]]; then
        cp -p "$src" "$dst"
    fi
}

prepare_validation_sample() {
    local sample="$1"
    local base_out="${BASE_SPECHLA_ROOT}/${sample}"
    local base_asm="${BASE_ASM_ROOT}/${sample}"
    local valid_out="${VALID_SPECHLA_ROOT}/${sample}"
    local valid_asm="${VALID_ASM_ROOT}/${sample}"
    if [[ ! -d "$base_out" ]]; then
        echo "[ERROR] missing baseline output: $base_out" >&2
        return 1
    fi
    if [[ "$RESET_SAMPLE" == "1" ]]; then
        rm -rf "$valid_out" "$valid_asm"
    fi
    mkdir -p "$valid_out"
    if [[ -d "$base_asm" && ! -d "$valid_asm" ]]; then
        mkdir -p "$valid_asm"
        cp -a "${base_asm}/." "$valid_asm/"
    fi

    copy_if_present "${base_out}/${sample}.uniq.R1.fq.gz" "${valid_out}/${sample}.uniq.R1.fq.gz"
    copy_if_present "${base_out}/${sample}.uniq.R2.fq.gz" "${valid_out}/${sample}.uniq.R2.fq.gz"
    copy_if_present "${base_out}/${sample}.map_database.bam" "${valid_out}/${sample}.map_database.bam"
    copy_if_present "${base_out}/${sample}.map_database.bam.bai" "${valid_out}/${sample}.map_database.bam.bai"

    for short in A B C DRB1 DPB1 DQB1; do
        copy_if_present "${base_out}/${short}.R1.fq.gz" "${valid_out}/${short}.R1.fq.gz"
        copy_if_present "${base_out}/${short}.R2.fq.gz" "${valid_out}/${short}.R2.fq.gz"
    done
}

write_validation_gene_bed() {
    local sample="$1"
    local manifest="${DIAG_DIR}/${sample}.read_bin_rescue.tsv"
    local bed="${DIAG_DIR}/${sample}.written_genes.bed"
    local budget="${DIAG_DIR}/${sample}.gene_budget.tsv"
    "$PYBIN" - "$manifest" "$SCRIPTS_DIR/gene.spechla.bed" "$bed" "$budget" \
        "$READ_BIN_RESCUE_MAX_RERUN_RESCUED_PAIRS" \
        "$READ_BIN_RESCUE_MAX_RERUN_ORIGINAL_PAIRS" \
        "$READ_BIN_RESCUE_DPB1_MAX_RESCUED_PAIRS" "$BASE_SPECHLA_ROOT" "$sample" \
        "$DPB1_ASSEMBLE_MAX_BLOCKS" "$DPB1_ASSEMBLE_MAX_BLOCK_HAPS" <<'PY'
import csv
import sys
from pathlib import Path

import pysam

manifest = Path(sys.argv[1])
source_bed = Path(sys.argv[2])
out_bed = Path(sys.argv[3])
budget_path = Path(sys.argv[4])
max_rerun_rescued_pairs = int(sys.argv[5])
max_rerun_original_pairs = int(sys.argv[6])
dpb1_max_rescued_pairs = int(sys.argv[7])
base_spechla_root = Path(sys.argv[8])
sample_name = sys.argv[9]
dpb1_max_blocks = int(sys.argv[10])
dpb1_max_block_haps = int(sys.argv[11])

PLOIDY = 4


def load_gene_regions(path):
    regions = {}
    with path.open() as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 4:
                regions[parts[3]] = (parts[0], int(parts[1]), int(parts[2]))
    return regions


def count_phased_blocks(vcf_path, region):
    if not vcf_path.exists() or not Path(str(vcf_path) + ".tbi").exists():
        return None
    chrom, start, end = region
    block_ids = set()
    with pysam.VariantFile(str(vcf_path)) as vcf:
        for rec in vcf.fetch(chrom, start, end):
            call = rec.samples[0]
            gt = call.get("GT")
            if gt is None or len(gt) != PLOIDY or any(allele is None for allele in gt):
                continue
            if not call.phased:
                continue
            ps = call.get("PS")
            block_ids.add(ps if ps is not None else rec.pos)
    return len(block_ids)

written = set()
budget_rows = []
regions = load_gene_regions(source_bed)
dpb1_baseline_blocks = None
if "HLA-DPB1" in regions:
    dpb1_vcf = base_spechla_root / sample_name / f"{sample_name}.phased.hla-dpb1.vcf.gz"
    dpb1_baseline_blocks = count_phased_blocks(dpb1_vcf, regions["HLA-DPB1"])
with manifest.open() as handle:
    reader = csv.DictReader((line for line in handle if not line.startswith("#")), delimiter="\t")
    for row in reader:
        gene = row.get("gene", "")
        status = row.get("status", "")
        original_pairs = int(row.get("original_pairs") or 0)
        full_gene_pairs = int(row.get("full_gene_pairs") or 0)
        retained_gene_pairs = int(row.get("retained_gene_pairs") or 0)
        rescued_pairs = int(row.get("rescued_pairs") or 0)
        decision = status
        include = status == "written"
        if include and gene == "HLA-DPB1" and dpb1_max_rescued_pairs > 0 \
                and rescued_pairs > dpb1_max_rescued_pairs:
            include = False
            decision = f"skipped_dpb1_rescued_pairs>{dpb1_max_rescued_pairs}"
        if include and gene == "HLA-DPB1" and dpb1_baseline_blocks is not None:
            block_haps = dpb1_baseline_blocks * PLOIDY
            if dpb1_max_blocks > 0 and dpb1_baseline_blocks > dpb1_max_blocks:
                include = False
                decision = f"skipped_dpb1_blocks>{dpb1_max_blocks}"
            elif dpb1_max_block_haps > 0 and block_haps > dpb1_max_block_haps:
                include = False
                decision = f"skipped_dpb1_block_haps>{dpb1_max_block_haps}"
        if include and max_rerun_original_pairs > 0 and original_pairs > max_rerun_original_pairs:
            include = False
            decision = f"skipped_original_pairs>{max_rerun_original_pairs}"
        if include and max_rerun_rescued_pairs > 0 and rescued_pairs > max_rerun_rescued_pairs:
            include = False
            decision = f"skipped_rescued_pairs>{max_rerun_rescued_pairs}"
        if include:
            written.add(gene)
            decision = "included"
        budget_rows.append({
            "gene": gene,
            "original_pairs": original_pairs,
            "full_gene_pairs": full_gene_pairs,
            "retained_gene_pairs": retained_gene_pairs,
            "rescued_pairs": rescued_pairs,
            "baseline_blocks": dpb1_baseline_blocks if gene == "HLA-DPB1" and dpb1_baseline_blocks is not None else "",
            "manifest_status": status,
            "budget_decision": decision,
        })

with budget_path.open("w", newline="") as handle:
    fields = [
        "gene", "original_pairs", "full_gene_pairs", "retained_gene_pairs",
        "rescued_pairs", "baseline_blocks", "manifest_status", "budget_decision",
    ]
    writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
    writer.writeheader()
    writer.writerows(budget_rows)

with source_bed.open() as src, out_bed.open("w") as out:
    for line in src:
        parts = line.rstrip("\n").split("\t")
        if len(parts) >= 4 and parts[3] in written:
            out.write(line)
PY
    if [[ ! -s "$bed" ]]; then
        echo "[sample] $sample no written rescued genes; keep copied baseline calls"
        return 1
    fi
    echo "$bed"
}

short_to_gene() {
    case "$1" in
        A) echo "HLA-A" ;;
        B) echo "HLA-B" ;;
        C) echo "HLA-C" ;;
        DRB1) echo "HLA-DRB1" ;;
        DPB1) echo "HLA-DPB1" ;;
        DQB1) echo "HLA-DQB1" ;;
        *) echo "$1" ;;
    esac
}

copy_unchanged_gene_bams() {
    local sample="$1"
    local gene_bed="$2"
    local base_out="${BASE_SPECHLA_ROOT}/${sample}"
    local valid_out="${VALID_SPECHLA_ROOT}/${sample}"
    local valid_asm="${VALID_ASM_ROOT}/${sample}"
    local rerun_genes
    rerun_genes="$(awk '{print $4}' "$gene_bed" | tr '\n' ' ')"
    for short in A B C DRB1 DPB1 DQB1; do
        local gene
        gene="$(short_to_gene "$short")"
        if [[ " $rerun_genes " == *" $gene "* ]]; then
            rm -rf "${valid_asm}/$(echo "$gene" | tr '[:upper:]' '[:lower:]')"
            continue
        fi
        copy_if_present "${base_out}/${short}.bam" "${valid_out}/${short}.bam"
        copy_if_present "${base_out}/${short}.bam.bai" "${valid_out}/${short}.bam.bai"
    done
}

run_rescue_for_sample() {
    local sample="$1"
    local fq1="${FQ_ROOT}/${sample}_R1_001.fastq.gz"
    local fq2="${FQ_ROOT}/${sample}_R2_001.fastq.gz"
    local valid_out="${VALID_SPECHLA_ROOT}/${sample}"
    local manifest="${DIAG_DIR}/${sample}.read_bin_rescue.tsv"
    local log="${DIAG_DIR}/${sample}.read_bin_rescue.log"
    if [[ "$REUSE_RESCUE_MANIFEST" == "1" && -s "$manifest" ]]; then
        echo "[sample] $sample reuse existing rescue manifest: $manifest"
        return 0
    fi
    local gene_args=()
    for gene in $GENES; do
        gene_args+=(--gene "$gene")
    done

    "$PYBIN" "$SCRIPTS_DIR/rescue_gene_binned_reads.py" \
        --fq1 "$fq1" \
        --fq2 "$fq2" \
        --fq-dir "$valid_out" \
        --k "$READ_BIN_RESCUE_K" \
        --min-hits "$READ_BIN_RESCUE_MIN_HITS" \
        --min-margin "$READ_BIN_RESCUE_MIN_MARGIN" \
        --require-both-mates \
        --retention-gate \
        --retention-min-full-pairs "$READ_BIN_RESCUE_RETENTION_MIN_FULL_PAIRS" \
        --retention-max-retained-fraction "$READ_BIN_RESCUE_RETENTION_MAX_RETAINED_FRACTION" \
        --retention-min-missing-fraction "$READ_BIN_RESCUE_RETENTION_MIN_MISSING_FRACTION" \
        --retention-min-rescue-pairs "$READ_BIN_RESCUE_RETENTION_MIN_RESCUE_PAIRS" \
        --max-rescue-fraction "$READ_BIN_RESCUE_MAX_FRACTION" \
        --max-rescue-pairs "$READ_BIN_RESCUE_MAX_PAIRS" \
        --manifest "$manifest" \
        "${gene_args[@]}" > "$log" 2>&1
}

summarize_manifests() {
    local summary="${DIAG_DIR}/read_bin_rescue_manifest_summary.tsv"
    local budget_summary="${DIAG_DIR}/read_bin_rescue_budget_summary.tsv"
    "$PYBIN" - "${DIAG_DIR}" "$summary" <<'PY'
import csv
import sys
from pathlib import Path

diag_dir = Path(sys.argv[1])
summary = Path(sys.argv[2])
fields = [
    "sample", "gene", "original_pairs", "full_gene_pairs", "retained_gene_pairs",
    "retained_fraction", "missing_fraction", "rescued_pairs", "allowed_pairs",
    "status", "retention_reason",
]
with summary.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
    writer.writeheader()
    for manifest in sorted(diag_dir.glob("*.read_bin_rescue.tsv")):
        sample = manifest.name.removesuffix(".read_bin_rescue.tsv")
        with manifest.open(newline="") as source:
            reader = csv.DictReader((line for line in source if not line.startswith("#")), delimiter="\t")
            for row in reader:
                writer.writerow({field: sample if field == "sample" else row.get(field, "") for field in fields})
PY
    echo "[summary] $summary"
    "$PYBIN" - "${DIAG_DIR}" "$budget_summary" <<'PY'
import csv
import sys
from pathlib import Path

diag_dir = Path(sys.argv[1])
summary = Path(sys.argv[2])
fields = [
    "sample", "gene", "original_pairs", "full_gene_pairs", "retained_gene_pairs",
    "rescued_pairs", "baseline_blocks", "manifest_status", "budget_decision",
]
with summary.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
    writer.writeheader()
    for budget in sorted(diag_dir.glob("*.gene_budget.tsv")):
        sample = budget.name.removesuffix(".gene_budget.tsv")
        with budget.open(newline="") as source:
            reader = csv.DictReader(source, delimiter="\t")
            for row in reader:
                out = {field: row.get(field, "") for field in fields}
                out["sample"] = sample
                writer.writerow(out)
PY
    echo "[summary] $budget_summary"
}

run_pipeline_for_sample() {
    local sample="$1"
    local gene_bed="$2"
    local fq1="${FQ_ROOT}/${sample}_R1_001.fastq.gz"
    local fq2="${FQ_ROOT}/${sample}_R2_001.fastq.gz"
    local log="${LOG_ROOT}/${sample}.polyphase_v2.log"

    if [[ "$STREAM_LOG" == "1" ]]; then
        env \
            FQ1="$fq1" \
            FQ2="$fq2" \
            SAMPLE="$sample" \
            OUT_ROOT="$VALID_SPECHLA_ROOT" \
            ASM_ROOT="$VALID_ASM_ROOT" \
            GENE_BED="$gene_bed" \
            THREADS="$THREADS" \
            SKIP_DONE=1 \
            EXON_TYPING=0 \
                DPB1_ASSEMBLE_PREFILTER_TOP="$DPB1_ASSEMBLE_PREFILTER_TOP" \
                DPB1_ASSEMBLE_TOP_N_PER_BLOCK="$DPB1_ASSEMBLE_TOP_N_PER_BLOCK" \
                DPB1_ASSEMBLE_GLOBAL_POOL_CAP="$DPB1_ASSEMBLE_GLOBAL_POOL_CAP" \
                DPB1_ASSEMBLE_MAX_BLOCKS="$DPB1_ASSEMBLE_MAX_BLOCKS" \
                DPB1_ASSEMBLE_MAX_BLOCK_HAPS="$DPB1_ASSEMBLE_MAX_BLOCK_HAPS" \
                DPB1_ASSEMBLE_MAX_PAIRED_COMBOS="$DPB1_ASSEMBLE_MAX_PAIRED_COMBOS" \
                DPB1_ASSEMBLE_MAX_WALL_SECONDS="$DPB1_ASSEMBLE_MAX_WALL_SECONDS" \
            PYBIN="$PYBIN" \
            bash "$SCRIPTS_DIR/polyphase_v2.sh" 2>&1 | tee "$log"
    else
        env \
            FQ1="$fq1" \
            FQ2="$fq2" \
            SAMPLE="$sample" \
            OUT_ROOT="$VALID_SPECHLA_ROOT" \
            ASM_ROOT="$VALID_ASM_ROOT" \
            GENE_BED="$gene_bed" \
            THREADS="$THREADS" \
            SKIP_DONE=1 \
            EXON_TYPING=0 \
                DPB1_ASSEMBLE_PREFILTER_TOP="$DPB1_ASSEMBLE_PREFILTER_TOP" \
                DPB1_ASSEMBLE_TOP_N_PER_BLOCK="$DPB1_ASSEMBLE_TOP_N_PER_BLOCK" \
                DPB1_ASSEMBLE_GLOBAL_POOL_CAP="$DPB1_ASSEMBLE_GLOBAL_POOL_CAP" \
                DPB1_ASSEMBLE_MAX_BLOCKS="$DPB1_ASSEMBLE_MAX_BLOCKS" \
                DPB1_ASSEMBLE_MAX_BLOCK_HAPS="$DPB1_ASSEMBLE_MAX_BLOCK_HAPS" \
                DPB1_ASSEMBLE_MAX_PAIRED_COMBOS="$DPB1_ASSEMBLE_MAX_PAIRED_COMBOS" \
                DPB1_ASSEMBLE_MAX_WALL_SECONDS="$DPB1_ASSEMBLE_MAX_WALL_SECONDS" \
            PYBIN="$PYBIN" \
            bash "$SCRIPTS_DIR/polyphase_v2.sh" > "$log" 2>&1
    fi
}

samples_file="${DIAG_DIR}/samples.tsv"
discover_samples > "$samples_file"
echo "[samples] $(wc -l < "$samples_file") samples -> $samples_file"

while read -r sample; do
    [[ -n "$sample" ]] || continue
    echo "[sample] $sample prepare"
    prepare_validation_sample "$sample"
    echo "[sample] $sample rescue"
    run_rescue_for_sample "$sample"
    if gene_bed=$(write_validation_gene_bed "$sample"); then
        echo "[sample] $sample downstream rerun genes: $(awk '{print $4}' "$gene_bed" | paste -sd, -)"
        copy_unchanged_gene_bams "$sample" "$gene_bed"
        run_pipeline_for_sample "$sample" "$gene_bed"
    fi
done < "$samples_file"

summarize_manifests

"$PYBIN" "$SCRIPTS_DIR/score_gendx_validation_root.py" \
    --quartet-summary "$QUARTET_SUMMARY" \
    --asm-root "$VALID_ASM_ROOT" \
    --out "${DIAG_DIR}/read_bin_rescue_score.tsv" \
    --summary "${DIAG_DIR}/read_bin_rescue_score.summary.tsv" \
    --samples-file "$samples_file"

echo "[score] ${DIAG_DIR}/read_bin_rescue_score.summary.tsv"
echo "[done] validation spechla root: $VALID_SPECHLA_ROOT"
echo "[done] validation asm root: $VALID_ASM_ROOT"