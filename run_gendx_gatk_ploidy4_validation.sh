#!/bin/bash
set -euo pipefail

# Independent GenDx noA4 validation branch for GATK variant callers.
#
# This does not modify the baseline SpecHLA/ASM roots or polyphase_v2.sh.
# It reuses the existing per-sample HLA BAMs/chi estimates, replaces only the
# variant-calling step with GATK HaplotypeCaller --sample-ploidy 4 by default,
# or Mutect2 when GATK_CALLER=Mutect2, then applies the existing chi_R-aware
# count reassignment before whatshap/polyphase and hla_polyphase_assemble.py.

ROOT=${ROOT:-/data6/wangxuedong/polyploid_hla}
REALSETS=${REALSETS:-/data2/wangxuedong/polyploid-hla-realsets}
BASE_SPECHLA_ROOT=${BASE_SPECHLA_ROOT:-${REALSETS}/spechla_out_gendx_amp_abc_readbin_rescue_full_conda_20260522}
VALID_SPECHLA_ROOT=${VALID_SPECHLA_ROOT:-${REALSETS}/spechla_out_gendx_amp_abc_gatk_ploidy4_chi_count_20260526}
VALID_ASM_ROOT=${VALID_ASM_ROOT:-${REALSETS}/asm_v2_gendx_amp_abc_gatk_ploidy4_chi_count_20260526}
DIAG_DIR=${DIAG_DIR:-${ROOT}/diagnostics/gendx_gatk_ploidy4_chi_count_20260526}
QUARTET_SUMMARY=${QUARTET_SUMMARY:-${ROOT}/diagnostics/gendx_quartet_summary_noA4.tsv}
SAMPLES_FILE=${SAMPLES_FILE:-${ROOT}/diagnostics/gendx_read_bin_rescue_full_conda_20260522/noA4_samples.tsv}
SAMPLES=${SAMPLES:-}
GENES=${GENES:-"HLA-A HLA-B HLA-C HLA-DRB1 HLA-DQB1 HLA-DPB1"}
SAMPLE_JOBS=${SAMPLE_JOBS:-1}
SKIP_SCORE=${SKIP_SCORE:-0}
SCORE_ONLY=${SCORE_ONLY:-0}

PYBIN=${PYBIN:-/data3/wangxuedong/app/miniconda3/envs/polyploid-hla/bin/python}
WHATSHAP=${WHATSHAP:-$(command -v whatshap || true)}
GATK=${GATK:-$(command -v gatk || true)}
GATK_JAR=${GATK_JAR:-}
FREEBAYES=${FREEBAYES:-$(command -v freebayes || true)}
PICARD=${PICARD:-$(command -v picard || true)}
THREADS=${THREADS:-8}
SAMTOOLS_THREADS=${SAMTOOLS_THREADS:-$THREADS}
GATK_JAVA_OPTIONS=${GATK_JAVA_OPTIONS:--Xmx24g}
PLOIDY=${PLOIDY:-4}
GATK_CALLER=${GATK_CALLER:-HaplotypeCaller}
CALLER_TAG=${CALLER_TAG:-gatk}
SCORE_PREFIX=${SCORE_PREFIX:-gatk_ploidy4_chi_count}

GATK_STAND_CALL_CONF=${GATK_STAND_CALL_CONF:-10}
GATK_MIN_BQ=${GATK_MIN_BQ:-13}
MUTECT2_MIN_AF=${MUTECT2_MIN_AF:-0.005}
MUTECT2_INITIAL_TUMOR_LOD=${MUTECT2_INITIAL_TUMOR_LOD:-0.0}
MUTECT2_TUMOR_LOD_TO_EMIT=${MUTECT2_TUMOR_LOD_TO_EMIT:-0.0}
MUTECT2_MAX_READS_PER_ALIGNMENT_START=${MUTECT2_MAX_READS_PER_ALIGNMENT_START:-0}
FB_MIN_BQ=${FB_MIN_BQ:-13}
FB_MIN_MQ=${FB_MIN_MQ:-20}
GT_MIN_DEPTH=${GT_MIN_DEPTH:-10}
GT_DROP_FP_AF=${GT_DROP_FP_AF:-0.05}
GT_DROP_UNASSIGNABLE=${GT_DROP_UNASSIGNABLE:-0}
MASK_MIN_DEPTH=${MASK_MIN_DEPTH:-5}

ASSEMBLE_ALIGNER=${ASSEMBLE_ALIGNER:-parasail}
ASSEMBLE_PREFILTER_TOP=${ASSEMBLE_PREFILTER_TOP:-200}
ASSEMBLE_TOP_N_PER_BLOCK=${ASSEMBLE_TOP_N_PER_BLOCK:-10}
ASSEMBLE_GLOBAL_POOL_CAP=${ASSEMBLE_GLOBAL_POOL_CAP:-30}
ASSEMBLE_MAX_BLOCKS=${ASSEMBLE_MAX_BLOCKS:-0}
ASSEMBLE_MAX_BLOCK_HAPS=${ASSEMBLE_MAX_BLOCK_HAPS:-0}
ASSEMBLE_MAX_PAIRED_COMBOS=${ASSEMBLE_MAX_PAIRED_COMBOS:-0}
ASSEMBLE_MAX_WALL_SECONDS=${ASSEMBLE_MAX_WALL_SECONDS:-0}

# DPB1 can become pathological when a caller emits many small phase blocks.
# These defaults are bounded for validation; set to 0 to remove a cap.
DPB1_ASSEMBLE_PREFILTER_TOP=${DPB1_ASSEMBLE_PREFILTER_TOP:-80}
DPB1_ASSEMBLE_TOP_N_PER_BLOCK=${DPB1_ASSEMBLE_TOP_N_PER_BLOCK:-6}
DPB1_ASSEMBLE_GLOBAL_POOL_CAP=${DPB1_ASSEMBLE_GLOBAL_POOL_CAP:-20}
DPB1_ASSEMBLE_MAX_BLOCKS=${DPB1_ASSEMBLE_MAX_BLOCKS:-30}
DPB1_ASSEMBLE_MAX_BLOCK_HAPS=${DPB1_ASSEMBLE_MAX_BLOCK_HAPS:-120}
DPB1_ASSEMBLE_MAX_PAIRED_COMBOS=${DPB1_ASSEMBLE_MAX_PAIRED_COMBOS:-80000}
DPB1_ASSEMBLE_MAX_WALL_SECONDS=${DPB1_ASSEMBLE_MAX_WALL_SECONDS:-900}

# GATK can fragment DRB1 similarly in chimeric samples, making uncapped noA4
# validation impractically slow. Keep this branch bounded without changing the
# production pipeline.
DRB1_ASSEMBLE_PREFILTER_TOP=${DRB1_ASSEMBLE_PREFILTER_TOP:-80}
DRB1_ASSEMBLE_TOP_N_PER_BLOCK=${DRB1_ASSEMBLE_TOP_N_PER_BLOCK:-6}
DRB1_ASSEMBLE_GLOBAL_POOL_CAP=${DRB1_ASSEMBLE_GLOBAL_POOL_CAP:-20}
DRB1_ASSEMBLE_MAX_BLOCKS=${DRB1_ASSEMBLE_MAX_BLOCKS:-30}
DRB1_ASSEMBLE_MAX_BLOCK_HAPS=${DRB1_ASSEMBLE_MAX_BLOCK_HAPS:-120}
DRB1_ASSEMBLE_MAX_PAIRED_COMBOS=${DRB1_ASSEMBLE_MAX_PAIRED_COMBOS:-120000}
DRB1_ASSEMBLE_MAX_WALL_SECONDS=${DRB1_ASSEMBLE_MAX_WALL_SECONDS:-900}

SCRIPTS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUNDLED_SPECHLA=${BUNDLED_SPECHLA:-${SCRIPTS_DIR}/resources/spechla}
SPECHLA=${SPECHLA:-$BUNDLED_SPECHLA}
SPECHLA_DB=${SPECHLA_DB:-${SPECHLA}/db}
HLA_REF=${HLA_REF:-${SPECHLA_DB}/ref/hla.ref.extend.fa}
DB_PREFIX=${DB_PREFIX:-${SPECHLA_DB}/ref/hla_gen.format.filter.extend.DRB.no26789.v2.fasta}
GENE_BED=${GENE_BED:-${SCRIPTS_DIR}/gene.spechla.bed}
ASSEMBLE_PY=${ASSEMBLE_PY:-${SCRIPTS_DIR}/hla_polyphase_assemble.py}
G_GROUP=${G_GROUP:-${SPECHLA_DB}/HLA/hla_nom_g.txt}

mkdir -p "$VALID_SPECHLA_ROOT" "$VALID_ASM_ROOT" "$DIAG_DIR"

die() {
    echo "[ERROR] $*" >&2
    exit 1
}

vcf_ready() {
    local vcf="$1"
    [[ -s "$vcf" && -s "${vcf}.tbi" ]]
}

asm_ready() {
    local asm_out="$1" gene="$2"
    [[ -s "${asm_out}/${gene}/calls.tsv" ]]
}

run_gatk() {
    if [[ -n "$GATK_JAR" ]]; then
        java $GATK_JAVA_OPTIONS -jar "$GATK_JAR" "$@"
    elif [[ -n "$GATK" ]]; then
        "$GATK" --java-options "$GATK_JAVA_OPTIONS" "$@"
    else
        die "GATK not found. Set GATK=/path/gatk or GATK_JAR=/path/gatk-package.jar."
    fi
}

require_tools() {
    [[ -x "$PYBIN" ]] || die "PYBIN is not executable: $PYBIN"
    [[ -n "$WHATSHAP" ]] || die "whatshap not found; set WHATSHAP=/path/whatshap"
    command -v bwa >/dev/null 2>&1 || die "bwa not found"
    [[ -n "$FREEBAYES" ]] || die "freebayes not found; set FREEBAYES=/path/freebayes"
    command -v bcftools >/dev/null 2>&1 || die "bcftools not found"
    command -v tabix >/dev/null 2>&1 || die "tabix not found"
    command -v samtools >/dev/null 2>&1 || die "samtools not found"
    [[ -n "$GATK" || -n "$GATK_JAR" ]] || die "GATK not found on PATH and GATK_JAR not set"
    [[ "$GATK_CALLER" == "HaplotypeCaller" || "$GATK_CALLER" == "Mutect2" ]] || die "unsupported GATK_CALLER=$GATK_CALLER"
}

ensure_reference_indexes() {
    [[ -s "$HLA_REF" ]] || die "missing reference: $HLA_REF"
    if [[ ! -s "${HLA_REF}.fai" ]]; then
        samtools faidx "$HLA_REF"
    fi
    local dict="${HLA_REF%.*}.dict"
    if [[ ! -s "$dict" ]]; then
        [[ -n "$PICARD" ]] || die "missing $dict and picard not found to create it"
        "$PICARD" CreateSequenceDictionary R="$HLA_REF" O="$dict"
    fi
}

discover_samples() {
    if [[ -n "$SAMPLES" ]]; then
        for sample in $SAMPLES; do printf '%s\n' "$sample"; done
    else
        [[ -s "$SAMPLES_FILE" ]] || die "missing SAMPLES_FILE: $SAMPLES_FILE"
        grep -v '^[[:space:]]*$' "$SAMPLES_FILE"
    fi
}

copy_if_present() {
    local src="$1" dst="$2"
    if [[ -e "$src" && ! -e "$dst" ]]; then
        cp -Pp "$src" "$dst"
    fi
}

prepare_sample() {
    local sample="$1"
    local base_out="${BASE_SPECHLA_ROOT}/${sample}"
    local out="${VALID_SPECHLA_ROOT}/${sample}"
    [[ -d "$base_out" ]] || die "missing baseline sample output: $base_out"
    mkdir -p "$out" "${VALID_ASM_ROOT}/${sample}"

    for suffix in \
        "uniq.R1.fq.gz" "uniq.R2.fq.gz" \
        "map_database.bam" "map_database.bam.bai" \
        "merge.bam" "merge.bam.bai" \
        "chimerism.txt" "chi_pooled.txt"; do
        copy_if_present "${base_out}/${sample}.${suffix}" "${out}/${sample}.${suffix}"
    done
    for short in A B C DPB1 DQB1 DRB1; do
        copy_if_present "${base_out}/${short}.R1.fq.gz" "${out}/${short}.R1.fq.gz"
        copy_if_present "${base_out}/${short}.R2.fq.gz" "${out}/${short}.R2.fq.gz"
        copy_if_present "${base_out}/${short}.bam" "${out}/${short}.bam"
        copy_if_present "${base_out}/${short}.bam.bai" "${out}/${short}.bam.bai"
    done
}

ensure_merged_bam() {
    local sample="$1" out="$2" merged_bam="$3"
    local hla fq1 fq2 group header per_bam per_ref

    if [[ -s "$merged_bam" ]]; then
        [[ -s "${merged_bam}.bai" ]] || samtools index -@ "$SAMTOOLS_THREADS" "$merged_bam"
        return
    fi

    fq1="${out}/${sample}.uniq.R1.fq.gz"
    fq2="${out}/${sample}.uniq.R2.fq.gz"
    [[ -s "$fq1" && -s "$fq2" ]] || die "$sample: missing uniq FASTQ needed to rebuild $merged_bam"

    group="@RG\tID:${sample}\tSM:${sample}"
    header="${out}/header.sam"
    echo "[$sample] rebuild missing merged BAM from per-gene FASTQs"
    bwa mem -U 10000 -L 10000,10000 -R "$group" \
        "$HLA_REF" "$fq1" "$fq2" 2>/dev/null \
        | samtools view -H - > "$header" || true

    local to_merge=()
    for hla in A B C DPB1 DQB1 DRB1; do
        per_bam="${out}/${hla}.bam"
        if [[ ! -s "$per_bam" ]]; then
            fq1="${out}/${hla}.R1.fq.gz"
            fq2="${out}/${hla}.R2.fq.gz"
            if [[ ! -s "$fq1" || ! -s "$fq2" ]]; then
                echo "[$sample] no reads for HLA-$hla; skip BAM rebuild"
                continue
            fi
            per_ref="${SPECHLA_DB}/HLA/HLA_${hla}/HLA_${hla}.fa"
            [[ -s "$per_ref" ]] || die "$sample: missing per-gene reference $per_ref"
            echo "[$sample] rebuild HLA-$hla BAM"
            bwa mem -t "$THREADS" -U 10000 -L 10000,10000 -R "$group" \
                "$per_ref" "$fq1" "$fq2" \
                | samtools view -@ "$SAMTOOLS_THREADS" -bS -F 0x800 - \
                | samtools sort -@ "$SAMTOOLS_THREADS" -o "$per_bam" -
        fi
        [[ -s "${per_bam}.bai" ]] || samtools index -@ "$SAMTOOLS_THREADS" "$per_bam"
        to_merge+=("$per_bam")
    done

    (( ${#to_merge[@]} > 0 )) || die "$sample: no per-gene BAMs available to merge"
    samtools merge -@ "$SAMTOOLS_THREADS" -f -h "$header" "$merged_bam" "${to_merge[@]}"
    samtools index -@ "$SAMTOOLS_THREADS" "$merged_bam"
}

read_chi_r() {
    local sample="$1" out="${VALID_SPECHLA_ROOT}/${sample}"
    local chi=""
    if [[ -s "${out}/${sample}.chi_pooled.txt" ]]; then
        chi=$(awk '/^GLOBAL[[:space:]]+chi_R=/{for(i=1;i<=NF;i++)if($i~/^chi_R=/){split($i,a,"=");print a[2];exit}}' "${out}/${sample}.chi_pooled.txt")
        if [[ -n "$chi" ]] && awk -v x="$chi" 'BEGIN{exit !(x>0 && x<0.5)}'; then
            printf '%s\n' "$chi"
            return
        fi
    fi
    if [[ -s "${out}/${sample}.chimerism.txt" ]]; then
        chi=$(awk '/chi_R=/{for(i=1;i<=NF;i++)if($i~/^chi_R=/){split($i,a,"=");print a[2];exit}}' "${out}/${sample}.chimerism.txt")
    fi
    [[ -n "$chi" ]] || return 1
    printf '%s\n' "$chi"
}

estimate_chi_r_pooled_fallback() {
    local sample="$1" out="$2" merged_bam="$3"
    local pooled_vcf="${out}/${sample}.pooled_continuous.vcf.gz"
    local pooled_log="${out}/${sample}.chi_pooled.txt"
    local chi=""

    if ! vcf_ready "$pooled_vcf"; then
        echo "[$sample] estimate missing chi_R with freebayes --pooled-continuous" >&2
        "$FREEBAYES" \
            --pooled-continuous \
            --min-alternate-fraction 0.005 \
            --min-alternate-count 2 \
            --min-base-quality "$FB_MIN_BQ" \
            --min-mapping-quality "$FB_MIN_MQ" \
            --min-coverage 30 \
            --haplotype-length 0 \
            --use-best-n-alleles 2 \
            -f "$HLA_REF" "$merged_bam" \
        | bcftools norm -f "$HLA_REF" -a -m -any -Oz -o "$pooled_vcf"
        tabix -f -p vcf "$pooled_vcf"
    fi

    "$PYBIN" "${SCRIPTS_DIR}/estimate_chi_pooled.py" "$pooled_vcf" > "$pooled_log" 2>&1 || true
    chi=$(awk '/^GLOBAL[[:space:]]+chi_R=/{for(i=1;i<=NF;i++)if($i~/^chi_R=/){split($i,a,"=");print a[2];exit}}' "$pooled_log")
    if [[ -n "$chi" ]] && awk -v x="$chi" 'BEGIN{exit !(x>0 && x<0.5)}'; then
        printf '%s\n' "$chi"
        return
    fi
    die "$sample: no usable chi_R found after pooled-continuous fallback"
}

gene_region() {
    local gene="$1"
    awk -v g="$gene" '$4==g{printf "%s:%d-%d\n", $1, $2+1, $3; exit}' "$GENE_BED"
}

gene_tag_lc() {
    echo "$1" | tr '[:upper:]' '[:lower:]'
}

run_sample() {
    local sample="$1"
    local out="${VALID_SPECHLA_ROOT}/${sample}"
    local asm_sample="${VALID_ASM_ROOT}/${sample}"
    local merged_bam="${out}/${sample}.merge.bam"
    local chi_r

    echo "===================================================="
    echo "[$sample] $GATK_CALLER + chi_R/count validation"
    prepare_sample "$sample"
    ensure_merged_bam "$sample" "$out" "$merged_bam"
    chi_r=$(read_chi_r "$sample" || estimate_chi_r_pooled_fallback "$sample" "$out" "$merged_bam")
    echo "[$sample] chi_R=$chi_r"

    for gene in $GENES; do
        local tag_lc region raw_gene_vcf nopl_gene_vcf gene_vcf regt_vcf phased_vcf asm_out
        tag_lc=$(gene_tag_lc "$gene")
        region=$(gene_region "$gene")
        [[ -n "$region" ]] || die "missing region for $gene in $GENE_BED"
        raw_gene_vcf="${out}/${sample}.${CALLER_TAG}.raw.${tag_lc}.vcf.gz"
        nopl_gene_vcf="${out}/${sample}.${CALLER_TAG}.raw_nopl.${tag_lc}.vcf.gz"
        gene_vcf="${out}/${sample}.${CALLER_TAG}.${tag_lc}.vcf.gz"
        regt_vcf="${out}/${sample}.${CALLER_TAG}_regt.${tag_lc}.vcf.gz"
        phased_vcf="${out}/${sample}.${CALLER_TAG}_phased.${tag_lc}.vcf.gz"
        asm_out="${asm_sample}/${tag_lc}"
        mkdir -p "$asm_out"

        echo "[$sample] $gene slice -> chi_R count reassignment -> whatshap -> assemble"
        if ! vcf_ready "$gene_vcf"; then
            if ! vcf_ready "$raw_gene_vcf"; then
                if [[ "$GATK_CALLER" == "Mutect2" ]]; then
                    echo "[$sample] $gene run GATK Mutect2 region=$region min_af=$MUTECT2_MIN_AF"
                    run_gatk Mutect2 \
                        -R "$HLA_REF" \
                        -I "$merged_bam" \
                        -L "$region" \
                        -O "$raw_gene_vcf" \
                        --native-pair-hmm-threads "$THREADS" \
                        --min-base-quality-score "$GATK_MIN_BQ" \
                        --minimum-allele-fraction "$MUTECT2_MIN_AF" \
                        --initial-tumor-lod "$MUTECT2_INITIAL_TUMOR_LOD" \
                        --tumor-lod-to-emit "$MUTECT2_TUMOR_LOD_TO_EMIT" \
                        --max-reads-per-alignment-start "$MUTECT2_MAX_READS_PER_ALIGNMENT_START"
                else
                    echo "[$sample] $gene run GATK HaplotypeCaller ploidy=$PLOIDY region=$region"
                    run_gatk HaplotypeCaller \
                        -R "$HLA_REF" \
                        -I "$merged_bam" \
                        -L "$region" \
                        -O "$raw_gene_vcf" \
                        --sample-ploidy "$PLOIDY" \
                        --native-pair-hmm-threads "$THREADS" \
                        --min-base-quality-score "$GATK_MIN_BQ" \
                        --standard-min-confidence-threshold-for-calling "$GATK_STAND_CALL_CONF"
                fi
            fi
            if ! vcf_ready "$nopl_gene_vcf"; then
                bcftools annotate -x FORMAT/PL -Oz -o "$nopl_gene_vcf" "$raw_gene_vcf"
                tabix -f -p vcf "$nopl_gene_vcf"
            fi
            bcftools norm -f "$HLA_REF" -a -m -any -Oz -o "$gene_vcf" "$nopl_gene_vcf"
            tabix -f -p vcf "$gene_vcf"
        fi
        if ! vcf_ready "$regt_vcf"; then
            local reassign_extra=()
            if [[ "$GT_DROP_UNASSIGNABLE" == "1" ]]; then
                reassign_extra+=(--drop-unassignable)
            fi
            "$PYBIN" "${SCRIPTS_DIR}/reassign_gt_chimeric.py" \
                --vcf "$gene_vcf" \
                --out "$regt_vcf" \
                --chi-r "$chi_r" \
                --min-depth "$GT_MIN_DEPTH" \
                --drop-fp-af "$GT_DROP_FP_AF" \
                "${reassign_extra[@]}"
        fi
        if ! vcf_ready "$phased_vcf"; then
            "$WHATSHAP" polyphase \
                --ploidy "$PLOIDY" \
                --reference "$HLA_REF" \
                --threads "$THREADS" \
                --ignore-read-groups \
                --output "$phased_vcf" \
                "$regt_vcf" "$merged_bam"
            tabix -f -p vcf "$phased_vcf"
        fi

        local prefilter_top="$ASSEMBLE_PREFILTER_TOP"
        local top_n="$ASSEMBLE_TOP_N_PER_BLOCK"
        local pool_cap="$ASSEMBLE_GLOBAL_POOL_CAP"
        local max_blocks="$ASSEMBLE_MAX_BLOCKS"
        local max_block_haps="$ASSEMBLE_MAX_BLOCK_HAPS"
        local max_combos="$ASSEMBLE_MAX_PAIRED_COMBOS"
        local max_wall="$ASSEMBLE_MAX_WALL_SECONDS"
        if [[ "$gene" == "HLA-DPB1" ]]; then
            prefilter_top="$DPB1_ASSEMBLE_PREFILTER_TOP"
            top_n="$DPB1_ASSEMBLE_TOP_N_PER_BLOCK"
            pool_cap="$DPB1_ASSEMBLE_GLOBAL_POOL_CAP"
            max_blocks="$DPB1_ASSEMBLE_MAX_BLOCKS"
            max_block_haps="$DPB1_ASSEMBLE_MAX_BLOCK_HAPS"
            max_combos="$DPB1_ASSEMBLE_MAX_PAIRED_COMBOS"
            max_wall="$DPB1_ASSEMBLE_MAX_WALL_SECONDS"
        elif [[ "$gene" == "HLA-DRB1" ]]; then
            prefilter_top="$DRB1_ASSEMBLE_PREFILTER_TOP"
            top_n="$DRB1_ASSEMBLE_TOP_N_PER_BLOCK"
            pool_cap="$DRB1_ASSEMBLE_GLOBAL_POOL_CAP"
            max_blocks="$DRB1_ASSEMBLE_MAX_BLOCKS"
            max_block_haps="$DRB1_ASSEMBLE_MAX_BLOCK_HAPS"
            max_combos="$DRB1_ASSEMBLE_MAX_PAIRED_COMBOS"
            max_wall="$DRB1_ASSEMBLE_MAX_WALL_SECONDS"
        fi

        if asm_ready "$asm_out" "$gene"; then
            echo "[$sample] $gene skip existing assembly $asm_out/$gene/calls.tsv"
        else
            "$PYBIN" -u "$ASSEMBLE_PY" \
                --vcf "$phased_vcf" \
                --ref "$HLA_REF" \
                --gene-bed "$GENE_BED" \
                --genes "$gene" \
                --out "$asm_out" \
                --imgt "$DB_PREFIX" \
                --paired-diploids \
                --bam "$merged_bam" \
                --mask-min-depth "$MASK_MIN_DEPTH" \
                --aligner "$ASSEMBLE_ALIGNER" \
                --prefilter-top "$prefilter_top" \
                --top-n-per-block "$top_n" \
                --global-pool-cap "$pool_cap" \
                --max-blocks "$max_blocks" \
                --max-block-haps "$max_block_haps" \
                --max-paired-combos "$max_combos" \
                --max-wall-seconds "$max_wall" \
                --chimerism "$chi_r" \
                --dump-block-fa
        fi
    done

    local final="${asm_sample}/${sample}.final_calls.tsv"
    local compact="${asm_sample}/${sample}.final_calls.compact.tsv"
    "$PYBIN" "${SCRIPTS_DIR}/aggregate_calls.py" \
        --asm-root "$VALID_ASM_ROOT" \
        --sample "$sample" \
        --out "$final" \
        --spechla-root "$VALID_SPECHLA_ROOT" \
        --compact-out "$compact" \
        --g-group "$G_GROUP"
    echo "[$sample] final $final"
}

score_all() {
    local detail="${DIAG_DIR}/${SCORE_PREFIX}.noA4_score.tsv"
    local summary="${DIAG_DIR}/${SCORE_PREFIX}.noA4_score.summary.tsv"
    "$PYBIN" "${SCRIPTS_DIR}/score_gendx_validation_root.py" \
        --quartet-summary "$QUARTET_SUMMARY" \
        --asm-root "$VALID_ASM_ROOT" \
        --out "$detail" \
        --summary "$summary" \
        --samples-file "$SAMPLES_FILE"
    echo "[score] $summary"
    cat "$summary"
}

run_all_samples() {
    local sample failures
    failures=0

    if (( SAMPLE_JOBS <= 1 )); then
        while read -r sample; do
            [[ -z "$sample" ]] && continue
            run_sample "$sample"
        done < <(discover_samples)
        return
    fi

    while read -r sample; do
        [[ -z "$sample" ]] && continue
        while (( $(jobs -rp | wc -l) >= SAMPLE_JOBS )); do
            wait -n || failures=$((failures + 1))
        done
        run_sample "$sample" &
    done < <(discover_samples)

    while (( $(jobs -rp | wc -l) > 0 )); do
        wait -n || failures=$((failures + 1))
    done
    (( failures == 0 )) || die "$failures sample job(s) failed"
}

require_tools
ensure_reference_indexes
if [[ "$SCORE_ONLY" == "1" ]]; then
    score_all
    echo "[INFO] $GATK_CALLER noA4 scoring complete."
    exit 0
fi

run_all_samples
if [[ "$SKIP_SCORE" != "1" ]]; then
    score_all
fi

echo "[INFO] $GATK_CALLER noA4 validation complete."