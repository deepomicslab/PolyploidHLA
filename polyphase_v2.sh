#!/bin/bash
#SBATCH --job-name='polyphase_v2'
#SBATCH --cpus-per-task=8
#SBATCH --ntasks=1
#SBATCH --output=log.polyphase_v2.log
#SBATCH --mem=50G
#SBATCH --time=14-00:00:00
#
# End-to-end polyploid (k=4) HLA typing for chimeric / mixed samples.
#
# Pipeline (per sample):
#   0. dedupe read names (SpecHLA uniq_read_name.py)
#   1. competitive bowtie2/bwa map to IMGT db -> assign reads to gene
#      (SpecHLA assign_reads_to_genes.py)
#   2. per-gene bwa mem to SpecHLA per-gene ref (HLA_X.fa, ~5-13kb)
#   3. merge per-gene BAMs -> single sample BAM on hla.ref.extend.fa
#   4. freebayes -p 4 with low MAF / low alt-count (chimerism + low-donor friendly)
#   5. estimate chimerism chi from VCF AD (whole hla VCF)
#   6. per-gene: subset VCF -> whatshap polyphase (ploidy=4)
#                          -> hla_polyphase_assemble.py with chimerism + N-mask
#
# Configurable by env vars; sensible defaults. Run:
#   FQ1=... FQ2=... SAMPLE=... bash polyphase_v2.sh
# Or edit the SAMPLES_FQ map below for batch runs.
#
set -euo pipefail

# ---------------- env / binaries (override-friendly) ----------------
# Resolve python / whatshap from the active env by default. Set PYBIN /
# WHATSHAP explicitly to pin a specific conda env.
PYBIN=${PYBIN:-$(command -v python)}
WHATSHAP=${WHATSHAP:-$(command -v whatshap)}
FREEBAYES=${FREEBAYES:-$(command -v freebayes)}

SCRIPTS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Bundled SpecHLA-derived resources (override with: SPECHLA=/path bash polyphase_v2.sh).
BUNDLED_SPECHLA=${BUNDLED_SPECHLA:-${SCRIPTS_DIR}/resources/spechla}
LEGACY_SPECHLA=${LEGACY_SPECHLA:-"$(cd "${SCRIPTS_DIR}/.." && pwd)/SpecHLA"}
if [[ -z "${SPECHLA:-}" ]]; then
    if [[ -d "$BUNDLED_SPECHLA" ]]; then
        SPECHLA="$BUNDLED_SPECHLA"
    else
        SPECHLA="$LEGACY_SPECHLA"
    fi
fi
SPECHLA_SCRIPT=${SPECHLA_SCRIPT:-${SPECHLA}/script}
SPECHLA_DB=${SPECHLA_DB:-${SPECHLA}/db}
HLA_REF=${HLA_REF:-${SPECHLA_DB}/ref/hla.ref.extend.fa}
DB_PREFIX=${DB_PREFIX:-${SPECHLA_DB}/ref/hla_gen.format.filter.extend.DRB.no26789.v2.fasta}

WORK_DIR=${WORK_DIR:-"$(cd "${SCRIPTS_DIR}/.." && pwd)"}
OUT_ROOT=${OUT_ROOT:-${WORK_DIR}/spechla_out}
ASM_ROOT=${ASM_ROOT:-${WORK_DIR}/asm_v2}
GENE_BED=${GENE_BED:-${SCRIPTS_DIR}/gene.spechla.bed}
ASSEMBLE_PY=${ASSEMBLE_PY:-${SCRIPTS_DIR}/hla_polyphase_assemble.py}

PLOIDY=${PLOIDY:-4}
THREADS=${THREADS:-8}
SAMTOOLS_THREADS=${SAMTOOLS_THREADS:-$THREADS}
SKIP_DONE=${SKIP_DONE:-1}
REUSE_BINNING_ROOT=${REUSE_BINNING_ROOT:-}
REUSE_BINNING_CLEAN_DOWNSTREAM=${REUSE_BINNING_CLEAN_DOWNSTREAM:-0}

# Speed/accuracy tradeoff knobs. Defaults preserve the validated workflow.
BOWTIE2_MODE=${BOWTIE2_MODE:-very-sensitive}
BOWTIE2_K=${BOWTIE2_K:-30}

# Truth-free read-binning rescue. Keep the default conservative: DPB1-only,
# paired-end evidence, and retention-gated so broad all-gene rescue remains
# an explicit validation/diagnostic mode.
READ_BIN_RESCUE=${READ_BIN_RESCUE:-1}
READ_BIN_RESCUE_PY=${READ_BIN_RESCUE_PY:-${SCRIPTS_DIR}/rescue_gene_binned_reads.py}
READ_BIN_RESCUE_GENES=${READ_BIN_RESCUE_GENES:-"HLA-DPB1"}
READ_BIN_RESCUE_BACKGROUND_GENES=${READ_BIN_RESCUE_BACKGROUND_GENES:-"HLA-A HLA-B HLA-C HLA-DRB1 HLA-DPB1 HLA-DQB1"}
READ_BIN_RESCUE_K=${READ_BIN_RESCUE_K:-31}
READ_BIN_RESCUE_MIN_HITS=${READ_BIN_RESCUE_MIN_HITS:-10}
READ_BIN_RESCUE_MIN_MARGIN=${READ_BIN_RESCUE_MIN_MARGIN:-5}
READ_BIN_RESCUE_REQUIRE_BOTH_MATES=${READ_BIN_RESCUE_REQUIRE_BOTH_MATES:-1}
READ_BIN_RESCUE_MIN_MATE_HITS=${READ_BIN_RESCUE_MIN_MATE_HITS:-1}
READ_BIN_RESCUE_MAX_FRAC=${READ_BIN_RESCUE_MAX_FRAC:-0.25}
READ_BIN_RESCUE_MAX_PAIRS=${READ_BIN_RESCUE_MAX_PAIRS:-100000}
READ_BIN_RESCUE_RETENTION_GATE=${READ_BIN_RESCUE_RETENTION_GATE:-1}
READ_BIN_RESCUE_RETENTION_MIN_FULL_PAIRS=${READ_BIN_RESCUE_RETENTION_MIN_FULL_PAIRS:-50}
READ_BIN_RESCUE_RETENTION_MAX_RETAINED_FRAC=${READ_BIN_RESCUE_RETENTION_MAX_RETAINED_FRAC:-0.10}
READ_BIN_RESCUE_RETENTION_MIN_MISSING_FRAC=${READ_BIN_RESCUE_RETENTION_MIN_MISSING_FRAC:-0.30}
READ_BIN_RESCUE_RETENTION_MIN_RESCUE_PAIRS=${READ_BIN_RESCUE_RETENTION_MIN_RESCUE_PAIRS:-50}

# freebayes params: low-donor friendly. Override via env if needed.
FB_MIN_AF=${FB_MIN_AF:-0.03}
FB_MIN_AC=${FB_MIN_AC:-2}
FB_MIN_BQ=${FB_MIN_BQ:-13}
FB_MIN_MQ=${FB_MIN_MQ:-20}
FB_MIN_COV=${FB_MIN_COV:-10}

# GT reassignment under chimerism model
GT_REASSIGN=${GT_REASSIGN:-1}      # set 0 to skip
GT_MIN_DEPTH=${GT_MIN_DEPTH:-10}
GT_DROP_FP_AF=${GT_DROP_FP_AF:-0.05}

# masking
MASK_MIN_DEPTH=${MASK_MIN_DEPTH:-5}
ASSEMBLE_ALIGNER=${ASSEMBLE_ALIGNER:-parasail}
ASSEMBLE_PREFILTER_TOP=${ASSEMBLE_PREFILTER_TOP:-200}

# EM-based iterative remap refinement (post-baseline). 0 = off.
EM_REFINE=${EM_REFINE:-1}
EM_REFINE_MAX_DIFF=${EM_REFINE_MAX_DIFF:-0.7}  # accept EM only if sumAbsDiff < this
EM_REFINE_GENES=${EM_REFINE_GENES:-"HLA-A HLA-B HLA-C HLA-DRB1 HLA-DPB1 HLA-DQB1"}
# Per-gene chi_R re-fitting (Salmon-like) is experimental. In 267016, fixed
# pooled/global chi_R recovers more recipient-minor alleles (21/24) than
# unconstrained per-gene chi (19/24), so keep it off by default.
EM_REFINE_PER_GENE_CHI=${EM_REFINE_PER_GENE_CHI:-0}
EM_REFINE_CHI_PRIOR=${EM_REFINE_CHI_PRIOR:-0.05}
EM_REFINE_TOP_N=${EM_REFINE_TOP_N:-25}
EM_REFINE_MIN_FRAC=${EM_REFINE_MIN_FRAC:-0.001}
EM_REFINE_RECIPIENT_MINOR_RESCUE=${EM_REFINE_RECIPIENT_MINOR_RESCUE:-1}
EM_REFINE_RESCUE_MIN_FRAC=${EM_REFINE_RESCUE_MIN_FRAC:-0.001}
EM_REFINE_RESCUE_MIN_COUNT=${EM_REFINE_RESCUE_MIN_COUNT:-20}
EM_REFINE_RESCUE_MAX_FRAC=${EM_REFINE_RESCUE_MAX_FRAC:-0.08}
EM_REFINE_LOW_RECIPIENT_PRIVATE_RESCUE=${EM_REFINE_LOW_RECIPIENT_PRIVATE_RESCUE:-1}
EM_REFINE_LOW_RECIPIENT_PRIVATE_GENES=${EM_REFINE_LOW_RECIPIENT_PRIVATE_GENES:-HLA-C}
EM_REFINE_LOW_RECIPIENT_PRIVATE_MAX_FRAC=${EM_REFINE_LOW_RECIPIENT_PRIVATE_MAX_FRAC:-0.02}
EM_REFINE_LOW_RECIPIENT_PRIVATE_DOSE_RATIO=${EM_REFINE_LOW_RECIPIENT_PRIVATE_DOSE_RATIO:-0.20}
EM_REFINE_PY=${EM_REFINE_PY:-${SCRIPTS_DIR}/iterative_remap_em.py}
EM_REFINE_GATE_PY=${EM_REFINE_GATE_PY:-${SCRIPTS_DIR}/em_refine_gate.py}
EM_REFINE_GATE_HIGH_MASK=${EM_REFINE_GATE_HIGH_MASK:-0.40}
EM_REFINE_GATE_AMBIG_DIFF=${EM_REFINE_GATE_AMBIG_DIFF:-0.35}
EM_REFINE_GATE_AMBIG_TOP=${EM_REFINE_GATE_AMBIG_TOP:-0.42}
EM_REFINE_GATE_SELECTED_MIN_FRAC=${EM_REFINE_GATE_SELECTED_MIN_FRAC:-0.0005}
EM_REFINE_GATE_BASELINE_MIN_FRAC=${EM_REFINE_GATE_BASELINE_MIN_FRAC:-0.005}
EM_REFINE_GATE_BASELINE_NEAR_TIE=${EM_REFINE_GATE_BASELINE_NEAR_TIE:-0.001}
EM_REFINE_GATE_CLASS_I_BASELINE_GENES=${EM_REFINE_GATE_CLASS_I_BASELINE_GENES:-HLA-A}
EM_REFINE_GATE_CLASS_I_BASELINE_DIFF=${EM_REFINE_GATE_CLASS_I_BASELINE_DIFF:-0.20}
EM_REFINE_GATE_CLASS_I_BASELINE_TOP=${EM_REFINE_GATE_CLASS_I_BASELINE_TOP:-0.50}
# Class-II hap assemblies are often partially masked in these short-read data;
# use stricter EM override gates there so lack of support is not overread as
# absence. A/B/C keep the global 0.7 gate unless explicitly overridden.
EM_REFINE_GATE_GENE_MAX_DIFF=${EM_REFINE_GATE_GENE_MAX_DIFF:-HLA-DRB1=0.35,HLA-DPB1=0.20,HLA-DQB1=0.50}
EM_REFINE_GATE_GENE_HIGH_MASK=${EM_REFINE_GATE_GENE_HIGH_MASK:-}
EM_REFINE_GATE_GENE_AMBIG_DIFF=${EM_REFINE_GATE_GENE_AMBIG_DIFF:-}
EM_REFINE_GATE_GENE_AMBIG_TOP=${EM_REFINE_GATE_GENE_AMBIG_TOP:-}
EM_REFINE_GATE_GENE_SELECTED_MIN_FRAC=${EM_REFINE_GATE_GENE_SELECTED_MIN_FRAC:-}
EM_REFINE_GATE_GENE_BASELINE_MIN_FRAC=${EM_REFINE_GATE_GENE_BASELINE_MIN_FRAC:-}
EM_REFINE_GATE_GENE_BASELINE_NEAR_TIE=${EM_REFINE_GATE_GENE_BASELINE_NEAR_TIE:-}

# Optional exon-level G group fallback/diagnostic for high-mask genes. This
# writes <SAMPLE>.exon_calls.tsv but does not override final calls.
EXON_TYPING=${EXON_TYPING:-1}
EXON_TYPING_GENES=${EXON_TYPING_GENES:-"HLA-DRB1 HLA-DPB1 HLA-DQB1"}
EXON_TYPING_PY=${EXON_TYPING_PY:-${SCRIPTS_DIR}/exon_typing_from_haps.py}

# Optional DRB345 add-on typing. DRB345 is a reporting group for DRB3/DRB4/DRB5,
# constrained by the final DRB1 haplotypes and emitted as an extra diagnostic row
# in <sample>.final_calls.tsv after the main six-gene calls are finalized.
DRB345_TYPING=${DRB345_TYPING:-1}
DRB345_TYPING_PY=${DRB345_TYPING_PY:-${SCRIPTS_DIR}/type_drb345.py}
DRB345_SUBS_PER_2FIELD=${DRB345_SUBS_PER_2FIELD:-5}
DRB345_TOP_PER_LOCUS=${DRB345_TOP_PER_LOCUS:-8}
DRB345_DB_MIN_AS_FRAC=${DRB345_DB_MIN_AS_FRAC:-0.90}
DRB345_DB_MIN_AS=${DRB345_DB_MIN_AS:--100000000}
DRB345_REMAP_MIN_AS_FRAC=${DRB345_REMAP_MIN_AS_FRAC:-0.95}
DRB345_EVIDENCE_K=${DRB345_EVIDENCE_K:-71}
DRB345_MIN_LOCUS_UNIQUE_FRAC=${DRB345_MIN_LOCUS_UNIQUE_FRAC:--1}
DRB345_DRB1_UNTRUSTED_MASK=${DRB345_DRB1_UNTRUSTED_MASK:-0.50}

# Optional post-aggregate constrained direct-likelihood gate. This is off by
# default; when enabled it only searches quartets made from current/baseline
# alleles and applies changes with a high likelihood-gap threshold.
DIRECT_CONSTRAINED_GATE=${DIRECT_CONSTRAINED_GATE:-0}
DIRECT_CONSTRAINED_GATE_GENES=${DIRECT_CONSTRAINED_GATE_GENES:-"HLA-DRB1 HLA-DPB1 HLA-DQB1"}
DIRECT_CONSTRAINED_GATE_GAP=${DIRECT_CONSTRAINED_GATE_GAP:-150}
DIRECT_CONSTRAINED_GATE_TOP_N=${DIRECT_CONSTRAINED_GATE_TOP_N:-0}
DIRECT_CONSTRAINED_GATE_MIN_FRAC=${DIRECT_CONSTRAINED_GATE_MIN_FRAC:-0.002}
DIRECT_CONSTRAINED_GATE_MIN_AS_FRAC=${DIRECT_CONSTRAINED_GATE_MIN_AS_FRAC:-0.95}
DIRECT_CONSTRAINED_GATE_FAMILY_AGG=${DIRECT_CONSTRAINED_GATE_FAMILY_AGG:-max}
DIRECT_CONSTRAINED_BATCH_PY=${DIRECT_CONSTRAINED_BATCH_PY:-${SCRIPTS_DIR}/direct_quartet_likelihood_batch.py}
DIRECT_CONSTRAINED_APPLY_PY=${DIRECT_CONSTRAINED_APPLY_PY:-${SCRIPTS_DIR}/apply_direct_constrained_gate.py}

# Class-II joint rescue. Enabled by default because the rules are guarded by
# high-mask gates and only touch DRB1/DPB1 failure modes validated on the real
# set. Set CLASS2_JOINT_RESCUE=0 to disable. It is truth-free and writes
# backups as calls.class2_joint_input.tsv.
CLASS2_JOINT_RESCUE=${CLASS2_JOINT_RESCUE:-1}
CLASS2_JOINT_RESCUE_PY=${CLASS2_JOINT_RESCUE_PY:-${SCRIPTS_DIR}/apply_class2_joint_rescue.py}
CLASS2_JOINT_DRB1_MIN_MASK=${CLASS2_JOINT_DRB1_MIN_MASK:-0.40}
CLASS2_JOINT_DPB1_MIN_MASK=${CLASS2_JOINT_DPB1_MIN_MASK:-0.40}
CLASS2_JOINT_DPB1_RARE_CUTOFF=${CLASS2_JOINT_DPB1_RARE_CUTOFF:-100}
CLASS2_JOINT_DPB1_MIN_FRAC=${CLASS2_JOINT_DPB1_MIN_FRAC:-0.02}
CLASS2_JOINT_DPB1_TOP_COMMON=${CLASS2_JOINT_DPB1_TOP_COMMON:-6}
CLASS2_JOINT_DPB1_COMMON_MINOR=${CLASS2_JOINT_DPB1_COMMON_MINOR:-1}
CLASS2_JOINT_DPB1_COMMON_MINOR_MAX_NUMBER=${CLASS2_JOINT_DPB1_COMMON_MINOR_MAX_NUMBER:-10}
CLASS2_JOINT_DPB1_COMMON_MINOR_MIN_FRAC=${CLASS2_JOINT_DPB1_COMMON_MINOR_MIN_FRAC:-0.005}
CLASS2_JOINT_DPB1_COMMON_MINOR_MAX_FRAC=${CLASS2_JOINT_DPB1_COMMON_MINOR_MAX_FRAC:-0.09}
CLASS2_JOINT_DPB1_COMMON_MINOR_MIN_WEIGHT=${CLASS2_JOINT_DPB1_COMMON_MINOR_MIN_WEIGHT:-50}
CLASS2_JOINT_DPB1_COMMON_MINOR_TOP=${CLASS2_JOINT_DPB1_COMMON_MINOR_TOP:-12}
CLASS2_JOINT_DPB1_ABSOLUTE_COMMON=${CLASS2_JOINT_DPB1_ABSOLUTE_COMMON:-1}
CLASS2_JOINT_DPB1_ABSOLUTE_COMMON_MAX_NUMBER=${CLASS2_JOINT_DPB1_ABSOLUTE_COMMON_MAX_NUMBER:-10}
CLASS2_JOINT_DPB1_ABSOLUTE_COMMON_MIN_WEIGHT=${CLASS2_JOINT_DPB1_ABSOLUTE_COMMON_MIN_WEIGHT:-150}
CLASS2_JOINT_DPB1_ABSOLUTE_COMMON_MIN_FRAC=${CLASS2_JOINT_DPB1_ABSOLUTE_COMMON_MIN_FRAC:-0.01}
CLASS2_JOINT_DPB1_ABSOLUTE_COMMON_MIN_RATIO=${CLASS2_JOINT_DPB1_ABSOLUTE_COMMON_MIN_RATIO:-2.0}

# chimerism prior: 0 = donor major (allo-HSCT recipient blood, default).
# Set to 1 for solid-organ tx etc.
RECIPIENT_MAJOR=${RECIPIENT_MAJOR:-0}

# ---------------- samples ----------------
# Either set env: FQ1, FQ2, SAMPLE for one sample,
# or fill SAMPLES_FQ below: "sample_id|fq1|fq2"
SAMPLES_FQ=(
    # "267015-HLA-20260415EM_S1_L001|/path/R1.fq.gz|/path/R2.fq.gz"
)
if [[ -n "${FQ1:-}" && -n "${FQ2:-}" && -n "${SAMPLE:-}" ]]; then
    SAMPLES_FQ=("${SAMPLE}|${FQ1}|${FQ2}")
fi
if [[ ${#SAMPLES_FQ[@]} -eq 0 ]]; then
    echo "[ERROR] no samples configured. Set FQ1/FQ2/SAMPLE env or edit SAMPLES_FQ." >&2
    exit 1
fi

mkdir -p "$OUT_ROOT" "$ASM_ROOT"

vcf_ready () {
    local vcf="$1"
    [[ -s "$vcf" && -s "${vcf}.tbi" ]]
}

clean_reuse_downstream_outputs () {
    local SPEC="$1"
    local OUT="$2"
    echo "[step] clean downstream outputs in target run before recomputing"
    rm -f "${OUT}"/A.bam "${OUT}"/A.bam.bai "${OUT}"/B.bam "${OUT}"/B.bam.bai "${OUT}"/C.bam "${OUT}"/C.bam.bai
    rm -f "${OUT}"/DPB1.bam "${OUT}"/DPB1.bam.bai "${OUT}"/DQB1.bam "${OUT}"/DQB1.bam.bai "${OUT}"/DRB1.bam "${OUT}"/DRB1.bam.bai
    rm -f "${OUT}/${SPEC}.merge.bam" "${OUT}/${SPEC}.merge.bam.bai"
    rm -f "${OUT}"/*.freebayes*.vcf.gz "${OUT}"/*.freebayes*.vcf.gz.tbi
    rm -f "${OUT}"/*.phased*.vcf.gz "${OUT}"/*.phased*.vcf.gz.tbi "${OUT}"/*.pooled_continuous.vcf.gz "${OUT}"/*.pooled_continuous.vcf.gz.tbi
    rm -f "${OUT}"/*.chimerism.txt "${OUT}"/*.chi_pooled.txt "${OUT}"/read_bin_rescue_manifest.tsv
    rm -rf "${OUT}"/em_refine "${ASM_ROOT}/${SPEC}"
}

seed_binning_from_reuse_root () {
    local SPEC="$1"
    local OUT="$2"
    [[ -n "$REUSE_BINNING_ROOT" ]] || return 0

    local SRC="${REUSE_BINNING_ROOT}/${SPEC}"
    local MARKER="${OUT}/.reuse_binning_seeded"
    if [[ -f "$MARKER" && $SKIP_DONE -eq 1 \
        && -s "${OUT}/${SPEC}.uniq.R1.fq.gz" \
        && -s "${OUT}/${SPEC}.uniq.R2.fq.gz" \
        && -s "${OUT}/${SPEC}.map_database.bam" \
        && -s "${OUT}/A.R1.fq.gz" ]]; then
        echo "[skip] reused binning cache already seeded from $REUSE_BINNING_ROOT"
        if [[ "$REUSE_BINNING_CLEAN_DOWNSTREAM" == "1" ]]; then
            clean_reuse_downstream_outputs "$SPEC" "$OUT"
            date > "$MARKER"
        fi
        return 0
    fi
    [[ -d "$SRC" ]] || { echo "[ERROR] reuse source missing: $SRC" >&2; exit 1; }
    [[ -s "${SRC}/${SPEC}.uniq.R1.fq.gz" ]] || { echo "[ERROR] source uniq R1 missing: ${SRC}/${SPEC}.uniq.R1.fq.gz" >&2; exit 1; }
    [[ -s "${SRC}/${SPEC}.uniq.R2.fq.gz" ]] || { echo "[ERROR] source uniq R2 missing: ${SRC}/${SPEC}.uniq.R2.fq.gz" >&2; exit 1; }
    [[ -s "${SRC}/${SPEC}.map_database.bam" ]] || { echo "[ERROR] source DB BAM missing: ${SRC}/${SPEC}.map_database.bam" >&2; exit 1; }
    [[ -s "${SRC}/A.R1.fq.gz" ]] || { echo "[ERROR] source per-gene FASTQs missing under $SRC" >&2; exit 1; }

    echo "[step] seed dedup/db-map/gene-bin cache from $SRC"
    mkdir -p "$OUT"
    cp -a "${SRC}/${SPEC}.uniq.R1.fq.gz" "${SRC}/${SPEC}.uniq.R2.fq.gz" "$OUT/"
    cp -a "${SRC}/${SPEC}.map_database.bam" "$OUT/"
    [[ -s "${SRC}/${SPEC}.map_database.bam.bai" ]] && cp -a "${SRC}/${SPEC}.map_database.bam.bai" "$OUT/"
    [[ -s "${SRC}/header.sam" ]] && cp -a "${SRC}/header.sam" "$OUT/"

    local fq
    for fq in "${SRC}"/*.R1.fq.gz "${SRC}"/*.R2.fq.gz; do
        [[ -e "$fq" ]] || continue
        case "$(basename "$fq")" in
            "${SPEC}.uniq.R1.fq.gz"|"${SPEC}.uniq.R2.fq.gz") continue ;;
        esac
        cp -a "$fq" "$OUT/"
    done

    if [[ "$REUSE_BINNING_CLEAN_DOWNSTREAM" == "1" ]]; then
        clean_reuse_downstream_outputs "$SPEC" "$OUT"
    fi
    date > "$MARKER"
}

# Read gene.bed; suffix dup gene names with _2,_3
declare -a CHROMS STARTS ENDS GENES TAGS
declare -A SEEN
while read -r chrom start end gene rest; do
    [[ -z "${chrom:-}" || "${chrom:0:1}" == "#" ]] && continue
    CHROMS+=("$chrom"); STARTS+=("$start"); ENDS+=("$end"); GENES+=("$gene")
    n=${SEEN[$gene]:-0}; n=$((n+1)); SEEN[$gene]=$n
    if [[ $n -gt 1 ]]; then TAGS+=("${gene}_${n}"); else TAGS+=("$gene"); fi
done < "$GENE_BED"

run_one_sample () {
    local SPEC FQ1 FQ2 OUT VCF MERGED_BAM
    SPEC="$1"; FQ1="$2"; FQ2="$3"
    OUT="${OUT_ROOT}/${SPEC}"
    MERGED_BAM="${OUT}/${SPEC}.merge.bam"
    VCF="${OUT}/${SPEC}.freebayes.vcf.gz"
    mkdir -p "$OUT"

    echo "===================================================="
    echo "[$SPEC] starting full pipeline"
    echo "  FQ1=$FQ1"
    echo "  FQ2=$FQ2"
    echo "  OUT=$OUT"
    echo "===================================================="

    seed_binning_from_reuse_root "$SPEC" "$OUT"

    # ---- 0. dedupe read names ----
    local UFQ1="${OUT}/${SPEC}.uniq.R1.fq.gz"
    local UFQ2="${OUT}/${SPEC}.uniq.R2.fq.gz"
    if [[ -f "$UFQ1" && -f "$UFQ2" && $SKIP_DONE -eq 1 ]]; then
        echo "[skip] uniq fastq exist"
    else
        echo "[step] dedupe read names"
        "$PYBIN" "${SPECHLA_SCRIPT}/uniq_read_name.py" "$FQ1" "$UFQ1"
        "$PYBIN" "${SPECHLA_SCRIPT}/uniq_read_name.py" "$FQ2" "$UFQ2"
    fi

    # ---- 1. assign reads to genes (competitive map to IMGT db) ----
    local DB_BAM="${OUT}/${SPEC}.map_database.bam"
    if [[ -f "$DB_BAM" && $SKIP_DONE -eq 1 ]]; then
        echo "[skip] $DB_BAM exists (bowtie2 -> IMGT db)"
    else
        echo "[step] competitive map to IMGT db (bowtie2)"
        bowtie2 "--${BOWTIE2_MODE}" -p "$THREADS" -k "$BOWTIE2_K" \
            -x "$DB_PREFIX" -1 "$UFQ1" -2 "$UFQ2" \
            | samtools view -@ "$SAMTOOLS_THREADS" -bS - \
            | samtools sort -@ "$SAMTOOLS_THREADS" -o "$DB_BAM" -
        samtools index -@ "$SAMTOOLS_THREADS" "$DB_BAM"
    fi
    [[ -f "${DB_BAM}.bai" ]] || samtools index -@ "$SAMTOOLS_THREADS" "$DB_BAM"

    if [[ -f "${OUT}/A.R1.fq.gz" && $SKIP_DONE -eq 1 ]]; then
        echo "[skip] per-gene binned fastqs exist"
    else
        echo "[step] assign_reads_to_genes.py -> per-gene fastqs"
        "$PYBIN" "${SPECHLA_SCRIPT}/assign_reads_to_genes.py" \
            -1 "$UFQ1" -2 "$UFQ2" \
            -n "${SPECHLA}/bin" -o "$OUT" -d 0.1 \
            -b "$DB_BAM" -nm 2
    fi

    if [[ "$READ_BIN_RESCUE" == "1" ]]; then
        local RB_MANIFEST="${OUT}/read_bin_rescue_manifest.tsv"
        if [[ -f "$RB_MANIFEST" && $SKIP_DONE -eq 1 ]]; then
            echo "[skip] read-bin rescue already ran"
        else
            echo "[step] rescue full-FASTQ reads missed by strict gene binning"
            local RB_GENE_ARGS=()
            for gx in $READ_BIN_RESCUE_GENES; do RB_GENE_ARGS+=(--gene "$gx"); done
            local RB_BACKGROUND_GENE_ARGS=()
            for gx in $READ_BIN_RESCUE_BACKGROUND_GENES; do RB_BACKGROUND_GENE_ARGS+=(--background-gene "$gx"); done
            local RB_MATE_ARGS=()
            if [[ "$READ_BIN_RESCUE_REQUIRE_BOTH_MATES" == "1" ]]; then
                RB_MATE_ARGS+=(--require-both-mates)
            fi
            local RB_RETENTION_ARGS=()
            if [[ "$READ_BIN_RESCUE_RETENTION_GATE" == "1" ]]; then
                RB_RETENTION_ARGS+=(
                    --retention-gate
                    --retention-min-full-pairs "$READ_BIN_RESCUE_RETENTION_MIN_FULL_PAIRS"
                    --retention-max-retained-fraction "$READ_BIN_RESCUE_RETENTION_MAX_RETAINED_FRAC"
                    --retention-min-missing-fraction "$READ_BIN_RESCUE_RETENTION_MIN_MISSING_FRAC"
                    --retention-min-rescue-pairs "$READ_BIN_RESCUE_RETENTION_MIN_RESCUE_PAIRS"
                )
            fi
            "$PYBIN" -u "$READ_BIN_RESCUE_PY" \
                --fq1 "$UFQ1" --fq2 "$UFQ2" \
                --fq-dir "$OUT" \
                --exon-dir "${SPECHLA_DB}/HLA/exon" \
                --k "$READ_BIN_RESCUE_K" \
                --min-hits "$READ_BIN_RESCUE_MIN_HITS" \
                --min-margin "$READ_BIN_RESCUE_MIN_MARGIN" \
                --min-mate-hits "$READ_BIN_RESCUE_MIN_MATE_HITS" \
                --max-rescue-fraction "$READ_BIN_RESCUE_MAX_FRAC" \
                --max-rescue-pairs "$READ_BIN_RESCUE_MAX_PAIRS" \
                --manifest "$RB_MANIFEST" \
                "${RB_GENE_ARGS[@]}" \
                "${RB_BACKGROUND_GENE_ARGS[@]}" \
                "${RB_MATE_ARGS[@]}" \
                "${RB_RETENTION_ARGS[@]}"
        fi
    fi

    # ---- 2. per-gene bwa to SpecHLA per-gene refs ----
    local HLAS=(A B C DPB1 DQB1 DRB1)
    local GROUP="@RG\tID:${SPEC}\tSM:${SPEC}"
    if [[ -s "${OUT}/header.sam" && $SKIP_DONE -eq 1 ]]; then
        echo "[skip] ${OUT}/header.sam exists"
    else
        bwa mem -U 10000 -L 10000,10000 -R "$GROUP" \
            "$HLA_REF" "$UFQ1" "$UFQ2" 2>/dev/null \
            | samtools view -H - > "${OUT}/header.sam" || true
    fi
    for hla in "${HLAS[@]}"; do
        local PER_BAM="${OUT}/${hla}.bam"
        if [[ -f "$PER_BAM" && $SKIP_DONE -eq 1 ]]; then
            echo "[skip] $PER_BAM exists"
            continue
        fi
        local R1="${OUT}/${hla}.R1.fq.gz"
        local R2="${OUT}/${hla}.R2.fq.gz"
        if [[ ! -s "$R1" || ! -s "$R2" ]]; then
            echo "[warn] no reads for HLA-${hla}; skipping"
            : > "${PER_BAM}.empty"
            continue
        fi
        local PER_REF="${SPECHLA_DB}/HLA/HLA_${hla}/HLA_${hla}.fa"
        echo "[step] bwa mem HLA-${hla}"
        bwa mem -t "$THREADS" -U 10000 -L 10000,10000 -R "$GROUP" \
            "$PER_REF" "$R1" "$R2" \
            | samtools view -@ "$SAMTOOLS_THREADS" -bS -F 0x800 - \
            | samtools sort -@ "$SAMTOOLS_THREADS" -o "$PER_BAM" -
        samtools index -@ "$SAMTOOLS_THREADS" "$PER_BAM"
    done

    # ---- 3. merge into one BAM on combined HLA ref ----
    if [[ -f "$MERGED_BAM" && $SKIP_DONE -eq 1 ]]; then
        echo "[skip] $MERGED_BAM exists"
    else
        local TO_MERGE=()
        for hla in "${HLAS[@]}"; do
            [[ -f "${OUT}/${hla}.bam" ]] && TO_MERGE+=("${OUT}/${hla}.bam")
        done
        echo "[step] samtools merge -> $MERGED_BAM"
        samtools merge -@ "$SAMTOOLS_THREADS" -f -h "${OUT}/header.sam" "$MERGED_BAM" "${TO_MERGE[@]}"
        samtools index -@ "$SAMTOOLS_THREADS" "$MERGED_BAM"
    fi

    # ---- 4. freebayes -p 4 ----
    if [[ $SKIP_DONE -eq 1 ]] && vcf_ready "$VCF"; then
        echo "[skip] $VCF exists"
    else
        echo "[step] freebayes ploidy=4 (low-MAF / low-AC; haplotype-length 0 to suppress combinatorial FP)"
        "$FREEBAYES" \
            -p "$PLOIDY" \
            --min-alternate-fraction "$FB_MIN_AF" \
            --min-alternate-count "$FB_MIN_AC" \
            --min-base-quality "$FB_MIN_BQ" \
            --min-mapping-quality "$FB_MIN_MQ" \
            --min-coverage "$FB_MIN_COV" \
            --haplotype-length 0 \
            --use-best-n-alleles 4 \
            -f "$HLA_REF" "$MERGED_BAM" \
        | bcftools norm -f "$HLA_REF" -a -m -any -Oz -o "$VCF"
        tabix -f -p vcf "$VCF"
    fi

    # ---- 5. estimate chimerism (whole VCF) ----
    local CHIM_LOG="${OUT}/${SPEC}.chimerism.txt"
    if [[ ! -f "$CHIM_LOG" || $SKIP_DONE -eq 0 ]]; then
        echo "[step] estimate chimerism"
        PYTHONPATH="${SCRIPTS_DIR}:${PYTHONPATH:-}" "$PYBIN" - "$VCF" "$RECIPIENT_MAJOR" > "$CHIM_LOG" <<'PYEOF'
import sys, pysam
from hla_polyphase_assemble import estimate_chimerism_from_vcf
vcf_path, recip_major = sys.argv[1], sys.argv[2] == "1"
v = pysam.VariantFile(vcf_path)
all_obs = [[]]
for rec in v:
    s = rec.samples[0]
    gt = s.get('GT')
    if gt is None or len(gt) != 4 or any(a is None for a in gt):
        continue
    if not rec.alts or len(rec.alts) != 1 or not (set(gt) <= {0, 1}):
        continue
    ad = s.get('AD')
    if not ad or len(ad) < 2 or any(x is None for x in ad) or sum(ad) < 10:
        continue
    all_obs[0].append((rec.pos - 1, tuple(gt), ad[1] / sum(ad)))
est, n, info = estimate_chimerism_from_vcf(all_obs)
if est is None:
    print(f"FAIL n={n}")
    sys.exit(0)
chi_r = est if recip_major else (1.0 - est)
print(f"major={est:.4f} n={n} info={info} recipient_major={recip_major} chi_R={chi_r:.4f}")
PYEOF
    fi
    cat "$CHIM_LOG"
    local CHI_R
    CHI_R=$(awk '/chi_R=/{for(i=1;i<=NF;i++)if($i~/^chi_R=/){split($i,a,"=");print a[2]}}' "$CHIM_LOG")

    # ---- 5b. pooled-continuous freebayes -> AF-based chi_R (Method 1) ----
    # Bypasses ploidy=4 GT decisions which inflate chi_R when true chi is small
    # (low-AF recipient sites get called 0/0/0/0 and dropped). Falls back to
    # the GT-based estimate if anything fails.
    if [[ "${USE_POOLED_CHI:-1}" -eq 1 ]]; then
        local PC_VCF="${OUT}/${SPEC}.pooled_continuous.vcf.gz"
        local PC_LOG="${OUT}/${SPEC}.chi_pooled.txt"
        if [[ $SKIP_DONE -eq 0 ]] || ! vcf_ready "$PC_VCF"; then
            echo "[step] freebayes --pooled-continuous (chi_R estimation only)"
            "$FREEBAYES" \
                --pooled-continuous \
                --min-alternate-fraction 0.005 \
                --min-alternate-count 2 \
                --min-base-quality "$FB_MIN_BQ" \
                --min-mapping-quality "$FB_MIN_MQ" \
                --min-coverage 30 \
                --haplotype-length 0 \
                --use-best-n-alleles 2 \
                -f "$HLA_REF" "$MERGED_BAM" \
            | bcftools norm -f "$HLA_REF" -a -m -any -Oz -o "$PC_VCF" \
            && tabix -f -p vcf "$PC_VCF" || echo "[warn] pooled-continuous freebayes failed"
        fi
        if [[ -f "$PC_VCF" ]]; then
            "$PYBIN" "${SCRIPTS_DIR}/estimate_chi_pooled.py" "$PC_VCF" > "$PC_LOG" 2>&1 || true
            cat "$PC_LOG"
            local CHI_R_PC
            CHI_R_PC=$(awk '/^GLOBAL[[:space:]]+chi_R=/{for(i=1;i<=NF;i++)if($i~/^chi_R=/){split($i,a,"=");print a[2];exit}}' "$PC_LOG")
            if [[ -n "$CHI_R_PC" ]] && awk -v x="$CHI_R_PC" 'BEGIN{exit !(x>0 && x<0.5)}'; then
                echo "[chi] pooled-continuous chi_R=$CHI_R_PC (overrides GT-based $CHI_R)"
                CHI_R="$CHI_R_PC"
            else
                echo "[chi] pooled-continuous chi_R unusable ('$CHI_R_PC'); keeping GT-based $CHI_R"
            fi
        fi
    fi

    # ---- 6. per-gene polyphase + assemble + mask ----
    for i in "${!TAGS[@]}"; do
        local TAG="${TAGS[$i]}" GENE="${GENES[$i]}"
        local CHROM="${CHROMS[$i]}" START="${STARTS[$i]}" END="${ENDS[$i]}"
        local TAG_LC; TAG_LC=$(echo "$TAG" | tr '[:upper:]' '[:lower:]')
        local REGION="${CHROM}:$((START + 1))-${END}"

        local FB_GENE_VCF="${OUT}/${SPEC}.freebayes.${TAG_LC}.vcf.gz"
        local PH_VCF="${OUT}/${SPEC}.phased.${TAG_LC}.vcf.gz"
        local ASM_OUT="${ASM_ROOT}/${SPEC}/${TAG_LC}"

        echo "==== [$SPEC] $TAG ($REGION) ===="

        if [[ -f "${ASM_OUT}/${GENE}/calls.tsv" && $SKIP_DONE -eq 1 ]]; then
            echo "[skip] typing already done: ${ASM_OUT}/${GENE}/calls.tsv"
            continue
        fi
        mkdir -p "$ASM_OUT"

        if [[ $SKIP_DONE -eq 1 ]] && vcf_ready "$FB_GENE_VCF"; then
            echo "[skip] $FB_GENE_VCF exists"
        else
            echo "[step] slice VCF -> $FB_GENE_VCF"
            bcftools view -r "$REGION" -Oz -o "$FB_GENE_VCF" "$VCF"
            tabix -f -p vcf "$FB_GENE_VCF"
        fi

        # GT reassignment under chimerism model (before whatshap)
        local IN_FOR_PHASE="$FB_GENE_VCF"
        if [[ $GT_REASSIGN -eq 1 && -n "${CHI_R:-}" ]]; then
            local FB_REGT_VCF="${OUT}/${SPEC}.freebayes_regt.${TAG_LC}.vcf.gz"
            if [[ $SKIP_DONE -eq 1 ]] && vcf_ready "$FB_REGT_VCF"; then
                echo "[skip] $FB_REGT_VCF exists"
            else
                echo "[step] reassign GT under chi_R=$CHI_R"
                "$PYBIN" "${SCRIPTS_DIR}/reassign_gt_chimeric.py" \
                    --vcf "$FB_GENE_VCF" \
                    --out "$FB_REGT_VCF" \
                    --chi-r "$CHI_R" \
                    --min-depth "$GT_MIN_DEPTH" \
                    --drop-fp-af "$GT_DROP_FP_AF"
            fi
            IN_FOR_PHASE="$FB_REGT_VCF"
        fi

        if [[ $SKIP_DONE -eq 1 ]] && vcf_ready "$PH_VCF"; then
            echo "[skip] $PH_VCF exists"
        else
            echo "[step] whatshap polyphase -> $PH_VCF"
            "$WHATSHAP" polyphase \
                --ploidy "$PLOIDY" \
                --reference "$HLA_REF" \
                --threads "$THREADS" \
                --ignore-read-groups \
                --output "$PH_VCF" \
                "$IN_FOR_PHASE" "$MERGED_BAM"
            tabix -f -p vcf "$PH_VCF"
        fi

        local CHIM_ARGS=()
        if [[ -n "${CHI_R:-}" ]]; then
            CHIM_ARGS=(--chimerism "$CHI_R")
        else
            CHIM_ARGS=(--chimerism auto)
            [[ $RECIPIENT_MAJOR -eq 1 ]] && CHIM_ARGS+=(--recipient-major)
        fi

        echo "[step] typing -> $ASM_OUT"
        "$PYBIN" -u "$ASSEMBLE_PY" \
            --vcf "$PH_VCF" \
            --ref "$HLA_REF" \
            --gene-bed "$GENE_BED" \
            --genes "$GENE" \
            --out "$ASM_OUT" \
            --imgt "$DB_PREFIX" \
            --paired-diploids \
            --bam "$MERGED_BAM" \
            --mask-min-depth "$MASK_MIN_DEPTH" \
            --aligner "$ASSEMBLE_ALIGNER" \
            --prefilter-top "$ASSEMBLE_PREFILTER_TOP" \
            "${CHIM_ARGS[@]}" \
            --dump-block-fa
    done

    echo "[$SPEC] done."
}

run_em_refine () {
    # Run iterative_remap_em.py for this sample, then per gene override the
    # baseline calls.tsv when sumAbsDiff < EM_REFINE_MAX_DIFF. Saves the
    # original as calls.baseline.tsv. No-op if EM_REFINE != 1.
    local SPEC="$1"
    [[ $EM_REFINE -ne 1 ]] && return 0
    local OUT="${OUT_ROOT}/${SPEC}"
    local CHIM_LOG="${OUT}/${SPEC}.chimerism.txt"
    local PC_LOG="${OUT}/${SPEC}.chi_pooled.txt"
    local CHI_R
    CHI_R=$(awk '/chi_R=/{for(i=1;i<=NF;i++)if($i~/^chi_R=/){split($i,a,"=");print a[2]}}' "$CHIM_LOG" 2>/dev/null || echo "")
    if [[ "${USE_POOLED_CHI:-1}" -eq 1 && -f "$PC_LOG" ]]; then
        local CHI_R_PC
        CHI_R_PC=$(awk '/^GLOBAL[[:space:]]+chi_R=/{for(i=1;i<=NF;i++)if($i~/^chi_R=/){split($i,a,"=");print a[2];exit}}' "$PC_LOG" 2>/dev/null || echo "")
        if [[ -n "$CHI_R_PC" ]] && awk -v x="$CHI_R_PC" 'BEGIN{exit !(x>0 && x<0.5)}'; then
            echo "[em-refine] using pooled-continuous chi_R=$CHI_R_PC (was $CHI_R)"
            CHI_R="$CHI_R_PC"
        fi
    fi
    if [[ -z "$CHI_R" ]]; then
        echo "[em-refine] no chi_R found; skipping"
        return 0
    fi
    local EM_OUT="${OUT}/em_refine"
    mkdir -p "$EM_OUT"
    local GENE_ARGS=()
    for gx in $EM_REFINE_GENES; do GENE_ARGS+=(--gene "$gx"); done
    local CHI_ARGS=()
    if [[ "$EM_REFINE_PER_GENE_CHI" == "1" ]]; then
        CHI_ARGS+=(--per-gene-chi --chi-prior "$EM_REFINE_CHI_PRIOR")
    fi
    local RESCUE_ARGS=()
    if [[ "$EM_REFINE_RECIPIENT_MINOR_RESCUE" == "1" ]]; then
        RESCUE_ARGS+=(--recipient-minor-rescue
            --rescue-min-frac "$EM_REFINE_RESCUE_MIN_FRAC"
            --rescue-min-count "$EM_REFINE_RESCUE_MIN_COUNT"
            --rescue-max-frac "$EM_REFINE_RESCUE_MAX_FRAC")
    fi
    if [[ "$EM_REFINE_LOW_RECIPIENT_PRIVATE_RESCUE" == "1" ]]; then
        RESCUE_ARGS+=(--low-recipient-private-rescue
            --low-recipient-private-genes "$EM_REFINE_LOW_RECIPIENT_PRIVATE_GENES"
            --low-recipient-private-max-frac "$EM_REFINE_LOW_RECIPIENT_PRIVATE_MAX_FRAC"
            --low-recipient-private-dose-ratio "$EM_REFINE_LOW_RECIPIENT_PRIVATE_DOSE_RATIO")
    fi
    echo "[step] EM iterative remap (chi_R=$CHI_R, max-diff=$EM_REFINE_MAX_DIFF)"
    "$PYBIN" -u "$EM_REFINE_PY" \
        --sample "$SPEC" \
        --fq-dir "$OUT" \
        --chi-r "$CHI_R" \
        --imgt "$DB_PREFIX" \
        --out-dir "$EM_OUT" \
        --baseline-root "${ASM_ROOT}/${SPEC}" \
        --threads "$THREADS" \
        "${CHI_ARGS[@]}" \
        "${RESCUE_ARGS[@]}" \
        ${EM_REFINE_TOP_N:+--top-n "$EM_REFINE_TOP_N"} \
        ${EM_REFINE_MIN_FRAC:+--min-frac "$EM_REFINE_MIN_FRAC"} \
        "${GENE_ARGS[@]}" 2>&1 | tee "${EM_OUT}/em_refine.log" || {
            echo "[em-refine] python failed; leaving baseline calls intact"; return 0
        }
    # per-gene override
    for gx in $EM_REFINE_GENES; do
        local TAG="$gx"
        local TAG_LC; TAG_LC=$(echo "$TAG" | tr '[:upper:]' '[:lower:]')
        local SUMM="${EM_OUT}/${gx}.summary.tsv"
        local EM_CALLS="${EM_OUT}/${gx}.calls.tsv"
        local TF_COUNTS="${EM_OUT}/${gx}.tf_counts.tsv"
        local DST="${ASM_ROOT}/${SPEC}/${TAG_LC}/${TAG}/calls.tsv"
        if [[ ! -f "$SUMM" || ! -f "$EM_CALLS" || ! -f "$DST" ]]; then
            echo "[em-refine] $gx: missing inputs (em or baseline); skip"
            continue
        fi
        local GATE GATE_ACTION GATE_REASON
        GATE=$("$PYBIN" "$EM_REFINE_GATE_PY" \
            --gene "$gx" \
            --gene-dir "${ASM_ROOT}/${SPEC}/${TAG_LC}/${TAG}" \
            --summary "$SUMM" \
            --em-calls "$EM_CALLS" \
            --baseline-calls "$DST" \
            --tf-counts "$TF_COUNTS" \
            --max-diff "$EM_REFINE_MAX_DIFF" \
            --high-mask "$EM_REFINE_GATE_HIGH_MASK" \
            --ambiguous-diff "$EM_REFINE_GATE_AMBIG_DIFF" \
            --ambiguous-top-frac "$EM_REFINE_GATE_AMBIG_TOP" \
            --selected-min-frac "$EM_REFINE_GATE_SELECTED_MIN_FRAC" \
            --baseline-min-frac "$EM_REFINE_GATE_BASELINE_MIN_FRAC" \
            --baseline-near-tie "$EM_REFINE_GATE_BASELINE_NEAR_TIE" \
            --class-i-baseline-genes "$EM_REFINE_GATE_CLASS_I_BASELINE_GENES" \
            --class-i-baseline-diff "$EM_REFINE_GATE_CLASS_I_BASELINE_DIFF" \
            --class-i-baseline-top-frac "$EM_REFINE_GATE_CLASS_I_BASELINE_TOP" \
            --gene-max-diff "$EM_REFINE_GATE_GENE_MAX_DIFF" \
            --gene-high-mask "$EM_REFINE_GATE_GENE_HIGH_MASK" \
            --gene-ambiguous-diff "$EM_REFINE_GATE_GENE_AMBIG_DIFF" \
            --gene-ambiguous-top-frac "$EM_REFINE_GATE_GENE_AMBIG_TOP" \
            --gene-selected-min-frac "$EM_REFINE_GATE_GENE_SELECTED_MIN_FRAC" \
            --gene-baseline-min-frac "$EM_REFINE_GATE_GENE_BASELINE_MIN_FRAC" \
            --gene-baseline-near-tie "$EM_REFINE_GATE_GENE_BASELINE_NEAR_TIE")
        GATE_ACTION=${GATE%%$'\t'*}
        GATE_REASON=${GATE#*$'\t'}
        if [[ "$GATE_ACTION" == "OVERRIDE" ]]; then
            if [[ ! -f "${DST%.tsv}.baseline.tsv" ]]; then
                cp "$DST" "${DST%.tsv}.baseline.tsv"
            fi
            cp "$EM_CALLS" "${DST}.tmp.$$"
            mv "${DST}.tmp.$$" "$DST"
            echo "[em-refine] $gx: $GATE_REASON -> OVERRIDE (baseline saved)"
        else
            echo "[em-refine] $gx: $GATE_REASON -> KEEP baseline"
        fi
    done
}

run_direct_constrained_gate () {
    local SPEC="$1"
    [[ "$DIRECT_CONSTRAINED_GATE" == "1" ]] || return 0
    if [[ "$EM_REFINE" != "1" ]]; then
        echo "[direct-gate] EM_REFINE=0; constrained direct gate needs em_refine outputs, skip"
        return 0
    fi

    local OUT="${OUT_ROOT}/${SPEC}"
    local FINAL="${ASM_ROOT}/${SPEC}/${SPEC}.final_calls.tsv"
    if [[ ! -f "$FINAL" ]]; then
        echo "[direct-gate] missing final calls before gate: $FINAL; skip"
        return 0
    fi

    local DIRECT_TSV="${OUT}/direct_constrained_gate.tsv"
    local DIRECT_FAIL="${OUT}/direct_constrained_gate.fail.tsv"
    local DIRECT_MANIFEST="${OUT}/direct_constrained_gate_manifest.tsv"
    local SAM_CACHE="${OUT}/direct_constrained_sam_cache"
    local GENE_ARGS=()
    for gx in $DIRECT_CONSTRAINED_GATE_GENES; do GENE_ARGS+=(--gene "$gx"); done

    if [[ $SKIP_DONE -eq 1 && -s "$DIRECT_TSV" ]]; then
        echo "[direct-gate] skip direct likelihood; $DIRECT_TSV exists"
    else
        echo "[step] constrained direct-likelihood gate scan (gap>=$DIRECT_CONSTRAINED_GATE_GAP)"
        "$PYBIN" -u "$DIRECT_CONSTRAINED_BATCH_PY" \
            --sample "$SPEC" \
            --spechla-root "$OUT_ROOT" \
            --asm-root "$ASM_ROOT" \
            --out-tsv "$DIRECT_TSV" \
            --fail-tsv "$DIRECT_FAIL" \
            --sam-cache-dir "$SAM_CACHE" \
            --threads "$THREADS" \
            --min-as-frac "$DIRECT_CONSTRAINED_GATE_MIN_AS_FRAC" \
            --direct-top-n "$DIRECT_CONSTRAINED_GATE_TOP_N" \
            --direct-min-frac "$DIRECT_CONSTRAINED_GATE_MIN_FRAC" \
            --direct-family-agg "$DIRECT_CONSTRAINED_GATE_FAMILY_AGG" \
            "${GENE_ARGS[@]}" || {
                echo "[direct-gate] likelihood scan failed; leaving current calls intact"
                return 0
            }
    fi

    echo "[step] apply constrained direct-likelihood gate"
    "$PYBIN" -u "$DIRECT_CONSTRAINED_APPLY_PY" \
        --direct-tsv "$DIRECT_TSV" \
        --in-asm-root "$ASM_ROOT" \
        --in-place \
        --spechla-root "$OUT_ROOT" \
        --g-group "${SPECHLA_DB}/HLA/hla_nom_g.txt" \
        --gap-threshold "$DIRECT_CONSTRAINED_GATE_GAP" \
        --manifest "$DIRECT_MANIFEST" || {
            echo "[direct-gate] apply failed; leaving current calls intact"
            return 0
        }
}

run_class2_joint_rescue () {
    local SPEC="$1"
    [[ "$CLASS2_JOINT_RESCUE" == "1" ]] || return 0
    if [[ "$EM_REFINE" != "1" ]]; then
        echo "[class2-joint] EM_REFINE=0; class-II joint rescue needs em_refine outputs, skip"
        return 0
    fi

    local OUT="${OUT_ROOT}/${SPEC}"
    local FINAL="${ASM_ROOT}/${SPEC}/${SPEC}.final_calls.tsv"
    if [[ ! -f "$FINAL" ]]; then
        echo "[class2-joint] missing final calls before rescue: $FINAL; skip"
        return 0
    fi

    local MANIFEST="${OUT}/class2_joint_rescue_manifest.tsv"
    local COMMON_MINOR_ARGS=()
    if [[ "$CLASS2_JOINT_DPB1_COMMON_MINOR" == "1" ]]; then
        COMMON_MINOR_ARGS+=(--dpb1-common-minor
            --dpb1-common-minor-max-number "$CLASS2_JOINT_DPB1_COMMON_MINOR_MAX_NUMBER"
            --dpb1-common-minor-min-fraction "$CLASS2_JOINT_DPB1_COMMON_MINOR_MIN_FRAC"
            --dpb1-common-minor-max-fraction "$CLASS2_JOINT_DPB1_COMMON_MINOR_MAX_FRAC"
            --dpb1-common-minor-min-weight "$CLASS2_JOINT_DPB1_COMMON_MINOR_MIN_WEIGHT"
            --dpb1-common-minor-top "$CLASS2_JOINT_DPB1_COMMON_MINOR_TOP")
    fi
    local ABSOLUTE_COMMON_ARGS=()
    if [[ "$CLASS2_JOINT_DPB1_ABSOLUTE_COMMON" == "1" ]]; then
        ABSOLUTE_COMMON_ARGS+=(--dpb1-absolute-common
            --dpb1-absolute-common-max-number "$CLASS2_JOINT_DPB1_ABSOLUTE_COMMON_MAX_NUMBER"
            --dpb1-absolute-common-min-weight "$CLASS2_JOINT_DPB1_ABSOLUTE_COMMON_MIN_WEIGHT"
            --dpb1-absolute-common-min-fraction "$CLASS2_JOINT_DPB1_ABSOLUTE_COMMON_MIN_FRAC"
            --dpb1-absolute-common-min-ratio "$CLASS2_JOINT_DPB1_ABSOLUTE_COMMON_MIN_RATIO")
    fi
    echo "[step] class-II joint rescue (DRB1-DQB1 LD + DPB1 rare/common-minor rescue)"
    "$PYBIN" -u "$CLASS2_JOINT_RESCUE_PY" \
        --sample "$SPEC" \
        --in-asm-root "$ASM_ROOT" \
        --in-place \
        --spechla-root "$OUT_ROOT" \
        --g-group "${SPECHLA_DB}/HLA/hla_nom_g.txt" \
        --manifest "$MANIFEST" \
        --drb1-min-mask "$CLASS2_JOINT_DRB1_MIN_MASK" \
        --dpb1-min-mask "$CLASS2_JOINT_DPB1_MIN_MASK" \
        --dpb1-rare-cutoff "$CLASS2_JOINT_DPB1_RARE_CUTOFF" \
        --dpb1-min-fraction "$CLASS2_JOINT_DPB1_MIN_FRAC" \
        --dpb1-top-common "$CLASS2_JOINT_DPB1_TOP_COMMON" \
        "${COMMON_MINOR_ARGS[@]}" \
        "${ABSOLUTE_COMMON_ARGS[@]}" || {
            echo "[class2-joint] rescue failed; leaving current calls intact"
            return 0
        }
}

run_drb345_typing () {
    local SPEC="$1"
    [[ "$DRB345_TYPING" == "1" ]] || return 0
    local OUT="${OUT_ROOT}/${SPEC}"
    local FINAL="${ASM_ROOT}/${SPEC}/${SPEC}.final_calls.tsv"
    local DB_BAM="${OUT}/${SPEC}.map_database.bam"
    if [[ ! -s "$FINAL" ]]; then
        echo "[drb345] missing final calls: $FINAL; skip"
        return 0
    fi
    if [[ ! -s "$DB_BAM" ]]; then
        echo "[drb345] missing DB BAM: $DB_BAM; skip"
        return 0
    fi
    if [[ -f "${OUT}/drb345/HLA-DRB345.manifest.tsv" && $SKIP_DONE -eq 1 ]] \
        && awk -F'\t' 'NR>1 && $2=="HLA-DRB345" {found=1} END{exit !found}' "$FINAL"; then
        echo "[skip] DRB345 typing already ran and final row exists"
        return 0
    fi
    echo "[step] DRB345 linked add-on typing"
    "$PYBIN" -u "$DRB345_TYPING_PY" \
        --sample "$SPEC" \
        --fq-dir "$OUT" \
        --db-bam "$DB_BAM" \
        --asm-root "$ASM_ROOT" \
        --final-calls "$FINAL" \
        --imgt "$DB_PREFIX" \
        --g-group "${SPECHLA_DB}/HLA/hla_nom_g.txt" \
        --threads "$THREADS" \
        --subs-per-2field "$DRB345_SUBS_PER_2FIELD" \
        --top-per-locus "$DRB345_TOP_PER_LOCUS" \
        --db-min-as-frac "$DRB345_DB_MIN_AS_FRAC" \
        --db-min-as "$DRB345_DB_MIN_AS" \
        --remap-min-as-frac "$DRB345_REMAP_MIN_AS_FRAC" \
        --evidence-k "$DRB345_EVIDENCE_K" \
        --min-locus-unique-frac "$DRB345_MIN_LOCUS_UNIQUE_FRAC" \
        --drb1-untrusted-mask "$DRB345_DRB1_UNTRUSTED_MASK" \
        || echo "[drb345] typing failed; continuing with main calls"
}

for entry in "${SAMPLES_FQ[@]}"; do
    IFS='|' read -r S F1 F2 <<<"$entry"
    run_one_sample "$S" "$F1" "$F2"
    run_em_refine "$S"

    if [[ "$EXON_TYPING" == "1" ]]; then
        if command -v blastn >/dev/null 2>&1; then
            EXON_OUT="${ASM_ROOT}/${S}/${S}.exon_calls.tsv"
            EXON_ARGS=(--genes)
            for gx in $EXON_TYPING_GENES; do EXON_ARGS+=("$gx"); done
            "$PYBIN" "${EXON_TYPING_PY}" \
                --asm-root "$ASM_ROOT" --sample "$S" --out "$EXON_OUT" \
                --g-group "${SPECHLA_DB}/HLA/hla_nom_g.txt" \
                "${EXON_ARGS[@]}" \
                || echo "[exon-typing] failed; continuing with primary calls"
        else
            echo "[exon-typing] blastn not found; skip exon fallback diagnostics"
        fi
    fi

    # ---- 8. aggregate per-gene calls.tsv into one final summary ----
    FINAL="${ASM_ROOT}/${S}/${S}.final_calls.tsv"
    "$PYBIN" "${SCRIPTS_DIR}/aggregate_calls.py" \
        --asm-root "$ASM_ROOT" --sample "$S" --out "$FINAL" \
        --spechla-root "$OUT_ROOT" \
        --g-group "${SPECHLA_DB}/HLA/hla_nom_g.txt" \
        && echo "[FINAL] ${S}: ${FINAL}"

    run_direct_constrained_gate "$S"
    run_class2_joint_rescue "$S"
    run_drb345_typing "$S"
done

echo "[INFO] All samples processed."
