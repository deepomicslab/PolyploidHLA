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

# Speed/accuracy tradeoff knobs. Defaults preserve the validated workflow.
BOWTIE2_MODE=${BOWTIE2_MODE:-very-sensitive}
BOWTIE2_K=${BOWTIE2_K:-30}

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
EM_REFINE_PY=${EM_REFINE_PY:-${SCRIPTS_DIR}/iterative_remap_em.py}

# Optional exon-level G group fallback/diagnostic for high-mask genes. This
# writes <SAMPLE>.exon_calls.tsv but does not override final calls.
EXON_TYPING=${EXON_TYPING:-1}
EXON_TYPING_GENES=${EXON_TYPING_GENES:-"HLA-DRB1 HLA-DPB1 HLA-DQB1"}
EXON_TYPING_PY=${EXON_TYPING_PY:-${SCRIPTS_DIR}/exon_typing_from_haps.py}

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

    # ---- 2. per-gene bwa to SpecHLA per-gene refs ----
    local HLAS=(A B C DPB1 DQB1 DRB1)
    local GROUP="@RG\tID:${SPEC}\tSM:${SPEC}"
    bwa mem -U 10000 -L 10000,10000 -R "$GROUP" \
        "$HLA_REF" "$UFQ1" "$UFQ2" 2>/dev/null \
        | samtools view -H - > "${OUT}/header.sam" || true
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
    if [[ -f "$VCF" && $SKIP_DONE -eq 1 ]]; then
        echo "[skip] $VCF exists"
    else
        echo "[step] freebayes ploidy=4 (low-MAF / low-AC; haplotype-length 0 to suppress combinatorial FP)"
        freebayes \
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
        if [[ ! -f "$PC_VCF" || $SKIP_DONE -eq 0 ]]; then
            echo "[step] freebayes --pooled-continuous (chi_R estimation only)"
            freebayes \
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
            CHI_R_PC=$(awk '/^GLOBAL chi_R=/{split($2,a,"="); print a[2]; exit}' "$PC_LOG")
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

        if [[ -f "$FB_GENE_VCF" && $SKIP_DONE -eq 1 ]]; then
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
            if [[ -f "$FB_REGT_VCF" && $SKIP_DONE -eq 1 ]]; then
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

        if [[ -f "$PH_VCF" && $SKIP_DONE -eq 1 ]]; then
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
        CHI_R_PC=$(awk '/^GLOBAL chi_R=/{split($2,a,"="); print a[2]; exit}' "$PC_LOG" 2>/dev/null || echo "")
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
    echo "[step] EM iterative remap (chi_R=$CHI_R, max-diff=$EM_REFINE_MAX_DIFF)"
    "$PYBIN" -u "$EM_REFINE_PY" \
        --sample "$SPEC" \
        --fq-dir "$OUT" \
        --chi-r "$CHI_R" \
        --imgt "$DB_PREFIX" \
        --out-dir "$EM_OUT" \
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
        local DST="${ASM_ROOT}/${SPEC}/${TAG_LC}/${TAG}/calls.tsv"
        if [[ ! -f "$SUMM" || ! -f "$EM_CALLS" || ! -f "$DST" ]]; then
            echo "[em-refine] $gx: missing inputs (em or baseline); skip"
            continue
        fi
        local DIFF
        DIFF=$(awk 'NR==2{print $2}' "$SUMM")
        local KEEP
        KEEP=$(awk -v d="$DIFF" -v t="$EM_REFINE_MAX_DIFF" 'BEGIN{print (d+0 < t+0) ? 1 : 0}')
        if [[ "$KEEP" == "1" ]]; then
            if [[ ! -f "${DST%.tsv}.baseline.tsv" ]]; then
                cp "$DST" "${DST%.tsv}.baseline.tsv"
            fi
            cp "$EM_CALLS" "$DST"
            echo "[em-refine] $gx: sumAbsDiff=$DIFF < $EM_REFINE_MAX_DIFF -> OVERRIDE (baseline saved)"
        else
            echo "[em-refine] $gx: sumAbsDiff=$DIFF >= $EM_REFINE_MAX_DIFF -> KEEP baseline"
        fi
    done
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
        --g-group "${SPECHLA_DB}/HLA/hla_nom_g.txt" \
        && echo "[FINAL] ${S}: ${FINAL}"
done

echo "[INFO] All samples processed."
