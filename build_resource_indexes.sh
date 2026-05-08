#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'USAGE'
Usage: bash scripts/build_resource_indexes.sh [--resources PATH] [--force]

Build or repair indexes for the bundled HLA resource directory.

Options:
  --resources PATH   Resource root laid out like scripts/resources/spechla
                     (default: scripts/resources/spechla)
  --force            Rebuild indexes even if output files already exist
  -h, --help         Show this help message

Required tools on PATH: samtools, bwa, bowtie2-build.
Optional tools: makeblastdb for exon/per-gene BLAST databases.
USAGE
}

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
RESOURCE_ROOT="${SCRIPT_DIR}/resources/spechla"
FORCE=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --resources)
            [[ $# -ge 2 ]] || { echo "[ERROR] --resources needs a path" >&2; exit 2; }
            RESOURCE_ROOT="$2"
            shift 2
            ;;
        --force)
            FORCE=1
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "[ERROR] Unknown argument: $1" >&2
            usage >&2
            exit 2
            ;;
    esac
done

[[ -d "$RESOURCE_ROOT" ]] || { echo "[ERROR] Resource root not found: $RESOURCE_ROOT" >&2; exit 1; }
RESOURCE_ROOT=$(cd "$RESOURCE_ROOT" && pwd)
REF_DIR="${RESOURCE_ROOT}/db/ref"
HLA_DIR="${RESOURCE_ROOT}/db/HLA"
EXON_DIR="${HLA_DIR}/exon"
COMBINED_REF="${REF_DIR}/hla.ref.extend.fa"
IMGT_FASTA="${REF_DIR}/hla_gen.format.filter.extend.DRB.no26789.v2.fasta"
GENES=(HLA_A HLA_B HLA_C HLA_DPB1 HLA_DQB1 HLA_DRB1)

log() {
    echo "[INFO] $*"
}

warn() {
    echo "[WARN] $*" >&2
}

need_tool() {
    command -v "$1" >/dev/null 2>&1 || { echo "[ERROR] Missing required tool: $1" >&2; exit 1; }
}

has_all() {
    local file
    for file in "$@"; do
        [[ -s "$file" ]] || return 1
    done
    return 0
}

run_if_missing() {
    local label="$1"
    local outputs_csv="$2"
    shift 2
    IFS=',' read -r -a outputs <<< "$outputs_csv"
    if [[ "$FORCE" -eq 0 ]] && has_all "${outputs[@]}"; then
        log "skip ${label}"
        return 0
    fi
    log "build ${label}"
    "$@"
}

build_faidx() {
    local fasta="$1"
    [[ -s "$fasta" ]] || { warn "missing FASTA: $fasta"; return 0; }
    run_if_missing "faidx $(basename "$fasta")" "${fasta}.fai" samtools faidx "$fasta"
}

build_bwa() {
    local fasta="$1"
    [[ -s "$fasta" ]] || { warn "missing FASTA: $fasta"; return 0; }
    run_if_missing "bwa index $(basename "$fasta")" "${fasta}.amb,${fasta}.ann,${fasta}.bwt,${fasta}.pac,${fasta}.sa" bwa index "$fasta"
}

build_dict() {
    local fasta="$1"
    local dict_out="${fasta%.fa}.dict"
    local dict_help
    [[ -s "$fasta" ]] || { warn "missing FASTA: $fasta"; return 0; }
    dict_help=$(samtools dict 2>&1 || true)
    if [[ "$dict_help" != *"Create a sequence dictionary"* ]]; then
        warn "samtools dict is unavailable; skip $(basename "$dict_out")"
        return 0
    fi
    run_if_missing "sequence dictionary $(basename "$dict_out")" "$dict_out" samtools dict -o "$dict_out" "$fasta"
}

build_bowtie2() {
    local fasta="$1"
    [[ -s "$fasta" ]] || { warn "missing FASTA: $fasta"; return 0; }
    run_if_missing "bowtie2 index $(basename "$fasta")" \
        "${fasta}.1.bt2,${fasta}.2.bt2,${fasta}.3.bt2,${fasta}.4.bt2,${fasta}.rev.1.bt2,${fasta}.rev.2.bt2" \
        bowtie2-build "$fasta" "$fasta"
}

build_blastdb() {
    local fasta="$1"
    local out_prefix="$2"
    [[ -s "$fasta" ]] || { warn "missing FASTA: $fasta"; return 0; }
    if ! command -v makeblastdb >/dev/null 2>&1; then
        warn "makeblastdb is unavailable; skip BLAST DB for $(basename "$fasta")"
        return 0
    fi
    run_if_missing "blast db $(basename "$out_prefix")" \
        "${out_prefix}.nhr,${out_prefix}.nin,${out_prefix}.nsq" \
        makeblastdb -in "$fasta" -dbtype nucl -out "$out_prefix" -parse_seqids
}

need_tool samtools
need_tool bwa
need_tool bowtie2-build

log "resource root: $RESOURCE_ROOT"

build_faidx "$COMBINED_REF"
build_bwa "$COMBINED_REF"
build_dict "$COMBINED_REF"
build_bowtie2 "$IMGT_FASTA"

for gene in "${GENES[@]}"; do
    gene_fasta="${HLA_DIR}/${gene}/${gene}.fa"
    build_faidx "$gene_fasta"
    build_bwa "$gene_fasta"
    build_blastdb "$gene_fasta" "${HLA_DIR}/${gene}/${gene}"
done

if [[ -d "$EXON_DIR" ]]; then
    for gene in "${GENES[@]}"; do
        exon_fasta="${EXON_DIR}/${gene}.fasta"
        build_faidx "$exon_fasta"
        build_blastdb "$exon_fasta" "$exon_fasta"
    done
else
    warn "missing exon directory: $EXON_DIR"
fi

log "resource index build/check complete"
