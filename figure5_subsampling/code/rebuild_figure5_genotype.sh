#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

# =========================================================
# Rebuild Figure 5 genotype-result input from one sample
# =========================================================

# -----------------------------
# repository-local defaults
# -----------------------------
READS_FASTQ="${READS_FASTQ:-${REPO_DIR}/../output/01_demultiplex/merged_by_sample/sample4.fastq.gz}"
ASSEMBLY_FASTA="${ASSEMBLY_FASTA:-${REPO_DIR}/../output/02_assembly/sample4/sample4_assembly_polish_consensus.fasta}"
SAMPLE_NAME="${SAMPLE_NAME:-}"
READ_COUNT="${READ_COUNT:-}"
REPLICATE_ID="${REPLICATE_ID:-1}"
SAMPLE_SOURCE="${SAMPLE_SOURCE:-}"
RESULT_DIR="${RESULT_DIR:-.}"
OUTPUT_NAME="${OUTPUT_NAME:-sample4_patch_genotype_result.txt}"
BAM_DIR="${BAM_DIR:-${RESULT_DIR}/bam_for_igv}"

# -----------------------------
# MLST reference
# -----------------------------
MLST_FASTA="${MLST_FASTA:-}"
MLST_DB_DIR="${MLST_DB_DIR:-${REPO_DIR}/../db}"
COMBINED_MLST_FASTA_NAME="${COMBINED_MLST_FASTA_NAME:-all_MLST_local.fas}"

# -----------------------------
# tools
# -----------------------------
MINIMAP2="${MINIMAP2:-}"
SAMTOOLS="${SAMTOOLS:-}"
VERBOSE="${VERBOSE:-1}"

# 1 = also write sorted BAM/BAI file for IGV inspection
WRITE_BAM="${WRITE_BAM:-0}"

GENO_OUT="${RESULT_DIR}/${OUTPUT_NAME}"

log() {
  if [[ "${VERBOSE}" -eq 1 ]]; then
    printf "%s\n" "$1"
  fi
}

uppercase() {
  printf "%s" "$1" | tr '[:lower:]' '[:upper:]'
}

require_file() {
  local path="$1"
  local label="$2"

  [[ -f "${path}" ]] || {
    echo "ERROR: ${label} not found: ${path}" >&2
    exit 1
  }
}

resolve_tool() {
  local requested="$1"
  local tool_name="$2"

  if [[ -n "${requested}" ]]; then
    if [[ "${requested}" == */* ]]; then
      [[ -x "${requested}" ]] || {
        echo "ERROR: ${tool_name} is not executable: ${requested}" >&2
        exit 1
      }
      printf "%s\n" "${requested}"
      return
    fi

    command -v "${requested}" >/dev/null 2>&1 || {
      echo "ERROR: ${tool_name} not found in PATH: ${requested}" >&2
      exit 1
    }
    command -v "${requested}"
    return
  fi

  if command -v "${tool_name}" >/dev/null 2>&1; then
    command -v "${tool_name}"
    return
  fi

  echo "ERROR: ${tool_name} not found. Set $(uppercase "${tool_name}") or add it to PATH." >&2
  exit 1
}

require_dir() {
  local path="$1"
  local label="$2"

  [[ -d "${path}" ]] || {
    echo "ERROR: ${label} not found: ${path}" >&2
    exit 1
  }
}

infer_sample_name() {
  local path="$1"
  local base

  base="$(basename "${path}")"
  base="${base%.fastq.gz}"
  base="${base%.fq.gz}"
  base="${base%.fastq}"
  base="${base%.fq}"
  base="${base%_assembly_polish_consensus.fasta}"
  base="${base%_assembly_consensus.fasta}"
  base="${base%_assembly.fasta}"
  base="${base%.fasta}"
  base="${base%.fa}"

  printf "%s\n" "${base}"
}

count_fastq_reads() {
  local fastq="$1"

  if [[ "${fastq}" == *.gz ]]; then
    gzip -cd -- "${fastq}" | awk 'END {print NR/4}'
  else
    awk 'END {print NR/4}' "${fastq}"
  fi
}

build_mlst_fasta() {
  local output_fasta="$1"
  local gene
  local -a genes=(arcC aroE glpF gmk pta tpi yqiL)

  : > "${output_fasta}"
  for gene in "${genes[@]}"; do
    require_file "${MLST_DB_DIR}/${gene}.fas" "MLST allele fasta"
    cat "${MLST_DB_DIR}/${gene}.fas" >> "${output_fasta}"
  done
}

mkdir -p "${RESULT_DIR}"

require_file "${ASSEMBLY_FASTA}" "assembly fasta"

if [[ -n "${READS_FASTQ}" ]]; then
  require_file "${READS_FASTQ}" "reads fastq"
fi

require_dir "${MLST_DB_DIR}" "MLST database directory"
MINIMAP2="$(resolve_tool "${MINIMAP2}" "minimap2")"
SAMTOOLS="$(resolve_tool "${SAMTOOLS}" "samtools")"

if [[ -z "${MLST_FASTA}" ]]; then
  MLST_FASTA="${RESULT_DIR}/${COMBINED_MLST_FASTA_NAME}"
  build_mlst_fasta "${MLST_FASTA}"
else
  require_file "${MLST_FASTA}" "MLST reference fasta"
fi

if [[ -z "${SAMPLE_NAME}" ]]; then
  if [[ -n "${READS_FASTQ}" ]]; then
    SAMPLE_NAME="$(infer_sample_name "${READS_FASTQ}")"
  else
    SAMPLE_NAME="$(infer_sample_name "${ASSEMBLY_FASTA}")"
  fi
fi

if [[ -z "${READ_COUNT}" ]]; then
  if [[ -n "${READS_FASTQ}" ]]; then
    READ_COUNT="$(count_fastq_reads "${READS_FASTQ}")"
  else
    echo "ERROR: READ_COUNT is required when READS_FASTQ is not provided." >&2
    exit 1
  fi
fi

if [[ -z "${SAMPLE_SOURCE}" ]]; then
  SAMPLE_SOURCE="${SAMPLE_NAME}_${READ_COUNT}_${REPLICATE_ID}"
fi

if [[ "${WRITE_BAM}" -eq 1 ]]; then
  mkdir -p "${BAM_DIR}"
fi

log "Assembly fasta : ${ASSEMBLY_FASTA}"
log "Reads fastq    : ${READS_FASTQ:-<not provided>}"
log "MLST fasta     : ${MLST_FASTA}"
log "Sample source  : ${SAMPLE_SOURCE}"
log "Writing output : ${GENO_OUT}"

"${MINIMAP2}" -a -x sr "${ASSEMBLY_FASTA}" "${MLST_FASTA}" \
  | "${SAMTOOLS}" view -d NM:0 \
  | awk -v sample="${SAMPLE_SOURCE}" '
      BEGIN { OFS = "\t" }
      {
        nm = ms = as = nn = tp = ""
        for (i = 12; i <= NF; i++) {
          if ($i ~ /^NM:/) nm = $i
          else if ($i ~ /^ms:/) ms = $i
          else if ($i ~ /^AS:/) as = $i
          else if ($i ~ /^nn:/) nn = $i
          else if ($i ~ /^tp:/) tp = $i
        }
        print $1, sample, $6, nm, ms, as, nn, tp
      }
    ' > "${GENO_OUT}"

if [[ "${WRITE_BAM}" -eq 1 ]]; then
  bam="${BAM_DIR}/${SAMPLE_SOURCE}_MLST_type.bam"

  "${MINIMAP2}" -a -x sr "${ASSEMBLY_FASTA}" "${MLST_FASTA}" \
    | "${SAMTOOLS}" view -S -b \
    | "${SAMTOOLS}" sort -o "${bam}"

  "${SAMTOOLS}" index "${bam}"
  log "BAM written    : ${bam}"
fi

log "Done."
