#!/usr/bin/env bash
set -euo pipefail

# =========================================================
# Rebuild full Figure 5 batch workflow for sample4
# =========================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

# -----------------------------
# repository-local defaults
# -----------------------------
INPUT_FASTQ="${INPUT_FASTQ:-${REPO_DIR}/../output/01_demultiplex/merged_by_sample/sample4.fastq.gz}"
SAMPLE_NAME="${SAMPLE_NAME:-sample4}"
RESULT_DIR="${RESULT_DIR:-${REPO_DIR}/output/figure5_rebuild}"
WORK_DIR="${WORK_DIR:-${RESULT_DIR}/work}"
MASTER_OUTPUT_NAME="${MASTER_OUTPUT_NAME:-sample4_patch_genotype_result.txt}"
MANIFEST_NAME="${MANIFEST_NAME:-sample4_patch_manifest.tsv}"
GENOTYPE_SCRIPT="${GENOTYPE_SCRIPT:-${SCRIPT_DIR}/rebuild_figure5_genotype.sh}"
BAM_DIR="${BAM_DIR:-${RESULT_DIR}/bam_for_igv}"

# -----------------------------
# MLST reference
# -----------------------------
MLST_FASTA="${MLST_FASTA:-}"
MLST_DB_DIR="${MLST_DB_DIR:-${REPO_DIR}/../db}"
COMBINED_MLST_FASTA_NAME="${COMBINED_MLST_FASTA_NAME:-all_MLST_local.fas}"

# -----------------------------
# batch parameters
# -----------------------------
READ_VALUES="${READ_VALUES:-}"
READ_START="${READ_START:-10}"
READ_END="${READ_END:-60}"
READ_STEP="${READ_STEP:-1}"
ITER_START="${ITER_START:-1}"
ITER_END="${ITER_END:-100}"

# -----------------------------
# assembly / polishing parameters
# -----------------------------
FLYE_THREADS="${FLYE_THREADS:-4}"
MEDAKA_THREADS="${MEDAKA_THREADS:-1}"
GENOME_SIZE="${GENOME_SIZE:-4k}"
ASM_COVERAGE="${ASM_COVERAGE:-40}"
MIN_OVERLAP="${MIN_OVERLAP:-1000}"
MEDAKA_MODEL="${MEDAKA_MODEL:-r941_min_sup_g507}"
WRITE_BAM="${WRITE_BAM:-0}"
FAIL_FAST="${FAIL_FAST:-0}"

# -----------------------------
# tools
# -----------------------------
SEQKIT="${SEQKIT:-}"
FLYE="${FLYE:-}"
MEDAKA_CONSENSUS="${MEDAKA_CONSENSUS:-}"
MINIMAP2="${MINIMAP2:-}"
SAMTOOLS="${SAMTOOLS:-}"

MASTER_OUTPUT="${RESULT_DIR}/${MASTER_OUTPUT_NAME}"
MANIFEST_PATH="${RESULT_DIR}/${MANIFEST_NAME}"

log() {
  printf "%s\n" "$1"
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

run_with_logs() {
  local stdout_log="$1"
  local stderr_log="$2"
  shift 2

  "$@" >"${stdout_log}" 2>"${stderr_log}"
}

sample_reads() {
  local reads="$1"
  local replicate="$2"
  local output_fastq="$3"

  "${SEQKIT}" sample2 \
    -2 \
    -n "${reads}" \
    -s "${replicate}" \
    -o "${output_fastq}" \
    "${INPUT_FASTQ}"
}

handle_failure() {
  local sample_tag="$1"
  local status="$2"
  local genotype_tsv="$3"

  printf "%s\n" "${status}" > "${genotype_tsv}.status"
  : > "${genotype_tsv}"
  log "[${sample_tag}] ${status}"

  if [[ "${FAIL_FAST}" -eq 1 ]]; then
    echo "ERROR: ${sample_tag} failed with status ${status}." >&2
    exit 1
  fi
}

declare -a READ_COUNTS=()

if [[ -n "${READ_VALUES}" ]]; then
  read -r -a READ_COUNTS <<< "${READ_VALUES}"
else
  for reads in $(seq "${READ_START}" "${READ_STEP}" "${READ_END}"); do
    READ_COUNTS+=("${reads}")
  done
fi

[[ "${#READ_COUNTS[@]}" -gt 0 ]] || {
  echo "ERROR: no read counts selected." >&2
  exit 1
}

require_file "${INPUT_FASTQ}" "input fastq"
require_file "${GENOTYPE_SCRIPT}" "genotype script"

require_dir "${MLST_DB_DIR}" "MLST database directory"
SEQKIT="$(resolve_tool "${SEQKIT}" "seqkit")"
FLYE="$(resolve_tool "${FLYE}" "flye")"
MEDAKA_CONSENSUS="$(resolve_tool "${MEDAKA_CONSENSUS}" "medaka_consensus")"
MINIMAP2="$(resolve_tool "${MINIMAP2}" "minimap2")"
SAMTOOLS="$(resolve_tool "${SAMTOOLS}" "samtools")"

export PATH="$(dirname "${SEQKIT}")":"$(dirname "${FLYE}")":"$(dirname "${MEDAKA_CONSENSUS}")":"$(dirname "${MINIMAP2}")":"$(dirname "${SAMTOOLS}")":"${PATH}"

mkdir -p "${RESULT_DIR}"
mkdir -p "${WORK_DIR}"

if [[ -z "${MLST_FASTA}" ]]; then
  MLST_FASTA="${RESULT_DIR}/${COMBINED_MLST_FASTA_NAME}"
  if [[ ! -f "${MLST_FASTA}" ]]; then
    build_mlst_fasta "${MLST_FASTA}"
  fi
else
  require_file "${MLST_FASTA}" "MLST reference fasta"
fi

total_reads="$(count_fastq_reads "${INPUT_FASTQ}")"

for reads in "${READ_COUNTS[@]}"; do
  if (( reads > total_reads )); then
    echo "ERROR: requested reads (${reads}) exceed total reads in ${INPUT_FASTQ} (${total_reads})." >&2
    exit 1
  fi
done

log "Input fastq     : ${INPUT_FASTQ}"
log "Total reads     : ${total_reads}"
log "Read counts     : ${READ_COUNTS[*]}"
log "Replicate range : ${ITER_START}-${ITER_END}"
log "Result dir      : ${RESULT_DIR}"
log "Work dir        : ${WORK_DIR}"
log "MLST fasta      : ${MLST_FASTA}"
log "Fail fast       : ${FAIL_FAST}"

for reads in "${READ_COUNTS[@]}"; do
  fastq_dir="${WORK_DIR}/${SAMPLE_NAME}_fastq/${SAMPLE_NAME}_${reads}"
  assembly_root="${WORK_DIR}/${SAMPLE_NAME}_assembly/${SAMPLE_NAME}_${reads}_assembly"
  mkdir -p "${fastq_dir}"
  mkdir -p "${assembly_root}"

  for replicate in $(seq "${ITER_START}" "${ITER_END}"); do
    sample_tag="${SAMPLE_NAME}_${reads}_${replicate}"
    run_dir="${assembly_root}/${SAMPLE_NAME}_${reads}_assembly_${replicate}"
    flye_dir="${run_dir}/flye"
    medaka_dir="${run_dir}/medaka"
    log_dir="${run_dir}/logs"
    sample_fastq="${fastq_dir}/${sample_tag}.fastq.gz"
    assembly_fasta="${run_dir}/${sample_tag}.assembly.fasta"
    consensus_fasta="${run_dir}/${sample_tag}_polish_consensus.fasta"
    genotype_tsv="${run_dir}/${sample_tag}.genotype.tsv"
    status_file="${genotype_tsv}.status"
    current_status=""

    mkdir -p "${run_dir}"
    mkdir -p "${log_dir}"

    log "[${sample_tag}] preparing subsample"
    sampled_n=""
    if [[ -f "${sample_fastq}" ]]; then
      sampled_n="$(count_fastq_reads "${sample_fastq}")"
    fi

    if [[ ! -f "${sample_fastq}" || "${sampled_n}" != "${reads}" ]]; then
      if [[ -f "${sample_fastq}" ]]; then
        log "[${sample_tag}] existing subsample has ${sampled_n} reads; rebuilding"
        rm -f "${sample_fastq}"
      fi
      sample_reads "${reads}" "${replicate}" "${sample_fastq}"
      sampled_n="$(count_fastq_reads "${sample_fastq}")"
      [[ "${sampled_n}" == "${reads}" ]] || {
        echo "ERROR: subsampling for ${sample_tag} produced ${sampled_n} reads instead of ${reads}." >&2
        exit 1
      }
    fi

    log "[${sample_tag}] assembling with flye"
    if [[ ! -f "${assembly_fasta}" ]]; then
      mkdir -p "${flye_dir}"
      if run_with_logs \
        "${log_dir}/flye.stdout.log" \
        "${log_dir}/flye.stderr.log" \
        "${FLYE}" \
        --nano-raw "${sample_fastq}" \
        --threads "${FLYE_THREADS}" \
        --min-overlap "${MIN_OVERLAP}" \
        --asm-coverage "${ASM_COVERAGE}" \
        --genome-size "${GENOME_SIZE}" \
        --out-dir "${flye_dir}" \
        && [[ -f "${flye_dir}/assembly.fasta" ]]; then
        cp "${flye_dir}/assembly.fasta" "${assembly_fasta}"
      else
        handle_failure "${sample_tag}" "assembly_failed" "${genotype_tsv}"
        continue
      fi
    fi

    log "[${sample_tag}] polishing with medaka"
    if [[ ! -f "${consensus_fasta}" ]]; then
      mkdir -p "${medaka_dir}"
      if run_with_logs \
        "${log_dir}/medaka.stdout.log" \
        "${log_dir}/medaka.stderr.log" \
        "${MEDAKA_CONSENSUS}" \
        -i "${sample_fastq}" \
        -d "${assembly_fasta}" \
        -o "${medaka_dir}" \
        -m "${MEDAKA_MODEL}" \
        -t "${MEDAKA_THREADS}" \
        && [[ -f "${medaka_dir}/consensus.fasta" ]]; then
        cp "${medaka_dir}/consensus.fasta" "${consensus_fasta}"
      else
        handle_failure "${sample_tag}" "medaka_failed" "${genotype_tsv}"
        continue
      fi
    fi

    log "[${sample_tag}] genotyping"
    if [[ -f "${status_file}" ]]; then
      current_status="$(cat "${status_file}")"
    fi

    if [[ ! -f "${genotype_tsv}" || ! -f "${status_file}" || "${current_status}" != "ok" ]]; then
      if run_with_logs \
        "${log_dir}/genotype.stdout.log" \
        "${log_dir}/genotype.stderr.log" \
        env \
        VERBOSE=0 \
        READS_FASTQ="${sample_fastq}" \
        ASSEMBLY_FASTA="${consensus_fasta}" \
        SAMPLE_NAME="${SAMPLE_NAME}" \
        READ_COUNT="${reads}" \
        REPLICATE_ID="${replicate}" \
        RESULT_DIR="${run_dir}" \
        OUTPUT_NAME="$(basename "${genotype_tsv}")" \
        MLST_FASTA="${MLST_FASTA}" \
        WRITE_BAM="${WRITE_BAM}" \
        BAM_DIR="${BAM_DIR}" \
        MINIMAP2="${MINIMAP2}" \
        SAMTOOLS="${SAMTOOLS}" \
        bash "${GENOTYPE_SCRIPT}"; then
        printf "%s\n" "ok" > "${status_file}"
      else
        handle_failure "${sample_tag}" "genotype_failed" "${genotype_tsv}"
        continue
      fi
    fi
  done
done

: > "${MASTER_OUTPUT}"
printf "Sample\tReads\tReplicate\tStatus\tFastq\tAssembly\tConsensus\tGenotype\n" > "${MANIFEST_PATH}"

for reads in "${READ_COUNTS[@]}"; do
  assembly_root="${WORK_DIR}/${SAMPLE_NAME}_assembly/${SAMPLE_NAME}_${reads}_assembly"

  for replicate in $(seq "${ITER_START}" "${ITER_END}"); do
    sample_tag="${SAMPLE_NAME}_${reads}_${replicate}"
    run_dir="${assembly_root}/${SAMPLE_NAME}_${reads}_assembly_${replicate}"
    sample_fastq="${WORK_DIR}/${SAMPLE_NAME}_fastq/${SAMPLE_NAME}_${reads}/${sample_tag}.fastq.gz"
    assembly_fasta="${run_dir}/${sample_tag}.assembly.fasta"
    consensus_fasta="${run_dir}/${sample_tag}_polish_consensus.fasta"
    genotype_tsv="${run_dir}/${sample_tag}.genotype.tsv"
    status_file="${genotype_tsv}.status"
    status="missing"

    require_file "${genotype_tsv}" "per-replicate genotype table"
    if [[ -f "${status_file}" ]]; then
      status="$(cat "${status_file}")"
    fi
    cat "${genotype_tsv}" >> "${MASTER_OUTPUT}"
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
      "${SAMPLE_NAME}" \
      "${reads}" \
      "${replicate}" \
      "${status}" \
      "${sample_fastq}" \
      "${assembly_fasta}" \
      "${consensus_fasta}" \
      "${genotype_tsv}" \
      >> "${MANIFEST_PATH}"
  done
done

log "Master genotype table : ${MASTER_OUTPUT}"
log "Manifest              : ${MANIFEST_PATH}"
log "Done."
