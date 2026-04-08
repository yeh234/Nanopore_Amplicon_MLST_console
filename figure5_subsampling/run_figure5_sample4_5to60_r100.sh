#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="${SCRIPT_DIR}"
cd "${REPO_DIR}"

INPUT_FASTQ="${INPUT_FASTQ:-${REPO_DIR}/../output/01_demultiplex/merged_by_sample/sample4.fastq.gz}"
SAMPLE_NAME="${SAMPLE_NAME:-sample4}"
OUT_DIR="${OUT_DIR:-output}"
WORK_DIR="${WORK_DIR:-${OUT_DIR}/work}"

FLYE_THREADS="${FLYE_THREADS:-1}"
MEDAKA_THREADS="${MEDAKA_THREADS:-1}"
WRITE_BAM="${WRITE_BAM:-0}"
FAIL_FAST="${FAIL_FAST:-0}"

GUIDE_READ_LABEL="${GUIDE_READ_LABEL:-Approximate transition into high-correctness range}"
GUIDE_READ_SECONDARY_LABEL="${GUIDE_READ_SECONDARY_LABEL:-All loci >=99% correctness}"
BATCH_SCRIPT="${SCRIPT_DIR}/code/rebuild_figure5_batch.sh"
PLOT_SCRIPT="${SCRIPT_DIR}/code/MAF_revised.R"

[[ -f "${INPUT_FASTQ}" ]] || {
  echo "ERROR: input FASTQ not found: ${INPUT_FASTQ}" >&2
  exit 1
}

mkdir -p "${OUT_DIR}"

RESULT_DIR="${OUT_DIR}" \
WORK_DIR="${WORK_DIR}" \
INPUT_FASTQ="${INPUT_FASTQ}" \
SAMPLE_NAME="${SAMPLE_NAME}" \
MASTER_OUTPUT_NAME="sample4_patch_genotype_result_full_r100.txt" \
MANIFEST_NAME="sample4_patch_manifest_full_r100.tsv" \
READ_START=5 \
READ_END=60 \
READ_STEP=1 \
ITER_START=1 \
ITER_END=100 \
FLYE_THREADS="${FLYE_THREADS}" \
MEDAKA_THREADS="${MEDAKA_THREADS}" \
WRITE_BAM="${WRITE_BAM}" \
FAIL_FAST="${FAIL_FAST}" \
bash "${BATCH_SCRIPT}"

INPUT_FILE="${OUT_DIR}/sample4_patch_genotype_result_full_r100.txt" \
OUTPUT_PREFIX="${OUT_DIR}/Figure5_sup_5to60_full_r100" \
SUMMARY_FILE="${OUT_DIR}/Figure5_sup_5to60_full_r100_summary.tsv" \
LOW_READ_ZOOM_MIN=5 \
LOW_READ_ZOOM_MAX=60 \
Y_AXIS_MIN=0 \
Y_AXIS_MAX=1.02 \
LABEL_READ_MIN=5 \
LABEL_READ_MAX=15 \
FORCE_LABEL_READS=35,45,60 \
GUIDE_READ=35 \
GUIDE_READ_LABEL="${GUIDE_READ_LABEL}" \
GUIDE_READ_SECONDARY=45 \
GUIDE_READ_SECONDARY_LABEL="${GUIDE_READ_SECONDARY_LABEL}" \
Rscript "${PLOT_SCRIPT}"

printf "Output dir  : %s\n" "${OUT_DIR}"
printf "Figure PDF  : %s\n" "${OUT_DIR}/Figure5_sup_5to60_full_r100.pdf"
printf "Figure PNG  : %s\n" "${OUT_DIR}/Figure5_sup_5to60_full_r100.png"
printf "Summary TSV : %s\n" "${OUT_DIR}/Figure5_sup_5to60_full_r100_summary.tsv"
