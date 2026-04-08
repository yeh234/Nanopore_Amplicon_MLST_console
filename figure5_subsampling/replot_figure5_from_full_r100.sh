#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="${SCRIPT_DIR}"
cd "${REPO_DIR}"

INPUT_FILE="${INPUT_FILE:-output/sample4_patch_genotype_result_full_r100.txt}"
OUTPUT_PREFIX="${OUTPUT_PREFIX:-output/Figure5_sup_5to60_full_r100}"
SUMMARY_FILE="${SUMMARY_FILE:-output/Figure5_sup_5to60_full_r100_summary.tsv}"
PLOT_SCRIPT="${SCRIPT_DIR}/code/MAF_revised.R"

LOW_READ_ZOOM_MIN="${LOW_READ_ZOOM_MIN:-5}"
LOW_READ_ZOOM_MAX="${LOW_READ_ZOOM_MAX:-60}"
Y_AXIS_MIN="${Y_AXIS_MIN:-0}"
Y_AXIS_MAX="${Y_AXIS_MAX:-1.02}"
LABEL_READ_MIN="${LABEL_READ_MIN:-5}"
LABEL_READ_MAX="${LABEL_READ_MAX:-15}"
FORCE_LABEL_READS="${FORCE_LABEL_READS:-35,45,60}"
GUIDE_READ="${GUIDE_READ:-35}"
GUIDE_READ_LABEL="${GUIDE_READ_LABEL:-Approximate transition into high-correctness range}"
GUIDE_READ_SECONDARY="${GUIDE_READ_SECONDARY:-45}"
GUIDE_READ_SECONDARY_LABEL="${GUIDE_READ_SECONDARY_LABEL:-All loci >=99% correctness}"

[[ -f "${INPUT_FILE}" ]] || {
  echo "ERROR: genotype table not found: ${INPUT_FILE}" >&2
  exit 1
}

[[ -f "${PLOT_SCRIPT}" ]] || {
  echo "ERROR: plotting script not found: ${PLOT_SCRIPT}" >&2
  exit 1
}

mkdir -p "$(dirname "${OUTPUT_PREFIX}")" "$(dirname "${SUMMARY_FILE}")"

INPUT_FILE="${INPUT_FILE}" \
OUTPUT_PREFIX="${OUTPUT_PREFIX}" \
SUMMARY_FILE="${SUMMARY_FILE}" \
LOW_READ_ZOOM_MIN="${LOW_READ_ZOOM_MIN}" \
LOW_READ_ZOOM_MAX="${LOW_READ_ZOOM_MAX}" \
Y_AXIS_MIN="${Y_AXIS_MIN}" \
Y_AXIS_MAX="${Y_AXIS_MAX}" \
LABEL_READ_MIN="${LABEL_READ_MIN}" \
LABEL_READ_MAX="${LABEL_READ_MAX}" \
FORCE_LABEL_READS="${FORCE_LABEL_READS}" \
GUIDE_READ="${GUIDE_READ}" \
GUIDE_READ_LABEL="${GUIDE_READ_LABEL}" \
GUIDE_READ_SECONDARY="${GUIDE_READ_SECONDARY}" \
GUIDE_READ_SECONDARY_LABEL="${GUIDE_READ_SECONDARY_LABEL}" \
Rscript "${PLOT_SCRIPT}"

printf "Input file  : %s\n" "${INPUT_FILE}"
printf "Figure PDF  : %s\n" "${OUTPUT_PREFIX}.pdf"
printf "Figure PNG  : %s\n" "${OUTPUT_PREFIX}.png"
printf "Summary TSV : %s\n" "${SUMMARY_FILE}"
