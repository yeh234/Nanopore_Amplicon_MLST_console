#!/bin/bash
set -euo pipefail

sample_fastq_dir="output/01_demultiplex/merged_by_sample"
assembly_output_dir="output/02_assembly"
mlst_output_dir="output/03_mlst"
assembly_env=""
mlst_env=""
threads=4
mlst_threads=1
genome_size="4k"
asm_coverage=40
min_overlap=1000
samples_csv=""
skip_polish=0
skip_assembly=0
skip_mlst=0
medaka_model="r941_min_sup_g507"
medaka_bacteria=0
scheme="saureus"

usage() {
cat <<'EOF'
Usage:
  bash code/04_assembly_polish_pipeline.sh \
    [--sample-fastq-dir=output/01_demultiplex/merged_by_sample] \
    [--assembly-output-dir=output/02_assembly] \
    [--mlst-output-dir=output/03_mlst] \
    [--assembly-env=nanopore-assembly-mlst-macos-arm64] \
    [--mlst-env=mlst-only] \
    [--threads=4] \
    [--mlst-threads=1] \
    [--genome-size=4k] \
    [--asm-coverage=40] \
    [--min-overlap=1000] \
    [--samples=sample1,sample2] \
    [--skip-polish] \
    [--skip-assembly] \
    [--skip-mlst] \
    [--medaka-bacteria] \
    [--medaka-model=r941_min_sup_g507] \
    [--scheme=saureus]

This wrapper runs:
  Step 4-5: assembly + polish
  Step 6: standard MLST calling with mlst
EOF
}

die() {
  echo "Error: $*" >&2
  exit 1
}

run_with_optional_env() {
  local env_name="$1"
  shift
  if [[ -n "$env_name" ]]; then
    conda run -n "$env_name" "$@"
  else
    "$@"
  fi
}

main() {
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --sample-fastq-dir=*)
        sample_fastq_dir="${1#*=}"
        ;;
      --assembly-output-dir=*)
        assembly_output_dir="${1#*=}"
        ;;
      --mlst-output-dir=*)
        mlst_output_dir="${1#*=}"
        ;;
      --assembly-env=*)
        assembly_env="${1#*=}"
        ;;
      --mlst-env=*)
        mlst_env="${1#*=}"
        ;;
      --threads=*)
        threads="${1#*=}"
        ;;
      --mlst-threads=*)
        mlst_threads="${1#*=}"
        ;;
      --genome-size=*)
        genome_size="${1#*=}"
        ;;
      --asm-coverage=*)
        asm_coverage="${1#*=}"
        ;;
      --min-overlap=*)
        min_overlap="${1#*=}"
        ;;
      --samples=*)
        samples_csv="${1#*=}"
        ;;
      --skip-polish)
        skip_polish=1
        ;;
      --skip-assembly)
        skip_assembly=1
        ;;
      --skip-mlst)
        skip_mlst=1
        ;;
      --medaka-bacteria)
        medaka_bacteria=1
        ;;
      --medaka-model=*)
        medaka_model="${1#*=}"
        ;;
      --scheme=*)
        scheme="${1#*=}"
        ;;
      --help|-h)
        usage
        exit 0
        ;;
      *)
        usage >&2
        die "Unknown argument: $1"
        ;;
    esac
    shift
  done

  if [[ "$skip_assembly" -eq 0 ]]; then
    [[ -d "$sample_fastq_dir" ]] || die "Sample FASTQ directory not found: $sample_fastq_dir"
    assembly_cmd=(
      bash code/02_assembly_polish_pipeline.sh
      "--sample-fastq-dir=$sample_fastq_dir"
      "--output-dir=$assembly_output_dir"
      "--threads=$threads"
      "--genome-size=$genome_size"
      "--asm-coverage=$asm_coverage"
      "--min-overlap=$min_overlap"
    )
    if [[ -n "$samples_csv" ]]; then
      assembly_cmd+=("--samples=$samples_csv")
    fi
    if [[ "$skip_polish" -eq 1 ]]; then
      assembly_cmd+=("--skip-polish")
    fi
    if [[ "$medaka_bacteria" -eq 1 ]]; then
      assembly_cmd+=("--medaka-bacteria")
    fi
    if [[ -n "$medaka_model" ]]; then
      assembly_cmd+=("--medaka-model=$medaka_model")
    fi

    echo "[progress] step2 stage assembly"
    echo "Running Step 4-5 assembly + polish"
    run_with_optional_env "$assembly_env" "${assembly_cmd[@]}"
  fi

  if [[ "$skip_mlst" -eq 0 ]]; then
    [[ -d "$assembly_output_dir" ]] || die "Assembly output directory not found: $assembly_output_dir"
    mlst_cmd=(
      python3 code/03_mlst_standard_caller.py
      "--assembly-dir=$assembly_output_dir"
      "--output-dir=$mlst_output_dir"
      "--scheme=$scheme"
      "--threads=$mlst_threads"
      "--prefer-polished"
    )
    if [[ -n "$samples_csv" ]]; then
      mlst_cmd+=("--samples=$samples_csv")
    fi

    echo "[progress] step2 stage mlst"
    echo "Running Step 6 standard mlst caller"
    run_with_optional_env "$mlst_env" "${mlst_cmd[@]}"
  fi

  echo "[progress] step2 stage complete"
  echo "Final SUP Step 4-6 wrapper finished"
  echo "Assembly output: $assembly_output_dir"
  echo "MLST output: $mlst_output_dir"
}

main "$@"
