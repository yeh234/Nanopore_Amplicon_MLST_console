#!/bin/bash
set -euo pipefail

sample_fastq_dir="output/01_demultiplex/merged_by_sample"
output_dir="output/02_assembly"
threads=8
genome_size="4k"
asm_coverage=40
min_overlap=1000
samples_csv=""
skip_polish=0
medaka_model=""
medaka_bacteria=0

usage() {
cat <<'EOF'
Usage:
  bash code/02_assembly_polish_pipeline.sh \
    [--sample-fastq-dir=output/01_demultiplex/merged_by_sample] \
    [--output-dir=output/02_assembly] \
    [--threads=8] \
    [--genome-size=4k] \
    [--asm-coverage=40] \
    [--min-overlap=1000] \
    [--samples=sample1,sample2] \
    [--skip-polish] \
    [--medaka-bacteria] \
    [--medaka-model=r1041_e82_400bps_sup_v5.0.0]

Inputs:
  - merged sample FASTQ files such as sample1.fastq.gz, sample2.fastq.gz

Outputs:
  - per-sample Flye assembly
  - read-to-assembly BAM
  - per-sample polished consensus FASTA
  - read-to-polished BAM
  - assembly_polish_summary.tsv
EOF
}

die() {
  echo "Error: $*" >&2
  exit 1
}

need_cmd() {
  command -v "$1" >/dev/null 2>&1 || die "Required command not found: $1"
}

count_fastq_reads() {
  local path="$1"
  if [[ "$path" == *.gz ]]; then
    gzip -dc "$path" | awk 'END{print NR/4}'
  else
    awk 'END{print NR/4}' "$path"
  fi
}

sample_selected() {
  local sample_name="$1"
  if [[ -z "$samples_csv" ]]; then
    return 0
  fi

  IFS=',' read -r -a wanted <<< "$samples_csv"
  local item
  for item in "${wanted[@]}"; do
    if [[ "$item" == "$sample_name" ]]; then
      return 0
    fi
  done
  return 1
}

discover_sample_fastqs() {
  local dir="$1"
  find "$dir" -maxdepth 1 -type f \( -name 'sample*.fastq.gz' -o -name 'sample*.fastq' \) | sort
}

run_flye() {
  local reads_path="$1"
  local flye_dir="$2"

  flye \
    --nano-raw "$reads_path" \
    --threads "$threads" \
    --min-overlap "$min_overlap" \
    --asm-coverage "$asm_coverage" \
    --genome-size "$genome_size" \
    --out-dir "$flye_dir"
}

run_bwa_alignment() {
  local reference_fasta="$1"
  local reads_path="$2"
  local output_bam="$3"

  bwa index "$reference_fasta"
  bwa mem "$reference_fasta" "$reads_path" | samtools sort -o "$output_bam"
  samtools index "$output_bam"
}

run_medaka_polish() {
  local reads_path="$1"
  local draft_fasta="$2"
  local draft_bam="$3"
  local sample_dir="$4"
  local sample_name="$5"

  local polished_fasta="$sample_dir/${sample_name}_assembly_polish_consensus.fasta"

  # Prefer the original repo's medaka consensus/stitch workflow when the installed
  # medaka version still supports those subcommands. On medaka 2.x this legacy CLI
  # is gone, so we fall back to medaka_consensus.
  if command -v medaka >/dev/null 2>&1 && medaka consensus -h >/dev/null 2>&1 && medaka stitch -h >/dev/null 2>&1; then
    local polish_hdf="$sample_dir/${sample_name}_assembly_polish.hdf"
    medaka consensus "$draft_bam" "$polish_hdf" >&2
    medaka stitch "$polish_hdf" "$draft_fasta" "$polished_fasta" >&2
  elif command -v medaka_consensus >/dev/null 2>&1; then
    local medaka_dir="$sample_dir/medaka_consensus"
    mkdir -p "$medaka_dir"

    local -a medaka_cmd=(medaka_consensus -i "$reads_path" -d "$draft_fasta" -o "$medaka_dir" -t "$threads")
    if [[ "$medaka_bacteria" -eq 1 ]]; then
      medaka_cmd+=(--bacteria)
    fi
    if [[ -n "$medaka_model" ]]; then
      medaka_cmd+=(-m "$medaka_model")
    fi

    "${medaka_cmd[@]}" >&2
    cp "$medaka_dir/consensus.fasta" "$polished_fasta"
  else
    die "Neither medaka_consensus nor medaka is available. Install one of them or use --skip-polish."
  fi

  printf '%s\n' "$polished_fasta"
}

main() {
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --sample-fastq-dir=*)
        sample_fastq_dir="${1#*=}"
        ;;
      --output-dir=*)
        output_dir="${1#*=}"
        ;;
      --threads=*)
        threads="${1#*=}"
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
      --medaka-bacteria)
        medaka_bacteria=1
        ;;
      --medaka-model=*)
        medaka_model="${1#*=}"
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

  [[ -d "$sample_fastq_dir" ]] || die "Sample FASTQ directory not found: $sample_fastq_dir"

  need_cmd flye
  need_cmd bwa
  need_cmd samtools
  if [[ "$skip_polish" -eq 0 ]]; then
    if ! command -v medaka_consensus >/dev/null 2>&1 && ! command -v medaka >/dev/null 2>&1; then
      die "Polish requested, but medaka_consensus and medaka are both unavailable."
    fi
  fi

  mkdir -p "$output_dir"
  local summary_path="$output_dir/assembly_polish_summary.tsv"
  printf 'sample\treads_fastq\tread_count\tdraft_fasta\tdraft_bam\tpolished_fasta\tpolished_bam\tstatus\n' > "$summary_path"

  local processed=0
  local reads_path sample_name read_count sample_dir flye_dir draft_fasta draft_bam polished_fasta polished_bam

  while IFS= read -r reads_path; do
    sample_name="$(basename "$reads_path")"
    sample_name="${sample_name%.gz}"
    sample_name="${sample_name%.fastq}"

    if ! sample_selected "$sample_name"; then
      continue
    fi

    read_count="$(count_fastq_reads "$reads_path")"
    if [[ "$read_count" == "0" ]]; then
      printf '%s\t%s\t%s\t\t\t\t\tskipped_empty\n' "$sample_name" "$reads_path" "$read_count" >> "$summary_path"
      continue
    fi

    processed=$((processed + 1))
    sample_dir="$output_dir/$sample_name"
    flye_dir="$sample_dir/flye"
    mkdir -p "$sample_dir"

    echo "Processing $sample_name"
    run_flye "$reads_path" "$flye_dir"

    draft_fasta="$sample_dir/${sample_name}_assembly.fasta"
    cp "$flye_dir/assembly.fasta" "$draft_fasta"

    draft_bam="$sample_dir/${sample_name}_assembly_output.bam"
    run_bwa_alignment "$draft_fasta" "$reads_path" "$draft_bam"

    polished_fasta=""
    polished_bam=""
    if [[ "$skip_polish" -eq 0 ]]; then
      polished_fasta="$(run_medaka_polish "$reads_path" "$draft_fasta" "$draft_bam" "$sample_dir" "$sample_name")"
      polished_bam="$sample_dir/${sample_name}_assembly_polish_output.bam"
      run_bwa_alignment "$polished_fasta" "$reads_path" "$polished_bam"
      printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\tassembled_and_polished\n' \
        "$sample_name" "$reads_path" "$read_count" "$draft_fasta" "$draft_bam" "$polished_fasta" "$polished_bam" >> "$summary_path"
    else
      printf '%s\t%s\t%s\t%s\t%s\t\t\tassembled_only\n' \
        "$sample_name" "$reads_path" "$read_count" "$draft_fasta" "$draft_bam" >> "$summary_path"
    fi
  done < <(discover_sample_fastqs "$sample_fastq_dir")

  [[ "$processed" -gt 0 ]] || die "No matching sample FASTQ files found in $sample_fastq_dir"
  echo "Wrote $summary_path"
}

main "$@"
