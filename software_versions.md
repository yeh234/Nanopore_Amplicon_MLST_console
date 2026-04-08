# Software Versions

This file records the software versions associated with the curated final workflow under `code/`.

It separates:
- the experiment-specific historical basecalling chain used to generate the FASTQ input
- the currently validated local workflow environments used for demultiplexing, assembly, polishing, and MLST calling in this repository

## 1. Experiment-Specific Sequencing and Basecalling

These settings describe the nanopore run and the historical basecalling chain used for the study dataset.

| Component | Version / model | Notes |
| --- | --- | --- |
| Nanopore platform | ONT MK1C | Sequencing instrument |
| Flow cell | FLO-MIN106 (R9.4.1) | Historical study chemistry |
| Raw signal format | FAST5 | Raw signal exported from the run |
| FAST5 to POD5 conversion | `pod5` `0.3.23` | Historical conversion step before basecalling |
| Basecaller | Dorado `v0.9.6+0949eb8d` | Run on `NVIDIA GeForce RTX 2080 Ti` |
| Dorado model | `dna_r9.4.1_e8_sup@v3.6` | Generated SUP FASTQ output |

## 2. Demultiplex and Read Processing

These tools are used by `code/01_sup_demultiplex.R`.

| Tool | Version | Role |
| --- | --- | --- |
| R | `4.4.3` | Demultiplex workflow runtime |
| Biostrings | `2.74.1` | Read and sequence processing in R |
| seqkit | `v2.13.0` | `seqkit amplicon` barcode-region extraction |
| bioawk | optional; version not recorded | Preferred FASTA export helper when available |

## 3. Assembly and Polish Environment

Environment file:
- `envs/macos-arm64-nanopore-assembly-mlst.yml`

Conda environment name:
- `nanopore-assembly-mlst-macos-arm64`

These tools are used by `code/02_assembly_polish_pipeline.sh` and `code/04_assembly_polish_pipeline.sh`.

| Tool | Version |
| --- | --- |
| Python | `3.12.13` |
| Flye | `2.9.6-b1802` |
| bwa | `0.7.19-r1273` |
| samtools | `1.23+htslib-1.23` |
| minimap2 | `2.30-r1287` |
| Medaka | `2.2.0` |

## 4. Standard MLST Environment

Environment file:
- `envs/macos-arm64-mlst-only.yml`

Conda environment name:
- `mlst-only`

These tools are used by `code/03_mlst_standard_caller.py` and the MLST stage in `code/04_assembly_polish_pipeline.sh`.

| Tool | Version |
| --- | --- |
| Python | `3.12.13` |
| mlst | `2.32.2` |

## 5. Figure 3 Barcode Screening

These tools are used by `figure3_barcode_screening/run_all.R` and the scripts
under `figure3_barcode_screening/scripts/`.

| Tool | Version | Role |
| --- | --- | --- |
| R | `4.4.3` | Figure 3 workflow runtime |
| dplyr | `1.1.4` | Table filtering, joins, and summaries |
| ggplot2 | `3.5.2` | Heatmap plotting and export |
| readr | `2.1.5` | CSV and TSV input/output |
| readxl | `1.4.5` | Default `barcode108_panel.xlsx` import |
| tidyr | `1.3.1` | Long/wide table reshaping |

## 6. Figure 5 Subsampling and Replot

These tools are used by `figure5_subsampling/run_figure5_sample4_5to60_r100.sh`,
`figure5_subsampling/replot_figure5_from_full_r100.sh`,
`figure5_subsampling/code/rebuild_figure5_batch.sh`,
`figure5_subsampling/code/rebuild_figure5_genotype.sh`, and
`figure5_subsampling/code/MAF_revised.R`.

| Tool | Version | Role |
| --- | --- | --- |
| R | `4.4.3` | Figure 5 plotting and replot runtime |
| seqkit | `v2.13.0` | Read-depth subsampling |
| Flye | `2.9.6-b1802` | Subsample assembly |
| Medaka / `medaka_consensus` | `2.2.0` | Consensus polishing |
| minimap2 | `2.30-r1287` | MLST allele alignment in the genotype helper |
| samtools | `1.23+htslib-1.23` | Alignment filtering, BAM generation, and sorting |
| dplyr | `1.1.4` | Plot summary transforms |
| ggplot2 | `3.5.2` | Figure 5 rendering |
| readr | `2.1.5` | Genotype-table input and summary output |
| tidyr | `1.3.1` | Table parsing and grid completion |
| patchwork | `1.3.2` | Optional combined-panel layout when available |

## 7. Workflow Outputs

The current staged output layout used by the final workflow is:
- `output/01_demultiplex`
- `output/02_assembly`
- `output/03_mlst`
- `output/04_st_guess`
- `output/Result_Table`

## 8. Notes

- The `macos-arm64` environment files are intended for Apple Silicon Macs, including M1, M2, M3, and M4 systems.
- The basecalling chain listed above is experiment-specific and may differ from newer ONT-recommended workflows.
- For manuscript writing, it is usually sufficient to keep the key platform, basecalling, and core analysis-tool versions in the Methods section and refer readers to this file for the full version summary.
