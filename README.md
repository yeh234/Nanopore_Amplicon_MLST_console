# Nanopore Amplicon *S. aureus* MLST Workflow

This repository contains the final workflow used to process **Oxford Nanopore Dorado FASTQ** data for within-read barcode demultiplexing, per-sample assembly and polish, standard MLST calling, optional ST inference, and final publication-style MLST table generation.

## Repository Scope

This repository provides the command-line workflow, a local Shiny wrapper,
lightweight public test assets, and figure-specific reproductions used in this
study.

- Full study FASTQ files are not expected to be committed to the public repository.
- PubMLST reference files may be supplied locally under `db/` or regenerated from the official PubMLST sources.
- Large local-only files are excluded through the root `.gitignore`.

## Repository layout

- `code/`: workflow scripts
- `app.R`: local Shiny wrapper for Step 1-4
- `envs/`: conda environment YAML files
- `demo/`: smaller demo FASTQ for lightweight testing
- `fixtures/`: downstream Step 2 fixture package for exact public testing
- `fastq/`: local raw FASTQ working directory
- `db/`: local PubMLST allele/profile snapshot used during local validation
- `figure3_barcode_screening/`: standalone Figure 3 barcode-screening workflow and exports
- `figure5_subsampling/`: standalone Figure 5 subsampling workflow; reuses the repository-root `db/` by default
- `OT-2_protocols/`: curated OT-2 notebooks and layout files used by this study; see [OT-2_protocols/README.md](OT-2_protocols/README.md)

## Data availability

### Full study FASTQ

The 128,431 Dorado SUP reads used for the study were deposited in the NCBI
Sequence Read Archive under BioProject
[PRJNA1561677](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1561677).
They were re-basecalled from the retained pass FAST5 files and are the input
before length filtering and demultiplexing.

Retrieve the FASTQ from SRA and place it under `fastq/` before starting Option A.
See [the FASTQ input notes](fastq/README.md) for the expected layout.

### Public test assets

Two public test assets are included for lightweight verification:

- `demo/Samlst_SUP.demo.fastq.gz`: subset for a demultiplexing smoke test
- `fixtures/downstream_step2_fixture.tar.gz`: fixed per-sample reads and expected downstream outputs

The fixture's validator compares results against its own reference outputs.
Table 1 reports read counts from the full study analysis.

## PubMLST reference data

The workflow expects local PubMLST allele and profile files under `db/` when
running the ST-guess step and figure-specific rebuilds that depend on MLST
references.

Default workflow note:

- Step 3 uses `db/saureus_Profile.txt`

Expected local filenames:

- `db/arcC.fas`
- `db/aroE.fas`
- `db/glpF.fas`
- `db/gmk.fas`
- `db/pta.fas`
- `db/tpi.fas`
- `db/yqiL.fas`
- `db/saureus_Profile.txt`

If these files are not already present locally, regenerate them from the
official PubMLST sources. The committed `db/README.md` documents the expected
layout.

Practical reconstruction guidance:

- open the official PubMLST MLST download page
- in the `Staphylococcus aureus` row, download `profiles`, `arcC`, `aroE`, `glpF`, `gmk`, `pta`, `tpi`, and `yqiL`
- store them locally as:
  - `db/saureus_Profile.txt`
  - `db/arcC.fas`
  - `db/aroE.fas`
  - `db/glpF.fas`
  - `db/gmk.fas`
  - `db/pta.fas`
  - `db/tpi.fas`
  - `db/yqiL.fas`
- for programmatic profile download, the PubMLST REST endpoint for the *S. aureus* MLST profile table is:
  - `https://rest.pubmlst.org/db/pubmlst_saureus_seqdef/schemes/1/profiles_csv`

Official PubMLST references:

- API overview: <https://pubmlst.org/api>
- API base URL: <https://rest.pubmlst.org/>
- general MLST data downloads page: <https://dev.pubmlst.org/data>
- terms and conditions: <https://pubmlst.org/terms-conditions>

Important licensing note:

- PubMLST states that data submitted on or after 1 January 2025 are subject to additional restrictions on download and redistribution.
- Review the current PubMLST terms before redistributing any database snapshot in a public repository.

Required acknowledgement for publications using PubMLST data:

- cite the PubMLST website and Jolley et al. 2018 as requested on the PubMLST site

## Conda environments

The repository uses the following environment files:

- `envs/macos-arm64-nanopore-assembly-mlst.yml` -> `nanopore-assembly-mlst-macos-arm64`
- `envs/macos-arm64-mlst-only.yml` -> `mlst-only`

Create them once:

```bash
conda env create -f envs/macos-arm64-nanopore-assembly-mlst.yml
conda env create -f envs/macos-arm64-mlst-only.yml
```

If you update the YAML files later:

```bash
conda env update -f envs/macos-arm64-nanopore-assembly-mlst.yml --prune
conda env update -f envs/macos-arm64-mlst-only.yml --prune
```

Quick checks:

```bash
conda run -n nanopore-assembly-mlst-macos-arm64 flye --version
conda run -n nanopore-assembly-mlst-macos-arm64 bwa 2>&1 | head -n 1
conda run -n nanopore-assembly-mlst-macos-arm64 samtools --version | head -n 1
conda run -n nanopore-assembly-mlst-macos-arm64 medaka --help > /dev/null
conda run -n mlst-only mlst --version
```

When using `app.R`, keep the Step 2 fields `assembly conda env` and `MLST conda env` aligned with these exact environment names. The environment diagnostics card checks:

- the current shell for `Rscript`, `python3`, `bash`, `conda`, `seqkit`, and `bioawk`
- the assembly env for `flye`, `bwa`, `samtools`, `medaka_consensus`, and `medaka`
- the MLST env for `mlst`

## Start options

### Option A. Full local run from a Dorado FASTQ

Step 1. Demultiplex from Dorado FASTQ

```bash
Rscript code/01_sup_demultiplex.R \
  --input-dir=fastq \
  --output-dir=output/01_demultiplex
```

Step 2. Assembly, polish, and standard MLST

```bash
bash code/04_assembly_polish_pipeline.sh \
  --sample-fastq-dir=output/01_demultiplex/merged_by_sample \
  --assembly-output-dir=output/02_assembly \
  --mlst-output-dir=output/03_mlst \
  --assembly-env=nanopore-assembly-mlst-macos-arm64 \
  --mlst-env=mlst-only \
  --threads=4 \
  --mlst-threads=1 \
  --asm-coverage=40 \
  --min-overlap=1000 \
  --medaka-model=r941_min_sup_g507
```

Step 3. Optional ST guess table

```bash
python3 code/05_infer_st_from_partial_mlst.py \
  --mlst-summary=output/03_mlst/mlst_sample_summary.tsv \
  --profile-table=db/saureus_Profile.txt \
  --output-dir=output/04_st_guess
```

Step 4. Build publication-style MLST table

```bash
python3 code/06_build_publication_style_mlst_table.py \
  --assembly-summary=output/02_assembly/assembly_polish_summary.tsv \
  --mlst-summary=output/03_mlst/mlst_sample_summary.tsv \
  --st-guess-summary=output/04_st_guess/st_guess_summary.tsv \
  --output-path=output/Result_Table/SUP_final_publication_style_mlst.tsv
```

### Option B. Demo from Step 1

If you want a smaller Step 1 test input instead of the full raw FASTQ, use:

- input: `demo/Samlst_SUP.demo.fastq.gz`
- notes: `demo/README.md`

Example Step 1-only smoke run:

```bash
mkdir -p output/demo_step1_input
cp demo/Samlst_SUP.demo.fastq.gz output/demo_step1_input/

Rscript code/01_sup_demultiplex.R \
  --input-dir=output/demo_step1_input \
  --output-dir=output/demo_step1_demultiplex
```

You can then continue with the same Step 2-4 commands shown above, replacing the Step 1 output directory with `output/demo_step1_demultiplex`.

### Option C. Public Step 2 fixture run

Extract the downstream fixture from the repository root:

```bash
mkdir -p output/fixture_inputs
tar -xzf fixtures/downstream_step2_fixture.tar.gz -C output/fixture_inputs
```

Follow `output/fixture_inputs/downstream_step2_fixture/README.md`. This route
starts from fixed `merged_by_sample/*.fastq.gz` inputs and checks the Step 2-4
outputs with the included validator.

## Shiny app

If you want a local GUI wrapper for the workflow:

```bash
Rscript app.R
```

This direct entrypoint attempts to open your default browser automatically. If browser launch is blocked, use the local URL printed in the terminal.

or:

```bash
R -e "install.packages('shiny')"
R -e "shiny::runApp('.')"
```

The Shiny app is only a wrapper. All CLI tools and conda environments still need to be installed first.

## Expected outputs

- `output/01_demultiplex/file_summary.tsv`
- `output/02_assembly/assembly_polish_summary.tsv`
- `output/03_mlst/mlst_sample_summary.tsv`
- `output/04_st_guess/st_guess_summary.tsv`
- `output/Result_Table/SUP_final_publication_style_mlst.tsv`

## Validated downstream settings

For the current *S. aureus* study dataset, the validated downstream settings were:

- `Flye --asm-coverage=40`
- `Flye --min-overlap=1000`
- `medaka model = r941_min_sup_g507`
