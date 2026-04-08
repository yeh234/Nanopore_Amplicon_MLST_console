# Figure3 barcode screening

This repository provides a stepwise R workflow for:

- preparing published 4-bp ligation-fidelity data
- screening barcode candidates against the barcode24 set
- plotting a barcode24 interaction heatmap
- formatting the barcode108 panel workbook
- exporting primer-barcode design tables

This repository supports Golden Gate barcode design. It does not include the
downstream nanopore basecalling or MLST calling pipeline.

## Repository Layout

- `run_all.R`
  Wrapper that sources steps 1-5 from the repository root.
- `scripts/`
  Main executable scripts for each step.
- `src/`
  Shared constants and helper functions.
- `data/public_reference/`
  Location for `06_matrix.csv` and `02_overhangs.csv`.
- `data/reference/`
  Barcode and primer mapping tables used in steps 4-5, including
  `barcode108_panel.xlsx`.
- `data/processed/`
  Processed panel and ligation tables.
- `output/`
  Exported figures and summary tables.

## Required Inputs

This repository does not redistribute the published ligation reference files.
Before running steps 1-3, prepare:

- `data/public_reference/06_matrix.csv`
- `data/public_reference/02_overhangs.csv`

See `data/public_reference/README.md` for file origin and placement notes.

If the public reference files are stored elsewhere, set environment variables
first:

```r
Sys.setenv(
  BARCODE_LIGATION_MATRIX_PATH = "/path/to/06_matrix.csv",
  BARCODE_OVERHANGS_PATH = "/path/to/02_overhangs.csv"
)
```

Step 4 uses `data/reference/barcode108_panel.xlsx` by default. If that file is
not present, it falls back to `data/processed/barcode108_panel.tsv` or
`data/processed/barcode108_panel.csv`.

If you want to use a different Excel workbook, set:

```r
Sys.setenv(BARCODE_PANEL_WORKBOOK_PATH = "/path/to/barcode108_panel.xlsx")
```

If you need to specify a different sheet in the workbook, set:

```r
Sys.setenv(BARCODE_PANEL_WORKBOOK_SHEET = "sheet_name")
```

## Run

From the repository root, source the scripts in order:

```r
source("scripts/01_prepare_public_ligation_data.R")
source("scripts/02_screen_candidate_barcodes.R")
source("scripts/03_plot_barcode24_heatmap.R")
source("scripts/04_generate_or_format_barcode108_panel.R")
source("scripts/05_export_primer_barcode_design.R")
```

Or run the wrapper:

```r
source("run_all.R")
```

Required R packages:

```r
install.packages(c("dplyr", "ggplot2", "readr", "readxl", "tidyr"))
```

## Script Outputs

`01_prepare_public_ligation_data.R`
- writes `data/processed/public_ligation_matrix.csv`
- writes `data/processed/public_ligation_long.csv`
- writes `data/processed/barcode24_candidate_long.csv`
- writes `output/tables/ligation_summary.tsv`

`02_screen_candidate_barcodes.R`
- writes `output/tables/barcode24_candidates.tsv`
- writes `output/tables/low_cross_ligation_candidates.tsv`

`03_plot_barcode24_heatmap.R`
- writes `output/figures/barcode24_heatmap_clustered.png`
- writes `output/figures/barcode24_heatmap_clustered.pdf`

`04_generate_or_format_barcode108_panel.R`
- writes `data/processed/barcode108_panel.tsv`
- writes `data/processed/barcode108_panel.csv`
- writes `output/tables/barcode108_panel_long.tsv`
- writes `output/tables/barcode15_reference.tsv`

`05_export_primer_barcode_design.R`
- writes `output/tables/barcode108_primer_design.tsv`
- writes `output/tables/primer_slot_map.tsv`
