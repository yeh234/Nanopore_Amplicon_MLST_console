# Input / Output specification

## 01_prepare_public_ligation_data.R
Input:
- data/public_reference/06_matrix.csv
- data/public_reference/02_overhangs.csv

Output:
- data/processed/public_ligation_matrix.csv
- data/processed/public_ligation_long.csv
- data/processed/barcode24_candidate_long.csv
- output/tables/ligation_summary.tsv

## 02_screen_candidate_barcodes.R
Input:
- data/processed/public_ligation_long.csv
- output/tables/ligation_summary.tsv

Output:
- output/tables/barcode24_candidates.tsv
- output/tables/low_cross_ligation_candidates.tsv

## 03_plot_barcode24_heatmap.R
Input:
- data/processed/public_ligation_matrix.csv
- data/processed/barcode24_candidate_long.csv

Output:
- output/figures/barcode24_heatmap_clustered.png
- output/figures/barcode24_heatmap_clustered.pdf

## 04_generate_or_format_barcode108_panel.R
Input:
- data/reference/barcode108_panel.xlsx, or
- data/processed/barcode108_panel.tsv, or
- data/processed/barcode108_panel.csv, or
- `BARCODE_PANEL_WORKBOOK_PATH` to an Excel workbook
- data/reference/barcode15_reference.tsv

Output:
- data/processed/barcode108_panel.tsv
- data/processed/barcode108_panel.csv
- output/tables/barcode108_panel_long.tsv
- output/tables/barcode15_reference.tsv

## 05_export_primer_barcode_design.R
Input:
- data/processed/barcode108_panel.tsv
- data/reference/barcode15_reference.tsv
- data/reference/primer_slot_map.tsv

Output:
- output/tables/barcode108_primer_design.tsv
- output/tables/primer_slot_map.tsv
