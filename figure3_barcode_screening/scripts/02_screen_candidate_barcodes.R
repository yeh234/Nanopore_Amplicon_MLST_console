source("src/barcode_design_common.R")

screen_candidate_barcodes <- function() {
  barcode_design_require_packages(c("dplyr", "readr"))
  barcode_design_ensure_output_dirs()

  long_path <- barcode_design_path("data", "processed", "public_ligation_long.csv")
  summary_path <- barcode_design_path("output", "tables", "ligation_summary.tsv")

  if (!file.exists(long_path) || !file.exists(summary_path)) {
    stop("Run scripts/01_prepare_public_ligation_data.R before screening candidate barcodes.", call. = FALSE)
  }

  long_df <- readr::read_csv(long_path, show_col_types = FALSE)
  summary_df <- readr::read_tsv(summary_path, show_col_types = FALSE)
  barcode24_reference <- barcode_design_default_barcode24_reference()
  barcode24_reference$legacy_order <- seq_len(nrow(barcode24_reference))

  barcode24_candidates <- barcode24_reference |>
    dplyr::left_join(summary_df, by = c("overhang" = "Overhang"))

  low_cross_ligation_candidates <- barcode_design_score_candidate_cross_ligation(
    long_df = long_df,
    summary_df = summary_df,
    reference_overhangs = barcode24_reference$overhang
  )

  readr::write_tsv(
    barcode24_candidates,
    barcode_design_path("output", "tables", "barcode24_candidates.tsv")
  )
  readr::write_tsv(
    low_cross_ligation_candidates,
    barcode_design_path("output", "tables", "low_cross_ligation_candidates.tsv")
  )

  invisible(
    list(
      barcode24_candidates = barcode24_candidates,
      low_cross_ligation_candidates = low_cross_ligation_candidates
    )
  )
}

if (barcode_design_should_autorun()) {
  screen_candidate_barcodes()
}
