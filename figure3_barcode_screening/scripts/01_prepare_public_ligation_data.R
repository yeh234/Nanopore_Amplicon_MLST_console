source("src/barcode_design_common.R")

prepare_public_ligation_data <- function() {
  barcode_design_require_packages(c("dplyr", "readr", "tidyr"))
  barcode_design_ensure_output_dirs()

  matrix_df <- barcode_design_read_public_ligation_matrix()
  barcode_design_validate_overhang_reference(matrix_df)

  long_df <- barcode_design_make_ligation_long(matrix_df)
  summary_df <- barcode_design_summarise_fidelity(long_df)
  barcode24_reference <- barcode_design_default_barcode24_reference()

  barcode24_long_df <- long_df |>
    dplyr::filter(
      .data$Overhang %in% barcode24_reference$overhang,
      .data$PairedSeq %in% barcode24_reference$overhang
    )

  readr::write_csv(
    matrix_df,
    barcode_design_path("data", "processed", "public_ligation_matrix.csv")
  )
  readr::write_csv(
    long_df,
    barcode_design_path("data", "processed", "public_ligation_long.csv")
  )
  readr::write_csv(
    barcode24_long_df,
    barcode_design_path("data", "processed", "barcode24_candidate_long.csv")
  )
  readr::write_tsv(
    summary_df,
    barcode_design_path("output", "tables", "ligation_summary.tsv")
  )

  invisible(
    list(
      matrix = matrix_df,
      long = long_df,
      summary = summary_df,
      barcode24_long = barcode24_long_df
    )
  )
}

if (barcode_design_should_autorun()) {
  prepare_public_ligation_data()
}
