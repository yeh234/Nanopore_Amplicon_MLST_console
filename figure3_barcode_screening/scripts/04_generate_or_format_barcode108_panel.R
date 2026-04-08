source("src/barcode_design_common.R")

generate_or_format_barcode108_panel <- function() {
  barcode_design_require_packages(c("readr"))
  barcode_design_ensure_output_dirs()

  panel_df <- barcode_design_read_panel_source() |>
    barcode_design_standardize_panel_table()

  barcode15_reference <- barcode_design_load_barcode15_reference()
  panel_long_df <- barcode_design_panel_to_long(
    panel_df = panel_df,
    barcode_reference_df = barcode15_reference
  )

  readr::write_tsv(
    panel_df,
    barcode_design_path("data", "processed", "barcode108_panel.tsv")
  )
  readr::write_csv(
    panel_df,
    barcode_design_path("data", "processed", "barcode108_panel.csv")
  )
  readr::write_tsv(
    panel_long_df,
    barcode_design_path("output", "tables", "barcode108_panel_long.tsv")
  )
  readr::write_tsv(
    barcode15_reference,
    barcode_design_path("output", "tables", "barcode15_reference.tsv")
  )

  invisible(
    list(
      panel = panel_df,
      panel_long = panel_long_df,
      barcode15_reference = barcode15_reference
    )
  )
}

if (barcode_design_should_autorun()) {
  generate_or_format_barcode108_panel()
}
