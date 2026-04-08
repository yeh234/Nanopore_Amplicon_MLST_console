source("src/barcode_design_common.R")

export_primer_barcode_design <- function() {
  barcode_design_require_packages(c("readr"))
  barcode_design_ensure_output_dirs()

  panel_path <- barcode_design_path("data", "processed", "barcode108_panel.tsv")
  if (!file.exists(panel_path)) {
    stop("Run scripts/04_generate_or_format_barcode108_panel.R before exporting primer-barcode designs.", call. = FALSE)
  }

  panel_df <- readr::read_tsv(panel_path, show_col_types = FALSE)
  barcode15_reference <- barcode_design_load_barcode15_reference()
  primer_slot_map <- barcode_design_load_primer_slot_map()

  panel_long_df <- barcode_design_panel_to_long(
    panel_df = panel_df,
    barcode_reference_df = barcode15_reference
  )

  primer_design_df <- barcode_design_build_primer_design(
    panel_long_df = panel_long_df,
    primer_slot_map_df = primer_slot_map
  )

  readr::write_tsv(
    primer_design_df,
    barcode_design_path("output", "tables", "barcode108_primer_design.tsv")
  )
  readr::write_tsv(
    primer_slot_map,
    barcode_design_path("output", "tables", "primer_slot_map.tsv")
  )

  invisible(primer_design_df)
}

if (barcode_design_should_autorun()) {
  export_primer_barcode_design()
}
