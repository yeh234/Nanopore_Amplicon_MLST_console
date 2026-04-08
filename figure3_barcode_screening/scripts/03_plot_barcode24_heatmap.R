source("src/barcode_design_common.R")

build_barcode24_heatmap_plot <- function(plot_df, plot_title) {
  ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = .data$PairedSeq, y = .data$Overhang, fill = .data$Count)
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.35) +
    ggplot2::scale_fill_gradient(low = "white", high = "darkred") +
    ggplot2::labs(
      title = plot_title,
      x = "Paired Sequence",
      y = "Overhang Sequence",
      fill = "Count"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
      axis.text.y = ggplot2::element_text(size = 8),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    )
}

plot_barcode24_heatmap <- function() {
  barcode_design_require_packages(c("dplyr", "ggplot2", "readr"))
  barcode_design_ensure_output_dirs()

  matrix_path <- barcode_design_path("data", "processed", "public_ligation_matrix.csv")
  long_path <- barcode_design_path("data", "processed", "barcode24_candidate_long.csv")

  if (!file.exists(matrix_path) || !file.exists(long_path)) {
    stop("Run scripts/01_prepare_public_ligation_data.R before plotting the barcode24 heatmap.", call. = FALSE)
  }

  matrix_df <- readr::read_csv(matrix_path, show_col_types = FALSE)
  long_df <- readr::read_csv(long_path, show_col_types = FALSE)
  barcode24_reference <- barcode_design_default_barcode24_reference()

  clustered_order <- barcode_design_compute_heatmap_order(
    matrix_df = matrix_df,
    reference_overhangs = barcode24_reference$overhang
  )

  clustered_plot_df <- long_df |>
    dplyr::mutate(
      Overhang = factor(.data$Overhang, levels = clustered_order$row_order),
      PairedSeq = factor(.data$PairedSeq, levels = clustered_order$col_order)
    )

  clustered_heatmap_plot <- build_barcode24_heatmap_plot(
    plot_df = clustered_plot_df,
    plot_title = "Barcode24 Ligase Interaction Map (Clustered)"
  )

  ggplot2::ggsave(
    filename = barcode_design_path("output", "figures", "barcode24_heatmap_clustered.png"),
    plot = clustered_heatmap_plot,
    width = 10,
    height = 8,
    units = "in",
    dpi = 300
  )
  ggplot2::ggsave(
    filename = barcode_design_path("output", "figures", "barcode24_heatmap_clustered.pdf"),
    plot = clustered_heatmap_plot,
    width = 10,
    height = 8,
    units = "in"
  )

  invisible(clustered_heatmap_plot)
}

if (barcode_design_should_autorun()) {
  plot_barcode24_heatmap()
}
