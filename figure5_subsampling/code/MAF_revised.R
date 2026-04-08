suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
})

# =========================
# 1) configurable inputs
# =========================
input_file <- Sys.getenv("INPUT_FILE", unset = "sample4_patch_genotype_result.txt")
output_prefix <- Sys.getenv("OUTPUT_PREFIX", unset = "Figure5_redrawn")
summary_file <- Sys.getenv("SUMMARY_FILE", unset = paste0(output_prefix, "_summary.tsv"))
expected_replicates_env <- Sys.getenv("EXPECTED_REPLICATES", unset = "")
zoom_read_min_env <- Sys.getenv("LOW_READ_ZOOM_MIN", unset = "5")
zoom_read_max_env <- Sys.getenv("LOW_READ_ZOOM_MAX", unset = "15")
guide_read_env <- Sys.getenv("GUIDE_READ", unset = "39")
guide_read_label <- Sys.getenv("GUIDE_READ_LABEL", unset = "")
guide_read_secondary_env <- Sys.getenv("GUIDE_READ_SECONDARY", unset = "")
guide_read_secondary_label <- Sys.getenv("GUIDE_READ_SECONDARY_LABEL", unset = "")
force_label_reads_env <- Sys.getenv("FORCE_LABEL_READS", unset = "39,60")
label_read_min_env <- Sys.getenv("LABEL_READ_MIN", unset = "5")
label_read_max_env <- Sys.getenv("LABEL_READ_MAX", unset = "15")
x_axis_label <- Sys.getenv("X_AXIS_LABEL", unset = "Subsampled reads")
y_axis_label <- Sys.getenv("Y_AXIS_LABEL", unset = "Locus-level correctness")
y_axis_min_env <- Sys.getenv("Y_AXIS_MIN", unset = "0")
y_axis_max_env <- Sys.getenv("Y_AXIS_MAX", unset = "1")

# =========================
# 2) expected allele profile for sample 4
# =========================
expected_profile <- tribble(
  ~Gene,   ~Type,
  "arcC",  "2",
  "aroE",  "3",
  "glpF",  "1",
  "gmk",   "1",
  "pta",   "4",
  "tpi",   "4",
  "yqiL",  "3"
)

gene_order <- expected_profile$Gene
gene_palette <- c(
  arcC = "#dc6b61",
  aroE = "#b7941f",
  glpF = "#71b63a",
  gmk = "#34b7a0",
  pta  = "#3fb7e8",
  tpi  = "#b0a0f5",
  yqiL = "#ef75d1"
)

stop_with_message <- function(...) {
  stop(paste0(...), call. = FALSE)
}

parse_int_list <- function(value) {
  if (!nzchar(value)) {
    return(integer())
  }

  out <- stringr::str_split(value, "[^0-9]+", simplify = FALSE)[[1]] %>%
    .[nzchar(.)] %>%
    as.integer()

  out[!is.na(out)]
}

parse_gene_type <- function(df) {
  df %>%
    tidyr::extract(
      col = Gene_type,
      into = c("Gene", "Type"),
      regex = "^(.*)_([^_]+)$",
      remove = FALSE
    )
}

parse_sample_source <- function(df) {
  df %>%
    tidyr::extract(
      col = SampleSource,
      into = c("Sample", "Reads", "Num"),
      regex = "^(.*)_([0-9]+)_([0-9]+)$",
      remove = FALSE
    )
}

if (!file.exists(input_file)) {
  stop_with_message("Input file not found: ", input_file)
}

zoom_read_min <- suppressWarnings(as.integer(zoom_read_min_env))
zoom_read_max <- suppressWarnings(as.integer(zoom_read_max_env))
guide_read <- suppressWarnings(as.integer(guide_read_env))
guide_read_secondary <- suppressWarnings(as.integer(guide_read_secondary_env))
force_label_reads <- parse_int_list(force_label_reads_env)
label_read_min <- suppressWarnings(as.integer(label_read_min_env))
label_read_max <- suppressWarnings(as.integer(label_read_max_env))
y_axis_min <- suppressWarnings(as.numeric(y_axis_min_env))
y_axis_max <- suppressWarnings(as.numeric(y_axis_max_env))

if (is.na(zoom_read_min) || is.na(zoom_read_max) || zoom_read_min >= zoom_read_max) {
  stop_with_message("LOW_READ_ZOOM_MIN and LOW_READ_ZOOM_MAX must be integers with MIN < MAX.")
}

if (nzchar(guide_read_env) && is.na(guide_read)) {
  stop_with_message("GUIDE_READ must be an integer.")
}

if (nzchar(guide_read_secondary_env) && is.na(guide_read_secondary)) {
  stop_with_message("GUIDE_READ_SECONDARY must be an integer.")
}

if (is.na(label_read_min) || is.na(label_read_max) || label_read_min > label_read_max) {
  stop_with_message("LABEL_READ_MIN and LABEL_READ_MAX must be integers with MIN <= MAX.")
}

if (is.na(y_axis_min) || is.na(y_axis_max) || y_axis_min >= y_axis_max) {
  stop_with_message("Y_AXIS_MIN and Y_AXIS_MAX must be numeric values with MIN < MAX.")
}

# =========================
# 3) read genotype result table
# =========================
geno_raw <- readr::read_tsv(
  file = input_file,
  col_names = FALSE,
  show_col_types = FALSE,
  progress = FALSE
) %>%
  rlang::set_names(c(
    "Gene_type",
    "SampleSource",
    "Cigar",
    "Tag_NM",
    "Tag_MS",
    "Tag_AS",
    "Tag_NN",
    "Tag_TP"
  )) %>%
  parse_gene_type() %>%
  parse_sample_source() %>%
  mutate(
    Reads = as.integer(Reads),
    Num = as.integer(Num),
    Cigar_check = !str_detect(Cigar, "S")
  )

if (nrow(geno_raw) == 0) {
  stop_with_message("Input file is empty: ", input_file)
}

bad_rows <- geno_raw %>%
  filter(is.na(Gene) | is.na(Type) | is.na(Reads) | is.na(Num))

if (nrow(bad_rows) > 0) {
  stop_with_message(
    "Failed to parse some rows. ",
    "Expected Gene_type like 'arcC_2' and SampleSource like 'sample4_25_1'."
  )
}

if (nzchar(expected_replicates_env)) {
  expected_replicates <- suppressWarnings(as.integer(expected_replicates_env))
  if (is.na(expected_replicates) || expected_replicates <= 0) {
    stop_with_message("EXPECTED_REPLICATES must be a positive integer.")
  }

  denominator_df <- geno_raw %>%
    distinct(Reads) %>%
    mutate(replicate_n = expected_replicates)
} else {
  denominator_df <- geno_raw %>%
    group_by(Reads) %>%
    summarise(replicate_n = max(Num), .groups = "drop")
}

read_levels <- denominator_df %>%
  pull(Reads) %>%
  sort()

n_read_levels <- length(read_levels)
min_reads <- min(read_levels)
max_reads <- max(read_levels)
dense_x_axis <- n_read_levels >= 50 && min_reads <= 5
zoom_levels <- read_levels[read_levels >= zoom_read_min & read_levels <= zoom_read_max]
has_zoom_panel <- length(zoom_levels) >= 4 && length(zoom_levels) < n_read_levels

if (n_read_levels <= 1) {
  x_breaks_major <- read_levels
} else if (dense_x_axis) {
  x_breaks_major <- unique(c(
    min_reads,
    read_levels[read_levels %% 5 == 0],
    max_reads
  )) %>% sort()
} else {
  x_breaks_major <- read_levels
}

x_breaks_minor <- if (n_read_levels > 1) read_levels else NULL
plot_width <- if (dense_x_axis) 15.5 else 13.5
x_axis_text_size <- if (dense_x_axis) 6 else 7
label_text_size <- if (dense_x_axis) 2.5 else 3.1
zoom_label_text_size <- 2.7
guide_label_text_size <- 3.4
main_plot_height <- 6.8
combined_plot_height <- if (has_zoom_panel) 10.2 else main_plot_height

# =========================
# 4) retain expected allele calls only
# =========================
correct_calls <- geno_raw %>%
  filter(Cigar_check) %>%
  inner_join(expected_profile, by = c("Gene", "Type")) %>%
  distinct(Reads, Num, Gene, Type)

# =========================
# 5) calculate locus-level correctness
# =========================
fig5_df <- tidyr::crossing(
  Reads = read_levels,
  expected_profile
) %>%
  left_join(
    correct_calls %>%
      count(Reads, Gene, Type, name = "correct_n"),
    by = c("Reads", "Gene", "Type")
  ) %>%
  left_join(denominator_df, by = "Reads") %>%
  mutate(
    correct_n = replace_na(correct_n, 0L),
    correctness = correct_n / replicate_n,
    Gene = factor(Gene, levels = gene_order)
  ) %>%
  group_by(Gene) %>%
  mutate(
    read_index = row_number(),
    label = percent(correctness, accuracy = 1),
    show_label = (Reads >= label_read_min & Reads <= label_read_max) | Reads %in% force_label_reads,
    label_y = pmin(
      correctness + if_else(read_index %% 2L == 0L, 0.035, 0.075),
      1.03
    )
  ) %>%
  ungroup() %>%
  arrange(Gene, Reads)

readr::write_tsv(fig5_df, summary_file)
print(fig5_df, n = nrow(fig5_df))

label_df <- fig5_df %>%
  filter(show_label)

zoom_df <- fig5_df %>%
  filter(Reads >= zoom_read_min, Reads <= zoom_read_max) %>%
  group_by(Gene) %>%
  mutate(
    zoom_index = row_number(),
    zoom_label_y = pmin(
      correctness + if_else(zoom_index %% 2L == 0L, 0.04, 0.08),
      1.04
    )
  ) %>%
  ungroup()

if (length(read_levels) == 1) {
  x_limits <- c(read_levels[1] - 0.5, read_levels[1] + 0.5)
} else {
  x_limits <- range(read_levels)
}

if (length(zoom_levels) == 1) {
  zoom_x_limits <- c(zoom_levels[1] - 0.5, zoom_levels[1] + 0.5)
} else if (length(zoom_levels) > 1) {
  zoom_x_limits <- range(zoom_levels)
} else {
  zoom_x_limits <- NULL
}

primary_guide_df <- tibble::tibble(
  x = if (!is.na(guide_read) && guide_read >= min_reads && guide_read <= max_reads) guide_read else numeric(),
  label = if (!is.na(guide_read) && nzchar(guide_read_label)) guide_read_label else character()
)

secondary_guide_df <- tibble::tibble(
  x = if (!is.na(guide_read_secondary) && guide_read_secondary >= min_reads && guide_read_secondary <= max_reads) guide_read_secondary else numeric(),
  label = if (!is.na(guide_read_secondary) && nzchar(guide_read_secondary_label)) guide_read_secondary_label else character()
)

top_facet_gene <- factor(gene_order[[1]], levels = gene_order)

primary_guide_label_df <- tibble::tibble(
  Reads = if (nrow(primary_guide_df) > 0) primary_guide_df$x else numeric(),
  correctness = if (nrow(primary_guide_df) > 0) y_axis_max + 0.17 else numeric(),
  label = if (nrow(primary_guide_df) > 0) primary_guide_df$label else character(),
  Gene = top_facet_gene,
  hjust = 0.5
)

secondary_guide_label_df <- tibble::tibble(
  Reads = if (nrow(secondary_guide_df) > 0) secondary_guide_df$x else numeric(),
  correctness = if (nrow(secondary_guide_df) > 0) y_axis_max + 0.17 else numeric(),
  label = if (nrow(secondary_guide_df) > 0) secondary_guide_df$label else character(),
  Gene = top_facet_gene,
  hjust = 0.5
)

# =========================
# 6) plot
# =========================
p_fig5 <- ggplot(
  fig5_df,
  aes(x = Reads, y = correctness, group = Gene, color = Gene)
) +
  geom_vline(
    data = secondary_guide_df,
    aes(xintercept = x),
    linetype = "dashed",
    color = "grey70",
    linewidth = 0.35
  ) +
  geom_vline(
    data = primary_guide_df,
    aes(xintercept = x),
    linetype = "dashed",
    color = "grey35",
    linewidth = 0.4
  ) +
  geom_point(size = 2.1) +
  geom_text(
    data = label_df,
    aes(y = label_y, label = label),
    color = "grey20",
    size = label_text_size,
    vjust = 0,
    check_overlap = TRUE,
    show.legend = FALSE
  ) +
  facet_grid(rows = vars(Gene)) +
  scale_color_manual(values = gene_palette, drop = FALSE) +
  scale_y_continuous(
    breaks = seq(y_axis_min, y_axis_max, by = 0.25),
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0.02, 0.12))
  ) +
  scale_x_continuous(
    breaks = x_breaks_major,
    minor_breaks = x_breaks_minor,
    limits = x_limits,
    expand = expansion(mult = c(0.005, 0.01))
  ) +
  labs(
    x = x_axis_label,
    y = y_axis_label
  ) +
  coord_cartesian(ylim = c(y_axis_min, y_axis_max), clip = "off") +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "grey88", color = NA),
    strip.text.y = element_text(face = "bold", size = 9.6),
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_line(
      color = if (dense_x_axis) "grey97" else NA,
      linewidth = if (dense_x_axis) 0.15 else 0
    ),
    panel.grid.major.x = element_line(color = "grey92", linewidth = 0.22),
    panel.grid.major.y = element_line(color = "grey92", linewidth = 0.22),
    panel.spacing.y = unit(1.2, "mm"),
    legend.position = "none",
    axis.text.x = element_text(size = x_axis_text_size),
    axis.text.y = element_text(size = 8),
    axis.title = element_text(face = "plain"),
    plot.margin = margin(34, 18, 8, 8)
  )

if (length(read_levels) > 1) {
  p_fig5 <- p_fig5 + geom_line(linewidth = 0.95)
}

if (nrow(primary_guide_label_df) > 0 && nzchar(primary_guide_label_df$label[[1]])) {
  p_fig5 <- p_fig5 +
    geom_text(
      data = primary_guide_label_df,
      aes(x = Reads, y = correctness, label = label, hjust = hjust),
      inherit.aes = FALSE,
      vjust = 0,
      size = guide_label_text_size,
      color = "grey10"
    )
}

if (nrow(secondary_guide_label_df) > 0 && nzchar(secondary_guide_label_df$label[[1]])) {
  p_fig5 <- p_fig5 +
    geom_text(
      data = secondary_guide_label_df,
      aes(x = Reads, y = correctness, label = label, hjust = hjust),
      inherit.aes = FALSE,
      vjust = 0,
      size = guide_label_text_size,
      color = "grey15"
    )
}

p_zoom <- NULL
if (has_zoom_panel) {
  p_zoom <- ggplot(
    zoom_df,
    aes(x = Reads, y = correctness, group = Gene, color = Gene)
  ) +
    geom_vline(
      data = secondary_guide_df %>%
        dplyr::filter(x >= zoom_x_limits[1], x <= zoom_x_limits[2]),
      aes(xintercept = x),
      linetype = "dashed",
      color = "grey70",
      linewidth = 0.35
    ) +
    geom_vline(
      data = primary_guide_df %>%
        dplyr::filter(x >= zoom_x_limits[1], x <= zoom_x_limits[2]),
      aes(xintercept = x),
      linetype = "dashed",
      color = "grey45",
      linewidth = 0.4
    ) +
    geom_point(size = 2.1) +
    geom_text(
      aes(y = zoom_label_y, label = label),
      color = "grey20",
      size = zoom_label_text_size,
      vjust = 0,
      check_overlap = TRUE,
      show.legend = FALSE
    ) +
    facet_grid(rows = vars(Gene)) +
    scale_color_manual(values = gene_palette, drop = FALSE) +
    scale_y_continuous(
      breaks = seq(y_axis_min, y_axis_max, by = 0.25),
      labels = percent_format(accuracy = 1),
      expand = expansion(mult = c(0.02, 0.12))
    ) +
    scale_x_continuous(
      breaks = zoom_levels,
      minor_breaks = zoom_levels,
      limits = zoom_x_limits,
      expand = expansion(mult = c(0.01, 0.02))
    ) +
    labs(
      x = x_axis_label,
      y = y_axis_label,
      title = sprintf("Low-read zoom: %d-%d reads", zoom_read_min, zoom_read_max)
    ) +
    coord_cartesian(ylim = c(y_axis_min, y_axis_max), clip = "off") +
    theme_bw(base_size = 11) +
    theme(
      strip.background = element_rect(fill = "grey88", color = NA),
      strip.text.y = element_text(face = "bold", size = 9.6),
      panel.grid.minor.y = element_blank(),
      panel.grid.minor.x = element_line(color = "grey97", linewidth = 0.15),
      panel.grid.major.x = element_line(color = "grey92", linewidth = 0.22),
      panel.grid.major.y = element_line(color = "grey92", linewidth = 0.22),
      panel.spacing.y = unit(1.2, "mm"),
      legend.position = "none",
      plot.title = element_text(size = 10, face = "plain"),
      axis.text.x = element_text(size = 7),
      axis.text.y = element_text(size = 8),
      axis.title = element_text(face = "plain"),
      plot.margin = margin(4, 18, 8, 8)
    )

  if (length(zoom_levels) > 1) {
    p_zoom <- p_zoom + geom_line(linewidth = 0.95)
  }
}

plot_object <- p_fig5
if (has_zoom_panel && requireNamespace("patchwork", quietly = TRUE)) {
  plot_object <- patchwork::wrap_plots(
    p_fig5 + theme(plot.margin = margin(8, 18, 2, 8)),
    p_zoom,
    ncol = 1,
    heights = c(3.4, 2.2)
  )
}

print(plot_object)

# =========================
# 7) save
# =========================
ggsave(
  filename = paste0(output_prefix, ".png"),
  plot = plot_object,
  width = plot_width,
  height = combined_plot_height,
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = paste0(output_prefix, ".pdf"),
  plot = plot_object,
  width = plot_width,
  height = combined_plot_height,
  bg = "white"
)

if (has_zoom_panel) {
  ggsave(
    filename = paste0(output_prefix, "_lowread_zoom.png"),
    plot = p_zoom,
    width = plot_width,
    height = 5.8,
    dpi = 600,
    bg = "white"
  )

  ggsave(
    filename = paste0(output_prefix, "_lowread_zoom.pdf"),
    plot = p_zoom,
    width = plot_width,
    height = 5.8,
    bg = "white"
  )
}

message("Summary written to: ", summary_file)
message("Figure written to: ", paste0(output_prefix, ".png"))
message("Figure written to: ", paste0(output_prefix, ".pdf"))
if (has_zoom_panel) {
  message("Low-read zoom written to: ", paste0(output_prefix, "_lowread_zoom.png"))
  message("Low-read zoom written to: ", paste0(output_prefix, "_lowread_zoom.pdf"))
}
