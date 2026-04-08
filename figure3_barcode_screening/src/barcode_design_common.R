barcode_design_best_barcode_space <- c(
  "AATG", "AACA", "GACC", "GGGA", "CGCC", "ACCT",
  "AGAC", "CAGA", "GTTA", "ACTA", "CTAA", "CCGC"
)

barcode_design_current_barcode_mlst <- c(
  "CATT", "TGTT", "GGTC", "GGTC", "TCCC", "GGCG",
  "AGGT", "GTCT", "TCTG", "TAAC", "TAGT", "TTAG", "GCGG"
)

barcode_design_bsaI_recognition_sequence <- "GGTCTC"
barcode_design_left_primer_prefix <- "A"
barcode_design_right_primer_prefix <- "T"

barcode_design_default_barcode24_reference <- function() {
  overhang <- unique(c(
    barcode_design_best_barcode_space,
    barcode_design_current_barcode_mlst
  ))

  source <- ifelse(
    overhang %in% barcode_design_best_barcode_space &
      overhang %in% unique(barcode_design_current_barcode_mlst),
    "legacy_best_barcode_space_and_current_mlst",
    ifelse(
      overhang %in% barcode_design_best_barcode_space,
      "legacy_best_barcode_space",
      "legacy_current_mlst"
    )
  )

  data.frame(
    overhang = overhang,
    source = source,
    stringsAsFactors = FALSE
  )
}

barcode_design_default_barcode15_reference <- function() {
  data.frame(
    barcode_id = 1:15,
    barcode_sequence = c(
      barcode_design_best_barcode_space,
      "ATGA", "ATGC", "CCTA"
    ),
    slot = c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 1, 2, 5),
    source = c(
      rep("legacy_best_barcode_space", 12),
      rep("legacy_screened_extension", 3)
    ),
    note = c(
      rep(
        "Barcode IDs 1-12 follow the best.barcode.space order from the archive script.",
        12
      ),
      "Barcode 13 is explicit in the archive primer block; slot assignment is inferred from the 108-panel workbook.",
      "Barcode 14 is explicit in the archive primer block; slot assignment is inferred from the 108-panel workbook.",
      "Barcode 15 is explicit in the archive primer block; slot assignment is inferred from the 108-panel workbook."
    ),
    stringsAsFactors = FALSE
  )
}

barcode_design_default_primer_slot_map <- function() {
  data.frame(
    slot = 1:6,
    left_primer_name = c("arcC_R", "aroE_R", NA, NA, "pta_R", NA),
    left_primer_seq = c(
      "AGGTATCTGCTTCAATCAGCG",
      "GGTGTTGTATTAATAACGATATC",
      NA, NA,
      "GACCCTTTTGTTGAAAAGCTTAA",
      NA
    ),
    right_primer_name = c("aroE_F", "glpF_F", NA, NA, "tpi_F", NA),
    right_primer_seq = c(
      "ATCGGAAATCCTATTTCACATTC",
      "CTAGGAACTGCAATCTTAATCC",
      NA, NA,
      "TCGTTCATTCTGAACGTCGTGAA",
      NA
    ),
    status = c("confirmed", "confirmed", "missing", "missing", "confirmed", "missing"),
    note = c(
      "Explicitly shown in the archive primer-design example for barcode 13.",
      "Explicitly shown in the archive primer-design example for barcode 14.",
      "Primer anchors for slot 3 are not present in the archive script and must be filled manually.",
      "Primer anchors for slot 4 are not present in the archive script and must be filled manually.",
      "Explicitly shown in the archive primer-design example for barcode 15.",
      "Primer anchors for slot 6 are not present in the archive script and must be filled manually."
    ),
    stringsAsFactors = FALSE
  )
}

barcode_design_collect_root_candidates <- function() {
  frame_paths <- vapply(
    sys.frames(),
    function(frame) {
      ofile <- frame$ofile
      if (is.null(ofile) || !nzchar(ofile)) {
        return("")
      }
      dirname(normalizePath(ofile, winslash = "/", mustWork = FALSE))
    },
    FUN.VALUE = character(1)
  )

  unique(c(
    Sys.getenv("BARCODE_DESIGN_ROOT", unset = ""),
    frame_paths,
    getwd()
  ))
}

barcode_design_find_root <- function() {
  markers <- c("README.md", "run_all.R", file.path("src", "barcode_design_common.R"))

  for (candidate in barcode_design_collect_root_candidates()) {
    if (!nzchar(candidate)) {
      next
    }

    current <- normalizePath(candidate, winslash = "/", mustWork = FALSE)
    if (file.exists(current) && !dir.exists(current)) {
      current <- dirname(current)
    }

    repeat {
      if (all(file.exists(file.path(current, markers)))) {
        return(current)
      }

      parent <- dirname(current)
      if (identical(parent, current)) {
        break
      }
      current <- parent
    }
  }

  stop(
    "Could not locate the project root. Set BARCODE_DESIGN_ROOT before sourcing the scripts.",
    call. = FALSE
  )
}

barcode_design_root <- barcode_design_find_root()

barcode_design_path <- function(...) {
  file.path(barcode_design_root, ...)
}

barcode_design_should_autorun <- function() {
  !tolower(Sys.getenv("BARCODE_DESIGN_AUTORUN", unset = "true")) %in% c("0", "false", "no")
}

barcode_design_require_packages <- function(packages) {
  missing_packages <- packages[
    !vapply(packages, requireNamespace, FUN.VALUE = logical(1), quietly = TRUE)
  ]

  if (length(missing_packages) > 0) {
    stop(
      sprintf(
        "Missing required package(s): %s",
        paste(missing_packages, collapse = ", ")
      ),
      call. = FALSE
    )
  }
}

barcode_design_ensure_output_dirs <- function() {
  dirs <- c(
    barcode_design_path("data", "processed"),
    barcode_design_path("output", "figures"),
    barcode_design_path("output", "tables")
  )

  invisible(vapply(
    dirs,
    function(dir_path) {
      dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
      TRUE
    },
    FUN.VALUE = logical(1)
  ))
}

barcode_design_get_public_matrix_path <- function() {
  override <- Sys.getenv("BARCODE_LIGATION_MATRIX_PATH", unset = "")
  if (nzchar(override)) {
    return(normalizePath(override, winslash = "/", mustWork = FALSE))
  }
  barcode_design_path("data", "public_reference", "06_matrix.csv")
}

barcode_design_get_overhang_reference_path <- function() {
  override <- Sys.getenv("BARCODE_OVERHANGS_PATH", unset = "")
  if (nzchar(override)) {
    return(normalizePath(override, winslash = "/", mustWork = FALSE))
  }
  barcode_design_path("data", "public_reference", "02_overhangs.csv")
}

barcode_design_get_panel_workbook_path <- function() {
  override <- Sys.getenv("BARCODE_PANEL_WORKBOOK_PATH", unset = "")
  if (nzchar(override)) {
    return(normalizePath(override, winslash = "/", mustWork = FALSE))
  }
  default_path <- barcode_design_path("data", "reference", "barcode108_panel.xlsx")
  if (file.exists(default_path)) {
    return(default_path)
  }
  ""
}

barcode_design_get_panel_sheet <- function() {
  override <- Sys.getenv("BARCODE_PANEL_WORKBOOK_SHEET", unset = "")
  if (nzchar(override)) {
    return(override)
  }
  "20220723_barcode108_combination"
}

barcode_design_sanitize_overhang <- function(x) {
  toupper(gsub("[^ACGT]", "", x))
}

barcode_design_reverse_complement <- function(x) {
  x <- as.character(x)

  vapply(
    x,
    function(seq_value) {
      if (is.na(seq_value) || !nzchar(seq_value)) {
        return(NA_character_)
      }

      bases <- strsplit(toupper(seq_value), "", fixed = TRUE)[[1]]
      paste0(rev(chartr("ACGT", "TGCA", bases)), collapse = "")
    },
    FUN.VALUE = character(1)
  )
}

barcode_design_read_public_ligation_matrix <- function(matrix_path = barcode_design_get_public_matrix_path()) {
  barcode_design_require_packages(c("dplyr", "readr"))

  if (!file.exists(matrix_path)) {
    stop(
      sprintf(
        "Public ligation matrix not found at %s. Place 06_matrix.csv under data/public_reference/ or set BARCODE_LIGATION_MATRIX_PATH.",
        matrix_path
      ),
      call. = FALSE
    )
  }

  matrix_df <- readr::read_csv(matrix_path, show_col_types = FALSE)
  if (ncol(matrix_df) < 2) {
    stop("The ligation matrix must contain one overhang column and at least one paired-sequence column.", call. = FALSE)
  }

  colnames(matrix_df)[1] <- "Overhang"
  matrix_df$Overhang <- barcode_design_sanitize_overhang(as.character(matrix_df$Overhang))
  matrix_df <- matrix_df[matrix_df$Overhang != "", , drop = FALSE]

  paired_sequence_names <- barcode_design_sanitize_overhang(colnames(matrix_df)[-1])
  if (any(!nzchar(paired_sequence_names))) {
    stop("One or more paired-sequence columns could not be sanitized into 4-bp overhang names.", call. = FALSE)
  }
  if (anyDuplicated(paired_sequence_names) > 0) {
    stop("Sanitized paired-sequence column names are duplicated. Please inspect 06_matrix.csv.", call. = FALSE)
  }

  colnames(matrix_df)[-1] <- paired_sequence_names
  matrix_df[-1] <- lapply(matrix_df[-1], as.numeric)
  matrix_df
}

barcode_design_validate_overhang_reference <- function(
  matrix_df,
  overhang_path = barcode_design_get_overhang_reference_path()
) {
  barcode_design_require_packages(c("readr"))

  if (!file.exists(overhang_path)) {
    warning(
      sprintf(
        "Optional overhang reference file not found at %s. Continuing without cross-checking 02_overhangs.csv.",
        overhang_path
      ),
      call. = FALSE
    )
    return(invisible(NULL))
  }

  reference_df <- readr::read_csv(overhang_path, show_col_types = FALSE)
  reference_overhangs <- unique(unlist(lapply(reference_df, function(column) {
    if (!is.character(column) && !is.factor(column)) {
      return(character(0))
    }

    values <- barcode_design_sanitize_overhang(as.character(column))
    values[nchar(values) == 4]
  }), use.names = FALSE))

  if (length(reference_overhangs) == 0) {
    warning("02_overhangs.csv was found but no 4-bp overhang strings were detected in it.", call. = FALSE)
    return(invisible(NULL))
  }

  matrix_overhangs <- unique(c(matrix_df$Overhang, colnames(matrix_df)[-1]))
  missing_overhangs <- setdiff(reference_overhangs, matrix_overhangs)
  if (length(missing_overhangs) > 0) {
    warning(
      sprintf(
        "Detected %d overhang(s) in 02_overhangs.csv that are absent from the matrix, e.g. %s",
        length(missing_overhangs),
        paste(utils::head(missing_overhangs, 5), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  invisible(reference_overhangs)
}

barcode_design_make_ligation_long <- function(matrix_df) {
  barcode_design_require_packages(c("dplyr", "tidyr"))

  long_df <- tidyr::pivot_longer(
    data = matrix_df,
    cols = -Overhang,
    names_to = "PairedSeq",
    values_to = "Count"
  )

  long_df$ComplementSeq <- barcode_design_reverse_complement(long_df$Overhang)
  long_df$Correct <- long_df$PairedSeq == long_df$ComplementSeq
  long_df$Count <- as.numeric(long_df$Count)
  long_df
}

barcode_design_summarise_fidelity <- function(long_df) {
  barcode_design_require_packages(c("dplyr", "tidyr"))

  summary_df <- long_df |>
    dplyr::group_by(.data$Overhang, .data$Correct) |>
    dplyr::summarise(count_sum = sum(.data$Count, na.rm = TRUE), .groups = "drop") |>
    tidyr::pivot_wider(
      names_from = "Correct",
      values_from = "count_sum",
      values_fill = 0
    )

  if (!"TRUE" %in% colnames(summary_df)) {
    summary_df[["TRUE"]] <- 0
  }
  if (!"FALSE" %in% colnames(summary_df)) {
    summary_df[["FALSE"]] <- 0
  }

  summary_df |>
    dplyr::transmute(
      Overhang = .data$Overhang,
      correct_count = .data$`TRUE`,
      incorrect_count = .data$`FALSE`,
      total_count = .data$correct_count + .data$incorrect_count,
      fidelity_rate = dplyr::if_else(
        .data$total_count > 0,
        .data$correct_count / .data$total_count,
        NA_real_
      )
    ) |>
    dplyr::arrange(dplyr::desc(.data$fidelity_rate), dplyr::desc(.data$total_count), .data$Overhang) |>
    dplyr::mutate(rank = dplyr::row_number())
}

barcode_design_score_candidate_cross_ligation <- function(long_df, summary_df, reference_overhangs) {
  barcode_design_require_packages(c("dplyr", "tidyr"))

  candidate_summary <- summary_df |>
    dplyr::filter(!.data$Overhang %in% reference_overhangs)

  if (nrow(candidate_summary) == 0) {
    return(candidate_summary)
  }

  cross_to_reference <- long_df |>
    dplyr::filter(.data$Overhang %in% candidate_summary$Overhang, .data$PairedSeq %in% reference_overhangs) |>
    dplyr::group_by(.data$Overhang) |>
    dplyr::summarise(
      cross_to_reference_count = sum(.data$Count, na.rm = TRUE),
      max_cross_to_reference_count = max(.data$Count, na.rm = TRUE),
      .groups = "drop"
    )

  cross_from_reference <- long_df |>
    dplyr::filter(.data$Overhang %in% reference_overhangs, .data$PairedSeq %in% candidate_summary$Overhang) |>
    dplyr::group_by(.data$PairedSeq) |>
    dplyr::summarise(
      cross_from_reference_count = sum(.data$Count, na.rm = TRUE),
      max_cross_from_reference_count = max(.data$Count, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::rename(Overhang = PairedSeq)

  candidate_summary |>
    dplyr::left_join(cross_to_reference, by = "Overhang") |>
    dplyr::left_join(cross_from_reference, by = "Overhang") |>
    dplyr::mutate(
      dplyr::across(
        c(
          "cross_to_reference_count",
          "max_cross_to_reference_count",
          "cross_from_reference_count",
          "max_cross_from_reference_count"
        ),
        ~ tidyr::replace_na(.x, 0)
      ),
      cross_total_count = .data$cross_to_reference_count + .data$cross_from_reference_count,
      cross_vs_self_ratio = dplyr::if_else(
        .data$correct_count > 0,
        .data$cross_total_count / .data$correct_count,
        NA_real_
      )
    ) |>
    dplyr::arrange(.data$cross_total_count, dplyr::desc(.data$fidelity_rate), dplyr::desc(.data$total_count), .data$Overhang)
}

barcode_design_compute_heatmap_order <- function(matrix_df, reference_overhangs) {
  numeric_matrix <- as.matrix(matrix_df[, -1, drop = FALSE])
  rownames(numeric_matrix) <- matrix_df$Overhang
  mode(numeric_matrix) <- "numeric"

  row_order <- rownames(numeric_matrix)
  col_order <- colnames(numeric_matrix)

  if (nrow(numeric_matrix) > 1) {
    row_order <- rownames(numeric_matrix)[stats::hclust(stats::dist(numeric_matrix))$order]
  }

  if (ncol(numeric_matrix) > 1) {
    col_order <- colnames(numeric_matrix)[stats::hclust(stats::dist(t(numeric_matrix)))$order]
  }

  list(
    row_order = intersect(row_order, reference_overhangs),
    col_order = intersect(col_order, reference_overhangs)
  )
}

barcode_design_read_panel_source <- function(sheet = barcode_design_get_panel_sheet()) {
  barcode_design_require_packages(c("readr"))

  workbook_path <- barcode_design_get_panel_workbook_path()
  if (nzchar(workbook_path)) {
    if (!file.exists(workbook_path)) {
      stop(
        sprintf(
          paste(
            "Barcode 108 workbook not found at %s.",
            "Check BARCODE_PANEL_WORKBOOK_PATH,",
            "restore data/reference/barcode108_panel.xlsx,",
            "or remove the override to use data/processed/barcode108_panel.tsv."
          ),
          workbook_path
        ),
        call. = FALSE
      )
    }

    barcode_design_require_packages(c("readxl"))
    return(readxl::read_excel(
      workbook_path,
      sheet = sheet,
      .name_repair = "unique_quiet"
    ))
  }

  panel_tsv_path <- barcode_design_path("data", "processed", "barcode108_panel.tsv")
  if (file.exists(panel_tsv_path)) {
    return(readr::read_tsv(panel_tsv_path, show_col_types = FALSE))
  }

  panel_csv_path <- barcode_design_path("data", "processed", "barcode108_panel.csv")
  if (file.exists(panel_csv_path)) {
    return(readr::read_csv(panel_csv_path, show_col_types = FALSE))
  }

  stop(
    paste(
      "No barcode108 panel source found.",
      "Provide data/processed/barcode108_panel.tsv or data/processed/barcode108_panel.csv,",
      "or set BARCODE_PANEL_WORKBOOK_PATH to an Excel workbook."
    ),
    call. = FALSE
  )
}

barcode_design_standardize_panel_table <- function(panel_df) {
  barcode_design_require_packages(c("dplyr"))

  panel_df <- as.data.frame(panel_df, stringsAsFactors = FALSE)

  pos_columns <- paste0("pos", 1:6)
  missing_pos_columns <- setdiff(pos_columns, colnames(panel_df))
  if (length(missing_pos_columns) > 0) {
    stop(
      sprintf("The barcode panel workbook is missing expected column(s): %s", paste(missing_pos_columns, collapse = ", ")),
      call. = FALSE
    )
  }

  if (!"code" %in% colnames(panel_df)) {
    panel_df$code <- apply(panel_df[, pos_columns, drop = FALSE], 1, paste, collapse = "-")
  }

  if (!"num" %in% colnames(panel_df)) {
    panel_df$num <- seq_len(nrow(panel_df))
  }

  extra_columns <- setdiff(colnames(panel_df), c(pos_columns, "code", "num"))
  if ("annotation" %in% colnames(panel_df)) {
    annotation <- panel_df$annotation
  } else if (length(extra_columns) > 0) {
    annotation <- panel_df[[extra_columns[1]]]
  } else {
    annotation <- NA_character_
  }

  panel_df$annotation <- as.character(annotation)
  panel_df[, pos_columns] <- lapply(panel_df[, pos_columns, drop = FALSE], as.integer)
  panel_df$num <- as.integer(panel_df$num)

  panel_df |>
    dplyr::select(dplyr::all_of(pos_columns), "code", "num", "annotation")
}

barcode_design_load_barcode15_reference <- function() {
  barcode_design_require_packages(c("readr"))

  reference_path <- barcode_design_path("data", "reference", "barcode15_reference.tsv")
  if (!file.exists(reference_path)) {
    return(barcode_design_default_barcode15_reference())
  }

  readr::read_tsv(reference_path, show_col_types = FALSE)
}

barcode_design_load_primer_slot_map <- function() {
  barcode_design_require_packages(c("readr"))

  reference_path <- barcode_design_path("data", "reference", "primer_slot_map.tsv")
  if (!file.exists(reference_path)) {
    return(barcode_design_default_primer_slot_map())
  }

  readr::read_tsv(reference_path, show_col_types = FALSE)
}

barcode_design_panel_to_long <- function(panel_df, barcode_reference_df) {
  barcode_design_require_packages(c("dplyr", "tidyr"))

  pos_columns <- paste0("pos", 1:6)
  colnames(barcode_reference_df)[colnames(barcode_reference_df) == "source"] <- "barcode_source"
  colnames(barcode_reference_df)[colnames(barcode_reference_df) == "note"] <- "barcode_note"

  panel_long_df <- panel_df |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(pos_columns),
      names_to = "slot_label",
      values_to = "barcode_id"
    ) |>
    dplyr::mutate(slot = as.integer(gsub("^pos", "", .data$slot_label))) |>
    dplyr::left_join(barcode_reference_df, by = "barcode_id")

  colnames(panel_long_df)[colnames(panel_long_df) == "slot.x"] <- "slot"
  colnames(panel_long_df)[colnames(panel_long_df) == "slot.y"] <- "reference_slot"
  panel_long_df
}

barcode_design_build_primer_design <- function(panel_long_df, primer_slot_map_df) {
  barcode_design_require_packages(c("dplyr"))

  colnames(primer_slot_map_df)[colnames(primer_slot_map_df) == "status"] <- "primer_status"
  colnames(primer_slot_map_df)[colnames(primer_slot_map_df) == "note"] <- "primer_note"

  panel_long_df |>
    dplyr::left_join(primer_slot_map_df, by = "slot") |>
    dplyr::mutate(
      slot_matches_reference = dplyr::if_else(
        is.na(.data$reference_slot),
        NA,
        .data$slot == .data$reference_slot
      ),
      barcode_reverse_complement = barcode_design_reverse_complement(.data$barcode_sequence),
      left_primer_construct = dplyr::if_else(
        !is.na(.data$left_primer_seq),
        paste0(
          barcode_design_bsaI_recognition_sequence,
          barcode_design_left_primer_prefix,
          .data$barcode_reverse_complement,
          .data$left_primer_seq
        ),
        NA_character_
      ),
      right_primer_construct = dplyr::if_else(
        !is.na(.data$right_primer_seq),
        paste0(
          barcode_design_bsaI_recognition_sequence,
          barcode_design_right_primer_prefix,
          .data$barcode_sequence,
          .data$right_primer_seq
        ),
        NA_character_
      )
    )
}
