#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Biostrings)
})

default_config <- list(
  input_dir = "fastq",
  output_dir = file.path("output", "01_demultiplex"),
  min_len = 3000L,
  max_len = 4000L,
  primer_mismatch = 7L,
  barcode_mismatch = 1L,
  pair_gap_slack = 0L,
  rescue_partial = FALSE
)

barcode_definitions <- list(
  barcode1 = list(
    primers = c("CGCTGATTGAAGCAGATACCT", "GAATGTGAAATAGGATTTCCGAT"),
    values = c("1" = "CATT", "2" = "TGTT"),
    region = "22:-24"
  ),
  barcode2 = list(
    primers = c("GATATCGTTATTAATACAACACC", "GGATTAAGATTGCAGTTCCTAG"),
    values = c("3" = "GGTC", "4" = "TCCC"),
    region = "24:-23"
  ),
  barcode3 = list(
    primers = c("GAATTGGACATGCGATTTTACCA", "GATGGTCCCGATAAAACGAT"),
    values = c("5" = "GGCG", "6" = "AGGT"),
    region = "24:-21"
  ),
  barcode4 = list(
    primers = c("TACGATTACGTTGTAGTTAATGA", "CCTTCAGGTAATACGATTTTAAC"),
    values = c("7" = "GTCT", "8" = "TCTG"),
    region = "24:-24"
  ),
  barcode5 = list(
    primers = c("TTAAGCTTTTCAACAAAAGGGTC", "TTCACGACGTTCAGAATGAACGA"),
    values = c("9" = "TAAC", "10" = "TAGT"),
    region = "24:-24"
  ),
  barcode6 = list(
    primers = c("GTACAATTGTTAGAAGGTGCAAA", "GCCAATAGGTGTCCTGTATGCTG"),
    values = c("11" = "TTAG", "12" = "GCGG"),
    region = "24:-24"
  )
)

fixed_barcodes <- c(
  barcode3 = "5",
  barcode4 = "7"
)

sample_barcode_positions <- c("barcode1", "barcode2", "barcode5", "barcode6")

sample_map <- c(
  "1-3-9-11" = "sample1",
  "2-3-9-11" = "sample2",
  "1-4-9-11" = "sample3",
  "2-4-9-11" = "sample4",
  "1-3-10-11" = "sample5",
  "2-3-10-11" = "sample6",
  "1-4-10-11" = "sample7",
  "2-4-10-11" = "sample8",
  "1-3-9-12" = "sample9",
  "2-3-9-12" = "sample10"
)

stop_with_usage <- function(message) {
  stop(
    paste(
      message,
      "",
      "Usage:",
      "Rscript code/01_sup_demultiplex.R [--input-dir=fastq] [--output-dir=output/01_demultiplex] [--min-len=3000] [--max-len=4000] [--primer-mismatch=7] [--barcode-mismatch=1] [--pair-gap-slack=0] [--rescue-partial=false]",
      sep = "\n"
    ),
    call. = FALSE
  )
}

parse_args <- function(args) {
  config <- default_config

  for (arg in args) {
    if (!grepl("^--[^=]+=.+$", arg)) {
      stop_with_usage(paste("Invalid argument:", arg))
    }

    key <- sub("^--([^=]+)=.*$", "\\1", arg)
    value <- sub("^--[^=]+=(.*)$", "\\1", arg)

    if (key == "input-dir") {
      config$input_dir <- value
    } else if (key == "output-dir") {
      config$output_dir <- value
    } else if (key == "min-len") {
      config$min_len <- as.integer(value)
    } else if (key == "max-len") {
      config$max_len <- as.integer(value)
    } else if (key == "primer-mismatch") {
      config$primer_mismatch <- as.integer(value)
    } else if (key == "barcode-mismatch") {
      config$barcode_mismatch <- as.integer(value)
    } else if (key == "pair-gap-slack") {
      config$pair_gap_slack <- as.integer(value)
    } else if (key == "rescue-partial") {
      config$rescue_partial <- tolower(value) %in% c("1", "true", "yes")
    } else {
      stop_with_usage(paste("Unknown argument:", arg))
    }
  }

  if (
    is.na(config$min_len) ||
    is.na(config$max_len) ||
    is.na(config$primer_mismatch) ||
    is.na(config$barcode_mismatch) ||
    is.na(config$pair_gap_slack)
  ) {
    stop_with_usage("Length and distance arguments must be integers.")
  }

  if (config$min_len > config$max_len) {
    stop_with_usage("--min-len must be <= --max-len.")
  }

  config
}

ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

find_fastq_files <- function(input_dir) {
  files <- unique(c(
    Sys.glob(file.path(input_dir, "*.fastq.gz")),
    Sys.glob(file.path(input_dir, "*.fastq"))
  ))

  sort(files)
}

file_stem <- function(path) {
  stem <- basename(path)
  if (grepl("\\.gz$", stem)) {
    stem <- tools::file_path_sans_ext(stem)
  }
  if (grepl("\\.fastq$", stem)) {
    stem <- tools::file_path_sans_ext(stem)
  }
  stem
}

read_fastq <- function(path) {
  is_gzip <- grepl("\\.gz$", path)
  con <- if (is_gzip) gzfile(path, open = "rt") else file(path, open = "rt")
  on.exit(close(con))

  lines <- readLines(con, warn = FALSE)

  if (length(lines) %% 4L != 0L) {
    stop(sprintf("FASTQ file does not contain a multiple of 4 lines: %s", path), call. = FALSE)
  }

  data.frame(
    read_id = sub("^@", "", sub("[[:space:]].*$", "", lines[seq(1L, length(lines), 4L)])),
    sequence = lines[seq(2L, length(lines), 4L)],
    quality = lines[seq(4L, length(lines), 4L)],
    stringsAsFactors = FALSE
  )
}

write_fastq <- function(path, fastq_df) {
  con <- gzfile(path, open = "wt")
  on.exit(close(con))

  lines <- character(nrow(fastq_df) * 4L)
  lines[seq(1L, length(lines), 4L)] <- paste0("@", fastq_df$read_id)
  lines[seq(2L, length(lines), 4L)] <- fastq_df$sequence
  lines[seq(3L, length(lines), 4L)] <- "+"
  lines[seq(4L, length(lines), 4L)] <- fastq_df$quality

  writeLines(lines, con = con, sep = "\n")
}

write_fasta <- function(path, read_ids, sequences) {
  con <- file(path, open = "wt")
  on.exit(close(con))

  if (!length(read_ids)) {
    return(invisible(NULL))
  }

  lines <- character(length(read_ids) * 2L)
  lines[seq(1L, length(lines), 2L)] <- paste0(">", read_ids)
  lines[seq(2L, length(lines), 2L)] <- sequences

  writeLines(lines, con = con, sep = "\n")
}

empty_fastq_df <- function() {
  data.frame(
    read_id = character(),
    sequence = character(),
    quality = character(),
    stringsAsFactors = FALSE
  )
}

empty_read_calls <- function(n = 0L) {
  data.frame(
    read_id = rep(NA_character_, n),
    read_length = rep(NA_integer_, n),
    barcode1 = rep(NA_character_, n),
    barcode2 = rep(NA_character_, n),
    barcode3 = rep(NA_character_, n),
    barcode4 = rep(NA_character_, n),
    barcode5 = rep(NA_character_, n),
    barcode6 = rep(NA_character_, n),
    fixed_barcode_ok = rep(FALSE, n),
    fixed_barcode_complete = rep(FALSE, n),
    fixed_barcode_strict_ok = rep(FALSE, n),
    sample_code = rep(NA_character_, n),
    sample = rep(NA_character_, n),
    rescued_partial = rep(FALSE, n),
    stringsAsFactors = FALSE
  )
}

check_dependencies <- function() {
  if (!nzchar(Sys.which("seqkit"))) {
    stop("seqkit is required but was not found in PATH.", call. = FALSE)
  }
}

bioawk_available <- function() {
  nzchar(Sys.which("bioawk"))
}

write_seqid_file <- function(path, read_ids) {
  con <- file(path, open = "wt")
  on.exit(close(con))
  if (!length(read_ids)) {
    return(invisible(NULL))
  }
  writeLines(read_ids, con = con, sep = "\n")
}

filter_fastq_to_fasta <- function(input_fastq, filtered_fastq_df, output_fasta, min_len, max_len) {
  if (bioawk_available()) {
    bioawk_program <- "length($seq)>min_len && length($seq)<max_len {print \">\"$name; print $seq}"
    command <- sprintf(
      "bioawk -c fastx -v min_len=%d -v max_len=%d '%s' %s > %s",
      min_len,
      max_len,
      bioawk_program,
      shQuote(input_fastq),
      shQuote(output_fasta)
    )
    status <- suppressWarnings(system(command, ignore.stderr = TRUE))
    if (is.null(status) || status == 0L) {
      return(invisible(NULL))
    }
    warning(
      sprintf("bioawk filtering failed for %s; falling back to internal FASTA writer.", basename(input_fastq)),
      call. = FALSE
    )
  }

  write_fasta(output_fasta, filtered_fastq_df$read_id, filtered_fastq_df$sequence)
}

write_barcode_table <- function(path, barcode_df) {
  utils::write.table(
    barcode_df,
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
}

seqkit_grep_fastq <- function(input_fastq, seqid_file, output_fastq) {
  result <- system2(
    "seqkit",
    args = c("grep", "-f", seqid_file, "-o", output_fastq, input_fastq),
    stdout = TRUE,
    stderr = TRUE
  )
  status <- attr(result, "status")
  if (!is.null(status) && status != 0L) {
    stop(
      sprintf("seqkit grep failed for %s:\n%s", basename(output_fastq), paste(result, collapse = "\n")),
      call. = FALSE
    )
  }
}

concatenate_fastq_files <- function(input_paths, output_path) {
  if (!length(input_paths)) {
    write_fastq(output_path, empty_fastq_df())
    return(invisible(NULL))
  }

  if (length(input_paths) == 1L) {
    ok <- file.copy(input_paths, output_path, overwrite = TRUE)
    if (!isTRUE(ok)) {
      stop(
        sprintf("file copy failed while creating %s", basename(output_path)),
        call. = FALSE
      )
    }
    return(invisible(NULL))
  }

  out_con <- file(output_path, open = "wb")
  on.exit(close(out_con), add = TRUE)

  for (input_path in input_paths) {
    in_con <- file(input_path, open = "rb")
    tryCatch({
      repeat {
        chunk <- readBin(in_con, what = "raw", n = 1048576L)
        if (!length(chunk)) {
          break
        }
        writeBin(chunk, out_con)
      }
    }, finally = {
      close(in_con)
    })
  }
}

run_seqkit_amplicon <- function(input_fasta, output_fasta, primers, mismatch, region) {
  args <- c(
    "amplicon",
    "-F", primers[1],
    "-R", primers[2],
    "-m", as.character(mismatch),
    "-r", region,
    input_fasta,
    "-o", output_fasta
  )

  result <- system2("seqkit", args = args, stdout = TRUE, stderr = TRUE)
  status <- attr(result, "status")
  if (!is.null(status) && status != 0L) {
    stop(
      sprintf("seqkit amplicon failed for %s:\n%s", basename(output_fasta), paste(result, collapse = "\n")),
      call. = FALSE
    )
  }
}

read_amplicon_fasta <- function(path, barcode_type) {
  if (!file.exists(path) || !file.info(path)$size) {
    return(data.frame(
      seqID = character(),
      Amplicon = character(),
      BarcodeType = character(),
      stringsAsFactors = FALSE
    ))
  }

  sequences <- readDNAStringSet(path)
  if (!length(sequences)) {
    return(data.frame(
      seqID = character(),
      Amplicon = character(),
      BarcodeType = character(),
      stringsAsFactors = FALSE
    ))
  }

  data.frame(
    seqID = names(sequences),
    Amplicon = as.character(sequences),
    BarcodeType = barcode_type,
    stringsAsFactors = FALSE
  )
}

barcode_call_from_amplicon <- function(seq, values, max_mismatch) {
  if (is.na(seq) || !nzchar(seq) || nchar(seq) < 4L) {
    return(NA_character_)
  }

  starts <- seq_len(nchar(seq) - 3L)
  observed_kmers <- vapply(starts, function(start_pos) {
    substr(seq, start_pos, start_pos + 3L)
  }, character(1))

  distances <- vapply(values, function(barcode_seq) {
    min(adist(observed_kmers, barcode_seq))
  }, numeric(1))

  best <- which(distances == min(distances))
  if (length(best) != 1L || distances[best] > max_mismatch) {
    return(NA_character_)
  }

  names(values)[best]
}

extract_barcode_amplicons <- function(filtered_fasta, tmp_dir, primer_mismatch) {
  barcode_tables <- lapply(names(barcode_definitions), function(barcode_type) {
    definition <- barcode_definitions[[barcode_type]]
    output_fasta <- file.path(tmp_dir, paste0(barcode_type, ".fasta"))

    run_seqkit_amplicon(
      input_fasta = filtered_fasta,
      output_fasta = output_fasta,
      primers = definition$primers,
      mismatch = primer_mismatch,
      region = definition$region
    )

    read_amplicon_fasta(output_fasta, barcode_type = barcode_type)
  })

  do.call(rbind, barcode_tables)
}

annotate_barcode_df <- function(barcode_df, barcode_mismatch) {
  if (!nrow(barcode_df)) {
    barcode_df$BarcodeResult <- character()
    return(barcode_df)
  }

  barcode_df$BarcodeResult <- vapply(seq_len(nrow(barcode_df)), function(i) {
    barcode_call_from_amplicon(
      seq = barcode_df$Amplicon[i],
      values = barcode_definitions[[barcode_df$BarcodeType[i]]]$values,
      max_mismatch = barcode_mismatch
    )
  }, character(1))

  barcode_df
}

collapse_barcode_calls <- function(annotated_barcode_df) {
  if (!nrow(annotated_barcode_df)) {
    return(data.frame(
      seqID = character(),
      BarcodeType = character(),
      BarcodeResult = character(),
      stringsAsFactors = FALSE
    ))
  }

  unique_calls <- unique(annotated_barcode_df[, c("seqID", "BarcodeType", "BarcodeResult"), drop = FALSE])
  call_counts <- aggregate(
    BarcodeResult ~ seqID + BarcodeType,
    data = unique_calls,
    FUN = length
  )
  ambiguous <- call_counts$BarcodeResult > 1L
  unique_calls <- merge(unique_calls, call_counts[!ambiguous, , drop = FALSE], by = c("seqID", "BarcodeType"))
  unique_calls[, c("seqID", "BarcodeType", "BarcodeResult.x"), drop = FALSE]
}

format_votes <- function(call_vector) {
  call_vector <- na.omit(call_vector)
  if (!length(call_vector)) {
    return("")
  }

  vote_table <- sort(table(call_vector), decreasing = TRUE)
  paste(sprintf("%s:%d", names(vote_table), as.integer(vote_table)), collapse = ";")
}

majority_call <- function(call_vector) {
  call_vector <- na.omit(call_vector)
  if (!length(call_vector)) {
    return(NA_character_)
  }

  vote_table <- sort(table(call_vector), decreasing = TRUE)
  if (length(vote_table) > 1L && vote_table[1] == vote_table[2]) {
    return(NA_character_)
  }

  names(vote_table)[1]
}

sample_code_from_positions <- function(position_df) {
  paste(
    position_df[[sample_barcode_positions[1]]],
    position_df[[sample_barcode_positions[2]]],
    position_df[[sample_barcode_positions[3]]],
    position_df[[sample_barcode_positions[4]]],
    sep = "-"
  )
}

sample_code_candidates_from_partial <- function(position_df) {
  sample_codes <- names(sample_map)
  candidate_list <- vector("list", nrow(position_df))

  for (i in seq_len(nrow(position_df))) {
    row <- position_df[i, , drop = FALSE]
    known_positions <- sample_barcode_positions[!is.na(row[1, sample_barcode_positions])]

    if (!length(known_positions)) {
      candidate_list[[i]] <- character()
      next
    }

    candidates <- sample_codes
    for (position_name in known_positions) {
      value <- row[[position_name]]
      part_index <- match(position_name, sample_barcode_positions)
      candidates <- candidates[vapply(strsplit(candidates, "-", fixed = TRUE), function(parts) {
        parts[[part_index]] == value
      }, logical(1))]
    }

    candidate_list[[i]] <- candidates
  }

  candidate_list
}

build_read_calls <- function(filtered_fastq, barcode_df, barcode_mismatch, rescue_partial) {
  if (!nrow(filtered_fastq)) {
    return(empty_read_calls())
  }

  read_calls <- empty_read_calls(nrow(filtered_fastq))
  read_calls$read_id <- filtered_fastq$read_id
  read_calls$read_length <- nchar(filtered_fastq$sequence)

  if (nrow(barcode_df)) {
    barcode_df <- annotate_barcode_df(barcode_df, barcode_mismatch = barcode_mismatch)
    collapsed_calls <- collapse_barcode_calls(barcode_df)

    if (nrow(collapsed_calls)) {
      colnames(collapsed_calls)[3] <- "BarcodeResult"
      split_calls <- split(collapsed_calls, collapsed_calls$BarcodeType)
      for (barcode_type in names(split_calls)) {
        call_table <- split_calls[[barcode_type]][, c("seqID", "BarcodeResult"), drop = FALSE]
        idx <- match(call_table$seqID, read_calls$read_id)
        read_calls[idx, barcode_type] <- call_table$BarcodeResult
      }
    }
  }

  read_calls$fixed_barcode_ok <- (is.na(read_calls$barcode3) | read_calls$barcode3 == fixed_barcodes["barcode3"]) &
    (is.na(read_calls$barcode4) | read_calls$barcode4 == fixed_barcodes["barcode4"])
  read_calls$fixed_barcode_complete <- !is.na(read_calls$barcode3) & !is.na(read_calls$barcode4)
  read_calls$fixed_barcode_strict_ok <- read_calls$barcode3 == fixed_barcodes["barcode3"] &
    read_calls$barcode4 == fixed_barcodes["barcode4"]

  variable_complete <- complete.cases(read_calls[, sample_barcode_positions])
  sample_ready <- variable_complete & read_calls$fixed_barcode_ok

  if (any(sample_ready)) {
    read_calls$sample_code[sample_ready] <- sample_code_from_positions(read_calls[sample_ready, , drop = FALSE])
    read_calls$sample[sample_ready] <- unname(sample_map[read_calls$sample_code[sample_ready]])
  }

  read_calls$rescued_partial <- FALSE

  if (isTRUE(rescue_partial)) {
    partial_ready <- rowSums(!is.na(read_calls[, sample_barcode_positions])) >= 3L &
      !sample_ready &
      read_calls$fixed_barcode_ok

    if (any(partial_ready)) {
      candidates <- sample_code_candidates_from_partial(read_calls[partial_ready, , drop = FALSE])
      unique_candidate <- vapply(candidates, length, integer(1)) == 1L

      if (any(unique_candidate)) {
        rescued_codes <- vapply(candidates[unique_candidate], `[`, character(1), 1L)
        rescued_rows <- which(partial_ready)[unique_candidate]
        read_calls$sample_code[rescued_rows] <- rescued_codes
        read_calls$sample[rescued_rows] <- unname(sample_map[rescued_codes])
        read_calls$rescued_partial[rescued_rows] <- TRUE
      }
    }
  }

  read_calls
}

summarize_file <- function(file_name, total_reads, filtered_fastq, read_calls) {
  majority <- c(
    barcode1 = majority_call(read_calls$barcode1),
    barcode2 = majority_call(read_calls$barcode2),
    barcode3 = majority_call(read_calls$barcode3),
    barcode4 = majority_call(read_calls$barcode4),
    barcode5 = majority_call(read_calls$barcode5),
    barcode6 = majority_call(read_calls$barcode6)
  )

  majority_sample_code <- if (any(is.na(majority[sample_barcode_positions]))) {
    NA_character_
  } else {
    paste(
      majority["barcode1"],
      majority["barcode2"],
      majority["barcode5"],
      majority["barcode6"],
      sep = "-"
    )
  }

  file_sample_from_reads <- majority_call(read_calls$sample)

  data.frame(
    file = file_name,
    total_reads = total_reads,
    filtered_reads = nrow(filtered_fastq),
    variable_complete_reads = sum(complete.cases(read_calls[, sample_barcode_positions])),
    fixed_barcode_ok_reads = sum(read_calls$fixed_barcode_ok, na.rm = TRUE),
    fixed_barcode_strict_reads = sum(read_calls$fixed_barcode_strict_ok, na.rm = TRUE),
    read_level_sample_reads = sum(!is.na(read_calls$sample)),
    rescued_partial_reads = sum(read_calls$rescued_partial, na.rm = TRUE),
    barcode1_majority = majority["barcode1"],
    barcode2_majority = majority["barcode2"],
    barcode3_majority = majority["barcode3"],
    barcode4_majority = majority["barcode4"],
    barcode5_majority = majority["barcode5"],
    barcode6_majority = majority["barcode6"],
    barcode1_votes = format_votes(read_calls$barcode1),
    barcode2_votes = format_votes(read_calls$barcode2),
    barcode3_votes = format_votes(read_calls$barcode3),
    barcode4_votes = format_votes(read_calls$barcode4),
    barcode5_votes = format_votes(read_calls$barcode5),
    barcode6_votes = format_votes(read_calls$barcode6),
    sample_votes = format_votes(read_calls$sample),
    majority_sample_code = majority_sample_code,
    file_sample_from_barcode_majority = unname(sample_map[majority_sample_code]),
    file_sample_from_assigned_reads = file_sample_from_reads,
    stringsAsFactors = FALSE
  )
}

main <- function() {
  check_dependencies()
  config <- parse_args(commandArgs(trailingOnly = TRUE))
  files <- find_fastq_files(config$input_dir)

  if (!length(files)) {
    stop(sprintf("No FASTQ files found in %s", config$input_dir), call. = FALSE)
  }

  ensure_dir(config$output_dir)
  ensure_dir(file.path(config$output_dir, "filtered_fastq"))
  ensure_dir(file.path(config$output_dir, "filtered_fasta"))
  ensure_dir(file.path(config$output_dir, "amplicon_barcode"))
  ensure_dir(file.path(config$output_dir, "sample_id"))
  ensure_dir(file.path(config$output_dir, "read_calls"))
  ensure_dir(file.path(config$output_dir, "sample_fastq_by_file"))
  ensure_dir(file.path(config$output_dir, "merged_by_sample"))

  sample_names <- sort(unique(unname(sample_map)))
  merged_fastq_paths <- setNames(vector("list", length(sample_names)), sample_names)

  summary_rows <- list()

  total_files <- length(files)

  for (file_index in seq_along(files)) {
    path <- files[[file_index]]
    message(sprintf("[progress] step1 %d/%d %s", file_index, total_files, basename(path)))
    message(sprintf("Processing %s", basename(path)))
    fastq <- read_fastq(path)
    total_reads <- nrow(fastq)
    read_length <- nchar(fastq$sequence)
    keep <- read_length >= config$min_len & read_length <= config$max_len
    filtered_fastq <- fastq[keep, , drop = FALSE]
    stem <- file_stem(path)

    filtered_fastq_path <- file.path(config$output_dir, "filtered_fastq", paste0(stem, ".filtered.fastq.gz"))
    write_fastq(filtered_fastq_path, filtered_fastq)

    barcode_df <- data.frame(
      seqID = character(),
      Amplicon = character(),
      BarcodeType = character(),
      stringsAsFactors = FALSE
    )

    if (nrow(filtered_fastq)) {
      filtered_fasta_path <- file.path(config$output_dir, "filtered_fasta", paste0(stem, ".filtered.fasta"))
      filter_fastq_to_fasta(
        input_fastq = path,
        filtered_fastq_df = filtered_fastq,
        output_fasta = filtered_fasta_path,
        min_len = config$min_len,
        max_len = config$max_len
      )

      tmp_amplicon_dir <- file.path(config$output_dir, "amplicon_barcode", stem)
      ensure_dir(tmp_amplicon_dir)

      barcode_df <- extract_barcode_amplicons(
        filtered_fasta = filtered_fasta_path,
        tmp_dir = tmp_amplicon_dir,
        primer_mismatch = config$primer_mismatch
      )

      write_barcode_table(
        file.path(config$output_dir, "amplicon_barcode", paste0(stem, "_barcode.txt")),
        barcode_df
      )
    }

    read_calls <- build_read_calls(
      filtered_fastq = filtered_fastq,
      barcode_df = barcode_df,
      barcode_mismatch = config$barcode_mismatch,
      rescue_partial = config$rescue_partial
    )

    utils::write.table(
      read_calls,
      file = file.path(config$output_dir, "read_calls", paste0(stem, ".tsv")),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE,
      col.names = TRUE
    )

    summary_rows[[length(summary_rows) + 1L]] <- summarize_file(
      file_name = basename(path),
      total_reads = total_reads,
      filtered_fastq = filtered_fastq,
      read_calls = read_calls
    )

    for (sample_name in sample_names) {
      sample_ids <- read_calls$read_id[read_calls$sample == sample_name]
      seqid_path <- file.path(
        config$output_dir,
        "sample_id",
        paste0(stem, "_", sample_name, "_SeqID.txt")
      )
      write_seqid_file(seqid_path, sample_ids)

      if (length(sample_ids)) {
        by_file_path <- file.path(
          config$output_dir,
          "sample_fastq_by_file",
          paste0(stem, ".", sample_name, ".fastq.gz")
        )
        seqkit_grep_fastq(
          input_fastq = path,
          seqid_file = seqid_path,
          output_fastq = by_file_path
        )
        merged_fastq_paths[[sample_name]] <- c(merged_fastq_paths[[sample_name]], by_file_path)
      }
    }
  }

  summary_df <- do.call(rbind, summary_rows)
  utils::write.table(
    summary_df,
    file = file.path(config$output_dir, "file_summary.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )

  for (sample_name in sample_names) {
    merged_path <- file.path(config$output_dir, "merged_by_sample", paste0(sample_name, ".fastq.gz"))
    concatenate_fastq_files(merged_fastq_paths[[sample_name]], merged_path)
  }

  message(sprintf("[progress] step1 %d/%d done", total_files, total_files))
}

main()
