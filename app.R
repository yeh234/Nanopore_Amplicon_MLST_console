suppressPackageStartupMessages({
  library(shiny)
})

repo_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
state_dir <- file.path(repo_root, "output", "_shiny_state")
dir.create(state_dir, recursive = TRUE, showWarnings = FALSE)

`%||%` <- function(x, y) {
  if (is.null(x) || !length(x)) {
    y
  } else {
    x
  }
}

is_blank <- function(x) {
  !nzchar(trimws(x %||% ""))
}

resolve_path <- function(path) {
  path <- trimws(path %||% "")
  if (!nzchar(path)) {
    return("")
  }

  expanded <- path.expand(path)
  if (!grepl("^(/|[A-Za-z]:[/\\\\])", expanded)) {
    expanded <- file.path(repo_root, expanded)
  }
  normalizePath(expanded, winslash = "/", mustWork = FALSE)
}

cli_kv <- function(flag, value) {
  sprintf("--%s=%s", flag, shQuote(as.character(value)))
}

cli_flag <- function(flag, enabled) {
  if (isTRUE(enabled)) {
    paste0("--", flag)
  } else {
    ""
  }
}

compact_parts <- function(...) {
  parts <- unlist(list(...), use.names = FALSE)
  parts[nzchar(parts)]
}

build_step1_command <- function(values) {
  paste(
    compact_parts(
      "Rscript",
      "code/01_sup_demultiplex.R",
      cli_kv("input-dir", values$step1_input_dir),
      cli_kv("output-dir", values$step1_output_dir),
      cli_kv("min-len", values$step1_min_len),
      cli_kv("max-len", values$step1_max_len),
      cli_kv("primer-mismatch", values$step1_primer_mismatch),
      cli_kv("barcode-mismatch", values$step1_barcode_mismatch),
      cli_kv("pair-gap-slack", values$step1_pair_gap_slack),
      cli_kv("rescue-partial", tolower(as.character(isTRUE(values$step1_rescue_partial))))
    ),
    collapse = " "
  )
}

build_step2_command <- function(values) {
  paste(
    compact_parts(
      "bash",
      "code/04_assembly_polish_pipeline.sh",
      cli_kv("sample-fastq-dir", values$step2_sample_fastq_dir),
      cli_kv("assembly-output-dir", values$step2_assembly_output_dir),
      cli_kv("mlst-output-dir", values$step2_mlst_output_dir),
      if (!is_blank(values$step2_assembly_env)) cli_kv("assembly-env", values$step2_assembly_env) else "",
      if (!is_blank(values$step2_mlst_env)) cli_kv("mlst-env", values$step2_mlst_env) else "",
      cli_kv("threads", values$step2_threads),
      cli_kv("mlst-threads", values$step2_mlst_threads),
      cli_kv("genome-size", values$step2_genome_size),
      cli_kv("asm-coverage", values$step2_asm_coverage),
      cli_kv("min-overlap", values$step2_min_overlap),
      if (!is_blank(values$step2_samples)) cli_kv("samples", values$step2_samples) else "",
      cli_flag("skip-polish", values$step2_skip_polish),
      cli_flag("skip-assembly", values$step2_skip_assembly),
      cli_flag("skip-mlst", values$step2_skip_mlst),
      cli_flag("medaka-bacteria", values$step2_medaka_bacteria),
      if (!is_blank(values$step2_medaka_model)) cli_kv("medaka-model", values$step2_medaka_model) else "",
      cli_kv("scheme", values$step2_scheme)
    ),
    collapse = " "
  )
}

build_step3_command <- function(values) {
  paste(
    compact_parts(
      "python3",
      "code/05_infer_st_from_partial_mlst.py",
      cli_kv("mlst-summary", values$step3_mlst_summary),
      cli_kv("profile-table", values$step3_profile_table),
      cli_kv("output-dir", values$step3_output_dir),
      cli_kv("top-n", values$step3_top_n)
    ),
    collapse = " "
  )
}

build_step4_command <- function(values) {
  paste(
    compact_parts(
      "python3",
      "code/06_build_publication_style_mlst_table.py",
      cli_kv("assembly-summary", values$step4_assembly_summary),
      cli_kv("mlst-summary", values$step4_mlst_summary),
      cli_kv("st-guess-summary", values$step4_st_guess_summary),
      cli_kv("output-path", values$step4_output_path)
    ),
    collapse = " "
  )
}

build_full_command <- function(values, pretty = FALSE) {
  separator <- if (isTRUE(pretty)) " && \\\n" else " && "
  paste(
    c(
      build_step1_command(values),
      build_step2_command(values),
      build_step3_command(values),
      build_step4_command(values)
    ),
    collapse = separator
  )
}

demo_values_from <- function(values) {
  demo_root <- "output/demo_run"
  step1_output_dir <- file.path(demo_root, "01_demultiplex")
  step2_assembly_output_dir <- file.path(demo_root, "02_assembly")
  step2_mlst_output_dir <- file.path(demo_root, "03_mlst")
  step3_output_dir <- file.path(demo_root, "04_st_guess")

  modifyList(
    values,
    list(
      step1_input_dir = "demo",
      step1_output_dir = step1_output_dir,
      step2_sample_fastq_dir = file.path(step1_output_dir, "merged_by_sample"),
      step2_assembly_output_dir = step2_assembly_output_dir,
      step2_mlst_output_dir = step2_mlst_output_dir,
      step3_mlst_summary = file.path(step2_mlst_output_dir, "mlst_sample_summary.tsv"),
      step3_output_dir = step3_output_dir,
      step4_assembly_summary = file.path(step2_assembly_output_dir, "assembly_polish_summary.tsv"),
      step4_mlst_summary = file.path(step2_mlst_output_dir, "mlst_sample_summary.tsv"),
      step4_st_guess_summary = file.path(step3_output_dir, "st_guess_summary.tsv"),
      step4_output_path = file.path(demo_root, "Result_Table", "SUP_demo_publication_style_mlst.tsv")
    )
  )
}

make_job_paths <- function(label) {
  safe_label <- gsub("[^A-Za-z0-9]+", "_", tolower(label))
  stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  base <- file.path(state_dir, paste0(stamp, "_", safe_label))
  list(
    log = paste0(base, ".log"),
    status = paste0(base, ".status"),
    script = paste0(base, ".sh")
  )
}

launch_background_job <- function(shell_command, label) {
  paths <- make_job_paths(label)
  script_lines <- c(
    "#!/bin/bash",
    "set -u",
    "set -o pipefail",
    sprintf("status_file=%s", shQuote(paths$status)),
    "trap 'status=$?; printf \"%s\\n\" \"$status\" > \"$status_file\"' EXIT",
    sprintf("cd %s", shQuote(repo_root)),
    shell_command
  )
  writeLines(script_lines, con = paths$script, sep = "\n")
  Sys.chmod(paths$script, mode = "0755")

  launcher <- if (nzchar(Sys.which("nohup"))) "nohup" else "bash"
  launch_args <- if (identical(launcher, "nohup")) c("bash", paths$script) else paths$script
  system2(launcher, launch_args, wait = FALSE, stdout = paths$log, stderr = paths$log)

  list(
    pid = NA_integer_,
    log = paths$log,
    status = paths$status,
    script = paths$script
  )
}

read_status_code <- function(path) {
  if (!file.exists(path)) {
    return(NA_integer_)
  }

  value <- trimws(readLines(path, warn = FALSE, n = 1L))
  suppressWarnings(as.integer(value))
}

read_log_tail <- function(path, n = 200L) {
  if (is_blank(path) || !file.exists(path)) {
    return("No log yet.")
  }

  lines <- readLines(path, warn = FALSE)
  if (!length(lines)) {
    return("Log file is empty.")
  }

  paste(tail(lines, n), collapse = "\n")
}

read_preview_tsv <- function(path, n_max = 8L) {
  if (is_blank(path) || !file.exists(path)) {
    return(data.frame(status = sprintf("Not found: %s", path), check.names = FALSE))
  }

  info <- file.info(path)
  if (isTRUE(info$isdir)) {
    return(data.frame(status = sprintf("Directory exists: %s", path), check.names = FALSE))
  }

  tryCatch(
    utils::read.delim(
      path,
      sep = "\t",
      nrows = n_max,
      check.names = FALSE,
      stringsAsFactors = FALSE
    ),
    error = function(error) {
      data.frame(error = conditionMessage(error), check.names = FALSE)
    }
  )
}

read_full_tsv <- function(path) {
  if (is_blank(path) || !file.exists(path)) {
    return(data.frame(status = sprintf("Not found: %s", path), check.names = FALSE))
  }

  info <- file.info(path)
  if (isTRUE(info$isdir)) {
    return(data.frame(status = sprintf("Directory exists: %s", path), check.names = FALSE))
  }

  tryCatch(
    utils::read.delim(
      path,
      sep = "\t",
      check.names = FALSE,
      stringsAsFactors = FALSE
    ),
    error = function(error) {
      data.frame(error = conditionMessage(error), check.names = FALSE)
    }
  )
}

list_fastq_files <- function(input_dir) {
  dir_path <- resolve_path(input_dir)
  if (!dir.exists(dir_path)) {
    return(character())
  }

  sort(unique(c(
    Sys.glob(file.path(dir_path, "*.fastq.gz")),
    Sys.glob(file.path(dir_path, "*.fastq"))
  )))
}

normalize_sample_name <- function(path) {
  sample_name <- basename(path)
  sample_name <- sub("\\.gz$", "", sample_name)
  sub("\\.fastq$", "", sample_name)
}

parse_csv_items <- function(x) {
  items <- trimws(strsplit(x %||% "", ",", fixed = TRUE)[[1]])
  items[nzchar(items)]
}

list_step2_samples <- function(sample_fastq_dir, samples_csv = "") {
  dir_path <- resolve_path(sample_fastq_dir)
  if (!dir.exists(dir_path)) {
    return(character())
  }

  files <- sort(unique(c(
    Sys.glob(file.path(dir_path, "sample*.fastq.gz")),
    Sys.glob(file.path(dir_path, "sample*.fastq"))
  )))
  samples <- unique(vapply(files, normalize_sample_name, character(1)))
  wanted <- parse_csv_items(samples_csv)
  if (length(wanted)) {
    samples <- samples[samples %in% wanted]
  }
  sort(samples)
}

read_log_lines <- function(path) {
  if (is_blank(path) || !file.exists(path)) {
    return(character())
  }
  readLines(path, warn = FALSE)
}

path_ready_after <- function(path, started_at = NULL) {
  if (is_blank(path) || !(file.exists(path) || dir.exists(path))) {
    return(FALSE)
  }
  if (is.null(started_at)) {
    return(TRUE)
  }

  info <- file.info(path)
  if (!nrow(info) || is.na(info$mtime[1])) {
    return(FALSE)
  }

  as.numeric(difftime(info$mtime[1], started_at, units = "secs")) >= -1
}

extract_progress_counter <- function(lines, step_token) {
  pattern <- sprintf("^\\[progress\\] %s ([0-9]+)/([0-9]+)", step_token)
  matched <- regexec(pattern, lines)
  captures <- regmatches(lines, matched)
  captures <- captures[lengths(captures) > 0]
  if (!length(captures)) {
    return(NULL)
  }

  last <- captures[[length(captures)]]
  list(
    current = as.integer(last[2]),
    total = as.integer(last[3])
  )
}

count_step2_processed_samples <- function(lines) {
  matched <- regexec("^Processing (sample[^[:space:]]+)$", lines)
  captures <- regmatches(lines, matched)
  captures <- captures[lengths(captures) > 0]
  if (!length(captures)) {
    return(0L)
  }

  samples <- unique(vapply(captures, `[`, character(1), 2))
  length(samples)
}

progress_state <- function(status, percent, summary, detail = "", indeterminate = FALSE) {
  list(
    status = status,
    percent = max(0, min(100, round(percent))),
    summary = summary,
    detail = detail,
    indeterminate = indeterminate
  )
}

detect_failed_step <- function(values, rv) {
  if (isTRUE(rv$active) || is.na(rv$last_exit) || identical(rv$last_exit, 0L)) {
    return("")
  }

  if (identical(rv$label, "Step 1 demultiplex")) {
    return("step1")
  }
  if (identical(rv$label, "Step 2 assembly_mlst")) {
    return("step2")
  }
  if (identical(rv$label, "Step 3 st_guess")) {
    return("step3")
  }
  if (identical(rv$label, "Step 4 result_table")) {
    return("step4")
  }
  if (!identical(rv$label, "Full workflow")) {
    return("")
  }

  step1_summary <- resolve_path(file.path(values$step1_output_dir, "file_summary.tsv"))
  step2_summary <- resolve_path(file.path(values$step2_mlst_output_dir, "mlst_sample_summary.tsv"))
  step3_summary <- resolve_path(file.path(values$step3_output_dir, "st_guess_summary.tsv"))
  step4_output <- resolve_path(values$step4_output_path)

  if (!path_ready_after(step1_summary, rv$started_at)) {
    return("step1")
  }
  if (!path_ready_after(step2_summary, rv$started_at)) {
    return("step2")
  }
  if (!path_ready_after(step3_summary, rv$started_at)) {
    return("step3")
  }
  if (!path_ready_after(step4_output, rv$started_at)) {
    return("step4")
  }

  ""
}

build_progress_states <- function(values, rv) {
  lines <- read_log_lines(rv$log_path)
  failed_step <- detect_failed_step(values, rv)
  current_label <- rv$label %||% ""
  active <- isTRUE(rv$active)
  started_at <- rv$started_at
  prefer_current_run <- !active && identical(current_label, "Full workflow") && !is.na(rv$last_exit)

  step1_summary <- resolve_path(file.path(values$step1_output_dir, "file_summary.tsv"))
  step2_assembly_summary <- resolve_path(file.path(values$step2_assembly_output_dir, "assembly_polish_summary.tsv"))
  step2_mlst_summary <- resolve_path(file.path(values$step2_mlst_output_dir, "mlst_sample_summary.tsv"))
  step3_summary <- resolve_path(file.path(values$step3_output_dir, "st_guess_summary.tsv"))
  step4_output <- resolve_path(values$step4_output_path)

  step1_done_current <- path_ready_after(step1_summary, started_at)
  step2_done_current <- path_ready_after(
    if (isTRUE(values$step2_skip_mlst)) step2_assembly_summary else step2_mlst_summary,
    started_at
  )
  step3_done_current <- path_ready_after(step3_summary, started_at)
  step4_done_current <- path_ready_after(step4_output, started_at)

  step1_total <- length(list_fastq_files(values$step1_input_dir))
  step2_total <- length(list_step2_samples(values$step2_sample_fastq_dir, values$step2_samples))

  step1_counter <- extract_progress_counter(lines, "step1")
  step3_counter <- extract_progress_counter(lines, "step3")
  step4_counter <- extract_progress_counter(lines, "step4")

  step1_processed <- if (!is.null(step1_counter)) {
    step1_counter$current
  } else {
    sum(grepl("^Processing .*\\.fastq(\\.gz)?$", lines))
  }
  step2_processed <- count_step2_processed_samples(lines)
  step2_stage_assembly <- any(grepl("^\\[progress\\] step2 stage assembly$", lines)) ||
    any(grepl("^Running Step 4-5 assembly \\+ polish$", lines))
  step2_stage_mlst <- any(grepl("^\\[progress\\] step2 stage mlst$", lines)) ||
    any(grepl("^Running Step 6 standard mlst caller$", lines))

  states <- list()

  step1_active <- active && current_label %in% c("Step 1 demultiplex", "Full workflow")
  if (step1_active) {
    if (step1_done_current) {
      states$step1 <- progress_state("complete", 100, "Step 1 complete", sprintf("%d/%d FASTQ files", step1_total, step1_total))
    } else {
      percent <- if (step1_total > 0) min(97, (step1_processed / step1_total) * 100) else 10
      detail <- if (step1_total > 0) sprintf("%d/%d FASTQ files", min(step1_processed, step1_total), step1_total) else "Scanning FASTQ inputs"
      states$step1 <- progress_state("running", percent, "Running demultiplex", detail, indeterminate = step1_total == 0)
    }
  } else if (identical(failed_step, "step1")) {
    states$step1 <- progress_state("failed", if (step1_total > 0) min(95, (step1_processed / step1_total) * 100) else 0, "Step 1 failed", "See log for the failing FASTQ")
  } else if (if (prefer_current_run) path_ready_after(step1_summary, started_at) else file.exists(step1_summary)) {
    states$step1 <- progress_state("complete", 100, "Step 1 output ready", basename(step1_summary))
  } else {
    states$step1 <- progress_state("pending", 0, "Step 1 not started", "Demultiplex")
  }

  step2_active <- active && current_label %in% c("Step 2 assembly_mlst", "Full workflow")
  if (step2_active && identical(current_label, "Full workflow") && !step1_done_current) {
    states$step2 <- progress_state("pending", 0, "Waiting for Step 1", "Assembly + MLST")
  } else if (step2_active) {
    if (step2_done_current) {
      states$step2 <- progress_state("complete", 100, "Step 2 complete", if (isTRUE(values$step2_skip_mlst)) "Assembly summary ready" else "MLST summary ready")
    } else if (isTRUE(values$step2_skip_assembly) && !isTRUE(values$step2_skip_mlst)) {
      states$step2 <- progress_state("running", 88, "Running MLST", "Assembly skipped for this run", indeterminate = TRUE)
    } else if (isTRUE(values$step2_skip_assembly) && isTRUE(values$step2_skip_mlst)) {
      states$step2 <- progress_state("complete", 100, "Step 2 complete", "Assembly and MLST both skipped")
    } else if (!step2_stage_assembly && !step2_stage_mlst) {
      states$step2 <- progress_state("running", 5, "Initializing Step 2", "Preparing assembly inputs", indeterminate = TRUE)
    } else if (!step2_stage_mlst) {
      percent <- if (step2_total > 0) 5 + min(78, (step2_processed / step2_total) * 78) else 40
      detail <- if (step2_total > 0) sprintf("Assembly %d/%d samples", min(step2_processed, step2_total), step2_total) else "Running assembly"
      states$step2 <- progress_state("running", percent, "Assembly + polish", detail)
    } else {
      states$step2 <- progress_state("running", 90, "Running MLST", "Assembly finished; waiting for mlst output", indeterminate = TRUE)
    }
  } else if (identical(failed_step, "step2")) {
    states$step2 <- progress_state("failed", 0, "Step 2 failed", "Assembly or MLST stage did not finish")
  } else if (if (prefer_current_run) {
    path_ready_after(if (isTRUE(values$step2_skip_mlst)) step2_assembly_summary else step2_mlst_summary, started_at)
  } else {
    file.exists(if (isTRUE(values$step2_skip_mlst)) step2_assembly_summary else step2_mlst_summary)
  }) {
    states$step2 <- progress_state("complete", 100, "Step 2 output ready", if (isTRUE(values$step2_skip_mlst)) basename(step2_assembly_summary) else basename(step2_mlst_summary))
  } else {
    states$step2 <- progress_state("pending", 0, "Step 2 not started", "Assembly + MLST")
  }

  step3_active <- active && current_label %in% c("Step 3 st_guess", "Full workflow")
  if (step3_active && identical(current_label, "Full workflow") && !step2_done_current) {
    states$step3 <- progress_state("pending", 0, "Waiting for Step 2", "ST guess")
  } else if (step3_active) {
    if (step3_done_current) {
      states$step3 <- progress_state("complete", 100, "Step 3 complete", basename(step3_summary))
    } else if (!is.null(step3_counter)) {
      percent <- if (!is.na(step3_counter$total) && step3_counter$total > 0) min(98, (step3_counter$current / step3_counter$total) * 100) else 20
      detail <- if (!is.na(step3_counter$total) && step3_counter$total > 0) {
        sprintf("%d/%d samples", step3_counter$current, step3_counter$total)
      } else {
        "Inferring ST candidates"
      }
      states$step3 <- progress_state("running", percent, "Running ST guess", detail)
    } else {
      states$step3 <- progress_state("running", 15, "Initializing Step 3", "Loading MLST summary", indeterminate = TRUE)
    }
  } else if (identical(failed_step, "step3")) {
    states$step3 <- progress_state("failed", 0, "Step 3 failed", "Could not write st_guess_summary.tsv")
  } else if (if (prefer_current_run) path_ready_after(step3_summary, started_at) else file.exists(step3_summary)) {
    states$step3 <- progress_state("complete", 100, "Step 3 output ready", basename(step3_summary))
  } else {
    states$step3 <- progress_state("pending", 0, "Step 3 not started", "ST guess")
  }

  step4_active <- active && current_label %in% c("Step 4 result_table", "Full workflow")
  if (step4_active && identical(current_label, "Full workflow") && !step3_done_current) {
    states$step4 <- progress_state("pending", 0, "Waiting for Step 3", "Result Table")
  } else if (step4_active) {
    if (step4_done_current) {
      states$step4 <- progress_state("complete", 100, "Step 4 complete", basename(step4_output))
    } else if (!is.null(step4_counter)) {
      percent <- if (!is.na(step4_counter$total) && step4_counter$total > 0) min(98, (step4_counter$current / step4_counter$total) * 100) else 20
      detail <- if (!is.na(step4_counter$total) && step4_counter$total > 0) {
        sprintf("%d/%d samples", step4_counter$current, step4_counter$total)
      } else {
        "Building final TSV"
      }
      states$step4 <- progress_state("running", percent, "Building result table", detail)
    } else {
      states$step4 <- progress_state("running", 15, "Initializing Step 4", "Collecting summary tables", indeterminate = TRUE)
    }
  } else if (identical(failed_step, "step4")) {
    states$step4 <- progress_state("failed", 0, "Step 4 failed", "Could not write final result table")
  } else if (if (prefer_current_run) path_ready_after(step4_output, started_at) else file.exists(step4_output)) {
    states$step4 <- progress_state("complete", 100, "Step 4 output ready", basename(step4_output))
  } else {
    states$step4 <- progress_state("pending", 0, "Step 4 not started", "Result Table")
  }

  states
}

run_command_capture <- function(command, args = character()) {
  tryCatch(
    {
      output <- suppressWarnings(system2(command, args, stdout = TRUE, stderr = TRUE))
      status <- attr(output, "status")
      if (is.null(status)) {
        status <- 0L
      }
      list(
        status = as.integer(status),
        output = trimws(output)
      )
    },
    error = function(error) {
      list(
        status = 1L,
        output = trimws(conditionMessage(error))
      )
    }
  )
}

list_conda_envs <- function(conda_bin = Sys.which("conda")) {
  empty <- data.frame(
    name = character(),
    path = character(),
    stringsAsFactors = FALSE
  )

  if (!nzchar(conda_bin)) {
    return(empty)
  }

  result <- run_command_capture(conda_bin, c("env", "list"))
  lines <- result$output
  lines <- lines[nzchar(lines)]
  lines <- lines[!grepl("^(#|sysctl:)", lines)]

  matched <- regexec("^([^[:space:]]+)\\s+\\*?\\s+(.+)$", lines)
  captures <- regmatches(lines, matched)
  captures <- captures[lengths(captures) > 0]
  if (!length(captures)) {
    return(empty)
  }

  data.frame(
    name = vapply(captures, `[`, character(1), 2),
    path = vapply(captures, `[`, character(1), 3),
    stringsAsFactors = FALSE
  )
}

conda_env_exists <- function(envs, env_name) {
  nzchar(env_name) && any(envs$name == env_name)
}

first_nonempty <- function(x, fallback = "") {
  x <- trimws(x)
  x <- x[nzchar(x)]
  if (!length(x)) {
    fallback
  } else {
    x[[1]]
  }
}

conda_tool_paths <- function(conda_bin, env_name, tools) {
  paths <- setNames(rep("", length(tools)), tools)

  if (!nzchar(conda_bin) || !nzchar(env_name) || !length(tools)) {
    return(list(paths = paths, error = ""))
  }

  shell_command <- paste0(
    "for tool in ", paste(shQuote(tools), collapse = " "), "; do ",
    "path=$(command -v \"$tool\" 2>/dev/null || true); ",
    "printf '%s\\t%s\\n' \"$tool\" \"$path\"; ",
    "done"
  )
  result <- run_command_capture(conda_bin, c("run", "-n", env_name, "sh", "-lc", shQuote(shell_command)))

  tab_lines <- result$output[grepl("\t", result$output, fixed = TRUE)]
  for (line in tab_lines) {
    parts <- strsplit(line, "\t", fixed = TRUE)[[1]]
    tool_name <- trimws(parts[[1]])
    tool_path <- if (length(parts) >= 2) trimws(parts[[2]]) else ""
    if (tool_name %in% names(paths)) {
      paths[[tool_name]] <- tool_path
    }
  }

  error_detail <- ""
  if (!identical(result$status, 0L) && !any(nzchar(paths))) {
    error_detail <- first_nonempty(result$output, sprintf("conda run failed for env: %s", env_name))
  }

  list(paths = paths, error = error_detail)
}

diagnostic_rows <- function(step, context, tools, paths, missing_detail) {
  data.frame(
    Step = step,
    Context = context,
    Tool = tools,
    Available = ifelse(nzchar(paths), "yes", "no"),
    Detail = ifelse(nzchar(paths), unname(paths), missing_detail),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

tool_diagnostics <- function(assembly_env = "", mlst_env = "") {
  shell_tools_general <- c("Rscript", "python3", "bash", "conda")
  shell_tools_step1 <- c("seqkit", "bioawk")
  assembly_tools <- c("flye", "bwa", "samtools", "medaka_consensus", "medaka")
  mlst_tools <- c("mlst")

  shell_general_paths <- Sys.which(shell_tools_general)
  shell_step1_paths <- Sys.which(shell_tools_step1)
  conda_bin <- unname(Sys.which("conda")[["conda"]])
  conda_envs <- list_conda_envs(conda_bin)

  shell_rows <- rbind(
    diagnostic_rows(
      "General",
      "current shell",
      shell_tools_general,
      shell_general_paths,
      "Not found in current shell PATH"
    ),
    diagnostic_rows(
      "Step 1",
      "current shell",
      shell_tools_step1,
      shell_step1_paths,
      "Not found in current shell PATH"
    )
  )

  build_env_rows <- function(step, env_name, tools) {
    context <- if (nzchar(env_name)) {
      sprintf("conda env: %s", env_name)
    } else {
      "conda env not set"
    }

    if (!nzchar(conda_bin)) {
      return(diagnostic_rows(step, context, tools, setNames(rep("", length(tools)), tools), "conda not found in current shell PATH"))
    }
    if (!nzchar(env_name)) {
      return(diagnostic_rows(step, context, tools, setNames(rep("", length(tools)), tools), "Environment name is blank"))
    }
    if (!conda_env_exists(conda_envs, env_name)) {
      return(diagnostic_rows(step, context, tools, setNames(rep("", length(tools)), tools), sprintf("Conda env not found: %s", env_name)))
    }

    env_probe <- conda_tool_paths(conda_bin, env_name, tools)
    fallback_detail <- if (nzchar(env_probe$error)) {
      env_probe$error
    } else {
      sprintf("Not found inside conda env: %s", env_name)
    }

    diagnostic_rows(step, context, tools, env_probe$paths, fallback_detail)
  }

  rbind(
    shell_rows,
    build_env_rows("Step 2 assembly", assembly_env %||% "", assembly_tools),
    build_env_rows("Step 2 MLST", mlst_env %||% "", mlst_tools)
  )
}

card_panel <- function(title, ..., subtitle = NULL, class = "") {
  div(
    class = paste("card-panel", class),
    div(
      class = "card-head",
      h3(title),
      if (!is.null(subtitle)) p(class = "card-subtitle", subtitle)
    ),
    ...
  )
}

status_chip <- function(label, path, ready) {
  div(
    class = paste("status-chip", if (ready) "ready" else "missing"),
    div(class = "status-top", span(class = "status-label", label), span(class = "status-state", if (ready) "ready" else "missing")),
    div(class = "status-path", path)
  )
}

table_preview_output <- function(output_id) {
  div(class = "table-wrap", tableOutput(output_id))
}

progress_panel_ui <- function(title, state) {
  badge_label <- switch(
    state$status,
    pending = "Pending",
    running = "Running",
    complete = "Complete",
    failed = "Failed",
    state$status
  )
  fill_class <- paste(
    "progress-fill",
    state$status,
    if (isTRUE(state$indeterminate) && identical(state$status, "running")) "indeterminate" else ""
  )
  fill_width <- if (isTRUE(state$indeterminate) && identical(state$status, "running")) 100 else state$percent

  div(
    class = "progress-card",
    div(
      class = "progress-head",
      span(class = "progress-title", title),
      span(class = paste("progress-badge", state$status), badge_label)
    ),
    div(class = "progress-summary", state$summary),
    div(
      class = "progress-track",
      div(class = fill_class, style = sprintf("width: %s%%;", fill_width))
    ),
    div(
      class = "progress-foot",
      span(state$detail),
      span(sprintf("%d%%", state$percent))
    )
  )
}

ui <- fluidPage(
  tags$head(
    tags$title("Nanopore Amplicon MLST Console"),
    tags$style(
      HTML(
        "
        body {
          background:
            radial-gradient(circle at top left, rgba(185, 138, 74, 0.25), transparent 32%),
            radial-gradient(circle at right top, rgba(77, 115, 92, 0.18), transparent 28%),
            linear-gradient(180deg, #f6efe4 0%, #efe5d4 100%);
          color: #2c271f;
          font-family: 'Avenir Next', 'Segoe UI', sans-serif;
        }
        .container-fluid {
          max-width: 1480px;
          padding: 28px 28px 40px;
        }
        .hero {
          background: rgba(255, 251, 244, 0.78);
          border: 1px solid rgba(67, 53, 31, 0.12);
          border-radius: 24px;
          box-shadow: 0 18px 40px rgba(74, 56, 31, 0.08);
          padding: 28px 30px 22px;
          margin-bottom: 20px;
          backdrop-filter: blur(6px);
        }
        .hero-kicker {
          text-transform: uppercase;
          letter-spacing: 0.18em;
          font-size: 12px;
          color: #87663a;
          margin-bottom: 10px;
        }
        .hero h1 {
          font-family: Georgia, 'Times New Roman', serif;
          font-size: 42px;
          margin: 0 0 10px;
          color: #241d17;
        }
        .hero p {
          max-width: 920px;
          font-size: 15px;
          line-height: 1.6;
          margin-bottom: 10px;
        }
        .hero-meta {
          display: flex;
          flex-wrap: wrap;
          gap: 10px;
          margin-top: 12px;
        }
        .hero-meta span {
          display: inline-flex;
          align-items: center;
          gap: 8px;
          padding: 8px 12px;
          border-radius: 999px;
          background: #eee1cb;
          color: #4a3b28;
          font-size: 12px;
        }
        .nav-tabs {
          border-bottom: none;
          margin-bottom: 16px;
        }
        .nav-tabs > li > a {
          border: none;
          border-radius: 999px;
          background: rgba(255, 250, 241, 0.85);
          color: #5d4a31;
          margin-right: 8px;
          padding: 10px 16px;
        }
        .nav-tabs > li.active > a,
        .nav-tabs > li.active > a:focus,
        .nav-tabs > li.active > a:hover {
          background: #3d5b46;
          color: #fdfaf4;
          border: none;
        }
        .card-panel {
          background: rgba(255, 252, 247, 0.9);
          border: 1px solid rgba(67, 53, 31, 0.12);
          border-radius: 22px;
          box-shadow: 0 14px 28px rgba(74, 56, 31, 0.07);
          padding: 22px;
          margin-bottom: 18px;
        }
        .card-head h3 {
          font-family: Georgia, 'Times New Roman', serif;
          font-size: 26px;
          margin-top: 0;
          margin-bottom: 6px;
          color: #241d17;
        }
        .card-subtitle {
          color: #685743;
          margin-bottom: 16px;
        }
        .help-block, .form-group label, .checkbox label {
          color: #4d4337;
        }
        .btn {
          border-radius: 999px;
          border: none;
          padding: 10px 16px;
          font-weight: 600;
          box-shadow: none;
        }
        .btn-primary {
          background: #8a5a31;
          color: #fffaf5;
        }
        .btn-default {
          background: #e8dac0;
          color: #4f412f;
        }
        .btn-primary:hover,
        .btn-primary:focus {
          background: #734726;
        }
        .btn-default:hover,
        .btn-default:focus {
          background: #dbc8a7;
        }
        .section-note {
          font-size: 13px;
          color: #675846;
          margin-bottom: 12px;
        }
        .status-grid {
          display: grid;
          gap: 12px;
        }
        .status-chip {
          padding: 14px 15px;
          border-radius: 18px;
          border: 1px solid rgba(67, 53, 31, 0.12);
          background: #fbf5ec;
        }
        .status-chip.ready {
          border-color: rgba(58, 111, 83, 0.28);
          background: rgba(226, 243, 233, 0.8);
        }
        .status-chip.missing {
          border-color: rgba(138, 90, 49, 0.22);
          background: rgba(251, 243, 232, 0.92);
        }
        .status-top {
          display: flex;
          justify-content: space-between;
          gap: 12px;
          margin-bottom: 6px;
        }
        .status-label {
          font-weight: 700;
        }
        .status-state {
          text-transform: uppercase;
          letter-spacing: 0.08em;
          font-size: 11px;
          color: #6c5b46;
        }
        .status-path {
          font-family: Menlo, Consolas, monospace;
          font-size: 12px;
          color: #544635;
          word-break: break-all;
        }
        .job-banner {
          padding: 18px;
          border-radius: 18px;
          border: 1px solid rgba(67, 53, 31, 0.12);
        }
        .job-banner h4 {
          margin-top: 0;
          margin-bottom: 8px;
          font-size: 22px;
          font-family: Georgia, 'Times New Roman', serif;
        }
        .job-banner.live {
          background: rgba(247, 232, 206, 0.92);
          border-color: rgba(138, 90, 49, 0.22);
        }
        .job-banner.success {
          background: rgba(223, 242, 230, 0.92);
          border-color: rgba(58, 111, 83, 0.24);
        }
        .job-banner.failed {
          background: rgba(251, 229, 221, 0.94);
          border-color: rgba(148, 69, 47, 0.24);
        }
        .job-banner.idle {
          background: rgba(244, 236, 223, 0.92);
        }
        .log-panel pre,
        .command-box pre {
          min-height: 240px;
          max-height: 520px;
          overflow: auto;
          background: #231d17;
          color: #f7f0e7;
          border: none;
          border-radius: 18px;
          padding: 18px;
          font-family: Menlo, Consolas, monospace;
          font-size: 12px;
          line-height: 1.55;
          white-space: pre-wrap;
          word-break: break-word;
        }
        .command-box pre {
          min-height: 180px;
          background: #2d271f;
        }
        .table-wrap {
          overflow-x: auto;
          border-radius: 16px;
          border: 1px solid rgba(67, 53, 31, 0.08);
          background: #fffaf3;
          padding: 8px;
        }
        .table > thead > tr > th {
          background: #e9dbc0;
          color: #3b3024;
          border-bottom: none;
        }
        .table > tbody > tr > td {
          background: #fffdf8;
        }
        .inline-actions {
          display: flex;
          flex-wrap: wrap;
          gap: 10px;
          margin-top: 12px;
          margin-bottom: 10px;
        }
        .small-code {
          font-family: Menlo, Consolas, monospace;
          font-size: 12px;
          word-break: break-all;
        }
        .progress-stack {
          display: grid;
          gap: 12px;
        }
        .progress-card {
          background: #fffaf3;
          border: 1px solid rgba(67, 53, 31, 0.08);
          border-radius: 16px;
          padding: 14px 15px;
        }
        .progress-head,
        .progress-foot {
          display: flex;
          justify-content: space-between;
          gap: 12px;
          align-items: center;
        }
        .progress-title {
          font-weight: 700;
          color: #30271f;
        }
        .progress-summary {
          margin: 8px 0 10px;
          color: #5f503f;
        }
        .progress-track {
          position: relative;
          overflow: hidden;
          height: 10px;
          border-radius: 999px;
          background: #eadfcf;
        }
        .progress-fill {
          height: 100%;
          border-radius: 999px;
          transition: width 0.35s ease;
          background: #b58a4a;
        }
        .progress-fill.complete {
          background: #3c7356;
        }
        .progress-fill.failed {
          background: #9b4a36;
        }
        .progress-fill.pending {
          background: #d7c6aa;
        }
        .progress-fill.indeterminate {
          position: relative;
          background:
            linear-gradient(90deg, rgba(181, 138, 74, 0.18) 0%, rgba(181, 138, 74, 0.9) 48%, rgba(181, 138, 74, 0.18) 100%);
          background-size: 220px 100%;
          animation: shimmer 1.3s linear infinite;
        }
        .progress-badge {
          display: inline-flex;
          align-items: center;
          justify-content: center;
          min-width: 72px;
          padding: 4px 10px;
          border-radius: 999px;
          font-size: 11px;
          font-weight: 700;
          letter-spacing: 0.06em;
          text-transform: uppercase;
          background: #eadfcf;
          color: #5a4934;
        }
        .progress-badge.running {
          background: #f1dfbf;
          color: #7a522e;
        }
        .progress-badge.complete {
          background: #dff0e4;
          color: #315944;
        }
        .progress-badge.failed {
          background: #f3d7cf;
          color: #803e2d;
        }
        .progress-badge.pending {
          background: #ece2d2;
          color: #6d5c48;
        }
        @keyframes shimmer {
          from {
            background-position: -220px 0;
          }
          to {
            background-position: 220px 0;
          }
        }
        "
      )
    )
  ),
  div(
    class = "hero",
    div(class = "hero-kicker", "OT-2 Golden Gate | Nanopore | Dorado SUP"),
    h1("Nanopore Amplicon MLST Console"),
    p("Downstream analysis console for Golden Gate-barcoded amplicon libraries generated by the OT-2 workflow, following Oxford Nanopore sequencing and Dorado SUP basecalling."),
    div(
      class = "hero-meta",
      span(paste("Repo:", repo_root)),
      span("Large FASTQ files are referenced by path, not uploaded through the browser"),
      span("One workflow job at a time")
    )
  ),
  tabsetPanel(
    id = "main_tabs",
    tabPanel(
      "Overview",
      fluidRow(
        column(
          width = 4,
          card_panel(
            "Execution Panel",
            p(class = "section-note", "Run Current Workflow uses the parameters and paths currently shown in each tab to execute Step 1-4 in sequence."),
            p(class = "section-note", "Run Demo Workflow switches Step 1-4 to the bundled demo FASTQ and writes outputs under `output/demo_run/`."),
            div(
              class = "inline-actions",
              actionButton("sync_paths", "Sync Suggested Downstream Paths", class = "btn-default"),
              actionButton("run_demo_pipeline", "Run Demo Workflow (Step 1-4)", class = "btn-default"),
              actionButton("run_full_pipeline", "Run Current Workflow (Step 1-4)", class = "btn-primary"),
              actionButton("clear_job_state", "Clear Job State", class = "btn-default")
            ),
            div(class = "command-box", verbatimTextOutput("full_command_preview"))
          ),
          card_panel(
            "Environment Diagnostics",
            p(class = "section-note", "This table checks both the current shell and the Step 2 conda environments. `flye`, `bwa`, `samtools`, and `medaka*` are checked inside the assembly env, and `mlst` is checked inside the MLST env."),
            table_preview_output("tool_diagnostics_table")
          ),
          card_panel(
            "Step Progress",
            uiOutput("workflow_progress_ui")
          )
        ),
        column(
          width = 8,
          card_panel(
            "Job Status",
            uiOutput("job_status_ui")
          ),
          card_panel(
            "Live Log",
            p(class = "section-note", "The log tail refreshes every second while a job is running."),
            div(class = "log-panel", verbatimTextOutput("job_log"))
          )
        )
      ),
      fluidRow(
        column(
          width = 12,
          card_panel(
            "Primary Output Status",
            uiOutput("output_status_ui"),
            div(
              class = "inline-actions",
              downloadButton("download_final_table", "Download Final TSV", class = "btn-default")
            )
          )
        )
      ),
      fluidRow(
        column(
          width = 12,
          card_panel(
            "Result Preview",
            p(class = "section-note", "Only the first few rows are loaded to avoid reading the entire table at once."),
            numericInput("preview_rows", "Preview rows", value = 8, min = 1, max = 50, step = 1),
            tabsetPanel(
              tabPanel("Step 1 summary", table_preview_output("preview_step1")),
              tabPanel("Assembly summary", table_preview_output("preview_step2_assembly")),
              tabPanel("MLST summary", table_preview_output("preview_step2_mlst")),
              tabPanel("ST guess", table_preview_output("preview_step3")),
              tabPanel("Final table", table_preview_output("preview_step4"))
            )
          )
        )
      )
    ),
    tabPanel(
      "Step 1 Demultiplex",
      fluidRow(
        column(
          width = 4,
          card_panel(
            "Parameters",
            textInput("step1_input_dir", "FASTQ input directory", value = "fastq"),
            textInput("step1_output_dir", "Output directory", value = "output/01_demultiplex"),
            numericInput("step1_min_len", "Minimum read length", value = 3000, min = 0, step = 100),
            numericInput("step1_max_len", "Maximum read length", value = 4000, min = 0, step = 100),
            numericInput("step1_primer_mismatch", "primer mismatch", value = 7, min = 0, step = 1),
            numericInput("step1_barcode_mismatch", "barcode mismatch", value = 1, min = 0, step = 1),
            numericInput("step1_pair_gap_slack", "pair gap slack", value = 0, min = 0, step = 1),
            checkboxInput("step1_rescue_partial", "Enable partial barcode rescue", value = FALSE),
            p(class = "section-note", "`pair_gap_slack` is still passed to the script, but the current production version of `01_sup_demultiplex.R` does not actively use this parameter."),
            div(class = "inline-actions", actionButton("run_step1", "Run Step 1", class = "btn-primary"))
          )
        ),
        column(
          width = 8,
          card_panel(
            "Progress",
            uiOutput("step1_progress_ui")
          ),
          card_panel(
            "Command Preview",
            div(class = "command-box", verbatimTextOutput("step1_command_preview"))
          ),
          card_panel(
            "Output Preview",
            table_preview_output("preview_step1_tab")
          )
        )
      )
    ),
    tabPanel(
      "Step 2 Assembly + MLST",
      fluidRow(
        column(
          width = 4,
          card_panel(
            "Parameters",
            textInput("step2_sample_fastq_dir", "Sample FASTQ directory", value = "output/01_demultiplex/merged_by_sample"),
            textInput("step2_assembly_output_dir", "Assembly output directory", value = "output/02_assembly"),
            textInput("step2_mlst_output_dir", "MLST output directory", value = "output/03_mlst"),
            textInput("step2_assembly_env", "assembly conda env", value = "nanopore-assembly-mlst-macos-arm64"),
            textInput("step2_mlst_env", "MLST conda env", value = "mlst-only"),
            numericInput("step2_threads", "assembly threads", value = 4, min = 1, step = 1),
            numericInput("step2_mlst_threads", "MLST threads", value = 1, min = 1, step = 1),
            textInput("step2_genome_size", "genome size", value = "4k"),
            numericInput("step2_asm_coverage", "Flye asm coverage", value = 40, min = 1, step = 1),
            numericInput("step2_min_overlap", "Flye min overlap", value = 1000, min = 1, step = 100),
            textInput("step2_medaka_model", "medaka model", value = "r941_min_sup_g507"),
            textInput("step2_samples", "Restrict sample names, comma-separated", value = ""),
            textInput("step2_scheme", "MLST scheme", value = "saureus"),
            checkboxInput("step2_skip_polish", "Skip polish", value = FALSE),
            checkboxInput("step2_skip_assembly", "Skip assembly", value = FALSE),
            checkboxInput("step2_skip_mlst", "Skip MLST", value = FALSE),
            checkboxInput("step2_medaka_bacteria", "Use medaka --bacteria", value = FALSE),
            div(class = "inline-actions", actionButton("run_step2", "Run Step 2", class = "btn-primary"))
          )
        ),
        column(
          width = 8,
          card_panel(
            "Progress",
            uiOutput("step2_progress_ui")
          ),
          card_panel(
            "Command Preview",
            div(class = "command-box", verbatimTextOutput("step2_command_preview"))
          ),
          card_panel(
            "Assembly summary",
            table_preview_output("preview_step2_assembly_tab")
          ),
          card_panel(
            "MLST summary",
            table_preview_output("preview_step2_mlst_tab")
          )
        )
      )
    ),
    tabPanel(
      "Step 3 ST Guess",
      fluidRow(
        column(
          width = 4,
          card_panel(
            "Parameters",
            textInput("step3_mlst_summary", "MLST summary path", value = "output/03_mlst/mlst_sample_summary.tsv"),
            textInput("step3_profile_table", "Profile table", value = "db/saureus_Profile.txt"),
            textInput("step3_output_dir", "Output directory", value = "output/04_st_guess"),
            numericInput("step3_top_n", "Top N candidates", value = 10, min = 1, step = 1),
            div(class = "inline-actions", actionButton("run_step3", "Run Step 3", class = "btn-primary"))
          )
        ),
        column(
          width = 8,
          card_panel(
            "Progress",
            uiOutput("step3_progress_ui")
          ),
          card_panel(
            "Command Preview",
            div(class = "command-box", verbatimTextOutput("step3_command_preview"))
          ),
          card_panel(
            "ST guess summary",
            table_preview_output("preview_step3_tab")
          )
        )
      )
    ),
    tabPanel(
      "Step 4 Result Table",
      fluidRow(
        column(
          width = 4,
          card_panel(
            "Parameters",
            textInput("step4_assembly_summary", "Assembly summary path", value = "output/02_assembly/assembly_polish_summary.tsv"),
            textInput("step4_mlst_summary", "MLST summary path", value = "output/03_mlst/mlst_sample_summary.tsv"),
            textInput("step4_st_guess_summary", "ST guess summary path", value = "output/04_st_guess/st_guess_summary.tsv"),
            textInput("step4_output_path", "Output TSV", value = "output/Result_Table/SUP_final_publication_style_mlst.tsv"),
            div(class = "inline-actions", actionButton("run_step4", "Run Step 4", class = "btn-primary"))
          )
        ),
        column(
          width = 8,
          card_panel(
            "Progress",
            uiOutput("step4_progress_ui")
          ),
          card_panel(
            "Command Preview",
            div(class = "command-box", verbatimTextOutput("step4_command_preview"))
          ),
          card_panel(
            "Result Table",
            p(class = "section-note", "This tab loads the complete final TSV instead of the preview row limit used elsewhere."),
            table_preview_output("preview_step4_tab")
          )
        )
      )
    )
  )
)

server <- function(input, output, session) {
  rv <- reactiveValues(
    active = FALSE,
    label = NULL,
    pid = NA_integer_,
    log_path = "",
    status_path = "",
    started_at = NULL,
    finished_at = NULL,
    last_exit = NA_integer_
  )

  current_values <- reactive({
    list(
      step1_input_dir = input$step1_input_dir,
      step1_output_dir = input$step1_output_dir,
      step1_min_len = input$step1_min_len,
      step1_max_len = input$step1_max_len,
      step1_primer_mismatch = input$step1_primer_mismatch,
      step1_barcode_mismatch = input$step1_barcode_mismatch,
      step1_pair_gap_slack = input$step1_pair_gap_slack,
      step1_rescue_partial = input$step1_rescue_partial,
      step2_sample_fastq_dir = input$step2_sample_fastq_dir,
      step2_assembly_output_dir = input$step2_assembly_output_dir,
      step2_mlst_output_dir = input$step2_mlst_output_dir,
      step2_assembly_env = input$step2_assembly_env,
      step2_mlst_env = input$step2_mlst_env,
      step2_threads = input$step2_threads,
      step2_mlst_threads = input$step2_mlst_threads,
      step2_genome_size = input$step2_genome_size,
      step2_asm_coverage = input$step2_asm_coverage,
      step2_min_overlap = input$step2_min_overlap,
      step2_medaka_model = input$step2_medaka_model,
      step2_samples = input$step2_samples,
      step2_scheme = input$step2_scheme,
      step2_skip_polish = input$step2_skip_polish,
      step2_skip_assembly = input$step2_skip_assembly,
      step2_skip_mlst = input$step2_skip_mlst,
      step2_medaka_bacteria = input$step2_medaka_bacteria,
      step3_mlst_summary = input$step3_mlst_summary,
      step3_profile_table = input$step3_profile_table,
      step3_output_dir = input$step3_output_dir,
      step3_top_n = input$step3_top_n,
      step4_assembly_summary = input$step4_assembly_summary,
      step4_mlst_summary = input$step4_mlst_summary,
      step4_st_guess_summary = input$step4_st_guess_summary,
      step4_output_path = input$step4_output_path
    )
  })

  progress_states <- reactive({
    invalidateLater(1000, session)
    build_progress_states(current_values(), reactiveValuesToList(rv))
  })

  output_paths <- reactive({
    list(
      "Step 1 file summary" = resolve_path(file.path(input$step1_output_dir, "file_summary.tsv")),
      "Merged sample FASTQ directory" = resolve_path(file.path(input$step1_output_dir, "merged_by_sample")),
      "Assembly summary" = resolve_path(file.path(input$step2_assembly_output_dir, "assembly_polish_summary.tsv")),
      "MLST summary" = resolve_path(file.path(input$step2_mlst_output_dir, "mlst_sample_summary.tsv")),
      "ST guess summary" = resolve_path(file.path(input$step3_output_dir, "st_guess_summary.tsv")),
      "Final result table" = resolve_path(input$step4_output_path)
    )
  })

  launch_job <- function(label, shell_command) {
    if (isTRUE(rv$active)) {
      showNotification("A job is already running. Wait for the current job to finish first.", type = "warning")
      return(invisible(FALSE))
    }

    job <- tryCatch(
      launch_background_job(shell_command, label),
      error = function(error) {
        showNotification(conditionMessage(error), type = "error", duration = NULL)
        NULL
      }
    )

    if (is.null(job)) {
      return(invisible(FALSE))
    }

    rv$active <- TRUE
    rv$label <- label
    rv$pid <- job$pid
    rv$log_path <- job$log
    rv$status_path <- job$status
    rv$started_at <- Sys.time()
    rv$finished_at <- NULL
    rv$last_exit <- NA_integer_

    if (is.na(job$pid)) {
      showNotification(sprintf("%s started.", label), type = "message")
    } else {
      showNotification(sprintf("%s started, PID %s", label, job$pid), type = "message")
    }
    invisible(TRUE)
  }

  observe({
    invalidateLater(1000, session)
    if (!isTRUE(rv$active)) {
      return()
    }

    if (!file.exists(rv$status_path)) {
      return()
    }

    status_code <- read_status_code(rv$status_path)
    rv$active <- FALSE
    rv$finished_at <- Sys.time()
    rv$last_exit <- status_code

    if (identical(status_code, 0L)) {
      showNotification(sprintf("%s completed.", rv$label), type = "message")
    } else {
      showNotification(sprintf("%s failed. Check the log.", rv$label), type = "error", duration = NULL)
    }
  })

  observeEvent(input$sync_paths, {
    updateTextInput(
      session,
      "step2_sample_fastq_dir",
      value = file.path(input$step1_output_dir, "merged_by_sample")
    )
    updateTextInput(
      session,
      "step3_mlst_summary",
      value = file.path(input$step2_mlst_output_dir, "mlst_sample_summary.tsv")
    )
    updateTextInput(
      session,
      "step4_assembly_summary",
      value = file.path(input$step2_assembly_output_dir, "assembly_polish_summary.tsv")
    )
    updateTextInput(
      session,
      "step4_mlst_summary",
      value = file.path(input$step2_mlst_output_dir, "mlst_sample_summary.tsv")
    )
    updateTextInput(
      session,
      "step4_st_guess_summary",
      value = file.path(input$step3_output_dir, "st_guess_summary.tsv")
    )
  })

  observeEvent(input$clear_job_state, {
    rv$active <- FALSE
    rv$label <- NULL
    rv$pid <- NA_integer_
    rv$log_path <- ""
    rv$status_path <- ""
    rv$started_at <- NULL
    rv$finished_at <- NULL
    rv$last_exit <- NA_integer_
  })

  observeEvent(input$run_step1, {
    launch_job("Step 1 demultiplex", build_step1_command(current_values()))
  })

  observeEvent(input$run_step2, {
    launch_job("Step 2 assembly_mlst", build_step2_command(current_values()))
  })

  observeEvent(input$run_step3, {
    launch_job("Step 3 st_guess", build_step3_command(current_values()))
  })

  observeEvent(input$run_step4, {
    launch_job("Step 4 result_table", build_step4_command(current_values()))
  })

  observeEvent(input$run_full_pipeline, {
    if (isTRUE(input$step2_skip_mlst)) {
      showNotification("Full workflow mode requires Step 2 to produce the MLST summary. Disable `Skip MLST` first.", type = "warning")
      return()
    }
    launch_job("Full workflow", build_full_command(current_values(), pretty = FALSE))
  })

  observeEvent(input$run_demo_pipeline, {
    demo_fastq <- resolve_path(file.path("demo", "Samlst_SUP.demo.fastq.gz"))
    if (!file.exists(demo_fastq)) {
      showNotification(sprintf("Demo FASTQ not found: %s", demo_fastq), type = "error", duration = NULL)
      return()
    }
    if (isTRUE(input$step2_skip_mlst)) {
      showNotification("Demo mode requires Step 2 to produce the MLST summary. Disable `Skip MLST` first.", type = "warning")
      return()
    }

    demo_values <- demo_values_from(current_values())

    updateTextInput(session, "step1_input_dir", value = demo_values$step1_input_dir)
    updateTextInput(session, "step1_output_dir", value = demo_values$step1_output_dir)
    updateTextInput(session, "step2_sample_fastq_dir", value = demo_values$step2_sample_fastq_dir)
    updateTextInput(session, "step2_assembly_output_dir", value = demo_values$step2_assembly_output_dir)
    updateTextInput(session, "step2_mlst_output_dir", value = demo_values$step2_mlst_output_dir)
    updateTextInput(session, "step3_mlst_summary", value = demo_values$step3_mlst_summary)
    updateTextInput(session, "step3_output_dir", value = demo_values$step3_output_dir)
    updateTextInput(session, "step4_assembly_summary", value = demo_values$step4_assembly_summary)
    updateTextInput(session, "step4_mlst_summary", value = demo_values$step4_mlst_summary)
    updateTextInput(session, "step4_st_guess_summary", value = demo_values$step4_st_guess_summary)
    updateTextInput(session, "step4_output_path", value = demo_values$step4_output_path)

    launch_job("Demo full workflow", build_full_command(demo_values, pretty = FALSE))
  })

  output$full_command_preview <- renderText({
    build_full_command(current_values(), pretty = TRUE)
  })

  output$step1_command_preview <- renderText({
    build_step1_command(current_values())
  })

  output$step2_command_preview <- renderText({
    build_step2_command(current_values())
  })

  output$step3_command_preview <- renderText({
    build_step3_command(current_values())
  })

  output$step4_command_preview <- renderText({
    build_step4_command(current_values())
  })

  output$workflow_progress_ui <- renderUI({
    states <- progress_states()
    div(
      class = "progress-stack",
      progress_panel_ui("Step 1 Demultiplex", states$step1),
      progress_panel_ui("Step 2 Assembly + MLST", states$step2),
      progress_panel_ui("Step 3 ST Guess", states$step3),
      progress_panel_ui("Step 4 Result Table", states$step4)
    )
  })

  output$step1_progress_ui <- renderUI({
    progress_panel_ui("Step 1 Demultiplex", progress_states()$step1)
  })

  output$step2_progress_ui <- renderUI({
    progress_panel_ui("Step 2 Assembly + MLST", progress_states()$step2)
  })

  output$step3_progress_ui <- renderUI({
    progress_panel_ui("Step 3 ST Guess", progress_states()$step3)
  })

  output$step4_progress_ui <- renderUI({
    progress_panel_ui("Step 4 Result Table", progress_states()$step4)
  })

  output$job_status_ui <- renderUI({
    banner_class <- "idle"
    title <- "No job is currently running"
    lines <- c(
      sprintf("Repo root: %s", repo_root)
    )

    if (isTRUE(rv$active)) {
      banner_class <- "live"
      title <- sprintf("%s running", rv$label %||% "Workflow")
      lines <- c(sprintf("Started: %s", format(rv$started_at, "%Y-%m-%d %H:%M:%S")), sprintf("Log: %s", rv$log_path))
    } else if (!is.na(rv$last_exit)) {
      banner_class <- if (identical(rv$last_exit, 0L)) "success" else "failed"
      title <- if (identical(rv$last_exit, 0L)) {
        sprintf("%s completed", rv$label %||% "Workflow")
      } else {
        sprintf("%s failed", rv$label %||% "Workflow")
      }
      lines <- c(
        sprintf("Exit status: %s", rv$last_exit),
        sprintf("Started: %s", format(rv$started_at, "%Y-%m-%d %H:%M:%S")),
        sprintf("Finished: %s", format(rv$finished_at, "%Y-%m-%d %H:%M:%S")),
        sprintf("Log: %s", rv$log_path)
      )
    }

    div(
      class = paste("job-banner", banner_class),
      h4(title),
      tags$pre(class = "small-code", paste(lines, collapse = "\n"))
    )
  })

  output$job_log <- renderText({
    invalidateLater(1000, session)
    read_log_tail(rv$log_path)
  })

  output$output_status_ui <- renderUI({
    invalidateLater(1500, session)
    tags$div(
      class = "status-grid",
      lapply(names(output_paths()), function(label) {
        path <- output_paths()[[label]]
        ready <- file.exists(path) || dir.exists(path)
        status_chip(label, path, ready)
      })
    )
  })

  output$tool_diagnostics_table <- renderTable({
    tool_diagnostics(input$step2_assembly_env, input$step2_mlst_env)
  }, striped = TRUE, bordered = TRUE, spacing = "xs", rownames = FALSE)

  output$preview_step1 <- renderTable({
    invalidateLater(1500, session)
    read_preview_tsv(resolve_path(file.path(input$step1_output_dir, "file_summary.tsv")), input$preview_rows)
  }, striped = TRUE, bordered = TRUE, spacing = "xs", rownames = FALSE)

  output$preview_step2_assembly <- renderTable({
    invalidateLater(1500, session)
    read_preview_tsv(resolve_path(file.path(input$step2_assembly_output_dir, "assembly_polish_summary.tsv")), input$preview_rows)
  }, striped = TRUE, bordered = TRUE, spacing = "xs", rownames = FALSE)

  output$preview_step2_mlst <- renderTable({
    invalidateLater(1500, session)
    read_preview_tsv(resolve_path(file.path(input$step2_mlst_output_dir, "mlst_sample_summary.tsv")), input$preview_rows)
  }, striped = TRUE, bordered = TRUE, spacing = "xs", rownames = FALSE)

  output$preview_step3 <- renderTable({
    invalidateLater(1500, session)
    read_preview_tsv(resolve_path(file.path(input$step3_output_dir, "st_guess_summary.tsv")), input$preview_rows)
  }, striped = TRUE, bordered = TRUE, spacing = "xs", rownames = FALSE)

  output$preview_step4 <- renderTable({
    invalidateLater(1500, session)
    read_preview_tsv(resolve_path(input$step4_output_path), input$preview_rows)
  }, striped = TRUE, bordered = TRUE, spacing = "xs", rownames = FALSE)

  output$preview_step1_tab <- renderTable({
    invalidateLater(1500, session)
    read_preview_tsv(resolve_path(file.path(input$step1_output_dir, "file_summary.tsv")), input$preview_rows)
  }, striped = TRUE, bordered = TRUE, spacing = "xs", rownames = FALSE)

  output$preview_step2_assembly_tab <- renderTable({
    invalidateLater(1500, session)
    read_preview_tsv(resolve_path(file.path(input$step2_assembly_output_dir, "assembly_polish_summary.tsv")), input$preview_rows)
  }, striped = TRUE, bordered = TRUE, spacing = "xs", rownames = FALSE)

  output$preview_step2_mlst_tab <- renderTable({
    invalidateLater(1500, session)
    read_preview_tsv(resolve_path(file.path(input$step2_mlst_output_dir, "mlst_sample_summary.tsv")), input$preview_rows)
  }, striped = TRUE, bordered = TRUE, spacing = "xs", rownames = FALSE)

  output$preview_step3_tab <- renderTable({
    invalidateLater(1500, session)
    read_preview_tsv(resolve_path(file.path(input$step3_output_dir, "st_guess_summary.tsv")), input$preview_rows)
  }, striped = TRUE, bordered = TRUE, spacing = "xs", rownames = FALSE)

  output$preview_step4_tab <- renderTable({
    invalidateLater(1500, session)
    read_full_tsv(resolve_path(input$step4_output_path))
  }, striped = TRUE, bordered = TRUE, spacing = "xs", rownames = FALSE)

  output$download_final_table <- downloadHandler(
    filename = function() {
      basename(resolve_path(input$step4_output_path))
    },
    content = function(file) {
      final_path <- resolve_path(input$step4_output_path)
      if (!file.exists(final_path)) {
        stop("Final table does not exist yet.", call. = FALSE)
      }
      file.copy(final_path, file, overwrite = TRUE)
    }
  )
}

app <- shinyApp(ui = ui, server = server)

if (sys.nframe() == 0) {
  shiny::runApp(app, launch.browser = TRUE)
} else {
  app
}
