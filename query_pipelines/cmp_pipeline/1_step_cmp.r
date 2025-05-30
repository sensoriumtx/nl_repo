#!/usr/bin/env Rscript

# ------------------------- Initialization -------------------------
suppressMessages(library(tidyverse))
suppressMessages(library(lubridate))
suppressMessages(library(parallel))

cat("Script starting...\n", flush = TRUE)

args <- commandArgs(trailingOnly = TRUE)
log_file <- NULL
log <- function(message) {
  cat(message, "\n", flush = TRUE)
  if (!is.null(log_file)) {
    cat(paste0(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), message, "\n"), file = log_file, append = TRUE)
  }
}

# ------------------------- Argument Parsing -------------------------
parseArgs <- function(args) {
  out <- list(endpoint = "dev")
  i <- 1
  while (i <= length(args)) {
    arg <- args[i]
    val <- if (i + 1 <= length(args)) args[i + 1] else NA
    if (is.na(val)) stop(paste("Missing value for", arg))
    if (arg == "--cmp") out$cmp <- val
    else if (arg == "--cmp_in_file") out$cmp_in_file <- val
    else if (arg == "--search") out$search <- val
    else if (arg == "--scoring") out$scoring <- val
    else if (arg == "--endpoint") out$endpoint <- val
    else if (arg == "--outdir") out$outdir <- val
    else stop(paste("Unknown argument:", arg))
    i <- i + 2
  }
  return(out)
}

params <- parseArgs(args)
if (is.null(params$outdir)) stop("❌ --outdir must be specified")

dir.create(params$outdir, showWarnings = FALSE, recursive = TRUE)
log_file <- file.path(params$outdir, "pipeline_log.txt")
cat("Pipeline Log\n==============\n", file = log_file)
log(paste("Args:", paste(args, collapse = " | ")))

# ------------------------- Step 0: Resolve CMPs -------------------------
log("[Step 0] Resolving input compound identifiers")
cmp_ids <- NULL
resolved_cmp_file <- NULL

if (!is.null(params$cmp_in_file)) {
  log("→ Using cmp values from --cmp_in_file")
  cmp_df <- read_csv(params$cmp_in_file, show_col_types = FALSE)
  if (!"cmp" %in% names(cmp_df)) stop("❌ The file provided to --cmp_in_file must contain a 'cmp' column.")
  cmp_ids <- cmp_df$cmp %>% unique() %>% na.omit() %>% trimws() %>% sort()
  resolved_cmp_file <- params$cmp_in_file

} else if (!is.null(params$cmp)) {
  log("→ Using cmp values from --cmp string")
  cmp_ids <- str_split(params$cmp, "\\|")[[1]] %>% unique() %>% na.omit() %>% trimws() %>% sort()
  resolved_cmp_file <- file.path(params$outdir, "step0_cmp_from_arg.csv")
  write_csv(tibble(cmp = cmp_ids), resolved_cmp_file)

} else if (!is.null(params$search) && !is.null(params$scoring)) {
  log("→ Performing grep-style search using --search across --scoring file")
  search_term <- tolower(params$search)
  master_df <- read_csv(params$scoring, show_col_types = FALSE)
  char_cols <- master_df %>% select(where(is.character)) %>% names()
  matches <- master_df %>%
    filter(if_any(all_of(char_cols), ~ str_detect(tolower(.), fixed(search_term, ignore_case = TRUE))))
  if (nrow(matches) == 0) stop(paste("❌ No matches found for search term:", search_term))
  resolved_cmp_file <- file.path(params$outdir, "step0_grep_matched_rows.csv")
  write_csv(matches, resolved_cmp_file)
  cmp_ids <- matches$cmp %>% unique() %>% na.omit() %>% trimws() %>% sort()
  log(paste("✓ Found", length(cmp_ids), "unique cmp matches from grep search"))

} else {
  stop("❌ No valid input provided. Use --cmp_in_file, --cmp, or --search + --scoring.")
}

if (length(cmp_ids) == 0) stop("❌ No compound identifiers resolved.")
log(paste("✓ Final resolved cmp count:", length(cmp_ids)))

# ------------------------- Utility: Chunking Function -------------------------
run_chunked_step <- function(input_ids, chunk_size, step_prefix, outdir, step_script, arg_flag, endpoint) {
  chunk_dir <- file.path(outdir, paste0(step_prefix, "_chunks"))
  dir.create(chunk_dir, showWarnings = FALSE, recursive = TRUE)
  input_chunks <- split(input_ids, ceiling(seq_along(input_ids) / chunk_size))
  chunk_files <- list()

  for (i in seq_along(input_chunks)) {
    chunk_vals <- input_chunks[[i]]
    chunk_str <- paste(chunk_vals, collapse = "|")
    chunk_out <- file.path(chunk_dir, sprintf("%s_chunk_%03d.csv", step_prefix, i))
    cmd <- paste(
      "Rscript", step_script,
      "--endpoint", endpoint,
      arg_flag, shQuote(chunk_str),
      "--out", shQuote(chunk_out)
    )
    log(paste("→ Running", step_prefix, "chunk", i, "→", chunk_out))
    exit_code <- system(cmd)
    if (exit_code != 0 || !file.exists(chunk_out)) {
      stop(paste("❌ Chunk", i, "failed or output not found at:", chunk_out))
    }
    chunk_files[[i]] <- chunk_out
  }

  final_out <- file.path(outdir, paste0(step_prefix, "_final.csv"))
  merged_df <- bind_rows(lapply(chunk_files, read_csv, show_col_types = FALSE))
  write_csv(merged_df, final_out)
  log(paste("✓", step_prefix, "complete. Final output →", final_out))
  return(final_out)
}

# ------------------------- Step 1: Pull Plants for CMP -------------------------
log("[Step 1] Pulling plants for compounds (chunked)")
step1_out <- run_chunked_step(
  input_ids = cmp_ids,
  chunk_size = 250,
  step_prefix = "step1_plants",
  outdir = params$outdir,
  step_script = "scripts/pull_plant_for_compound_ids.r",
  arg_flag = "--compound",
  endpoint = params$endpoint
)

# ------------------------- Step 2: Pull Acts for CMP -------------------------
log("[Step 2] Pulling activities for compounds (chunked)")
step2_out <- run_chunked_step(
  input_ids = cmp_ids,
  chunk_size = 250,
  step_prefix = "step2_acts",
  outdir = params$outdir,
  step_script = "scripts/pull_acts_for_compound_ids.r",
  arg_flag = "--compound",
  endpoint = params$endpoint
)

# ------------------------- Step 3: Pull Acts for PLNs -------------------------
log("[Step 3] Pulling activities for plants (chunked)")
pln_df <- read_csv(step1_out, show_col_types = FALSE)
if (!"pln_label" %in% names(pln_df)) stop("❌ step1 output missing 'pln_label' column.")
pln_ids <- pln_df$pln_label %>% unique() %>% na.omit()
if (length(pln_ids) == 0) stop("❌ No plant identifiers found from Step 1.")

step3_out <- run_chunked_step(
  input_ids = pln_ids,
  chunk_size = 250,
  step_prefix = "step3_acts",
  outdir = params$outdir,
  step_script = "scripts/pull_acts_for_specific_pln.r",
  arg_flag = "--plants",
  endpoint = params$endpoint
)

# ------------------------- Completion -------------------------
log("[Pipeline] ✅ Success. Outputs written to:")
log(paste("  - Resolved CMPs:", resolved_cmp_file))
log(paste("  - Plants for CMPs:", step1_out))
log(paste("  - Activities for CMPs:", step2_out))
log(paste("  - Activities for Plants:", step3_out))
