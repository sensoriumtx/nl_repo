#!/usr/bin/env Rscript

suppressMessages(library(tidyverse))
suppressMessages(library(lubridate))
suppressMessages(library(parallel))

cat("Script starting...\n", flush = TRUE)

# ------------------------- Logging -------------------------
args <- commandArgs(trailingOnly = TRUE)
log_file <- NULL
log <- function(message) {
  cat(message, "\n", flush = TRUE)
  if (!is.null(log_file)) {
    cat(message, "\n", file = log_file, append = TRUE)
  }
}
log(paste("Args:", paste(args, collapse = " | ")))

# ------------------------- Argument Parsing -------------------------
parseArgs <- function(args) {
  out <- list(endpoint = "dev")
  i <- 1
  while (i <= length(args)) {
    arg <- args[i]
    val <- if (i + 1 <= length(args)) args[i + 1] else NA
    if (is.na(val)) stop(paste("Missing value for", arg))
    if (arg == "--cmp") out$cmp <- val
    else if (arg == "--smiles") out$smiles <- val
    else if (arg == "--inchi") out$inchi <- val
    else if (arg == "--inchi_key") out$inchi_key <- val
    else if (arg == "--iupac") out$iupac <- val
    else if (arg == "--cmp_label") out$cmp_label <- val
    else if (arg == "--in_file") out$in_file <- val
    else if (arg == "--endpoint") out$endpoint <- val
    else if (arg == "--outdir") out$outdir <- val
    else stop(paste("Unknown argument:", arg))
    i <- i + 2
  }
  return(out)
}

params <- parseArgs(args)
if (is.null(params$outdir)) stop("--outdir must be specified")

dir.create(params$outdir, showWarnings = FALSE, recursive = TRUE)
log_file <- file.path(params$outdir, "pipeline_log.txt")
cat("Pipeline Log\n==============\n", file = log_file)

# ------------------------- Step 0: Flexible Input Matching -------------------------
if (!is.null(params$in_file)) {
  log("[Step 0] Resolving input via flexible matching from master file")
  master_df <- read_csv(params$in_file, show_col_types = FALSE)
  matches <- master_df

  if (!is.null(params$cmp_label)) {
    matches <- matches %>% filter(str_detect(tolower(cmp_label), tolower(params$cmp_label)))
  }
  if (!is.null(params$smiles)) {
    matches <- matches %>% filter(str_detect(tolower(primarySmiles), tolower(params$smiles)))
  }
  if (!is.null(params$inchi)) {
    matches <- matches %>% filter(str_detect(tolower(inchi), tolower(params$inchi)))
  }
  if (!is.null(params$inchi_key)) {
    matches <- matches %>% filter(str_detect(tolower(inchi_key), tolower(params$inchi_key)))
  }
  if (!is.null(params$cmp)) {
    matches <- matches %>% filter(str_detect(tolower(cmp), tolower(params$cmp)))
  }
  if (!is.null(params$inpac)) {
    matches <- matches %>% filter(str_detect(tolower(iupac_name), tolower(params$iupac)))
  }

  resolved_cmp_table <- matches %>% distinct(cmp, .keep_all = TRUE)
  resolved_cmp_file <- file.path(params$outdir, "step0_resolved_from_input.csv")
  write_csv(resolved_cmp_table, resolved_cmp_file)
  log(paste("[Step 0] Matching complete. Found", nrow(resolved_cmp_table), "matches"))
  log(paste("Output:", resolved_cmp_file))

  input_file_path <- resolved_cmp_file
} else {
  log("[Step 0] Skipped: no --in_file provided")
  resolve_input <- tibble(
    cmp = params$cmp,
    primarySmiles = params$smiles,
    inchi = params$inchi,
    inchi_key = params$inchi_key,
    iupac_name = params$iupac,
    cmp_label = params$cmp_label
  )
  input_file_path <- file.path(params$outdir, "step1_identifier_input.csv")
  write_csv(resolve_input, input_file_path)
}

# ------------------------- Step 1: Resolve Identifiers to CMP -------------------------
log("[Step 1] Resolving compound identifiers to `cmp` ID")

resolved_cmp_file <- file.path(params$outdir, "step1_resolved_cmp.csv")
resolve_cmd <- paste(
  "Rscript scripts/resolve_cmp_id.r",
  "--in", shQuote(input_file_path),
  "--out", shQuote(resolved_cmp_file)
)
system(resolve_cmd)
if (!file.exists(resolved_cmp_file)) stop("Step 1 failed: cmp resolution output not found.")
log("[Step 1] Complete")

cmp_ids <- read_csv(resolved_cmp_file, show_col_types = FALSE)$cmp %>% unique() %>% paste(collapse = "|")
if (is.null(cmp_ids) || cmp_ids == "") stop("No resolved cmp IDs")

# ------------------------- Step 2: Pull Plants for CMP -------------------------
log("[Step 2] Pulling plants associated with resolved compounds")

plants_file <- file.path(params$outdir, "step2_plants_for_cmp.csv")
plant_cmd <- paste(
  "Rscript scripts/pull_plant_for_cmp_id.r",
  "--endpoint", params$endpoint,
  "--compound_activity_file", shQuote(resolved_cmp_file),
  "--cmp_id_column", "cmp",
  "--out", shQuote(plants_file)
)
system(plant_cmd)
if (!file.exists(plants_file)) stop("Step 2 failed: plant output not found.")
log("[Step 2] Complete")

# ------------------------- Step 3: Pull Acts for CMP -------------------------
log("[Step 3] Pulling activities associated with compounds")

cmp_acts_file <- file.path(params$outdir, "step3_acts_for_cmp.csv")
act_cmd <- paste(
  "Rscript scripts/pull_acts_for_cmp_id.r",
  "--endpoint", params$endpoint,
  "--compound", shQuote(cmp_ids),
  "--out", shQuote(cmp_acts_file)
)
system(act_cmd)
if (!file.exists(cmp_acts_file)) stop("Step 3 failed: cmp activities output not found.")
log("[Step 3] Complete")

# ------------------------- Step 4: Pull Acts for PLN -------------------------
log("[Step 4] Pulling activities associated with plants")

plant_labels <- read_csv(plants_file, show_col_types = FALSE)$pln_label %>%
  unique() %>% na.omit() %>% paste(collapse = "|")

plant_acts_file <- file.path(params$outdir, "step4_acts_for_pln.csv")
pln_act_cmd <- paste(
  "Rscript scripts/pull_acts_for_pln_id.r",
  "--endpoint", params$endpoint,
  "--plants", shQuote(plant_labels),
  "--out", shQuote(plant_acts_file)
)
system(pln_act_cmd)
if (!file.exists(plant_acts_file)) stop("Step 4 failed: plant activities output not found.")
log("[Step 4] Complete")

# ------------------------- Completion -------------------------
log("[Pipeline] Success. Outputs written to:")
log(paste("  - Resolved CMPs:", resolved_cmp_file))
log(paste("  - Plants for CMPs:", plants_file))
log(paste("  - Activities for CMPs:", cmp_acts_file))
log(paste("  - Activities for Plants:", plant_acts_file))
