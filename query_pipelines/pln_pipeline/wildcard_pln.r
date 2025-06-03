#!/usr/bin/env Rscript

# Nick Laskowski
# Version 1.0
# Start from PLN values and resolve downstream data

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
    if (arg == "--plants") out$plants <- val
    else if (arg == "--endpoint") out$endpoint <- val
    else if (arg == "--outdir") out$outdir <- val
    else stop(paste("Unknown argument:", arg))
    i <- i + 2
  }
  return(out)
}

params <- parseArgs(args)
if (is.null(params$plants)) stop("--plants must be provided (e.g., 'Panax ginseng|Withania somnifera')")
if (is.null(params$outdir)) stop("--outdir must be specified")

dir.create(params$outdir, showWarnings = FALSE, recursive = TRUE)
log_file <- file.path(params$outdir, "pipeline_log.txt")
cat("Pipeline Log\n==============\n", file = log_file)

# ------------------------- Step 1: Pull Compounds for PLN -------------------------
log("[Step 1] Pulling compounds associated with provided plants")

compounds_file <- file.path(params$outdir, "step1_compounds_for_pln.csv")
cmp_cmd <- paste(
  "Rscript scripts/pull_cmp_for_pln.r",
  "--endpoint", params$endpoint,
  "--plants", shQuote(params$plants),
  "--out", shQuote(compounds_file)
)
system(cmp_cmd)
if (!file.exists(compounds_file)) stop("Step 1 failed: compound output not found.")
log("[Step 1] Complete")

# ------------------------- Step 2: Pull Activities for CMP -------------------------
log("[Step 2] Pulling activities associated with compounds")

cmp_df <- read_csv(compounds_file, show_col_types = FALSE)
if (!"cmp" %in% colnames(cmp_df)) stop("The column 'cmp' must be present in compounds file.")
cmp_ids <- cmp_df$cmp %>% unique() %>% na.omit() %>% sort() %>% paste(collapse = "|")
if (cmp_ids == "") stop("No valid cmp IDs found from compound list.")

cmp_acts_file <- file.path(params$outdir, "step2_acts_for_cmp.csv")
act_cmd <- paste(
  "Rscript scripts/pull_acts_for_specific_cmp_ids.r",
  "--endpoint", params$endpoint,
  "--compound", shQuote(cmp_ids),
  "--out", shQuote(cmp_acts_file)
)
system(act_cmd)
if (!file.exists(cmp_acts_file)) stop("Step 2 failed: compound activity output not found.")
log("[Step 2] Complete")

# ------------------------- Step 3: Pull Activities for PLN -------------------------
log("[Step 3] Pulling activities associated with provided plants")

plant_acts_file <- file.path(params$outdir, "step3_acts_for_pln.csv")
pln_act_cmd <- paste(
  "Rscript scripts/pull_acts_for_specific_pln.r",
  "--endpoint", params$endpoint,
  "--plants", shQuote(params$plants),
  "--out", shQuote(plant_acts_file)
)
system(pln_act_cmd)
if (!file.exists(plant_acts_file)) stop("Step 3 failed: plant activity output not found.")
log("[Step 3] Complete")

# ------------------------- Completion -------------------------
log("[Pipeline] Success. Outputs written to:")
log(paste("  - Compounds for PLNs:", compounds_file))
log(paste("  - Activities for CMPs:", cmp_acts_file))
log(paste("  - Activities for PLNs:", plant_acts_file))
