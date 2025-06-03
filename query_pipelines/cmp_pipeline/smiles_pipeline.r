# Nick Laskowski
# Version 1.4
# SMILES-Based Semantic Association Pipeline — With Proper Stepwise Output Management and Internal Join Handling

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
    if (arg == "--smiles") out$smiles <- val
    else if (arg == "--smiles_file") out$smiles_file <- val
    else if (arg == "--smiles_column") out$smiles_column <- val
    else if (arg == "--endpoint") out$endpoint <- val
    else if (arg == "--outdir") out$outdir <- val
    else stop(paste("Unknown argument:", arg))
    i <- i + 2
  }
  return(out)
}

params <- parseArgs(args)
if (is.null(params$outdir)) stop("--outdir must be specified")

# ------------------------- Directory Setup -------------------------
dir.create(params$outdir, showWarnings = FALSE, recursive = TRUE)
log_file <- file.path(params$outdir, "pipeline_log.txt")
cat("Pipeline Log\n==============\n", file = log_file)

# ------------------------- Step 0: SMILES to Activities -------------------------
log("[Step 0] Pulling activities directly for SMILES")

resolved_cmp_file <- file.path(params$outdir, "step0_acts_for_smiles.csv")

# Determine whether to use direct input or file input
if (!is.null(params$smiles_file)) {
  smiles_data <- read_csv(params$smiles_file, show_col_types = FALSE)
  if (is.null(params$smiles_column)) stop("--smiles_column must be specified if --smiles_file is used")
  if (!params$smiles_column %in% names(smiles_data)) stop(paste("Column", params$smiles_column, "not found in SMILES file"))
  params$smiles <- smiles_data[[params$smiles_column]] %>% unique() %>% na.omit() %>% paste(collapse = "|")
}

if (is.null(params$smiles)) stop("You must provide either --smiles or a valid --smiles_file and --smiles_column")

pull_smiles_cmd <- paste(
  "Rscript Nick_dev/sensgit/scripts/pull_act_for_smiles.r",
  "--endpoint", params$endpoint,
  "--smiles", shQuote(params$smiles),
  "--out", shQuote(resolved_cmp_file)
)

system(pull_smiles_cmd)
if (!file.exists(resolved_cmp_file)) stop("Step 0 failed: SMILES activities output not found.")
log(paste("[Step 0] Complete. Output written to:", resolved_cmp_file))

# ------------------------- Step 1: Pull Plants for CMP -------------------------
log("[Step 1] Pulling plants associated with resolved compounds")

plants_file <- file.path(params$outdir, "step1_plants_for_cmp.csv")
plant_cmd <- paste(
  "Rscript scripts/pull_plant_for_compound_ids.r",
  "--endpoint", params$endpoint,
  "--compound_activity_file", shQuote(resolved_cmp_file),
  "--cmp_id_column", "cmp",
  "--out", shQuote(plants_file)
)
system(plant_cmd)
if (!file.exists(plants_file)) stop("Step 1 failed: plant output not found.")
log(paste("[Step 1] Complete. Output written to:", plants_file))

# ------------------------- Step 2: Pull Acts for CMP -------------------------
log("[Step 2] Pulling activities associated with compounds")

resolved_df <- read_csv(resolved_cmp_file, show_col_types = FALSE)
if (!"cmp" %in% colnames(resolved_df)) stop("Column 'cmp' not found in SMILES output file.")

cmp_ids <- resolved_df$cmp %>% unique() %>% na.omit() %>% paste(collapse = "|")

cmp_acts_file <- file.path(params$outdir, "step2_acts_for_cmp.csv")
act_cmd <- paste(
  "Rscript scripts/pull_acts_for_specific_cmp_ids.r",
  "--endpoint", params$endpoint,
  "--compound", shQuote(cmp_ids),
  "--out", shQuote(cmp_acts_file)
)
system(act_cmd)
if (!file.exists(cmp_acts_file)) stop("Step 2 failed: cmp activities output not found.")
log(paste("[Step 2] Complete. Output written to:", cmp_acts_file))

# ------------------------- Step 3: Pull Acts for PLN -------------------------
log("[Step 3] Pulling activities associated with plants")

plant_labels <- read_csv(plants_file, show_col_types = FALSE)$pln_label %>%
  unique() %>% na.omit() %>% paste(collapse = "|")

plant_acts_file <- file.path(params$outdir, "step3_acts_for_pln.csv")
pln_act_cmd <- paste(
  "Rscript scripts/pull_acts_for_specific_pln.r",
  "--endpoint", params$endpoint,
  "--plants", shQuote(plant_labels),
  "--out", shQuote(plant_acts_file)
)
system(pln_act_cmd)
if (!file.exists(plant_acts_file)) stop("Step 3 failed: plant activities output not found.")
log(paste("[Step 3] Complete. Output written to:", plant_acts_file))

# ------------------------- Completion -------------------------
log("[Pipeline] Success. Outputs written to:")
log(paste("  - SMILES activities:", resolved_cmp_file))
log(paste("  - Plants for CMPs:", plants_file))
log(paste("  - Activities for CMPs:", cmp_acts_file))
log(paste("  - Activities for Plants:", plant_acts_file))
