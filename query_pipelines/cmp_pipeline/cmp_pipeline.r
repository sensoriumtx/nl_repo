# Nick Laskowski
# Version 1.1
# CMP-Based Semantic Association Pipeline 

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
    else if (arg == "--endpoint") out$endpoint <- val
    else if (arg == "--outdir") out$outdir <- val
    else stop(paste("Unknown argument:", arg))
    i <- i + 2
  }
  return(out)
}

params <- parseArgs(args)
if (is.null(params$outdir)) stop("--outdir must be specified")
if (is.null(params$cmp)) stop("--cmp must be provided")

# ------------------------- Directory Setup -------------------------
dir.create(params$outdir, showWarnings = FALSE, recursive = TRUE)
log_file <- file.path(params$outdir, "pipeline_log.txt")
cat("Pipeline Log\n==============\n", file = log_file)

# ------------------------- Step 0: Pull Activities for CMP -------------------------
log("[Step 0] Pulling activities for CMP")

cmp_acts_file <- file.path(params$outdir, "step0_acts_for_cmp.csv")

pull_cmp_acts_cmd <- paste(
  "Rscript scripts/pull_acts_for_specific_cmp_ids.r",
  "--endpoint", params$endpoint,
  "--compound", shQuote(params$cmp),
  "--out", shQuote(cmp_acts_file)
)

system(pull_cmp_acts_cmd)
if (!file.exists(cmp_acts_file)) stop("Step 0 failed: CMP activities output not found.")
log(paste("[Step 0] Complete. Output written to:", cmp_acts_file))

# ------------------------- Step 1: Pull Plants for CMP -------------------------
log("[Step 1] Pulling plants associated with resolved compounds")

plants_file <- file.path(params$outdir, "step1_plants_for_cmp.csv")
plant_cmd <- paste(
  "Rscript scripts/pull_plant_for_compound_ids.r",
  "--endpoint", params$endpoint,
  "--compound_activity_file", shQuote(cmp_acts_file),
  "--cmp_id_column", "cmp",
  "--out", shQuote(plants_file)
)
system(plant_cmd)
if (!file.exists(plants_file)) stop("Step 1 failed: plant output not found.")
log(paste("[Step 1] Complete. Output written to:", plants_file))

# ------------------------- Step 2: Pull Acts for PLN -------------------------
log("[Step 2] Pulling activities associated with plants")

plant_labels <- read_csv(plants_file, show_col_types = FALSE)$pln_label %>%
  unique() %>% na.omit() %>% paste(collapse = "|")

plant_acts_file <- file.path(params$outdir, "step2_acts_for_pln.csv")
pln_act_cmd <- paste(
  "Rscript scripts/pull_acts_for_specific_pln.r",
  "--endpoint", params$endpoint,
  "--plants", shQuote(plant_labels),
  "--out", shQuote(plant_acts_file)
)
system(pln_act_cmd)
if (!file.exists(plant_acts_file)) stop("Step 2 failed: plant activities output not found.")
log(paste("[Step 2] Complete. Output written to:", plant_acts_file))

# ------------------------- Completion -------------------------
log("[Pipeline] Success. Outputs written to:")
log(paste("  - CMP activities:", cmp_acts_file))
log(paste("  - Plants for CMPs:", plants_file))
log(paste("  - Activities for Plants:", plant_acts_file))
