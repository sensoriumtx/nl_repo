#!/usr/bin/env Rscript

suppressMessages(library(tidyverse))
suppressMessages(library(lubridate))
suppressMessages(library(parallel))

cat("Script starting...\n", flush = TRUE)

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
    if (arg == "--plants") {
      out$plants <- val
    } else if (arg == "--in_file") {
      out$in_file <- val
    } else if (arg == "--endpoint") {
      out$endpoint <- val
    } else if (arg == "--outdir") {
      out$outdir <- val
    } else {
      stop(paste("Unknown argument:", arg))
    }
    i <- i + 2
  }
  return(out)
}

params <- parseArgs(args)
if (is.null(params$outdir)) stop("--outdir must be specified")
dir.create(params$outdir, showWarnings = FALSE, recursive = TRUE)
log_file <- file.path(params$outdir, "pipeline_log.txt")
cat("Pipeline Log\n==============\n", file = log_file)

# ------------------------- Step 1: Resolve or Use PLN -------------------------
log("[Step 1] Preparing plant identifiers")

# Step 1a - Get user input list from --plants or --in_file
plant_input <- NULL
if (!is.null(params$plants)) {
  plant_input <- str_split(params$plants, "\\|")[[1]] %>% unique()
} else if (!is.null(params$in_file)) {
  df <- read_csv(params$in_file, show_col_types = FALSE)
  cols <- tolower(names(df))
  key_col <- intersect(cols, c("plant", "plants", "pln", "pln_label"))[1]
  if (is.null(key_col)) stop("No usable plant column found in input file")
  plant_input <- df[[key_col]] %>% unique() %>% na.omit()
} else {
  stop("You must provide either --plants or --in_file")
}

# Step 1b - If already in 'sen:SENPLN...' form, skip resolution
resolved_pln_file <- file.path(params$outdir, "step1_resolved_pln.csv")
if (all(grepl("^sen:SENPLN[0-9]+$", plant_input))) {
  log("[Step 1] Skipping resolution; PLN IDs provided")
  tibble(pln = plant_input) %>% write_csv(resolved_pln_file)
} else {
  log("[Step 1] Resolving plant names to `pln` URIs via resolve_pln_id.r")
  plant_str <- paste(plant_input, collapse = "|")
  resolve_cmd <- paste(
    "Rscript scripts/resolve_pln_id.r",
    "--endpoint", params$endpoint,
    "--plants", shQuote(plant_str),
    "--out", shQuote(resolved_pln_file)
  )
  system(resolve_cmd)
  if (!file.exists(resolved_pln_file)) stop("Step 1 failed: resolved plant file not found")
}
log("[Step 1] Complete")

pln_ids <- read_csv(resolved_pln_file, show_col_types = FALSE)$pln %>%
  unique() %>% paste(collapse = "|")
if (is.null(pln_ids) || pln_ids == "") stop("No valid `pln` identifiers found")

# ------------------------- Step 2: Pull Compounds for Plants -------------------------
log("[Step 2] Pulling compounds associated with plants")

plant_cmp_file <- file.path(params$outdir, "step2_cmp_for_plants.csv")
cmp_cmd <- paste(
  "Rscript scripts/pull_cmp_for_pln.r",
  "--endpoint", params$endpoint,
  "--plants", shQuote(pln_ids),
  "--out", shQuote(plant_cmp_file)
)
system(cmp_cmd)
if (!file.exists(plant_cmp_file)) stop("Step 2 failed: compound file not found")
log("[Step 2] Complete")

# ------------------------- Step 3: Pull Activities for Plants -------------------------
log("[Step 3] Pulling activities associated with plants")

plant_acts_file <- file.path(params$outdir, "step3_acts_for_pln.csv")
acts_cmd <- paste(
  "Rscript scripts/pull_acts_for_pln_id.r",
  "--endpoint", params$endpoint,
  "--plants", shQuote(pln_ids),
  "--out", shQuote(plant_acts_file)
)
system(acts_cmd)
if (!file.exists(plant_acts_file)) stop("Step 3 failed: plant activity file not found")
log("[Step 3] Complete")

log("[Pipeline] Success. Outputs written to:")
log(paste("  - Resolved Plants:", resolved_pln_file))
log(paste("  - Compounds for Plants:", plant_cmp_file))
log(paste("  - Activities for Plants:", plant_acts_file))
