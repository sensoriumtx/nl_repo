# Nick Laskowski
# Version 1.0
# Wildcard PLN Search Pipeline

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
    if (arg == "--search") out$search <- val
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

# ------------------------- Step 0: Search in PLN Master -------------------------
log("[Step 0] Performing wildcard search across all string fields in PLN master file")

if (is.null(params$in_file)) stop("You must provide --in_file to search from a master CSV.")
if (is.null(params$search)) stop("You must provide a search string using --search")

search_term <- tolower(params$search) %>% str_trim()
log(paste("Search string (normalized):", search_term))

master_df <- read_csv(params$in_file, show_col_types = FALSE, guess_max = 10000)

# Identify all character columns
char_cols <- master_df %>%
  select(where(is.character)) %>%
  names()

# Filter rows where any string column contains the search term
matches <- master_df %>%
  filter(if_any(all_of(char_cols), ~ str_detect(tolower(.), fixed(search_term, ignore_case = TRUE))))

if (nrow(matches) == 0) stop(paste("No matches found for search term:", search_term))

resolved_pln_file <- file.path(params$outdir, "step0_wildcard_matched_rows.csv")
write_csv(matches, resolved_pln_file)

# Ensure required column
if (!"pln_label" %in% colnames(matches)) stop("The column 'pln_label' must be present in matched data.")
pln_labels <- matches$pln_label %>% unique() %>% na.omit() %>% paste(collapse = "|")

log(paste("[Step 0] Match complete. Found", nrow(matches), "rows matching:", search_term))
log(paste("Unique PLN labels:", str_count(pln_labels, "\\|") + 1))
log(paste("Matched rows saved to:", resolved_pln_file))

# ------------------------- Step 1: Pull Compounds for PLN -------------------------
log("[Step 1] Pulling compounds associated with matched plants")

compounds_file <- file.path(params$outdir, "step1_compounds_for_pln.csv")
cmp_cmd <- paste(
  "Rscript scripts/pull_cmp_for_pln.r",
  "--endpoint", params$endpoint,
  "--plants", shQuote(pln_labels),
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
log("[Step 3] Pulling activities associated with matched plants")

plant_acts_file <- file.path(params$outdir, "step3_acts_for_pln.csv")
pln_act_cmd <- paste(
  "Rscript scripts/pull_acts_for_specific_pln.r",
  "--endpoint", params$endpoint,
  "--plants", shQuote(pln_labels),
  "--out", shQuote(plant_acts_file)
)
system(pln_act_cmd)
if (!file.exists(plant_acts_file)) stop("Step 3 failed: plant activity output not found.")
log("[Step 3] Complete")

# ------------------------- Completion -------------------------
log("[Pipeline] Success. Outputs written to:")
log(paste("  - Matched PLN rows:", resolved_pln_file))
log(paste("  - Compounds for PLNs:", compounds_file))
log(paste("  - Activities for CMPs:", cmp_acts_file))
log(paste("  - Activities for PLNs:", plant_acts_file))
