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

# ------------------------- Step 0: Open Grep-style CMP Search -------------------------
log("[Step 0] Performing grep-like search across all string fields in master file")

if (is.null(params$in_file)) stop("You must provide --in_file to search from a master CSV.")
if (is.null(params$search)) stop("You must provide a search string using --search")

search_term <- tolower(params$search)
master_df <- read_csv(params$in_file, show_col_types = FALSE)

# Identify all character/text columns
char_cols <- master_df %>%
  select(where(is.character)) %>%
  names()

# Filter rows where any string column contains the search term
matches <- master_df %>%
  filter(if_any(all_of(char_cols), ~ str_detect(tolower(.), fixed(search_term, ignore_case = TRUE))))

if (nrow(matches) == 0) stop("No matches found for search term:", search_term)

# Output full matched rows and extract unique cmp
resolved_cmp_file <- file.path(params$outdir, "step0_grep_matched_rows.csv")
write_csv(matches, resolved_cmp_file)

cmp_ids <- matches$cmp %>% unique() %>% na.omit() %>% sort() %>% paste(collapse = "|")
if (cmp_ids == "") stop("No valid cmp IDs found in matched rows.")

log(paste("[Step 0] Match complete. Found", nrow(matches), "rows matching:", search_term))
log(paste("Unique cmp IDs:", str_count(cmp_ids, "sen:")))
log(paste("Matched rows saved to:", resolved_cmp_file))

# ------------------------- Step 1: Pull Plants for CMP -------------------------
log("[Step 1] Pulling plants associated with resolved compounds")

plants_file <- file.path(params$outdir, "step1_plants_for_cmp.csv")
plant_cmd <- paste(
  "Rscript scripts/pull_plant_for_cmp_id.r",
  "--endpoint", params$endpoint,
  "--compound_activity_file", shQuote(resolved_cmp_file),
  "--cmp_id_column", "cmp",
  "--out", shQuote(plants_file)
)
system(plant_cmd)
if (!file.exists(plants_file)) stop("Step 1 failed: plant output not found.")
log("[Step 1] Complete")

# ------------------------- Step 2: Pull Acts for CMP -------------------------
log("[Step 2] Pulling activities associated with compounds")

cmp_acts_file <- file.path(params$outdir, "step2_acts_for_cmp.csv")
act_cmd <- paste(
  "Rscript scripts/pull_acts_for_cmp_id.r",
  "--endpoint", params$endpoint,
  "--compound", shQuote(cmp_ids),
  "--out", shQuote(cmp_acts_file)
)
system(act_cmd)
if (!file.exists(cmp_acts_file)) stop("Step 2 failed: cmp activities output not found.")
log("[Step 2] Complete")

# ------------------------- Step 3: Pull Acts for PLN -------------------------
log("[Step 3] Pulling activities associated with plants")

plant_labels <- read_csv(plants_file, show_col_types = FALSE)$pln_label %>%
  unique() %>% na.omit() %>% paste(collapse = "|")

plant_acts_file <- file.path(params$outdir, "step3_acts_for_pln.csv")
pln_act_cmd <- paste(
  "Rscript scripts/pull_acts_for_pln_id.r",
  "--endpoint", params$endpoint,
  "--plants", shQuote(plant_labels),
  "--out", shQuote(plant_acts_file)
)
system(pln_act_cmd)
if (!file.exists(plant_acts_file)) stop("Step 3 failed: plant activities output not found.")
log("[Step 3] Complete")

# ------------------------- Completion -------------------------
log("[Pipeline] Success. Outputs written to:")
log(paste("  - Matched CMP rows:", resolved_cmp_file))
log(paste("  - Plants for CMPs:", plants_file))
log(paste("  - Activities for CMPs:", cmp_acts_file))
log(paste("  - Activities for Plants:", plant_acts_file))
