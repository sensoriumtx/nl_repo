# Nick Laskowski
# Version 1.1
# Wildcard Search to Semantic Association (PLN version) — Parallel Multi-Search and Fault Tolerance

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
if (is.null(params$search)) stop("--search must be provided (e.g., 'Panax ginseng|Withania somnifera')")
if (is.null(params$outdir)) stop("--outdir must be specified")

dir.create(params$outdir, showWarnings = FALSE, recursive = TRUE)
log_file <- file.path(params$outdir, "pipeline_log.txt")
cat("Pipeline Log\n==============\n", file = log_file)

# ------------------------- Step 0: Search PLNs -------------------------
log("[Step 0] Performing parallel wildcard search across plant name fields")

if (is.null(params$in_file)) stop("You must provide --in_file to search from a master CSV.")
search_terms <- str_split(params$search, "\\|")[[1]] %>% tolower() %>% str_trim()
log(paste("Search terms:", paste(search_terms, collapse = " | ")))

master_df <- read_csv(params$in_file, show_col_types = FALSE, guess_max = 10000)

char_cols <- master_df %>%
  select(where(is.character)) %>%
  names()

match_term <- function(term) {
  master_df %>%
    filter(if_any(all_of(char_cols), ~ str_detect(tolower(.), fixed(term, ignore_case = TRUE))))
}

matched_list <- mclapply(search_terms, match_term, mc.cores = min(length(search_terms), detectCores()))
all_matches <- bind_rows(matched_list) %>% distinct()

resolved_pln_file <- file.path(params$outdir, "step0_wildcard_matched_plants.csv")
if (nrow(all_matches) == 0) {
  log("[Step 0] No matching rows found. Writing blank file.")
  write_csv(tibble(pln_label = character()), resolved_pln_file)
  cmp_ids <- ""
} else {
  write_csv(all_matches, resolved_pln_file)
  if (!"pln_label" %in% colnames(all_matches)) stop("The column 'pln_label' must be present.")
  log(paste("[Step 0] Match complete. Found", nrow(all_matches), "matched rows"))
}

# ------------------------- Step 1: Pull Compounds for PLN -------------------------
log("[Step 1] Pulling compounds associated with matched plants")

compounds_file <- file.path(params$outdir, "step1_compounds_for_pln.csv")
step1_has_data <- FALSE

if (nrow(all_matches) > 0) {
  plants_joined <- all_matches$pln_label %>% unique() %>% na.omit() %>% paste(collapse = "|")
  cmp_cmd <- paste(
    "Rscript scripts/pull_cmp_for_pln.r",
    "--endpoint", params$endpoint,
    "--plants", shQuote(plants_joined),
    "--out", shQuote(compounds_file)
  )
  system(cmp_cmd)
  if (file.exists(compounds_file)) {
    cmp_df <- read_csv(compounds_file, show_col_types = FALSE)
    if (nrow(cmp_df) > 0) step1_has_data <- TRUE
    else {
      log("[Step 1] Compound file is empty. Writing blank.")
      write_csv(tibble(), compounds_file)
    }
  } else {
    log("[Step 1] Compound file not found. Writing blank.")
    write_csv(tibble(), compounds_file)
  }
} else {
  write_csv(tibble(), compounds_file)
  log("[Step 1] Skipped due to no plant matches")
}

log("[Step 1] Complete")

# ------------------------- Step 2: Pull Activities for CMP -------------------------
log("[Step 2] Pulling activities associated with compounds")

cmp_acts_file <- file.path(params$outdir, "step2_acts_for_cmp.csv")
step2_has_data <- FALSE

if (step1_has_data) {
  cmp_ids <- cmp_df$cmp %>% unique() %>% na.omit() %>% sort() %>% paste(collapse = "|")
  if (cmp_ids != "") {
    act_cmd <- paste(
      "Rscript scripts/pull_acts_for_specific_cmp_ids.r",
      "--endpoint", params$endpoint,
      "--compound", shQuote(cmp_ids),
      "--out", shQuote(cmp_acts_file)
    )
    system(act_cmd)

    if (file.exists(cmp_acts_file)) {
      cmp_acts <- read_csv(cmp_acts_file, show_col_types = FALSE)
      if (nrow(cmp_acts) > 0) step2_has_data <- TRUE
      else {
        log("[Step 2] CMP activity file is empty. Writing blank.")
        write_csv(tibble(), cmp_acts_file)
      }
    } else {
      log("[Step 2] CMP activity file not found. Writing blank.")
      write_csv(tibble(), cmp_acts_file)
    }
  } else {
    write_csv(tibble(), cmp_acts_file)
    log("[Step 2] No cmp IDs resolved.")
  }
} else {
  write_csv(tibble(), cmp_acts_file)
  log("[Step 2] Skipped due to no compound output from Step 1")
}

log("[Step 2] Complete")

# ------------------------- Step 3: Pull Activities for PLN -------------------------
log("[Step 3] Pulling activities directly associated with matched plants")

plant_acts_file <- file.path(params$outdir, "step3_acts_for_pln.csv")

if (nrow(all_matches) > 0) {
  plants_joined <- all_matches$pln_label %>% unique() %>% na.omit() %>% paste(collapse = "|")
  pln_act_cmd <- paste(
    "Rscript scripts/pull_acts_for_specific_pln.r",
    "--endpoint", params$endpoint,
    "--plants", shQuote(plants_joined),
    "--out", shQuote(plant_acts_file)
  )
  system(pln_act_cmd)

  if (!file.exists(plant_acts_file)) {
    log("[Step 3] Plant activity file not found. Writing blank.")
    write_csv(tibble(), plant_acts_file)
  } else {
    pln_acts <- read_csv(plant_acts_file, show_col_types = FALSE)
    if (nrow(pln_acts) == 0) {
      log("[Step 3] Plant activity file is empty. Writing blank.")
      write_csv(tibble(), plant_acts_file)
    }
  }
} else {
  write_csv(tibble(), plant_acts_file)
  log("[Step 3] Skipped due to no plant matches")
}

log("[Step 3] Complete")

# ------------------------- Completion -------------------------
log("[Pipeline] Success. Outputs written to:")
log(paste("  - Matched PLN rows:", resolved_pln_file))
log(paste("  - Compounds for PLNs:", compounds_file))
log(paste("  - Activities for CMPs:", cmp_acts_file))
log(paste("  - Activities for PLNs:", plant_acts_file))
