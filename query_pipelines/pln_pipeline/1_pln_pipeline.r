#!/usr/bin/env Rscript

suppressMessages(library(tidyverse))
suppressMessages(library(parallel))
suppressMessages(source("scripts/SPARQL.R"))
suppressMessages(source("scripts/common.r"))

options(scipen = 100, digits = 4)
args <- commandArgs(TRUE)

# ------------------ Argument Parsing ------------------
endpoint <- "dev"
ids <- NULL
outFile <- NULL
filter_out_act <- NULL

while (length(args) > 0) {
  if (args[1] == "--endpoint") {
    endpoint <- if (args[2] == "dev") endpoint_dev else args[2]
    args <- args[-c(1, 2)]
  } else if (args[1] == "--plants") {
    ids <- args[2]
    args <- args[-c(1, 2)]
  } else if (args[1] == "--out") {
    outFile <- args[2]
    args <- args[-c(1, 2)]
  } else if (args[1] == "--filter_out_act") {
    filter_out_act <- args[2]
    args <- args[-c(1, 2)]
  } else {
    stop(paste("Unknown argument:", args[1]))
  }
}

if (is.null(ids) || is.null(outFile)) {
  stop("Both --plants and --out must be provided.")
}

dir.create(dirname(outFile), recursive = TRUE, showWarnings = FALSE)

# ------------------ Step 1 ------------------
message("[Step 1] Using provided plant URIs (no resolution)")
plant_input <- str_split(ids, "\\|")[[1]] %>% unique()
pln_ids <- paste(plant_input, collapse = "|")
if (is.null(pln_ids) || pln_ids == "") stop("No valid PLN identifiers found")

# ------------------ Step 2: Pull Compounds ------------------
message("[Step 2] Pulling compounds associated with plants")
cmp_cmd <- paste(
  "Rscript scripts/pull_cmp_for_pln.r",
  "--endpoint", endpoint,
  "--plants", shQuote(pln_ids),
  "--out", shQuote(outFile)
)
system(cmp_cmd)

if (!file.exists(outFile)) {
  warning("Step 2 failed: compound file not found")
} else {
  message("✓ Compounds successfully retrieved")
}

# ------------------ Step 3: Pull Activities for Plants ------------------
message("[Step 3] Pulling activities associated with plants")
acts_outfile <- file.path(dirname(outFile), "step3_acts_for_pln.csv")
acts_cmd <- paste(
  "Rscript scripts/pull_act_for_plant_enrich.r",
  "--endpoint", endpoint,
  "--plants", shQuote(pln_ids),
  "--out", shQuote(acts_outfile)
)
if (!is.null(filter_out_act)) {
  acts_cmd <- paste(acts_cmd, "--filter_out_act", shQuote(filter_out_act))
}
system(acts_cmd)

if (!file.exists(acts_outfile)) {
  warning("Step 3 failed: plant activity file not found")
} else {
  message("✓ Activities for plants successfully retrieved")
}

# ------------------ Step 4: Pull Activities for Compounds ------------------
message("[Step 4] Pulling activities associated with compounds")
if (file.exists(outFile)) {
  cmp_df <- tryCatch({
    read_csv(outFile, show_col_types = FALSE)
  }, error = function(e) {
    warning("Failed to read compound file: ", e$message)
    return(NULL)
  })

  if (!is.null(cmp_df) && "cmp" %in% names(cmp_df)) {
    cmp_ids <- cmp_df$cmp %>% unique() %>% paste(collapse = "|")
    cmp_acts_outfile <- file.path(dirname(outFile), "step4_acts_for_cmp.csv")

    cmp_acts_cmd <- paste(
      "Rscript scripts/pull_acts_for_cmp_id.r",
      "--endpoint", endpoint,
      "--compounds", shQuote(cmp_ids),
      "--out", shQuote(cmp_acts_outfile)
    )
    system(cmp_acts_cmd)
    if (!file.exists(cmp_acts_outfile)) {
      warning("Step 4 failed: compound activity file not found")
    } else {
      message("✓ Activities for compounds successfully retrieved")
    }
  } else {
    warning("No compound IDs found in compound file.")
  }
} else {
  warning("Compound file missing. Skipping Step 4.")
}

message("[Pipeline] Execution complete. All outputs saved in:", dirname(outFile))
