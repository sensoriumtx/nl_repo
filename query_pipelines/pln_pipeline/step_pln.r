#!/usr/bin/env Rscript

# ------------------ Load Libraries and Sources ------------------
suppressMessages(library(tidyverse))
suppressMessages(library(parallel))
suppressMessages(source("scripts/SPARQL.R"))
suppressMessages(source("scripts/common.r"))

options(scipen = 100, digits = 4)

# ------------------ Argument Parsing ------------------
args <- commandArgs(TRUE)

endpoint <- "dev"
ids <- NULL
cmp_outdir <- NULL

while (length(args) > 0) {
  if (args[1] == "--endpoint") {
    endpoint <- if (args[2] == "dev") endpoint_dev else args[2]
    args <- args[-c(1, 2)]
  } else if (args[1] == "--plants") {
    ids <- args[2]
    args <- args[-c(1, 2)]
  } else if (args[1] == "--outdir") {
    cmp_outdir <- args[2]
    args <- args[-c(1, 2)]
  } else {
    stop(paste("Unknown argument:", args[1]))
  }
}

# ------------------ Output Directory Handling ------------------
if (is.null(ids) || is.null(cmp_outdir)) {
  stop("Both --plants and --outdir must be provided.")
}

if (!dir.exists(cmp_outdir)) {
  dir.create(cmp_outdir, recursive = TRUE, showWarnings = FALSE)
}
message("✓ Output directory created or confirmed at: ", cmp_outdir)

# ------------------ Step 1: Pull Compounds for Plants ------------------
message("[Step 1] Pulling compounds associated with plants")

plant_input <- str_split(ids, "\\|")[[1]] %>% unique()
pln_ids <- paste(plant_input, collapse = "|")
if (is.null(pln_ids) || pln_ids == "") stop("No valid PLN identifiers found")

cmp_outfile <- file.path(cmp_outdir, "step1_cmp_for_pln.csv")

cmp_cmd <- paste(
  "Rscript scripts/pull_cmp_for_pln.r",
  "--endpoint", endpoint,
  "--plants", shQuote(pln_ids),
  "--out", shQuote(cmp_outfile)
)
system(cmp_cmd)
message("✓ Step 1 completed. Compounds saved to:", cmp_outfile)

# ------------------ Step 2: Pull Activities for Plants ------------------
message("[Step 2] Pulling activities associated with plants")

acts_outfile <- file.path(cmp_outdir, "step2_acts_for_pln.csv")

acts_cmd <- paste(
  "Rscript scripts/pull_acts_for_specific_pln.r",
  "--endpoint", endpoint,
  "--plants", shQuote(pln_ids),
  "--out", shQuote(acts_outfile)
)
system(acts_cmd)
message("✓ Step 2 completed. Activities for plants saved to:", acts_outfile)

# ------------------ Step 3: Pull Activities for Compounds ------------------
message("[Step 3] Pulling activities associated with compounds")

step3_outfile <- file.path(cmp_outdir, "step3_acts_for_cmp.csv")

# Read Step 1 output and extract unique cmp IDs
cmp_df <- tryCatch(read_csv(cmp_outfile, show_col_types = FALSE), error = function(e) NULL)
if (is.null(cmp_df) || !"cmp" %in% colnames(cmp_df)) {
  stop("Step 1 output missing or 'cmp' column not found. Cannot proceed to Step 3.")
}

# Extract unique cmp values and collapse into a single string
cmp_ids_vector <- cmp_df$cmp %>% unique() %>% sort()
cmp_ids_string <- paste(cmp_ids_vector, collapse = "|")
if (cmp_ids_string == "") stop("No valid compound identifiers found for Step 3.")

# Write collapsed string to a file
cmp_ids_file <- file.path(cmp_outdir, "step3_cmp_ids.txt")
writeLines(cmp_ids_string, cmp_ids_file)

# Construct command using file-based compound input
cmp_acts_cmd <- paste(
  "Rscript scripts/pull_acts_for_specific_cmp_ids.r",
  "--endpoint", endpoint,
  "--compound_file", shQuote(cmp_ids_file),
  "--out", shQuote(step3_outfile)
)
system(cmp_acts_cmd)

message("✓ Step 3 completed. Activities for compounds saved to:", step3_outfile)
