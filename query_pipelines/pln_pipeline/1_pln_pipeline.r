#!/usr/bin/env Rscript

suppressMessages(library(tidyverse))
suppressMessages(library(parallel))
suppressMessages(source("scripts/SPARQL.R"))
suppressMessages(source("scripts/common.r"))

options(scipen = 100, digits = 4)
args <- commandArgs(TRUE)

# ------------------ Argument Defaults ------------------
endpoint <- "dev"
ids <- NULL
outdir <- NULL
loop <- TRUE

# ------------------ Argument Parsing ------------------
while (loop) {
  if (args[1] == "--endpoint") {
    if (args[2] == "dev") {
      endpoint <- endpoint_dev
    }
  } else if (args[1] == "--plants") {
    ids <- args[2]
  } else if (args[1] == "--outdir") {
    outdir <- args[2]
  }

  if (length(args) > 1) {
    args <- args[2:length(args)]
  } else {
    loop <- FALSE
  }
}

if (is.null(ids) || is.null(outdir)) {
  stop("Both --plants and --outdir must be provided.")
}

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ------------------ Step 1: Accept plant URIs ------------------
message("[Step 1] Using provided plant URIs (no resolution)")

plant_input <- str_split(ids, "\\|")[[1]] %>% unique()
pln_ids <- plant_input %>% paste(collapse = "|")
if (is.null(pln_ids) || pln_ids == "") stop("No valid PLN identifiers found")

# ------------------ Step 2: Pull Compounds ------------------
message("[Step 2] Pulling compounds associated with plants")

cmp_cmd <- paste(
  "Rscript scripts/pull_cmp_for_pln.r",
  "--endpoint", endpoint,
  "--plants", shQuote(pln_ids),
  "--outdir", shQuote(outdir)
)
system(cmp_cmd)

cmp_outfile <- file.path(outdir, "step2_cmp_for_plants.csv")
if (!file.exists(cmp_outfile)) stop("Step 2 failed: compound file not found")

# ------------------ Step 3: Pull Activities for Plants ------------------
message("[Step 3] Pulling activities associated with plants")

acts_cmd <- paste(
  "Rscript scripts/pull_act_for_plant_enrich.r",
  "--endpoint", endpoint,
  "--plants", shQuote(pln_ids),
  "--out", shQuote(file.path(outdir, "step3_acts_for_pln.csv"))
)
system(acts_cmd)

# ------------------ Step 4: Pull Activities for Compounds ------------------
message("[Step 4] Pulling activities associated with compounds")

cmp_df <- read_csv(cmp_outfile, show_col_types = FALSE)
if (!"cmp" %in% names(cmp_df)) stop("Compound file missing 'cmp' column")

cmp_ids <- cmp_df$cmp %>% unique() %>% paste(collapse = "|")

cmp_acts_cmd <- paste(
  "Rscript scripts/pull_acts_for_cmp_id.r",
  "--endpoint", endpoint,
  "--compounds", shQuote(cmp_ids),
  "--out", shQuote(file.path(outdir, "step4_acts_for_cmp.csv"))
)
system(cmp_acts_cmd)

message("[Pipeline] Execution complete. All outputs saved in:", outdir)
