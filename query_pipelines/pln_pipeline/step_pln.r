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
cmp_outfile <- NULL
acts_outfile <- NULL

while (length(args) > 0) {
  if (args[1] == "--endpoint") {
    endpoint <- if (args[2] == "dev") endpoint_dev else args[2]
    args <- args[-c(1, 2)]
  } else if (args[1] == "--plants") {
    ids <- args[2]
    args <- args[-c(1, 2)]
  } else if (args[1] == "--out") {
    cmp_outfile <- args[2]
    args <- args[-c(1, 2)]
  } else {
    stop(paste("Unknown argument:", args[1]))
  }
}

if (is.null(ids) || is.null(cmp_outfile)) {
  stop("Both --plants and --out must be provided.")
}

dir.create(dirname(cmp_outfile), recursive = TRUE, showWarnings = FALSE)

# ------------------ Step 1: Pull Compounds for Plants ------------------
message("[Step 1] Pulling compounds associated with plants")

plant_input <- str_split(ids, "\\|")[[1]] %>% unique()
pln_ids <- paste(plant_input, collapse = "|")
if (is.null(pln_ids) || pln_ids == "") stop("No valid PLN identifiers found")

cmp_cmd <- paste(
  "Rscript scripts/pull_cmp_for_pln.r",
  "--endpoint", endpoint,
  "--plants", shQuote(pln_ids),
  "--out", shQuote(cmp_outfile)
)
system(cmp_cmd)

# ------------------ Step 2: Pull Activities for Plants ------------------
message("[Step 2] Pulling activities associated with plants")

acts_outfile <- file.path(dirname(cmp_outfile), "step2_acts_for_pln.csv")

acts_cmd <- paste(
  "Rscript scripts/pull_act_for_plant_enrich.r",
  "--endpoint", endpoint,
  "--plants", shQuote(pln_ids),
  "--out", shQuote(acts_outfile)
)
system(acts_cmd)

message("✓ Step 2 completed. Activities for plants saved to:", acts_outfile)
