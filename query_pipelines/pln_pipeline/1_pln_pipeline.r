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

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

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

# ------------------ Step 3: Pull Activities ------------------
message("[Step 3] Pulling activities associated with plants")

act_cmd <- paste(
  "Rscript scripts/pull_acts_for_plant_enrich.r",
  "--endpoint", endpoint,
  "--plants", shQuote(pln_ids),
  "--out", file.path(outdir, "step3_acts_for_pln.csv")
)
system(act_cmd)

message("[Pipeline] Execution complete.")
