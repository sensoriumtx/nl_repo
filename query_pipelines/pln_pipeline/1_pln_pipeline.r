#!/usr/bin/env Rscript

suppressMessages(library(tidyverse))
suppressMessages(library(parallel))
suppressMessages(source("scripts/SPARQL.R"))
suppressMessages(source("scripts/common.r"))

options(scipen = 100, digits = 4)
args <- commandArgs(TRUE)

# ------------------ Argument Defaults ------------------
ids <- NULL
endpoint <- NULL
outdir <- NULL
loop <- TRUE

# ------------------ Argument Parsing ------------------
while (loop) {
  if (args[1] == "--endpoint") {
    if (args[2] == "dev") {
      endpoint <- endpoint_dev
    }
  }

  if (args[1] == "--plants") {
    ids <- args[2]
  }

  if (args[1] == "--outdir") {
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

# Construct output path
outFile <- file.path(outdir, "step2_cmp_for_plants.csv")
dir.create(dirname(outFile), recursive = TRUE, showWarnings = FALSE)

# ------------------ SPARQL: Get PLN Metadata ------------------
log("[Step 1] Using provided plant URIs (no resolution)")

plant_input <- NULL
if (!is.null(params$plants)) {
  plant_input <- str_split(params$plants, "\\|")[[1]] %>% unique()
} else if (!is.null(params$in_file)) {
  df <- read_csv(params$in_file, show_col_types = FALSE)
  cols <- tolower(names(df))
  key_col <- intersect(cols, c("pln", "plant", "plants"))[1]
  if (is.null(key_col)) stop("No usable column found in input file")
  plant_input <- df[[key_col]] %>% unique() %>% na.omit()
} else {
  stop("You must provide either --plants or --in_file")
}

pln_ids <- plant_input %>% unique() %>% paste(collapse = "|")
if (is.null(pln_ids) || pln_ids == "") stop("No valid PLN identifiers found")



# ------------------ SPARQL: Get Compounds for Plants ------------------
log("[Step 2] Pulling compounds associated with plants")

cmp_cmd <- paste(
  "Rscript scripts/pull_cmp_for_pln.r",
  "--endpoint", params$endpoint,
  "--plants", shQuote(pln_ids),
  "--outdir", shQuote(params$outdir)
)
system(cmp_cmd)


# ------------------ Safe Join ------------------
if ("pln" %in% names(combined) && "pln" %in% names(df.id)) {
  message("Joining plant metadata with compound results")
  combined <- combined %>%
    left_join(df.id, by = "pln")
} else {
  warning("Join skipped: 'pln' column missing in one or both data frames.")
}

# ------------------ Save Output ------------------
if (nrow(combined) == 0) {
  message("No compound mappings found.")
} else {
  message(paste("Total rows retrieved:", nrow(combined)))
}
write_csv(combined, outFile)
