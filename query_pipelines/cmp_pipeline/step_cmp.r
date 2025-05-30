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
scoring_file <- NULL
cmp_string <- NULL
cmp_file <- NULL
outdir <- NULL
chunk_size <- 50

while (length(args) > 0) {
  if (args[1] == "--endpoint") {
    endpoint <- if (args[2] == "dev") endpoint_dev else args[2]
    args <- args[-c(1, 2)]
  } else if (args[1] == "--cmp") {
    cmp_string <- args[2]
    args <- args[-c(1, 2)]
  } else if (args[1] == "--cmp_in_file") {
    cmp_file <- args[2]
    args <- args[-c(1, 2)]
  } else if (args[1] == "--scoring") {
    scoring_file <- args[2]
    args <- args[-c(1, 2)]
  } else if (args[1] == "--outdir") {
    outdir <- args[2]
    args <- args[-c(1, 2)]
  } else if (args[1] == "--chunk") {
    chunk_size <- as.integer(args[2])
    args <- args[-c(1, 2)]
  } else {
    stop(paste("Unknown argument:", args[1]))
  }
}

if (is.null(outdir) || is.null(scoring_file)) stop("--scoring and --outdir are required")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
message("✓ Output directory confirmed: ", outdir)

# ------------------ Step 0: Resolve cmp list ------------------
message("[Step 0] Resolving input compound identifiers")
cmp_list <- character()

if (!is.null(cmp_string)) {
  message("→ Using cmp values from --cmp string")
  cmp_list <- unique(str_split(cmp_string, "\\|")[[1]])
} else if (!is.null(cmp_file)) {
  message("→ Reading cmp values from file: ", cmp_file)
  cmp_df <- read_csv(cmp_file, show_col_types = FALSE)
  if (!"cmp" %in% colnames(cmp_df)) stop("Input file must contain 'cmp' column")
  cmp_list <- unique(cmp_df$cmp)
} else {
  message("→ Searching master scoring file using partial matches")
  scoring_df <- read_csv(scoring_file, show_col_types = FALSE)
  message("Available sample rows:")
  print(head(scoring_df, 3))
  search_term <- readline(prompt = "Enter search string for cmp filtering: ")
  matches <- scoring_df %>% filter(if_any(everything(), ~ str_detect(as.character(.), fixed(search_term, ignore_case = TRUE))))
  if (nrow(matches) == 0) stop("No matches found in scoring file.")
  cmp_list <- unique(matches$cmp)
  write_csv(matches, file.path(outdir, "matched_cmp_rows.csv"))
}

cmp_list <- unique(cmp_list[!is.na(cmp_list)])
message("✓ Final resolved cmp count: ", length(cmp_list))

# ------------------ Step 1: Pull Plants for Compounds ------------------
message("[Step 1] Pulling plants for compounds (chunked)")
step1_dir <- file.path(outdir, "step1_plants_chunks")
if (!dir.exists(step1_dir)) dir.create(step1_dir, recursive = TRUE)

chunks <- split(cmp_list, ceiling(seq_along(cmp_list)/chunk_size))
chunk_files <- c()

for (i in seq_along(chunks)) {
  chunk <- chunks[[i]]
  chunk_file <- file.path(step1_dir, sprintf("step1_plants_chunk_%03d.csv", i))
  write_csv(tibble(cmp = chunk), chunk_file)
  chunk_files <- c(chunk_files, chunk_file)

  cmd <- paste(
    "Rscript scripts/step1_plants_chunked.r",
    shQuote(chunk_file),
    shQuote(chunk_file)  # Reusing for output inside script
  )
  message("Running step1_plants chunk ", i, " → ", chunk_file)
  system(cmd)
}

# Merge all chunked outputs
step1_merged <- map_dfr(chunk_files, ~read_csv(.x, show_col_types = FALSE))
write_csv(step1_merged, file.path(outdir, "step1_plants_complete.csv"))
message("✓ Step 1 complete → All plant data merged.")

# ------------------ Future Steps Go Here ------------------
# Example placeholder
message("→ Next steps (e.g., activity pulling) can now proceed using merged compound data")
