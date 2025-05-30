#!/usr/bin/env Rscript

# ------------------ Load Libraries and Sources ------------------
suppressMessages(library(tidyverse))
suppressMessages(library(parallel))
suppressMessages(source("scripts/SPARQL.R"))
suppressMessages(source("scripts/common.r"))

options(scipen = 100, digits = 4)

# ------------------ Chunked Step Function ------------------
run_chunked_step <- function(
  step_name,
  input_ids,
  step_script,
  arg_flag,
  chunk_size = 50,
  endpoint = NULL,
  output_dir
) {
  message("[", toupper(step_name), "] Starting chunked execution...")

  chunk_dir <- file.path(output_dir, paste0(step_name, "_chunks"))
  dir.create(chunk_dir, recursive = TRUE, showWarnings = FALSE)

  input_chunks <- split(input_ids, ceiling(seq_along(input_ids) / chunk_size))
  chunk_files <- list()

  for (i in seq_along(input_chunks)) {
    chunk_vals <- input_chunks[[i]]
    chunk_str <- paste(chunk_vals, collapse = "|")
    chunk_out <- file.path(chunk_dir, sprintf("%s_chunk_%03d.csv", step_name, i))

    cmd <- paste(
      "Rscript", step_script,
      if (!is.null(endpoint)) paste("--endpoint", endpoint) else NULL,
      arg_flag, shQuote(chunk_str),
      "--out", shQuote(chunk_out)
    )

    message("Running chunk ", i, " (", length(chunk_vals), " IDs): ", chunk_out)
    system(cmd)
    chunk_files[[i]] <- chunk_out
  }

  # Merge chunk results
  final_outfile <- file.path(output_dir, paste0(step_name, "_final.csv"))
  final_df <- bind_rows(lapply(chunk_files, read_csv, show_col_types = FALSE))
  write_csv(final_df, final_outfile)
  message("[", toupper(step_name), "] Completed. Merged file saved to: ", final_outfile)

  return(final_outfile)
}

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

if (is.null(ids) || is.null(cmp_outdir)) {
  stop("Both --plants and --outdir must be provided.")
}

# ------------------ Output Directory Setup ------------------
if (!dir.exists(cmp_outdir)) {
  dir.create(cmp_outdir, recursive = TRUE, showWarnings = FALSE)
}
message("✓ Output directory created or confirmed at: ", cmp_outdir)

# ------------------ Step 1: Pull Compounds for Plants ------------------
message("[Step 1] Pulling compounds associated with plants")

plant_input <- str_split(ids, "\\|")[[1]] %>% unique()
if (length(plant_input) == 0) stop("No valid plant identifiers provided.")

step1_outfile <- run_chunked_step(
  step_name = "step1",
  input_ids = plant_input,
  step_script = "scripts/pull_cmp_for_pln.r",
  arg_flag = "--plants",
  chunk_size = 300,
  endpoint = endpoint,
  output_dir = cmp_outdir
)

# ------------------ Step 2: Pull Activities for Plants ------------------
message("[Step 2] Pulling activities associated with plants")

step2_outfile <- run_chunked_step(
  step_name = "step2",
  input_ids = plant_input,
  step_script = "scripts/pull_acts_for_specific_pln.r",
  arg_flag = "--plants",
  chunk_size = 300,
  endpoint = endpoint,
  output_dir = cmp_outdir
)

# ------------------ Step 3: Pull Activities for Compounds ------------------
message("[Step 3] Pulling activities associated with compounds")

cmp_df <- tryCatch(read_csv(step1_outfile, show_col_types = FALSE), error = function(e) NULL)
if (is.null(cmp_df) || !"cmp" %in% colnames(cmp_df)) {
  stop("Step 1 output missing or 'cmp' column not found. Cannot proceed to Step 3.")
}

cmp_ids <- cmp_df$cmp %>% unique() %>% sort()
if (length(cmp_ids) == 0) stop("No valid compound identifiers found for Step 3.")

step3_outfile <- run_chunked_step(
  step_name = "step3",
  input_ids = cmp_ids,
  step_script = "scripts/pull_acts_for_specific_cmp_ids.r",
  arg_flag = "--compound",
  chunk_size = 300,
  endpoint = endpoint,
  output_dir = cmp_outdir
)

message("All steps completed successfully.")
