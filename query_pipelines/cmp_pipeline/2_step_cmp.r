#!/usr/bin/env Rscript

suppressMessages(library(tidyverse))
suppressMessages(library(lubridate))

cat("Step 1 starting...\n", flush = TRUE)

# ------------------ Parse Args ------------------
args <- commandArgs(trailingOnly = TRUE)

params <- list()
i <- 1
while (i <= length(args)) {
  key <- args[i]
  val <- if (i + 1 <= length(args)) args[i + 1] else NA
  if (is.na(val)) stop(paste("Missing value for", key))
  if (key %in% c("--cmp", "--cmp_file", "--smiles", "--smiles_file", "--endpoint", "--outdir")) {
    params[[substr(key, 3, nchar(key))]] <- val
  } else {
    stop(paste("Unknown argument:", key))
  }
  i <- i + 2
}

# ------------------ Validate ------------------
if (is.null(params$outdir)) stop("--outdir must be specified")
dir.create(params$outdir, recursive = TRUE, showWarnings = FALSE)

if (is.null(params$endpoint)) stop("--endpoint must be specified")

if (is.null(params$cmp) && is.null(params$cmp_file) && is.null(params$smiles) && is.null(params$smiles_file)) {
  stop("You must provide one of: --cmp, --cmp_file, --smiles, or --smiles_file")
}

# ------------------ Determine Input ------------------
if (!is.null(params$cmp)) {
  cmp_vals <- str_split(params$cmp, "\\|")[[1]] %>% unique() %>% na.omit() %>% trimws()
  input_type <- "cmp"
  target_script <- "scripts/pull_acts_for_specific_cmp_ids.r"
  input_arg_flag <- "--cmp"
  input_str <- paste(cmp_vals, collapse = "|")

} else if (!is.null(params$cmp_file)) {
  cmp_df <- read_csv(params$cmp_file, show_col_types = FALSE)
  if (!"cmp" %in% names(cmp_df)) stop("❌ cmp_file must have a column named 'cmp'")
  cmp_vals <- cmp_df$cmp %>% unique() %>% na.omit() %>% trimws()
  input_type <- "cmp"
  target_script <- "scripts/pull_acts_for_specific_cmp_ids.r"
  input_arg_flag <- "--cmp"
  input_str <- paste(cmp_vals, collapse = "|")

} else if (!is.null(params$smiles)) {
  smiles_vals <- str_split(params$smiles, "\\|")[[1]] %>% unique() %>% na.omit() %>% trimws()
  input_type <- "smiles"
  target_script <- "/sensorium-research-kb/dev/data/query_output/testing/for_nick/Nick_dev/R_Script_Dev/pull_act_for_specific_smiles.r"
  input_arg_flag <- "--smiles"
  input_str <- paste(smiles_vals, collapse = "|")

} else if (!is.null(params$smiles_file)) {
  smiles_df <- read_csv(params$smiles_file, show_col_types = FALSE)
  if (!"isoSmiles" %in% names(smiles_df)) stop("❌ smiles_file must have a column named 'isoSmiles'")
  smiles_vals <- smiles_df$isoSmiles %>% unique() %>% na.omit() %>% trimws()
  input_type <- "smiles"
  target_script <- "/sensorium-research-kb/dev/data/query_output/testing/for_nick/Nick_dev/R_Script_Dev/pull_act_for_specific_smiles.r"
  input_arg_flag <- "--smiles"
  input_str <- paste(smiles_vals, collapse = "|")
}


# ------------------ Build & Run Command ------------------
out_file <- file.path(params$outdir, paste0("step1_", input_type, "_acts.csv"))

cmd <- paste(
  "Rscript", target_script,
  "--endpoint", shQuote(params$endpoint),
  input_arg_flag, shQuote(input_str),
  "--out", shQuote(out_file)
)

cat("→ Running command:\n", cmd, "\n")
exit_code <- system(cmd)

if (exit_code != 0 || !file.exists(out_file)) {
  stop("Script failed or output file not generated.")
}

cat("Step 1 complete. Output saved to:\n", out_file, "\n")
