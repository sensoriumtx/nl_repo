#!/usr/bin/env Rscript

suppressMessages(library(tidyverse))
suppressMessages(library(readr))
suppressMessages(library(stringr))

args <- commandArgs(trailingOnly = TRUE)

params <- list()
i <- 1
while (i <= length(args)) {
  if (args[i] == "--cmp") {
    params$cmp <- args[i + 1]; i <- i + 1
  } else if (args[i] == "--cmp_file") {
    params$cmp_file <- args[i + 1]; i <- i + 1
  } else if (args[i] == "--smiles") {
    params$smiles <- args[i + 1]; i <- i + 1
  } else if (args[i] == "--smiles_file") {
    params$smiles_file <- args[i + 1]; i <- i + 1
  } else if (args[i] == "--endpoint") {
    params$endpoint <- args[i + 1]; i <- i + 1
  } else if (args[i] == "--outdir") {
    params$outdir <- args[i + 1]; i <- i + 1
  }
  i <- i + 1
}

# -------------- Helper to read file and extract column values --------------
read_compound_file_column <- function(file_path, column_name) {
  message(paste0("📄 Reading column '", column_name, "' from: ", file_path))
  if (grepl("\\.tsv$", file_path, ignore.case = TRUE)) {
    df <- read_tsv(file_path, locale = locale(encoding = "UTF-8"), show_col_types = FALSE)
  } else {
    df <- read_csv(file_path, locale = locale(encoding = "UTF-8"), show_col_types = FALSE)
  }
  if (!(column_name %in% colnames(df))) {
    stop(paste0("❌ Column '", column_name, "' not found in ", file_path))
  }
  values <- df[[column_name]] %>% unique() %>% na.omit() %>% trimws() %>% unlist()
  return(values)
}

# ------------------ Determine Input ------------------
input_type <- NULL
input_str <- NULL
target_script <- NULL
input_arg_flag <- NULL

if (!is.null(params$cmp)) {
  cmp_vals <- str_split(params$cmp, "\\|")[[1]] %>% unique() %>% na.omit() %>% trimws()
  input_type <- "cmp"
  target_script <- "scripts/pull_acts_for_specific_cmp_ids.r"
  input_arg_flag <- "--compound"
  input_str <- paste(cmp_vals, collapse = "|")

} else if (!is.null(params$cmp_file)) {
  cmp_vals <- read_compound_file_column(params$cmp_file, "cmp")
  input_type <- "cmp"
  target_script <- "scripts/pull_acts_for_specific_cmp_ids.r"
  input_arg_flag <- "--compound"
  input_str <- paste(cmp_vals, collapse = "|")

} else if (!is.null(params$smiles)) {
  smiles_vals <- str_split(params$smiles, "\\|")[[1]] %>% unique() %>% na.omit() %>% trimws()
  input_type <- "smiles"
  target_script <- "/sensorium-research-kb/dev/data/query_output/testing/for_nick/Nick_dev/R_Script_Dev/pull_act_for_specific_smiles.r"
  input_arg_flag <- "--smiles"
  input_str <- paste(smiles_vals, collapse = "|")

} else if (!is.null(params$smiles_file)) {
  smiles_vals <- read_compound_file_column(params$smiles_file, "isoSmiles")
  input_type <- "smiles"
  target_script <- "/sensorium-research-kb/dev/data/query_output/testing/for_nick/Nick_dev/R_Script_Dev/pull_act_for_specific_smiles.r"
  input_arg_flag <- "--smiles"
  input_str <- paste(smiles_vals, collapse = "|")
}

if (is.null(input_str) || is.null(target_script)) {
  stop("❌ No valid input provided. Please use --cmp, --cmp_file, --smiles, or --smiles_file.")
}

# 🧪 DEBUG: Print input string being passed
cat("✅ Collapsed input string:\n", input_str, "\n")

# ------------------ Build & Run Command ------------------
out_file <- file.path(params$outdir, paste0("step1_", input_type, "_acts.csv"))

cmd <- paste(
  "Rscript", target_script,
  "--endpoint", shQuote(params$endpoint),
  input_arg_flag, shQuote(input_str),
  "--out", shQuote(out_file)
)

cat("TRUE→ Running command:\n", cmd, "\n")
exit_code <- system(cmd)

if (exit_code != 0 || !file.exists(out_file)) {
  stop("❌ Script failed or output file not generated.")
}

cat("✅ Step 1 complete. Output saved to:\n", out_file, "\n")
