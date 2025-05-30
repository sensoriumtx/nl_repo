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

# ------------------ Validate Required Args ------------------
if (is.null(params$outdir)) stop("❌ --outdir must be specified")
dir.create(params$outdir, recursive = TRUE, showWarnings = FALSE)

if (is.null(params$endpoint)) stop("❌ --endpoint must be specified")

if (is.null(params$cmp) && is.null(params$cmp_file) &&
    is.null(params$smiles) && is.null(params$smiles_file)) {
  stop("❌ You must provide one of: --cmp, --cmp_file, --smiles, or --smiles_file")
}

# ------------------ File Reading Helper ------------------
read_compound_file_column <- function(file_path, column_name) {
  message(paste0("📄 Reading column '", column_name, "' from: ", file_path))
  
  if (grepl("\\.tsv$", file_path, ignore.case = TRUE)) {
    df <- read_tsv(file_path, locale = locale(encoding = "UTF-8"), show_col_types = FALSE)
  } else {
    df <- read_csv(file_path, locale = locale(encoding = "UTF-8"), show_col_types = FALSE)
  }

  colnames(df) <- trimws(colnames(df))  # Clean column names

  if (!(column_name %in% colnames(df))) {
    stop(paste0("❌ Column '", column_name, "' not found in ", file_path, "\nAvailable columns: ", paste(colnames(df), collapse = ", ")))
  }

  values <- df[[column_name]] %>%
    as.character() %>%
    unique() %>%
    na.omit() %>%
    trimws()

  if (length(values) == 0) {
    stop(paste0("❌ No valid entries found in column '", column_name, "'"))
  }

  return(values)
}

# ------------------ Determine Input Source ------------------
if (!is.null(params$cmp)) {
  compound_vals <- str_split(params$cmp, "\\|")[[1]] %>% unique() %>% na.omit() %>% trimws()
  input_type <- "cmp"
  target_script <- "scripts/pull_acts_for_specific_cmp_ids.r"

} else if (!is.null(params$cmp_file)) {
  compound_vals <- read_compound_file_column(params$cmp_file, "cmp")
  input_type <- "cmp"
  target_script <- "scripts/pull_acts_for_specific_cmp_ids.r"

} else if (!is.null(params$smiles)) {
  compound_vals <- str_split(params$smiles, "\\|")[[1]] %>% unique() %>% na.omit() %>% trimws()
  input_type <- "smiles"
  target_script <- "/sensorium-research-kb/dev/data/query_output/testing/for_nick/Nick_dev/R_Script_Dev/pull_act_for_specific_smiles.r"

} else if (!is.null(params$smiles_file)) {
  compound_vals <- read_compound_file_column(params$smiles_file, "isoSmiles")
  input_type <- "smiles"
  target_script <- "/sensorium-research-kb/dev/data/query_output/testing/for_nick/Nick_dev/R_Script_Dev/pull_act_for_specific_smiles.r"
}

# ------------------ Prepare Input and Run ------------------
input_str <- paste(compound_vals, collapse = "|")
input_arg_flag <- "--compound"
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
  stop("❌ Script failed or output file not generated.")
}

cat("✔ Step 1 complete. Output saved to:\n", out_file, "\n")
