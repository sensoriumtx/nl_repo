#!/usr/bin/env Rscript

# --- Install Required Packages ---
packages <- c("argparse", "rcdk", "rcdklibs", "rJava", "dplyr", "readr")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
    library(pkg, character.only = TRUE)
  }
}

# --- Load Libraries Silently ---
suppressMessages({
  library(argparse)
  library(rcdk)
  library(rcdklibs)
  library(rJava)
  library(dplyr)
  library(readr)
})

# --- Functions ---

# Check if SMILES is already isomeric
is_isomeric_smiles <- function(smiles) {
  tryCatch({
    mol <- parse.smiles(smiles)[[1]]
    if (is.null(mol)) return(list(is_isomeric = FALSE, status = "Invalid SMILES"))
    canonical <- get.smiles(mol, flavor = smiles.flavors(c("Canonical")))
    isomeric <- get.smiles(mol, flavor = smiles.flavors(c("Canonical", "Isomeric")))
    return(list(is_isomeric = canonical != isomeric, status = "Valid"))
  }, error = function(e) list(is_isomeric = FALSE, status = paste("Error:", e$message)))
}

# Convert SMILES to isomeric format
standardize_smiles <- function(smiles) {
  tryCatch({
    mol <- parse.smiles(smiles)[[1]]
    if (is.null(mol)) return(list(smiles = NA, status = "Invalid SMILES"))
    iso <- get.smiles(mol, flavor = smiles.flavors(c("Canonical", "Isomeric")))
    return(list(smiles = iso, status = "Converted"))
  }, error = function(e) list(smiles = NA, status = paste("Error:", e$message)))
}

# Process SMILES from CSV input
process_smiles_csv <- function(input_file, output_file, smiles_column = "primary_SMILES") {
  df <- tryCatch(read_csv(input_file, show_col_types = FALSE), 
                 error = function(e) stop("Error reading CSV: ", e$message))

  if (!(smiles_column %in% colnames(df))) {
    stop("Column ", smiles_column, " not found in the input CSV.")
  }

  df <- df %>%
    mutate(
      isomeric_SMILES = NA_character_,
      is_already_isomeric = FALSE,
      status = NA_character_
    )

  for (i in seq_len(nrow(df))) {
    smiles <- df[[i, smiles_column]]
    if (is.na(smiles) || !is.character(smiles)) {
      df$status[i] <- "Invalid or missing SMILES"
      next
    }

    iso_check <- is_isomeric_smiles(smiles)
    df$is_already_isomeric[i] <- iso_check$is_isomeric
    df$status[i] <- iso_check$status

    if (!iso_check$is_isomeric) {
      result <- standardize_smiles(smiles)
      df$isomeric_SMILES[i] <- result$smiles
      df$status[i] <- result$status
    } else {
      df$isomeric_SMILES[i] <- smiles
    }
  }

  # Prevent writing to directory instead of file
  if (dir.exists(output_file)) {
    stop("ERROR: The --out path appears to be a directory. Please specify a full file name like 'output.csv'")
  }

  write_csv(df, output_file)
  cat("Output saved to:", output_file, "\n")
}

# Process inline SMILES string(s)
process_inline_smiles <- function(smiles_input) {
  smiles_list <- unlist(strsplit(smiles_input, "\\|"))
  results <- data.frame(
    input_SMILES = smiles_list,
    isomeric_SMILES = NA_character_,
    is_already_isomeric = FALSE,
    status = NA_character_,
    stringsAsFactors = FALSE
  )

  for (i in seq_along(smiles_list)) {
    smi <- smiles_list[i]
    check <- is_isomeric_smiles(smi)
    results$is_already_isomeric[i] <- check$is_isomeric
    results$status[i] <- check$status

    if (!check$is_isomeric) {
      conv <- standardize_smiles(smi)
      results$isomeric_SMILES[i] <- conv$smiles
      results$status[i] <- conv$status
    } else {
      results$isomeric_SMILES[i] <- smi
    }
  }

  print(results)
}

# --- Argument Parser Setup ---
parser <- ArgumentParser(description = "Standardize SMILES to isomeric form")

parser$add_argument("--in", dest = "in_file", help = "Path to input CSV file")
parser$add_argument("--out", dest = "out_file", help = "Path to output CSV file")
parser$add_argument("--smiles_column", default = "primary_SMILES", help = "SMILES column name in input CSV")
parser$add_argument("--smiles", help = "SMILES string(s) separated by pipe (|)")

args <- parser$parse_args()

# --- Execution Logic ---
if (!is.null(args$smiles)) {
  process_inline_smiles(args$smiles)
} else if (!is.null(args$in_file) && !is.null(args$out_file)) {
  process_smiles_csv(args$in_file, args$out_file, args$smiles_column)
} else {
  cat("ERROR: Provide either --smiles or both --in and --out paths.\n")
  parser$print_help()
}
