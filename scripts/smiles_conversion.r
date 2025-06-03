#!/usr/bin/env Rscript
install.packages(argparse)
install.packages(rcdk)
install.packages(dplyr)
install.packages(readr)
library(argparse)
library(rcdk)
library(dplyr)
library(readr)

# Function to check if SMILES is isomeric
is_isomeric_smiles <- function(smiles) {
  tryCatch({
    mol <- parse.smiles(smiles)[[1]]
    if (is.null(mol)) {
      return(list(is_isomeric = FALSE, status = "Invalid SMILES"))
    }
    # Convert to SMILES with stereochemistry to check if it differs
    smiles_canonical <- generate.smiles(mol, canonical = TRUE)
    smiles_isomeric <- generate.smiles(mol, canonical = TRUE, do.isomeric = TRUE)
    is_isomeric <- smiles_canonical != smiles_isomeric
    return(list(is_isomeric = is_isomeric, status = "Valid"))
  }, error = function(e) {
    return(list(is_isomeric = FALSE, status = paste("Error:", e$message)))
  })
}

# Function to standardize SMILES to isomeric form
standardize_smiles <- function(smiles) {
  tryCatch({
    mol <- parse.smiles(smiles)[[1]]
    if (is.null(mol)) {
      return(list(smiles = NA, status = "Invalid SMILES"))
    }
    # Generate isomeric SMILES
    isomeric_smiles <- generate.smiles(mol, canonical = TRUE, do.isomeric = TRUE)
    return(list(smiles = isomeric_smiles, status = "Converted"))
  }, error = function(e) {
    return(list(smiles = NA, status = paste("Error:", e$message)))
  })
}

# Function to process CSV file
process_smiles_csv <- function(input_file, output_file, smiles_column = "primary_SMILES") {
  # Read input CSV
  tryCatch({
    df <- read_csv(input_file, show_col_types = FALSE)
    if (!(smiles_column %in% colnames(df))) {
      stop(paste("Column", smiles_column, "not found in input CSV"))
    }
  }, error = function(e) {
    cat(paste("Error reading input CSV:", e$message, "\n"))
    return()
  })

  # Initialize result columns
  df$isomeric_SMILES <- NA_character_
  df$is_already_isomeric <- FALSE
  df$status <- NA_character_

  # Process each SMILES
  for (i in 1:nrow(df)) {
    smiles <- df[[i, smiles_column]]
    if (is.na(smiles) || !is.character(smiles)) {
      df$status[i] <- "Invalid or missing SMILES"
      next
    }

    # Check if SMILES is already isomeric
    isomeric_result <- is_isomeric_smiles(smiles)
    df$is_already_isomeric[i] <- isomeric_result$is_isomeric
    df$status[i] <- isomeric_result$status

    # If not isomeric, convert to isomeric SMILES
    if (!isomeric_result$is_isomeric) {
      standardized <- standardize_smiles(smiles)
      df$isomeric_SMILES[i] <- standardized$smiles
      df$status[i] <- standardized$status
    } else {
      df$isomeric_SMILES[i] <- smiles  # Keep original if already isomeric
    }
  }

  # Save results to output CSV
  tryCatch({
    write_csv(df, output_file)
    cat(paste("Output saved to", output_file, "\n"))
  }, error = function(e) {
    cat(paste("Error saving output CSV:", e$message, "\n"))
  })
}

# Parse command-line arguments
parser <- ArgumentParser(description = "Convert canonical SMILES to isomeric SMILES in a CSV file")
parser$add_argument("--in", required = TRUE, help = "Path to input CSV file with SMILES column")
parser$add_argument("--out", required = TRUE, help = "Path to output CSV file")
parser$add_argument("--smiles_column", default = "primary_SMILES", help = "Name of the SMILES column in the input CSV")
args <- parser$parse_args()

# Run the pipeline
process_smiles_csv(args$in, args$out, args$smiles_column)