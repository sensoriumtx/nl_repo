suppressMessages(library(tidyverse))
suppressMessages(library(readr))
suppressMessages(library(stringr))
suppressMessages(library(optparse))

# ---- Parse Args ----
option_list <- list(
  make_option("--cmp", type = "character", help = "Compound IDs, |-delimited"),
  make_option("--cmp_file", type = "character", help = "Path to CSV or TSV with a 'cmp' column"),
  make_option("--smiles", type = "character", help = "SMILES strings, |-delimited"),
  make_option("--smiles_file", type = "character", help = "Path to CSV/TSV with a 'isoSmiles' column"),
  make_option("--endpoint", type = "character", help = "SPARQL endpoint"),
  make_option("--outdir", type = "character", help = "Output directory")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$endpoint) || is.null(opt$outdir)) {
  stop("❌ Must provide --endpoint and --outdir")
}

# ---- Utility Function ----
read_compound_file_column <- function(file_path, column_name) {
  message(paste0("📄 Reading column '", column_name, "' from: ", file_path))
  if (grepl("\\.tsv$", file_path, ignore.case = TRUE)) {
    df <- read_tsv(file_path, locale = locale(encoding = "UTF-8"), show_col_types = FALSE)
  } else {
    df <- read_csv(file_path, locale = locale(encoding = "UTF-8"), show_col_types = FALSE)
  }
  colnames(df) <- trimws(colnames(df))
  if (!(column_name %in% colnames(df))) {
    stop(paste0("❌ Column '", column_name, "' not found in ", file_path, ". Found: ", paste(colnames(df), collapse = ", ")))
  }
  values <- df[[column_name]] %>% as.character() %>% unique() %>% na.omit() %>% trimws()
  return(values)
}

# ---- Determine Input ----
input_type <- NULL
input_str <- NULL
input_arg_flag <- NULL
target_script <- NULL

if (!is.null(opt$cmp)) {
  cmp_vals <- str_split(opt$cmp, "\\|")[[1]] %>% unique() %>% na.omit() %>% trimws()
  input_type <- "cmp"
  target_script <- "scripts/pull_acts_for_specific_cmp_ids.r"
  input_arg_flag <- "--compound"
  input_str <- paste(cmp_vals, collapse = "|")

} else if (!is.null(opt$cmp_file)) {
  cmp_vals <- read_compound_file_column(opt$cmp_file, "cmp")
  input_type <- "cmp"
  target_script <- "scripts/pull_acts_for_specific_cmp_ids.r"
  input_arg_flag <- "--compound"
  input_str <- paste(cmp_vals, collapse = "|")

} else if (!is.null(opt$smiles)) {
  smiles_vals <- str_split(opt$smiles, "\\|")[[1]] %>% unique() %>% na.omit() %>% trimws()
  input_type <- "smiles"
  target_script <- "scripts/pull_act_for_specific_smiles.r"
  input_arg_flag <- "--smiles"
  input_str <- paste(smiles_vals, collapse = "|")

} else if (!is.null(opt$smiles_file)) {
  smiles_vals <- read_compound_file_column(opt$smiles_file, "isoSmiles")
  input_type <- "smiles"
  target_script <- "scripts/pull_act_for_specific_smiles.r"
  input_arg_flag <- "--smiles"
  input_str <- paste(smiles_vals, collapse = "|")
} else {
  stop("❌ No valid compound or SMILES input provided.")
}

# ---- Build & Run Command ----
out_file <- file.path(opt$outdir, paste0("step1_", input_type, "_acts.csv"))

cmd <- paste(
  "Rscript", target_script,
  "--endpoint", shQuote(opt$endpoint),
  input_arg_flag, shQuote(input_str),
  "--out", shQuote(out_file)
)

cat("\n→ Running command:\n", cmd, "\n")
exit_code <- system(cmd)

if (exit_code != 0 || !file.exists(out_file)) {
  stop("❌ Script failed or output file not generated.")
}

cat("✅ Step 1 complete. Output saved to:\n", out_file, "\n")
