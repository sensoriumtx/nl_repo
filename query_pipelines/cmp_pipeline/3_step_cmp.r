# ------------------ Determine Input ------------------
if (!is.null(params$cmp)) {
  compound_vals <- str_split(params$cmp, "\\|")[[1]] %>% unique() %>% na.omit() %>% trimws()
  input_type <- "cmp"
  target_script <- "scripts/pull_acts_for_specific_cmp_ids.r"

} else if (!is.null(params$cmp_file)) {
  df <- read_csv(params$cmp_file, show_col_types = FALSE)
  if (!"cmp" %in% names(df)) stop("❌ cmp_file must have a column named 'cmp'")
  compound_vals <- df$cmp %>% unique() %>% na.omit() %>% trimws()
  input_type <- "cmp"
  target_script <- "scripts/pull_acts_for_specific_cmp_ids.r"

} else if (!is.null(params$smiles)) {
  compound_vals <- str_split(params$smiles, "\\|")[[1]] %>% unique() %>% na.omit() %>% trimws()
  input_type <- "smiles"
  target_script <- "/sensorium-research-kb/dev/data/query_output/testing/for_nick/Nick_dev/R_Script_Dev/pull_act_for_specific_smiles.r"

} else if (!is.null(params$smiles_file)) {
  df <- read_csv(params$smiles_file, show_col_types = FALSE)
  if (!"isoSmiles" %in% names(df)) stop("❌ smiles_file must have a column named 'isoSmiles'")
  compound_vals <- df$isoSmiles %>% unique() %>% na.omit() %>% trimws()
  input_type <- "smiles"
  target_script <- "/sensorium-research-kb/dev/data/query_output/testing/for_nick/Nick_dev/R_Script_Dev/pull_act_for_specific_smiles.r"

} else {
  stop("❌ One of --cmp, --cmp_file, --smiles, or --smiles_file must be provided.")
}

# Collapse the values to be passed to --compound
input_str <- paste(compound_vals, collapse = "|")
input_arg_flag <- "--compound"

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
  stop("❌ Script failed or output file not generated.")
}

cat("✔ Step 1 complete. Output saved to:\n", out_file, "\n")
