suppressMessages(library(tidyverse))
suppressMessages(library(lubridate))
suppressMessages(library(parallel))

cat("Script starting...\n", flush = TRUE)

args <- commandArgs(trailingOnly = TRUE)
log_file <- NULL
log <- function(message) {
  cat(message, "\n", flush = TRUE)
  if (!is.null(log_file)) {
    cat(message, "\n", file = log_file, append = TRUE)
  }
}
log(paste("Args:", paste(args, collapse = " | ")))

# ------------------------- Argument Parsing -------------------------
parseArgs <- function(args) {
  out <- list(endpoint = "dev")  # default
  i <- 1
  while (i <= length(args)) {
    arg <- args[i]
    val <- if (i + 1 <= length(args)) args[i + 1] else NA
    if (is.na(val)) stop(paste("Missing value for", arg))
    if (arg == "--cmp") {
      out$cmp <- val
    } else if (arg == "--smiles") {
      out$smiles <- val
    } else if (arg == "--inchi") {
      out$inchi <- val
    } else if (arg == "--inchi_key") {
      out$inchi_key <- val
    } else if (arg == "--iupac") {
      out$iupac <- val
    } else if (arg == "--cmp_label") {
      out$cmp_label <- val
    } else if (arg == "--in_file") {
      out$in_file <- val
    } else if (arg == "--endpoint") {
      out$endpoint <- val
    } else if (arg == "--outdir") {
      out$outdir <- val
    } else {
      stop(paste("Unknown argument:", arg))
    }
    i <- i + 2
  }
  return(out)
}

params <- parseArgs(args)
if (is.null(params$outdir)) stop("--outdir must be specified")

dir.create(params$outdir, showWarnings = FALSE, recursive = TRUE)
log_file <- file.path(params$outdir, "pipeline_log.txt")
cat("Pipeline Log\n==============\n", file = log_file)

# ------------------------- Step 1: Resolve Identifiers to CMP -------------------------
log("[Step 1] Resolving compound identifiers to `cmp` ID")

resolve_input <- tibble(
  cmp = params$cmp,
  primarySmiles = params$smiles,
  inchi = params$inchi,
  inchi_key = params$inchi_key,
  iupac_name = params$iupac,
  cmp_label = params$cmp_label
)

input_file_path <- file.path(params$outdir, "step1_identifier_input.csv")
write_csv(resolve_input, input_file_path)

resolved_cmp_file <- file.path(params$outdir, "step1_resolved_cmp.csv")
resolve_cmd <- paste(
  "Rscript scripts/resolve_cmp_id.r",
  "--in", shQuote(input_file_path),
  "--out", shQuote(resolved_cmp_file)
)
system(resolve_cmd)
if (!file.exists(resolved_cmp_file)) stop("Step 1 failed: cmp resolution output not found.")
log("[Step 1] Complete")

cmp_ids <- read_csv(resolved_cmp_file, show_col_types = FALSE)$cmp %>% unique() %>% paste(collapse = "|")
if (is.null(cmp_ids) || cmp_ids == "") stop("No resolved cmp IDs")

# ------------------------- Step 2: Pull Plants for CMP -------------------------
log("[Step 2] Pulling plants associated with resolved compounds")

plants_file <- file.path(params$outdir, "step2_plants_for_cmp.csv")
plant_cmd <- paste(
  "Rscript scripts/pull_plant_for_cmp_id.r",
  "--endpoint", params$endpoint,
  "--compound_activity_file", shQuote(resolved_cmp_file),
  "--cmp_id_column", "cmp",
  "--out", shQuote(plants_file)
)
system(plant_cmd)
if (!file.exists(plants_file)) stop("Step 2 failed: plant output not found.")
log("[Step 2] Complete")

# ------------------------- Step 3: Pull Acts for CMP -------------------------
log("[Step 3] Pulling activities associated with compounds")

cmp_acts_file <- file.path(params$outdir, "step3_acts_for_cmp.csv")
act_cmd <- paste(
  "Rscript scripts/pull_acts_for_cmp_id.r",
  "--endpoint", params$endpoint,
  "--compound", shQuote(cmp_ids),
  "--out", shQuote(cmp_acts_file)
)
system(act_cmd)
if (!file.exists(cmp_acts_file)) stop("Step 3 failed: cmp activities output not found.")
log("[Step 3] Complete")


### Add Direct References *********


# ------------------------- Step 4: Pull Acts for PLN -------------------------
log("[Step 4] Pulling activities associated with plants")

plant_labels <- read_csv(plants_file, show_col_types = FALSE)$pln_label %>%
  unique() %>% na.omit() %>% paste(collapse = "|")

plant_acts_file <- file.path(params$outdir, "step4_acts_for_pln.csv")
pln_act_cmd <- paste(
  "Rscript scripts/pull_acts_for_specific_pln.r",
  "--endpoint", params$endpoint,
  "--plants", shQuote(plant_labels),
  "--out", shQuote(plant_acts_file)
)

### Add Direct References *********



##### Should I enrich here? what to enrich? all cmp? all pln? panacea filter? p-value?

# ------------------------- Step 5: Final Output -------------------------
system(pln_act_cmd)
if (!file.exists(plant_acts_file)) stop("Step 4 failed: plant activities output not found.")
log("[Step 4] Complete")

log("[Pipeline] Success. Outputs written to:")
log(paste("  - Resolved CMPs:", resolved_cmp_file))
log(paste("  - Plants for CMPs:", plants_file))
log(paste("  - Activities for CMPs:", cmp_acts_file))
log(paste("  - Activities for Plants:", plant_acts_file))
