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
    if (arg == "--plants") {
      out$plants <- val
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

# ------------------------- Step 1: Resolve Plants -------------------------
log("[Step 1] Resolving plant identifiers to `pln` IDs")

plant_input <- if (!is.null(params$plants)) {
  str_split(params$plants, "\\|")[[1]] %>% unique()
} else if (!is.null(params$in_file)) {
  df <- read_csv(params$in_file, show_col_types = FALSE)
  cols <- names(df)
  key_col <- intersect(tolower(cols), c("plant", "plants", "pln", "pln_label"))[1]
  if (is.null(key_col)) stop("No usable plant column found in input file")
  df[[key_col]] %>% unique() %>% na.omit()
} else {
  stop("You must provide --plants or --in_file")
}

plant_str <- paste(plant_input, collapse = "|")
plant_resolved_file <- file.path(params$outdir, "step1_resolved_pln.csv")
resolve_cmd <- paste(
  "Rscript scripts/resolve_pln_id.r",  # you can implement or rename this
  "--endpoint", params$endpoint,
  "--plants", shQuote(plant_str),
  "--out", shQuote(plant_resolved_file)
)
system(resolve_cmd)
if (!file.exists(plant_resolved_file)) stop("Step 1 failed: resolved plant file not found")
log("[Step 1] Complete")

pln_ids <- read_csv(plant_resolved_file, show_col_types = FALSE)$pln %>%
  unique() %>% paste(collapse = "|")
if (is.null(pln_ids) || pln_ids == "") stop("No valid `pln` identifiers found")

# ------------------------- Step 2: Pull Compounds for Plants -------------------------
log("[Step 2] Pulling compounds associated with plants")

plant_cmp_file <- file.path(params$outdir, "step2_cmp_for_plants.csv")
cmp_cmd <- paste(
  "Rscript scripts/pull_cmp_for_pln_id.r",
  "--endpoint", params$endpoint,
  "--plants", shQuote(pln_ids),
  "--out", shQuote(plant_cmp_file)
)
system(cmp_cmd)
if (!file.exists(plant_cmp_file)) stop("Step 2 failed: compound file not found")
log("[Step 2] Complete")

# ------------------------- Step 3: Pull Activities for Plants -------------------------
log("[Step 3] Pulling activities associated with plants")

plant_acts_file <- file.path(params$outdir, "step3_acts_for_pln.csv")
acts_cmd <- paste(
  "Rscript scripts/pull_acts_for_pln_id.r",
  "--endpoint", params$endpoint,
  "--plants", shQuote(pln_ids),
  "--out", shQuote(plant_acts_file)
)
system(acts_cmd)
if (!file.exists(plant_acts_file)) stop("Step 3 failed: plant activity file not found")
log("[Step 3] Complete")

log("[Pipeline] Success. Outputs written to:")
log(paste("  - Resolved Plants:", plant_resolved_file))
log(paste("  - Compounds for Plants:", plant_cmp_file))
log(paste("  - Activities for Plants:", plant_acts_file))
