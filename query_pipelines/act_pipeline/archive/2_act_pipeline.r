suppressMessages(library(tidyverse))
suppressMessages(library(lubridate))
suppressMessages(library(parallel))
suppressMessages(library(furrr))
suppressMessages(library(future))

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

parseArgs <- function(args) {
  out <- list()
  i <- 1
  while (i <= length(args)) {
    arg <- args[i]
    val <- if (i + 1 <= length(args)) args[i + 1] else NA
    if (is.na(val)) stop(paste("Missing value for", arg))
    switch(arg,
           "--endpoint" = { out$endpoint <- val },
           "--acts" = { out$acts <- val },
           "--act_file" = { out$act_file <- val },
           "--filter_out_act" = { out$filter <- val },
           "--outdir" = { out$outdir <- val },
           "--scoring" = { out$scoring <- val },
           stop(paste("Unknown argument:", arg))
    )
    i <- i + 2
  }
  return(out)
}

params <- parseArgs(args)

dir.create(params$outdir, showWarnings = FALSE, recursive = TRUE)
log_file <- file.path(params$outdir, "pipeline_log.txt")
cat("Pipeline Log\n==============\n", file = log_file)

# Set optimized worker count
workers <- max(1, floor(parallel::detectCores(logical = FALSE) * 0.75))
log(paste("Detected cores:", detectCores(), "| Using workers:", workers))

# -------- Load Activity Terms --------
if (!is.null(params$act_file)) {
  if (!file.exists(params$act_file)) stop("Provided --act_file does not exist.")
  act_table <- read_csv(params$act_file, show_col_types = FALSE)
  names(act_table) <- tolower(names(act_table))
  term_col <- intersect(names(act_table), c("term", "terms"))
  if (length(term_col) == 0) stop("--act_file must contain a column named 'term', 'terms', 'Term', or 'Terms'")
  params$acts <- paste(act_table[[term_col[1]]], collapse = "|")
} else if (is.null(params$acts)) {
  stop("You must provide either --acts or --act_file")
}

# -------- Future Plan --------
plan(multisession, workers = workers)

# -------- Chunk and Run Functions --------
batch_execute <- function(items, batch_size, fn) {
  results <- list()
  for (i in seq(1, length(items), by = batch_size)) {
    range <- i:min(i + batch_size - 1, length(items))
    log(paste("Processing batch", range[1], "to", range[length(range)]))
    batch_results <- future_map(range, fn, .progress = TRUE)
    results <- c(results, batch_results)
  }
  return(results)
}

# -------- Step 1 --------
log("[Step 1] Fetching Plants for Activity")
all_acts <- strsplit(params$acts, "\\|")[[1]]
chunk_size <- 25
act_chunks <- split(all_acts, ceiling(seq_along(all_acts) / chunk_size))

step1_chunk_dir <- file.path(params$outdir, "step1_chunks")
dir.create(step1_chunk_dir, showWarnings = FALSE, recursive = TRUE)

run_act_chunk <- function(i) {
  chunk_file <- file.path(step1_chunk_dir, paste0("step1_acts_chunk_", i, ".csv"))
  chunk_acts <- act_chunks[[i]]
  cmd <- paste(
    "Rscript scripts/pull_act_for_plant_enrich.r",
    "--endpoint", params$endpoint,
    "--acts", shQuote(paste(chunk_acts, collapse = "|")),
    "--filter_out_act", shQuote(params$filter),
    "--out", shQuote(chunk_file)
  )
  log(paste("Executing chunk", i))
  system(cmd)
  return(chunk_file)
}

step1_files <- batch_execute(seq_along(act_chunks), 3, run_act_chunk)
valid_step1_files <- step1_files[file.exists(unlist(step1_files))]
plants_df <- map_dfr(valid_step1_files, read_csv, show_col_types = FALSE) %>%
  distinct(pln, pln_label, act_pln = act, act_label_pln = act_label) %>%
  drop_na(pln, pln_label)
write_csv(plants_df, file.path(params$outdir, "step1_plants.csv"))

# -------- Step 2 --------
log("[Step 2] Fetching Compounds for Plants")
plant_chunks <- split(plants_df$pln_label, ceiling(seq_along(plants_df$pln_label) / 100))
step2_chunk_dir <- file.path(params$outdir, "step2_chunks")
dir.create(step2_chunk_dir, showWarnings = FALSE, recursive = TRUE)

run_cmp_chunk <- function(i) {
  chunk_file <- file.path(step2_chunk_dir, paste0("step2_cmp_chunk_", i, ".csv"))
  cmd <- paste(
    "Rscript scripts/pull_cmp_for_pln.r",
    "--endpoint", params$endpoint,
    "--plants", shQuote(paste(plant_chunks[[i]], collapse = "|")),
    "--out", shQuote(chunk_file)
  )
  system(cmd)
  return(chunk_file)
}

step2_files <- batch_execute(seq_along(plant_chunks), 3, run_cmp_chunk)
valid_step2_files <- step2_files[file.exists(unlist(step2_files))]
cmp_df <- map_dfr(valid_step2_files, read_csv, show_col_types = FALSE) %>%
  distinct(pln, cmp, cmp_labels) %>%
  inner_join(plants_df, by = "pln")
write_csv(cmp_df, file.path(params$outdir, "step2_cmp.csv"))

# -------- Step 3 --------
log("[Step 3] Fetching Activities for Compounds")
act_chunks_cmp <- split(all_acts, ceiling(seq_along(all_acts) / chunk_size))
step3_chunk_dir <- file.path(params$outdir, "step3_chunks")
dir.create(step3_chunk_dir, showWarnings = FALSE, recursive = TRUE)

run_cmp_act_chunk <- function(i) {
  chunk_file <- file.path(step3_chunk_dir, paste0("step3_cmp_chunk_", i, ".csv"))
  cmd <- paste(
    "Rscript scripts/pull_act_for_cmp_enrich.r",
    "--endpoint", params$endpoint,
    "--acts", shQuote(paste(act_chunks_cmp[[i]], collapse = "|")),
    "--filter_out_act", shQuote(params$filter),
    "--out", shQuote(chunk_file)
  )
  system(cmd)
  return(chunk_file)
}

step3_files <- batch_execute(seq_along(act_chunks_cmp), 3, run_cmp_act_chunk)
valid_step3_files <- step3_files[file.exists(unlist(step3_files))]
cmp_act_df <- map_dfr(valid_step3_files, read_csv, show_col_types = FALSE) %>%
  distinct(cmp, act_cmp = act.x, act_label_cmp = act_label) %>%
  drop_na(cmp)

# -------- Step 4 --------
log("[Step 4] Merging Dataframes")
final_df <- cmp_df %>%
  full_join(plants_df, by = "pln") %>%
  full_join(cmp_act_df, by = "cmp") %>%
  left_join(cmp_df %>% select(cmp, cmp_labels_source = cmp_labels), by = "cmp") %>%
  mutate(cmp_labels = coalesce(cmp_labels, cmp_labels_source)) %>%
  select(-cmp_labels_source) %>%
  distinct(cmp, cmp_labels, pln, pln_label, act_cmp, act_label_cmp, act_pln, act_label_pln)
write_csv(final_df, file.path(params$outdir, "final_merged_output.csv"))

# -------- Step 5 --------
log("[Step 5] Aggregating and Scoring")
deliverable_df <- final_df %>%
  group_by(cmp) %>%
  summarize(
    cmp_label = first(cmp_labels),
    pln = paste(unique(na.omit(pln)), collapse = "|"),
    pln_label = paste(unique(na.omit(pln_label)), collapse = "|"),
    act_cmp = paste(unique(na.omit(act_cmp)), collapse = "|"),
    act_label_cmp = paste(unique(na.omit(act_label_cmp)), collapse = "|"),
    act_pln = paste(unique(na.omit(act_pln)), collapse = "|"),
    act_label_pln = paste(unique(na.omit(act_label_pln)), collapse = "|"),
    .groups = "drop"
  )
if (!is.null(params$scoring)) {
  scoring_df <- read_csv(params$scoring, show_col_types = FALSE)
  deliverable_df <- deliverable_df %>% left_join(scoring_df, by = "cmp")
}
write_csv(deliverable_df, file.path(params$outdir, "deliverable.csv"))

# -------- Step 6 --------
log("[Step 6] Classification of Compounds")
cmd6 <- paste(
  "Rscript scripts/pull_class_for_cmp.r",
  "--endpoint", params$endpoint,
  "--in", shQuote(file.path(params$outdir, "deliverable.csv")),
  "--out", shQuote(file.path(params$outdir, "final_deliverable.csv"))
)
system(cmd6)

if (file.exists(file.path(params$outdir, "final_deliverable.csv"))) {
  final_classified_df <- read_csv(file.path(params$outdir, "final_deliverable.csv"), show_col_types = FALSE)
  if ("cmpLabel" %in% names(final_classified_df)) {
    final_classified_df <- final_classified_df %>% mutate(cmp_label = coalesce(cmp_label, cmpLabel))
    write_csv(final_classified_df, file.path(params$outdir, "final_deliverable.csv"))
  }
}

log("[Pipeline] Execution complete.")
quit(save = "no")
