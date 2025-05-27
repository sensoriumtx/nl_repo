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

# Ensure output directory and initialize log file
dir.create(params$outdir, showWarnings = FALSE, recursive = TRUE)
log_file <- file.path(params$outdir, "pipeline_log.txt")
cat("Pipeline Log\n==============\n", file = log_file)

workers <- 10
log(paste("Using", workers, "parallel workers"))

# -------- Validate or Load Activities --------
if (!is.null(params$act_file)) {
  if (!file.exists(params$act_file)) stop("Provided --act_file does not exist.")

  act_table <- read_csv(params$act_file, show_col_types = FALSE)
  names(act_table) <- tolower(names(act_table))  # Normalize all column names to lowercase

  # Acceptable column names (case insensitive): "term", "terms"
  col_options <- c("term", "terms")
  term_col <- col_options[col_options %in% names(act_table)]

  if (length(term_col) == 0) stop("--act_file must contain a column named 'term', 'terms', 'Term', or 'Terms'")

  # Use the first matching column
  params$acts <- paste(act_table[[term_col[1]]], collapse = "|")
} else if (is.null(params$acts)) {
  stop("You must provide either --acts or --act_file")
}


# --------- Step 1: Fetching Activities for Plants (Parallelized) ---------------

if (!is.null(params$act_file)) {
  if (!file.exists(params$act_file)) stop("Provided --act_file does not exist.")
  act_table <- read_csv(params$act_file, show_col_types = FALSE)
  if (!"terms" %in% names(act_table)) stop("--act_file must contain a 'terms' column")
  params$acts <- paste(act_table$terms, collapse = "|")
} else if (is.null(params$acts)) {
  stop("You must provide either --acts or --act_file")
}
log("[Step 1] Fetching Plants for Activity")

all_acts <- strsplit(params$acts, "\\|")[[1]]
chunk_size <- 50
act_chunks <- split(all_acts, ceiling(seq_along(all_acts) / chunk_size))

cl <- makeCluster(workers)
clusterExport(cl, list("params", "act_chunks"))

chunk_output_files <- parLapply(cl, seq_along(act_chunks), function(i) {
  chunk_acts <- act_chunks[[i]]
  chunk_file <- file.path(params$outdir, paste0("step1_acts_chunk_", i, ".csv"))
  cmd <- paste(
    "Rscript scripts/pull_act_for_plant_enrich.r",
    "--endpoint", params$endpoint,
    "--acts", shQuote(paste(chunk_acts, collapse = "|")),
    "--filter_out_act", shQuote(params$filter),
    "--out", shQuote(chunk_file)
  )
  system(cmd)
  return(chunk_file)
})
stopCluster(cl)

plants_df <- map_dfr(chunk_output_files, read_csv, show_col_types = FALSE) %>%
  distinct(pln, pln_label, act_pln = act, act_label_pln = act_label) %>%
  drop_na(pln, pln_label)
step1_out <- file.path(params$outdir, "step1_plants.csv")
write_csv(plants_df, step1_out)
log(paste("[Step 1] Complete: Total Plants:", nrow(plants_df)))

# --------- Step 2: Fetching Compounds for Plants (Fully Parallelized) ----------
log("[Step 2] Fetching Compounds for Plants")

plant_label_map <- plants_df %>% distinct(pln, pln_label)
plant_chunks <- split(plant_label_map$pln_label, ceiling(seq_along(plant_label_map$pln_label) / 100))
log(paste("Total Plant Labels:", nrow(plant_label_map)))
log(paste("Total Chunks:", length(plant_chunks)))

step2_chunk_dir <- file.path(params$outdir, "step2_chunks")
dir.create(step2_chunk_dir, showWarnings = FALSE, recursive = TRUE)

run_chunk <- function(i) {
  chunk_file <- file.path(step2_chunk_dir, paste0("step2_cmp_chunk_", i, ".csv"))
  chunk_plants <- paste(plant_chunks[[i]], collapse = "|")
  cmd2 <- paste(
    "Rscript scripts/pull_cmp_for_pln.r",
    "--endpoint", params$endpoint,
    "--plants", shQuote(chunk_plants),
    "--out", shQuote(chunk_file)
  )
  system(cmd2)
  return(chunk_file)
}

chunk_output_files <- parallel::mclapply(seq_along(plant_chunks), run_chunk, mc.cores = workers)
valid_chunks <- chunk_output_files[sapply(chunk_output_files, file.exists)]

log(paste("[Step 2] Total Valid Chunks:", length(valid_chunks)))
cmp_df <- map_dfr(valid_chunks, read_csv, show_col_types = FALSE) %>%
  distinct(pln, cmp, cmp_labels) %>%
  inner_join(plant_label_map, by = "pln")
step2_out <- file.path(params$outdir, "step2_cmp.csv")
write_csv(cmp_df, step2_out)
log(paste("[Step 2] Complete: Total Compound Mappings:", nrow(cmp_df)))

# --------- Step 3: Fetching Activities for Compounds ---------------------------
log("[Step 3] Fetching Compounds for Activities")

act_chunks_cmp <- split(all_acts, ceiling(seq_along(all_acts) / chunk_size))

step3_chunk_dir <- file.path(params$outdir, "step3_chunks")
dir.create(step3_chunk_dir, showWarnings = FALSE, recursive = TRUE)

run_cmp_act_chunk <- function(i) {
  chunk_file <- file.path(step3_chunk_dir, paste0("step3_cmp_chunk_", i, ".csv"))
  chunk_acts <- act_chunks_cmp[[i]]
  cmd3 <- paste(
    "Rscript scripts/pull_act_for_cmp_enrich.r",
    "--endpoint", params$endpoint,
    "--acts", shQuote(paste(chunk_acts, collapse = "|")),
    "--filter_out_act", shQuote(params$filter),
    "--out", shQuote(chunk_file)
  )
  system(cmd3)
  return(chunk_file)
}

cmp_act_chunk_files <- parallel::mclapply(seq_along(act_chunks_cmp), run_cmp_act_chunk, mc.cores = workers)
valid_cmp_chunks <- cmp_act_chunk_files[sapply(cmp_act_chunk_files, file.exists)]

cmp_act_df <- map_dfr(valid_cmp_chunks, read_csv, show_col_types = FALSE) %>%
  distinct(cmp, act_cmp = act.x, act_label_cmp = act_label) %>%
  drop_na(cmp)
log(paste("[Step 3] Complete: Total Compounds:", nrow(cmp_act_df)))

# --------- Step 4: Merging Dataframes (Expanded Merge) -------------------------
log("[Step 4] Merging Plant, Compound, and Activity Data")

final_df <- cmp_df %>%
  full_join(plants_df, by = "pln") %>%
  full_join(cmp_act_df, by = "cmp")

cmp_unmapped <- cmp_act_df %>%
  filter(!(cmp %in% cmp_df$cmp)) %>%
  mutate(pln = NA, pln_label = NA, act_pln = NA, act_label_pln = NA)

final_df <- bind_rows(final_df, cmp_unmapped) %>%
  left_join(cmp_df %>% select(cmp, cmp_labels_source = cmp_labels) %>% distinct(), by = "cmp") %>%
  mutate(cmp_labels = coalesce(cmp_labels, cmp_labels_source)) %>%
  select(-cmp_labels_source) %>%
  left_join(plants_df %>% select(pln, pln_label) %>% distinct(), by = "pln", suffix = c("", ".plndf")) %>%
  mutate(pln_label = coalesce(pln_label, pln_label.plndf)) %>%
  select(-pln_label.plndf) %>%
  distinct(cmp, cmp_labels, pln, pln_label, act_cmp, act_label_cmp, act_pln, act_label_pln)

final_out <- file.path(params$outdir, "final_merged_output.csv")
write_csv(final_df, final_out)
log(paste("[Step 4] Complete: Final Merged File Rows:", nrow(final_df)))

# --------- Step 5: Generating Clean Deliverable File with Scoring --------------
log("[Step 5] Aggregating and Scoring Final Output")

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

step5_out <- file.path(params$outdir, "deliverable.csv")
write_csv(deliverable_df, step5_out)
log("[Step 5] Complete: Deliverable File with Scoring Saved")

# --------- Step 6: Run Classification on Compounds -----------------------------
log("[Step 6] Running Classification Query for Compounds")

step6_out <- file.path(params$outdir, "final_deliverable.csv")
cmd6 <- paste(
  "Rscript /sensorium-research-kb/dev/data/query_output/testing/for_nick/scripts/pull_class_for_cmp.r",
  "--endpoint", params$endpoint,
  "--in", shQuote(step5_out),
  "--out", shQuote(step6_out)
)
system(cmd6)

log("[Step 6] Complete: Final Deliverable with Classifications Saved")

# Patch cmp_label if missing using cmpLabel from classification output
if (file.exists(step6_out)) {
  final_classified_df <- read_csv(step6_out, show_col_types = FALSE)
  if ("cmpLabel" %in% names(final_classified_df)) {
    final_classified_df <- final_classified_df %>%
      mutate(cmp_label = coalesce(cmp_label, cmpLabel))
    write_csv(final_classified_df, step6_out)
    log("[Step 6] Patched: Missing cmp_label Values Mapped")
  }
}
log("[Pipeline] Execution complete.")
quit(save = "no")
