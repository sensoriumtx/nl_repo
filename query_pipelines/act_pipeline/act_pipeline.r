#!/usr/bin/env Rscript

suppressMessages(library(tidyverse))
suppressMessages(library(lubridate))
suppressMessages(library(parallel))

cat("Script starting...\n", flush = TRUE)

# ------------------------- Argument Parsing -------------------------
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
    if (arg == "--endpoint") {
      out$endpoint <- val
    } else if (arg == "--acts") {
      out$acts <- val
    } else if (arg == "--act_file") {
      out$act_file <- val
    } else if (arg == "--filter_out_act") {
      out$filter <- val
    } else if (arg == "--outdir") {
      out$outdir <- val
    } else if (arg == "--scoring") {
      out$scoring <- val
    } else {
      stop(paste("Unknown argument:", arg))
    }
    i <- i + 2
  }
  return(out)
}

params <- parseArgs(args)

# ------------------------- Output Setup -------------------------
dir.create(params$outdir, showWarnings = FALSE, recursive = TRUE)
log_file <- file.path(params$outdir, "pipeline_log.txt")
cat("Pipeline Log\n==============\n", file = log_file)

workers <- 10
log(paste("[Setup] Using", workers, "parallel workers"))

# ------------------------- Activity Loading -------------------------
if (!is.null(params$act_file)) {
  if (!file.exists(params$act_file)) stop("Provided --act_file does not exist.")
  act_table <- read_csv(params$act_file, show_col_types = FALSE)
  names_lower <- tolower(names(act_table))
  col_match <- which(names_lower %in% c("term", "terms"))
  if (length(col_match) == 0) stop("--act_file must contain a 'term' or 'terms' column")
  selected_col <- names(act_table)[col_match[1]]
  params$acts <- paste(act_table[[selected_col]], collapse = "|")
} else if (is.null(params$acts)) {
  stop("You must provide either --acts or --act_file")
}

# ------------------------- Step 1: Pull Plants for Activities -------------------------
log("[Step 1] Fetching Plants Associated with Activity")
all_acts <- strsplit(params$acts, "\\|")[[1]]
chunk_size <- 100
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
  rename(act_pln = act, act_label_pln = act_label) %>%
  distinct(pln, pln_label, act_pln, act_label_pln) %>%
  drop_na(pln, pln_label)

write_csv(plants_df, file.path(params$outdir, "step1_plants.csv"))
log(paste("[Step 1] Complete: Total Plants:", nrow(plants_df)))

# ------------------------- Step 2: Pull Compounds for Plants -------------------------
log("[Step 2] Fetching Compounds Associated with Identified Plants")

# Read the output of Step 1 to get plant labels
step1_out <- file.path(params$outdir, "step1_plants.csv")
if (!file.exists(step1_out)) stop("[Step 2] Missing step1_plants.csv")

plants_df <- read_csv(step1_out, show_col_types = FALSE) %>%
  drop_na(pln, pln_label)

unique_pln_labels <- unique(plants_df$pln_label)
if (length(unique_pln_labels) == 0) {
  stop("[Step 2] No valid plant labels found from Step 1 output.")
}

joined_plants <- paste(unique_pln_labels, collapse = "|")
step2_outfile <- file.path(params$outdir, "step2_cmp.csv")

cmd2 <- paste(
  "Rscript scripts/pull_cmp_for_pln.r",
  "--endpoint", params$endpoint,
  "--plants", shQuote(joined_plants),
  "--out", shQuote(step2_outfile)
)
log(paste("[Step 2]", cmd2))
system(cmd2)

if (!file.exists(step2_outfile)) {
  stop("[Step 2] Failed to generate output file")
}

cmp_df <- read_csv(step2_outfile, show_col_types = FALSE) %>%
  distinct(pln, cmp, cmp_labels) %>%
  inner_join(plants_df %>% distinct(pln, pln_label), by = "pln")

write_csv(cmp_df, step2_outfile)
log(paste("[Step 2] Complete: Total Compound Mappings:", nrow(cmp_df)))

# ------------------------- Step 3: Pull Activities for Compounds -------------------------
log("[Step 3] Fetching Compound-Level Activity Annotations")

if (is.null(params$acts) || nchar(params$acts) == 0) {
  stop("[Step 3] No valid acts string provided")
}

step3_chunk_dir <- file.path(params$outdir, "step3_chunks")
dir.create(step3_chunk_dir, showWarnings = FALSE, recursive = TRUE)

step3_outfile <- file.path(step3_chunk_dir, "step3_cmp_enrichment.csv")

cmd3 <- paste(
  "Rscript scripts/pull_act_for_cmp_enrich.r",
  "--endpoint", params$endpoint,
  "--acts", shQuote(params$acts),
  "--filter_out_act", shQuote(params$filter),
  "--out", shQuote(step3_outfile)
)

log(paste("[Step 3]", cmd3))
system(cmd3)

if (!file.exists(step3_outfile)) {
  log("[Step 3] Warning: Compound activity enrichment file not found")
  cmp_act_df <- tibble(cmp = character(), act_cmp = character(), act_label_cmp = character())
} else {
  cmp_act_df <- read_csv(step3_outfile, show_col_types = FALSE) %>%
    rename(act_cmp = act.x, act_label_cmp = act_label) %>%
    distinct(cmp, act_cmp, act_label_cmp) %>%
    drop_na(cmp)
  log(paste("[Step 3] Complete: Total Compound Activities:", nrow(cmp_act_df)))
}


# ------------------------- Step 4: Merge All -------------------------
log("[Step 4] Merging Plant, Compound, and Activity Data")

# Read cmp_df and plants_df if not in memory
if (!exists("cmp_df")) {
  cmp_df <- read_csv(file.path(params$outdir, "step2_cmp.csv"), show_col_types = FALSE)
}
if (!exists("plants_df")) {
  plants_df <- read_csv(file.path(params$outdir, "step1_plants.csv"), show_col_types = FALSE)
}

# Merge blocks
final_df <- cmp_df %>%
  full_join(plants_df, by = "pln") %>%
  full_join(cmp_act_df, by = "cmp")

# Patch to avoid logical() issue
cmp_unmapped <- cmp_act_df %>%
  filter(!(cmp %in% cmp_df$cmp)) %>%
  mutate(
    pln = as.character(NA),
    pln_label = as.character(NA),
    act_pln = as.character(NA),
    act_label_pln = as.character(NA)
  )

final_df <- bind_rows(final_df, cmp_unmapped) %>%
  left_join(cmp_df %>% select(cmp, cmp_labels_source = cmp_labels) %>% distinct(), by = "cmp") %>%
  mutate(cmp_labels = coalesce(cmp_labels, cmp_labels_source)) %>%
  select(-cmp_labels_source) %>%
  left_join(plants_df %>% select(pln, pln_label) %>% distinct(), by = "pln", suffix = c("", ".plndf")) %>%
  mutate(pln_label = coalesce(pln_label, pln_label.plndf)) %>%
  select(-pln_label.plndf) %>%
  distinct(cmp, cmp_labels, pln, pln_label, act_cmp, act_label_cmp, act_pln, act_label_pln)

write_csv(final_df, file.path(params$outdir, "final_merged_output.csv"))
log(paste("[Step 4] Complete: Final Merged File Rows:", nrow(final_df)))

# ------------------------- Step 5: Deliverable File + Scoring -------------------------
log("[Step 5] Aggregating and Merging Final Output with Scoring")

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
  log("[Step 5] Reading and merging scoring file on 'cmp'")
  scoring_df <- read_csv(params$scoring, show_col_types = FALSE)

  if ("cmp" %in% names(scoring_df)) {
    scoring_df <- scoring_df %>% mutate(cmp = as.character(cmp))
    deliverable_df <- deliverable_df %>% mutate(cmp = as.character(cmp)) %>%
      left_join(scoring_df, by = "cmp")
    log("[Step 5] Merge with scoring file successful")
  } else {
    log("[Step 5] Scoring file missing 'cmp' column — skipping merge")
  }
}

write_csv(deliverable_df, file.path(params$outdir, "final_deliverable.csv"))
log("[Step 5] Complete: Final Deliverable File Saved")

# ------------------------- Done -------------------------
log("[Pipeline] Execution complete.")
quit(save = "no")
