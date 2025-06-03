#!/usr/bin/env Rscript

suppressMessages(library(tidyverse))
suppressMessages(library(stringi))
suppressMessages(library(parallel))
suppressMessages(library(lubridate))
suppressMessages(source("scripts/common.r"))
suppressMessages(source("scripts/SPARQL.R"))

options(scipen = 100, digits = 4)

# -------------------- Argument Parsing --------------------
args <- commandArgs(trailingOnly = TRUE)
parsed_args <- list(
  endpoint = "dev", 
  cmp_file = NULL, 
  cmp_column = "cmp", 
  chunk_size = 50, 
  outdir = NULL,
  out = NULL
)

i <- 1
while (i <= length(args)) {
  key <- args[i]
  val <- args[i + 1]
  if (key == "--endpoint") parsed_args$endpoint <- val
  if (key == "--cmp_file") parsed_args$cmp_file <- val
  if (key == "--cmp_column") parsed_args$cmp_column <- val
  if (key == "--chunk_size") parsed_args$chunk_size <- as.integer(val)
  if (key == "--outdir") parsed_args$outdir <- val
  if (key == "--out") parsed_args$out <- val
  i <- i + 2
}

# -------------------- Validate Inputs --------------------
if (is.null(parsed_args$cmp_file) || is.null(parsed_args$outdir) || is.null(parsed_args$out)) {
  stop("Missing required arguments: --cmp_file, --outdir, and/or --out", call. = FALSE)
}
if (!file.exists(parsed_args$cmp_file)) stop("File not found: ", parsed_args$cmp_file)
if (!dir.exists(parsed_args$outdir)) dir.create(parsed_args$outdir, recursive = TRUE)

chunk_dir <- file.path(parsed_args$outdir, "chunked")
if (!dir.exists(chunk_dir)) dir.create(chunk_dir, recursive = TRUE)

# -------------------- Load & Chunk CMPs --------------------
message("[INFO] Reading input file...")
df <- read_csv(parsed_args$cmp_file, show_col_types = FALSE)
if (!(parsed_args$cmp_column %in% colnames(df))) stop("Column not found: ", parsed_args$cmp_column)

cmp_list <- df[[parsed_args$cmp_column]] %>% unique() %>% discard(is.na)
message("[INFO] Loaded ", length(cmp_list), " unique cmp values")

chunks <- split(cmp_list, ceiling(seq_along(cmp_list) / parsed_args$chunk_size))
message("[INFO] Split into ", length(chunks), " chunks (chunk size = ", parsed_args$chunk_size, ")")

# -------------------- SPARQL Utility Functions --------------------
outcomes <- c(
  "\"Agonist\"", "\"Inhibitor\"", "\"Antagonist\"", "\"Activator\"",
  "\"Inhibition\"", "\"Blocker\"", "\"Channel blocker\"", "\"Blocker (channel blocker)\"",
  "\"Full agonist\"", "\"Partial agonist\"", "\"Activation\"", "\"Inhibitor (gating inhibitor)\"",
  "\"Stimulator\"", "\"Inducer\"", "\"Positive\"", "\"Pore blocker\"", "\"Suppressor\"",
  "\"Potentiation\"", "\"Gating inhibitor\"", "\"Opener\"", "\"Active\"", "\"Unspecified\""
)

get_cmp_metadata <- function(cmp_pipe_joined) {
  ids_space_sep <- gsub("\\|", " ", cmp_pipe_joined)
  q <- paste0(sparql_prefix, "
    select distinct ?cmp (sample(?label) as ?cmp_label) where {
      values ?cmp { ", ids_space_sep, " }
      {
        ?cmp (^sen:hasCompound|(sen:maps_to+/^sen:hasCompound))/rdf:type sen:use .
      } UNION {
        values ?outcome { ", paste(outcomes, collapse = " "), " }
        ?assayResult a sen:assayResult .
        ?assayResult sen:measureOutcome ?outcome .
        ?cmp ^sen:hasCompound ?assayResult .
        ?assayResult sen:targetGene|sen:targetProtein ?target .
        ?target ^sen:hasTarget/rdf:type sen:use .
      }
      ?cmp sen:lcLabel ?label
    }
    group by ?cmp
  ")
  tryCatch({
    SPARQL(endpoint, q, ns = prefix, extra = query_options, format = 'json')$results
  }, error = function(e) {
    message("[ERROR] Metadata SPARQL failed: ", e$message)
    tibble()
  })
}

get_cmp_activities <- function(cmp_id) {
  q <- paste0(sparql_prefix, "
    select distinct ?cmp ?isoSmiles ?inchi ?act ?act_label (group_concat(distinct ?path; separator=\"|\") as ?path)
           (count(distinct ?tgt) as ?target_count_for_this_compound_activity)
           (group_concat(distinct ?tgt; separator=\"|\") as ?target) ?assay_source where {
      select ?assay_source ?cmp ?isoSmiles ?inchi ?act ?act_label ?path ?tgt ?primary_symbol where {
        ?use a sen:use .
        values ?cmp { ", cmp_id, " } .
        {
          ?cmp (^sen:hasCompound|(sen:maps_to+/^sen:hasCompound)) ?use .
        } UNION {
          values ?outcome { ", paste(outcomes, collapse = " "), " }
          ?assayResult a sen:assayResult .
          ?assayResult sen:measureOutcome ?outcome .
          ?cmp ^sen:hasCompound ?assayResult .
          ?assayResult (sen:targetGene/sen:maps_to*)|(sen:targetProtein/sen:maps_to*) ?tgt_entity .
          ?use sen:hasTarget/sen:maps_to* ?tgt_entity .
          ?tgt_entity sen:primary_symbol ?tgt .
          ?tgt_entity sen:tox21 \"false\" .
          ?assayResult sen:hasSource ?src .
          ?src rdfs:label ?assay_source .
        }
        ?use sen:hasActivity ?act .
        ?act sen:lcLabel ?act_label .
        OPTIONAL { ?cmp sen:IsomericSMILES ?isoSmiles . }
        OPTIONAL { ?cmp sen:InChI ?inchi . }
        BIND(IF(bound(?tgt), \"indirect\", \"direct\") as ?path)
      } order by ?path ?cmp ?act
    } group by ?assay_source ?cmp ?isoSmiles ?inchi ?act ?act_label
  ")
  tryCatch({
    SPARQL(endpoint, q, ns = prefix, extra = query_options, format = 'json')$results
  }, error = function(e) {
    message("[ERROR] Activity SPARQL failed for ", cmp_id, ": ", e$message)
    tibble()
  })
}

# -------------------- Execute Each Chunk --------------------
endpoint <- if (parsed_args$endpoint == "dev") endpoint_dev else stop("Invalid endpoint")

chunk_files <- c()
for (i in seq_along(chunks)) {
  chunk_vals <- chunks[[i]]
  cmp_arg <- paste(chunk_vals, collapse = "|")
  out_file <- file.path(chunk_dir, paste0("chunk_", i, ".csv"))
  chunk_files <- c(chunk_files, out_file)

  message("[INFO] Running chunk ", i, " with ", length(chunk_vals), " compounds")

  meta_df <- get_cmp_metadata(cmp_arg)
  if (nrow(meta_df) == 0) {
    message("[WARN] No metadata returned for chunk ", i)
    next
  }

  activity_df <- bind_rows(lapply(meta_df$cmp, get_cmp_activities))
  final_df <- fix_sparql_ids(activity_df) %>%
    left_join(meta_df, by = "cmp")

  write_csv(final_df, out_file)
  message("[INFO] Chunk ", i, " saved to: ", out_file)
}

# -------------------- Final Merge --------------------
message("[INFO] Merging ", length(chunk_files), " chunk files into final output")

merged_df <- bind_rows(lapply(chunk_files, read_csv, show_col_types = FALSE))
write_csv(merged_df, parsed_args$out)

message("✅ All chunks complete. Final file written to: ", parsed_args$out)
